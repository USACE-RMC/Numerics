using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Net.Sockets;
using System.Threading;
using System.Threading.Tasks;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;

namespace Data.TimeSeriesAnalysis
{
    /// <summary>
    /// Provides integration and validation tests for the <see cref="TimeSeriesDownload"/> class,
    /// including downloads from the Canadian Hydrometric Monitoring Network (CHMN),
    /// the United States Geological Survey (USGS), the Global Historical Climatology Network (GHCN),
    /// and the Australian Bureau of Meteorology (BOM).
    /// </summary>
    /// <remarks>
    /// <para>
    /// Live provider checks carry the <c>LiveServiceIntegration</c> category and are excluded from
    /// PR integration. Run them explicitly with
    /// <c>dotnet test -c Release --filter "TestCategory=LiveServiceIntegration"</c>.
    /// Offline request, parsing, and input-validation tests remain in the PR gate.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    [DoNotParallelize]
    public class Test_TimeSeriesDownload
    {
        #region Station Lists

        /// <summary>
        /// CHMN station: Cold River at Merritt.
        /// </summary>
        private const string CHMN_1 = "08LG010";

        /// <summary>
        /// CHMN station: Lillooet River near Pemberton.
        /// </summary>
        private const string CHMN_2 = "08MG005";

        /// <summary>
        /// CHMN station: Capilano River above intake.
        /// </summary>
        private const string CHMN_3 = "08GA010";

        /// <summary>
        /// USGS site number for Little River near Durham, NH (example).
        /// </summary>
        private const string USGS_1 = "01134500";

        /// <summary>
        /// USGS site number for Kaweah River near Three Rivers, CA.
        /// </summary>
        private const string USGS_2 = "11274500";

        /// <summary>
        /// USGS site number for Potomac River near Washington, DC.
        /// </summary>
        private const string USGS_3 = "01614000";

        /// <summary>
        /// USGS site number for Mississippi River at St. Louis, MO.
        /// </summary>
        private const string USGS_4 = "07010000";

        /// <summary>
        /// USGS site number for Susquehanna River at Harrisburg, PA.
        /// </summary>
        private const string USGS_5 = "01570500";

        /// <summary>
        /// USGS site number for Potomac River at Little Falls near Washington, DC.
        /// </summary>
        private const string USGS_6 = "01646500";

        /// <summary>
        /// GHCN station: USC00040741.
        /// </summary>
        private const string GHCN_1 = "USC00040741";

        /// <summary>
        /// GHCN station: USC00042402.
        /// </summary>
        private const string GHCN_2 = "USC00042402";

        /// <summary>
        /// GHCN station: USC00046685.
        /// </summary>
        private const string GHCN_3 = "USC00046685";

        /// <summary>
        /// BOM station: Cotter River at Gingera.
        /// </summary>
        private const string BOM_1 = "410730";

        /// <summary>
        /// BOM station: Goodradigbee River at Wee Jasper.
        /// </summary>
        private const string BOM_2 = "410761";

        /// <summary>
        /// BOM station: Murray River at Doctors Point.
        /// </summary>
        private const string BOM_3 = "409202";

        /// <summary>
        /// Start date for deterministic test windows.
        /// </summary>
        private static readonly DateTime WinStart = new DateTime(2021, 11, 01);

        /// <summary>
        /// End date for deterministic test windows.
        /// </summary>
        private static readonly DateTime WinEnd = new DateTime(2021, 11, 30);

        #endregion

        #region Helper Methods

        /// <summary>
        /// Checks whether the system currently has an active Internet connection.
        /// </summary>
        /// <returns>True if an Internet connection is available; otherwise false.</returns>
        private static async Task<bool> Online() => await TimeSeriesDownload.IsConnectedToInternet();

        /// <summary>
        /// Per-service availability cache. Successful probes are cached for the rest of
        /// the run. Failed probes are cached only for <see cref="ProbeFailureTtl"/>,
        /// so a single transient flake (e.g. BOM's WDP backend returning 500 /
        /// DatasourceError once) does not silently disable every test for that
        /// provider for the entire run.
        /// </summary>
        private static readonly ConcurrentDictionary<string, (bool ok, DateTime at)> _serviceState = new();

        /// <summary>
        /// How long a "service unreachable" verdict sticks before the next probe.
        /// </summary>
        private static readonly TimeSpan ProbeFailureTtl = TimeSpan.FromSeconds(30);

        /// <summary>
        /// "Is this service reachable?" probe with TTL on negative results. Returns false
        /// when the host is offline or when the supplied probe throws (any exception is
        /// treated as "service down" for the next 30 seconds).
        /// </summary>
        private static async Task<bool> AvailableAsync(string serviceName, Func<Task> probe)
        {
            if (_serviceState.TryGetValue(serviceName, out var prior))
            {
                if (prior.ok) return true;
                if (DateTime.UtcNow - prior.at < ProbeFailureTtl) return false;
            }
            if (!await Online())
            {
                _serviceState[serviceName] = (false, DateTime.UtcNow);
                return false;
            }
            try
            {
                await probe();
                _serviceState[serviceName] = (true, DateTime.UtcNow);
                return true;
            }
            catch
            {
                _serviceState[serviceName] = (false, DateTime.UtcNow);
                return false;
            }
        }

        /// <summary>
        /// Fake USGS Water Data API that answers requests through a routing function.
        /// </summary>
        /// <remarks>
        /// Requests are recorded with the API key header they carried. A route that returns null answers
        /// 404, and connectivity probes fail, so tests also prove no preflight request is made.
        /// </remarks>
        private sealed class FakeUsgsApiHandler : HttpMessageHandler
        {
            /// <summary>
            /// Maps a request URL to a response, or null for 404.
            /// </summary>
            private readonly Func<string, HttpResponseMessage> _route;

            /// <summary>
            /// Whether responses complete on the thread pool.
            /// </summary>
            private readonly bool _completeAsynchronously;

            /// <summary>
            /// Guards the request lists; continuous downloads send requests in parallel.
            /// </summary>
            private readonly object _gate = new object();

            /// <summary>
            /// Initializes a new instance of the <see cref="FakeUsgsApiHandler"/> class.
            /// </summary>
            /// <param name="route">Maps a request URL to a response, or null for 404.</param>
            /// <param name="completeAsynchronously">True to complete responses on the thread pool.</param>
            internal FakeUsgsApiHandler(Func<string, HttpResponseMessage> route, bool completeAsynchronously = false)
            {
                _route = route;
                _completeAsynchronously = completeAsynchronously;
            }

            /// <summary>
            /// Gets the absolute URLs requested through this handler.
            /// </summary>
            internal List<string> Requests { get; } = new List<string>();

            /// <summary>
            /// Gets the X-Api-Key header of each request, or null.
            /// </summary>
            internal List<string> ApiKeys { get; } = new List<string>();

            /// <summary>
            /// Returns a 200 response with a text body.
            /// </summary>
            /// <param name="body">The body.</param>
            internal static HttpResponseMessage Text(string body) => new HttpResponseMessage(HttpStatusCode.OK) { Content = new StringContent(body) };

            /// <summary>
            /// Determines whether a URL is a connectivity probe.
            /// </summary>
            /// <param name="url">The absolute URL.</param>
            /// <returns>True for a probe of a provider root or a generic connectivity check.</returns>
            internal static bool IsBlockedPreflight(string url)
            {
                string normalized = url.TrimEnd('/');
                return url.Contains("google", StringComparison.OrdinalIgnoreCase) ||
                       url.Contains("generate_204", StringComparison.OrdinalIgnoreCase) ||
                       url.StartsWith("https://waterservices.usgs.gov", StringComparison.OrdinalIgnoreCase) ||
                       string.Equals(normalized, "https://api.waterdata.usgs.gov", StringComparison.OrdinalIgnoreCase);
            }

            /// <inheritdoc/>
            protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
            {
                string url = request.RequestUri?.AbsoluteUri ?? "";
                lock (_gate)
                {
                    Requests.Add(url);
                    ApiKeys.Add(request.Headers.TryGetValues("X-Api-Key", out var keys) ? keys.FirstOrDefault() : null);
                }

                if (IsBlockedPreflight(url))
                    return Task.FromException<HttpResponseMessage>(new HttpRequestException($"Unexpected connectivity preflight request: {url}"));

                HttpResponseMessage response;
                try
                {
                    response = _route(url) ?? new HttpResponseMessage(HttpStatusCode.NotFound) { Content = new StringContent($"Unexpected URL: {url}") };
                }
                catch (Exception ex)
                {
                    return Task.FromException<HttpResponseMessage>(ex);
                }
                return _completeAsynchronously ? Task.Run(() => response) : Task.FromResult(response);
            }
        }

        /// <summary>
        /// Synchronization context that records posted continuations and executes them on the thread pool.
        /// </summary>
        /// <remarks>
        /// Tests use this context to prove downloader internals do not marshal long-running parsing
        /// continuations back to a WPF-like caller context.
        /// </remarks>
        private sealed class CountingSynchronizationContext : SynchronizationContext
        {
            /// <summary>
            /// Count of posted callbacks.
            /// </summary>
            private int _postCount;

            /// <summary>
            /// Gets the number of callbacks posted to this context.
            /// </summary>
            internal int PostCount => _postCount;

            /// <summary>
            /// Records and dispatches an asynchronous callback.
            /// </summary>
            /// <param name="d">The callback to execute.</param>
            /// <param name="state">The callback state.</param>
            /// <remarks>
            /// The callback is still executed so tests do not deadlock if a regression posts to the
            /// context; the post count then fails the assertion.
            /// </remarks>
            public override void Post(SendOrPostCallback d, object state)
            {
                Interlocked.Increment(ref _postCount);
                ThreadPool.QueueUserWorkItem(_ => d(state));
            }
        }

        /// <summary>
        /// HTTP handler that waits until the supplied cancellation token is canceled.
        /// </summary>
        /// <remarks>
        /// This handler verifies that downloader cancellation reaches the HTTP request rather than
        /// only canceling caller-side waiting.
        /// </remarks>
        private sealed class CancellableUsgsPeakHandler : HttpMessageHandler
        {
            /// <summary>
            /// Gets the absolute URLs requested through this handler.
            /// </summary>
            /// <remarks>
            /// The test uses this collection to confirm the request was started before cancellation.
            /// </remarks>
            internal List<string> Requests { get; } = new List<string>();

            /// <inheritdoc/>
            /// <remarks>
            /// The handler deliberately waits for cancellation and never returns a successful response.
            /// </remarks>
            protected override async Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
            {
                Requests.Add(request.RequestUri?.AbsoluteUri ?? "");
                await Task.Delay(TimeSpan.FromSeconds(30), cancellationToken);
                return new HttpResponseMessage(HttpStatusCode.OK)
                {
                    Content = new StringContent("")
                };
            }
        }

        /// <summary>
        /// HTTP handler that serves deterministic GHCN daily station-file responses.
        /// </summary>
        /// <remarks>
        /// The handler allows only the expected NOAA NCEI <c>USC00040741.dly</c> station-file URL
        /// and rejects generic connectivity probes. This keeps GHCN regression tests deterministic.
        /// </remarks>
        private sealed class GuardedGhcnDailyHandler : HttpMessageHandler
        {
            /// <summary>
            /// Whether the handler should simulate an HTTP timeout after URL validation.
            /// </summary>
            private readonly bool _simulateTimeout;

            /// <summary>
            /// Whether the handler should include a duplicate precipitation row for the same station month.
            /// </summary>
            private readonly bool _includeDuplicatePrecipitationRow;

            /// <summary>
            /// Gets the absolute URLs requested through this handler.
            /// </summary>
            internal List<string> Requests { get; } = new List<string>();

            /// <summary>
            /// Initializes a new instance of the <see cref="GuardedGhcnDailyHandler"/> class.
            /// </summary>
            /// <param name="simulateTimeout">Whether to throw a timeout after validating the request URL.</param>
            /// <param name="includeDuplicatePrecipitationRow">Whether to include duplicate daily precipitation records for the same dates.</param>
            /// <remarks>
            /// Timeout simulation is used to verify that the downloader reports the GHCN-specific
            /// timeout window without waiting for the full timeout duration. Duplicate-row simulation
            /// verifies inline last-value-wins normalization without a post-download pass.
            /// </remarks>
            internal GuardedGhcnDailyHandler(
                bool simulateTimeout = false,
                bool includeDuplicatePrecipitationRow = false)
            {
                _simulateTimeout = simulateTimeout;
                _includeDuplicatePrecipitationRow = includeDuplicatePrecipitationRow;
            }

            /// <inheritdoc/>
            /// <remarks>
            /// Returns a one-month GHCN daily file containing precipitation and snow rows only for
            /// the expected station-file URL.
            /// </remarks>
            protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
            {
                string url = request.RequestUri?.AbsoluteUri ?? "";
                Requests.Add(url);

                if (IsBlockedPreflight(url))
                {
                    return Task.FromException<HttpResponseMessage>(
                        new HttpRequestException($"Unexpected connectivity preflight request: {url}"));
                }

                if (!IsExpectedGhcnStationUrl(url))
                {
                    return Task.FromResult(new HttpResponseMessage(HttpStatusCode.NotFound)
                    {
                        Content = new StringContent($"Unexpected URL: {url}")
                    });
                }

                if (_simulateTimeout)
                {
                    return Task.FromException<HttpResponseMessage>(
                        new TaskCanceledException($"Simulated GHCN timeout: {url}"));
                }

                var rows = new List<string>
                {
                    BuildGhcnDailyLine("USC00040741", 2024, 1, "PRCP", new Dictionary<int, int>
                    {
                        [1] = 254,
                        [2] = -9999,
                        [3] = 508
                    }),
                    BuildGhcnDailyLine("USC00040741", 2024, 1, "SNOW", new Dictionary<int, int>
                    {
                        [1] = 0
                    })
                };

                if (_includeDuplicatePrecipitationRow)
                {
                    rows.Add(BuildGhcnDailyLine("USC00040741", 2024, 1, "PRCP", new Dictionary<int, int>
                    {
                        [1] = 762,
                        [2] = 1016,
                        [3] = 1270
                    }));
                }

                string body = string.Join(Environment.NewLine, rows);

                return Task.FromResult(new HttpResponseMessage(HttpStatusCode.OK)
                {
                    Content = new StringContent(body)
                });
            }

            /// <summary>
            /// Determines whether a URL is one of the connectivity probes this regression forbids.
            /// </summary>
            /// <param name="url">The absolute URL requested through the handler.</param>
            /// <returns><c>true</c> when the URL is a forbidden preflight; otherwise, <c>false</c>.</returns>
            /// <remarks>
            /// GHCN daily data should request the station file directly; generic root checks are not
            /// required before attempting the provider request.
            /// </remarks>
            internal static bool IsBlockedPreflight(string url)
            {
                string normalized = url.TrimEnd('/');
                return url.Contains("google", StringComparison.OrdinalIgnoreCase) ||
                       url.Contains("generate_204", StringComparison.OrdinalIgnoreCase) ||
                       string.Equals(normalized, "https://www.ncei.noaa.gov", StringComparison.OrdinalIgnoreCase) ||
                       string.Equals(normalized, "https://www.ncei.noaa.gov/pub/data/ghcn/daily/all", StringComparison.OrdinalIgnoreCase);
            }

            /// <summary>
            /// Determines whether a URL matches the expected GHCN station-file request.
            /// </summary>
            /// <param name="url">The absolute URL requested through the handler.</param>
            /// <returns><c>true</c> when the URL targets the expected fixture endpoint; otherwise, <c>false</c>.</returns>
            /// <remarks>
            /// The check pins the request to the NOAA NCEI daily station-file path and the exact
            /// station that failed in BestFit.
            /// </remarks>
            private static bool IsExpectedGhcnStationUrl(string url)
            {
                return string.Equals(
                    url,
                    "https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/USC00040741.dly",
                    StringComparison.OrdinalIgnoreCase);
            }

            /// <summary>
            /// Builds one fixed-width GHCN daily station-file line.
            /// </summary>
            /// <param name="stationId">The 11-character GHCN station identifier.</param>
            /// <param name="year">The four-digit year.</param>
            /// <param name="month">The one-based month number.</param>
            /// <param name="element">The four-character GHCN element code.</param>
            /// <param name="valuesByDay">Raw GHCN values keyed by one-based day of month.</param>
            /// <returns>A fixed-width daily line in GHCN <c>.dly</c> format.</returns>
            /// <remarks>
            /// Each daily value is represented by a five-character integer followed by three flag
            /// characters. Missing days are filled with <c>-9999</c>, matching the provider format.
            /// </remarks>
            private static string BuildGhcnDailyLine(
                string stationId,
                int year,
                int month,
                string element,
                IReadOnlyDictionary<int, int> valuesByDay)
            {
                var builder = new System.Text.StringBuilder();
                builder.Append(stationId);
                builder.Append(year.ToString("0000", CultureInfo.InvariantCulture));
                builder.Append(month.ToString("00", CultureInfo.InvariantCulture));
                builder.Append(element);

                for (int day = 1; day <= 31; day++)
                {
                    int value = valuesByDay.TryGetValue(day, out int parsedValue) ? parsedValue : -9999;
                    builder.Append(value.ToString(CultureInfo.InvariantCulture).PadLeft(5));
                    builder.Append("   ");
                }

                return builder.ToString();
            }
        }

        /// <summary>
        /// Probes a small windowed CHMN download to confirm Water Survey of Canada services are responsive.
        /// </summary>
        private static Task<bool> ChmnAvailable() => AvailableAsync("CHMN", () =>
            TimeSeriesDownload.FromCHMN(CHMN_1,
                TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                startDate: WinStart, endDate: WinStart.AddDays(1)));

        /// <summary>
        /// Probes a small USGS peak-discharge download (annual maxes, a small payload) to confirm USGS is responsive.
        /// </summary>
        private static Task<bool> UsgsAvailable() => AvailableAsync("USGS", async () =>
        {
            await TimeSeriesDownload.FromUSGS(USGS_3, TimeSeriesDownload.TimeSeriesType.PeakDischarge);
        });

        /// <summary>
        /// Probes a GHCN station file fetch to confirm NOAA NCEI services are responsive.
        /// </summary>
        private static Task<bool> GhcnAvailable() => AvailableAsync("GHCN", () =>
            TimeSeriesDownload.FromGHCN(GHCN_1));

        /// <summary>
        /// Probes a small windowed BOM download to confirm the KiWIS values endpoint (and the
        /// underlying WDP backend) is responsive. False when BOM returns 500 / DatasourceError.
        /// </summary>
        private static Task<bool> BomAvailable() => AvailableAsync("BOM", () =>
            TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                startDate: WinStart, endDate: WinStart.AddDays(1)));

        /// <summary>
        /// Verifies that a time series has monotonically increasing, non-duplicated date indices.
        /// </summary>
        /// <param name="ts">The time series to validate.</param>
        private static void AssertDailySeriesMonotonic(TimeSeries ts)
        {
            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);

            DateTime? prev = null;
            foreach (var pt in ts)
            {
                if (prev.HasValue)
                    Assert.IsGreaterThanOrEqualTo(prev.Value, pt.Index, "Dates not sorted chronologically.");
                prev = pt.Index;
            }

            // Ensure all indices are unique
            Assert.AreEqual(ts.Count, ts.Select(o => o.Index).Distinct().Count(), "Duplicate date indices detected.");
        }

        /// <summary>
        /// Compares two numeric values within specified relative and absolute tolerances.
        /// </summary>
        /// <param name="a">First value.</param>
        /// <param name="b">Second value.</param>
        /// <param name="relTol">Relative tolerance (default = 1e-6).</param>
        /// <param name="absTol">Absolute tolerance (default = 1e-9).</param>
        private static void AssertRoughlyEqual(double a, double b, double relTol = 1e-6, double absTol = 1e-9)
        {
            if (double.IsNaN(a) && double.IsNaN(b)) return;
            double diff = Math.Abs(a - b);
            if (diff <= absTol) return;

            double denom = Math.Max(Math.Abs(a), Math.Abs(b));
            if (denom == 0)
                Assert.IsLessThanOrEqualTo(absTol,diff);
            else
                Assert.IsLessThanOrEqualTo(relTol, diff / denom);
        }

        #endregion

        #region Connection Handling Tests

        /// <summary>
        /// Verifies connect-address ordering puts IPv4 before IPv6 while preserving relative order.
        /// </summary>
        /// <remarks>
        /// Providers such as NOAA NCEI publish several AAAA records. On networks that silently
        /// drop IPv6, connecting in DNS order costs a full TCP timeout per dead address, so the
        /// downloader must try IPv4 first.
        /// </remarks>
        [TestMethod]
        public void OrderAddressesForConnect_PutsIPv4First_PreservingRelativeOrder()
        {
            var v6First = IPAddress.Parse("2610:20:8040:2::178");
            var v4First = IPAddress.Parse("205.167.25.177");
            var v6Second = IPAddress.Parse("2610:20:8040:2::177");
            var v4Second = IPAddress.Parse("205.167.25.178");

            var ordered = TimeSeriesDownload.OrderAddressesForConnect(
                new[] { v6First, v4First, v6Second, v4Second });

            CollectionAssert.AreEqual(new[] { v4First, v4Second, v6First, v6Second }, ordered);
        }

        /// <summary>
        /// Verifies connect-address ordering passes an all-IPv6 result through unchanged.
        /// </summary>
        /// <remarks>
        /// IPv6-only hosts must remain reachable; preferring IPv4 must never drop IPv6 addresses.
        /// </remarks>
        [TestMethod]
        public void OrderAddressesForConnect_AllIPv6_PassesThroughUnchanged()
        {
            var addresses = new[]
            {
                IPAddress.Parse("2610:20:8040:2::178"),
                IPAddress.Parse("2610:20:8040:2::177")
            };

            var ordered = TimeSeriesDownload.OrderAddressesForConnect(addresses);

            CollectionAssert.AreEqual(addresses, ordered);
        }

        /// <summary>
        /// Verifies connect-address ordering tolerates an empty DNS result.
        /// </summary>
        [TestMethod]
        public void OrderAddressesForConnect_Empty_ReturnsEmpty()
        {
            Assert.IsEmpty(TimeSeriesDownload.OrderAddressesForConnect(Array.Empty<IPAddress>()));
        }

        /// <summary>
        /// Verifies the .NET Framework bind delegate rejects IPv6 remote endpoints.
        /// </summary>
        /// <remarks>
        /// On .NET Framework the service point tries resolved addresses in DNS order with no
        /// per-address bound. Throwing a socket exception from the bind delegate makes a blocked
        /// IPv6 attempt fail immediately so the service point moves on to an IPv4 address.
        /// </remarks>
        [TestMethod]
        public void RejectIPv6BindEndPoint_ThrowsSocketExceptionForIPv6()
        {
            var remote = new IPEndPoint(IPAddress.Parse("2610:20:8040:2::178"), 443);
            Assert.Throws<SocketException>(() => TimeSeriesDownload.RejectIPv6BindEndPoint(remote));
        }

        /// <summary>
        /// Verifies the .NET Framework bind delegate lets IPv4 remote endpoints bind normally.
        /// </summary>
        [TestMethod]
        public void RejectIPv6BindEndPoint_ReturnsNullForIPv4()
        {
            var remote = new IPEndPoint(IPAddress.Parse("205.167.25.177"), 443);
            Assert.IsNull(TimeSeriesDownload.RejectIPv6BindEndPoint(remote));
        }

        /// <summary>
        /// The User-Agent names the library and its major and minor version.
        /// </summary>
        [TestMethod]
        public void UserAgent_NamesLibraryVersion()
        {
            var version = typeof(TimeSeriesDownload).Assembly.GetName().Version!;
            Assert.AreEqual($"USACE-Numerics/{version.Major}.{version.Minor}", TimeSeriesDownload.UserAgent);
            Assert.IsTrue(System.Net.Http.Headers.ProductInfoHeaderValue.TryParse(TimeSeriesDownload.UserAgent, out _));
        }

        #endregion

        #region CHMN (Canada) Tests

        /// <summary>
        /// Validates a full-period-of-record download for the CHMN Cold River station (flow).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_FullPor_ColdRiver_Flow()
        {
            if (!await ChmnAvailable()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the CHMN Lillooet River station (flow).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_FullPor_Lillooet_Flow()
        {
            if (!await ChmnAvailable()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_2);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the CHMN Capilano River station (flow).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_FullPor_Capilano_Flow()
        {
            if (!await ChmnAvailable()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_3);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests flow unit conversions (cms ↔ cfs) for CHMN data.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_UnitConversion_Flow_CmsCfs()
        {
            if (!await ChmnAvailable()) return;

            var tsCms = await TimeSeriesDownload.FromCHMN(CHMN_1,
                TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                TimeSeriesDownload.DischargeUnit.CubicMetersPerSecond,
                startDate: WinStart, endDate: WinEnd);

            var tsCfs = await TimeSeriesDownload.FromCHMN(CHMN_1,
                TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                TimeSeriesDownload.DischargeUnit.CubicFeetPerSecond,
                startDate: WinStart, endDate: WinEnd);

            const double factor = 35.3146667;
            for (int i = 0; i < tsCms.Count; i++)
            {
                if (double.IsNaN(tsCms[i].Value) || double.IsNaN(tsCfs[i].Value)) continue;
                AssertRoughlyEqual(tsCfs[i].Value, tsCms[i].Value * factor);
            }
        }

        /// <summary>
        /// Tests stage unit conversions (m ↔ ft) for CHMN data.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_UnitConversion_Stage_MFt()
        {
            if (!await ChmnAvailable()) return;

            var tsM = await TimeSeriesDownload.FromCHMN(CHMN_1,
                TimeSeriesDownload.TimeSeriesType.DailyStage,
                heightUnit: TimeSeriesDownload.HeightUnit.Meters,
                startDate: WinStart, endDate: WinEnd);

            var tsFt = await TimeSeriesDownload.FromCHMN(CHMN_1,
                TimeSeriesDownload.TimeSeriesType.DailyStage,
                heightUnit: TimeSeriesDownload.HeightUnit.Feet,
                startDate: WinStart, endDate: WinEnd);

            const double factor = 3.280839895;
            for (int i = 0; i < tsM.Count; i++)
            {
                if (double.IsNaN(tsM[i].Value) || double.IsNaN(tsFt[i].Value)) continue;
                AssertRoughlyEqual(tsFt[i].Value, tsM[i].Value * factor);
            }
        }

        /// <summary>
        /// Ensures CHMN rejects invalid station identifiers.
        /// </summary>
        [TestMethod]
        public async Task CHMN_InvalidStation_Throws()
        {
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromCHMN("08LG01"));
        }

        /// <summary>
        /// Ensures CHMN rejects unsupported time series types (field measurements not available).
        /// </summary>
        [TestMethod]
        public async Task CHMN_UnsupportedType_Throws()
        {
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromCHMN(CHMN_1,
                    TimeSeriesDownload.TimeSeriesType.MeasuredDischarge));
        }

        /// <summary>
        /// Tests CHMN instantaneous discharge download (real-time 5-minute data).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_InstantaneousDischarge_Works()
        {
            if (!await ChmnAvailable()) return;

            var ts = await TimeSeriesDownload.FromCHMN(CHMN_2,
                TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        /// <summary>
        /// Tests CHMN instantaneous stage download (real-time 5-minute data).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_InstantaneousStage_Works()
        {
            if (!await ChmnAvailable()) return;

            var ts = await TimeSeriesDownload.FromCHMN(CHMN_2,
                TimeSeriesDownload.TimeSeriesType.InstantaneousStage);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        /// <summary>
        /// Tests CHMN peak discharge download (annual instantaneous maximums).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_PeakDischarge_Works()
        {
            if (!await ChmnAvailable()) return;

            var ts = await TimeSeriesDownload.FromCHMN(CHMN_1,
                TimeSeriesDownload.TimeSeriesType.PeakDischarge);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        /// <summary>
        /// Tests CHMN peak stage download (annual instantaneous maximums).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task CHMN_PeakStage_Works()
        {
            if (!await ChmnAvailable()) return;

            var ts = await TimeSeriesDownload.FromCHMN(CHMN_1,
                TimeSeriesDownload.TimeSeriesType.PeakStage);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        #endregion

        #region USGS Tests

        /// <summary>
        /// Runs a USGS download against a fake Water Data API.
        /// </summary>
        /// <param name="handler">The fake API.</param>
        /// <param name="download">The download to run.</param>
        /// <returns>The download result.</returns>
        private static async Task<(TimeSeries TimeSeries, string RawText)> WithFakeUsgs(FakeUsgsApiHandler handler, Func<Task<(TimeSeries TimeSeries, string RawText)>> download)
        {
            using var defaultClient = new HttpClient(handler) { Timeout = TimeSpan.FromSeconds(30) };
            using var decompressClient = new HttpClient(handler) { Timeout = TimeSpan.FromSeconds(30) };
            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                return await download();
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Single-series daily metadata response.
        /// </summary>
        private const string UsgsDailyMetadata = "x,y,id,begin,end,primary\n,,daily1,1930-03-01 05:00:00+00:00,2026-11-03 05:00:00+00:00,Primary\n";

        /// <summary>
        /// Daily values: CSV columns, a missing day, an empty value and a quoted qualifier list. No preflight request is made.
        /// </summary>
        [TestMethod]
        public async Task USGS_DailyDischarge_ParsesCsvAndFillsMissingDays()
        {
            var handler = new FakeUsgsApiHandler(url =>
            {
                if (url.Contains("/time-series-metadata/")) return FakeUsgsApiHandler.Text(UsgsDailyMetadata);
                if (url.Contains("/daily/items?")) return FakeUsgsApiHandler.Text(
                    "x,y,time,value,approval_status,qualifier\r\n" +
                    ",,2020-01-01,101.5,Approved,\r\n" +
                    ",,2020-01-02,,Approved,\"ICE,EQUIP\"\r\n" +
                    ",,2020-01-04,104,Provisional,\r\n");
                return null;
            });
            var (ts, raw) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6));

            Assert.AreEqual(TimeInterval.OneDay, ts.TimeInterval);
            Assert.AreEqual(4, ts.Count);
            Assert.AreEqual(new DateTime(2020, 1, 1), ts[0].Index);
            Assert.AreEqual(DateTimeKind.Unspecified, ts[0].Index.Kind);
            Assert.AreEqual(101.5, ts[0].Value);
            Assert.IsTrue(double.IsNaN(ts[1].Value), "An empty value is missing.");
            Assert.AreEqual(new DateTime(2020, 1, 3), ts[2].Index);
            Assert.IsTrue(double.IsNaN(ts[2].Value), "A missing day is filled with NaN.");
            Assert.AreEqual(104d, ts[3].Value);

            string data = handler.Requests.Single(u => u.Contains("/daily/items?"));
            StringAssert.StartsWith(data, "https://api.waterdata.usgs.gov/ogcapi/v1/collections/daily/items?");
            foreach (var part in new[] { "monitoring_location_id=USGS-01646500", "parameter_code=00060", "statistic_id=00003", "f=csv", "sortby=time", "limit=50000", "skipGeometry=true" })
                Assert.Contains(part, data);
            Assert.DoesNotContain("datetime=", data);
            Assert.IsFalse(handler.Requests.Any(FakeUsgsApiHandler.IsBlockedPreflight));
            StringAssert.StartsWith(raw, "# U.S. Geological Survey, Water Data API");
            Assert.Contains("\"ICE,EQUIP\"", raw);
        }

        /// <summary>
        /// Daily stage defaults to the daily mean; options choose the statistic and the date range.
        /// </summary>
        [TestMethod]
        public async Task USGS_DailyStage_OptionsSetStatisticAndDates()
        {
            var handler = new FakeUsgsApiHandler(url =>
                url.Contains("/time-series-metadata/") ? FakeUsgsApiHandler.Text(UsgsDailyMetadata) :
                url.Contains("/daily/items?") ? FakeUsgsApiHandler.Text("x,y,time,value\n,,2020-06-01,6.1\n") : null);

            await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.DailyStage));
            string meanUrl = handler.Requests.Single(u => u.Contains("/daily/items?"));
            Assert.Contains("parameter_code=00065", meanUrl);
            Assert.Contains("statistic_id=00003", meanUrl);

            handler.Requests.Clear();
            var options = new TimeSeriesDownload.UsgsDownloadOptions
            {
                Statistic = TimeSeriesDownload.UsgsStatistic.Maximum,
                StartDate = new DateTime(2020, 1, 1),
                EndDate = new DateTime(2020, 12, 31)
            };
            await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.DailyStage, options));
            string maxUrl = handler.Requests.Single(u => u.Contains("/daily/items?"));
            Assert.Contains("statistic_id=00001", maxUrl);
            Assert.Contains("datetime=2020-01-01/2020-12-31", maxUrl);

            await Assert.ThrowsAsync<ArgumentException>(async () => await TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.DailyStage,
                new TimeSeriesDownload.UsgsDownloadOptions { StartDate = new DateTime(2021, 1, 1), EndDate = new DateTime(2020, 1, 1) }));
            await Assert.ThrowsAsync<ArgumentNullException>(async () => await TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.DailyStage, null));
        }

        /// <summary>
        /// When a site has two daily series for the statistic, the values are requested again for the primary one.
        /// </summary>
        [TestMethod]
        public async Task USGS_Daily_MultipleSeries_UsesPrimary()
        {
            var handler = new FakeUsgsApiHandler(url =>
            {
                if (url.Contains("/time-series-metadata/")) return FakeUsgsApiHandler.Text(
                    "x,y,id,begin,end,primary\n,,aux,2001-01-01 05:00:00+00:00,2002-01-01 05:00:00+00:00,\n,,main,2001-01-01 05:00:00+00:00,2002-01-01 05:00:00+00:00,Primary\n");
                if (url.Contains("time_series_id=main")) return FakeUsgsApiHandler.Text("x,y,time,value\n,,2001-01-01,7\n,,2001-01-02,8\n");
                if (url.Contains("/daily/items?")) return FakeUsgsApiHandler.Text("x,y,time,value\n,,2001-01-01,7\n,,2001-01-01,99\n,,2001-01-02,8\n");
                return null;
            });
            var (ts, _) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6));
            Assert.AreEqual(2, ts.Count);
            Assert.AreEqual(7d, ts[0].Value);
            Assert.IsTrue(handler.Requests.Any(u => u.Contains("time_series_id=main")));
        }

        /// <summary>
        /// A full page of 50,000 records is followed by a query that starts the day after the last record.
        /// </summary>
        [TestMethod]
        public async Task USGS_Daily_FullPage_PagesForwardByDate()
        {
            var start = new DateTime(1880, 1, 1);
            var page = new System.Text.StringBuilder("x,y,time,value\n");
            for (int i = 0; i < 50000; i++) page.Append(",,").Append(start.AddDays(i).ToString("yyyy-MM-dd", CultureInfo.InvariantCulture)).Append(",1\n");
            string next = start.AddDays(50000).ToString("yyyy-MM-dd", CultureInfo.InvariantCulture);
            var handler = new FakeUsgsApiHandler(url =>
            {
                if (url.Contains("/time-series-metadata/")) return FakeUsgsApiHandler.Text(UsgsDailyMetadata);
                if (url.Contains("datetime=" + next + "/..")) return FakeUsgsApiHandler.Text($"x,y,time,value\n,,{next},2\n");
                if (url.Contains("/daily/items?")) return FakeUsgsApiHandler.Text(page.ToString());
                return null;
            });
            var (ts, raw) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6));
            Assert.AreEqual(50001, ts.Count);
            Assert.AreEqual(2d, ts[50000].Value);
            Assert.AreEqual(2, handler.Requests.Count(u => u.Contains("/daily/items?")));
            Assert.HasCount(1, System.Text.RegularExpressions.Regex.Matches(raw, "^x,y,time,value", System.Text.RegularExpressions.RegexOptions.Multiline), "The CSV header appears once.");
        }

        /// <summary>
        /// Peaks with a time of day are dated in the site's local standard time; missing month or day is taken as 1;
        /// stage peaks use parameter 00065.
        /// </summary>
        [TestMethod]
        public async Task USGS_Peaks_UseLocalStandardDate()
        {
            var handler = new FakeUsgsApiHandler(url =>
            {
                if (url.Contains("/monitoring-locations/items/USGS-01646500")) return FakeUsgsApiHandler.Text("{\"type\":\"Feature\",\"properties\":{\"time_zone_abbreviation\":\"EST\"}}");
                if (url.Contains("/peaks/items?")) return FakeUsgsApiHandler.Text(
                    "x,y,time,value,year,month,day,time_of_day,qualifier\n" +
                    ",,1889-01-01,484000,1889,,,,HISTORIC\n" +
                    ",,1931-04-03,37900,1931,4,3,,\n" +
                    ",,1942-10-01,,1942,10,,,\n" +
                    ",,1998-03-23,152000,1998,3,23,02:15:00+00:00,\n" +
                    ",,2006-06-28,103000,2006,6,28,15:00:00+00:00,\n");
                return null;
            });
            var (ts, raw) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.PeakDischarge));

            CollectionAssert.AreEqual(new[] { new DateTime(1889, 1, 1), new DateTime(1931, 4, 3), new DateTime(1998, 3, 22), new DateTime(2006, 6, 28) },
                ts.Select(o => o.Index).ToArray());
            CollectionAssert.AreEqual(new[] { 484000d, 37900d, 152000d, 103000d }, ts.Select(o => o.Value).ToArray());
            Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval, "Peaks fall on irregular dates.");
            Assert.Contains("parameter_code=00060", handler.Requests.Single(u => u.Contains("/peaks/")));
            Assert.Contains("HISTORIC", raw);

            handler.Requests.Clear();
            await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.PeakStage));
            Assert.Contains("parameter_code=00065", handler.Requests.Single(u => u.Contains("/peaks/")));
        }

        /// <summary>
        /// Verifies USGS downloads honor caller cancellation while waiting on HTTP.
        /// </summary>
        /// <returns>A task that completes when the cancellation regression check finishes.</returns>
        /// <remarks>
        /// BestFit relies on this behavior to stop waiting when an external provider is blocked or
        /// unavailable from the user's network.
        /// </remarks>
        [TestMethod]
        public async Task USGS_PeakDischarge_HonorsCancellationToken()
        {
            var defaultHandler = new CancellableUsgsPeakHandler();
            var decompressHandler = new CancellableUsgsPeakHandler();

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(60) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(60) };
            using var cts = new CancellationTokenSource(TimeSpan.FromMilliseconds(50));

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                await Assert.ThrowsAsync<OperationCanceledException>(async () =>
                    await TimeSeriesDownload.FromUSGS(USGS_3, TimeSeriesDownload.TimeSeriesType.PeakDischarge, cts.Token));

                Assert.IsTrue(decompressHandler.Requests.Any(u => u.Contains("/peaks/items?")), "Expected a started USGS peak request.");
                Assert.IsEmpty(defaultHandler.Requests);
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies duplicate USGS measured-stage timestamps keep the last source value.
        /// </summary>
        /// <remarks>
        /// The fixture is based on USGS 01570500, where the field-measurements API returns two approved gage-height
        /// records for 1936-03-15 05:00:00 with values 14.51 and 14.52.
        /// </remarks>
        [TestMethod]
        public async Task USGS_MeasuredStage_DuplicateDateTimes_LastValueWins()
        {
            var handler = new FakeUsgsApiHandler(url => url.Contains("/field-measurements/items?") ? FakeUsgsApiHandler.Text(
                "x,y,time,time_of_day,value,approval_status,qualifier\n" +
                ",,1935-01-01,17:30:00+00:00,9.5,Approved,\n" +
                ",,1935-01-02,17:30:00+00:00,,Approved,\n" +
                ",,1935-06-01,,7.25,Approved,\n" +
                ",,1936-03-15,05:00:00+00:00,14.51,Approved,\n" +
                ",,1936-03-15,05:00:00+00:00,14.52,Approved,\n") : null);
            var (ts, raw) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_5, TimeSeriesDownload.TimeSeriesType.MeasuredStage));

            var duplicateTime = new DateTime(1936, 3, 15, 5, 0, 0, DateTimeKind.Utc);
            Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
            Assert.AreEqual(3, ts.Count);
            Assert.AreEqual(new DateTime(1935, 1, 1, 17, 30, 0, DateTimeKind.Utc), ts[0].Index, "Sorted by time.");
            Assert.AreEqual(new DateTime(1935, 6, 1, 0, 0, 0, DateTimeKind.Utc), ts[1].Index, "A bare date is midnight UTC on every machine.");
            Assert.AreEqual(DateTimeKind.Utc, ts[1].Index.Kind);
            Assert.AreEqual(14.52d, ts.Single(o => o.Index == duplicateTime).Value);
            string url = handler.Requests.Single();
            StringAssert.StartsWith(url, "https://api.waterdata.usgs.gov/ogcapi/v1/collections/field-measurements/items?");
            Assert.Contains("monitoring_location_id=USGS-01570500", url);
            Assert.Contains("parameter_code=00065", url);
            Assert.Contains("f=csv", url);
            Assert.Contains("time_of_day", url);
            Assert.DoesNotContain("sortby", url);
            Assert.Contains("14.51", raw);
        }

        /// <summary>
        /// Builds a fake API for continuous values over 2023-11-15 to 2024-01-02, one window per year.
        /// </summary>
        /// <param name="parameterCode">The expected parameter code.</param>
        /// <param name="simulateTimeout">True to time out every continuous request.</param>
        /// <param name="completeAsynchronously">True to complete responses on the thread pool.</param>
        private static FakeUsgsApiHandler ContinuousFake(string parameterCode, bool simulateTimeout = false, bool completeAsynchronously = false)
        {
            return new FakeUsgsApiHandler(url =>
            {
                if (url.Contains("/time-series-metadata/") && url.Contains("parameter_code=" + parameterCode) && url.Contains("computation_period_identifier=Points"))
                    return FakeUsgsApiHandler.Text("x,y,id,begin,end,primary\n,,iv1,2023-11-15 00:00:00+00:00,2024-01-02 00:00:00+00:00,Primary\n");
                if (url.Contains("/continuous/items?time_series_id=iv1"))
                {
                    if (simulateTimeout) throw new TaskCanceledException("Simulated instantaneous timeout");
                    if (url.Contains("datetime=2023-11-15T00:00:00Z/2023-12-31T23:59:59Z"))
                        return FakeUsgsApiHandler.Text("x,y,time,value\n,,2023-11-15 00:00:00+00:00,101.5\n,,2023-11-15 00:15:00+00:00,102\n,,2023-11-15 00:45:00+00:00,104.5\n");
                    if (url.Contains("datetime=2024-01-01T00:00:00Z/2024-01-02T00:00:00Z"))
                        return FakeUsgsApiHandler.Text("x,y,time,value\n,,2024-01-01 00:00:00+00:00,110\n");
                }
                return null;
            }, completeAsynchronously);
        }

        /// <summary>
        /// Instantaneous discharge and stage: the primary series, one request per calendar year, irregular UTC dates in time order.
        /// </summary>
        [TestMethod]
        public async Task USGS_Instantaneous_YearlyChunks_UtcIrregular()
        {
            foreach (var (type, code) in new[] { (TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge, "00060"), (TimeSeriesDownload.TimeSeriesType.InstantaneousStage, "00065") })
            {
                var handler = ContinuousFake(code);
                var (ts, raw) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, type));
                Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
                CollectionAssert.AreEqual(new[]
                {
                    new DateTime(2023, 11, 15, 0, 0, 0, DateTimeKind.Utc), new DateTime(2023, 11, 15, 0, 15, 0, DateTimeKind.Utc),
                    new DateTime(2023, 11, 15, 0, 45, 0, DateTimeKind.Utc), new DateTime(2024, 1, 1, 0, 0, 0, DateTimeKind.Utc)
                }, ts.Select(o => o.Index).ToArray());
                Assert.IsTrue(ts.All(o => o.Index.Kind == DateTimeKind.Utc));
                Assert.AreEqual(110d, ts[3].Value);
                Assert.AreEqual(2, handler.Requests.Count(u => u.Contains("/continuous/items?")));
                Assert.Contains("time series iv1", raw);
                Assert.IsFalse(handler.Requests.Any(FakeUsgsApiHandler.IsBlockedPreflight));
            }
        }

        /// <summary>
        /// The date range in the options clips the yearly windows.
        /// </summary>
        [TestMethod]
        public async Task USGS_Instantaneous_OptionsClipWindows()
        {
            var handler = ContinuousFake("00060");
            var options = new TimeSeriesDownload.UsgsDownloadOptions { StartDate = new DateTime(2024, 1, 1, 0, 0, 0, DateTimeKind.Utc) };
            var (ts, _) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge, options));
            Assert.AreEqual(1, ts.Count);
            Assert.AreEqual(1, handler.Requests.Count(u => u.Contains("/continuous/items?")));
        }

        /// <summary>
        /// Verifies instantaneous USGS timeout messages use the longer instantaneous limit.
        /// </summary>
        [TestMethod]
        public async Task USGS_InstantaneousDischarge_TimeoutMessage_UsesInstantaneousLimit()
        {
            var handler = ContinuousFake("00060", simulateTimeout: true);
            var exception = await Assert.ThrowsAsync<TimeoutException>(async () =>
                await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge)));
            Assert.Contains("120 seconds", exception.Message);
        }

        /// <summary>
        /// Verifies USGS instantaneous parsing does not resume on the caller synchronization context.
        /// </summary>
        /// <remarks>
        /// BestFit starts downloads from the WPF UI thread. The fake responses complete asynchronously while a
        /// synchronization context is installed; the downloader must not post continuations back to it.
        /// </remarks>
        [TestMethod]
        public async Task USGS_InstantaneousDischarge_DoesNotPostBackToSynchronizationContext()
        {
            var handler = ContinuousFake("00060", completeAsynchronously: true);
            using var defaultClient = new HttpClient(handler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(handler) { Timeout = TimeSpan.FromSeconds(5) };
            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var context = new CountingSynchronizationContext();
                SynchronizationContext previousContext = SynchronizationContext.Current;
                Task<(TimeSeries TimeSeries, string RawText)> downloadTask;
                try
                {
                    SynchronizationContext.SetSynchronizationContext(context);
                    downloadTask = TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge);
                }
                finally
                {
                    SynchronizationContext.SetSynchronizationContext(previousContext);
                }

                var (ts, _) = await downloadTask;
                Assert.AreEqual(0, context.PostCount, "USGS instantaneous download should not post continuations to the caller synchronization context.");
                Assert.AreEqual(4, ts.Count);
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies duplicate-heavy USGS instantaneous responses are normalized quickly, last value winning.
        /// </summary>
        [TestMethod]
        [Timeout(15000, CooperativeCancellation = true)]
        public async Task USGS_InstantaneousDischarge_DuplicateHeavyResponse_NormalizesQuickly()
        {
            const int pairs = 20000;
            var body = new System.Text.StringBuilder("x,y,time,value\n");
            for (int i = 0; i < pairs; i++)
            {
                string time = new DateTime(2024, 1, 1).AddMinutes(i).ToString("yyyy-MM-dd HH:mm:ss", CultureInfo.InvariantCulture) + "+00:00";
                body.Append(",,").Append(time).Append(',').Append((1000 + i).ToString(CultureInfo.InvariantCulture)).Append('\n');
                body.Append(",,").Append(time).Append(',').Append((2000 + i).ToString(CultureInfo.InvariantCulture)).Append('\n');
            }
            var handler = new FakeUsgsApiHandler(url =>
                url.Contains("/time-series-metadata/") ? FakeUsgsApiHandler.Text("x,y,id,begin,end,primary\n,,iv1,2024-01-01 00:00:00+00:00,2024-02-01 00:00:00+00:00,Primary\n") :
                url.Contains("/continuous/items?") ? FakeUsgsApiHandler.Text(body.ToString()) : null);
            var (ts, _) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge));
            Assert.AreEqual(pairs, ts.Count);
            Assert.AreEqual(2000d + pairs - 1, ts[pairs - 1].Value);
        }

        /// <summary>
        /// The API key is sent in the X-Api-Key header of every request and never written to the raw text; the
        /// per-download key overrides the default key.
        /// </summary>
        [TestMethod]
        public async Task USGS_ApiKey_SentAsHeaderOnly()
        {
            var handler = new FakeUsgsApiHandler(url =>
                url.Contains("/time-series-metadata/") ? FakeUsgsApiHandler.Text(UsgsDailyMetadata) :
                url.Contains("/daily/items?") ? FakeUsgsApiHandler.Text("x,y,time,value\n,,2020-06-01,6.1\n") : null);
            string previous = TimeSeriesDownload.UsgsApiKey;
            try
            {
                TimeSeriesDownload.UsgsApiKey = "default-key";
                var (_, raw) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6));
                Assert.IsTrue(handler.ApiKeys.All(k => k == "default-key"));
                Assert.DoesNotContain("default-key", raw);
                Assert.IsFalse(handler.Requests.Any(u => u.Contains("api_key")));

                handler.ApiKeys.Clear();
                await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                    new TimeSeriesDownload.UsgsDownloadOptions { ApiKey = "my-key" }));
                Assert.IsTrue(handler.ApiKeys.All(k => k == "my-key"));
            }
            finally
            {
                TimeSeriesDownload.UsgsApiKey = previous;
            }
        }

        /// <summary>
        /// A rate-limit response with a long Retry-After fails at once, with guidance about the API key.
        /// </summary>
        [TestMethod]
        public async Task USGS_RateLimit_FailsFastWithGuidance()
        {
            var handler = new FakeUsgsApiHandler(url =>
            {
                var response = new HttpResponseMessage((HttpStatusCode)429) { Content = new StringContent("{\"error\":{\"code\":\"OVER_RATE_LIMIT\"}}") };
                response.Headers.RetryAfter = new System.Net.Http.Headers.RetryConditionHeaderValue(TimeSpan.FromSeconds(835));
                return response;
            });
            var exception = await Assert.ThrowsAsync<Exception>(async () =>
                await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.PeakDischarge)));
            Assert.Contains("rate limit", exception.Message);
            Assert.Contains("UsgsApiKey", exception.Message);
            Assert.AreEqual(1, handler.Requests.Count(u => u.Contains("/peaks/")), "No retries while the limit is spent.");
        }

        /// <summary>
        /// A site with no daily data reports that nothing matched.
        /// </summary>
        [TestMethod]
        public async Task USGS_Daily_NoData_Throws()
        {
            var handler = new FakeUsgsApiHandler(url =>
                url.Contains("/time-series-metadata/") ? FakeUsgsApiHandler.Text("x,y,id,begin,end,primary\n") :
                url.Contains("/daily/items?") ? FakeUsgsApiHandler.Text("x,y,time,value\n") : null);
            var exception = await Assert.ThrowsAsync<Exception>(async () => await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6)));
            Assert.Contains("No data found", exception.Message);
        }

        #region Resilience

        /// <summary>
        /// Response content whose body stream fails part way, like a connection that drops mid-body.
        /// </summary>
        private sealed class DroppedBodyContent : HttpContent
        {
            /// <summary>The bytes sent before the connection drops.</summary>
            private readonly byte[] _prefix;

            /// <summary>
            /// Initializes a new instance of the <see cref="DroppedBodyContent"/> class.
            /// </summary>
            /// <param name="prefix">The text sent before the connection drops.</param>
            internal DroppedBodyContent(string prefix) { _prefix = System.Text.Encoding.UTF8.GetBytes(prefix); }

            /// <inheritdoc/>
            protected override Task SerializeToStreamAsync(System.IO.Stream stream, TransportContext context) => throw new System.IO.IOException("The response ended prematurely.");

            /// <inheritdoc/>
            protected override Task<System.IO.Stream> CreateContentReadStreamAsync() => Task.FromResult<System.IO.Stream>(new DroppingStream(_prefix));

            /// <inheritdoc/>
            protected override bool TryComputeLength(out long length) { length = 0; return false; }
        }

        /// <summary>
        /// A read-only stream that returns a prefix and then throws <see cref="System.IO.IOException"/>.
        /// </summary>
        private sealed class DroppingStream : System.IO.MemoryStream
        {
            /// <summary>
            /// Initializes a new instance of the <see cref="DroppingStream"/> class.
            /// </summary>
            /// <param name="prefix">The bytes returned before the failure.</param>
            internal DroppingStream(byte[] prefix) : base(prefix) { }

            /// <inheritdoc/>
            public override int Read(byte[] buffer, int offset, int count)
            {
                int n = base.Read(buffer, offset, count);
                if (n == 0) throw new System.IO.IOException("The response ended prematurely.");
                return n;
            }

            /// <inheritdoc/>
            public override Task<int> ReadAsync(byte[] buffer, int offset, int count, CancellationToken cancellationToken) => Task.FromResult(Read(buffer, offset, count));

            /// <inheritdoc/>
            public override Task CopyToAsync(System.IO.Stream destination, int bufferSize, CancellationToken cancellationToken)
            {
                var buffer = new byte[bufferSize];
                int n;
                while ((n = Read(buffer, 0, buffer.Length)) > 0) destination.Write(buffer, 0, n);
                return Task.CompletedTask;
            }
        }

        /// <summary>
        /// A rate-limit response with a Retry-After delay.
        /// </summary>
        /// <param name="seconds">The Retry-After delay in seconds.</param>
        private static HttpResponseMessage RateLimited(int seconds)
        {
            var response = new HttpResponseMessage((HttpStatusCode)429) { Content = new StringContent("rate") };
            response.Headers.RetryAfter = new System.Net.Http.Headers.RetryConditionHeaderValue(TimeSpan.FromSeconds(seconds));
            return response;
        }

        /// <summary>
        /// A USGS page whose connection drops part way is retried (an IOException is transient), for every provider
        /// that uses the shared retry helper.
        /// </summary>
        [TestMethod]
        public async Task USGS_DroppedBody_IsRetried()
        {
            int calls = 0;
            var handler = new FakeUsgsApiHandler(url =>
            {
                if (url.Contains("/time-series-metadata/")) return FakeUsgsApiHandler.Text(UsgsDailyMetadata);
                if (!url.Contains("/daily/items?")) return null;
                if (Interlocked.Increment(ref calls) == 1) return new HttpResponseMessage(HttpStatusCode.OK) { Content = new DroppedBodyContent("x,y,time,value\n,,2020-01") };
                return FakeUsgsApiHandler.Text("x,y,time,value\n,,2020-01-01,1\n,,2020-01-02,2\n");
            });
            var (ts, _) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6));
            Assert.HasCount(2, ts);
            Assert.AreEqual(2, calls);
        }

        /// <summary>
        /// A USGS page cut off inside its last record (no final line break) is retried; one that stays cut off fails
        /// instead of returning partial data.
        /// </summary>
        [TestMethod]
        public async Task USGS_TruncatedPage_IsRetriedThenFails()
        {
            int calls = 0;
            var handler = new FakeUsgsApiHandler(url =>
            {
                if (url.Contains("/time-series-metadata/")) return FakeUsgsApiHandler.Text(UsgsDailyMetadata);
                if (!url.Contains("/daily/items?")) return null;
                return Interlocked.Increment(ref calls) == 1
                    ? FakeUsgsApiHandler.Text("x,y,time,value\n,,2020-01-01,1\n,,2020-01-0")
                    : FakeUsgsApiHandler.Text("x,y,time,value\n,,2020-01-01,1\n,,2020-01-02,2\n");
            });
            var (ts, _) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6));
            Assert.HasCount(2, ts);
            Assert.AreEqual(2, calls);

            calls = 0;
            var cut = new FakeUsgsApiHandler(url =>
            {
                if (url.Contains("/time-series-metadata/")) return FakeUsgsApiHandler.Text(UsgsDailyMetadata);
                if (!url.Contains("/daily/items?")) return null;
                Interlocked.Increment(ref calls);
                return FakeUsgsApiHandler.Text("x,y,time,value\n,,2020-01-01,1\n,,2020-01-0");
            });
            var e = await Assert.ThrowsAsync<System.IO.IOException>(async () => await WithFakeUsgs(cut, () => TimeSeriesDownload.FromUSGS(USGS_6)));
            Assert.Contains("ended before its last record", e.Message);
            Assert.AreEqual(4, calls, "Four attempts, then a clear failure.");
        }

        /// <summary>
        /// A timed-out request is retried at most once: two attempts instead of four.
        /// </summary>
        [TestMethod]
        public async Task USGS_Timeout_RetriedAtMostOnce()
        {
            var handler = ContinuousFake("00060", simulateTimeout: true);
            var options = new TimeSeriesDownload.UsgsDownloadOptions { StartDate = new DateTime(2024, 1, 1, 0, 0, 0, DateTimeKind.Utc) };
            await Assert.ThrowsAsync<TimeoutException>(async () =>
                await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge, options)));
            Assert.AreEqual(2, handler.Requests.Count(u => u.Contains("/continuous/items?")));
        }

        /// <summary>
        /// GHCN honours a short Retry-After (waits it, then succeeds) and stops at once on a long one, with a
        /// rate-limit message.
        /// </summary>
        [TestMethod]
        public async Task GHCN_RetryAfter_IsHonoured()
        {
            string dly = "USC00040741202401PRCP" + string.Concat(Enumerable.Repeat("   10  0", 31)) + "\n";
            int calls = 0;
            var shortWait = new FakeUsgsApiHandler(_ => Interlocked.Increment(ref calls) == 1 ? RateLimited(1) : FakeUsgsApiHandler.Text(dly));
            var watch = System.Diagnostics.Stopwatch.StartNew();
            var (ts, _) = await WithFakeUsgs(shortWait, async () => (await TimeSeriesDownload.FromGHCN(GHCN_1), ""));
            Assert.HasCount(31, ts);
            Assert.IsGreaterThanOrEqualTo(0.95, watch.Elapsed.TotalSeconds, "Waited the 1 s Retry-After, not the 0.5 s backoff.");

            var longWait = new FakeUsgsApiHandler(_ => RateLimited(3600));
            var e = await Assert.ThrowsAsync<Exception>(async () => await WithFakeUsgs(longWait, async () => (await TimeSeriesDownload.FromGHCN(GHCN_1), "")));
            Assert.Contains("rate limit", e.Message);
            Assert.HasCount(1, longWait.Requests, "A long Retry-After is final.");

            var notFound = new FakeUsgsApiHandler(_ => new HttpResponseMessage(HttpStatusCode.NotFound) { Content = new StringContent("missing") });
            await Assert.ThrowsAsync<HttpRequestException>(async () => await WithFakeUsgs(notFound, async () => (await TimeSeriesDownload.FromGHCN(GHCN_1), "")));
            Assert.HasCount(1, notFound.Requests, "A 404 is not retried.");
        }

        /// <summary>
        /// A GHCN station file whose connection drops part way is downloaded again.
        /// </summary>
        [TestMethod]
        public async Task GHCN_DroppedBody_IsRetried()
        {
            string dly = "USC00040741202401PRCP" + string.Concat(Enumerable.Repeat("   10  0", 31)) + "\n";
            int calls = 0;
            var handler = new FakeUsgsApiHandler(_ => Interlocked.Increment(ref calls) == 1
                ? new HttpResponseMessage(HttpStatusCode.OK) { Content = new DroppedBodyContent(dly.Substring(0, 100)) }
                : FakeUsgsApiHandler.Text(dly));
            var (ts, _) = await WithFakeUsgs(handler, async () => (await TimeSeriesDownload.FromGHCN(GHCN_1), ""));
            Assert.HasCount(31, ts);
            Assert.AreEqual(2, calls);
        }

        /// <summary>
        /// BOM requests go through the swappable client; the qualifier overload asks for the Quality Code field and a
        /// malformed values response gives a clear error.
        /// </summary>
        [TestMethod]
        public async Task BOM_TestClient_QualityCodesAndMalformedValues()
        {
            string list = "[[\"ts_id\",\"ts_name\"],[\"1573010\",\"DMQaQc.Merged.DailyMean.24HR\"]]";
            string values = "[{\"ts_id\":\"1573010\",\"columns\":\"Timestamp,Value,Quality Code\",\"data\":[[\"2021-11-01T00:00:00.000+10:00\",1.5,10],[\"2021-11-02T00:00:00.000+10:00\",2.5,90]]}]";
            var handler = new FakeUsgsApiHandler(url => FakeUsgsApiHandler.Text(url.Contains("getTimeseriesList") ? list : values));
            using var client = new HttpClient(handler);
            TimeSeriesDownload.SetHttpClientsForTesting(client, client, client);
            try
            {
                var plain = await TimeSeriesDownload.FromABOM(BOM_1, startDate: WinStart, endDate: WinEnd);
                Assert.HasCount(2, plain);
                Assert.DoesNotContain("returnfields", handler.Requests.Last());
                var q = new List<string>();
                var withQ = await TimeSeriesDownload.FromABOM(BOM_1, TimeSeriesDownload.TimeSeriesType.DailyDischarge, TimeSeriesDownload.DischargeUnit.CubicMetersPerSecond,
                    TimeSeriesDownload.HeightUnit.Meters, TimeSeriesDownload.DepthUnit.Millimeters, WinStart, WinEnd, q);
                Assert.Contains("returnfields=Timestamp%2CValue%2CQuality+Code", handler.Requests.Last());
                CollectionAssert.AreEqual(new[] { "10", "90" }, q);
                CollectionAssert.AreEqual(plain.Select(o => o.Value).ToArray(), withQ.Select(o => o.Value).ToArray());
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }

            var bad = new FakeUsgsApiHandler(url => FakeUsgsApiHandler.Text(url.Contains("getTimeseriesList") ? list : "[{\"data\":[[\"2021-11-01"));
            using var badClient = new HttpClient(bad);
            TimeSeriesDownload.SetHttpClientsForTesting(badClient, badClient, badClient);
            try
            {
                var e = await Assert.ThrowsAsync<Exception>(async () => await TimeSeriesDownload.FromABOM(BOM_1, startDate: WinStart, endDate: WinEnd));
                Assert.Contains("Failed to parse BOM time series values as JSON", e.Message);
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        #endregion

        /// <summary>
        /// The qualifier overloads return the same series plus one qualifier per ordinate; only instantaneous USGS
        /// requests change, by asking for the qualifier field.
        /// </summary>
        [TestMethod]
        public async Task Download_QualifierOverloads_SameSeries()
        {
            var handler = ContinuousFake("00060");
            var options = new TimeSeriesDownload.UsgsDownloadOptions();
            var (plain, _) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge, options));
            Assert.IsTrue(handler.Requests.Where(u => u.Contains("/continuous/")).All(u => u.Contains("properties=time,value&")));
            handler.Requests.Clear();
            var q = new List<string>();
            var (withQ, _) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge, options, q));
            Assert.IsTrue(handler.Requests.Where(u => u.Contains("/continuous/")).All(u => u.Contains("properties=time,value,qualifier&")));
            Assert.HasCount(withQ.Count, q);
            CollectionAssert.AreEqual(plain.Select(o => o.Value).ToArray(), withQ.Select(o => o.Value).ToArray());

            var daily = new FakeUsgsApiHandler(url =>
                url.Contains("/time-series-metadata/") ? FakeUsgsApiHandler.Text(UsgsDailyMetadata) :
                url.Contains("/daily/items?") ? FakeUsgsApiHandler.Text("x,y,time,value,approval_status,qualifier\n,,2020-01-01,1.5,Approved,\n,,2020-01-03,,Provisional,\"ICE,EQUIP\"\n") : null);
            var dq = new List<string>();
            var (dts, _) = await WithFakeUsgs(daily, () => TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.DailyDischarge, options, dq));
            CollectionAssert.AreEqual(new[] { "", "", "ICE,EQUIP" }, dq);
            Assert.HasCount(3, dts);

            var ghcn = new FakeUsgsApiHandler(_ => FakeUsgsApiHandler.Text("USC00040741202401PRCP" + "   10T 0" + string.Concat(Enumerable.Repeat("    0  0", 30)) + "\n"));
            var gq = new List<string>();
            var (gts, _) = await WithFakeUsgs(ghcn, async () => (await TimeSeriesDownload.FromGHCN(GHCN_1, TimeSeriesDownload.TimeSeriesType.DailyPrecipitation, TimeSeriesDownload.DepthUnit.Millimeters, gq), ""));
            Assert.AreEqual(1d, gts[0].Value);
            Assert.AreEqual("T 0", gq[0]);

            var chmn = new FakeUsgsApiHandler(_ => FakeUsgsApiHandler.Text("﻿ ID,Date,Parameter/Paramètre,Value/Valeur,Symbol/Symbole\n08LG010,2021-11-01,discharge/débit,1.5,B\n"));
            var cq = new List<string>();
            await WithFakeUsgs(chmn, async () => (await TimeSeriesDownload.FromCHMN(CHMN_1, TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                TimeSeriesDownload.DischargeUnit.CubicMetersPerSecond, TimeSeriesDownload.HeightUnit.Meters, null, null, cq), ""));
            CollectionAssert.AreEqual(new[] { "B" }, cq);

            await Assert.ThrowsAsync<ArgumentNullException>(async () => await TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.DailyDischarge, options, null!));
            await Assert.ThrowsAsync<ArgumentNullException>(async () => await TimeSeriesDownload.FromGHCN(GHCN_1, TimeSeriesDownload.TimeSeriesType.DailySnow, TimeSeriesDownload.DepthUnit.Inches, null!));
        }

        /// <summary>
        /// GHCN and CHMN report an empty, header-only or HTML response with status 200 as "No data found", like USGS.
        /// </summary>
        [TestMethod]
        public async Task GHCN_CHMN_EmptyOrHtmlResponse_ThrowsNoData()
        {
            const string chmnHeader = "﻿ ID,Date,Parameter/Paramètre,Value/Valeur,Symbol/Symbole\n";
            foreach (string body in new[] { "", "  \r\n", "<html><body>Access denied</body></html>" })
            {
                var ghcn = new FakeUsgsApiHandler(_ => FakeUsgsApiHandler.Text(body));
                var e = await Assert.ThrowsAsync<Exception>(async () => await WithFakeUsgs(ghcn, async () => (await TimeSeriesDownload.FromGHCN(GHCN_1), "")));
                Assert.Contains("No data found", e.Message);
            }
            foreach (string body in new[] { "", chmnHeader, chmnHeader + "\n", "<!DOCTYPE html><html></html>" })
            {
                var chmn = new FakeUsgsApiHandler(_ => FakeUsgsApiHandler.Text(body));
                var e = await Assert.ThrowsAsync<Exception>(async () => await WithFakeUsgs(chmn, async () => (await TimeSeriesDownload.FromCHMN("99ZZ999"), "")));
                Assert.Contains("No data found", e.Message);
            }

            // A station file without the requested element returns an empty series.
            var noSnow = new FakeUsgsApiHandler(_ => FakeUsgsApiHandler.Text("USC00040741202401PRCP" + string.Concat(Enumerable.Repeat("    0   ", 31)) + "\n"));
            var (snow, _) = await WithFakeUsgs(noSnow, async () => (await TimeSeriesDownload.FromGHCN(GHCN_1, TimeSeriesDownload.TimeSeriesType.DailySnow), ""));
            Assert.IsEmpty(snow);
            var data = new FakeUsgsApiHandler(_ => FakeUsgsApiHandler.Text(chmnHeader + "08LG010,2021-11-01,discharge/débit,1.5,\n"));
            var (one, _) = await WithFakeUsgs(data, async () => (await TimeSeriesDownload.FromCHMN(CHMN_1), ""));
            Assert.HasCount(1, one);
        }

        /// <summary>
        /// ParseUSGS applied to the raw text of a download gives the downloaded series, for daily, instantaneous, peak and
        /// measured types.
        /// </summary>
        [TestMethod]
        public async Task USGS_ParseRawText_ReproducesDownload()
        {
            var dailyHandler = new FakeUsgsApiHandler(url =>
                url.Contains("/time-series-metadata/") ? FakeUsgsApiHandler.Text(UsgsDailyMetadata) :
                url.Contains("/daily/items?") ? FakeUsgsApiHandler.Text("x,y,time,value,approval_status,qualifier\n,,2020-01-01,1.5,Approved,\n,,2020-01-03,,Approved,\"ICE,EQUIP\"\n") : null);
            var peakHandler = new FakeUsgsApiHandler(url =>
                url.Contains("/monitoring-locations/") ? FakeUsgsApiHandler.Text("{\"properties\":{\"time_zone_abbreviation\":\"EST\"}}") :
                url.Contains("/peaks/items?") ? FakeUsgsApiHandler.Text("x,y,time,value,year,month,day,time_of_day,qualifier\n,,1998-03-23,152000,1998,3,23,02:15:00+00:00,\n,,1889-01-01,484000,1889,,,,\n") : null);
            var measuredHandler = new FakeUsgsApiHandler(url => url.Contains("/field-measurements/items?") ? FakeUsgsApiHandler.Text(
                "x,y,time,time_of_day,value\n,,1936-03-15,05:00:00+00:00,14.5\n,,1935-06-01,,7.25\n") : null);
            var cases = new (FakeUsgsApiHandler Handler, TimeSeriesDownload.TimeSeriesType Type, TimeSpan? Offset)[]
            {
                (dailyHandler, TimeSeriesDownload.TimeSeriesType.DailyDischarge, null),
                (ContinuousFake("00060"), TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge, null),
                (peakHandler, TimeSeriesDownload.TimeSeriesType.PeakDischarge, TimeSpan.FromHours(-5)),
                (measuredHandler, TimeSeriesDownload.TimeSeriesType.MeasuredStage, null)
            };
            foreach (var (handler, type, offset) in cases)
            {
                var (ts, raw) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS(USGS_6, type));
                var parsed = TimeSeriesDownload.ParseUSGS(raw, type, offset);
                Assert.AreEqual(ts.TimeInterval, parsed.TimeInterval, type.ToString());
                CollectionAssert.AreEqual(ts.Select(o => o.Index).ToArray(), parsed.Select(o => o.Index).ToArray(), type.ToString());
                CollectionAssert.AreEqual(ts.Select(o => o.Index.Kind).ToArray(), parsed.Select(o => o.Index.Kind).ToArray(), type.ToString());
                CollectionAssert.AreEqual(ts.Select(o => BitConverter.DoubleToInt64Bits(o.Value)).ToArray(),
                    parsed.Select(o => BitConverter.DoubleToInt64Bits(o.Value)).ToArray(), type.ToString());
            }
        }

        /// <summary>
        /// Tests full-period-of-record USGS daily discharge download.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task USGS_FullPor_DailyDischarge()
        {
            if (!await UsgsAvailable()) return;

            var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests USGS daily stage download for correctness and continuity.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task USGS_FullPor_DailyStage()
        {
            if (!await UsgsAvailable()) return;

            var (ts, _) = await TimeSeriesDownload.FromUSGS(USGS_4, TimeSeriesDownload.TimeSeriesType.DailyStage);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests USGS peak discharge data retrieval for non-daily data.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task USGS_PeakDischarge_Works()
        {
            if (!await UsgsAvailable()) return;

            var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_3,
                TimeSeriesDownload.TimeSeriesType.PeakDischarge);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsNotNull(raw, "Raw text is null.");
        }

        /// <summary>
        /// Ensures USGS rejects invalid station identifiers.
        /// </summary>
        [TestMethod]
        public async Task USGS_InvalidStation_Throws()
        {
            foreach (string bad in new[] { "1134500", "1234567890123456", "0113450A", "USGS-01134500", "" })
                await Assert.ThrowsAsync<ArgumentException>(async () => await TimeSeriesDownload.FromUSGS(bad), bad);
            await Assert.ThrowsAsync<ArgumentException>(async () => await TimeSeriesDownload.FromUSGS(null!));
        }

        /// <summary>
        /// A 15-digit latitude-longitude site number is accepted and requested as USGS-{site}.
        /// </summary>
        [TestMethod]
        public async Task USGS_FifteenDigitSite_IsRequested()
        {
            var handler = new FakeUsgsApiHandler(url =>
                url.Contains("/time-series-metadata/") ? FakeUsgsApiHandler.Text(UsgsDailyMetadata) :
                url.Contains("/daily/items?") ? FakeUsgsApiHandler.Text("x,y,time,value\n,,2026-06-01,6.1\n") : null);
            var (ts, _) = await WithFakeUsgs(handler, () => TimeSeriesDownload.FromUSGS("393109104464500"));
            Assert.HasCount(1, ts);
            Assert.IsTrue(handler.Requests.All(u => u.Contains("monitoring_location_id=USGS-393109104464500")));
        }

        /// <summary>
        /// Ensures USGS rejects unsupported precipitation and snow types.
        /// </summary>
        [TestMethod]
        public async Task USGS_UnsupportedType_Throws()
        {
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromUSGS(USGS_1,
                    TimeSeriesDownload.TimeSeriesType.DailyPrecipitation));
        }

        /// <summary>
        /// Tests USGS field measurement discharge data retrieval from the OGC API.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task USGS_MeasuredDischarge_Works()
        {
            if (!await UsgsAvailable()) return;

            var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_5,
                TimeSeriesDownload.TimeSeriesType.MeasuredDischarge);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsNotNull(raw, "Raw text is null.");
            Assert.IsGreaterThan(0, ts.Count);
            Assert.IsFalse(string.IsNullOrEmpty(raw), "Raw text should contain JSON response.");

            // Verify all values are positive (discharge must be > 0)
            foreach (var pt in ts)
            {
                Assert.IsGreaterThan(0, pt.Value, $"Discharge value {pt.Value} at {pt.Index} should be positive.");
            }
        }

        /// <summary>
        /// Tests USGS field measurement stage (gage height) data retrieval from the OGC API.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task USGS_MeasuredStage_Works()
        {
            if (!await UsgsAvailable()) return;

            var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_5,
                TimeSeriesDownload.TimeSeriesType.MeasuredStage);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsNotNull(raw, "Raw text is null.");
            Assert.IsGreaterThan(0, ts.Count);
            Assert.IsFalse(string.IsNullOrEmpty(raw), "Raw text should contain JSON response.");

            // Verify all values are positive (gage height must be > 0)
            foreach (var pt in ts)
            {
                Assert.IsGreaterThan(0, pt.Value, $"Gage height value {pt.Value} at {pt.Index} should be positive.");
            }
        }

        /// <summary>
        /// Tests USGS peak stage data retrieval.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task USGS_PeakStage_Works()
        {
            if (!await UsgsAvailable()) return;

            var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_3,
                TimeSeriesDownload.TimeSeriesType.PeakStage);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsNotNull(raw, "Raw text is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        /// <summary>
        /// Tests USGS instantaneous discharge download.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task USGS_InstantaneousDischarge_Works()
        {
            if (!await UsgsAvailable()) return;

            var (ts, _) = await TimeSeriesDownload.FromUSGS(USGS_6,
                TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        /// <summary>
        /// Tests USGS instantaneous stage download.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task USGS_InstantaneousStage_Works()
        {
            if (!await UsgsAvailable()) return;

            var (ts, _) = await TimeSeriesDownload.FromUSGS(USGS_6,
                TimeSeriesDownload.TimeSeriesType.InstantaneousStage);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        #endregion

        #region GHCN Tests

        /// <summary>
        /// Verifies deterministic GHCN daily precipitation downloads parse station files directly.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// The fixture pins the request to NOAA NCEI station <c>USC00040741</c>, confirms no
        /// connectivity preflight is made, and verifies GHCN tenths-of-millimeters are converted to
        /// inches.
        /// </remarks>
        [TestMethod]
        public async Task GHCN_DailyPrecipitation_DeterministicDownload_ParsesAndNoPreflight()
        {
            var defaultHandler = new GuardedGhcnDailyHandler();
            var decompressHandler = new GuardedGhcnDailyHandler();

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var ts = await TimeSeriesDownload.FromGHCN(
                    GHCN_1,
                    TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                    TimeSeriesDownload.DepthUnit.Inches);

                Assert.IsNotNull(ts, "Time series is null.");
                Assert.AreEqual(TimeInterval.OneDay, ts.TimeInterval);
                Assert.AreEqual(31, ts.Count);
                Assert.AreEqual(new DateTime(2024, 1, 1), ts.StartDate);
                Assert.AreEqual(new DateTime(2024, 1, 31), ts.EndDate);
                Assert.AreEqual(1d, ts[0].Value);
                Assert.IsTrue(double.IsNaN(ts[1].Value));
                Assert.AreEqual(2d, ts[2].Value);
                Assert.HasCount(1, defaultHandler.Requests, "Expected exactly one GHCN station-file request.");
                Assert.IsEmpty(decompressHandler.Requests, "GHCN should not use the compressed client.");
                Assert.IsFalse(defaultHandler.Requests.Any(GuardedGhcnDailyHandler.IsBlockedPreflight));
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies GHCN duplicate daily timestamps are replaced while station rows are parsed.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// This pins the intended single-pass behavior: the later GHCN row for a duplicate date wins
        /// at ordinate creation time, avoiding the previous post-download duplicate scan.
        /// </remarks>
        [TestMethod]
        public async Task GHCN_DailyPrecipitation_DuplicateDateTimes_LastValueWinsInline()
        {
            var defaultHandler = new GuardedGhcnDailyHandler(includeDuplicatePrecipitationRow: true);
            var decompressHandler = new GuardedGhcnDailyHandler(includeDuplicatePrecipitationRow: true);

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var ts = await TimeSeriesDownload.FromGHCN(
                    GHCN_1,
                    TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                    TimeSeriesDownload.DepthUnit.Inches);

                Assert.AreEqual(31, ts.Count);
                Assert.AreEqual(ts.Count, ts.Select(o => o.Index).Distinct().Count(), "Duplicate date indices detected.");
                Assert.AreEqual(3d, ts[0].Value);
                Assert.AreEqual(4d, ts[1].Value);
                Assert.AreEqual(5d, ts[2].Value);
                Assert.HasCount(1, defaultHandler.Requests, "Expected exactly one GHCN station-file request.");
                Assert.IsEmpty(decompressHandler.Requests, "GHCN should not use the compressed client.");
                Assert.IsFalse(defaultHandler.Requests.Any(GuardedGhcnDailyHandler.IsBlockedPreflight));
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies GHCN timeout messages use the dedicated GHCN provider limit.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// Full station files can be larger than other provider payloads, so GHCN uses a limit
        /// above the generic provider timeout. The handler throws a timeout immediately after URL
        /// validation so this test does not wait for the full provider timeout.
        /// </remarks>
        [TestMethod]
        public async Task GHCN_DailyPrecipitation_TimeoutMessage_UsesGhcnLimit()
        {
            var defaultHandler = new GuardedGhcnDailyHandler(simulateTimeout: true);
            var decompressHandler = new GuardedGhcnDailyHandler(simulateTimeout: true);

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var exception = await Assert.ThrowsAsync<TimeoutException>(async () =>
                    await TimeSeriesDownload.FromGHCN(
                        GHCN_1,
                        TimeSeriesDownload.TimeSeriesType.DailyPrecipitation));

                Assert.Contains("60 seconds", exception.Message);
                Assert.HasCount(2, defaultHandler.Requests, "Expected the retry helper to make two GHCN attempts.");
                Assert.IsEmpty(decompressHandler.Requests, "GHCN should not use the compressed client.");
                Assert.IsFalse(defaultHandler.Requests.Any(GuardedGhcnDailyHandler.IsBlockedPreflight));
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Tests full-period-of-record GHCN daily precipitation download.
        /// </summary>
        /// <remarks>
        /// The ceiling guards against connection-establishment regressions: walking dead IPv6
        /// addresses before IPv4 once made this download take minutes instead of seconds.
        /// </remarks>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        [Timeout(60000, CooperativeCancellation = true)]
        public async Task GHCN_FullPor_Precipitation()
        {
            if (!await GhcnAvailable()) return;
            var ts = await TimeSeriesDownload.FromGHCN(GHCN_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests full-period-of-record GHCN daily snow download.
        /// </summary>
        /// <remarks>
        /// The ceiling guards against connection-establishment regressions: walking dead IPv6
        /// addresses before IPv4 once made this download take minutes instead of seconds.
        /// </remarks>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        [Timeout(60000, CooperativeCancellation = true)]
        public async Task GHCN_FullPor_Snow()
        {
            if (!await GhcnAvailable()) return;
            var ts = await TimeSeriesDownload.FromGHCN(GHCN_3,
                TimeSeriesDownload.TimeSeriesType.DailySnow);
            if (ts.Count > 0) AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests precipitation unit conversion between millimeters and inches.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task GHCN_UnitConversion_Mm_In()
        {
            if (!await GhcnAvailable()) return;

            var tsMm = await TimeSeriesDownload.FromGHCN(GHCN_2,
                TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                TimeSeriesDownload.DepthUnit.Millimeters);

            var tsIn = await TimeSeriesDownload.FromGHCN(GHCN_2,
                TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                TimeSeriesDownload.DepthUnit.Inches);

            const double factor = 25.4;
            for (int i = 0; i < Math.Min(tsMm.Count, tsIn.Count); i++)
            {
                if (double.IsNaN(tsMm[i].Value) || double.IsNaN(tsIn[i].Value)) continue;
                AssertRoughlyEqual(tsIn[i].Value, tsMm[i].Value / factor);
            }
        }

        /// <summary>
        /// Tests precipitation unit conversion between millimeters and centimeters.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration")]
        public async Task GHCN_UnitConversion_Mm_Cm()
        {
            if (!await GhcnAvailable()) return;

            var tsMm = await TimeSeriesDownload.FromGHCN(GHCN_1,
                TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                TimeSeriesDownload.DepthUnit.Millimeters);

            var tsCm = await TimeSeriesDownload.FromGHCN(GHCN_1,
                TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                TimeSeriesDownload.DepthUnit.Centimeters);

            for (int i = 0; i < Math.Min(tsMm.Count, tsCm.Count); i++)
            {
                if (double.IsNaN(tsMm[i].Value) || double.IsNaN(tsCm[i].Value)) continue;
                AssertRoughlyEqual(tsCm[i].Value, tsMm[i].Value / 10.0);
            }
        }

        /// <summary>
        /// Ensures GHCN rejects invalid station identifiers.
        /// </summary>
        [TestMethod]
        public async Task GHCN_InvalidStation_Throws()
        {
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromGHCN("USC0004074"));
        }

        /// <summary>
        /// Ensures GHCN rejects unsupported stage-type requests.
        /// </summary>
        [TestMethod]
        public async Task GHCN_UnsupportedType_Throws()
        {
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromGHCN(GHCN_1,
                    TimeSeriesDownload.TimeSeriesType.DailyStage));
        }

        #endregion

        #region BOM (Australia) Tests

        // Live BOM checks remain available locally but are excluded from PR integration.
        // Run them with: dotnet test -c Release --filter "TestCategory=BOMIntegration"

        /// <summary>
        /// Validates a full-period-of-record download for the BOM Cotter River station (discharge).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_FullPor_CotterRiver_Discharge()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the BOM Goodradigbee River station (discharge).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_FullPor_Goodradigbee_Discharge()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_2);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the BOM Murray River station (stage).
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_FullPor_MurrayRiver_Stage()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_3, TimeSeriesDownload.TimeSeriesType.DailyStage);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests discharge unit conversions (cms ↔ cfs) for BOM data.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_UnitConversion_Discharge_CmsCfs()
        {
            if (!await BomAvailable()) return;

            var tsCms = await TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                TimeSeriesDownload.DischargeUnit.CubicMetersPerSecond,
                startDate: WinStart, endDate: WinEnd);

            var tsCfs = await TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                TimeSeriesDownload.DischargeUnit.CubicFeetPerSecond,
                startDate: WinStart, endDate: WinEnd);

            const double factor = 35.3146667;
            for (int i = 0; i < tsCms.Count; i++)
            {
                if (double.IsNaN(tsCms[i].Value) || double.IsNaN(tsCfs[i].Value)) continue;
                AssertRoughlyEqual(tsCfs[i].Value, tsCms[i].Value * factor);
            }
        }

        /// <summary>
        /// Tests stage unit conversions (m ↔ ft) for BOM data.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_UnitConversion_Stage_MFt()
        {
            if (!await BomAvailable()) return;

            var tsM = await TimeSeriesDownload.FromABOM(BOM_3,
                TimeSeriesDownload.TimeSeriesType.DailyStage,
                heightUnit: TimeSeriesDownload.HeightUnit.Meters,
                startDate: WinStart, endDate: WinEnd);

            var tsFt = await TimeSeriesDownload.FromABOM(BOM_3,
                TimeSeriesDownload.TimeSeriesType.DailyStage,
                heightUnit: TimeSeriesDownload.HeightUnit.Feet,
                startDate: WinStart, endDate: WinEnd);

            const double factor = 3.280839895;
            for (int i = 0; i < tsM.Count; i++)
            {
                if (double.IsNaN(tsM[i].Value) || double.IsNaN(tsFt[i].Value)) continue;
                AssertRoughlyEqual(tsFt[i].Value, tsM[i].Value * factor);
            }
        }

        /// <summary>
        /// Ensures BOM rejects invalid station identifiers.
        /// </summary>
        [TestMethod]
        public async Task BOM_InvalidStation_Throws()
        {
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromABOM("41073"));
        }

        /// <summary>
        /// Ensures BOM rejects unsupported time series types.
        /// </summary>
        [TestMethod]
        public async Task BOM_UnsupportedType_Throws()
        {
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromABOM(BOM_1,
                    TimeSeriesDownload.TimeSeriesType.PeakDischarge));
        }

        /// <summary>
        /// Tests BOM with a windowed date range to verify date filtering works correctly.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_WindowedDownload_Works()
        {
            if (!await BomAvailable()) return;

            var ts = await TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                startDate: WinStart, endDate: WinEnd);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);

            // Verify data is within requested window (allowing for some timezone flexibility)
            var firstDate = ts.First().Index;
            var lastDate = ts.Last().Index;

            Assert.IsGreaterThanOrEqualTo(WinStart.AddDays(-1), firstDate,
                $"First date {firstDate} is before window start {WinStart}");
            Assert.IsLessThanOrEqualTo(WinEnd.AddDays(1), lastDate,
                $"Last date {lastDate} is after window end {WinEnd}");
        }

        /// <summary>
        /// Validates instantaneous discharge download from BOM.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_InstantaneousDischarge_Works()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge,
                startDate: WinStart, endDate: WinEnd);
            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        /// <summary>
        /// Validates instantaneous stage download from BOM.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_InstantaneousStage_Works()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_3,
                TimeSeriesDownload.TimeSeriesType.InstantaneousStage,
                startDate: WinStart, endDate: WinEnd);
            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        /// <summary>
        /// Validates daily precipitation download from BOM.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_DailyPrecipitation_Works()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                startDate: WinStart, endDate: WinEnd);
            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);
        }

        /// <summary>
        /// Tests precipitation unit conversions (mm ↔ inches) for BOM data.
        /// </summary>
        [TestMethod, TestCategory("Integration"), TestCategory("LiveServiceIntegration"), TestCategory("BOMIntegration")]
        public async Task BOM_UnitConversion_Precip_MmIn()
        {
            if (!await BomAvailable()) return;

            var tsMm = await TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                depthUnit: TimeSeriesDownload.DepthUnit.Millimeters,
                startDate: WinStart, endDate: WinEnd);

            var tsIn = await TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.DailyPrecipitation,
                depthUnit: TimeSeriesDownload.DepthUnit.Inches,
                startDate: WinStart, endDate: WinEnd);

            const double factor = 25.4;
            for (int i = 0; i < tsMm.Count; i++)
            {
                if (double.IsNaN(tsMm[i].Value) || double.IsNaN(tsIn[i].Value)) continue;
                if (tsMm[i].Value == 0 && tsIn[i].Value == 0) continue;
                AssertRoughlyEqual(tsMm[i].Value, tsIn[i].Value * factor);
            }
        }

        #endregion
    }
}

