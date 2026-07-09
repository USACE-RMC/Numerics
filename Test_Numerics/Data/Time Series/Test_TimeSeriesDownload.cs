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
    /// the United States Geological Survey (USGS), and the Global Historical Climatology Network (GHCN).
    /// </summary>
    /// <remarks>
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
        /// How long a "service unreachable" verdict sticks before we re-probe.
        /// </summary>
        private static readonly TimeSpan ProbeFailureTtl = TimeSpan.FromSeconds(30);

        /// <summary>
        /// "Can we hit this service?" probe with TTL on negative results. Returns false
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
        /// HTTP handler that serves one deterministic USGS peak-flow response and rejects connectivity preflight requests.
        /// </summary>
        /// <remarks>
        /// The handler makes the regression test independent of live internet access while still
        /// verifying the downloader requests the expected USGS endpoint directly.
        /// </remarks>
        private sealed class GuardedUsgsPeakHandler : HttpMessageHandler
        {
            /// <summary>
            /// Gets the absolute URLs requested through this handler.
            /// </summary>
            /// <remarks>
            /// The test uses this collection to verify no hidden preflight request was issued.
            /// </remarks>
            internal List<string> Requests { get; } = new List<string>();

            /// <inheritdoc/>
            /// <remarks>
            /// Returns a small USGS peak-flow RDB body only for the expected data URL and fails
            /// connectivity preflight requests explicitly.
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

                if (!IsExpectedUsgsPeakUrl(url))
                {
                    return Task.FromResult(new HttpResponseMessage(HttpStatusCode.NotFound)
                    {
                        Content = new StringContent($"Unexpected URL: {url}")
                    });
                }

                string body = string.Join(Environment.NewLine, new[]
                {
                    "USGS\t01614000\t2020-04-21\t1\t12345"
                });

                return Task.FromResult(new HttpResponseMessage(HttpStatusCode.OK)
                {
                    Content = new StringContent(body)
                });
            }

            /// <summary>
            /// Determines whether a URL is one of the connectivity probes this regression test forbids.
            /// </summary>
            /// <param name="url">The absolute URL requested through the handler.</param>
            /// <returns>True when the URL is a forbidden connectivity preflight; otherwise false.</returns>
            /// <remarks>
            /// The blocked list covers the old Google probe and the older generic USGS root probe.
            /// </remarks>
            internal static bool IsBlockedPreflight(string url)
            {
                return url.Contains("google", StringComparison.OrdinalIgnoreCase) ||
                       url.Contains("generate_204", StringComparison.OrdinalIgnoreCase) ||
                       url.StartsWith("https://waterservices.usgs.gov/nwis", StringComparison.OrdinalIgnoreCase);
            }

            /// <summary>
            /// Determines whether a URL matches the expected USGS peak-flow request for this test fixture.
            /// </summary>
            /// <param name="url">The absolute URL requested through the handler.</param>
            /// <returns>True when the URL targets the expected USGS peak-flow endpoint; otherwise false.</returns>
            /// <remarks>
            /// The check keeps the test tied to the production URL builder, including site number,
            /// agency code, and RDB output format.
            /// </remarks>
            private static bool IsExpectedUsgsPeakUrl(string url)
            {
                return url.StartsWith("https://nwis.waterdata.usgs.gov/nwis/peak?", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("site_no=01614000", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("agency_cd=USGS", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("format=rdb", StringComparison.OrdinalIgnoreCase);
            }
        }

        /// <summary>
        /// HTTP handler that serves deterministic USGS instantaneous RDB responses.
        /// </summary>
        /// <remarks>
        /// The handler allows only the expected <c>/nwis/iv</c> request for USGS 01646500 and
        /// rejects generic connectivity probes. This keeps the test deterministic and independent
        /// from the live USGS service.
        /// </remarks>
        private sealed class GuardedUsgsInstantaneousHandler : HttpMessageHandler
        {
            /// <summary>
            /// USGS parameter code expected by this handler.
            /// </summary>
            private readonly string _parameterCode;

            /// <summary>
            /// Whether the handler should simulate an HTTP timeout after URL validation.
            /// </summary>
            private readonly bool _simulateTimeout;

            /// <summary>
            /// Whether successful responses should complete asynchronously.
            /// </summary>
            private readonly bool _completeAsynchronously;

            /// <summary>
            /// Number of duplicate timestamp pairs to append to the deterministic response.
            /// </summary>
            private readonly int _duplicateTimestampPairs;

            /// <summary>
            /// Gets the absolute URLs requested through this handler.
            /// </summary>
            internal List<string> Requests { get; } = new List<string>();

            /// <summary>
            /// Initializes a new instance of the <see cref="GuardedUsgsInstantaneousHandler"/> class.
            /// </summary>
            /// <param name="parameterCode">The expected USGS parameter code.</param>
            /// <param name="simulateTimeout">Whether to throw a timeout after validating the request URL.</param>
            /// <param name="completeAsynchronously">Whether to complete successful responses on the thread pool.</param>
            /// <param name="duplicateTimestampPairs">The number of duplicate timestamp pairs to include in the response.</param>
            /// <exception cref="ArgumentException">Thrown when <paramref name="parameterCode"/> is empty.</exception>
            /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="duplicateTimestampPairs"/> is negative.</exception>
            /// <remarks>
            /// Parameter <c>00060</c> represents discharge and <c>00065</c> represents gage height.
            /// Asynchronous completion is used by synchronization-context regression tests, while
            /// duplicate timestamp pairs exercise the downloader's normalization cleanup.
            /// </remarks>
            internal GuardedUsgsInstantaneousHandler(
                string parameterCode,
                bool simulateTimeout = false,
                bool completeAsynchronously = false,
                int duplicateTimestampPairs = 0)
            {
                if (string.IsNullOrWhiteSpace(parameterCode))
                    throw new ArgumentException("A parameter code is required.", nameof(parameterCode));
                if (duplicateTimestampPairs < 0)
                    throw new ArgumentOutOfRangeException(nameof(duplicateTimestampPairs), "Duplicate timestamp pair count cannot be negative.");

                _parameterCode = parameterCode;
                _simulateTimeout = simulateTimeout;
                _completeAsynchronously = completeAsynchronously;
                _duplicateTimestampPairs = duplicateTimestampPairs;
            }

            /// <inheritdoc/>
            /// <remarks>
            /// Returns a small nonuniform instantaneous RDB body only for the expected USGS data URL.
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

                if (!IsExpectedUsgsInstantaneousUrl(url))
                {
                    return Task.FromResult(new HttpResponseMessage(HttpStatusCode.NotFound)
                    {
                        Content = new StringContent($"Unexpected URL: {url}")
                    });
                }

                if (_simulateTimeout)
                {
                    return Task.FromException<HttpResponseMessage>(
                        new TaskCanceledException($"Simulated instantaneous timeout: {url}"));
                }

                string valueOne = _parameterCode == "00060" ? "101.5" : "6.12";
                string valueTwo = _parameterCode == "00060" ? "102.0" : "6.13";
                string valueThree = _parameterCode == "00060" ? "104.5" : "6.17";
                var lines = new List<string>
                {
                    $"agency_cd\tsite_no\tdatetime\ttz_cd\t{_parameterCode}_00000\t{_parameterCode}_00000_cd",
                    "5s\t15s\t20d\t6s\t14n\t10s",
                    $"USGS\t01646500\t2024-01-01 00:00\tEST\t{valueOne}\tA",
                    $"USGS\t01646500\t2024-01-01 00:15\tEST\t{valueTwo}\tA",
                    $"USGS\t01646500\t2024-01-01 00:45\tEST\t{valueThree}\tA"
                };

                for (int i = 0; i < _duplicateTimestampPairs; i++)
                {
                    string timestamp = new DateTime(2024, 1, 2).AddMinutes(i).ToString("yyyy-MM-dd HH:mm", CultureInfo.InvariantCulture);
                    lines.Add($"USGS\t01646500\t{timestamp}\tEST\t{1000d + i:0.0}\tP");
                    lines.Add($"USGS\t01646500\t{timestamp}\tEST\t{2000d + i:0.0}\tA");
                }

                string body = string.Join(Environment.NewLine, lines);

                var response = new HttpResponseMessage(HttpStatusCode.OK)
                {
                    Content = new StringContent(body)
                };

                return _completeAsynchronously
                    ? Task.Run(() => response)
                    : Task.FromResult(response);
            }

            /// <summary>
            /// Determines whether a URL is one of the connectivity probes this regression forbids.
            /// </summary>
            /// <param name="url">The absolute URL requested through the handler.</param>
            /// <returns><c>true</c> when the URL is a forbidden preflight; otherwise, <c>false</c>.</returns>
            /// <remarks>
            /// The instantaneous data endpoint is under <c>/nwis/iv</c>; only generic roots and the
            /// old Google probe are treated as preflight requests.
            /// </remarks>
            internal static bool IsBlockedPreflight(string url)
            {
                string normalized = url.TrimEnd('/');
                return url.Contains("google", StringComparison.OrdinalIgnoreCase) ||
                       url.Contains("generate_204", StringComparison.OrdinalIgnoreCase) ||
                       string.Equals(normalized, "https://waterservices.usgs.gov/nwis", StringComparison.OrdinalIgnoreCase);
            }

            /// <summary>
            /// Determines whether a URL matches the expected USGS instantaneous request.
            /// </summary>
            /// <param name="url">The absolute URL requested through the handler.</param>
            /// <returns><c>true</c> when the URL targets the expected fixture endpoint; otherwise, <c>false</c>.</returns>
            /// <remarks>
            /// The check pins the request to the instantaneous-value API, the 01646500 gage, the
            /// requested parameter code, and RDB formatting.
            /// </remarks>
            private bool IsExpectedUsgsInstantaneousUrl(string url)
            {
                return url.StartsWith("https://waterservices.usgs.gov/nwis/iv/?", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("sites=01646500", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains($"parameterCd={_parameterCode}", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("startDT=1900-01-01", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("siteStatus=all", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("format=rdb", StringComparison.OrdinalIgnoreCase);
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
        /// HTTP handler that serves a deterministic USGS measured-stage response with duplicate timestamps.
        /// </summary>
        /// <remarks>
        /// The fixture is based on USGS 01570500, where the live field-measurements API returns two
        /// approved gage-height records for 1936-03-15 05:00:00 with values 14.51 and 14.52.
        /// </remarks>
        private sealed class DuplicateUsgsMeasuredStageHandler : HttpMessageHandler
        {
            /// <summary>
            /// Gets the absolute URLs requested through this handler.
            /// </summary>
            internal List<string> Requests { get; } = new List<string>();

            /// <inheritdoc/>
            protected override Task<HttpResponseMessage> SendAsync(HttpRequestMessage request, CancellationToken cancellationToken)
            {
                string url = request.RequestUri?.AbsoluteUri ?? "";
                Requests.Add(url);

                if (GuardedUsgsPeakHandler.IsBlockedPreflight(url))
                {
                    return Task.FromException<HttpResponseMessage>(
                        new HttpRequestException($"Unexpected connectivity preflight request: {url}"));
                }

                if (!IsExpectedUsgsMeasuredStageUrl(url))
                {
                    return Task.FromResult(new HttpResponseMessage(HttpStatusCode.NotFound)
                    {
                        Content = new StringContent($"Unexpected URL: {url}")
                    });
                }

                string body = string.Join("", new[]
                {
                    "{",
                    "\"type\":\"FeatureCollection\",",
                    "\"features\":[",
                    "{\"type\":\"Feature\",\"id\":\"75f11a56-5796-4dc9-b22a-054da446a164\",\"properties\":{\"time\":\"1936-03-15T05:00:00+00:00\",\"value\":\"14.51\"}},",
                    "{\"type\":\"Feature\",\"id\":\"8dccaa2d-8480-4b5d-b19b-18122749b6c6\",\"properties\":{\"time\":\"1936-03-15T05:00:00+00:00\",\"value\":\"14.52\"}}",
                    "],",
                    "\"links\":[]",
                    "}"
                });

                return Task.FromResult(new HttpResponseMessage(HttpStatusCode.OK)
                {
                    Content = new StringContent(body)
                });
            }

            /// <summary>
            /// Determines whether a URL matches the expected USGS measured-stage request.
            /// </summary>
            private static bool IsExpectedUsgsMeasuredStageUrl(string url)
            {
                return url.StartsWith("https://api.waterdata.usgs.gov/ogcapi/v0/collections/field-measurements/items?",
                           StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("monitoring_location_id=USGS-01570500", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("parameter_code=00065", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("limit=10000", StringComparison.OrdinalIgnoreCase) &&
                       url.Contains("f=json", StringComparison.OrdinalIgnoreCase);
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
        /// Probes a small USGS peak-discharge download (annual maxes — a small payload) to confirm USGS is responsive.
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

        #endregion

        #region CHMN (Canada) Tests

        /// <summary>
        /// Validates a full-period-of-record download for the CHMN Cold River station (flow).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task CHMN_FullPor_ColdRiver_Flow()
        {
            if (!await ChmnAvailable()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the CHMN Lillooet River station (flow).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task CHMN_FullPor_Lillooet_Flow()
        {
            if (!await ChmnAvailable()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_2);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the CHMN Capilano River station (flow).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task CHMN_FullPor_Capilano_Flow()
        {
            if (!await ChmnAvailable()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_3);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests flow unit conversions (cms ↔ cfs) for CHMN data.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        /// Verifies USGS downloads do not run a generic connectivity preflight before requesting the data.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// The test uses a deterministic in-memory HTTP handler so it can prove request routing
        /// without relying on the live USGS service or the host machine's internet connection.
        /// </remarks>
        [TestMethod]
        public async Task USGS_PeakDischarge_DoesNotRunConnectivityPreflight()
        {
            var defaultHandler = new GuardedUsgsPeakHandler();
            var decompressHandler = new GuardedUsgsPeakHandler();

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_3, TimeSeriesDownload.TimeSeriesType.PeakDischarge);

                Assert.IsNotNull(ts, "Time series is null.");
                Assert.AreEqual(1, ts.Count);
                Assert.AreEqual(12345d, ts[0].Value);
                Assert.IsFalse(string.IsNullOrWhiteSpace(raw), "Raw text should contain the fake USGS response.");
                Assert.HasCount(1, defaultHandler.Requests, "Expected exactly one USGS data request.");
                Assert.IsFalse(defaultHandler.Requests.Any(GuardedUsgsPeakHandler.IsBlockedPreflight));
                Assert.IsEmpty(decompressHandler.Requests, "Peak download should not use the decompression client.");
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
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

                Assert.HasCount(1, defaultHandler.Requests, "Expected one started USGS peak request.");
                Assert.IsEmpty(decompressHandler.Requests, "Peak download should not use the decompression client.");
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies duplicate USGS measured-stage timestamps keep the last source value.
        /// </summary>
        [TestMethod]
        public async Task USGS_MeasuredStage_DuplicateDateTimes_LastValueWins()
        {
            var defaultHandler = new DuplicateUsgsMeasuredStageHandler();
            var decompressHandler = new DuplicateUsgsMeasuredStageHandler();

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_5,
                    TimeSeriesDownload.TimeSeriesType.MeasuredStage);

                var duplicateTime = new DateTime(1936, 3, 15, 5, 0, 0, DateTimeKind.Utc);

                Assert.IsNotNull(ts, "Time series is null.");
                Assert.IsFalse(string.IsNullOrWhiteSpace(raw), "Raw text should contain the fake USGS response.");
                Assert.AreEqual(ts.Count, ts.Select(o => o.Index).Distinct().Count(), "Duplicate date indices detected.");
                Assert.AreEqual(1, ts.Count(o => o.Index == duplicateTime), "Expected one retained record for the duplicate timestamp.");
                Assert.AreEqual(14.52d, ts.Single(o => o.Index == duplicateTime).Value);
                Assert.IsEmpty(defaultHandler.Requests, "Measured stage should not use the default client.");
                Assert.HasCount(1, decompressHandler.Requests, "Expected exactly one USGS OGC data request.");
                Assert.IsTrue(decompressHandler.Requests[0].Contains("monitoring_location_id=USGS-01570500", StringComparison.OrdinalIgnoreCase));
                Assert.IsTrue(decompressHandler.Requests[0].Contains("parameter_code=00065", StringComparison.OrdinalIgnoreCase));
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies deterministic USGS instantaneous discharge downloads remain irregular.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// USGS instantaneous discharge uses parameter <c>00060</c> on the <c>/nwis/iv</c> endpoint.
        /// </remarks>
        [TestMethod]
        public async Task USGS_InstantaneousDischarge_DeterministicDownload_IsIrregular()
        {
            await VerifyUsgsInstantaneousDownloadIsIrregular(
                TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge,
                "00060");
        }

        /// <summary>
        /// Verifies deterministic USGS instantaneous stage downloads remain irregular.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// USGS instantaneous stage uses parameter <c>00065</c> on the <c>/nwis/iv</c> endpoint.
        /// </remarks>
        [TestMethod]
        public async Task USGS_InstantaneousStage_DeterministicDownload_IsIrregular()
        {
            await VerifyUsgsInstantaneousDownloadIsIrregular(
                TimeSeriesDownload.TimeSeriesType.InstantaneousStage,
                "00065");
        }

        /// <summary>
        /// Verifies instantaneous USGS timeout messages use the longer instantaneous limit.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// The handler throws a timeout immediately after URL validation so the test does not wait
        /// for the full timeout duration.
        /// </remarks>
        [TestMethod]
        public async Task USGS_InstantaneousDischarge_TimeoutMessage_UsesInstantaneousLimit()
        {
            var defaultHandler = new GuardedUsgsInstantaneousHandler("00060", simulateTimeout: true);
            var decompressHandler = new GuardedUsgsInstantaneousHandler("00060", simulateTimeout: true);

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var exception = await Assert.ThrowsAsync<TimeoutException>(async () =>
                    await TimeSeriesDownload.FromUSGS(USGS_6, TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge));

                Assert.Contains("120 seconds", exception.Message);
                Assert.IsEmpty(defaultHandler.Requests, "Instantaneous download should not use the default client.");
                Assert.HasCount(2, decompressHandler.Requests, "Expected the retry helper to make two instantaneous attempts.");
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies USGS instantaneous parsing does not resume on the caller synchronization context.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// BestFit starts downloads from the WPF UI thread. This test forces the fake HTTP response
        /// to complete asynchronously while a synchronization context is installed, then verifies
        /// the downloader did not post its large-response continuation back to that context.
        /// </remarks>
        [TestMethod]
        public async Task USGS_InstantaneousDischarge_DoesNotPostBackToSynchronizationContext()
        {
            var defaultHandler = new GuardedUsgsInstantaneousHandler("00060", completeAsynchronously: true);
            var decompressHandler = new GuardedUsgsInstantaneousHandler("00060", completeAsynchronously: true);

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var context = new CountingSynchronizationContext();
                SynchronizationContext previousContext = SynchronizationContext.Current;
                Task<(TimeSeries TimeSeries, string RawText)> downloadTask;

                try
                {
                    SynchronizationContext.SetSynchronizationContext(context);
                    downloadTask = TimeSeriesDownload.FromUSGS(
                        USGS_6,
                        TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge);
                }
                finally
                {
                    SynchronizationContext.SetSynchronizationContext(previousContext);
                }

                var (ts, raw) = await downloadTask;

                Assert.AreEqual(0, context.PostCount, "USGS instantaneous download should not post continuations to the caller synchronization context.");
                Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
                Assert.AreEqual(3, ts.Count);
                Assert.IsFalse(string.IsNullOrWhiteSpace(raw), "Raw text should contain the fake USGS response.");
                Assert.IsEmpty(defaultHandler.Requests, "Instantaneous download should not use the default client.");
                Assert.HasCount(1, decompressHandler.Requests, "Expected exactly one USGS instantaneous data request.");
                Assert.IsFalse(decompressHandler.Requests.Any(GuardedUsgsInstantaneousHandler.IsBlockedPreflight));
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies duplicate-heavy USGS instantaneous responses are normalized without a long rebuild.
        /// </summary>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// The fixture mimics provisional/approved duplicate timestamps in instantaneous RDB data.
        /// A previous cleanup path cleared the million-row series through per-item removal, which
        /// made full-period BestFit downloads appear to spin indefinitely after the response arrived.
        /// </remarks>
        [TestMethod]
        [Timeout(15000, CooperativeCancellation = true)]
        public async Task USGS_InstantaneousDischarge_DuplicateHeavyResponse_NormalizesQuickly()
        {
            const int duplicateTimestampPairs = 50000;
            var defaultHandler = new GuardedUsgsInstantaneousHandler("00060", duplicateTimestampPairs: duplicateTimestampPairs);
            var decompressHandler = new GuardedUsgsInstantaneousHandler("00060", duplicateTimestampPairs: duplicateTimestampPairs);

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var (ts, raw) = await TimeSeriesDownload.FromUSGS(
                    USGS_6,
                    TimeSeriesDownload.TimeSeriesType.InstantaneousDischarge);
                var finalDuplicateTimestamp = new DateTime(2024, 1, 2).AddMinutes(duplicateTimestampPairs - 1);

                Assert.IsNotNull(ts, "Time series is null.");
                Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
                Assert.AreEqual(3 + duplicateTimestampPairs, ts.Count);
                Assert.AreEqual(ts.Count, ts.Select(o => o.Index).Distinct().Count(), "Duplicate date indices detected.");
                Assert.AreEqual(2000d + duplicateTimestampPairs - 1, ts.Single(o => o.Index == finalDuplicateTimestamp).Value);
                Assert.IsFalse(string.IsNullOrWhiteSpace(raw), "Raw text should contain the fake USGS response.");
                Assert.IsEmpty(defaultHandler.Requests, "Instantaneous download should not use the default client.");
                Assert.HasCount(1, decompressHandler.Requests, "Expected exactly one USGS instantaneous data request.");
                Assert.IsFalse(decompressHandler.Requests.Any(GuardedUsgsInstantaneousHandler.IsBlockedPreflight));
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Verifies a deterministic USGS instantaneous download preserves nonuniform timestamps.
        /// </summary>
        /// <param name="seriesType">The instantaneous USGS series type to request.</param>
        /// <param name="parameterCode">The expected USGS parameter code.</param>
        /// <returns>A task that completes when the regression check finishes.</returns>
        /// <remarks>
        /// The fixture has a 15-minute step followed by a 30-minute step. The downloader must keep
        /// those timestamps exactly and must not convert the series into a regular interval.
        /// </remarks>
        private static async Task VerifyUsgsInstantaneousDownloadIsIrregular(
            TimeSeriesDownload.TimeSeriesType seriesType,
            string parameterCode)
        {
            var defaultHandler = new GuardedUsgsInstantaneousHandler(parameterCode);
            var decompressHandler = new GuardedUsgsInstantaneousHandler(parameterCode);

            using var defaultClient = new HttpClient(defaultHandler) { Timeout = TimeSpan.FromSeconds(5) };
            using var decompressClient = new HttpClient(decompressHandler) { Timeout = TimeSpan.FromSeconds(5) };

            TimeSeriesDownload.SetHttpClientsForTesting(defaultClient, decompressClient);
            try
            {
                var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_6, seriesType);
                var expectedDates = new[]
                {
                    new DateTime(2024, 1, 1, 0, 0, 0),
                    new DateTime(2024, 1, 1, 0, 15, 0),
                    new DateTime(2024, 1, 1, 0, 45, 0)
                };

                Assert.IsNotNull(ts, "Time series is null.");
                Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
                Assert.AreEqual(expectedDates.Length, ts.Count);
                CollectionAssert.AreEqual(expectedDates, ts.Select(o => o.Index).ToArray());
                Assert.IsFalse(string.IsNullOrWhiteSpace(raw), "Raw text should contain the fake USGS response.");
                Assert.IsEmpty(defaultHandler.Requests, "Instantaneous download should not use the default client.");
                Assert.HasCount(1, decompressHandler.Requests, "Expected exactly one USGS instantaneous data request.");
                Assert.IsTrue(decompressHandler.Requests[0].Contains("/nwis/iv/", StringComparison.OrdinalIgnoreCase));
                Assert.IsTrue(decompressHandler.Requests[0].Contains($"parameterCd={parameterCode}", StringComparison.OrdinalIgnoreCase));
                Assert.IsTrue(decompressHandler.Requests[0].Contains("startDT=1900-01-01", StringComparison.OrdinalIgnoreCase));
                Assert.IsFalse(decompressHandler.Requests.Any(GuardedUsgsInstantaneousHandler.IsBlockedPreflight));
            }
            finally
            {
                TimeSeriesDownload.ResetHttpClientsForTesting();
            }
        }

        /// <summary>
        /// Tests full-period-of-record USGS daily discharge download.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task USGS_FullPor_DailyDischarge()
        {
            if (!await UsgsAvailable()) return;

            var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests USGS daily stage download for correctness and continuity.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task USGS_FullPor_DailyStage()
        {
            if (!await UsgsAvailable()) return;

            var (ts, _) = await TimeSeriesDownload.FromUSGS(USGS_4, TimeSeriesDownload.TimeSeriesType.DailyStage);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests USGS peak discharge data retrieval for non-daily data.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
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
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromUSGS("1134500"));
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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

        /// <summary>
        /// Validates a full-period-of-record download for the BOM Cotter River station (discharge).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task BOM_FullPor_CotterRiver_Discharge()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the BOM Goodradigbee River station (discharge).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task BOM_FullPor_Goodradigbee_Discharge()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_2);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the BOM Murray River station (stage).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task BOM_FullPor_MurrayRiver_Stage()
        {
            if (!await BomAvailable()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_3, TimeSeriesDownload.TimeSeriesType.DailyStage);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests discharge unit conversions (cms ↔ cfs) for BOM data.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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
        [TestMethod, TestCategory("Integration")]
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

