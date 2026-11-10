using System.IO.Compression;
using System.IO;
using System.Linq;
using System.Net;
using System.Net.Sockets;
using System.Xml.Linq;
using System;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using System.Net.Http;
using System.Globalization;
using System.Collections.Generic;
using System.Text;

namespace Numerics.Data
{

    /// <summary>
    /// Download time series data from the Internet.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class TimeSeriesDownload
    {
        private static HttpClient _defaultClient = CreateDefaultClient();
        private static HttpClient _decompressClient = CreateDecompressClient();
        private static HttpClient _bomClient = CreateBomClient();

        /// <summary>
        /// Gets the User-Agent product token sent with every provider request, for example USACE-Numerics/2.3.
        /// </summary>
        /// <remarks>
        /// The version is the major and minor assembly version, so providers can tell library releases apart.
        /// </remarks>
        internal static string UserAgent
        {
            get
            {
                Version? version = typeof(TimeSeriesDownload).Assembly.GetName().Version;
                return "USACE-Numerics/" + (version == null ? "2" : version.ToString(2));
            }
        }

        /// <summary>
        /// Maximum time allowed for one provider request attempt before retry/failure handling runs.
        /// </summary>
        /// <remarks>
        /// A short per-attempt timeout prevents blocked government-network routes from making UI
        /// callers wait for the default minute-long HTTP timeout on each retry.
        /// </remarks>
        private static readonly TimeSpan ProviderRequestTimeout = TimeSpan.FromSeconds(10);

        /// <summary>
        /// Maximum time allowed for one instantaneous provider request attempt.
        /// </summary>
        /// <remarks>
        /// Instantaneous series can be substantially larger than daily, peak, or field-measurement
        /// responses, especially when a provider accepts a long period-of-record query.
        /// </remarks>
        private static readonly TimeSpan InstantaneousProviderRequestTimeout = TimeSpan.FromSeconds(120);

        /// <summary>
        /// Maximum time allowed for one GHCN station-file request attempt.
        /// </summary>
        /// <remarks>
        /// NCEI serves full station files in seconds over IPv4. Multi-minute
        /// header waits come from sequential TCP timeouts against NCEI's six IPv6 addresses on
        /// networks that drop IPv6 traffic, which the downloader avoids by preferring IPv4.
        /// This margin only covers genuinely slow links and large station files.
        /// </remarks>
        private static readonly TimeSpan GhcnProviderRequestTimeout = TimeSpan.FromSeconds(60);

        /// <summary>
        /// Maximum time allowed for one explicit connectivity probe.
        /// </summary>
        /// <remarks>
        /// Connectivity probes are advisory only and must remain faster than real provider downloads.
        /// </remarks>
        private static readonly TimeSpan ConnectivityProbeTimeout = TimeSpan.FromSeconds(2);

        /// <summary>
        /// Default number of attempts for provider requests made through the retry helper.
        /// </summary>
        /// <remarks>
        /// Two attempts are enough to smooth over one transient outage without turning blocked
        /// network paths into multi-minute waits.
        /// </remarks>
        private const int DefaultMaxDownloadAttempts = 2;

        /// <summary>
        /// Provider root endpoints used for the public Internet connectivity check.
        /// </summary>
        /// <remarks>
        /// These URLs are not a precondition for data requests; they are used only by callers that
        /// explicitly ask whether at least one supported provider endpoint is reachable.
        /// </remarks>
        private static readonly string[] InternetProbeUrls =
        {
            "https://api.waterdata.usgs.gov/",
            "https://www.ncei.noaa.gov/",
            "https://wateroffice.ec.gc.ca/",
            "https://www.bom.gov.au/waterdata/"
        };

#if !NETFRAMEWORK
        /// <summary>
        /// Maximum time allowed for one TCP connect attempt to a single resolved address.
        /// </summary>
        /// <remarks>
        /// Bounding each address separately keeps one unreachable address family from consuming
        /// the whole request timeout while still letting later addresses in the list be tried.
        /// </remarks>
        private static readonly TimeSpan ConnectAttemptTimeout = TimeSpan.FromSeconds(5);
#endif

        /// <summary>
        /// Creates the HTTP handler used by the downloader's shared clients.
        /// </summary>
        /// <param name="automaticDecompression">True to enable automatic GZip and Deflate response decompression.</param>
        /// <returns>A configured HTTP handler.</returns>
        /// <remarks>
        /// Some providers publish IPv6 (AAAA) DNS records while many restricted networks silently
        /// drop IPv6 traffic. A default handler tries every unreachable IPv6 address before
        /// falling back to IPv4 (~21 s of TCP timeout per address; measured at ~126 s to first
        /// byte for NOAA NCEI's six AAAA records). On modern targets the handler connects with
        /// <c>ConnectPreferringIPv4Async</c> so IPv4 is tried first and IPv6 remains a
        /// fallback. On .NET Framework the equivalent ordering is applied per request through
        /// <see cref="PreferIPv4ForNetFramework"/>.
        /// </remarks>
        private static HttpMessageHandler CreateHandler(bool automaticDecompression)
        {
#if NETFRAMEWORK
            var handler = new HttpClientHandler();
#else
            var handler = new SocketsHttpHandler { ConnectCallback = ConnectPreferringIPv4Async };
#endif
            if (automaticDecompression)
            {
                handler.AutomaticDecompression = DecompressionMethods.GZip | DecompressionMethods.Deflate;
            }
            return handler;
        }

#if !NETFRAMEWORK
        /// <summary>
        /// Opens a TCP connection for an HTTP request, preferring IPv4 addresses over IPv6.
        /// </summary>
        /// <param name="context">The connection context supplied by the HTTP handler.</param>
        /// <param name="cancellationToken">Token used to cancel the connection attempt.</param>
        /// <returns>A connected network stream.</returns>
        /// <exception cref="HttpRequestException">Thrown when no resolved address accepts a connection.</exception>
        /// <exception cref="OperationCanceledException">Thrown when <paramref name="cancellationToken"/> is canceled.</exception>
        /// <remarks>
        /// The default connect logic tries DNS results in returned order, which lists IPv6 first.
        /// On networks that silently drop IPv6 traffic every dead address costs a full TCP
        /// timeout before IPv4 is reached. Trying IPv4 first keeps downloads fast there while the
        /// IPv6 fallback keeps IPv6-only networks working.
        /// </remarks>
        private static async ValueTask<Stream> ConnectPreferringIPv4Async(SocketsHttpConnectionContext context, CancellationToken cancellationToken)
        {
            IPAddress[] addresses = await Dns.GetHostAddressesAsync(context.DnsEndPoint.Host, cancellationToken).ConfigureAwait(false);
            IPAddress[] ordered = OrderAddressesForConnect(addresses);
            Exception? lastFailure = null;
            foreach (IPAddress address in ordered)
            {
                var socket = new Socket(address.AddressFamily, SocketType.Stream, ProtocolType.Tcp) { NoDelay = true };
                try
                {
                    using var attemptCts = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
                    attemptCts.CancelAfter(ConnectAttemptTimeout);
                    await socket.ConnectAsync(new IPEndPoint(address, context.DnsEndPoint.Port), attemptCts.Token).ConfigureAwait(false);
                    return new NetworkStream(socket, ownsSocket: true);
                }
                catch (Exception ex)
                {
                    socket.Dispose();
                    cancellationToken.ThrowIfCancellationRequested();
                    lastFailure = ex;
                }
            }

            throw new HttpRequestException(
                $"Could not connect to '{context.DnsEndPoint.Host}:{context.DnsEndPoint.Port}' on any resolved address.",
                lastFailure);
        }
#endif

        /// <summary>
        /// Orders resolved addresses so IPv4 is attempted before IPv6.
        /// </summary>
        /// <param name="addresses">The addresses returned by DNS resolution.</param>
        /// <returns>The addresses with IPv4 first; relative order within each family is preserved.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="addresses"/> is null.</exception>
        internal static IPAddress[] OrderAddressesForConnect(IPAddress[] addresses)
        {
            if (addresses == null) throw new ArgumentNullException(nameof(addresses));
            return addresses.OrderBy(address => address.AddressFamily == AddressFamily.InterNetworkV6 ? 1 : 0).ToArray();
        }

        /// <summary>
        /// Local-endpoint binder that rejects IPv6 remote endpoints so connects fall through to IPv4.
        /// </summary>
        /// <param name="remoteEndPoint">The remote endpoint the connection is being bound for.</param>
        /// <returns>Null for IPv4 remote endpoints, letting the OS choose the local endpoint.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="remoteEndPoint"/> is null.</exception>
        /// <exception cref="SocketException">Thrown for IPv6 remote endpoints to make the connection attempt skip to the next resolved address.</exception>
        /// <remarks>
        /// On .NET Framework a <see cref="ServicePoint"/> tries resolved addresses in DNS order
        /// (IPv6 first) with no per-address bound, so blocked IPv6 routes cost a full TCP timeout
        /// per address. Throwing from the bind delegate fails an IPv6 attempt immediately and the
        /// service point moves on to the next (IPv4) address.
        /// </remarks>
        internal static IPEndPoint? RejectIPv6BindEndPoint(IPEndPoint remoteEndPoint)
        {
            if (remoteEndPoint == null) throw new ArgumentNullException(nameof(remoteEndPoint));
            if (remoteEndPoint.AddressFamily == AddressFamily.InterNetworkV6)
            {
                throw new SocketException((int)SocketError.AddressFamilyNotSupported);
            }
            return null;
        }

        /// <summary>
        /// Configures the .NET Framework service point for a request URI to prefer IPv4 connections.
        /// </summary>
        /// <param name="requestUri">The absolute request URI about to be downloaded.</param>
        /// <remarks>
        /// No-op on modern targets, where <see cref="CreateHandler"/> installs an IPv4-preferring
        /// connect callback instead. Skipped when the OS has no IPv4 stack so IPv6-only machines
        /// keep working.
        /// </remarks>
        private static void PreferIPv4ForNetFramework(Uri requestUri)
        {
#if NETFRAMEWORK
            if (!Socket.OSSupportsIPv4) return;
            ServicePoint servicePoint = ServicePointManager.FindServicePoint(requestUri);
            servicePoint.BindIPEndPointDelegate = (servicePointArg, remoteEndPoint, retryCount) => RejectIPv6BindEndPoint(remoteEndPoint);
#endif
        }

        /// <summary>
        /// Creates the default HTTP client used for uncompressed provider requests.
        /// </summary>
        /// <returns>A configured HTTP client whose timeout is controlled by the retry helper.</returns>
        /// <remarks>
        /// The client is held as a static instance in production to avoid socket churn. Request
        /// timeouts are applied by <see cref="GetWithRetryAsync(HttpClient, string, int, TimeSpan?, CancellationToken, string, Func{string, bool})"/>
        /// so each provider call can choose the right timeout for its payload size.
        /// </remarks>
        private static HttpClient CreateDefaultClient()
        {
            var client = new HttpClient(CreateHandler(automaticDecompression: false)) { Timeout = Timeout.InfiniteTimeSpan };
            client.DefaultRequestHeaders.UserAgent.ParseAdd(UserAgent);
            return client;
        }

        /// <summary>
        /// Creates the HTTP client used for provider requests that may return compressed payloads.
        /// </summary>
        /// <returns>A configured HTTP client with automatic GZip and Deflate decompression.</returns>
        /// <remarks>
        /// The client is held as a static instance in production to avoid socket churn. Request
        /// timeouts are applied by <see cref="GetWithRetryAsync(HttpClient, string, int, TimeSpan?, CancellationToken, string, Func{string, bool})"/>.
        /// </remarks>
        private static HttpClient CreateDecompressClient()
        {
            var client = new HttpClient(CreateHandler(automaticDecompression: true)) { Timeout = Timeout.InfiniteTimeSpan };
            client.DefaultRequestHeaders.UserAgent.ParseAdd(UserAgent);
            return client;
        }

        /// <summary>
        /// Creates the HTTP client used for BOM KiWIS requests.
        /// </summary>
        /// <returns>A configured HTTP client with browser-like default headers.</returns>
        /// <remarks>
        /// BOM sits behind a security appliance that rejects bare programmatic requests, so the
        /// client sends browser-like headers. The client is held as a static instance so the
        /// station lookup and value download in one call, and repeated calls, reuse pooled
        /// connections instead of paying a fresh TCP and TLS handshake each time. The
        /// Accept-Encoding header is supplied by the handler's automatic decompression.
        /// </remarks>
        private static HttpClient CreateBomClient()
        {
            var client = new HttpClient(CreateHandler(automaticDecompression: true)) { Timeout = Timeout.InfiniteTimeSpan };
            client.DefaultRequestHeaders.Add("User-Agent", UserAgent);
            client.DefaultRequestHeaders.Add("Accept", "application/json,text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
            client.DefaultRequestHeaders.Add("Accept-Language", "en-US,en;q=0.9");
            client.DefaultRequestHeaders.Add("DNT", "1");
            client.DefaultRequestHeaders.Add("Connection", "keep-alive");
            return client;
        }

        /// <summary>
        /// Replaces the static HTTP clients used by the downloader for deterministic tests.
        /// </summary>
        /// <param name="defaultClient">The client to use for ordinary provider requests.</param>
        /// <param name="decompressClient">The client to use for compressed provider requests.</param>
        /// <param name="bomClient">The client to use for BOM requests, or null to keep the current one.</param>
        /// <exception cref="ArgumentNullException">Thrown when the default or decompressing client is null.</exception>
        /// <remarks>
        /// This method is internal so tests can verify request routing without adding public API.
        /// </remarks>
        internal static void SetHttpClientsForTesting(HttpClient defaultClient, HttpClient decompressClient, HttpClient? bomClient = null)
        {
            _defaultClient = defaultClient ?? throw new ArgumentNullException(nameof(defaultClient));
            _decompressClient = decompressClient ?? throw new ArgumentNullException(nameof(decompressClient));
            if (bomClient != null) _bomClient = bomClient;
        }

        /// <summary>
        /// Restores the production HTTP clients after a deterministic test.
        /// </summary>
        /// <remarks>
        /// Tests call this method in a finally block so later tests do not inherit fake handlers.
        /// </remarks>
        internal static void ResetHttpClientsForTesting()
        {
            _defaultClient = CreateDefaultClient();
            _decompressClient = CreateDecompressClient();
            _bomClient = CreateBomClient();
        }

        /// <summary>
        /// GET with bounded retry and exponential backoff on transient failures. Retries 5xx, 408, 429, network
        /// exceptions, and KiWIS DatasourceError bodies. Does NOT retry 4xx:
        /// those mean the request itself is wrong.
        /// </summary>
        /// <param name="client">The HTTP client used to issue the request.</param>
        /// <param name="url">The request URL.</param>
        /// <param name="maxAttempts">The maximum number of request attempts.</param>
        /// <param name="requestTimeout">The maximum duration allowed for each request attempt.</param>
        /// <param name="cancellationToken">Token used to cancel the request and retry delay.</param>
        /// <param name="apiKey">An API key sent in the <c>X-Api-Key</c> header, or null for none.</param>
        /// <param name="isComplete">
        /// Optional check of a successful body; a body it rejects (for example one cut off mid-record) is retried like a
        /// transient failure, and on the last attempt an <see cref="IOException"/> is thrown.
        /// </param>
        /// <returns>The final HTTP status code and response body.</returns>
        /// <exception cref="HttpRequestException">Thrown when all attempts fail with a network error.</exception>
        /// <exception cref="IOException">Thrown when the last attempt's body ends early or fails <paramref name="isComplete"/>.</exception>
        /// <exception cref="OperationCanceledException">Thrown when <paramref name="cancellationToken"/> is canceled.</exception>
        /// <exception cref="TimeoutException">Thrown when attempts exceed the per-attempt timeout twice, or on the last attempt.</exception>
        /// <remarks>
        /// The helper bounds each attempt separately so blocked network paths fail quickly while
        /// still allowing retries for transient service errors. A transient response with a
        /// <c>Retry-After</c> delay of up to 30 seconds waits that long before the next attempt.
        /// A connection that drops while the body is read (an <see cref="IOException"/>) is transient.
        /// A timed-out attempt is retried at most once, so a stalled service costs at most two
        /// timeouts instead of one per attempt.
        /// </remarks>
        private static async Task<(HttpStatusCode status, string body)> GetWithRetryAsync(
            HttpClient client,
            string url,
            int maxAttempts = DefaultMaxDownloadAttempts,
            TimeSpan? requestTimeout = null,
            CancellationToken cancellationToken = default,
            string? apiKey = null,
            Func<string, bool>? isComplete = null)
        {
            TimeSpan timeout = requestTimeout ?? ProviderRequestTimeout;
            PreferIPv4ForNetFramework(new Uri(url));
            Exception? lastNet = null;
            int timeouts = 0;
            for (int attempt = 1; attempt <= maxAttempts; attempt++)
            {
                TimeSpan? retryAfter = null;
                try
                {
                    using var attemptCts = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
                    attemptCts.CancelAfter(timeout);
                    using var message = new HttpRequestMessage(HttpMethod.Get, url);
                    if (apiKey != null) message.Headers.TryAddWithoutValidation("X-Api-Key", apiKey);
                    using var resp = await client.SendAsync(
                        message,
                        HttpCompletionOption.ResponseHeadersRead,
                        attemptCts.Token).ConfigureAwait(false);
                    string body = await ReadResponseBodyAsync(resp, attemptCts.Token).ConfigureAwait(false);
                    bool transient =
                        (int)resp.StatusCode >= 500 ||
                        resp.StatusCode == HttpStatusCode.RequestTimeout ||
                        (int)resp.StatusCode == 429 ||
                        body.Contains("DatasourceError") ||
                        body.Contains("Error connecting to WDP");
                    retryAfter = resp.Headers.RetryAfter?.Delta;
                    if (!transient && resp.IsSuccessStatusCode && isComplete != null && !isComplete(body))
                    {
                        var truncated = new IOException($"The response from '{url}' ended before its last record.");
                        if (attempt == maxAttempts) throw truncated;
                        lastNet = truncated;
                    }
                    else if (!transient || attempt == maxAttempts || IsLongRetryAfter(retryAfter))
                    {
                        // A rate limit that resets later than a short wait is final for this call.
                        return (resp.StatusCode, body);
                    }
                }
                catch (HttpRequestException ex) { lastNet = ex; if (attempt == maxAttempts) throw; }
                catch (IOException ex) when (attempt < maxAttempts) { lastNet = ex; }
                catch (TaskCanceledException) when (cancellationToken.IsCancellationRequested) { throw; }
                catch (TaskCanceledException ex)
                {
                    lastNet = CreateDownloadTimeoutException(url, timeout, ex);
                    if (attempt == maxAttempts || ++timeouts >= 2) throw lastNet;
                }

                await Task.Delay(RetryDelay(retryAfter, attempt), cancellationToken).ConfigureAwait(false);
            }
            throw lastNet ?? new Exception("GetWithRetryAsync exhausted retries");
        }

        /// <summary>
        /// The longest <c>Retry-After</c> delay a download waits before retrying.
        /// </summary>
        private static readonly TimeSpan MaxRetryAfter = TimeSpan.FromSeconds(30);

        /// <summary>
        /// Determines whether a <c>Retry-After</c> delay is too long to wait for, which makes the response final.
        /// </summary>
        /// <param name="retryAfter">The delay from the response, or null.</param>
        /// <returns>True when the delay is longer than <see cref="MaxRetryAfter"/>.</returns>
        private static bool IsLongRetryAfter(TimeSpan? retryAfter) => retryAfter.HasValue && retryAfter.Value > MaxRetryAfter;

        /// <summary>
        /// Returns the wait before the next attempt: a short <c>Retry-After</c> delay when the service gave one,
        /// otherwise exponential backoff of 0.5, 1, 2, ... seconds.
        /// </summary>
        /// <param name="retryAfter">The delay from the response, or null.</param>
        /// <param name="attempt">The one-based attempt that just failed.</param>
        /// <returns>The delay.</returns>
        private static TimeSpan RetryDelay(TimeSpan? retryAfter, int attempt)
        {
            return retryAfter.HasValue && retryAfter.Value > TimeSpan.Zero && retryAfter.Value <= MaxRetryAfter
                ? retryAfter.Value
                : TimeSpan.FromMilliseconds(500 * Math.Pow(2, attempt - 1));
        }

        /// <summary>
        /// Determines whether a USGS CSV page is complete: empty, or ending with a line break.
        /// </summary>
        /// <param name="body">The page text.</param>
        /// <returns>False for a page cut off inside its last record.</returns>
        /// <remarks>Every USGS Water Data API CSV page ends with a line break.</remarks>
        private static bool IsCompleteUsgsCsv(string body) => body.Length == 0 || body[body.Length - 1] == '\n';

        /// <summary>
        /// Reads an HTTP response body with cancellation support.
        /// </summary>
        /// <param name="response">The HTTP response whose content should be read.</param>
        /// <param name="cancellationToken">Token used to cancel the body read.</param>
        /// <returns>The decoded response text.</returns>
        /// <exception cref="OperationCanceledException">Thrown when <paramref name="cancellationToken"/> is canceled.</exception>
        /// <remarks>
        /// <see cref="HttpContent.ReadAsStringAsync()"/> cannot be canceled on older target
        /// frameworks. Full-period instantaneous USGS responses are large enough that body reads
        /// must honor the same bounded timeout as the request headers.
        /// </remarks>
        private static async Task<string> ReadResponseBodyAsync(HttpResponseMessage response, CancellationToken cancellationToken)
        {
            using var stream = await response.Content.ReadAsStreamAsync().ConfigureAwait(false);
            using var memory = new MemoryStream();
            byte[] buffer = new byte[81920];
            int bytesRead;
            while ((bytesRead = await stream.ReadAsync(buffer, 0, buffer.Length, cancellationToken).ConfigureAwait(false)) > 0)
            {
                memory.Write(buffer, 0, bytesRead);
            }

            return GetResponseEncoding(response).GetString(memory.ToArray());
        }

        /// <summary>
        /// Gets the text encoding declared by an HTTP response, defaulting to UTF-8.
        /// </summary>
        /// <param name="response">The HTTP response to inspect.</param>
        /// <returns>The declared response encoding, or UTF-8 when none is declared or recognized.</returns>
        /// <remarks>
        /// USGS RDB and most provider text responses are UTF-8-compatible. The fallback keeps parsing
        /// stable when a provider omits or misstates the charset.
        /// </remarks>
        private static Encoding GetResponseEncoding(HttpResponseMessage response)
        {
            string charset = response.Content.Headers.ContentType?.CharSet ?? "";
            if (string.IsNullOrWhiteSpace(charset)) return Encoding.UTF8;

            try
            {
                return Encoding.GetEncoding(charset.Trim('"'));
            }
            catch (ArgumentException)
            {
                return Encoding.UTF8;
            }
        }

        /// <summary>
        /// Creates a timeout exception for a provider request that did not complete quickly enough.
        /// </summary>
        /// <param name="url">The URL that timed out.</param>
        /// <param name="requestTimeout">The timeout used for the request attempt.</param>
        /// <param name="innerException">The cancellation exception raised by the HTTP stack.</param>
        /// <returns>An exception with a user-facing timeout message.</returns>
        /// <remarks>
        /// The message names the endpoint so application callers can distinguish a blocked provider
        /// request from a general Internet connectivity failure.
        /// </remarks>
        private static TimeoutException CreateDownloadTimeoutException(string url, TimeSpan requestTimeout, Exception innerException)
        {
            return new TimeoutException(
                $"The request to '{url}' did not complete within {requestTimeout.TotalSeconds:0} seconds.",
                innerException);
        }

        /// <summary>
        /// Checks if there is an Internet connection by probing provider endpoints used by the downloader.
        /// </summary>
        /// <returns>True when at least one supported provider endpoint completes an HTTP request; otherwise false.</returns>
        /// <remarks>
        /// This method is not used as a preflight gate for provider downloads because restricted
        /// networks may allow one data provider while blocking unrelated endpoints.
        /// </remarks>
        public static async Task<bool> IsConnectedToInternet()
        {
            using var cts = new CancellationTokenSource();
            var probes = InternetProbeUrls.Select(url => CanReachEndpoint(url, cts.Token)).ToList();

            while (probes.Count > 0)
            {
                Task<bool> completed = await Task.WhenAny(probes);
                probes.Remove(completed);
                if (await completed)
                {
                    cts.Cancel();
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Determines whether an endpoint can complete an HTTP request.
        /// </summary>
        /// <param name="url">The endpoint URL to probe.</param>
        /// <param name="cancellationToken">Token used to cancel the probe.</param>
        /// <returns>True when the request completes with any HTTP response; otherwise false.</returns>
        /// <remarks>
        /// Any HTTP response means the network path exists. HTTP status codes are intentionally not
        /// interpreted because security appliances and provider roots may return redirects or errors
        /// while still proving transport reachability.
        /// </remarks>
        private static async Task<bool> CanReachEndpoint(string url, CancellationToken cancellationToken = default)
        {
            try
            {
                PreferIPv4ForNetFramework(new Uri(url));
                using var cts = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
                cts.CancelAfter(ConnectivityProbeTimeout);
                using (await _defaultClient.GetAsync(
                    url,
                    HttpCompletionOption.ResponseHeadersRead,
                    cts.Token))
                {
                }

                return true;
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        /// Writes the qualifier of each ordinate, in series order, into a caller's list.
        /// </summary>
        /// <param name="timeSeries">The parsed series.</param>
        /// <param name="qualifierByDate">The qualifier of each parsed record by date; the last record for a date wins, as for values.</param>
        /// <param name="qualifiers">The caller's list; cleared, then given one entry per ordinate (empty for filled dates).</param>
        private static void WriteQualifiers(TimeSeries timeSeries, Dictionary<DateTime, string> qualifierByDate, IList<string> qualifiers)
        {
            qualifiers.Clear();
            foreach (var ordinate in timeSeries)
                qualifiers.Add(qualifierByDate.TryGetValue(ordinate.Index, out string? qualifier) ? qualifier : "");
        }

        /// <summary>
        /// The message of the exception thrown when a provider returns no data.
        /// </summary>
        private const string NoDataMessage = "No data found matching all criteria.";

        /// <summary>
        /// Determines whether a successful text response holds no data.
        /// </summary>
        /// <param name="body">The response text.</param>
        /// <param name="hasHeader">True when the format starts with a header line (CSV).</param>
        /// <returns>True for an empty body, an HTML page, or, with a header, a body with no line after the header.</returns>
        /// <remarks>
        /// Some providers answer an unknown station with an empty file or a header-only file, and proxies or security
        /// appliances can answer with an HTML page and status 200. These are reported as missing data instead of
        /// returning an empty series.
        /// </remarks>
        private static bool IsNoDataResponse(string body, bool hasHeader)
        {
            string text = body.TrimStart('\uFEFF', ' ', '\t', '\r', '\n');
            if (text.Length == 0 || text[0] == '<') return true;
            if (!hasHeader) return false;
            int lineEnd = text.IndexOf('\n');
            return lineEnd < 0 || text.Substring(lineEnd + 1).Trim().Length == 0;
        }

        /// <summary>
        /// Adds a downloaded ordinate or replaces the previously parsed ordinate at the same date/time.
        /// </summary>
        /// <param name="timeSeries">The time series being populated.</param>
        /// <param name="indexesByDateTime">Lookup from parsed date/time values to their current series index.</param>
        /// <param name="index">The ordinate date/time.</param>
        /// <param name="value">The ordinate value.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="timeSeries"/> or <paramref name="indexesByDateTime"/> is null.</exception>
        /// <remarks>
        /// Provider data can contain distinct source records with the same timestamp. Handling that
        /// as records are parsed keeps the "last value wins" behavior without a
        /// second full-series pass after download completion.
        /// </remarks>
        private static void AddOrReplaceByDateTime(
            TimeSeries timeSeries,
            Dictionary<DateTime, int> indexesByDateTime,
            DateTime index,
            double value)
        {
            if (timeSeries == null) throw new ArgumentNullException(nameof(timeSeries));
            if (indexesByDateTime == null) throw new ArgumentNullException(nameof(indexesByDateTime));

            var ordinate = new SeriesOrdinate<DateTime, double>(index, value);
            if (indexesByDateTime.TryGetValue(index, out int existingIndex))
            {
                timeSeries[existingIndex] = ordinate;
                return;
            }

            indexesByDateTime.Add(index, timeSeries.Count);
            timeSeries.Add(ordinate);
        }

        /// <summary>
        /// Enumeration of time series options.
        /// </summary>
        public enum TimeSeriesType
        {
            /// <summary>
            /// Daily mean discharge.
            /// </summary>
            DailyDischarge,
            /// <summary>
            /// Daily mean stage.
            /// </summary>
            DailyStage,
            /// <summary>
            /// Daily total precipitation.
            /// </summary>
            DailyPrecipitation,
            /// <summary>
            /// Daily total snow.
            /// </summary>
            DailySnow,
            /// <summary>
            /// Instantaneous discharge, typically recorded at a 15-minute interval.
            /// </summary>
            InstantaneousDischarge,
            /// <summary>
            /// Instantaneous stage, typically recorded at a 15-minute interval.
            /// </summary>
            InstantaneousStage,
            /// <summary>
            /// Annual max peak discharge.
            /// </summary>
            PeakDischarge,
            /// <summary>
            /// Annual max peak stage.
            /// </summary>
            PeakStage,
            /// <summary>
            /// Discrete field measurement of discharge, used for rating curves.
            /// </summary>
            MeasuredDischarge,
            /// <summary>
            /// Discrete field measurement of stage (gage height), used for rating curves.
            /// </summary>
            MeasuredStage
        }

        /// <summary>
        /// Enumeration of depth unit options.
        /// </summary>
        public enum DepthUnit
        {
            /// <summary>
            /// Millimeters (mm)
            /// </summary>
            Millimeters,
            /// <summary>
            /// Centimeters (cm).
            /// </summary>
            Centimeters,
            /// <summary>
            /// Inches (in).
            /// </summary>
            Inches
        }

        /// <summary>
        /// Enumeration of discharge unit options.
        /// </summary>
        public enum DischargeUnit
        {
            /// <summary>
            /// Cubic feet per second (cfs)
            /// </summary>
            CubicFeetPerSecond,
            /// <summary>
            /// Cubic meters per second (cms)
            /// </summary>
            CubicMetersPerSecond,
        }

        /// <summary>
        /// Enumeration of gage height unit options.
        /// </summary>
        public enum HeightUnit
        {
            /// <summary>
            /// Feet (ft)
            /// </summary>
            Feet,
            /// <summary>
            /// Meters (m)
            /// </summary>
            Meters,
        }

        #region GHCN

        /// <summary>
        /// Download data from the Global Historical Climatology Network (GHCN). 
        /// </summary>
        /// <param name="siteNumber">The station identification code.</param>
        /// <param name="timeSeriesType">The time series type. Default = Daily precipitation.</param>
        /// <param name="unit">The depth unit. Default = inches.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>A downloaded time series.</returns>
        /// <exception cref="Exception">Thrown when the service returns an empty file or an HTML page.</exception>
        public static Task<TimeSeries> FromGHCN(string siteNumber, TimeSeriesType timeSeriesType = TimeSeriesType.DailyPrecipitation, DepthUnit unit = DepthUnit.Inches, CancellationToken cancellationToken = default)
        {
            return FromGhcnCore(siteNumber, timeSeriesType, unit, null, cancellationToken);
        }

        /// <summary>
        /// Download data from the Global Historical Climatology Network (GHCN), with the flags of each day.
        /// </summary>
        /// <param name="siteNumber">The station identification code.</param>
        /// <param name="timeSeriesType">The time series type: daily precipitation or daily snow.</param>
        /// <param name="unit">The depth unit.</param>
        /// <param name="qualifiers">
        /// Receives one entry per ordinate of the returned series, in the same order: the three flag characters of the
        /// day (measurement, quality and source flags, as in the station file, for example "T 1" for a trace or " I0"
        /// for a failed internal-consistency check), or an empty string for a day not in the file. The list is cleared first.
        /// </param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The same series as the overload without flags.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="qualifiers"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for an invalid station code or type.</exception>
        /// <exception cref="Exception">Thrown when the service returns an empty file or an HTML page.</exception>
        public static Task<TimeSeries> FromGHCN(string siteNumber, TimeSeriesType timeSeriesType, DepthUnit unit, IList<string> qualifiers, CancellationToken cancellationToken = default)
        {
            if (qualifiers == null) throw new ArgumentNullException(nameof(qualifiers));
            return FromGhcnCore(siteNumber, timeSeriesType, unit, qualifiers, cancellationToken);
        }

        /// <summary>
        /// Validates and runs a GHCN download.
        /// </summary>
        /// <param name="siteNumber">The station identification code.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="unit">The depth unit.</param>
        /// <param name="qualifiers">Receives the flags of each ordinate, or null.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The downloaded series.</returns>
        private static async Task<TimeSeries> FromGhcnCore(string siteNumber, TimeSeriesType timeSeriesType, DepthUnit unit, IList<string>? qualifiers, CancellationToken cancellationToken)
        {
            // Check site number
            if (siteNumber.Length != 11 || !Regex.IsMatch(siteNumber, @"^[A-Za-z0-9]{11}$"))
            {
                throw new ArgumentException("The GHCN site number must be exactly 11 alphanumeric characters.", nameof(siteNumber));
            }

            ValidateGhcnType(timeSeriesType);

            string tempFilePath = Path.Combine(Path.GetTempPath(), $"{Path.GetRandomFileName()}.dly");

            try
            {
 
                // Download the GHCN file
                string ghcnBaseUrl = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/";
                string stationFileUrl = $"{ghcnBaseUrl}{Uri.EscapeDataString(siteNumber)}.dly";
                await DownloadGhcnFileAsync(stationFileUrl, tempFilePath, cancellationToken).ConfigureAwait(false);

                if (!File.Exists(tempFilePath))
                {
                    throw new Exception("Downloaded file does not exist.");
                }

                string dlyText = File.ReadAllText(tempFilePath);
                if (IsNoDataResponse(dlyText, hasHeader: false)) throw new Exception(NoDataMessage);
                return ParseGhcnText(dlyText, timeSeriesType, unit, DateTime.Today, qualifiers);
            }
            finally
            {
                // Ensure temporary file is deleted even if an exception occurs
                if (File.Exists(tempFilePath))
                {
                    File.Delete(tempFilePath);
                }
            }
        }

        /// <summary>
        /// Downloads a GHCN station file to a local file, with bounded retries.
        /// </summary>
        /// <param name="stationFileUrl">The station file URL.</param>
        /// <param name="filePath">The local file to write.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>A task that completes when the whole file is written.</returns>
        /// <exception cref="HttpRequestException">Thrown for an error status, or when the last attempt fails with a network error.</exception>
        /// <exception cref="IOException">Thrown when the connection drops while the last attempt's body is read.</exception>
        /// <exception cref="TimeoutException">Thrown when the last attempt's headers or body take longer than the GHCN limit.</exception>
        /// <exception cref="Exception">Thrown when the service rate limit is spent.</exception>
        /// <remarks>
        /// <para>
        /// Station files are streamed to disk rather than buffered. The headers and then the body each get the GHCN
        /// request timeout, so a response that stalls part way through fails instead of waiting indefinitely. A
        /// connection that drops while the body is read is retried like a network error.
        /// </para>
        /// <para>
        /// 5xx, 408 and 429 responses are retried with the same rule as the other providers: a <c>Retry-After</c> delay
        /// of up to 30 seconds is honoured, and a longer one makes the response final. A final 429 reports that the
        /// rate limit is spent.
        /// </para>
        /// </remarks>
        private static async Task DownloadGhcnFileAsync(string stationFileUrl, string filePath, CancellationToken cancellationToken)
        {
            PreferIPv4ForNetFramework(new Uri(stationFileUrl));
            Exception? lastNet = null;
            for (int attempt = 1; ; attempt++)
            {
                bool last = attempt == DefaultMaxDownloadAttempts;
                bool final = false;
                TimeSpan? retryAfter = null;
                try
                {
                    using var attemptCts = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
                    attemptCts.CancelAfter(GhcnProviderRequestTimeout);
                    using var response = await _defaultClient.GetAsync(
                        stationFileUrl,
                        HttpCompletionOption.ResponseHeadersRead,
                        attemptCts.Token).ConfigureAwait(false);
                    int code = (int)response.StatusCode;
                    bool transient = code >= 500 || response.StatusCode == HttpStatusCode.RequestTimeout || code == 429;
                    retryAfter = response.Headers.RetryAfter?.Delta;
                    if (!transient || last || IsLongRetryAfter(retryAfter))
                    {
                        final = true;
                        if (code == 429)
                        {
                            throw new Exception("Failed to retrieve GHCN data: the NOAA NCEI service rate limit was exceeded. " +
                                $"Try again later. Status=429. Url={stationFileUrl}", lastNet);
                        }
                        response.EnsureSuccessStatusCode();

                        // The body gets its own timeout, started when the headers arrive.
                        final = false;
                        attemptCts.CancelAfter(GhcnProviderRequestTimeout);
                        using (Stream stream = await response.Content.ReadAsStreamAsync().ConfigureAwait(false))
                        using (FileStream fs = new FileStream(filePath, FileMode.Create, FileAccess.Write, FileShare.None, bufferSize: 81920, useAsync: true))
                        {
                            await stream.CopyToAsync(fs, 81920, attemptCts.Token).ConfigureAwait(false);
                        }
                        return;
                    }
                }
                catch (HttpRequestException ex) when (!last && !final) { lastNet = ex; }
                catch (IOException ex) when (!last) { lastNet = ex; }
                catch (OperationCanceledException) when (cancellationToken.IsCancellationRequested) { throw; }
                catch (OperationCanceledException ex)
                {
                    if (last) throw CreateDownloadTimeoutException(stationFileUrl, GhcnProviderRequestTimeout, ex);
                }
                await Task.Delay(RetryDelay(retryAfter, attempt), cancellationToken).ConfigureAwait(false);
            }
        }

        /// <summary>
        /// Parses the text of a GHCN-Daily station file (<c>.dly</c>) into a daily time series.
        /// </summary>
        /// <param name="dlyText">The station file text, one fixed-width record per line.</param>
        /// <param name="timeSeriesType">The element to read: daily precipitation (PRCP) or daily snow (SNOW). Default = Daily precipitation.</param>
        /// <param name="unit">The depth unit. Default = inches.</param>
        /// <returns>A regular one-day series; days without a valid value are NaN.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="dlyText"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the type is not daily precipitation or daily snow.</exception>
        /// <exception cref="FormatException">Thrown for a record too short to hold its station, date, element and days, for example a file cut off mid-line, or a record with an invalid year or month.</exception>
        /// <remarks>
        /// <para>
        /// This is the parser <see cref="FromGHCN(string, TimeSeriesType, DepthUnit, CancellationToken)"/> applies to the
        /// downloaded file, so a saved station file gives the same series offline.
        /// </para>
        /// <para>
        /// Each record holds one station, year, month and element, followed by 31 daily values of five characters plus
        /// three flag characters (NOAA NCEI GHCN-Daily readme, section III). Precipitation is stored in tenths of
        /// millimeters and snowfall in millimeters; both are converted to <paramref name="unit"/>. The missing value -9999 becomes NaN, every
        /// day between two records is filled with NaN, and a later record for the same date replaces an earlier one.
        /// Days after today (the rest of the current month, which the file pads with -9999) are left out. The flags
        /// are not read.
        /// </para>
        /// </remarks>
        public static TimeSeries ParseGHCN(string dlyText, TimeSeriesType timeSeriesType = TimeSeriesType.DailyPrecipitation, DepthUnit unit = DepthUnit.Inches)
        {
            return ParseGhcnText(dlyText, timeSeriesType, unit, DateTime.Today, null);
        }

        /// <summary>
        /// Parses the text of a GHCN-Daily station file into a daily time series, with the flags of each day.
        /// </summary>
        /// <param name="dlyText">The station file text.</param>
        /// <param name="timeSeriesType">Daily precipitation (PRCP) or daily snow (SNOW).</param>
        /// <param name="unit">The depth unit.</param>
        /// <param name="qualifiers">
        /// Receives one entry per ordinate of the returned series, in the same order: the day's three flag characters
        /// (MFLAG, QFLAG, SFLAG), or an empty string for a day not in the file. The list is cleared first.
        /// </param>
        /// <returns>The same series as the overload without flags.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="dlyText"/> or <paramref name="qualifiers"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the type is not daily precipitation or daily snow.</exception>
        /// <exception cref="FormatException">Thrown for a record too short to hold its station, date, element and days, or a record with an invalid year or month.</exception>
        public static TimeSeries ParseGHCN(string dlyText, TimeSeriesType timeSeriesType, DepthUnit unit, IList<string> qualifiers)
        {
            if (qualifiers == null) throw new ArgumentNullException(nameof(qualifiers));
            return ParseGhcnText(dlyText, timeSeriesType, unit, DateTime.Today, qualifiers);
        }

        /// <summary>
        /// Parses a GHCN-Daily station file, leaving out days after a given date.
        /// </summary>
        /// <param name="dlyText">The station file text.</param>
        /// <param name="timeSeriesType">Daily precipitation or daily snow.</param>
        /// <param name="unit">The depth unit.</param>
        /// <param name="today">The last date kept; later days of the current month are not observations yet.</param>
        /// <param name="qualifiers">Receives the flags of each ordinate, or null.</param>
        /// <returns>A regular one-day series.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="dlyText"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the type is not daily precipitation or daily snow.</exception>
        internal static TimeSeries ParseGhcnText(string dlyText, TimeSeriesType timeSeriesType, DepthUnit unit, DateTime today, IList<string>? qualifiers = null)
        {
            if (dlyText == null) throw new ArgumentNullException(nameof(dlyText));
            ValidateGhcnType(timeSeriesType);
            var qualifierByDate = qualifiers != null ? new Dictionary<DateTime, string>() : null;

            var timeSeries = new TimeSeries(TimeInterval.OneDay);
            DateTime? previousDate = null;
            bool previousSuppressCollectionChanged = timeSeries.SuppressCollectionChanged;
            var indexesByDateTime = new Dictionary<DateTime, int>();
            timeSeries.SuppressCollectionChanged = true;
            try
            {
                using var reader = new StringReader(dlyText);
                string? line;
                int lineNumber = 0;
                while ((line = reader.ReadLine()) != null)
                {
                    lineNumber++;
                    if (line.Length < 21)
                    {
                        if (line.Trim().Length == 0) continue;
                        throw GhcnRecordTooShort(lineNumber, line.Length, 21);
                    }

                    // Extract the type (PRCP, SNOW, etc.)
                    var typeString = line.Substring(17, 4);
                    if ((typeString == "PRCP" && timeSeriesType == TimeSeriesType.DailyPrecipitation) ||
                        (typeString == "SNOW" && timeSeriesType == TimeSeriesType.DailySnow))
                    {

                        // Extract year, month, and parse days
                        if (!int.TryParse(line.Substring(11, 4), NumberStyles.Integer, CultureInfo.InvariantCulture, out var year) ||
                            !int.TryParse(line.Substring(15, 2), NumberStyles.Integer, CultureInfo.InvariantCulture, out var month) ||
                            year < 1 || year > 9999 || month < 1 || month > 12)
                            throw new FormatException($"GHCN station file line {lineNumber} has an invalid year or month.");
                        int daysInMonth = DateTime.DaysInMonth(year, month);
                        int needed = 21 + (daysInMonth - 1) * 8 + 5;
                        if (line.Length < needed) throw GhcnRecordTooShort(lineNumber, line.Length, needed);

                        for (int i = 0; i < daysInMonth; i++)
                        {
                            int offset = 21 + (i * 8);
                            string strgValue = line.Substring(offset, 5).Trim();
                            var currentDate = new DateTime(year, month, i + 1);
                            if (currentDate > today) break;

                            // Fill missing dates with NaN if there is a gap
                            if (previousDate.HasValue && (currentDate - previousDate.Value).Days > 1)
                            {
                                FillMissingDates(timeSeries, previousDate.Value, currentDate, indexesByDateTime);
                            }

                            // Parse the precipitation or snow value
                            if (int.TryParse(strgValue, NumberStyles.Integer, CultureInfo.InvariantCulture, out int rawValue) && rawValue != -9999)  // Valid precipitation value
                            {
                                double convertedValue = ConvertToDesiredUnit(rawValue, unit, timeSeriesType == TimeSeriesType.DailySnow);
                                AddOrReplaceByDateTime(timeSeries, indexesByDateTime, currentDate, convertedValue);
                            }
                            else
                            {
                                // Add missing value if invalid
                                AddOrReplaceByDateTime(timeSeries, indexesByDateTime, currentDate, double.NaN);
                            }
                            if (qualifierByDate != null) qualifierByDate[currentDate] = offset + 8 <= line.Length ? line.Substring(offset + 5, 3) : "";

                            previousDate = currentDate;
                        }
                    }
                }
            }
            finally
            {
                timeSeries.SuppressCollectionChanged = previousSuppressCollectionChanged;
            }

            if (qualifiers != null) WriteQualifiers(timeSeries, qualifierByDate!, qualifiers);
            return timeSeries;
        }

        /// <summary>
        /// Creates the exception for a GHCN record that is too short to parse.
        /// </summary>
        /// <param name="lineNumber">The one-based line number.</param>
        /// <param name="length">The line length.</param>
        /// <param name="needed">The length the record needs.</param>
        /// <returns>The exception to throw.</returns>
        /// <remarks>A complete GHCN-Daily record is 269 characters; a shorter one is usually a file cut off mid-line.</remarks>
        private static FormatException GhcnRecordTooShort(int lineNumber, int length, int needed)
        {
            return new FormatException($"GHCN station file line {lineNumber} has {length} characters; the record needs at least {needed} " +
                "(a complete record has 269). The file may have been cut off.");
        }

        /// <summary>
        /// Checks that a GHCN request asks for daily precipitation or daily snow.
        /// </summary>
        /// <param name="timeSeriesType">The requested type.</param>
        /// <exception cref="ArgumentException">Thrown for any other type.</exception>
        private static void ValidateGhcnType(TimeSeriesType timeSeriesType)
        {
            if (timeSeriesType != TimeSeriesType.DailyPrecipitation && timeSeriesType != TimeSeriesType.DailySnow)
            {
                throw new ArgumentException("The time series type must be either daily precipitation or daily snow.", nameof(timeSeriesType));
            }
        }

        /// <summary>
        /// Converts raw GHCN values to the desired unit.
        /// </summary>
        /// <param name="rawValue">The raw value from the GHCN dataset.</param>
        /// <param name="unit">The unit to convert to.</param>
        /// <param name="isSnow">True for snowfall (SNOW), false for precipitation (PRCP).</param>
        /// <returns>The converted value.</returns>
        /// <remarks>
        /// The GHCN-Daily readme gives PRCP in tenths of millimeters and SNOW in millimeters, so each element
        /// is divided once into its native millimeters before the unit conversion.
        /// </remarks>
        private static double ConvertToDesiredUnit(int rawValue, DepthUnit unit, bool isSnow)
        {
            if (isSnow)
            {
                return unit switch
                {
                    DepthUnit.Millimeters => rawValue,
                    DepthUnit.Centimeters => rawValue / 10.0,
                    DepthUnit.Inches => rawValue / 25.4,
                    _ => throw new ArgumentOutOfRangeException(nameof(unit), "Unsupported depth unit.")
                };
            }
            return unit switch
            {
                DepthUnit.Millimeters => rawValue / 10.0,
                DepthUnit.Centimeters => rawValue / 100.0,
                DepthUnit.Inches => rawValue / 254.0,
                _ => throw new ArgumentOutOfRangeException(nameof(unit), "Unsupported depth unit.")
            };
        }

        /// <summary>
        /// Fills missing dates between the last available date and the current date.
        /// </summary>
        /// <param name="timeSeries">The time series to fill.</param>
        /// <param name="previousDate">The last available date.</param>
        /// <param name="currentDate">The current date to reach.</param>
        /// <param name="indexesByDateTime">Optional parsed-date lookup used when duplicate replacement is active.</param>
        /// <remarks>
        /// When a lookup is supplied, inserted gap ordinates participate in the same inline
        /// duplicate handling as provider records. Without a lookup the method keeps the
        /// append-only behavior.
        /// </remarks>
        private static void FillMissingDates(
            TimeSeries timeSeries,
            DateTime previousDate,
            DateTime currentDate,
            Dictionary<DateTime, int>? indexesByDateTime = null)
        {
            DateTime missingDate = previousDate.AddDays(1);
            while (missingDate < currentDate)
            {
                if (indexesByDateTime == null)
                {
                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(missingDate, double.NaN));
                }
                else
                {
                    AddOrReplaceByDateTime(timeSeries, indexesByDateTime, missingDate, double.NaN);
                }

                missingDate = missingDate.AddDays(1);
            }
        }

        #endregion

        #region USGS

        /// <summary>
        /// Base URL of the USGS Water Data OGC API collections.
        /// </summary>
        private const string UsgsApiBaseUrl = "https://api.waterdata.usgs.gov/ogcapi/v1/collections/";

        /// <summary>
        /// Largest page the USGS Water Data API returns for one request.
        /// </summary>
        private const int UsgsPageLimit = 50000;

        /// <summary>
        /// Number of yearly continuous-value requests run at the same time.
        /// </summary>
        /// <remarks>
        /// The continuous collection accepts at most three years per query and the service recommends
        /// one query per year. Four concurrent yearly requests cut a period-of-record download to
        /// about a quarter of the sequential time while keeping the request rate modest.
        /// </remarks>
        private const int UsgsParallelRequests = 4;

        /// <summary>
        /// Gets or sets the default USGS Water Data API key, used when <see cref="UsgsDownloadOptions.ApiKey"/> is not set.
        /// </summary>
        /// <remarks>
        /// The USGS Water Data API rate-limits anonymous requests. A free key from
        /// https://api.waterdata.usgs.gov/signup/ raises the limit. The key is sent in the <c>X-Api-Key</c>
        /// header and is never written to the raw text.
        /// </remarks>
        public static string? UsgsApiKey { get; set; }

        /// <summary>
        /// The daily statistic to download from USGS.
        /// </summary>
        public enum UsgsStatistic
        {
            /// <summary>
            /// Daily mean (USGS statistic code 00003).
            /// </summary>
            Mean,
            /// <summary>
            /// Daily maximum (USGS statistic code 00001).
            /// </summary>
            Maximum,
            /// <summary>
            /// Daily minimum (USGS statistic code 00002).
            /// </summary>
            Minimum
        }

        /// <summary>
        /// Options for a USGS download.
        /// </summary>
        public sealed class UsgsDownloadOptions
        {
            /// <summary>
            /// Gets or sets the inclusive start date, or null for the start of record.
            /// </summary>
            /// <remarks>Interpreted as a calendar date for daily values and peaks, and as UTC for instantaneous values and field measurements.</remarks>
            public DateTime? StartDate { get; set; }

            /// <summary>
            /// Gets or sets the inclusive end date, or null for the end of record.
            /// </summary>
            /// <remarks>Interpreted as a calendar date for daily values and peaks, and as UTC for instantaneous values and field measurements.</remarks>
            public DateTime? EndDate { get; set; }

            /// <summary>
            /// Gets or sets the daily statistic for daily discharge and stage. Default = mean. Ignored for other types.
            /// </summary>
            public UsgsStatistic Statistic { get; set; } = UsgsStatistic.Mean;

            /// <summary>
            /// Gets or sets the USGS Water Data API key for this download, or null to use <see cref="UsgsApiKey"/>.
            /// </summary>
            public string? ApiKey { get; set; }
        }

        /// <summary>
        /// Download time series data from USGS.
        /// </summary>
        /// <param name="siteNumber">USGS site number.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The downloaded time series and the raw provider text.</returns>
        /// <remarks>
        /// Downloads the full period of record with the default <see cref="UsgsDownloadOptions"/>. See
        /// <see cref="FromUSGS(string, TimeSeriesType, UsgsDownloadOptions, CancellationToken)"/>.
        /// </remarks>
        public static Task<(TimeSeries TimeSeries, string RawText)> FromUSGS(string siteNumber, TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge, CancellationToken cancellationToken = default)
        {
            return FromUSGS(siteNumber, timeSeriesType, new UsgsDownloadOptions(), cancellationToken);
        }

        /// <summary>
        /// Download time series data from USGS.
        /// </summary>
        /// <param name="siteNumber">USGS site number: 8 to 15 digits, for example 01646500 or 393109104464500.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="options">The date range, daily statistic and API key.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The downloaded time series and the raw provider text.</returns>
        /// <exception cref="ArgumentException">Thrown for a site number that is not 8 to 15 digits, an unsupported type, or a start date after the end date.</exception>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="options"/> is null.</exception>
        /// <exception cref="Exception">Thrown when the service returns an error or no data.</exception>
        /// <remarks>
        /// <para>
        /// Data come from the USGS Water Data API (https://api.waterdata.usgs.gov), which replaces the legacy
        /// NWIS WaterServices (scheduled to be retired on 22 February 2027). Each request asks for CSV with only the fields
        /// used, sorted by time, in pages of up to 50,000 records.
        /// </para>
        /// <list type="bullet">
        /// <item><b>Daily discharge and stage:</b> one value per day for the chosen statistic (mean by default), as a
        /// regular one-day series; days missing between records are NaN.</item>
        /// <item><b>Instantaneous discharge and stage:</b> the primary continuous series, as an irregular series of
        /// UTC dates (<see cref="DateTimeKind.Utc"/>). The record is requested one calendar year at a time, four
        /// years at once.</item>
        /// <item><b>Peak discharge and stage:</b> the annual peak record, as an irregular series. A peak with a time of
        /// day is dated in the site's local standard time, as in the official peak-flow file; a missing month or day is
        /// taken as 1.</item>
        /// <item><b>Measured discharge and stage:</b> field measurements as an irregular series in UTC; a measurement
        /// recorded without a time of day is dated at midnight UTC.</item>
        /// </list>
        /// <para>
        /// The raw text is the CSV the service returned, after comment lines (starting with #) that name the source,
        /// the request and the retrieval time.
        /// </para>
        /// </remarks>
        public static Task<(TimeSeries TimeSeries, string RawText)> FromUSGS(string siteNumber, TimeSeriesType timeSeriesType, UsgsDownloadOptions options, CancellationToken cancellationToken = default)
        {
            return FromUsgsCore(siteNumber, timeSeriesType, options, null, cancellationToken);
        }

        /// <summary>
        /// Download time series data from USGS, with the qualifier of each value.
        /// </summary>
        /// <param name="siteNumber">USGS site number: 8 to 15 digits.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="options">The date range, daily statistic and API key.</param>
        /// <param name="qualifiers">
        /// Receives one entry per ordinate of the returned series, in the same order: the service's qualifier field
        /// (for example "ESTIMATED,ICE"), or an empty string when there is none or the day was filled as missing.
        /// The list is cleared first.
        /// </param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The downloaded time series and the raw provider text.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="options"/> or <paramref name="qualifiers"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for a site number that is not 8 to 15 digits, an unsupported type, or a start date after the end date.</exception>
        /// <exception cref="Exception">Thrown when the service returns an error or no data.</exception>
        /// <remarks>
        /// The series is the same as <see cref="FromUSGS(string, TimeSeriesType, UsgsDownloadOptions, CancellationToken)"/>
        /// returns. Instantaneous requests also ask for the qualifier field, so their raw text has that column.
        /// </remarks>
        public static Task<(TimeSeries TimeSeries, string RawText)> FromUSGS(string siteNumber, TimeSeriesType timeSeriesType, UsgsDownloadOptions options, IList<string> qualifiers, CancellationToken cancellationToken = default)
        {
            if (qualifiers == null) throw new ArgumentNullException(nameof(qualifiers));
            return FromUsgsCore(siteNumber, timeSeriesType, options, qualifiers, cancellationToken);
        }

        /// <summary>
        /// Validates and runs a USGS download.
        /// </summary>
        /// <param name="siteNumber">USGS site number.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="options">The download options.</param>
        /// <param name="qualifiers">Receives the qualifier of each ordinate, or null.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The downloaded time series and the raw provider text.</returns>
        private static async Task<(TimeSeries TimeSeries, string RawText)> FromUsgsCore(string siteNumber, TimeSeriesType timeSeriesType, UsgsDownloadOptions options, IList<string>? qualifiers, CancellationToken cancellationToken)
        {
            if (options == null) throw new ArgumentNullException(nameof(options));

            // Check site number: 8 to 15 digits (surface-water sites have 8, latitude-longitude sites have 15)
            if (siteNumber == null || !Regex.IsMatch(siteNumber, @"^[0-9]{8,15}$"))
            {
                throw new ArgumentException("The USGS site number must be 8 to 15 digits long", nameof(siteNumber));
            }

            // Check time series type
            if (timeSeriesType == TimeSeriesType.DailyPrecipitation || timeSeriesType == TimeSeriesType.DailySnow)
            {
                throw new ArgumentException("The time series type cannot be daily precipitation or daily snow.", nameof(timeSeriesType));
            }

            if (options.StartDate.HasValue && options.EndDate.HasValue && options.StartDate.Value > options.EndDate.Value)
            {
                throw new ArgumentException("The start date must be on or before the end date.", nameof(options));
            }

            string parameterCode = timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.InstantaneousDischarge ||
                                   timeSeriesType == TimeSeriesType.PeakDischarge || timeSeriesType == TimeSeriesType.MeasuredDischarge ? "00060" : "00065";
            var request = new UsgsRequest("USGS-" + siteNumber, parameterCode, options.ApiKey ?? UsgsApiKey, cancellationToken) { Qualifiers = qualifiers };

            switch (timeSeriesType)
            {
                case TimeSeriesType.DailyDischarge:
                case TimeSeriesType.DailyStage:
                    return await DownloadUsgsDailyAsync(request, UsgsStatisticCode(options.Statistic), options.StartDate, options.EndDate).ConfigureAwait(false);
                case TimeSeriesType.InstantaneousDischarge:
                case TimeSeriesType.InstantaneousStage:
                    return await DownloadUsgsContinuousAsync(request, options.StartDate, options.EndDate).ConfigureAwait(false);
                case TimeSeriesType.PeakDischarge:
                case TimeSeriesType.PeakStage:
                    return await DownloadUsgsPeaksAsync(request, options.StartDate, options.EndDate).ConfigureAwait(false);
                default:
                    return await DownloadUsgsMeasurementsAsync(request, options.StartDate, options.EndDate).ConfigureAwait(false);
            }
        }

        /// <summary>
        /// The fixed parts of one USGS download.
        /// </summary>
        private sealed class UsgsRequest
        {
            /// <summary>
            /// Creates the request description.
            /// </summary>
            /// <param name="monitoringLocationId">The monitoring location id, for example USGS-01646500.</param>
            /// <param name="parameterCode">The USGS parameter code.</param>
            /// <param name="apiKey">The API key, or null.</param>
            /// <param name="cancellationToken">Token used to cancel the download.</param>
            internal UsgsRequest(string monitoringLocationId, string parameterCode, string? apiKey, CancellationToken cancellationToken)
            {
                MonitoringLocationId = monitoringLocationId;
                ParameterCode = parameterCode;
                ApiKey = string.IsNullOrWhiteSpace(apiKey) ? null : apiKey;
                CancellationToken = cancellationToken;
            }

            /// <summary>The monitoring location id.</summary>
            internal string MonitoringLocationId { get; }

            /// <summary>The USGS parameter code.</summary>
            internal string ParameterCode { get; }

            /// <summary>The API key, or null.</summary>
            internal string? ApiKey { get; }

            /// <summary>Token used to cancel the download.</summary>
            internal CancellationToken CancellationToken { get; }

            /// <summary>Receives the qualifier of each ordinate, or null when not requested.</summary>
            internal IList<string>? Qualifiers { get; set; }
        }

        /// <summary>
        /// Returns the USGS statistic code of a daily statistic.
        /// </summary>
        /// <param name="statistic">The daily statistic.</param>
        private static string UsgsStatisticCode(UsgsStatistic statistic)
        {
            return statistic == UsgsStatistic.Maximum ? "00001" : statistic == UsgsStatistic.Minimum ? "00002" : "00003";
        }

        /// <summary>
        /// Downloads daily values.
        /// </summary>
        /// <param name="request">The download description.</param>
        /// <param name="statisticCode">The USGS statistic code.</param>
        /// <param name="startDate">The inclusive start date, or null.</param>
        /// <param name="endDate">The inclusive end date, or null.</param>
        /// <returns>A regular daily series and the raw text.</returns>
        private static async Task<(TimeSeries TimeSeries, string RawText)> DownloadUsgsDailyAsync(UsgsRequest request, string statisticCode, DateTime? startDate, DateTime? endDate)
        {
            // The series metadata is read at the same time as the values. A site can have more than one daily series
            // for a parameter and statistic; then the values are requested again for the primary series only.
            var metadataTask = GetUsgsSeriesMetadataAsync(request, "Daily", statisticCode, InstantaneousProviderRequestTimeout);
            string query = $"daily/items?monitoring_location_id={request.MonitoringLocationId}&parameter_code={request.ParameterCode}" +
                           $"&statistic_id={statisticCode}&properties=time,value,approval_status,qualifier";
            var raw = UsgsRawTextHeader(request, $"daily values, parameter {request.ParameterCode}, statistic {statisticCode}");
            int records = await GetUsgsCsvAsync(request, query, startDate?.Date, endDate?.Date, true, InstantaneousProviderRequestTimeout, raw, _ => { }).ConfigureAwait(false);

            var metadata = await metadataTask.ConfigureAwait(false);
            if (metadata.Count > 1)
            {
                var primary = metadata.FirstOrDefault(m => m.Primary) ?? metadata[0];
                query = $"daily/items?time_series_id={primary.Id}&properties=time,value,approval_status,qualifier";
                raw = UsgsRawTextHeader(request, $"daily values, parameter {request.ParameterCode}, statistic {statisticCode}, time series {primary.Id}");
                records = await GetUsgsCsvAsync(request, query, startDate?.Date, endDate?.Date, true, InstantaneousProviderRequestTimeout, raw, _ => { }).ConfigureAwait(false);
            }
            if (records == 0) throw new Exception(NoDataMessage);

            string rawText = raw.ToString();
            return (ParseUsgsText(rawText, TimeSeriesType.DailyDischarge, null, null, null, request.Qualifiers), rawText);
        }

        /// <summary>
        /// Downloads continuous (instantaneous) values.
        /// </summary>
        /// <param name="request">The download description.</param>
        /// <param name="startDate">The inclusive start, as UTC, or null.</param>
        /// <param name="endDate">The inclusive end, as UTC, or null.</param>
        /// <returns>An irregular UTC series and the raw text.</returns>
        /// <remarks>
        /// The period of record of the primary continuous series is read from the time-series metadata, then the
        /// record is requested one calendar year at a time with up to <see cref="UsgsParallelRequests"/> requests
        /// running at once. The yearly pieces are joined in time order.
        /// </remarks>
        private static async Task<(TimeSeries TimeSeries, string RawText)> DownloadUsgsContinuousAsync(UsgsRequest request, DateTime? startDate, DateTime? endDate)
        {
            var metadata = await GetUsgsSeriesMetadataAsync(request, "Points", null).ConfigureAwait(false);
            var primary = metadata.FirstOrDefault(m => m.Primary) ?? metadata.FirstOrDefault();
            if (primary == null || primary.Begin == null || primary.End == null) throw new Exception(NoDataMessage);

            DateTime begin = startDate.HasValue && ToUtc(startDate.Value) > primary.Begin.Value ? ToUtc(startDate.Value) : primary.Begin.Value;
            DateTime end = endDate.HasValue && ToUtc(endDate.Value) < primary.End.Value ? ToUtc(endDate.Value) : primary.End.Value;
            var raw = UsgsRawTextHeader(request, $"continuous values, parameter {request.ParameterCode}, time series {primary.Id}");
            var timeSeries = new TimeSeries(TimeInterval.Irregular);
            if (begin > end)
            {
                request.Qualifiers?.Clear();
                return (timeSeries, raw.ToString());
            }

            // One window per calendar year, clipped to the requested range.
            var windows = new List<(DateTime From, DateTime To)>();
            for (int year = begin.Year; year <= end.Year; year++)
            {
                var from = new DateTime(year, 1, 1, 0, 0, 0, DateTimeKind.Utc);
                var to = new DateTime(year, 12, 31, 23, 59, 59, DateTimeKind.Utc);
                windows.Add((from < begin ? begin : from, to > end ? end : to));
            }

            string query = $"continuous/items?time_series_id={primary.Id}&properties=time,value" + (request.Qualifiers != null ? ",qualifier" : "");
            var texts = new StringBuilder[windows.Count];
            using (var throttle = new SemaphoreSlim(UsgsParallelRequests))
            {
                var tasks = new Task[windows.Count];
                for (int i = 0; i < windows.Count; i++)
                {
                    int k = i;
                    tasks[k] = Task.Run(async () =>
                    {
                        await throttle.WaitAsync(request.CancellationToken).ConfigureAwait(false);
                        try
                        {
                            texts[k] = new StringBuilder();
                            await GetUsgsCsvAsync(request, query, windows[k].From, windows[k].To, false, InstantaneousProviderRequestTimeout, texts[k], _ => { }).ConfigureAwait(false);
                        }
                        finally
                        {
                            throttle.Release();
                        }
                    }, request.CancellationToken);
                }
                await Task.WhenAll(tasks).ConfigureAwait(false);
            }

            // The yearly pieces are joined in time order; each piece keeps its own CSV header line.
            for (int k = 0; k < texts.Length; k++) raw.Append(texts[k]);
            string rawText = raw.ToString();
            return (ParseUsgsText(rawText, TimeSeriesType.InstantaneousDischarge, null, null, null, request.Qualifiers), rawText);
        }

        /// <summary>
        /// Downloads the annual peak record.
        /// </summary>
        /// <param name="request">The download description.</param>
        /// <param name="startDate">The inclusive start date, or null.</param>
        /// <param name="endDate">The inclusive end date, or null.</param>
        /// <returns>The peak series and the raw text.</returns>
        /// <remarks>
        /// The service dates a peak that has a time of day in UTC. The official peak-flow file dates it in the
        /// site's local standard time, so the date is converted back with the standard offset of the site's time
        /// zone, read from the monitoring-location record at the same time as the peaks.
        /// </remarks>
        private static async Task<(TimeSeries TimeSeries, string RawText)> DownloadUsgsPeaksAsync(UsgsRequest request, DateTime? startDate, DateTime? endDate)
        {
            string query = $"peaks/items?monitoring_location_id={request.MonitoringLocationId}&parameter_code={request.ParameterCode}" +
                           "&properties=time,value,year,month,day,time_of_day,qualifier";
            var raw = UsgsRawTextHeader(request, $"annual peaks, parameter {request.ParameterCode}");
            var zoneTask = GetUsgsStandardOffsetAsync(request);
            bool timed = false;
            await GetUsgsCsvAsync(request, query, null, null, false, ProviderRequestTimeout, raw, row => timed |= TryGetUsgsPeakParts(row, out _, out _, out _, out TimeSpan? t) && t.HasValue).ConfigureAwait(false);

            // The site's time zone is needed only when a peak has a time of day.
            TimeSpan? offset = null;
            if (timed) offset = await zoneTask.ConfigureAwait(false);
            else ObserveTask(zoneTask);
            string rawText = raw.ToString();
            return (ParseUsgsText(rawText, TimeSeriesType.PeakDischarge, offset, startDate, endDate, request.Qualifiers), rawText);
        }

        /// <summary>
        /// Reads the date parts of one annual peak record.
        /// </summary>
        /// <param name="row">The record.</param>
        /// <param name="year">The peak year.</param>
        /// <param name="month">The month, or 0 when missing or invalid.</param>
        /// <param name="day">The day, or 0 when missing or invalid.</param>
        /// <param name="timeOfDay">The UTC time of day when the record has a month, a day and a time; otherwise null.</param>
        /// <returns>False when the record has no value or no year and is skipped.</returns>
        private static bool TryGetUsgsPeakParts(UsgsCsvRow row, out int year, out int month, out int day, out TimeSpan? timeOfDay)
        {
            month = 0; day = 0; timeOfDay = null;
            if (string.IsNullOrWhiteSpace(row.Value)) { year = 0; return false; }
            if (!int.TryParse(row.Year, NumberStyles.Integer, CultureInfo.InvariantCulture, out year)) return false;
            month = int.TryParse(row.Month, NumberStyles.Integer, CultureInfo.InvariantCulture, out int m) && m >= 1 && m <= 12 ? m : 0;
            day = int.TryParse(row.Day, NumberStyles.Integer, CultureInfo.InvariantCulture, out int d) && d >= 1 && d <= 31 ? d : 0;
            if (month > 0 && day > 0 && !string.IsNullOrWhiteSpace(row.TimeOfDay) &&
                TimeSpan.TryParse(row.TimeOfDay!.Split('+', '-')[0], CultureInfo.InvariantCulture, out TimeSpan t))
            {
                timeOfDay = t;
            }
            return true;
        }

        /// <summary>
        /// Downloads field measurements.
        /// </summary>
        /// <param name="request">The download description.</param>
        /// <param name="startDate">The inclusive start, as UTC, or null.</param>
        /// <param name="endDate">The inclusive end, as UTC, or null.</param>
        /// <returns>An irregular series in time order and the raw text.</returns>
        private static async Task<(TimeSeries TimeSeries, string RawText)> DownloadUsgsMeasurementsAsync(UsgsRequest request, DateTime? startDate, DateTime? endDate)
        {
            string query = $"field-measurements/items?monitoring_location_id={request.MonitoringLocationId}&parameter_code={request.ParameterCode}" +
                           "&properties=time,time_of_day,value,approval_status,qualifier";
            var raw = UsgsRawTextHeader(request, $"field measurements, parameter {request.ParameterCode}");
            await GetUsgsCsvAsync(request, query, startDate.HasValue ? ToUtc(startDate.Value) : (DateTime?)null,
                endDate.HasValue ? ToUtc(endDate.Value) : (DateTime?)null, false, ProviderRequestTimeout, raw, _ => { }, sorted: false).ConfigureAwait(false);
            string rawText = raw.ToString();
            return (ParseUsgsText(rawText, TimeSeriesType.MeasuredDischarge, null, null, null, request.Qualifiers), rawText);
        }

        /// <summary>
        /// Parses USGS Water Data API CSV text into a time series.
        /// </summary>
        /// <param name="rawText">
        /// The CSV text of one or more responses, as <see cref="FromUSGS(string, TimeSeriesType, UsgsDownloadOptions, CancellationToken)"/>
        /// returns it: lines starting with # are skipped, and a repeated header line starts a new piece.
        /// </param>
        /// <param name="timeSeriesType">
        /// The time series type the text holds. Discharge and stage parse alike; daily, instantaneous, peak and measured
        /// types differ.
        /// </param>
        /// <param name="peakStandardOffset">
        /// For peak types only: the standard UTC offset of the site's time zone (for example -5 hours for EST), used to
        /// date a peak that has a time of day in local standard time. Not needed when no peak has a time of day.
        /// </param>
        /// <returns>The parsed time series.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="rawText"/> is null.</exception>
        /// <exception cref="ArgumentException">
        /// Thrown for daily precipitation or daily snow, or for a peak with a time of day when
        /// <paramref name="peakStandardOffset"/> is null.
        /// </exception>
        /// <remarks>
        /// <para>
        /// This is the parser <c>FromUSGS</c> applies to the downloaded text, so saved raw text gives the same series
        /// offline. Fields are located by header name, and quoted fields (RFC 4180) are supported.
        /// </para>
        /// <list type="bullet">
        /// <item><b>Daily:</b> a regular one-day series; days missing between records, and empty values, are NaN.</item>
        /// <item><b>Instantaneous:</b> an irregular series of UTC dates (<see cref="DateTimeKind.Utc"/>).</item>
        /// <item><b>Peak:</b> an irregular series of one value per peak; a missing month or day is taken as 1.</item>
        /// <item><b>Measured:</b> an irregular UTC series sorted by time; a record without a time of day is dated at midnight UTC.</item>
        /// </list>
        /// <para>A later record at the same date and time replaces an earlier one.</para>
        /// </remarks>
        public static TimeSeries ParseUSGS(string rawText, TimeSeriesType timeSeriesType, TimeSpan? peakStandardOffset = null)
        {
            if (rawText == null) throw new ArgumentNullException(nameof(rawText));
            if (timeSeriesType == TimeSeriesType.DailyPrecipitation || timeSeriesType == TimeSeriesType.DailySnow)
            {
                throw new ArgumentException("The time series type cannot be daily precipitation or daily snow.", nameof(timeSeriesType));
            }
            return ParseUsgsText(rawText, timeSeriesType, peakStandardOffset, null, null, null);
        }

        /// <summary>
        /// Parses USGS Water Data API CSV text into a time series, with the qualifier of each value.
        /// </summary>
        /// <param name="rawText">The CSV text, as for <see cref="ParseUSGS(string, TimeSeriesType, TimeSpan?)"/>.</param>
        /// <param name="timeSeriesType">The time series type the text holds.</param>
        /// <param name="peakStandardOffset">For timed peaks: the standard UTC offset of the site's time zone, or null.</param>
        /// <param name="qualifiers">
        /// Receives one entry per ordinate of the returned series, in the same order: the qualifier field, or an empty
        /// string when the text has none or the day was filled as missing. The list is cleared first.
        /// </param>
        /// <returns>The parsed time series, the same as the overload without qualifiers returns.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="rawText"/> or <paramref name="qualifiers"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for daily precipitation or daily snow, or for a timed peak without an offset.</exception>
        public static TimeSeries ParseUSGS(string rawText, TimeSeriesType timeSeriesType, TimeSpan? peakStandardOffset, IList<string> qualifiers)
        {
            if (qualifiers == null) throw new ArgumentNullException(nameof(qualifiers));
            if (rawText == null) throw new ArgumentNullException(nameof(rawText));
            if (timeSeriesType == TimeSeriesType.DailyPrecipitation || timeSeriesType == TimeSeriesType.DailySnow)
            {
                throw new ArgumentException("The time series type cannot be daily precipitation or daily snow.", nameof(timeSeriesType));
            }
            return ParseUsgsText(rawText, timeSeriesType, peakStandardOffset, null, null, qualifiers);
        }

        /// <summary>
        /// Parses USGS CSV text, keeping only peaks dated within an optional range.
        /// </summary>
        /// <param name="rawText">The CSV text.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="peakStandardOffset">The site's standard UTC offset, for timed peaks.</param>
        /// <param name="startDate">For peaks: the inclusive start date, or null.</param>
        /// <param name="endDate">For peaks: the inclusive end date, or null.</param>
        /// <param name="qualifiers">Receives the qualifier of each ordinate, or null.</param>
        /// <returns>The parsed time series.</returns>
        /// <exception cref="ArgumentException">Thrown for a timed peak when <paramref name="peakStandardOffset"/> is null.</exception>
        private static TimeSeries ParseUsgsText(string rawText, TimeSeriesType timeSeriesType, TimeSpan? peakStandardOffset, DateTime? startDate, DateTime? endDate, IList<string>? qualifiers)
        {
            var qualifierByDate = qualifiers != null ? new Dictionary<DateTime, string>() : null;
            var timeSeries = ParseUsgsRecords(rawText, timeSeriesType, peakStandardOffset, startDate, endDate, qualifierByDate);
            if (qualifiers != null) WriteQualifiers(timeSeries, qualifierByDate!, qualifiers);
            return timeSeries;
        }

        /// <summary>
        /// Parses the records of USGS CSV text by type.
        /// </summary>
        /// <param name="rawText">The CSV text.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="peakStandardOffset">The site's standard UTC offset, for timed peaks.</param>
        /// <param name="startDate">For peaks: the inclusive start date, or null.</param>
        /// <param name="endDate">For peaks: the inclusive end date, or null.</param>
        /// <param name="qualifierByDate">Receives the qualifier field of each record by date, or null.</param>
        /// <returns>The parsed time series.</returns>
        private static TimeSeries ParseUsgsRecords(string rawText, TimeSeriesType timeSeriesType, TimeSpan? peakStandardOffset, DateTime? startDate, DateTime? endDate, Dictionary<DateTime, string>? qualifierByDate)
        {
            var indexesByDateTime = new Dictionary<DateTime, int>();
            switch (timeSeriesType)
            {
                case TimeSeriesType.DailyDischarge:
                case TimeSeriesType.DailyStage:
                    {
                        var timeSeries = new TimeSeries();
                        timeSeries.SuppressCollectionChanged = true;
                        ParseUsgsCsv(rawText, row =>
                        {
                            if (!DateTime.TryParseExact(row.Time, "yyyy-MM-dd", CultureInfo.InvariantCulture, DateTimeStyles.None, out DateTime index)) return;

                            // Fill in missing days
                            if (timeSeries.Count > 0 && index != TimeSeries.AddTimeInterval(timeSeries.Last().Index, TimeInterval.OneDay))
                            {
                                while (timeSeries.Last().Index < TimeSeries.SubtractTimeInterval(index, TimeInterval.OneDay))
                                {
                                    DateTime missingIndex = TimeSeries.AddTimeInterval(timeSeries.Last().Index, TimeInterval.OneDay);
                                    AddOrReplaceByDateTime(timeSeries, indexesByDateTime, missingIndex, double.NaN);
                                }
                            }
                            AddOrReplaceByDateTime(timeSeries, indexesByDateTime, index, ParseUsgsValue(row.Value));
                            if (qualifierByDate != null) qualifierByDate[index] = row.Qualifier ?? "";
                        });
                        timeSeries.SuppressCollectionChanged = false;
                        return timeSeries;
                    }
                case TimeSeriesType.InstantaneousDischarge:
                case TimeSeriesType.InstantaneousStage:
                    {
                        var timeSeries = new TimeSeries(TimeInterval.Irregular);
                        timeSeries.SuppressCollectionChanged = true;
                        ParseUsgsCsv(rawText, row =>
                        {
                            if (!TryParseUsgsTimestamp(row.Time, out DateTime time)) return;
                            AddOrReplaceByDateTime(timeSeries, indexesByDateTime, time, ParseUsgsValue(row.Value));
                            if (qualifierByDate != null) qualifierByDate[time] = row.Qualifier ?? "";
                        });
                        timeSeries.SuppressCollectionChanged = false;
                        return timeSeries;
                    }
                case TimeSeriesType.PeakDischarge:
                case TimeSeriesType.PeakStage:
                    {
                        // Peaks fall on irregular dates, one or more per water year, as Canadian peaks do.
                        var timeSeries = new TimeSeries(TimeInterval.Irregular);
                        ParseUsgsCsv(rawText, row =>
                        {
                            if (!TryGetUsgsPeakParts(row, out int year, out int month, out int day, out TimeSpan? timeOfDay)) return;
                            DateTime index;
                            if (timeOfDay.HasValue)
                            {
                                if (!peakStandardOffset.HasValue)
                                    throw new ArgumentException("A peak with a time of day needs the site's standard UTC offset.", nameof(peakStandardOffset));
                                index = new DateTime(year, month, day).Add(timeOfDay.Value).Add(peakStandardOffset.Value).Date;
                            }
                            else
                            {
                                // A missing month or day is taken as 1, as in the USGS peak-flow file convention.
                                index = new DateTime(year, month > 0 ? month : 1, day > 0 ? day : 1);
                            }

                            if (startDate.HasValue && index < startDate.Value.Date) return;
                            if (endDate.HasValue && index > endDate.Value.Date) return;
                            if (!double.TryParse(row.Value, NumberStyles.Float, CultureInfo.InvariantCulture, out double value)) return;
                            AddOrReplaceByDateTime(timeSeries, indexesByDateTime, index, value);
                            if (qualifierByDate != null) qualifierByDate[index] = row.Qualifier ?? "";
                        });
                        return timeSeries;
                    }
                default:
                    {
                        // Records are read in service order, so a later record at the same time replaces an earlier one, and
                        // the series is sorted by time at the end.
                        var timeSeries = new TimeSeries(TimeInterval.Irregular);
                        ParseUsgsCsv(rawText, row =>
                        {
                            // The service gives the UTC date and time of day separately. A measurement without a time of day is
                            // dated at midnight UTC, independent of the machine's time zone.
                            string text = string.IsNullOrWhiteSpace(row.TimeOfDay) ? row.Time : row.Time + "T" + row.TimeOfDay;
                            if (!DateTime.TryParse(text, CultureInfo.InvariantCulture, DateTimeStyles.AdjustToUniversal | DateTimeStyles.AssumeUniversal, out DateTime index)) return;
                            if (string.IsNullOrWhiteSpace(row.Value) ||
                                !double.TryParse(row.Value, NumberStyles.Float, CultureInfo.InvariantCulture, out double value)) return;
                            AddOrReplaceByDateTime(timeSeries, indexesByDateTime, index, value);
                            if (qualifierByDate != null) qualifierByDate[index] = row.Qualifier ?? "";
                        });
                        timeSeries.SortByTime();
                        return timeSeries;
                    }
            }
        }

        /// <summary>
        /// One CSV record from the USGS Water Data API, holding the fields the readers use.
        /// </summary>
        private sealed class UsgsCsvRow
        {
            /// <summary>The time or date.</summary>
            internal string Time = "";
            /// <summary>The value text.</summary>
            internal string Value = "";
            /// <summary>The peak year, or null.</summary>
            internal string? Year;
            /// <summary>The peak month, or null.</summary>
            internal string? Month;
            /// <summary>The peak day, or null.</summary>
            internal string? Day;
            /// <summary>The peak time of day, or null.</summary>
            internal string? TimeOfDay;
            /// <summary>The qualifier field, or null.</summary>
            internal string? Qualifier;
        }

        /// <summary>
        /// Requests a sorted CSV query and pages forward in time until every record is read.
        /// </summary>
        /// <param name="request">The download description.</param>
        /// <param name="query">The collection path and filters, without paging, sorting, format or time arguments.</param>
        /// <param name="from">The inclusive start, or null.</param>
        /// <param name="to">The inclusive end, or null.</param>
        /// <param name="dateOnly">True when the records are dated by day.</param>
        /// <param name="timeout">The per-request timeout.</param>
        /// <param name="raw">Receives the raw CSV text.</param>
        /// <param name="onRow">Called for each record. The record object is reused; copy it to keep it.</param>
        /// <param name="sorted">True to sort by time; false to keep the service order.</param>
        /// <returns>The number of records read.</returns>
        /// <remarks>
        /// The service returns only the first page of a sorted query, so a full sorted page is followed by a new
        /// query that starts just after the last record read. Unsorted queries page by offset.
        /// </remarks>
        private static async Task<int> GetUsgsCsvAsync(UsgsRequest request, string query, DateTime? from, DateTime? to, bool dateOnly, TimeSpan timeout, StringBuilder raw, Action<UsgsCsvRow> onRow, bool sorted = true)
        {
            int total = 0;
            bool first = true;
            while (true)
            {
                string url = UsgsApiBaseUrl + query + "&skipGeometry=true" + (sorted ? "&sortby=time" : "") + $"&limit={UsgsPageLimit}&f=csv";
                if (!sorted && total > 0)
                    url += "&offset=" + total.ToString(CultureInfo.InvariantCulture);
                if (from.HasValue || to.HasValue)
                    url += "&datetime=" + FormatUsgsTime(from, dateOnly) + "/" + FormatUsgsTime(to, dateOnly);

                var (status, body) = await GetWithRetryAsync(_decompressClient, url, maxAttempts: 4, requestTimeout: timeout,
                    cancellationToken: request.CancellationToken, apiKey: request.ApiKey, isComplete: IsCompleteUsgsCsv).ConfigureAwait(false);
                if ((int)status >= 400) throw UsgsRequestException("USGS data", status, body);
                raw.Append("# Request: ").Append(url).Append('\n');
                AppendUsgsCsvText(raw, body, first);
                first = false;

                string lastTime = "";
                int count = ParseUsgsCsv(body, row => { lastTime = row.Time; onRow(row); });
                total += count;
                if (count < UsgsPageLimit) return total;
                if (!sorted) continue;

                // Page forward from just after the last record read.
                if (dateOnly)
                {
                    if (!DateTime.TryParseExact(lastTime, "yyyy-MM-dd", CultureInfo.InvariantCulture, DateTimeStyles.None, out DateTime lastDate)) return total;
                    from = lastDate.AddDays(1);
                }
                else
                {
                    if (!TryParseUsgsTimestamp(lastTime, out DateTime last)) return total;
                    from = last.AddSeconds(1);
                }
                if (to.HasValue && from > to) return total;
            }
        }

        /// <summary>
        /// Creates the exception for a failed USGS request.
        /// </summary>
        /// <param name="what">What was requested.</param>
        /// <param name="status">The HTTP status.</param>
        /// <param name="body">The response body.</param>
        /// <returns>The exception to throw.</returns>
        /// <remarks>
        /// HTTP 429 means the per-hour request limit is spent (the USGS documents 1,000 requests per hour for
        /// anonymous use); the message says how to raise it with an API key.
        /// </remarks>
        private static Exception UsgsRequestException(string what, HttpStatusCode status, string body)
        {
            string snippet = body.Length > 500 ? body.Substring(0, 500) + "…" : body;
            if ((int)status == 429)
            {
                return new Exception($"Failed to retrieve {what}: the USGS Water Data API rate limit was exceeded. " +
                    "Set TimeSeriesDownload.UsgsApiKey (or UsgsDownloadOptions.ApiKey) to a free key from https://api.waterdata.usgs.gov/signup/, " +
                    $"or try again later. Status=429. Body={snippet}");
            }
            return new Exception($"Failed to retrieve {what}. Status={(int)status} ({status}). Body={snippet}");
        }

        /// <summary>
        /// Formats one end of a datetime filter, or ".." for an open end.
        /// </summary>
        /// <param name="time">The time, or null.</param>
        /// <param name="dateOnly">True to format a calendar date.</param>
        private static string FormatUsgsTime(DateTime? time, bool dateOnly)
        {
            if (!time.HasValue) return "..";
            return dateOnly
                ? time.Value.ToString("yyyy-MM-dd", CultureInfo.InvariantCulture)
                : ToUtc(time.Value).ToString("yyyy-MM-dd'T'HH:mm:ss'Z'", CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Returns a date as UTC, treating an unspecified kind as UTC.
        /// </summary>
        /// <param name="time">The date.</param>
        private static DateTime ToUtc(DateTime time)
        {
            if (time.Kind == DateTimeKind.Local) return time.ToUniversalTime();
            return DateTime.SpecifyKind(time, DateTimeKind.Utc);
        }

        /// <summary>
        /// Parses a USGS timestamp such as "2024-03-10 06:15:00+00:00" to a UTC date.
        /// </summary>
        /// <param name="text">The timestamp.</param>
        /// <param name="time">The UTC date.</param>
        /// <returns>True when the text is a timestamp.</returns>
        private static bool TryParseUsgsTimestamp(string text, out DateTime time)
        {
            return DateTime.TryParse(text, CultureInfo.InvariantCulture, DateTimeStyles.AdjustToUniversal | DateTimeStyles.AssumeUniversal, out time);
        }

        /// <summary>
        /// Parses a USGS value; an empty or non-numeric value is NaN.
        /// </summary>
        /// <param name="text">The value text.</param>
        private static double ParseUsgsValue(string text)
        {
            return !string.IsNullOrWhiteSpace(text) && double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out double value) ? value : double.NaN;
        }

        /// <summary>
        /// Parses the records of a USGS CSV response.
        /// </summary>
        /// <param name="body">The CSV text, with a header line.</param>
        /// <param name="onRow">Called for each record; the record object is reused.</param>
        /// <returns>The number of records in the response.</returns>
        /// <remarks>
        /// Fields are located by header name. Quoted fields (RFC 4180) are supported; the qualifier field, for
        /// example, can hold a quoted list. Lines starting with # are comments. The first record is the header, and a
        /// later record whose fields include both "time" and "value" is a repeated header that starts a new piece of
        /// joined text.
        /// </remarks>
        private static int ParseUsgsCsv(string body, Action<UsgsCsvRow> onRow)
        {
            int count = 0;
            var row = new UsgsCsvRow();
            int position = 0;
            int[]? map = null;
            var fields = new List<string>();
            while (position < body.Length)
            {
                if (body[position] == '#')
                {
                    int lineEnd = body.IndexOf('\n', position);
                    position = lineEnd < 0 ? body.Length : lineEnd + 1;
                    continue;
                }
                position = ReadCsvRecord(body, position, fields);
                if (fields.Count == 0 || (fields.Count == 1 && fields[0].Length == 0)) continue;
                if (map == null || (fields.Contains("time") && fields.Contains("value")))
                {
                    map = new int[8];
                    string[] names = { "", "time", "value", "year", "month", "day", "time_of_day", "qualifier" };
                    for (int j = 0; j < names.Length; j++) map[j] = fields.IndexOf(names[j]);
                    continue;
                }
                count++;
                row.Time = Field(fields, map[1]) ?? "";
                row.Value = Field(fields, map[2]) ?? "";
                row.Year = Field(fields, map[3]);
                row.Month = Field(fields, map[4]);
                row.Day = Field(fields, map[5]);
                row.TimeOfDay = Field(fields, map[6]);
                row.Qualifier = Field(fields, map[7]);
                onRow(row);
            }
            return count;
        }

        /// <summary>
        /// Returns a field by position, or null when the position is absent.
        /// </summary>
        /// <param name="fields">The record fields.</param>
        /// <param name="index">The field position, or -1.</param>
        private static string? Field(List<string> fields, int index)
        {
            return index >= 0 && index < fields.Count ? fields[index] : null;
        }

        /// <summary>
        /// Reads one CSV record starting at a position.
        /// </summary>
        /// <param name="text">The CSV text.</param>
        /// <param name="position">The position of the record.</param>
        /// <param name="fields">Receives the fields of the record.</param>
        /// <returns>The position after the record and its line break.</returns>
        private static int ReadCsvRecord(string text, int position, List<string> fields)
        {
            fields.Clear();
            var quoted = new StringBuilder();
            int i = position;
            while (true)
            {
                if (i < text.Length && text[i] == '"')
                {
                    // Quoted field: "" is an escaped quote.
                    quoted.Clear();
                    i++;
                    while (i < text.Length)
                    {
                        if (text[i] == '"')
                        {
                            if (i + 1 < text.Length && text[i + 1] == '"') { quoted.Append('"'); i += 2; continue; }
                            i++;
                            break;
                        }
                        quoted.Append(text[i++]);
                    }
                    fields.Add(quoted.ToString());
                }
                else
                {
                    int start = i;
                    while (i < text.Length && text[i] != ',' && text[i] != '\n' && text[i] != '\r') i++;
                    fields.Add(text.Substring(start, i - start));
                }

                if (i < text.Length && text[i] == ',') { i++; continue; }
                if (i < text.Length && text[i] == '\r') i++;
                if (i < text.Length && text[i] == '\n') i++;
                return i;
            }
        }

        /// <summary>
        /// Appends CSV text to the raw text, without the header line except for the first piece.
        /// </summary>
        /// <param name="raw">The raw text.</param>
        /// <param name="text">The CSV text of one response.</param>
        /// <param name="includeHeader">True to keep the header line.</param>
        private static void AppendUsgsCsvText(StringBuilder raw, string text, bool includeHeader)
        {
            int start = 0;
            if (!includeHeader)
            {
                int lineEnd = text.IndexOf('\n');
                start = lineEnd < 0 ? text.Length : lineEnd + 1;
            }
            raw.Append(text, start, text.Length - start);
            if (text.Length > start && text[text.Length - 1] != '\n') raw.Append('\n');
        }

        /// <summary>
        /// Starts the raw text with comment lines naming the source and retrieval.
        /// </summary>
        /// <param name="request">The download description.</param>
        /// <param name="description">What was requested.</param>
        private static StringBuilder UsgsRawTextHeader(UsgsRequest request, string description)
        {
            var raw = new StringBuilder();
            raw.Append("# U.S. Geological Survey, Water Data API (https://api.waterdata.usgs.gov)\n");
            raw.Append("# Monitoring location: ").Append(request.MonitoringLocationId).Append('\n');
            raw.Append("# Data: ").Append(description).Append('\n');
            raw.Append("# Retrieved: ").Append(DateTime.UtcNow.ToString("yyyy-MM-dd'T'HH:mm:ss'Z'", CultureInfo.InvariantCulture)).Append('\n');
            raw.Append("# Recent data are provisional and subject to revision (see approval_status where listed).\n");
            return raw;
        }

        /// <summary>
        /// Time-series metadata of one USGS series.
        /// </summary>
        private sealed class UsgsSeriesMetadata
        {
            /// <summary>The time series id.</summary>
            internal string Id = "";
            /// <summary>True for the primary series.</summary>
            internal bool Primary;
            /// <summary>The first record time, UTC.</summary>
            internal DateTime? Begin;
            /// <summary>The last record time, UTC.</summary>
            internal DateTime? End;
        }

        /// <summary>
        /// Reads the time-series metadata of a site and parameter.
        /// </summary>
        /// <param name="request">The download description.</param>
        /// <param name="computationPeriod">The computation period, for example Points or Daily.</param>
        /// <param name="statisticCode">The statistic code, or null for any.</param>
        /// <param name="timeout">The per-request timeout, or null for the default provider timeout.</param>
        /// <remarks>
        /// The daily download reads the metadata while the values are requested, so it gives the metadata request the
        /// values request's timeout; a slow metadata response then does not fail a download whose values arrive.
        /// </remarks>
        private static async Task<List<UsgsSeriesMetadata>> GetUsgsSeriesMetadataAsync(UsgsRequest request, string computationPeriod, string? statisticCode, TimeSpan? timeout = null)
        {
            string url = UsgsApiBaseUrl + $"time-series-metadata/items?monitoring_location_id={request.MonitoringLocationId}&parameter_code={request.ParameterCode}" +
                         $"&computation_period_identifier={computationPeriod}" + (statisticCode != null ? $"&statistic_id={statisticCode}" : "") +
                         "&properties=id,begin,end,primary&skipGeometry=true&limit=100&f=csv";
            var (status, body) = await GetWithRetryAsync(_decompressClient, url, maxAttempts: 4, requestTimeout: timeout, cancellationToken: request.CancellationToken,
                apiKey: request.ApiKey, isComplete: IsCompleteUsgsCsv).ConfigureAwait(false);
            if ((int)status >= 400) throw UsgsRequestException("USGS time-series metadata", status, body);

            var result = new List<UsgsSeriesMetadata>();
            var fields = new List<string>();
            int position = 0;
            int idCol = -1, beginCol = -1, endCol = -1, primaryCol = -1;
            bool header = true;
            while (position < body.Length)
            {
                position = ReadCsvRecord(body, position, fields);
                if (fields.Count <= 1) continue;
                if (header)
                {
                    idCol = fields.IndexOf("id"); beginCol = fields.IndexOf("begin"); endCol = fields.IndexOf("end"); primaryCol = fields.IndexOf("primary");
                    header = false;
                    continue;
                }
                var item = new UsgsSeriesMetadata
                {
                    Id = Field(fields, idCol) ?? "",
                    Primary = string.Equals(Field(fields, primaryCol), "Primary", StringComparison.OrdinalIgnoreCase)
                };
                if (TryParseUsgsTimestamp(Field(fields, beginCol) ?? "", out DateTime b)) item.Begin = b;
                if (TryParseUsgsTimestamp(Field(fields, endCol) ?? "", out DateTime e)) item.End = e;
                if (item.Id.Length > 0) result.Add(item);
            }
            return result;
        }

        /// <summary>
        /// Reads the standard UTC offset of a site's time zone.
        /// </summary>
        /// <param name="request">The download description.</param>
        /// <returns>The standard offset, for example -5 hours for EST.</returns>
        /// <exception cref="Exception">Thrown when the site's time zone is unknown.</exception>
        private static async Task<TimeSpan> GetUsgsStandardOffsetAsync(UsgsRequest request)
        {
            string url = UsgsApiBaseUrl + $"monitoring-locations/items/{request.MonitoringLocationId}?properties=time_zone_abbreviation&f=json";
            var (status, body) = await GetWithRetryAsync(_decompressClient, url, maxAttempts: 4, cancellationToken: request.CancellationToken, apiKey: request.ApiKey).ConfigureAwait(false);
            if ((int)status >= 400) throw UsgsRequestException("the USGS monitoring location", status, body);
            using var doc = System.Text.Json.JsonDocument.Parse(body);
            string? zone = doc.RootElement.TryGetProperty("properties", out var props) && props.TryGetProperty("time_zone_abbreviation", out var z) ? z.GetString() : null;
            if (zone != null && UsgsStandardOffsets.TryGetValue(zone.Trim(), out double hours)) return TimeSpan.FromHours(hours);
            throw new Exception($"Unknown time zone '{zone}' for {request.MonitoringLocationId}.");
        }

        /// <summary>
        /// Standard UTC offsets, in hours, of the time zone abbreviations used by USGS monitoring locations.
        /// </summary>
        private static readonly Dictionary<string, double> UsgsStandardOffsets = new Dictionary<string, double>(StringComparer.OrdinalIgnoreCase)
        {
            { "UTC", 0 }, { "GMT", 0 }, { "ZULU", 0 },
            { "AST", -4 }, { "ADT", -4 }, { "EST", -5 }, { "EDT", -5 }, { "CST", -6 }, { "CDT", -6 },
            { "MST", -7 }, { "MDT", -7 }, { "PST", -8 }, { "PDT", -8 }, { "AKST", -9 }, { "AKDT", -9 },
            { "HST", -10 }, { "SST", -11 }, { "ChST", 10 }, { "GST", 10 }, { "WAKT", 12 }
        };

        /// <summary>
        /// Observes a task that is no longer awaited so a later fault is not unobserved.
        /// </summary>
        /// <param name="task">The task.</param>
        private static void ObserveTask(Task task)
        {
            task.ContinueWith(t => _ = t.Exception, TaskContinuationOptions.OnlyOnFaulted | TaskContinuationOptions.ExecuteSynchronously);
        }

        #endregion

        #region Canadian Hydrometric Monitoring Network (CHMN)

        /// <summary>
        /// Download time series data from the Water Survey of Canada
        /// (Environment and Climate Change Canada Hydrometric Monitoring Network).
        /// </summary>
        /// <param name="stationNumber">WSC 7-character station id, e.g., "08LG010".</param>
        /// <param name="timeSeriesType">
        ///     DailyDischarge / DailyStage → daily_data endpoint (parameters[]=flow/level).
        ///     InstantaneousDischarge / InstantaneousStage → real_time_data endpoint (parameters[]=46/47), dated in UTC (DateTimeKind.Utc).
        ///     PeakDischarge / PeakStage → peak_data endpoint (parameters[]=flow/level).
        /// </param>
        /// <param name="dischargeUnit">Desired discharge unit for discharge types.</param>
        /// <param name="heightUnit">Desired stage unit for stage types.</param>
        /// <param name="startDate">
        ///     Optional inclusive start date. If null, defaults vary by type:
        ///     Daily → 1800-01-01, Instantaneous → 18 months ago, Peak → year 1800.
        /// </param>
        /// <param name="endDate">
        ///     Optional inclusive end date. If null, defaults to DateTime.Today (or Now for instantaneous).
        /// </param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>TimeSeries of values.</returns>
        /// <exception cref="Exception">Thrown when the service returns an error, or an empty, header-only or HTML response (for example for an unknown station).</exception>
        public static Task<TimeSeries> FromCHMN(
            string stationNumber,
            TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge,
            DischargeUnit dischargeUnit = DischargeUnit.CubicMetersPerSecond,
            HeightUnit heightUnit = HeightUnit.Meters,
            DateTime? startDate = null,
            DateTime? endDate = null,
            CancellationToken cancellationToken = default)
        {
            return FromChmnCore(stationNumber, timeSeriesType, dischargeUnit, heightUnit, startDate, endDate, null, cancellationToken);
        }

        /// <summary>
        /// Download time series data from the Water Survey of Canada, with the data symbol of each value.
        /// </summary>
        /// <param name="stationNumber">WSC 7-character station id, e.g., "08LG010".</param>
        /// <param name="timeSeriesType">Daily, instantaneous or peak discharge or stage.</param>
        /// <param name="dischargeUnit">Desired discharge unit for discharge types.</param>
        /// <param name="heightUnit">Desired stage unit for stage types.</param>
        /// <param name="startDate">Optional inclusive start date, as for the overload without symbols.</param>
        /// <param name="endDate">Optional inclusive end date, as for the overload without symbols.</param>
        /// <param name="qualifiers">
        /// Receives one entry per ordinate of the returned series, in the same order: the Symbol/Symbole field
        /// (for example B for ice conditions, E for estimated, A for partial day), or an empty string when there is
        /// none or the day was filled as missing. The list is cleared first.
        /// </param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The same series as the overload without symbols.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="qualifiers"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for a station id that is not 7 characters or an unsupported type.</exception>
        /// <exception cref="Exception">Thrown when the service returns an error, or an empty, header-only or HTML response.</exception>
        public static Task<TimeSeries> FromCHMN(
            string stationNumber,
            TimeSeriesType timeSeriesType,
            DischargeUnit dischargeUnit,
            HeightUnit heightUnit,
            DateTime? startDate,
            DateTime? endDate,
            IList<string> qualifiers,
            CancellationToken cancellationToken = default)
        {
            if (qualifiers == null) throw new ArgumentNullException(nameof(qualifiers));
            return FromChmnCore(stationNumber, timeSeriesType, dischargeUnit, heightUnit, startDate, endDate, qualifiers, cancellationToken);
        }

        /// <summary>
        /// Validates and runs a CHMN download.
        /// </summary>
        /// <param name="stationNumber">WSC station id.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="dischargeUnit">Desired discharge unit.</param>
        /// <param name="heightUnit">Desired stage unit.</param>
        /// <param name="startDate">Optional inclusive start date.</param>
        /// <param name="endDate">Optional inclusive end date.</param>
        /// <param name="qualifiers">Receives the data symbol of each ordinate, or null.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The downloaded series.</returns>
        private static async Task<TimeSeries> FromChmnCore(
            string stationNumber,
            TimeSeriesType timeSeriesType,
            DischargeUnit dischargeUnit,
            HeightUnit heightUnit,
            DateTime? startDate,
            DateTime? endDate,
            IList<string>? qualifiers,
            CancellationToken cancellationToken)
        {
            ValidateChmnRequest(stationNumber, timeSeriesType);

            // Determine if this is a discharge or stage request
            bool isDischarge = timeSeriesType == TimeSeriesType.DailyDischarge ||
                               timeSeriesType == TimeSeriesType.InstantaneousDischarge ||
                               timeSeriesType == TimeSeriesType.PeakDischarge;

            // Build the URL based on the time series type
            string url;

            if (timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.DailyStage)
            {
                // Daily data endpoint
                string parameter = isDischarge ? "flow" : "level";
                DateTime sd = startDate ?? new DateTime(1800, 1, 1);
                DateTime ed = endDate ?? DateTime.Today;
                url = "https://wateroffice.ec.gc.ca/services/daily_data/csv/inline" +
                    $"?stations[]={Uri.EscapeDataString(stationNumber)}" +
                    $"&parameters[]={parameter}" +
                    $"&start_date={sd:yyyy-MM-dd}" +
                    $"&end_date={ed:yyyy-MM-dd}";
            }
            else if (timeSeriesType == TimeSeriesType.InstantaneousDischarge || timeSeriesType == TimeSeriesType.InstantaneousStage)
            {
                // Real-time data endpoint (parameter codes: 46=level, 47=discharge)
                string paramCode = isDischarge ? "47" : "46";
                DateTime sd = startDate ?? DateTime.UtcNow.AddMonths(-18);
                DateTime ed = endDate ?? DateTime.UtcNow;
                // Use '+' as the URL-safe space separator instead of a literal space
                // in the query string. ECCC accepts the literal space, but some
                // intermediate proxies and CDNs reject it.
                url = "https://wateroffice.ec.gc.ca/services/real_time_data/csv/inline" +
                    $"?stations[]={Uri.EscapeDataString(stationNumber)}" +
                    $"&parameters[]={paramCode}" +
                    $"&start_date={sd:yyyy-MM-dd}+{sd:HH:mm:ss}" +
                    $"&end_date={ed:yyyy-MM-dd}+{ed:HH:mm:ss}";
            }
            else // PeakDischarge or PeakStage
            {
                // Peak data endpoint
                string parameter = isDischarge ? "flow" : "level";
                int startYear = startDate?.Year ?? 1800;
                int endYear = endDate?.Year ?? DateTime.Today.Year;
                url = "https://wateroffice.ec.gc.ca/services/peak_data/csv/inline" +
                    $"?stations[]={Uri.EscapeDataString(stationNumber)}" +
                    $"&parameters[]={parameter}" +
                    $"&start_year={startYear}" +
                    $"&end_year={endYear}";
            }

            var client = _decompressClient;

            // The endpoint returns CSV (automatically decompressed by HttpClientHandler).
            TimeSpan requestTimeout = timeSeriesType == TimeSeriesType.InstantaneousDischarge ||
                                      timeSeriesType == TimeSeriesType.InstantaneousStage
                ? InstantaneousProviderRequestTimeout
                : ProviderRequestTimeout;
            var (status, csv) = await GetWithRetryAsync(
                client,
                url,
                requestTimeout: requestTimeout,
                cancellationToken: cancellationToken);
            if ((int)status >= 400)
            {
                string snippet = csv.Length > 500 ? csv.Substring(0, 500) + "…" : csv;
                throw new Exception(
                    $"Failed to retrieve CHMN data. Status={(int)status} ({status}). Body={snippet}");
            }

            if (IsNoDataResponse(csv, hasHeader: true)) throw new Exception(NoDataMessage);
            return ParseChmnText(csv, stationNumber, timeSeriesType, dischargeUnit, heightUnit, qualifiers);
        }

        /// <summary>
        /// Parses a Water Survey of Canada (wateroffice.ec.gc.ca services) CSV response into a time series.
        /// </summary>
        /// <param name="csvText">The CSV text of a daily_data, real_time_data or peak_data response.</param>
        /// <param name="stationNumber">The WSC 7-character station id; rows for other stations are skipped.</param>
        /// <param name="timeSeriesType">
        /// The time series type the text holds: daily, instantaneous (real-time) or peak discharge or stage.
        /// </param>
        /// <param name="dischargeUnit">Desired discharge unit for discharge types. Default = cubic meters per second.</param>
        /// <param name="heightUnit">Desired stage unit for stage types. Default = meters.</param>
        /// <returns>The parsed time series.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="csvText"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for a station id that is not 7 characters or an unsupported type.</exception>
        /// <remarks>
        /// <para>
        /// This is the parser <see cref="FromCHMN(string, TimeSeriesType, DischargeUnit, HeightUnit, DateTime?, DateTime?, CancellationToken)"/>
        /// applies to the downloaded text, so a saved response gives the same series offline.
        /// </para>
        /// <para>
        /// Daily types give a regular one-day series with NaN for missing days and empty values; real-time types give
        /// a five-minute series of the reported times in UTC (<see cref="DateTimeKind.Utc"/>); peak types give an irregular series of the annual maximums
        /// (minimum rows are skipped). Values are converted from the native m³/s and meters. A later row for the same
        /// date and time replaces an earlier one. The data symbols are returned only by the overload that takes a
        /// qualifier list.
        /// </para>
        /// </remarks>
        public static TimeSeries ParseCHMN(
            string csvText,
            string stationNumber,
            TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge,
            DischargeUnit dischargeUnit = DischargeUnit.CubicMetersPerSecond,
            HeightUnit heightUnit = HeightUnit.Meters)
        {
            return ParseChmnText(csvText, stationNumber, timeSeriesType, dischargeUnit, heightUnit, null);
        }

        /// <summary>
        /// Parses a Water Survey of Canada CSV response into a time series, with the data symbol of each value.
        /// </summary>
        /// <param name="csvText">The CSV text of a daily_data, real_time_data or peak_data response.</param>
        /// <param name="stationNumber">The WSC 7-character station id; rows for other stations are skipped.</param>
        /// <param name="timeSeriesType">Daily, instantaneous (real-time) or peak discharge or stage.</param>
        /// <param name="dischargeUnit">Desired discharge unit for discharge types.</param>
        /// <param name="heightUnit">Desired stage unit for stage types.</param>
        /// <param name="qualifiers">
        /// Receives one entry per ordinate of the returned series, in the same order: the Symbol/Symbole field, or an
        /// empty string when there is none or the day was filled as missing. The list is cleared first.
        /// </param>
        /// <returns>The same series as the overload without symbols.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="csvText"/> or <paramref name="qualifiers"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for a station id that is not 7 characters or an unsupported type.</exception>
        public static TimeSeries ParseCHMN(
            string csvText,
            string stationNumber,
            TimeSeriesType timeSeriesType,
            DischargeUnit dischargeUnit,
            HeightUnit heightUnit,
            IList<string> qualifiers)
        {
            if (qualifiers == null) throw new ArgumentNullException(nameof(qualifiers));
            return ParseChmnText(csvText, stationNumber, timeSeriesType, dischargeUnit, heightUnit, qualifiers);
        }

        /// <summary>
        /// Parses a CHMN CSV response, optionally collecting the data symbols.
        /// </summary>
        /// <param name="csvText">The CSV text.</param>
        /// <param name="stationNumber">The station id.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="dischargeUnit">Desired discharge unit.</param>
        /// <param name="heightUnit">Desired stage unit.</param>
        /// <param name="qualifiers">Receives the data symbol of each ordinate, or null.</param>
        /// <returns>The parsed series.</returns>
        private static TimeSeries ParseChmnText(
            string csvText,
            string stationNumber,
            TimeSeriesType timeSeriesType,
            DischargeUnit dischargeUnit,
            HeightUnit heightUnit,
            IList<string>? qualifiers)
        {
            if (csvText == null) throw new ArgumentNullException(nameof(csvText));
            ValidateChmnRequest(stationNumber, timeSeriesType);
            var qualifierByDate = qualifiers != null ? new Dictionary<DateTime, string>() : null;

            bool isDischarge = timeSeriesType == TimeSeriesType.DailyDischarge ||
                               timeSeriesType == TimeSeriesType.InstantaneousDischarge ||
                               timeSeriesType == TimeSeriesType.PeakDischarge;
            bool isPeak = timeSeriesType == TimeSeriesType.PeakDischarge || timeSeriesType == TimeSeriesType.PeakStage;
            TimeInterval interval = isPeak ? TimeInterval.Irregular
                : timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.DailyStage ? TimeInterval.OneDay
                : TimeInterval.FiveMinute;

            var ts = new TimeSeries(interval);
            var indexesByDateTime = new Dictionary<DateTime, int>();
            bool previousSuppressCollectionChanged = ts.SuppressCollectionChanged;
            ts.SuppressCollectionChanged = true;
            try
            {
                if (isPeak)
                {
                    // Peak data CSV format:
                    // ID,Parameter/Paramètre,Date,Timezone/Fuseau horaire,Type/Catégorie,Value/Valeur,Symbol/Symbole
                    ParseCHMNPeakCsv(csvText, stationNumber, isDischarge, dischargeUnit, heightUnit, ts, indexesByDateTime, qualifierByDate);
                }
                else
                {
                    // Daily and real-time CSV share the same column layout:
                    // ID,Date,Parameter/Paramètre,Value/Valeur,...
                    ParseCHMNDailyCsv(csvText, stationNumber, timeSeriesType, isDischarge, dischargeUnit, heightUnit, ts, indexesByDateTime, qualifierByDate);
                }
            }
            finally
            {
                ts.SuppressCollectionChanged = previousSuppressCollectionChanged;
            }

            if (qualifiers != null) WriteQualifiers(ts, qualifierByDate!, qualifiers);
            return ts;
        }

        /// <summary>
        /// Checks a CHMN station id and time series type.
        /// </summary>
        /// <param name="stationNumber">The WSC station id.</param>
        /// <param name="timeSeriesType">The requested type.</param>
        /// <exception cref="ArgumentException">Thrown for an id that is not 7 characters or an unsupported type.</exception>
        private static void ValidateChmnRequest(string stationNumber, TimeSeriesType timeSeriesType)
        {
            // Basic validation (WSC station numbers are 7 chars including letters)
            if (string.IsNullOrWhiteSpace(stationNumber) || stationNumber.Length != 7)
                throw new ArgumentException("The WSC station number must be 7 characters, e.g., 08LG010.", nameof(stationNumber));

            // Supported series types
            if (timeSeriesType == TimeSeriesType.DailyPrecipitation || timeSeriesType == TimeSeriesType.DailySnow ||
                timeSeriesType == TimeSeriesType.MeasuredDischarge || timeSeriesType == TimeSeriesType.MeasuredStage)
                throw new ArgumentException("Canadian API supports Daily, Instantaneous, and Peak discharge/stage types. Field measurements are not available through the API.", nameof(timeSeriesType));
        }

        /// <summary>
        /// Parse CHMN daily or real-time CSV data.
        /// </summary>
        /// <param name="csv">The CSV response body returned by the CHMN endpoint.</param>
        /// <param name="stationNumber">The CHMN station identifier to retain while parsing.</param>
        /// <param name="timeSeriesType">The requested time series type.</param>
        /// <param name="isDischarge">A value indicating whether discharge values are being parsed.</param>
        /// <param name="dischargeUnit">The discharge unit requested by the caller.</param>
        /// <param name="heightUnit">The height unit requested by the caller.</param>
        /// <param name="ts">The time series being populated.</param>
        /// <param name="indexesByDateTime">Lookup used to replace duplicate dates inline while parsing.</param>
        /// <param name="qualifierByDate">Receives the Symbol/Symbole field of each row by date, or null.</param>
        /// <remarks>
        /// CSV format: ID,Date,Parameter/Paramètre,Value/Valeur,... (additional columns vary by endpoint)
        /// </remarks>
        private static void ParseCHMNDailyCsv(string csv, string stationNumber, TimeSeriesType timeSeriesType,
            bool isDischarge, DischargeUnit dischargeUnit, HeightUnit heightUnit, TimeSeries ts,
            Dictionary<DateTime, int> indexesByDateTime, Dictionary<DateTime, string>? qualifierByDate)
        {
            bool isDailyType = timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.DailyStage;

            using (var sr = new StringReader(csv))
            {
                string? line;
                int idCol = 0, dateCol = 1, paramCol = 2, valueCol = 3, symbolCol = -1;
                bool headerParsed = false;
                DateTime? prevDate = null;

                while ((line = sr.ReadLine()) != null)
                {
                    if (string.IsNullOrWhiteSpace(line)) continue;

                    // Detect header line
                    if (!headerParsed)
                    {
                        if (line.Contains("ID") && line.Contains("Date") && line.Contains("Value"))
                        {
                            headerParsed = true;
                            symbolCol = ChmnSymbolColumn(line);
                        }
                        continue;
                    }

                    // Parse data rows
                    string[] parts = line.Split(',');
                    if (parts.Length <= Math.Max(Math.Max(idCol, dateCol), valueCol))
                        continue;

                    // Check station match
                    if (!parts[idCol].Trim().Equals(stationNumber, StringComparison.OrdinalIgnoreCase))
                        continue;

                    // Check parameter match for daily data (parameter column contains text like "discharge/débit" or "level")
                    if (isDailyType && paramCol >= 0 && parts.Length > paramCol)
                    {
                        string paramValue = parts[paramCol].Trim().ToLowerInvariant();
                        bool isFlowParam = paramValue.Contains("discharge") || paramValue.Contains("débit") || paramValue.Contains("flow");
                        bool isLevelParam = paramValue.Contains("level") || paramValue.Contains("stage") || paramValue.Contains("niveau");

                        if (isDischarge && !isFlowParam) continue;
                        if (!isDischarge && !isLevelParam) continue;
                    }

                    // Real-time timestamps are UTC with a Z suffix, e.g. "2026-11-01T08:00:00Z" (wateroffice
                    // services index). DateTime.TryParse with DateTimeStyles.None would convert them to the
                    // machine's local time, so they are read as DateTimeOffset and returned as UTC
                    // (DateTimeKind.Utc). Daily rows are calendar dates.
                    if (!DateTimeOffset.TryParse(parts[dateCol].Trim(), CultureInfo.InvariantCulture,
                            DateTimeStyles.None, out DateTimeOffset dto))
                        continue;
                    DateTime date = isDailyType ? dto.Date : dto.UtcDateTime;

                    // Parse value
                    double val = double.NaN;
                    string valueStr = parts[valueCol].Trim();

                    if (!string.IsNullOrWhiteSpace(valueStr))
                    {
                        if (double.TryParse(valueStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double raw))
                        {
                            val = ConvertCHMNValue(raw, isDischarge, dischargeUnit, heightUnit);
                        }
                    }

                    // Fill gaps with NaN for daily series only
                    if (isDailyType && prevDate.HasValue && (date - prevDate.Value).Days > 1)
                        FillMissingDates(ts, prevDate.Value, date, indexesByDateTime);

                    AddOrReplaceByDateTime(ts, indexesByDateTime, date, val);
                    if (qualifierByDate != null) qualifierByDate[date] = ChmnSymbol(parts, symbolCol);
                    prevDate = date;
                }
            }
        }

        /// <summary>
        /// Finds the Symbol/Symbole column of a CHMN header line.
        /// </summary>
        /// <param name="header">The header line.</param>
        /// <returns>The column position, or -1.</returns>
        private static int ChmnSymbolColumn(string header)
        {
            string[] names = header.Split(',');
            for (int i = 0; i < names.Length; i++)
                if (names[i].IndexOf("Symbol", StringComparison.OrdinalIgnoreCase) >= 0) return i;
            return -1;
        }

        /// <summary>
        /// Returns the data symbol of a CHMN row, or an empty string.
        /// </summary>
        /// <param name="parts">The row fields.</param>
        /// <param name="symbolCol">The Symbol/Symbole column, or -1.</param>
        private static string ChmnSymbol(string[] parts, int symbolCol)
        {
            return symbolCol >= 0 && symbolCol < parts.Length ? parts[symbolCol].Trim() : "";
        }

        /// <summary>
        /// Parse CHMN peak data CSV.
        /// </summary>
        /// <param name="csv">The CSV response body returned by the CHMN endpoint.</param>
        /// <param name="stationNumber">The CHMN station identifier to retain while parsing.</param>
        /// <param name="isDischarge">A value indicating whether discharge values are being parsed.</param>
        /// <param name="dischargeUnit">The discharge unit requested by the caller.</param>
        /// <param name="heightUnit">The height unit requested by the caller.</param>
        /// <param name="ts">The time series being populated.</param>
        /// <param name="indexesByDateTime">Lookup used to replace duplicate dates inline while parsing.</param>
        /// <param name="qualifierByDate">Receives the Symbol/Symbole field of each row by date, or null.</param>
        /// <remarks>
        /// CSV format: ID,Parameter/Paramètre,Date,Timezone/Fuseau horaire,Type/Catégorie,Value/Valeur,Symbol/Symbole.
        /// Only "maximum" rows are included (minimums are skipped).
        /// </remarks>
        private static void ParseCHMNPeakCsv(string csv, string stationNumber, bool isDischarge,
            DischargeUnit dischargeUnit, HeightUnit heightUnit, TimeSeries ts,
            Dictionary<DateTime, int> indexesByDateTime, Dictionary<DateTime, string>? qualifierByDate)
        {
            using (var sr = new StringReader(csv))
            {
                string? line;
                bool headerParsed = false;

                // Peak CSV columns: ID(0), Parameter(1), Date(2), Timezone(3), Type(4), Value(5), Symbol(6)
                int idCol = 0, dateCol = 2, typeCol = 4, valueCol = 5, symbolCol = -1;

                while ((line = sr.ReadLine()) != null)
                {
                    if (string.IsNullOrWhiteSpace(line)) continue;

                    // Detect header line
                    if (!headerParsed)
                    {
                        if (line.Contains("ID") && line.Contains("Date") && line.Contains("Value"))
                        {
                            headerParsed = true;
                            symbolCol = ChmnSymbolColumn(line);
                        }
                        continue;
                    }

                    string[] parts = line.Split(',');
                    if (parts.Length <= valueCol) continue;

                    // Check station match
                    if (!parts[idCol].Trim().Equals(stationNumber, StringComparison.OrdinalIgnoreCase))
                        continue;

                    // Only include maximum values (skip minimums)
                    string typeValue = parts[typeCol].Trim().ToLowerInvariant();
                    if (typeValue != "maximum") continue;

                    // Parse date
                    if (!DateTime.TryParse(parts[dateCol].Trim(), CultureInfo.InvariantCulture, DateTimeStyles.None, out DateTime date))
                        continue;

                    // Parse value
                    string valueStr = parts[valueCol].Trim();
                    if (string.IsNullOrWhiteSpace(valueStr)) continue;
                    if (!double.TryParse(valueStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double raw))
                        continue;

                    double val = ConvertCHMNValue(raw, isDischarge, dischargeUnit, heightUnit);
                    AddOrReplaceByDateTime(ts, indexesByDateTime, date, val);
                    if (qualifierByDate != null) qualifierByDate[date] = ChmnSymbol(parts, symbolCol);
                }
            }
        }

        /// <summary>
        /// Convert a CHMN raw value to the desired unit.
        /// </summary>
        /// <remarks>
        /// WSC native units: discharge in m³/s, stage in meters.
        /// </remarks>
        private static double ConvertCHMNValue(double raw, bool isDischarge, DischargeUnit dischargeUnit, HeightUnit heightUnit)
        {
            if (isDischarge)
            {
                return dischargeUnit == DischargeUnit.CubicFeetPerSecond
                    ? raw * 35.3146667   // cms → cfs
                    : raw;               // cms
            }
            else
            {
                return heightUnit == HeightUnit.Feet
                    ? raw * 3.280839895  // m → ft
                    : raw;               // m
            }
        }

        #endregion

        #region Australian Bureau of Meteorology (BOM)

        /// <summary>
        /// Download time series data from the Australian Bureau of Meteorology (BOM)
        /// via the KiWIS API.
        /// </summary>
        /// <param name="stationNumber">BOM station number (e.g., "410730").</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="dischargeUnit">Desired discharge unit for discharge types.</param>
        /// <param name="heightUnit">Desired stage unit for stage types.</param>
        /// <param name="depthUnit">Desired depth unit for precipitation types.</param>
        /// <param name="startDate">Optional start date. If null, attempts to retrieve full period of record.</param>
        /// <param name="endDate">Optional end date. If null, defaults to today.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>TimeSeries of values.</returns>
        public static Task<TimeSeries> FromABOM(
            string stationNumber,
            TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge,
            DischargeUnit dischargeUnit = DischargeUnit.CubicMetersPerSecond,
            HeightUnit heightUnit = HeightUnit.Meters,
            DepthUnit depthUnit = DepthUnit.Millimeters,
            DateTime? startDate = null,
            DateTime? endDate = null,
            CancellationToken cancellationToken = default)
        {
            return FromAbomCore(stationNumber, timeSeriesType, dischargeUnit, heightUnit, depthUnit, startDate, endDate, null, cancellationToken);
        }

        /// <summary>
        /// Download time series data from the Australian Bureau of Meteorology (BOM), with the quality code of each value.
        /// </summary>
        /// <param name="stationNumber">BOM station number (e.g., "410730").</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="dischargeUnit">Desired discharge unit for discharge types.</param>
        /// <param name="heightUnit">Desired stage unit for stage types.</param>
        /// <param name="depthUnit">Desired depth unit for precipitation types.</param>
        /// <param name="startDate">Optional start date. If null, attempts to retrieve full period of record.</param>
        /// <param name="endDate">Optional end date. If null, defaults to today.</param>
        /// <param name="qualifiers">
        /// Receives one entry per ordinate of the returned series, in the same order: the KiWIS quality code (for
        /// example 10, 90, 110, 130 or 140), or an empty string when there is none or the day was filled as missing.
        /// The list is cleared first.
        /// </param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The same series as the overload without quality codes.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="qualifiers"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for an invalid station number or an unsupported type.</exception>
        /// <exception cref="Exception">Thrown when the service returns an error or no matching series.</exception>
        /// <remarks>The values request also asks for the Quality Code field (returnfields=Timestamp,Value,Quality Code).</remarks>
        public static Task<TimeSeries> FromABOM(
            string stationNumber,
            TimeSeriesType timeSeriesType,
            DischargeUnit dischargeUnit,
            HeightUnit heightUnit,
            DepthUnit depthUnit,
            DateTime? startDate,
            DateTime? endDate,
            IList<string> qualifiers,
            CancellationToken cancellationToken = default)
        {
            if (qualifiers == null) throw new ArgumentNullException(nameof(qualifiers));
            return FromAbomCore(stationNumber, timeSeriesType, dischargeUnit, heightUnit, depthUnit, startDate, endDate, qualifiers, cancellationToken);
        }

        /// <summary>
        /// Validates and runs a BOM download.
        /// </summary>
        /// <param name="stationNumber">BOM station number.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="dischargeUnit">Desired discharge unit.</param>
        /// <param name="heightUnit">Desired stage unit.</param>
        /// <param name="depthUnit">Desired depth unit.</param>
        /// <param name="startDate">Optional start date.</param>
        /// <param name="endDate">Optional end date.</param>
        /// <param name="qualifiers">Receives the quality code of each ordinate, or null.</param>
        /// <param name="cancellationToken">Token used to cancel the download.</param>
        /// <returns>The downloaded series.</returns>
        private static async Task<TimeSeries> FromAbomCore(
            string stationNumber,
            TimeSeriesType timeSeriesType,
            DischargeUnit dischargeUnit,
            HeightUnit heightUnit,
            DepthUnit depthUnit,
            DateTime? startDate,
            DateTime? endDate,
            IList<string>? qualifiers,
            CancellationToken cancellationToken)
        {
            // Validate station number (BOM station numbers are typically 6 digits)
            if (string.IsNullOrWhiteSpace(stationNumber) || stationNumber.Length < 6)
                throw new ArgumentException("BOM station number must be at least 6 digits.", nameof(stationNumber));

            ValidateAbomType(timeSeriesType);

            // Set default dates
            DateTime sd = startDate ?? new DateTime(1800, 1, 1);
            DateTime ed = endDate ?? DateTime.Today;

            string parameterType = AbomParameterType(timeSeriesType);
            bool isInstantaneous = timeSeriesType == TimeSeriesType.InstantaneousDischarge || timeSeriesType == TimeSeriesType.InstantaneousStage;

            // Step 1: Get timeseries list to find the appropriate ts_id
            string tsListUrl = "https://www.bom.gov.au/waterdata/services" +
                $"?service=kisters&type=QueryServices&datasource=0&format=json" +
                $"&request=getTimeseriesList" +
                $"&station_no={Uri.EscapeDataString(stationNumber)}" +
                $"&parametertype_name={Uri.EscapeDataString(parameterType)}";

            // The shared BOM client carries the browser-like headers the BOM security
            // appliance expects and reuses pooled connections across both KiWIS requests.
            var client = _bomClient;

            HttpStatusCode listStatus;
            string listResponse;
            try
            {
                (listStatus, listResponse) = await GetWithRetryAsync(client, tsListUrl, cancellationToken: cancellationToken);
            }
            catch (HttpRequestException ex)
            {
                throw new Exception($"Failed to connect to BOM API: {ex.Message}", ex);
            }
            catch (TaskCanceledException ex)
            {
                throw new Exception("Request to BOM API timed out.", ex);
            }

            // Surface the response body in the exception so transient upstream failures
            // (KiWIS 5xx with a JSON error payload) are self-documenting.
            if ((int)listStatus >= 400)
            {
                string snippet = listResponse.Length > 500 ? listResponse.Substring(0, 500) + "…" : listResponse;
                throw new Exception(
                    $"Failed to retrieve time series list from BOM API. Status={(int)listStatus} ({listStatus}). Body={snippet}");
            }

            string tsId = SelectABOMTimeSeriesId(listResponse, stationNumber, timeSeriesType);

            // Step 2: Get time series values using the ts_id
            string valuesUrl = "https://www.bom.gov.au/waterdata/services" +
                $"?service=kisters&type=QueryServices&datasource=0&format=json" +
                $"&request=getTimeseriesValues" +
                $"&ts_id={tsId}" +
                $"&from={sd:yyyy-MM-dd}" +
                $"&to={ed:yyyy-MM-dd}" +
                (qualifiers != null ? "&returnfields=Timestamp%2CValue%2CQuality+Code" : "");

            // Reuse the shared BOM client so the values request rides the connection the
            // station lookup already established.
            HttpStatusCode valuesStatus;
            string valuesResponse;
            try
            {
                (valuesStatus, valuesResponse) = await GetWithRetryAsync(
                    client,
                    valuesUrl,
                    requestTimeout: isInstantaneous ? InstantaneousProviderRequestTimeout : ProviderRequestTimeout,
                    cancellationToken: cancellationToken);
            }
            catch (HttpRequestException ex)
            {
                throw new Exception($"Failed to retrieve time series values from BOM API: {ex.Message}", ex);
            }
            catch (TaskCanceledException ex)
            {
                throw new Exception("Request to BOM API timed out while retrieving values.", ex);
            }

            // Surface the response body in the exception so transient upstream failures
            // (e.g. KiWIS returning 500 with `{"code":"DatasourceError","message":"Error connecting to WDP."}`)
            // are self-documenting after the retry helper has exhausted its attempts.
            if ((int)valuesStatus >= 400)
            {
                string snippet = valuesResponse.Length > 500 ? valuesResponse.Substring(0, 500) + "…" : valuesResponse;
                throw new Exception(
                    $"Failed to retrieve time series values from BOM API. Status={(int)valuesStatus} ({valuesStatus}). Body={snippet}");
            }

            return ParseAbomText(valuesResponse, timeSeriesType, dischargeUnit, heightUnit, depthUnit, qualifiers);
        }

        /// <summary>
        /// Checks that a BOM request asks for a supported time series type.
        /// </summary>
        /// <param name="timeSeriesType">The requested type.</param>
        /// <exception cref="ArgumentException">Thrown for peak, field-measurement and snow types.</exception>
        private static void ValidateAbomType(TimeSeriesType timeSeriesType)
        {
            // BOM supports daily/instantaneous discharge and stage, plus daily precipitation
            var supportedTypes = new[] {
                TimeSeriesType.DailyDischarge, TimeSeriesType.DailyStage,
                TimeSeriesType.InstantaneousDischarge, TimeSeriesType.InstantaneousStage,
                TimeSeriesType.DailyPrecipitation
            };
            if (!supportedTypes.Contains(timeSeriesType))
                throw new ArgumentException(
                    "BOM API supports DailyDischarge, DailyStage, InstantaneousDischarge, InstantaneousStage, and DailyPrecipitation.",
                    nameof(timeSeriesType));
        }

        /// <summary>
        /// Returns the KiWIS parameter type name of a time series type.
        /// </summary>
        /// <param name="timeSeriesType">A supported BOM type.</param>
        /// <returns>The parametertype_name used in the station lookup.</returns>
        private static string AbomParameterType(TimeSeriesType timeSeriesType)
        {
            if (timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.InstantaneousDischarge) return "Water Course Discharge";
            if (timeSeriesType == TimeSeriesType.DailyStage || timeSeriesType == TimeSeriesType.InstantaneousStage) return "Water Course Level";
            return "Rainfall";
        }

        /// <summary>
        /// Chooses the KiWIS time series id for a requested type from a getTimeseriesList response.
        /// </summary>
        /// <param name="listJson">The getTimeseriesList JSON: a header row, then one row per series.</param>
        /// <param name="stationNumber">The station number, for error messages.</param>
        /// <param name="timeSeriesType">The requested type.</param>
        /// <returns>The ts_id of the chosen series.</returns>
        /// <exception cref="Exception">Thrown for an error response, malformed JSON, or no matching series.</exception>
        /// <remarks>
        /// Quality-controlled merged series (DMQaQc.Merged) are preferred over other matches. A daily type takes the
        /// midnight-to-midnight series (ts_name ending .24HR) whatever the service order, so a 9am-to-9am series
        /// (.09HR) is never picked when a .24HR series exists; a station without a .24HR series keeps the first match.
        /// </remarks>
        internal static string SelectABOMTimeSeriesId(string listJson, string stationNumber, TimeSeriesType timeSeriesType)
        {
            string parameterType = AbomParameterType(timeSeriesType);
            bool isInstantaneous = timeSeriesType == TimeSeriesType.InstantaneousDischarge || timeSeriesType == TimeSeriesType.InstantaneousStage;
            bool isPrecip = timeSeriesType == TimeSeriesType.DailyPrecipitation;

            // Check for error response from KiWIS
            if (listJson.Contains("\"type\":\"error\"") || listJson.StartsWith("{\"type\":\"error\""))
            {
                throw new Exception($"BOM API returned an error. Response: {listJson}");
            }

            // Parse JSON response to find appropriate timeseries
            // BOM returns array where first element is column headers, subsequent elements are data rows
            using var listData = ParseAbomJson(listJson, "BOM API response");
            var root = listData.RootElement;

            if (root.ValueKind != System.Text.Json.JsonValueKind.Array || root.GetArrayLength() < 2)
                throw new Exception($"No time series found for station {stationNumber} with parameter {parameterType}. " +
                    $"API returned {root.GetArrayLength()} items.");

            // Find ts_id column index from header row
            int tsIdIndex = -1;
            int tsNameIndex = -1;
            var headers = root[0];

            if (headers.ValueKind != System.Text.Json.JsonValueKind.Array)
                throw new Exception("Unexpected response format: first element is not an array of headers");

            for (int i = 0; i < headers.GetArrayLength(); i++)
            {
                string? header = headers[i].GetString();
                if (header == "ts_id") tsIdIndex = i;
                if (header == "ts_name") tsNameIndex = i;
            }

            if (tsIdIndex == -1)
                throw new Exception("Could not find ts_id in response");

            // Determine the ts_name pattern to search for.
            // Prioritize DMQaQc.Merged (quality-controlled merged data).
            string tsNamePattern;
            if (isInstantaneous)
                tsNamePattern = "AsStored";
            else if (isPrecip)
                tsNamePattern = "DailyTotal";
            else
                tsNamePattern = "DailyMean";

            // Search passes: a daily type first looks for the .24HR series (merged, then any), then every type looks
            // for a DMQaQc.Merged match, then any match.
            string? tsId = null;
            for (int pass = isInstantaneous ? 2 : 0; pass < 4 && tsId == null; pass++)
            {
                for (int i = 1; i < root.GetArrayLength(); i++)
                {
                    var row = root[i];
                    string tsName = tsNameIndex >= 0 ? row[tsNameIndex].GetString() ?? "" : "";

                    if (!tsName.Contains(tsNamePattern)) continue;
                    if (pass < 2 && !tsName.EndsWith(".24HR", StringComparison.Ordinal)) continue;
                    if ((pass == 0 || pass == 2) && !tsName.StartsWith("DMQaQc.Merged")) continue;

                    tsId = row[tsIdIndex].GetString();
                    break;
                }
            }

            // No silent fallback: if the requested cadence is not present, throw rather than
            // quietly hand back hourly data when the caller asked for daily.
            if (tsId == null)
            {
                var available = new List<string>();
                for (int i = 1; i < root.GetArrayLength() && available.Count < 10; i++)
                    if (tsNameIndex >= 0)
                        available.Add(root[i][tsNameIndex].GetString() ?? "");
                throw new Exception(
                    $"No suitable time series found for station {stationNumber} matching '{tsNamePattern}'. " +
                    $"Available ts_name values: [{string.Join(", ", available)}]");
            }
            return tsId;
        }

        /// <summary>
        /// Parses a BOM Water Data Online (KiWIS) getTimeseriesValues JSON response into a time series.
        /// </summary>
        /// <param name="valuesJson">The JSON response: an array whose first object holds a "data" array of [timestamp, value] rows.</param>
        /// <param name="timeSeriesType">The time series type the response holds. Default = Daily discharge.</param>
        /// <param name="dischargeUnit">Desired discharge unit for discharge types. Default = cubic meters per second.</param>
        /// <param name="heightUnit">Desired stage unit for stage types. Default = meters.</param>
        /// <param name="depthUnit">Desired depth unit for precipitation types. Default = millimeters.</param>
        /// <returns>The parsed time series.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="valuesJson"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for an unsupported type.</exception>
        /// <exception cref="Exception">Thrown for a KiWIS error response, or a response with no data array.</exception>
        /// <remarks>
        /// <para>
        /// This is the parser <see cref="FromABOM(string, TimeSeriesType, DischargeUnit, HeightUnit, DepthUnit, DateTime?, DateTime?, CancellationToken)"/>
        /// applies to the downloaded values, so a saved response gives the same series offline.
        /// </para>
        /// <para>
        /// Timestamps carry an explicit offset (for example +10:00). Daily types are dated by the reported wall-clock date
        /// and form a regular one-day series with NaN for missing days and non-numeric values; instantaneous types keep
        /// the reported wall-clock time in an irregular series. Values are converted from the native m³/s, meters and
        /// millimeters. A later row for the same date and time replaces an earlier one. Quality codes are not read.
        /// </para>
        /// </remarks>
        public static TimeSeries ParseABOM(
            string valuesJson,
            TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge,
            DischargeUnit dischargeUnit = DischargeUnit.CubicMetersPerSecond,
            HeightUnit heightUnit = HeightUnit.Meters,
            DepthUnit depthUnit = DepthUnit.Millimeters)
        {
            return ParseAbomText(valuesJson, timeSeriesType, dischargeUnit, heightUnit, depthUnit, null);
        }

        /// <summary>
        /// Parses a BOM KiWIS getTimeseriesValues JSON response into a time series, with the quality code of each value.
        /// </summary>
        /// <param name="valuesJson">The JSON response; its "columns" field names a Quality Code column when the request asked for one.</param>
        /// <param name="timeSeriesType">The time series type the response holds.</param>
        /// <param name="dischargeUnit">Desired discharge unit for discharge types.</param>
        /// <param name="heightUnit">Desired stage unit for stage types.</param>
        /// <param name="depthUnit">Desired depth unit for precipitation types.</param>
        /// <param name="qualifiers">
        /// Receives one entry per ordinate of the returned series, in the same order: the quality code, or an empty
        /// string when the response has none or the day was filled as missing. The list is cleared first.
        /// </param>
        /// <returns>The same series as the overload without quality codes.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="valuesJson"/> or <paramref name="qualifiers"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown for an unsupported type.</exception>
        /// <exception cref="Exception">Thrown for a KiWIS error response, or a response with no data array.</exception>
        public static TimeSeries ParseABOM(
            string valuesJson,
            TimeSeriesType timeSeriesType,
            DischargeUnit dischargeUnit,
            HeightUnit heightUnit,
            DepthUnit depthUnit,
            IList<string> qualifiers)
        {
            if (qualifiers == null) throw new ArgumentNullException(nameof(qualifiers));
            return ParseAbomText(valuesJson, timeSeriesType, dischargeUnit, heightUnit, depthUnit, qualifiers);
        }

        /// <summary>
        /// Parses a BOM values response, optionally collecting the quality codes.
        /// </summary>
        /// <param name="valuesJson">The JSON response.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <param name="dischargeUnit">Desired discharge unit.</param>
        /// <param name="heightUnit">Desired stage unit.</param>
        /// <param name="depthUnit">Desired depth unit.</param>
        /// <param name="qualifiers">Receives the quality code of each ordinate, or null.</param>
        /// <returns>The parsed series.</returns>
        private static TimeSeries ParseAbomText(
            string valuesJson,
            TimeSeriesType timeSeriesType,
            DischargeUnit dischargeUnit,
            HeightUnit heightUnit,
            DepthUnit depthUnit,
            IList<string>? qualifiers)
        {
            if (valuesJson == null) throw new ArgumentNullException(nameof(valuesJson));
            ValidateAbomType(timeSeriesType);
            var qualifierByDate = qualifiers != null ? new Dictionary<DateTime, string>() : null;

            bool isDischarge = timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.InstantaneousDischarge;
            bool isStage = timeSeriesType == TimeSeriesType.DailyStage || timeSeriesType == TimeSeriesType.InstantaneousStage;
            bool isPrecip = timeSeriesType == TimeSeriesType.DailyPrecipitation;
            bool isInstantaneous = timeSeriesType == TimeSeriesType.InstantaneousDischarge || timeSeriesType == TimeSeriesType.InstantaneousStage;

            var ts = isInstantaneous ? new TimeSeries(TimeInterval.Irregular) : new TimeSeries(TimeInterval.OneDay);
            DateTime? prevDate = null;
            var indexesByDateTime = new Dictionary<DateTime, int>();

            // Check for error response
            if (valuesJson.Contains("\"type\":\"error\"") || valuesJson.StartsWith("{\"type\":\"error\""))
            {
                throw new Exception($"BOM API returned an error for time series values. Response: {valuesJson}");
            }

            using var valuesData = ParseAbomJson(valuesJson, "BOM time series values");
            var root = valuesData.RootElement;

            // Response structure: array with single object containing 'data' array
            if (root.ValueKind != System.Text.Json.JsonValueKind.Array || root.GetArrayLength() == 0)
                throw new Exception("No data returned from BOM");

            var dataObj = root[0];
            if (dataObj.ValueKind != System.Text.Json.JsonValueKind.Object || !dataObj.TryGetProperty("data", out var dataArray) ||
                dataArray.ValueKind != System.Text.Json.JsonValueKind.Array)
                throw new Exception("Invalid response format from BOM");
            int qualityCol = qualifierByDate != null ? AbomQualityColumn(dataObj) : -1;

            bool previousSuppressCollectionChanged = ts.SuppressCollectionChanged;
            ts.SuppressCollectionChanged = true;
            try
            {
                // Parse the data array - each element is [timestamp, value, ...]
                foreach (var point in dataArray.EnumerateArray())
                {
                    if (point.GetArrayLength() < 2) continue;

                    // BOM returns ISO-8601 with explicit offset, e.g. "2021-11-01T00:00:00.000+10:00".
                    // DateTime.TryParse converts these to LOCAL time, which on a non-AEST machine
                    // shifts daily series labels by a calendar day. DateTimeOffset preserves the
                    // offset; .Date gives the wall-clock date BOM intended.
                    string? timestampStr = point[0].GetString();
                    if (!DateTimeOffset.TryParse(timestampStr, CultureInfo.InvariantCulture,
                            DateTimeStyles.None, out DateTimeOffset dto))
                        continue;
                    DateTime date = isInstantaneous ? dto.DateTime : dto.Date;

                    // Parse value
                    double val = double.NaN;
                    var valueElement = point[1];

                    if (valueElement.ValueKind == System.Text.Json.JsonValueKind.Number)
                    {
                        double raw = valueElement.GetDouble();

                        // Apply unit conversion
                        if (isDischarge)
                        {
                            // BOM returns discharge in m³/s
                            val = dischargeUnit == DischargeUnit.CubicFeetPerSecond
                                ? raw * 35.3146667
                                : raw;
                        }
                        else if (isStage)
                        {
                            // BOM returns stage in meters
                            val = heightUnit == HeightUnit.Feet
                                ? raw * 3.280839895
                                : raw;
                        }
                        else if (isPrecip)
                        {
                            // BOM returns rainfall in millimeters
                            val = depthUnit switch
                            {
                                DepthUnit.Millimeters => raw,
                                DepthUnit.Centimeters => raw / 10.0,
                                DepthUnit.Inches => raw / 25.4,
                                _ => raw
                            };
                        }
                    }

                    // Fill gaps with NaN for daily series only
                    if (!isInstantaneous && prevDate.HasValue && (date - prevDate.Value).Days > 1)
                        FillMissingDates(ts, prevDate.Value, date, indexesByDateTime);

                    AddOrReplaceByDateTime(ts, indexesByDateTime, date, val);
                    if (qualifierByDate != null)
                        qualifierByDate[date] = qualityCol >= 0 && point.GetArrayLength() > qualityCol && point[qualityCol].ValueKind != System.Text.Json.JsonValueKind.Null
                            ? point[qualityCol].ToString()
                            : "";
                    prevDate = date;
                }
            }
            finally
            {
                ts.SuppressCollectionChanged = previousSuppressCollectionChanged;
            }

            if (qualifiers != null) WriteQualifiers(ts, qualifierByDate!, qualifiers);
            return ts;
        }

        /// <summary>
        /// Parses a BOM KiWIS JSON response, reporting malformed JSON with the start of the response.
        /// </summary>
        /// <param name="json">The response text.</param>
        /// <param name="what">What the response holds, for the message.</param>
        /// <returns>The parsed document; the caller disposes it.</returns>
        /// <exception cref="Exception">Thrown when the text is not valid JSON.</exception>
        private static System.Text.Json.JsonDocument ParseAbomJson(string json, string what)
        {
            try
            {
                return System.Text.Json.JsonDocument.Parse(json);
            }
            catch (System.Text.Json.JsonException ex)
            {
                throw new Exception($"Failed to parse {what} as JSON. Response: {json.Substring(0, Math.Min(500, json.Length))}", ex);
            }
        }

        /// <summary>
        /// Finds the Quality Code position in a KiWIS values object from its "columns" field.
        /// </summary>
        /// <param name="dataObj">The values object.</param>
        /// <returns>The position in each data row, or -1.</returns>
        private static int AbomQualityColumn(System.Text.Json.JsonElement dataObj)
        {
            if (!dataObj.TryGetProperty("columns", out var columns) || columns.ValueKind != System.Text.Json.JsonValueKind.String) return -1;
            string[] names = (columns.GetString() ?? "").Split(',');
            for (int i = 0; i < names.Length; i++)
                if (string.Equals(names[i].Trim(), "Quality Code", StringComparison.OrdinalIgnoreCase)) return i;
            return -1;
        }

        #endregion

    }
}
