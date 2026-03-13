/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using System.IO.Compression;
using System.IO;
using System.Linq;
using System.Net;
using System.Xml.Linq;
using System;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using System.Net.Http;
using System.Globalization;
using System.Collections.Generic;

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
        private static readonly HttpClient _defaultClient = new HttpClient() { Timeout = TimeSpan.FromSeconds(60) };
        private static readonly HttpClient _decompressClient = new HttpClient(new HttpClientHandler
        {
            AutomaticDecompression = DecompressionMethods.GZip | DecompressionMethods.Deflate
        })
        { Timeout = TimeSpan.FromSeconds(60) };

        private const string UserAgent = "USACE-Numerics/2.0";

        /// <summary>
        /// Checks if there is an Internet connection.
        /// </summary>
        public static async Task<bool> IsConnectedToInternet()
        {
            try
            {
                using (var cts = new CancellationTokenSource(TimeSpan.FromSeconds(5)))
                {
                    await _defaultClient.GetAsync("https://waterservices.usgs.gov/nwis/", cts.Token);
                    return true;
                }
            }
            catch
            {
                return false;
            }
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
        /// <returns>A downloaded time series.</returns>
        public static async Task<TimeSeries> FromGHCN(string siteNumber, TimeSeriesType timeSeriesType = TimeSeriesType.DailyPrecipitation, DepthUnit unit = DepthUnit.Inches)
        {
            

            // Check site number
            if (siteNumber.Length != 11 || !Regex.IsMatch(siteNumber, @"^[A-Za-z0-9]{11}$"))
            {
                throw new ArgumentException("The GHCN site number must be exactly 11 alphanumeric characters.", nameof(siteNumber));
            }

            // Check time series type
            if (timeSeriesType != TimeSeriesType.DailyPrecipitation && timeSeriesType != TimeSeriesType.DailySnow)
            {
                throw new ArgumentException("The time series type must be either daily precipitation or daily snow.", nameof(timeSeriesType));
            }

            var timeSeries = new TimeSeries(TimeInterval.OneDay);
            DateTime? previousDate = null;
            string tempFilePath = Path.Combine(Path.GetTempPath(), $"{Path.GetRandomFileName()}.dly");

            // Check internet connection
            if (!await IsConnectedToInternet())
            {
                throw new InvalidOperationException("No internet connection.");
            }

            try
            {
 
                // Download the GHCN file
                string ghcnBaseUrl = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/";
                string stationFileUrl = $"{ghcnBaseUrl}{Uri.EscapeDataString(siteNumber)}.dly";
                {
                    var client = _defaultClient;
                    // Request the file and ensure that response headers are read first.
                    using (HttpResponseMessage response = await client.GetAsync(stationFileUrl, HttpCompletionOption.ResponseHeadersRead))
                    {
                        response.EnsureSuccessStatusCode();
                        using (Stream stream = await response.Content.ReadAsStreamAsync())
                        using (FileStream fs = new FileStream(tempFilePath, FileMode.Create, FileAccess.Write, FileShare.None, bufferSize: 81920, useAsync: true))
                        {
                            await stream.CopyToAsync(fs);
                        }
                    }
                }

                // Read and parse the file
                if (!File.Exists(tempFilePath))
                {
                    throw new Exception("Downloaded file does not exist.");
                }

                // Read and parse the file
                string[] lines = File.ReadAllLines(tempFilePath);
                foreach (string line in lines)
                {
                    // Extract the type (PRCP, SNOW, etc.)
                    var typeString = line.Substring(17, 4);
                    if ((typeString == "PRCP" && timeSeriesType == TimeSeriesType.DailyPrecipitation) || 
                        (typeString == "SNOW" && timeSeriesType == TimeSeriesType.DailySnow))
                    {

                        // Extract year, month, and parse days
                        int.TryParse(line.Substring(11, 4), out var year);
                        int.TryParse(line.Substring(15, 2), out var month);
                        int daysInMonth = DateTime.DaysInMonth(year, month);

                        for (int i = 0; i < daysInMonth; i++)
                        {
                            int offset = 21 + (i * 8);
                            string strgValue = line.Substring(offset, 5).Trim();
                            var currentDate = new DateTime(year, month, i + 1);

                            // Fill missing dates with NaN if there is a gap
                            if (previousDate.HasValue && (currentDate - previousDate.Value).Days > 1)
                            {
                                FillMissingDates(timeSeries, previousDate.Value, currentDate);
                            }

                            // Parse the precipitation or snow value
                            if (int.TryParse(strgValue, out int rawValue) && rawValue != -9999)  // Valid precipitation value
                            {
                                double convertedValue = ConvertToDesiredUnit(rawValue, unit);
                                timeSeries.Add(new SeriesOrdinate<DateTime, double>(currentDate, convertedValue));
                            }
                            else  
                            {
                                // Add missing value if invalid
                                timeSeries.Add(new SeriesOrdinate<DateTime, double>(currentDate, double.NaN));
                            }

                            previousDate = currentDate;
                        }                      
                    }
                }

            }
            finally
            {
                // Ensure temporary file is deleted even if an exception occurs
                if (File.Exists(tempFilePath))
                {
                    File.Delete(tempFilePath);
                }
            }

            return timeSeries;
        }

        /// <summary>
        /// Converts raw GHCN values to the desired unit.
        /// </summary>
        /// <param name="rawValue">The raw value from the GHCN dataset.</param>
        /// <param name="unit">The unit to convert to.</param>
        /// <returns>The converted value.</returns>
        private static double ConvertToDesiredUnit(int rawValue, DepthUnit unit)
        {
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
        private static void FillMissingDates(TimeSeries timeSeries, DateTime previousDate, DateTime currentDate)
        {
            DateTime missingDate = previousDate.AddDays(1);
            while (missingDate < currentDate)
            {
                timeSeries.Add(new SeriesOrdinate<DateTime, double>(missingDate, double.NaN));
                missingDate = missingDate.AddDays(1);
            }
        }

        #endregion

        #region USGS

        /// <summary>
        /// Create the URL string for USGS download.
        /// </summary>
        /// <param name="siteNumber">USGS site number.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        /// <returns>The URL string.</returns>
        private static string CreateURLForUSGSDownload(string siteNumber, TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge)
        {
            // For annual max data (peak values)
            if (timeSeriesType == TimeSeriesType.PeakDischarge || timeSeriesType == TimeSeriesType.PeakStage)
            {
                return $"https://nwis.waterdata.usgs.gov/nwis/peak?site_no={siteNumber}&agency_cd=USGS&format=rdb";
            }

            // For field measurements (modernized OGC API)
            if (timeSeriesType == TimeSeriesType.MeasuredDischarge || timeSeriesType == TimeSeriesType.MeasuredStage)
            {
                string paramCode = timeSeriesType == TimeSeriesType.MeasuredDischarge ? "00060" : "00065";
                return $"https://api.waterdata.usgs.gov/ogcapi/v0/collections/field-measurements/items?monitoring_location_id=USGS-{siteNumber}&parameter_code={paramCode}&limit=10000&f=json";
            }

            // Daily and instantaneous data (legacy NWIS API)
            string dataTypePart = (timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.DailyStage) ? "dv/?" : "iv/?";
            string siteNumberPart = $"&sites={siteNumber}";
            string parameterCodePart = $"&parameterCd={(timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.InstantaneousDischarge ? "00060" : "00065")}";
            string startDatePart = $"&startDT={(timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.DailyStage ? "1800-01-01" : "1900-01-01")}";
            string endDatePart = $"&endDT={DateTime.Now:yyyy-MM-dd}";
            string statCodePart = (timeSeriesType == TimeSeriesType.DailyDischarge) ? "&statCd=00003" : "";
            string siteStatusPart = "&siteStatus=all";
            string formatPart = "&format=rdb";

            return $"https://waterservices.usgs.gov/nwis/{dataTypePart}{siteNumberPart}{parameterCodePart}{startDatePart}{endDatePart}{statCodePart}{siteStatusPart}{formatPart}";
        }

        /// <summary>
        /// Download time series data from USGS
        /// </summary>
        /// <param name="siteNumber">USGS site number.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        public static async Task<(TimeSeries TimeSeries, string RawText)> FromUSGS(string siteNumber, TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge)
        {
            

            // Check site number
            if (siteNumber.Length != 8)
            {
                throw new ArgumentException("The USGS site number must be 8-digits long", nameof(siteNumber));
            }

            // Check time series type
            if (timeSeriesType == TimeSeriesType.DailyPrecipitation || timeSeriesType == TimeSeriesType.DailySnow)
            {
                throw new ArgumentException("The time series type cannot be daily precipitation or daily snow.", nameof(timeSeriesType));
            }

            // Check internet connection
            if (!await IsConnectedToInternet())
            {
                throw new InvalidOperationException("No internet connection.");
            }

            var timeSeries = (timeSeriesType == TimeSeriesType.MeasuredDischarge || timeSeriesType == TimeSeriesType.MeasuredStage ||
                              timeSeriesType == TimeSeriesType.InstantaneousDischarge || timeSeriesType == TimeSeriesType.InstantaneousStage)
                ? new TimeSeries(TimeInterval.Irregular)
                : new TimeSeries();
            string textDownload = "";

            // Get URL
            string url = CreateURLForUSGSDownload(siteNumber, timeSeriesType);

            // For daily or instantaneous data, use HttpClient with gzip support
            {
                if (timeSeriesType == TimeSeriesType.DailyDischarge ||
                    timeSeriesType == TimeSeriesType.DailyStage ||
                    timeSeriesType == TimeSeriesType.InstantaneousDischarge ||
                    timeSeriesType == TimeSeriesType.InstantaneousStage)
                {
                    bool isInstantaneous = timeSeriesType == TimeSeriesType.InstantaneousDischarge ||
                                           timeSeriesType == TimeSeriesType.InstantaneousStage;
                    // Instantaneous RDB has an extra tz_cd column at index 3, pushing the value to index 4
                    int valueIndex = isInstantaneous ? 4 : 3;

                    {
                        var client = _decompressClient;

                        HttpResponseMessage response = await client.GetAsync(url);
                        response.EnsureSuccessStatusCode();

                        using (Stream contentStream = await response.Content.ReadAsStreamAsync())
                        using (StreamReader reader = new StreamReader(contentStream))
                        {
                            string? line;
                            bool isHeader = true;

                            while ((line = await reader.ReadLineAsync()) != null)
                            {
                                // Skip header row
                                if (isHeader)
                                {
                                    isHeader = false;
                                    continue;
                                }

                                // Skip comment lines
                                if (line.StartsWith("#"))
                                {
                                    if (line.Trim() == "#  No sites found matching all criteria")
                                    {
                                        throw new Exception("No data found matching all criteria.");
                                    }
                                    continue;
                                }

                                string[] fields = line.Split('\t');
                                // Validate expected number of fields and record type
                                if (fields.Length <= valueIndex || fields[0] != "USGS")
                                    continue;

                                // Parse date (assume fields[2] contains the date)
                                if (!DateTime.TryParse(fields[2], out DateTime index))
                                {
                                    continue;
                                }

                                // Fill in missing days only for daily series (not instantaneous)
                                if (!isInstantaneous && timeSeries.Count > 0 && index != TimeSeries.AddTimeInterval(timeSeries.Last().Index, TimeInterval.OneDay))
                                {
                                    while (timeSeries.Last().Index < TimeSeries.SubtractTimeInterval(index, TimeInterval.OneDay))
                                    {
                                        timeSeries.Add(new SeriesOrdinate<DateTime, double>(
                                            TimeSeries.AddTimeInterval(timeSeries.Last().Index, TimeInterval.OneDay), double.NaN));
                                    }
                                }

                                // Get and parse value
                                string valueStr = fields[valueIndex];
                                double value = string.IsNullOrWhiteSpace(valueStr) ? double.NaN : double.TryParse(valueStr, out double tempVal) ? tempVal : double.NaN;
                                timeSeries.Add(new SeriesOrdinate<DateTime, double>(index, value));
                            }
                        }
                    }
                }
                // For peak data (annual max values)
                else if (timeSeriesType == TimeSeriesType.PeakDischarge || timeSeriesType == TimeSeriesType.PeakStage)
                {
                    {
                        // Download data as string (assumes USGS peak data is not compressed)
                        textDownload = await _defaultClient.GetStringAsync(url);

                        var lines = textDownload.Split(new[] { Environment.NewLine }, StringSplitOptions.RemoveEmptyEntries);
                        foreach (string line in lines)
                        {
                            var segments = line.Split('\t');
                            if (segments.First() == "USGS" && segments.Count() >= 5)
                            {
                                // Get date
                                DateTime index = DateTime.Now;
                                int year = 2000;
                                int month = 1;
                                int day = 1;
                                var dateString = segments[2].Split('-');
                                DateTime.TryParse(segments[2], out index);

                                if (index == DateTime.MinValue)
                                {
                                    // The date parsing failed, so try to manually parse it
                                    if (dateString[1] == "00" && dateString[2] != "00")
                                    {
                                        int.TryParse(dateString[0], out year);
                                        int.TryParse(dateString[2], out day);
                                        index = new DateTime(year, month, day, 0, 0, 0);
                                    }
                                    else if (dateString[1] != "00" && dateString[2] == "00")
                                    {
                                        int.TryParse(dateString[0], out year);
                                        int.TryParse(dateString[1], out month);
                                        index = new DateTime(year, month, day, 0, 0, 0);
                                    }
                                    else if (dateString[1] == "00" && dateString[2] == "00")
                                    {
                                        int.TryParse(dateString[0], out year);
                                        index = new DateTime(year, month, day, 0, 0, 0);
                                    }
                                }
                                // Get value
                                double value = 0;
                                int idx = timeSeriesType == TimeSeriesType.PeakDischarge ? 4 : 6;
                                if (segments[idx] != "" && segments[idx] != " " && segments[idx] != "  " && !string.IsNullOrEmpty(segments[idx]))
                                {
                                    double.TryParse(segments[idx], out value);
                                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(index, value));
                                }
                            }
                        }
                    }
                }
                // For field measurements (modernized OGC API)
                else if (timeSeriesType == TimeSeriesType.MeasuredDischarge || timeSeriesType == TimeSeriesType.MeasuredStage)
                {
                    {
                        var rawText = new System.Text.StringBuilder();
                        await ParseUSGSOgcApiPages(_decompressClient, url, timeSeries, rawText);
                        textDownload = rawText.ToString();
                    }

                    timeSeries.SortByTime();
                }
            }

            return (timeSeries, textDownload);
        }

        /// <summary>
        /// Parse paginated USGS OGC API responses (GeoJSON FeatureCollection).
        /// Follows cursor-based pagination via links[rel=next].
        /// </summary>
        /// <param name="client">HttpClient to use for requests.</param>
        /// <param name="initialUrl">The initial API URL.</param>
        /// <param name="timeSeries">TimeSeries to populate with parsed data.</param>
        /// <param name="rawText">Optional StringBuilder to accumulate raw JSON responses.</param>
        private static async Task ParseUSGSOgcApiPages(HttpClient client, string initialUrl, TimeSeries timeSeries, System.Text.StringBuilder? rawText = null)
        {
            string? nextUrl = initialUrl;

            while (nextUrl != null)
            {
                // Retry with exponential backoff for rate limiting (429)
                HttpResponseMessage response = await client.GetAsync(nextUrl);
                for (int attempt = 1; attempt < 6; attempt++)
                {
                    if ((int)response.StatusCode != 429) break;
                    int delayMs = Math.Min((int)Math.Pow(2, attempt) * 2000, 60000);
                    await System.Threading.Tasks.Task.Delay(delayMs);
                    response = await client.GetAsync(nextUrl);
                }
                response.EnsureSuccessStatusCode();
                var json = await response.Content.ReadAsStringAsync();

                if (rawText != null)
                {
                    if (rawText.Length > 0) rawText.AppendLine();
                    rawText.Append(json);
                }

                var doc = System.Text.Json.JsonDocument.Parse(json);
                var root = doc.RootElement;

                if (root.TryGetProperty("features", out var features))
                {
                    foreach (var feature in features.EnumerateArray())
                    {
                        var props = feature.GetProperty("properties");
                        string? timeStr = props.GetProperty("time").GetString();
                        if (!DateTime.TryParse(timeStr, CultureInfo.InvariantCulture,
                            DateTimeStyles.AdjustToUniversal, out DateTime dt))
                            continue;

                        string? valueStr = props.GetProperty("value").GetString();
                        if (string.IsNullOrWhiteSpace(valueStr) ||
                            !double.TryParse(valueStr, NumberStyles.Float,
                            CultureInfo.InvariantCulture, out double val))
                            continue;

                        timeSeries.Add(new SeriesOrdinate<DateTime, double>(dt, val));
                    }
                }

                // Follow cursor-based pagination
                nextUrl = null;
                if (root.TryGetProperty("links", out var links))
                {
                    foreach (var link in links.EnumerateArray())
                    {
                        if (link.GetProperty("rel").GetString() is "next")
                        {
                            var candidateUrl = link.GetProperty("href").GetString();
                            // Validate the URL domain to prevent open redirect
                            if (candidateUrl != null && Uri.TryCreate(candidateUrl, UriKind.Absolute, out var uri) &&
                                uri.Host.EndsWith(".usgs.gov", StringComparison.OrdinalIgnoreCase))
                            {
                                nextUrl = candidateUrl;
                            }
                            // Small delay between pages to avoid rate limiting
                            await System.Threading.Tasks.Task.Delay(200);
                            break;
                        }
                    }
                }
            }
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
        ///     InstantaneousDischarge / InstantaneousStage → real_time_data endpoint (parameters[]=46/47).
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
        /// <returns>TimeSeries of values.</returns>
        public static async Task<TimeSeries> FromCHMN(
            string stationNumber,
            TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge,
            DischargeUnit dischargeUnit = DischargeUnit.CubicMetersPerSecond,
            HeightUnit heightUnit = HeightUnit.Meters,
            DateTime? startDate = null,
            DateTime? endDate = null)
        {
            // Connectivity
            if (!await IsConnectedToInternet())
                throw new InvalidOperationException("No internet connection.");

            // Basic validation (WSC station numbers are 7 chars including letters)
            if (string.IsNullOrWhiteSpace(stationNumber) || stationNumber.Length != 7)
                throw new ArgumentException("The WSC station number must be 7 characters, e.g., 08LG010.", nameof(stationNumber));

            // Supported series types
            if (timeSeriesType == TimeSeriesType.DailyPrecipitation || timeSeriesType == TimeSeriesType.DailySnow ||
                timeSeriesType == TimeSeriesType.MeasuredDischarge || timeSeriesType == TimeSeriesType.MeasuredStage)
                throw new ArgumentException("Canadian API supports Daily, Instantaneous, and Peak discharge/stage types. Field measurements are not available through the API.", nameof(timeSeriesType));

            // Determine if this is a discharge or stage request
            bool isDischarge = timeSeriesType == TimeSeriesType.DailyDischarge ||
                               timeSeriesType == TimeSeriesType.InstantaneousDischarge ||
                               timeSeriesType == TimeSeriesType.PeakDischarge;

            // Build the URL based on the time series type
            string url;
            TimeInterval interval;

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
                interval = TimeInterval.OneDay;
            }
            else if (timeSeriesType == TimeSeriesType.InstantaneousDischarge || timeSeriesType == TimeSeriesType.InstantaneousStage)
            {
                // Real-time data endpoint (parameter codes: 46=level, 47=discharge)
                string paramCode = isDischarge ? "47" : "46";
                DateTime sd = startDate ?? DateTime.UtcNow.AddMonths(-18);
                DateTime ed = endDate ?? DateTime.UtcNow;
                url = "https://wateroffice.ec.gc.ca/services/real_time_data/csv/inline" +
                    $"?stations[]={Uri.EscapeDataString(stationNumber)}" +
                    $"&parameters[]={paramCode}" +
                    $"&start_date={sd:yyyy-MM-dd HH:mm:ss}" +
                    $"&end_date={ed:yyyy-MM-dd HH:mm:ss}";
                interval = TimeInterval.FiveMinute;
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
                interval = TimeInterval.Irregular;
            }

            var ts = new TimeSeries(interval);

            {
                var client = _decompressClient;

                // The endpoint returns CSV (automatically decompressed by HttpClientHandler)
                using HttpResponseMessage resp = await client.GetAsync(url);
                resp.EnsureSuccessStatusCode();

                string csv;
                {
                    csv = await resp.Content.ReadAsStringAsync();
                }

                // Parse the CSV response based on endpoint type
                if (timeSeriesType == TimeSeriesType.PeakDischarge || timeSeriesType == TimeSeriesType.PeakStage)
                {
                    // Peak data CSV format:
                    // ID,Parameter/Paramètre,Date,Timezone/Fuseau horaire,Type/Catégorie,Value/Valeur,Symbol/Symbole
                    ParseCHMNPeakCsv(csv, stationNumber, isDischarge, dischargeUnit, heightUnit, ts);
                }
                else
                {
                    // Daily and real-time CSV share the same column layout:
                    // ID,Date,Parameter/Paramètre,Value/Valeur,...
                    ParseCHMNDailyCsv(csv, stationNumber, timeSeriesType, isDischarge, dischargeUnit, heightUnit, ts);
                }
            }

            return ts;
        }

        /// <summary>
        /// Parse CHMN daily or real-time CSV data.
        /// </summary>
        /// <remarks>
        /// CSV format: ID,Date,Parameter/Paramètre,Value/Valeur,... (additional columns vary by endpoint)
        /// </remarks>
        private static void ParseCHMNDailyCsv(string csv, string stationNumber, TimeSeriesType timeSeriesType,
            bool isDischarge, DischargeUnit dischargeUnit, HeightUnit heightUnit, TimeSeries ts)
        {
            bool isDailyType = timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.DailyStage;

            using (var sr = new StringReader(csv))
            {
                string? line;
                int idCol = 0, dateCol = 1, paramCol = 2, valueCol = 3;
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

                    // Parse date
                    if (!DateTime.TryParse(parts[dateCol].Trim(), CultureInfo.InvariantCulture, DateTimeStyles.None, out DateTime date))
                        continue;

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
                        FillMissingDates(ts, prevDate.Value, date);

                    ts.Add(new SeriesOrdinate<DateTime, double>(date, val));
                    prevDate = date;
                }
            }
        }

        /// <summary>
        /// Parse CHMN peak data CSV.
        /// </summary>
        /// <remarks>
        /// CSV format: ID,Parameter/Paramètre,Date,Timezone/Fuseau horaire,Type/Catégorie,Value/Valeur,Symbol/Symbole.
        /// Only "maximum" rows are included (minimums are skipped).
        /// </remarks>
        private static void ParseCHMNPeakCsv(string csv, string stationNumber, bool isDischarge,
            DischargeUnit dischargeUnit, HeightUnit heightUnit, TimeSeries ts)
        {
            using (var sr = new StringReader(csv))
            {
                string? line;
                bool headerParsed = false;

                // Peak CSV columns: ID(0), Parameter(1), Date(2), Timezone(3), Type(4), Value(5), Symbol(6)
                int idCol = 0, dateCol = 2, typeCol = 4, valueCol = 5;

                while ((line = sr.ReadLine()) != null)
                {
                    if (string.IsNullOrWhiteSpace(line)) continue;

                    // Detect header line
                    if (!headerParsed)
                    {
                        if (line.Contains("ID") && line.Contains("Date") && line.Contains("Value"))
                        {
                            headerParsed = true;
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
                    ts.Add(new SeriesOrdinate<DateTime, double>(date, val));
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
        /// <returns>TimeSeries of values.</returns>
        public static async Task<TimeSeries> FromABOM(
            string stationNumber,
            TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge,
            DischargeUnit dischargeUnit = DischargeUnit.CubicMetersPerSecond,
            HeightUnit heightUnit = HeightUnit.Meters,
            DepthUnit depthUnit = DepthUnit.Millimeters,
            DateTime? startDate = null,
            DateTime? endDate = null)
        {
            

            // Validate station number (BOM station numbers are typically 6 digits)
            if (string.IsNullOrWhiteSpace(stationNumber) || stationNumber.Length < 6)
                throw new ArgumentException("BOM station number must be at least 6 digits.", nameof(stationNumber));

            // Validate time series type - BOM supports daily/instantaneous discharge and stage, plus daily precipitation
            var supportedTypes = new[] {
                TimeSeriesType.DailyDischarge, TimeSeriesType.DailyStage,
                TimeSeriesType.InstantaneousDischarge, TimeSeriesType.InstantaneousStage,
                TimeSeriesType.DailyPrecipitation
            };
            if (!supportedTypes.Contains(timeSeriesType))
                throw new ArgumentException(
                    "BOM API supports DailyDischarge, DailyStage, InstantaneousDischarge, InstantaneousStage, and DailyPrecipitation.",
                    nameof(timeSeriesType));

            // Check connectivity
            if (!await IsConnectedToInternet())
                throw new InvalidOperationException("No internet connection.");

            // Set default dates
            DateTime sd = startDate ?? new DateTime(1800, 1, 1);
            DateTime ed = endDate ?? DateTime.Today;

            // Determine parameter type and ts_name search pattern based on requested series
            bool isDischarge = timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.InstantaneousDischarge;
            bool isStage = timeSeriesType == TimeSeriesType.DailyStage || timeSeriesType == TimeSeriesType.InstantaneousStage;
            bool isPrecip = timeSeriesType == TimeSeriesType.DailyPrecipitation;
            bool isInstantaneous = timeSeriesType == TimeSeriesType.InstantaneousDischarge || timeSeriesType == TimeSeriesType.InstantaneousStage;

            string parameterType;
            if (isDischarge) parameterType = "Water Course Discharge";
            else if (isStage) parameterType = "Water Course Level";
            else parameterType = "Rainfall";

            // Step 1: Get timeseries list to find the appropriate ts_id
            string tsListUrl = "http://www.bom.gov.au/waterdata/services" +
                $"?service=kisters&type=QueryServices&datasource=0&format=json" +
                $"&request=getTimeseriesList" +
                $"&station_no={Uri.EscapeDataString(stationNumber)}" +
                $"&parametertype_name={Uri.EscapeDataString(parameterType)}";

            string? tsId = null;

            // Create HttpClientHandler with automatic decompression
            var handler = new HttpClientHandler
            {
                AutomaticDecompression = System.Net.DecompressionMethods.GZip | System.Net.DecompressionMethods.Deflate
            };

            using (var client = new HttpClient(handler))
            {
                // Add browser-like headers to avoid security proxy issues
                client.DefaultRequestHeaders.Add("User-Agent", UserAgent);
                client.DefaultRequestHeaders.Add("Accept", "application/json,text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
                client.DefaultRequestHeaders.Add("Accept-Language", "en-US,en;q=0.9");
                client.DefaultRequestHeaders.Add("Accept-Encoding", "gzip, deflate");
                client.DefaultRequestHeaders.Add("DNT", "1");
                client.DefaultRequestHeaders.Add("Connection", "keep-alive");
                client.Timeout = TimeSpan.FromSeconds(30);

                HttpResponseMessage response;
                try
                {
                    response = await client.GetAsync(tsListUrl);
                    response.EnsureSuccessStatusCode();
                }
                catch (HttpRequestException ex)
                {
                    throw new Exception($"Failed to connect to BOM API: {ex.Message}", ex);
                }
                catch (TaskCanceledException ex)
                {
                    throw new Exception("Request to BOM API timed out.", ex);
                }

                var listResponse = await response.Content.ReadAsStringAsync();

                // Check for error response from KiWIS
                if (listResponse.Contains("\"type\":\"error\"") || listResponse.StartsWith("{\"type\":\"error\""))
                {
                    throw new Exception($"BOM API returned an error. Response: {listResponse}");
                }

                // Parse JSON response to find appropriate timeseries
                // BOM returns array where first element is column headers, subsequent elements are data rows
                System.Text.Json.JsonDocument listData;
                try
                {
                    listData = System.Text.Json.JsonDocument.Parse(listResponse);
                }
                catch (System.Text.Json.JsonException ex)
                {
                    throw new Exception($"Failed to parse BOM API response as JSON. Response: {listResponse.Substring(0, Math.Min(500, listResponse.Length))}", ex);
                }

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

                // Two-pass search: first look for DMQaQc.Merged match, then any match
                for (int pass = 0; pass < 2 && tsId == null; pass++)
                {
                    for (int i = 1; i < root.GetArrayLength(); i++)
                    {
                        var row = root[i];
                        string tsName = tsNameIndex >= 0 ? row[tsNameIndex].GetString() ?? "" : "";

                        if (!tsName.Contains(tsNamePattern)) continue;
                        if (pass == 0 && !tsName.StartsWith("DMQaQc.Merged")) continue;

                        tsId = row[tsIdIndex].GetString();
                        break;
                    }
                }

                // If no match found, take the first available series
                if (tsId == null && root.GetArrayLength() > 1)
                {
                    tsId = root[1][tsIdIndex].GetString();
                }

                if (tsId == null)
                    throw new Exception($"No suitable time series found for station {stationNumber}");
            }

            // Step 2: Get time series values using the ts_id
            string valuesUrl = "http://www.bom.gov.au/waterdata/services" +
                $"?service=kisters&type=QueryServices&datasource=0&format=json" +
                $"&request=getTimeseriesValues" +
                $"&ts_id={tsId}" +
                $"&from={sd:yyyy-MM-dd}" +
                $"&to={ed:yyyy-MM-dd}";

            var ts = isInstantaneous ? new TimeSeries(TimeInterval.Irregular) : new TimeSeries(TimeInterval.OneDay);
            DateTime? prevDate = null;

            // Create HttpClientHandler with automatic decompression
            var handler2 = new HttpClientHandler
            {
                AutomaticDecompression = System.Net.DecompressionMethods.GZip | System.Net.DecompressionMethods.Deflate
            };

            using (var client = new HttpClient(handler2))
            {
                // Add browser-like headers to avoid security proxy issues
                client.DefaultRequestHeaders.Add("User-Agent", UserAgent);
                client.DefaultRequestHeaders.Add("Accept", "application/json,text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8");
                client.DefaultRequestHeaders.Add("Accept-Language", "en-US,en;q=0.9");
                client.DefaultRequestHeaders.Add("Accept-Encoding", "gzip, deflate");
                client.DefaultRequestHeaders.Add("DNT", "1");
                client.DefaultRequestHeaders.Add("Connection", "keep-alive");
                client.Timeout = TimeSpan.FromSeconds(60);

                HttpResponseMessage response;
                try
                {
                    response = await client.GetAsync(valuesUrl);
                    response.EnsureSuccessStatusCode();
                }
                catch (HttpRequestException ex)
                {
                    throw new Exception($"Failed to retrieve time series values from BOM API: {ex.Message}", ex);
                }
                catch (TaskCanceledException ex)
                {
                    throw new Exception("Request to BOM API timed out while retrieving values.", ex);
                }

                var valuesResponse = await response.Content.ReadAsStringAsync();

                // Check for error response
                if (valuesResponse.Contains("\"type\":\"error\"") || valuesResponse.StartsWith("{\"type\":\"error\""))
                {
                    throw new Exception($"BOM API returned an error for time series values. Response: {valuesResponse}");
                }

                var valuesData = System.Text.Json.JsonDocument.Parse(valuesResponse);
                var root = valuesData.RootElement;

                // Response structure: array with single object containing 'data' array
                if (root.GetArrayLength() == 0)
                    throw new Exception("No data returned from BOM");

                var dataObj = root[0];
                if (!dataObj.TryGetProperty("data", out var dataArray))
                    throw new Exception("Invalid response format from BOM");

                // Parse the data array - each element is [timestamp, value, ...]
                foreach (var point in dataArray.EnumerateArray())
                {
                    if (point.GetArrayLength() < 2) continue;

                    // Parse timestamp
                    string? timestampStr = point[0].GetString();
                    if (!DateTime.TryParse(timestampStr, out DateTime date))
                        continue;

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
                        FillMissingDates(ts, prevDate.Value, date);

                    ts.Add(new SeriesOrdinate<DateTime, double>(date, val));
                    prevDate = date;
                }
            }

            return ts;
        }

        #endregion

    }
}
