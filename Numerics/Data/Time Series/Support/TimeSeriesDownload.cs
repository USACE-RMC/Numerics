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

        /// <summary>
        /// Checks if there is an Internet connection.
        /// </summary>
        public static async Task<bool> IsConnectedToInternet()
        {
            try
            {
                using (HttpClient client = new HttpClient())
                {
                    client.Timeout = TimeSpan.FromSeconds(5);
                    await client.GetAsync("https://www.google.com");
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
            PeakStage
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
            if (siteNumber.Length != 11)
            {
                throw new ArgumentException("The GHCN site number must be 11-digits long", nameof(siteNumber));
            }

            // Check time series type
            if (timeSeriesType != TimeSeriesType.DailyPrecipitation && timeSeriesType != TimeSeriesType.DailySnow)
            {
                throw new ArgumentException("The time series type must be either daily precipitation or daily snow.", nameof(timeSeriesType));
            }

            var timeSeries = new TimeSeries(TimeInterval.OneDay);
            DateTime? previousDate = null;
            string tempFilePath = Path.Combine(Path.GetTempPath(), $"{siteNumber}.dly");

            // Check internet connection
            if (!await IsConnectedToInternet())
            {
                throw new InvalidOperationException("No internet connection.");
            }

            try
            {
 
                // Download the GHCN file
                string ghcnBaseUrl = "https://www.ncei.noaa.gov/pub/data/ghcn/daily/all/";
                string stationFileUrl = $"{ghcnBaseUrl}{siteNumber}.dly";
                using (var client = new HttpClient())
                {
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
            catch (Exception)
            {
                throw;
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

            // Determine URL parts for daily and instantaneous data
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

            var timeSeries = new TimeSeries();
            string textDownload = "";

            // Get URL
            string url = CreateURLForUSGSDownload(siteNumber, timeSeriesType);

            try
            {
                // For daily or instantaneous data, use HttpClient with gzip support
                if (timeSeriesType == TimeSeriesType.DailyDischarge || 
                    timeSeriesType == TimeSeriesType.DailyStage || 
                    timeSeriesType == TimeSeriesType.InstantaneousDischarge || 
                    timeSeriesType == TimeSeriesType.InstantaneousStage)
                {
                    using (HttpClient client = new HttpClient())
                    {
                        // Set request headers to accept gzip encoding
                        client.DefaultRequestHeaders.AcceptEncoding.Add(new System.Net.Http.Headers.StringWithQualityHeaderValue("gzip"));

                        HttpResponseMessage response = await client.GetAsync(url);
                        response.EnsureSuccessStatusCode();

                        // Ensure the response content is compressed as expected
                        if (!response.Content.Headers.ContentEncoding.Contains("gzip"))
                        {
                            throw new Exception("Response is not compressed as expected.");
                        }

                        using (Stream compressedStream = await response.Content.ReadAsStreamAsync())
                        using (GZipStream decompressionStream = new GZipStream(compressedStream, CompressionMode.Decompress))
                        using (StreamReader reader = new StreamReader(decompressionStream))
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
                                if (fields.Length < 5 || fields[0] != "USGS")
                                    continue;

                                // Parse date (assume fields[2] contains the date)
                                if (!DateTime.TryParse(fields[2], out DateTime index))
                                {
                                    // Optionally log or handle date parse failures
                                    continue;
                                }

                                // Fill in missing days if the time series is not continuous
                                if (timeSeries.Count > 0 && index != TimeSeries.AddTimeInterval(timeSeries.Last().Index, TimeInterval.OneDay))
                                {
                                    while (timeSeries.Last().Index < TimeSeries.SubtractTimeInterval(index, TimeInterval.OneDay))
                                    {
                                        timeSeries.Add(new SeriesOrdinate<DateTime, double>(
                                            TimeSeries.AddTimeInterval(timeSeries.Last().Index, TimeInterval.OneDay), double.NaN));
                                    }
                                }

                                // Get and parse value
                                string valueStr = fields[3];
                                double value = string.IsNullOrWhiteSpace(valueStr) ? double.NaN : double.TryParse(valueStr, out double tempVal) ? tempVal : double.NaN;
                                timeSeries.Add(new SeriesOrdinate<DateTime, double>(index, value));
                            }
                        }
                    }
                }
                // For peak data (annual max values)
                else if (timeSeriesType == TimeSeriesType.PeakDischarge || timeSeriesType == TimeSeriesType.PeakStage)
                {
                    using (HttpClient client = new HttpClient())
                    {
                        // Download data as string (assumes USGS peak data is not compressed)
                        textDownload = await client.GetStringAsync(url);

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
            }
            catch (Exception)
            {
                throw;
            }

            return (timeSeries, textDownload);
        }

        #endregion

        #region Canadian Hydrometric Monitoring Network (CHMN)

        /// <summary>
        /// Download historical daily means from the Water Survey of Canada
        /// (Environment and Climate Change Canada Hydrometric Monitoring Network).
        /// </summary>
        /// <param name="stationNumber">WSC 7-character station id, e.g., "08LG010".</param>
        /// <param name="timeSeriesType">
        ///     DailyDischarge → parameters[]=flow
        ///     DailyStage     → parameters[]=level
        /// </param>
        /// <param name="dischargeUnit">Desired discharge unit if DailyDischarge.</param>
        /// <param name="heightUnit">Desired stage unit if DailyStage.</param>
        /// <param name="startDate">
        ///     Optional inclusive start date. If null, defaults to 1800-01-01
        ///     to request full POR. The service returns only available days.
        /// </param>
        /// <param name="endDate">
        ///     Optional inclusive end date. If null, defaults to DateTime.Today.
        /// </param>
        /// <returns>TimeSeries of daily means.</returns>
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

            // Supported series types for this endpoint
            if (timeSeriesType != TimeSeriesType.DailyDischarge && timeSeriesType != TimeSeriesType.DailyStage)
                throw new ArgumentException("Canadian daily_data API supports DailyDischarge or DailyStage.", nameof(timeSeriesType));

            // Parameter mapping per WSC docs: parameters[]=flow or parameters[]=level
            string parameter = timeSeriesType == TimeSeriesType.DailyDischarge ? "flow" : "level";

            // Date range: request wide window to effectively get full POR
            DateTime sd = startDate ?? new DateTime(1800, 1, 1);
            DateTime ed = endDate ?? DateTime.Today;

            string url =
                "https://wateroffice.ec.gc.ca/services/daily_data/csv/inline" +
                $"?stations[]={Uri.EscapeDataString(stationNumber)}" +
                $"&parameters[]={parameter}" +
                $"&start_date={sd:yyyy-MM-dd}" +
                $"&end_date={ed:yyyy-MM-dd}";

            var ts = new TimeSeries(TimeInterval.OneDay);

            // Create HttpClientHandler to bypass proxy if needed
            var handler = new HttpClientHandler
            {
                UseProxy = false,  // Bypass corporate proxy
                AutomaticDecompression = System.Net.DecompressionMethods.GZip | System.Net.DecompressionMethods.Deflate
            };

            using (var client = new HttpClient(handler))
            {
                // Add browser-like headers to avoid security proxy issues
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36");
                client.DefaultRequestHeaders.Add("Accept", "text/csv,application/csv,text/plain,*/*");
                client.DefaultRequestHeaders.Add("Accept-Language", "en-US,en;q=0.9");
                client.DefaultRequestHeaders.Add("Accept-Encoding", "gzip, deflate, br");
                client.DefaultRequestHeaders.Add("DNT", "1");
                client.DefaultRequestHeaders.Add("Connection", "keep-alive");
                client.DefaultRequestHeaders.Add("Upgrade-Insecure-Requests", "1");

                // The endpoint returns CSV (possibly gzip-compressed)
                using HttpResponseMessage resp = await client.GetAsync(url);
                resp.EnsureSuccessStatusCode();

                string csv;
                byte[] responseBytes = await resp.Content.ReadAsByteArrayAsync();

                // Check if the response is gzip-compressed by checking the magic number
                // Gzip magic number is 0x1f 0x8b
                bool isGzipped = responseBytes.Length >= 2 && responseBytes[0] == 0x1f && responseBytes[1] == 0x8b;

                if (isGzipped)
                {
                    // Decompress gzip data
                    using (var compressedStream = new MemoryStream(responseBytes))
                    using (var gzipStream = new GZipStream(compressedStream, CompressionMode.Decompress))
                    using (var reader = new StreamReader(gzipStream))
                    {
                        csv = await reader.ReadToEndAsync();
                    }
                }
                else
                {
                    // Read as plain text
                    csv = System.Text.Encoding.UTF8.GetString(responseBytes);
                }

                // Parse the CSV response
                using (var sr = new StringReader(csv))
                {
                    string? line;
                    int idCol = 0, dateCol = 1, paramCol = 2, valueCol = 3;
                    bool headerParsed = false;

                    // For filling missing dates
                    DateTime? prevDate = null;

                    while ((line = sr.ReadLine()) != null)
                    {
                        // Skip empty lines
                        if (string.IsNullOrWhiteSpace(line)) continue;

                        // Try to detect and parse the header
                        if (!headerParsed)
                        {

                            // If this line contains the key header terms, treat it as the header
                            if (line == "﻿ ID,Date,Parameter/Paramètre,Value/Valeur,Symbol/Symbole")
                            {
                                headerParsed = true;
                            }

                            // If we haven't found a header yet, skip this line (could be metadata)
                            continue;
                        }

                        // We're past the header - parse data rows
                        string[] parts = line.Split(',');

                        // Validate row has enough columns
                        if (parts.Length <= Math.Max(Math.Max(idCol, dateCol), valueCol))
                            continue;

                        // Check station match
                        if (!parts[idCol].Trim().Equals(stationNumber, StringComparison.OrdinalIgnoreCase))
                            continue;

                        // Check parameter match (parameter column contains "discharge/débit" or "level")
                        if (paramCol >= 0 && parts.Length > paramCol)
                        {
                            string paramValue = parts[paramCol].Trim().ToLowerInvariant();
                            bool isFlowParam = paramValue.Contains("discharge") || paramValue.Contains("débit") || paramValue.Contains("flow");
                            bool isLevelParam = paramValue.Contains("level") || paramValue.Contains("stage") || paramValue.Contains("niveau");

                            // Skip if parameter doesn't match requested type
                            if (timeSeriesType == TimeSeriesType.DailyDischarge && !isFlowParam)
                                continue;
                            if (timeSeriesType == TimeSeriesType.DailyStage && !isLevelParam)
                                continue;
                        }

                        // Parse date (handle both M/d/yyyy and yyyy-MM-dd formats)
                        if (!DateTime.TryParse(parts[dateCol].Trim(), out DateTime date))
                            continue;

                        // Parse value
                        double val = double.NaN;
                        string valueStr = parts[valueCol].Trim();

                        if (!string.IsNullOrWhiteSpace(valueStr))
                        {
                            // WSC uses dot decimal, values in m^3/s for flow and m for level
                            if (double.TryParse(valueStr, NumberStyles.Float, CultureInfo.InvariantCulture, out double raw))
                            {
                                if (timeSeriesType == TimeSeriesType.DailyDischarge)
                                {
                                    val = dischargeUnit == DischargeUnit.CubicFeetPerSecond
                                        ? raw * 35.3146667   // cms → cfs
                                        : raw;               // cms
                                }
                                else
                                {
                                    val = heightUnit == HeightUnit.Feet
                                        ? raw * 3.280839895 // m → ft
                                        : raw;              // m
                                }
                            }
                        }

                        // Fill gaps with NaN to keep a continuous daily series
                        if (prevDate.HasValue && (date - prevDate.Value).Days > 1)
                            FillMissingDates(ts, prevDate.Value, date);

                        ts.Add(new SeriesOrdinate<DateTime, double>(date, val));
                        prevDate = date;
                    }
                }
            }

            return ts;
        }


        #endregion

        #region Australian Bureau of Meteorology (BOM)

        /// <summary>
        /// Download daily discharge or stage data from the Australian Bureau of Meteorology (BOM) 
        /// via the KiWIS API.
        /// </summary>
        /// <param name="stationNumber">BOM station number (e.g., "410730").</param>
        /// <param name="timeSeriesType">The time series type (DailyDischarge or DailyStage).</param>
        /// <param name="dischargeUnit">Desired discharge unit if DailyDischarge.</param>
        /// <param name="heightUnit">Desired stage unit if DailyStage.</param>
        /// <param name="startDate">Optional start date. If null, attempts to retrieve full period of record.</param>
        /// <param name="endDate">Optional end date. If null, defaults to today.</param>
        /// <returns>TimeSeries of daily values.</returns>
        public static async Task<TimeSeries> FromABOM(
            string stationNumber,
            TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge,
            DischargeUnit dischargeUnit = DischargeUnit.CubicMetersPerSecond,
            HeightUnit heightUnit = HeightUnit.Meters,
            DateTime? startDate = null,
            DateTime? endDate = null)
        {
            

            // Validate station number (BOM station numbers are typically 6 digits)
            if (string.IsNullOrWhiteSpace(stationNumber) || stationNumber.Length < 6)
                throw new ArgumentException("BOM station number must be at least 6 digits.", nameof(stationNumber));

            // Validate time series type
            if (timeSeriesType != TimeSeriesType.DailyDischarge && timeSeriesType != TimeSeriesType.DailyStage)
                throw new ArgumentException("BOM API supports DailyDischarge or DailyStage only.", nameof(timeSeriesType));

            // Check connectivity
            if (!await IsConnectedToInternet())
                throw new InvalidOperationException("No internet connection.");

            // Set default dates
            DateTime sd = startDate ?? new DateTime(1800, 1, 1);
            DateTime ed = endDate ?? DateTime.Today;

            // Determine parameter type based on requested series
            string parameterType = timeSeriesType == TimeSeriesType.DailyDischarge
                ? "Water Course Discharge"
                : "Water Course Level";

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
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36");
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

                // Look for daily mean time series (prioritize merged/quality-checked data)
                for (int i = 1; i < root.GetArrayLength(); i++)
                {
                    var row = root[i];
                    string? tsName = tsNameIndex >= 0 ? row[tsNameIndex].GetString() : "";

                    if (tsName == null) continue;
                    // Prioritize: DMQaQc.Merged.DailyMean.24HR or similar daily mean series
                    if (tsName.Contains("DailyMean") || tsName.Contains("Daily Mean"))
                    {
                        tsId = row[tsIdIndex].GetString();
                        break;
                    }
                }

                // If no daily mean found, take the first available series
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

            var ts = new TimeSeries(TimeInterval.OneDay);
            DateTime? prevDate = null;

            // Create HttpClientHandler with automatic decompression
            var handler2 = new HttpClientHandler
            {
                AutomaticDecompression = System.Net.DecompressionMethods.GZip | System.Net.DecompressionMethods.Deflate
            };

            using (var client = new HttpClient(handler2))
            {
                // Add browser-like headers to avoid security proxy issues
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36");
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
                        if (timeSeriesType == TimeSeriesType.DailyDischarge)
                        {
                            // BOM returns discharge in m³/s
                            val = dischargeUnit == DischargeUnit.CubicFeetPerSecond
                                ? raw * 35.3146667
                                : raw;
                        }
                        else
                        {
                            // BOM returns stage in meters
                            val = heightUnit == HeightUnit.Feet
                                ? raw * 3.280839895
                                : raw;
                        }
                    }

                    // Fill gaps with NaN
                    if (prevDate.HasValue && (date - prevDate.Value).Days > 1)
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
