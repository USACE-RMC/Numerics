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
            // Check internet connection
            if (!await IsConnectedToInternet())
            {
                throw new InvalidOperationException("No internet connection.");
            }

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
            catch (Exception ex)
            {
                throw ex;
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
            string statCodePart = (timeSeriesType == TimeSeriesType.DailyDischarge || timeSeriesType == TimeSeriesType.DailyStage) ? "&statCd=00003" : "";
            string siteStatusPart = "&siteStatus=all";
            string formatPart = "&format=rdb";

            return $"https://waterservices.usgs.gov/nwis/{dataTypePart}{siteNumberPart}{parameterCodePart}{startDatePart}{endDatePart}{statCodePart}{siteStatusPart}{formatPart}";
        }

        /// <summary>
        /// Download time series data from USGS
        /// </summary>
        /// <param name="siteNumber">USGS site number.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        public static async Task<TimeSeries> FromUSGS(string siteNumber, TimeSeriesType timeSeriesType = TimeSeriesType.DailyDischarge)
        {
            // Check internet connection
            if (!await IsConnectedToInternet())
            {
                throw new InvalidOperationException("No internet connection.");
            }

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

            var timeSeries = new TimeSeries();

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
                            string line;
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
                        string textDownload = await client.GetStringAsync(url);

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
            catch (Exception ex)
            {
                throw ex;
            }

            return timeSeries;
        }

        #endregion

    }

}
