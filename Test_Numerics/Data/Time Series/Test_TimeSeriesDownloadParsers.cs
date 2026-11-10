using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using static Numerics.Data.TimeSeriesDownload;

namespace Data.TimeSeriesAnalysis
{
    /// <summary>
    /// Offline tests of the <see cref="TimeSeriesDownload"/> parsers, driven by recorded provider responses.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Each fixture in <c>Fixtures/</c> is a few kilobytes cut from a real response of the USGS Water Data API, NOAA NCEI
    /// GHCN-Daily, the Water Survey of Canada or the BOM Water Data Online service. The paired <c>.expected.csv</c>
    /// holds the oracle series, computed with independent Python parsers that follow each provider's format (the
    /// csv and json modules and the GHCN-Daily readme layout). The USGS peak oracle is the official NWIS peak-flow
    /// file for the site.
    /// </para>
    /// <para>
    /// These tests run without network access and cover the parsing that <c>FromUSGS</c>, <c>FromGHCN</c>,
    /// <c>FromCHMN</c> and <c>FromABOM</c> apply to a download.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_TimeSeriesDownloadParsers
    {
        /// <summary>
        /// Reads an embedded fixture as text.
        /// </summary>
        /// <param name="name">The fixture file name.</param>
        /// <returns>The fixture text.</returns>
        internal static string Fixture(string name)
        {
            var assembly = Assembly.GetExecutingAssembly();
            string resource = assembly.GetManifestResourceNames().Single(n => n.EndsWith("Fixtures." + name, StringComparison.Ordinal));
            using var reader = new StreamReader(assembly.GetManifestResourceStream(resource)!);
            return reader.ReadToEnd();
        }

        /// <summary>
        /// Reads an oracle series: one ISO date and the hexadecimal bits of the value per line.
        /// </summary>
        /// <param name="name">The fixture name without the .expected.csv suffix.</param>
        /// <returns>The oracle dates and values.</returns>
        internal static List<(DateTime Index, double Value)> Expected(string name)
        {
            var result = new List<(DateTime, double)>();
            foreach (string line in Fixture(name + ".expected.csv").Split('\n').Skip(1))
            {
                if (line.Trim().Length == 0) continue;
                string[] parts = line.Trim().Split(',');
                result.Add((DateTime.ParseExact(parts[0], "yyyy-MM-ddTHH:mm:ss", CultureInfo.InvariantCulture),
                    BitConverter.Int64BitsToDouble(long.Parse(parts[1], NumberStyles.HexNumber, CultureInfo.InvariantCulture))));
            }
            return result;
        }

        /// <summary>
        /// Asserts that a series equals an oracle series: the same dates in order, NaN where the oracle is NaN, and
        /// otherwise the same bits.
        /// </summary>
        /// <param name="expected">The oracle series.</param>
        /// <param name="actual">The parsed series.</param>
        /// <param name="scale">A factor applied to the oracle values before comparing, or 1.</param>
        internal static void AssertSeries(List<(DateTime Index, double Value)> expected, TimeSeries actual, double scale = 1)
        {
            Assert.HasCount(expected.Count, actual);
            for (int i = 0; i < expected.Count; i++)
            {
                Assert.AreEqual(expected[i].Index, actual[i].Index, $"Date at {i}.");
                if (double.IsNaN(expected[i].Value))
                    Assert.IsTrue(double.IsNaN(actual[i].Value), $"NaN expected at {expected[i].Index:yyyy-MM-dd HH:mm}.");
                else
                    Assert.AreEqual(expected[i].Value * scale, actual[i].Value, 0d, $"Value at {expected[i].Index:yyyy-MM-dd HH:mm}.");
            }
        }

        #region USGS

        /// <summary>
        /// Daily values over a page boundary: the second page's request line and missing header are joined in order.
        /// </summary>
        [TestMethod]
        public void ParseUSGS_Daily_TwoPages_MatchesOracle()
        {
            var ts = ParseUSGS(Fixture("usgs-daily-05082500-paging.csv"), TimeSeriesType.DailyDischarge);
            AssertSeries(Expected("usgs-daily-05082500-paging"), ts);
            Assert.AreEqual(TimeInterval.OneDay, ts.TimeInterval);
            Assert.IsTrue(ts.All(o => o.Index.Kind == DateTimeKind.Unspecified));
        }

        /// <summary>
        /// Ice-affected days without a value (qualifier "ESTIMATED,ICE" in a quoted field) are NaN.
        /// </summary>
        [TestMethod]
        public void ParseUSGS_Daily_IceDays_AreNaN()
        {
            var ts = ParseUSGS(Fixture("usgs-daily-05082500-ice.csv"), TimeSeriesType.DailyStage);
            AssertSeries(Expected("usgs-daily-05082500-ice"), ts);
            Assert.AreEqual(50, ts.Count(o => double.IsNaN(o.Value)));
        }

        /// <summary>
        /// Instantaneous values across the 2024 spring daylight-saving change stay in UTC.
        /// </summary>
        [TestMethod]
        public void ParseUSGS_Instantaneous_UtcAcrossDaylightSaving()
        {
            var ts = ParseUSGS(Fixture("usgs-instantaneous-01646500.csv"), TimeSeriesType.InstantaneousDischarge);
            AssertSeries(Expected("usgs-instantaneous-01646500"), ts);
            Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
            Assert.IsTrue(ts.All(o => o.Index.Kind == DateTimeKind.Utc));
        }

        /// <summary>
        /// Annual peaks match the official NWIS peak-flow file of the site, dated in local standard time (EST).
        /// </summary>
        [TestMethod]
        public void ParseUSGS_Peaks_MatchNwisPeakFile()
        {
            var ts = ParseUSGS(Fixture("usgs-peaks-01646500.csv"), TimeSeriesType.PeakDischarge, TimeSpan.FromHours(-5));
            AssertSeries(Expected("usgs-peaks-01646500"), ts);
            Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
        }

        /// <summary>
        /// A peak with a time of day cannot be dated without the site's standard offset.
        /// </summary>
        [TestMethod]
        public void ParseUSGS_TimedPeakWithoutOffset_Throws()
        {
            Assert.Throws<ArgumentException>(() => ParseUSGS(Fixture("usgs-peaks-01646500.csv"), TimeSeriesType.PeakStage));
            string untimed = "x,y,time,value,year,month,day,time_of_day,qualifier\n,,1931-04-03,37900,1931,4,3,,\n";
            Assert.HasCount(1, ParseUSGS(untimed, TimeSeriesType.PeakDischarge));
        }

        /// <summary>
        /// Field measurements are read in UTC and sorted by time.
        /// </summary>
        [TestMethod]
        public void ParseUSGS_Measurements_SortedUtc()
        {
            var ts = ParseUSGS(Fixture("usgs-measurements-01646500.csv"), TimeSeriesType.MeasuredDischarge);
            AssertSeries(Expected("usgs-measurements-01646500"), ts);
            Assert.IsTrue(ts.All(o => o.Index.Kind == DateTimeKind.Utc));
        }

        /// <summary>
        /// Invalid arguments and empty text.
        /// </summary>
        [TestMethod]
        public void ParseUSGS_InvalidInputs()
        {
            Assert.Throws<ArgumentNullException>(() => ParseUSGS(null!, TimeSeriesType.DailyDischarge));
            Assert.Throws<ArgumentException>(() => ParseUSGS("x,y,time,value\n", TimeSeriesType.DailySnow));
            Assert.IsEmpty(ParseUSGS("", TimeSeriesType.DailyDischarge));
            Assert.IsEmpty(ParseUSGS("# comment only\nx,y,time,value\n", TimeSeriesType.InstantaneousStage));
        }

        #endregion

        #region GHCN

        /// <summary>
        /// Precipitation from a station file with other elements (SNOW, SNWD, TMAX) present, in millimeters and inches.
        /// </summary>
        [TestMethod]
        public void ParseGHCN_Precipitation_MatchesReadmeOracle()
        {
            string text = Fixture("ghcn-USW00094728.dly");
            var mm = ParseGHCN(text, TimeSeriesType.DailyPrecipitation, DepthUnit.Millimeters);
            AssertSeries(Expected("ghcn-USW00094728-prcp-mm"), mm);
            Assert.AreEqual(TimeInterval.OneDay, mm.TimeInterval);
            AssertSeries(Expected("ghcn-USW00094728-prcp-in"), ParseGHCN(text));
            var cm = ParseGHCN(text, TimeSeriesType.DailyPrecipitation, DepthUnit.Centimeters);
            for (int i = 0; i < cm.Count; i++)
                if (!double.IsNaN(cm[i].Value)) Assert.AreEqual(mm[i].Value / 10, cm[i].Value, 1e-12);
        }

        /// <summary>
        /// Snowfall is read in its native millimeters (readme: "SNOW = Snowfall (mm)"), not tenths: the 30 Dec 2000
        /// Central Park storm is 305 mm, 30.5 cm and about 12 inches.
        /// </summary>
        [TestMethod]
        public void ParseGHCN_Snow_IsMillimeters()
        {
            string text = Fixture("ghcn-USW00094728.dly");
            var mm = ParseGHCN(text, TimeSeriesType.DailySnow, DepthUnit.Millimeters);
            AssertSeries(Expected("ghcn-USW00094728-snow-mm"), mm);
            AssertSeries(Expected("ghcn-USW00094728-snow-in"), ParseGHCN(text, TimeSeriesType.DailySnow));
            var cm = ParseGHCN(text, TimeSeriesType.DailySnow, DepthUnit.Centimeters);
            var storm = new DateTime(2000, 12, 30);
            Assert.AreEqual(305d, mm.Single(o => o.Index == storm).Value);
            Assert.AreEqual(30.5, cm.Single(o => o.Index == storm).Value);
            Assert.AreEqual(305 / 25.4, ParseGHCN(text, TimeSeriesType.DailySnow).Single(o => o.Index == storm).Value);
        }

        /// <summary>
        /// The current month is padded with -9999 to its end; days after today are left out, while days up to today
        /// that are not reported yet stay NaN.
        /// </summary>
        [TestMethod]
        public void ParseGHCN_DaysAfterToday_AreLeftOut()
        {
            string text = Fixture("ghcn-USW00094728-recent.dly");
            var ts = ParseGhcnText(text, TimeSeriesType.DailyPrecipitation, DepthUnit.Millimeters, new DateTime(2026, 11, 6));
            Assert.HasCount(31 + 6, ts);
            Assert.AreEqual(new DateTime(2026, 11, 6), ts.EndDate);
            Assert.AreEqual(3.8, ts.Single(o => o.Index == new DateTime(2026, 10, 4)).Value);
            Assert.AreEqual(0d, ts.Single(o => o.Index == new DateTime(2026, 11, 4)).Value);
            Assert.IsTrue(double.IsNaN(ts.Single(o => o.Index == new DateTime(2026, 11, 5)).Value), "An unreported day up to today stays NaN.");
            Assert.HasCount(61, ParseGhcnText(text, TimeSeriesType.DailySnow, DepthUnit.Millimeters, new DateTime(2026, 12, 1)), "A finished month keeps all its days.");
            Assert.IsTrue(ParseGHCN(text).All(o => o.Index <= DateTime.Today));
        }

        /// <summary>
        /// A month missing from the station file is filled with NaN days.
        /// </summary>
        [TestMethod]
        public void ParseGHCN_MissingMonth_FilledWithNaN()
        {
            var ts = ParseGHCN(Fixture("ghcn-USC00042319.dly"), TimeSeriesType.DailyPrecipitation, DepthUnit.Millimeters);
            AssertSeries(Expected("ghcn-USC00042319-prcp-mm"), ts);
            Assert.IsTrue(ts.Where(o => o.Index.Year == 1984 && o.Index.Month == 12).All(o => double.IsNaN(o.Value)));
        }

        /// <summary>
        /// A record cut off mid-line gives a clear FormatException naming the line; blank lines are skipped.
        /// </summary>
        [TestMethod]
        public void ParseGHCN_CutRecord_ThrowsFormatException()
        {
            string[] lines = Fixture("ghcn-USW00094728.dly").Split('\n');
            string cut = lines[0] + "\n" + lines[1].Substring(0, 150);
            var e = Assert.Throws<FormatException>(() => ParseGHCN(cut, TimeSeriesType.DailyPrecipitation));
            Assert.Contains("line 2", e.Message);
            Assert.IsEmpty(ParseGHCN(lines[0] + "\n" + lines[3].Substring(0, 150), TimeSeriesType.DailySnow), "Other elements need only their first 21 characters.");
            Assert.Throws<FormatException>(() => ParseGHCN("USW0009472820", TimeSeriesType.DailyPrecipitation));
            string withBlank = lines[0] + "\n\n   \n" + lines[1] + "\n";
            Assert.AreEqual(ParseGHCN(lines[0] + "\n" + lines[1]).Count, ParseGHCN(withBlank).Count);
        }

        /// <summary>
        /// A record with a non-numeric or out-of-range year or month gives a FormatException naming the line.
        /// </summary>
        [TestMethod]
        public void ParseGHCN_InvalidYearOrMonth_ThrowsFormatException()
        {
            string record = Fixture("ghcn-USW00094728.dly").Split('\n')[1].TrimEnd('\r');
            Assert.AreEqual("PRCP", record.Substring(17, 4));
            string badYear = record.Substring(0, 11) + "XXXX" + record.Substring(15);
            string badMonth = record.Substring(0, 15) + "13" + record.Substring(17);
            Assert.Contains("line 1", Assert.Throws<FormatException>(() => ParseGHCN(badYear, TimeSeriesType.DailyPrecipitation)).Message);
            Assert.Throws<FormatException>(() => ParseGHCN(badMonth, TimeSeriesType.DailyPrecipitation));
        }

        /// <summary>
        /// A peak with a non-numeric value is skipped rather than stored as zero.
        /// </summary>
        [TestMethod]
        public void ParseUSGS_NonNumericPeak_Skipped()
        {
            var ts = ParseUSGS("time,value,year,month,day,time_of_day\n2020-01-05,abc,2020,1,5,\n2021-02-06,125,2021,2,6,\n", TimeSeriesType.PeakDischarge);
            Assert.HasCount(1, ts);
            Assert.AreEqual(125d, ts[0].Value);
        }

        /// <summary>
        /// Invalid arguments and empty text.
        /// </summary>
        [TestMethod]
        public void ParseGHCN_InvalidInputs()
        {
            Assert.Throws<ArgumentNullException>(() => ParseGHCN(null!));
            Assert.Throws<ArgumentException>(() => ParseGHCN("", TimeSeriesType.DailyStage));
            Assert.IsEmpty(ParseGHCN(""));
        }

        #endregion

        #region CHMN

        /// <summary>
        /// Daily discharge with a 16-day gap and ice symbols (B); the byte-order mark before the header is tolerated.
        /// </summary>
        [TestMethod]
        public void ParseCHMN_Daily_MatchesOracle()
        {
            string text = Fixture("chmn-05BH004-daily.csv");
            var ts = ParseCHMN(text, "05BH004");
            AssertSeries(Expected("chmn-05BH004-daily"), ts);
            Assert.AreEqual(TimeInterval.OneDay, ts.TimeInterval);
            AssertSeries(Expected("chmn-05BH004-daily"), ParseCHMN(text, "05BH004", TimeSeriesType.DailyDischarge, DischargeUnit.CubicFeetPerSecond), 35.3146667);
            Assert.IsEmpty(ParseCHMN(text, "05BH004", TimeSeriesType.DailyStage), "Discharge rows are not stage.");
            Assert.IsEmpty(ParseCHMN(text, "08MF005"), "Rows of another station are skipped.");
        }

        /// <summary>
        /// Annual peaks keep only the maximum rows.
        /// </summary>
        [TestMethod]
        public void ParseCHMN_Peaks_MatchOracle()
        {
            var ts = ParseCHMN(Fixture("chmn-08MF005-peaks.csv"), "08MF005", TimeSeriesType.PeakDischarge);
            AssertSeries(Expected("chmn-08MF005-peaks"), ts);
            Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
        }

        /// <summary>
        /// Real-time values keep the reported UTC wall-clock time across the Nov 1 daylight-saving change.
        /// </summary>
        [TestMethod]
        public void ParseCHMN_RealTime_MatchesOracle()
        {
            var ts = ParseCHMN(Fixture("chmn-08MF005-realtime.csv"), "08MF005", TimeSeriesType.InstantaneousDischarge);
            AssertSeries(Expected("chmn-08MF005-realtime"), ts);
            Assert.AreEqual(TimeInterval.FiveMinute, ts.TimeInterval);
            Assert.IsTrue(ts.All(o => o.Index.Kind == DateTimeKind.Utc), "Real-time dates are UTC.");
            Assert.AreEqual(new DateTime(2026, 11, 1, 7, 0, 0, DateTimeKind.Utc), ts[0].Index);
            var offset = ParseCHMN("ID,Date,Parameter/Paramètre,Value/Valeur\n08MF005,2024-11-01T00:00:00-08:00,47,5\n", "08MF005", TimeSeriesType.InstantaneousStage);
            Assert.AreEqual(new DateTime(2024, 11, 1, 8, 0, 0, DateTimeKind.Utc), offset[0].Index, "An explicit offset is converted to UTC.");
            Assert.IsTrue(ParseCHMN(Fixture("chmn-05BH004-daily.csv"), "05BH004").All(o => o.Index.Kind == DateTimeKind.Unspecified), "Daily dates are calendar dates.");
        }

        /// <summary>
        /// Invalid arguments.
        /// </summary>
        [TestMethod]
        public void ParseCHMN_InvalidInputs()
        {
            Assert.Throws<ArgumentNullException>(() => ParseCHMN(null!, "05BH004"));
            Assert.Throws<ArgumentException>(() => ParseCHMN("", "05BH04"));
            Assert.Throws<ArgumentException>(() => ParseCHMN("", "05BH004", TimeSeriesType.MeasuredStage));
        }

        #endregion

        #region BOM

        /// <summary>
        /// Daily discharge with a leading null value; dates are the reported +10:00 wall-clock dates.
        /// </summary>
        [TestMethod]
        public void ParseABOM_Daily_MatchesOracle()
        {
            string text = Fixture("bom-410730-daily.json");
            var ts = ParseABOM(text);
            AssertSeries(Expected("bom-410730-daily"), ts);
            Assert.AreEqual(TimeInterval.OneDay, ts.TimeInterval);
            AssertSeries(Expected("bom-410730-daily"), ParseABOM(text, TimeSeriesType.DailyDischarge, DischargeUnit.CubicFeetPerSecond), 35.3146667);
        }

        /// <summary>
        /// Instantaneous discharge across the April 5 daylight-saving change keeps the reported wall-clock times.
        /// </summary>
        [TestMethod]
        public void ParseABOM_Instantaneous_MatchesOracle()
        {
            var ts = ParseABOM(Fixture("bom-410730-instantaneous.json"), TimeSeriesType.InstantaneousDischarge);
            AssertSeries(Expected("bom-410730-instantaneous"), ts);
            Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
        }

        /// <summary>
        /// Error responses and invalid arguments.
        /// </summary>
        [TestMethod]
        public void ParseABOM_InvalidInputs()
        {
            Assert.Throws<ArgumentNullException>(() => ParseABOM(null!));
            Assert.Throws<ArgumentException>(() => ParseABOM("[]", TimeSeriesType.PeakDischarge));
            var error = Assert.Throws<Exception>(() => ParseABOM("{\"type\":\"error\",\"code\":\"InvalidParameterValue\"}"));
            Assert.Contains("returned an error", error.Message);
            Assert.Throws<Exception>(() => ParseABOM("[]"));
            Assert.Throws<Exception>(() => ParseABOM("[{\"ts_id\":\"1\"}]"));
        }

        /// <summary>
        /// The station lookup prefers quality-controlled merged series for each type.
        /// </summary>
        [TestMethod]
        public void SelectABOMTimeSeriesId_PrefersMergedSeries()
        {
            string discharge = Fixture("bom-410730-discharge-list.json");
            Assert.AreEqual("1573010", SelectABOMTimeSeriesId(discharge, "410730", TimeSeriesType.DailyDischarge));
            Assert.AreEqual("1570010", SelectABOMTimeSeriesId(discharge, "410730", TimeSeriesType.InstantaneousDischarge));
            Assert.AreEqual("1554010", SelectABOMTimeSeriesId(Fixture("bom-410730-rainfall-list.json"), "410730", TimeSeriesType.DailyPrecipitation));
            var none = Assert.Throws<Exception>(() => SelectABOMTimeSeriesId("[[\"ts_id\",\"ts_name\"],[\"1\",\"X.AsStored.1\"]]", "410730", TimeSeriesType.DailyStage));
            Assert.Contains("No suitable time series", none.Message);
        }

        /// <summary>
        /// Daily discharge, stage and rainfall take the .24HR series even when a .09HR series is listed first; without
        /// a .24HR series the first match is kept.
        /// </summary>
        [TestMethod]
        public void SelectABOMTimeSeriesId_DailyTakes24HourSeries()
        {
            const string head = "[[\"ts_id\",\"ts_name\"],";
            Assert.AreEqual("24", SelectABOMTimeSeriesId(head + "[\"9\",\"DMQaQc.Merged.DailyMean.09HR\"],[\"24\",\"DMQaQc.Merged.DailyMean.24HR\"]]", "410730", TimeSeriesType.DailyDischarge));
            Assert.AreEqual("24", SelectABOMTimeSeriesId(head + "[\"9\",\"DMQaQc.Merged.DailyMean.09HR\"],[\"24\",\"PR01QaQc.Merged.DailyMean.24HR\"]]", "410730", TimeSeriesType.DailyStage));
            Assert.AreEqual("24", SelectABOMTimeSeriesId(head + "[\"9\",\"DMQaQc.Merged.DailyTotal.09HR\"],[\"24\",\"DMQaQc.Merged.DailyTotal.24HR\"]]", "410730", TimeSeriesType.DailyPrecipitation));
            Assert.AreEqual("9", SelectABOMTimeSeriesId(head + "[\"8\",\"PR01QaQc.Merged.DailyMean.09HR\"],[\"9\",\"DMQaQc.Merged.DailyMean.09HR\"]]", "410730", TimeSeriesType.DailyDischarge));
            Assert.AreEqual("a", SelectABOMTimeSeriesId(head + "[\"a\",\"DMQaQc.Merged.AsStored.1\"],[\"b\",\"X.AsStored.24HR\"]]", "410730", TimeSeriesType.InstantaneousStage), "Instantaneous types are unchanged.");
        }

        #endregion

        #region Qualifiers

        /// <summary>
        /// Asserts that the qualifier overload returns the same series bits as the default parser.
        /// </summary>
        /// <param name="plain">The series from the default parser.</param>
        /// <param name="withQualifiers">The series from the qualifier overload.</param>
        /// <param name="qualifiers">The qualifier list.</param>
        private static void AssertSameSeries(TimeSeries plain, TimeSeries withQualifiers, IList<string> qualifiers)
        {
            Assert.HasCount(plain.Count, qualifiers);
            CollectionAssert.AreEqual(plain.Select(o => o.Index).ToArray(), withQualifiers.Select(o => o.Index).ToArray());
            CollectionAssert.AreEqual(plain.Select(o => BitConverter.DoubleToInt64Bits(o.Value)).ToArray(),
                withQualifiers.Select(o => BitConverter.DoubleToInt64Bits(o.Value)).ToArray());
            Assert.AreEqual(plain.TimeInterval, withQualifiers.TimeInterval);
        }

        /// <summary>
        /// USGS qualifiers are the qualifier field; the series is unchanged.
        /// </summary>
        [TestMethod]
        public void ParseUSGS_Qualifiers_ParallelToSeries()
        {
            string text = Fixture("usgs-daily-05082500-ice.csv");
            var q = new List<string> { "stale" };
            var ts = ParseUSGS(text, TimeSeriesType.DailyDischarge, null, q);
            AssertSameSeries(ParseUSGS(text, TimeSeriesType.DailyDischarge), ts, q);
            Assert.AreEqual("ESTIMATED,ICE", q[ts.ToList().FindIndex(o => o.Index == new DateTime(2026, 2, 14))]);
            Assert.AreEqual(17, q.Count(s => s == "ESTIMATED,ICE"));
            Assert.AreEqual(33, q.Count(s => s == "ICE"));
            Assert.DoesNotContain("stale", q, "The list is cleared first.");

            var peaks = new List<string>();
            ParseUSGS(Fixture("usgs-peaks-01646500.csv"), TimeSeriesType.PeakDischarge, TimeSpan.FromHours(-5), peaks);
            Assert.HasCount(95, peaks);
            Assert.Throws<ArgumentNullException>(() => ParseUSGS(text, TimeSeriesType.DailyDischarge, null, null!));
        }

        /// <summary>
        /// GHCN qualifiers are the day's three flag characters (MFLAG, QFLAG, SFLAG).
        /// </summary>
        [TestMethod]
        public void ParseGHCN_Qualifiers_AreDayFlags()
        {
            string text = Fixture("ghcn-USW00094728.dly");
            var q = new List<string>();
            var ts = ParseGHCN(text, TimeSeriesType.DailySnow, DepthUnit.Millimeters, q);
            AssertSameSeries(ParseGHCN(text, TimeSeriesType.DailySnow, DepthUnit.Millimeters), ts, q);
            Assert.AreEqual(" I0", q[ts.ToList().FindIndex(o => o.Index == new DateTime(2001, 1, 23))], "Failed internal-consistency check.");
            Assert.IsTrue(q.All(s => s.Length == 3));

            var gap = new List<string>();
            var gapped = ParseGHCN(Fixture("ghcn-USC00042319.dly"), TimeSeriesType.DailyPrecipitation, DepthUnit.Millimeters, gap);
            Assert.IsTrue(gapped.Select((o, i) => (o, i)).Where(x => x.o.Index.Year == 1984 && x.o.Index.Month == 12).All(x => gap[x.i] == ""), "Filled days have no flags.");
        }

        /// <summary>
        /// CHMN qualifiers are the Symbol/Symbole field (B for ice conditions).
        /// </summary>
        [TestMethod]
        public void ParseCHMN_Qualifiers_AreSymbols()
        {
            string text = Fixture("chmn-05BH004-daily.csv");
            var q = new List<string>();
            var ts = ParseCHMN(text, "05BH004", TimeSeriesType.DailyDischarge, DischargeUnit.CubicMetersPerSecond, HeightUnit.Meters, q);
            AssertSameSeries(ParseCHMN(text, "05BH004"), ts, q);
            Assert.AreEqual("B", q[0]);
            Assert.AreEqual("", q[ts.ToList().FindIndex(o => o.Index == new DateTime(1911, 1, 10))], "A filled gap day has no symbol.");

            var rt = new List<string>();
            ParseCHMN(Fixture("chmn-08MF005-realtime.csv"), "08MF005", TimeSeriesType.InstantaneousDischarge, DischargeUnit.CubicMetersPerSecond, HeightUnit.Meters, rt);
            Assert.HasCount(30, rt);
            var pk = new List<string>();
            ParseCHMN(Fixture("chmn-08MF005-peaks.csv"), "08MF005", TimeSeriesType.PeakDischarge, DischargeUnit.CubicMetersPerSecond, HeightUnit.Meters, pk);
            Assert.HasCount(74, pk);
        }

        /// <summary>
        /// BOM qualifiers are the KiWIS Quality Code when the response has that column, and empty otherwise.
        /// </summary>
        [TestMethod]
        public void ParseABOM_Qualifiers_AreQualityCodes()
        {
            string text = Fixture("bom-410730-daily-quality.json");
            var q = new List<string>();
            var ts = ParseABOM(text, TimeSeriesType.DailyDischarge, DischargeUnit.CubicMetersPerSecond, HeightUnit.Meters, DepthUnit.Millimeters, q);
            AssertSameSeries(ParseABOM(text), ts, q);
            CollectionAssert.AreEqual(new[] { "10", "10", "10", "90" }, q.Take(4).ToArray());
            Assert.AreEqual("-1", q[39]);
            Assert.IsTrue(double.IsNaN(ts[39].Value));

            var none = new List<string>();
            ParseABOM(Fixture("bom-410730-daily.json"), TimeSeriesType.DailyDischarge, DischargeUnit.CubicMetersPerSecond, HeightUnit.Meters, DepthUnit.Millimeters, none);
            Assert.IsTrue(none.All(s => s == ""), "A response without a Quality Code column has no codes.");
        }

        #endregion
    }
}
