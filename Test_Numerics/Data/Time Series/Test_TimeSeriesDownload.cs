using System;
using System.Collections.Concurrent;
using System.Linq;
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
        /// Per-service availability cache. Each entry runs its probe at most once per test
        /// run; subsequent calls reuse the cached <see cref="Lazy{T}"/> task. Used so that
        /// integration tests skip cleanly when the upstream provider is unreachable
        /// (e.g. BOM's WDP backend returning 500 / DatasourceError) instead of failing CI.
        /// </summary>
        private static readonly ConcurrentDictionary<string, Lazy<Task<bool>>> _serviceAvailability = new();

        /// <summary>
        /// Memoized "can we hit this service?" probe. Returns false when the host is offline
        /// or when the supplied probe throws (any exception is treated as "service down").
        /// </summary>
        private static Task<bool> AvailableAsync(string serviceName, Func<Task> probe) =>
            _serviceAvailability.GetOrAdd(serviceName, _ => new Lazy<Task<bool>>(async () =>
            {
                if (!await Online()) return false;
                try { await probe(); return true; }
                catch { return false; }
            })).Value;

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
        /// Tests full-period-of-record GHCN daily precipitation download.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task GHCN_FullPor_Precipitation()
        {
            if (!await GhcnAvailable()) return;
            var ts = await TimeSeriesDownload.FromGHCN(GHCN_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests full-period-of-record GHCN daily snow download.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
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

