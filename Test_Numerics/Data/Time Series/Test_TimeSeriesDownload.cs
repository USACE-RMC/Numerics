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

using System;
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
        /// Verifies that a time series has monotonically increasing, non-duplicated date indices.
        /// </summary>
        /// <param name="ts">The time series to validate.</param>
        private static void AssertDailySeriesMonotonic(TimeSeries ts)
        {
            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan( 0, ts.Count);

            DateTime? prev = null;
            foreach (var pt in ts)
            {
                if (prev.HasValue)
                    Assert.IsTrue(pt.Index >= prev.Value, "Dates not sorted chronologically.");
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
            if (!await Online()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the CHMN Lillooet River station (flow).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task CHMN_FullPor_Lillooet_Flow()
        {
            if (!await Online()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_2);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the CHMN Capilano River station (flow).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task CHMN_FullPor_Capilano_Flow()
        {
            if (!await Online()) return;
            var ts = await TimeSeriesDownload.FromCHMN(CHMN_3);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests flow unit conversions (cms ↔ cfs) for CHMN data.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task CHMN_UnitConversion_Flow_CmsCfs()
        {
            if (!await Online()) return;

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
            if (!await Online()) return;

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
        /// Ensures CHMN rejects unsupported time series types.
        /// </summary>
        [TestMethod]
        public async Task CHMN_UnsupportedType_Throws()
        {
            await Assert.ThrowsAsync<ArgumentException>(async () =>
                await TimeSeriesDownload.FromCHMN(CHMN_1,
                    TimeSeriesDownload.TimeSeriesType.PeakDischarge));
        }

        #endregion

        #region USGS Tests

        /// <summary>
        /// Tests full-period-of-record USGS daily discharge download.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task USGS_FullPor_DailyDischarge()
        {
            if (!await Online()) return;

            var (ts, raw) = await TimeSeriesDownload.FromUSGS(USGS_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests USGS daily stage download for correctness and continuity.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task USGS_FullPor_DailyStage()
        {
            if (!await Online()) return;

            var (ts, _) = await TimeSeriesDownload.FromUSGS(USGS_4, TimeSeriesDownload.TimeSeriesType.DailyStage);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests USGS peak discharge data retrieval for non-daily data.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task USGS_PeakDischarge_Works()
        {
            if (!await Online()) return;

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

        #endregion

        #region GHCN Tests

        /// <summary>
        /// Tests full-period-of-record GHCN daily precipitation download.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task GHCN_FullPor_Precipitation()
        {
            if (!await Online()) return;
            var ts = await TimeSeriesDownload.FromGHCN(GHCN_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests full-period-of-record GHCN daily snow download.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task GHCN_FullPor_Snow()
        {
            if (!await Online()) return;
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
            if (!await Online()) return;

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
            if (!await Online()) return;

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
            if (!await Online()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_1);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the BOM Goodradigbee River station (discharge).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task BOM_FullPor_Goodradigbee_Discharge()
        {
            if (!await Online()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_2);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Validates a full-period-of-record download for the BOM Murray River station (stage).
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task BOM_FullPor_MurrayRiver_Stage()
        {
            if (!await Online()) return;
            var ts = await TimeSeriesDownload.FromABOM(BOM_3, TimeSeriesDownload.TimeSeriesType.DailyStage);
            AssertDailySeriesMonotonic(ts);
        }

        /// <summary>
        /// Tests discharge unit conversions (cms ↔ cfs) for BOM data.
        /// </summary>
        [TestMethod, TestCategory("Integration")]
        public async Task BOM_UnitConversion_Discharge_CmsCfs()
        {
            if (!await Online()) return;

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
            if (!await Online()) return;

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
            if (!await Online()) return;

            var ts = await TimeSeriesDownload.FromABOM(BOM_1,
                TimeSeriesDownload.TimeSeriesType.DailyDischarge,
                startDate: WinStart, endDate: WinEnd);

            Assert.IsNotNull(ts, "Time series is null.");
            Assert.IsGreaterThan(0, ts.Count);

            // Verify data is within requested window (allowing for some timezone flexibility)
            var firstDate = ts.First().Index;
            var lastDate = ts.Last().Index;

            Assert.IsTrue(firstDate >= WinStart.AddDays(-1),
                $"First date {firstDate} is before window start {WinStart}");
            Assert.IsTrue(lastDate <= WinEnd.AddDays(1),
                $"Last date {lastDate} is after window end {WinEnd}");
        }

        #endregion
    }
}

