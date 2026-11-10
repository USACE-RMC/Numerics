using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Data.TimeSeriesAnalysis
{
    /// <summary>
    /// Unit tests for the one-pass block series methods of the TimeSeries class.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Each block method is compared bit for bit with a verbatim copy of the 2.2.0 implementation, kept as a bitwise
    /// reference, which clones the series through XML, shifts the dates month by month and scans the series once per
    /// block. The one intended difference is that
    /// <see cref="TimeSeries.CustomYearSeries(int, BlockFunctionType, SmoothingFunctionType, int)"/> reports the
    /// selected ordinate's own date, where the reference reports it after a month shift that clamps the day.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [TestClass]
    public partial class Test_TimeSeriesBlockSeries
    {
        /// <summary>
        /// The water-year maximum reports the ordinate's own date: 30 November of a non-leap year, which the reference reports
        /// as 28 November.
        /// </summary>
        [TestMethod]
        public void Test_CustomYearSeries_ReportsOrdinateDate()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2021, 10, 1), new DateTime(2022, 9, 30), 1.0);
            for (int i = 0; i < ts.Count; i++)
                if (ts[i].Index == new DateTime(2021, 11, 30)) ts[i].Value = 100;
            var result = ts.CustomYearSeries(10);
            Assert.AreEqual(1, result.Count);
            Assert.AreEqual(new DateTime(2021, 11, 30), result[0].Index);
            Assert.AreEqual(100d, result[0].Value);
            Assert.AreEqual(new DateTime(2021, 11, 28), LegacyCustomYearSeries(ts, 10, BlockFunctionType.Maximum, SmoothingFunctionType.None, 1)[0].Index);
        }

        /// <summary>
        /// Small known values: hourly series over two water years, with a missing value in the second.
        /// </summary>
        [TestMethod]
        public void Test_CustomYearSeries_KnownValues()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2020, 9, 29), new double[] { 1, 5, 2, 7, double.NaN, 3 });
            // Water year 2020: 29 and 30 September. Water year 2021: 1 to 4 October.
            var max = ts.CustomYearSeries(10, BlockFunctionType.Maximum);
            Assert.AreEqual(2, max.Count);
            Assert.AreEqual(new DateTime(2020, 9, 30), max[0].Index);
            Assert.AreEqual(5d, max[0].Value);
            Assert.AreEqual(new DateTime(2020, 10, 2), max[1].Index);
            Assert.AreEqual(7d, max[1].Value);
            var sum = ts.CustomYearSeries(10, BlockFunctionType.Sum);
            // The NaN makes the 2021 sum NaN, so that year is omitted.
            Assert.AreEqual(1, sum.Count);
            Assert.AreEqual(6d, sum[0].Value);
            Assert.AreEqual(new DateTime(2020, 9, 30), sum[0].Index);
        }

        /// <summary>
        /// All block methods match the 2.2.0 implementation bit for bit on randomized series with ties, missing values,
        /// signed zeros, gaps, duplicate and unsorted dates, clamped days and every smoothing option: twelve trials, one
        /// for each combination of daily, hourly or irregular series with the four smoothing options. The 60-trial sweep
        /// with longer series runs with <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep BlockSeries</c>.
        /// </summary>
        [TestMethod]
        public void Test_BlockSeries_MatchLegacy()
        {
            MatchLegacy(12, 300, 600);
        }

        /// <summary>
        /// Block methods leave the source series as the 2.2.0 reference does: sorted only by a moving average or moving sum.
        /// </summary>
        [TestMethod]
        public void Test_BlockSeries_SourceOrderUnchanged()
        {
            var ts = new TimeSeries(TimeInterval.OneDay);
            ts.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2021, 3, 1), 2));
            ts.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2020, 3, 1), 1));
            ts.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2022, 3, 1), 3));
            ts.CustomYearSeries(10);
            ts.CalendarYearSeries();
            Assert.AreEqual(new DateTime(2021, 3, 1), ts[0].Index);
            ts.CustomYearSeries(10, BlockFunctionType.Maximum, SmoothingFunctionType.MovingAverage, 2);
            Assert.AreEqual(new DateTime(2020, 3, 1), ts[0].Index);
        }

        /// <summary>
        /// An empty series gives an empty block series, and invalid months still throw.
        /// </summary>
        [TestMethod]
        public void Test_BlockSeries_EdgeCases()
        {
            var empty = new TimeSeries(TimeInterval.OneDay);
            Assert.AreEqual(0, empty.CustomYearSeries(10).Count);
            Assert.AreEqual(0, empty.CalendarYearSeries().Count);
            Assert.AreEqual(0, empty.MonthlySeries().Count);
            Assert.AreEqual(0, empty.QuarterlySeries().Count);
            Assert.AreEqual(0, empty.CustomYearSeries(11, 3).Count);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => empty.CustomYearSeries(13));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => empty.CustomYearSeries(1, 0));
        }
    }
}
