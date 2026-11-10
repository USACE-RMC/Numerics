#nullable disable
using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;

namespace Data.TimeSeriesAnalysis
{
    /// <summary>
    /// The comparison of the block series methods with verbatim copies of the 2.2.0 implementations, shared by the unit
    /// test and the full sweep in Benchmarks_Numerics.
    /// </summary>
    public partial class Test_TimeSeriesBlockSeries
    {
        private static readonly BlockFunctionType[] BlockFunctions = { BlockFunctionType.Minimum, BlockFunctionType.Maximum, BlockFunctionType.Average, BlockFunctionType.Sum };

        /// <summary>
        /// Compares every block method with the 2.2.0 implementation bit for bit on randomized series.
        /// </summary>
        /// <param name="trials">The number of random series; trials 0 to 11 cover each kind of series with each
        /// smoothing option.</param>
        /// <param name="dailyLength">The bound on the length of daily and irregular series, less 20.</param>
        /// <param name="hourlyLength">The bound on the length of hourly series, less 20.</param>
        public static void MatchLegacy(int trials, int dailyLength, int hourlyLength)
        {
            var rnd = new Random(20261105);
            for (int trial = 0; trial < trials; trial++)
            {
                var ts = RandomSeries(rnd, trial, dailyLength, hourlyLength);
                var smoothing = (SmoothingFunctionType)(trial % 4);
                int period = trial % 3 == 0 ? 1 : 1 + rnd.Next(4);
                if (period >= ts.Count) period = 1;
                foreach (var f in BlockFunctions)
                {
                    AssertSame(() => LegacyCalendarYearSeries(ts.Clone(), f, smoothing, period), () => ts.Clone().CalendarYearSeries(f, smoothing, period));
                    AssertSame(() => LegacyMonthlySeries(ts.Clone(), f, smoothing, period), () => ts.Clone().MonthlySeries(f, smoothing, period));
                    AssertSame(() => LegacyQuarterlySeries(ts.Clone(), f, smoothing, period), () => ts.Clone().QuarterlySeries(f, smoothing, period));
                    int start = 1 + rnd.Next(12), end = 1 + rnd.Next(12);
                    AssertSame(() => LegacyCustomYearSeries(ts.Clone(), start, f, smoothing, period), () => ts.Clone().CustomYearSeries(start, f, smoothing, period), start);
                    AssertSame(() => LegacyCustomYearSeries(ts.Clone(), 10, f, smoothing, period), () => ts.Clone().CustomYearSeries(10, f, smoothing, period), 10);
                    AssertSame(() => LegacyCustomYearSeries(ts.Clone(), start, end, f, smoothing, period), () => ts.Clone().CustomYearSeries(start, end, f, smoothing, period));
                }
            }
        }

        /// <summary>
        /// Builds a random daily, hourly or irregular series with the cases the legacy comparison must cover.
        /// </summary>
        /// <param name="rnd">The random generator.</param>
        /// <param name="trial">The trial number, which selects the kind of series.</param>
        /// <param name="dailyLength">The bound on the length of daily and irregular series, less 20.</param>
        /// <param name="hourlyLength">The bound on the length of hourly series, less 20.</param>
        private static TimeSeries RandomSeries(Random rnd, int trial, int dailyLength, int hourlyLength)
        {
            int kind = trial % 3;
            int n = 20 + rnd.Next(kind == 1 ? hourlyLength : dailyLength);
            var start = new DateTime(1999 + rnd.Next(3), 1 + rnd.Next(12), 1 + rnd.Next(28), rnd.Next(24), 0, 0);
            var values = new double[n];
            for (int i = 0; i < n; i++)
            {
                double u = rnd.NextDouble();
                values[i] = u < 0.05 ? double.NaN : u < 0.08 ? -0.0 : u < 0.11 ? 0.0 : u < 0.3 ? Math.Round(rnd.NextDouble() * 4) : rnd.NextDouble() * 100 - 20;
            }
            if (kind == 0) return new TimeSeries(TimeInterval.OneDay, start, values);
            if (kind == 1) return new TimeSeries(TimeInterval.OneHour, start, values);

            // Irregular: month-end dates (clamped by the shift), gaps, duplicates and some out-of-order ordinates.
            var ts = new TimeSeries(TimeInterval.Irregular);
            var date = new DateTime(2000 + rnd.Next(4), 1 + rnd.Next(12), 1).AddDays(-1);
            for (int i = 0; i < n; i++)
            {
                double u = rnd.NextDouble();
                if (u < 0.4) date = date.AddMonths(1);
                else if (u < 0.9) date = date.AddDays(1 + rnd.Next(20));
                ts.Add(new SeriesOrdinate<DateTime, double>(date, values[i]));
            }
            for (int s = 0; s < n / 50; s++)
            {
                int a = rnd.Next(n), b = rnd.Next(n);
                var tmp = ts[a];
                ts[a] = ts[b];
                ts[b] = tmp;
            }
            return ts;
        }

        /// <summary>
        /// Asserts that two block series computations give the same series, or throw the same exception type.
        /// </summary>
        /// <param name="expected">The 2.2.0 computation.</param>
        /// <param name="actual">The computation under test.</param>
        /// <param name="startMonth">For the custom year with a start month, the start month: the reference date is then the
        /// tested date shifted forward and back by the month shift, as the reference clamps it. Zero for exact dates.</param>
        private static void AssertSame(Func<TimeSeries> expected, Func<TimeSeries> actual, int startMonth = 0)
        {
            TimeSeries e = null, a = null;
            Type eError = null, aError = null;
            try { e = expected(); } catch (Exception ex) { eError = ex.GetType(); }
            try { a = actual(); } catch (Exception ex) { aError = ex.GetType(); }
            Assert.AreEqual(eError, aError);
            if (e != null && a != null) AssertSame(e, a, startMonth);
        }

        /// <summary>
        /// Asserts that two block series have the same dates, kinds and value bits, in the same order.
        /// </summary>
        /// <param name="expected">The 2.2.0 result.</param>
        /// <param name="actual">The result under test.</param>
        /// <param name="startMonth">See <see cref="AssertSame(Func{TimeSeries}, Func{TimeSeries}, int)"/>.</param>
        private static void AssertSame(TimeSeries expected, TimeSeries actual, int startMonth = 0)
        {
            int shift = startMonth > 1 ? 13 - startMonth : 0;
            Assert.AreEqual(expected.TimeInterval, actual.TimeInterval);
            Assert.AreEqual(expected.Count, actual.Count);
            for (int i = 0; i < expected.Count; i++)
            {
                var date = shift != 0 ? actual[i].Index.AddMonths(shift).AddMonths(-shift) : actual[i].Index;
                Assert.AreEqual(expected[i].Index.Ticks, date.Ticks);
                Assert.AreEqual(expected[i].Index.Kind, date.Kind);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i].Value), BitConverter.DoubleToInt64Bits(actual[i].Value));
            }
        }

        #region 2.2.0 implementations

        private static TimeSeries LegacyClone(TimeSeries ts) => new TimeSeries(ts.ToXElement());

        private static TimeSeries LegacySmoothed(TimeSeries ts, SmoothingFunctionType smoothingFunction, int period)
        {
            TimeSeries smoothedSeries = LegacyClone(ts);
            if (smoothingFunction == SmoothingFunctionType.MovingAverage)
                smoothedSeries = period == 1 ? LegacyClone(ts) : ts.MovingAverage(period);
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
                smoothedSeries = period == 1 ? LegacyClone(ts) : ts.MovingSum(period);
            else if (smoothingFunction == SmoothingFunctionType.Difference)
                smoothedSeries = ts.Difference(period);
            return smoothedSeries;
        }

        private static SeriesOrdinate<DateTime, double> LegacyBlock(List<SeriesOrdinate<DateTime, double>> blockData, BlockFunctionType blockFunction)
        {
            var ordinate = new SeriesOrdinate<DateTime, double>() { Value = double.NaN };
            if (blockFunction == BlockFunctionType.Minimum)
            {
                double min = double.MaxValue;
                for (int j = 0; j < blockData.Count; j++)
                {
                    if (blockData[j].Value < min)
                    {
                        min = blockData[j].Value;
                        ordinate.Index = blockData[j].Index;
                        ordinate.Value = blockData[j].Value;
                    }
                }
            }
            else if (blockFunction == BlockFunctionType.Maximum)
            {
                double max = double.MinValue;
                for (int j = 0; j < blockData.Count; j++)
                {
                    if (blockData[j].Value > max)
                    {
                        max = blockData[j].Value;
                        ordinate.Index = blockData[j].Index;
                        ordinate.Value = blockData[j].Value;
                    }
                }
            }
            else if (blockFunction == BlockFunctionType.Sum)
            {
                double sum = 0;
                for (int j = 0; j < blockData.Count; j++)
                    sum += blockData[j].Value;
                if (blockData.Count > 0)
                {
                    ordinate.Index = blockData.Last().Index;
                    ordinate.Value = sum;
                }
            }
            else if (blockFunction == BlockFunctionType.Average)
            {
                double sum = 0;
                for (int j = 0; j < blockData.Count; j++)
                    sum += blockData[j].Value;
                if (blockData.Count > 0)
                {
                    ordinate.Index = blockData.Last().Index;
                    ordinate.Value = sum / blockData.Count;
                }
            }
            return ordinate;
        }

        private static TimeSeries LegacyCalendarYearSeries(TimeSeries ts, BlockFunctionType blockFunction, SmoothingFunctionType smoothingFunction, int period)
        {
            var result = new TimeSeries(TimeInterval.Irregular);
            var smoothedSeries = LegacySmoothed(ts, smoothingFunction, period);
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {
                var ordinate = LegacyBlock(smoothedSeries.Where(x => x.Index.Year == i).ToList(), blockFunction);
                if (!double.IsNaN(ordinate.Value)) result.Add(ordinate);
            }
            return result;
        }

        private static TimeSeries LegacyCustomYearSeries(TimeSeries ts, int startMonth, BlockFunctionType blockFunction, SmoothingFunctionType smoothingFunction, int period)
        {
            var result = new TimeSeries(TimeInterval.Irregular);
            var smoothedSeries = LegacySmoothed(ts, smoothingFunction, period);
            int shift = startMonth != 1 ? 12 - startMonth + 1 : 0;
            smoothedSeries = startMonth != 1 ? smoothedSeries.ShiftDatesByMonth(shift) : smoothedSeries;
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {
                var ordinate = LegacyBlock(smoothedSeries.Where(x => x.Index.Year == i).ToList(), blockFunction);
                if (!double.IsNaN(ordinate.Value)) result.Add(ordinate);
            }
            result = startMonth != 1 ? result.ShiftDatesByMonth(-shift) : result;
            return result;
        }

        private static TimeSeries LegacyCustomYearSeries(TimeSeries ts, int startMonth, int endMonth, BlockFunctionType blockFunction, SmoothingFunctionType smoothingFunction, int period)
        {
            var result = new TimeSeries(TimeInterval.Irregular);
            var smoothedSeries = LegacySmoothed(ts, smoothingFunction, period);
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {
                var blockData = new List<SeriesOrdinate<DateTime, double>>();
                if (startMonth <= endMonth)
                {
                    for (int j = startMonth; j <= endMonth; j++)
                        blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == j).ToList());
                }
                else
                {
                    for (int j = startMonth; j <= 12; j++)
                        blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i - 1 && x.Index.Month == j).ToList());
                    for (int j = 1; j <= endMonth; j++)
                        blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == j).ToList());
                }
                var ordinate = LegacyBlock(blockData, blockFunction);
                if (!double.IsNaN(ordinate.Value)) result.Add(ordinate);
            }
            return result;
        }

        private static TimeSeries LegacyMonthlySeries(TimeSeries ts, BlockFunctionType blockFunction, SmoothingFunctionType smoothingFunction, int period)
        {
            var result = new TimeSeries(TimeInterval.Irregular);
            var smoothedSeries = LegacySmoothed(ts, smoothingFunction, period);
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {
                for (int k = 1; k <= 12; k++)
                {
                    var ordinate = LegacyBlock(smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == k).ToList(), blockFunction);
                    if (!double.IsNaN(ordinate.Value)) result.Add(ordinate);
                }
            }
            return result;
        }

        private static TimeSeries LegacyQuarterlySeries(TimeSeries ts, BlockFunctionType blockFunction, SmoothingFunctionType smoothingFunction, int period)
        {
            var qStart = new int[] { 1, 4, 7, 10 };
            var qEnd = new int[] { 3, 6, 9, 12 };
            var result = new TimeSeries(TimeInterval.Irregular);
            var smoothedSeries = LegacySmoothed(ts, smoothingFunction, period);
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {
                for (int q = 0; q < qEnd.Length; q++)
                {
                    var blockData = new List<SeriesOrdinate<DateTime, double>>();
                    for (int j = qStart[q]; j <= qEnd[q]; j++)
                        blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == j).ToList());
                    var ordinate = LegacyBlock(blockData, blockFunction);
                    if (!double.IsNaN(ordinate.Value)) result.Add(ordinate);
                }
            }
            return result;
        }

        #endregion
    }
}
