using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using System;
using System.Xml.Linq;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Threading.Tasks;
using System.Text.RegularExpressions;
using System.IO;
using System.Net;

namespace Data.TimeSeriesAnalysis
{
    /// <summary>
    /// Unit tests for the TimeSeries class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_TimeSeries
    {
        /// <summary>
        /// Test the construction methods for the TimeSeries class
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var ts1 = new TimeSeries();

            var ts2 = new TimeSeries(TimeInterval.OneYear);

            var ts3 = new TimeSeries(TimeInterval.OneYear, new DateTime(2023, 01, 01), new DateTime(2023, 12, 31));

            var ts4 = new TimeSeries(TimeInterval.OneYear, new DateTime(2023, 01, 01), new DateTime(2023, 12, 31), 5);

            var ts5 = new TimeSeries(TimeInterval.OneYear, new DateTime(2023, 01, 01), new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 });

            var xElement = new XElement("xElement",
                new XAttribute("TimeInterval", TimeInterval.OneYear),
                new XElement("SeriesOrdinate", new XAttribute("Index", new DateTime(2023, 01, 01).ToString("o", CultureInfo.InvariantCulture)), new XAttribute("Value", 5.0.ToString("G17", CultureInfo.InvariantCulture))),
                new XElement("SeriesOrdinate", new XAttribute("Index", new DateTime(2023, 12, 31).ToString("o", CultureInfo.InvariantCulture)), new XAttribute("Value", 5.0.ToString("G17", CultureInfo.InvariantCulture))),
                new XElement("SeriesOrdinate", new XAttribute("Index", new DateTime(2024, 01, 01).ToString("o", CultureInfo.InvariantCulture)), new XAttribute("Value", 5.0.ToString("G17", CultureInfo.InvariantCulture))));

            var ts6 = new TimeSeries(xElement);
        }

        /// <summary>
        /// Test the getter for the TimeSeries class
        /// </summary>
        [TestMethod]
        public void Test_Getters()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 });

            var ti = ts.TimeInterval;
            Assert.AreEqual(TimeInterval.OneMonth, ti);

            var missing = ts.HasMissingValues;
            Assert.IsFalse(missing);

            var start = ts.StartDate;
            Assert.AreEqual(new DateTime(2023, 01, 01), start);

            var end = ts.EndDate;
            Assert.AreEqual(new DateTime(2023, 12, 01), end);
        }

        /// <summary>
        /// Test the method that converts a TimeSeries object to an XElement object
        /// </summary>
        [TestMethod]
        public void Test_ToXElement()
        {
            var ts = new TimeSeries(TimeInterval.OneYear, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });
            var xElement = ts.ToXElement();

            var name = xElement.Name;

            // will be overwritten, but must initialize to something;
            TimeInterval time = TimeInterval.OneDay;
            List<SeriesOrdinate<DateTime, double>> ordinates = new List<SeriesOrdinate<DateTime, double>>();

            if (xElement.Attribute(nameof(TimeInterval)) != null)
                Enum.TryParse(xElement.Attribute(nameof(TimeInterval)).Value, out time);

            foreach (XElement ordinate in xElement.Elements("SeriesOrdinate"))
            {
                var index = DateTime.ParseExact(ordinate.Attribute("Index").Value, "o", CultureInfo.InvariantCulture, DateTimeStyles.RoundtripKind);
                double.TryParse(ordinate.Attribute("Value").Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var value);
                ordinates.Add(new SeriesOrdinate<DateTime, double>(index, value));
            }

            for (int i = 0; i < ts.Count; i++)
            {
                Assert.AreEqual(ts[i].Index, ordinates[i].Index);
                Assert.AreEqual(ts[i].Value, ordinates[i].Value);
            }

            Assert.AreEqual("TimeSeries", name);
            Assert.AreEqual(ts.TimeInterval, time);
        }

        /// <summary>
        /// Test the clone method
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });
            var newTS = ts.Clone();

            for (int i = 0; i < ts.Count; i++)
            {
                Assert.AreEqual(ts[i].Index, newTS[i].Index);
                Assert.AreEqual(ts[i].Value, newTS[i].Value);
            }

            Assert.AreEqual(ts.HasMissingValues, newTS.HasMissingValues);
            Assert.AreEqual(ts.TimeInterval, newTS.TimeInterval);
            Assert.AreEqual(ts.StartDate, newTS.StartDate);
            Assert.AreEqual(ts.EndDate, newTS.EndDate);
        }

        /// <summary>
        /// Test the sort method
        /// </summary>
        [TestMethod]
        public void Test_Sort()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);
            ts.SortByTime();

            for (int i = 0; i < ts.Count; i++)
            {
                // values are in original order b/c that's how they were added and time stamped
                Assert.AreEqual(values[i], ts[i].Value);
            }

            ts.SortByValue();
            // now values must be sorted
            Array.Sort(values);

            for (int i = 0; i < ts.Count; i++)
            {
                Assert.AreEqual(values[i], ts[i].Value);
            }
        }

        /// <summary>
        /// Helper function for Test_Math
        /// </summary>
        /// <param name="ts">The TimeSeries object being evaluated.</param>
        /// <param name="values">The validation values to compare against.</param>
        /// <param name="tolerance">The tolerance level for the comparison. Default: 1E-6</param>
        public void Equal(TimeSeries ts, double[] values, double tolerance = 1E-6)
        {
            for (int i = 0; i < ts.Count; i++)
            {
                Assert.AreEqual(values[i], ts[i].Value, tolerance);
            }
        }

        /// <summary>
        /// Test the methods that manipulate the TimeSeries values 
        /// </summary>
        [TestMethod]
        public void Test_Math()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            // Test Add method
            for (int i = 0; i < values.Length; i++) { values[i] += 10; }
            ts.Add(10);
            Equal(ts, values);

            // Test Subtract method
            for (int i = 0; i < values.Length; i++) { values[i] -= 5; }
            ts.Subtract(5);
            Equal(ts, values);

            // Test Multiply method
            for (int i = 0; i < values.Length; i++) { values[i] *= 4; }
            ts.Multiply(4);
            Equal(ts, values);

            // Test Divide method
            for (int i = 0; i < values.Length; i++) { values[i] /= -2; }
            ts.Divide(-2);
            Equal(ts, values);

            // Test Equal method
            for (int i = 0; i < values.Length; i++) { values[i] = Math.Abs(values[i]); }
            ts.AbsoluteValue();
            Equal(ts, values);

            // Test Exponentiate method
            for (int i = 0; i < values.Length; i++) { values[i] = Math.Pow(values[i], 2); }
            ts.Exponentiate(2);
            Equal(ts, values);

            // Test LogTransform method
            for (int i = 0; i < values.Length; i++) { values[i] = Math.Log10(values[i]); }
            ts.LogTransform();
            Equal(ts, values);

            // Test Inverse method
            for (int i = 0; i < values.Length; i++) { values[i] = 1 / values[i]; }
            ts.Inverse();
            Equal(ts, values);

            // Test Standardize method
            var mean = Numerics.Data.Statistics.Statistics.Mean(values);
            var sd = Numerics.Data.Statistics.Statistics.StandardDeviation(values);
            for (int i = 0; i < values.Length; i++) { values[i] = (values[i] - mean) / sd; }
            ts.Standardize();
            Equal(ts, values);
        }

        /// <summary>
        /// The indexed divide overload rejects a zero divisor with the same exception contract as
        /// its all-values twin and leaves the selected values unchanged.
        /// </summary>
        [TestMethod]
        public void Test_Divide_Indexed_RejectsZeroLikeTwin()
        {
            var plain = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 2d, 4d });
            var indexed = plain.Clone();

            var expected = Assert.Throws<ArgumentException>(() => plain.Divide(0d));
            var actual = Assert.Throws<ArgumentException>(() => indexed.Divide(0d, new[] { 1 }));

            Assert.AreEqual(expected.ParamName, actual.ParamName);
            Assert.AreEqual(expected.Message, actual.Message);
            Assert.AreEqual(2d, indexed[0].Value, 0d);
            Assert.AreEqual(4d, indexed[1].Value, 0d);
        }

        /// <summary>
        /// Verifies that the indexed log transform skips out-of-range indexes like its sibling overloads.
        /// </summary>
        [TestMethod]
        public void Test_LogTransform_Indexed_SkipsOutOfRangeIndexes()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 10, 100, 1000 });
            ts.LogTransform(new[] { 1, 7, -1 });
            Assert.AreEqual(10d, ts[0].Value, 1E-6);
            Assert.AreEqual(2d, ts[1].Value, 1E-6);
            Assert.AreEqual(1000d, ts[2].Value, 1E-6);
        }

        /// <summary>
        /// Verifies that the indexed log transform still marks in-range non-positive values as missing.
        /// </summary>
        [TestMethod]
        public void Test_LogTransform_Indexed_MarksNonPositiveInRangeValuesMissing()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { -5, 100, 1000 });
            ts.LogTransform(new[] { 0, 1 });
            Assert.IsTrue(double.IsNaN(ts[0].Value));
            Assert.AreEqual(2d, ts[1].Value, 1E-6);
        }

        /// <summary>
        /// Verifies that the indexed inverse skips out-of-range indexes like its sibling overloads.
        /// </summary>
        [TestMethod]
        public void Test_Inverse_Indexed_SkipsOutOfRangeIndexes()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 2, 4, 8 });
            ts.Inverse(new[] { 0, 5, -2 });
            Assert.AreEqual(0.5, ts[0].Value, 1E-6);
            Assert.AreEqual(4d, ts[1].Value, 1E-6);
            Assert.AreEqual(8d, ts[2].Value, 1E-6);
        }

        /// <summary>
        /// Verifies that the indexed inverse still marks in-range zero and missing values as missing.
        /// </summary>
        [TestMethod]
        public void Test_Inverse_Indexed_MarksZeroInRangeValuesMissing()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 0, 4, 8 });
            ts.Inverse(new[] { 0, 2 });
            Assert.IsTrue(double.IsNaN(ts[0].Value));
            Assert.AreEqual(0.125, ts[2].Value, 1E-6);
        }

        /// <summary>
        /// Test the CumulativeSum method
        /// </summary>
        [TestMethod]
        public void Test_Cumulative()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);
            var newValues = new double[12];

            for (int i = 1; i < values.Length; i++) { values[i] = values[i - 1] + values[i]; }
            var newTS = ts.CumulativeSum();
            Equal(newTS, values);

        }

        /// <summary>
        /// Verifies that the cumulative sum preserves the source series' time interval.
        /// </summary>
        [TestMethod]
        public void Test_Cumulative_PreservesTimeInterval()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });
            var newTS = ts.CumulativeSum();
            Assert.AreEqual(ts.TimeInterval, newTS.TimeInterval);
        }

        /// <summary>
        /// Test the successive Difference method
        /// </summary>
        [TestMethod]
        public void Test_Difference()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);
            var newValues = new double[12];

            for (int i = 1; i < newValues.Length; i++) { newValues[i] = values[i] - values[i - 1]; }
            var newTS = ts.Difference();
            // don't want first value
            Equal(newTS, newValues.Skip(1).ToArray());
            for (int i = 0; i < newTS.Count; i++)
                Assert.AreEqual(ts[i].Index, newTS[i].Index);
        }

        /// <summary>
        /// Test the methods that handle missing values
        /// </summary>
        [TestMethod]
        public void Test_Missing()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, double.NaN, double.NaN };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            var missing = ts.NumberOfMissingValues();
            Assert.AreEqual(2, missing);

            ts.ReplaceMissingData(-1);
            missing = ts.NumberOfMissingValues();
            Assert.AreEqual(0, missing);
            Assert.AreEqual(-1, ts[10].Value);
            Assert.AreEqual(-1, ts[11].Value);

            values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, double.NaN, double.NaN };
            ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            ts.InterpolateMissingData(2);
            missing = ts.NumberOfMissingValues();
            Assert.AreEqual(0, missing);
            Assert.AreEqual(11.9, ts[10].Value, 1E-6);
            Assert.AreEqual(8.9, ts[11].Value, 1E-6);
        }

        /// <summary>
        /// Verifies that the indexed interpolation matches its non-indexed twin at the series start instead of reading before the first ordinate.
        /// </summary>
        [TestMethod]
        public void Test_InterpolateMissingData_Indexed_MatchesTwinAtSeriesStart()
        {
            var indexed = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 1, double.NaN, double.NaN });
            var plain = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 1, double.NaN, double.NaN });
            indexed.InterpolateMissingData(1, new[] { 1 });
            plain.InterpolateMissingData(1);
            for (int i = 0; i < plain.Count; i++)
            {
                Assert.AreEqual(plain[i].Value, indexed[i].Value);
            }
        }

        /// <summary>
        /// Test the method that adds different interval times
        /// </summary>
        [TestMethod]
        public void Test_AddInterval()
        {
            // DateTime(year, month, day, hour, minute, second)
            DateTime original = new DateTime(2023, 01, 01, 09, 30, 25);

            var test = TimeSeries.AddTimeInterval(original, TimeInterval.OneMinute);
            Assert.AreEqual(original.AddMinutes(1), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.FiveMinute);
            Assert.AreEqual(original.AddMinutes(5), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.FifteenMinute);
            Assert.AreEqual(original.AddMinutes(15), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.ThirtyMinute);
            Assert.AreEqual(original.AddMinutes(30), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.OneHour);
            Assert.AreEqual(original.AddHours(1), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.SixHour);
            Assert.AreEqual(original.AddHours(6), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.TwelveHour);
            Assert.AreEqual(original.AddHours(12), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.OneDay);
            Assert.AreEqual(original.AddDays(1), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.SevenDay);
            Assert.AreEqual(original.AddDays(7), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.OneMonth);
            Assert.AreEqual(original.AddMonths(1), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.OneQuarter);
            Assert.AreEqual(original.AddMonths(3), test);

            test = TimeSeries.AddTimeInterval(original, TimeInterval.OneYear);
            Assert.AreEqual(original.AddYears(1), test);
        }

        /// <summary>
        /// Test the method that subtracts different interval times
        /// </summary>
        [TestMethod]
        public void Test_SubtractInterval()
        {
            // DateTime(year, month, day, hour, minute, second)
            DateTime original = new DateTime(2023, 01, 01, 09, 30, 25);

            var test = TimeSeries.SubtractTimeInterval(original, TimeInterval.OneMinute);
            Assert.AreEqual(original.AddMinutes(-1), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.FiveMinute);
            Assert.AreEqual(original.AddMinutes(-5), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.FifteenMinute);
            Assert.AreEqual(original.AddMinutes(-15), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.ThirtyMinute);
            Assert.AreEqual(original.AddMinutes(-30), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.OneHour);
            Assert.AreEqual(original.AddHours(-1), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.SixHour);
            Assert.AreEqual(original.AddHours(-6), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.TwelveHour);
            Assert.AreEqual(original.AddHours(-12), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.OneDay);
            Assert.AreEqual(original.AddDays(-1), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.SevenDay);
            Assert.AreEqual(original.AddDays(-7), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.OneMonth);
            Assert.AreEqual(original.AddMonths(-1), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.OneQuarter);
            Assert.AreEqual(original.AddMonths(-3), test);

            test = TimeSeries.SubtractTimeInterval(original, TimeInterval.OneYear);
            Assert.AreEqual(original.AddYears(-1), test);
        }

        /// <summary>
        /// Test the method that converts the time interval to hours
        /// </summary>
        [TestMethod]
        public void Test_TimeIntervalInHours()
        {
            var valid = new double[] { 1d / 60d, 1 / 12, 3 / 12, 1 / 2, 1, 6, 12, 24, 168 };

            double test = TimeSeries.TimeIntervalInHours(TimeInterval.OneMinute);
            Assert.AreEqual(1d / 60d, test);

            test = TimeSeries.TimeIntervalInHours(TimeInterval.FiveMinute);
            Assert.AreEqual(1d / 12d, test);

            test = TimeSeries.TimeIntervalInHours(TimeInterval.FifteenMinute);
            Assert.AreEqual(3d / 12d, test);

            test = TimeSeries.TimeIntervalInHours(TimeInterval.ThirtyMinute);
            Assert.AreEqual(1d / 2d, test);

            test = TimeSeries.TimeIntervalInHours(TimeInterval.OneHour);
            Assert.AreEqual(1, test);

            test = TimeSeries.TimeIntervalInHours(TimeInterval.SixHour);
            Assert.AreEqual(6, test);

            test = TimeSeries.TimeIntervalInHours(TimeInterval.TwelveHour);
            Assert.AreEqual(12, test);

            test = TimeSeries.TimeIntervalInHours(TimeInterval.OneDay);
            Assert.AreEqual(24, test);

            test = TimeSeries.TimeIntervalInHours(TimeInterval.SevenDay);
            Assert.AreEqual(168, test);
        }

        /// <summary>
        /// Test the moving average method. Validation values were attained from MS Excel.
        /// </summary>
        [TestMethod]
        public void Test_MovingAverage()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });

            var newTS = ts.MovingAverage(5);
            var valid = new double[] { 17.6, 20.4, 26.8, 22.2, 24.8, 25.4, 22.6, 15.6 };
            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(valid[i], newTS[i].Value);
            }
        }

        /// <summary>
        /// Test the moving sum method. Validation values were attained from MS Excel.
        /// </summary>
        [TestMethod]
        public void Test_MovingSum()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });

            var newTS = ts.MovingSum(5);
            var valid = new double[] { 88, 102, 134, 111, 124, 127, 113, 78 };
            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(valid[i], newTS[i].Value);
            }
        }

        /// <summary>
        /// Test that MovingSum propagates NaN strictly by default (any NaN in window -> NaN out).
        /// Matches pandas Series.rolling(w).sum() default min_periods=w behavior.
        /// </summary>
        [TestMethod]
        public void Test_MovingSum_NaN_Strict()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 1.0, double.NaN, 3.0, 4.0, 5.0 });

            var newTS = ts.MovingSum(2);
            // window positions: [1,NaN], [NaN,3], [3,4], [4,5] -> [NaN, NaN, 7, 9]
            var valid = new double[] { double.NaN, double.NaN, 7.0, 9.0 };
            Assert.AreEqual(valid.Length, newTS.Count);
            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(valid[i], newTS[i].Value);
            }
        }

        /// <summary>
        /// Test that MovingSum with minValidCount=1 skips NaN values (pandas min_periods=1 semantics).
        /// Sum is over observed values only; no scaling is applied.
        /// </summary>
        [TestMethod]
        public void Test_MovingSum_NaN_SkipMode()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 1.0, double.NaN, 3.0, 4.0, 5.0 });

            var newTS = ts.MovingSum(2, minValidCount: 1);
            // window positions: [1,NaN]->1, [NaN,3]->3, [3,4]->7, [4,5]->9
            var valid = new double[] { 1.0, 3.0, 7.0, 9.0 };
            Assert.AreEqual(valid.Length, newTS.Count);
            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(valid[i], newTS[i].Value);
            }
        }

        /// <summary>
        /// Test that MovingAverage propagates NaN strictly by default.
        /// </summary>
        [TestMethod]
        public void Test_MovingAverage_NaN_Strict()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 1.0, double.NaN, 3.0, 4.0, 5.0 });

            var newTS = ts.MovingAverage(2);
            var valid = new double[] { double.NaN, double.NaN, 3.5, 4.5 };
            Assert.AreEqual(valid.Length, newTS.Count);
            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(valid[i], newTS[i].Value);
            }
        }

        /// <summary>
        /// Test that MovingAverage with minValidCount=1 averages over non-NaN values only.
        /// Denominator is the count of observed values, not the window size, so the average is not biased low.
        /// </summary>
        [TestMethod]
        public void Test_MovingAverage_NaN_SkipMode()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 1.0, double.NaN, 3.0, 4.0, 5.0 });

            var newTS = ts.MovingAverage(2, minValidCount: 1);
            // [1,NaN]->1/1=1, [NaN,3]->3/1=3, [3,4]->7/2=3.5, [4,5]->9/2=4.5
            var valid = new double[] { 1.0, 3.0, 3.5, 4.5 };
            Assert.AreEqual(valid.Length, newTS.Count);
            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(valid[i], newTS[i].Value);
            }
        }

        /// <summary>
        /// Regression guard: Difference relies on IEEE-754 subtraction to propagate NaN automatically.
        /// </summary>
        [TestMethod]
        public void Test_Difference_NaN()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), new double[] { 1.0, double.NaN, 3.0, 7.0 });

            var newTS = ts.Difference();
            // diffs: NaN-1=NaN, 3-NaN=NaN, 7-3=4
            Assert.AreEqual(3, newTS.Count);
            Assert.AreEqual(double.NaN, newTS[0].Value);
            Assert.AreEqual(double.NaN, newTS[1].Value);
            Assert.AreEqual(4.0, newTS[2].Value);
        }

        /// <summary>
        /// Test that CalendarYearSeries with Sum/Average drops years containing any NaN
        /// (rather than silently treating NaN as zero).
        /// </summary>
        [TestMethod]
        public void Test_CalendarYearSeries_NaN()
        {
            // Two years of monthly data; year 2 has one NaN month.
            var values = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                        1, double.NaN, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            // Sum: year 1 = 78, year 2 = NaN (dropped from output)
            var sumTS = ts.CalendarYearSeries(BlockFunctionType.Sum);
            Assert.AreEqual(1, sumTS.Count);
            Assert.AreEqual(78.0, sumTS[0].Value);
            Assert.AreEqual(2023, sumTS[0].Index.Year);

            // Average: same; year 2 dropped because NaN propagates
            var avgTS = ts.CalendarYearSeries(BlockFunctionType.Average);
            Assert.AreEqual(1, avgTS.Count);
            Assert.AreEqual(78.0 / 12.0, avgTS[0].Value, 1E-9);

            // Maximum: NaN values are skipped via the < / > comparisons (false on NaN), so both years survive.
            var maxTS = ts.CalendarYearSeries(BlockFunctionType.Maximum);
            Assert.AreEqual(2, maxTS.Count);
            Assert.AreEqual(12.0, maxTS[0].Value);
            Assert.AreEqual(12.0, maxTS[1].Value);
        }

        /// <summary>
        /// Regression guard: MonthlySeries already propagated NaN naively before this fix; verify it still does.
        /// </summary>
        [TestMethod]
        public void Test_MonthlySeries_NaN()
        {
            // Two months of weekly data; month 2 has one NaN.
            var values = new double[] { 10, 20, 30, 40, 50, 60, double.NaN, 80 };
            var ts = new TimeSeries(TimeInterval.SevenDay, new DateTime(2023, 01, 01), values);

            var sumTS = ts.MonthlySeries(BlockFunctionType.Sum);
            // Jan 2023 has 5 weekly values: 10+20+30+40+50 = 150
            // Feb 2023 has 3 weekly values, one NaN -> sum = NaN, dropped
            // Mar-Dec 2023 are empty -> skipped (no spurious zeros)
            Assert.AreEqual(1, sumTS.Count);
            Assert.AreEqual(150.0, sumTS[0].Value, 1E-9);
        }

        /// <summary>
        /// Test that PeaksOverThresholdSeries with a multi-day MovingSum smoothing on data with NaN
        /// does not produce spurious exceedances from windows that touched a missing day.
        /// This is the user-reported bug: 2-day moving sum on precipitation with gaps was returning
        /// the value of the single non-missing day instead of NaN, causing zero-baseline exceedances.
        /// </summary>
        [TestMethod]
        public void Test_PeaksOverThreshold_MovingSum_NaN()
        {
            // One legitimate 2-day exceedance (5 + 6 = 11). The trailing 8.0 sits at the very end
            // of the series with a NaN before it, so the only window touching it is [NaN, 8.0]:
            // a NaN window carries no event, so nothing above threshold 7 may come from it and
            // only [5, 6] = 11 remains.
            var values = new double[] { 0.1, 0.2, 5.0, 6.0, 0.1, 0.0, 0.1, 0.2, 0.0, double.NaN, 8.0 };
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 01, 01), values);

            var pot = ts.PeaksOverThresholdSeries(threshold: 7.0, minStepsBetweenEvents: 1,
                                                  smoothingFunction: SmoothingFunctionType.MovingSum, period: 2);

            Assert.AreEqual(1, pot.Count);
            Assert.AreEqual(11.0, pot[0].Value, 1E-9);
        }

        /// <summary>
        /// Test the method that shifts all dates to a new starting date
        /// </summary>
        [TestMethod]
        public void Test_ShiftAllDates()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });
            var newTS = ts.Clone();
            newTS.ShiftAllDates(new DateTime(2024, 02, 04));

            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(ts[i].Index.AddYears(1).AddMonths(1).AddDays(3), newTS[i].Index);
                Assert.AreEqual(ts[i].Value, newTS[i].Value);
            }
        }

        /// <summary>
        /// Test the method that shifts all dates by a specified number of days.
        /// </summary>
        [TestMethod]
        public void Test_ShiftDatesByDays()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });
            var newTS = ts.ShiftDatesByDay(5);

            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(ts[i].Index.AddDays(5), newTS[i].Index);
                Assert.AreEqual(ts[i].Value, newTS[i].Value);
            }
        }

        /// <summary>
        /// Test the method that shifts all dates by a specified number of months
        /// </summary>
        [TestMethod]
        public void Test_ShiftDatesByMonth()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });
            var newTS = ts.ShiftDatesByMonth(5);

            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(ts[i].Index.AddMonths(5), newTS[i].Index);
                Assert.AreEqual(ts[i].Value, newTS[i].Value);
            }
        }

        /// <summary>
        /// Test the method that shifts all dates by a specified number of years
        /// </summary>
        [TestMethod]
        public void Test_ShiftDatesByYear()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });
            var newTS = ts.ShiftDatesByYear(5);

            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(ts[i].Index.AddYears(5), newTS[i].Index);
                Assert.AreEqual(ts[i].Value, newTS[i].Value);
            }
        }

        /// <summary>
        /// Test the method that clips the TimeSeries object, from both the beginning and end.
        /// </summary>
        [TestMethod]
        public void Test_ClipTimeSeries()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 });
            var start = new DateTime(2023, 11, 01);
            var end = new DateTime(2023, 12, 01);
            var newTS = ts.ClipTimeSeries(start, end);

            Assert.AreEqual(2, newTS.Count);
            Assert.AreEqual(start, newTS.StartDate);
            Assert.AreEqual(end, newTS.EndDate);

            for (int i = 0; i < newTS.Count; i++)
            {
                Assert.AreEqual(ts[10 + i].Index, newTS[i].Index);
                Assert.AreEqual(ts[10 + i].Value, newTS[i].Value);
            }
        }

        /// <summary>
        /// Convert from 1-hr to 15-min assuming block averages. 
        /// </summary>
        [TestMethod]
        public void Test_ConvertTimeInterval_1hr_to_15min()
        {
            var data = new double[] { 43, 42, 42, 46, 51, 56, 58, 60, 67, 278, 278, 218, 256, 243, 225, 199, 262, 616, 1440, 2180, 3140, 4520, 4700, 4340, 4190, 3720, 3120, 2550, 2110, 1720, 1350, 1110, 950, 830, 755, 696, 664, 612, 572, 520, 458, 410, 362, 315, 273, 270, 277, 246, 247, 248, 261, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162 };
            var trueValues = new double[] { 43, 42.75, 42.5, 42.25, 42, 42, 42, 42, 42, 43, 44, 45, 46, 47.25, 48.5, 49.75, 51, 52.25, 53.5, 54.75, 56, 56.5, 57, 57.5, 58, 58.5, 59, 59.5, 60, 61.75, 63.5, 65.25, 67, 119.75, 172.5, 225.25, 278, 278, 278, 278, 278, 263, 248, 233, 218, 227.5, 237, 246.5, 256, 252.75, 249.5, 246.25, 243, 238.5, 234, 229.5, 225, 218.5, 212, 205.5, 199, 214.75, 230.5, 246.25, 262, 350.5, 439, 527.5, 616, 822, 1028, 1234, 1440, 1625, 1810, 1995, 2180, 2420, 2660, 2900, 3140, 3485, 3830, 4175, 4520, 4565, 4610, 4655, 4700, 4610, 4520, 4430, 4340, 4302.5, 4265, 4227.5, 4190, 4072.5, 3955, 3837.5, 3720, 3570, 3420, 3270, 3120, 2977.5, 2835, 2692.5, 2550, 2440, 2330, 2220, 2110, 2012.5, 1915, 1817.5, 1720, 1627.5, 1535, 1442.5, 1350, 1290, 1230, 1170, 1110, 1070, 1030, 990, 950, 920, 890, 860, 830, 811.25, 792.5, 773.75, 755, 740.25, 725.5, 710.75, 696, 688, 680, 672, 664, 651, 638, 625, 612, 602, 592, 582, 572, 559, 546, 533, 520, 504.5, 489, 473.5, 458, 446, 434, 422, 410, 398, 386, 374, 362, 350.25, 338.5, 326.75, 315, 304.5, 294, 283.5, 273, 272.25, 271.5, 270.75, 270, 271.75, 273.5, 275.25, 277, 269.25, 261.5, 253.75, 246, 246.25, 246.5, 246.75, 247, 247.25, 247.5, 247.75, 248, 251.25, 254.5, 257.75, 261, 254, 247, 240, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 215.25, 197.5, 179.75, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162 };

            var ts = new TimeSeries(TimeInterval.OneHour, new DateTime(1973, 5, 1), data);
            var newTS = ts.ConvertTimeInterval(TimeInterval.FifteenMinute);

            for (int i = 0; i < trueValues.Length; i++)
                Assert.AreEqual(trueValues[i], newTS[i].Value, 1E-2);
        }

        /// <summary>
        /// Convert from 1-hr to 6-hr block averages. 
        /// </summary>
        [TestMethod]
        public void Test_ConvertTimeInterval_1hr_to_6hr()
        {
            var data = new double[] { 43, 42, 42, 46, 51, 56, 58, 60, 67, 278, 278, 218, 256, 243, 225, 199, 262, 616, 1440, 2180, 3140, 4520, 4700, 4340, 4190, 3720, 3120, 2550, 2110, 1720, 1350, 1110, 950, 830, 755, 696, 664, 612, 572, 520, 458, 410, 362, 315, 273, 270, 277, 246, 247, 248, 261, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162 };
            var trueValues = new double[] { 46.67, 159.83, 300.17, 3386.67, 2901.67, 948.5, 539.33, 290.5, 242.5, 233, 233, 233, 209.33, 162, 162, 162, 162 };

            var ts = new TimeSeries(TimeInterval.OneHour, new DateTime(1973, 5, 1), data);
            var newTS = ts.ConvertTimeInterval(TimeInterval.SixHour);

            for (int i = 0; i < trueValues.Length; i++)
                Assert.AreEqual(trueValues[i], newTS[i].Value, 1E-2);
        }

        /// <summary>
        /// Convert from 1-hr to 1-day block averages. 
        /// </summary>
        [TestMethod]
        public void Test_ConvertTimeInterval_1hr_to_1Day()
        {
            var data = new double[] { 43, 42, 42, 46, 51, 56, 58, 60, 67, 278, 278, 218, 256, 243, 225, 199, 262, 616, 1440, 2180, 3140, 4520, 4700, 4340, 4190, 3720, 3120, 2550, 2110, 1720, 1350, 1110, 950, 830, 755, 696, 664, 612, 572, 520, 458, 410, 362, 315, 273, 270, 277, 246, 247, 248, 261, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162 };
            var trueValues = new double[] { 973.33, 1170.00, 235.38, 173.83 };

            var ts = new TimeSeries(TimeInterval.OneHour, new DateTime(1973, 5, 1), data);
            var newTS = ts.ConvertTimeInterval(TimeInterval.OneDay);

            for (int i = 0; i < trueValues.Length; i++)
                Assert.AreEqual(trueValues[i], newTS[i].Value, 1E-2);
        }


        /// <summary>
        /// Convert from 1-hr to 1-day block sums, then back to 6-hr blocks. 
        /// </summary>
        [TestMethod]
        public void Test_ConvertTimeInterval_1hr_to_1Day_Sum()
        {
            var data = new double[] { 43, 42, 42, 46, 51, 56, 58, 60, 67, 278, 278, 218, 256, 243, 225, 199, 262, 616, 1440, 2180, 3140, 4520, 4700, 4340, 4190, 3720, 3120, 2550, 2110, 1720, 1350, 1110, 950, 830, 755, 696, 664, 612, 572, 520, 458, 410, 362, 315, 273, 270, 277, 246, 247, 248, 261, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 233, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162 };
            var ts = new TimeSeries(TimeInterval.OneHour, new DateTime(1973, 5, 1), data);

            // Convert to 1-Day
            var newTS = ts.ConvertTimeInterval(TimeInterval.OneDay, false);
            var trueValues = new double[] { 23360, 28080, 5649, 4172 };
            for (int i = 0; i < trueValues.Length; i++)
                Assert.AreEqual(trueValues[i], newTS[i].Value, 1E-2);

            // Convert to 6-hrs
            newTS = newTS.ConvertTimeInterval(TimeInterval.SixHour, false);
            trueValues = new double[] { 5840, 5840, 5840, 5840, 7020, 7020, 7020, 7020, 1412.25, 1412.25, 1412.25, 1412.25, 1043, 1043, 1043, 1043 };
            for (int i = 0; i < trueValues.Length; i++)
                Assert.AreEqual(trueValues[i], newTS[i].Value, 1E-2);
        }

        /// <summary>
        /// The standard deviation of a series with missing values must equal the sample standard
        /// deviation of its observed values.
        /// </summary>
        /// <remarks>
        /// Reference values verified against Python pandas 2.3.3 (Series.std(), which skips NaN)
        /// and exact rational arithmetic. The observed values of { 1, NaN, 2, 3 } have sample
        /// standard deviation exactly 1, and those of { NaN, 3.1, NaN, 4.7, 2.2, NaN, 5.9, 1.4 }
        /// have 1.8338484124921557, with pandas within two units in the last place of the exact
        /// rational value.
        /// </remarks>
        [TestMethod]
        public void Test_StandardDeviation_WithMissingValues()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 01, 01), new double[] { 1, double.NaN, 2, 3 });
            Assert.AreEqual(1.0, ts.StandardDeviation(), 1E-14);

            var gappy = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 01, 01), new double[] { double.NaN, 3.1, double.NaN, 4.7, 2.2, double.NaN, 5.9, 1.4 });
            Assert.AreEqual(1.8338484124921557, gappy.StandardDeviation(), 1E-14);
        }

        /// <summary>
        /// Test the statistics methods of the TimeSeries class
        /// </summary>
        [TestMethod]
        public void Test_Stats()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            // Test MinValue
            double min = ts.MinValue();
            Assert.AreEqual(5, min);

            // Test MaxValue
            double max = ts.MaxValue();
            Assert.AreEqual(48, max);

            // Test MeanValue
            double mean = ts.MeanValue();
            Assert.AreEqual(20.83333333333, mean, 1E-10);

            // Test StandardDeviation
            double sd = ts.StandardDeviation();
            Assert.AreEqual(12.40112, sd, 1E-3);

            // Test Duration
            var duration = ts.Duration();
            var dValues = new double[] { 48, 36, 33, 22, 22, 18, 16, 15, 13, 12, 10, 5 };
            var dValid = new double[] { 7.69230769230769, 15.3846153846154, 23.0769230769231, 30.7692307692308, 38.4615384615385, 46.1538461538462, 53.8461538461538, 61.5384615384615, 69.2307692307692, 76.9230769230769, 84.6153846153846, 92.3076923076923 };

            for(int i = 0; i < dValid.Length; i++)
            {
                Assert.AreEqual(dValid[i], duration[i, 0], 1E-10);
                Assert.AreEqual(dValues[i], duration[i, 1], 1E-10);
            }
        }

        /// <summary>
        /// Test the methods that return the summary statistics about a TimeSeries object. These method include SummaryPercentiles(), Percentiles(), and
        /// SummaryStatistics()
        /// </summary>
        [TestMethod]
        public void Test_SummaryStats()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            // Test SummaryPercecntiles
            var summaryP = ts.SummaryPercentiles();
            var spValid = new double[] { 7.75, 12.75, 17.00, 24.75, 41.40 };
            for (int i = 0; i < summaryP.Length; i++)
            {
                Assert.AreEqual(spValid[i], summaryP[i], 1E-10);
            }

            // Test Percentiles
            var percentiles = ts.Percentiles(new double[] { 0.10, 0.20, 0.30, 0.40, 0.60, 0.70, 0.80, 0.90 });
            var pValid = new double[] { 10.2, 12.2, 13.6, 15.4, 20.4, 22.0, 30.8, 35.7 };
            for (int i = 0; i < percentiles.Length; i++)
            {
                Assert.AreEqual(pValid[i], percentiles[i], 1E-10);
            }

            // Test SummaryStatistics
            var summaryStats = ts.SummaryStatistics();
            var testStats = new double[13];

            summaryStats.TryGetValue("Record Length", out testStats[0]);
            summaryStats.TryGetValue("Missing Values", out testStats[1]);
            summaryStats.TryGetValue("Minimum", out testStats[2]);
            summaryStats.TryGetValue("Maximum", out testStats[3]);
            summaryStats.TryGetValue("Mean", out testStats[4]);
            summaryStats.TryGetValue("Std Dev", out testStats[5]);
            summaryStats.TryGetValue("Skewness", out testStats[6]);
            summaryStats.TryGetValue("Kurtosis", out testStats[7]);
            summaryStats.TryGetValue("5%", out testStats[8]);
            summaryStats.TryGetValue("25%", out testStats[9]);
            summaryStats.TryGetValue("50%", out testStats[10]);
            summaryStats.TryGetValue("75%", out testStats[11]);
            summaryStats.TryGetValue("95%", out testStats[12]);

            var validStats = new double[] { 12, 0, 5, 48, 20.83333333333, 12.4011240937215, 1.06382220138287, 0.682181791356259, 7.75, 12.75, 17.00, 24.75, 41.40 };
            for (int i = 0; i < validStats.Length; i++)
            {
                Assert.AreEqual(validStats[i], testStats[i], 1E-10);
            }
        }

        /// <summary>
        /// Test the method that returns a summary of different hypothesis tests on the TimeSeries object. These methods were validated against
        /// their corresponding R commands. 
        /// </summary>
        /// <remarks>
        /// <list type="bullet">
        /// <item>"jarque.test()" from the "moments" package</item>
        /// <item>"Box.test()", "wilcox.test()", "t.test()", and "var.test()" from the "stats" package</item>
        /// <item>"ww.test()" and "mk.test()" from the "trend" package</item>
        /// </list>
        /// <para>
        /// <b> References: </b>
        /// <list type="bullet">
        /// <item> Lukasz Komsta (2005). moments: Moments, Cumulants, Skewness, Kurtosis and Related Tests. R package version 0.14.1, https://cran.r-project.org/web/packages/moments</item>
        /// <item> R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, 
        /// Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org/. </item>
        /// <item> Pohlert T (2023). trend: Non-Parametric Trend Tests and Change-Point Detection. R package version 1.1.6, https://CRAN.R-project.org/package=trend</item>
        /// </list>
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_SummaryHypothesisTest()
        {
            var values = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);
            var hypothesis = ts.SummaryHypothesisTest();
            var testH = new double[7];

            hypothesis.TryGetValue("Jarque-Bera test for normality", out testH[0]);
            hypothesis.TryGetValue("Ljung-Box test for independence", out testH[1]);
            hypothesis.TryGetValue("Wald-Wolfowitz test for independence and stationarity (trend)", out testH[2]);
            hypothesis.TryGetValue("Mann-Whitney test for homogeneity and stationarity (jump)", out testH[3]);
            hypothesis.TryGetValue("Mann-Kendall test for homogeneity and stationarity (trend)", out testH[4]);
            hypothesis.TryGetValue("t-test for differences in the means of two samples", out testH[5]);
            hypothesis.TryGetValue("F-test for differences in the variances of two samples", out testH[6]);

            var validH = new double[] { 0.0012195d, 0.8204, 0.2436, 0.5956, 0.7757, 0.4691, 0.2287 };
            for (int i = 0; i < validH.Length; i++)
            {
                Assert.AreEqual(validH[i], testH[i], 1E-1);
            }
        }

        /// <summary>
        /// Test the methods that return monthly statistics. These methods include MonthlyPercentiles() and MonthlySummaryStatistics().
        /// </summary>
        [TestMethod]
        public void Test_MonthlyStats()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13, 38, 17, 3, 6, 27, 11, 2, 41, 44, 37, 50, 8 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            var percentiles = ts.MonthlyPercentiles(new double[] { 0.10, 0.20, 0.30, 0.40, 0.60, 0.70, 0.80, 0.90 });
            var validP = new double[,] { { 23.6, 25.2, 26.8, 28.4, 31.6, 33.2, 34.8, 36.4 },
                                         { 16.1, 16.2, 16.3, 16.4, 16.6, 16.7, 16.8, 16.9 },
                                         { 6, 9, 12, 15, 21, 24, 27, 30 },
                                         { 5.1, 5.2, 5.3, 5.4, 5.6, 5.7, 5.8, 5.9 },
                                         { 13.5, 15.0, 16.5, 18.0, 21.0, 22.5, 24.0, 25.5 },
                                         { 13.5, 16.0, 18.5, 21.0, 26.0, 28.5, 31.0, 33.5 },
                                         { 6.6, 11.2, 15.8, 20.4, 29.6, 34.2, 38.8, 43.4 },
                                         { 13.1, 16.2, 19.3, 22.4, 28.6, 31.7, 34.8, 37.9 },
                                         { 20.6, 23.2, 25.8, 28.4, 33.6, 36.2, 38.8, 41.4 },
                                         { 17.2, 19.4, 21.6, 23.8, 28.2, 30.4, 32.6, 34.8 },
                                         { 24.8, 27.6, 30.4, 33.2, 38.8, 41.6, 44.4, 47.2 },
                                         { 8.5,  9.0,  9.5, 10.0, 11.0, 11.5, 12.0, 12.5 } };

            var summary = ts.MonthlySummaryStatistics();
            var validS = new double[,] { { 22, 22.8, 26.0, 30.0, 34.0, 37.2, 38, 30 },
                                         { 16, 16.05, 16.25, 16.50, 16.75, 16.95, 17, 16.5 },
                                         { 3, 4.5, 10.5, 18.0, 25.5, 31.5, 33, 18 },
                                         { 5, 5.05, 5.25, 5.50, 5.75, 5.95, 6, 5.5 },
                                         { 12, 12.75, 15.75, 19.50, 23.25, 26.25, 27, 19.5},
                                         { 11, 12.25, 17.25, 23.50, 29.75, 34.75 , 36, 23.5 },
                                         { 2, 4.3, 13.5, 25.0, 36.5, 45.7, 48, 25},
                                         { 10, 11.55, 17.75, 25.50, 33.25, 39.45, 41, 25.5},
                                         { 18, 19.3, 24.5, 31.0, 37.5, 42.7, 44, 31 },
                                         { 15, 16.1, 20.5, 26.0, 31.5, 35.9, 37, 26 },
                                         { 22, 23.4, 29.0, 36.0, 43.0, 48.6, 50, 36 },
                                         { 8, 8.25, 9.25, 10.50, 11.75, 12.75, 13, 10.5 } };

            for (int i = 0; i < 11; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    Assert.AreEqual(validP[i, j], percentiles[i, j], 1E-10);
                    Assert.AreEqual(validS[i, j], summary[i, j], 1E-10);
                }
            }
        }

        /// <summary>
        /// Test that the mean column of MonthlySummaryStatistics() is the sequential arithmetic
        /// mean of each month's values, reproduced bit-for-bit by Statistics.Mean.
        /// </summary>
        /// <remarks>
        /// The summary sorts each month's values ascending before reducing them, so the oracle
        /// applies the same order. Sequential accumulation makes the reduction deterministic on
        /// every host regardless of processor count, and the delta of zero detects any
        /// reassociation of the summation order.
        /// </remarks>
        [TestMethod]
        public void Test_MonthlySummaryStats_MeanColumnIsSequentialMean()
        {
            // Three years of monthly values spanning several orders of magnitude so a change in
            // summation order would alter the last bits of the mean.
            var values = new double[36];
            for (int i = 0; i < values.Length; i++)
            {
                values[i] = Math.Pow(10d, (i % 5) - 2) * (1d + i / 35d);
            }
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2021, 01, 01), values);

            var summary = ts.MonthlySummaryStatistics();
            for (int month = 1; month <= 12; month++)
            {
                var monthlyData = new List<double>();
                for (int j = 0; j < ts.Count; j++)
                {
                    if (ts[j].Index.Month == month) { monthlyData.Add(ts[j].Value); }
                }
                monthlyData.Sort();
                Assert.AreEqual(Numerics.Data.Statistics.Statistics.Mean(monthlyData), summary[month - 1, 7], 0d);
            }
        }

        /// <summary>
        /// Test the monthly frequency method
        /// </summary>
        [TestMethod]
        public void Test_MonthlyFrequency()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10, 18, 15, 22, 13, 38, 17, 3, 6, 27, 11, 2, 41, 44, 37, 50, 8 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            var frequencies = ts.MonthlyFrequency();
            for (int i = 0; i < frequencies.Length; i++)
            {
                Assert.AreEqual(2, frequencies[i]);
            }
        }

        /// <summary>
        /// Test the calendar year series method
        /// </summary>
        [TestMethod]
        public void Test_CalendarYearSeries()
        {
            var values = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            var annualTS = ts.CalendarYearSeries();

            // same values as annual max series test
            var maxVals = new double[] { 263, 331, 317, 337, 262, 239 };

            for (int i = 0; i < maxVals.Length; i++)
            {
                Assert.AreEqual(maxVals[i], annualTS[i].Value);
            }
        }

        /// <summary>
        /// Test the water year series method
        /// </summary>
        [TestMethod]
        public void Test_WaterYearSeries()
        {
            var values = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            var customTS = ts.CustomYearSeries();
            var maxVals = new double[] { 263, 331, 317, 204, 337, 250 };
            for (int i = 0; i < maxVals.Length; i++)
            {
                Assert.AreEqual(maxVals[i], customTS[i].Value);
            }
        }

        /// <summary>
        /// Test the custom year series method
        /// </summary>
        [TestMethod]
        public void Test_CustomYearSeries()
        {
            var values = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            // First test calendar year
            var customTS = ts.CustomYearSeries(1, 12);
            var maxVals = new double[] { 263, 331, 317, 337, 262, 239 };
            for (int i = 0; i < maxVals.Length; i++)
            {
                Assert.AreEqual(maxVals[i], customTS[i].Value);
            }

            // Next test custom overlapping year
            customTS = ts.CustomYearSeries(10, 3);
            maxVals = new double[] { 244, 331, 240, 204, 337, 250 };
            for (int i = 0; i < maxVals.Length; i++)
            {
                Assert.AreEqual(maxVals[i], customTS[i].Value);
            }
        }

        /// <summary>
        /// Test the monthly series method
        /// </summary>
        [TestMethod]
        public void Test_MonthlySeries()
        {
            var values = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);
            var monthlyTS = ts.MonthlySeries();

            for (int i = 0; i < monthlyTS.Count; i++)
            {
                // since the time interval is one month, all of the original values are the max of their own month
                Assert.AreEqual(values[i], monthlyTS[i].Value);
            }

            // with a weekly time interval the maxes differ
            var ts2 = new TimeSeries(TimeInterval.SevenDay, new DateTime(2023, 01, 01), values);
            var monthlyTS2 = ts2.MonthlySeries();
            var maxVals2 = new double[] { 244, 263, 205, 331, 262, 240, 195, 317, 204, 185, 337, 262, 182, 250, 239, 180 };
            for (int i = 0; i < maxVals2.Length; i++)
            {
                Assert.AreEqual(maxVals2[i], monthlyTS2[i].Value);
            }
        }

        /// <summary>
        /// Test the quarterly series method
        /// </summary>
        [TestMethod]
        public void Test_QuarterlySeries()
        {
            var values = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            var quarterlyTS = ts.QuarterlySeries();
            var maxVals = new double[] { 244, 229, 263, 205, 331, 174, 262, 240, 230, 173, 317, 204, 194, 179, 192, 337, 230, 262, 182, 250, 237, 239, 180 };
            for (int i = 0; i < maxVals.Length; i++)
            {
                Assert.AreEqual(maxVals[i], quarterlyTS[i].Value);
            }
        }

        /// <summary>
        /// Test using the 'clust' method included in the POT R package (https://cran.r-project.org/web/packages/POT/index.html)
        /// </summary>
        [TestMethod]
        public void Test_PeaksOverThreshold()
        {
            var values = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166, 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var ts = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 01, 01), values);

            // Threshold of 100, 2 min steps between
            var truePOT = new double[] { 331, 337, 262 };
            var POT = ts.PeaksOverThresholdSeries(100, 2);
            for (int i = 0; i < POT.Count; i++)
            {
                Assert.AreEqual(truePOT[i], POT[i].Value);
            }

            // Threshold of 90, 1 min steps between
            truePOT = new double[] { 337 };
            POT = ts.PeaksOverThresholdSeries(90, 1);
            for (int i = 0; i < POT.Count; i++)
            {
                Assert.AreEqual(truePOT[i], POT[i].Value);
            }

            // Threshold of 150, 5 min steps between
            truePOT = new double[] { 331, 240, 317, 337 };
            POT = ts.PeaksOverThresholdSeries(150, 5);
            for (int i = 0; i < POT.Count; i++)
            {
                Assert.AreEqual(truePOT[i], POT[i].Value);
            }

            // Threshold of 200, 2 min steps between
            truePOT = new double[] { 263, 331, 262, 317, 337, 250 };
            POT = ts.PeaksOverThresholdSeries(200, 2);
            for (int i = 0; i < POT.Count; i++)
            {
                Assert.AreEqual(truePOT[i], POT[i].Value);
            }

        }

        /// <summary>
        /// Test the SeasonalDecompose method with a known sinusoidal + linear trend signal.
        /// </summary>
        [TestMethod]
        public void Test_SeasonalDecompose()
        {
            // Create a time series: linear trend + sinusoidal seasonal + small noise
            int period = 12;
            int nYears = 5;
            int n = period * nYears; // 60 data points
            var ts = new TimeSeries(TimeInterval.OneMonth);
            var startDate = new DateTime(2000, 1, 1);

            double[] truetrend = new double[n];
            double[] trueSeasonal = new double[n];

            for (int i = 0; i < n; i++)
            {
                // Linear trend: 100 + 0.5 * i
                truetrend[i] = 100.0 + 0.5 * i;
                // Seasonal: 10 * sin(2π * i / 12)
                trueSeasonal[i] = 10.0 * Math.Sin(2.0 * Math.PI * i / period);
                double value = truetrend[i] + trueSeasonal[i];
                ts.Add(new SeriesOrdinate<DateTime, double>(startDate.AddMonths(i), value));
            }

            // Decompose
            var (trend, seasonal, residual) = ts.SeasonalDecompose(period);

            // Verify trend is returned
            Assert.IsGreaterThan(0, trend.Count, "Trend should have values.");
            Assert.IsLessThanOrEqualTo(n, trend.Count, "Trend should not exceed original length.");

            // Verify seasonal has correct length
            // HasCount does not accept double[] (not IEnumerable), so suppress MSTEST0037
#pragma warning disable MSTEST0037
            Assert.AreEqual(n, seasonal.Length, "Seasonal should have same length as original.");
#pragma warning restore MSTEST0037

            // Verify residual is returned
            Assert.IsGreaterThan(0, residual.Count, "Residual should have values.");

            // Verify the decomposition is additive: original = trend + seasonal + residual
            // This is an exact mathematical identity by construction
            for (int i = 0; i < residual.Count; i++)
            {
                var residOrd = residual[i];
                var trendOrd = trend.FirstOrDefault(t => t.Index == residOrd.Index);
                int origIdx = -1;
                for (int k = 0; k < n; k++)
                {
                    if (ts[k].Index == residOrd.Index) { origIdx = k; break; }
                }
                if (trendOrd != null && origIdx >= 0)
                {
                    double reconstructed = trendOrd.Value + seasonal[origIdx] + residOrd.Value;
                    Assert.AreEqual(ts[origIdx].Value, reconstructed, 1E-6,
                        $"Additive decomposition should hold at index {origIdx}.");
                }
            }

            // The seasonal component should capture the sinusoidal signal
            // Verify seasonal values are not all zero (they should have nonzero amplitude)
            double maxSeasonal = seasonal.Max();
            double minSeasonal = seasonal.Min();
            Assert.IsGreaterThan(1.0, maxSeasonal - minSeasonal,
                "Seasonal component should have nonzero amplitude for sinusoidal input.");
        }

        /// <summary>
        /// Test that SeasonalDecompose throws for invalid inputs.
        /// </summary>
        [TestMethod]
        public void Test_SeasonalDecompose_InvalidInputs()
        {
            var ts = new TimeSeries(TimeInterval.OneMonth);
            var startDate = new DateTime(2000, 1, 1);
            for (int i = 0; i < 10; i++)
                ts.Add(new SeriesOrdinate<DateTime, double>(startDate.AddMonths(i), i));

            // Period too small
            Assert.Throws<ArgumentException>(() => ts.SeasonalDecompose(1));

            // Too few data points for 2 complete periods
            Assert.Throws<ArgumentException>(() => ts.SeasonalDecompose(12));
        }

        #region ResampleWithKNN tests

        /// <summary>
        /// Build a stationary AR(1) series x[t] = phi*x[t-1] + N(0, sigma^2) using a fixed seed.
        /// Used by the ResampleWithKNN/ResampleWithBlockBootstrap statistical tests.
        /// </summary>
        private static TimeSeries MakeAR1(int n, double phi, double sigma, int seed)
        {
            var rng = new Random(seed);
            var values = new double[n];
            values[0] = 0.0;
            for (int i = 1; i < n; i++)
            {
                // Box-Muller for a standard normal draw.
                double u1 = 1.0 - rng.NextDouble();
                double u2 = 1.0 - rng.NextDouble();
                double z = Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
                values[i] = phi * values[i - 1] + sigma * z;
            }
            return new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1), values);
        }

        /// <summary>
        /// Same seed must produce identical KNN resamples.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithKNN_DeterministicWithSeed()
        {
            var ts = MakeAR1(n: 100, phi: 0.5, sigma: 1.0, seed: 1);
            var a = ts.ResampleWithKNN(timeSteps: 50, k: 10, seed: 42);
            var b = ts.ResampleWithKNN(timeSteps: 50, k: 10, seed: 42);
            Assert.AreEqual(a.Count, b.Count);
            for (int i = 0; i < a.Count; i++)
                Assert.AreEqual(a[i].Value, b[i].Value, 1e-12,
                    $"KNN must be deterministic with a fixed seed; index {i} differs.");
        }

        /// <summary>
        /// Output length must equal the requested timeSteps exactly.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithKNN_OutputLength()
        {
            var ts = MakeAR1(n: 100, phi: 0.5, sigma: 1.0, seed: 1);
            foreach (int len in new[] { 1, 7, 50, 250 })
            {
                var r = ts.ResampleWithKNN(timeSteps: len, k: 10, seed: 42);
                Assert.AreEqual(len, r.Count,
                    $"KNN output length must match the requested timeSteps; expected {len}, got {r.Count}.");
            }
        }

        /// <summary>
        /// On a strongly trended series x[t] = t, the conditional KNN bootstrap (Lall-Sharma)
        /// advances through the trend, accumulating mean drift of ≈ +1 per step: each step
        /// must take the neighbor's successor value x[j+1], never the neighbor's own value,
        /// or the trajectory hovers near the starting value (zero net drift). Drift is
        /// averaged over multiple seeds to suppress the random-walk noise.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithKNN_AdvancesThroughTime()
        {
            // x[t] = t, t in [0, 199]
            var values = Enumerable.Range(0, 200).Select(i => (double)i).ToArray();
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1), values);

            const int steps = 50;
            const int trials = 20;
            double avgDrift = 0;
            for (int seed = 1; seed <= trials; seed++)
            {
                var r = ts.ResampleWithKNN(timeSteps: steps, k: 10, seed: seed);
                avgDrift += (r[steps - 1].Value - r[0].Value);
            }
            avgDrift /= trials;

            // Successor sampling gives avgDrift ≈ +50 (+1 per step over 50 steps);
            // neighbor-value sampling gives avgDrift ~ N(0, ~5) and fails clearly.
            Assert.IsGreaterThan(25.0, avgDrift,
                $"Expected KNN trajectory to advance through the trend (avgDrift > 25 over {steps} steps); " +
                $"observed avgDrift = {avgDrift:F2}. Sampling the neighbor's own value keeps the trajectory near the starting value.");
        }

        /// <summary>
        /// The marginal mean of resampled tails should match the input mean within sampling error.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithKNN_PreservesMarginalMean()
        {
            var ts = MakeAR1(n: 200, phi: 0.5, sigma: 1.0, seed: 7);
            double inputMean = ts.MeanValue();
            double inputStd = ts.StandardDeviation();

            const int trials = 200;
            const int steps = 100;
            double sumOfTailMeans = 0;
            for (int seed = 1; seed <= trials; seed++)
            {
                var r = ts.ResampleWithKNN(timeSteps: steps, k: 10, seed: seed);
                sumOfTailMeans += r.MeanValue();
            }
            double grandMean = sumOfTailMeans / trials;
            double standardError = inputStd / Math.Sqrt(steps * trials);

            // Allow 4x SE, generous to absorb short-series boundary effects in KNN.
            Assert.IsLessThan(4.0 * standardError + 0.05, Math.Abs(grandMean - inputMean),
                $"KNN tails should preserve the marginal mean; input={inputMean:F4}, grand-mean={grandMean:F4}, " +
                $"4*SE={4.0 * standardError:F4}.");
        }

        /// <summary>
        /// The marginal variance of resampled tails should be in the same regime as the input variance.
        /// KNN can shrink variance because the conditional sampling pulls toward the local mean, so the test
        /// only requires the resample variance to remain within 50% of the input variance.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithKNN_PreservesMarginalVariance()
        {
            var ts = MakeAR1(n: 200, phi: 0.5, sigma: 1.0, seed: 7);
            double inputVar = ts.StandardDeviation() * ts.StandardDeviation();

            const int trials = 200;
            const int steps = 100;
            double sumOfTailVars = 0;
            for (int seed = 1; seed <= trials; seed++)
            {
                var r = ts.ResampleWithKNN(timeSteps: steps, k: 10, seed: seed);
                double rStd = r.StandardDeviation();
                sumOfTailVars += rStd * rStd;
            }
            double avgTailVar = sumOfTailVars / trials;
            double ratio = avgTailVar / inputVar;

            Assert.IsTrue(ratio > 0.5 && ratio < 2.0,
                $"KNN tail variance should be within 0.5x..2x of input variance; input={inputVar:F4}, " +
                $"avg-tail={avgTailVar:F4}, ratio={ratio:F2}.");
        }

        /// <summary>
        /// On a synthetic AR(1) with phi=0.7, the resampled tails should recover lag-1 autocorrelation
        /// in the same regime. The conditional bootstrap is what gives KNN this property.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithKNN_RecoversAR1Autocorrelation()
        {
            const double phiTrue = 0.7;
            var ts = MakeAR1(n: 300, phi: phiTrue, sigma: 1.0, seed: 7);

            const int trials = 200;
            const int steps = 100;
            double sumLag1 = 0;
            for (int seed = 1; seed <= trials; seed++)
            {
                var r = ts.ResampleWithKNN(timeSteps: steps, k: 10, seed: seed);
                sumLag1 += Lag1Autocorrelation(r);
            }
            double avgLag1 = sumLag1 / trials;

            // The resampler must recover the AR(1) persistence: avgLag1 above 0.4 (near
            // phiTrue = 0.7 with shrinkage), not the near-zero autocorrelation of a
            // KNN-neighbor random walk about a fixed point.
            Assert.IsTrue(avgLag1 > 0.4 && avgLag1 < 0.95,
                $"KNN should recover lag-1 autocorrelation in the AR(1) regime; expected ~{phiTrue}, got {avgLag1:F3}.");
        }

        private static double Lag1Autocorrelation(TimeSeries ts)
        {
            int n = ts.Count;
            if (n < 3) return double.NaN;
            double mean = ts.MeanValue();
            double num = 0, den = 0;
            for (int i = 1; i < n; i++) num += (ts[i].Value - mean) * (ts[i - 1].Value - mean);
            for (int i = 0; i < n; i++) den += (ts[i].Value - mean) * (ts[i].Value - mean);
            return den > 0 ? num / den : double.NaN;
        }

        /// <summary>
        /// KNN must reject series that are too short for a sensible conditional bootstrap.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithKNN_ThrowsForShortSeries()
        {
            // Count = 10: one less than required since the last index is excluded from candidates
            // and the underlying KNearestNeighbors needs ≥10 training points.
            var ts = MakeAR1(n: 10, phi: 0.5, sigma: 1.0, seed: 1);
            Assert.Throws<ArgumentException>(() => ts.ResampleWithKNN(timeSteps: 5, k: 5, seed: 42));
        }

        #endregion

        #region ResampleWithBlockBootstrap tests

        /// <summary>
        /// Same seed must produce identical block-bootstrap resamples.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithBlockBootstrap_DeterministicWithSeed()
        {
            var ts = MakeAR1(n: 100, phi: 0.5, sigma: 1.0, seed: 1);
            var a = ts.ResampleWithBlockBootstrap(timeSteps: 50, blockSize: 5, seed: 42);
            var b = ts.ResampleWithBlockBootstrap(timeSteps: 50, blockSize: 5, seed: 42);
            Assert.AreEqual(a.Count, b.Count);
            for (int i = 0; i < a.Count; i++)
                Assert.AreEqual(a[i].Value, b[i].Value, 1e-12,
                    $"BlockBootstrap must be deterministic with a fixed seed; index {i} differs.");
        }

        /// <summary>
        /// Output length must equal the requested timeSteps exactly even when not divisible by blockSize.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithBlockBootstrap_OutputLength()
        {
            var ts = MakeAR1(n: 100, phi: 0.5, sigma: 1.0, seed: 1);
            // Try lengths that are NOT multiples of blockSize=7 to verify the truncation logic.
            foreach (int len in new[] { 1, 7, 13, 50, 251 })
            {
                var r = ts.ResampleWithBlockBootstrap(timeSteps: len, blockSize: 7, seed: 42);
                Assert.AreEqual(len, r.Count,
                    $"BlockBootstrap output length must match requested timeSteps; expected {len}, got {r.Count}.");
            }
        }

        /// <summary>
        /// The marginal mean of block-bootstrap tails should match the input mean within sampling error.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithBlockBootstrap_PreservesMarginalMean()
        {
            var ts = MakeAR1(n: 200, phi: 0.5, sigma: 1.0, seed: 7);
            double inputMean = ts.MeanValue();
            double inputStd = ts.StandardDeviation();

            const int trials = 200;
            const int steps = 100;
            double sumOfTailMeans = 0;
            for (int seed = 1; seed <= trials; seed++)
            {
                var r = ts.ResampleWithBlockBootstrap(timeSteps: steps, blockSize: 5, seed: seed);
                sumOfTailMeans += r.MeanValue();
            }
            double grandMean = sumOfTailMeans / trials;
            double standardError = inputStd / Math.Sqrt(steps * trials);

            Assert.IsLessThan(4.0 * standardError + 0.05, Math.Abs(grandMean - inputMean),
                $"BlockBootstrap tails should preserve the marginal mean; input={inputMean:F4}, " +
                $"grand-mean={grandMean:F4}, 4*SE={4.0 * standardError:F4}.");
        }

        /// <summary>
        /// The marginal variance of block-bootstrap tails should be close to the input variance.
        /// Block bootstrap preserves marginal variance well; the test requires within ±20%.
        /// </summary>
        [TestMethod]
        public void Test_ResampleWithBlockBootstrap_PreservesMarginalVariance()
        {
            var ts = MakeAR1(n: 200, phi: 0.5, sigma: 1.0, seed: 7);
            double inputVar = ts.StandardDeviation() * ts.StandardDeviation();

            const int trials = 200;
            const int steps = 100;
            double sumOfTailVars = 0;
            for (int seed = 1; seed <= trials; seed++)
            {
                var r = ts.ResampleWithBlockBootstrap(timeSteps: steps, blockSize: 5, seed: seed);
                double rStd = r.StandardDeviation();
                sumOfTailVars += rStd * rStd;
            }
            double avgTailVar = sumOfTailVars / trials;
            double ratio = avgTailVar / inputVar;

            Assert.IsTrue(ratio > 0.8 && ratio < 1.2,
                $"BlockBootstrap should preserve marginal variance within ±20%; input={inputVar:F4}, " +
                $"avg-tail={avgTailVar:F4}, ratio={ratio:F2}.");
        }

        #endregion

        #region Collection Behavior

        /// <summary>
        /// Clear should empty the series and raise a single reset event.
        /// </summary>
        /// <remarks>
        /// Clearing must raise a single Reset rather than per-element notifications;
        /// element-by-element removal costs two full equality scans per element
        /// (O(n²) on large downloads).
        /// </remarks>
        [TestMethod]
        public void Test_Clear_EmptiesSeriesAndRaisesSingleReset()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 1, 1), new double[] { 1, 2, 3, 4, 5 });
            var actions = new List<System.Collections.Specialized.NotifyCollectionChangedAction>();
            ts.CollectionChanged += (s, e) => actions.Add(e.Action);

            ts.Clear();

            Assert.AreEqual(0, ts.Count);
            Assert.HasCount(1, actions);
            Assert.AreEqual(System.Collections.Specialized.NotifyCollectionChangedAction.Reset, actions[0]);
        }

        /// <summary>
        /// Clearing a large series should complete quickly.
        /// </summary>
        /// <remarks>
        /// Clearing a century of daily data must complete immediately; element-by-element
        /// removal takes seconds to minutes at this size.
        /// </remarks>
        [TestMethod]
        [Timeout(5000, CooperativeCancellation = true)]
        public void Test_Clear_LargeSeries_CompletesQuickly()
        {
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(1900, 1, 1), new double[200000]);

            ts.Clear();

            Assert.AreEqual(0, ts.Count);
        }

        /// <summary>
        /// RemoveAt should remove the ordinate at the requested index even when an equal-valued
        /// ordinate appears earlier in the series.
        /// </summary>
        /// <remarks>
        /// RemoveAt must remove by position: delegating to Remove(item) removes the first equal
        /// element, silently deleting the wrong ordinate when duplicates exist.
        /// </remarks>
        [TestMethod]
        public void Test_RemoveAt_WithDuplicateOrdinates_RemovesRequestedIndex()
        {
            var first = new SeriesOrdinate<DateTime, double>(new DateTime(2024, 1, 1), 1.0);
            var middle = new SeriesOrdinate<DateTime, double>(new DateTime(2024, 1, 2), 2.0);
            var duplicateOfFirst = new SeriesOrdinate<DateTime, double>(new DateTime(2024, 1, 1), 1.0);
            var ts = new TimeSeries(TimeInterval.Irregular) { first, middle, duplicateOfFirst };

            var removals = new List<(object item, int index)>();
            ts.CollectionChanged += (s, e) =>
            {
                if (e.Action == System.Collections.Specialized.NotifyCollectionChangedAction.Remove)
                    removals.Add((e.OldItems?[0], e.OldStartingIndex));
            };

            ts.RemoveAt(2);

            Assert.AreEqual(2, ts.Count);
            Assert.IsTrue(ReferenceEquals(first, ts[0]), "The first ordinate should remain at index 0.");
            Assert.IsTrue(ReferenceEquals(middle, ts[1]), "The middle ordinate should remain at index 1.");
            Assert.HasCount(1, removals);
            Assert.IsTrue(ReferenceEquals(duplicateOfFirst, removals[0].item), "The event should carry the removed ordinate.");
            Assert.AreEqual(2, removals[0].index, "The event should report the requested index.");
        }

        #endregion

        /// <summary>
        /// Verify SmoothedSeries reproduces the exact preprocessing PeaksOverThresholdSeries applies.
        /// </summary>
        /// <remarks>
        /// The smoothing branch was extracted from PeaksOverThresholdSeries so threshold-selection
        /// diagnostics can operate on the same series the extraction thresholds; this test pins the
        /// branch behavior (moving average and moving sum smooth for periods above 1 and clone at
        /// a period of 1, differencing applies at every period, and None clones) and that the
        /// smoothed values genuinely differ from the raw values when smoothing is configured.
        /// </remarks>
        [TestMethod]
        public void Test_SmoothedSeries_MatchesPeaksOverThresholdPreprocessing()
        {
            var series = new TimeSeries(TimeInterval.OneDay, new DateTime(2020, 1, 1), new double[] { 5, 9, 2, 14, 7, 11, 3, 16, 8, 12 });

            // Moving average over 3 steps matches the direct transform and differs from the raw values.
            var smoothed = series.SmoothedSeries(SmoothingFunctionType.MovingAverage, 3);
            var direct = series.MovingAverage(3);
            Assert.AreEqual(direct.Count, smoothed.Count);
            bool anyDifferent = false;
            for (int i = 0; i < direct.Count; i++)
            {
                Assert.AreEqual(direct[i].Value, smoothed[i].Value, 0d);
                if (!double.IsNaN(smoothed[i].Value) && smoothed[i].Value != series[i].Value)
                    anyDifferent = true;
            }
            Assert.IsTrue(anyDifferent, "Smoothing must change the diagnostic value scale.");

            // A period of 1 is an identity clone for moving average and moving sum.
            var identity = series.SmoothedSeries(SmoothingFunctionType.MovingAverage, 1);
            for (int i = 0; i < series.Count; i++)
                Assert.AreEqual(series[i].Value, identity[i].Value, 0d);

            // Moving sum and differencing route to their transforms.
            var movingSum = series.SmoothedSeries(SmoothingFunctionType.MovingSum, 3);
            var directSum = series.MovingSum(3);
            for (int i = 0; i < directSum.Count; i++)
                Assert.AreEqual(directSum[i].Value, movingSum[i].Value, 0d);
            var difference = series.SmoothedSeries(SmoothingFunctionType.Difference, 1);
            var directDifference = series.Difference(1);
            for (int i = 0; i < directDifference.Count; i++)
                Assert.AreEqual(directDifference[i].Value, difference[i].Value, 0d);

            // None clones the series.
            var none = series.SmoothedSeries(SmoothingFunctionType.None);
            for (int i = 0; i < series.Count; i++)
                Assert.AreEqual(series[i].Value, none[i].Value, 0d);
        }

        /// <summary>
        /// The direct clone matches an XML round trip (the 2.2.0 clone, kept as a bitwise reference) bit for bit,
        /// and is independent of the source.
        /// </summary>
        [TestMethod]
        public void Test_Clone_MatchesXmlRoundTrip()
        {
            var ts = new TimeSeries(TimeInterval.OneHour);
            var start = new DateTime(2021, 3, 14, 0, 0, 0);
            var values = new double[] { 1.5, -0.0, 0.0, double.NaN, BitConverter.Int64BitsToDouble(0x7FF8000000000001), double.PositiveInfinity,
                double.NegativeInfinity, double.Epsilon, double.MaxValue, 0.1 + 0.2, Math.PI * 1e-300, -123456.789 };
            var kinds = new[] { DateTimeKind.Unspecified, DateTimeKind.Utc, DateTimeKind.Local };
            for (int i = 0; i < values.Length; i++)
                ts.Add(new SeriesOrdinate<DateTime, double>(DateTime.SpecifyKind(start.AddMinutes(37 * i).AddTicks(i), kinds[i % 3]), values[i]));
            // A local time that may fall in a daylight-saving gap.
            ts.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2021, 3, 14, 2, 30, 0, DateTimeKind.Local), 2));

            var expected = new TimeSeries(ts.ToXElement());
            var clone = ts.Clone();
            Assert.AreEqual(expected.TimeInterval, clone.TimeInterval);
            Assert.AreEqual(expected.Count, clone.Count);
            for (int i = 0; i < expected.Count; i++)
            {
                Assert.AreEqual(expected[i].Index.Ticks, clone[i].Index.Ticks);
                Assert.AreEqual(expected[i].Index.Kind, clone[i].Index.Kind);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i].Value), BitConverter.DoubleToInt64Bits(clone[i].Value));
                Assert.AreNotSame(ts[i], clone[i]);
            }

            clone[0].Value = 99;
            Assert.AreEqual(1.5, ts[0].Value);
            Assert.AreEqual(0, new TimeSeries(TimeInterval.Irregular).Clone().Count);
            Assert.AreEqual(TimeInterval.Irregular, new TimeSeries(TimeInterval.Irregular).Clone().TimeInterval);
        }

        /// <summary>
        /// The direct clone matches the XML round trip on random values of every magnitude.
        /// </summary>
        [TestMethod]
        public void Test_Clone_MatchesXmlRoundTrip_Random()
        {
            var rnd = new Random(42);
            var bytes = new byte[8];
            var data = new double[5000];
            for (int i = 0; i < data.Length; i++)
            {
                rnd.NextBytes(bytes);
                data[i] = BitConverter.ToDouble(bytes, 0);
            }
            var ts = new TimeSeries(TimeInterval.FifteenMinute, new DateTime(1950, 1, 1), data);
            var expected = new TimeSeries(ts.ToXElement());
            var clone = ts.Clone();
            for (int i = 0; i < data.Length; i++)
            {
                Assert.AreEqual(expected[i].Index, clone[i].Index);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i].Value), BitConverter.DoubleToInt64Bits(clone[i].Value));
            }
        }

        /// <summary>
        /// Missing values are interpolated on exact tick differences. Expected values are the exact rational
        /// interpolants (Python fractions) rounded to double; interpolation on OLE Automation days would miss
        /// the first two by about 7E-14 relative.
        /// </summary>
        [TestMethod]
        public void Test_InterpolateMissingData_ExactTimeDifferences()
        {
            var hourly = new TimeSeries(TimeInterval.OneHour, new DateTime(2020, 6, 1, 12, 0, 0), new[] { 1000.123456789, double.NaN, double.NaN, 1003.987654321 });
            hourly.InterpolateMissingData(5);
            Assert.AreEqual(1001.411522633, hourly[1].Value, 1001.411522633 * 4E-16);
            Assert.AreEqual(1002.699588477, hourly[2].Value, 1002.699588477 * 4E-16);

            // A gap that runs to the end is extrapolated from the two preceding ordinates.
            var tail = new TimeSeries(TimeInterval.OneHour, new DateTime(2020, 6, 1), new[] { 0.1, 0.7, double.NaN, double.NaN });
            tail.InterpolateMissingData(5);
            Assert.AreEqual(1.2999999999999998, tail[2].Value, 1.3 * 4E-16);
            Assert.AreEqual(1.9, tail[3].Value, 1.9 * 4E-16);

            // Irregular spacing: 7 of 30 minutes.
            var irregular = new TimeSeries(TimeInterval.Irregular);
            irregular.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2020, 6, 1, 0, 0, 0), 2.0));
            irregular.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2020, 6, 1, 0, 7, 0), double.NaN));
            irregular.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2020, 6, 1, 0, 30, 0), 9.0));
            irregular.InterpolateMissingData(1);
            Assert.AreEqual(3.6333333333333333, irregular[1].Value, 3.6333333333333333 * 4E-16);

            // The indexed overload uses the same interpolation.
            var indexed = new TimeSeries(TimeInterval.OneHour, new DateTime(2020, 6, 1, 12, 0, 0), new[] { 1000.123456789, double.NaN, double.NaN, 1003.987654321 });
            indexed.InterpolateMissingData(5, new[] { 1, 2 });
            Assert.AreEqual(hourly[1].Value, indexed[1].Value, 0d);
            Assert.AreEqual(hourly[2].Value, indexed[2].Value, 0d);
        }

        /// <summary>
        /// Irregular readings to hourly blocks match pandas 2.3.3 <c>resample('h', label='left', closed='left')</c>
        /// for mean, sum with min_count 1 and 2, min and max, including an empty hour and a missing value.
        /// </summary>
        [TestMethod]
        public void Test_AggregateByInterval_MatchesPandas()
        {
            var irregular = new TimeSeries(TimeInterval.Irregular);
            var times = new[] { "00:05", "00:20", "00:50", "01:10", "03:40", "03:55", "03:59" };
            var values = new[] { 1.5, 2.25, double.NaN, 4.0, 10.0, 11.5, 0.125 };
            // Added out of order; the method sorts a copy.
            foreach (int i in new[] { 4, 0, 2, 1, 3, 6, 5 })
                irregular.Add(new SeriesOrdinate<DateTime, double>(DateTime.Parse("2024-03-10 " + times[i], CultureInfo.InvariantCulture), values[i]));

            var mean = irregular.AggregateByInterval(TimeInterval.OneHour);
            Assert.AreEqual(TimeInterval.OneHour, mean.TimeInterval);
            CollectionAssert.AreEqual(new[] { new DateTime(2024, 3, 10, 0, 0, 0), new DateTime(2024, 3, 10, 1, 0, 0), new DateTime(2024, 3, 10, 2, 0, 0), new DateTime(2024, 3, 10, 3, 0, 0) },
                mean.IndexesToArray());
            AssertValues(new[] { 1.875, 4.0, double.NaN, 7.208333333333333 }, mean);
            AssertValues(new[] { 3.75, 4.0, double.NaN, 21.625 }, irregular.AggregateByInterval(TimeInterval.OneHour, BlockFunctionType.Sum));
            AssertValues(new[] { 3.75, double.NaN, double.NaN, 21.625 }, irregular.AggregateByInterval(TimeInterval.OneHour, BlockFunctionType.Sum, 2));
            AssertValues(new[] { 1.5, 4.0, double.NaN, 0.125 }, irregular.AggregateByInterval(TimeInterval.OneHour, BlockFunctionType.Minimum));
            AssertValues(new[] { 2.25, 4.0, double.NaN, 11.5 }, irregular.AggregateByInterval(TimeInterval.OneHour, BlockFunctionType.Maximum));
            Assert.AreEqual(new DateTime(2024, 3, 10, 3, 40, 0), irregular[0].Index, "The series is not modified.");
        }

        /// <summary>
        /// An hourly series starting at 05:00 aggregates to calendar days (pandas <c>resample('D').mean()</c>: 10, 25),
        /// where ConvertTimeInterval counts 24-step blocks from 05:00. Quarters match pandas <c>resample('QS').sum(min_count=1)</c>.
        /// </summary>
        [TestMethod]
        public void Test_AggregateByInterval_CalendarAlignment()
        {
            var hourly = new TimeSeries(TimeInterval.OneHour, new DateTime(2020, 1, 1, 5, 0, 0), Enumerable.Range(1, 30).Select(i => (double)i).ToArray());
            var daily = hourly.AggregateByInterval(TimeInterval.OneDay);
            CollectionAssert.AreEqual(new[] { new DateTime(2020, 1, 1), new DateTime(2020, 1, 2) }, daily.IndexesToArray());
            AssertValues(new[] { 10.0, 25.0 }, daily);

            var sixHour = hourly.AggregateByInterval(TimeInterval.SixHour, BlockFunctionType.Sum);
            Assert.AreEqual(new DateTime(2020, 1, 1, 0, 0, 0), sixHour[0].Index);
            Assert.AreEqual(1d, sixHour[0].Value, "Only 05:00 falls in the 00:00 block.");
            Assert.AreEqual(new DateTime(2020, 1, 1, 6, 0, 0), sixHour[1].Index);
            Assert.AreEqual(2d + 3d + 4d + 5d + 6d + 7d, sixHour[1].Value, "06:00 to 11:00 fall in the 06:00 block.");

            var points = new TimeSeries(TimeInterval.Irregular);
            foreach (var (d, v) in new[] { ("2021-01-31", 1.0), ("2021-02-01", 2.0), ("2021-05-15", 3.0), ("2021-12-31", 4.0) })
                points.Add(new SeriesOrdinate<DateTime, double>(DateTime.Parse(d, CultureInfo.InvariantCulture), v));
            var quarters = points.AggregateByInterval(TimeInterval.OneQuarter, BlockFunctionType.Sum);
            CollectionAssert.AreEqual(new[] { new DateTime(2021, 1, 1), new DateTime(2021, 4, 1), new DateTime(2021, 7, 1), new DateTime(2021, 10, 1) }, quarters.IndexesToArray());
            AssertValues(new[] { 3.0, 3.0, double.NaN, 4.0 }, quarters);
            AssertValues(new[] { 10.0 }, points.AggregateByInterval(TimeInterval.OneYear, BlockFunctionType.Sum));
            Assert.AreEqual(12, points.AggregateByInterval(TimeInterval.OneMonth).Count);

            var weekly = points.AggregateByInterval(TimeInterval.SevenDay, BlockFunctionType.Sum);
            Assert.AreEqual(new DateTime(2021, 1, 31), weekly[0].Index, "Seven-day blocks start at the first ordinate's date.");
            Assert.AreEqual(3d, weekly[0].Value);
        }

        /// <summary>
        /// Edge cases: an empty series, invalid arguments, infinities and UTC kind.
        /// </summary>
        [TestMethod]
        public void Test_AggregateByInterval_EdgeCases()
        {
            Assert.AreEqual(0, new TimeSeries().AggregateByInterval(TimeInterval.OneDay).Count);
            var ts = new TimeSeries(TimeInterval.Irregular);
            ts.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2024, 1, 1, 1, 0, 0, DateTimeKind.Utc), double.NegativeInfinity));
            Assert.ThrowsExactly<ArgumentException>(() => ts.AggregateByInterval(TimeInterval.Irregular));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => ts.AggregateByInterval(TimeInterval.OneDay, BlockFunctionType.Average, 0));
            var max = ts.AggregateByInterval(TimeInterval.OneDay, BlockFunctionType.Maximum);
            Assert.AreEqual(double.NegativeInfinity, max[0].Value);
            Assert.AreEqual(DateTimeKind.Utc, max[0].Index.Kind);
        }

        /// <summary>
        /// Asserts value bits equal, NaN included.
        /// </summary>
        private static void AssertValues(double[] expected, TimeSeries actual)
        {
            Assert.HasCount(expected.Length, actual.ValuesToArray());
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i].Value, 0d, $"Block {i}");
        }

        /// <summary>
        /// The ticks accessor returns each ordinate date's ticks in collection order.
        /// </summary>
        [TestMethod]
        public void Test_IndexTicksToArray()
        {
            var ts = new TimeSeries(TimeInterval.OneHour, new DateTime(2020, 1, 1, 5, 0, 0), new double[] { 1, 2, 3 });
            ts.Add(new SeriesOrdinate<DateTime, double>(new DateTime(2019, 12, 31, 0, 0, 0, DateTimeKind.Utc), 4));
            var ticks = ts.IndexTicksToArray();
            var dates = ts.IndexesToArray();
            Assert.HasCount(dates.Length, ticks);
            for (int i = 0; i < dates.Length; i++)
                Assert.AreEqual(dates[i].Ticks, ticks[i]);
            Assert.AreEqual(new DateTime(2020, 1, 1, 6, 0, 0).Ticks, ticks[1]);
            Assert.IsEmpty(new TimeSeries().IndexTicksToArray());
        }

        /// <summary>
        /// Converting to a coarser step returns one ordinate per block that holds data; no empty block follows the data.
        /// </summary>
        /// <remarks>
        /// Expected values match pandas <c>Series.resample('D').mean()</c> and <c>.sum()</c> of hourly values 1, 2, ..., n
        /// starting at midnight.
        /// </remarks>
        [TestMethod]
        public void Test_ConvertTimeInterval_Coarser_NoEmptyTrailingBlock()
        {
            var start = new DateTime(2020, 1, 1);
            var exact = new TimeSeries(TimeInterval.OneHour, start, Enumerable.Range(1, 72).Select(i => (double)i).ToArray());
            var mean = exact.ConvertTimeInterval(TimeInterval.OneDay)!;
            var sum = exact.ConvertTimeInterval(TimeInterval.OneDay, false)!;
            AssertValues(new[] { 12.5, 36.5, 60.5 }, mean);
            AssertValues(new[] { 300.0, 876.0, 1452.0 }, sum);
            Assert.AreEqual(new DateTime(2020, 1, 3), mean.EndDate);

            // A partial last block keeps the steps it has.
            var partial = new TimeSeries(TimeInterval.OneHour, start, Enumerable.Range(1, 50).Select(i => (double)i).ToArray());
            AssertValues(new[] { 12.5, 36.5, 49.5 }, partial.ConvertTimeInterval(TimeInterval.OneDay)!);
            AssertValues(new[] { 300.0, 876.0, 99.0 }, partial.ConvertTimeInterval(TimeInterval.OneDay, false)!);

            // One value into the next day gives a fourth block of one step.
            var single = new TimeSeries(TimeInterval.OneHour, start, Enumerable.Range(1, 73).Select(i => (double)i).ToArray());
            AssertValues(new[] { 12.5, 36.5, 60.5, 73.0 }, single.ConvertTimeInterval(TimeInterval.OneDay)!);

            // A series shorter than one block gives one partial block.
            var shortSeries = new TimeSeries(TimeInterval.OneHour, start, new double[] { 2, 4 });
            AssertValues(new[] { 3.0 }, shortSeries.ConvertTimeInterval(TimeInterval.SixHour)!);
        }

        /// <summary>
        /// By default each difference is labelled with the date of its earlier ordinate.
        /// </summary>
        /// <remarks>
        /// Values are numpy <c>np.diff(x, n=differences)</c> for lag 1 and <c>x[lag:] - x[:-lag]</c> for a seasonal lag;
        /// the labels are the first <c>Count - lag * differences</c> dates of the series.
        /// </remarks>
        [TestMethod]
        public void Test_Difference_LabelsAtEarlierDate()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10 };
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 1, 1), values);

            var first = ts.Difference();
            Assert.AreEqual(7, first.Count);
            Assert.AreEqual(ts.StartDate, first.StartDate);
            Assert.AreEqual(new DateTime(2023, 1, 7), first.EndDate);
            Assert.AreEqual(-6.0, first[0].Value);
            Assert.AreEqual(TimeInterval.OneDay, first.TimeInterval);

            var seasonal = ts.Difference(3);
            AssertValues(new[] { -17.0, -4.0, 3.0, 43.0, -2.0 }, seasonal);
            for (int i = 0; i < seasonal.Count; i++)
                Assert.AreEqual(ts[i].Index, seasonal[i].Index);

            var second = ts.Difference(1, 2);
            AssertValues(new[] { 23.0, -45.0, 35.0, 17.0, -12.0, -50.0 }, second);
            for (int i = 0; i < second.Count; i++)
                Assert.AreEqual(ts[i].Index, second[i].Index);

            // The two-argument call and the explicit earlier label are the same.
            var explicitEarlier = ts.Difference(3, 1, false);
            AssertDates(seasonal.IndexesToArray(), explicitEarlier);
            AssertValues(seasonal.ValuesToArray(), explicitEarlier);
        }

        /// <summary>
        /// A series with gaps, or an irregular one, keeps its own dates as the difference labels.
        /// </summary>
        [TestMethod]
        public void Test_Difference_GappedAndIrregularKeepRealDates()
        {
            var dates = new[] { new DateTime(2020, 1, 1), new DateTime(2020, 1, 2), new DateTime(2020, 1, 5), new DateTime(2020, 1, 6), new DateTime(2020, 1, 20) };
            var values = new double[] { 1, 4, 9, 16, 25 };
            foreach (var interval in new[] { TimeInterval.OneDay, TimeInterval.Irregular })
            {
                var ts = new TimeSeries(interval, dates.Select(d => d.Ticks).ToArray(), values);

                var earlier = ts.Difference();
                AssertDates(dates.Take(4).ToArray(), earlier);
                AssertValues(new double[] { 3, 5, 7, 9 }, earlier);
                Assert.AreEqual(interval, earlier.TimeInterval);

                var later = ts.Difference(1, 1, true);
                AssertDates(dates.Skip(1).ToArray(), later);
                AssertValues(new double[] { 3, 5, 7, 9 }, later);

                var second = ts.Difference(1, 2);
                AssertDates(dates.Take(3).ToArray(), second);
                AssertValues(new double[] { 2, 2, 2 }, second);
            }
        }

        /// <summary>
        /// With the later label, each difference is labelled with the date of its later ordinate.
        /// </summary>
        /// <remarks>
        /// pandas: <c>s.diff(lag)</c> keeps the index of x[t], so after <c>dropna()</c> the first label is the
        /// ordinate at <c>lag * differences</c>; <c>s.diff(1).diff(1)</c> for the second order.
        /// </remarks>
        [TestMethod]
        public void Test_Difference_LabelsAtLaterDate()
        {
            var values = new double[] { 22, 16, 33, 5, 12, 36, 48, 10 };
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2023, 1, 1), values);

            var first = ts.Difference(1, 1, true);
            Assert.AreEqual(7, first.Count);
            Assert.AreEqual(new DateTime(2023, 1, 2), first.StartDate);
            Assert.AreEqual(ts.EndDate, first.EndDate);
            Assert.AreEqual(-6.0, first[0].Value);
            Assert.AreEqual(TimeInterval.OneDay, first.TimeInterval);

            var seasonal = ts.Difference(3, 1, true);
            AssertValues(new[] { -17.0, -4.0, 3.0, 43.0, -2.0 }, seasonal);
            for (int i = 0; i < seasonal.Count; i++)
                Assert.AreEqual(ts[i + 3].Index, seasonal[i].Index);

            var second = ts.Difference(1, 2, true);
            AssertValues(new[] { 23.0, -45.0, 35.0, 17.0, -12.0, -50.0 }, second);
            for (int i = 0; i < second.Count; i++)
                Assert.AreEqual(ts[i + 2].Index, second[i].Index);
        }

        /// <summary>
        /// Difference rejects an empty series and a series no longer than the lag, with either label.
        /// </summary>
        [TestMethod]
        public void Test_Difference_InvalidInputs()
        {
            var shortSeries = new TimeSeries(TimeInterval.OneDay, new DateTime(2020, 1, 1), new double[] { 1, 2 });
            Assert.Throws<ArgumentException>(() => shortSeries.Difference(2));
            Assert.Throws<ArgumentException>(() => shortSeries.Difference(1, 2));
            Assert.Throws<ArgumentException>(() => shortSeries.Difference(2, 1, true));
            Assert.Throws<ArgumentException>(() => new TimeSeries(TimeInterval.OneDay).Difference());
            Assert.Throws<ArgumentException>(() => new TimeSeries(TimeInterval.Irregular).Difference(1, 1, true));
        }

        /// <summary>
        /// Month, quarter and year series chain <see cref="TimeSeries.AddTimeInterval"/> from the previous date, so a
        /// start after day 28 settles on the shortest day of month it meets.
        /// </summary>
        /// <remarks>
        /// Expected dates are .NET <c>DateTime.AddMonths(1)</c> (or <c>AddYears(1)</c>) applied to the previous date.
        /// </remarks>
        [TestMethod]
        public void Test_CalendarSteps_ChainFromPreviousDate()
        {
            var chained = new[] { new DateTime(2020, 1, 31), new DateTime(2020, 2, 29), new DateTime(2020, 3, 29), new DateTime(2020, 4, 29), new DateTime(2020, 5, 29) };
            AssertDates(chained, new TimeSeries(TimeInterval.OneMonth, chained[0], new double[5]));
            AssertDates(chained, new TimeSeries(TimeInterval.OneMonth, chained[0], chained[4]));
            AssertDates(chained, new TimeSeries(TimeInterval.OneMonth, chained[0], chained[4], 2.5));

            var quarters = new[] { new DateTime(2019, 11, 30), new DateTime(2020, 2, 29), new DateTime(2020, 5, 29), new DateTime(2020, 8, 29) };
            AssertDates(quarters, new TimeSeries(TimeInterval.OneQuarter, quarters[0], new double[4]));

            var years = new[] { new DateTime(2020, 2, 29), new DateTime(2021, 2, 28), new DateTime(2022, 2, 28), new DateTime(2023, 2, 28), new DateTime(2024, 2, 28) };
            AssertDates(years, new TimeSeries(TimeInterval.OneYear, years[0], new double[5]));

            // A start on day 28 or earlier keeps its day of month.
            var fifteenths = new TimeSeries(TimeInterval.OneMonth, new DateTime(2023, 1, 15), new double[3]);
            AssertDates(new[] { new DateTime(2023, 1, 15), new DateTime(2023, 2, 15), new DateTime(2023, 3, 15) }, fifteenths);

            var hours = new TimeSeries(TimeInterval.SixHour, new DateTime(2020, 1, 31, 18, 0, 0), new double[3]);
            AssertDates(new[] { new DateTime(2020, 1, 31, 18, 0, 0), new DateTime(2020, 2, 1), new DateTime(2020, 2, 1, 6, 0, 0) }, hours);
        }

        /// <summary>
        /// Shifting, filling and resampling chain <see cref="TimeSeries.AddTimeInterval"/> from the previous date.
        /// </summary>
        [TestMethod]
        public void Test_CalendarSteps_ShiftFillAndResampleChain()
        {
            var chained = new[] { new DateTime(2020, 1, 31), new DateTime(2020, 2, 29), new DateTime(2020, 3, 29), new DateTime(2020, 4, 29), new DateTime(2020, 5, 29) };

            var shifted = new TimeSeries(TimeInterval.OneMonth, new DateTime(2019, 1, 1), new double[] { 1, 2, 3, 4, 5 });
            shifted.ShiftAllDates(chained[0]);
            AssertDates(chained, shifted);
            AssertValues(new double[] { 1, 2, 3, 4, 5 }, shifted);

            var sparse = new TimeSeries(TimeInterval.OneMonth);
            sparse.Add(new SeriesOrdinate<DateTime, double>(chained[2], 7));
            var filled = TimeSeries.FillMissingDates(sparse, chained[0], chained[4], -1);
            AssertDates(chained, filled);
            AssertValues(new double[] { -1, -1, 7, -1, -1 }, filled);

            var source = new TimeSeries(TimeInterval.OneMonth, chained[0], Enumerable.Range(0, 24).Select(i => Math.Sin(i)).ToArray());
            AssertDates(chained, source.ResampleWithBlockBootstrap(5, 3, 7));
            AssertDates(chained, source.ResampleWithKNN(5, 3, 7));
        }

        /// <summary>
        /// The start and end date constructors stop at the last stepped date on or before the end date.
        /// </summary>
        [TestMethod]
        public void Test_StartEndConstructor_NoDateAfterEnd()
        {
            // Chained monthly dates from January 31 fall on the 29th, so June 30 is between two of them.
            var chained = new[] { new DateTime(2020, 1, 31), new DateTime(2020, 2, 29), new DateTime(2020, 3, 29), new DateTime(2020, 4, 29), new DateTime(2020, 5, 29), new DateTime(2020, 6, 29) };
            var monthly = new TimeSeries(TimeInterval.OneMonth, chained[0], new DateTime(2020, 6, 30));
            AssertDates(chained, monthly);
            Assert.IsTrue(monthly.ValuesToArray().All(double.IsNaN));
            var fixedMonthly = new TimeSeries(TimeInterval.OneMonth, chained[0], new DateTime(2020, 6, 30), 4.5);
            AssertDates(chained, fixedMonthly);
            AssertValues(Enumerable.Repeat(4.5, 6).ToArray(), fixedMonthly);

            // An end between two fixed steps.
            var hourly = new TimeSeries(TimeInterval.SixHour, new DateTime(2020, 1, 1), new DateTime(2020, 1, 1, 13, 0, 0));
            AssertDates(new[] { new DateTime(2020, 1, 1), new DateTime(2020, 1, 1, 6, 0, 0), new DateTime(2020, 1, 1, 12, 0, 0) }, hourly);

            // An end on a stepped date is the last ordinate, and equal start and end dates give one ordinate.
            AssertDates(chained.Take(5).ToArray(), new TimeSeries(TimeInterval.OneMonth, chained[0], chained[4]));
            AssertDates(new[] { chained[0] }, new TimeSeries(TimeInterval.OneMonth, chained[0], chained[0], 1.0));
            AssertDates(new[] { chained[0] }, new TimeSeries(TimeInterval.OneYear, chained[0], new DateTime(2020, 12, 31)));

            Assert.Throws<ArgumentException>(() => new TimeSeries(TimeInterval.OneDay, chained[1], chained[0]));
            Assert.Throws<ArgumentException>(() => new TimeSeries(TimeInterval.Irregular, chained[0], chained[1], 0.0));
        }

        /// <summary>
        /// An end date found by chaining <see cref="TimeSeries.AddTimeInterval"/> n - 1 times from the start gives
        /// exactly n ordinates, for every regular interval.
        /// </summary>
        [TestMethod]
        public void Test_StartEndConstructor_ChainedEndGivesRequestedCount()
        {
            var starts = new[] { new DateTime(2000, 1, 1), new DateTime(2020, 1, 31), new DateTime(2020, 2, 29), new DateTime(2019, 8, 30), new DateTime(2021, 12, 31, 23, 45, 0) };
            foreach (TimeInterval interval in Enum.GetValues(typeof(TimeInterval)))
            {
                if (interval == TimeInterval.Irregular) continue;
                foreach (var start in starts)
                {
                    var end = start;
                    for (int n = 1; n <= 30; n++)
                    {
                        var ts = new TimeSeries(interval, start, end);
                        Assert.AreEqual(n, ts.Count, $"{interval} from {start:O}");
                        Assert.AreEqual(end, ts.EndDate, $"{interval} from {start:O}");
                        end = TimeSeries.AddTimeInterval(end, interval);
                    }
                }
            }
        }

        /// <summary>
        /// <see cref="TimeSeries.AnchoredDates"/> counts month, quarter and year steps from the start date.
        /// </summary>
        /// <remarks>
        /// Expected dates are pandas <c>Timestamp(start) + DateOffset(months=k)</c> (or <c>years=k</c>) for k = 0, 1, ...
        /// </remarks>
        [TestMethod]
        public void Test_AnchoredDates_KeepDayOfMonth()
        {
            var monthEnds = new[] { new DateTime(2020, 1, 31), new DateTime(2020, 2, 29), new DateTime(2020, 3, 31), new DateTime(2020, 4, 30), new DateTime(2020, 5, 31) };
            CollectionAssert.AreEqual(monthEnds, TimeSeries.AnchoredDates(monthEnds[0], 5, TimeInterval.OneMonth));

            var quarters = new[] { new DateTime(2019, 11, 30), new DateTime(2020, 2, 29), new DateTime(2020, 5, 30), new DateTime(2020, 8, 30) };
            CollectionAssert.AreEqual(quarters, TimeSeries.AnchoredDates(quarters[0], 4, TimeInterval.OneQuarter));

            var leapYears = new[] { new DateTime(2020, 2, 29), new DateTime(2021, 2, 28), new DateTime(2022, 2, 28), new DateTime(2023, 2, 28), new DateTime(2024, 2, 29) };
            CollectionAssert.AreEqual(leapYears, TimeSeries.AnchoredDates(leapYears[0], 5, TimeInterval.OneYear));

            // A series on the anchored dates, through the ticks constructor.
            var ticks = TimeSeries.AnchoredDates(monthEnds[0], 5, TimeInterval.OneMonth).Select(d => d.Ticks).ToArray();
            var series = new TimeSeries(TimeInterval.OneMonth, ticks, new double[] { 1, 2, 3, 4, 5 });
            AssertDates(monthEnds, series);
        }

        /// <summary>
        /// <see cref="TimeSeries.AnchoredDates"/> equals the chained series dates for fixed intervals and for a start on
        /// day 28 or earlier.
        /// </summary>
        [TestMethod]
        public void Test_AnchoredDates_MatchChainedWhereTheyAgree()
        {
            var fixedIntervals = new[] { TimeInterval.OneMinute, TimeInterval.FiveMinute, TimeInterval.FifteenMinute, TimeInterval.ThirtyMinute, TimeInterval.OneHour, TimeInterval.SixHour, TimeInterval.TwelveHour, TimeInterval.OneDay, TimeInterval.SevenDay };
            var start = new DateTime(2020, 1, 31, 18, 30, 0);
            foreach (var interval in fixedIntervals)
                CollectionAssert.AreEqual(new TimeSeries(interval, start, new double[40]).IndexesToArray(), TimeSeries.AnchoredDates(start, 40, interval), interval.ToString());

            foreach (var interval in new[] { TimeInterval.OneMonth, TimeInterval.OneQuarter, TimeInterval.OneYear })
            {
                for (int day = 1; day <= 28; day++)
                {
                    var first = new DateTime(2019, 1, day);
                    CollectionAssert.AreEqual(new TimeSeries(interval, first, new double[30]).IndexesToArray(), TimeSeries.AnchoredDates(first, 30, interval), $"{interval} day {day}");
                }
            }
        }

        /// <summary>
        /// <see cref="TimeSeries.AnchoredDates"/> handles an empty and a single-date sequence and rejects invalid input.
        /// </summary>
        [TestMethod]
        public void Test_AnchoredDates_EdgeCasesAndInvalidInputs()
        {
            var start = new DateTime(2020, 1, 31);
            Assert.IsEmpty(TimeSeries.AnchoredDates(start, 0, TimeInterval.OneMonth));
            CollectionAssert.AreEqual(new[] { start }, TimeSeries.AnchoredDates(start, 1, TimeInterval.OneYear));
            Assert.Throws<ArgumentException>(() => TimeSeries.AnchoredDates(start, 3, TimeInterval.Irregular));
            Assert.Throws<ArgumentOutOfRangeException>(() => TimeSeries.AnchoredDates(start, -1, TimeInterval.OneMonth));
            Assert.Throws<ArgumentOutOfRangeException>(() => TimeSeries.AnchoredDates(DateTime.MaxValue.AddMonths(-1), 3, TimeInterval.OneMonth));
        }

        /// <summary>
        /// Asserts the ordinate dates equal, in order.
        /// </summary>
        private static void AssertDates(DateTime[] expected, TimeSeries actual)
        {
            Assert.HasCount(expected.Length, actual.IndexesToArray());
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i].Index, $"Ordinate {i}");
        }

        /// <summary>
        /// The ticks constructor round-trips <see cref="TimeSeries.IndexTicksToArray"/> and the values bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_TicksConstructor_RoundTrip()
        {
            var nanPayload = BitConverter.Int64BitsToDouble(unchecked((long)0xFFF8_0000_0000_1234UL));
            var ticks = new long[] { new DateTime(2020, 2, 29, 23, 59, 59).Ticks + 7, new DateTime(1850, 6, 1).Ticks, DateTime.MinValue.Ticks, DateTime.MaxValue.Ticks };
            var values = new double[] { nanPayload, -0.0, double.PositiveInfinity, 1.5 };

            var ts = new TimeSeries(TimeInterval.Irregular, ticks, values);

            Assert.AreEqual(TimeInterval.Irregular, ts.TimeInterval);
            CollectionAssert.AreEqual(ticks, ts.IndexTicksToArray());
            var stored = ts.ValuesToArray();
            for (int i = 0; i < values.Length; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(values[i]), BitConverter.DoubleToInt64Bits(stored[i]));
            Assert.AreEqual(DateTimeKind.Unspecified, ts[0].Index.Kind);

            // The same dates and values as the start-date constructor.
            var regular = new TimeSeries(TimeInterval.OneMonth, new DateTime(2020, 1, 31), new double[] { 1, 2, 3 });
            var copy = new TimeSeries(TimeInterval.OneMonth, regular.IndexTicksToArray(), regular.ValuesToArray());
            AssertDates(regular.IndexesToArray(), copy);
            AssertValues(regular.ValuesToArray(), copy);
            Assert.AreEqual(0, new TimeSeries(TimeInterval.OneDay, new long[0], new double[0]).Count);
        }

        /// <summary>
        /// The ticks constructor rejects null, mismatched and out-of-range input.
        /// </summary>
        [TestMethod]
        public void Test_TicksConstructor_InvalidInputs()
        {
            Assert.Throws<ArgumentNullException>(() => new TimeSeries(TimeInterval.OneDay, (IList<long>)null!, new double[0]));
            Assert.Throws<ArgumentNullException>(() => new TimeSeries(TimeInterval.OneDay, new long[0], null!));
            Assert.Throws<ArgumentException>(() => new TimeSeries(TimeInterval.OneDay, new long[] { 0, 1 }, new double[] { 1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new TimeSeries(TimeInterval.OneDay, new long[] { -1 }, new double[] { 1 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new TimeSeries(TimeInterval.OneDay, new long[] { DateTime.MaxValue.Ticks + 1 }, new double[] { 1 }));
        }

    }
}
