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

            // now the time interval is a week, so we will get different maxes
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

    }
}
