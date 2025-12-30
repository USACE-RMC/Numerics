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

using Numerics.Data.Statistics;
using Numerics.MachineLearning;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Globalization;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Linq;

namespace Numerics.Data
{
    /// <summary>
    /// A time-series class, which is a collection of time-series ordinates.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class TimeSeries : Series<DateTime, double>
    {
        /// <summary>
        /// Constructs an empty time-series.
        /// </summary>
        public TimeSeries() { }

        /// <summary>
        /// Constructs an empty time-series with a specified time interval.
        /// </summary>
        /// <param name="timeInterval">The time interval for the series.</param>
        public TimeSeries(TimeInterval timeInterval)
        {
            _timeInterval = timeInterval;
        }

        /// <summary>
        /// Constructs an empty time-series with a specified start and end date.
        /// </summary>
        /// <param name="timeInterval">The time interval for the series.</param>
        /// <param name="startDate">The start date/time of the series.</param>
        /// <param name="endDate">The end date/time of the series.</param>
        public TimeSeries(TimeInterval timeInterval, DateTime startDate, DateTime endDate)
        {
            if (timeInterval == TimeInterval.Irregular)
                throw new ArgumentException("The time interval cannot be irregular with this constructor.");

            if (startDate > endDate)
                throw new ArgumentException("Start date must be less than or equal to end date.");

            _timeInterval = timeInterval;
            Add(new SeriesOrdinate<DateTime, double>(startDate, double.NaN));
            while (_seriesOrdinates.Last().Index < endDate)
            {
                Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(_seriesOrdinates.Last().Index, _timeInterval), double.NaN));
            }
        }

        /// <summary>
        /// Constructs a time-series with a specified start and end date, and a constant fixed value.
        /// </summary>
        /// <param name="timeInterval">The time interval for the series.</param>
        /// <param name="startDate">The start date/time of the series.</param>
        /// <param name="endDate">The end date/time of the series.</param>
        /// <param name="fixedValue">A fixed value to be assigned to each ordinate.</param>
        public TimeSeries(TimeInterval timeInterval, DateTime startDate, DateTime endDate, double fixedValue)
        {
            if (timeInterval == TimeInterval.Irregular)
                throw new ArgumentException("The time interval cannot be irregular with this constructor.");

            if (startDate > endDate)
                throw new ArgumentException("Start date must be less than or equal to end date.");

            _timeInterval = timeInterval;
            Add(new SeriesOrdinate<DateTime, double>(startDate, fixedValue));
            while (_seriesOrdinates.Last().Index < endDate)
            {
                Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(_seriesOrdinates.Last().Index, _timeInterval), fixedValue));
            }
        }

        /// <summary>
        /// Constructs a time-series based on the start date and a list of data values.
        /// </summary>
        /// <param name="timeInterval">The time interval for the series.</param>
        /// <param name="startDate">The start date/time of the series.</param>
        /// <param name="data">A list of data values.</param>
        public TimeSeries(TimeInterval timeInterval, DateTime startDate, IList<double> data)
        {
            if (timeInterval == TimeInterval.Irregular)
                throw new ArgumentException("The time interval cannot be irregular with this constructor.");

            _timeInterval = timeInterval;
            if (data.Count == 0) return;
            Add(new SeriesOrdinate<DateTime, double>(startDate, data[0]));
            for (int i = 1; i < data.Count; i++)
            {
                Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(this[i - 1].Index, _timeInterval), data[i]));
            }
        }

        /// <summary>
        /// Constructs a time-series based on XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        public TimeSeries(XElement xElement)
        {
            // Get time interval
            var timeIntervalAttr = xElement.Attribute(nameof(TimeInterval));
            if (timeIntervalAttr != null)
                Enum.TryParse(timeIntervalAttr.Value, out _timeInterval);

            // Get Ordinates
            foreach (XElement ordinate in xElement.Elements("SeriesOrdinate"))
            {
                // Try to parse the invariant date string using TryParseExact
                // If it fails, do a regular try parse.
                DateTime index = default;
                var ordAttr = ordinate.Attribute("Index");
                if (ordAttr != null)
                {
                    if (!DateTime.TryParseExact(ordAttr.Value, "o", CultureInfo.InvariantCulture, DateTimeStyles.RoundtripKind, out index))
                    {
                        DateTime.TryParse(ordAttr.Value, out index);
                    }
                }

                double value = 0.0;
                var ordVal = ordinate.Attribute("Value");
                if (ordVal != null)
                {
                    double.TryParse(ordVal.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out value);
                }
                
                Add(new SeriesOrdinate<DateTime, double>(index, value));
            }
        }

        private TimeInterval _timeInterval = TimeInterval.OneDay;

        /// <summary>
        /// Returns the time interval of the time-series.
        /// </summary>
        public TimeInterval TimeInterval
        {
            get { return _timeInterval; }
        }

        /// <summary>
        /// Gets whether there are missing values.
        /// </summary>
        public bool HasMissingValues
        {
            get { return NumberOfMissingValues() > 0 ? true : false; }
        }

        /// <summary>
        /// Gets the start date of the time-series.
        /// </summary>
        public DateTime StartDate
        {
            get { return _seriesOrdinates.Min(x => x.Index); }
        }

        /// <summary>
        /// Gets the end date of the time-series.
        /// </summary>
        public DateTime EndDate
        {
            get { return _seriesOrdinates.Max(x => x.Index); }
        }

        /// <summary>
        /// Sorts the elements in the entire collection by the time ordinate given a specified sort direction.
        /// </summary>
        /// <param name="Order">Optional. Ascending or descending order. Default = Ascending.</param>
        public void SortByTime(ListSortDirection Order = ListSortDirection.Ascending)
        {
            if (Order == ListSortDirection.Ascending)
            {
                _seriesOrdinates.Sort((x, y) => x.Index.CompareTo(y.Index));
            }
            else
            {
                _seriesOrdinates.Sort((x, y) => -1 * x.Index.CompareTo(y.Index));
            }
        }

        /// <summary>
        /// Sorts the elements in the entire collection by value given a specified sort direction.
        /// </summary>
        /// <param name="Order">Optional. Ascending or descending order. Default = Ascending.</param>
        public void SortByValue(ListSortDirection Order = ListSortDirection.Ascending)
        {
            if (Order == ListSortDirection.Ascending)
            {
                _seriesOrdinates.Sort((x, y) => x.Value.CompareTo(y.Value));
            }
            else
            {
                _seriesOrdinates.Sort((x, y) => -1 * x.Value.CompareTo(y.Value));
            }
        }

        /// <summary>
        /// Add a constant to each value in the time-series. Missing values are kept as missing.
        /// </summary>
        /// <param name="constant">Factor to add to each value in the series.</param>
        public void Add(double constant)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i <= Count - 1; i++)
            {
                if (!double.IsNaN(this[i].Value)) { this[i].Value += constant; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Add a constant to specified values in the time-series. Missing values are kept as missing.
        /// </summary>
        /// <param name="constant">Factor to add to each value in the series.</param>
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the calculation to.</param>
        public void Add(double constant, IList<int> indexes)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && (!double.IsNaN(this[indexes[i]].Value))) { this[indexes[i]].Value += constant; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Subtract a constant from each value in the time-series. Missing values are kept as missing.
        /// </summary>
        /// <param name="constant">Factor to subtract each value in the series.</param>
        public void Subtract(double constant)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i <= Count - 1; i++)
            {
                if (!double.IsNaN(this[i].Value)) { this[i].Value -= constant; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Subtract a constant from specified values in the time-series. Missing values are kept as missing.
        /// </summary>
        /// <param name="constant">Factor to subtract each value in the series.</param>
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the calculation to.</param>
        public void Subtract(double constant, IList<int> indexes)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && (!double.IsNaN(this[indexes[i]].Value))) { this[indexes[i]].Value -= constant; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Multiply each value in the time-series by a constant. Missing values are kept as missing.
        /// </summary>
        /// <param name="constant">Factor to multiply each value by in the series.</param>
        public void Multiply(double constant)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i <= Count - 1; i++)
            {
                if (!double.IsNaN(this[i].Value)) { this[i].Value *= constant; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Multiply specified values in the time-series by a constant. Missing values are kept as missing.
        /// </summary>
        /// <param name="constant">Factor to multiply each value by in the series.</param>
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the calculation to.</param>
        public void Multiply(double constant, IList<int> indexes)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && (!double.IsNaN(this[indexes[i]].Value))) { this[indexes[i]].Value *= constant; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Divide each value in the time-series by a constant. Missing values are kept as missing.
        /// </summary>
        /// <param name="constant">Factor to divide each value by in the series.</param>
        public void Divide(double constant)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i <= Count - 1; i++)
            {
                if (!double.IsNaN(this[i].Value)) { this[i].Value /= constant; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Divide specified values in the time-series by a constant. Missing values are kept as missing.
        /// </summary>
        /// <param name="constant">Factor to divide each value by in the series.</param>
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the calculation to.</param>
        public void Divide(double constant, IList<int> indexes)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && (!double.IsNaN(this[indexes[i]].Value))) { this[indexes[i]].Value /= constant; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Set each value in the time-series to its absolute value. Missing values are kept as missing.
        /// </summary>
        public void AbsoluteValue()
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i <= Count - 1; i++)
            {
                if (!double.IsNaN(this[i].Value)) { this[i].Value = Math.Abs(this[i].Value); }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Set specified values in the time-series to its absolute value. Missing values are kept as missing.
        /// </summary>
        public void AbsoluteValue(IList<int> indexes)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && (!double.IsNaN(this[indexes[i]].Value))) { this[indexes[i]].Value = Math.Abs(this[indexes[i]].Value); }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Raise each value in the time-series by the specified power or exponent. Missing values are kept as missing.
        /// </summary>
        /// <param name="power">Power or exponent.</param>
        public void Exponentiate(double power)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i <= Count - 1; i++)
                if (!double.IsNaN(this[i].Value))
                    this[i].Value = Math.Pow(this[i].Value, power);
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Raise specified values in the time-series by the specified power or exponent. Missing values are kept as missing.
        /// </summary>
        /// <param name="power">Power or exponent.</param>
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the calculation to.</param>
        public void Exponentiate(double power, IList<int> indexes)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && (!double.IsNaN(this[indexes[i]].Value))) { this[indexes[i]].Value = Math.Pow(this[indexes[i]].Value, power); }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Log transform values in the time-series. Missing values are kept as missing.
        /// </summary>
        /// <param name="baseValue">The log base value.</param>
        public void LogTransform(double baseValue = 10)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < Count; i++)
            {
                if (this[i].Value > 0 && !double.IsNaN(this[i].Value))
                    this[i].Value = Math.Log(this[i].Value, baseValue);
                else
                    this[i].Value = double.NaN;
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Log transform specified values in the time-series. Missing values are kept as missing.
        /// </summary>
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the calculation to.</param>
        /// <param name="baseValue">The log base value.</param>
        public void LogTransform(IList<int> indexes, double baseValue = 10)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && this[indexes[i]].Value > 0 && (!double.IsNaN(this[indexes[i]].Value))) { this[indexes[i]].Value = Math.Log(this[indexes[i]].Value, baseValue); }
                else { this[indexes[i]].Value = double.NaN; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Standardize the time series values. This action is not reversible. 
        /// </summary>
        public void Standardize()
        {
            double mean = MeanValue();
            double stdDev = StandardDeviation();

            if (stdDev == 0)
                throw new InvalidOperationException("Standard deviation is zero, standardization cannot be performed.");

            SuppressCollectionChanged = true;
            for (int i = 0; i < Count; i++)
            {
                this[i].Value = (this[i].Value - mean) / stdDev;
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Each value in the time-series is replaced by its inverse (1/x). Missing values are kept as missing. If the value is 0.0, the value is set to Double.NaN.
        /// </summary>
        public void Inverse()
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < Count; i++)
            {
                if (this[i].Value != 0 && !double.IsNaN(this[i].Value)) { this[i].Value = 1d / this[i].Value; }
                else if (this[i].Value == 0 || double.IsNaN(this[i].Value)) { this[i].Value = double.NaN; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Specified values in the time-series are replaced by their inverse (1/x). Missing values are kept as missing. If the value is 0.0, the value is set to Double.NaN.
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the inverse calculation to.</param>
        /// </summary>
        public void Inverse(IList<int> indexes)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && this[indexes[i]].Value != 0 && !double.IsNaN(this[indexes[i]].Value)) { this[indexes[i]].Value = 1d / this[indexes[i]].Value; }
                else if (this[indexes[i]].Value == 0 || double.IsNaN(this[indexes[i]].Value)) { this[indexes[i]].Value = double.NaN; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Returns the cumulative sum of the time-series. Missing values are treated as zero when accumulating values.
        /// </summary>
        public TimeSeries CumulativeSum()
        {
            var timeSeries = new TimeSeries();
            double sum = 0d;
            for (int i = 0; i < Count; i++)
            {
                if (this[i].Value != default && double.IsNaN(this[i].Value) == false)
                    sum += this[i].Value;
                timeSeries.Add(this[i].Clone());
                timeSeries.Last().Value = sum;
            }
            return timeSeries;
        }

        /// <summary>
        /// Returns a time-series of the successive differences.
        /// </summary>
        /// <param name="lag">The lag between elements to subtract (default is 1).</param>
        /// <param name="differences">The order of differences (default is 1).</param>
        /// <returns>A time series of differences.</returns>
        public TimeSeries Difference(int lag = 1, int differences = 1)
        {
            double[] result = ValuesToArray();
            for (int d = 0; d < differences; d++)
            {
                if (result.Length <= lag)
                    throw new ArgumentException("The length of the array must be greater than the lag.");

                // Create a new array that is shorter by lag.
                double[] temp = new double[result.Length - lag];
                for (int i = 0; i < result.Length - lag; i++)
                {
                    temp[i] = result[i + lag] - result[i];
                }
                result = temp;
            }
            // Adjust the start date if needed, depending on how TimeSeries is defined.
            var timeSeries = new TimeSeries(TimeInterval, StartDate, result);
            return timeSeries;
        }


        /// <summary>
        /// Returns the number of missing values.
        /// </summary>
        public int NumberOfMissingValues()
        {
            return _seriesOrdinates.Where(x => double.IsNaN(x.Value)).Count();
        }

        /// <summary>
        /// Replaces all missing data (Double.NaN) with the specified value.
        /// </summary>
        /// <param name="value">Value for missing data.</param>
        public void ReplaceMissingData(double value)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < Count; i++)
            {
                if (double.IsNaN(this[i].Value)) { this[i].Value = value; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Replaces missing data (Double.NaN) for specified indexes with the specified value.
        /// </summary>
        /// <param name="value">Value for missing data.</param>
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the calculation to.</param>
        public void ReplaceMissingData(IList<int> indexes, double value)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i < indexes.Count; i++)
            {
                if (indexes[i] >= 0 && indexes[i] < Count && double.IsNaN(this[indexes[i]].Value)) { this[indexes[i]].Value = value; }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Interpolate missing data. Data will only be interpolated if the number of consecutive missing value is less than the specified limit.
        /// </summary>
        /// <param name="maxNumberOfMissing">The maximum number of consecutive missing values.</param>
        public void InterpolateMissingData(int maxNumberOfMissing)
        {
            SuppressCollectionChanged = true;
            SortByTime();
            double x;
            //double y;
            double x1;
            double x2;
            double y1;
            double y2;
            int upper;
            // 
            for (int i = 1; i < Count; i++)
            {
                // Find missing value
                if (double.IsNaN(this[i].Value))
                {
                    // ok we found one
                    x = this[i].Index.ToOADate();
                    x1 = this[i - 1].Index.ToOADate();
                    y1 = this[i - 1].Value;
                    upper = i + maxNumberOfMissing;
                    // Find the next non-missing value
                    for (int j = i; j <= Math.Min(Count - 1, upper); j++)
                    {
                        // ok we found one
                        // the interpolation case
                        if (!double.IsNaN(this[j].Value))
                        {
                            x2 = this[j].Index.ToOADate();
                            y2 = this[j].Value;
                            this[i].Value = y1 + (x - x1) / (x2 - x1) * (y2 - y1);
                            break;
                        }
                        // the extrapolation case
                        if (j == Count - 1)
                        {
                            x1 = this[i - 2].Index.ToOADate();
                            x2 = this[i - 1].Index.ToOADate();
                            y1 = this[i - 2].Value;
                            y2 = this[i - 1].Value;
                            this[i].Value = y1 - (x1 - x) * (y2 - y1) / (x2 - x1);
                        }
                    }
                }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Interpolate missing data. Data will only be interpolated if the number of consecutive missing value is less than the specified limit.
        /// </summary>
        /// <param name="maxNumberOfMissing">The maximum number of consecutive missing values.</param>
        /// <param name="indexes">List of integer index values (0 based) for each ordinate in the time series to apply the calculation to.</param>

        public void InterpolateMissingData(int maxNumberOfMissing, IList<int> indexes)
        {
            SuppressCollectionChanged = true;
            SortByTime();
            double x;
            //double y;
            double x1;
            double x2;
            double y1;
            double y2;
            int upper;
            // 
            for (int i = 0; i < indexes.Count; i++)
            {
                // Find missing value
                if(indexes[i] >= 1 && indexes[i] < Count && double.IsNaN(this[indexes[i]].Value))
                {
                    // ok we found one
                    int idx = indexes[i];
                    x = this[idx].Index.ToOADate();
                    x1 = this[idx - 1].Index.ToOADate();
                    y1 = this[idx - 1].Value;
                    upper = idx + maxNumberOfMissing;
                    // Find the next non-missing value
                    for (int j = idx; j <= Math.Min(Count - 1, upper); j++)
                    {
                        // ok we found one
                        // the interpolation case
                        if (!double.IsNaN(this[j].Value))
                        {
                            x2 = this[j].Index.ToOADate();
                            y2 = this[j].Value;
                            this[idx].Value = y1 + (x - x1) / (x2 - x1) * (y2 - y1);
                            break;
                        }
                        // the extrapolation case
                        if (j == Count - 1)
                        {
                            x1 = this[idx - 2].Index.ToOADate();
                            x2 = this[idx - 1].Index.ToOADate();
                            y1 = this[idx - 2].Value;
                            y2 = this[idx - 1].Value;
                            this[idx].Value = y1 - (x1 - x) * (y2 - y1) / (x2 - x1);
                        }
                    }
                }
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        public static TimeSeries FillMissingDates(TimeSeries timeSeries, DateTime startDate, DateTime endDate, double value = 0.0)
        {
            if (timeSeries.TimeInterval == TimeInterval.Irregular) throw new Exception("This method does not work with irregular data.");

            // Create a dictionary for fast lookup
            // var lookup = this.ToList().ToDictionary(p => p.Index, p => p.Value);

            var lookup = timeSeries.IndexesToList();

            var result = new TimeSeries(timeSeries.TimeInterval);

            // Loop over and add values
            for (DateTime date = startDate; date <= endDate; date = AddTimeInterval(date, timeSeries.TimeInterval))
            {
                var idx = lookup.IndexOf(date);
                if (idx == -1)
                {
                    result.Add(new SeriesOrdinate<DateTime, double>(date, value));
                }
                else
                {
                    result.Add(new SeriesOrdinate<DateTime, double>(date, timeSeries[idx].Value));
                }
            }

            return result;

        }

        /// <summary>
        /// Returns a new date/time that adds the time interval to the specified data/time value.
        /// </summary>
        /// <param name="time">Time to increase.</param>
        /// <param name="timeInterval">The time interval.</param>
        public static DateTime AddTimeInterval(DateTime time, TimeInterval timeInterval)
        {
            switch (timeInterval)
            {
                case TimeInterval.OneMinute:
                    {
                        return time.AddMinutes(1d);
                    }

                case TimeInterval.FiveMinute:
                    {
                        return time.AddMinutes(5d);
                    }

                case TimeInterval.FifteenMinute:
                    {
                        return time.AddMinutes(15d);
                    }

                case TimeInterval.ThirtyMinute:
                    {
                        return time.AddMinutes(30d);
                    }

                case TimeInterval.OneHour:
                    {
                        return time.AddHours(1d);
                    }

                case TimeInterval.SixHour:
                    {
                        return time.AddHours(6d);
                    }

                case TimeInterval.TwelveHour:
                    {
                        return time.AddHours(12d);
                    }

                case TimeInterval.OneDay:
                    {
                        return time.AddDays(1d);
                    }

                case TimeInterval.SevenDay:
                    {
                        return time.AddDays(7d);
                    }

                case TimeInterval.OneMonth:
                    {
                        return time.AddMonths(1);
                    }

                case TimeInterval.OneQuarter:
                    {
                        return time.AddMonths(3);
                    }

                case TimeInterval.OneYear:
                    {
                        return time.AddYears(1);
                    }
            }

            return time;
        }

        /// <summary>
        /// Returns a new date/time that subtracts the time interval to the specified data/time value.
        /// </summary>
        /// <param name="time">Time to decrease.</param>
        /// <param name="timeInterval">The time interval.</param>
        public static DateTime SubtractTimeInterval(DateTime time, TimeInterval timeInterval)
        {
            switch (timeInterval)
            {
                case TimeInterval.OneMinute:
                    {
                        return time.AddMinutes(-1d);
                    }

                case TimeInterval.FiveMinute:
                    {
                        return time.AddMinutes(-5d);
                    }

                case TimeInterval.FifteenMinute:
                    {
                        return time.AddMinutes(-15d);
                    }

                case TimeInterval.ThirtyMinute:
                    {
                        return time.AddMinutes(-30d);
                    }

                case TimeInterval.OneHour:
                    {
                        return time.AddHours(-1d);
                    }

                case TimeInterval.SixHour:
                    {
                        return time.AddHours(-6d);
                    }

                case TimeInterval.TwelveHour:
                    {
                        return time.AddHours(-12d);
                    }

                case TimeInterval.OneDay:
                    {
                        return time.AddDays(-1d);
                    }

                case TimeInterval.SevenDay:
                    {
                        return time.AddDays(-7d);
                    }

                case TimeInterval.OneMonth:
                    {
                        return time.AddMonths(-1);
                    }

                case TimeInterval.OneQuarter:
                    {
                        return time.AddMonths(-3);
                    }

                case TimeInterval.OneYear:
                    {
                        return time.AddYears(-1);
                    }
            }

            return time;
        }

        /// <summary>
        /// Converts the time interval to hours.
        /// </summary>
        /// <returns>The time interval in hours.</returns>
        public static double TimeIntervalInHours(TimeInterval timeInterval)
        {
            if (timeInterval == TimeInterval.OneMinute) return 1d / 60d;
            if (timeInterval == TimeInterval.FiveMinute) return 5d / 60d;
            if (timeInterval == TimeInterval.FifteenMinute) return 15d / 60d;
            if (timeInterval == TimeInterval.ThirtyMinute) return 30d / 60d;
            if (timeInterval == TimeInterval.OneHour) return 1d;
            if (timeInterval == TimeInterval.SixHour) return 6d;
            if (timeInterval == TimeInterval.TwelveHour) return 12d;
            if (timeInterval == TimeInterval.OneDay) return 24d;
            if (timeInterval == TimeInterval.SevenDay) return 24d * 7d;
            return double.NaN;
        }

        /// <summary>
        /// Determines if the minimum step between events has been exceeded. 
        /// </summary>
        /// <param name="startTime">Start time of the starting event.</param>
        /// <param name="endTime">End time of the ending event.</param>
        /// <param name="minStepsBetweenEvents">Minimum time steps between events.</param>
        private bool CheckIfMinStepsExceeded(DateTime startTime, DateTime endTime, int minStepsBetweenEvents)
        {
            switch (TimeInterval)
            {
                case TimeInterval.OneMinute:
                    {
                        return endTime > startTime.AddMinutes(1 * minStepsBetweenEvents);
                    }

                case TimeInterval.FiveMinute:
                    {
                        return endTime > startTime.AddMinutes(5 * minStepsBetweenEvents);
                    }

                case TimeInterval.FifteenMinute:
                    {
                        return endTime > startTime.AddMinutes(15 * minStepsBetweenEvents);
                    }

                case TimeInterval.ThirtyMinute:
                    {
                        return endTime > startTime.AddMinutes(30 * minStepsBetweenEvents);
                    }

                case TimeInterval.OneHour:
                    {
                        return endTime > startTime.AddHours(1 * minStepsBetweenEvents);
                    }

                case TimeInterval.SixHour:
                    {
                        return endTime > startTime.AddHours(6 * minStepsBetweenEvents);
                    }

                case TimeInterval.TwelveHour:
                    {
                        return endTime > startTime.AddHours(12 * minStepsBetweenEvents);
                    }

                case TimeInterval.OneDay:
                    {
                        return endTime > startTime.AddDays(minStepsBetweenEvents);
                    }

                case TimeInterval.SevenDay:
                    {
                        return endTime > startTime.AddDays(7 * minStepsBetweenEvents);
                    }

                case TimeInterval.OneMonth:
                    {
                        return endTime > startTime.AddMonths(1 * minStepsBetweenEvents);
                    }

                case TimeInterval.OneQuarter:
                    {
                        return endTime > startTime.AddYears(1 * minStepsBetweenEvents);
                    }

                case TimeInterval.OneYear:
                    {
                        return endTime > startTime.AddYears(1 * minStepsBetweenEvents);
                    }
            }

            return false;
        }

        /// <summary>
        /// Returns a moving average time-series based on the specified time period. The average is computed based on the previous n=period ordinates.
        /// </summary>
        /// <param name="period">The time period to average over. If time interval is 1-hour, and period is 12, the moving average will be computed over a moving 12 hour block.</param>
        public TimeSeries MovingAverage(int period)
        {
            if (period >= Count)
                throw new ArgumentException(nameof(period), "The period must be less than the length of the time-series.");
            SortByTime();
            var timeSeries = new TimeSeries(TimeInterval);
            double sum = 0d;
            double avg = 0d;
            for (int i = 1; i <= Count; i++)
            {
                sum += !double.IsNaN(this[i - 1].Value) ? this[i - 1].Value : 0;
                if (i > period)
                {
                    sum -= !double.IsNaN(this[i - period - 1].Value) ? this[i - period - 1].Value : 0;
                    avg = sum / period;
                }
                else
                {
                    avg = sum / i;
                }
                if (i >= period)
                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(this[i - 1].Index, avg));
            }
            return timeSeries;
        }

        /// <summary>
        /// Returns a moving sum time-series based on the specified time period. The sum is computed based on the previous n=period ordinates.
        /// </summary>
        /// <param name="period">The time period to sum over. If time interval is 1-hour, and period is 12, the moving sum will be computed over a moving 12 hour block.</param>
        public TimeSeries MovingSum(int period)
        {
            if (period >= Count)
                throw new ArgumentException(nameof(period), "The period must be less than the length of the time-series.");
            SortByTime();
            var timeSeries = new TimeSeries(TimeInterval);
            double sum = 0d;
            for (int i = 1; i <= Count; i++)
            {
                sum += !double.IsNaN(this[i - 1].Value) ? this[i - 1].Value : 0;
                if (i > period)
                    sum -= !double.IsNaN(this[i - period - 1].Value) ? this[i - period - 1].Value : 0;
                if (i >= period)
                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(this[i - 1].Index, sum));
            }
            return timeSeries;
        }

        /// <summary>
        /// Shift all of the dates to match the new start date.
        /// </summary>
        /// <param name="newStartDate">The new start date.</param>
        public void ShiftAllDates(DateTime newStartDate)
        {
            if (Count == 0) return;
            bool wasSuppressed = SuppressCollectionChanged;
            SuppressCollectionChanged = true;
            this[0].Index = newStartDate;
            if (TimeInterval != TimeInterval.Irregular)
            {
                for (int i = 1; i < Count; i++) 
                { 
                    this[i].Index = AddTimeInterval(this[i - 1].Index, _timeInterval); 
                }
            }
            
            SuppressCollectionChanged = wasSuppressed;
            if (SuppressCollectionChanged == false) { RaiseCollectionChangedReset(); }
        }

        /// <summary>
        /// Shift the dates by a specified number of days.
        /// </summary>
        /// <param name="numberOfDays">The number of days to shift by.</param>
        /// <returns> A new TimeSeries object with the dates shifted</returns>
        public TimeSeries ShiftDatesByDay(int numberOfDays)
        {
            SortByTime();
            var timeSeries = new TimeSeries(TimeInterval);
            for (int i = 0; i < Count; i++)
            {
                var ordinate = new SeriesOrdinate<DateTime, double>();
                ordinate.Index = this[i].Index.AddDays(numberOfDays);
                ordinate.Value = this[i].Value;
                timeSeries.Add(ordinate);
            }
            return timeSeries;
        }

        /// <summary>
        /// Shift the dates by a specified number of months.
        /// </summary>
        /// <param name="numberOfMonths">The number of months to shift by.</param>
        /// <returns> A new TimeSeries object with the dates shifted</returns>
        public TimeSeries ShiftDatesByMonth(int numberOfMonths)
        {
            SortByTime();
            var timeSeries = new TimeSeries(TimeInterval);
            for (int i = 0; i < Count; i++)
            {
                var ordinate = new SeriesOrdinate<DateTime, double>();
                ordinate.Index = this[i].Index.AddMonths(numberOfMonths);
                ordinate.Value = this[i].Value;
                timeSeries.Add(ordinate);
            }
            return timeSeries;
        }

        /// <summary>
        /// Shift the dates by a specified number of years. 
        /// </summary>
        /// <param name="numberOfYears">The number of years to shift by.</param>
        /// <returns> A new TimeSeries object with the dates shifted</returns>
        public TimeSeries ShiftDatesByYear(int numberOfYears)
        {
            SortByTime();
            var timeSeries = new TimeSeries(TimeInterval);
            for (int i = 0; i < Count; i++)
            {
                var ordinate = new SeriesOrdinate<DateTime, double>();
                ordinate.Index = this[i].Index.AddYears(numberOfYears);
                ordinate.Value = this[i].Value;
                timeSeries.Add(ordinate);
            }
            return timeSeries;
        }

        /// <summary>
        /// Clip the time-series.
        /// </summary>
        /// <param name="startDate">The new start date/time of the series.</param>
        /// <param name="endDate">The new end date/time of the series.</param>
        /// <returns> A new TimeSeries object with the dates clipped</returns>
        public TimeSeries ClipTimeSeries(DateTime startDate, DateTime endDate)
        {
            if (startDate < StartDate)
                throw new ArgumentOutOfRangeException(nameof(startDate), "The start date is earlier than the start date of the time-series.");
            if (endDate > EndDate)
                throw new ArgumentOutOfRangeException(nameof(endDate), "The end date is later than the end date of the time-series.");
            var timeSeries = new TimeSeries(TimeInterval);
            for (int i = 0; i < Count; i++)
            {
                if (this[i].Index >= startDate && this[i].Index <= endDate)
                {
                    timeSeries.Add(this[i].Clone());
                }
            }
            return timeSeries;
        }

        /// <summary>
        /// Convert the current time-series to a new time interval.
        /// </summary>
        /// <param name="timeInterval">The new time interval.</param>
        /// <param name="average">Optional. Determines if values should be averaged (true) or cumulated (false) for larger time steps. Default = true.</param>
        /// <returns>A new TimeSeries object with the new interval.</returns>
        public TimeSeries? ConvertTimeInterval(TimeInterval timeInterval, bool average = true)
        {
            var TS = TimeSeries.TimeIntervalInHours(TimeInterval); // The time step in hours
            var newTS = TimeSeries.TimeIntervalInHours(timeInterval); // The new time step in hours
            int blockDuration = (int)Math.Floor(newTS / TS);
            double N = EndDate.Subtract(StartDate).TotalHours / newTS + 1;

            if (newTS == TS)
            {
                return Clone();
            }
            else if (newTS < TS && average == true)
            {
                // Create interpolater with existing data set.
                double t = 0, value = 0;
                var x = new double[Count];
                var y = new double[Count];
                for (int i = 0; i < Count; i++)
                {
                    x[i] = t;
                    y[i] = this[i].Value;
                    t += TS;
                }
                var linInt = new Linear(x, y);

                // Now Interpolate values for a smaller time step
                t = 0;
                var timeSeries = new TimeSeries(timeInterval);
                timeSeries.Add(new SeriesOrdinate<DateTime, double>(StartDate, this[0].Value));
                for (int i = 1; i < N; i++)
                {
                    t += newTS;
                    value = linInt.Interpolate(t);
                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(timeSeries[i - 1].Index, timeInterval), value));
                }
                return timeSeries;

            }
            else if (newTS > TS && average == true)
            {
                // Calculate block average
                int t = 0;
                var timeSeries = new TimeSeries(timeInterval);
                for (int i = 0; i < N; i++)
                {
                    double avg = 0;
                    for (int j = t; j < Math.Min(t + blockDuration, Count); j++)
                        avg += this[j].Value;
                    avg /= Math.Min(t + blockDuration, Count) - t;
                    t += blockDuration;
                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(i == 0 ? StartDate : AddTimeInterval(timeSeries[i - 1].Index, timeInterval), avg));
                }
                return timeSeries;
            }
            else if (newTS < TS && average == false) 
            {
                // Create interpolater with existing data set.
                double t = 0, value = 0, rate = TS / newTS;
                var x = new double[Count];
                var y = new double[Count];
                for (int i = 0; i < Count; i++)
                {
                    x[i] = t;
                    y[i] = this[i].Value;
                    t += TS;
                }
                var linInt = new Linear(x, y);

                // Now disaggregate values for a smaller time step
                t = 0;
                var timeSeries = new TimeSeries(timeInterval);
                timeSeries.Add(new SeriesOrdinate<DateTime, double>(StartDate, y[0] / rate));
                for (int i = 1; i < N; i++)
                {
                    t += newTS;
                    int idx = linInt.Search(t);
                    value = y[idx] / rate;
                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(timeSeries[i - 1].Index, timeInterval), value));
                }
                return timeSeries;
            }        
            else if (newTS > TS && average == false)
            {
                // Calculate block sum
                int t = 0;
                var timeSeries = new TimeSeries(timeInterval);
                for (int i = 0; i < N; i++)
                {
                    double sum = 0;
                    for (int j = t; j < Math.Min(t + blockDuration, Count); j++)
                        sum += this[j].Value;
                    t += blockDuration;
                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(i == 0 ? StartDate : AddTimeInterval(timeSeries[i - 1].Index, timeInterval), sum));
                }
                return timeSeries;
            }

            return null;
        }

        #region Summary Statistics

        /// <summary>
        /// Gets the min value of the time-series.
        /// </summary>
        public double MinValue()
        {
            double min = double.MaxValue;
            for (int i = 0; i < Count; i++)
            {
                if (!double.IsNaN(this[i].Value) && this[i].Value < min)
                {
                    min = this[i].Value;
                }
            }
            return min;
        }

        /// <summary>
        /// Gets the max value of the time-series.
        /// </summary>
        public double MaxValue()
        {
            double max = double.MinValue;
            for (int i = 0; i < Count; i++)
            {
                if (!double.IsNaN(this[i].Value) && this[i].Value > max)
                {
                    max = this[i].Value;
                }
            }
            return max;
        }

        /// <summary>
        /// Gets the mean of the time-series values.
        /// </summary>
        public double MeanValue()
        {
            if (Count == 0) return double.NaN;
            double mean = 0d;
            int n = 0;
            for (int i = 0; i < Count; i++)
            {
                if (!double.IsNaN(this[i].Value))
                {
                    mean += this[i].Value;
                    n += 1;
                }
            }
            return mean / n;
        }

        /// <summary>
        /// Gets the standard deviation of the time-series values.
        /// </summary>
        public double StandardDeviation()
        {
            if (Count < 2) return double.NaN;
            double variance = 0d;
            double t = this[0].Value;
            double n = 1;
            for (int i = 1; i < Count; i++)
            {
                if (!double.IsNaN(this[i].Value))
                {
                    t += this[i].Value;
                    double diff = (i + 1) * this[i].Value - t;
                    variance += diff * diff / ((i + 1.0d) * i);
                    n += 1;
                }
            }
            return Math.Sqrt(variance / (n - 1));
        }

        /// <summary>
        /// Returns summary percentile stats for the 5th, 25th, 50th, 75th, and 95th percentiles.
        /// </summary>
        public double[] SummaryPercentiles()
        {
            return Percentiles(new[] { 0.05d, 0.25d, 0.5d, 0.75d, 0.95d });
        }

        /// <summary>
        /// Returns an array of percentiles given a list of k-th percentile values.
        /// </summary>
        /// <param name="kValues">A list of k-th percentile values.</param>
        public double[] Percentiles(IList<double> kValues)
        {
            var data = ValuesToArray();
            Array.Sort(data);
            var result = Statistics.Statistics.Percentile(data, kValues, true);
            return result;
        }

        /// <summary>
        /// Returns the duration (percent of time exceedance curve)
        /// </summary>
        public double[,] Duration()
        {
            var result = new double[Count, 2];
            var pp = PlottingPositions.Weibull(Count);
            var data = ValuesToArray();
            Array.Sort(data);
            Array.Reverse(data);
            for (int i = 0; i < data.Length; i++)
            {
                result[i, 0] = pp[i] * 100;
                result[i, 1] = data[i];
            }
            return result;
        }

        /// <summary>
        /// Returns an array of percentiles given a list of k-th percentile values, for each month of the year.
        /// Number of rows = 12. Number of columns = length of the k-value list.
        /// </summary>
        /// <param name="kValues">A list of k-th percentile values.</param>
        public double[,] MonthlyPercentiles(IList<double> kValues)
        {
            var result = new double[12, kValues.Count];
            if (kValues == null || kValues.Count == 0) { return result; }

            Parallel.For(1, 13, index =>
            {
                // Filter data by month
                var monthlyData = new List<double>();
                for (int j = 0; j < Count; j++)
                {
                    if (this[j].Index.Month == index) { monthlyData.Add(this[j].Value); }
                }
                // Compute percentiles
                monthlyData.Sort();
                if (monthlyData.Count > 0)
                {
                    for (int j = 0; j < kValues.Count; j++)
                    {
                        result[index - 1, j] = Statistics.Statistics.Percentile(monthlyData, kValues[j], true);
                    }
                }
            });
            return result;
        }

        /// <summary>
        /// Returns an array of summary statistics for each month of the year. 
        /// Number of rows = 12. Number of columns = 8 {min, 5%, 25%, 50%, 75%, 95%, max, mean}.
        /// </summary>
        public double[,] MonthlySummaryStatistics()
        {
            var monthlySummary = new double[12, 8];
            Parallel.For(1, 13, index =>
            {
                // Filter data by month
                var monthlyData = new List<double>();
                for (int j = 0; j < Count; j++)
                {
                    if (this[j].Index.Month == index && !double.IsNaN(this[j].Value)) { monthlyData.Add(this[j].Value); }
                }
                if (monthlyData.Count == 0) { return; }
                // Compute percentiles
                monthlyData.Sort();
                monthlySummary[index - 1, 0] = monthlyData[0];
                monthlySummary[index - 1, 1] = Statistics.Statistics.Percentile(monthlyData, 0.05, true);
                monthlySummary[index - 1, 2] = Statistics.Statistics.Percentile(monthlyData, 0.25, true);
                monthlySummary[index - 1, 3] = Statistics.Statistics.Percentile(monthlyData, 0.5, true);
                monthlySummary[index - 1, 4] = Statistics.Statistics.Percentile(monthlyData, 0.75, true);
                monthlySummary[index - 1, 5] = Statistics.Statistics.Percentile(monthlyData, 0.95, true);
                monthlySummary[index - 1, 6] = monthlyData[monthlyData.Count - 1];
                monthlySummary[index - 1, 7] = Statistics.Statistics.ParallelMean(monthlyData);
            });
            return monthlySummary;
        }

        /// <summary>
        /// Returns a dictionary of the time series summary statistics.
        /// </summary>
        public Dictionary<string, double> SummaryStatistics()
        {
            var values = _seriesOrdinates.Where(y => !double.IsNaN(y.Value)).Select(x => x.Value).ToArray();
            var moments = Count <= 2 ? new double[] { double.NaN, double.NaN, double.NaN, double.NaN } : Statistics.Statistics.ProductMoments(values);
            var percentiles = Count <= 2 ? new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN, double.NaN } : Percentiles(new[] { 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99 });

            var result = new Dictionary<string, double>();
            result.Add("Record Length", Count);
            result.Add("Missing Values", NumberOfMissingValues());
            result.Add("Minimum", Statistics.Statistics.Minimum(values));
            result.Add("Maximum", Statistics.Statistics.Maximum(values));
            result.Add("Mean", moments[0]);
            result.Add("Std Dev", moments[1]);
            result.Add("Skewness", moments[2]);
            result.Add("Kurtosis", moments[3]);
            result.Add("1%", percentiles[0]);
            result.Add("5%", percentiles[1]);
            result.Add("25%", percentiles[2]);
            result.Add("50%", percentiles[3]);
            result.Add("75%", percentiles[4]);
            result.Add("95%", percentiles[5]);
            result.Add("99%", percentiles[6]);

            return result;
        }

        /// <summary>
        /// Returns a dictionary of the hypothesis test results.
        /// </summary>
        /// <param name="splitLocation">The location in the series to split the data samples.</param>
        public Dictionary<string, double> SummaryHypothesisTest(int splitLocation = -1)
        {
            var values = _seriesOrdinates.Where(y => !double.IsNaN(y.Value)).Select(x => x.Value).ToArray();
            splitLocation = splitLocation < 0 ? (int)((double)values.Length / 2) : splitLocation;
            var v1 = values.Subset(0, splitLocation);
            var v2 = values.Subset(splitLocation + 1, values.Length - 1);

            var result = new Dictionary<string, double>();
            result.Add("Jarque-Bera test for normality", HypothesisTests.JarqueBeraTest(values));
            result.Add("Ljung-Box test for independence", HypothesisTests.LjungBoxTest(values));
            result.Add("Wald-Wolfowitz test for independence and stationarity (trend)", HypothesisTests.WaldWolfowitzTest(values));
            result.Add("Mann-Whitney test for homogeneity and stationarity (jump)", HypothesisTests.MannWhitneyTest(v1.Length <= v2.Length ? v1 : v2, v1.Length > v2.Length ? v1 : v2));
            result.Add("Mann-Kendall test for homogeneity and stationarity (trend)", HypothesisTests.MannKendallTest(values));
            result.Add("t-test for differences in the means of two samples", HypothesisTests.UnequalVarianceTtest(v1, v2));
            result.Add("F-test for differences in the variances of two samples", HypothesisTests.Ftest(v1, v2));

            return result;
        }

        #endregion

        #region Frequency Analysis Methods

        /// <summary>
        /// Compute the monthly frequency of occurrence.
        /// </summary>
        public double[] MonthlyFrequency()
        {
            var frequencies = new double[12];
            for (int i = 1; i <= 12; i++)
                frequencies[i - 1] = (double)_seriesOrdinates.Where(x => x.Index.Month == i).ToList().Count;
            return frequencies;
        }

        /// <summary>
        /// Returns an annual (irregular) block series.  
        /// </summary>
        /// <param name="blockFunction">Optional. The block function type; e.g. min, max, sum, or average. Default = Maximum.</param>
        /// <param name="smoothingFunction">Optional. The smoothing function type. Default = None.</param>
        /// <param name="period">Optional. The time period to perform smoothing over. If time interval is 1-hour, and period = 12, the smoothing will be computed over a moving 12 hour block. Default = 1.</param>
        public TimeSeries CalendarYearSeries(BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            var result = new TimeSeries(TimeInterval.Irregular);

            // First, perform smoothing function
            TimeSeries? smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? Clone() : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? Clone() : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            // Then, perform block function
            if (smoothedSeries == null) return result;
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {
                var blockData = smoothedSeries.Where(x => x.Index.Year == i).ToList();
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
                    {
                        sum += !double.IsNaN(blockData[j].Value) ? blockData[j].Value : 0;
                    }
                    ordinate.Index = blockData.Last().Index;
                    ordinate.Value = sum;
                }
                else if (blockFunction == BlockFunctionType.Average)
                {
                    double sum = 0;
                    for (int j = 0; j < blockData.Count; j++)
                    {
                        sum += !double.IsNaN(blockData[j].Value) ? blockData[j].Value : 0;
                    }
                    ordinate.Index = blockData.Last().Index;
                    ordinate.Value = sum / blockData.Count;
                }

                if (!double.IsNaN(ordinate.Value))
                    result.Add(ordinate);
            }
            return result;
        }

        /// <summary>
        /// Returns an annual (irregular) block series based on a 12-month year with a customized starting month. 
        /// </summary>
        /// <param name="startMonth">Optional. The month when the custom year begins. Default = 10 (or October).</param>
        /// <param name="blockFunction">Optional. The block function type; e.g. min, max, sum, or average. Default = Maximum.</param>
        /// <param name="smoothingFunction">Optional. The smoothing function type. Default = None.</param>
        /// <param name="period">Optional. The time period to perform smoothing over. If time interval is 1-hour, and period = 12, the smoothing will be computed over a moving 12 hour block. Default = 1.</param>
        public TimeSeries CustomYearSeries(int startMonth = 10, BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            if (startMonth < 1 || startMonth > 12)
            {
                throw new ArgumentOutOfRangeException(nameof(startMonth), "The start month be between 1 and 12.");
            }

            var result = new TimeSeries(TimeInterval.Irregular);

            // First, perform smoothing function
            TimeSeries? smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? Clone() : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? Clone() : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            // Then, shift the dates
            if( smoothedSeries == null) return result;
            int shift = startMonth != 1 ? 12 - startMonth + 1 : 0;
            smoothedSeries = startMonth != 1 ? smoothedSeries.ShiftDatesByMonth(shift) : smoothedSeries;

            // Then, perform block function
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {
                var blockData = smoothedSeries.Where(x => x.Index.Year == i).ToList();
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
                    {
                        sum += !double.IsNaN(blockData[j].Value) ? blockData[j].Value : 0;
                    }
                    ordinate.Index = blockData.Last().Index;
                    ordinate.Value = sum;
                }
                else if (blockFunction == BlockFunctionType.Average)
                {
                    double sum = 0;
                    for (int j = 0; j < blockData.Count; j++)
                    {
                        sum += !double.IsNaN(blockData[j].Value) ? blockData[j].Value : 0;
                    }
                    ordinate.Index = blockData.Last().Index;
                    ordinate.Value = sum / blockData.Count;
                }

                if (!double.IsNaN(ordinate.Value))
                    result.Add(ordinate);
            }

            // Finally, shift the dates back
            result = startMonth != 1 ? result.ShiftDatesByMonth(-shift) : result;
            return result;
        }

        /// <summary>
        /// Returns a custom annual (irregular) block series based on a customized window. 
        /// </summary>
        /// <param name="startMonth">The month when the custom year begins.</param>
        /// <param name="endMonth">The month when the custom year ends.</param>
        /// <param name="blockFunction">Optional. The block function type; e.g. min, max, sum, or average. Default = Maximum.</param>
        /// <param name="smoothingFunction">Optional. The smoothing function type. Default = None.</param>
        /// <param name="period">Optional. The time period to perform smoothing over. If time interval is 1-hour, and period = 12, the smoothing will be computed over a moving 12 hour block. Default = 1.</param>
        public TimeSeries CustomYearSeries(int startMonth, int endMonth, BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            if (startMonth < 1 || startMonth > 12)
            {
                throw new ArgumentOutOfRangeException(nameof(startMonth), "The start month be between 1 and 12.");
            }
            if (endMonth < 1 || endMonth > 12)
            {
                throw new ArgumentOutOfRangeException(nameof(endMonth), "The start month be between 1 and 12.");
            }

            var result = new TimeSeries(TimeInterval.Irregular);

            // First, perform smoothing function
            TimeSeries? smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? Clone() : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? Clone() : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            // Then, perform block function
            if(smoothedSeries == null) return result;
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {

                var blockData = new List<SeriesOrdinate<DateTime, double>>();
                var ordinate = new SeriesOrdinate<DateTime, double>() { Value = double.NaN };

                if (startMonth <= endMonth)
                {
                    // Months are within a single calendar year
                    for (int j = startMonth; j <= endMonth; j++)
                    {
                        blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == j).ToList());
                    }
                }
                else
                {
                    // Months overlap two calendar years
                    // First get the previous year's months
                    for (int j = startMonth; j <= 12; j++)
                    {
                        blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i - 1 && x.Index.Month == j).ToList());
                    }
                    // Then the current year's months
                    for (int j = 1; j <= endMonth; j++)
                    {
                        blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == j).ToList());
                    }
                }

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
                    {
                        sum += !double.IsNaN(blockData[j].Value) ? blockData[j].Value : 0;
                    }
                    ordinate.Index = blockData.Last().Index;
                    ordinate.Value = sum;
                }
                else if (blockFunction == BlockFunctionType.Average)
                {
                    double sum = 0;
                    for (int j = 0; j < blockData.Count; j++)
                    {
                        sum += !double.IsNaN(blockData[j].Value) ? blockData[j].Value : 0;
                    }
                    ordinate.Index = blockData.Last().Index;
                    ordinate.Value = sum / blockData.Count;
                }

                if (!double.IsNaN(ordinate.Value))
                    result.Add(ordinate);
            }
            return result;
        }

        /// <summary>
        /// Returns an monthly (irregular) block series.  
        /// </summary>
        /// <param name="blockFunction">Optional. The block function type; e.g. min, max, sum, or average. Default = Maximum.</param>
        /// <param name="smoothingFunction">Optional. The smoothing function type. Default = None.</param>
        /// <param name="period">Optional. The time period to perform smoothing over. If time interval is 1-hour, and period = 12, the smoothing will be computed over a moving 12 hour block. Default = 1.</param>
        public TimeSeries MonthlySeries(BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            var result = new TimeSeries(TimeInterval.Irregular);

            // Create smoothed series
            TimeSeries? smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? Clone() : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? Clone() : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            if(smoothedSeries == null) return result;
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {

                for (int k = 1; k <= 12; k++)
                {
                    var blockData = smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == k).ToList();
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
                        {
                            sum += blockData[j].Value;
                        }
                        ordinate.Index = blockData.Last().Index;
                        ordinate.Value = sum;
                    }
                    else if (blockFunction == BlockFunctionType.Average)
                    {
                        double sum = 0;
                        for (int j = 0; j < blockData.Count; j++)
                        {
                            sum += blockData[j].Value;
                        }
                        ordinate.Index = blockData.Last().Index;
                        ordinate.Value = sum / blockData.Count;
                    }

                    if (!double.IsNaN(ordinate.Value))
                        result.Add(ordinate);
                }


            }
            return result;
        }

        /// <summary>
        /// Returns an quarterly (irregular) block series.  
        /// </summary>
        /// <param name="blockFunction">Optional. The block function type; e.g. min, max, sum, or average. Default = Maximum.</param>
        /// <param name="smoothingFunction">Optional. The smoothing function type. Default = None.</param>
        /// <param name="period">Optional. The time period to perform smoothing over. If time interval is 1-hour, and period = 12, the smoothing will be computed over a moving 12 hour block. Default = 1.</param>
        public TimeSeries QuarterlySeries(BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            var qStart = new int[] { 1, 4, 7, 10 };
            var qEnd = new int[] { 3, 6, 9, 12 };

            var result = new TimeSeries(TimeInterval.Irregular);

            // Create smoothed series
            TimeSeries? smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? Clone() : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? Clone() : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            if (smoothedSeries == null) return result;
            for (int i = smoothedSeries.StartDate.Year; i <= smoothedSeries.EndDate.Year; i++)
            {

                for (int q = 0; q < qEnd.Length; q++)
                {
                    var blockData = new List<SeriesOrdinate<DateTime, double>>();
                    for (int j = qStart[q]; j <= qEnd[q]; j++)
                    {
                        blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == j).ToList());
                    }

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
                        {
                            sum += blockData[j].Value;
                        }
                        ordinate.Index = blockData.Last().Index;
                        ordinate.Value = sum;
                    }
                    else if (blockFunction == BlockFunctionType.Average)
                    {
                        double sum = 0;
                        for (int j = 0; j < blockData.Count; j++)
                        {
                            sum += blockData[j].Value;
                        }
                        ordinate.Index = blockData.Last().Index;
                        ordinate.Value = sum / blockData.Count;
                    }

                    if (!double.IsNaN(ordinate.Value))
                        result.Add(ordinate);

                }

            }
            return result;
        }

        /// <summary>
        /// Returns a peaks-over-threshold (POT) series.
        /// </summary>
        /// <param name="threshold">The threshold value.</param>
        /// <param name="minStepsBetweenEvents">The minimum number of time steps between independent peak events. This time condition ensures independence between events. Default = 1.</param>
        /// <param name="smoothingFunction">The smoothing function type. Smoothing is performed before the peaks-over-threshold analysis.</param>
        /// <param name="period">The time period to perform smoothing over. If time interval is 1-hour, and period is 12. The smoothing will be computed over a moving 12 hour block.</param>
        /// <remarks>
        /// This routine is based on the 'clust' method included in the POT R package (https://cran.r-project.org/web/packages/POT/index.html).
        /// The clusters of exceedances are defined as follows:
        /// <list type="bullet">
        /// <item>
        /// The first exceedance initiates the first cluster;
        /// </item>
        /// <item>
        /// The first observation under the threshold u “ends” the current cluster unless the minimum steps between events does not hold;
        /// </item>
        /// <item>
        /// The next exceedance initiates a new cluster;
        /// </item>
        /// </list>
        /// </remarks>
        public TimeSeries PeaksOverThresholdSeries(double threshold, int minStepsBetweenEvents = 1, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            // Create smoothed time series
            TimeSeries? smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? Clone() : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? Clone() : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            // First, create the cluster indexes. 
            int i = 0, idx, idxMax;
            var clusters = new List<int[]>();

            if(smoothedSeries == null) return new TimeSeries(TimeInterval.Irregular);
            while (i < smoothedSeries.Count)
            {
                if (!double.IsNaN(smoothedSeries[i].Value) && smoothedSeries[i].Value > threshold)
                {
                    // Set the start of the cluster
                    clusters.Add(new int[2]);
                    clusters.Last()[0] = i;
                    idx = i + 1;
                    idxMax = idx;

                    while ((!double.IsNaN(smoothedSeries[idx].Value) && smoothedSeries[idx].Value > threshold) ||
                        CheckIfMinStepsExceeded(smoothedSeries[idxMax].Index, smoothedSeries[idx].Index, minStepsBetweenEvents) == false)
                    {
                        if (!double.IsNaN(smoothedSeries[idx].Value) && smoothedSeries[idx].Value >= smoothedSeries[idxMax].Value)
                        {
                            idxMax = idx;
                        }

                        // Increment inner loop
                        idx++;
                        if (idx >= smoothedSeries.Count)
                        {
                            idx--;
                            break;
                        }
                    }

                    // Set the end of the cluster
                    clusters.Last()[1] = idx - 1;
                    // Increment outer loop
                    i = idx + 1;
                }
                else
                {
                    i++;
                }
            }

            // Next get the max values within each cluster
            var result = new TimeSeries(TimeInterval.Irregular);
            for (i = 0; i < clusters.Count; i++)
            {
                var max = smoothedSeries[clusters[i][0]].Clone();
                for (int j = clusters[i][0] + 1; j <= clusters[i][1]; j++)
                {
                    if (smoothedSeries[j].Value >= max.Value)
                    {
                        max = smoothedSeries[j].Clone();
                    }
                }
                result.Add(max);
            }
            return result;
        }

        #endregion

        /// <summary>
        /// Resample the time series using k-Nearest neighbors.
        /// </summary>
        /// <param name="timeSteps">The number of time steps to resample.</param>
        /// <param name="k">The number of nearest neighbors.</param>
        /// <param name="seed">The prng seed. Default = 12345.</param>
        /// <returns></returns>
        public TimeSeries ResampleWithKNN(int timeSteps, int k, int seed = 12345)
        {
            if (Count < 2) throw new ArgumentException("Need at least 2 points for resampling.");
            if (timeSteps < 1) throw new ArgumentException("The number of time steps must be at least 1.");

            // Initialize
            var timeSeries = new TimeSeries(TimeInterval);
            var prng = new MersenneTwister(seed);
            int startIdx = prng.Next(Count - 1);
            double currentValue = this[startIdx].Value;
            DateTime currentDate = StartDate;
            timeSeries.Add(new SeriesOrdinate<DateTime, double>(currentDate, currentValue));

            // Standardize the values
            double mean = MeanValue();
            double stdDev = StandardDeviation();
            var stdData = this.Clone();
            stdData.Standardize();
            var stdValues = stdData.ValuesToArray();

            // Perform k-NN
            var kNN = new KNearestNeighbors(stdValues, stdValues, k);
            for (int i = 1; i < timeSteps; i++)
            {
                double val = (currentValue - mean) / stdDev;
                var kNearest = kNN.GetNeighbors([val]);
                int selectedIdx = kNearest[prng.Next(k)];
                currentValue = this[selectedIdx].Value;
                currentDate = TimeSeries.AddTimeInterval(currentDate, TimeInterval);
                timeSeries.Add(new SeriesOrdinate<DateTime, double>(currentDate, currentValue));
            }

            return timeSeries;
        }

        /// <summary>
        /// Resample the time series using the block bootstrap.
        /// </summary>
        /// <param name="timeSteps">The number of time steps to resample.</param>
        /// <param name="blockSize">The resampling block size.</param>
        /// <param name="seed">The prng seed. Default = 12345.</param>
        /// <returns></returns>
        public TimeSeries ResampleWithBlockBootstrap(int timeSteps, int blockSize, int seed = 12345)
        {
            if (Count < 2) throw new ArgumentException("Need at least 2 points for resampling.");
            if (Count < blockSize) throw new ArgumentException("Block size is too large for this time series data.");
            if (timeSteps < 1) throw new ArgumentException("The number of time steps must be at least 1.");

            // Step 1: collect all blocks of length B
            var data = ValuesToArray();
            List<double[]> blocks = new List<double[]>();
            for (int i = 0; i <= Count - blockSize; i++)
            {
                double[] block = new double[blockSize];
                Array.Copy(data, i, block, 0, blockSize);
                blocks.Add(block);
            }

            // Step 2: sample blocks and build synthetic series
            var prng = new MersenneTwister(seed);
            List<double> syntheticSeries = new List<double>();
            while (syntheticSeries.Count < timeSteps)
            {
                int index = prng.Next(blocks.Count);
                syntheticSeries.AddRange(blocks[index]);
            }

            // Step 3: trim to desired length
            var timeSeries = new TimeSeries(TimeInterval, StartDate, syntheticSeries.Take(timeSteps).ToArray());
            return timeSeries;
        }

        /// <summary>
        /// Returns an XElement of a series ordinate.
        /// </summary>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(TimeSeries));
            result.SetAttributeValue(nameof(TimeInterval), TimeInterval.ToString());
            for (int i = 0; i < Count; i++)
            {
                var ordinate = new XElement("SeriesOrdinate");
                ordinate.SetAttributeValue("Index", this[i].Index.ToString("o", CultureInfo.InvariantCulture));
                ordinate.SetAttributeValue("Value", this[i].Value.ToString("G17", CultureInfo.InvariantCulture));
                result.Add(ordinate);
            }
            return result;
        }

        /// <summary>
        /// Creates a copy of the time series.
        /// </summary>
        public TimeSeries Clone()
        {
            return new TimeSeries(ToXElement());
        }

    }
}