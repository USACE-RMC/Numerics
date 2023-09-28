﻿using Numerics.Data.Statistics;
using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Threading.Tasks;
using System.Xml.Linq;

namespace Numerics.Data
{

    /// <summary>
    /// A time-series class, which is a collection of time-series ordinates.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
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

            _timeInterval = timeInterval;
            Add(new SeriesOrdinate<DateTime, double>(startDate, double.NaN));
            while (_seriesOrdinates.Last().Index < endDate)
                Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(_seriesOrdinates.Last().Index, _timeInterval), double.NaN));                   
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

            _timeInterval = timeInterval;
            Add(new SeriesOrdinate<DateTime, double>(startDate, fixedValue));
            while (_seriesOrdinates.Last().Index < endDate)
                Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(_seriesOrdinates.Last().Index, _timeInterval), fixedValue));
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
                Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(this[i - 1].Index, _timeInterval), data[i]));
        }

        /// <summary>
        /// Constructs a time-series based on XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        public TimeSeries(XElement xElement)
        {
            // Get time interval
            if (xElement.Attribute(nameof(TimeInterval)) != null)         
                Enum.TryParse(xElement.Attribute(nameof(TimeInterval)).Value, out _timeInterval);

            // Get Ordinates
            foreach (XElement ordinate in xElement.Elements("SeriesOrdinate"))
            {
                DateTime index = DateTime.Now;
                DateTime.TryParse(ordinate.Attribute("Index").Value, out index);
                double value = 0;
                double.TryParse(ordinate.Attribute("Value").Value, out value);
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
                if (!double.IsNaN(this[i].Value))
                    this[i].Value += constant;
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
                if (!double.IsNaN(this[i].Value))
                    this[i].Value -= constant;
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
                if (!double.IsNaN(this[i].Value))
                    this[i].Value *= constant;
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
                if (!double.IsNaN(this[i].Value))
                    this[i].Value /= constant;
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
                if (!double.IsNaN(this[i].Value))
                    this[i].Value = Math.Abs(this[i].Value);
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
        /// Log transform value in the time-series. Missing values are kept as missing.
        /// </summary>
        /// <param name="baseValue">The log base value.</param>
        public void LogTransform(double baseValue = 10)
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i <= Count - 1; i++)
            {
                if (this[i].Value > 0 && !double.IsNaN(this[i].Value))
                    this[i].Value = Math.Log(this[i].Value, baseValue);
                else if (this[i].Value <=0 || double.IsNaN(this[i].Value))
                    this[i].Value = double.NaN;
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        public void Standardize()
        {
            SuppressCollectionChanged = true;
            double mean = MeanValue();
            double stdDev = StandardDeviation();
            for (int i = 0; i <= Count - 1; i++)
            {
                this[i].Value = (this[i].Value - mean) / stdDev;
            }
            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Raise each value in the time-series is replaced by its inverse (1/x). Missing values are kept as missing. If the value is 0.0, the value is set to Double.NaN.
        /// </summary>
        public void Inverse()
        {
            SuppressCollectionChanged = true;
            for (int i = 0; i <= Count - 1; i++)
            {
                if (this[i].Value != 0 && !double.IsNaN(this[i].Value))
                    this[i].Value = 1d / this[i].Value;
                else if (this[i].Value == 0 || double.IsNaN(this[i].Value))
                    this[i].Value = double.NaN;
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
        /// Returns a time-series of the successive differences per time period.
        /// </summary>
        /// <param name="period">Time period for taking differences. If time interval is 1-hour, and period is 12, the difference will be computed over a moving 12 hour block.</param>
        public TimeSeries Difference(int period = 1)
        {
            var timeSeries = new TimeSeries();
            for (int i = period; i < Count; i++)
            {
                timeSeries.Add(this[i].Clone());
                timeSeries.Last().Value = this[i].Value - this[i - period].Value;
            }
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
            for (int i = 0; i <= Count - 1; i++)
                if (double.IsNaN(this[i].Value))
                    this[i].Value = value;
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
            double y;
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
        /// Determines if the minimum step between events has been exceeded. 
        /// </summary>
        /// <param name="startTime">Start time of the starting event.</param>
        /// <param name="endTime">End time of the ending event.</param>
        /// <param name="minStepsBetweenEvents">Minimum time steps between events.</param>
        private bool CheckIfMinStepsExceeded(DateTime startTime, DateTime endTime, int minStepsBetweenEvents)
        {
            DateTime _endTime;
            switch (TimeInterval)
            {
                case TimeInterval.OneMinute:
                    {
                        _endTime = startTime.AddMinutes(1 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.FiveMinute:
                    {
                        _endTime = startTime.AddMinutes(5 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.FifteenMinute:
                    {
                        _endTime = startTime.AddMinutes(15 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.ThirtyMinute:
                    {
                        _endTime = startTime.AddMinutes(30 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.OneHour:
                    {
                        _endTime = startTime.AddHours(1 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.SixHour:
                    {
                        _endTime = startTime.AddHours(6 * minStepsBetweenEvents);
                        return endTime <= _endTime;
                    }

                case TimeInterval.TwelveHour:
                    {
                        _endTime = startTime.AddHours(12 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.OneDay:
                    {
                        _endTime = startTime.AddDays(minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.SevenDay:
                    {
                        _endTime = startTime.AddDays(7 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.OneMonth:
                    {
                        _endTime = startTime.AddMonths(1 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.OneQuarter:
                    {
                        _endTime = startTime.AddMonths(3 * minStepsBetweenEvents);
                        return endTime > _endTime;
                    }

                case TimeInterval.OneYear:
                    {
                        _endTime = startTime.AddYears(1 * minStepsBetweenEvents);
                        return endTime > _endTime;
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
        /// <param name="period">The time period to sum over. If time interval is 1-hour, and period is 12. The moving sum will be computed over a moving 12 hour block.</param>
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
            SuppressCollectionChanged = true;
            this[0].Index = newStartDate;
            for (int i = 1; i < Count; i++)
                this[i].Index = AddTimeInterval(this[i - 1].Index, _timeInterval);

            SuppressCollectionChanged = false;
            RaiseCollectionChangedReset();
        }

        /// <summary>
        /// Shift the dates by a specified number of months.
        /// </summary>
        /// <param name="numberOfMonths">The number of months to shift by.</param>
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
        public TimeSeries ConvertTimeInterval(TimeInterval timeInterval)
        {
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
            if (Count == 0) return double.NaN;
            double variance_ = 0d;
            double t = this[0].Value;
            int n = 0;
            for (int i = 1; i < Count; i++)
            {
                if (!double.IsNaN(this[i].Value))
                {
                    t += this[i].Value;
                    double diff = (i + 1) * this[i].Value - t;
                    variance_ += diff * diff / ((i + 1.0d) * i);
                    n += 1;
                }
            }
            return Math.Sqrt(variance_ / (n - 1));
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
            var perc = new double[kValues.Count];
            var data = ValuesToArray();
            Array.Sort(data);
            for (int i = 0; i < kValues.Count; i++)
                perc[i] = Statistics.Statistics.Percentile(data, kValues[i], true);
            return perc;
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
            var monthlyPercValues = new double[12, kValues.Count];
            Parallel.For(1, 13, index =>
            {
                // Filter data by month
                var monthlyData = new List<double>();
                for (int j = 0; j < Count; j++)
                {
                    if (this[j].Index.Month == index)
                        monthlyData.Add(this[j].Value);
                }
                // Compute percentiles
                monthlyData.Sort();
                for (int j = 0; j < kValues.Count; j++)
                    monthlyPercValues[index - 1, j] = Statistics.Statistics.Percentile(monthlyData, kValues[j], true);
            });
            return monthlyPercValues;
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
                    if (this[j].Index.Month == index)
                        monthlyData.Add(this[j].Value);
                }
                // Compute percentiles
                monthlyData.Sort();
                monthlySummary[index - 1, 0] = monthlyData.First();
                monthlySummary[index - 1, 1] = Statistics.Statistics.Percentile(monthlyData, 0.05, true);
                monthlySummary[index - 1, 2] = Statistics.Statistics.Percentile(monthlyData, 0.25, true);
                monthlySummary[index - 1, 3] = Statistics.Statistics.Percentile(monthlyData, 0.5, true);
                monthlySummary[index - 1, 4] = Statistics.Statistics.Percentile(monthlyData, 0.75, true);
                monthlySummary[index - 1, 5] = Statistics.Statistics.Percentile(monthlyData, 0.95, true);
                monthlySummary[index - 1, 6] = monthlyData.Last();
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
            var percentiles = Count <= 2 ? new double[] { double.NaN, double.NaN, double.NaN, double.NaN, double.NaN } : Percentiles(new[] { 0.05, 0.25, 0.5, 0.75, 0.95 });

            var result = new Dictionary<string, double>();
            result.Add("Record Length", Count);
            result.Add("Missing Values", NumberOfMissingValues());
            result.Add("Minimum", Statistics.Statistics.Minimum(values));
            result.Add("Maximum", Statistics.Statistics.Maximum(values));
            result.Add("Mean", moments[0]);
            result.Add("Std Dev", moments[1]);
            result.Add("Skewness", moments[2]);
            result.Add("Kurtosis", moments[3]);
            result.Add("5%", percentiles[0]);
            result.Add("25%", percentiles[1]);
            result.Add("50%", percentiles[2]);
            result.Add("75%", percentiles[3]);
            result.Add("95%", percentiles[4]);

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
        public int[] MonthlyFrequency()
        {
            var frequencies = new int[12];
            for (int i = 1; i <= 12; i++)
                frequencies[i - 1] = _seriesOrdinates.Where(x => x.Index.Month == i).ToList().Count;
            return frequencies;
        }

        /// <summary>
        /// Creates an annual max series.
        /// </summary>
        /// <param name="startMonth">The month when the year begins. If not 1, dates are shifted.</param>
        /// <param name="period">The time period to average or sum over.</param>
        /// <param name="isMovingAverage">If true, a moving average is performed, if false, a moving sum is performed.</param>
        public TimeSeries AnnualMaxSeries(int startMonth = 1, int period = 1, bool isMovingAverage = true)
        {
            var result = new TimeSeries(TimeInterval.Irregular);
            // Shift series
            int shift = startMonth != 1 ? 12 - startMonth + 1 : startMonth;
            var shiftedSeries = startMonth != 1 ? ShiftDatesByMonth(shift) : this;
            // Get moving series
            var movingSeries = period == 1 ? shiftedSeries : isMovingAverage ? shiftedSeries.MovingAverage(period) : shiftedSeries.MovingSum(period);
            // Create block max series
            for (int i = movingSeries.StartDate.Year; i <= movingSeries.EndDate.Year; i++)
            {
                int y = i;
                var blockData = movingSeries.Where(x => x.Index.Year == y).ToList();
                double max = double.MinValue;
                var maxOrdinate = new SeriesOrdinate<DateTime, double>();
                for (int j = 0; j < blockData.Count; j++)
                {
                    if (blockData[j].Value > max)
                    {
                        max = blockData[j].Value;
                        maxOrdinate.Index = blockData[j].Index;
                        maxOrdinate.Value = blockData[j].Value;
                    }
                }
                result.Add(maxOrdinate);
            }
            return result;
        }

        /// <summary>
        /// Creates a seasonal annual max series. 
        /// </summary>
        /// <param name="months">The months that define the season.</param>
        /// <param name="period">The time period to average or sum over.</param>
        /// <param name="isMovingAverage">If true, a moving average is performed, if false, a moving sum is performed.</param>
        public TimeSeries AnnualMaxSeries(int[] months, int period = 1, bool isMovingAverage = true)
        {
            var result = new TimeSeries(TimeInterval.Irregular);
            // Get moving series
            var movingSeries = isMovingAverage ? MovingAverage(period) : MovingSum(period);
            // Create block max series
            for (int i = movingSeries.StartDate.Year; i <= movingSeries.EndDate.Year; i++)
            {
                int y = i;
                var blockData = new List<SeriesOrdinate<DateTime, double>>();
                for (int j = 0; j < months.Length; j++)
                {
                    blockData.AddRange(movingSeries.Where(x => x.Index.Year == y && x.Index.Month == months[j]).ToList());
                }
                double max = double.MinValue;
                var maxOrdinate = new SeriesOrdinate<DateTime, double>();
                for (int j = 0; j < blockData.Count; j++)
                {
                    if (blockData[j].Value > max)
                    {
                        max = blockData[j].Value;
                        maxOrdinate.Index = blockData[j].Index;
                        maxOrdinate.Value = blockData[j].Value;
                    }
                }
                result.Add(maxOrdinate);
            }
            return result;
        }

        /// <summary>
        /// Returns an annual (irregular) block series.  
        /// </summary>
        /// <param name="blockFunction">The block function type; e.g. min, max, sum, or average.</param>
        /// <param name="smoothingFunction">The smoothing function type.</param>
        /// <param name="period">The time period to perform smoothing over. If time interval is 1-hour, and period is 12. The smoothing will be computed over a moving 12 hour block.</param>
        public TimeSeries CalendarYearSeries(BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            var result = new TimeSeries(TimeInterval.Irregular);

            // Create smoothed series
            TimeSeries smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = this.Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? this : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? this : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            for (int i = StartDate.Year; i <= EndDate.Year; i++)
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
        /// Returns an annual (irregular) block series based on the water year. 
        /// </summary>
        /// <param name="startMonth">The month when the water year begins. If not 1, dates are shifted.</param>
        /// <param name="blockFunction">The block function type; e.g. min, max, sum, or average.</param>
        /// <param name="smoothingFunction">The smoothing function type.</param>
        /// <param name="period">The time period to perform smoothing over. If time interval is 1-hour, and period is 12. The smoothing will be computed over a moving 12 hour block.</param>
        public TimeSeries WaterYearSeries(int startMonth = 10, BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            // Shift series
            int shift = startMonth != 1 ? 12 - startMonth + 1 : startMonth;
            var shiftedSeries = startMonth != 1 ? ShiftDatesByMonth(shift) : this;
            return shiftedSeries.CalendarYearSeries(blockFunction, smoothingFunction, period);
        }

        /// <summary>
        /// Returns a custom annual (irregular) block series. 
        /// </summary>
        /// <param name="startMonth">The month when the season begins.</param>
        /// <param name="endMonth">The month when the season ends.</param>
        /// <param name="blockFunction">The block function type; e.g. min, max, sum, or average.</param>
        /// <param name="smoothingFunction">The smoothing function type.</param>
        /// <param name="period">The time period to perform smoothing over. If time interval is 1-hour, and period is 12. The smoothing will be computed over a moving 12 hour block.</param>
        public TimeSeries CustomYearSeries(int startMonth = 1, int endMonth = 12, BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            var months = new List<int>();
            if (startMonth <= endMonth)
            {
                for (int i = startMonth; i <= endMonth; i++)
                    months.Add(i);
            }
            else
            {
                for (int i = startMonth; i <= 12; i++)
                    months.Add(i);
                for (int i = 1; i <= endMonth; i++)
                    months.Add(i);
            }

            var result = new TimeSeries(TimeInterval.Irregular);

            // Create smoothed series
            TimeSeries smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = this.Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? this : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? this : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            for (int i = StartDate.Year; i <= EndDate.Year; i++)
            {

                var blockData = new List<SeriesOrdinate<DateTime, double>>();
                for (int j = 0; j < months.Count; j++)
                {
                    blockData.AddRange(smoothedSeries.Where(x => x.Index.Year == i && x.Index.Month == months[j]).ToList());
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
        /// <param name="blockFunction">The block function type; e.g. min, max, sum, or average.</param>
        /// <param name="smoothingFunction">The smoothing function type.</param>
        /// <param name="period">The time period to perform smoothing over. If time interval is 1-hour, and period is 12. The smoothing will be computed over a moving 12 hour block.</param>
        public TimeSeries MonthlySeries(BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            var result = new TimeSeries(TimeInterval.Irregular);

            // Create smoothed series
            TimeSeries smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = this.Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? this : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? this : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            for (int i = StartDate.Year; i <= EndDate.Year; i++)
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
        /// <param name="blockFunction">The block function type; e.g. min, max, sum, or average.</param>
        /// <param name="smoothingFunction">The smoothing function type.</param>
        /// <param name="period">The time period to perform smoothing over. If time interval is 1-hour, and period is 12. The smoothing will be computed over a moving 12 hour block.</param>
        public TimeSeries QuarterlySeries(BlockFunctionType blockFunction = BlockFunctionType.Maximum, SmoothingFunctionType smoothingFunction = SmoothingFunctionType.None, int period = 1)
        {
            var qStart = new int[] { 1, 4, 7, 10 };
            var qEnd = new int[] { 3, 6, 9, 12 };

            var result = new TimeSeries(TimeInterval.Irregular);

            // Create smoothed series
            TimeSeries smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = this.Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? this : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? this : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }

            for (int i = StartDate.Year; i <= EndDate.Year; i++)
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
        /// This routine is based on the "clust" method included in the POT R package (https://cran.r-project.org/web/packages/POT/index.html).
        /// The clusters of exceedances are defines as follows:
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
            TimeSeries smoothedSeries = null;
            if (smoothingFunction == SmoothingFunctionType.None)
            {
                smoothedSeries = this.Clone();
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingAverage)
            {
                smoothedSeries = period == 1 ? this : MovingAverage(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.MovingSum)
            {
                smoothedSeries = period == 1 ? this : MovingSum(period);
            }
            else if (smoothingFunction == SmoothingFunctionType.Difference)
            {
                smoothedSeries = Difference(period);
            }


            var result = new TimeSeries(TimeInterval.Irregular);       
            for (int i = 0; i < Count; i++)
            {
                if (!double.IsNaN(smoothedSeries[i].Value) && smoothedSeries[i].Value > threshold)
                {
                    result.Add(smoothedSeries[i].Clone());

                    for (int j = i + 1; j < Count; j++)
                    {
                        if (!double.IsNaN(smoothedSeries[i].Value) && smoothedSeries[j].Value <= threshold && CheckIfMinStepsExceeded(result.Last().Index, smoothedSeries[j].Index, minStepsBetweenEvents))
                        {
                            i = j;
                            break;
                        }
                        if (!double.IsNaN(smoothedSeries[i].Value) && smoothedSeries[j].Value >= result.Last().Value)
                        {
                            result[result.Count - 1] = smoothedSeries[j].Clone();
                        }
                    }
                }
            }

            return result;
        }

        #endregion



        /// <summary>
        /// Download time series data from USGS
        /// </summary>
        /// <param name="siteNumber">USGS site number.</param>
        /// <param name="timeSeriesType">The time series type.</param>
        public static TimeSeries DownloadfromUSGS(string siteNumber, USGSTimeSeriesType timeSeriesType = USGSTimeSeriesType.DailyDischarge)
        {
            var timeSeries = new TimeSeries();

            try
            {
                // Check if there is an internet connection
                if (IsConnectedToInternet() == false)
                    throw new Exception("There is no internet connection!");

                // Setup url parameters
                string timeInterval = "dv";
                string startDate = "1800-01-01";
                string endDate = DateTime.Now.ToString("yyyy-MM-dd");
                string statCode = "&statCd=00003"; 
                string parameterCode = "00060";      
                string url = "";

                if (timeSeriesType == USGSTimeSeriesType.DailyDischarge || timeSeriesType == USGSTimeSeriesType.DailyStage)
                {
                    timeInterval = "dv";
                    statCode = "&statCd=00003";
                    parameterCode = timeSeriesType == USGSTimeSeriesType.DailyDischarge ? "00060" : "00065";
                    timeSeries = new TimeSeries(TimeInterval.OneDay);
                    url = "https://waterservices.usgs.gov/nwis/" + timeInterval + "/?format=waterml,2.0&sites=" + siteNumber + "&startDT=" + startDate + "&endDT=" + endDate + statCode + "&parameterCd=" + parameterCode + "&siteStatus=all";
                }
                else if (timeSeriesType == USGSTimeSeriesType.InstantaneousDischarge || timeSeriesType == USGSTimeSeriesType.InstantaneousStage)
                {
                    timeInterval = "iv";
                    statCode = "";
                    parameterCode = timeSeriesType == USGSTimeSeriesType.InstantaneousDischarge ? "00060" : "00065";
                    timeSeries = new TimeSeries(TimeInterval.FifteenMinute);
                    url = "https://waterservices.usgs.gov/nwis/" + timeInterval + "/?format=waterml,2.0&sites=" + siteNumber + "&startDT=" + startDate + "&endDT=" + endDate + statCode + "&parameterCd=" + parameterCode + "&siteStatus=all";
                }
                else if (timeSeriesType == USGSTimeSeriesType.PeakDischarge || timeSeriesType == USGSTimeSeriesType.PeakStage)
                {
                    timeInterval = "peak";
                    statCode = "";
                    timeSeries = new TimeSeries(TimeInterval.Irregular);
                    url = "https://nwis.waterdata.usgs.gov/nwis/peak?site_no="+ siteNumber + "&agency_cd=USGS&format=rdb";
                }
                            
                // Download data
                string textDownload;
                using (var client = new WebClient())
                    textDownload = client.DownloadString(url);

                if (timeSeriesType == USGSTimeSeriesType.DailyDischarge ||  timeSeriesType == USGSTimeSeriesType.DailyStage )
                {
                    // Convert to XElement and check if the download was valid
                    var xElement = XElement.Parse(textDownload);
                    var points = xElement.Descendants("{http://www.opengis.net/waterml/2.0}point");

                    // Create time series
                    foreach (XElement point in points.Elements())
                    {
                        if (point.Element("{http://www.opengis.net/waterml/2.0}time") != null && point.Element("{http://www.opengis.net/waterml/2.0}value") != null)
                        {
                            // Get date
                            DateTime index = DateTime.Now;
                            DateTime.TryParse(point.Element("{http://www.opengis.net/waterml/2.0}time").Value, out index);

                            // See if this date is 
                            if (timeSeries.Count > 0 && index != AddTimeInterval(timeSeries.Last().Index, TimeInterval.OneDay))
                            {
                                while (timeSeries.Last().Index < SubtractTimeInterval(index, TimeInterval.OneDay))
                                    timeSeries.Add(new SeriesOrdinate<DateTime, double>(AddTimeInterval(timeSeries.Last().Index, TimeInterval.OneDay), double.NaN));
                            }

                            // Get value
                            double value = 0;
                            string valueStg = point.Element("{http://www.opengis.net/waterml/2.0}value").Value;
                            if (valueStg == "" || valueStg == " " || valueStg == "  " || string.IsNullOrEmpty(valueStg))
                                value = double.NaN;
                            else
                            {
                                double.TryParse(valueStg, out value);
                            }
                            timeSeries.Add(new SeriesOrdinate<DateTime, double>(index, value));
                        }
                    }
                }
                else if (timeSeriesType == USGSTimeSeriesType.PeakDischarge || timeSeriesType == USGSTimeSeriesType.PeakStage)
                {
                    var delimiters = new char[] { '\t' };
                    var lines = textDownload.Split(new string[] { Environment.NewLine }, StringSplitOptions.RemoveEmptyEntries);        
                    foreach (string line in lines)
                    {
                        var segments = line.Split(delimiters, StringSplitOptions.RemoveEmptyEntries);
                        if (segments.Count() >= 4 && segments.First() == "USGS")
                        {
                            // Get date
                            DateTime index = DateTime.Now;
                            var dateString = segments[2].Split('-');
                            if (dateString[1] == "00")
                            {
                                int year = 2000;
                                int.TryParse(dateString[0], out year);
                                index = new DateTime(year, 1, 1, 0, 0, 0);
                            }
                            else
                            {
                                DateTime.TryParse(segments[2], out index);
                            }
                            

                            // Get value
                            double value = 0;

                            // see if the 4th column has a time
                            int offset = 0;
                            var timeString = segments[3].Split(':');
                            if (timeString.Length == 2)
                                offset = 1;

                            int idx = timeSeriesType == USGSTimeSeriesType.PeakDischarge ? 3 + offset : 4 + offset;
                            if (segments[idx] == "" || segments[idx] == " " || segments[idx] == "  " || string.IsNullOrEmpty(segments[idx]))
                                value = double.NaN;
                            else
                            {
                                double.TryParse(segments[idx], out value);
                            }
                            timeSeries.Add(new SeriesOrdinate<DateTime, double>(index, value));
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

        /// <summary>
        /// Checks if there is an Internet connection.
        /// </summary>
        public static bool IsConnectedToInternet()
        {
            try
            {
                using (var client = new WebClient())
                using (client.OpenRead("http://google.com/generate_204"))
                    return true;
            }
            catch { }
            return false;
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
                ordinate.SetAttributeValue("Index", this[i].Index.ToString());
                ordinate.SetAttributeValue("Value", this[i].Value.ToString());
                result.Add(ordinate);
            }
            return result;
        }

        /// <summary>
        /// Creates a copy of the time series.
        /// </summary>
        public TimeSeries Clone()
        {
            return new TimeSeries(TimeInterval, StartDate, ValuesToArray());          
        }
    }
}