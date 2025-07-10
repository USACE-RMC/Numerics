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

// <copyright file="Statistics.cs" company="Math.NET">
// Math.NET Numerics, part of the Math.NET Project
// http://numerics.mathdotnet.com
// http://github.com/mathnet/mathnet-numerics
// 
// Copyright (c) 2009-2015 Math.NET
// 
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// </copyright>

using System;
using System.Collections.Generic;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Contains functions for computing running statistics of a sample of data.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// This class is copied from the Math.NET Numerics library, <see href="http://numerics.mathdotnet.com"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public class RunningStatistics
    {

        /// <summary>
        /// Constructs an empty running statistics class.
        /// </summary>
        public RunningStatistics()
        {
        }

        /// <summary>
        /// Constructs running statistics based on a list of values.
        /// </summary>
        /// <param name="values">List of data values.</param>
        public RunningStatistics(IList<double> values)
        {
            Push(values);
        }

        private long _n;
        private double _min = double.PositiveInfinity;
        private double _max = double.NegativeInfinity;
        private double _m1;
        private double _m2;
        private double _m3;
        private double _m4;

        /// <summary>
        /// Gets the total number of samples.
        /// </summary>
        public long Count
        {
            get
            {
                return _n;
            }
        }

        /// <summary>
        /// Returns the minimum value in the sample data.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        public double Minimum
        {
            get
            {
                return _n > 0L ? _min : double.NaN;
            }
        }

        /// <summary>
        /// Returns the maximum value in the sample data.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        public double Maximum
        {
            get
            {
                return _n > 0L ? _max : double.NaN;
            }
        }

        /// <summary>
        /// Returns the sample mean, an estimate of the population mean.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        public double Mean
        {
            get
            {
                return _n > 0L ? _m1 : double.NaN;
            }
        }

        /// <summary>
        /// Returns the sample variance.
        /// On a dataset of size N, it will use an N-1 normalizer (Bessel's correction).
        /// Returns NaN if data has less than two entries or if any entry is NaN.
        /// </summary>
        public double Variance
        {
            get
            {
                return _n < 2L ? double.NaN : _m2 / (_n - 1L);
            }
        }

        /// <summary>
        /// Returns the population variance.
        /// On a dataset of size N, it will use an N normalizer and would thus be biased if applied to a subset.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        public double PopulationVariance
        {
            get
            {
                return _n < 2L ? double.NaN : _m2 / _n;
            }
        }

        /// <summary>
        /// Returns the sample standard deviation.
        /// On a dataset of size N, it will use an N-1 normalizer (Bessel's correction).
        /// Returns NaN if data has less than two entries or if any entry is NaN.
        /// </summary>
        public double StandardDeviation
        {
            get
            {
                return _n < 2L ? double.NaN : Math.Sqrt(_m2 / (_n - 1L));
            }
        }

        /// <summary>
        /// Returns the population standard deviation.
        /// On a dataset of size N, it will use an N normalizer and would thus be biased if applied to a subset.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        public double PopulationStandardDeviation
        {
            get
            {
                return _n < 2L ? double.NaN : Math.Sqrt(_m2 / _n);
            }
        }

        /// <summary>
        /// Returns the coefficient of variation of the sample.
        /// </summary>
        public double CoefficientOfVariation
        {
            get
            {
                return StandardDeviation / Mean;
            }
        }

        /// <summary>
        /// Returns the sample skewness.
        /// Uses a normalizer (Bessel's correction; type 2).
        /// Returns NaN if data has less than three entries or if any entry is NaN.
        /// </summary>
        public double Skewness
        {
            get
            {
                return _n < 3L ? double.NaN : _n * _m3 * Math.Sqrt(_m2 / (_n - 1L)) / (_m2 * _m2 * (_n - 2L)) * (_n - 1L);
            }
        }

        /// <summary>
        /// Returns the population skewness.
        /// Does not use a normalizer and would thus be biased if applied to a subset (type 1).
        /// Returns NaN if data has less than two entries or if any entry is NaN.
        /// </summary>
        public double PopulationSkewness
        {
            get
            {
                return _n < 2L ? double.NaN : Math.Sqrt(_n) * _m3 / Math.Pow(_m2, 1.5d);
            }
        }

        /// <summary>
        /// Returns the sample kurtosis.
        /// Uses a normalizer (Bessel's correction; type 2).
        /// Returns NaN if data has less than four entries or if any entry is NaN.
        /// </summary>
        public double Kurtosis
        {
            get
            {
                return _n < 4L ? double.NaN : (_n * (double)_n - 1d) / ((_n - 2L) * (_n - 3L)) * (_n * _m4 / (_m2 * _m2) - 3d + 6.0d / (_n + 1L));
            }
        }

        /// <summary>
        /// Returns the population kurtosis.
        /// Does not use a normalizer and would thus be biased if applied to a subset (type 1).
        /// Returns NaN if data has less than three entries or if any entry is NaN.
        /// </summary>
        public double PopulationKurtosis
        {
            get
            {
                return _n < 3L ? double.NaN : _n * _m4 / (_m2 * _m2) - 3.0d;
            }
        }

        /// <summary>
        /// Updates the running statistics by adding another data value (in-place).
        /// </summary>
        /// <param name="value">Data value.</param>
        public void Push(double value)
        {
            _n += 1L;
            // Update moments
            double d = value - _m1;
            double s = d / _n;
            double s2 = s * s;
            double t = d * s * (_n - 1L);
            _m1 += s;
            _m4 += t * s2 * (_n * _n - 3L * _n + 3L) + 6d * s2 * _m2 - 4d * s * _m3;
            _m3 += t * s * (_n - 2L) - 3d * s * _m2;
            _m2 += t;
            // Update min
            if (value < _min || double.IsNaN(value))
            {
                _min = value;
            }
            // Update max
            if (value > _max || double.IsNaN(value))
            {
                _max = value;
            }
        }

        /// <summary>
        /// Updates the running statistics by adding a sequence of data values (in-place).
        /// </summary>
        /// <param name="values">List of data values.</param>
        public void Push(IList<double> values)
        {
            foreach (double value in values)
                Push(value);
        }

        /// <summary>
        /// Create a new running statistics over the combined samples of two existing running statistics.
        /// </summary>
        public static RunningStatistics Combine(RunningStatistics a, RunningStatistics b)
        {
            if (a._n == 0L)
            {
                return b;
            }
            else if (b._n == 0L)
            {
                return a;
            }

            long n = a._n + b._n;
            double d = b._m1 - a._m1;
            double d2 = d * d;
            double d3 = d2 * d;
            double d4 = d2 * d2;
            double m1 = (a._n * a._m1 + b._n * b._m1) / n;
            double m2 = a._m2 + b._m2 + d2 * a._n * b._n / n;
            double m3 = a._m3 + b._m3 + d3 * a._n * b._n * (a._n - b._n) / (n * n) + 3d * d * (a._n * b._m2 - b._n * a._m2) / n;
            double m4 = a._m4 + b._m4 + d4 * a._n * b._n * (a._n * a._n - a._n * b._n + b._n * b._n) / (n * n * n) + 6d * d2 * (a._n * a._n * b._m2 + b._n * b._n * a._m2) / (n * n) + 4d * d * (a._n * b._m3 - b._n * a._m3) / n;
            double min = Math.Min(a._min, b._min);
            double max = Math.Max(a._max, b._max);
            return new RunningStatistics()
            {
                _n = n,
                _m1 = m1,
                _m2 = m2,
                _m3 = m3,
                _m4 = m4,
                _min = min,
                _max = max
            };
        }

        /// <summary>
        /// Create a new running statistics over the combined samples of two existing running statistics.
        /// </summary>
        public static RunningStatistics operator +(RunningStatistics a, RunningStatistics b)
        {
            return Combine(a, b);
        }

    }
}