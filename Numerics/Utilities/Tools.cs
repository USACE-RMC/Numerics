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
using System;
using System.Collections.Generic; 
using System.IO;
using System.IO.Compression;
using System.Threading;

namespace Numerics
{
    /// <summary>
    /// A class of public utility functions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public static class Tools
    {

        /// <summary>
        /// Value used for comparing double precision equivalence. Represented as 2^(-53) = 1.11022302462516E-16.
        /// </summary>
        public const double DoubleMachineEpsilon = 1.11022302462516E-16;

        /// <summary>
        /// Value used for comparing single precision equivalence. Represented as 2^(-24) = 5.96046447753906E-07
        /// </summary>
        public const float SingleMachineEpsilon = (float)5.96046447753906E-07;

        /// <summary>
        /// Euler constant.
        /// </summary>
        public const double Euler = 0.5772156649015328606065120;

        /// <summary>
        /// The golden ratio. 
        /// </summary>
        public static double GoldenRatio = 1.61803398874989484820; // 0.5d * (1d + Math.Sqrt(5));

        /// <summary>
        /// Log of square root of 2*pi
        /// </summary>
        public const double LogSqrt2PI = 0.91893853320467274178032973640562;

        /// <summary>
        /// Square root of 2*pi
        /// </summary>
        public const double Sqrt2PI = 2.50662827463100050242E0;

        /// <summary>
        /// Square root of 2
        /// </summary>
        public const double Sqrt2 = 1.4142135623730950488016887;

        /// <summary>
        /// Log of 2.
        /// </summary>
        public const double Log2 = 0.69314718055994530941;

        /// <summary>
        /// Returns a value with the same magnitude as a and the same sign as b.
        /// </summary>
        /// <param name="a">Value a.</param>
        /// <param name="b">Value b.</param>
        public static double Sign(double a, double b)
        {
            return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
        }

        /// <summary>
        /// Swap value a with value b.
        /// </summary>
        /// <param name="a">Value a.</param>
        /// <param name="b">Value b.</param>
        public static void Swap(ref double a, ref double b)
        {
            double temp = a;
            a = b;
            b = temp;
        }

        /// <summary>
        /// Returns the squared value of a.
        /// </summary>
        /// <param name="a">Value a.</param>
        public static double Sqr(double a)
        {
            return a * a;
        }

        /// <summary>
        /// Returns the value of a raised to the power of b.
        /// </summary>
        /// <param name="a">The value a.</param>
        /// <param name="b">The exponent b.</param>
        public static double Pow(double a, int b)
        {
            if (b == 0d) return 1d;
            double e = a; int n = Math.Abs(b);
            for (int i = 2; i <= n; i++) e *= a;
            if (b > 0d) return e;
            return 1d / e;
        }

        /// <summary>
        /// Returns the base 10 logarithm of a specified number. 
        /// </summary>
        /// <param name="x">The number whose logarithm is to be found.</param>
        public static double Log10(double x)
        {
            if (x < 1E-16 && Math.Sign(x) != -1) x = 1E-16;
            return Math.Log10(x);
        }

        /// <summary>
        /// Returns the natural (base e) logarithm of a specified number. 
        /// </summary>
        /// <param name="x">The number whose logarithm is to be found.</param>
        public static double Log(double x)
        {
            if (x < 1E-16 && Math.Sign(x) != -1) x = 1E-16;
            return Math.Log(x);
        }

        /// <summary>
        /// Returns the Euclidean distance between two points ||x - y||.
        /// </summary>
        /// <param name="x1">X of point 1.</param>
        /// <param name="y1">Y of point 1.</param>
        /// <param name="x2">X of point 2.</param>
        /// <param name="y2">Y of point 2.</param>
        public static double Distance(double x1, double y1, double x2, double y2)
        {       
            return Math.Sqrt(Math.Pow(x2 - x1, 2d) + Math.Pow(y2 - y1, 2d));
        }

        /// <summary>
        /// Returns the Euclidean distance between two points ||x - y||.
        /// </summary>
        /// <param name="x">First point.</param>
        /// <param name="y">Second point.</param>
        /// <returns>The Euclidean distance ||x - y||.</returns>
        public static double Distance(IList<double> x, IList<double> y)
        {
            double d = 0;
            for (int i = 0; i < x.Count; i++)
            {
                double dx = x[i] - y[i];
                d += dx * dx;
            }
            return Math.Sqrt(d);
        }

        /// <summary>
        /// Returns the Euclidean distance between two points after applying min-max normalization to each dimension.
        /// Normalization ensures each parameter dimension contributes equally to the distance metric,
        /// based on the specified lower and upper bounds.
        /// </summary>
        /// <param name="x">First point.</param>
        /// <param name="y">Second point.</param>
        /// <param name="lower">The lower bounds for each parameter dimension.</param>
        /// <param name="upper">The upper bounds for each parameter dimension.</param>
        /// <returns>The normalized Euclidean distance ||x - y||.</returns>
        public static double NormalizedDistance(IList<double> x, IList<double> y, IList<double> lower, IList<double> upper)
        {
            double d = 0;
            for (int i = 0; i < x.Count; i++)
            {
                double range = upper[i] - lower[i];
                double xi = (x[i] - lower[i]) / range;
                double yi = (y[i] - lower[i]) / range;
                double dx = xi - yi;
                d += dx * dx;
            }
            return Math.Sqrt(d);
        }


        /// <summary>
        /// Returns the weighted Euclidean distance between two points ||x - y||.
        /// </summary>
        /// <param name="x">First point.</param>
        /// <param name="y">Second point.</param>
        /// <param name="weights">Weights.</param>
        /// <returns>The weighted Euclidean distance ||x - y||.</returns>
        public static double WeightedDistance(IList<double> x, IList<double> y, IList<double> weights)
        {
            double d = 0;
            for (int i = 0; i < x.Count; i++)
            {
                double dx = x[i] - y[i];
                d += (dx * dx) * weights[i];
            }
            return Math.Sqrt(d);
        }

        /// <summary>
        /// Returns the normalized values ranging between 0 and 1.
        /// <para> 
        /// <see href="https://www.codecademy.com/article/normalization"/>
        /// </para>
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double[] Normalize(IList<double> values)
        {
            var result = new double[values.Count];
            double min = 0, max = 0;
            MinMax(values, out min, out max);
            for (int i = 0; i < values.Count; i++)
                result[i] = (values[i] - min) / (max - min);
            return result;
        }

        /// <summary>
        /// Returns the denormalized values. 
        /// </summary>
        /// <param name="values">The list of normalized values to denormalize.</param>
        /// <param name="min">The minimum value from the original data.</param>
        /// <param name="max">The maximum value from the original data.</param>
        public static double[] Denormalize(IList<double> values, double min, double max)
        {
            var result = new double[values.Count];
            for (int i = 0; i < values.Count; i++)
                result[i] = values[i] * (max - min) + min;
            return result;
        }

        /// <summary>
        /// Returns the standardized values. 
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double[] Standardize(IList<double> values)
        {
            var result = new double[values.Count];
            var meanSD = Statistics.MeanStandardDeviation(values);
            for (int i = 0; i < values.Count; i++)
                result[i] = (values[i] - meanSD.Item1) / meanSD.Item2;
            return result;
        }

        /// <summary>
        /// Returns the destandardized values. 
        /// </summary>
        /// <param name="values">The list of standardized values to destandardize.</param>
        /// <param name="mean">The mean of the original data.</param>
        /// <param name="standardDeviation">The standard deviation of the original data.</param>
        public static double[] Destandardize(IList<double> values, double mean, double standardDeviation)
        {
            var result = new double[values.Count];
            for (int i = 0; i < values.Count; i++)
                result[i] = values[i] * standardDeviation + mean;
            return result;
        }

        /// <summary>
        /// Estimates the sum of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Sum(IList<int> values)
        {
            if (values.Count == 0) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values.Count; i++)
                sum += values[i];
            return sum;
        }

        /// <summary>
        /// Estimates the sum of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Sum(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values.Count; i++)
                sum += values[i];
            return sum;
        }

        /// <summary>
        /// Estimate the sum of a list of values. 
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param>
        public static double Sum(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values.Count; i++)
                if (indicators[i] == (useComplement ? 0 : 1)) sum += values[i];      
            return sum;
        }

        /// <summary>
        /// Returns the sum product of two lists of values.
        /// </summary>
        /// <param name="values1">The first list of values.</param>
        /// <param name="values2">The second list of values.</param>
        public static double SumProduct(IList<double> values1, IList<double> values2)
        {
            if (values1.Count == 0) return double.NaN;
            if (values2.Count != values1.Count) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values1.Count; i++)
                sum += values1[i] * values2[i];
            return sum;
        }


        /// <summary>
        /// Estimates the mean of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Mean(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values.Count; i++)
                sum += values[i];
            return sum / values.Count;
        }

        /// <summary>
        /// Estimates the mean of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param>
        public static double Mean(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double sum = 0d;
            int count = 0;
            for (int i = 0; i < values.Count; i++)
            {
                if (indicators[i] == (useComplement ? 0 : 1))
                {
                    sum += values[i];
                    count += 1;
                }
            }
            return sum == 0 || count == 0 ? 0 : sum / count;
        }

        /// <summary>
        /// Estimates the product of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Product(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double product = values[0];
            for (int i = 1; i < values.Count; i++)
                product *= values[i];
            return product;
        }

        /// <summary>
        /// Estimates the product of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param> 
        public static double Product(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            int i, j, idx = 0;
            double product = 0;
            for (i = 0; i < values.Count; i++)
            {
                if (indicators[i] == (useComplement ? 0 : 1))
                {
                    product = values[i];
                    idx = i;
                    break;
                }
            }
            for (j = idx + 1; j < values.Count; j++)
                if (indicators[j] == (useComplement ? 0 : 1)) product *= values[j];
            return product;
        }

        /// <summary>
        /// Returns the minimum and maximum values from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="min">Output. Minimum value.</param>
        /// <param name="max">Output. Maximum value.</param>
        public static void MinMax(IList<double> values, out double min, out double max)
        {
            min = double.MaxValue;
            max = double.MinValue;
            if (values.Count == 0) return;
            for (int i = 0; i < values.Count; i++)
            {
                if (values[i] < min || double.IsNaN(values[i]))
                    min = values[i];

                if (values[i] > max || double.IsNaN(values[i]))
                    max = values[i];
            }

        }

        /// <summary>
        /// Returns the index of the minimum value.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static int ArgMin(IList<double> values)
        {
            double min = double.MaxValue;
            int index = -1;
            for (int i = 0; i < values.Count; i++)
            {
                if (values[i] < min)
                {
                    min = values[i];
                    index = i;
                }
            }
            return index;
        }

        /// <summary>
        /// Returns the index of the maximum value.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static int ArgMax(IList<double> values)
        {
            double max = double.MinValue;
            int index = -1;
            for (int i = 0; i < values.Count; i++)
            {
                if (values[i] > max)
                {
                    max = values[i];
                    index = i;
                }
            }
            return index;
        }

        /// <summary>
        /// Returns the smallest value from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Min(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double min = double.MaxValue;
            for (int i = 0; i < values.Count; i++)
                if (values[i] < min || double.IsNaN(values[i]))
                    min = values[i];
            return min;
        }

        /// <summary>
        /// Returns the smallest value from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param> 
        public static double Min(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double min = double.MaxValue;
            for (int i = 0; i < values.Count; i++)
                if (indicators[i] == (useComplement ? 0 : 1) && values[i] < min || double.IsNaN(values[i]))
                    min = values[i];
            return min;
        }

        /// <summary>
        /// Returns the largest value from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Max(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double max = double.MinValue;
            for (int i = 0; i < values.Count; i++)
                if (values[i] > max || double.IsNaN(values[i]))
                    max = values[i];
            return max;
        }

        /// <summary>
        /// Returns the largest value from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param> 
        public static double Max(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double max = double.MinValue;
            for (int i = 0; i < values.Count; i++)
                if(indicators[i] == (useComplement ? 0 : 1) && values[i] > max|| double.IsNaN(values[i]))
                    max = values[i];
            return max;
        }

        /// <summary>
        /// Supporting function used to add doubles from an interlocked parallel loop.
        /// </summary>
        /// <param name="valueToAddTo">The value to add to.</param>
        /// <param name="valueToAdd">The value to add.</param>
        public static double ParallelAdd(ref double valueToAddTo, double valueToAdd)
        {
            double newCurrentValue = valueToAddTo;
            while (true)
            {
                double currentValue = newCurrentValue;
                double newValue = currentValue + valueToAdd;
                newCurrentValue = Interlocked.CompareExchange(ref valueToAddTo, newValue, currentValue);
                if (newCurrentValue == currentValue) return newValue;
                if (double.IsNaN(newCurrentValue) && double.IsNaN(currentValue)) return newValue;
            }
        }

        /// <summary>
        /// The log-sum-exponential function.
        /// </summary>
        /// <param name="u">Log of u.</param>
        /// <param name="v">Log of v.</param>
        public static double LogSumExp(double u, double v)
        {
            double max = Math.Max(u, v);
            return max + Math.Log(Math.Exp(u - max) + Math.Exp(v - max));
        }

        /// <summary>
        /// The log-sum-exponential function. 
        /// </summary>
        /// <param name="values">The list of values of log.</param>
        public static double LogSumExp(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double max = Max(values);
            double exp = 0;
            for (int i = 0; i < values.Count; i++)
            {
                exp += Math.Exp(values[i] - max);
            }
            return max + Math.Log(exp);
        }

        /// <summary>
        /// Returns a sequence of integers.
        /// </summary>
        /// <param name="start">The starting value.</param>
        /// <param name="end">The ending value.</param>
        /// <param name="step">The step size.</param>
        public static int[] Sequence(int start, int end, int step = 1)
        {
            var result = new List<int>();
            int val = start;
            while (val <= end)
            {
                result.Add(val);
                val += step;
            };
            return result.ToArray();
        }

        /// <summary>
        /// Returns a sequence of doubles.
        /// </summary>
        /// <param name="start">The starting value.</param>
        /// <param name="end">The ending value.</param>
        /// <param name="step">The step size.</param>
        public static double[] Sequence(double start, double end, double step = 1)
        {
            var result = new List<double>();
            double val = start;
            while (val <= end)
            {
                result.Add(val);
                val += step;
            };
            return result.ToArray();
        }


        /// <summary>
        /// Returns a compressed a byte array.
        /// </summary>
        /// <param name="data">An array of bytes.</param>
        public static byte[] Compress(byte[] data)
        {
            if (data is null) return null;
            var output = new MemoryStream();
            using (var dstream = new DeflateStream(output, CompressionLevel.Optimal))
            {
                dstream.Write(data, 0, data.Length);
            }
            return output.ToArray();
        }

        /// <summary>
        /// Returns a decompressed byte array.
        /// </summary>
        /// <param name="data">An array of bytes.</param>
        public static byte[] Decompress(byte[] data)
        {
            if (data is null) return null;
            var input = new MemoryStream(data);
            var output = new MemoryStream();
            using (var dstream = new DeflateStream(input, CompressionMode.Decompress))
            {
                dstream.CopyTo(output);
            }
            return output.ToArray();
        }

    }
}
