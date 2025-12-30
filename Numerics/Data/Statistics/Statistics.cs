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
using System.Linq;
using System.Threading.Tasks;

namespace Numerics.Data.Statistics
{

    /// <summary>
    /// Contains functions for computing descriptive statistics of a sample of data.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// <see href = "https://en.wikipedia.org/wiki/Summary_statistics" />
    /// </item>
    /// <item>
    /// <see href = "https://en.wikipedia.org/wiki/Descriptive_statistics" />
    /// </item>
    /// <item>
    /// This class contains some functions from the Math.NET Numerics library, <see href="http://numerics.mathdotnet.com"/>
    /// </item>
    /// </list>
    /// </para>
    /// </remarks>
    public static class Statistics
    {
        /// <summary>
        /// Returns the smallest value from the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double Minimum(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            double min = double.PositiveInfinity;
            for (int i = 0; i < data.Count; i++)
            {
                if (double.IsNaN(data[i]))
                    return double.NaN;
                if (data[i] < min)
                    min = data[i];
            }

            return double.IsPositiveInfinity(min) ? double.NaN : min;
        }

        /// <summary>
        /// Returns the largest value from the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double Maximum(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            double max = double.NegativeInfinity;
            for (int i = 0; i < data.Count; i++)
            {
                if (double.IsNaN(data[i]))
                    return double.NaN;
                if (data[i] > max)
                    max = data[i];
            }

            return double.IsNegativeInfinity(max) ? double.NaN : max;
        }

        /// <summary>
        /// Estimates the sum of the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double Sum(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            double sum = 0;
            for (int i = 0; i < data.Count; i++)
                sum += data[i];
            return sum;
        }

        /// <summary>
        /// Estimates the arithmetic sample mean from the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double Mean(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;
            double sum = 0;
            for (int i = 0; i < data.Count; i++)
                sum += data[i];
            return sum / data.Count;
        }

        /// <summary>
        /// Computes the arithmetic sample mean from the unsorted data array by first enabling parallelization of the array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double ParallelMean(IList<double> data)
        {
            if (data.Count == 0) return double.NaN;

            double sum = data.AsParallel().Sum();
            return sum / data.Count;
        }

        /// <summary>
        /// Evaluates the geometric mean of the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double GeometricMean(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            double sum = 0;
            for (int i = 0; i < data.Count; i++)
            {
                if (data[i] <= 0) return double.NaN;
                sum += Math.Log(data[i]);
            }

            return Math.Exp(sum / data.Count);
        }

        /// <summary>
        /// Evaluates the harmonic mean of the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double HarmonicMean(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            double sum = 0;
            for (int i = 0; i < data.Count; i++)
                sum += 1.0 / data[i];
            return data.Count / sum;
        }

        /// <summary>
        /// Estimates the unbiased population variance from the provided samples as unsorted array.
        /// On a dataset of size N will use an N-1 normalizer (Bessel's correction).
        /// Returns NaN if data has less than two entries or if any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double Variance(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count <= 1) return double.NaN;

            double variance_ = 0;
            double t = data[0];
            for (int i = 1; i < data.Count; i++)
            {
                t += data[i];
                double diff = (i + 1) * data[i] - t;
                variance_ += diff * diff / ((i + 1.0d) * i);
            }
            return variance_ / (data.Count - 1);
        }

        /// <summary>
        /// Evaluates the population variance from the full population provided as unsorted array.
        /// On a dataset of size N will use an N normalizer and would thus be biased if applied to a subset.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double PopulationVariance(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            double variance = 0;
            double t = data[0];
            for (int i = 1; i < data.Count; i++)
            {
                t += data[i];
                double diff = (i + 1) * data[i] - t;
                variance += diff * diff / ((i + 1.0) * i);
            }
            return variance / data.Count;
        }

        /// <summary>
        /// Estimates the unbiased population standard deviation from the provided samples as unsorted array.
        /// On a dataset of size N will use an N-1 normalizer (Bessel's correction).
        /// Returns NaN if data has less than two entries or if any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double StandardDeviation(IList<double> data)
        {
            return Math.Sqrt(Variance(data));
        }

        /// <summary>
        /// Evaluates the population standard deviation from the full population provided as unsorted array.
        /// On a dataset of size N will use an N normalizer and would thus be biased if applied to a subset.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double PopulationStandardDeviation(IList<double> data)
        {
            return Math.Sqrt(PopulationVariance(data));
        }

        /// <summary>
        /// Estimates the arithmetic sample mean and the unbiased population variance from the provided samples as unsorted array.
        /// On a dataset of size N will use an N-1 normalizer (Bessel's correction).
        /// Returns NaN for mean if data is empty or any entry is NaN and NaN for variance if data has less than two entries or if any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static (double mean, double variance) MeanVariance(IList<double> data)
        {
            return (Mean(data), Variance(data));
        }

        /// <summary>
        /// Estimates the arithmetic sample mean and the unbiased population standard deviation from the provided samples as unsorted array.
        /// On a dataset of size N will use an N-1 normalizer (Bessel's correction).
        /// Returns NaN for mean if data is empty or any entry is NaN and NaN for standard deviation if data has less than two entries or if any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static (double mean, double standardDeviation) MeanStandardDeviation(IList<double> data)
        {
            return (Mean(data), StandardDeviation(data));
        }

        /// <summary>
        /// Estimates the coefficient of variation from the provided sample of data.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double CoefficientOfVariation(IList<double> data)
        {
            return StandardDeviation(data) / Mean(data);
        }

        /// <summary>
        /// Estimates the skewness coefficient from the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double Skewness(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            double mean = Mean(data);
            int n = data.Count;
            double s2 = 0, s3 = 0;
            for (int i = 0; i < n; i++)
            {
                double xm = data[i] - mean;
                s2 += xm * xm;
                s3 += xm * xm * xm;
            }
            double m2 = s2 / n;
            double m3 = s3 / n;
            double g = m3 / Math.Pow(m2, 3.0d / 2.0d);
            double a = Math.Sqrt(n * (n - 1.0));
            double b = n - 2;
            return a / b * g;
        }

        /// <summary>
        /// Computes the standard error of the statistic, using the jackknife method.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <param name="statistic">The statistic for estimating standard error.</param>
        public static double JackKnifeStandardError(IList<double> data, Func<IList<double>, double> statistic)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            int N = data.Count;
            double theta = statistic(data);
            double I = 0;
            Parallel.For(0, N, () => 0d, (i, loop, subI) =>
            {
                // Remove data point
                var jackSample = new List<double>(data);
                jackSample.RemoveAt(i);
                // Compute statistic
                subI += Tools.Sqr(statistic(jackSample) - theta);
                return subI;
            }, z => Tools.ParallelAdd(ref I, z));
            return Math.Sqrt((N - 1) / (double)N * I);
        }

        /// <summary>
        /// Returns a jackknifed sample.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <param name="statistic">The statistic for estimating a sample.</param>
        public static double[]? JackKnifeSample(IList<double> data, Func<IList<double>, double> statistic)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return null;

            int N = data.Count;
            var thetaJack = new double[N];
            // Perform Jackknife
            Parallel.For(0, N, i =>
            {
                // Remove data point
                var jackSample = new List<double>(data);
                jackSample.RemoveAt(i);
                // Compute statistic
                thetaJack[i] = statistic(jackSample);
            });
            return thetaJack;
        }

        /// <summary>
        /// Estimates the kurtosis from the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double Kurtosis(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;

            double mean = Mean(data);
            int n = data.Count;
            double s2 = 0, s4 = 0;
            for (int i = 0; i < n; i++)
            {
                double xm = data[i] - mean;
                s2 += xm * xm;
                s4 += xm * xm * xm * xm;
            }
            double m2 = s2 / n;
            double m4 = s4 / n;
            double v = s2 / (n - 1);
            double a = n * (n + 1) / (double)((n - 1) * (n - 2) * (n - 3));
            double b = s4 / (v * v);
            double c = (n - 1) * (n - 1) / (double)((n - 2) * (n - 3));
            return a * b - 3 * c;
        }

        /// <summary>
        /// Estimates the unbiased population covariance from the provided two sample arrays.
        /// On a dataset of size N will use an N-1 normalizer (Bessel's correction).
        /// Returns NaN if data has less than two entries or if any entry is NaN.
        /// </summary>
        /// <param name="data1">First sample of data, no sorting is assumed.</param>
        /// <param name="data2">Second sample of data, no sorting is assumed.</param>
        public static double Covariance(IList<double> data1, IList<double> data2)
        {
            if (data1 == null) throw new ArgumentNullException(nameof(data1));
            if (data2 == null) throw new ArgumentNullException(nameof(data2));
            if (data1.Count != data2.Count)
            {
                throw new ArgumentException("All vectors must have the same dimensionality.");
            }

            if (data1.Count <= 1) return double.NaN;

            double mean1 = Mean(data1);
            double mean2 = Mean(data2);
            double covariance = 0.0;
            for (int i = 0; i < data1.Count; i++)
                covariance += (data1[i] - mean1) * (data2[i] - mean2);
            return covariance / (data1.Count - 1);
        }

        /// <summary>
        /// Evaluates the population covariance from the full population provided as two arrays.
        /// On a dataset of size N will use an N normalizer and would thus be biased if applied to a subset.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        /// <param name="data1">First sample of data, no sorting is assumed.</param>
        /// <param name="data2">Second sample of data, no sorting is assumed.</param>
        public static double PopulationCovariance(IList<double> data1, IList<double> data2)
        {
            if (data1 == null) throw new ArgumentNullException(nameof(data1));
            if (data2 == null) throw new ArgumentNullException(nameof(data2));
            if (data1.Count != data2.Count)
            {
                throw new ArgumentException("All vectors must have the same dimensionality.");
            }

            if (data1.Count == 0) return double.NaN;

            double mean1 = Mean(data1);
            double mean2 = Mean(data2);
            double covariance = 0.0;
            for (int i = 0; i < data1.Count; i++)
                covariance += (data1[i] - mean1) * (data2[i] - mean2);
            return covariance / data1.Count;
        }

        /// <summary>
        /// Returns the first four product moments of a sample {Mean, Standard Deviation, Skew, and Kurtosis}, or returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double[] ProductMoments(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            double N = data.Count;
            if (N < 4) return [double.NaN, double.NaN, double.NaN, double.NaN];

            // sums of powers
            double X1 = 0, X2 = 0, X3 = 0, X4 = 0;
            foreach (var x in data)
            {
                double x2 = x * x;
                X1 += x;
                X2 += x2;
                X3 += x2 * x;
                X4 += x2 * x2;
            }

            // raw moments
            double U1 = X1 / N;
            double U2 = X2 / N;
            double U3 = X3 / N;
            double U4 = X4 / N;

            // central moments
            double m2 = (U2 - U1 * U1) * (N / (N - 1));  // sample variance
            double S = Math.Sqrt(m2);

            // pre-compute powers
            double U1_2 = U1 * U1;
            double U1_3 = U1_2 * U1;
            double U1_4 = U1_3 * U1;
            double S3 = S * S * S;
            double S4 = S3 * S;

            // third central moment
            double c3 = U3 - 3 * U1 * U2 + 2 * U1_3;
            // fourth central moment
            double c4 = U4 - 4 * U1 * U3 + 6 * U2 * U1_2 - 3 * U1_4;

            // bias-corrected skewness
            double G = (N * N) / ((N - 1) * (N - 2)) * (c3 / S3);

            // bias-corrected excess kurtosis
            double K = ((N * N) * (N + 1)) / ((N - 1) * (N - 2) * (N - 3)) * (c4 / S4) - 3d * (N - 1) * (N - 1) / ((N - 2) * (N - 3));

            return [U1, S, G, K];

        }

        /// <summary>
        /// Returns the linear moments of a sample {L-Mean (λ1), L-Scale (λ2), L-Skewness (τ3), and L-Kurtosis (τ4)}, or returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double[] LinearMoments(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            double N = data.Count;
            if (N < 4) return [double.NaN, double.NaN, double.NaN, double.NaN];

            // Copy and sort data
            var sortedData = data.ToArray();
            Array.Sort(sortedData);

            double B0 = 0, B1 = 0, B2 = 0, B3 = 0;
            for (int i = 1; i <= N; i++)
            {
                B0 += sortedData[i - 1];
                if (i > 1)
                    B1 += (i - 1) / (N - 1) * sortedData[i - 1];
                if (i > 2)
                    B2 += (i - 2) * (i - 1) / ((N - 2) * (N - 1)) * sortedData[i - 1];
                if (i > 3)
                    B3 += (i - 3) * (i - 2) * (i - 1) / ((N - 3) * (N - 2) * (N - 1)) * sortedData[i - 1];
            }
 
            B0 /= N;
            B1 /= N;
            B2 /= N;
            B3 /= N;
            // L-Mean (λ1)
            // L-Scale (λ2)
            // L-Skewness (τ3)
            // L-Kurtosis (τ4)
            double L1 = B0;
            double L2 = 2 * B1 - B0;
            double T3 = 2 * (3 * B2 - B0) / (2 * B1 - B0) - 3;
            double T4 = 5 * (2 * (2 * B3 - 3 * B2) + B0) / (2 * B1 - B0) + 6;
            return [L1, L2, T3, T4];
        }

        /// <summary>
        /// Returns the k-th percentile of values in a sample.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="k">The k-th percentile to find.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The k-th percentile.</returns>
        public static double Percentile(IList<double> data, double k, bool dataIsSorted = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            int n = data.Count;
            if (n == 0) throw new ArgumentException("Sequence contains no elements.", nameof(data));
            if (k < 0.0 || k > 1.0) throw new ArgumentOutOfRangeException(nameof(k), "k must be in [0,1].");

            // Copy & sort if needed
            var sortedData = dataIsSorted ? data: data.OrderBy(x => x).ToArray();

            // Trivial cases
            if (n == 1 || k == 0.0) return sortedData[0];
            if (k == 1.0) return sortedData[n - 1];

            // Zero-based linear interpolation (Type 7)
            double h = (n - 1) * k;
            int lower = (int)Math.Floor(h);
            int upper = (int)Math.Ceiling(h);
            double w = h - lower;
            return sortedData[lower] + w * (sortedData[upper] - sortedData[lower]);
        }

        /// <summary>
        /// Returns an array of percentile values of a sample.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="k">The list of k-th percentiles to find.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The k-th percentile.</returns>
        public static double[] Percentile(IList<double> data, IList<double> k, bool dataIsSorted = false)
        {
            // Copy & sort if needed
            var sortedData = dataIsSorted ? data : data.OrderBy(x => x).ToArray();
            var result = new double[k.Count];
            for (int i = 0; i < k.Count; i++)
            {
                result[i] = Percentile(sortedData, k[i], true);
            }
            return result;
        }

        /// <summary>
        /// Estimates the 5-number summary {min, 25th-percentile, 50th-percentile, 75th-percentile, max} from a sample of data.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <returns>5-number summary statistics.</returns>
        public static double[] FiveNumberSummary(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            // Copy & sort
            var sortedData = data.ToArray();
            Array.Sort(sortedData);
            double min = sortedData[0];
            double max = sortedData[sortedData.Count() - 1];
            double p25 = Percentile(sortedData, 0.25, true);
            double p50 = Percentile(sortedData, 0.5, true);
            double p75 = Percentile(sortedData, 0.75, true);
            return [min, p25, p50, p75, max];
        }

        /// <summary>
        /// Estimates the 7-number summary {min, 5th percentile, 25th-percentile, 50th-percentile, 75th-percentile, 95th-percentile, max} from a sample of data.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <returns>7-number summary statistics.</returns>
        public static double[] SevenNumberSummary(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            // Copy & sort
            var sortedData = data.ToArray();
            Array.Sort(sortedData);
            double min = sortedData[0];
            double max = sortedData[sortedData.Count() - 1];
            double p5 = Percentile(sortedData, 0.05, true);
            double p25 = Percentile(sortedData, 0.25, true);
            double p50 = Percentile(sortedData, 0.5, true);
            double p75 = Percentile(sortedData, 0.75, true);
            double p95 = Percentile(sortedData, 0.95, true);
            return [min, p5, p25, p50, p75, p95, max];
        }

        /// <summary>
        /// Returns the rank of each entry of the unsorted data array.
        /// </summary>
        /// <param name="data">The array of sample of data, no sorting is assumed.</param>
        public static double[] RanksInPlace(double[] data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));

            var ranks = new double[data.Length];
            var index = new int[data.Length];
            for (int i = 0; i < index.Length; i++)
            {
                index[i] = i;
            }

            // Copy and sort array
            var work = (double[])data.Clone();
            Array.Sort(work, index);

            int previousIndex = 0;
            for (int i = 1; i < work.Length; i++)
            {

                if (Math.Abs(work[i] - work[previousIndex]) <= 0)
                {
                    continue;
                }

                if (i == previousIndex + 1)
                {
                    ranks[index[previousIndex]] = i;
                }
                else
                {
                    RanksTies(ranks, index, previousIndex, i);
                }

                previousIndex = i;
            }

            RanksTies(ranks, index, previousIndex, work.Length);
            return ranks;
        }

        /// <summary>
        /// Returns the rank of each entry of the unsorted data array.
        /// </summary>
        /// <param name="data">The array of sample of data, no sorting is assumed.</param>
        /// <param name="ties">Output. The number of ties in the data.</param>
        public static double[] RanksInPlace(double[] data, out double [] ties)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));

            var ranks = new double[data.Length];
            ties = new double[data.Length];
            var index = new int[data.Length];
            for (int i = 0; i < index.Length; i++)
            {
                index[i] = i;
            }

            // Copy and sort array
            var work = (double[])data.Clone();
            Array.Sort(work, index);

            int previousIndex = 0;
            int t = 0;
            for (int i = 1; i < work.Length; i++)
            {
                if (work[i].AlmostEquals(work[previousIndex], Tools.DoubleMachineEpsilon))
                {
                    t += 1;
                    continue;
                }

                if (i == previousIndex + 1)
                {
                    ranks[index[previousIndex]] = i;
                    t = 0;      
                }
                else
                {
                    RanksTies(ranks, index, previousIndex, i);
                    ties[i - 1] = t;
                    t = 0;
                }

                previousIndex = i;
            }

            RanksTies(ranks, index, previousIndex, work.Length);
            return ranks;
        }

        /// <summary>
        /// Helper function for RanksInplace(double[], out double[])
        /// </summary>
        private static void RanksTies(double[] ranks, int[] index, int a, int b)
        {       
            double rank = (b + a - 1) / 2d + 1;
            for (int k = a; k < b; k++)
            {
                ranks[index[k]] = rank;
            }
        }

        /// <summary>
        /// Computes the entropy function for a set of numerical values in a given Probability Density Function (pdf).
        /// </summary>
        /// <param name="data">The array of values.</param>
        /// <param name="pdf">A probability distribution function.</param>
        public static double Entropy(double[] data, Func<double, double> pdf)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));

            double sum = 0;
            for (int i = 0; i < data.Length; i++)
            {
                double p = pdf(data[i]);
                if (p > 0)
                {
                    sum += p * Math.Log(p);
                }               
            }
            return -sum;
        }
    }
}