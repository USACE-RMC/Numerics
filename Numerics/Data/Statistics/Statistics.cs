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
        /// The number of accumulation chunks used by the parallel mean reduction. Fixed, not derived
        /// from the processor count, so the summation order — and therefore the mean's last bits —
        /// does not vary with the machine or the thread count. Matches the bootstrap's jackknife
        /// reduction.
        /// </summary>
        private const int ParallelMeanChunks = 64;

        /// <summary>
        /// The sample size below which the parallel mean falls through to the sequential mean:
        /// scheduling a parallel loop costs more than summing this few values, and the sequential
        /// result is then bitwise identical to <see cref="Mean(IList{double})"/>.
        /// </summary>
        private const int ParallelMeanSequentialThreshold = 8192;

        /// <summary>
        /// Computes the arithmetic sample mean from the unsorted data array using a parallel
        /// reduction with a deterministic summation order.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <remarks>
        /// Large samples are split into a fixed number of chunks with balanced ranges, each chunk is
        /// summed sequentially into its own slot in parallel, and the chunk sums are combined
        /// serially in chunk order — the same deterministic reduction the bootstrap's jackknife
        /// accumulation uses. The summation tree therefore depends only on the sample length, so the
        /// result is bit-identical on every machine, core count, and scheduler. A processor-count
        /// partition would produce different last bits on different machines, and a shared
        /// accumulator such as <c>Tools.ParallelAdd</c> would be race-free but commit its additions
        /// in thread-scheduler order — the same non-reproducibility — so both are excluded from this
        /// reduction. Samples below the threshold fall through to the sequential
        /// <see cref="Mean(IList{double})"/>.
        /// </remarks>
        public static double ParallelMean(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            int n = data.Count;
            if (n == 0) return double.NaN;
            if (n < ParallelMeanSequentialThreshold) return Mean(data);

            int chunks = Math.Min(ParallelMeanChunks, n);
            var chunkSums = new double[chunks];
            Parallel.For(0, chunks, c =>
            {
                int start = (int)((long)c * n / chunks);
                int end = (int)((long)(c + 1) * n / chunks);
                double sum = 0d;
                for (int i = start; i < end; i++)
                    sum += data[i];
                chunkSums[c] = sum;
            });

            double total = 0d;
            for (int c = 0; c < chunks; c++)
                total += chunkSums[c];
            return total / n;
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
            {
                if (data[i] <= 0) return double.NaN;
                sum += 1.0 / data[i];
            }
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
        /// Computes the jackknife standard error of a statistic.
        /// </summary>
        /// <param name="data">Sample data.</param>
        /// <param name="statistic">The statistic evaluated on isolated sample arrays.</param>
        /// <returns>The jackknife standard error.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="statistic"/> is null.</exception>
        /// <remarks>A single-element sample returns zero without invoking <paramref name="statistic"/>. Parallel callbacks receive independent arrays.</remarks>
        public static double JackKnifeStandardError(IList<double> data, Func<IList<double>, double> statistic)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (statistic == null) throw new ArgumentNullException(nameof(statistic));
            if (data.Count == 0) return double.NaN;

            int N = data.Count;
            if (N == 1) return 0d;
            double theta = statistic(data.ToArray());

            int chunks = Math.Min(JackKnifeChunks, N);
            var chunkSums = new double[chunks];
            Parallel.For(0, chunks, c =>
            {
                double sum = 0d;
                int start = (int)((long)c * N / chunks);
                int end = (int)((long)(c + 1) * N / chunks);
                for (int i = start; i < end; i++)
                {
                    var jackSample = new double[N - 1];
                    for (int k = 0; k < i; k++) jackSample[k] = data[k];
                    for (int k = i + 1; k < N; k++) jackSample[k - 1] = data[k];
                    sum += Tools.Sqr(statistic(jackSample) - theta);
                }
                chunkSums[c] = sum;
            });

            double sumOfSquares = 0d;
            for (int c = 0; c < chunks; c++) sumOfSquares += chunkSums[c];
            return Math.Sqrt((N - 1) / (double)N * sumOfSquares);
        }
        /// <summary>
        /// The number of accumulation chunks used by the jackknife reductions — fixed, so the
        /// floating-point association order does not vary with the machine or the thread count.
        /// </summary>
        private const int JackKnifeChunks = 64;

        /// <summary>
        /// Evaluates a statistic for every leave-one-out sample.
        /// </summary>
        /// <param name="data">Sample data.</param>
        /// <param name="statistic">The statistic evaluated on isolated leave-one-out arrays.</param>
        /// <returns>The statistic values, or <see langword="null"/> for an empty input sample.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="statistic"/> is null.</exception>
        /// <remarks>Each callback receives an independent leave-one-out array; a single-element sample therefore supplies one empty array.</remarks>
        public static double[]? JackKnifeSample(IList<double> data, Func<IList<double>, double> statistic)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (statistic == null) throw new ArgumentNullException(nameof(statistic));
            if (data.Count == 0) return null;

            int N = data.Count;
            var thetaJack = new double[N];
            int chunks = Math.Min(JackKnifeChunks, N);
            Parallel.For(0, chunks, c =>
            {
                int start = (int)((long)c * N / chunks);
                int end = (int)((long)(c + 1) * N / chunks);
                for (int i = start; i < end; i++)
                {
                    var jackSample = new double[N - 1];
                    for (int k = 0; k < i; k++) jackSample[k] = data[k];
                    for (int k = i + 1; k < N; k++) jackSample[k - 1] = data[k];
                    thetaJack[i] = statistic(jackSample);
                }
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
        /// Returns the first four product moments of a sample {Mean, Standard Deviation, bias-corrected Skew, and bias-corrected Excess Kurtosis},
        /// or four NaN values when the sample has fewer than four entries or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        public static double[] ProductMoments(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            double N = data.Count;
            if (N < 4) return [double.NaN, double.NaN, double.NaN, double.NaN];

            // Sums of powers accumulated about the first value. The shift keeps the higher moments
            // conditioned on the sample spread rather than on the distance from zero, so the
            // central moments below are not dominated by cancellation when the mean is large
            // relative to the spread.
            double shift = data[0];
            double X1 = 0, X2 = 0, X3 = 0, X4 = 0;
            foreach (var x in data)
            {
                double y = x - shift;
                double y2 = y * y;
                X1 += y;
                X2 += y2;
                X3 += y2 * y;
                X4 += y2 * y2;
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

            return [shift + U1, S, G, K];

        }

        /// <summary>
        /// Returns the linear moments of a sample {L-Mean (λ1), L-Scale (λ2), L-Skewness (τ3), and L-Kurtosis (τ4)}, or returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <returns>The linear moments {λ1, λ2, τ3, τ4}, or four NaN values when the sample has fewer than four entries.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// The probability weighted moment numerators are accumulated in double precision. Evaluating
        /// them in integer arithmetic overflows silently under the unchecked default: the b₃ numerator
        /// (i-3)(i-2)(i-1) exceeds <see cref="int.MaxValue"/> at i = 1,293 and the b₂ numerator
        /// (i-2)(i-1) exceeds it at i = 46,343, which corrupts τ₃ and τ₄ for large samples. Below those
        /// thresholds the products are exact integers well under 2⁵³, so the double accumulation
        /// reproduces the integer results exactly. See USACE-RMC/Numerics#146.
        /// </remarks>
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
                // Form the b2 and b3 numerators in double so that large samples do not overflow the int products.
                double di = i;
                B0 += sortedData[i - 1];
                if (i > 1)
                    B1 += (i - 1) / (N - 1) * sortedData[i - 1];
                if (i > 2)
                    B2 += (di - 2) * (di - 1) / ((N - 2) * (N - 1)) * sortedData[i - 1];
                if (i > 3)
                    B3 += (di - 3) * (di - 2) * (di - 1) / ((N - 3) * (N - 2) * (N - 1)) * sortedData[i - 1];
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
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="data"/> is empty.</exception>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when <paramref name="k"/> is not a finite value in [0,1], with <c>ParamName</c> equal
        /// to <c>k</c>. Every comparison against NaN is false, so NaN is rejected explicitly; without
        /// that test it reaches the interpolation index, where the float-to-int conversion produces an
        /// out-of-range index and the failure surfaces as an indexer exception rather than this contract
        /// exception.
        /// </exception>
        public static double Percentile(IList<double> data, double k, bool dataIsSorted = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            int n = data.Count;
            if (n == 0) throw new ArgumentException("Sequence contains no elements.", nameof(data));
            if (double.IsNaN(k) || k < 0.0 || k > 1.0) throw new ArgumentOutOfRangeException(nameof(k), "k must be in [0,1].");

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
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="k"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="k"/> has entries and <paramref name="data"/> is empty. An empty <paramref name="k"/> returns an empty array without inspecting <paramref name="data"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any entry of <paramref name="k"/> is not a finite value in [0,1].</exception>
        public static double[] Percentile(IList<double> data, IList<double> k, bool dataIsSorted = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (k == null) throw new ArgumentNullException(nameof(k));

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
        /// <param name="ties">Output. A sparse array of tie-run lengths: the entry at a tie run's last
        /// position in the sorted order holds the run length minus one, and every other entry is zero.
        /// A group of k equal values therefore reports k - 1, so the number of tied groups is the count
        /// of entries greater than zero, not greater than one.</param>
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
            // The loop records a run's length only when the run closes at a later, distinct value, so a
            // tie run containing the largest values never closes inside the loop and its length would be
            // silently dropped without this trailing write.
            if (t > 0)
            {
                ties[work.Length - 1] = t;
            }
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

        #region Weighted Statistics

        /// <summary>
        /// Validates a weight vector against its data vector and returns the weight total.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry.</param>
        /// <returns>The sum of the weights.</returns>
        /// <exception cref="ArgumentNullException">Thrown when either vector is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vectors differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite.</exception>
        private static double ValidateWeights(IList<double> data, IList<double> weights)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            if (data.Count != weights.Count) throw new ArgumentException("The data and weights must have the same length.", nameof(weights));
            double total = 0;
            for (int i = 0; i < weights.Count; i++)
            {
                if (!Tools.IsFinite(weights[i]) || weights[i] < 0d)
                    throw new ArgumentOutOfRangeException(nameof(weights), "Weights must be finite and non-negative.");
                total += weights[i];
            }
            return total;
        }

        /// <summary>
        /// Computes the weighted mean and the weighted second, third, and fourth central power sums
        /// with a deterministic, sequential two-pass sweep.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry.</param>
        /// <param name="total">The sum of the weights.</param>
        /// <param name="mean">Output. The weighted mean.</param>
        /// <param name="s2">Output. The weighted sum of squared deviations.</param>
        /// <param name="s3">Output. The weighted sum of cubed deviations.</param>
        /// <param name="s4">Output. The weighted sum of fourth-power deviations.</param>
        private static void WeightedCentralSums(IList<double> data, IList<double> weights, double total,
            out double mean, out double s2, out double s3, out double s4)
        {
            double sum = 0;
            for (int i = 0; i < data.Count; i++)
                sum += weights[i] * data[i];
            mean = sum / total;
            s2 = 0;
            s3 = 0;
            s4 = 0;
            for (int i = 0; i < data.Count; i++)
            {
                double xm = data[i] - mean;
                double w = weights[i];
                double xm2 = xm * xm;
                s2 += w * xm2;
                s3 += w * xm2 * xm;
                s4 += w * xm2 * xm2;
            }
        }

        /// <summary>
        /// Returns the effective sample size implied by the weights under the given interpretation:
        /// the weight total for frequency weights, or W^2 / sum(w^2) for reliability weights.
        /// </summary>
        /// <param name="weights">The non-negative, finite weight of each data entry.</param>
        /// <param name="total">The sum of the weights.</param>
        /// <param name="weightType">The weight interpretation.</param>
        /// <returns>The effective sample size.</returns>
        private static double EffectiveSampleSize(IList<double> weights, double total, WeightType weightType)
        {
            if (weightType == WeightType.Frequency) return total;
            double sumSq = 0;
            for (int i = 0; i < weights.Count; i++)
                sumSq += weights[i] * weights[i];
            return total * total / sumSq;
        }

        /// <summary>
        /// Estimates the weighted arithmetic mean from the unsorted data array.
        /// Returns NaN if data is empty or any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry. Equal weights reproduce the unweighted mean.</param>
        /// <exception cref="ArgumentNullException">Thrown when either vector is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vectors differ in length, or when every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite.</exception>
        public static double Mean(IList<double> data, IList<double> weights)
        {
            double total = ValidateWeights(data, weights);
            if (data.Count == 0) return double.NaN;
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
            double sum = 0;
            for (int i = 0; i < data.Count; i++)
                sum += weights[i] * data[i];
            return sum / total;
        }

        /// <summary>
        /// Estimates the unbiased weighted variance from the unsorted data array.
        /// Returns NaN if data has less than two entries, if any entry is NaN, or if the
        /// bias-correction denominator implied by the weights is not positive.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry.</param>
        /// <param name="weightType">Optional. The weight interpretation. Default = Frequency.</param>
        /// <remarks>
        /// Frequency weights use the denominator W - 1, so integer weights reproduce the unweighted
        /// variance of the replicated sample exactly; reliability weights use W - sum(w^2)/W, which
        /// reduces to N - 1 at equal weights.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either vector is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vectors differ in length, or when every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite.</exception>
        public static double Variance(IList<double> data, IList<double> weights, WeightType weightType = WeightType.Frequency)
        {
            double total = ValidateWeights(data, weights);
            if (data.Count <= 1) return double.NaN;
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
            WeightedCentralSums(data, weights, total, out _, out double s2, out _, out _);
            double denominator;
            if (weightType == WeightType.Frequency)
            {
                denominator = total - 1d;
            }
            else
            {
                double sumSq = 0;
                for (int i = 0; i < weights.Count; i++)
                    sumSq += weights[i] * weights[i];
                denominator = total - sumSq / total;
            }
            if (denominator <= 0d) return double.NaN;
            return s2 / denominator;
        }

        /// <summary>
        /// Estimates the unbiased weighted standard deviation from the unsorted data array.
        /// Returns NaN if data has less than two entries, if any entry is NaN, or if the
        /// bias-correction denominator implied by the weights is not positive.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry.</param>
        /// <param name="weightType">Optional. The weight interpretation. Default = Frequency.</param>
        /// <exception cref="ArgumentNullException">Thrown when either vector is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vectors differ in length, or when every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite.</exception>
        public static double StandardDeviation(IList<double> data, IList<double> weights, WeightType weightType = WeightType.Frequency)
        {
            return Math.Sqrt(Variance(data, weights, weightType));
        }

        /// <summary>
        /// Estimates the weighted skewness coefficient from the unsorted data array.
        /// Returns NaN if data is empty, if any entry is NaN, or if the effective sample size
        /// implied by the weights does not exceed two, where the bias correction has a pole.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry.</param>
        /// <param name="weightType">Optional. The weight interpretation. Default = Frequency.</param>
        /// <remarks>
        /// The adjusted Fisher-Pearson correction of the unweighted estimator is applied with the
        /// effective sample size in place of the count: the weight total for frequency weights
        /// (integer weights reproduce the replicated sample exactly), or W^2 / sum(w^2) for
        /// reliability weights.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either vector is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vectors differ in length, or when every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite.</exception>
        public static double Skewness(IList<double> data, IList<double> weights, WeightType weightType = WeightType.Frequency)
        {
            double total = ValidateWeights(data, weights);
            if (data.Count == 0) return double.NaN;
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
            WeightedCentralSums(data, weights, total, out _, out double s2, out double s3, out _);
            double n = EffectiveSampleSize(weights, total, weightType);
            // The effective sample size is continuous, so the correction's pole at n = 2 is reachable
            // from ordinary weights, e.g. normalized frequency weights; past it the correction flips
            // sign. The unweighted estimator needs n > 2 for the same reason.
            if (n <= 2d) return double.NaN;
            double m2 = s2 / total;
            double m3 = s3 / total;
            double g = m3 / Math.Pow(m2, 3.0d / 2.0d);
            double a = Math.Sqrt(n * (n - 1.0));
            double b = n - 2;
            return a / b * g;
        }

        /// <summary>
        /// Estimates the weighted excess kurtosis from the unsorted data array.
        /// Returns NaN if data is empty, if any entry is NaN, or if the effective sample size
        /// implied by the weights does not exceed three, where the small-sample correction has a pole.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry.</param>
        /// <param name="weightType">Optional. The weight interpretation. Default = Frequency.</param>
        /// <remarks>
        /// The unweighted small-sample correction is applied with the effective sample size in place
        /// of the count: the weight total for frequency weights, or W^2 / sum(w^2) for reliability
        /// weights. The moment ratio is computed from weight-normalized central moments, which makes
        /// reliability weights scale-invariant; integer frequency weights reproduce the replicated
        /// sample, and equal weights the unweighted estimator, to floating-point rounding.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either vector is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vectors differ in length, or when every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite.</exception>
        public static double Kurtosis(IList<double> data, IList<double> weights, WeightType weightType = WeightType.Frequency)
        {
            double total = ValidateWeights(data, weights);
            if (data.Count == 0) return double.NaN;
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
            WeightedCentralSums(data, weights, total, out _, out double s2, out _, out double s4);
            double n = EffectiveSampleSize(weights, total, weightType);
            // The effective sample size is continuous, so the correction's poles at n = 2 and n = 3 are
            // reachable from ordinary weights, e.g. reliability weights concentrated on two entries;
            // past them the correction changes sign. The unweighted estimator needs n > 3 for the same
            // reason.
            if (n <= 3d) return double.NaN;
            double m2 = s2 / total;
            double m4 = s4 / total;
            double a = n * (n + 1) / ((n - 1) * (n - 2) * (n - 3));
            double b = m4 / (m2 * m2) * ((n - 1) * (n - 1) / n);
            double c = (n - 1) * (n - 1) / ((n - 2) * (n - 3));
            return a * b - 3 * c;
        }

        /// <summary>
        /// Returns the weighted k-th percentile of values in a sample.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="k">The k-th percentile to find.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry, aligned with the data order.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The weighted k-th percentile.</returns>
        /// <remarks>
        /// <para>
        /// With the positive-weight points sorted ascending, A(i) the weight strictly below point i
        /// and B(i) the weight strictly above it, each point sits at the plotting position
        /// p(i) = A(i) / (A(i) + B(i)) and the percentile interpolates linearly between adjacent
        /// positions. The construction is reflection-symmetric and strictly monotone, the first and
        /// last points sit at 0 and 1, and equal weights reproduce the unweighted Type 7 positions
        /// i / (n - 1) - with unit weights the interpolation is arithmetically identical to
        /// <see cref="Percentile(IList{double}, double, bool)"/>.
        /// </para>
        /// <para>
        /// Zero-weight entries carry no mass and are excluded. Note that no percentile estimator can
        /// reproduce the unweighted result at equal weights and, simultaneously, the replicated
        /// sample at integer weights - replication changes the Type 7 interior positions - so integer
        /// weights agree with the replicated sample exactly at the plotting positions and to within
        /// the replicated interpolation gaps between them.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="data"/> is empty, the vectors differ in length, or every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not a finite value in [0,1], or a weight is negative or not finite.</exception>
        public static double Percentile(IList<double> data, double k, IList<double> weights, bool dataIsSorted = false)
        {
            double total = ValidateWeights(data, weights);
            if (data.Count == 0) throw new ArgumentException("Sequence contains no elements.", nameof(data));
            if (double.IsNaN(k) || k < 0.0 || k > 1.0) throw new ArgumentOutOfRangeException(nameof(k), "k must be in [0,1].");
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
            PrepareWeightedSample(data, weights, dataIsSorted, out var x, out var w, out var prefix, out double sum);
            return WeightedPercentile(x, w, prefix, sum, k);
        }

        /// <summary>
        /// Returns an array of weighted percentile values of a sample.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="k">The list of k-th percentiles to find.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry, aligned with the data order.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The weighted k-th percentiles.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/>, <paramref name="k"/>, or <paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="k"/> has entries and <paramref name="data"/> is empty, the vectors differ in length, or every weight is zero. An empty <paramref name="k"/> returns an empty array.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any entry of <paramref name="k"/> is not a finite value in [0,1], or a weight is negative or not finite.</exception>
        public static double[] Percentile(IList<double> data, IList<double> k, IList<double> weights, bool dataIsSorted = false)
        {
            double total = ValidateWeights(data, weights);
            if (k == null) throw new ArgumentNullException(nameof(k));
            if (k.Count == 0) return Array.Empty<double>();
            if (data.Count == 0) throw new ArgumentException("Sequence contains no elements.", nameof(data));
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
            PrepareWeightedSample(data, weights, dataIsSorted, out var x, out var w, out var prefix, out double sum);
            var result = new double[k.Count];
            for (int i = 0; i < k.Count; i++)
            {
                if (double.IsNaN(k[i]) || k[i] < 0.0 || k[i] > 1.0) throw new ArgumentOutOfRangeException(nameof(k), "k must be in [0,1].");
                result[i] = WeightedPercentile(x, w, prefix, sum, k[i]);
            }
            return result;
        }

        /// <summary>
        /// Copies the positive-weight sample points, sorts them by value when required, and builds
        /// the exclusive prefix-weight vector.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry.</param>
        /// <param name="dataIsSorted">True when the data (and aligned weights) are already sorted ascending.</param>
        /// <param name="x">Output. The positive-weight data values, sorted ascending.</param>
        /// <param name="w">Output. The aligned positive weights.</param>
        /// <param name="prefix">Output. The weight strictly below each point.</param>
        /// <param name="total">Output. The positive-weight total.</param>
        private static void PrepareWeightedSample(IList<double> data, IList<double> weights, bool dataIsSorted,
            out double[] x, out double[] w, out double[] prefix, out double total)
        {
            int m = 0;
            for (int i = 0; i < weights.Count; i++)
                if (weights[i] > 0d) m++;
            x = new double[m];
            w = new double[m];
            int j = 0;
            for (int i = 0; i < data.Count; i++)
            {
                if (weights[i] > 0d)
                {
                    x[j] = data[i];
                    w[j] = weights[i];
                    j++;
                }
            }
            // Ties among equal data values leave the quantile function unchanged, so an unstable
            // pair sort is safe.
            if (!dataIsSorted) Array.Sort(x, w);
            prefix = new double[m];
            double run = 0;
            for (int i = 0; i < m; i++)
            {
                prefix[i] = run;
                run += w[i];
            }
            total = run;
        }

        /// <summary>
        /// Evaluates the weighted percentile from a prepared sorted sample.
        /// </summary>
        /// <param name="x">The positive-weight data values, sorted ascending.</param>
        /// <param name="w">The aligned positive weights.</param>
        /// <param name="prefix">The weight strictly below each point.</param>
        /// <param name="total">The positive-weight total.</param>
        /// <param name="k">The k-th percentile to find.</param>
        /// <returns>The weighted k-th percentile.</returns>
        private static double WeightedPercentile(double[] x, double[] w, double[] prefix, double total, double k)
        {
            int m = x.Length;
            if (m == 1 || k == 0.0) return x[0];
            if (k == 1.0) return x[m - 1];

            // Plotting positions p(i) = A(i) / (total - w(i)) are strictly increasing with
            // p(0) = 0 and p(m-1) = 1, so a bracket always exists. Binary search for it.
            // In exact arithmetic the denominator is the sum of the other positive weights, so it is
            // positive for m >= 2; in floating point a weight that dominates the total cancels it to
            // zero, which would poison the search with 0/0. A zero denominator means the point holds
            // essentially all the mass, so its position collapses to the boundary it sits against.
            double Position(int i)
            {
                double d = total - w[i];
                return d > 0d ? prefix[i] / d : (i == 0 ? 0d : 1d);
            }

            int lo = 0, hi = m - 1;
            while (hi - lo > 1)
            {
                int mid = (lo + hi) >> 1;
                if (Position(mid) <= k) lo = mid;
                else hi = mid;
            }
            double pLo = Position(lo);
            double pHi = Position(hi);
            if (k <= pLo) return x[lo];
            if (k >= pHi) return x[hi];
            double theta = (k - pLo) / (pHi - pLo);
            return x[lo] + theta * (x[hi] - x[lo]);
        }

        #endregion

    }
}