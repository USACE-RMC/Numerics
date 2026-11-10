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

using Numerics.MachineLearning;
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

            return min;
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

            return max;
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
            if (TryGetArray(data, out var array, out int offset))
            {
                int n = data.Count;
                for (int i = offset; i < offset + n; i++)
                    sum += array[i];
                return sum / n;
            }
            for (int i = 0; i < data.Count; i++)
                sum += data[i];
            return sum / data.Count;
        }

        /// <summary>
        /// Returns the array behind a list that is a <c>double[]</c> or an <see cref="ArraySegment{T}"/> of one, so that
        /// a sequential loop can index the array directly instead of calling the list indexer.
        /// </summary>
        /// <param name="data">The list.</param>
        /// <param name="array">On return, the array behind <paramref name="data"/>, or null.</param>
        /// <param name="offset">On return, the index in <paramref name="array"/> of the first element of
        /// <paramref name="data"/>.</param>
        /// <returns>True when <paramref name="data"/> is backed by an array.</returns>
        /// <remarks>
        /// The direct loop visits the same elements in the same order as the indexer loop, so the result is identical;
        /// it avoids an interface call per element.
        /// </remarks>
        internal static bool TryGetArray(IList<double> data, out double[] array, out int offset)
        {
            if (data is double[] values)
            {
                array = values;
                offset = 0;
                return true;
            }
            if (data is ArraySegment<double> segment && segment.Array != null)
            {
                array = segment.Array;
                offset = segment.Offset;
                return true;
            }
            array = null!;
            offset = 0;
            return false;
        }

        /// <summary>
        /// The number of accumulation chunks used by the parallel mean reduction. Fixed, not derived
        /// from the processor count, so the summation order, and therefore the mean's last bits,
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
        /// serially in chunk order, the same deterministic reduction the bootstrap's jackknife
        /// accumulation uses. The summation tree therefore depends only on the sample length, so the
        /// result is bit-identical on every machine, core count, and scheduler. A processor-count
        /// partition would produce different last bits on different machines, and a shared
        /// accumulator such as <c>Tools.ParallelAdd</c> would be race-free but commit its additions
        /// in thread-scheduler order (the same non-reproducibility), so both are excluded from this
        /// reduction. Samples below the threshold fall through to the sequential
        /// <see cref="Mean(IList{double})"/>.
        /// For samples at or above the threshold, indexed elements are read concurrently. The
        /// caller must not mutate <paramref name="data"/> during the call, and custom
        /// <see cref="IList{T}"/> implementations must support concurrent indexed reads. This
        /// method owns its internal <see cref="Parallel.For(int, int, Action{int})"/> scheduling;
        /// parallel options from an outer operation are not inherited.
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
        /// <remarks>
        /// The sum of squared deviations comes from <see cref="CentralSumOfSquares(IList{double})"/>, the corrected
        /// two-pass algorithm, which keeps its accuracy when the mean is large relative to the spread.
        /// </remarks>
        public static double Variance(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count <= 1) return double.NaN;
            return CentralSumOfSquares(data) / (data.Count - 1);
        }

        /// <summary>
        /// Evaluates the population variance from the full population provided as unsorted array.
        /// On a dataset of size N will use an N normalizer and would thus be biased if applied to a subset.
        /// Returns NaN if data is empty or if any entry is NaN.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <remarks>
        /// The sum of squared deviations comes from <see cref="CentralSumOfSquares(IList{double})"/>, the corrected
        /// two-pass algorithm, which keeps its accuracy when the mean is large relative to the spread.
        /// </remarks>
        public static double PopulationVariance(IList<double> data)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Count == 0) return double.NaN;
            return CentralSumOfSquares(data) / data.Count;
        }

        /// <summary>
        /// Returns the sum of squared deviations from the mean by the corrected two-pass algorithm.
        /// </summary>
        /// <param name="data">Sample of data, not empty.</param>
        /// <returns>Σ(x - m)² - (Σ(x - m))² / n with m = <see cref="Mean(IList{double})"/>, or NaN if any entry is NaN.</returns>
        /// <remarks>
        /// The second term removes, to first order, the rounding error of the computed mean, so the result does not
        /// cancel when the mean is large relative to the spread, as a running or one-pass sum of squares does (Chan,
        /// Golub and LeVeque, 1983, "Algorithms for computing the sample variance: analysis and recommendations",
        /// The American Statistician 37(3), 242-247).
        /// </remarks>
        private static double CentralSumOfSquares(IList<double> data)
        {
            return CentralCrossProductSum(data, data);
        }

        /// <summary>
        /// Returns the sum of cross products of deviations from the means by the corrected two-pass algorithm.
        /// </summary>
        /// <param name="data1">First sample of data, not empty.</param>
        /// <param name="data2">Second sample of data, of the same length.</param>
        /// <returns>Σ(x - mx)(y - my) - Σ(x - mx) Σ(y - my) / n, or NaN if any entry is NaN.</returns>
        /// <remarks>
        /// The second term removes, to first order, the rounding errors of the computed means
        /// (see <see cref="CentralSumOfSquares(IList{double})"/>).
        /// </remarks>
        private static double CentralCrossProductSum(IList<double> data1, IList<double> data2)
        {
            int n = data1.Count;
            double mean1 = Mean(data1);
            double mean2 = ReferenceEquals(data1, data2) ? mean1 : Mean(data2);
            double sum = 0d, sum1 = 0d, sum2 = 0d;
            if (ReferenceEquals(data1, data2) && TryGetArray(data1, out var array, out int offset))
            {
                for (int i = offset; i < offset + n; i++)
                {
                    double d1 = array[i] - mean1;
                    double d2 = array[i] - mean2;
                    sum += d1 * d2;
                    sum1 += d1;
                    sum2 += d2;
                }
                return sum - sum1 * sum2 / n;
            }
            for (int i = 0; i < n; i++)
            {
                double d1 = data1[i] - mean1;
                double d2 = data2[i] - mean2;
                sum += d1 * d2;
                sum1 += d1;
                sum2 += d2;
            }
            return sum - sum1 * sum2 / n;
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
        /// The number of accumulation chunks used by the jackknife reductions. Fixed, so the
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
            // The products are formed in double: in int, (n-1)(n-2)(n-3) overflows silently at n = 1,293.
            // Below that every product is an exact integer, so the values are unchanged.
            double a = n * (n + 1d) / ((n - 1d) * (n - 2d) * (n - 3d));
            double b = s4 / (v * v);
            double c = (n - 1d) * (n - 1d) / ((n - 2d) * (n - 3d));
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

            return CentralCrossProductSum(data1, data2) / (data1.Count - 1);
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

            return CentralCrossProductSum(data1, data2) / data1.Count;
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
        /// Returns the linear moments of a sample {L-Mean (λ1), L-Scale (λ2), L-Skewness (τ3), and L-Kurtosis (τ4)}, or four NaN values when the sample has fewer than four entries or any entry is NaN.
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
            if (data.Count < 4) return [double.NaN, double.NaN, double.NaN, double.NaN];
            return SampleLinearMoments(data, 4);
        }

        /// <summary>
        /// Returns the first <paramref name="nmom"/> sample L-moments {λ1, λ2, τ3, τ4, τ5, ...}: the L-location, the
        /// L-scale and the L-moment ratios τr = λr / λ2 for r ≥ 3.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <param name="nmom">The number of L-moments to return, at least 1; 5 gives τ5 for regional frequency analysis.</param>
        /// <returns>The <paramref name="nmom"/> values {λ1, λ2, τ3, ..., τ_nmom}; all NaN when the sample has fewer than
        /// <paramref name="nmom"/> entries, and NaN when any entry is NaN.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="nmom"/> is less than 1.</exception>
        /// <remarks>
        /// <para>
        /// The unbiased estimator of Hosking (1990), the estimator of the SAMLMR routine of Hosking's LMOMENTS Fortran
        /// package with its unbiased option, evaluated through the probability weighted moments. With the data sorted
        /// ascending, the unbiased probability weighted moments are
        /// b_r = n⁻¹ Σ_(i > r) x(i) (i − 1)(i − 2)...(i − r) / ((n − 1)(n − 2)...(n − r)), and the L-moments are
        /// λ_(r+1) = Σ_k p*_(r,k) b_k with the shifted Legendre coefficients
        /// p*_(r,k) = (−1)^(r−k) (r + k)! / ((k!)² (r − k)!), formed by the recursion
        /// p*_(r,k) = −p*_(r,k−1) (r + k)(r − k + 1) / k² of SAMLMR.
        /// </para>
        /// <para>
        /// The first four values are bitwise those of <see cref="LinearMoments(IList{double})"/>, which shares this
        /// evaluation: λ2 = 2b1 − b0, τ3 = 2(3b2 − b0)/(2b1 − b0) − 3 and τ4 = 5(2(2b3 − 3b2) + b0)/(2b1 − b0) + 6. Each λr
        /// for r ≥ 2 is a difference of probability weighted moments, so it loses relative accuracy when the mean is
        /// large relative to λr; the loss grows with the order because the coefficients p*_(r,k) grow. Against exact
        /// rational evaluation, τ5 is within about 1E-13 of the exact value for lognormal samples, and within about
        /// 2E-11 for samples uniform on [1000, 1015] (n = 6 to 200), whose mean is 305 to 1262 times λ2, with a median
        /// of 401. A constant sample has λ2 = 0, so its ratios are not finite.
        /// </para>
        /// <para>
        /// Practical bound on <paramref name="nmom"/>: the error of τr grows by roughly a factor of 10 to 30 per order.
        /// For lognormal samples of 60 and 300 values, the largest absolute error was about 6E-14 for τ5, 5E-12 for τ8,
        /// 1E-10 for τ10, 3E-9 for τ12, 7E-7 for τ15 and 4E-3 for τ20; by τ25 the ratios are meaningless. Keep
        /// <paramref name="nmom"/> at 10 or below when the ratios must be accurate to 1E-10. For very high orders the
        /// coefficient products overflow double, and the ratios become infinite or NaN: the first non-finite value
        /// appeared at nmom = 105 for n = 1,000, 79 for n = 10,000, 63 for n = 100,000 and 53 for n = 1,000,000.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations
        /// of order statistics. Journal of the Royal Statistical Society: Series B, 52(1), 105-124.
        /// Hosking, J. R. M. and Wallis, J. R. (1997). Regional Frequency Analysis: An Approach Based on L-Moments.
        /// Cambridge University Press.
        /// </para>
        /// </remarks>
        public static double[] LinearMoments(IList<double> data, int nmom)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (nmom < 1) throw new ArgumentOutOfRangeException(nameof(nmom), "The number of L-moments must be at least 1.");
            if (data.Count < nmom)
            {
                var missing = new double[nmom];
                for (int j = 0; j < nmom; j++) missing[j] = double.NaN;
                return missing;
            }
            return SampleLinearMoments(data, nmom);
        }

        /// <summary>
        /// Returns the first <paramref name="nmom"/> unbiased sample L-moments {λ1, λ2, τ3, ...} of at least
        /// <paramref name="nmom"/> values.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed, with at least <paramref name="nmom"/> entries.</param>
        /// <param name="nmom">The number of L-moments, at least 1.</param>
        /// <returns>The L-moments, as documented on <see cref="LinearMoments(IList{double}, int)"/>.</returns>
        /// <remarks>
        /// The probability weighted moment numerators are accumulated in double precision. Evaluating
        /// them in integer arithmetic overflows silently under the unchecked default: the b₃ numerator
        /// (i-3)(i-2)(i-1) exceeds <see cref="int.MaxValue"/> at i = 1,293 and the b₂ numerator
        /// (i-2)(i-1) exceeds it at i = 46,343, which corrupts τ₃ and τ₄ for large samples. Below those
        /// thresholds the products are exact integers well under 2⁵³, so the double accumulation
        /// reproduces the integer results exactly. See USACE-RMC/Numerics#146. Each coefficient
        /// (i − r)...(i − 1) / ((n − r)...(n − 1)) is formed with its products in that order, then divided,
        /// then multiplied by x(i).
        /// </remarks>
        private static double[] SampleLinearMoments(IList<double> data, int nmom)
        {
            double N = data.Count;

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
            var b = new double[Math.Max(nmom, 4)];
            b[0] = B0;
            b[1] = B1;
            b[2] = B2;
            b[3] = B3;

            // Higher orders in a second pass, with the coefficient products in the same order.
            for (int r = 4; r < nmom; r++)
            {
                double denominator = N - r;
                for (int j = r - 1; j >= 1; j--) denominator *= N - j;
                double sum = 0;
                for (int i = r + 1; i <= N; i++)
                {
                    double di = i;
                    double product = di - r;
                    for (int j = r - 1; j >= 1; j--) product *= di - j;
                    sum += product / denominator * sortedData[i - 1];
                }
                b[r] = sum;
            }
            for (int r = 0; r < nmom; r++) b[r] /= N;

            // L-Mean (λ1), L-Scale (λ2), L-Skewness (τ3) and L-Kurtosis (τ4)
            var moments = new double[nmom];
            moments[0] = b[0];
            if (nmom > 1) moments[1] = 2 * b[1] - b[0];
            if (nmom > 2) moments[2] = 2 * (3 * b[2] - b[0]) / (2 * b[1] - b[0]) - 3;
            if (nmom > 3) moments[3] = 5 * (2 * (2 * b[3] - 3 * b[2]) + b[0]) / (2 * b[1] - b[0]) + 6;

            // Higher orders: λ_(r+1) = Σ_k p*_(r,k) b_k, divided by λ2.
            for (int r = 4; r < nmom; r++)
            {
                double ak = r + 1;
                double p = r % 2 == 0 ? 1d : -1d;
                double lambda = p * b[0];
                for (int k = 1; k <= r; k++)
                {
                    double ai = k;
                    p = -p * (ak + ai - 1d) * (ak - ai) / (ai * ai);
                    lambda += p * b[k];
                }
                moments[r] = lambda / moments[1];
            }
            return moments;
        }

        /// <summary>
        /// Returns the sample L-coefficient of variation (L-CV), τ = λ2 / λ1.
        /// </summary>
        /// <param name="data">Sample of data, no sorting is assumed.</param>
        /// <returns>The L-CV, or NaN when the sample has fewer than four entries or any entry is NaN.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// The ratio of the L-scale to the L-location of <see cref="LinearMoments(IList{double})"/>, bitwise
        /// <c>LinearMoments(data)[1] / LinearMoments(data)[0]</c>. It is undefined (infinite or NaN) when λ1 is zero.
        /// <para>
        /// <b> References: </b>
        /// Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations
        /// of order statistics. Journal of the Royal Statistical Society: Series B, 52(1), 105-124.
        /// </para>
        /// </remarks>
        public static double LinearCoefficientOfVariation(IList<double> data)
        {
            var moments = LinearMoments(data);
            return moments[1] / moments[0];
        }

        /// <summary>
        /// Returns the median of a sample.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The middle order statistic for an odd count, or the midpoint of the two middle order statistics for
        /// an even count; NaN if any entry is NaN.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="data"/> is empty.</exception>
        /// <remarks>
        /// For an even count the midpoint is computed as <c>0.5 * (x[lo] + x[hi])</c>, which is correctly rounded
        /// unless the sum overflows; when it does, the midpoint is <c>0.5 * x[lo] + 0.5 * x[hi]</c>, so the median of
        /// two finite values is always finite. <see cref="Percentile(IList{double}, double, bool)"/> at k = 0.5
        /// computes the same midpoint as <c>x[lo] + 0.5 * (x[hi] - x[lo])</c> (Hyndman and Fan Type 7), which has two
        /// roundings and can differ from it in the last bit. Any NaN entry makes the median NaN, as in
        /// <c>numpy.median</c>; with <paramref name="dataIsSorted"/> true every entry is still checked for NaN.
        /// </remarks>
        public static double Median(IList<double> data, bool dataIsSorted = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            int n = data.Count;
            if (n == 0) throw new ArgumentException("Sequence contains no elements.", nameof(data));
            if (dataIsSorted)
            {
                for (int i = 0; i < n; i++)
                    if (double.IsNaN(data[i])) return double.NaN;
                return SortedMedian(data, n);
            }
            var sortedData = SortedCopy(data);
            // The sort places every NaN first.
            if (double.IsNaN(sortedData[0])) return double.NaN;
            return SortedMedian(sortedData, n);
        }

        /// <summary>
        /// Sorts the first <paramref name="count"/> values of a buffer in place and returns their median.
        /// </summary>
        /// <param name="values">The buffer, of length at least <paramref name="count"/>; its first values are reordered.</param>
        /// <param name="count">The number of values, at least one.</param>
        /// <returns>The median, computed as in <see cref="Median(IList{double}, bool)"/>; NaN if any of the values is NaN.</returns>
        /// <remarks>Used by callers that reuse a scratch buffer, so no array is allocated.</remarks>
        internal static double MedianInPlace(double[] values, int count)
        {
            Array.Sort(values, 0, count);
            // The sort places every NaN first.
            if (double.IsNaN(values[0])) return double.NaN;
            return SortedMedian(values, count);
        }

        /// <summary>
        /// Returns the median of ascending values without NaN.
        /// </summary>
        /// <param name="sortedData">The ascending values.</param>
        /// <param name="n">The number of values, at least one.</param>
        /// <returns>The middle value, or the midpoint of the two middle values, finite whenever they are.</returns>
        private static double SortedMedian(IList<double> sortedData, int n)
        {
            if (n % 2 == 1) return sortedData[n / 2];
            double lo = sortedData[n / 2 - 1], hi = sortedData[n / 2];
            double mid = 0.5 * (lo + hi);
            // Halving each value first avoids the overflow of the sum; without overflow the halved sum is kept, since
            // halving a subnormal value first can lose its last bit.
            if (double.IsInfinity(mid) && !double.IsInfinity(lo) && !double.IsInfinity(hi)) mid = 0.5 * lo + 0.5 * hi;
            return mid;
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
            var sortedData = dataIsSorted ? data : SortedCopy(data);

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
            var sortedData = dataIsSorted ? data : SortedCopy(data);
            var result = new double[k.Count];
            for (int i = 0; i < k.Count; i++)
            {
                result[i] = Percentile(sortedData, k[i], true);
            }
            return result;
        }

        /// <summary>
        /// Returns the k-th percentile of a sample by linear interpolation between the order statistics placed at
        /// the given plotting positions, for example the Weibull positions i / (n + 1) (Hyndman and Fan type 6).
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="k">The k-th percentile to find, in [0,1].</param>
        /// <param name="plottingPositionType">The plotting positions of the order statistics.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The k-th percentile.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="data"/> is empty.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not a finite value in
        /// [0,1], or <paramref name="plottingPositionType"/> is not a defined value.</exception>
        /// <remarks>
        /// <para>
        /// The i-th order statistic x(i), i = 1..n, sits at the position p_i of
        /// <see cref="PlottingPositions.Function(int, PlottingPositions.PlottingPostionType)"/>: (i − α)/(n + 1 − 2α)
        /// for the alpha formulas and i / n for Kaplan-Meier. Inverting the piecewise-linear interpolation gives the
        /// one-based rank h = (n + 1 − 2α) k + α, or h = n k for Kaplan-Meier; with j = floor(h) and w = h − j the
        /// percentile is x(j) + w (x(j+1) − x(j)), the same interpolation form as
        /// <see cref="Percentile(IList{double}, double, bool)"/>, and it is held at x(1) for h ≤ 1 and at x(n) for
        /// h ≥ n. For Weibull, Blom and Hazen the same weight is evaluated from k and the two bracketing positions, as
        /// w = (k − p_j) / (p_(j+1) − p_j) with p_j = (j − α)/(n + 1 − 2α), which keeps the bits of k that the rounding of h
        /// loses. Median, Cunnane and Gringorten keep the weight h − j: their alphas (0.3175, 0.4 and 0.44) are not exact
        /// in binary, so their positions carry an extra rounding, and the ratio form measured worse for them against exact
        /// rational evaluation. Weibull is Hyndman and Fan type 6, Hazen type 5 and Kaplan-Meier type 4; the zero-based Type 7
        /// of <see cref="Percentile(IList{double}, double, bool)"/> is not a plotting-position formula.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Hyndman, R. J. and Fan, Y. (1996). Sample quantiles in statistical packages. The American
        /// Statistician, 50(4), 361-365.
        /// </para>
        /// </remarks>
        public static double Percentile(IList<double> data, double k, PlottingPositions.PlottingPostionType plottingPositionType, bool dataIsSorted = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            int n = data.Count;
            if (n == 0) throw new ArgumentException("Sequence contains no elements.", nameof(data));
            if (double.IsNaN(k) || k < 0.0 || k > 1.0) throw new ArgumentOutOfRangeException(nameof(k), "k must be in [0,1].");
            double h = PlottingPositionRank(n, k, plottingPositionType, out double alpha);
            var sortedData = dataIsSorted ? data : SortedCopy(data);
            if (plottingPositionType == PlottingPositions.PlottingPostionType.Weibull ||
                plottingPositionType == PlottingPositions.PlottingPostionType.Blom ||
                plottingPositionType == PlottingPositions.PlottingPostionType.Hazen)
                return InterpolatePlottingPosition(sortedData, k, h, alpha);
            return InterpolateRank(sortedData, h);
        }

        /// <summary>
        /// Returns an array of percentiles of a sample by linear interpolation between the order statistics placed
        /// at the given plotting positions.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="k">The list of k-th percentiles to find, each in [0,1].</param>
        /// <param name="plottingPositionType">The plotting positions of the order statistics.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The k-th percentiles.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="k"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="k"/> has entries and <paramref name="data"/> is empty.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any entry of <paramref name="k"/> is not a finite
        /// value in [0,1], or <paramref name="plottingPositionType"/> is not a defined value.</exception>
        /// <remarks>
        /// Each value equals <see cref="Percentile(IList{double}, double, PlottingPositions.PlottingPostionType, bool)"/>;
        /// the data are sorted once.
        /// </remarks>
        public static double[] Percentile(IList<double> data, IList<double> k, PlottingPositions.PlottingPostionType plottingPositionType, bool dataIsSorted = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (k == null) throw new ArgumentNullException(nameof(k));
            var result = new double[k.Count];
            if (k.Count == 0) return result;
            var sortedData = dataIsSorted ? data : SortedCopy(data);
            for (int i = 0; i < k.Count; i++)
                result[i] = Percentile(sortedData, k[i], plottingPositionType, true);
            return result;
        }

        /// <summary>
        /// Returns an ascending copy of the data, ordered exactly as a stable sort orders it.
        /// </summary>
        /// <param name="data">The data.</param>
        /// <returns>The values in ascending order, NaN first.</returns>
        /// <remarks>
        /// <para>
        /// The result is identical, bit for bit, to a stable sort (the <c>OrderBy</c> of <c>System.Linq</c> under
        /// <see cref="Comparer{T}.Default"/>).
        /// <see cref="Array.Sort{T}(T[])"/> orders the values the same way (NaN first, then ascending) but is not
        /// stable, and the only values that compare equal while differing in their bits are the zeros (−0 and +0)
        /// and the NaNs (their payloads). Each of those groups occupies one contiguous run of the sorted array, so
        /// writing the group's members back into its run in their original order reproduces the stable order
        /// exactly. Every other run of equal values is bitwise uniform, so its order cannot be observed.
        /// </para>
        /// </remarks>
        private static double[] SortedCopy(IList<double> data)
        {
            var sorted = new double[data.Count];
            data.CopyTo(sorted, 0);
            Array.Sort(sorted);

            int nanCount = 0, zeroCount = 0;
            for (int i = 0; i < sorted.Length; i++)
            {
                if (double.IsNaN(sorted[i])) nanCount++;
                else if (sorted[i] == 0d) zeroCount++;
            }
            if (nanCount < 2 && zeroCount < 2) return sorted;

            // The zeros sit after the NaNs and every negative value.
            int zeroStart = nanCount;
            while (zeroStart < sorted.Length && sorted[zeroStart] < 0d) zeroStart++;
            int nanIndex = 0, zeroIndex = zeroStart;
            for (int i = 0; i < data.Count; i++)
            {
                double value = data[i];
                if (double.IsNaN(value)) sorted[nanIndex++] = value;
                else if (value == 0d) sorted[zeroIndex++] = value;
            }
            return sorted;
        }

        /// <summary>
        /// Returns the one-based rank h at which the plotting-position interpolation reaches probability k.
        /// </summary>
        /// <param name="n">The sample size.</param>
        /// <param name="k">The probability, in [0,1].</param>
        /// <param name="plottingPositionType">The plotting positions of the order statistics.</param>
        /// <param name="alpha">On return, the alpha coefficient of the formula; 0 for Kaplan-Meier, where it is unused.</param>
        /// <returns>(n + 1 − 2α) k + α for the alpha formulas; n k for Kaplan-Meier.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="plottingPositionType"/> is not a defined value.</exception>
        private static double PlottingPositionRank(int n, double k, PlottingPositions.PlottingPostionType plottingPositionType, out double alpha)
        {
            alpha = 0d;
            if (plottingPositionType == PlottingPositions.PlottingPostionType.KaplanMeier) return n * k;
            if (!Enum.IsDefined(typeof(PlottingPositions.PlottingPostionType), plottingPositionType))
                throw new ArgumentOutOfRangeException(nameof(plottingPositionType), "Unknown plotting position type.");
            alpha = PlottingPositions.Alpha(plottingPositionType);
            return (n + 1 - 2d * alpha) * k + alpha;
        }

        /// <summary>
        /// Interpolates sorted data at a one-based rank, holding the end values outside [1, n].
        /// </summary>
        /// <param name="sortedData">The ascending data.</param>
        /// <param name="h">The one-based rank.</param>
        /// <returns>x(j) + w (x(j+1) − x(j)) with j = floor(h) and w = h − j.</returns>
        private static double InterpolateRank(IList<double> sortedData, double h)
        {
            int n = sortedData.Count;
            if (h <= 1d) return sortedData[0];
            if (h >= n) return sortedData[n - 1];
            int lower = (int)Math.Floor(h);
            double w = h - lower;
            return sortedData[lower - 1] + w * (sortedData[lower] - sortedData[lower - 1]);
        }

        /// <summary>
        /// Interpolates sorted data at probability k between the plotting positions (i − α) / (n + 1 − 2α), holding the
        /// end values outside [1, n].
        /// </summary>
        /// <param name="sortedData">The ascending data.</param>
        /// <param name="k">The probability, in [0,1].</param>
        /// <param name="h">The one-based rank (n + 1 − 2α) k + α.</param>
        /// <param name="alpha">The alpha coefficient of the plotting positions.</param>
        /// <returns>x(r) + w (x(r+1) − x(r)) with r = floor(h) and w = (k − p_r) / (p_(r+1) − p_r), where
        /// p_r = (r − α)/d and d = n + 1 − 2α.</returns>
        /// <remarks>
        /// The weight is formed from k and the two bracketing plotting positions rather than as h − r, the same weight in
        /// exact arithmetic. It is more accurate than h − r, because k − p_r keeps the significant bits of k that the
        /// rounding of h discards: against exact rational evaluation of the weight at the return-period probabilities
        /// 1 − 1/T (T = 2 to 100, n = 5 to 120), at uniform random probabilities (n = 2 to 1000) and at k = j/1000, both
        /// the mean and the maximum error are lower for Weibull (α = 0), Blom (α = 3/8) and Hazen (α = 1/2). For
        /// Weibull d = n + 1 and p_r = r/d, so the expression reduces exactly to (k − r/d) / ((r + 1)/d − r/d). The
        /// Median, Cunnane and Gringorten formulas keep the weight h − r: their alphas (0.3175, 0.4 and 0.44) are not exact
        /// in binary, so the positions p_r carry an extra rounding, and on the same grids the ratio form raised the mean
        /// result error for Median and Cunnane and the maximum for Gringorten. When h rounds across an integer the weight
        /// can fall just outside [0, 1]; it is clamped, so the result stays between x(r) and x(r+1).
        /// </remarks>
        private static double InterpolatePlottingPosition(IList<double> sortedData, double k, double h, double alpha)
        {
            int n = sortedData.Count;
            if (h <= 1d) return sortedData[0];
            if (h >= n) return sortedData[n - 1];
            int lower = (int)Math.Floor(h);
            double d = n + 1 - 2d * alpha;
            double pLower = (lower - alpha) / d;
            double w = Tools.Clamp((k - pLower) / ((lower + 1 - alpha) / d - pLower), 0d, 1d);
            return sortedData[lower - 1] + w * (sortedData[lower] - sortedData[lower - 1]);
        }

        /// <summary>
        /// Returns the percentiles of each column of a matrix over its finite entries.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable, for example one row
        /// per realization of an ensemble and one column per node of a curve.</param>
        /// <param name="k">The list of k-th percentiles to find, each in [0,1].</param>
        /// <returns>
        /// A 2-D array with one row per entry of <paramref name="k"/> and one column per column of
        /// <paramref name="data"/>; a column without a finite entry gives <see cref="double.NaN"/>.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="k"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any entry of <paramref name="k"/> is not a finite value in [0,1].</exception>
        /// <remarks>
        /// See <see cref="ColumnPercentiles(double[,], IList{double}, out int[])"/>.
        /// </remarks>
        public static double[,] ColumnPercentiles(double[,] data, IList<double> k)
        {
            return ColumnPercentiles(data, k, out _);
        }

        /// <summary>
        /// Returns the percentiles of each column of a matrix over its finite entries, and the number of
        /// non-finite entries excluded from each column.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable, for example one row
        /// per realization of an ensemble and one column per node of a curve.</param>
        /// <param name="k">The list of k-th percentiles to find, each in [0,1].</param>
        /// <param name="excluded">On return, the number of NaN or infinite entries excluded from each column.</param>
        /// <returns>
        /// A 2-D array with one row per entry of <paramref name="k"/> and one column per column of
        /// <paramref name="data"/>; a column without a finite entry gives <see cref="double.NaN"/>.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="k"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any entry of <paramref name="k"/> is not a finite value in [0,1].</exception>
        /// <remarks>
        /// <para>
        /// Each column's finite entries are gathered in row order and passed to
        /// <see cref="Percentile(IList{double}, IList{double}, bool)"/>, so every value is bitwise the
        /// zero-based linear interpolation (Type 7) that member returns for that column. NaN and infinite entries
        /// are excluded rather than sorted in: left in, NaN would sort first and move the low percentiles. With
        /// every entry of a column finite, the column's result is the percentile over all its entries.
        /// </para>
        /// <para>
        /// The columns are evaluated in parallel above a work threshold (sequentially below 1,000 entries), each into
        /// its own result column from its own gathered copy,
        /// so the result does not depend on the machine or the thread count. One call replaces one call per
        /// column.
        /// </para>
        /// </remarks>
        public static double[,] ColumnPercentiles(double[,] data, IList<double> k, out int[] excluded)
        {
            return ColumnPercentiles(data, k, column => Percentile(column, k, false), out excluded);
        }

        /// <summary>
        /// Returns the plotting-position percentiles of each column of a matrix over its finite entries.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable, for example one row
        /// per year and one column per grid cell.</param>
        /// <param name="k">The list of k-th percentiles to find, each in [0,1].</param>
        /// <param name="plottingPositionType">The plotting positions of the order statistics.</param>
        /// <returns>
        /// A 2-D array with one row per entry of <paramref name="k"/> and one column per column of
        /// <paramref name="data"/>; a column without a finite entry gives <see cref="double.NaN"/>.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="k"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any entry of <paramref name="k"/> is not a finite
        /// value in [0,1], or <paramref name="plottingPositionType"/> is not a defined value.</exception>
        /// <remarks>
        /// See <see cref="ColumnPercentiles(double[,], IList{double}, PlottingPositions.PlottingPostionType, out int[])"/>.
        /// </remarks>
        public static double[,] ColumnPercentiles(double[,] data, IList<double> k, PlottingPositions.PlottingPostionType plottingPositionType)
        {
            return ColumnPercentiles(data, k, plottingPositionType, out _);
        }

        /// <summary>
        /// Returns the plotting-position percentiles of each column of a matrix over its finite entries, and the
        /// number of non-finite entries excluded from each column.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable, for example one row
        /// per year and one column per grid cell.</param>
        /// <param name="k">The list of k-th percentiles to find, each in [0,1].</param>
        /// <param name="plottingPositionType">The plotting positions of the order statistics.</param>
        /// <param name="excluded">On return, the number of NaN or infinite entries excluded from each column.</param>
        /// <returns>
        /// A 2-D array with one row per entry of <paramref name="k"/> and one column per column of
        /// <paramref name="data"/>; a column without a finite entry gives <see cref="double.NaN"/>.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="k"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any entry of <paramref name="k"/> is not a finite
        /// value in [0,1], or <paramref name="plottingPositionType"/> is not a defined value.</exception>
        /// <remarks>
        /// <para>
        /// The plotting-position counterpart of <see cref="ColumnPercentiles(double[,], IList{double}, out int[])"/>:
        /// each column's finite entries are gathered in row order and passed to
        /// <see cref="Percentile(IList{double}, IList{double}, PlottingPositions.PlottingPostionType, bool)"/>, so every
        /// value is bitwise the percentile that member returns for that column (Weibull is Hyndman and Fan type 6).
        /// NaN and infinite entries are excluded, and the columns are evaluated in parallel above a work threshold, each into its own result
        /// column, so the result does not depend on the machine or the thread count.
        /// </para>
        /// <para>
        /// One call replaces one call per column, which matters when the caller crosses a language boundary for
        /// each call (for example the annual maxima of every cell of a grid).
        /// </para>
        /// </remarks>
        public static double[,] ColumnPercentiles(double[,] data, IList<double> k, PlottingPositions.PlottingPostionType plottingPositionType, out int[] excluded)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (k == null) throw new ArgumentNullException(nameof(k));
            if (!Enum.IsDefined(typeof(PlottingPositions.PlottingPostionType), plottingPositionType))
                throw new ArgumentOutOfRangeException(nameof(plottingPositionType), "Unknown plotting position type.");
            return ColumnPercentiles(data, k, column => Percentile(column, k, plottingPositionType, false), out excluded);
        }

        /// <summary>
        /// Evaluates a percentile member over the finite entries of each column of a matrix.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable.</param>
        /// <param name="k">The list of k-th percentiles to find, each in [0,1].</param>
        /// <param name="columnPercentiles">Returns the percentiles of one column's gathered finite entries, one per
        /// entry of <paramref name="k"/>.</param>
        /// <param name="excluded">On return, the number of NaN or infinite entries excluded from each column.</param>
        /// <returns>One row per entry of <paramref name="k"/> and one column per column of <paramref name="data"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="k"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any entry of <paramref name="k"/> is not a finite value in [0,1].</exception>
        private static double[,] ColumnPercentiles(double[,] data, IList<double> k, Func<IList<double>, double[]> columnPercentiles, out int[] excluded)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (k == null) throw new ArgumentNullException(nameof(k));
            for (int i = 0; i < k.Count; i++)
            {
                if (double.IsNaN(k[i]) || k[i] < 0.0 || k[i] > 1.0) throw new ArgumentOutOfRangeException(nameof(k), "k must be in [0,1].");
            }
            return FiniteColumnStatistics(data, k.Count, columnPercentiles, SortWorkPerEntry, out excluded);
        }

        /// <summary>
        /// Returns the arithmetic mean of each column of a matrix over its finite entries.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable, for example one row
        /// per bootstrap replicate and one column per node of a curve.</param>
        /// <returns>One mean per column of <paramref name="data"/>; a column without a finite entry gives
        /// <see cref="double.NaN"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// See <see cref="ColumnMeans(double[,], out int[])"/>.
        /// </remarks>
        public static double[] ColumnMeans(double[,] data)
        {
            return ColumnMeans(data, out _);
        }

        /// <summary>
        /// Returns the arithmetic mean of each column of a matrix over its finite entries, and the number of
        /// non-finite entries excluded from each column.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable, for example one row
        /// per bootstrap replicate and one column per node of a curve.</param>
        /// <param name="excluded">On return, the number of NaN or infinite entries excluded from each column.</param>
        /// <returns>One mean per column of <paramref name="data"/>; a column without a finite entry gives
        /// <see cref="double.NaN"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// <para>
        /// The mean counterpart of <see cref="ColumnPercentiles(double[,], IList{double}, out int[])"/>, with the same
        /// rule for non-finite entries: each column's finite entries are gathered in row order and passed to
        /// <see cref="Mean(IList{double})"/>, so every value is bitwise the sequential mean of those entries, and
        /// NaN and infinite entries are excluded. With every entry of a column finite, the value is
        /// <see cref="Mean(IList{double})"/> of the whole column.
        /// </para>
        /// <para>
        /// The columns are evaluated in parallel above a work threshold (sequentially below 10,000 entries), each
        /// into its own result entry, so the result does not depend
        /// on the machine or the thread count. Evaluates all columns in one call.
        /// </para>
        /// </remarks>
        public static double[] ColumnMeans(double[,] data, out int[] excluded)
        {
            return FiniteColumnStatistics(data, 1, column => new[] { Mean(column) }, SumWorkPerEntry, out excluded).GetRow(0);
        }

        /// <summary>
        /// Returns the unbiased sample variance of each column of a matrix over its finite entries.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable.</param>
        /// <returns>One variance per column of <paramref name="data"/>; NaN for a column with fewer than two finite
        /// entries.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// See <see cref="ColumnVariances(double[,], out int[])"/>.
        /// </remarks>
        public static double[] ColumnVariances(double[,] data)
        {
            return ColumnVariances(data, out _);
        }

        /// <summary>
        /// Returns the unbiased sample variance of each column of a matrix over its finite entries, and the number
        /// of non-finite entries excluded from each column.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable.</param>
        /// <param name="excluded">On return, the number of NaN or infinite entries excluded from each column.</param>
        /// <returns>One variance per column of <paramref name="data"/>; NaN for a column with fewer than two finite
        /// entries.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// The variance counterpart of <see cref="ColumnMeans(double[,], out int[])"/>, with the same rule for
        /// non-finite entries and the same parallel evaluation: every value is bitwise
        /// <see cref="Variance(IList{double})"/> (n − 1 divisor, corrected two-pass) of the column's finite entries
        /// in row order.
        /// </remarks>
        public static double[] ColumnVariances(double[,] data, out int[] excluded)
        {
            return FiniteColumnStatistics(data, 1, column => new[] { Variance(column) }, SumWorkPerEntry, out excluded).GetRow(0);
        }

        /// <summary>
        /// Returns the unbiased sample standard deviation of each column of a matrix over its finite entries.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable.</param>
        /// <returns>One standard deviation per column of <paramref name="data"/>; NaN for a column with fewer than
        /// two finite entries.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// See <see cref="ColumnStandardDeviations(double[,], out int[])"/>.
        /// </remarks>
        public static double[] ColumnStandardDeviations(double[,] data)
        {
            return ColumnStandardDeviations(data, out _);
        }

        /// <summary>
        /// Returns the unbiased sample standard deviation of each column of a matrix over its finite entries, and the
        /// number of non-finite entries excluded from each column.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable.</param>
        /// <param name="excluded">On return, the number of NaN or infinite entries excluded from each column.</param>
        /// <returns>One standard deviation per column of <paramref name="data"/>; NaN for a column with fewer than
        /// two finite entries.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// Every value is bitwise <see cref="StandardDeviation(IList{double})"/> of the column's finite entries, the
        /// square root of the corresponding <see cref="ColumnVariances(double[,], out int[])"/> value.
        /// </remarks>
        public static double[] ColumnStandardDeviations(double[,] data, out int[] excluded)
        {
            return FiniteColumnStatistics(data, 1, column => new[] { StandardDeviation(column) }, SumWorkPerEntry, out excluded).GetRow(0);
        }

        /// <summary>
        /// Returns the linear moments of each column of a matrix over its finite entries.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable.</param>
        /// <returns>A 4 × p array with rows λ1, λ2, τ3 and τ4 and one column per column of <paramref name="data"/>;
        /// NaN for a column with fewer than four finite entries.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// See <see cref="ColumnLinearMoments(double[,], out int[])"/>.
        /// </remarks>
        public static double[,] ColumnLinearMoments(double[,] data)
        {
            return ColumnLinearMoments(data, out _);
        }

        /// <summary>
        /// Returns the linear moments of each column of a matrix over its finite entries, and the number of
        /// non-finite entries excluded from each column.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable.</param>
        /// <param name="excluded">On return, the number of NaN or infinite entries excluded from each column.</param>
        /// <returns>A 4 × p array with rows λ1, λ2, τ3 and τ4 and one column per column of <paramref name="data"/>;
        /// NaN for a column with fewer than four finite entries.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// Every column is bitwise <see cref="LinearMoments(IList{double})"/> of the column's finite entries, with the
        /// same rule for non-finite entries and the same parallel evaluation as
        /// <see cref="ColumnMeans(double[,], out int[])"/>. The L-CV of a column, row 1 divided by row 0, is bitwise
        /// <see cref="LinearCoefficientOfVariation(IList{double})"/>.
        /// </remarks>
        public static double[,] ColumnLinearMoments(double[,] data, out int[] excluded)
        {
            return FiniteColumnStatistics(data, 4, LinearMoments, SortWorkPerEntry, out excluded);
        }


        /// <summary>
        /// The work per matrix entry passed to <see cref="RangeParallel.For(int, long, Action{int})"/> by the column
        /// means, variances and standard deviations, whose per-column cost is a few passes over the entries.
        /// </summary>
        /// <remarks>
        /// With the work threshold of <see cref="RangeParallel.WorkThreshold"/> the columns run sequentially below
        /// 10,000 entries, where scheduling the parallel loop costs more than it saves (measured on 22 cores: at
        /// 100 × 100 the two paths take the same time, and at 10 × 8 the sequential path is faster).
        /// </remarks>
        private const long SumWorkPerEntry = 20;

        /// <summary>
        /// The work per matrix entry passed to <see cref="RangeParallel.For(int, long, Action{int})"/> by the column
        /// percentiles and L-moments, whose per-column cost is dominated by a sort.
        /// </summary>
        /// <remarks>
        /// The columns run sequentially below 1,000 entries; above that the parallel loop is faster (measured on
        /// 22 cores: 0.04 ms against 0.27 ms for the L-moments of a 100 × 100 matrix).
        /// </remarks>
        private const long SortWorkPerEntry = 200;

        /// <summary>
        /// Evaluates a statistic over the finite entries of each column of a matrix, in parallel over the columns
        /// when the matrix is large enough.
        /// </summary>
        /// <param name="data">The samples, one row per member and one column per variable.</param>
        /// <param name="outputs">The number of values the statistic returns for one column.</param>
        /// <param name="statistic">Returns the <paramref name="outputs"/> values of one column's gathered finite entries.</param>
        /// <param name="workPerEntry">The approximate cost of the statistic per entry, which sets the matrix size below
        /// which the columns run sequentially (see <see cref="RangeParallel"/>).</param>
        /// <param name="excluded">On return, the number of NaN or infinite entries excluded from each column.</param>
        /// <returns>One row per output and one column per column of <paramref name="data"/>; NaN for a column
        /// without a finite entry.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <remarks>
        /// Each column is gathered and evaluated independently into its own result entries, so the result is
        /// bitwise the same whether the columns run sequentially or in parallel, on any machine.
        /// </remarks>
        private static double[,] FiniteColumnStatistics(double[,] data, int outputs, Func<IList<double>, double[]> statistic, long workPerEntry, out int[] excluded)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            int rows = data.GetLength(0);
            int columns = data.GetLength(1);
            var result = new double[outputs, columns];
            var counts = new int[columns];
            RangeParallel.For(columns, (long)rows * columns * workPerEntry, j =>
            {
                var column = new double[rows];
                int finite = 0;
                for (int r = 0; r < rows; r++)
                {
                    double value = data[r, j];
                    if (!double.IsNaN(value) && !double.IsInfinity(value)) column[finite++] = value;
                }
                counts[j] = rows - finite;
                if (finite == 0)
                {
                    for (int i = 0; i < outputs; i++) result[i, j] = double.NaN;
                    return;
                }
                var values = statistic(new ArraySegment<double>(column, 0, finite));
                for (int i = 0; i < outputs; i++) result[i, j] = values[i];
            });
            excluded = counts;
            return result;
        }

        /// <summary>
        /// Returns the Wilson score interval for a binomial proportion.
        /// </summary>
        /// <param name="successes">The number of successes, in [0, <paramref name="trials"/>].</param>
        /// <param name="trials">The number of trials; must be positive.</param>
        /// <param name="confidenceLevel">The two-sided confidence level, in (0, 1); for example 0.9.</param>
        /// <returns>The lower and upper bounds, clamped to [0, 1].</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="trials"/> is not positive and
        /// finite, <paramref name="successes"/> is outside [0, <paramref name="trials"/>], or
        /// <paramref name="confidenceLevel"/> is not strictly between zero and one.</exception>
        /// <remarks>
        /// <para>
        /// With p̂ = s / n and z the standard normal quantile at 0.5 + level / 2, the center is
        /// (p̂ + z²/2n) / (1 + z²/n) and the half-width is z √(p̂(1 − p̂)/n + z²/4n²) / (1 + z²/n), evaluated in that
        /// order with z² / n formed once. The counts may be fractional, for example weighted counts.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. Journal
        /// of the American Statistical Association, 22(158), 209-212.
        /// </para>
        /// </remarks>
        public static (double Lower, double Upper) WilsonScoreInterval(double successes, double trials, double confidenceLevel)
        {
            if (!Tools.IsFinite(trials) || trials <= 0d)
                throw new ArgumentOutOfRangeException(nameof(trials), "The number of trials must be positive and finite.");
            if (double.IsNaN(successes) || successes < 0d || successes > trials)
                throw new ArgumentOutOfRangeException(nameof(successes), "The number of successes must be between zero and the number of trials.");
            if (double.IsNaN(confidenceLevel) || confidenceLevel <= 0d || confidenceLevel >= 1d)
                throw new ArgumentOutOfRangeException(nameof(confidenceLevel), "The confidence level must be between 0 and 1.");

            double z = new Distributions.Normal(0.0, 1.0).InverseCDF(0.5 + confidenceLevel / 2.0);
            double pHat = successes / trials;
            double zSquaredOverN = z * z / trials;
            double denominator = 1.0 + zSquaredOverN;
            double center = (pHat + zSquaredOverN / 2.0) / denominator;
            double halfWidth = z * Math.Sqrt(pHat * (1.0 - pHat) / trials + zSquaredOverN / (4.0 * trials)) / denominator;
            return (Math.Max(0.0, center - halfWidth), Math.Min(1.0, center + halfWidth));
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
        /// <returns>The one-based ranks; exactly equal values, infinities included, share their average rank.</returns>
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

                // Exact equality: a difference test would split equal infinities (inf - inf is NaN).
                if (work[i] == work[previousIndex])
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
        /// <returns>The one-based ranks; exactly equal values, infinities included, share their average rank, the same
        /// ties as <see cref="RanksInPlace(double[])"/>.</returns>
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
                // Exact equality, the tie rule of RanksInPlace(double[]).
                if (work[i] == work[previousIndex])
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
        /// <remarks>
        /// The second power sum is corrected by (Σw(x - m))² / W, which removes the rounding error of the computed mean to
        /// first order, as in the corrected two-pass algorithm of <see cref="CentralSumOfSquares(IList{double})"/>.
        /// </remarks>
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
            double s1 = 0;
            for (int i = 0; i < data.Count; i++)
            {
                double xm = data[i] - mean;
                double w = weights[i];
                double xm2 = xm * xm;
                s1 += w * xm;
                s2 += w * xm2;
                s3 += w * xm2 * xm;
                s4 += w * xm2 * xm2;
            }
            s2 -= s1 * s1 / total;
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
        /// Computes the Kish effective sample size of a weight vector, (sum w)^2 / sum(w^2).
        /// Returns NaN if the weight vector is empty.
        /// </summary>
        /// <param name="weights">The non-negative, finite weights.</param>
        /// <returns>The effective sample size, between one and the number of positive weights.</returns>
        /// <remarks>
        /// <para>
        /// This is the effective sample size the reliability-weighted statistics use in their bias
        /// corrections (see <see cref="WeightType.Reliability"/>), computed by the same code. Equal
        /// positive weights give the count, and a single dominant weight gives a value near one.
        /// </para>
        /// <para>
        /// Summation order: the weight total W is accumulated sequentially from zero in index order,
        /// the square sum S is then accumulated sequentially from zero in index order as
        /// S += w[i] * w[i], and the result is (W * W) / S. No pairwise or compensated summation is
        /// used, so the bits are reproducible from the same left-to-right loop in any language. The
        /// result can differ in the last bits from a pairwise sum such as numpy.sum for more than
        /// eight weights.
        /// </para>
        /// <para>
        /// The weights are not rescaled: when W * W or S underflows to zero (for example ten weights
        /// of 1e-200) the result is NaN, and weights above about 1e154 overflow the squares.
        /// Rescale such weights by their maximum first; the ratio is scale invariant.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Kish, L. (1965). Survey Sampling. Wiley, New York.
        /// Kish, L. (1992). Weighting for unequal Pi. Journal of Official Statistics, 8(2), 183-200.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when the weight vector is null.</exception>
        /// <exception cref="ArgumentException">Thrown when every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite.</exception>
        public static double EffectiveSampleSize(IList<double> weights)
        {
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            double total = ValidateWeights(weights, weights);
            if (weights.Count == 0) return double.NaN;
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
            return EffectiveSampleSize(weights, total, WeightType.Reliability);
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
        /// Returns the weighted k-th percentile of a sample by linear interpolation between the order statistics placed
        /// at weighted plotting positions, for example the weighted Hazen positions (c_(i−1) + w_i / 2) / W.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="k">The k-th percentile to find, in [0,1].</param>
        /// <param name="weights">The non-negative, finite weight of each data entry, aligned with the data order.</param>
        /// <param name="plottingPositionType">The plotting positions of the order statistics.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The weighted k-th percentile.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> or <paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="data"/> is empty or holds NaN, the vectors differ
        /// in length, or every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="k"/> is not a finite value in [0,1], a
        /// weight is negative or not finite, the total weight overflows, or <paramref name="plottingPositionType"/> is
        /// not a defined value.</exception>
        /// <remarks>
        /// <para>
        /// With the points in ascending order, c_(i−1) the sequential sum of the weights before point i, W the total
        /// weight and w̄ = W / n the mean weight over all n points, point i sits at
        /// F_i = (c_(i−1) + (1 − α) w_i) / (W + (1 − 2α) w̄) for the alpha formulas, the position of
        /// <see cref="PlottingPositions.Function(IList{double}, PlottingPositions.PlottingPostionType)"/> in exact
        /// arithmetic, and at F_i = (c_(i−1) + w_i) / W for Kaplan-Meier. For Hazen (α = 1/2) this is
        /// (c_(i−1) + w_i / 2) / W, evaluated with exactly those operations. The percentile interpolates linearly
        /// between successive (F, x) points, as x(i−1) + (x(i) − x(i−1)) (k − F_(i−1)) / (F_i − F_(i−1)); it returns x(i)
        /// where two positions coincide, and it is held at the first value for k at or below F_1 and at the last value
        /// for k above F_n.
        /// </para>
        /// <para>
        /// Zero-weight entries carry no mass but take a position. This differs from
        /// <see cref="Percentile(IList{double}, double, IList{double}, bool)"/>, whose positions follow the weight below
        /// and above each positive-weight point. Unsorted data are sorted by value with ties kept in their original
        /// order. With <paramref name="dataIsSorted"/> true no array is allocated: the weights are summed once for W and
        /// the positions are formed during the interpolation walk.
        /// </para>
        /// </remarks>
        public static double Percentile(IList<double> data, double k, IList<double> weights, PlottingPositions.PlottingPostionType plottingPositionType, bool dataIsSorted = false)
        {
            ValidateWeightedPlottingPositionSample(data, weights, out var x, out var w, out double total, out bool hasNaN);
            if (x.Length == 0) throw new ArgumentException("Sequence contains no elements.", nameof(data));
            if (double.IsNaN(k) || k < 0.0 || k > 1.0) throw new ArgumentOutOfRangeException(nameof(k), "k must be in [0,1].");
            double share = WeightedPositionShare(plottingPositionType, out double alpha);
            CheckWeightedPlottingPositionTotal(total, hasNaN);
            bool kaplanMeier = plottingPositionType == PlottingPositions.PlottingPostionType.KaplanMeier;
            if (!dataIsSorted) total = SortWeightedSample(ref x, ref w);
            double denominator = WeightedPositionDenominator(total, x.Length, alpha, kaplanMeier);
            int index = 0;
            double before = 0d, previousPosition = double.NaN, previousValue = double.NaN;
            return WeightedPlottingPositionPercentile(x, w, denominator, k, share, ref index, ref before, ref previousPosition, ref previousValue);
        }

        /// <summary>
        /// Returns an array of weighted percentiles of a sample by linear interpolation between the order statistics
        /// placed at weighted plotting positions.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="k">The list of k-th percentiles to find, each in [0,1], in any order.</param>
        /// <param name="weights">The non-negative, finite weight of each data entry, aligned with the data order.</param>
        /// <param name="plottingPositionType">The plotting positions of the order statistics.</param>
        /// <param name="dataIsSorted">Boolean value indicating if the data is sorted or not. Assumed false, not sorted, by default.</param>
        /// <returns>The weighted k-th percentiles, in the order of <paramref name="k"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/>, <paramref name="k"/> or
        /// <paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vectors differ in length, or when <paramref name="k"/> has
        /// entries and <paramref name="data"/> is empty or holds NaN, or every weight is zero. An empty
        /// <paramref name="k"/> returns an empty array.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite, or when
        /// <paramref name="k"/> has entries and any of them is not a finite value in [0,1], the total weight overflows,
        /// or <paramref name="plottingPositionType"/> is not a defined value.</exception>
        /// <remarks>
        /// Each value equals, bit for bit,
        /// <see cref="Percentile(IList{double}, double, IList{double}, PlottingPositions.PlottingPostionType, bool)"/>
        /// at that k. The sample is validated once, sorted once and summed once for the total weight W, and the
        /// probabilities are served in ascending order by one walk over the positions, so the cost is
        /// O(n log n + m log m) for n points and m probabilities rather than O(m n log n). Probabilities that are not
        /// ascending are visited in sorted order and their results written back to their own entries.
        /// </remarks>
        public static double[] Percentile(IList<double> data, IList<double> k, IList<double> weights, PlottingPositions.PlottingPostionType plottingPositionType, bool dataIsSorted = false)
        {
            ValidateWeightedPlottingPositionSample(data, weights, out var x, out var w, out double total, out bool hasNaN);
            if (k == null) throw new ArgumentNullException(nameof(k));
            int m = k.Count;
            if (m == 0) return Array.Empty<double>();
            if (x.Length == 0) throw new ArgumentException("Sequence contains no elements.", nameof(data));
            bool ascending = true;
            for (int j = 0; j < m; j++)
            {
                double kj = k[j];
                if (double.IsNaN(kj) || kj < 0.0 || kj > 1.0) throw new ArgumentOutOfRangeException(nameof(k), "k must be in [0,1].");
                if (j > 0 && kj < k[j - 1]) ascending = false;
            }
            double share = WeightedPositionShare(plottingPositionType, out double alpha);
            CheckWeightedPlottingPositionTotal(total, hasNaN);
            bool kaplanMeier = plottingPositionType == PlottingPositions.PlottingPostionType.KaplanMeier;
            if (!dataIsSorted) total = SortWeightedSample(ref x, ref w);
            double denominator = WeightedPositionDenominator(total, x.Length, alpha, kaplanMeier);

            int[]? order = null;
            if (!ascending)
            {
                var keys = k.ToArray();
                order = new int[m];
                for (int j = 0; j < m; j++) order[j] = j;
                // Equal probabilities give equal results, so an unstable sort is safe.
                Array.Sort(keys, order);
            }
            var result = new double[m];
            int index = 0;
            double before = 0d, previousPosition = double.NaN, previousValue = double.NaN;
            for (int j = 0; j < m; j++)
            {
                int target = order == null ? j : order[j];
                result[target] = WeightedPlottingPositionPercentile(x, w, denominator, k[target], share, ref index, ref before, ref previousPosition, ref previousValue);
            }
            return result;
        }

        /// <summary>
        /// Validates a weighted sample for the weighted plotting-position percentile in one pass.
        /// </summary>
        /// <param name="data">Sample of data.</param>
        /// <param name="weights">The weight of each data entry.</param>
        /// <param name="x">On return, the data as an array: <paramref name="data"/> itself when it is one, otherwise a copy.</param>
        /// <param name="w">On return, the weights as an array: <paramref name="weights"/> itself when it is one, otherwise a copy.</param>
        /// <param name="total">On return, the sequential sum of the weights in data order.</param>
        /// <param name="hasNaN">On return, true when a data entry is NaN.</param>
        /// <exception cref="ArgumentNullException">Thrown when either vector is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vectors differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite.</exception>
        private static void ValidateWeightedPlottingPositionSample(IList<double> data, IList<double> weights, out double[] x, out double[] w, out double total, out bool hasNaN)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            int n = data.Count;
            if (n != weights.Count) throw new ArgumentException("The data and weights must have the same length.", nameof(weights));
            // Arrays are read directly; other lists are copied once.
            x = data as double[] ?? data.ToArray();
            w = weights as double[] ?? weights.ToArray();

            // One pass validates the weights, forms their sequential total and notes NaN data.
            double sum = 0d, probe = 0d;
            for (int i = 0; i < n; i++)
            {
                double weight = w[i];
                // False for NaN, a negative weight and +∞.
                if (!(weight >= 0d && weight <= double.MaxValue))
                    throw new ArgumentOutOfRangeException(nameof(weights), "Weights must be finite and non-negative.");
                sum += weight;
                // x − x is NaN exactly for NaN and infinite values, which the scan below tells apart.
                probe += x[i] - x[i];
            }
            total = sum;
            hasNaN = false;
            if (double.IsNaN(probe))
                for (int i = 0; i < n && !hasNaN; i++) hasNaN = double.IsNaN(x[i]);
        }

        /// <summary>
        /// Refuses a weight total that overflowed or is zero, and data that hold NaN.
        /// </summary>
        /// <param name="total">The sequential sum of the weights.</param>
        /// <param name="hasNaN">True when a data entry is NaN.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="total"/> is infinite.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="total"/> is zero or <paramref name="hasNaN"/> is true.</exception>
        private static void CheckWeightedPlottingPositionTotal(double total, bool hasNaN)
        {
            if (double.IsInfinity(total)) throw new ArgumentOutOfRangeException("weights", "The total weight must be finite.");
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", "weights");
            if (hasNaN) throw new ArgumentException("The data must not contain NaN.", "data");
        }

        /// <summary>
        /// Replaces validated data and weights with copies sorted by value, with ties in their original order, and returns
        /// the sequential total of the weights in that order.
        /// </summary>
        /// <param name="x">The validated data, unsorted; on return, a sorted copy.</param>
        /// <param name="w">The validated weights, aligned with <paramref name="x"/>; on return, the copy aligned with the sorted data.</param>
        /// <returns>The sequential sum of the weights in sorted order.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the total weight in sorted order overflows.</exception>
        /// <remarks>Kept apart from the sorted path so that the comparison's closure is allocated only here.</remarks>
        private static double SortWeightedSample(ref double[] x, ref double[] w)
        {
            var values = x;
            int n = values.Length;
            var order = new int[n];
            for (int i = 0; i < n; i++) order[i] = i;
            Array.Sort(order, (a, c) => { int byValue = values[a].CompareTo(values[c]); return byValue != 0 ? byValue : a.CompareTo(c); });
            var xs = new double[n];
            var ws = new double[n];
            double sortedTotal = 0d;
            for (int i = 0; i < n; i++)
            {
                xs[i] = values[order[i]];
                ws[i] = w[order[i]];
                sortedTotal += ws[i];
            }
            // Summed in a different order, the total can overflow where the original order did not.
            if (double.IsInfinity(sortedTotal)) throw new ArgumentOutOfRangeException("weights", "The total weight must be finite.");
            x = xs;
            w = ws;
            return sortedTotal;
        }

        /// <summary>
        /// Returns the share (1 − α) of a point's own weight below its weighted plotting position, and α, for a
        /// plotting position type; Kaplan-Meier has share 1 and α = 0.
        /// </summary>
        /// <param name="plottingPositionType">The plotting position type.</param>
        /// <param name="alpha">On return, the alpha coefficient.</param>
        /// <returns>The share 1 − α.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="plottingPositionType"/> is not a defined value.</exception>
        /// <remarks>A switch over the defined values, so no value is boxed.</remarks>
        private static double WeightedPositionShare(PlottingPositions.PlottingPostionType plottingPositionType, out double alpha)
        {
            switch (plottingPositionType)
            {
                case PlottingPositions.PlottingPostionType.KaplanMeier:
                    alpha = 0d;
                    return 1d;
                case PlottingPositions.PlottingPostionType.Weibull:
                case PlottingPositions.PlottingPostionType.Median:
                case PlottingPositions.PlottingPostionType.Blom:
                case PlottingPositions.PlottingPostionType.Cunnane:
                case PlottingPositions.PlottingPostionType.Gringorten:
                case PlottingPositions.PlottingPostionType.Hazen:
                    alpha = PlottingPositions.Alpha(plottingPositionType);
                    return 1d - alpha;
                default:
                    throw new ArgumentOutOfRangeException(nameof(plottingPositionType), "Unknown plotting position type.");
            }
        }

        /// <summary>
        /// Returns the denominator of the weighted plotting positions: W + (1 − 2α) W / n, or W for Kaplan-Meier.
        /// </summary>
        /// <param name="total">The total weight W, positive and finite.</param>
        /// <param name="n">The number of points.</param>
        /// <param name="alpha">The alpha coefficient (unused for Kaplan-Meier).</param>
        /// <param name="kaplanMeier">True for the Kaplan-Meier positions c_i / W.</param>
        /// <returns>The denominator.</returns>
        private static double WeightedPositionDenominator(double total, int n, double alpha, bool kaplanMeier)
        {
            return kaplanMeier ? total : total + (1d - 2d * alpha) * (total / n);
        }

        /// <summary>
        /// Interpolates ascending data at probability k between weighted plotting positions formed during the walk,
        /// resuming from a cursor so that ascending probabilities share one walk.
        /// </summary>
        /// <param name="x">The ascending data.</param>
        /// <param name="w">The weights, aligned with <paramref name="x"/>.</param>
        /// <param name="denominator">The position denominator of <see cref="WeightedPositionDenominator(double, int, double, bool)"/>.</param>
        /// <param name="k">The probability, in [0,1], at least the probability of the previous call on the same cursor.</param>
        /// <param name="share">The share 1 − α of a point's own weight below its position (1 for Kaplan-Meier).</param>
        /// <param name="index">The cursor: the first point whose position has not been passed; 0 to start.</param>
        /// <param name="before">The cursor: the sequential sum of the weights before <paramref name="index"/>; 0 to start.</param>
        /// <param name="previousPosition">The cursor: the position of the point before <paramref name="index"/>; NaN to start.</param>
        /// <param name="previousValue">The cursor: the value of the point before <paramref name="index"/>; NaN to start.</param>
        /// <returns>The percentile, as documented on
        /// <see cref="Percentile(IList{double}, double, IList{double}, PlottingPositions.PlottingPostionType, bool)"/>.</returns>
        /// <remarks>
        /// The walk stops at the first point whose position is at least k and leaves the cursor on it, so a later call
        /// with a larger k recomputes that position with the same operations and continues from there. Every result
        /// is therefore bitwise the one a walk from the first point gives.
        /// </remarks>
        private static double WeightedPlottingPositionPercentile(double[] x, double[] w, double denominator, double k, double share,
            ref int index, ref double before, ref double previousPosition, ref double previousValue)
        {
            int n = x.Length;
            for (; index < n; index++)
            {
                double position = (before + share * w[index]) / denominator;
                if (k <= position)
                {
                    if (index == 0) return x[0];
                    double span = position - previousPosition;
                    if (span <= 0d) return x[index];
                    return previousValue + (x[index] - previousValue) * (k - previousPosition) / span;
                }
                previousPosition = position;
                previousValue = x[index];
                before += w[index];
            }
            return x[n - 1];
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
