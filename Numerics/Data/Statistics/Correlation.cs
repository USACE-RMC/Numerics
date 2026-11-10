using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Contains methods for determining the correlation of two data sets.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This class contains a shared functions for Pearson and Spearman Correlation Coefficients. The Pearson correlation coefficient measures linear correlation
    /// between two sets of data. The Spearman correlation coefficient is a nonparametric measure of rank correlation (i.e. measures the strength and direction of
    /// association between two ranked variables).
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </item>
    /// <item>
    /// <see href="https://en.wikipedia.org/wiki/Pearson_correlation_coefficient"/>
    /// </item>
    /// <item>
    /// <see href="https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient"/>
    /// </item>
    /// </list>
    /// </remarks>
    public sealed class Correlation
    {
        /// <summary>
        /// From this number of columns on, the correlation matrix sums are vectorized across columns.
        /// </summary>
        private const int VectorizedColumns = 8;


        /// <summary>
        /// Computes the Pearson correlation coefficient.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2.</param>
        /// <returns>The Pearson correlation coefficient.</returns>
        /// <remarks>
        /// The sums of squares and cross products are accumulated by the corrected two-pass algorithm: each sum of
        /// products of deviations is reduced by the product of the deviation sums divided by n, which removes, to
        /// first order, the rounding errors of the computed means (Chan, T. F., Golub, G. H. and LeVeque, R. J., 1983,
        /// "Algorithms for computing the sample variance: analysis and recommendations", The American Statistician
        /// 37(3), 242-247).
        /// The covariance and variance members of <see cref="Statistics"/> use the same correction.
        /// </remarks>
        public static double Pearson(IList<double> sample1, IList<double> sample2)
        {
            // Check if the arrays are the same length
            if (sample2.Count != sample1.Count)
            {
                throw new ArgumentOutOfRangeException(nameof(sample1), "The samples arrays must be the same length.");
            }
            // 
            int i, n = sample1.Count;
            // Find means.
            double ax = 0d;
            double ay = 0d;
            for (i = 0; i <= n - 1; i++)
            {
                ax += sample1[i];
                ay += sample2[i];
            }
            ax /= n;
            ay /= n;
            // Compute the correlation coefficient
            double sxx = 0d, syy = 0d, sxy = 0d, sx = 0d, sy = 0d;
            double xt, yt;
            for (i = 0; i <= n - 1; i++)
            {
                xt = sample1[i] - ax;
                yt = sample2[i] - ay;
                sxx += xt * xt;
                syy += yt * yt;
                sxy += xt * yt;
                sx += xt;
                sy += yt;
            }
            sxx -= sx * sx / n;
            syy -= sy * sy / n;
            sxy -= sx * sy / n;
            return sxy / Math.Sqrt(sxx * syy);
        }

        /// <summary>
        /// Computes the Pearson correlation coefficient matrix for the variables in a 2D array.
        /// </summary>
        /// <param name="samples">
        /// A 2D array of size [n, p] where each row is one observation and each column is one variable.
        /// </param>
        /// <returns>
        /// A p×p matrix C where C[j,k] is the Pearson correlation between column j and column k of <paramref name="samples"/>.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown if <paramref name="samples"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown if <paramref name="samples"/> has zero columns.</exception>
        /// <remarks>
        /// From eight columns on, the sums are vectorized across columns; every entry is accumulated over the rows in
        /// the same order as the scalar loops, so the result does not depend on the number of columns. The sums are
        /// corrected as in <see cref="Pearson(IList{double}, IList{double})"/>.
        /// </remarks>
        public static double[,] Pearson(double[,] samples)
        {
            if (samples == null)
                throw new ArgumentNullException(nameof(samples));

            int n = samples.GetLength(0);
            int p = samples.GetLength(1);
            if (p == 0)
                throw new ArgumentException("Input must have at least one column.", nameof(samples));

            // 1) Compute means for each column
            double[] means = new double[p];
            // 2) Compute sum of squares (ss) and cross-products (cov)
            double[] ss = new double[p];
            double[] deviationSums = new double[p];
            double[,] cov = new double[p, p];
            if (p < VectorizedColumns)
            {
                for (int j = 0; j < p; j++)
                    for (int i = 0; i < n; i++)
                        means[j] += samples[i, j];
                for (int j = 0; j < p; j++)
                    means[j] /= n;

                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < p; j++)
                    {
                        double dx = samples[i, j] - means[j];
                        ss[j] += dx * dx;
                        deviationSums[j] += dx;
                        for (int k = j; k < p; k++)
                        {
                            double dy = samples[i, k] - means[k];
                            cov[j, k] += dx * dy;
                        }
                    }
                }
            }
            else
            {
                // The same sums in the same order for every entry, vectorized across columns: the column sums add
                // the rows in turn, and each row is centered once and added to the rows of the cross-product matrix.
                var row = new double[p];
                for (int i = 0; i < n; i++)
                {
                    Buffer.BlockCopy(samples, i * p * sizeof(double), row, 0, p * sizeof(double));
                    CrossProductKernel.Axpy(1d, row, 0, means, 0, p);
                }
                for (int j = 0; j < p; j++)
                    means[j] /= n;

                var flat = new double[p * p];
                for (int i = 0; i < n; i++)
                {
                    for (int k = 0; k < p; k++)
                        row[k] = samples[i, k] - means[k];
                    CrossProductKernel.Axpy(1d, row, 0, deviationSums, 0, p);
                    for (int j = 0; j < p; j++)
                    {
                        double dx = row[j];
                        ss[j] += dx * dx;
                        CrossProductKernel.Axpy(dx, row, j, flat, j * p + j, p - j);
                    }
                }
                Buffer.BlockCopy(flat, 0, cov, 0, flat.Length * sizeof(double));
            }
            // Remove the first-order rounding error of the means, then mirror the upper triangle to the lower.
            for (int j = 0; j < p; j++)
            {
                ss[j] -= deviationSums[j] * deviationSums[j] / n;
                for (int k = j; k < p; k++)
                    cov[j, k] -= deviationSums[j] * deviationSums[k] / n;
            }
            for (int j = 0; j < p; j++)
                for (int k = j + 1; k < p; k++)
                    cov[k, j] = cov[j, k];

            // 3) Build correlation matrix
            double[,] corr = new double[p, p];
            for (int j = 0; j < p; j++)
            {
                for (int k = 0; k < p; k++)
                {
                    corr[j, k] = cov[j, k] / Math.Sqrt(ss[j] * ss[k]);
                }
            }

            return corr;
        }

        /// <summary>
        /// Computes the Spearman ranked correlation coefficient.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2.</param>
        /// <returns>The Spearman ranked correlation coefficient.</returns>
        public static double Spearman(IList<double> sample1, IList<double> sample2)
        {
            // Check if the arrays are the same length
            if (sample2.Count != sample1.Count)
            {
                throw new ArgumentOutOfRangeException(nameof(sample1), "The sample arrays must be the same length.");
            }

            var rank1 = Statistics.RanksInPlace(sample1.ToArray());
            var rank2 = Statistics.RanksInPlace(sample2.ToArray());
            return Pearson(rank1, rank2);
        }

        /// <summary>
        /// Computes the Spearman correlation coefficient matrix for the variables in a 2D array.
        /// </summary>
        /// <param name="samples">
        /// A 2D array of size [n, p] where each row is one observation and each column is one variable.
        /// </param>
        /// <returns>
        /// A p×p matrix C where C[j,k] is the Spearman correlation between column j and column k of <paramref name="samples"/>.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown if <paramref name="samples"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown if <paramref name="samples"/> has zero columns.</exception>
        public static double[,] Spearman(double[,] samples)
        {
            if (samples == null)
                throw new ArgumentNullException(nameof(samples));

            int n = samples.GetLength(0);
            int p = samples.GetLength(1);
            if (p == 0)
                throw new ArgumentException("Input must have at least one column.", nameof(samples));

            var ranks = new double[n, p];
            for (int i = 0; i < p; i++)
            {
                ranks.SetColumn(i, Statistics.RanksInPlace(samples.GetColumn(i)));
            }
            return Pearson(ranks);
        }

        /// <summary>
        /// Computes Kendall's Tau ranked correlation coefficient.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2.</param>
        /// <returns>The Kendall's Tau-b ranked correlation coefficient, or NaN when either sample contains NaN or
        /// either sample is constant.</returns>
        /// <remarks>
        /// <para>
        /// Tau-b is (n_c − n_d) / (√n_1 √n_2), with n_c − n_d the concordant minus discordant pairs among the pairs
        /// tied in neither sample, and n_1 and n_2 the pairs not tied in the first and second sample. Two values are
        /// tied when they compare equal, and a pair's concordance is decided by comparing the values, so infinite
        /// values and values of any magnitude are ranked exactly. A NaN in either sample gives NaN, as
        /// scipy.stats.kendalltau does by default.
        /// </para>
        /// <para>
        /// Samples of 48 or more values obtain the three integer counts in O(n log n) by Knight's (1966) algorithm:
        /// sort the pairs by the first sample (ties by the second), then count the inversions of the second sample
        /// with a merge sort. Smaller samples use the O(n²) pair loop, which gives the same counts. The counts are
        /// 64-bit, so any sample size is supported.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Knight, W. R. (1966). A computer method for calculating Kendall's tau with ungrouped data. Journal of
        /// the American Statistical Association, 61(314), 436-439.
        /// </para>
        /// </remarks>
        public static double KendallsTau(IList<double> sample1, IList<double> sample2)
        {
            // Check if the arrays are the same length
            if (sample2.Count != sample1.Count)
            {
                throw new ArgumentOutOfRangeException(nameof(sample1), "The sample arrays must be the same length.");
            }

            int n = sample1.Count;
            for (int r = 0; r < n; r++)
            {
                if (double.IsNaN(sample1[r]) || double.IsNaN(sample2[r])) return double.NaN;
            }
            long i, n1, n2;
            if (n >= KnightMinimumSize)
                KnightCounts(sample1, sample2, out i, out n1, out n2);
            else
                PairCounts(sample1, sample2, out i, out n1, out n2);
            return i / (Math.Sqrt(n1) * Math.Sqrt(n2));
        }

        /// <summary>
        /// The smallest sample size for which <see cref="KendallsTau"/> uses Knight's algorithm. Below it the
        /// pair loop is as fast; both give the same counts.
        /// </summary>
        internal const int KnightMinimumSize = 48;

        /// <summary>
        /// Computes Kendall's tau-b by the O(n²) pair loop, for samples without NaN.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2, of the same length.</param>
        /// <returns>The Kendall's Tau ranked correlation coefficient.</returns>
        internal static double KendallsTauPairs(IList<double> sample1, IList<double> sample2)
        {
            PairCounts(sample1, sample2, out long i, out long n1, out long n2);
            return i / (Math.Sqrt(n1) * Math.Sqrt(n2));
        }

        /// <summary>
        /// Counts the pairs of tau-b by comparing every pair of observations.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2, of the same length, without NaN.</param>
        /// <param name="score">Output. Concordant minus discordant pairs among the pairs tied in neither sample.</param>
        /// <param name="untied1">Output. The pairs not tied in sample 1.</param>
        /// <param name="untied2">Output. The pairs not tied in sample 2.</param>
        private static void PairCounts(IList<double> sample1, IList<double> sample2, out long score, out long untied1, out long untied2)
        {
            int n = sample1.Count;
            long i = 0, n1 = 0, n2 = 0;
            // Loop over first member of pair and second member.
            for (int j = 0; j <= n - 2; j++)
            {
                for (int k = j + 1; k <= n - 1; k++)
                {
                    int a1 = CompareSign(sample1[j], sample1[k]);
                    int a2 = CompareSign(sample2[j], sample2[k]);
                    if (a1 != 0 && a2 != 0)
                    {
                        // Neither has a tie
                        n1 += 1;
                        n2 += 1;
                        i += a1 == a2 ? 1 : -1;
                    }
                    else
                    {
                        // One or both arrays have ties.
                        if (a1 != 0) n1 += 1;
                        if (a2 != 0) n2 += 1;
                    }
                }
            }
            score = i;
            untied1 = n1;
            untied2 = n2;
        }

        /// <summary>
        /// Returns the sign of a − b by comparison: −1, 0 or +1.
        /// </summary>
        /// <param name="a">The first value, not NaN.</param>
        /// <param name="b">The second value, not NaN.</param>
        /// <returns>The sign of a − b, exact for infinite values and for values whose difference underflows or overflows.</returns>
        private static int CompareSign(double a, double b)
        {
            return a > b ? 1 : a < b ? -1 : 0;
        }

        /// <summary>
        /// Counts the pairs of tau-b by Knight's algorithm.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2, of the same length, without NaN.</param>
        /// <param name="score">Output. Concordant minus discordant pairs among the pairs tied in neither sample.</param>
        /// <param name="untied1">Output. The pairs not tied in sample 1.</param>
        /// <param name="untied2">Output. The pairs not tied in sample 2.</param>
        /// <remarks>
        /// With n_0 = n (n − 1) / 2, T_1 and T_2 the pairs tied in each sample, T_12 the pairs tied in both, and S
        /// the strict inversions of sample 2 once the pairs are ordered by sample 1 then sample 2, the counts are
        /// n_1 = n_0 − T_1, n_2 = n_0 − T_2 and n_c − n_d = n_0 − T_1 − T_2 + T_12 − 2S.
        /// </remarks>
        private static void KnightCounts(IList<double> sample1, IList<double> sample2, out long score, out long untied1, out long untied2)
        {
            int n = sample1.Count;
            var x = new double[n];
            var y = new double[n];
            for (int r = 0; r < n; r++)
            {
                x[r] = sample1[r];
                y[r] = sample2[r];
            }

            // Order the pairs by sample 1, then by sample 2 within each run of ties in sample 1.
            Array.Sort(x, y);
            long tied1 = 0, tied12 = 0;
            int start = 0;
            while (start < n)
            {
                int end = start + 1;
                while (end < n && x[end] == x[start]) end++;
                long run = end - start;
                tied1 += run * (run - 1) / 2;
                if (run > 1)
                {
                    Array.Sort(y, start, end - start);
                    tied12 += TiedPairs(y, start, end);
                }
                start = end;
            }

            long swaps = MergeSortInversions(y, new double[n]);
            long tied2 = TiedPairs(y, 0, n);
            long total = (long)n * (n - 1) / 2;
            untied1 = total - tied1;
            untied2 = total - tied2;
            score = total - tied1 - tied2 + tied12 - 2 * swaps;
        }

        /// <summary>
        /// Counts the pairs of equal values in an ascending range.
        /// </summary>
        /// <param name="values">The ascending values.</param>
        /// <param name="start">The first index of the range.</param>
        /// <param name="end">One past the last index of the range.</param>
        /// <returns>The sum of t (t − 1) / 2 over the runs of t equal values.</returns>
        private static long TiedPairs(double[] values, int start, int end)
        {
            long tied = 0;
            int r0 = start;
            while (r0 < end)
            {
                int r1 = r0 + 1;
                while (r1 < end && values[r1] == values[r0]) r1++;
                long run = r1 - r0;
                tied += run * (run - 1) / 2;
                r0 = r1;
            }
            return tied;
        }

        /// <summary>
        /// Sorts the values ascending with a bottom-up merge sort and counts the strict inversions removed.
        /// </summary>
        /// <param name="values">The values, sorted in place.</param>
        /// <param name="work">A work array of the same length.</param>
        /// <returns>The number of pairs (a, b), a before b, with values[b] &lt; values[a].</returns>
        private static long MergeSortInversions(double[] values, double[] work)
        {
            int n = values.Length;
            long swaps = 0;
            for (int width = 1; width < n; width *= 2)
            {
                for (int lo = 0; lo < n - width; lo += 2 * width)
                {
                    int mid = lo + width;
                    int hi = Math.Min(lo + 2 * width, n);
                    int a = lo, b = mid, k = lo;
                    while (a < mid && b < hi)
                    {
                        if (values[b] < values[a])
                        {
                            work[k++] = values[b++];
                            swaps += mid - a;
                        }
                        else
                        {
                            work[k++] = values[a++];
                        }
                    }
                    while (a < mid) work[k++] = values[a++];
                    while (b < hi) work[k++] = values[b++];
                    Array.Copy(work, lo, values, lo, hi - lo);
                }
            }
            return swaps;
        }
      
    }
}