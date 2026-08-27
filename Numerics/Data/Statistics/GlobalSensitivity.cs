using System;
using System.Collections.Generic;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Given-data global sensitivity estimators over paired input-output samples.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// These estimators require no special sampling design: they operate on any stored sample of
    /// one input against the output, conditioning through equal-frequency bins of the input. The
    /// first-order Sobol index measures the share of output variance explained by the input alone,
    /// the PAWN indices measure conditional-versus-unconditional distribution shifts through
    /// Kolmogorov-Smirnov statistics (sensitive to changes a variance ratio misses), and the
    /// Borgonovo delta is a moment-independent total-variation measure suited to tail-driven
    /// outputs. All three are deterministic: binning is by rank with ties resolved by the sort's
    /// deterministic order, and no randomness is used anywhere.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>Plischke, E., Borgonovo, E., and Smith, C. L. (2013). Global sensitivity measures from given data. European Journal of Operational Research, 226(3), 536-550.</item>
    /// <item>Pianosi, F., and Wagener, T. (2015). A simple and efficient method for global sensitivity analysis based on cumulative distribution functions. Environmental Modelling and Software, 67, 1-11.</item>
    /// <item>Borgonovo, E. (2007). A new uncertainty importance measure. Reliability Engineering and System Safety, 92(6), 771-784.</item>
    /// </list>
    /// </remarks>
    public static class GlobalSensitivity
    {

        /// <summary>
        /// Estimates the first-order Sobol index of an input from given data: the variance of the
        /// conditional output mean across equal-frequency input bins, over the output variance.
        /// </summary>
        /// <param name="x">The input sample.</param>
        /// <param name="y">The output sample, aligned with the input.</param>
        /// <param name="bins">Optional. The number of equal-frequency input bins. Default = 20.</param>
        /// <returns>The first-order index, clamped to [0, 1].</returns>
        /// <remarks>
        /// The estimator carries a positive finite-sample bias on the order of the bin count over
        /// the sample size, so an independent input reads near zero only when the sample is much
        /// larger than the bin count. A degenerate output with zero variance returns zero.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either sample is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the samples differ in length or are shorter than the bin count.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the bin count is less than two, or a sample value is not finite.</exception>
        public static double FirstOrderSobol(IList<double> x, IList<double> y, int bins = 20)
        {
            ValidateSamples(x, y, bins);
            int n = x.Count;
            var order = SortIndicesBy(x);

            double mean = 0;
            for (int i = 0; i < n; i++)
                mean += y[i];
            mean /= n;
            double totalVariance = 0;
            for (int i = 0; i < n; i++)
            {
                double d = y[i] - mean;
                totalVariance += d * d;
            }
            totalVariance /= n;
            if (totalVariance <= 0d) return 0d;

            double betweenVariance = 0;
            for (int b = 0; b < bins; b++)
            {
                int start = (int)((long)b * n / bins);
                int end = (int)((long)(b + 1) * n / bins);
                double binMean = 0;
                for (int i = start; i < end; i++)
                    binMean += y[order[i]];
                binMean /= end - start;
                double d = binMean - mean;
                betweenVariance += (end - start) * d * d;
            }
            betweenVariance /= n;
            return Tools.Clamp(betweenVariance / totalVariance, 0d, 1d);
        }

        /// <summary>
        /// Estimates the PAWN sensitivity statistics of an input from given data: for each
        /// equal-frequency input bin, the two-sample Kolmogorov-Smirnov statistic between the
        /// conditional output distribution in the bin and the unconditional output distribution.
        /// </summary>
        /// <param name="x">The input sample.</param>
        /// <param name="y">The output sample, aligned with the input.</param>
        /// <param name="bins">Optional. The number of equal-frequency input bins. Default = 20.</param>
        /// <returns>The per-bin Kolmogorov-Smirnov statistics, in bin order.</returns>
        /// <remarks>
        /// The caller summarizes the per-bin statistics as the application demands - the median is
        /// the customary index and <see cref="PawnMedian"/> provides it directly. An independent
        /// input yields statistics at the two-sample sampling-noise scale rather than exactly zero.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either sample is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the samples differ in length or are shorter than the bin count.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the bin count is less than two, or a sample value is not finite.</exception>
        public static double[] Pawn(IList<double> x, IList<double> y, int bins = 20)
        {
            ValidateSamples(x, y, bins);
            int n = x.Count;
            var order = SortIndicesBy(x);

            var sortedY = new double[n];
            for (int i = 0; i < n; i++)
                sortedY[i] = y[i];
            Array.Sort(sortedY);

            var result = new double[bins];
            for (int b = 0; b < bins; b++)
            {
                int start = (int)((long)b * n / bins);
                int end = (int)((long)(b + 1) * n / bins);
                int size = end - start;
                var binY = new double[size];
                for (int i = start; i < end; i++)
                    binY[i - start] = y[order[i]];
                Array.Sort(binY);
                result[b] = TwoSampleKolmogorovSmirnov(binY, sortedY);
            }
            return result;
        }

        /// <summary>
        /// Estimates the customary PAWN index of an input from given data: the median of the
        /// per-bin Kolmogorov-Smirnov statistics of <see cref="Pawn"/>.
        /// </summary>
        /// <param name="x">The input sample.</param>
        /// <param name="y">The output sample, aligned with the input.</param>
        /// <param name="bins">Optional. The number of equal-frequency input bins. Default = 20.</param>
        /// <returns>The median Kolmogorov-Smirnov statistic.</returns>
        /// <exception cref="ArgumentNullException">Thrown when either sample is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the samples differ in length or are shorter than the bin count.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the bin count is less than two, or a sample value is not finite.</exception>
        public static double PawnMedian(IList<double> x, IList<double> y, int bins = 20)
        {
            return Statistics.Percentile(Pawn(x, y, bins), 0.5d);
        }

        /// <summary>
        /// Estimates the Borgonovo delta of an input from given data with the double-histogram
        /// total-variation estimator: the output is partitioned into equal-frequency rank classes,
        /// and delta is half the bin-weighted total variation between the conditional and
        /// unconditional class frequencies.
        /// </summary>
        /// <param name="x">The input sample.</param>
        /// <param name="y">The output sample, aligned with the input.</param>
        /// <param name="xBins">Optional. The number of equal-frequency input bins. Default = 20.</param>
        /// <param name="yBins">Optional. The number of equal-frequency output classes. Default = 20.</param>
        /// <returns>The delta estimate, within [0, 1 - 1/yBins].</returns>
        /// <remarks>
        /// The class partition is rank-based, so the estimate is exactly invariant to strictly
        /// increasing transforms of the output - the defining property of the delta measure. A
        /// functionally dependent output saturates at 1 - 1/yBins, and an independent input carries
        /// a positive finite-sample bias on the order of the square root of the class count over the
        /// bin size, so small deltas require samples much larger than the histogram resolution.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either sample is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the samples differ in length or are shorter than either bin count.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when either bin count is less than two, or a sample value is not finite.</exception>
        public static double BorgonovoDelta(IList<double> x, IList<double> y, int xBins = 20, int yBins = 20)
        {
            ValidateSamples(x, y, xBins);
            if (yBins < 2) throw new ArgumentOutOfRangeException(nameof(yBins), "The class count must be at least two.");
            int n = x.Count;
            if (n < yBins) throw new ArgumentException("The sample must be at least as long as the class count.", nameof(yBins));
            var xOrder = SortIndicesBy(x);
            var yOrder = SortIndicesBy(y);

            // Assign each observation its equal-frequency output class by rank.
            var yClass = new int[n];
            var classCounts = new int[yBins];
            for (int c = 0; c < yBins; c++)
            {
                int start = (int)((long)c * n / yBins);
                int end = (int)((long)(c + 1) * n / yBins);
                classCounts[c] = end - start;
                for (int i = start; i < end; i++)
                    yClass[yOrder[i]] = c;
            }

            double delta = 0;
            var conditionalCounts = new int[yBins];
            for (int b = 0; b < xBins; b++)
            {
                int start = (int)((long)b * n / xBins);
                int end = (int)((long)(b + 1) * n / xBins);
                int size = end - start;
                Array.Clear(conditionalCounts, 0, yBins);
                for (int i = start; i < end; i++)
                    conditionalCounts[yClass[xOrder[i]]]++;
                double totalVariation = 0;
                for (int c = 0; c < yBins; c++)
                    totalVariation += Math.Abs((double)conditionalCounts[c] / size - (double)classCounts[c] / n);
                delta += (double)size / n * totalVariation;
            }
            return 0.5d * delta;
        }

        /// <summary>
        /// Validates a paired sample and bin count.
        /// </summary>
        /// <param name="x">The input sample.</param>
        /// <param name="y">The output sample.</param>
        /// <param name="bins">The bin count.</param>
        private static void ValidateSamples(IList<double> x, IList<double> y, int bins)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (x.Count != y.Count) throw new ArgumentException("The input and output samples must have the same length.", nameof(y));
            if (bins < 2) throw new ArgumentOutOfRangeException(nameof(bins), "The bin count must be at least two.");
            if (x.Count < bins) throw new ArgumentException("The sample must be at least as long as the bin count.", nameof(x));
            for (int i = 0; i < x.Count; i++)
            {
                if (!Tools.IsFinite(x[i]))
                    throw new ArgumentOutOfRangeException(nameof(x), "Sample values must be finite.");
                if (!Tools.IsFinite(y[i]))
                    throw new ArgumentOutOfRangeException(nameof(y), "Sample values must be finite.");
            }
        }

        /// <summary>
        /// Returns the sample indices sorted ascending by value. Ties keep the sort's deterministic
        /// order, so identical inputs always produce identical bin assignments.
        /// </summary>
        /// <param name="values">The values to rank.</param>
        /// <returns>The sorted index array.</returns>
        private static int[] SortIndicesBy(IList<double> values)
        {
            int n = values.Count;
            var keys = new double[n];
            var order = new int[n];
            for (int i = 0; i < n; i++)
            {
                keys[i] = values[i];
                order[i] = i;
            }
            Array.Sort(keys, order);
            return order;
        }

        /// <summary>
        /// Computes the exact two-sample Kolmogorov-Smirnov statistic between two sorted samples,
        /// evaluating the gap only after advancing both empirical distributions through tied values.
        /// </summary>
        /// <param name="first">The first sorted sample.</param>
        /// <param name="second">The second sorted sample.</param>
        /// <returns>The supremum absolute difference of the empirical distributions.</returns>
        private static double TwoSampleKolmogorovSmirnov(double[] first, double[] second)
        {
            int n1 = first.Length, n2 = second.Length;
            int i = 0, j = 0;
            double supremum = 0;
            while (i < n1 && j < n2)
            {
                double value = Math.Min(first[i], second[j]);
                while (i < n1 && first[i] == value) i++;
                while (j < n2 && second[j] == value) j++;
                double gap = Math.Abs((double)i / n1 - (double)j / n2);
                if (gap > supremum) supremum = gap;
            }
            // The remaining tail of either sample only shrinks toward equality at one, and the gap
            // at the crossover was already recorded, so the sweep is complete.
            return supremum;
        }

    }
}
