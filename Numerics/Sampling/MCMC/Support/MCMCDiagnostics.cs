using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Sampling.MCMC
{
    /// <summary>
    /// A class for assessing Bayesian MCMC convergence diagnostics.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class MCMCDiagnostics
    {
        /// <summary>
        /// Blom offset used by the rank-normal inverse transformation.
        /// </summary>
        private const double RankOffset = 3d / 8d;

        /// <summary>
        /// IEEE-754 double-precision machine epsilon used by R <c>posterior</c>
        /// to identify constant diagnostic inputs.
        /// </summary>
        private const double MachineEpsilon = 2.2204460492503131E-16d;

        /// <summary>
        /// Computes a conservative rank-normalized effective sample size.
        /// </summary>
        /// <param name="series">The series of posterior samples to evaluate.</param>
        /// <returns>
        /// The minimum of bulk, lower-tail, and upper-tail effective sample size,
        /// or <see cref="double.NaN"/> when the series is insufficient or degenerate.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="series"/> is null.</exception>
        /// <remarks>
        /// The series is split into two half-chains. Bulk ESS uses pooled rank-normalized
        /// draws; tail ESS uses indicators at the pooled 5th and 95th percentiles. Each
        /// component uses Geyer's multi-chain initial-positive and initial-monotone sequence.
        /// Reference: Vehtari et al. (2021), Bayesian Analysis 16(2), 667-718.
        /// </remarks>
        public static double EffectiveSampleSize(IList<double> series)
        {
            if (series == null)
                throw new ArgumentNullException(nameof(series));

            return ComputeConservativeEffectiveSampleSize(new[] { series.ToArray() });
        }

        /// <summary>
        /// Computes a conservative rank-normalized effective sample size for each model parameter.
        /// </summary>
        /// <param name="markovChains">
        /// The Markov chains to evaluate. When lengths differ, only the common leading
        /// length of every chain is used.
        /// </param>
        /// <param name="averageACF">Output. A jagged array of averaged autocorrelation functions, one for each parameter.</param>
        /// <returns>
        /// One conservative scalar ESS per parameter: the minimum of bulk, lower-tail,
        /// and upper-tail ESS. Unusable inputs return <see cref="double.NaN"/>.
        /// </returns>
        /// <exception cref="ArgumentException">Thrown when no chains are provided or a chain is empty.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when no model parameters are present.</exception>
        /// <remarks>
        /// The 51-row averaged original-scale autocorrelation output is retained for plotting.
        /// It is not used to compute the modern scalar ESS.
        /// </remarks>
        public static double[] EffectiveSampleSize(IList<List<ParameterSet>> markovChains, out double[][,] averageACF)
        {
            if (markovChains == null || markovChains.Count == 0)
                throw new ArgumentException("No chains provided.", nameof(markovChains));
            if (markovChains.Any(chain => chain == null || chain.Count == 0))
                throw new ArgumentException("Every chain must contain at least one iteration.", nameof(markovChains));

            int chainCount = markovChains.Count;
            int parameterCount = markovChains[0][0].Values.Length;
            int commonLength = markovChains.Min(chain => chain.Count);
            if (parameterCount < 1)
                throw new ArgumentOutOfRangeException(nameof(markovChains), "There must be at least one parameter to evaluate.");

            var effectiveSampleSizes = new double[parameterCount];
            averageACF = new double[parameterCount][,];
            for (int parameterIndex = 0; parameterIndex < parameterCount; parameterIndex++)
            {
                var valuesByChain = new double[chainCount][];
                averageACF[parameterIndex] = new double[51, 2];
                for (int chainIndex = 0; chainIndex < chainCount; chainIndex++)
                {
                    var values = new double[commonLength];
                    for (int iteration = 0; iteration < commonLength; iteration++)
                        values[iteration] = markovChains[chainIndex][iteration].Values[parameterIndex];
                    valuesByChain[chainIndex] = values;
                    AccumulateAverageAutocorrelation(values, averageACF[parameterIndex], chainCount);
                }

                effectiveSampleSizes[parameterIndex] =
                    ComputeConservativeEffectiveSampleSize(valuesByChain);
            }

            return effectiveSampleSizes;
        }

        /// <summary>
        /// Computes rank-normalized split and folded split R-hat.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Each retained chain is split in half and pooled ranks are transformed to normal
        /// scores. The returned value is the maximum of rank-normalized split R-hat and
        /// folded rank-normalized split R-hat, making it sensitive to both location and scale.
        /// </para>
        /// <para>
        /// A single original chain returns <see cref="double.NaN"/> for compatibility and
        /// because independent-chain comparison is unavailable. Reference: Vehtari et al.
        /// (2021), Bayesian Analysis 16(2), 667-718.
        /// </para>
        /// </remarks>
        /// <param name="markovChains">The list of Markov Chains to be evaluated. The chains must be of equal length.</param>
        /// <param name="warmupIterations">The number of warm up MCMC iterations to discard at the beginning of the chains.</param>
        /// <returns>One modern R-hat value per model parameter.</returns>
        /// <exception cref="ArgumentException">Thrown when no chains are provided or chain lengths differ.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown for an invalid warmup count or parameter count.</exception>
        public static double[] GelmanRubin(IList<List<ParameterSet>> markovChains, int warmupIterations = 0)
        {
            if (markovChains == null || markovChains.Count == 0)
                throw new ArgumentException("No chains provided.", nameof(markovChains));
            if (markovChains.Any(chain => chain == null || chain.Count == 0))
                throw new ArgumentException("Every chain must contain at least one iteration.", nameof(markovChains));

            int chainCount = markovChains.Count;
            int parameterCount = markovChains[0][0].Values.Length;
            int iterationCount = markovChains[0].Count;
            foreach (var chain in markovChains)
            {
                if (chain.Count != iterationCount)
                    throw new ArgumentException("All chains must have the same length.");
            }

            if (parameterCount < 1)
                throw new ArgumentOutOfRangeException(nameof(markovChains), "There must be at least one parameter to evaluate.");
            if (warmupIterations < 0 || warmupIterations >= iterationCount)
                throw new ArgumentOutOfRangeException(nameof(warmupIterations), "Warmup iterations must leave at least one retained iteration.");

            var rhat = Enumerable.Repeat(double.NaN, parameterCount).ToArray();
            int retainedCount = iterationCount - warmupIterations;
            if (chainCount < 2 || retainedCount < 4)
                return rhat;

            for (int parameterIndex = 0; parameterIndex < parameterCount; parameterIndex++)
            {
                var retained = new double[chainCount][];
                for (int chainIndex = 0; chainIndex < chainCount; chainIndex++)
                {
                    retained[chainIndex] = new double[retainedCount];
                    for (int iteration = 0; iteration < retainedCount; iteration++)
                    {
                        retained[chainIndex][iteration] =
                            markovChains[chainIndex][warmupIterations + iteration].Values[parameterIndex];
                    }
                }

                if (ShouldReturnNaN(retained))
                    continue;

                double rankRhat = ComputeBasicRhat(RankNormalize(SplitChains(retained)));
                double foldedRhat = ComputeBasicRhat(
                    RankNormalize(SplitChains(FoldAroundMedian(retained))));
                if (Tools.IsFinite(rankRhat) && Tools.IsFinite(foldedRhat))
                    rhat[parameterIndex] = Math.Max(rankRhat, foldedRhat);
            }

            return rhat;
        }

        /// <summary>
        /// Adds one chain's original-scale autocorrelation to the retained plotting output.
        /// </summary>
        /// <param name="values">Common-length parameter draws for one chain.</param>
        /// <param name="averageAutocorrelation">The 51-row output accumulator.</param>
        /// <param name="chainCount">Number of chains contributing to the average.</param>
        private static void AccumulateAverageAutocorrelation(
            double[] values,
            double[,] averageAutocorrelation,
            int chainCount)
        {
            if (values.Length < 2)
                return;

            int maximumLag = Math.Min(50, values.Length - 1);
            double[,]? autocorrelation = Fourier.Autocorrelation(values, maximumLag);
            if (autocorrelation == null)
                return;

            for (int lag = 0; lag < autocorrelation.GetLength(0); lag++)
                averageAutocorrelation[lag, 1] += autocorrelation[lag, 1] / chainCount;
        }

        /// <summary>
        /// Computes the minimum of rank-normalized bulk and two quantile-indicator ESS values.
        /// </summary>
        /// <param name="chains">Original-scale chains with a common length.</param>
        /// <returns>The conservative scalar ESS, or <see cref="double.NaN"/>.</returns>
        private static double ComputeConservativeEffectiveSampleSize(double[][] chains)
        {
            if (chains.Length == 0 || chains[0].Length < 6 || ShouldReturnNaN(chains))
                return double.NaN;

            double bulk = ComputeGeyerEffectiveSampleSize(RankNormalize(SplitChains(chains)));
            double lowerTail = ComputeQuantileEffectiveSampleSize(chains, 0.05d);
            double upperTail = ComputeQuantileEffectiveSampleSize(chains, 0.95d);
            if (!Tools.IsFinite(bulk) || !Tools.IsFinite(lowerTail) || !Tools.IsFinite(upperTail))
                return double.NaN;

            return Math.Min(bulk, Math.Min(lowerTail, upperTail));
        }

        /// <summary>
        /// Computes ESS for an indicator at a pooled sample quantile.
        /// </summary>
        /// <param name="chains">Original-scale chains.</param>
        /// <param name="probability">Quantile probability.</param>
        /// <returns>The split-chain indicator ESS.</returns>
        private static double ComputeQuantileEffectiveSampleSize(double[][] chains, double probability)
        {
            double threshold = Quantile(Flatten(chains), probability);
            var indicators = new double[chains.Length][];
            for (int chainIndex = 0; chainIndex < chains.Length; chainIndex++)
            {
                indicators[chainIndex] = new double[chains[chainIndex].Length];
                for (int iteration = 0; iteration < chains[chainIndex].Length; iteration++)
                    indicators[chainIndex][iteration] = chains[chainIndex][iteration] <= threshold ? 1d : 0d;
            }

            return ComputeGeyerEffectiveSampleSize(SplitChains(indicators));
        }

        /// <summary>
        /// Computes multi-chain ESS using posterior's FFT autocovariance and Geyer sequence.
        /// </summary>
        /// <param name="chains">Transformed chains with a common length.</param>
        /// <returns>The effective sample size, which may exceed the draw count for negative correlation.</returns>
        private static double ComputeGeyerEffectiveSampleSize(double[][] chains)
        {
            int chainCount = chains.Length;
            int iterationCount = chains[0].Length;
            if (iterationCount < 3 || ShouldReturnNaN(chains))
                return double.NaN;

            var averageAutocovariance = new double[iterationCount];
            var chainMeans = new double[chainCount];
            for (int chainIndex = 0; chainIndex < chainCount; chainIndex++)
            {
                chainMeans[chainIndex] = chains[chainIndex].Average();
                double[] autocovariance = Autocovariance(chains[chainIndex], chainMeans[chainIndex]);
                for (int lag = 0; lag < iterationCount; lag++)
                    averageAutocovariance[lag] += autocovariance[lag] / chainCount;
            }

            double meanVariance = averageAutocovariance[0] * iterationCount / (iterationCount - 1d);
            double variancePlus = meanVariance * (iterationCount - 1d) / iterationCount;
            if (chainCount > 1)
                variancePlus += SampleVariance(chainMeans);
            if (!Tools.IsFinite(variancePlus) || variancePlus <= 0d)
                return double.NaN;

            var rho = new double[iterationCount];
            int sequenceIndex = 0;
            double rhoEven = 1d;
            rho[0] = rhoEven;
            double rhoOdd = 1d - (meanVariance - averageAutocovariance[1]) / variancePlus;
            rho[1] = rhoOdd;
            while (sequenceIndex < iterationCount - 5 &&
                   Tools.IsFinite(rhoEven + rhoOdd) &&
                   rhoEven + rhoOdd > 0d)
            {
                sequenceIndex += 2;
                rhoEven = 1d - (meanVariance - averageAutocovariance[sequenceIndex]) / variancePlus;
                rhoOdd = 1d - (meanVariance - averageAutocovariance[sequenceIndex + 1]) / variancePlus;
                if (rhoEven + rhoOdd >= 0d)
                {
                    rho[sequenceIndex] = rhoEven;
                    rho[sequenceIndex + 1] = rhoOdd;
                }
            }

            int maximumIndex = sequenceIndex;
            if (rhoEven > 0d)
                rho[maximumIndex] = rhoEven;

            sequenceIndex = 0;
            while (sequenceIndex <= maximumIndex - 4)
            {
                sequenceIndex += 2;
                double previousPair = rho[sequenceIndex - 2] + rho[sequenceIndex - 1];
                double currentPair = rho[sequenceIndex] + rho[sequenceIndex + 1];
                if (currentPair > previousPair)
                {
                    rho[sequenceIndex] = previousPair / 2d;
                    rho[sequenceIndex + 1] = rho[sequenceIndex];
                }
            }

            int sumLength = Math.Max(maximumIndex, 1);
            double rhoSum = 0d;
            for (int lag = 0; lag < sumLength; lag++)
                rhoSum += rho[lag];
            double integratedAutocorrelationTime =
                -1d + 2d * rhoSum + rho[maximumIndex];
            double totalDraws = chainCount * (double)iterationCount;
            double lowerBound = 1d / Math.Log10(totalDraws);
            if (integratedAutocorrelationTime < lowerBound)
                integratedAutocorrelationTime = lowerBound;

            return totalDraws / integratedAutocorrelationTime;
        }

        /// <summary>
        /// Computes biased lag autocovariances using a zero-padded FFT.
        /// </summary>
        /// <param name="values">One transformed chain.</param>
        /// <param name="mean">Precomputed chain mean.</param>
        /// <returns>Autocovariances for every nonnegative lag, divided by chain length.</returns>
        private static double[] Autocovariance(double[] values, double mean)
        {
            int transformLength = 1;
            while (transformLength < 2 * values.Length)
                transformLength <<= 1;

            var centered = new double[transformLength];
            for (int index = 0; index < values.Length; index++)
                centered[index] = values[index] - mean;
            double[] correlation = Fourier.Correlation(centered, centered);
            var autocovariance = new double[values.Length];
            for (int lag = 0; lag < values.Length; lag++)
                autocovariance[lag] = correlation[lag] / values.Length;
            return autocovariance;
        }

        /// <summary>
        /// Splits every chain in half, discarding the middle draw when its length is odd.
        /// </summary>
        /// <param name="chains">Original chains with a common length.</param>
        /// <returns>Twice as many chains, each with half the original length.</returns>
        private static double[][] SplitChains(double[][] chains)
        {
            int halfLength = chains[0].Length / 2;
            var split = new double[chains.Length * 2][];
            for (int chainIndex = 0; chainIndex < chains.Length; chainIndex++)
            {
                split[chainIndex] = new double[halfLength];
                split[chainIndex + chains.Length] = new double[halfLength];
                Array.Copy(chains[chainIndex], 0, split[chainIndex], 0, halfLength);
                Array.Copy(
                    chains[chainIndex],
                    chains[chainIndex].Length - halfLength,
                    split[chainIndex + chains.Length],
                    0,
                    halfLength);
            }
            return split;
        }

        /// <summary>
        /// Rank-normalizes pooled chain values using average ranks and Blom's offset.
        /// </summary>
        /// <param name="chains">Chains to rank as one pooled collection.</param>
        /// <returns>Chains with values replaced by standard-normal rank scores.</returns>
        private static double[][] RankNormalize(double[][] chains)
        {
            double[] flattened = Flatten(chains);
            var order = Enumerable.Range(0, flattened.Length)
                .OrderBy(index => flattened[index])
                .ToArray();
            var ranks = new double[flattened.Length];
            int start = 0;
            while (start < order.Length)
            {
                int end = start + 1;
                while (end < order.Length && flattened[order[end]].Equals(flattened[order[start]]))
                    end++;
                double averageRank = (start + 1d + end) / 2d;
                for (int index = start; index < end; index++)
                    ranks[order[index]] = averageRank;
                start = end;
            }

            var normalized = new double[chains.Length][];
            int flatIndex = 0;
            double denominator = flattened.Length - 2d * RankOffset + 1d;
            for (int chainIndex = 0; chainIndex < chains.Length; chainIndex++)
            {
                normalized[chainIndex] = new double[chains[chainIndex].Length];
                for (int iteration = 0; iteration < chains[chainIndex].Length; iteration++)
                {
                    double probability = (ranks[flatIndex] - RankOffset) / denominator;
                    normalized[chainIndex][iteration] = Normal.StandardZ(probability);
                    flatIndex++;
                }
            }
            return normalized;
        }

        /// <summary>
        /// Folds pooled draws around their median.
        /// </summary>
        /// <param name="chains">Original-scale chains.</param>
        /// <returns>Absolute deviations from the pooled median.</returns>
        private static double[][] FoldAroundMedian(double[][] chains)
        {
            double median = Quantile(Flatten(chains), 0.5d);
            var folded = new double[chains.Length][];
            for (int chainIndex = 0; chainIndex < chains.Length; chainIndex++)
            {
                folded[chainIndex] = new double[chains[chainIndex].Length];
                for (int iteration = 0; iteration < chains[chainIndex].Length; iteration++)
                    folded[chainIndex][iteration] = Math.Abs(chains[chainIndex][iteration] - median);
            }
            return folded;
        }

        /// <summary>
        /// Computes the basic between/within-chain R-hat for transformed split chains.
        /// </summary>
        /// <param name="chains">Transformed split chains.</param>
        /// <returns>The basic R-hat, or <see cref="double.NaN"/> for unusable input.</returns>
        private static double ComputeBasicRhat(double[][] chains)
        {
            if (chains.Length < 2 || chains[0].Length < 2 || ShouldReturnNaN(chains))
                return double.NaN;

            int iterationCount = chains[0].Length;
            var means = chains.Select(chain => chain.Average()).ToArray();
            double within = chains.Average(SampleVariance);
            if (!Tools.IsFinite(within) || within <= 0d)
                return double.NaN;
            double between = iterationCount * SampleVariance(means);
            return Math.Sqrt((between / within + iterationCount - 1d) / iterationCount);
        }

        /// <summary>
        /// Computes sample variance with denominator <c>n-1</c>.
        /// </summary>
        /// <param name="values">Values to evaluate.</param>
        /// <returns>The sample variance, or <see cref="double.NaN"/> for fewer than two values.</returns>
        private static double SampleVariance(double[] values)
        {
            if (values.Length < 2)
                return double.NaN;
            double mean = values.Average();
            double sum = 0d;
            for (int index = 0; index < values.Length; index++)
                sum += Tools.Sqr(values[index] - mean);
            return sum / (values.Length - 1d);
        }

        /// <summary>
        /// Flattens chains in column-major chain order matching an iterations-by-chains R matrix.
        /// </summary>
        /// <param name="chains">Chains to flatten.</param>
        /// <returns>All values with each chain contiguous.</returns>
        private static double[] Flatten(double[][] chains)
        {
            int total = chains.Sum(chain => chain.Length);
            var flattened = new double[total];
            int offset = 0;
            foreach (double[] chain in chains)
            {
                Array.Copy(chain, 0, flattened, offset, chain.Length);
                offset += chain.Length;
            }
            return flattened;
        }

        /// <summary>
        /// Computes the R type-7 sample quantile.
        /// </summary>
        /// <param name="values">Pooled sample values.</param>
        /// <param name="probability">Probability in the closed unit interval.</param>
        /// <returns>The interpolated sample quantile.</returns>
        private static double Quantile(double[] values, double probability)
        {
            var sorted = values.OrderBy(value => value).ToArray();
            double position = (sorted.Length - 1d) * probability;
            int lower = (int)Math.Floor(position);
            int upper = (int)Math.Ceiling(position);
            if (lower == upper)
                return sorted[lower];
            return sorted[lower] + (position - lower) * (sorted[upper] - sorted[lower]);
        }

        /// <summary>
        /// Determines whether pooled diagnostic input is nonfinite or effectively constant.
        /// </summary>
        /// <param name="chains">Chains to inspect.</param>
        /// <returns><see langword="true"/> when the diagnostic must return NaN.</returns>
        private static bool ShouldReturnNaN(double[][] chains)
        {
            double minimum = double.PositiveInfinity;
            double maximum = double.NegativeInfinity;
            foreach (double[] chain in chains)
            {
                foreach (double value in chain)
                {
                    if (!Tools.IsFinite(value))
                        return true;
                    minimum = Math.Min(minimum, value);
                    maximum = Math.Max(maximum, value);
                }
            }
            return Math.Abs(maximum - minimum) < MachineEpsilon;
        }

        /// <summary>
        /// Computes the minimum sample size rounded to the nearest 100 based on the Raftery-Lewis method.
        /// </summary>
        /// <param name="quantile">The posterior quantile of interest; e.g., 0.975.</param>
        /// <param name="tolerance">The acceptable tolerance for this quantile; e.g., ±0.005.</param>
        /// <param name="probability">Probability of being within the range of tolerance; e.g., 0.95.</param>
        /// <returns>The minimum sample size as an integer, rounded to the nearest thousands place.</returns>
        public static int MinimumSampleSize(double quantile, double tolerance, double probability)
        {
            double q = quantile;
            double r = tolerance;
            double s = probability;
            double N = (q * (1d - q) * Math.Pow(Normal.StandardZ(0.5d * (s + 1d)), 2d)) / Math.Pow(r, 2d);
            int Nmin = (int)Math.Round(N / 100d) * 100;
            return Nmin;
        }

    }

   
}
