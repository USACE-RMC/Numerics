using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// A class for performing the bootstrap uncertainty analysis.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Bootstrapping_(statistics)" />
    /// </para>
    /// </remarks>
    public class BootstrapAnalysis
    {

        #region Construction

        /// <summary>
        /// Construct a new Bootstrap Analysis.
        /// </summary>
        /// <param name="distribution">The univariate distribution to bootstrap.</param>
        /// <param name="estimationMethod">The parameter estimation method.</param>
        /// <param name="sampleSize">Size of the bootstrap sample to generate.</param>
        /// <param name="replications">The number of bootstrap replications to be sampled.</param>
        /// <param name="seed">Optional. Seed for random number generator. Default = 12345.</param>
        public BootstrapAnalysis(IUnivariateDistribution distribution, ParameterEstimationMethod estimationMethod, int sampleSize, int replications = 10000, int seed = 12345)
        {
            if (distribution as IBootstrappable == null) throw new ArgumentException("The distribution must implement IBootstrappable.", nameof(Distribution));
            if (sampleSize < 10) throw new ArgumentOutOfRangeException(nameof(SampleSize), "The sample size must at least 10.");
            if (replications < 100) throw new ArgumentOutOfRangeException(nameof(Replications), "The number of bootstrap replications must be at least 100.");

            Distribution = (IBootstrappable)distribution;
            EstimationMethod = estimationMethod;
            SampleSize = sampleSize;
            Replications = replications;
            PRNGSeed = seed;
        }

        #endregion

        #region Members

        private int _retries = 20;

        /// <summary>
        /// The number of accumulation chunks used by the reductions over replications. Fixed, not
        /// derived from the processor count, so the summation order — and therefore the result —
        /// does not vary with the machine or the thread count.
        /// </summary>
        private const int ReductionChunks = 64;

        /// <summary>
        /// The distribution parameter estimation method.
        /// </summary>
        public ParameterEstimationMethod EstimationMethod { get; private set; }

        /// <summary>
        /// The number of replications whose parameter fit failed on the most recent call to
        /// <see cref="Distributions()"/> or <see cref="Distributions(ParameterSet[])"/>. Failed
        /// replications are recorded as null and excluded from every summary, so a non-zero count
        /// means the results rest on fewer samples than requested.
        /// </summary>
        public int FailedReplications { get; private set; }

        /// <summary>
        /// The univariate distribution to bootstrap.
        /// </summary>
        public IBootstrappable Distribution { get; private set; }

        /// <summary>
        /// Size of the bootstrap sample to generate.
        /// </summary>
        public int SampleSize { get; private set; }

        /// <summary>
        /// The number of bootstrap replications to be sampled.
        /// </summary>
        public int Replications { get; private set; }

        /// <summary>
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; private set; }

        #endregion

        #region Methods

        /// <summary>
        /// Bootstrap a list of fitted distributions.
        /// </summary>
        public IUnivariateDistribution[] Distributions()
        {
            var bootDistributions = new IUnivariateDistribution[Replications];
            var r = new MersenneTwister(PRNGSeed);
            var seeds = r.NextIntegers(Replications);
            int failures = 0;
            Parallel.For(0, Replications, idx =>
            {
                bool failed = false;
                for (int m = 0; m < _retries; m++)
                {
                    try
                    {
                        bootDistributions[idx] = Distribution.Bootstrap(EstimationMethod, SampleSize, seeds[idx] + 10 * m);
                        failed = false;
                    }
                    catch (Exception)
                    {
                        failed = true;
                    };

                    if (failed == false) break;
                }

                // MLE and certain L-moments methods can fail to find a solution
                // On fail, set to null
                if (failed == true)
                {
                    bootDistributions[idx] = null!;
                    Interlocked.Increment(ref failures);
                }

            });
            FailedReplications = failures;
            return bootDistributions;
        }

        /// <summary>
        /// Return a list of distributions given an array of parameter sets.
        /// </summary>
        /// <param name="parameterSets">An array of parameter sets.</param>
        public IUnivariateDistribution[] Distributions(ParameterSet[] parameterSets)
        {
            var bootDistributions = new IUnivariateDistribution[parameterSets.Length];
            int failures = 0;
            Parallel.For(0, parameterSets.Length, idx =>
            {
                bool failed = false;

                try
                {
                    var dist = ((UnivariateDistributionBase)Distribution).Clone();
                    dist.SetParameters(parameterSets[idx].Values);
                    bootDistributions[idx] = dist;
                    failed = false;
                }
                catch (Exception)
                {
                    failed = true;
                };

                // On fail, set to null
                if (failed == true)
                {
                    bootDistributions[idx] = null!;
                    Interlocked.Increment(ref failures);
                }


            });
            FailedReplications = failures;
            return bootDistributions;
        }


        /// <summary>
        /// Bootstrap an array of distribution parameters.
        /// </summary>
        public double[,] Parameters(IUnivariateDistribution[]? distributions = null)
        {
            var bootDistributions = distributions != null ? distributions : Distributions();
            var bootParameters = new double[bootDistributions.Length, Distribution.NumberOfParameters];
            Parallel.For(0, bootDistributions.Length, idx =>
            {
                if (bootDistributions[idx] != null)
                {
                    var parameters = bootDistributions[idx].GetParameters;
                    for (int i = 0; i < Distribution.NumberOfParameters; i++)
                        bootParameters[idx, i] = parameters[i];
                }
                else
                {
                    for (int i = 0; i < Distribution.NumberOfParameters; i++)
                        bootParameters[idx, i] = double.NaN;
                }

            });
            return bootParameters;
        }

        /// <summary>
        /// Bootstrap an array of distribution parameter sets.
        /// </summary>
        public ParameterSet[] ParameterSets(IUnivariateDistribution[]? distributions = null)
        {
            var bootDistributions = distributions != null ? distributions : Distributions();
            var bootParameters = new ParameterSet[bootDistributions.Length];
            Parallel.For(0, bootDistributions.Length, idx =>
            {
                if (bootDistributions[idx] != null)
                {
                    bootParameters[idx] = new ParameterSet(bootDistributions[idx].GetParameters, double.NaN);
                }
                else
                {
                    var parameters = new double[Distribution.NumberOfParameters];
                    for (int i = 0; i < Distribution.NumberOfParameters; i++)
                        parameters[i] = double.NaN;
                    bootParameters[idx] = new ParameterSet(parameters, double.NaN);
                }

            });
            return bootParameters;
        }

        /// <summary>
        /// Bootstrap a list of product moments for each bootstrapped sample.
        /// </summary>
        public double[,] ProductMoments()
        {
            var bootMoments = new double[Replications, 4];
            var r = new MersenneTwister(PRNGSeed);
            var seeds = r.NextIntegers(Replications);
            Parallel.For(0, Replications, idx =>
            {
                var moments = Statistics.ProductMoments(Distribution.GenerateRandomValues(SampleSize, seeds[idx]));
                for (int i = 0; i < moments.Length; i++) bootMoments[idx, i] = moments[i];
            });
            return bootMoments;
        }

        /// <summary>
        /// Bootstrap a list of linear moments for each bootstrapped sample.
        /// </summary>
        public double[,] LinearMoments()
        {
            var bootMoments = new double[Replications, 4];
            var r = new MersenneTwister(PRNGSeed);
            var seeds = r.NextIntegers(Replications);
            Parallel.For(0, Replications, idx =>
            {
                var moments = Statistics.LinearMoments(Distribution.GenerateRandomValues(SampleSize, seeds[idx]));
                for (int i = 0; i < moments.Length; i++) bootMoments[idx, i] = moments[i];
            });
            return bootMoments;
        }

        /// <summary>
        /// Bootstrap a list of quantiles given the input non-exceedance probabilities.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[,] Quantiles(IList<double> probabilities, IUnivariateDistribution[]? distributions = null)
        {
            var bootDistributions = distributions != null ? distributions : Distributions();
            var Output = new double[bootDistributions.Length, probabilities.Count];
            Parallel.For(0, bootDistributions.Length, idx =>
            {
                var distribution = bootDistributions[idx];
                for (int i = 0; i < probabilities.Count; i++)
                    Output[idx, i] = distribution != null ? distribution.InverseCDF(probabilities[i]) : double.NaN;
            });
            return Output;
        }

        /// <summary>
        /// Bootstrap a list of non-exceedance probabilities given the input quantile values.
        /// </summary>
        /// <param name="quantiles">List quantile values.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[,] Probabilities(IList<double> quantiles, IUnivariateDistribution[]? distributions = null)
        {
            var bootDistributions = distributions != null ? distributions : Distributions();
            var Output = new double[bootDistributions.Length, quantiles.Count];
            Parallel.For(0, bootDistributions.Length, idx =>
            {
                var distribution = bootDistributions[idx];
                for (int i = 0; i < quantiles.Count; i++)
                    Output[idx, i] = distribution != null ? distribution.CDF(quantiles[i]) : double.NaN;
            });
            return Output;
        }

        /// <summary>
        /// Bootstrap full uncertainty analysis results using the percentile method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        /// <param name="recordParameterSets">Optional. Determines whether to record parameter sets. Default = true.</param>
        public UncertaintyAnalysisResults Estimate(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[]? distributions = null, bool recordParameterSets = true)
        {
            var results = new UncertaintyAnalysisResults();
            results.ParentDistribution = (UnivariateDistributionBase)Distribution;

            // get mode curve
            results.ModeCurve = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                results.ModeCurve[i] = Distribution.InverseCDF(probabilities[i]);

            // get bootstrapped list of distributions
            var bootDistributions = distributions != null ? distributions : Distributions();

            // get parameter sets
            if (recordParameterSets == true)
                results.ParameterSets = ParameterSets(bootDistributions);

            // get confidence intervals
            results.ConfidenceIntervals = PercentileQuantileCI(probabilities, alpha, bootDistributions);
     
            // create list of quantiles
            var minMax = ComputeMinMaxQuantiles(0.001, 1 - 1E-9, bootDistributions);
            double shift = 0;
            if (minMax[0] <= 0) shift = Math.Abs(minMax[0]) + 1d;
            double min = minMax[0] + shift;
            double max = minMax[1] + shift;
            double logMin = Math.Log10(min);
            int order = (int)Math.Floor(Math.Log10(max) - logMin);
            int bins = Math.Max(200, Math.Min(1000, 100 * order));
            double delta = (Math.Log10(max) - logMin) / (bins - 1);

            // Each ordinate is computed from the origin, not accumulated from its predecessor, so
            // rounding does not compound across the ladder.
            var quantiles = new List<double>(bins);
            for (int i = 0; i < bins; i++)
            {
                quantiles.Add(Math.Pow(10, logMin + i * delta) - shift);
            }

            // get mean curve
            results.MeanCurve = ExpectedProbabilities(quantiles, probabilities, bootDistributions);

            return results;
        }

        /// <summary>
        /// Bootstrap the expected non-exceedance probabilities given the input quantile values. Returns the x-values interpolated from the list of desired non-exceedance probabilities.
        /// </summary>
        /// <param name="quantiles">List quantile values.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[] ExpectedProbabilities(IList<double> quantiles, IList<double> probabilities, IUnivariateDistribution[]? distributions = null)
        {
            var quants = quantiles.ToArray();
            var probs = probabilities.ToArray();
            Array.Sort(quants);
            var bootDistributions = distributions != null ? distributions : Distributions();
            var expected = MeanCDFs(quants, bootDistributions);

            double minY = double.MaxValue;
            double maxY = double.MinValue;
            var yVals = new List<double>();
            var xVals = new List<double>();
            yVals.Add(quantiles[0]);
            xVals.Add(expected[0]);
            for (int i = 1; i < quantiles.Count; i++)
            {
                if (expected[i] > xVals.Last())
                {
                    minY = Math.Min(minY, quantiles[i]);
                    maxY = Math.Max(maxY, quantiles[i]);
                    yVals.Add(quantiles[i]);
                    xVals.Add(expected[i]);
                }
            }
            bool useLogTransform = false;
            if (minY > 0 && (Math.Log10(maxY) - Math.Log10(minY)) > 1)
                useLogTransform = true;

            Linear linint = new Linear(xVals, yVals) { XTransform = Transform.NormalZ, YTransform = useLogTransform ?  Transform.Logarithmic : Transform.None };
            return linint.Interpolate(probs);
        }

        /// <summary>
        /// The mean CDF across the bootstrapped distributions at each quantile.
        /// </summary>
        /// <param name="quantiles">The quantile values to evaluate, ascending.</param>
        /// <param name="distributions">The bootstrapped distributions; null entries are failed fits.</param>
        /// <returns>The expected non-exceedance probability at each quantile.</returns>
        /// <remarks>
        /// Replications are split into <see cref="ReductionChunks"/> chunks, each summed
        /// sequentially and merged in chunk order, so the result does not depend on the thread
        /// count. Failed fits are excluded from both the sum and the divisor.
        /// </remarks>
        private static double[] MeanCDFs(double[] quantiles, IUnivariateDistribution[] distributions)
        {
            int replications = distributions.Length;
            int quantileCount = quantiles.Length;
            var expected = new double[quantileCount];
            if (replications == 0 || quantileCount == 0) return expected;

            int chunks = Math.Min(ReductionChunks, replications);
            var chunkSums = new double[chunks][];
            var chunkValid = new int[chunks];
            for (int c = 0; c < chunks; c++) chunkSums[c] = new double[quantileCount];

            Parallel.For(0, chunks, c =>
            {
                var accumulator = chunkSums[c];
                int start = (int)((long)c * replications / chunks);
                int end = (int)((long)(c + 1) * replications / chunks);
                int valid = 0;
                for (int j = start; j < end; j++)
                {
                    var distribution = distributions[j];
                    if (distribution == null) continue;
                    valid++;
                    for (int i = 0; i < quantileCount; i++)
                    {
                        accumulator[i] += distribution.CDF(quantiles[i]);
                    }
                }
                chunkValid[c] = valid;
            });

            int validCount = 0;
            for (int c = 0; c < chunks; c++) validCount += chunkValid[c];
            if (validCount == 0) return expected;

            for (int i = 0; i < quantileCount; i++)
            {
                double total = 0d;
                for (int c = 0; c < chunks; c++) total += chunkSums[c][i];
                expected[i] = total / validCount;
            }
            return expected;
        }

        /// <summary>
        /// Bootstrap the expected non-exceedance probabilities given the input quantile values.
        /// </summary>
        /// <param name="quantiles">List quantile values.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[] ExpectedProbabilities(IList<double> quantiles, IUnivariateDistribution[]? distributions = null)
        {
            var quants = quantiles.ToArray();
            Array.Sort(quants);
            var bootDistributions = distributions != null ? distributions : Distributions();
            return MeanCDFs(quants, bootDistributions);
        }

        /// <summary>
        /// Returns the min and max quantiles from a bootstrap analysis.
        /// </summary>
        /// <param name="minProbability">The minimum probability to compute quantiles.</param>
        /// <param name="maxProbability">The maximum probability to compute quantiles.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[] ComputeMinMaxQuantiles(double minProbability, double maxProbability, IUnivariateDistribution[] distributions)
        {
            // Thread-local extremes merged once per partition rather than a lock per distribution.
            var output = new double[] { double.MaxValue, double.MinValue };
            object lockObject = new object();
            int count = distributions.Length;
            Parallel.For(0, count, () => (Min: double.MaxValue, Max: double.MinValue), (j, loop, local) =>
            {
                var distribution = distributions[j];
                if (distribution == null) return local;
                double minX = distribution.InverseCDF(minProbability);
                double maxX = distribution.InverseCDF(maxProbability);
                return (minX < local.Min ? minX : local.Min, maxX > local.Max ? maxX : local.Max);
            },
            local =>
            {
                lock (lockObject)
                {
                    if (local.Min < output[0]) output[0] = local.Min;
                    if (local.Max > output[1]) output[1] = local.Max;
                }
            });
            return output;
        }

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the percentile method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[,] PercentileQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[]? distributions = null)
        {
            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];
            var bootDistributions = distributions != null ? distributions : Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var XValues = new double[bootDistributions.Length];
                Parallel.For(0, bootDistributions.Length, idx => { XValues[idx] = bootDistributions[idx] != null ? bootDistributions[idx].InverseCDF(probabilities[i]) : double.NaN; });

                // Filter valid values and sort
                int validCount = 0;
                for (int k = 0; k < XValues.Length; k++)
                {
                    if (!double.IsNaN(XValues[k])) validCount++;
                }
                var validValues = new double[validCount];
                int writeIdx = 0;
                for (int k = 0; k < XValues.Length; k++)
                {
                    if (!double.IsNaN(XValues[k]))
                        validValues[writeIdx++] = XValues[k];
                }
                Array.Sort(validValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                    Output[i, j] = Statistics.Percentile(validValues, CIs[j], true);
            }
            return Output;
        }

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the bias-corrected percentile method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[,] BiasCorrectedQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[]? distributions = null)
        {
            // Create list of original X values given probability values
            var populationXValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationXValues[i] = Distribution.InverseCDF(probabilities[i]);

            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];
            var bootDistributions = distributions != null ? distributions : Distributions();
            int replications = bootDistributions.Length;
            for (int i = 0; i < probabilities.Count; i++)
            {
                var XValues = new double[replications];
                Parallel.For(0, replications, idx =>
                {
                    XValues[idx] = bootDistributions[idx] != null ? bootDistributions[idx].InverseCDF(probabilities[i]) : double.NaN;
                });

                // Counted sequentially so the proportion does not depend on the thread count.
                double P0 = 0d; // proportions of values less than population
                for (int idx = 0; idx < replications; idx++)
                {
                    if (!double.IsNaN(XValues[idx]) && XValues[idx] <= populationXValues[i]) P0 += 1d;
                }

                // get proportion
                P0 = P0 / (replications + 1);

                // Filter valid values and sort
                int validCount = 0;
                for (int k = 0; k < XValues.Length; k++)
                {
                    if (!double.IsNaN(XValues[k])) validCount++;
                }
                var validValues = new double[validCount];
                int writeIdx = 0;
                for (int k = 0; k < XValues.Length; k++)
                {
                    if (!double.IsNaN(XValues[k]))
                        validValues[writeIdx++] = XValues[k];
                }
                Array.Sort(validValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                {
                    double Z0 = Normal.StandardZ(P0);
                    double Z = Normal.StandardZ(CIs[j]);
                    double BC = Normal.StandardCDF(2d * Z0 + Z);
                    Output[i, j] = Statistics.Percentile(validValues, BC, true);
                }
            }
            return Output;
        }

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the Normal, or standard method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[,] NormalQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[]? distributions = null)
        {

            // Create list of original X values given probability values
            // Use a cube-root transform to make results transformation invariant
            var populationXValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationXValues[i] = Math.Pow(Distribution.InverseCDF(probabilities[i]), 1d / 3d);

            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];
            var bootDistributions = distributions != null ? distributions : Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var XValues = new double[bootDistributions.Length];
                Parallel.For(0, bootDistributions.Length, idx => { XValues[idx] = bootDistributions[idx] != null ? Math.Pow(bootDistributions[idx].InverseCDF(probabilities[i]), 1d / 3d) : double.NaN; });

                // Filter valid values
                int validCount = 0;
                for (int k = 0; k < XValues.Length; k++)
                {
                    if (!double.IsNaN(XValues[k])) validCount++;
                }
                var validValues = new double[validCount];
                int writeIdx = 0;
                for (int k = 0; k < XValues.Length; k++)
                {
                    if (!double.IsNaN(XValues[k]))
                        validValues[writeIdx++] = XValues[k];
                }

                // Get Standard error
                double SE = Statistics.StandardDeviation(validValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                {
                    double Z = Normal.StandardZ(CIs[j]);
                    Output[i, j] = Math.Pow(populationXValues[i] + SE * Z, 3d);
                }
            }
            return Output;
        }

        #region Bias-Corrected and Accelerated

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the bias-corrected and accelerated (BCa) percentile method.
        /// </summary>
        /// <param name="sampleData">Sample of data.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,] BCaQuantileCI(IList<double> sampleData, IList<double> probabilities, double alpha = 0.1)
        {
            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];

            // Estimate distribution
            SampleSize = sampleData.Count;
            ((IEstimation)Distribution).Estimate(sampleData, EstimationMethod);

            // Create list of original X values given probability values
            var populationXValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationXValues[i] = Distribution.InverseCDF(probabilities[i]);

            // Get acceleration constants
            var a = AccelerationConstants(sampleData, probabilities, populationXValues);

            // Get bootstrapped distributions
            var bootDistributions = Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var XValues = new double[Replications];
                Parallel.For(0, Replications, idx =>
                {
                    XValues[idx] = bootDistributions[idx] != null ? bootDistributions[idx].InverseCDF(probabilities[i]) : double.NaN;
                });

                // Counted sequentially so the proportion does not depend on the thread count.
                double P0 = 0d; // proportions of values less than population
                for (int idx = 0; idx < Replications; idx++)
                {
                    if (!double.IsNaN(XValues[idx]) && XValues[idx] <= populationXValues[i]) P0 += 1d;
                }

                // get proportion
                P0 = (P0 + 1) / (Replications + 1);

                // Filter valid values and sort
                int validCount = 0;
                for (int k = 0; k < XValues.Length; k++)
                {
                    if (!double.IsNaN(XValues[k])) validCount++;
                }
                var validValues = new double[validCount];
                int writeIdx = 0;
                for (int k = 0; k < XValues.Length; k++)
                {
                    if (!double.IsNaN(XValues[k]))
                        validValues[writeIdx++] = XValues[k];
                }
                Array.Sort(validValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                {
                    double Z0 = Normal.StandardZ(P0);
                    double Z = Normal.StandardZ(CIs[j]);
                    double num = Z0 + Z;
                    double den = 1 - a[i] * (Z0 + Z);
                    double BC = Normal.StandardCDF(Z0 + num / den);
                    Output[i, j] = Statistics.Percentile(validValues, BC, true);
                }
            }
            return Output;
        }

        /// <summary>
        /// Estimates the acceleration constants for each probability.
        /// </summary>
        /// <param name="sampleData">Sample of data.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="thetaHats">The list of best-estimate quantiles.</param>
        /// <remarks>
        /// Chunked so the moment sums merge in a fixed order, independent of the thread count.
        /// Each chunk refills one leave-one-out buffer rather than copying the sample per point.
        /// </remarks>
        private double[] AccelerationConstants(IList<double> sampleData, IList<double> probabilities, IList<double> thetaHats)
        {
            var N = sampleData.Count;
            int probabilityCount = probabilities.Count;
            var a = new double[probabilityCount];
            if (N == 0) return a;

            int chunks = Math.Min(ReductionChunks, N);
            var chunkI2 = new double[chunks][];
            var chunkI3 = new double[chunks][];
            for (int c = 0; c < chunks; c++)
            {
                chunkI2[c] = new double[probabilityCount];
                chunkI3[c] = new double[probabilityCount];
            }

            Parallel.For(0, chunks, c =>
            {
                var i2 = chunkI2[c];
                var i3 = chunkI3[c];
                var jackSample = new double[N - 1];
                int start = (int)((long)c * N / chunks);
                int end = (int)((long)(c + 1) * N / chunks);
                for (int idx = start; idx < end; idx++)
                {
                    for (int k = 0; k < idx; k++) jackSample[k] = sampleData[k];
                    for (int k = idx + 1; k < N; k++) jackSample[k - 1] = sampleData[k];

                    // Cloned per point: a failed Estimate can leave the instance partially set.
                    var newDistribution = ((UnivariateDistributionBase)Distribution).Clone();
                    try
                    {
                        ((IEstimation)newDistribution).Estimate(jackSample, EstimationMethod);
                        for (int i = 0; i < probabilityCount; i++)
                        {
                            double thetaJack = newDistribution.InverseCDF(probabilities[i]);
                            i2[i] += Math.Pow(thetaHats[i] - thetaJack, 2);
                            i3[i] += Math.Pow(thetaHats[i] - thetaJack, 3);
                        }
                    }
                    catch (Exception)
                    {
                        // MLE and certain L-moments methods can fail to find a solution
                    };
                }
            });

            // Get acceleration constant
            for (int i = 0; i < probabilityCount; i++)
            {
                double I2 = 0d, I3 = 0d;
                for (int c = 0; c < chunks; c++)
                {
                    I2 += chunkI2[c][i];
                    I3 += chunkI3[c][i];
                }
                a[i] = I3 / (Math.Pow(I2, 1.5) * 6);
            }

            return a;
        }

        #endregion

        #region Bootstrap-t (aka Student-t Bootstrap)

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the Bootstrap-t method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,] BootstrapTQuantileCI(IList<double> probabilities, double alpha = 0.1)
        {
            // Create list of original X values given probability values
            // Use a cube-root transform to make results transformation invariant
            var populationXValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationXValues[i] = Math.Pow(Distribution.InverseCDF(probabilities[i]), 1d / 3d);

            var xValues = new double[Replications, probabilities.Count];
            var studentT = new double[Replications, probabilities.Count];
            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];

            // First create list of bootstrap distributions, 
            // and estimate standard error for each quantiles          
            var bootDistributions = new IUnivariateDistribution[Replications];
            var r = new MersenneTwister(PRNGSeed);
            var seeds = r.NextIntegers(Replications);
            Parallel.For(0, Replications, i =>
            {
                try
                {
                    var newDistribution = ((UnivariateDistributionBase)Distribution).Clone();
                    var sample = newDistribution.GenerateRandomValues(SampleSize, seeds[i]);
                    ((IEstimation)newDistribution).Estimate(sample, EstimationMethod);
                    bootDistributions[i] = newDistribution;

                    // Record inner boot thetas
                    var bootXValues = new double[probabilities.Count];
                    for (int j = 0; j < probabilities.Count; j++)
                        bootXValues[j] = Math.Pow(bootDistributions[i].InverseCDF(probabilities[j]), 1d / 3d);

                    // Now estimate the standard error at each quantile using the jackknife method
                    //var bootSE = StandardError(sample, probabilities, bootXValues);
                    var bootSE = BootstrapStandardError(newDistribution, probabilities, 300, seeds[i]);
                    for (int j = 0; j < probabilities.Count; j++)
                    {
                        xValues[i, j] = bootXValues[j];
                        studentT[i, j] = (populationXValues[j] - bootXValues[j]) / bootSE[j];
                    }

                }
                catch (Exception)
                {
                    // MLE and certain L-moments methods can fail to find a solution
                    // On fail, set to null
                    bootDistributions[i] = null!;
                    for (int j = 0; j < probabilities.Count; j++)
                    {
                        xValues[i, j] = double.NaN;
                        studentT[i, j] = double.NaN;
                    }
                };

            });


            for (int i = 0; i < probabilities.Count; i++)
            {
                var rawX = xValues.GetColumn(i);
                var rawT = studentT.GetColumn(i);
                int validCount = 0;
                for (int k = 0; k < rawX.Length; k++)
                {
                    if (!double.IsNaN(rawX[k])) validCount++;
                }
                var XValues = new double[validCount];
                var TValues = new double[validCount];
                int writeIdx = 0;
                for (int k = 0; k < rawX.Length; k++)
                {
                    if (!double.IsNaN(rawX[k]))
                    {
                        XValues[writeIdx] = rawX[k];
                        TValues[writeIdx] = rawT[k];
                        writeIdx++;
                    }
                }

                // Get Standard error
                double SE = Statistics.StandardDeviation(XValues);
                Array.Sort(TValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                {
                    double T = Statistics.Percentile(TValues, CIs[j], true);
                    Output[i, j] = Math.Pow(populationXValues[i] + SE * T, 3d);
                }
            }

            return Output;
        }

        /// <summary>
        /// Estimates the standard error for each probability using the parametric bootstrap.
        /// </summary>
        ///<param name="parentDist">The parent distribution.</param>
        ///<param name="probabilities">The list of probabilities where the standard error is calculated.</param>
        ///<param name="replications">The number of bootstrap replications. Default = 300.</param>
        ///<param name="seed">The PRNG seed. Default = 12345.</param>
        private double[] BootstrapStandardError(UnivariateDistributionBase parentDist, IList<double> probabilities, int replications = 300, int seed = 12345)
        {
            int B = replications;
            var r = new MersenneTwister(seed);
            var seeds = r.NextIntegers(replications);
            var xValues = new double[B, probabilities.Count];
            var se = new double[probabilities.Count];
            Parallel.For(0, replications, i =>
            {
                try
                {
                    var bootDist = parentDist.Clone();
                    var sample = bootDist.GenerateRandomValues(SampleSize, seeds[i]);
                    ((IEstimation)bootDist).Estimate(sample, EstimationMethod);

                    // Record inner boot thetas
                    for (int j = 0; j < probabilities.Count; j++)
                        xValues[i, j] = Math.Pow(bootDist.InverseCDF(probabilities[j]), 1d / 3d);

                }
                catch (Exception)
                {
                    // MLE and certain L-moments methods can fail to find a solution
                    // On fail, set to null

                };

            });

            // Get standard error
            for (int i = 0; i < probabilities.Count; i++)
                se[i] = Statistics.StandardDeviation(xValues.GetColumn(i));
            return se;
        }

        /// <summary>
        /// Estimates the standard error for each probability.
        /// </summary>
        /// <param name="sampleData">Sample of data.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="thetaHats">The list of best-estimate quantiles.</param>
        /// <remarks>
        /// Chunked as <see cref="AccelerationConstants"/>.
        /// </remarks>
        private double[] StandardError(IList<double> sampleData, IList<double> probabilities, IList<double> thetaHats)
        {
            var N = sampleData.Count;
            int probabilityCount = probabilities.Count;
            var se = new double[probabilityCount];
            if (N == 0) return se;

            int chunks = Math.Min(ReductionChunks, N);
            var chunkI2 = new double[chunks][];
            for (int c = 0; c < chunks; c++) chunkI2[c] = new double[probabilityCount];

            // Perform Jackknife
            Parallel.For(0, chunks, c =>
            {
                var i2 = chunkI2[c];
                var jackSample = new double[N - 1];
                int start = (int)((long)c * N / chunks);
                int end = (int)((long)(c + 1) * N / chunks);
                for (int idx = start; idx < end; idx++)
                {
                    for (int k = 0; k < idx; k++) jackSample[k] = sampleData[k];
                    for (int k = idx + 1; k < N; k++) jackSample[k - 1] = sampleData[k];

                    var newDistribution = ((UnivariateDistributionBase)Distribution).Clone();
                    try
                    {
                        ((IEstimation)newDistribution).Estimate(jackSample, EstimationMethod);
                        for (int i = 0; i < probabilityCount; i++)
                        {
                            double thetaJack = Math.Pow(newDistribution.InverseCDF(probabilities[i]), 1d / 3d);
                            i2[i] += Math.Pow(thetaHats[i] - thetaJack, 2);
                        }
                    }
                    catch (Exception)
                    {
                        // MLE and certain L-moments methods can fail to find a solution
                    };
                }
            });

            // Get standard error
            for (int i = 0; i < probabilityCount; i++)
            {
                double I2 = 0d;
                for (int c = 0; c < chunks; c++) I2 += chunkI2[c][i];
                se[i] = Math.Sqrt((N - 1) / (double)N * I2);
            }

            return se;
        }

        #endregion

        #endregion

    }
}