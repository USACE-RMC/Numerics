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
        /// Generates fitted bootstrap distributions.
        /// </summary>
        /// <returns>The fitted distributions; isolated failed replications are represented by null entries.</returns>
        /// <exception cref="AggregateException">Thrown when every replication fails after all retries.</exception>
        public IUnivariateDistribution[] Distributions()
        {
            var bootDistributions = new IUnivariateDistribution[Replications];
            var failuresByReplication = new Exception?[Replications];
            var random = new MersenneTwister(PRNGSeed);
            var seeds = random.NextIntegers(Replications);
            int failures = 0;

            Parallel.For(0, Replications, index =>
            {
                Exception? lastFailure = null;
                for (int attempt = 0; attempt < _retries; attempt++)
                {
                    try
                    {
                        bootDistributions[index] = Distribution.Bootstrap(EstimationMethod, SampleSize, seeds[index] + 10 * attempt);
                        lastFailure = null;
                        break;
                    }
                    catch (Exception exception)
                    {
                        lastFailure = exception;
                    }
                }

                if (lastFailure != null || bootDistributions[index] == null)
                {
                    failuresByReplication[index] = lastFailure
                        ?? new InvalidOperationException("The bootstrap fit returned no distribution.");
                    bootDistributions[index] = null!;
                    Interlocked.Increment(ref failures);
                }
            });

            FailedReplications = failures;
            if (failures == Replications)
                throw new AggregateException("Every bootstrap distribution fit failed.", failuresByReplication.Where(exception => exception != null).Cast<Exception>());
            return bootDistributions;
        }

        /// <summary>
        /// Creates fitted distributions from parameter sets.
        /// </summary>
        /// <param name="parameterSets">The parameter sets.</param>
        /// <returns>The distributions; isolated invalid sets are represented by null entries.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="parameterSets"/> is null.</exception>
        /// <exception cref="AggregateException">Thrown when every parameter set is invalid.</exception>
        public IUnivariateDistribution[] Distributions(ParameterSet[] parameterSets)
        {
            if (parameterSets == null) throw new ArgumentNullException(nameof(parameterSets));
            var bootDistributions = new IUnivariateDistribution[parameterSets.Length];
            if (parameterSets.Length == 0) return bootDistributions;

            var failuresByReplication = new Exception?[parameterSets.Length];
            int failures = 0;
            Parallel.For(0, parameterSets.Length, index =>
            {
                try
                {
                    var distribution = ((UnivariateDistributionBase)Distribution).Clone();
                    distribution.ValidateParameters(parameterSets[index].Values, true);
                    distribution.SetParameters(parameterSets[index].Values);
                    if (!distribution.ParametersValid)
                        throw new ArgumentException("The parameter set does not define a valid distribution.", nameof(parameterSets));
                    bootDistributions[index] = distribution;
                }
                catch (Exception exception)
                {
                    failuresByReplication[index] = exception;
                    bootDistributions[index] = null!;
                    Interlocked.Increment(ref failures);
                }
            });

            FailedReplications = failures;
            if (failures == parameterSets.Length)
                throw new AggregateException("Every bootstrap parameter set was invalid.", failuresByReplication.Where(exception => exception != null).Cast<Exception>());
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
        /// Interpolates quantiles at requested probabilities from the mean bootstrap CDF.
        /// </summary>
        /// <param name="quantiles">Quantile ordinates; ordering is not required.</param>
        /// <param name="probabilities">The probabilities to interpolate.</param>
        /// <param name="distributions">Optional precomputed bootstrap distributions.</param>
        /// <returns>The interpolated quantiles.</returns>
        public double[] ExpectedProbabilities(IList<double> quantiles, IList<double> probabilities, IUnivariateDistribution[]? distributions = null)
        {
            if (quantiles == null) throw new ArgumentNullException(nameof(quantiles));
            if (probabilities == null) throw new ArgumentNullException(nameof(probabilities));
            if (quantiles.Count < 2) throw new ArgumentException("At least two quantiles are required.", nameof(quantiles));

            var sortedQuantiles = quantiles.ToArray();
            Array.Sort(sortedQuantiles);
            var targetProbabilities = probabilities.ToArray();
            var bootDistributions = distributions ?? Distributions();
            var expected = MeanCDFs(sortedQuantiles, bootDistributions);

            var cdfValues = new List<double> { expected[0] };
            var ordinateValues = new List<double> { sortedQuantiles[0] };
            double minimumOrdinate = sortedQuantiles[0];
            double maximumOrdinate = sortedQuantiles[0];
            for (int i = 1; i < sortedQuantiles.Length; i++)
            {
                if (expected[i] > cdfValues[cdfValues.Count - 1])
                {
                    minimumOrdinate = Math.Min(minimumOrdinate, sortedQuantiles[i]);
                    maximumOrdinate = Math.Max(maximumOrdinate, sortedQuantiles[i]);
                    ordinateValues.Add(sortedQuantiles[i]);
                    cdfValues.Add(expected[i]);
                }
            }
            if (cdfValues.Count < 2)
                throw new InvalidOperationException("The mean bootstrap CDF does not contain two distinct probabilities.");

            bool useLogTransform = minimumOrdinate > 0d
                && Math.Log10(maximumOrdinate) - Math.Log10(minimumOrdinate) > 1d;
            var interpolation = new Linear(cdfValues, ordinateValues)
            {
                XTransform = Transform.NormalZ,
                YTransform = useLogTransform ? Transform.Logarithmic : Transform.None
            };
            return interpolation.Interpolate(targetProbabilities);
        }

        /// <summary>
        /// Computes the mean CDF across successful bootstrap fits at each quantile.
        /// </summary>
        /// <param name="quantiles">The quantiles to evaluate.</param>
        /// <param name="distributions">The bootstrap distributions; null entries represent failed fits.</param>
        /// <returns>The mean CDF values.</returns>
        private static double[] MeanCDFs(double[] quantiles, IUnivariateDistribution[] distributions)
        {
            int replications = distributions.Length;
            int quantileCount = quantiles.Length;
            var expected = new double[quantileCount];
            if (quantileCount == 0) return expected;
            if (replications == 0)
                throw new InvalidOperationException("No bootstrap distributions were supplied.");

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
                        accumulator[i] += distribution.CDF(quantiles[i]);
                }
                chunkValid[c] = valid;
            });

            int validCount = 0;
            for (int c = 0; c < chunks; c++) validCount += chunkValid[c];
            if (validCount == 0)
                throw new InvalidOperationException("Every bootstrap distribution fit failed.");

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
            if (output[0] == double.MaxValue || output[1] == double.MinValue)
                throw new InvalidOperationException("Every bootstrap distribution fit failed.");
            return output;
        }

        /// <summary>
        /// Computes percentile bootstrap confidence intervals for quantiles.
        /// </summary>
        /// <param name="probabilities">The non-exceedance probabilities.</param>
        /// <param name="alpha">The excluded two-sided probability.</param>
        /// <param name="distributions">Optional precomputed bootstrap distributions.</param>
        /// <returns>The lower and upper confidence limits for each probability.</returns>
        public double[,] PercentileQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[]? distributions = null)
        {
            var confidenceProbabilities = new[] { alpha / 2d, 1d - alpha / 2d };
            var output = new double[probabilities.Count, 2];
            var bootDistributions = distributions ?? Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var values = new double[bootDistributions.Length];
                Parallel.For(0, bootDistributions.Length, index =>
                {
                    values[index] = bootDistributions[index] != null
                        ? bootDistributions[index].InverseCDF(probabilities[i])
                        : double.NaN;
                });

                var successfulValues = FiniteValuesOrThrow(values, "percentile confidence intervals");
                Array.Sort(successfulValues);
                for (int j = 0; j < confidenceProbabilities.Length; j++)
                    output[i, j] = Statistics.Percentile(successfulValues, confidenceProbabilities[j], true);
            }
            return output;
        }
        /// <summary>
        /// Computes bias-corrected percentile confidence intervals for quantiles.
        /// </summary>
        /// <param name="probabilities">The non-exceedance probabilities.</param>
        /// <param name="alpha">The excluded two-sided probability.</param>
        /// <param name="distributions">Optional precomputed bootstrap distributions.</param>
        /// <returns>The lower and upper confidence limits for each probability.</returns>
        public double[,] BiasCorrectedQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[]? distributions = null)
        {
            var populationValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationValues[i] = Distribution.InverseCDF(probabilities[i]);

            var confidenceProbabilities = new[] { alpha / 2d, 1d - alpha / 2d };
            var output = new double[probabilities.Count, 2];
            var bootDistributions = distributions ?? Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var values = new double[bootDistributions.Length];
                Parallel.For(0, bootDistributions.Length, index =>
                {
                    values[index] = bootDistributions[index] != null
                        ? bootDistributions[index].InverseCDF(probabilities[i])
                        : double.NaN;
                });

                var successfulValues = FiniteValuesOrThrow(values, "bias-corrected confidence intervals");
                int lessOrEqual = 0;
                for (int index = 0; index < successfulValues.Length; index++)
                    if (successfulValues[index] <= populationValues[i]) lessOrEqual++;

                double proportion = (lessOrEqual + 1d) / (successfulValues.Length + 1d);
                Array.Sort(successfulValues);
                double bias = Normal.StandardZ(proportion);
                for (int j = 0; j < confidenceProbabilities.Length; j++)
                {
                    double adjusted = Normal.StandardCDF(2d * bias + Normal.StandardZ(confidenceProbabilities[j]));
                    output[i, j] = Statistics.Percentile(successfulValues, adjusted, true);
                }
            }
            return output;
        }
        /// <summary>
        /// Computes normal-approximation bootstrap confidence intervals for quantiles.
        /// </summary>
        /// <param name="probabilities">The non-exceedance probabilities.</param>
        /// <param name="alpha">The excluded two-sided probability.</param>
        /// <param name="distributions">Optional precomputed bootstrap distributions.</param>
        /// <returns>The lower and upper confidence limits for each probability.</returns>
        public double[,] NormalQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[]? distributions = null)
        {
            var transformedPopulationValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                transformedPopulationValues[i] = CubeRoot(Distribution.InverseCDF(probabilities[i]));

            var confidenceProbabilities = new[] { alpha / 2d, 1d - alpha / 2d };
            var output = new double[probabilities.Count, 2];
            var bootDistributions = distributions ?? Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var transformedValues = new double[bootDistributions.Length];
                Parallel.For(0, bootDistributions.Length, index =>
                {
                    transformedValues[index] = bootDistributions[index] != null
                        ? CubeRoot(bootDistributions[index].InverseCDF(probabilities[i]))
                        : double.NaN;
                });

                var successfulValues = FiniteValuesOrThrow(transformedValues, "normal confidence intervals", 2);
                double standardError = Statistics.StandardDeviation(successfulValues);
                for (int j = 0; j < confidenceProbabilities.Length; j++)
                {
                    double transformedLimit = transformedPopulationValues[i]
                        + standardError * Normal.StandardZ(confidenceProbabilities[j]);
                    output[i, j] = transformedLimit * transformedLimit * transformedLimit;
                }
            }
            return output;
        }
        #region Bias-Corrected and Accelerated

        /// <summary>
        /// Computes bias-corrected and accelerated percentile confidence intervals.
        /// </summary>
        /// <param name="sampleData">The sample used to estimate the parent distribution.</param>
        /// <param name="probabilities">The non-exceedance probabilities.</param>
        /// <param name="alpha">The excluded two-sided probability.</param>
        /// <returns>The lower and upper confidence limits for each probability.</returns>
        public double[,] BCaQuantileCI(IList<double> sampleData, IList<double> probabilities, double alpha = 0.1)
        {
            var confidenceProbabilities = new[] { alpha / 2d, 1d - alpha / 2d };
            var output = new double[probabilities.Count, 2];

            SampleSize = sampleData.Count;
            ((IEstimation)Distribution).Estimate(sampleData, EstimationMethod);
            var populationValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationValues[i] = Distribution.InverseCDF(probabilities[i]);

            var acceleration = AccelerationConstants(sampleData, probabilities, populationValues);
            var bootDistributions = Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var values = new double[bootDistributions.Length];
                Parallel.For(0, bootDistributions.Length, index =>
                {
                    values[index] = bootDistributions[index] != null
                        ? bootDistributions[index].InverseCDF(probabilities[i])
                        : double.NaN;
                });

                var successfulValues = FiniteValuesOrThrow(values, "BCa confidence intervals");
                int lessOrEqual = 0;
                for (int index = 0; index < successfulValues.Length; index++)
                    if (successfulValues[index] <= populationValues[i]) lessOrEqual++;

                double proportion = (lessOrEqual + 1d) / (successfulValues.Length + 1d);
                double bias = Normal.StandardZ(proportion);
                Array.Sort(successfulValues);
                for (int j = 0; j < confidenceProbabilities.Length; j++)
                {
                    double normalQuantile = Normal.StandardZ(confidenceProbabilities[j]);
                    double denominator = 1d - acceleration[i] * (bias + normalQuantile);
                    double adjusted = Normal.StandardCDF(bias + (bias + normalQuantile) / denominator);
                    output[i, j] = Statistics.Percentile(successfulValues, adjusted, true);
                }
            }
            return output;
        }
        /// <summary>
        /// Estimates acceleration constants from successful leave-one-out fits.
        /// </summary>
        /// <param name="sampleData">The observed sample.</param>
        /// <param name="probabilities">The non-exceedance probabilities.</param>
        /// <param name="thetaHats">The fitted population quantiles.</param>
        /// <returns>One acceleration constant per probability.</returns>
        /// <exception cref="AggregateException">Thrown when every leave-one-out fit fails.</exception>
        private double[] AccelerationConstants(IList<double> sampleData, IList<double> probabilities, IList<double> thetaHats)
        {
            int sampleCount = sampleData.Count;
            int probabilityCount = probabilities.Count;
            var acceleration = new double[probabilityCount];
            if (sampleCount == 0) return acceleration;

            int chunks = Math.Min(ReductionChunks, sampleCount);
            var chunkSecondMoments = new double[chunks][];
            var chunkThirdMoments = new double[chunks][];
            var chunkSuccesses = new int[chunks];
            var failures = new Exception?[sampleCount];
            for (int chunk = 0; chunk < chunks; chunk++)
            {
                chunkSecondMoments[chunk] = new double[probabilityCount];
                chunkThirdMoments[chunk] = new double[probabilityCount];
            }

            Parallel.For(0, chunks, chunk =>
            {
                var secondMoments = chunkSecondMoments[chunk];
                var thirdMoments = chunkThirdMoments[chunk];
                int start = (int)((long)chunk * sampleCount / chunks);
                int end = (int)((long)(chunk + 1) * sampleCount / chunks);
                int successes = 0;
                for (int index = start; index < end; index++)
                {
                    var jackknifeSample = new double[sampleCount - 1];
                    for (int k = 0; k < index; k++) jackknifeSample[k] = sampleData[k];
                    for (int k = index + 1; k < sampleCount; k++) jackknifeSample[k - 1] = sampleData[k];

                    var distribution = ((UnivariateDistributionBase)Distribution).Clone();
                    try
                    {
                        ((IEstimation)distribution).Estimate(jackknifeSample, EstimationMethod);
                        for (int i = 0; i < probabilityCount; i++)
                        {
                            double difference = thetaHats[i] - distribution.InverseCDF(probabilities[i]);
                            secondMoments[i] += difference * difference;
                            thirdMoments[i] += difference * difference * difference;
                        }
                        successes++;
                    }
                    catch (Exception exception)
                    {
                        failures[index] = exception;
                    }
                }
                chunkSuccesses[chunk] = successes;
            });

            int successfulFits = 0;
            for (int chunk = 0; chunk < chunks; chunk++) successfulFits += chunkSuccesses[chunk];
            if (successfulFits == 0)
                throw new AggregateException("Every jackknife acceleration fit failed.", failures.Where(exception => exception != null).Cast<Exception>());

            for (int i = 0; i < probabilityCount; i++)
            {
                double secondMoment = 0d;
                double thirdMoment = 0d;
                for (int chunk = 0; chunk < chunks; chunk++)
                {
                    secondMoment += chunkSecondMoments[chunk][i];
                    thirdMoment += chunkThirdMoments[chunk][i];
                }
                acceleration[i] = secondMoment > 0d && Tools.IsFinite(secondMoment)
                    ? thirdMoment / (6d * Math.Pow(secondMoment, 1.5d))
                    : 0d;
            }
            return acceleration;
        }
        #endregion

        #region Bootstrap-t (aka Student-t Bootstrap)

        /// <summary>
        /// Computes studentized bootstrap confidence intervals for quantiles.
        /// </summary>
        /// <param name="probabilities">The non-exceedance probabilities.</param>
        /// <param name="alpha">The excluded two-sided probability.</param>
        /// <returns>The lower and upper confidence limits for each probability.</returns>
        public double[,] BootstrapTQuantileCI(IList<double> probabilities, double alpha = 0.1)
        {
            var populationValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationValues[i] = CubeRoot(Distribution.InverseCDF(probabilities[i]));

            var transformedValues = new double[Replications, probabilities.Count];
            var studentizedValues = new double[Replications, probabilities.Count];
            var failures = new Exception?[Replications];
            var confidenceProbabilities = new[] { alpha / 2d, 1d - alpha / 2d };
            var output = new double[probabilities.Count, 2];
            var random = new MersenneTwister(PRNGSeed);
            var seeds = random.NextIntegers(Replications);
            int failedFits = 0;

            Parallel.For(0, Replications, index =>
            {
                try
                {
                    var distribution = ((UnivariateDistributionBase)Distribution).Clone();
                    var sample = distribution.GenerateRandomValues(SampleSize, seeds[index]);
                    ((IEstimation)distribution).Estimate(sample, EstimationMethod);

                    var bootstrapValues = new double[probabilities.Count];
                    for (int j = 0; j < probabilities.Count; j++)
                        bootstrapValues[j] = CubeRoot(distribution.InverseCDF(probabilities[j]));

                    var standardErrors = BootstrapStandardError(distribution, probabilities, 300, seeds[index]);
                    for (int j = 0; j < probabilities.Count; j++)
                    {
                        transformedValues[index, j] = bootstrapValues[j];
                        studentizedValues[index, j] = (populationValues[j] - bootstrapValues[j]) / standardErrors[j];
                    }
                }
                catch (Exception exception)
                {
                    failures[index] = exception;
                    Interlocked.Increment(ref failedFits);
                    for (int j = 0; j < probabilities.Count; j++)
                    {
                        transformedValues[index, j] = double.NaN;
                        studentizedValues[index, j] = double.NaN;
                    }
                }
            });

            if (failedFits == Replications)
                throw new AggregateException("Every studentized bootstrap fit failed.", failures.Where(exception => exception != null).Cast<Exception>());

            for (int i = 0; i < probabilities.Count; i++)
            {
                var rawValues = transformedValues.GetColumn(i);
                var rawStudentized = studentizedValues.GetColumn(i);
                var values = new List<double>();
                var studentized = new List<double>();
                for (int index = 0; index < rawValues.Length; index++)
                {
                    if (Tools.IsFinite(rawValues[index]) && Tools.IsFinite(rawStudentized[index]))
                    {
                        values.Add(rawValues[index]);
                        studentized.Add(rawStudentized[index]);
                    }
                }
                if (values.Count < 2)
                    throw new InvalidOperationException("Insufficient successful fits are available for studentized confidence intervals.");

                double standardError = Statistics.StandardDeviation(values);
                var sortedStudentized = studentized.ToArray();
                Array.Sort(sortedStudentized);
                for (int j = 0; j < confidenceProbabilities.Length; j++)
                {
                    double studentizedQuantile = Statistics.Percentile(sortedStudentized, confidenceProbabilities[j], true);
                    double transformedLimit = populationValues[i] + standardError * studentizedQuantile;
                    output[i, j] = transformedLimit * transformedLimit * transformedLimit;
                }
            }
            return output;
        }
        /// <summary>
        /// Estimates quantile standard errors using successful inner bootstrap fits.
        /// </summary>
        /// <param name="parentDist">The fitted parent distribution used to generate inner samples.</param>
        /// <param name="probabilities">The non-exceedance probabilities.</param>
        /// <param name="replications">The number of inner bootstrap fits.</param>
        /// <param name="seed">The deterministic random-number seed.</param>
        /// <returns>One standard error per probability.</returns>
        /// <exception cref="AggregateException">Thrown when every inner bootstrap fit fails.</exception>
        /// <exception cref="InvalidOperationException">Thrown when fewer than two finite fits are available for a requested probability.</exception>
        private double[] BootstrapStandardError(UnivariateDistributionBase parentDist, IList<double> probabilities, int replications = 300, int seed = 12345)
        {
            var random = new MersenneTwister(seed);
            var seeds = random.NextIntegers(replications);
            var values = new double[replications, probabilities.Count];
            var failures = new Exception?[replications];
            int failedFits = 0;

            Parallel.For(0, replications, index =>
            {
                try
                {
                    var distribution = parentDist.Clone();
                    var sample = distribution.GenerateRandomValues(SampleSize, seeds[index]);
                    ((IEstimation)distribution).Estimate(sample, EstimationMethod);
                    for (int j = 0; j < probabilities.Count; j++)
                        values[index, j] = CubeRoot(distribution.InverseCDF(probabilities[j]));
                }
                catch (Exception exception)
                {
                    failures[index] = exception;
                    Interlocked.Increment(ref failedFits);
                    for (int j = 0; j < probabilities.Count; j++) values[index, j] = double.NaN;
                }
            });

            if (failedFits == replications)
                throw new AggregateException("Every inner bootstrap fit failed.", failures.Where(exception => exception != null).Cast<Exception>());

            var standardErrors = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
            {
                var successfulValues = FiniteValuesOrThrow(values.GetColumn(i), "inner bootstrap standard errors", 2);
                standardErrors[i] = Statistics.StandardDeviation(successfulValues);
            }
            return standardErrors;
        }
        /// <summary>
        /// Estimates jackknife standard errors from successful leave-one-out fits.
        /// </summary>
        /// <param name="sampleData">The observed sample.</param>
        /// <param name="probabilities">The non-exceedance probabilities.</param>
        /// <param name="thetaHats">The cube-root-transformed fitted population quantiles.</param>
        /// <returns>One jackknife standard error per probability.</returns>
        /// <exception cref="AggregateException">Thrown when every leave-one-out fit fails.</exception>
        private double[] StandardError(IList<double> sampleData, IList<double> probabilities, IList<double> thetaHats)
        {
            int sampleCount = sampleData.Count;
            int probabilityCount = probabilities.Count;
            var standardErrors = new double[probabilityCount];
            if (sampleCount == 0) return standardErrors;

            int chunks = Math.Min(ReductionChunks, sampleCount);
            var chunkSecondMoments = new double[chunks][];
            var chunkSuccesses = new int[chunks];
            var failures = new Exception?[sampleCount];
            for (int chunk = 0; chunk < chunks; chunk++)
                chunkSecondMoments[chunk] = new double[probabilityCount];

            Parallel.For(0, chunks, chunk =>
            {
                var secondMoments = chunkSecondMoments[chunk];
                int start = (int)((long)chunk * sampleCount / chunks);
                int end = (int)((long)(chunk + 1) * sampleCount / chunks);
                int successes = 0;
                for (int index = start; index < end; index++)
                {
                    var jackknifeSample = new double[sampleCount - 1];
                    for (int k = 0; k < index; k++) jackknifeSample[k] = sampleData[k];
                    for (int k = index + 1; k < sampleCount; k++) jackknifeSample[k - 1] = sampleData[k];

                    var distribution = ((UnivariateDistributionBase)Distribution).Clone();
                    try
                    {
                        ((IEstimation)distribution).Estimate(jackknifeSample, EstimationMethod);
                        for (int i = 0; i < probabilityCount; i++)
                        {
                            double difference = thetaHats[i] - CubeRoot(distribution.InverseCDF(probabilities[i]));
                            secondMoments[i] += difference * difference;
                        }
                        successes++;
                    }
                    catch (Exception exception)
                    {
                        failures[index] = exception;
                    }
                }
                chunkSuccesses[chunk] = successes;
            });

            int successfulFits = 0;
            for (int chunk = 0; chunk < chunks; chunk++) successfulFits += chunkSuccesses[chunk];
            if (successfulFits == 0)
                throw new AggregateException("Every jackknife standard-error fit failed.", failures.Where(exception => exception != null).Cast<Exception>());

            for (int i = 0; i < probabilityCount; i++)
            {
                double secondMoment = 0d;
                for (int chunk = 0; chunk < chunks; chunk++) secondMoment += chunkSecondMoments[chunk][i];
                standardErrors[i] = successfulFits > 1
                    ? Math.Sqrt((successfulFits - 1d) / successfulFits * secondMoment)
                    : 0d;
            }
            return standardErrors;
        }
        /// <summary>
        /// Returns finite results from successful fits and enforces a minimum sample count.
        /// </summary>
        /// <param name="values">The fit results, with non-finite entries representing failed fits.</param>
        /// <param name="operation">The operation name included in failure messages.</param>
        /// <param name="minimumCount">The minimum number of finite results required.</param>
        /// <returns>A new array containing only finite fit results.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        /// <exception cref="InvalidOperationException">Thrown when fewer than <paramref name="minimumCount"/> fits succeeded.</exception>
        private static double[] FiniteValuesOrThrow(double[] values, string operation, int minimumCount = 1)
        {
            var successfulValues = values.Where(Tools.IsFinite).ToArray();
            if (successfulValues.Length < minimumCount)
                throw new InvalidOperationException("Insufficient successful fits are available for " + operation + ".");
            return successfulValues;
        }

        /// <summary>
        /// Computes the real cube root while preserving the sign of negative values.
        /// </summary>
        /// <param name="value">The value whose real cube root is required.</param>
        /// <returns>The sign-preserving real cube root.</returns>
        private static double CubeRoot(double value)
        {
            if (value == 0d) return value;
            return value < 0d ? -Math.Pow(-value, 1d / 3d) : Math.Pow(value, 1d / 3d);
        }
        #endregion

        #endregion

    }
}