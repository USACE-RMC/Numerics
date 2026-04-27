using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics.Optimization;
using System;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Numerics.Sampling
{

    /// <summary>
    /// A general-purpose bootstrap class for parametric or non-parametric bootstrap analysis.
    /// Supports all major confidence interval methods: Percentile, Bias-Corrected, BCa, Normal, and Bootstrap-t.
    /// </summary>
    /// <typeparam name="TData">The type of data being bootstrapped.</typeparam>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href="https://en.wikipedia.org/wiki/Bootstrapping_(statistics)" />
    /// </para>
    /// </remarks>
    public class Bootstrap<TData>
    {

        #region Construction

        /// <summary>
        /// Constructs a new bootstrap analysis.
        /// </summary>
        /// <param name="originalData">The original data.</param>
        /// <param name="originalParameters">The original fitted parameter set.</param>
        public Bootstrap(TData originalData, ParameterSet originalParameters)
        {
            _originalData = originalData;
            _originalParameters = originalParameters;
        }

        #endregion

        #region Members

        private TData _originalData;
        private ParameterSet _originalParameters;
        private ParameterSet[] _bootstrapParameterSets = null!;
        private double[,] _bootstrapStatistics = null!;
        private int _numStats;
        private int _numParams;
        private int _failedCount;
        private bool[] _validFlags = null!;
        private double[,]? _studentizedValues;
        private double[,]? _transformedStatistics;

        #endregion

        #region Properties

        /// <summary>
        /// Delegate function for resampling the original data given the current parameters and a random number generator.
        /// </summary>
        public Func<TData, ParameterSet, Random, TData>? ResampleFunction { get; set; }

        /// <summary>
        /// Delegate function for fitting a model to data and returning a parameter set.
        /// </summary>
        public Func<TData, ParameterSet>? FitFunction { get; set; }

        /// <summary>
        /// Delegate function for extracting statistics from a fitted parameter set.
        /// Returns an array of statistic values (e.g., quantiles at multiple probabilities).
        /// </summary>
        public Func<ParameterSet, double[]>? StatisticFunction { get; set; }

        /// <summary>
        /// Optional delegate for computing a leave-one-out jackknife sample.
        /// Takes the original data and the index of the observation to remove.
        /// Required for the BCa confidence interval method.
        /// </summary>
        public Func<TData, int, TData>? JackknifeFunction { get; set; }

        /// <summary>
        /// Optional delegate that returns the number of observations in the data.
        /// Required for the BCa confidence interval method.
        /// </summary>
        public Func<TData, int>? SampleSizeFunction { get; set; }

        /// <summary>
        /// Optional transform applied to statistic values before Normal and Bootstrap-t CI computation.
        /// Default is cube-root: x → x^(1/3). Set to null for no transform.
        /// </summary>
        public Func<double, double> Transform { get; set; } = x => Math.Pow(x, 1d / 3d);

        /// <summary>
        /// Optional inverse transform corresponding to Transform.
        /// Default is cube: x → x^3. Set to null for no transform.
        /// </summary>
        public Func<double, double> InverseTransform { get; set; } = x => Math.Pow(x, 3d);

        /// <summary>
        /// Number of bootstrap replicates. Default = 10,000.
        /// </summary>
        public int Replicates { get; set; } = 10000;

        /// <summary>
        /// Gets and sets the PRNG seed for reproducibility. Default = 12345.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// The maximum number of retries for a failed bootstrap replicate. Default = 20.
        /// </summary>
        public int MaxRetries { get; set; } = 20;

        /// <summary>
        /// The number of inner bootstrap replicates for Bootstrap-t standard error estimation. Default = 300.
        /// </summary>
        public int InnerReplicates { get; set; } = 300;

        /// <summary>
        /// Gets the bootstrapped model parameter sets.
        /// </summary>
        public ParameterSet[] BootstrapParameterSets => _bootstrapParameterSets;

        /// <summary>
        /// Gets the bootstrapped statistics as a 2D array [replicates, statistics].
        /// </summary>
        public double[,] BootstrapStatistics => _bootstrapStatistics;

        /// <summary>
        /// Gets the number of replicates that failed after all retries.
        /// </summary>
        public int FailedReplicates => _failedCount;

        #endregion

        #region Run Methods

        /// <summary>
        /// Runs the basic bootstrap procedure with error handling and retry logic.
        /// </summary>
        public void Run()
        {
            ValidateCoreDelegates();
            var resample = ResampleFunction!;
            var fit = FitFunction!;
            var statistic = StatisticFunction!;
            InitializeState();

            var prng = new MersenneTwister(PRNGSeed);
            var seeds = prng.NextIntegers(Replicates);

            Parallel.For(0, Replicates, idx =>
            {
                bool succeeded = false;
                for (int retry = 0; retry < MaxRetries; retry++)
                {
                    try
                    {
                        var rng = new MersenneTwister(seeds[idx] + 10 * retry);
                        var sample = resample(_originalData, _originalParameters, rng);
                        var fitResult = fit(sample);
                        var stat = statistic(fitResult);

                        // Validate: check for NaN in statistics
                        bool hasNaN = false;
                        for (int k = 0; k < _numStats; k++)
                        {
                            if (double.IsNaN(stat[k]) || double.IsInfinity(stat[k]))
                            {
                                hasNaN = true;
                                break;
                            }
                        }
                        if (hasNaN) continue;

                        _bootstrapParameterSets[idx] = fitResult;
                        for (int k = 0; k < _numStats; k++)
                            _bootstrapStatistics[idx, k] = stat[k];
                        _validFlags[idx] = true;
                        succeeded = true;
                    }
                    catch (Exception)
                    {
                        // retry
                    }
                    if (succeeded) break;
                }

                if (!succeeded)
                {
                    MarkFailed(idx);
                }
            });
        }

        /// <summary>
        /// Runs the double bootstrap procedure with bias correction.
        /// </summary>
        /// <param name="innerReplicates">Number of inner bootstrap replicates. Default = 300.</param>
        public void RunDoubleBootstrap(int innerReplicates = 300)
        {
            ValidateCoreDelegates();
            var resample = ResampleFunction!;
            var fit = FitFunction!;
            var statistic = StatisticFunction!;
            InitializeState();

            var prng = new MersenneTwister(PRNGSeed);
            var seeds = prng.NextIntegers(Replicates);

            Parallel.For(0, Replicates, idx =>
            {
                bool succeeded = false;
                for (int retry = 0; retry < MaxRetries; retry++)
                {
                    try
                    {
                        var rng = new MersenneTwister(seeds[idx] + 10 * retry);

                        // Outer bootstrap
                        var outerSample = resample(_originalData, _originalParameters, rng);
                        var outerFit = fit(outerSample);
                        var outerStat = statistic(outerFit);

                        // Inner bootstrap for bias estimation
                        int p = outerFit.Values.Length;
                        var parmsInnerSum = new double[p];
                        var statsInnerSum = new double[_numStats];
                        int validInner = 0;

                        for (int k = 0; k < innerReplicates; k++)
                        {
                            try
                            {
                                var innerSample = resample(outerSample, outerFit, rng);
                                var innerFit = fit(innerSample);
                                var innerStat = statistic(innerFit);

                                for (int i = 0; i < p; i++)
                                    parmsInnerSum[i] += innerFit.Values[i];
                                for (int i = 0; i < _numStats; i++)
                                    statsInnerSum[i] += innerStat[i];
                                validInner++;
                            }
                            catch (Exception)
                            {
                                // skip failed inner replicate
                            }
                        }

                        if (validInner == 0) continue;

                        // Bias-correct parameters
                        var biasCorrectedParms = new double[p];
                        for (int i = 0; i < p; i++)
                        {
                            double innerMean = parmsInnerSum[i] / validInner;
                            biasCorrectedParms[i] = outerFit.Values[i] - (innerMean - outerFit.Values[i]);
                        }

                        // Bias-correct statistics
                        var biasCorrectedStats = new double[_numStats];
                        for (int i = 0; i < _numStats; i++)
                        {
                            double innerMean = statsInnerSum[i] / validInner;
                            biasCorrectedStats[i] = outerStat[i] - (innerMean - outerStat[i]);
                        }

                        _bootstrapParameterSets[idx] = new ParameterSet(biasCorrectedParms, outerFit.Fitness);
                        for (int k = 0; k < _numStats; k++)
                            _bootstrapStatistics[idx, k] = biasCorrectedStats[k];
                        _validFlags[idx] = true;
                        succeeded = true;
                    }
                    catch (Exception)
                    {
                        // retry
                    }
                    if (succeeded) break;
                }

                if (!succeeded)
                {
                    MarkFailed(idx);
                }
            });
        }

        /// <summary>
        /// Runs the bootstrap procedure with nested inner bootstrap for studentized (Bootstrap-t) confidence intervals.
        /// Must be called before requesting Bootstrap-t CIs.
        /// </summary>
        public void RunWithStudentizedBootstrap()
        {
            ValidateCoreDelegates();
            var resample = ResampleFunction!;
            var fitFunc = FitFunction!;
            var statistic = StatisticFunction!;

            var originalStats = statistic(_originalParameters);
            _numStats = originalStats.Length;
            _numParams = _originalParameters.Values.Length;

            // Apply transform to population statistics
            var popTransformed = new double[_numStats];
            for (int i = 0; i < _numStats; i++)
                popTransformed[i] = ApplyTransform(originalStats[i]);

            _bootstrapParameterSets = new ParameterSet[Replicates];
            _bootstrapStatistics = new double[Replicates, _numStats];
            _validFlags = new bool[Replicates];
            _failedCount = 0;
            _studentizedValues = new double[Replicates, _numStats];
            _transformedStatistics = new double[Replicates, _numStats];
            var studentizedValues = _studentizedValues;
            var transformedStatistics = _transformedStatistics;

            var prng = new MersenneTwister(PRNGSeed);
            var seeds = prng.NextIntegers(Replicates);

            Parallel.For(0, Replicates, idx =>
            {
                bool succeeded = false;
                for (int retry = 0; retry < MaxRetries; retry++)
                {
                    try
                    {
                        var rng = new MersenneTwister(seeds[idx] + 10 * retry);
                        var sample = resample(_originalData, _originalParameters, rng);
                        var outerFit = fitFunc(sample);
                        var outerStats = statistic(outerFit);

                        _bootstrapParameterSets[idx] = outerFit;
                        for (int k = 0; k < _numStats; k++)
                            _bootstrapStatistics[idx, k] = outerStats[k];

                        // Transform outer statistics
                        var outerTransformed = new double[_numStats];
                        for (int j = 0; j < _numStats; j++)
                            outerTransformed[j] = ApplyTransform(outerStats[j]);

                        // Inner bootstrap for SE estimation
                        var innerPrng = new MersenneTwister(seeds[idx]);
                        var innerSeeds = innerPrng.NextIntegers(InnerReplicates);
                        var innerTransformed = new double[InnerReplicates, _numStats];
                        int validInner = 0;

                        for (int k = 0; k < InnerReplicates; k++)
                        {
                            try
                            {
                                var innerSample = resample(sample, outerFit, new MersenneTwister(innerSeeds[k]));
                                var innerFit = fitFunc(innerSample);
                                var innerStats = statistic(innerFit);
                                for (int j = 0; j < _numStats; j++)
                                    innerTransformed[k, j] = ApplyTransform(innerStats[j]);
                                validInner++;
                            }
                            catch (Exception)
                            {
                                for (int j = 0; j < _numStats; j++)
                                    innerTransformed[k, j] = double.NaN;
                            }
                        }

                        if (validInner < 2) continue;

                        // Compute inner SE per statistic and studentized values
                        for (int j = 0; j < _numStats; j++)
                        {
                            var col = innerTransformed.GetColumn(j);
                            var validCol = col.Where(x => !double.IsNaN(x) && !double.IsInfinity(x)).ToArray();
                            double se = validCol.Length > 1 ? Statistics.StandardDeviation(validCol) : double.NaN;
                            transformedStatistics[idx, j] = outerTransformed[j];
                            studentizedValues[idx, j] = se > 0 ? (popTransformed[j] - outerTransformed[j]) / se : double.NaN;
                        }

                        _validFlags[idx] = true;
                        succeeded = true;
                    }
                    catch (Exception)
                    {
                        // retry
                    }
                    if (succeeded) break;
                }

                if (!succeeded)
                {
                    MarkFailed(idx);
                    for (int j = 0; j < _numStats; j++)
                    {
                        transformedStatistics[idx, j] = double.NaN;
                        studentizedValues[idx, j] = double.NaN;
                    }
                }
            });
        }

        #endregion

        #region Confidence Intervals

        /// <summary>
        /// Computes bootstrap confidence intervals using the specified method.
        /// </summary>
        /// <param name="method">The confidence interval method.</param>
        /// <param name="alpha">The confidence level. Default = 0.1, resulting in 90% CIs.</param>
        /// <returns>A BootstrapResults object containing CIs for both parameters and statistics.</returns>
        public BootstrapResults GetConfidenceIntervals(BootstrapCIMethod method, double alpha = 0.1)
        {
            if (_bootstrapStatistics == null)
                throw new InvalidOperationException("Run() or RunWithStudentizedBootstrap() must be called first.");
            if (alpha <= 0 || alpha >= 1)
                throw new ArgumentOutOfRangeException(nameof(alpha), "Alpha must be between 0 and 1.");
            if (method == BootstrapCIMethod.BCa && (JackknifeFunction == null || SampleSizeFunction == null))
                throw new InvalidOperationException("JackknifeFunction and SampleSizeFunction must be set for BCa method.");
            if (method == BootstrapCIMethod.BootstrapT && _studentizedValues == null)
                throw new InvalidOperationException("RunWithStudentizedBootstrap() must be called before requesting Bootstrap-t CIs.");

            if (StatisticFunction == null)
                throw new InvalidOperationException("StatisticFunction must be set.");
            var originalStats = StatisticFunction(_originalParameters);

            // Compute acceleration constants once for BCa
            double[]? accelConstants = null;
            if (method == BootstrapCIMethod.BCa)
                accelConstants = ComputeAccelerationConstants(originalStats);

            var results = new BootstrapResults
            {
                Method = method,
                Alpha = alpha,
                StatisticResults = new BootstrapStatisticResult[_numStats],
                ParameterResults = new BootstrapStatisticResult[_numParams],
                FailedReplicates = _failedCount
            };

            // Compute CIs for each statistic
            for (int i = 0; i < _numStats; i++)
            {
                var values = _bootstrapStatistics.GetColumn(i);
                switch (method)
                {
                    case BootstrapCIMethod.Percentile:
                        results.StatisticResults[i] = ComputePercentileCI(values, originalStats[i], alpha);
                        break;
                    case BootstrapCIMethod.BiasCorrected:
                        results.StatisticResults[i] = ComputeBiasCorrectedCI(values, originalStats[i], alpha);
                        break;
                    case BootstrapCIMethod.BCa:
                        results.StatisticResults[i] = ComputeBCaCI(values, originalStats[i], alpha, accelConstants![i]);
                        break;
                    case BootstrapCIMethod.Normal:
                        results.StatisticResults[i] = ComputeNormalCI(values, originalStats[i], alpha);
                        break;
                    case BootstrapCIMethod.BootstrapT:
                        results.StatisticResults[i] = ComputeBootstrapTCI(i, originalStats[i], alpha);
                        break;
                }
            }

            // Compute percentile CIs for each parameter
            for (int i = 0; i < _numParams; i++)
            {
                var values = _bootstrapParameterSets.Select(ps => ps.Values[i]).ToArray();
                results.ParameterResults[i] = ComputePercentileCI(values, _originalParameters.Values[i], alpha);
            }

            return results;
        }

        #endregion

        #region CI Methods

        /// <summary>
        /// Computes percentile confidence intervals for a single statistic.
        /// </summary>
        private BootstrapStatisticResult ComputePercentileCI(double[] values, double populationEstimate, double alpha)
        {
            var validValues = values.Where(x => !double.IsNaN(x) && !double.IsInfinity(x)).ToArray();
            Array.Sort(validValues);

            double lowerP = alpha / 2d;
            double upperP = 1d - alpha / 2d;

            return new BootstrapStatisticResult
            {
                PopulationEstimate = populationEstimate,
                LowerCI = validValues.Length > 0 ? Statistics.Percentile(validValues, lowerP, true) : double.NaN,
                UpperCI = validValues.Length > 0 ? Statistics.Percentile(validValues, upperP, true) : double.NaN,
                ValidCount = validValues.Length,
                TotalCount = values.Length,
                StandardError = validValues.Length > 1 ? Statistics.StandardDeviation(validValues) : double.NaN,
                Mean = validValues.Length > 0 ? Statistics.Mean(validValues) : double.NaN
            };
        }

        /// <summary>
        /// Computes bias-corrected (BC) confidence intervals for a single statistic.
        /// </summary>
        private BootstrapStatisticResult ComputeBiasCorrectedCI(double[] values, double populationEstimate, double alpha)
        {
            var validValues = values.Where(x => !double.IsNaN(x) && !double.IsInfinity(x)).ToArray();
            int validN = validValues.Length;
            if (validN == 0) return EmptyResult(populationEstimate, values.Length);

            // Count proportion <= population estimate
            int countLeq = 0;
            for (int i = 0; i < validN; i++)
                if (validValues[i] <= populationEstimate) countLeq++;
            double P0 = (double)countLeq / (validN + 1);

            Array.Sort(validValues);

            double Z0 = Normal.StandardZ(P0);
            double ZLower = Normal.StandardZ(alpha / 2d);
            double ZUpper = Normal.StandardZ(1d - alpha / 2d);
            double bcLower = Normal.StandardCDF(2d * Z0 + ZLower);
            double bcUpper = Normal.StandardCDF(2d * Z0 + ZUpper);

            return new BootstrapStatisticResult
            {
                PopulationEstimate = populationEstimate,
                LowerCI = Statistics.Percentile(validValues, bcLower, true),
                UpperCI = Statistics.Percentile(validValues, bcUpper, true),
                ValidCount = validN,
                TotalCount = values.Length,
                StandardError = validN > 1 ? Statistics.StandardDeviation(validValues) : double.NaN,
                Mean = Statistics.Mean(validValues)
            };
        }

        /// <summary>
        /// Computes bias-corrected and accelerated (BCa) confidence intervals for a single statistic.
        /// </summary>
        private BootstrapStatisticResult ComputeBCaCI(double[] values, double populationEstimate, double alpha, double acceleration)
        {
            var validValues = values.Where(x => !double.IsNaN(x) && !double.IsInfinity(x)).ToArray();
            int validN = validValues.Length;
            if (validN == 0) return EmptyResult(populationEstimate, values.Length);

            // Count proportion <= population estimate (matching BCaQuantileCI line 593)
            int countLeq = 0;
            for (int i = 0; i < validN; i++)
                if (validValues[i] <= populationEstimate) countLeq++;
            double P0 = (double)(countLeq + 1) / (validN + 1);

            Array.Sort(validValues);

            double Z0 = Normal.StandardZ(P0);
            double ZLower = Normal.StandardZ(alpha / 2d);
            double ZUpper = Normal.StandardZ(1d - alpha / 2d);

            double numLower = Z0 + ZLower;
            double denLower = 1d - acceleration * numLower;
            double bcLower = Normal.StandardCDF(Z0 + numLower / denLower);

            double numUpper = Z0 + ZUpper;
            double denUpper = 1d - acceleration * numUpper;
            double bcUpper = Normal.StandardCDF(Z0 + numUpper / denUpper);

            return new BootstrapStatisticResult
            {
                PopulationEstimate = populationEstimate,
                LowerCI = Statistics.Percentile(validValues, bcLower, true),
                UpperCI = Statistics.Percentile(validValues, bcUpper, true),
                ValidCount = validN,
                TotalCount = values.Length,
                StandardError = validN > 1 ? Statistics.StandardDeviation(validValues) : double.NaN,
                Mean = Statistics.Mean(validValues)
            };
        }

        /// <summary>
        /// Computes Normal (standard) confidence intervals for a single statistic.
        /// Uses a configurable transform for transformation invariance.
        /// </summary>
        private BootstrapStatisticResult ComputeNormalCI(double[] values, double populationEstimate, double alpha)
        {
            double popTransformed = ApplyTransform(populationEstimate);

            var transformedValid = new double[values.Length];
            int validCount = 0;
            for (int i = 0; i < values.Length; i++)
            {
                if (!double.IsNaN(values[i]) && !double.IsInfinity(values[i]))
                {
                    transformedValid[validCount] = ApplyTransform(values[i]);
                    validCount++;
                }
            }

            if (validCount < 2) return EmptyResult(populationEstimate, values.Length);

            var validSlice = new double[validCount];
            Array.Copy(transformedValid, validSlice, validCount);

            double SE = Statistics.StandardDeviation(validSlice);
            double ZLower = Normal.StandardZ(alpha / 2d);
            double ZUpper = Normal.StandardZ(1d - alpha / 2d);

            double lowerTransformed = popTransformed + SE * ZLower;
            double upperTransformed = popTransformed + SE * ZUpper;

            return new BootstrapStatisticResult
            {
                PopulationEstimate = populationEstimate,
                LowerCI = ApplyInverseTransform(lowerTransformed),
                UpperCI = ApplyInverseTransform(upperTransformed),
                ValidCount = validCount,
                TotalCount = values.Length,
                StandardError = SE,
                Mean = Statistics.Mean(validSlice)
            };
        }

        /// <summary>
        /// Computes Bootstrap-t (studentized) confidence intervals for a single statistic.
        /// Requires RunWithStudentizedBootstrap() to have been called first.
        /// </summary>
        private BootstrapStatisticResult ComputeBootstrapTCI(int statisticIndex, double populationEstimate, double alpha)
        {
            double popTransformed = ApplyTransform(populationEstimate);

            // GetConfidenceIntervals() validates _studentizedValues is non-null before calling this method
            var xCol = _transformedStatistics!.GetColumn(statisticIndex);
            var tCol = _studentizedValues!.GetColumn(statisticIndex);

            var validX = xCol.Where(x => !double.IsNaN(x) && !double.IsInfinity(x)).ToArray();
            var validT = tCol.Where(x => !double.IsNaN(x) && !double.IsInfinity(x)).ToArray();

            if (validT.Length < 2) return EmptyResult(populationEstimate, Replicates);

            double SE = Statistics.StandardDeviation(validX);
            Array.Sort(validT);

            double tLower = Statistics.Percentile(validT, alpha / 2d, true);
            double tUpper = Statistics.Percentile(validT, 1d - alpha / 2d, true);

            return new BootstrapStatisticResult
            {
                PopulationEstimate = populationEstimate,
                LowerCI = ApplyInverseTransform(popTransformed + SE * tLower),
                UpperCI = ApplyInverseTransform(popTransformed + SE * tUpper),
                ValidCount = validT.Length,
                TotalCount = Replicates,
                StandardError = SE,
                Mean = Statistics.Mean(validX)
            };
        }

        #endregion

        #region BCa Support

        /// <summary>
        /// Computes the acceleration constants for each statistic using jackknife leave-one-out.
        /// </summary>
        private double[] ComputeAccelerationConstants(double[] populationEstimates)
        {
            // Caller (GetConfidenceIntervals) validates these delegates before calling this method
            var sampleSize = SampleSizeFunction!;
            var jackknife = JackknifeFunction!;
            var fitFunc = FitFunction!;
            var statistic = StatisticFunction!;

            int N = sampleSize(_originalData);
            var I2 = new double[_numStats];
            var I3 = new double[_numStats];
            var a = new double[_numStats];

            Parallel.For(0, N, idx =>
            {
                try
                {
                    var jackData = jackknife(_originalData, idx);
                    var jackFit = fitFunc(jackData);
                    var jackStats = statistic(jackFit);

                    for (int i = 0; i < _numStats; i++)
                    {
                        double diff = populationEstimates[i] - jackStats[i];
                        Tools.ParallelAdd(ref I2[i], diff * diff);
                        Tools.ParallelAdd(ref I3[i], diff * diff * diff);
                    }
                }
                catch (Exception)
                {
                    // Skip failed jackknife samples
                }
            });

            for (int i = 0; i < _numStats; i++)
                a[i] = I3[i] / (Math.Pow(I2[i], 1.5) * 6d);

            return a;
        }

        #endregion

        #region Private Helpers

        /// <summary>
        /// Validates that the core delegate functions are set.
        /// </summary>
        private void ValidateCoreDelegates()
        {
            if (ResampleFunction == null)
                throw new InvalidOperationException("ResampleFunction must be set.");
            if (FitFunction == null)
                throw new InvalidOperationException("FitFunction must be set.");
            if (StatisticFunction == null)
                throw new InvalidOperationException("StatisticFunction must be set.");
        }

        /// <summary>
        /// Initializes internal state arrays before a bootstrap run.
        /// </summary>
        private void InitializeState()
        {
            // ValidateCoreDelegates() is always called before this method
            var originalStats = StatisticFunction!(_originalParameters);
            _numStats = originalStats.Length;
            _numParams = _originalParameters.Values.Length;

            _bootstrapParameterSets = new ParameterSet[Replicates];
            _bootstrapStatistics = new double[Replicates, _numStats];
            _validFlags = new bool[Replicates];
            _failedCount = 0;
            _studentizedValues = null;
            _transformedStatistics = null;
        }

        /// <summary>
        /// Marks a replicate as failed with NaN values.
        /// </summary>
        private void MarkFailed(int idx)
        {
            var nanParams = new double[_numParams];
            for (int k = 0; k < _numParams; k++) nanParams[k] = double.NaN;
            _bootstrapParameterSets[idx] = new ParameterSet(nanParams, double.NaN);
            for (int k = 0; k < _numStats; k++)
                _bootstrapStatistics[idx, k] = double.NaN;
            _validFlags[idx] = false;
            Interlocked.Increment(ref _failedCount);
        }

        /// <summary>
        /// Applies the Transform function, or returns the value unchanged if Transform is null.
        /// </summary>
        private double ApplyTransform(double value)
        {
            return Transform != null ? Transform(value) : value;
        }

        /// <summary>
        /// Applies the InverseTransform function, or returns the value unchanged if InverseTransform is null.
        /// </summary>
        private double ApplyInverseTransform(double value)
        {
            return InverseTransform != null ? InverseTransform(value) : value;
        }

        /// <summary>
        /// Creates an empty result for when there are insufficient valid values.
        /// </summary>
        private BootstrapStatisticResult EmptyResult(double populationEstimate, int totalCount)
        {
            return new BootstrapStatisticResult
            {
                PopulationEstimate = populationEstimate,
                LowerCI = double.NaN,
                UpperCI = double.NaN,
                ValidCount = 0,
                TotalCount = totalCount,
                StandardError = double.NaN,
                Mean = double.NaN
            };
        }

        #endregion
    }
}
