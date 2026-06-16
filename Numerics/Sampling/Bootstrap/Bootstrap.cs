using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Functions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Numerics.Sampling
{

    /// <summary>
    /// A general-purpose bootstrap class for parametric or non-parametric bootstrap analysis.
    /// Supports Percentile, Bias-Corrected, BCa, Normal, Bootstrap-t, and covariance-aware pivotal bootstrap workflows.
    /// </summary>
    /// <typeparam name="TData">The type of data being bootstrapped.</typeparam>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    ///     The regular bootstrap methods use <see cref="FitFunction"/> and require <see cref="StatisticFunction"/>.
    ///     The pivotal bootstrap is a separate covariance-aware run mode that uses
    ///     <see cref="FitWithCovarianceFunction"/> and <see cref="OriginalCovariance"/>. The two modes share
    ///     result storage, but their required delegates and supported confidence interval methods are validated separately.
    /// </para>
    /// <para>
    /// <see href="https://en.wikipedia.org/wiki/Bootstrapping_(statistics)" />
    /// </para>
    /// </remarks>
    public class Bootstrap<TData>
    {
        #region Construction

        /// <summary>
        /// Constructs a new regular bootstrap analysis.
        /// </summary>
        /// <param name="originalData">The original data.</param>
        /// <param name="originalParameters">The original fitted parameter set.</param>
        public Bootstrap(TData originalData, ParameterSet originalParameters)
        {
            _originalData = originalData;
            _originalParameters = originalParameters;
        }

        /// <summary>
        /// Constructs a new covariance-aware bootstrap analysis.
        /// </summary>
        /// <param name="originalData">The original data.</param>
        /// <param name="originalFit">The original fitted parameter set and covariance matrix.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="originalFit"/> is null.</exception>
        public Bootstrap(TData originalData, BootstrapFit originalFit)
        {
            if (originalFit == null)
                throw new ArgumentNullException(nameof(originalFit));

            _originalData = originalData;
            _originalParameters = originalFit.Parameters.Clone();
            _originalCovariance = originalFit.Covariance.Clone();
        }

        #endregion

        #region Members

        /// <summary>
        /// Identifies the workflow that produced the active bootstrap results.
        /// </summary>
        private enum BootstrapRunType
        {
            /// <summary>
            /// No bootstrap workflow has completed.
            /// </summary>
            None,

            /// <summary>
            /// Results came from <see cref="Run"/>.
            /// </summary>
            Regular,

            /// <summary>
            /// Results came from <see cref="RunDoubleBootstrap(int)"/>.
            /// </summary>
            DoubleBootstrap,

            /// <summary>
            /// Results came from <see cref="RunWithStudentizedBootstrap"/>.
            /// </summary>
            Studentized,

            /// <summary>
            /// Results came from <see cref="RunPivotalBootstrap"/> or <see cref="TransformPivotalBootstrap(IEnumerable{BootstrapFit})"/>.
            /// </summary>
            Pivotal
        }

        private readonly TData _originalData;
        private ParameterSet _originalParameters;
        private Matrix? _originalCovariance;
        private ParameterSet[] _bootstrapParameterSets = null!;
        private double[,] _bootstrapStatistics = null!;
        private ParameterSet[] _rawBootstrapParameterSets = Array.Empty<ParameterSet>();
        private double[,]? _rawBootstrapStatistics;
        private BootstrapFit[] _rawBootstrapFits = Array.Empty<BootstrapFit>();
        private ILinkFunction?[] _pivotalLinks = Array.Empty<ILinkFunction?>();
        private PivotalBootstrapDiagnostics? _pivotalDiagnostics;
        private double[]? _parentStatistics;
        private int _numStats;
        private int _numParams;
        private int _failedCount;
        private bool[] _validFlags = null!;
        private double[,]? _studentizedValues;
        private double[,]? _transformedStatistics;
        private BootstrapRunType _runType;

        #endregion

        #region Properties

        /// <summary>
        /// Delegate function for resampling the original data given the current parameters and a random number generator.
        /// </summary>
        public Func<TData, ParameterSet, Random, TData>? ResampleFunction { get; set; }

        /// <summary>
        /// Delegate function for fitting a model to data and returning a parameter set.
        /// Required by regular bootstrap methods.
        /// </summary>
        public Func<TData, ParameterSet>? FitFunction { get; set; }

        /// <summary>
        /// Delegate function for fitting a model to data and returning a parameter set plus covariance matrix.
        /// Required by <see cref="RunPivotalBootstrap"/>.
        /// </summary>
        public Func<TData, BootstrapFit>? FitWithCovarianceFunction { get; set; }

        /// <summary>
        /// Delegate function for extracting statistics from a fitted parameter set.
        /// Required by regular bootstrap methods and optional for pivotal bootstrap.
        /// </summary>
        public Func<ParameterSet, double[]>? StatisticFunction { get; set; }

        /// <summary>
        /// Optional delegate for computing a leave-one-out jackknife sample.
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
        /// Default is cube-root.
        /// </summary>
        public Func<double, double> Transform { get; set; } = x => Math.Pow(x, 1d / 3d);

        /// <summary>
        /// Optional inverse transform corresponding to <see cref="Transform"/>.
        /// Default is cube.
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
        /// Gets or sets the original covariance matrix used by the pivotal bootstrap.
        /// </summary>
        public Matrix? OriginalCovariance
        {
            get => _originalCovariance?.Clone();
            set => _originalCovariance = value?.Clone();
        }

        /// <summary>
        /// Gets or sets a factory that returns one optional <see cref="ILinkFunction"/> per fitted parameter for pivotal bootstrap.
        /// Null link entries are treated as identity links by <see cref="LinkController"/>.
        /// </summary>
        public Func<PivotalBootstrapContext, ILinkFunction?[]>? PivotalLinkFactory { get; set; }

        /// <summary>
        /// Gets or sets an optional filter applied to valid raw covariance-aware bootstrap fits before pivotal transformation.
        /// </summary>
        public Func<BootstrapFit, bool>? PivotalReplicateFilter { get; set; }

        /// <summary>
        /// Gets or sets an optional validator applied to transformed pivotal parameter values.
        /// </summary>
        public Func<double[], bool>? PivotalParameterValidator { get; set; }

        /// <summary>
        /// Gets or sets the policy used when a pivotal draw is invalid. Default is <see cref="PivotalBootstrapInvalidDrawPolicy.Drop"/>.
        /// </summary>
        public PivotalBootstrapInvalidDrawPolicy PivotalInvalidDrawPolicy { get; set; } = PivotalBootstrapInvalidDrawPolicy.Drop;

        /// <summary>
        /// Gets or sets a value indicating whether pivotal bootstrap covariance matrices should be regularized before Cholesky decomposition.
        /// </summary>
        public bool RegularizePivotalCovariances { get; set; } = true;

        /// <summary>
        /// Gets or sets an optional absolute component limit for the standardized pivotal vector.
        /// </summary>
        public double? PivotalZLimit { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether Gaussian jitter should be added to pivotal vectors before the optional limit check.
        /// </summary>
        public bool AddPivotalJitter { get; set; }

        /// <summary>
        /// Gets or sets the base standard deviation of optional Gaussian jitter. The applied scale is this value divided by sqrt(parameter count).
        /// </summary>
        public double PivotalJitterScale { get; set; } = 0.01d;

        /// <summary>
        /// Gets the active bootstrapped model parameter sets.
        /// For pivotal bootstrap, these are the retained pivotal parameter draws.
        /// </summary>
        public ParameterSet[] BootstrapParameterSets => _bootstrapParameterSets;

        /// <summary>
        /// Gets the active bootstrapped statistics as a 2D array [replicate, statistic].
        /// For pivotal bootstrap, these are pivotal statistics when <see cref="StatisticFunction"/> was supplied.
        /// </summary>
        public double[,] BootstrapStatistics => _bootstrapStatistics;

        /// <summary>
        /// Gets the accepted raw covariance-aware bootstrap fits from the most recent pivotal bootstrap run.
        /// </summary>
        public BootstrapFit[] RawBootstrapFits => _rawBootstrapFits;

        /// <summary>
        /// Gets the accepted raw covariance-aware bootstrap parameter sets from the most recent pivotal bootstrap run.
        /// </summary>
        public ParameterSet[] RawBootstrapParameterSets => _rawBootstrapParameterSets;

        /// <summary>
        /// Gets the accepted raw bootstrap statistics from the most recent pivotal bootstrap run, when <see cref="StatisticFunction"/> was supplied.
        /// </summary>
        public double[,]? RawBootstrapStatistics => _rawBootstrapStatistics;

        /// <summary>
        /// Gets the pivotal link functions used by the most recent pivotal bootstrap run.
        /// </summary>
        public ILinkFunction?[] PivotalLinks => (ILinkFunction?[])_pivotalLinks.Clone();

        /// <summary>
        /// Gets diagnostics from the most recent pivotal bootstrap run.
        /// </summary>
        public PivotalBootstrapDiagnostics? PivotalDiagnostics => _pivotalDiagnostics;

        /// <summary>
        /// Gets the number of replicates that failed after all retries.
        /// For pivotal bootstrap, this is the raw covariance-aware fit failure count.
        /// </summary>
        public int FailedReplicates => _failedCount;

        #endregion

        #region Run Methods

        /// <summary>
        /// Runs the regular bootstrap procedure with error handling and retry logic.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when a required regular-bootstrap delegate has not been set.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <see cref="Replicates"/> or <see cref="MaxRetries"/> is not positive.</exception>
        public void Run()
        {
            ValidateCoreDelegates();
            ValidateReplicationSettings();
            var resample = ResampleFunction!;
            var fit = FitFunction!;
            var statistic = StatisticFunction!;
            InitializeState();
            ResetPivotalState();

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
                        if (!HasExpectedFiniteParameterValues(fitResult, _numParams))
                            continue;
                        var stat = ValidateStatistics(statistic(fitResult), _numStats);

                        _bootstrapParameterSets[idx] = fitResult;
                        for (int k = 0; k < _numStats; k++)
                            _bootstrapStatistics[idx, k] = stat[k];
                        _validFlags[idx] = true;
                        succeeded = true;
                    }
                    catch (Exception)
                    {
                        // Retry failed regular bootstrap replicates.
                    }
                    if (succeeded) break;
                }

                if (!succeeded)
                    MarkFailed(idx);
            });

            _runType = BootstrapRunType.Regular;
        }

        /// <summary>
        /// Runs the double bootstrap procedure with bias correction.
        /// </summary>
        /// <param name="innerReplicates">Number of inner bootstrap replicates. Default = 300.</param>
        /// <exception cref="InvalidOperationException">Thrown when a required regular-bootstrap delegate has not been set.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when replicate counts are not positive.</exception>
        public void RunDoubleBootstrap(int innerReplicates = 300)
        {
            ValidateCoreDelegates();
            ValidateReplicationSettings();
            if (innerReplicates < 1)
                throw new ArgumentOutOfRangeException(nameof(innerReplicates), "The number of inner replicates must be positive.");

            var resample = ResampleFunction!;
            var fit = FitFunction!;
            var statistic = StatisticFunction!;
            InitializeState();
            ResetPivotalState();

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

                        var outerSample = resample(_originalData, _originalParameters, rng);
                        var outerFit = fit(outerSample);
                        if (!HasExpectedFiniteParameterValues(outerFit, _numParams))
                            continue;
                        var outerStat = ValidateStatistics(statistic(outerFit), _numStats);

                        int p = _numParams;
                        var parmsInnerSum = new double[p];
                        var statsInnerSum = new double[_numStats];
                        int validInner = 0;

                        for (int k = 0; k < innerReplicates; k++)
                        {
                            try
                            {
                                var innerSample = resample(outerSample, outerFit, rng);
                                var innerFit = fit(innerSample);
                                if (!HasExpectedFiniteParameterValues(innerFit, p))
                                    continue;
                                var innerStat = ValidateStatistics(statistic(innerFit), _numStats);

                                for (int i = 0; i < p; i++)
                                    parmsInnerSum[i] += innerFit.Values[i];
                                for (int i = 0; i < _numStats; i++)
                                    statsInnerSum[i] += innerStat[i];
                                validInner++;
                            }
                            catch (Exception)
                            {
                                // Skip failed inner replicate.
                            }
                        }

                        if (validInner == 0) continue;

                        var biasCorrectedParms = new double[p];
                        for (int i = 0; i < p; i++)
                        {
                            double innerMean = parmsInnerSum[i] / validInner;
                            biasCorrectedParms[i] = outerFit.Values[i] - (innerMean - outerFit.Values[i]);
                        }

                        var biasCorrectedStats = new double[_numStats];
                        for (int i = 0; i < _numStats; i++)
                        {
                            double innerMean = statsInnerSum[i] / validInner;
                            biasCorrectedStats[i] = outerStat[i] - (innerMean - outerStat[i]);
                        }

                        _bootstrapParameterSets[idx] = new ParameterSet(biasCorrectedParms, outerFit.Fitness, outerFit.Weight);
                        for (int k = 0; k < _numStats; k++)
                            _bootstrapStatistics[idx, k] = biasCorrectedStats[k];
                        _validFlags[idx] = true;
                        succeeded = true;
                    }
                    catch (Exception)
                    {
                        // Retry failed outer replicate.
                    }
                    if (succeeded) break;
                }

                if (!succeeded)
                    MarkFailed(idx);
            });

            _runType = BootstrapRunType.DoubleBootstrap;
        }

        /// <summary>
        /// Runs the regular bootstrap procedure with nested inner bootstrap for studentized Bootstrap-t confidence intervals.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when a required regular-bootstrap delegate has not been set.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when replicate counts are not positive.</exception>
        public void RunWithStudentizedBootstrap()
        {
            ValidateCoreDelegates();
            ValidateReplicationSettings();
            if (InnerReplicates < 1)
                throw new ArgumentOutOfRangeException(nameof(InnerReplicates), "The number of inner replicates must be positive.");

            var resample = ResampleFunction!;
            var fitFunc = FitFunction!;
            var statistic = StatisticFunction!;
            ResetPivotalState();

            var originalStats = statistic(_originalParameters);
            ValidateStatistics(originalStats);
            _parentStatistics = (double[])originalStats.Clone();
            _numStats = originalStats.Length;
            _numParams = _originalParameters.Values.Length;

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
                        if (!HasExpectedFiniteParameterValues(outerFit, _numParams))
                            continue;
                        var outerStats = ValidateStatistics(statistic(outerFit), _numStats);

                        _bootstrapParameterSets[idx] = outerFit;
                        for (int k = 0; k < _numStats; k++)
                            _bootstrapStatistics[idx, k] = outerStats[k];

                        var outerTransformed = new double[_numStats];
                        for (int j = 0; j < _numStats; j++)
                            outerTransformed[j] = ApplyTransform(outerStats[j]);

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
                                if (!HasExpectedFiniteParameterValues(innerFit, _numParams))
                                {
                                    for (int j = 0; j < _numStats; j++)
                                        innerTransformed[k, j] = double.NaN;
                                    continue;
                                }
                                var innerStats = ValidateStatistics(statistic(innerFit), _numStats);

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

                        for (int j = 0; j < _numStats; j++)
                        {
                            var col = innerTransformed.GetColumn(j);
                            var validCol = col.Where(Tools.IsFinite).ToArray();
                            double se = validCol.Length > 1 ? Statistics.StandardDeviation(validCol) : double.NaN;
                            transformedStatistics[idx, j] = outerTransformed[j];
                            studentizedValues[idx, j] = se > 0 ? (popTransformed[j] - outerTransformed[j]) / se : double.NaN;
                        }

                        _validFlags[idx] = true;
                        succeeded = true;
                    }
                    catch (Exception)
                    {
                        // Retry failed outer replicate.
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

            _runType = BootstrapRunType.Studentized;
        }

        #endregion

        #region Pivotal Bootstrap

        /// <summary>
        /// Runs the covariance-aware pivotal bootstrap as a distinct bootstrap mode.
        /// </summary>
        /// <exception cref="InvalidOperationException">
        /// Thrown when <see cref="ResampleFunction"/>, <see cref="FitWithCovarianceFunction"/>, or
        /// <see cref="OriginalCovariance"/> has not been supplied.
        /// </exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <see cref="Replicates"/> or <see cref="MaxRetries"/> is not positive.</exception>
        public void RunPivotalBootstrap()
        {
            ValidatePivotalDelegates();
            ValidateReplicationSettings();

            var parentFit = CreateOriginalFit();
            var resample = ResampleFunction!;
            var fit = FitWithCovarianceFunction!;

            var stopwatch = Stopwatch.StartNew();
            var prng = new MersenneTwister(PRNGSeed);
            var seeds = prng.NextIntegers(Replicates);
            var rawFits = new BootstrapFit?[Replicates];
            int failed = 0;

            Parallel.For(0, Replicates, index =>
            {
                bool succeeded = false;
                for (int retry = 0; retry < MaxRetries; retry++)
                {
                    try
                    {
                        var rng = new MersenneTwister(seeds[index] + 10 * retry);
                        TData sample = resample(_originalData, parentFit.Parameters, rng);
                        BootstrapFit rawFit = fit(sample);
                        if (!IsValidFit(rawFit, parentFit.ParameterCount))
                            continue;

                        rawFits[index] = rawFit;
                        succeeded = true;
                    }
                    catch (Exception)
                    {
                        // Retry failed pivotal bootstrap resamples or covariance-aware fits.
                    }

                    if (succeeded)
                        break;
                }

                if (!succeeded)
                    Interlocked.Increment(ref failed);
            });

            stopwatch.Stop();

            var acceptedRawFits = rawFits.Where(f => f != null).Select(f => f!).ToArray();
            TransformPivotalBootstrap(acceptedRawFits, Replicates, failed, stopwatch.Elapsed);
        }

        /// <summary>
        /// Transforms precomputed raw covariance-aware bootstrap fits into pivotal bootstrap draws.
        /// </summary>
        /// <param name="rawFits">Raw bootstrap fits and covariances to transform.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="rawFits"/> is null.</exception>
        /// <exception cref="InvalidOperationException">Thrown when no valid raw fits are accepted.</exception>
        public void TransformPivotalBootstrap(IEnumerable<BootstrapFit> rawFits)
        {
            if (rawFits == null)
                throw new ArgumentNullException(nameof(rawFits));

            BootstrapFit[] suppliedFits = rawFits.ToArray();
            TransformPivotalBootstrap(suppliedFits, suppliedFits.Length, 0, TimeSpan.Zero);
        }

        /// <summary>
        /// Applies the pivotal transformation to raw fits and stores the resulting active bootstrap ensemble.
        /// </summary>
        /// <param name="rawFits">Raw bootstrap fits supplied to the transformation.</param>
        /// <param name="requestedReplicates">The number of raw replicates requested or supplied.</param>
        /// <param name="failedRawReplicates">The number of raw replicates that failed before transformation.</param>
        /// <param name="resamplingTime">The elapsed raw resampling and fitting time.</param>
        /// <exception cref="InvalidOperationException">Thrown when no valid raw fits are accepted.</exception>
        private void TransformPivotalBootstrap(BootstrapFit[] rawFits, int requestedReplicates, int failedRawReplicates, TimeSpan resamplingTime)
        {
            BootstrapFit parentFit = CreateOriginalFit();
            int p = parentFit.ParameterCount;
            var stopwatch = Stopwatch.StartNew();

            BootstrapFit[] acceptedRawFits = rawFits
                .Where(fit => IsAcceptedRawFit(fit, p))
                .ToArray();

            if (acceptedRawFits.Length == 0)
                throw new InvalidOperationException("No valid raw bootstrap fits were supplied for pivotal transformation.");

            ILinkFunction?[] links = CreatePivotalLinks(parentFit, acceptedRawFits);
            var linkController = new LinkController(links);
            double[] parentEta = linkController.Link(parentFit.Parameters.Values);
            ValidateTransformedValues(parentEta, "The parent link transformation produced a non-finite value.");

            Matrix parentLinkCovariance = LinkCovariance(parentFit, linkController);
            var parentCholesky = new CholeskyDecomposition(parentLinkCovariance);
            var pivotalParameterSets = new List<ParameterSet>(acceptedRawFits.Length);
            var jitterRng = new MersenneTwister(PRNGSeed);
            int invalid = 0;

            for (int i = 0; i < acceptedRawFits.Length; i++)
            {
                BootstrapFit rawFit = acceptedRawFits[i];
                if (TryCreatePivotalDraw(rawFit, parentFit, linkController, parentEta, parentCholesky, jitterRng, out ParameterSet pivotalParameters))
                {
                    pivotalParameterSets.Add(pivotalParameters);
                }
                else
                {
                    invalid++;
                    if (TryApplyInvalidPolicy(rawFit, parentFit, out ParameterSet fallbackParameters))
                        pivotalParameterSets.Add(fallbackParameters);
                }
            }

            stopwatch.Stop();

            _rawBootstrapFits = acceptedRawFits;
            _rawBootstrapParameterSets = acceptedRawFits.Select(f => f.Parameters.Clone()).ToArray();
            _bootstrapParameterSets = pivotalParameterSets.ToArray();
            _numParams = p;
            _failedCount = failedRawReplicates;
            _validFlags = Enumerable.Repeat(true, _bootstrapParameterSets.Length).ToArray();
            _studentizedValues = null;
            _transformedStatistics = null;
            _pivotalLinks = links;
            _pivotalDiagnostics = new PivotalBootstrapDiagnostics
            {
                RequestedReplicates = requestedReplicates,
                FailedRawReplicates = failedRawReplicates,
                RejectedRawReplicates = rawFits.Length - acceptedRawFits.Length,
                AcceptedRawReplicates = acceptedRawFits.Length,
                InvalidPivotalReplicates = invalid,
                RetainedPivotalReplicates = _bootstrapParameterSets.Length,
                ResamplingTime = resamplingTime,
                TransformationTime = stopwatch.Elapsed
            };

            if (StatisticFunction != null)
            {
                _parentStatistics = ValidateStatistics(StatisticFunction(parentFit.Parameters));
                _numStats = _parentStatistics.Length;
                _rawBootstrapStatistics = ComputeStatistics(_rawBootstrapParameterSets, StatisticFunction);
                _bootstrapStatistics = ComputeStatistics(_bootstrapParameterSets, StatisticFunction);
            }
            else
            {
                _parentStatistics = Array.Empty<double>();
                _numStats = 0;
                _rawBootstrapStatistics = new double[_rawBootstrapParameterSets.Length, 0];
                _bootstrapStatistics = new double[_bootstrapParameterSets.Length, 0];
            }

            _runType = BootstrapRunType.Pivotal;
        }

        /// <summary>
        /// Attempts to create one pivotal draw from an accepted raw fit.
        /// </summary>
        /// <param name="rawFit">The accepted raw bootstrap fit.</param>
        /// <param name="parentFit">The parent bootstrap fit.</param>
        /// <param name="linkController">The parameter link controller.</param>
        /// <param name="parentEta">The parent parameters in link space.</param>
        /// <param name="parentCholesky">The Cholesky decomposition of the parent link-space covariance.</param>
        /// <param name="jitterRng">The random-number generator used for optional jitter.</param>
        /// <param name="pivotalParameters">The created pivotal parameters when the method returns true.</param>
        /// <returns>True when a finite and valid pivotal draw is created; otherwise false.</returns>
        private bool TryCreatePivotalDraw(
            BootstrapFit rawFit,
            BootstrapFit parentFit,
            LinkController linkController,
            double[] parentEta,
            CholeskyDecomposition parentCholesky,
            Random jitterRng,
            out ParameterSet pivotalParameters)
        {
            pivotalParameters = default;
            try
            {
                double[] rawEta = linkController.Link(rawFit.Parameters.Values);
                ValidateTransformedValues(rawEta, "The raw link transformation produced a non-finite value.");

                Matrix rawLinkCovariance = LinkCovariance(rawFit, linkController);
                var rawCholesky = new CholeskyDecomposition(rawLinkCovariance);
                var difference = new double[parentEta.Length];
                for (int j = 0; j < difference.Length; j++)
                    difference[j] = parentEta[j] - rawEta[j];

                double[] z = rawCholesky.Forward(new Vector(difference)).ToArray();
                if (!ContainsOnlyFiniteValues(z))
                    return false;

                if (AddPivotalJitter && PivotalJitterScale > 0d)
                {
                    double jitterScale = PivotalJitterScale / Math.Sqrt(Math.Max(1, z.Length));
                    for (int j = 0; j < z.Length; j++)
                    {
                        double u = Tools.Clamp(jitterRng.NextDouble(), 1e-16, 1d - 1e-16);
                        z[j] += jitterScale * Normal.StandardZ(u);
                    }
                }

                if (PivotalZLimit.HasValue && z.Any(value => Math.Abs(value) > PivotalZLimit.Value))
                    return false;

                double[] reinflated = parentCholesky.L * z;
                var pivotalEta = new double[parentEta.Length];
                for (int j = 0; j < pivotalEta.Length; j++)
                    pivotalEta[j] = parentEta[j] + reinflated[j];

                double[] pivotalValues = linkController.InverseLink(pivotalEta);
                if (!ContainsOnlyFiniteValues(pivotalValues))
                    return false;
                if (PivotalParameterValidator != null && !PivotalParameterValidator(pivotalValues))
                    return false;

                pivotalParameters = new ParameterSet(pivotalValues, rawFit.Parameters.Fitness, rawFit.Parameters.Weight);
                return true;
            }
            catch
            {
                return false;
            }
        }

        /// <summary>
        /// Applies the configured invalid pivotal draw policy.
        /// </summary>
        /// <param name="rawFit">The raw bootstrap fit associated with the invalid pivotal draw.</param>
        /// <param name="parentFit">The parent bootstrap fit.</param>
        /// <param name="parameterSet">The fallback parameter set when the method returns true.</param>
        /// <returns>True when a fallback parameter set should be retained; otherwise false.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <see cref="PivotalInvalidDrawPolicy"/> is unknown.</exception>
        private bool TryApplyInvalidPolicy(BootstrapFit rawFit, BootstrapFit parentFit, out ParameterSet parameterSet)
        {
            parameterSet = default;
            switch (PivotalInvalidDrawPolicy)
            {
                case PivotalBootstrapInvalidDrawPolicy.Drop:
                    return false;
                case PivotalBootstrapInvalidDrawPolicy.UseRaw:
                    parameterSet = rawFit.Parameters.Clone();
                    return true;
                case PivotalBootstrapInvalidDrawPolicy.UseParent:
                    parameterSet = parentFit.Parameters.Clone();
                    return true;
                default:
                    throw new ArgumentOutOfRangeException(nameof(PivotalInvalidDrawPolicy), $"Unknown invalid-draw policy: {PivotalInvalidDrawPolicy}.");
            }
        }

        /// <summary>
        /// Computes the link-space covariance for a bootstrap fit using the supplied link controller.
        /// </summary>
        /// <param name="fit">The fit containing native parameter covariance.</param>
        /// <param name="linkController">The link controller used to compute the diagonal Jacobian.</param>
        /// <returns>The covariance matrix in link space.</returns>
        private Matrix LinkCovariance(BootstrapFit fit, LinkController linkController)
        {
            Matrix jacobian = linkController.LinkJacobian(fit.Parameters.Values);
            for (int i = 0; i < jacobian.NumberOfRows; i++)
            {
                if (!Tools.IsFinite(jacobian[i, i]))
                    throw new ArgumentException("The link derivative produced a non-finite value.");
            }

            Matrix linkCovariance = jacobian * fit.Covariance * jacobian.Transpose();
            return RegularizePivotalCovariances
                ? MatrixRegularization.MakeSymmetricPositiveDefinite(linkCovariance)
                : linkCovariance;
        }

        /// <summary>
        /// Builds the pivotal link array from <see cref="PivotalLinkFactory"/> or the identity default.
        /// </summary>
        /// <param name="parentFit">The original parent fit.</param>
        /// <param name="acceptedRawFits">The accepted raw bootstrap fits.</param>
        /// <returns>One optional link function per fitted parameter.</returns>
        /// <exception cref="ArgumentException">Thrown when the link factory returns the wrong number of links.</exception>
        private ILinkFunction?[] CreatePivotalLinks(BootstrapFit parentFit, BootstrapFit[] acceptedRawFits)
        {
            ILinkFunction?[] links = PivotalLinkFactory != null
                ? PivotalLinkFactory(new PivotalBootstrapContext(parentFit, acceptedRawFits))
                : new ILinkFunction?[parentFit.ParameterCount];

            if (links == null)
                throw new ArgumentException("The pivotal link factory must not return null.", nameof(PivotalLinkFactory));
            if (links.Length != parentFit.ParameterCount)
                throw new ArgumentException("The pivotal link factory must return one link per parameter.", nameof(PivotalLinkFactory));

            return (ILinkFunction?[])links.Clone();
        }

        #endregion

        #region Confidence Intervals

        /// <summary>
        /// Computes bootstrap confidence intervals using the specified method.
        /// </summary>
        /// <param name="method">The confidence interval method.</param>
        /// <param name="alpha">The two-sided alpha level. Default = 0.1, resulting in 90% confidence intervals.</param>
        /// <returns>A <see cref="BootstrapResults"/> object containing confidence intervals for parameters and statistics.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the requested interval method is incompatible with the last run mode.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="alpha"/> is not between zero and one.</exception>
        public BootstrapResults GetConfidenceIntervals(BootstrapCIMethod method, double alpha = 0.1)
        {
            ValidateConfidenceIntervalRequest(method, alpha);

            if (_runType == BootstrapRunType.Pivotal)
            {
                return CreatePercentileResults(
                    _bootstrapParameterSets,
                    _bootstrapStatistics,
                    _originalParameters.Values,
                    _parentStatistics,
                    alpha,
                    _failedCount,
                    method);
            }

            if (StatisticFunction == null)
                throw new InvalidOperationException("StatisticFunction must be set.");
            var originalStats = ValidateStatistics(StatisticFunction(_originalParameters));

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

            for (int i = 0; i < _numParams; i++)
            {
                var values = _bootstrapParameterSets.Select(ps => ps.Values[i]).ToArray();
                results.ParameterResults[i] = ComputePercentileCI(values, _originalParameters.Values[i], alpha);
            }

            return results;
        }

        /// <summary>
        /// Computes raw-bootstrap percentile confidence intervals after a pivotal bootstrap run.
        /// </summary>
        /// <param name="alpha">The two-sided alpha level. Default = 0.1, resulting in 90% confidence intervals.</param>
        /// <returns>Percentile confidence intervals for the accepted raw bootstrap ensemble.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the last run was not pivotal bootstrap.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="alpha"/> is not between zero and one.</exception>
        public BootstrapResults GetRawPivotalConfidenceIntervals(double alpha = 0.1)
        {
            if (_runType != BootstrapRunType.Pivotal)
                throw new InvalidOperationException("RunPivotalBootstrap() or TransformPivotalBootstrap() must be called before requesting raw pivotal-bootstrap confidence intervals.");
            if (alpha <= 0d || alpha >= 1d)
                throw new ArgumentOutOfRangeException(nameof(alpha), "Alpha must be between 0 and 1.");

            return CreatePercentileResults(
                _rawBootstrapParameterSets,
                _rawBootstrapStatistics,
                _originalParameters.Values,
                _parentStatistics,
                alpha,
                _failedCount,
                BootstrapCIMethod.Percentile);
        }

        #endregion

        #region CI Methods

        /// <summary>
        /// Computes percentile confidence intervals for a single statistic or parameter.
        /// </summary>
        /// <param name="values">Bootstrap values for one statistic or parameter.</param>
        /// <param name="populationEstimate">The original population estimate.</param>
        /// <param name="alpha">The two-sided alpha level.</param>
        /// <returns>The percentile confidence interval summary.</returns>
        private BootstrapStatisticResult ComputePercentileCI(double[] values, double populationEstimate, double alpha)
        {
            var validValues = values.Where(Tools.IsFinite).ToArray();
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
        /// Computes bias-corrected confidence intervals for a single statistic.
        /// </summary>
        /// <param name="values">Bootstrap statistic values.</param>
        /// <param name="populationEstimate">The original statistic estimate.</param>
        /// <param name="alpha">The two-sided alpha level.</param>
        /// <returns>The bias-corrected confidence interval summary.</returns>
        private BootstrapStatisticResult ComputeBiasCorrectedCI(double[] values, double populationEstimate, double alpha)
        {
            var validValues = values.Where(Tools.IsFinite).ToArray();
            int validN = validValues.Length;
            if (validN == 0) return EmptyResult(populationEstimate, values.Length);

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
        /// Computes bias-corrected and accelerated confidence intervals for a single statistic.
        /// </summary>
        /// <param name="values">Bootstrap statistic values.</param>
        /// <param name="populationEstimate">The original statistic estimate.</param>
        /// <param name="alpha">The two-sided alpha level.</param>
        /// <param name="acceleration">The jackknife acceleration constant.</param>
        /// <returns>The BCa confidence interval summary.</returns>
        private BootstrapStatisticResult ComputeBCaCI(double[] values, double populationEstimate, double alpha, double acceleration)
        {
            var validValues = values.Where(Tools.IsFinite).ToArray();
            int validN = validValues.Length;
            if (validN == 0) return EmptyResult(populationEstimate, values.Length);

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
        /// Computes Normal confidence intervals for a single statistic using the configured transform pair.
        /// </summary>
        /// <param name="values">Bootstrap statistic values.</param>
        /// <param name="populationEstimate">The original statistic estimate.</param>
        /// <param name="alpha">The two-sided alpha level.</param>
        /// <returns>The Normal confidence interval summary.</returns>
        private BootstrapStatisticResult ComputeNormalCI(double[] values, double populationEstimate, double alpha)
        {
            double popTransformed = ApplyTransform(populationEstimate);
            var transformedValid = new double[values.Length];
            int validCount = 0;
            for (int i = 0; i < values.Length; i++)
            {
                if (Tools.IsFinite(values[i]))
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
        /// Computes Bootstrap-t confidence intervals for a single statistic.
        /// </summary>
        /// <param name="statisticIndex">The zero-based statistic index.</param>
        /// <param name="populationEstimate">The original statistic estimate.</param>
        /// <param name="alpha">The two-sided alpha level.</param>
        /// <returns>The Bootstrap-t confidence interval summary.</returns>
        private BootstrapStatisticResult ComputeBootstrapTCI(int statisticIndex, double populationEstimate, double alpha)
        {
            double popTransformed = ApplyTransform(populationEstimate);
            var xCol = _transformedStatistics!.GetColumn(statisticIndex);
            var tCol = _studentizedValues!.GetColumn(statisticIndex);

            var validX = xCol.Where(Tools.IsFinite).ToArray();
            var validT = tCol.Where(Tools.IsFinite).ToArray();

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
        /// Computes acceleration constants for each statistic using leave-one-out jackknife samples.
        /// </summary>
        /// <param name="populationEstimates">The original statistic estimates.</param>
        /// <returns>The acceleration constant for each statistic.</returns>
        private double[] ComputeAccelerationConstants(double[] populationEstimates)
        {
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
                    // Skip failed jackknife samples.
                }
            });

            for (int i = 0; i < _numStats; i++)
                a[i] = I3[i] / (Math.Pow(I2[i], 1.5) * 6d);

            return a;
        }

        #endregion

        #region Private Helpers

        /// <summary>
        /// Validates that the regular bootstrap delegates are set.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when a required regular-bootstrap delegate is missing.</exception>
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
        /// Validates that the pivotal bootstrap delegates and parent covariance are set.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when a required pivotal-bootstrap input is missing.</exception>
        private void ValidatePivotalDelegates()
        {
            if (ResampleFunction == null)
                throw new InvalidOperationException("ResampleFunction must be set before running the pivotal bootstrap.");
            if (FitWithCovarianceFunction == null)
                throw new InvalidOperationException("FitWithCovarianceFunction must be set before running the pivotal bootstrap.");
            if (_originalCovariance == null)
                throw new InvalidOperationException("OriginalCovariance must be set before running the pivotal bootstrap.");
        }

        /// <summary>
        /// Validates the configured replicate and retry counts.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <see cref="Replicates"/> or <see cref="MaxRetries"/> is not positive.</exception>
        private void ValidateReplicationSettings()
        {
            if (Replicates < 1)
                throw new ArgumentOutOfRangeException(nameof(Replicates), "The number of replicates must be positive.");
            if (MaxRetries < 1)
                throw new ArgumentOutOfRangeException(nameof(MaxRetries), "The maximum retry count must be positive.");
        }

        /// <summary>
        /// Initializes regular-bootstrap state arrays before a regular run.
        /// </summary>
        private void InitializeState()
        {
            var originalStats = ValidateStatistics(StatisticFunction!(_originalParameters));
            _parentStatistics = (double[])originalStats.Clone();
            _numStats = originalStats.Length;
            _numParams = _originalParameters.Values.Length;

            _bootstrapParameterSets = new ParameterSet[Replicates];
            _bootstrapStatistics = new double[Replicates, _numStats];
            _validFlags = new bool[Replicates];
            _failedCount = 0;
            _studentizedValues = null;
            _transformedStatistics = null;
            _runType = BootstrapRunType.None;
        }

        /// <summary>
        /// Clears stored pivotal-bootstrap state before a regular bootstrap run.
        /// </summary>
        private void ResetPivotalState()
        {
            _rawBootstrapFits = Array.Empty<BootstrapFit>();
            _rawBootstrapParameterSets = Array.Empty<ParameterSet>();
            _rawBootstrapStatistics = null;
            _pivotalLinks = Array.Empty<ILinkFunction?>();
            _pivotalDiagnostics = null;
        }

        /// <summary>
        /// Marks a replicate as failed with NaN parameter and statistic values.
        /// </summary>
        /// <param name="idx">The zero-based replicate index.</param>
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
        /// Applies the configured regular-bootstrap statistic transform.
        /// </summary>
        /// <param name="value">The value to transform.</param>
        /// <returns>The transformed value, or the original value when <see cref="Transform"/> is null.</returns>
        private double ApplyTransform(double value)
        {
            return Transform != null ? Transform(value) : value;
        }

        /// <summary>
        /// Applies the configured inverse regular-bootstrap statistic transform.
        /// </summary>
        /// <param name="value">The value to inverse transform.</param>
        /// <returns>The inverse-transformed value, or the original value when <see cref="InverseTransform"/> is null.</returns>
        private double ApplyInverseTransform(double value)
        {
            return InverseTransform != null ? InverseTransform(value) : value;
        }

        /// <summary>
        /// Creates an empty confidence interval result for insufficient valid values.
        /// </summary>
        /// <param name="populationEstimate">The original estimate.</param>
        /// <param name="totalCount">The total number of attempted values.</param>
        /// <returns>A confidence interval result containing NaN interval bounds.</returns>
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

        /// <summary>
        /// Creates percentile results for parameter and statistic ensembles.
        /// </summary>
        /// <param name="parameterSets">The parameter ensemble.</param>
        /// <param name="statistics">The statistic ensemble, or null when no statistic function was supplied.</param>
        /// <param name="parameterEstimates">The original parameter estimates.</param>
        /// <param name="statisticEstimates">The original statistic estimates.</param>
        /// <param name="alpha">The two-sided alpha level.</param>
        /// <param name="failedReplicates">The number of failed replicates.</param>
        /// <param name="method">The confidence interval method label to store in the result.</param>
        /// <returns>A bootstrap results object containing percentile confidence intervals.</returns>
        private BootstrapResults CreatePercentileResults(
            ParameterSet[] parameterSets,
            double[,]? statistics,
            double[] parameterEstimates,
            double[]? statisticEstimates,
            double alpha,
            int failedReplicates,
            BootstrapCIMethod method)
        {
            int statisticCount = statistics?.GetLength(1) ?? 0;
            var results = new BootstrapResults
            {
                Method = method,
                Alpha = alpha,
                StatisticResults = new BootstrapStatisticResult[statisticCount],
                ParameterResults = new BootstrapStatisticResult[parameterEstimates.Length],
                FailedReplicates = failedReplicates
            };

            for (int i = 0; i < statisticCount; i++)
            {
                double estimate = statisticEstimates != null && i < statisticEstimates.Length ? statisticEstimates[i] : double.NaN;
                results.StatisticResults[i] = ComputePercentileCI(statistics!.GetColumn(i), estimate, alpha);
            }

            for (int i = 0; i < parameterEstimates.Length; i++)
                results.ParameterResults[i] = ComputePercentileCI(GetParameterColumn(parameterSets, i), parameterEstimates[i], alpha);

            return results;
        }

        /// <summary>
        /// Extracts one parameter column from an ensemble of parameter sets.
        /// </summary>
        /// <param name="parameterSets">The parameter ensemble.</param>
        /// <param name="parameterIndex">The zero-based parameter index.</param>
        /// <returns>The values for the requested parameter.</returns>
        private static double[] GetParameterColumn(ParameterSet[] parameterSets, int parameterIndex)
        {
            var values = new double[parameterSets.Length];
            for (int i = 0; i < parameterSets.Length; i++)
                values[i] = parameterSets[i].Values[parameterIndex];
            return values;
        }

        /// <summary>
        /// Computes statistic values for each parameter set.
        /// </summary>
        /// <param name="parameterSets">The parameter sets to evaluate.</param>
        /// <param name="statistic">The statistic function.</param>
        /// <returns>A two-dimensional statistic array indexed by replicate and statistic.</returns>
        private static double[,] ComputeStatistics(ParameterSet[] parameterSets, Func<ParameterSet, double[]> statistic)
        {
            if (parameterSets.Length == 0)
                return new double[0, 0];

            double[] first = ValidateStatistics(statistic(parameterSets[0]));
            var values = new double[parameterSets.Length, first.Length];
            for (int j = 0; j < first.Length; j++)
                values[0, j] = first[j];

            for (int i = 1; i < parameterSets.Length; i++)
            {
                double[] row = ValidateStatistics(statistic(parameterSets[i]), first.Length);
                for (int j = 0; j < first.Length; j++)
                    values[i, j] = row[j];
            }

            return values;
        }

        /// <summary>
        /// Validates statistic values returned by <see cref="StatisticFunction"/>.
        /// </summary>
        /// <param name="statistics">The statistic values to validate.</param>
        /// <param name="expectedLength">The expected statistic count, when known.</param>
        /// <returns>The validated statistic values.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the statistic values are null, empty, non-finite, or have the wrong length.</exception>
        private static double[] ValidateStatistics(double[] statistics, int? expectedLength = null)
        {
            if (statistics == null)
                throw new InvalidOperationException("The statistic function returned null.");
            if (statistics.Length == 0)
                throw new InvalidOperationException("The statistic function must return at least one statistic.");
            if (expectedLength.HasValue && statistics.Length != expectedLength.Value)
                throw new InvalidOperationException("The statistic function must return the same number of statistics for every draw.");
            if (!ContainsOnlyFiniteValues(statistics))
                throw new InvalidOperationException("The statistic function returned a non-finite value.");

            return statistics;
        }

        /// <summary>
        /// Determines whether a fitted parameter set has the expected finite parameter vector.
        /// </summary>
        /// <param name="parameters">The fitted parameter set to evaluate.</param>
        /// <param name="expectedLength">The expected number of parameters.</param>
        /// <returns>True when the parameter values are non-null, have the expected length, and are finite; otherwise false.</returns>
        private static bool HasExpectedFiniteParameterValues(ParameterSet parameters, int expectedLength)
        {
            return parameters.Values != null
                && parameters.Values.Length == expectedLength
                && ContainsOnlyFiniteValues(parameters.Values);
        }

        /// <summary>
        /// Creates the original covariance-aware fit for pivotal bootstrap.
        /// </summary>
        /// <returns>The original fit and covariance.</returns>
        /// <exception cref="InvalidOperationException">Thrown when <see cref="OriginalCovariance"/> is not set.</exception>
        private BootstrapFit CreateOriginalFit()
        {
            if (_originalCovariance == null)
                throw new InvalidOperationException("OriginalCovariance must be set before running or transforming a pivotal bootstrap.");

            var fit = new BootstrapFit(_originalParameters, _originalCovariance);
            if (!IsValidFit(fit, null))
                throw new ArgumentException("The original fit must contain finite parameters and a finite square covariance matrix.");

            return fit;
        }

        /// <summary>
        /// Determines whether a raw covariance-aware fit should be accepted before pivotal transformation.
        /// </summary>
        /// <param name="fit">The fit to evaluate.</param>
        /// <param name="parameterCount">The required parameter count.</param>
        /// <returns>True when the fit is valid and passes the optional replicate filter; otherwise false.</returns>
        private bool IsAcceptedRawFit(BootstrapFit? fit, int parameterCount)
        {
            if (!IsValidFit(fit, parameterCount))
                return false;
            if (PivotalReplicateFilter != null && !PivotalReplicateFilter(fit!))
                return false;
            return true;
        }

        /// <summary>
        /// Determines whether a covariance-aware fit has finite parameters and a finite square covariance matrix.
        /// </summary>
        /// <param name="fit">The fit to evaluate.</param>
        /// <param name="parameterCount">The required parameter count, when known.</param>
        /// <returns>True when the fit is valid; otherwise false.</returns>
        private static bool IsValidFit(BootstrapFit? fit, int? parameterCount)
        {
            if (fit == null)
                return false;
            if (fit.Parameters.Values == null || fit.Parameters.Values.Length == 0)
                return false;
            if (parameterCount.HasValue && fit.Parameters.Values.Length != parameterCount.Value)
                return false;
            if (!ContainsOnlyFiniteValues(fit.Parameters.Values))
                return false;
            if (fit.Covariance == null || !fit.Covariance.IsSquare || fit.Covariance.NumberOfRows != fit.Parameters.Values.Length)
                return false;

            for (int i = 0; i < fit.Covariance.NumberOfRows; i++)
            {
                for (int j = 0; j < fit.Covariance.NumberOfColumns; j++)
                {
                    if (!Tools.IsFinite(fit.Covariance[i, j]))
                        return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Validates that transformed values are finite.
        /// </summary>
        /// <param name="values">The transformed values.</param>
        /// <param name="message">The exception message to use when validation fails.</param>
        /// <exception cref="ArgumentException">Thrown when any value is non-finite.</exception>
        private static void ValidateTransformedValues(double[] values, string message)
        {
            if (!ContainsOnlyFiniteValues(values))
                throw new ArgumentException(message);
        }

        /// <summary>
        /// Determines whether every value in an array is finite.
        /// </summary>
        /// <param name="values">The values to evaluate.</param>
        /// <returns>True when every value is finite; otherwise false.</returns>
        private static bool ContainsOnlyFiniteValues(double[] values)
        {
            if (values == null)
                return false;

            for (int i = 0; i < values.Length; i++)
            {
                if (!Tools.IsFinite(values[i]))
                    return false;
            }

            return true;
        }

        /// <summary>
        /// Validates a confidence interval request against the last bootstrap run mode.
        /// </summary>
        /// <param name="method">The requested confidence interval method.</param>
        /// <param name="alpha">The requested alpha level.</param>
        /// <exception cref="InvalidOperationException">Thrown when no run has completed or when the method is incompatible with the run mode.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="alpha"/> is not between zero and one.</exception>
        private void ValidateConfidenceIntervalRequest(BootstrapCIMethod method, double alpha)
        {
            if (_runType == BootstrapRunType.None || _bootstrapStatistics == null)
                throw new InvalidOperationException("A bootstrap run must be completed before requesting confidence intervals.");
            if (alpha <= 0d || alpha >= 1d)
                throw new ArgumentOutOfRangeException(nameof(alpha), "Alpha must be between 0 and 1.");

            if (_runType == BootstrapRunType.Pivotal)
            {
                if (method != BootstrapCIMethod.Percentile)
                    throw new InvalidOperationException("Only percentile confidence intervals are supported after a pivotal bootstrap run.");
                return;
            }

            if (method == BootstrapCIMethod.BCa && (JackknifeFunction == null || SampleSizeFunction == null))
                throw new InvalidOperationException("JackknifeFunction and SampleSizeFunction must be set for BCa method.");
            if (method == BootstrapCIMethod.BootstrapT && _studentizedValues == null)
                throw new InvalidOperationException("RunWithStudentizedBootstrap() must be called before requesting Bootstrap-t confidence intervals.");
        }

        #endregion
    }
}
