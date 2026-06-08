using Numerics.Data.Statistics;
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
    /// Specifies how invalid pivotal bootstrap draws should be handled.
    /// </summary>
    public enum PivotalBootstrapInvalidDrawPolicy
    {
        /// <summary>
        /// Drop invalid pivotal draws from the returned ensemble.
        /// </summary>
        Drop,

        /// <summary>
        /// Replace invalid pivotal draws with the corresponding raw bootstrap fit.
        /// </summary>
        UseRaw,

        /// <summary>
        /// Replace invalid pivotal draws with the parent fit.
        /// </summary>
        UseParent
    }

    /// <summary>
    /// Stores a fitted parameter vector and its covariance matrix.
    /// </summary>
    [Serializable]
    public sealed class PivotalBootstrapFit
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="PivotalBootstrapFit"/> class.
        /// </summary>
        /// <param name="parameters">The fitted parameter set.</param>
        /// <param name="covariance">The covariance matrix for the fitted parameters.</param>
        public PivotalBootstrapFit(ParameterSet parameters, Matrix covariance)
        {
            if (parameters.Values == null)
                throw new ArgumentException("The parameter set must contain values.", nameof(parameters));
            if (covariance == null)
                throw new ArgumentNullException(nameof(covariance));
            if (!covariance.IsSquare)
                throw new ArgumentException("The covariance matrix must be square.", nameof(covariance));
            if (covariance.NumberOfRows != parameters.Values.Length)
                throw new ArgumentException("The covariance dimension must match the parameter count.", nameof(covariance));

            Parameters = parameters.Clone();
            Covariance = covariance.Clone();
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="PivotalBootstrapFit"/> class.
        /// </summary>
        /// <param name="parameters">The fitted parameter values.</param>
        /// <param name="covariance">The covariance matrix for the fitted parameters.</param>
        public PivotalBootstrapFit(double[] parameters, Matrix covariance)
            : this(new ParameterSet(parameters != null ? (double[])parameters.Clone() : throw new ArgumentNullException(nameof(parameters)), double.NaN), covariance)
        {
        }

        /// <summary>
        /// Gets the fitted parameter set.
        /// </summary>
        public ParameterSet Parameters { get; }

        /// <summary>
        /// Gets the covariance matrix for <see cref="Parameters"/>.
        /// </summary>
        public Matrix Covariance { get; }

        /// <summary>
        /// Gets the number of fitted parameters.
        /// </summary>
        public int ParameterCount => Parameters.Values.Length;
    }

    /// <summary>
    /// Delegate-backed scalar link used by <see cref="PivotalBootstrap"/>.
    /// </summary>
    /// <remarks>
    /// This type gives callers complete control over the link, inverse link, and
    /// derivative used for each parameter. Use <see cref="FromLinkFunction"/> to adapt
    /// Numerics link-function classes such as <see cref="LogLink"/>, <see cref="FisherZLink"/>,
    /// or <see cref="YeoJohnsonLink"/>.
    /// </remarks>
    [Serializable]
    public sealed class PivotalBootstrapLink
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="PivotalBootstrapLink"/> class.
        /// </summary>
        /// <param name="link">The link transformation.</param>
        /// <param name="inverseLink">The inverse link transformation.</param>
        /// <param name="derivative">The first derivative of the link with respect to the raw parameter.</param>
        /// <param name="name">An optional descriptive name.</param>
        public PivotalBootstrapLink(
            Func<double, double> link,
            Func<double, double> inverseLink,
            Func<double, double> derivative,
            string? name = null)
        {
            Link = link ?? throw new ArgumentNullException(nameof(link));
            InverseLink = inverseLink ?? throw new ArgumentNullException(nameof(inverseLink));
            Derivative = derivative ?? throw new ArgumentNullException(nameof(derivative));
            Name = name ?? "Custom";
        }

        /// <summary>
        /// Gets an identity link.
        /// </summary>
        public static PivotalBootstrapLink Identity { get; } = FromLinkFunction(new IdentityLink());

        /// <summary>
        /// Gets the link transformation.
        /// </summary>
        public Func<double, double> Link { get; }

        /// <summary>
        /// Gets the inverse link transformation.
        /// </summary>
        public Func<double, double> InverseLink { get; }

        /// <summary>
        /// Gets the first derivative of the link with respect to the raw parameter.
        /// </summary>
        public Func<double, double> Derivative { get; }

        /// <summary>
        /// Gets the descriptive link name.
        /// </summary>
        public string Name { get; }

        /// <summary>
        /// Adapts an <see cref="ILinkFunction"/> to a pivotal bootstrap link.
        /// </summary>
        /// <param name="linkFunction">The Numerics link function.</param>
        /// <param name="name">An optional descriptive name.</param>
        public static PivotalBootstrapLink FromLinkFunction(ILinkFunction linkFunction, string? name = null)
        {
            if (linkFunction == null)
                throw new ArgumentNullException(nameof(linkFunction));

            return new PivotalBootstrapLink(
                linkFunction.Link,
                linkFunction.InverseLink,
                linkFunction.DLink,
                name ?? linkFunction.GetType().Name);
        }
    }

    /// <summary>
    /// Context supplied to a pivotal bootstrap link factory.
    /// </summary>
    public sealed class PivotalBootstrapLinkContext
    {
        private readonly double[,] _rawParameterValues;

        /// <summary>
        /// Initializes a new instance of the <see cref="PivotalBootstrapLinkContext"/> class.
        /// </summary>
        /// <param name="parentFit">The parent fit.</param>
        /// <param name="rawBootstrapFits">The accepted raw bootstrap fits.</param>
        public PivotalBootstrapLinkContext(PivotalBootstrapFit parentFit, PivotalBootstrapFit[] rawBootstrapFits)
        {
            ParentFit = parentFit ?? throw new ArgumentNullException(nameof(parentFit));
            RawBootstrapFits = rawBootstrapFits ?? throw new ArgumentNullException(nameof(rawBootstrapFits));

            int p = parentFit.ParameterCount;
            _rawParameterValues = new double[rawBootstrapFits.Length, p];
            for (int i = 0; i < rawBootstrapFits.Length; i++)
            {
                for (int j = 0; j < p; j++)
                    _rawParameterValues[i, j] = rawBootstrapFits[i].Parameters.Values[j];
            }
        }

        /// <summary>
        /// Gets the parent fit.
        /// </summary>
        public PivotalBootstrapFit ParentFit { get; }

        /// <summary>
        /// Gets the accepted raw bootstrap fits.
        /// </summary>
        public PivotalBootstrapFit[] RawBootstrapFits { get; }

        /// <summary>
        /// Gets the number of parameters.
        /// </summary>
        public int ParameterCount => ParentFit.ParameterCount;

        /// <summary>
        /// Gets a copy of the raw bootstrap values for one parameter.
        /// </summary>
        /// <param name="parameterIndex">The zero-based parameter index.</param>
        public double[] GetRawParameterValues(int parameterIndex)
        {
            if (parameterIndex < 0 || parameterIndex >= ParameterCount)
                throw new ArgumentOutOfRangeException(nameof(parameterIndex));

            var values = new double[RawBootstrapFits.Length];
            for (int i = 0; i < values.Length; i++)
                values[i] = _rawParameterValues[i, parameterIndex];
            return values;
        }
    }

    /// <summary>
    /// Options for transforming raw bootstrap fits into pivotal bootstrap draws.
    /// </summary>
    public class PivotalBootstrapTransformOptions
    {
        /// <summary>
        /// Gets or sets a factory that returns one scalar link per parameter.
        /// </summary>
        /// <remarks>
        /// The factory is called after raw bootstrap fits have been accepted, so it can
        /// estimate data-adaptive links from the parent and raw bootstrap ensembles.
        /// If omitted, identity links are used for every parameter.
        /// </remarks>
        public Func<PivotalBootstrapLinkContext, PivotalBootstrapLink[]>? LinkFactory { get; set; }

        /// <summary>
        /// Gets or sets an optional statistic function applied to each returned parameter set.
        /// </summary>
        public Func<ParameterSet, double[]>? StatisticFunction { get; set; }

        /// <summary>
        /// Gets or sets an optional filter applied to raw bootstrap fits before the pivotal transform.
        /// </summary>
        public Func<PivotalBootstrapFit, bool>? ReplicateFilter { get; set; }

        /// <summary>
        /// Gets or sets an optional validator applied to transformed parameter values.
        /// </summary>
        public Func<double[], bool>? ParameterValidator { get; set; }

        /// <summary>
        /// Gets or sets the policy used when a pivotal draw is invalid. Default is <see cref="PivotalBootstrapInvalidDrawPolicy.Drop"/>.
        /// </summary>
        public PivotalBootstrapInvalidDrawPolicy InvalidDrawPolicy { get; set; } = PivotalBootstrapInvalidDrawPolicy.Drop;

        /// <summary>
        /// Gets or sets a value indicating whether covariance matrices should be repaired to be symmetric positive definite.
        /// </summary>
        public bool RegularizeCovariances { get; set; } = true;

        /// <summary>
        /// Gets or sets an optional absolute component limit for the standardized pivot vector.
        /// </summary>
        public double? ZLimit { get; set; }

        /// <summary>
        /// Gets or sets a value indicating whether small Gaussian jitter should be added in pivot space.
        /// </summary>
        public bool AddJitter { get; set; }

        /// <summary>
        /// Gets or sets the standard deviation of optional Gaussian jitter in pivot space.
        /// </summary>
        public double JitterScale { get; set; } = 0.01d;

        /// <summary>
        /// Gets or sets the pseudo-random-number generator seed. Default = 12345.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;
    }

    /// <summary>
    /// Options for running an end-to-end pivotal bootstrap.
    /// </summary>
    /// <typeparam name="TData">The data type being bootstrapped.</typeparam>
    public sealed class PivotalBootstrapOptions<TData> : PivotalBootstrapTransformOptions
    {
        /// <summary>
        /// Gets or sets the function used to resample data from the parent fit.
        /// </summary>
        public Func<TData, PivotalBootstrapFit, Random, TData>? ResampleFunction { get; set; }

        /// <summary>
        /// Gets or sets the function used to fit one bootstrap data set and compute its covariance.
        /// </summary>
        public Func<TData, PivotalBootstrapFit>? FitFunction { get; set; }

        /// <summary>
        /// Gets or sets the number of raw bootstrap replicates to request. Default = 10,000.
        /// </summary>
        public int Replicates { get; set; } = 10000;

        /// <summary>
        /// Gets or sets the maximum number of retries for a failed raw bootstrap fit. Default = 20.
        /// </summary>
        public int MaxRetries { get; set; } = 20;

        /// <summary>
        /// Gets or sets the maximum number of parallel workers. Values less than 1 use the default scheduler setting.
        /// </summary>
        public int MaxDegreeOfParallelism { get; set; } = 0;
    }

    /// <summary>
    /// Diagnostics returned by a pivotal bootstrap run.
    /// </summary>
    [Serializable]
    public sealed class PivotalBootstrapDiagnostics
    {
        /// <summary>
        /// Gets or sets the number of raw bootstrap replicates requested.
        /// </summary>
        public int RequestedReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of raw bootstrap fits rejected before transformation.
        /// </summary>
        public int RejectedRawReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of raw bootstrap fits that failed during end-to-end resampling.
        /// </summary>
        public int FailedRawReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of raw bootstrap fits accepted for transformation.
        /// </summary>
        public int AcceptedRawReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of transformed pivotal draws retained.
        /// </summary>
        public int RetainedPivotalReplicates { get; set; }

        /// <summary>
        /// Gets or sets the number of transformed pivotal draws dropped or replaced according to policy.
        /// </summary>
        public int InvalidPivotalReplicates { get; set; }

        /// <summary>
        /// Gets or sets the raw bootstrap resampling time.
        /// </summary>
        public TimeSpan ResamplingTime { get; set; }

        /// <summary>
        /// Gets or sets the pivotal transformation time.
        /// </summary>
        public TimeSpan TransformationTime { get; set; }
    }

    /// <summary>
    /// Results returned by a pivotal bootstrap run.
    /// </summary>
    [Serializable]
    public sealed class PivotalBootstrapResults
    {
        internal PivotalBootstrapResults(
            PivotalBootstrapFit parentFit,
            PivotalBootstrapFit[] rawFits,
            ParameterSet[] pivotalParameterSets,
            PivotalBootstrapLink[] links,
            double[]? parentStatistics,
            double[,]? rawStatistics,
            double[,]? pivotalStatistics,
            PivotalBootstrapDiagnostics diagnostics)
        {
            ParentFit = parentFit;
            RawFits = rawFits;
            RawParameterSets = rawFits.Select(f => f.Parameters.Clone()).ToArray();
            PivotalParameterSets = pivotalParameterSets;
            Links = links;
            ParentStatistics = parentStatistics;
            RawStatistics = rawStatistics;
            PivotalStatistics = pivotalStatistics;
            Diagnostics = diagnostics;
        }

        /// <summary>
        /// Gets the parent fit.
        /// </summary>
        public PivotalBootstrapFit ParentFit { get; }

        /// <summary>
        /// Gets the accepted raw bootstrap fits and their covariances.
        /// </summary>
        public PivotalBootstrapFit[] RawFits { get; }

        /// <summary>
        /// Gets the accepted raw bootstrap parameter sets.
        /// </summary>
        public ParameterSet[] RawParameterSets { get; }

        /// <summary>
        /// Gets the retained pivotal bootstrap parameter sets.
        /// </summary>
        public ParameterSet[] PivotalParameterSets { get; }

        /// <summary>
        /// Gets the links used for the pivotal transformation.
        /// </summary>
        public PivotalBootstrapLink[] Links { get; }

        /// <summary>
        /// Gets the statistic values for the parent fit, when a statistic function was supplied.
        /// </summary>
        public double[]? ParentStatistics { get; }

        /// <summary>
        /// Gets raw bootstrap statistics as [replicate, statistic], when a statistic function was supplied.
        /// </summary>
        public double[,]? RawStatistics { get; }

        /// <summary>
        /// Gets pivotal bootstrap statistics as [replicate, statistic], when a statistic function was supplied.
        /// </summary>
        public double[,]? PivotalStatistics { get; }

        /// <summary>
        /// Gets run diagnostics.
        /// </summary>
        public PivotalBootstrapDiagnostics Diagnostics { get; }

        /// <summary>
        /// Gets percentile confidence intervals for pivotal bootstrap parameters.
        /// </summary>
        /// <param name="alpha">The two-sided alpha level, e.g. 0.1 for a 90 percent interval.</param>
        public BootstrapStatisticResult[] GetPivotalParameterConfidenceIntervals(double alpha)
        {
            return CreateIntervalResults(PivotalParameterSets, ParentFit.Parameters.Values, alpha);
        }

        /// <summary>
        /// Gets percentile confidence intervals for raw bootstrap parameters.
        /// </summary>
        /// <param name="alpha">The two-sided alpha level, e.g. 0.1 for a 90 percent interval.</param>
        public BootstrapStatisticResult[] GetRawParameterConfidenceIntervals(double alpha)
        {
            return CreateIntervalResults(RawParameterSets, ParentFit.Parameters.Values, alpha);
        }

        /// <summary>
        /// Gets percentile confidence intervals for pivotal bootstrap statistics.
        /// </summary>
        /// <param name="alpha">The two-sided alpha level, e.g. 0.1 for a 90 percent interval.</param>
        public BootstrapStatisticResult[] GetPivotalStatisticConfidenceIntervals(double alpha)
        {
            if (PivotalStatistics == null || ParentStatistics == null)
                throw new InvalidOperationException("No statistic function was supplied.");

            return CreateIntervalResults(PivotalStatistics, ParentStatistics, alpha);
        }

        /// <summary>
        /// Gets percentile confidence intervals for raw bootstrap statistics.
        /// </summary>
        /// <param name="alpha">The two-sided alpha level, e.g. 0.1 for a 90 percent interval.</param>
        public BootstrapStatisticResult[] GetRawStatisticConfidenceIntervals(double alpha)
        {
            if (RawStatistics == null || ParentStatistics == null)
                throw new InvalidOperationException("No statistic function was supplied.");

            return CreateIntervalResults(RawStatistics, ParentStatistics, alpha);
        }

        private static BootstrapStatisticResult[] CreateIntervalResults(ParameterSet[] parameterSets, double[] estimates, double alpha)
        {
            if (parameterSets.Length == 0)
                return Array.Empty<BootstrapStatisticResult>();

            int columns = parameterSets[0].Values.Length;
            var values = new double[parameterSets.Length, columns];
            for (int i = 0; i < parameterSets.Length; i++)
            {
                for (int j = 0; j < columns; j++)
                    values[i, j] = parameterSets[i].Values[j];
            }

            return CreateIntervalResults(values, estimates, alpha);
        }

        private static BootstrapStatisticResult[] CreateIntervalResults(double[,] values, double[] estimates, double alpha)
        {
            if (alpha <= 0d || alpha >= 1d)
                throw new ArgumentOutOfRangeException(nameof(alpha), "Alpha must be between 0 and 1.");

            int rows = values.GetLength(0);
            int columns = values.GetLength(1);
            var results = new BootstrapStatisticResult[columns];
            for (int j = 0; j < columns; j++)
            {
                var column = new List<double>(rows);
                for (int i = 0; i < rows; i++)
                {
                    double value = values[i, j];
                    if (IsFinite(value))
                        column.Add(value);
                }

                column.Sort();
                results[j] = new BootstrapStatisticResult
                {
                    PopulationEstimate = estimates[j],
                    LowerCI = column.Count > 0 ? Statistics.Percentile(column, alpha / 2d, true) : double.NaN,
                    UpperCI = column.Count > 0 ? Statistics.Percentile(column, 1d - alpha / 2d, true) : double.NaN,
                    ValidCount = column.Count,
                    TotalCount = rows,
                    StandardError = column.Count > 1 ? Statistics.StandardDeviation(column) : double.NaN,
                    Mean = column.Count > 0 ? column.Average() : double.NaN
                };
            }

            return results;
        }

        private static bool IsFinite(double value)
        {
            return !double.IsNaN(value) && !double.IsInfinity(value);
        }
    }

    /// <summary>
    /// Runs the bias-corrected pivotal bootstrap.
    /// </summary>
    public static class PivotalBootstrap
    {
        /// <summary>
        /// Runs an end-to-end pivotal bootstrap.
        /// </summary>
        /// <typeparam name="TData">The data type being bootstrapped.</typeparam>
        /// <param name="originalData">The original data.</param>
        /// <param name="parentFit">The parent fit and covariance.</param>
        /// <param name="options">Bootstrap options.</param>
        public static PivotalBootstrapResults Run<TData>(
            TData originalData,
            PivotalBootstrapFit parentFit,
            PivotalBootstrapOptions<TData> options)
        {
            if (options == null)
                throw new ArgumentNullException(nameof(options));
            if (options.ResampleFunction == null)
                throw new InvalidOperationException("A resample function is required.");
            if (options.FitFunction == null)
                throw new InvalidOperationException("A fit function is required.");
            if (options.Replicates < 1)
                throw new ArgumentOutOfRangeException(nameof(options.Replicates), "The number of replicates must be positive.");
            if (options.MaxRetries < 1)
                throw new ArgumentOutOfRangeException(nameof(options.MaxRetries), "The maximum retry count must be positive.");

            ValidateFit(parentFit, null);
            var stopwatch = Stopwatch.StartNew();
            var prng = new MersenneTwister(options.PRNGSeed);
            int[] seeds = prng.NextIntegers(options.Replicates);
            var rawFits = new PivotalBootstrapFit?[options.Replicates];
            int failed = 0;

            var parallelOptions = new ParallelOptions();
            if (options.MaxDegreeOfParallelism > 0)
                parallelOptions.MaxDegreeOfParallelism = options.MaxDegreeOfParallelism;

            Parallel.For(0, options.Replicates, parallelOptions, index =>
            {
                bool accepted = false;
                for (int retry = 0; retry < options.MaxRetries; retry++)
                {
                    try
                    {
                        var rng = new MersenneTwister(seeds[index] + 10 * retry);
                        TData sample = options.ResampleFunction(originalData, parentFit, rng);
                        PivotalBootstrapFit fit = options.FitFunction(sample);
                        if (!IsAcceptedRawFit(fit, parentFit.ParameterCount, options))
                            continue;

                        rawFits[index] = fit;
                        accepted = true;
                    }
                    catch
                    {
                        // Retry failed resamples or fits.
                    }

                    if (accepted)
                        break;
                }

                if (!accepted)
                    Interlocked.Increment(ref failed);
            });

            stopwatch.Stop();
            PivotalBootstrapFit[] acceptedFits = rawFits.Where(f => f != null).Select(f => f!).ToArray();
            PivotalBootstrapResults results = Transform(parentFit, acceptedFits, options);
            results.Diagnostics.RequestedReplicates = options.Replicates;
            results.Diagnostics.FailedRawReplicates = failed;
            results.Diagnostics.ResamplingTime = stopwatch.Elapsed;
            return results;
        }

        /// <summary>
        /// Transforms raw bootstrap fits into pivotal bootstrap draws.
        /// </summary>
        /// <param name="parentFit">The parent fit and covariance.</param>
        /// <param name="rawBootstrapFits">Raw bootstrap fits and covariances.</param>
        /// <param name="options">Transformation options. If omitted, identity links are used.</param>
        public static PivotalBootstrapResults Transform(
            PivotalBootstrapFit parentFit,
            IEnumerable<PivotalBootstrapFit> rawBootstrapFits,
            PivotalBootstrapTransformOptions? options = null)
        {
            if (rawBootstrapFits == null)
                throw new ArgumentNullException(nameof(rawBootstrapFits));

            options ??= new PivotalBootstrapTransformOptions();
            ValidateFit(parentFit, null);
            int p = parentFit.ParameterCount;
            var stopwatch = Stopwatch.StartNew();
            PivotalBootstrapFit[] suppliedRaw = rawBootstrapFits.ToArray();
            var acceptedRaw = suppliedRaw
                .Where(f => IsAcceptedRawFit(f, p, options))
                .ToArray();

            if (acceptedRaw.Length == 0)
                throw new InvalidOperationException("No valid raw bootstrap fits were supplied.");

            var diagnostics = new PivotalBootstrapDiagnostics
            {
                RequestedReplicates = acceptedRaw.Length,
                RejectedRawReplicates = suppliedRaw.Length - acceptedRaw.Length,
                AcceptedRawReplicates = acceptedRaw.Length
            };

            var context = new PivotalBootstrapLinkContext(parentFit, acceptedRaw);
            PivotalBootstrapLink[] links = options.LinkFactory != null
                ? options.LinkFactory(context)
                : Enumerable.Repeat(PivotalBootstrapLink.Identity, p).ToArray();
            ValidateLinks(links, p);

            double[] parentEta = LinkParameters(parentFit.Parameters.Values, links);
            Matrix parentLinkCovariance = LinkCovariance(parentFit, links, options.RegularizeCovariances);
            var parentCholesky = new CholeskyDecomposition(parentLinkCovariance);

            var pivotalParameterSets = new List<ParameterSet>(acceptedRaw.Length);
            var jitterRng = new MersenneTwister(options.PRNGSeed);
            int invalid = 0;

            for (int i = 0; i < acceptedRaw.Length; i++)
            {
                PivotalBootstrapFit rawFit = acceptedRaw[i];
                if (TryCreatePivotalDraw(rawFit, parentFit, links, parentEta, parentCholesky, options, jitterRng, out ParameterSet pivotalParameters))
                {
                    pivotalParameterSets.Add(pivotalParameters);
                }
                else
                {
                    invalid++;
                    if (TryApplyInvalidPolicy(rawFit, parentFit, options.InvalidDrawPolicy, out ParameterSet fallbackParameters))
                        pivotalParameterSets.Add(fallbackParameters);
                }
            }

            stopwatch.Stop();
            diagnostics.InvalidPivotalReplicates = invalid;
            diagnostics.RetainedPivotalReplicates = pivotalParameterSets.Count;
            diagnostics.TransformationTime = stopwatch.Elapsed;

            Func<ParameterSet, double[]>? statistic = options.StatisticFunction;
            double[]? parentStatistics = statistic != null ? ValidateStatistics(statistic(parentFit.Parameters)) : null;
            double[,]? rawStatistics = statistic != null ? ComputeStatistics(acceptedRaw.Select(f => f.Parameters).ToArray(), statistic) : null;
            double[,]? pivotalStatistics = statistic != null ? ComputeStatistics(pivotalParameterSets.ToArray(), statistic) : null;

            return new PivotalBootstrapResults(
                parentFit,
                acceptedRaw,
                pivotalParameterSets.ToArray(),
                links,
                parentStatistics,
                rawStatistics,
                pivotalStatistics,
                diagnostics);
        }

        private static bool TryCreatePivotalDraw(
            PivotalBootstrapFit rawFit,
            PivotalBootstrapFit parentFit,
            PivotalBootstrapLink[] links,
            double[] parentEta,
            CholeskyDecomposition parentCholesky,
            PivotalBootstrapTransformOptions options,
            Random jitterRng,
            out ParameterSet pivotalParameters)
        {
            pivotalParameters = default;
            try
            {
                double[] rawEta = LinkParameters(rawFit.Parameters.Values, links);
                Matrix rawLinkCovariance = LinkCovariance(rawFit, links, options.RegularizeCovariances);
                var rawCholesky = new CholeskyDecomposition(rawLinkCovariance);
                var difference = new double[parentEta.Length];
                for (int j = 0; j < difference.Length; j++)
                    difference[j] = parentEta[j] - rawEta[j];

                double[] z = rawCholesky.Forward(new Vector(difference)).ToArray();
                if (options.ZLimit.HasValue && z.Any(value => Math.Abs(value) > options.ZLimit.Value))
                    return false;

                if (options.AddJitter && options.JitterScale > 0d)
                {
                    for (int j = 0; j < z.Length; j++)
                        z[j] += options.JitterScale * Numerics.Distributions.Normal.StandardZ(Math.Min(1d - 1e-16, Math.Max(1e-16, jitterRng.NextDouble())));
                }

                double[] reinflated = parentCholesky.L * z;
                double[] pivotalEta = new double[parentEta.Length];
                double[] pivotalValues = new double[parentEta.Length];
                for (int j = 0; j < pivotalValues.Length; j++)
                {
                    pivotalEta[j] = parentEta[j] + reinflated[j];
                    pivotalValues[j] = links[j].InverseLink(pivotalEta[j]);
                }

                if (!AllFinite(pivotalValues))
                    return false;
                if (options.ParameterValidator != null && !options.ParameterValidator(pivotalValues))
                    return false;

                pivotalParameters = new ParameterSet(pivotalValues, rawFit.Parameters.Fitness, rawFit.Parameters.Weight);
                return true;
            }
            catch
            {
                return false;
            }
        }

        private static bool TryApplyInvalidPolicy(
            PivotalBootstrapFit rawFit,
            PivotalBootstrapFit parentFit,
            PivotalBootstrapInvalidDrawPolicy policy,
            out ParameterSet parameterSet)
        {
            parameterSet = default;
            switch (policy)
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
                    throw new ArgumentOutOfRangeException(nameof(policy), $"Unknown invalid-draw policy: {policy}.");
            }
        }

        private static bool IsAcceptedRawFit(PivotalBootstrapFit? fit, int parameterCount, PivotalBootstrapTransformOptions options)
        {
            if (fit == null || !IsValidFit(fit, parameterCount))
                return false;
            if (options.ReplicateFilter != null && !options.ReplicateFilter(fit))
                return false;
            return true;
        }

        private static void ValidateFit(PivotalBootstrapFit fit, int? parameterCount)
        {
            if (!IsValidFit(fit, parameterCount))
                throw new ArgumentException("The fit must contain finite parameters and a finite square covariance matrix.", nameof(fit));
        }

        private static bool IsValidFit(PivotalBootstrapFit? fit, int? parameterCount)
        {
            if (fit == null)
                return false;
            if (fit.Parameters.Values == null || fit.Parameters.Values.Length == 0)
                return false;
            if (parameterCount.HasValue && fit.Parameters.Values.Length != parameterCount.Value)
                return false;
            if (!AllFinite(fit.Parameters.Values))
                return false;
            if (fit.Covariance == null || !fit.Covariance.IsSquare || fit.Covariance.NumberOfRows != fit.Parameters.Values.Length)
                return false;

            for (int i = 0; i < fit.Covariance.NumberOfRows; i++)
            {
                for (int j = 0; j < fit.Covariance.NumberOfColumns; j++)
                {
                    if (!IsFinite(fit.Covariance[i, j]))
                        return false;
                }
            }

            return true;
        }

        private static void ValidateLinks(PivotalBootstrapLink[] links, int parameterCount)
        {
            if (links == null)
                throw new ArgumentNullException(nameof(links));
            if (links.Length != parameterCount)
                throw new ArgumentException("The link factory must return one link per parameter.", nameof(links));
            if (links.Any(link => link == null))
                throw new ArgumentException("The link factory returned a null link.", nameof(links));
        }

        private static double[] LinkParameters(double[] parameters, PivotalBootstrapLink[] links)
        {
            var eta = new double[parameters.Length];
            for (int i = 0; i < parameters.Length; i++)
                eta[i] = links[i].Link(parameters[i]);

            if (!AllFinite(eta))
                throw new ArgumentException("The link transformation produced a non-finite value.");

            return eta;
        }

        private static Matrix LinkCovariance(PivotalBootstrapFit fit, PivotalBootstrapLink[] links, bool regularize)
        {
            int p = fit.ParameterCount;
            var jacobian = new Matrix(p);
            for (int i = 0; i < p; i++)
            {
                double derivative = links[i].Derivative(fit.Parameters.Values[i]);
                if (!IsFinite(derivative))
                    throw new ArgumentException("The link derivative produced a non-finite value.");
                jacobian[i, i] = derivative;
            }

            Matrix linkCovariance = jacobian * fit.Covariance * jacobian.Transpose();
            if (!regularize)
                return linkCovariance;

            return MatrixRegularization.MakeSymmetricPositiveDefinite(linkCovariance);
        }

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

        private static double[] ValidateStatistics(double[] statistics, int? expectedLength = null)
        {
            if (statistics == null)
                throw new InvalidOperationException("The statistic function returned null.");
            if (statistics.Length == 0)
                throw new InvalidOperationException("The statistic function must return at least one statistic.");
            if (expectedLength.HasValue && statistics.Length != expectedLength.Value)
                throw new InvalidOperationException("The statistic function must return the same number of statistics for every draw.");
            if (!AllFinite(statistics))
                throw new InvalidOperationException("The statistic function returned a non-finite value.");

            return statistics;
        }

        private static bool AllFinite(double[] values)
        {
            for (int i = 0; i < values.Length; i++)
            {
                if (!IsFinite(values[i]))
                    return false;
            }

            return true;
        }

        private static bool IsFinite(double value)
        {
            return !double.IsNaN(value) && !double.IsInfinity(value);
        }

    }
}
