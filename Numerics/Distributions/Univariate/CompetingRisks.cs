using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Threading;
using System.Xml.Linq;

namespace Numerics.Distributions
{
    /// <summary>
    /// A competing risks distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://reliability.readthedocs.io/en/latest/Competing%20risk%20models.html" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class CompetingRisks : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {
        /// <summary>
        /// Construct new competing risks distribution.
        /// </summary>
        /// <param name="distributions">The competing distributions.</param>
        public CompetingRisks(UnivariateDistributionBase[] distributions)
        {
            SetParameters(distributions);
        }

        /// <summary>
        /// Construct new competing risks distribution.
        /// </summary>
        /// <param name="distributions">The competing distributions.</param>
        public CompetingRisks(IUnivariateDistribution[] distributions)
        {
            SetParameters(distributions);
        }

        private UnivariateDistributionBase[] _distributions = null!;
        private EmpiricalDistribution _empiricalCDF = null!;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;
        private bool _empiricalCDFCreated = false;
        private double[,] _correlationMatrix = null!;
        private bool _mvnCreated = false;
        private Probability.DependencyType _dependency = Probability.DependencyType.Independent;
        private MultivariateNormal _mvn = null!;
        private int _prngSeed = MultivariateNormal.DefaultMVNUNISeed;

        private string? _cachedConfiguration;
        [NonSerialized] private WeibullConfiguration? _weibullConfiguration;

        /// <summary>Immutable scalar state whose equality implies identical built-in Weibull configuration XML.</summary>
        private sealed class WeibullConfiguration
        {
            private readonly bool _minimum;
            private readonly Probability.DependencyType _dependency;
            private readonly int _seed;
            private readonly Transform _xTransform;
            private readonly Transform _probabilityTransform;
            private readonly long[] _parameterBits;
            private readonly long[]? _matrixBits;
            private readonly int _matrixRows;
            private readonly int _matrixColumns;
            private DensityStep? _densityStep;

            /// <summary>Publishes the lazily computed step atomically for concurrent readers of unchanged state.</summary>
            private sealed class DensityStep
            {
                internal readonly double Value;
                internal DensityStep(double value) { Value = value; }
            }

            /// <summary>Returns a previously computed component-derived step, excluding observation-dependent fallbacks.</summary>
            internal bool TryGetDensityStep(out double step)
            {
                var cached = Volatile.Read(ref _densityStep);
                step = cached is null ? 0d : cached.Value;
                return cached is not null;
            }

            /// <summary>Stores the unchanged derivative-step expression for this exact component configuration.</summary>
            internal void CacheDensityStep(double step) => Volatile.Write(ref _densityStep, new DensityStep(step));

            private WeibullConfiguration(CompetingRisks owner)
            {
                _minimum = owner.MinimumOfRandomVariables;
                _dependency = owner.Dependency;
                _seed = owner.PRNGSeed;
                _xTransform = owner.XTransform;
                _probabilityTransform = owner.ProbabilityTransform;
                _parameterBits = new long[2 * owner._distributions.Length];
                for (int i = 0; i < owner._distributions.Length; i++)
                {
                    var weibull = (Weibull)owner._distributions[i];
                    _parameterBits[2 * i] = BitConverter.DoubleToInt64Bits(weibull.Lambda);
                    _parameterBits[2 * i + 1] = BitConverter.DoubleToInt64Bits(weibull.Kappa);
                }
                var matrix = owner._correlationMatrix;
                if (matrix is not null)
                {
                    _matrixRows = matrix.GetLength(0);
                    _matrixColumns = matrix.GetLength(1);
                    _matrixBits = new long[matrix.Length];
                    int rowStart = matrix.GetLowerBound(0), columnStart = matrix.GetLowerBound(1), index = 0;
                    for (int row = 0; row < _matrixRows; row++)
                        for (int column = 0; column < _matrixColumns; column++)
                            _matrixBits[index++] = BitConverter.DoubleToInt64Bits(matrix[rowStart + row, columnStart + column]);
                }
            }

            /// <summary>Captures only exact built-in Weibulls; derived and custom XML callbacks retain the generic path.</summary>
            internal static WeibullConfiguration? Capture(CompetingRisks owner)
            {
                if (owner._distributions is null) return null;
                foreach (var distribution in owner._distributions)
                    if (distribution is null || distribution.GetType() != typeof(Weibull)) return null;
                return new WeibullConfiguration(owner);
            }

            /// <summary>Compares live values without allocating wrappers, parameter arrays, or XML.</summary>
            internal bool Matches(CompetingRisks owner)
            {
                if (owner._distributions is null || owner._distributions.Length != _parameterBits.Length / 2
                    || owner.MinimumOfRandomVariables != _minimum || owner.Dependency != _dependency
                    || owner.PRNGSeed != _seed || owner.XTransform != _xTransform
                    || owner.ProbabilityTransform != _probabilityTransform) return false;
                for (int i = 0; i < owner._distributions.Length; i++)
                {
                    var distribution = owner._distributions[i];
                    if (distribution is null || distribution.GetType() != typeof(Weibull)) return false;
                    var weibull = (Weibull)distribution;
                    if (BitConverter.DoubleToInt64Bits(weibull.Lambda) != _parameterBits[2 * i]
                        || BitConverter.DoubleToInt64Bits(weibull.Kappa) != _parameterBits[2 * i + 1]) return false;
                }
                var matrix = owner._correlationMatrix;
                if (matrix is null) return _matrixBits is null;
                if (_matrixBits is null || matrix.GetLength(0) != _matrixRows || matrix.GetLength(1) != _matrixColumns) return false;
                int rowStart = matrix.GetLowerBound(0), columnStart = matrix.GetLowerBound(1), index = 0;
                for (int row = 0; row < _matrixRows; row++)
                    for (int column = 0; column < _matrixColumns; column++)
                        if (BitConverter.DoubleToInt64Bits(matrix[rowStart + row, columnStart + column]) != _matrixBits[index++]) return false;
                return true;
            }
        }

        /// <summary>Invalidates derived caches when mutable components or configuration change.</summary>
        private void RefreshCachedConfiguration()
        {
            var previous = Volatile.Read(ref _weibullConfiguration);
            if (previous is not null && previous.Matches(this)) return;
            // Capture before canonical serialization: fallback callbacks may mutate their configuration.
            var next = WeibullConfiguration.Capture(this);
            string configuration = DistributionNumerics.ConfigurationState(this);
            if (configuration != _cachedConfiguration)
            {
                _cachedConfiguration = configuration;
                _momentsComputed = false;
                _empiricalCDFCreated = false;
                _mvnCreated = false;
            }
            Volatile.Write(ref _weibullConfiguration, next);
        }

        /// <summary>
        /// Returns the array of univariate probability distributions.
        /// </summary>
        public ReadOnlyCollection<UnivariateDistributionBase> Distributions => new(_distributions);

        /// <summary>
        /// The seed for the multivariate normal's quadrature randomizer, used by the dependent
        /// (perfectly negative and correlation-matrix) branches.
        /// </summary>
        /// <remarks>
        /// Those branches evaluate a randomized-lattice rectangle integral drawing from
        /// <see cref="MultivariateNormal.MVNUNI"/>, so their results reproduce only when the seed
        /// is fixed. Applied when the multivariate normal is built — set it before the first
        /// dependent evaluation.
        /// </remarks>
        public int PRNGSeed
        {
            get { return _prngSeed; }
            set
            {
                if (_prngSeed == value) return;
                _prngSeed = value;
                _mvnCreated = false;
                _empiricalCDFCreated = false;
            }
        }

        /// <summary>
        /// Determines the interpolation transform for the X-values.
        /// </summary>
        public Transform XTransform { get; set; } = Transform.None;

        /// <summary>
        /// Determines the interpolation transform for the Probability-values.
        /// </summary>
        public Transform ProbabilityTransform { get; set; } = Transform.NormalZ;

        /// <summary>
        /// If true, the competing risks model computes the minimum of the random variables. If false, it computes the maximum of random variables. 
        /// </summary>
        public bool MinimumOfRandomVariables { get; set; } = true;

        /// <summary>
        /// The dependency between random variables. 
        /// </summary>
        public Probability.DependencyType Dependency
        {
            get { return _dependency; }
            set
            {
                if (_dependency != value)
                {
                    _dependency = value;
                    _mvnCreated = false;
                }
            }
        }

        /// <summary>
        /// The correlation matrix used for modeling dependency between the marginal distributions.
        /// This is only used when the Dependency Type = CorrelationMatrix.
        /// </summary>
        public double[,] CorrelationMatrix 
        { 
            get {  return _correlationMatrix; } 
            set
            {
                _correlationMatrix = value;
                _mvnCreated = false;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get
            {
                int sum = 0;
                for (int i = 0; i < Distributions.Count; i++)
                    sum += Distributions[i].NumberOfParameters;
                return sum;
            }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type => UnivariateDistributionType.CompetingRisks;

        /// <inheritdoc/>
        public override string DisplayName => "Competing Risks";

        /// <inheritdoc/>
        public override string ShortDisplayName => "CR";

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[1, 2];
                string Dstring = "{";
                for (int i = 0; i < Distributions.Count; i++)
                {
                    Dstring += Distributions[i].DisplayName;
                    if (i < Distributions.Count - 1)
                    {
                        Dstring += ",";
                    }
                }
                Dstring += "}";
                parmString[0, 0] = "Distributions";
                parmString[0, 1] = Dstring;
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNames
        {
            get
            {
                var result = new List<string>();
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    for (int j = 0; j < Distributions[i].ParameterNames.Length; j++)
                    {
                        result.Add("D" + (i + 1).ToString() + " " + Distributions[i].ParameterNames[j]);
                    }
                }
                return result.ToArray();
            }
        }


        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get
            {
                var result = new List<string>();
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    for (int j = 0; j < Distributions[i].ParameterNamesShortForm.Length; j++)
                    {
                        result.Add("D" + (i + 1).ToString() + " " + Distributions[i].ParameterNamesShortForm[j]);
                    }
                }
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override double[] GetParameters 
        {
            get
            {
                var result = new List<double>();
                for (int i = 0; i < Distributions.Count; i++)
                    result.AddRange(Distributions[i].GetParameters);
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Distributions)]; }
        }

        /// <summary>
        /// Compute central moments of the distribution.
        /// </summary>
        private void ComputeMoments()
        {
            double center = InverseCDF(.5), scale = InverseCDF(.75) - InverseCDF(.25);
            var mom = DistributionMomentIntegration.Compute(LogPDF, Minimum, Maximum, center, scale);
            u1 = mom[0];
            u2 = mom[1];
            u3 = mom[2];
            u4 = mom[3];
            _momentsComputed = true;
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                RefreshCachedConfiguration();
                if (!_momentsComputed) ComputeMoments();
                return u1;
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get 
            {
                var brent = new BrentSearch(PDF, InverseCDF(0.001), InverseCDF(0.999));
                brent.Maximize();
                return brent.BestParameterSet.Values[0];
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                RefreshCachedConfiguration();
                if (!_momentsComputed) ComputeMoments();
                return u2;
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                RefreshCachedConfiguration();
                if (!_momentsComputed) ComputeMoments();
                return u3;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                RefreshCachedConfiguration();
                if (!_momentsComputed) ComputeMoments();
                return u4;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return MinimumOfRandomVariables ? _distributions.Min(p => p.Minimum) : _distributions.Max(p => p.Minimum); }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return MinimumOfRandomVariables ? _distributions.Min(p => p.Maximum) : _distributions.Max(p => p.Maximum); }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get
            {
                var result = new List<double>();
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.AddRange(Distributions[i].MinimumOfParameters);
                }
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get
            {
                var result = new List<double>();
                for (int i = 0; i < Distributions.Count(); i++)
                {
                    result.AddRange(Distributions[i].MaximumOfParameters);
                }
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                SetParameters(MLE(sample));
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <inheritdoc/>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var newDistribution = (CompetingRisks)Clone();
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="distributions">The competing distributions.</param>
        public void SetParameters(UnivariateDistributionBase[] distributions)
        {
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            _distributions = distributions;
            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
            _mvnCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="distributions">The competing distributions.</param>
        public void SetParameters(IUnivariateDistribution[] distributions)
        {
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            _distributions = new UnivariateDistributionBase[distributions.Length];
            for (int i = 0; i < distributions.Length; i++)
            {
                _distributions[i] = (UnivariateDistributionBase)distributions[i];
            }
            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
            _mvnCreated = false;
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentException">Thrown when the flattened parameter count does not match the component distributions.</exception>
        public override void SetParameters(IList<double> parameters)
        {
            if (Distributions == null || Distributions.Count == 0)
            {
                _parametersValid = false;
                return;
            }
            if (parameters.Count != NumberOfParameters)
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }

            int t = 0;
            for (int i = 0; i < Distributions.Count; i++)
            {
                var parms = new List<double>();
                for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                {
                    parms.Add(parameters[j]);
                }
                Distributions[i].SetParameters(parms);
                t += Distributions[i].NumberOfParameters;
            }

            _parametersValid = ValidateParameters(parameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
            _mvnCreated = false;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            ArgumentOutOfRangeException? error = null;
            if (_distributions == null || _distributions.Length == 0 || _distributions.Any(d => d is null))
                error = new ArgumentOutOfRangeException(nameof(Distributions), "There must be at least one non-null distribution.");
            else if (parameters == null || parameters.Count != _distributions.Sum(d => d.GetParameters.Length))
                error = new ArgumentOutOfRangeException(nameof(parameters), "The flattened parameter count must match the component distributions.");
            else
            {
                int offset = 0;
                foreach (var distribution in Distributions)
                {
                    var candidate = new double[distribution.GetParameters.Length];
                    for (int j = 0; j < candidate.Length; j++) candidate[j] = parameters[offset++];
                    error = distribution.ValidateParameters(candidate, false);
                    if (error != null) break;
                }
            }
            if (throwException && error != null) throw error;
            return error;
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample);
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];

            int t = 0;
            for (int i = 0; i < Distributions.Count; i++)
            {
                var tuple = ((IMaximumLikelihoodEstimation)Distributions[i]).GetParameterConstraints(sample);
                var initials = tuple.Item1;
                var lowers = tuple.Item2;
                var uppers = tuple.Item3;

                for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                {
                    initialVals[j] = initials[j - t];
                    lowerVals[j] = lowers[j - t];
                    upperVals[j] = uppers[j - t];
                }
                t += Distributions[i].NumberOfParameters;
            }

            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        public double[] MLE(IList<double> sample)
        {
            // Set constraints
            var tuple = GetParameterConstraints(sample);
            var Initials = tuple.Item1;
            var Lowers = tuple.Item2;
            var Uppers = tuple.Item3;

            // Solve using Nelder-Mead (Downhill Simplex)
            double logLH(double[] x)
            {
                var dist = (CompetingRisks)Clone();
                dist.SetParameters(x);
                return dist.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <inheritdoc/>
        public override double PDF(double x) => Math.Exp(LogPDF(x));

        /// <inheritdoc/>
        /// <remarks>Dependent densities use a support-bounded CDF derivative. An interior
        /// derivative that does not exceed 1E-300 returns negative infinity, retaining the
        /// established rejection of unresolved estimation candidates.</remarks>
        public override double LogPDF(double x)
        {
            ValidateEvaluation();
            if (double.IsNaN(x)) return double.NaN;
            double minimum = Minimum;
            if (x < minimum) return double.NegativeInfinity;
            double maximum = Maximum;
            if (x > maximum || double.IsInfinity(x)) return double.NegativeInfinity;
            if (_distributions.Length == 1) return _distributions[0].LogPDF(x);
            if (Dependency != Probability.DependencyType.Independent)
            {
                double density = DependentDensity(x, minimum, maximum, out bool reuseBounds);
                if (x > (reuseBounds ? minimum : Minimum) && x < (reuseBounds ? maximum : Maximum))
                    return density > 1E-300 ? Math.Log(density) : double.NegativeInfinity;
                return Math.Log(density);
            }

            // Sum f_i times the other factors, without dividing by possibly zero tails.
            // Computing each excluded product also avoids infinity-minus-infinity in the log sum.
            double total = double.NegativeInfinity;
            for (int i = 0; i < Distributions.Count; i++)
            {
                double product = 0;
                for (int j = 0; j < Distributions.Count; j++)
                    if (j != i) product += MinimumOfRandomVariables ? Distributions[j].LogCCDF(x) : Distributions[j].LogCDF(x);
                double logDensity = Distributions[i].LogPDF(x);
                if (double.IsNegativeInfinity(product))
                {
                    if (double.IsPositiveInfinity(logDensity)) return IndependentEndpointLogDensity(x);
                    continue;
                }
                double term = logDensity + product;
                total = DistributionNumerics.LogSum(total, term);
            }
            return total;
        }

        /// <summary>Combines endpoint tail exponents before evaluating the one-sided density limit.</summary>
        /// <remarks>For a product tail c*t^a*log(1/t)^b, its density tends to zero for a&gt;1,
        /// infinity for a&lt;1, and c times the logarithmic limit for a=1.</remarks>
        private double IndependentEndpointLogDensity(double x)
        {
            double power = 0, logPower = 0, logCoefficient = 0;
            bool lower = !MinimumOfRandomVariables;
            foreach (var distribution in Distributions)
            {
                double logTail = lower ? distribution.LogCDF(x) : distribution.LogCCDF(x);
                if (!double.IsNegativeInfinity(logTail)) { logCoefficient += logTail; continue; }
                if (!DistributionEndpointTail.TryExpansion(distribution, lower, out double componentPower,
                    out double componentLogPower, out double componentCoefficient))
                    throw new InvalidOperationException("The independent endpoint density limit is unavailable for this component family.");
                power += componentPower;
                logPower += componentLogPower;
                logCoefficient += componentCoefficient;
            }
            if (power > 1 || (power == 1 && logPower < 0)) return double.NegativeInfinity;
            if (power < 1 || logPower > 0) return double.PositiveInfinity;
            return logCoefficient;
        }

        /// <summary>Checks current component validity, including mutations through public component references.</summary>
        private void ValidateEvaluation()
        {
            if (!_parametersValid || Array.Exists(_distributions, d => !d.ParametersValid)) ValidateParameters(GetParameters, true);
        }

        /// <summary>Numerically differentiates the dependent CDF with a support-bounded stencil.</summary>
        /// <param name="x">Observation at which to evaluate the density.</param>
        /// <param name="minimum">Lower support bound already read by the caller.</param>
        /// <param name="maximum">Upper support bound already read by the caller.</param>
        /// <param name="reuseBounds">Whether exact built-in Weibulls permit reuse of the caller's support bounds.</param>
        /// <returns>The resolved finite CDF derivative.</returns>
        /// <remarks>A centered local step is used in the interior; endpoints use a one-sided
        /// step. Finite negative interior slopes are returned for candidate rejection by
        /// <see cref="LogPDF(double)"/>. Negative boundary or nonfinite density remains a failure.</remarks>
        private double DependentDensity(double x, double minimum, double maximum, out bool reuseBounds)
        {
            var configuration = Volatile.Read(ref _weibullConfiguration);
            if (configuration is not null && !configuration.Matches(this)) configuration = null;
            // Exact built-in Weibulls have fixed support and no custom evaluation callbacks.
            // Derived owners and generic components retain every live support read.
            reuseBounds = configuration is not null && GetType() == typeof(CompetingRisks);
            double step;
            if (configuration is null || !configuration.TryGetDensityStep(out step))
            {
                double scale = double.PositiveInfinity;
                foreach (var distribution in Distributions)
                {
                    double width = distribution.InverseCDF(.75) - distribution.InverseCDF(.25);
                    if (width > 0 && DistributionNumerics.IsFinite(width)) scale = Math.Min(scale, width);
                }
                bool componentScale = DistributionNumerics.IsFinite(scale);
                if (!componentScale) scale = Math.Max(1, Math.Abs(x));
                step = Math.Pow(Tools.DoubleMachineEpsilon, 1.0 / 3) * scale;
                if (componentScale && configuration is not null) configuration.CacheDensityStep(step);
            }
            double left = Math.Max(reuseBounds ? minimum : Minimum, x - step), right = Math.Min(reuseBounds ? maximum : Maximum, x + step);
            if (!(right > left)) throw new InvalidOperationException("The dependent density cannot be resolved at this floating-point scale.");
            // The matching built-in case has already validated this fixed configuration and
            // clamped both endpoints to support. Generic and derived cases retain virtual calls.
            double density = reuseBounds
                ? (DependentCDFCore(right) - DependentCDFCore(left)) / (right - left)
                : (CDF(right) - CDF(left)) / (right - left);
            if (!DistributionNumerics.IsFinite(density) || (density < 0 && (x == (reuseBounds ? minimum : Minimum) || x == (reuseBounds ? maximum : Maximum))))
                throw new InvalidOperationException("Numerical differentiation of the dependent CDF did not produce a nonnegative finite density.");
            return density;
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            ValidateEvaluation();
            if (Dependency != Probability.DependencyType.Independent) return Math.Log(CDF(x));
            if (!MinimumOfRandomVariables) return Distributions.Sum(d => d.LogCDF(x));
            double union = double.NegativeInfinity, precedingSurvival = 0;
            foreach (var distribution in Distributions)
            {
                union = DistributionNumerics.LogSum(union, precedingSurvival + distribution.LogCDF(x));
                precedingSurvival += distribution.LogCCDF(x);
            }
            return union;
        }

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            ValidateEvaluation();
            if (Dependency != Probability.DependencyType.Independent) return DistributionNumerics.Log1mExp(Math.Log(CDF(x)));
            if (MinimumOfRandomVariables) return Distributions.Sum(d => d.LogCCDF(x));
            double union = double.NegativeInfinity, precedingCDF = 0;
            foreach (var distribution in Distributions)
            {
                union = DistributionNumerics.LogSum(union, precedingCDF + distribution.LogCCDF(x));
                precedingCDF += distribution.LogCDF(x);
            }
            return union;
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            ValidateEvaluation();
            if (Dependency == Probability.DependencyType.Independent) return Math.Exp(LogCDF(x));
            RefreshCachedConfiguration();
            if (x < Minimum) return 0;
            if (x > Maximum) return 1;
            if (_distributions.Length == 1)
            {
                return _distributions[0].CDF(x);
            }

            return DependentCDFCore(x);
        }

        /// <summary>Combines component CDFs after validation, configuration refresh, and support checks.</summary>
        /// <param name="x">Observation within the current support.</param>
        /// <returns>The dependent composite probability with the existing probability bounds.</returns>
        /// <remarks>The density reuses this body only for an already validated, unchanged exact
        /// built-in Weibull configuration. Public and generic evaluation retain their live guards.</remarks>
        private double DependentCDFCore(double x)
        {
            double p = double.NaN;
            var ind = new int[_distributions.Length];
            var cdf = new double[_distributions.Length];
            for (int i = 0; i < _distributions.Length; i++)
            {
                ind[i] = 1;
                cdf[i] = _distributions[i].CDF(x);
            }

            if (MinimumOfRandomVariables == true)
            {
                
                if (Dependency == Probability.DependencyType.PerfectlyNegative || Dependency == Probability.DependencyType.CorrelationMatrix)
                {
                    if (_mvnCreated == false)
                        CreateMultivariateNormal();
                    p = Probability.UnionPCM(cdf, _mvn.Covariance);
                }
                else
                {
                    p = Probability.Union(cdf, Dependency);
                }
            }
            else
            {
                if (Dependency == Probability.DependencyType.PerfectlyNegative || Dependency == Probability.DependencyType.CorrelationMatrix)
                {
                    if (_mvnCreated == false)
                        CreateMultivariateNormal();
                    p = Probability.JointProbability(cdf, ind, _mvn.Covariance);
                }
                else
                {
                    p = Probability.JointProbability(cdf, Dependency);
                }
                
            }
            return p < 0d ? 0d : p > 1d ? 1d : p;
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            ValidateEvaluation();
            RefreshCachedConfiguration();
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);

            // If there is only one distribution, return its inverse CDF
            if (Distributions.Count() == 1)
            {
                return Distributions[0].InverseCDF(probability);
            }

            double x = 0;
            if (_empiricalCDFCreated == true)
            {
                x = _empiricalCDF.InverseCDF(probability);
            }
            else
            {
                // use a root finder to solve the inverse CDF
                var xVals = Distributions.Select(d => d.InverseCDF(probability));
                double minX = xVals.Min();
                double maxX = xVals.Max();
                try
                {
                    double reference = minX / 2 + maxX / 2;
                    double scale = maxX / 2 - minX / 2;
                    if (!(scale > 0) || !DistributionNumerics.IsFinite(scale))
                        scale = Distributions.Select(d => d.InverseCDF(.75) - d.InverseCDF(.25))
                            .Where(width => width > 0 && DistributionNumerics.IsFinite(width)).DefaultIfEmpty(1).Min();
                    double Argument(double t)
                    {
                        double value = reference + scale * t;
                        if (double.IsInfinity(value) && DistributionNumerics.IsFinite(t)) value = scale * (reference / scale + t);
                        return Math.Max(Minimum, Math.Min(Maximum, value));
                    }
                    double Residual(double t) => probability <= .5 ? LogCDF(Argument(t)) - Math.Log(probability)
                        : LogCCDF(Argument(t)) - Tools.Log1p(-probability);
                    double lower = -1, upper = 1;
                    Brent.Bracket(Residual, ref lower, ref upper, out _, out _);
                    x = Argument(Brent.Solve(Residual, lower, upper, 1E-6 / Math.Max(1, scale), 100, true));
                }
                catch (Exception)
                {
                    // If the root finder fails, create an empirical CDF
                    if (_empiricalCDFCreated == false)
                        CreateEmpiricalCDF();
                    x = _empiricalCDF.InverseCDF(probability);
                }
            }
            double min = Minimum;
            double max = Maximum;
            return x < min ? min : x > max ? max : x;
        }

        /// <summary>
        /// Returns a list of cumulative incidence functions. 
        /// </summary>
        /// <param name="bins">Optional. The stratification bins to integrate over. Default is 200 bins.</param>
        public List<EmpiricalDistribution> CumulativeIncidenceFunctions(List<StratificationBin>? bins = null)
        {
            RefreshCachedConfiguration();
            // Get stratification bins
            if (bins == null)
            {
                double minP = 1E-16;
                double maxP = 1 - 1E-16;
                double minX = Distributions.Min(d => d.InverseCDF(minP));
                double maxX = Distributions.Max(d => d.InverseCDF(maxP));
                bins = Stratify.XValues(new StratificationOptions(minX, maxX, 200, false), true);
            }

            var D = Distributions.Count();
            var x = new List<double[]>();
            var p = new List<double[]>();
            var dF = new List<double[]>();

            if (Dependency == Probability.DependencyType.PerfectlyNegative || Dependency == Probability.DependencyType.CorrelationMatrix)
            {
                /* 
                * For perfect negative dependency or a custom correlation matrix,
                * use the Genz Method. This method is slow but accurate.
                */

                if (_mvnCreated == false)
                    CreateMultivariateNormal();

                var lower = new double[D];
                var upper = new double[D];
                for (int i = 0; i < D; i++)
                {
                    x.Add(new double[bins.Count + 1]);
                    p.Add(new double[bins.Count + 1]);
                    dF.Add(new double[bins.Count + 1]);

                    // Record the first bin
                    x[i][0] = bins[0].LowerBound;
                    for (int k = 0; k < D; k++)
                    {
                        if (MinimumOfRandomVariables == true)
                        {
                            lower[k] = k == i ? Normal.StandardZ(1E-16) : Normal.StandardZ(Distributions[k].CDF(bins[0].LowerBound));
                            upper[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[0].LowerBound)) : Normal.StandardZ(1 - 1E-16);
                        }
                        else
                        {
                            lower[k] = Normal.StandardZ(1E-16);
                            upper[k] = Normal.StandardZ(Distributions[i].CDF(bins[0].LowerBound));
                        }
                    }
                    dF[i][0] = _mvn.Interval(lower, upper);
                    if (double.IsNaN(dF[i][0])) dF[i][0] = 0;
                    dF[i][0] = Math.Max(0, Math.Min(1, dF[i][0]));
                    p[i][0] = dF[i][0];

                    // Record the remaining bins
                    for (int j = 0; j < bins.Count; j++)
                    {
                        x[i][j + 1] = bins[j].UpperBound;
                        for (int k = 0; k < D; k++)
                        {
                            if (MinimumOfRandomVariables == true)
                            {
                                lower[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[j].LowerBound)) : Normal.StandardZ(Distributions[k].CDF(bins[j].Midpoint));
                                upper[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[j].UpperBound)) : Normal.StandardZ(1 - 1E-16);
                            }
                            else
                            {
                                lower[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[j].LowerBound)) : Normal.StandardZ(1E-16);
                                upper[k] = k == i ? Normal.StandardZ(Distributions[i].CDF(bins[j].UpperBound)) : Normal.StandardZ(Distributions[k].CDF(bins[j].Midpoint));
                            }
                        }
                        dF[i][j + 1] = _mvn.Interval(lower, upper);
                        if (double.IsNaN(dF[i][j + 1])) dF[i][j + 1] = 0;
                        dF[i][j + 1] = Math.Max(0, Math.Min(1, dF[i][j + 1]));
                    }
                }
            }
            else if (Dependency == Probability.DependencyType.Independent || Dependency == Probability.DependencyType.PerfectlyPositive)
            {
                /* 
                * For perfect independent or perfectly positive,
                * use the "Delta Method" developed by Haden Smith and Dave Margo.
                * It is fast and accurate.
                */

                double F1 = 0, F2 = 0;
                var pm = new double[D];
                var pl = new double[D];
                var pu = new double[D];
                var ind = new int[D];
                ind.Fill(1);

                for (int i = 0; i < D; i++)
                {
                    x.Add(new double[bins.Count + 1]);
                    p.Add(new double[bins.Count + 1]);
                    dF.Add(new double[bins.Count + 1]);

                    // Record first bin
                    x[i][0] = bins[0].LowerBound;
                    for (int k = 0; k < D; k++)
                    {
                        if (MinimumOfRandomVariables == true)
                        {
                            pl[k] = k == i ? Distributions[k].CCDF(bins[0].LowerBound) : Distributions[k].CCDF(bins[0].LowerBound);
                            pu[k] = k == i ? 1.0 : Distributions[k].CCDF(bins[0].LowerBound);
                        }
                        else
                        {
                            pu[k] = Distributions[k].CDF(bins[0].LowerBound);
                        }
                    }
                    F1 = Probability.JointProbability(pl, Dependency);
                    F2 = Probability.JointProbability(pu, Dependency);
                    dF[i][0] = F2 - F1;
                    if (double.IsNaN(dF[i][0])) 
                        dF[i][0] = 0;
                    dF[i][0] = Math.Max(0, Math.Min(1, dF[i][0]));
                    p[i][0] = dF[i][0];

                    // Record remaining bins
                    for (int j = 0; j < bins.Count; j++)
                    {
                        x[i][j + 1] = bins[j].UpperBound;
                        for (int k = 0; k < D; k++)
                        {
                            if (MinimumOfRandomVariables == true)
                            {
                                pl[k] = k == i ? Distributions[k].CCDF(bins[j].UpperBound) : Distributions[k].CCDF(bins[j].Midpoint);
                                pu[k] = k == i ? Distributions[k].CCDF(bins[j].LowerBound) : Distributions[k].CCDF(bins[j].Midpoint);
                            }
                            else
                            {
                                pl[k] = k == i ? Distributions[k].CDF(bins[j].LowerBound) : Distributions[k].CDF(bins[j].Midpoint);
                                pu[k] = k == i ? Distributions[k].CDF(bins[j].UpperBound) : Distributions[k].CDF(bins[j].Midpoint);
                            }

                        }
                        F1 = Probability.JointProbability(pl, Dependency);
                        F2 = Probability.JointProbability(pu, Dependency);
                        dF[i][j + 1] = F2 - F1;
                        if (double.IsNaN(dF[i][j + 1])) 
                            dF[i][j + 1] = 0;
                        dF[i][j + 1] = Math.Max(0, Math.Min(1, dF[i][j + 1]));
                    }
                }
            }

            // Get cumulative probabilities and make sure they sum <= 1 across D
            bool fixDF = false;
            var sum = new double[bins.Count + 1];
            for (int j = 1; j <= bins.Count; j++)
            {
                for (int i = 0; i < D; i++)
                {
                    sum[j] += p[i][j - 1] + dF[i][j];
                    p[i][j] = Math.Max(0, Math.Min(1, p[i][j - 1] + dF[i][j]));
                }
                if (sum[j] > 1 && sum[j] != sum[j - 1] && fixDF == false)
                {
                    double s = 0;
                    for (int i = 0; i < D; i++)
                    {
                        dF[i][j] *= (1 - sum[j - 1]) / (sum[j] - sum[j - 1]);
                        s += p[i][j - 1] + dF[i][j];

                        p[i][j] = Math.Max(0, Math.Min(1, p[i][j - 1] + dF[i][j]));
                    }
                    sum[j] = s;
                    fixDF = true;
                }
                else if (fixDF == true)
                {
                    for (int i = 0; i < D; i++)
                    {
                        dF[i][j] = 0;
                        p[i][j] = Math.Max(0, Math.Min(1, p[i][j - 1] + dF[i][j]));
                    }
                }
            }

            // Return CIFs
            var CIFs = new List<EmpiricalDistribution>();
            for (int i = 0; i < D; i++)
                CIFs.Add(new EmpiricalDistribution(x[i], p[i]));

            return CIFs;

        }

        /// <summary>
        /// Validates the user-supplied correlation matrix used by the Gaussian copula.
        /// </summary>
        /// <exception cref="ArgumentException">
        /// Thrown when the matrix is missing, has the wrong dimensions, contains invalid
        /// entries, is not symmetric with unit diagonal, or is not positive definite.
        /// </exception>
        private void ValidateCorrelationMatrix()
        {
            int dimension = Distributions.Count;
            if (CorrelationMatrix == null)
                throw new ArgumentException("A correlation-matrix dependency requires a correlation matrix.", nameof(CorrelationMatrix));
            if (CorrelationMatrix.GetLength(0) != dimension || CorrelationMatrix.GetLength(1) != dimension)
                throw new ArgumentException("The correlation matrix dimensions must match the number of distributions.", nameof(CorrelationMatrix));

            for (int i = 0; i < dimension; i++)
            {
                if (Math.Abs(CorrelationMatrix[i, i] - 1d) > 1E-12)
                    throw new ArgumentException("The correlation matrix must have unit diagonal entries.", nameof(CorrelationMatrix));

                for (int j = 0; j < dimension; j++)
                {
                    double value = CorrelationMatrix[i, j];
                    if (!Tools.IsFinite(value) || value < -1d || value > 1d)
                        throw new ArgumentException("The correlation matrix must contain finite values between -1 and 1.", nameof(CorrelationMatrix));
                    if (j > i && Math.Abs(value - CorrelationMatrix[j, i]) > 1E-12)
                        throw new ArgumentException("The correlation matrix must be symmetric.", nameof(CorrelationMatrix));
                }
            }

            try
            {
                var cholesky = new CholeskyDecomposition(new Matrix(CorrelationMatrix));
                if (!cholesky.IsPositiveDefinite)
                    throw new ArgumentException("The correlation matrix must be positive definite.", nameof(CorrelationMatrix));
            }
            catch (ArgumentException)
            {
                throw;
            }
            catch (Exception exception)
            {
                throw new ArgumentException("The correlation matrix must be positive definite.", nameof(CorrelationMatrix), exception);
            }
        }

        /// <summary>
        /// Create a Multivariate Normal distribution used for modeling dependency between the marginal distributions.
        /// </summary>
        private void CreateMultivariateNormal()
        {
            var D = Distributions.Count();
            var mu = new double[D];
            var sigma = new double[D, D];
            if (Dependency == Probability.DependencyType.PerfectlyNegative)
            {
                double rho = -1d / (D - 1d) + Math.Sqrt(Tools.DoubleMachineEpsilon);
                for (int i = 0; i < D; i++)
                {
                    mu[i] = 0d;
                    for (int j = 0; j < D; j++)
                        sigma[i, j] = i == j ? 1d : rho;
                }
            }
            else
            {
                ValidateCorrelationMatrix();
                for (int i = 0; i < D; i++)
                {
                    mu[i] = 0d;
                    for (int j = 0; j < D; j++)
                        sigma[i, j] = CorrelationMatrix[i, j];
                }
            }
            _mvn = new MultivariateNormal(mu, sigma) { MVNUNI = new MersenneTwister(PRNGSeed) };
            _mvnCreated = true;
        }

        /// <summary>
        /// Create empirical distribution for the CDF.
        /// </summary>
        public void CreateEmpiricalCDF()
        {
            RefreshCachedConfiguration();
            // Get min & max
            double minP = 1E-16;
            double maxP = 1 - 1E-16;
            double minX = Distributions.Min(d => d.InverseCDF(minP));
            double maxX = Distributions.Max(d => d.InverseCDF(maxP));
            // Get number of bins
            double shift = 0;
            if (minX <= 0) shift = Math.Abs(minX) + 1d;
            double min = minX + shift;
            double max = maxX + shift;
            int order = (int)Math.Floor(Math.Log10(max) - Math.Log10(min));
            int binN = Math.Max(200, 100 * order) - 1;
            // Create bins
            var bins = Stratify.XValues(new StratificationOptions(minX, maxX, binN, false), true);
            var xValues = new List<double>();
            var pValues = new List<double>();
            var x = bins.First().LowerBound;
            var p = CDF(bins.First().LowerBound);
            xValues.Add(x);
            pValues.Add(p);
            for (int i = 1; i < bins.Count; i++)
            {
                x = bins[i].LowerBound;
                p = CDF(x);
                if (x > xValues.Last() && p > pValues.Last())
                {
                    xValues.Add(x);
                    pValues.Add(p);
                }
            }
            x = maxX;
            p = CDF(x);
            if (x > xValues.Last() && p > pValues.Last())
            {
                xValues.Add(x);
                pValues.Add(p);
            }
            _empiricalCDF = new EmpiricalDistribution(xValues, pValues) { XTransform = XTransform, ProbabilityTransform = ProbabilityTransform };
            _empiricalCDFCreated = true;
            _momentsComputed = false;
        }

        /// <inheritdoc/>
        public override double[] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            return GenerateRandomValuesWithDependency(sampleSize, seed);
        }

        /// <summary>
        /// Generates random values accounting for dependency structure.
        /// </summary>
        /// <param name="sampleSize"> Size of random sample to generate. </param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <returns>Array of random samples.</returns>
        public double[] GenerateRandomValuesWithDependency(int sampleSize, int seed = -1)
        {
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize];

            if (Dependency == Probability.DependencyType.Independent)
            {
                // Original implementation is correct for independent case
                for (int i = 0; i < sampleSize; i++)
                {
                    double xMin = double.MaxValue;
                    double xMax = double.MinValue;
                    for (int j = 0; j < Distributions.Count; j++)
                    {
                        var x = Distributions[j].InverseCDF(rnd.NextDouble());
                        if (x < xMin) xMin = x;
                        if (x > xMax) xMax = x;
                    }
                    sample[i] = MinimumOfRandomVariables ? xMin : xMax;
                }
            }
            else if (Dependency == Probability.DependencyType.PerfectlyPositive)
            {
                // For perfectly positive dependency, all variables share the same quantile
                for (int i = 0; i < sampleSize; i++)
                {
                    double u = rnd.NextDouble();
                    double xMin = double.MaxValue;
                    double xMax = double.MinValue;
                    for (int j = 0; j < Distributions.Count; j++)
                    {
                        var x = Distributions[j].InverseCDF(u); // Same u for all
                        if (x < xMin) xMin = x;
                        if (x > xMax) xMax = x;
                    }
                    sample[i] = MinimumOfRandomVariables ? xMin : xMax;
                }
            }
            else if (Dependency == Probability.DependencyType.PerfectlyNegative ||
                     Dependency == Probability.DependencyType.CorrelationMatrix)
            {
                // Use Gaussian copula for correlation structure
                if (_mvnCreated == false)
                    CreateMultivariateNormal();

                // Generate correlated standard normal samples
                var mvnSamples = _mvn.GenerateRandomValues(sampleSize, seed);

                for (int i = 0; i < sampleSize; i++)
                {
                    double xMin = double.MaxValue;
                    double xMax = double.MinValue;

                    for (int j = 0; j < Distributions.Count; j++)
                    {
                        // Transform standard normal to uniform via Phi, then to marginal via inverse CDF
                        double z = mvnSamples[i, j];
                        double u = Normal.StandardCDF(z);
                        var x = Distributions[j].InverseCDF(u);

                        if (x < xMin) xMin = x;
                        if (x > xMax) xMax = x;
                    }
                    sample[i] = MinimumOfRandomVariables ? xMin : xMax;
                }
            }

            return sample;
        }


        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            var dists = new UnivariateDistributionBase[Distributions.Count];
            for (int i = 0; i < Distributions.Count; i++)
                dists[i] = Distributions[i].Clone();

            var cr = new CompetingRisks(dists)
            {
                MinimumOfRandomVariables = MinimumOfRandomVariables,
                Dependency = Dependency,
                XTransform = XTransform,
                ProbabilityTransform = ProbabilityTransform,
                PRNGSeed = PRNGSeed
            };
            if (CorrelationMatrix != null)
                cr.CorrelationMatrix = (double[,])CorrelationMatrix.Clone();

            return cr;
        }

        /// <inheritdoc/>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue(nameof(XTransform), XTransform.ToString());
            result.SetAttributeValue(nameof(ProbabilityTransform), ProbabilityTransform.ToString());
            result.SetAttributeValue(nameof(MinimumOfRandomVariables), MinimumOfRandomVariables.ToString());
            result.SetAttributeValue(nameof(Dependency), Dependency.ToString());
            result.SetAttributeValue(nameof(PRNGSeed), PRNGSeed.ToString(CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Distributions), String.Join("|", Distributions.Select(x => x.Type)));
            // Parameters
            var parms = GetParameters;
            var parmStrings = new string[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("Parameters", String.Join("|", parmStrings));

            // Correlation matrix
            var corrMatrixElement = new XElement(nameof(CorrelationMatrix));
            if (CorrelationMatrix != null 
                && CorrelationMatrix.GetLength(0) == Distributions.Count 
                && CorrelationMatrix.GetLength(1) == Distributions.Count)
            {
                int rows = Distributions.Count;
                int cols = Distributions.Count;
                var row = new double[cols];

                for (int i = 0; i < rows; i++)
                {
                    var corrRowElement = new XElement("Correlation_Row");

                    // collect one row of the 2D array
                    for (int j = 0; j < cols; j++)
                    {
                        row[j] = _correlationMatrix[i, j];
                    }

                    // format each double to "G17" with invariant culture
                    var formatted = row.Select(v => v.ToString("G17", CultureInfo.InvariantCulture));
                    // join with '|' and set as the element's text
                    corrRowElement.Value = string.Join("|", formatted);
                    corrMatrixElement.Add(corrRowElement);
                }
            }
            result.Add(corrMatrixElement);

            return result;
        }

        /// <summary>
        /// Creates a competing-risks distribution from its serialized representation.
        /// </summary>
        /// <param name="xElement">The element to deserialize.</param>
        /// <returns>A deserialized competing-risks distribution, or <see langword="null"/> when the element identifies another distribution type.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when serialized configuration, parameters, or correlation data is malformed.</exception>
        /// <remarks>
        /// Deserialization preserves the saved dependency mode even when optional correlation data
        /// is absent. Empty correlation elements and complete component-sized all-zero matrices
        /// written by earlier applications are treated as an unconfigured matrix. This permits
        /// legacy import without declaring the configuration numerically ready: correlation-based
        /// evaluation still calls the strict matrix validator and fails until a valid matrix is set.
        /// </remarks>
        public static CompetingRisks? FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));

            var typeAttribute = xElement.Attribute(nameof(UnivariateDistributionBase.Type));
            if (typeAttribute == null
                || !Enum.TryParse(typeAttribute.Value, out UnivariateDistributionType type)
                || !Enum.IsDefined(typeof(UnivariateDistributionType), type))
                throw new ArgumentException("The serialized distribution type is missing or invalid.", nameof(xElement));
            if (type != UnivariateDistributionType.CompetingRisks) return null;

            var distributionsAttribute = xElement.Attribute(nameof(Distributions));
            if (distributionsAttribute == null || string.IsNullOrWhiteSpace(distributionsAttribute.Value))
                throw new ArgumentException("The serialized competing-risks distribution has no component distributions.", nameof(xElement));

            string[] typeTokens = distributionsAttribute.Value.Split('|');
            var distributions = new UnivariateDistributionBase[typeTokens.Length];
            for (int i = 0; i < typeTokens.Length; i++)
            {
                if (!Enum.TryParse(typeTokens[i], out UnivariateDistributionType componentType)
                    || !Enum.IsDefined(typeof(UnivariateDistributionType), componentType))
                    throw new ArgumentException("The serialized competing-risks distribution contains an invalid component type.", nameof(xElement));
                distributions[i] = UnivariateDistributionFactory.CreateDistribution(componentType);
            }

            var competingRisks = new CompetingRisks(distributions);

            var xTransformAttribute = xElement.Attribute(nameof(XTransform));
            if (xTransformAttribute != null)
            {
                if (!Enum.TryParse(xTransformAttribute.Value, out Transform xTransform)
                    || !Enum.IsDefined(typeof(Transform), xTransform))
                    throw new ArgumentException("The serialized X transform is invalid.", nameof(xElement));
                competingRisks.XTransform = xTransform;
            }

            var probabilityTransformAttribute = xElement.Attribute(nameof(ProbabilityTransform));
            if (probabilityTransformAttribute != null)
            {
                if (!Enum.TryParse(probabilityTransformAttribute.Value, out Transform probabilityTransform)
                    || !Enum.IsDefined(typeof(Transform), probabilityTransform))
                    throw new ArgumentException("The serialized probability transform is invalid.", nameof(xElement));
                competingRisks.ProbabilityTransform = probabilityTransform;
            }

            var minimumAttribute = xElement.Attribute(nameof(MinimumOfRandomVariables));
            if (minimumAttribute != null)
            {
                if (!bool.TryParse(minimumAttribute.Value, out bool minimumOfRandomVariables))
                    throw new ArgumentException("The serialized minimum-selection flag is invalid.", nameof(xElement));
                competingRisks.MinimumOfRandomVariables = minimumOfRandomVariables;
            }

            var dependencyAttribute = xElement.Attribute(nameof(Dependency));
            if (dependencyAttribute != null)
            {
                if (!Enum.TryParse(dependencyAttribute.Value, out Probability.DependencyType dependency)
                    || !Enum.IsDefined(typeof(Probability.DependencyType), dependency))
                    throw new ArgumentException("The serialized dependency type is invalid.", nameof(xElement));
                competingRisks.Dependency = dependency;
            }

            var seedAttribute = xElement.Attribute(nameof(PRNGSeed));
            if (seedAttribute != null)
            {
                if (!int.TryParse(seedAttribute.Value, NumberStyles.Integer, CultureInfo.InvariantCulture, out int seed))
                    throw new ArgumentException("The serialized competing-risks seed is invalid.", nameof(xElement));
                competingRisks.PRNGSeed = seed;
            }

            var parametersAttribute = xElement.Attribute("Parameters");
            if (parametersAttribute == null)
                throw new ArgumentException("The serialized competing-risks parameters are missing.", nameof(xElement));
            string[] parameterTokens = parametersAttribute.Value.Split('|');
            if (parameterTokens.Length != competingRisks.NumberOfParameters)
                throw new ArgumentException("The serialized competing-risks parameter count is invalid.", nameof(xElement));
            var parameters = new double[parameterTokens.Length];
            for (int i = 0; i < parameters.Length; i++)
            {
                if (!double.TryParse(parameterTokens[i], NumberStyles.Any, CultureInfo.InvariantCulture, out parameters[i])
                    || !Tools.IsFinite(parameters[i]))
                    throw new ArgumentException("The serialized competing-risks parameters contain an invalid value.", nameof(xElement));
            }

            int offset = 0;
            for (int i = 0; i < distributions.Length; i++)
            {
                int count = distributions[i].NumberOfParameters;
                var componentParameters = new double[count];
                Array.Copy(parameters, offset, componentParameters, 0, count);
                distributions[i].ValidateParameters(componentParameters, true);
                offset += count;
            }
            competingRisks.SetParameters(parameters);
            if (!competingRisks.ParametersValid)
                throw new ArgumentException("The serialized competing-risks parameters are invalid.", nameof(xElement));

            var correlationElement = xElement.Element(nameof(CorrelationMatrix));
            if (correlationElement != null)
            {
                var correlationRows = correlationElement.Elements("Correlation_Row").ToArray();
                bool containsUnsupportedContent = correlationElement.Elements().Count() != correlationRows.Length
                    || correlationRows.Any(row => row.HasElements)
                    || correlationElement.Nodes().OfType<XText>().Any(text => !string.IsNullOrWhiteSpace(text.Value));
                if (containsUnsupportedContent)
                    throw new ArgumentException("The serialized correlation matrix contains unsupported content.", nameof(xElement));

                if (correlationRows.Length == 0)
                    return competingRisks;

                int dimension = distributions.Length;
                if (correlationRows.Length != dimension)
                    throw new ArgumentException("The serialized correlation matrix has an invalid row count.", nameof(xElement));

                var correlation = new double[dimension, dimension];
                bool allZero = true;
                for (int i = 0; i < dimension; i++)
                {
                    string[] entries = correlationRows[i].Value.Split('|');
                    if (entries.Length != dimension)
                        throw new ArgumentException("The serialized correlation matrix has an invalid column count.", nameof(xElement));
                    for (int j = 0; j < dimension; j++)
                    {
                        if (!double.TryParse(entries[j], NumberStyles.Any, CultureInfo.InvariantCulture, out correlation[i, j])
                            || !Tools.IsFinite(correlation[i, j])
                            || correlation[i, j] < -1d
                            || correlation[i, j] > 1d)
                            throw new ArgumentException("The serialized correlation matrix contains an invalid value.", nameof(xElement));
                        if (correlation[i, j] != 0d)
                            allZero = false;
                    }
                }

                if (allZero)
                    return competingRisks;

                for (int i = 0; i < dimension; i++)
                {
                    if (Math.Abs(correlation[i, i] - 1d) > 1E-12)
                        throw new ArgumentException("The serialized correlation matrix must have unit diagonal entries.", nameof(xElement));
                    for (int j = i + 1; j < dimension; j++)
                    {
                        if (Math.Abs(correlation[i, j] - correlation[j, i]) > 1E-12)
                            throw new ArgumentException("The serialized correlation matrix must be symmetric.", nameof(xElement));
                    }
                }
                competingRisks.CorrelationMatrix = correlation;
            }

            return competingRisks;
        }

    }
}
