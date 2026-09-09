using Numerics.Data;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading;
using System.Xml.Linq;

namespace Numerics.Distributions
{
    /// <summary>
    /// A Mixture distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    /// When zero inflation is enabled, the distribution is a positive-hurdle mixture:
    /// the zero atom has probability <see cref="ZeroWeight"/>, component weights sum to
    /// the remaining mass, and every component is conditioned on a value strictly greater
    /// than zero. <see cref="PDF(double)"/> at zero reports the atom probability under the
    /// mixed Lebesgue-plus-Dirac reference measure.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Mixture_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Mixture : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {

        /// <summary>
        /// Construct new mixture distribution.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public Mixture(double[] weights, UnivariateDistributionBase[] distributions)
        {
            SetParameters(weights, distributions);
        }

        /// <summary>
        /// Construct new mixture distribution.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public Mixture(double[] weights, IUnivariateDistribution[] distributions)
        {
            SetParameters(weights, distributions);
        }

        private double[] _weights = null!;
        private bool _isZeroInflated;
        private double _zeroWeight;
        private UnivariateDistributionBase[] _distributions = null!;
        private EmpiricalDistribution _empiricalCDF = null!;
        private bool _momentsComputed = false;
        private double u1, u2, u3, u4;
        private bool _empiricalCDFCreated = false;
        [NonSerialized] private WeightLogEntry?[]? _logWeightCache;

        /// <summary>An immutable weight/log pair published atomically to concurrent density readers.</summary>
        private sealed class WeightLogEntry
        {
            /// <summary>The exact binary representation of the cached weight.</summary>
            internal readonly long WeightBits;

            /// <summary>The natural logarithm associated with <see cref="WeightBits"/>.</summary>
            internal readonly double LogValue;

            /// <summary>Initializes an immutable cache entry for one exact weight value.</summary>
            /// <param name="weightBits">The binary64 bits of the weight read by the caller.</param>
            /// <param name="logValue">The natural logarithm of that weight.</param>
            internal WeightLogEntry(long weightBits, double logValue)
            {
                WeightBits = weightBits;
                LogValue = logValue;
            }
        }

        /// <summary>
        /// Returns the array of distribution weights.
        /// </summary>
        public double[] Weights => _weights;

        /// <summary>
        /// Returns the array of univariate probability distributions.
        /// </summary>
        public UnivariateDistributionBase[] Distributions => _distributions;

        /// <summary>
        /// Gets or sets whether a separate probability weight is assigned to values less than or equal to zero.
        /// </summary>
        /// <remarks>
        /// Enabling zero inflation rescales finite, nonnegative component weights so their sum
        /// equals <c>1 - <see cref="ZeroWeight"/></c>.
        /// </remarks>
        public bool IsZeroInflated
        {
            get { return _isZeroInflated; }
            set
            {
                _isZeroInflated = value;
                if (_isZeroInflated)
                {
                    NormalizeComponentWeights();
                }
                RefreshConfigurationState();
            }
        }

        /// <summary>
        /// Gets or sets the zero-value probability weight used when the mixture is zero-inflated.
        /// </summary>
        /// <remarks>
        /// When zero inflation is enabled, assigning this property rescales finite, nonnegative
        /// component weights to the remaining probability mass.
        /// </remarks>
        public double ZeroWeight
        {
            get { return _zeroWeight; }
            set
            {
                _zeroWeight = value;
                if (IsZeroInflated)
                {
                    NormalizeComponentWeights();
                }
                RefreshConfigurationState();
            }
        }

        /// <summary>
        /// Rescales valid component weights to the probability mass remaining after zero inflation.
        /// </summary>
        /// <remarks>
        /// Invalid zero weights and invalid component weights are left unchanged so parameter
        /// validation can report the original configuration error.
        /// </remarks>
        private void NormalizeComponentWeights()
        {
            if (_weights is null || _weights.Length == 0 ||
                double.IsNaN(ZeroWeight) || double.IsInfinity(ZeroWeight) ||
                ZeroWeight < 0.0 || ZeroWeight > 1.0)
            {
                return;
            }

            double sum = 0.0;
            for (int i = 0; i < _weights.Length; i++)
            {
                if (double.IsNaN(_weights[i]) || double.IsInfinity(_weights[i]) || _weights[i] < 0.0)
                {
                    return;
                }
                sum += _weights[i];
            }

            if (sum <= 0.0 || double.IsInfinity(sum))
            {
                return;
            }

            double scale = (1.0 - ZeroWeight) / sum;
            for (int i = 0; i < _weights.Length; i++)
            {
                _weights[i] *= scale;
            }
        }

        /// <summary>
        /// Returns the next representable value greater than the supplied value.
        /// </summary>
        /// <param name="value">The starting value.</param>
        /// <returns>The adjacent representable value toward positive infinity.</returns>
        private static double BitIncrement(double value)
        {
            if (double.IsNaN(value) || value == double.PositiveInfinity) return value;
            if (value == 0.0) return double.Epsilon;
            long bits = BitConverter.DoubleToInt64Bits(value);
            return BitConverter.Int64BitsToDouble(value > 0.0 ? bits + 1 : bits - 1);
        }

        /// <summary>
        /// Returns the next representable value less than the supplied value.
        /// </summary>
        /// <param name="value">The starting value.</param>
        /// <returns>The adjacent representable value toward negative infinity.</returns>
        private static double BitDecrement(double value)
        {
            if (double.IsNaN(value) || value == double.NegativeInfinity) return value;
            if (value == 0.0) return -double.Epsilon;
            long bits = BitConverter.DoubleToInt64Bits(value);
            return BitConverter.Int64BitsToDouble(value > 0.0 ? bits - 1 : bits + 1);
        }

        /// <summary>Gets the component log probability above zero without requiring representable probability mass.</summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <returns>The component log survival probability at zero.</returns>
        /// <exception cref="InvalidOperationException">The component does not have a finite, nonpositive log probability above zero.</exception>
        private double PositiveLogMass(int componentIndex)
        {
            double log = Distributions[componentIndex] is Normal normal
                ? normal.LogCCDFAtZero() : Distributions[componentIndex].LogCCDF(0);
            if (!Tools.IsFinite(log) || log > 0) throw new InvalidOperationException("The active component must have positive probability above zero.");
            return log;
        }

        /// <summary>Evaluates the positive-conditional component log density.</summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="x">The value in the component's physical coordinates.</param>
        /// <returns>The component log density conditional on a value above zero, or negative infinity when <paramref name="x"/> is not positive.</returns>
        /// <exception cref="InvalidOperationException">The component has no valid positive log mass.</exception>
        private double PositiveConditionalLogPDF(int componentIndex, double x) => x > 0
            ? Distributions[componentIndex].LogPDF(x) - PositiveLogMass(componentIndex) : double.NegativeInfinity;

        /// <summary>Evaluates the positive-conditional component log CDF through an interval probability.</summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="x">The positive upper evaluation endpoint.</param>
        /// <returns>The log conditional probability in the interval from zero through <paramref name="x"/>, or negative infinity when <paramref name="x"/> is not positive.</returns>
        /// <exception cref="InvalidOperationException">The component has no valid positive log mass.</exception>
        private double PositiveConditionalLogCDF(int componentIndex, double x) => x <= 0 ? double.NegativeInfinity
            : Distributions[componentIndex].LogLikelihood_Intervals(0, x) - PositiveLogMass(componentIndex);

        /// <summary>Evaluates the positive-conditional component log survival directly.</summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="x">The physical-coordinate survival threshold.</param>
        /// <returns>Zero when <paramref name="x"/> is not positive; otherwise, the component log survival probability conditional on a value above zero.</returns>
        /// <exception cref="InvalidOperationException">The component has no valid positive log mass.</exception>
        private double PositiveConditionalLogCCDF(int componentIndex, double x) => x <= 0 ? 0
            : Math.Min(0, Distributions[componentIndex].LogCCDF(x) - PositiveLogMass(componentIndex));

        /// <summary>Inverts a positive-conditional component using a direct log-survival equation.</summary>
        /// <param name="componentIndex">The zero-based component index.</param>
        /// <param name="probability">The conditional cumulative probability in the closed unit interval.</param>
        /// <returns>The component quantile conditional on a value above zero, including its endpoint limits.</returns>
        /// <exception cref="InvalidOperationException">The component has no valid positive mass or a finite quantile bracket cannot be formed.</exception>
        /// <remarks>Retains the existing root tolerance and iteration limit; it avoids constructing
        /// an unconditional probability that can round to one when positive mass is very small.</remarks>
        private double PositiveConditionalQuantile(int componentIndex, double probability)
        {
            var distribution = Distributions[componentIndex];
            double minimum = Math.Max(0, distribution.Minimum);
            if (probability == 0) return minimum;
            if (probability == 1) return distribution.Maximum;
            double target = PositiveLogMass(componentIndex) + Tools.Log1p(-probability);
            double scale = distribution.InverseCDF(.75) - distribution.InverseCDF(.25);
            if (!(scale > 0) || !Tools.IsFinite(scale)) scale = Math.Max(1, Math.Abs(minimum));
            double upper = Math.Min(distribution.Maximum, minimum + scale);
            for (int i = 0; distribution.LogCCDF(upper) > target && i < 1024; i++)
            {
                scale *= 2;
                double next = minimum + scale;
                upper = Math.Min(distribution.Maximum, Tools.IsFinite(next) ? next : double.MaxValue);
            }
            if (!(upper > minimum) || distribution.LogCCDF(upper) > target)
                throw new InvalidOperationException("The positive-conditional quantile could not be bracketed.");
            // Solve in a unit interval so the existing tolerance does not erase a tiny physical scale.
            double width = upper - minimum;
            return minimum + width * Brent.Solve(t => distribution.LogCCDF(minimum + width * t) - target,
                0, 1, 1E-6 / Math.Max(1, width), 100, true);
        }

        private string? _cachedConfiguration;
        [NonSerialized] private DistributionSnapshot? _configurationCache;

        /// <summary>Refreshes cached moments and interpolation when public arrays or nested components change.</summary>
        /// <remarks>A published bitwise snapshot short-circuits the canonical-string serialization on the
        /// unchanged path; a mismatch or an uncapturable component tree falls back to the string
        /// comparison, which remains the deciding authority for cache invalidation.</remarks>
        private void RefreshCachedConfiguration()
        {
            var previous = Volatile.Read(ref _configurationCache);
            if (previous is not null && previous.Matches(this)) return;
            // Capture before canonical serialization: fallback callbacks may mutate their configuration.
            var next = DistributionSnapshot.TryCapture(this);
            string configuration = DistributionNumerics.ConfigurationState(this);
            if (configuration != _cachedConfiguration)
            {
                _cachedConfiguration = configuration;
                _momentsComputed = false;
                _empiricalCDFCreated = false;
            }
            Volatile.Write(ref _configurationCache, next);
        }

        [NonSerialized] private ValidationCertificate? _validationCertificate;

        /// <summary>Publishes the exact bitwise state that has already passed full evaluation validation.</summary>
        /// <remarks>A certificate exists only for capturable component trees, so a bitwise match
        /// proves the full validator - including the zero-inflated positive-mass checks, which are
        /// pure functions of the captured component parameters - already accepted exactly this
        /// state. Any mutation changes the bits and routes the next evaluation back through the
        /// full validator, preserving every exception and its precedence.</remarks>
        private sealed class ValidationCertificate
        {
            /// <summary>The captured configuration whose bits passed the full validator.</summary>
            internal readonly DistributionSnapshot State;

            /// <summary>The lazily published support minimum for the certified state.</summary>
            private CachedValue? _minimum;

            /// <summary>The lazily published support maximum for the certified state.</summary>
            private CachedValue? _maximum;

            /// <summary>An immutable value publication for concurrent readers of unchanged state.</summary>
            private sealed class CachedValue
            {
                /// <summary>The cached value.</summary>
                internal readonly double Value;

                /// <summary>Initializes an immutable value publication.</summary>
                /// <param name="value">The value to publish.</param>
                internal CachedValue(double value) { Value = value; }
            }

            /// <summary>Initializes a certificate for a validated state.</summary>
            /// <param name="state">The captured configuration that passed validation.</param>
            internal ValidationCertificate(DistributionSnapshot state) { State = state; }

            /// <summary>Returns the cached support minimum for the certified state.</summary>
            /// <param name="value">The cached minimum, or zero when none has been published.</param>
            /// <returns><see langword="true"/> when a minimum is available; otherwise, <see langword="false"/>.</returns>
            internal bool TryGetMinimum(out double value)
            {
                var cached = Volatile.Read(ref _minimum);
                value = cached is null ? 0d : cached.Value;
                return cached is not null;
            }

            /// <summary>Publishes the support minimum for the certified state.</summary>
            /// <param name="value">The computed minimum.</param>
            internal void CacheMinimum(double value) => Volatile.Write(ref _minimum, new CachedValue(value));

            /// <summary>Returns the cached support maximum for the certified state.</summary>
            /// <param name="value">The cached maximum, or zero when none has been published.</param>
            /// <returns><see langword="true"/> when a maximum is available; otherwise, <see langword="false"/>.</returns>
            internal bool TryGetMaximum(out double value)
            {
                var cached = Volatile.Read(ref _maximum);
                value = cached is null ? 0d : cached.Value;
                return cached is not null;
            }

            /// <summary>Publishes the support maximum for the certified state.</summary>
            /// <param name="value">The computed maximum.</param>
            internal void CacheMaximum(double value) => Volatile.Write(ref _maximum, new CachedValue(value));
        }

        /// <summary>Checks mutable weights and current component validity before evaluation.</summary>
        /// <exception cref="ArgumentOutOfRangeException">The component collection, weights, zero weight, total mass, or a component's parameters are invalid.</exception>
        /// <remarks>A published certificate for the bitwise-identical state skips re-validation
        /// without allocating; any mutation, and every uncapturable component tree, runs the full
        /// validator exactly as before.</remarks>
        private void ValidateEvaluation()
        {
            var certificate = Volatile.Read(ref _validationCertificate);
            if (certificate is not null && certificate.State.Matches(this)) return;
            ValidateEvaluationSlow();
            // Prefer the refresh-published snapshot instance so steady-state callers can unify
            // the two checks by reference identity instead of walking twice.
            var state = Volatile.Read(ref _configurationCache);
            if (state is null || !state.Matches(this)) state = DistributionSnapshot.TryCapture(this);
            if (state is not null) Volatile.Write(ref _validationCertificate, new ValidationCertificate(state));
        }

        /// <summary>Runs the full evaluation validator against live state.</summary>
        /// <exception cref="ArgumentOutOfRangeException">The component collection, weights, zero weight, total mass, or a component's parameters are invalid.</exception>
        /// <remarks>Validates live component parameters directly so evaluation does not flatten and re-slice
        /// the same state. Sealed Normal components delegate to their existing scalar validator without
        /// allocating parameter arrays; other components retain their list validator.</remarks>
        private void ValidateEvaluationSlow()
        {
            ArgumentOutOfRangeException? error = null;
            if (_distributions is null || _weights is null || _distributions.Length == 0
                || _distributions.Length != _weights.Length || Array.Exists(_distributions, d => d is null))
                error = new ArgumentOutOfRangeException(nameof(Distributions), "At least one non-null component and a matching weight vector are required.");
            else if (IsZeroInflated && (!Tools.IsFinite(ZeroWeight) || ZeroWeight < 0 || ZeroWeight >= 1))
                error = new ArgumentOutOfRangeException(nameof(ZeroWeight), "The zero weight must be finite and in [0,1).");
            else
            {
                int count = Distributions.Length;
                double mass = IsZeroInflated ? ZeroWeight : 0;
                for (int i = 0; i < count; i++)
                {
                    double weight = Weights[i];
                    if (!Tools.IsFinite(weight) || weight < 0 || weight > 1)
                    { error = new ArgumentOutOfRangeException(nameof(Weights), "Weights must be finite and between zero and one."); break; }
                    mass += weight;
                }
                if (error is null && (!Tools.IsFinite(mass) || !mass.AlmostEquals(1, 1E-8)))
                    error = new ArgumentOutOfRangeException(nameof(Weights), "Component and zero weights must sum to one.");
                for (int i = 0; i < count && error is null; i++)
                {
                    if (Distributions[i] is Normal normal)
                        error = normal.ValidateParameters(normal.Mu, normal.Sigma, false);
                    else
                    {
                        double[] parameters = Distributions[i].GetParameters;
                        error = Distributions[i].ValidateParameters(parameters, false);
                    }
                    if (error is null && IsZeroInflated && Weights[i] > 0)
                    {
                        double logMass = Distributions[i] is Normal zeroNormal
                            ? zeroNormal.LogCCDFAtZero() : Distributions[i].LogCCDF(0);
                        if (!Tools.IsFinite(logMass) || logMass > 0)
                            error = new ArgumentOutOfRangeException(nameof(Distributions), "Each active component must have positive probability above zero.");
                    }
                }
            }
            if (error != null) throw error;
        }

        /// <summary>
        /// Refreshes validity and cached results after zero-inflation configuration changes.
        /// </summary>
        private void RefreshConfigurationState()
        {
            if (_weights is null || _distributions is null)
            {
                _parametersValid = false;
            }
            else
            {
                _parametersValid = ValidateParameters(GetParameters, false) is null;
            }
            _momentsComputed = false;
            _empiricalCDFCreated = false;
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
        /// The maximum iterations in the Expectation Maximization algorithm. Default = 1,000. 
        /// </summary>
        public int MaxIterations { get; set; } = 1000;

        /// <summary>
        /// The relative tolerance for convergence. Default = 1E-8.
        /// </summary>
        public double Tolerance { get; set; } = 1E-8;

        /// <summary>
        /// The total number of iterations required to find the MLE.
        /// </summary>
        public int Iterations { get; private set; }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get
            {
                int sum = 0;
                sum += Distributions.Count();
                for (int i = 0; i < Distributions.Length; i++)
                    sum += Distributions[i].NumberOfParameters;
                return sum;
            }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type => UnivariateDistributionType.Mixture;

        /// <inheritdoc/>
        public override string DisplayName => "Mixture";

        /// <inheritdoc/>
        public override string ShortDisplayName => "MIX";

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                string Wstring = "{";
                string Dstring = "{";
                for (int i = 0; i < Weights.Count(); i++)
                {
                    Wstring += Weights[i].ToString();
                    Dstring += Distributions[i].DisplayName;
                    if (i < Weights.Count() - 1)
                    {
                        Wstring += ",";
                        Dstring += ",";
                    }
                }
                Wstring += "}";
                Dstring += "}";
                parmString[0, 0] = "Weights";
                parmString[1, 0] = "Distributions";
                parmString[0, 1] = Wstring;
                parmString[1, 1] = Dstring;
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNames
        {
            get
            {
                var result = new List<string>();
                for (int i = 1; i <= Distributions.Count(); i++)
                {
                    result.Add("Weight " + i.ToString());
                }
                for (int i = 0; i < Distributions.Length; i++)
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
                for (int i = 1; i <= Distributions.Count(); i++)
                {
                    result.Add("W" + i.ToString());
                }
                for (int i = 0; i < Distributions.Length; i++)
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
                result.AddRange(Weights);
                for (int i = 0; i < Distributions.Length; i++)
                {
                    result.AddRange(Distributions[i].GetParameters);
                }                  
                return result.ToArray();
            }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Weights), nameof(Distributions)]; }
        }

        /// <summary>
        /// Compute central moments of the distribution.
        /// </summary>
        private void ComputeMoments()
        {
            ValidateEvaluation();
            var components = new List<(double weight, double mean, double sd, double skew, double kurt)>();
            if (IsZeroInflated && ZeroWeight > 0) components.Add((ZeroWeight, 0, 0, 0, 0));
            for (int i = 0; i < Distributions.Length; i++)
            {
                if (Weights[i] == 0) continue;
                var distribution = Distributions[i];
                if (IsZeroInflated && distribution.LogCDF(0) != double.NegativeInfinity)
                {
                    int index = i;
                    double center = PositiveConditionalQuantile(i, .5);
                    double scale = PositiveConditionalQuantile(i, .75) - PositiveConditionalQuantile(i, .25);
                    var moments = DistributionMomentIntegration.Compute(x => PositiveConditionalLogPDF(index, x),
                        Math.Max(0, distribution.Minimum), distribution.Maximum, center, scale);
                    components.Add((Weights[i], moments[0], moments[1], moments[2], moments[3]));
                }
                else components.Add((Weights[i], distribution.Mean, distribution.StandardDeviation, distribution.Skewness, distribution.Kurtosis));
            }
            double reference = components[0].mean, offset = 0, totalWeight = 0;
            foreach (var component in components) { offset += component.weight * (component.mean - reference); totalWeight += component.weight; }
            u1 = reference * totalWeight + offset;
            if (!Tools.IsFinite(u1))
            {
                double magnitude = components.Max(c => Math.Abs(c.mean));
                u1 = magnitude * components.Sum(c => c.weight * (c.mean / magnitude));
            }
            double scaleMoment = components.Max(c => Math.Max(c.sd, double.IsInfinity(c.mean - u1)
                && Tools.IsFinite(c.mean) && Tools.IsFinite(u1) ? Math.Max(Math.Abs(c.mean), Math.Abs(u1)) : Math.Abs(c.mean - u1)));
            if (!Tools.IsFinite(scaleMoment)) { u2 = scaleMoment; u3 = u4 = double.NaN; _momentsComputed = true; return; }
            if (scaleMoment == 0) { u2 = 0; u3 = u4 = double.NaN; _momentsComputed = true; return; }
            double m2 = 0, m3 = 0, m4 = 0;
            foreach (var component in components)
            {
                double d = DistributionNumerics.Standardize(component.mean, u1, scaleMoment), sd = component.sd / scaleMoment;
                double v = sd * sd, d2 = d * d;
                double third = sd == 0 ? 0 : component.skew * v * sd;
                double fourth = sd == 0 ? 0 : component.kurt * v * v;
                m2 += component.weight * (v + d2);
                m3 += component.weight * (third + 3 * d * v + d * d2);
                m4 += component.weight * (fourth + 4 * d * third + 6 * d2 * v + d2 * d2);
            }
            u2 = scaleMoment * Math.Sqrt(m2);
            u3 = m3 / m2 / Math.Sqrt(m2);
            u4 = m4 / m2 / m2;
            _momentsComputed = true;
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                RefreshCachedConfiguration();
                if (!_momentsComputed)
                    ComputeMoments();
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
                if (!_momentsComputed)
                    ComputeMoments();
                return u2;
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                RefreshCachedConfiguration();
                if (!_momentsComputed)
                    ComputeMoments();
                return u3;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                RefreshCachedConfiguration();
                if (!_momentsComputed)
                    ComputeMoments();
                return u4;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get
            {
                var certificate = Volatile.Read(ref _validationCertificate);
                if (certificate is not null && certificate.TryGetMinimum(out double cached)
                    && certificate.State.Matches(this)) return cached;
                ValidateEvaluation();
                double minimum;
                if (IsZeroInflated && ZeroWeight > 0) minimum = 0;
                else
                {
                    minimum = Distributions.Where((d, i) => Weights[i] > 0).Min(d => d.Minimum);
                    if (IsZeroInflated) minimum = Math.Max(0, minimum);
                }
                certificate = Volatile.Read(ref _validationCertificate);
                if (certificate is not null && certificate.State.Matches(this)) certificate.CacheMinimum(minimum);
                return minimum;
            }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get
            {
                var certificate = Volatile.Read(ref _validationCertificate);
                if (certificate is not null && certificate.TryGetMaximum(out double cached)
                    && certificate.State.Matches(this)) return cached;
                ValidateEvaluation();
                double maximum = Distributions.Where((d, i) => Weights[i] > 0).Max(d => d.Maximum);
                certificate = Volatile.Read(ref _validationCertificate);
                if (certificate is not null && certificate.State.Matches(this)) certificate.CacheMaximum(maximum);
                return maximum;
            }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        { 
            get 
            {
                var result = new List<double>();
                if (IsZeroInflated) { result.Add(0.0); }
                for (int i = 0; i < Distributions.Length; i++)
                {
                    result.Add(0.0);
                }
                for (int i = 0; i < Distributions.Length; i++)
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
                if (IsZeroInflated) { result.Add(1.0); }
                for (int i = 0; i < Distributions.Length; i++)
                {
                    result.Add(1.0);
                }
                for (int i = 0; i < Distributions.Length; i++)
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
            var newDistribution = (Mixture)Clone();
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public void SetParameters(double[] weights, UnivariateDistributionBase[] distributions)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            if (weights.Length != distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));

            _weights = weights.ToArray();
            _distributions = distributions.ToArray();
            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="distributions">The mixture distributions.</param>
        public void SetParameters(double[] weights, IUnivariateDistribution[] distributions)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (distributions == null) throw new ArgumentNullException(nameof(Distributions));
            if (weights.Length != distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));

            _weights = weights.ToArray();
            _distributions = new UnivariateDistributionBase[distributions.Length];
            for (int i = 0; i < distributions.Length; i++)
            {
                _distributions[i] = (UnivariateDistributionBase)distributions[i];
            }
            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="weights">The mixture weights.</param>
        /// <param name="parameters">The mixture distribution parameters.</param>
        public void SetParameters(double[] weights, double[] parameters)
        {
            if (weights == null) throw new ArgumentNullException(nameof(Weights));
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (weights.Length != Distributions.Length)
                throw new ArgumentException("The weight and distribution arrays must have the same length.", nameof(Weights));
            if (parameters.Length != Distributions.Sum(x => x.NumberOfParameters))
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }

            double[] parameterCopy = parameters.ToArray();

            // Set weights
            _weights = weights.ToArray();
            // Set distribution parameters
            int t = 0;
            for (int i = 0; i < Distributions.Length; i++)
            {
                var parms = new List<double>();
                for (int j = t; j < t + Distributions[i].NumberOfParameters; j++)
                {
                    parms.Add(parameterCopy[j]);
                }
                Distributions[i].SetParameters(parms);
                t += Distributions[i].NumberOfParameters;
            }
            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count != NumberOfParameters)
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }

            double[] parameterCopy = parameters.ToArray();

            // Set the weights.
            int parameterIndex = 0;
            for (int i = 0; i < Distributions.Length; i++)
            {
                Weights[i] = parameterCopy[parameterIndex++];
            }

            // Set the distribution parameters.
            for (int i = 0; i < Distributions.Length; i++)
            {
                double[] distributionParameters = parameterCopy
                    .Skip(parameterIndex)
                    .Take(Distributions[i].NumberOfParameters)
                    .ToArray();
                Distributions[i].SetParameters(distributionParameters);
                parameterIndex += Distributions[i].NumberOfParameters;
            }

            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <summary>
        /// Set the distribution parameters from a referenced array. Weights are normalized to the configured simplex.
        /// </summary>
        /// <param name="parameters">The array of parameters. The caller's array is not modified.</param>
        public void SetParameters(ref double[] parameters)
        {
            if (parameters == null) return;
            if (Weights == null || Weights.Length == 0) return;
            if (Distributions == null || Distributions.Count() == 0) return;

            double[] parameterCopy = parameters.ToArray();
            if (Distributions.Count() == 1 && parameterCopy.Length == Distributions[0].NumberOfParameters)
            {
                Weights[0] = IsZeroInflated ? 1.0 - ZeroWeight : 1.0;
                Distributions[0].SetParameters(parameterCopy);
            }
            else
            {
                int componentCount = Distributions.Count();
                int parameterIndex = componentCount;
                double weightSum = 0.0;

                for (int i = 0; i < componentCount; i++)
                {
                    Weights[i] = parameterCopy[i];
                    weightSum += Weights[i];
                }

                double componentMass = IsZeroInflated ? 1.0 - ZeroWeight : 1.0;
                if (weightSum <= 0.0 || !Tools.IsFinite(weightSum))
                {
                    double uniformWeight = componentMass / componentCount;
                    for (int i = 0; i < componentCount; i++) Weights[i] = uniformWeight;
                }
                else
                {
                    double scale = componentMass / weightSum;
                    for (int i = 0; i < componentCount; i++) Weights[i] *= scale;
                }

                for (int i = 0; i < componentCount; i++)
                {
                    double[] distributionParameters = parameterCopy
                        .Skip(parameterIndex)
                        .Take(Distributions[i].NumberOfParameters)
                        .ToArray();
                    Distributions[i].SetParameters(distributionParameters);
                    parameterIndex += Distributions[i].NumberOfParameters;
                }
            }

            _parametersValid = ValidateParameters(GetParameters, false) is null;
            _momentsComputed = false;
            _empiricalCDFCreated = false;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            ArgumentOutOfRangeException? error = null;
            if (_distributions is null || _weights is null || _distributions.Length == 0
                || _distributions.Length != _weights.Length || _distributions.Any(d => d is null))
                error = new ArgumentOutOfRangeException(nameof(Distributions), "At least one non-null component and a matching weight vector are required.");
            else if (parameters is null || parameters.Count != _weights.Length + _distributions.Sum(d => d.GetParameters.Length))
                error = new ArgumentOutOfRangeException(nameof(parameters), "The flattened parameter count must match the mixture.");
            else if (IsZeroInflated && (!Tools.IsFinite(ZeroWeight) || ZeroWeight < 0 || ZeroWeight >= 1))
                error = new ArgumentOutOfRangeException(nameof(ZeroWeight), "The zero weight must be finite and in [0,1).");
            else
            {
                int count = Distributions.Length;
                double mass = IsZeroInflated ? ZeroWeight : 0;
                for (int i = 0; i < count; i++)
                {
                    double weight = parameters[i];
                    if (!Tools.IsFinite(weight) || weight < 0 || weight > 1)
                    { error = new ArgumentOutOfRangeException(nameof(Weights), "Weights must be finite and between zero and one."); break; }
                    mass += weight;
                }
                if (error is null && (!Tools.IsFinite(mass) || !mass.AlmostEquals(1, 1E-8)))
                    error = new ArgumentOutOfRangeException(nameof(Weights), "Component and zero weights must sum to one.");
                int offset = count;
                for (int i = 0; i < count && error is null; i++)
                {
                    // Nonparametric components expose no flattened scalar parameters.
                    var candidate = new double[Distributions[i].GetParameters.Length];
                    for (int j = 0; j < candidate.Length; j++) candidate[j] = parameters[offset++];
                    error = Distributions[i].ValidateParameters(candidate, false);
                    if (error is null && IsZeroInflated && parameters[i] > 0)
                    {
                        var distribution = Distributions[i];
                        if (!candidate.SequenceEqual(distribution.GetParameters))
                        { distribution = distribution.Clone(); distribution.SetParameters(candidate); }
                        double logMass = distribution.LogCCDF(0);
                        if (!Tools.IsFinite(logMass) || logMass > 0)
                            error = new ArgumentOutOfRangeException(nameof(Distributions), "Each active component must have positive probability above zero.");
                    }
                }
            }
            if (throwException && error != null) throw error;
            return error;
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];

            // Weights are first
            int t = 0;
            for (int i = 0; i < Distributions.Length; i++)
            {
                initialVals[i] = IsZeroInflated ? (1d - ZeroWeight) / Distributions.Count() : 1d / Distributions.Count();
                lowerVals[i] = 0.0;
                upperVals[i] = 1.0;
                t += 1;
            }

            for (int i = 0; i < Distributions.Length; i++)
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
            ValidateParameters(GetParameters, true);

            int observationCount = sample.Count;
            int distributionParameterCount = Distributions.Sum(x => x.NumberOfParameters);
            int componentCount = Distributions.Count();

            for (int rowIndex = 0; rowIndex < observationCount; rowIndex++)
                if (sample[rowIndex] < 0 && Distributions.All(d => d is GammaDistribution || d is Weibull
                    || d is LnNormal || d is LogNormal || d is LogPearsonTypeIII))
                    throw CreateImpossibleRowException(rowIndex, sample[rowIndex]);

            if (IsZeroInflated)
            {
                for (int rowIndex = 0; rowIndex < observationCount; rowIndex++)
                {
                    if (sample[rowIndex] < 0.0)
                    {
                        throw new InvalidOperationException(
                            "Mixture EM row " + rowIndex.ToString(CultureInfo.InvariantCulture) +
                            " has negative exact value " + sample[rowIndex].ToString("R", CultureInfo.InvariantCulture) +
                            " in a zero-inflated model.");
                    }
                    if (sample[rowIndex] == 0.0 && ZeroWeight == 0.0)
                        throw CreateImpossibleRowException(rowIndex, sample[rowIndex]);
                }
            }

            Tuple<double[], double[], double[]> constraints = GetParameterConstraints(sample);
            double[] initialParameters = constraints.Item1.Subset(componentCount);
            double[] lowerParameters = constraints.Item2.Subset(componentCount);
            double[] upperParameters = constraints.Item3.Subset(componentCount);

            double[] mleWeights = constraints.Item1.Subset(0, componentCount - 1);
            double[] mleParameters = initialParameters;
            var responsibilities = new double[observationCount, componentCount];
            double oldLogLikelihood = double.MinValue;
            double newLogLikelihood = double.MinValue;

            double EStep(double[] parameters)
            {
                var distribution = (Mixture)Clone();
                distribution.SetParameters(mleWeights, parameters);
                double logLikelihood = 0.0;

                for (int rowIndex = 0; rowIndex < observationCount; rowIndex++)
                {
                    double value = sample[rowIndex];
                    if (IsZeroInflated && value == 0.0)
                    {
                        if (!Tools.IsFinite(ZeroWeight) || ZeroWeight <= 0.0)
                        {
                            throw CreateImpossibleRowException(rowIndex, value);
                        }

                        for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                        {
                            responsibilities[rowIndex, componentIndex] = 0.0;
                        }
                        logLikelihood += Math.Log(ZeroWeight);
                        continue;
                    }

                    var componentLogs = new double[componentCount];
                    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                    {
                        if (mleWeights[componentIndex] == 0)
                        {
                            componentLogs[componentIndex] = double.NegativeInfinity;
                            continue;
                        }
                        componentLogs[componentIndex] = IsZeroInflated
                            ? distribution.PositiveConditionalLogPDF(componentIndex, value)
                            : distribution.Distributions[componentIndex].LogPDF(value);
                    }
                    var rowResponsibilities = new double[componentCount];
                    double rowLogProbability = MixtureLogWeights.Normalize(componentLogs, mleWeights, rowResponsibilities);
                    if (!Tools.IsFinite(rowLogProbability))
                    {
                        throw CreateImpossibleRowException(rowIndex, value);
                    }

                    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                    {
                        responsibilities[rowIndex, componentIndex] = rowResponsibilities[componentIndex];
                    }
                    logLikelihood += rowLogProbability;
                }

                return logLikelihood;
            }

            double[] MStep(double[] parameters)
            {
                for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                {
                    double weight = 0.0;
                    for (int rowIndex = 0; rowIndex < observationCount; rowIndex++)
                    {
                        if (!IsZeroInflated || sample[rowIndex] > 0.0)
                        {
                            weight += responsibilities[rowIndex, componentIndex];
                        }
                    }
                    mleWeights[componentIndex] = weight;
                }

                double componentWeightSum = mleWeights.Sum();
                double componentWeightTarget = IsZeroInflated ? 1.0 - ZeroWeight : 1.0;
                if (!Tools.IsFinite(componentWeightSum) || componentWeightSum <= 0.0)
                {
                    throw new InvalidOperationException("Mixture EM cannot update component weights because no finite positive responsibility mass is available.");
                }

                double scale = componentWeightTarget / componentWeightSum;
                for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                {
                    mleWeights[componentIndex] *= scale;
                }

                var solver = new NelderMead(Objective, distributionParameterCount, parameters, lowerParameters, upperParameters);
                solver.Maximize();
                return solver.BestParameterSet.Values;
            }

            double Objective(double[] parameters)
            {
                var distribution = (Mixture)Clone();
                distribution.SetParameters(mleWeights, parameters);
                double logLikelihood = distribution.LogLikelihood(sample);
                return Tools.IsFinite(logLikelihood) ? logLikelihood : double.NegativeInfinity;
            }

            InvalidOperationException CreateImpossibleRowException(int rowIndex, double value)
            {
                return new InvalidOperationException(
                    "Mixture EM row " + rowIndex.ToString(CultureInfo.InvariantCulture) +
                    " with value " + value.ToString("R", CultureInfo.InvariantCulture) +
                    " has zero or nonfinite total probability.");
            }

            for (Iterations = 1; Iterations <= MaxIterations; Iterations++)
            {
                newLogLikelihood = EStep(mleParameters);
                if (Math.Abs((oldLogLikelihood - newLogLikelihood) / oldLogLikelihood) < Tolerance) break;
                mleParameters = MStep(mleParameters);
                oldLogLikelihood = newLogLikelihood;
            }

            var result = new List<double>();
            result.AddRange(mleWeights);
            result.AddRange(mleParameters);
            return result.ToArray();
        }

        /// <inheritdoc/>
        public override double PDF(double x) => Math.Exp(LogPDF(x));

        /// <summary>Reuses the exact logarithm of a weight read at its existing post-callback evaluation point.</summary>
        /// <param name="index">The component index.</param>
        /// <param name="weight">The live weight already read after the component callback.</param>
        /// <returns>The value of <see cref="Math.Log(double)"/> for the supplied weight.</returns>
        /// <remarks>Entries are immutable and published with release/acquire semantics. No other live weight
        /// is read, and field-based deserialization starts with an empty transient cache.</remarks>
        private double LogWeight(int index, double weight)
        {
            var cache = Volatile.Read(ref _logWeightCache);
            if (cache is null || index >= cache.Length)
            {
                var expanded = new WeightLogEntry?[_weights.Length];
                if (cache is not null) Array.Copy(cache, expanded, cache.Length);
                Volatile.Write(ref _logWeightCache, expanded);
                cache = expanded;
            }
            long bits = BitConverter.DoubleToInt64Bits(weight);
            var entry = Volatile.Read(ref cache[index]);
            if (entry is not null && entry.WeightBits == bits) return entry.LogValue;
            double log = Math.Log(weight);
            Volatile.Write(ref cache[index], new WeightLogEntry(bits, log));
            return log;
        }

        /// <inheritdoc/>
        public override double LogPDF(double x)
        {
            ValidateEvaluation();
            if (IsZeroInflated && x <= 0) return x == 0 ? Math.Log(ZeroWeight) : double.NegativeInfinity;
            double total = double.NegativeInfinity;
            for (int i = 0; i < Distributions.Length; i++)
            {
                if (Weights[i] == 0) continue;
                double log = IsZeroInflated ? PositiveConditionalLogPDF(i, x) : Distributions[i].LogPDF(x);
                total = DistributionNumerics.LogSum(total, LogWeight(i, Weights[i]) + log);
            }
            return total;
        }

        /// <inheritdoc/>
        public override double CDF(double x) => Math.Exp(LogCDF(x));

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            ValidateEvaluation();
            if (IsZeroInflated && x < 0) return double.NegativeInfinity;
            double total = IsZeroInflated ? Math.Log(ZeroWeight) : double.NegativeInfinity;
            for (int i = 0; i < Distributions.Length; i++)
            {
                if (Weights[i] == 0) continue;
                double log = IsZeroInflated ? PositiveConditionalLogCDF(i, x) : Distributions[i].LogCDF(x);
                total = DistributionNumerics.LogSum(total, LogWeight(i, Weights[i]) + log);
            }
            return Math.Min(0, total);
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            ValidateEvaluation();
            if (IsZeroInflated && x < 0) return 0;
            double total = double.NegativeInfinity;
            for (int i = 0; i < Distributions.Length; i++)
            {
                if (Weights[i] == 0) continue;
                double log = IsZeroInflated ? PositiveConditionalLogCCDF(i, x) : Distributions[i].LogCCDF(x);
                total = DistributionNumerics.LogSum(total, LogWeight(i, Weights[i]) + log);
            }
            return Math.Min(0, total);
        }

        /// <summary>Combines component interval log probabilities while retaining the hurdle atom's endpoint convention.</summary>
        /// <param name="lower">The open lower interval endpoint.</param>
        /// <param name="upper">The closed upper interval endpoint.</param>
        /// <returns>The logarithm of the total mixture probability in <c>(lower, upper]</c>, including the zero hurdle mass when applicable.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The live mixture configuration or a component interval is invalid.</exception>
        internal double LogIntervalProbability(double lower, double upper)
        {
            ValidateEvaluation();
            double total = IsZeroInflated && lower < 0 && upper >= 0 ? Math.Log(ZeroWeight) : double.NegativeInfinity;
            if (IsZeroInflated && upper <= 0) return total;
            for (int i = 0; i < Distributions.Length; i++)
            {
                if (Weights[i] == 0) continue;
                double log = Distributions[i].LogLikelihood_Intervals(IsZeroInflated ? Math.Max(0, lower) : lower, upper);
                if (IsZeroInflated) log -= Distributions[i].LogCCDF(0);
                total = DistributionNumerics.LogSum(total, LogWeight(i, Weights[i]) + log);
            }
            return total;
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // One snapshot walk covers refresh and validation on the unchanged path: the
            // certificate holding the refresh-published instance proves both checks at once.
            var certificate = Volatile.Read(ref _validationCertificate);
            bool certified = certificate is not null
                && ReferenceEquals(Volatile.Read(ref _configurationCache), certificate.State)
                && certificate.State.Matches(this);
            if (!certified)
            {
                RefreshCachedConfiguration();
                certificate = null;
            }
            if (!(probability >= 0.0 && probability <= 1.0))
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between 0 and 1.");
            if (!certified) ValidateEvaluation();
            if (probability == 0.0) return Minimum;
            if (probability == 1.0) return Maximum;
            if (IsZeroInflated && probability <= ZeroWeight) return 0.0;
            if (Distributions.Length == 1 && !IsZeroInflated)
            {
                return Distributions[0].InverseCDF(probability);
            }

            if (_empiricalCDFCreated)
            {
                double empiricalValue = _empiricalCDF.InverseCDF(probability);
                double clampMinimum = certificate is not null && certificate.TryGetMinimum(out double cachedMinimum)
                    ? cachedMinimum : Minimum;
                double clampMaximum = certificate is not null && certificate.TryGetMaximum(out double cachedMaximum)
                    ? cachedMaximum : Maximum;
                return Tools.Clamp(empiricalValue, clampMinimum, clampMaximum);
            }

            double componentProbability = IsZeroInflated
                ? (probability - ZeroWeight) / (1.0 - ZeroWeight)
                : probability;
            var componentQuantiles = new List<double>();
            for (int i = 0; i < Distributions.Length; i++)
            {
                if (Weights[i] == 0) continue;
                componentQuantiles.Add(IsZeroInflated ? PositiveConditionalQuantile(i, componentProbability)
                    : Distributions[i].InverseCDF(componentProbability));
            }

            double lowerBound = componentQuantiles.Min();
            double upperBound = componentQuantiles.Max();
            double value;
            try
            {
                if (lowerBound == upperBound) return Tools.Clamp(lowerBound, Minimum, Maximum);
                double width = upperBound - lowerBound;
                double Argument(double t) => Tools.IsFinite(width) ? lowerBound + width * t : (1 - t) * lowerBound + t * upperBound;
                double Residual(double t) => probability <= .5 ? LogCDF(Argument(t)) - Math.Log(probability)
                    : LogCCDF(Argument(t)) - Tools.Log1p(-probability);
                double scale = Tools.IsFinite(width) ? width : Math.Max(Math.Abs(lowerBound), Math.Abs(upperBound));
                value = Argument(Brent.Solve(Residual, 0, 1, 1E-6 / Math.Max(1, scale), 100, true));
            }
            catch (Exception)
            {
                if (!_empiricalCDFCreated) CreateEmpiricalCDF();
                value = _empiricalCDF.InverseCDF(probability);
            }

            return Tools.Clamp(value, Minimum, Maximum);
        }

        /// <summary>
        /// Create empirical distribution for the CDF.
        /// </summary>
        public void CreateEmpiricalCDF()
        {
            // Get min & max
            double minP = 1E-16;
            double maxP = 1 - 1E-16;
            RefreshCachedConfiguration();
            var activeIndices = Enumerable.Range(0, Distributions.Length).Where(i => Weights[i] > 0).ToArray();
            double minX = IsZeroInflated ? 0 : activeIndices.Min(i => Distributions[i].InverseCDF(minP));
            double maxX = activeIndices.Max(i => IsZeroInflated ? PositiveConditionalQuantile(i, maxP) : Distributions[i].InverseCDF(maxP));
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
            ValidateEvaluation();

            var random = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize];
            int lastActiveComponent = Array.FindLastIndex(Weights, weight => weight > 0);
            for (int sampleIndex = 0; sampleIndex < sampleSize; sampleIndex++)
            {
                double mixtureProbability = random.NextDouble();
                double componentProbability = random.NextDouble();
                if (IsZeroInflated && mixtureProbability <= ZeroWeight)
                {
                    sample[sampleIndex] = 0.0;
                    continue;
                }

                double cumulativeWeight = IsZeroInflated ? ZeroWeight : 0.0;
                for (int componentIndex = 0; componentIndex < Distributions.Count(); componentIndex++)
                {
                    if (Weights[componentIndex] == 0) continue;
                    cumulativeWeight += Weights[componentIndex];
                    if (mixtureProbability <= cumulativeWeight || componentIndex == lastActiveComponent)
                    {
                        double probability = componentProbability;
                        if (IsZeroInflated)
                        {
                            sample[sampleIndex] = PositiveConditionalQuantile(componentIndex, componentProbability);
                            break;
                        }
                        sample[sampleIndex] = Distributions[componentIndex].InverseCDF(probability);
                        break;
                    }
                }
            }

            return sample;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            var dists = new UnivariateDistributionBase[Distributions.Count()];
            for (int i = 0; i < Distributions.Length; i++)
                dists[i] = Distributions[i].Clone();

            return new Mixture(Weights.ToArray(), dists)
            {
                IsZeroInflated = IsZeroInflated,
                ZeroWeight = ZeroWeight,
                XTransform = XTransform,
                ProbabilityTransform = ProbabilityTransform
            };
        }

        /// <inheritdoc/>
        public override XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            result.SetAttributeValue(nameof(IsZeroInflated), IsZeroInflated.ToString());
            result.SetAttributeValue(nameof(ZeroWeight), ZeroWeight.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(XTransform), XTransform.ToString());
            result.SetAttributeValue(nameof(ProbabilityTransform), ProbabilityTransform.ToString());
            result.SetAttributeValue(nameof(Distributions), String.Join("|", Distributions.Select(x => x.Type)));
            // Weights
            var weights = Weights;
            var weightStrings = new string[Weights.Length];
            for (int i = 0; i < Weights.Length; i++)
            {
                weightStrings[i] = weights[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue(nameof(Weights), String.Join("|", weightStrings));
            // Parameters
            var parms = GetParameters;
            var parmStrings = new string[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                parmStrings[i] = parms[i].ToString("G17", CultureInfo.InvariantCulture);
            }
            result.SetAttributeValue("Parameters", String.Join("|", parmStrings));
            return result;
        }

        /// <summary>
        /// Create a mixture distribution from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new mixture distribution.</returns>
        public static Mixture? FromXElement(XElement xElement)
        {
            UnivariateDistributionType type = UnivariateDistributionType.Deterministic;
            var typeAttr = xElement.Attribute(nameof(UnivariateDistributionBase.Type));
            if (typeAttr != null)
            {
                Enum.TryParse(typeAttr.Value, out type);

            }
            if (type == UnivariateDistributionType.Mixture)
            {
                var weights = new List<double>();
                var distributions = new List<UnivariateDistributionBase>();
                var weightsAttr = xElement.Attribute(nameof(Weights));
                if (weightsAttr != null)
                {
                    var w = weightsAttr.Value.Split('|');
                    for (int i = 0; i < w.Length; i++)
                    {
                        double.TryParse(w[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var weight);
                        weights.Add(weight);
                    }
                }
                var distsAttr = xElement.Attribute(nameof(Distributions));
                if (distsAttr != null)
                {
                    var types = distsAttr.Value.Split('|');
                    for (int i = 0; i < types.Length; i++)
                    {
                        Enum.TryParse(types[i], out UnivariateDistributionType distType);
                        distributions.Add(UnivariateDistributionFactory.CreateDistribution(distType));
                    }
                }
                var mixture = new Mixture(weights.ToArray(), distributions.ToArray());

                var zeroInflatedAttr = xElement.Attribute(nameof(IsZeroInflated));
                if (zeroInflatedAttr != null)
                {
                    bool.TryParse(zeroInflatedAttr.Value, out var isZeroInflated);
                    mixture.IsZeroInflated = isZeroInflated;
                }
                var zeroWeightAttr = xElement.Attribute(nameof(ZeroWeight));
                if (zeroWeightAttr != null)
                {
                    double.TryParse(zeroWeightAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out var zeroWeight);
                    mixture.ZeroWeight = zeroWeight;
                }
                var xTransformAttr = xElement.Attribute(nameof(XTransform));
                if (xTransformAttr != null)
                {
                    Enum.TryParse(xTransformAttr.Value, out Transform xTransform);
                    mixture.XTransform = xTransform;
                }
                var probTransformAttr = xElement.Attribute(nameof(ProbabilityTransform));
                if (probTransformAttr != null)
                {
                    Enum.TryParse(probTransformAttr.Value, out Transform probabilityTransform);
                    mixture.ProbabilityTransform = probabilityTransform;
                }
                var paramsAttr = xElement.Attribute("Parameters");
                if (paramsAttr != null)
                {
                    var vals = paramsAttr.Value.Split('|');
                    var parameters = new List<double>();
                    for (int i = 0; i < vals.Length; i++)
                    {
                        double.TryParse(vals[i], NumberStyles.Any, CultureInfo.InvariantCulture, out var parm);
                        parameters.Add(parm);
                    }
                    mixture.SetParameters(parameters);
                }

                return mixture;
            }
            else
            {
                return null;
            }
        }

    }
}
