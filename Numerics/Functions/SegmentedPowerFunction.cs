using Numerics.Distributions;
using Numerics.Mathematics.RootFinding;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// A segmented power function in the BaRatin matrix-of-controls ADDITION mode:
    /// Q(h) = Σₖ 10^(log₁₀αₖ) · (h − hₖ)^βₖ · 𝟙{h &gt; hₖ}, with a log₁₀-space Gaussian residual σ.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// The parameter-vector layout is exactly the RMC-BestFit rating-curve layout
    /// (<c>RMC.BestFit/Models/RatingCurve/RatingCurve.cs</c>), so a fitted posterior parameter
    /// set applies directly through <see cref="SetParameters(IList{double})"/>:
    /// <c>[h₁, log₁₀α₁, β₁, h₂, log₁₀α₂, β₂, …, σ]</c> with length 3·segments + 1. Under
    /// addition mode the BaRatin continuity derivation collapses each control's offset to its
    /// activation stage (b_k = κ_k), the breakpoints must be strictly ordered
    /// (h₁ &lt; h₂ &lt; …), and discharge is zero at and below the main-channel cease-to-flow
    /// stage h₁. One segment degenerates to the plain power law
    /// <see cref="PowerFunction"/> with α = 10^(log₁₀α₁), β = β₁, ξ = h₁.
    /// </para>
    /// <para>
    /// The residual is log₁₀-space, matching the BestFit stochastic prediction: with a
    /// confidence level u, Q(h) is multiplied by 10^(z) where z ~ N(0, σ) evaluated at u. The
    /// numeric <see cref="InverseFunction(double)"/> uses monotone bracketing with Brent's
    /// method — the addition-mode sum is strictly increasing above h₁ for non-negative
    /// exponents.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// Le Coz, J., Renard, B., Bonnifait, L., Branger, F., Le Boursicaud, R. (2014). Combining
    /// hydraulic knowledge and uncertain gaugings in the estimation of hydrometric rating
    /// curves: a Bayesian approach. J. Hydrol. 509:573–587.
    /// </para>
    /// </remarks>
    [Serializable]
    public class SegmentedPowerFunction : IUnivariateFunction
    {

        /// <summary>
        /// Construct a new deterministic single-segment power function with
        /// h₁ = 0, log₁₀α₁ = 0, β₁ = 1.5, and σ = 0.1.
        /// </summary>
        public SegmentedPowerFunction() : this(1)
        {
        }

        /// <summary>
        /// Construct a new deterministic segmented power function with default, strictly
        /// ordered breakpoints (hₖ = k − 1), log₁₀αₖ = 0, βₖ = 1.5, and σ = 0.1.
        /// </summary>
        /// <param name="numberOfSegments">The number of segments (≥ 1).</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="numberOfSegments"/> is less than one.</exception>
        public SegmentedPowerFunction(int numberOfSegments)
        {
            if (numberOfSegments < 1)
                throw new ArgumentOutOfRangeException(nameof(numberOfSegments), "The number of segments must be at least one.");
            _numberOfSegments = numberOfSegments;
            _parameters = new double[3 * numberOfSegments + 1];
            for (int k = 0; k < numberOfSegments; k++)
            {
                _parameters[3 * k] = k;
                _parameters[3 * k + 1] = 0d;
                _parameters[3 * k + 2] = 1.5d;
            }
            _parameters[_parameters.Length - 1] = 0.1d;
            _normal.SetParameters(0d, _parameters[_parameters.Length - 1]);
            IsDeterministic = true;
        }

        /// <summary>
        /// Construct a new segmented power function directly from a BestFit-layout parameter
        /// vector. The segment count is inferred from the vector length (3·segments + 1), and
        /// the function is stochastic (σ is the last entry).
        /// </summary>
        /// <param name="parameters">The parameter vector [h₁, log₁₀α₁, β₁, …, σ].</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="parameters"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the vector length is not 3·segments + 1 for a positive segment count.</exception>
        public SegmentedPowerFunction(IList<double> parameters)
        {
            if (parameters == null) throw new ArgumentNullException(nameof(parameters));
            if (parameters.Count < 4 || (parameters.Count - 1) % 3 != 0)
                throw new ArgumentException("The parameter vector length must be 3·segments + 1.", nameof(parameters));
            _numberOfSegments = (parameters.Count - 1) / 3;
            _parameters = new double[parameters.Count];
            IsDeterministic = false;
            SetParameters(parameters);
        }

        private bool _parametersValid = true;
        private int _numberOfSegments;
        private double[] _parameters;
        private Normal _normal = new Normal();

        /// <summary>
        /// The number of power-law segments (controls). Fixed at construction so the parameter
        /// vector length stays coherent.
        /// </summary>
        public int NumberOfSegments => _numberOfSegments;

        /// <summary>
        /// The log₁₀-space standard error σ (the last parameter-vector entry).
        /// </summary>
        public double Sigma => _parameters[_parameters.Length - 1];

        /// <inheritdoc/>
        public int NumberOfParameters => 3 * _numberOfSegments + 1;

        /// <inheritdoc/>
        public bool ParametersValid => _parametersValid;

        /// <inheritdoc/>
        public double Minimum
        {
            get { return _parameters[0]; }
            set { throw new NotSupportedException("Minimum is derived from the first breakpoint h₁ and cannot be set directly."); }
        }

        /// <inheritdoc/>
        public double Maximum { get; set; } = double.MaxValue;

        /// <inheritdoc/>
        public double[] MinimumOfParameters
        {
            get
            {
                var result = new double[NumberOfParameters];
                for (int k = 0; k < _numberOfSegments; k++)
                {
                    result[3 * k] = double.MinValue;
                    result[3 * k + 1] = double.MinValue;
                    result[3 * k + 2] = 0d;
                }
                result[result.Length - 1] = 0d;
                return result;
            }
        }

        /// <inheritdoc/>
        public double[] MaximumOfParameters
        {
            get
            {
                var result = new double[NumberOfParameters];
                for (int i = 0; i < result.Length; i++) result[i] = double.MaxValue;
                return result;
            }
        }

        /// <inheritdoc/>
        public bool IsDeterministic { get; set; }

        /// <inheritdoc/>
        public double ConfidenceLevel { get; set; } = -1;

        /// <summary>
        /// Gets the breakpoint (activation stage) hₖ of a one-based segment.
        /// </summary>
        /// <param name="segmentOneBased">The segment index (1..NumberOfSegments).</param>
        /// <returns>The breakpoint hₖ.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the segment index is out of range.</exception>
        public double GetBreakpoint(int segmentOneBased)
        {
            ValidateSegmentIndex(segmentOneBased);
            return _parameters[3 * (segmentOneBased - 1)];
        }

        /// <summary>
        /// Gets log₁₀(αₖ) of a one-based segment.
        /// </summary>
        /// <param name="segmentOneBased">The segment index (1..NumberOfSegments).</param>
        /// <returns>log₁₀(αₖ).</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the segment index is out of range.</exception>
        public double GetLog10Alpha(int segmentOneBased)
        {
            ValidateSegmentIndex(segmentOneBased);
            return _parameters[3 * (segmentOneBased - 1) + 1];
        }

        /// <summary>
        /// Gets the exponent βₖ of a one-based segment.
        /// </summary>
        /// <param name="segmentOneBased">The segment index (1..NumberOfSegments).</param>
        /// <returns>βₖ.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the segment index is out of range.</exception>
        public double GetBeta(int segmentOneBased)
        {
            ValidateSegmentIndex(segmentOneBased);
            return _parameters[3 * (segmentOneBased - 1) + 2];
        }

        /// <inheritdoc/>
        public void SetParameters(IList<double> parameters)
        {
            // Validate parameters
            _parametersValid = ValidateParameters(parameters, false) is null;
            // Set parameters
            for (int i = 0; i < _parameters.Length && i < parameters.Count; i++)
                _parameters[i] = parameters[i];
            _normal.SetParameters(0d, _parameters[_parameters.Length - 1]);
        }

        /// <inheritdoc/>
        public ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
            {
                var error = new ArgumentOutOfRangeException(nameof(parameters), "The parameter vector length must be 3·segments + 1.");
                if (throwException) throw error;
                return error;
            }
            for (int i = 0; i < parameters.Count; i++)
            {
                if (!Tools.IsFinite(parameters[i]))
                {
                    var error = new ArgumentOutOfRangeException(nameof(parameters), "All parameters must be finite.");
                    if (throwException) throw error;
                    return error;
                }
            }
            for (int k = 0; k < (parameters.Count - 1) / 3; k++)
            {
                if (parameters[3 * k + 2] < 0)
                {
                    var error = new ArgumentOutOfRangeException(nameof(parameters), "Exponents must be non-negative for a monotone rating.");
                    if (throwException) throw error;
                    return error;
                }
                if (k > 0 && parameters[3 * k] <= parameters[3 * (k - 1)])
                {
                    var error = new ArgumentOutOfRangeException(nameof(parameters), "Breakpoints must be strictly increasing (h₁ < h₂ < …).");
                    if (throwException) throw error;
                    return error;
                }
            }
            if (IsDeterministic == false && parameters[parameters.Count - 1] <= 0)
            {
                var error = new ArgumentOutOfRangeException(nameof(Sigma), "Standard error must be greater than zero.");
                if (throwException) throw error;
                return error;
            }
            return null!;
        }

        /// <inheritdoc/>
        public double Function(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_parameters, true);

            // Check support
            if (x >= Maximum) x = Maximum;
            double q = DeterministicFunction(x);
            if (q <= 0d) return 0d;

            if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                return q;
            // Log₁₀-space residual, matching the BestFit stochastic prediction.
            return q * Math.Pow(10d, _normal.InverseCDF(ConfidenceLevel));
        }

        /// <inheritdoc/>
        public double InverseFunction(double y)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_parameters, true);

            // Fold the residual out first: the stochastic curve is the deterministic curve
            // scaled by 10^z, so the inverse divides before the monotone root find.
            if (IsDeterministic == false && ConfidenceLevel >= 0 && ConfidenceLevel <= 1)
                y /= Math.Pow(10d, _normal.InverseCDF(ConfidenceLevel));

            double h1 = _parameters[0];
            if (y <= 0d) return h1;

            // Bracket the root: the addition-mode sum is strictly increasing above h₁, so
            // expand the upper bound geometrically until it covers y (or the support cap).
            double lower = h1;
            double upper = _parameters[3 * (_numberOfSegments - 1)] + 1d;
            if (upper <= lower) upper = lower + 1d;
            while (DeterministicFunction(upper) < y)
            {
                if (upper >= Maximum) return Maximum;
                double span = upper - h1;
                upper = h1 + span * 2d;
                if (upper > Maximum) upper = Maximum;
            }

            double x = Brent.Solve(h => DeterministicFunction(h) - y, lower, upper);
            if (x < h1) return h1;
            if (x > Maximum) return Maximum;
            return x;
        }

        /// <summary>
        /// Serializes the function's configuration to an XElement: the segment count, the
        /// BestFit-layout parameter vector, the deterministic flag, and the upper support bound.
        /// <see cref="Minimum"/> derives from h₁ and <see cref="ConfidenceLevel"/> is runtime
        /// sampling state — neither is serialized.
        /// </summary>
        /// <returns>An XElement representation of the segmented power function.</returns>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(SegmentedPowerFunction));
            result.SetAttributeValue(nameof(NumberOfSegments), NumberOfSegments.ToString(CultureInfo.InvariantCulture));
            var values = new string[_parameters.Length];
            for (int i = 0; i < _parameters.Length; i++)
                values[i] = _parameters[i].ToString("G17", CultureInfo.InvariantCulture);
            result.SetAttributeValue("Parameters", string.Join("|", values));
            result.SetAttributeValue(nameof(IsDeterministic), IsDeterministic.ToString());
            result.SetAttributeValue(nameof(Maximum), Maximum.ToString("G17", CultureInfo.InvariantCulture));
            return result;
        }

        /// <summary>
        /// Deserializes a segmented power function from an XElement produced by
        /// <see cref="ToXElement"/>.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new <see cref="SegmentedPowerFunction"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the element carries no parseable parameter vector.</exception>
        public static SegmentedPowerFunction FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            string? text = xElement.Attribute("Parameters")?.Value;
            if (string.IsNullOrEmpty(text))
                throw new ArgumentException("The serialized segmented power function is missing its parameter vector.", nameof(xElement));

            string[] tokens = text!.Split('|');
            var parameters = new double[tokens.Length];
            for (int i = 0; i < tokens.Length; i++)
            {
                if (!double.TryParse(tokens[i], NumberStyles.Any, CultureInfo.InvariantCulture, out parameters[i]))
                    throw new ArgumentException("The serialized segmented power function carries an unparseable parameter value.", nameof(xElement));
            }

            var function = new SegmentedPowerFunction(parameters);
            if (bool.TryParse(xElement.Attribute(nameof(IsDeterministic))?.Value, out bool isDeterministic))
                function.IsDeterministic = isDeterministic;
            if (double.TryParse(xElement.Attribute(nameof(Maximum))?.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out double maximum))
                function.Maximum = maximum;
            return function;
        }

        /// <summary>
        /// The deterministic addition-mode discharge: the sum of the active controls' power
        /// laws, zero at and below the main-channel cease-to-flow stage h₁ (the BestFit
        /// <c>Predict</c> body).
        /// </summary>
        /// <param name="x">The stage.</param>
        /// <returns>The deterministic discharge.</returns>
        private double DeterministicFunction(double x)
        {
            double depth1 = x - _parameters[0];
            if (depth1 <= 0d) return 0d;

            double q = Math.Pow(10d, _parameters[1]) * Math.Pow(depth1, _parameters[2]);
            for (int k = 1; k < _numberOfSegments; k++)
            {
                double depth = x - _parameters[3 * k];
                if (depth > 0d)
                    q += Math.Pow(10d, _parameters[3 * k + 1]) * Math.Pow(depth, _parameters[3 * k + 2]);
            }
            return q;
        }

        /// <summary>
        /// Throws when a one-based segment index is out of range.
        /// </summary>
        /// <param name="segmentOneBased">The segment index.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the segment index is out of range.</exception>
        private void ValidateSegmentIndex(int segmentOneBased)
        {
            if (segmentOneBased < 1 || segmentOneBased > _numberOfSegments)
                throw new ArgumentOutOfRangeException(nameof(segmentOneBased), "Segment index must be between 1 and " + _numberOfSegments + ".");
        }

    }
}
