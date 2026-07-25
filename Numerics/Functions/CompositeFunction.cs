using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// A weighted combination of univariate functions with two modes: the pointwise weighted
    /// average F(x) = Σ wᵢ·fᵢ(x), and a mixture where a single confidence-level draw selects a
    /// child by cumulative weight and re-scales the remainder as the child's own draw.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// This is the shared math behind composite (e.g., day/night exposure) risk input
    /// functions. Weights must be non-negative and sum to one. The mixture mode is driven by
    /// one uniform: with confidence level u, the child i with cumulative weight bracketing u is
    /// selected and evaluated at the re-scaled remainder (u − Σ w₍&lt;i₎)/wᵢ — deterministic
    /// composition sampling with no internal random source, so the same u always reproduces the
    /// same curve. Outside [0, 1] (the mean convention shared by the other function types),
    /// both modes evaluate the weighted average of the children's own mean evaluations.
    /// </para>
    /// <para>
    /// The numeric <see cref="InverseFunction(double)"/> assumes the composed function is
    /// monotone in x (true when every child is monotone non-decreasing, the risk-function use);
    /// the mixture inverse inverts the selected child directly. Serialization embeds each child
    /// through its own serialized form and reconstructs through
    /// <see cref="UnivariateFunctionFactory.CreateFromXElement(XElement)"/>, the
    /// <c>Mixture</c>/<c>CompetingRisks</c> distribution idiom.
    /// </para>
    /// </remarks>
    [Serializable]
    public class CompositeFunction : IUnivariateFunction
    {

        /// <summary>
        /// Construct a new composite function over child functions with equal weights, in the
        /// weighted-average mode.
        /// </summary>
        /// <param name="functions">The child functions.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="functions"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when no child functions are supplied.</exception>
        public CompositeFunction(IList<IUnivariateFunction> functions)
        {
            if (functions == null) throw new ArgumentNullException(nameof(functions));
            if (functions.Count == 0) throw new ArgumentException("At least one child function is required.", nameof(functions));
            _functions = new IUnivariateFunction[functions.Count];
            _weights = new double[functions.Count];
            for (int i = 0; i < functions.Count; i++)
            {
                _functions[i] = functions[i];
                _weights[i] = 1d / functions.Count;
            }
        }

        /// <summary>
        /// Construct a new composite function over child functions and weights, in the
        /// weighted-average mode.
        /// </summary>
        /// <param name="functions">The child functions.</param>
        /// <param name="weights">The non-negative weights; must sum to one.</param>
        /// <exception cref="ArgumentNullException">Thrown when either argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the collections are empty or mismatched in length.</exception>
        public CompositeFunction(IList<IUnivariateFunction> functions, IList<double> weights)
        {
            if (functions == null) throw new ArgumentNullException(nameof(functions));
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            if (functions.Count == 0) throw new ArgumentException("At least one child function is required.", nameof(functions));
            if (functions.Count != weights.Count) throw new ArgumentException("The weight count must match the function count.", nameof(weights));
            _functions = new IUnivariateFunction[functions.Count];
            _weights = new double[functions.Count];
            for (int i = 0; i < functions.Count; i++) _functions[i] = functions[i];
            SetParameters(weights);
        }

        private bool _parametersValid = true;
        private IUnivariateFunction[] _functions;
        private double[] _weights;

        /// <summary>
        /// The combination mode. Default = weighted average.
        /// </summary>
        public CompositeFunctionMode Mode { get; set; } = CompositeFunctionMode.WeightedAverage;

        /// <summary>
        /// The child functions.
        /// </summary>
        public IReadOnlyList<IUnivariateFunction> Functions => _functions;

        /// <summary>
        /// The child weights (non-negative, summing to one).
        /// </summary>
        public IReadOnlyList<double> Weights => _weights;

        /// <inheritdoc/>
        /// <remarks>The parameters are the child weights; children own their own parameters.</remarks>
        public int NumberOfParameters => _functions.Length;

        /// <inheritdoc/>
        /// <remarks>
        /// Reflects the composite's own parameters (the weights) only: children validate their
        /// own parameters at evaluation time, exactly as they do standalone.
        /// </remarks>
        public bool ParametersValid => _parametersValid;

        /// <inheritdoc/>
        /// <remarks>Derived: the smallest child minimum. The setter is not supported.</remarks>
        public double Minimum
        {
            get
            {
                double minimum = double.MaxValue;
                for (int i = 0; i < _functions.Length; i++)
                {
                    if (_functions[i].Minimum < minimum) minimum = _functions[i].Minimum;
                }
                return minimum;
            }
            set { throw new NotSupportedException("Minimum is derived from the child functions and cannot be set directly."); }
        }

        /// <inheritdoc/>
        /// <remarks>Derived: the largest child maximum. The setter is not supported.</remarks>
        public double Maximum
        {
            get
            {
                double maximum = double.MinValue;
                for (int i = 0; i < _functions.Length; i++)
                {
                    if (_functions[i].Maximum > maximum) maximum = _functions[i].Maximum;
                }
                return maximum;
            }
            set { throw new NotSupportedException("Maximum is derived from the child functions and cannot be set directly."); }
        }

        /// <inheritdoc/>
        public double[] MinimumOfParameters
        {
            get
            {
                var result = new double[_functions.Length];
                for (int i = 0; i < result.Length; i++) result[i] = 0d;
                return result;
            }
        }

        /// <inheritdoc/>
        public double[] MaximumOfParameters
        {
            get
            {
                var result = new double[_functions.Length];
                for (int i = 0; i < result.Length; i++) result[i] = 1d;
                return result;
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// True only when every child is deterministic and, in the mixture mode, when at most
        /// one child carries positive weight (a multi-branch mixture varies with the draw even
        /// over deterministic children). The setter propagates the flag to every child.
        /// </remarks>
        public bool IsDeterministic
        {
            get
            {
                for (int i = 0; i < _functions.Length; i++)
                {
                    if (!_functions[i].IsDeterministic) return false;
                }
                if (Mode == CompositeFunctionMode.Mixture)
                {
                    int reachable = 0;
                    for (int i = 0; i < _weights.Length; i++)
                    {
                        if (_weights[i] > 0d) reachable++;
                    }
                    if (reachable > 1) return false;
                }
                return true;
            }
            set
            {
                for (int i = 0; i < _functions.Length; i++) _functions[i].IsDeterministic = value;
            }
        }

        /// <inheritdoc/>
        public double ConfidenceLevel { get; set; } = -1;

        /// <inheritdoc/>
        /// <remarks>The parameters are the child weights.</remarks>
        public void SetParameters(IList<double> parameters)
        {
            // Validate parameters
            _parametersValid = ValidateParameters(parameters, false) is null;
            // Set parameters
            for (int i = 0; i < _weights.Length && i < parameters.Count; i++)
                _weights[i] = parameters[i];
        }

        /// <inheritdoc/>
        public ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters == null || parameters.Count != _functions.Length)
            {
                var error = new ArgumentOutOfRangeException(nameof(parameters), "The weight count must match the child function count.");
                if (throwException) throw error;
                return error;
            }
            double sum = 0d;
            for (int i = 0; i < parameters.Count; i++)
            {
                if (!Tools.IsFinite(parameters[i]) || parameters[i] < 0d)
                {
                    var error = new ArgumentOutOfRangeException(nameof(parameters), "Weights must be non-negative and finite.");
                    if (throwException) throw error;
                    return error;
                }
                sum += parameters[i];
            }
            if (Math.Abs(sum - 1d) > 1E-8)
            {
                var error = new ArgumentOutOfRangeException(nameof(parameters), "Weights must sum to one.");
                if (throwException) throw error;
                return error;
            }
            return null!;
        }

        /// <inheritdoc/>
        public double Function(double x)
        {
            // Validate parameters
            if (ParametersValid == false)
                ValidateParameters(_weights, true);

            if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
            {
                // The mean convention: the weighted average of the children's own mean
                // evaluations (children left at their configured state).
                double mean = 0d;
                for (int i = 0; i < _functions.Length; i++)
                    mean += _weights[i] * _functions[i].Function(x);
                return mean;
            }

            if (Mode == CompositeFunctionMode.Mixture)
            {
                int index = SelectMixtureChild(ConfidenceLevel, out double remainder);
                return EvaluateChildAt(index, remainder, x);
            }

            // Weighted average: the composite's draw drives every child co-monotonically.
            double result = 0d;
            for (int i = 0; i < _functions.Length; i++)
                result += _weights[i] * EvaluateChildAt(i, ConfidenceLevel, x);
            return result;
        }

        /// <inheritdoc/>
        public double InverseFunction(double y)
        {
            // Validate parameters
            if (ParametersValid == false)
                ValidateParameters(_weights, true);

            if (Mode == CompositeFunctionMode.Mixture && IsDeterministic == false && ConfidenceLevel >= 0 && ConfidenceLevel <= 1)
            {
                int index = SelectMixtureChild(ConfidenceLevel, out double remainder);
                return InvertChildAt(index, remainder, y);
            }

            // Weighted average (or the mean convention): numeric monotone bracketing over the
            // composed forward evaluation.
            double lower = Minimum;
            double upper = Maximum;
            if (double.IsInfinity(lower) || lower == double.MinValue) lower = -1d;
            if (double.IsInfinity(upper) || upper == double.MaxValue) upper = 1d;
            while (Function(lower) > y && Tools.IsFinite(lower)) lower = lower < 0 ? lower * 2d : (lower - 1d) * 2d;
            while (Function(upper) < y && Tools.IsFinite(upper)) upper = upper > 0 ? upper * 2d : (upper + 1d) * 2d;
            return Mathematics.RootFinding.Brent.Solve(h => Function(h) - y, lower, upper);
        }

        /// <summary>
        /// Serializes the composite's configuration to an XElement: the mode, the weights, and
        /// each child's own serialized form. <see cref="ConfidenceLevel"/> is runtime sampling
        /// state and is deliberately not serialized.
        /// </summary>
        /// <returns>An XElement representation of the composite function.</returns>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(CompositeFunction));
            result.SetAttributeValue(nameof(Mode), Mode.ToString());
            var values = new string[_weights.Length];
            for (int i = 0; i < _weights.Length; i++)
                values[i] = _weights[i].ToString("G17", CultureInfo.InvariantCulture);
            result.SetAttributeValue(nameof(Weights), string.Join("|", values));
            var children = new XElement("Functions");
            for (int i = 0; i < _functions.Length; i++)
            {
                children.Add(SerializeChild(_functions[i]));
            }
            result.Add(children);
            return result;
        }

        /// <summary>
        /// Deserializes a composite function from an XElement produced by
        /// <see cref="ToXElement"/>; children reconstruct through the function factory.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new <see cref="CompositeFunction"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the element carries no child functions or mismatched weights.</exception>
        public static CompositeFunction FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            var childrenElement = xElement.Element("Functions");
            if (childrenElement == null)
                throw new ArgumentException("The serialized composite function is missing its child functions.", nameof(xElement));

            var functions = new List<IUnivariateFunction>();
            foreach (var child in childrenElement.Elements())
            {
                functions.Add(UnivariateFunctionFactory.CreateFromXElement(child));
            }
            if (functions.Count == 0)
                throw new ArgumentException("The serialized composite function carries no child functions.", nameof(xElement));

            var composite = new CompositeFunction(functions);
            string? weightsText = xElement.Attribute(nameof(Weights))?.Value;
            if (!string.IsNullOrEmpty(weightsText))
            {
                string[] tokens = weightsText!.Split('|');
                if (tokens.Length != functions.Count)
                    throw new ArgumentException("The serialized composite function's weight count does not match its child count.", nameof(xElement));
                var weights = new double[tokens.Length];
                for (int i = 0; i < tokens.Length; i++)
                {
                    if (!double.TryParse(tokens[i], NumberStyles.Any, CultureInfo.InvariantCulture, out weights[i]))
                        throw new ArgumentException("The serialized composite function carries an unparseable weight.", nameof(xElement));
                }
                composite.SetParameters(weights);
            }
            if (Enum.TryParse(xElement.Attribute(nameof(Mode))?.Value, out CompositeFunctionMode mode))
                composite.Mode = mode;
            return composite;
        }

        /// <summary>
        /// Selects the mixture child whose cumulative weight brackets the draw, and re-scales
        /// the remainder as the child's own confidence level.
        /// </summary>
        /// <param name="u">The composite's uniform draw in [0, 1].</param>
        /// <param name="remainder">The re-scaled child draw.</param>
        /// <returns>The selected child index.</returns>
        private int SelectMixtureChild(double u, out double remainder)
        {
            double cumulative = 0d;
            for (int i = 0; i < _weights.Length; i++)
            {
                if (_weights[i] <= 0d) continue;
                if (u <= cumulative + _weights[i] || i == _weights.Length - 1)
                {
                    remainder = (u - cumulative) / _weights[i];
                    if (remainder < 0d) remainder = 0d;
                    if (remainder > 1d) remainder = 1d;
                    return i;
                }
                cumulative += _weights[i];
            }
            // All-zero weights cannot validate; fall back to the last child defensively.
            remainder = u;
            return _weights.Length - 1;
        }

        /// <summary>
        /// Evaluates a child at a given confidence level without disturbing the child's
        /// configured state (the level is restored after the call).
        /// </summary>
        /// <param name="index">The child index.</param>
        /// <param name="confidenceLevel">The confidence level to evaluate at.</param>
        /// <param name="x">The evaluation point.</param>
        /// <returns>The child's value.</returns>
        private double EvaluateChildAt(int index, double confidenceLevel, double x)
        {
            var child = _functions[index];
            double restore = child.ConfidenceLevel;
            child.ConfidenceLevel = confidenceLevel;
            double value = child.Function(x);
            child.ConfidenceLevel = restore;
            return value;
        }

        /// <summary>
        /// Inverts a child at a given confidence level without disturbing the child's
        /// configured state (the level is restored after the call).
        /// </summary>
        /// <param name="index">The child index.</param>
        /// <param name="confidenceLevel">The confidence level to invert at.</param>
        /// <param name="y">The value to invert.</param>
        /// <returns>The child's inverse.</returns>
        private double InvertChildAt(int index, double confidenceLevel, double y)
        {
            var child = _functions[index];
            double restore = child.ConfidenceLevel;
            child.ConfidenceLevel = confidenceLevel;
            double value = child.InverseFunction(y);
            child.ConfidenceLevel = restore;
            return value;
        }

        /// <summary>
        /// Serializes a child through its concrete ToXElement method.
        /// </summary>
        /// <param name="function">The child function.</param>
        /// <returns>The child's serialized form.</returns>
        /// <exception cref="NotSupportedException">Thrown when the child is not a library function type.</exception>
        private static XElement SerializeChild(IUnivariateFunction function)
        {
            if (function is LinearFunction linear) return linear.ToXElement();
            if (function is PowerFunction power) return power.ToXElement();
            if (function is TabularFunction tabular) return tabular.ToXElement();
            if (function is SegmentedPowerFunction segmented) return segmented.ToXElement();
            if (function is CompositeFunction composite) return composite.ToXElement();
            throw new NotSupportedException("The child function type '" + function.GetType().Name + "' does not support serialization.");
        }

    }
}
