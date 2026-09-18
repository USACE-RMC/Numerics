using Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;

namespace Numerics.Functions
{

    /// <summary>
    /// A class for a simple linear function, with a single predictor and a slope and intercept coefficient, and normally distributed noise.
    /// Y = α + βX + ϵ, where ϵ ~ N(0,σ) 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class LinearFunction : IUnivariateFunction
    {

        /// <summary>
        /// Construct a new linear function with an intercept of 0 and slope of 1. 
        /// </summary>
        public LinearFunction()
        {
            Alpha = 0;
            Beta = 1;
            IsDeterministic = true;
        }

        /// <summary>
        /// Construct a new linear function with a given intercept and slope. 
        /// </summary>
        /// <param name="alpha">The intercept parameter.</param>
        /// <param name="beta">The slope parameter.</param>
        public LinearFunction(double alpha, double beta)
        {
            Alpha = alpha;
            Beta = beta;
            IsDeterministic = true;
        }

        /// <summary>
        /// Construct a new linear function with a given intercept, slope and standard error.
        /// </summary>
        /// <param name="alpha">The intercept parameter.</param>
        /// <param name="beta">The slope parameter.</param>
        /// <param name="sigma">The standard error.</param>
        public LinearFunction(double alpha, double beta, double sigma)
        {
            Alpha = alpha;
            Beta = beta;
            Sigma = sigma;
            IsDeterministic = false;
        }

        private bool _parametersValid = true;
        private double _alpha, _beta, _sigma;
        private Normal _normal = new Normal();

        /// <summary>
        /// The intercept parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters(new[] { value, Beta, Sigma }, false) is null;
                _alpha = value;
            }
        }

        /// <summary>
        /// The slope parameter β (beta).
        /// </summary>
        public double Beta
        {
            get { return _beta; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, value , Sigma}, false) is null;
                _beta = value;
            }
        }

        /// <summary>
        /// The standard error parameter σ (sigma).
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, Beta, value}, false) is null;
                _sigma = value;
                _normal.SetParameters(0, _sigma);
            }
        }

        /// <inheritdoc/>
        public int NumberOfParameters => 3;

        /// <inheritdoc/>
        public bool ParametersValid => _parametersValid;

        /// <inheritdoc/>
        public double Minimum { get; set; } = double.MinValue;

        /// <inheritdoc/>
        public double Maximum { get; set; } = double.MaxValue;

        /// <inheritdoc/>
        public double[] MinimumOfParameters => new double[] { double.MinValue, double.MinValue, 0 };

        /// <inheritdoc/>
        public double[] MaximumOfParameters => new double[] { double.MaxValue, double.MaxValue, double.MaxValue };

        /// <inheritdoc/>
        public bool IsDeterministic { get; set; }

        /// <inheritdoc/>
        public double ConfidenceLevel { get; set; } = -1;

        /// <inheritdoc/>
        public void SetParameters(IList<double> parameters)
        {
            // Validate parameters
            _parametersValid = ValidateParameters(parameters, false) is null;
            // Set parameters
            _alpha = parameters[0];
            _beta = parameters[1];
            _sigma = parameters[2];
            _normal.SetParameters(0, _sigma);
        }

        /// <inheritdoc/>
        public ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (IsDeterministic == false && parameters[2] <= 0)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Standard error must be greater than zero.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Standard error must be greater than zero.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public double Function(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Alpha, Beta, Sigma }, true);
        
            // Check support
            if (x <= Minimum)
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1) return Alpha + Beta * Minimum;
                return Alpha + Beta * Minimum + _normal.InverseCDF(ConfidenceLevel);
            }
            if (x >= Maximum)
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1) return Alpha + Beta * Maximum;
                return Alpha + Beta * Maximum + _normal.InverseCDF(ConfidenceLevel);
            }

            if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1) return Alpha + Beta * x;
            return Alpha + Beta * x + _normal.InverseCDF(ConfidenceLevel);
        }

        /// <inheritdoc/>
        public double InverseFunction(double y)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Alpha, Beta, Sigma }, true);
            if (Math.Abs(Beta) < double.Epsilon)
                throw new InvalidOperationException("Cannot compute inverse function when Beta is zero.");

            double x = 0;
            if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
            {
                x = (y - Alpha) / Beta;
            }
            else
            {
                x = (y - Alpha - _normal.InverseCDF(ConfidenceLevel)) / Beta;
            }
            if (x < Minimum) return Minimum;
            if (x > Maximum) return Maximum;
            return x;
        }

        /// <summary>
        /// Serializes the function's configuration to an XElement: the parameters α, β, and σ,
        /// the deterministic flag, and the support bounds. <see cref="ConfidenceLevel"/> is
        /// runtime sampling state and is deliberately not serialized.
        /// </summary>
        /// <returns>An XElement representation of the linear function.</returns>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(LinearFunction));
            result.SetAttributeValue(nameof(Alpha), Alpha.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Beta), Beta.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Sigma), Sigma.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(IsDeterministic), IsDeterministic.ToString());
            result.SetAttributeValue(nameof(Minimum), Minimum.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Maximum), Maximum.ToString("G17", CultureInfo.InvariantCulture));
            return result;
        }

        /// <summary>
        /// Deserializes a linear function from an XElement produced by <see cref="ToXElement"/>.
        /// Missing attributes keep the default-constructed values.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new <see cref="LinearFunction"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when a present attribute is malformed or non-finite.</exception>
        public static LinearFunction FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            var function = new LinearFunction();
            // Set the deterministic flag first: parameter validation is gated on it.
            var deterministicAttribute = xElement.Attribute(nameof(IsDeterministic));
            if (deterministicAttribute != null)
            {
                if (!bool.TryParse(deterministicAttribute.Value, out bool isDeterministic))
                    throw new ArgumentException("The serialized deterministic flag is invalid.", nameof(xElement));
                function.IsDeterministic = isDeterministic;
            }
            if (TryReadFiniteDouble(xElement, nameof(Alpha), out double alpha))
                function.Alpha = alpha;
            if (TryReadFiniteDouble(xElement, nameof(Beta), out double beta))
                function.Beta = beta;
            if (TryReadFiniteDouble(xElement, nameof(Sigma), out double sigma))
                function.Sigma = sigma;
            if (TryReadFiniteDouble(xElement, nameof(Minimum), out double minimum))
                function.Minimum = minimum;
            if (TryReadFiniteDouble(xElement, nameof(Maximum), out double maximum))
                function.Maximum = maximum;
            return function;
        }

        /// <summary>Reads an optional finite double attribute.</summary>
        private static bool TryReadFiniteDouble(XElement xElement, string attributeName, out double value)
        {
            var attribute = xElement.Attribute(attributeName);
            if (attribute == null)
            {
                value = 0d;
                return false;
            }
            if (!double.TryParse(attribute.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out value)
                || !Tools.IsFinite(value))
            {
                throw new ArgumentException("The serialized " + attributeName + " value is invalid.", nameof(xElement));
            }
            return true;
        }

    }
}
