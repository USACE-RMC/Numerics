using Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// A class for a power function with normally distributed noise.
    /// Y = [α * (X - ξ)^β] * ϵ, where ϵ ~ N(0,σ) 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class PowerFunction : IUnivariateFunction
    {

        /// <summary>
        /// Construct a new deterministic power function with α=1 and β=1.5 and ξ=0. 
        /// </summary>
        public PowerFunction()
        {
            Alpha = 1;
            Beta = 1.5;
            Xi = 0;
            IsDeterministic = true;
        }

        /// <summary>
        /// Construct a new power function with a given α, β, and ξ. 
        /// </summary>
        /// <param name="alpha">The coefficient parameter α.</param>
        /// <param name="beta">The exponent parameter β.</param>
        /// <param name="xi">The location parameter ξ.</param>
        public PowerFunction(double alpha, double beta, double xi = 0)
        {
            Alpha = alpha;
            Beta = beta;
            Xi = xi;
            IsDeterministic = true;
        }

        /// <summary>
        /// Construct a new power function with a given α, β, ξ and standard error σ.
        /// </summary>
        /// <param name="alpha">The coefficient parameter α.</param>
        /// <param name="beta">The exponent parameter β.</param>
        /// <param name="xi">The location parameter ξ.</param>
        /// <param name="sigma">The log-space standard error σ.</param>
        public PowerFunction(double alpha, double beta, double xi, double sigma)
        {
            Alpha = alpha;
            Beta = beta;
            Sigma = sigma;
            Xi = xi;
            IsDeterministic = false;
        }

        private bool _parametersValid = true;
        private double _alpha, _beta, _xi, _sigma;
        private Normal _normal = new Normal();

        /// <summary>
        /// The coefficient α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters(new[] { value, Beta, Xi, Sigma }, false) is null;
                _alpha = value;
            }
        }

        /// <summary>
        /// The exponent β (beta).
        /// </summary>
        public double Beta
        {
            get { return _beta; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, value, Xi, Sigma }, false) is null;
                _beta = value;
            }
        }

        /// <summary>
        /// The location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return _xi; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, Beta, value, Sigma }, false) is null;
                _xi = value;
            }
        }

        /// <summary>
        /// The standard error parameter σ (sigma) in log-space.
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, Beta, Xi, value }, false) is null;
                _sigma = value;
                _normal.SetParameters(0, _sigma);
            }
        }

        /// <inheritdoc/>
        public int NumberOfParameters => 4;

        /// <inheritdoc/>
        public bool ParametersValid => _parametersValid;

        /// <inheritdoc/>
        public double Minimum
        {
            get { return Xi; }
            set { throw new NotSupportedException("Minimum is derived from Xi and cannot be set directly."); }
        }

        /// <inheritdoc/>
        public double Maximum { get; set; } = double.MaxValue;

        /// <inheritdoc/>
        public double[] MinimumOfParameters => new double[] { 0, -10, 0, 0 };

        /// <inheritdoc/>
        public double[] MaximumOfParameters => new double[] { double.MaxValue, 10, double.MaxValue, double.MaxValue };

        /// <inheritdoc/>
        public bool IsDeterministic { get; set; }

        /// <summary>
        /// Determines if the power function should be inverted.
        /// </summary>
        public bool IsInverse { get; set; }

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
            _xi = parameters[2];
            _sigma = parameters[3];
            _normal.SetParameters(0, _sigma);
        }

        /// <inheritdoc/>
        public ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (IsDeterministic == false && parameters[3] <= 0)
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
                ValidateParameters(new[] { Alpha, Beta,  Xi, Sigma }, true);

            // Check support
            if (x <= Xi)
                x = Xi + Tools.DoubleMachineEpsilon;
            else if (x >= Maximum)
                x = Maximum;

            double y = 0;
            if (IsInverse)
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                {
                    y = Math.Exp((Math.Log(x) - Math.Log(Alpha)) / Beta) + Xi;
                }
                else
                {
                    y = Math.Exp((Math.Log(x) - Math.Log(Alpha) - _normal.InverseCDF(ConfidenceLevel)) / Beta) + Xi;
                }
            }
            else
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                {
                    y = Math.Exp(Math.Log(Alpha) + Beta * Math.Log(x - Xi));
                }
                else
                {
                    y = Math.Exp(Math.Log(Alpha) + Beta * Math.Log(x - Xi) + _normal.InverseCDF(ConfidenceLevel));
                }
            }
            return y;
        }

        /// <inheritdoc/>
        public double InverseFunction(double y)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Alpha, Beta, Xi, Sigma }, true);

            double x = 0;
            if (IsInverse)
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                {
                    x = Math.Exp(Math.Log(Alpha) + Beta * Math.Log(y - Xi));
                }
                else
                {
                    x = Math.Exp(Math.Log(Alpha) + Beta * Math.Log(y - Xi) + _normal.InverseCDF(ConfidenceLevel));
                }
            }
            else
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                {
                    x = Math.Exp((Math.Log(y) - Math.Log(Alpha)) / Beta) + Xi;
                }
                else
                {
                    x = Math.Exp((Math.Log(y) - Math.Log(Alpha) - _normal.InverseCDF(ConfidenceLevel)) / Beta) + Xi;
                }
            }
            if (x < Minimum) return Minimum;
            if (x > Maximum) return Maximum;
            return x;
        }

        /// <summary>
        /// Serializes the function's configuration to an XElement: the parameters α, β, ξ, and σ,
        /// the deterministic and inverse flags, and the upper support bound.
        /// <see cref="Minimum"/> is derived from ξ and <see cref="ConfidenceLevel"/> is runtime
        /// sampling state — neither is serialized.
        /// </summary>
        /// <returns>An XElement representation of the power function.</returns>
        public XElement ToXElement()
        {
            var result = new XElement(nameof(PowerFunction));
            result.SetAttributeValue(nameof(Alpha), Alpha.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Beta), Beta.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Xi), Xi.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(Sigma), Sigma.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(IsDeterministic), IsDeterministic.ToString());
            result.SetAttributeValue(nameof(IsInverse), IsInverse.ToString());
            result.SetAttributeValue(nameof(Maximum), Maximum.ToString("G17", CultureInfo.InvariantCulture));
            return result;
        }

        /// <summary>
        /// Deserializes a power function from an XElement produced by <see cref="ToXElement"/>.
        /// Missing attributes keep the default-constructed values.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        /// <returns>A new <see cref="PowerFunction"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when a present attribute is malformed or non-finite.</exception>
        public static PowerFunction FromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            var function = new PowerFunction();
            // Set the deterministic flag first: parameter validation is gated on it.
            var deterministicAttribute = xElement.Attribute(nameof(IsDeterministic));
            if (deterministicAttribute != null)
            {
                if (!bool.TryParse(deterministicAttribute.Value, out bool isDeterministic))
                    throw new ArgumentException("The serialized deterministic flag is invalid.", nameof(xElement));
                function.IsDeterministic = isDeterministic;
            }
            var inverseAttribute = xElement.Attribute(nameof(IsInverse));
            if (inverseAttribute != null)
            {
                if (!bool.TryParse(inverseAttribute.Value, out bool isInverse))
                    throw new ArgumentException("The serialized inverse flag is invalid.", nameof(xElement));
                function.IsInverse = isInverse;
            }
            if (TryReadFiniteDouble(xElement, nameof(Alpha), out double alpha))
                function.Alpha = alpha;
            if (TryReadFiniteDouble(xElement, nameof(Beta), out double beta))
                function.Beta = beta;
            if (TryReadFiniteDouble(xElement, nameof(Xi), out double xi))
                function.Xi = xi;
            if (TryReadFiniteDouble(xElement, nameof(Sigma), out double sigma))
                function.Sigma = sigma;
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
