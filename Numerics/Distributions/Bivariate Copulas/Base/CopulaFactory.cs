using System;
using System.Globalization;
using System.Xml.Linq;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// A bivariate copula factory class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public sealed class CopulaFactory
    {

        /// <summary>
        /// Create a bivariate copula based on the copula type, with its default parameters.
        /// </summary>
        /// <param name="copulaType">Copula type.</param>
        /// <returns>
        /// A bivariate copula.
        /// </returns>
        /// <exception cref="NotSupportedException">
        /// <paramref name="copulaType"/> is not a defined <see cref="CopulaType"/> value.
        /// </exception>
        public static BivariateCopula CreateCopula(CopulaType copulaType)
        {
            switch (copulaType)
            {
                case CopulaType.AliMikhailHaq:
                    return new AMHCopula();
                case CopulaType.Clayton:
                    return new ClaytonCopula();
                case CopulaType.Frank:
                    return new FrankCopula();
                case CopulaType.Gumbel:
                    return new GumbelCopula();
                case CopulaType.Joe:
                    return new JoeCopula();
                case CopulaType.Normal:
                    return new NormalCopula();
                case CopulaType.StudentT:
                    return new StudentTCopula();
                case CopulaType.Independence:
                    return new IndependenceCopula();
                default:
                    throw new NotSupportedException($"The copula type {copulaType} is not supported.");
            }
        }

        /// <summary>
        /// Creates a bivariate copula from its serialized representation.
        /// </summary>
        /// <param name="xElement">The element to deserialize, as written by
        /// <see cref="BivariateCopula.ToXElement"/>.</param>
        /// <returns>A validated bivariate copula. Marginal distributions are not part of the
        /// serialized form and are left unattached.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the type or parameter data is missing or malformed.</exception>
        public static BivariateCopula CreateCopula(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));

            var typeAttribute = xElement.Attribute(nameof(BivariateCopula.Type));
            if (typeAttribute == null
                || !Enum.TryParse(typeAttribute.Value, out CopulaType type)
                || !Enum.IsDefined(typeof(CopulaType), type))
                throw new ArgumentException("The serialized copula type is missing or invalid.", nameof(xElement));

            var copula = CreateCopula(type);

            // A missing or empty "Parameters" attribute carries zero values, so the count
            // check below accepts it for a zero-parameter copula and rejects it otherwise.
            var parametersAttribute = xElement.Attribute("Parameters");
            string parameterText = parametersAttribute == null ? string.Empty : parametersAttribute.Value;
            string[] parts = parameterText.Length == 0 ? Array.Empty<string>() : parameterText.Split('|');
            if (parts.Length != copula.NumberOfCopulaParameters)
                throw new ArgumentException("The serialized copula parameter count does not match the copula type.", nameof(xElement));

            if (parts.Length > 0)
            {
                var values = new double[parts.Length];
                for (int i = 0; i < parts.Length; i++)
                {
                    if (!double.TryParse(parts[i], NumberStyles.Any, CultureInfo.InvariantCulture, out values[i])
                        || !Tools.IsFinite(values[i]))
                        throw new ArgumentException("The serialized copula parameter at position " + i + " is missing or invalid.", nameof(xElement));
                }
                copula.SetCopulaParameters(values);
                if (!copula.ParametersValid)
                    throw new ArgumentException("The serialized parameters do not define a valid copula.", nameof(xElement));
            }
            return copula;
        }

    }
}
