using System;
using System.Globalization;
using System.Linq;
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
        /// <see cref="BivariateCopula.ToXElement()"/>.</param>
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

        /// <summary>
        /// Creates a bivariate copula from its serialized representation, optionally restoring its
        /// marginal distributions.
        /// </summary>
        /// <param name="xElement">The element to deserialize, as written by
        /// <see cref="BivariateCopula.ToXElement(bool)"/>.</param>
        /// <param name="includeMarginals">True to restore the marginal distributions stored in the
        /// "MarginalDistributionX" and "MarginalDistributionY" child elements.</param>
        /// <returns>A validated bivariate copula. A marginal is attached only when
        /// <paramref name="includeMarginals"/> is true and its child element is present.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the type, parameter or marginal data is missing or malformed.</exception>
        /// <exception cref="NotSupportedException">Thrown when a marginal type requires a user-provided implementation.</exception>
        /// <remarks>
        /// Each marginal is rebuilt through <see cref="UnivariateDistributionFactory.CreateDistribution(XElement)"/>.
        /// </remarks>
        public static BivariateCopula CreateCopula(XElement xElement, bool includeMarginals)
        {
            var copula = CreateCopula(xElement);
            if (!includeMarginals) return copula;
            copula.MarginalDistributionX = ReadMarginal(xElement, nameof(BivariateCopula.MarginalDistributionX));
            copula.MarginalDistributionY = ReadMarginal(xElement, nameof(BivariateCopula.MarginalDistributionY));
            return copula;
        }

        /// <summary>
        /// Reads one serialized marginal distribution.
        /// </summary>
        /// <param name="xElement">The copula element.</param>
        /// <param name="name">The name of the child element that holds the marginal.</param>
        /// <returns>The marginal distribution, or null when the child element is absent.</returns>
        /// <exception cref="ArgumentException">Thrown when the child element holds no distribution.</exception>
        private static UnivariateDistributionBase? ReadMarginal(XElement xElement, string name)
        {
            var container = xElement.Element(name);
            if (container == null) return null;
            var payload = container.Elements().FirstOrDefault();
            if (payload == null)
                throw new ArgumentException("The serialized marginal distribution " + name + " is empty.", nameof(xElement));
            return UnivariateDistributionFactory.CreateDistribution(payload);
        }

    }
}
