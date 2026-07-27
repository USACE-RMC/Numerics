using System;
using System.Xml.Linq;

namespace Numerics.Distributions
{

    /// <summary>
    /// A univariate distribution factory class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public sealed class UnivariateDistributionFactory
    {

        /// <summary>
        /// Create a distribution based on the distribution type.
        /// </summary>
        /// <param name="distributionType">Distribution type.</param>
        /// <returns>
        /// A univariate distribution.
        /// </returns>
        /// <exception cref="NotSupportedException">
        /// The requested distribution type requires additional component distributions or a user implementation.
        /// </exception>
        /// <exception cref="ArgumentOutOfRangeException">
        /// <paramref name="distributionType"/> is not a defined <see cref="UnivariateDistributionType"/> value.
        /// </exception>
        public static UnivariateDistributionBase CreateDistribution(UnivariateDistributionType distributionType)
        {
            switch (distributionType)
            {
                case UnivariateDistributionType.ChiSquared:
                    return new ChiSquared();
                case UnivariateDistributionType.Bernoulli:
                    return new Bernoulli();
                case UnivariateDistributionType.Beta:
                    return new BetaDistribution();
                case UnivariateDistributionType.Binomial:
                    return new Binomial();
                case UnivariateDistributionType.Cauchy:
                    return new Cauchy();
                case UnivariateDistributionType.Deterministic:
                    return new Deterministic();
                case UnivariateDistributionType.Empirical:
                    return new EmpiricalDistribution();
                case UnivariateDistributionType.Exponential:
                    return new Exponential();
                case UnivariateDistributionType.GammaDistribution:
                    return new GammaDistribution();
                case UnivariateDistributionType.GeneralizedBeta:
                    return new GeneralizedBeta();
                case UnivariateDistributionType.GeneralizedExtremeValue:
                    return new GeneralizedExtremeValue();
                case UnivariateDistributionType.GeneralizedLogistic:
                    return new GeneralizedLogistic();
                case UnivariateDistributionType.GeneralizedNormal:
                    return new GeneralizedNormal();
                case UnivariateDistributionType.GeneralizedPareto:
                    return new GeneralizedPareto();
                case UnivariateDistributionType.Geometric:
                    return new Geometric();
                case UnivariateDistributionType.Gumbel:
                    return new Gumbel();
                case UnivariateDistributionType.InverseChiSquared:
                    return new InverseChiSquared();
                case UnivariateDistributionType.InverseGamma:
                    return new InverseGamma();
                case UnivariateDistributionType.KappaFour:
                    return new KappaFour();
                case UnivariateDistributionType.KernelDensity:
                    return new KernelDensity();
                case UnivariateDistributionType.LnNormal:
                    return new LnNormal();
                case UnivariateDistributionType.Logistic:
                    return new Logistic();
                case UnivariateDistributionType.LogNormal:
                    return new LogNormal();
                case UnivariateDistributionType.LogPearsonTypeIII:
                    return new LogPearsonTypeIII();
                case UnivariateDistributionType.NoncentralT:
                    return new NoncentralT();
                case UnivariateDistributionType.Normal:
                    return new Normal();
                case UnivariateDistributionType.Pareto:
                    return new Pareto();
                case UnivariateDistributionType.PearsonTypeIII:
                    return new PearsonTypeIII();
                case UnivariateDistributionType.Pert:
                    return new Pert();
                case UnivariateDistributionType.PertPercentile:
                    return new PertPercentile();
                case UnivariateDistributionType.PertPercentileZ:
                    return new PertPercentileZ();
                case UnivariateDistributionType.Poisson:
                    return new Poisson();
                case UnivariateDistributionType.Rayleigh:
                    return new Rayleigh();
                case UnivariateDistributionType.StudentT:
                    return new StudentT();
                case UnivariateDistributionType.Triangular:
                    return new Triangular();
                case UnivariateDistributionType.TruncatedNormal:
                    return new TruncatedNormal();
                case UnivariateDistributionType.Uniform:
                    return new Uniform();
                case UnivariateDistributionType.UniformDiscrete:
                    return new UniformDiscrete();
                case UnivariateDistributionType.VonMises:
                    return new VonMises();
                case UnivariateDistributionType.Weibull:
                    return new Weibull();
                case UnivariateDistributionType.CompetingRisks:
                case UnivariateDistributionType.Mixture:
                case UnivariateDistributionType.UserDefined:
                    throw new NotSupportedException(
                        "Distribution type " + distributionType + " requires external components and cannot be created without parameters.");
                default:
                    throw new ArgumentOutOfRangeException(
                        nameof(distributionType), distributionType, "The distribution type is not defined.");
            }
        }

        /// <summary>
        /// Attempts to create a distribution that can be initialized without external component distributions.
        /// </summary>
        /// <param name="distributionType">Distribution type.</param>
        /// <param name="distribution">
        /// When this method returns <see langword="true"/>, contains a distribution whose
        /// <see cref="UnivariateDistributionBase.Type"/> matches <paramref name="distributionType"/>;
        /// otherwise, <see langword="null"/>.
        /// </param>
        /// <returns>
        /// <see langword="true"/> when the requested distribution can be created without external components;
        /// otherwise, <see langword="false"/>.
        /// </returns>
        /// <remarks>
        /// Composite distributions such as <see cref="CompetingRisks"/> and <see cref="Mixture"/>,
        /// user-defined distributions, and undefined enumeration values cannot be created by this method.
        /// </remarks>
        public static bool TryCreateDistribution(UnivariateDistributionType distributionType, out UnivariateDistributionBase? distribution)
        {
            if (!Enum.IsDefined(typeof(UnivariateDistributionType), distributionType) ||
                distributionType == UnivariateDistributionType.CompetingRisks ||
                distributionType == UnivariateDistributionType.Mixture ||
                distributionType == UnivariateDistributionType.UserDefined)
            {
                distribution = null;
                return false;
            }

            distribution = CreateDistribution(distributionType);
            return true;
        }

        /// <summary>
        /// Creates a distribution from its serialized representation.
        /// </summary>
        /// <param name="xElement">The element to deserialize.</param>
        /// <returns>A validated univariate distribution.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the type or parameter data is missing or malformed.</exception>
        /// <exception cref="NotSupportedException">Thrown when the serialized type requires a user-provided implementation.</exception>
        public static UnivariateDistributionBase CreateDistribution(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));

            var typeAttribute = xElement.Attribute(nameof(UnivariateDistributionBase.Type));
            if (typeAttribute == null
                || !Enum.TryParse(typeAttribute.Value, out UnivariateDistributionType type)
                || !Enum.IsDefined(typeof(UnivariateDistributionType), type))
                throw new ArgumentException("The serialized distribution type is missing or invalid.", nameof(xElement));

            if (type == UnivariateDistributionType.Mixture)
                return Mixture.FromXElement(xElement)
                    ?? throw new ArgumentException("The serialized mixture is invalid.", nameof(xElement));
            if (type == UnivariateDistributionType.CompetingRisks)
                return CompetingRisks.FromXElement(xElement)
                    ?? throw new ArgumentException("The serialized competing-risks distribution is invalid.", nameof(xElement));
            if (type == UnivariateDistributionType.PertPercentile)
                return PertPercentile.FromXElement(xElement)
                    ?? throw new ArgumentException("The serialized percentile PERT distribution is invalid.", nameof(xElement));
            if (type == UnivariateDistributionType.PertPercentileZ)
                return PertPercentileZ.FromXElement(xElement)
                    ?? throw new ArgumentException("The serialized transformed percentile PERT distribution is invalid.", nameof(xElement));
            if (type == UnivariateDistributionType.Empirical)
                return EmpiricalDistribution.FromXElement(xElement);
            if (type == UnivariateDistributionType.KernelDensity)
                return KernelDensity.FromXElement(xElement);

            var distribution = CreateDistribution(type);
            var names = distribution.GetParameterPropertyNames;
            var values = new double[distribution.NumberOfParameters];
            for (int i = 0; i < values.Length; i++)
            {
                var parameterAttribute = xElement.Attribute(names[i]);
                if (parameterAttribute == null
                    || !double.TryParse(parameterAttribute.Value, System.Globalization.NumberStyles.Any, System.Globalization.CultureInfo.InvariantCulture, out values[i])
                    || !Tools.IsFinite(values[i]))
                    throw new ArgumentException("The serialized distribution parameter '" + names[i] + "' is missing or invalid.", nameof(xElement));
            }

            distribution.ValidateParameters(values, true);
            distribution.SetParameters(values);
            if (!distribution.ParametersValid)
                throw new ArgumentException("The serialized parameters do not define a valid distribution.", nameof(xElement));
            return distribution;
        }

    }
}