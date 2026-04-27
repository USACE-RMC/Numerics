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
        public static UnivariateDistributionBase CreateDistribution(UnivariateDistributionType distributionType)
        {
            if (distributionType == UnivariateDistributionType.Bernoulli)
            {
                return new Bernoulli();
            }
            else if (distributionType == UnivariateDistributionType.Beta)
            {
                return new BetaDistribution();
            }
            else if (distributionType == UnivariateDistributionType.Binomial)
            {
                return new Binomial();
            }
            else if (distributionType == UnivariateDistributionType.Cauchy)
            {
                return new Cauchy();
            }
            else if (distributionType == UnivariateDistributionType.ChiSquared)
            {
                return new ChiSquared();
            }
            else if (distributionType == UnivariateDistributionType.Exponential)
            {
                return new Exponential();
            }
            else if (distributionType == UnivariateDistributionType.GammaDistribution)
            {
                return new GammaDistribution();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedBeta)
            {
                return new GeneralizedBeta();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedExtremeValue)
            {
                return new GeneralizedExtremeValue();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedLogistic)
            {
                return new GeneralizedLogistic();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedNormal)
            {
                return new GeneralizedNormal();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedPareto)
            {
                return new GeneralizedPareto();
            }
            else if (distributionType == UnivariateDistributionType.Geometric)
            {
                return new Geometric();
            }
            else if (distributionType == UnivariateDistributionType.Gumbel)
            {
                return new Gumbel();
            }
            else if (distributionType == UnivariateDistributionType.InverseChiSquared)
            {
                return new InverseChiSquared();
            }
            else if (distributionType == UnivariateDistributionType.InverseGamma)
            {
                return new InverseGamma();
            }
            else if (distributionType == UnivariateDistributionType.KappaFour)
            {
                return new KappaFour();
            }
            else if (distributionType == UnivariateDistributionType.LnNormal)
            {
                return new LnNormal();
            }
            else if (distributionType == UnivariateDistributionType.Logistic)
            {
                return new Logistic();
            }
            else if (distributionType == UnivariateDistributionType.LogNormal)
            {
                return new LogNormal();
            }
            else if (distributionType == UnivariateDistributionType.LogPearsonTypeIII)
            {
                return new LogPearsonTypeIII();
            }
            else if (distributionType == UnivariateDistributionType.NoncentralT)
            {
                return new NoncentralT();
            }
            else if (distributionType == UnivariateDistributionType.Normal)
            {
                return new Normal();
            }
            else if (distributionType == UnivariateDistributionType.Pareto)
            {
                return new Pareto();
            }
            else if (distributionType == UnivariateDistributionType.PearsonTypeIII)
            {
                return new PearsonTypeIII();
            }
            else if (distributionType == UnivariateDistributionType.Pert)
            {
                return new Pert();
            }
            else if (distributionType == UnivariateDistributionType.PertPercentile)
            {
                return new PertPercentile();
            }
            else if (distributionType == UnivariateDistributionType.PertPercentileZ)
            {
                return new PertPercentileZ();
            }
            else if (distributionType == UnivariateDistributionType.Poisson)
            {
                return new Poisson();
            }
            else if (distributionType == UnivariateDistributionType.Rayleigh)
            {
                return new Rayleigh();
            }
            else if (distributionType == UnivariateDistributionType.StudentT)
            {
                return new StudentT();
            }
            else if (distributionType == UnivariateDistributionType.Triangular)
            {
                return new Triangular();
            }
            else if (distributionType == UnivariateDistributionType.TruncatedNormal)
            {
                return new TruncatedNormal();
            }
            else if (distributionType == UnivariateDistributionType.Uniform)
            {
                return new Uniform();
            }
            else if (distributionType == UnivariateDistributionType.UniformDiscrete)
            {
                return new UniformDiscrete();
            }
            else if (distributionType == UnivariateDistributionType.Weibull)
            {
                return new Weibull();
            }

            // Default to Deterministic for unrecognized types
            return new Deterministic();
        }

        /// <summary>
        /// Create a distribution from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize into a univariate distribution.</param>
        /// <returns>
        /// A univariate distribution.
        /// </returns>
        public static UnivariateDistributionBase CreateDistribution(XElement xElement)
        {
            UnivariateDistributionType type = UnivariateDistributionType.Deterministic;
            var typeAttr = xElement.Attribute(nameof(UnivariateDistributionBase.Type));
            if (typeAttr != null)
            {
                Enum.TryParse(typeAttr.Value, out type);

                if (type == UnivariateDistributionType.Mixture)
                {
                    return Mixture.FromXElement(xElement)!;
                }
                else if (type == UnivariateDistributionType.CompetingRisks)
                {
                    return CompetingRisks.FromXElement(xElement)!;
                }
                else if (type == UnivariateDistributionType.PertPercentile)
                {
                    return PertPercentile.FromXElement(xElement)!;
                }
                else if (type == UnivariateDistributionType.PertPercentileZ)
                {
                    return PertPercentileZ.FromXElement(xElement)!;
                }
            }

            var dist = CreateDistribution(type);
            var names = dist.GetParameterPropertyNames;
            var parms = dist.GetParameters;
            var vals = new double[dist.NumberOfParameters];
            for (int i = 0; i < dist.NumberOfParameters; i++)
            {
                var paramAttr = xElement.Attribute(names[i]);
                if (paramAttr != null)
                {
                    double.TryParse(paramAttr.Value, System.Globalization.NumberStyles.Any, System.Globalization.CultureInfo.InvariantCulture, out vals[i]);
                }
            }
            dist.SetParameters(vals);
            return dist;
        }


    }
}