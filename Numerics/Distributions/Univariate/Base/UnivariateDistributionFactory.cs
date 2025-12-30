/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
            var distribution = default(UnivariateDistributionBase);
            if (distributionType == UnivariateDistributionType.Bernoulli)
            {
                distribution = new Bernoulli();
            }
            else if (distributionType == UnivariateDistributionType.Beta)
            {
                distribution = new BetaDistribution();
            }
            else if (distributionType == UnivariateDistributionType.Binomial)
            {
                distribution = new Binomial();
            }
            else if (distributionType == UnivariateDistributionType.Cauchy)
            {
                distribution = new Cauchy();
            }
            else if (distributionType == UnivariateDistributionType.ChiSquared)
            {
                distribution = new ChiSquared();
            }
            else if (distributionType == UnivariateDistributionType.Deterministic)
            {
                distribution = new Deterministic();
            }
            else if (distributionType == UnivariateDistributionType.Exponential)
            {
                distribution = new Exponential();
            }
            else if (distributionType == UnivariateDistributionType.GammaDistribution)
            {
                distribution = new GammaDistribution();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedBeta)
            {
                distribution = new GeneralizedBeta();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedExtremeValue)
            {
                distribution = new GeneralizedExtremeValue();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedLogistic)
            {
                distribution = new GeneralizedLogistic();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedNormal)
            {
                distribution = new GeneralizedNormal();
            }
            else if (distributionType == UnivariateDistributionType.GeneralizedPareto)
            {
                distribution = new GeneralizedPareto();
            }
            else if (distributionType == UnivariateDistributionType.Geometric)
            {
                distribution = new Geometric();
            }
            else if (distributionType == UnivariateDistributionType.Gumbel)
            {
                distribution = new Gumbel();
            }
            else if (distributionType == UnivariateDistributionType.InverseChiSquared)
            {
                distribution = new InverseChiSquared();
            }
            else if (distributionType == UnivariateDistributionType.InverseGamma)
            {
                distribution = new InverseGamma();
            }
            else if (distributionType == UnivariateDistributionType.KappaFour)
            {
                distribution = new KappaFour();
            }
            else if (distributionType == UnivariateDistributionType.LnNormal)
            {
                distribution = new LnNormal();
            }
            else if (distributionType == UnivariateDistributionType.Logistic)
            {
                distribution = new Logistic();
            }
            else if (distributionType == UnivariateDistributionType.LogNormal)
            {
                distribution = new LogNormal();
            }
            else if (distributionType == UnivariateDistributionType.LogPearsonTypeIII)
            {
                distribution = new LogPearsonTypeIII();
            }
            else if (distributionType == UnivariateDistributionType.NoncentralT)
            {
                distribution = new NoncentralT();
            }
            else if (distributionType == UnivariateDistributionType.Normal)
            {
                distribution = new Normal();
            }
            else if (distributionType == UnivariateDistributionType.Pareto)
            {
                distribution = new Pareto();
            }
            else if (distributionType == UnivariateDistributionType.PearsonTypeIII)
            {
                distribution = new PearsonTypeIII();
            }
            else if (distributionType == UnivariateDistributionType.Pert)
            {
                distribution = new Pert();
            }
            else if (distributionType == UnivariateDistributionType.PertPercentile)
            {
                distribution = new PertPercentile();
            }
            else if (distributionType == UnivariateDistributionType.PertPercentileZ)
            {
                distribution = new PertPercentileZ();
            }
            else if (distributionType == UnivariateDistributionType.Poisson)
            {
                distribution = new Poisson();
            }
            else if (distributionType == UnivariateDistributionType.Rayleigh)
            {
                distribution = new Rayleigh();
            }
            else if (distributionType == UnivariateDistributionType.StudentT)
            {
                distribution = new StudentT();
            }
            else if (distributionType == UnivariateDistributionType.Triangular)
            {
                distribution = new Triangular();
            }
            else if (distributionType == UnivariateDistributionType.TruncatedNormal)
            {
                distribution = new TruncatedNormal();
            }
            else if (distributionType == UnivariateDistributionType.Uniform)
            {
                distribution = new Uniform();
            }
            else if (distributionType == UnivariateDistributionType.UniformDiscrete)
            {
                distribution = new UniformDiscrete();
            }
            else if (distributionType == UnivariateDistributionType.Weibull)
            {
                distribution = new Weibull();
            }
            if (distribution is null)
            {
                throw new ArgumentException("Distribution is not found.");
            }
            return distribution;
        }

        /// <summary>
        /// Create a distribution from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize into a univariate distribution.</param>
        /// <returns>
        /// A univariate distribution.
        /// </returns>
        public static UnivariateDistributionBase? CreateDistribution(XElement xElement)
        {
            UnivariateDistributionType type;
            UnivariateDistributionBase? dist = null;
            var xAttr = xElement.Attribute(nameof(UnivariateDistributionBase.Type));
            if (xAttr != null)
            {
                Enum.TryParse(xAttr.Value, out type);

                if (type == UnivariateDistributionType.Mixture)
                {
                    dist = Mixture.FromXElement(xElement);
                    return dist;
                }
                else if (type == UnivariateDistributionType.CompetingRisks)
                {
                    dist = CompetingRisks.FromXElement(xElement);
                    return dist;
                }
                else if (type == UnivariateDistributionType.PertPercentile)
                {
                    dist = PertPercentile.FromXElement(xElement);
                    return dist;
                }
                else if (type == UnivariateDistributionType.PertPercentileZ)
                {
                    dist = PertPercentileZ.FromXElement(xElement);
                    return dist;
                }
                else
                {
                    dist = CreateDistribution(type);
                }

            }
            if (dist is null)
            {
                throw new ArgumentException("Distribution is not found.");
            }
            var names = dist.GetParameterPropertyNames;
            var parms = dist.GetParameters;
            var vals = new double[dist.NumberOfParameters];
            for (int i = 0; i < dist.NumberOfParameters; i++)
            {
                var xAttrParm = xElement.Attribute(names[i]);
                if (xAttrParm != null)
                {
                    double.TryParse(xAttrParm.Value, System.Globalization.NumberStyles.Any, System.Globalization.CultureInfo.InvariantCulture, out vals[i]);
                }
            }
            dist.SetParameters(vals);
            return dist;
        }


    }
}