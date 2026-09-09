using System;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{
    /// <summary>One-sided endpoint expansions for independent products of zero tails and infinite densities.</summary>
    internal static class DistributionEndpointTail
    {
        /// <summary>Returns tail ~ exp(logCoefficient)*distance^power*log(1/distance)^logPower.</summary>
        /// <param name="distribution">The distribution whose finite endpoint behavior is requested.</param>
        /// <param name="lower"><see langword="true"/> to describe the lower CDF tail; <see langword="false"/> to describe the upper survival tail.</param>
        /// <param name="power">The exponent applied to distance from the finite endpoint.</param>
        /// <param name="logPower">The exponent applied to the logarithm of the reciprocal endpoint distance.</param>
        /// <param name="logCoefficient">The logarithm of the expansion coefficient.</param>
        /// <returns><see langword="true"/> when the distribution has a recognized expansion; otherwise, <see langword="false"/>.</returns>
        /// <remarks>These are finite lower CDF or upper survival endpoint limits. Infinite power denotes
        /// faster-than-polynomial decay. No numerical endpoint offset or density floor is used.</remarks>
        internal static bool TryExpansion(UnivariateDistributionBase distribution, bool lower,
            out double power, out double logPower, out double logCoefficient)
        {
            power = logPower = logCoefficient = 0;
            switch (distribution)
            {
                case GammaDistribution gamma when lower:
                    power = gamma.Kappa;
                    logCoefficient = -Gamma.LogGamma(power + 1) - power * Math.Log(gamma.Theta);
                    return true;
                case Weibull weibull when lower:
                    power = weibull.Kappa;
                    logCoefficient = -power * Math.Log(weibull.Lambda);
                    return true;
                case Exponential exponential when lower:
                    power = 1; logCoefficient = -Math.Log(exponential.Alpha);
                    return true;
                case Uniform uniform:
                    power = 1; logCoefficient = -Math.Log(uniform.Max - uniform.Min);
                    return true;
                case GeneralizedPareto pareto:
                    power = lower ? 1 : 1 / pareto.Kappa;
                    logCoefficient = lower ? -Math.Log(pareto.Alpha) : power * (Math.Log(pareto.Kappa) - Math.Log(pareto.Alpha));
                    return true;
                case GeneralizedExtremeValue extreme:
                    power = lower ? double.PositiveInfinity : 1 / extreme.Kappa;
                    logCoefficient = lower ? 0 : power * (Math.Log(extreme.Kappa) - Math.Log(extreme.Alpha));
                    return true;
                case GeneralizedLogistic logistic:
                    power = 1 / Math.Abs(logistic.Kappa);
                    logCoefficient = power * (Math.Log(Math.Abs(logistic.Kappa)) - Math.Log(logistic.Alpha));
                    return true;
                case GeneralizedNormal _:
                case LnNormal _:
                case LogNormal _:
                    power = double.PositiveInfinity;
                    return true;
                case PearsonTypeIII pearson:
                    power = pearson.Alpha;
                    logCoefficient = -Gamma.LogGamma(power + 1) - power * Math.Log(Math.Abs(pearson.Beta));
                    return true;
                case LogPearsonTypeIII pearson:
                    if (pearson.Gamma == 0) { power = double.PositiveInfinity; return true; }
                    if (lower && pearson.Gamma < 0)
                    {
                        // A reflected gamma survival becomes an algebraic-logarithmic lower tail after exponentiation.
                        power = 1 / (Math.Abs(pearson.Beta) * Math.Log(pearson.Base));
                        logPower = pearson.Alpha - 1;
                        logCoefficient = -power * pearson.Xi * Math.Log(pearson.Base)
                            + logPower * Math.Log(power) - Gamma.LogGamma(pearson.Alpha);
                    }
                    else
                    {
                        power = pearson.Alpha;
                        double logEndpoint = pearson.Xi * Math.Log(pearson.Base);
                        logCoefficient = -Gamma.LogGamma(power + 1)
                            - power * (Math.Log(Math.Abs(pearson.Beta)) + Math.Log(Math.Log(pearson.Base)) + logEndpoint);
                    }
                    return true;
                case KappaFour kappa:
                    if (!lower)
                    {
                        power = 1 / kappa.Kappa;
                        logCoefficient = power * (Math.Log(kappa.Kappa) - Math.Log(kappa.Alpha));
                    }
                    else if (kappa.Hondo > 0)
                    {
                        power = 1 / kappa.Hondo;
                        logCoefficient = power * (kappa.Kappa * Math.Log(kappa.Hondo) - Math.Log(kappa.Alpha));
                    }
                    else if (kappa.Hondo < 0)
                    {
                        power = 1 / (kappa.Kappa * kappa.Hondo);
                        logCoefficient = Math.Log(-kappa.Hondo) / kappa.Hondo
                            + power * (Math.Log(-kappa.Kappa) - Math.Log(kappa.Alpha));
                    }
                    else power = double.PositiveInfinity;
                    return true;
                default:
                    return false;
            }
        }
    }
}
