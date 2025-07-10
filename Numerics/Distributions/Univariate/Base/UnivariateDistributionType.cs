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

namespace Numerics.Distributions
{

    /// <summary>
    /// Univariate Distribution Types.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum UnivariateDistributionType
    {
        /// <summary>
        /// Chi-Squared distribution.
        /// </summary>
        ChiSquared,
        /// <summary>
        /// Bernoulli distribution.
        /// </summary>
        Bernoulli,
        /// <summary>
        /// Beta distribution.
        /// </summary>
        Beta,
        /// <summary>
        /// Binomial distribution.
        /// </summary>
        Binomial,
        /// <summary>
        /// Cauchy distribution.
        /// </summary>
        Cauchy,
        /// <summary>
        /// Competing risks distribution.
        /// </summary>
        CompetingRisks,
        /// <summary>
        /// Deterministic distribution.
        /// </summary>
        Deterministic,
        /// <summary>
        /// Empirical distribution.
        /// </summary>
        Empirical,
        /// <summary>
        /// Exponential distribution.
        /// </summary>
        Exponential,
        /// <summary>
        /// Gamma distribution.
        /// </summary>
        GammaDistribution,
        /// <summary>
        /// Generalized Beta distribution.
        /// </summary>
        GeneralizedBeta,
        /// <summary>
        /// Generalized Extreme Value (GEV) distribution.
        /// </summary>
        GeneralizedExtremeValue,
        /// <summary>
        /// Generalized Logistic (GLO) distribution.
        /// </summary>
        GeneralizedLogistic,
        /// <summary>
        /// Generalized Normal (GNO) distribution.
        /// </summary>
        GeneralizedNormal,
        /// <summary>
        /// Generalized Pareto (GPA) distribution.
        /// </summary>
        GeneralizedPareto,
        /// <summary>
        /// Geometric distribution.
        /// </summary>
        Geometric,
        /// <summary>
        /// Gumbel (EV1) distribution.
        /// </summary>
        Gumbel,
        /// <summary>
        /// Inverse Chi-Squared distribution.
        /// </summary>
        InverseChiSquared,
        /// <summary>
        /// Inverse Gamma distribution.
        /// </summary>
        InverseGamma,
        /// <summary>
        /// Kappa-4 distribution.
        /// </summary>
        KappaFour,
        /// <summary>
        /// Kernel Density distribution.
        /// </summary>
        KernelDensity,
        /// <summary>
        /// Log-Normal (base e) distribution.
        /// </summary>
        LnNormal,
        /// <summary>
        /// Logistic distribution.
        /// </summary>
        Logistic,
        /// <summary>
        /// Log-Normal distribution.
        /// </summary>
        LogNormal,
        /// <summary>
        /// Log-Pearson Type III (LPIII) distribution.
        /// </summary>
        LogPearsonTypeIII,
        /// <summary>
        /// Mixture distribution.
        /// </summary>
        Mixture,
        /// <summary>
        /// Non-central t distribution.
        /// </summary>
        NoncentralT,
        /// <summary>
        /// Normal distribution.
        /// </summary>
        Normal,
        /// <summary>
        /// Pareto distribution.
        /// </summary>
        Pareto,
        /// <summary>
        /// Pearson Type III (PIII) distribution.
        /// </summary>
        PearsonTypeIII,
        /// <summary>
        /// PERT distribution.
        /// </summary>
        Pert,
        /// <summary>
        /// PERT-Percentile distribution.
        /// </summary>
        PertPercentile,
        /// <summary>
        /// PERT-Percentile Z distribution.
        /// </summary>
        PertPercentileZ,
        /// <summary>
        /// Poisson distribution.
        /// </summary>
        Poisson,
        /// <summary>
        /// Rayleigh distribution.
        /// </summary>
        Rayleigh,
        /// <summary>
        /// Student t distribution.
        /// </summary>
        StudentT,
        /// <summary>
        /// Triangular distribution.
        /// </summary>
        Triangular,
        /// <summary>
        /// Truncated Normal distribution.
        /// </summary>
        TruncatedNormal,
        /// <summary>
        /// Uniform distribution.
        /// </summary>
        Uniform,
        /// <summary>
        /// Uniform-Discrete distribution.
        /// </summary>
        UniformDiscrete,
        /// <summary>
        /// User-defined distribution.
        /// </summary>
        UserDefined,
        /// <summary>
        /// Weibull distribution.
        /// </summary>
        Weibull
    }
}