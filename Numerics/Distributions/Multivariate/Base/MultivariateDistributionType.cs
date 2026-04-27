using System;

namespace Numerics.Distributions
{

    /// <summary>
    /// Multivariate Distribution Types.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum MultivariateDistributionType
    {
        /// <summary>
        /// Bivariate empirical distribution.
        /// </summary>
        BivariateEmpiricalDistribution,
        /// <summary>
        /// Multivariate Normal (MVN) distribution.
        /// </summary>
        MultivariateNormal,
        /// <summary>
        /// Dirichlet distribution.
        /// </summary>
        Dirichlet,
        /// <summary>
        /// Multinomial distribution.
        /// </summary>
        Multinomial,
        /// <summary>
        /// Multivariate Student's t-distribution.
        /// </summary>
        MultivariateStudentT
    }
}