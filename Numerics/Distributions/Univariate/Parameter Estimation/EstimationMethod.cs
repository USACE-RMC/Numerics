namespace Numerics.Distributions
{

    /// <summary>
    /// Probability Distribution Estimation Methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public enum ParameterEstimationMethod
    {
        /// <summary>
        /// Method of maximum likelihood
        /// </summary>
        MaximumLikelihood,
        /// <summary>
        /// Method of moments (or product moments)
        /// </summary>
        MethodOfMoments,
        /// <summary>
        /// Method of linear moments (or L-moments)
        /// </summary>
        MethodOfLinearMoments,
        /// <summary>
        /// Method of percentiles (or least squares)
        /// </summary>
        MethodOfPercentiles
    }
}