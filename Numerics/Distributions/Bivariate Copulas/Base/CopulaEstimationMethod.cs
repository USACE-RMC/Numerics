namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// Copula Estimation Methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public enum CopulaEstimationMethod
    {
        /// <summary>
        /// The full maximum likelihood estimation method requires parametric marginal distributions to predict the data. 
        /// The marginal distributions and the copula dependency are estimated simultaneously. 
        /// </summary>
        FullLikelihood,

        /// <summary>
        /// A semi-parametric approach that uses the plotting positions of the data to 
        /// estimate the copula dependency using maximum likelihood.
        /// </summary>
        PseudoLikelihood,

        /// <summary>
        /// The inference from margins (IFM) method includes two procedures:
        /// 1) marginal distributions are independently estimated from the observed values;
        /// 2) the copula dependency is estimate through the maximization of the likelihood function
        /// given the marginal distributions.
        /// </summary>
        InferenceFromMargins,

    }
}