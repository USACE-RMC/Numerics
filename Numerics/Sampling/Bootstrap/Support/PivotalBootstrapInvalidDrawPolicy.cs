namespace Numerics.Sampling
{
    /// <summary>
    /// Specifies how invalid pivotal bootstrap draws are handled after the two-covariance transform.
    /// </summary>
    public enum PivotalBootstrapInvalidDrawPolicy
    {
        /// <summary>
        /// Drop invalid pivotal draws from the retained pivotal ensemble.
        /// </summary>
        Drop,

        /// <summary>
        /// Replace invalid pivotal draws with the corresponding accepted raw bootstrap fit.
        /// </summary>
        UseRaw,

        /// <summary>
        /// Replace invalid pivotal draws with the original parent fit.
        /// </summary>
        UseParent
    }
}
