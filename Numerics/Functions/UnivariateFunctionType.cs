namespace Numerics.Functions
{
    /// <summary>
    /// The univariate function types.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// Members are append-only, like <see cref="Numerics.Distributions.UnivariateDistributionType"/>:
    /// downstream serialization and factory dispatch key on these values, so existing members are
    /// never renamed, reordered, or removed.
    /// </para>
    /// </remarks>
    public enum UnivariateFunctionType
    {
        /// <summary>
        /// The linear function Y = α + βX with optional additive Gaussian noise (<see cref="LinearFunction"/>).
        /// </summary>
        Linear,

        /// <summary>
        /// The power function Y = α(X − ξ)^β with optional log-space Gaussian noise and an optional
        /// inverse form (<see cref="PowerFunction"/>).
        /// </summary>
        Power,

        /// <summary>
        /// The tabular (nonparametric) function over uncertain ordered paired data (<see cref="TabularFunction"/>).
        /// </summary>
        Tabular,
    }
}
