namespace Numerics.Functions
{
    /// <summary>
    /// Enumeration of standard link function types for generalized linear models.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Each value corresponds to a canonical link for a specific GLM family.
    ///     Use <see cref="LinkFunctionFactory.Create(LinkFunctionType)"/> to obtain
    ///     an <see cref="ILinkFunction"/> instance from an enum value.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public enum LinkFunctionType
    {
        /// <summary>
        /// Identity link: &#951; = x. Canonical link for the Normal (Gaussian) family.
        /// </summary>
        Identity,

        /// <summary>
        /// Log link: &#951; = log(x). Canonical link for the Poisson and Exponential families.
        /// </summary>
        Log,

        /// <summary>
        /// Logit link: &#951; = log(x / (1 &#8722; x)). Canonical link for the Binomial family.
        /// </summary>
        Logit,

        /// <summary>
        /// Probit link: &#951; = &#934;&#8315;&#185;(x). Alternative link for the Binomial family using the standard normal quantile function.
        /// </summary>
        Probit,

        /// <summary>
        /// Complementary log-log link: &#951; = log(&#8722;log(1 &#8722; x)). Used for asymmetric binary response models.
        /// </summary>
        ComplementaryLogLog
    }
}
