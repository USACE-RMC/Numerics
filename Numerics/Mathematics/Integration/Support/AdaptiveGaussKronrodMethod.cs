namespace Numerics.Mathematics.Integration
{
    /// <summary>
    /// The adaptive strategies available to <see cref="AdaptiveGaussKronrod"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Both strategies use the 10-point Gauss, 21-point Kronrod rule. Local bisection is the cheaper choice for
    /// smooth, peaked and oscillatory integrands. QAGS is the better choice for endpoint singularities and for integrands that local bisection
    /// fails to converge on: on x^-0.5, log x, x^-0.9 and log(x)/sqrt(x) over [0, 1] it needs 231 to 315
    /// evaluations and is accurate to about 1E-13 or better, where local bisection needs 819 to 2247 and leaves
    /// errors of 1.6E-9 to 0.12.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W., and Kahaner, D. K. (1983). QUADPACK: A Subroutine Package for Automatic Integration. Springer-Verlag. <see href="https://doi.org/10.1007/978-3-642-61786-7"/></item>
    /// <item>Wynn, P. (1956). On a device for computing the e_m(S_n) transformation. Mathematical Tables and Other Aids to Computation, 10(54), 91-96. <see href="https://doi.org/10.2307/2002183"/></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public enum AdaptiveGaussKronrodMethod
    {
        /// <summary>
        /// Recursive local bisection: an interval is accepted when its Kronrod-Gauss difference meets the tolerance scaled
        /// by its share of the range. This is the default.
        /// </summary>
        LocalBisection,
        /// <summary>
        /// The QUADPACK QAGS algorithm: global bisection of the interval with the largest QUADPACK error estimate, with
        /// Wynn epsilon-algorithm extrapolation. It is suited to endpoint singularities and integrands the local method
        /// fails to converge on.
        /// </summary>
        Qags
    }
}
