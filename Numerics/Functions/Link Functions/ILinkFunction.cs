using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// Interface for link functions that transform parameters between real-space and link-space.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     A link function h provides a bijective, differentiable mapping from the parameter's
    ///     natural domain to an unconstrained link-space. This is used in generalized linear models (GLMs),
    ///     Bayesian MCMC estimation, and parametric bootstrap procedures.
    /// </para>
    /// <para>
    ///     Implementations must satisfy the round-trip identity: InverseLink(Link(x)) = x
    ///     for all x in the valid domain, and the derivative consistency: DLink(x) = dLink(x)/dx.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface ILinkFunction
    {
        /// <summary>
        /// Evaluates the link function mapping real-space to link-space: &#951; = h(x).
        /// </summary>
        /// <param name="x">The real-space value to transform.</param>
        /// <returns>The link-space value &#951;.</returns>
        double Link(double x);

        /// <summary>
        /// Evaluates the inverse link function mapping link-space back to real-space: x = h&#8315;&#185;(&#951;).
        /// </summary>
        /// <param name="eta">The link-space value to transform.</param>
        /// <returns>The real-space value x.</returns>
        double InverseLink(double eta);

        /// <summary>
        /// Evaluates the derivative of the link function with respect to x: h&#8242;(x) = d&#951;/dx.
        /// </summary>
        /// <param name="x">The real-space value at which to evaluate the derivative.</param>
        /// <returns>The derivative d&#951;/dx.</returns>
        double DLink(double x);

        /// <summary>
        /// Serializes the link function to an <see cref="XElement"/>.
        /// </summary>
        /// <returns>An XElement representing the link function and its parameters.</returns>
        XElement ToXElement();
    }
}
