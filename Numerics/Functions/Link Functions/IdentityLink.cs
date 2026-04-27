using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// Identity link function: h(x) = x. No transformation is applied.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     The identity link maps parameters directly without transformation.
    ///     It is the canonical link for the Normal (Gaussian) GLM family.
    /// </para>
    /// <para>
    ///     Domain: all real numbers. h(x) = x, h&#8315;&#185;(&#951;) = &#951;, h&#8242;(x) = 1.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class IdentityLink : ILinkFunction
    {
        /// <inheritdoc/>
        public double Link(double x) => x;

        /// <inheritdoc/>
        public double InverseLink(double eta) => eta;

        /// <inheritdoc/>
        public double DLink(double x) => 1.0;

        /// <inheritdoc/>
        public XElement ToXElement() => new XElement(nameof(IdentityLink));
    }
}
