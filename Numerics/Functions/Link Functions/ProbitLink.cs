using System;
using System.Xml.Linq;
using Numerics.Distributions;

namespace Numerics.Functions
{
    /// <summary>
    /// Probit link function mapping the unit interval (0, 1) to the unconstrained real line
    /// using the standard normal quantile function.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Domain: x &#8712; (0, 1). h(x) = &#934;&#8315;&#185;(x), h&#8315;&#185;(&#951;) = &#934;(&#951;), h&#8242;(x) = 1 / &#966;(&#934;&#8315;&#185;(x)),
    ///     where &#934; is the standard normal CDF and &#966; is the standard normal PDF.
    /// </para>
    /// <para>
    ///     The probit link is an alternative link for the Binomial GLM family.
    ///     It assumes an underlying normally distributed latent variable.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class ProbitLink : ILinkFunction
    {
        /// <summary>
        /// Smallest admissible x to avoid quantile blow-ups near the boundary.
        /// </summary>
        private const double MinX = 1e-12;

        /// <summary>
        /// Largest admissible x to avoid quantile blow-ups near the boundary.
        /// </summary>
        private const double MaxX = 1.0 - 1e-12;

        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is not in (0, 1).</exception>
        public double Link(double x)
        {
            if (x <= 0.0 || x >= 1.0)
                throw new ArgumentOutOfRangeException(nameof(x), "ProbitLink requires x in (0, 1).");
            x = Math.Max(MinX, Math.Min(MaxX, x));
            return Normal.StandardZ(x);
        }

        /// <inheritdoc/>
        public double InverseLink(double eta)
        {
            return Normal.StandardCDF(eta);
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is not in (0, 1).</exception>
        public double DLink(double x)
        {
            if (x <= 0.0 || x >= 1.0)
                throw new ArgumentOutOfRangeException(nameof(x), "ProbitLink derivative requires x in (0, 1).");
            x = Math.Max(MinX, Math.Min(MaxX, x));
            double z = Normal.StandardZ(x);
            double pdf = Normal.StandardPDF(z);
            return 1.0 / Math.Max(pdf, 1e-16);
        }

        /// <inheritdoc/>
        public XElement ToXElement() => new XElement(nameof(ProbitLink));
    }
}
