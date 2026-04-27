using System;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// Log link function mapping positive reals to the unconstrained real line.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Domain: x &gt; 0. h(x) = log(x), h&#8315;&#185;(&#951;) = exp(&#951;), h&#8242;(x) = 1/x.
    /// </para>
    /// <para>
    ///     The log link is the canonical link for the Poisson and Exponential GLM families.
    ///     It is commonly used for scale parameters that must remain positive.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class LogLink : ILinkFunction
    {
        /// <summary>
        /// Smallest admissible x to avoid log(0) and 1/x blow-ups.
        /// </summary>
        private const double MinX = 1e-12;

        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is less than or equal to zero.</exception>
        public double Link(double x)
        {
            if (x <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(x), "LogLink requires x > 0.");
            if (x < MinX) x = MinX;
            return Math.Log(x);
        }

        /// <inheritdoc/>
        public double InverseLink(double eta)
        {
            return Math.Exp(eta);
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is less than or equal to zero.</exception>
        public double DLink(double x)
        {
            if (x <= 0.0)
                throw new ArgumentOutOfRangeException(nameof(x), "LogLink derivative requires x > 0.");
            if (x < MinX) x = MinX;
            return 1.0 / x;
        }

        /// <inheritdoc/>
        public XElement ToXElement() => new XElement(nameof(LogLink));
    }
}
