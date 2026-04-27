using System;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// Logit link function mapping the unit interval (0, 1) to the unconstrained real line.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Domain: x &#8712; (0, 1). h(x) = log(x / (1 &#8722; x)), h&#8315;&#185;(&#951;) = 1 / (1 + exp(&#8722;&#951;)), h&#8242;(x) = 1 / (x(1 &#8722; x)).
    /// </para>
    /// <para>
    ///     The logit link is the canonical link for the Binomial GLM family.
    ///     It is used for probability parameters, mixing proportions, and any parameter
    ///     naturally bounded to the unit interval.
    /// </para>
    /// <para>
    ///     The inverse link (sigmoid/logistic function) uses numerically stable formulation
    ///     to avoid overflow for large |&#951;|.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class LogitLink : ILinkFunction
    {
        /// <summary>
        /// Smallest admissible x to avoid log(0) blow-ups near the boundary.
        /// </summary>
        private const double MinX = 1e-12;

        /// <summary>
        /// Largest admissible x to avoid log(0) blow-ups near the boundary.
        /// </summary>
        private const double MaxX = 1.0 - 1e-12;

        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is not in (0, 1).</exception>
        public double Link(double x)
        {
            if (x <= 0.0 || x >= 1.0)
                throw new ArgumentOutOfRangeException(nameof(x), "LogitLink requires x in (0, 1).");
            x = Math.Max(MinX, Math.Min(MaxX, x));
            return Math.Log(x / (1.0 - x));
        }

        /// <inheritdoc/>
        public double InverseLink(double eta)
        {
            // Numerically stable sigmoid to avoid overflow for large |eta|.
            if (eta >= 0)
            {
                double e = Math.Exp(-eta);
                return 1.0 / (1.0 + e);
            }
            else
            {
                double e = Math.Exp(eta);
                return e / (1.0 + e);
            }
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is not in (0, 1).</exception>
        public double DLink(double x)
        {
            if (x <= 0.0 || x >= 1.0)
                throw new ArgumentOutOfRangeException(nameof(x), "LogitLink derivative requires x in (0, 1).");
            x = Math.Max(MinX, Math.Min(MaxX, x));
            return 1.0 / (x * (1.0 - x));
        }

        /// <inheritdoc/>
        public XElement ToXElement() => new XElement(nameof(LogitLink));
    }
}
