using System;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// Complementary log-log link function mapping the unit interval (0, 1) to the unconstrained real line.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Domain: x &#8712; (0, 1). h(x) = log(&#8722;log(1 &#8722; x)), h&#8315;&#185;(&#951;) = 1 &#8722; exp(&#8722;exp(&#951;)),
    ///     h&#8242;(x) = 1 / ((1 &#8722; x) &#183; (&#8722;log(1 &#8722; x))).
    /// </para>
    /// <para>
    ///     The complementary log-log link is used for asymmetric binary response models,
    ///     particularly when the probability of the event is small (rare events). It arises
    ///     naturally from extreme value (Gumbel) latent variable models.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class ComplementaryLogLogLink : ILinkFunction
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
                throw new ArgumentOutOfRangeException(nameof(x), "ComplementaryLogLogLink requires x in (0, 1).");
            x = Math.Max(MinX, Math.Min(MaxX, x));
            return Math.Log(-Math.Log(1.0 - x));
        }

        /// <inheritdoc/>
        public double InverseLink(double eta)
        {
            return 1.0 - Math.Exp(-Math.Exp(eta));
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is not in (0, 1).</exception>
        public double DLink(double x)
        {
            if (x <= 0.0 || x >= 1.0)
                throw new ArgumentOutOfRangeException(nameof(x), "ComplementaryLogLogLink derivative requires x in (0, 1).");
            x = Math.Max(MinX, Math.Min(MaxX, x));
            // h'(x) = 1 / ((1 - x) * (-log(1 - x)))
            double oneMinusX = 1.0 - x;
            double negLogOneMinusX = -Math.Log(oneMinusX);
            return 1.0 / Math.Max(oneMinusX * negLogOneMinusX, 1e-16);
        }

        /// <inheritdoc/>
        public XElement ToXElement() => new XElement(nameof(ComplementaryLogLogLink));
    }
}
