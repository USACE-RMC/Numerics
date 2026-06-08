using System;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// Fisher z link function mapping correlations from (-1, 1) to the unconstrained real line.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Domain: -1 &lt; x &lt; 1. h(x) = atanh(x), h&#8315;&#185;(&#951;) = tanh(&#951;),
    ///     h&#8242;(x) = 1 / (1 - x^2).
    /// </para>
    /// <para>
    ///     This link is commonly used for correlation parameters and other signed bounded parameters.
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class FisherZLink : ILinkFunction
    {
        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is not in (-1, 1).</exception>
        public double Link(double x)
        {
            if (x <= -1d || x >= 1d)
                throw new ArgumentOutOfRangeException(nameof(x), "FisherZLink requires -1 < x < 1.");

            return 0.5d * Math.Log((1d + x) / (1d - x));
        }

        /// <inheritdoc/>
        public double InverseLink(double eta)
        {
            if (eta >= 0d)
            {
                double e = Math.Exp(-2d * eta);
                return (1d - e) / (1d + e);
            }

            double expEta = Math.Exp(2d * eta);
            return (expEta - 1d) / (expEta + 1d);
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="x"/> is not in (-1, 1).</exception>
        public double DLink(double x)
        {
            if (x <= -1d || x >= 1d)
                throw new ArgumentOutOfRangeException(nameof(x), "FisherZLink derivative requires -1 < x < 1.");

            return 1d / (1d - x * x);
        }

        /// <inheritdoc/>
        public XElement ToXElement() => new XElement(nameof(FisherZLink));
    }
}
