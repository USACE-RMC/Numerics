/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
