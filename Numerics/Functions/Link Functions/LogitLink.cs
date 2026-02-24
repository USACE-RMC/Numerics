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
