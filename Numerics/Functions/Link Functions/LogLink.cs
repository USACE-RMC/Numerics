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
