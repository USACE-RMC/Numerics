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
