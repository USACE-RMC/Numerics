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

namespace Numerics.Functions
{
    /// <summary>
    /// Enumeration of standard link function types for generalized linear models.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Each value corresponds to a canonical link for a specific GLM family.
    ///     Use <see cref="LinkFunctionFactory.Create(LinkFunctionType)"/> to obtain
    ///     an <see cref="ILinkFunction"/> instance from an enum value.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public enum LinkFunctionType
    {
        /// <summary>
        /// Identity link: &#951; = x. Canonical link for the Normal (Gaussian) family.
        /// </summary>
        Identity,

        /// <summary>
        /// Log link: &#951; = log(x). Canonical link for the Poisson and Exponential families.
        /// </summary>
        Log,

        /// <summary>
        /// Logit link: &#951; = log(x / (1 &#8722; x)). Canonical link for the Binomial family.
        /// </summary>
        Logit,

        /// <summary>
        /// Probit link: &#951; = &#934;&#8315;&#185;(x). Alternative link for the Binomial family using the standard normal quantile function.
        /// </summary>
        Probit,

        /// <summary>
        /// Complementary log-log link: &#951; = log(&#8722;log(1 &#8722; x)). Used for asymmetric binary response models.
        /// </summary>
        ComplementaryLogLog
    }
}
