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
    /// Factory for creating <see cref="ILinkFunction"/> instances from <see cref="LinkFunctionType"/> enum values
    /// or from serialized <see cref="XElement"/> representations.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public static class LinkFunctionFactory
    {
        /// <summary>
        /// Creates an <see cref="ILinkFunction"/> instance corresponding to the specified link function type.
        /// </summary>
        /// <param name="type">The link function type.</param>
        /// <returns>A new <see cref="ILinkFunction"/> instance.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="type"/> is not a recognized link function type.</exception>
        public static ILinkFunction Create(LinkFunctionType type)
        {
            switch (type)
            {
                case LinkFunctionType.Identity:
                    return new IdentityLink();
                case LinkFunctionType.Log:
                    return new LogLink();
                case LinkFunctionType.Logit:
                    return new LogitLink();
                case LinkFunctionType.Probit:
                    return new ProbitLink();
                case LinkFunctionType.ComplementaryLogLog:
                    return new ComplementaryLogLogLink();
                default:
                    throw new ArgumentOutOfRangeException(nameof(type), $"Unknown link function type: {type}.");
            }
        }

        /// <summary>
        /// Creates an <see cref="ILinkFunction"/> instance from a serialized <see cref="XElement"/>.
        /// The element name determines the link function type.
        /// </summary>
        /// <param name="xElement">The XElement representing the link function. The element name must match a known link function class name.</param>
        /// <returns>A new <see cref="ILinkFunction"/> instance.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="NotSupportedException">Thrown when the element name does not correspond to a known link function type.</exception>
        public static ILinkFunction CreateFromXElement(XElement xElement)
        {
            if (xElement == null) throw new ArgumentNullException(nameof(xElement));
            switch (xElement.Name.LocalName)
            {
                case nameof(IdentityLink):
                    return new IdentityLink();
                case nameof(LogLink):
                    return new LogLink();
                case nameof(LogitLink):
                    return new LogitLink();
                case nameof(ProbitLink):
                    return new ProbitLink();
                case nameof(ComplementaryLogLogLink):
                    return new ComplementaryLogLogLink();
                default:
                    throw new NotSupportedException($"Unknown link function type: '{xElement.Name.LocalName}'.");
            }
        }
    }
}
