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
