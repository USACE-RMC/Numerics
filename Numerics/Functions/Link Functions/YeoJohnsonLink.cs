using Numerics.Data.Statistics;
using System;
using System.Globalization;
using System.Xml.Linq;

namespace Numerics.Functions
{
    /// <summary>
    /// Yeo-Johnson power-transformation link function.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     The Yeo-Johnson transformation is defined on the full real line. It is useful
    ///     for skew or shape parameters that do not have a closed-form
    ///     variance-stabilizing link.
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class YeoJohnsonLink : ILinkFunction
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="YeoJohnsonLink"/> class with lambda equal to 1.
        /// </summary>
        public YeoJohnsonLink()
        {
            Lambda = 1d;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="YeoJohnsonLink"/> class.
        /// </summary>
        /// <param name="lambda">The Yeo-Johnson transformation parameter.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="lambda"/> is not finite.</exception>
        public YeoJohnsonLink(double lambda)
        {
            if (double.IsNaN(lambda) || double.IsInfinity(lambda))
                throw new ArgumentOutOfRangeException(nameof(lambda), "Lambda must be finite.");

            Lambda = lambda;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="YeoJohnsonLink"/> class by estimating lambda from representative values.
        /// </summary>
        /// <param name="values">Representative values used to estimate lambda.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        public YeoJohnsonLink(double[] values)
        {
            if (values == null)
                throw new ArgumentNullException(nameof(values));
            if (values.Length < 2)
                throw new ArgumentException("At least 2 values are required to fit lambda.", nameof(values));

            Lambda = YeoJohnson.FitLambda(values);
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="YeoJohnsonLink"/> class from its XML representation.
        /// </summary>
        /// <param name="xElement">The XML element to read.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xElement"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the element does not contain a Lambda attribute.</exception>
        public YeoJohnsonLink(XElement xElement)
        {
            if (xElement == null)
                throw new ArgumentNullException(nameof(xElement));

            XAttribute? lambdaAttribute = xElement.Attribute("Lambda");
            if (lambdaAttribute == null)
                throw new ArgumentException("The YeoJohnsonLink element must contain a Lambda attribute.", nameof(xElement));

            Lambda = double.Parse(lambdaAttribute.Value, NumberStyles.Any, CultureInfo.InvariantCulture);
        }

        /// <summary>
        /// Gets the Yeo-Johnson transformation parameter.
        /// </summary>
        public double Lambda { get; }

        /// <inheritdoc/>
        public double Link(double x)
        {
            return YeoJohnson.Transform(x, Lambda);
        }

        /// <inheritdoc/>
        public double InverseLink(double eta)
        {
            return YeoJohnson.InverseTransform(eta, Lambda);
        }

        /// <inheritdoc/>
        public double DLink(double x)
        {
            return YeoJohnson.Derivative(x, Lambda);
        }

        /// <inheritdoc/>
        public XElement ToXElement() => new XElement(nameof(YeoJohnsonLink), new XAttribute("Lambda", Lambda.ToString("G17", CultureInfo.InvariantCulture)));
    }
}
