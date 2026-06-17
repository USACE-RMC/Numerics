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
    /// <para>
    ///     The transformation parameter lambda is restricted to the range used by
    ///     <see cref="YeoJohnson.FitLambda(System.Collections.Generic.IList{double})"/>, -5 &lt;= lambda &lt;= 5.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class YeoJohnsonLink : ILinkFunction
    {
        private const double MinimumLambda = -5d;
        private const double MaximumLambda = 5d;

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
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="lambda"/> is not finite or is outside [-5, 5].</exception>
        public YeoJohnsonLink(double lambda)
        {
            Lambda = ValidateLambda(lambda, nameof(lambda));
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="YeoJohnsonLink"/> class by estimating lambda from representative values.
        /// </summary>
        /// <param name="values">Representative values used to estimate lambda.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when fewer than two finite values are supplied or lambda fitting fails.</exception>
        public YeoJohnsonLink(double[] values)
        {
            if (values == null)
                throw new ArgumentNullException(nameof(values));
            if (values.Length < 2)
                throw new ArgumentException("At least 2 values are required to fit lambda.", nameof(values));
            for (int i = 0; i < values.Length; i++)
                if (!Tools.IsFinite(values[i]))
                    throw new ArgumentException("Every representative value must be finite.", nameof(values));

            Lambda = ValidateLambda(YeoJohnson.FitLambda(values), nameof(values));
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

            Lambda = ValidateLambda(double.Parse(lambdaAttribute.Value, NumberStyles.Any, CultureInfo.InvariantCulture), "Lambda");
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

        /// <summary>
        /// Validates a Yeo-Johnson lambda value.
        /// </summary>
        /// <param name="lambda">The lambda value to validate.</param>
        /// <param name="paramName">The parameter name to use in the exception.</param>
        /// <returns>The validated lambda value.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="lambda"/> is not finite or is outside [-5, 5].</exception>
        private static double ValidateLambda(double lambda, string paramName)
        {
            if (!Tools.IsFinite(lambda) || lambda < MinimumLambda || lambda > MaximumLambda)
                throw new ArgumentOutOfRangeException(paramName, "Lambda must be finite and in the range [-5, 5].");

            return lambda;
        }
    }
}
