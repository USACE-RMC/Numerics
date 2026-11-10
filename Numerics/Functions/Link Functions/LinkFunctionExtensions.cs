using System;
using System.Collections.Generic;

namespace Numerics.Functions
{
    /// <summary>
    /// List overloads of the <see cref="ILinkFunction"/> members.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Each overload calls the scalar member once per element, in order, so every result has the same bits as the
    /// scalar call. They are extension methods, rather than interface members, so that every implementation of
    /// <see cref="ILinkFunction"/>, including one outside this library, gets them without a change: a default
    /// interface implementation is not supported by .NET Framework 4.8.1, one of the target frameworks.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public static class LinkFunctionExtensions
    {
        /// <summary>
        /// Evaluates the link function at each real-space value.
        /// </summary>
        /// <param name="link">The link function.</param>
        /// <param name="x">The real-space values.</param>
        /// <returns>One link-space value per input, each equal to <see cref="ILinkFunction.Link(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="link"/> or <paramref name="x"/> is null.</exception>
        public static double[] Link(this ILinkFunction link, IList<double> x)
        {
            if (link == null) throw new ArgumentNullException(nameof(link));
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = link.Link(x[i]);
            return result;
        }

        /// <summary>
        /// Evaluates the inverse link function at each link-space value.
        /// </summary>
        /// <param name="link">The link function.</param>
        /// <param name="eta">The link-space values.</param>
        /// <returns>One real-space value per input, each equal to <see cref="ILinkFunction.InverseLink(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="link"/> or <paramref name="eta"/> is null.</exception>
        public static double[] InverseLink(this ILinkFunction link, IList<double> eta)
        {
            if (link == null) throw new ArgumentNullException(nameof(link));
            if (eta == null) throw new ArgumentNullException(nameof(eta));
            var result = new double[eta.Count];
            for (int i = 0; i < eta.Count; i++)
                result[i] = link.InverseLink(eta[i]);
            return result;
        }

        /// <summary>
        /// Evaluates the derivative of the link function at each real-space value.
        /// </summary>
        /// <param name="link">The link function.</param>
        /// <param name="x">The real-space values.</param>
        /// <returns>One derivative per input, each equal to <see cref="ILinkFunction.DLink(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="link"/> or <paramref name="x"/> is null.</exception>
        public static double[] DLink(this ILinkFunction link, IList<double> x)
        {
            if (link == null) throw new ArgumentNullException(nameof(link));
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = link.DLink(x[i]);
            return result;
        }
    }
}
