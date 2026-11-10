using System;
using System.Collections.Generic;

namespace Numerics.Functions
{
    /// <summary>
    /// List overloads of the <see cref="IUnivariateFunction"/> members.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Each overload calls the scalar member once per element, in order, so every result has the same bits as the
    /// scalar call. They are extension methods, rather than interface members, so that every implementation of
    /// <see cref="IUnivariateFunction"/>, including one outside this library, gets them without a change: a default
    /// interface implementation is not supported by .NET Framework 4.8.1, one of the target frameworks.
    /// </para>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public static class UnivariateFunctionExtensions
    {
        /// <summary>
        /// Returns the function evaluated at each x-value.
        /// </summary>
        /// <param name="function">The function to evaluate.</param>
        /// <param name="x">The x-values.</param>
        /// <returns>One value per x-value, each equal to <see cref="IUnivariateFunction.Function(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="function"/> or <paramref name="x"/> is null.</exception>
        public static double[] Function(this IUnivariateFunction function, IList<double> x)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = function.Function(x[i]);
            return result;
        }

        /// <summary>
        /// Returns the inverse function evaluated at each y-value.
        /// </summary>
        /// <param name="function">The function to invert.</param>
        /// <param name="y">The y-values.</param>
        /// <returns>One value per y-value, each equal to <see cref="IUnivariateFunction.InverseFunction(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="function"/> or <paramref name="y"/> is null.</exception>
        public static double[] InverseFunction(this IUnivariateFunction function, IList<double> y)
        {
            if (function == null) throw new ArgumentNullException(nameof(function));
            if (y == null) throw new ArgumentNullException(nameof(y));
            var result = new double[y.Count];
            for (int i = 0; i < y.Count; i++)
                result[i] = function.InverseFunction(y[i]);
            return result;
        }
    }
}
