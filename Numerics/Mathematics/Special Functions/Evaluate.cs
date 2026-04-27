using System;

namespace Numerics.Mathematics.SpecialFunctions
{

    /// <summary>
    /// Evaluation functions useful for computing polynomials. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class Evaluate
    {

        /// <summary>
        /// Evaluates a double precision polynomial.
        /// </summary>
        /// <remarks>
        /// For sanity's sake, the value of N indicates the NUMBER of
        /// coefficients, or more precisely, the ORDER of the polynomial,
        /// rather than the DEGREE Of the polynomial. The two quantities
        /// differ by 1, but cause a great deal of confusion.
        /// <para>
        /// For example, a polynomial of order 4 in this function would be:
        /// coefficients[3]*x^3 + coefficients[2]*x^2 + coefficients[1]*x + coefficients[0]
        /// </para>
        /// <para>
        /// References: 
        /// <list type="bullet">
        /// <item><description>
        /// Based on a function contained in algorithm AS241, Applied Statistics, 1988, Vol. 37, No. 3.
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <param name="coefficients">The coefficients of the polynomial. Item[0] is the constant term. </param>
        /// <param name="x">The point at which the polynomial is to be evaluated.</param>
        /// <returns>
        /// The polynomial with the given coefficients evaluated at x
        /// </returns>
        public static double Polynomial(double[] coefficients, double x)
        {
            int n = coefficients.Length;
            double value = coefficients[n - 1];
            for (int i = n - 2; i >= 0; i -= 1)
            {
                value *= x;
                value += coefficients[i];
            }

            return value;
        }

        /// <summary>
        /// Evaluates a double precision polynomial. Coefficients are in reverse order.
        /// </summary>
        /// <remarks>
        /// For example, a polynomial of order 4 in this function would be:
        /// coefficients[0]*x^3 + coefficients[1]*x^2 + coefficients[2]*x + coefficients[3]
        /// </remarks>
        /// <param name="coefficients">The coefficients of the polynomial. The last element in the list is the constant term.</param>
        /// <param name="x">The point at which the polynomial is to be evaluated.</param>
        /// <param name="n"> Optional parameter to redefine the order of the polynomial to be n+1 </param>
        /// <returns>
        /// The polynomial with the given coefficients in reverse order evaluated at x
        /// </returns>
        public static double PolynomialRev(double[] coefficients, double x, int n = -1)
        {
            if (n > coefficients.Length)
            {
                throw new ArgumentOutOfRangeException("n cannot be greater than the number of coefficients");
            }
            else if (n == -1 || n == coefficients.Length)
            {
                n = coefficients.Length - 1;
            }

            double value = coefficients[0];
            for (int i = 1; i <= n; i++)
            {
                value *= x;
                value += coefficients[i];
            }

            return value;
        }

        /// <summary>
        /// Evaluates a double precision polynomial. Coefficients are in reverse order, and coefficient(N) = 1.0.
        /// </summary>
        /// <remarks>
        /// For example, a polynomial of order 4 in this function would be:
        /// x^3 + coefficients[0]*x^2 + coefficients[1]*x + coefficients[2]
        /// </remarks>
        /// <param name="coefficients"> The coefficients of the polynomial. The last element in the list is the constant term, and the first element 
        /// is the coefficient for the second term, the coefficient for the first term is always 1 </param>
        /// <param name="x"> The point at which the polynomial is to be evaluated </param>
        /// <returns>
        /// The polynomial with the given coefficients in reverse order, with 1 as the coefficient for the first term, evaluated at x
        /// </returns>
        public static double PolynomialRev_1(double[] coefficients, double x)
        {
            int n = coefficients.Length;
            double value = x + coefficients[0];
            for (int i = 1; i < n; i++)
            {
                value *= x;
                value += coefficients[i];
            }
            return value;
        }
    }
}