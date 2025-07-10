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