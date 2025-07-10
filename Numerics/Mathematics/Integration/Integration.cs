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

namespace Numerics.Mathematics.Integration
{

    /// <summary>
    /// Contains methods for numerical integration.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <b> Description: </b>
    /// The basic problem of numerical integration is to approximate a solution to f(x) integrated over the interval [a,b].
    /// <para>
    /// <b> References: </b>
    /// <see href = "https://en.wikipedia.org/wiki/Numerical_integration" />
    /// </para>
    /// </remarks>
    public sealed class Integration
    {
        /// <summary>
        /// Returns the integral of a function between a and b by ten-point Gauss-Legendre integration. 
        /// </summary>
        /// <param name="f">The function to integrate.</param>
        /// <param name="a">Start point for integration.</param>
        /// <param name="b">End point for integration.</param>
        /// <returns>The value of a definite integral.</returns>
        public static double GaussLegendre(Func<double, double> f, double a, double b)
        {
            var x = new double[] { 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717 };
            var w = new double[] { 0.2955242247147529, 0.2692667193099963, 0.2190863625159821, 0.1494513491505806, 0.0666713443086881 };
            double xm = 0.5 * (b + a);
            double xr = 0.5 * (b - a);
            double s = 0;
            for (int j = 0; j < 5; j++)
            {
                double dx = xr * x[j];
                s += w[j] * (f(xm + dx) + f(xm - dx));
            }
            return s *= xr;
        }

        /// <summary>
        /// Numerical integration using the Trapezoidal Rule.
        /// </summary>
        /// <param name="f">The function to integrate.</param>
        /// <param name="a">Start point for integration.</param>
        /// <param name="b">End point for integration.</param>
        /// <param name="steps">Number of integration steps. Default = 2.</param>
        /// <returns>The value of a definite integral.</returns>
        public static double TrapezoidalRule(Func<double, double> f, double a, double b, int steps = 2)
        {
            double h = (b - a) / steps;
            double x = a;
            double sum = 0.5d * (f(a) + f(b));
            for (int i = 1; i <= steps - 1; i++)
            {
                x += h;
                sum += f(x);
            }
            return h * sum;
        }


        /// <summary>
        /// Numerical integration using Simpson's Rule.
        /// </summary>
        /// <param name="f">The function to integrate.</param>
        /// <param name="a">Start point for integration.</param>
        /// <param name="b">End point for integration.</param>
        /// <param name="steps">Number of integration steps. Default = 2.</param>
        /// <returns>The value of a definite integral.</returns>
        public static double SimpsonsRule(Func<double, double> f, double a, double b, int steps = 2)
        {
            double h = (b - a) / steps;
            double sum1 = f(a  + h / 2d);
            double sum2 = 0d;
            for (int i = 1; i <= steps - 1; i++)
            {
                sum1 += f(a + h * i + h / 2d);
                sum2 += f(a + h * i);
            }          
            return h / 6d * (f(a) + f(b) + 4d * sum1 + 2d * sum2);
        }

        /// <summary>
        /// Numerical integration using the Midpoint method.
        /// </summary>
        /// <param name="f">The function to integrate.</param>
        /// <param name="a">Start point for integration.</param>
        /// <param name="b">End point for integration.</param>
        /// <param name="steps">Number of integration steps. Default = 2.</param>
        /// <returns>The value of a definite integral.</returns>
        public static double Midpoint(Func<double, double> f, double a, double b, int steps = 2)
        {
            double h = (b - a) / steps;
            double x = a + h / 2d;
            double sum = 0;
            for (int i = 1; i <= steps; i++)
            {
                sum += f(x);
                x += h;
            }
            return h * sum;
        }
  
    }
}