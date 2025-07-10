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

namespace Numerics.Mathematics
{

    /// <summary>
    /// Contains methods for numerical differentiation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <b> Description: </b>
    /// These methods estimate the derivative of a function using values and other knowledge about the function.
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Numerical_differentiation" />
    /// </description></item>
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant" />
    /// </description></item>
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Hessian_matrix" />
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public sealed class NumericalDerivative
    {
        /// <summary>
        /// Computes the derivative of a function.
        /// </summary>
        /// <param name="f">Function for which the derivative is to be evaluated.</param>
        /// <param name="point">The location where the derivative is to be evaluated.</param>
        /// <param name="stepSize">Optional. The finite difference step size.</param>
        /// <remarks>The most common three point method is an average of a forward and backward difference derivative.</remarks>
        /// <returns>
        /// The derivative of the function f, evaluated at the given point
        /// </returns>
        public static double Derivative(Func<double, double> f, double point, double stepSize = -1)
        {
            double h = stepSize <= 0 ? CalculateStepSize(point) : stepSize;
            return (f(point + h) - f(point - h)) / (2d * h);
        }

        /// <summary>
        /// Computes the derivative of a function using Ridders' method of polynomial extrapolation.
        /// </summary>
        /// <param name="f">Function for which the derivative is to be evaluated.</param>
        /// <param name="point">The location where the derivative is to be evaluated.</param>
        /// <param name="err">Output. An estimate of the error in the derivative.</param>
        /// <param name="stepSize">Optional. The finite difference step size.</param>
        /// <remarks>
        /// References: Taken from "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017.
        /// </remarks>
        /// <returns>
        /// The derivative of the function f, evaluated at the given point
        /// </returns>
        public static double RiddersMethod(Func<double, double> f, double point, out double err, double stepSize = -1)
        {
            int ntab = 10;        // Sets maximum size of tableau.
            double con = 1.4d;         // Step size decreased by CON as each iteration
            double con2 = Math.Pow(con, 2d);
            double big = double.MaxValue;
            double safe = 2.0d;        // Return when error SAFE is worse than the best so far
            double errt;
            double fac;
            double hh = stepSize <= 0 ? CalculateStepSize(point) : stepSize;
            var ans = default(double);
            var a = new double[ntab + 1, ntab + 1];
            a[0, 0] = (f(point + hh) - f(point - hh)) / (2d * hh);
            err = big;
            for (int i = 1; i < ntab; i++)
            {
                // Successive columns in the Neville tableau will go to smaller step sizes and higher order of extrapolation.
                hh /= con;
                a[0, i] = (f(point + hh) - f(point - hh)) / (2d * hh); // Try new, smaller step size.
                fac = con2;
                // Compute extrapolation of various orders, requiring no new function evaluations.
                for (int j = 1; j <= i; j++)
                {
                    a[j, i] = (a[j - 1, i] * fac - a[j - 1, i - 1]) / (fac - 1.0d);
                    fac = con2 * fac;
                    errt = Math.Max(Math.Abs(a[j, i] - a[j - 1, i]), Math.Abs(a[j, i] - a[j - 1, i - 1]));
                    // The error strategy is to compare each new extrapolation to one order lower, both
                    // at the present step size and the previous one. 
                    if (errt <= err) // If error is decreased, save the improved answer.
                    {
                        err = errt;
                        ans = a[j, i];
                    }
                }

                if (Math.Abs(a[i, i] - a[i - 1, i - 1]) >= safe * err)
                    break;
                // If higher order is worse by a significant factor SAFE, then quit early. 
            }

            return ans;
        }

        /// <summary>
        /// Computes the gradient of a function using the symmetric difference quotient method.
        /// </summary>
        /// <param name="f">Function for which the derivative is to be evaluated.</param>
        /// <param name="point">The location where the derivative is to be evaluated.</param>
        /// <param name="stepSize">Optional. The finite difference step size.</param>
        /// <remarks>The most common three point method is an average of a forward and backward difference derivative.</remarks>
        /// <returns>
        /// The gradient of the function f, evaluated at the given points
        /// </returns>
        public static double[] Gradient(Func<double[], double> f, double[] point, double stepSize = -1)
        {
            double h;
            var grad = new double[point.Length];
            var hi = new double[point.Length];
            var lo = new double[point.Length];

            point.CopyTo(hi, 0);
            point.CopyTo(lo, 0);
            for (int i = 0; i < point.Length; i++)
            {
                h = stepSize <= 0 ? CalculateStepSize(point[i]) : stepSize;
                hi[i] += h;
                lo[i] -= h;
                grad[i] = (f(hi) - f(lo)) / (2d * h);
                hi[i] = point[i];
                lo[i] = point[i];
            }
            return grad;
        }

        /// <summary>
        /// Computes the Jacobian matrix at given function locations for each parameter point.
        /// </summary>
        /// <param name="f">Function for which the derivative is to be evaluate.</param>
        /// <param name="x">Functional value locations. Determines the rows of the Jacobian.</param>
        /// <param name="point">The location where the derivative is to be evaluated. Determines the columns of the Jacobian.</param>
        /// <param name="stepSize">Optional. The finite difference step size.</param>
        /// <returns>The Jacobian matrix.</returns>
        public static double[,] Jacobian(Func<double, double[], double> f, double[] x, double[] point, double stepSize = -1)
        {
            double h;
            var jac = new double[x.Length, point.Length];
            var hi = new double[point.Length];
            var lo = new double[point.Length];
            point.CopyTo(hi, 0);
            point.CopyTo(lo, 0);
            for (int i = 0; i < x.Length; i++)
            {
                for (int j = 0; j < point.Length; j++)
                {
                    h = stepSize <= 0 ? CalculateStepSize(point[i]) : stepSize;
                    hi[j] += h;
                    lo[j] -= h;
                    jac[i,j] = (f(x[i], hi) - f(x[i], lo)) / (2d * h);
                    hi[j] = point[j];
                    lo[j] = point[j];
                }
            }
            return jac;
        }

        /// <summary>
        /// Computes the Hessian matrix at a given point.
        /// </summary>
        /// <param name="f">Function which the derivative is to be evaluated.</param>
        /// <param name="point">The location where the derivative is to be evaluated.</param>
        /// <param name="stepSize">Optional. The finite difference step size..</param>
        /// <returns>The Hessian matrix.</returns>
        public static double[,] Hessian(Func<double[], double> f, double[] point, double stepSize = -1)
        {

            var hess = new double[point.Length, point.Length];
            var x = new double[point.Length];
            point.CopyTo(x, 0);
            double f1, f2, f3, f4;
            double hi, hj;

            for (int i = 0; i < point.Length; i++)
            {
                hi = stepSize <= 0 ? CalculateStepSize(point[i], 3) : stepSize;

                for (int j = 0; j < point.Length; j++)
                {
                    hj = stepSize <= 0 ? CalculateStepSize(point[j], 3) : stepSize;

                    if (i == j)
                    {
                        x[i] += hi;
                        f1 = f(x);
                        f2 = f(point);
                        x[i] -= 2 * hi;
                        f3 = f(x);
                        hess[i, j] = (f1 - 2 * f2 + f3) / (hi * hi);
                    }
                    else
                    {
                        x[i] += hi;
                        x[j] += hj;
                        f1 = f(x);
                        x[j] -= 2 * hj;
                        f2 = f(x);
                        x[i] -= 2 * hi;
                        x[j] += 2 * hj;
                        f3 = f(x);
                        x[j] -= 2 * hj;
                        f4 = f(x);
                        hess[i, j] = (f1 - f2 - f3 + f4) / (4 * hi * hj);
                    }
                    x[i] = point[i];
                    x[j] = point[j];
                }

            }

            return hess;
        }

        /// <summary>
        /// A base step size value, h, will be scaled according to the function input parameter.
        /// </summary>
        /// <param name="x">The input parameter.</param>
        /// <param name="order">The order of the derivative.</param>
        /// <returns> 
        /// Scaled h
        /// </returns>
        public static double CalculateStepSize(double x, int order = 1)
        {
            return Math.Pow(Tools.DoubleMachineEpsilon, 1d / (1d + order)) * (1.0 + Math.Abs(x));
        }
    }
}