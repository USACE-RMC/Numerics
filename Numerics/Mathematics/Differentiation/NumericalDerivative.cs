﻿using Numerics.Mathematics.LinearAlgebra;
using System;

namespace Numerics.Mathematics
{

    /// <summary>
    /// Contains methods for numerical differentiation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Numerical_differentiation" />
    /// <see href = "https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant" />
    /// <see href = "https://en.wikipedia.org/wiki/Hessian_matrix" />
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
        public static double Derivative(Func<double, double> f, double point, double stepSize = -1)
        {
            double h = stepSize <= 0 ? CalculateStepSize(point) : stepSize;
            //double h = stepSize;
            return (f(point + h) - f(point - h)) / (2d * h);
        }

        /// <summary>
        /// Computes the derivative of a function using Ridders' method of polynomial extrapolation.
        /// </summary>
        /// <param name="f">Function for which the derivative is to be evaluated.</param>
        /// <param name="point">The location where the derivative is to be evaluated.</param>
        /// <param name="stepSize">The initial estimate of the finite difference step size.</param>
        /// <param name="err">An estimate of the error in the derivative is returned by reference.</param>
        /// <remarks>
        /// Taken from "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017.
        /// </remarks>
        public static double RiddersMethod(Func<double, double> f, double point, double stepSize = -1, OptionalOut<double> err = null)
        {

            if (err == null) { err = new OptionalOut<double>() { Result = 0 }; }
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
            err.Result = big;
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
                    if (errt <= err.Result) // If error is decreased, save the improved answer.
                    {
                        err.Result = errt;
                        ans = a[j, i];
                    }
                }

                if (Math.Abs(a[i, i] - a[i - 1, i - 1]) >= safe * err.Result)
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
            double h = stepSize <= 0 ? CalculateStepSize(Tools.Min(point), 2) : stepSize;
            var hess= new double[point.Length, point.Length];
            var hi = new double[point.Length];
            var lo = new double[point.Length];
            point.CopyTo(hi, 0);
            point.CopyTo(lo, 0);
            for (int i = 0; i < point.Length; i++)
            {
                //double h = stepSize <= 0 ? CalculateStepSize(point[i], 2) : stepSize;
                hi[i] += h;
                lo[i] -= h;
                var grad1 = Gradient(f, hi, h);
                var grad2 = Gradient(f, lo, h);
                for (int j = 0; j < point.Length; j++)
                    hess[i, j] = (grad1[j] - grad2[j]) / (2d * h);
                hi[i] = point[i];
                lo[i] = point[i];
            }
            return hess;
        }

        /// <summary>
        /// A base step size value, h, will be scaled according to the function input parameter.
        /// </summary>
        /// <param name="x">The input parameter.</param>
        /// <param name="order">The order of the derivative.</param>
        public static double CalculateStepSize(double x, int order = 1)
        {
            return x != 0 ? Math.Pow(Tools.DoubleMachineEpsilon, 1d / (1d + order)) * Math.Abs(x) : Math.Pow(Tools.DoubleMachineEpsilon, 1d / (1d + order));
        }
    }
}