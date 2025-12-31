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

using Numerics.Mathematics.LinearAlgebra;
using System;

namespace Numerics.Mathematics
{
    /// <summary>
    /// Contains methods for numerical differentiation using finite difference approximations.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This class provides various finite difference methods for estimating derivatives of functions.
    /// Methods include forward, backward, and central differences for first and second derivatives,
    /// as well as advanced techniques like Ridders' method and adaptive step sizing.
    /// The class also supports multivariable calculus operations including Jacobian matrices,
    /// gradients, and Hessian matrices with boundary constraint handling.
    /// </para>
    /// <para>
    /// <b> Key Features: </b>
    /// </para>
    /// <list type="bullet">
    /// <item><description>Multiple finite difference schemes (forward, backward, central)</description></item>
    /// <item><description>Adaptive step sizing with boundary awareness</description></item>
    /// <item><description>Ridders' method for high-accuracy derivative estimation</description></item>
    /// <item><description>Support for constrained optimization (with bounds)</description></item>
    /// <item><description>Robust error handling and numerical stability checks</description></item>
    /// </list>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <list type="bullet">
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Numerical_differentiation"/>
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Finite_difference"/>
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant"/>
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Hessian_matrix"/>
    /// </description></item>
    /// <item><description>
    /// Numerical Recipes: The Art of Scientific Computing, Third Edition. Press et al. 2017.
    /// </description></item>
    /// </list>
    /// </remarks>
    public sealed class NumericalDerivative
    {
        #region First Derivatives - Single Variable

        /// <summary>
        /// Computes the first derivative using forward difference approximation.
        /// </summary>
        /// <param name="f">Function to differentiate.</param>
        /// <param name="point">Point at which to evaluate the derivative.</param>
        /// <param name="stepSize">Finite difference step size. If negative or zero, an optimal step is calculated.</param>
        /// <returns>Approximation of f'(point).</returns>
        /// <remarks>
        /// <para>
        /// Uses the formula: f'(x) ≈ [f(x+h) - f(x)] / h
        /// </para>
        /// <para>
        /// This is a first-order accurate method with truncation error O(h).
        /// It's useful when evaluations at x-h are problematic or when only forward
        /// information is available.
        /// </para>
        /// </remarks>
        public static double ForwardDifference(Func<double, double> f, double point, double stepSize = -1)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));

            double h = stepSize <= 0 ? CalculateStepSize(point) : stepSize;
            return (f(point + h) - f(point)) / h;
        }

        /// <summary>
        /// Computes the first derivative using backward difference approximation.
        /// </summary>
        /// <param name="f">Function to differentiate.</param>
        /// <param name="point">Point at which to evaluate the derivative.</param>
        /// <param name="stepSize">Finite difference step size. If negative or zero, an optimal step is calculated.</param>
        /// <returns>Approximation of f'(point).</returns>
        /// <remarks>
        /// <para>
        /// Uses the formula: f'(x) ≈ [f(x) - f(x-h)] / h
        /// </para>
        /// <para>
        /// This is a first-order accurate method with truncation error O(h).
        /// It's useful when evaluations at x+h are problematic or when only backward
        /// information is available.
        /// </para>
        /// </remarks>
        public static double BackwardDifference(Func<double, double> f, double point, double stepSize = -1)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));

            double h = stepSize <= 0 ? CalculateStepSize(point) : stepSize;
            return (f(point) - f(point - h)) / h;
        }

        /// <summary>
        /// Computes the first derivative using central difference approximation.
        /// </summary>
        /// <param name="f">Function to differentiate.</param>
        /// <param name="point">Point at which to evaluate the derivative.</param>
        /// <param name="stepSize">Finite difference step size. If negative or zero, an optimal step is calculated.</param>
        /// <returns>Approximation of f'(point).</returns>
        /// <remarks>
        /// <para>
        /// Uses the formula: f'(x) ≈ [f(x+h) - f(x-h)] / (2h)
        /// </para>
        /// <para>
        /// This is a second-order accurate method with truncation error O(h²), making it
        /// significantly more accurate than forward or backward differences for the same step size.
        /// This is the recommended method when function evaluations are available on both sides of the point.
        /// </para>
        /// </remarks>
        public static double CentralDifference(Func<double, double> f, double point, double stepSize = -1)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));

            double h = stepSize <= 0 ? CalculateStepSize(point) : stepSize;
            return (f(point + h) - f(point - h)) / (2.0 * h);
        }

        /// <summary>
        /// Computes the first derivative of a function (alias for CentralDifference).
        /// </summary>
        /// <param name="f">Function for which the derivative is to be evaluated.</param>
        /// <param name="point">The location where the derivative is to be evaluated.</param>
        /// <param name="stepSize">Optional. The finite difference step size.</param>
        /// <returns>The derivative of the function f, evaluated at the given point.</returns>
        /// <remarks>
        /// This method uses the central difference approximation, which is an average of 
        /// forward and backward differences and provides second-order accuracy O(h²).
        /// </remarks>
        public static double Derivative(Func<double, double> f, double point, double stepSize = -1)
        {
            return CentralDifference(f, point, stepSize);
        }

        /// <summary>
        /// Computes the derivative of a function using Ridders' method of polynomial extrapolation.
        /// </summary>
        /// <param name="f">Function for which the derivative is to be evaluated.</param>
        /// <param name="point">The location where the derivative is to be evaluated.</param>
        /// <param name="err">Output. An estimate of the error in the derivative.</param>
        /// <param name="stepSize">Optional. The initial finite difference step size.</param>
        /// <returns>The derivative of the function f, evaluated at the given point.</returns>
        /// <remarks>
        /// <para>
        /// Ridders' method uses polynomial extrapolation with successively smaller step sizes
        /// to achieve high accuracy. It automatically estimates the error and stops when
        /// further refinement no longer improves the result.
        /// </para>
        /// <para>
        /// This method is particularly useful when:
        /// </para>
        /// <list type="bullet">
        /// <item><description>High accuracy is required</description></item>
        /// <item><description>The function is smooth</description></item>
        /// <item><description>You need an error estimate</description></item>
        /// </list>
        /// <para>
        /// Reference: Numerical Recipes: The Art of Scientific Computing, Third Edition. Press et al. 2017.
        /// </para>
        /// </remarks>
        public static double RiddersMethod(Func<double, double> f, double point, out double err, double stepSize = -1)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));

            // Algorithm parameters
            const int ntab = 10;              // Maximum size of extrapolation tableau
            const double con = 1.4;           // Step size reduction factor at each iteration
            double con2 = con * con;          // Square of the reduction factor
            const double big = double.MaxValue;
            const double safe = 2.0;          // Safety factor for early termination

            double hh = stepSize <= 0 ? CalculateStepSize(point) : stepSize;
            var ans = 0.0;
            var a = new double[ntab, ntab];

            // Initialize with central difference at largest step size
            a[0, 0] = (f(point + hh) - f(point - hh)) / (2.0 * hh);
            err = big;

            // Build tableau of successively smaller step sizes and higher-order extrapolations
            for (int i = 1; i < ntab; i++)
            {
                // Reduce step size
                hh /= con;

                // Compute derivative at new, smaller step size
                a[0, i] = (f(point + hh) - f(point - hh)) / (2.0 * hh);

                double fac = con2;

                // Perform polynomial extrapolation using Neville's algorithm
                for (int j = 1; j <= i; j++)
                {
                    // Extrapolate to h → 0 using previous estimates
                    a[j, i] = (a[j - 1, i] * fac - a[j - 1, i - 1]) / (fac - 1.0);
                    fac *= con2;

                    // Estimate error by comparing neighboring extrapolations
                    double errt = Math.Max(Math.Abs(a[j, i] - a[j - 1, i]),
                                          Math.Abs(a[j, i] - a[j - 1, i - 1]));

                    // Keep the best estimate
                    if (errt <= err)
                    {
                        err = errt;
                        ans = a[j, i];
                    }
                }

                // Terminate if higher order extrapolation is significantly worse
                if (Math.Abs(a[i, i] - a[i - 1, i - 1]) >= safe * err)
                    break;
            }

            return ans;
        }

        #endregion

        #region Second Derivatives - Single Variable

        /// <summary>
        /// Computes the second derivative using central difference approximation.
        /// </summary>
        /// <param name="f">Function to differentiate twice.</param>
        /// <param name="point">Point at which to evaluate the second derivative.</param>
        /// <param name="stepSize">Finite difference step size. If negative or zero, an optimal step is calculated.</param>
        /// <returns>Approximation of f''(point).</returns>
        /// <remarks>
        /// <para>
        /// Uses the formula: f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h²
        /// </para>
        /// <para>
        /// This is a second-order accurate method with truncation error O(h²).
        /// For second derivatives, a larger step size (order 2) is used by default to balance
        /// truncation error and rounding error.
        /// </para>
        /// </remarks>
        public static double SecondDerivative(Func<double, double> f, double point, double stepSize = -1)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));

            double h = stepSize <= 0 ? CalculateStepSize(point, order: 2) : stepSize;
            return (f(point + h) - 2.0 * f(point) + f(point - h)) / (h * h);
        }

        /// <summary>
        /// Computes the second derivative using forward difference approximation.
        /// </summary>
        /// <param name="f">Function to differentiate twice.</param>
        /// <param name="point">Point at which to evaluate the second derivative.</param>
        /// <param name="stepSize">Finite difference step size. If negative or zero, an optimal step is calculated.</param>
        /// <returns>Approximation of f''(point).</returns>
        /// <remarks>
        /// <para>
        /// Uses the three-point forward formula: f''(x) ≈ [f(x) - 2f(x+h) + f(x+2h)] / h²
        /// </para>
        /// <para>
        /// This is a second-order accurate method useful when backward evaluations are not available.
        /// </para>
        /// </remarks>
        public static double SecondDerivativeForward(Func<double, double> f, double point, double stepSize = -1)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));

            double h = stepSize <= 0 ? CalculateStepSize(point, order: 2) : stepSize;
            return (f(point) - 2.0 * f(point + h) + f(point + 2.0 * h)) / (h * h);
        }

        /// <summary>
        /// Computes the second derivative using backward difference approximation.
        /// </summary>
        /// <param name="f">Function to differentiate twice.</param>
        /// <param name="point">Point at which to evaluate the second derivative.</param>
        /// <param name="stepSize">Finite difference step size. If negative or zero, an optimal step is calculated.</param>
        /// <returns>Approximation of f''(point).</returns>
        /// <remarks>
        /// <para>
        /// Uses the three-point backward formula: f''(x) ≈ [f(x-2h) - 2f(x-h) + f(x)] / h²
        /// </para>
        /// <para>
        /// This is a second-order accurate method useful when forward evaluations are not available.
        /// </para>
        /// </remarks>
        public static double SecondDerivativeBackward(Func<double, double> f, double point, double stepSize = -1)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));

            double h = stepSize <= 0 ? CalculateStepSize(point, order: 2) : stepSize;
            return (f(point - 2.0 * h) - 2.0 * f(point - h) + f(point)) / (h * h);
        }

        #endregion

        #region Jacobian Matrices

        /// <summary>
        /// Computes the Jacobian matrix at given function locations for each parameter point.
        /// </summary>
        /// <param name="f">Multivariate function f(x_i, θ) where x_i is a fixed value and θ is the parameter vector.</param>
        /// <param name="x">Function value locations (determines the rows of the Jacobian).</param>
        /// <param name="point">Parameter vector location (determines the columns of the Jacobian).</param>
        /// <param name="stepSize">Optional. The finite difference step size.</param>
        /// <returns>The Jacobian matrix J where J[i,j] = ∂f(x_i, θ)/∂θ_j evaluated at θ = point.</returns>
        /// <remarks>
        /// <para>
        /// The Jacobian matrix contains all first-order partial derivatives of a multivariate function.
        /// For a function f: R^n → R^m, the Jacobian is an m×n matrix where the (i,j)-th entry is ∂f_i/∂x_j.
        /// </para>
        /// <para>
        /// This overload is designed for functions where some variables (x) are held constant while
        /// computing derivatives with respect to parameters (point).
        /// </para>
        /// <para>
        /// Uses central differences for improved accuracy O(h²).
        /// </para>
        /// </remarks>
        public static double[,] Jacobian(Func<double, double[], double> f, double[] x, double[] point, double stepSize = -1)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (point == null) throw new ArgumentNullException(nameof(point));

            var jac = new double[x.Length, point.Length];
            var hi = new double[point.Length];
            var lo = new double[point.Length];

            // Initialize parameter vectors for forward and backward differences
            point.CopyTo(hi, 0);
            point.CopyTo(lo, 0);

            // Compute each element of the Jacobian
            for (int i = 0; i < x.Length; i++)
            {
                for (int j = 0; j < point.Length; j++)
                {
                    // BUG FIX: Use point[j] instead of point[i] for step size calculation
                    double h = stepSize <= 0 ? CalculateStepSize(point[j]) : stepSize;

                    // Perturb the j-th parameter
                    hi[j] += h;
                    lo[j] -= h;

                    // Central difference approximation
                    jac[i, j] = (f(x[i], hi) - f(x[i], lo)) / (2.0 * h);

                    // Reset for next iteration
                    hi[j] = point[j];
                    lo[j] = point[j];
                }
            }

            return jac;
        }

        /// <summary>
        /// Computes the Jacobian matrix of a vector-valued function g: R^p → R^q at theta, 
        /// with adaptive, bound-aware steps.
        /// </summary>
        /// <param name="g">Vector-valued function mapping p parameters to q outputs.</param>
        /// <param name="theta">Parameter vector at which to evaluate the Jacobian.</param>
        /// <param name="lowerBounds">Optional. Lower bounds for parameters.</param>
        /// <param name="upperBounds">Optional. Upper bounds for parameters.</param>
        /// <param name="relStep">Relative step size for finite differences (default 1e-5).</param>
        /// <param name="absStep">Absolute minimum step size (default 1e-7).</param>
        /// <param name="maxBacktrack">Maximum number of step size reductions if function evaluations fail (default 5).</param>
        /// <returns>The q×p Jacobian matrix J where J[i,j] = ∂g_i/∂θ_j.</returns>
        /// <remarks>
        /// <para>
        /// This method automatically adapts the finite difference scheme based on parameter bounds:
        /// </para>
        /// <list type="bullet">
        /// <item><description>Uses central differences when sufficient space is available on both sides</description></item>
        /// <item><description>Switches to forward differences near upper bounds</description></item>
        /// <item><description>Switches to backward differences near lower bounds</description></item>
        /// <item><description>Reduces step size if function evaluations fail (return NaN or Infinity)</description></item>
        /// </list>
        /// <para>
        /// This adaptive approach ensures robust derivative estimation even for constrained optimization problems.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown if g or theta is null.</exception>
        /// <exception cref="InvalidOperationException">Thrown if g(theta) returns null.</exception>
        public static double[,] Jacobian(
            Func<double[], double[]> g,
            double[] theta,
            double[] lowerBounds = null!,
            double[] upperBounds = null!,
            double relStep = 1e-5,
            double absStep = 1e-7,
            int maxBacktrack = 5)
        {
            if (g == null) throw new ArgumentNullException(nameof(g));
            if (theta == null) throw new ArgumentNullException(nameof(theta));

            int p = theta.Length;
            var g0 = g(theta);
            if (g0 == null) throw new InvalidOperationException("g(theta) returned null.");
            int q = g0.Length;

            var J = new double[q, p];
            var theta0 = (double[])theta.Clone();

            // Compute each column of the Jacobian (derivatives with respect to each parameter)
            for (int j = 0; j < p; j++)
            {
                var col = new double[q];
                bool success = false;

                // Adaptive step sizing
                double scale = Math.Abs(theta0[j]) + 1.0;
                double hjBase = Math.Max(absStep, relStep * scale);
                double hj = hjBase;

                // Try different step sizes and difference schemes until successful
                for (int backtrack = 0; backtrack <= maxBacktrack && !success; backtrack++)
                {
                    // Determine available space relative to bounds
                    double roomLeft = AvailableLeft(theta0, j, lowerBounds);
                    double roomRight = AvailableRight(theta0, j, upperBounds);

                    // Choose finite difference scheme based on available space
                    bool canCentral = roomLeft >= hj && roomRight >= hj;
                    bool useForward = !canCentral && roomRight >= hj;
                    bool useBackward = !canCentral && !useForward && roomLeft >= hj;

                    try
                    {
                        if (canCentral)
                        {
                            // Central difference: most accurate when possible
                            var thetaPlus = (double[])theta0.Clone();
                            var thetaMinus = (double[])theta0.Clone();
                            thetaPlus[j] += hj;
                            thetaMinus[j] -= hj;

                            ClampInPlace(thetaPlus, lowerBounds, upperBounds);
                            ClampInPlace(thetaMinus, lowerBounds, upperBounds);

                            double actualStep = thetaPlus[j] - thetaMinus[j];
                            if (Math.Abs(actualStep) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            var gPlus = g(thetaPlus);
                            var gMinus = g(thetaMinus);

                            if (IsBad(gPlus, q) || IsBad(gMinus, q))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            for (int i = 0; i < q; i++)
                            {
                                col[i] = (gPlus[i] - gMinus[i]) / actualStep;
                            }
                            success = true;
                        }
                        else if (useForward)
                        {
                            // Forward difference: when at or near lower bound
                            var thetaPlus = (double[])theta0.Clone();
                            thetaPlus[j] += hj;

                            ClampInPlace(thetaPlus, lowerBounds, upperBounds);

                            double actualStep = thetaPlus[j] - theta0[j];
                            if (Math.Abs(actualStep) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            var gPlus = g(thetaPlus);

                            if (IsBad(gPlus, q))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            for (int i = 0; i < q; i++)
                            {
                                col[i] = (gPlus[i] - g0[i]) / actualStep;
                            }
                            success = true;
                        }
                        else if (useBackward)
                        {
                            // Backward difference: when at or near upper bound
                            var thetaMinus = (double[])theta0.Clone();
                            thetaMinus[j] -= hj;

                            ClampInPlace(thetaMinus, lowerBounds, upperBounds);

                            double actualStep = theta0[j] - thetaMinus[j];
                            if (Math.Abs(actualStep) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            var gMinus = g(thetaMinus);

                            if (IsBad(gMinus, q))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            for (int i = 0; i < q; i++)
                            {
                                col[i] = (g0[i] - gMinus[i]) / actualStep;
                            }
                            success = true;
                        }
                        else
                        {
                            // No room available: reduce step and retry
                            hj *= 0.5;
                        }
                    }
                    catch
                    {
                        // Function evaluation failed: reduce step and retry
                        hj *= 0.5;
                    }
                }

                // If all attempts failed, set derivatives to zero
                if (!success)
                {
                    col.Fill(0.0);
                }

                // Store column in Jacobian matrix
                for (int i = 0; i < q; i++)
                {
                    J[i, j] = col[i];
                }
            }

            return J;
        }

        #endregion

        #region Gradient Vectors

        /// <summary>
        /// Computes the gradient of a scalar function f: R^p → R at theta, with adaptive, bound-aware steps.
        /// </summary>
        /// <param name="f">Scalar-valued function mapping p parameters to a single output.</param>
        /// <param name="theta">Parameter vector at which to evaluate the gradient.</param>
        /// <param name="lowerBounds">Optional. Lower bounds for parameters.</param>
        /// <param name="upperBounds">Optional. Upper bounds for parameters.</param>
        /// <param name="relStep">Relative step size for finite differences (default 1e-5).</param>
        /// <param name="absStep">Absolute minimum step size (default 1e-7).</param>
        /// <param name="maxBacktrack">Maximum number of step size reductions if function evaluations fail (default 5).</param>
        /// <returns>The gradient vector ∇f where grad[j] = ∂f/∂θ_j.</returns>
        /// <remarks>
        /// <para>
        /// The gradient vector points in the direction of steepest ascent of the function.
        /// For optimization, the negative gradient points toward the minimum.
        /// </para>
        /// <para>
        /// This method uses the same adaptive finite difference strategy as the Jacobian computation:
        /// central differences when possible, switching to one-sided differences near boundaries.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown if f or theta is null.</exception>
        /// <exception cref="ArithmeticException">Thrown if f(theta) is not finite.</exception>
        public static double[] Gradient(
            Func<double[], double> f,
            double[] theta,
            double[] lowerBounds = null!,
            double[] upperBounds = null!,
            double relStep = 1e-5,
            double absStep = 1e-7,
            int maxBacktrack = 5)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));
            if (theta == null) throw new ArgumentNullException(nameof(theta));

            int p = theta.Length;
            var grad = new double[p];
            var theta0 = (double[])theta.Clone();
            double f0 = f(theta0);

            if (!Tools.IsFinite(f0))
            {
                throw new ArithmeticException("f(theta) is not finite.");
            }

            // Compute each component of the gradient
            for (int j = 0; j < p; j++)
            {
                bool success = false;
                double scale = Math.Abs(theta0[j]) + 1.0;
                double hjBase = Math.Max(absStep, relStep * scale);
                double hj = hjBase;

                for (int backtrack = 0; backtrack <= maxBacktrack && !success; backtrack++)
                {
                    double roomLeft = AvailableLeft(theta0, j, lowerBounds);
                    double roomRight = AvailableRight(theta0, j, upperBounds);

                    bool canCentral = roomLeft >= hj && roomRight >= hj;
                    bool useForward = !canCentral && roomRight >= hj;
                    bool useBackward = !canCentral && !useForward && roomLeft >= hj;

                    try
                    {
                        if (canCentral)
                        {
                            var thetaPlus = (double[])theta0.Clone();
                            var thetaMinus = (double[])theta0.Clone();
                            thetaPlus[j] += hj;
                            thetaMinus[j] -= hj;

                            ClampInPlace(thetaPlus, lowerBounds, upperBounds);
                            ClampInPlace(thetaMinus, lowerBounds, upperBounds);

                            double actualStep = thetaPlus[j] - thetaMinus[j];
                            if (Math.Abs(actualStep) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            double fPlus = f(thetaPlus);
                            double fMinus = f(thetaMinus);

                            if (!Tools.IsFinite(fPlus) || !Tools.IsFinite(fMinus))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            grad[j] = (fPlus - fMinus) / actualStep;
                            success = true;
                        }
                        else if (useForward)
                        {
                            var thetaPlus = (double[])theta0.Clone();
                            thetaPlus[j] += hj;

                            ClampInPlace(thetaPlus, lowerBounds, upperBounds);

                            double actualStep = thetaPlus[j] - theta0[j];
                            if (Math.Abs(actualStep) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            double fPlus = f(thetaPlus);

                            if (!Tools.IsFinite(fPlus))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            grad[j] = (fPlus - f0) / actualStep;
                            success = true;
                        }
                        else if (useBackward)
                        {
                            var thetaMinus = (double[])theta0.Clone();
                            thetaMinus[j] -= hj;

                            ClampInPlace(thetaMinus, lowerBounds, upperBounds);

                            double actualStep = theta0[j] - thetaMinus[j];
                            if (Math.Abs(actualStep) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            double fMinus = f(thetaMinus);

                            if (!Tools.IsFinite(fMinus))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            grad[j] = (f0 - fMinus) / actualStep;
                            success = true;
                        }
                        else
                        {
                            hj *= 0.5;
                        }
                    }
                    catch
                    {
                        hj *= 0.5;
                    }
                }

                if (!success)
                {
                    grad[j] = 0.0;
                }
            }

            return grad;
        }

        #endregion

        #region Hessian Matrices

        /// <summary>
        /// Computes the Hessian matrix of a scalar function f: R^p → R at theta, 
        /// with adaptive, bound-aware steps.
        /// </summary>
        /// <param name="f">Scalar-valued function mapping p parameters to a single output.</param>
        /// <param name="theta">Parameter vector at which to evaluate the Hessian.</param>
        /// <param name="lowerBounds">Optional. Lower bounds for parameters.</param>
        /// <param name="upperBounds">Optional. Upper bounds for parameters.</param>
        /// <param name="relStep">Relative step size for finite differences (default 1e-4, larger for second derivatives).</param>
        /// <param name="absStep">Absolute minimum step size (default 1e-6).</param>
        /// <param name="maxBacktrack">Maximum number of step size reductions if function evaluations fail (default 6).</param>
        /// <returns>The p×p symmetric Hessian matrix H where H[i,j] = ∂²f/(∂θ_i ∂θ_j).</returns>
        /// <remarks>
        /// <para>
        /// The Hessian matrix contains all second-order partial derivatives of a scalar function.
        /// It characterizes the local curvature of the function and is used in:
        /// </para>
        /// <list type="bullet">
        /// <item><description>Newton's method for optimization</description></item>
        /// <item><description>Determining if a critical point is a minimum, maximum, or saddle point</description></item>
        /// <item><description>Quadratic approximations of functions</description></item>
        /// <item><description>Uncertainty quantification in parameter estimation</description></item>
        /// </list>
        /// <para>
        /// The Hessian is symmetric for smooth functions (by Schwarz's theorem), so H[i,j] = H[j,i].
        /// This implementation enforces symmetry by averaging the upper and lower triangular elements.
        /// </para>
        /// <para>
        /// Uses central differences when interior; otherwise switches to one-sided formulas near boundaries.
        /// Diagonal elements use three-point formulas; off-diagonals use four-point stencils (central) or
        /// three-point one-sided formulas when necessary.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown if f or theta is null.</exception>
        /// <exception cref="ArithmeticException">Thrown if f(theta) is not finite.</exception>
        public static double[,] Hessian(
            Func<double[], double> f,
            double[] theta,
            double[] lowerBounds = null!,
            double[] upperBounds = null!,
            double relStep = 1e-4,
            double absStep = 1e-6,
            int maxBacktrack = 6)
        {
            if (f == null) throw new ArgumentNullException(nameof(f));
            if (theta == null) throw new ArgumentNullException(nameof(theta));

            int p = theta.Length;
            var H = new double[p, p];
            var x0 = (double[])theta.Clone();
            double f0 = f(x0);

            if (!Tools.IsFinite(f0))
            {
                throw new ArithmeticException("f(theta) is not finite.");
            }

            // Compute diagonal elements (∂²f/∂θ_j²)
            for (int j = 0; j < p; j++)
            {
                bool success = false;
                double scale = Math.Abs(x0[j]) + 1.0;
                double hjBase = Math.Max(absStep, relStep * scale);
                double hj = hjBase;

                for (int backtrack = 0; backtrack <= maxBacktrack && !success; backtrack++)
                {
                    double roomLeft = AvailableLeft(x0, j, lowerBounds);
                    double roomRight = AvailableRight(x0, j, upperBounds);

                    bool canCentral = roomLeft >= hj && roomRight >= hj;
                    bool canForward = roomRight >= 2.0 * hj;
                    bool canBackward = roomLeft >= 2.0 * hj;

                    try
                    {
                        if (canCentral)
                        {
                            // Central difference: f''(x) ≈ [f(x+h) - 2f(x) + f(x-h)] / h²
                            var xPlus = (double[])x0.Clone();
                            var xMinus = (double[])x0.Clone();
                            xPlus[j] += hj;
                            xMinus[j] -= hj;

                            ClampInPlace(xPlus, lowerBounds, upperBounds);
                            ClampInPlace(xMinus, lowerBounds, upperBounds);

                            double fPlus = f(xPlus);
                            double fMinus = f(xMinus);

                            if (!Tools.IsFinite(fPlus) || !Tools.IsFinite(fMinus))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            double actualStep = xPlus[j] - xMinus[j];
                            if (Math.Abs(actualStep) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            double hEff = 0.5 * actualStep;
                            H[j, j] = (fPlus - 2.0 * f0 + fMinus) / (hEff * hEff);
                            success = true;
                        }
                        else if (canForward)
                        {
                            // Forward difference: f''(x) ≈ [f(x) - 2f(x+h) + f(x+2h)] / h²
                            var x1 = (double[])x0.Clone();
                            var x2 = (double[])x0.Clone();
                            x1[j] += hj;
                            x2[j] += 2.0 * hj;

                            ClampInPlace(x1, lowerBounds, upperBounds);
                            ClampInPlace(x2, lowerBounds, upperBounds);

                            double f1 = f(x1);
                            double f2 = f(x2);

                            if (!Tools.IsFinite(f1) || !Tools.IsFinite(f2))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            double hEff = x1[j] - x0[j];
                            if (Math.Abs(hEff) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            H[j, j] = (f0 - 2.0 * f1 + f2) / (hEff * hEff);
                            success = true;
                        }
                        else if (canBackward)
                        {
                            // Backward difference: f''(x) ≈ [f(x-2h) - 2f(x-h) + f(x)] / h²
                            var x1 = (double[])x0.Clone();
                            var x2 = (double[])x0.Clone();
                            x1[j] -= hj;
                            x2[j] -= 2.0 * hj;

                            ClampInPlace(x1, lowerBounds, upperBounds);
                            ClampInPlace(x2, lowerBounds, upperBounds);

                            double f1 = f(x1);
                            double f2 = f(x2);

                            if (!Tools.IsFinite(f1) || !Tools.IsFinite(f2))
                            {
                                hj *= 0.5;
                                continue;
                            }

                            double hEff = x0[j] - x1[j];
                            if (Math.Abs(hEff) < 1e-16 * scale)
                            {
                                hj *= 0.5;
                                continue;
                            }

                            H[j, j] = (f2 - 2.0 * f1 + f0) / (hEff * hEff);
                            success = true;
                        }
                        else
                        {
                            hj *= 0.5;
                        }
                    }
                    catch
                    {
                        hj *= 0.5;
                    }
                }

                if (!success)
                {
                    H[j, j] = 0.0;
                }
            }

            // Compute off-diagonal elements (mixed partials ∂²f/(∂θ_i ∂θ_j), i<j)
            for (int i = 0; i < p; i++)
            {
                for (int j = i + 1; j < p; j++)
                {
                    bool success = false;

                    double scaleI = Math.Abs(x0[i]) + 1.0;
                    double scaleJ = Math.Abs(x0[j]) + 1.0;
                    double hiBase = Math.Max(absStep, relStep * scaleI);
                    double hjBase = Math.Max(absStep, relStep * scaleJ);
                    double hi = hiBase;
                    double hj = hjBase;

                    for (int backtrack = 0; backtrack <= maxBacktrack && !success; backtrack++)
                    {
                        double leftI = AvailableLeft(x0, i, lowerBounds);
                        double rightI = AvailableRight(x0, i, upperBounds);
                        double leftJ = AvailableLeft(x0, j, lowerBounds);
                        double rightJ = AvailableRight(x0, j, upperBounds);

                        bool centralI = leftI >= hi && rightI >= hi;
                        bool centralJ = leftJ >= hj && rightJ >= hj;

                        try
                        {
                            if (centralI && centralJ)
                            {
                                // Central 4-point stencil for mixed partial:
                                // ∂²f/(∂i∂j) ≈ [f(x+hi,y+hj) - f(x+hi,y-hj) - f(x-hi,y+hj) + f(x-hi,y-hj)] / (4·hi·hj)
                                var xPlusPlus = (double[])x0.Clone();
                                var xPlusMinus = (double[])x0.Clone();
                                var xMinusPlus = (double[])x0.Clone();
                                var xMinusMinus = (double[])x0.Clone();

                                xPlusPlus[i] += hi; xPlusPlus[j] += hj;
                                xPlusMinus[i] += hi; xPlusMinus[j] -= hj;
                                xMinusPlus[i] -= hi; xMinusPlus[j] += hj;
                                xMinusMinus[i] -= hi; xMinusMinus[j] -= hj;

                                ClampInPlace(xPlusPlus, lowerBounds, upperBounds);
                                ClampInPlace(xPlusMinus, lowerBounds, upperBounds);
                                ClampInPlace(xMinusPlus, lowerBounds, upperBounds);
                                ClampInPlace(xMinusMinus, lowerBounds, upperBounds);

                                double fPlusPlus = f(xPlusPlus);
                                double fPlusMinus = f(xPlusMinus);
                                double fMinusPlus = f(xMinusPlus);
                                double fMinusMinus = f(xMinusMinus);

                                if (!Tools.IsFinite(fPlusPlus) || !Tools.IsFinite(fPlusMinus) ||
                                    !Tools.IsFinite(fMinusPlus) || !Tools.IsFinite(fMinusMinus))
                                {
                                    hi *= 0.5;
                                    hj *= 0.5;
                                    continue;
                                }

                                double actualStepI = xPlusPlus[i] - xMinusMinus[i];
                                double actualStepJ = xPlusPlus[j] - xMinusMinus[j];

                                if (Math.Abs(actualStepI) < 1e-16 * scaleI ||
                                    Math.Abs(actualStepJ) < 1e-16 * scaleJ)
                                {
                                    hi *= 0.5;
                                    hj *= 0.5;
                                    continue;
                                }

                                double mixedPartial = (fPlusPlus - fPlusMinus - fMinusPlus + fMinusMinus) /
                                                     (actualStepI * actualStepJ);

                                H[i, j] = mixedPartial;
                                H[j, i] = mixedPartial;
                                success = true;
                            }
                            else
                            {
                                // One-sided mixed partial formula
                                // Prefer forward-forward when possible
                                bool forwardForward = (rightI >= hi && rightJ >= hj);
                                bool backwardBackward = (leftI >= hi && leftJ >= hj);

                                if (forwardForward)
                                {
                                    // ∂²f/(∂i∂j) ≈ [f(x+hi,y+hj) - f(x+hi,y) - f(x,y+hj) + f(x,y)] / (hi·hj)
                                    var xI = (double[])x0.Clone();
                                    var xJ = (double[])x0.Clone();
                                    var xIJ = (double[])x0.Clone();

                                    xI[i] += hi;
                                    xJ[j] += hj;
                                    xIJ[i] += hi; xIJ[j] += hj;

                                    ClampInPlace(xI, lowerBounds, upperBounds);
                                    ClampInPlace(xJ, lowerBounds, upperBounds);
                                    ClampInPlace(xIJ, lowerBounds, upperBounds);

                                    double fI = f(xI);
                                    double fJ = f(xJ);
                                    double fIJ = f(xIJ);

                                    if (!Tools.IsFinite(fI) || !Tools.IsFinite(fJ) || !Tools.IsFinite(fIJ))
                                    {
                                        hi *= 0.5;
                                        hj *= 0.5;
                                        continue;
                                    }

                                    double actualStepI = xI[i] - x0[i];
                                    double actualStepJ = xJ[j] - x0[j];

                                    if (Math.Abs(actualStepI) < 1e-16 * scaleI ||
                                        Math.Abs(actualStepJ) < 1e-16 * scaleJ)
                                    {
                                        hi *= 0.5;
                                        hj *= 0.5;
                                        continue;
                                    }

                                    double mixedPartial = (fIJ - fI - fJ + f0) / (actualStepI * actualStepJ);

                                    H[i, j] = mixedPartial;
                                    H[j, i] = mixedPartial;
                                    success = true;
                                }
                                else if (backwardBackward)
                                {
                                    // Backward-backward version
                                    var xI = (double[])x0.Clone();
                                    var xJ = (double[])x0.Clone();
                                    var xIJ = (double[])x0.Clone();

                                    xI[i] -= hi;
                                    xJ[j] -= hj;
                                    xIJ[i] -= hi; xIJ[j] -= hj;

                                    ClampInPlace(xI, lowerBounds, upperBounds);
                                    ClampInPlace(xJ, lowerBounds, upperBounds);
                                    ClampInPlace(xIJ, lowerBounds, upperBounds);

                                    double fI = f(xI);
                                    double fJ = f(xJ);
                                    double fIJ = f(xIJ);

                                    if (!Tools.IsFinite(fI) || !Tools.IsFinite(fJ) || !Tools.IsFinite(fIJ))
                                    {
                                        hi *= 0.5;
                                        hj *= 0.5;
                                        continue;
                                    }

                                    double actualStepI = x0[i] - xI[i];
                                    double actualStepJ = x0[j] - xJ[j];

                                    if (Math.Abs(actualStepI) < 1e-16 * scaleI ||
                                        Math.Abs(actualStepJ) < 1e-16 * scaleJ)
                                    {
                                        hi *= 0.5;
                                        hj *= 0.5;
                                        continue;
                                    }

                                    double mixedPartial = (fIJ - fI - fJ + f0) / (actualStepI * actualStepJ);

                                    H[i, j] = mixedPartial;
                                    H[j, i] = mixedPartial;
                                    success = true;
                                }
                                else
                                {
                                    // Insufficient room: reduce step sizes
                                    hi *= 0.5;
                                    hj *= 0.5;
                                }
                            }
                        }
                        catch
                        {
                            hi *= 0.5;
                            hj *= 0.5;
                        }
                    }

                    if (!success)
                    {
                        H[i, j] = 0.0;
                        H[j, i] = 0.0;
                    }
                }
            }

            // Final symmetry enforcement: average upper and lower triangular parts
            // This reduces numerical noise from finite difference approximations
            for (int i = 0; i < p; i++)
            {
                for (int j = i + 1; j < p; j++)
                {
                    double average = 0.5 * (H[i, j] + H[j, i]);
                    H[i, j] = average;
                    H[j, i] = average;
                }
            }

            return H;
        }

        #endregion

        #region Step Size Calculation

        /// <summary>
        /// Calculates an optimal step size for numerical differentiation based on the input parameter value.
        /// </summary>
        /// <param name="x">The input parameter value.</param>
        /// <param name="order">The order of the derivative (1 for first derivative, 2 for second derivative, etc.).</param>
        /// <returns>An adaptive step size scaled to the magnitude of x and appropriate for the derivative order.</returns>
        /// <remarks>
        /// <para>
        /// The step size h is chosen to balance two competing sources of error:
        /// </para>
        /// <list type="number">
        /// <item><description>
        /// <b>Truncation error:</b> Decreases as h gets smaller. For an n-th order method, 
        /// truncation error is O(h^n).
        /// </description></item>
        /// <item><description>
        /// <b>Rounding error:</b> Increases as h gets smaller because we're subtracting 
        /// nearly equal numbers (catastrophic cancellation).
        /// </description></item>
        /// </list>
        /// <para>
        /// The optimal step size minimizes the total error. For first derivatives with central differences,
        /// this gives h ≈ ε^(1/3) where ε is machine epsilon (~2.22e-16 for double precision).
        /// For second derivatives, h ≈ ε^(1/4).
        /// </para>
        /// <para>
        /// The formula used is: h = ε^(1/(1+order)) · (1 + |x|)
        /// </para>
        /// <para>
        /// The (1 + |x|) term scales the step size relative to the magnitude of the input,
        /// ensuring that the relative step size is appropriate regardless of the scale of x.
        /// </para>
        /// </remarks>
        public static double CalculateStepSize(double x, int order = 1)
        {
            // Machine epsilon for double precision: approximately 2.22e-16
            // For first derivatives: ε^(1/2) ≈ 1.5e-8
            // For second derivatives: ε^(1/3) ≈ 6.1e-6
            return Math.Pow(Tools.DoubleMachineEpsilon, 1.0 / (1.0 + order)) * (1.0 + Math.Abs(x));
        }

        #endregion

        #region Helper Methods

        /// <summary>
        /// Clamps all elements of a vector to stay within specified bounds.
        /// </summary>
        /// <param name="v">Vector to clamp (modified in place).</param>
        /// <param name="lowerBounds">Optional lower bounds for each element.</param>
        /// <param name="upperBounds">Optional upper bounds for each element.</param>
        private static void ClampInPlace(double[] v, double[]? lowerBounds, double[]? upperBounds)
        {
            if (lowerBounds == null && upperBounds == null)
                return;

            for (int k = 0; k < v.Length; k++)
            {
                if (lowerBounds != null && v[k] < lowerBounds[k])
                    v[k] = lowerBounds[k];

                if (upperBounds != null && v[k] > upperBounds[k])
                    v[k] = upperBounds[k];
            }
        }

        /// <summary>
        /// Calculates available space to the left of a parameter relative to its lower bound.
        /// </summary>
        /// <param name="x">Parameter vector.</param>
        /// <param name="j">Index of the parameter to check.</param>
        /// <param name="lowerBounds">Optional lower bounds.</param>
        /// <returns>Distance to lower bound, or positive infinity if unbounded.</returns>
        private static double AvailableLeft(double[] x, int j, double[]? lowerBounds)
        {
            return lowerBounds == null ? double.PositiveInfinity : x[j] - lowerBounds[j];
        }

        /// <summary>
        /// Calculates available space to the right of a parameter relative to its upper bound.
        /// </summary>
        /// <param name="x">Parameter vector.</param>
        /// <param name="j">Index of the parameter to check.</param>
        /// <param name="upperBounds">Optional upper bounds.</param>
        /// <returns>Distance to upper bound, or positive infinity if unbounded.</returns>
        private static double AvailableRight(double[] x, int j, double[]? upperBounds)
        {
            return upperBounds == null ? double.PositiveInfinity : upperBounds[j] - x[j];
        }

        /// <summary>
        /// Checks if a vector contains any non-finite values (NaN or Infinity).
        /// </summary>
        /// <param name="v">Vector to check.</param>
        /// <param name="expectedLength">Expected length of the vector.</param>
        /// <returns>True if the vector is null, has wrong length, or contains non-finite values.</returns>
        private static bool IsBad(double[] v, int expectedLength)
        {
            if (v == null || v.Length != expectedLength)
                return true;

            for (int i = 0; i < expectedLength; i++)
            {
                if (double.IsNaN(v[i]) || double.IsInfinity(v[i]))
                    return true;
            }

            return false;
        }

        #endregion
    }
}