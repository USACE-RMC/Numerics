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

namespace Numerics.Mathematics.RootFinding
{

    /// <summary>
    /// This class contains a shared function for finding the solution to the equation f(x)=0 using the Newton-Raphson method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// In numerical analysis, Newton's method (also known as the Newton–Raphson method), named after Isaac Newton and Joseph Raphson,
    /// is a method for finding successively better approximations to the roots (or zeros) of a real-valued function. In this method, 
    /// one takes an initial guess at the root, then a tangent can be extended from the point evaluated at the inital value. The point
    /// where this tangent crosses the x axis, usually represents an improved estimate of the root. 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.</item>
    /// <item> "Applied Numerical Methods with MATLAB for Engineers and Scientists, Third Edition.", Steven C. Chapra, McGraw-Hill, 2012. </item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Root-finding_algorithm"/>
    /// </description></item>
    /// <item> <description> 
    /// <see href="https://en.wikipedia.org/wiki/Newton%27s_method"/>
    /// </description></item>
    /// <item> <description> 
    /// <see href="http://www.m-hikari.com/ams/ams-2017/ams-53-56-2017/p/hahmAMS53-56-2017.pdf"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class NewtonRaphson
    {

        /// <summary>
        /// Use the basic Newton-Raphson method to find a solution of the equation f(x)=0.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Remarks:
        /// </para>
        /// The basic algorithm aborts immediately if the root leaves the bound interval.
        /// </remarks>
        /// <param name="f">The function to find roots from.</param>
        /// <param name="df">The first derivative of the function to find roots from.</param>
        /// <param name="firstGuess">This is the first guess at the root (c).</param>
        /// <param name="tolerance">Optional. Desired tolerance for both the root and the function value at the root.
        /// The root will be refined until the tolerance is achieved or the maximum number of iterations is reached. Default = 1e-8.</param>
        /// <param name="maxIterations">Optional. Maximum number of iterations. Default = 1000.</param>
        /// <param name="reportFailure">Optional. If set to true, an exception will be thrown if the routine fails to converge.
        /// If set to false, the root from the last iteration will be returned if the routine fails to converge. Default = True.</param>
        /// <returns>
        /// The root to the equation f(x)=0 given the specified tolerance.
        /// </returns>
        public static double Solve(Func<double, double> f, Func<double, double> df, double firstGuess, double tolerance = 1E-8, int maxIterations = 1000, bool reportFailure = true)
        {

            // Define variables
            double root = firstGuess;
            bool solutionFound = false;
            double y;
            double yPrime;
            double x0;
            double x1;
            double eps = Tools.DoubleMachineEpsilon; // Don't want to divide by a number smaller than this

            // Newton-Raphson loop
            for (int i = 1; i <= maxIterations; i++)
            {
                // do Newton's computation
                x0 = root;
                y = f(x0);
                yPrime = df(x0);
                if (Math.Abs(yPrime) < eps)
                {
                    // the denominator is too small
                    solutionFound = false;
                    break;
                }
                x1 = x0 - y / yPrime;
                root = x1;
                // check if the result is within the desired tolerance
                if (Math.Abs(x1 - x0) < tolerance)
                {
                    solutionFound = true;
                    break;
                }
            }

            // return results of solver
            if (solutionFound == false && reportFailure == true)
            {
                throw new ArgumentException("Newton-Raphson method failed to find root.");
            }
            else
            {
                return root;
            }

        }


        /// <summary>
        /// Use the robust Newton-Raphson method to find a solution of the equation f(x)=0.
        /// </summary>
        /// <remarks>
        /// Robust Newton-Raphson method falls back to bisection when over or undershooting the bounds.
        /// </remarks>
        /// <param name="f">The function to find roots from.</param>
        /// <param name="df">The first derivative of the function to find roots from.</param>
        /// <param name="firstGuess">This is the first guess at the root.</param>
        /// <param name="lowerBound">The lower bound of the interval containing the root. Aborts if it leaves the interval.</param>
        /// <param name="upperBound">The upper bound of the interval containing the root. Aborts if it leaves the interval.</param>
        /// <param name="tolerance">Optional. Desired tolerance for both the root and the function value at the root.
        /// The root will be refined until the tolerance is achieved or the maximum number of iterations is reached. Default = 1e-8.</param>
        /// <param name="maxIterations">Optional. Maximum number of iterations. Default = 1000.</param>
        /// <param name="reportFailure">Optional. If set to true, an exception will be thrown if the routine fails to converge.
        /// If set to false, the root from the last iteration will be returned if the routine fails to converge. Default = True.</param>
        /// <returns>
        /// The root to the equation f(x)=0 given the specified tolerance.
        /// </returns>
        public static double RobustSolve(Func<double, double> f, Func<double, double> df, double firstGuess, double lowerBound, double upperBound, double tolerance = 1E-8, int maxIterations = 1000, bool reportFailure = true)
        {

            // Define variables
            double root = firstGuess;
            bool solutionFound = false;
            double xl = lowerBound;
            double xh = upperBound;
            double fl = f(xl);
            double fh = f(xh);

            // validate inputs
            if (upperBound < lowerBound)
            {
                throw new ArgumentOutOfRangeException(nameof(upperBound), "The upper bound (b) cannot be less than the lower bound (a).");
            }
            if (root < lowerBound || root > upperBound)
            {
                throw new ArgumentOutOfRangeException(nameof(firstGuess), "The first guess must be between the upper and lower bound.");
            }
            if (fl * fh >= 0d)
            {
                throw new ArgumentException("Robust Newton-Raphson method failed because the root is not bracketed.");
            }

            if (fl < 0d)
            {
                // Orient the search so that f(lowerBound)<0
                xl = lowerBound;
                xh = upperBound;
            }
            else
            {
                xh = lowerBound;
                xl = upperBound;
            }
            // The "step-size before last"
            double dxold = Math.Abs(upperBound - lowerBound);
            // The last step
            double dx = dxold;
            double F = f(root);
            double dF = df(root);
            // Robust Newton-Raphson loop
            for (int i = 1; i <= maxIterations; i++)
            {
                // Bisect if Newton is out of range, or not decreasing fast enough.
                if ((((root - xh) * dF - F) * ((root - xl) * dF - F) > 0.0)
                    || (Math.Abs(2.0 * F) > Math.Abs(dxold * dF)))
                {
                    dxold = dx;
                    dx = 0.5 * (xh - xl);
                    root = xl + dx;
                    if (xl == root)
                    {
                        // Change in root is negligible.
                        // Bisection step acceptable. 
                        solutionFound = true;
                        break;
                    }
                }
                else
                {
                    dxold = dx;
                    dx = F / dF;
                    double temp = root;
                    root -= dx;
                    if (temp == root)
                    {
                        // Change in root is negligible.
                        // Newton step acceptable. 
                    }
                }
                // Check convergence
                if (Math.Abs(dx) < tolerance)
                {
                    solutionFound = true;
                    break;
                }
                // The one new function evaluation per iteration.
                F = f(root);
                dF = df(root);
                // Maintain the bracket on the root.
                if (F < 0.0)
                    xl = root;
                else
                    xh = root;
            }

            // return results of solver
            if (solutionFound == false && reportFailure == true)
            {
                throw new ArgumentException("Robust Newton-Raphson method failed to find root.");
            }
            else
            {
                return root;
            }

        }

        /// <summary>
        /// Solves a system of nonlinear equations using the Newton-Raphson method.
        /// </summary>
        /// <param name="f">
        /// The system of equations to solve, provided as a function that takes a vector 
        /// and returns a vector of the same dimension. The root is where f(x) = 0.
        /// </param>
        /// <param name="df">
        /// The Jacobian matrix function. Takes a vector x and returns the matrix of 
        /// partial derivatives ∂fᵢ/∂xⱼ evaluated at x.
        /// </param>
        /// <param name="firstGuess">Initial guess for the solution vector.</param>
        /// <param name="tolerance">Optional. Desired tolerance for both the root and the function value at the root.
        /// The root will be refined until the tolerance is achieved or the maximum number of iterations is reached. Default = 1e-8.</param>
        /// <param name="maxIterations"> Maximum number of iterations allowed. Default is 1000. </param>
        /// <param name="reportFailure"> If true, throws an exception on failure. If false, returns the best estimate found. Default is true.
        /// </param>
        /// <returns>The solution vector where f(x) ≈ 0.</returns>
        /// <exception cref="ArgumentException">
        /// Thrown when the Jacobian is singular or the method fails to converge 
        /// (only if reportFailure is true).
        /// </exception>
        /// <remarks>
        /// <para>
        /// The multivariate Newton-Raphson method uses the iterative formula:
        /// x_{n+1} = x_n - J(x_n)^{-1} * f(x_n)
        /// where J is the Jacobian matrix of partial derivatives.
        /// </para>
        /// <para>
        /// <b>Performance:</b> For an n-dimensional system, each iteration requires:
        /// - Evaluation of n function values
        /// - Evaluation of n² Jacobian elements
        /// - Solution of an n×n linear system
        /// Time complexity per iteration: O(n³)
        /// </para>
        /// <para>
        /// <b>Convergence:</b> Quadratic convergence when started near the solution 
        /// and the Jacobian is non-singular. Convergence is not guaranteed for poor 
        /// initial guesses.
        /// </para>
        /// </remarks>
        /// <example>
        /// <code>
        /// // Solve the system:
        /// // f1(x,y) = x² + y² - 4 = 0
        /// // f2(x,y) = x - y = 0
        /// 
        /// Func&lt;Vector, Vector&gt; f = (x) => new Vector(new[] {
        ///     x[0]*x[0] + x[1]*x[1] - 4,
        ///     x[0] - x[1]
        /// });
        /// 
        /// Func&lt;Vector, Matrix&gt; df = (x) => new Matrix(new[,] {
        ///     { 2*x[0], 2*x[1] },
        ///     { 1, -1 }
        /// });
        /// 
        /// var initial = new Vector(new[] { 1.0, 1.0 });
        /// var solution = NewtonRaphson.Solve(f, df, initial);
        /// // solution ≈ (√2, √2)
        /// </code>
        /// </example>
        public static Vector Solve(Func<Vector, Vector> f, Func<Vector, Matrix> df, Vector firstGuess, double tolerance = 1E-8, int maxIterations = 1000, bool reportFailure = true) 
        {

            Vector x = firstGuess.Clone();
            bool converged = false;
            double eps = Tools.DoubleMachineEpsilon;

            for (int iter = 1; iter <= maxIterations; iter++)
            {
                Vector fx = f(x);
                Matrix J = df(x);

                // Check for near-singular Jacobian
                if (Math.Abs(J.Determinant()) < eps)
                {
                    if (reportFailure)
                        throw new ArgumentException("Jacobian is singular or nearly singular.");
                    break;
                }

                // Solve J·Δx = F(x)  →  Δx = J⁻¹·F(x)
                Vector deltaX = J.Inverse().Multiply(fx);

                // Update
                Vector xNew = x - deltaX;

                // Check convergence: step size and residual
                if (deltaX.Norm() < tolerance && f(xNew).Norm() < tolerance)
                {
                    x = xNew;
                    converged = true;
                    break;
                }

                x = xNew;
            }

            if (!converged && reportFailure)
                throw new ArgumentException($"Newton–Raphson did not converge after {maxIterations} iterations.");

            return x;

        }
    
    
    }
}