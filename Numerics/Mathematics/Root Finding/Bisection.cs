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

namespace Numerics.Mathematics.RootFinding
{

    /// <summary>
    /// Contains the bisection root-finding algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// This class contains a shared function for finding the solution to the equation f(x)=0 using the bisection method.
    /// Bisection method is a variation of the incremental search method in which the interval is always divided in half.  
    /// If a function changes sign over an interval, the function value at the midpoint is evaluated. The location 
    /// of the root is then determined as lying within the interval for the next iteration. 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes: The art of Scientific Computing, Third Edition." Press et al. 2017. </item>
    /// <item>"Applied Numerical Methods with MATLAB for Engineers and Scientists, Third Edition.", Steven C. Chapra, McGraw-Hill, 2012. </item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Root-finding_algorithm"/>
    /// </description></item>
    /// <item> <description> 
    /// <see href="https://en.wikipedia.org/wiki/Bisection_method"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class Bisection
    {

        /// <summary>
        /// Use the bisection method to find the solution of the equation where f(x)=0. 
        /// </summary>
        /// <param name="f">The function to solve.</param>
        /// <param name="firstGuess">This is the first guess at the root.</param>
        /// <param name="lowerBound">The lower bound of the interval containing the root.</param>
        /// <param name="upperBound">The upper bound of the interval containing the root.</param>
        /// <param name="tolerance">Optional. Desired tolerance for both the root and the function value at the root.
        /// The root will be refined until the tolerance is achieved or the maximum number of iterations is reached. Default = 1e-8.</param>
        /// <param name="maxIterations">Optional. Maximum number of iterations. Default = 1000.</param>
        /// <param name="reportFailure">Optional. If set to true, an exception will be thrown if the routine fails to converge.
        /// If set to false, the root from the last iteration will be returned if the routine fails to converge. Default = True.</param>
        /// <returns>
        /// The root to the equation f(x)=0 given the specified tolerance.
        /// </returns>
        public static double Solve(Func<double, double> f, double firstGuess, double lowerBound, double upperBound, double tolerance = 1E-8, int maxIterations = 1000, bool reportFailure = true)
        {

            // INPUT: Function f(), endpoint values a, b, tolerance TOL, maximum iterations NMAX
            // CONDITIONS: a < b, either f(a) < 0 and f(b) > 0 or f(a) > 0 and f(b) <0
            // OUTPUT: value which differs from a root of f(x)=0 by less than TOL

            // Define variables
            double root = firstGuess;
            bool solutionFound = false;
            double xl = lowerBound;
            double xh = upperBound;
            double fl = f(xl);
            double fmid = f(xh);
            double dx, xmid;

            // validate inputs
            if (upperBound < lowerBound)
            {
                throw new ArgumentOutOfRangeException(nameof(upperBound), "The upper bound (b) cannot be less than the lower bound (a).");
            }
            if (root < lowerBound || root > upperBound)
            {
                throw new ArgumentOutOfRangeException(nameof(firstGuess), "The first guess must be between the upper and lower bound.");
            }
            if (fl * fmid >=0.0)
            {
                throw new ArgumentException("Bisection method failed because the root is not bracketed.");
            }
            if (fl < 0d)
            {
                // Orient the search so that f>0 lies at X+DX
                dx = xh - xl;
            }
            else
            {
                dx = xl - xh;
            }
            // Bisection loop
            for (int i = 1; i <= maxIterations; i++)
            {
                fmid = f(xmid = root + (dx *= 0.5));
                if (fmid <= 0) root = xmid;
                // check if the solution meets required tolerance
                if (Math.Abs(dx) <= tolerance || fmid == 0.0)
                {
                    // a solution has been achieved, so exit loop
                    solutionFound = true;
                    break;
                }
            }

            // return results of solver
            if (solutionFound == false && reportFailure == true)
            {
                throw new ArgumentException("Bisection method failed to find root.");
            }
            else
            {
                return root;
            }
        }
    }
}