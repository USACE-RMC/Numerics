using System;

namespace Numerics.Mathematics.RootFinding
{

    /// <summary>
    /// The Brent root-finding algorithm.
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
    /// This class contains a shared function for finding the solution to the equation f(x)=0 using Brent's method. Brent's method is a 
    /// hybrid root-finding method that combines the use of bracketing methods, namely, the bisection method, secant method, and inverse 
    /// quadratic interpolation. 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.</item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Root-finding_algorithm"/>
    /// </description></item>
    /// <item> <description> 
    /// <see href="https://en.wikipedia.org/wiki/Brent%27s_method"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class Brent
    {

        /// <summary>
        /// Use the Brent method to find a solution of the equation where f(x)=0.
        /// </summary>
        /// <param name="f">The function to solve.</param>
        /// <param name="lowerBound">The lower bound (a) of the interval containing the root.</param>
        /// <param name="upperBound">The upper bound (b) of the interval containing the root.</param>
        /// <param name="tolerance">Optional. Desired tolerance for both the root and the function value at the root.
        /// The root will be refined until the tolerance is achieved or the maximum number of iterations is reached. Default = 1e-8.</param>
        /// <param name="maxIterations">Optional. Maximum number of iterations. Default = 1000.</param>
        /// <param name="reportFailure">Optional. If set to true, an exception will be thrown if the routine fails to converge.
        /// If set to false, the root from the last iteration will be returned if the routine fails to converge. Default = True.</param>
        /// <returns>
        /// The root to the equation f(x)=0 given the specified tolerance.
        /// </returns>
        public static double Solve(Func<double, double> f, double lowerBound, double upperBound, double tolerance = 1E-8, int maxIterations = 1000, bool reportFailure = true)
        {

            // validate inputs
            if (upperBound < lowerBound)
            {
                throw new ArgumentOutOfRangeException(nameof(upperBound), "The upper bound (b) cannot be less than the lower bound (a).");
            }

            // This algorithm is implemented as shown in:
            // "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.

            // Define variables
            bool solutionFound = false;
            double EPS = Tools.DoubleMachineEpsilon;
            double a = lowerBound, b = upperBound, c = upperBound, d = 0, e = 0, fa = f(a), fb = f(b), fc, p, q, r, s, tol1, xm, root = 0;
            //Verifying opposite sign change in interval
            if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
            {
                throw new ArgumentException("Brent's method failed because the root is not bracketed.");
            }
            fc = fb; 

            // Brent's loop
            for (int i = 1; i <= maxIterations; i++)
            {
                //If both upper bounds are same sign switch c lower bound then rearrange points
                if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
                {
                    c = a;
                    fc = fa;
                    e = d = b - a;
                }

                //So that c can be regarded as the current approximate solution
                if (Math.Abs(fc) < Math.Abs(fb))
                {
                    a = b;
                    b = c;
                    c = a;
                    fa = fb;
                    fb = fc;
                    fc = fa;
                }
                tol1 = 2.0 * EPS * Math.Abs(b) + 0.5 * tolerance; 
                xm = 0.5 * (c - b); //Termination test and possible exit

                //Breaks if b is proven to be root
                if (Math.Abs(xm) <= tol1 || fb == 0.0)
                {
                    root = b;
                    solutionFound = true;
                    break;
                }
                //Chooses open methods or bisection
                if (Math.Abs(e) >= tol1 && Math.Abs(fa) > Math.Abs(fb))
                {
                    s = fb / fa;
                    // Secant Method
                    if (a == c)
                    {
                        p = 2.0 * xm * s;
                        q = 1.0 - s;
                    }
                    // Inverse quadratic interpolation
                    else 
                    {
                        q = fa / fc;
                        r = fb / fc;
                        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                    }
                    if (p > 0.0) q = -q;
                    p = Math.Abs(p);
                    double min1 = 3.0 * xm * q - Math.Abs(tol1 * q);
                    double min2 = Math.Abs(e * q);
                    if (2.0 * p < (min1 < min2 ? min1 : min2))
                    {
                        e = d;
                        d = p / q;
                    }
                    else
                    {
                        d = xm;
                        e = d;
                    }
                }
                else
                {
                    d = xm;
                    e = d;
                }
                //Bisection Method
                a = b;
                fa = fb;
                if (Math.Abs(d) > tol1)
                    b += d;
                else
                    b += Tools.Sign(tol1, xm);
                fb = f(b);
            }

            // return results of solver
            if (solutionFound == false && reportFailure == true)
            {
                throw new ArgumentException("Brent's method failed to find root.");
            }
            else
            {
                return solutionFound ? root : b;
            }
        }

        /// <summary>
        /// Bracket the root by expanding outward.
        /// </summary>
        /// <param name="f">The function to solve.</param>
        /// <param name="lowerBound">The lower bound (a) of the interval containing the root.</param>
        /// <param name="upperBound">The upper bound (b) of the interval containing the root.</param>
        /// <param name="f1">Output. The function value at the lower bound.</param>
        /// <param name="f2">Output. The function value at the upper bound.</param>
        /// <param name="maxIterations">Optional. Maximum number of iterations. Default = 1000.</param>
        public static bool Bracket(Func<double, double> f, ref double lowerBound, ref double upperBound, out double f1, out double f2, int maxIterations = 10)
        {
            double FACTOR = 1.6;
            if (lowerBound == upperBound) throw new Exception("Bad initial range in bracket.");
            f1 = f(lowerBound);
            f2 = f(upperBound);
            
            for (int j = 0; j < maxIterations; j++)
            {
                if (f1 * f2 < 0.0) return true;
                if (Math.Abs(f1) < Math.Abs(f2))
                    f1 = f(lowerBound += FACTOR * (lowerBound - upperBound));
                else
                    f2 = f(upperBound += FACTOR * (upperBound - lowerBound));
            }
            return false;

        }



    }
}