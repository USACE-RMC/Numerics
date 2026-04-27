using System;

namespace Numerics.Mathematics.RootFinding
{

    /// <summary>
    /// The secant root-finding algorithm.
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
    /// This class contains a shared function for finding the solution to the equation, where f(x)=0 using the Secant method.
    /// The secant method is a root-finding procedure in numerical analysis that uses a series of roots of secant lines to
    /// better approximate a root of a function.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.</item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Root-finding_algorithm"/>
    /// </description></item>
    /// <item> <description> 
    /// <see href="https://en.wikipedia.org/wiki/Secant_method"/>
    /// </description></item>
    /// <item> <description> 
    /// <see href="http://www.m-hikari.com/ams/ams-2017/ams-53-56-2017/p/hahmAMS53-56-2017.pdf"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class Secant
    {

        /// <summary>
        /// Use the secant method to find a solution of the equation where f(x)=0. 
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

            // Define variables
            bool solutionFound = false;
            double xl, root;
            double x1 = lowerBound;
            double x2 = upperBound;
            double fl = f(x1);
            double fh = f(x2);

            // Pick the bound with the smaller function value as the most recent guess
            if (Math.Abs(fl) < Math.Abs(fh))
            {
                root = x1;
                xl = x2;
                Tools.Swap(ref fl, ref fh);
            }
            else
            {
                xl = x1;
                root = x2;
            }

            // Secant loop
            for (int i = 1; i <= maxIterations; i++)
            {
                double dx = (xl - root) * fh / (fh - fl);
                xl = root;
                fl = fh;
                root += dx;
                fh = f(root);
                if (Math.Abs(dx) < tolerance || fh == 0d)
                {
                    solutionFound = true;
                    break;
                }
            }

            // return results of solver
            if (solutionFound == false && reportFailure == true)
            {
                throw new ArgumentException("Secant method failed to find root.");
            }
            else
            {
                return root;
            }

        }
    }
}