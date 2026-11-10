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
    /// quadratic interpolation. The solver is a literal port of brentq.c in SciPy and returns the same doubles as
    /// scipy.optimize.brentq.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> Brent, R. P. (1973). Algorithms for Minimization without Derivatives. Prentice-Hall.</item>
    /// <item> SciPy brentq.c (scipy/optimize/Zeros/brentq.c); Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods, 17, 261-272.</item>
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
        /// The smallest relative tolerance accepted by <see cref="Solve(Func{double, double}, double, double, double, double, int, bool)"/>:
        /// four times the double-precision machine epsilon 2^-52, or 8.881784197001252E-16.
        /// </summary>
        /// <remarks>
        /// This is the default and the minimum relative tolerance of scipy.optimize.brentq. A smaller relative tolerance
        /// would ask for a bracket narrower than the spacing of doubles near the root.
        /// </remarks>
        public const double MinimumRelativeTolerance = 8.8817841970012523E-16;

        /// <summary>
        /// Use the Brent method to find a solution of the equation where f(x)=0.
        /// </summary>
        /// <param name="f">The function to solve.</param>
        /// <param name="lowerBound">The lower bound (a) of the interval containing the root.</param>
        /// <param name="upperBound">The upper bound (b) of the interval containing the root.</param>
        /// <param name="tolerance">Optional. Absolute tolerance on the root: the search stops when the half-width of the bracket is less than
        /// tolerance / 2 + 2 * 2^-52 * |x|, where x is the current estimate, or the function is exactly zero. Default = 1e-8.</param>
        /// <param name="maxIterations">Optional. Maximum number of iterations. Default = 1000.</param>
        /// <param name="reportFailure">Optional. If set to true, an exception will be thrown if the routine fails to converge.
        /// If set to false, the root from the last iteration will be returned if the routine fails to converge. Default = True.</param>
        /// <returns>
        /// The root to the equation f(x)=0 given the specified tolerance.
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the upper bound is less than the lower bound, or the tolerance is negative or NaN.</exception>
        /// <exception cref="ArgumentException">Thrown when the root is not bracketed, or when the method fails to converge and
        /// <paramref name="reportFailure"/> is true.</exception>
        /// <remarks>
        /// This overload calls <see cref="Solve(Func{double, double}, double, double, double, double, int, bool)"/> with the relative tolerance
        /// <see cref="MinimumRelativeTolerance"/>, the default of scipy.optimize.brentq, so the result is bitwise equal to
        /// <c>scipy.optimize.brentq(f, lowerBound, upperBound, xtol=tolerance, maxiter=maxIterations)</c>.
        /// </remarks>
        public static double Solve(Func<double, double> f, double lowerBound, double upperBound, double tolerance = 1E-8, int maxIterations = 1000, bool reportFailure = true)
        {
            return Solve(f, lowerBound, upperBound, tolerance, MinimumRelativeTolerance, maxIterations, reportFailure);
        }

        /// <summary>
        /// Use the Brent method to find a solution of the equation where f(x)=0, with an absolute and a relative tolerance on the root.
        /// </summary>
        /// <param name="f">The function to solve.</param>
        /// <param name="lowerBound">The lower bound (a) of the interval containing the root.</param>
        /// <param name="upperBound">The upper bound (b) of the interval containing the root.</param>
        /// <param name="absoluteTolerance">The absolute tolerance on the root (xtol in SciPy). Must be non-negative.</param>
        /// <param name="relativeTolerance">The relative tolerance on the root (rtol in SciPy). Must be at least <see cref="MinimumRelativeTolerance"/>.</param>
        /// <param name="maxIterations">Optional. Maximum number of iterations. Default = 1000.</param>
        /// <param name="reportFailure">Optional. If set to true, an exception will be thrown if the routine fails to converge.
        /// If set to false, the estimate from the last iteration will be returned if the routine fails to converge. Default = True.</param>
        /// <returns>
        /// The root to the equation f(x)=0 given the specified tolerances.
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the upper bound is less than the lower bound, the absolute tolerance
        /// is negative or NaN, or the relative tolerance is less than <see cref="MinimumRelativeTolerance"/> or NaN.</exception>
        /// <exception cref="ArgumentException">Thrown when the root is not bracketed, or when the method fails to converge and
        /// <paramref name="reportFailure"/> is true.</exception>
        /// <remarks>
        /// <para>
        /// This is a literal port of brentq.c in SciPy (scipy/optimize/Zeros/brentq.c, written by Charles Harris), so the
        /// result is bitwise equal to <c>scipy.optimize.brentq(f, lowerBound, upperBound, xtol=absoluteTolerance,
        /// rtol=relativeTolerance, maxiter=maxIterations)</c>. The bracket [xblk, xcur] always holds a sign change, with
        /// |f(xcur)| no greater than |f(xblk)|. Each step tries inverse linear interpolation (secant) when only two points are
        /// distinct, or inverse quadratic extrapolation otherwise, and falls back to bisection unless twice the trial step is
        /// less than both the step before last and 3 * |sbis| - delta, where sbis is the bisection step. A step no longer than
        /// delta is lengthened to delta toward the other end of the bracket.
        /// </para>
        /// <para>
        /// The search stops when f(xcur) is exactly zero, or when the half-width of the bracket is less than
        /// delta = (absoluteTolerance + relativeTolerance * |xcur|) / 2. One iteration costs one function evaluation, after the two at
        /// the bounds. A zero at either bound is returned before the sign test, which throws when f has the same strict sign at
        /// both bounds. Where f is NaN at a bound the C code compares sign bits and may reject the bracket; this method does not
        /// throw then, and the result is unspecified, so give it a bracket on which f is defined. scipy.optimize.brentq also raises
        /// when f returns NaN anywhere; this method does not check. Two further inputs that SciPy rejects are
        /// accepted: an absolute tolerance of zero (SciPy raises "xtol too small"), which leaves only the relative
        /// term, and a negative maximum iteration count (SciPy raises "maxiter must be >= 0"), which runs no iterations, as zero
        /// does in both, so the call throws, or returns the upper bound when <paramref name="reportFailure"/> is false.
        /// </para>
        /// <para>
        /// References: Brent, R. P. (1973), Algorithms for Minimization without Derivatives, Prentice-Hall, chapter 4; and Virtanen,
        /// P., et al. (2020), SciPy 1.0: fundamental algorithms for scientific computing in Python, Nature Methods, 17, 261-272.
        /// </para>
        /// </remarks>
        public static double Solve(Func<double, double> f, double lowerBound, double upperBound, double absoluteTolerance, double relativeTolerance, int maxIterations = 1000, bool reportFailure = true)
        {
            // validate inputs
            if (upperBound < lowerBound)
            {
                throw new ArgumentOutOfRangeException(nameof(upperBound), "The upper bound (b) cannot be less than the lower bound (a).");
            }
            if (!(absoluteTolerance >= 0d))
            {
                throw new ArgumentOutOfRangeException(nameof(absoluteTolerance), "The absolute tolerance must be non-negative.");
            }
            if (!(relativeTolerance >= MinimumRelativeTolerance))
            {
                throw new ArgumentOutOfRangeException(nameof(relativeTolerance), "The relative tolerance must be at least 4 times the machine epsilon (8.881784197001252E-16).");
            }

            double xpre = lowerBound, xcur = upperBound;
            double xblk = 0d, fpre, fcur, fblk = 0d, spre = 0d, scur = 0d, sbis;
            // The tolerance is 2 * delta.
            double delta;
            double stry, dpre, dblk;

            fpre = f(xpre);
            fcur = f(xcur);
            if (fpre == 0d) return xpre;
            if (fcur == 0d) return xcur;
            // brentq.c rejects equal sign bits; this test agrees with it for every non-NaN pair and
            // lets a NaN at a bound through.
            if ((fpre > 0d && fcur > 0d) || (fpre < 0d && fcur < 0d))
            {
                throw new ArgumentException("Brent's method failed because the root is not bracketed.");
            }

            for (int i = 0; i < maxIterations; i++)
            {
                if (fpre != 0d && fcur != 0d && SignBit(fpre) != SignBit(fcur))
                {
                    xblk = xpre;
                    fblk = fpre;
                    spre = scur = xcur - xpre;
                }
                if (Math.Abs(fblk) < Math.Abs(fcur))
                {
                    xpre = xcur;
                    xcur = xblk;
                    xblk = xpre;

                    fpre = fcur;
                    fcur = fblk;
                    fblk = fpre;
                }

                delta = (absoluteTolerance + relativeTolerance * Math.Abs(xcur)) / 2d;
                sbis = (xblk - xcur) / 2d;
                if (fcur == 0d || Math.Abs(sbis) < delta)
                {
                    return xcur;
                }

                if (Math.Abs(spre) > delta && Math.Abs(fcur) < Math.Abs(fpre))
                {
                    if (xpre == xblk)
                    {
                        // Interpolate.
                        stry = -fcur * (xcur - xpre) / (fcur - fpre);
                    }
                    else
                    {
                        // Extrapolate.
                        dpre = (fpre - fcur) / (xpre - xcur);
                        dblk = (fblk - fcur) / (xblk - xcur);
                        stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre));
                    }
                    double bisectionBound = 3d * Math.Abs(sbis) - delta;
                    if (2d * Math.Abs(stry) < (Math.Abs(spre) < bisectionBound ? Math.Abs(spre) : bisectionBound))
                    {
                        // Good short step.
                        spre = scur;
                        scur = stry;
                    }
                    else
                    {
                        // Bisect.
                        spre = sbis;
                        scur = sbis;
                    }
                }
                else
                {
                    // Bisect.
                    spre = sbis;
                    scur = sbis;
                }

                xpre = xcur;
                fpre = fcur;
                if (Math.Abs(scur) > delta)
                    xcur += scur;
                else
                    xcur += (sbis > 0d ? delta : -delta);

                fcur = f(xcur);
            }

            if (reportFailure)
            {
                throw new ArgumentException("Brent's method failed to find root.");
            }
            return xcur;
        }

        /// <summary>
        /// Returns true when the sign bit of the value is set, as the C function signbit does, including for -0 and a negative NaN.
        /// </summary>
        /// <param name="value">The value to test.</param>
        /// <returns>True when the sign bit is set.</returns>
        private static bool SignBit(double value)
        {
            return BitConverter.DoubleToInt64Bits(value) < 0L;
        }

        /// <summary>
        /// Bracket the root by expanding outward.
        /// </summary>
        /// <param name="f">The function to solve.</param>
        /// <param name="lowerBound">The lower bound (a) of the interval containing the root.</param>
        /// <param name="upperBound">The upper bound (b) of the interval containing the root.</param>
        /// <param name="f1">Output. The function value at the lower bound.</param>
        /// <param name="f2">Output. The function value at the upper bound.</param>
        /// <param name="maxIterations">Optional. Maximum number of expansion steps. Default = 10.</param>
        /// <returns>
        /// True when the updated bounds bracket a sign change, that is <c>f1 * f2 &lt; 0</c>; otherwise false.
        /// </returns>
        /// <exception cref="ArgumentException">Thrown when the lower and upper bounds are equal.</exception>
        /// <remarks>
        /// Each step moves the endpoint with the smaller absolute function value outward by 1.6 times the
        /// current width, and both bounds are updated in place. The sign test is applied before each step and
        /// once more after the final step.
        /// </remarks>
        public static bool Bracket(Func<double, double> f, ref double lowerBound, ref double upperBound, out double f1, out double f2, int maxIterations = 10)
        {
            double FACTOR = 1.6;
            if (lowerBound == upperBound) throw new ArgumentException("Bad initial range in bracket.");
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
            // Test the endpoint produced by the final step as well.
            return f1 * f2 < 0.0;

        }



    }
}