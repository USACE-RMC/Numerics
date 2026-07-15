using System;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The Brent optimization algorithm. The function need not be differentiable, and no derivatives are taken.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This class contains methods for finding the minimum or maximum of a function using Brent's method.
    /// Brent's method is a hybrid root-finding algorithm that combines the bisection method, secant method, and inverse
    /// quadratic interpolation.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Brent%27s_method"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class BrentSearch : Optimizer
    {
        /// <summary>
        /// Construct a new Brent optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="lowerBound">The lower bound (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBound">The upper bound (inclusive) of the interval containing the optimal point.</param>
        public BrentSearch(Func<double, double> objectiveFunction, double lowerBound, double upperBound) : base((x) => objectiveFunction(x[0]), 1)
        {
            // validate inputs
            if (upperBound < lowerBound)
            {
                throw new ArgumentOutOfRangeException("upperBound", "The upper bound cannot be less than the lower bound.");
            }
            LowerBound = lowerBound;
            UpperBound = upperBound;
        }

        /// <summary>
        /// The lower bound (inclusive) of the interval containing the optimal point. 
        /// </summary>
        public double LowerBound { get; private set; }

        /// <summary>
        /// The upper bound (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double UpperBound { get; private set; }

        /// <inheritdoc/>
        protected override void Optimize()
        { 
            // Define variables
            bool cancel = false;
            double ax = LowerBound, bx = 0.5 * (UpperBound + LowerBound), cx = UpperBound;
            // Golder ratio and a small number which protects against trying to achieve 
            // fractional accuracy for a minimum that happens to be exactly zero.
            double CGOLD = 0.381966d;
            double ZEPS = Tools.DoubleMachineEpsilon * 1.0e-3;
            double a, b, d = 0.0, etemp, fu, fv, fw, fx;
            double p, q, r, tol1, tol2, u, v, w, x, xm;
            double e = 0.0;

            a = (ax < cx ? ax : cx);
            b = (ax > cx ? ax : cx);
            x = w = v = bx;
            fw = fv = fx = Evaluate(new double[] { x }, ref cancel);
            for (int i = 1; i <= MaxIterations; i++)
            {
                xm = 0.5 * (a + b);
                tol2 = 2.0 * (tol1 = RelativeTolerance * Math.Abs(x) + ZEPS);
                // Test for done here
                if (Math.Abs(x - xm) <= tol2 - 0.5d * (b - a))
                {
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }
                // Construct a trial parabolic fit
                if (Math.Abs(e) > tol1)
                {
                    r = (x - w) * (fx - fv);
                    q = (x - v) * (fx - fw);
                    p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);
                    if (q > 0.0) p = -p;
                    q = Math.Abs(q);
                    etemp = e;
                    e = d;
                    if (Math.Abs(p) >= Math.Abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
                        d = CGOLD * (e = (x >= xm ? a - x : b - x));
                    // The above conditions determine the acceptability of the parabolic fit. 
                    // here we take the golden section step into the larger of the two segments. 
                    else
                    {
                        d = p / q;
                        u = x + d;
                        if (u - a < tol2 || b - u < tol2)
                            d = Tools.Sign(tol1, xm - x);
                    }
                }
                else
                {
                    d = CGOLD * (e = (x >= xm ? a - x : b - x));
                }

                u = (Math.Abs(d) >= tol1 ? x + d : x + Tools.Sign(tol1, d));
                fu = Evaluate(new double[] { u }, ref cancel);
                if (cancel == true) return;
                // This is the one function evaluation per iteration.
                if (fu <= fx)
                {
                    if (u >= x) a = x; else b = x;
                    v = w;
                    w = x;
                    x = u;
                    fv = fw;
                    fw = fx;
                    fx = fu;
                }
                else
                {
                    if (u < x) a = u; else b = u;
                    if (fu <= fw || w == x)
                    {
                        v = w;
                        w = u;
                        fv = fw;
                        fw = fu;
                    }
                    else if (fu <= fv || v == x || v == w)
                    {
                        v = u;
                        fv = fu;
                    }
                }
                // Done with housekeeping. Back for another iteration
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

        /// <summary>
        /// Bracket the objective function minimum.
        /// </summary>
        /// <param name="s">The finite, nonzero starting step size. A negative value searches toward decreasing coordinates first. Default = 1E-2.</param>
        /// <param name="k">The finite geometric expansion factor, which must be greater than one. Default = 2.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="s"/> is zero or non-finite, or <paramref name="k"/> is non-finite or not greater than one.</exception>
        /// <exception cref="ArgumentException">The bracket is not found within <see cref="Optimizer.MaxIterations"/> and <see cref="Optimizer.ReportFailure"/> is <see langword="true"/>.</exception>
        /// <exception cref="ArithmeticException">The search leaves the finite range of <see cref="double"/> and <see cref="Optimizer.ReportFailure"/> is <see langword="true"/>.</exception>
        /// <exception cref="InvalidOperationException">The objective function returns <see cref="double.NaN"/> and <see cref="Optimizer.ReportFailure"/> is <see langword="true"/>.</exception>
        /// <remarks>
        /// The first three trial coordinates are equally spaced. After each unsuccessful downhill step,
        /// the signed step is multiplied by <paramref name="k"/> until the middle value is no greater
        /// than either endpoint. Failed searches leave the existing bounds unchanged.
        /// </remarks>
        public void Bracket(double s = 1E-2, double k = 2d)
        {
            if (!Tools.IsFinite(s) || s == 0d)
                throw new ArgumentOutOfRangeException(nameof(s), "The starting step size must be finite and nonzero.");
            if (!Tools.IsFinite(k) || k <= 1d)
                throw new ArgumentOutOfRangeException(nameof(k), "The expansion factor must be finite and greater than one.");

            double a = LowerBound, b = a + s;
            if (!Tools.IsFinite(a) || !Tools.IsFinite(b))
            {
                UpdateStatus(OptimizationStatus.Failure, new ArithmeticException("The initial bracketing coordinates must be finite."));
                return;
            }

            double fa = ObjectiveFunction(new double[] { a });
            double fb = ObjectiveFunction(new double[] { b });
            if (double.IsNaN(fa) || double.IsNaN(fb))
            {
                UpdateStatus(OptimizationStatus.Failure, new InvalidOperationException("The objective function returned NaN while initializing the bracket."));
                return;
            }

            if (fb > fa)
            {
                double temp = a;
                a = b;
                b = temp;
                fb = fa;
                s *= -1;
            }

            for (int iteration = 0; iteration < MaxIterations; iteration++)
            {
                double c = b + s;
                if (!Tools.IsFinite(c))
                {
                    UpdateStatus(OptimizationStatus.Failure, new ArithmeticException("The bracketing search exceeded the finite range of double-precision coordinates."));
                    return;
                }

                double fc = ObjectiveFunction(new double[] { c });
                if (double.IsNaN(fc))
                {
                    UpdateStatus(OptimizationStatus.Failure, new InvalidOperationException("The objective function returned NaN while expanding the bracket."));
                    return;
                }

                if (fc >= fb)
                {
                    LowerBound = Math.Min(a, c);
                    UpperBound = Math.Max(a, c);
                    return;
                }

                a = b;
                b = c;
                fb = fc;
                s *= k;
                if (!Tools.IsFinite(s))
                {
                    UpdateStatus(OptimizationStatus.Failure, new ArithmeticException("The geometric bracketing step exceeded the finite range of double precision."));
                    return;
                }
            }

            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

    }
}