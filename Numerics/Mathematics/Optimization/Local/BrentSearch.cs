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
    /// This class contains methods for finding the minimum or maximum of a function of one variable on a bounded interval
    /// using Brent's method, which combines golden-section search with successive parabolic interpolation. The search is a
    /// port of the bounded minimizer of SciPy, scipy.optimize.minimize_scalar(method='bounded').
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Brent, R. P. (1973). Algorithms for Minimization without Derivatives. Prentice-Hall.
    /// </description></item>
    /// <item><description>
    /// Forsythe, G. E., Malcolm, M. A., and Moler, C. B. (1977). Computer Methods for Mathematical Computations. Prentice-Hall.
    /// </description></item>
    /// <item><description>
    /// Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nature Methods, 17, 261-272.
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
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the upper bound is less than the lower bound.</exception>
        /// <remarks>
        /// The constructor sets <see cref="Optimizer.AbsoluteTolerance"/> to 3 * (1E-3 * <see cref="Tools.DoubleMachineEpsilon"/>),
        /// about 3.3E-19, so the absolute term of the step tolerance, <see cref="Optimizer.AbsoluteTolerance"/> / 3, is
        /// about 1.1E-19: a guard against seeking relative accuracy at a minimum that is exactly zero.
        /// </remarks>
        public BrentSearch(Func<double, double> objectiveFunction, double lowerBound, double upperBound) : base((x) => objectiveFunction(x[0]), 1)
        {
            // validate inputs
            if (upperBound < lowerBound)
            {
                throw new ArgumentOutOfRangeException("upperBound", "The upper bound cannot be less than the lower bound.");
            }
            LowerBound = lowerBound;
            UpperBound = upperBound;
            AbsoluteTolerance = 3d * (Tools.DoubleMachineEpsilon * 1E-3);
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
        /// <remarks>
        /// <para>
        /// The search is a port of the bounded Brent minimizer of SciPy (<c>_minimize_scalar_bounded</c> in
        /// scipy/optimize/_optimize.py, used by <c>scipy.optimize.minimize_scalar(method='bounded')</c> and
        /// <c>scipy.optimize.fminbound</c>), after Forsythe, Malcolm and Moler (1977) and Brent (1973). It starts at the
        /// golden-section point a + (3 - sqrt(5)) / 2 * (b - a) and combines parabolic interpolation with golden-section
        /// steps, never evaluating closer than tol1 to the best point or the bounds.
        /// </para>
        /// <para>
        /// The tolerances map onto SciPy as follows: tol1 = <see cref="Optimizer.RelativeTolerance"/> * |x| +
        /// <see cref="Optimizer.AbsoluteTolerance"/> / 3, where SciPy uses tol1 = sqrt(2.2E-16) * |x| + xatol / 3. The
        /// search stops with <see cref="OptimizationStatus.Success"/> when |x - m| &lt;= 2 * tol1 - (b - a) / 2, with m the
        /// midpoint of the current interval [a, b]. SciPy's maxiter counts function evaluations and maps onto
        /// <see cref="Optimizer.MaxFunctionEvaluations"/>; <see cref="Optimizer.MaxIterations"/> bounds the evaluations after the
        /// first. The defaults, <see cref="Optimizer.RelativeTolerance"/> = 1E-8 and the absolute tolerance set by the
        /// constructor, give tol1 = 1E-8 * |x| + 1.1E-19. With
        /// <see cref="Optimizer.RelativeTolerance"/> = Math.Sqrt(2.2E-16), <see cref="Optimizer.AbsoluteTolerance"/> =
        /// 1E-5 and <see cref="Optimizer.MaxFunctionEvaluations"/> = 500, the evaluated points and the minimizer are bitwise equal
        /// to SciPy's default settings. A NaN objective value is ranked as positive infinity, as everywhere in
        /// <see cref="Optimizer"/>.
        /// </para>
        /// </remarks>
        protected override void Optimize()
        {
            bool cancel = false;
            double relativeTolerance = RelativeTolerance;
            double absoluteTolerance = AbsoluteTolerance;
            double goldenMean = 0.5 * (3.0 - Math.Sqrt(5.0));
            double a = LowerBound, b = UpperBound;
            // fulc, nfc and xf are the third best, second best and best points so far, with values ffulc, fnfc and fx.
            double fulc = a + goldenMean * (b - a);
            double nfc = fulc, xf = fulc;
            double rat = 0.0, e = 0.0;
            double x = xf;
            // One trial buffer for the search: Evaluate copies any point it keeps.
            var point = new double[] { x };
            double fx = Evaluate(point, ref cancel);
            if (cancel) return;
            double ffulc = fx, fnfc = fx;
            double xm = 0.5 * (a + b);
            double tol1 = relativeTolerance * Math.Abs(xf) + absoluteTolerance / 3.0;
            double tol2 = 2.0 * tol1;

            while (Math.Abs(xf - xm) > (tol2 - 0.5 * (b - a)))
            {
                bool golden = true;
                // Check for a parabolic fit.
                if (Math.Abs(e) > tol1)
                {
                    golden = false;
                    double r = (xf - nfc) * (fx - ffulc);
                    double q = (xf - fulc) * (fx - fnfc);
                    double p = (xf - fulc) * q - (xf - nfc) * r;
                    q = 2.0 * (q - r);
                    if (q > 0.0) p = -p;
                    q = Math.Abs(q);
                    r = e;
                    e = rat;

                    // Check for acceptability of the parabola.
                    if (Math.Abs(p) < Math.Abs(0.5 * q * r) && p > q * (a - xf) && p < q * (b - xf))
                    {
                        rat = (p + 0.0) / q;
                        x = xf + rat;
                        if ((x - a) < tol2 || (b - x) < tol2)
                        {
                            // A step of tol1 toward the midpoint, upward when x is at the midpoint.
                            rat = Tools.Sign(tol1, xm - xf);
                        }
                    }
                    else
                    {
                        golden = true;
                    }
                }

                if (golden)
                {
                    // Golden-section step into the larger of the two segments.
                    e = xf >= xm ? a - xf : b - xf;
                    rat = goldenMean * e;
                }

                // The step is at least tol1 long, and upward when rat is zero.
                x = xf + Tools.Sign(Math.Max(Math.Abs(rat), tol1), rat);
                point[0] = x;
                double fu = Evaluate(point, ref cancel);
                Iterations += 1;
                if (cancel) return;

                if (fu <= fx)
                {
                    if (x >= xf) a = xf; else b = xf;
                    fulc = nfc; ffulc = fnfc;
                    nfc = xf; fnfc = fx;
                    xf = x; fx = fu;
                }
                else
                {
                    if (x < xf) a = x; else b = x;
                    if (fu <= fnfc || nfc == xf)
                    {
                        fulc = nfc; ffulc = fnfc;
                        nfc = x; fnfc = fu;
                    }
                    else if (fu <= ffulc || fulc == xf || fulc == nfc)
                    {
                        fulc = x; ffulc = fu;
                    }
                }

                xm = 0.5 * (a + b);
                tol1 = relativeTolerance * Math.Abs(xf) + absoluteTolerance / 3.0;
                tol2 = 2.0 * tol1;

                if (Iterations >= MaxIterations)
                {
                    UpdateStatus(OptimizationStatus.MaximumIterationsReached);
                    return;
                }
            }

            UpdateStatus(OptimizationStatus.Success);
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