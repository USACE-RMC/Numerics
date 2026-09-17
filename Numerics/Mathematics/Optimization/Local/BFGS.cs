using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// Contains the Broyden-Fletcher-Goldfarb-Shanno (BFGS) optimization algorithm. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This is an iterative method for solving unconstrained nonlinear optimization problems. It gradually improves
    /// an approximation to the Hessian matrix of the loss function, obtained from gradient evaluations via a 
    /// generalized secant method.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// <item><description>
    /// </description></item>
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017.
    /// <item><description>
    /// </description></item>
    /// <see href="https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm"/>
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class BFGS : Optimizer
    {

        /// <summary>
        /// Construct a new BFGS optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values to evaluate.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="gradient">Optional. Function to evaluate the gradient. Default uses finite difference.</param>
        public BFGS(Func<double[], double> objectiveFunction, int numberOfParameters, 
                    IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds, 
                    Func<double[], double[]>? gradient = null) : base(objectiveFunction, numberOfParameters)
        {
            // Check if the length of the initial, lower and upper bounds equal the number of parameters
            if (initialValues.Count != numberOfParameters || lowerBounds.Count != numberOfParameters || upperBounds.Count != numberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(lowerBounds), "The initial values and lower and upper bounds must be the same length as the number of parameters.");
            }
            // Check if the initial values are between the lower and upper values
            for (int j = 0; j < initialValues.Count; j++)
            {
                if (upperBounds[j] < lowerBounds[j])
                {
                    throw new ArgumentOutOfRangeException(nameof(upperBounds), "The upper bound cannot be less than the lower bound.");
                }
                if (initialValues[j] < lowerBounds[j] || initialValues[j] > upperBounds[j])
                {
                    throw new ArgumentOutOfRangeException(nameof(initialValues), "The initial values must be between the upper and lower bounds.");
                }
            }
            InitialValues = initialValues.ToArray();
            LowerBounds = lowerBounds.ToArray();
            UpperBounds = upperBounds.ToArray();
            Gradient = gradient;
        }

        /// <summary>
        /// An array of initial values to evaluate. 
        /// </summary>
        public double[] InitialValues { get; private set; }

        /// <summary>
        /// An array of lower bounds (inclusive) of the interval containing the optimal point. 
        /// </summary>
        public double[] LowerBounds { get; private set; }

        /// <summary>
        /// An array of upper bounds (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double[] UpperBounds { get; private set; }

        /// <inheritdoc />
        protected override double[]? ParameterLowerBounds => LowerBounds;

        /// <inheritdoc />
        protected override double[]? ParameterUpperBounds => UpperBounds;

        /// <summary>
        /// The function for evaluating the gradient of the objective function.
        /// </summary>
        public Func<double[], double[]>? Gradient;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            int n = NumberOfParameters;
            bool cancel = false;
            var x = (double[])InitialValues.Clone();
            double f = EvaluateObjective(x, ref cancel);
            if (cancel) return;
            if (!Tools.IsFinite(f))
                throw new ArgumentException("The initial objective value must be finite.", nameof(ObjectiveFunction));
            var g = EvaluateGradient(x, ref cancel);
            if (cancel) return;

            var inverseHessian = Matrix.Identity(n);
            var projected = new double[n];
            var direction = new double[n];
            double stpmax = 100 * Math.Max(Math.Sqrt(Tools.SumProduct(x, x)), n);

            while (true)
            {
                if (ProjectedGradient(x, g, projected) <= AbsoluteTolerance)
                {
                    // Objective probes (including finite differences) can have a slightly lower rounded
                    // value. Return the accepted point whose convergence was actually established.
                    BestParameterSet = new ParameterSet((double[])x.Clone(), f);
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }
                if (Iterations >= MaxIterations)
                {
                    UpdateStatus(OptimizationStatus.MaximumIterationsReached);
                    return;
                }

                for (int i = 0; i < n; i++)
                {
                    direction[i] = 0;
                    for (int j = 0; j < n; j++) direction[i] -= inverseHessian[i, j] * projected[j];
                }
                MakeFeasible(x, direction);
                double slope = Tools.SumProduct(g, direction);
                if (!Tools.IsFinite(slope) || slope >= 0)
                {
                    inverseHessian = Matrix.Identity(n);
                    for (int i = 0; i < n; i++) direction[i] = -projected[i];
                }

                bool canRestart = false;
                for (int i = 0; i < n; i++) canRestart |= direction[i] != -projected[i];
                bool accepted = LineSearch(x, f, g, direction, stpmax, out var nextX, out double nextF, out var nextG, ref cancel);
                if (!accepted && !cancel && canRestart)
                {
                    // A stale metric can exhaust the search even with a negative slope. As in
                    // classical BFGS implementations (for example R's vmmin), restart the metric
                    // once at this point. Both searches must satisfy the same Wolfe conditions.
                    inverseHessian = Matrix.Identity(n);
                    for (int i = 0; i < n; i++) direction[i] = -projected[i];
                    accepted = LineSearch(x, f, g, direction, stpmax, out nextX, out nextF, out nextG, ref cancel);
                }
                if (cancel) return;
                if (!accepted)
                {
                    UpdateStatus(OptimizationStatus.LineSearchFailed);
                    return;
                }
                Iterations++;

                var step = new double[n];
                var change = new double[n];
                for (int i = 0; i < n; i++)
                {
                    step[i] = nextX[i] - x[i];
                    change[i] = nextG[i] - g[i];
                }
                UpdateInverseHessian(ref inverseHessian, step, change);
                x = nextX;
                f = nextF;
                g = nextG;
            }
        }

        /// <summary>Evaluates an objective trial while retaining only finite incumbents.</summary>
        /// <param name="x">The trial point.</param>
        /// <param name="cancel">The evaluation-budget cancellation flag.</param>
        /// <returns>The scaled objective value, including a non-finite rejection value.</returns>
        private double EvaluateObjective(double[] x, ref bool cancel)
        {
            var incumbent = BestParameterSet;
            double f = Evaluate(x, ref cancel);
            if (!Tools.IsFinite(f)) BestParameterSet = incumbent;
            return f;
        }

        /// <summary>Evaluates and validates a gradient in minimization coordinates.</summary>
        /// <param name="x">The point at which to differentiate.</param>
        /// <param name="cancel">The evaluation-budget cancellation flag.</param>
        /// <returns>A private copy of the scaled gradient.</returns>
        /// <exception cref="ArgumentException">The gradient has an invalid dimension or non-finite component.</exception>
        private double[] EvaluateGradient(double[] x, ref bool cancel)
        {
            double[] g;
            if (Gradient != null)
            {
                var supplied = Gradient(x);
                if (supplied == null || supplied.Length != NumberOfParameters)
                    throw new ArgumentException("The gradient must contain one value per parameter.", nameof(Gradient));
                g = (double[])supplied.Clone();
                for (int i = 0; i < g.Length; i++) g[i] *= functionScale;
            }
            else
            {
                bool stopped = cancel;
                // Finite differences may request more probes after cancellation. Do not spend beyond
                // the budget, and leave its status intact instead of misclassifying the partial gradient.
                g = NumericalDerivative.Gradient(p => stopped ? double.NaN : EvaluateObjective(p, ref stopped),
                    x, LowerBounds, UpperBounds);
                cancel = stopped;
                if (cancel) return g;
            }
            for (int i = 0; i < g.Length; i++)
                if (!Tools.IsFinite(g[i]))
                    throw new ArgumentException("The gradient must contain only finite values.", nameof(Gradient));
            return g;
        }

        /// <summary>Projects the gradient onto feasible descent coordinates and returns its infinity norm.</summary>
        /// <param name="x">The feasible point.</param>
        /// <param name="g">Its objective gradient.</param>
        /// <param name="projected">The projected gradient buffer.</param>
        /// <returns>The largest absolute projected component.</returns>
        private double ProjectedGradient(double[] x, double[] g, double[] projected)
        {
            double norm = 0;
            for (int i = 0; i < x.Length; i++)
            {
                projected[i] = (x[i] <= LowerBounds[i] && g[i] > 0) ||
                               (x[i] >= UpperBounds[i] && g[i] < 0) || LowerBounds[i] == UpperBounds[i] ? 0 : g[i];
                norm = Math.Max(norm, Math.Abs(projected[i]));
            }
            return norm;
        }

        /// <summary>Removes direction components that would immediately leave the feasible box.</summary>
        /// <param name="x">The current point.</param>
        /// <param name="direction">The search direction, modified in place.</param>
        private void MakeFeasible(double[] x, double[] direction)
        {
            for (int i = 0; i < x.Length; i++)
                if ((x[i] <= LowerBounds[i] && direction[i] < 0) || (x[i] >= UpperBounds[i] && direction[i] > 0))
                    direction[i] = 0;
        }

        /// <summary>Applies the inverse BFGS update only when its curvature denominators are reliable.</summary>
        /// <param name="h">The inverse Hessian, reset if arithmetic becomes non-finite.</param>
        /// <param name="s">The accepted parameter step.</param>
        /// <param name="y">The change in gradients.</param>
        /// <remarks>Uses the existing Numerical Recipes symmetric BFGS formula with positive-curvature guards.</remarks>
        private static void UpdateInverseHessian(ref Matrix h, double[] s, double[] y)
        {
            int n = s.Length;
            var hy = new double[n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) hy[i] += h[i, j] * y[j];
            double ys = Tools.SumProduct(y, s), yhy = Tools.SumProduct(y, hy);
            double floor = Math.Sqrt(Tools.DoubleMachineEpsilon) * Math.Sqrt(Tools.SumProduct(y, y)) * Math.Sqrt(Tools.SumProduct(s, s));
            if (!Tools.IsFinite(ys) || !Tools.IsFinite(yhy) || yhy <= 0 || ys <= floor) return;
            var v = new double[n];
            for (int i = 0; i < n; i++) v[i] = s[i] / ys - hy[i] / yhy;
            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                {
                    double value = h[i, j] + s[i] * s[j] / ys - hy[i] * hy[j] / yhy + yhy * v[i] * v[j];
                    if (!Tools.IsFinite(value))
                    {
                        h = Matrix.Identity(n);
                        return;
                    }
                    h[i, j] = h[j, i] = value;
                }
        }

        /// <summary>Searches a feasible ray using strong Wolfe conditions, or sufficient decrease at its bound.</summary>
        /// <param name="x0">The current point.</param>
        /// <param name="f0">Its scaled objective.</param>
        /// <param name="g0">Its scaled gradient.</param>
        /// <param name="p">The feasible direction, scaled in place.</param>
        /// <param name="stpmax">The maximum direction length.</param>
        /// <param name="x">The accepted point, or the starting point on failure.</param>
        /// <param name="f">The objective at the returned point.</param>
        /// <param name="g">The gradient at the returned point.</param>
        /// <param name="cancel">The evaluation-budget cancellation flag.</param>
        /// <returns>Whether a step was accepted.</returns>
        /// <remarks>A bound can truncate a descending ray before Wolfe curvature is attainable; subsequent
        /// convergence still requires the projected gradient tolerance. Interior searches use c1=1e-4 and c2=0.9.</remarks>
        private bool LineSearch(double[] x0, double f0, double[] g0, double[] p, double stpmax,
            out double[] x, out double f, out double[] g, ref bool cancel)
        {
            x = x0; f = f0; g = g0;
            double norm = Math.Sqrt(Tools.SumProduct(p, p));
            if (norm > stpmax)
                for (int i = 0; i < p.Length; i++) p[i] *= stpmax / norm;
            double slope0 = Tools.SumProduct(g0, p);
            if (!Tools.IsFinite(slope0) || slope0 >= 0) return false;

            double limit = double.PositiveInfinity;
            for (int i = 0; i < p.Length; i++)
            {
                if (p[i] > 0) limit = Math.Min(limit, (UpperBounds[i] - x0[i]) / p[i]);
                else if (p[i] < 0) limit = Math.Min(limit, (LowerBounds[i] - x0[i]) / p[i]);
            }
            double alpha = Math.Min(1, limit), previous = 0, fPrevious = f0, slopePrevious = slope0;
            for (int iteration = 0; iteration < 20; iteration++)
            {
                var trial = TrialPoint(x0, p, alpha);
                if (alpha <= 0 || trial.SequenceEqual(x0)) return false;
                double value = EvaluateObjective(trial, ref cancel);
                if (cancel) return false;
                if (!Tools.IsFinite(value) || value > f0 + 1e-4 * alpha * slope0)
                    return Zoom(x0, f0, g0, p, slope0, previous, fPrevious, slopePrevious,
                        alpha, value, double.NaN, out x, out f, out g, ref cancel);

                var gradient = EvaluateGradient(trial, ref cancel);
                if (cancel) return false;
                double slope = Tools.SumProduct(gradient, p);
                if (Math.Abs(slope) <= -0.9 * slope0 || (alpha == limit && slope < 0))
                {
                    x = trial; f = value; g = gradient;
                    return true;
                }
                if (!Tools.IsFinite(slope)) return false;
                // Equal rounded values can still satisfy both Wolfe conditions. Check their
                // gradients before reducing the bracket, without relaxing either condition.
                if (iteration > 0 && value >= fPrevious)
                    return Zoom(x0, f0, g0, p, slope0, previous, fPrevious, slopePrevious,
                        alpha, value, slope, out x, out f, out g, ref cancel);
                if (slope >= 0)
                    return Zoom(x0, f0, g0, p, slope0, alpha, value, slope,
                        previous, fPrevious, slopePrevious, out x, out f, out g, ref cancel);
                previous = alpha;
                fPrevious = value;
                slopePrevious = slope;
                alpha = Math.Min(2 * alpha, limit);
                if (alpha == previous) return false;
            }
            return false;
        }

        /// <summary>Constructs a point on a feasible ray, correcting only boundary roundoff.</summary>
        /// <param name="x0">The ray origin.</param>
        /// <param name="p">The feasible direction.</param>
        /// <param name="alpha">A step no greater than the feasible limit.</param>
        /// <returns>The trial coordinates.</returns>
        private double[] TrialPoint(double[] x0, double[] p, double alpha)
        {
            var x = new double[x0.Length];
            for (int i = 0; i < x.Length; i++)
                x[i] = RepairParameter(x0[i] + alpha * p[i], LowerBounds[i], UpperBounds[i]);
            return x;
        }

        /// <summary>Refines a Wolfe bracket while retaining both endpoint values and available slopes.</summary>
        /// <param name="x0">The initial point.</param>
        /// <param name="f0">Its objective.</param>
        /// <param name="g0">Its gradient.</param>
        /// <param name="p">The feasible search direction.</param>
        /// <param name="slope0">The initial directional derivative.</param>
        /// <param name="low">The endpoint with sufficient decrease.</param>
        /// <param name="fLow">Its function value.</param>
        /// <param name="slopeLow">Its directional derivative.</param>
        /// <param name="high">The other bracket endpoint, which may precede low.</param>
        /// <param name="fHigh">Its function value.</param>
        /// <param name="slopeHigh">Its derivative, or NaN if not evaluated.</param>
        /// <param name="x">The accepted point, or initial point on failure.</param>
        /// <param name="f">The returned point's objective.</param>
        /// <param name="g">The returned point's gradient.</param>
        /// <param name="cancel">The evaluation-budget cancellation flag.</param>
        /// <returns>Whether a Wolfe step was found.</returns>
        /// <remarks>Uses the bracket logic of Nocedal and Wright, Numerical Optimization, algorithm 3.6;
        /// compare SciPy 1.16.2 optimize/_linesearch.py. Interpolation is safeguarded away from both endpoints.</remarks>
        private bool Zoom(double[] x0, double f0, double[] g0, double[] p, double slope0,
            double low, double fLow, double slopeLow, double high, double fHigh, double slopeHigh,
            out double[] x, out double f, out double[] g, ref bool cancel)
        {
            x = x0; f = f0; g = g0;
            double[]? previousTrial = null;
            for (int iteration = 0; iteration < 20; iteration++)
            {
                double alpha = Interpolate(low, fLow, slopeLow, high, fHigh, slopeHigh);
                if (alpha == low || alpha == high) return false;
                var trial = TrialPoint(x0, p, alpha);
                if (trial.SequenceEqual(x0) || (previousTrial != null && trial.SequenceEqual(previousTrial))) return false;
                previousTrial = trial;
                double value = EvaluateObjective(trial, ref cancel);
                if (cancel) return false;
                if (!Tools.IsFinite(value) || value > f0 + 1e-4 * alpha * slope0)
                {
                    high = alpha; fHigh = value; slopeHigh = double.NaN;
                }
                else
                {
                    var gradient = EvaluateGradient(trial, ref cancel);
                    if (cancel) return false;
                    double slope = Tools.SumProduct(gradient, p);
                    if (Math.Abs(slope) <= -0.9 * slope0)
                    {
                        x = trial; f = value; g = gradient;
                        return true;
                    }
                    if (!Tools.IsFinite(slope)) return false;
                    if (value >= fLow)
                    {
                        high = alpha; fHigh = value; slopeHigh = slope;
                        continue;
                    }
                    if (slope * (high - low) >= 0)
                    {
                        high = low; fHigh = fLow; slopeHigh = slopeLow;
                    }
                    low = alpha; fLow = value; slopeLow = slope;
                }
            }
            return false;
        }

        /// <summary>Chooses a safeguarded cubic or quadratic interpolant, falling back to bisection.</summary>
        /// <param name="a">The first bracket endpoint.</param>
        /// <param name="fa">Its function value.</param>
        /// <param name="ga">Its slope.</param>
        /// <param name="b">The second endpoint.</param>
        /// <param name="fb">Its function value.</param>
        /// <param name="gb">Its slope, or NaN when unavailable.</param>
        /// <returns>An interior trial step.</returns>
        private static double Interpolate(double a, double fa, double ga, double b, double fb, double gb)
        {
            double width = b - a;
            double left = Math.Min(a, b) + 0.1 * Math.Abs(width);
            double right = Math.Max(a, b) - 0.1 * Math.Abs(width);
            double candidate = double.NaN;
            if (Tools.IsFinite(gb) && Tools.IsFinite(fb))
            {
                double d1 = ga + gb - 3 * (fb - fa) / width;
                double radicand = d1 * d1 - ga * gb;
                if (radicand >= 0)
                {
                    double d2 = Math.Sign(width) * Math.Sqrt(radicand);
                    candidate = b - width * (gb + d2 - d1) / (gb - ga + 2 * d2);
                }
            }
            if (!Tools.IsFinite(candidate) || candidate <= left || candidate >= right)
                candidate = a - ga * width * width / (2 * (fb - fa - ga * width));
            if (!Tools.IsFinite(candidate) || candidate <= Math.Min(a, b) || candidate >= Math.Max(a, b))
                return a + 0.5 * width;
            // Preserve useful interpolation on very steep objectives while guaranteeing a
            // contraction of at least ten percent. Repeated bisection can exhaust the bracket
            // budget before reaching a perfectly representable, very short Wolfe step.
            return Math.Max(left, Math.Min(right, candidate));
        }

    }
}
