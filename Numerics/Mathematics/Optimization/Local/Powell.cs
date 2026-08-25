using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// Contains the Powell optimization algorithm. The function need not be differentiable, and no derivatives are taken.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors:</b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This method minimizes the function by bi-directionally searching along search vectors via Brent's method. The minima of
    /// these search vectors are recorded and used to create new search vectors, as the current ones are deleted after use.
    /// The algorithm iterates until no significant improvement is made.
    /// </para>
    /// <para>
    /// <b> Feasible region: </b>
    /// The search is confined to the box defined by <see cref="LowerBounds"/> and <see cref="UpperBounds"/>, and the objective
    /// function is never evaluated outside it. Each line search is restricted to the feasible step interval, which is the set
    /// of step lengths along the search direction that keep every coordinate between its bounds. That interval always contains
    /// zero because the current iterate is feasible, so restricting the search never excludes the current point. The
    /// extrapolated point used by the direction-set update is scored only when the reflection that produces it is itself
    /// feasible; an infeasible extrapolation is treated as no improvement and the direction set is left alone for that
    /// iteration. Projecting the reflection back onto the box instead was measured and rejected, because the projected point
    /// sits on a face where the objective can be far smaller than at the reflection, which makes the acceptance test fire on a
    /// point the test was not derived for and admits direction replacements that slow convergence badly on smooth problems.
    /// These restrictions matter most when this class is the local solver inside a global optimizer, because the objective
    /// function is then the global solver's own evaluation routine and any point it scores can be recorded as the reported
    /// solution.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description> 
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </description></item>
    /// <item><description> 
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017
    /// </description></item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Powell%27s_method"/>
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class Powell : Optimizer
    {
        /// <summary>
        /// Construct a new Powell optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values to evaluate.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public Powell(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
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

        /// <inheritdoc/>
        protected override void Optimize()
        {
            // Set variables
            int i, j, D = NumberOfParameters, ibig;
            bool cancel = false;
            double t, fret, fp, fptt, delta;
            var p = InitialValues.ToArray();
            var pt = new double[D];
            var ptt = new double[D];
            var xi = new double[D]; // Direction vector
            // Set the initial matrix for directions
            // and save the initial point
            var ximat = new double[D, D];
            for (i = 0; i < D; i++)
            {
                ximat[i, i] = 1d;
                pt[i] = p[i];
            }
            // initial function evaluation
            fret = Evaluate(p, ref cancel);
            while (Iterations < MaxIterations)
            {
                fp = fret;
                ibig = 0;
                delta = 0.0; // Will be the biggest function decrease.
                // In each iteration, loop over all directions in the set.
                for (i = 0; i < D; i++)
                {
                    // Copy the direction
                    for (j = 0; j < D; j++) xi[j] = ximat[j, i];
                    fptt = fret;
                    fret = LineMinimization(p, xi, ref cancel);
                    if (cancel == true) return;
                    // And record it if it is the larges decrease so far.
                    if (fptt - fret > delta)
                    {
                        delta = fptt - fret;
                        ibig = i + 1;
                    }
                }
                // Check convergence
                if (CheckConvergence(fp, fret))
                {
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }
                // Construct the extrapolated point and save the average direction moved.
                // Save the old starting point.
                // The reflection through the current point readily leaves the box even when both points are
                // inside it, and this point is scored through the objective function, so it has to be kept
                // feasible. It is skipped rather than projected back: projection moves the point onto a face
                // where the objective can be much smaller than at the reflection, which turns the acceptance
                // test below into a test on a different point and admits direction replacements the test was
                // never meant to admit. See the remarks on this class.
                bool extrapolationIsFeasible = true;
                for (j = 0; j < D; j++)
                {
                    ptt[j] = 2.0 * p[j] - pt[j];
                    if (ptt[j] < LowerBounds[j] || ptt[j] > UpperBounds[j]) extrapolationIsFeasible = false;
                    xi[j] = p[j] - pt[j];
                    pt[j] = p[j];
                }
                // Function evaluated at the extrapolated point. An infeasible extrapolation is treated as no
                // improvement, which is the same outcome the acceptance test reaches for a point that does
                // not beat the value at the start of this iteration.
                fptt = extrapolationIsFeasible ? Evaluate(ptt, ref cancel) : double.PositiveInfinity;
                if (cancel == true) return;
                if (fptt < fp)
                {
                    t = 2.0 * (fp - 2.0 * fret + fptt) * Tools.Sqr(fp - fret - delta) - delta * Tools.Sqr(fp - fptt);
                    if (t < 0.0)
                    {
                        // Move to the minimum of the new direction and save the new direction
                        fret = LineMinimization(p, xi, ref cancel);
                        if (cancel == true) return;
                        for (j = 0; j < D; j++)
                        {
                            ximat[j, ibig - 1] = ximat[j, D - 1];
                            ximat[j, D - 1] = xi[j];
                        }
                    }
                }

                Iterations += 1;
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);

        }

        /// <summary>
        /// Auxiliary line minimization routine.
        /// </summary>
        /// <param name="startPoint">The initial point. Updated in place to the minimizing point along the direction.</param>
        /// <param name="direction">The initial direction. Updated in place to the vector displacement actually taken.</param>
        /// <param name="cancel">Determines if the solver should be canceled.</param>
        /// <returns>
        /// The fitness at the returned point, or <see cref="double.NaN"/> when the solver was canceled.
        /// </returns>
        /// <remarks>
        /// <para>
        /// The search is restricted to the feasible step interval returned by
        /// <see cref="GetFeasibleStepInterval"/>. That interval always contains zero, so restricting the
        /// search never excludes the current point. Three guards enforce the restriction, because
        /// <see cref="BrentSearch.Bracket"/> expands geometrically and overwrites its own bounds, so
        /// constructing the search over a narrower interval would not by itself confine it.
        /// </para>
        /// <list type="number">
        /// <item><description>
        /// The step length is clamped to the feasible interval before the trial point is formed, and each
        /// coordinate of that point is then repaired. The objective seen by <see cref="BrentSearch"/> is
        /// therefore the restriction of the function to the feasible segment, extended as a constant past
        /// each end. No trial point can lie outside the box however far the bracket expands, and the
        /// constant tails stop the bracketing search within a step or two of leaving the interval instead
        /// of letting it run away.
        /// </description></item>
        /// <item><description>
        /// The first bracketing step comes from <see cref="GetBracketingStep"/> and always lands inside the
        /// feasible interval, so the single comparison that chooses the bracketing direction is made
        /// between two genuinely different points rather than on a constant tail.
        /// </description></item>
        /// <item><description>
        /// The interval reported by the bracketing search is trimmed back to the feasible interval before
        /// the minimization runs over it, so the minimization never sees the constant tails at all. This
        /// is required for correctness, not just for tidiness. <see cref="BrentSearch"/> accepts a trial
        /// point whose value merely ties the incumbent, and every point on a constant tail ties, so a
        /// bracket that included one would be collapsed into the tail and the half of it that holds the
        /// minimum would be discarded. Trimming rather than widening keeps the search local: the interval
        /// still comes from the outward expansion that started at the current point. It also means the
        /// accepted step is already inside the feasible interval, so the clamp inside the objective is
        /// inactive for it and the returned fitness is the objective at the returned point.
        /// </description></item>
        /// </list>
        /// <para>
        /// An end of the trimmed interval that is a bound of the feasible interval is evaluated after the
        /// minimization and taken when it is strictly better. A line minimum that is cut off by a bound sits
        /// exactly on that end, and <see cref="BrentSearch"/> stops short of an end by its own tolerance and
        /// never evaluates it, so without this a solution on a face is reported a little inside the face
        /// instead of on it. The comparison is strict so that a flat stretch reaching the boundary cannot
        /// pull the iterate out to the boundary for no improvement.
        /// </para>
        /// <para>
        /// The zero step is evaluated before the search and is kept unless the search strictly improves on
        /// it, so this routine is non-increasing. The enclosing algorithm relies on that: it identifies the
        /// direction of largest decrease by index, and that index is only defined when some direction did
        /// decrease.
        /// </para>
        /// <para>
        /// When the feasible interval collapses to the single point zero, the line search cannot move. That
        /// happens when the direction is identically zero and when the iterate sits on a bound with the
        /// direction pointing out of the box in every constrained coordinate. In that case the value at the
        /// current point is returned and the direction is left unchanged, so a blocked direction is retained
        /// in the direction set and can be retried from a later iterate.
        /// </para>
        /// </remarks>
        private double LineMinimization(double[] startPoint, double[] direction, ref bool cancel)
        {
            // Line-minimization routine, Given an n-dimensional point p[0..n-1] and an n-dimension
            // direction xi[0..n-1], moves and resets p to where the function of functor func(p) takes on
            // a minimum along the direction xi from p, and replaces xi by the actual vector displacement
            // that p was moved. Also returns the value of func at the return location p. This is actually
            // all accomplished by calling the Brent minimize routine.
            int D = NumberOfParameters;
            bool c = cancel;
            GetFeasibleStepInterval(startPoint, direction, out double alphaMin, out double alphaMax);

            // The value at the current point, which is the zero step. Brent starts from the middle of the
            // interval it is given rather than from an endpoint, so it need not evaluate the zero step at
            // all. Taking it here makes it the fallback below, which is what keeps this routine from ever
            // returning a point worse than the one it was given.
            double zeroStep = Evaluate(startPoint, ref c);
            cancel = c;
            if (cancel) return double.NaN;

            // The interval always contains zero, so this is the degenerate case in which no step is
            // feasible in either sense. Report the value at the current point instead of searching.
            if (alphaMin == 0d && alphaMax == 0d) return zeroStep;

            double func(double alpha)
            {
                double step = alpha < alphaMin ? alphaMin : (alpha > alphaMax ? alphaMax : alpha);
                var x = new double[D];
                for (int i = 0; i < D; i++)
                    x[i] = RepairParameter(startPoint[i] + step * direction[i], LowerBounds[i], UpperBounds[i]);
                return Evaluate(x, ref c);
            }

            // Bracket the minimum first. The bracketing search steps outward from the current point and can
            // run past an endpoint of the feasible interval, where the objective above is constant. That
            // constant is what stops the expansion, but the interval it reports must be trimmed back to the
            // feasible interval before the minimization runs over it. See the remarks on this method.
            var bracketing = new BrentSearch(func, 0d, 1d) { RelativeTolerance = RelativeTolerance, AbsoluteTolerance = AbsoluteTolerance };
            bracketing.Bracket(GetBracketingStep(alphaMin, alphaMax));
            cancel = c;
            if (cancel) return double.NaN;
            double lower = bracketing.LowerBound < alphaMin ? alphaMin : (bracketing.LowerBound > alphaMax ? alphaMax : bracketing.LowerBound);
            double upper = bracketing.UpperBound < alphaMin ? alphaMin : (bracketing.UpperBound > alphaMax ? alphaMax : bracketing.UpperBound);

            var brent = new BrentSearch(func, lower, upper) { RelativeTolerance = RelativeTolerance, AbsoluteTolerance = AbsoluteTolerance };
            brent.Minimize();
            cancel = c;
            if (cancel) return double.NaN;
            double xmin = brent.BestParameterSet.Values[0];
            double fmin = brent.BestParameterSet.Fitness;

            // A line minimum that is cut off by a bound sits exactly on the end of the feasible interval, and
            // the minimization above stops short of an end by its own tolerance and never evaluates it. Try
            // an end that is a bound of the feasible interval and take it only when it is strictly better, so
            // that a step limited by a bound lands on the bound rather than just inside it, without a flat
            // stretch being able to pull the iterate all the way out to the boundary for nothing.
            if (upper == alphaMax && upper != 0d)
            {
                double atEnd = func(alphaMax);
                if (atEnd < fmin) { xmin = alphaMax; fmin = atEnd; }
            }
            if (lower == alphaMin && lower != 0d)
            {
                double atEnd = func(alphaMin);
                if (atEnd < fmin) { xmin = alphaMin; fmin = atEnd; }
            }
            cancel = c;
            if (cancel) return double.NaN;

            // Fall back on the zero step unless the search strictly improved on it. Written this way so a
            // search that returned NaN keeps the current point rather than moving to it.
            if (!(fmin < zeroStep))
            {
                xmin = 0d;
                fmin = zeroStep;
            }
            for (int j = 0; j < NumberOfParameters; j++)
            {
                direction[j] *= xmin;
                startPoint[j] += direction[j];
                // Make sure the parameter is within bounds
                startPoint[j] = RepairParameter(startPoint[j], LowerBounds[j], UpperBounds[j]);
            }
            return fmin;
        }

        /// <summary>
        /// Returns the first bracketing step for a line search over the given feasible step interval.
        /// </summary>
        /// <param name="alphaMin">The most negative feasible step length. Must not be positive.</param>
        /// <param name="alphaMax">The most positive feasible step length. Must not be negative.</param>
        /// <returns>A nonzero step length that lies inside the feasible interval.</returns>
        /// <remarks>
        /// <para>
        /// <see cref="BrentSearch.Bracket"/> compares the objective at the origin against the objective one
        /// step away, and only that single comparison decides which way it then expands. The objective handed
        /// to it by <see cref="LineMinimization"/> is constant outside the feasible interval, so a first step
        /// that lands outside carries no information: the comparison sees two equal values, the search never
        /// turns around, and the whole feasible side can be missed. Returning a step that lies inside the
        /// interval is what prevents that.
        /// </para>
        /// <para>
        /// The default magnitude is the same 0.1 the routine has always used, and it is kept whenever the
        /// interval is wide enough to hold it in either sense, so an interior iterate brackets exactly as
        /// before. Only when the interval is narrower than the default in both senses is the magnitude
        /// reduced, to half of the wider side, which is inside the interval and nonzero because the caller
        /// has already handled the interval that collapses to a point.
        /// </para>
        /// </remarks>
        private static double GetBracketingStep(double alphaMin, double alphaMax)
        {
            const double defaultStep = 0.1;
            if (alphaMax >= defaultStep) return defaultStep;
            if (alphaMin <= -defaultStep) return -defaultStep;
            double wider = alphaMax >= -alphaMin ? alphaMax : alphaMin;
            double half = 0.5 * wider;
            // Halving underflows to zero once the wider side is the smallest subnormal, and a zero step is
            // rejected by the bracketing routine. The endpoint itself is nonzero there, and the caller has
            // already returned for the interval that collapses to a point, so it is a usable step.
            return half == 0d ? wider : half;
        }

        /// <summary>
        /// Returns the interval of step lengths along a direction that keeps every coordinate within its bounds.
        /// </summary>
        /// <param name="startPoint">The feasible point the step is taken from.</param>
        /// <param name="direction">The search direction.</param>
        /// <param name="alphaMin">When this method returns, the most negative feasible step length.</param>
        /// <param name="alphaMax">When this method returns, the most positive feasible step length.</param>
        /// <remarks>
        /// <para>
        /// For each coordinate the two ratios <c>(UpperBounds[i] - startPoint[i]) / direction[i]</c> and
        /// <c>(LowerBounds[i] - startPoint[i]) / direction[i]</c> are the step lengths that put that coordinate
        /// exactly on a bound. The larger of the two is an upper limit and the smaller is a lower limit when the
        /// component is positive, and the roles swap when it is negative. A zero component imposes no limit, and
        /// an infinite bound produces an infinite limit, which is also no limit.
        /// </para>
        /// <para>
        /// A direction with no nonzero component cannot move the point, so the interval is reported as the single
        /// point zero rather than as the whole line. Because <paramref name="startPoint"/> is feasible the interval
        /// always contains zero; that is enforced explicitly at the end so that a rounding error in one of the
        /// ratios, or a non-finite ratio, cannot produce an interval that excludes the current point.
        /// </para>
        /// </remarks>
        private void GetFeasibleStepInterval(double[] startPoint, double[] direction, out double alphaMin, out double alphaMax)
        {
            alphaMin = double.NegativeInfinity;
            alphaMax = double.PositiveInfinity;
            bool moves = false;
            for (int i = 0; i < NumberOfParameters; i++)
            {
                double d = direction[i];
                if (d == 0d) continue;
                moves = true;
                double toUpper = (UpperBounds[i] - startPoint[i]) / d;
                double toLower = (LowerBounds[i] - startPoint[i]) / d;
                double high = d > 0d ? toUpper : toLower;
                double low = d > 0d ? toLower : toUpper;
                if (high < alphaMax) alphaMax = high;
                if (low > alphaMin) alphaMin = low;
            }

            if (!moves)
            {
                alphaMin = 0d;
                alphaMax = 0d;
                return;
            }

            // The start point is feasible, so zero is always a feasible step. These comparisons are written
            // so that a NaN limit collapses that side of the interval to zero rather than propagating.
            if (!(alphaMax > 0d)) alphaMax = 0d;
            if (!(alphaMin < 0d)) alphaMin = 0d;
        }

    }
}
