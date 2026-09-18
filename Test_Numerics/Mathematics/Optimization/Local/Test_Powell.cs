using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the Powell optimization algorithm.
    /// </summary>
    [TestClass]
    public class Test_Powell
    {
        /// <summary>
        /// Test the Powell algorithm with a simple 3-dimensional test function.
        /// </summary>
        [TestMethod]
        public void Test_FXYZ()
        {
            var initial = new double[] { 0.2d, 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new Powell(TestFunctions.FXYZ, 3, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            double x = solution[0];
            double y = solution[1];
            double z = solution[2];
            double validX = 0.125d;
            double validY = 0.2d;
            double validZ = 0.35d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
            Assert.AreEqual(z, validZ, 1E-4);
        }

        /// <summary>
        /// Test the Powell algorithm with the De Jong Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            var initial = new double[] { 1.0d, -1.0d, 2.0d, -2.0d, 1.0d };
            var lower = new double[] { -5.12d, -5.12d, -5.12d, -5.12d, -5.12d };
            var upper = new double[] { 5.12d, 5.12d, 5.12d, 5.12d, 5.12d };
            var solver = new Powell(TestFunctions.DeJong, 5, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0d, 0.0d, 0.0d, 0.0d, 0.0d };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the Powell algorithm with the Sum of Power functions in 3-D.
        /// </summary>
        [TestMethod]
        public void Test_SumOfPowerFunctions()
        {
            var initial = new double[] { 0.5d, -0.5d, 0.5d };
            var lower = new double[] { -1d, -1d, -1d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new Powell(TestFunctions.SumOfPowerFunctions, 3, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0d, 0.0d, 0.0d };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the Powell algorithm with the Rosenbrock Function in 2-D.
        /// </summary>
        [TestMethod]
        public void Test_Rosenbrock()
        {
            var initial = new double[] { 0, 0 };
            var lower = new double[] { -2.048, -2.048 };
            var upper = new double[] { 2.048, 2.048 };
            var solver = new Powell(TestFunctions.Rosenbrock, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 1.0d, 1.0d };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the Powell algorithm with the Booth Function
        /// </summary>
        [TestMethod]
        public void Test_Booth()
        {
            var initial = new double[] { 0.0d, 0.0d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new Powell(TestFunctions.Booth, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 1.0d;
            var validY = 3.0d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the Powell algorithm with the Matyas Function
        /// </summary>
        [TestMethod]
        public void Test_Matyas()
        {
            var initial = new double[] { 1.0d, -1.0d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new Powell(TestFunctions.Matyas, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = 0.0d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the Powell algorithm with the McCormick Function
        /// </summary>
        [TestMethod]
        public void Test_McCormick()
        {
            var initial = new double[] { 0.0d, 0.0d };
            var lower = new double[] { -1.5d, -3d };
            var upper = new double[] { 4d, 4d };
            var solver = new Powell(TestFunctions.McCormick, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -1.9133;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = -0.54719d;
            var validY = -1.54719d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the Powell algorithm with the Beale Function
        /// </summary>
        [TestMethod]
        public void Test_Beale()
        {
            var initial = new double[] { 0.0d, 0.0d };
            var lower = new double[] { -4.5d, -4.5d };
            var upper = new double[] { 4.5d, 4.5d };
            var solver = new Powell(TestFunctions.Beale, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 3.0d;
            var validY = 0.5d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Verifies that Powell's Brent line search reaches a distant minimum without linear bracketing work.
        /// </summary>
        [TestMethod]
        public void Test_DistantLineMinimumUsesGeometricBracket()
        {
            var solver = new Powell(
                x => (x[0] - 50d) * (x[0] - 50d),
                1,
                new[] { 0d },
                new[] { -100d },
                new[] { 100d })
            {
                ComputeHessian = false,
                RecordTraces = false
            };

            solver.Minimize();

            Assert.AreEqual(50d, solver.BestParameterSet.Values[0], 1E-4);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-8);
            Assert.IsLessThan(200, solver.FunctionEvaluations, $"Expected fewer than 200 objective evaluations, but observed {solver.FunctionEvaluations}.");
        }

        /// <summary>
        /// The quadratic used by the feasible-region tests below. Its unconstrained minimum is the point
        /// (20, 20), which is far outside the unit square those tests declare, so the constrained solution
        /// is the corner (1, 1) and the objective falls monotonically toward the centre from there.
        /// </summary>
        /// <param name="x">The point to evaluate.</param>
        /// <returns>The value of the quadratic at <paramref name="x"/>.</returns>
        private static double CornerQuadratic(double[] x)
        {
            return (x[0] - 20d) * (x[0] - 20d) + (x[1] - 20d) * (x[1] - 20d);
        }

        /// <summary>
        /// Test that no point outside the declared bounds is ever passed to the objective function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// This covers both paths that can produce an out-of-box evaluation: the line search drives a
        /// bracketing routine that expands geometrically without regard for the bounds and must be
        /// confined to the feasible step interval, and the extrapolated point built for the direction-set
        /// update is a reflection through the current point, which leaves the box readily even when both
        /// points forming it are inside.
        /// </para>
        /// <para>
        /// Both matter because this class is used as the local solver inside <see cref="MultiStart"/> and
        /// <see cref="MLSL"/>, where the objective function is the global solver's own evaluation routine and
        /// any point it scores can be recorded as the reported solution.
        /// </para>
        /// <para>
        /// The unconstrained minimum lies outside the box in the direction of the corner (1, 1), so an
        /// unrestricted line search runs straight out of the box and the reflection through a corner-bound
        /// iterate lands outside it as well. The iteration count is asserted so that a future change that
        /// converges before the extrapolation is ever built cannot leave that second path untested.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_NoPointOutsideTheBoundsIsEvaluated()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var offending = new List<double[]>();

            double recording(double[] x)
            {
                for (int i = 0; i < x.Length; i++)
                {
                    if (x[i] < lower[i] || x[i] > upper[i])
                    {
                        offending.Add((double[])x.Clone());
                        break;
                    }
                }
                return CornerQuadratic(x);
            }

            var solver = new Powell(recording, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false };
            solver.Minimize();

            Assert.IsGreaterThanOrEqualTo(1, solver.Iterations, "The extrapolated point is only built once a full iteration completes, so this test must run at least one iteration to cover it.");
            string report = offending.Count == 0
                ? "No point was evaluated outside the bounds."
                : $"{offending.Count} points were evaluated outside the bounds, the first being ({offending[0][0]}, {offending[0][1]}).";
            Assert.IsEmpty(offending, report);
        }

        /// <summary>
        /// Test that the reported solution is the constrained corner and that its fitness is the value of the
        /// objective function at that corner.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The unconstrained minimum of this objective is the point (20, 20) with value zero, and the
        /// constrained minimum over the unit square is the corner (1, 1) with value 722. An unrestricted
        /// search reaches the unconstrained minimum, scores it as zero and reports it, which is both
        /// infeasible and wrong by 722.
        /// </para>
        /// <para>
        /// Both halves are asserted exactly. A caller that checks feasibility and then trusts the fitness has
        /// to be able to rely on the two describing the same point.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_ConstrainedCornerIsReportedWithItsOwnFitness()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new Powell(CornerQuadratic, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false };
            solver.Minimize();

            var solution = solver.BestParameterSet.Values;
            Assert.AreEqual(1d, solution[0]);
            Assert.AreEqual(1d, solution[1]);
            Assert.AreEqual(722d, solver.BestParameterSet.Fitness);
            Assert.AreEqual(CornerQuadratic(solution), solver.BestParameterSet.Fitness, "The reported fitness must be the objective function evaluated at the reported point.");
        }

        /// <summary>
        /// Test that maximization stores the negated objective, so that the same consistency check applies.
        /// </summary>
        /// <remarks>
        /// <see cref="Optimizer.Evaluate"/> multiplies the objective by a scale that is minus one while
        /// maximizing, and stores that scaled value as the fitness alongside the point. The stored fitness is
        /// therefore the negated objective on this path, which is what the assertion below pins.
        /// </remarks>
        [TestMethod]
        public void Test_MaximizeStoresTheNegatedObjective()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            double peak(double[] x) => 5d - (x[0] - 0.3d) * (x[0] - 0.3d) - (x[1] - 0.4d) * (x[1] - 0.4d);
            var solver = new Powell(peak, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Maximize();

            var solution = solver.BestParameterSet.Values;
            Assert.AreEqual(0.3d, solution[0], 1E-6);
            Assert.AreEqual(0.4d, solution[1], 1E-6);
            Assert.AreEqual(-peak(solution), solver.BestParameterSet.Fitness, "While maximizing, the stored fitness is the negated objective at the stored point.");
        }

        /// <summary>
        /// Test that a box with no width returns its single feasible point.
        /// </summary>
        /// <remarks>
        /// Every coordinate is pinned, so the feasible step interval collapses to the single step length zero
        /// in every direction and no line search can move. The solver must report the pinned point and the
        /// value there, rather than dividing by a zero interval width or spending its whole iteration budget
        /// discovering that it cannot move.
        /// </remarks>
        [TestMethod]
        public void Test_ZeroWidthBoxReturnsItsOnlyFeasiblePoint()
        {
            var initial = new double[] { 1d, 1d };
            var lower = new double[] { 1d, 1d };
            var upper = new double[] { 1d, 1d };
            double f(double[] x) => (x[0] - 3d) * (x[0] - 3d) + (x[1] + 2d) * (x[1] + 2d);
            var solver = new Powell(f, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Minimize();

            var solution = solver.BestParameterSet.Values;
            Assert.AreEqual(1d, solution[0]);
            Assert.AreEqual(1d, solution[1]);
            Assert.AreEqual(13d, solver.BestParameterSet.Fitness);
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        }

        /// <summary>
        /// Test that an iterate sitting on a bound with the minimum outside it stays on that bound.
        /// </summary>
        /// <remarks>
        /// The starting point is the corner (1, 1) and the unconstrained minimum is (5, 5), so every
        /// coordinate direction that improves the objective points out of the box and the feasible step
        /// interval collapses to a point in each of them. The corner is the constrained minimum, with value
        /// 32, and it must be reported with that value rather than with the value at some point outside.
        /// </remarks>
        [TestMethod]
        public void Test_StartOnABoundWithTheMinimumOutsideStaysOnTheBound()
        {
            var initial = new double[] { 1d, 1d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            double f(double[] x) => (x[0] - 5d) * (x[0] - 5d) + (x[1] - 5d) * (x[1] - 5d);
            var solver = new Powell(f, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Minimize();

            var solution = solver.BestParameterSet.Values;
            Assert.AreEqual(1d, solution[0]);
            Assert.AreEqual(1d, solution[1]);
            Assert.AreEqual(32d, solver.BestParameterSet.Fitness);
        }

        /// <summary>
        /// Test that a minimum lying just inside a bound is found from a start on that bound.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The objective is strictly convex and smooth with its unique minimum at (0.04, 0.5), which is
        /// inside the box. The start is on the lower bound of the first coordinate, so the first bracketing
        /// probe is uphill and the bracketing search reverses back through the start and past the bound.
        /// </para>
        /// <para>
        /// Past a bound the objective handed to the line search is constant, and the minimization accepts a
        /// trial point that merely ties its incumbent, so a bracket carrying that constant region is
        /// collapsed into it and the half holding the minimum is thrown away. The solver then stops on the
        /// bound and reports success at a point whose gradient is nowhere near zero. Trimming the bracket
        /// back to the feasible step interval before minimizing is what prevents that.
        /// </para>
        /// <para>
        /// The assertions are on the solution rather than on the status, because the failure reported
        /// success. A stop that leaves a feasible coordinate step still reducing the objective is the
        /// property being ruled out.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MinimumJustInsideABoundIsFoundFromThatBound()
        {
            var initial = new double[] { 0d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            double f(double[] x) => (x[0] - 0.04d) * (x[0] - 0.04d) + (x[1] - 0.5d) * (x[1] - 0.5d);
            var solver = new Powell(f, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Minimize();

            var solution = solver.BestParameterSet.Values;
            Assert.AreEqual(0.04d, solution[0], 1E-6);
            Assert.AreEqual(0.5d, solution[1], 1E-6);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-10);

            // No feasible coordinate step may still reduce the objective at the reported point.
            double at = f(solution);
            foreach (double h in new[] { 1E-2, 1E-3, 1E-4 })
                for (int i = 0; i < solution.Length; i++)
                    foreach (int sign in new[] { 1, -1 })
                    {
                        var probe = (double[])solution.Clone();
                        probe[i] = Math.Max(lower[i], Math.Min(upper[i], solution[i] + sign * h));
                        Assert.IsGreaterThanOrEqualTo(at - 1E-12, f(probe), $"A feasible step of {sign * h} in parameter {i} still reduces the objective, so the reported point is not a constrained minimum.");
                    }
        }

        /// <summary>
        /// Test that an iterate on the far corner searches back into the box.
        /// </summary>
        /// <remarks>
        /// The start is the upper corner and the minimum is interior, so every improving direction points
        /// back into the box and the feasible step interval is entirely negative. The bracketing routine
        /// picks its direction from a single comparison against a first step of a fixed sign, and a first
        /// step pointing out of the box is clamped back onto the start, which makes that comparison a tie
        /// and leaves the search facing the wrong way. Choosing the sign of the first step from the feasible
        /// interval is what prevents that.
        /// </remarks>
        [TestMethod]
        public void Test_IterateOnTheUpperCornerSearchesBackIntoTheBox()
        {
            var initial = new double[] { 1d, 1d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            double f(double[] x) => (x[0] - 0.4d) * (x[0] - 0.4d) + (x[1] - 0.4d) * (x[1] - 0.4d);
            var solver = new Powell(f, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Minimize();

            var solution = solver.BestParameterSet.Values;
            Assert.AreEqual(0.4d, solution[0], 1E-6);
            Assert.AreEqual(0.4d, solution[1], 1E-6);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-10);
        }

        /// <summary>
        /// Test that a line search which finds nothing better than its starting point keeps that point.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The bracketing search hands the minimization an interval that need not contain the zero step, and
        /// the minimization starts from the middle of the interval it is given, so it can return a point
        /// worse than the one the line search was asked to improve on. Evaluating the zero step first and
        /// keeping it unless the search strictly improves on it makes each line search non-increasing.
        /// </para>
        /// <para>
        /// Goldstein-Price over this box is used because it is strongly multimodal, which is what makes a
        /// line search settle on a point worse than its start. From this start the solver reaches the global
        /// minimum of 3; without the fallback it stops in a basin four orders of magnitude worse.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_LineSearchNeverKeepsAPointWorseThanItsStart()
        {
            var initial = new double[] { -10d, -8d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new Powell(TestFunctions.GoldsteinPrice, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Minimize();

            Assert.AreEqual(3d, solver.BestParameterSet.Fitness, 1E-6);
            var solution = solver.BestParameterSet.Values;
            Assert.AreEqual(0d, solution[0], 1E-4);
            Assert.AreEqual(-1d, solution[1], 1E-4);
        }

        /// <summary>
        /// Test that a line search that fails to improve on its start does not disable a direction of the
        /// direction set.
        /// </summary>
        /// <remarks>
        /// When the search along the average direction fails to strictly improve on the current point, the
        /// fallback keeps the current point. Were the direction scaled by that zero step, the zero vector
        /// would be stored into the direction set by the enclosing direction-set update, its feasible step
        /// interval would collapse to the single point zero on every later iteration, and the direction
        /// would be permanently lost. Eggholder from this corner start reaches that state and then stops
        /// at f = -683.29 with one direction dead, while a direction set that keeps the unscaled direction
        /// continues to f = -715.98. The assertion sits between the two regimes rather than pinning the
        /// trajectory, because the objective is chaotic and its trigonometry varies in the last bits
        /// across target frameworks.
        /// </remarks>
        [TestMethod]
        public void Test_ZeroStepFallbackKeepsTheDirectionSetAlive()
        {
            var initial = new double[] { 512d, -128d };
            var lower = new double[] { -512d, -512d };
            var upper = new double[] { 512d, 512d };
            var solver = new Powell(TestFunctions.Eggholder, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Minimize();

            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsLessThan(-700d, solver.BestParameterSet.Fitness,
                "The run stopped early with a dead direction in the direction set.");
        }

        /// <summary>
        /// Test that a box whose width is the smallest representable number still runs a line search.
        /// </summary>
        /// <remarks>
        /// The feasible step interval here is a single subnormal wide. Halving it to pick a first bracketing
        /// step underflows to zero, and a zero step is rejected by the bracketing routine, which fails the
        /// run on a box that is perfectly legal. Falling back on the interval endpoint itself keeps the step
        /// nonzero. The interval that collapses to a single point is handled before this and does not reach
        /// the step selection.
        /// </remarks>
        [TestMethod]
        public void Test_SubnormalWidthBoxStillRunsALineSearch()
        {
            var initial = new double[] { 0d };
            var lower = new double[] { 0d };
            var upper = new double[] { double.Epsilon };
            double f(double[] x) => (x[0] - 1d) * (x[0] - 1d);
            var solver = new Powell(f, 1, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Minimize();

            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            var solution = solver.BestParameterSet.Values;
            Assert.IsGreaterThanOrEqualTo(lower[0], solution[0], "The solution is below its lower bound.");
            Assert.IsLessThanOrEqualTo(upper[0], solution[0], "The solution is above its upper bound.");
            Assert.AreEqual(1d, solver.BestParameterSet.Fitness, 1E-300);
        }

        /// <summary>
        /// Test that infinite bounds leave the search unconstrained.
        /// </summary>
        /// <remarks>
        /// A caller that does not want to bound a coordinate passes an infinite bound for it. The feasible
        /// step interval for such a coordinate is unbounded on that side, which is no constraint at all, so
        /// the solver must still reach a minimum that lies far from its starting point.
        /// </remarks>
        [TestMethod]
        public void Test_InfiniteBoundsLeaveTheSearchUnconstrained()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { double.NegativeInfinity, double.NegativeInfinity };
            var upper = new double[] { double.PositiveInfinity, double.PositiveInfinity };
            double f(double[] x) => (x[0] - 30d) * (x[0] - 30d) + (x[1] + 40d) * (x[1] + 40d);
            var solver = new Powell(f, 2, initial, lower, upper) { ReportFailure = false, RecordTraces = false, ComputeHessian = false };
            solver.Minimize();

            var solution = solver.BestParameterSet.Values;
            Assert.AreEqual(30d, solution[0], 1E-4);
            Assert.AreEqual(-40d, solution[1], 1E-4);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-8);
        }

    }
}
