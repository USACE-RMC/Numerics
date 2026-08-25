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
        /// This covers both places the solver used to leave the box. The line search used to hand the step
        /// length straight to a bracketing routine that expands geometrically without regard for the bounds,
        /// and the extrapolated point built for the direction-set update is a reflection through the current
        /// point, which leaves the box readily even when both points forming it are inside.
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
