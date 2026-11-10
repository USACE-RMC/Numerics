using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the Multi-Start (MS) optimization algorithm
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_MultiStart
    {
        /// <summary>
        /// Test the MS algorithm with a simple 3-dimensional test function.
        /// </summary>
        [TestMethod]
        public void Test_FXYZ()
        {
            var initial = new double[] { 0.2d, 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new MultiStart(TestFunctions.FXYZ, 3, initial, lower, upper);
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
        /// Test the MS algorithm with the De Jong Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            var initial = new double[] { 1.0d, -1.0d, 2.0d, -2.0d, 1.0d };
            var lower = new double[] { -5.12d, -5.12d, -5.12d, -5.12d, -5.12d };
            var upper = new double[] { 5.12d, 5.12d, 5.12d, 5.12d, 5.12d };
            var solver = new MultiStart(TestFunctions.DeJong, 5, initial, lower, upper);
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
        /// Test the MS algorithm with the Sum of Power functions in 3-D.
        /// </summary>
        [TestMethod]
        public void Test_SumOfPowerFunctions()
        {
            var initial = new double[] { 0.5d, -0.5d, 0.5d };
            var lower = new double[] { -1d, -1d, -1d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new MultiStart(TestFunctions.SumOfPowerFunctions, 3, initial, lower, upper);
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
        /// Test the MS algorithm with the Rosenbrock Function
        /// </summary>
        [TestMethod]
        public void Test_Rosenbrock()
        {
            var initial = new double[] { 0, 0, 0, 0, 0 };
            var lower = new double[] { -2.048, -2.048, -2.048, -2.048, -2.048 };
            var upper = new double[] { 2.048, 2.048, 2.048, 2.048, 2.048 };
            var solver = new MultiStart(TestFunctions.Rosenbrock, 5, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the MS algorithm with the Booth Function
        /// </summary>
        [TestMethod]
        public void Test_Booth()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new MultiStart(TestFunctions.Booth, 2, initial, lower, upper);
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
        /// Test the MS algorithm with the Matyas Function
        /// </summary>
        [TestMethod]
        public void Test_Matyas()
        {
            var initial = new double[] { 2d, -2d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new MultiStart(TestFunctions.Matyas, 2, initial, lower, upper);
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
        /// Test the MS algorithm with the McCormick Function
        /// </summary>
        [TestMethod]
        public void Test_McCormick()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -1.5d, -3d };
            var upper = new double[] { 4d, 4d };
            var solver = new MultiStart(TestFunctions.McCormick, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -1.9133;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = -0.54719d;
            var validY = -1.54719d;
            Assert.AreEqual(x, validX, 1E-3);
            Assert.AreEqual(y, validY, 1E-3);
        }

        /// <summary>
        /// Test the MS algorithm with the Beale Function
        /// </summary>
        [TestMethod]
        public void Test_Beale()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -4.5d, -4.5d };
            var upper = new double[] { 4.5d, 4.5d };
            var solver = new MultiStart(TestFunctions.Beale, 2, initial, lower, upper);
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
        /// Test the MS algorithm with the Goldstein-Price Function
        /// </summary>
        [TestMethod]
        public void Test_GoldsteinPrice()
        {
            var initial = new double[] { -1d, 1d };
            var lower = new double[] { -2d, -2d };
            var upper = new double[] { 2d, 2d };
            var solver = new MultiStart(TestFunctions.GoldsteinPrice, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 3.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = -1.0d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the MS algorithm with the Rastrigin Function
        /// </summary>
        [TestMethod]
        public void Test_Rastrigin()
        {
            var initial = new double[] { 1, 1, 1, 1, 1 };
            var lower = new double[] { -5.12, -5.12, -5.12, -5.12, -5.12 };
            var upper = new double[] { 5.12, 5.12, 5.12, 5.12, 5.12 };
            // Need to run a lot of starts
            var solver = new MultiStart(TestFunctions.Rastrigin, 5, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0 };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the MS algorithm with the Ackley Function
        /// </summary>
        [TestMethod]
        public void Test_Ackley()
        {
            var initial = new double[] { 1d, 1d };
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5d, 5d };
            var solver = new MultiStart(TestFunctions.Ackley, 2, initial, lower, upper);
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
        /// Test the MS algorithm with the three hump camel Function
        /// </summary>
        [TestMethod]
        public void Test_ThreeHumpCamel()
        {
            var initial = new double[] { 2d, -2d };
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5, 5d };
            var solver = new MultiStart(TestFunctions.ThreeHumpCamel, 2, initial, lower, upper);
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
        /// Test the MS algorithm with the Eggholder Function
        /// </summary>
        [TestMethod]
        public void Test_Eggholder()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -512d, -512d };
            var upper = new double[] { 512d, 512d };
            var solver = new MultiStart(TestFunctions.Eggholder, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -959.6407;
            Assert.AreEqual(F, trueF, 0.5);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 512d;
            var validY = 404.2319d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-1);

            // The optimum sits on the upper bound, so the local solver's finite-difference gradient is
            // taken at a point on that bound. The local solver scores its probes through this solver's
            // objective delegate, so an unclamped perturbation can be reported as the solution. Guard the
            // declared feasible region explicitly rather than relying on the tolerances above.
            for (int i = 0; i < solution.Length; i++)
            {
                Assert.IsGreaterThanOrEqualTo(lower[i], solution[i], "Parameter " + i + " is below its lower bound.");
                Assert.IsLessThanOrEqualTo(upper[i], solution[i], "Parameter " + i + " is above its upper bound.");
            }
        }

        /// <summary>
        /// Test the MS algorithm with the tp2 Function
        /// </summary>
        [TestMethod]
        public void Test_TP2()
        {
            var initial = new double[] { 2d, 2d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 2d, 2d };
            var solver = new MultiStart(TestFunctions.tp2, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 1d;
            var validY = 0.666667d;

            bool match1 = Math.Abs(x - validX) < 1E-4 && Math.Abs(y - validY) < 1E-4;
            bool match2 = Math.Abs(x - validY) < 1E-4 && Math.Abs(y - validX) < 1E-4;
            Assert.IsTrue(match1 || match2);
        }

        /// <summary>
        /// Test that the reported solution lies inside the declared bounds when the unconstrained optimum
        /// lies outside them.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The objective is a smooth quadratic centred at (20, 20) restricted to the unit square, so the
        /// constrained solution is the corner (1, 1) and every finite difference taken there steps outside
        /// the feasible region unless it is clamped.
        /// </para>
        /// <para>
        /// The local solvers are built over this solver's own objective delegate, so their probes are scored
        /// through <see cref="Optimizer.Evaluate"/> and can be recorded as the incumbent. Both finite
        /// difference paths matter: the local gradient and the local solver's end-of-run Hessian. The
        /// objective falls monotonically toward the centre, so any unclamped probe scores better than the
        /// corner and would be reported as the solution.
        /// </para>
        /// <para>
        /// <see cref="Optimizer.ComputeHessian"/> is deliberately left at its default of true. Switching it
        /// off would hide the Hessian path, which is the half of this guard that a bounded gradient alone
        /// does not cover.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_SolutionStaysWithinBounds()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new MultiStart(x => Math.Pow(x[0] - 20d, 2d) + Math.Pow(x[1] - 20d, 2d), 2, initial, lower, upper)
            {
                ReportFailure = false,
                // Every local start converges on the same corner, so a handful of starts exercises the
                // guarded path as thoroughly as the default hundred and keeps the test cheap.
                MaxIterations = 10,
                // Bookkeeping only, and it dominates the runtime of this test.
                RecordTraces = false
            };
            solver.Minimize();

            var solution = solver.BestParameterSet.Values;
            Assert.IsNotNull(solution);
            for (int i = 0; i < solution.Length; i++)
            {
                Assert.IsGreaterThanOrEqualTo(lower[i], solution[i], "Parameter " + i + " is below its lower bound.");
                Assert.IsLessThanOrEqualTo(upper[i], solution[i], "Parameter " + i + " is above its upper bound.");
            }

            // The constrained minimum is the corner nearest the unconstrained optimum.
            Assert.AreEqual(1d, solution[0], 1E-6);
            Assert.AreEqual(1d, solution[1], 1E-6);
        }

        /// <summary>
        /// Test that the public <see cref="MultiStart.InitialValues"/> property still equals what was
        /// passed to the constructor after a run, and that the caller's own array is left untouched.
        /// </summary>
        /// <remarks>
        /// The solver's working array is a private copy of <see cref="MultiStart.InitialValues"/>:
        /// uniform draws overwrite only that copy, no sampled point may alias the public array, and
        /// the local-search entry point must not modify its argument. This pins both the public
        /// <see cref="MultiStart.InitialValues"/> array and the caller's own array.
        /// </remarks>
        [TestMethod]
        public void Test_InitialValuesAreNotMutatedByARun()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var callerSnapshot = (double[])initial.Clone();
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new MultiStart(TestFunctions.Booth, 2, initial, lower, upper);
            solver.Minimize();

            CollectionAssert.AreEqual(callerSnapshot, initial, "The caller's own array must not be modified by a run.");
            CollectionAssert.AreEqual(callerSnapshot, solver.InitialValues, "InitialValues must still equal what was passed to the constructor after a run.");
        }

        /// <summary>
        /// An objective with a quadratic well around (0.5, 0.5) and a high plateau everywhere else, so a
        /// point perturbed out of the well cannot re-attain a fitness recorded inside it.
        /// </summary>
        /// <param name="x">The point to evaluate.</param>
        /// <returns>The squared distance to the well centre inside the well; 1000 outside.</returns>
        private static double WellInPlateau(double[] x)
        {
            double dx = x[0] - 0.5d, dy = x[1] - 0.5d;
            double r2 = dx * dx + dy * dy;
            return r2 <= 0.09d ? r2 : 1000d;
        }

        /// <summary>
        /// A MultiStart whose bounds repair perturbs every value it is given, which makes any
        /// write-through of the repair observable as a values/fitness mismatch on the recorded best set.
        /// </summary>
        private sealed class PerturbingRepairMultiStart : MultiStart
        {
            public PerturbingRepairMultiStart(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds, LocalMethod method)
                : base(objectiveFunction, numberOfParameters, initialValues, lowerBounds, upperBounds, method) { }

            protected override double RepairParameter(double value, double lowerBound, double upperBound)
            {
                return base.RepairParameter(value + 0.25d, lowerBound, upperBound);
            }
        }

        /// <summary>
        /// Test that the bounds repair before the polish search cannot rewrite the recorded best set.
        /// </summary>
        /// <remarks>
        /// The polish path hands the live array inside <see cref="Optimizer.BestParameterSet"/> to the
        /// local-search entry point, whose bounds repair must not write through it, because the recorded
        /// fitness would keep describing the unrepaired point. A repair that changes its input is
        /// unreachable for a legally constructed solver, whose samples are already inside the box, so the
        /// repair is forced here through an override that perturbs every value it is given. The
        /// perturbation kicks the well-bottom best point onto the plateau, where the polish search cannot
        /// re-attain the recorded fitness, so a write-through cannot be healed by re-recording. The
        /// invariant is that the best set still satisfies fitness = f(values) exactly.
        /// </remarks>
        [TestMethod]
        public void Test_BoundsRepairDoesNotRewriteTheBestParameterSet()
        {
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new PerturbingRepairMultiStart(WellInPlateau, 2, new double[] { 0.5d, 0.5d }, lower, upper, LocalMethod.NelderMead)
            { MaxIterations = 25, ReportFailure = false, ComputeHessian = false };
            solver.Minimize();

            Assert.AreEqual(solver.BestParameterSet.Fitness, WellInPlateau(solver.BestParameterSet.Values), 0d,
                "The best parameter set's recorded fitness no longer matches its recorded values; the bounds repair wrote through the recorded array.");
        }

        /// <summary>
        /// A quadratic whose unconstrained minimum lies far outside the unit square used by the local method
        /// tests below, so that the constrained solution is the corner (1, 1) with value 722.
        /// </summary>
        /// <param name="x">The point to evaluate.</param>
        /// <returns>The value of the quadratic at <paramref name="x"/>.</returns>
        private static double CornerQuadratic(double[] x)
        {
            return (x[0] - 20d) * (x[0] - 20d) + (x[1] - 20d) * (x[1] - 20d);
        }

        /// <summary>
        /// Test that every local method this solver supports returns a point inside the declared bounds, and
        /// that the reported fitness is the objective function evaluated at that same point.
        /// </summary>
        /// <param name="method">The local search method to drive.</param>
        /// <param name="polish">Whether a final local search polishes the best member.</param>
        /// <remarks>
        /// <para>
        /// The local solvers are built over this solver's own objective delegate, so their probes are scored
        /// through <see cref="Optimizer.Evaluate"/> and any of them can be recorded as the incumbent. The
        /// objective here falls monotonically toward its unconstrained minimum outside the box, so any probe
        /// that leaves the box scores better than the constrained corner and would be reported.
        /// </para>
        /// <para>
        /// Both settings of <see cref="MultiStart.Polish"/> are covered because the two can fail in
        /// different ways: without polishing, an unguarded solver reports the infeasible probe itself;
        /// with polishing, the probe's point is repaired into the box while its fitness is kept,
        /// producing a feasible-looking point carrying a fitness from somewhere else.
        /// </para>
        /// <para>
        /// The fitness comparison is exact rather than approximate. <see cref="Optimizer.Evaluate"/> stores
        /// the point and the value it computed there in the same operation, so the two can only disagree if
        /// something later replaced one of them.
        /// </para>
        /// </remarks>
        [TestMethod]
        [DataRow(LocalMethod.BFGS, true)]
        [DataRow(LocalMethod.BFGS, false)]
        [DataRow(LocalMethod.NelderMead, true)]
        [DataRow(LocalMethod.NelderMead, false)]
        [DataRow(LocalMethod.Powell, true)]
        [DataRow(LocalMethod.Powell, false)]
        [DataRow(LocalMethod.LBFGSB, true)]
        [DataRow(LocalMethod.LBFGSB, false)]
        public void Test_SupportedLocalMethodReturnsAFeasibleAndConsistentSolution(LocalMethod method, bool polish)
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new MultiStart(CornerQuadratic, 2, initial, lower, upper, method)
            {
                ReportFailure = false,
                RecordTraces = false,
                // Every local start converges on the same corner, so a handful of starts covers the guarded
                // path as thoroughly as the default hundred and keeps the test cheap.
                MaxIterations = 10,
                Polish = polish
            };
            solver.Minimize();

            var solution = solver.BestParameterSet.Values;
            Assert.IsNotNull(solution);
            for (int i = 0; i < solution.Length; i++)
            {
                Assert.IsGreaterThanOrEqualTo(lower[i], solution[i], "Parameter " + i + " is below its lower bound.");
                Assert.IsLessThanOrEqualTo(upper[i], solution[i], "Parameter " + i + " is above its upper bound.");
            }

            // The constrained minimum is the corner nearest the unconstrained optimum.
            Assert.AreEqual(1d, solution[0], 1E-6);
            Assert.AreEqual(1d, solution[1], 1E-6);
            Assert.AreEqual(CornerQuadratic(solution), solver.BestParameterSet.Fitness, "The reported fitness must be the objective function evaluated at the reported point.");
        }

        /// <summary>
        /// Test that the local methods this solver does not build report that rather than returning a point.
        /// </summary>
        /// <param name="method">The unsupported local search method.</param>
        /// <remarks>
        /// <see cref="LocalMethod"/> also names the Adam and gradient descent algorithms, which this solver
        /// has never constructed. This pins the behaviour so that the coverage above is known to be complete
        /// for the enumeration as it stands, and so that adding support later is a deliberate change rather
        /// than a silent one.
        /// </remarks>
        [TestMethod]
        [DataRow(LocalMethod.ADAM)]
        [DataRow(LocalMethod.GradientDescent)]
        public void Test_UnsupportedLocalMethodIsReported(LocalMethod method)
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new MultiStart(CornerQuadratic, 2, initial, lower, upper, method) { RecordTraces = false, MaxIterations = 10 };

            Assert.ThrowsExactly<NotSupportedException>(() => solver.Minimize());
        }

        /// <summary>
        /// Test that maximization stores the negated objective at the reported point.
        /// </summary>
        /// <param name="method">The local search method to drive.</param>
        /// <remarks>
        /// <see cref="Optimizer.Evaluate"/> scales the objective by minus one while maximizing and stores
        /// that scaled value as the fitness. The local solvers are minimized over this solver's evaluation
        /// routine, which has already applied the scale, so the convention is the same for every method: the
        /// stored fitness is the negated objective at the stored point.
        /// </remarks>
        [TestMethod]
        [DataRow(LocalMethod.BFGS)]
        [DataRow(LocalMethod.NelderMead)]
        [DataRow(LocalMethod.Powell)]
        [DataRow(LocalMethod.LBFGSB)]
        public void Test_MaximizeStoresTheNegatedObjective(LocalMethod method)
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            double peak(double[] x) => 5d - (x[0] - 0.3d) * (x[0] - 0.3d) - (x[1] - 0.4d) * (x[1] - 0.4d);
            var solver = new MultiStart(peak, 2, initial, lower, upper, method) { ReportFailure = false, RecordTraces = false, MaxIterations = 10, ComputeHessian = false };
            solver.Maximize();

            var solution = solver.BestParameterSet.Values;
            Assert.IsNotNull(solution);
            // The location tolerance is loose because the point of this test is the sign convention on the
            // stored fitness, not the accuracy of any one local method. The simplex method in particular
            // stops on its own simplex size rather than on a gradient.
            Assert.AreEqual(0.3d, solution[0], 1E-3);
            Assert.AreEqual(0.4d, solution[1], 1E-3);
            Assert.AreEqual(-peak(solution), solver.BestParameterSet.Fitness, "While maximizing, the stored fitness is the negated objective at the stored point.");
        }

        /// <summary>
        /// Seeded runs with Nelder-Mead and Powell local searches are reproducible bit for bit.
        /// </summary>
        /// <remarks>
        /// The fingerprint pins the exact search path; see <see cref="OptimizerFingerprint"/>. A different
        /// fingerprint means different results.
        /// </remarks>
        [TestMethod]
        public void Test_SeededRunFingerprint()
        {
            var a = new MultiStart(OptimizerFingerprint.SixHumpCamel, 2, new[] { 1d, 1d }, new[] { -3d, -2d }, new[] { 3d, 2d }, LocalMethod.NelderMead) { MaxIterations = 10, ReportFailure = false, ComputeHessian = false };
            a.Minimize();
            Assert.AreEqual("x=BFB7019491425294,3FE6CE0AC7BDBD28|f=BFF0818CD58D6B55|it=10|fe=884|st=Success|trace=884:2A20AAD904B2F617", OptimizerFingerprint.Of(a));
            var b = new MultiStart(OptimizerFingerprint.SixHumpCamel, 2, new[] { 1d, 1d }, new[] { -3d, -2d }, new[] { 3d, 2d }, LocalMethod.Powell) { MaxIterations = 10, ReportFailure = false, ComputeHessian = false };
            b.Minimize();
            Assert.AreEqual("x=3FB6FFE2D4D130E2,BFE6CE14CD6483CF|f=BFF0818CD655CB15|it=10|fe=2347|st=Success|trace=2347:3357740C55648D6E", OptimizerFingerprint.Of(b));
        }

        /// <summary>
        /// A completed run reports the polish status and returns a Hessian, and the search
        /// stops at the evaluation budget.
        /// </summary>
        [TestMethod]
        public void Test_StatusHessianAndBudget()
        {
            var solver = new MultiStart(OptimizerFingerprint.SixHumpCamel, 2, new[] { 1d, 1d }, new[] { -3d, -2d }, new[] { 3d, 2d }) { MaxIterations = 10 };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsNotNull(solver.Hessian);
            Assert.AreEqual(-1.031628453489877, solver.BestParameterSet.Fitness, 1E-8);

            // Nelder-Mead local searches stop exactly at the budget. (With BFGS, the finite-difference
            // gradient retries a few evaluations after the budget is reached.)
            var budget = new MultiStart(OptimizerFingerprint.SixHumpCamel, 2, new[] { 1d, 1d }, new[] { -3d, -2d }, new[] { 3d, 2d }, LocalMethod.NelderMead) { MaxFunctionEvaluations = 200 };
            budget.Minimize();
            Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, budget.Status);
            Assert.AreEqual(200, budget.FunctionEvaluations);
        }

        /// <summary>
        /// L-BFGS-B local searches keep a bound-active coordinate exactly on its bound and are reproducible bit for bit.
        /// </summary>
        /// <remarks>
        /// The minimum of (x - 2)² + (y - 0.3)² on the unit square is (1, 0.3): x is held at its upper bound and y is
        /// interior. The fingerprint pins the seeded search path; see <see cref="OptimizerFingerprint"/>.
        /// </remarks>
        [TestMethod]
        public void Test_LBFGSBLocalSearchWithActiveBound()
        {
            double f(double[] p) => (p[0] - 2) * (p[0] - 2) + (p[1] - 0.3) * (p[1] - 0.3);
            var solver = new MultiStart(f, 2, new[] { 0.5, 0.5 }, new[] { 0d, 0d }, new[] { 1d, 1d }, LocalMethod.LBFGSB) { MaxIterations = 10, ReportFailure = false, ComputeHessian = false };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(1d, solver.BestParameterSet.Values[0]);
            Assert.AreEqual(0.3, solver.BestParameterSet.Values[1], 1E-6);
            Assert.AreEqual(f(solver.BestParameterSet.Values), solver.BestParameterSet.Fitness);
            Assert.AreEqual("x=3FF0000000000000,3FD3333333337FDE|f=3FF0000000000000|it=10|fe=211|st=Success|trace=211:D280F2AC3DA103F4", OptimizerFingerprint.Of(solver));
        }

        /// <summary>
        /// When fewer than 10 evaluations remain for a local search, the local solver still gets its minimum
        /// budget and this optimizer's evaluation cap stops the run, which ends with
        /// <see cref="OptimizationStatus.MaximumFunctionEvaluationsReached"/>, never with <see cref="OptimizationStatus.None"/>.
        /// The budgets are ones where the last local search starts with fewer than 10 evaluations left. The cap is
        /// exact whether or not <see cref="Optimizer.ReportFailure"/> is set.
        /// </summary>
        [TestMethod]
        public void Test_LocalSearchWithNearlySpentBudget()
        {
            static double Rosenbrock(double[] x) => 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
            foreach (bool reportFailure in new[] { true, false })
            {
                foreach (int budget in new[] { 190, 194, 198 })
                {
                    var solver = new MultiStart(Rosenbrock, 2, new[] { 0.1, 0.1 }, new[] { -5d, -5d }, new[] { 5d, 5d }, LocalMethod.NelderMead) { MaxFunctionEvaluations = budget, ReportFailure = reportFailure };
                    solver.Minimize();
                    Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, solver.Status, $"MaxFunctionEvaluations {budget}, ReportFailure {reportFailure}");
                    Assert.AreEqual(budget, solver.FunctionEvaluations, $"MaxFunctionEvaluations {budget}, ReportFailure {reportFailure}");
                }
            }
        }

        /// <summary>
        /// A Powell local search whose first line-search evaluation meets this optimizer's evaluation cap ends
        /// the run as a budget stop rather than failing with a missing line minimum.
        /// </summary>
        [TestMethod]
        public void Test_PowellLocalSearchAtEvaluationCap()
        {
            static double Rosenbrock(double[] x) => 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
            foreach (bool reportFailure in new[] { true, false })
            {
                var solver = new MultiStart(Rosenbrock, 2, new[] { 0.1, 0.1 }, new[] { -5d, -5d }, new[] { 5d, 5d }, LocalMethod.Powell) { MaxFunctionEvaluations = 272, ReportFailure = reportFailure };
                solver.Minimize();
                Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, solver.Status, $"ReportFailure {reportFailure}");
                Assert.IsLessThanOrEqualTo(272, solver.FunctionEvaluations, $"ReportFailure {reportFailure}");
            }
        }
    }
}
