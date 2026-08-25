using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the MLSL optimization algorithm
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
    public class Test_MLSL
    {
        /// <summary>
        /// Test the MLSL algorithm with a simple 3-dimensional test function.
        /// </summary>
        [TestMethod]
        public void Test_FXYZ()
        {
            var initial = new double[] { 0.2d, 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new MLSL(TestFunctions.FXYZ, 3, initial, lower, upper);
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
        /// Test the MLSL algorithm with the De Jong Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            var initial = new double[] { 1.0d, -1.0d, 2.0d, -2.0d, 1.0d };
            var lower = new double[] { -5.12d, -5.12d, -5.12d, -5.12d, -5.12d };
            var upper = new double[] { 5.12d, 5.12d, 5.12d, 5.12d, 5.12d };
            var solver = new MLSL(TestFunctions.DeJong, 5, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Sum of Power functions in 3-D.
        /// </summary>
        [TestMethod]
        public void Test_SumOfPowerFunctions()
        {
            var initial = new double[] { 0.5d, -0.5d, 0.5d };
            var lower = new double[] { -1d, -1d, -1d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new MLSL(TestFunctions.SumOfPowerFunctions, 3, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Rosenbrock Function
        /// </summary>
        [TestMethod]
        public void Test_Rosenbrock()
        {
            var initial = new double[] { 0, 0, 0, 0, 0 };
            var lower = new double[] { -2.048, -2.048, -2.048, -2.048, -2.048 };
            var upper = new double[] { 2.048, 2.048, 2.048, 2.048, 2.048 };
            var solver = new MLSL(TestFunctions.Rosenbrock, 5, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Booth Function
        /// </summary>
        [TestMethod]
        public void Test_Booth()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new MLSL(TestFunctions.Booth, 2, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Matyas Function
        /// </summary>
        [TestMethod]
        public void Test_Matyas()
        {
            var initial = new double[] { 2d, -2d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new MLSL(TestFunctions.Matyas, 2, initial, lower, upper);
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
        /// Test the MLSL algorithm with the McCormick Function
        /// </summary>
        [TestMethod]
        public void Test_McCormick()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -1.5d, -3d };
            var upper = new double[] { 4d, 4d };
            var solver = new MLSL(TestFunctions.McCormick, 2, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Beale Function
        /// </summary>
        [TestMethod]
        public void Test_Beale()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -4.5d, -4.5d };
            var upper = new double[] { 4.5d, 4.5d };
            var solver = new MLSL(TestFunctions.Beale, 2, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Goldstein-Price Function
        /// </summary>
        [TestMethod]
        public void Test_GoldsteinPrice()
        {
            var initial = new double[] { -1d, 1d };
            var lower = new double[] { -2d, -2d };
            var upper = new double[] { 2d, 2d };
            var solver = new MLSL(TestFunctions.GoldsteinPrice, 2, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Rastrigin Function
        /// </summary>
        [TestMethod]
        public void Test_Rastrigin()
        {
            var initial = new double[] { 1, 1, 1, 1, 1 };
            var lower = new double[] { -5.12, -5.12, -5.12, -5.12, -5.12 };
            var upper = new double[] { 5.12, 5.12, 5.12, 5.12, 5.12 };
            // Need to run a lot of starts
            var solver = new MLSL(TestFunctions.Rastrigin, 5, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Ackley Function
        /// </summary>
        [TestMethod]
        public void Test_Ackley()
        {
            var initial = new double[] { 1d, 1d };
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5d, 5d };
            var solver = new MLSL(TestFunctions.Ackley, 2, initial, lower, upper);
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
        /// Test the MLSL algorithm with the three hump camel Function
        /// </summary>
        [TestMethod]
        public void Test_ThreeHumpCamel()
        {
            var initial = new double[] { 2d, -2d };
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5, 5d };
            var solver = new MLSL(TestFunctions.ThreeHumpCamel, 2, initial, lower, upper);
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
        /// Test the MLSL algorithm with the Eggholder Function
        /// </summary>
        [TestMethod]
        public void Test_Eggholder()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -512d, -512d };
            var upper = new double[] { 512d, 512d };
            // Need to run a lot of starts
            var solver = new MLSL(TestFunctions.Eggholder, 2, initial, lower, upper);
            solver.SampleSize = 200;
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -959.6407;
            Assert.AreEqual(F, trueF, 0.5);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 512d;
            var validY = 404.2319d;
            Assert.AreEqual(x, validX, 1E-1);
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
        /// Test the MLSL algorithm with the tp2 Function
        /// </summary>
        [TestMethod]
        public void Test_TP2()
        {
            var initial = new double[] { 2d, 2d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 2d, 2d };
            var solver = new MLSL(TestFunctions.tp2, 2, initial, lower, upper);
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
        /// Test that the sample reduction sort orders by fitness and keeps equally fit sample points in the order they were generated.
        /// </summary>
        /// <remarks>
        /// The objective function returns the same value everywhere except in a narrow strip of the
        /// search space, where it is strictly lower. Most sampled points therefore tie exactly, and a
        /// few are distinctly better. The reduced sample is truncated at the sort boundary, so under an
        /// unstable sort the ties would decide which points start local searches, and therefore which
        /// optimum is returned. The strictly better points are not the first points generated, so the
        /// test pins the primary fitness key as well as the tie-break and fails if the sort is removed.
        /// </remarks>
        [TestMethod]
        public void Test_TiedFitnessPreservesSampleOrder()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };

            // Record the evaluated parameter arrays by reference. A sample point stores the same array
            // instance that was handed to the objective function, so the first occurrence of an array
            // in this list is the position at which that sample point was added.
            var evaluated = new List<double[]>();
            var sync = new object();

            // The objective is constant except in a narrow strip near the upper bound of the first
            // parameter. The initial point is outside the strip, so it belongs to the tied group.
            var solver = new MLSL(x => { lock (sync) { evaluated.Add(x); } return x[0] > 0.9d ? 0d : 1d; }, 2, initial, lower, upper)
            {
                // A small reduction parameter keeps the number of local searches down.
                Gamma = 0.02,
                ReportFailure = false
            };
            solver.Minimize();

            // Each iteration appends exactly one generation of sample points.
            Assert.AreEqual(0, solver.SampledPoints.Count % solver.SampleSize);
            Assert.IsGreaterThan(16, solver.SampledPoints.Count);

            // Resolve the generation order of every sample point.
            var generation = new int[solver.SampledPoints.Count];
            for (int i = 0; i < solver.SampledPoints.Count; i++)
            {
                var values = solver.SampledPoints[i].ParameterSet.Values;
                generation[i] = evaluated.FindIndex(v => ReferenceEquals(v, values));
                Assert.IsGreaterThanOrEqualTo(0, generation[i], "Every sample point must have been evaluated.");
            }

            // The objective takes exactly two values, and both groups must be populated for the test to
            // pin the primary sort key as well as the tie-break.
            int better = solver.SampledPoints.Count(p => p.ParameterSet.Fitness == 0d);
            int tied = solver.SampledPoints.Count(p => p.ParameterSet.Fitness == 1d);
            Assert.AreEqual(solver.SampledPoints.Count, better + tied);
            Assert.IsGreaterThan(0, better, "At least one point must fall in the strictly better strip.");

            // The tied group is larger than the insertion-sort threshold of the framework sort, so an
            // unstable sort is free to permute it.
            Assert.IsGreaterThan(16, tied);

            // The strictly better points must sort ahead of the tied points.
            for (int i = 0; i < better; i++)
                Assert.AreEqual(0d, solver.SampledPoints[i].ParameterSet.Fitness, $"Sample point {i} should be in the better group.");
            for (int i = better; i < solver.SampledPoints.Count; i++)
                Assert.AreEqual(1d, solver.SampledPoints[i].ParameterSet.Fitness, $"Sample point {i} should be in the tied group.");

            // No point of the better group was generated first, so leaving the sample in generation
            // order, or sorting on anything other than fitness, fails here.
            Assert.IsGreaterThan(0, generation[0], "The best sample point must have moved to the front of the sample.");

            // Within each fitness group the sample must still be in generation order.
            for (int i = 1; i < better; i++)
                Assert.IsGreaterThan(generation[i - 1], generation[i], $"Better sample point {i} is out of generation order.");
            for (int i = better + 1; i < solver.SampledPoints.Count; i++)
                Assert.IsGreaterThan(generation[i - 1], generation[i], $"Tied sample point {i} is out of generation order.");
        }

        /// <summary>
        /// Test that the sampled point list keeps the same instance for the whole optimization run.
        /// </summary>
        /// <remarks>
        /// The sample reduction sort runs on every iteration, and <see cref="MLSL.SampledPoints"/> is a
        /// public property. Sorting into a new list instance would leave a caller that captured the list
        /// during a run holding a detached copy that stops growing, so the ordered result must be copied
        /// back into the existing list.
        /// </remarks>
        [TestMethod]
        public void Test_SampledPointsKeepSameListInstance()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };

            // Record every distinct list instance observed from inside the run. A run that replaces the
            // list on each iteration records one instance per iteration.
            var instances = new List<List<MLSL.SamplePoint>>();
            var sync = new object();
            MLSL solver = null;
            solver = new MLSL(x =>
            {
                lock (sync)
                {
                    var points = solver?.SampledPoints;
                    if (points != null && (instances.Count == 0 || !ReferenceEquals(instances[instances.Count - 1], points)))
                        instances.Add(points);
                }
                return Math.Pow(x[0] - 0.3d, 2d) + Math.Pow(x[1] - 0.7d, 2d);
            }, 2, initial, lower, upper)
            {
                ReportFailure = false
            };
            solver.Minimize();

            // More than one generation of sample points was drawn, so the sort ran at least once.
            Assert.IsGreaterThan(solver.SampleSize, solver.SampledPoints.Count);

            // The list instance never changed, and it is the one the solver still exposes.
            Assert.HasCount(1, instances);
            Assert.IsTrue(ReferenceEquals(instances[0], solver.SampledPoints), "The sampled point list instance must not be replaced during a run.");
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
            var solver = new MLSL(x => Math.Pow(x[0] - 20d, 2d) + Math.Pow(x[1] - 20d, 2d), 2, initial, lower, upper)
            {
                ReportFailure = false
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
        /// Test that the public <see cref="MLSL.InitialValues"/> property still equals what was passed to
        /// the constructor after a run, and that the caller's own array is left untouched.
        /// </summary>
        /// <remarks>
        /// On the first iteration, <see cref="MLSL.InitialValues"/> itself was stored by reference into a
        /// <see cref="ParameterSet"/> and handed straight to the local solver's in-place bounds repair, so
        /// either path could silently corrupt the public property's own array. This pins both the public
        /// <see cref="MLSL.InitialValues"/> array and the caller's own array against that regression.
        /// </remarks>
        [TestMethod]
        public void Test_InitialValuesAreNotMutatedByARun()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var callerSnapshot = (double[])initial.Clone();
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new MLSL(TestFunctions.Booth, 2, initial, lower, upper);
            solver.Minimize();

            CollectionAssert.AreEqual(callerSnapshot, initial, "The caller's own array must not be modified by a run.");
            CollectionAssert.AreEqual(callerSnapshot, solver.InitialValues, "InitialValues must still equal what was passed to the constructor after a run.");
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
        /// <param name="polish">Whether a final local search polishes the best sampled point.</param>
        /// <remarks>
        /// <para>
        /// The local solvers are built over this solver's own objective delegate, so their probes are scored
        /// through <see cref="Optimizer.Evaluate"/> and any of them can be recorded as the incumbent. The
        /// objective here falls monotonically toward its unconstrained minimum outside the box, so any probe
        /// that leaves the box scores better than the constrained corner and would be reported.
        /// </para>
        /// <para>
        /// Both settings of <see cref="MLSL.Polish"/> are covered because the two used to fail differently.
        /// Without polishing, the reported point was the infeasible probe itself. With polishing, the
        /// reported point was repaired back to the corner while the fitness recorded at the infeasible probe
        /// was kept, which produced a feasible-looking point carrying a fitness from somewhere else.
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
        public void Test_SupportedLocalMethodReturnsAFeasibleAndConsistentSolution(LocalMethod method, bool polish)
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new MLSL(CornerQuadratic, 2, initial, lower, upper, method)
            {
                ReportFailure = false,
                RecordTraces = false,
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
            var solver = new MLSL(CornerQuadratic, 2, initial, lower, upper, method) { RecordTraces = false, MaxIterations = 10 };

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
        public void Test_MaximizeStoresTheNegatedObjective(LocalMethod method)
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            double peak(double[] x) => 5d - (x[0] - 0.3d) * (x[0] - 0.3d) - (x[1] - 0.4d) * (x[1] - 0.4d);
            var solver = new MLSL(peak, 2, initial, lower, upper, method) { ReportFailure = false, RecordTraces = false, MaxIterations = 10, ComputeHessian = false };
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
    }
}
