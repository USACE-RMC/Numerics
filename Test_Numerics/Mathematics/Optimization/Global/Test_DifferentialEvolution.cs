using System.Collections.Generic;
using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;
using System.Linq;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the Differential Evolution optimization algorithm
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
    public class Test_DifferentialEvolution
    {
        /// <summary>
        /// Test the DE algorithm with a simple 3-dimensional test function.
        /// </summary>
        [TestMethod]
        public void Test_FXYZ()
        {
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new DifferentialEvolution(TestFunctions.FXYZ, 3, lower, upper);
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
        /// Test the DE algorithm with the De Jong Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            var lower = new double[] { -5.12d, -5.12d, -5.12d, -5.12d, -5.12d };
            var upper = new double[] { 5.12d, 5.12d, 5.12d, 5.12d, 5.12d };
            var solver = new DifferentialEvolution(TestFunctions.DeJong, 5, lower, upper);
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
        /// Test the DE algorithm with the Sum of Power functions in 3-D.
        /// </summary>
        [TestMethod]
        public void Test_SumOfPowerFunctions()
        {
            var lower = new double[] { -1d, -1d, -1d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new DifferentialEvolution(TestFunctions.SumOfPowerFunctions, 3, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0d, 0.0d, 0.0d };
            // This one is challenging to get great precision on the parameters.
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-2);
        }

        /// <summary>
        /// Test the DE algorithm with the Rosenbrock Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_Rosenbrock()
        {
            var lower = new double[] { -2.048, -2.048, -2.048, -2.048, -2.048 };
            var upper = new double[] { 2.048, 2.048, 2.048, 2.048, 2.048 };
            var solver = new DifferentialEvolution(TestFunctions.Rosenbrock, 5, lower, upper);
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
        /// Test the DE algorithm with the Booth Function
        /// </summary>
        [TestMethod]
        public void Test_Booth()
        {
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new DifferentialEvolution(TestFunctions.Booth, 2, lower, upper);
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
        /// Test the DE algorithm with the Matyas Function
        /// </summary>
        [TestMethod]
        public void Test_Matyas()
        {
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new DifferentialEvolution(TestFunctions.Matyas, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = 0.0d;
            Assert.AreEqual(x, validX, 1E-3);
            Assert.AreEqual(y, validY, 1E-3);
        }

        /// <summary>
        /// Test the DE algorithm with the McCormick Function
        /// </summary>
        [TestMethod]
        public void Test_McCormick()
        {
            var lower = new double[] { -1.5d, -3d };
            var upper = new double[] { 4d, 4d };
            var solver = new DifferentialEvolution(TestFunctions.McCormick, 2, lower, upper);
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
        /// Test the DE algorithm with the Beale Function
        /// </summary>
        [TestMethod]
        public void Test_Beale()
        {
            var lower = new double[] { -4.5d, -4.5d };
            var upper = new double[] { 4.5d, 4.5d };
            var solver = new DifferentialEvolution(TestFunctions.Beale, 2, lower, upper);
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
        /// Test the DE algorithm with the Goldstein-Price Function
        /// </summary>
        [TestMethod]
        public void Test_GoldsteinPrice()
        {
            var lower = new double[] { -2d, -2d };
            var upper = new double[] { 2d, 2d };
            var solver = new DifferentialEvolution(TestFunctions.GoldsteinPrice, 2, lower, upper);
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
        /// Test the DE algorithm with the Rastrigin Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_Rastrigin()
        {
            var lower = new double[] { -5.12, -5.12, -5.12, -5.12, -5.12 };
            var upper = new double[] { 5.12, 5.12, 5.12, 5.12, 5.12 };
            var solver = new DifferentialEvolution(TestFunctions.Rastrigin, 5, lower, upper);
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
        /// Test the DE algorithm with the Ackley Function
        /// </summary>
        [TestMethod]
        public void Test_Ackley()
        {
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5d, 5d };
            var solver = new DifferentialEvolution(TestFunctions.Ackley, 2, lower, upper);
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
        /// Test the DE algorithm with the three hump camel Function
        /// </summary>
        [TestMethod]
        public void Test_ThreeHumpCamel()
        {
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5, 5d };
            var solver = new DifferentialEvolution(TestFunctions.ThreeHumpCamel, 2, lower, upper);
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
        /// Test the DE algorithm with the Eggholder Function
        /// </summary>
        [TestMethod]
        public void Test_Eggholder()
        {
            var lower = new double[] { -512d, -512d };
            var upper = new double[] { 512d, 512d };
            var solver = new DifferentialEvolution(TestFunctions.Eggholder, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -959.6407;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 512d;
            var validY = 404.2319d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-3);
        }

        /// <summary>
        /// Test the DE algorithm with the tp2 Function
        /// </summary>
        [TestMethod]
        public void Test_TP2()
        {
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 2d, 2d };
            var solver = new DifferentialEvolution(TestFunctions.tp2, 2, lower, upper);
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
        /// Verifies infeasible trial coordinates are repaired halfway toward the target instead of clamped to a boundary.
        /// </summary>
        /// <param name="seed">The deterministic pseudo-random seed used for the differential-evolution run.</param>
        /// <param name="expectedFirstTrial">The independently derived midpoint repair for the first trial vector.</param>
        [TestMethod]
        [DataRow(1, 0.8196490542613901d)]
        [DataRow(2, 0.7595016034028959d)]
        [DataRow(3, 0.7520687940414064d)]
        [DataRow(4, 0.48277057497762144d)]
        [DataRow(5, 0.8296836465888191d)]
        [DataRow(12345, 0.1893519861914683d)]
        public void Test_InfeasibleTrialCoordinates_AreRepairedHalfwayToTarget(int seed, double expectedFirstTrial)
        {
            var evaluatedCoordinates = new List<double>();
            double Objective(double[] values)
            {
                evaluatedCoordinates.Add(values[0]);
                double difference = values[0] - 0.5d;
                return difference * difference;
            }

            var solver = new DifferentialEvolution(Objective, 1, new double[] { 0d }, new double[] { 1d })
            {
                PopulationSize = 4,
                PRNGSeed = seed,
                Mutation = 2d,
                DitherRate = 0d,
                CrossoverProbability = 1d,
                MaxIterations = 11,
                ComputeHessian = false,
                ReportFailure = false
            };

            solver.Minimize();

            Assert.IsGreaterThan(solver.PopulationSize, evaluatedCoordinates.Count);
            Assert.AreEqual(expectedFirstTrial, evaluatedCoordinates[solver.PopulationSize], 1E-15);
            foreach (double coordinate in evaluatedCoordinates.Skip(solver.PopulationSize))
            {
                Assert.IsTrue(
                    coordinate > 0d && coordinate < 1d,
                    $"Expected an interior repaired coordinate, but evaluated {coordinate:R}.");
            }
        }

        /// <summary>
        /// The one-pass donor selection draws the same keys and returns the same three donors as the stable
        /// random-key shuffle it replaces, leaving the generator in the same state.
        /// </summary>
        [TestMethod]
        public void Test_SelectDonors_MatchesStableShuffle()
        {
            foreach (int populationSize in new[] { 4, 5, 10, 37 })
            {
                var linq = new MersenneTwister(91);
                var fast = new MersenneTwister(91);
                for (int trial = 0; trial < 200; trial++)
                {
                    int target = trial % populationSize;
                    var expected = Enumerable.Range(0, populationSize).Where(x => x != target).OrderBy(_ => linq.NextDouble()).Take(3).ToArray();
                    DifferentialEvolution.SelectDonors(fast, populationSize, target, out int r0, out int r1, out int r2);
                    CollectionAssert.AreEqual(expected, new[] { r0, r1, r2 });
                }
                Assert.AreEqual(linq.NextDouble(), fast.NextDouble(), 0d);
            }
        }

        /// <summary>
        /// Seeded runs are reproducible bit for bit: the default, a second seed, a four-member population under maximization, and an evaluation budget.
        /// </summary>
        /// <remarks>
        /// The fingerprint pins the exact search path; see <see cref="OptimizerFingerprint"/>. A different
        /// fingerprint means different results.
        /// </remarks>
        [TestMethod]
        public void Test_SeededRunFingerprint()
        {
            var lower = new[] { -5d, -5d, -5d };
            var upper = new[] { 5d, 5d, 5d };
            var a = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { ComputeHessian = false };
            a.Minimize();
            Assert.AreEqual("x=3FEFFFF4FD66B384,3FEFFFE80FC16821,3FEFFFC56F5F7E6A|f=3E28631D229539D4|it=135|fe=4080|st=Success|trace=4080:5647497C593DC47A", OptimizerFingerprint.Of(a));
            var b = new DifferentialEvolution(OptimizerFingerprint.SixHumpCamel, 2, new[] { -3d, -2d }, new[] { 3d, 2d }) { PRNGSeed = 7, ComputeHessian = false };
            b.Minimize();
            Assert.AreEqual("x=BFB6FFEBC551E82B,3FE6CE0357BCBCE5|f=BFF0818CD62E6A4C|it=43|fe=880|st=Success|trace=880:98F35F8755E70BE0", OptimizerFingerprint.Of(b));
            var c = new DifferentialEvolution(x => -OptimizerFingerprint.Rosenbrock(x), 3, lower, upper) { PopulationSize = 4, MaxIterations = 200, ReportFailure = false, ComputeHessian = false };
            c.Maximize();
            Assert.AreEqual("x=3FB74776CF1C104C,BFD258DC0C2ED8EF,3FA88959C0DA1627|f=402698C91529DB9D|it=106|fe=428|st=Success|trace=428:692D3206F9C961FA", OptimizerFingerprint.Of(c));
            var d = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { MaxFunctionEvaluations = 500, ReportFailure = false };
            d.Minimize();
            Assert.AreEqual("x=3FF63FBEF6B5E002,3FFE5AE60D002B70,400DAEB7D588CECD|f=4002938BDDDF2B7B|it=16|fe=500|st=MaximumFunctionEvaluationsReached|trace=500:AC198443EE0908A4", OptimizerFingerprint.Of(d));
        }

        /// <summary>
        /// Polish is off by default; when on, the best point is refined by the
        /// chosen local method and the evaluations are counted.
        /// </summary>
        [TestMethod]
        public void Test_Polish()
        {
            var lower = new[] { -5d, -5d, -5d };
            var upper = new[] { 5d, 5d, 5d };
            var plain = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { ComputeHessian = false };
            plain.Minimize();
            foreach (var method in new[] { LocalMethod.NelderMead, LocalMethod.BFGS, LocalMethod.Powell })
            {
                var polished = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { Polish = true, PolishMethod = method };
                polished.Minimize();
                Assert.AreEqual(OptimizationStatus.Success, polished.Status);
                Assert.IsLessThanOrEqualTo(plain.BestParameterSet.Fitness, polished.BestParameterSet.Fitness);
                Assert.IsGreaterThan(plain.FunctionEvaluations, polished.FunctionEvaluations);
                Assert.IsNotNull(polished.Hessian);
            }
            var invalid = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { Polish = true, PolishMethod = LocalMethod.ADAM };
            Assert.ThrowsExactly<NotSupportedException>(() => invalid.Minimize());
        }

        /// <summary>
        /// Integer parameters are rounded before every evaluation, so a pure-integer and a mixed problem reach
        /// the optimum found by enumerating every integer point (the oracle).
        /// </summary>
        [TestMethod]
        public void Test_IntegerParameters()
        {
            // f = (x - 2.4)^2 + (y + 1.6)^2 + (x - y - 3.7)^2 / 4 over integers in [-5, 5]
            Func<double[], double> f = v => (v[0] - 2.4) * (v[0] - 2.4) + (v[1] + 1.6) * (v[1] + 1.6) + (v[0] - v[1] - 3.7) * (v[0] - v[1] - 3.7) / 4;
            double best = double.PositiveInfinity; double bx = 0, by = 0;
            for (int x = -5; x <= 5; x++)
                for (int y = -5; y <= 5; y++)
                    if (f(new double[] { x, y }) < best) { best = f(new double[] { x, y }); bx = x; by = y; }

            var solver = new DifferentialEvolution(f, 2, new[] { -5.3, -5d }, new[] { 5d, 5.7 }) { IntegerParameters = new[] { true, true }, ComputeHessian = false };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(bx, solver.BestParameterSet.Values[0], 0d);
            Assert.AreEqual(by, solver.BestParameterSet.Values[1], 0d);
            Assert.AreEqual(best, solver.BestParameterSet.Fitness, 0d);
            foreach (var p in solver.ParameterSetTrace)
                foreach (var v in p.Values) Assert.AreEqual(Math.Round(v), v, 0d);

            // Mixed: x integer, y continuous; for each integer x the best y is (x - 3.7 - 4 * 1.6) / 5 in closed form.
            double mixedBest = double.PositiveInfinity;
            for (int x = -5; x <= 5; x++)
            {
                double y = (x - 3.7 - 6.4) / 5d;
                mixedBest = Math.Min(mixedBest, f(new double[] { x, y }));
            }
            var mixed = new DifferentialEvolution(f, 2, new[] { -5d, -5d }, new[] { 5d, 5d }) { IntegerParameters = new[] { true, false }, ComputeHessian = false, Polish = true };
            mixed.Minimize();
            Assert.AreEqual(Math.Round(mixed.BestParameterSet.Values[0]), mixed.BestParameterSet.Values[0], 0d);
            Assert.AreEqual(mixedBest, mixed.BestParameterSet.Fitness, 1E-6);

            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new DifferentialEvolution(f, 2, new[] { 0.2, 0d }, new[] { 0.8, 1d }) { IntegerParameters = new[] { true, false } }.Minimize());
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new DifferentialEvolution(f, 2, new[] { 0d, 0d }, new[] { 1d, 1d }) { IntegerParameters = new[] { true } }.Minimize());
        }

        /// <summary>
        /// Initial members replace the first Latin hypercube members: a run seeded with the optimum keeps it,
        /// and invalid rows are rejected.
        /// </summary>
        [TestMethod]
        public void Test_InitialPopulation()
        {
            var lower = new[] { -5d, -5d, -5d };
            var upper = new[] { 5d, 5d, 5d };
            var plain = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { ComputeHessian = false };
            plain.Minimize();
            var warm = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { ComputeHessian = false, InitialPopulation = new List<double[]> { new[] { 1d, 1d, 1d } } };
            warm.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, warm.Status);
            Assert.AreEqual(0d, warm.BestParameterSet.Fitness, 0d);
            Assert.AreEqual(0d, warm.ParameterSetTrace[0].Fitness, 0d);
            Assert.IsLessThanOrEqualTo(plain.BestParameterSet.Fitness, warm.BestParameterSet.Fitness);

            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { InitialPopulation = new List<double[]> { new[] { 9d, 0d, 0d } } }.Minimize());
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { InitialPopulation = new List<double[]> { new[] { 0d, 0d } } }.Minimize());
        }

        /// <summary>
        /// Parallel evaluation gives the same result for any number of threads, converges, stops exactly at
        /// the evaluation budget, and surfaces an objective exception unwrapped.
        /// </summary>
        [TestMethod]
        public void Test_UseParallel()
        {
            var lower = new[] { -5d, -5d, -5d, -5d };
            var upper = new[] { 5d, 5d, 5d, 5d };
            string first = null;
            foreach (int threads in new[] { 1, 2, -1 })
            {
                var solver = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 4, lower, upper) { UseParallel = true, MaxDegreeOfParallelism = threads, ComputeHessian = false };
                solver.Minimize();
                string fp = OptimizerFingerprint.Of(solver);
                if (first == null) first = fp; else Assert.AreEqual(first, fp);
                Assert.IsLessThan(1E-2, solver.BestParameterSet.Fitness);
            }

            var budget = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 4, lower, upper) { UseParallel = true, MaxFunctionEvaluations = 500 };
            budget.Minimize();
            Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, budget.Status);
            Assert.AreEqual(500, budget.FunctionEvaluations);

            var failing = new DifferentialEvolution(x => x[0] > 4.0 ? throw new InvalidOperationException("bad point") : OptimizerFingerprint.Rosenbrock(x), 4, lower, upper) { UseParallel = true };
            Assert.ThrowsExactly<InvalidOperationException>(() => failing.Minimize());
        }
    }
}
