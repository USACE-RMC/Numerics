using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using System;
using System.Linq;
using Utilities;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the CMA-ES optimization algorithm with IPOP restarts.
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
    public class Test_CMAES
    {
        /// <summary>
        /// The sphere function shifted to x* = (1, ..., 1).
        /// </summary>
        private static double Sphere(double[] x)
        {
            double s = 0;
            for (int i = 0; i < x.Length; i++) s += (x[i] - 1) * (x[i] - 1);
            return s;
        }

        /// <summary>
        /// The ellipsoid with condition number 1e6, f = sum c_i z_i^2, optionally in a rotated basis R.
        /// By default z = R (x - 1), so x* = (1, ..., 1); about the origin z = R x - 1, which rotates the
        /// axis-aligned problem about the box center.
        /// </summary>
        private static Func<double[], double> Ellipsoid(int n, bool rotate, bool aboutOrigin = false)
        {
            // R is a product of Givens rotations G(p, q, 0.5 + p q), so it is orthogonal and mixes every coordinate
            var R = new double[n, n];
            for (int i = 0; i < n; i++) R[i, i] = 1;
            for (int p = 0; p < n && rotate; p++)
                for (int q = p + 1; q < n; q++)
                {
                    double t = 0.5 + p * q, cs = Math.Cos(t), sn = Math.Sin(t);
                    for (int j = 0; j < n; j++)
                    {
                        double ap = R[p, j], aq = R[q, j];
                        R[p, j] = cs * ap - sn * aq;
                        R[q, j] = sn * ap + cs * aq;
                    }
                }
            var c = new double[n];
            for (int i = 0; i < n; i++) c[i] = Math.Pow(10, 6.0 * i / (n - 1));
            return x =>
            {
                double f = 0;
                for (int i = 0; i < n; i++)
                {
                    double z = 0;
                    if (aboutOrigin)
                    {
                        for (int j = 0; j < n; j++) z += R[i, j] * x[j];
                        z -= 1;
                    }
                    else
                    {
                        for (int j = 0; j < n; j++) z += R[i, j] * (x[j] - 1);
                    }
                    f += c[i] * z * z;
                }
                return f;
            };
        }

        /// <summary>
        /// The Rastrigin function shifted to x* = (1, ..., 1), so the box center is not the optimum.
        /// </summary>
        private static double ShiftedRastrigin(double[] x)
        {
            double f = 10 * x.Length;
            for (int i = 0; i < x.Length; i++)
            {
                double z = x[i] - 1;
                f += z * z - 10 * Math.Cos(2 * Math.PI * z);
            }
            return f;
        }

        /// <summary>
        /// Returns n copies of a value.
        /// </summary>
        private static double[] Fill(int n, double v) => Enumerable.Repeat(v, n).ToArray();

        /// <summary>
        /// A single run converges on the shifted 10-D sphere.
        /// </summary>
        [TestMethod]
        public void Test_Sphere()
        {
            var solver = new CMAES(Sphere, 10, Fill(10, -5), Fill(10, 5)) { Restarts = 0 };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-8);
            foreach (var v in solver.BestParameterSet.Values) Assert.AreEqual(1d, v, 1E-4);
        }

        /// <summary>
        /// A single run converges on the rotated 10-D ellipsoid with condition number 1e6.
        /// </summary>
        [TestMethod]
        public void Test_RotatedEllipsoid()
        {
            var solver = new CMAES(Ellipsoid(10, true), 10, Fill(10, -5), Fill(10, 5)) { Restarts = 0 };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-8);
            foreach (var v in solver.BestParameterSet.Values) Assert.AreEqual(1d, v, 1E-3);
        }

        /// <summary>
        /// The Rosenbrock function in 2, 5 and 10 dimensions, with two restarts.
        /// </summary>
        [TestMethod]
        public void Test_Rosenbrock()
        {
            foreach (int n in new[] { 2, 5, 10 })
            {
                var solver = new CMAES(TestFunctions.Rosenbrock, n, Fill(n, -5), Fill(n, 5)) { Restarts = 2 };
                solver.Minimize();
                Assert.AreEqual(OptimizationStatus.Success, solver.Status);
                Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-8);
                foreach (var v in solver.BestParameterSet.Values) Assert.AreEqual(1d, v, 1E-3);
            }
        }

        /// <summary>
        /// The first run ends in a local minimum of the shifted 5-D Rastrigin function; the IPOP restarts with
        /// larger populations find the global minimum.
        /// </summary>
        [TestMethod]
        public void Test_RastriginIPOP()
        {
            var single = new CMAES(ShiftedRastrigin, 5, Fill(5, -5.12), Fill(5, 5.12)) { Restarts = 0 };
            single.Minimize();
            Assert.IsGreaterThan(0.5, single.BestParameterSet.Fitness);

            var solver = new CMAES(ShiftedRastrigin, 5, Fill(5, -5.12), Fill(5, 5.12));
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsGreaterThan(0, solver.RestartsPerformed);
            Assert.HasCount(solver.RestartsPerformed + 1, solver.RunTerminations);
            // Each restart doubles the population size: 8, 16, 32, ...
            for (int r = 0; r < solver.RunTerminations.Count; r++)
                Assert.AreEqual((8 << r).ToString(), solver.RunTerminations[r].Split('|')[1]);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-8);
            foreach (var v in solver.BestParameterSet.Values) Assert.AreEqual(1d, v, 1E-4);
        }

        /// <summary>
        /// The six-hump camel function has two global minima f = -1.0316 at (0.0898, -0.7126) and (-0.0898, 0.7126).
        /// </summary>
        [TestMethod]
        public void Test_SixHumpCamel()
        {
            var solver = new CMAES(OptimizerFingerprint.SixHumpCamel, 2, new[] { -3d, -2d }, new[] { 3d, 2d }) { Restarts = 0 };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(-1.0316284534898774, solver.BestParameterSet.Fitness, 1E-8);
            double x = solver.BestParameterSet.Values[0], y = solver.BestParameterSet.Values[1];
            Assert.AreEqual(0.08984201368301331, Math.Abs(x), 1E-4);
            Assert.AreEqual(0.7126564032704135, Math.Abs(y), 1E-4);
            Assert.IsLessThan(0, x * y);
        }

        /// <summary>
        /// CMA-ES is invariant to rotations of the search space: the axis-aligned ellipsoid and the same ellipsoid
        /// rotated about the initial mean need a similar number of evaluations to reach f &lt; 1e-8.
        /// </summary>
        [TestMethod]
        public void Test_RotationInvariance()
        {
            int n = 10;
            double Mean(bool rotate)
            {
                double total = 0;
                for (int seed = 1; seed <= 10; seed++)
                {
                    var f = Ellipsoid(n, rotate, aboutOrigin: true);
                    int count = 0, hit = -1;
                    Func<double[], double> g = x => { count++; double v = f(x); if (hit < 0 && v < 1E-8) hit = count; return v; };
                    // A wide box and a small step size keep the samples far from the bounds, where the
                    // bound transformation is the identity and does not break the invariance.
                    var solver = new CMAES(g, n, Fill(n, -100), Fill(n, 100)) { Restarts = 0, InitialStepSize = 0.01, PRNGSeed = seed, ComputeHessian = false, RecordTraces = false };
                    solver.Minimize();
                    Assert.IsGreaterThan(0, hit);
                    total += hit;
                }
                return total / 10;
            }
            double rotated = Mean(true), aligned = Mean(false);
            Assert.AreEqual(1d, rotated / aligned, 0.15);
        }

        /// <summary>
        /// Runs with the same seed are identical, and different seeds give different paths.
        /// </summary>
        [TestMethod]
        public void Test_Determinism()
        {
            string Run(int seed)
            {
                var solver = new CMAES(OptimizerFingerprint.Rosenbrock, 3, Fill(3, -5), Fill(3, 5)) { PRNGSeed = seed, Restarts = 1, ComputeHessian = false };
                solver.Minimize();
                return OptimizerFingerprint.Of(solver);
            }
            Assert.AreEqual(Run(7), Run(7));
            Assert.AreNotEqual(Run(7), Run(8));
        }

        /// <summary>
        /// Every evaluated point is inside the bounds, also when the optimum lies outside the box.
        /// </summary>
        [TestMethod]
        public void Test_BoundsRespected()
        {
            var lower = new[] { -5d, -1d, 2d };
            var upper = new[] { 5d, 1d, 3d };
            bool inside = true;
            Func<double[], double> f = x =>
            {
                for (int i = 0; i < 3; i++) inside &= x[i] >= lower[i] && x[i] <= upper[i];
                return (x[0] - 10) * (x[0] - 10) + x[1] * x[1] + (x[2] - 2.5) * (x[2] - 2.5);
            };
            var solver = new CMAES(f, 3, lower, upper) { Restarts = 1, ComputeHessian = false };
            solver.Minimize();
            Assert.IsTrue(inside);
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(5d, solver.BestParameterSet.Values[0], 1E-4);
            Assert.AreEqual(0d, solver.BestParameterSet.Values[1], 1E-4);
            Assert.AreEqual(2.5, solver.BestParameterSet.Values[2], 1E-4);
        }

        /// <summary>
        /// Single runs that start near a bound converge on the sphere for every seed. The cap on the coordinate
        /// standard deviations keeps the step size from drifting beyond the box, where the periodic bound
        /// transformation would turn the search into a random search.
        /// </summary>
        [TestMethod]
        public void Test_StartNearBound()
        {
            for (int seed = 1; seed <= 20; seed++)
            {
                var solver = new CMAES(Sphere, 10, Fill(10, -4), Fill(10, -5), Fill(10, 5)) { Restarts = 0, PRNGSeed = seed, ComputeHessian = false, RecordTraces = false };
                solver.Minimize();
                Assert.AreEqual(OptimizationStatus.Success, solver.Status);
                Assert.IsLessThan(1E-8, solver.BestParameterSet.Fitness, $"seed {seed}");
            }
        }

        /// <summary>
        /// The bound transformation maps every real number into [0, 1], is the identity in the interior,
        /// has period 2.2, and is inverted by InverseTransform.
        /// </summary>
        [TestMethod]
        public void Test_BoundTransform()
        {
            for (double y = -10; y <= 10; y += 0.013)
            {
                double t = CMAES.Transform(y);
                Assert.IsTrue(t >= 0 && t <= 1);
                Assert.AreEqual(t, CMAES.Transform(y + 2.2), 1E-12);
            }
            Assert.AreEqual(0.3, CMAES.Transform(0.3), 0);
            Assert.AreEqual(0d, CMAES.Transform(-0.05), 1E-15);
            Assert.AreEqual(1d, CMAES.Transform(1.05), 1E-15);
            // Values from pycma 4.5.0: BoxConstraintsLinQuadTransformation([[0, 1]]) at -0.2, 0.01, 0.99, 1.3, 3.0
            Assert.AreEqual(0.1, CMAES.Transform(-0.2), 1E-14);
            Assert.AreEqual(0.018, CMAES.Transform(0.01), 1E-14);
            Assert.AreEqual(0.982, CMAES.Transform(0.99), 1E-14);
            Assert.AreEqual(0.8, CMAES.Transform(1.3), 1E-14);
            Assert.AreEqual(0.8, CMAES.Transform(3.0), 1E-12);
            for (double u = 0; u <= 1; u += 0.01)
                Assert.AreEqual(u, CMAES.Transform(CMAES.InverseTransform(u)), 1E-14);
        }

        /// <summary>
        /// The initial point overload starts the first run at the given point.
        /// </summary>
        [TestMethod]
        public void Test_InitialValues()
        {
            var solver = new CMAES(TestFunctions.Rosenbrock, 2, new[] { -1.2, 1d }, new[] { -5d, -5d }, new[] { 5d, 5d }) { Restarts = 0, InitialStepSize = 0.05 };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(1d, solver.BestParameterSet.Values[0], 1E-4);
            Assert.AreEqual(1d, solver.BestParameterSet.Values[1], 1E-4);
            Assert.AreEqual(-1.2, solver.InitialValues![0]);
        }

        /// <summary>
        /// Maximization negates the objective.
        /// </summary>
        [TestMethod]
        public void Test_Maximize()
        {
            var solver = new CMAES(x => -Sphere(x), 3, Fill(3, -5), Fill(3, 5)) { Restarts = 0 };
            solver.Maximize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            foreach (var v in solver.BestParameterSet.Values) Assert.AreEqual(1d, v, 1E-4);
        }

        /// <summary>
        /// After a converged run, reaching the evaluation limit during a restart ends with success.
        /// Without a converged run, the limit is reported as usual.
        /// </summary>
        [TestMethod]
        public void Test_BudgetEndsRestarts()
        {
            var solver = new CMAES(ShiftedRastrigin, 5, Fill(5, -5.12), Fill(5, 5.12)) { MaxFunctionEvaluations = 20000, ComputeHessian = false };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(20000, solver.FunctionEvaluations);

            var cut = new CMAES(Sphere, 10, Fill(10, -5), Fill(10, 5)) { MaxFunctionEvaluations = 100, ReportFailure = false };
            cut.Minimize();
            Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, cut.Status);
        }

        /// <summary>
        /// Invalid constructor arguments and settings are rejected.
        /// </summary>
        [TestMethod]
        public void Test_InvalidInputs()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new CMAES(Sphere, 2, new[] { 0d }, new[] { 1d, 1d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new CMAES(Sphere, 2, new[] { 0d, 1d }, new[] { 1d, 1d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new CMAES(Sphere, 2, new[] { 2d, 0.5 }, new[] { 0d, 0d }, new[] { 1d, 1d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new CMAES(Sphere, 2, new[] { 0.5 }, new[] { 0d, 0d }, new[] { 1d, 1d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new CMAES(Sphere, 2, new[] { 0d, 0d }, new[] { 1d, 1d }) { PopulationSize = 1 }.Minimize());
            Assert.Throws<ArgumentOutOfRangeException>(() => new CMAES(Sphere, 2, new[] { 0d, 0d }, new[] { 1d, 1d }) { InitialStepSize = 0 }.Minimize());
            Assert.Throws<ArgumentOutOfRangeException>(() => new CMAES(Sphere, 2, new[] { 0d, 0d }, new[] { 1d, 1d }) { Restarts = -1 }.Minimize());
        }

        /// <summary>
        /// The default parameters follow Table 1 of Hansen (2023): lambda = 4 + floor(3 ln n).
        /// </summary>
        [TestMethod]
        public void Test_DefaultParameters()
        {
            Assert.AreEqual(4, new CMAES(Sphere, 1, Fill(1, 0), Fill(1, 1)).PopulationSize);
            Assert.AreEqual(10, new CMAES(Sphere, 10, Fill(10, 0), Fill(10, 1)).PopulationSize);
            Assert.AreEqual(12, new CMAES(Sphere, 20, Fill(20, 0), Fill(20, 1)).PopulationSize);
            var s = new CMAES(Sphere, 3, Fill(3, 0), Fill(3, 1));
            Assert.AreEqual(9, s.Restarts);
            Assert.AreEqual(0.3, s.InitialStepSize);
            Assert.AreEqual(12345, s.PRNGSeed);
        }

        /// <summary>
        /// Pins a seeded run on the products-only Rosenbrock function; see <see cref="OptimizerFingerprint"/>.
        /// </summary>
        /// <remarks>
        /// The objective uses products only, but the algorithm itself calls Math.Exp (step size), Math.Log (weights and
        /// the inverse normal CDF) and Math.Atan2, Math.Cos and Math.Sin (Jacobi eigensolver). On .NET Framework 4.8.1
        /// these differ from .NET 10 in the last bits, and on Linux the C library rounds them differently again, so the
        /// bitwise fingerprint is pinned on .NET 10 on Windows only, where it was recorded, and every target framework and
        /// platform checks the run on a tolerance (see <see cref="CLibraryRounding"/>).
        /// </remarks>
        [TestMethod]
        public void Test_SeededRunFingerprint()
        {
            var solver = new CMAES(OptimizerFingerprint.Rosenbrock, 4, Fill(4, -5), Fill(4, 5)) { Restarts = 1, ComputeHessian = false };
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-10);
            foreach (var v in solver.BestParameterSet.Values) Assert.AreEqual(1d, v, 1E-5);
#if NET
            if (CLibraryRounding.IsReferencePlatform)
                Assert.AreEqual("x=3FF0000055E9F204,3FF000007A6222F0,3FF00000DB7BF5BC,3FF00001DDC86F14|f=3DA0295D6888784E|it=343|fe=3792|st=Success|trace=3792:7B020A14680A142D", OptimizerFingerprint.Of(solver));
#endif
        }

        /// <summary>
        /// With Polish on, the best point after the restarts is refined by Nelder-Mead and is no worse.
        /// </summary>
        [TestMethod]
        public void Test_Polish()
        {
            var lower = new[] { -5d, -5d, -5d };
            var upper = new[] { 5d, 5d, 5d };
            var plain = new CMAES(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { Restarts = 0, ComputeHessian = false };
            plain.Minimize();
            var polished = new CMAES(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { Restarts = 0, Polish = true };
            polished.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, polished.Status);
            Assert.IsLessThanOrEqualTo(plain.BestParameterSet.Fitness, polished.BestParameterSet.Fitness);
            Assert.IsGreaterThan(plain.FunctionEvaluations, polished.FunctionEvaluations);
            Assert.IsNotNull(polished.Hessian);
        }

        /// <summary>
        /// The default evaluation budget is 10,000 times the number of parameters; it bounds the restart
        /// sequence, and a run that converged before the budget ran out reports Success.
        /// </summary>
        [TestMethod]
        public void Test_DefaultEvaluationBudget()
        {
            var lower = new[] { -5d, -5d, -5d, -5d, -5d, -5d };
            var upper = new[] { 5d, 5d, 5d, 5d, 5d, 5d };
            var solver = new CMAES(OptimizerFingerprint.Rosenbrock, 6, lower, upper) { ComputeHessian = false };
            Assert.AreEqual(60000, solver.MaxFunctionEvaluations);
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsLessThanOrEqualTo(60000, solver.FunctionEvaluations);
            Assert.IsLessThan(1E-8, solver.BestParameterSet.Fitness);

            var larger = new CMAES(OptimizerFingerprint.Rosenbrock, 6, lower, upper) { MaxFunctionEvaluations = 1000000 };
            Assert.AreEqual(1000000, larger.MaxFunctionEvaluations);
        }
    }
}
