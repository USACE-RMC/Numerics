using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the Adaptive Movement (Adam) optimization algorithm. The objective function must be differentiable and convex.
    /// </summary>
    [TestClass]
    public class Test_Adam
    {
        /// <summary>
        /// Test the ADAM algorithm with a simple 3-dimensional test function.
        /// </summary>
        [TestMethod]
        public void Test_FXYZ()
        {
            var initial = new double[] { 0.2d, 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new ADAM(TestFunctions.FXYZ, 3, initial, lower, upper);
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
        /// Test the ADAM algorithm with the De Jong Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            var initial = new double[] { 1.0d, -1.0d, 2.0d, -2.0d, 1.0d };
            var lower = new double[] { -5.12d, -5.12d, -5.12d, -5.12d, -5.12d };
            var upper = new double[] { 5.12d, 5.12d, 5.12d, 5.12d, 5.12d };
            var solver = new ADAM(TestFunctions.DeJong, 5, initial, lower, upper);
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
        /// Test the ADAM algorithm with the Sum of Power functions in 2-D.
        /// </summary>
        [TestMethod]
        public void Test_SumOfPowerFunctions()
        {
            var initial = new double[] { 0.5d, -0.5d };
            var lower = new double[] { -1d, -1d };
            var upper = new double[] { 1d, 1d };
            var solver = new ADAM(TestFunctions.SumOfPowerFunctions, 2, initial, lower, upper);
            solver.MaxIterations = 100000; // Requires a lot of iterations
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0d, 0.0d };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the ADAM algorithm with the Rosenbrock Function in 2-D.
        /// </summary>
        [TestMethod]
        public void Test_Rosenbrock()
        {
            var initial = new double[] { 0, 0 };
            var lower = new double[] { -2.048, -2.048 };
            var upper = new double[] { 2.048, 2.048 };
            var solver = new ADAM(TestFunctions.Rosenbrock, 2, initial, lower, upper);
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
        /// Test the ADAM algorithm with the Booth Function
        /// </summary>
        [TestMethod]
        public void Test_Booth()
        {
            var initial = new double[] { 0.0d, 0.0d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new ADAM(TestFunctions.Booth, 2, initial, lower, upper);
            solver.MaxIterations = 100000; // Requires a lot of iterations
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
        /// Test the ADAM algorithm with the Matyas Function
        /// </summary>
        [TestMethod]
        public void Test_Matyas()
        {
            var initial = new double[] { 1.0d, -1.0d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new ADAM(TestFunctions.Matyas, 2, initial, lower, upper);
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
        /// Test the ADAM algorithm with the McCormick Function
        /// </summary>
        [TestMethod]
        public void Test_McCormick()
        {
            var initial = new double[] { 0.0d, 0.0d };
            var lower = new double[] { -1.5d, -3d };
            var upper = new double[] { 4d, 4d };
            var solver = new ADAM(TestFunctions.McCormick, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -1.9133;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = -0.54719d;
            var validY = -1.54719d;
            // Level of precision in parameters is not great
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
        }

        /// <summary>
        /// The finite-difference gradient must not evaluate the objective outside the declared bounds.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The unconstrained minimum of this objective is at (2, -2), so the solver parks on the corner
        /// (1, 0) and every subsequent gradient is taken at a point on both bounds. A perturbation that is
        /// not clamped evaluates the objective where it may be undefined, and the resulting point can
        /// displace the incumbent.
        /// </para>
        /// <para>
        /// The end-of-run Hessian pass in <see cref="Optimizer.Minimize"/> is a separate finite-difference
        /// path that is not given the bounds, so it is switched off here to keep this test scoped to the
        /// gradient.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_GradientDoesNotProbeOutsideBounds()
        {
            var initial = new double[] { 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var violations = new List<string>();

            var solver = new ADAM(x =>
            {
                for (int i = 0; i < x.Length; i++)
                    if (x[i] < lower[i] || x[i] > upper[i]) violations.Add("p" + i + "=" + x[i]);
                return Math.Pow(x[0] - 2d, 2d) + Math.Pow(x[1] + 2d, 2d);
            }, 2, initial, lower, upper) { ComputeHessian = false };
            solver.Minimize();

            Assert.IsEmpty(violations, "The objective was evaluated outside the declared bounds: " + string.Join(", ", violations));
            Assert.AreEqual(1d, solver.BestParameterSet.Values[0], 1E-6);
            Assert.AreEqual(0d, solver.BestParameterSet.Values[1], 1E-6);
        }

    }
}
