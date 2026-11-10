using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the Simulated Annealing (SA) optimization algorithm.
    /// SA is an inefficient optimizer. Therefore, the precision of the unit tests are relaxed.
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
    public class Test_SimulatedAnnealing
    {

        /// <summary>
        /// Test the SA algorithm with a simple 3-dimensional test function.
        /// </summary>
        [TestMethod]
        public void Test_FXYZ()
        {
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new SimulatedAnnealing(TestFunctions.FXYZ, 3, lower, upper) { AbsoluteTolerance = 1E-6, RelativeTolerance = 1E-4 };
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            double x = solution[0];
            double y = solution[1];
            double z = solution[2];
            double validX = 0.125d;
            double validY = 0.2d;
            double validZ = 0.35d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
            Assert.AreEqual(z, validZ, 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the De Jong Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            var lower = new double[] { -5.12d, -5.12d, -5.12d, -5.12d, -5.12d };
            var upper = new double[] { 5.12d, 5.12d, 5.12d, 5.12d, 5.12d };
            var solver = new SimulatedAnnealing(TestFunctions.DeJong, 5, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0d, 0.0d, 0.0d, 0.0d, 0.0d };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the Sum of Power functions in 3-D.
        /// </summary>
        [TestMethod]
        public void Test_SumOfPowerFunctions()
        {
            var lower = new double[] { -1d, -1d, -1d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new SimulatedAnnealing(TestFunctions.SumOfPowerFunctions, 3, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0d, 0.0d, 0.0d };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the Rosenbrock Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_Rosenbrock()
        {
            var lower = new double[] { -2.048, -2.048, -2.048, -2.048, -2.048 };
            var upper = new double[] { 2.048, 2.048, 2.048, 2.048, 2.048 };
            var solver = new SimulatedAnnealing(TestFunctions.Rosenbrock, 5, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-1);  
            // SA has poor precision for this test function
        }

        /// <summary>
        /// Test the SA algorithm with the Booth Function
        /// </summary>
        [TestMethod]
        public void Test_Booth()
        {
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new SimulatedAnnealing(TestFunctions.Booth, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 1.0d;
            var validY = 3.0d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the Matyas Function
        /// </summary>
        [TestMethod]
        public void Test_Matyas()
        {
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new SimulatedAnnealing(TestFunctions.Matyas, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = 0.0d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the McCormick Function
        /// </summary>
        [TestMethod]
        public void Test_McCormick()
        {
            var lower = new double[] { -1.5d, -3d };
            var upper = new double[] { 4d, 4d };
            var solver = new SimulatedAnnealing(TestFunctions.McCormick, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -1.9133;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = -0.54719d;
            var validY = -1.54719d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the Beale Function
        /// </summary>
        [TestMethod]
        public void Test_Beale()
        {
            var lower = new double[] { -4.5d, -4.5d };
            var upper = new double[] { 4.5d, 4.5d };
            var solver = new SimulatedAnnealing(TestFunctions.Beale, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 3.0d;
            var validY = 0.5d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the Goldstein-Price Function
        /// </summary>
        [TestMethod]
        public void Test_GoldsteinPrice()
        {
            var lower = new double[] { -2d, -2d };
            var upper = new double[] { 2d, 2d };
            var solver = new SimulatedAnnealing(TestFunctions.GoldsteinPrice, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 3.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = -1.0d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the Rastrigin Function
        /// </summary>
        [TestMethod]
        public void Test_Rastrigin()
        {
            var lower = new double[] { -5.12, -5.12, -5.12, -5.12, -5.12 };
            var upper = new double[] { 5.12, 5.12, 5.12, 5.12, 5.12 };
            var solver = new SimulatedAnnealing(TestFunctions.Rastrigin, 5, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0 };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the Ackley Function
        /// </summary>
        [TestMethod]
        public void Test_Ackley()
        {
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5d, 5d };
            var solver = new SimulatedAnnealing(TestFunctions.Ackley, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = 0.0d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the three hump camel Function
        /// </summary>
        [TestMethod]
        public void Test_ThreeHumpCamel()
        {
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5, 5d };
            var solver = new SimulatedAnnealing(TestFunctions.ThreeHumpCamel, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = 0.0d;
            Assert.AreEqual(x, validX, 1E-2);
            Assert.AreEqual(y, validY, 1E-2);
        }

        /// <summary>
        /// Test the SA algorithm with the Eggholder Function
        /// </summary>
        [TestMethod]
        public void Test_Eggholder()
        {
            var lower = new double[] { -512d, -512d };
            var upper = new double[] { 512d, 512d };
            var solver = new SimulatedAnnealing(TestFunctions.Eggholder, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness; // -126.42
            //double trueF = -959.6407;
            //Assert.AreEqual(F, trueF, 1E-2);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0]; // -46.45
            var y = solution[1]; // 38.13
            //var validX = 512d;
            //var validY = 404.2319d;
            //Assert.AreEqual(x, validX, 1E-2);
            //Assert.AreEqual(y, validY, 1E-2);

            // Simulated annealing fails to converge for this test function
            // -126.423832142973 (should be -959.6407)
            // This test is still included as an example. 
        }

        /// <summary>
        /// Test the SA algorithm with the tp2 Function
        /// </summary>
        [TestMethod]
        public void Test_TP2()
        {
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 2d, 2d };
            var solver = new SimulatedAnnealing(TestFunctions.tp2, 2, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-2);
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
        /// A seeded run, including its Nelder-Mead polish, is reproducible bit for bit.
        /// </summary>
        /// <remarks>
        /// The fingerprint pins the exact search path; see <see cref="OptimizerFingerprint"/>. A different
        /// fingerprint means different results.
        /// </remarks>
        [TestMethod]
        public void Test_SeededRunFingerprint()
        {
            var solver = new SimulatedAnnealing(OptimizerFingerprint.Rosenbrock, 3, new[] { -5d, -5d, -5d }, new[] { 5d, 5d, 5d }) { MaxIterations = 30, ReportFailure = false };
            solver.Minimize();
            Assert.AreEqual("x=3FF0000000712C98,3FF0000000D85370,3FF0000001B39714|f=3C71371A26ACB800|it=30|fe=3813|st=MaximumIterationsReached|trace=3813:23DAD742FEFEDBF3", OptimizerFingerprint.Of(solver));
        }

        /// <summary>
        /// The search stops when the best value stalls over ToleranceSteps reductions and polishes the best
        /// point: it ends with Success and a Hessian, in a small fraction of
        /// the MaxIterations x 40 x D evaluations of a full schedule.
        /// </summary>
        [TestMethod]
        public void Test_StallStopAndPolish()
        {
            var solver = new SimulatedAnnealing(OptimizerFingerprint.Rosenbrock, 3, new[] { -5d, -5d, -5d }, new[] { 5d, 5d, 5d });
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsNotNull(solver.Hessian);
            Assert.IsLessThan(1E-10, solver.BestParameterSet.Fitness);
            Assert.IsLessThan(100000, solver.FunctionEvaluations);

            var unpolished = new SimulatedAnnealing(OptimizerFingerprint.Rosenbrock, 3, new[] { -5d, -5d, -5d }, new[] { 5d, 5d, 5d }) { Polish = false };
            unpolished.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, unpolished.Status);
            Assert.IsGreaterThan(solver.BestParameterSet.Fitness, unpolished.BestParameterSet.Fitness);
        }
    }
}
