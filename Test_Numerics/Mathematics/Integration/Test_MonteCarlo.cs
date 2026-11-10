using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the Monte Carlo method of integration
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_MonteCarlo
    {
        /// <summary>
        /// Test the Monte Carlo method with the Pi function
        /// </summary>
        [TestMethod()]
        public void Test_PI()
        {
            var mc = new MonteCarloIntegration(Integrands.PI, 2, new double[] { -1, -1 }, new double[] { 1, 1 });
            mc.Random = new MersenneTwister(12345);
            mc.MinIterations = 1000;
            mc.MaxIterations = 100000000;
            mc.AbsoluteTolerance = 1E-8;
            mc.RelativeTolerance = 1E-4;
            mc.Integrate();
            var result = mc.Result;
            double trueResult = 3.14;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Monte Carlo method with the GSL function
        /// </summary>
        [TestMethod()]
        public void Test_GSL()
        {
            var mc = new MonteCarloIntegration(Integrands.GSL, 3, new double[] { 0, 0, 0 }, new double[] { Math.PI, Math.PI, Math.PI });
            mc.Random = new MersenneTwister(12345);
            mc.MinIterations = 1000;
            mc.MaxIterations = 100000000;
            mc.MaxFunctionEvaluations = 100000000;
            mc.AbsoluteTolerance = 1E-8;
            mc.RelativeTolerance = 1E-4;
            mc.Integrate();
            var result = mc.Result;
            double trueResult = 1.3932039296856768591842462603255;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Monte Carlo method with the sum of 3 normal distributions
        /// </summary>
        [TestMethod()]
        public void Test_SumOfThreeNormals()
        {
            var min = new double[3];
            var max = new double[3];
            for (int i = 0; i < 3; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            var mc = new MonteCarloIntegration(Integrands.SumOfNormals, 3, min, max);
            mc.Random = new MersenneTwister(12345);
            mc.MinIterations = 1000;
            mc.MaxIterations = 100000000;
            mc.AbsoluteTolerance = 1E-8;
            mc.RelativeTolerance = 1E-4;
            mc.Integrate();
            var result = mc.Result;
            double trueResult = 57;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Monte Carlo method with the sum of 5 normal distributions
        /// </summary>
        [TestMethod()]
        public void Test_SumOfFiveNormals()
        {
            var min = new double[5];
            var max = new double[5];
            for (int i = 0; i < 5; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            var mc = new MonteCarloIntegration(Integrands.SumOfNormals, 5, min, max);
            mc.Random = new MersenneTwister(12345);
            mc.MinIterations = 1000;
            mc.MaxIterations = 100000000;
            mc.AbsoluteTolerance = 1E-8;
            mc.RelativeTolerance = 1E-4;
            mc.Integrate();
            var result = mc.Result;
            double trueResult = 224;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Monte Carlo method with the sum of 20 normal distributions
        /// </summary>
        [TestMethod()]
        public void Test_SumOfTwentyNormals()
        {
            var min = new double[20];
            var max = new double[20];
            for (int i = 0; i < 20; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            var mc = new MonteCarloIntegration(Integrands.SumOfNormals, 20, min, max);
            mc.Random = new MersenneTwister(12345);
            mc.MinIterations = 1000;
            mc.MaxIterations = 100000000;
            mc.AbsoluteTolerance = 1E-8;
            mc.RelativeTolerance = 1E-4;
            mc.Integrate();
            var result = mc.Result;
            double trueResult = 837;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);

        }

        /// <summary>
        /// MaxFunctionEvaluations caps the sample count, and stopping on a limit without convergence is
        /// not reported as Success. The capped estimate equals a hand count of the first 100 points
        /// drawn from the same generator.
        /// </summary>
        [TestMethod]
        public void Test_Budget_And_Status()
        {
            var a = new double[] { -1, -1 };
            var b = new double[] { 1, 1 };
            var capped = new MonteCarloIntegration(Integrands.PI, 2, a, b) { Random = new MersenneTwister(12345), MaxIterations = 1000, MaxFunctionEvaluations = 100 };
            capped.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumFunctionEvaluationsReached, capped.Status);
            Assert.AreEqual(100, capped.FunctionEvaluations);
            var prng = new MersenneTwister(12345);
            int inside = 0;
            for (int i = 0; i < 100; i++)
            {
                double x = -1 + prng.NextDouble() * 2, y = -1 + prng.NextDouble() * 2;
                if (x * x + y * y <= 1) inside++;
            }
            Assert.AreEqual(4.0 * inside / 100, capped.Result, 1E-15);

            // The iteration limit result keeps its bits (3.15512) and is reported as MaximumIterationsReached.
            var iterations = new MonteCarloIntegration(Integrands.PI, 2, a, b) { Random = new MersenneTwister(12345), MaxIterations = 100000 };
            iterations.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, iterations.Status);
            Assert.AreEqual(0x40093DAF8DF7A4E8, BitConverter.DoubleToInt64Bits(iterations.Result));

            // A constant integrand has zero standard error and converges on the second sample.
            var constant = new MonteCarloIntegration(x => 2.0, 2, a, b) { Random = new MersenneTwister(12345) };
            constant.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, constant.Status);
            Assert.AreEqual(2, constant.FunctionEvaluations);
            Assert.AreEqual(8.0, constant.Result, 0.0);
        }

        /// <summary>
        /// UseSobolSequence defaults to false (the pseudo-random path). When true, each call samples a
        /// Sobol sequence from its first point, so repeated calls agree and match a hand count.
        /// </summary>
        [TestMethod]
        public void Test_SobolSequence()
        {
            var a = new double[] { -1, -1 };
            var b = new double[] { 1, 1 };
            var mc = new MonteCarloIntegration(Integrands.PI, 2, a, b) { MinIterations = 4096, MaxIterations = 4096 };
            Assert.IsFalse(mc.UseSobolSequence);
            mc.UseSobolSequence = true;
            mc.Integrate();
            double first = mc.Result;
            mc.Integrate();
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(first), BitConverter.DoubleToInt64Bits(mc.Result));

            var sobol = new SobolSequence(2);
            int inside = 0;
            for (int i = 0; i < 4096; i++)
            {
                var u = sobol.NextDouble();
                double x = -1 + u[0] * 2, y = -1 + u[1] * 2;
                if (x * x + y * y <= 1) inside++;
            }
            Assert.AreEqual(4.0 * inside / 4096, first, 1E-15);
            Assert.AreEqual(Math.PI, first, 5E-2);
        }
    }
}
