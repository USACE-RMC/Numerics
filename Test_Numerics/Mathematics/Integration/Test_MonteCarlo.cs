using Microsoft.VisualStudio.TestTools.UnitTesting;
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
    }
}
