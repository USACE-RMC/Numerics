using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the Miser algorithm
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_Miser
    {
        /// <summary>
        /// Test the Miser algorithm with the Pi function
        /// </summary>
        [TestMethod()]
        public void Test_PI()
        {
            var miser = new Miser(Integrands.PI, 2, new double[] { -1, -1 }, new double[] { 1, 1 });
            miser.Random = new MersenneTwister(12345);
            miser.Integrate();
            var result = miser.Result;
            double trueResult = 3.14;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }


        /// <summary>
        /// Test the Miser algorithm with the GSL function
        /// </summary>
        [TestMethod()]
        public void Test_GSL()
        {

            var miser = new Miser(Integrands.GSL, 3, new double[] { 0, 0, 0 }, new double[] { Math.PI, Math.PI, Math.PI });
            miser.Random = new MersenneTwister(12345);
            miser.Integrate();
            var result = miser.Result;
            double trueResult = 1.3932039296856768591842462603255;
            Assert.AreEqual(trueResult, result, 1E-2 * trueResult);
        }

        /// <summary>
        /// Test the Miser algorithm with the sum of three normal distributions
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

            var miser = new Miser(Integrands.SumOfNormals, 3, min, max);
            miser.Random = new MersenneTwister(12345);
            miser.Integrate();
            var result = miser.Result;
            double trueResult = 57;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Miser algorithm with the sum of five normal distributions
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

            var miser = new Miser(Integrands.SumOfNormals, 5, min, max);
            miser.Random = new MersenneTwister(12345);
            miser.Integrate();
            var result = miser.Result;
            double trueResult = 224;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Miser algorithm with the sum of 20 normal distributions
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

            var miser = new Miser(Integrands.SumOfNormals, 20, min, max);
            miser.Random = new MersenneTwister(12345);
            miser.Integrate();
            var result = miser.Result;
            double trueResult = 837;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }
    }
}
