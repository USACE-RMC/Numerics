using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics;
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

        /// <summary>
        /// Miser reports Success, and repeated calls on one instance reproduce the first call bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_Status_And_RepeatableCalls()
        {
            var miser = new Miser(Integrands.PI, 2, new double[] { -1, -1 }, new double[] { 1, 1 }) { MaxFunctionEvaluations = 100000 };
            miser.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, miser.Status);
            Assert.AreEqual(100000, miser.FunctionEvaluations);
            Assert.AreEqual(0x40092938B4938B49, BitConverter.DoubleToInt64Bits(miser.Result));
            miser.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, miser.Status);
            Assert.AreEqual(0x40092938B4938B49, BitConverter.DoubleToInt64Bits(miser.Result));
        }

        /// <summary>
        /// Per-depth work arrays reproduce the allocating implementation bit for bit on both sampling
        /// paths, including on a second call that reuses the arrays.
        /// </summary>
        [TestMethod]
        public void Test_PooledWorkArrays_BitsUnchanged()
        {
            var mt = new Miser(Integrands.PI, 2, new double[] { -1, -1 }, new double[] { 1, 1 })
            {
                MaxFunctionEvaluations = 100000,
                UseSobolSequence = false,
                Random = new MersenneTwister(12345)
            };
            mt.Integrate();
            Assert.AreEqual(0x400910DC370DC371, BitConverter.DoubleToInt64Bits(mt.Result));

            var gsl = new Miser(Integrands.GSL, 3, new double[] { 0, 0, 0 }, new double[] { Math.PI, Math.PI, Math.PI }) { MaxFunctionEvaluations = 200000 };
            for (int run = 0; run < 2; run++)
            {
                gsl.Integrate();
                Assert.AreEqual(0x3FF64657C1CDAFC3, BitConverter.DoubleToInt64Bits(gsl.Result));
                Assert.AreEqual(200000, gsl.FunctionEvaluations);
            }
        }
    }
}
