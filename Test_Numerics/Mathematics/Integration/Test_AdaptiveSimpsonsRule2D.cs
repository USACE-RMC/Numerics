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
    public class Test_AdaptiveSimpsonsRule2D
    {

        /// <summary>
        /// Test the 2D Adaptive Simpson's Rule algorithm with the Pi function
        /// </summary>
        [TestMethod()]
        public void Test_PI()
        {
            var asr2D = new AdaptiveSimpsonsRule2D(Integrands.PI2D, -1, 1, -1, 1 );
            asr2D.Integrate();
            var result = asr2D.Result;
            double trueResult = 3.14;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the 2D Adaptive Simpson's Rule algorithm with the sum of two normal distributions.
        /// </summary>
        [TestMethod()]
        public void Test_SumOfTwoNormals()
        {
            var asr2D = new AdaptiveSimpsonsRule2D(Integrands.SumOfNormals2D, 1E-15, 1- 1E-15, 1E-15, 1- 1E-15);
            asr2D.Integrate();
            var result = asr2D.Result;
            double trueResult = 40;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the 2D Adaptive Simpson's Rule algorithm with X + Y integrand.
        /// </summary>
        [TestMethod()]
        public void Test_XPlusY()
        {
            // Define the function f(x, y) = x + y
            Func<double, double, double> func = (x, y) => x + y;

            // Set up the 2D Adaptive Simpson's Rule with the bounds [-1, 1] for both x and y
            var asr2D = new AdaptiveSimpsonsRule2D(func, -1, 1, -1, 1);

            // Run the integration process
            asr2D.Integrate();

            // Get the result of the integration
            var result = asr2D.Result;

            // The exact result is 0
            double trueResult = 0;

            // Assert that the result is within a small error margin of 0
            Assert.AreEqual(trueResult, result, 1E-5);
        }

        /// <summary>
        /// Test the 2D Adaptive Simpson's Rule algorithm with X^2 + Y^2 integrand.
        /// </summary>
        [TestMethod()]
        public void Test_XSquaredPlusYSquared()
        {
            // Define the function f(x, y) = x^2 + y^2
            Func<double, double, double> func = (x, y) => x * x + y * y;

            // Set up the 2D Adaptive Simpson's Rule with the bounds [-1, 1] for both x and y
            var asr2D = new AdaptiveSimpsonsRule2D(func, -1, 1, -1, 1);

            // Run the integration process
            asr2D.Integrate();

            // Get the result of the integration
            var result = asr2D.Result;

            // The exact result is 2.666...
            double trueResult = 2.6666667;

            // Assert that the result is within a small error margin of 2
            Assert.AreEqual(trueResult, result, 1E-6);
        }

    }
}
