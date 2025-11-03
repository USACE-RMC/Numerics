/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
