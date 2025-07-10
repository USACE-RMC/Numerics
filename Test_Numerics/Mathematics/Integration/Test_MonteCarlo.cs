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
