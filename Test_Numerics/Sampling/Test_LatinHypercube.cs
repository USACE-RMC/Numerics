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
using Numerics.Sampling;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using System.Diagnostics;

namespace Sampling
{
    /// <summary>
    /// Unit tests for the Latin Hypercube class. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_LatinHypercube
    {

        /// <summary>
        /// Test for Latin Hypercube Sampling
        /// </summary>
        [TestMethod]
        public void Test_LHS()
        {
            // There is not a great way to test LHS against a known theoretical solution,
            // or other code solutions. So here I am just testing that I get back the correct
            // sampling statistics.

            int N = 10000;
            var norm1 = new Normal(100, 15);
            var norm2 = new Normal(200, 30);
            var x1 = new double[N];
            var x2 = new double[N];

            // LHS
            var lhs = LatinHypercube.Random(N, 2, 45678);
            for (int i = 0; i < N; i++)
            {
                x1[i] = norm1.InverseCDF(lhs[i, 0]);
                x2[i] = norm2.InverseCDF(lhs[i, 1]);
            }

            // Test mean
            Assert.AreEqual(norm1.Mean, Statistics.Mean(x1), 1E-2);
            Assert.AreEqual(norm2.Mean, Statistics.Mean(x2), 1E-2);

            // Test standard deviation
            Assert.AreEqual(norm1.StandardDeviation, Statistics.StandardDeviation(x1), 1E-2);
            Assert.AreEqual(norm2.StandardDeviation, Statistics.StandardDeviation(x2), 1E-2);

            // Test correlation. Should be close to zero.
            Assert.AreEqual(0.0, Correlation.Pearson(x1, x2), 1E-2);

        }

        /// <summary>
        /// Test for Latin Hypercube Sampling using the median value per bin.
        /// </summary>
        [TestMethod]
        public void Test_LHS_Median()
        {
            // There is not a great way to test LHS against a known theoretical solution,
            // or other code solutions. So here I am just testing that I get back the correct
            // sampling statistics.

            int N = 10000;
            var norm1 = new Normal(100, 15);
            var norm2 = new Normal(200, 30);
            var x1 = new double[N];
            var x2 = new double[N];

            // LHS
            var lhs = LatinHypercube.Median(N, 2, 45678);
            for (int i = 0; i < N; i++)
            {
                x1[i] = norm1.InverseCDF(lhs[i, 0]);
                x2[i] = norm2.InverseCDF(lhs[i, 1]);
            }

            // Test mean
            Assert.AreEqual(norm1.Mean, Statistics.Mean(x1), 1E-2);
            Assert.AreEqual(norm2.Mean, Statistics.Mean(x2), 1E-2);

            // Test standard deviation
            Assert.AreEqual(norm1.StandardDeviation, Statistics.StandardDeviation(x1), 1E-2);
            Assert.AreEqual(norm2.StandardDeviation, Statistics.StandardDeviation(x2), 1E-2);

            // Test correlation. Should be close to zero.
            // Correlation is higher with the median option
            Assert.AreEqual(0.0, Correlation.Pearson(x1, x2), 0.05);

        }

    }
}
