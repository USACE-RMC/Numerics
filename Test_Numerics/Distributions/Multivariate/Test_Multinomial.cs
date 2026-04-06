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

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Multivariate
{
    /// <summary>
    /// Unit tests for the Multinomial distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// Reference values verified analytically and against R dmultinom.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Multinomial
    {

        /// <summary>
        /// Test basic construction.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            Assert.AreEqual(3, M.Dimension);
            Assert.AreEqual(10, M.NumberOfTrials);
            Assert.IsTrue(M.ParametersValid);
        }

        /// <summary>
        /// Test invalid parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(0, new double[] { 0.5, 0.5 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(10, new double[] { 0.5 })); // too few categories
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(10, new double[] { -0.1, 0.6, 0.5 })); // negative prob
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(10, new double[] { 0.3, 0.3, 0.3 })); // don't sum to 1
        }

        /// <summary>
        /// Test type and display names.
        /// </summary>
        [TestMethod]
        public void Test_TypeAndName()
        {
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });
            Assert.AreEqual(MultivariateDistributionType.Multinomial, M.Type);
            Assert.AreEqual("Multinomial", M.DisplayName);
            Assert.AreEqual("Mult", M.ShortDisplayName);
        }

        /// <summary>
        /// Test mean vector. Mean[i] = N * p[i].
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            var mean = M.Mean;
            Assert.AreEqual(2.0, mean[0], 1e-10);
            Assert.AreEqual(3.0, mean[1], 1e-10);
            Assert.AreEqual(5.0, mean[2], 1e-10);
        }

        /// <summary>
        /// Test variance vector. Var[i] = N * p[i] * (1 - p[i]).
        /// </summary>
        [TestMethod]
        public void Test_Variance()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            var variance = M.Variance;
            Assert.AreEqual(10 * 0.2 * 0.8, variance[0], 1e-10);
            Assert.AreEqual(10 * 0.3 * 0.7, variance[1], 1e-10);
            Assert.AreEqual(10 * 0.5 * 0.5, variance[2], 1e-10);
        }

        /// <summary>
        /// Test covariance. Cov(Xi, Xj) = -N * pi * pj.
        /// </summary>
        [TestMethod]
        public void Test_Covariance()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            Assert.AreEqual(-10 * 0.2 * 0.3, M.Covariance(0, 1), 1e-10);
            Assert.AreEqual(M.Variance[0], M.Covariance(0, 0), 1e-10);
        }

        /// <summary>
        /// Test PMF for a simple fair coin case: Mult(10, (0.5, 0.5)) = Binomial(10, 0.5).
        /// </summary>
        [TestMethod]
        public void Test_PMF_Binomial()
        {
            // Mult(10, (0.5, 0.5)) at x = (5, 5) should equal C(10,5) * 0.5^10
            // = 252 * 0.0009765625 = 0.24609375
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });
            double pmf = M.PDF(new double[] { 5, 5 });
            Assert.AreEqual(0.24609375, pmf, 1e-8);
        }

        /// <summary>
        /// Test PMF for a 3-category case.
        /// </summary>
        [TestMethod]
        public void Test_PMF_ThreeCategory()
        {
            // Mult(4, (0.2, 0.3, 0.5)) at x = (1, 1, 2)
            // = 4! / (1! 1! 2!) * 0.2^1 * 0.3^1 * 0.5^2
            // = 24 / (1*1*2) * 0.2 * 0.3 * 0.25
            // = 12 * 0.015 = 0.18
            var M = new Multinomial(4, new double[] { 0.2, 0.3, 0.5 });
            double pmf = M.PDF(new double[] { 1, 1, 2 });
            Assert.AreEqual(0.18, pmf, 1e-10);
        }

        /// <summary>
        /// Test that PMF returns 0 for invalid count vectors.
        /// </summary>
        [TestMethod]
        public void Test_PMF_Invalid()
        {
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });
            // Counts don't sum to N
            Assert.AreEqual(0.0, M.PDF(new double[] { 3, 3 }), 1e-10);
            // Negative count
            Assert.AreEqual(0.0, M.PDF(new double[] { -1, 11 }), 1e-10);
            // Wrong dimension
            Assert.AreEqual(0.0, M.PDF(new double[] { 5, 3, 2 }), 1e-10);
        }

        /// <summary>
        /// Test LogPMF consistency with PMF.
        /// </summary>
        [TestMethod]
        public void Test_LogPMF()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            var x = new double[] { 2, 3, 5 };
            Assert.AreEqual(Math.Log(M.PDF(x)), M.LogPMF(x), 1e-10);
        }

        /// <summary>
        /// Test random sampling: each sample sums to N and moments converge.
        /// </summary>
        [TestMethod]
        public void Test_Sampling()
        {
            var M = new Multinomial(100, new double[] { 0.2, 0.3, 0.5 });
            var samples = M.GenerateRandomValues(5000, seed: 42);

            // All samples should sum to N
            for (int i = 0; i < 5000; i++)
            {
                double sum = 0;
                for (int j = 0; j < 3; j++)
                {
                    Assert.IsGreaterThanOrEqualTo(0, samples[i, j], $"Sample [{i},{j}] is negative");
                    sum += samples[i, j];
                }
                Assert.AreEqual(100.0, sum, 1e-10, $"Sample {i} does not sum to N");
            }

            // Sample means should converge to theoretical means
            var sampleMeans = new double[3];
            for (int i = 0; i < 5000; i++)
                for (int j = 0; j < 3; j++)
                    sampleMeans[j] += samples[i, j];

            var theoreticalMeans = M.Mean;
            for (int j = 0; j < 3; j++)
            {
                sampleMeans[j] /= 5000;
                Assert.AreEqual(theoreticalMeans[j], sampleMeans[j], 1.0); // within 1 count of N*p
            }
        }

        /// <summary>
        /// Test the static Sample method for weighted categorical sampling.
        /// </summary>
        [TestMethod]
        public void Test_WeightedSample()
        {
            var rng = new Random(42);
            var weights = new double[] { 1.0, 3.0, 6.0 }; // 10%, 30%, 60%
            var counts = new int[3];

            int n = 10000;
            for (int i = 0; i < n; i++)
            {
                int idx = Multinomial.Sample(weights, rng);
                Assert.IsGreaterThanOrEqualTo(0, idx);
                Assert.IsLessThan(3, idx);
                counts[idx]++;
            }

            // Check proportions
            Assert.AreEqual(0.1, counts[0] / (double)n, 0.02);
            Assert.AreEqual(0.3, counts[1] / (double)n, 0.02);
            Assert.AreEqual(0.6, counts[2] / (double)n, 0.02);
        }

        /// <summary>
        /// Test the static Sample method with edge cases.
        /// </summary>
        [TestMethod]
        public void Test_WeightedSample_EdgeCases()
        {
            var rng = new Random(42);
            // All weight on one category
            Assert.AreEqual(1, Multinomial.Sample(new double[] { 0.0, 1.0, 0.0 }, rng));

            // Single positive weight
            Assert.AreEqual(0, Multinomial.Sample(new double[] { 5.0, 0.0 }, rng));

            // Invalid: all zeros
            Assert.Throws<ArgumentException>(() => Multinomial.Sample(new double[] { 0.0, 0.0 }, rng));
        }

        /// <summary>
        /// Test Clone produces an independent copy.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            var M2 = M.Clone() as Multinomial;
            Assert.IsNotNull(M2);
            Assert.AreEqual(M.Dimension, M2.Dimension);
            Assert.AreEqual(M.NumberOfTrials, M2.NumberOfTrials);
            Assert.IsTrue(M2.ParametersValid);
        }

        /// <summary>
        /// Test CDF throws NotImplementedException.
        /// </summary>
        [TestMethod]
        public void Test_CDF_Throws()
        {
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });
            Assert.Throws<NotImplementedException>(() => M.CDF(new double[] { 5, 5 }));
        }

    }
}
