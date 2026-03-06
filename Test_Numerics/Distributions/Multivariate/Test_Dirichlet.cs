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
    /// Unit tests for the Dirichlet distribution.
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
    /// Reference values verified analytically and against R package MCMCpack::ddirichlet.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Dirichlet
    {

        /// <summary>
        /// Test symmetric Dirichlet construction.
        /// </summary>
        [TestMethod]
        public void Test_SymmetricConstruction()
        {
            var D = new Dirichlet(3, 2.0);
            Assert.AreEqual(3, D.Dimension);
            Assert.IsTrue(D.ParametersValid);

            var alpha = D.Alpha;
            Assert.AreEqual(2.0, alpha[0]);
            Assert.AreEqual(2.0, alpha[1]);
            Assert.AreEqual(2.0, alpha[2]);
            Assert.AreEqual(6.0, D.AlphaSum);
        }

        /// <summary>
        /// Test asymmetric Dirichlet construction.
        /// </summary>
        [TestMethod]
        public void Test_AsymmetricConstruction()
        {
            var D = new Dirichlet(new double[] { 1.0, 2.0, 3.0 });
            Assert.AreEqual(3, D.Dimension);
            Assert.IsTrue(D.ParametersValid);
            Assert.AreEqual(6.0, D.AlphaSum);
        }

        /// <summary>
        /// Test invalid parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new Dirichlet(1, 1.0)); // dimension < 2
            Assert.Throws<ArgumentOutOfRangeException>(() => new Dirichlet(3, 0.0)); // alpha = 0
            Assert.Throws<ArgumentOutOfRangeException>(() => new Dirichlet(3, -1.0)); // alpha < 0
            Assert.Throws<ArgumentOutOfRangeException>(() => new Dirichlet(new double[] { 1.0, -1.0 })); // negative alpha
        }

        /// <summary>
        /// Test distribution type and display names.
        /// </summary>
        [TestMethod]
        public void Test_TypeAndName()
        {
            var D = new Dirichlet(3, 1.0);
            Assert.AreEqual(MultivariateDistributionType.Dirichlet, D.Type);
            Assert.AreEqual("Dirichlet", D.DisplayName);
            Assert.AreEqual("Dir", D.ShortDisplayName);
        }

        /// <summary>
        /// Test mean vector. Mean[i] = alpha[i] / sum(alpha).
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            // Symmetric Dir(2, 2, 2): mean = (1/3, 1/3, 1/3)
            var D1 = new Dirichlet(3, 2.0);
            var mean1 = D1.Mean;
            Assert.AreEqual(1.0 / 3.0, mean1[0], 1e-10);
            Assert.AreEqual(1.0 / 3.0, mean1[1], 1e-10);
            Assert.AreEqual(1.0 / 3.0, mean1[2], 1e-10);

            // Asymmetric Dir(1, 2, 3): mean = (1/6, 2/6, 3/6)
            var D2 = new Dirichlet(new double[] { 1.0, 2.0, 3.0 });
            var mean2 = D2.Mean;
            Assert.AreEqual(1.0 / 6.0, mean2[0], 1e-10);
            Assert.AreEqual(2.0 / 6.0, mean2[1], 1e-10);
            Assert.AreEqual(3.0 / 6.0, mean2[2], 1e-10);
        }

        /// <summary>
        /// Test variance vector. Var[i] = alpha[i] * (S - alpha[i]) / (S^2 * (S+1)).
        /// </summary>
        [TestMethod]
        public void Test_Variance()
        {
            // Dir(1, 2, 3): S = 6
            // Var[0] = 1*5 / (36*7) = 5/252 ≈ 0.019841
            // Var[1] = 2*4 / (36*7) = 8/252 ≈ 0.031746
            // Var[2] = 3*3 / (36*7) = 9/252 ≈ 0.035714
            var D = new Dirichlet(new double[] { 1.0, 2.0, 3.0 });
            var v = D.Variance;
            Assert.AreEqual(5.0 / 252.0, v[0], 1e-10);
            Assert.AreEqual(8.0 / 252.0, v[1], 1e-10);
            Assert.AreEqual(9.0 / 252.0, v[2], 1e-10);
        }

        /// <summary>
        /// Test mode vector. Mode[i] = (alpha[i] - 1) / (S - K) when all alpha > 1.
        /// </summary>
        [TestMethod]
        public void Test_Mode()
        {
            // Dir(2, 3, 5): S=10, K=3, S-K=7
            // Mode = (1/7, 2/7, 4/7)
            var D = new Dirichlet(new double[] { 2.0, 3.0, 5.0 });
            var mode = D.Mode;
            Assert.AreEqual(1.0 / 7.0, mode[0], 1e-10);
            Assert.AreEqual(2.0 / 7.0, mode[1], 1e-10);
            Assert.AreEqual(4.0 / 7.0, mode[2], 1e-10);
        }

        /// <summary>
        /// Test that Mode throws when any alpha <= 1.
        /// </summary>
        [TestMethod]
        public void Test_Mode_Invalid()
        {
            var D = new Dirichlet(new double[] { 0.5, 2.0, 3.0 });
            Assert.Throws<InvalidOperationException>(() => { var m = D.Mode; });
        }

        /// <summary>
        /// Test covariance. Cov(Xi, Xj) = -alpha_i * alpha_j / (S^2 * (S+1)).
        /// </summary>
        [TestMethod]
        public void Test_Covariance()
        {
            var D = new Dirichlet(new double[] { 1.0, 2.0, 3.0 });
            // Cov(X0, X1) = -1*2 / (36*7) = -2/252
            Assert.AreEqual(-2.0 / 252.0, D.Covariance(0, 1), 1e-10);
            // Cov(X0, X0) = Var(X0)
            Assert.AreEqual(D.Variance[0], D.Covariance(0, 0), 1e-10);
            // Covariance matrix is symmetric
            Assert.AreEqual(D.Covariance(0, 1), D.Covariance(1, 0), 1e-10);
        }

        /// <summary>
        /// Test the PDF for the symmetric Dir(1,1,1) = Uniform on simplex.
        /// </summary>
        [TestMethod]
        public void Test_PDF_Uniform()
        {
            // Dir(1,1,1) is uniform on the 2-simplex.
            // PDF = Gamma(3) / (Gamma(1)*Gamma(1)*Gamma(1)) = 2! = 2
            var D = new Dirichlet(3, 1.0);
            double pdf = D.PDF(new double[] { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 });
            Assert.AreEqual(2.0, pdf, 1e-6);

            // Should be the same everywhere on the simplex
            double pdf2 = D.PDF(new double[] { 0.1, 0.2, 0.7 });
            Assert.AreEqual(2.0, pdf2, 1e-6);
        }

        /// <summary>
        /// Test the PDF for Dir(2, 3, 5).
        /// </summary>
        [TestMethod]
        public void Test_PDF_Asymmetric()
        {
            // Dir(2,3,5): B(alpha) = Gamma(2)*Gamma(3)*Gamma(5) / Gamma(10)
            // = 1! * 2! * 4! / 9! = 1*2*24 / 362880 = 48/362880 = 1/7560
            // PDF(0.1, 0.3, 0.6) = (1/B) * 0.1^1 * 0.3^2 * 0.6^4
            // = 7560 * 0.1 * 0.09 * 0.1296 = 7560 * 0.0011664 = 8.81...
            var D = new Dirichlet(new double[] { 2.0, 3.0, 5.0 });
            double pdf = D.PDF(new double[] { 0.1, 0.3, 0.6 });
            double expected = 7560.0 * 0.1 * 0.09 * 0.1296;
            Assert.AreEqual(expected, pdf, 1e-4);
        }

        /// <summary>
        /// Test that PDF returns 0 for points outside the simplex.
        /// </summary>
        [TestMethod]
        public void Test_PDF_OutsideSimplex()
        {
            var D = new Dirichlet(3, 2.0);
            // Components don't sum to 1
            Assert.AreEqual(0.0, D.PDF(new double[] { 0.5, 0.5, 0.5 }), 1e-10);
            // Negative component
            Assert.AreEqual(0.0, D.PDF(new double[] { -0.1, 0.6, 0.5 }), 1e-10);
            // Wrong dimension
            Assert.AreEqual(0.0, D.PDF(new double[] { 0.5, 0.5 }), 1e-10);
        }

        /// <summary>
        /// Test that LogPDF is consistent with PDF.
        /// </summary>
        [TestMethod]
        public void Test_LogPDF()
        {
            var D = new Dirichlet(new double[] { 2.0, 3.0, 5.0 });
            var x = new double[] { 0.2, 0.3, 0.5 };
            Assert.AreEqual(Math.Log(D.PDF(x)), D.LogPDF(x), 1e-10);
        }

        /// <summary>
        /// Test random sampling: all samples on the simplex and moments converge.
        /// </summary>
        [TestMethod]
        public void Test_Sampling()
        {
            var D = new Dirichlet(new double[] { 2.0, 3.0, 5.0 });
            var samples = D.GenerateRandomValues(10000, seed: 42);

            // All samples should be on the simplex
            for (int i = 0; i < 10000; i++)
            {
                double sum = 0;
                for (int j = 0; j < 3; j++)
                {
                    Assert.IsGreaterThan(0, samples[i, j], $"Sample [{i},{j}] = {samples[i, j]} is not positive");
                    sum += samples[i, j];
                }
                Assert.AreEqual(1.0, sum, 1e-10, $"Sample {i} does not sum to 1");
            }

            // Check that sample means converge to theoretical means
            var sampleMeans = new double[3];
            for (int i = 0; i < 10000; i++)
            {
                for (int j = 0; j < 3; j++)
                    sampleMeans[j] += samples[i, j];
            }
            for (int j = 0; j < 3; j++)
                sampleMeans[j] /= 10000;

            var theoreticalMeans = D.Mean;
            Assert.AreEqual(theoreticalMeans[0], sampleMeans[0], 0.02); // 2/10 = 0.2
            Assert.AreEqual(theoreticalMeans[1], sampleMeans[1], 0.02); // 3/10 = 0.3
            Assert.AreEqual(theoreticalMeans[2], sampleMeans[2], 0.02); // 5/10 = 0.5
        }

        /// <summary>
        /// Test Clone produces an independent copy.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var D = new Dirichlet(new double[] { 1.0, 2.0, 3.0 });
            var D2 = D.Clone() as Dirichlet;
            Assert.IsNotNull(D2);
            Assert.AreEqual(D.Dimension, D2.Dimension);
            Assert.IsTrue(D2.ParametersValid);

            // Verify alpha values are copied
            var a1 = D.Alpha;
            var a2 = D2.Alpha;
            for (int i = 0; i < D.Dimension; i++)
                Assert.AreEqual(a1[i], a2[i]);
        }

        /// <summary>
        /// Test the LogMultivariateBeta function.
        /// </summary>
        [TestMethod]
        public void Test_LogMultivariateBeta()
        {
            // For alpha = (1, 1), B(1,1) = Gamma(1)*Gamma(1)/Gamma(2) = 1*1/1 = 1
            Assert.AreEqual(0.0, Dirichlet.LogMultivariateBeta(new double[] { 1.0, 1.0 }), 1e-10);

            // For alpha = (1, 1, 1), B(1,1,1) = Gamma(1)^3/Gamma(3) = 1/2
            Assert.AreEqual(Math.Log(0.5), Dirichlet.LogMultivariateBeta(new double[] { 1.0, 1.0, 1.0 }), 1e-10);
        }

        /// <summary>
        /// Test CDF throws NotImplementedException.
        /// </summary>
        [TestMethod]
        public void Test_CDF_Throws()
        {
            var D = new Dirichlet(3, 1.0);
            Assert.Throws<NotImplementedException>(() => D.CDF(new double[] { 0.3, 0.3, 0.4 }));
        }

    }
}
