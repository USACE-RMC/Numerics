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

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Ln-Normal (Galton) distribution algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil</item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests" />
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_LnNormal
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 1.8.1 Wabash River at Lafayette, IN
        private double[] sample = new double[] { 41500d, 57000d, 44000d, 49000d, 31000d, 45900d, 19000d, 41100d, 37300d, 76000d, 33200d, 61200d, 76000d, 59800d, 44400d, 58400d, 53600d, 59800d, 63300d, 57700d, 64000d, 63500d, 38000d, 74600d, 13100d, 37600d, 67500d, 21700d, 37000d, 93500d, 58500d, 63300d, 74400d, 34200d, 14600d, 44200d, 13100d, 73300d, 46600d, 39400d, 41200d, 41300d, 62000d, 90000d, 50600d, 41900d, 35000d, 16500d, 35300d, 30000d, 52600d, 99000d, 89000d, 39500d, 55400d, 46000d, 63000d, 58300d, 36500d, 14600d, 64900d, 68500d, 69100d, 42600d, 31000d, 39400d, 40700d, 53400d, 36000d, 43900d, 23600d, 50500d, 49700d, 48100d, 44500d, 56400d, 60800d, 40400d, 80400d, 41600d, 14700d, 33300d, 40700d, 53300d, 77400d };

        /// <summary>
        /// Verification of Ln-Normal Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.2.1 page 99.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_LnNormal_MOM_Fit()
        {

            var LN = new LnNormal();
            LN.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double u1 = LN.Mu;
            double u2 = LN.Sigma;
            double true_u1 = 10.7676d;
            double true_u2 = 0.4544d;
            Assert.AreEqual((u1 - true_u1) / true_u1 < 0.01d, true);
            Assert.AreEqual((u2 - true_u2) / true_u2 < 0.01d, true);
        }

        /// <summary>
        /// Verification of Ln-Normal Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.2.1 page 99.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_LnNormal_MLE_Fit()
        {
            var LN = new LnNormal();
            LN.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double u1 = LN.Mu;
            double u2 = LN.Sigma;
            double true_u1 = 10.7711d;
            double true_u2 = 0.4562d;
            Assert.AreEqual((u1 - true_u1) / true_u1 < 0.01d, true);
            Assert.AreEqual((u2 - true_u2) / true_u2 < 0.01d, true);
        }

        /// <summary>
        /// Test the quantile function for the Ln-Normal Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.2.2 page 102.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_LnNormal_Quantile()
        {
            var LN = new LnNormal() { Mu = 10.7676d, Sigma = 0.4544d };
            double q100 = LN.InverseCDF(0.99d);
            double true_q100 = 136611d;
            Assert.AreEqual((q100 - true_q100) / true_q100 < 0.01d, true);
            double p = LN.CDF(q100);
            double true_p = 0.99d;
            Assert.AreEqual((p - true_p) / true_p < 0.01d, true);
        }

        /// <summary>
        /// Test the standard error for the Ln-Normal Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.1.3 page 94.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_LnNormal_StandardError()
        {
            // Maximum Likelihood
            var LN = new LnNormal() { Mu = 10.7711d, Sigma = 0.4562d };
            LN.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double qVar99 = Math.Sqrt(LN.QuantileVariance(0.99d, 85, ParameterEstimationMethod.MaximumLikelihood));
            double true_qVar99 = 13113d;
            Assert.AreEqual((qVar99 - true_qVar99) / true_qVar99 < 0.01d, true);
        }

        /// <summary>
        /// Test the conditional expectation using textbook solution.
        /// </summary>
        [TestMethod()]
        public void Test_ConditionalExpectation()
        {
            double mu = 10d;
            double sigma = 5d;
            double alpha = 0.999;
            var LN = new LnNormal(mu, sigma);

            double CE = LN.ConditionalExpectedValue(alpha);
            double true_CE = Math.Exp(LN.Mu + 0.5 * LN.Sigma * LN.Sigma) / (1d - alpha) * (1 - Normal.StandardCDF(Normal.StandardZ(alpha) - LN.Sigma));
            Assert.AreEqual(CE, true_CE, 1E-4);

        }

        /// <summary>
        /// Verifying that inputs can create an Ln Normal distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.Mean, 10, 1E-4);
            Assert.AreEqual(LN.StandardDeviation, 10, 1E-4);

            var LN2 = new LnNormal(1, 1);
            Assert.AreEqual(LN2.Mean, 1, 1E-4);
            Assert.AreEqual(LN2.StandardDeviation, 1, 1E-4);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var LN = new LnNormal(double.NaN, 0);
            Assert.IsFalse(LN.ParametersValid);

            var LN2 = new LnNormal(double.PositiveInfinity, 1);
            Assert.IsFalse(LN2.ParametersValid);

            var LN3 = new LnNormal(0,double.PositiveInfinity);
            Assert.IsFalse(LN3.ParametersValid);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new LnNormal(10, 5);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.Mean, 1.142e26, 1e30);

            var LN2 = new LnNormal(1, 1);
            Assert.AreEqual(LN2.Mean, 1,1e-04);
        }

        /// <summary>
        /// Testing median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.Median, 7.07106,1e-05);

            var LN2 = new LnNormal(1, 1);
            Assert.AreEqual(LN2.Median, 0.707106, 1e-05);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.Mode, 3.5355, 1e-04);

            var LN2 = new LnNormal(1, 1);
            Assert.AreEqual(LN2.Mode, 0.35355,1e-04);
        }

        /// <summary>
        /// Testing standard deviation function.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.StandardDeviation, 5.92e47, 1e49);

            var LN2 = new LnNormal(1, 1);
            Assert.AreEqual(LN2.StandardDeviation, 1,1e-4);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.Skewness, 1.39e65, 1e67);

            var LN2 = new LnNormal(1, 1);
            Assert.AreEqual(LN2.Skewness, 4, 1e-04);
        }

        /// <summary>
        /// Testing Kurtosis
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var LN = new LnNormal(1,1);
            Assert.AreEqual(LN.Kurtosis, 41, 1e-04);
        }

        /// <summary>
        /// Testing Minimum and Maximum functions are 0 and infinity, respectively.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.Minimum, 0);
            Assert.AreEqual(LN.Maximum,double.PositiveInfinity);

            var LN2 = new LnNormal(1, 1);
            Assert.AreEqual(LN2.Minimum, 0);
            Assert.AreEqual(LN2.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.PDF(1), 0.03033, 1e-04);
            Assert.AreEqual(LN.PDF(-1), 0);

            var LN2 = new LnNormal(2.5, 2.5);
            Assert.AreEqual(LN2.PDF(0.5), 0.303322, 1e-04);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var LN = new LnNormal(2.5,2.5);
            Assert.AreEqual(LN.CDF(0.5), 0.06465, 1e-05);
            Assert.AreEqual(LN.CDF(0.8), 0.17046, 1e-05);
        }

        /// <summary>
        /// Testing inverse cdf method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var LN = new LnNormal();
            Assert.AreEqual(LN.InverseCDF(0), 0);
            Assert.AreEqual(LN.InverseCDF(1),double.PositiveInfinity);
            Assert.AreEqual(LN.InverseCDF(0.5), 7.07106,1e-04);
        }
    }
}
