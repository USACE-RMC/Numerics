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
    /// Testing the Log-Normal distribution algorithm.
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
    public class Test_LogNormal
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 1.8.1 Wabash River at Lafayette, IN
        private double[] sample = new double[] { 41500d, 57000d, 44000d, 49000d, 31000d, 45900d, 19000d, 41100d, 37300d, 76000d, 33200d, 61200d, 76000d, 59800d, 44400d, 58400d, 53600d, 59800d, 63300d, 57700d, 64000d, 63500d, 38000d, 74600d, 13100d, 37600d, 67500d, 21700d, 37000d, 93500d, 58500d, 63300d, 74400d, 34200d, 14600d, 44200d, 13100d, 73300d, 46600d, 39400d, 41200d, 41300d, 62000d, 90000d, 50600d, 41900d, 35000d, 16500d, 35300d, 30000d, 52600d, 99000d, 89000d, 39500d, 55400d, 46000d, 63000d, 58300d, 36500d, 14600d, 64900d, 68500d, 69100d, 42600d, 31000d, 39400d, 40700d, 53400d, 36000d, 43900d, 23600d, 50500d, 49700d, 48100d, 44500d, 56400d, 60800d, 40400d, 80400d, 41600d, 14700d, 33300d, 40700d, 53300d, 77400d };

        /// <summary>
        /// Verification of Log-Normal Distribution fit with method of moments.
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
        public void Test_LogNormal_MOM_Fit()
        {
            var LogN = new LogNormal() { Base = Math.E };
            LogN.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double u1 = LogN.Mu;
            double u2 = LogN.Sigma;
            double true_u1 = 10.716952223744224d;
            double true_u2 = 0.45007398831588075d;
            Assert.AreEqual((u1 - true_u1) / true_u1 < 0.01d, true);
            Assert.AreEqual((u2 - true_u2) / true_u2 < 0.01d, true);
        }

        /// <summary>
        /// Verification of Log-Normal Distribution fit with method of linear moments.
        /// </summary>
        [TestMethod()]
        public void Test_LogNormal_LMOM_Fit()
        {
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var norm = new LogNormal();
            norm.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
            double u1 = norm.Mu;
            double u2 = norm.Sigma;
            double true_u1 = 0.9672391d;
            double true_u2 = 0.1675344d;
            Assert.AreEqual(u1, true_u1, 0.0001d);
            Assert.AreEqual(u2, true_u2, 0.0001d);
            var lmom = norm.LinearMomentsFromParameters(norm.GetParameters);
            Assert.AreEqual(lmom[0], 0.96723909d, 0.0001d);
            Assert.AreEqual(lmom[1], 0.09452119d, 0.0001d);
            Assert.AreEqual(lmom[2], 0.00000000d, 0.0001d);
            Assert.AreEqual(lmom[3], 0.12260172d, 0.0001d);
        }

        /// <summary>
        /// Verification of Log-Normal Distribution fit with method of maximum likelihood.
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
        public void Test_LogNormal_MLE_Fit()
        {
            var LogN = new LogNormal() { Base = Math.E };
            LogN.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double u1 = LogN.Mu;
            double u2 = LogN.Sigma;
            double true_u1 = 10.716950857801747d;
            double true_u2 = 0.44742859657407796d;
            Assert.AreEqual((u1 - true_u1) / true_u1 < 0.01d, true);
            Assert.AreEqual((u2 - true_u2) / true_u2 < 0.01d, true);
            var LN = new LnNormal();
            LN.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
        }

        /// <summary>
        /// Test the quantile function for the Log-Normal Distribution.
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
        public void Test_LogNormal_Quantile()
        {
            var LogN = new LogNormal() { Mu = 10.7676d, Sigma = 0.4544d, Base = Math.E };
            double q100 = LogN.InverseCDF(0.99d);
            double true_q100 = 136611d;
            Assert.AreEqual((q100 - true_q100) / true_q100 < 0.01d, true);
            double p = LogN.CDF(q100);
            double true_p = 0.99d;
            Assert.AreEqual((p - true_p) / true_p < 0.01d, true);
        }

        /// <summary>
        /// Test the standard error for the Log-Normal Distribution.
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
        public void Test_LogNormal_StandardError()
        {

            // Maximum Likelihood
            var LogN = new LogNormal() { Mu = 10.7711d, Sigma = 0.4562d, Base = Math.E };
            double qVar99 = Math.Sqrt(LogN.QuantileVariance(0.99d, 85, ParameterEstimationMethod.MaximumLikelihood));
            double true_qVar99 = 13113d;
            Assert.AreEqual((qVar99 - true_qVar99) / true_qVar99 < 0.01d, true);
        }

        /// <summary>
        /// Verifying Log-Normal distribution is being created with input parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var LogN = new LogNormal();
            Assert.AreEqual(LogN.Mu, 3);
            Assert.AreEqual(LogN.Sigma, 0.5);

            var LogN2 = new LogNormal(1, 1);
            Assert.AreEqual(LogN2.Mu, 1);
            Assert.AreEqual(LogN2.Sigma, 1);
        }

        /// <summary>
        /// Testing Log-Normal with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var LogN = new LogNormal(double.NaN,double.NaN);
            Assert.IsFalse(LogN.ParametersValid);
            
            var LogN2 = new LogNormal(double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(LogN2.ParametersValid);

            var LogN3 = new LogNormal(1, -1);
            Assert.IsFalse(LogN3.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString()
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var LogN = new LogNormal();
            Assert.AreEqual(LogN.ParametersToString[0, 0], "Mean (of log) (µ)");
            Assert.AreEqual(LogN.ParametersToString[1, 0], "Std Dev (of log) (σ)");
            Assert.AreEqual(LogN.ParametersToString[0, 1], "3");
            Assert.AreEqual(LogN.ParametersToString[1, 1], "0.5");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var trueDist = new LnNormal(10, 5);
            var dist = new LogNormal(trueDist.Mu, trueDist.Sigma) { Base = Math.E };
            Assert.AreEqual(trueDist.Mean, dist.Mean, 1E-2);
            Assert.AreEqual(trueDist.StandardDeviation, dist.StandardDeviation, 1E-2);
            Assert.AreEqual(trueDist.Skewness, dist.Skewness, 1E-2);
            Assert.AreEqual(trueDist.Kurtosis, dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing minimum and maximum functions respectively.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var LogN = new LogNormal();
            Assert.AreEqual(LogN.Minimum, 0);
            Assert.AreEqual(LogN.Maximum, double.PositiveInfinity);

            var LogN2 = new LogNormal(1, 1);
            Assert.AreEqual(LogN2.Minimum, 0);
            Assert.AreEqual(LogN2.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var LogN = new LogNormal(1.5,0.1);
            Assert.AreEqual(LogN.PDF(0.1), 3.32e-135,1e-04);

            var LogN2 = new LogNormal(-0.1, 0.1);
            Assert.AreEqual(LogN.PDF(0.8), 9.12888e-56, 1e-04);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var LogN = new LogNormal(1.5, 0.1);
            Assert.AreEqual(LogN.CDF(0.1), 0);

            var LogN2 = new LogNormal(1.5, 1.5);
            Assert.AreEqual(LogN2.CDF(0.5), 0.11493, 1e-05);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var LogN = new LogNormal(2.5, 2.5);
            Assert.AreEqual(LogN.InverseCDF(0.8), 40183.99248, 1e-04);

            var LogN2 = new LogNormal(1.5, 2.5);
            Assert.AreEqual(LogN.InverseCDF(0.8), 40183.99248, 1e-05);
        }
    }
}
