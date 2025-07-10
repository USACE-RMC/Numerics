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
    /// Testing the Normal distribution algorithm.
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
    public class Test_Normal
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 5.1.1 Tippecanoe River Near Delphi, Indiana (Station 43) Data
        private double[] sample = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

        /// <summary>
        /// Verification of Normal Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.1.1 page 88.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Normal_MOM_Fit()
        {
            var z = Normal.StandardZ(1 - 1E-16);

            var norm = new Normal();
            norm.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double u1 = norm.Mu;
            double u2 = norm.Sigma;
            double true_u1 = 12665d;
            double true_u2 = 4710d;
            Assert.AreEqual((u1 - true_u1) / true_u1 < 0.01d, true);
            Assert.AreEqual((u2 - true_u2) / true_u2 < 0.01d, true);

        }

        /// <summary>
        /// Verification of Normal Distribution fit with method of linear moments.
        /// </summary>
        [TestMethod()]
        public void Test_Normal_LMOM_Fit()
        {
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var norm = new Normal();
            norm.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
            double u1 = norm.Mu;
            double u2 = norm.Sigma;
            double true_u1 = 9.957516d;
            double true_u2 = 3.513431d;
            Assert.AreEqual(u1, true_u1, 0.0001d);
            Assert.AreEqual(u2, true_u2, 0.0001d);
            var lmom = norm.LinearMomentsFromParameters(norm.GetParameters);
            Assert.AreEqual(lmom[0], 9.9575163d, 0.0001d);
            Assert.AreEqual(lmom[1], 1.9822411d, 0.0001d);
            Assert.AreEqual(lmom[2], 0.0000000d, 0.0001d);
            Assert.AreEqual(lmom[3], 0.1226017d, 0.0001d);
        }

        /// <summary>
        /// Verification of Normal Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// Sigma should be the population standard deviation.
        /// </remarks>
        [TestMethod()]
        public void Test_Normal_MLE_Fit()
        {
            var norm = new Normal();
            norm.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double u1 = norm.Mu;
            double u2 = norm.Sigma;
            double true_u1 = 12665d;
            double true_u2 = 4660d;
            Assert.AreEqual((u1 - true_u1) / true_u1 < 0.01d, true);
            Assert.AreEqual((u2 - true_u2) / true_u2 < 0.01d, true);

            var lh = norm.LogLikelihood(sample);
        }

        /// <summary>
        /// Test the quantile function for the Normal Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.1.2 page 91.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Normal_Quantile()
        {
            var N = new Normal(12665d, 4710d);
            double q100 = N.InverseCDF(0.99d);
            double true_q100 = 23624d;
            Assert.AreEqual((q100 - true_q100) / true_q100 < 0.01d, true);
            double p = N.CDF(q100);
            double true_p = 0.99d;
            Assert.AreEqual((p - true_p) / true_p < 0.01d, true);
        }

        /// <summary>
        /// Test the standard error for the Normal Distribution.
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
        public void Test_Normal_StandardError()
        {

            // Maximum Likelihood
            var N = new Normal(12665d, 4710d);
            double qVar99 = Math.Sqrt(N.QuantileVariance(0.99d, 48, ParameterEstimationMethod.MaximumLikelihood));
            double true_qVar99 = 1309d;
            Assert.AreEqual((qVar99 - true_qVar99) / true_qVar99 < 0.01d, true);
        }

        /// <summary>
        /// Verifying that input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var N = new Normal();
            Assert.AreEqual(N.Mu, 0);
            Assert.AreEqual(N.Sigma, 1);

            var N2 = new Normal(1, 1);
            Assert.AreEqual(N2.Mu, 1);
            Assert.AreEqual(N2.Sigma,1);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var N = new Normal(double.NaN, 1);
            Assert.IsFalse(N.ParametersValid);

            var N2 = new Normal(1,double.NaN);
            Assert.IsFalse(N2.ParametersValid);

            var N3 = new Normal(double.PositiveInfinity, 1);
            Assert.IsFalse(N3.ParametersValid);

            var N4 = new Normal(1,double.PositiveInfinity);
            Assert.IsFalse(N4.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var N = new Normal();
            Assert.AreEqual(N.ParametersToString[0, 0], "Mean (µ)");
            Assert.AreEqual(N.ParametersToString[1, 0], "Std Dev (σ)");
            Assert.AreEqual(N.ParametersToString[0, 1], "0");
            Assert.AreEqual(N.ParametersToString[1,1], "1");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Normal(5, 9);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var N = new Normal();
            Assert.AreEqual(N.Mean, 0);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(N2.Mean, 5);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var N = new Normal();
            Assert.AreEqual(N.Median, 0);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(N2.Median, 5);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var N = new Normal();
            Assert.AreEqual(N.Mode, 0);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(N2.Mode, 5);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var N = new Normal();
            Assert.AreEqual(N.StandardDeviation, 1);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(N2.StandardDeviation, 9);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var N = new Normal();
            Assert.AreEqual(N.Skewness, 0);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(N2.Skewness, 0);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var N = new Normal();
            Assert.AreEqual(N.Kurtosis, 3);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(N2.Kurtosis, 3);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var N = new Normal();
            Assert.AreEqual(N.Minimum, double.NegativeInfinity);
            Assert.AreEqual(N.Maximum, double.PositiveInfinity);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(N2.Minimum,double.NegativeInfinity);
            Assert.AreEqual(N2.Maximum,double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var N = new Normal();
            Assert.AreEqual(N.PDF(0), 0.39894, 1e-04);
            Assert.AreEqual(N.PDF(1), 0.24197, 1e-04);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(N2.PDF(-1), 0.03549, 1e-04);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var N = new Normal(5,2);
            Assert.AreEqual(N.CDF(0), 0.006209, 1e-04);
            Assert.AreEqual(N.CDF(4), 0.30853, 1e-04);
            Assert.AreEqual(N.CDF(5), 0.5);
        }
    }
}
