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
    /// Testing the Exponential distribution algorithm.
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
    public class Test_Exponential
    {
        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 1.8.1 Wabash River at Lafayette, IN
        private double[] sample = new double[] { 41500d, 57000d, 44000d, 49000d, 31000d, 45900d, 19000d, 41100d, 37300d, 76000d, 33200d, 61200d, 76000d, 59800d, 44400d, 58400d, 53600d, 59800d, 63300d, 57700d, 64000d, 63500d, 38000d, 74600d, 13100d, 37600d, 67500d, 21700d, 37000d, 93500d, 58500d, 63300d, 74400d, 34200d, 14600d, 44200d, 13100d, 73300d, 46600d, 39400d, 41200d, 41300d, 62000d, 90000d, 50600d, 41900d, 35000d, 16500d, 35300d, 30000d, 52600d, 99000d, 89000d, 39500d, 55400d, 46000d, 63000d, 58300d, 36500d, 14600d, 64900d, 68500d, 69100d, 42600d, 31000d, 39400d, 40700d, 53400d, 36000d, 43900d, 23600d, 50500d, 49700d, 48100d, 44500d, 56400d, 60800d, 40400d, 80400d, 41600d, 14700d, 33300d, 40700d, 53300d, 77400d };

        /// <summary>
        /// Verification of Exponential Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 6.1.1 page 132.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_EXP_MOM_Fit()
        {
            var EXP = new Exponential();
            EXP.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double x = EXP.Xi;
            double a = EXP.Alpha;
            double true_x = 30314.48d;
            double true_a = 18907.87d;
            Assert.AreEqual((x - true_x) / true_x < 0.01d, true);
            Assert.AreEqual((a - true_a) / true_a < 0.01d, true);
        }

        /// <summary>
        /// Test the fit against the R 'lmom' package.
        /// </summary>
        [TestMethod()]
        public void Test_EXP_LMOM_Fit()
        {
            var sample = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var EXP = new Exponential();
            EXP.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
            double x = EXP.Xi;
            double a = EXP.Alpha;
            double true_x = 1372.333d;
            double true_a = 276.4731d;
            Assert.AreEqual(x, true_x, 0.001d);
            Assert.AreEqual(a, true_a, 0.001d);
            var lmom = EXP.LinearMomentsFromParameters(EXP.GetParameters);
            Assert.AreEqual(lmom[0], 1648.806d, 0.001d);
            Assert.AreEqual(lmom[1], 138.2366d, 0.001d);
            Assert.AreEqual(lmom[2], 0.3333333d, 0.001d);
            Assert.AreEqual(lmom[3], 0.1666667d, 0.001d);
        }

        /// <summary>
        /// Verification of Exponential Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 6.1.1 page 132.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_EXP_MLE_Fit()
        {
            var EXP = new Exponential();
            EXP.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double x = EXP.Xi;
            double a = EXP.Alpha;
            double true_x = 13100d;
            double true_a = 36122.35d;
            Assert.AreEqual((x - true_x) / true_x < 0.01d, true);
            Assert.AreEqual((a - true_a) / true_a < 0.01d, true);
        }

        /// <summary>
        /// Test the quantile function for the Exponential Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 6.1.2 page 134.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_EXP_Quantile()
        {
            var EXP = new Exponential(27421d, 25200d);
            double q100 = EXP.InverseCDF(0.99d);
            double true_100 = 143471d;
            Assert.AreEqual((q100 - true_100) / true_100 < 0.01d, true);
            double p = EXP.CDF(q100);
            double true_p = 0.99d;
            Assert.AreEqual(p == true_p, true);
        }

        /// <summary>
        /// Test the standard error for the Exponential Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 6.1.3 page 138.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_EXP_StandardError()
        {

            // Method of Moments
            var EXP = new Exponential(27421d, 25200d);
            double se100 = Math.Sqrt(EXP.QuantileVariance(0.99d, 85, ParameterEstimationMethod.MethodOfMoments));
            double true_se100 = 15986d;
            Assert.AreEqual((se100 - true_se100) / true_se100 < 0.01d, true);

            // Maximum Likelihood
            EXP = new Exponential(12629d, 39991d);
            se100 = Math.Sqrt(EXP.QuantileVariance(0.99d, 85, ParameterEstimationMethod.MaximumLikelihood));
            true_se100 = 20048d;
            Assert.AreEqual((se100 - true_se100) / true_se100 < 0.01d, true);
        }

        /// <summary>
        /// Test the partial derivatives for the Exponential Distribution.
        /// </summary>
        [TestMethod()]
        public void Test_EXP_Partials()
        {
            var EXP = new Exponential(30314.48d, 18907.87d);
            double dQdLocation = EXP.QuantileGradient(0.99d)[0];
            double dQdScale = EXP.QuantileGradient(0.99d)[1];
            double true_dLocation = 1.0d;
            double true_dScale = 4.60517d;
            Assert.AreEqual((dQdLocation - true_dLocation) / true_dLocation < 0.01d, true);
            Assert.AreEqual((dQdScale - true_dScale) / true_dScale < 0.01d, true);
        }

        /// <summary>
        /// Validating parameters can create exponential function
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var EXP = new Exponential(-5, 100);
            Assert.AreEqual(EXP.Xi, -5);
            Assert.AreEqual(EXP.Alpha, 100);

            var EXP2 = new Exponential(0, 1);
            Assert.AreEqual(EXP2.Xi, 0);
            Assert.AreEqual(EXP2.Alpha, 1);
        }

        /// <summary>
        /// Testing Exponential distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var EXP = new Exponential(double.NaN, 5);
            Assert.IsFalse(EXP.ParametersValid);

            var EXP2 = new Exponential(double.PositiveInfinity, 100);
            Assert.IsFalse(EXP2.ParametersValid);

            var EXP3 = new Exponential(0, double.NegativeInfinity);
            Assert.IsFalse(EXP3.ParametersValid);

            var EXP4 = new Exponential(1,double.NaN);
            Assert.IsFalse(EXP4.ParametersValid);
        }

        /// <summary>
        /// Checking parameters to string function.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var EXP = new Exponential(1, 1);
            Assert.AreEqual(EXP.ParametersToString[0, 0], "Location (ξ)");
            Assert.AreEqual(EXP.ParametersToString[1, 0], "Scale (α)");
            Assert.AreEqual(EXP.ParametersToString[0, 1], "1");
            Assert.AreEqual(EXP.ParametersToString[1, 1], "1");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Exponential(1, 1);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Checking Mean of exponential distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var EXP = new Exponential(1, 1);
            Assert.AreEqual(EXP.Mean, 2);

            var EXP2 = new Exponential(-100, 4);
            Assert.AreEqual(EXP2.Mean, -96);
        }

        /// <summary>
        /// Checking median of distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var EXP = new Exponential(0,1);
            Assert.AreEqual(EXP.Median, 0.693147, 1e-04);

            var EXP2 = new Exponential(-100, 1);
            Assert.AreEqual(EXP2.Median, -99.306852, 1e-04);
        }

        /// <summary>
        /// Checking mode of Exponential is equal to the location parameter.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var EXP = new Exponential(0,1);
            Assert.AreEqual(EXP.Mode, 0);

            var EXP2 = new Exponential(-100,1);
            Assert.AreEqual(EXP2.Mode, -100);
        }

        /// <summary>
        /// Checking Standard deviation is equal to the scale parameter.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var EXP = new Exponential(0,1);
            Assert.AreEqual(EXP.StandardDeviation, 1);

            var EXP2 = new Exponential(-100, 1);
            Assert.AreEqual(EXP2.StandardDeviation, 1);
        }

        /// <summary>
        /// Checking skew is equal to 2 regardless of parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var EXP = new Exponential(0, 1);
            Assert.AreEqual(EXP.Skewness, 2);
        }

        /// <summary>
        /// Checking Kurtosis is equal to 9 regardless of parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var EXP = new Exponential(0, 1);
            Assert.AreEqual(EXP.Kurtosis, 9);
        }

        /// <summary>
        /// Checking minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var EXP = new Exponential(0, 1);
            Assert.AreEqual(EXP.Minimum, 0);

            var EXP2 = new Exponential(-100, 1);
            Assert.AreEqual(EXP2.Minimum, -100);
        }

        /// <summary>
        /// Checking maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var EXP = new Exponential(0, 1);
            Assert.AreEqual(EXP.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Checking PDF function with different points in the distribution range.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var EXP = new Exponential(0, 3);
            Assert.AreEqual(0.33333, EXP.PDF(0),1e-04);
            Assert.AreEqual(0.23884377, EXP.PDF(1),1e-04);
            Assert.AreEqual(0.171139, EXP.PDF(2),1e-04);
            Assert.AreEqual(0.122626, EXP.PDF(3),1e-04);
        }

        /// <summary>
        /// Checking CDF function with different points in the distribution range.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var EXP = new Exponential(0, 3);
            Assert.AreEqual(0, EXP.CDF(0));
            Assert.AreEqual(0.2834686, EXP.CDF(1), 1e-04);
            Assert.AreEqual(0.4865828, EXP.CDF(2), 1e-04);
            Assert.AreEqual(0.632120, EXP.CDF(3), 1e-04);
        }

        /// <summary>
        /// Test the inverse CDF.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var EXP = new Exponential(0, 1);
            Assert.AreEqual(0,EXP.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity, EXP.InverseCDF(1));
            Assert.AreEqual(0.693147, EXP.InverseCDF(0.5), 1e-04);

            var EXP2 = new Exponential();
            Assert.AreEqual(103.5667494, EXP2.InverseCDF(0.3), 1e-04);
            Assert.AreEqual(106.93147, EXP2.InverseCDF(0.5), 1e-04);
            Assert.AreEqual(112.039728, EXP2.InverseCDF(0.7), 1e-04);
        }
    }
}
