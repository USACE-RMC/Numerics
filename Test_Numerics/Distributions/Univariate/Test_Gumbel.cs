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
    /// Testing the Gumbel distribution algorithm.
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
    public class Test_Gumbel
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 7.2.1 Sugar Creek at Crawfordsville, IN
        private double[] sample = new double[] { 17600d, 3660d, 903d, 5050d, 24000d, 11400d, 9470d, 8970d, 7710d, 14800d, 13900d, 20800d, 9470d, 7860d, 7860d, 2730d, 6480d, 18200d, 26300d, 15100d, 14600d, 7300d, 8580d, 15100d, 15100d, 21800d, 6200d, 2130d, 11100d, 14300d, 11200d, 6670d, 5440d, 9370d, 6900d, 9680d, 6810d, 7730d, 5290d, 12200d, 9750d, 7390d, 13100d, 7190d, 8850d, 6290d, 18800d, 9740d, 2990d, 6950d, 9390d, 12400d, 21200d };

        /// <summary>
        /// Verification of Gumbel Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.2.1 page 234.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GUM_MOM_Fit()
        {
            var GUM = new Gumbel();
            GUM.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double x = GUM.Xi;
            double a = GUM.Alpha;
            double true_x = 8074.4d;
            double true_a = 4441.4d;
            Assert.AreEqual((x - true_x) / true_x < 0.01d, true);
            Assert.AreEqual((a - true_a) / true_a < 0.01d, true);
        }

        /// <summary>
        /// Verification of Gumbel Distribution fit with method of linear moments.
        /// </summary>
        [TestMethod()]
        public void Test_GUM_LMOM_Fit()
        {
            var sample = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var GUM = new Gumbel();
            GUM.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
            double x = GUM.Xi;
            double a = GUM.Alpha;
            double true_x = 1533.69d;
            double true_a = 199.4332d;
            Assert.AreEqual(x, true_x, 0.001d);
            Assert.AreEqual(a, true_a, 0.001d);
            var lmom = GUM.LinearMomentsFromParameters(GUM.GetParameters);
            Assert.AreEqual(lmom[0], 1648.806d, 0.001d);
            Assert.AreEqual(lmom[1], 138.2366d, 0.001d);
            Assert.AreEqual(lmom[2], 0.169925d, 0.001d);
            Assert.AreEqual(lmom[3], 0.150375d, 0.001d);
        }

        /// <summary>
        /// Verification of Gumbel Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.2.1 page 234.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GUM_MLE_Fit()
        {
            var GUM = new Gumbel();
            GUM.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double x = GUM.Xi;
            double a = GUM.Alpha;
            double true_x = 8049.6d;
            double true_a = 4478.6d;
            Assert.AreEqual((x - true_x) / true_x < 0.01d, true);
            Assert.AreEqual((a - true_a) / true_a < 0.01d, true);
        }

        /// <summary>
        /// Test the quantile function for the Gumbel Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.7.2 page 237.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Gumbel_Quantile()
        {
            var GUM = new Gumbel(8049.6d, 4478.6d);
            double q100 = GUM.InverseCDF(0.99d);
            double true_q100 = 28652d;
            Assert.AreEqual((q100 - true_q100) / true_q100 < 0.01d, true);
            double p = GUM.CDF(q100);
            double true_p = 0.99d;
            Assert.AreEqual((p - true_p) / true_p < 0.01d, true);

        }

        /// <summary>
        /// Test the standard error for the Gumbel Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.2.3 page 240.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Gumbel_StandardError()
        {

            // Maximum Likelihood
            var GUM = new Gumbel(8049.6d, 4478.6d);
            double qVar99 = Math.Sqrt(GUM.QuantileVariance(0.99d, 53, ParameterEstimationMethod.MaximumLikelihood));
            double true_qVar99 = 2486.5d;
            Assert.AreEqual((qVar99 - true_qVar99) / true_qVar99 < 0.01d, true);
        }

        /// <summary>
        /// Checking Gumbel is created with inputs.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.Xi, 100);
            Assert.AreEqual(GUM.Alpha, 10);

            var GUM2 = new Gumbel(-100, 1);
            Assert.AreEqual(GUM2.Xi, -100);
            Assert.AreEqual(GUM2.Alpha, 1);
        }

        /// <summary>
        /// Testing Gumbel with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var GUM = new Gumbel(double.NaN,double.NaN);
            Assert.IsFalse(GUM.ParametersValid);

            var GUM2 = new Gumbel(double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(GUM2.ParametersValid);

            var GUM3 = new Gumbel(100, 0);
            Assert.IsFalse(GUM3.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString()
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.ParametersToString[0, 0], "Location (ξ)");
            Assert.AreEqual(GUM.ParametersToString[1, 0], "Scale (α)");
            Assert.AreEqual(GUM.ParametersToString[0, 1], "100");
            Assert.AreEqual(GUM.ParametersToString[1, 1], "10");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Gumbel(10, 1);
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
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.Mean, 105.77215, 1e-04);

            var GUM2 = new Gumbel(10, 1);
            Assert.AreEqual(GUM2.Mean, 10.577215, 1e-04);
        }

        /// <summary>
        /// Testing median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.Median, 103.66512, 1e-05);

            var GUM2 = new Gumbel(10, 1);
            Assert.AreEqual(GUM2.Median, 10.366512, 1e-04);
        }

        /// <summary>
        /// Testing Standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.StandardDeviation, 12.82549, 1e-04);

            var GUM2 = new Gumbel(10, 1);
            Assert.AreEqual(GUM2.StandardDeviation, 1.28254, 1e-04);
        }

        /// <summary>
        /// Testing skew is 1.1396.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.Skewness, 1.1396);

            var GUM2 = new Gumbel(10, 1);
            Assert.AreEqual(GUM2.Skewness, 1.1396);
        }

        /// <summary>
        /// Testing kurtosis 5.4.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.Kurtosis, 5.4);

            var GUM2 = new Gumbel(10, 1);
            Assert.AreEqual(GUM2.Kurtosis, 5.4);
        }

        /// <summary>
        /// Testing minimum and maximum of distribution.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.Minimum, double.NegativeInfinity);
            Assert.AreEqual(GUM.Maximum,double.PositiveInfinity);

            var GUM2 = new Gumbel(10, 1);
            Assert.AreEqual(GUM2.Minimum, double.NegativeInfinity);
            Assert.AreEqual(GUM2.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method at different locations and with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.PDF(100), 0.0367879, 1e-04);
            Assert.AreEqual(GUM.PDF(0), 0);
            Assert.AreEqual(GUM.PDF(200), 4.5397e-06, 1e-10);

            var GUM2 = new Gumbel(10, 1);
            Assert.AreEqual(GUM2.PDF(17), 9.1105e-04, 1e-09);
        }

        /// <summary>
        /// Testing CDF method at different locations and with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.CDF(100), 0.36787, 1e-04);
            Assert.AreEqual(GUM.CDF(50), 3.5073e-65, 1e-68);
            Assert.AreEqual(GUM.CDF(-10), 0);

            var GUM2 = new Gumbel(10, 2);
            Assert.AreEqual(GUM2.CDF(5), 5.11929e-06, 1e-10);
        }

        /// <summary>
        /// Testing InverseCDF with different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var GUM = new Gumbel();
            Assert.AreEqual(GUM.InverseCDF(0), double.NegativeInfinity);
            Assert.AreEqual(GUM.InverseCDF(1), double.PositiveInfinity);
            Assert.AreEqual(GUM.InverseCDF(0.3), 98.14373, 1e-04);
            Assert.AreEqual(GUM.InverseCDF(0.7), 110.309304, 1e-04);
        }
    }
}
