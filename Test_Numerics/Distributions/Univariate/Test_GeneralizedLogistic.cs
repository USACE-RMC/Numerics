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
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Generalized Logistic distribution algorithm.
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
    public class Test_GeneralizedLogistic
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 9.2.1 East Fork White River at Seymour, IN
        private double[] sample = new double[] { 24000d, 7920d, 21900d, 47100d, 30400d, 36100d, 67100d, 7030d, 28200d, 40100d, 10300d, 11100d, 17000d, 65600d, 32600d, 36200d, 46400d, 3650d, 16800d, 44800d, 37100d, 42900d, 15000d, 33000d, 28000d, 78500d, 54000d, 28600d, 44000d, 13300d, 6120d, 11100d, 42100d, 33400d, 30100d, 32100d, 28100d, 59400d, 23800d, 52000d, 54900d, 25600d, 10900d, 33700d, 60200d, 39200d, 26300d, 27900d, 27000d, 22700d, 17500d, 46400d, 19300d, 12700d, 36000d, 39900d, 25400d, 30200d, 47000d, 39800d, 23800d, 29600d, 33400d, 15400d, 28400d, 26700d, 46500d, 61200d };

        /// <summary>
        /// Verification of Generalized Logistic Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 9.2.1 page 311.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GLO_MOM_Fit()
        {

            // The data is different than book, so I am using book moments.
            var GLO = new GeneralizedLogistic();
            // GLO.Estimate(sample, ParameterEstimationMethod.MethodOfMoments)
            GLO.SetParameters(GLO.DirectMethodOfMoments(new[] { 32714d, 16560.2848d, 0.49051d }));
            double x = GLO.Xi;
            double a = GLO.Alpha;
            double k = GLO.Kappa;
            double true_x = 31892d;
            double true_a = 9030d;
            double true_k = -0.05515d;
            Assert.AreEqual((x - true_x) / true_x < 0.01d, true);
            Assert.AreEqual((a - true_a) / true_a < 0.01d, true);
            Assert.AreEqual((k - true_k) / true_k < 0.01d, true);
        }

        /// <summary>
        /// Verification of Generalized Logistic Distribution fit with method of linear moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 9.2.1 page 311.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GLO_LMOM_Fit()
        {
            var sample = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var GLO = new GeneralizedLogistic();
            GLO.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
            double x = GLO.Xi;
            double a = GLO.Alpha;
            double k = GLO.Kappa;
            double true_x = 1625.42d;
            double true_a = 135.8186d;
            double true_k = -0.1033903d;
            Assert.AreEqual(x, true_x, 0.001d);
            Assert.AreEqual(a, true_a, 0.001d);
            Assert.AreEqual(k, true_k, 0.001d);
            var lmom = GLO.LinearMomentsFromParameters(GLO.GetParameters);
            Assert.AreEqual(lmom[0], 1648.806d, 0.001d);
            Assert.AreEqual(lmom[1], 138.2366d, 0.001d);
            Assert.AreEqual(lmom[2], 0.1033903d, 0.001d);
            Assert.AreEqual(lmom[3], 0.1755746d, 0.001d);
        }

        /// <summary>
        /// Verification of Generalized Logistic Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 9.1.1 page 313.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GLO_MLE_Fit()
        {
            var GLO = new GeneralizedLogistic();
            GLO.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double x = GLO.Xi;
            double a = GLO.Alpha;
            double k = GLO.Kappa;
            double true_x = 30911.83d;
            double true_a = 9305.0205d;
            double true_k = -0.144152d;
            Assert.AreEqual((x - true_x) / true_x < 0.01d, true);
            Assert.AreEqual((a - true_a) / true_a < 0.01d, true);
            Assert.AreEqual((k - true_k) / true_k < 0.01d, true);
        }

        /// <summary>
        /// Test the quantile function for the Generalized Logistic Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 9.2.2 page 315.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GLO_Quantile()
        {
            var GLO = new GeneralizedLogistic(31892d, 9030d, -0.05515d);
            double q100 = GLO.InverseCDF(0.99d);
            double true_100 = 79117d;
            Assert.AreEqual((q100 - true_100) / true_100 < 0.01d, true);
            double p = GLO.CDF(q100);
            double true_p = 0.99d;
            Assert.AreEqual(p == true_p, true);
        }

        /// <summary>
        /// Test the partial derivatives for the Generalized Logistic Distribution.
        /// </summary>
        [TestMethod()]
        public void Test_GLO_Partials()
        {
            var GLO = new GeneralizedLogistic(30911.83d, 9305.0205d, -0.144152d);
            double dQdLocation = GLO.QuantileGradient(0.99d)[0];
            double dQdScale = GLO.QuantileGradient(0.99d)[1];
            double dQdShape = GLO.QuantileGradient(0.99d)[2];
            double true_dLocation = 1.0d;
            double true_dScale = 6.51695d;
            double true_dShape = -154595.08d;
            Assert.AreEqual((dQdLocation - true_dLocation) / true_dLocation < 0.01d, true);
            Assert.AreEqual((dQdScale - true_dScale) / true_dScale < 0.01d, true);
            Assert.AreEqual((dQdShape - true_dShape) / true_dShape < 0.01d, true);
        }

        /// <summary>
        /// Testing parameters can create a generalized logistic distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.Xi, 100);
            Assert.AreEqual(l.Alpha, 10);
            Assert.AreEqual(l.Kappa, 0);

            var l2 = new GeneralizedLogistic(-100, 10, 1);
            Assert.AreEqual(l2.Xi, -100);
            Assert.AreEqual(l2.Alpha, 10);
            Assert.AreEqual(l2.Kappa, 1);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var l = new GeneralizedLogistic(double.NaN,double.NaN,double.NaN);
            Assert.IsFalse(l.ParametersValid);

            var l2 = new GeneralizedLogistic(double.PositiveInfinity,double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(l2.ParametersValid);

            var l3 = new GeneralizedLogistic(100, 0, 1);
            Assert.IsFalse(l3.ParametersValid);
        }

        /// <summary>
        /// Testing ParameterToString()
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.ParametersToString[0, 0], "Location (ξ)");
            Assert.AreEqual(l.ParametersToString[1, 0], "Scale (α)");
            Assert.AreEqual(l.ParametersToString[2, 0], "Shape (κ)");
            Assert.AreEqual(l.ParametersToString[0, 1], "100");
            Assert.AreEqual(l.ParametersToString[1, 1], "10");
            Assert.AreEqual(l.ParametersToString[2, 1], "0");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new GeneralizedLogistic(100, 10, -0.1);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Checking mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.Mean, 100);

            var l2 = new GeneralizedLogistic(100, 10, 0.9);
            Assert.AreEqual(l2.Mean, 9.44703, 1e-04);

            var l3 = new GeneralizedLogistic(100, 10, 1);
            Assert.AreEqual(l3.Mean, double.NaN);
        }

        /// <summary>
        /// Testing Median is Xi.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.Median, 100);

            var l2 = new GeneralizedLogistic(10, 10, 1);
            Assert.AreEqual(l2.Median, 10);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.Mode, 100);

            var l2 = new GeneralizedLogistic(100, 10, 1);
            Assert.AreEqual(l2.Mode, 95);
        }
        
        /// <summary>
        /// Checking Standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.StandardDeviation, 18.13799, 1e-04);

            var l2 = new GeneralizedLogistic(100, 10, 0.4);
            Assert.AreEqual(l2.StandardDeviation, 39.76482, 1e-04);

            var l3 = new GeneralizedLogistic(100, 10, 1);
            Assert.AreEqual(l3.StandardDeviation, double.NaN);
        }

        /// <summary>
        /// Checking skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.Skewness, 0);

            var l2 = new GeneralizedLogistic(100, 10, 0.3);
            Assert.AreEqual(l2.Skewness, -10.90354, 1e-04);

            var l3 = new GeneralizedLogistic(100, 10, 1);
            Assert.AreEqual(l3.Skewness, double.NaN);
        }

        /// <summary>
        /// Testing Kurtosis method.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.Kurtosis, 21d / 5d);

            var l2 = new GeneralizedLogistic(100, 10, 0.24);
            Assert.AreEqual(l2.Kurtosis, 199.733369,1e-04);

            var l3 = new GeneralizedLogistic(100, 10, 0.25);
            Assert.AreEqual(l3.Kurtosis, double.NaN);
        }

        /// <summary>
        /// Testing minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.Minimum, double.NegativeInfinity);

            var l2 = new GeneralizedLogistic(100, 10, -5);
            Assert.AreEqual(l2.Minimum, 98);
        }

        /// <summary>
        /// Testing maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum() 
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.Maximum, double.PositiveInfinity);

            var l2 = new GeneralizedLogistic(100, 10, 1);
            Assert.AreEqual(l2.Maximum, 110);
        }

        /// <summary>
        /// Testing PDF method
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.PDF(100), 0.025);
            Assert.AreEqual(l.PDF(0), 4.5395e-06,1e-10);

            var l2 = new GeneralizedLogistic(100, 10, 1);
            Assert.AreEqual(l2.PDF(100), 0.025);
            Assert.AreEqual(l2.PDF(0),6.9444e-04,1e-08);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.CDF(100), 0.5);
            Assert.AreEqual(l.CDF(0), 4.5397e-05, 1e-8);

            var l2 = new GeneralizedLogistic(100, 10, 1);
            Assert.AreEqual(l2.CDF(0), 0.083333, 1e-04);
        }

        /// <summary>
        /// Testing the inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var l = new GeneralizedLogistic();
            Assert.AreEqual(l.InverseCDF(0), double.NegativeInfinity);
            Assert.AreEqual(l.InverseCDF(0.5), 100);
            Assert.AreEqual(l.InverseCDF(1),double.PositiveInfinity);

            var l2 = new GeneralizedLogistic(100, 10, 1);
            Assert.AreEqual(l2.InverseCDF(0.5), 100);
            Assert.AreEqual(l2.InverseCDF(0.7), 105.714285,1e-04);
        }
    }
}
