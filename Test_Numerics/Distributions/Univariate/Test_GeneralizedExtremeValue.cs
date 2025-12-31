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
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Generalized Extreme Value distribution algorithm.
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
    public class Test_GeneralizedExtremeValue
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 7.1.2 White River near Nora, IN
        private double[] sample = new double[] { 23200d, 2950d, 10300d, 23200d, 4540d, 9960d, 10800d, 26900d, 23300d, 20400d, 8480d, 3150d, 9380d, 32400d, 20800d, 11100d, 7270d, 9600d, 14600d, 14300d, 22500d, 14700d, 12700d, 9740d, 3050d, 8830d, 12000d, 30400d, 27000d, 15200d, 8040d, 11700d, 20300d, 22700d, 30400d, 9180d, 4870d, 14700d, 12800d, 13700d, 7960d, 9830d, 12500d, 10700d, 13200d, 14700d, 14300d, 4050d, 14600d, 14400d, 19200d, 7160d, 12100d, 8650d, 10600d, 24500d, 14400d, 6300d, 9560d, 15800d, 14300d, 28700d };

        /// <summary>
        /// Verification of GEV Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.1 page 218.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_MOM_Fit()
        {
            var GEV = new GeneralizedExtremeValue();
            GEV.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double x = GEV.Xi;
            double a = GEV.Alpha;
            double k = GEV.Kappa;
            double true_x = 11012d;
            double true_a = 6209.4d;
            double true_k = 0.0736d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x );
            Assert.IsLessThan(0.01d, (a - true_a) / true_a );
            Assert.IsLessThan(0.01d, (k - true_k) / true_k);
        }

        /// <summary>
        /// Verification of GEV Distribution fit with method of linear moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.1 page 218.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_LMOM_Fit()
        {
            var sample = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var GEV = new GeneralizedExtremeValue();
            GEV.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
            double x = GEV.Xi;
            double a = GEV.Alpha;
            double k = GEV.Kappa;
            double true_x = 1543.933d;
            double true_a = 218.1148d;
            double true_k = 0.1068473d;
            Assert.AreEqual(x, true_x, 0.001d);
            Assert.AreEqual(a, true_a, 0.001d);
            Assert.AreEqual(k, true_k, 0.001d);
            var lmom = GEV.LinearMomentsFromParameters(GEV.GetParameters);
            Assert.AreEqual(1648.806d, lmom[0],  0.001d);
            Assert.AreEqual(138.2366d, lmom[1],  0.001d);
            Assert.AreEqual(0.1030703d, lmom[2],  0.001d);
            Assert.AreEqual(0.1277244d, lmom[3],  0.001d);
        }

        /// <summary>
        /// Verification of GEV Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.1 page 219.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_MLE_Fit()
        {
            var GEV = new GeneralizedExtremeValue();
            GEV.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double x = GEV.Xi;
            double a = GEV.Alpha;
            double k = GEV.Kappa;
            double true_x = 10849d;
            double true_a = 5745.6d;
            double true_k = 0.005d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x);
            Assert.IsLessThan(0.01d, (a - true_a) / true_a);
            Assert.IsLessThan(0.01d, (k - true_k) / true_k);
        }

        /// <summary>
        /// Test the quantile function for the GEV Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.2 page 221.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_Quantile()
        {
            var GEV = new GeneralizedExtremeValue(10849d, 5745.6d, 0.005d);
            double q100 = GEV.InverseCDF(0.99d);
            double true_q100 = 36977d;
            Assert.IsLessThan(0.01d, (q100 - true_q100) / true_q100);
            double p = GEV.CDF(q100);
            double true_p = 0.99d;
            Assert.IsLessThan(0.01d, (p - true_p) / true_p);
        }

        /// <summary>
        /// Test the standard error for the GEV Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.3 page 226.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_StandardError()
        {
            double u = 10849d;
            double a = 5745.6d;
            double k = 0.005d;
            double true_dXdU = 1.0d;
            double true_dxdA = 4.5472d;
            double true_dxdK = -59861;
            double true_VarU = 664669d;
            double true_VarA = 346400d;
            double true_VarK = 0.007655d;
            double true_CovarUA = 176180d;
            double true_CovarUK = 23.977d;
            double true_CovarAK = 13.8574d;
            double true_QVar = 26445364d;
            double true_QSigma = 5142d;
            var GEV = new GeneralizedExtremeValue(u, a, k);
            var partials = GEV.QuantileGradient(0.99d).ToArray();
            var covar = GEV.ParameterCovariance(sample.Length, ParameterEstimationMethod.MaximumLikelihood);
            double qVar = GEV.QuantileVariance(0.99d, sample.Length, ParameterEstimationMethod.MaximumLikelihood);
            double qSigma = Math.Sqrt(qVar);
            Assert.IsLessThan(0.01d,(partials[0] - true_dXdU) / true_dXdU);
            Assert.IsLessThan(0.01d, (partials[1] - true_dxdA) / true_dxdA);
            Assert.IsLessThan(0.01d, (partials[2] - true_dxdK) / true_dxdK);
            Assert.IsLessThan(0.01d, (covar[0, 0] - true_VarU) / true_VarU);
            Assert.IsLessThan(0.01d, (covar[1, 1] - true_VarA) / true_VarA);
            Assert.IsLessThan(0.01d, (covar[2, 2] - true_VarK) / true_VarK);
            Assert.IsLessThan(0.01d, (covar[0, 1] - true_CovarUA) / true_CovarUA);
            Assert.IsLessThan(0.01d, (covar[0, 2] - true_CovarUK) / true_CovarUK);
            Assert.IsLessThan(0.01d, (covar[1, 2] - true_CovarAK) / true_CovarAK);
            Assert.IsLessThan(0.01d, (qVar - true_QVar) / true_QVar);
            Assert.IsLessThan(0.01d, (qSigma - true_QSigma) / true_QSigma);
        }

        /// <summary>
        /// Testing GEV constructor
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(100,GEV.Xi);
            Assert.AreEqual(10,GEV.Alpha);
            Assert.AreEqual(0,GEV.Kappa);

            var GEV2 = new GeneralizedExtremeValue(-100, 1, 1);
            Assert.AreEqual(-100,GEV2.Xi);
            Assert.AreEqual(1, GEV2.Alpha);
            Assert.AreEqual(1, GEV2.Kappa);
        }

        /// <summary>
        /// Testing bad parameters on GEV distribution.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var GEV = new GeneralizedExtremeValue(double.NaN, double.NaN, double.NaN);
            Assert.IsFalse(GEV.ParametersValid);

            var GEV2 = new GeneralizedExtremeValue(double.PositiveInfinity,double.PositiveInfinity, double.PositiveInfinity);
            Assert.IsFalse(GEV2.ParametersValid);

            var GEV3 = new GeneralizedExtremeValue(100, 0, 1);
            Assert.IsFalse(GEV3.ParametersValid);
        }

        /// <summary>
        /// Testing Parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual("Location (ξ)", GEV.ParametersToString[0, 0]);
            Assert.AreEqual("Scale (α)", GEV.ParametersToString[1, 0]);
            Assert.AreEqual("Shape (κ)", GEV.ParametersToString[2, 0]);
            Assert.AreEqual("100", GEV.ParametersToString[0, 1]);
            Assert.AreEqual("10", GEV.ParametersToString[1, 1]);
            Assert.AreEqual("0", GEV.ParametersToString[2, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new GeneralizedExtremeValue(100, 10, -0.1);
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
            var GEV = new GeneralizedExtremeValue();
            var true_val = 100 + 10 * Tools.Euler;
            Assert.AreEqual(true_val, GEV.Mean);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.9);
            Assert.AreEqual(100.42482, GEV2.Mean, 1e-04);

            var GEV3 = new GeneralizedExtremeValue(100, 10, 10);
            Assert.AreEqual(double.NaN, GEV3.Mean);
        }

        /// <summary>
        /// Testing Median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(103.66512, GEV.Median,  1e-04);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.9);
            Assert.AreEqual(104.3419519, GEV2.Median,  1e-04);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(100, GEV.Mode);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(95, GEV2.Mode);
        }

        /// <summary>
        /// Testing standard deviation method.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(12.825498, GEV.StandardDeviation,  1e-05);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.49);
            Assert.AreEqual(9.280898, GEV2.StandardDeviation, 1e-04);

            var GEV3 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(double.NaN, GEV3.StandardDeviation);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(1.1396, GEV.Skewness);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.3);
            Assert.AreEqual(-0.0690175, GEV2.Skewness,  1e-03);

            var GEV3 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(double.NaN, GEV3.Skewness);
        }

        /// <summary>
        /// Testing kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(3 + 12d / 5d, GEV.Kurtosis);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.24);
            Assert.AreEqual(2.7659607, GEV2.Kurtosis,  1e-04);

            var GEV3 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(double.NaN, GEV3.Kurtosis);
        }

        /// <summary>
        /// Testing minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(double.NegativeInfinity, GEV.Minimum);

            var GEV2 = new GeneralizedExtremeValue(100, 10, -5);
            Assert.AreEqual(98, GEV2.Minimum);
        }

        /// <summary>
        /// Testing maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(double.PositiveInfinity, GEV.Maximum);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(110, GEV2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(0,GEV.PDF(0));
            Assert.AreEqual(0,GEV.PDF(1));

            var GEV2 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(1.67017007902456E-06, GEV2.PDF(0), 1e-10);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(0.367879, GEV.CDF(100),  1e-04);
            Assert.AreEqual(0.9999546, GEV.CDF(200),  1e-07);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(0.367879, GEV2.CDF(100),  1e-05);
            Assert.AreEqual(1, GEV2.CDF(200));
        }

        /// <summary>
        /// Testing InverseCDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(double.NegativeInfinity, GEV.InverseCDF(0));
            Assert.AreEqual(103.66512, GEV.InverseCDF(0.5),  1e-05);
            Assert.AreEqual(double.PositiveInfinity,GEV.InverseCDF(1));
        }
    }
}
