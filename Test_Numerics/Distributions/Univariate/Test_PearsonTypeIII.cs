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
    /// Testing the Pearson Type III distribution algorithm.
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
    public class Test_PearsonTypeIII
    {

        // Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee & F. Ashkar, Water Resources Publications, 1991.
        // Table 1.2 Maximum annual peak discharge values in cms, observed at the Harricana River at Amos (Quebec, Canada)
        private double[] sample = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

        /// <summary>
        /// Verification of PIII fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee & F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 6.3 page 70.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_P3_MOM()
        {
            var P3 = new PearsonTypeIII();
            P3.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double xi = P3.Xi;
            double beta = P3.Beta;
            double alpha = P3.Alpha;
            double mu = P3.Mu;
            double sigma = P3.Sigma;
            double gamma = P3.Gamma;
            double mean = P3.Mean;
            double stDev = P3.StandardDeviation;
            double skew = P3.Skewness;
            double true_xi = 79.84941d;
            double true_beta = 1d / 0.04846d;
            double true_alpha = 5.40148d;
            double true_mu = 191.31739d;
            double true_sigma = 47.96161d;
            double true_gamma = 0.86055d;
            double true_mean = 191.31739d;
            double true_stDev = 47.96161d;
            double true_skew = 0.86055d;
            Assert.AreEqual((xi - true_xi) / true_xi < 0.01d, true);
            Assert.AreEqual((beta - true_beta) / true_beta < 0.01d, true);
            Assert.AreEqual((alpha - true_alpha) / true_alpha < 0.01d, true);
            Assert.AreEqual((mu - true_mu) / true_mu < 0.01d, true);
            Assert.AreEqual((sigma - true_sigma) / true_sigma < 0.01d, true);
            Assert.AreEqual((gamma - true_gamma) / true_gamma < 0.01d, true);
            Assert.AreEqual((mean - true_mean) / true_mean < 0.01d, true);
            Assert.AreEqual((stDev - true_stDev) / true_stDev < 0.01d, true);
            Assert.AreEqual((skew - true_skew) / true_skew < 0.01d, true);
        }

        /// <summary>
        /// Verification of PIII fit with method of linear moments.
        /// </summary>
        [TestMethod()]
        public void Test_P3_LMOM_Fit()
        {
            var sample = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var P3 = new PearsonTypeIII();
            P3.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
            double x = P3.Xi;
            double a = P3.Alpha;
            double b = P3.Beta;
            double true_x = 863.4104d;
            double true_a = 10.02196d;
            double true_b = 78.36751d;
            Assert.AreEqual(x, true_x, 0.001d);
            Assert.AreEqual(a, true_a, 0.001d);
            Assert.AreEqual(b, true_b, 0.001d);
            var lmom = P3.LinearMomentsFromParameters(P3.GetParameters);
            Assert.AreEqual(lmom[0], 1648.806d, 0.001d);
            Assert.AreEqual(lmom[1], 138.2366d, 0.001d);
            Assert.AreEqual(lmom[2], 0.1033889d, 0.001d);
            Assert.AreEqual(lmom[3], 0.1258521d, 0.001d);
        }

        /// <summary>
        /// Verification of PIII fit with maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee & F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 6.1 page 64.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_P3_MLE()
        {
            var P3 = new PearsonTypeIII();
            P3.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double xi = P3.Xi;
            double beta = P3.Beta;
            double alpha = P3.Alpha;
            double mu = P3.Mu;
            double sigma = P3.Sigma;
            double gamma = P3.Gamma;
            double mean = P3.Mean;
            double stDev = P3.StandardDeviation;
            double skew = P3.Skewness;
            double true_xi = 39.38903d;
            double true_beta = 1d / 0.06872d;
            double true_alpha = 10.44062d;
            double true_mu = 191.31739d;
            double true_sigma = 47.01925d;
            double true_gamma = 0.61897d;
            double true_mean = 191.31739d;
            double true_stDev = 47.01925d;
            double true_skew = 0.61897d;
            Assert.AreEqual((xi - true_xi) / true_xi < 0.01d, true);
            Assert.AreEqual((beta - true_beta) / true_beta < 0.01d, true);
            Assert.AreEqual((alpha - true_alpha) / true_alpha < 0.01d, true);
            Assert.AreEqual((mu - true_mu) / true_mu < 0.01d, true);
            Assert.AreEqual((sigma - true_sigma) / true_sigma < 0.01d, true);
            Assert.AreEqual((gamma - true_gamma) / true_gamma < 0.01d, true);
            Assert.AreEqual((mean - true_mean) / true_mean < 0.01d, true);
            Assert.AreEqual((stDev - true_stDev) / true_stDev < 0.01d, true);
            Assert.AreEqual((skew - true_skew) / true_skew < 0.01d, true);
        }

        /// <summary>
        /// Test the quantile function for the Pearson Type III Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee & F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 6.1 page 64.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_P3_Quantile()
        {
            var P3 = new PearsonTypeIII(191.31739d, 47.01925d, -0.61897d);
            double q999 = P3.InverseCDF(0.99d);
            double true_q999 = 321.48d;
            Assert.AreEqual((q999 - true_q999) / true_q999 < 0.01d, true);
        }

        /// <summary>
        /// Test the standard error for the Pearson Type III Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee & F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 6.1 & 6.3 page 64-70.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_P3_StandardError()
        {

            // Method of Moments
            var P3 = new PearsonTypeIII(191.31739d, 47.96161d, 0.86055d);
            double qVar999 = Math.Sqrt(P3.QuantileVariance(0.99d, 69, ParameterEstimationMethod.MethodOfMoments));
            double true_qVar999 = 27.175d;
            Assert.AreEqual((qVar999 - true_qVar999) / true_qVar999 < 0.01d, true);

            // Maximum Likelihood
            P3 = new PearsonTypeIII(191.31739d, 47.01925d, 0.61897d);
            qVar999 = Math.Sqrt(P3.QuantileVariance(0.99d, 69, ParameterEstimationMethod.MaximumLikelihood));
            true_qVar999 = 20.045d;
            Assert.AreEqual((qVar999 - true_qVar999) / true_qVar999 < 0.01d, true);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.Mu, 100);
            Assert.AreEqual(P3.Sigma, 10);
            Assert.AreEqual(P3.Gamma, 0);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(P3ii.Mu, 1);
            Assert.AreEqual(P3ii.Sigma, 1);
            Assert.AreEqual(P3ii.Gamma, 1);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var P3 = new PearsonTypeIII(double.PositiveInfinity, double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(P3.ParametersValid);

            var P3ii = new PearsonTypeIII(double.NaN, double.NaN, double.NaN);
            Assert.IsFalse(P3ii.ParametersValid);

            var P3iii = new PearsonTypeIII(1, 0, 1);
            Assert.IsFalse(P3iii.ParametersValid);
        }

        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.ParametersToString[0, 0], "Mean (µ)");
            Assert.AreEqual(P3.ParametersToString[1, 0], "Std Dev (σ)");
            Assert.AreEqual(P3.ParametersToString[2, 0], "Skew (γ)");
            Assert.AreEqual(P3.ParametersToString[0, 1], "100");
            Assert.AreEqual(P3.ParametersToString[1, 1], "10");
            Assert.AreEqual(P3.ParametersToString[2, 1], "0");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new PearsonTypeIII();
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
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.Mean, 100);

            var P3ii = new PearsonTypeIII(100, 1, 1);
            Assert.AreEqual(P3ii.Mean, 100);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.Median, 100);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.Mode, 100);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(P3ii.Mode, 0.5);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.StandardDeviation, 10);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(P3ii.StandardDeviation, 1);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.Skewness, 0);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(P3.Skewness, 0);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.Kurtosis, 3);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(P3ii.Kurtosis, 4.5);
        }

        /// <summary>
        /// Testing minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.Minimum,double.NegativeInfinity);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(P3ii.Minimum, -1);
        }

        /// <summary>
        /// Testing maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum() 
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(P3.Maximum,double.PositiveInfinity);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(P3ii.Maximum, double.PositiveInfinity);
        }

    }


}
