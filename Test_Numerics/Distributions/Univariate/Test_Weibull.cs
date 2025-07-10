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
using Numerics.Data.Statistics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Weibull distribution algorithm.
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
    public class Test_Weibull
    {

        private double[] sample = new double[] { 5.85217239831041d, 11.7689913217945d, 7.92234431846412d, 2.38759244506478d, 15.5436696441499d, 6.08101380255694d, 2.86541835011654d, 14.3272883381316d, 15.2293040674263d, 9.97119823770777d, 15.0078313002315d, 7.61635445751015d, 15.1579888433448d, 1.10494899327883d, 5.16175715794861d, 1.57293533830851d, 2.63953242550359d, 10.2188413398386d, 17.9577682621499d, 12.4224983994072d, 11.1290575697936d, 4.58559057415731d, 8.36394145136537d, 7.73020853953012d, 4.3003409576186d, 13.6952916728348d, 5.66874359549936d, 3.90607944288712d, 6.87181072557784d, 6.39556271370504d, 8.34155604676012d, 1.29993054120872d, 8.10693236597578d, 1.69364361593805d, 11.2322364035615d, 10.5062782641941d, 10.6177306237836d, 8.22440512843937d, 4.92470777947916d, 10.9335442619582d, 9.14151844388434d, 10.4858916284884d, 18.6061687709787d, 11.6935338563956d, 13.9506199815183d, 4.17482853658164d, 14.7373704498786d, 2.36384889353382d, 5.44157070547486d, 5.96907635426314d };

        /// <summary>
        /// Verification of Weibull Distribution fit with method of maximum likelihood.
        /// </summary>
        [TestMethod()]
        public void Test_Weibull_MLE_Fit()
        {
            var W = new Weibull();
            W.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double lamda = W.Lambda;
            double kappa = W.Kappa;
            double true_L = 9.589d;
            double true_k = 1.907d;
            Assert.AreEqual((lamda - true_L) / true_L < 0.01d, true);
            Assert.AreEqual((kappa - true_k) / true_k < 0.01d, true);
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
        public void Test_Weibull_Quantile()
        {
            var W = new Weibull(9.589d, 1.907d);
            double q100 = W.InverseCDF(0.99d);
            double true_q100 = 21.358d;
            Assert.AreEqual((q100 - true_q100) / true_q100 < 0.01d, true);
            double p = W.CDF(q100);
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
        public void Test_Weibull_StandardError()
        {

            // Maximum Likelihood
            var GUM = new Gumbel(8049.6d, 4478.6d);
            double qVar99 = Math.Sqrt(GUM.QuantileVariance(0.99d, 53, ParameterEstimationMethod.MaximumLikelihood));
            double true_qVar99 = 2486.5d;
            Assert.AreEqual((qVar99 - true_qVar99) / true_qVar99 < 0.01d, true);
        }

        /// <summary>
        /// Tested against Palisade's @Risk.
        /// </summary>
        [TestMethod()]
        public void Test_Weibull_GOF()
        {
            var W = new Weibull(9.589d, 1.907d);
            double logL = W.LogLikelihood(sample);
            double AIC = GoodnessOfFit.AIC(W.NumberOfParameters, logL);
            double BIC = GoodnessOfFit.BIC(sample.Length, W.NumberOfParameters, logL);
            var pp = PlottingPositions.Weibull(sample.Length);
            var modeled = new double[(sample.Length)];
            Array.Sort(sample);
            for (int i = 0; i < sample.Length; i++)
                modeled[i] = W.CDF(sample[i]);
            double true_AIC = 294.5878d;
            double true_BIC = 298.1566d;
            Assert.AreEqual((AIC - true_AIC) / true_AIC < 0.01d, true);
            Assert.AreEqual((BIC - true_BIC) / true_BIC < 0.01d, true);
        }

        /// <summary>
        /// Verifying input parameters can create Weibull.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var W = new Weibull();
            Assert.AreEqual(W.Lambda, 10);
            Assert.AreEqual(W.Kappa, 2);

            var W2 = new Weibull(1, 1);
            Assert.AreEqual(W2.Lambda, 1);
            Assert.AreEqual(W2.Kappa, 1);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var W = new Weibull(double.NaN, double.NaN);
            Assert.IsFalse(W.ParametersValid);

            var W2 = new Weibull(double.NegativeInfinity,double.PositiveInfinity);
            Assert.IsFalse(W2.ParametersValid);

            var W3 = new Weibull(0, 1);
            Assert.IsFalse(W3.ParametersValid);

            var W4 = new Weibull(1, 0);
            Assert.IsFalse(W4.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var W = new Weibull();
            Assert.AreEqual(W.ParametersToString[0, 0], "Scale (λ)");
            Assert.AreEqual(W.ParametersToString[1, 0], "Shape (κ)");
            Assert.AreEqual(W.ParametersToString[0, 1], "10");
            Assert.AreEqual(W.ParametersToString[1, 1], "2");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Weibull(1, 1);
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
            var W = new Weibull(0.1, 1);
            Assert.AreEqual(W.Mean, 0.1);

            var W2 = new Weibull(1, 1);
            Assert.AreEqual(W2.Mean, 1);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var W = new Weibull(0.1, 1);
            Assert.AreEqual(W.Median, 0.06931, 1e-04);

            var W2 = new Weibull(1, 1);
            Assert.AreEqual(W2.Median, 0.69314, 1e-04);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var W = new Weibull(0.1, 1);
            Assert.AreEqual(W.Mode, 0);

            var W2 = new Weibull(10, 10);
            Assert.AreEqual(W2.Mode, 9.89519, 1e-05);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var W = new Weibull(0.1, 1);
            Assert.AreEqual(W.StandardDeviation, 0.1);

            var W2 = new Weibull(1, 1);
            Assert.AreEqual(W2.StandardDeviation, 1);
        }

        /// <summary>
        /// Testing skew
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var W = new Weibull(0.1, 1);
            Assert.AreEqual(W.Skewness, 2,1e-04);

            var W2 = new Weibull(1, 1);
            Assert.AreEqual(W2.Skewness, 2);
        }

        /// <summary>
        /// Testing kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var W = new Weibull();
            Assert.AreEqual(W.Kurtosis, 3.24508,1e-04);

            var W2 = new Weibull(1, 1);
            Assert.AreEqual(W2.Kurtosis, 9);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var W = new Weibull();
            Assert.AreEqual(W.Minimum, 0);
            Assert.AreEqual(W.Maximum,double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var W = new Weibull(1, 1);
            Assert.AreEqual(W.PDF(0), 1);
            Assert.AreEqual(W.PDF(1), 0.36787, 1e-05);
            Assert.AreEqual(W.PDF(10), 0.00004539, 1e-08);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var W = new Weibull(1, 1);
            Assert.AreEqual(W.CDF(0), 0);
            Assert.AreEqual(W.CDF(1), 0.63212, 1e-05);
            Assert.AreEqual(W.CDF(10), 0.99995, 1e-05);
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var W = new Weibull();
            Assert.AreEqual(W.InverseCDF(0),0);
            Assert.AreEqual(W.InverseCDF(1),double.PositiveInfinity);
            Assert.AreEqual(W.InverseCDF(0.4), 7.1472, 1e-04);
        }
    }
}
