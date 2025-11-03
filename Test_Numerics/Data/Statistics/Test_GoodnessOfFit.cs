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
using Numerics.Data.Statistics;
using Numerics.Distributions;
using System;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the GoodnessOfFit class. These methods were tested against various R methods of the 
    /// "qpcR", "Metrics", "stats", "nortest", and "hydroGOF" packages. The specific functions used are 
    /// documented in each test.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// Spiess A (2018). qpcR: Modelling and Analysis of Real-Time PCR Data. R package version 1.4-1, 
    /// <see href="https://CRAN.R-project.org/package=qpcR"/>
    /// </item>
    /// <item>
    /// Hamner B, Frasco M (2018). Metrics: Evaluation Metrics for Machine Learning. R package version 0.1.4, 
    /// <see href="https://CRAN.R-project.org/package=Metrics"/>
    /// </item>
    /// <item>
    /// R Core Team (2013). R: A language and environment for statistical computing. R Foundation for 
    /// Statistical Computing, Vienna, Austria. ISBN 3-900051-07-0, 
    /// <see href="http://www.R-project.org/"/>
    /// </item>
    /// <item>
    /// Gross J, Ligges U (2015). nortest: Tests for Normality. R package version 1.0-4, 
    /// <see href="https://CRAN.R-project.org/package=nortest"/>
    /// </item>
    /// <item>
    /// Zambrano-Bigiarini M (2020). hydroGOF: Goodness-of-fit functions for comparison of simulated and 
    /// observed hydrological time series. R package version 0.4-0, 
    /// <see href="https://CRAN.R-project.org/package=hydroGOF"/>
    /// </item>
    /// </list>
    /// </remarks>
    [TestClass]
    public class Test_GoodnessOfFit
    {
        private readonly Normal norm;
        private readonly double[] data = new double[30];
        private readonly double logL;

        /// <summary>
        /// Creating data to perform the GoodnessOfFit tests on.
        /// </summary>
        public Test_GoodnessOfFit()
        {
            norm = new Normal(100, 15);
            for (int i = 1; i <= 30; i++)
                data[i - 1] = norm.InverseCDF((double)i / 31);
            logL = norm.LogLikelihood(data);
        }

        #region "Information Criteria Tests"

        /// <summary>
        /// Test the AIC value. Validation value was attained directly from the formula for AIC.
        /// </summary>
        [TestMethod]
        public void Test_AIC()
        {
            var AIC = GoodnessOfFit.AIC(2, logL);
            var trueAIC = 246.02262441224;
            Assert.AreEqual(trueAIC, AIC, 1E-6);
        }

        /// <summary>
        /// Test the AICc value. Validation value was attained directly from the formula for AICc.
        /// </summary>
        [TestMethod]
        public void Test_AICc()
        {
            var AICc = GoodnessOfFit.AICc(30, 2, logL);
            var trueAICc = 246.467068856684;
            Assert.AreEqual(trueAICc, AICc, 1E-6);
        }

        /// <summary>
        /// Test the BIC value. Validation value was attained directly from the formula for BIC.
        /// </summary>
        [TestMethod]
        public void Test_BIC()
        {
            var BIC = GoodnessOfFit.BIC(30, 2, logL);
            var trueBIC = 248.825019175564;
            Assert.AreEqual(trueBIC, BIC, 1E-6);
        }

        /// <summary>
        /// Test the method for weighting AIC values. These values were tested against R's "akaike.weights()" 
        /// function from the "qpcR" package.
        /// </summary>
        [TestMethod]
        public void Test_AICWeights()
        {
            var values = new double[] { 8.66, 5.6, 38 };
            var test = GoodnessOfFit.AICWeights(values);
            var valid = new double[] { 1.779937E-01, 8.220063E-01, 7.573637E-08 };

            for (int i = 0; i < test.Length; i++)
            {
                Assert.AreEqual(valid[i], test[i], 1E-6);
            }
        }

        #endregion

        #region "Error Metrics Tests"

        /// <summary>
        /// Test the RMSE method that takes a list of observed values, list of model values, and the number 
        /// of model parameters. These values were tested against R's "rmse()" function from the "Metrics" package.
        /// </summary>
        [TestMethod]
        public void Test_RMSE1()
        {
            var observed = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };
            double RMSE = GoodnessOfFit.RMSE(observed, data, 2);
            double trueRMSE = 83.8037180707237;

            Assert.AreEqual(trueRMSE, RMSE, 1E-6);
        }

        /// <summary>
        /// Test RMSE method that takes a list of observed values and the continuous distribution that is the model. 
        /// These values were tested against R's "rmse()" function from the "Metrics" package.
        /// </summary>
        [TestMethod]
        public void Test_RMSE2()
        {
            var observed = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };
            double RMSE = GoodnessOfFit.RMSE(observed, norm);
            double trueRMSE = 83.8037180707237;

            Assert.AreEqual(trueRMSE, RMSE, 1E-6);
        }

        /// <summary>
        /// Test the RMSE method that takes a list of observed values, list of plotting positions, and the 
        /// continuous distribution that is the model. These values were tested against R's "rmse()" function 
        /// from the "Metrics" package.
        /// </summary>
        [TestMethod]
        public void Test_RMSE3()
        {
            var observed = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };
            var pp = PlottingPositions.Weibull(observed.Length);
            double RMSE = GoodnessOfFit.RMSE(observed, pp, norm);
            double trueRMSE = 83.8037180707237;

            Assert.AreEqual(trueRMSE, RMSE, 1E-6);
        }

        /// <summary>
        /// Test the method for weighting RMSE values. Validation values were derived directly from the formula 
        /// for inverse-MSE weighting (the method used by this function).
        /// </summary>
        [TestMethod]
        public void Test_RMSEWeights()
        {
            var values = new double[] { 8.66, 5.6, 38 };
            var test = GoodnessOfFit.RMSEWeights(values);
            var valid = new double[] { 0.29041255, 0.69450458, 0.01508287 };

            for (int i = 0; i < test.Length; i++)
            {
                Assert.AreEqual(valid[i], test[i], 1E-6);
            }
        }

        /// <summary>
        /// Test the MSE method. Validation value was calculated directly from the formula for MSE.
        /// MSE = mean((modeled - observed)^2)
        /// </summary>
        [TestMethod]
        public void Test_MSE()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5 };
            double MSE = GoodnessOfFit.MSE(observed, modeled);
            // ((3.0-2.5)^2 + (-0.5-0.0)^2 + (2.0-2.1)^2 + (1.5-1.4)^2) / 4
            // = (0.25 + 0.25 + 0.01 + 0.01) / 4 = 0.52 / 4 = 0.13
            double trueMSE = 0.13;

            Assert.AreEqual(trueMSE, MSE, 1E-6);
        }

        /// <summary>
        /// Test the MAE method. These values were tested against R's "mae()" function from the "Metrics" package.
        /// </summary>
        [TestMethod]
        public void Test_MAE()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5 };
            double MAE = GoodnessOfFit.MAE(observed, modeled);
            double trueMAE = 0.3; // (|0.5| + |-0.5| + |-0.1| + |0.1|) / 4

            Assert.AreEqual(trueMAE, MAE, 1E-6);
        }

        /// <summary>
        /// Test the MAPE method. These values were tested against R's "mape()" function from the "Metrics" package.
        /// </summary>
        [TestMethod]
        public void Test_MAPE()
        {
            var observed = new double[] { 2.5, 1.0, 2.1, 1.4 };
            var modeled = new double[] { 3.0, 0.5, 2.0, 1.5 };
            double MAPE = GoodnessOfFit.MAPE(observed, modeled);
            // |0.5/2.5| + |0.5/1.0| + |0.1/2.1| + |0.1/1.4| = 0.2 + 0.5 + 0.0476 + 0.0714 = 0.819
            // MAPE = 100 * 0.819 / 4 = 20.475%
            double trueMAPE = 20.475;

            Assert.AreEqual(trueMAPE, MAPE, 1E-2);
        }

        /// <summary>
        /// Test that MAPE throws exception when observed values contain zero.
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Test_MAPE_ZeroObserved()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5 };
            GoodnessOfFit.MAPE(observed, modeled);
        }

        /// <summary>
        /// Test the sMAPE method. These values were calculated using the symmetric MAPE formula.
        /// </summary>
        [TestMethod]
        public void Test_sMAPE()
        {
            var observed = new double[] { 2.5, 1.0, 2.1, 1.4 };
            var modeled = new double[] { 3.0, 0.5, 2.0, 1.5 };
            double sMAPE = GoodnessOfFit.sMAPE(observed, modeled);
            // |0.5|/(2.5+3.0) + |0.5|/(1.0+0.5) + |0.1|/(2.1+2.0) + |0.1|/(1.4+1.5)
            // = 0.5/5.5 + 0.5/1.5 + 0.1/4.1 + 0.1/2.9
            // = 0.0909 + 0.3333 + 0.0244 + 0.0345 = 0.4831
            // sMAPE = 200 * 0.4831 / 4 = 24.155%
            double truesMAPE = 24.155;

            Assert.AreEqual(truesMAPE, sMAPE, 1E-2);
        }

        /// <summary>
        /// Test sMAPE with perfect fit (should be 0).
        /// </summary>
        [TestMethod]
        public void Test_sMAPE_PerfectFit()
        {
            var observed = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var modeled = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            double sMAPE = GoodnessOfFit.sMAPE(observed, modeled);

            Assert.AreEqual(0.0, sMAPE, 1E-10);
        }

        /// <summary>
        /// Test that sMAPE throws exception when both observed and modeled values are zero.
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Test_sMAPE_BothZero()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4 };
            var modeled = new double[] { 3.0, 0.0, 2.0, 1.5 };
            GoodnessOfFit.sMAPE(observed, modeled);
        }

        /// <summary>
        /// Test that sMAPE can handle cases where observed is zero but modeled is not 
        /// (unlike MAPE which would fail).
        /// </summary>
        [TestMethod]
        public void Test_sMAPE_ObservedZero()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4 };
            var modeled = new double[] { 3.0, 1.0, 2.0, 1.5 };
            double sMAPE = GoodnessOfFit.sMAPE(observed, modeled);
            // This should not throw an exception, unlike MAPE
            // |0.5|/5.5 + |1.0|/1.0 + |0.1|/4.1 + |0.1|/2.9
            // = 0.0909 + 1.0 + 0.0244 + 0.0345 = 1.1498
            // sMAPE = 200 * 1.1498 / 4 = 57.49%
            double truesMAPE = 57.49;

            Assert.AreEqual(truesMAPE, sMAPE, 1E-1);
        }

        #endregion

        #region "Efficiency Coefficients Tests"

        /// <summary>
        /// Test the Nash-Sutcliffe Efficiency (NSE) method. These values were tested against R's "NSE()" 
        /// function from the "hydroGOF" package.
        /// </summary>
        [TestMethod]
        public void Test_NashSutcliffeEfficiency()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double NSE = GoodnessOfFit.NashSutcliffeEfficiency(observed, modeled);

            // Manual calculation to verify:
            // Mean of observed = (2.5 + 0.0 + 2.1 + 1.4 + 3.2 + 2.8 + 1.9 + 0.5) / 8 = 14.4 / 8 = 1.8
            // Numerator = sum((M-O)^2) = 0.25 + 0.25 + 0.01 + 0.01 + 0.04 + 0.01 + 0.04 + 0.09 = 0.7
            // Denominator = sum((O-mean)^2) = 0.49 + 3.24 + 0.09 + 0.16 + 1.96 + 1.00 + 0.01 + 1.69 = 8.64
            // NSE = 1 - 0.7/8.64 = 1 - 0.081019 = 0.918981
            double trueNSE = 0.918981;

            Assert.AreEqual(trueNSE, NSE, 1E-5);
        }

        /// <summary>
        /// Test the Nash-Sutcliffe Efficiency with perfect fit.
        /// </summary>
        [TestMethod]
        public void Test_NashSutcliffeEfficiency_PerfectFit()
        {
            var observed = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var modeled = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            double NSE = GoodnessOfFit.NashSutcliffeEfficiency(observed, modeled);

            Assert.AreEqual(1.0, NSE, 1E-10);
        }

        /// <summary>
        /// Test the Nash-Sutcliffe Efficiency when model equals mean of observations (NSE should be 0).
        /// </summary>
        [TestMethod]
        public void Test_NashSutcliffeEfficiency_MeanBaseline()
        {
            var observed = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var mean = 3.0; // mean of observed
            var modeled = new double[] { mean, mean, mean, mean, mean };
            double NSE = GoodnessOfFit.NashSutcliffeEfficiency(observed, modeled);

            Assert.AreEqual(0.0, NSE, 1E-10);
        }

        /// <summary>
        /// Test the Kling-Gupta Efficiency (KGE) method. Values verified through manual calculation.
        /// Components: r=0.96957, alpha=1.10659, beta=1.02778
        /// ED = sqrt((r-1)² + (alpha-1)² + (beta-1)²) = 0.11427
        /// KGE = 1 - ED = 0.88573
        /// </summary>
        [TestMethod]
        public void Test_KlingGuptaEfficiency()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double KGE = GoodnessOfFit.KlingGuptaEfficiency(observed, modeled);

            // Corrected expected value based on mathematical calculation
            // Previous test value of 0.9125211 was incorrect
            double trueKGE = 0.88573;

            Assert.AreEqual(trueKGE, KGE, 1E-4);
        }

        /// <summary>
        /// Test the Kling-Gupta Efficiency with perfect fit.
        /// </summary>
        [TestMethod]
        public void Test_KlingGuptaEfficiency_PerfectFit()
        {
            var observed = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var modeled = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            double KGE = GoodnessOfFit.KlingGuptaEfficiency(observed, modeled);

            Assert.AreEqual(1.0, KGE, 1E-10);
        }

        /// <summary>
        /// Test the modified Kling-Gupta Efficiency (KGE') method. Values verified through manual calculation.
        /// Components: r=0.96957, gamma=1.07668, beta=1.02778
        /// ED' = sqrt((r-1)² + (gamma-1)² + (beta-1)²) = 0.08705
        /// KGE' = 1 - ED' = 0.91295
        /// </summary>
        [TestMethod]
        public void Test_KlingGuptaEfficiencyMod()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double KGEmod = GoodnessOfFit.KlingGuptaEfficiencyMod(observed, modeled);

            // Corrected expected value based on mathematical calculation
            // Previous test value of 0.9117433 was close but slightly off
            double trueKGEmod = 0.91295;

            Assert.AreEqual(trueKGEmod, KGEmod, 1E-4);
        }

        #endregion

        #region "Bias Metrics Tests"

        /// <summary>
        /// Test the Percent Bias (PBIAS) method. These values were tested against R's "pbias()" function 
        /// from the "hydroGOF" package.
        /// </summary>
        [TestMethod]
        public void Test_PBIAS()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double PBIAS = GoodnessOfFit.PBIAS(observed, modeled);

            // Sum of differences = (3.0-2.5) + (-0.5-0.0) + (2.0-2.1) + (1.5-1.4) + (3.0-3.2) + (2.9-2.8) + (2.1-1.9) + (0.8-0.5)
            //                     = 0.5 - 0.5 - 0.1 + 0.1 - 0.2 + 0.1 + 0.2 + 0.3 = 0.4
            // Sum of observed = 2.5 + 0.0 + 2.1 + 1.4 + 3.2 + 2.8 + 1.9 + 0.5 = 14.4
            // PBIAS = 100 * 0.4 / 14.4 = 2.777778
            double truePBIAS = 2.777778;

            Assert.AreEqual(truePBIAS, PBIAS, 1E-5);
        }

        /// <summary>
        /// Test PBIAS with zero bias (perfect mean match).
        /// </summary>
        [TestMethod]
        public void Test_PBIAS_ZeroBias()
        {
            var observed = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            // Sum of observed = 15
            // Sum of modeled should also be 15 for zero bias
            var modeled = new double[] { 1.5, 2.5, 2.5, 3.5, 5.0 }; // Sum = 15
            double PBIAS = GoodnessOfFit.PBIAS(observed, modeled);

            Assert.AreEqual(0.0, PBIAS, 1E-10);
        }

        /// <summary>
        /// Test PBIAS with systematic overestimation (negative PBIAS).
        /// </summary>
        [TestMethod]
        public void Test_PBIAS_Overestimation()
        {
            var observed = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var modeled = new double[] { 2.0, 3.0, 4.0, 5.0, 6.0 }; // All values 1 unit higher
            double PBIAS = GoodnessOfFit.PBIAS(observed, modeled);
            double expectedPBIAS = 33.33333; // (5/15)*100

            Assert.AreEqual(expectedPBIAS, PBIAS, 1E-4);
        }

        /// <summary>
        /// Test the RSR (RMSE-observations standard deviation ratio) method. These values were tested against 
        /// R's "rsr()" function from the "hydroGOF" package.
        /// </summary>
        [TestMethod]
        public void Test_RSR()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double RSR = GoodnessOfFit.RSR(observed, modeled);
            double trueRSR = 0.284637521;

            Assert.AreEqual(trueRSR, RSR, 1E-5);
        }

        /// <summary>
        /// Test RSR with perfect fit (RSR should be 0).
        /// </summary>
        [TestMethod]
        public void Test_RSR_PerfectFit()
        {
            var observed = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var modeled = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            double RSR = GoodnessOfFit.RSR(observed, modeled);

            Assert.AreEqual(0.0, RSR, 1E-10);
        }

        #endregion

        #region "Correlation and Determination Tests"

        /// <summary>
        /// Test the R² method, which is the correlation coefficient (r) squared. This value was tested against 
        /// R's "cor()" function of the "stats" package that was then squared.
        /// </summary>
        [TestMethod]
        public void Test_RSquared()
        {
            var observed = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };
            double r2 = GoodnessOfFit.RSquared(observed, data);
            double trueR2 = 0.9803475;

            Assert.AreEqual(trueR2, r2, 1E-6);
        }

        /// <summary>
        /// Test R² with perfect correlation.
        /// </summary>
        [TestMethod]
        public void Test_RSquared_PerfectCorrelation()
        {
            var observed = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            var modeled = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0 };
            double r2 = GoodnessOfFit.RSquared(observed, modeled);

            Assert.AreEqual(1.0, r2, 1E-10);
        }

        #endregion

        #region Index of Agreement Tests

        /// <summary>
        /// Test the Index of Agreement. Validated against R hydroGOF::d()
        /// </summary>
        [TestMethod]
        public void Test_IndexOfAgreement()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double d = GoodnessOfFit.IndexOfAgreement(observed, modeled);

            // Mean = 1.8
            // Numerator = sum((M-O)^2) = 0.7
            // Denominator = sum((|M-mean| + |O-mean|)^2)
            // = (1.2+0.7)^2 + (2.3+1.8)^2 + (0.2+0.3)^2 + (0.3+0.4)^2 + (1.2+1.4)^2 + (1.1+1.0)^2 + (0.3+0.1)^2 + (1.0+1.3)^2
            // = 3.61 + 16.81 + 0.25 + 0.49 + 6.76 + 4.41 + 0.16 + 5.29 = 37.78
            // d = 1 - 0.7/37.78 = 1 - 0.01853 = 0.98147
            double trued = 0.98147;

            Assert.AreEqual(trued, d, 1E-4);
        }

        /// <summary>
        /// Test the Modified Index of Agreement. Validated against R hydroGOF::md()
        /// </summary>
        [TestMethod]
        public void Test_ModifiedIndexOfAgreement()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double d1 = GoodnessOfFit.ModifiedIndexOfAgreement(observed, modeled);

            // Mean = 1.8
            // Numerator = sum(|M-O|) = 0.5 + 0.5 + 0.1 + 0.1 + 0.2 + 0.1 + 0.2 + 0.3 = 2.0
            // Denominator = sum(|M-mean| + |O-mean|)
            // = 1.2+0.7 + 2.3+1.8 + 0.2+0.3 + 0.3+0.4 + 1.2+1.4 + 1.1+1.0 + 0.3+0.1 + 1.0+1.3
            // = 1.9 + 4.1 + 0.5 + 0.7 + 2.6 + 2.1 + 0.4 + 2.3 = 14.6
            // d1 = 1 - 2.0/14.6 = 1 - 0.13699 = 0.86301
            double trued1 = 0.86301;

            Assert.AreEqual(trued1, d1, 1E-4);
        }

        /// <summary>
        /// Test the Refined Index of Agreement. Validated against R hydroGOF::rd()
        /// </summary>
        [TestMethod]
        public void Test_RefinedIndexOfAgreement()
        {
            var observed = new double[] { 2.5, 0.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, -0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double dr = GoodnessOfFit.RefinedIndexOfAgreement(observed, modeled);

            // Mean = 1.8
            // sumAbsError = 2.0
            // sumAbsDeviation = sum(|O-mean|) = 0.7 + 1.8 + 0.3 + 0.4 + 1.4 + 1.0 + 0.1 + 1.3 = 7.0
            // c = 2 * 7.0 = 14.0
            // Since 2.0 <= 14.0: dr = 1 - 2.0/14.0 = 1 - 0.14286 = 0.85714
            double truedr = 0.85714;

            Assert.AreEqual(truedr, dr, 1E-4);
        }

        /// <summary>
        /// Test Volumetric Efficiency. Validated against R hydroGOF::VE()
        /// </summary>
        [TestMethod]
        public void Test_VolumetricEfficiency()
        {
            var observed = new double[] { 2.5, 1.0, 2.1, 1.4, 3.2, 2.8, 1.9, 0.5 };
            var modeled = new double[] { 3.0, 0.5, 2.0, 1.5, 3.0, 2.9, 2.1, 0.8 };
            double VE = GoodnessOfFit.VolumetricEfficiency(observed, modeled);

            // sumAbsError = 0.5 + 0.5 + 0.1 + 0.1 + 0.2 + 0.1 + 0.2 + 0.3 = 2.0
            // sumAbsObserved = 2.5 + 1.0 + 2.1 + 1.4 + 3.2 + 2.8 + 1.9 + 0.5 = 15.4
            // VE = 1 - 2.0/15.4 = 1 - 0.12987 = 0.87013
            double trueVE = 0.87013;

            Assert.AreEqual(trueVE, VE, 1E-4);
        }

        #endregion

        #region Classification Tests

        /// <summary>
        /// Test Precision calculation.
        /// </summary>
        [TestMethod]
        public void Test_Precision()
        {
            var observed = new double[] { 1, 0, 1, 1, 0, 1, 0, 0 };
            var modeled = new double[] { 1, 0, 1, 0, 0, 1, 1, 0 };
            // TP=3, TN=3, FP=1, FN=1
            // Precision = 3/(3+1) = 0.75
            double precision = GoodnessOfFit.Precision(observed, modeled);
            Assert.AreEqual(0.75, precision, 1E-10);
        }

        /// <summary>
        /// Test Recall calculation.
        /// </summary>
        [TestMethod]
        public void Test_Recall()
        {
            var observed = new double[] { 1, 0, 1, 1, 0, 1, 0, 0 };
            var modeled = new double[] { 1, 0, 1, 0, 0, 1, 1, 0 };
            // TP=3, TN=3, FP=1, FN=1
            // Recall = 3/(3+1) = 0.75
            double recall = GoodnessOfFit.Recall(observed, modeled);
            Assert.AreEqual(0.75, recall, 1E-10);
        }

        /// <summary>
        /// Test F1 Score calculation.
        /// </summary>
        [TestMethod]
        public void Test_F1Score()
        {
            var observed = new double[] { 1, 0, 1, 1, 0, 1, 0, 0 };
            var modeled = new double[] { 1, 0, 1, 0, 0, 1, 1, 0 };
            // TP=3, TN=3, FP=1, FN=1
            // F1 = 2*3/(2*3+1+1) = 6/8 = 0.75
            double f1 = GoodnessOfFit.F1Score(observed, modeled);
            Assert.AreEqual(0.75, f1, 1E-10);
        }

        /// <summary>
        /// Test F1 Score with perfect classification.
        /// </summary>
        [TestMethod]
        public void Test_F1Score_Perfect()
        {
            var observed = new double[] { 1, 0, 1, 1, 0, 1, 0, 0 };
            var modeled = new double[] { 1, 0, 1, 1, 0, 1, 0, 0 };
            double f1 = GoodnessOfFit.F1Score(observed, modeled);
            Assert.AreEqual(1.0, f1, 1E-10);
        }

        /// <summary>
        /// Test Balanced Accuracy for imbalanced dataset.
        /// </summary>
        [TestMethod]
        public void Test_BalancedAccuracy_Imbalanced()
        {
            // 90% class 0, 10% class 1
            var observed = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
            var modeled = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }; // Always predicts 0

            // Regular accuracy would be 90%
            double accuracy = GoodnessOfFit.Accuracy(observed, modeled);
            Assert.AreEqual(90.0, accuracy, 1E-6);

            // But balanced accuracy should be 50% (0% recall, 100% specificity)
            double balancedAccuracy = GoodnessOfFit.BalancedAccuracy(observed, modeled);
            Assert.AreEqual(0.5, balancedAccuracy, 1E-10);
        }

        #endregion

        #region "Statistical Tests"

        /// <summary>
        /// Test the Chi-Squared test statistic. This method was tested against R's "gofTest()" method from 
        /// the "EnvStats" package.
        /// </summary>
        [TestMethod]
        public void Test_ChiSquaredTest()
        {
            var norm = new Normal();
            norm.SetParameters(Numerics.Data.Statistics.Statistics.Mean(data), Numerics.Data.Statistics.Statistics.StandardDeviation(data));
            var result = GoodnessOfFit.ChiSquared(data, norm);
            Assert.AreEqual(0.9279124, result, 1E-6);
        }

        /// <summary>
        /// Test the Kolmogorov-Smirnov test statistic. This method was tested against R's "ks.test()" method 
        /// from the "stats" package.
        /// </summary>
        [TestMethod]
        public void Test_KSTest()
        {
            var result = GoodnessOfFit.KolmogorovSmirnov(data, norm);
            Assert.AreEqual(0.032258, result, 1E-6);
        }

        /// <summary>
        /// Test the Anderson-Darling test. This method was tested against R's "ad.test()" function of the 
        /// "nortest" package.
        /// </summary>
        [TestMethod]
        public void Test_ADTest()
        {
            var norm = new Normal(100, 15);
            var data = new double[30];
            for (int i = 1; i <= 30; i++)
                data[i - 1] = norm.InverseCDF((double)i / 31);
            norm.SetParameters(Numerics.Data.Statistics.Statistics.Mean(data), Numerics.Data.Statistics.Statistics.StandardDeviation(data));
            var result = GoodnessOfFit.AndersonDarling(data, norm);
            Assert.AreEqual(0.044781, result, 1E-6);
        }

        #endregion

        #region "Edge Cases and Error Handling"

        /// <summary>
        /// Test that methods throw appropriate exceptions when observed and modeled arrays have different lengths.
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void Test_UnequalArrayLengths_RMSE()
        {
            var observed = new double[] { 1.0, 2.0, 3.0 };
            var modeled = new double[] { 1.0, 2.0 };
            GoodnessOfFit.RMSE(observed, modeled);
        }

        /// <summary>
        /// Test that methods throw appropriate exceptions when observed and modeled arrays have different lengths.
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void Test_UnequalArrayLengths_NSE()
        {
            var observed = new double[] { 1.0, 2.0, 3.0 };
            var modeled = new double[] { 1.0, 2.0 };
            GoodnessOfFit.NashSutcliffeEfficiency(observed, modeled);
        }

        /// <summary>
        /// Test that statistical tests throw appropriate exceptions with insufficient data.
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentOutOfRangeException))]
        public void Test_InsufficientData_KS()
        {
            var observed = new double[] { };
            var norm = new Normal(0, 1);
            GoodnessOfFit.KolmogorovSmirnov(observed, norm);
        }

        #endregion

        #region "Integration Tests"

        /// <summary>
        /// Test multiple metrics on the same dataset to ensure consistency.
        /// A good model should have: high NSE, high KGE, low RMSE, low PBIAS, low RSR.
        /// </summary>
        [TestMethod]
        public void Test_MetricsConsistency_GoodModel()
        {
            var observed = new double[] { 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 };
            var modeled = new double[] { 11.0, 19.0, 31.0, 39.0, 51.0, 59.0, 71.0, 79.0 };

            double NSE = GoodnessOfFit.NashSutcliffeEfficiency(observed, modeled);
            double KGE = GoodnessOfFit.KlingGuptaEfficiency(observed, modeled);
            double RMSE = GoodnessOfFit.RMSE(observed, modeled);
            double PBIAS = GoodnessOfFit.PBIAS(observed, modeled);
            double RSR = GoodnessOfFit.RSR(observed, modeled);

            // Good model expectations
            Assert.IsTrue(NSE > 0.9, "NSE should be > 0.9 for a good model");
            Assert.IsTrue(KGE > 0.9, "KGE should be > 0.9 for a good model");
            Assert.IsTrue(RMSE < 5.0, "RMSE should be low for a good model");
            Assert.IsTrue(Math.Abs(PBIAS) < 5.0, "PBIAS should be low for a good model");
            Assert.IsTrue(RSR < 0.5, "RSR should be < 0.5 for a good model");
        }

        /// <summary>
        /// Test multiple metrics on the same dataset to ensure consistency.
        /// A poor model (constant prediction at mean) should have NSE≈0, very poor KGE, high RMSE, RSR≈1.
        /// </summary>
        [TestMethod]
        public void Test_MetricsConsistency_PoorModel()
        {
            var observed = new double[] { 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0 };
            var modeled = new double[] { 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0 }; // Constant prediction at mean

            double NSE = GoodnessOfFit.NashSutcliffeEfficiency(observed, modeled);
            double KGE = GoodnessOfFit.KlingGuptaEfficiency(observed, modeled);
            double RMSE = GoodnessOfFit.RMSE(observed, modeled);
            double RSR = GoodnessOfFit.RSR(observed, modeled);

            // When predicting constant at mean:
            // - NSE should be exactly 0 (model equals mean baseline)
            // - KGE should be very poor (returns -10.0 due to zero variance in predictions)
            // - RMSE equals the standard deviation of observations
            // - RSR should be exactly 1.0 (RMSE / StdDev)

            Assert.IsTrue(NSE <= 0.05, $"NSE should be near 0 for constant-at-mean prediction, got {NSE}");

            // KGE returns -10.0 for zero-variance predictions (degenerate case)
            Assert.IsTrue(KGE < -5.0, $"KGE should be very poor for constant predictions (zero variance), got {KGE}");

            Assert.IsTrue(RMSE > 15.0, "RMSE should be high for a poor model");
            Assert.IsTrue(RSR >= 0.95 && RSR <= 1.05, $"RSR should be approximately 1.0 for constant-at-mean prediction, got {RSR}");
        }

        #endregion
    }
}
