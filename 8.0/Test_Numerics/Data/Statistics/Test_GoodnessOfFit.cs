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

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the GoodnessOfFit class. These methods were tested against various R methods of the "qpcR", "Metrics", "stats", and "nortest" packages. The 
    /// specific functions used are documented in each test.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// Spiess A (2018). qpcR: Modelling and Analysis of Real-Time PCR Data. R package version 1.4-1, <see href="https://CRAN.R-project.org/package=qpcR."/>
    /// </item>
    /// <item>
    /// Hamner B, Frasco M (2018). Metrics: Evaluation Metrics for Machine Learning. R package version 0.1.4, <see href="https://CRAN.R-project.org/package=Metrics"/>
    /// </item>
    /// <item>
    /// R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, 
    /// Vienna, Austria. ISBN 3-900051-07-0, URL <see href="http://www.R-project.org/"/>
    /// </item>
    /// <item>
    /// Gross J, Ligges U (2015). nortest: Tests for Normality. R package version 1.0-4, <see href="https://CRAN.R-project.org/package=nortest"/>
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
        /// Creating data to perform the GoodnessOfFit tests on
        /// </summary>
        public Test_GoodnessOfFit()
        {
            norm = new Normal(100, 15);
            for (int i = 1; i <= 30; i++)
                data[i - 1] = norm.InverseCDF((double)i / 31);
            logL = norm.LogLikelihood(data);
        }

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
        /// Test the method for weighting AIC values. These values were tested against R's "akaike.weights()" function
        /// from the "qpcR" package.
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

        /// <summary>
        /// Test the RMSE method that takes a list of observed values, list of model values, and the number of model parameters. These values were tested 
        /// against R's "rmse()" function from the "Metrics" package.
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
        /// Test RMSE method that takes a list of observed values and the continuous distribution that is the model. These values were tested 
        /// against R's "rmse()" function from the "Metrics" package.
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
        /// Test the RMSE method that takes a list of observed values,list of plotting positions, and the continues distribution that is the model. These values were tested
        /// against R's "rmse()" function from the "Metrics" package.
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
        /// Test the method for weighting RMSE values.  Validation values was were directly from the formula for inverse-MSE weighting (the method used by this function).
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
        /// Test the R^2 method, which is the correlation coefficient (r) squared. This value was tested against R's "cor()" function of the
        /// "stats" package that was then squared.
        /// </summary>
        [TestMethod]
        public void Test_RSqaured()
        {
            var observed = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };
            double r2 = GoodnessOfFit.RSquared(observed, data);
            double trueR2 = 0.9803475;

            Assert.AreEqual(trueR2, r2, 1E-6);
        }

        /// <summary>
        /// Test the Chi-Squared test statistic. This method was tested against R's "gofTest()" method from the "EnvStats" package.
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
        /// Test the Kolmogrov Smirnov test statistic. This method was tested against R's "ks.test()" method from the "stats" package.
        /// </summary>
        [TestMethod]
        public void Test_KSTest()
        {
            var result = GoodnessOfFit.KolmogorovSmirnov(data, norm);
            Assert.AreEqual(0.032258, result, 1E-6);
        }

        /// <summary>
        /// Tested the Anderson Darling test. This method was tested against R's "ad.test()" function of the "nortest" package.
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
    }
}