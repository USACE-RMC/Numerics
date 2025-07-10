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
using Numerics.Sampling.MCMC;
using System.Collections.Generic;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit test for the Differential Evolution Monte Carlo with Snooker Update (DEMCzs) sampler.  
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_DEMCzs
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 5.1.1 Tippecanoe River Near Delphi, Indiana (Station 43) Data
        private double[] sample1 = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 7.2.1 Sugar Creek at Crawfordsville, IN
        private double[] sample2 = new double[] { 17600d, 3660d, 903d, 5050d, 24000d, 11400d, 9470d, 8970d, 7710d, 14800d, 13900d, 20800d, 9470d, 7860d, 7860d, 2730d, 6480d, 18200d, 26300d, 15100d, 14600d, 7300d, 8580d, 15100d, 15100d, 21800d, 6200d, 2130d, 11100d, 14300d, 11200d, 6670d, 5440d, 9370d, 6900d, 9680d, 6810d, 7730d, 5290d, 12200d, 9750d, 7390d, 13100d, 7190d, 8850d, 6290d, 18800d, 9740d, 2990d, 6950d, 9390d, 12400d, 21200d };

        /// <summary>
        /// This test compares the results obtained using DEMCzs for the Normal distribution with those from the 'rstan' package. 
        /// </summary>
        [TestMethod]
        public void Test_DEMCzs_NormalDist_RStan()
        {

            // Create uniform priors
            var normDist = new Normal();
            var constraints = normDist.GetParameterConstraints(sample1);
            var muPrior = new Uniform(constraints.Item2[0], constraints.Item3[0]);
            var sigmaPrior = new Uniform(constraints.Item2[1], constraints.Item3[1]);
            var priors = new List<IUnivariateDistribution> { muPrior, sigmaPrior };

            // Create log-likelihood function
            double logLH(double[] x)
            {
                var dist = new Normal(x[0], x[1]);
                return dist.LogLikelihood(sample1);
            }

            // Create and run MCMC sampler
            var sampler = new DEMCzs(priors, logLH);
            sampler.Sample();
            var results = new MCMCResults(sampler);

            /* Below are the results from 'rstan' using comparable MCMC settings:
            *            mean se_mean     sd       5%      50%      95% n_eff Rhat
            *  mu    12663.69    7.10 706.60 11488.50 12671.08 13801.45  9897    1
            *  sigma  4844.09    5.22 519.08  4077.80  4796.63  5771.81  9880    1
            *  lp__   -466.13    0.01   1.03  -468.17  -465.81  -465.15  9958    1
            * 
            *  Since MCMC methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 5% of 'rstan' results. 
            */

            // Mu 
            Assert.AreEqual(12663.69, results.ParameterResults[0].SummaryStatistics.Mean, 0.05 * 12663.69);
            Assert.AreEqual(706.60, results.ParameterResults[0].SummaryStatistics.StandardDeviation, 0.05 * 706.60);
            Assert.AreEqual(11488.50, results.ParameterResults[0].SummaryStatistics.LowerCI, 0.05 * 11488.50);
            Assert.AreEqual(12671.08, results.ParameterResults[0].SummaryStatistics.Median, 0.05 * 12671.08);
            Assert.AreEqual(13801.45, results.ParameterResults[0].SummaryStatistics.UpperCI, 0.05 * 13801.45);
            // Sigma 
            Assert.AreEqual(4844.09, results.ParameterResults[1].SummaryStatistics.Mean, 0.05 * 4844.09);
            Assert.AreEqual(519.08, results.ParameterResults[1].SummaryStatistics.StandardDeviation, 0.05 * 519.08);
            Assert.AreEqual(4077.80, results.ParameterResults[1].SummaryStatistics.LowerCI, 0.05 * 4077.80);
            Assert.AreEqual(4796.63, results.ParameterResults[1].SummaryStatistics.Median, 0.05 * 4796.63);
            Assert.AreEqual(5771.81, results.ParameterResults[1].SummaryStatistics.UpperCI, 0.05 * 5771.81);
        }

        /// <summary>
        /// This test compares the results obtained using DEMCzs for the Logistic distribution with those from the 'rstan' package. 
        /// </summary>
        [TestMethod]
        public void Test_DEMCzs_LogisticDist_RStan()
        {

            // Create uniform priors
            var logDist = new Logistic();
            var constraints = logDist.GetParameterConstraints(sample1);
            var xiPrior = new Uniform(constraints.Item2[0], constraints.Item3[0]);
            var alphaPrior = new Uniform(constraints.Item2[1], constraints.Item3[1]);
            var priors = new List<IUnivariateDistribution> { xiPrior, alphaPrior };

            // Create log-likelihood function
            double logLH(double[] x)
            {
                var dist = new Logistic(x[0], x[1]);
                return dist.LogLikelihood(sample1);
            }

            // Create and run sampler
            var sampler = new DEMCzs(priors, logLH);
            sampler.Sample();
            var results = new MCMCResults(sampler);

            /* Below are the results from 'rstan' using comparable MCMC settings:
            *            mean se_mean     sd       5%      50%      95% n_eff Rhat
            *  mu    12631.74    7.19 719.36 11458.52 12622.01 13830.74 10002    1
            *  sigma  2823.47    3.49 348.26  2307.77  2794.40  3441.31  9947    1
            *  lp__   -467.66    0.01   1.05  -469.75  -467.33  -466.69 10128    
            * 
            *  Since MCMC methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 5% of 'rstan' results. 
            */

            // Mu in R is equal to Xi in Numerics
            Assert.AreEqual(12631.74, results.ParameterResults[0].SummaryStatistics.Mean, 0.05 * 12631.74);
            Assert.AreEqual(719.36, results.ParameterResults[0].SummaryStatistics.StandardDeviation, 0.05 * 719.36);
            Assert.AreEqual(11458.52, results.ParameterResults[0].SummaryStatistics.LowerCI, 0.05 * 11458.52);
            Assert.AreEqual(12622.01, results.ParameterResults[0].SummaryStatistics.Median, 0.05 * 12622.01);
            Assert.AreEqual(13830.74, results.ParameterResults[0].SummaryStatistics.UpperCI, 0.05 * 13830.74);
            // Sigma in R is equal to Alpha in Numerics
            Assert.AreEqual(2823.47, results.ParameterResults[1].SummaryStatistics.Mean, 0.05 * 2823.47);
            Assert.AreEqual(348.26, results.ParameterResults[1].SummaryStatistics.StandardDeviation, 0.05 * 348.26);
            Assert.AreEqual(2307.77, results.ParameterResults[1].SummaryStatistics.LowerCI, 0.05 * 2307.77);
            Assert.AreEqual(2794.40, results.ParameterResults[1].SummaryStatistics.Median, 0.05 * 2794.40);
            Assert.AreEqual(3441.31, results.ParameterResults[1].SummaryStatistics.UpperCI, 0.05 * 3441.31);
        }

        /// <summary>
        /// This test compares the results obtained using DEMCzs for the Gumbel distribution with those from the 'rstan' package. 
        /// </summary>
        [TestMethod]
        public void Test_DEMCzs_GumbelDist_RStan()
        {

            // Create uniform priors
            var gumDist = new Gumbel();
            var constraints = gumDist.GetParameterConstraints(sample2);
            var xiPrior = new Uniform(constraints.Item2[0], constraints.Item3[0]);
            var alphaPrior = new Uniform(constraints.Item2[1], constraints.Item3[1]);
            var priors = new List<IUnivariateDistribution> { xiPrior, alphaPrior };

            // Create log-likelihood function
            double logLH(double[] x)
            {
                var dist = new Gumbel(x[0], x[1]);
                return dist.LogLikelihood(sample2);
            }

            // Create and run sampler
            var sampler = new DEMCzs(priors, logLH);
            sampler.Sample();
            var results = new MCMCResults(sampler);

            /* Below are the results from 'rstan' using comparable MCMC settings:
            *            mean se_mean     sd      5%      50%    95% n_eff Rhat
            *  mu     8067.90    6.78 680.85 6973.92 8055.82 9201.14 10095    1
            *  beta   4658.04    5.32 522.33 3885.53 4612.71 5595.40  9623    1
            *  lp__   -521.82    0.01   1.02 -523.81 -521.50 -520.84  9860    1    
            * 
            *  Since MCMC methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 5% of 'rstan' results. 
            */

            // Mu in R is equal to Xi in Numerics
            Assert.AreEqual(8067.90, results.ParameterResults[0].SummaryStatistics.Mean, 0.05 * 8067.90);
            Assert.AreEqual(680.85, results.ParameterResults[0].SummaryStatistics.StandardDeviation, 0.05 * 680.85);
            Assert.AreEqual(6973.92, results.ParameterResults[0].SummaryStatistics.LowerCI, 0.05 * 6973.92);
            Assert.AreEqual(8055.82, results.ParameterResults[0].SummaryStatistics.Median, 0.05 * 8055.82);
            Assert.AreEqual(9201.14, results.ParameterResults[0].SummaryStatistics.UpperCI, 0.05 * 9201.14);
            // Beta in R is equal to Alpha in Numerics
            Assert.AreEqual(4658.04, results.ParameterResults[1].SummaryStatistics.Mean, 0.05 * 4658.04);
            Assert.AreEqual(522.33, results.ParameterResults[1].SummaryStatistics.StandardDeviation, 0.05 * 522.33);
            Assert.AreEqual(3885.53, results.ParameterResults[1].SummaryStatistics.LowerCI, 0.05 * 3885.53);
            Assert.AreEqual(4612.71, results.ParameterResults[1].SummaryStatistics.Median, 0.05 * 4612.71);
            Assert.AreEqual(5595.40, results.ParameterResults[1].SummaryStatistics.UpperCI, 0.05 * 5595.40);
        }

        /// <summary>
        /// This test compares the results obtained using DEMCzs for the Weibull distribution with those from the 'rstan' package. 
        /// </summary>
        [TestMethod]
        public void Test_DEMCzs_WeibullDist_RStan()
        {

            // Create uniform priors
            var weiDist = new Weibull();
            var constraints = weiDist.GetParameterConstraints(sample1);
            var lambdaPrior = new Uniform(constraints.Item2[0], constraints.Item3[0]);
            var kappaPrior = new Uniform(constraints.Item2[1], constraints.Item3[1]);
            var priors = new List<IUnivariateDistribution> { lambdaPrior, kappaPrior };

            // Create log-likelihood function
            double logLH(double[] x)
            {
                var dist = new Weibull(x[0], x[1]);
                return dist.LogLikelihood(sample1);
            }

            // Create and run sampler
            var sampler = new DEMCzs(priors, logLH);
            sampler.Sample();
            var results = new MCMCResults(sampler);

            /* Below are the results from 'rstan' using comparable MCMC settings:
            *            mean se_mean     sd      5%       50%      95% n_eff Rhat
            *  alpha     2.98    0.00   0.34     2.44     2.97     3.55  9702    1
            *  sigma 14310.25    7.46 745.27 13103.65 14309.95 15534.76  9990    1
            *  lp__   -463.40    0.01   1.04  -465.48  -463.09  -462.44  9975    1    
            * 
            *  Since MCMC methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 5% of 'rstan' results. 
            */

            // Alpha from R is equal to Kappa in Numerics
            Assert.AreEqual(2.98, results.ParameterResults[1].SummaryStatistics.Mean, 0.05 * 2.98);
            Assert.AreEqual(0.34, results.ParameterResults[1].SummaryStatistics.StandardDeviation, 0.05 * 0.34);
            Assert.AreEqual(2.44, results.ParameterResults[1].SummaryStatistics.LowerCI, 0.05 * 2.44);
            Assert.AreEqual(2.97, results.ParameterResults[1].SummaryStatistics.Median, 0.05 * 2.97);
            Assert.AreEqual(3.55, results.ParameterResults[1].SummaryStatistics.UpperCI, 0.05 * 3.55);
            // Sigma from R is equal to Lambda in Numerics
            Assert.AreEqual(14310.25, results.ParameterResults[0].SummaryStatistics.Mean, 0.05 * 14310.25);
            Assert.AreEqual(745.27, results.ParameterResults[0].SummaryStatistics.StandardDeviation, 0.05 * 745.27);
            Assert.AreEqual(13103.65, results.ParameterResults[0].SummaryStatistics.LowerCI, 0.05 * 13103.65);
            Assert.AreEqual(14309.95, results.ParameterResults[0].SummaryStatistics.Median, 0.05 * 14309.95);
            Assert.AreEqual(15534.76, results.ParameterResults[0].SummaryStatistics.UpperCI, 0.05 * 15534.76);
        }

    }
}
