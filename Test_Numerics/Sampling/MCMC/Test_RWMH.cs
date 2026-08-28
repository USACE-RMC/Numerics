using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit test for the Random Walk Metropolis-Hastings (RWMH) sampler. 
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
    public class Test_RWMH
    {

        /// <summary>
        /// This test compares the results obtained using RWMH with those from the 'rstan' package. 
        /// </summary>
        [TestMethod]
        public void Test_RWMH_NormalDist_RStan()
        {

            // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
            // Table 5.1.1 Tippecanoe River Near Delphi, Indiana (Station 43) Data
            double[] sample = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

            // Create uniform priors
            var normDist = new Normal();
            var constraints = normDist.GetParameterConstraints(sample);
            var muPrior = new Uniform(constraints.Item2[0], constraints.Item3[0]);
            var sigmaPrior = new Uniform(constraints.Item2[1], constraints.Item3[1]);
            var priors = new List<IUnivariateDistribution> { muPrior, sigmaPrior };

            // Create log-likelihood function
            double logLH(double[] x)
            {
                var dist = new Normal(x[0], x[1]);
                return dist.LogLikelihood(sample);
            }

            // Create and run sampler
            var sampler = new RWMH(priors, logLH, new Matrix(2));
            sampler.Initialize = MCMCSampler.InitializationType.MAP;
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
        /// Pins seeded RWMH draws bitwise: a fixed seed and configuration must reproduce every drawn
        /// value, fitness, and accept count exactly, on every platform.
        /// </summary>
        /// <remarks>
        /// The proposal covariance is fixed for a whole RWMH run, so each chain's proposal is
        /// factorized once at initialization and only its mean moves with the chain state — a
        /// translation changes nothing the factorization derives from the covariance, so the drawn
        /// values are identical to re-parameterizing the proposal at every step. The expected
        /// literals below are a golden master captured at seed 12345, 2 chains, 200 iterations,
        /// 100 warmup, Randomize initialization, and every one must reproduce exactly — no
        /// tolerance.
        /// </remarks>
        [TestMethod]
        public void Test_RWMH_MeanOnlyProposalUpdate_ReproducesReferenceDrawsExactly()
        {
            double[] sample = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

            var normDist = new Normal();
            var constraints = normDist.GetParameterConstraints(sample);
            var muPrior = new Uniform(constraints.Item2[0], constraints.Item3[0]);
            var sigmaPrior = new Uniform(constraints.Item2[1], constraints.Item3[1]);
            var priors = new List<IUnivariateDistribution> { muPrior, sigmaPrior };

            double logLH(double[] x)
            {
                var dist = new Normal(x[0], x[1]);
                return dist.LogLikelihood(sample);
            }

            var proposal = new Matrix(2);
            proposal[0, 0] = 500d * 500d;
            proposal[1, 1] = 300d * 300d;

            var sampler = new RWMH(priors, logLH, proposal)
            {
                Initialize = MCMCSampler.InitializationType.Randomize,
                PRNGSeed = 12345,
                NumberOfChains = 2,
                Iterations = 200,
                WarmupIterations = 100,
                ThinningInterval = 1,
                OutputLength = 100,
            };
            sampler.Sample();

            Assert.AreEqual(11934.435149791721, sampler.Output[0][0].Values[0], 0d);
            Assert.AreEqual(4526.5732644023383, sampler.Output[0][0].Values[1], 0d);
            Assert.AreEqual(-474.22548961374116, sampler.Output[0][0].Fitness, 0d);
            Assert.AreEqual(12095.33865486568, sampler.Output[0][1].Values[0], 0d);
            Assert.AreEqual(4982.8206043145146, sampler.Output[0][1].Values[1], 0d);
            Assert.AreEqual(12317.827258667343, sampler.Output[0][2].Values[0], 0d);
            Assert.AreEqual(5032.5525299089159, sampler.Output[0][2].Values[1], 0d);
            Assert.AreEqual(12178.391385589634, sampler.Output[0][sampler.Output[0].Count - 1].Values[0], 0d);
            Assert.AreEqual(4339.4907777225499, sampler.Output[0][sampler.Output[0].Count - 1].Values[1], 0d);
            Assert.AreEqual(178, sampler.AcceptCount[0]);

            Assert.AreEqual(20109.085873154781, sampler.Output[1][0].Values[0], 0d);
            Assert.AreEqual(28890.963429995856, sampler.Output[1][0].Values[1], 0d);
            Assert.AreEqual(21498.176703183646, sampler.Output[1][2].Values[0], 0d);
            Assert.AreEqual(28440.722339269822, sampler.Output[1][2].Values[1], 0d);
            Assert.AreEqual(19608.601940250355, sampler.Output[1][sampler.Output[1].Count - 1].Values[0], 0d);
            Assert.AreEqual(25961.674646968706, sampler.Output[1][sampler.Output[1].Count - 1].Values[1], 0d);
            Assert.AreEqual(209, sampler.AcceptCount[1]);

            Assert.AreEqual(12729.094200262518, sampler.MAP.Values[0], 0d);
            Assert.AreEqual(4543.2396366299008, sampler.MAP.Values[1], 0d);
            Assert.AreEqual(-473.59481939858011, sampler.MAP.Fitness, 0d);
        }

    }
}
