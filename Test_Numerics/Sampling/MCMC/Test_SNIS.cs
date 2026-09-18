using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Distributions;
using Numerics.Sampling;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit test for the Self-Normalizing Importance Sampler (SNIS). 
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
    public class Test_SNIS
    {
        /// <summary>
        /// This test compares the results obtained using SNIS with those from the 'rstan' package. 
        /// </summary>
        [TestMethod]
        public void Test_SNIS_NormalDist_RStan()
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
            var sampler = new SNIS(priors, logLH) { PRNGSeed = 45678 };
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
        /// This test verifies that SNIS fails fast when every importance weight is invalid.
        /// </summary>
        [TestMethod]
        public void Test_SNIS_AllInvalidWeights_ThrowsInvalidOperationException()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(0d, 1d) };
            double logLH(double[] x) => double.NegativeInfinity;

            var sampler = new SNIS(priors, logLH)
            {
                Iterations = 100,
                OutputLength = 100,
                PRNGSeed = 12345
            };

            Assert.Throws<System.InvalidOperationException>(() => sampler.Sample());
        }

        /// <summary>
        /// This test verifies that finite SNIS weights normalize cleanly when invalid samples are present.
        /// </summary>
        [TestMethod]
        public void Test_SNIS_MixedFiniteAndInvalidWeights_NormalizesFiniteWeights()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(0d, 1d) };
            double logLH(double[] x) => x[0] < 0.5d ? 0d : double.NegativeInfinity;

            var sampler = new SNIS(priors, logLH)
            {
                Iterations = 100,
                OutputLength = 100,
                PRNGSeed = 12345
            };

            sampler.Sample();

            double sum = 0d;
            bool foundZeroWeight = false;
            bool foundPositiveWeight = false;
            foreach (var parameterSet in sampler.MarkovChains[0])
            {
                Assert.IsFalse(double.IsNaN(parameterSet.Weight));
                Assert.IsFalse(double.IsInfinity(parameterSet.Weight));
                sum += parameterSet.Weight;
                foundZeroWeight |= parameterSet.Weight == 0d;
                foundPositiveWeight |= parameterSet.Weight > 0d;
            }

            Assert.IsTrue(foundZeroWeight);
            Assert.IsTrue(foundPositiveWeight);
            Assert.AreEqual(1d, sum, 1e-12);
        }

        /// <summary>
        /// SNIS identifies NumberOfChains when rejecting configurations with more than one chain.
        /// </summary>
        [TestMethod]
        public void Test_SNIS_MultipleChains_ReportsNumberOfChainsParameterName()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(0d, 1d) };
            var sampler = new SNIS(priors, x => 0d)
            {
                NumberOfChains = 2,
                InitialIterations = 2,
                Iterations = 100,
                OutputLength = 100
            };

            var exception = Assert.Throws<System.ArgumentException>(() => sampler.Sample());
            Assert.AreEqual(nameof(SNIS.NumberOfChains), exception.ParamName);
        }

        /// <summary>
        /// This test verifies that the resampling sort is stable — draws tied at the same fitness
        /// keep their original draw order — and that an identically-seeded run reproduces the
        /// output exactly, element for element.
        /// </summary>
        /// <remarks>
        /// The fixture yields two tied-fitness runs: half the draws sit at a log-likelihood of
        /// exactly negative infinity and the rest at exactly zero, so the sorted chain is fully
        /// determined by the tie-breaking rule. The oracle reproduces the sampler's parameter
        /// draws from the same master PRNG seed and requires the sorted chain to be the
        /// negative-infinity run followed by the finite run, each in original draw order. An
        /// unstable sort would reorder the tied draws and change which samples the resampling
        /// plotting positions select.
        /// </remarks>
        [TestMethod]
        public void Test_SNIS_TiedFitnessDraws_KeepDrawOrderAndReproduceExactly()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(0d, 1d) };
            double logLH(double[] x) => x[0] < 0.5d ? 0d : double.NegativeInfinity;

            var first = new SNIS(priors, logLH) { Iterations = 100, OutputLength = 100, PRNGSeed = 12345 };
            var second = new SNIS(priors, logLH) { Iterations = 100, OutputLength = 100, PRNGSeed = 12345 };
            first.Sample();
            second.Sample();

            Assert.HasCount(first.Output[0].Count, second.Output[0]);
            for (int i = 0; i < first.Output[0].Count; i++)
            {
                Assert.AreEqual(first.Output[0][i].Values[0], second.Output[0][i].Values[0], 0d);
            }

            // Reproduce the sampler's parameter draws: the master PRNG is seeded with PRNGSeed and
            // each draw is the prior inverse CDF of its uniform matrix, in draw-index order.
            var prior = new Uniform(0d, 1d);
            var rnds = new MersenneTwister(12345).NextDoubles(100, 1);
            var draws = new double[100];
            for (int i = 0; i < draws.Length; i++)
                draws[i] = prior.InverseCDF(rnds[i, 0]);

            // The stable ascending sort on Fitness places the negative-infinity run first and the
            // zero-fitness run second, each preserving its original draw order.
            var expected = draws.Where(d => !(d < 0.5d)).Concat(draws.Where(d => d < 0.5d)).ToArray();
            Assert.IsTrue(draws.Any(d => !(d < 0.5d)));
            Assert.IsTrue(draws.Any(d => d < 0.5d));
            for (int i = 0; i < expected.Length; i++)
            {
                Assert.AreEqual(expected[i], first.MarkovChains[0][i].Values[0], 0d);
            }
        }

    }
}
