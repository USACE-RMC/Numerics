using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit test for the Gibbs sampler. 
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
    public class Test_Gibbs
    {
        /// <summary>
        /// This test compares the results obtained using Gibbs with those from the 'rstan' package. 
        /// </summary>
        [TestMethod]
        public void Test_Gibbs_NormalDist_RStan()
        {

            // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
            // Table 5.1.1 Tippecanoe River Near Delphi, Indiana (Station 43) Data
            double[] sample = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

            // Create non-informative priors
            int n = sample.Length;
            var mu = Statistics.Mean(sample);
            double mu0 = 0, sigma0 = 5E5;
            // Preserve the reference model's limiting non-informative inverse-gamma prior.
            double variancePriorShape = 0d, variancePriorScale = 0d;

            // Use proper distributions solely to initialize the sampler state.
            double initializationShape = 2d, initializationScale = 0.001d;
            var muInitializationPrior = new Normal(mu0, sigma0);
            var sigmaInitializationPrior = new InverseGamma(initializationScale, initializationShape);
            var conditionalMean = new Normal();
            var conditionalVariance = new InverseGamma();
            var priors = new List<IUnivariateDistribution> { muInitializationPrior, sigmaInitializationPrior };

            // Create log-likelihood function
            double logLH(double[] x)
            {
                var dist = new Normal(x[0], x[1]);
                return dist.LogLikelihood(sample);
            }

            // Create proposal function
            double[] proposal(double[] x, Random random)
            {
                // Sample the conditional mean given the current standard deviation.
                var meanParameters = ConditionalMeanParameters(n, mu, x[1], mu0, sigma0);
                conditionalMean.SetParameters(meanParameters.Mean, meanParameters.StandardDeviation);
                double mup = conditionalMean.InverseCDF(random.NextDouble());

                // Sample the conditional variance, then return its square root as sigma.
                var varianceParameters = ConditionalVarianceParameters(sample, mup, variancePriorShape, variancePriorScale);
                conditionalVariance.SetParameters(new[] { varianceParameters.Scale, varianceParameters.Shape });
                double sig2p = conditionalVariance.InverseCDF(random.NextDouble());

                // return proposal vector
                return new double[] { mup, Math.Sqrt(sig2p) };
            }

            // Create and run sampler
            var sampler = new Gibbs(priors, logLH, proposal);
            sampler.Sample();
            var results = new MCMCResults(sampler);
            Assert.AreEqual(mu0, muInitializationPrior.Mu);
            Assert.AreEqual(sigma0, muInitializationPrior.Sigma);
            Assert.AreEqual(initializationScale, sigmaInitializationPrior.Beta);
            Assert.AreEqual(initializationShape, sigmaInitializationPrior.Alpha);

            /* Below are the results from 'rstan' using comparable MCMC settings:
            *            mean se_mean     sd       5%      50%      95% n_eff Rhat
            *  mu    12663.69    7.10 706.60 11488.50 12671.08 13801.45  9897    1
            *  sigma  4844.09    5.22 519.08  4077.80  4796.63  5771.81  9880    1
            *  lp__   -466.13    0.01   1.03  -468.17  -465.81  -465.15  9958    1
            * 
            *  Since MCMC methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these
            *  comparisons verify that the results remain within 5% of the 'rstan' results.
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
        /// Verifies the Normal and inverse-gamma full-conditionals with an informative, nonzero prior.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalParameters_InformativePrior()
        {
            var meanParameters = ConditionalMeanParameters(4, 10d, 2d, 2d, 3d);
            Assert.AreEqual(9.2d, meanParameters.Mean, 1E-12);
            Assert.AreEqual(Math.Sqrt(0.9d), meanParameters.StandardDeviation, 1E-12);

            var varianceParameters = ConditionalVarianceParameters(new[] { 8d, 10d, 12d, 14d }, 10d, 2d, 0.5d);
            Assert.AreEqual(12.5d, varianceParameters.Scale, 1E-12);
            Assert.AreEqual(4d, varianceParameters.Shape, 1E-12);
        }

        /// <summary>
        /// Computes the Normal full-conditional parameters for the population mean when the
        /// likelihood variance is fixed at its current Gibbs state.
        /// </summary>
        /// <param name="sampleSize">The number of observations.</param>
        /// <param name="sampleMean">The arithmetic mean of the observations.</param>
        /// <param name="currentStandardDeviation">The current likelihood standard deviation.</param>
        /// <param name="priorMean">The Normal prior mean.</param>
        /// <param name="priorStandardDeviation">The Normal prior standard deviation.</param>
        /// <returns>The conditional Normal mean and standard deviation.</returns>
        private static (double Mean, double StandardDeviation) ConditionalMeanParameters(
            int sampleSize,
            double sampleMean,
            double currentStandardDeviation,
            double priorMean,
            double priorStandardDeviation)
        {
            double likelihoodVariance = currentStandardDeviation * currentStandardDeviation;
            double priorVariance = priorStandardDeviation * priorStandardDeviation;
            double posteriorVariance = 1d / (sampleSize / likelihoodVariance + 1d / priorVariance);
            double posteriorMean = posteriorVariance *
                (sampleSize * sampleMean / likelihoodVariance + priorMean / priorVariance);
            return (posteriorMean, Math.Sqrt(posteriorVariance));
        }

        /// <summary>
        /// Computes the inverse-gamma full-conditional parameters for the likelihood variance.
        /// </summary>
        /// <param name="sample">The observed sample.</param>
        /// <param name="conditionalMean">The mean drawn during the current Gibbs iteration.</param>
        /// <param name="priorShape">The inverse-gamma prior shape.</param>
        /// <param name="priorScale">The inverse-gamma prior scale.</param>
        /// <returns>The conditional inverse-gamma scale and shape, in the order expected by <see cref="InverseGamma"/>.</returns>
        private static (double Scale, double Shape) ConditionalVarianceParameters(
            IList<double> sample,
            double conditionalMean,
            double priorShape,
            double priorScale)
        {
            double sumOfSquaredErrors = 0d;
            for (int i = 0; i < sample.Count; i++)
            {
                double residual = sample[i] - conditionalMean;
                sumOfSquaredErrors += residual * residual;
            }
            return (priorScale + sumOfSquaredErrors / 2d, priorShape + sample.Count / 2d);
        }
      
    }
}
