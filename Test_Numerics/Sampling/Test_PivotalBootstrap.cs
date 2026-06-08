using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Functions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;
using System;
using System.Linq;

namespace Sampling
{
    /// <summary>
    /// Unit tests for the bias-corrected pivotal bootstrap.
    /// </summary>
    [TestClass]
    public class Test_PivotalBootstrap
    {
        /// <summary>
        /// Verifies the pivotal transform uses the replicate covariance to standardize
        /// and the parent covariance to re-inflate.
        /// </summary>
        [TestMethod]
        public void Transform_IdentityLink_UsesBothCovariances()
        {
            var parent = Fit(new[] { 10d, 20d }, new double[,] { { 4d, 0d }, { 0d, 9d } });
            var raw = Fit(new[] { 8d, 17d }, new double[,] { { 1d, 0d }, { 0d, 1d } });

            var options = new PivotalBootstrapTransformOptions
            {
                RegularizeCovariances = false
            };

            PivotalBootstrapResults result = PivotalBootstrap.Transform(parent, new[] { raw }, options);

            Assert.HasCount(1, result.PivotalParameterSets);
            Assert.AreEqual(14d, result.PivotalParameterSets[0].Values[0], 1e-12);
            Assert.AreEqual(29d, result.PivotalParameterSets[0].Values[1], 1e-12);
        }

        /// <summary>
        /// Verifies callers have full control over the links through the link factory.
        /// </summary>
        [TestMethod]
        public void Transform_CustomLinkFactory_ControlsLinkDefinition()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var raw = Fit(new[] { 8d }, new double[,] { { 1d } });
            bool factoryWasCalled = false;

            var options = new PivotalBootstrapTransformOptions
            {
                RegularizeCovariances = false,
                LinkFactory = context =>
                {
                    factoryWasCalled = true;
                    Assert.AreEqual(1, context.ParameterCount);
                    CollectionAssert.AreEqual(new[] { 8d }, context.GetRawParameterValues(0));
                    return new[]
                    {
                        new PivotalBootstrapLink(
                            x => 2d * x,
                            eta => eta / 2d,
                            x => 2d,
                            "Double")
                    };
                }
            };

            PivotalBootstrapResults result = PivotalBootstrap.Transform(parent, new[] { raw }, options);

            Assert.IsTrue(factoryWasCalled);
            Assert.AreEqual("Double", result.Links[0].Name);
            Assert.AreEqual(14d, result.PivotalParameterSets[0].Values[0], 1e-12);
        }

        /// <summary>
        /// Verifies the normal location-scale pivotal bootstrap recovers objective-Bayes marginals
        /// and quantile intervals to Monte-Carlo tolerance.
        /// </summary>
        [TestMethod]
        public void Run_NormalLocationScale_MatchesObjectiveBayesMarginalsAndQuantileIntervals()
        {
            const int sampleSize = 100;
            const int replicates = 8000;
            var parentDistribution = new Normal(3d, 0.7d);
            var parent = new PivotalBootstrapFit(
                new ParameterSet(parentDistribution.GetParameters, double.NaN),
                new Matrix(parentDistribution.ParameterCovariance(sampleSize, ParameterEstimationMethod.MaximumLikelihood)));
            double[] originalData = new double[sampleSize];
            double[] probabilities = { 0.5d, 0.95d };

            var options = new PivotalBootstrapOptions<double[]>
            {
                Replicates = replicates,
                PRNGSeed = 12345,
                ResampleFunction = (data, fit, rng) =>
                {
                    var distribution = new Normal(fit.Parameters.Values[0], fit.Parameters.Values[1]);
                    return distribution.GenerateRandomValues(sampleSize, rng.Next());
                },
                FitFunction = sample => FitNormal(sample),
                LinkFactory = context => new[]
                {
                    PivotalBootstrapLink.FromLinkFunction(new IdentityLink()),
                    PivotalBootstrapLink.FromLinkFunction(new LogLink())
                },
                ParameterValidator = values => values[1] > 0d,
                StatisticFunction = ps =>
                {
                    var distribution = new Normal(ps.Values[0], ps.Values[1]);
                    return probabilities.Select(distribution.InverseCDF).ToArray();
                }
            };

            PivotalBootstrapResults result = PivotalBootstrap.Run(originalData, parent, options);

            Assert.IsGreaterThan(0.99d * replicates, result.PivotalParameterSets.Length);
            double[] mus = result.PivotalParameterSets.Select(ps => ps.Values[0]).ToArray();
            double[] sigmas = result.PivotalParameterSets.Select(ps => ps.Values[1]).ToArray();
            var muPosterior = new Normal(parentDistribution.Mu, parentDistribution.Sigma / Math.Sqrt(sampleSize));
            var sigmaSquaredPosterior = new InverseChiSquared(sampleSize - 1, parentDistribution.Sigma * parentDistribution.Sigma);

            Array.Sort(mus);
            double[] sigmaSquares = sigmas.Select(sigma => sigma * sigma).ToArray();
            Array.Sort(sigmaSquares);

            double muKs = GoodnessOfFit.KolmogorovSmirnov(mus, muPosterior);
            double sigmaKs = GoodnessOfFit.KolmogorovSmirnov(sigmaSquares, sigmaSquaredPosterior);

            Assert.IsLessThan(0.04d, muKs, $"Mu marginal KS distance was {muKs}.");
            Assert.IsLessThan(0.08d, sigmaKs, $"Sigma marginal KS distance was {sigmaKs}.");

            BootstrapStatisticResult[] quantileIntervals = result.GetPivotalStatisticConfidenceIntervals(0.1d);
            double[,] objectiveBayesIntervals = parentDistribution.MonteCarloConfidenceIntervals(
                sampleSize,
                20000,
                probabilities,
                new[] { 0.05d, 0.95d });

            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(objectiveBayesIntervals[i, 0], quantileIntervals[i].LowerCI, 0.10d * parentDistribution.Sigma);
                Assert.AreEqual(objectiveBayesIntervals[i, 1], quantileIntervals[i].UpperCI, 0.10d * parentDistribution.Sigma);
            }
        }

        private static PivotalBootstrapFit Fit(double[] values, double[,] covariance)
        {
            return new PivotalBootstrapFit(new ParameterSet(values, double.NaN), new Matrix(covariance));
        }

        private static PivotalBootstrapFit FitNormal(double[] sample)
        {
            var distribution = new Normal();
            ((IEstimation)distribution).Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            if (!distribution.ParametersValid)
                throw new InvalidOperationException("The normal fit produced invalid parameters.");

            return new PivotalBootstrapFit(
                new ParameterSet(distribution.GetParameters, double.NaN),
                new Matrix(distribution.ParameterCovariance(sample.Length, ParameterEstimationMethod.MaximumLikelihood)));
        }

    }
}
