using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
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
    /// Unit tests for the covariance-aware pivotal bootstrap workflow on <see cref="Bootstrap{TData}"/>.
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
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;

            boot.TransformPivotalBootstrap(new[] { raw });

            Assert.HasCount(1, boot.BootstrapParameterSets);
            Assert.AreEqual(14d, boot.BootstrapParameterSets[0].Values[0], 1e-12);
            Assert.AreEqual(29d, boot.BootstrapParameterSets[0].Values[1], 1e-12);
        }

        /// <summary>
        /// Verifies a full covariance matrix is used in the two-covariance transform.
        /// </summary>
        [TestMethod]
        public void Transform_FullCovariance_UsesCholeskyAlgebra()
        {
            var parent = Fit(new[] { 10d, 20d }, new double[,] { { 4d, 1.2d }, { 1.2d, 9d } });
            var raw = Fit(new[] { 8d, 17d }, new double[,] { { 1d, 0.25d }, { 0.25d, 2d } });
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;

            boot.TransformPivotalBootstrap(new[] { raw });

            double[] expected = ExpectedIdentityPivotalValues(parent, raw);
            Assert.AreEqual(expected[0], boot.BootstrapParameterSets[0].Values[0], 1e-12);
            Assert.AreEqual(expected[1], boot.BootstrapParameterSets[0].Values[1], 1e-12);
        }

        /// <summary>
        /// Verifies callers have full control over pivotal links through the link factory.
        /// </summary>
        [TestMethod]
        public void Transform_CustomLinkFactory_ControlsLinkDefinition()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var raw = Fit(new[] { 8d }, new double[,] { { 1d } });
            var boot = CreatePivotalBootstrap(parent);
            bool factoryWasCalled = false;

            boot.RegularizePivotalCovariances = false;
            boot.PivotalLinkFactory = context =>
            {
                factoryWasCalled = true;
                Assert.AreEqual(1, context.ParameterCount);
                CollectionAssert.AreEqual(new[] { 8d }, context.GetRawParameterValues(0));
                return new ILinkFunction[] { new LinearScaleLink(2d) };
            };

            boot.TransformPivotalBootstrap(new[] { raw });

            Assert.IsTrue(factoryWasCalled);
            Assert.IsInstanceOfType(boot.PivotalLinks[0], typeof(LinearScaleLink));
            Assert.AreEqual(14d, boot.BootstrapParameterSets[0].Values[0], 1e-12);
        }

        /// <summary>
        /// Verifies standard Numerics links can be mixed in a pivotal transformation.
        /// </summary>
        [TestMethod]
        public void Transform_LogFisherZAndYeoJohnsonLinks_RetainsFiniteDraws()
        {
            var parent = Fit(
                new[] { 1d, 2d, 0.25d, -0.5d },
                new double[,]
                {
                    { 0.04d, 0d, 0d, 0d },
                    { 0d, 0.09d, 0d, 0d },
                    { 0d, 0d, 0.01d, 0d },
                    { 0d, 0d, 0d, 0.16d }
                });
            var rawFits = new[]
            {
                Fit(
                    new[] { 0.9d, 1.8d, 0.10d, -0.7d },
                    new double[,]
                    {
                        { 0.01d, 0d, 0d, 0d },
                        { 0d, 0.04d, 0d, 0d },
                        { 0d, 0d, 0.0025d, 0d },
                        { 0d, 0d, 0d, 0.09d }
                    }),
                Fit(
                    new[] { 1.1d, 2.1d, 0.35d, -0.3d },
                    new double[,]
                    {
                        { 0.01d, 0d, 0d, 0d },
                        { 0d, 0.04d, 0d, 0d },
                        { 0d, 0d, 0.0025d, 0d },
                        { 0d, 0d, 0d, 0.09d }
                    })
            };
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;
            boot.PivotalLinkFactory = context => new ILinkFunction[]
            {
                null,
                new LogLink(),
                new FisherZLink(),
                new YeoJohnsonLink(context.GetRawParameterValues(3))
            };
            boot.PivotalParameterValidator = values => values[1] > 0d && values[2] > -1d && values[2] < 1d;

            boot.TransformPivotalBootstrap(rawFits);

            Assert.HasCount(2, boot.BootstrapParameterSets);
            Assert.IsNull(boot.PivotalLinks[0]);
            Assert.IsInstanceOfType(boot.PivotalLinks[1], typeof(LogLink));
            Assert.IsInstanceOfType(boot.PivotalLinks[2], typeof(FisherZLink));
            Assert.IsInstanceOfType(boot.PivotalLinks[3], typeof(YeoJohnsonLink));
            foreach (double value in boot.BootstrapParameterSets[0].Values)
                Assert.IsTrue(Tools.IsFinite(value));
        }

        /// <summary>
        /// Verifies raw and pivotal confidence intervals are exposed separately after pivotal transformation.
        /// </summary>
        [TestMethod]
        public void Transform_WithStatisticFunction_StoresRawAndPivotalStatistics()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var rawFits = new[]
            {
                Fit(new[] { 8d }, new double[,] { { 1d } }),
                Fit(new[] { 9d }, new double[,] { { 1d } }),
                Fit(new[] { 10d }, new double[,] { { 1d } }),
                Fit(new[] { 11d }, new double[,] { { 1d } }),
                Fit(new[] { 12d }, new double[,] { { 1d } })
            };
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;
            boot.StatisticFunction = ps => new[] { ps.Values[0] * 2d };

            boot.TransformPivotalBootstrap(rawFits);
            BootstrapResults pivotal = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile, 0.2d);
            BootstrapResults raw = boot.GetRawPivotalConfidenceIntervals(0.2d);

            Assert.HasCount(1, pivotal.StatisticResults);
            Assert.HasCount(1, raw.StatisticResults);
            Assert.AreEqual(20d, pivotal.StatisticResults[0].PopulationEstimate, 1e-12);
            Assert.AreEqual(20d, raw.StatisticResults[0].PopulationEstimate, 1e-12);
            Assert.IsGreaterThan(pivotal.ParameterResults[0].LowerCI, raw.ParameterResults[0].LowerCI);
        }

        /// <summary>
        /// Verifies pivotal parameter ensembles can be generated without a statistic function.
        /// </summary>
        [TestMethod]
        public void Transform_WithoutStatisticFunction_ProducesParameterIntervalsOnly()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var raw = Fit(new[] { 8d }, new double[,] { { 1d } });
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;

            boot.TransformPivotalBootstrap(new[] { raw });
            BootstrapResults results = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile);

            Assert.HasCount(1, results.ParameterResults);
            Assert.HasCount(0, results.StatisticResults);
            Assert.AreEqual(0, boot.BootstrapStatistics.GetLength(1));
        }

        /// <summary>
        /// Verifies invalid pivotal draws can be dropped.
        /// </summary>
        [TestMethod]
        public void Transform_InvalidDrawPolicyDrop_DropsInvalidPivotalDraws()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var raw = Fit(new[] { 8d }, new double[,] { { 1d } });
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;
            boot.PivotalParameterValidator = _ => false;
            boot.PivotalInvalidDrawPolicy = PivotalBootstrapInvalidDrawPolicy.Drop;

            boot.TransformPivotalBootstrap(new[] { raw });

            Assert.HasCount(0, boot.BootstrapParameterSets);
            Assert.AreEqual(1, boot.PivotalDiagnostics!.InvalidPivotalReplicates);
            Assert.AreEqual(0, boot.PivotalDiagnostics.RetainedPivotalReplicates);
        }

        /// <summary>
        /// Verifies invalid pivotal draws can fall back to raw bootstrap fits.
        /// </summary>
        [TestMethod]
        public void Transform_InvalidDrawPolicyUseRaw_UsesRawParameters()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var raw = Fit(new[] { 8d }, new double[,] { { 1d } });
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;
            boot.PivotalParameterValidator = _ => false;
            boot.PivotalInvalidDrawPolicy = PivotalBootstrapInvalidDrawPolicy.UseRaw;

            boot.TransformPivotalBootstrap(new[] { raw });

            Assert.HasCount(1, boot.BootstrapParameterSets);
            Assert.AreEqual(8d, boot.BootstrapParameterSets[0].Values[0], 1e-12);
        }

        /// <summary>
        /// Verifies invalid pivotal draws can fall back to the parent fit.
        /// </summary>
        [TestMethod]
        public void Transform_InvalidDrawPolicyUseParent_UsesParentParameters()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var raw = Fit(new[] { 8d }, new double[,] { { 1d } });
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;
            boot.PivotalParameterValidator = _ => false;
            boot.PivotalInvalidDrawPolicy = PivotalBootstrapInvalidDrawPolicy.UseParent;

            boot.TransformPivotalBootstrap(new[] { raw });

            Assert.HasCount(1, boot.BootstrapParameterSets);
            Assert.AreEqual(10d, boot.BootstrapParameterSets[0].Values[0], 1e-12);
        }

        /// <summary>
        /// Verifies raw replicate filtering occurs before pivotal transformation and is reported in diagnostics.
        /// </summary>
        [TestMethod]
        public void Transform_ReplicateFilter_RejectsRawFitsAndUpdatesDiagnostics()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var rawFits = new[]
            {
                Fit(new[] { 8d }, new double[,] { { 1d } }),
                Fit(new[] { 9d }, new double[,] { { 1d } }),
                Fit(new[] { 12d }, new double[,] { { 1d } })
            };
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;
            boot.PivotalReplicateFilter = fit => fit.Parameters.Values[0] >= 9d;

            boot.TransformPivotalBootstrap(rawFits);

            Assert.HasCount(2, boot.RawBootstrapParameterSets);
            Assert.HasCount(2, boot.BootstrapParameterSets);
            Assert.AreEqual(1, boot.PivotalDiagnostics!.RejectedRawReplicates);
        }

        /// <summary>
        /// Verifies the pivotal run method uses the covariance-aware fit delegate, not the regular fit delegate.
        /// </summary>
        [TestMethod]
        public void RunPivotalBootstrap_IgnoresRegularFitFunction()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var samples = new[] { 8d, 9d, 10d, 11d };
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;
            boot.Replicates = samples.Length;
            boot.MaxRetries = 1;
            boot.ResampleFunction = (data, ps, rng) => data;
            boot.FitFunction = _ => throw new InvalidOperationException("Regular fit should not be used by pivotal bootstrap.");
            int index = -1;
            boot.FitWithCovarianceFunction = _ =>
            {
                int next = System.Threading.Interlocked.Increment(ref index);
                return Fit(new[] { samples[next] }, new double[,] { { 1d } });
            };

            boot.RunPivotalBootstrap();

            Assert.HasCount(samples.Length, boot.RawBootstrapFits);
            Assert.HasCount(samples.Length, boot.BootstrapParameterSets);
        }

        /// <summary>
        /// Verifies regular bootstrap ignores pivotal-only properties and delegates.
        /// </summary>
        [TestMethod]
        public void Run_RegularBootstrap_IgnoresPivotalOnlyProperties()
        {
            var boot = new Bootstrap<double[]>(new[] { 1d }, new ParameterSet(new[] { 10d }, double.NaN));
            boot.Replicates = 3;
            boot.MaxRetries = 1;
            boot.ResampleFunction = (data, ps, rng) => data;
            boot.FitFunction = _ => new ParameterSet(new[] { 11d }, double.NaN);
            boot.StatisticFunction = ps => new[] { ps.Values[0] };
            boot.FitWithCovarianceFunction = _ => throw new InvalidOperationException("Pivotal fit should not be used by regular bootstrap.");
            boot.PivotalLinkFactory = _ => throw new InvalidOperationException("Pivotal links should not be used by regular bootstrap.");
            boot.PivotalReplicateFilter = _ => throw new InvalidOperationException("Pivotal filters should not be used by regular bootstrap.");
            boot.OriginalCovariance = new Matrix(new[,] { { 1d } });

            boot.Run();

            Assert.HasCount(3, boot.BootstrapParameterSets);
            Assert.HasCount(0, boot.RawBootstrapParameterSets);
        }

        /// <summary>
        /// Verifies missing covariance-aware fit delegate produces a clear failure.
        /// </summary>
        [TestMethod]
        public void RunPivotalBootstrap_WithoutFitWithCovarianceFunction_Throws()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var boot = CreatePivotalBootstrap(parent);
            boot.ResampleFunction = (data, ps, rng) => data;
            boot.FitWithCovarianceFunction = null;

            var exception = Assert.Throws<InvalidOperationException>(() => boot.RunPivotalBootstrap());
            StringAssert.Contains(exception.Message, nameof(Bootstrap<double[]>.FitWithCovarianceFunction));
        }

        /// <summary>
        /// Verifies missing original covariance produces a clear failure.
        /// </summary>
        [TestMethod]
        public void RunPivotalBootstrap_WithoutOriginalCovariance_Throws()
        {
            var boot = new Bootstrap<double[]>(new[] { 1d }, new ParameterSet(new[] { 10d }, double.NaN));
            boot.ResampleFunction = (data, ps, rng) => data;
            boot.FitWithCovarianceFunction = _ => Fit(new[] { 10d }, new double[,] { { 1d } });

            var exception = Assert.Throws<InvalidOperationException>(() => boot.RunPivotalBootstrap());
            StringAssert.Contains(exception.Message, nameof(Bootstrap<double[]>.OriginalCovariance));
        }

        /// <summary>
        /// Verifies incompatible confidence interval methods are rejected after a pivotal run.
        /// </summary>
        [TestMethod]
        public void GetConfidenceIntervals_BCaAfterPivotalRun_Throws()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var raw = Fit(new[] { 8d }, new double[,] { { 1d } });
            var boot = CreatePivotalBootstrap(parent);
            boot.RegularizePivotalCovariances = false;
            boot.TransformPivotalBootstrap(new[] { raw });

            var exception = Assert.Throws<InvalidOperationException>(() => boot.GetConfidenceIntervals(BootstrapCIMethod.BCa));
            StringAssert.Contains(exception.Message, "Only percentile");
        }

        /// <summary>
        /// Verifies raw pivotal confidence intervals require a pivotal run first.
        /// </summary>
        [TestMethod]
        public void GetRawPivotalConfidenceIntervals_BeforePivotalRun_Throws()
        {
            var parent = Fit(new[] { 10d }, new double[,] { { 4d } });
            var boot = CreatePivotalBootstrap(parent);

            Assert.Throws<InvalidOperationException>(() => boot.GetRawPivotalConfidenceIntervals());
        }

        /// <summary>
        /// Verifies repeated pivotal runs with the same seed retain reproducibility.
        /// </summary>
        [TestMethod]
        public void RunPivotalBootstrap_WithSameSeed_IsReproducible()
        {
            Bootstrap<double[]> first = CreateSeededNormalPivotalBootstrap();
            Bootstrap<double[]> second = CreateSeededNormalPivotalBootstrap();

            first.RunPivotalBootstrap();
            second.RunPivotalBootstrap();

            Assert.HasCount(first.BootstrapParameterSets.Length, second.BootstrapParameterSets);
            for (int i = 0; i < first.BootstrapParameterSets.Length; i++)
            {
                Assert.AreEqual(first.BootstrapParameterSets[i].Values[0], second.BootstrapParameterSets[i].Values[0], 1e-12);
                Assert.AreEqual(first.BootstrapParameterSets[i].Values[1], second.BootstrapParameterSets[i].Values[1], 1e-12);
            }
        }

        /// <summary>
        /// Verifies the normal location-scale pivotal bootstrap produces plausible objective-Bayes-style quantile intervals.
        /// </summary>
        [TestMethod]
        public void Run_NormalLocationScale_MatchesObjectiveBayesQuantileIntervals()
        {
            const int sampleSize = 100;
            const int replicates = 2500;
            var parentDistribution = new Normal(3d, 0.7d);
            var parent = new BootstrapFit(
                new ParameterSet(parentDistribution.GetParameters, double.NaN),
                new Matrix(parentDistribution.ParameterCovariance(sampleSize, ParameterEstimationMethod.MaximumLikelihood)));
            double[] originalData = new double[sampleSize];
            double[] probabilities = { 0.5d, 0.95d };
            var boot = new Bootstrap<double[]>(originalData, parent);

            boot.Replicates = replicates;
            boot.PRNGSeed = 12345;
            boot.ResampleFunction = (data, ps, rng) =>
            {
                var distribution = new Normal(ps.Values[0], ps.Values[1]);
                return distribution.GenerateRandomValues(sampleSize, rng.Next());
            };
            boot.FitWithCovarianceFunction = FitNormal;
            boot.PivotalLinkFactory = _ => new ILinkFunction[] { new IdentityLink(), new LogLink() };
            boot.PivotalParameterValidator = values => values[1] > 0d;
            boot.StatisticFunction = ps =>
            {
                var distribution = new Normal(ps.Values[0], ps.Values[1]);
                return probabilities.Select(distribution.InverseCDF).ToArray();
            };

            boot.RunPivotalBootstrap();
            BootstrapResults quantileIntervals = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile, 0.1d);
            double[,] objectiveBayesIntervals = parentDistribution.MonteCarloConfidenceIntervals(
                sampleSize,
                12000,
                probabilities,
                new[] { 0.05d, 0.95d });

            Assert.IsGreaterThan(0.98d * replicates, boot.BootstrapParameterSets.Length);
            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(objectiveBayesIntervals[i, 0], quantileIntervals.StatisticResults[i].LowerCI, 0.12d * parentDistribution.Sigma);
                Assert.AreEqual(objectiveBayesIntervals[i, 1], quantileIntervals.StatisticResults[i].UpperCI, 0.12d * parentDistribution.Sigma);
            }
        }

        /// <summary>
        /// Creates a covariance-aware bootstrap configured with a fixed parent fit.
        /// </summary>
        /// <param name="parent">The parent fit.</param>
        /// <returns>A bootstrap configured for pivotal transformations.</returns>
        private static Bootstrap<double[]> CreatePivotalBootstrap(BootstrapFit parent)
        {
            return new Bootstrap<double[]>(Array.Empty<double>(), parent);
        }

        /// <summary>
        /// Creates a covariance-aware fit from parameter values and covariance entries.
        /// </summary>
        /// <param name="values">The fitted parameter values.</param>
        /// <param name="covariance">The covariance matrix.</param>
        /// <returns>The covariance-aware fit.</returns>
        private static BootstrapFit Fit(double[] values, double[,] covariance)
        {
            return new BootstrapFit(new ParameterSet(values, double.NaN), new Matrix(covariance));
        }

        /// <summary>
        /// Fits a normal distribution by maximum likelihood and returns parameter covariance.
        /// </summary>
        /// <param name="sample">The sample to fit.</param>
        /// <returns>The covariance-aware normal fit.</returns>
        private static BootstrapFit FitNormal(double[] sample)
        {
            var distribution = new Normal();
            ((IEstimation)distribution).Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            if (!distribution.ParametersValid)
                throw new InvalidOperationException("The normal fit produced invalid parameters.");

            return new BootstrapFit(
                new ParameterSet(distribution.GetParameters, double.NaN),
                new Matrix(distribution.ParameterCovariance(sample.Length, ParameterEstimationMethod.MaximumLikelihood)));
        }

        /// <summary>
        /// Creates a reproducible normal location-scale pivotal bootstrap for reproducibility checks.
        /// </summary>
        /// <returns>A configured pivotal bootstrap instance.</returns>
        private static Bootstrap<double[]> CreateSeededNormalPivotalBootstrap()
        {
            const int sampleSize = 25;
            var distribution = new Normal(3d, 0.7d);
            var parent = new BootstrapFit(
                new ParameterSet(distribution.GetParameters, double.NaN),
                new Matrix(distribution.ParameterCovariance(sampleSize, ParameterEstimationMethod.MaximumLikelihood)));
            var boot = new Bootstrap<double[]>(new double[sampleSize], parent);
            boot.Replicates = 50;
            boot.PRNGSeed = 8675309;
            boot.ResampleFunction = (data, ps, rng) => new Normal(ps.Values[0], ps.Values[1]).GenerateRandomValues(sampleSize, rng.Next());
            boot.FitWithCovarianceFunction = FitNormal;
            boot.PivotalLinkFactory = _ => new ILinkFunction[] { null, new LogLink() };
            boot.PivotalParameterValidator = values => values[1] > 0d;
            return boot;
        }

        /// <summary>
        /// Computes the expected identity-link pivotal value for a single raw fit.
        /// </summary>
        /// <param name="parent">The parent fit.</param>
        /// <param name="raw">The raw bootstrap fit.</param>
        /// <returns>The expected pivotal parameter values.</returns>
        private static double[] ExpectedIdentityPivotalValues(BootstrapFit parent, BootstrapFit raw)
        {
            var parentCholesky = new CholeskyDecomposition(parent.Covariance);
            var rawCholesky = new CholeskyDecomposition(raw.Covariance);
            var difference = new double[parent.ParameterCount];
            for (int i = 0; i < difference.Length; i++)
                difference[i] = parent.Parameters.Values[i] - raw.Parameters.Values[i];

            double[] z = rawCholesky.Forward(new Vector(difference)).ToArray();
            double[] reinflated = parentCholesky.L * z;
            var expected = new double[difference.Length];
            for (int i = 0; i < expected.Length; i++)
                expected[i] = parent.Parameters.Values[i] + reinflated[i];
            return expected;
        }

        /// <summary>
        /// A linear scaling link used to verify the pivotal link factory accepts custom <see cref="ILinkFunction"/> instances.
        /// </summary>
        [Serializable]
        private sealed class LinearScaleLink : ILinkFunction
        {
            private readonly double _scale;

            /// <summary>
            /// Initializes a new instance of the <see cref="LinearScaleLink"/> class.
            /// </summary>
            /// <param name="scale">The nonzero linear scale.</param>
            public LinearScaleLink(double scale)
            {
                _scale = scale;
            }

            /// <inheritdoc/>
            public double Link(double x)
            {
                return _scale * x;
            }

            /// <inheritdoc/>
            public double InverseLink(double eta)
            {
                return eta / _scale;
            }

            /// <inheritdoc/>
            public double DLink(double x)
            {
                return _scale;
            }

            /// <inheritdoc/>
            public System.Xml.Linq.XElement ToXElement()
            {
                return new System.Xml.Linq.XElement(nameof(LinearScaleLink));
            }
        }
    }
}
