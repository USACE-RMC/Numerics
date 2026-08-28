using System.Reflection;
using Numerics;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit tests for adaptive-random-walk covariance updates and NUTS diagnostics.
    /// </summary>
    [TestClass]
    public class Test_MCMCSamplerDiagnostics
    {
        /// <summary>
        /// Confirms that rejected transitions before the warmup boundary advance the
        /// ARWMH running covariance sample count.
        /// </summary>
        /// <remarks>
        /// Adaptive Metropolis covariance is defined from the realized chain history, so every
        /// realized transition, accepted or rejected, contributes one observation to the running
        /// covariance; a rejected proposal contributes the repeated retained state.
        /// </remarks>
        [TestMethod]
        public void ARWMH_RejectedWarmupTransitionsEnterCovariance()
        {
            InspectableARWMH sampler = CreateRejectingArwmh(warmupIterations: 50);
            ParameterSet state = new ParameterSet(new[] { 0d, 0d }, 0d);
            sampler.InitializeAt(state);

            for (int i = 0; i < 12; i++)
                state = sampler.Iterate(state);

            Assert.AreEqual(12, sampler.SampleCount[0]);
            Assert.AreEqual(0, sampler.AcceptCount[0]);
            Assert.AreEqual(
                12,
                GetArwmhCovarianceSampleCount(sampler),
                "Every rejected warmup transition must contribute its repeated retained state.");
        }

        /// <summary>
        /// Confirms that ARWMH continues its original continual-adaptation schedule after
        /// the configured warmup boundary.
        /// </summary>
        /// <remarks>
        /// Continued Adaptive Metropolis updating is consistent with the original
        /// Haario-Saksman-Tamminen construction. Together with the pre-warmup test above, this
        /// pins the schedule on both sides of the boundary: every realized transition enters the
        /// running covariance before the boundary and after it alike.
        /// </remarks>
        [TestMethod]
        public void ARWMH_RejectedPostWarmupTransitionsContinueEnteringCovariance()
        {
            InspectableARWMH sampler = CreateRejectingArwmh(warmupIterations: 5);
            ParameterSet state = new ParameterSet(new[] { 0d, 0d }, 0d);
            sampler.InitializeAt(state);

            for (int i = 0; i < 12; i++)
                state = sampler.Iterate(state);

            Assert.AreEqual(12, sampler.SampleCount[0]);
            Assert.AreEqual(0, sampler.AcceptCount[0]);
            Assert.AreEqual(
                12,
                GetArwmhCovarianceSampleCount(sampler),
                "Continual adaptation must keep every realized transition before and after the boundary.");
        }

        /// <summary>
        /// Confirms that NUTS diagnostic arrays are non-null and empty before sampling initializes them.
        /// </summary>
        [TestMethod]
        public void NUTS_DiagnosticArraysAreNonNullAndEmptyBeforeSampling()
        {
            SeedableNUTS sampler = CreateNutsSampler();

            Assert.IsEmpty(sampler.HamiltonianAcceptanceRates);
            Assert.IsEmpty(sampler.DiagnosticSampleCounts);
            Assert.IsEmpty(sampler.DivergenceCounts);
            Assert.IsEmpty(sampler.MaxTreeDepthHitCounts);
            Assert.IsEmpty(sampler.MeanTreeDepths);
            Assert.IsEmpty(sampler.MeanLeapfrogSteps);
            Assert.IsEmpty(sampler.StepSizes);
            Assert.IsEmpty(sampler.EnergyBayesianFractionOfMissingInformation);
        }

        /// <summary>
        /// Confirms that generic acceptance retains accepted-transition semantics while
        /// NUTS results persist the Hamiltonian acceptance statistic used by BestFit.
        /// </summary>
        [TestMethod]
        public void NUTS_AcceptanceContractsRemainSeparatedAndResultsPersistHamiltonianRates()
        {
            SeedableNUTS sampler = CreateNutsSampler();
            sampler.Seed(new ParameterSet(new[] { 0d, 0d }, 0d));

            sampler.Sample();

            for (int chainIndex = 0; chainIndex < sampler.NumberOfChains; chainIndex++)
            {
                double expectedTransitionRate = (double)sampler.AcceptCount[chainIndex] /
                    sampler.SampleCount[chainIndex];
                Assert.AreEqual(sampler.SampleCount[chainIndex], sampler.AcceptCount[chainIndex]);
                Assert.AreEqual(expectedTransitionRate, sampler.AcceptanceRates[chainIndex], 0d);
                Assert.AreEqual(1d, sampler.AcceptanceRates[chainIndex], 0d);
                Assert.AreEqual(
                    sampler.SampleCount[chainIndex] - sampler.WarmupIterations,
                    sampler.DiagnosticSampleCounts[chainIndex]);
                Assert.IsGreaterThan(0d, sampler.HamiltonianAcceptanceRates[chainIndex]);
                Assert.IsLessThan(1d, sampler.HamiltonianAcceptanceRates[chainIndex]);
                Assert.AreNotEqual(
                    sampler.AcceptanceRates[chainIndex],
                    sampler.HamiltonianAcceptanceRates[chainIndex]);
                Assert.IsGreaterThanOrEqualTo(1d, sampler.MeanTreeDepths[chainIndex]);
                Assert.IsLessThanOrEqualTo(4d, sampler.MeanTreeDepths[chainIndex]);
                Assert.IsGreaterThanOrEqualTo(1d, sampler.MeanLeapfrogSteps[chainIndex]);
                Assert.IsGreaterThan(0d, sampler.StepSizes[chainIndex]);
                Assert.IsTrue(Tools.IsFinite(sampler.EnergyBayesianFractionOfMissingInformation[chainIndex]));
            }

            var results = new MCMCResults(sampler);
            CollectionAssert.AreEqual(sampler.HamiltonianAcceptanceRates, results.AcceptanceRates);

            byte[] serialized = MCMCResults.ToByteArray(results);
            string currentJson = System.Text.Encoding.UTF8.GetString(serialized);
            Assert.IsFalse(currentJson.Contains("NUTSDivergenceCounts", StringComparison.Ordinal));
            string staleJson = currentJson.Insert(1,
                "\"NUTSDiagnosticSampleCounts\":[100,100]," +
                "\"NUTSDivergenceCounts\":[0,0]," +
                "\"NUTSMaxTreeDepthHitCounts\":[0,0]," +
                "\"NUTSMeanTreeDepths\":[2.0,2.0]," +
                "\"NUTSMeanLeapfrogSteps\":[4.0,4.0]," +
                "\"NUTSStepSizes\":[0.1,0.1]," +
                "\"NUTSEnergyBayesianFractionOfMissingInformation\":[0.8,0.8],");
            MCMCResults restored = MCMCResults.FromByteArray(System.Text.Encoding.UTF8.GetBytes(staleJson));
            Assert.IsNotNull(restored);
            CollectionAssert.AreEqual(results.AcceptanceRates, restored.AcceptanceRates);
        }

        /// <summary>
        /// Confirms that MCMCResults retains its baseline nullability metadata and has no
        /// sampler-specific NUTS result properties.
        /// </summary>
        [TestMethod]
        public void MCMCResults_RetainsBaselineNullabilityAndOmitsNutsDiagnostics()
        {
#if NET6_0_OR_GREATER
            var nullabilityContext = new NullabilityInfoContext();
            (string Name, NullabilityState State)[] expectedStates =
            {
                (nameof(MCMCResults.MarkovChains), NullabilityState.Nullable),
                (nameof(MCMCResults.MeanLogLikelihood), NullabilityState.Nullable),
                (nameof(MCMCResults.Output), NullabilityState.NotNull),
                (nameof(MCMCResults.AcceptanceRates), NullabilityState.NotNull),
                (nameof(MCMCResults.ParameterResults), NullabilityState.NotNull),
                (nameof(MCMCResults.MAP), NullabilityState.NotNull),
                (nameof(MCMCResults.PosteriorMean), NullabilityState.NotNull)
            };

            foreach ((string propertyName, NullabilityState expectedState) in expectedStates)
            {
                PropertyInfo property = typeof(MCMCResults).GetProperty(propertyName);
                Assert.IsNotNull(property, $"Missing baseline MCMCResults property {propertyName}.");
                Assert.AreEqual(
                    expectedState,
                    nullabilityContext.Create(property).ReadState,
                    $"Unexpected nullability metadata for MCMCResults.{propertyName}.");
            }

#endif

            string[] removedProperties =
            {
                "NUTSDiagnosticSampleCounts",
                "NUTSDivergenceCounts",
                "NUTSMaxTreeDepthHitCounts",
                "NUTSMeanTreeDepths",
                "NUTSMeanLeapfrogSteps",
                "NUTSStepSizes",
                "NUTSEnergyBayesianFractionOfMissingInformation"
            };
            foreach (string propertyName in removedProperties)
            {
                Assert.IsNull(
                    typeof(MCMCResults).GetProperty(propertyName),
                    $"NUTS diagnostic {propertyName} must remain on NUTS rather than MCMCResults.");
            }
        }
        /// <summary>
        /// Confirms the streaming E-BFMI calculation matches Stan's
        /// <c>mean(diff(E)^2) / var(E)</c> convention.
        /// </summary>
        [TestMethod]
        public void NUTS_EnergyBayesianFractionOfMissingInformationMatchesStanFormula()
        {
            double actual = NUTS.ComputeEnergyBayesianFractionOfMissingInformation(
                sampleCount: 4,
                energyM2: 5d,
                squaredDifferenceSum: 6d);

            Assert.AreEqual(0.9d, actual, 1e-15);
            Assert.IsTrue(double.IsNaN(
                NUTS.ComputeEnergyBayesianFractionOfMissingInformation(1, 0d, 0d)));
        }

        /// <summary>
        /// Creates the deterministic two-chain NUTS fixture used by acceptance-contract tests.
        /// </summary>
        /// <returns>The configured sampler before sampling or user-defined seeding.</returns>
        private static SeedableNUTS CreateNutsSampler()
        {
            var priors = new List<IUnivariateDistribution>
            {
                new Uniform(-20d, 20d),
                new Uniform(-20d, 20d)
            };
            double LogTarget(double[] values) =>
                -0.5d * (values[0] * values[0] + values[1] * values[1]);
            Vector Gradient(IList<double> values) => new Vector(new[] { -values[0], -values[1] });
            return new SeedableNUTS(priors, LogTarget, Gradient)
            {
                NumberOfChains = 2,
                InitialIterations = 2,
                WarmupIterations = 50,
                Iterations = 100,
                OutputLength = 100,
                ThinningInterval = 1,
                PRNGSeed = 8675309,
                ParallelizeChains = false
            };
        }
        /// <summary>
        /// Creates a deterministic ARWMH sampler whose only finite-density state is the origin.
        /// </summary>
        /// <param name="warmupIterations">Configured warmup boundary in raw transitions.</param>
        /// <returns>The configured inspectable sampler.</returns>
        private static InspectableARWMH CreateRejectingArwmh(int warmupIterations)
        {
            var priors = new List<IUnivariateDistribution>
            {
                new Uniform(-10d, 10d),
                new Uniform(-10d, 10d)
            };
            double SpikeTarget(double[] values) =>
                values[0] == 0d && values[1] == 0d ? 0d : double.NegativeInfinity;
            return new InspectableARWMH(priors, SpikeTarget)
            {
                NumberOfChains = 1,
                WarmupIterations = warmupIterations,
                ThinningInterval = 1,
                PRNGSeed = 314159,
                Beta = 1d
            };
        }

        /// <summary>
        /// Reads the current running-covariance observation count for the first ARWMH chain.
        /// </summary>
        /// <param name="sampler">Sampler to inspect.</param>
        /// <returns>The number of states pushed into the running covariance.</returns>
        private static int GetArwmhCovarianceSampleCount(ARWMH sampler)
        {
            FieldInfo field = typeof(ARWMH).GetField("sigma", BindingFlags.Instance | BindingFlags.NonPublic);
            Assert.IsNotNull(field, "The ARWMH running covariance field must be accessible.");
            var covariance = field.GetValue(sampler) as RunningCovarianceMatrix[];
            Assert.IsNotNull(covariance);
            return covariance[0].N;
        }


        /// <summary>
        /// Exposes deterministic initialization and one ARWMH transition for focused tests.
        /// </summary>
        private sealed class InspectableARWMH : ARWMH
        {
            /// <summary>
            /// Initializes a test sampler.
            /// </summary>
            /// <param name="priors">Parameter prior bounds.</param>
            /// <param name="target">Log target.</param>
            public InspectableARWMH(List<IUnivariateDistribution> priors, LogLikelihood target)
                : base(priors, target)
            {
            }

            /// <summary>
            /// Initializes sampler-specific state at a deterministic chain state.
            /// </summary>
            /// <param name="state">Initial state.</param>
            internal void InitializeAt(ParameterSet state)
            {
                Reset();
                _chainStates = new[] { state.Clone() };
                InitializeCustomSettings();
            }

            /// <summary>
            /// Executes one raw transition.
            /// </summary>
            /// <param name="state">Current state.</param>
            /// <returns>The retained state.</returns>
            internal ParameterSet Iterate(ParameterSet state) => ChainIteration(0, state);
        }

        /// <summary>
        /// Allows a deterministic user-defined starting state for a focused NUTS run.
        /// </summary>
        private sealed class SeedableNUTS : NUTS
        {
            /// <summary>
            /// Initializes a test sampler with an analytic gradient.
            /// </summary>
            /// <param name="priors">Parameter prior bounds.</param>
            /// <param name="target">Log target.</param>
            /// <param name="gradient">Analytic log-target gradient.</param>
            public SeedableNUTS(List<IUnivariateDistribution> priors, LogLikelihood target, HMC.Gradient gradient)
                : base(priors, target, maxTreeDepth: 4, gradientFunction: gradient)
            {
            }

            /// <summary>
            /// Seeds every chain for user-defined initialization.
            /// </summary>
            /// <param name="state">Initial state.</param>
            internal void Seed(ParameterSet state)
            {
                Reset();
                for (int chainIndex = 0; chainIndex < NumberOfChains; chainIndex++)
                    MarkovChains[chainIndex].Add(state.Clone());
                Initialize = InitializationType.UserDefined;
            }
        }
    }
}
