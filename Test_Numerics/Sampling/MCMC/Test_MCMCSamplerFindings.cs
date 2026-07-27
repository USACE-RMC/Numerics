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
    /// Characterizes the adaptive-random-walk and NUTS diagnostic findings used by
    /// the RMC.BestFit verification program.
    /// </summary>
    /// <remarks>
    /// These tests intentionally pin the current public behavior before any sampler
    /// correction. They use deterministic inline targets and do not require R or Python.
    /// </remarks>
    [TestClass]
    public class Test_MCMCSamplerFindings
    {
        /// <summary>
        /// Confirms that rejected transitions before the warmup boundary advance the
        /// running covariance sample count in the corrected ARWMH implementation.
        /// </summary>
        /// <remarks>
        /// Adaptive Metropolis covariance is defined from the realized chain history;
        /// a rejected proposal therefore contributes the repeated retained state.
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
        /// Haario-Saksman-Tamminen construction. This test separates that intended
        /// behavior from the omitted repeated states before the boundary.
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
        /// Confirms that NUTS exposes post-warmup Hamiltonian acceptance and additive sampler diagnostics.
        /// </summary>
        [TestMethod]
        public void NUTS_AcceptanceRatesAndDiagnosticsUsePostWarmupHamiltonianStatistics()
        {
            var priors = new List<IUnivariateDistribution>
            {
                new Uniform(-20d, 20d),
                new Uniform(-20d, 20d)
            };
            double LogTarget(double[] values) => -0.5d * (values[0] * values[0] + values[1] * values[1]);
            Vector Gradient(IList<double> values) => new Vector(new[] { -values[0], -values[1] });
            var sampler = new SeedableNUTS(priors, LogTarget, Gradient)
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
            sampler.Seed(new ParameterSet(new[] { 0d, 0d }, 0d));

            sampler.Sample();

            for (int chainIndex = 0; chainIndex < sampler.NumberOfChains; chainIndex++)
            {
                Assert.AreEqual(sampler.SampleCount[chainIndex], sampler.AcceptCount[chainIndex]);
                Assert.AreEqual(
                    sampler.SampleCount[chainIndex] - sampler.WarmupIterations,
                    sampler.DiagnosticSampleCounts[chainIndex]);
                Assert.AreEqual(
                    sampler.HamiltonianAcceptanceRates[chainIndex],
                    sampler.AcceptanceRates[chainIndex],
                    0d);
                Assert.IsGreaterThan(0d, sampler.AcceptanceRates[chainIndex]);
                Assert.IsLessThan(1d, sampler.AcceptanceRates[chainIndex]);
                Assert.IsGreaterThanOrEqualTo(1d, sampler.MeanTreeDepths[chainIndex]);
                Assert.IsLessThanOrEqualTo(4d, sampler.MeanTreeDepths[chainIndex]);
                Assert.IsGreaterThanOrEqualTo(1d, sampler.MeanLeapfrogSteps[chainIndex]);
                Assert.IsGreaterThan(0d, sampler.StepSizes[chainIndex]);
                Assert.IsTrue(Tools.IsFinite(sampler.EnergyBayesianFractionOfMissingInformation[chainIndex]));
            }

            var results = new MCMCResults(sampler);
            CollectionAssert.AreEqual(sampler.AcceptanceRates, results.AcceptanceRates);
            CollectionAssert.AreEqual(sampler.DiagnosticSampleCounts, results.NUTSDiagnosticSampleCounts);
            CollectionAssert.AreEqual(sampler.DivergenceCounts, results.NUTSDivergenceCounts);
            CollectionAssert.AreEqual(sampler.MaxTreeDepthHitCounts, results.NUTSMaxTreeDepthHitCounts);
            CollectionAssert.AreEqual(sampler.MeanTreeDepths, results.NUTSMeanTreeDepths);
            CollectionAssert.AreEqual(sampler.MeanLeapfrogSteps, results.NUTSMeanLeapfrogSteps);
            CollectionAssert.AreEqual(sampler.StepSizes, results.NUTSStepSizes);
            CollectionAssert.AreEqual(
                sampler.EnergyBayesianFractionOfMissingInformation,
                results.NUTSEnergyBayesianFractionOfMissingInformation);

            byte[] serialized = MCMCResults.ToByteArray(results);
            MCMCResults restored = MCMCResults.FromByteArray(serialized);
            Assert.IsNotNull(restored);
            CollectionAssert.AreEqual(results.AcceptanceRates, restored.AcceptanceRates);
            CollectionAssert.AreEqual(results.NUTSDiagnosticSampleCounts, restored.NUTSDiagnosticSampleCounts);
            CollectionAssert.AreEqual(results.NUTSDivergenceCounts, restored.NUTSDivergenceCounts);
            CollectionAssert.AreEqual(results.NUTSMaxTreeDepthHitCounts, restored.NUTSMaxTreeDepthHitCounts);
            CollectionAssert.AreEqual(results.NUTSMeanTreeDepths, restored.NUTSMeanTreeDepths);
            CollectionAssert.AreEqual(results.NUTSMeanLeapfrogSteps, restored.NUTSMeanLeapfrogSteps);
            CollectionAssert.AreEqual(results.NUTSStepSizes, restored.NUTSStepSizes);
            CollectionAssert.AreEqual(
                results.NUTSEnergyBayesianFractionOfMissingInformation,
                restored.NUTSEnergyBayesianFractionOfMissingInformation);
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
            Assert.IsNotNull(field, "The ARWMH running covariance field must be available for characterization.");
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
