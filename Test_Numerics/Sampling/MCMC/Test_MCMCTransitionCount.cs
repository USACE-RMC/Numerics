using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit tests for the derived work-estimate properties on <see cref="MCMCSampler"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <see cref="MCMCSampler.Iterations"/> does not describe the work a run performs. Two multipliers sit
    /// between the two figures: <c>Sample()</c> runs ceil(<c>OutputLength</c> / <c>NumberOfChains</c>)
    /// recorded iterations beyond <c>Iterations</c> to collect the posterior output, and each recorded
    /// iteration advances the chain <c>ThinningInterval</c> times. <see cref="MCMCSampler.TransitionCount"/>
    /// and <see cref="MCMCSampler.TotalTransitionCount"/> report the product.
    /// </para>
    /// <para>
    /// The expected values below are computed by hand rather than from the properties' own expressions, and
    /// <see cref="TransitionCount_MatchesTheTransitionsSampleActuallyPerforms"/> ties the arithmetic to the
    /// real <c>Sample()</c> loop so that the two cannot drift apart.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_MCMCTransitionCount
    {
        /// <summary>
        /// Builds a minimal, cheap RWMH sampler. The transition counts are pure functions of the settings,
        /// so the priors and likelihood only have to be well formed.
        /// </summary>
        /// <returns>A two-parameter RWMH sampler left at its default settings.</returns>
        private static RWMH CreateSampler()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(-10d, 10d), new Uniform(0.1d, 10d) };
            double logLH(double[] x) => new Normal(x[0], x[1]).LogPDF(0d);
            return new RWMH(priors, logLH, Matrix.Identity(2));
        }

        /// <summary>
        /// At the shipped defaults a single chain performs 120,000 transitions and the four chains together
        /// perform 480,000, against an <see cref="MCMCSampler.Iterations"/> that reads 3,500.
        /// </summary>
        [TestMethod]
        public void TransitionCount_AtTheDefaults()
        {
            var sampler = CreateSampler();

            // Guard the defaults the hand-computed values below depend on.
            Assert.AreEqual(3500, sampler.Iterations);
            Assert.AreEqual(10000, sampler.OutputLength);
            Assert.AreEqual(4, sampler.NumberOfChains);
            Assert.AreEqual(20, sampler.ThinningInterval);
            Assert.AreEqual(1750, sampler.WarmupIterations);

            // (3500 + ceil(10000 / 4)) * 20 = (3500 + 2500) * 20 = 120,000.
            Assert.AreEqual(120000L, sampler.TransitionCount);
            Assert.AreEqual(480000L, sampler.TotalTransitionCount);

            // Warmup is a subset of Iterations, not an addition to it, so it is already inside the figure.
            Assert.IsLessThanOrEqualTo((int)(0.5 * sampler.Iterations), sampler.WarmupIterations);
        }

        /// <summary>
        /// A configuration in which <c>OutputLength / NumberOfChains</c> does not divide evenly, so the
        /// ceiling in the output-iteration count is exercised.
        /// </summary>
        [TestMethod]
        public void TransitionCount_WhenOutputDoesNotDivideEvenly()
        {
            var sampler = CreateSampler();
            sampler.Iterations = 1000;
            sampler.OutputLength = 10001;
            sampler.NumberOfChains = 3;
            sampler.ThinningInterval = 7;

            // 10001 / 3 = 3333.667, so ceil gives 3334 — one more than truncation would.
            // (1000 + 3334) * 7 = 4334 * 7 = 30,338.
            Assert.AreEqual(30338L, sampler.TransitionCount);
            Assert.AreEqual(91014L, sampler.TotalTransitionCount);
        }

        /// <summary>
        /// A second non-default configuration with a different uneven division and a different thinning
        /// interval.
        /// </summary>
        [TestMethod]
        public void TransitionCount_AtASecondNonDefaultConfiguration()
        {
            var sampler = CreateSampler();
            sampler.Iterations = 250;
            sampler.OutputLength = 1000;
            sampler.NumberOfChains = 7;
            sampler.ThinningInterval = 3;

            // 1000 / 7 = 142.857, so ceil gives 143.
            // (250 + 143) * 3 = 393 * 3 = 1179.
            Assert.AreEqual(1179L, sampler.TransitionCount);
            Assert.AreEqual(8253L, sampler.TotalTransitionCount);

            // A thinning interval of 1 removes that multiplier entirely, leaving the recorded iterations.
            sampler.ThinningInterval = 1;
            Assert.AreEqual(393L, sampler.TransitionCount);
            Assert.AreEqual(2751L, sampler.TotalTransitionCount);
        }

        /// <summary>
        /// The counts are <see cref="long"/> so that settings the sampler accepts cannot overflow them.
        /// </summary>
        /// <remarks>
        /// Nothing bounds <see cref="MCMCSampler.Iterations"/> from above, so the product leaves
        /// <see cref="int"/> range long before the settings become unreasonable. In <see cref="int"/>
        /// arithmetic the per-chain figure below would have wrapped to a negative number.
        /// </remarks>
        [TestMethod]
        public void TransitionCount_DoesNotOverflowForLargeSettings()
        {
            var sampler = CreateSampler();
            sampler.Iterations = 200000000;
            sampler.OutputLength = 10000;
            sampler.NumberOfChains = 4;
            sampler.ThinningInterval = 20;

            // (200,000,000 + 2,500) * 20 = 4,000,050,000, which exceeds int.MaxValue (2,147,483,647).
            Assert.AreEqual(4000050000L, sampler.TransitionCount);
            Assert.AreEqual(16000200000L, sampler.TotalTransitionCount);
            Assert.IsGreaterThan(int.MaxValue, sampler.TransitionCount);
        }

        /// <summary>
        /// Zero chains is not a samplable configuration, so both counts report zero rather than
        /// dividing by zero. The properties are readable at any time; the configuration itself is
        /// only rejected by <c>ValidateSettings</c> when <c>Sample()</c> runs.
        /// </summary>
        [TestMethod]
        public void TransitionCount_AtZeroChains_ReportsZero()
        {
            var sampler = CreateSampler();
            sampler.NumberOfChains = 0;

            Assert.AreEqual(0L, sampler.TransitionCount);
            Assert.AreEqual(0L, sampler.TotalTransitionCount);
        }

        /// <summary>
        /// Ties <see cref="MCMCSampler.TransitionCount"/> to the transitions <c>Sample()</c> actually
        /// performs, so the property cannot drift away from the loop it describes.
        /// </summary>
        /// <remarks>
        /// <c>SampleChain</c> is invoked once per recorded iteration per chain and advances the chain
        /// <c>ThinningInterval</c> times, so the observed call count times the thinning interval must equal
        /// the advertised per-chain transition count. Chains are sampled serially here so the counters are
        /// deterministic. The configuration is the smallest one <c>ValidateSettings</c> accepts, to keep the
        /// run cheap.
        /// </remarks>
        [TestMethod]
        public void TransitionCount_MatchesTheTransitionsSampleActuallyPerforms()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(-5d, 5d), new Uniform(0.5d, 5d) };
            double logLH(double[] x) => new Normal(x[0], x[1]).LogPDF(0d);
            var sampler = new CountingRWMH(priors, logLH, Matrix.Identity(2))
            {
                Iterations = 100,
                WarmupIterations = 10,
                OutputLength = 100,
                NumberOfChains = 2,
                ThinningInterval = 3,
                InitialIterations = 10,
                ParallelizeChains = false,
                PRNGSeed = 12345
            };

            // ceil(100 / 2) = 50 output iterations, so 150 recorded iterations, each of 3 transitions.
            Assert.AreEqual(450L, sampler.TransitionCount);
            Assert.AreEqual(900L, sampler.TotalTransitionCount);

            sampler.Sample();

            // Compare against the COUNTED transitions directly. Nothing here multiplies by
            // ThinningInterval: TransitionCallCount is incremented once per ChainIteration, which is one
            // real transition, so if the inner loop in SampleChain ever ran a different number of times
            // than ThinningInterval the counts would diverge and this would fail. Multiplying an observed
            // SampleChain count by ThinningInterval would instead assume the very thing under test.
            Assert.AreEqual(sampler.TransitionCount, sampler.TransitionCallCount[0],
                "TransitionCount must equal the transitions Sample() actually performs on one chain.");

            long observedTotal = 0;
            for (int i = 0; i < sampler.NumberOfChains; i++)
            {
                Assert.AreEqual(sampler.TransitionCallCount[0], sampler.TransitionCallCount[i],
                    "Every chain is advanced the same number of times.");
                observedTotal += sampler.TransitionCallCount[i];
            }
            Assert.AreEqual(sampler.TotalTransitionCount, observedTotal,
                "TotalTransitionCount must equal the transitions Sample() performs across all chains.");

            // The recorded-iteration count is a separate fact worth pinning: it is the outer loop,
            // Iterations + ceil(OutputLength / NumberOfChains) = 150, and the transitions are that times
            // the thinning interval. Asserting both separately locates a future break to one loop or the
            // other rather than just reporting a mismatched product.
            Assert.AreEqual(150L, sampler.ChainCallCount[0],
                "SampleChain must be invoked once per recorded iteration.");
            Assert.AreEqual(sampler.ChainCallCount[0] * sampler.ThinningInterval, sampler.TransitionCallCount[0],
                "Each recorded iteration must advance the chain exactly ThinningInterval times.");
        }

        /// <summary>
        /// An RWMH sampler that records both how many times each chain is advanced (<c>SampleChain</c>, the
        /// outer recorded-iteration loop) and how many transitions actually occur (<c>ChainIteration</c>,
        /// the inner thinning loop).
        /// </summary>
        /// <remarks>
        /// Counting <c>ChainIteration</c> is what makes the anti-drift assertion real. Counting only
        /// <c>SampleChain</c> and multiplying by <see cref="MCMCSampler.ThinningInterval"/> would assume the
        /// inner loop performs exactly that many transitions, which is one of the two things
        /// <see cref="MCMCSampler.TransitionCount"/> claims.
        /// </remarks>
        private sealed class CountingRWMH : RWMH
        {
            /// <summary>
            /// Initializes a counting sampler.
            /// </summary>
            /// <param name="priors">Parameter priors.</param>
            /// <param name="target">Log-likelihood function.</param>
            /// <param name="proposalSigma">Proposal covariance.</param>
            public CountingRWMH(List<IUnivariateDistribution> priors, LogLikelihood target, Matrix proposalSigma)
                : base(priors, target, proposalSigma)
            {
            }

            /// <summary>
            /// The number of times each chain's recorded iteration ran, indexed by chain.
            /// </summary>
            public long[] ChainCallCount { get; private set; } = [];

            /// <summary>
            /// The number of transitions each chain actually performed, indexed by chain.
            /// </summary>
            public long[] TransitionCallCount { get; private set; } = [];

            /// <inheritdoc/>
            public override void Sample()
            {
                // Sized here rather than at construction because NumberOfChains is assigned by the caller
                // after the constructor runs.
                ChainCallCount = new long[NumberOfChains];
                TransitionCallCount = new long[NumberOfChains];
                base.Sample();
            }

            /// <inheritdoc/>
            protected override ParameterSet SampleChain(int index, ParameterSet state)
            {
                ChainCallCount[index]++;
                return base.SampleChain(index, state);
            }

            /// <inheritdoc/>
            protected override ParameterSet ChainIteration(int index, ParameterSet state)
            {
                TransitionCallCount[index]++;
                return base.ChainIteration(index, state);
            }
        }
    }
}
