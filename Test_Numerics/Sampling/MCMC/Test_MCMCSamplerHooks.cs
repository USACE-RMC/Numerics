using System;
using System.Collections.Generic;
using System.Threading;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit tests for the chain index and external cancellation token on <see cref="MCMCSampler"/>.
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
    public class Test_MCMCSamplerHooks
    {
        /// <summary>
        /// Builds a seeded two-parameter RWMH sampler on a standard normal target.
        /// </summary>
        /// <param name="logLikelihood">The log-likelihood delegate.</param>
        /// <returns>An unsampled sampler.</returns>
        private static RWMH BuildSampler(LogLikelihood logLikelihood)
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(-10d, 10d), new Uniform(-10d, 10d) };
            return new RWMH(priors, logLikelihood, Matrix.Identity(2) * 0.5)
            {
                PRNGSeed = 12345,
                NumberOfChains = 4,
                Iterations = 200,
                WarmupIterations = 50,
                ThinningInterval = 2,
                OutputLength = 100
            };
        }

        /// <summary>
        /// The standard normal log-density, up to an additive constant.
        /// </summary>
        /// <param name="x">The parameter vector.</param>
        /// <returns>The log-density.</returns>
        private static double StandardNormal(double[] x) => -0.5 * (x[0] * x[0] + x[1] * x[1]);

        /// <summary>
        /// Runs a sampler whose likelihood holds a per-chain resource that fails if two threads use it at
        /// once, and checks that the chain index lets every chain own its resource.
        /// </summary>
        /// <param name="parallel">Whether chains run in parallel.</param>
        private static void AssertChainIndexSelectsExclusiveResource(bool parallel)
        {
            // Slot 0 is the resource for index -1 (initialization); slots 1-4 belong to chains 0-3.
            var busy = new int[5];
            var calls = new int[5];
            int collisions = 0;
            int outOfRange = 0;

            double LogLikelihood(double[] x)
            {
                int slot = MCMCSampler.CurrentChainIndex + 1;
                if (slot < 0 || slot >= busy.Length)
                {
                    Interlocked.Increment(ref outOfRange);
                    return StandardNormal(x);
                }
                if (Interlocked.Exchange(ref busy[slot], 1) != 0)
                    Interlocked.Increment(ref collisions);
                Interlocked.Increment(ref calls[slot]);
                Thread.SpinWait(50);
                double value = StandardNormal(x);
                Volatile.Write(ref busy[slot], 0);
                return value;
            }

            var sampler = BuildSampler(LogLikelihood);
            sampler.ParallelizeChains = parallel;
            sampler.Sample();

            Assert.AreEqual(0, outOfRange, "CurrentChainIndex left the range -1 to NumberOfChains - 1.");
            Assert.AreEqual(0, collisions, "Two threads used the same per-chain resource at once.");
            Assert.IsGreaterThan(0, calls[0], "Initialization must run with CurrentChainIndex = -1.");
            for (int c = 1; c < calls.Length; c++)
                Assert.IsGreaterThan(0, calls[c], $"Chain {c - 1} never evaluated the likelihood under its own index.");
            Assert.AreEqual(-1, MCMCSampler.CurrentChainIndex, "The index must be restored after sampling.");
        }

        /// <summary>
        /// With parallel chains, each chain's likelihood calls see that chain's index, so a per-chain
        /// resource is never used by two threads at once.
        /// </summary>
        [TestMethod]
        public void Test_CurrentChainIndex_Parallel_SelectsPerChainResource()
        {
            AssertChainIndexSelectsExclusiveResource(true);
        }

        /// <summary>
        /// With serial chains, the chain index is set in the same way.
        /// </summary>
        [TestMethod]
        public void Test_CurrentChainIndex_Serial_SelectsPerChainResource()
        {
            AssertChainIndexSelectsExclusiveResource(false);
        }

        /// <summary>
        /// Outside a sampler the chain index is -1.
        /// </summary>
        [TestMethod]
        public void Test_CurrentChainIndex_OutsideSampling_IsMinusOne()
        {
            Assert.AreEqual(-1, MCMCSampler.CurrentChainIndex);
        }

        /// <summary>
        /// The NUTS gradient delegate sees the chain index too.
        /// </summary>
        [TestMethod]
        public void Test_CurrentChainIndex_VisibleToNutsGradient()
        {
            var seen = new int[3];
            Vector Gradient(IList<double> x)
            {
                int index = MCMCSampler.CurrentChainIndex;
                if (index >= -1 && index <= 1) Interlocked.Increment(ref seen[index + 1]);
                return new Vector(new[] { -x[0] });
            }

            var priors = new List<IUnivariateDistribution> { new Uniform(-10d, 10d) };
            var sampler = new NUTS(priors, x => -0.5 * x[0] * x[0], gradientFunction: Gradient)
            {
                PRNGSeed = 12345,
                NumberOfChains = 2,
                Iterations = 100,
                WarmupIterations = 50,
                ThinningInterval = 1,
                OutputLength = 100
            };
            sampler.Sample();

            Assert.IsGreaterThan(0, seen[1], "Chain 0 gradient calls must see index 0.");
            Assert.IsGreaterThan(0, seen[2], "Chain 1 gradient calls must see index 1.");
        }

        /// <summary>
        /// Reading the chain index does not change the sampled values.
        /// </summary>
        [TestMethod]
        public void Test_CurrentChainIndex_DoesNotChangeDraws()
        {
            var plain = BuildSampler(StandardNormal);
            plain.Sample();

            int sink = 0;
            var reading = BuildSampler(x => { Interlocked.Add(ref sink, MCMCSampler.CurrentChainIndex); return StandardNormal(x); });
            reading.Sample();

            for (int c = 0; c < plain.NumberOfChains; c++)
            {
                Assert.HasCount(plain.Output[c].Count, reading.Output[c]);
                for (int i = 0; i < plain.Output[c].Count; i++)
                    CollectionAssert.AreEqual(plain.Output[c][i].Values, reading.Output[c][i].Values);
            }
        }

        /// <summary>
        /// A token cancelled before sampling stops the run after one recorded iteration, without throwing.
        /// </summary>
        [TestMethod]
        public void Test_CancellationToken_AlreadyCancelled_StopsAfterOneIteration()
        {
            using var source = new CancellationTokenSource();
            source.Cancel();
            var sampler = BuildSampler(StandardNormal);
            sampler.CancellationToken = source.Token;
            sampler.Sample();

            for (int c = 0; c < sampler.NumberOfChains; c++)
                Assert.HasCount(1, sampler.MarkovChains[c]);
        }

        /// <summary>
        /// A token cancelled during sampling stops the run early.
        /// </summary>
        [TestMethod]
        public void Test_CancellationToken_CancelledDuringRun_StopsEarly()
        {
            using var source = new CancellationTokenSource();
            int calls = 0;
            var sampler = BuildSampler(x =>
            {
                if (Interlocked.Increment(ref calls) == 400) source.Cancel();
                return StandardNormal(x);
            });
            sampler.CancellationToken = source.Token;
            sampler.Sample();

            Assert.IsLessThan(sampler.Iterations, sampler.MarkovChains[0].Count);
            Assert.IsEmpty(sampler.Output[0]);
        }

        /// <summary>
        /// A token that is never cancelled runs the full sample and changes no draw.
        /// </summary>
        [TestMethod]
        public void Test_CancellationToken_NotCancelled_DoesNotChangeDraws()
        {
            var plain = BuildSampler(StandardNormal);
            plain.Sample();

            using var source = new CancellationTokenSource();
            var withToken = BuildSampler(StandardNormal);
            withToken.CancellationToken = source.Token;
            withToken.Sample();

            for (int c = 0; c < plain.NumberOfChains; c++)
            {
                Assert.HasCount(plain.Iterations, withToken.MarkovChains[c]);
                Assert.HasCount(plain.Output[c].Count, withToken.Output[c]);
                for (int i = 0; i < plain.Output[c].Count; i++)
                    CollectionAssert.AreEqual(plain.Output[c][i].Values, withToken.Output[c][i].Values);
            }
        }
    }
}
