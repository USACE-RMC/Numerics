using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling.MCMC;
using Utilities;

namespace Sampling.MCMC
{
    /// <summary>
    /// Reference-stream tests that pin the exact seeded output of the NUTS and ARWMH samplers.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// Each test hashes the bit pattern of every recorded chain state, every output draw and its fitness
    /// with 64-bit FNV-1a and compares it with a reference value for the same seeded configuration. The
    /// NUTS hashes use the step size fixed at the smoothed dual-averaging value from the first sampling
    /// transition; the ARWMH hash does not depend on that rule.
    /// </para>
    /// <para>
    /// The reference hashes were recorded on Windows and are checked there. The samplers call Math.Exp and Math.Log
    /// (acceptance tests, proposals and the normal deviates), which the C library on Linux rounds differently for some
    /// arguments. On Linux the NUTS chains first differ by a few ulp at the 32nd (analytic gradient) or 141st (numerical
    /// gradient) recorded state, and the trajectories then diverge, so on other platforms the NUTS tests check instead
    /// that the output reproduces the mean and second moment of every coordinate of the target within Monte Carlo error.
    /// The ARWMH chains differ in six recorded states of one chain, by at most 64 ulp in one coordinate, and then rejoin;
    /// its output is the same bit for bit, so on other platforms the ARWMH test pins a hash of the output alone. The
    /// moment check runs on Windows as well (see <see cref="CLibraryRounding"/>).
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_MCMCReferenceStreams
    {
        /// <summary>The per-coordinate standard deviations of the target.</summary>
        private static readonly double[] Sd = { 0.5, 1d, 2d, 4d, 8d };

        /// <summary>
        /// The log-density of a zero-mean diagonal Gaussian, up to an additive constant.
        /// </summary>
        /// <param name="x">The parameter vector.</param>
        /// <returns>The log-density.</returns>
        private static double LogLikelihood(double[] x)
        {
            double sum = 0;
            for (int i = 0; i < x.Length; i++)
                sum += -0.5 * x[i] * x[i] / (Sd[i] * Sd[i]);
            return sum;
        }

        /// <summary>
        /// The exact gradient of <see cref="LogLikelihood(double[])"/>.
        /// </summary>
        /// <param name="x">The parameter vector.</param>
        /// <returns>The gradient.</returns>
        private static Vector Gradient(IList<double> x)
        {
            var g = new Vector(x.Count);
            for (int i = 0; i < x.Count; i++)
                g[i] = -x[i] / (Sd[i] * Sd[i]);
            return g;
        }

        /// <summary>
        /// Builds the seeded configuration the reference hashes were captured from.
        /// </summary>
        /// <param name="sampler">The unsampled sampler.</param>
        /// <param name="thinning">The thinning interval.</param>
        /// <returns>The configured sampler.</returns>
        private static MCMCSampler Configure(MCMCSampler sampler, int thinning)
        {
            sampler.PRNGSeed = 2468;
            sampler.NumberOfChains = 4;
            sampler.Iterations = 600;
            sampler.WarmupIterations = 150;
            sampler.ThinningInterval = thinning;
            sampler.OutputLength = 400;
            return sampler;
        }

        /// <summary>
        /// Builds bounded uniform priors for the target.
        /// </summary>
        /// <returns>The prior list.</returns>
        private static List<IUnivariateDistribution> Priors()
        {
            var priors = new List<IUnivariateDistribution>();
            for (int i = 0; i < Sd.Length; i++) priors.Add(new Uniform(-50d, 50d));
            return priors;
        }

        /// <summary>
        /// Hashes every recorded chain state and output draw with 64-bit FNV-1a over the IEEE bit patterns.
        /// </summary>
        /// <param name="sampler">A sampled sampler.</param>
        /// <param name="includeChains">Whether the recorded chain states are hashed before the output draws.</param>
        /// <returns>The hash.</returns>
        private static ulong HashStreams(MCMCSampler sampler, bool includeChains = true)
        {
            ulong hash = 14695981039346656037UL;
            void Add(double value)
            {
                ulong bits = (ulong)BitConverter.DoubleToInt64Bits(value);
                for (int b = 0; b < 8; b++)
                {
                    hash ^= (bits >> (8 * b)) & 0xFF;
                    hash *= 1099511628211UL;
                }
            }

            foreach (var lists in includeChains ? new[] { sampler.MarkovChains, sampler.Output } : new[] { sampler.Output })
            {
                foreach (var chain in lists)
                {
                    foreach (var state in chain)
                    {
                        foreach (var value in state.Values) Add(value);
                        Add(state.Fitness);
                    }
                }
            }
            return hash;
        }

        /// <summary>
        /// Asserts the reference hash on Windows and, on every platform, that the output reproduces the target moments.
        /// </summary>
        /// <param name="sampler">A sampled sampler.</param>
        /// <param name="referenceHash">The hash of the chains and the output, recorded on Windows.</param>
        /// <param name="outputHash">The hash of the output alone, recorded on Windows and pinned on other platforms;
        /// null when the chains diverge there.</param>
        private static void AssertReferenceStream(MCMCSampler sampler, ulong referenceHash, ulong? outputHash = null)
        {
            if (CLibraryRounding.IsReferencePlatform) Assert.AreEqual(referenceHash, HashStreams(sampler));
            if (outputHash.HasValue) Assert.AreEqual(outputHash.Value, HashStreams(sampler, includeChains: false));
            for (int j = 0; j < Sd.Length; j++)
            {
                int k = j;
                AssertWithinMonteCarloError(sampler, s => s.Values[k], 0d, $"mean of coordinate {k}");
                AssertWithinMonteCarloError(sampler, s => s.Values[k] * s.Values[k], Sd[k] * Sd[k], $"second moment of coordinate {k}");
            }
        }

        /// <summary>
        /// Asserts that the output mean of a function of the draws lies within four Monte Carlo standard errors of its
        /// exact value under the target.
        /// </summary>
        /// <param name="sampler">A sampled sampler.</param>
        /// <param name="f">The function of a draw.</param>
        /// <param name="exact">The expectation of <paramref name="f"/> under the target.</param>
        /// <param name="what">The quantity named in the failure message.</param>
        /// <remarks>
        /// The standard error is the batch-means estimate: each chain's output is split into ten consecutive batches, and
        /// the error of the overall mean is the standard deviation of the batch means divided by the square root of their
        /// number, which accounts for the autocorrelation within a chain.
        /// </remarks>
        private static void AssertWithinMonteCarloError(MCMCSampler sampler, Func<ParameterSet, double> f, double exact, string what)
        {
            const int batchesPerChain = 10;
            var batchMeans = new List<double>();
            foreach (var chain in sampler.Output)
            {
                int size = chain.Count / batchesPerChain;
                for (int b = 0; b < batchesPerChain; b++)
                {
                    double sum = 0;
                    for (int i = b * size; i < (b + 1) * size; i++) sum += f(chain[i]);
                    batchMeans.Add(sum / size);
                }
            }
            double mean = 0, squares = 0;
            foreach (double m in batchMeans) mean += m;
            mean /= batchMeans.Count;
            foreach (double m in batchMeans) squares += (m - mean) * (m - mean);
            double standardError = Math.Sqrt(squares / (batchMeans.Count - 1) / batchMeans.Count);
            Assert.IsLessThanOrEqualTo(4d * standardError, Math.Abs(mean - exact), $"{what}: {mean:R} against {exact:R}, standard error {standardError:R}");
        }

        /// <summary>
        /// NUTS with an analytic gradient reproduces the reference stream.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_AnalyticGradient_ReproducesReferenceStream()
        {
            var sampler = Configure(new NUTS(Priors(), LogLikelihood, gradientFunction: Gradient), 1);
            sampler.Sample();
            AssertReferenceStream(sampler, 0xB3E23D0BE2FB1F5AUL);
        }

        /// <summary>
        /// A combined log-likelihood-and-gradient function that returns the same values reproduces the
        /// analytic-gradient reference stream bit for bit, while the separate log-likelihood is called only
        /// for chain initialization and the separate gradient not at all.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_CombinedValueAndGradient_ReproducesStreamWithOneEvaluationPerStep()
        {
            int likelihoodCalls = 0, gradientCalls = 0, combinedCalls = 0;
            var nuts = new NUTS(Priors(), x => { System.Threading.Interlocked.Increment(ref likelihoodCalls); return LogLikelihood(x); },
                gradientFunction: x => { System.Threading.Interlocked.Increment(ref gradientCalls); return Gradient(x); })
            {
                LogLikelihoodAndGradientFunction = x =>
                {
                    System.Threading.Interlocked.Increment(ref combinedCalls);
                    var values = new double[x.Count];
                    x.CopyTo(values, 0);
                    return (LogLikelihood(values), Gradient(x));
                }
            };
            var sampler = Configure(nuts, 1);
            sampler.Sample();

            AssertReferenceStream(sampler, 0xB3E23D0BE2FB1F5AUL);
            Assert.AreEqual(sampler.InitialIterations, likelihoodCalls, "Only chain initialization may call the separate log-likelihood.");
            Assert.AreEqual(0, gradientCalls);
            Assert.IsGreaterThan(0, combinedCalls);
        }

        /// <summary>
        /// NUTS with the default finite-difference gradient reproduces the reference stream.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_NumericalGradient_ReproducesReferenceStream()
        {
            var sampler = Configure(new NUTS(Priors(), LogLikelihood), 1);
            sampler.Sample();
            AssertReferenceStream(sampler, 0x028968BE5D65BD19UL);
        }

        /// <summary>
        /// A gradient of the wrong length fails the first leapfrog step, in the step-size heuristic, with an
        /// <see cref="ArgumentException"/> that reaches the caller directly.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_WrongLengthGradient_ThrowsArgumentException()
        {
            var sampler = Configure(new NUTS(Priors(), LogLikelihood, gradientFunction: x => new Vector(x.Count - 1)), 1);
            Assert.ThrowsExactly<ArgumentException>(() => sampler.Sample());
        }

        /// <summary>
        /// An exception from the gradient delegate propagates unwrapped from initialization; the step-size
        /// heuristic does not replace it with the configured step size.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_ThrowingGradient_PropagatesFromInitialization()
        {
            int calls = 0;
            var sampler = Configure(new NUTS(Priors(), LogLikelihood, gradientFunction: x =>
            {
                System.Threading.Interlocked.Increment(ref calls);
                throw new InvalidOperationException("gradient failed");
            }), 1);

            var exception = Assert.ThrowsExactly<InvalidOperationException>(() => sampler.Sample());
            Assert.AreEqual("gradient failed", exception.Message);
            Assert.AreEqual(1, calls);
        }

        /// <summary>
        /// ARWMH reproduces the reference stream through both the identity and the adaptive proposal.
        /// </summary>
        /// <remarks>
        /// With five parameters the identity proposal is used for the first 500 transitions of each chain
        /// and on about 5% of the later ones, and the thinning interval of 10 gives 10,000 transitions per
        /// chain, so both proposal branches are exercised.
        /// </remarks>
        [TestMethod]
        public void Test_ARWMH_ReproducesReferenceStream()
        {
            var sampler = Configure(new ARWMH(Priors(), LogLikelihood), 10);
            sampler.Sample();
            AssertReferenceStream(sampler, 0x5FA307CAFB91C514UL, 0x98C5098077CDD1DDUL);
        }
    }
}
