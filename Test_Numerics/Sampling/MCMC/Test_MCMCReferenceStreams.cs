using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling.MCMC;

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
        /// <returns>The hash.</returns>
        private static ulong HashStreams(MCMCSampler sampler)
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

            foreach (var lists in new[] { sampler.MarkovChains, sampler.Output })
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
        /// NUTS with an analytic gradient reproduces the reference stream.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_AnalyticGradient_ReproducesReferenceStream()
        {
            var sampler = Configure(new NUTS(Priors(), LogLikelihood, gradientFunction: Gradient), 1);
            sampler.Sample();
            Assert.AreEqual(0xB3E23D0BE2FB1F5AUL, HashStreams(sampler));
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

            Assert.AreEqual(0xB3E23D0BE2FB1F5AUL, HashStreams(sampler));
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
            Assert.AreEqual(0x028968BE5D65BD19UL, HashStreams(sampler));
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
            Assert.AreEqual(0x5FA307CAFB91C514UL, HashStreams(sampler));
        }
    }
}
