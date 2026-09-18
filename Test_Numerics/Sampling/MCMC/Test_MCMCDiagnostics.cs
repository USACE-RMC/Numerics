using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit tests for MCMCDiagnostics, particularly GelmanRubin R-hat with warmup.
    /// </summary>
    /// <remarks>Reference values use R <c>posterior</c> 1.7.0.</remarks>
    [TestClass]
    public class Test_MCMCDiagnostics
    {
        /// <summary>
        /// Verify Gelman-Rubin R-hat computation with warmup correctly divides by
        /// (N - warmupIterations) instead of N. For chains from the same distribution,
        /// R-hat should be close to 1.0.
        /// </summary>
        [TestMethod]
        public void Test_GelmanRubin_WithWarmup()
        {
            // Create 3 chains of 200 samples each from similar distributions.
            // Use deterministic sequences for reproducibility.
            var rng1 = new MersenneTwister(42);
            var rng2 = new MersenneTwister(123);
            var rng3 = new MersenneTwister(456);
            int chainLength = 200;

            var chain1 = new List<ParameterSet>();
            var chain2 = new List<ParameterSet>();
            var chain3 = new List<ParameterSet>();

            for (int i = 0; i < chainLength; i++)
            {
                // Simulate chains that start dispersed but converge
                double drift1 = i < 50 ? 5.0 : 0.0;
                double drift2 = i < 50 ? -5.0 : 0.0;
                double drift3 = i < 50 ? 3.0 : 0.0;

                chain1.Add(new ParameterSet(new[] { rng1.NextDouble() * 2 - 1 + drift1 }, 0));
                chain2.Add(new ParameterSet(new[] { rng2.NextDouble() * 2 - 1 + drift2 }, 0));
                chain3.Add(new ParameterSet(new[] { rng3.NextDouble() * 2 - 1 + drift3 }, 0));
            }

            var chains = new List<List<ParameterSet>> { chain1, chain2, chain3 };

            // Without warmup, R-hat should be high (chains have different initial distributions)
            var rhatNoWarmup = MCMCDiagnostics.GelmanRubin(chains, 0);
            Assert.IsGreaterThan(1.1, rhatNoWarmup[0], $"R-hat without warmup should be > 1.1, got {rhatNoWarmup[0]}");

            // With warmup=50, R-hat should be close to 1.0 (converged portion only)
            var rhatWithWarmup = MCMCDiagnostics.GelmanRubin(chains, 50);
            Assert.IsLessThan(1.1, rhatWithWarmup[0], $"R-hat with warmup=50 should be < 1.1, got {rhatWithWarmup[0]}");

            // Warmup should improve R-hat (make it closer to 1.0)
            Assert.IsLessThan(rhatNoWarmup[0], rhatWithWarmup[0],
                "R-hat with warmup should be closer to 1.0 than without warmup");
        }

        /// <summary>
        /// Verifies rank-normalized R-hat and conservative ESS against R posterior 1.7.0.
        /// </summary>
        /// <remarks>
        /// R reference: rhat=0.98937039497892809, bulk ESS=616.50943111983349,
        /// lower-tail ESS=319.75259831061169, and upper-tail/scalar ESS=306.3117099147583.
        /// </remarks>
        [TestMethod]
        public void Test_ModernDiagnostics_MatchRPosteriorReference()
        {
            List<List<ParameterSet>> chains = CreateModuloFixture();

            double actualRhat = MCMCDiagnostics.GelmanRubin(chains)[0];
            double actualEss = MCMCDiagnostics.EffectiveSampleSize(chains, out double[][,] averageAcf)[0];

            Assert.AreEqual(0.98937039497892809, actualRhat, 1E-10,
                "Rank-normalized split/folded R-hat must match R posterior 1.7.0.");
            Assert.AreEqual(306.3117099147583, actualEss, 1E-8,
                "The existing scalar ESS must equal min(bulk, lower-tail, upper-tail) from R posterior.");
            Assert.AreEqual(51, averageAcf[0].GetLength(0));
            Assert.AreEqual(2, averageAcf[0].GetLength(1));
            Assert.AreEqual(1d, averageAcf[0][0, 1], 1E-12,
                "The existing original-scale ACF plotting output must be preserved.");
        }

        /// <summary>
        /// Verifies folded rank-normalized R-hat detects equal-center scale disagreement.
        /// </summary>
        /// <remarks>R posterior 1.7.0 returns 1.2810140532084646 for this fixture.</remarks>
        [TestMethod]
        public void Test_GelmanRubin_FoldedRanksDetectScaleMismatch()
        {
            double[] scales = { 0.5d, 0.5d, 2d, 2d };
            var chains = new List<List<ParameterSet>>();
            for (int chainIndex = 0; chainIndex < scales.Length; chainIndex++)
            {
                var chain = new List<ParameterSet>();
                for (int iteration = 1; iteration <= 100; iteration++)
                {
                    double centered = ((iteration * 37 + chainIndex * 13) % 101 - 50) / 10d;
                    chain.Add(new ParameterSet(new[] { scales[chainIndex] * centered }, 0d));
                }
                chains.Add(chain);
            }

            double actual = MCMCDiagnostics.GelmanRubin(chains)[0];

            Assert.AreEqual(1.2810140532084646, actual, 1E-10);
            Assert.IsGreaterThan(1.01, actual,
                "Folded R-hat must flag scale disagreement even when chain centers agree.");
        }

        /// <summary>
        /// Verifies modern diagnostic edge and compatibility behavior.
        /// </summary>
        [TestMethod]
        public void Test_ModernDiagnostics_EdgeCases()
        {
            var constantChains = new List<List<ParameterSet>>();
            for (int chainIndex = 0; chainIndex < 4; chainIndex++)
            {
                var chain = new List<ParameterSet>();
                for (int iteration = 0; iteration < 20; iteration++)
                    chain.Add(new ParameterSet(new[] { 1d }, 0d));
                constantChains.Add(chain);
            }

            Assert.IsTrue(double.IsNaN(MCMCDiagnostics.GelmanRubin(constantChains)[0]));
            Assert.IsTrue(double.IsNaN(
                MCMCDiagnostics.EffectiveSampleSize(constantChains, out _)[0]));
            Assert.IsTrue(double.IsNaN(MCMCDiagnostics.EffectiveSampleSize(new[] { 1d, 2d, 3d })));

            List<List<ParameterSet>> fixture = CreateModuloFixture();
            double original = MCMCDiagnostics.EffectiveSampleSize(fixture, out _)[0];
            fixture.Reverse();
            double permuted = MCMCDiagnostics.EffectiveSampleSize(fixture, out _)[0];
            Assert.AreEqual(original, permuted, 1E-10,
                "ESS must be invariant to chain ordering.");
        }

        /// <summary>
        /// Verify GelmanRubin handles edge cases.
        /// </summary>
        [TestMethod]
        public void Test_GelmanRubin_EdgeCases()
        {
            // Single chain should return NaN
            var chain = new List<ParameterSet>();
            for (int i = 0; i < 10; i++)
                chain.Add(new ParameterSet(new[] { 1.0 }, 0));
            var singleChain = new List<List<ParameterSet>> { chain };
            var result = MCMCDiagnostics.GelmanRubin(singleChain);
            Assert.IsTrue(double.IsNaN(result[0]), "Single chain should return NaN");
        }

        /// <summary>
        /// Creates the compact deterministic fixture evaluated by R posterior 1.7.0.
        /// </summary>
        /// <returns>Four 64-draw chains containing one parameter.</returns>
        private static List<List<ParameterSet>> CreateModuloFixture()
        {
            var chains = new List<List<ParameterSet>>();
            for (int chainIndex = 0; chainIndex < 4; chainIndex++)
            {
                var chain = new List<ParameterSet>();
                for (int iteration = 1; iteration <= 64; iteration++)
                {
                    double value = ((iteration * 17 + chainIndex * 11) % 31) / 10d +
                        chainIndex * 0.05d;
                    chain.Add(new ParameterSet(new[] { value }, 0d));
                }
                chains.Add(chain);
            }
            return chains;
        }
    }
}
