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
    }
}
