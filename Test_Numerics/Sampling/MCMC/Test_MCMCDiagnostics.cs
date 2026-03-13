/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
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
            var rng1 = new Random(42);
            var rng2 = new Random(123);
            var rng3 = new Random(456);
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
            Assert.IsGreaterThan(rhatNoWarmup[0], 1.1, $"R-hat without warmup should be > 1.1, got {rhatNoWarmup[0]}");

            // With warmup=50, R-hat should be close to 1.0 (converged portion only)
            var rhatWithWarmup = MCMCDiagnostics.GelmanRubin(chains, 50);
            Assert.IsLessThan(rhatWithWarmup[0], 1.1, $"R-hat with warmup=50 should be < 1.1, got {rhatWithWarmup[0]}");

            // Warmup should improve R-hat (make it closer to 1.0)
            Assert.IsLessThan(rhatWithWarmup[0], rhatNoWarmup[0],
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
