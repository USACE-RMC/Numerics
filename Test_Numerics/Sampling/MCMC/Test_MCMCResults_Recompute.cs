using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit tests for <see cref="MCMCResults.RecomputeParameterResults"/>. The contract:
    /// recomputing parameter results at a new alpha must update the credible-interval
    /// percentiles (LowerCI/UpperCI) while preserving alpha-independent diagnostics
    /// (Rhat, ESS, Autocorrelation) AND the underlying chain output (Output, MAP).
    /// </summary>
    /// <remarks>
    /// This is the cornerstone of the RMC.BestFit Phase 2 reprocess-don't-clear refactor.
    /// A regression that drops one of the three snapshot/restore lines silently corrupts
    /// every convergence-diagnostic display the moment a user changes CredibleIntervalWidth
    /// on an estimated analysis.
    /// </remarks>
    [TestClass]
    public class Test_MCMCResults_Recompute
    {
        /// <summary>
        /// Builds a minimal <see cref="MCMCResults"/> populated via the
        /// <see cref="MCMCResults(ParameterSet, IList{ParameterSet}, double)"/> constructor.
        /// MarkovChains/AcceptanceRates/MeanLogLikelihood remain null per the constructor's
        /// documented behavior.
        /// </summary>
        private static MCMCResults BuildResults(double alpha)
        {
            // Three parameters (mu, sigma, tau) sampled from synthetic distributions.
            // 1000 samples is enough for smooth percentile estimates.
            var rng = new MersenneTwister(2026);
            var output = new List<ParameterSet>(capacity: 1000);
            for (int i = 0; i < 1000; i++)
            {
                // Param 0: ~ Normal(10, 1)
                double p0 = 10.0 + StandardNormal(rng);
                // Param 1: ~ Normal(2, 0.5)
                double p1 = 2.0 + 0.5 * StandardNormal(rng);
                // Param 2: ~ Uniform(-1, 1)
                double p2 = -1.0 + 2.0 * rng.NextDouble();
                output.Add(new ParameterSet(new[] { p0, p1, p2 }, 0));
            }
            var map = new ParameterSet(new[] { 10.0, 2.0, 0.0 }, 0);
            return new MCMCResults(map, output, alpha);
        }

        /// <summary>Box-Muller transform, single sample.</summary>
        private static double StandardNormal(MersenneTwister rng)
        {
            double u1 = rng.NextDouble();
            double u2 = rng.NextDouble();
            if (u1 < 1e-300) u1 = 1e-300;
            return System.Math.Sqrt(-2.0 * System.Math.Log(u1)) * System.Math.Cos(2.0 * System.Math.PI * u2);
        }

        /// <summary>
        /// Seed each parameter's alpha-independent diagnostics with sentinel values.
        /// These would normally be populated by the <see cref="MCMCResults(MCMCSampler, double)"/>
        /// path; for test purposes we set them directly to verify the snapshot/restore.
        /// </summary>
        private static (double[] rhats, double[] esss, double[][,] acfs) SeedDiagnostics(MCMCResults results)
        {
            int n = results.ParameterResults.Length;
            var rhats = new double[n];
            var esss = new double[n];
            var acfs = new double[n][,];
            for (int i = 0; i < n; i++)
            {
                rhats[i] = 1.0 + 0.01 * (i + 1);                // 1.01, 1.02, 1.03
                esss[i] = 800.0 + 10.0 * i;                     // 800, 810, 820
                acfs[i] = new double[1, 5] { { 1.0, 0.5 - 0.1 * i, 0.25, 0.125, 0.0625 } };
                results.ParameterResults[i].SummaryStatistics.Rhat = rhats[i];
                results.ParameterResults[i].SummaryStatistics.ESS = esss[i];
                results.ParameterResults[i].Autocorrelation = acfs[i];
            }
            return (rhats, esss, acfs);
        }

        [TestMethod]
        public void Recompute_PreservesOutputReference()
        {
            var results = BuildResults(alpha: 0.10);
            var outputRef = results.Output;

            results.RecomputeParameterResults(alpha: 0.05);

            // The Output list must be the same reference — the chain itself is not touched.
            Assert.AreSame(outputRef, results.Output, "Output reference must not change.");
            Assert.AreEqual(1000, results.Output.Count, "Output count must be preserved.");
        }

        [TestMethod]
        public void Recompute_PreservesMAP()
        {
            var results = BuildResults(alpha: 0.10);
            var valuesBefore = (double[])results.MAP.Values.Clone();
            var fitnessBefore = results.MAP.Fitness;

            results.RecomputeParameterResults(alpha: 0.05);

            Assert.AreEqual(valuesBefore.Length, results.MAP.Values.Length, "MAP value count must not change.");
            for (int i = 0; i < valuesBefore.Length; i++)
                Assert.AreEqual(valuesBefore[i], results.MAP.Values[i], 1e-12, $"MAP value [{i}] must not change.");
            Assert.AreEqual(fitnessBefore, results.MAP.Fitness, 1e-12, "MAP fitness must not change.");
        }

        [TestMethod]
        public void Recompute_PreservesRhat()
        {
            var results = BuildResults(alpha: 0.10);
            var (rhats, _, _) = SeedDiagnostics(results);

            results.RecomputeParameterResults(alpha: 0.05);

            for (int i = 0; i < results.ParameterResults.Length; i++)
            {
                Assert.AreEqual(rhats[i], results.ParameterResults[i].SummaryStatistics.Rhat, 1e-12,
                    $"Parameter {i}: Rhat must be preserved across alpha change.");
            }
        }

        [TestMethod]
        public void Recompute_PreservesESS()
        {
            var results = BuildResults(alpha: 0.10);
            var (_, esss, _) = SeedDiagnostics(results);

            results.RecomputeParameterResults(alpha: 0.05);

            for (int i = 0; i < results.ParameterResults.Length; i++)
            {
                Assert.AreEqual(esss[i], results.ParameterResults[i].SummaryStatistics.ESS, 1e-12,
                    $"Parameter {i}: ESS must be preserved across alpha change.");
            }
        }

        [TestMethod]
        public void Recompute_PreservesAutocorrelation()
        {
            var results = BuildResults(alpha: 0.10);
            var (_, _, acfs) = SeedDiagnostics(results);

            results.RecomputeParameterResults(alpha: 0.05);

            for (int i = 0; i < results.ParameterResults.Length; i++)
            {
                var acfAfter = results.ParameterResults[i].Autocorrelation;
                Assert.AreEqual(acfs[i].GetLength(0), acfAfter.GetLength(0));
                Assert.AreEqual(acfs[i].GetLength(1), acfAfter.GetLength(1));
                for (int r = 0; r < acfs[i].GetLength(0); r++)
                {
                    for (int c = 0; c < acfs[i].GetLength(1); c++)
                    {
                        Assert.AreEqual(acfs[i][r, c], acfAfter[r, c], 1e-12,
                            $"Parameter {i} ACF[{r},{c}]: must be preserved across alpha change.");
                    }
                }
            }
        }

        [TestMethod]
        public void Recompute_NarrowsCIWhenAlphaIncreases()
        {
            // alpha=0.10 → 90% CI; alpha=0.20 → 80% CI (narrower band).
            var results = BuildResults(alpha: 0.10);
            var lower90 = results.ParameterResults.Select(p => p.SummaryStatistics.LowerCI).ToArray();
            var upper90 = results.ParameterResults.Select(p => p.SummaryStatistics.UpperCI).ToArray();

            results.RecomputeParameterResults(alpha: 0.20);

            for (int i = 0; i < results.ParameterResults.Length; i++)
            {
                double lower80 = results.ParameterResults[i].SummaryStatistics.LowerCI;
                double upper80 = results.ParameterResults[i].SummaryStatistics.UpperCI;
                Assert.IsTrue(lower80 > lower90[i],
                    $"Parameter {i}: 80% LowerCI ({lower80}) must be > 90% LowerCI ({lower90[i]}) — band narrows when alpha increases.");
                Assert.IsTrue(upper80 < upper90[i],
                    $"Parameter {i}: 80% UpperCI ({upper80}) must be < 90% UpperCI ({upper90[i]}) — band narrows when alpha increases.");
            }
        }

        [TestMethod]
        public void Recompute_WidensCIWhenAlphaDecreases()
        {
            // alpha=0.10 → 90% CI; alpha=0.05 → 95% CI (wider band).
            var results = BuildResults(alpha: 0.10);
            var lower90 = results.ParameterResults.Select(p => p.SummaryStatistics.LowerCI).ToArray();
            var upper90 = results.ParameterResults.Select(p => p.SummaryStatistics.UpperCI).ToArray();

            results.RecomputeParameterResults(alpha: 0.05);

            for (int i = 0; i < results.ParameterResults.Length; i++)
            {
                double lower95 = results.ParameterResults[i].SummaryStatistics.LowerCI;
                double upper95 = results.ParameterResults[i].SummaryStatistics.UpperCI;
                Assert.IsTrue(lower95 < lower90[i],
                    $"Parameter {i}: 95% LowerCI ({lower95}) must be < 90% LowerCI ({lower90[i]}) — band widens when alpha decreases.");
                Assert.IsTrue(upper95 > upper90[i],
                    $"Parameter {i}: 95% UpperCI ({upper95}) must be > 90% UpperCI ({upper90[i]}) — band widens when alpha decreases.");
            }
        }

        [TestMethod]
        public void Recompute_PreservesMeanAndMedian()
        {
            // Mean and Median are alpha-independent — should be preserved exactly.
            var results = BuildResults(alpha: 0.10);
            var meansBefore = results.ParameterResults.Select(p => p.SummaryStatistics.Mean).ToArray();
            var mediansBefore = results.ParameterResults.Select(p => p.SummaryStatistics.Median).ToArray();

            results.RecomputeParameterResults(alpha: 0.05);

            for (int i = 0; i < results.ParameterResults.Length; i++)
            {
                Assert.AreEqual(meansBefore[i], results.ParameterResults[i].SummaryStatistics.Mean, 1e-10,
                    $"Parameter {i}: Mean must be preserved.");
                Assert.AreEqual(mediansBefore[i], results.ParameterResults[i].SummaryStatistics.Median, 1e-10,
                    $"Parameter {i}: Median must be preserved.");
            }
        }

        [TestMethod]
        public void Recompute_OnEmptyResults_DoesNotThrow()
        {
            // ParameterResults == null guard.
            var results = new MCMCResults();

            results.RecomputeParameterResults(alpha: 0.05);

            // Should be a no-op — ParameterResults stays null, no exception.
            Assert.IsNull(results.ParameterResults);
        }
    }
}
