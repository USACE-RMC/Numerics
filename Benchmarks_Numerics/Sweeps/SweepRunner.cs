using System.Diagnostics;
using Data.Statistics;
using Data.TimeSeriesAnalysis;
using Mathematics.LinearAlgebra;
using Utilities;

namespace Benchmarks_Numerics.Sweeps
{
    /// <summary>
    /// Runs the full-size randomized and accuracy sweeps whose reduced forms are the Test_Numerics unit tests. The sweep
    /// bodies are the unit-test methods themselves, compiled into this project from the <c>*.Sweep.cs</c> files of
    /// Test_Numerics, so the unit test and the sweep cannot drift apart.
    /// </summary>
    /// <remarks>
    /// Run with <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep [name ...]</c>; with no name every sweep
    /// runs. A failed assertion stops the run with a nonzero exit code.
    /// </remarks>
    public static class SweepRunner
    {
        private static readonly (string Name, Action Run)[] Sweeps =
        {
            ("BlockSeries", () => Test_TimeSeriesBlockSeries.MatchLegacy(60, 900, 3000)),
            ("CyclicJacobi", () => Test_EigenValueDecomposition.CyclicJacobiDocumentedAccuracy(10, 50, 120, 300)),
            ("Unimodality", () => Test_HypothesisTests.UnimodalityTestBootstrap(199)),
            ("Cholesky", () => Test_VectorizedLinearAlgebra.CholeskyBlockedErrorsMatchScalarLoop(16, 64, 200, 400)),
            ("PlottingPosition", () =>
            {
                Test_PlottingPositionPercentile.WeightedListOfProbabilitiesMatchesScalarBitwiseSweep(600);
                Test_PlottingPositionPercentile.BlomAndHazenWeightCloserOnProbabilityGridsSweep(20000);
                Test_PlottingPositionPercentile.OtherFormulasKeepTheirWeightsBitwiseSweep(3000);
                Test_PlottingPositionPercentile.BlomAndHazenWeightCloserToExactRationalSweep(3000);
                Test_PlottingPositionPercentile.WeibullWeightCloserToExactRationalSweep(3000);
                Test_PlottingPositionPercentile.WeightedHazenMatchesSequentialReferenceBitwiseSweep(2000);
            }),
            ("ToolsBitwise", () =>
            {
                Test_ToolsBitwiseReference.Expm1MatchesReferenceBitwiseSweep(1000, 4096, 1000000);
                Test_ToolsBitwiseReference.Log1pMatchesReferenceBitwiseSweep(1000, 4096, 1000000);
                Test_ToolsBitwiseReference.Log1mExpMatchesReferenceBitwiseSweep(1000, 4096, 1000000);
            }),
        };

        /// <summary>
        /// Runs the named sweeps, or every sweep when no name is given.
        /// </summary>
        /// <param name="names">The sweep names, matched without regard to case.</param>
        /// <returns>Zero when every selected sweep passes, otherwise one.</returns>
        public static int Run(IReadOnlyCollection<string> names)
        {
            var unknown = names.Where(n => !Sweeps.Any(s => string.Equals(s.Name, n, StringComparison.OrdinalIgnoreCase))).ToList();
            if (unknown.Count > 0)
            {
                Console.Error.WriteLine($"Unknown sweep: {string.Join(", ", unknown)}. Sweeps: {string.Join(", ", Sweeps.Select(s => s.Name))}.");
                return 1;
            }
            int failures = 0;
            foreach (var (name, run) in Sweeps)
            {
                if (names.Count > 0 && !names.Contains(name, StringComparer.OrdinalIgnoreCase)) continue;
                var watch = Stopwatch.StartNew();
                try
                {
                    run();
                    Console.WriteLine($"{name}: passed in {watch.Elapsed.TotalSeconds:F1} s");
                }
                catch (Exception ex)
                {
                    failures++;
                    Console.WriteLine($"{name}: FAILED after {watch.Elapsed.TotalSeconds:F1} s");
                    Console.WriteLine(ex);
                }
            }
            return failures == 0 ? 0 : 1;
        }
    }
}
