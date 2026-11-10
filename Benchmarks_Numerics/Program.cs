using BenchmarkDotNet.Running;

namespace Benchmarks_Numerics
{
    /// <summary>
    /// Entry point for the Numerics BenchmarkDotNet suite. Run with
    /// <c>dotnet run -c Release --project Benchmarks_Numerics -- --filter *</c>.
    /// <c>-- agk-compare &lt;path&gt;</c> writes the AdaptiveGaussKronrod old-versus-new comparison table instead, and
    /// <c>-- sweep [name ...]</c> runs the full-size forms of the reduced Test_Numerics sweeps, and <c>-- compare [name ...]</c>
    /// times the machine-learning methods against the scikit-learn, R and DiceKriging reference timings.
    /// </summary>
    public static class Program
    {
        /// <summary>
        /// Runs the benchmarks selected on the command line.
        /// </summary>
        /// <param name="args">BenchmarkDotNet command-line arguments.</param>
        /// <returns>The process exit code: nonzero when a sweep or comparison fails.</returns>
        public static int Main(string[] args)
        {
            if (args.Length == 2 && args[0] == "agk-compare")
            {
                Integration.AdaptiveGaussKronrodCases.WriteComparison(args[1]);
                return 0;
            }
            if (args.Length >= 1 && args[0] == "compare")
                return MachineLearning.ReferenceComparisons.Run(args.Skip(1).ToArray());
            if (args.Length >= 1 && args[0] == "sweep")
                return Sweeps.SweepRunner.Run(args.Skip(1).ToArray());
            BenchmarkSwitcher.FromAssembly(typeof(Program).Assembly).Run(args);
            return 0;
        }
    }
}
