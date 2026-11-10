using BenchmarkDotNet.Attributes;

namespace Benchmarks_Numerics.Integration
{
    /// <summary>
    /// The QUADPACK QAGS body of AdaptiveGaussKronrod against its default local-bisection method on the
    /// comparison cases: Test_Numerics integrands, internal Numerics callers, representative TotalRisk
    /// integrands and classic hard cases. Each run constructs a new integrator, as the call sites do.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class AdaptiveGaussKronrodBenchmarks
    {
        private AdaptiveGaussKronrodCase _case = null!;

        /// <summary>The case names, or the comma-separated subset named by the AGK_CASES environment variable.</summary>
        public static IEnumerable<string> CaseNames
        {
            get
            {
                string? subset = Environment.GetEnvironmentVariable("AGK_CASES");
                var names = AdaptiveGaussKronrodCases.All().Select(c => c.Name);
                return string.IsNullOrEmpty(subset) ? names : names.Where(n => subset.Split(',').Contains(n));
            }
        }

        /// <summary>The case under test.</summary>
        [ParamsSource(nameof(CaseNames))]
        public string Case { get; set; } = "";

        /// <summary>
        /// Selects the case.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            _case = AdaptiveGaussKronrodCases.All().Single(c => c.Name == Case);
        }

        /// <summary>The default local-bisection method.</summary>
        /// <returns>The integral.</returns>
        [Benchmark(Baseline = true)]
        public double Old() => _case.RunOld().Result;

        /// <summary>The QAGS method.</summary>
        /// <returns>The integral.</returns>
        [Benchmark]
        public double New() => _case.RunNew().Result;
    }
}
