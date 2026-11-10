using BenchmarkDotNet.Attributes;
using Numerics.Distributions;

namespace Benchmarks_Numerics.Distributions
{
    /// <summary>
    /// The per-element-parameter static Poisson and binomial members against a loop that creates one distribution per
    /// element, on 100,000 elements.
    /// </summary>
    /// <remarks>
    /// The inputs match the Python comparison with scipy.stats: Poisson rates uniform on [0.5, 50] with uniform
    /// probabilities and counts drawn near each rate, and binomial success probabilities uniform on [0.01, 0.99] with
    /// 1 to 100 trials.
    /// </remarks>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class DiscreteElementwiseBenchmarks
    {
        private const int N = 100_000;
        private double[] _rates = null!, _probabilities = null!, _counts = null!, _success = null!, _successes = null!;
        private int[] _trials = null!;

        /// <summary>
        /// Prepares the seeded inputs.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(1);
            _rates = new double[N];
            _probabilities = new double[N];
            _counts = new double[N];
            _success = new double[N];
            _trials = new int[N];
            _successes = new double[N];
            for (int i = 0; i < N; i++)
            {
                _rates[i] = 0.5 + 49.5 * rng.NextDouble();
                _probabilities[i] = rng.NextDouble();
                _counts[i] = Math.Floor(_rates[i] * 2 * rng.NextDouble());
                _success[i] = 0.01 + 0.98 * rng.NextDouble();
                _trials[i] = 1 + rng.Next(100);
                _successes[i] = rng.Next(_trials[i] + 1);
            }
        }

        /// <summary>One <see cref="Poisson"/> instance per element.</summary>
        /// <returns>The quantiles.</returns>
        [Benchmark(Baseline = true)]
        public double[] PoissonInverseCDFInstanceLoop()
        {
            var result = new double[N];
            for (int i = 0; i < N; i++) result[i] = new Poisson(_rates[i]).InverseCDF(_probabilities[i]);
            return result;
        }

        /// <summary>The static per-element quantile.</summary>
        /// <returns>The quantiles.</returns>
        [Benchmark]
        public double[] PoissonInverseCDFStatic() => Poisson.InverseCDF(_rates, _probabilities);

        /// <summary>The static per-element log mass.</summary>
        /// <returns>The log masses.</returns>
        [Benchmark]
        public double[] PoissonLogPDFStatic() => Poisson.LogPDF(_rates, _counts);

        /// <summary>The static per-element CDF.</summary>
        /// <returns>The probabilities.</returns>
        [Benchmark]
        public double[] PoissonCDFStatic() => Poisson.CDF(_rates, _counts);

        /// <summary>The static per-element binomial log mass.</summary>
        /// <returns>The log masses.</returns>
        [Benchmark]
        public double[] BinomialLogPDFStatic() => Binomial.LogPDF(_success, _trials, _successes);

        /// <summary>The static per-element binomial CDF.</summary>
        /// <returns>The probabilities.</returns>
        [Benchmark]
        public double[] BinomialCDFStatic() => Binomial.CDF(_success, _trials, _successes);
    }
}
