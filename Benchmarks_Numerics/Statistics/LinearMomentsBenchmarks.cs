using BenchmarkDotNet.Attributes;

namespace Benchmarks_Numerics.Statistics
{
    /// <summary>
    /// The sample L-moments up to τ5 and the four-value member against the sequential four-moment accumulation.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class LinearMomentsBenchmarks
    {
        private double[] _x = null!;

        /// <summary>The sample size.</summary>
        [Params(50, 1000, 100000)]
        public int N { get; set; }

        /// <summary>
        /// Prepares lognormal-like values.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(1);
            _x = new double[N];
            for (int i = 0; i < N; i++) _x[i] = Math.Exp(3d + rng.NextDouble());
        }

        /// <summary>The sequential four-moment accumulation written out, without the shared evaluation.</summary>
        /// <returns>The L-moments.</returns>
        [Benchmark(Baseline = true)]
        public double[] SequentialFourMoments()
        {
            var sorted = (double[])_x.Clone();
            Array.Sort(sorted);
            double N = sorted.Length;
            double B0 = 0, B1 = 0, B2 = 0, B3 = 0;
            for (int i = 1; i <= N; i++)
            {
                double di = i;
                B0 += sorted[i - 1];
                if (i > 1) B1 += (i - 1) / (N - 1) * sorted[i - 1];
                if (i > 2) B2 += (di - 2) * (di - 1) / ((N - 2) * (N - 1)) * sorted[i - 1];
                if (i > 3) B3 += (di - 3) * (di - 2) * (di - 1) / ((N - 3) * (N - 2) * (N - 1)) * sorted[i - 1];
            }
            B0 /= N; B1 /= N; B2 /= N; B3 /= N;
            return new[] { B0, 2 * B1 - B0, 2 * (3 * B2 - B0) / (2 * B1 - B0) - 3, 5 * (2 * (2 * B3 - 3 * B2) + B0) / (2 * B1 - B0) + 6 };
        }

        /// <summary>λ1, λ2, τ3 and τ4.</summary>
        /// <returns>The L-moments.</returns>
        [Benchmark]
        public double[] FourMoments()
        {
            return Numerics.Data.Statistics.Statistics.LinearMoments(_x);
        }

        /// <summary>λ1, λ2 and τ3 to τ5.</summary>
        /// <returns>The L-moments.</returns>
        [Benchmark]
        public double[] FiveMoments()
        {
            return Numerics.Data.Statistics.Statistics.LinearMoments(_x, 5);
        }
    }
}
