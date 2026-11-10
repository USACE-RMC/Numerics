using BenchmarkDotNet.Attributes;
using Numerics.Data.Statistics;

namespace Benchmarks_Numerics.Statistics
{
    /// <summary>
    /// The weighted Hazen percentile of unsorted data at 99 unordered probabilities: one call of the list overload against
    /// one scalar call per probability.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class WeightedPercentileListBenchmarks
    {
        private double[] _x = null!, _w = null!, _k = null!;

        /// <summary>The number of weighted points.</summary>
        [Params(100, 2000)]
        public int N { get; set; }

        /// <summary>
        /// Prepares unsorted values with random non-negative weights and the probabilities 0.01 to 0.99 in shuffled order.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(1);
            _x = new double[N];
            _w = new double[N];
            for (int i = 0; i < N; i++)
            {
                _x[i] = rng.NextDouble() * 100d;
                _w[i] = rng.NextDouble();
            }
            _k = new double[99];
            for (int j = 0; j < 99; j++) _k[j] = (j + 1) / 100d;
            for (int j = 98; j > 0; j--)
            {
                int swap = rng.Next(j + 1);
                (_k[j], _k[swap]) = (_k[swap], _k[j]);
            }
        }

        /// <summary>One scalar call per probability, each validating and sorting the sample.</summary>
        /// <returns>The percentiles.</returns>
        [Benchmark(Baseline = true)]
        public double[] ScalarPerProbability()
        {
            var result = new double[_k.Length];
            for (int j = 0; j < _k.Length; j++)
                result[j] = Numerics.Data.Statistics.Statistics.Percentile(_x, _k[j], _w, PlottingPositions.PlottingPostionType.Hazen);
            return result;
        }

        /// <summary>The list overload: one validation, one sort and one walk.</summary>
        /// <returns>The percentiles.</returns>
        [Benchmark]
        public double[] ListOverload()
        {
            return Numerics.Data.Statistics.Statistics.Percentile(_x, _k, _w, PlottingPositions.PlottingPostionType.Hazen);
        }
    }
}
