using BenchmarkDotNet.Attributes;
using Numerics.Data.Statistics;

namespace Benchmarks_Numerics.Statistics
{
    /// <summary>
    /// The weighted Hazen percentile on sorted arrays against a plain sequential loop that performs the same operations in
    /// the same order.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class WeightedPlottingPositionPercentileBenchmarks
    {
        private double[] _x = null!, _w = null!;

        /// <summary>The number of weighted points.</summary>
        [Params(50, 500)]
        public int N { get; set; }

        /// <summary>
        /// Prepares ascending values with random non-negative weights.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(1);
            _x = new double[N];
            _w = new double[N];
            for (int i = 0; i < N; i++)
            {
                _x[i] = i + rng.NextDouble();
                _w[i] = rng.NextDouble();
            }
        }

        /// <summary>The plain sequential loop: positions (c_(i−1) + w_i / 2) / W and linear interpolation.</summary>
        /// <returns>The 90th percentile.</returns>
        [Benchmark(Baseline = true)]
        public double SequentialLoop()
        {
            double k = 0.9, total = 0d;
            for (int i = 0; i < _w.Length; i++) total += _w[i];
            double before = 0d, previousPosition = double.NaN, previousValue = double.NaN;
            for (int i = 0; i < _x.Length; i++)
            {
                double position = (before + 0.5 * _w[i]) / total;
                if (k <= position)
                {
                    if (i == 0) return _x[0];
                    double span = position - previousPosition;
                    if (span <= 0d) return _x[i];
                    return previousValue + (_x[i] - previousValue) * (k - previousPosition) / span;
                }
                previousPosition = position;
                previousValue = _x[i];
                before += _w[i];
            }
            return _x[_x.Length - 1];
        }

        /// <summary>The library member on sorted arrays.</summary>
        /// <returns>The 90th percentile.</returns>
        [Benchmark]
        public double WeightedHazenPercentile()
        {
            return Numerics.Data.Statistics.Statistics.Percentile(_x, 0.9, _w, PlottingPositions.PlottingPostionType.Hazen, true);
        }
    }
}
