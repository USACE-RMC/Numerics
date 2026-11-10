using BenchmarkDotNet.Attributes;
using Numerics;

namespace Benchmarks_Numerics.Utilities
{
    /// <summary>
    /// <see cref="Tools.LogSumExp(IList{double})"/> and <see cref="Tools.LogSumExp(double[,], int)"/> against the direct
    /// evaluation max + log(Σ exp(x − max)) over the same values.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class LogSumExpBenchmarks
    {
        private double[] _values = null!;
        private double[,] _matrix = null!;

        /// <summary>The number of rows of the matrix; the list holds rows × 100 values.</summary>
        [Params(100, 1000)]
        public int Rows { get; set; }

        /// <summary>
        /// Prepares normal log values with standard deviation 5, as a list and as a matrix of 100 columns.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(1);
            _matrix = new double[Rows, 100];
            _values = new double[Rows * 100];
            for (int r = 0; r < Rows; r++)
                for (int c = 0; c < 100; c++)
                {
                    double u1 = 1.0 - rng.NextDouble(), u2 = rng.NextDouble();
                    double z = 5 * Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
                    _matrix[r, c] = z;
                    _values[r * 100 + c] = z;
                }
        }

        /// <summary>The direct evaluation max + log(Σ exp(x − max)).</summary>
        /// <returns>The log-sum-exp of the list.</returns>
        [Benchmark(Baseline = true)]
        public double DirectList()
        {
            double max = double.NegativeInfinity;
            for (int i = 0; i < _values.Length; i++) if (_values[i] > max) max = _values[i];
            double sum = 0;
            for (int i = 0; i < _values.Length; i++) sum += Math.Exp(_values[i] - max);
            return max + Math.Log(sum);
        }

        /// <summary>The list form.</summary>
        /// <returns>The log-sum-exp of the list.</returns>
        [Benchmark]
        public double List() => Tools.LogSumExp(_values);

        /// <summary>One value per column.</summary>
        /// <returns>The column log-sum-exps.</returns>
        [Benchmark]
        public double[] MatrixAxis0() => Tools.LogSumExp(_matrix, 0);

        /// <summary>One value per row.</summary>
        /// <returns>The row log-sum-exps.</returns>
        [Benchmark]
        public double[] MatrixAxis1() => Tools.LogSumExp(_matrix, 1);
    }
}
