using BenchmarkDotNet.Attributes;
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;

namespace Benchmarks_Numerics.LinearAlgebra
{
    /// <summary>
    /// Vectorized matrix products, Gram matrix, transpose, Cholesky and LU against scalar reference loops that perform
    /// the same operations in the same order.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class MatrixProductBenchmarks
    {
        private Matrix _a = null!, _b = null!, _f = null!, _spd = null!, _general = null!, _square = null!;

        /// <summary>
        /// Prepares the matrices: 200 × 300 by 300 × 200 products, a 2,000 × 300 matrix for its Gram matrix, a 1,000 ×
        /// 1,000 square matrix, and 300 × 300 matrices for the factorizations.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(1);
            Matrix Make(int r, int c)
            {
                var m = new Matrix(r, c);
                for (int i = 0; i < r; i++) for (int j = 0; j < c; j++) m[i, j] = rng.NextDouble() - 0.5;
                return m;
            }
            _a = Make(200, 300);
            _b = Make(300, 200);
            _f = Make(2000, 300);
            _square = Make(1000, 1000);
            _general = Make(300, 300);
            var x = Make(400, 300);
            _spd = x.TransposeMultiply(x);
        }

        /// <summary>Scalar matrix product.</summary>
        [Benchmark(Baseline = true)]
        public double MultiplyScalar()
        {
            int rows = _a.NumberOfRows, inner = _a.NumberOfColumns, cols = _b.NumberOfColumns;
            var a = _a.Array;
            var b = _b.Array;
            var c = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                {
                    double sum = 0.0d;
                    for (int k = 0; k < inner; k++) sum += a[i, k] * b[k, j];
                    c[i, j] = sum;
                }
            return c[0, 0];
        }

        /// <summary>Vectorized matrix product.</summary>
        [Benchmark]
        public double Multiply() => _a.Multiply(_b)[0, 0];

        /// <summary>FᵀF as an explicit transpose followed by the scalar product.</summary>
        [Benchmark]
        public double GramScalar()
        {
            var t = new double[_f.NumberOfColumns, _f.NumberOfRows];
            var f = _f.Array;
            for (int i = 0; i < _f.NumberOfRows; i++) for (int j = 0; j < _f.NumberOfColumns; j++) t[j, i] = f[i, j];
            int rows = t.GetLength(0), inner = t.GetLength(1), cols = f.GetLength(1);
            var c = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                {
                    double sum = 0.0d;
                    for (int k = 0; k < inner; k++) sum += t[i, k] * f[k, j];
                    c[i, j] = sum;
                }
            return c[0, 0];
        }

        /// <summary>FᵀF with <see cref="Matrix.TransposeMultiply(Matrix)"/>.</summary>
        [Benchmark]
        public double Gram() => _f.TransposeMultiply(_f)[0, 0];

        /// <summary>Element-by-element transpose.</summary>
        [Benchmark]
        public double TransposeScalar()
        {
            var t = new Matrix(_square.NumberOfColumns, _square.NumberOfRows);
            for (int i = 0; i < _square.NumberOfRows; i++)
                for (int j = 0; j < _square.NumberOfColumns; j++) t[j, i] = _square[i, j];
            return t[0, 1];
        }

        /// <summary>Tiled transpose.</summary>
        [Benchmark]
        public double Transpose() => _square.Transpose()[0, 1];

        /// <summary>Element-by-element copy constructor.</summary>
        [Benchmark]
        public double CopyScalar()
        {
            var a = _square.Array;
            var c = new double[a.GetLength(0), a.GetLength(1)];
            for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) c[i, j] = a[i, j];
            return c[0, 0];
        }

        /// <summary>Block copy constructor.</summary>
        [Benchmark]
        public double Copy() => _square.Clone()[0, 0];

        /// <summary>Scalar Cholesky loop (flat row-major, inner products over k descending).</summary>
        [Benchmark]
        public double CholeskyScalar()
        {
            int n = _spd.NumberOfRows;
            var w = new double[n * n];
            var s = _spd.Array;
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) w[i * n + j] = s[i, j];
            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                {
                    double sum = w[i * n + j];
                    for (int k = i - 1; k >= 0; k--) sum -= w[i * n + k] * w[j * n + k];
                    if (i == j) w[i * n + i] = Math.Sqrt(sum);
                    else w[j * n + i] = sum / w[i * n + i];
                }
            return w[n * n - 1];
        }

        /// <summary>Cholesky factorization.</summary>
        [Benchmark]
        public double Cholesky() => new CholeskyDecomposition(_spd).L[0, 0];

        /// <summary>Scalar LU elimination with partial pivoting.</summary>
        [Benchmark]
        public double LUScalar()
        {
            int n = _general.NumberOfRows;
            var lu = (double[,])_general.Array.Clone();
            var vv = new double[n];
            for (int i = 0; i < n; i++)
            {
                double big = 0;
                for (int j = 0; j < n; j++) big = Math.Max(big, Math.Abs(lu[i, j]));
                vv[i] = 1.0 / big;
            }
            for (int k = 0; k < n; k++)
            {
                double big = 0; int imax = k;
                for (int i = k; i < n; i++) { double t = vv[i] * Math.Abs(lu[i, k]); if (t > big) { big = t; imax = i; } }
                if (k != imax) { for (int j = 0; j < n; j++) (lu[imax, j], lu[k, j]) = (lu[k, j], lu[imax, j]); vv[imax] = vv[k]; }
                if (lu[k, k] == 0.0) lu[k, k] = 1e-40;
                for (int i = k + 1; i < n; i++)
                {
                    lu[i, k] /= lu[k, k];
                    for (int j = k + 1; j < n; j++) lu[i, j] -= lu[i, k] * lu[k, j];
                }
            }
            return lu[n - 1, n - 1];
        }

        /// <summary>LU decomposition.</summary>
        [Benchmark]
        public double LU() => new LUDecomposition(_general).LU[0, 0];
    }

    /// <summary>
    /// Running covariance updates and correlation matrices against matrix-expression and scalar reference versions.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class StatisticsBenchmarks
    {
        private double[][] _vectors = null!;
        private double[,] _samples = null!;

        /// <summary>
        /// The number of variables in each pushed vector.
        /// </summary>
        [Params(10, 50)]
        public int Size { get; set; }

        /// <summary>
        /// Prepares 1,000 vectors to push and a 2,000 × 20 sample.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(2);
            _vectors = new double[1000][];
            for (int i = 0; i < _vectors.Length; i++)
            {
                _vectors[i] = new double[Size];
                for (int j = 0; j < Size; j++) _vectors[i][j] = rng.NextDouble();
            }
            _samples = new double[2000, 20];
            for (int i = 0; i < 2000; i++) for (int j = 0; j < 20; j++) _samples[i, j] = rng.NextDouble() + j;
        }

        /// <summary>The matrix-expression Welford update for 1,000 vectors.</summary>
        [Benchmark(Baseline = true)]
        public double PushExpression()
        {
            var mean = new Matrix(Size, 1);
            var covariance = Matrix.Identity(Size);
            for (int n = 1; n <= _vectors.Length; n++)
            {
                var x = new Matrix(_vectors[n - 1]);
                var oldMean = mean;
                mean += n == 1 ? x : (x - mean) * (1d / n);
                covariance += (x - oldMean) * Matrix.Transpose(x - mean);
            }
            return covariance[0, 0];
        }

        /// <summary><see cref="RunningCovarianceMatrix.Push"/> for 1,000 vectors.</summary>
        [Benchmark]
        public double Push()
        {
            var r = new RunningCovarianceMatrix(Size);
            foreach (var v in _vectors) r.Push(v);
            return r.Covariance[0, 0];
        }

        /// <summary>Scalar Pearson correlation matrix of the 2,000 × 20 sample.</summary>
        [Benchmark]
        public double PearsonScalar()
        {
            var s = _samples;
            int n = s.GetLength(0), p = s.GetLength(1);
            var means = new double[p];
            for (int j = 0; j < p; j++) for (int i = 0; i < n; i++) means[j] += s[i, j];
            for (int j = 0; j < p; j++) means[j] /= n;
            var ss = new double[p];
            var cov = new double[p, p];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < p; j++)
                {
                    double dx = s[i, j] - means[j];
                    ss[j] += dx * dx;
                    for (int k = j; k < p; k++) cov[j, k] += dx * (s[i, k] - means[k]);
                }
            return cov[0, p - 1] / Math.Sqrt(ss[0] * ss[p - 1]);
        }

        /// <summary><see cref="Correlation.Pearson(double[,])"/> of the 2,000 × 20 sample.</summary>
        [Benchmark]
        public double Pearson() => Correlation.Pearson(_samples)[0, 19];
    }
}
