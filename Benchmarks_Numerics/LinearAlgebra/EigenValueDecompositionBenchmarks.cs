using BenchmarkDotNet.Attributes;
using Numerics.Mathematics.LinearAlgebra;

namespace Benchmarks_Numerics.LinearAlgebra
{
    /// <summary>
    /// The symmetric eigen-decomposition by Householder tridiagonalization and implicit QL iterations, on random Gram
    /// matrices of the orders timed against numpy.linalg.eigh in docs/mathematics/linear-algebra.md.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class EigenValueDecompositionBenchmarks
    {
        private Matrix _a = null!;

        /// <summary>
        /// The order of the matrix.
        /// </summary>
        [Params(10, 50, 100, 200, 300, 500)]
        public int N { get; set; }

        /// <summary>
        /// Forms the Gram matrix B·Bᵀ of an N × N matrix B with entries uniform on (−0.5, 0.5) from a seeded generator.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(N);
            var b = new double[N, N];
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++) b[i, j] = rng.NextDouble() - 0.5;
            _a = new Matrix(N);
            for (int i = 0; i < N; i++)
                for (int j = 0; j <= i; j++)
                {
                    double s = 0d;
                    for (int k = 0; k < N; k++) s += b[i, k] * b[j, k];
                    _a[i, j] = s;
                    _a[j, i] = s;
                }
        }

        /// <summary>
        /// Decomposes the matrix with <see cref="EigenSolverMethod.QL"/>.
        /// </summary>
        [Benchmark]
        public double QL() => new EigenValueDecomposition(_a, EigenSolverMethod.QL).EigenValues[0];
    }
}
