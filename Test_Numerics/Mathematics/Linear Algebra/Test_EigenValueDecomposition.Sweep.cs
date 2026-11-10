#nullable disable
using System;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// The documented-accuracy measurement of the cyclic Jacobi eigenvalue decomposition, shared by the unit test and
    /// the full sweep in Benchmarks_Numerics.
    /// </summary>
    public partial class Test_EigenValueDecomposition
    {
        /// <summary>
        /// Asserts the documented residual, orthonormality and eigenvalue bounds of the cyclic Jacobi decomposition of a
        /// random Gram matrix and of evenly spaced and graded known spectra at each order.
        /// </summary>
        /// <param name="orders">The matrix orders.</param>
        public static void CyclicJacobiDocumentedAccuracy(params int[] orders)
        {
            foreach (int n in orders)
            {
                CheckAccuracy(RandomGram(n, 7 + n), null, 1E-14);
                var spaced = new double[n];
                var graded = new double[n];
                for (int i = 0; i < n; i++)
                {
                    spaced[i] = 1d + i;
                    graded[i] = System.Math.Pow(10d, -12d * i / (n - 1));
                }
                CheckAccuracy(KnownSpectrum(spaced, 100 + n), spaced, 2E-14);
                CheckAccuracy(KnownSpectrum(graded, 200 + n), graded, 2E-12);
            }
        }

        /// <summary>
        /// Returns the symmetric matrix B * B^T for B with uniform(-0.5, 0.5) entries from a seeded generator.
        /// </summary>
        private static Matrix RandomGram(int n, int seed)
        {
            var rng = new System.Random(seed);
            var B = new Matrix(n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    B[i, j] = rng.NextDouble() - 0.5;
            var S = B * Matrix.Transpose(B);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < i; j++)
                    S[i, j] = S[j, i];
            return S;
        }

        /// <summary>
        /// Asserts the residual, orthonormality and (when a spectrum is given) eigenvalue bounds of the decomposition of
        /// <paramref name="A"/>, and writes the measured values to the test output.
        /// </summary>
        private static void CheckAccuracy(Matrix A, double[] spectrum, double eigenvalueTolerance)
        {
            int n = A.NumberOfRows;
            var eig = new EigenValueDecomposition(A);
            Assert.IsTrue(eig.IsConverged);
            var V = eig.EigenVectors;
            double frobenius = 0d, largest = 0d;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) frobenius += A[i, j] * A[i, j];
            frobenius = System.Math.Sqrt(frobenius);
            for (int i = 0; i < n; i++) largest = System.Math.Max(largest, System.Math.Abs(eig.EigenValues[i]));
            double residual = 0d, orthogonality = 0d;
            for (int c = 0; c < n; c++)
            {
                double sum = 0d;
                for (int i = 0; i < n; i++)
                {
                    double av = 0d;
                    for (int k = 0; k < n; k++) av += A[i, k] * V[k, c];
                    double r = av - eig.EigenValues[c] * V[i, c];
                    sum += r * r;
                }
                residual = System.Math.Max(residual, System.Math.Sqrt(sum) / frobenius);
                for (int d = 0; d <= c; d++)
                {
                    double dot = 0d;
                    for (int k = 0; k < n; k++) dot += V[k, c] * V[k, d];
                    orthogonality = System.Math.Max(orthogonality, System.Math.Abs(dot - (c == d ? 1d : 0d)));
                }
            }
            double eigenvalueError = 0d;
            if (spectrum != null)
            {
                var computed = eig.EigenValues.ToArray();
                var expected = (double[])spectrum.Clone();
                System.Array.Sort(computed);
                System.Array.Sort(expected);
                for (int i = 0; i < n; i++)
                    eigenvalueError = System.Math.Max(eigenvalueError, System.Math.Abs(computed[i] - expected[i]) / largest);
            }
            System.Console.WriteLine($"n = {n}: residual {residual:E2}, orthogonality {orthogonality:E2}, eigenvalues {eigenvalueError:E2}");
            Assert.IsLessThan(1E-11, residual, $"n = {n}");
            Assert.IsLessThan(5E-14, orthogonality, $"n = {n}");
            Assert.IsLessThan(eigenvalueTolerance, eigenvalueError, $"n = {n}");
        }

        /// <summary>
        /// Returns Q Λ Qᵀ with Q the product of three Householder reflections I − 2vvᵀ/(vᵀv), v uniform(−0.5, 0.5) from a
        /// seeded generator, symmetrized exactly.
        /// </summary>
        private static Matrix KnownSpectrum(double[] spectrum, int seed)
        {
            int n = spectrum.Length;
            var rng = new System.Random(seed);
            var Q = new Matrix(n);
            for (int i = 0; i < n; i++) Q[i, i] = 1d;
            for (int h = 0; h < 3; h++)
            {
                var v = new double[n];
                double norm = 0d;
                for (int i = 0; i < n; i++) { v[i] = rng.NextDouble() - 0.5; norm += v[i] * v[i]; }
                // Q ← Q (I − 2vvᵀ/norm)
                for (int i = 0; i < n; i++)
                {
                    double qv = 0d;
                    for (int k = 0; k < n; k++) qv += Q[i, k] * v[k];
                    for (int k = 0; k < n; k++) Q[i, k] -= 2d * qv * v[k] / norm;
                }
            }
            var A = new Matrix(n);
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++)
                {
                    double s = 0d;
                    for (int k = 0; k < n; k++) s += Q[i, k] * spectrum[k] * Q[j, k];
                    A[i, j] = s;
                    A[j, i] = s;
                }
            return A;
        }
    }
}
