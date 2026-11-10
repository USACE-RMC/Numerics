#nullable disable
using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// The scalar reference loops and test matrices of the vectorized linear algebra tests, and the blocked Cholesky
    /// error comparison shared by the unit test and the full sweep in Benchmarks_Numerics.
    /// </summary>
    public partial class Test_VectorizedLinearAlgebra
    {
        /// <summary>
        /// Asserts that the blocked Cholesky factor has a backward error no larger than the scalar loop and an inverse
        /// residual within twice that of substitution with the scalar factor, at each order.
        /// </summary>
        /// <param name="orders">The matrix orders.</param>
        public static void CholeskyBlockedErrorsMatchScalarLoop(params int[] orders)
        {
            foreach (int n in orders)
            {
                var a = Spd(n, 3 * n, 1e-3);
                Assert.IsTrue(ReferenceCholesky(a, 0d, out var reference, out _, out _));
                var decomposition = new CholeskyDecomposition(new Matrix(a));
                var blocked = decomposition.L.Array;
                double referenceBackward = BackwardError(a, reference), blockedBackward = BackwardError(a, blocked);

                var inverse = decomposition.InverseA().Array;
                var referenceInverse = new double[n, n];
                var column = new double[n];
                for (int c = 0; c < n; c++)
                {
                    for (int i = 0; i < n; i++)
                    {
                        double sum = i == c ? 1d : 0d;
                        for (int k = 0; k < i; k++) sum -= reference[i, k] * column[k];
                        column[i] = sum / reference[i, i];
                    }
                    for (int i = n - 1; i >= 0; i--)
                    {
                        double sum = column[i];
                        for (int k = i + 1; k < n; k++) sum -= reference[k, i] * column[k];
                        column[i] = sum / reference[i, i];
                    }
                    for (int i = 0; i < n; i++) referenceInverse[i, c] = column[i];
                }
                double referenceResidual = InverseResidual(a, referenceInverse), blockedResidual = InverseResidual(a, inverse);
                Console.WriteLine($"n = {n}: backward {blockedBackward:E2} vs {referenceBackward:E2}, inverse residual {blockedResidual:E2} vs {referenceResidual:E2}");
                Assert.IsLessThanOrEqualTo(referenceBackward, blockedBackward, $"n = {n}");
                Assert.IsLessThanOrEqualTo(2d * referenceResidual, blockedResidual, $"n = {n}");
            }
        }

        /// <summary>
        /// Random values with special entries sprinkled in.
        /// </summary>
        private static double[,] Values(int rows, int cols, int seed, bool special = true)
        {
            var rng = new Random(seed);
            var a = new double[rows, cols];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    double v = (rng.NextDouble() - 0.5) * Math.Pow(10, rng.Next(-3, 4));
                    if (special)
                    {
                        int r = rng.Next(40);
                        if (r == 0) v = -0.0;
                        else if (r == 1) v = 0.0;
                        else if (r == 2 && rng.Next(4) == 0) v = double.NaN;
                        else if (r == 3 && rng.Next(4) == 0) v = rng.Next(2) == 0 ? double.PositiveInfinity : double.NegativeInfinity;
                    }
                    a[i, j] = v;
                }
            }
            return a;
        }

        /// <summary>
        /// The scalar triple loop: each entry summed from 0 over k in increasing order.
        /// </summary>
        private static double[,] ReferenceMultiply(double[,] a, double[,] b)
        {
            int rows = a.GetLength(0), inner = a.GetLength(1), cols = b.GetLength(1);
            var c = new double[rows, cols];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < cols; j++)
                {
                    double sum = 0.0d;
                    for (int k = 0; k < inner; k++) sum += a[i, k] * b[k, j];
                    c[i, j] = sum;
                }
            return c;
        }

        /// <summary>
        /// A scalar Cholesky loop: row-major, inner products over k from i − 1 down to 0.
        /// </summary>
        private static bool ReferenceCholesky(double[,] source, double tol, out double[,] lower, out int failedRow, out double failedPivot)
        {
            int n = source.GetLength(0);
            var w = (double[,])source.Clone();
            failedRow = -1;
            failedPivot = double.NaN;
            for (int i = 0; i < n; i++)
            {
                for (int j = i; j < n; j++)
                {
                    double sum = w[i, j];
                    for (int k = i - 1; k >= 0; k--) sum -= w[i, k] * w[j, k];
                    if (i == j)
                    {
                        double diagonal = source[i, i];
                        double threshold = diagonal > 0d && !double.IsInfinity(diagonal) ? tol * diagonal : 0d;
                        if (double.IsNaN(sum) || sum <= 0d || sum <= threshold)
                        {
                            failedRow = i;
                            failedPivot = sum;
                            lower = w;
                            return false;
                        }
                        w[i, i] = Math.Sqrt(sum);
                    }
                    else w[j, i] = sum / w[i, i];
                }
            }
            for (int i = 0; i < n; i++) for (int j = 0; j < i; j++) w[j, i] = 0.0d;
            lower = w;
            return true;
        }

        /// <summary>
        /// Symmetric positive definite test matrix XᵀX + shift·I.
        /// </summary>
        private static double[,] Spd(int n, int seed, double shift)
        {
            var x = Values(n + 5, n, seed, special: false);
            var a = ReferenceMultiply(Transpose(x), x);
            for (int i = 0; i < n; i++) a[i, i] += shift;
            return a;
        }

        private static double[,] Transpose(double[,] x)
        {
            var t = new double[x.GetLength(1), x.GetLength(0)];
            for (int i = 0; i < x.GetLength(0); i++) for (int j = 0; j < x.GetLength(1); j++) t[j, i] = x[i, j];
            return t;
        }

        /// <summary>
        /// Returns ‖A − LLᵀ‖_F / ‖A‖_F for the lower triangle of <paramref name="lower"/>.
        /// </summary>
        private static double BackwardError(double[,] a, double[,] lower)
        {
            int n = a.GetLength(0);
            double residual = 0d, norm = 0d;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    double sum = 0d;
                    for (int k = 0; k <= Math.Min(i, j); k++) sum += lower[i, k] * lower[j, k];
                    residual += (a[i, j] - sum) * (a[i, j] - sum);
                    norm += a[i, j] * a[i, j];
                }
            return Math.Sqrt(residual / norm);
        }

        /// <summary>
        /// Returns max |A X − I| over all entries.
        /// </summary>
        private static double InverseResidual(double[,] a, double[,] x)
        {
            int n = a.GetLength(0);
            double max = 0d;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    double sum = 0d;
                    for (int k = 0; k < n; k++) sum += a[i, k] * x[k, j];
                    max = Math.Max(max, Math.Abs(sum - (i == j ? 1d : 0d)));
                }
            return max;
        }
    }
}
