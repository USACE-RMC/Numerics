using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// Unit tests for the internal blocked Cholesky kernels.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// The tests cover a known factor, the backward error below and above the block width and the parallel
    /// threshold, a block addressed by an offset and leading dimension inside a larger array, rejected pivots,
    /// multiple right-hand sides with leading zeros, and the inverse.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_CholeskyKernel
    {
        /// <summary>
        /// A row-major symmetric positive-definite matrix XᵀX / m + shift·I.
        /// </summary>
        private static double[] Spd(int n, int seed, double shift)
        {
            var rng = new Random(seed);
            int m = n + 4;
            var x = new double[m * n];
            for (int i = 0; i < x.Length; i++) x[i] = rng.NextDouble() - 0.5;
            var a = new double[n * n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j <= i; j++)
                {
                    double s = 0;
                    for (int k = 0; k < m; k++) s += x[k * n + i] * x[k * n + j];
                    a[i * n + j] = a[j * n + i] = s / m + (i == j ? shift : 0d);
                }
            return a;
        }

        /// <summary>
        /// The exact factor of a matrix with an integer Cholesky factor.
        /// </summary>
        [TestMethod]
        public void Test_TryFactor_KnownFactor()
        {
            var a = new double[] { 16, 4, 8, 4, 5, -4, 8, -4, 22 };
            Assert.IsTrue(CholeskyKernel.TryFactor(a, 0, 3, 3, 0d, null, out int row, out double pivot));
            Assert.AreEqual(-1, row);
            Assert.IsTrue(double.IsNaN(pivot));
            var expected = new double[] { 4, 4, 8, 1, 2, -4, 2, -3, 3 };
            for (int i = 0; i < 9; i++) Assert.AreEqual(expected[i], a[i], 0d, $"entry {i}");
        }

        /// <summary>
        /// The factor satisfies L Lᵀ = A to the Cholesky backward-error bound below and above the block width and the
        /// parallel threshold, and the solves and the inverse have small residuals.
        /// </summary>
        [TestMethod]
        public void Test_FactorSolveInverse_Residuals()
        {
            foreach (int n in new[] { 1, 2, 5, 31, 128, 129, 300, 400 })
            {
                var a = Spd(n, n, 0.01);
                var l = (double[])a.Clone();
                Assert.IsTrue(CholeskyKernel.TryFactor(l, 0, n, n, CholeskyDecomposition.DefaultRelativeTolerance(n), new double[n], out _, out _));
                double scale = 0;
                for (int i = 0; i < n * n; i++) scale = Math.Max(scale, Math.Abs(a[i]));
                for (int i = 0; i < n; i++)
                    for (int j = 0; j <= i; j++)
                    {
                        double s = 0;
                        for (int k = 0; k <= j; k++) s += l[i * n + k] * l[j * n + k];
                        Assert.AreEqual(a[i * n + j], s, 2 * n * 2.22E-16 * scale, $"n {n} [{i},{j}]");
                    }

                // Two right-hand sides solved together: A x = b through L y = b and Lᵀ x = y.
                var rng = new Random(n);
                var b = new double[2 * n];
                for (int i = 0; i < b.Length; i++) b[i] = rng.NextDouble() - 0.5;
                var x = (double[])b.Clone();
                CholeskyKernel.SolveLowerRows(l, 0, n, n, x, 0, 2, null);
                CholeskyKernel.SolveUpper(l, 0, n, n, x, 0);
                CholeskyKernel.SolveUpper(l, 0, n, n, x, n);
                for (int p = 0; p < 2; p++)
                    for (int i = 0; i < n; i++)
                    {
                        double s = 0;
                        for (int k = 0; k < n; k++) s += a[i * n + k] * x[p * n + k];
                        Assert.AreEqual(b[p * n + i], s, 1E-10, $"solve n {n} rhs {p} row {i}");
                    }

                if (n > 300) continue;
                var w = (double[])l.Clone();
                var t = new double[n * n];
                var diagonal = new double[n];
                CholeskyKernel.Invert(w, n, t, diagonal);
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        double s = 0;
                        for (int k = 0; k < n; k++)
                        {
                            double inv = k == j ? diagonal[k] : k < j ? w[k * n + j] : w[j * n + k];
                            s += a[i * n + k] * inv;
                        }
                        Assert.AreEqual(i == j ? 1d : 0d, s, 1E-9, $"inverse n {n} [{i},{j}]");
                    }
            }
        }

        /// <summary>
        /// A block addressed by an offset and a leading dimension inside a larger array gives the same factor and
        /// solutions bit for bit as the same block stored on its own, and the surrounding storage is untouched.
        /// </summary>
        [TestMethod]
        public void Test_OffsetAndLeadingDimension_MatchCompactStorage()
        {
            foreach (int n in new[] { 3, 17, 130 })
            {
                var a = Spd(n, 7 * n, 0.1);
                int ld = n + 5, rowOffset = 2, colOffset = 3;
                var big = new double[(n + rowOffset + 1) * ld];
                for (int i = 0; i < big.Length; i++) big[i] = -999d;
                int offset = rowOffset * ld + colOffset;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) big[offset + i * ld + j] = a[i * n + j];
                var sentinel = (double[])big.Clone();

                var compact = (double[])a.Clone();
                Assert.IsTrue(CholeskyKernel.TryFactor(compact, 0, n, n, 0d, null, out _, out _));
                Assert.IsTrue(CholeskyKernel.TryFactor(big, offset, ld, n, 0d, null, out _, out _));
                for (int i = 0; i < big.Length; i++)
                {
                    int r = i / ld - rowOffset, c = i % ld - colOffset;
                    bool inLower = r >= 0 && r < n && c >= 0 && c <= r;
                    if (inLower) Assert.AreEqual(BitConverter.DoubleToInt64Bits(compact[r * n + c]), BitConverter.DoubleToInt64Bits(big[i]));
                    else Assert.AreEqual(sentinel[i], big[i], 0d, $"storage {i} was written");
                }

                var rhs = new double[n + 4];
                var rng = new Random(n);
                for (int i = 0; i < rhs.Length; i++) rhs[i] = rng.NextDouble();
                var x1 = new double[n];
                Array.Copy(rhs, 4, x1, 0, n);
                var x2 = (double[])rhs.Clone();
                CholeskyKernel.SolveLowerRows(compact, 0, n, n, x1, 0, 1, null);
                CholeskyKernel.SolveUpper(compact, 0, n, n, x1, 0);
                CholeskyKernel.SolveLowerRows(big, offset, ld, n, x2, 4, 1, null);
                CholeskyKernel.SolveUpper(big, offset, ld, n, x2, 4);
                for (int i = 0; i < 4; i++) Assert.AreEqual(rhs[i], x2[i], 0d);
                for (int i = 0; i < n; i++) Assert.AreEqual(BitConverter.DoubleToInt64Bits(x1[i]), BitConverter.DoubleToInt64Bits(x2[4 + i]));
            }
        }

        /// <summary>
        /// A rejected pivot reports its row and value; the relative test rejects a pivot the absolute test accepts,
        /// and NaN is rejected.
        /// </summary>
        [TestMethod]
        public void Test_TryFactor_RejectedPivots()
        {
            var indefinite = new double[] { 1, 2, 2, 1 };
            Assert.IsFalse(CholeskyKernel.TryFactor(indefinite, 0, 2, 2, 0d, null, out int row, out double pivot));
            Assert.AreEqual(1, row);
            Assert.AreEqual(-3d, pivot, 0d);

            double rho = 0.99999999999999978d;
            var nearSingular = new double[] { 1, rho, rho, 1 };
            Assert.IsFalse(CholeskyKernel.TryFactor((double[])nearSingular.Clone(), 0, 2, 2, CholeskyDecomposition.DefaultRelativeTolerance(2), null, out row, out pivot));
            Assert.AreEqual(1, row);
            Assert.AreEqual(4.440892098500626E-16d, pivot, 0d);
            Assert.IsTrue(CholeskyKernel.TryFactor((double[])nearSingular.Clone(), 0, 2, 2, 0d, null, out _, out _));

            var nan = new double[] { 4, double.NaN, double.NaN, 4 };
            Assert.IsFalse(CholeskyKernel.TryFactor(nan, 0, 2, 2, 0d, null, out row, out pivot));
            Assert.AreEqual(1, row);
            Assert.IsTrue(double.IsNaN(pivot));
        }

        /// <summary>
        /// Right-hand sides with known leading zeros, the form used to build the inverse factor, give the same
        /// solutions as dense right-hand sides.
        /// </summary>
        [TestMethod]
        public void Test_SolveLowerRows_StartIndices()
        {
            int n = 150, m = 20;
            var l = Spd(n, 3, 0.05);
            Assert.IsTrue(CholeskyKernel.TryFactor(l, 0, n, n, 0d, null, out _, out _));
            var rng = new Random(5);
            var x = new double[m * n];
            var start = new int[m];
            for (int p = 0; p < m; p++)
            {
                start[p] = rng.Next(n);
                for (int i = start[p]; i < n; i++) x[p * n + i] = rng.NextDouble() - 0.5;
            }
            var dense = (double[])x.Clone();
            CholeskyKernel.SolveLowerRows(l, 0, n, n, x, 0, m, start);
            CholeskyKernel.SolveLowerRows(l, 0, n, n, dense, 0, m, null);
            for (int p = 0; p < m; p++)
                for (int i = 0; i < n; i++)
                {
                    if (i < start[p]) Assert.AreEqual(0d, x[p * n + i], 0d);
                    Assert.AreEqual(dense[p * n + i], x[p * n + i], 1E-12 * Math.Max(1d, Math.Abs(dense[p * n + i])));
                }
        }
    }
}
