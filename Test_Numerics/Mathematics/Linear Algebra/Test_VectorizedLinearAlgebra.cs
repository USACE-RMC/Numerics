using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// Bit-equality tests for the vectorized matrix products, copies and LU elimination, and accuracy tests for the
    /// SIMD Gram matrix and the blocked Cholesky factorization.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// Each test compares the library result with a scalar reference loop that performs the same operations in the
    /// same order, on sizes below and above the vectorization thresholds and on data containing −0.0, NaN, infinities,
    /// exact zeros and mixed magnitudes. Results must agree bit for bit; NaN results only need to be NaN in both, since
    /// subtracting a product and adding its negation can differ in the sign bit of a NaN.
    /// </para>
    /// </remarks>
    [TestClass]
    public partial class Test_VectorizedLinearAlgebra
    {
        /// <summary>
        /// Asserts bitwise equality, with any NaN matching any NaN.
        /// </summary>
        internal static void AssertSame(double expected, double actual, string where)
        {
            if (double.IsNaN(expected) && double.IsNaN(actual)) return;
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(actual), $"{where}: expected {expected:R}, actual {actual:R}");
        }

        private static void AssertSame(double[,] expected, double[,] actual, string where)
        {
            Assert.AreEqual(expected.GetLength(0), actual.GetLength(0));
            Assert.AreEqual(expected.GetLength(1), actual.GetLength(1));
            for (int i = 0; i < expected.GetLength(0); i++)
                for (int j = 0; j < expected.GetLength(1); j++)
                    AssertSame(expected[i, j], actual[i, j], $"{where}[{i},{j}]");
        }

        /// <summary>
        /// Matrix products match the scalar triple loop bit for bit, below and above the vectorization width.
        /// </summary>
        [TestMethod]
        public void Test_Multiply_BitIdentical()
        {
            foreach (var (rows, inner, cols) in new[] { (5, 7, 3), (6, 9, 8), (13, 31, 17), (40, 70, 33), (3, 1, 12) })
            {
                var a = Values(rows, inner, rows * 31 + cols);
                var b = Values(inner, cols, inner * 17 + rows);
                var expected = ReferenceMultiply(a, b);
                AssertSame(expected, new Matrix(a).Multiply(new Matrix(b)).Array, $"Multiply {rows}x{inner}x{cols}");
                AssertSame(expected, (new Matrix(a) * new Matrix(b)).Array, "operator *");
            }
        }

        /// <summary>
        /// The transpose product matches the product with an explicit transpose bit for bit; the Gram matrix agrees with
        /// it to rounding and is exactly symmetric.
        /// </summary>
        [TestMethod]
        public void Test_TransposeMultiply_BitIdentical()
        {
            foreach (var (inner, rows, cols) in new[] { (9, 4, 5), (50, 12, 20), (200, 33, 33), (7, 1, 9) })
            {
                var a = Values(inner, rows, inner + rows);
                var b = Values(inner, cols, inner + cols + 1);
                var A = new Matrix(a);
                var at = new double[rows, inner];
                for (int i = 0; i < inner; i++) for (int j = 0; j < rows; j++) at[j, i] = a[i, j];
                AssertSame(ReferenceMultiply(at, b), A.TransposeMultiply(new Matrix(b)).Array, "TransposeMultiply");
                // The Gram matrix uses the SIMD kernel: within rounding of the row-by-row sum, and exactly symmetric.
                var gram = A.TransposeMultiply(A);
                var exact = ReferenceMultiply(at, a);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < rows; j++)
                    {
                        AssertSame(gram[i, j], gram[j, i], "Gram symmetry");
                        double e = exact[i, j], g = gram[i, j];
                        if (double.IsNaN(e) || double.IsInfinity(e))
                        {
                            Assert.AreEqual(double.IsNaN(e), double.IsNaN(g), "Gram NaN");
                            if (!double.IsNaN(e)) Assert.AreEqual(e, g, "Gram infinity");
                            continue;
                        }
                        double scale = 0;
                        for (int k = 0; k < inner; k++) scale += Math.Abs(a[k, i] * a[k, j]);
                        Assert.AreEqual(e, g, 1E-14 * inner * scale, "Gram");
                    }
            }
            Assert.ThrowsExactly<ArgumentException>(() => new Matrix(3, 2).TransposeMultiply(new Matrix(4, 2)));
            Assert.ThrowsExactly<ArgumentNullException>(() => new Matrix(3, 2).TransposeMultiply(null!));
        }

        /// <summary>
        /// Copies and transposes reproduce every bit, including −0.0 and NaN.
        /// </summary>
        [TestMethod]
        public void Test_Copies_BitIdentical()
        {
            var a = Values(37, 45, 5);
            var m = new Matrix(a);
            AssertSame(a, m.Array, "ctor");
            Assert.AreNotSame(a, m.Array);
            AssertSame(a, m.Clone().Array, "Clone");
            var t = m.Transpose();
            var ts = Matrix.Transpose(m);
            for (int i = 0; i < 37; i++)
                for (int j = 0; j < 45; j++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(a[i, j]), BitConverter.DoubleToInt64Bits(t[j, i]));
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(a[i, j]), BitConverter.DoubleToInt64Bits(ts[j, i]));
                }
            var v = new double[] { -0.0, 1, double.NaN, 3 };
            var column = new Matrix(v);
            for (int i = 0; i < v.Length; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(v[i]), BitConverter.DoubleToInt64Bits(column[i, 0]));
        }

        /// <summary>
        /// Checks the incomplete factor of a rejected matrix against the scalar loop: rows up to the rejected one,
        /// in the columns before it, agree within <paramref name="tolerance"/> (bitwise when it is zero), NaN matching
        /// NaN, and the strict upper triangle is zero. Entries the scalar loop computes below the rejected row are not
        /// compared, because the blocked kernel has not reached them.
        /// </summary>
        private static void AssertIncompleteFactor(double[,] expected, double[,] actual, int failedRow, double tolerance, string where)
        {
            int n = expected.GetLength(0);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (j > i) Assert.AreEqual(0d, actual[i, j], 0d, $"{where} upper [{i},{j}]");
                    else if (i > failedRow || j >= failedRow) continue;
                    else if (tolerance == 0d || double.IsNaN(expected[i, j])) AssertSame(expected[i, j], actual[i, j], $"{where} [{i},{j}]");
                    else Assert.AreEqual(expected[i, j], actual[i, j], tolerance, $"{where} [{i},{j}]");
                }
            }
        }

        /// <summary>
        /// From order 16 the blocked factor has a backward error ‖A − LLᵀ‖_F / ‖A‖_F no larger than the scalar loop
        /// (measured 1.05 to 1.7 times smaller), and the blocked inverse a residual max |A A⁻¹ − I| within twice that of
        /// column-by-column substitution with the scalar factor (measured 0.8 to 1.3 times it), on XᵀX + 10⁻³ I test
        /// matrices of order 16 to 400. This test checks orders 16, 64 and 200; the sweep that adds order 400 runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep Cholesky</c>.
        /// </summary>
        [TestMethod]
        public void Test_Cholesky_BlockedErrorsMatchScalarLoop()
        {
            CholeskyBlockedErrorsMatchScalarLoop(16, 64, 200);
        }

        /// <summary>
        /// The blocked Cholesky factor agrees with the scalar loop to rounding error and has a backward error within the
        /// Cholesky bound, and a rejected pivot is reported at the same row with the same value to rounding error, below
        /// and above the block width, including a NaN entry and an indefinite matrix.
        /// </summary>
        [TestMethod]
        public void Test_Cholesky_MatchesScalarLoop()
        {
            foreach (int n in new[] { 10, 47, 48, 64, 97, 130 })
            {
                var a = Spd(n, n, 1e-3);
                Assert.IsTrue(ReferenceCholesky(a, 0d, out var expected, out _, out _));
                Assert.IsTrue(CholeskyDecomposition.TryFactorizeFlat(new Matrix(a), 0d, out var factor, out int row, out double pivot));
                Assert.AreEqual(-1, row);
                Assert.IsTrue(double.IsNaN(pivot));
                double scaleL = 0, scaleA = 0;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        scaleL = Math.Max(scaleL, Math.Abs(expected[i, j]));
                        scaleA = Math.Max(scaleA, Math.Abs(a[i, j]));
                    }
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        Assert.AreEqual(expected[i, j], factor[i * n + j], 1E-12 * scaleL, $"factor {n} [{i},{j}]");
                        if (j > i) Assert.AreEqual(0d, factor[i * n + j]);
                        if (j > i) continue;
                        double product = 0;
                        for (int k = 0; k <= j; k++) product += factor[i * n + k] * factor[j * n + k];
                        Assert.AreEqual(a[i, j], product, 4 * n * 2.22E-16 * scaleA, $"residual {n} [{i},{j}]");
                    }
                }
                var ctor = new CholeskyDecomposition(new Matrix(a)).L.Array;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) AssertSame(factor[i * n + j], ctor[i, j], $"Cholesky ctor {n}");

                // Indefinite: a negative diagonal entry part way through.
                var bad = (double[,])a.Clone();
                bad[n / 2, n / 2] = -bad[n / 2, n / 2];
                Assert.IsFalse(ReferenceCholesky(bad, 0d, out var expectedBad, out int expectedRow, out double expectedPivot));
                Assert.IsFalse(CholeskyDecomposition.TryFactorize(new Matrix(bad), 0d, out var lowerBad, out int badRow, out double badPivot));
                Assert.AreEqual(expectedRow, badRow);
                Assert.AreEqual(expectedPivot, badPivot, 1E-12 * Math.Abs(expectedPivot), "pivot");
                AssertIncompleteFactor(expectedBad, lowerBad.Array, badRow, n < CholeskyDecomposition.BlockedDimension ? 0d : 1E-12 * scaleL, $"Cholesky failure {n}");

                // A NaN below the first row.
                var nan = (double[,])a.Clone();
                nan[n - 1, 1] = double.NaN;
                nan[1, n - 1] = double.NaN;
                Assert.AreEqual(ReferenceCholesky(nan, 0d, out var expectedNan, out int nanRow, out _),
                    CholeskyDecomposition.TryFactorize(new Matrix(nan), 0d, out var lowerNan, out int nanRow2, out _));
                Assert.AreEqual(nanRow, nanRow2);
                AssertIncompleteFactor(expectedNan, lowerNan.Array, nanRow2, n < CholeskyDecomposition.BlockedDimension ? 0d : 1E-12 * scaleL, $"Cholesky NaN {n}");
            }
        }

        /// <summary>
        /// The scalar elimination of <see cref="LUDecomposition"/>.
        /// </summary>
        private static double[,] ReferenceLU(double[,] source, out int[] index)
        {
            int n = source.GetLength(0);
            var lu = (double[,])source.Clone();
            index = new int[n];
            var vv = new double[n];
            for (int i = 0; i < n; i++)
            {
                double big = 0.0d;
                for (int j = 0; j < n; j++)
                {
                    double temp = Math.Abs(lu[i, j]);
                    if (temp > big) big = temp;
                }
                vv[i] = 1.0d / big;
            }
            for (int k = 0; k < n; k++)
            {
                double big = 0.0d;
                int imax = k;
                for (int i = k; i < n; i++)
                {
                    double temp = vv[i] * Math.Abs(lu[i, k]);
                    if (temp > big) { big = temp; imax = i; }
                }
                if (k != imax)
                {
                    for (int j = 0; j < n; j++) { double temp = lu[imax, j]; lu[imax, j] = lu[k, j]; lu[k, j] = temp; }
                    vv[imax] = vv[k];
                }
                index[k] = imax;
                if (lu[k, k] == 0.0d) lu[k, k] = 1.0E-40d;
                for (int i = k + 1; i < n; i++)
                {
                    lu[i, k] /= lu[k, k];
                    for (int j = k + 1; j < n; j++) lu[i, j] -= lu[i, k] * lu[k, j];
                }
            }
            return lu;
        }

        /// <summary>
        /// The LU factors and determinant match the scalar elimination bit for bit below and above the vectorization
        /// order, including a matrix with −0.0 entries and an exactly singular column that takes the tiny-pivot path.
        /// </summary>
        [TestMethod]
        public void Test_LU_BitIdentical()
        {
            foreach (int n in new[] { 8, 31, 32, 57 })
            {
                var a = Values(n, n, 100 + n, special: false);
                a[1, 2] = -0.0;
                a[n - 1, 0] = -0.0;
                var expected = ReferenceLU(a, out var index);
                var lu = new LUDecomposition(new Matrix(a));
                AssertSame(expected, lu.LU.Array, $"LU {n}");

                // Two equal rows make the last pivot exactly zero.
                var singular = (double[,])a.Clone();
                for (int j = 0; j < n; j++) singular[n - 1, j] = singular[0, j];
                AssertSame(ReferenceLU(singular, out _), new LUDecomposition(new Matrix(singular)).LU.Array, $"LU singular {n}");
            }
        }
    }
}
