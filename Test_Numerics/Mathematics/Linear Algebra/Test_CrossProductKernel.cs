using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// Unit tests for the internal SIMD cross-product kernel.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// The kernel uses <see cref="System.Numerics.Vector{T}"/> loads through unchecked references on .NET 8 and later
    /// and through the array constructor on .NET Framework, so these tests cover both paths across the target
    /// frameworks. Reference values are plain scalar loops; the tolerance of 1e-13 relative allows for the lane-wise
    /// summation order.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_CrossProductKernel
    {
        /// <summary>
        /// Column-major random data with lengths that leave SIMD and tile remainders.
        /// </summary>
        private static double[] Data(int n, int p, int seed)
        {
            var rng = new Random(seed);
            var z = new double[n * p];
            for (int i = 0; i < z.Length; i++) z[i] = rng.NextDouble() * 2 - 1 + 0.1 * (i % p);
            return z;
        }

        private static double Naive(double[] z, int n, int a, int b, int r0, int r1)
        {
            double s = 0;
            for (int i = r0; i < r1; i++) s += z[a * n + i] * z[b * n + i];
            return s;
        }

        /// <summary>
        /// General left and right column sets over several row segments, including a multi-block segment, match scalar
        /// sums and are reproducible from call to call.
        /// </summary>
        [TestMethod]
        public void Test_Compute_Segments()
        {
            int n = 9001, p = 11;
            var z = Data(n, p, 1);
            var left = new[] { 0, 3, 4, 5, 7, 9, 10 };
            var right = new[] { 2, 3, 8 };
            var segments = new[] { 0, 13, 4500, 9001 };
            var results = new[] { new double[21], new double[21], new double[21] };
            CrossProductKernel.Compute(z, n, left, right, segments, results);
            for (int s = 0; s < 3; s++)
                for (int r = 0; r < right.Length; r++)
                    for (int l = 0; l < left.Length; l++)
                    {
                        double expected = Naive(z, n, left[l], right[r], segments[s], segments[s + 1]);
                        Assert.AreEqual(expected, results[s][r * left.Length + l], 1E-13 * Math.Max(1, Math.Abs(expected)) * (segments[s + 1] - segments[s]));
                    }
            var again = new[] { new double[21], new double[21], new double[21] };
            CrossProductKernel.Compute(z, n, left, right, segments, again);
            for (int s = 0; s < 3; s++) CollectionAssert.AreEqual(results[s], again[s]);
        }

        /// <summary>
        /// The symmetric mode mirrors the upper triangle into an exactly symmetric matrix that matches scalar sums.
        /// </summary>
        [TestMethod]
        public void Test_Compute_Symmetric()
        {
            int n = 5003, p = 9;
            var z = Data(n, p, 2);
            var all = new int[p];
            for (int j = 0; j < p; j++) all[j] = j;
            var g = new double[p * p];
            CrossProductKernel.Compute(z, n, all, all, new[] { 0, n }, new[] { g }, symmetric: true);
            for (int a = 0; a < p; a++)
                for (int b = 0; b < p; b++)
                {
                    Assert.AreEqual(g[a * p + b], g[b * p + a], 0.0);
                    double expected = Naive(z, n, a, b, 0, n);
                    Assert.AreEqual(expected, g[b * p + a], 1E-13 * n * Math.Max(1, Math.Abs(expected)));
                }
        }

        /// <summary>
        /// The row-major Gram of shifted data and its column sums match scalar loops, and the result is exactly
        /// symmetric.
        /// </summary>
        [TestMethod]
        public void Test_RowMajorGram()
        {
            int n = 4099, p = 7;
            var rng = new Random(3);
            var x = new double[n, p];
            for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) x[i, j] = 100 + j + rng.NextDouble();
            var shift = new double[p];
            for (int j = 0; j < p; j++) shift[j] = x[0, j];
            var sums = new double[p];
            var g = new double[p * p];
            CrossProductKernel.RowMajorGram(x, shift, sums, g);
            for (int a = 0; a < p; a++)
            {
                double s = 0;
                for (int i = 0; i < n; i++) s += x[i, a] - shift[a];
                Assert.AreEqual(s, sums[a], 1E-12 * Math.Max(1, Math.Abs(s)));
                for (int b = 0; b < p; b++)
                {
                    double e = 0;
                    for (int i = 0; i < n; i++) e += (x[i, a] - shift[a]) * (x[i, b] - shift[b]);
                    Assert.AreEqual(e, g[a * p + b], 1E-12 * Math.Max(1, Math.Abs(e)));
                    Assert.AreEqual(g[a * p + b], g[b * p + a], 0.0);
                }
            }
        }

        /// <summary>
        /// Dot and Axpy match scalar loops for lengths around the SIMD width, and propagate NaN.
        /// </summary>
        [TestMethod]
        public void Test_Dot_Axpy()
        {
            var rng = new Random(4);
            for (int len = 0; len <= 19; len++)
            {
                var x = new double[len + 3];
                var y = new double[len + 5];
                for (int i = 0; i < x.Length; i++) x[i] = rng.NextDouble() - 0.5;
                for (int i = 0; i < y.Length; i++) y[i] = rng.NextDouble() - 0.5;
                double expected = 0;
                for (int i = 0; i < len; i++) expected += x[3 + i] * y[2 + i];
                Assert.AreEqual(expected, CrossProductKernel.Dot(x, 3, y, 2, len), 1E-14);

                var copy = (double[])y.Clone();
                CrossProductKernel.Axpy(-1.5, x, 1, y, 4, len);
                for (int i = 0; i < y.Length; i++)
                {
                    double e = i >= 4 && i < 4 + len ? copy[i] + -1.5 * x[i - 3] : copy[i];
                    Assert.AreEqual(e, y[i], 0.0);
                }
            }
            var nan = new double[] { 1, 2, double.NaN, 4, 5, 6, 7, 8, 9 };
            Assert.IsTrue(double.IsNaN(CrossProductKernel.Dot(nan, 0, nan, 0, nan.Length)));
        }
    

        /// <summary>
        /// Dot4x2 matches the eight scalar inner products for lengths around the SIMD width and with unaligned offsets,
        /// is exact on small integers, and propagates NaN.
        /// </summary>
        [TestMethod]
        public void Test_Dot4x2()
        {
            var rng = new Random(11);
            var s = new double[8];
            for (int len = 0; len <= 37; len++)
            {
                var a = new double[4 * len + 9];
                var b = new double[2 * len + 7];
                for (int i = 0; i < a.Length; i++) a[i] = rng.NextDouble() - 0.5;
                for (int i = 0; i < b.Length; i++) b[i] = rng.NextDouble() - 0.5;
                int[] ao = { 1, len + 3, 2 * len + 5, 3 * len + 9 - 1 };
                int[] bo = { 2, len + 7 - 1 };
                CrossProductKernel.Dot4x2(a, ao[0], ao[1], ao[2], ao[3], b, bo[0], bo[1], len, s);
                for (int r = 0; r < 4; r++)
                {
                    for (int c = 0; c < 2; c++)
                    {
                        double expected = 0, scale = 0;
                        for (int i = 0; i < len; i++)
                        {
                            expected += a[ao[r] + i] * b[bo[c] + i];
                            scale += Math.Abs(a[ao[r] + i] * b[bo[c] + i]);
                        }
                        Assert.AreEqual(expected, s[2 * r + c], 1E-15 * Math.Max(1d, scale), $"len {len} ({r},{c})");
                    }
                }
            }

            var ints = new double[32];
            for (int i = 0; i < ints.Length; i++) ints[i] = i % 7 - 3;
            CrossProductKernel.Dot4x2(ints, 0, 4, 8, 12, ints, 16, 20, 11, s);
            for (int r = 0; r < 4; r++)
            {
                for (int c = 0; c < 2; c++)
                {
                    double expected = 0;
                    for (int i = 0; i < 11; i++) expected += ints[4 * r + i] * ints[16 + 4 * c + i];
                    Assert.AreEqual(expected, s[2 * r + c], 0d);
                }
            }

            ints[5] = double.NaN;
            CrossProductKernel.Dot4x2(ints, 0, 4, 8, 12, ints, 16, 20, 8, s);
            Assert.IsTrue(double.IsNaN(s[0]) && double.IsNaN(s[2]));
        }
    }
}
