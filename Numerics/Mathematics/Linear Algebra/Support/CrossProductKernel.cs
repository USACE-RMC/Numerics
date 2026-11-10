using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
#if NET8_0_OR_GREATER
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
#endif

namespace Numerics.Mathematics.LinearAlgebra
{
    /// <summary>
    /// Deterministic, cache-blocked and SIMD-accelerated column cross products for column-major data.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// The data are stored column-major in a flat array, so column j of an n-row matrix occupies
    /// elements [j·n, (j + 1)·n). Each row segment is split into row blocks whose number depends only on the
    /// segment length, never on the processor count. Blocks accumulate their partial products in parallel, and the
    /// partial products are added in block order, so the same request gives the same result from run to run on the
    /// same machine. Products in four-column groups are accumulated lane-wise in SIMD registers and the remaining
    /// columns in scalar loops, so a product's last bits can depend on the position of its column in the request. The
    /// four-column groups use 256-bit fused multiply-adds on .NET 8 and later where the hardware supports them (see
    /// <see cref="Dot4x2"/>), and the SIMD width of <see cref="System.Numerics.Vector{T}"/> sets the lane-wise summation
    /// order elsewhere, so results can differ in the last bits between processor architectures and target frameworks.
    /// </para>
    /// </remarks>
    internal static class CrossProductKernel
    {
        /// <summary>
        /// Rows per inner tile, small enough that a tile of the requested columns stays in the core-local cache.
        /// </summary>
        private const int TileRows = 256;

        /// <summary>
        /// Minimum rows per parallel row block.
        /// </summary>
        private const int MinBlockRows = 2048;

        /// <summary>
        /// Maximum number of row blocks per segment.
        /// </summary>
        private const int MaxBlocks = 16;

        /// <summary>
        /// Computes C[l, r] = Σ data[left[l]·n + i]·data[right[r]·n + i] over each row segment.
        /// </summary>
        /// <param name="data">Column-major data with column stride <paramref name="n"/>.</param>
        /// <param name="n">The number of rows (the column stride).</param>
        /// <param name="left">The left column indexes.</param>
        /// <param name="right">The right column indexes.</param>
        /// <param name="segments">Row segment boundaries: segment s covers rows [segments[s], segments[s + 1]).</param>
        /// <param name="results">One output per segment, each of length left.Length·right.Length, stored as
        /// result[r·left.Length + l].</param>
        /// <param name="symmetric">If true, <paramref name="left"/> and <paramref name="right"/> must be the same list;
        /// only l ≤ r + 1 is computed and the rest is mirrored.</param>
        public static void Compute(double[] data, int n, int[] left, int[] right, int[] segments, double[][] results, bool symmetric = false)
        {
            int L = left.Length, R = right.Length, S = segments.Length - 1;
            var blockCount = new int[S];
            var taskStart = new int[S + 1];
            for (int s = 0; s < S; s++)
            {
                int rows = segments[s + 1] - segments[s];
                blockCount[s] = rows <= 0 ? 0 : Math.Max(1, Math.Min(MaxBlocks, rows / MinBlockRows));
                taskStart[s + 1] = taskStart[s] + blockCount[s];
                Array.Clear(results[s], 0, L * R);
            }
            int tasks = taskStart[S];
            if (tasks == 0 || L == 0 || R == 0) return;

            var partial = new double[tasks][];
            Action<int> run = t =>
            {
                int s = 0;
                while (taskStart[s + 1] <= t) s++;
                int b = t - taskStart[s];
                int rows = segments[s + 1] - segments[s];
                int size = rows / blockCount[s];
                int r0 = segments[s] + b * size;
                int r1 = b == blockCount[s] - 1 ? segments[s + 1] : r0 + size;
                var buffer = blockCount[s] == 1 ? results[s] : new double[L * R];
                Accumulate(data, n, left, right, r0, r1, buffer, symmetric);
                partial[t] = buffer;
            };
            if (tasks == 1) run(0);
            else Parallel.For(0, tasks, run);

            for (int s = 0; s < S; s++)
            {
                var result = results[s];
                if (blockCount[s] > 1)
                {
                    for (int t = taskStart[s]; t < taskStart[s + 1]; t++)
                    {
                        var buffer = partial[t];
                        for (int i = 0; i < result.Length; i++) result[i] += buffer[i];
                    }
                }
                if (symmetric)
                {
                    for (int r = 0; r < R; r++)
                        for (int l = r + 1; l < L; l++)
                            result[r * L + l] = result[l * L + r];
                }
            }
        }

        /// <summary>
        /// Computes the cross-product matrix of row-major data shifted by a constant row, Σᵢ (xᵢ − s)(xᵢ − s)ᵀ, and the
        /// column sums of the shifted data, without transposing the data.
        /// </summary>
        /// <param name="x">The n × p data, observations in rows.</param>
        /// <param name="shift">The p values subtracted from every row.</param>
        /// <param name="sums">Output: the p column sums of the shifted data.</param>
        /// <param name="result">Output: the p × p symmetric cross-product matrix, row-major.</param>
        /// <param name="rowStart">The first row of the data to include.</param>
        /// <remarks>
        /// Rows are processed eight at a time as a rank-8 update of the upper triangle, so each vector of the result is
        /// loaded and stored once per eight rows. Row blocks are fixed by n, as in <see cref="Compute"/>, so the result
        /// is reproducible.
        /// </remarks>
        public static void RowMajorGram(double[,] x, double[] shift, double[] sums, double[] result, int rowStart = 0)
        {
            int n = x.GetLength(0), p = x.GetLength(1);
            int w = System.Numerics.Vector<double>.Count;
            int stride = (p + w - 1) / w * w;
            int rows = n - rowStart;
            int blocks = Math.Max(1, Math.Min(MaxBlocks, rows / MinBlockRows));
            var partialGram = new double[blocks][];
            var partialSums = new double[blocks][];
            Action<int> run = b =>
            {
                int size = rows / blocks;
                int r0 = rowStart + b * size, r1 = b == blocks - 1 ? n : r0 + size;
                var g = new double[p * stride];
                var sum = new double[stride];
                var t = new double[8 * stride];
                for (int i0 = r0; i0 < r1; i0 += 8)
                {
                    int m = Math.Min(8, r1 - i0);
                    for (int k = 0; k < 8; k++)
                    {
                        int o = k * stride;
                        if (k < m)
                        {
                            Buffer.BlockCopy(x, (i0 + k) * p * sizeof(double), t, o * sizeof(double), p * sizeof(double));
                            for (int c = 0; c < p; c++)
                            {
                                double v = t[o + c] - shift[c];
                                t[o + c] = v;
                                sum[c] += v;
                            }
                        }
                        else Array.Clear(t, o, stride);
                    }
                    for (int a = 0; a < p; a++)
                    {
                        var v0 = new System.Numerics.Vector<double>(t[a]);
                        var v1 = new System.Numerics.Vector<double>(t[stride + a]);
                        var v2 = new System.Numerics.Vector<double>(t[2 * stride + a]);
                        var v3 = new System.Numerics.Vector<double>(t[3 * stride + a]);
                        var v4 = new System.Numerics.Vector<double>(t[4 * stride + a]);
                        var v5 = new System.Numerics.Vector<double>(t[5 * stride + a]);
                        var v6 = new System.Numerics.Vector<double>(t[6 * stride + a]);
                        var v7 = new System.Numerics.Vector<double>(t[7 * stride + a]);
                        int ga = a * stride;
                        for (int j = a / w * w; j < stride; j += w)
                        {
                            var acc = v0 * Load(t, j) + v1 * Load(t, stride + j) + v2 * Load(t, 2 * stride + j) + v3 * Load(t, 3 * stride + j)
                                + v4 * Load(t, 4 * stride + j) + v5 * Load(t, 5 * stride + j) + v6 * Load(t, 6 * stride + j) + v7 * Load(t, 7 * stride + j);
                            (Load(g, ga + j) + acc).CopyTo(g, ga + j);
                        }
                    }
                }
                partialGram[b] = g;
                partialSums[b] = sum;
            };
            if (blocks == 1) run(0);
            else Parallel.For(0, blocks, run);

            Array.Clear(sums, 0, p);
            Array.Clear(result, 0, p * p);
            for (int b = 0; b < blocks; b++)
            {
                var g = partialGram[b];
                for (int a = 0; a < p; a++)
                    for (int j = a; j < p; j++) result[a * p + j] += g[a * stride + j];
                for (int c = 0; c < p; c++) sums[c] += partialSums[b][c];
            }
            for (int a = 0; a < p; a++)
                for (int j = a + 1; j < p; j++) result[j * p + a] = result[a * p + j];
        }

        /// <summary>
        /// Sequentially accumulates the cross products over one row block, tile by tile.
        /// </summary>
        private static void Accumulate(double[] data, int n, int[] left, int[] right, int rowStart, int rowEnd, double[] result, bool symmetric)
        {
            int L = left.Length, R = right.Length;
            var s = new double[8];
            for (int t0 = rowStart; t0 < rowEnd; t0 += TileRows)
            {
                int len = Math.Min(TileRows, rowEnd - t0);
                int r = 0;
                for (; r + 1 < R; r += 2)
                {
                    int ob0 = right[r] * n + t0, ob1 = right[r + 1] * n + t0;
                    int lEnd = symmetric ? r + 2 : L;
                    int l = 0;
                    for (; l + 3 < lEnd; l += 4)
                    {
                        Dot4x2(data, left[l] * n + t0, left[l + 1] * n + t0, left[l + 2] * n + t0, left[l + 3] * n + t0,
                            data, ob0, ob1, len, s);
                        int r0 = r * L + l, r1 = (r + 1) * L + l;
                        result[r0] += s[0]; result[r0 + 1] += s[2]; result[r0 + 2] += s[4]; result[r0 + 3] += s[6];
                        result[r1] += s[1]; result[r1 + 1] += s[3]; result[r1 + 2] += s[5]; result[r1 + 3] += s[7];
                    }
                    for (; l < lEnd; l++)
                    {
                        int oa = left[l] * n + t0;
                        double s0 = 0d, s1 = 0d;
                        for (int i = 0; i < len; i++)
                        {
                            double v = data[oa + i];
                            s0 += v * data[ob0 + i];
                            s1 += v * data[ob1 + i];
                        }
                        result[r * L + l] += s0;
                        result[(r + 1) * L + l] += s1;
                    }
                }
                if (r < R)
                {
                    int ob0 = right[r] * n + t0;
                    int lEnd = symmetric ? r + 1 : L;
                    for (int l = 0; l < lEnd; l++)
                    {
                        int oa = left[l] * n + t0;
                        double s0 = 0d;
                        for (int i = 0; i < len; i++) s0 += data[oa + i] * data[ob0 + i];
                        result[r * L + l] += s0;
                    }
                }
            }
        }

        /// <summary>
        /// Computes the eight inner products between four segments of <paramref name="a"/> and two segments of
        /// <paramref name="b"/>, all of the same length.
        /// </summary>
        /// <param name="a">The array holding the four left segments.</param>
        /// <param name="a0">The start of the first left segment.</param>
        /// <param name="a1">The start of the second left segment.</param>
        /// <param name="a2">The start of the third left segment.</param>
        /// <param name="a3">The start of the fourth left segment.</param>
        /// <param name="b">The array holding the two right segments.</param>
        /// <param name="b0">The start of the first right segment.</param>
        /// <param name="b1">The start of the second right segment.</param>
        /// <param name="length">The segment length.</param>
        /// <param name="s">Receives the products in the order (a0,b0), (a0,b1), (a1,b0), ..., (a3,b1).</param>
        /// <remarks>
        /// On .NET 8 and later, where 256-bit vectors are hardware accelerated, the products are accumulated lane-wise in
        /// 256-bit registers with fused multiply-adds where supported; otherwise, on every target framework, they are
        /// accumulated lane-wise in <see cref="System.Numerics.Vector{T}"/> registers. The remaining elements are added by
        /// a scalar loop.
        /// </remarks>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        internal static void Dot4x2(double[] a, int a0, int a1, int a2, int a3, double[] b, int b0, int b1, int length, double[] s)
        {
            int k = 0;
            double s00 = 0, s01 = 0, s10 = 0, s11 = 0, s20 = 0, s21 = 0, s30 = 0, s31 = 0;
#if NET8_0_OR_GREATER
            if (Vector256.IsHardwareAccelerated && length >= 4)
            {
                ref double ra = ref MemoryMarshal.GetArrayDataReference(a);
                ref double rb = ref MemoryMarshal.GetArrayDataReference(b);
                var c00 = Vector256<double>.Zero; var c01 = Vector256<double>.Zero;
                var c10 = Vector256<double>.Zero; var c11 = Vector256<double>.Zero;
                var c20 = Vector256<double>.Zero; var c21 = Vector256<double>.Zero;
                var c30 = Vector256<double>.Zero; var c31 = Vector256<double>.Zero;
                for (; k <= length - 4; k += 4)
                {
                    var y0 = Vector256.LoadUnsafe(ref rb, (nuint)(b0 + k));
                    var y1 = Vector256.LoadUnsafe(ref rb, (nuint)(b1 + k));
                    var x = Vector256.LoadUnsafe(ref ra, (nuint)(a0 + k));
                    c00 = MultiplyAdd(x, y0, c00); c01 = MultiplyAdd(x, y1, c01);
                    x = Vector256.LoadUnsafe(ref ra, (nuint)(a1 + k));
                    c10 = MultiplyAdd(x, y0, c10); c11 = MultiplyAdd(x, y1, c11);
                    x = Vector256.LoadUnsafe(ref ra, (nuint)(a2 + k));
                    c20 = MultiplyAdd(x, y0, c20); c21 = MultiplyAdd(x, y1, c21);
                    x = Vector256.LoadUnsafe(ref ra, (nuint)(a3 + k));
                    c30 = MultiplyAdd(x, y0, c30); c31 = MultiplyAdd(x, y1, c31);
                }
                s00 = Vector256.Sum(c00); s01 = Vector256.Sum(c01);
                s10 = Vector256.Sum(c10); s11 = Vector256.Sum(c11);
                s20 = Vector256.Sum(c20); s21 = Vector256.Sum(c21);
                s30 = Vector256.Sum(c30); s31 = Vector256.Sum(c31);
            }
            else
#endif
            {
                int w = System.Numerics.Vector<double>.Count;
                if (length >= w)
                {
                    var c00 = System.Numerics.Vector<double>.Zero; var c01 = c00;
                    var c10 = c00; var c11 = c00;
                    var c20 = c00; var c21 = c00;
                    var c30 = c00; var c31 = c00;
                    for (; k + w <= length; k += w)
                    {
                        var y0 = new System.Numerics.Vector<double>(b, b0 + k);
                        var y1 = new System.Numerics.Vector<double>(b, b1 + k);
                        var x = new System.Numerics.Vector<double>(a, a0 + k);
                        c00 += x * y0; c01 += x * y1;
                        x = new System.Numerics.Vector<double>(a, a1 + k);
                        c10 += x * y0; c11 += x * y1;
                        x = new System.Numerics.Vector<double>(a, a2 + k);
                        c20 += x * y0; c21 += x * y1;
                        x = new System.Numerics.Vector<double>(a, a3 + k);
                        c30 += x * y0; c31 += x * y1;
                    }
                    var one = System.Numerics.Vector<double>.One;
                    s00 = System.Numerics.Vector.Dot(c00, one); s01 = System.Numerics.Vector.Dot(c01, one);
                    s10 = System.Numerics.Vector.Dot(c10, one); s11 = System.Numerics.Vector.Dot(c11, one);
                    s20 = System.Numerics.Vector.Dot(c20, one); s21 = System.Numerics.Vector.Dot(c21, one);
                    s30 = System.Numerics.Vector.Dot(c30, one); s31 = System.Numerics.Vector.Dot(c31, one);
                }
            }
            for (; k < length; k++)
            {
                double y0 = b[b0 + k], y1 = b[b1 + k];
                double x0 = a[a0 + k], x1 = a[a1 + k], x2 = a[a2 + k], x3 = a[a3 + k];
                s00 += x0 * y0; s01 += x0 * y1;
                s10 += x1 * y0; s11 += x1 * y1;
                s20 += x2 * y0; s21 += x2 * y1;
                s30 += x3 * y0; s31 += x3 * y1;
            }
            s[0] = s00; s[1] = s01; s[2] = s10; s[3] = s11;
            s[4] = s20; s[5] = s21; s[6] = s30; s[7] = s31;
        }

#if NET8_0_OR_GREATER
        /// <summary>Returns <c>x * y + c</c>, fused when the hardware supports it.</summary>
        /// <param name="x">The first factor.</param>
        /// <param name="y">The second factor.</param>
        /// <param name="c">The addend.</param>
        /// <returns>The product plus the addend.</returns>
        private static Vector256<double> MultiplyAdd(Vector256<double> x, Vector256<double> y, Vector256<double> c)
        {
            return Fma.IsSupported ? Fma.MultiplyAdd(x, y, c) : c + x * y;
        }
#endif

        /// <summary>
        /// Loads a SIMD vector from <paramref name="a"/> starting at index <paramref name="i"/>; the caller guarantees the range.
        /// </summary>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        private static System.Numerics.Vector<double> Load(double[] a, int i)
        {
#if NET
            return System.Runtime.CompilerServices.Unsafe.ReadUnaligned<System.Numerics.Vector<double>>(
                ref System.Runtime.CompilerServices.Unsafe.As<double, byte>(
                    ref System.Runtime.CompilerServices.Unsafe.Add(ref System.Runtime.InteropServices.MemoryMarshal.GetArrayDataReference(a), i)));
#else
            return new System.Numerics.Vector<double>(a, i);
#endif
        }

        /// <summary>
        /// Dot product of two contiguous ranges, accumulated lane-wise in one SIMD register and then by a scalar tail.
        /// </summary>
        /// <param name="x">The first array.</param>
        /// <param name="ox">Offset into <paramref name="x"/>.</param>
        /// <param name="y">The second array.</param>
        /// <param name="oy">Offset into <paramref name="y"/>.</param>
        /// <param name="len">The range length.</param>
        public static double Dot(double[] x, int ox, double[] y, int oy, int len)
        {
            int w = System.Numerics.Vector<double>.Count;
            if (len < w)
            {
                // No full vector: the scalar tail alone, which is the same sum as below.
                double t = 0d;
                for (int k = 0; k < len; k++) t += x[ox + k] * y[oy + k];
                return t;
            }
            var a = System.Numerics.Vector<double>.Zero;
            int i = 0;
            for (; i + w <= len; i += w)
                a += Load(x, ox + i) * Load(y, oy + i);
            double s = System.Numerics.Vector.Dot(a, System.Numerics.Vector<double>.One);
            for (; i < len; i++) s += x[ox + i] * y[oy + i];
            return s;
        }

        /// <summary>
        /// Computes y[oy + i] += a·x[ox + i] for i in [0, len).
        /// </summary>
        /// <param name="a">The scalar multiplier.</param>
        /// <param name="x">The source array.</param>
        /// <param name="ox">Offset into <paramref name="x"/>.</param>
        /// <param name="y">The destination array.</param>
        /// <param name="oy">Offset into <paramref name="y"/>.</param>
        /// <param name="len">The range length.</param>
        public static void Axpy(double a, double[] x, int ox, double[] y, int oy, int len)
        {
            int w = System.Numerics.Vector<double>.Count;
            var va = new System.Numerics.Vector<double>(a);
            int i = 0;
            for (; i + w <= len; i += w)
            {
                var r = Load(y, oy + i) + va * Load(x, ox + i);
                r.CopyTo(y, oy + i);
            }
            for (; i < len; i++) y[oy + i] += a * x[ox + i];
        }
    }
}
