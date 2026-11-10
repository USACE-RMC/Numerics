using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
#if NET
using System.Numerics;
using System.Runtime.InteropServices;
#endif
using SimdVector = System.Numerics.Vector<double>;

namespace Numerics.Mathematics.LinearAlgebra
{
    /// <summary>
    /// Element-wise SIMD kernels for the Householder tridiagonalization (tred2) and the implicit QL iterations (tql2) of
    /// <see cref="EigenValueDecomposition"/>, on flat row-major arrays.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// Every kernel is element-wise: each SIMD lane carries one output element through the same sequence of IEEE
    /// operations, in the same order, as the scalar loop it replaces. There are no horizontal sums and no fused
    /// multiply-adds, so the results do not depend on the SIMD width or the target framework. The rotation kernel works on a
    /// range of columns and gives every column its operations in order, so splitting the columns among parallel tasks does
    /// not change any bit either.
    /// </para>
    /// </remarks>
    internal static class TridiagonalKernel
    {
        /// <summary>
        /// The smallest number of columns given to one parallel task.
        /// </summary>
        internal const int ColumnsPerTask = 64;

        /// <summary>
        /// The largest number of parallel tasks.
        /// </summary>
        internal const int MaxTasks = 16;

        /// <summary>
        /// Runs <paramref name="body"/> on [0, <paramref name="columns"/>), split into blocks of whole 16-column groups that
        /// run in parallel when <paramref name="parallel"/> is true and there are at least two blocks.
        /// </summary>
        /// <param name="columns">The number of columns.</param>
        /// <param name="parallel">Whether the blocks may run in parallel.</param>
        /// <param name="body">The work on the column range [start, end).</param>
        /// <remarks>
        /// There are at most <see cref="MaxTasks"/> blocks, one per processor at most, each of at least
        /// <see cref="ColumnsPerTask"/> columns.
        /// </remarks>
        internal static void ForColumnBlocks(int columns, bool parallel, Action<int, int> body)
        {
            int tasks = parallel ? Math.Min(Math.Min(MaxTasks, Environment.ProcessorCount), columns / ColumnsPerTask) : 1;
            if (tasks < 2)
            {
                body(0, columns);
                return;
            }
            int groups = (columns + 15) / 16;
            Parallel.For(0, tasks, t =>
            {
                int start = Math.Min(columns, (int)((long)groups * t / tasks) * 16);
                int end = Math.Min(columns, (int)((long)groups * (t + 1) / tasks) * 16);
                if (end > start) body(start, end);
            });
        }
        /// <summary>
        /// Forms <c>y[c] = Σ_m coef[m]·M[m·ld + c]</c> for c in [0, len), each sum accumulated from zero over m in
        /// increasing order.
        /// </summary>
        /// <param name="M">The row-major matrix.</param>
        /// <param name="ld">The leading dimension of <paramref name="M"/>.</param>
        /// <param name="coef">The coefficients, one per row.</param>
        /// <param name="y">Output: the len sums.</param>
        /// <param name="len">The number of rows and of columns.</param>
        /// <remarks>
        /// The product is formed as a sequence of row axpys, so every lane sums its column in row order. Four vectors of
        /// columns are held in registers while the rows are streamed.
        /// </remarks>
#if NET
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        internal static void ColumnSums(double[] M, int ld, double[] coef, double[] y, int len)
        {
            int w = SimdVector.Count;
            int c = 0;
#if NET
            nuint w1 = (nuint)w, stride = (nuint)ld;
            ref double rm = ref MemoryMarshal.GetArrayDataReference(M);
            ref double ry = ref MemoryMarshal.GetArrayDataReference(y);
            for (; c + 4 * w <= len; c += 4 * w)
            {
                SimdVector a0 = SimdVector.Zero, a1 = a0, a2 = a0, a3 = a0;
                nuint o = (nuint)c;
                for (int m = 0; m < len; m++, o += stride)
                {
                    var u = new SimdVector(coef[m]);
                    a0 += System.Numerics.Vector.LoadUnsafe(ref rm, o) * u;
                    a1 += System.Numerics.Vector.LoadUnsafe(ref rm, o + w1) * u;
                    a2 += System.Numerics.Vector.LoadUnsafe(ref rm, o + 2 * w1) * u;
                    a3 += System.Numerics.Vector.LoadUnsafe(ref rm, o + 3 * w1) * u;
                }
                nuint q = (nuint)c;
                a0.StoreUnsafe(ref ry, q);
                a1.StoreUnsafe(ref ry, q + w1);
                a2.StoreUnsafe(ref ry, q + 2 * w1);
                a3.StoreUnsafe(ref ry, q + 3 * w1);
            }
#endif
            for (; c + w <= len; c += w)
            {
                SimdVector a0 = SimdVector.Zero;
                for (int m = 0, o = c; m < len; m++, o += ld) a0 += new SimdVector(M, o) * new SimdVector(coef[m]);
                a0.CopyTo(y, c);
            }
            for (; c < len; c++)
            {
                double s = 0d;
                for (int m = 0, o = c; m < len; m++, o += ld) s += M[o] * coef[m];
                y[c] = s;
            }
        }

        /// <summary>
        /// Applies the symmetric rank-two update of tred2 to one row segment: <c>V[ov + k] −= a·x[k] + b·y[oy + k]</c>.
        /// </summary>
        /// <param name="V">The working matrix, updated in place.</param>
        /// <param name="ov">The offset of the row segment in <paramref name="V"/>.</param>
        /// <param name="a">The multiplier of <paramref name="x"/>.</param>
        /// <param name="x">The first update vector, from index zero.</param>
        /// <param name="b">The multiplier of <paramref name="y"/>.</param>
        /// <param name="y">The second update vector.</param>
        /// <param name="oy">The offset into <paramref name="y"/>.</param>
        /// <param name="len">The number of entries.</param>
#if NET
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        internal static void SubtractRankTwo(double[] V, int ov, double a, double[] x, double b, double[] y, int oy, int len)
        {
            int w = SimdVector.Count;
            var va = new SimdVector(a);
            var vb = new SimdVector(b);
            int k = 0;
#if NET
            if (len >= 4 * w)
            {
                ref double rv = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(V), ov);
                ref double rx = ref MemoryMarshal.GetArrayDataReference(x);
                ref double ry = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(y), oy);
                nuint i = 0, count = (nuint)len, w1 = (nuint)w, w4 = 4 * w1;
                for (; i + w4 <= count; i += w4)
                {
                    var v0 = System.Numerics.Vector.LoadUnsafe(ref rv, i) - (va * System.Numerics.Vector.LoadUnsafe(ref rx, i) + vb * System.Numerics.Vector.LoadUnsafe(ref ry, i));
                    var v1 = System.Numerics.Vector.LoadUnsafe(ref rv, i + w1) - (va * System.Numerics.Vector.LoadUnsafe(ref rx, i + w1) + vb * System.Numerics.Vector.LoadUnsafe(ref ry, i + w1));
                    var v2 = System.Numerics.Vector.LoadUnsafe(ref rv, i + 2 * w1) - (va * System.Numerics.Vector.LoadUnsafe(ref rx, i + 2 * w1) + vb * System.Numerics.Vector.LoadUnsafe(ref ry, i + 2 * w1));
                    var v3 = System.Numerics.Vector.LoadUnsafe(ref rv, i + 3 * w1) - (va * System.Numerics.Vector.LoadUnsafe(ref rx, i + 3 * w1) + vb * System.Numerics.Vector.LoadUnsafe(ref ry, i + 3 * w1));
                    v0.StoreUnsafe(ref rv, i);
                    v1.StoreUnsafe(ref rv, i + w1);
                    v2.StoreUnsafe(ref rv, i + 2 * w1);
                    v3.StoreUnsafe(ref rv, i + 3 * w1);
                }
                k = (int)i;
            }
#endif
            for (; k + w <= len; k += w)
                (new SimdVector(V, ov + k) - (va * new SimdVector(x, k) + vb * new SimdVector(y, oy + k))).CopyTo(V, ov + k);
            for (; k < len; k++) V[ov + k] -= a * x[k] + b * y[oy + k];
        }

        /// <summary>
        /// Adds a multiple of a vector to one row segment: <c>T[o + k] = T[o + k] + g[k]·s</c>.
        /// </summary>
        /// <param name="T">The working matrix, updated in place.</param>
        /// <param name="o">The offset of the row segment in <paramref name="T"/>.</param>
        /// <param name="g">The vector, from index zero.</param>
        /// <param name="s">The multiplier.</param>
        /// <param name="len">The number of entries.</param>
#if NET
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        internal static void AddScaled(double[] T, int o, double[] g, double s, int len)
        {
            int w = SimdVector.Count;
            var vs = new SimdVector(s);
            int k = 0;
#if NET
            if (len >= 4 * w)
            {
                ref double rt = ref Unsafe.Add(ref MemoryMarshal.GetArrayDataReference(T), o);
                ref double rg = ref MemoryMarshal.GetArrayDataReference(g);
                nuint i = 0, count = (nuint)len, w1 = (nuint)w, w4 = 4 * w1;
                for (; i + w4 <= count; i += w4)
                {
                    var t0 = System.Numerics.Vector.LoadUnsafe(ref rt, i) + System.Numerics.Vector.LoadUnsafe(ref rg, i) * vs;
                    var t1 = System.Numerics.Vector.LoadUnsafe(ref rt, i + w1) + System.Numerics.Vector.LoadUnsafe(ref rg, i + w1) * vs;
                    var t2 = System.Numerics.Vector.LoadUnsafe(ref rt, i + 2 * w1) + System.Numerics.Vector.LoadUnsafe(ref rg, i + 2 * w1) * vs;
                    var t3 = System.Numerics.Vector.LoadUnsafe(ref rt, i + 3 * w1) + System.Numerics.Vector.LoadUnsafe(ref rg, i + 3 * w1) * vs;
                    t0.StoreUnsafe(ref rt, i);
                    t1.StoreUnsafe(ref rt, i + w1);
                    t2.StoreUnsafe(ref rt, i + 2 * w1);
                    t3.StoreUnsafe(ref rt, i + 3 * w1);
                }
                k = (int)i;
            }
#endif
            for (; k + w <= len; k += w) (new SimdVector(T, o + k) + new SimdVector(g, k) * vs).CopyTo(T, o + k);
            for (; k < len; k++) T[o + k] += g[k] * s;
        }

        /// <summary>
        /// Writes the transpose of the n × n row-major matrix <paramref name="source"/> to <paramref name="target"/>.
        /// </summary>
        /// <param name="source">The matrix.</param>
        /// <param name="target">Output: the transpose. It must not be <paramref name="source"/>.</param>
        /// <param name="n">The order.</param>
        internal static void Transpose(double[] source, double[] target, int n)
        {
            const int block = 16;
            for (int i0 = 0; i0 < n; i0 += block)
            {
                int i1 = Math.Min(n, i0 + block);
                for (int j0 = 0; j0 < n; j0 += block)
                {
                    int j1 = Math.Min(n, j0 + block);
                    for (int i = i0; i < i1; i++)
                        for (int j = j0; j < j1; j++) target[j * n + i] = source[i * n + j];
                }
            }
        }

        /// <summary>
        /// Applies a log of QL sweeps to columns [c0, c1) of the rows of <paramref name="W"/>, sweep by sweep in log order.
        /// </summary>
        /// <param name="W">The n × n row-major matrix whose rows are rotated.</param>
        /// <param name="n">The order.</param>
        /// <param name="cosines">The cosine of each logged rotation.</param>
        /// <param name="sines">The sine of each logged rotation.</param>
        /// <param name="first">The first row l of each sweep.</param>
        /// <param name="last">The last row m of each sweep.</param>
        /// <param name="sweeps">The number of logged sweeps.</param>
        /// <param name="c0">The first column.</param>
        /// <param name="c1">One past the last column.</param>
        /// <remarks>
        /// Sweep q holds the rotations for i = m − 1 down to l, in that order. Rotation i replaces rows i and i + 1 by
        /// <c>row i+1 ← s·row i + c·row i+1</c> and <c>row i ← c·row i − s·row i+1</c>, the update of procedure tql2. The
        /// matrix is processed in blocks of columns; within a block the sweeps and rotations run in log order, and row
        /// i + 1 of the block is carried in registers from one rotation to the next. Each element therefore receives the
        /// same operations in the same order as when the rotations are applied one at a time to whole rows, while a block
        /// of columns stays in cache across all the logged sweeps.
        /// </remarks>
#if NET
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        internal static void ApplySweeps(double[] W, int n, double[] cosines, double[] sines, int[] first, int[] last, int sweeps, int c0, int c1)
        {
            int w = SimdVector.Count;
            int k = c0;
#if NET
            nuint w1 = (nuint)w, stride = (nuint)n;
            ref double rw = ref MemoryMarshal.GetArrayDataReference(W);
            for (; k + 4 * w <= c1; k += 4 * w)
            {
                int pos = 0;
                for (int q = 0; q < sweeps; q++)
                {
                    int l = first[q], m = last[q];
                    nuint o = (nuint)(m * n + k);
                    var y0 = System.Numerics.Vector.LoadUnsafe(ref rw, o);
                    var y1 = System.Numerics.Vector.LoadUnsafe(ref rw, o + w1);
                    var y2 = System.Numerics.Vector.LoadUnsafe(ref rw, o + 2 * w1);
                    var y3 = System.Numerics.Vector.LoadUnsafe(ref rw, o + 3 * w1);
                    for (int i = m - 1; i >= l; i--, pos++)
                    {
                        var vc = new SimdVector(cosines[pos]);
                        var vs = new SimdVector(sines[pos]);
                        nuint a = o - stride;
                        var x0 = System.Numerics.Vector.LoadUnsafe(ref rw, a);
                        var x1 = System.Numerics.Vector.LoadUnsafe(ref rw, a + w1);
                        var x2 = System.Numerics.Vector.LoadUnsafe(ref rw, a + 2 * w1);
                        var x3 = System.Numerics.Vector.LoadUnsafe(ref rw, a + 3 * w1);
                        (vs * x0 + vc * y0).StoreUnsafe(ref rw, o);
                        (vs * x1 + vc * y1).StoreUnsafe(ref rw, o + w1);
                        (vs * x2 + vc * y2).StoreUnsafe(ref rw, o + 2 * w1);
                        (vs * x3 + vc * y3).StoreUnsafe(ref rw, o + 3 * w1);
                        y0 = vc * x0 - vs * y0;
                        y1 = vc * x1 - vs * y1;
                        y2 = vc * x2 - vs * y2;
                        y3 = vc * x3 - vs * y3;
                        o = a;
                    }
                    y0.StoreUnsafe(ref rw, o);
                    y1.StoreUnsafe(ref rw, o + w1);
                    y2.StoreUnsafe(ref rw, o + 2 * w1);
                    y3.StoreUnsafe(ref rw, o + 3 * w1);
                }
            }
#endif
            for (; k + w <= c1; k += w)
            {
                int pos = 0;
                for (int q = 0; q < sweeps; q++)
                {
                    int l = first[q], m = last[q];
                    int o = m * n + k;
                    var y = new SimdVector(W, o);
                    for (int i = m - 1; i >= l; i--, pos++)
                    {
                        var vc = new SimdVector(cosines[pos]);
                        var vs = new SimdVector(sines[pos]);
                        var x = new SimdVector(W, o - n);
                        (vs * x + vc * y).CopyTo(W, o);
                        y = vc * x - vs * y;
                        o -= n;
                    }
                    y.CopyTo(W, o);
                }
            }
            for (; k < c1; k++)
            {
                int pos = 0;
                for (int q = 0; q < sweeps; q++)
                {
                    int l = first[q], m = last[q];
                    int o = m * n + k;
                    double y = W[o];
                    for (int i = m - 1; i >= l; i--, pos++)
                    {
                        double c = cosines[pos], s = sines[pos];
                        double x = W[o - n];
                        W[o] = s * x + c * y;
                        y = c * x - s * y;
                        o -= n;
                    }
                    W[o] = y;
                }
            }
        }
    }
}
