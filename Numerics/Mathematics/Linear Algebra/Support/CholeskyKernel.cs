using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;

namespace Numerics.Mathematics.LinearAlgebra
{
    /// <summary>
    /// Blocked dense kernels for the Cholesky factorization of symmetric positive-definite matrices stored as flat
    /// row-major arrays, with the triangular solves and the inverse that use the factor.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// Every routine is a right-looking blocked form of the textbook algorithm (Golub and Van Loan, Matrix
    /// Computations, 4th ed., sections 4.2 and 3.1), written as inner products of contiguous row segments so
    /// that the inner loops vectorize with 256-bit fused multiply-adds where the hardware supports them.
    /// A matrix is addressed by an offset and a leading dimension, so a block of a larger array can be factored
    /// and solved in place.
    /// </para>
    /// <para>
    /// Work is divided into tiles whose boundaries depend only on the matrix dimension, and every inner
    /// product is accumulated in a fixed order, so the results do not depend on the number of threads.
    /// They may differ in the last bits between machines with and without fused multiply-add.
    /// The factor has the usual backward error of Cholesky factorization: <c>L Lᵀ = A + E</c> with
    /// <c>|E| &lt;= c n ε |L| |Lᵀ|</c> for a modest constant <c>c</c>.
    /// </para>
    /// </remarks>
    internal static class CholeskyKernel
    {
        /// <summary>The block width of the panel operations.</summary>
        private const int BlockSize = 128;

        /// <summary>The dimension at and above which the blocked updates run in parallel.</summary>
        internal const int ParallelThreshold = 384;

        /// <summary>The length of the inner-product chunks of the triangular solve.</summary>
        private const int KChunk = 256;

        /// <summary>The number of rows in one parallel work item.</summary>
        private const int RowsPerTask = 8;

        #region Rank-k Update

        /// <summary>
        /// Accumulates the inner products of row segments for every row pair (i, j) with
        /// <c>colStart &lt;= j &lt;= i</c> and <c>rowStart &lt;= i &lt; rowEnd</c>.
        /// </summary>
        /// <param name="src">The row-major source.</param>
        /// <param name="srcOffset">The index of the source's element (0, 0).</param>
        /// <param name="dest">The row-major destination.</param>
        /// <param name="destOffset">The index of the destination's element (0, 0).</param>
        /// <param name="ld">The leading dimension of the source and destination.</param>
        /// <param name="rowStart">The first row i.</param>
        /// <param name="rowEnd">One past the last row i.</param>
        /// <param name="colStart">The first row j.</param>
        /// <param name="k0">The first column of the segments.</param>
        /// <param name="k1">One past the last column of the segments.</param>
        /// <param name="sign">The sign applied to each inner product before it is added.</param>
        /// <param name="transposed">
        /// When true, the product for (i, j) is added to <c>dest[j, i]</c> and diagonal products to
        /// <paramref name="diagonal"/>; otherwise it is added to <c>dest[i, j]</c>.
        /// </param>
        /// <param name="diagonal">The diagonal accumulator used when <paramref name="transposed"/> is true.</param>
        /// <param name="parallel">Whether the rows may be processed in parallel.</param>
        private static void RankUpdate(double[] src, int srcOffset, double[] dest, int destOffset, int ld, int rowStart, int rowEnd,
            int colStart, int k0, int k1, double sign, bool transposed, double[]? diagonal, bool parallel)
        {
            int rows = rowEnd - rowStart;
            if (rows <= 0 || k1 <= k0) return;
            int tasks = (rows + RowsPerTask - 1) / RowsPerTask;
            if (parallel && tasks > 1)
                RankUpdateParallel(src, srcOffset, dest, destOffset, ld, rowStart, rowEnd, colStart, k0, k1, sign, transposed, diagonal, tasks);
            else
                RankUpdateRows(src, srcOffset, dest, destOffset, ld, rowStart, rowEnd, colStart, k0, k1, sign, transposed, diagonal, new double[8]);
        }

        /// <summary>Applies <see cref="RankUpdate"/> in parallel bands of <see cref="RowsPerTask"/> rows.</summary>
        /// <param name="src">The source.</param>
        /// <param name="srcOffset">The index of the source's element (0, 0).</param>
        /// <param name="dest">The destination.</param>
        /// <param name="destOffset">The index of the destination's element (0, 0).</param>
        /// <param name="ld">The leading dimension.</param>
        /// <param name="rowStart">The first row i.</param>
        /// <param name="rowEnd">One past the last row i.</param>
        /// <param name="colStart">The first row j.</param>
        /// <param name="k0">The first segment column.</param>
        /// <param name="k1">One past the last segment column.</param>
        /// <param name="sign">The sign applied to each product.</param>
        /// <param name="transposed">Whether products are stored at the transposed position.</param>
        /// <param name="diagonal">The diagonal accumulator for the transposed form.</param>
        /// <param name="tasks">The number of bands.</param>
        /// <remarks>Kept apart from the serial path so that the serial path captures no variables in a closure.</remarks>
        private static void RankUpdateParallel(double[] src, int srcOffset, double[] dest, int destOffset, int ld, int rowStart, int rowEnd,
            int colStart, int k0, int k1, double sign, bool transposed, double[]? diagonal, int tasks)
        {
            Parallel.For(0, tasks, () => new double[8], (t, _, s) =>
            {
                int r0 = rowStart + t * RowsPerTask;
                RankUpdateRows(src, srcOffset, dest, destOffset, ld, r0, Math.Min(rowEnd, r0 + RowsPerTask), colStart, k0, k1, sign, transposed, diagonal, s);
                return s;
            }, _ => { });
        }

        /// <summary>Applies <see cref="RankUpdate"/> to one band of rows.</summary>
        /// <param name="src">The source.</param>
        /// <param name="srcOffset">The index of the source's element (0, 0).</param>
        /// <param name="dest">The destination.</param>
        /// <param name="destOffset">The index of the destination's element (0, 0).</param>
        /// <param name="ld">The leading dimension.</param>
        /// <param name="r0">The first row of the band.</param>
        /// <param name="r1">One past the last row of the band.</param>
        /// <param name="colStart">The first row j.</param>
        /// <param name="k0">The first segment column.</param>
        /// <param name="k1">One past the last segment column.</param>
        /// <param name="sign">The sign applied to each product.</param>
        /// <param name="transposed">Whether products are stored at the transposed position.</param>
        /// <param name="diagonal">The diagonal accumulator for the transposed form.</param>
        /// <param name="s">An eight-element scratch buffer.</param>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private static void RankUpdateRows(double[] src, int srcOffset, double[] dest, int destOffset, int ld, int r0, int r1,
            int colStart, int k0, int k1, double sign, bool transposed, double[]? diagonal, double[] s)
        {
            int len = k1 - k0;
            int i = r0;
            for (; i + 3 < r1; i += 4)
            {
                int ia = srcOffset + i * ld + k0;
                // Columns j < i are shared by all four rows; the triangle at the end is done one row at a time.
                int j = colStart;
                for (; j + 1 < i; j += 2)
                {
                    int ja = srcOffset + j * ld + k0;
                    CrossProductKernel.Dot4x2(src, ia, ia + ld, ia + 2 * ld, ia + 3 * ld, src, ja, ja + ld, len, s);
                    for (int r = 0; r < 4; r++)
                    {
                        Add(dest, destOffset, ld, i + r, j, sign * s[2 * r], transposed, diagonal);
                        Add(dest, destOffset, ld, i + r, j + 1, sign * s[2 * r + 1], transposed, diagonal);
                    }
                }
                for (int r = 0; r < 4; r++)
                {
                    int row = i + r;
                    for (int jj = j; jj <= row; jj++)
                        Add(dest, destOffset, ld, row, jj, sign * CrossProductKernel.Dot(src, srcOffset + row * ld + k0, src, srcOffset + jj * ld + k0, len), transposed, diagonal);
                }
            }
            for (; i < r1; i++)
            {
                for (int j = colStart; j <= i; j++)
                    Add(dest, destOffset, ld, i, j, sign * CrossProductKernel.Dot(src, srcOffset + i * ld + k0, src, srcOffset + j * ld + k0, len), transposed, diagonal);
            }
        }

        /// <summary>Adds a value at (i, j), or at its transposed position.</summary>
        /// <param name="dest">The destination.</param>
        /// <param name="offset">The index of the destination's element (0, 0).</param>
        /// <param name="ld">The leading dimension.</param>
        /// <param name="i">The row, with <c>j &lt;= i</c>.</param>
        /// <param name="j">The column.</param>
        /// <param name="value">The value to add.</param>
        /// <param name="transposed">Whether to add at (j, i), with the diagonal held separately.</param>
        /// <param name="diagonal">The diagonal accumulator for the transposed form.</param>
        private static void Add(double[] dest, int offset, int ld, int i, int j, double value, bool transposed, double[]? diagonal)
        {
            if (!transposed) dest[offset + i * ld + j] += value;
            else if (i == j) diagonal![i] += value;
            else dest[offset + j * ld + i] += value;
        }

        #endregion

        #region Factorization and Solves

        /// <summary>
        /// Overwrites the lower triangle of a symmetric matrix with its Cholesky factor.
        /// </summary>
        /// <param name="a">
        /// The storage of the row-major n-by-n matrix. Its lower triangle and diagonal are read and replaced by L;
        /// the strictly upper triangle is neither read nor written.
        /// </param>
        /// <param name="offset">The index of the matrix element (0, 0).</param>
        /// <param name="ld">The leading dimension, at least <paramref name="n"/>.</param>
        /// <param name="n">The dimension.</param>
        /// <param name="relativeTolerance">
        /// A pivot is rejected when it is NaN, not positive, or at most this fraction of the original diagonal
        /// entry, the test of <see cref="CholeskyDecomposition"/>. Zero gives the purely absolute test.
        /// </param>
        /// <param name="originalDiagonal">
        /// An n-element buffer that receives the original diagonal, or null to allocate one when
        /// <paramref name="relativeTolerance"/> is positive.
        /// </param>
        /// <param name="failedRow">The row of the rejected pivot, or -1 on success.</param>
        /// <param name="failedPivot">The rejected pivot before its square root, or NaN on success.</param>
        /// <returns>True when every pivot passes; otherwise the factor is incomplete.</returns>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        internal static bool TryFactor(double[] a, int offset, int ld, int n, double relativeTolerance, double[]? originalDiagonal,
            out int failedRow, out double failedPivot)
        {
            failedRow = -1;
            failedPivot = double.NaN;
            if (relativeTolerance > 0d)
            {
                originalDiagonal ??= new double[n];
                for (int i = 0; i < n; i++) originalDiagonal[i] = a[offset + i * ld + i];
            }
            else if (originalDiagonal != null)
            {
                for (int i = 0; i < n; i++) originalDiagonal[i] = a[offset + i * ld + i];
            }
            bool parallel = n >= ParallelThreshold;
            for (int k0 = 0; k0 < n; k0 += BlockSize)
            {
                int k1 = Math.Min(n, k0 + BlockSize);
                // Diagonal block.
                for (int i = k0; i < k1; i++)
                {
                    int ri = offset + i * ld;
                    for (int j = k0; j <= i; j++)
                    {
                        int rj = offset + j * ld;
                        double sum = a[ri + j] - CrossProductKernel.Dot(a, ri + k0, a, rj + k0, j - k0);
                        if (i == j)
                        {
                            double threshold = 0d;
                            if (relativeTolerance > 0d)
                            {
                                double d = originalDiagonal![i];
                                threshold = d > 0d && !double.IsInfinity(d) ? relativeTolerance * d : 0d;
                            }
                            if (double.IsNaN(sum) || sum <= 0d || sum <= threshold)
                            {
                                failedRow = i;
                                failedPivot = sum;
                                return false;
                            }
                            a[ri + i] = Math.Sqrt(sum);
                        }
                        else
                        {
                            a[ri + j] = sum / a[rj + j];
                        }
                    }
                }
                if (k1 == n) break;
                // Panel below the diagonal block; its rows are independent.
                int tasks = (n - k1 + RowsPerTask - 1) / RowsPerTask;
                if (parallel && tasks > 1) PanelParallel(a, offset, ld, n, k0, k1, tasks);
                else PanelRows(a, offset, ld, k0, k1, k1, n);
                // Trailing update of the lower triangle.
                RankUpdate(a, offset, a, offset, ld, k1, n, k1, k0, k1, -1d, false, null, parallel);
            }
            return true;
        }

        /// <summary>Computes rows [r0, r1) of the panel below the diagonal block of columns [k0, k1).</summary>
        /// <param name="a">The matrix storage.</param>
        /// <param name="offset">The index of the matrix element (0, 0).</param>
        /// <param name="ld">The leading dimension.</param>
        /// <param name="k0">The first column of the block.</param>
        /// <param name="k1">One past the last column of the block.</param>
        /// <param name="r0">The first row.</param>
        /// <param name="r1">One past the last row.</param>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private static void PanelRows(double[] a, int offset, int ld, int k0, int k1, int r0, int r1)
        {
            for (int i = r0; i < r1; i++)
            {
                int ri = offset + i * ld;
                for (int j = k0; j < k1; j++)
                {
                    int rj = offset + j * ld;
                    a[ri + j] = (a[ri + j] - CrossProductKernel.Dot(a, ri + k0, a, rj + k0, j - k0)) / a[rj + j];
                }
            }
        }

        /// <summary>Computes the panel below the diagonal block in parallel bands of <see cref="RowsPerTask"/> rows.</summary>
        /// <param name="a">The matrix storage.</param>
        /// <param name="offset">The index of the matrix element (0, 0).</param>
        /// <param name="ld">The leading dimension.</param>
        /// <param name="n">The dimension.</param>
        /// <param name="k0">The first column of the block.</param>
        /// <param name="k1">One past the last column of the block.</param>
        /// <param name="tasks">The number of bands.</param>
        private static void PanelParallel(double[] a, int offset, int ld, int n, int k0, int k1, int tasks)
        {
            Parallel.For(0, tasks, t =>
            {
                int r0 = k1 + t * RowsPerTask;
                PanelRows(a, offset, ld, k0, k1, r0, Math.Min(n, r0 + RowsPerTask));
            });
        }

        /// <summary>
        /// Solves <c>L x = b</c> in place for each row of <paramref name="x"/>.
        /// </summary>
        /// <param name="l">The storage of the row-major lower triangular factor.</param>
        /// <param name="offset">The index of the factor element (0, 0).</param>
        /// <param name="ld">The leading dimension of the factor.</param>
        /// <param name="n">The dimension of the factor.</param>
        /// <param name="x">
        /// The storage of the m-by-n row-major right-hand sides, with consecutive rows n elements apart, replaced by
        /// the solutions.
        /// </param>
        /// <param name="xOffset">The index of the first element of the first right-hand side.</param>
        /// <param name="m">The number of right-hand sides.</param>
        /// <param name="start">
        /// Optional. For each right-hand side, the index of its first possibly nonzero entry; entries before
        /// it must be zero and stay zero in the solution.
        /// </param>
        internal static void SolveLowerRows(double[] l, int offset, int ld, int n, double[] x, int xOffset, int m, int[]? start)
        {
            if (m <= 0) return;
            if (m == 1 && n <= BlockSize)
            {
                // One diagonal block: plain forward substitution, as the general path below computes it.
                SubstituteWithinBlock(l, offset, ld, n, x, xOffset, 0, start, 0, n);
                return;
            }
            // Parallel when the dimension is large, or when many right-hand sides make the few barriers cheap.
            bool parallel = m >= 2 * RowsPerTask && (n >= ParallelThreshold || (long)m * n * n >= 1L << 22);
            int tasks = (m + RowsPerTask - 1) / RowsPerTask;
            if (parallel)
            {
                SolveLowerRowsParallel(l, offset, ld, n, x, xOffset, m, start, tasks);
                return;
            }
            // The scratch buffer is used only by blocks after the first.
            var s = n > BlockSize ? new double[8] : Array.Empty<double>();
            for (int i0 = 0; i0 < n; i0 += BlockSize)
            {
                int i1 = Math.Min(n, i0 + BlockSize);
                for (int t = 0; t < tasks; t++)
                    SolveLowerBlock(l, offset, ld, n, x, xOffset, t * RowsPerTask, Math.Min(m, t * RowsPerTask + RowsPerTask), start, i0, i1, s);
            }
        }

        /// <summary>
        /// Runs <see cref="SolveLowerRows"/> with the right-hand sides in parallel bands of <see cref="RowsPerTask"/>,
        /// one diagonal block at a time.
        /// </summary>
        /// <param name="l">The storage of the row-major lower triangular factor.</param>
        /// <param name="offset">The index of the factor element (0, 0).</param>
        /// <param name="ld">The leading dimension of the factor.</param>
        /// <param name="n">The dimension of the factor.</param>
        /// <param name="x">The right-hand sides, replaced by the solutions.</param>
        /// <param name="xOffset">The index of the first element of the first right-hand side.</param>
        /// <param name="m">The number of right-hand sides.</param>
        /// <param name="start">The optional first nonzero index of each right-hand side.</param>
        /// <param name="tasks">The number of bands.</param>
        private static void SolveLowerRowsParallel(double[] l, int offset, int ld, int n, double[] x, int xOffset, int m, int[]? start, int tasks)
        {
            for (int i0 = 0; i0 < n; i0 += BlockSize)
            {
                int block = i0, blockEnd = Math.Min(n, i0 + BlockSize);
                Parallel.For(0, tasks, () => new double[8], (t, _, s) =>
                {
                    int p0 = t * RowsPerTask, p1 = Math.Min(m, p0 + RowsPerTask);
                    SolveLowerBlock(l, offset, ld, n, x, xOffset, p0, p1, start, block, blockEnd, s);
                    return s;
                }, _ => { });
            }
        }

        /// <summary>Advances the forward substitution of rows [p0, p1) through one diagonal block.</summary>
        /// <param name="l">The lower triangular factor.</param>
        /// <param name="offset">The index of the factor element (0, 0).</param>
        /// <param name="ld">The leading dimension of the factor.</param>
        /// <param name="n">The dimension, which is also the distance between right-hand sides.</param>
        /// <param name="x">The right-hand sides.</param>
        /// <param name="xOffset">The index of the first element of the first right-hand side.</param>
        /// <param name="p0">The first right-hand side.</param>
        /// <param name="p1">One past the last right-hand side.</param>
        /// <param name="start">The optional first nonzero index of each right-hand side.</param>
        /// <param name="i0">The first row of the block.</param>
        /// <param name="i1">One past the last row of the block.</param>
        /// <param name="s">An eight-element scratch buffer.</param>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private static void SolveLowerBlock(double[] l, int offset, int ld, int n, double[] x, int xOffset, int p0, int p1, int[]? start,
            int i0, int i1, double[] s)
        {
            // Update with the solved part: x[p, i] -= L[i, k] x[p, k] over k < i0. The k range is taken in
            // chunks so that four rows of L and the right-hand sides of this band stay in cache while every
            // pair of right-hand sides is updated against them. Entries before a right-hand side's start are
            // zero, so a pair may begin at the smaller of its two starts.
            int pairs = (p1 - p0) / 2;
            for (int kc = 0; kc < i0; kc += KChunk)
            {
                int ke = Math.Min(i0, kc + KChunk);
                int i = i0;
                for (; i + 3 < i1; i += 4)
                {
                    int li = offset + i * ld;
                    for (int t = 0; t < pairs; t++)
                    {
                        int p = p0 + 2 * t;
                        int from = Math.Max(kc, start == null ? 0 : Math.Min(start[p], start[p + 1]));
                        if (from >= ke) continue;
                        int xp = xOffset + p * n, xq = xp + n;
                        CrossProductKernel.Dot4x2(l, li + from, li + ld + from, li + 2 * ld + from, li + 3 * ld + from, x, xp + from, xq + from, ke - from, s);
                        for (int r = 0; r < 4; r++)
                        {
                            x[xp + i + r] -= s[2 * r];
                            x[xq + i + r] -= s[2 * r + 1];
                        }
                    }
                }
                for (int p = p0; p < p1; p++)
                {
                    int from = Math.Max(kc, start == null ? 0 : start[p]);
                    if (from >= ke) continue;
                    int xp = xOffset + p * n;
                    // Rows not covered by the groups of four, and every row for an unpaired right-hand side.
                    int first = p < p0 + 2 * pairs ? i : i0;
                    for (int r = first; r < i1; r++)
                        x[xp + r] -= CrossProductKernel.Dot(l, offset + r * ld + from, x, xp + from, ke - from);
                }
            }
            for (int p = p0; p < p1; p++)
                SubstituteWithinBlock(l, offset, ld, n, x, xOffset, p, start, i0, i1);
        }

        /// <summary>Completes the forward substitution of one right-hand side inside a diagonal block.</summary>
        /// <param name="l">The lower triangular factor.</param>
        /// <param name="offset">The index of the factor element (0, 0).</param>
        /// <param name="ld">The leading dimension of the factor.</param>
        /// <param name="n">The dimension, which is also the distance between right-hand sides.</param>
        /// <param name="x">The right-hand sides.</param>
        /// <param name="xOffset">The index of the first element of the first right-hand side.</param>
        /// <param name="p">The right-hand side.</param>
        /// <param name="start">The optional first nonzero index of each right-hand side.</param>
        /// <param name="i0">The first row of the block.</param>
        /// <param name="i1">One past the last row of the block.</param>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private static void SubstituteWithinBlock(double[] l, int offset, int ld, int n, double[] x, int xOffset, int p, int[]? start, int i0, int i1)
        {
            int xp = xOffset + p * n;
            int first = start == null ? i0 : Math.Max(i0, start[p]);
            for (int i = first; i < i1; i++)
            {
                int li = offset + i * ld;
                x[xp + i] = (x[xp + i] - CrossProductKernel.Dot(l, li + first, x, xp + first, i - first)) / l[li + i];
            }
        }

        /// <summary>
        /// Solves <c>Lᵀ x = b</c> in place for one right-hand side.
        /// </summary>
        /// <param name="l">The storage of the row-major lower triangular factor.</param>
        /// <param name="offset">The index of the factor element (0, 0).</param>
        /// <param name="ld">The leading dimension of the factor.</param>
        /// <param name="n">The dimension.</param>
        /// <param name="x">The storage of the right-hand side, replaced by the solution.</param>
        /// <param name="xOffset">The index of the first element of the right-hand side.</param>
        /// <remarks>
        /// The solve is column oriented: once x[i] is known, its multiple of row i of L is subtracted from
        /// x[0..i), so every update runs over a contiguous row.
        /// </remarks>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        internal static void SolveUpper(double[] l, int offset, int ld, int n, double[] x, int xOffset)
        {
            for (int i = n - 1; i >= 0; i--)
            {
                int ri = offset + i * ld;
                double xi = x[xOffset + i] / l[ri + i];
                x[xOffset + i] = xi;
                CrossProductKernel.Axpy(-xi, l, ri, x, xOffset, i);
            }
        }

        /// <summary>
        /// Computes the inverse of the factored matrix, <c>A⁻¹ = L⁻ᵀ L⁻¹</c>.
        /// </summary>
        /// <param name="l">
        /// The row-major n-by-n factor in the lower triangle. The strictly upper triangle receives the off-diagonal
        /// entries of the inverse, with <c>A⁻¹[i, j]</c> for <c>i &lt; j</c> stored at <c>l[i, j]</c>.
        /// </param>
        /// <param name="n">The dimension.</param>
        /// <param name="t">An n-by-n buffer that receives <c>L⁻ᵀ</c> in row-major order (upper triangular).</param>
        /// <param name="inverseDiagonal">An n-element buffer that receives the diagonal of the inverse.</param>
        internal static void Invert(double[] l, int n, double[] t, double[] inverseDiagonal)
        {
            InvertFactor(l, n, t);
            // The inverse entry (a, b) is the inner product of rows a and b of L⁻ᵀ over k >= max(a, b),
            // accumulated one column block at a time in increasing order.
            Array.Clear(inverseDiagonal, 0, n);
            for (int i = 0; i < n; i++)
                Array.Clear(l, i * n + i + 1, n - i - 1);
            bool parallel = n >= ParallelThreshold;
            for (int k0 = 0; k0 < n; k0 += BlockSize)
            {
                int k1 = Math.Min(n, k0 + BlockSize);
                RankUpdate(t, 0, l, 0, n, 0, k1, 0, k0, k1, 1d, true, inverseDiagonal, parallel);
            }
        }

        /// <summary>
        /// Computes <c>L⁻ᵀ</c>, whose row j is the solution of <c>L x = e_j</c>.
        /// </summary>
        /// <param name="l">The row-major n-by-n lower triangular factor.</param>
        /// <param name="n">The dimension.</param>
        /// <param name="t">An n-by-n buffer that receives <c>L⁻ᵀ</c> in row-major order.</param>
        internal static void InvertFactor(double[] l, int n, double[] t)
        {
            Array.Clear(t, 0, n * n);
            var start = new int[n];
            for (int j = 0; j < n; j++)
            {
                t[j * n + j] = 1d;
                start[j] = j;
            }
            SolveLowerRows(l, 0, n, n, t, 0, n, start);
        }

        #endregion

    }
}
