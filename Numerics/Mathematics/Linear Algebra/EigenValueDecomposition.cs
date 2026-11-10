namespace Numerics.Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class for computing all eigenvalues and eigenvectors of a real <b>symmetric</b> matrix using the cyclic Jacobi rotation
    /// method or Householder tridiagonalization with the implicit QL algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b> Authors: </b>
    /// <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// Computes the eigen-decomposition A = V * D * V^T for a real symmetric matrix A, where D contains the eigenvalues and
    /// the columns of V are the corresponding orthonormal eigenvectors. The implementation uses the cyclic Jacobi method:
    /// each sweep rotates, row by row, every off-diagonal element above the convergence threshold, until a sweep finds
    /// none. Each sweep costs O(n³); ill-conditioned matrices, such as spectra graded down to 1E-12 of the largest
    /// eigenvalue, take more sweeps than well-separated ones. The threshold is 1E-12 times the largest element magnitude.
    /// On random Gram matrices and on matrices with evenly spaced or graded known spectra, for n = 10 to 300, the
    /// eigenpair residual ‖A v − λ v‖₂ / ‖A‖_F is 2E-13 to 3E-12, the eigenvectors are orthonormal to 1.5E-14, and the
    /// eigenvalues agree with the known spectrum to 1E-14 of the largest eigenvalue, or 1.2E-12 for a spectrum graded
    /// down to 1E-12. The eigenvalues are left in diagonal order and are not sorted.
    /// </para>
    /// <para>
    /// With <see cref="EigenSolverMethod.QL"/> the matrix is reduced to tridiagonal form by Householder reflections (tred2)
    /// and diagonalized by implicit QL iterations (tql2), the EISPACK procedures as given in JAMA, at about 9n³ flops;
    /// the eigenvalues are returned in ascending order, as LAPACK's dsyevd and numpy.linalg.eigh return them, with the
    /// eigenvectors in the matching columns. The eigenvalues are accurate to a small multiple of ε‖A‖, and the
    /// orthogonality error ‖VᵀV − I‖ of the eigenvectors grows with n, for random Gram matrices from 1.4E-15 at n = 10 to
    /// 6.8E-14 at n = 500 (LAPACK's divide-and-conquer dsyevd: 2.7E-15 at n = 500). That is the accuracy of the cited
    /// EISPACK procedures, which are implemented without modification.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes: The Art of Scientific Computing, Third Edition." Press et al., 2017. </item>
    /// <item> <description><see href = "https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm"/></description> </item>
    /// <item> Martin, R. S., Reinsch, C., and Wilkinson, J. H. (1968). Householder's tridiagonalization of a symmetric
    /// matrix. Numerische Mathematik, 11(3), 181-195. </item>
    /// <item> Bowdler, H., Martin, R. S., Reinsch, C., and Wilkinson, J. H. (1968). The QR and QL algorithms for
    /// symmetric matrices. Numerische Mathematik, 11(4), 293-306. </item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class EigenValueDecomposition
    {
        /// <summary>
        /// Constructs a new symmetric eigen-decomposition by the cyclic Jacobi method.
        /// </summary>
        /// <param name="A">The <b>symmetric</b> input matrix A [0..n-1][0..n-1].</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when A is not square.</exception>
        /// <exception cref="ArgumentException">Thrown when A contains a NaN or infinite entry, or is not symmetric.</exception>
        public EigenValueDecomposition(Matrix A) : this(A, EigenSolverMethod.Jacobi)
        {
        }

        /// <summary>
        /// Constructs a new symmetric eigen-decomposition by the given method.
        /// </summary>
        /// <param name="A">The <b>symmetric</b> input matrix A [0..n-1][0..n-1].</param>
        /// <param name="method">The algorithm: <see cref="EigenSolverMethod.Jacobi"/>, with the eigenvalues in diagonal
        /// order, or <see cref="EigenSolverMethod.QL"/>, with the eigenvalues in ascending order.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when A is not square, or <paramref name="method"/> is not a
        /// defined value.</exception>
        /// <exception cref="ArgumentException">Thrown when A contains a NaN or infinite entry, or is not symmetric.</exception>
        /// <exception cref="ArithmeticException">Thrown by <see cref="EigenSolverMethod.QL"/> when an eigenvalue does not
        /// converge within 60 QL iterations.</exception>
        /// <remarks>
        /// <see cref="EigenSolverMethod.Jacobi"/> gives results bitwise identical to <see cref="EigenValueDecomposition(Matrix)"/>.
        /// For <see cref="EigenSolverMethod.QL"/>, column i of <see cref="EigenVectors"/> is the unit eigenvector of the i-th
        /// smallest eigenvalue, and <see cref="IsConverged"/> is true on return.
        /// </remarks>
        public EigenValueDecomposition(Matrix A, EigenSolverMethod method)
        {
            if (!Enum.IsDefined(typeof(EigenSolverMethod), method))
                throw new ArgumentOutOfRangeException(nameof(method), "Unknown eigen solver method.");
            if (A.NumberOfRows != A.NumberOfColumns)
                throw new ArgumentOutOfRangeException(nameof(A), "The matrix A must be square.");
            for (int i = 0; i < A.NumberOfRows; i++)
            {
                for (int j = 0; j < A.NumberOfColumns; j++)
                {
                    if (!Tools.IsFinite(A[i, j]))
                        throw new ArgumentException("The matrix A must contain only finite entries.", nameof(A));
                }
            }
            if (!A.IsSymmetric)
                throw new ArgumentException("The matrix A must be symmetric for this decomposition.", nameof(A));

            n = A.NumberOfRows;
            this.A = new Matrix(A.Array);
            EigenVectors = Matrix.Identity(n);
            EigenValues = new Vector(n);
            if (method == EigenSolverMethod.QL)
            {
                DecomposeAscending(A);
                return;
            }

            // Rotate a working copy so the public input matrix keeps the values that were decomposed
            var work = new Matrix(A.Array);
            var a = work.Array; // same storage as work

            // The convergence threshold follows the scale of the matrix, so a decomposition of c*A
            // stops at the same relative accuracy as a decomposition of A. The largest element
            // magnitude is the reference: it dominates every off-diagonal, and rotations can always
            // drive the off-diagonals below 1E-12 of it before reaching roundoff.
            double elementMax = 0.0;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    elementMax = Math.Max(elementMax, Math.Abs(a[i, j]));
            double tol = 1e-12 * elementMax;
            // Cyclic Jacobi: each sweep visits the off-diagonal elements row by row and rotates every one above the
            // threshold. A sweep is at most n(n - 1)/2 rotations. Ill-conditioned matrices take up to about 18 sweeps
            // (about 9 n^2 rotations), so the limit of 50 n^2, about 100 sweeps with a floor of 2000, leaves a wide
            // margin; a matrix that reaches it is reported through IsConverged.
            long maxIter = Math.Max(2000L, 50L * n * n);
            long rotations = 0;
            var v = EigenVectors.Array;
            while (rotations < maxIter)
            {
                bool rotated = false;
                for (int p = 0; p < n - 1 && rotations < maxIter; p++)
                {
                    for (int q = p + 1; q < n && rotations < maxIter; q++)
                    {
                        double apq = a[p, q];
                        if (Math.Abs(apq) <= tol) continue;
                        rotated = true;
                        rotations++;

                        double app = a[p, p];
                        double aqq = a[q, q];

                        // Robust angle computation: θ = 0.5 * atan2(2*apq, aqq - app)
                        double theta = 0.5 * Math.Atan2(2.0 * apq, aqq - app);
                        double c = Math.Cos(theta);
                        double s = Math.Sin(theta);

                        // Rotate the 2x2 pivot block exactly
                        double appNew = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                        double aqqNew = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                        a[p, p] = appNew;
                        a[q, q] = aqqNew;
                        a[p, q] = a[q, p] = 0.0;

                        // Update the rest of the rows/cols (preserve symmetry)
                        for (int k = 0; k < n; k++)
                        {
                            if (k == p || k == q) continue;
                            double aik = a[k, p];
                            double akq = a[k, q];
                            double akpNew = c * aik - s * akq;
                            double akqNew = s * aik + c * akq;
                            a[k, p] = a[p, k] = akpNew;
                            a[k, q] = a[q, k] = akqNew;
                        }

                        // Accumulate eigenvectors: V = V * J(p,q)
                        for (int k = 0; k < n; k++)
                        {
                            double vip = v[k, p];
                            double viq = v[k, q];
                            v[k, p] = c * vip - s * viq;
                            v[k, q] = s * vip + c * viq;
                        }
                    }
                }
                if (!rotated)
                {
                    IsConverged = true;
                    break;
                }
            }

            // Extract eigenvalues from the diagonal of the rotated working copy
            for (int i = 0; i < n; i++)
                EigenValues[i] = a[i, i];
        }

        private readonly int n; // Size of A

        /// <summary>
        /// Fills the eigenvalues in ascending order and the matching eigenvector columns by tred2 and tql2.
        /// </summary>
        /// <param name="A">The validated symmetric matrix.</param>
        /// <exception cref="ArithmeticException">Thrown when the QL iteration does not converge.</exception>
        private void DecomposeAscending(Matrix A)
        {
            if (n > 0)
            {
                var flat = new double[n * n];
                Buffer.BlockCopy(A.Array, 0, flat, 0, flat.Length * sizeof(double));
                DecomposeTridiagonal(flat, n, out var values, out var vectors);
                var v = EigenVectors.Array;
                for (int i = 0; i < n; i++)
                {
                    // The tridiagonal solver returns decreasing values with one eigenvector per row.
                    int source = n - 1 - i;
                    EigenValues[i] = values[source];
                    for (int k = 0; k < n; k++) v[k, i] = vectors[source * n + k];
                }
            }
            IsConverged = true;
        }

        /// <summary>
        /// Stores the input matrix A that was decomposed (copied from the constructor input).
        /// </summary>
        public Matrix A { get; private set; }

        /// <summary>
        /// Returns the vector of eigenvalues (length n).
        /// </summary>
        public Vector EigenValues { get; private set; }

        /// <summary>
        /// Gets whether the Jacobi iteration converged, that is, every off-diagonal element fell to at most
        /// 1E-12 times the largest element magnitude within the rotation limit. Always true after a decomposition by
        /// <see cref="EigenSolverMethod.QL"/>, which throws instead of returning unconverged values.
        /// </summary>
        /// <remarks>
        /// The decomposition stops after max(2000, 50 n^2) rotations, about 100 cyclic sweeps. Random Gram matrices and
        /// spectra graded down to 1E-12, for n = 10 to 300, converge within the limit. When this property is false, the
        /// eigenvalues and eigenvectors are those of the last rotation and may be inaccurate. Each sweep costs O(n^3).
        /// </remarks>
        public bool IsConverged { get; private set; }

        /// <summary>
        /// Returns the matrix of eigenvectors (n x n). Each column corresponds to an eigenvalue.
        /// </summary>
        public Matrix EigenVectors { get; private set; }

        /// <summary>
        /// Returns the effective sample size based on Dutilleul's method (1993). 
        /// </summary>
        /// <returns>The scale-invariant effective sample size, or zero for an all-zero spectrum.</returns>
        /// <remarks>
        /// Eigenvalues are normalized by the largest absolute eigenvalue before evaluating
        /// <c>(Σλ)² / Σλ²</c>. Negative eigenvalues no larger than <c>1E-10</c> of that spectral
        /// scale are treated as numerical roundoff; materially negative eigenvalues retain the
        /// established formula behavior.
        /// </remarks>
        public double EffectiveSampleSize()
        {
            double spectralScale = 0d;
            for (int i = 0; i < EigenValues.Length; i++)
                spectralScale = Math.Max(spectralScale, Math.Abs(EigenValues[i]));
            if (spectralScale == 0d) return 0d;

            double sum = 0;
            double sumsq = 0;
            for (int i = 0; i < EigenValues.Length; i++)
            {
                // Normalize first so the roundoff threshold and ESS are independent of matrix scale.
                double lambda = EigenValues[i] / spectralScale;
                if (lambda < 0.0 && Math.Abs(lambda) <= 1e-10) lambda = 0.0;
                sum += lambda;
                sumsq += lambda * lambda;
            }
            if (sumsq == 0d) return 0d;
            double neff = (sum * sum) / sumsq;
            return neff;
        }

        /// <summary>
        /// Decomposes the symmetric matrix A = V·diag(λ)·Vᵀ on a flat row-major array by Householder tridiagonalization
        /// and the implicit QL algorithm.
        /// </summary>
        /// <param name="a">The p × p symmetric matrix, row-major. It is overwritten.</param>
        /// <param name="p">The order of the matrix.</param>
        /// <param name="values">Output: the p eigenvalues in decreasing order.</param>
        /// <param name="vectors">Output: p × p, row i holds the unit eigenvector of eigenvalue i.</param>
        /// <exception cref="ArithmeticException">Thrown when the QL iteration does not converge.</exception>
        /// <remarks>
        /// <para>
        /// The matrix is reduced to tridiagonal form by Householder reflections with the transformations accumulated
        /// (procedure tred2), and the tridiagonal matrix is diagonalized by QL iterations with implicit shifts (procedure
        /// tql2), the EISPACK procedures in the form given by JAMA. The cost is about 9p³ flops on flat arrays. The loops
        /// run lane-wise in SIMD registers and the QL rotations are applied in batches of sweeps, with every element
        /// keeping the operations and the order of operations of the sequential procedures, so the result is bitwise
        /// that of the sequential procedures and does not depend on the SIMD width. From p = 128 the logged rotations are
        /// applied by parallel tasks, in independent blocks of columns, while the next sweeps are computed; this changes no
        /// bits either. Eigenvalues are accurate to a small multiple of ε·‖A‖; the orthogonality error of the eigenvectors
        /// grows with p (6.8E-14 at p = 500 on random Gram matrices), the accuracy of the cited procedures.
        /// </para>
        /// <para>
        /// References: Bowdler, H., Martin, R. S., Reinsch, C., and Wilkinson, J. H. (1968). The QR and QL algorithms for
        /// symmetric matrices. Numerische Mathematik, 11(4), 293-306. Wilkinson, J. H., and Reinsch, C. (1971). Handbook
        /// for Automatic Computation, Volume II: Linear Algebra. Springer, Berlin.
        /// </para>
        /// </remarks>
        internal static void DecomposeTridiagonal(double[] a, int p, out double[] values, out double[] vectors)
        {
            var V = a;
            var d = new double[p];
            var e = new double[p];
            // The reduction leaves the transformation transposed, so each eigenvector is a contiguous row for the QL rotations.
            Tridiagonalize(V, d, e, p);
            Diagonalize(V, d, e, p);

            // Sort into decreasing order, keeping the original order of ties.
            var order = new int[p];
            for (int i = 0; i < p; i++) order[i] = i;
            Array.Sort(order, (x, y) => d[y] != d[x] ? d[y].CompareTo(d[x]) : x.CompareTo(y));
            values = new double[p];
            vectors = new double[p * p];
            for (int i = 0; i < p; i++)
            {
                values[i] = d[order[i]];
                Array.Copy(V, order[i] * p, vectors, i * p, p);
            }
        }

        /// <summary>
        /// Householder reduction to tridiagonal form with accumulated transformations (tred2), on the transposed layout.
        /// On exit row i of V is column i of the orthogonal transformation, d is the diagonal and e[1..p-1] the
        /// subdiagonal.
        /// </summary>
        /// <param name="V">The symmetric matrix, row-major, overwritten by the transposed transformation.</param>
        /// <param name="d">Output: the diagonal.</param>
        /// <param name="e">Output: the subdiagonal in e[1..n-1].</param>
        /// <param name="n">The order of the matrix.</param>
        /// <remarks>
        /// <para>
        /// The procedure of Martin, Reinsch and Wilkinson (1968) reads the lower triangle by columns. Here the same
        /// arithmetic in the same order runs on the transpose, where those columns are contiguous rows, so the input is
        /// read from its upper triangle; the transformation then comes out transposed, which is the row layout the QL
        /// rotations use. For an exactly symmetric input the result is bitwise that of the column-indexed procedure
        /// followed by a transpose; for an input symmetric only to rounding, the two read different triangles and agree
        /// to rounding.
        /// </para>
        /// <para>
        /// Every inner product of the procedure is a sum in a fixed order, and every update is element-wise, so the
        /// loops are rearranged to run lane-wise in SIMD registers with each element keeping its order of operations.
        /// The reduction keeps the whole active block symmetric (the upper triangle is first mirrored into the lower
        /// one, and the rank-two update is applied to full rows), so the product A·u is formed as row axpys
        /// (<see cref="TridiagonalKernel.ColumnSums"/>) with each entry summed in the procedure's order. The two terms of
        /// the update are added in either order, which is exact. The accumulation of the transformation runs on the
        /// transpose of the working array, where its inner products are also row axpys. These rearrangements change no
        /// bits.
        /// </para>
        /// </remarks>
        private static void Tridiagonalize(double[] V, double[] d, double[] e, int n)
        {
            for (int r = 1; r < n; r++)
                for (int c = 0; c < r; c++) V[r * n + c] = V[c * n + r];
            for (int j = 0; j < n; j++) d[j] = V[j * n + n - 1];
            for (int i = n - 1; i > 0; i--)
            {
                double scale = 0d, h = 0d;
                for (int k = 0; k < i; k++) scale += Math.Abs(d[k]);
                if (scale == 0d)
                {
                    e[i] = d[i - 1];
                    for (int j = 0; j < i; j++)
                    {
                        d[j] = V[j * n + i - 1];
                        V[j * n + i] = 0d;
                        V[i * n + j] = 0d;
                    }
                }
                else
                {
                    for (int k = 0; k < i; k++)
                    {
                        d[k] /= scale;
                        h += d[k] * d[k];
                    }
                    double f = d[i - 1];
                    double g = Math.Sqrt(h);
                    if (f > 0) g = -g;
                    e[i] = scale * g;
                    h -= f * g;
                    d[i - 1] = f - g;
                    // Row i keeps the Householder vector u; e = A·u over the leading i × i block.
                    int u = i * n;
                    Array.Copy(d, 0, V, u, i);
                    TridiagonalKernel.ColumnSums(V, n, d, e, i);
                    f = 0d;
                    for (int j = 0; j < i; j++)
                    {
                        e[j] /= h;
                        f += e[j] * d[j];
                    }
                    double hh = f / (h + h);
                    for (int j = 0; j < i; j++) e[j] -= hh * d[j];
                    for (int j = 0; j < i; j++)
                    {
                        // A[j, k] -= u[j]·e[k] + e[j]·u[k] on the whole row keeps the block symmetric.
                        TridiagonalKernel.SubtractRankTwo(V, j * n, d[j], e, e[j], V, u, i);
                        d[j] = V[j * n + i - 1];
                        V[j * n + i] = 0d;
                    }
                }
                d[i] = h;
            }

            // Accumulate the transformation on the transpose T, where row k of T is column k of the working array.
            var T = new double[n * n];
            TridiagonalKernel.Transpose(V, T, n);
            var x = new double[n];
            var g2 = new double[n];
            for (int i = 0; i < n - 1; i++)
            {
                T[(n - 1) * n + i] = T[i * n + i];
                T[i * n + i] = 1d;
                double h = d[i + 1];
                if (h != 0d)
                {
                    for (int k = 0; k <= i; k++)
                    {
                        x[k] = T[k * n + i + 1];
                        d[k] = x[k] / h;
                    }
                    TridiagonalKernel.ColumnSums(T, n, x, g2, i + 1);
                    for (int j = 0; j <= i; j++) g2[j] = -g2[j];
                    for (int k = 0; k <= i; k++) TridiagonalKernel.AddScaled(T, k * n, g2, d[k], i + 1);
                }
                for (int k = 0; k <= i; k++) T[k * n + i + 1] = 0d;
            }
            for (int j = 0; j < n; j++)
            {
                d[j] = T[(n - 1) * n + j];
                T[(n - 1) * n + j] = 0d;
            }
            T[(n - 1) * n + n - 1] = 1d;
            e[0] = 0d;
            TridiagonalKernel.Transpose(T, V, n);
        }

        /// <summary>
        /// Implicit QL iterations on the tridiagonal matrix (tql2). W holds the eigenvectors as rows.
        /// </summary>
        /// <param name="W">The transposed transformation from <see cref="Tridiagonalize"/>; on exit row i is the
        /// eigenvector of d[i].</param>
        /// <param name="d">The diagonal; on exit the eigenvalues, unsorted.</param>
        /// <param name="e">The subdiagonal in e[1..n-1]; destroyed.</param>
        /// <param name="n">The order of the matrix.</param>
        /// <exception cref="ArithmeticException">Thrown when an eigenvalue does not converge within 60 iterations.</exception>
        /// <remarks>
        /// The scalar recurrence is that of procedure tql2 (Bowdler et al., 1968) as given in JAMA. The plane rotations
        /// of each sweep depend only on d and e, so they are logged and applied to W in batches by
        /// <see cref="TridiagonalKernel.ApplySweeps"/>, which gives every element of W the same operations in the same
        /// order as rotating whole rows one at a time; the result is bitwise unchanged. From order
        /// <see cref="ParallelOrder"/> two logs alternate: a full log is applied by parallel tasks, each on its own block
        /// of columns, while the scalar recurrence fills the other, and each log waits for the previous one, so the sweeps
        /// still reach every column in order.
        /// </remarks>
        private static void Diagonalize(double[] W, double[] d, double[] e, int n)
        {
            bool parallel = n >= ParallelOrder;
            // The log holds a few sweeps; a full log is applied and reused, which changes no bits. For large orders two
            // logs alternate, so one is applied by parallel tasks while the next is being computed.
            int capacity = parallel ? Math.Max(n, ParallelLogCapacity) : Math.Max(n, Math.Min(1024, 2 * n * n));
            var log = new RotationLog(capacity);
            var spare = parallel ? new RotationLog(capacity) : null;
            Task? pending = null;

            for (int i = 1; i < n; i++) e[i - 1] = e[i];
            e[n - 1] = 0d;
            double f = 0d, tst1 = 0d;
            double eps = Math.Pow(2d, -52d);
            for (int l = 0; l < n; l++)
            {
                tst1 = Math.Max(tst1, Math.Abs(d[l]) + Math.Abs(e[l]));
                int m = l;
                while (m < n)
                {
                    if (Math.Abs(e[m]) <= eps * tst1) break;
                    m++;
                }
                if (m > l)
                {
                    int iteration = 0;
                    do
                    {
                        if (++iteration > 60)
                        {
                            pending?.Wait();
                            throw new ArithmeticException("The QL iteration did not converge.");
                        }
                        double g = d[l];
                        double p = (d[l + 1] - g) / (2d * e[l]);
                        double r = Tools.Hypot(p, 1d);
                        if (p < 0) r = -r;
                        d[l] = e[l] / (p + r);
                        d[l + 1] = e[l] * (p + r);
                        double dl1 = d[l + 1];
                        double h = g - d[l];
                        for (int i = l + 2; i < n; i++) d[i] -= h;
                        f += h;

                        p = d[m];
                        double c = 1d, c2 = c, c3 = c;
                        double el1 = e[l + 1];
                        double s = 0d, s2 = 0d;
                        if (log.Rotations + (m - l) > capacity)
                        {
                            if (spare == null)
                            {
                                log.Apply(W, n, false);
                            }
                            else
                            {
                                pending?.Wait();
                                var full = log;
                                pending = Task.Run(() => full.Apply(W, n, true));
                                log = spare;
                                spare = full;
                            }
                            log.Rotations = 0;
                            log.Sweeps = 0;
                        }
                        double[] cosines = log.Cosines, sines = log.Sines;
                        int logged = log.Rotations;
                        for (int i = m - 1; i >= l; i--)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * e[i];
                            h = c * p;
                            r = Tools.Hypot(p, e[i]);
                            e[i + 1] = s * r;
                            s = e[i] / r;
                            c = p / r;
                            p = c * d[i] - s * g;
                            d[i + 1] = h + s * (c * g + s * d[i]);
                            cosines[logged] = c;
                            sines[logged] = s;
                            logged++;
                        }
                        log.Rotations = logged;
                        log.First[log.Sweeps] = l;
                        log.Last[log.Sweeps] = m;
                        log.Sweeps++;
                        p = -s * s2 * c3 * el1 * e[l] / dl1;
                        e[l] = s * p;
                        d[l] = c * p;
                    } while (Math.Abs(e[l]) > eps * tst1);
                }
                d[l] += f;
                e[l] = 0d;
            }
            pending?.Wait();
            if (log.Sweeps > 0) log.Apply(W, n, parallel);
        }

        /// <summary>
        /// The order at and above which the QL rotations are applied by parallel tasks.
        /// </summary>
        private const int ParallelOrder = 128;

        /// <summary>
        /// The number of rotations logged before the log is applied, for parallel orders.
        /// </summary>
        private const int ParallelLogCapacity = 4096;

        /// <summary>
        /// A log of QL sweeps: the cosine and sine of each rotation, and the rows of each sweep.
        /// </summary>
        private sealed class RotationLog
        {
            /// <summary>
            /// Creates a log for the given number of rotations.
            /// </summary>
            /// <param name="capacity">The number of rotations, at least the largest sweep.</param>
            public RotationLog(int capacity)
            {
                Cosines = new double[capacity];
                Sines = new double[capacity];
                First = new int[capacity];
                Last = new int[capacity];
            }

            /// <summary>The cosine of each rotation.</summary>
            public double[] Cosines { get; }

            /// <summary>The sine of each rotation.</summary>
            public double[] Sines { get; }

            /// <summary>The first row l of each sweep.</summary>
            public int[] First { get; }

            /// <summary>The last row m of each sweep.</summary>
            public int[] Last { get; }

            /// <summary>The number of logged rotations.</summary>
            public int Rotations { get; set; }

            /// <summary>The number of logged sweeps.</summary>
            public int Sweeps { get; set; }

            /// <summary>
            /// Applies the logged sweeps to the rows of W, in column blocks that may run in parallel.
            /// </summary>
            /// <param name="W">The n × n row-major matrix.</param>
            /// <param name="n">The order.</param>
            /// <param name="parallel">Whether the column blocks may run in parallel.</param>
            public void Apply(double[] W, int n, bool parallel)
            {
                int sweeps = Sweeps;
                TridiagonalKernel.ForColumnBlocks(n, parallel,
                    (c0, c1) => TridiagonalKernel.ApplySweeps(W, n, Cosines, Sines, First, Last, sweeps, c0, c1));
            }
        }
    }
}
