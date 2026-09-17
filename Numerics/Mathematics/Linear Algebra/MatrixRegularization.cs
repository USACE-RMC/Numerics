namespace Numerics.Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class for performing Matrix regularization.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list>
    /// </para>
    /// </remarks>
    public static class MatrixRegularization
    {
        /// <summary>
        /// Eigen-regularizes a symmetric matrix Vb and returns a PSD matrix suitable for Cholesky.
        /// Vb is (approximately) symmetrized, eigen-decomposed Vb = Q diag(λ) Qᵀ, eigenvalues are floored/capped,
        /// and the matrix is reconstructed as Q diag(λ_reg) Qᵀ.
        /// </summary>
        /// <param name="Vb">Symmetric matrix.</param>
        /// <param name="eps">Floor expressed as eps * (trace(Vb)/p), falling back to eps if trace ≤ 0.</param>
        /// <param name="capMult">Cap is capMult * median(λ).</param>
        public static Matrix Regularize(Matrix Vb, double eps = 1e-6, double capMult = 50.0)
        {
            if (!Vb.IsSquare) throw new ArgumentException("Vb must be square.", nameof(Vb));
            int p = Vb.NumberOfRows;

            // Ensure exact symmetry: A := (A + Aᵀ)/2 (nice to have before Jacobi)
            var VbSym = Symmetrize(Vb);

            // Eigen-decomposition (symmetric): Vb = Q diag(λ) Qᵀ
            var eig = new EigenValueDecomposition(VbSym);
            var Q = eig.EigenVectors;
            var w = eig.EigenValues; // length p

            // Compute trace and robust floor
            double trace = 0.0; for (int i = 0; i < p; i++) trace += w[i];
            double floor = eps * (trace > 0.0 ? trace / p : 1.0);

            // Median of eigenvalues (copy -> sort)
            double median = MedianFromVector(w);
            double cap = capMult * median;

            // Floor and cap eigenvalues
            var D = new Matrix(p);
            for (int i = 0; i < p; i++)
            {
                double li = w[i];
                if (li < floor) li = floor;
                if (li > cap) li = cap;
                D[i, i] = li;
            }

            // Recompose: Q * D * Qᵀ
            return Q * D * Matrix.Transpose(Q);
        }

        /// <summary>
        /// Eigen-regularizes a symmetric matrix Vb and returns a PSD matrix suitable for Cholesky.
        /// Vb is (approximately) symmetrized, eigen-decomposed Vb = Q diag(λ) Qᵀ, eigenvalues are floored/capped,
        /// and the matrix is reconstructed as Q diag(λ_reg) Qᵀ.
        /// </summary>
        /// <param name="Vb">Symmetric matrix.</param>
        /// <param name="eps">Floor expressed as eps * (trace(Vb)/p), falling back to eps if trace ≤ 0.</param>
        /// <param name="capMult">Cap is capMult * median(λ).</param>
        public static double[,] Regularize(double[,] Vb, double eps = 1e-6, double capMult = 50.0)
        {
            var A = new Matrix(Vb);
            var R = Regularize(A, eps, capMult);
            return R.ToArray();
        }

        /// <summary>
        /// Symmetrizes the matrix.
        /// </summary>
        /// <param name="A">The matrix to evaluate.</param>
        private static Matrix Symmetrize(Matrix A)
        {
            int n = A.NumberOfRows, m = A.NumberOfColumns;
            var S = new Matrix(n, m);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    S[i, j] = 0.5 * (A[i, j] + A[j, i]);
                }
            }
            return S;
        }

        /// <summary>
        /// Returns the median value from the vector.
        /// </summary>
        /// <param name="v">The vector to evaluate.</param>
        private static double MedianFromVector(Vector v)
        {
            int n = v.Length;
            var arr = new double[n];
            for (int i = 0; i < n; i++) arr[i] = v[i];
            Array.Sort(arr);
            if ((n & 1) == 1) return arr[n / 2];
            return 0.5 * (arr[n / 2 - 1] + arr[n / 2]);
        }

        /// <summary>
        /// Determines whether the Cholesky decomposition accepts a matrix as positive definite.
        /// </summary>
        /// <param name="matrix">The symmetric matrix to test.</param>
        /// <returns><see langword="true"/> when the decomposition succeeds; otherwise, <see langword="false"/>.</returns>
        private static bool CholeskyAccepts(Matrix matrix)
        {
            return CholeskyDecomposition.TryFactorize(matrix, CholeskyDecomposition.DefaultRelativeTolerance(matrix.NumberOfRows), out _, out _, out _);
        }

        /// <summary>
        /// Makes the matrix symmetric and, when necessary, adds a ridge until Cholesky accepts it as
        /// positive definite.
        /// </summary>
        /// <param name="M">The finite square matrix to adjust.</param>
        /// <returns>A finite symmetric matrix accepted by the default Cholesky pivot test.</returns>
        /// <exception cref="ArgumentNullException">Thrown when M is null.</exception>
        /// <exception cref="ArgumentException">Thrown when M is not square or contains non-finite entries.</exception>
        /// <exception cref="InvalidOperationException">
        /// Thrown when symmetrization or the ridge scale is not representable, or no finite candidate
        /// on the ridge ladder passes the Cholesky test before a ridge or diagonal overflows.
        /// </exception>
        /// <remarks>
        /// The symmetric input is returned without a ridge when its Cholesky decomposition succeeds.
        /// Otherwise the base ridge remains 1E-10 times the mean diagonal for a positive trace, or
        /// 1E-10 for a non-positive trace. The original first eight candidates retain their arithmetic;
        /// subsequent candidates multiply the previous ridge by ten. Every returned candidate is checked.
        /// The maximum permitted ridge is the largest finite value on this decade ladder for which all
        /// adjusted diagonal entries remain finite. Exhaustion throws rather than returning an unchecked
        /// or smaller ridge. No pivot tolerance or eigenvalue floor is changed.
        /// </remarks>
        public static Matrix MakeSymmetricPositiveDefinite(Matrix M)
        {
            if (M == null) throw new ArgumentNullException(nameof(M));
            if (M.NumberOfRows != M.NumberOfColumns)
                throw new ArgumentException("The matrix must be square.", nameof(M));
            for (int i = 0; i < M.NumberOfRows; i++)
            {
                for (int j = 0; j < M.NumberOfColumns; j++)
                {
                    if (!Tools.IsFinite(M[i, j]))
                        throw new ArgumentException("The matrix must contain only finite entries.", nameof(M));
                }
            }

            var S = 0.5 * (M + M.Transpose());
            for (int i = 0; i < S.NumberOfRows; i++)
            {
                for (int j = 0; j < S.NumberOfColumns; j++)
                {
                    if (!Tools.IsFinite(S[i, j]))
                        throw new InvalidOperationException("Matrix symmetrization produced a non-finite entry.");
                }
            }
            if (CholeskyAccepts(S)) return S;

            double tr = 0.0;
            for (int i = 0; i < S.NumberOfRows; i++) tr += S[i, i];
            double baseRidge = (tr > 0 ? 1e-10 * tr / S.NumberOfRows : 1e-10);
            if (!Tools.IsFinite(tr) || !Tools.IsFinite(baseRidge) || baseRidge <= 0d)
                throw new InvalidOperationException("The trace-scaled matrix ridge is not a positive finite number.");

            double ridge = baseRidge;
            for (int k = 0; Tools.IsFinite(ridge); k++)
            {
                var T = S.Clone();
                for (int i = 0; i < T.NumberOfRows; i++)
                {
                    T[i, i] += ridge;
                    if (!Tools.IsFinite(T[i, i]))
                        throw new InvalidOperationException("Matrix ridge escalation overflowed a diagonal before positive definiteness was achieved.");
                }
                if (CholeskyAccepts(T)) return T;

                // Preserve established candidates exactly, then continue monotonically. Multiplying the
                // ridge itself avoids overflowing an unscaled power of ten for very small trace scales.
                ridge = k < 7 ? baseRidge * Math.Pow(10.0, k + 1) : ridge * 10.0;
            }
            throw new InvalidOperationException("Matrix ridge escalation exhausted finite candidates before positive definiteness was achieved.");
        }
    }
}
