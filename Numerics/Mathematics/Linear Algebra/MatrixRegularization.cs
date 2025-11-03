/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
        /// Makes the matrix symmetric and positive definite.
        /// </summary>
        /// <param name="M">The matrix to adjust.</param>
        /// <returns>A symmetric and positive definite matrix.</returns>
        public static Matrix MakeSymmetricPositiveDefinite(Matrix M)
        {
            // Symmetrize
            var S = 0.5 * (M + M.Transpose());
            // Tiny trace-scaled ridge
            double tr = 0.0;
            for (int i = 0; i < S.NumberOfRows; i++) tr += S[i, i];
            double baseRidge = (tr > 0 ? 1e-10 * tr / S.NumberOfRows : 1e-10);

            // Try increasing ridge until Cholesky succeeds
            for (int k = 0; k < 8; k++)
            {
                var T = S.Clone();
                double ridge = baseRidge * Math.Pow(10.0, k);
                for (int i = 0; i < T.NumberOfRows; i++) T[i, i] += ridge;
                try { var _ = new CholeskyDecomposition(T); return T; } catch { /* retry bigger ridge */ }
            }
            // Last resort: add a biggish ridge
            var U = S.Clone();
            double big = (tr > 0 ? 1e-4 * tr / S.NumberOfRows : 1e-4);
            for (int i = 0; i < U.NumberOfRows; i++) U[i, i] += big;
            return U;
        }   
    }
}
