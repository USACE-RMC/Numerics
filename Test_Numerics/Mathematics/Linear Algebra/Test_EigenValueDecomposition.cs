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


using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class testing individual components of the Eigen Value Decomposition Method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_EigenValueDecomposition
    {

        /// <summary>
        /// 2x2 symmetric matrix with known analytic eigenvalues (3 and 1).
        /// Verifies eigenvalues and that A ≈ V D V^T and V is orthonormal.
        /// </summary>
        [TestMethod]
        public void SymEig_2x2_KnownEigenvalues()
        {
            // A = [2 1; 1 2] has eigenvalues {3, 1}
            var A = new Matrix(new double[,] { { 2.0, 1.0 }, { 1.0, 2.0 } });

            var eig = new EigenValueDecomposition(A);
            var V = eig.EigenVectors;
            var w = eig.EigenValues; // length 2

            // Check eigenvalues (order may differ): sort both and compare
            var computed = new double[] { w[0], w[1] };
            Array.Sort(computed);
            var expected = new double[] { 1.0, 3.0 };
            Array.Sort(expected);

            Assert.AreEqual(expected[0], computed[0], 1e-12);
            Assert.AreEqual(expected[1], computed[1], 1e-12);

            // Build diagonal D
            var D = new Matrix(2);
            D[0, 0] = w[0]; D[1, 1] = w[1];

            // Check reconstruction: A ≈ V * D * V^T
            var Arec = V * D * Matrix.Transpose(V);
            AssertMatrixAlmostEqual(A, Arec);

            // Check orthonormality: V^T V ≈ I
            var I = Matrix.Identity(2);
            var VtV = Matrix.Transpose(V) * V;
            AssertMatrixAlmostEqual(I, VtV);
        }

        /// <summary>
        /// 3x3 symmetric matrix. Verifies A ≈ V D V^T and V orthonormality.
        /// </summary>
        [TestMethod]
        public void SymEig_3x3_ReconstructAndOrthonormal()
        {
            // Symmetric test matrix
            var A = new Matrix(new double[,]
            {
                { 4, 1, 1 },
                { 1, 3, 0 },
                { 1, 0, 2 }
            });

            var eig = new EigenValueDecomposition(A);
            var V = eig.EigenVectors;
            var w = eig.EigenValues;

            var D = new Matrix(3);
            for (int i = 0; i < 3; i++) 
                D[i, i] = w[i];

            // A ≈ V D V^T
            var Arec = V * D * Matrix.Transpose(V);
            AssertMatrixAlmostEqual(A, Arec);

            // V^T V ≈ I
            var I = Matrix.Identity(3);
            var VtV = Matrix.Transpose(V) * V;
            AssertMatrixAlmostEqual(I, VtV);
        }

        /// <summary>
        /// 5x5 symmetric matrix (constructed deterministically). Verifies reconstruction and orthonormality.
        /// </summary>
        [TestMethod]
        public void SymEig_5x5_SymmetricRandomLike()
        {
            // Build a deterministic symmetric 5x5
            var baseM = new double[,]
            {
                {  2, -1,  0,  2, -3 },
                {  1,  4,  5, -2,  0 },
                { -2,  3,  6,  1,  4 },
                {  0, -1,  2,  3, -2 },
                {  1,  2, -4,  0,  5 }
            };
            // Symmetrize: A = (M + M^T)/2
            var A = new Matrix(5);
            for (int i = 0; i < 5; i++)
                for (int j = 0; j < 5; j++)
                    A[i, j] = 0.5 * (baseM[i, j] + baseM[j, i]);

            var eig = new EigenValueDecomposition(A);
            var V = eig.EigenVectors;
            var w = eig.EigenValues;

            var D = new Matrix(5);
            for (int i = 0; i < 5; i++) 
                D[i, i] = w[i];

            // A ≈ V D V^T
            var Arec = V * D * Matrix.Transpose(V);
            AssertMatrixAlmostEqual(A, Arec);

            // V^T V ≈ I
            var I = Matrix.Identity(5);
            var VtV = Matrix.Transpose(V) * V;
            AssertMatrixAlmostEqual(I, VtV);
        }

        /// <summary>
        /// Degenerate spectrum: A = 2 * I (3x3). All eigenvalues are 2 (multiplicity 3).
        /// V may be any orthonormal basis. Checks values, reconstruction, orthonormality, and residuals.
        /// </summary>
        [TestMethod]
        public void SymEig_3x3_RepeatedEigenvalues_AllTwos()
        {
            var A = new Matrix(3);
            for (int i = 0; i < 3; i++) A[i, i] = 2.0;

            var eig = new EigenValueDecomposition(A);
            var V = eig.EigenVectors;
            var w = eig.EigenValues;

            // All eigenvalues = 2 (order arbitrary)
            for (int i = 0; i < 3; i++)
                Assert.AreEqual(2.0, w[i], 1e-12);

            // Orthonormal V
            AssertOrthonormal(V);

            // Reconstruction A ≈ V D V^T
            var D = new Matrix(3);
            for (int i = 0; i < 3; i++) D[i, i] = w[i];
            var Arec = V * D * Matrix.Transpose(V);
            AssertMatrixAlmostEqual(A, Arec);

            // Max eigen residual
            var maxRes = MaxEigenResidual(A, V, w);
            Assert.IsLessThan(1e-12, maxRes, $"Max eigen residual too large: {maxRes}");
        }

        /// <summary>
        /// Structured Toeplitz tridiagonal (n=8): diag=2, offdiag=-1. 
        /// Eigenvalues are λ_k = 2 - 2 cos(kπ/(n+1)), k = 1..n. 
        /// Checks spectrum, reconstruction, orthonormality, and residuals.
        /// </summary>
        [TestMethod]
        public void SymEig_8x8_TridiagonalToeplitz_KnownSpectrum()
        {
            int n = 8;
            var A = new Matrix(n);
            for (int i = 0; i < n; i++)
            {
                A[i, i] = 2.0;
                if (i + 1 < n) { A[i, i + 1] = -1.0; A[i + 1, i] = -1.0; }
            }

            // Expected eigenvalues (sorted ascending)
            var expected = new double[n];
            for (int k = 1; k <= n; k++)
                expected[k - 1] = 2.0 - 2.0 * Math.Cos(k * Math.PI / (n + 1));

            var eig = new EigenValueDecomposition(A);
            var V = eig.EigenVectors;
            var w = eig.EigenValues;

            // Sort both
            var computed = (double[])w.ToArray().Clone();
            Array.Sort(expected);
            Array.Sort(computed);

            for (int i = 0; i < n; i++)
                Assert.AreEqual(expected[i], computed[i], 1e-8, $"Eigenvalue mismatch at i={i}");

            // Orthonormal V
            AssertOrthonormal(V, 1e-9);

            // Reconstruction A ≈ V D V^T
            var D = new Matrix(n);
            for (int i = 0; i < n; i++) D[i, i] = w[i];
            var Arec = V * D * Matrix.Transpose(V);
            AssertMatrixAlmostEqual(A, Arec, 1e-8);

            // Max eigen residual
            var maxRes = MaxEigenResidual(A, V, w);
            Assert.IsLessThan(1e-8, maxRes, $"Max eigen residual too large: {maxRes}");
        }

        /// <summary>
        /// Nearly-diagonal symmetric with tiny off-diagonal couplings. 
        /// Stresses convergence when off-diagonals are very small.
        /// Checks reconstruction, orthonormality, and residuals.
        /// </summary>
        [TestMethod]
        public void SymEig_5x5_NearlyDiagonal_SmallCoupling()
        {
            int n = 5;
            var A = new Matrix(n);
            // Start diagonal with separated values
            for (int i = 0; i < n; i++) A[i, i] = i + 1; // 1..5
            // Add very small symmetric couplings
            double eps = 1e-9;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    A[i, j] = A[j, i] = ((i + j) % 2 == 0 ? eps : -eps);

            var eig = new EigenValueDecomposition(A);
            var V = eig.EigenVectors;
            var w = eig.EigenValues;

            // Orthonormal V
            AssertOrthonormal(V, 1e-9);

            // Reconstruction A ≈ V D V^T
            var D = new Matrix(n);
            for (int i = 0; i < n; i++) D[i, i] = w[i];
            var Arec = V * D * Matrix.Transpose(V);
            AssertMatrixAlmostEqual(A, Arec, 1e-8);

            // Max eigen residual
            var maxRes = MaxEigenResidual(A, V, w);
            Assert.IsLessThan(1e-8, maxRes, $"Max eigen residual too large: {maxRes}");
        }

        // ---------- Helpers ----------

        private static void AssertOrthonormal(Matrix V, double tol = 1E-10)
        {
            var I = Matrix.Identity(V.NumberOfColumns);
            var VtV = Matrix.Transpose(V) * V;
            AssertMatrixAlmostEqual(I, VtV, tol);
        }

        private static double MaxEigenResidual(Matrix A, Matrix V, Vector w)
        {
            int n = A.NumberOfRows;
            double max = 0.0;
            for (int j = 0; j < n; j++)
            {
                // r_j = A v_j - λ_j v_j
                var vj = new Vector(V.Column(j));
                var Av = A * vj;
                double lambda = w[j];
                double norm = 0.0;
                for (int i = 0; i < n; i++)
                {
                    double r = Av[i] - lambda * vj[i];
                    norm += r * r;
                }
                norm = Math.Sqrt(norm);
                if (norm > max) max = norm;
            }
            return max;
        }

        private static void AssertMatrixAlmostEqual(Matrix expected, Matrix actual, double tol = 1E-10)
        {
            Assert.AreEqual(expected.NumberOfRows, actual.NumberOfRows);
            Assert.AreEqual(expected.NumberOfColumns, actual.NumberOfColumns);
            for (int i = 0; i < expected.NumberOfRows; i++)
                for (int j = 0; j < expected.NumberOfColumns; j++)
                    Assert.AreEqual(expected[i, j], actual[i, j], tol);
        }
    }
}
