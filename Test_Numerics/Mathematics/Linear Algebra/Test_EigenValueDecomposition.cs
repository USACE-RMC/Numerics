
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
        /// <summary>
        /// The decomposition converges to the same relative accuracy at any matrix scale, and the
        /// public input matrix keeps the values that were decomposed.
        /// </summary>
        /// <remarks>
        /// Reference eigenvalues computed with Python numpy 2.4.2 numpy.linalg.eigh on the 3x3
        /// fixture scaled by 1E-9 and by 1E+6. The assertions demand 1E-12 relative accuracy, so a
        /// convergence threshold that ignores the matrix scale is observable at either extreme.
        /// </remarks>
        [TestMethod]
        public void Test_ScaleInvariantConvergence_And_InputPreserved()
        {
            var baseA = new double[,] { { 4, 1, 1 }, { 1, 3, 0 }, { 1, 0, 2 } };
            var expected = new double[] { 1.4679111137620429, 2.6527036446661385, 4.879385241571816 };

            foreach (double scale in new double[] { 1E-9, 1.0, 1E+6 })
            {
                var A = new Matrix(3, 3);
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        A[i, j] = baseA[i, j] * scale;

                var evd = new EigenValueDecomposition(A);

                // The input matrix is preserved bit for bit
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        Assert.AreEqual(A[i, j], evd.A[i, j], $"A[{i},{j}] at scale {scale}");

                // Eigenvalues scale linearly and hold 1E-12 relative accuracy
                var w = new double[] { evd.EigenValues[0], evd.EigenValues[1], evd.EigenValues[2] };
                Array.Sort(w);
                for (int i = 0; i < 3; i++)
                    Assert.AreEqual(expected[i] * scale, w[i], Math.Abs(expected[i] * scale) * 1E-12);
            }
        }

        /// <summary>
        /// Verifies Dutilleul effective sample size is invariant to covariance scale, including
        /// scales below the former absolute degeneracy threshold.
        /// </summary>
        [TestMethod]
        public void Test_EffectiveSampleSize_IsScaleInvariant()
        {
            const double expected = 1.8d; // (1 + 2)^2 / (1^2 + 2^2)
            foreach (double scale in new double[] { 1E-9d, 1d, 1E+6d })
            {
                var matrix = new Matrix(new double[,]
                {
                    { scale, 0d },
                    { 0d, 2d * scale },
                });

                double actual = new EigenValueDecomposition(matrix).EffectiveSampleSize();
                Assert.AreEqual(expected, actual, 1E-12d, $"scale={scale:R}");
            }
        }

        /// <summary>
        /// Verifies the zero spectrum remains degenerate and tiny negative roundoff is clipped
        /// relative to the spectrum rather than by an absolute eigenvalue threshold.
        /// </summary>
        [TestMethod]
        public void Test_EffectiveSampleSize_HandlesZeroAndRelativeNegativeRoundoff()
        {
            Assert.AreEqual(0d,
                new EigenValueDecomposition(new Matrix(new double[2, 2])).EffectiveSampleSize(), 0d);

            foreach (double scale in new double[] { 1E-9d, 1d, 1E+6d })
            {
                var matrix = new Matrix(new double[,]
                {
                    { scale, 0d },
                    { 0d, -5E-11d * scale },
                });

                double actual = new EigenValueDecomposition(matrix).EffectiveSampleSize();
                Assert.AreEqual(1d, actual, 1E-12d, $"scale={scale:R}");
            }
        }

    }
}
