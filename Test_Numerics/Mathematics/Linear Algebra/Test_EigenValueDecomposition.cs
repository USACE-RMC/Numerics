
using System;
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
    public partial class Test_EigenValueDecomposition
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

        /// <summary>
        /// IsConverged reports whether the rotation limit was reached, without changing any result.
        /// </summary>
        [TestMethod]
        public void Test_IsConvergedReportsRotationLimit()
        {
            // Small and trivial matrices converge.
            Assert.IsTrue(new EigenValueDecomposition(new Matrix(new double[,] { { 2d, 1d }, { 1d, 3d } })).IsConverged);
            Assert.IsTrue(new EigenValueDecomposition(new Matrix(new double[,] { { 5d } })).IsConverged);
            Assert.IsTrue(new EigenValueDecomposition(new Matrix(0)).IsConverged);
            Assert.IsTrue(new EigenValueDecomposition(RandomGram(10, 12345)).IsConverged);

            // A fixed 2000-rotation limit stops n = 60 early (max |AV - VD| of 7.2E-2); the limit scales with n^2,
            // so the decomposition converges and the residual is at roundoff.
            var S = RandomGram(60, 12345);
            var eig = new EigenValueDecomposition(S);
            Assert.IsTrue(eig.IsConverged);
            var residual = S * eig.EigenVectors - eig.EigenVectors * Matrix.Diagonal(eig.EigenValues);
            double max = 0d;
            for (int i = 0; i < 60; i++)
                for (int j = 0; j < 60; j++)
                    max = System.Math.Max(max, System.Math.Abs(residual[i, j]));
            Assert.IsLessThan(1E-10, max);

            // Non-finite entries are rejected rather than returning eigenvalues such as (0, 2) for NaN.
            Assert.ThrowsExactly<System.ArgumentException>(() => new EigenValueDecomposition(new Matrix(new double[,] { { 1d, double.NaN }, { double.NaN, 2d } })));
            Assert.ThrowsExactly<System.ArgumentException>(() => new EigenValueDecomposition(new Matrix(new double[,] { { double.PositiveInfinity, 0d }, { 0d, 2d } })));
        }

        /// <summary>
        /// The 50 x 50 tridiagonal Toeplitz matrix with 2 on the diagonal and -1 beside it has the analytic
        /// eigenvalues 2 - 2 cos(k pi / 51), k = 1..50. It needs more than 2000 rotations, so it also pins
        /// the scaled rotation limit.
        /// </summary>
        [TestMethod]
        public void Test_LargeMatrixMatchesAnalyticEigenvalues()
        {
            const int n = 50;
            var A = new Matrix(n);
            for (int i = 0; i < n; i++)
            {
                A[i, i] = 2d;
                if (i > 0) A[i, i - 1] = A[i - 1, i] = -1d;
            }
            var eig = new EigenValueDecomposition(A);
            Assert.IsTrue(eig.IsConverged);
            var computed = eig.EigenValues.ToArray();
            System.Array.Sort(computed);
            for (int k = 1; k <= n; k++)
            {
                double exact = 2d - 2d * System.Math.Cos(k * System.Math.PI / (n + 1));
                // The stopping threshold is 1E-12 of the largest element (2), so eigenvalues agree to about 1E-11.
                Assert.AreEqual(exact, computed[k - 1], 1E-10);
            }
        }

    
        /// <summary>
        /// A 120 × 120 Gram matrix converges within the rotation limit; the eigenvectors are orthonormal, A·V = V·Λ to
        /// rounding, and the eigenvalues match an independent Householder/QL decomposition.
        /// </summary>
        [TestMethod]
        public void Test_CyclicJacobi_LargeGram()
        {
            int n = 120;
            var rng = new System.Random(9);
            var F = new Matrix(400, n);
            for (int i = 0; i < 400; i++)
                for (int j = 0; j < n; j++) F[i, j] = rng.NextDouble() - 0.5 + 0.01 * j;
            var A = F.TransposeMultiply(F);
            var eigen = new EigenValueDecomposition(A);
            Assert.IsTrue(eigen.IsConverged);
            var V = eigen.EigenVectors;
            double scale = 0;
            for (int i = 0; i < n; i++) scale = System.Math.Max(scale, System.Math.Abs(eigen.EigenValues[i]));
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double dot = 0, av = 0;
                    for (int k = 0; k < n; k++)
                    {
                        dot += V[k, i] * V[k, j];
                        av += A[j, k] * V[k, i];
                    }
                    Assert.AreEqual(i == j ? 1.0 : 0.0, dot, 1E-12);
                    Assert.AreEqual(eigen.EigenValues[i] * V[j, i], av, 1E-10 * scale);
                }
            }

            var flat = new double[n * n];
            for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) flat[i * n + j] = A[i, j];
            EigenValueDecomposition.DecomposeTridiagonal(flat, n, out var reference, out _);
            var sorted = eigen.EigenValues.ToArray();
            System.Array.Sort(sorted);
            System.Array.Reverse(sorted);
            for (int i = 0; i < n; i++) Assert.AreEqual(reference[i], sorted[i], 1E-12 * scale);
        }

        /// <summary>
        /// Measures the accuracy figures stated in the documentation on random Gram matrices and on matrices
        /// Q Λ Qᵀ with a known spectrum, Q a product of three Householder reflections: Λ evenly spaced in [1, n], and
        /// Λ graded geometrically from 1 down to 1E-12. For n = 10 to 300 the measured eigenpair residuals
        /// ‖Av − λv‖₂ / ‖A‖_F are 2.2E-13 to 2.7E-12 (asserted below 1E-11), the eigenvectors are orthonormal to
        /// 1.5E-14 (asserted below 5E-14), and the eigenvalues agree with Λ to 9.3E-15 of the largest eigenvalue for
        /// the evenly spaced spectrum (asserted below 2E-14) and 1.2E-12 for the graded one (asserted below 2E-12).
        /// This test checks n = 10, 50 and 120; the sweep that adds n = 300 runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep CyclicJacobi</c>.
        /// </summary>
        [TestMethod]
        public void Test_CyclicJacobi_DocumentedAccuracy()
        {
            CyclicJacobiDocumentedAccuracy(10, 50, 120);
        }

        /// <summary>
        /// The QL method returns ascending eigenvalues equal to numpy.linalg.eigh and 50-digit mpmath.eigsy, with unit
        /// eigenvector columns that satisfy A v = λ v.
        /// </summary>
        [TestMethod]
        public void Test_QL_KnownValues()
        {
            var A = new Matrix(new double[,] { { 4.0, 1.0, -2.0, 0.5 }, { 1.0, 3.0, 0.0, 1.5 }, { -2.0, 0.0, 5.0, -1.0 }, { 0.5, 1.5, -1.0, 2.0 } });
            double[] exact = { 0.641278407914022, 2.4510237937334756, 3.858181525316084, 7.049516273036418 };
            var evd = new EigenValueDecomposition(A, EigenSolverMethod.QL);
            Assert.IsTrue(evd.IsConverged);
            for (int i = 0; i < 4; i++)
            {
                Assert.AreEqual(exact[i], evd.EigenValues[i], 4e-15 * 7.05, $"eigenvalue {i}");
                double norm = 0;
                for (int k = 0; k < 4; k++)
                {
                    double av = 0;
                    for (int m = 0; m < 4; m++) av += A[k, m] * evd.EigenVectors[m, i];
                    Assert.AreEqual(evd.EigenValues[i] * evd.EigenVectors[k, i], av, 1e-14, $"residual {i},{k}");
                    norm += evd.EigenVectors[k, i] * evd.EigenVectors[k, i];
                }
                Assert.AreEqual(1d, norm, 1e-15);
            }
            Assert.AreEqual(4.0, evd.A[0, 0]);
        }

        /// <summary>
        /// On random Gram matrices the QL eigenvalues are ascending and agree with the Jacobi eigenvalues, the eigenpair
        /// residual ‖A v − λ v‖ is a small multiple of ε‖A‖, and the eigenvectors are orthonormal; the default constructor
        /// and the Jacobi method give identical bits; a 0 × 0 matrix and an undefined method are handled.
        /// </summary>
        [TestMethod]
        public void Test_QL_RandomMatricesAndJacobiIdentity()
        {
            foreach (int n in new[] { 1, 2, 5, 17, 60 })
            {
                var rng = new Random(n);
                var B = new double[n, n];
                for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) B[i, j] = rng.NextDouble() - 0.5;
                var a = new double[n, n];
                for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { double s = 0; for (int k = 0; k < n; k++) s += B[i, k] * B[j, k]; a[i, j] = s; }
                var A = new Matrix(a);
                var ql = new EigenValueDecomposition(A, EigenSolverMethod.QL);
                var jacobi = new EigenValueDecomposition(A, EigenSolverMethod.Jacobi);
                var byDefault = new EigenValueDecomposition(A);

                double scale = 0;
                for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) scale = Math.Max(scale, Math.Abs(a[i, j]));
                var sortedJacobi = jacobi.EigenValues.ToArray();
                Array.Sort(sortedJacobi);
                for (int i = 0; i < n; i++)
                {
                    if (i > 0) Assert.IsLessThanOrEqualTo(ql.EigenValues[i], ql.EigenValues[i - 1]);
                    Assert.AreEqual(sortedJacobi[i], ql.EigenValues[i], 1e-12 * scale);
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(byDefault.EigenValues[i]), BitConverter.DoubleToInt64Bits(jacobi.EigenValues[i]));
                    for (int k = 0; k < n; k++)
                    {
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(byDefault.EigenVectors[k, i]), BitConverter.DoubleToInt64Bits(jacobi.EigenVectors[k, i]));
                        double av = 0;
                        for (int m = 0; m < n; m++) av += a[k, m] * ql.EigenVectors[m, i];
                        Assert.AreEqual(ql.EigenValues[i] * ql.EigenVectors[k, i], av, 50 * n * 2.2e-16 * scale, $"n = {n} residual");
                    }
                    for (int j = 0; j <= i; j++)
                    {
                        double dot = 0;
                        for (int k = 0; k < n; k++) dot += ql.EigenVectors[k, i] * ql.EigenVectors[k, j];
                        Assert.AreEqual(i == j ? 1d : 0d, dot, 1e-13, $"n = {n} orthonormality");
                    }
                }
            }

            var empty = new EigenValueDecomposition(new Matrix(0), EigenSolverMethod.QL);
            Assert.AreEqual(0, empty.EigenValues.Length);
            Assert.IsTrue(empty.IsConverged);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new EigenValueDecomposition(Matrix.Identity(2), (EigenSolverMethod)7));
        }

        /// <summary>
        /// On exactly symmetric matrices the decomposition is bitwise the column-indexed tred2 followed by a transpose,
        /// tql2 and the decreasing sort, evaluated by the sequential reference below.
        /// </summary>
        [TestMethod]
        public void Test_DecomposeTridiagonal_MatchesColumnIndexedReferenceBitwise()
        {
            int cases = 0;
            foreach (int n in new[] { 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 33, 64, 100 })
            {
                for (int seed = 0; seed < 4; seed++)
                {
                    var rng = new Random(1000 * n + seed);
                    var B = new double[n * n];
                    for (int i = 0; i < n * n; i++) B[i] = rng.NextDouble() - 0.5;
                    var A = new double[n * n];
                    for (int i = 0; i < n; i++)
                        for (int j = 0; j <= i; j++)
                        {
                            double s = 0;
                            for (int k = 0; k < n; k++) s += B[i * n + k] * B[j * n + k];
                            if (seed == 3 && (i + j) % 4 == 1) s = 0;
                            A[i * n + j] = A[j * n + i] = s;
                        }

                    var V = (double[])A.Clone();
                    var d = new double[n];
                    var e = new double[n];
                    ReferenceTridiagonalize(V, d, e, n);
                    for (int i = 0; i < n; i++)
                        for (int j = i + 1; j < n; j++)
                            (V[i * n + j], V[j * n + i]) = (V[j * n + i], V[i * n + j]);
                    ReferenceDiagonalize(V, d, e, n);
                    var order = new int[n];
                    for (int i = 0; i < n; i++) order[i] = i;
                    Array.Sort(order, (x, y) => d[y] != d[x] ? d[y].CompareTo(d[x]) : x.CompareTo(y));

                    EigenValueDecomposition.DecomposeTridiagonal((double[])A.Clone(), n, out var values, out var vectors);
                    for (int i = 0; i < n; i++)
                    {
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(d[order[i]]), BitConverter.DoubleToInt64Bits(values[i]), $"n = {n} seed {seed} value {i}");
                        for (int k = 0; k < n; k++)
                            Assert.AreEqual(BitConverter.DoubleToInt64Bits(V[order[i] * n + k]), BitConverter.DoubleToInt64Bits(vectors[i * n + k]), $"n = {n} seed {seed} vector {i}");
                    }
                    cases++;
                }
            }
            Assert.AreEqual(60, cases);
        }

        /// <summary>
        /// The SIMD and batched-rotation paths stay bitwise equal to the sequential reference at orders that leave every
        /// kind of vector tail and flush the rotation log several times, on random Gram matrices, on a matrix with zero
        /// columns (the scale = 0 branch of tred2), on a diagonal matrix with repeated entries, and on an input symmetric
        /// only to rounding, which is decomposed as the matrix mirrored from its upper triangle.
        /// </summary>
        [TestMethod]
        public void Test_DecomposeTridiagonal_LargeAndStructuredMatchReferenceBitwise()
        {
            var inputs = new System.Collections.Generic.List<(string Name, int N, double[] A)>();
            foreach (int n in new[] { 63, 129, 211 })
            {
                var rng = new Random(n);
                var B = new double[n * n];
                for (int i = 0; i < n * n; i++) B[i] = rng.NextDouble() - 0.5;
                var A = new double[n * n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j <= i; j++)
                    {
                        double s = 0;
                        for (int k = 0; k < n; k++) s += B[i * n + k] * B[j * n + k];
                        A[i * n + j] = A[j * n + i] = s;
                    }
                inputs.Add(($"gram {n}", n, A));

                var Z = (double[])A.Clone();
                foreach (int z in new[] { 0, n / 3, n - 2 })
                    for (int k = 0; k < n; k++) Z[z * n + k] = Z[k * n + z] = 0d;
                inputs.Add(($"zero columns {n}", n, Z));
            }
            {
                int n = 37;
                var D = new double[n * n];
                for (int i = 0; i < n; i++) D[i * n + i] = 1d + (i % 3);
                inputs.Add(("diagonal", n, D));
            }

            foreach (var (name, n, A) in inputs)
            {
                AssertMatchesReference(A, A, n, name);

                // Perturb the lower triangle by one ulp; the decomposition reads the upper triangle only.
                var perturbed = (double[])A.Clone();
                for (int i = 1; i < n; i++)
                    for (int j = 0; j < i; j += 3) perturbed[i * n + j] = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(perturbed[i * n + j]) + 1);
                AssertMatchesReference(perturbed, A, n, name + " perturbed");
            }
        }

        /// <summary>
        /// Asserts that DecomposeTridiagonal of <paramref name="input"/> is bitwise the sequential reference applied to the
        /// exactly symmetric <paramref name="symmetric"/>.
        /// </summary>
        private static void AssertMatchesReference(double[] input, double[] symmetric, int n, string name)
        {
            var V = (double[])symmetric.Clone();
            var d = new double[n];
            var e = new double[n];
            ReferenceTridiagonalize(V, d, e, n);
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    (V[i * n + j], V[j * n + i]) = (V[j * n + i], V[i * n + j]);
            ReferenceDiagonalize(V, d, e, n);
            var order = new int[n];
            for (int i = 0; i < n; i++) order[i] = i;
            Array.Sort(order, (x, y) => d[y] != d[x] ? d[y].CompareTo(d[x]) : x.CompareTo(y));

            EigenValueDecomposition.DecomposeTridiagonal((double[])input.Clone(), n, out var values, out var vectors);
            for (int i = 0; i < n; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(d[order[i]]), BitConverter.DoubleToInt64Bits(values[i]), $"{name} value {i}");
                for (int k = 0; k < n; k++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(V[order[i] * n + k]), BitConverter.DoubleToInt64Bits(vectors[i * n + k]), $"{name} vector {i}");
            }
        }

        /// <summary>The column-indexed Householder tridiagonalization (tred2), the sequential reference.</summary>
        private static void ReferenceTridiagonalize(double[] V, double[] d, double[] e, int n)
        {
            for (int j = 0; j < n; j++) d[j] = V[(n - 1) * n + j];
            for (int i = n - 1; i > 0; i--)
            {
                double scale = 0d, h = 0d;
                for (int k = 0; k < i; k++) scale += Math.Abs(d[k]);
                if (scale == 0d)
                {
                    e[i] = d[i - 1];
                    for (int j = 0; j < i; j++)
                    {
                        d[j] = V[(i - 1) * n + j];
                        V[i * n + j] = 0d;
                        V[j * n + i] = 0d;
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
                    for (int j = 0; j < i; j++) e[j] = 0d;
                    for (int j = 0; j < i; j++)
                    {
                        f = d[j];
                        V[j * n + i] = f;
                        g = e[j] + V[j * n + j] * f;
                        for (int k = j + 1; k <= i - 1; k++)
                        {
                            g += V[k * n + j] * d[k];
                            e[k] += V[k * n + j] * f;
                        }
                        e[j] = g;
                    }
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
                        f = d[j];
                        g = e[j];
                        for (int k = j; k <= i - 1; k++) V[k * n + j] -= f * e[k] + g * d[k];
                        d[j] = V[(i - 1) * n + j];
                        V[i * n + j] = 0d;
                    }
                }
                d[i] = h;
            }

            // Accumulate the transformations.
            for (int i = 0; i < n - 1; i++)
            {
                V[(n - 1) * n + i] = V[i * n + i];
                V[i * n + i] = 1d;
                double h = d[i + 1];
                if (h != 0d)
                {
                    for (int k = 0; k <= i; k++) d[k] = V[k * n + i + 1] / h;
                    for (int j = 0; j <= i; j++)
                    {
                        double g = 0d;
                        for (int k = 0; k <= i; k++) g += V[k * n + i + 1] * V[k * n + j];
                        for (int k = 0; k <= i; k++) V[k * n + j] -= g * d[k];
                    }
                }
                for (int k = 0; k <= i; k++) V[k * n + i + 1] = 0d;
            }
            for (int j = 0; j < n; j++)
            {
                d[j] = V[(n - 1) * n + j];
                V[(n - 1) * n + j] = 0d;
            }
            V[(n - 1) * n + n - 1] = 1d;
            e[0] = 0d;
        }

        /// <summary>Implicit QL iterations (tql2) on rows of W, the sequential reference.</summary>
        private static void ReferenceDiagonalize(double[] W, double[] d, double[] e, int n)
        {
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
                        if (++iteration > 60) throw new ArithmeticException("The QL iteration did not converge.");
                        double g = d[l];
                        double p = (d[l + 1] - g) / (2d * e[l]);
                        double r = Numerics.Tools.Hypot(p, 1d);
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
                        for (int i = m - 1; i >= l; i--)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * e[i];
                            h = c * p;
                            r = Numerics.Tools.Hypot(p, e[i]);
                            e[i + 1] = s * r;
                            s = e[i] / r;
                            c = p / r;
                            p = c * d[i] - s * g;
                            d[i + 1] = h + s * (c * g + s * d[i]);
                            ReferenceRotate(W, i * n, (i + 1) * n, n, c, s);
                        }
                        p = -s * s2 * c3 * el1 * e[l] / dl1;
                        e[l] = s * p;
                        d[l] = c * p;
                    } while (Math.Abs(e[l]) > eps * tst1);
                }
                d[l] += f;
                e[l] = 0d;
            }
        }

        /// <summary>The plane rotation of rows a and b.</summary>
        private static void ReferenceRotate(double[] W, int a, int b, int n, double c, double s)
        {
            int w = System.Numerics.Vector<double>.Count;
            var vc = new System.Numerics.Vector<double>(c);
            var vs = new System.Numerics.Vector<double>(s);
            int k = 0;
            for (; k + w <= n; k += w)
            {
                var x = new System.Numerics.Vector<double>(W, a + k);
                var y = new System.Numerics.Vector<double>(W, b + k);
                (vs * x + vc * y).CopyTo(W, b + k);
                (vc * x - vs * y).CopyTo(W, a + k);
            }
            for (; k < n; k++)
            {
                double x = W[a + k], y = W[b + k];
                W[b + k] = s * x + c * y;
                W[a + k] = c * x - s * y;
            }
        }

    }
}
