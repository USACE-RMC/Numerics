using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class characterizing <see cref="MatrixRegularization"/>, whose ridge-escalation loop is driven by
    /// whether <see cref="CholeskyDecomposition"/> accepts a candidate matrix.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list>
    /// </para>
    /// <para>
    /// <see cref="MatrixRegularization.MakeSymmetricPositiveDefinite"/> symmetrizes its input and first tests
    /// that un-ridged candidate. It returns the candidate unchanged when Cholesky accepts it. Only a rejected
    /// candidate enters the trace-scaled ridge loop, beginning at <c>1E-10 * trace / p</c> and multiplying the
    /// ridge by ten after each rejection until a finite candidate is accepted. These tests pin both the no-ridge and fallback
    /// paths so that a conditioning-policy change cannot silently alter downstream fitted results.
    /// </para>
    /// <para>
    /// The loop is structurally immune to the scale-relative pivot test at any realistic dimension. For a
    /// positive semi-definite input the ridged matrix has a smallest eigenvalue of at least the ridge, so
    /// every pivot is at least <c>1E-10 * trace / p</c> while no diagonal entry exceeds the trace. The pivot
    /// ratio is therefore at least <c>1E-10 / p</c>, against a tolerance of <c>p * 2^-52</c>. Those two cross
    /// only near <c>p = 671</c>; below that the first attempt always succeeds.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_MatrixRegularization
    {
        /// <summary>
        /// Asserts two matrices agree entry for entry.
        /// </summary>
        /// <param name="expected">The expected matrix.</param>
        /// <param name="actual">The matrix produced by the method under test.</param>
        /// <param name="delta">The permitted absolute difference per entry.</param>
        private static void AssertMatricesEqual(Matrix expected, Matrix actual, double delta)
        {
            Assert.AreEqual(expected.NumberOfRows, actual.NumberOfRows);
            Assert.AreEqual(expected.NumberOfColumns, actual.NumberOfColumns);
            for (int i = 0; i < expected.NumberOfRows; i++)
            {
                for (int j = 0; j < expected.NumberOfColumns; j++)
                    Assert.AreEqual(expected[i, j], actual[i, j], delta, "entry [" + i + "," + j + "]");
            }
        }

        /// <summary>
        /// Builds the matrix the first ridge attempt produces: the symmetrized input plus
        /// <c>1E-10 * trace / p</c> on the diagonal.
        /// </summary>
        /// <param name="M">The input matrix.</param>
        /// <returns>The candidate the loop tests first.</returns>
        private static Matrix FirstRidgeCandidate(Matrix M)
        {
            int p = M.NumberOfRows;
            var S = new Matrix(p);
            for (int i = 0; i < p; i++)
            {
                for (int j = 0; j < p; j++)
                    S[i, j] = 0.5d * (M[i, j] + M[j, i]);
            }
            double trace = 0d;
            for (int i = 0; i < p; i++) trace += S[i, i];
            double ridge = trace > 0d ? 1E-10d * trace / p : 1E-10d;
            for (int i = 0; i < p; i++) S[i, i] += ridge;
            return S;
        }

        /// <summary>
        /// Verifies that a well-conditioned symmetric matrix is returned without a ridge.
        /// </summary>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_WellConditionedReturnsWithoutRidge()
        {
            var M = new Matrix(new[,] { { 4d, 1d, 0.5d }, { 1d, 3d, 0.25d }, { 0.5d, 0.25d, 2d } });
            AssertMatricesEqual(M, MatrixRegularization.MakeSymmetricPositiveDefinite(M), 0d);
        }

        /// <summary>
        /// Verifies that a positive-definite matrix with widely separated coordinate scales is not changed.
        /// </summary>
        /// <remarks>
        /// The diagonal scales mirror real-space moment covariances. An unconditional trace-scaled ridge is
        /// dominated by the largest coordinate and materially changes the smallest coordinate even though
        /// each Cholesky pivot is healthy relative to its own diagonal.
        /// </remarks>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_ScaleSeparatedReturnsWithoutRidge()
        {
            var M = new Matrix(new[,] { { 1d, 0d, 0d }, { 0d, 1E4d, 0d }, { 0d, 0d, 1E8d } });
            AssertMatricesEqual(M, MatrixRegularization.MakeSymmetricPositiveDefinite(M), 0d);
        }

        /// <summary>
        /// Verifies that an exactly rank-deficient input is still resolved by the base ridge.
        /// </summary>
        /// <remarks>
        /// This is the motivating case for the scale-relative pivot test: the third row of the input
        /// equals the first, so the raw matrix is exactly rank two. The base ridge of
        /// <c>1E-10 * 5 / 3 = 1.666667E-10</c> lifts the smallest eigenvalue clear of the tolerance —
        /// the final pivot ratio is 1.666668E-10 against a tolerance of 6.66E-16 — so the first attempt
        /// succeeds under both the absolute and the scale-relative pivot tests, and the returned matrix
        /// is the same under either.
        /// </remarks>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_RankDeficientTakesTheBaseRidge()
        {
            var M = new Matrix(new[,] { { 2d, 0.5d, 2d }, { 0.5d, 1d, 0.5d }, { 2d, 0.5d, 2d } });
            var regularized = MatrixRegularization.MakeSymmetricPositiveDefinite(M);
            AssertMatricesEqual(FirstRidgeCandidate(M), regularized, 0d);

            // The candidate the loop accepted really is accepted by the tightened test.
            var chol = new CholeskyDecomposition(regularized);
            Assert.IsTrue(chol.IsPositiveDefinite);
            Assert.IsGreaterThan(chol.RelativeTolerance, chol.L[2, 2] * chol.L[2, 2] / regularized[2, 2]);
        }

        /// <summary>
        /// Verifies that an asymmetric input whose symmetric part is positive definite is only symmetrized.
        /// </summary>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_SymmetrizesFirst()
        {
            var M = new Matrix(new[,] { { 2d, 0.8d }, { 0.2d, 2d } });
            var regularized = MatrixRegularization.MakeSymmetricPositiveDefinite(M);
            var expected = new Matrix(new[,] { { 2d, 0.5d }, { 0.5d, 2d } });
            AssertMatricesEqual(expected, regularized, 0d);
            Assert.AreEqual(0.5d, regularized[0, 1], 0d);
            Assert.AreEqual(0.5d, regularized[1, 0], 0d);
        }

        /// <summary>
        /// Continues the established ridge ladder after the eighth rejected candidate, returning an SPD matrix.
        /// </summary>
        /// <remarks>
        /// The eigenvalues of [[1,2],[2,1]] are -1 and 3. A ridge of 1 leaves a zero eigenvalue;
        /// the first accepted decade is 10, giving eigenvalues 9 and 13 and determinant 117.
        /// The former unchecked 1E-4 ridge left the matrix indefinite.
        /// </remarks>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_IndefiniteContinuesTheCheckedRidgeLadder()
        {
            var original = new Matrix(new[,] { { 1d, 2d }, { 2d, 1d } });
            var regularized = MatrixRegularization.MakeSymmetricPositiveDefinite(original);
            AssertMatricesEqual(new Matrix(new[,] { { 11d, 2d }, { 2d, 11d } }), regularized, 0d);
            Assert.AreEqual(117d, new CholeskyDecomposition(regularized).Determinant(), 1E-12d);
            AssertMatricesEqual(new Matrix(new[,] { { 1d, 2d }, { 2d, 1d } }), original, 0d);
        }

        /// <summary>
        /// Each of the original eight ridge candidates retains its exact value and selection boundary.
        /// </summary>
        /// <param name="decade">The zero-based index of the first ridge larger than the negative eigenvalue.</param>
        /// <remarks>The diagonal fixture has known eigenvalues, so its required shift is analytical.</remarks>
        [TestMethod]
        [DataRow(0)]
        [DataRow(1)]
        [DataRow(2)]
        [DataRow(3)]
        [DataRow(4)]
        [DataRow(5)]
        [DataRow(6)]
        [DataRow(7)]
        public void Test_MakeSymmetricPositiveDefinite_OriginalEightCandidatesAreUnchanged(int decade)
        {
            double negativeEigenvalue = -5E-11d * Math.Pow(10d, decade);
            var matrix = new Matrix(new[,] { { negativeEigenvalue, 0d }, { 0d, 2d } });
            double expectedRidge = (1E-10d * (negativeEigenvalue + 2d) / 2d) * Math.Pow(10d, decade);
            var actual = MatrixRegularization.MakeSymmetricPositiveDefinite(matrix);
            AssertMatricesEqual(new Matrix(new[,] { { negativeEigenvalue + expectedRidge, 0d }, { 0d, 2d + expectedRidge } }), actual, 0d);
            Assert.IsTrue(actual[0, 0] > 0d && actual[1, 1] > 0d);
        }

        /// <summary>
        /// The captured B17C indefinite moment matrix requires the ninth trace-scaled ridge.
        /// </summary>
        /// <remarks>
        /// Independent NumPy eigvalsh results give a smallest eigenvalue of -9.092331071132015E-6
        /// before repair and 4.71784861949334E-5 after adding the ninth ridge. A positive determinant
        /// alone is insufficient, so every leading principal minor is also checked (Sylvester's criterion).
        /// </remarks>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_CapturedB17CWeightHasPositivePrincipalMinors()
        {
            var original = new Matrix(new[,]
            {
                { 0.01634642985063478d, -0.0005398365277961114d, 0.0007743723735985425d },
                { -0.0005398365277961114d, 0.0004912013373976332d, -0.00011353771970719928d },
                { 0.0007743723735985425d, -0.00011353771970719928d, 4.361399178720463e-05d }
            });
            double baseRidge = 1E-10d * (original[0, 0] + original[1, 1] + original[2, 2]) / 3d;
            double expectedRidge = (baseRidge * 1E7d) * 10d;
            var actual = MatrixRegularization.MakeSymmetricPositiveDefinite(original);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(original[i, j] + (i == j ? expectedRidge : 0d), actual[i, j], 0d);
            }
            double leadingTwo = actual[0, 0] * actual[1, 1] - actual[0, 1] * actual[1, 0];
            double determinant = actual[0, 0] * (actual[1, 1] * actual[2, 2] - actual[1, 2] * actual[2, 1])
                - actual[0, 1] * (actual[1, 0] * actual[2, 2] - actual[1, 2] * actual[2, 0])
                + actual[0, 2] * (actual[1, 0] * actual[2, 1] - actual[1, 1] * actual[2, 0]);
            Assert.IsTrue(actual[0, 0] > 0d && leadingTwo > 0d && determinant > 0d);
            Assert.IsTrue(new CholeskyDecomposition(actual).IsPositiveDefinite);
        }

        /// <summary>
        /// Zero and negative traces retain the existing absolute base ridge and receive checked escalation.
        /// </summary>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_NonPositiveTraceUsesExistingBaseRidge()
        {
            AssertMatricesEqual(new Matrix(new[,] { { 1E-10d, 0d }, { 0d, 1E-10d } }),
                MatrixRegularization.MakeSymmetricPositiveDefinite(new Matrix(2)), 0d);
            var actual = MatrixRegularization.MakeSymmetricPositiveDefinite(new Matrix(new[,] { { -1d, 0d }, { 0d, -2d } }));
            AssertMatricesEqual(new Matrix(new[,] { { 9d, 0d }, { 0d, 8d } }), actual, 0d);
        }

        /// <summary>
        /// Invalid matrix arguments are distinguished from expected rejected positive-definiteness probes.
        /// </summary>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_InvalidArgumentsAreRejected()
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(null));
            Assert.ThrowsExactly<ArgumentException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(new Matrix(2, 3)));
            foreach (double value in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            {
                Assert.ThrowsExactly<ArgumentException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(
                    new Matrix(new[,] { { value, 0d }, { 0d, 1d } })));
                Assert.ThrowsExactly<ArgumentException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(
                    new Matrix(new[,] { { 1d, value }, { value, 1d } })));
            }
        }

        /// <summary>
        /// Unrepresentable symmetrization, trace scale, and exhausted finite ridge candidates fail explicitly.
        /// </summary>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_UnrepresentableRepairThrows()
        {
            // Symmetrization overflows before a finite candidate can be formed.
            Assert.ThrowsExactly<InvalidOperationException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(
                new Matrix(new[,] { { double.MaxValue, 0d }, { 0d, 1d } })));
            // The finite zero-trace input requires a shift exceeding 8E307. The next decade overflows the positive diagonal.
            Assert.ThrowsExactly<InvalidOperationException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(
                new Matrix(new[,] { { -8E307d, 0d }, { 0d, 8E307d } })));
            // All finite decade shifts remain indefinite; the next ridge itself overflows.
            Assert.ThrowsExactly<InvalidOperationException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(
                new Matrix(new[,] { { -8E307d, 8E307d }, { 8E307d, -8E307d } })));
            // A finite indefinite matrix can still overflow when its diagonal entries are summed.
            var overflowingTrace = new Matrix(new[,]
            {
                { 8E307d, 0d, 0d, 0d }, { 0d, 8E307d, 0d, 0d },
                { 0d, 0d, 8E307d, 0d }, { 0d, 0d, 0d, -8E307d }
            });
            Assert.ThrowsExactly<InvalidOperationException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(overflowingTrace));
            // A trace-scaled starting ridge underflows to zero; do not silently select a different scale.
            Assert.ThrowsExactly<InvalidOperationException>(() => MatrixRegularization.MakeSymmetricPositiveDefinite(
                new Matrix(new[,] { { 0d, 0d }, { 0d, 1E-320d } })));
        }

        /// <summary>
        /// Expected rejected pivots must not construct or throw exceptions during ridge selection.
        /// </summary>
        /// <remarks>
        /// First-chance observation detects exceptions even when the regularizer catches them internally.
        /// The indefinite fixture is the captured B17C realization 213 moment covariance (seed 12345).
        /// </remarks>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_ExpectedRejectionsDoNotThrow()
        {
            int rejectedPivotExceptions = 0;
            int threadId = System.Threading.Thread.CurrentThread.ManagedThreadId;
            EventHandler<System.Runtime.ExceptionServices.FirstChanceExceptionEventArgs> handler = (sender, args) =>
            {
                if (System.Threading.Thread.CurrentThread.ManagedThreadId == threadId &&
                    args.Exception.Message.StartsWith("Cholesky Decomposition failed.", StringComparison.Ordinal))
                    rejectedPivotExceptions++;
            };
            AppDomain.CurrentDomain.FirstChanceException += handler;
            try
            {
                MatrixRegularization.MakeSymmetricPositiveDefinite(
                    new Matrix(new[,] { { 2d, 0.5d, 2d }, { 0.5d, 1d, 0.5d }, { 2d, 0.5d, 2d } }));
                MatrixRegularization.MakeSymmetricPositiveDefinite(new Matrix(new[,]
                {
                    { 0.01634642985063478d, -0.0005398365277961114d, 0.0007743723735985425d },
                    { -0.0005398365277961114d, 0.0004912013373976332d, -0.00011353771970719928d },
                    { 0.0007743723735985425d, -0.00011353771970719928d, 4.361399178720463e-05d }
                }));
            }
            finally
            {
                AppDomain.CurrentDomain.FirstChanceException -= handler;
            }
            Assert.AreEqual(0, rejectedPivotExceptions);
        }

        /// <summary>
        /// Verifies that <see cref="MatrixRegularization.Regularize(Matrix, double, double)"/> is untouched,
        /// since it floors eigenvalues directly and never consults the Cholesky factorization.
        /// </summary>
        [TestMethod]
        public void Test_Regularize_FloorsAndCapsEigenvalues()
        {
            var M = new Matrix(new[,] { { 1d, 1d }, { 1d, 1d } });
            var regularized = MatrixRegularization.Regularize(M);

            // Eigenvalues are 2 and 0; the floor is eps * trace / p = 1E-6 * 2 / 2 = 1E-6, and the cap is
            // 50 * median(0, 2) = 50, which binds on neither.
            Assert.AreEqual(1.0000005d, regularized[0, 0], 1E-12d);
            Assert.AreEqual(1.0000005d, regularized[1, 1], 1E-12d);
            Assert.AreEqual(0.9999995d, regularized[0, 1], 1E-12d);
            Assert.AreEqual(0.9999995d, regularized[1, 0], 1E-12d);
        }
    }
}
