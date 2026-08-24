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
    /// <see cref="MatrixRegularization.MakeSymmetricPositiveDefinite"/> symmetrizes its input, then adds a
    /// trace-scaled ridge of <c>1E-10 * trace / p</c> and retries with the ridge multiplied by ten each time
    /// the factorization is rejected, up to eight attempts. Tightening the Cholesky pivot test could in
    /// principle make the loop reject a matrix it used to accept, escalate to a larger ridge, and return a
    /// different matrix to every downstream consumer. These tests pin the returned matrix so that any such
    /// escalation shows up as a failure rather than as a silent change in a fitted result.
    /// </para>
    /// <para>
    /// The loop is structurally immune to the scale-relative pivot test at any realistic dimension. For a
    /// positive semi-definite input the ridged matrix has a smallest eigenvalue of at least the ridge, so
    /// every pivot is at least <c>1E-10 * trace / p</c> while no diagonal entry exceeds the trace. The pivot
    /// ratio is therefore at least <c>1E-10 / p</c>, against a tolerance of <c>p * 2^-52</c>. Those two cross
    /// only near <c>p = 671</c>; below that the first attempt always succeeds, exactly as it does today.
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
        /// Verifies that a well-conditioned symmetric matrix is returned with the base ridge only.
        /// </summary>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_WellConditionedTakesTheBaseRidge()
        {
            var M = new Matrix(new[,] { { 4d, 1d, 0.5d }, { 1d, 3d, 0.25d }, { 0.5d, 0.25d, 2d } });
            AssertMatricesEqual(FirstRidgeCandidate(M), MatrixRegularization.MakeSymmetricPositiveDefinite(M), 0d);
        }

        /// <summary>
        /// Verifies that an exactly rank-deficient input is still resolved by the base ridge.
        /// </summary>
        /// <remarks>
        /// This is the case the scale-relative pivot test was introduced for: the third row of the input
        /// equals the first, so the raw matrix is exactly rank two. The base ridge of
        /// <c>1E-10 * 5 / 3 = 1.667E-10</c> lifts the smallest eigenvalue clear of the tolerance —
        /// the final pivot ratio is about 8.3E-11 against a tolerance of 6.66E-16 — so the first attempt
        /// succeeds under both the old absolute test and the new relative one, and the returned matrix is
        /// unchanged.
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
        /// Verifies that an asymmetric input is symmetrized before the ridge is applied.
        /// </summary>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_SymmetrizesFirst()
        {
            var M = new Matrix(new[,] { { 2d, 0.8d }, { 0.2d, 2d } });
            var regularized = MatrixRegularization.MakeSymmetricPositiveDefinite(M);
            AssertMatricesEqual(FirstRidgeCandidate(M), regularized, 0d);
            Assert.AreEqual(0.5d, regularized[0, 1], 0d);
            Assert.AreEqual(0.5d, regularized[1, 0], 0d);
        }

        /// <summary>
        /// Verifies the indefinite path, where the ridge must escalate, is reached identically.
        /// </summary>
        /// <remarks>
        /// The input has an eigenvalue of -1, which no ridge in the loop's range can lift, so the method
        /// falls through to the last-resort ridge of <c>1E-4 * trace / p</c>. Every rejection along the way
        /// comes from a strictly negative pivot, which both the absolute and the relative test reject
        /// identically, so this path cannot move.
        /// </remarks>
        [TestMethod]
        public void Test_MakeSymmetricPositiveDefinite_IndefiniteFallsThroughToTheLastResortRidge()
        {
            var M = new Matrix(new[,] { { 1d, 2d }, { 2d, 1d } });
            var regularized = MatrixRegularization.MakeSymmetricPositiveDefinite(M);
            double expectedRidge = 1E-4d * 2d / 2d;
            Assert.AreEqual(1d + expectedRidge, regularized[0, 0], 0d);
            Assert.AreEqual(1d + expectedRidge, regularized[1, 1], 0d);
            Assert.AreEqual(2d, regularized[0, 1], 0d);
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
