using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class testing individual components of the Cholesky Decomposition Method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_CholeskyDecomposition
    {
        /// <summary>
        /// Test Cholesky Decomposition methods with LUDecomposition.cs.
        /// </summary>
        [TestMethod()]
        public void Test_CholeskyDecomp()
        {
            var A = new Matrix(3);
            A[0, 0] = 25d;
            A[0, 1] = 15d;
            A[0, 2] = -5;
            A[1, 0] = 15d;
            A[1, 1] = 18d;
            A[1, 2] = 0d;
            A[2, 0] = -5;
            A[2, 1] = 0d;
            A[2, 2] = 11d;
            var chol = new CholeskyDecomposition(A);
            var lu = new LUDecomposition(A);

            // Test Determinant
            double true_det = Math.Log(lu.Determinant());
            double det = chol.LogDeterminant();
            Assert.AreEqual(det, true_det, 0.0001d);

            // Test Inverse
            var true_invA = lu.InverseA();
            var invA = chol.InverseA();
            for (int i = 0; i < invA.NumberOfRows; i++)
            {
                for (int j = 0; j < invA.NumberOfColumns; j++)
                    Assert.AreEqual(invA[i, j], true_invA[i, j], 0.0001d);
            }

            // Test Solve x given B
            var B = new Vector(new[]{ 6d, -4, 27d });
            var true_x = lu.Solve(B);
            var x = chol.Solve(B);
            for (int i = 0; i < x.Length; i++)
                Assert.AreEqual(x[i], true_x[i], 0.0001d);
        }

        /// <summary>
        /// Test Cholesky Decomposition method with known values from textbook.
        /// </summary>


        [TestMethod()]
        public void Test_CholeskyDecompVals()
        {

            var A = new Matrix(3);
            A[0, 0] = 16d;
            A[0, 1] = 4d;
            A[0, 2] = 8d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = -4d;
            A[2, 0] = 8d;
            A[2, 1] = -4d;
            A[2, 2] = 22d;

            var chol = new CholeskyDecomposition(A);

            var L = new Matrix(3);
            L[0, 0] = 4d;
            L[0, 1] = 0d;
            L[0, 2] = 0d;
            L[1, 0] = 1d;
            L[1, 1] = 2d;
            L[1, 2] = 0d;
            L[2, 0] = 2d;
            L[2, 1] = -3d;
            L[2, 2] = 3d;

            //Test Decomposition
            var prod = L;
            for (int i = 0; i < prod.NumberOfRows; i++)
            {
                for (int j = 0; j < prod.NumberOfColumns; j++)
                    Assert.AreEqual(L[i,j], chol.L[i, j], 0.0001d);
            }
        }
        /// <summary>
        /// Tests Log determinant of matrix A
        /// </summary>
        [TestMethod()]
        public void Test_LogDeterminant()
        {
            var A = new Matrix(3);
            A[0, 0] = 16d;
            A[0, 1] = 4d;
            A[0, 2] = 8d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = -4d;
            A[2, 0] = 8d;
            A[2, 1] = -4d;
            A[2, 2] = 22d;
            var chol = new CholeskyDecomposition(A);

            // Test Determinant
            double true_det = 6.356108; //Math.Log(576)
            double det = chol.LogDeterminant();
            Assert.AreEqual(det, true_det, 0.0001d);
        }

        /// <summary>
        /// Tests invertible matrix function
        /// </summary>
        [TestMethod()]
        public void Test_Inverse()
        {
            var A = new Matrix(3);
            A[0, 0] = 16d;
            A[0, 1] = 4d;
            A[0, 2] = 8d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = -4d;
            A[2, 0] = 8d;
            A[2, 1] = -4d;
            A[2, 2] = 22d;
            var chol = new CholeskyDecomposition(A);

            // Test Inverse
            var true_invA = new Matrix(3);
            true_invA[0, 0] = 0.1631944444d;
            true_invA[0, 1] = -0.2083333333d;
            true_invA[0, 2] = -0.09722222222d;
            true_invA[1, 0] = -0.2083333333d;
            true_invA[1, 1] = 0.5d;
            true_invA[1, 2] = 0.1666666667d;
            true_invA[2, 0] = -0.09722222222d;
            true_invA[2, 1] = 0.1666666667d;
            true_invA[2, 2] = 0.1111111111d;

            var invA = chol.InverseA();
            for (int i = 0; i < invA.NumberOfRows; i++)
            {
                for (int j = 0; j < invA.NumberOfColumns; j++)
                    Assert.AreEqual(invA[i, j], true_invA[i, j], 0.0001d);
            }
        }

        /// <summary>
        /// Tests solve function for the set of n linear equations A*x=b
        /// </summary>
        [TestMethod()]
        public void Test_Solve()
        {
            var A = new Matrix(3);
            A[0, 0] = 16d;
            A[0, 1] = 4d;
            A[0, 2] = 8d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = -4d;
            A[2, 0] = 8d;
            A[2, 1] = -4d;
            A[2, 2] = 22d;
            var chol = new CholeskyDecomposition(A);

            //Test Solve
            var B = new Vector(new[] { 16d, 18d, -22d });
            var trueX = new Vector(new[] { 1d, 2d, -1d });
            var X = chol.Solve(B);
            for (int i = 0; i < X.Length; i++)
                Assert.AreEqual(X[i], trueX[i], 0.0001d);
        }

        /// <summary>
        /// Test Forward-substitution method
        /// </summary>
        [TestMethod()]
        public void Test_Forward()
        {
            var A = new Matrix(3);
            A[0, 0] = 16d;
            A[0, 1] = 4d;
            A[0, 2] = 8d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = -4d;
            A[2, 0] = 8d;
            A[2, 1] = -4d;
            A[2, 2] = 22d;
            var chol = new CholeskyDecomposition(A);

            var b = new Vector(new double[] { 16d, 18d, -22d });
            var true_y = new double[] { 4d, 7d, -3d };
            var y = chol.Forward(b);
            for (int i = 0; i < y.Length; i++)
                Assert.AreEqual(y[i], true_y[i], 0.0001d);
        }

        /// <summary>
        /// Tests Back-substitution method
        /// </summary>
        [TestMethod()]
        public void Test_Back()
        {
            var A = new Matrix(3);
            A[0, 0] = 16d;
            A[0, 1] = 4d;
            A[0, 2] = 8d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = -4d;
            A[2, 0] = 8d;
            A[2, 1] = -4d;
            A[2, 2] = 22d;
            var chol = new CholeskyDecomposition(A);

            var Y = new Vector(new double[] { 4d, 7d, -3d });
            var right_x = new double[] { 1d, 2d, -1d };
            var x = chol.Backward(Y);
            for (int i = 0; i < x.Length; i++)
                Assert.AreEqual(x[i], right_x[i], 0.0001d);
        }

        /// <summary>
        /// The exactly rank-two covariance from USACE-RMC/Numerics#145, whose third row equals its first.
        /// </summary>
        private static readonly double[,] RankTwoCovariance =
        {
            { 2d, 0.5d, 2d },
            { 0.5d, 1d, 0.5d },
            { 2d, 0.5d, 2d }
        };

        /// <summary>
        /// Asserts the action throws, and returns the exception, without requiring a specific type
        /// (net481-compatible).
        /// </summary>
        /// <param name="action">The action expected to throw.</param>
        /// <returns>The exception that was thrown.</returns>
        private static Exception AssertThrowsAny(Action action)
        {
            try
            {
                action();
            }
            catch (Exception ex)
            {
                return ex;
            }
            Assert.Fail("The action was expected to throw.");
            throw new InvalidOperationException("unreachable");
        }

        /// <summary>
        /// Verifies that an exactly rank-deficient matrix is rejected rather than silently factorized.
        /// </summary>
        /// <remarks>
        /// Row three of <see cref="RankTwoCovariance"/> is identical to row one, so the matrix is exactly
        /// rank two and its third pivot is exactly zero in exact arithmetic. In floating point the pivot
        /// evaluates to 4.440892E-16 — a positive number — so the purely absolute <c>pivot &lt;= 0</c> test
        /// let the factorization complete and reported the matrix positive-definite. The pivot ratio is
        /// 2.220446E-16, which the default tolerance of <c>3 * 2^-52 = 6.661338E-16</c> rejects with a
        /// factor of three to spare.
        /// </remarks>
        [TestMethod]
        public void Test_RejectsExactlyRankDeficientMatrix()
        {
            var exception = AssertThrowsAny(() => new CholeskyDecomposition(new Matrix(RankTwoCovariance)));
            Assert.IsGreaterThanOrEqualTo(0, exception.Message.IndexOf("positive-definite", StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Pins the pre-existing behaviour of the absolute pivot test, reachable by passing a zero tolerance.
        /// </summary>
        /// <remarks>
        /// With <c>relativeTolerance = 0</c> the test reduces to <c>pivot &lt;= 0</c> exactly, so the rank-two
        /// covariance still factorizes and still yields the wrong answers it always did: a final factor entry
        /// of 2.107342425544702E-08 and a log determinant of -34.790890420621793, against a true log
        /// pseudo-determinant of 1.2527629684953678. This test exists so that the zero-tolerance escape hatch
        /// is demonstrably identical to the old behaviour, not so that the old behaviour is endorsed.
        /// </remarks>
        [TestMethod]
        public void Test_ZeroToleranceReproducesTheAbsolutePivotTest()
        {
            var chol = new CholeskyDecomposition(new Matrix(RankTwoCovariance), 0d);
            Assert.IsTrue(chol.IsPositiveDefinite);
            Assert.AreEqual(0d, chol.RelativeTolerance, 0d);
            Assert.AreEqual(2.107342425544702E-08d, chol.L[2, 2], 1E-14d);
            Assert.AreEqual(-34.790890420621793d, chol.LogDeterminant(), 1E-5d);
        }

        /// <summary>
        /// Verifies that a genuinely positive-definite but severely ill-conditioned correlation matrix is
        /// still accepted.
        /// </summary>
        /// <remarks>
        /// For two variables with correlation ρ the second pivot ratio is exactly <c>1 - ρ²</c>. At
        /// ρ = 0.9999 the ratio is 1.999900E-04, some 4.5E+11 times the tolerance of 4.440892E-16 that
        /// applies at n = 2. At ρ = 1 − 1E-12 — the tightest legitimately positive-definite matrix anywhere
        /// in this test suite — the ratio is 1.999956E-12, still 4,503 times the tolerance. A false
        /// rejection would require ρ within roughly two ulp of one, at which point the matrix is
        /// numerically singular in any case.
        /// </remarks>
        [TestMethod]
        public void Test_IllConditionedCorrelationIsStillAccepted()
        {
            foreach (double rho in new[] { 0.9999d, 1d - 1E-12d })
            {
                var A = new Matrix(new[,] { { 1d, rho }, { rho, 1d } });
                var chol = new CholeskyDecomposition(A);
                Assert.IsTrue(chol.IsPositiveDefinite, "rho = " + rho.ToString("R"));

                // The second pivot is L[1,1]^2 and must match 1 - rho^2 to rounding.
                Assert.AreEqual(1d - rho * rho, chol.L[1, 1] * chol.L[1, 1], 1E-16d);
                Assert.IsGreaterThan(chol.RelativeTolerance, chol.L[1, 1] * chol.L[1, 1]);
            }
        }

        /// <summary>
        /// Verifies that a covariance mixing a very large variance with a very small one is still accepted,
        /// which is the reason the tolerance is relative to each diagonal entry rather than to the largest.
        /// </summary>
        /// <remarks>
        /// For <c>[[1E12, 1], [1, 1E-11]]</c> the determinant is 9, so the matrix is comfortably positive
        /// definite, and the second pivot ratio is 0.9 — entirely healthy. Measured against the largest
        /// diagonal instead, the same pivot would read 9E-24 and any sane tolerance would reject it. The
        /// pivot is the conditional variance of the second variable, so its own diagonal is the only
        /// meaningful scale.
        /// </remarks>
        [TestMethod]
        public void Test_TinyVarianceBesideLargeVarianceIsStillAccepted()
        {
            var A = new Matrix(new[,] { { 1E12d, 1d }, { 1d, 1E-11d } });
            var chol = new CholeskyDecomposition(A);
            Assert.IsTrue(chol.IsPositiveDefinite);

            double pivotRatio = chol.L[1, 1] * chol.L[1, 1] / A[1, 1];
            Assert.AreEqual(0.9d, pivotRatio, 1E-12d);
            Assert.IsLessThan(1E-23d, chol.L[1, 1] * chol.L[1, 1] / A[0, 0],
                "The same pivot measured against the largest diagonal is negligible, which is why that scale is not used.");
        }

        /// <summary>
        /// Verifies the default tolerance value and that it is reported on the instance.
        /// </summary>
        /// <remarks>
        /// The default is <c>n * 2^-52</c>, expressed as <c>2 * n * Tools.DoubleMachineEpsilon</c> because
        /// that constant is the unit roundoff <c>2^-53</c>.
        /// </remarks>
        [TestMethod]
        public void Test_DefaultRelativeTolerance()
        {
            Assert.AreEqual(0d, CholeskyDecomposition.DefaultRelativeTolerance(0), 0d);
            Assert.AreEqual(0d, CholeskyDecomposition.DefaultRelativeTolerance(-4), 0d);
            Assert.AreEqual(2d * 3d * Tools.DoubleMachineEpsilon, CholeskyDecomposition.DefaultRelativeTolerance(3), 0d);
            Assert.AreEqual(6.661338147750960E-16d, CholeskyDecomposition.DefaultRelativeTolerance(3), 1E-30d);

            var A = new Matrix(new[,] { { 4d, 1d }, { 1d, 3d } });
            Assert.AreEqual(CholeskyDecomposition.DefaultRelativeTolerance(2), new CholeskyDecomposition(A).RelativeTolerance, 0d);
            Assert.AreEqual(1E-8d, new CholeskyDecomposition(A, 1E-8d).RelativeTolerance, 0d);
        }

        /// <summary>
        /// Verifies that an out-of-range or non-finite tolerance is rejected.
        /// </summary>
        [TestMethod]
        public void Test_InvalidRelativeToleranceIsRejected()
        {
            var A = new Matrix(new[,] { { 4d, 1d }, { 1d, 3d } });
            foreach (double tolerance in new[] { -1E-16d, 1d, 2d, double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            {
                var exception = AssertThrowsAny(() => new CholeskyDecomposition(A, tolerance));
                Assert.IsInstanceOfType(exception, typeof(ArgumentOutOfRangeException), "tolerance = " + tolerance.ToString("R"));
            }
        }

        /// <summary>
        /// Verifies that a matrix rejected by the absolute test is still rejected, with the same message, and
        /// that a well-conditioned matrix factorizes identically under either tolerance.
        /// </summary>
        /// <remarks>
        /// Only pivots that are strictly positive yet at or below the relative threshold change outcome. A
        /// negative or NaN pivot takes the original branch and therefore carries the original message
        /// verbatim, which the callers that catch and re-wrap this exception rely on.
        /// </remarks>
        [TestMethod]
        public void Test_NonPositiveDefiniteRejectionIsUnchanged()
        {
            var indefinite = new Matrix(new[,] { { 1d, 2d }, { 2d, 1d } });
            Assert.AreEqual(
                "Cholesky Decomposition failed. The input matrix is not positive-definite.",
                AssertThrowsAny(() => new CholeskyDecomposition(indefinite)).Message);
            Assert.AreEqual(
                "Cholesky Decomposition failed. The input matrix is not positive-definite.",
                AssertThrowsAny(() => new CholeskyDecomposition(indefinite, 0d)).Message);

            var A = new Matrix(new[,] { { 16d, 4d, 8d }, { 4d, 5d, -4d }, { 8d, -4d, 22d } });
            var withDefault = new CholeskyDecomposition(A);
            var withZero = new CholeskyDecomposition(A, 0d);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(withZero.L[i, j], withDefault.L[i, j], 0d);
            }
        }

    }
}

