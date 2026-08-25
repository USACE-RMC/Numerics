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
        /// <para>
        /// This case being caught must not be read as a guarantee that rank deficiency is caught in general.
        /// The test rejects most numerically rank-deficient matrices but not all of them: over 1,000 trials
        /// per dimension on a covariance estimated from <c>m = n - 1</c> observations of <c>n</c> variables,
        /// it removes about two thirds of the matrices the absolute test accepted, and 11% to 18% still
        /// factorize. A successful factorization is not a rank certificate;
        /// <c>DecompositionMethod.SingularValue</c> on <c>MultivariateNormal</c> is the only reliable rank
        /// test in this library.
        /// </para>
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

        /// <summary>
        /// A correlation whose second pivot ratio is exactly <c>3 * 2^-52</c>, straddling the tolerance
        /// between <c>n = 2</c> and <c>n = 4</c>.
        /// </summary>
        /// <remarks>
        /// The first factor row is exact — <c>L[0,0] = 1</c> and <c>L[1,0] = ρ</c> — and
        /// <c>1 - fl(ρ²)</c> is exact by Sterbenz, so the pivot is bit-reproducible at
        /// 6.661338147750939E-16 rather than being an artifact of rounding order.
        /// </remarks>
        private const double ThreeUlpCorrelation = 0.99999999999999967d;

        /// <summary>
        /// Verifies that the tolerance really does scale with the dimension, by presenting the same
        /// near-singular two-variable block at two different dimensions and getting opposite verdicts.
        /// </summary>
        /// <remarks>
        /// The block's pivot ratio is exactly <c>3 * 2^-52 = 6.661338E-16</c>. At <c>n = 2</c> the tolerance
        /// is <c>2 * 2^-52 = 4.440892E-16</c> and the block is accepted; padded with two independent unit
        /// variances to <c>n = 4</c> the tolerance is <c>4 * 2^-52 = 8.881784E-16</c> and the identical block
        /// is rejected. Nothing about the block changed, only the dimension. Were the <c>n</c> factor dropped
        /// from <see cref="CholeskyDecomposition.DefaultRelativeTolerance"/> the tolerance would be
        /// 2.220446E-16 at both sizes and the 4x4 assertion here would fail, so this test pins the scaling as
        /// behaviour and not merely as a formula.
        /// </remarks>
        [TestMethod]
        public void Test_ToleranceScalesWithDimension()
        {
            double rho = ThreeUlpCorrelation;
            double expectedPivot = 3d * 2d * Tools.DoubleMachineEpsilon;

            // n = 2: tolerance 2 * 2^-52, below the pivot, so the block factorizes.
            var small = new CholeskyDecomposition(new Matrix(new[,] { { 1d, rho }, { rho, 1d } }));
            Assert.IsTrue(small.IsPositiveDefinite);
            Assert.AreEqual(expectedPivot, small.L[1, 1] * small.L[1, 1], 1E-24d);
            Assert.IsGreaterThan(small.RelativeTolerance, small.L[1, 1] * small.L[1, 1]);

            // n = 4: the same block, padded. Tolerance 4 * 2^-52 now exceeds the pivot.
            var padded = new Matrix(new[,]
            {
                { 1d, rho, 0d, 0d },
                { rho, 1d, 0d, 0d },
                { 0d, 0d, 1d, 0d },
                { 0d, 0d, 0d, 1d }
            });
            Assert.IsLessThan(CholeskyDecomposition.DefaultRelativeTolerance(4), expectedPivot);
            var exception = AssertThrowsAny(() => new CholeskyDecomposition(padded));
            Assert.IsGreaterThanOrEqualTo(0, exception.Message.IndexOf("positive-definite", StringComparison.OrdinalIgnoreCase));

            // The same 4x4 is accepted once the tolerance is lowered below the pivot, confirming the pivot
            // itself — not some other property of the padded matrix — is what the dimension changed.
            Assert.IsTrue(new CholeskyDecomposition(padded, 2d * Tools.DoubleMachineEpsilon).IsPositiveDefinite);
        }

        /// <summary>
        /// Verifies the diagonal guard: when <c>A[i,i]</c> is not a positive finite number the threshold falls
        /// back to zero, so the outcome is identical to the purely absolute test.
        /// </summary>
        /// <remarks>
        /// The class is public and can be handed anything, so the guard is exercised with a zero matrix, a
        /// zero on the diagonal behind a healthy leading entry, a negative diagonal, an infinite diagonal and
        /// a NaN diagonal. Each is compared against the same matrix at <c>relativeTolerance = 0</c> rather
        /// than against a hard-coded verdict, which is exactly the property the guard is there to provide:
        /// a non-positive or non-finite diagonal cannot make the relative test behave differently from the
        /// absolute one. The infinite case is accepted by both, which is pre-existing behaviour, not an
        /// endorsement.
        /// </remarks>
        [TestMethod]
        public void Test_NonPositiveOrNonFiniteDiagonalFallsBackToTheAbsoluteTest()
        {
            var cases = new[]
            {
                new[,] { { 0d, 0d }, { 0d, 0d } },
                new[,] { { 1d, 0d }, { 0d, 0d } },
                new[,] { { -1d, 0d }, { 0d, 1d } },
                new[,] { { double.PositiveInfinity, 0d }, { 0d, 1d } },
                new[,] { { double.NaN, 0d }, { 0d, 1d } }
            };

            foreach (var entries in cases)
            {
                var A = new Matrix(entries);
                string label = "[[" + entries[0, 0].ToString("R") + ",...]]";

                CholeskyDecomposition withDefault = null, withZero = null;
                Exception defaultException = null, zeroException = null;
                try { withDefault = new CholeskyDecomposition(A); }
                catch (Exception ex) { defaultException = ex; }
                try { withZero = new CholeskyDecomposition(A, 0d); }
                catch (Exception ex) { zeroException = ex; }

                Assert.AreEqual(zeroException == null, defaultException == null,
                    "The relative and absolute tests must agree on " + label);
                if (defaultException != null)
                {
                    // A guarded diagonal takes the original branch, so the message is the original one.
                    Assert.AreEqual(zeroException.Message, defaultException.Message, label);
                    Assert.AreEqual("Cholesky Decomposition failed. The input matrix is not positive-definite.",
                        defaultException.Message, label);
                }
                else
                {
                    Assert.AreEqual(withZero.L[0, 0], withDefault.L[0, 0], 0d, label);
                    Assert.AreEqual(withZero.L[1, 1], withDefault.L[1, 1], 0d, label);
                }
            }
        }

    }
}

