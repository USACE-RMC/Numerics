using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class testing individual components of the QR Decomposition Method.
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
    public class Test_QRDecomposition
    {
        /// <summary>
        /// Test QR Decomposition with original matrix.
        /// </summary>
        [TestMethod]
        public void Test_QRDecomp()
        {
            var A = new Matrix(3);
            A[0, 0] = 1d; A[0, 1] = 1d; A[0, 2] = 1d;
            A[1, 0] = 0d; A[1, 1] = 2d; A[1, 2] = 5d;
            A[2, 0] = 2d; A[2, 1] = 5d; A[2, 2] = -1d;

            var qr = new QRDecomposition(A);
            var Q = qr.Q;
            var R = qr.RMatrix;
            var QR = Q * R;

            for (int i = 0; i < A.NumberOfRows; ++i)
                for (int j = 0; j < A.NumberOfColumns; ++j)
                    Assert.AreEqual(A[i, j], QR[i, j], 1e-10);
        }

        /// <summary>
        /// Testing Solve with vector input.
        /// </summary>
        [TestMethod]
        public void Test_SolveVector()
        {
            var A = new Matrix(3);
            A[0, 0] = 1d; A[0, 1] = 1d; A[0, 2] = 1d;
            A[1, 0] = 0d; A[1, 1] = 2d; A[1, 2] = 5d;
            A[2, 0] = 2d; A[2, 1] = 5d; A[2, 2] = -1d;

            var B = new Vector(new[] { 6d, -4, 27d });
            var qr = new QRDecomposition(A);
            var x = qr.Solve(B);

            // Verify that A * x ≈ B
            var Ax = A * x;
            for (int i = 0; i < B.Length; ++i)
                Assert.AreEqual(B[i], Ax[i], 1e-10);
        }

        /// <summary>
        /// Testing Solve with matrix input.
        /// </summary>
        [TestMethod]
        public void Test_SolveMatrix()
        {
            var A = new Matrix(3);
            A[0, 0] = 1d; A[0, 1] = 1d; A[0, 2] = 1d;
            A[1, 0] = 0d; A[1, 1] = 2d; A[1, 2] = 5d;
            A[2, 0] = 2d; A[2, 1] = 5d; A[2, 2] = -1d;

            var matB = new Matrix(3, 1);
            matB[0, 0] = 6d;
            matB[1, 0] = -4d;
            matB[2, 0] = 27d;

            var qr = new QRDecomposition(A);
            var matX = qr.Solve(matB);

            var AX = A * matX;
            for (int i = 0; i < AX.NumberOfRows; ++i)
                Assert.AreEqual(matB[i, 0], AX[i, 0], 1e-10);
        }

        /// <summary>
        /// Tests QR decomposition with a square non-symmetric matrix using Solve(Vector).
        /// </summary>
        [TestMethod]
        public void Test_QR_Solve_Vector_NonSymmetric_Square()
        {
            var A = new Matrix(new double[,]
            {
                { 1, 2, 3 },
                { 0, 1, 4 },
                { 5, 6, 0 }
            });

            var b = new Vector(new double[] { 1, 2, 3 });

            var qr = new QRDecomposition(A);
            var x = qr.Solve(b);
            var Ax = A * x;

            for (int i = 0; i < b.Length; i++)
                Assert.AreEqual(b[i], Ax[i], 1e-10);
        }

        /// <summary>
        /// Tests QR decomposition with a square non-symmetric matrix using Solve(Matrix).
        /// </summary>
        [TestMethod]
        public void Test_QR_Solve_Matrix_NonSymmetric_Square()
        {
            var A = new Matrix(new double[,]
            {
                { 1, 2, 3 },
                { 0, 1, 4 },
                { 5, 6, 0 }
            });

            var B = new Matrix(3, 1);
            B[0, 0] = 1;
            B[1, 0] = 2;
            B[2, 0] = 3;

            var qr = new QRDecomposition(A);
            var X = qr.Solve(B);
            var AX = A * X;

            for (int i = 0; i < A.NumberOfRows; i++)
                Assert.AreEqual(B[i, 0], AX[i, 0], 1e-10);
        }

        /// <summary>
        /// Tests QR decomposition with an overdetermined matrix using Solve(Vector).
        /// </summary>
        [TestMethod]
        public void Test_QR_Solve_Vector_NonSymmetric_Overdetermined()
        {
            var A = new Matrix(new double[,]
            {
                { 1, 1, 1 },
                { 1, 2, 3 },
                { 1, 3, 6 },
                { 1, 4, 10 }
            });

            var b = new Vector(new double[] { 6, 0, 0, 6 });

            var qr = new QRDecomposition(A);
            var x = qr.Solve(b);
            var Ax = A * x;

            for (int i = 0; i < b.Length; i++)
                Assert.AreEqual(b[i], Ax[i], 1e-10);
        }

        /// <summary>
        /// Tests QR decomposition with an overdetermined matrix using Solve(Matrix).
        /// </summary>
        [TestMethod]
        public void Test_QR_Solve_Matrix_NonSymmetric_Overdetermined()
        {
            var A = new Matrix(new double[,]
            {
                { 1, 1, 1 },
                { 1, 2, 3 },
                { 1, 3, 6 },
                { 1, 4, 10 }
            });

            var B = new Matrix(4, 1);
            B[0, 0] = 6;
            B[1, 0] = 0;
            B[2, 0] = 0;
            B[3, 0] = 6;

            var qr = new QRDecomposition(A);
            var X = qr.Solve(B);
            var AX = A * X;

            for (int i = 0; i < A.NumberOfRows; i++)
                Assert.AreEqual(B[i, 0], AX[i, 0], 1e-10);
        }

        /// <summary>
        /// Tests QR decomposition with an underdetermined matrix using Solve(Vector).
        /// </summary>
        [TestMethod]
        public void Test_QR_Solve_Vector_NonSymmetric_Underdetermined()
        {
            var A = new Matrix(new double[,]
            {
                { 2, 3, 5, 1 },
                { 1, 0, 2, 3 },
                { 0, 1, 4, 2 }
            });

            var b = new Vector(new double[] { 1, 2, 3 });

            var qr = new QRDecomposition(A);
            var x = qr.Solve(b);
            var Ax = A * x;

            for (int i = 0; i < b.Length; i++)
                Assert.AreEqual(b[i], Ax[i], 1e-10);
        }

        /// <summary>
        /// Tests QR decomposition with an underdetermined matrix using Solve(Matrix).
        /// </summary>
        [TestMethod]
        public void Test_QR_Solve_Matrix_NonSymmetric_Underdetermined()
        {
            var A = new Matrix(new double[,]
            {
                { 2, 3, 5, 1 },
                { 1, 0, 2, 3 },
                { 0, 1, 4, 2 }
            });

            var B = new Matrix(3, 1);
            B[0, 0] = 1;
            B[1, 0] = 2;
            B[2, 0] = 3;

            var qr = new QRDecomposition(A);
            var X = qr.Solve(B);
            var AX = A * X;

            for (int i = 0; i < A.NumberOfRows; i++)
                Assert.AreEqual(B[i, 0], AX[i, 0], 1e-10);
        }

        /// <summary>
        /// A numerically rank-deficient matrix is rejected by Solve instead of dividing by roundoff, while a
        /// full-rank least-squares problem is unaffected.
        /// </summary>
        [TestMethod()]
        public void Test_SolveRejectsRankDeficientMatrix()
        {
            // The second column is twice the first.
            var deficient = new QRDecomposition(new Matrix(new double[,] { { 1, 2 }, { 2, 4 }, { 3, 6 } }));
            Assert.ThrowsExactly<System.ArgumentException>(() => deficient.Solve(new Vector(new double[] { 1, 2, 3 })));
            Assert.ThrowsExactly<System.ArgumentException>(() => deficient.Solve(new Matrix(3, 1)));

            // y = 1 + 2x exactly at x = 0, 1, 2, 3: the least-squares fit recovers (1, 2).
            var A = new Matrix(new double[,] { { 1, 0 }, { 1, 1 }, { 1, 2 }, { 1, 3 } });
            var x = new QRDecomposition(A).Solve(new Vector(new double[] { 1, 3, 5, 7 }));
            Assert.AreEqual(1d, x[0], 1E-13);
            Assert.AreEqual(2d, x[1], 1E-13);
        }

        /// <summary>
        /// R matches numpy.linalg.qr on the textbook matrix [[12, -51, 4], [6, 167, -68], [-4, 24, -41]], and
        /// columns near the overflow and underflow limits are handled by the scaled norm: a column of two 1E200
        /// entries gives a finite R[0,0], and a column of two 1E-200 entries is triangularized.
        /// </summary>
        [TestMethod()]
        public void Test_RMatchesNumpyAndScaledNorm()
        {
            var A = new Matrix(new double[,] { { 12, -51, 4 }, { 6, 167, -68 }, { -4, 24, -41 } });
            var R = new QRDecomposition(A).RMatrix;
            var expected = new double[,] { { -14, -21, 14 }, { 0, -175, 70 }, { 0, 0, -35 } };
            // R is unique up to the sign of each row; the last 1 x 1 reflection differs from LAPACK's.
            for (int i = 0; i < 3; i++)
            {
                double sign = System.Math.Sign(R[i, i]) * System.Math.Sign(expected[i, i]);
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(expected[i, j], sign * R[i, j], 1E-12);
            }

            var big = new QRDecomposition(new Matrix(new double[,] { { 1E200 }, { 1E200 } })).RMatrix;
            Assert.AreEqual(-System.Math.Sqrt(2d) * 1E200, big[0, 0], 1E-15 * 1E200);
            Assert.AreEqual(0d, big[1, 0], 1E-15 * 1E200);
            var small = new QRDecomposition(new Matrix(new double[,] { { 1E-200 }, { 1E-200 } })).RMatrix;
            Assert.AreEqual(-System.Math.Sqrt(2d) * 1E-200, small[0, 0], 1E-15 * 1E-200);
            Assert.AreEqual(0d, small[1, 0], 1E-15 * 1E-200);
        }
    }
}
