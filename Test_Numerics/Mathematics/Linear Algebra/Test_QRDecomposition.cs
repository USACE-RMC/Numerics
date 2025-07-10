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
    }
}
