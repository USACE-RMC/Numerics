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
    /// A class unit testing all Matrix Operations.
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
    public class Test_Matrix
    {
        /// <summary>
        /// Test structural methods for the matrix class.
        /// </summary>
        [TestMethod()]
        public void Test_StructuralMethods()
        {
            var A = new Matrix(3);
            A[0, 0] = 1d;
            A[0, 1] = 2d;
            A[0, 2] = 3d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = 6d;
            A[2, 0] = 7d;
            A[2, 1] = 8d;
            A[2, 2] = 9d;

            ///<summary>
            /// Testing .Row
            /// </summary>
            var true_row = new[] { 4d, 5d, 6d };
            var row = A.Row(1);
            for (int i = 0; i < row.Length; i++)
                Assert.AreEqual(row[i], true_row[i]);

            ///<summary>
            /// Testing .Column
            /// </summary>
            var true_col = new[] { 2d, 5d, 8d };
            var col = A.Column(1);
            for (int i = 0; i < col.Length; i++)
                Assert.AreEqual(col[i], true_col[i]);

            ///<summary> 
            ///Testing .Diagonal
            /// </summary>
            var true_diag = new[] { 1d, 5d, 9d };
            var diag = A.Diagonal();
            for (int i = 0; i < diag.Length; i++)
                Assert.AreEqual(diag[i], true_diag[i]);

            ///<summary>
            ///Testing upper triangular
            /// </summary>
            var true_upT = new Matrix(3);
            true_upT[0, 0] = 1d;
            true_upT[0, 1] = 2d;
            true_upT[0, 2] = 3d;
            true_upT[1, 1] = 5d;
            true_upT[1, 2] = 6d;
            true_upT[2, 2] = 9d;
            var upT = A.UpperTriangle();
            for (int i = 0; i < upT.NumberOfRows; i++)
            {
                for (int j = 0; j < upT.NumberOfColumns - 1; j++)
                    Assert.AreEqual(upT[i, j], true_upT[i, j]);
            }

            ///<summary>
            ///Testing lower triangular
            /// </summary>
            var true_lowT = new Matrix(3);
            true_lowT[0, 0] = 1d;
            true_lowT[1, 0] = 4d;
            true_lowT[1, 1] = 5d;
            true_lowT[2, 0] = 7d;
            true_lowT[2, 1] = 8d;
            true_lowT[2, 2] = 9d;
            var lowT = A.LowerTriangle();
            for (int i = 0; i < lowT.NumberOfRows; i++)
            {
                for (int j = 0; j < lowT.NumberOfColumns; j++)
                    Assert.AreEqual(lowT[i, j], true_lowT[i, j]);
            }
        }

        /// <summary>
        /// Test is symmetric.
        /// </summary>
        [TestMethod()]
        public void Test_IsSymmetric()
        {
            var A = new Matrix(3);
            A[0, 0] = 1d;
            A[0, 1] = 0d;
            A[0, 2] = 0d;
            A[1, 0] = 0d;
            A[1, 1] = 1d;
            A[1, 2] = 0d;
            A[2, 0] = 0d;
            A[2, 1] = 0d;
            A[2, 2] = 1d;
            bool true_result = true;
            bool result = A.IsSymmetric;
            Assert.AreEqual(result, true_result);
            //
            var B = new Matrix(3);
            B[0, 0] = 1d;
            B[0, 1] = 1d;
            B[0, 2] = 1d;
            B[1, 0] = 0d;
            B[1, 1] = 1d;
            B[1, 2] = 0d;
            B[2, 0] = 0d;
            B[2, 1] = 0d;
            B[2, 2] = 1d;
            true_result = false;
            result = B.IsSymmetric;
            Assert.AreEqual(result, true_result);
        }

        ///<summary>
        /// Testing if matrix is square
        /// </summary>
        [TestMethod()]
        public void Test_IsSquare()
        {
            var A = new Matrix(3,4);
            A[0, 0] = 1d;
            A[0, 1] = 0d;
            A[0, 2] = 0d;
            A[0, 3] = 1d;
            A[1, 0] = 0d;
            A[1, 1] = 1d;
            A[1, 2] = 0d;
            A[1, 3] = 1d;
            A[2, 0] = 0d;
            A[2, 1] = 0d;
            A[2, 2] = 1d;
            A[2, 3] = 1d;
            Assert.IsFalse(A.IsSquare);
        }

        ///<summary>
        /// Testing if matrix is square
        /// </summary>
        [TestMethod()]
        public void Test_Clone()
        {
            var A = new Matrix(3, 4);
            A[0, 0] = 1d;
            A[0, 1] = 0d;
            A[0, 2] = 0d;
            A[0, 3] = 1d;
            A[1, 0] = 0d;
            A[1, 1] = 1d;
            A[1, 2] = 0d;
            A[1, 3] = 1d;
            A[2, 0] = 0d;
            A[2, 1] = 0d;
            A[2, 2] = 1d;
            A[2, 3] = 1d;
            var clone = A.Clone();
            for (int i = 0; i < A.NumberOfRows; i++)
            {
                for (int j = 0; j < A.NumberOfColumns; j++)
                    Assert.AreEqual(A[i, j], clone[i, j], 0.0001d);
            }
        }

        ///<summary>
        /// Test Inverse
        /// </summary>
        [TestMethod()]
        public void Test_Inverse()
        {
            var A = new Matrix(3);
            A[0, 0] = 1d;
            A[0, 1] = 1d;
            A[0, 2] = 1d;
            A[1, 0] = 0d;
            A[1, 1] = 2d;
            A[1, 2] = 5d;
            A[2, 0] = 2d;
            A[2, 1] = 5d;
            A[2, 2] = -1d;

            //Testing Inverse
            var true_invA = new Matrix(3);
            true_invA[0, 0] = -27d / -21d;
            true_invA[0, 1] = 6d / -21d;
            true_invA[0, 2] = 3d / -21d;
            true_invA[1, 0] = 10d / -21d;
            true_invA[1, 1] = -3d / -21d;
            true_invA[1, 2] = -5d / -21d;
            true_invA[2, 0] = -4d / -21d;
            true_invA[2, 1] = -3d / -21d;
            true_invA[2, 2] = 2d / -21d;
            var invA = A.Inverse();
            for (int i = 0; i < invA.NumberOfRows; i++)
            {
                for (int j = 0; j < invA.NumberOfColumns; j++)
                    Assert.AreEqual(invA[i, j], true_invA[i, j], 0.0001d);
            }
        }

        /// <summary>
        /// Test identity.
        /// </summary>
        [TestMethod()]
        public void Test_Identity()
        {
            var true_result = new Matrix(3);
            true_result[0, 0] = 1d;
            true_result[0, 1] = 0d;
            true_result[0, 2] = 0d;
            true_result[1, 0] = 0d;
            true_result[1, 1] = 1d;
            true_result[1, 2] = 0d;
            true_result[2, 0] = 0d;
            true_result[2, 1] = 0d;
            true_result[2, 2] = 1d;
            var result = Matrix.Identity(3);
            for (int i = 0; i < result.NumberOfRows ; i++)
            {
                for (int j = 0; j < result.NumberOfColumns; j++)
                    Assert.AreEqual(result[i, j], true_result[i, j]);
            }
        }

        /// <summary>
        /// Test transpose.
        /// </summary>
        [TestMethod()]
        public void Test_Transpose()
        {
            var A = new Matrix(3, 2);
            A[0, 0] = 1d;
            A[0, 1] = 2d;
            A[1, 0] = 3d;
            A[1, 1] = 4d;
            A[2, 0] = 5d;
            A[2, 1] = 6d;

            var true_result = new Matrix(2, 3);
            true_result[0, 0] = 1d;
            true_result[0, 1] = 3d;
            true_result[0, 2] = 5d;
            true_result[1, 0] = 2d;
            true_result[1, 1] = 4d;
            true_result[1, 2] = 6d;

            var result = Matrix.Transpose(A);
            for (int i = 0; i < result.NumberOfRows; i++)
            {
                for (int j = 0; j < result.NumberOfColumns; j++)
                    Assert.AreEqual(result[i, j], true_result[i, j]);
            }
        }

        /// <summary>
        /// Testing diagonal function
        /// </summary>
        [TestMethod()]
        public void Test_Diagonal()
        {
            var A = new Matrix(3);
            A[0, 0] = 1;
            A[0, 1] = 2;
            A[0, 2] = 3;
            A[1,0] = 4;
            A[1,1] = 5;
            A[1,2] = 6;
            A[2,0] = 7;
            A[2,1] = 8;
            A[2,2] = 9;
            var result = Matrix.Diagonal(A);

            var true_result = new Matrix(3);
            true_result[0, 0] = 1;
            true_result[0, 1] = 0;
            true_result[0, 2] = 0;
            true_result[1,0] = 0;
            true_result[1, 1] = 5;
            true_result[1, 2] = 0;
            true_result[2, 0] = 0;
            true_result[2, 1] = 0;
            true_result[2, 2] = 9;

            for (int i = 0; i < result.NumberOfRows; i++)
            {
                for (int j = 0; j < result.NumberOfColumns; j++)
                    Assert.AreEqual(result[i, j], true_result[i, j]);
            }
        }
        
        /// <summary>
        /// Sum all elements in the matrix
        /// </summary>
        [TestMethod()]
        public void Test_Sum()
        {
            var A = new Matrix(3);
            A[0,0] = 1;
            A[0,1] = 2;
            A[0,2] = 3;
            A[1,0] = 4;
            A[1,1] = 5;
            A[1,2] = 6;
            A[2,0] = 7;
            A[2,1] = 8;
            A[2,2] = 9;
            var sum = A.Sum();
            Assert.AreEqual(45, sum);
        }

        /// <summary>
        /// Add matrices
        /// </summary>
        [TestMethod()]
        public void Test_Add()
        {
            var A = new Matrix(3, 2);
            A[0, 0] = 5d;
            A[0, 1] = -2;
            A[1, 0] = -1;
            A[1, 1] = 3d;
            A[2, 0] = 1d;
            A[2, 1] = 0d;
            var B = new Matrix(3, 2);
            B[0, 0] = 2d;
            B[0, 1] = -2;
            B[1, 0] = 0d;
            B[1, 1] = 1d;
            B[2, 0] = 4d;
            B[2, 1] = -1;
            var true_result = new Matrix(3, 2);
            true_result[0, 0] = 7d;
            true_result[0, 1] = -4;
            true_result[1, 0] = -1;
            true_result[1, 1] = 4d;
            true_result[2, 0] = 5d;
            true_result[2, 1] = -1;
            var result = A + B;
            for (int i = 0; i < result.NumberOfRows; i++)
            {
                for (int j = 0; j < result.NumberOfColumns; j++)
                    Assert.AreEqual(result[i, j], true_result[i, j]);
            }
        }

        /// <summary>
        /// Subtract matrices.
        /// </summary>
        [TestMethod()]
        public void Test_Subtract()
        {
            var A = new Matrix(3, 2);
            A[0, 0] = 5d;
            A[0, 1] = -2;
            A[1, 0] = -1;
            A[1, 1] = 3d;
            A[2, 0] = 1d;
            A[2, 1] = 0d;
            var B = new Matrix(3, 2);
            B[0, 0] = 2d;
            B[0, 1] = -2;
            B[1, 0] = 0d;
            B[1, 1] = 1d;
            B[2, 0] = 4d;
            B[2, 1] = -1;
            var true_result = new Matrix(3, 2);
            true_result[0, 0] = 3d;
            true_result[0, 1] = 0d;
            true_result[1, 0] = -1;
            true_result[1, 1] = 2d;
            true_result[2, 0] = -3;
            true_result[2, 1] = 1d;
            var result = A - B;
            for (int i = 0; i < result.NumberOfRows; i++)
            {
                for (int j = 0; j < result.NumberOfColumns; j++)
                    Assert.AreEqual(result[i, j], true_result[i, j]);
            }
        }

        /// <summary>
        /// Multiply matrices.
        /// </summary>
        [TestMethod()]
        public void Test_MultiplyByMatrix()
        {
            var A = new Matrix(2, 3);
            A[0, 0] = 1d;
            A[0, 1] = 2d;
            A[0, 2] = 3d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = 6d;
            var B = new Matrix(3, 2);
            B[0, 0] = 7d;
            B[0, 1] = 8d;
            B[1, 0] = 9d;
            B[1, 1] = 10d;
            B[2, 0] = 11d;
            B[2, 1] = 12d;
            var true_result = new Matrix(2);
            true_result[0, 0] = 58d;
            true_result[0, 1] = 64d;
            true_result[1, 0] = 139d;
            true_result[1, 1] = 154d;
            var result = A * B;
            for (int i = 0; i < result.NumberOfRows; i++)
            {
                for (int j = 0; j < result.NumberOfColumns; j++)
                    Assert.AreEqual(result[i, j], true_result[i, j]);
            }

        }

        /// <summary>
        /// Multiple matrix by vector.
        /// </summary>
        [TestMethod()]
        public void Test_MultiplyByVector()
        {
            var vector = new[] { 2d, 1d, 3d };
            var A = new Matrix(3);
            A[0, 0] = 1d;
            A[0, 1] = 2d;
            A[0, 2] = 3d;
            A[1, 0] = 4d;
            A[1, 1] = 5d;
            A[1, 2] = 6d;
            A[2, 0] = 7d;
            A[2, 1] = 8d;
            A[2, 2] = 9d;
            var true_result = new[] { 13d, 31d, 49d };
            var result = A * vector;
            for (int i = 0; i < result.Length; i++)
                Assert.AreEqual(result[i], true_result[i]);
        }

        /// <summary>
        /// Multiply matrix by scalar. 
        /// </summary>
        [TestMethod()]
        public void Test_MultiplyByScalar()
        {
            double scalar = 5d;
            var A = new Matrix(3, 2);
            A[0, 0] = 5d;
            A[0, 1] = -2;
            A[1, 0] = -1;
            A[1, 1] = 3d;
            A[2, 0] = 1d;
            A[2, 1] = 0d;
            var true_result = new Matrix(3, 2);
            true_result[0, 0] = 25d;
            true_result[0, 1] = -10;
            true_result[1, 0] = -5;
            true_result[1, 1] = 15d;
            true_result[2, 0] = 5d;
            true_result[2, 1] = 0d;
            var result = A * scalar;
            for (int i = 0; i < result.NumberOfRows; i++)
            {
                for (int j = 0; j < result.NumberOfColumns; j++)
                    Assert.AreEqual(result[i, j], true_result[i, j]);
            }
        }

        /// <summary>
        /// Test determinant.
        /// </summary>
        [TestMethod()]
        public void Test_Determinant()
        {
            var A = new Matrix(2);
            A[0, 0] = 1d;
            A[0, 1] = 2d;
            A[1, 0] = 3d;
            A[1, 1] = 4d;
            // determinant
            double true_result = -2.0d;
            double result = A.Determinant();
            Assert.AreEqual(result, true_result);
        }
    }
}
