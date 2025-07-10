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

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class testing individual components of the Singular Value Decomposition Method.
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
    public class Test_SingularValueDecomp
    {
        /// <summary>
        /// Test Singular Value Decomposition methods with LUDecomposition.cs.
        /// </summary>
        [TestMethod()]
        public void Test_SVDWithLU()
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
            var lu = new LUDecomposition(A);
            var svd = new SingularValueDecomposition(A);

            // Test Solve x given B
            var B = new Vector(new[] { 6d, -4, 27d });
            var true_x = lu.Solve(B);
            var x = svd.Solve(B);
            for (int i = 0; i < x.Length; i++)
                Assert.AreEqual(x[i], true_x[i], 0.0001d);
        }

        /// <summary>
        /// Tests rank function.
        /// </summary>
        [TestMethod()]
        public void Test_Rank()
        {
            var A = new Matrix(2);
            A[0, 0] = 3;
            A[0, 1] = 0;
            A[1, 0] = 4;
            A[1, 1] = 5;

            var svd = new SingularValueDecomposition(A);
            var r = svd.Rank();
            Assert.AreEqual(2,r);
        }

        /// <summary>
        /// Tests nullity function. 
        /// </summary>
        [TestMethod()]
        public void Test_Nullity()
        {
            var A = new Matrix(3);
            A[0, 0] = 10;
            A[0, 1] = 20;
            A[0, 2] = 10;
            A[1, 0] = 20;
            A[1, 1] = 40;
            A[1, 2] = 20;
            A[2, 0] = 30;
            A[2, 1] = 50;
            A[2, 2] = 0;

            var svd = new SingularValueDecomposition(A);
            var nn = svd.Nullity();
            Assert.AreEqual(1, nn);
        }

        /// <summary>
        /// Tests known Range of matrix. 
        /// The function should return the column space of the matrix.
        /// </summary>
        [TestMethod()]
        public void Test_Range()
        {
            var A = new Matrix(2);
            A[0, 0] = 3;
            A[0, 1] = 0;
            A[1, 0] = 4;
            A[1, 1] = 5;
            var svd = new SingularValueDecomposition(A); ;

            var ranB = new Matrix(4, 3);
            ranB[0, 0] = 0;
            ranB[0, 1] = 1;
            ranB[0, 2] = 1;
            ranB[1,0] = 0;
            ranB[1, 1] = 0;
            ranB[1, 2] = 1;
            ranB[2, 0] = 1;
            ranB[2, 1] = 0;
            ranB[2, 2] = 0;
            ranB[3, 0] = 1;
            ranB[3, 1] = 1;
            ranB[3, 2] = 0;

            var true_range = new Matrix(3, 2);
            true_range[0, 0] = 1/Math.Sqrt(10);
            true_range[0, 1] = -3 / Math.Sqrt(10);
            true_range[1, 0] = 3/Math.Sqrt(10);
            true_range[1, 1] = 1 / Math.Sqrt(10);

            var range = svd.Range();

            for (int i = 0; i < range.NumberOfRows; i++)
            {
                for (int j = 0; j < range.NumberOfColumns; j++)
                    Assert.AreEqual(range[i, j], true_range[i, j], 0.0001d);
            }
        }

        /// <summary>
        /// Tests known Nullspace of a matrix 
        /// </summary>
        [TestMethod()]
        public void Test_Nullspace()
        {
            var B = new Matrix(2);
            B[0, 0] = 2;
            B[0, 1] = 0;
            B[1, 0] = 2;
            B[1, 1] = 0;

            var svd = new SingularValueDecomposition(B);
            var nullspace = svd.Nullspace();

            var bNull = new Matrix(2, 1);
            bNull[0, 0] = 0;
            bNull[1, 0] = 1;

            for (int i = 0; i < nullspace.NumberOfRows; i++)
            {
                Assert.AreEqual(nullspace[i, 0], bNull[i, 0], 0.0001d);

            }
        }

        /// <summary>
        /// Evaluating A*x=b with a b vector input
        /// </summary>
        [TestMethod()]
        public void Test_SolveVector()
        {
            var A = new Matrix(3);
            A[0, 0] = 2;
            A[0, 1] = -3;
            A[0, 2] = 1;
            A[1, 0] = 1;
            A[1, 1] = -1;
            A[1, 2] = -2;
            A[2, 0] = 3;
            A[2, 1] = 1;
            A[2, 2] = -1;

            var svd = new SingularValueDecomposition(A);

            var B = new Vector(new[] { 7d, -2d, 0 });

            //Testing Solve with vector input 
            var true_x = new[] { 1d, -1d, 2 };
            var x = svd.Solve(B);
            for (int i = 0; i < x.Length; i++)
                Assert.AreEqual(x[i], true_x[i], 0.0001d);
        }

        /// <summary>
        /// Evaluating A*x = B with a B matrix input 
        /// </summary>
        [TestMethod()]
        public void Test_SolveMatrix()
        {
            var A = new Matrix(3);
            A[0, 0] = 3;
            A[0, 1] = -2;
            A[0, 2] = 1;
            A[1, 0] = -2;
            A[1, 1] = 3;
            A[1, 2] = 2;
            A[2, 0] = 1;
            A[2, 1] = 2;
            A[2, 2] = 2;

            var svd = new SingularValueDecomposition(A);
            var B = new Vector(new[] { 3d, -3d, 2d });
            var matB = new Matrix(3, 1);
            matB[0, 0] = 3;
            matB[1, 0] = -3;
            matB[2,0] = 2;

            //Testing Solve with Matrix B input
            var true_matX = new Matrix(3, 1);
            true_matX[0, 0] = 2d;
            true_matX[1, 0] = 1d;
            true_matX[2, 0] = -1d;

            var matX = svd.Solve(matB);

            for (int i = 0; i < matX.NumberOfRows; i++)
            {
                Assert.AreEqual(matX[i, 0], true_matX[i, 0], 0.0001d);

            }
        }

        /// <summary>
        /// Testing known square matrix with positive singular values to check correct decomposition
        /// </summary>
        [TestMethod()]
        public void Test_Decompose()
        {
            var A = new Matrix(2);
            A[0, 0] = 3;
            A[0, 1] = 0;
            A[1, 0] = 4;
            A[1, 1] = 5;

            var svd = new SingularValueDecomposition(A);

            var U = new Matrix(2);
            U[0, 0] = 1;
            U[0, 1] = -3;
            U[1, 0] = 3;
            U[1, 1] = 1;
            U = (1 / Math.Sqrt(10)) * U;

            var W = new Vector(new[] { Math.Sqrt(45), Math.Sqrt(5) });

            var V = new Matrix(2);
            V[0, 0] = 1;
            V[0, 1] = -1;
            V[1, 0] = 1;
            V[1, 1] = 1;
            V = (1 / Math.Sqrt(2)) * V;

            for (int i = 0; i < U.NumberOfRows; i++)
            {
                for (int j = 0; j < U.NumberOfColumns; j++)
                    Assert.AreEqual(U[i, j], svd.U[i, j], 0.0001d);
            }
            for (int i = 0; i < W.Length; i++)
                Assert.AreEqual(W[i], svd.W[i], 0.0001d);

            for (int i = 0; i < V.NumberOfRows; i++)
            {
                for (int j = 0; j < V.NumberOfColumns; j++)
                    Assert.AreEqual(V[i, j], svd.V[i, j], 0.0001d);
            }

            ///Tested with the Test_Decompose() matrix to see if U and V were equal
            ///svd() function assigns switched left and right singular vectors, however 
            ///the overall decomposition is the same for both functions. Test passed.
        }

        /// <summary>
        /// Testing log determinant of decomposed matrix 
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
            var svd = new SingularValueDecomposition(A);

            // Test Determinant
            double true_det = 6.356108; //Math.Log(576) 
            double det = svd.LogDeterminant();
            Assert.AreEqual(det, true_det, 0.0001d);
        }

        /// <summary>
        /// Testing log of determinant of the Pseudo A Matrix
        /// </summary>
        [TestMethod()]
        public void Test_LogPseudoDeterminant()
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
            var svd = new SingularValueDecomposition(A);

            // Test Determinant
            double true_det = 6.356108; //Math.Log(576) 
            double det = svd.LogPseudoDeterminant();
            Assert.AreEqual(det, true_det, 0.0001d);
        }

    }
}
