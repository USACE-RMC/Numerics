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
    public class Test_CholeskyDecomp
    {
        /// <summary>
        /// Test Cholesky Decomposition methods with LUDecomposition.cs.
        /// </summary>
        [TestMethod()]
        public void Test_CholeskyDecompMethods()
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
        /// <para>
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
  
    }
}

