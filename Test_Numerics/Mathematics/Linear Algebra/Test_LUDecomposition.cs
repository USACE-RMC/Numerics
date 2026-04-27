using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class testing individual components of the LU Decomposition Method.
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
    public class Test_LUDecomposition
    {
        /// <summary>
        /// Test LU Decomposition with original matrix
        /// </summary>
        [TestMethod()]
        public void Test_LUDecomp()
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

            var B = new Vector(new[] { 6d, -4, 27d });
            var lu = new LUDecomposition(A);

            var L = new Matrix(3);
            L[0, 0] = 1d;
            L[0, 1] = 0d;
            L[0, 2] = 0d;
            L[1,0] = 0d;
            L[1, 1] = 1d;
            L[1, 2] = 0d;
            L[2, 0] = 2d;
            L[2, 1] = 1.5d;
            L[2, 2] = 1d;

            var U = new Matrix(3);
            U[0, 0] = 1d;
            U[0, 1] = 1d;
            U[0, 2] = 1d;
            U[1, 0] = 0d;
            U[1, 1] = 2d;
            U[1, 2] = 5d;
            U[2, 0] = 0d;
            U[2, 1] = 0d;
            U[2, 2] = -10.5d;
            for (int i = 0; i < A.NumberOfRows; ++i)
            {
                for (int j = 0; j < A.NumberOfColumns; ++j)
                {
                    Assert.AreEqual(A[i, j], lu.A[i, j], 0.0001d);
                }
            }

            /// Replicated LU decomposition in R and compared with decomposition results.
            /// Used lu() function for recreation. Test passed. 
        }
        /// <summary>
        /// Testing Solve with vector input 
        /// </summary>
        [TestMethod()]
        public void Test_SolveVector() {
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

            var B = new Vector(new[] { 6d, -4, 27d });
            var lu = new LUDecomposition(A);

            
            var true_x = new[] { 5d, 3d, -2 };
            var x = lu.Solve(B);
            for (int i = 0; i < x.Length; i++)
                Assert.AreEqual(x[i], true_x[i], 0.0001d);
        }

        /// <summary>
        /// Testing Solve with matrix input 
        /// </summary>
        [TestMethod()]
        public void Test_SolveMatrix() {
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

            var lu = new LUDecomposition(A);

            var matB = new Matrix(3, 1);
            matB[0, 0] = 6d;
            matB[1, 0] = -4d;
            matB[2, 0] = 27d;

            //Testing Solve with Matrix B
            var true_matX = new Matrix(3, 1);
            true_matX[0, 0] = 5d;
            true_matX[1, 0] = 3d;
            true_matX[2, 0] = -2d;

            var matX = lu.Solve(matB);

            for (int i = 0; i < matX.NumberOfRows; i++)
            {
                    Assert.AreEqual(matX[i,0], true_matX[i,0], 0.0001d);
                
            }
                    
        }

        /// <summary>
        /// Test determinant function on decomposed A
        /// </summary>
        [TestMethod()]
        public void Test_Determinant() {
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

            var B = new Vector(new[] { 6d, -4, 27d });
            var lu = new LUDecomposition(A);

            //Testing determinant
            double true_det = -21;
            double det = lu.Determinant();
            Assert.AreEqual(det, true_det, 0.0001d);

        }

        /// <summary>
        /// Testing inverse function on matrix A
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

            var B = new Vector(new[] { 6d, -4, 27d });
            var lu = new LUDecomposition(A);

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
            var invA = lu.InverseA();
            for (int i = 0; i < invA.NumberOfRows; i++)
            {
                for (int j = 0; j < invA.NumberOfColumns; j++)
                        Assert.AreEqual(invA[i, j], true_invA[i, j], 0.0001d);
            }
        }

    }
}
