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

            //<summary>
            // Testing .Row
            // </summary>
            var true_row = new[] { 4d, 5d, 6d };
            var row = A.Row(1);
            for (int i = 0; i < row.Length; i++)
                Assert.AreEqual(row[i], true_row[i]);

            //<summary>
            // Testing .Column
            // </summary>
            var true_col = new[] { 2d, 5d, 8d };
            var col = A.Column(1);
            for (int i = 0; i < col.Length; i++)
                Assert.AreEqual(col[i], true_col[i]);

            //<summary>
            //Testing .Diagonal
            // </summary>
            var true_diag = new[] { 1d, 5d, 9d };
            var diag = A.Diagonal();
            for (int i = 0; i < diag.Length; i++)
                Assert.AreEqual(diag[i], true_diag[i]);

            //<summary>
            //Testing upper triangular
            // </summary>
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

            //<summary>
            //Testing lower triangular
            // </summary>
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
        /// Testing Clone
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
        /// Testing Clone
        /// </summary>
        [TestMethod()]
        public void Test_XElement()
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
            var clone = new Matrix(A.ToXElement());
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

        /// <summary>
        /// Matrix multiplication is bit-identical to the direct i-j-k triple loop it replaced,
        /// including rectangular shapes, empty inner dimensions and special values.
        /// </summary>
        [TestMethod()]
        public void Test_MultiplyMatchesDirectTripleLoop()
        {
            var rng = new System.Random(42);
            int[][] shapes = { new[] { 1, 1, 1 }, new[] { 3, 4, 2 }, new[] { 7, 1, 5 }, new[] { 5, 0, 3 }, new[] { 17, 23, 11 }, new[] { 64, 64, 64 } };
            foreach (var shape in shapes)
            {
                var A = new Matrix(shape[0], shape[1]);
                var B = new Matrix(shape[1], shape[2]);
                for (int i = 0; i < A.NumberOfRows; i++)
                    for (int j = 0; j < A.NumberOfColumns; j++)
                        A[i, j] = (rng.NextDouble() - 0.5) * System.Math.Pow(10d, rng.Next(-8, 9));
                for (int i = 0; i < B.NumberOfRows; i++)
                    for (int j = 0; j < B.NumberOfColumns; j++)
                        B[i, j] = (rng.NextDouble() - 0.5) * System.Math.Pow(10d, rng.Next(-8, 9));
                AssertBitIdentical(DirectProduct(A, B), A * B);
            }

            var special = new Matrix(new[,] { { double.NaN, 1d, -0d }, { double.PositiveInfinity, double.Epsilon, 1E308 }, { -0d, -0d, 2d } });
            var other = new Matrix(new[,] { { 1d, -0d, 3d }, { 1E308, 0d, -1d }, { -0d, double.NegativeInfinity, 0.5d } });
            AssertBitIdentical(DirectProduct(special, other), special * other);
            AssertBitIdentical(DirectProduct(other, special), other.Multiply(special));
            AssertBitIdentical(DirectProduct(special, special), special * special);
        }

        /// <summary>
        /// The direct i-j-k matrix product, the reference for the vectorized product.
        /// </summary>
        private static Matrix DirectProduct(Matrix A, Matrix B)
        {
            var result = new Matrix(A.NumberOfRows, B.NumberOfColumns);
            for (int i = 0; i < A.NumberOfRows; i++)
            {
                for (int j = 0; j < B.NumberOfColumns; j++)
                {
                    double sum = 0.0d;
                    for (int k = 0; k < A.NumberOfColumns; k++)
                        sum += A[i, k] * B[k, j];
                    result[i, j] = sum;
                }
            }
            return result;
        }

        /// <summary>
        /// Asserts two matrices have the same shape and the same bit pattern in every entry.
        /// </summary>
        private static void AssertBitIdentical(Matrix expected, Matrix actual)
        {
            Assert.AreEqual(expected.NumberOfRows, actual.NumberOfRows);
            Assert.AreEqual(expected.NumberOfColumns, actual.NumberOfColumns);
            for (int i = 0; i < expected.NumberOfRows; i++)
                for (int j = 0; j < expected.NumberOfColumns; j++)
                    Assert.AreEqual(System.BitConverter.DoubleToInt64Bits(expected[i, j]), System.BitConverter.DoubleToInt64Bits(actual[i, j]));
        }

        /// <summary>
        /// IsSymmetric uses a tolerance relative to the largest entry, so it rejects [[1E-13, 0.5E-13], [0.6E-13, 1E-13]]
        /// (which an absolute 1E-12 tolerance accepts) and accepts large matrices with tiny relative asymmetry.
        /// </summary>
        [TestMethod()]
        public void Test_IsSymmetricIsScaleInvariant()
        {
            Assert.IsFalse(new Matrix(new double[,] { { 1E-13, 0.5E-13 }, { 0.6E-13, 1E-13 } }).IsSymmetric);
            Assert.IsTrue(new Matrix(new double[,] { { 1E6, 1E6 + 1E-7 }, { 1E6, 2E6 } }).IsSymmetric);
            Assert.IsFalse(new Matrix(new double[,] { { 1, 0.5 }, { 0.5 + 1E-9, 1 } }).IsSymmetric);
            Assert.IsTrue(new Matrix(new double[,] { { 1, 0.5 }, { 0.5 + 1E-13, 1 } }).IsSymmetric);
            Assert.IsTrue(new Matrix(2).IsSymmetric);
            Assert.IsFalse(new Matrix(new double[,] { { 1, double.NaN }, { double.NaN, 1 } }).IsSymmetric);
            Assert.IsFalse(new Matrix(2, 3).IsSymmetric);
        }
    }
}
