using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class testing individual components of Gauss-Jordan Elimination.
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
    public class Test_GaussJordanElimination
    {
        /// <summary>
        /// Tests Gauss-Jordan elimination against a known solution.
        /// </summary>
        [TestMethod()]
        public void Test_GaussJordanElim()
        {
            var _matrix = new double[,] { { 1d, 3d, 3d }, { 1d, 4d, 3d }, { 1d, 3d, 4d } };
            var true_IA = new double[,] { { 7d, -3, -3 }, { -1, 1d, 0d }, { -1, 0d, 1d } };
            var A = new Matrix(_matrix);
            Matrix argB = null;
            GaussJordanElimination.Solve(ref A, B: ref argB);
            for (int i = 0; i < A.NumberOfRows; i++)
            {
                for (int j = 0; j < A.NumberOfColumns - 1; j++)
                    Assert.AreEqual(true_IA[i, j], A[i, j]);
            }

            // Recreated Gauss Jordan test in R to compare the inverted A matrices.
            // Uses the library(matlib) gaussianElimination() and inv() functions. Test passed.
        }

        /// <summary>
        /// Invalid inputs are rejected: a NaN entry (which never wins the pivot comparison), a non-square A, and a
        /// right-hand side with the wrong number of rows.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidInputs()
        {
            var A = new Matrix(new double[,] { { double.NaN, 1 }, { 1, 2 } });
            Matrix B = null!;
            Assert.ThrowsExactly<System.ArgumentException>(() => GaussJordanElimination.Solve(ref A, ref B));
            var square = new Matrix(new double[,] { { 2, 1 }, { 1, 2 } });
            var shortB = new Matrix(1, 1);
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => GaussJordanElimination.Solve(ref square, ref shortB));
            var wide = new Matrix(2, 3);
            Matrix none = null!;
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => GaussJordanElimination.Solve(ref wide, ref none));
        }
    }
}
