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

            /// Recreated Gauss Jordan test in R to compare the inverted A matrices.
            /// I utilized library(matlib), gaussianElimination(), and inv() functions. Test passed.
        }

    }
}
