using System;

namespace Numerics.Mathematics.LinearAlgebra
{

    /// <summary>
    /// A class for solving a set of linear equations Gauss-Jordan elimination.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// Gauss-Jordan Elimination is a row reduction algorithm that correctly formats matrices
    /// so that linear systems in the form A*x = b can be solved. 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes: The art of Scientific Computing, Third Edition." Press et al. 2017. </item>
    /// <item>"Numerical Methods for Engineers, Second Edition.", D.V. Griffiths and I.M. Smith, Taylor and Francis Group, 2006. </item>
    /// <item><description> 
    /// <see href = "https://en.wikipedia.org/wiki/Gaussian_elimination" />
    /// </description></item>
    ///</list>
    ///</para>
    /// </remarks>
    public class GaussJordanElimination
    {

        /// <summary>
        /// Perform Gauss-Jordan elimination.
        /// </summary>
        /// <param name="A">The input matrix A.</param>
        /// <param name="B">The input matrix B containing the M right-hand side vectors.</param>
        /// <remarks>
        /// On output, A is replaced by its matrix inverse, and B is replaced by the corresponding set of solution vectors.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when A is not square or B does not have as many rows as A.</exception>
        /// <exception cref="ArgumentException">Thrown when A contains a NaN or infinite entry.</exception>
        public static void Solve(ref Matrix A, ref Matrix B)
        {
            // No right-hand sides. Replaces A by its inverse.
            if (B == null)
            {
                // Dummy vector with zero columns
                B = new Matrix(A.NumberOfRows, 0);
            }
            if (A.NumberOfRows != A.NumberOfColumns)
                throw new ArgumentOutOfRangeException(nameof(A), "The matrix A must be square.");
            if (B.NumberOfRows != A.NumberOfRows)
                throw new ArgumentOutOfRangeException(nameof(B), "The matrix B must have the same number of rows as A.");
            // A NaN never wins the pivot comparison, so it would leave a stale pivot; reject it up front.
            for (int r = 0; r < A.NumberOfRows; r++)
            {
                for (int c = 0; c < A.NumberOfColumns; c++)
                {
                    if (!Tools.IsFinite(A[r, c]))
                        throw new ArgumentException("The matrix A must contain only finite entries.", nameof(A));
                }
            }

            int n = A.NumberOfRows;
            int m = B.NumberOfColumns;
            var icol = default(int);
            var irow = default(int);
            double big;
            double dum;
            double pivinv;
            // These integer arrays are used for bookkeeping on the pivoting.
            var indxc = new int[n];
            var indxr = new int[n];
            var ipiv = new int[n];
            for (int j = 0; j < n; j++)
                ipiv[j] = 0;
            // This is the main loop over the columns to be reduced. 
            for (int i = 0; i < n; i++)
            {
                big = 0.0d;
                // This is the outer loop of the search for a pivot element.
                for (int j = 0; j < n; j++)
                {
                    if (ipiv[j] != 1.0d)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            if (ipiv[k] == 0.0d)
                            {
                                if (Math.Abs(A[j, k]) >= big)
                                {
                                    big = Math.Abs(A[j, k]);
                                    irow = j;
                                    icol = k;
                                }
                            }
                            else if (ipiv[k] > 1)
                            {
                                throw new Exception("Singular matrix");
                            }
                        }
                    }
                }

                ipiv[icol] += 1;
                // With the pivot element found, rows are interchanged, if needed, to put the pivot
                // element on the diagonal. The columns are not physically interchanged, only relabeled:
                // index-c[i], the column of the (i+1)th pivot element, is the (i+1)th column that is
                // reduced, while index-r[i] is the row in which that pivot element was originally located.
                // If index-r[i] <> index-c[i], there is an implied column interchange. With the form of 
                // bookkeeping, the solution B's will end up in the correct order, and the inverse matrix
                // with be scrambled by columns. 
                if (irow != icol)
                {
                    for (int l = 0; l < n; l++)
                    {
                        dum = A[irow, l];
                        A[irow, l] = A[icol, l];
                        A[icol, l] = dum;
                    }

                    for (int l = 0; l < m; l++)
                    {
                        dum = B[irow, l];
                        B[irow, l] = B[icol, l];
                        B[icol, l] = dum;
                    }
                }
                // Divide the pivot row by the pivot element, located at i-row and i-col.
                indxr[i] = irow;
                indxc[i] = icol;
                if (A[icol, icol] == 0.0d)
                {
                    throw new Exception("Singular matrix");
                }

                pivinv = 1.0d / A[icol, icol];
                A[icol, icol] = 1.0d;
                for (int l = 0; l < n; l++)
                    A[icol, l] *= pivinv;
                for (int l = 0; l < m; l++)
                    B[icol, l] *= pivinv;
                // Reduce the rows except for the pivot one.
                for (int ll = 0; ll < n; ll++)
                {
                    if (ll != icol)
                    {
                        dum = A[ll, icol];
                        A[ll, icol] = 0.0d;
                        for (int l = 0; l < n ; l++)
                            A[ll, l] -= A[icol, l] * dum;
                        for (int l = 0; l < m; l++)
                            B[ll, l] -= B[icol, l] * dum;
                    }
                }
            }
            // This is the of the main loop over columns of the reduction. It only remains to unscramble
            // the solution in view of the column interchanges. This is done by interchanging pairs of 
            // columns in the reverse order that the permutation was built up. 
            for (int l = n - 1; l >= 0; l -= 1)
            {
                if (indxr[l] != indxc[l])
                {
                    for (int k = 0; k < n ; k++)
                    {
                        dum = A[k, indxr[l]];
                        A[k, indxr[l]] = A[k, indxc[l]];
                        A[k, indxc[l]] = dum;
                    }
                }
            }
            // Done.
        }
    }
}