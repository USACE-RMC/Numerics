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
    /// so that we can solve linear systems in the form A*x = b. 
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
        public static void Solve(ref Matrix A, ref Matrix B)
        {
            // No right-hand sides. Replaces A by its inverse.
            if (B == null)
            {
                // Dummy vector with zero columns
                B = new Matrix(A.NumberOfRows, 0);
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
                // We now have the pivot element, so we interchange rows, if needed, to put the pivot
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
                // We are now ready to divide the pivot row by the pivot element, located at i-row and i-col.
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
                // Now we reduce the rows except for the pivot one, of course.
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
            // the solution in view of the column interchanges. We do this by interchanging pairs of 
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
            // And we are done.
        }
    }
}