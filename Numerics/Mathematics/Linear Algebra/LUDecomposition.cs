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
    /// A class for solving a set of linear equations using LU Decomposition.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b> Authors: </b>
    ///    <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// LU decomposition is a useful method utilizing Gaussian Elimination to decompose a matrix into 
    /// a lower and upper triangular matrix. 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes: The art of Scientific Computing, Third Edition." Press et al. 2017. </item>
    /// <item>"Applied Numerical Methods with MATLAB for Engineers and Scientists, Third Edition.", Steven C. Chapra, McGraw-Hill, 2012. </item>
    /// <item><description> 
    /// <see href = "https://en.wikipedia.org/wiki/LU_decomposition" />
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class LUDecomposition
    {
      
        /// <summary>
        /// Constructs new LU Decomposition.
        /// </summary>
        /// <param name="A">The input matrix A that is to be LU decomposed.</param>
        /// <remarks>
        /// The input matrix is not altered; a copy is made, on which outer product
        /// Gaussian elimination is then done in-place.
        /// </remarks>
        public LUDecomposition(Matrix A)
        {
            // Given a square matrix A[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a 
            // row-wise permutation of itself. Index[0..n-1] is an output vector that records the row permutation effected by the
            // partial pivoting; d is output as +-1 depending on whether the number of row interchanges was even or odd, respectively.
            // This routine is used in combination with Solve to solve linear equations or inverse a matrix. 

            this.A = new Matrix(A.ToArray());
            LU = new Matrix(A.ToArray());
            n = A.NumberOfRows;
            index = new int[n];
            var vv = new double[n];
            double TINY = 1.0E-40d;
            int imax, i, j, k;
            double big, temp;
            d = 1.0d;

            // Loop over rows to get the implicit scaling information.
            for (i = 0; i < n; i++)
            {
                big = 0.0d;
                for (j = 0; j < n; j++)
                {
                    temp = Math.Abs(LU[i, j]);
                    if (temp > big) big = temp;
                }

                if (big == 0.0d) throw new ArgumentException("Singular matrix in LU decomposition.");
                vv[i] = 1.0d / big;
            }
            // This is the outermost k-i-j loop. 
            for (k = 0; k < n; k++)
            {
                big = 0.0d;
                imax = k;
                for (i = k; i < n; i++)
                {
                    temp = vv[i] * Math.Abs(LU[i, k]);
                    // Is the figure of merit for the pivot better than the best so far?
                    if (temp > big)
                    {
                        big = temp;
                        imax = i;
                    }
                }

                if (k != imax)
                {
                    // Do we need to interchange rows? Yes, do so...
                    for (j = 0; j < n; j++)
                    {
                        temp = LU[imax, j];
                        LU[imax, j] = LU[k, j];
                        LU[k, j] = temp;
                    }
                    d = -d;
                    vv[imax] = vv[k];
                }

                index[k] = imax;

                // If the pivot element is zer0, the matrix is singular (at least to the precision of the algorithm.)
                // For some applications on singular matrices, it is desirable to substitute TINY for zero.
                if (LU[k, k] == 0.0d) LU[k, k] = TINY;
                
                //Perform Gaussian Elimination step
                for (i = k + 1; i < n; i++)
                {
                    LU[i, k] /= LU[k, k]; // True Lower triangular matrix
                    for (j = k + 1; j < n; j++)
                        LU[i, j] -= LU[i, k] * LU[k, j]; // True Upper triangular matrix
                }
            }
        }

        private readonly int n; // Number of rows in A
        private readonly double d; // used by determinant routine
        private readonly int[] index; // stores the permutation

        /// <summary>
        /// Stores the decomposition.
        /// </summary>
        public Matrix LU { get; private set; }

        /// <summary>
        /// Stores the input matrix A that was LU decomposed.
        /// </summary>
        public Matrix A { get; private set; }

        /// <summary>
        /// Solves the set of n linear equations A*x=b using the stored LU decomposition of A.
        /// </summary>
        /// <param name="B">Right-hand side vector b [0..n-1].</param>   
        public Vector Solve(Vector B)
        {
            // Solves the set of n linear equations A*x=b using the stored LU decomposition of A.
            // b[0..n-1] is input as the right-hand side vector b, while x returns the solution vector x;
            // This routine takes into account the possibility that b will begin with many zero elements, 
            // so it is efficient for use in matrix inversion.
            int i, ii = 0, ip, j;
            double sum;
            var x = new Vector(n);
            if (B.Length != n)
            {
                throw new ArgumentOutOfRangeException(nameof(B), "The vector b must have the same number of rows as the matrix A.");
            }
            for (i = 0; i < n; i++) x[i] = B[i];
            for (i = 0; i < n; i++)
            {
                ip = index[i];
                sum = x[ip];
                x[ip] = x[i];
                if (ii != 0)
                    for (j = ii - 1; j < i; j++) sum -= LU[i,j] * x[j];
                else if (sum != 0.0)
                    ii = i + 1;
                x[i] = sum;
            }
            for (i = n - 1; i >= 0; i--)
            {
                sum = x[i];
                for (j = i + 1; j < n; j++) sum -= LU[i,j] * x[j];
                x[i] = sum / LU[i,i];
            }
            return x;
        }

        /// <summary>
        /// Solves m sets of n linear equations A*X=B using LU decomposition of A.
        /// </summary>
        /// <param name="B">The right-hand side matrix B [0..n-1][0..n-1].</param>
        public Matrix Solve(Matrix B)
        {
            int i, j, m = B.NumberOfColumns;
            var x = new Matrix(n);
            var xx = new Vector(n);
            if (B.NumberOfRows != n)
            {
                throw new ArgumentOutOfRangeException(nameof(B), "The vector b must have the same number of rows as the matrix A.");
            }
            for (j = 0; j < m; j++)
            {
                for (i = 0; i < n; i++) { xx[i] = B[i, j]; }
                xx = Solve(xx);
                for (i = 0; i < n; i++) { x[i, j] = xx[i]; }
            }
            return x;
        }

        /// <summary>
        /// Returns the matrix inverse A^-1 using the stored LU decomposition.
        /// </summary>
        public Matrix InverseA()
        {
            int i, j;
            var Ainv = new Matrix(n);
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++) Ainv[i, j] = 0d;
                Ainv[i, i] = 1d;
            }
            return Solve(Ainv);
        }

        /// <summary>
        /// Using the stored LU decomposition, returns the determinant of the matrix A.
        /// </summary>
        public double Determinant()
        {
            double dd = d;
            for (int i = 0; i < n; i++) 
                dd *= LU[i, i];
            return dd;
        }

    }
}