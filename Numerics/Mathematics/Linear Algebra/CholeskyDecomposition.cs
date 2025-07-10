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
    /// A class for solving a set of linear equations using Cholesky Decomposition.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b> Authors: </b>
    /// <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The Cholesky decomposition or Cholesky factorization is a decomposition of a Hermitian, positive-definite matrix into the product of 
    /// a lower triangular matrix and its conjugate transpose. The basic idea of this method includes decomposing a square positive-definite
    /// matrix by symmetrically applying column-clearing operations from Gaussian Elimination. This method is very useful for Monte Carlo 
    /// simulations that is utilized in TotalRisk. This method is also roughly twice as efficient as the LU decomposition for solving 
    /// systems of linear equations. 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes: The art of Scientific Computing, Third Edition." Press et al. 2017. </item>
    /// <item>"Numerical Methods for Engineers, Second Edition.", D.V. Griffiths and I.M. Smith, Taylor and Francis Group, 2006. </item>
    /// <item><description> 
    /// <see href = "https://en.wikipedia.org/wiki/Cholesky_decomposition" /> 
    /// </description></item>
    /// <item><description> 
    ///  <see href = "https://www.geeksforgeeks.org/cholesky-decomposition-matrix-decomposition/"/> 
    /// </description></item> 
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class CholeskyDecomposition
    {
     
        /// <summary>
        /// Constructs new Cholesky Decomposition.
        /// </summary>
        /// <param name="A">The positive-definite symmetric input matrix A [0..n-1][0..n-1] that is to be Cholesky decomposed.</param>
        public CholeskyDecomposition(Matrix A)
        {

            IsPositiveDefinite = false;
            int i, j, k;
            n = A.NumberOfRows;
            this.A = new Matrix(A.ToArray());
            L = new Matrix(A.ToArray()); // Lower triangular matrix
            double sum;
            if (A.NumberOfColumns != A.NumberOfRows)
            {
                throw new ArgumentOutOfRangeException(nameof(A), "The matrix A must be square.");
            }
            
            //Decomposing a matrix into Lower triangular
            for (i = 0; i < n; i++)
            {
                for (j = i; j < n; j++)
                {
                    sum = L[i, j];

                    for (k = i - 1; k >= 0; k -= 1)
                        sum -= L[i, k] * L[j, k]; // Cholesky formula 
                    if (i == j)
                    {
                        if (double.IsNaN(sum) || sum <= 0d)
                            throw new Exception("Cholesky Decomposition failed. The input matrix is not positive-definite.");
                        L[i, i] = Math.Sqrt(sum);
                    }
                    else
                    {
                        L[j, i] = sum / L[i, i]; // Upper Triangular matrix
                    }
                }
            }
            
            // Making sure 0 entries for upper triangular matrix
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < i; j++)
                    L[j, i] = 0.0d;
            }
            // Failure of the decomposition indicates that the matrix A is not positive-definite. 
            // Success, means it is. 
            IsPositiveDefinite = true;
        }

        private readonly int n; // Number of rows in A

        /// <summary>
        /// Stores the decomposition.
        /// </summary>
        public Matrix L { get; private set; }

        /// <summary>
        /// Stores the input matrix A that was decomposed.
        /// </summary>
        public Matrix A { get; private set; }

        /// <summary>
        /// Determines whether the input matrix A is positive definite.
        /// </summary>
        public bool IsPositiveDefinite { get; private set; }

        /// <summary>
        /// Solves the set of n linear equations A*x=b using the stored Cholesky decomposition of A=L*L^T.
        /// </summary>
        /// <returns> x vector from A*x=b </returns>
        /// <param name="b">Right-hand side vector b [0..n-1].</param>
        public Vector Solve(Vector b)
        {
            int i, k;
            double sum;
            var x = new Vector(n);
            if (b.Length != n)
            {
                throw new ArgumentOutOfRangeException(nameof(b), "The vector b must have the same number of rows as the matrix A.");
            }
            for (i = 0; i < n; i++)
            {
                for (sum = b[i], k = i - 1; k >= 0; k--) sum -= L[i,k] * x[k];
                x[i] = sum / L[i,i];
            }
            for (i = n - 1; i >= 0; i--)
            {
                for (sum = x[i], k = i + 1; k < n; k++) sum -= L[k,i] * x[k];
                x[i] = sum / L[i,i];
            }
            return x;
        }

        /// <summary>
        /// Solving the L^T * x = y equation with backward substitution
        /// </summary>
        /// <param name="y">The input vector y [0..n-1].</param>

        public Vector Backward(Vector y)
        {
            int i, j;
            var x = new Vector(n);

            if (y.Length != n)
            {
                throw new ArgumentOutOfRangeException(nameof(y), "The vector y must have the same number of rows as the matrix A.");
            }

            for (i = n - 1; i >= 0; --i)
            {

                var sum = y[i];
                for (j = n - 1; j > i; --j)
                {
                    sum -= x[j] * L[j, i];
                }
                x[i] = sum / L[i, i];
            }

            return x;
        }

        /// <summary>
        /// Solving the L * y = b equation using Forward substitution
        /// </summary>
        /// <param name="b">The right-hand side vector b [0..n-1].</param>
        public Vector Forward(Vector b)
        {
            int i, j;
            var y = new Vector(n);
            double sum;
            if (b.Length != n)
            {
                throw new ArgumentOutOfRangeException(nameof(b), "The vector b must have the same number of rows as the matrix A.");
            }

            for (i = 0; i < n; i++)
            {
                sum = b[i];
                for (j = 0; j < i; j++)
                    sum -= L[i, j] * y[j];
                y[i] = sum / L[i, i];
            }

            return y;
        }

        /// <returns>
        /// Matrix inverse A^-1 using the stored Cholesky decomposition.
        /// </returns>
        public Matrix InverseA()
        {
            int i, j, k;
            var Ainv = new Matrix(n);
            double sum;
            for (i = 0; i < n; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    sum = i == j ? 1.0d : 0.0d;
                    for (k = i - 1; k >= j; k -= 1)
                        sum -= L[i, k] * Ainv[j, k];
                    Ainv[j, i] = sum / L[i, i];
                }
            }

            for (i = n - 1; i >= 0; i -= 1)
            {
                for (j = 0; j <= i; j++)
                {
                    sum = i < j ? 0.0d : Ainv[j, i];
                    for (k = i + 1; k < n; k++)
                        sum -= L[k, i] * Ainv[j, k];
                    Ainv[j, i] = sum / L[i, i];
                    Ainv[i, j] = Ainv[j, i];
                }
            }

            return Ainv;
        }

        /// <summary>
        /// Using the stored Cholesky decomposition, returns the determinant of the matrix A.
        /// </summary>
        public double Determinant()
        {
            double d = 1d;
            for (int i = 0; i < n; i++)
                d *= L[i, i];
            return Math.Pow(d, 2d);
        }

        /// <summary>
        /// Using the stored Cholesky decomposition, returns the logarithm of the determinant of the matrix A.
        /// </summary>
        public double LogDeterminant()
        {
            double sum = 0d;
            for (int i = 0; i < n; i++)
                sum += Math.Log(L[i, i]);
            return 2.0d * sum;
        }

        
    }
}