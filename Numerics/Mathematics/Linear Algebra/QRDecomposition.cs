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
    /// A class for performing QR Decomposition using Householder reflections.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b> Authors: </b>
    ///    <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///    </list>
    /// </para>
    /// <para>
    /// This class implements the QR decomposition of a general real matrix A = QR
    /// using Householder reflections.
    /// </para>
    /// <para>
    /// <b>References:</b>
    /// <list type="bullet">
    /// <item>"Numerical Recipes: The Art of Scientific Computing, Third Edition." Press et al. 2017.</item>
    /// <item><see href="https://en.wikipedia.org/wiki/QR_decomposition" /></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class QRDecomposition
    {
        private readonly int m, n;
        private readonly Matrix _Qt, _R;

        /// <summary>
        /// Constructs a new QR decomposition of matrix A.
        /// </summary>
        /// <param name="A">The real m-by-n matrix A.</param>
        public QRDecomposition(Matrix A)
        {
            m = A.NumberOfRows;
            n = A.NumberOfColumns;
            _R = A.Clone();
            _Qt = Matrix.Identity(m);

            for (int k = 0; k < Math.Min(m, n); k++)
            {
                double normX = 0.0;
                for (int i = k; i < m; i++)
                    normX += _R[i, k] * _R[i, k];
                normX = Math.Sqrt(normX);
                if (normX == 0.0) continue;

                if (_R[k, k] >= 0.0)
                    normX = -normX;

                double[] v = new double[m];
                for (int i = 0; i < m; i++) v[i] = 0.0;
                v[k] = _R[k, k] - normX;
                for (int i = k + 1; i < m; i++)
                    v[i] = _R[i, k];

                double beta = 0.0;
                for (int i = k; i < m; i++)
                    beta += v[i] * v[i];
                beta = 2.0 / beta;

                // Apply to R
                for (int j = k; j < n; j++)
                {
                    double s = 0.0;
                    for (int i = k; i < m; i++)
                        s += v[i] * _R[i, j];
                    s *= beta;
                    for (int i = k; i < m; i++)
                        _R[i, j] -= s * v[i];
                }

                // Apply to Q
                for (int j = 0; j < m; j++)
                {
                    double s = 0.0;
                    for (int i = k; i < m; i++)
                        s += v[i] * _Qt[i, j];
                    s *= beta;
                    for (int i = k; i < m; i++)
                        _Qt[i, j] -= s * v[i];
                }
            }
        }

        /// <summary>
        /// Gets the orthogonal matrix Q.
        /// </summary>
        public Matrix Q => Matrix.Transpose(_Qt);

        /// <summary>
        /// Gets the upper triangular matrix R.
        /// </summary>
        public Matrix RMatrix => _R;

        /// <summary>
        /// Solves Ax = b using QR decomposition.
        /// </summary>
        /// <param name="b">The right-hand side vector.</param>
        public Vector Solve(Vector b)
        {
            if (b.Length != m)
                throw new ArgumentException("Vector b must match the number of rows in A.");

            int min = Math.Min(m, n);
            var Qtb = _Qt * b;
            var x = new double[n];
            for (int i = min - 1; i >= 0; i--)
            {
                x[i] = Qtb[i];
                for (int j = i + 1; j < min; j++)
                    x[i] -= _R[i, j] * x[j];
                x[i] /= _R[i, i];
            }
            return new Vector(x);
        }

        /// <summary>
        /// Solves AX = B using QR decomposition where B is a matrix.
        /// </summary>
        /// <param name="B">The right-hand side matrix.</param>
        public Matrix Solve(Matrix B)
        {
            if (B.NumberOfRows != m)
                throw new ArgumentException("Matrix B must have the same number of rows as A.");

            int min = Math.Min(m, n);
            var X = new Matrix(n, B.NumberOfColumns);
            for (int j = 0; j < B.NumberOfColumns; j++)
            {
                var col = new Vector(B.Column(j));
                var Qtb = _Qt * col;
                var xCol = new double[n];
                for (int i = min - 1; i >= 0; i--)
                {
                    xCol[i] = Qtb[i];
                    for (int k = i + 1; k < min; k++)
                        xCol[i] -= _R[i, k] * xCol[k];
                    xCol[i] /= _R[i, i];
                }
                for (int i = 0; i < min; i++)
                    X[i, j] = xCol[i];
            }
            return X;
        }
    }
}
