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
                // Column norm scaled by its largest entry, so squares cannot overflow or underflow.
                double normX = ScaledNorm(_R, k, m, k);
                if (normX == 0.0) continue;

                if (_R[k, k] >= 0.0)
                    normX = -normX;

                double[] v = new double[m];
                for (int i = 0; i < m; i++) v[i] = 0.0;
                v[k] = _R[k, k] - normX;
                for (int i = k + 1; i < m; i++)
                    v[i] = _R[i, k];

                // Normalize the Householder vector with a scaled norm, so H = I - 2 v v^T without forming v^T v.
                double vMax = 0.0;
                for (int i = k; i < m; i++)
                    vMax = Math.Max(vMax, Math.Abs(v[i]));
                double vSum = 0.0;
                for (int i = k; i < m; i++)
                {
                    double t = v[i] / vMax;
                    vSum += t * t;
                }
                double vNorm = vMax * Math.Sqrt(vSum);
                for (int i = k; i < m; i++)
                    v[i] /= vNorm;
                double beta = 2.0;

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
        /// <returns>The solution x of length n, in the least-squares sense when A has more rows than columns.</returns>
        /// <exception cref="ArgumentException">
        /// Thrown when b does not match the rows of A, or when A is numerically rank-deficient.
        /// </exception>
        public Vector Solve(Vector b)
        {
            if (b.Length != m)
                throw new ArgumentException("Vector b must match the number of rows in A.");
            CheckFullRank();

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
        /// <returns>The solution X, n by the number of columns of B, in the least-squares sense when A has more rows
        /// than columns.</returns>
        /// <exception cref="ArgumentException">
        /// Thrown when B does not match the rows of A, or when A is numerically rank-deficient.
        /// </exception>
        public Matrix Solve(Matrix B)
        {
            if (B.NumberOfRows != m)
                throw new ArgumentException("Matrix B must have the same number of rows as A.");
            CheckFullRank();

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

        /// <summary>
        /// Returns the Euclidean norm of rows [start, end) of a column, scaled by the column's largest entry.
        /// </summary>
        /// <param name="A">The matrix.</param>
        /// <param name="start">The first row.</param>
        /// <param name="end">The row after the last.</param>
        /// <param name="column">The column.</param>
        /// <returns>The norm, or zero for an all-zero segment.</returns>
        private static double ScaledNorm(Matrix A, int start, int end, int column)
        {
            double largest = 0.0;
            for (int i = start; i < end; i++)
                largest = Math.Max(largest, Math.Abs(A[i, column]));
            if (largest == 0.0) return 0.0;
            double sum = 0.0;
            for (int i = start; i < end; i++)
            {
                double t = A[i, column] / largest;
                sum += t * t;
            }
            return largest * Math.Sqrt(sum);
        }

        /// <summary>
        /// Throws when a diagonal entry of R is negligible, so back substitution would divide by roundoff.
        /// </summary>
        /// <exception cref="ArgumentException">Thrown when A is numerically rank-deficient.</exception>
        /// <remarks>
        /// A diagonal entry is negligible when |R[i,i]| &lt;= max(m, n) * 2^-52 * max |R[j,j]|, the rank tolerance
        /// that numpy.linalg.matrix_rank applies to singular values, here applied to the diagonal of R.
        /// </remarks>
        private void CheckFullRank()
        {
            int min = Math.Min(m, n);
            double largest = 0d;
            for (int i = 0; i < min; i++)
                largest = Math.Max(largest, Math.Abs(_R[i, i]));
            double tolerance = Math.Max(m, n) * 2d * Tools.DoubleMachineEpsilon * largest;
            for (int i = 0; i < min; i++)
            {
                if (!(Math.Abs(_R[i, i]) > tolerance))
                    throw new ArgumentException("The matrix A is numerically rank-deficient, so the least-squares solution is not unique.");
            }
        }
    }
}
