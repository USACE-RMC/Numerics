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
    /// A class for solving a set of linear equations using Singular Value Decomposition.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///<list type = "bullet" >
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// Singular value decomposition (SVD) is one of the most useful matrix decomposition methods. Used in:
    /// <list type="bullet">
    /// <item> Reducing dimension of large data sets (PCA) </item>
    /// <item>Separating "signal" from "noise" in real-world data </item>
    /// <item>Identifying line / low-dim plane of best-fit  to high dimensional data </item>
    /// <item>Solving/interpreting linear systems of equations </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017. </item>
    /// <item> <description> 
    /// <see href = "https://en.wikipedia.org/wiki/Singular_value_decomposition" />
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class SingularValueDecomposition
    {

        /// <summary>
        /// Constructs a Singular Value Decomposition. 
        /// </summary>
        /// <param name="A">The input matrix A that is to be SV decomposed.</param>
        public SingularValueDecomposition(Matrix A)
        {
            this.A = A;
            n = A.NumberOfColumns;
            m = A.NumberOfRows;
            U = new Matrix(A.ToArray());
            V = new Matrix(n, n);
            W = new Vector(n);
            Decompose();
            Reorder();
            Threshold = 0.5 * Math.Sqrt(m + n + 1d) * W[0] * eps; // Default threshold for nonzero singular values.
        }

        private readonly double eps = Tools.DoubleMachineEpsilon;      
        private readonly int n; // Number of columns in A
        private readonly int m; // Number of rows in A

        /// <summary>
        /// The default threshold value.
        /// </summary>
        public double Threshold { get; private set; } 

        /// <summary>
        /// Stores the input matrix A that was SV decomposed.
        /// </summary>
        public Matrix A { get; private set; }

        /// <summary>
        /// M x N column-orthogonal matrix U.
        /// </summary>
        public Matrix U { get; private set; }

        /// <summary>
        /// Transpose of an N x N orthogonal matrix V.
        /// </summary>
        public Matrix V { get; private set; }

        /// <summary>
        /// Diagonal matrix W, stored as a vector.
        /// </summary>
        public Vector W { get; private set; }

        /// <summary>
        /// Returns the reciprocal of the condition number of A.
        /// </summary>
        public double InverseCondition {
            get { return (W[0] <= 0d || W[n - 1] <= 0d) ? 0d : W[n - 1] / W[0]; }
        }

        /// <summary>
        /// Return the rank of A, after zeroing any singular values smaller than the threshold.  
        /// </summary>
        /// <param name="threshold">The threshold to evaluate.
        /// If the threshold is negative, a default value based on estimated roundoff is used.</param>
        public int Rank(double threshold = -1)
        {
            int j, nr = 0;
            int n = A.NumberOfColumns;
            int m = A.NumberOfRows;
            Threshold = (threshold >= 0d ? threshold : 0.5 * Math.Sqrt(m + n + 1d) * W[0] * eps);
            for (j = 0; j < n; j++) if (W[j] > Threshold) nr++;
            return nr;
        }

        /// <summary>
        /// Return the nullity of A, after zeroing any singular values smaller than the threshold. 
        /// </summary>
        /// <param name="threshold">The threshold to evaluate.
        /// If the threshold is negative, a default value based on estimated roundoff is used.</param>
        public int Nullity(double threshold = -1)
        {
            int j, nn = 0;
            Threshold = (threshold >= 0d ? threshold : 0.5 * Math.Sqrt(m + n + 1d) * W[0] * eps);
            for (j = 0; j < n; j++) if (W[j] <= Threshold) nn++;
            return nn;
        }

        /// <summary>
        /// Gives an orthonormal basis for the range of A as the columns of a return matrix. 
        /// </summary>
        /// <param name="threshold">The threshold to evaluate.
        /// If the threshold is negative, a default value based on estimated roundoff is used.</param>
        public Matrix Range(double threshold = -1)
        {
            int i, j, nr = 0;
            var rnge = new Matrix(m, Rank(threshold));
            for (j = 0; j < n; j++)
            {
                if (W[j] > Threshold)
                {
                    for (i = 0; i < m; i++) rnge[i,nr] = U[i,j];
                    nr++;
                }
            }
            return rnge;
        }

        /// <summary>
        /// Gives an orthonormal basis for the nullspace of A as the columns of a return matrix. 
        /// </summary>
        /// <param name="threshold">The threshold to evaluate.
        /// If the threshold is negative, a default value based on estimated roundoff is used.</param>
        public Matrix Nullspace(double threshold = -1)
        {
            int j, jj, nn = 0;
            var nullsp = new Matrix(n, Nullity(threshold));
            for (j = 0; j < n; j++)
            {
                if (W[j] <= Threshold)
                {
                    for (jj = 0; jj < n; jj++) nullsp[jj,nn] = V[jj,j];
                    nn++;
                }
            }
            return nullsp;
        }

        /// <summary>
        /// Solves A*x=b for a vector x using the pseudoinverse of A as obtained by SVD.
        /// </summary>
        /// <param name="B">Right-hand side vector b [0..n-1].</param>
        /// <param name="threshold">The threshold to evaluate. 
        /// If positive, the threshold is the threshold value below which singular values are considered zero. 
        /// If threshold is negative, a default based on expected roundoff error is used. 
        /// </param>
        public Vector Solve(Vector B, double threshold = -1)
        {
            // Solve A*x=b for a vector x using the pseudoinverse of A as obtained by SVD. If positive,
            // the threshold is the threshold value below which singular values are considered zero. If threshold is
            // negative, a default based on expected roundoff error is used. 
            int i, j, jj;
            double s;
            var x = new Vector(n);
            var tmp = new Vector(n);
            if (B.Length != m)
            {
                throw new ArgumentOutOfRangeException(nameof(B), "The vector b must have the same number of rows as the matrix A.");
            }
            Threshold = (threshold >= 0d ? threshold : 0.5 * Math.Sqrt(m + n + 1d) * W[0] * eps);
            for (j = 0; j < n; j++) // Calculate U^T*B
            {
                s = 0.0;
                if (W[j] > Threshold) // Nonzero result only if W-j is nonzero.
                {
                    for (i = 0; i < m; i++) s += U[i,j] * B[i];
                    s /= W[j]; // This s is then divided by W-j.
                }
                tmp[j] = s;
            }
            // Matrix multiply by V to get answer.
            for (j = 0; j < n; j++)
            {
                s = 0.0;
                for (jj = 0; jj < n; jj++) s += V[j,jj] * tmp[jj];
                x[j] = s;
            }
            return x;
        }

        /// <summary>
        /// Solve m sets of n equations A*X=B using the pseudoinverse of A. 
        /// </summary>
        /// <param name="B">The right-hand side matrix B [0..n-1][0..m-1].</param>
        /// <param name="threshold">The threshold to evaluate. 
        /// If positive, the threshold is the threshold value below which singular values are considered zero. 
        /// If threshold is negative, a default based on expected roundoff error is used. 
        /// </param>
        public Matrix Solve(Matrix B, double threshold = -1)
        {
            int i, j, p = B.NumberOfColumns;
            var x = new Matrix(n, p);
            var xx = new Vector(n);
            var bcol = new Vector(m);
            if (B.NumberOfRows != n)
            {
                throw new ArgumentOutOfRangeException(nameof(B), "The matrix B must have the same number of rows as columns in the matrix A.");
            }
            // Copy and solve each column in turn.
            for (j = 0; j < p; j++)
            {
                for (i = 0; i < m; i++) bcol[i] = B[i, j];
                xx = Solve(bcol, threshold);
                for (i = 0; i < n; i++) x[i, j] = xx[i];
            }
            return x;
        }

        /// <summary>
        /// Performs the singular value decomposition.
        /// </summary>
        private void Decompose()
        {
            bool flag;
            int i, its, j, jj, k, l=0, nm=0;
            double anorm, c, f, g, h, s, scale, x, y, z;
            var rv1 = new Vector(n);
            g = scale = anorm = 0.0;
            for (i = 0; i < n; i++)
            {
                l = i + 2;
                rv1[i] = scale * g;
                g = s = scale = 0.0;
                if (i < m)
                {
                    for (k = i; k < m; k++) scale += Math.Abs(U[k,i]);
                    if (scale != 0.0)
                    {
                        for (k = i; k < m; k++)
                        {
                            U[k,i] /= scale;
                            s += U[k,i] * U[k,i];
                        }
                        f = U[i,i];
                        g = -Tools.Sign(Math.Sqrt(s), f);
                        h = f * g - s;
                        U[i,i] = f - g;
                        for (j = l - 1; j < n; j++)
                        {
                            for (s = 0.0, k = i; k < m; k++) s += U[k,i] * U[k,j];
                            f = s / h;
                            for (k = i; k < m; k++) U[k,j] += f * U[k,i];
                        }
                        for (k = i; k < m; k++) U[k,i] *= scale;
                    }
                }
                W[i] = scale * g;
                g = s = scale = 0.0;
                if (i + 1 <= m && i + 1 != n)
                {
                    for (k = l - 1; k < n; k++) scale += Math.Abs(U[i,k]);
                    if (scale != 0.0)
                    {
                        for (k = l - 1; k < n; k++)
                        {
                            U[i,k] /= scale;
                            s += U[i,k] * U[i,k];
                        }
                        f = U[i,l - 1];
                        g = -Tools.Sign(Math.Sqrt(s), f);
                        h = f * g - s;
                        U[i,l - 1] = f - g;
                        for (k = l - 1; k < n; k++) rv1[k] = U[i,k] / h;
                        for (j = l - 1; j < m; j++)
                        {
                            for (s = 0.0, k = l - 1; k < n; k++) s += U[j,k] * U[i,k];
                            for (k = l - 1; k < n; k++) U[j,k] += s * rv1[k];
                        }
                        for (k = l - 1; k < n; k++) U[i,k] *= scale;
                    }
                }
                anorm = Math.Max(anorm, (Math.Abs(W[i]) + Math.Abs(rv1[i])));
            }
            for (i = n - 1; i >= 0; i--)
            {
                if (i < n - 1)
                {
                    if (g != 0.0)
                    {
                        for (j = l; j < n; j++)
                            V[j,i] = (U[i,j] / U[i,l]) / g;
                        for (j = l; j < n; j++)
                        {
                            for (s = 0.0, k = l; k < n; k++) s += U[i,k] * V[k,j];
                            for (k = l; k < n; k++) V[k,j] += s * V[k,i];
                        }
                    }
                    for (j = l; j < n; j++) V[i,j] = V[j,i] = 0.0;
                }
                V[i,i] = 1.0;
                g = rv1[i];
                l = i;
            }
            for (i = Math.Min(m, n) - 1; i >= 0; i--)
            {
                l = i + 1;
                g = W[i];
                for (j = l; j < n; j++) U[i,j] = 0.0;
                if (g != 0.0)
                {
                    g = 1.0 / g;
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = l; k < m; k++) s += U[k,i] * U[k,j];
                        f = (s / U[i,i]) * g;
                        for (k = i; k < m; k++) U[k,j] += f * U[k,i];
                    }
                    for (j = i; j < m; j++) U[j,i] *= g;
                }
                else for (j = i; j < m; j++) U[j,i] = 0.0;
                ++U[i,i];
            }
            for (k = n - 1; k >= 0; k--)
            {
                for (its = 0; its < 30; its++)
                {
                    flag = true;
                    for (l = k; l >= 0; l--)
                    {
                        nm = l - 1;
                        if (l == 0 || Math.Abs(rv1[l]) <= eps * anorm)
                        {
                            flag = false;
                            break;
                        }
                        if (Math.Abs(W[nm]) <= eps * anorm) break;
                    }
                    if (flag)
                    {
                        c = 0.0;
                        s = 1.0;
                        for (i = l; i < k + 1; i++)
                        {
                            f = s * rv1[i];
                            rv1[i] = c * rv1[i];
                            if (Math.Abs(f) <= eps * anorm) break;
                            g = W[i];
                            h = Pythag(f, g);
                            W[i] = h;
                            h = 1.0 / h;
                            c = g * h;
                            s = -f * h;
                            for (j = 0; j < m; j++)
                            {
                                y = U[j,nm];
                                z = U[j,i];
                                U[j,nm] = y * c + z * s;
                                U[j,i] = z * c - y * s;
                            }
                        }
                    }
                    z = W[k];
                    if (l == k)
                    {
                        if (z < 0.0)
                        {
                            W[k] = -z;
                            for (j = 0; j < n; j++) V[j,k] = -V[j,k];
                        }
                        break;
                    }
                    if (its == 99) throw new ArgumentException("There was no convergence in 100 iterations");
                    x = W[l];
                    nm = k - 1;
                    y = W[nm];
                    g = rv1[nm];
                    h = rv1[k];
                    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
                    g = Pythag(f, 1.0);
                    f = ((x - z) * (x + z) + h * ((y / (f + Tools.Sign(g, f))) - h)) / x;
                    c = s = 1.0;
                    for (j = l; j <= nm; j++)
                    {
                        i = j + 1;
                        g = rv1[i];
                        y = W[i];
                        h = s * g;
                        g = c * g;
                        z = Pythag(f, h);
                        rv1[j] = z;
                        c = f / z;
                        s = h / z;
                        f = x * c + g * s;
                        g = g * c - x * s;
                        h = y * s;
                        y *= c;
                        for (jj = 0; jj < n; jj++)
                        {
                            x = V[jj,j];
                            z = V[jj,i];
                            V[jj,j] = x * c + z * s;
                            V[jj,i] = z * c - x * s;
                        }
                        z = Pythag(f, h);
                        W[j] = z;
                        if (z != 0.0)
                        {
                            z = 1.0 / z;
                            c = f * z;
                            s = h * z;
                        }
                        f = c * g + s * y;
                        x = c * y - s * g;
                        for (jj = 0; jj < m; jj++)
                        {
                            y = U[jj,j];
                            z = U[jj,i];
                            U[jj,j] = y * c + z * s;
                            U[jj,i] = z * c - y * s;
                        }
                    }
                    rv1[l] = 0.0;
                    rv1[k] = f;
                    W[k] = x;
                }
            }
        }

        /// <summary>
        /// Axillary function to reorder.
        /// </summary>
        private void Reorder()
        {
            int i, j, k, s, inc = 1;
            double sw;
            var su = new Vector(m);
            var sv = new Vector(n);
            do { inc *= 3; inc++; } while (inc <= n);
            do
            {
                inc /= 3;
                for (i = inc; i < n; i++)
                {
                    sw = W[i];
                    for (k = 0; k < m; k++) su[k] = U[k,i];
                    for (k = 0; k < n; k++) sv[k] = V[k,i];
                    j = i;
                    while (W[j - inc] < sw)
                    {
                        W[j] = W[j - inc];
                        for (k = 0; k < m; k++) U[k,j] = U[k,j - inc];
                        for (k = 0; k < n; k++) V[k,j] = V[k,j - inc];
                        j -= inc;
                        if (j < inc) break;
                    }
                    W[j] = sw;
                    for (k = 0; k < m; k++) U[k,j] = su[k];
                    for (k = 0; k < n; k++) V[k,j] = sv[k];

                }
            } while (inc > 1);
            for (k = 0; k < n; k++)
            {
                s = 0;
                for (i = 0; i < m; i++) if (U[i,k] < 0d) s++;
                for (j = 0; j < n; j++) if (V[j,k] < 0d) s++;
                if (s > (m + n) / 2)
                {
                    for (i = 0; i < m; i++) U[i,k] = -U[i,k];
                    for (j = 0; j < n; j++) V[j,k] = -V[j,k];
                }
            }
        }

        /// <summary>
        /// Computes (a^2 + b^2)^1/2 without destructive underflow or overflow.
        /// </summary>
        private double Pythag(double a, double b) {
            double absa = Math.Abs(a), absb = Math.Abs(b);
            return (absa > absb ? absa * Math.Sqrt(1.0 + Math.Pow(absb / absa, 2)) :
                (absb == 0.0 ? 0.0 : absb * Math.Sqrt(1.0 + Math.Pow(absa / absb, 2))));
        }

        /// <summary>
        /// Takes Log determinant of the Matrix W
        /// </summary>
        public double LogDeterminant()
        {
            double det = 0;
            for (int i = 0; i < W.Length; i++)
                det += Math.Log((double)W[i]);
            return det;
        }

        /// <summary>
        /// Takes Log determinant of the Matrix W
        /// </summary>
        public double LogPseudoDeterminant()
        {
            double det = 0;
            for (int i = 0; i < W.Length; i++)
                if (W[i] != 0) det += Math.Log((double)W[i]);
            return det;
        }

    }
}