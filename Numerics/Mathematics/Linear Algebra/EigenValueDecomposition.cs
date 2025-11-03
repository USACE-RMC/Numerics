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

namespace Numerics.Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class for computing all eigenvalues and eigenvectors of a real <b>symmetric</b> matrix using the Jacobi rotation method.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b> Authors: </b>
    /// <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// Computes the eigen-decomposition A = V * D * V^T for a real symmetric matrix A, where D contains the eigenvalues and
    /// the columns of V are the corresponding orthonormal eigenvectors. The implementation uses the classic Jacobi iterative
    /// method, which is robust and accurate for small to medium-sized problems.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes: The Art of Scientific Computing, Third Edition." Press et al., 2017. </item>
    /// <item> <description><see href = "https://en.wikipedia.org/wiki/Jacobi_eigenvalue_algorithm"/></description> </item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class EigenValueDecomposition
    {
        /// <summary>
        /// Constructs a new symmetric eigen-decomposition.
        /// </summary>
        /// <param name="A">The <b>symmetric</b> input matrix A [0..n-1][0..n-1].</param>
        public EigenValueDecomposition(Matrix A)
        {
            if (A.NumberOfRows != A.NumberOfColumns)
                throw new ArgumentOutOfRangeException(nameof(A), "The matrix A must be square.");
            if (!A.IsSymmetric)
                throw new ArgumentException("The matrix A must be symmetric for this decomposition.", nameof(A));

            n = A.NumberOfRows;
            this.A = new Matrix(A.ToArray());
            EigenVectors = Matrix.Identity(n);
            EigenValues = new Vector(n);

            // Work on a local copy of A (array) for speed
            var a = this.A.Array; // same storage as this.A
            const double tol = 1e-12;
            const int maxIter = 2000;

            for (int iter = 0; iter < maxIter; iter++)
            {
                // Find largest off-diagonal element
                double max = 0.0; int p = 0, q = 0;
                for (int i = 0; i < n - 1; i++)
                {
                    for (int j = i + 1; j < n; j++)
                    {
                        double val = Math.Abs(a[i, j]);
                        if (val > max) { max = val; p = i; q = j; }
                    }
                }

                if (max <= tol) break; // Converged

                // Current largest off-diagonal pair (p, q)
                double app = a[p, p];
                double aqq = a[q, q];
                double apq = a[p, q];

                // Robust angle computation: θ = 0.5 * atan2(2*apq, aqq - app)
                double theta = 0.5 * Math.Atan2(2.0 * apq, aqq - app);
                double c = Math.Cos(theta);
                double s = Math.Sin(theta);

                // Rotate the 2x2 pivot block exactly
                double appNew = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                double aqqNew = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                a[p, p] = appNew;
                a[q, q] = aqqNew;
                a[p, q] = a[q, p] = 0.0;

                // Update the rest of the rows/cols (preserve symmetry)
                for (int k = 0; k < n; k++)
                {
                    if (k == p || k == q) continue;
                    double aik = a[k, p];
                    double akq = a[k, q];
                    double akpNew = c * aik - s * akq;
                    double akqNew = s * aik + c * akq;
                    a[k, p] = a[p, k] = akpNew;
                    a[k, q] = a[q, k] = akqNew;
                }

                // Accumulate eigenvectors: V = V * J(p,q)
                for (int k = 0; k < n; k++)
                {
                    double vip = EigenVectors[k, p];
                    double viq = EigenVectors[k, q];
                    EigenVectors[k, p] = c * vip - s * viq;
                    EigenVectors[k, q] = s * vip + c * viq;
                }
            }

            // Extract eigenvalues from diagonal of A
            for (int i = 0; i < n; i++) 
                EigenValues[i] = this.A[i, i];
        }

        private readonly int n; // Size of A

        /// <summary>
        /// Stores the input matrix A that was decomposed (copied from the constructor input).
        /// </summary>
        public Matrix A { get; private set; }

        /// <summary>
        /// Returns the vector of eigenvalues (length n).
        /// </summary>
        public Vector EigenValues { get; private set; }

        /// <summary>
        /// Returns the matrix of eigenvectors (n x n). Each column corresponds to an eigenvalue.
        /// </summary>
        public Matrix EigenVectors { get; private set; }

        /// <summary>
        /// Returns the effective sample size based on Dutilleul's method (1993). 
        /// </summary>
        public double EffectiveSampleSize()
        {
            double sum = 0;
            double sumsq = 0;
            for (int i = 0; i < EigenValues.Length; i++)
            {
                // Clip tiny negative eigenvalues that can appear from numerical error
                double lambda = EigenValues[i];
                if (lambda < 0.0 && Math.Abs(lambda) <= 1e-10) lambda = 0.0;
                sum += lambda;
                sumsq += lambda * lambda;
            }
            if (sumsq <= 1E-12) return 0.0; // degenerate case
            double neff = (sum * sum) / sumsq;
            return neff;
        }
    }
}
