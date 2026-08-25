using System;
using System.Globalization;


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
        /// Constructs new Cholesky Decomposition using the default scale-relative pivot tolerance.
        /// </summary>
        /// <param name="A">The positive-definite symmetric input matrix A [0..n-1][0..n-1] that is to be Cholesky decomposed.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the matrix A is not square.</exception>
        /// <exception cref="Exception">Thrown when the matrix A is not positive-definite.</exception>
        /// <remarks>
        /// <para>
        /// The pivot tolerance is <see cref="DefaultRelativeTolerance(int)"/> evaluated at the dimension of A.
        /// </para>
        /// <para>
        /// <b>Success is not a rank certificate.</b> The pivot test rejects most numerically rank-deficient
        /// matrices but not all of them — for a covariance estimated from <c>m = n - 1</c> observations, 11% to
        /// 18% still factorize — so <see cref="IsPositiveDefinite"/> being true does not prove full rank. Use a
        /// singular value decomposition when the rank genuinely has to be known.
        /// </para>
        /// </remarks>
        public CholeskyDecomposition(Matrix A)
            : this(A, DefaultRelativeTolerance(A.NumberOfRows))
        {
        }

        /// <summary>
        /// Constructs new Cholesky Decomposition with an explicit scale-relative pivot tolerance.
        /// </summary>
        /// <param name="A">The positive-definite symmetric input matrix A [0..n-1][0..n-1] that is to be Cholesky decomposed.</param>
        /// <param name="relativeTolerance">
        /// The pivot tolerance, expressed as a fraction of the corresponding diagonal entry of A. A pivot is
        /// rejected when it falls at or below <c>relativeTolerance * A[i,i]</c>. Pass zero to reproduce the
        /// purely absolute <c>pivot &lt;= 0</c> test exactly.
        /// </param>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when the matrix A is not square, or when <paramref name="relativeTolerance"/> is not a finite
        /// value in the interval [0, 1).
        /// </exception>
        /// <exception cref="Exception">Thrown when the matrix A is not positive-definite.</exception>
        /// <remarks>
        /// <para>
        /// The pivot at step i is the conditional variance of variable i given variables 0..i-1, so it is
        /// naturally measured against <c>A[i,i]</c> — the unconditional variance of that same variable — rather
        /// than against the largest diagonal of A. A test relative to the largest diagonal would falsely reject
        /// a covariance that legitimately mixes a very small variance with a very large one, because the small
        /// variable's pivot is small in absolute terms while being a perfectly healthy fraction of its own
        /// diagonal.
        /// </para>
        /// <para>
        /// The scale-relative test exists because an exactly rank-deficient matrix does not generally produce a
        /// non-positive pivot in floating point. Rounding leaves a small positive residue instead, the
        /// factorization completes, and the matrix is silently reported positive-definite with a wildly wrong
        /// determinant and inverse. For example, the exactly rank-two covariance
        /// <c>[[2, 0.5, 2], [0.5, 1, 0.5], [2, 0.5, 2]]</c> yields a final pivot of 4.44E-16 rather than zero;
        /// under the absolute test it factorizes, and its log determinant comes out near -34.8 instead of the
        /// log pseudo-determinant 1.2528. See <see href="https://github.com/USACE-RMC/Numerics/issues/145"/>.
        /// </para>
        /// <para>
        /// <b>This test is not a rank certificate.</b> It rejects most numerically rank-deficient matrices but
        /// not all of them, so <see cref="IsPositiveDefinite"/> being true does not prove the matrix has full
        /// rank. Measured over 1,000 trials per dimension on a covariance estimated from <c>m = n - 1</c>
        /// observations of <c>n</c> variables — exactly rank <c>m</c>, and the commonest way rank deficiency
        /// arises in practice — the test removes roughly two thirds of the matrices the absolute test had
        /// accepted, while 11% to 18% still factorize. When the rank of a matrix genuinely has to be known,
        /// use a singular value decomposition; for a multivariate normal that is
        /// <c>DecompositionMethod.SingularValue</c>, which is the only reliable rank test in this library.
        /// </para>
        /// <para>
        /// When <c>A[i,i]</c> is not a positive finite number the threshold falls back to zero, which is the
        /// absolute test. That case cannot weaken the result: a positive-definite matrix has a strictly positive
        /// finite diagonal, so a non-positive or non-finite diagonal is rejected on its own merits.
        /// </para>
        /// </remarks>
        public CholeskyDecomposition(Matrix A, double relativeTolerance)
        {

            IsPositiveDefinite = false;
            int i, j, k;
            if (double.IsNaN(relativeTolerance) || double.IsInfinity(relativeTolerance) || relativeTolerance < 0d || relativeTolerance >= 1d)
            {
                throw new ArgumentOutOfRangeException(nameof(relativeTolerance), "The relative tolerance must be a finite value in the interval [0, 1).");
            }
            RelativeTolerance = relativeTolerance;
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
                        // Reject a pivot that is negligible relative to its own diagonal entry. The diagonal
                        // guard keeps the threshold at zero — today's absolute test — whenever A[i,i] is not a
                        // positive finite number.
                        double diagonal = this.A[i, i];
                        double threshold = diagonal > 0d && !double.IsInfinity(diagonal) ? relativeTolerance * diagonal : 0d;
                        if (double.IsNaN(sum) || sum <= 0d)
                            throw new Exception("Cholesky Decomposition failed. The input matrix is not positive-definite.");
                        if (sum <= threshold)
                            throw new Exception("Cholesky Decomposition failed. The input matrix is not positive-definite. The pivot at row "
                                + i.ToString(CultureInfo.InvariantCulture) + " is "
                                + (sum / diagonal).ToString("E6", CultureInfo.InvariantCulture)
                                + " times its diagonal entry, at or below the relative tolerance "
                                + relativeTolerance.ToString("E6", CultureInfo.InvariantCulture)
                                + ", so the matrix is numerically rank-deficient.");
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

        /// <summary>
        /// Returns the default scale-relative pivot tolerance for a matrix of the given dimension.
        /// </summary>
        /// <param name="dimension">The number of rows in the matrix to be decomposed.</param>
        /// <returns>
        /// The tolerance, approximately <c>dimension * 2^-52</c>, or zero when the dimension is not positive.
        /// </returns>
        /// <remarks>
        /// <para>
        /// The value is approximately <c>n * 2^-52 ≈ n * 2.22E-16</c>, computed as <c>2 * n *</c>
        /// <see cref="Tools.DoubleMachineEpsilon"/> because that constant is the unit roundoff <c>2^-53</c>
        /// rather than the double-precision spacing <c>2^-52</c>. The library constant is a decimal-truncated
        /// <c>2^-53</c>, so the computed tolerance exceeds <c>n * 2^-52</c> by about 3.1E-15 relative.
        /// </para>
        /// <para>
        /// This tracks the standard backward-error bound for Cholesky factorization, in which the computed pivot
        /// differs from the exact one by at most <c>γ_i * A[i,i]</c> — of order <c>n</c> unit roundoffs times
        /// the corresponding diagonal entry. A pivot at or below that level carries no information beyond
        /// rounding noise. Because the noise floor scales with <c>A[i,i]</c> and not with the largest diagonal,
        /// the test is invariant under the rescaling <c>D*A*D</c> for a positive diagonal <c>D</c>, which is the
        /// correct invariance for a covariance: changing the units of one variable must not change whether the
        /// matrix is accepted.
        /// </para>
        /// <para>
        /// The margin against falsely rejecting a legitimate matrix is large. For two variables with correlation
        /// ρ the pivot ratio is <c>1 - ρ²</c>, so a false rejection at <c>n = 2</c> requires
        /// <c>1 - ρ² &lt;= 4.44E-16</c>, that is ρ within about two ulp of one. Measured across the Numerics
        /// test suite (29.8 million factorizations) the smallest ratio produced by a genuinely positive-definite
        /// matrix is 2.0E-12, some 4,500 times the tolerance that applies to it.
        /// </para>
        /// <para>
        /// The margin in the other direction is far smaller, and no tolerance of this form can close it. A
        /// rank-deficient matrix whose accumulated rounding leaves a pivot ratio above the tolerance is still
        /// accepted, and that is not rare: for a covariance estimated from <c>m = n - 1</c> observations of
        /// <c>n</c> variables, 11% to 18% of trials still factorize. Raising the tolerance would eat into the
        /// false-rejection margin without fixing this; determining rank requires a singular value decomposition.
        /// </para>
        /// </remarks>
        public static double DefaultRelativeTolerance(int dimension)
        {
            return dimension > 0 ? 2d * dimension * Tools.DoubleMachineEpsilon : 0d;
        }

        private readonly int n; // Number of rows in A

        /// <summary>
        /// The scale-relative pivot tolerance applied during the factorization.
        /// </summary>
        /// <remarks>
        /// A pivot was rejected when it fell at or below this fraction of the corresponding diagonal entry of
        /// <see cref="A"/>. Zero means the purely absolute test was used.
        /// </remarks>
        public double RelativeTolerance { get; private set; }

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