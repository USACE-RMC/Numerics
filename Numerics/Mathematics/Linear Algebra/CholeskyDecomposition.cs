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
        /// Success is not a rank certificate; see <see cref="DefaultRelativeTolerance(int)"/> for the
        /// limits of the pivot test.
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
        /// naturally measured against <c>A[i,i]</c>, the unconditional variance of that same variable, rather
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
        /// The test is not a rank certificate: it rejects most numerically rank-deficient matrices but
        /// not all of them, so <see cref="IsPositiveDefinite"/> being true does not prove the matrix has
        /// full rank. See <see cref="DefaultRelativeTolerance(int)"/> for the margins in both directions.
        /// When the rank genuinely has to be known, use a singular value decomposition; for a
        /// multivariate normal that is <c>DecompositionMethod.SingularValue</c>.
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
            ValidateRelativeTolerance(relativeTolerance);
            n = A.NumberOfRows;
            this.A = new Matrix(A.Array);
            RelativeTolerance = relativeTolerance;
            if (!TryFactorizeFlat(this.A, relativeTolerance, out double[] factor, out int failedRow, out double failedPivot))
            {
                if (double.IsNaN(failedPivot) || failedPivot <= 0d)
                    throw new Exception("Cholesky Decomposition failed. The input matrix is not positive-definite.");
                throw new Exception("Cholesky Decomposition failed. The input matrix is not positive-definite. The pivot at row "
                    + failedRow.ToString(CultureInfo.InvariantCulture) + " is "
                    + (failedPivot / this.A[failedRow, failedRow]).ToString("E6", CultureInfo.InvariantCulture)
                    + " times its diagonal entry, at or below the relative tolerance "
                    + relativeTolerance.ToString("E6", CultureInfo.InvariantCulture)
                    + ", so the matrix is numerically rank-deficient.");
            }
            _factor = factor;
            L = ToMatrix(factor, n);
            IsPositiveDefinite = true;
        }

        /// <summary>
        /// Factors a symmetric matrix using the same pivot decisions as the public constructor,
        /// without throwing for an expected rejected pivot.
        /// </summary>
        /// <param name="matrix">The symmetric input matrix, which is not modified.</param>
        /// <param name="relativeTolerance">The finite pivot tolerance in [0, 1).</param>
        /// <param name="lower">The lower triangular factor on success; an incomplete factor on failure. Either way the strict
        /// upper triangle is zero.</param>
        /// <param name="failedRow">The rejected pivot's row, or -1 on success.</param>
        /// <param name="failedPivot">The rejected pivot before its square root, or NaN on success.</param>
        /// <returns>True when every pivot passes the existing positive-definiteness checks.</returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when the matrix is not square or the tolerance is outside its permitted range.
        /// </exception>
        /// <remarks>The factor of <see cref="TryFactorizeFlat"/> as a matrix.</remarks>
        internal static bool TryFactorize(Matrix matrix, double relativeTolerance, out Matrix lower,
            out int failedRow, out double failedPivot)
        {
            bool success = TryFactorizeFlat(matrix, relativeTolerance, out double[] factor, out failedRow, out failedPivot);
            lower = ToMatrix(factor, matrix.NumberOfRows);
            return success;
        }

        /// <summary>
        /// Factors a symmetric matrix into a flat row-major lower triangle, using the same pivot decisions as the
        /// public constructor, without throwing for an expected rejected pivot.
        /// </summary>
        /// <param name="matrix">The symmetric input matrix, which is not modified.</param>
        /// <param name="relativeTolerance">The finite pivot tolerance in [0, 1).</param>
        /// <param name="factor">
        /// The n-by-n row-major lower triangular factor on success; an incomplete factor on failure. Either way the strict
        /// upper triangle is zero.
        /// </param>
        /// <param name="failedRow">The rejected pivot's row, or -1 on success.</param>
        /// <param name="failedPivot">The rejected pivot before its square root, or NaN on success.</param>
        /// <returns>True when every pivot passes the existing positive-definiteness checks.</returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when the matrix is not square or the tolerance is outside its permitted range.
        /// </exception>
        /// <remarks>
        /// <para>
        /// Preserves the constructor's upper-triangle input convention and diagonal-relative threshold: the upper
        /// triangle and diagonal of the input are read and transposed into a flat lower triangle. Invalid arguments
        /// remain distinct from rejected pivots.
        /// </para>
        /// <para>
        /// From order <see cref="BlockedDimension"/> on, the factor is computed by the blocked kernel
        /// <see cref="CholeskyKernel.TryFactor"/>. Below it, the scalar loop is used, with each inner product
        /// subtracted term by term for k from i − 1 down to 0, which is faster at those orders.
        /// </para>
        /// </remarks>
        internal static bool TryFactorizeFlat(Matrix matrix, double relativeTolerance, out double[] factor,
            out int failedRow, out double failedPivot)
        {
            ValidateRelativeTolerance(relativeTolerance);
            if (matrix.NumberOfColumns != matrix.NumberOfRows)
                throw new ArgumentOutOfRangeException("A", "The matrix A must be square.");

            int dimension = matrix.NumberOfRows;
            var source = matrix.Array;
            factor = new double[dimension * dimension];
            for (int i = 0; i < dimension; i++)
            {
                int row = i * dimension;
                for (int j = 0; j <= i; j++)
                    factor[row + j] = source[j, i];
            }
            if (dimension >= BlockedDimension)
                return CholeskyKernel.TryFactor(factor, 0, dimension, dimension, relativeTolerance, null, out failedRow, out failedPivot);

            failedRow = -1;
            failedPivot = double.NaN;
            for (int i = 0; i < dimension; i++)
            {
                int rowI = i * dimension;
                for (int j = i; j < dimension; j++)
                {
                    int rowJ = j * dimension;
                    double sum = factor[rowJ + i];
                    for (int k = i - 1; k >= 0; k--)
                        sum -= factor[rowI + k] * factor[rowJ + k];
                    if (i == j)
                    {
                        double diagonal = source[i, i];
                        double threshold = diagonal > 0d && !double.IsInfinity(diagonal) ? relativeTolerance * diagonal : 0d;
                        if (double.IsNaN(sum) || sum <= 0d || sum <= threshold)
                        {
                            failedRow = i;
                            failedPivot = sum;
                            return false;
                        }
                        factor[rowI + i] = Math.Sqrt(sum);
                    }
                    else
                    {
                        factor[rowJ + i] = sum / factor[rowI + i];
                    }
                }
            }
            return true;
        }

        /// <summary>
        /// The order from which the factorization and the solves use the blocked kernels of <see cref="CholeskyKernel"/>.
        /// </summary>
        /// <remarks>
        /// Below this order the scalar loops are faster, because the kernels' call and vector-setup overhead dominates;
        /// from it on the kernels are faster, and the backward error of the factor is no larger than that of the scalar loop
        /// on random symmetric positive-definite test matrices.
        /// </remarks>
        internal const int BlockedDimension = 16;

        /// <summary>Copies a flat row-major square array into a new matrix.</summary>
        /// <param name="work">The row-major values.</param>
        /// <param name="dimension">The number of rows and columns.</param>
        /// <returns>The matrix.</returns>
        private static Matrix ToMatrix(double[] work, int dimension)
        {
            var result = new Matrix(dimension);
            var target = result.Array;
            for (int i = 0; i < dimension; i++)
            {
                int row = i * dimension;
                for (int j = 0; j < dimension; j++)
                    target[i, j] = work[row + j];
            }
            return result;
        }

        /// <summary>Validates the pivot tolerance before accessing the input matrix.</summary>
        /// <param name="relativeTolerance">The finite relative pivot tolerance in [0, 1).</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the tolerance is invalid.</exception>
        /// <remarks>Preserves the public constructor's original argument-validation order.</remarks>
        private static void ValidateRelativeTolerance(double relativeTolerance)
        {
            if (double.IsNaN(relativeTolerance) || double.IsInfinity(relativeTolerance) || relativeTolerance < 0d || relativeTolerance >= 1d)
                throw new ArgumentOutOfRangeException(nameof(relativeTolerance), "The relative tolerance must be a finite value in the interval [0, 1).");
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
        /// rather than the double-precision spacing <c>2^-52</c>. The library constant is <c>2^-53</c>
        /// rounded to fifteen significant figures, which sits just above the exact value, so the computed
        /// tolerance exceeds <c>n * 2^-52</c> by about 3.1E-15 relative.
        /// </para>
        /// <para>
        /// This tracks the standard backward-error bound for Cholesky factorization, in which the computed pivot
        /// differs from the exact one by at most <c>γ_i * A[i,i]</c>, of order <c>n</c> unit roundoffs times
        /// the corresponding diagonal entry. A pivot at or below that level carries no information beyond
        /// rounding noise. Because the noise floor scales with <c>A[i,i]</c> and not with the largest diagonal,
        /// the test is invariant under the rescaling <c>D*A*D</c> for a positive diagonal <c>D</c>, which is the
        /// correct invariance for a covariance: changing the units of one variable must not change whether the
        /// matrix is accepted.
        /// </para>
        /// <para>
        /// The margin against falsely rejecting a legitimate matrix is large. For two variables with correlation
        /// ρ the pivot ratio is <c>1 - ρ²</c>, so a false rejection at <c>n = 2</c> requires
        /// <c>1 - ρ² &lt;= 4.44E-16</c>, that is ρ within about two ulp of one.
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

        /// <summary>
        /// Creates a decomposition from a factor already computed by <see cref="TryFactorize"/>.
        /// </summary>
        /// <param name="A">The matrix that was factored; stored without copying.</param>
        /// <param name="factor">The flat row-major lower triangular factor of A.</param>
        /// <param name="relativeTolerance">The pivot tolerance used for the factorization.</param>
        private CholeskyDecomposition(Matrix A, double[] factor, double relativeTolerance)
        {
            n = A.NumberOfRows;
            this.A = A;
            _factor = factor;
            L = ToMatrix(factor, n);
            RelativeTolerance = relativeTolerance;
            IsPositiveDefinite = true;
        }

        /// <summary>
        /// Attempts a Cholesky decomposition, retrying with an escalating diagonal jitter when the matrix is
        /// not numerically positive-definite. Uses relative jitters from 1E-10 to 1E-2.
        /// </summary>
        /// <param name="A">The symmetric input matrix A [0..n-1][0..n-1], which is not modified.</param>
        /// <param name="decomposition">
        /// The decomposition of <c>A + jitter * I</c> on success; otherwise <see langword="null"/>.
        /// </param>
        /// <param name="jitter">
        /// The absolute diagonal jitter that was added: zero when A factors as given, or NaN on failure.
        /// </param>
        /// <returns><see langword="true"/> when A or a jittered A was factored; otherwise <see langword="false"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when A is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when A is not square.</exception>
        /// <remarks>
        /// Equivalent to <see cref="TryCholesky(Matrix, double, double, out CholeskyDecomposition, out double)"/>
        /// with a minimum relative jitter of 1E-10 and a maximum of 1E-2.
        /// </remarks>
        public static bool TryCholesky(Matrix A, out CholeskyDecomposition? decomposition, out double jitter)
        {
            return TryCholesky(A, 1E-10, 1E-2, out decomposition, out jitter);
        }

        /// <summary>
        /// Attempts a Cholesky decomposition, retrying with an escalating diagonal jitter when the matrix is
        /// not numerically positive-definite.
        /// </summary>
        /// <param name="A">The symmetric input matrix A [0..n-1][0..n-1], which is not modified.</param>
        /// <param name="minimumRelativeJitter">The first jitter, as a fraction of the mean diagonal of A.</param>
        /// <param name="maximumRelativeJitter">The largest jitter allowed, as a fraction of the mean diagonal of A.</param>
        /// <param name="decomposition">
        /// The decomposition of <c>A + jitter * I</c> on success; otherwise <see langword="null"/>.
        /// </param>
        /// <param name="jitter">
        /// The absolute diagonal jitter that was added: zero when A factors as given, or NaN on failure.
        /// </param>
        /// <returns><see langword="true"/> when A or a jittered A was factored; otherwise <see langword="false"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when A is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when A is not square, or when the jitter bounds are not finite with
        /// <c>0 &lt; minimumRelativeJitter &lt;= maximumRelativeJitter</c>.
        /// </exception>
        /// <remarks>
        /// <para>
        /// A is first factored as given. If a pivot is rejected, the factorization is repeated on
        /// <c>A + δ_k * I</c> with <c>δ_k = s * minimumRelativeJitter * 10^k</c>, for k = 0, 1, ... while
        /// <c>δ_k</c> does not exceed <c>s * maximumRelativeJitter</c>. The scale <c>s</c> is the mean of the
        /// diagonal of A when that is a positive finite number, and one otherwise, so the jitter is invariant
        /// to the units of A. This is the diagonal "jitter" regularization commonly used for Gaussian process
        /// covariance matrices: the returned factor is the exact Cholesky factor of the jittered matrix, which
        /// is stored in <see cref="A"/>.
        /// </para>
        /// <para>
        /// Each attempt uses the same factorization and the same default pivot test,
        /// <see cref="DefaultRelativeTolerance(int)"/>, as <see cref="CholeskyDecomposition(Matrix)"/>, so a
        /// matrix the constructor accepts is returned with a zero jitter and a bit-identical factor. As in the
        /// constructor, only the upper triangle and diagonal of A are read. A matrix with a non-finite entry in
        /// that part is rejected without any attempt. A failed attempt costs a full O(n³) factorization.
        /// </para>
        /// </remarks>
        public static bool TryCholesky(Matrix A, double minimumRelativeJitter, double maximumRelativeJitter,
            out CholeskyDecomposition? decomposition, out double jitter)
        {
            if (A == null) throw new ArgumentNullException(nameof(A));
            if (A.NumberOfColumns != A.NumberOfRows)
                throw new ArgumentOutOfRangeException(nameof(A), "The matrix A must be square.");
            if (!Tools.IsFinite(minimumRelativeJitter) || minimumRelativeJitter <= 0d)
                throw new ArgumentOutOfRangeException(nameof(minimumRelativeJitter), "The minimum relative jitter must be a positive finite value.");
            if (!Tools.IsFinite(maximumRelativeJitter) || maximumRelativeJitter < minimumRelativeJitter)
                throw new ArgumentOutOfRangeException(nameof(maximumRelativeJitter), "The maximum relative jitter must be finite and no smaller than the minimum relative jitter.");

            decomposition = null;
            jitter = double.NaN;
            int dimension = A.NumberOfRows;
            double trace = 0d;
            for (int i = 0; i < dimension; i++)
            {
                for (int j = i; j < dimension; j++)
                {
                    if (!Tools.IsFinite(A[i, j])) return false;
                }
                trace += A[i, i];
            }

            double tolerance = DefaultRelativeTolerance(dimension);
            var candidate = new Matrix(A.Array);
            if (TryFactorizeFlat(candidate, tolerance, out double[] factor, out _, out _))
            {
                decomposition = new CholeskyDecomposition(candidate, factor, tolerance);
                jitter = 0d;
                return true;
            }

            double meanDiagonal = trace / dimension;
            double scale = meanDiagonal > 0d && Tools.IsFinite(meanDiagonal) ? meanDiagonal : 1d;
            // The small allowance keeps the maximum on the ladder when it is an exact decade above the minimum.
            int steps = (int)Math.Floor(Math.Log10(maximumRelativeJitter / minimumRelativeJitter) + 1E-9);
            for (int k = 0; k <= steps; k++)
            {
                double delta = scale * minimumRelativeJitter * Math.Pow(10d, k);
                if (!Tools.IsFinite(delta)) return false;
                candidate = new Matrix(A.Array);
                for (int i = 0; i < dimension; i++)
                    candidate[i, i] += delta;
                if (TryFactorizeFlat(candidate, tolerance, out factor, out _, out _))
                {
                    decomposition = new CholeskyDecomposition(candidate, factor, tolerance);
                    jitter = delta;
                    return true;
                }
            }
            return false;
        }

        private readonly int n; // Number of rows in A

        /// <summary>
        /// The factor as a flat row-major array with a zero strict upper triangle, used by the solve methods.
        /// </summary>
        private readonly double[] _factor;

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
        /// <remarks>
        /// A copy of the lower triangular factor. <see cref="Solve(Vector)"/>, <see cref="Forward(Vector)"/>,
        /// <see cref="Backward(Vector)"/> and <see cref="InverseA"/> use the factor computed by the decomposition, so
        /// changing entries of this matrix does not change their results.
        /// </remarks>
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
        /// <exception cref="ArgumentOutOfRangeException">Thrown when b does not have n entries.</exception>
        /// <remarks>
        /// Forward substitution with L followed by back substitution with L^T, by the blocked kernels from order
        /// <see cref="BlockedDimension"/> on and by scalar loops below it.
        /// </remarks>
        public Vector Solve(Vector b)
        {
            if (b.Length != n)
            {
                throw new ArgumentOutOfRangeException(nameof(b), "The vector b must have the same number of rows as the matrix A.");
            }
            var x = b.ToArray();
            if (n >= BlockedDimension)
            {
                CholeskyKernel.SolveLowerRows(_factor, 0, n, n, x, 0, 1, null);
                CholeskyKernel.SolveUpper(_factor, 0, n, n, x, 0);
                return new Vector(x);
            }
            var l = _factor;
            for (int i = 0; i < n; i++)
            {
                double sum = x[i];
                for (int k = i - 1; k >= 0; k--) sum -= l[i * n + k] * x[k];
                x[i] = sum / l[i * n + i];
            }
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = x[i];
                for (int k = i + 1; k < n; k++) sum -= l[k * n + i] * x[k];
                x[i] = sum / l[i * n + i];
            }
            return new Vector(x);
        }

        /// <summary>
        /// Solving the L^T * x = y equation with backward substitution
        /// </summary>
        /// <param name="y">The input vector y [0..n-1].</param>
        /// <returns>The solution x.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when y does not have n entries.</exception>
        /// <remarks>
        /// By the column-oriented kernel from order <see cref="BlockedDimension"/> on and by a scalar loop below it; both
        /// subtract the terms for j from n − 1 down to i + 1 in that order.
        /// </remarks>
        public Vector Backward(Vector y)
        {
            if (y.Length != n)
            {
                throw new ArgumentOutOfRangeException(nameof(y), "The vector y must have the same number of rows as the matrix A.");
            }
            var x = y.ToArray();
            if (n >= BlockedDimension)
            {
                CholeskyKernel.SolveUpper(_factor, 0, n, n, x, 0);
                return new Vector(x);
            }
            var l = _factor;
            for (int i = n - 1; i >= 0; i--)
            {
                double sum = x[i];
                for (int j = n - 1; j > i; j--) sum -= x[j] * l[j * n + i];
                x[i] = sum / l[i * n + i];
            }
            return new Vector(x);
        }

        /// <summary>
        /// Solving the L * y = b equation using Forward substitution
        /// </summary>
        /// <param name="b">The right-hand side vector b [0..n-1].</param>
        /// <returns>The solution y.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when b does not have n entries.</exception>
        /// <remarks>By the blocked kernel from order <see cref="BlockedDimension"/> on and by a scalar loop below it.</remarks>
        public Vector Forward(Vector b)
        {
            if (b.Length != n)
            {
                throw new ArgumentOutOfRangeException(nameof(b), "The vector b must have the same number of rows as the matrix A.");
            }
            var y = b.ToArray();
            if (n >= BlockedDimension)
            {
                CholeskyKernel.SolveLowerRows(_factor, 0, n, n, y, 0, 1, null);
                return new Vector(y);
            }
            var l = _factor;
            for (int i = 0; i < n; i++)
            {
                double sum = y[i];
                for (int j = 0; j < i; j++) sum -= l[i * n + j] * y[j];
                y[i] = sum / l[i * n + i];
            }
            return new Vector(y);
        }

        /// <returns>
        /// Matrix inverse A^-1 using the stored Cholesky decomposition.
        /// </returns>
        /// <remarks>
        /// From order <see cref="BlockedDimension"/> on, the inverse is formed as L^-T L^-1 by
        /// <see cref="CholeskyKernel.Invert"/>, whose residual max |A A^-1 - I| is within a factor of two of column-by-column
        /// substitution on random symmetric positive-definite test matrices. Below it, A X = I is
        /// solved column by column with forward and back substitution. The result is exactly symmetric.
        /// </remarks>
        public Matrix InverseA()
        {
            var Ainv = new Matrix(n);
            if (n >= BlockedDimension)
            {
                var work = (double[])_factor.Clone();
                var inverseFactor = new double[n * n];
                var diagonal = new double[n];
                CholeskyKernel.Invert(work, n, inverseFactor, diagonal);
                var target = Ainv.Array;
                for (int r = 0; r < n; r++)
                {
                    target[r, r] = diagonal[r];
                    int row = r * n;
                    for (int c = r + 1; c < n; c++)
                    {
                        target[r, c] = work[row + c];
                        target[c, r] = work[row + c];
                    }
                }
                return Ainv;
            }

            var l = _factor;
            int i, j, k;
            double sum;
            for (i = 0; i < n; i++)
            {
                for (j = 0; j <= i; j++)
                {
                    sum = i == j ? 1.0d : 0.0d;
                    for (k = i - 1; k >= j; k -= 1)
                        sum -= l[i * n + k] * Ainv[j, k];
                    Ainv[j, i] = sum / l[i * n + i];
                }
            }

            for (i = n - 1; i >= 0; i -= 1)
            {
                for (j = 0; j <= i; j++)
                {
                    sum = i < j ? 0.0d : Ainv[j, i];
                    for (k = i + 1; k < n; k++)
                        sum -= l[k * n + i] * Ainv[j, k];
                    Ainv[j, i] = sum / l[i * n + i];
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