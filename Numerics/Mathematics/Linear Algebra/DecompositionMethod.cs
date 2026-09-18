using System;

namespace Numerics.Mathematics.LinearAlgebra
{

    /// <summary>
    /// The matrix decomposition method used to factorize a symmetric covariance matrix.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// A covariance matrix Σ is factorized to obtain a matrix A satisfying A·Aᵀ = Σ, which correlates
    /// independent standard normal variates, and to obtain the log-determinant and the quadratic form
    /// required by a Gaussian density. The two members below trade computational speed against the range
    /// of supported covariance structures: <see cref="Cholesky"/> is the faster factorization but requires
    /// a strictly positive-definite matrix, while <see cref="SingularValue"/> also handles the singular
    /// (perfectly or near-perfectly collinear) covariance matrices that arise, for example, in gridded
    /// data where adjacent pixels are perfectly correlated.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> "Numerical Recipes: The Art of Scientific Computing, Third Edition." Press et al., 2017. </item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public enum DecompositionMethod
    {
        /// <summary>
        /// Cholesky decomposition Σ = L·Lᵀ. This is the faster factorization — roughly twice as efficient
        /// as an LU decomposition and considerably cheaper than a singular value decomposition — and is the
        /// default. It requires a strictly positive-definite matrix and fails when the matrix is singular,
        /// so choose it when the covariance is known to be well conditioned.
        /// </summary>
        Cholesky,
        /// <summary>
        /// Singular value decomposition Σ = U·W·Vᵀ. This is the slower factorization but supports a wider
        /// range of covariance structures: it tolerates singular (collinear) matrices by working on the
        /// affine support of the distribution, using the rank, the pseudo-determinant and the pseudo-inverse
        /// in place of the dimension, the determinant and the inverse. Choose it when perfect or
        /// near-perfect correlation can occur, such as gridded data with adjacent, highly correlated cells.
        /// </summary>
        SingularValue
    }
}
