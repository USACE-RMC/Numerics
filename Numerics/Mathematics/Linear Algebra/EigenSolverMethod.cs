using System;

namespace Numerics.Mathematics.LinearAlgebra
{

    /// <summary>
    /// The algorithm used by <see cref="EigenValueDecomposition"/> to decompose a real symmetric matrix.
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
    /// Both methods return the eigenvalues and orthonormal eigenvectors of A = V·diag(λ)·Vᵀ. The cyclic Jacobi method
    /// leaves the eigenvalues in diagonal order; Householder tridiagonalization followed by the implicit QL algorithm is
    /// several times faster for matrices beyond a few dozen rows and returns the eigenvalues in ascending order.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> Martin, R. S., Reinsch, C., and Wilkinson, J. H. (1968). Householder's tridiagonalization of a symmetric
    /// matrix. Numerische Mathematik, 11(3), 181-195. </item>
    /// <item> Bowdler, H., Martin, R. S., Reinsch, C., and Wilkinson, J. H. (1968). The QR and QL algorithms for
    /// symmetric matrices. Numerische Mathematik, 11(4), 293-306. </item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public enum EigenSolverMethod
    {
        /// <summary>
        /// The cyclic Jacobi rotation method. The eigenvalues are left in diagonal order and are not sorted.
        /// </summary>
        Jacobi,
        /// <summary>
        /// Householder reduction to tridiagonal form (tred2) followed by implicit QL iterations (tql2). The eigenvalues
        /// are in ascending order, with the eigenvectors in the matching columns.
        /// </summary>
        QL
    }
}
