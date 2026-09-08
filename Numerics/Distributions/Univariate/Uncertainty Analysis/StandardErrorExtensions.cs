using System.Collections.Generic;

namespace Numerics.Distributions
{
    /// <summary>Stable, additive quantile-uncertainty operations for existing distributions.</summary>
    public static class StandardErrorExtensions
    {
        /// <summary>Returns the logarithm of the absolute quantile Jacobian determinant.</summary>
        /// <param name="distribution">A distribution whose gradients use its public parameter coordinates.</param>
        /// <param name="probabilities">One finite, strictly interior probability per parameter, in row order.</param>
        /// <returns>The log absolute determinant, or negative infinity for an exactly singular Jacobian.</returns>
        /// <exception cref="System.ArgumentNullException">An argument is null.</exception>
        /// <exception cref="System.ArgumentOutOfRangeException">The probability count or an individual probability is invalid.</exception>
        /// <exception cref="System.InvalidOperationException">A derivative is not finite.</exception>
        /// <remarks>Equilibrates rows and columns before pivoting and sums log pivots. This avoids forming
        /// a determinant that can overflow or underflow while its logarithm remains representable.
        /// No interface members are added and no artificial singularity pivots are used.</remarks>
        public static double LogAbsQuantileJacobian(this IStandardError distribution, IList<double> probabilities)
        {
            var matrix = DistributionNumerics.QuantileGradientMatrix(distribution, probabilities);
            return DistributionNumerics.LogAbsDeterminant(matrix, out _);
        }
    }
}
