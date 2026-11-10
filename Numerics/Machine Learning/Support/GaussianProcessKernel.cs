namespace Numerics.MachineLearning
{
    /// <summary>
    /// Stationary covariance functions for <see cref="GaussianProcessRegression"/>.
    /// </summary>
    /// <remarks>
    /// Each covariance is <c>k(x, x') = σ_f² ρ(r)</c> with the scaled distance
    /// <c>r = sqrt(Σ_j (x_j - x'_j)² / ℓ_j²)</c>, where the length scales ℓ_j are shared or set per feature
    /// (automatic relevance determination). The forms are those of Rasmussen and Williams (2006), section 4.2,
    /// and of scikit-learn's <c>RBF</c> and <c>Matern</c> kernels.
    /// </remarks>
    public enum GaussianProcessKernel
    {
        /// <summary>
        /// Squared exponential (radial basis function), <c>ρ(r) = exp(-r²/2)</c>. Infinitely differentiable.
        /// </summary>
        SquaredExponential,

        /// <summary>
        /// Matérn with ν = 3/2, <c>ρ(r) = (1 + √3 r) exp(-√3 r)</c>. Once differentiable.
        /// </summary>
        Matern32,

        /// <summary>
        /// Matérn with ν = 5/2, <c>ρ(r) = (1 + √5 r + 5r²/3) exp(-√5 r)</c>. Twice differentiable.
        /// </summary>
        Matern52,

        /// <summary>
        /// Exponential (Matérn with ν = 1/2), <c>ρ(r) = exp(-r)</c>. Continuous but not differentiable.
        /// </summary>
        Exponential
    }

    /// <summary>
    /// Mean functions for <see cref="GaussianProcessRegression"/>.
    /// </summary>
    /// <remarks>
    /// A nonzero mean is a linear combination of basis functions whose coefficients are estimated by
    /// generalized least squares for the current covariance, as in universal kriging. The log marginal
    /// likelihood is evaluated at those coefficients, and predictive variances include their estimation
    /// uncertainty.
    /// </remarks>
    public enum GaussianProcessMean
    {
        /// <summary>
        /// Zero mean.
        /// </summary>
        Zero,

        /// <summary>
        /// Unknown constant mean, estimated by generalized least squares.
        /// </summary>
        Constant,

        /// <summary>
        /// Mean linear in the features, <c>β₀ + Σ_j β_j x_j</c>, estimated by generalized least squares.
        /// </summary>
        Linear
    }

    /// <summary>
    /// Covariance forms used internally, including those available only to <see cref="Kriging"/>.
    /// </summary>
    internal enum CovarianceForm
    {
        /// <summary>Squared exponential.</summary>
        SquaredExponential,

        /// <summary>Matérn 3/2.</summary>
        Matern32,

        /// <summary>Matérn 5/2.</summary>
        Matern52,

        /// <summary>Exponential.</summary>
        Exponential,

        /// <summary>Spherical, <c>ρ(r) = 1 - 1.5 r + 0.5 r³</c> for r &lt; 1 and zero beyond.</summary>
        Spherical
    }
}
