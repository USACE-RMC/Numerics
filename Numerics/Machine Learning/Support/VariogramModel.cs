namespace Numerics.MachineLearning
{
    /// <summary>
    /// Isotropic variogram models for <see cref="Kriging"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Each model has semivariance <c>γ(h) = c₀ + c (1 - ρ(h))</c> for h &gt; 0 and γ(0) = 0, with nugget c₀,
    /// partial sill c and range parameter a. The range parameter follows the convention of the R package gstat
    /// (<c>vgm</c>). Isaaks and Srivastava (1989) write the exponential and Gaussian models with the practical
    /// range, at which the correlation has fallen to about 0.05; it equals 3a and √3 a respectively.
    /// </para>
    /// </remarks>
    public enum VariogramModel
    {
        /// <summary>
        /// Spherical, <c>ρ(h) = 1 - 1.5 (h/a) + 0.5 (h/a)³</c> for h &lt; a and zero beyond. Valid in up to three dimensions.
        /// </summary>
        Spherical,

        /// <summary>
        /// Exponential, <c>ρ(h) = exp(-h/a)</c>.
        /// </summary>
        Exponential,

        /// <summary>
        /// Gaussian, <c>ρ(h) = exp(-(h/a)²)</c>.
        /// </summary>
        Gaussian
    }

    /// <summary>
    /// Kriging variants, distinguished by their treatment of the mean.
    /// </summary>
    public enum KrigingType
    {
        /// <summary>
        /// Simple kriging: the mean is a known constant.
        /// </summary>
        Simple,

        /// <summary>
        /// Ordinary kriging: the mean is an unknown constant.
        /// </summary>
        Ordinary,

        /// <summary>
        /// Universal kriging: the mean is an unknown linear function of the coordinates.
        /// </summary>
        Universal
    }
}
