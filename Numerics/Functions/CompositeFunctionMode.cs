namespace Numerics.Functions
{
    /// <summary>
    /// The combination modes of a <see cref="CompositeFunction"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public enum CompositeFunctionMode
    {
        /// <summary>
        /// The pointwise weighted average of the child functions: F(x) = Σ wᵢ·fᵢ(x).
        /// </summary>
        WeightedAverage,

        /// <summary>
        /// A mixture: one confidence-level draw selects a child by cumulative weight and
        /// re-scales the remainder as the child's own confidence level, so a single uniform
        /// drives the whole composition (no internal random source).
        /// </summary>
        Mixture,
    }
}
