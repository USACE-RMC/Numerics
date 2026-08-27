using System;

namespace Numerics.Data
{
    /// <summary>
    /// Enumeration of the sides of an ordered data range on which out-of-range lookups may
    /// extrapolate rather than hold the boundary value.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// Sides are defined in value space and are independent of the sort orientation:
    /// <see cref="Below"/> refers to lookups beyond the minimum value of the lookup axis and
    /// <see cref="Above"/> to lookups beyond its maximum. Extrapolation extends the boundary
    /// segment linearly in the configured transform space, so a logarithmic axis extends
    /// log-linearly and a normal-probability axis extends linearly in standard-normal z - the
    /// latter therefore always back-transforms into (0, 1).
    /// </para>
    /// </remarks>
    [Flags]
    [Serializable]
    public enum ExtrapolationSides
    {
        /// <summary>
        /// No extrapolation: out-of-range lookups hold the boundary value on both sides.
        /// </summary>
        None = 0,

        /// <summary>
        /// Extrapolate below the minimum value of the lookup axis.
        /// </summary>
        Below = 1,

        /// <summary>
        /// Extrapolate above the maximum value of the lookup axis.
        /// </summary>
        Above = 2,

        /// <summary>
        /// Extrapolate on both sides of the range.
        /// </summary>
        Both = Below | Above
    }
}
