using System;

namespace Numerics.Data
{
    /// <summary>
    /// Enumeration of time block function types.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum BlockFunctionType
    {
        /// <summary>
        /// Compute the minimum over each time block.
        /// </summary>
        Minimum,

        /// <summary>
        /// Compute the maximum over each time block.
        /// </summary>
        Maximum,

        /// <summary>
        /// Compute the average over each time block.
        /// </summary>
        Average,

        /// <summary>
        /// Compute the sum over each time block.
        /// </summary>
        Sum,
    }
}
