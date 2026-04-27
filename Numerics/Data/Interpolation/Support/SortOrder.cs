using System;

namespace Numerics.Data
{
    /// <summary>
    /// Enumeration defining the sort order of values.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum SortOrder
    {
        /// <summary>
        /// Ascending order.
        /// </summary>
        Ascending,

        /// <summary>
        /// Descending order.
        /// </summary>
        Descending,

        /// <summary>
        /// Not sorted. 
        /// </summary>
        None
    }
}
