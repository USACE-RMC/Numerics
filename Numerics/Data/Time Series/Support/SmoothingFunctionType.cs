using System;

namespace Numerics.Data
{
    /// <summary>
    /// Enumeration of time series smoothing function options.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>

    [Serializable]
    public enum SmoothingFunctionType
    {
        /// <summary>
        /// Perform successive differences of the data.
        /// </summary>
        Difference,

        /// <summary>
        /// Perform a forward moving average of the data. 
        /// </summary>
        MovingAverage,

        /// <summary>
        /// Performs a forward moving sum of the data.
        /// </summary>
        MovingSum,

        /// <summary>
        /// Do nothing.
        /// </summary>
        None,
    }
}
