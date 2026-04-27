using System;

namespace Numerics.Data
{
    /// <summary>
    /// Enumeration of math function types.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum MathFunctionType
    {
        /// <summary>
        /// Compute the sum over each time block.
        /// </summary>
        Add,

        /// <summary>
        /// Compute the difference over each time block.
        /// </summary>
        Subtract,

        /// <summary>
        /// Compute the product over each time block.
        /// </summary>
        Multiply,

        /// <summary>
        /// Compute the quotient over each time block.
        /// </summary>
        Divide,

        /// <summary>
        /// Compute the logarithmic transform over each time block.
        /// </summary>
        Logarithm,

        /// <summary>
        /// Raise each time block to some power.
        /// </summary>
        Exponentiate,

        /// <summary>
        /// Compute the inverse over each time block.
        /// </summary>
        Inverse,

        /// <summary>
        /// Replace missing values.
        /// </summary>
        Replace,

        /// <summary>
        /// Interpolate missing values.
        /// </summary>
        Interpolate
    }
}
