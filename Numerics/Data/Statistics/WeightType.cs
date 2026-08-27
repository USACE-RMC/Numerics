using System;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Enumeration of statistical weight interpretations for the weighted sample statistics.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public enum WeightType
    {
        /// <summary>
        /// Frequency (repeat-count) weights: a weight states how many observations the entry
        /// represents, and the weight total plays the role of the sample size in the bias
        /// corrections. Integer frequency weights reproduce the unweighted statistics of the
        /// correspondingly replicated sample exactly.
        /// </summary>
        Frequency,

        /// <summary>
        /// Reliability (precision or importance) weights: only relative weights carry meaning, and
        /// the bias corrections use the effective sample size W^2 / sum(w^2) in place of the count,
        /// with the variance denominator W - sum(w^2)/W.
        /// </summary>
        Reliability
    }
}
