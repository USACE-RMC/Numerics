using System;

namespace Numerics.Sampling.MCMC
{

    /// <summary>
    /// A class for storing parameter statistics.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class ParameterStatistics
    {

        /// <summary>
        /// Constructs an empty parameter summary statistics class.
        /// </summary>
        public ParameterStatistics() { }

        /// <summary>
        /// The Gelman-Rubin diagnostic.
        /// </summary>
        public double Rhat { get; set; } = double.NaN;

        /// <summary>
        /// The effective sample size.
        /// </summary>
        public double ESS { get; set; } = double.NaN;

        /// <summary>
        /// The total sample size.
        /// </summary>
        public int N { get; set; }

        /// <summary>
        /// The parameter mean.
        /// </summary>
        public double Mean { get; set; }

        /// <summary>
        /// The parameter median.
        /// </summary>
        public double Median { get; set; }

        /// <summary>
        /// The parameter standard deviation.
        /// </summary>
        public double StandardDeviation { get; set; }

        /// <summary>
        /// The lower confidence interval.
        /// </summary>
        public double LowerCI { get; set; }

        /// <summary>
        /// The upper confidence interval.
        /// </summary>
        public double UpperCI { get; set; }

    }
}
