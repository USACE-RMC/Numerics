using Numerics.Data.Statistics;
using Numerics.Distributions;
using System;
using System.Text.Json.Serialization;

namespace Numerics.Sampling.MCMC
{

    /// <summary>
    /// A class for saving Bayesian MCMC results for each parameter.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class ParameterResults
    {

        /// <summary>
        /// Parameterless constructor for JSON deserialization.
        /// </summary>
        [JsonConstructor]
        public ParameterResults() { }

        /// <summary>
        /// Constructs new parameter results.
        /// </summary>
        /// <param name="values">List of posterior parameter values, aggregated together from each chain.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="sorted">Determines if the values have been sorted. Default = false.</param>
        public ParameterResults(double[] values, double alpha = 0.1, bool sorted = false)
        {
            // Sort the values
            if (sorted == false)
            {
                Array.Sort(values);
            }

            // Create Kernel Density Estimate
            var kde = new KernelDensity(values);
            KernelDensity = kde.CreatePDFGraph();

            // Set summary statistics
            SummaryStatistics = new ParameterStatistics();
            SummaryStatistics.N = values.Length;
            SummaryStatistics.Mean = kde.Mean;
            SummaryStatistics.StandardDeviation = kde.StandardDeviation;
            SummaryStatistics.Median = Statistics.Percentile(values, 0.5, true);
            SummaryStatistics.LowerCI = Statistics.Percentile(values, alpha / 2d, true);
            SummaryStatistics.UpperCI = Statistics.Percentile(values, 1d - alpha / 2d, true);

            // Create Histogram
            Histogram = new Histogram(values);

        }

        /// <summary>
        /// Parameter summary statistics.
        /// </summary>
        [JsonInclude]
        public ParameterStatistics SummaryStatistics { get; private set; } = null!;

        /// <summary>
        /// The kernel density results.
        /// </summary>
        [JsonInclude]
        public double[,] KernelDensity { get; private set; } = new double[0, 0];

        /// <summary>
        /// The histogram results.
        /// </summary>
        [JsonInclude]
        public Histogram Histogram { get; private set; } = null!;

        /// <summary>
        /// The autocorrelation function for each parameter. This is averaged across each chain.
        /// </summary>
        public double[,] Autocorrelation { get; set; } = new double[0, 0];

    }
}
