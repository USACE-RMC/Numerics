/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
        public ParameterStatistics SummaryStatistics { get; private set; }

        /// <summary>
        /// The kernel density results.
        /// </summary>
        [JsonInclude]
        public double[,] KernelDensity { get; private set; }

        /// <summary>
        /// The histogram results.
        /// </summary>
        [JsonInclude]
        public Histogram Histogram { get; private set; }

        /// <summary>
        /// The autocorrelation function for each parameter. This is averaged across each chain.
        /// </summary>
        public double[,] Autocorrelation { get; set; } = null!;

    }
}
