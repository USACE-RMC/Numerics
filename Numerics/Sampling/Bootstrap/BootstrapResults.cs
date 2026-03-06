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

using System;

namespace Numerics.Sampling
{

    /// <summary>
    /// Enumeration of bootstrap confidence interval methods.
    /// </summary>
    public enum BootstrapCIMethod
    {
        /// <summary>
        /// Percentile method. Direct percentile extraction from sorted bootstrap statistics.
        /// </summary>
        Percentile,

        /// <summary>
        /// Bias-corrected (BC) percentile method. Adjusts for median bias using the standard normal transformation.
        /// </summary>
        BiasCorrected,

        /// <summary>
        /// Bias-corrected and accelerated (BCa) method. Adjusts for both bias and skewness using jackknife acceleration.
        /// </summary>
        BCa,

        /// <summary>
        /// Normal (standard) method. Uses a configurable transform for transformation invariance with normal approximation.
        /// </summary>
        Normal,

        /// <summary>
        /// Bootstrap-t (studentized) method. Uses nested bootstrap to estimate standard errors and studentized pivotal statistic.
        /// </summary>
        BootstrapT
    }

    /// <summary>
    /// Stores bootstrap confidence interval results for a single statistic or parameter.
    /// </summary>
    [Serializable]
    public class BootstrapStatisticResult
    {
        /// <summary>
        /// The population (original) estimate for this statistic.
        /// </summary>
        public double PopulationEstimate { get; set; }

        /// <summary>
        /// The lower confidence interval bound.
        /// </summary>
        public double LowerCI { get; set; }

        /// <summary>
        /// The upper confidence interval bound.
        /// </summary>
        public double UpperCI { get; set; }

        /// <summary>
        /// The number of valid (non-NaN) bootstrap replicates used.
        /// </summary>
        public int ValidCount { get; set; }

        /// <summary>
        /// The total number of bootstrap replicates attempted.
        /// </summary>
        public int TotalCount { get; set; }

        /// <summary>
        /// The bootstrap standard error (standard deviation of valid replicates).
        /// </summary>
        public double StandardError { get; set; }

        /// <summary>
        /// The bootstrap mean of the valid replicates.
        /// </summary>
        public double Mean { get; set; }
    }

    /// <summary>
    /// Stores complete bootstrap analysis results including confidence intervals for statistics and parameters.
    /// </summary>
    [Serializable]
    public class BootstrapResults
    {
        /// <summary>
        /// The confidence interval method used.
        /// </summary>
        public BootstrapCIMethod Method { get; set; }

        /// <summary>
        /// The alpha level used (e.g., 0.1 for 90% CI).
        /// </summary>
        public double Alpha { get; set; }

        /// <summary>
        /// Results for each statistic (indexed by statistic index from StatisticFunction output).
        /// </summary>
        public BootstrapStatisticResult[] StatisticResults { get; set; } = Array.Empty<BootstrapStatisticResult>();

        /// <summary>
        /// Results for each model parameter (indexed by parameter index from ParameterSet.Values). Uses percentile CIs.
        /// </summary>
        public BootstrapStatisticResult[] ParameterResults { get; set; } = Array.Empty<BootstrapStatisticResult>();

        /// <summary>
        /// The total number of replicates that failed after all retries.
        /// </summary>
        public int FailedReplicates { get; set; }
    }
}
