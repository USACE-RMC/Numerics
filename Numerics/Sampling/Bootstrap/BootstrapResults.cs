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
