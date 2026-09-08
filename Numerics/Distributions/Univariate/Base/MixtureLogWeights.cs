using System;

namespace Numerics.Distributions
{
    /// <summary>Normalizes component weights for a single logarithmic mixture observation.</summary>
    internal static class MixtureLogWeights
    {
        /// <summary>Returns the row log probability and corresponding responsibilities.</summary>
        /// <remarks>Nonfinite or impossible rows are returned as nonfinite log probabilities so the
        /// caller can retain its observation-specific error message and aggregate likelihood convention.</remarks>
        internal static double Normalize(double[] logDensities, double[] weights, double[] responsibilities)
        {
            double maximum = double.NegativeInfinity;
            for (int i = 0; i < weights.Length; i++)
                if (weights[i] > 0) maximum = Math.Max(maximum, logDensities[i]);
            if (!DistributionNumerics.IsFinite(maximum)) return maximum;
            double weightedMaximum = double.NegativeInfinity;
            for (int i = 0; i < weights.Length; i++)
            {
                responsibilities[i] = weights[i] > 0 ? (logDensities[i] - maximum) + Math.Log(weights[i]) : double.NegativeInfinity;
                weightedMaximum = Math.Max(weightedMaximum, responsibilities[i]);
            }
            double sum = 0;
            for (int i = 0; i < weights.Length; i++) sum += Math.Exp(responsibilities[i] - weightedMaximum);
            double centeredLogRow = weightedMaximum + Math.Log(sum);
            for (int i = 0; i < weights.Length; i++) responsibilities[i] = Math.Exp(responsibilities[i] - centeredLogRow);
            return maximum + centeredLogRow;
        }
    }
}
