using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{
    internal static partial class DistributionNumerics
    {
        /// <summary>Returns a finite magnitude used to evaluate the same initialization estimator in unit coordinates.</summary>
        internal static double InitializationScale(IList<double> sample)
        {
            double scale = 0;
            for (int i = 0; i < sample.Count; i++) scale = Math.Max(scale, Math.Abs(sample[i]));
            if (!(scale > 0) || !IsFinite(scale))
                throw new ArgumentOutOfRangeException(nameof(sample), "A finite nonzero sample magnitude is required for initialization.");
            return scale;
        }

        /// <summary>Constructs ordered finite positive bounds that contain a representable positive initial parameter.</summary>
        internal static void PositiveParameterBounds(double initial, out double lower, out double upper)
        {
            if (!(initial > 0) || !IsFinite(initial))
                throw new ArgumentOutOfRangeException(nameof(initial), "The initial parameter must be finite and positive.");
            lower = Math.Max(double.Epsilon, Math.Min(Tools.DoubleMachineEpsilon, initial / 10));
            double decade = Math.Pow(10, Math.Ceiling(Math.Log10(initial) + 1));
            upper = IsFinite(decade) ? Math.Max(initial, decade) : double.MaxValue;
            if (!(lower < upper))
                throw new ArgumentOutOfRangeException(nameof(initial), "Finite ordered parameter bounds cannot be represented.");
        }

        /// <summary>Constructs signed, scale-aware finite location bounds and a feasible initial location.</summary>
        /// <remarks>For a lower-endpoint family the upper bound is the actual sample minimum.</remarks>
        internal static void LocationParameterBounds(ref double initial, double scale, double dataMinimum,
            double dataMaximum, bool upperAtMinimum, out double lower, out double upper)
        {
            if (!IsFinite(initial) || !(scale > 0) || !IsFinite(scale))
                throw new ArgumentOutOfRangeException(nameof(initial), "Initialization requires a finite location and positive scale.");
            double magnitude = Math.Max(Math.Abs(initial), Math.Max(scale, Math.Max(Math.Abs(dataMinimum), Math.Abs(dataMaximum))));
            double radius = Math.Pow(10, Math.Ceiling(Math.Log10(magnitude) + 1));
            if (!IsFinite(radius)) radius = double.MaxValue;
            lower = -radius;
            upper = upperAtMinimum ? dataMinimum : radius;
            if (!(lower < upper))
                throw new ArgumentOutOfRangeException(nameof(initial), "Finite ordered location bounds cannot be represented for this sample.");
            if (initial < lower || initial > upper) initial = lower / 2 + upper / 2;
        }
    }
}
