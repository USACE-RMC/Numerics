using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{
    /// <summary>Provides scale-aware construction of finite initialization and parameter bounds for univariate distributions.</summary>
    internal static partial class DistributionNumerics
    {
        /// <summary>Returns a finite magnitude used to evaluate the same initialization estimator in unit coordinates.</summary>
        /// <param name="sample">The validated finite observations whose largest absolute magnitude defines the scale.</param>
        /// <returns>The largest absolute observation.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The sample has no finite positive magnitude.</exception>
        internal static double InitializationScale(IList<double> sample)
        {
            double scale = 0;
            for (int i = 0; i < sample.Count; i++) scale = Math.Max(scale, Math.Abs(sample[i]));
            if (!(scale > 0) || !Tools.IsFinite(scale))
                throw new ArgumentOutOfRangeException(nameof(sample), "A finite nonzero sample magnitude is required for initialization.");
            return scale;
        }

        /// <summary>Constructs ordered finite positive bounds that contain a representable positive initial parameter.</summary>
        /// <param name="initial">The finite positive initial parameter that the bounds must contain.</param>
        /// <param name="lower">The resulting finite positive lower bound.</param>
        /// <param name="upper">The resulting finite upper bound.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="initial"/> is not finite and positive, or ordered finite bounds cannot be represented.</exception>
        internal static void PositiveParameterBounds(double initial, out double lower, out double upper)
        {
            if (!(initial > 0) || !Tools.IsFinite(initial))
                throw new ArgumentOutOfRangeException(nameof(initial), "The initial parameter must be finite and positive.");
            lower = Math.Max(double.Epsilon, Math.Min(Tools.DoubleMachineEpsilon, initial / 10));
            double decade = Math.Pow(10, Math.Ceiling(Math.Log10(initial) + 1));
            upper = Tools.IsFinite(decade) ? Math.Max(initial, decade) : double.MaxValue;
            if (!(lower < upper))
                throw new ArgumentOutOfRangeException(nameof(initial), "Finite ordered parameter bounds cannot be represented.");
        }

        /// <summary>Constructs signed, scale-aware finite location bounds and a feasible initial location.</summary>
        /// <param name="initial">On input, the proposed finite location; on output, that value or a feasible midpoint when it was outside the constructed bounds.</param>
        /// <param name="scale">A finite positive characteristic sample scale.</param>
        /// <param name="dataMinimum">The minimum observation in physical coordinates.</param>
        /// <param name="dataMaximum">The maximum observation in physical coordinates.</param>
        /// <param name="upperAtMinimum"><see langword="true"/> to use <paramref name="dataMinimum"/> as the upper location bound; otherwise, <see langword="false"/>.</param>
        /// <param name="lower">The resulting finite lower location bound.</param>
        /// <param name="upper">The resulting finite upper location bound.</param>
        /// <exception cref="ArgumentOutOfRangeException">The initial location or scale is invalid, or ordered finite location bounds cannot be represented.</exception>
        /// <remarks>For a lower-endpoint family the upper bound is the actual sample minimum.</remarks>
        internal static void LocationParameterBounds(ref double initial, double scale, double dataMinimum,
            double dataMaximum, bool upperAtMinimum, out double lower, out double upper)
        {
            if (!Tools.IsFinite(initial) || !(scale > 0) || !Tools.IsFinite(scale))
                throw new ArgumentOutOfRangeException(nameof(initial), "Initialization requires a finite location and positive scale.");
            double magnitude = Math.Max(Math.Abs(initial), Math.Max(scale, Math.Max(Math.Abs(dataMinimum), Math.Abs(dataMaximum))));
            double radius = Math.Pow(10, Math.Ceiling(Math.Log10(magnitude) + 1));
            if (!Tools.IsFinite(radius)) radius = double.MaxValue;
            lower = -radius;
            upper = upperAtMinimum ? dataMinimum : radius;
            if (!(lower < upper))
                throw new ArgumentOutOfRangeException(nameof(initial), "Finite ordered location bounds cannot be represented for this sample.");
            if (initial < lower || initial > upper) initial = lower / 2 + upper / 2;
        }
    }
}
