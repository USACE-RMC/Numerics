using System;
using System.Collections.Generic;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Contains plotting position formulas.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The choice of quantiles from a theoretical distribution can depend upon context and purpose.
    /// Several different formulas have been used or proposed as affine symmetrical plotting positions.
    /// Such formulas have the form (k − a) / (n + 1 − 2a) for some value of a in the range from 0 to 1,
    /// which gives a range between k / (n + 1) and (k − 1) / (n - 1).
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot"/>
    /// </remarks>
    public class PlottingPositions
    {
        /// <summary>
        /// The general plotting position formula. This formula assumes the data is uncensored and complete.
        /// </summary>
        /// <param name="N">The sample size, at least one.</param>
        /// <param name="alpha">The alpha coefficient. Range [0,1].</param>
        /// <returns>An array of plotting positions of size N, (i − α) / (N + 1 − 2α) for i = 1..N.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="N"/> is less than one,
        /// <paramref name="alpha"/> is NaN or outside [0, 1], or N = 1 with α = 1, where the formula is 0 / 0.</exception>
        /// <remarks>
        /// The formula is defined for one and two values: a single value sits at the median position 1/2 for every
        /// α &lt; 1.
        /// </remarks>
        public static double[] Function(int N, double alpha)
        {
            if (N < 1)
            {
                throw new ArgumentOutOfRangeException(nameof(N), "The sample size N must be at least 1.");
            }

            if (double.IsNaN(alpha) || alpha < 0d || alpha > 1d)
            {
                throw new ArgumentOutOfRangeException(nameof(alpha), "The alpha coefficient must be between 0 and 1.");
            }

            if (N == 1 && alpha == 1d)
            {
                throw new ArgumentOutOfRangeException(nameof(alpha), "The alpha coefficient must be less than 1 for a single value.");
            }

            var PP = new double[N];
            for (int i = 1; i <= N; i++)
                PP[i - 1] = (i - alpha) / (N + 1 - 2d * alpha);
            return PP;
        }

        /// <summary>
        /// The general plotting position formula for weighted observations, for example importance-sampling
        /// weights.
        /// </summary>
        /// <param name="weights">The non-negative, finite weight of each observation, in the ascending order of the
        /// observations' values.</param>
        /// <param name="alpha">The alpha coefficient. Range [0,0.5], which holds every named formula from Weibull (0)
        /// to Hazen (0.5).</param>
        /// <returns>An array of non-decreasing plotting positions in [0,1], one per weight.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="weights"/> is empty or every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite, the total
        /// weight overflows, or <paramref name="alpha"/> is outside [0,0.5].</exception>
        /// <remarks>
        /// <para>
        /// With c_i the cumulative weight through observation i, W the total weight and w̄ = W / n the mean
        /// weight, the position of observation i is F_i = (c_i − α w_i) / (W + w̄ − 2α w̄), which is
        /// (c_i − α w_i) / (W + (1 − 2α) w̄). Unit weights give c_i = i, W = n and w̄ = 1, so every position is
        /// bitwise that of <see cref="Function(int, double)"/>, (i − α) / (n + 1 − 2α). A zero-weight observation
        /// takes the position of its predecessor.
        /// </para>
        /// <para>
        /// The cumulative and total weights are sequential sums in the order given, the total being the last
        /// cumulative weight.
        /// </para>
        /// <para>
        /// Successive positions differ by ((1 − α) w_i + α w_(i−1)) / denominator ≥ 0, and the first is ≥ 0. The
        /// last is at most one exactly when α w_n + (1 − 2α) w̄ ≥ 0, which holds for every weight vector exactly
        /// when α ≤ 1/2 (it equals one only at α = 1/2 with a zero last weight); above one half a small last weight
        /// gives a position above one, so alpha is limited to [0, 0.5] here. The unweighted formula accepts [0, 1]
        /// because its weights are equal.
        /// </para>
        /// </remarks>
        public static double[] Function(IList<double> weights, double alpha)
        {
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            if (weights.Count == 0) throw new ArgumentException("At least one weight is required.", nameof(weights));
            if (double.IsNaN(alpha) || alpha < 0d || alpha > 0.5d)
                throw new ArgumentOutOfRangeException(nameof(alpha), "The alpha coefficient must be between 0 and 0.5.");

            int n = weights.Count;
            var cumulative = new double[n];
            double run = 0d;
            for (int i = 0; i < n; i++)
            {
                double weight = weights[i];
                if (!Tools.IsFinite(weight) || weight < 0d)
                    throw new ArgumentOutOfRangeException(nameof(weights), "Each weight must be finite and non-negative.");
                run += weight;
                cumulative[i] = run;
            }
            double total = run;
            if (double.IsInfinity(total))
                throw new ArgumentOutOfRangeException(nameof(weights), "The total weight must be finite.");
            if (total <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));

            double mean = total / n;
            double denominator = total + mean - 2d * alpha * mean;
            var PP = new double[n];
            for (int i = 0; i < n; i++)
                PP[i] = (cumulative[i] - alpha * weights[i]) / denominator;
            return PP;
        }

        /// <summary>
        /// Returns the alpha coefficient of a named plotting position formula, (i - alpha) / (n + 1 - 2 alpha).
        /// </summary>
        /// <param name="plottingPostionType">The plotting position formula type.</param>
        /// <returns>0 for Weibull, 0.3175 for Median, 0.375 for Blom, 0.40 for Cunnane, 0.44 for Gringorten and 0.50 for Hazen.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the type is <see cref="PlottingPostionType.KaplanMeier"/>,
        /// which is the formula i / n rather than an alpha formula, or is not a defined value.</exception>
        /// <remarks>
        /// The single table of the alpha coefficients: the named formulas (<see cref="Weibull(int)"/> through
        /// <see cref="Hazen(int)"/>), the typed overloads <see cref="Function(int, PlottingPostionType)"/> and
        /// <see cref="Function(IList{double}, PlottingPostionType)"/>, and the plotting-position percentiles of
        /// <see cref="Statistics.Percentile(IList{double}, double, PlottingPostionType, bool)"/> all read it.
        /// </remarks>
        public static double Alpha(PlottingPostionType plottingPostionType)
        {
            switch (plottingPostionType)
            {
                case PlottingPostionType.Weibull: return 0d;
                case PlottingPostionType.Median: return 0.3175d;
                case PlottingPostionType.Blom: return 0.375d;
                case PlottingPostionType.Cunnane: return 0.4d;
                case PlottingPostionType.Gringorten: return 0.44d;
                case PlottingPostionType.Hazen: return 0.5d;
                case PlottingPostionType.KaplanMeier:
                    throw new ArgumentOutOfRangeException(nameof(plottingPostionType), "The Kaplan-Meier formula i / n has no alpha coefficient.");
                default:
                    throw new ArgumentOutOfRangeException(nameof(plottingPostionType), "Unknown plotting position type.");
            }
        }

        /// <summary>
        /// The named plotting position formula for weighted observations, for example importance-sampling weights.
        /// </summary>
        /// <param name="weights">The non-negative, finite weight of each observation, in the ascending order of the
        /// observations' values.</param>
        /// <param name="plottingPostionType">The plotting position formula type.</param>
        /// <returns>An array of non-decreasing plotting positions in [0,1], one per weight.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when <paramref name="weights"/> is empty or every weight is zero.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a weight is negative or not finite, the total
        /// weight overflows, or <paramref name="plottingPostionType"/> is not a defined value.</exception>
        /// <remarks>
        /// The alpha formulas are <see cref="Function(IList{double}, double)"/> at <see cref="Alpha(PlottingPostionType)"/>.
        /// Kaplan-Meier is F_i = c_i / W with c_i the sequential cumulative weight and W the total, so unit weights
        /// give bitwise i / n, the positions of <see cref="KaplanMeier(int)"/>. Every named alpha lies in [0, 0.5],
        /// where every weighted position lies in [0, 1].
        /// </remarks>
        public static double[] Function(IList<double> weights, PlottingPostionType plottingPostionType)
        {
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            if (plottingPostionType != PlottingPostionType.KaplanMeier)
                return Function(weights, Alpha(plottingPostionType));
            if (weights.Count == 0) throw new ArgumentException("At least one weight is required.", nameof(weights));
            int n = weights.Count;
            var cumulative = new double[n];
            double run = 0d;
            for (int i = 0; i < n; i++)
            {
                double weight = weights[i];
                if (!Tools.IsFinite(weight) || weight < 0d)
                    throw new ArgumentOutOfRangeException(nameof(weights), "Each weight must be finite and non-negative.");
                run += weight;
                cumulative[i] = run;
            }
            if (double.IsInfinity(run))
                throw new ArgumentOutOfRangeException(nameof(weights), "The total weight must be finite.");
            if (run <= 0d) throw new ArgumentException("The weights must not all be zero.", nameof(weights));
            var PP = new double[n];
            for (int i = 0; i < n; i++)
                PP[i] = cumulative[i] / run;
            return PP;
        }

        /// <summary>
        /// The general plotting position formula. This formula assumes the data is uncensored and complete.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <param name="plottingPostionType">The plotting position formula type.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[]? Function(int N, PlottingPositions.PlottingPostionType plottingPostionType)
        {
            if (plottingPostionType == PlottingPostionType.Blom)
            {
                return Blom(N);
            }
            else if (plottingPostionType == PlottingPostionType.Cunnane)
            {
                return Cunnane(N);
            }
            else if (plottingPostionType == PlottingPostionType.Gringorten)
            {
                return Gringorten(N);
            }
            else if (plottingPostionType == PlottingPostionType.Hazen)
            {
                return Hazen(N);
            }
            else if (plottingPostionType == PlottingPostionType.KaplanMeier)
            {
                return KaplanMeier(N);
            }
            else if (plottingPostionType == PlottingPostionType.Median)
            {
                return Median(N);
            }
            else if (plottingPostionType == PlottingPostionType.Weibull)
            {
                return Weibull(N);
            }
            return null;
        }

        /// <summary>
        /// The Weibull plotting position formula (alpha = 0.0). Recommended for uniform distribution.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Weibull(int N)
        {
            return Function(N, Alpha(PlottingPostionType.Weibull));
        }

        /// <summary>
        /// The Median plotting position formula (alpha = 0.3175). Provides median exceedance probabilities for all distributions.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Median(int N)
        {
            return Function(N, Alpha(PlottingPostionType.Median));
        }

        /// <summary>
        /// The Blom (1958) plotting position formula (alpha = 0.375). Recommended for Normal, Gamma, 2-parameter Log Normal,
        /// 3-parameter Log Normal, and Log Pearson Type III distributions.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Blom(int N)
        {
            return Function(N, Alpha(PlottingPostionType.Blom));
        }

        /// <summary>
        /// The Cunnane (1978) plotting position formula (alpha = 0.40). Recommended for GEV and Log-Gumbel distributions.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Cunnane(int N)
        {
            return Function(N, Alpha(PlottingPostionType.Cunnane));
        }

        /// <summary>
        /// The Gringorten (1963) plotting position formula (alpha = 0.44). Recommended for Exponential, Gumbel and Weibull distributions.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Gringorten(int N)
        {
            return Function(N, Alpha(PlottingPostionType.Gringorten));
        }

        /// <summary>
        /// The Hazen plotting position formula (alpha = 0.50). Recommended when the parameters of the parent distribution
        /// are unknown.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Hazen(int N)
        {
            return Function(N, Alpha(PlottingPostionType.Hazen));
        }

        /// <summary>
        /// The Kaplan-Meier plotting position formula for uncensored data (k/n).
        /// </summary>
        /// <param name="N">The sample size, at least one.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="N"/> is less than one.</exception>
        public static double[] KaplanMeier(int N)
        {
            if (N < 1)
            {
                throw new ArgumentOutOfRangeException(nameof(N), "The sample size N must be at least 1.");
            }
            var PP = new double[N];
            for (int i = 1; i <= N; i++)
                PP[i - 1] = i / (double)N;
            return PP;
        }

        /// <summary>
        /// Enumeration of plotting position types.
        /// </summary>
        public enum PlottingPostionType
        {
            /// <summary>
            /// The Weibull plotting position formula (alpha = 0.0).
            /// </summary>
            Weibull,

            /// <summary>
            /// The Median plotting position formula (alpha = 0.3175).
            /// </summary>
            Median,

            /// <summary>
            /// The Blom (1958) plotting position formula (alpha = 0.375).
            /// </summary>
            Blom,

            /// <summary>
            /// The Cunnane (1978) plotting position formula (alpha = 0.40).
            /// </summary>
            Cunnane,

            /// <summary>
            /// The Gringorten (1963) plotting position formula (alpha = 0.44).
            /// </summary>
            Gringorten,

            /// <summary>
            /// The Hazen plotting position formula (alpha = 0.50).
            /// </summary>
            Hazen,

            /// <summary>
            /// The Kaplan-Meier plotting position formula for uncensored data (k/n).
            /// </summary>
            KaplanMeier
        }
    }
}