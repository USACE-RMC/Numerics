using System;
using System.Collections.Generic;
using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{
    internal static partial class DistributionNumerics
    {
        /// <summary>
        /// Locates the global maximum of a density on [lower, upper] from its values on an equally spaced grid, then
        /// polishes the best grid peaks with a bounded Brent search of the density itself.
        /// </summary>
        /// <param name="pdf">The density.</param>
        /// <param name="lower">The left end of the grid.</param>
        /// <param name="upper">The right end of the grid; greater than <paramref name="lower"/>.</param>
        /// <param name="gridDensity">The density, or any approximation of it that ranks the peaks, at the grid points
        /// lower + i·(upper − lower)/(n − 1), i = 0, …, n − 1; at least two points.</param>
        /// <param name="exactGrid">Whether <paramref name="gridDensity"/> holds the values of <paramref name="pdf"/>
        /// itself at the grid points, so that they are not evaluated again.</param>
        /// <returns>The location of the largest density found.</returns>
        /// <remarks>
        /// The grid's local maxima whose value is at least 0.9 times the largest, at most three of them in decreasing
        /// order, are each polished with <see cref="BrentSearch"/> over the two grid cells around them, and the polished
        /// point with the largest exact density wins; a grid point beats its polish only if its exact density is larger.
        /// The polish brackets a local maximum because a grid peak is no lower than its neighbours. A grid point counts as
        /// a peak when it is no lower than its left neighbour and higher than its right one, so a plateau resolves to its
        /// right end; peaks of equal value are taken from left to right. An end point is polished only over the one cell
        /// inside the grid, so a maximum beyond [lower, upper] is not found. The search is global up to the grid
        /// resolution: a peak narrower than a grid cell, or one whose grid value is below 0.9 of the best, can be missed.
        /// The three highest peaks are kept during the single scan of the grid, so no peak list is built or sorted.
        /// </remarks>
        internal static double GridMode(Func<double, double> pdf, double lower, double upper, IList<double> gridDensity, bool exactGrid = false)
        {
            int n = gridDensity.Count;
            double step = (upper - lower) / (n - 1);

            // Grid peaks: points no lower than both neighbours (an end point compares with its one neighbour). The three
            // highest are kept in decreasing value, equal values in increasing index; best is the highest peak value.
            int first = -1, second = -1, third = -1;
            double firstValue = 0d, secondValue = 0d, thirdValue = 0d;
            for (int i = 0; i < n; i++)
            {
                double d = gridDensity[i];
                if (!(d >= 0d)) continue;
                bool left = i == 0 || d >= gridDensity[i - 1];
                bool right = i == n - 1 || d > gridDensity[i + 1];
                if (!(left && right)) continue;
                // Strictly greater, so an equal value keeps the earlier index ahead.
                if (first < 0 || d > firstValue)
                {
                    third = second; thirdValue = secondValue;
                    second = first; secondValue = firstValue;
                    first = i; firstValue = d;
                }
                else if (second < 0 || d > secondValue)
                {
                    third = second; thirdValue = secondValue;
                    second = i; secondValue = d;
                }
                else if (third < 0 || d > thirdValue)
                {
                    third = i; thirdValue = d;
                }
            }
            if (first < 0) return lower + 0.5d * (upper - lower);
            double best = firstValue;

            double mode = double.NaN, modeDensity = double.NegativeInfinity;
            for (int k = 0; k < 3; k++)
            {
                int i = k == 0 ? first : k == 1 ? second : third;
                if (i < 0) break;
                if (k > 0 && gridDensity[i] < 0.9d * best) break;
                double x = lower + i * step;
                double a = Math.Max(lower, x - step), b = Math.Min(upper, x + step);
                var polish = new BrentSearch(pdf, a, b) { ReportFailure = false, ComputeHessian = false, RecordTraces = false };
                polish.Maximize();
                double candidate = polish.BestParameterSet.Values[0];
                // Maximize minimizes the negated density, so the best fitness is the negated density at the candidate.
                double candidateDensity = -polish.BestParameterSet.Fitness;
                double gridPointDensity = exactGrid ? gridDensity[i] : pdf(x);
                if (gridPointDensity > candidateDensity)
                {
                    candidate = x;
                    candidateDensity = gridPointDensity;
                }
                if (candidateDensity > modeDensity)
                {
                    mode = candidate;
                    modeDensity = candidateDensity;
                }
            }
            return mode;
        }
    }
}
