using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// Exact 0/1 knapsack by dynamic programming over integer weights, for budget-constrained
    /// project selection and as an exact benchmark for heuristic portfolio optimizers.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The solver maximizes Σ v_i x_i subject to Σ w_i x_i ≤ W with x_i in {0, 1}, by the recursion
    /// f_i(c) = max(f_{i-1}(c), f_{i-1}(c - w_i) + v_i) over capacities c = 0..W. Time is O(n·W) and memory
    /// is O(n·W) bytes for the decision table. The multi-constraint overload (for example one budget per period)
    /// runs over the product state space of all capacities, O(n·m·Π(W_k + 1)), so it is meant for small instances.
    /// An item is selected only when it strictly increases the value, so items with non-positive value are
    /// never selected, and among optimal selections the later item is left out on a tie.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Bellman, R. (1957). Dynamic Programming. Princeton University Press.
    /// </description></item>
    /// </list>
    /// </remarks>
    public static class Knapsack
    {
        /// <summary>The largest decision table, in entries, the solver will allocate.</summary>
        private const long MAX_TABLE_ENTRIES = 1L << 30;

        /// <summary>
        /// Solves the single-constraint 0/1 knapsack problem exactly.
        /// </summary>
        /// <param name="values">The value of each item; must be finite.</param>
        /// <param name="weights">The non-negative integer weight (cost) of each item.</param>
        /// <param name="capacity">The non-negative integer capacity (budget).</param>
        /// <param name="selectedItems">The indices of the selected items, in increasing order.</param>
        /// <returns>The optimal total value, equal to the sum of the selected values in increasing index order.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="values"/> or <paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">The lengths differ, a value is not finite, a weight is negative, or the decision table is too large.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="capacity"/> is negative.</exception>
        public static double Solve(IList<double> values, IList<int> weights, int capacity, out int[] selectedItems)
        {
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            int n = ValidateValues(values, weights.Count);
            if (capacity < 0) throw new ArgumentOutOfRangeException(nameof(capacity), "The capacity must be non-negative.");
            for (int i = 0; i < n; i++)
                if (weights[i] < 0) throw new ArgumentException($"The weight of item {i} is negative.", nameof(weights));

            long tableWidth = (long)capacity + 1;
            if (Math.Max(n, 1) * tableWidth > MAX_TABLE_ENTRIES) throw new ArgumentException("The decision table n·(W + 1) is too large for the exact solver.", nameof(capacity));
            int width = (int)tableWidth;

            var best = new double[width];
            var keep = new bool[n * width];
            for (int i = 0; i < n; i++)
            {
                int w = weights[i];
                double v = values[i];
                if (w > capacity || !(v > 0)) continue;
                int row = i * width;
                for (int c = capacity; c >= w; c--)
                {
                    double candidate = best[c - w] + v;
                    if (candidate > best[c])
                    {
                        best[c] = candidate;
                        keep[row + c] = true;
                    }
                }
            }

            var selected = new List<int>();
            int remaining = capacity;
            for (int i = n - 1; i >= 0; i--)
            {
                if (keep[i * width + remaining])
                {
                    selected.Add(i);
                    remaining -= weights[i];
                }
            }
            selected.Reverse();
            selectedItems = selected.ToArray();
            return best[capacity];
        }

        /// <summary>
        /// Solves the multi-constraint (multidimensional) 0/1 knapsack problem exactly, for example one budget per period.
        /// </summary>
        /// <param name="values">The value of each item; must be finite.</param>
        /// <param name="weights">The non-negative integer weight of item i against constraint k, as weights[i, k].</param>
        /// <param name="capacities">The non-negative integer capacity of each constraint.</param>
        /// <param name="selectedItems">The indices of the selected items, in increasing order.</param>
        /// <returns>The optimal total value, equal to the sum of the selected values in increasing index order.</returns>
        /// <exception cref="ArgumentNullException">An argument is null.</exception>
        /// <exception cref="ArgumentException">The dimensions disagree, there is no constraint, a value is not finite, a weight or capacity is negative, or the state space is too large.</exception>
        public static double Solve(IList<double> values, int[,] weights, IList<int> capacities, out int[] selectedItems)
        {
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            if (capacities == null) throw new ArgumentNullException(nameof(capacities));
            int n = ValidateValues(values, weights.GetLength(0));
            int m = capacities.Count;
            if (m < 1) throw new ArgumentException("At least one constraint is required.", nameof(capacities));
            if (weights.GetLength(1) != m) throw new ArgumentException("The weights must have one column per constraint.", nameof(weights));

            var radix = new int[m];
            var stride = new int[m];
            long size = 1;
            for (int k = 0; k < m; k++)
            {
                if (capacities[k] < 0) throw new ArgumentException($"The capacity of constraint {k} is negative.", nameof(capacities));
                long digitCount = (long)capacities[k] + 1;
                if (digitCount * (size * Math.Max(n, 1)) > MAX_TABLE_ENTRIES) throw new ArgumentException("The state space n·Π(W_k + 1) is too large for the exact solver.", nameof(capacities));
                radix[k] = (int)digitCount;
                stride[k] = (int)size;
                size *= digitCount;
            }
            for (int i = 0; i < n; i++)
                for (int k = 0; k < m; k++)
                    if (weights[i, k] < 0) throw new ArgumentException($"The weight of item {i} against constraint {k} is negative.", nameof(weights));

            int P = (int)size;
            var best = new double[P];
            var keep = new bool[n * P];
            var digits = new int[m];
            for (int i = 0; i < n; i++)
            {
                double v = values[i];
                if (!(v > 0)) continue;
                int offset = 0;
                bool fits = true;
                for (int k = 0; k < m; k++)
                {
                    if (weights[i, k] > capacities[k]) { fits = false; break; }
                    offset += weights[i, k] * stride[k];
                }
                if (!fits) continue;

                // Visit states in decreasing flat index, tracking the mixed-radix digits of each index.
                for (int k = 0; k < m; k++) digits[k] = radix[k] - 1;
                int row = i * P;
                for (int idx = P - 1; idx >= 0; idx--)
                {
                    bool room = true;
                    for (int k = 0; k < m; k++)
                        if (digits[k] < weights[i, k]) { room = false; break; }
                    if (room)
                    {
                        double candidate = best[idx - offset] + v;
                        if (candidate > best[idx])
                        {
                            best[idx] = candidate;
                            keep[row + idx] = true;
                        }
                    }
                    for (int k = 0; k < m; k++)
                    {
                        if (digits[k] > 0) { digits[k]--; break; }
                        digits[k] = radix[k] - 1;
                    }
                }
            }

            var selected = new List<int>();
            int state = P - 1;
            for (int i = n - 1; i >= 0; i--)
            {
                if (keep[i * P + state])
                {
                    selected.Add(i);
                    for (int k = 0; k < m; k++) state -= weights[i, k] * stride[k];
                }
            }
            selected.Reverse();
            selectedItems = selected.ToArray();
            return best[P - 1];
        }

        /// <summary>
        /// Validates the item values against the item count.
        /// </summary>
        /// <param name="values">The item values.</param>
        /// <param name="count">The number of items implied by the weights.</param>
        /// <returns>The number of items.</returns>
        private static int ValidateValues(IList<double> values, int count)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Count != count) throw new ArgumentException("The values and weights must have one entry per item.", nameof(values));
            for (int i = 0; i < values.Count; i++)
                if (!Tools.IsFinite(values[i])) throw new ArgumentException($"The value of item {i} is not finite.", nameof(values));
            return values.Count;
        }
    }
}
