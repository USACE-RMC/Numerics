using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for pooled independent-exclusive probabilities, including agreement with the
    /// allocating overload, output-row reuse, and early-convergence reporting.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_ProbabilityPooledExclusive
    {
        /// <summary>The three-event full enumeration: sizes 1, 2, 3 in combination order.</summary>
        private static readonly int[,] ThreeEventIndicators = new int[,]
        {
            { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 },
            { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 },
            { 1, 1, 1 },
        };

        /// <summary>The binomial combination counts for three events: C(3,1), C(3,2), C(3,3).</summary>
        private static readonly int[] ThreeEventCombinations = { 3, 3, 1 };

        /// <summary>
        /// Verifies agreement with the allocating overload for a complete enumeration.
        /// </summary>
        [TestMethod]
        public void Test_Pooled_MatchesAllocatingOverload()
        {
            double[] probabilities = { 0.3d, 0.2d, 0.1d };
            Probability.IndependentExclusive(probabilities, ThreeEventCombinations, ThreeEventIndicators,
                out List<double> expectedProbabilities, out List<int[]> expectedIndicators);

            var pooledProbabilities = new List<double>();
            var pooledIndicators = new List<int[]>();
            bool truncated = Probability.IndependentExclusive(probabilities, ThreeEventCombinations, ThreeEventIndicators,
                pooledProbabilities, pooledIndicators);

            Assert.IsFalse(truncated, "A full three-event enumeration never truncates.");
            Assert.HasCount(expectedProbabilities.Count, pooledProbabilities);
            Assert.HasCount(expectedIndicators.Count, pooledIndicators);
            for (int i = 0; i < expectedProbabilities.Count; i++)
            {
                Assert.AreEqual(expectedProbabilities[i], pooledProbabilities[i], 0d, $"Row {i} probability mismatch.");
                CollectionAssert.AreEqual(expectedIndicators[i], pooledIndicators[i], $"Row {i} indicator mismatch.");
            }
        }

        /// <summary>
        /// Test the steady-state reuse contract: a second call with the same shape refills the
        /// same row arrays in place (no per-call output allocations), and a smaller subsequent
        /// call trims the stale tail.
        /// </summary>
        [TestMethod]
        public void Test_Pooled_ReusesRowArraysAcrossCalls()
        {
            double[] first = { 0.3d, 0.2d, 0.1d };
            double[] second = { 0.05d, 0.4d, 0.25d };
            var pooledProbabilities = new List<double>();
            var pooledIndicators = new List<int[]>();

            Probability.IndependentExclusive(first, ThreeEventCombinations, ThreeEventIndicators, pooledProbabilities, pooledIndicators);
            var firstCallRow0 = pooledIndicators[0];
            var firstCallRow6 = pooledIndicators[6];

            Probability.IndependentExclusive(second, ThreeEventCombinations, ThreeEventIndicators, pooledProbabilities, pooledIndicators);
            Assert.AreSame(firstCallRow0, pooledIndicators[0], "Row arrays must be refilled in place across same-shape calls.");
            Assert.AreSame(firstCallRow6, pooledIndicators[6], "Every pooled row must be reused.");
            Assert.HasCount(7, pooledIndicators);

            // Compare with a fresh allocating call on the second inputs.
            Probability.IndependentExclusive(second, ThreeEventCombinations, ThreeEventIndicators,
                out List<double> expectedProbabilities, out _);
            for (int i = 0; i < expectedProbabilities.Count; i++)
            {
                Assert.AreEqual(expectedProbabilities[i], pooledProbabilities[i], 0d);
            }
        }

        /// <summary>
        /// Verifies that early convergence is reported and appends one closing pseudo-row.
        /// </summary>
        [TestMethod]
        public void Test_Pooled_TruncationIsReported()
        {
            // Ten rare events: the union converges within the first few inclusion-exclusion
            // terms at the default 1e-4 tolerances, so the expansion exits early.
            int n = 10;
            var probabilities = new double[n];
            for (int i = 0; i < n; i++) probabilities[i] = 0.01d;
            var (combinations, indicators) = BuildFullEnumeration(n);

            var pooledProbabilities = new List<double>();
            var pooledIndicators = new List<int[]>();
            bool truncated = Probability.IndependentExclusive(probabilities, combinations, indicators, pooledProbabilities, pooledIndicators);

            Assert.IsTrue(truncated, "Rare-event expansions converge early and must report the truncation.");
            Assert.IsLessThan(indicators.GetLength(0), pooledProbabilities.Count, "The deepest combinations were not enumerated.");
            Assert.HasCount(pooledProbabilities.Count, pooledIndicators);

            // The allocating overload produces the same truncated outputs.
            Probability.IndependentExclusive(probabilities, combinations, indicators,
                out List<double> expectedProbabilities, out _);
            Assert.HasCount(expectedProbabilities.Count, pooledProbabilities);
            for (int i = 0; i < expectedProbabilities.Count; i++)
            {
                Assert.AreEqual(expectedProbabilities[i], pooledProbabilities[i], 0d);
            }
        }

        /// <summary>
        /// Test that structurally inconsistent combination metadata and non-finite tolerances
        /// are rejected: the per-size counts must sum to the indicator row count, match the
        /// indicator layout, and the tolerance must be a finite number.
        /// </summary>
        [TestMethod]
        public void Test_Pooled_RejectsMalformedMetadata()
        {
            double[] probabilities = { 0.2d, 0.3d };
            var output = new List<double>();
            var indicators = new List<int[]>();
            int[,] rows = { { 1, 0 }, { 0, 1 }, { 1, 1 } };

            Assert.Throws<ArgumentException>(() => Probability.IndependentExclusive(
                probabilities, new[] { 2 }, rows, output, indicators));
            Assert.Throws<ArgumentException>(() => Probability.IndependentExclusive(
                probabilities, new[] { 1, 2 }, rows, output, indicators));
            Assert.Throws<ArgumentOutOfRangeException>(() => Probability.IndependentExclusive(
                probabilities, new[] { 2, 1 }, rows, output, indicators, double.NaN));
        }

        /// <summary>
        /// Builds the full combination enumeration (sizes 1..n in combination order) for n
        /// events.
        /// </summary>
        /// <param name="n">The event count.</param>
        /// <returns>The per-size combination counts and the indicator matrix.</returns>
        private static (int[] Combinations, int[,] Indicators) BuildFullEnumeration(int n)
        {
            var rows = new List<int[]>();
            var counts = new List<int>();
            for (int size = 1; size <= n; size++)
            {
                int before = rows.Count;
                foreach (var combo in EnumerateCombinations(size, n))
                {
                    var row = new int[n];
                    for (int j = 0; j < combo.Length; j++) row[combo[j]] = 1;
                    rows.Add(row);
                }
                counts.Add(rows.Count - before);
            }
            var indicators = new int[rows.Count, n];
            for (int i = 0; i < rows.Count; i++)
            {
                for (int j = 0; j < n; j++) indicators[i, j] = rows[i][j];
            }
            return (counts.ToArray(), indicators);
        }

        /// <summary>
        /// Enumerates all index combinations of a given size from n items, lexicographically.
        /// </summary>
        /// <param name="size">The combination size.</param>
        /// <param name="n">The item count.</param>
        /// <returns>The index combinations.</returns>
        private static IEnumerable<int[]> EnumerateCombinations(int size, int n)
        {
            var indexes = new int[size];
            for (int i = 0; i < size; i++) indexes[i] = i;
            while (true)
            {
                var copy = new int[size];
                Array.Copy(indexes, copy, size);
                yield return copy;

                int position = size - 1;
                while (position >= 0 && indexes[position] == n - size + position) position--;
                if (position < 0) yield break;
                indexes[position]++;
                for (int i = position + 1; i < size; i++) indexes[i] = indexes[i - 1] + 1;
            }
        }
    }
}
