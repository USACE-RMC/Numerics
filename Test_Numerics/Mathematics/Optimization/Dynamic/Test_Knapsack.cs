using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Tests the exact 0/1 knapsack solvers against brute-force enumeration and an independent numpy oracle.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_Knapsack
    {
        /// <summary>
        /// Brute-force best subset over m constraints; the value is summed in increasing index order.
        /// </summary>
        private static double BruteForce(double[] values, int[,] weights, int[] capacities, out List<int> best)
        {
            int n = values.Length, m = capacities.Length;
            double bestValue = 0;
            best = new List<int>();
            for (int mask = 0; mask < (1 << n); mask++)
            {
                double v = 0;
                bool fits = true;
                for (int k = 0; k < m && fits; k++)
                {
                    int w = 0;
                    for (int i = 0; i < n; i++) if ((mask >> i & 1) == 1) w += weights[i, k];
                    fits = w <= capacities[k];
                }
                if (!fits) continue;
                var items = new List<int>();
                for (int i = 0; i < n; i++) if ((mask >> i & 1) == 1) { v += values[i]; items.Add(i); }
                if (v > bestValue) { bestValue = v; best = items; }
            }
            return bestValue;
        }

        /// <summary>
        /// Random single-constraint instances: the value equals brute force exactly and the selection is the
        /// brute-force selection.
        /// </summary>
        [TestMethod]
        public void SingleConstraintMatchesBruteForce()
        {
            var rng = new MersenneTwister(97);
            for (int trial = 0; trial < 40; trial++)
            {
                int n = rng.Next(1, 13);
                var values = new double[n];
                var weights = new int[n];
                var w2 = new int[n, 1];
                for (int i = 0; i < n; i++)
                {
                    values[i] = rng.NextDouble() * 20 - 2;
                    weights[i] = rng.Next(0, 15);
                    w2[i, 0] = weights[i];
                }
                int capacity = rng.Next(0, 40);
                double expected = BruteForce(values, w2, new[] { capacity }, out var expectedItems);
                double actual = Knapsack.Solve(values, weights, capacity, out int[] items);
                Assert.AreEqual(expected, actual, $"trial {trial}");
                CollectionAssert.AreEqual(expectedItems, items, $"trial {trial}");
                double sum = 0;
                int weight = 0;
                foreach (int i in items) { sum += values[i]; weight += weights[i]; }
                Assert.AreEqual(sum, actual);
                Assert.IsLessThanOrEqualTo(capacity, weight);
            }
        }

        /// <summary>
        /// Random two- and three-constraint (multi-period budget) instances match brute force exactly, and the
        /// one-constraint overload agrees with the single-constraint solver.
        /// </summary>
        [TestMethod]
        public void MultiConstraintMatchesBruteForce()
        {
            var rng = new MersenneTwister(2024);
            for (int trial = 0; trial < 30; trial++)
            {
                int n = rng.Next(1, 11);
                int m = 1 + trial % 3;
                var values = new double[n];
                var weights = new int[n, m];
                var capacities = new int[m];
                for (int i = 0; i < n; i++)
                {
                    values[i] = rng.NextDouble() * 10;
                    for (int k = 0; k < m; k++) weights[i, k] = rng.Next(0, 8);
                }
                for (int k = 0; k < m; k++) capacities[k] = rng.Next(0, 20);
                double expected = BruteForce(values, weights, capacities, out var expectedItems);
                double actual = Knapsack.Solve(values, weights, capacities, out int[] items);
                Assert.AreEqual(expected, actual, $"trial {trial}");
                CollectionAssert.AreEqual(expectedItems, items, $"trial {trial}");
                if (m == 1)
                {
                    var single = new int[n];
                    for (int i = 0; i < n; i++) single[i] = weights[i, 0];
                    double singleValue = Knapsack.Solve(values, single, capacities[0], out int[] singleItems);
                    Assert.AreEqual(expected, singleValue);
                    CollectionAssert.AreEqual(items, singleItems);
                }
            }
        }

        /// <summary>
        /// The formula fixture of the numpy oracle (n = 30, W = 400) gives the oracle's value and selection.
        /// </summary>
        [TestMethod]
        public void MatchesNumpyOracle()
        {
            int n = 30;
            var values = new double[n];
            var weights = new int[n];
            for (int i = 0; i < n; i++)
            {
                weights[i] = 1 + (53 * i) % 197;
                values[i] = 1.0 + (37 * i) % 101 + 0.001 * i;
            }
            double value = Knapsack.Solve(values, weights, 400, out int[] items);
            Assert.AreEqual(682.1410000000001, value, 1E-12);
            CollectionAssert.AreEqual(new[] { 0, 4, 5, 8, 15, 16, 19, 23, 24, 27 }, items);
        }

        /// <summary>
        /// Capacity zero, items all too heavy, ties, zero-weight items, non-positive values and no items.
        /// </summary>
        [TestMethod]
        public void EdgeCases()
        {
            Assert.AreEqual(0.0, Knapsack.Solve(new[] { 5.0, 3.0 }, new[] { 1, 2 }, 0, out int[] items));
            Assert.IsEmpty(items);

            Assert.AreEqual(0.0, Knapsack.Solve(new[] { 5.0, 3.0 }, new[] { 11, 12 }, 10, out items));
            Assert.IsEmpty(items);

            // Two identical items and room for one: the earlier item is kept.
            Assert.AreEqual(4.0, Knapsack.Solve(new[] { 4.0, 4.0 }, new[] { 3, 3 }, 5, out items));
            CollectionAssert.AreEqual(new[] { 0 }, items);

            // A zero-weight positive item is always taken; zero and negative values never are.
            Assert.AreEqual(9.0, Knapsack.Solve(new[] { 2.0, 0.0, -1.0, 7.0 }, new[] { 0, 0, 0, 4 }, 4, out items));
            CollectionAssert.AreEqual(new[] { 0, 3 }, items);

            Assert.AreEqual(0.0, Knapsack.Solve(Array.Empty<double>(), Array.Empty<int>(), 7, out items));
            Assert.IsEmpty(items);

            Assert.AreEqual(6.0, Knapsack.Solve(new[] { 6.0, 5.0 }, new[,] { { 2, 0 }, { 0, 3 } }, new[] { 2, 2 }, out items));
            CollectionAssert.AreEqual(new[] { 0 }, items);
        }

        /// <summary>
        /// Invalid inputs are rejected with clear exceptions.
        /// </summary>
        [TestMethod]
        public void InvalidInputs()
        {
            Assert.Throws<ArgumentNullException>(() => Knapsack.Solve(null!, new[] { 1 }, 1, out _));
            Assert.Throws<ArgumentNullException>(() => Knapsack.Solve(new[] { 1.0 }, (int[])null!, 1, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { 1.0, 2.0 }, new[] { 1 }, 1, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { double.NaN }, new[] { 1 }, 1, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { double.PositiveInfinity }, new[] { 1 }, 1, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { 1.0 }, new[] { -1 }, 1, out _));
            Assert.Throws<ArgumentOutOfRangeException>(() => Knapsack.Solve(new[] { 1.0 }, new[] { 1 }, -1, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new double[3000], new int[3000], int.MaxValue - 1, out _));

            Assert.Throws<ArgumentNullException>(() => Knapsack.Solve(new[] { 1.0 }, (int[,])null!, new[] { 1 }, out _));
            Assert.Throws<ArgumentNullException>(() => Knapsack.Solve(new[] { 1.0 }, new int[1, 1], null!, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { 1.0 }, new int[1, 0], Array.Empty<int>(), out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { 1.0 }, new int[1, 2], new[] { 1 }, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { 1.0 }, new int[1, 1], new[] { -1 }, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { 1.0 }, new[,] { { -1 } }, new[] { 1 }, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { 1.0, 2.0 }, new int[1, 1], new[] { 1 }, out _));
            Assert.Throws<ArgumentException>(() => Knapsack.Solve(new[] { 1.0 }, new int[1, 3], new[] { 5000, 5000, 5000 }, out _));
        }

        /// <summary>
        /// A capacity of int.MaxValue, whose table width overflows an int, is rejected as too large rather than
        /// overflowing, and an empty item list cannot bypass the table-size limit.
        /// </summary>
        [TestMethod]
        public void LargestCapacityIsRejected()
        {
            Assert.ThrowsExactly<ArgumentException>(() => Knapsack.Solve(new[] { 1.0 }, new[] { 1 }, int.MaxValue, out _));
            Assert.ThrowsExactly<ArgumentException>(() => Knapsack.Solve(new[] { 1.0 }, new[,] { { 1 } }, new[] { int.MaxValue }, out _));
            Assert.ThrowsExactly<ArgumentException>(() => Knapsack.Solve(Array.Empty<double>(), Array.Empty<int>(), int.MaxValue, out _));
            Assert.ThrowsExactly<ArgumentException>(() => Knapsack.Solve(Array.Empty<double>(), new int[0, 1], new[] { int.MaxValue }, out _));

            Assert.AreEqual(0d, Knapsack.Solve(Array.Empty<double>(), Array.Empty<int>(), 100, out var none));
            Assert.IsEmpty(none);
            Assert.AreEqual(0d, Knapsack.Solve(Array.Empty<double>(), new int[0, 2], new[] { 10, 10 }, out none));
            Assert.IsEmpty(none);
        }
    }
}
