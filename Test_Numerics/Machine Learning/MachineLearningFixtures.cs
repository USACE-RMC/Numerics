using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using System.Collections.Generic;

namespace MachineLearning
{
    /// <summary>
    /// Shared deterministic fixtures and exact-comparison helpers for the machine learning tests.
    /// </summary>
    internal static class MachineLearningFixtures
    {
        /// <summary>
        /// Builds a deterministic 30 by 2 fixture of three overlapping groups.
        /// </summary>
        public static double[,] ReferenceFixture()
        {
            var x = new double[30, 2];
            for (int i = 0; i < 30; i++)
            {
                x[i, 0] = 3 * Math.Sin(1.3 * i) + 4 * (i % 3);
                x[i, 1] = 2 * Math.Cos(0.7 * i) + 3 * (i % 3);
            }
            return x;
        }

        /// <summary>
        /// Builds a deterministic n by 4 fixture of five overlapping groups from integer arithmetic, so
        /// the values are identical on every runtime.
        /// </summary>
        public static double[,] LargeFixture(int n)
        {
            var x = new double[n, 4];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < 4; j++)
                    x[i, j] = 4 * ((i * 7 + j) % 5) + ((uint)(i * 2654435761u + (uint)j * 40503u) % 1000u) / 333.0 + 0.5 * j;
            return x;
        }

        /// <summary>
        /// Returns the FNV-1a hash of the bit patterns of a sequence of values.
        /// </summary>
        public static ulong BitHash(IEnumerable<double> values)
        {
            ulong hash = 14695981039346656037UL;
            foreach (var v in values)
            {
                hash ^= (ulong)BitConverter.DoubleToInt64Bits(v);
                hash *= 1099511628211UL;
            }
            return hash;
        }

        /// <summary>
        /// Asserts exact equality of a matrix against row-major expected values.
        /// </summary>
        public static void AssertExact(double[] expected, double[,] actual)
        {
            int c = actual.GetLength(1);
            Assert.HasCount(expected.Length, actual);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i / c, i % c], 0d);
        }
    }
}
