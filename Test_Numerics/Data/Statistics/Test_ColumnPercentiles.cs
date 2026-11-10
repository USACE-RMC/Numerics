using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Stats = Numerics.Data.Statistics.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Tests for <see cref="Numerics.Data.Statistics.Statistics.ColumnPercentiles(double[,], IList{double}, out int[])"/>:
    /// each column equals <see cref="Numerics.Data.Statistics.Statistics.Percentile(IList{double}, IList{double}, bool)"/>
    /// over its finite entries, bit for bit.
    /// </summary>
    [TestClass]
    public class Test_ColumnPercentiles
    {
        /// <summary>The percentiles evaluated in every case, including both ends.</summary>
        private static readonly double[] K = { 0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0 };

        /// <summary>Known Type 7 values (R quantile(type = 7)) with non-finite entries excluded.</summary>
        [TestMethod]
        public void KnownValuesExcludeNonFiniteEntries()
        {
            // Column 0: {3, 1, 4, 1, 5, 9, 2, 6}; column 1: the same with NaN and infinities mixed in;
            // column 2: one finite entry; column 3: none.
            double nan = double.NaN, inf = double.PositiveInfinity;
            var data = new double[,]
            {
                { 3, nan, nan, nan },
                { 1, 3, nan, inf },
                { 4, 1, 7, nan },
                { 1, -inf, nan, nan },
                { 5, 4, nan, -inf },
                { 9, 1, nan, nan },
                { 2, 5, nan, nan },
                { 6, 9, nan, nan },
                { nan, 2, nan, nan },
                { inf, 6, nan, nan },
            };
            var result = Stats.ColumnPercentiles(data, new[] { 0.0, 0.25, 0.5, 0.9, 1.0 }, out int[] excluded);

            double[] expected = { 1.0, 1.75, 3.5, 6.9, 9.0 };
            for (int i = 0; i < expected.Length; i++)
            {
                Assert.AreEqual(expected[i], result[i, 0], 1E-12);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(result[i, 0]), BitConverter.DoubleToInt64Bits(result[i, 1]));
                Assert.AreEqual(7.0, result[i, 2]);
                Assert.IsTrue(double.IsNaN(result[i, 3]));
            }
            CollectionAssert.AreEqual(new[] { 2, 2, 9, 10 }, excluded);
        }

        /// <summary>Each column is bitwise the scalar member over its finite entries, on a random ensemble.</summary>
        [TestMethod]
        public void ColumnsMatchPercentileBitwise()
        {
            var random = new Random(20261023);
            int rows = 257, columns = 41;
            var data = new double[rows, columns];
            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < columns; c++)
                {
                    double u = random.NextDouble();
                    data[r, c] = u < 0.03 ? double.NaN : u < 0.04 ? double.NegativeInfinity : 100 * Math.Exp(3 * random.NextDouble()) - 50;
                }
            }
            data[0, 5] = -0.0;
            data[1, 5] = 0.0;

            var result = Stats.ColumnPercentiles(data, K, out int[] excluded);
            var plain = Stats.ColumnPercentiles(data, K);

            for (int c = 0; c < columns; c++)
            {
                var finite = new List<double>();
                for (int r = 0; r < rows; r++)
                    if (!double.IsNaN(data[r, c]) && !double.IsInfinity(data[r, c])) finite.Add(data[r, c]);
                var reference = Stats.Percentile(finite, K);
                Assert.AreEqual(rows - finite.Count, excluded[c]);
                for (int i = 0; i < K.Length; i++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference[i]), BitConverter.DoubleToInt64Bits(result[i, c]), $"k = {K[i]}, column {c}.");
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference[i]), BitConverter.DoubleToInt64Bits(plain[i, c]));
                }
            }
        }

        /// <summary>
        /// A single row gives that row at every percentile; a fully finite column equals the percentile over all
        /// its entries; and a large matrix (many columns, evaluated in parallel) is bitwise reproducible.
        /// </summary>
        [TestMethod]
        public void SingleRowFullyFiniteAndReproducible()
        {
            var single = Stats.ColumnPercentiles(new double[,] { { 4.0, -2.5, double.NaN } }, K, out int[] singleExcluded);
            for (int i = 0; i < K.Length; i++)
            {
                Assert.AreEqual(4.0, single[i, 0]);
                Assert.AreEqual(-2.5, single[i, 1]);
                Assert.IsTrue(double.IsNaN(single[i, 2]));
            }
            CollectionAssert.AreEqual(new[] { 0, 0, 1 }, singleExcluded);

            var random = new Random(7);
            int rows = 500, columns = 300;
            var data = new double[rows, columns];
            for (int r = 0; r < rows; r++)
                for (int c = 0; c < columns; c++)
                    data[r, c] = random.NextDouble() * 1000 - 500;
            var first = Stats.ColumnPercentiles(data, K, out int[] firstExcluded);
            for (int repeat = 0; repeat < 3; repeat++)
            {
                var again = Stats.ColumnPercentiles(data, K);
                for (int i = 0; i < K.Length; i++)
                    for (int c = 0; c < columns; c++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(first[i, c]), BitConverter.DoubleToInt64Bits(again[i, c]));
            }
            foreach (int c in new[] { 0, 151, 299 })
            {
                var column = new double[rows];
                for (int r = 0; r < rows; r++) column[r] = data[r, c];
                var reference = Stats.Percentile(column, K);
                for (int i = 0; i < K.Length; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference[i]), BitConverter.DoubleToInt64Bits(first[i, c]));
            }
            Assert.AreEqual(0, firstExcluded[42]);
        }

        /// <summary>Empty shapes give empty results; null inputs and percentiles outside [0, 1] are refused.</summary>
        [TestMethod]
        public void EdgeCasesAndArgumentErrors()
        {
            var noRows = Stats.ColumnPercentiles(new double[0, 3], K, out int[] excluded);
            Assert.AreEqual(K.Length, noRows.GetLength(0));
            Assert.AreEqual(3, noRows.GetLength(1));
            Assert.IsTrue(double.IsNaN(noRows[0, 0]));
            CollectionAssert.AreEqual(new[] { 0, 0, 0 }, excluded);
            Assert.AreEqual(0, Stats.ColumnPercentiles(new double[4, 0], K).GetLength(1));
            Assert.AreEqual(0, Stats.ColumnPercentiles(new double[4, 2], new double[0]).GetLength(0));

            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.ColumnPercentiles(null!, K));
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.ColumnPercentiles(new double[2, 2], null!));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.ColumnPercentiles(new double[2, 2], new[] { 0.5, 1.5 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.ColumnPercentiles(new double[2, 2], new[] { -0.1 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.ColumnPercentiles(new double[0, 0], new[] { double.NaN }));
        }
    }
}
