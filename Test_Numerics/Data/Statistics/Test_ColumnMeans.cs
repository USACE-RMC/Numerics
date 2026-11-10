using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Stats = Numerics.Data.Statistics.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Tests for <see cref="Stats.ColumnMeans(double[,], out int[])"/>: each column equals
    /// <see cref="Stats.Mean(IList{double})"/> over its finite entries, bit for bit.
    /// </summary>
    /// <remarks>
    /// NumPy 2.5.4 <c>numpy.nanmean</c> over the finite entries is the independent oracle for the known values; NumPy
    /// sums pairwise rather than sequentially, so it is compared to 1E-14 relative error.
    /// </remarks>
    [TestClass]
    public class Test_ColumnMeans
    {
        /// <summary>Known values with non-finite entries excluded, by NumPy.</summary>
        [TestMethod]
        public void KnownValuesExcludeNonFiniteEntries()
        {
            double nan = double.NaN, inf = double.PositiveInfinity;
            var data = new double[,]
            {
                { 3, nan, nan, nan },
                { 1, 3, nan, inf },
                { 4, 1, 7, nan },
                { 1, -inf, nan, nan },
                { 5, 4, nan, -inf },
            };
            var means = Stats.ColumnMeans(data, out int[] excluded);
            Assert.AreEqual(2.8, means[0], 1E-14);
            Assert.AreEqual(8.0 / 3.0, means[1], 1E-14);
            Assert.AreEqual(7.0, means[2]);
            Assert.IsTrue(double.IsNaN(means[3]));
            CollectionAssert.AreEqual(new[] { 0, 2, 4, 5 }, excluded);
        }

        /// <summary>Each column is bitwise the sequential mean of its finite entries, on a random ensemble.</summary>
        [TestMethod]
        public void ColumnsMatchMeanBitwise()
        {
            var random = new Random(41);
            int rows = 1000, columns = 77;
            var data = new double[rows, columns];
            for (int r = 0; r < rows; r++)
                for (int c = 0; c < columns; c++)
                {
                    double u = random.NextDouble();
                    data[r, c] = u < 0.02 ? double.NaN : u < 0.03 ? double.NegativeInfinity : 1E3 * Math.Exp(4 * random.NextDouble()) - 500;
                }
            var means = Stats.ColumnMeans(data, out int[] excluded);
            var plain = Stats.ColumnMeans(data);
            for (int c = 0; c < columns; c++)
            {
                var finite = new List<double>();
                for (int r = 0; r < rows; r++)
                    if (!double.IsNaN(data[r, c]) && !double.IsInfinity(data[r, c])) finite.Add(data[r, c]);
                Assert.AreEqual(rows - finite.Count, excluded[c]);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Mean(finite)), BitConverter.DoubleToInt64Bits(means[c]), $"column {c}");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(means[c]), BitConverter.DoubleToInt64Bits(plain[c]));
            }
            for (int repeat = 0; repeat < 3; repeat++)
            {
                var again = Stats.ColumnMeans(data);
                for (int c = 0; c < columns; c++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(means[c]), BitConverter.DoubleToInt64Bits(again[c]));
            }
        }

        /// <summary>Empty shapes, and a null matrix.</summary>
        [TestMethod]
        public void EdgeCasesAndArgumentErrors()
        {
            var noRows = Stats.ColumnMeans(new double[0, 3], out int[] counts);
            Assert.HasCount(3, noRows);
            Assert.IsTrue(double.IsNaN(noRows[0]));
            CollectionAssert.AreEqual(new[] { 0, 0, 0 }, counts);
            Assert.IsEmpty(Stats.ColumnMeans(new double[4, 0]));
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.ColumnMeans(null!));
        }
    }
}
