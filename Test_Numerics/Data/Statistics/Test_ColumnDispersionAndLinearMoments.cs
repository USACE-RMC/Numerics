using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Stats = Numerics.Data.Statistics.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Tests for <see cref="Stats.ColumnVariances(double[,], out int[])"/>,
    /// <see cref="Stats.ColumnStandardDeviations(double[,], out int[])"/> and
    /// <see cref="Stats.ColumnLinearMoments(double[,], out int[])"/>: each column equals the one-sample member over its
    /// finite entries, bit for bit.
    /// </summary>
    /// <remarks>
    /// The known values are exact rational evaluations (fractions.Fraction) of the n − 1 variance and of Hosking's
    /// unbiased probability weighted moments, correctly rounded; the standard deviations are 50-digit mpmath square
    /// roots of the exact variances.
    /// </remarks>
    [TestClass]
    public class Test_ColumnDispersionAndLinearMoments
    {
        private static readonly double[,] Data =
        {
            { 3.2, 101.5, 0.25 },
            { 1.7, 99.25, double.NaN },
            { 4.4, 100.75, 0.5 },
            { 2.9, 98.5, 0.125 },
            { 5.1, 102.0, double.PositiveInfinity },
            { 3.8, 100.0, 0.75 },
            { 2.2, 101.25, 0.375 },
            { 4.9, 99.75, 1.5 },
        };

        /// <summary>Known variances and standard deviations, with non-finite entries excluded.</summary>
        [TestMethod]
        public void VariancesMatchExactRationals()
        {
            var variances = Stats.ColumnVariances(Data, out var excluded);
            var deviations = Stats.ColumnStandardDeviations(Data, out var excludedSd);
            double[] exactVariance = { 1.5421428571428573, 1.4464285714285714, 0.24791666666666667 };
            double[] exactSd = { 1.2418304462135148, 1.2026755886059097, 0.4979123082096552 };
            for (int j = 0; j < 3; j++)
            {
                Assert.AreEqual(exactVariance[j], variances[j], 2e-16 * exactVariance[j], $"variance {j}");
                Assert.AreEqual(exactSd[j], deviations[j], 2e-16 * exactSd[j], $"sd {j}");
            }
            CollectionAssert.AreEqual(new[] { 0, 0, 2 }, excluded);
            CollectionAssert.AreEqual(excluded, excludedSd);
            CollectionAssert.AreEqual(variances, Stats.ColumnVariances(Data));
            CollectionAssert.AreEqual(deviations, Stats.ColumnStandardDeviations(Data));
        }

        /// <summary>Known L-moments, with non-finite entries excluded.</summary>
        [TestMethod]
        public void LinearMomentsMatchExactRationals()
        {
            var moments = Stats.ColumnLinearMoments(Data, out var excluded);
            double[,] exact =
            {
                { 3.525, 100.375, 0.5833333333333334 },
                { 0.7571428571428571, 0.7321428571428571, 0.2833333333333333 },
                { -0.05188679245283026, -0.06097560975609756, 0.4264705882352941 },
                { -0.04716981132075484, 0.012195121951219513, 0.3382352941176471 },
            };
            Assert.AreEqual(4, moments.GetLength(0));
            // λ2 = 2 b1 − b0 cancels when the mean is large relative to the L-scale: column 1 (mean 100.375, λ2 0.73)
            // measures 1.3E-14 relative error in LinearMoments itself, which this member reproduces bitwise.
            for (int r = 0; r < 4; r++)
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(exact[r, j], moments[r, j], 1e-13 * Math.Max(1d, Math.Abs(exact[r, j])), $"row {r} column {j}");
            CollectionAssert.AreEqual(new[] { 0, 0, 2 }, excluded);
        }

        /// <summary>
        /// Every column is bitwise the one-sample member over its gathered finite entries, including columns too short for
        /// the statistic, and the L-CV is row 1 over row 0, bitwise <see cref="Stats.LinearCoefficientOfVariation(IList{double})"/>.
        /// </summary>
        [TestMethod]
        public void ColumnsEqualOneSampleMembersBitwise()
        {
            var random = new Random(31);
            int rows = 60, columns = 40;
            var data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    data[i, j] = random.NextDouble() < 0.05 ? double.NaN : 1e6 + 50d * random.NextDouble();
            // Columns with no, one and three finite entries.
            for (int i = 0; i < rows; i++) { data[i, 0] = double.NaN; data[i, 1] = i == 7 ? 4d : double.NaN; data[i, 2] = i < 3 ? i : double.NegativeInfinity; }

            var variances = Stats.ColumnVariances(data);
            var deviations = Stats.ColumnStandardDeviations(data);
            var moments = Stats.ColumnLinearMoments(data);
            for (int j = 0; j < columns; j++)
            {
                var column = new List<double>();
                for (int i = 0; i < rows; i++) if (!double.IsNaN(data[i, j]) && !double.IsInfinity(data[i, j])) column.Add(data[i, j]);
                double variance = column.Count == 0 ? double.NaN : Stats.Variance(column);
                double sd = column.Count == 0 ? double.NaN : Stats.StandardDeviation(column);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(variance), BitConverter.DoubleToInt64Bits(variances[j]), $"variance {j}");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(sd), BitConverter.DoubleToInt64Bits(deviations[j]), $"sd {j}");
                var expected = column.Count == 0 ? new[] { double.NaN, double.NaN, double.NaN, double.NaN } : Stats.LinearMoments(column);
                for (int r = 0; r < 4; r++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[r]), BitConverter.DoubleToInt64Bits(moments[r, j]), $"moment {r} column {j}");
                if (column.Count >= 4)
                    Assert.AreEqual(Stats.LinearCoefficientOfVariation(column), moments[1, j] / moments[0, j], 0d, $"L-CV {j}");
            }
            Assert.IsTrue(double.IsNaN(variances[0]) && double.IsNaN(variances[1]) && !double.IsNaN(variances[2]));
            Assert.IsTrue(double.IsNaN(moments[0, 2]));
        }

        /// <summary>
        /// The column members give each column bitwise the one-sample member over a <see cref="List{T}"/> of its finite
        /// entries both below the work threshold, where the columns run sequentially (6 × 4), and above it, where they run
        /// in parallel (300 × 60); and the one-sample mean and variance give the same bits for an array, an array segment
        /// and a list of the same values.
        /// </summary>
        [TestMethod]
        public void ColumnsMatchOneSampleMembersOnSequentialAndParallelPaths()
        {
            var random = new Random(41);
            foreach (var (rows, columns) in new[] { (6, 4), (300, 60) })
            {
                var data = new double[rows, columns];
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < columns; j++)
                        data[i, j] = random.NextDouble() < 0.05 ? double.NaN : 1e3 * Math.Exp(random.NextDouble() * 3d);
                var means = Stats.ColumnMeans(data);
                var variances = Stats.ColumnVariances(data);
                var deviations = Stats.ColumnStandardDeviations(data);
                var moments = Stats.ColumnLinearMoments(data);
                for (int j = 0; j < columns; j++)
                {
                    var column = new List<double>();
                    for (int i = 0; i < rows; i++) if (!double.IsNaN(data[i, j])) column.Add(data[i, j]);
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Mean(column)), BitConverter.DoubleToInt64Bits(means[j]), $"{rows} x {columns} mean {j}");
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Variance(column)), BitConverter.DoubleToInt64Bits(variances[j]), $"{rows} x {columns} variance {j}");
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.StandardDeviation(column)), BitConverter.DoubleToInt64Bits(deviations[j]), $"{rows} x {columns} sd {j}");
                    var expected = column.Count < 4 ? new[] { double.NaN, double.NaN, double.NaN, double.NaN } : Stats.LinearMoments(column);
                    for (int r = 0; r < 4; r++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[r]), BitConverter.DoubleToInt64Bits(moments[r, j]), $"{rows} x {columns} moment {r} column {j}");
                }
            }

            var values = new double[257];
            for (int i = 0; i < values.Length; i++) values[i] = 1e8 + random.NextDouble();
            var padded = new double[values.Length + 5];
            Array.Copy(values, 0, padded, 3, values.Length);
            var segment = new ArraySegment<double>(padded, 3, values.Length);
            var list = new List<double>(values);
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Mean(list)), BitConverter.DoubleToInt64Bits(Stats.Mean(values)));
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Mean(list)), BitConverter.DoubleToInt64Bits(Stats.Mean(segment)));
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Variance(list)), BitConverter.DoubleToInt64Bits(Stats.Variance(values)));
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Variance(list)), BitConverter.DoubleToInt64Bits(Stats.Variance(segment)));
        }

        /// <summary>Null input is rejected; an empty matrix gives empty results.</summary>
        [TestMethod]
        public void InvalidAndEmptyInput()
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.ColumnVariances(null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.ColumnStandardDeviations(null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.ColumnLinearMoments(null!));
            Assert.IsEmpty(Stats.ColumnVariances(new double[5, 0]));
            Assert.AreEqual(0, Stats.ColumnLinearMoments(new double[5, 0]).GetLength(1));
        }
    }
}
