using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Stats = Numerics.Data.Statistics.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Tests for <see cref="Stats.ColumnPercentiles(double[,], IList{double}, PlottingPositions.PlottingPostionType, out int[])"/>:
    /// each column equals <see cref="Stats.Percentile(IList{double}, IList{double}, PlottingPositions.PlottingPostionType, bool)"/>
    /// over its finite entries, bit for bit.
    /// </summary>
    /// <remarks>
    /// The known values are NumPy 2.5.4 <c>numpy.quantile</c>: <c>method="weibull"</c> (Hyndman and Fan type 6),
    /// <c>"hazen"</c> (type 5) and <c>"interpolated_inverted_cdf"</c> (type 4, Kaplan-Meier). NumPy forms the rank in
    /// another operation order, so the known values are compared to 1E-12; the column member is compared to the
    /// scalar member bit for bit.
    /// </remarks>
    [TestClass]
    public class Test_ColumnPlottingPositionPercentiles
    {
        /// <summary>The percentiles of the known-value cases.</summary>
        private static readonly double[] K = { 0.0, 0.1, 0.25, 0.5, 0.9, 1.0 };

        /// <summary>Every defined plotting-position formula.</summary>
        private static readonly PlottingPositions.PlottingPostionType[] Kinds =
        {
            PlottingPositions.PlottingPostionType.Weibull, PlottingPositions.PlottingPostionType.Median,
            PlottingPositions.PlottingPostionType.Blom, PlottingPositions.PlottingPostionType.Cunnane,
            PlottingPositions.PlottingPostionType.Gringorten, PlottingPositions.PlottingPostionType.Hazen,
            PlottingPositions.PlottingPostionType.KaplanMeier,
        };

        /// <summary>The two known-value columns: eight finite values, and six finite values among NaN and -infinity.</summary>
        private static double[,] KnownData() => new double[,]
        {
            { 3, 1.5 }, { 1, 2.5 }, { 4, 0.5 }, { 1, 7.0 }, { 5, double.NaN }, { 9, 3.0 }, { 2, double.NegativeInfinity }, { 6, 4.0 },
        };

        /// <summary>NumPy oracle values for Weibull, Hazen and Kaplan-Meier, with the non-finite entries excluded.</summary>
        [TestMethod]
        public void KnownValuesMatchNumPy()
        {
            var cases = new (PlottingPositions.PlottingPostionType Kind, double[] First, double[] Second)[]
            {
                (PlottingPositions.PlottingPostionType.Weibull, new[] { 1.0, 1.0, 1.25, 3.5, 9.0, 9.0 }, new[] { 0.5, 0.5, 1.25, 2.75, 7.0, 7.0 }),
                (PlottingPositions.PlottingPostionType.Hazen, new[] { 1.0, 1.0, 1.5, 3.5, 8.1, 9.0 }, new[] { 0.5, 0.6, 1.5, 2.75, 6.7, 7.0 }),
                (PlottingPositions.PlottingPostionType.KaplanMeier, new[] { 1.0, 1.0, 1.0, 3.0, 6.6, 9.0 }, new[] { 0.5, 0.5, 1.0, 2.5, 5.2, 7.0 }),
            };
            foreach (var c in cases)
            {
                var result = Stats.ColumnPercentiles(KnownData(), K, c.Kind, out int[] excluded);
                CollectionAssert.AreEqual(new[] { 0, 2 }, excluded);
                for (int i = 0; i < K.Length; i++)
                {
                    Assert.AreEqual(c.First[i], result[i, 0], 1E-12, $"{c.Kind} k={K[i]}");
                    Assert.AreEqual(c.Second[i], result[i, 1], 1E-12, $"{c.Kind} k={K[i]}");
                }
            }
        }

        /// <summary>Each column is bitwise the scalar member over its finite entries, for every formula, on a random ensemble.</summary>
        [TestMethod]
        public void ColumnsMatchPercentileBitwiseForEveryFormula()
        {
            var random = new Random(20261023);
            int rows = 123, columns = 37;
            var data = new double[rows, columns];
            for (int r = 0; r < rows; r++)
                for (int c = 0; c < columns; c++)
                {
                    double u = random.NextDouble();
                    data[r, c] = u < 0.03 ? double.NaN : u < 0.04 ? double.PositiveInfinity : 100 * Math.Exp(3 * random.NextDouble()) - 50;
                }
            data[0, 5] = -0.0;
            data[1, 5] = 0.0;
            data[2, 5] = -0.0;
            var k = new[] { 0.0, 1E-9, 0.01, 0.2, 0.5, 0.8, 0.96, 0.999, 1.0 };

            foreach (var kind in Kinds)
            {
                var result = Stats.ColumnPercentiles(data, k, kind, out int[] excluded);
                var plain = Stats.ColumnPercentiles(data, k, kind);
                for (int c = 0; c < columns; c++)
                {
                    var finite = new List<double>();
                    for (int r = 0; r < rows; r++)
                        if (!double.IsNaN(data[r, c]) && !double.IsInfinity(data[r, c])) finite.Add(data[r, c]);
                    var reference = Stats.Percentile(finite, k, kind);
                    Assert.AreEqual(rows - finite.Count, excluded[c]);
                    for (int i = 0; i < k.Length; i++)
                    {
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference[i]), BitConverter.DoubleToInt64Bits(result[i, c]), $"{kind} column {c} k={k[i]}");
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(result[i, c]), BitConverter.DoubleToInt64Bits(plain[i, c]));
                    }
                }
            }
        }

        /// <summary>Repeated calls over many columns (evaluated in parallel) are bitwise identical.</summary>
        [TestMethod]
        public void RepeatedCallsAreBitwiseIdentical()
        {
            var random = new Random(7);
            var data = new double[40, 500];
            for (int r = 0; r < 40; r++)
                for (int c = 0; c < 500; c++)
                    data[r, c] = 90 * random.NextDouble();
            var k = new[] { 0.5, 0.9, 0.96, 0.975 };
            var first = Stats.ColumnPercentiles(data, k, PlottingPositions.PlottingPostionType.Weibull);
            for (int repeat = 0; repeat < 5; repeat++)
            {
                var again = Stats.ColumnPercentiles(data, k, PlottingPositions.PlottingPostionType.Weibull);
                for (int i = 0; i < k.Length; i++)
                    for (int c = 0; c < 500; c++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(first[i, c]), BitConverter.DoubleToInt64Bits(again[i, c]));
            }
        }

        /// <summary>Columns without a finite entry give NaN; empty shapes give empty results; one row holds its value.</summary>
        [TestMethod]
        public void EdgeCases()
        {
            var nan = double.NaN;
            var result = Stats.ColumnPercentiles(new double[,] { { nan, 2.0 }, { double.NegativeInfinity, nan } }, K,
                PlottingPositions.PlottingPostionType.Weibull, out int[] excluded);
            CollectionAssert.AreEqual(new[] { 2, 1 }, excluded);
            for (int i = 0; i < K.Length; i++)
            {
                Assert.IsTrue(double.IsNaN(result[i, 0]));
                Assert.AreEqual(2.0, result[i, 1]);
            }

            var noRows = Stats.ColumnPercentiles(new double[0, 3], K, PlottingPositions.PlottingPostionType.Hazen, out int[] counts);
            Assert.AreEqual(K.Length, noRows.GetLength(0));
            Assert.AreEqual(3, noRows.GetLength(1));
            Assert.IsTrue(double.IsNaN(noRows[0, 0]));
            CollectionAssert.AreEqual(new[] { 0, 0, 0 }, counts);

            var noColumns = Stats.ColumnPercentiles(new double[4, 0], K, PlottingPositions.PlottingPostionType.Hazen);
            Assert.AreEqual(0, noColumns.GetLength(1));
            var noK = Stats.ColumnPercentiles(new double[,] { { 1.0 } }, new double[0], PlottingPositions.PlottingPostionType.Hazen);
            Assert.AreEqual(0, noK.GetLength(0));
        }

        /// <summary>Null arguments, percentiles outside [0, 1] or NaN, and an undefined formula are refused.</summary>
        [TestMethod]
        public void ArgumentErrors()
        {
            var weibull = PlottingPositions.PlottingPostionType.Weibull;
            var data = new double[,] { { 1.0 }, { 2.0 } };
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.ColumnPercentiles(null!, K, weibull));
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.ColumnPercentiles(data, null!, weibull));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.ColumnPercentiles(data, new[] { 1.5 }, weibull));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.ColumnPercentiles(data, new[] { -0.1 }, weibull));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.ColumnPercentiles(data, new[] { double.NaN }, weibull));
            var error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() =>
                Stats.ColumnPercentiles(data, K, (PlottingPositions.PlottingPostionType)99));
            Assert.AreEqual("plottingPositionType", error.ParamName);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() =>
                Stats.ColumnPercentiles(new double[2, 0], K, (PlottingPositions.PlottingPostionType)99));
        }
    }
}
