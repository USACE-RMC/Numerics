using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Sampling;

namespace Distributions.Univariate
{
    /// <summary>
    /// Tests for the stratified-range graph tables inherited from <see cref="UnivariateDistributionBase"/>:
    /// <see cref="UnivariateDistributionBase.CreatePDFGraph(List{StratificationBin}[])"/> and
    /// <see cref="UnivariateDistributionBase.CreateCDFGraph(List{StratificationBin}[])"/>.
    /// </summary>
    [TestClass]
    public class Test_UnivariateDistributionBaseGraphs
    {
        /// <summary>Checks that two tables have the same shape and the same bits in every cell.</summary>
        /// <param name="expected">The expected table.</param>
        /// <param name="actual">The actual table.</param>
        private static void AssertBitwiseEqual(double[,] expected, double[,] actual)
        {
            Assert.AreEqual(expected.GetLength(0), actual.GetLength(0));
            Assert.AreEqual(expected.GetLength(1), actual.GetLength(1));
            for (int i = 0; i < expected.GetLength(0); i++)
                for (int j = 0; j < expected.GetLength(1); j++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i, j]), BitConverter.DoubleToInt64Bits(actual[i, j]));
        }

        /// <summary>
        /// A single range gives the same table, bit for bit, as the single-list overload.
        /// </summary>
        [TestMethod]
        public void SingleRange_MatchesSingleListOverload()
        {
            var normal = new Normal(3d, 2d);
            var bins = Stratify.XValues(new StratificationOptions(-3d, 9d, 50));
            AssertBitwiseEqual(normal.CreatePDFGraph(bins), normal.CreatePDFGraph(new[] { bins }));
            AssertBitwiseEqual(normal.CreateCDFGraph(bins), normal.CreateCDFGraph(new[] { bins }));
        }

        /// <summary>
        /// Several ranges give one row per bin with the ranges concatenated in order and no gap rows.
        /// </summary>
        [TestMethod]
        public void MultipleRanges_ConcatenatesEveryBinWithoutGaps()
        {
            var gev = new GeneralizedExtremeValue(8d, 2d, 0.1d);
            var ranges = new[]
            {
                Stratify.XValues(new StratificationOptions(2d, 6d, 4)),
                Stratify.XValues(new StratificationOptions(6d, 12d, 7)),
                Stratify.XValues(new StratificationOptions(12d, 20d, 3)),
            };
            var pdf = gev.CreatePDFGraph(ranges);
            var cdf = gev.CreateCDFGraph(ranges);
            Assert.AreEqual(14, pdf.GetLength(0));
            Assert.AreEqual(14, cdf.GetLength(0));
            int row = 0;
            foreach (var range in ranges)
            {
                var pdfPart = gev.CreatePDFGraph(range);
                var cdfPart = gev.CreateCDFGraph(range);
                for (int j = 0; j < range.Count; j++, row++)
                {
                    Assert.AreEqual(range[j].Midpoint, pdf[row, 0]);
                    Assert.AreEqual(gev.PDF(range[j].Midpoint), pdf[row, 1]);
                    Assert.AreEqual(range[j].Midpoint, cdf[row, 0]);
                    Assert.AreEqual(gev.CDF(range[j].Midpoint), cdf[row, 1]);
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(pdfPart[j, 1]), BitConverter.DoubleToInt64Bits(pdf[row, 1]));
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(cdfPart[j, 1]), BitConverter.DoubleToInt64Bits(cdf[row, 1]));
                }
            }
            for (int i = 1; i < cdf.GetLength(0); i++)
            {
                Assert.IsGreaterThan(cdf[i - 1, 0], cdf[i, 0]);
                Assert.IsGreaterThan(cdf[i - 1, 1], cdf[i, 1]);
            }
        }

        /// <summary>
        /// Empty ranges contribute no rows, and an array of empty ranges gives an empty table.
        /// </summary>
        [TestMethod]
        public void EmptyRanges_ContributeNoRows()
        {
            var normal = new Normal();
            var bins = Stratify.XValues(new StratificationOptions(-1d, 1d, 5));
            var empty = new List<StratificationBin>();
            var pdf = normal.CreatePDFGraph(new[] { empty, bins, empty });
            var cdf = normal.CreateCDFGraph(new[] { empty, bins, empty });
            AssertBitwiseEqual(normal.CreatePDFGraph(bins), pdf);
            AssertBitwiseEqual(normal.CreateCDFGraph(bins), cdf);
            Assert.AreEqual(0, normal.CreatePDFGraph(new[] { empty, empty }).GetLength(0));
            Assert.AreEqual(0, normal.CreateCDFGraph(Array.Empty<List<StratificationBin>>()).GetLength(0));
        }
    }
}
