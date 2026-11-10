using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Stats = Numerics.Data.Statistics.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Regression tests for the sort inside <see cref="Stats.Percentile(IList{double}, double, bool)"/> and the
    /// plotting-position percentiles: every value must stay bitwise what the stable LINQ <c>OrderBy</c> sort the
    /// members used before produced, including the order of −0 and +0 and of NaN payloads.
    /// </summary>
    /// <remarks>
    /// The reference sorts with <c>data.OrderBy(x =&gt; x).ToArray()</c> and passes the sorted array with
    /// <c>dataIsSorted = true</c>, which is exactly the former unsorted path. NumPy's <c>numpy.percentile</c> (linear,
    /// Type 7) is the independent oracle for the known values.
    /// </remarks>
    [TestClass]
    public class Test_PercentileSortOrder
    {
        /// <summary>NaN with payloads and both signs, so that their order is observable in the bits.</summary>
        private static readonly double[] NaNs =
        {
            BitConverter.Int64BitsToDouble(0x7FF8000000000001), BitConverter.Int64BitsToDouble(unchecked((long)0xFFF8000000000000)),
            BitConverter.Int64BitsToDouble(0x7FF8000000000ABC), double.NaN,
        };

        /// <summary>Random samples rich in ties, signed zeros and NaN.</summary>
        private static IEnumerable<double[]> Samples()
        {
            var random = new Random(20261023);
            foreach (int n in new[] { 1, 2, 3, 7, 40, 41, 257, 1000 })
            {
                for (int repeat = 0; repeat < 6; repeat++)
                {
                    var data = new double[n];
                    for (int i = 0; i < n; i++)
                    {
                        double u = random.NextDouble();
                        data[i] = repeat switch
                        {
                            0 => 100 * random.NextDouble() - 50,
                            1 => u < 0.5 ? (random.Next(2) == 0 ? -0.0 : 0.0) : random.Next(-3, 4),
                            2 => u < 0.2 ? NaNs[random.Next(NaNs.Length)] : u < 0.6 ? (random.Next(2) == 0 ? -0.0 : 0.0) : random.Next(-2, 3),
                            3 => u < 0.1 ? double.NegativeInfinity : u < 0.2 ? double.PositiveInfinity : Math.Round(10 * u) / 10,
                            4 => -0.0,
                            _ => u < 0.3 ? 0.0 : u < 0.6 ? -0.0 : u < 0.8 ? -1E-300 : 1E-300,
                        };
                    }
                    yield return data;
                }
            }
        }

        /// <summary>The percentiles asked for, both ends included.</summary>
        private static readonly double[] K = { 0.0, 1E-9, 0.05, 0.1, 0.25, 0.333, 0.5, 0.75, 0.9, 0.95, 0.999, 1.0 };

        /// <summary>Type 7 percentiles (scalar and list) are bitwise those of the stable sort.</summary>
        [TestMethod]
        public void TypeSevenMatchesStableSortBitwise()
        {
            foreach (var data in Samples())
            {
                var stable = data.OrderBy(x => x).ToArray();
                var list = Stats.Percentile(data, K);
                var reference = Stats.Percentile(stable, K, true);
                for (int i = 0; i < K.Length; i++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference[i]), BitConverter.DoubleToInt64Bits(list[i]), $"n={data.Length} k={K[i]}");
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Percentile(stable, K[i], true)),
                        BitConverter.DoubleToInt64Bits(Stats.Percentile(data, K[i])), $"n={data.Length} k={K[i]}");
                }
            }
        }

        /// <summary>Plotting-position percentiles (scalar and list, every formula) are bitwise those of the stable sort.</summary>
        [TestMethod]
        public void PlottingPositionPercentilesMatchStableSortBitwise()
        {
            var kinds = (PlottingPositions.PlottingPostionType[])Enum.GetValues(typeof(PlottingPositions.PlottingPostionType));
            foreach (var data in Samples())
            {
                var stable = data.OrderBy(x => x).ToArray();
                foreach (var kind in kinds)
                {
                    var list = Stats.Percentile(data, K, kind);
                    var reference = Stats.Percentile(stable, K, kind, true);
                    for (int i = 0; i < K.Length; i++)
                    {
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference[i]), BitConverter.DoubleToInt64Bits(list[i]), $"{kind} n={data.Length} k={K[i]}");
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Percentile(stable, K[i], kind, true)),
                            BitConverter.DoubleToInt64Bits(Stats.Percentile(data, K[i], kind)), $"{kind} n={data.Length} k={K[i]}");
                    }
                }
            }
        }

        /// <summary>The input list is left unchanged by the unsorted path.</summary>
        [TestMethod]
        public void InputIsNotModified()
        {
            var data = new[] { 3.0, -0.0, double.NaN, 0.0, 1.0, -2.0 };
            var copy = (double[])data.Clone();
            Stats.Percentile(data, 0.5);
            Stats.Percentile(data, K);
            Stats.Percentile(data, K, PlottingPositions.PlottingPostionType.Weibull);
            for (int i = 0; i < data.Length; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(copy[i]), BitConverter.DoubleToInt64Bits(data[i]));
        }

        /// <summary>NumPy oracle: <c>numpy.percentile(data, [10, 25, 50, 90], method="linear")</c>.</summary>
        [TestMethod]
        public void KnownValuesMatchNumPy()
        {
            var data = new List<double> { 3, 1, 4, 1, 5, 9, 2, 6, 5, 3 };
            var result = Stats.Percentile(data, new[] { 0.1, 0.25, 0.5, 0.9 });
            // NumPy interpolates as b - (b - a)(1 - t) above t = 0.5, so its 90th percentile is one unit in the last place
            // below the a + t(b - a) of Type 7 here; the tolerance covers that ordering difference.
            var expected = new[] { 1.0, 2.25, 3.5, 6.299999999999999 };
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], result[i], 1E-12);
        }
    }
}
