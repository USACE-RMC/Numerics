using System;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Stats = Numerics.Data.Statistics.Statistics;
using PPType = Numerics.Data.Statistics.PlottingPositions.PlottingPostionType;

namespace Data.Statistics
{
    /// <summary>
    /// Tests for <see cref="Numerics.Data.Statistics.Statistics.Percentile(System.Collections.Generic.IList{double}, double, PlottingPositions.PlottingPostionType, bool)"/>.
    /// </summary>
    [TestClass]
    public partial class Test_PlottingPositionPercentile
    {
        /// <summary>The sample {3, 1, 4, 1, 5, 9, 2, 6}; sorted {1, 1, 2, 3, 4, 5, 6, 9}.</summary>
        private static readonly double[] Sample = { 3, 1, 4, 1, 5, 9, 2, 6 };

        /// <summary>
        /// Known values of R quantile(x, k, type = 6, 5, 4) for Weibull, Hazen and Kaplan-Meier, including the
        /// held ends.
        /// </summary>
        [TestMethod]
        public void KnownValuesMatchHyndmanFanTypes()
        {
            double[] k = { 0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };
            double[] type6 = { 1.0, 1.0, 1.25, 3.5, 5.75, 9.0, 9.0 };
            double[] type5 = { 1.0, 1.0, 1.5, 3.5, 5.5, 8.1, 9.0 };
            double[] type4 = { 1.0, 1.0, 1.0, 3.0, 5.0, 6.6, 9.0 };
            for (int i = 0; i < k.Length; i++)
            {
                Assert.AreEqual(type6[i], Stats.Percentile(Sample, k[i], PPType.Weibull), 1E-12, $"Weibull k = {k[i]}");
                Assert.AreEqual(type5[i], Stats.Percentile(Sample, k[i], PPType.Hazen), 1E-12, $"Hazen k = {k[i]}");
                Assert.AreEqual(type4[i], Stats.Percentile(Sample, k[i], PPType.KaplanMeier), 1E-12, $"Kaplan-Meier k = {k[i]}");
            }
        }

        /// <summary>
        /// At each plotting position the percentile returns the order statistic, for every formula, and the array
        /// overload equals the scalar one bit for bit with sorted or unsorted data.
        /// </summary>
        [TestMethod]
        public void PlottingPositionsMapToOrderStatistics()
        {
            var sorted = (double[])Sample.Clone();
            Array.Sort(sorted);
            foreach (PPType type in Enum.GetValues(typeof(PPType)))
            {
                var positions = PlottingPositions.Function(Sample.Length, type)!;
                for (int i = 0; i < positions.Length; i++)
                    Assert.AreEqual(sorted[i], Stats.Percentile(Sample, positions[i], type), 1E-12, $"{type}, i = {i}");

                var k = new[] { 0.0, 0.05, 0.3, 0.5, 0.77, 0.99, 1.0 };
                var array = Stats.Percentile(Sample, k, type);
                var arraySorted = Stats.Percentile(sorted, k, type, true);
                for (int i = 0; i < k.Length; i++)
                {
                    long scalar = BitConverter.DoubleToInt64Bits(Stats.Percentile(Sample, k[i], type));
                    Assert.AreEqual(scalar, BitConverter.DoubleToInt64Bits(array[i]));
                    Assert.AreEqual(scalar, BitConverter.DoubleToInt64Bits(arraySorted[i]));
                }
            }
        }

        /// <summary>
        /// Hand-computed values for the Gringorten and Filliben median formulas, the two-observation boundary, the
        /// extremes at k = 0 and 1 for every formula, and agreement of sorted and unsorted scalar calls.
        /// </summary>
        [TestMethod]
        public void HandValuesAndBoundaries()
        {
            // Gringorten: h = (9 − 0.88) k + 0.44; k = 0.5 → h = 4.5 → 3.5; k = 0.9 → h = 7.748 → 6 + 0.748 · 3.
            Assert.AreEqual(3.5, Stats.Percentile(Sample, 0.5, PPType.Gringorten), 1E-12);
            Assert.AreEqual(8.244, Stats.Percentile(Sample, 0.9, PPType.Gringorten), 1E-12);
            // Median (alpha = 0.3175): k = 0.25 → h = 8.365 · 0.25 + 0.3175 = 2.40875 → 1 + 0.40875 · 1.
            Assert.AreEqual(1.40875, Stats.Percentile(Sample, 0.25, PPType.Median), 1E-12);

            var pair = new[] { 10.0, 2.0 };
            Assert.AreEqual(6.0, Stats.Percentile(pair, 0.5, PPType.Weibull), 1E-12);
            Assert.AreEqual(2.0, Stats.Percentile(pair, 1.0 / 3.0, PPType.Weibull), 1E-12);
            Assert.AreEqual(10.0, Stats.Percentile(pair, 0.9, PPType.Weibull));

            var sorted = (double[])Sample.Clone();
            Array.Sort(sorted);
            foreach (PPType type in Enum.GetValues(typeof(PPType)))
            {
                Assert.AreEqual(1.0, Stats.Percentile(Sample, 0.0, type), $"{type}");
                Assert.AreEqual(9.0, Stats.Percentile(Sample, 1.0, type), $"{type}");
                foreach (double k in new[] { 0.13, 0.5, 0.87 })
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Percentile(Sample, k, type)),
                        BitConverter.DoubleToInt64Bits(Stats.Percentile(sorted, k, type, true)), $"{type}, k = {k}");
            }
        }

        /// <summary>A single value is returned for every k; invalid inputs are refused.</summary>
        [TestMethod]
        public void EdgeCasesAndArgumentErrors()
        {
            Assert.AreEqual(7.0, Stats.Percentile(new[] { 7.0 }, 0.3, PPType.Weibull));
            Assert.IsEmpty(Stats.Percentile(Sample, new double[0], PPType.Weibull));

            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.Percentile(null!, 0.5, PPType.Weibull));
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.Percentile(Sample, (double[])null!, PPType.Weibull));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(new double[0], 0.5, PPType.Weibull));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(new double[0], new[] { 0.5 }, PPType.Weibull));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(Sample, 1.5, PPType.Weibull));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(Sample, double.NaN, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(Sample, 0.5, (PPType)99));
        }

        /// <summary>
        /// The Weibull weight formed from k and the bracketing positions, (k − r/d) / ((r + 1)/d − r/d), is closer to
        /// exact rational evaluation than h − floor(h) at the return-period probabilities 1 − 1/T, over seeded
        /// samples with n = 5 to 120 and T = 2 to 100. The exact value is evaluated in integers (every double is an
        /// integer multiple of 2^-1074), and the error is in units in the last place of the computed value.
        /// This test runs 300 cases; the 3000-case sweep runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep PlottingPosition</c>.
        /// </summary>
        [TestMethod]
        public void WeibullWeightCloserToExactRational()
        {
            WeibullWeightCloserToExactRationalSweep(300);
        }

        /// <summary>
        /// The Weibull weight is clamped to [0, 1]: at n = 11 and k = 1 − 1/6, (n + 1) k rounds up to the integer 10 and
        /// the unclamped weight is −1.3E-15, so the percentile is exactly x(10) rather than just below it. At n = 4 the
        /// probabilities 0.2, 0.4, 0.6 and 0.8 return exactly x(1) to x(4).
        /// </summary>
        [TestMethod]
        public void WeibullWeightIsClampedAndExactAtThePositions()
        {
            var x = new double[11];
            for (int i = 0; i < 11; i++) x[i] = i + 1;
            Assert.AreEqual(10d, Stats.Percentile(x, 1d - 1d / 6d, PPType.Weibull, true), 0d);
            Assert.AreEqual(10d, Stats.Percentile(x, 1d - 1d / 6d, PPType.Weibull), 0d);

            double[] four = { 2.5, 7.25, 11.0, 19.75 };
            double[] k = { 0.2, 0.4, 0.6, 0.8 };
            for (int i = 0; i < 4; i++)
                Assert.AreEqual(four[i], Stats.Percentile(four, k[i], PPType.Weibull, true), 0d, $"k = {k[i]}");
        }

        /// <summary>
        /// For Blom and Hazen the weight formed from k and the bracketing positions, (k − p_r) / (p_(r+1) − p_r) with
        /// p_r = (r − α)/(n + 1 − 2α), gives percentiles at least as close to exact rational evaluation as h − floor(h), in
        /// mean and maximum error, at the return-period probabilities 1 − 1/T over 3000 seeded samples with n = 5 to 120
        /// and T = 2 to 100. The exact value is evaluated in integers, as in the Weibull test, with α = a/8 and the rank
        /// scaled by 8.
        /// This test runs 600 samples; the 3000-sample sweep runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep PlottingPosition</c>.
        /// </summary>
        [TestMethod]
        public void BlomAndHazenWeightCloserToExactRational()
        {
            BlomAndHazenWeightCloserToExactRationalSweep(600);
        }

        /// <summary>
        /// The Blom and Hazen weights themselves, against the exact weight h − floor(h) evaluated in integers: on 20000
        /// cases each of k uniform in (0, 1) with n = 2 to 1000, of k = j/1000 with n = 5 to 200, and of k = 1 − 1/T with
        /// n = 5 to 120, the bracket weight has a lower mean and a maximum no larger than h − floor(h). The library
        /// result is checked to be the interpolation with exactly that weight.
        /// This test runs 3000 cases of each grid; the 20000-case sweep runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep PlottingPosition</c>.
        /// </summary>
        [TestMethod]
        public void BlomAndHazenWeightCloserOnProbabilityGrids()
        {
            BlomAndHazenWeightCloserOnProbabilityGridsSweep(3000);
        }

        /// <summary>
        /// The Weibull weight is unchanged, bit for bit (k − r/d) / ((r + 1)/d − r/d) with d = n + 1, and the Median,
        /// Cunnane and Gringorten formulas keep the weight h − floor(h).
        /// This test runs 500 random cases; the 3000-case sweep runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep PlottingPosition</c>.
        /// </summary>
        [TestMethod]
        public void OtherFormulasKeepTheirWeightsBitwise()
        {
            OtherFormulasKeepTheirWeightsBitwiseSweep(500);
        }

        /// <summary>
        /// The weighted plotting-position percentile of a list of probabilities equals the scalar member bit for bit, for
        /// every formula, sorted and unsorted data, ascending and unordered probabilities with repeats and the ends, and
        /// list inputs that are not arrays.
        /// This test runs 100 random trials; the 600-trial sweep runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep PlottingPosition</c>.
        /// </summary>
        [TestMethod]
        public void WeightedListOfProbabilitiesMatchesScalarBitwise()
        {
            WeightedListOfProbabilitiesMatchesScalarBitwiseSweep(100);
        }

        /// <summary>
        /// Known values of the weighted Hazen list overload, the empty list, and the argument errors of the scalar member.
        /// </summary>
        [TestMethod]
        public void WeightedListOfProbabilitiesKnownValuesAndErrors()
        {
            double[] x = { 1, 2, 3, 4 };
            double[] unit = { 1, 1, 1, 1 };
            // Positions 0.125, 0.375, 0.625, 0.875; the k are not ascending.
            var values = Stats.Percentile(x, new[] { 0.5, 0.1, 0.25, 1.0, 0.0, 0.9, 0.5 }, unit, PPType.Hazen, true);
            double[] expected = { 2.5, 1.0, 1.5, 4.0, 1.0, 4.0, 2.5 };
            for (int j = 0; j < expected.Length; j++) Assert.AreEqual(expected[j], values[j], 1e-15, $"j = {j}");
            Assert.IsEmpty(Stats.Percentile(x, new double[0], unit, PPType.Hazen));
            Assert.IsEmpty(Stats.Percentile(new double[0], new double[0], new double[0], PPType.Hazen));

            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.Percentile(null!, new[] { 0.5 }, unit, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.Percentile(x, new[] { 0.5 }, null!, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.Percentile(x, (double[])null!, unit, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(new double[0], new[] { 0.5 }, new double[0], PPType.Hazen));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(x, new[] { 0.5 }, new double[] { 1, 1 }, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(x, new[] { 0.5 }, new double[4], PPType.Hazen));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(new[] { 1d, double.NaN, 3d, 4d }, new[] { 0.5 }, unit, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, new[] { 0.5, 1.5 }, unit, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, new[] { double.NaN }, unit, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, new[] { 0.5 }, new double[] { 1, -1, 1, 1 }, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, new[] { 0.5 }, unit, (PPType)99));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, new[] { 0.5 }, new[] { double.MaxValue, double.MaxValue, 1d, 1d }, PPType.Hazen));
        }

        /// <summary>
        /// The weighted Hazen percentile against known values: the positions (c_(i−1) + w_i / 2) / W, linear interpolation
        /// between them, the held ends, zero weights taking a position, and unit weights giving the unweighted Hazen
        /// percentile.
        /// </summary>
        [TestMethod]
        public void WeightedHazenKnownValues()
        {
            double[] x = { 1, 2, 3, 4 };
            double[] unit = { 1, 1, 1, 1 };
            // Positions 0.125, 0.375, 0.625, 0.875.
            Assert.AreEqual(2.5, Stats.Percentile(x, 0.5, unit, PPType.Hazen, true), 1e-15);
            Assert.AreEqual(1.5, Stats.Percentile(x, 0.25, unit, PPType.Hazen, true), 1e-15);
            Assert.AreEqual(1.0, Stats.Percentile(x, 0.1, unit, PPType.Hazen, true));
            Assert.AreEqual(4.0, Stats.Percentile(x, 0.9, unit, PPType.Hazen, true));
            Assert.AreEqual(1.0, Stats.Percentile(x, 0.0, unit, PPType.Hazen, true));
            Assert.AreEqual(4.0, Stats.Percentile(x, 1.0, unit, PPType.Hazen, true));
            // Weights {1, 3}: positions 0.125 and 0.625, so k = 0.5 is 10 + 10 (0.375 / 0.5).
            Assert.AreEqual(17.5, Stats.Percentile(new double[] { 20, 10 }, 0.5, new double[] { 3, 1 }, PPType.Hazen), 1e-14);
            // A zero weight takes the position c_(i-1) / W: positions 0.25, 0.5, 0.75.
            Assert.AreEqual(2.0, Stats.Percentile(new double[] { 1, 2, 3 }, 0.5, new double[] { 1, 0, 1 }, PPType.Hazen, true), 1e-15);
            foreach (double k in new[] { 0.05, 0.2, 0.37, 0.5, 0.81, 0.99 })
                Assert.AreEqual(Stats.Percentile(Sample, k, PPType.Hazen), Stats.Percentile(Sample, k, new double[] { 1, 1, 1, 1, 1, 1, 1, 1 }, PPType.Hazen), 1e-12, $"k = {k}");

            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.Percentile(null!, 0.5, unit, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(x, 0.5, new double[] { 1, 1 }, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(x, 0.5, new double[4], PPType.Hazen));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, 1.5, unit, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, 0.5, new double[] { 1, -1, 1, 1 }, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, 0.5, unit, (PPType)99));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(new[] { 1d, double.NaN, 3d, 4d }, 0.5, unit, PPType.Hazen));
            Assert.ThrowsExactly<ArgumentException>(() => Stats.Percentile(new[] { 1d, 2d, 3d, double.NaN }, 0.5, unit, PPType.Hazen, true));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.Percentile(x, 0.5, new[] { double.MaxValue, double.MaxValue, 1d, 1d }, PPType.Hazen));
        }

        /// <summary>
        /// The weighted Hazen percentile equals, bit for bit, the sequential reference below on random weighted sets with
        /// ties and zero weights, for sorted input and for unsorted input sorted stably by value; and the other formulas
        /// agree with the positions of
        /// <see cref="PlottingPositions.Function(System.Collections.Generic.IList{double}, PPType)"/>.
        /// This test runs 200 cases; the 2000-case sweep runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep PlottingPosition</c>.
        /// </summary>
        [TestMethod]
        public void WeightedHazenMatchesSequentialReferenceBitwise()
        {
            WeightedHazenMatchesSequentialReferenceBitwiseSweep(200);
        }

#if NET6_0_OR_GREATER
        /// <summary>
        /// With sorted arrays the weighted plotting-position percentile allocates nothing, for every formula.
        /// </summary>
        [TestMethod]
        public void WeightedPlottingPositionPercentileDoesNotAllocateOnSortedArrays()
        {
            double[] x = { 1, 2, 3, 4, 5, 6, 7, 8 };
            double[] w = { 0.5, 1, 0, 2, 1, 0.25, 1, 3 };
            var types = new[] { PPType.Weibull, PPType.Median, PPType.Blom, PPType.Cunnane, PPType.Gringorten, PPType.Hazen, PPType.KaplanMeier };
            double sink = 0;
            foreach (var type in types) sink += Stats.Percentile(x, 0.37, w, type, true);
            long before = GC.GetAllocatedBytesForCurrentThread();
            for (int r = 0; r < 100; r++)
                foreach (var type in types) sink += Stats.Percentile(x, 0.37, w, type, true);
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
            Assert.AreEqual(0L, allocated);
            Assert.IsGreaterThan(0d, sink);
        }
#endif

    }
}
