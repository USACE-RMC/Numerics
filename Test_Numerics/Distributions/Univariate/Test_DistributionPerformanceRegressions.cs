using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Allocation regressions for distribution numerical paths used in repeated fitting work.</summary>
    [TestClass]
    public class Test_DistributionPerformanceRegressions
    {
        /// <summary>Receives diagnostic allocation output from the focused tests.</summary>
        public TestContext TestContext { get; set; }

#if NET8_0_OR_GREATER
        /// <summary>The small-shape gamma series reuses its immutable zeta constants.</summary>
        [TestMethod]
        public void GammaLogTails_RepeatedCallsDoNotAllocateZetaTables()
        {
            double checksum = 0d;
            foreach (double shape in new[] { 0.25d, 0.5d, 2d })
                checksum += DistributionNumerics.GammaLogCDF(shape, 1d);

            long before = GC.GetAllocatedBytesForCurrentThread();
            for (int repetition = 0; repetition < 100; repetition++)
            {
                checksum += DistributionNumerics.GammaLogCDF(0.25d, 1d);
                checksum += DistributionNumerics.GammaLogCDF(0.5d, 1d);
                checksum += DistributionNumerics.GammaLogCDF(2d, 1d);
            }
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

            Assert.IsTrue(double.IsFinite(checksum));
            Assert.AreEqual(0L, allocated, $"Repeated gamma log-tail evaluation allocated {allocated} bytes.");
        }
#endif

        /// <summary>Preserves the original literals, ordered runtime sums, fallback boundary and invalid-index behavior.</summary>
        [TestMethod]
        public void ZetaInteger_PrecomputedRangeRetainsBitwiseValues()
        {
            double[] literals = { 1.6449340668482264365, 1.2020569031595942854, 1.0823232337111381915,
                1.0369277551433699263, 1.0173430619844491397, 1.0083492773819228268,
                1.0040773561979443394, 1.0020083928260822144, 1.0009945751278180853,
                1.0004941886041194646, 1.0002460865533080483, 1.0001227133475784891,
                1.0000612481350587048, 1.0000305882363070205, 1.0000152822594086519 };
            for (int n = 2; n <= 16; n++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(literals[n - 2]),
                    BitConverter.DoubleToInt64Bits(DistributionNumerics.ZetaInteger(n)), $"n={n}");
            }
            for (int n = 17; n <= 60; n++)
            {
                double expected = 1d;
                for (int k = 2; k <= 32; k++) expected += Math.Pow(k, -n);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected),
                    BitConverter.DoubleToInt64Bits(DistributionNumerics.ZetaInteger(n)), $"n={n}");
            }
            Assert.ThrowsExactly<IndexOutOfRangeException>(() => DistributionNumerics.ZetaInteger(1));
        }

#if NET8_0_OR_GREATER
        /// <summary>LP3 gradients solve each standardized quantile once across skew, base, tail and scale cases.</summary>
        [TestMethod]
        public void LogPearsonQuantileGradient_RepeatedCallsAvoidDuplicateInverseWork()
        {
            var cases = new[]
            {
                (Distribution: new LogPearsonTypeIII(0.3d, 0.2d, -0.8d) { Base = Math.E }, Probability: 1E-10),
                (Distribution: new LogPearsonTypeIII(0.3d, 0.2d, 0d) { Base = 2d }, Probability: 0.83d),
                (Distribution: new LogPearsonTypeIII(3d, 0.35d, 0.8d) { Base = 10d }, Probability: 1d - 1E-10),
                (Distribution: new LogPearsonTypeIII(-250d, 1E-6d, 4d) { Base = Math.E }, Probability: 0.5d)
            };
            double checksum = 0d;
            foreach (var item in cases)
                foreach (double value in item.Distribution.QuantileGradient(item.Probability)) checksum += value;

            long before = GC.GetAllocatedBytesForCurrentThread();
            const int repetitions = 20;
            for (int repetition = 0; repetition < repetitions; repetition++)
                foreach (var item in cases)
                    foreach (double value in item.Distribution.QuantileGradient(item.Probability)) checksum += value;
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

            Assert.IsTrue(double.IsFinite(checksum));
            TestContext.WriteLine($"Measured LP3 gradient allocation: {allocated} bytes for {repetitions * cases.Length} calls.");
            // Each repaired call retains one Pearson object and its returned gradient array (about 104 bytes).
            // A 128-byte ceiling allows modest runtime layout variation but rejects the duplicate path's extra Pearson object.
            Assert.IsLessThanOrEqualTo(128L * repetitions * cases.Length, allocated,
                $"Repeated LP3 gradient evaluation allocated {allocated} bytes.");
        }
#endif
    }
}
