using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Regression coverage for allocation-free LP3 log-density delegation.</summary>
    [TestClass]
    public class Test_LogPearsonPerformanceRegressions
    {
        /// <summary>Receives diagnostic allocation output from the focused tests.</summary>
        public TestContext TestContext { get; set; }

#if NET8_0_OR_GREATER
        /// <summary>Repeated valid LP3 log-density evaluation does not allocate Pearson wrappers.</summary>
        [TestMethod]
        public void LogPDF_RepeatedCallsDoNotAllocatePearsonWrappers()
        {
            var cases = new[]
            {
                (Distribution: new LogPearsonTypeIII(3d, 0.35d, 0d) { Base = 10d }, X: 700d),
                (Distribution: new LogPearsonTypeIII(3d, 0.35d, 1E-5d) { Base = 2d }, X: 8d),
                (Distribution: new LogPearsonTypeIII(3d, 0.35d, -0.5d) { Base = 10d }, X: 1200d),
                (Distribution: new LogPearsonTypeIII(0.3d, 0.2d, 1.5d) { Base = Math.E }, X: 2d),
                (Distribution: new LogPearsonTypeIII(3d, 0.35d, 4d) { Base = 10d }, X: 900d)
            };
            long checksum = 0L;
            foreach (var item in cases) checksum ^= BitConverter.DoubleToInt64Bits(item.Distribution.LogPDF(item.X));

            long before = GC.GetAllocatedBytesForCurrentThread();
            const int repetitions = 100;
            for (int repetition = 0; repetition < repetitions; repetition++)
                foreach (var item in cases)
                    checksum ^= BitConverter.DoubleToInt64Bits(item.Distribution.LogPDF(item.X));
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

            TestContext.WriteLine($"Measured LP3 LogPDF allocation: {allocated} bytes for {repetitions * cases.Length} calls; checksum={checksum}.");
            Assert.AreEqual(0L, allocated, $"Repeated LP3 LogPDF evaluation allocated {allocated} bytes.");
        }
#endif

    }
}
