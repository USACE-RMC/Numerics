using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for accepted-interval quadrature records, including domain mass, weighted
    /// values, stratified integration, and the disabled-recorder path.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_AdaptiveGaussKronrodRecorder
    {
        /// <summary>A sharp peak that forces deep adaptive subdivision on [0, 1].</summary>
        private static double SharpPeak(double x)
        {
            double t = (x - 0.3) * 60d;
            return Math.Exp(-t * t) + 0.1 * x;
        }

        /// <summary>
        /// Verifies that accepted-interval weights partition the domain and their weighted
        /// function values reproduce the integral after adaptive subdivision.
        /// </summary>
        [TestMethod]
        public void Test_Recorder_MassAndResultIdentities()
        {
            var records = new List<(double X, double Weight, double Value)>();
            var gk = new AdaptiveGaussKronrod(SharpPeak, 0d, 1d)
            {
                RelativeTolerance = 1E-10,
                MinDepth = 2,
                Recorder = (x, w, fx) => records.Add((x, w, fx)),
            };
            gk.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, gk.Status);

            // Records arrive as 21 nodes per accepted interval, and subdivision must have
            // happened for this peak at this tolerance.
            Assert.AreEqual(0, records.Count % 21, "Records contain complete 21-node intervals.");
            int intervals = records.Count / 21;
            Assert.IsGreaterThan(1, intervals, "The sharp peak must force subdivision.");

            double weightSum = 0d;
            double weightedValueSum = 0d;
            for (int i = 0; i < records.Count; i++)
            {
                Assert.IsGreaterThan(0d, records[i].X, "G10K21 nodes are strictly interior (lower end).");
                Assert.IsLessThan(1d, records[i].X, "G10K21 nodes are strictly interior (upper end).");
                Assert.IsGreaterThan(0d, records[i].Weight, "Kronrod weights are strictly positive.");
                weightSum += records[i].Weight;
                weightedValueSum += records[i].Weight * records[i].Value;
            }

            Assert.AreEqual(1d, weightSum, 1E-12, "The accepted-interval weights partition the domain exactly.");
            Assert.AreEqual(gk.Result, weightedValueSum, Math.Abs(gk.Result) * 1E-12,
                "The weighted node sum reproduces the returned integral (floating-point reassociation only).");

            // The function-value channel echoes the integrand exactly at each node.
            Assert.AreEqual(SharpPeak(records[0].X), records[0].Value, 0d);
        }

        /// <summary>
        /// Test the stratified-bin form: the records cover the union of the bins and reproduces
        /// the summed result, with per-bin acceptance.
        /// </summary>
        [TestMethod]
        public void Test_Recorder_StratifiedBins()
        {
            var bins = new List<StratificationBin>
            {
                new StratificationBin(0d, 0.25d),
                new StratificationBin(0.25d, 0.6d),
                new StratificationBin(0.6d, 1d),
            };
            var records = new List<(double X, double Weight, double Value)>();
            var gk = new AdaptiveGaussKronrod(SharpPeak, 0d, 1d)
            {
                RelativeTolerance = 1E-10,
                Recorder = (x, w, fx) => records.Add((x, w, fx)),
            };
            gk.Integrate(bins);
            Assert.AreEqual(IntegrationStatus.Success, gk.Status);

            double weightSum = 0d;
            double weightedValueSum = 0d;
            for (int i = 0; i < records.Count; i++)
            {
                weightSum += records[i].Weight;
                weightedValueSum += records[i].Weight * records[i].Value;
            }
            Assert.AreEqual(1d, weightSum, 1E-12, "The records cover the union of the stratification bins.");
            Assert.AreEqual(gk.Result, weightedValueSum, Math.Abs(gk.Result) * 1E-12);
        }

        /// <summary>
        /// Verifies that recording does not change the result or function-evaluation count.
        /// </summary>
        [TestMethod]
        public void Test_Recorder_OffIsByteIdentical()
        {
            var plain = new AdaptiveGaussKronrod(SharpPeak, 0d, 1d) { RelativeTolerance = 1E-10 };
            plain.Integrate();

            var recorded = new AdaptiveGaussKronrod(SharpPeak, 0d, 1d)
            {
                RelativeTolerance = 1E-10,
                Recorder = (x, w, fx) => { },
            };
            recorded.Integrate();

            Assert.AreEqual(BitConverter.DoubleToInt64Bits(plain.Result), BitConverter.DoubleToInt64Bits(recorded.Result),
                "Recording must not perturb the integration result.");
            Assert.AreEqual(plain.FunctionEvaluations, recorded.FunctionEvaluations,
                "Recording must not change the evaluation count.");
        }
    }
}
