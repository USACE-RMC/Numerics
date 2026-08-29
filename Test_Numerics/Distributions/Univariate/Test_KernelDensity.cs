using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing Kernel Density Estimation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list> 
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_KernelDensity
    {
        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 5.1.1 Tippecanoe River Near Delphi, Indiana (Station 43) Data
        private double[] sample = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

        /// <summary>
        /// Test the KDE PDF against the R 'stats' package
        /// </summary>
        [TestMethod()]
        public void Test_KernelDensity_PDF()
        {
            var KDE = new KernelDensity(sample);
           
            // To replicate this, set the bandwidth in R to be the same
            // Results from R 
            var x1 = 2328.878;
            var x2 = 12221.33;
            var x3 = 28708.74;
            //
            var f1 = 1.04475e-05;
            var f2 = 7.417907e-05;
            var f3 = 1.845702e-07;
            //
            Assert.AreEqual(f1, KDE.PDF(x1), 1E-6);
            Assert.AreEqual(f2, KDE.PDF(x2), 1E-6);
            Assert.AreEqual(f3, KDE.PDF(x3), 1E-6);

        }

        /// <summary>
        /// Test the KDE CDF against the R 'spatstat' package
        /// </summary>
        [TestMethod()]
        public void Test_KernelDensity_CDF()
        {
            var KDE = new KernelDensity(sample);
            KDE.BoundedByData = false;

            // To replicate this, set the bandwidth in R to be the same
            // Results from R 
            var x1 = 2328.878;
            var x2 = 12221.33;
            var x3 = 28708.74;
            //
            var f1 = 0.01734183;
            var f2 = 0.4669572;
            var f3 = 0.9999118;
            // Test CDF
            Assert.AreEqual(f1, KDE.CDF(x1), 1E-2);
            Assert.AreEqual(f2, KDE.CDF(x2), 1E-2);
            Assert.AreEqual(f3, KDE.CDF(x3), 1E-2);
            // Test inverse CDF
            Assert.AreEqual(2328.878, KDE.InverseCDF(KDE.CDF(x1)), 1E-2);
            Assert.AreEqual(12221.33, KDE.InverseCDF(KDE.CDF(x2)), 1E-2);
            Assert.AreEqual(28708.74, KDE.InverseCDF(KDE.CDF(x3)), 1E-2);

        }
        /// <summary>
        /// Verifies that a constant nonzero sample yields a near-point-mass bandwidth proportional
        /// to the constant's magnitude rather than a fabricated spread.
        /// </summary>
        [TestMethod]
        public void ConstantNonzeroSample_UsesNearPointMassBandwidth()
        {
            double[] constantSample = { 5d, 5d, 5d, 5d };
            double expected = 5d * KernelDensity.DegenerateRelativeBandwidth;

            var distribution = new KernelDensity(constantSample);

            Assert.AreEqual(expected, distribution.Bandwidth, 0d);
            Assert.IsTrue(Tools.IsFinite(distribution.PDF(5d)));
            Assert.IsGreaterThan(0d, distribution.PDF(5d));
        }
        /// <summary>
        /// Verifies that a single finite observation yields a near-point-mass bandwidth at the
        /// observed value: a one-point sample supports no dispersion estimate.
        /// </summary>
        [TestMethod]
        public void SingleObservation_UsesNearPointMassBandwidth()
        {
            double[] sample = { 7d };
            double expected = 7d * KernelDensity.DegenerateRelativeBandwidth;

            var distribution = new KernelDensity(sample);

            Assert.AreEqual(expected, distribution.Bandwidth, 0d);
            Assert.IsTrue(Tools.IsFinite(distribution.PDF(7d)));
            Assert.IsGreaterThan(0d, distribution.PDF(7d));
        }

        /// <summary>
        /// Verifies that an all-zero sample yields the absolute degenerate bandwidth: the constant
        /// supplies no magnitude, and no spread is invented for it.
        /// </summary>
        [TestMethod]
        public void ConstantZeroSample_UsesAbsoluteDegenerateBandwidth()
        {
            double[] constantSample = { 0d, 0d, 0d, 0d };

            var distribution = new KernelDensity(constantSample);

            Assert.AreEqual(KernelDensity.DegenerateAbsoluteBandwidth, distribution.Bandwidth, 0d);
            Assert.IsTrue(Tools.IsFinite(distribution.PDF(0d)));
            Assert.IsGreaterThan(0d, distribution.PDF(0d));
        }

        /// <summary>
        /// A constant whose relative automatic bandwidth is subnormal uses the absolute
        /// degenerate fallback, avoiding an infinite density caused by reciprocal overflow.
        /// </summary>
        [TestMethod]
        public void ConstantWithSubnormalDerivedBandwidth_UsesAbsoluteDegenerateBandwidth()
        {
            double[] constantSample = { 1E-300, 1E-300, 1E-300, 1E-300 };

            var distribution = new KernelDensity(constantSample);

            Assert.AreEqual(KernelDensity.DegenerateAbsoluteBandwidth, distribution.Bandwidth, 0d);
            Assert.IsTrue(Tools.IsFinite(distribution.PDF(1E-300)));
            Assert.IsGreaterThan(0d, distribution.PDF(1E-300));
        }

        /// <summary>
        /// Verifies that a weighted constant sample yields the near-point-mass bandwidth from the
        /// constant's magnitude, independent of the weights.
        /// </summary>
        [TestMethod]
        public void WeightedConstantSample_UsesNearPointMassBandwidth()
        {
            double[] constantSample = { -3d, -3d, -3d, -3d };
            double[] weights = { 1d, 2d, 3d, 4d };
            double expected = 3d * KernelDensity.DegenerateRelativeBandwidth;

            var distribution = new KernelDensity(constantSample, weights);

            Assert.AreEqual(expected, distribution.Bandwidth, 0d);
            Assert.IsTrue(Tools.IsFinite(distribution.PDF(-3d)));
            Assert.IsGreaterThan(0d, distribution.PDF(-3d));
        }

        /// <summary>
        /// Verifies that an explicitly supplied zero bandwidth remains invalid.
        /// </summary>
        [TestMethod]
        public void ExplicitZeroBandwidth_RemainsInvalid()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                new KernelDensity(new[] { -1d, 0d, 1d }, KernelDensity.KernelType.Gaussian, 0d));
        }

        /// <summary>
        /// The automatic-bandwidth floor does not alter a positive explicit bandwidth.
        /// </summary>
        [TestMethod]
        public void ExplicitSubnormalBandwidth_RemainsAsSupplied()
        {
            var distribution = new KernelDensity(
                new[] { -1d, 0d, 1d }, KernelDensity.KernelType.Gaussian, double.Epsilon);

            Assert.AreEqual(double.Epsilon, distribution.Bandwidth, 0d);
        }

        /// <summary>
        /// Creates a deterministic sample large enough for the density reduction to be partitioned
        /// across several threads.
        /// </summary>
        /// <returns>A fixed Normal sample.</returns>
        private static double[] CreateReductionSample()
        {
            return new Normal(10d, 3d).GenerateRandomValues(2000, 20250824);
        }

        /// <summary>
        /// Verifies that repeated density evaluations at the same ordinate return bit-identical values.
        /// </summary>
        /// <remarks>
        /// The density is a sum over every sample point, so a scheduler-dependent accumulation order lets
        /// the last bits of a public distribution function change between otherwise identical calls.
        /// </remarks>
        [TestMethod]
        public void PDF_IsBitReproducibleAcrossRepeatedCalls()
        {
            var reductionSample = CreateReductionSample();
            var weights = new double[reductionSample.Length];
            for (int i = 0; i < weights.Length; i++) weights[i] = 0.5d + i % 7;

            var unweighted = new KernelDensity(reductionSample);
            var weighted = new KernelDensity(reductionSample, weights);

            foreach (double x in new[] { 2.5d, 7d, 10d, 13.25d, 19d })
            {
                long expectedUnweighted = BitConverter.DoubleToInt64Bits(unweighted.PDF(x));
                long expectedWeighted = BitConverter.DoubleToInt64Bits(weighted.PDF(x));
                for (int trial = 0; trial < 25; trial++)
                {
                    Assert.AreEqual(expectedUnweighted, BitConverter.DoubleToInt64Bits(unweighted.PDF(x)),
                        $"The unweighted density is not bit-identical at x = {x}.");
                    Assert.AreEqual(expectedWeighted, BitConverter.DoubleToInt64Bits(weighted.PDF(x)),
                        $"The weighted density is not bit-identical at x = {x}.");
                }
            }
        }

        /// <summary>
        /// Verifies that the mode and the cached CDF are bit-identical across separate instances built
        /// from the same sample. Both are derived from the density, so a drifting density moves them too.
        /// </summary>
        [TestMethod]
        public void ModeAndCDF_AreBitReproducibleAcrossInstances()
        {
            var reductionSample = CreateReductionSample();
            long expectedMode = 0L;
            long expectedCDF = 0L;

            for (int trial = 0; trial < 3; trial++)
            {
                var distribution = new KernelDensity(reductionSample);
                long mode = BitConverter.DoubleToInt64Bits(distribution.Mode);
                long cdf = BitConverter.DoubleToInt64Bits(distribution.CDF(11.5d));

                if (trial == 0)
                {
                    expectedMode = mode;
                    expectedCDF = cdf;
                    continue;
                }

                Assert.AreEqual(expectedMode, mode, "The mode is not bit-identical across instances.");
                Assert.AreEqual(expectedCDF, cdf, "The CDF is not bit-identical across instances.");
            }
        }

        /// <summary>
        /// Kernel density does not compose the empirical distribution — it carries its own
        /// lookup table with hard support gates — so the empirical extrapolation property cannot
        /// reach it: the gates and the serialized form are unchanged.
        /// </summary>
        [TestMethod]
        public void Test_KernelDensity_EmpiricalUnderTheHood_NoRegression()
        {
            var sample = new double[] { 8d, 9d, 10d, 11d, 12d, 10.5d, 9.5d, 10.2d };
            var kde = new KernelDensity(sample);
            Assert.AreEqual(0d, kde.CDF(kde.Minimum - 1d));
            Assert.AreEqual(1d, kde.CDF(kde.Maximum + 1d));
            Assert.DoesNotContain("Extrapolation", kde.ToXElement().ToString());
        }





    }
}
