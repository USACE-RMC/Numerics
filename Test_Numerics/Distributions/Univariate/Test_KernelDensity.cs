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





    }
}
