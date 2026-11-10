#nullable disable
using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// The parametric-bootstrap unimodality test check, shared by the unit test and the full sweep in
    /// Benchmarks_Numerics.
    /// </summary>
    public partial class Test_HypothesisTests
    {
        /// <summary>
        /// Asserts that the parametric-bootstrap p-value of the Gaussian mixture likelihood-ratio test is reproducible
        /// for a fixed seed, large for the unimodal fixture and at its minimum 1/(B + 1) for the bimodal fixture, and
        /// that the chi-square overload returns its reference value.
        /// </summary>
        /// <param name="realizations">The number of bootstrap realizations B.</param>
        public static void UnimodalityTestBootstrap(int realizations)
        {
            var unimodalData = new double[] { 4.5, 5.2, 5.1, 4.9, 5.0, 5.3, 5.4, 4.8, 4.7, 5.2,
                   5.1, 4.6, 5.0, 5.3, 5.1, 4.9, 5.2, 5.0, 4.8, 5.3,
                   4.9, 5.1, 5.2, 4.8, 5.0, 5.1, 5.2, 4.7, 5.3, 5.0 };
            var bimodalData = new double[] { 3.8, 3.9, 4.0, 3.9, 4.0, 4.1, 4.0, 3.8, 3.9, 4.0,
                4.1, 4.0, 3.9, 4.0, 4.0, 4.3, 4.4, 4.4, 4.5, 4.4, 4.3,
                4.4, 4.5, 4.4, 4.3, 4.4, 4.4, 4.5, 4.4, 4.3 };

            double pUni = HypothesisTests.UnimodalityTest(unimodalData, realizations);
            Assert.AreEqual(pUni, HypothesisTests.UnimodalityTest(unimodalData, realizations), 0d);
            Assert.IsGreaterThan(0.2, pUni);
            Assert.AreEqual(1d / (realizations + 1), HypothesisTests.UnimodalityTest(bimodalData, realizations), 1E-15);

            Assert.AreEqual(0.4142441, HypothesisTests.UnimodalityTest(unimodalData), 1E-4);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => HypothesisTests.UnimodalityTest(unimodalData, 0));
            Assert.ThrowsExactly<ArgumentException>(() => HypothesisTests.UnimodalityTest(new double[] { 1, 2, 3 }, 99));
        }
    }
}
