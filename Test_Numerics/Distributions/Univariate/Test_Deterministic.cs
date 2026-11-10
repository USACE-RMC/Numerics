using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Deterministic distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Deterministic
    {
        /// <summary>
        /// A point mass has zero spread (closed form).
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var d = new Deterministic(7);
            Assert.AreEqual(0, d.StandardDeviation);
            Assert.AreEqual(0, d.Variance);
            Assert.AreEqual(7, d.Mean);
        }

        /// <summary>
        /// Setting the value revalidates the parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Value_SetterRevalidates()
        {
            var d = new Deterministic(7);
            Assert.IsTrue(d.ParametersValid);
            d.Value = double.NaN;
            Assert.IsFalse(d.ParametersValid);
            d.Value = 3;
            Assert.IsTrue(d.ParametersValid);
            Assert.AreEqual(3, d.InverseCDF(0.5));
        }

        /// <summary>
        /// InverseCDF rejects a NaN probability.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_NaNProbabilityThrows()
        {
            var d = new Deterministic(7);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => d.InverseCDF(double.NaN));
        }

        /// <summary>
        /// An invalid point value makes PDF, CDF and InverseCDF throw, as the other distributions do, and a valid value
        /// set afterwards restores them.
        /// </summary>
        [TestMethod]
        public void Test_InvalidValue_Throws()
        {
            foreach (double invalid in new[] { double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            {
                var d = new Deterministic(7) { Value = invalid };
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => d.PDF(0));
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => d.CDF(0));
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => d.InverseCDF(0.5));
                d.Value = 3;
                Assert.AreEqual(1d, d.CDF(3));
                Assert.AreEqual(3d, d.InverseCDF(0.5));
            }
        }
    }
}
