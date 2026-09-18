using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Regression contracts for shared censored log likelihoods.</summary>
    [TestClass]
    public class Test_DistributionLikelihoodRegression
    {
        /// <summary>R 4.4.3 pnorm(log.p=TRUE) oracle for a representable upper-tail interval.</summary>
        [TestMethod]
        public void NormalTailIntervalRemainsFinite()
        {
            var distribution = new Normal(0, 1);
            Assert.AreEqual(-43.628216632280818, distribution.LogLikelihood_Intervals(9, 10), 2E-13);
            Assert.AreEqual(-43.628216632280818, distribution.LogLikelihood_Intervals(-10, -9), 2E-13);
            Assert.AreEqual(-46.970640393085588, distribution.LogLikelihood_Intervals(0, 1E-20), 2E-13);
        }

        /// <summary>An empty censoring category contributes zero even at an impossible threshold.</summary>
        [TestMethod]
        public void ZeroCensoringCountsContributeZero()
        {
            var distribution = new Exponential(0, 1);
            Assert.AreEqual(0, distribution.LogLikelihood_LeftCensored(-1, 0));
            Assert.AreEqual(0, distribution.LogLikelihood_RightCensored(double.PositiveInfinity, 0));
        }

        /// <summary>Intervals are empty when their limits coincide, including infinities.</summary>
        [TestMethod]
        public void EqualIntervalBoundsHaveZeroProbability()
        {
            var distribution = new Normal(0, 1);
            Assert.AreEqual(double.NegativeInfinity, distribution.LogLikelihood_Intervals(0, 0));
            Assert.AreEqual(double.NegativeInfinity, distribution.LogLikelihood_Intervals(double.PositiveInfinity, double.PositiveInfinity));
        }

        /// <summary>Invalid censoring counts and reversed intervals are rejected explicitly.</summary>
        [TestMethod]
        public void InvalidCensoringInputsAreRejected()
        {
            var distribution = new Normal(0, 1);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.LogLikelihood_LeftCensored(0, -1));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.LogLikelihood_RightCensored(0, -1));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.LogLikelihood_Intervals(1, 0));
        }
    }
}
