using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Verifies lower-tail accuracy shared by the logarithmic normal-family distributions.
    /// </summary>
    [TestClass]
    public class Test_LogDistributionFarTails
    {
        private static readonly double[] Probabilities = { 1E-17d, 1E-40d, 1E-100d };
        private static readonly double[] BaseTenQuantiles =
        {
            3.207796253986819E-9d,
            4.887408373653465E-14d,
            5.327781911580459E-22d,
        };
        private static readonly double[] NaturalQuantiles =
        {
            2.0473517891357346E-4d,
            1.6563049483082078E-6d,
            5.768415132086823E-10d,
        };

        /// <summary>
        /// Verifies the supplied distribution against externally computed lower-tail quantiles and
        /// checks CDF/inverse-CDF round trips where the cancelling error-function form returns zero.
        /// </summary>
        /// <param name="distribution">The distribution under test.</param>
        /// <param name="quantiles">The externally computed quantiles.</param>
        private static void AssertFarLowerTail(UnivariateDistributionBase distribution, double[] quantiles)
        {
            for (int i = 0; i < Probabilities.Length; i++)
            {
                double probability = Probabilities[i];
                double quantile = quantiles[i];
                double computedQuantile = distribution.InverseCDF(probability);

                Assert.AreEqual(quantile, computedQuantile, quantile * 2E-8d,
                    $"InverseCDF at p={probability:R}.");
                Assert.AreEqual(probability, distribution.CDF(quantile), probability * 2E-9d,
                    $"CDF at the reference quantile for p={probability:R}.");
                Assert.AreEqual(probability, distribution.CDF(computedQuantile), probability * 2E-8d,
                    $"CDF/InverseCDF round trip at p={probability:R}.");
            }
        }

        /// <summary>
        /// Pins base-10 LogNormal lower-tail quantiles computed with mpmath 1.4.1 at 100 digits.
        /// </summary>
        [TestMethod]
        public void Test_LogNormal_FarLowerTail()
        {
            AssertFarLowerTail(new LogNormal(0d, 1d), BaseTenQuantiles);
        }

        /// <summary>
        /// Pins natural-log LnNormal lower-tail quantiles computed with mpmath 1.4.1 at 100 digits.
        /// </summary>
        [TestMethod]
        public void Test_LnNormal_FarLowerTail()
        {
            AssertFarLowerTail(new LnNormal { Mu = 0d, Sigma = 1d }, NaturalQuantiles);
        }

        /// <summary>
        /// Pins the zero-skew, base-10 LogPearson III limit against the same normal-tail oracle.
        /// </summary>
        [TestMethod]
        public void Test_LogPearsonTypeIII_ZeroSkewFarLowerTail()
        {
            AssertFarLowerTail(new LogPearsonTypeIII(0d, 1d, 0d), BaseTenQuantiles);
        }
    }
}
