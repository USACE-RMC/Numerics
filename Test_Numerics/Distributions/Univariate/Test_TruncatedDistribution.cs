using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the general Truncated distribution class.
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
    public class Test_TruncatedDistribution
    {

        /// <summary>
        /// This method was verified against the R package "truncnorm"
        /// </summary>
        [TestMethod]
        public void Test_TruncatedNormalDist()
        {
            var tn = new TruncatedDistribution(new Normal(2, 1), 1.10, 2.11);
            var d = tn.PDF(1.5);
            var p = tn.CDF(1.5);
            var q = tn.InverseCDF(p);

            Assert.AreEqual(0.9786791, d, 1E-5);
            Assert.AreEqual(0.3460251, p, 1E-5);
            Assert.AreEqual(1.5, q, 1E-5);

            tn = new TruncatedDistribution(new Normal(10, 3), 8, 25);
            d = tn.PDF(12.75);
            p = tn.CDF(12.75);
            q = tn.InverseCDF(p);

            Assert.AreEqual(0.1168717, d, 1E-5);
            Assert.AreEqual(0.7596566, p, 1E-5);
            Assert.AreEqual(12.75, q, 1E-5);

            tn = new TruncatedDistribution(new Normal(0, 3), 0, 9);
            d = tn.PDF(4.5);
            p = tn.CDF(4.5);
            q = tn.InverseCDF(p);

            Assert.AreEqual(0.08657881, d, 1E-5);
            Assert.AreEqual(0.868731, p, 1E-5);
            Assert.AreEqual(4.5, q, 1E-5);

        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(0.5, ((Normal)tn.BaseDistribution).Mu);
            Assert.AreEqual(0.2, ((Normal)tn.BaseDistribution).Sigma);
            Assert.AreEqual(0, tn.Min);
            Assert.AreEqual(1, tn.Max);

            var tn2 = new TruncatedDistribution(new Normal(1, 1), 1, 2);
            Assert.AreEqual(1, ((Normal)tn2.BaseDistribution).Mu);
            Assert.AreEqual(1, ((Normal)tn2.BaseDistribution).Sigma);
            Assert.AreEqual(1, tn2.Min);
            Assert.AreEqual(2, tn2.Max);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual("Mean (µ)", tn.ParametersToString[0, 0]);
            Assert.AreEqual("Std Dev (σ)", tn.ParametersToString[1, 0]);
            Assert.AreEqual("Min", tn.ParametersToString[2, 0]);
            Assert.AreEqual("Max", tn.ParametersToString[3, 0]);
            Assert.AreEqual("0.5", tn.ParametersToString[0, 1]);
            Assert.AreEqual("0.2", tn.ParametersToString[1, 1]);
            Assert.AreEqual("0", tn.ParametersToString[2, 1]);
            Assert.AreEqual("1", tn.ParametersToString[3, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new TruncatedDistribution(new Normal(10, 3), 8, 25);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(0.5, tn.Mean, 1e-4);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(0.5, tn.Median, 1e-4);
        }

        /// <summary>
        /// Testing mode
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(0.5, tn.Mode, 1e-4);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(0.19091, tn.StandardDeviation, 1e-4);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(0, tn.Skewness, 1E-4);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(2.62422, tn.Kurtosis, 1e-04);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(0, tn.Minimum);
            Assert.AreEqual(1, tn.Maximum);
        }

        /// <summary>
        /// Testing that a clone is independent of the original and reproduces it exactly.
        /// Reference quantiles from scipy.stats.truncnorm.ppf.
        /// </summary>
        [TestMethod()]
        public void Test_Clone_IsIndependent()
        {
            var tn = new TruncatedDistribution(new Normal(0, 1), -1, 2);
            double pdf = tn.PDF(0.3), cdf = tn.CDF(0.3), median = tn.InverseCDF(0.5);
            Assert.AreEqual(0.17116391801782482, median, 1E-10);

            // The clone reproduces the original bit for bit.
            var clone = (TruncatedDistribution)tn.Clone();
            Assert.AreNotSame(tn.BaseDistribution, clone.BaseDistribution);
            Assert.AreEqual(pdf, clone.PDF(0.3), 0.0);
            Assert.AreEqual(cdf, clone.CDF(0.3), 0.0);
            Assert.AreEqual(median, clone.InverseCDF(0.5), 0.0);
            Assert.AreEqual(tn.Min, clone.Min, 0.0);
            Assert.AreEqual(tn.Max, clone.Max, 0.0);

            // Evaluating another parameter set does not change the original.
            var q = tn.ParameterSetQuantiles(new double[,] { { 5, 2, -1, 10 } }, new[] { 0.5 }, out bool[] valid);
            Assert.IsTrue(valid[0]);
            Assert.AreEqual(4.987818294572313, q[0, 0], 1E-10);
            Assert.AreEqual(median, tn.InverseCDF(0.5), 0.0);
            Assert.AreEqual(0.0, ((Normal)tn.BaseDistribution).Mu, 0.0);
            Assert.AreEqual(1.0, ((Normal)tn.BaseDistribution).Sigma, 0.0);

            // Changing the parameters of the clone's base distribution does not change the original.
            clone.BaseDistribution.SetParameters(new[] { 5d, 2d });
            Assert.AreEqual(0.0, ((Normal)tn.BaseDistribution).Mu, 0.0);
            Assert.AreEqual(pdf, tn.PDF(0.3), 0.0);
            Assert.AreEqual(cdf, tn.CDF(0.3), 0.0);
            Assert.AreEqual(median, tn.InverseCDF(0.5), 0.0);
        }

        /// <summary>
        /// When the base CDF at the lower bound is at most one half, the density, distribution, survival and quantile
        /// functions equal the CDF-difference expressions bit for bit.
        /// </summary>
        [TestMethod()]
        public void Test_LowerHalfMatchesDirectExpressions()
        {
            var cases = new (UnivariateDistributionBase Base, double Min, double Max)[]
            {
                (new Normal(0, 1), -1, 2),
                (new Normal(2, 1), 1.10, 2.11),
                (new Normal(10, 3), 8, 25),
                (new Normal(0, 3), 0, 9),
                (new Normal(0, 1), -10, -7),
                (new Exponential(1), 0, 5),
                (new Gumbel(100, 10), 80, 160),
            };
            foreach (var (baseDistribution, min, max) in cases)
            {
                var tn = new TruncatedDistribution(baseDistribution, min, max);
                double fMin = baseDistribution.CDF(min), fMax = baseDistribution.CDF(max);
                Assert.IsLessThanOrEqualTo(0.5, fMin);
                for (int i = 1; i < 20; i++)
                {
                    double x = min + (max - min) * i / 20d;
                    Assert.AreEqual(baseDistribution.PDF(x) / (fMax - fMin), tn.PDF(x), 0d);
                    double cdf = (baseDistribution.CDF(x) - fMin) / (fMax - fMin);
                    Assert.AreEqual(cdf, tn.CDF(x), 0d);
                    Assert.AreEqual(1d - cdf, tn.CCDF(x), 0d);
                    double p = i / 20d;
                    Assert.AreEqual(baseDistribution.InverseCDF(p * (fMax - fMin) + fMin), tn.InverseCDF(p), 0d);
                }
            }
        }

        /// <summary>
        /// Changing the base distribution through <see cref="TruncatedDistribution.BaseDistribution"/> refreshes the cached
        /// truncation probabilities, moments and validity on the next evaluation.
        /// </summary>
        [TestMethod()]
        public void Test_BaseDistributionChangeRefreshesCache()
        {
            var tn = new TruncatedDistribution(new Normal(0, 1), -1, 1);
            double pdf = tn.PDF(0), mean = tn.Mean;
            var normal = (Normal)tn.BaseDistribution;

            normal.Mu = 5;
            // scipy.stats.truncnorm(-6, -4, loc=5).pdf(0), .cdf(0.9), .ppf(0.5), .mean() and .std()
            Assert.AreEqual(0.04694371747909183, tn.PDF(0), 1E-12 * 0.047);
            Assert.AreEqual(0.6522372543193549, tn.CDF(0.9), 1E-12);
            Assert.AreEqual(0.8389028343201401, tn.InverseCDF(0.5), 1E-12);
            Assert.AreEqual(0.7744530681937984, tn.Mean, 1E-9 * 0.774);
            Assert.AreEqual(0.2157710773954302, tn.StandardDeviation, 1E-9 * 0.216);
            Assert.IsTrue(tn.ParametersValid);

            // The interval [-1, 1] has no mass under Normal(100, 1).
            normal.Mu = 100;
            Assert.Throws<ArgumentOutOfRangeException>(() => tn.PDF(0));
            Assert.IsFalse(tn.ParametersValid);

            normal.Mu = 0;
            Assert.AreEqual(pdf, tn.PDF(0), 0d);
            Assert.AreEqual(mean, tn.Mean, 0d);
            Assert.IsTrue(tn.ParametersValid);
        }

        /// <summary>
        /// An interval in the upper tail of the base distribution is evaluated from survival differences.
        /// </summary>
        [TestMethod()]
        public void Test_UpperTailUsesSurvivalDifferences()
        {
            var tn = new TruncatedDistribution(new Normal(0, 1), 7, 10);
            // mpmath (50 digits): (ncdf(-7) - ncdf(-x)) / (ncdf(-7) - ncdf(-10)) and npdf(x) / (ncdf(-7) - ncdf(-10));
            // scipy.stats.truncnorm(7, 10).ppf(0.5) and .ppf(0.99).
            Assert.AreEqual(0.76476023809564087, tn.CDF(7.2), 1E-13);
            Assert.AreEqual(0.23523976190435913, tn.CCDF(7.2), 1E-13 * 0.235);
            Assert.AreEqual(8.8177929871187948E-8, tn.CCDF(9), 1E-11 * 8.8E-8);
            Assert.AreEqual(1.7252448209690691, tn.PDF(7.2), 1E-13 * 1.73);
            Assert.AreEqual(7.09647265564641, tn.InverseCDF(0.5), 1E-13 * 7.1);
            Assert.AreEqual(7.618841509072155, tn.InverseCDF(0.99), 1E-13 * 7.6);
            Assert.AreEqual(0.5, tn.CDF(tn.InverseCDF(0.5)), 1E-13);
            Assert.AreEqual(1d, tn.CCDF(7), 0d);
            Assert.AreEqual(0d, tn.CCDF(10), 0d);

            // The mass of [8, 10] is 6.2E-16 (mpmath: ncdf(-8) - ncdf(-10)), below the 1E-15 zero-mass threshold.
            Assert.IsFalse(new TruncatedDistribution(new Normal(0, 1), 8, 10).ParametersValid);
        }

        /// <summary>
        /// The moments are integrated adaptively over the support.
        /// </summary>
        [TestMethod()]
        public void Test_MomentsAgainstOracle()
        {
            // mpmath quadrature (50 digits) of the truncated standard normal on [-1, 2];
            // scipy.stats.truncnorm(-1, 2).stats('mvsk') agrees to within 1e-15 (scipy kurtosis is excess).
            var tn = new TruncatedDistribution(new Normal(0, 1), -1, 2);
            Assert.AreEqual(0.22963717909132897, tn.Mean, 1E-9 * 0.23);
            Assert.AreEqual(0.72094558685904579, tn.StandardDeviation, 1E-9 * 0.72);
            Assert.AreEqual(0.32254591283309152, tn.Skewness, 1E-9 * 0.32);
            Assert.AreEqual(2.2801724579228966, tn.Kurtosis, 1E-9 * 2.28);

            // mpmath quadrature of exp(-x) / (1 - exp(-5)) on [0, 5]; scipy.stats.truncexpon(5).stats('mvsk').
            tn = new TruncatedDistribution(new Exponential(1), 0, 5);
            Assert.AreEqual(0.96608172546847884, tn.Mean, 1E-9 * 0.97);
            Assert.AreEqual(0.91063613918798467, tn.StandardDeviation, 1E-9 * 0.91);
            Assert.AreEqual(1.5026215309903858, tn.Skewness, 1E-9 * 1.5);
            Assert.AreEqual(5.2634982040870563, tn.Kurtosis, 1E-9 * 5.26);

            // mpmath quadrature for Normal(10, 3) on [8, 25].
            tn = new TruncatedDistribution(new Normal(10, 3), 8, 25);
            Assert.AreEqual(11.28204729992183, tn.Mean, 1E-9 * 11.3);
            Assert.AreEqual(2.1890999718718675, tn.StandardDeviation, 1E-9 * 2.19);
            Assert.AreEqual(0.73067614427346896, tn.Skewness, 1E-9 * 0.73);
            Assert.AreEqual(3.2319259596961974, tn.Kurtosis, 1E-9 * 3.23);
        }

        /// <summary>
        /// A mode at a truncation bound is returned exactly; an interior mode is found by the density search.
        /// </summary>
        [TestMethod()]
        public void Test_ModeAtBounds()
        {
            // The truncated Exponential density exp(-x) is largest at the lower bound.
            Assert.AreEqual(0d, new TruncatedDistribution(new Exponential(1), 0, 5).Mode, 0d);
            // The Normal(0, 1) density is decreasing on [1, 3] and increasing on [-3, -1].
            Assert.AreEqual(1d, new TruncatedDistribution(new Normal(0, 1), 1, 3).Mode, 0d);
            Assert.AreEqual(-1d, new TruncatedDistribution(new Normal(0, 1), -3, -1).Mode, 0d);
            // Interior mode at the Normal mean.
            Assert.AreEqual(0.5, new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1).Mode, 1E-6);
            Assert.AreEqual(10, new TruncatedDistribution(new Normal(10, 3), 8, 25).Mode, 1E-6);
        }

        /// <summary>
        /// A unimodal base with a closed-form mode inside its support gives the truncated mode directly: the base mode
        /// when the truncation interval holds it, otherwise the nearer bound. The expected values are the textbook
        /// modes: (κ − 1)θ = 4 for Gamma(θ = 2, κ = 3), ξ = 5 for Gumbel, μ − σγ/2 = 9.5 for Pearson III(10, 2, 0.5),
        /// 3 = k − 2 for chi-squared with k = 5, and σ = 2 for Rayleigh. A base whose mode is a support bound, and a
        /// U-shaped Beta, take the global search.
        /// </summary>
        [TestMethod]
        public void Test_ModeOfUnimodalBase()
        {
            Assert.AreEqual(4d, new TruncatedDistribution(new GammaDistribution(2, 3), 0, 20).Mode, 0d);
            Assert.AreEqual(5d, new TruncatedDistribution(new GammaDistribution(2, 3), 5, 20).Mode, 0d);
            Assert.AreEqual(2d, new TruncatedDistribution(new GammaDistribution(2, 3), 1, 2).Mode, 0d);
            Assert.AreEqual(5d, new TruncatedDistribution(new Gumbel(5, 2), 0, 20).Mode, 0d);
            Assert.AreEqual(9.5, new TruncatedDistribution(new PearsonTypeIII(10, 2, 0.5), 0, 30).Mode, 0d);
            Assert.AreEqual(3d, new TruncatedDistribution(new ChiSquared(5), 0, 20).Mode, 0d);
            Assert.AreEqual(2d, new TruncatedDistribution(new Rayleigh(2), 0, 10).Mode, 0d);
            // Gamma shape 1/2: the density is decreasing from its support bound 0, the mode of the base, and the search
            // returns that bound.
            Assert.AreEqual(0d, new TruncatedDistribution(new GammaDistribution(1, 0.5), 0, 5).Mode, 0d);
            // Beta(0.5, 0.5) is U-shaped; truncated to [0.2, 0.9] its density is largest at the upper bound.
            Assert.AreEqual(0.9, new TruncatedDistribution(new BetaDistribution(0.5, 0.5), 0.2, 0.9).Mode, 0d);
        }

        /// <summary>
        /// Mode is global: an equal-weight mixture of N(-2, 0.6) and N(2, 1.2) truncated to [-3.8, 5.6] has its taller peak
        /// near -1.998. The oracle maximizes the mixture density with scipy.optimize.minimize_scalar(method='bounded',
        /// xatol=1e-12), SciPy 1.18.1: -1.9980575185329832.
        /// </summary>
        [TestMethod]
        public void Test_ModeIsGlobalOnBimodalBase()
        {
            var mixture = new Mixture(new[] { 0.5, 0.5 }, new IUnivariateDistribution[] { new Normal(-2d, 0.6), new Normal(2d, 1.2) });
            var truncated = new TruncatedDistribution(mixture, -3.8, 5.6);
            Assert.AreEqual(-1.9980575185329832, truncated.Mode, 1E-8);
        }

        /// <summary>
        /// A spike narrower than a grid cell is found through the base distribution's mode: 0.9·N(0, 10) + 0.1·N(5, 0.01)
        /// truncated to [-30, 30] has its global mode at the spike, 5 (the narrow component's peak, 0.1·φ(0)/0.01 ≈ 3.99,
        /// is far above the wide one's 0.036).
        /// </summary>
        [TestMethod]
        public void Test_ModeFindsSpikeThroughBaseMode()
        {
            var mixture = new Mixture(new[] { 0.9, 0.1 }, new IUnivariateDistribution[] { new Normal(0d, 10d), new Normal(5d, 0.01) });
            var truncated = new TruncatedDistribution(mixture, -30d, 30d);
            Assert.AreEqual(5d, truncated.Mode, 1E-6);
        }
    }
}
