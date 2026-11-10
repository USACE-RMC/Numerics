using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Truncated Normal distribution algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil</item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests" />
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_TruncatedNormal
    {

        /// <summary>
        /// This method was verified against the R package "truncnorm"
        /// </summary>
        [TestMethod]
        public void Test_TruncatedNormalDist()
        {
            var tn = new TruncatedNormal(2, 1, 1.10, 2.11);
            var d = tn.PDF(1.5);
            var p = tn.CDF(1.5);
            var q = tn.InverseCDF(p);

            Assert.AreEqual(0.9786791, d, 1E-5);
            Assert.AreEqual(0.3460251, p, 1E-5);
            Assert.AreEqual(1.5, q, 1E-5);

            tn = new TruncatedNormal(10, 3, 8, 25);
            d = tn.PDF(12.75);
            p = tn.CDF(12.75);
            q = tn.InverseCDF(p);

            Assert.AreEqual(0.1168717, d, 1E-5);
            Assert.AreEqual(0.7596566, p, 1E-5);
            Assert.AreEqual(12.75, q, 1E-5);

            tn = new TruncatedNormal(0, 3, 0, 9);
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
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.5, tn.Mu);
            Assert.AreEqual(0.2, tn.Sigma);
            Assert.AreEqual(0, tn.Min);
            Assert.AreEqual(1, tn.Max);

            var tn2 = new TruncatedNormal(1, 1, 1, 2);
            Assert.AreEqual(1, tn2.Mu);
            Assert.AreEqual(1, tn2.Sigma);
            Assert.AreEqual(1, tn2.Min);
            Assert.AreEqual(2, tn2.Max);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var tn = new TruncatedNormal(double.NaN, double.NaN, double.NaN, double.NaN);
            Assert.IsFalse(tn.ParametersValid);

            var tn2 = new TruncatedNormal(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
            Assert.IsFalse(tn2.ParametersValid);

            var tn3 = new TruncatedNormal(0, -1, -1, 0);
            Assert.IsFalse(tn3.ParametersValid);

            var tn4 = new TruncatedNormal(1, 1, 1, 0);
            Assert.IsFalse(tn4.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var tn = new TruncatedNormal();
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
            var dist = new TruncatedNormal(10, 3, 8, 25);
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
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.5, tn.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.5, tn.Median);
        }

        /// <summary>
        /// Testing mode
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.5, tn.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.19091, tn.StandardDeviation, 1e-05);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0, tn.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(2.62422, tn.Kurtosis, 1e-04);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0, tn.Minimum);
            Assert.AreEqual(1, tn.Maximum);
        }

        /// <summary>
        /// The cached standardized bounds reproduce the direct per-call expressions bit for bit,
        /// including after each parameter setter is used on its own.
        /// </summary>
        [TestMethod()]
        public void Test_CachedBoundsMatchDirectExpressions()
        {
            var cases = new[]
            {
                new[] { 2d, 1d, 1.10, 2.11 },
                new[] { 10d, 3d, 8d, 25d },
                new[] { 0d, 3d, 0d, 9d },
                new[] { 0.5, 0.2, 0d, 1d },
                new[] { 4d, 0.7, -3d, 6d },
            };
            foreach (var c in cases)
            {
                var tn = new TruncatedNormal();
                tn.Mu = c[0];
                tn.Sigma = c[1];
                tn.Min = c[2];
                tn.Max = c[3];
                AssertMatchesDirect(tn);
                AssertMatchesDirect(new TruncatedNormal(c[0], c[1], c[2], c[3]));
            }
        }

        /// <summary>
        /// Intervals that contain the parent mean (α ≤ 0 ≤ β) are evaluated with the Φ(β) − Φ(α) expressions,
        /// bit for bit, over a grid of locations, scales and one- or two-sided truncations.
        /// </summary>
        [TestMethod()]
        public void Test_InteriorTruncationsMatchDirectExpressions()
        {
            var multiples = new[] { 0d, 0.1, 1d, 3d, 9d };
            foreach (double mu in new[] { -2d, -0.5, 0d, 0.3, 1.7 })
            {
                foreach (double sigma in new[] { 0.3, 1d, 4d })
                {
                    foreach (double below in multiples)
                    {
                        foreach (double above in multiples)
                        {
                            if (below == 0 && above == 0) continue;
                            AssertMatchesDirect(new TruncatedNormal(mu, sigma, mu - below * sigma, mu + above * sigma));
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Both bounds above the parent mean: TN(0, 1, 8, 9) cancels in Φ(β) − Φ(α).
        /// </summary>
        [TestMethod()]
        public void Test_UpperTailInterval()
        {
            var tn = new TruncatedNormal(0, 1, 8, 9);
            // mpmath (80 digits): npdf(8.1) / (ncdf(-8) - ncdf(-9)) and the corresponding CDF and bisection quantiles.
            // scipy.stats.truncnorm(8, 9) agrees to within 1e-14.
            Assert.AreEqual(3.631624468775537, tn.PDF(8.1), 1E-12 * 3.631624468775537);
            Assert.AreEqual(0.5583754014201233, tn.CDF(8.1), 1E-12);
            Assert.AreEqual(0.9849406286168291, tn.CDF(8.5), 1E-12);
            Assert.AreEqual(8.084888899018166, tn.InverseCDF(0.5), 1E-12 * 8.08);
            Assert.AreEqual(8.001237199059592, tn.InverseCDF(0.01), 1E-12 * 8.0);
            Assert.AreEqual(8.546703681127418, tn.InverseCDF(0.99), 1E-12 * 8.5);
            Assert.AreEqual(8.084888899018166, tn.Median, 1E-12 * 8.08);
            // mpmath quadrature of the moments over [8, 9]. The closed-form central moments lose
            // digits to cancellation in far tails; scipy.stats.truncnorm(8, 9).stats('mvsk') shows the same loss.
            Assert.AreEqual(8.121188992979796, tn.Mean, 1E-13 * 8.12);
            Assert.AreEqual(0.11894764723502567, tn.StandardDeviation, 1E-9 * 0.119);
            Assert.AreEqual(1.8532542803795513, tn.Skewness, 1E-7 * 1.85);
            Assert.AreEqual(7.678063491684261, tn.Kurtosis, 1E-6 * 7.68);
            Assert.AreEqual(8, tn.Mode);
        }

        /// <summary>
        /// Both bounds far above the parent mean: TN(0, 1, 40, 41), where Φ(β) − Φ(α) is exactly zero in double precision.
        /// </summary>
        [TestMethod()]
        public void Test_FarUpperTailInterval()
        {
            var tn = new TruncatedNormal(0, 1, 40, 41);
            // mpmath (80 digits) as in Test_UpperTailInterval; scipy.stats.truncnorm(40, 41) agrees to within 1e-13.
            Assert.AreEqual(26.82819751682549, tn.PDF(40.01), 1E-11 * 26.8);
            Assert.AreEqual(7.280388487857421e-08, tn.PDF(40.5), 1E-10 * 7.28e-08);
            Assert.AreEqual(0.32988079019628447, tn.CDF(40.01), 1E-11);
            Assert.AreEqual(0.9999999982034672, tn.CDF(40.5), 1E-13);
            Assert.AreEqual(40.01731412676465, tn.InverseCDF(0.5), 1E-13 * 40);
            Assert.AreEqual(40.0000000000025, tn.InverseCDF(1E-10), 1E-13 * 40);
            Assert.AreEqual(40.17221588247908, tn.InverseCDF(0.999), 1E-13 * 40);
            Assert.AreEqual(40.02496884720726, tn.Mean, 1E-13 * 40);
            // The variance is a difference of terms near α² = 1600, so about six digits remain.
            Assert.AreEqual(0.02495332399884605, tn.StandardDeviation, 1E-5 * 0.025);
        }

        /// <summary>
        /// Both bounds far below the parent mean mirror the far upper tail: TN(0, 1, -41, -40).
        /// </summary>
        [TestMethod()]
        public void Test_FarLowerTailInterval()
        {
            var tn = new TruncatedNormal(0, 1, -41, -40);
            // mpmath (80 digits): npdf(x) / (ncdf(-40) - ncdf(-41)) and the corresponding CDF and quantiles.
            Assert.AreEqual(26.82819751682549, tn.PDF(-40.01), 1E-11 * 26.8);
            Assert.AreEqual(0.6701192098037155, tn.CDF(-40.01), 1E-11);
            Assert.AreEqual(1.7965328361726676e-09, tn.CDF(-40.5), 1E-10 * 1.8e-09);
            Assert.AreEqual(-40.01731412676465, tn.InverseCDF(0.5), 1E-13 * 40);
            Assert.AreEqual(-40.00002499689694, tn.InverseCDF(0.999), 1E-13 * 40);
            Assert.AreEqual(-40.02496884720726, tn.Mean, 1E-13 * 40);
            Assert.AreEqual(0.02495332399884605, tn.StandardDeviation, 1E-5 * 0.025);
            Assert.AreEqual(-40, tn.Mode);
        }

        /// <summary>
        /// Moderate one-tail intervals above and below the parent mean, against mpmath.
        /// </summary>
        [TestMethod()]
        public void Test_OneTailIntervals()
        {
            // mpmath (80 digits) for TN(-4, 0.7, -3, 6): α = 1.43, β = 14.3; scipy.stats.truncnorm agrees to within 1e-14.
            var tn = new TruncatedNormal(-4, 0.7, -3, 6);
            Assert.AreEqual(2.1655500739175264, tn.PDF(-2.9), 1E-13 * 2.17);
            Assert.AreEqual(6.043197873903228e-07, tn.PDF(0.0), 1E-12 * 6.04e-07);
            Assert.AreEqual(0.2419181997371335, tn.CDF(-2.9), 1E-13);
            Assert.AreEqual(0.7902102399423461, tn.CDF(-2.5), 1E-13);
            Assert.AreEqual(-2.760312728335569, tn.InverseCDF(0.5), 1E-13 * 2.76);
            Assert.AreEqual(-2.3025723413949675, tn.InverseCDF(0.9), 1E-13 * 2.3);
            Assert.AreEqual(-2.68529772742997, tn.Mean, 1E-13 * 2.69);
            Assert.AreEqual(0.27615250690375404, tn.StandardDeviation, 1E-12 * 0.276);
            Assert.AreEqual(1.4219454592837737, tn.Skewness, 1E-11 * 1.42);
            Assert.AreEqual(5.460721171387187, tn.Kurtosis, 1E-11 * 5.46);

            // mpmath (80 digits) for TN(5, 2, 1, 3): α = -2, β = -1.
            tn = new TruncatedNormal(5, 2, 1, 3);
            Assert.AreEqual(0.31741746583004915, tn.PDF(1.5), 1E-13 * 0.317);
            Assert.AreEqual(0.12736109326157433, tn.CDF(1.5), 1E-13);
            Assert.AreEqual(0.9132027006530883, tn.CDF(2.9), 1E-13);
            Assert.AreEqual(2.327119304509506, tn.InverseCDF(0.5), 1E-13 * 2.33);
            Assert.AreEqual(1.8342751518795535, tn.InverseCDF(0.25), 1E-13 * 1.83);
            Assert.AreEqual(2.2336619067368946, tn.Mean, 1E-13 * 2.23);
            Assert.AreEqual(0.5394177828014248, tn.StandardDeviation, 1E-12 * 0.539);
            Assert.AreEqual(-0.49096097265072586, tn.Skewness, 1E-11);
            Assert.AreEqual(2.1748993666482543, tn.Kurtosis, 1E-11 * 2.17);

            // mpmath (80 digits) for TN(0, 1, 2, 3).
            tn = new TruncatedNormal(0, 1, 2, 3);
            Assert.AreEqual(0.8190705093187625, tn.PDF(2.5), 1E-13 * 0.819);
            Assert.AreEqual(0.7729105526090576, tn.CDF(2.5), 1E-13);
            Assert.AreEqual(2.2555308204204385, tn.InverseCDF(0.5), 1E-13 * 2.26);
            Assert.AreEqual(2.315821326743782, tn.Mean, 1E-13 * 2.32);
            Assert.AreEqual(0.2480338274802249, tn.StandardDeviation, 1E-12 * 0.248);
            Assert.AreEqual(0.8015928244415931, tn.Skewness, 1E-10 * 0.80);
            Assert.AreEqual(2.7575863323156984, tn.Kurtosis, 1E-10 * 2.76);
        }

        /// <summary>
        /// Min = Max is a point mass at Min with the conventions of <see cref="Deterministic"/>.
        /// </summary>
        [TestMethod()]
        public void Test_PointMass()
        {
            foreach (double mu in new[] { -1d, 2.5, 7d })
            {
                var tn = new TruncatedNormal(mu, 1.5, 2.5, 2.5);
                var reference = new Deterministic(2.5);
                Assert.IsTrue(tn.ParametersValid);
                foreach (double x in new[] { 2.0, 2.5, 3.0 })
                {
                    Assert.AreEqual(reference.PDF(x), tn.PDF(x));
                    Assert.AreEqual(reference.CDF(x), tn.CDF(x));
                }
                foreach (double p in new[] { 0d, 0.3, 1d })
                    Assert.AreEqual(2.5, tn.InverseCDF(p));
                Assert.AreEqual(2.5, tn.Mean);
                Assert.AreEqual(2.5, tn.Median);
                Assert.AreEqual(2.5, tn.Mode);
                Assert.AreEqual(0d, tn.StandardDeviation);
                Assert.IsTrue(double.IsNaN(tn.Skewness));
                Assert.IsTrue(double.IsNaN(tn.Kurtosis));
                Assert.AreEqual(2.5, tn.Minimum);
                Assert.AreEqual(2.5, tn.Maximum);
            }
            Assert.IsFalse(new TruncatedNormal(0, 1, double.NegativeInfinity, double.NegativeInfinity).ParametersValid);
            Assert.IsFalse(new TruncatedNormal(0, 1, double.PositiveInfinity, double.PositiveInfinity).ParametersValid);
            Assert.IsFalse(new TruncatedNormal(0, 1, 0, double.PositiveInfinity).ParametersValid);
        }

        /// <summary>
        /// The method of moments keeps the sample extremes as bounds and reproduces the sample mean and standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_MethodOfMoments()
        {
            // mpmath.findroot on the closed-form truncated mean and standard deviation (50 digits), started from
            // scipy.optimize.least_squares on scipy.stats.truncnorm(...).mean() and .std(); both agree to 1e-15.
            var sample1 = new[] { 5.0, 5.1, 4.9, 5.2, 4.8, 5.05, 4.95, 6.5, 3.6, 5.0 };
            var sample2 = new[] { 2.962, 3.442, 3.709, 3.907, 4.07, 4.21, 4.337, 4.454, 4.564, 4.669, 4.77, 4.869,
                4.967, 5.064, 5.162, 5.262, 5.364, 5.47, 5.582, 5.701, 5.832, 5.979, 6.15, 6.364 };
            var cases = new[]
            {
                (sample1, 4.986691586233058, 0.8713918092253526),
                (sample2, 5.189149491017554, 1.4158498058612876),
            };
            foreach (var (sample, mu, sigma) in cases)
            {
                var tn = new TruncatedNormal();
                tn.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
                var mom = Numerics.Data.Statistics.Statistics.ProductMoments(sample);
                Assert.IsTrue(tn.ParametersValid);
                Assert.AreEqual(sample.Min(), tn.Min);
                Assert.AreEqual(sample.Max(), tn.Max);
                Assert.AreEqual(mom[0], tn.Mean, 1E-8 * Math.Abs(mom[0]));
                Assert.AreEqual(mom[1], tn.StandardDeviation, 1E-8 * mom[1]);
                Assert.AreEqual(mu, tn.Mu, 1E-8 * Math.Abs(mu));
                Assert.AreEqual(sigma, tn.Sigma, 1E-8 * sigma);
            }
        }

        /// <summary>
        /// A sample whose standard deviation exceeds what any truncated normal on its range can reach falls back to the
        /// least-squares fit at the largest scale, and a constant sample yields the point mass.
        /// </summary>
        [TestMethod()]
        public void Test_MethodOfMoments_InfeasibleAndConstantSamples()
        {
            // Sample SD 0.6746 exceeds the uniform SD (2.2 - 0.05) / sqrt(12) = 0.6207 on the sample range.
            var sample = new[] { 0.2, 0.5, 0.1, 1.4, 0.3, 0.9, 0.05, 0.7, 2.2, 0.4 };
            var tn = new TruncatedNormal();
            tn.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            Assert.IsTrue(tn.ParametersValid);
            Assert.AreEqual(0.05, tn.Min);
            Assert.AreEqual(2.2, tn.Max);
            Assert.AreEqual(5d * (2.2 - 0.05), tn.Sigma, 1E-12 * tn.Sigma);
            // mpmath (50 digits): golden-section minimum over µ of (mean - 0.675)^2 + (sd - 0.6746398216003032)^2
            // for the truncated normal with σ = 10.75 on [0.05, 2.2], using the closed-form truncated moments. The
            // objective is flat in µ (it changes by 1E-12 over 3E-4), so µ is resolved only to about 1E-5 relative.
            // The truncated mean and standard deviation change by about 2.5E-3 per unit of µ here, so that resolution
            // leaves them uncertain by about 3.3E-6.
            Assert.AreEqual(-128.99079286692825, tn.Mu, 1E-5 * 129);
            Assert.AreEqual(0.72893145379445987, tn.Mean, 5E-6);
            Assert.AreEqual(0.54209246706701137, tn.StandardDeviation, 5E-6);
            Assert.IsLessThanOrEqualTo((2.2 - 0.05) / Math.Sqrt(12d), tn.StandardDeviation);

            var constant = new TruncatedNormal();
            constant.Estimate(new[] { 3.25, 3.25, 3.25, 3.25, 3.25 }, ParameterEstimationMethod.MethodOfMoments);
            Assert.IsTrue(constant.ParametersValid);
            Assert.AreEqual(3.25, constant.Min);
            Assert.AreEqual(3.25, constant.Max);
            Assert.AreEqual(3.25, constant.Mean);
            Assert.AreEqual(0d, constant.StandardDeviation);
        }

        /// <summary>
        /// Invalid samples and unsupported methods are rejected.
        /// </summary>
        [TestMethod()]
        public void Test_MethodOfMoments_InvalidSamples()
        {
            var tn = new TruncatedNormal();
            Assert.Throws<ArgumentNullException>(() => tn.Estimate(null!, ParameterEstimationMethod.MethodOfMoments));
            Assert.Throws<ArgumentOutOfRangeException>(() => tn.Estimate(new[] { 1d, 2d, 3d }, ParameterEstimationMethod.MethodOfMoments));
            Assert.Throws<ArgumentOutOfRangeException>(() => tn.Estimate(new[] { 1d, 2d, double.NaN, 3d }, ParameterEstimationMethod.MethodOfMoments));
            Assert.Throws<ArgumentOutOfRangeException>(() => tn.Estimate(new[] { double.NaN, double.NaN, double.NaN, double.NaN }, ParameterEstimationMethod.MethodOfMoments));
            Assert.Throws<ArgumentOutOfRangeException>(() => tn.Estimate(new[] { 1d, 2d, double.PositiveInfinity, 3d }, ParameterEstimationMethod.MethodOfMoments));
            Assert.Throws<NotImplementedException>(() => tn.Estimate(new[] { 1d, 2d, 3d, 4d }, ParameterEstimationMethod.MaximumLikelihood));
        }

        /// <summary>
        /// The parametric bootstrap runs on the truncated normal with the method of moments.
        /// </summary>
        [TestMethod()]
        public void Test_Bootstrap()
        {
            var tn = new TruncatedNormal(5, 1, 2, 6.5);
            Assert.IsInstanceOfType(tn, typeof(IBootstrappable));
            var replicate = (TruncatedNormal)tn.Bootstrap(ParameterEstimationMethod.MethodOfMoments, 200, 12345);
            Assert.IsTrue(replicate.ParametersValid);
            Assert.IsTrue(replicate.Min >= 2 && replicate.Max <= 6.5);

            var boot = new BootstrapAnalysis(tn, ParameterEstimationMethod.MethodOfMoments, 100, 200, 12345);
            var distributions = boot.Distributions();
            Assert.HasCount(200, distributions);
            var results = boot.Estimate(new[] { 0.1, 0.5, 0.9 });
            var median = results.ModeCurve![1];
            Assert.IsTrue(median > 2 && median < 6.5);
        }

        /// <summary>
        /// Concurrent evaluations on one instance return the same values as sequential evaluations.
        /// </summary>
        [TestMethod()]
        public void Test_ConcurrentEvaluationsAreIndependent()
        {
            var tn = new TruncatedNormal(0.5, 0.2, 0, 1);
            var xs = new double[2000];
            for (int i = 0; i < xs.Length; i++) xs[i] = (i + 0.5) / xs.Length;
            var expectedPdf = new double[xs.Length];
            var expectedCdf = new double[xs.Length];
            for (int i = 0; i < xs.Length; i++)
            {
                expectedPdf[i] = tn.PDF(xs[i]);
                expectedCdf[i] = tn.CDF(xs[i]);
            }
            for (int repeat = 0; repeat < 20; repeat++)
            {
                var pdf = new double[xs.Length];
                var cdf = new double[xs.Length];
                System.Threading.Tasks.Parallel.For(0, xs.Length, i =>
                {
                    pdf[i] = tn.PDF(xs[i]);
                    cdf[i] = tn.CDF(xs[i]);
                    _ = tn.Mean;
                    _ = tn.StandardDeviation;
                });
                CollectionAssert.AreEqual(expectedPdf, pdf);
                CollectionAssert.AreEqual(expectedCdf, cdf);
            }
        }

        /// <summary>
        /// Asserts exact equality with the expressions evaluated directly from the parameters.
        /// </summary>
        private static void AssertMatchesDirect(TruncatedNormal tn)
        {
            double mu = tn.Mu, sigma = tn.Sigma, min = tn.Min, max = tn.Max;
            double alpha = (min - mu) / sigma;
            double beta = (max - mu) / sigma;
            double z = Normal.StandardCDF(beta) - Normal.StandardCDF(alpha);

            Assert.AreEqual(mu + sigma * ((Normal.StandardPDF(alpha) - Normal.StandardPDF(beta)) / z), tn.Mean, 0d);
            Assert.AreEqual(mu + sigma * Normal.StandardZ((Normal.StandardCDF(alpha) + Normal.StandardCDF(beta)) / 2d), tn.Median, 0d);
            Assert.AreEqual(Math.Sqrt(Math.Pow(sigma, 2d) * (1d + (alpha * Normal.StandardPDF(alpha) - beta * Normal.StandardPDF(beta)) / z - Math.Pow((Normal.StandardPDF(alpha) - Normal.StandardPDF(beta)) / z, 2d))), tn.StandardDeviation, 0d);
            double z0 = Normal.StandardPDF(alpha) / z, z1 = Normal.StandardPDF(beta) / z;
            double v = 1d - (beta * z1 - alpha * z0) - Math.Pow(z1 - z0, 2d);
            Assert.AreEqual(-(1d / Math.Pow(v, 3d / 2d)) * (2d * Math.Pow(z1 - z0, 3d) + (3d * beta * z1 - 3d * alpha * z0 - 1d) * (z1 - z0) + Math.Pow(beta, 2d) * z1 - Math.Pow(alpha, 2d) * z0), tn.Skewness, 0d);
            Assert.AreEqual((1d / Math.Pow(v, 2d)) * (-3 * Math.Pow(z1 - z0, 4d) - 6d * (beta * z1 - alpha * z0) * Math.Pow(z1 - z0, 2d) - 2d * Math.Pow(z1 - z0, 2d) - 4d * (Math.Pow(beta, 2d) * z1 - Math.Pow(alpha, 2d) * z0) * (z1 - z0) - 3d * (beta * z1 - alpha * z0) - (Math.Pow(beta, 3d) * z1 - Math.Pow(alpha, 3d) * z0) + 3d), tn.Kurtosis, 0d);
            for (int i = 1; i < 20; i++)
            {
                double x = min + (max - min) * i / 20d;
                double xi = (x - mu) / sigma;
                Assert.AreEqual(Normal.StandardPDF(xi) / (sigma * z), tn.PDF(x), 0d);
                Assert.AreEqual((Normal.StandardCDF(xi) - Normal.StandardCDF(alpha)) / z, tn.CDF(x), 0d);
                double p = i / 20d;
                Assert.AreEqual(mu + sigma * Normal.StandardZ(Normal.StandardCDF(alpha) + p * z), tn.InverseCDF(p), 0d);
            }
        }
    }
}
