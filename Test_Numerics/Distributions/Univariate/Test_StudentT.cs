using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Student T distribution algorithm.
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
    public class Test_StudentT
    {
        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_StudentT_PDF()
        {
            var t = new StudentT(4d);
            double pdf = t.PDF(1.4d);
            double result = 0.138377537135553d;
            Assert.AreEqual(result, pdf, 1E-10);
            t = new StudentT(2.5d, 0.5d, 4d);
            pdf = t.PDF(1.4d);
            result = 0.1032953042520084d;
            Assert.AreEqual(result, pdf, 1E-10);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_StudentT_CDF()
        {
            var t = new StudentT(4d);
            double cdf = t.CDF(1.4d);
            double result = 0.882949686336585d;
            Assert.AreEqual(result, cdf, 1E-10);
            t = new StudentT(2.5d, 0.5d, 4d);
            cdf = t.CDF(1.4d);
            result = 0.0463263350898173d;
            Assert.AreEqual(result, cdf, 1E-10);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_StudentT_InverseCDF()
        {
            var t = new StudentT(4d);
            double cdf = t.CDF(1.4d);
            double invcdf = t.InverseCDF(cdf);
            double result = 1.4d;
            Assert.AreEqual(result, invcdf, 1E-2);
            t = new StudentT(2.5d, 0.5d, 4d);
            cdf = t.CDF(1.4d);
            invcdf = t.InverseCDF(cdf);
            result = 1.4d;
            Assert.AreEqual(result, invcdf, 1E-2);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var t = new StudentT();
            Assert.AreEqual(0d, t.Mu);
            Assert.AreEqual(1d, t.Sigma);
            Assert.AreEqual(10d, t.DegreesOfFreedom);

            var t2 = new StudentT(10, 10, 10);
            Assert.AreEqual(10d, t2.Mu);
            Assert.AreEqual(10d, t2.Sigma);
            Assert.AreEqual(10d, t2.DegreesOfFreedom);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var t = new StudentT(double.NaN, double.NaN, 1);
            Assert.IsFalse(t.ParametersValid);

            var t2 = new StudentT(double.PositiveInfinity, double.PositiveInfinity, 1);
            Assert.IsFalse(t2.ParametersValid);

            var t3 = new StudentT(1, 1, 0);
            Assert.IsFalse(t3.ParametersValid);
        }
        
        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var t = new StudentT();
            Assert.AreEqual("Location (µ)", t.ParametersToString[0, 0]);
            Assert.AreEqual("Scale (σ)", t.ParametersToString[1, 0]);
            Assert.AreEqual("Degrees of Freedom (ν)", t.ParametersToString[2, 0]);
            Assert.AreEqual("0", t.ParametersToString[0, 1]);
            Assert.AreEqual("1", t.ParametersToString[1,1]);
            Assert.AreEqual("10", t.ParametersToString[2, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new StudentT(10, 1, 100);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var t = new StudentT();
            Assert.AreEqual(0, t.Mean);

            var t2 = new StudentT(1, 1, 1);
            Assert.AreEqual(double.NaN, t2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var t = new StudentT();
            Assert.AreEqual(0, t.Median);

            var t2 = new StudentT(1, 1, 1);
            Assert.AreEqual(1, t2.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var t = new StudentT();
            Assert.AreEqual(0, t.Mode);

            var t2 = new StudentT(1,1,1);
            Assert.AreEqual(1, t2.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var t = new StudentT();
            Assert.AreEqual(1.11803, t.StandardDeviation, 1e-04);

            var t2 = new StudentT(1, 1, 2);
            Assert.AreEqual(double.PositiveInfinity, t2.StandardDeviation);

            var t3 = new StudentT(1, 1, 1);
            Assert.AreEqual(double.NaN, t3.StandardDeviation);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var t = new StudentT();
            Assert.AreEqual(0, t.Skewness);

            var t2 = new StudentT(1, 1, 1);
            Assert.AreEqual(double.NaN, t2.Skewness);
        }

        /// <summary>
        /// Testing kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var t = new StudentT();
            Assert.AreEqual(4, t.Kurtosis);

            var t2 = new StudentT(1, 1, 4);
            Assert.AreEqual(double.PositiveInfinity, t2.Kurtosis);

            var t3 = new StudentT(1, 1, 2);
            Assert.AreEqual(double.NaN, t3.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var t = new StudentT();
            Assert.AreEqual(double.NegativeInfinity, t.Minimum);
            Assert.AreEqual(double.PositiveInfinity, t.Maximum);
        }

        /// <summary>
        /// Verifies that location-scale Student-t densities include the scale Jacobian.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_PDF_LocationScaleIdentityAndNormalization()
        {
            var standard = new StudentT(0.0d, 1.0d, 4.0d);
            var scaled = new StudentT(2.5d, 0.5d, 4.0d);
            const double x = 1.4d;
            double z = (x - scaled.Mu) / scaled.Sigma;

            Assert.AreEqual(standard.PDF(z) / scaled.Sigma, scaled.PDF(x), 1E-14);

            const int intervals = 100000;
            double lower = scaled.Mu - 100.0d * scaled.Sigma;
            double upper = scaled.Mu + 100.0d * scaled.Sigma;
            double width = (upper - lower) / intervals;
            double integral = 0.5d * (scaled.PDF(lower) + scaled.PDF(upper));
            for (int i = 1; i < intervals; i++)
            {
                integral += scaled.PDF(lower + i * width);
            }

            Assert.AreEqual(1.0d, integral * width, 1E-6);
        }

        /// <summary>
        /// Verifies that extreme finite Student-t quantiles retain location and scale.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_InverseCDF_ExtremeTailRetainsLocationAndScale()
        {
            const double probability = 1E-155;
            var distribution = new StudentT(2.0d, 3.0d, 1.0d);

            double actual = distribution.InverseCDF(probability);
            double expected = 2.0d - 3.0d / (Math.PI * probability);

            Assert.IsFalse(double.IsInfinity(actual));
            Assert.AreEqual(expected, actual, Math.Abs(expected) * 1E-12);
            Assert.IsLessThan(distribution.InverseCDF(1E-150), actual);
        }
    
        /// <summary>
        /// The cached ν-dependent density factors reproduce the direct expression bit for bit,
        /// whether ν is set through the constructor, the property or SetParameters.
        /// </summary>
        [TestMethod]
        public void Test_PDF_CachedConstantsMatchDirectExpression()
        {
            var byProperty = new StudentT(0.5d, 2.0d, 3.0d);
            foreach (double v in new[] { 1d, 1.5d, 3d, 10d, 37.25d, 1E4, 9.9E7 })
            {
                byProperty.DegreesOfFreedom = v;
                var bySetter = new StudentT();
                bySetter.SetParameters(0.5d, 2.0d, v);
                var byConstructor = new StudentT(0.5d, 2.0d, v);
                for (double x = -40d; x <= 40d; x += 0.37d)
                {
                    double inverseSigma = 1.0d / 2.0d;
                    double z = (x - 0.5d) * inverseSigma;
                    double expected = Math.Exp(Gamma.LogGamma(0.5d) - Beta.LogFunction(v / 2.0d, 0.5d)) * Math.Pow(1.0d + z * z / v, -0.5d * (v + 1.0d)) * inverseSigma / Math.Sqrt(v * Math.PI);
                    Assert.AreEqual(expected, byProperty.PDF(x), 0d);
                    Assert.AreEqual(expected, bySetter.PDF(x), 0d);
                    Assert.AreEqual(expected, byConstructor.PDF(x), 0d);
                }
            }
        }

        /// <summary>
        /// Infinite and NaN arguments: CDF(+∞) = 1, CDF(−∞) = 0, NaN propagates, and a NaN probability is rejected
        /// by InverseCDF like any other probability outside [0, 1].
        /// </summary>
        [TestMethod]
        public void Test_StudentT_NonFiniteArguments()
        {
            foreach (double v in new[] { 1d, 5d, 1E9 })
            {
                var t = new StudentT(0.5d, 2d, v);
                Assert.AreEqual(1d, t.CDF(double.PositiveInfinity), 0d);
                Assert.AreEqual(0d, t.CDF(double.NegativeInfinity), 0d);
                Assert.AreEqual(0d, t.CCDF(double.PositiveInfinity), 0d);
                Assert.AreEqual(1d, t.CCDF(double.NegativeInfinity), 0d);
                Assert.AreEqual(double.NegativeInfinity, t.LogCDF(double.NegativeInfinity));
                Assert.AreEqual(0d, t.LogCDF(double.PositiveInfinity), 0d);
                Assert.AreEqual(double.NegativeInfinity, t.LogPDF(double.PositiveInfinity));
                Assert.IsTrue(double.IsNaN(t.CDF(double.NaN)));
                Assert.IsTrue(double.IsNaN(t.CCDF(double.NaN)));
                Assert.IsTrue(double.IsNaN(t.LogCDF(double.NaN)));
                Assert.IsTrue(double.IsNaN(t.LogCCDF(double.NaN)));
                Assert.Throws<ArgumentOutOfRangeException>(() => t.InverseCDF(double.NaN));
            }
        }

        /// <summary>
        /// The lower tail is evaluated without cancellation, and the upper tail directly by symmetry.
        /// Reference values are the regularized incomplete beta form I_{ν/(ν+z²)}(ν/2, 1/2)/2 evaluated with
        /// mpmath at 50 digits; scipy.stats.t.cdf agrees where it is finite.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_Tails_MatchOracle()
        {
            // { location, scale, ν, x, CDF } with scipy.stats.t.cdf((x - location)/scale, ν)
            var table = new double[,]
            {
                { 0, 1, 5, -1E5, 9.4901672353943245E-25 },
                { 0, 1, 5, -1E10, 9.4901672455623608E-50 },
                { 0, 1, 5, -3, 0.015049623948731287 },
                { 0, 1, 5, 0.3, 0.61187547886836277 },
                { 0, 1, 1, -1E300, 3.1830988618379065E-301 },
                { 0, 1, 1, 1E10, 0.99999999996816901 },
                { 0, 1, 2.5, -10, 0.0022207478836537124 },
                { 2.5, 0.5, 4.3, -12.5, 1.8438781758720978E-6 },
                { 0, 1, 1000, -1E-8, 0.49999999601157443 },
                // Large ν near the center: scipy.stats.t.cdf(-0.3, 1E6) and scipy.stats.t.cdf(-1E-8, 3E7).
                { 0, 1, 1E6, -0.3, 0.38208860898949715 },
                { 0, 1, 3E7, -1E-8, 0.49999999601057723 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                var t = new StudentT(table[i, 0], table[i, 1], table[i, 2]);
                double expected = table[i, 4];
                // The exp(ln F) path beyond |z| = 1E150 √ν carries the absolute error of a log near -692.
                // At ν = 1E6 the incomplete beta function with parameter ν/2 limits the accuracy to about 1E-10.
                double tolerance = table[i, 3] == -1E300 ? 1E-12 : table[i, 2] == 1E6 ? 1E-10 : 1E-13;
                Assert.AreEqual(expected, t.CDF(table[i, 3]), tolerance * expected, $"CDF row {i}");
                // Reflection: CCDF at the mirrored point equals the lower tail.
                double mirrored = 2d * table[i, 0] - table[i, 3];
                Assert.AreEqual(expected, t.CCDF(mirrored), tolerance * expected, $"CCDF row {i}");
            }
        }

        /// <summary>
        /// Log tail probabilities stay finite after the probability underflows. Reference values from mpmath at
        /// 50 digits of ln(I_{ν/(ν+z²)}(ν/2, 1/2)/2); scipy.stats.t.logcdf(-1E10, 5) = -112.87899841388881 and
        /// scipy returns -inf at -1E200.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_LogTails_MatchOracle()
        {
            var t = new StudentT(5d);
            Assert.AreEqual(-112.87899841388881, t.LogCDF(-1E10), 1E-14 * 112.87899841388881);
            Assert.AreEqual(-2300.3348367582322, t.LogCDF(-1E200), 1E-14 * 2300.3348367582322);
            Assert.AreEqual(-2300.3348367582322, t.LogCCDF(1E200), 1E-14 * 2300.3348367582322);
            Assert.AreEqual(-55.314371090109094, t.LogCCDF(1E5), 1E-14 * 55.314371090109094);
            Assert.AreEqual(-4.1964022748946908, t.LogCDF(-3d), 1E-14 * 4.1964022748946908);
            // ln(1 − F) for the upper body: scipy.stats.t.logcdf(1E10, 1) = -3.1830988618885674E-11.
            Assert.AreEqual(-3.1830988618885673E-11, new StudentT(1d).LogCDF(1E10), 1E-13 * 3.1830988618885673E-11);
            Assert.AreEqual(-691.92025778406311, new StudentT(1d).LogCDF(-1E300), 1E-14 * 691.92025778406311);
        }

        /// <summary>
        /// The log density is evaluated in log space. Reference values are
        /// ln Γ((ν+1)/2) − ln Γ(ν/2) − ln √(νπ) − ((ν+1)/2) ln(1 + z²/ν) from mpmath at 50 digits;
        /// scipy.stats.t.logpdf(40, 3.7) = -15.254323253964012 and scipy.stats.t.logpdf(1E30, 1) = -139.29983546549212.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_LogPDF_MatchesOracle()
        {
            Assert.AreEqual(-2759.2424174446072, new StudentT(5d).LogPDF(-1E200), 1E-14 * 2759.2424174446072);
            Assert.AreEqual(-15.254323253964013, new StudentT(3.7d).LogPDF(40d), 1E-14 * 15.254323253964013);
            Assert.AreEqual(-139.29983546549214, new StudentT(1d).LogPDF(1E30), 1E-14 * 139.29983546549214);
            // Location-scale Jacobian and agreement with ln PDF in the body.
            var t = new StudentT(2.5d, 0.5d, 4.3d);
            for (double x = -5d; x <= 10d; x += 0.75d)
            {
                Assert.AreEqual(Math.Log(t.PDF(x)), t.LogPDF(x), 1E-13 * Math.Max(1d, Math.Abs(Math.Log(t.PDF(x)))));
            }
            Assert.AreEqual(Math.Log(new StudentT(1E9).PDF(1.3d)), new StudentT(1E9).LogPDF(1.3d), 1E-15);
        }

        /// <summary>
        /// At large ν the quantile is the Student's t quantile, not the Normal one, and it is consistent with the CDF.
        /// For ν = 10^9 scipy.stats.t.ppf(0.975, 1E9) = 1.959963986912325, which differs from the Normal quantile
        /// 1.959963984540054 by the O(1/ν) term (z³ + z)/(4ν) ≈ 2.4E-9.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_LargeDegreesOfFreedom_QuantileMatchesCDF()
        {
            var t = new StudentT(1E9);
            Assert.AreEqual(1.959963986912325, t.InverseCDF(0.975), 1E-14);
            var scaled = new StudentT(3d, 2d, 1E8);
            foreach (double p in new[] { 1E-10, 0.025, 0.4, 0.5, 0.975, 1 - 1E-10 })
            {
                double x = scaled.InverseCDF(p);
                Assert.AreEqual(p, scaled.CDF(x), 1E-14 * Math.Max(p, 1E-3));
            }
        }

        /// <summary>
        /// For 1/4 &lt; p &lt; 3/4 (and ν below 1E16), and in the tails below ν = 30, the quantile reproduces the incomplete
        /// beta inversion bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_InverseCDF_CentreAndSmallDegreesOfFreedomMatchIncompleteBetaInversion()
        {
            foreach (double v in new[] { 1d, 4.3d, 29.9d, 30d, 999d, 1E4, 1E7 })
            {
                var t = new StudentT(1.5d, 2d, v);
                foreach (double p in new[] { 1E-12, 0.1d, 0.3d, 0.5d, 0.7d, 0.9d, 1 - 1E-12 })
                {
                    bool centre = p > 0.25d && p < 0.75d;
                    if (!centre && v >= 30d) continue;
                    double expected;
                    if (centre)
                    {
                        if (p == 0.5d) expected = 1.5d;
                        else
                        {
                            double z = Beta.IncompleteInverse(0.5d, 0.5d * v, Math.Abs(1.0d - 2.0d * p));
                            double q = Math.Sqrt(v * z / (1.0d - z));
                            expected = 1.5d + 2d * (p < 0.5d ? -q : q);
                        }
                    }
                    else
                    {
                        double probability = p;
                        int rflg = -1;
                        if (probability >= 0.5d) { probability = 1.0d - probability; rflg = 1; }
                        double z = Beta.IncompleteInverse(0.5d * v, 0.5d, 2.0d * probability);
                        double q = Math.Sqrt(v) * Math.Sqrt(1.0d - z) / Math.Sqrt(z);
                        expected = 1.5d + 2d * (rflg * q);
                    }
                    Assert.AreEqual(expected, t.InverseCDF(p), 0d, $"v={v}, p={p}");
                }
            }
        }

        /// <summary>
        /// Far-tail quantiles for small ν. For ν = 1 below 1E-150 the quantile is the Cauchy closed form -1/(πq) (before,
        /// 1E-300 returned -∞ and 1E-160 was 7.6E-6 off); a quantile beyond the double range returns -∞. For 1 &lt; ν &lt; 30,
        /// where 2 min(p, 1 - p) is subnormal, the log CDF is solved (before, p = 5E-324 was 41% off at ν = 2 and 1% at
        /// ν = 29). Reference values are mpmath roots of F(t) = p at 45 digits; the largest error measured on these points
        /// is 3.7E-14 relative, at ν = 2 and p = 1E-310.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_InverseCDF_SmallDegreesOfFreedomFarTail_MpmathOracle()
        {
            var cases = new (double Nu, double P, double T)[]
            {
                (1d, 1E-300, -3.1830988618379066E+299),
                (1d, 1E-160, -3.1830988618379068E+159),
                (2d, 5E-324, -3.1812124520951962E+161),
                (2d, 1E-310, -7.071067811865486E+154),
                (5d, 5E-324, -7.1894859915199676E+64),
                (5d, 1E-310, -1.5683925590993388E+62),
                (10d, 5E-324, -5.4907110967913065E+32),
                (29d, 5E-324, -692762909739.36694),
                (29.9, 5E-324, -325541834931.01195),
            };
            foreach (var (nu, p, expected) in cases)
            {
                Assert.AreEqual(expected, new StudentT(nu).InverseCDF(p), 1E-13 * Math.Abs(expected), $"nu = {nu}, p = {p}");
            }
            Assert.AreEqual(double.NegativeInfinity, new StudentT(1d).InverseCDF(5E-324));
            // The upper tail by symmetry: mpmath cot(π 2^-53) = 2.8670805696113293E+15 at ν = 1, p = 1 - 2^-53 (1 - p exact).
            Assert.AreEqual(2.8670805696113293E+15, new StudentT(1d).InverseCDF(1d - 1.1102230246251565E-16), 1E-13 * 2.8670805696113293E+15);
        }

        /// <summary>
        /// The quantile keeps its relative accuracy near the median for large ν and in the far tails for ν ≥ 30. Reference
        /// values are mpmath roots at 45 digits of F(t) = p, with F(t) = I_y(ν/2, 1/2)/2, y = ν/(ν + t²), for t &lt; 0, and
        /// for ν = 1E305 the Cornish-Fisher series z + (z³ + z)/(4ν) + (5z⁵ + 16z³ + 3z)/(96ν²) about the Normal quantile z.
        /// The largest error measured on these points is 5.2E-15 relative, at ν = 31 and p = 5E-324.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_InverseCDF_CentreAndFarTail_MpmathOracle()
        {
            var cases = new (double Nu, double P, double T)[]
            {
                (1000d, 0.4999999999, -2.5072552173850608E-10),
                (1000d, 0.49, -0.025075180209466442),
                (1E4, 0.49999999, -2.5066909398016785E-08),
                (1E6, 0.499999, -2.5066289012237064E-06),
                (1E6, 0.49, -0.025068914529877549),
                (1E7, 0.4999, -0.00025066283635457872),
                (1E305, 0.4999999999, -2.5066284820303539E-10),
                (1E305, 0.26, -0.64334540539291694),
                (30d, 1E-12, -11.397217523311411),
                (30d, 1E-20, -22.658878371940183),
                (150d, 1E-10, -6.8281782136557930),
                (300d, 1E-10, -6.5878115764844004),
                (300d, 1E-20, -9.9750513703355011),
                (500d, 1E-20, -9.6794176375543989),
                (999d, 1E-20, -9.4672535155058546),
                (31d, 5E-324, -137363880149.62192),
                (1E305, 1E-300, -37.047096299361199),
                (300d, 0.9999999999, 6.5878115624166830),
            };
            foreach (var (nu, p, expected) in cases)
            {
                Assert.AreEqual(expected, new StudentT(nu).InverseCDF(p), 1E-14 * Math.Abs(expected), $"nu = {nu}, p = {p}");
            }
        }

        /// <summary>
        /// For ν above 2^997, where the double-double residual of z²/ν needs scaling, the CDF, CCDF and log CDF are finite and
        /// equal the Normal values (scipy.stats.norm: cdf(1) = 0.8413447460685429, cdf(-3) = 0.001349898031630093,
        /// sf(5) = 2.8665157187919333E-7, logcdf(-40) = -804.6084420137539).
        /// </summary>
        [TestMethod]
        public void Test_StudentT_HugeDegreesOfFreedom_MatchesNormal()
        {
            foreach (double nu in new[] { 1E305, 1.7E308 })
            {
                var t = new StudentT(nu);
                Assert.AreEqual(0.8413447460685429, t.CDF(1d), 1E-14, $"nu = {nu}");
                Assert.AreEqual(0.001349898031630093, t.CDF(-3d), 1E-14 * 0.001349898031630093, $"nu = {nu}");
                Assert.AreEqual(2.8665157187919333E-7, t.CCDF(5d), 1E-14 * 2.8665157187919333E-7, $"nu = {nu}");
                Assert.AreEqual(-804.6084420137539, t.LogCDF(-40d), 1E-14 * 804.6084420137539, $"nu = {nu}");
            }
        }

        /// <summary>
        /// Maximum likelihood for the location, scale and degrees of freedom attains at least the log-likelihood
        /// of scipy.stats.t.fit on the same samples.
        /// </summary>
        /// <remarks>
        /// Sample 1: scipy.stats.t.fit gives df = 1.1626019650189297, loc = 0.9910552845161906,
        /// scale = 0.6137294370282971, log-likelihood -37.65269862393729; a tight scipy Nelder-Mead refinement
        /// reaches -37.65269861397859 at (0.9910739847712444, 0.6137139112335337, 1.1626138184828332).
        /// Sample 2: scipy.stats.t.fit stops at a local optimum with log-likelihood -33.872165523590866; the global
        /// optimum is -28.61633037997114 at (10.24056126084135, 0.5564721050323235, 1.1581478079288554).
        /// </remarks>
        [TestMethod]
        public void Test_StudentT_MLE_MatchesScipyFit()
        {
            var sample1 = new[] { 2.1, -0.4, 1.3, 0.8, 5.9, 1.1, 0.2, 1.7, -2.6, 0.9, 1.4, 0.6, 3.2, 1.0, -0.1, 1.2, 0.7, 8.4, 1.5, 0.3 };
            var t1 = new StudentT();
            t1.Estimate(sample1, ParameterEstimationMethod.MaximumLikelihood);
            Assert.IsTrue(t1.ParametersValid);
            Assert.IsGreaterThanOrEqualTo(-37.65269862393729 - 1E-6, t1.LogLikelihood(sample1));
            Assert.AreEqual(0.9910739847712444, t1.Mu, 1E-3);
            Assert.AreEqual(0.6137139112335337, t1.Sigma, 1E-3);
            Assert.AreEqual(1.1626138184828332, t1.DegreesOfFreedom, 1E-3);

            var sample2 = new[] { 10.2, 11.5, 9.1, 10.8, 14.9, 10.4, 9.7, 10.1, 3.8, 10.6, 11.0, 9.9, 12.3, 10.3, 9.5, 10.0 };
            var t2 = new StudentT();
            t2.Estimate(sample2, ParameterEstimationMethod.MaximumLikelihood);
            Assert.IsGreaterThanOrEqualTo(-33.872165523590866 - 1E-6, t2.LogLikelihood(sample2));
        }

        /// <summary>
        /// Estimation input validation, the unsupported method of moments, and bootstrapping.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_Estimation_ValidationAndBootstrap()
        {
            var t = new StudentT(1d, 2d, 6d);
            Assert.Throws<NotImplementedException>(() => t.Estimate(new[] { 1d, 2d, 3d, 5d, 8d }, ParameterEstimationMethod.MethodOfMoments));
            Assert.Throws<ArgumentOutOfRangeException>(() => t.Estimate(new[] { 1d, 2d, 3d }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.Throws<ArgumentOutOfRangeException>(() => t.Estimate(new[] { 1d, double.NaN, 3d, 4d }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.Throws<ArgumentOutOfRangeException>(() => t.GetParameterConstraints(new[] { 2d, 2d, 2d, 2d }));
            var constraints = t.GetParameterConstraints(new[] { 2.1, -0.4, 1.3, 0.8, 5.9, 1.1, 0.2, 1.7 });
            for (int i = 0; i < 3; i++)
            {
                Assert.IsTrue(constraints.Item2[i] <= constraints.Item1[i] && constraints.Item1[i] <= constraints.Item3[i]);
            }
            var bootstrapped = (StudentT)t.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 200, 12345);
            Assert.IsTrue(bootstrapped.ParametersValid);
            Assert.AreEqual(1d, bootstrapped.Mu, 0.5d);
        }

        /// <summary>
        /// The two-sided tail 2 CCDF(|z|) = I_x(ν/2, 1/2), x = ν/(ν + z²), matches mpmath (the incomplete beta function at
        /// 120 digits) on a 97-point grid with ν from 1 to 1E10 and |z| from 1E-6 to 1E6, leaving out tails below the double
        /// range: within 3.4E-15 relative where the tail is at least 0.01, 1.1E-14 down to 1E-10, and 8E-13 below. The
        /// largest errors measured on these points are 6.5E-16, 3.4E-15 and 2.9E-13.
        /// </summary>
        [TestMethod]
        public void Test_StudentT_TwoSidedTail_MpmathGrid()
        {
            var cases = new (double nu, double z, double tail)[]
            {
                (1.0, 1e-06, 0.99999936338022766),
                (1.0, 0.01, 0.99363401447018351),
                (1.0, 0.3, 0.81445284184451527),
                (1.0, 1.0, 0.50000000000000000),
                (1.0, 2.0, 0.29516723530086653),
                (1.0, 4.0, 0.15595826075473865),
                (1.0, 8.0, 0.079166848321131078),
                (1.0, 16.0, 0.039737048611081675),
                (1.0, 40.0, 0.015912179824051628),
                (1.0, 1000.0, 0.00063661956016111785),
                (1.0, 1000000.0, 6.3661977236736911E-7),
                (3.0, 1e-06, 0.99999926489480606),
                (3.0, 0.01, 0.99264911141285195),
                (3.0, 0.3, 0.78376329203991901),
                (3.0, 1.0, 0.39100221895577064),
                (3.0, 2.0, 0.13932596855884319),
                (3.0, 4.0, 0.028008456010146166),
                (3.0, 8.0, 0.0040765775877854684),
                (3.0, 16.0, 0.00053092969702083740),
                (3.0, 40.0, 3.4380680789158526E-5),
                (3.0, 1000.0, 2.2053076425765927E-9),
                (3.0, 1000000.0, 2.2053155816792292E-18),
                (10.0, 1e-06, 0.99999922178323208),
                (10.0, 0.01, 0.99221797499097142),
                (10.0, 0.3, 0.77032060756579857),
                (10.0, 1.0, 0.34089313230205986),
                (10.0, 2.0, 0.073388034770740365),
                (10.0, 4.0, 0.0025183326247366924),
                (10.0, 8.0, 1.1774942789666160E-5),
                (10.0, 16.0, 1.8776253538036420E-8),
                (10.0, 40.0, 2.2808577430857547E-12),
                (10.0, 1000.0, 2.4608247101732710E-26),
                (10.0, 1000000.0, 2.4609374998872070E-56),
                (30.0, 1e-06, 0.99999920873563020),
                (30.0, 0.01, 0.99208749257317774),
                (30.0, 0.3, 0.76624610528435277),
                (30.0, 1.0, 0.32530861542602990),
                (30.0, 2.0, 0.054625044962983101),
                (30.0, 4.0, 0.00038184563608375686),
                (30.0, 8.0, 6.2658224757007590E-9),
                (30.0, 16.0, 3.1158405836886140E-16),
                (30.0, 40.0, 1.3726045194406403E-27),
                (30.0, 1000.0, 2.0720034831117330E-69),
                (30.0, 1000000.0, 2.0729069296087555E-159),
                (70.0, 1e-06, 0.99999920495984784),
                (70.0, 0.01, 0.99204973287615850),
                (70.0, 0.3, 0.76506623745549374),
                (70.0, 1.0, 0.32075485777547313),
                (70.0, 2.0, 0.049381610801980447),
                (70.0, 4.0, 0.00015513361399488923),
                (70.0, 8.0, 1.8268184528593240E-11),
                (70.0, 16.0, 4.4107232586832800E-25),
                (70.0, 40.0, 5.8915516533823502E-50),
                (70.0, 1000.0, 3.5910564845991770E-147),
                (300.0, 1e-06, 0.99999920278006482),
                (300.0, 0.01, 0.99202793395888178),
                (300.0, 0.3, 0.76438491882814130),
                (300.0, 1.0, 0.31811640443062628),
                (300.0, 2.0, 0.046401520983149092),
                (300.0, 4.0, 7.9840044837194648E-5),
                (300.0, 8.0, 2.7350153901538423E-14),
                (300.0, 16.0, 4.3324653679334651E-42),
                (300.0, 40.0, 2.8524180395311874E-122),
                (1000.0, 1e-06, 0.99999920231488537),
                (1000.0, 0.01, 0.99202328193218725),
                (1000.0, 0.3, 0.76423950416724407),
                (1000.0, 1.0, 0.31755241808467233),
                (1000.0, 2.0, 0.045770346493251637),
                (1000.0, 4.0, 6.8009919208781580E-5),
                (1000.0, 8.0, 3.4266614823914745E-15),
                (1000.0, 16.0, 1.7811877613675644E-51),
                (1000.0, 40.0, 1.0478852155173361E-209),
                (3000.0, 1e-06, 0.99999920218192684),
                (3000.0, 0.01, 0.99202195228006107),
                (3000.0, 0.3, 0.76419794032952204),
                (3000.0, 1.0, 0.31739115804930473),
                (3000.0, 2.0, 0.045590263070947962),
                (3000.0, 4.0, 6.4872169589387304E-5),
                (3000.0, 8.0, 1.7603665036110302E-15),
                (3000.0, 16.0, 2.3397750301300284E-55),
                (3000.0, 40.0, 8.6601595494364699E-281),
                (100000.0, 1e-06, 0.99999920211743387),
                (100000.0, 0.01, 0.99202130731882321),
                (100000.0, 0.3, 0.76417777919033769),
                (100000.0, 1.0, 0.31731292756411000),
                (100000.0, 2.0, 0.045502963457506464),
                (100000.0, 4.0, 6.3387997557195344E-5),
                (100000.0, 8.0, 1.2573919876256373E-15),
                (100000.0, 16.0, 1.5067247617721307E-57),
                (10000000000.0, 1e-06, 0.99999920211543925),
                (10000000000.0, 0.01, 0.99202128737093631),
                (10000000000.0, 0.3, 0.76417715562833044),
                (10000000000.0, 1.0, 0.31731050788711118),
                (10000000000.0, 2.0, 0.045500263923353899),
                (10000000000.0, 4.0, 6.3342484121262614E-5),
                (10000000000.0, 8.0, 1.2441922462134114E-15),
                (10000000000.0, 16.0, 1.2777529899006010E-57),
            };
            foreach (var (nu, z, expected) in cases)
            {
                double bound = expected >= 0.01 ? 3.4E-15 : expected >= 1E-10 ? 1.1E-14 : 8E-13;
                Assert.AreEqual(expected, 2d * new StudentT(nu).CCDF(z), bound * expected, $"nu = {nu}, z = {z}");
            }
        }

        /// <summary>
        /// The upper tail at large ν, including ν ≥ 10^8 where a Normal approximation is 3.8E-3 off at z = 35, matches
        /// mpmath (the incomplete beta continued fraction at 50 digits) to the documented accuracy, and the log tail
        /// stays finite.
        /// </summary>
        [TestMethod]
        public void Test_CCDF_LargeDegreesOfFreedom_MpmathOracle()
        {
            var cases = new (double nu, double z, double ccdf)[]
            {
                (1E8, 35d, 1.1291456838795503E-268),
                (1E9, 5d, 2.8665162019758142E-7),
                (1E4, 0.5d, 0.30854303966161671),
                (1E7, 1E-6, 0.49999960105772957),
                (1E12, 30d, 4.9067149229646632E-198),
                (1E30, 10d, 7.6198530241605261E-24),
                (5E3, 40d, 4.2074252359681858E-304),
                (1E6, 8d, 6.2275317166012598E-16),
            };
            foreach (var (nu, z, expected) in cases)
            {
                var t = new StudentT(nu);
                Assert.AreEqual(expected, t.CCDF(z), expected * 1E-12, $"nu = {nu}, z = {z}");
                Assert.AreEqual(expected, t.CDF(-z), expected * 1E-12, $"nu = {nu}, z = {z}");
                Assert.AreEqual(System.Math.Log(expected), t.LogCCDF(z), 1E-12 * System.Math.Abs(System.Math.Log(expected)), $"nu = {nu}, z = {z}");
            }
        }

        /// <summary>
        /// The quantile at ν = 10^7 and p = 1E-300 returns within a time limit (the inverse incomplete beta function did
        /// not terminate there) and inverts the CDF; the quantile solves the log CDF in the tails for ν ≥ 30.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_LargeDegreesOfFreedom_FarTail()
        {
            foreach (var (nu, p) in new[] { (1E7, 1E-300), (1E5, 1E-300), (5E7, 1E-200), (1E9, 1E-10), (1E7, 0.3), (1E7, 1d - 1E-12) })
            {
                var t = new StudentT(nu);
                var task = System.Threading.Tasks.Task.Factory.StartNew(() => t.InverseCDF(p), System.Threading.Tasks.TaskCreationOptions.LongRunning);
                Assert.IsTrue(task.Wait(10000), $"InverseCDF did not return for nu = {nu}, p = {p}");
                double x = task.Result;
                double q = p < 0.5 ? t.CDF(x) : t.CCDF(x);
                double target = p < 0.5 ? p : 1d - p;
                Assert.AreEqual(target, q, target * 1E-11, $"nu = {nu}, p = {p}");
            }
        }
    }
}
