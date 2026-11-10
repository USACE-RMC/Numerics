using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;
using System;
using System.Diagnostics;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Noncentral T distribution algorithm.
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
    public class Test_NoncentralT
    {
        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_PDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double pdf = t.PDF(1.4d);
            double result = 0.23552141805184526d;
            Assert.AreEqual(pdf, result, 1E-6);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_CDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double cdf = t.CDF(1.4d);
            double result = 0.15955740661144721d;
            Assert.AreEqual(cdf, result, 1E-6);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_InverseCDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double cdf = t.CDF(1.4d);
            double invcdf = t.InverseCDF(cdf);
            double result = 1.4d;
            Assert.AreEqual(invcdf, result, 1E-6);
            var table = new[,] { { 3.0d, 0.0d, 1d, 0.89758361765043326d }, { 3.0d, 0.0d, 2d, 0.9522670169d }, { 3.0d, 0.0d, 3d, 0.97116555718878128d }, { 3.0d, 0.5d, 1d, 0.8231218864d }, { 3.0d, 0.5d, 2d, 0.904902151d }, { 3.0d, 0.5d, 3d, 0.9363471834d }, { 3.0d, 1.0d, 1d, 0.7301025986d }, { 3.0d, 1.0d, 2d, 0.8335594263d }, { 3.0d, 1.0d, 3d, 0.8774010255d }, { 3.0d, 2.0d, 1d, 0.5248571617d }, { 3.0d, 2.0d, 2d, 0.6293856597d }, { 3.0d, 2.0d, 3d, 0.6800271741d }, { 3.0d, 4.0d, 1d, 0.20590131975d }, { 3.0d, 4.0d, 2d, 0.2112148916d }, { 3.0d, 4.0d, 3d, 0.2074730718d }, { 15.0d, 7.0d, 15d, 0.9981130072d }, { 15.0d, 7.0d, 20d, 0.999487385d }, { 15.0d, 7.0d, 25d, 0.9998391562d }, { 0.05d, 1.0d, 1d, 0.168610566972d }, { 0.05d, 1.0d, 2d, 0.16967950985d }, { 0.05d, 1.0d, 3d, 0.1701041003d }, { 4.0d, 2.0d, 10d, 0.9247683363d }, { 4.0d, 3.0d, 10d, 0.7483139269d }, { 4.0d, 4.0d, 10d, 0.4659802096d }, { 5.0d, 2.0d, 10d, 0.9761872541d }, { 5.0d, 3.0d, 10d, 0.8979689357d }, { 5.0d, 4.0d, 10d, 0.7181904627d }, { 6.0d, 2.0d, 10d, 0.9923658945d }, { 6.0d, 3.0d, 10d, 0.9610341649d }, { 6.0d, 4.0d, 10d, 0.868800735d } };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double x = table[i, 0];
                double delta = table[i, 1];
                double degF = table[i, 2];
                var target = new NoncentralT(degF, delta);
                double expected = table[i, 3];
                double actual = target.CDF(x);
                Assert.AreEqual(expected, actual, 0.000001d);
                double expectedX = target.InverseCDF(actual);
                Assert.AreEqual(expectedX, x, 0.000001d);
            }

        }

        /// <summary>
        /// Verifies that zero noncentrality reduces to the central Student-t distribution across both tails.
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_ZeroNoncentralityMatchesStudentT()
        {
            var noncentral = new NoncentralT(12d, 0d);
            var central = new StudentT(0d, 1d, 12d);
            foreach (double x in new[] { -3d, -1d, 0d, 1d, 3d })
            {
                Assert.AreEqual(central.CDF(x), noncentral.CDF(x), 1E-8, $"CDF mismatch at x={x}.");
                Assert.AreEqual(central.PDF(x), noncentral.PDF(x), 1E-7, $"PDF mismatch at x={x}.");
            }
        }

        /// <summary>
        /// Verifies negative-tail symmetry, the exact value at zero, and inverse-CDF round trips.
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_NegativeTailSymmetryAndRoundTrips()
        {
            const double degreesOfFreedom = 9d;
            const double noncentrality = 1.75d;
            var positive = new NoncentralT(degreesOfFreedom, noncentrality);
            var reflected = new NoncentralT(degreesOfFreedom, -noncentrality);

            foreach (double x in new[] { -4d, -2d, -0.5d, 0d, 1.5d })
            {
                Assert.AreEqual(
                    positive.CDF(x),
                    1d - reflected.CDF(-x),
                    1E-8,
                    $"Reflection identity failed at x={x}.");
            }
            Assert.AreEqual(Normal.StandardCDF(-noncentrality), positive.CDF(0d), 1E-12);

            foreach (double probability in new[] { 0.05d, 0.25d, 0.75d, 0.95d })
            {
                double quantile = positive.InverseCDF(probability);
                Assert.AreEqual(probability, positive.CDF(quantile), 1E-6, $"Round trip failed at p={probability}.");
            }
        }

        /// <summary>
        /// Verifying input parameters can create Noncentral T distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var t = new NoncentralT();
            Assert.AreEqual(10d, t.DegreesOfFreedom);
            Assert.AreEqual(0d, t.Noncentrality);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(1d, t2.DegreesOfFreedom);
            Assert.AreEqual(1d, t2.Noncentrality);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var t = new NoncentralT(0, 1);
            Assert.IsFalse(t.ParametersValid);

            var t2 = new NoncentralT(1,double.PositiveInfinity);
            Assert.IsFalse(t2.ParametersValid);

            var t3 = new NoncentralT(1,double.NaN);
            Assert.IsFalse(t3.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var t = new NoncentralT();
            Assert.AreEqual("Degrees of Freedom (ν)", t.ParametersToString[0, 0]);
            Assert.AreEqual("Noncentrality (μ)", t.ParametersToString[1, 0]);
            Assert.AreEqual("10", t.ParametersToString[0, 1]);
            Assert.AreEqual("0", t.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0, t.Mean);

            var t2 = new NoncentralT(0, 1);
            Assert.AreEqual(double.NaN, t2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0, t.Median, 1e-04);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(1.3202, t2.Median, 1e-04);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0, t.Mode,  1E-4);

            var t3 = new NoncentralT(10, 1);
            Assert.AreEqual(0.9329, t3.Mode, 1e-04);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var t = new NoncentralT();
            Assert.AreEqual(1.1180, t.StandardDeviation,1e-04);

            var t2 = new NoncentralT(1, 0);
            Assert.AreEqual(double.NaN,t2.StandardDeviation);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0.0, t.Skewness,  1E-4);
        }

        /// <summary>
        /// Testing Kurtosis
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var t = new NoncentralT();
            Assert.AreEqual(4.0, t.Kurtosis, 1E-4);
        }

        /// <summary>
        /// Testing the mean and standard deviation against high-precision values of the closed forms
        /// E[T] = μ √(ν/2) Γ((ν − 1)/2) / Γ(ν/2) and Var[T] = ν(1 + μ²)/(ν − 2) − E[T]².
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath 1.4.1 at 50 significant digits. The mean is within 1E-15 relative
        /// and the standard deviation, which subtracts E[T]² from a nearly equal term, within 2E-15.
        /// </remarks>
        [TestMethod()]
        public void Test_Mean_StandardDeviation_MatchHighPrecisionClosedForm()
        {
            // { v, delta, mean, standard deviation }
            var table = new double[,]
            {
                { 2.5, 0, 0, 2.23606797749979 },
                { 2.5, 1.5, 2.267299944151681, 3.3330693006971495 },
                { 2.5, -2, -3.0230665922022415, 3.9825956836102163 },
                { 10, 0, 0, 1.118033988749895 },
                { 10, 1.5, 1.6255834619087155, 1.1916284690996921 },
                { 10, -2, -2.167444615878287, 1.2458667011764235 },
                { 100, 0, 0, 1.0101525445522108 },
                { 100, 1.5, 1.5113684302350077, 1.015919287493457 },
                { 100, -2, -2.01515790698001, 1.020382000165857 },
                { 343, 0, 0, 1.0029282639547537 },
                { 343, 1.5, 1.5032898747352625, 1.0045800794840356 },
                { 343, -2, -2.004386499647017, 1.0058629495260576 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double v = table[i, 0], d = table[i, 1];
                var t = new NoncentralT(v, d);
                Assert.AreEqual(table[i, 2], t.Mean, 1E-15 * Math.Abs(table[i, 2]), $"Mean v={v}, d={d}");
                Assert.AreEqual(table[i, 3], t.StandardDeviation, 2E-15 * table[i, 3], $"SD v={v}, d={d}");
            }
        }

        /// <summary>
        /// Testing the mean and standard deviation for large degrees of freedom, where Γ(ν/2) overflows.
        /// Reference values from scipy.stats.nct.stats.
        /// </summary>
        [TestMethod()]
        public void Test_Mean_StandardDeviation_LargeDegreesOfFreedom()
        {
            // { v, delta, mean, standard deviation }
            var table = new double[,]
            {
                { 343, -2, -2.0043864996470164, 1.0058629495260585 },
                { 343.5, 0.5, 0.5010950238044641, 1.0031073827537582 },
                { 344, 1.5, 1.5032802820970828, 1.0045666933615367 },
                { 345, 0, 0.0, 1.0029112143106624 },
                { 345, 1.5, 1.5032707452368215, 1.0045533852578048 },
                { 345, -2, -2.004360993649095, 1.0058287757309636 },
                { 1000, 0, 0.0, 1.0010015025043828 },
                { 1000, 1.5, 1.5011261731066854, 1.001565394003897 },
                { 1000, -2, -2.001501564142247, 1.002003756887318 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                var t = new NoncentralT(table[i, 0], table[i, 1]);
                Assert.AreEqual(table[i, 2], t.Mean, 1E-12 * Math.Abs(table[i, 2]), $"Mean row {i}");
                Assert.AreEqual(table[i, 3], t.StandardDeviation, 1E-12 * table[i, 3], $"SD row {i}");
            }
        }

        /// <summary>
        /// Testing the skewness and kurtosis against the closed-form moments of the noncentral t distribution,
        /// evaluated with 40-digit arithmetic (mpmath). For nonzero noncentrality these agree with
        /// scipy.stats.nct.stats to about 1E-13 relative.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness_Kurtosis_ClosedForm()
        {
            // { v, delta, skewness, kurtosis, relative tolerance for skewness, relative tolerance for kurtosis }
            var table = new double[,]
            {
                { 5, 0, 0.0, 9.0, 0, 1E-14 },
                { 5, 1.5, 1.7323349900113402, 17.164195033641584, 1E-14, 1E-14 },
                { 5, -2, -2.0763577541042182, 20.843148665973196, 1E-14, 1E-14 },
                { 10, 0, 0.0, 4.0, 0, 1E-14 },
                { 10, 1.5, 0.57405824314791092, 4.5175595001275408, 1E-13, 1E-14 },
                { 10, -2, -0.72363329408286697, 4.8273295684844696, 1E-13, 1E-14 },
                { 100, 0, 0.0, 3.0625, 0, 1E-14 },
                { 100, 1.5, 0.045924745854224846, 3.0653440729924442, 1E-12, 1E-13 },
                { 100, -2, -0.060878015186453743, 3.0675012540385171, 1E-12, 1E-13 },
                { 345, 0, 0.0, 3.0175953079178886, 0, 1E-14 },
                { 345, 1.5, 0.01311968259623633, 3.0178255449638961, 5E-11, 1E-12 },
                { 345, -2, -0.01746338571611568, 3.0180033234821498, 5E-11, 1E-12 },
                { 1000, 0, 0.0, 3.0060240963855422, 0, 1E-14 },
                { 1000, 1.5, 0.0045090240910851196, 3.0060512345450807, 2E-9, 2E-11 },
                { 1000, -2, -0.006008527268037405, 3.0060722892744272, 2E-9, 2E-11 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                var t = new NoncentralT(table[i, 0], table[i, 1]);
                Assert.AreEqual(table[i, 2], t.Skewness, table[i, 4] * Math.Abs(table[i, 2]), $"Skewness row {i}");
                Assert.AreEqual(table[i, 3], t.Kurtosis, table[i, 5] * table[i, 3], $"Kurtosis row {i}");
            }
        }

        /// <summary>
        /// Testing that the skewness is NaN for ν ≤ 3 and the kurtosis is NaN for ν ≤ 4, where the moments do not exist.
        /// Reference values from the closed-form moments evaluated with 40-digit arithmetic (mpmath), which agree with
        /// scipy.stats.nct.stats.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness_Kurtosis_ExistenceThresholds()
        {
            Assert.AreEqual(double.NaN, new NoncentralT(2.5, 1.5).Skewness);
            Assert.AreEqual(double.NaN, new NoncentralT(3, 1.5).Skewness);
            Assert.AreEqual(5.8286740738679441, new NoncentralT(3.5, 1.5).Skewness, 1E-13 * 5.8286740738679441);

            Assert.AreEqual(double.NaN, new NoncentralT(3.5, 1.5).Kurtosis);
            Assert.AreEqual(double.NaN, new NoncentralT(4, 1.5).Kurtosis);
            Assert.AreEqual(35.160976085109604, new NoncentralT(4.5, 1.5).Kurtosis, 1E-13 * 35.160976085109604);
        }

        /// <summary>
        /// Testing min and max functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var t = new NoncentralT();
            Assert.AreEqual(double.NegativeInfinity,t.Minimum);
            Assert.AreEqual(double.PositiveInfinity, t.Maximum);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(double.NegativeInfinity, t2.Minimum);
            Assert.AreEqual(double.PositiveInfinity, t2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0.38910, t.PDF(0), 1e-04);
            Assert.AreEqual(0.23036, t.PDF(1),1e-04);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0.82955, t.CDF(1), 1e-04);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(double.NegativeInfinity,t.InverseCDF(0) );
            Assert.AreEqual(double.PositiveInfinity, t.InverseCDF(1));
            Assert.AreEqual(-0.26018, t.InverseCDF(0.4), 1e-04);
        }

        /// <summary>
        /// Validates NoncentralT InverseCDF against Stedinger (1983) Table 1.
        /// Table 1: Percentage points of ζ(0.90)-distribution for the 10-year event.
        /// Relationship: ζ_α(p) = NCT.InverseCDF(α; df=n-1, δ=z_p·√n) / √n
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_StedsingerTable1()
        {
            double zp = Normal.StandardZ(0.90);

            // Exact values from Stedinger (1983) Table 1 (p=0.90)
            // { sampleSize, alpha, expected zeta }
            var table = new[,]
            {
                { 10d, 0.005, 0.436 }, { 10d, 0.050, 0.712 }, { 10d, 0.250, 1.043 },
                { 10d, 0.750, 1.671 }, { 10d, 0.950, 2.355 }, { 10d, 0.995, 3.368 },
                { 20d, 0.005, 0.651 }, { 20d, 0.050, 0.858 }, { 20d, 0.250, 1.104 },
                { 20d, 0.750, 1.528 }, { 20d, 0.950, 1.926 }, { 20d, 0.995, 2.423 },
                { 50d, 0.005, 0.857 }, { 50d, 0.050, 1.000 }, { 50d, 0.250, 1.164 },
                { 50d, 0.750, 1.426 }, { 50d, 0.950, 1.646 }, { 50d, 0.995, 1.890 },
                { 100d, 0.005, 0.970 }, { 100d, 0.050, 1.077 }, { 100d, 0.250, 1.196 },
                { 100d, 0.750, 1.380 }, { 100d, 0.950, 1.527 }, { 100d, 0.995, 1.683 },
            };

            for (int i = 0; i < table.GetLength(0); i++)
            {
                int n = (int)table[i, 0];
                double alpha = table[i, 1];
                double expectedZeta = table[i, 2];

                int df = n - 1;
                double delta = zp * Math.Sqrt(n);
                var nct = new NoncentralT(df, delta);
                double computedZeta = nct.InverseCDF(alpha) / Math.Sqrt(n);

                Assert.AreEqual(expectedZeta, computedZeta, 0.002,
                    $"Table 1 mismatch: n={n}, α={alpha}, expected={expectedZeta}, computed={computedZeta:F4}");
            }
        }

        /// <summary>
        /// Validates NoncentralT InverseCDF against Stedinger (1983) Table 3.
        /// Table 3: Percentage points of ζ(0.99)-distribution for the 100-year event.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_StedsingerTable3()
        {
            double zp = Normal.StandardZ(0.99);

            // Exact values from Stedinger (1983) Table 3 (p=0.99)
            var table = new[,]
            {
                { 10d, 0.005, 1.232 }, { 10d, 0.050, 1.562 }, { 10d, 0.250, 2.008 },
                { 10d, 0.750, 2.927 }, { 10d, 0.950, 3.981 }, { 10d, 0.995, 5.582 },
                { 20d, 0.005, 1.483 }, { 20d, 0.050, 1.749 }, { 20d, 0.250, 2.085 },
                { 20d, 0.750, 2.697 }, { 20d, 0.950, 3.295 }, { 20d, 0.995, 4.059 },
                { 50d, 0.005, 1.744 }, { 50d, 0.050, 1.936 }, { 50d, 0.250, 2.163 },
                { 50d, 0.750, 2.538 }, { 50d, 0.950, 2.862 }, { 50d, 0.995, 3.230 },
                { 100d, 0.005, 1.894 }, { 100d, 0.050, 2.040 }, { 100d, 0.250, 2.207 },
                { 100d, 0.750, 2.470 }, { 100d, 0.950, 2.684 }, { 100d, 0.995, 2.915 },
            };

            for (int i = 0; i < table.GetLength(0); i++)
            {
                int n = (int)table[i, 0];
                double alpha = table[i, 1];
                double expectedZeta = table[i, 2];

                int df = n - 1;
                double delta = zp * Math.Sqrt(n);
                var nct = new NoncentralT(df, delta);
                double computedZeta = nct.InverseCDF(alpha) / Math.Sqrt(n);

                Assert.AreEqual(expectedZeta, computedZeta, 0.002,
                    $"Table 3 mismatch: n={n}, α={alpha}, expected={expectedZeta}, computed={computedZeta:F4}");
            }
        }

        /// <summary>
        /// Tests NoncentralT CDF and InverseCDF with large noncentrality parameters
        /// encountered in real confidence interval calculations (δ > 20).
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_LargeNoncentrality()
        {
            double[] probs = { 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99 };

            // δ = z_0.99 * √100 ≈ 23.3
            var nct1 = new NoncentralT(99, Normal.StandardZ(0.99) * Math.Sqrt(100));
            foreach (double p in probs)
            {
                double x = nct1.InverseCDF(p);
                Assert.IsFalse(double.IsNaN(x), $"InverseCDF({p}) returned NaN for δ≈23.3");
                Assert.IsFalse(double.IsInfinity(x), $"InverseCDF({p}) returned Infinity for δ≈23.3");
                double roundtrip = nct1.CDF(x);
                Assert.AreEqual(p, roundtrip, 1e-4, $"CDF(InverseCDF({p})) roundtrip failed for δ≈23.3");
            }

            // δ = z_0.99 * √200 ≈ 32.9
            var nct2 = new NoncentralT(199, Normal.StandardZ(0.99) * Math.Sqrt(200));
            foreach (double p in probs)
            {
                double x = nct2.InverseCDF(p);
                Assert.IsFalse(double.IsNaN(x), $"InverseCDF({p}) returned NaN for δ≈32.9");
                Assert.IsFalse(double.IsInfinity(x), $"InverseCDF({p}) returned Infinity for δ≈32.9");
                double roundtrip = nct2.CDF(x);
                Assert.AreEqual(p, roundtrip, 1e-4, $"CDF(InverseCDF({p})) roundtrip failed for δ≈32.9");
            }
        }

        /// <summary>
        /// Verifies that NoncentralTConfidenceIntervals produces results equivalent to
        /// MonteCarloConfidenceIntervals. The noncentral-t method is exact (Stedinger 1983);
        /// the Monte Carlo method is approximate.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_CI_vs_MonteCarlo_CI()
        {
            var dist = new Normal(100, 15);
            int sampleSize = 30;
            var quantiles = new double[] { 0.01, 0.10, 0.50, 0.90, 0.99 };
            var percentiles = new double[] { 0.05, 0.25, 0.50, 0.75, 0.95 };

            var nctCI = dist.NoncentralTConfidenceIntervals(sampleSize, quantiles, percentiles);
            var mcCI = dist.MonteCarloConfidenceIntervals(sampleSize, 100000, quantiles, percentiles);

            for (int i = 0; i < quantiles.Length; i++)
            {
                for (int j = 0; j < percentiles.Length; j++)
                {
                    double nctVal = nctCI[i, j];
                    double mcVal = mcCI[i, j];
                    double tol = Math.Abs(nctVal) > 1.0 ? 0.02 * Math.Abs(nctVal) : 0.5;
                    Assert.AreEqual(nctVal, mcVal, tol,
                        $"NCT vs MC: quantile={quantiles[i]}, percentile={percentiles[j]}: NCT={nctVal:F4}, MC={mcVal:F4}");
                }
            }
        }

        /// <summary>
        /// Relative closeness helper for oracle comparisons.
        /// </summary>
        private static void AssertRelative(double expected, double actual, double tolerance, string message)
        {
            Assert.AreEqual(expected, actual, tolerance * Math.Abs(expected), message);
        }

        /// <summary>
        /// Noncentrality above about 38.6, where e^(−μ²/2) underflows, keeps the mass of the Poisson-weighted series.
        /// Reference values: the Guenther series summed with mpmath at 80 digits (CDF), the confluent hypergeometric
        /// closed form of the density with mpmath (PDF), and its root (quantile); scipy.stats.nct.cdf(40, 5, 40) =
        /// 0.41626, scipy.stats.nct.pdf(40, 5, 40) = 0.030415 and scipy.stats.nct.ppf(0.5, 5, 40) = 42.87261857692019 agree.
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_VeryLargeNoncentrality()
        {
            var t = new NoncentralT(5d, 40d);
            AssertRelative(0.41625955025322483, t.CDF(40d), 1E-12, "CDF");
            AssertRelative(0.58374044974677517, t.CCDF(40d), 1E-12, "CCDF");
            AssertRelative(0.030415436072516136, t.PDF(40d), 1E-12, "PDF");
            AssertRelative(42.872618576920202, t.InverseCDF(0.5d), 1E-12, "Median");
        }

        /// <summary>
        /// The noncentral t confidence intervals of the Normal distribution for a large sample, whose noncentralities
        /// (40.5 and 73.6) exceed the range of e^(−μ²/2). Reference values from
        /// 10 + scipy.stats.nct.ppf(c, 999, scipy.stats.norm.ppf(q) * sqrt(1000)) * 2 / sqrt(1000).
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_ConfidenceIntervals_LargeSample()
        {
            var ci = new Normal(10d, 2d).NoncentralTConfidenceIntervals(1000, new[] { 0.9d, 0.99d }, new[] { 0.05d, 0.5d, 0.95d });
            var expected = new double[,]
            {
                { 12.42647633461422, 12.563841213942627, 12.707634942450525 },
                { 14.459090820301334, 14.65414405301809, 14.860280306483387 },
            };
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 3; j++)
                    AssertRelative(expected[i, j], ci[i, j], 1E-11, $"CI [{i},{j}]");
        }

        /// <summary>
        /// Infinite and NaN arguments.
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_NonFiniteArguments()
        {
            foreach (double delta in new[] { -3d, 0d, 2d, 45d })
            {
                var t = new NoncentralT(7d, delta);
                Assert.AreEqual(1d, t.CDF(double.PositiveInfinity), 0d);
                Assert.AreEqual(0d, t.CDF(double.NegativeInfinity), 0d);
                Assert.AreEqual(0d, t.CCDF(double.PositiveInfinity), 0d);
                Assert.AreEqual(1d, t.CCDF(double.NegativeInfinity), 0d);
                Assert.AreEqual(0d, t.PDF(double.PositiveInfinity), 0d);
                Assert.AreEqual(0d, t.PDF(double.NegativeInfinity), 0d);
                Assert.IsTrue(double.IsNaN(t.CDF(double.NaN)));
                Assert.IsTrue(double.IsNaN(t.CCDF(double.NaN)));
                Assert.IsTrue(double.IsNaN(t.PDF(double.NaN)));
                Assert.Throws<ArgumentOutOfRangeException>(() => t.InverseCDF(double.NaN));
            }
        }

        /// <summary>
        /// The density at zero is finite for large degrees of freedom and matches high-precision values of
        /// Γ((ν + 1)/2) e^(−μ²/2) / (√(πν) Γ(ν/2)). Reference: scipy.stats.nct.pdf(0, 400, 1) = 0.24181954022381516.
        /// </summary>
        /// <remarks>
        /// The tabulated values were computed with Python mpmath 1.4.1 at 50 significant digits; the tolerance is 1E-15
        /// relative.
        /// </remarks>
        [TestMethod]
        public void Test_NoncentralT_PDFAtZero()
        {
            AssertRelative(0.24181954022381516, new NoncentralT(400d, 1d).PDF(0d), 1E-13, "v=400");
            Assert.AreEqual(new StudentT(5000d).PDF(0d) * Math.Exp(-0.5d), new NoncentralT(5000d, 1d).PDF(0d), 1E-13);
            // { v, u, f(0) }
            var table = new double[,]
            {
                { 1, -2, 0.04307855860369726 }, { 1, 0, 0.3183098861837907 }, { 1, 1.5, 0.10334008993412312 },
                { 2.5, -2, 0.04896548614401862 }, { 2.5, 0, 0.36180872402956504 }, { 2.5, 1.5, 0.11746209496797452 },
                { 10, -2, 0.05266009335378346 }, { 10, 0, 0.38910838396603104 }, { 10, 1.5, 0.1263249969243921 },
                { 100, -2, 0.05385615992189921 }, { 100, 0, 0.3979461869358938 }, { 100, 1.5, 0.129194211464585 },
                { 340, -2, 0.0539512819220309 }, { 340, 0, 0.39864904873110935 }, { 340, 1.5, 0.12942239728061364 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double v = table[i, 0], u = table[i, 1];
                Assert.AreEqual(table[i, 2], new NoncentralT(v, u).PDF(0d), 1E-15 * table[i, 2], $"v={v}, u={u}");
            }
        }

        /// <summary>
        /// Tail and body values to near machine precision. Reference CDF values are the Guenther series summed with
        /// mpmath at 80 digits (for x &lt; 0 as 1 − F(−x; ν, −μ) at 400 digits); reference densities are the confluent
        /// hypergeometric closed form evaluated with mpmath. scipy.stats.nct agrees to its own accuracy, for example
        /// scipy.stats.nct.cdf(-5, 10, 2) = 2.5188462768697946E-07 (relative error 1E-11).
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_MatchesHighPrecisionOracle()
        {
            // { x, ν, μ, PDF, CDF, CCDF }
            var table = new double[,]
            {
                { -5, 10, 2, 4.2222480713063658E-7, 2.5188462768958445E-7, 0.99999974811537231 },
                { 60, 10, 2, 8.598197149800068E-13, 0.99999999999481749, 5.1825136361448232E-12 },
                { 1E-9, 10, 2, 0.05266009346176539, 0.022750132000839301, 0.9772498679991607 },
                { -1E-9, 10, 2, 0.052660093245801524, 0.022750131895519114, 0.97724986810448089 },
                { 1.4, 4, 2.42, 0.23552141801722211, 0.15955740663227623, 0.84044259336772377 },
                { -2, 9, 1.75, 0.00083062247363953004, 0.00035703726877226363, 0.99964296273122774 },
                { 3, 1000, 1.5, 0.12970229767045064, 0.93265888363915856, 0.067341116360841441 },
                { -0.5, 2.5, 3, 0.00083470115575142798, 0.00035293391112280933, 0.99964706608887719 },
                { -1E-9, 1, 1, 0.19306470501813709, 0.15865525373839235, 0.84134474626160765 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                var t = new NoncentralT(table[i, 1], table[i, 2]);
                double x = table[i, 0];
                AssertRelative(table[i, 3], t.PDF(x), 1E-12, $"PDF row {i}");
                AssertRelative(Math.Log(table[i, 3]), t.LogPDF(x), 1E-13, $"LogPDF row {i}");
                AssertRelative(table[i, 4], t.CDF(x), 1E-12, $"CDF row {i}");
                AssertRelative(table[i, 5], t.CCDF(x), 1E-12, $"CCDF row {i}");
            }
        }

        /// <summary>
        /// Extreme quantiles. Reference values are roots of the mpmath Guenther series; scipy.stats.nct.ppf(1E-10, 10, 2)
        /// = -11.768545467517564 and scipy.stats.nct.ppf(0.999999, 10, 2) = 17.336750141870535; R qt(1E-10, 10, 2)
        /// returns -11.7712 because its lower tail is formed by subtraction.
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_ExtremeQuantiles()
        {
            var t = new NoncentralT(10d, 2d);
            AssertRelative(-11.768545467626908, t.InverseCDF(1E-10), 1E-12, "p = 1E-10");
            AssertRelative(17.336750141870533, t.InverseCDF(0.999999), 1E-12, "p = 0.999999");
            foreach (double p in new[] { 1E-30, 1E-5, 0.3, 0.9, 1 - 1E-12 })
            {
                double x = t.InverseCDF(p);
                double actual = p < 0.5 ? t.CDF(x) : t.CCDF(x);
                double target = p < 0.5 ? p : 1 - p;
                AssertRelative(target, actual, 1E-9, $"Round trip p = {p}");
            }
        }

        /// <summary>
        /// Maximum likelihood for the degrees of freedom and noncentrality attains at least the log-likelihood of
        /// scipy.stats.nct.fit(sample, floc=0, fscale=1), which gives df = 8.172735310413508, nc = 2.3277100431320847
        /// and log-likelihood -31.906983722436404 (a tight refinement reaches -31.906983722423817 at
        /// df = 8.172751328417844, nc = 2.327709174283645).
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_MLE_MatchesScipyFit()
        {
            var sample = new[] { 1.9, 2.8, 1.2, 3.5, 2.2, 0.7, 2.9, 4.8, 1.6, 2.4, 3.1, 1.8, 2.0, 6.3, 2.6, 1.4, 2.3, 3.9, 0.9, 2.7 };
            var t = new NoncentralT();
            t.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            Assert.IsTrue(t.ParametersValid);
            Assert.IsGreaterThanOrEqualTo(-31.906983722436404 - 1E-6, t.LogLikelihood(sample));
            Assert.AreEqual(8.172751328417844, t.DegreesOfFreedom, 1E-2);
            Assert.AreEqual(2.327709174283645, t.Noncentrality, 1E-3);
        }

        /// <summary>
        /// Estimation input validation, the unsupported methods, and bootstrapping.
        /// </summary>
        [TestMethod]
        public void Test_NoncentralT_Estimation_ValidationAndBootstrap()
        {
            var t = new NoncentralT(12d, 1.5d);
            Assert.Throws<NotImplementedException>(() => t.Estimate(new[] { 1d, 2d, 3d, 5d, 8d }, ParameterEstimationMethod.MethodOfMoments));
            Assert.Throws<ArgumentOutOfRangeException>(() => t.Estimate(new[] { 1d, 2d, 3d }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.Throws<ArgumentOutOfRangeException>(() => t.GetParameterConstraints(new[] { 1d, double.PositiveInfinity, 3d, 4d }));
            var constraints = t.GetParameterConstraints(new[] { 1.9, 2.8, 1.2, 3.5, 2.2, 0.7, 2.9, 4.8 });
            for (int i = 0; i < 2; i++)
            {
                Assert.IsTrue(constraints.Item2[i] <= constraints.Item1[i] && constraints.Item1[i] <= constraints.Item3[i]);
            }
            var bootstrapped = (NoncentralT)t.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 100, 12345);
            Assert.IsTrue(bootstrapped.ParametersValid);
            Assert.AreEqual(1.5d, bootstrapped.Noncentrality, 0.75d);
        }
    }
}
