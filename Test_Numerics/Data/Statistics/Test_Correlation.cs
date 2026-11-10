using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the Correlation class. These methods were tested against values attained from R's "cor()" function from the "stats" package.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, 
    /// Vienna, Austria. ISBN 3-900051-07-0, URL <see href=" http://www.R-project.org/."/>
    /// </remarks>
    [TestClass]
    public class Test_Correlation
    {
        /// <summary>
        /// Test the Pearson correlation coefficient on a small dataset
        /// </summary>
        [TestMethod()]
        public void Test_Pearson()
        {
            var XArray = new double[] { 14d, 8d, 32d, 7d, 3d, 15d };
            var YArray = new double[] { 10d, 5d, 7d, 4d, 3d, 8d };
            double R = Correlation.Pearson(XArray, YArray);
            double trueR = 0.54502739907793d;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// The corrected two-pass sums keep the Pearson coefficient accurate when the means are large relative to the
        /// spread. The exact value is from fractions.Fraction sums and a 60-digit mpmath square root; the uncorrected
        /// two-pass sums are about 1e3 ulps off on this fixture. The matrix overload, scalar and vectorized, agrees with it.
        /// </summary>
        [TestMethod]
        public void Test_Pearson_LargeOffset()
        {
            int n = 500;
            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = 1e9 + (i % 10) * 0.1;
                y[i] = 2e9 - (i % 7) * 0.3 + (i % 10) * 0.05;
            }
            double exact = 0.22979238343446176;
            double tolerance = 2 * 2.7755575615628914e-17; // two ulps of the exact value, 2^-55 each
            Assert.AreEqual(exact, Correlation.Pearson(x, y), tolerance);

            foreach (int p in new[] { 2, 9 })
            {
                var samples = new double[n, p];
                for (int i = 0; i < n; i++)
                {
                    samples[i, 0] = x[i];
                    samples[i, 1] = y[i];
                    for (int j = 2; j < p; j++) samples[i, j] = 1e6 * j + ((i * j) % 13) * 0.25;
                }
                var r = Correlation.Pearson(samples);
                Assert.AreEqual(exact, r[0, 1], tolerance, $"p = {p}");
                Assert.AreEqual(r[0, 1], r[1, 0], 0d);
                Assert.AreEqual(1d, r[0, 0], 1e-15);
            }
        }

        /// <summary>
        /// Test the Pearson correlation coefficient on a larger dataset
        /// </summary>
        [TestMethod]
        public void Test_Pearson_Big()
        {
            var XArray = new double[] { 230408, 288010, 345611, 403213, 460815, 518417, 576019, 633612, 691223, 748825, 806427, 864029, 921631, 1036834, 1152038 };
            var YArray = new double[] { 1519.7, 1520.5, 1520.9, 1521.7, 1523.5, 1525.9, 1528.4, 1530.9, 1533.2, 1534.7, 1535.9, 1538, 1541.3, 1547.7, 1552.7 };
            double R = Correlation.Pearson(XArray, YArray);
            double trueR = 0.988054377242161;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Test the Spearman ranked correlation coefficient on a small dataset
        /// </summary>
        [TestMethod()]
        public void Test_Spearman()
        {
            var XArray = new double[] { 14d, 8d, 32d, 7d, 3d, 15d };
            var YArray = new double[] { 10d, 5d, 7d, 4d, 3d, 8d };
            double R = Correlation.Spearman(XArray, YArray);
            double trueR = 0.771428571428571d;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Test the Spearman ranked correlation coefficient on a larger dataset
        /// </summary>
        [TestMethod]
        public void Test_Spearman_Big()
        {
            var XArray = new double[] { 230408, 288010, 345611, 403213, 460815, 518417, 576019, 633612, 691223, 748825, 806427, 864029, 921631, 1036834, 1152038 };
            var YArray = new double[] { 1519.7, 1520.5, 1520.9, 1521.7, 1523.5, 1525.9, 1528.4, 1530.9, 1533.2, 1534.7, 1535.9, 1538, 1541.3, 1547.7, 1552.7 };
            double R = Correlation.Spearman(XArray, YArray);
            double trueR = 1;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Test the Kendall's Tau ranked correlation coefficient on a small dataset
        /// </summary>
        [TestMethod()]
        public void Test_KendallsTau()
        {
            var XArray = new double[] { 14d, 8d, 32d, 7d, 3d, 15d };
            var YArray = new double[] { 10d, 5d, 7d, 4d, 3d, 8d };
            double R = Correlation.KendallsTau(XArray, YArray);
            double trueR = 0.6d;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Knight's O(n log n) path reproduces the O(n²) pair loop bit for bit on untied, tied, heavily tied,
        /// signed-zero and constant samples at and around the size thresholds.
        /// </summary>
        [TestMethod]
        public void Test_KendallsTau_KnightMatchesPairLoop()
        {
            var rng = new System.Random(12345);
            foreach (int n in new[] { Correlation.KnightMinimumSize, Correlation.KnightMinimumSize + 1, 100, 257, 1000, 3001 })
            {
                for (int scenario = 0; scenario < 6; scenario++)
                {
                    var x = new double[n];
                    var y = new double[n];
                    for (int i = 0; i < n; i++)
                    {
                        double u = rng.NextDouble(), v = rng.NextDouble();
                        switch (scenario)
                        {
                            case 0: x[i] = u; y[i] = u + v; break;                                  // untied, correlated
                            case 1: x[i] = System.Math.Floor(10 * u); y[i] = System.Math.Floor(5 * v); break; // ties in both
                            case 2: x[i] = System.Math.Floor(2 * u); y[i] = System.Math.Floor(3 * (u + v)); break; // heavy ties
                            case 3: x[i] = u < 0.5 ? 0d : -0d; y[i] = v < 0.3 ? -0d : v; break;    // signed zeros
                            case 4: x[i] = 7d; y[i] = v; break;                                       // constant sample
                            case 5: x[i] = 1E6 * (u - 0.5); y[i] = -System.Math.Round(x[i], 1) + 1E-3 * v; break; // negative association
                        }
                    }
                    double fast = Correlation.KendallsTau(x, y);
                    double pairs = Correlation.KendallsTauPairs(x, y);
                    Assert.AreEqual(System.BitConverter.DoubleToInt64Bits(pairs), System.BitConverter.DoubleToInt64Bits(fast),
                        $"n = {n}, scenario = {scenario}: {fast} vs {pairs}");
                }
            }
        }

        /// <summary>
        /// Tau-b on a tied sample of 40 pairs (Knight path) against scipy.stats.kendalltau 1.17.1 (variant 'b'),
        /// -0.014064711520113382.
        /// </summary>
        [TestMethod]
        public void Test_KendallsTau_TiedScipy()
        {
            var x = new double[] { 0, 7, 5, 3, 1, 8, 6, 4, 2, 0, 7, 5, 3, 1, 8, 6, 4, 2, 0, 7, 5, 3, 1, 8, 6, 4, 2, 0, 7, 5, 3, 1, 8, 6, 4, 2, 0, 7, 5, 3 };
            var y = new double[] { 0, 2, 6, 9, 6, 5, 3, 6, 11, 4, 2, 2, 1, 5, 11, 5, 4, 5, 5, 10, 6, 1, 1, 3, 4, 10, 7, 3, 4, 7, 9, 5, 3, 0, 2, 6, 9, 6, 5, 3 };
            Assert.AreEqual(-0.014064711520113382, Correlation.KendallsTau(x, y), 1E-15);
        }

        /// <summary>
        /// Edge cases against scipy.stats.kendalltau 1.17.1 and analytic values: a NaN gives NaN (scipy's default
        /// propagation), infinite values rank exactly, values whose differences underflow or overflow keep their
        /// order, and a sample beyond 65,536 values, whose pair count exceeds a 32-bit counter, is still exact.
        /// </summary>
        [TestMethod]
        public void Test_KendallsTau_EdgeCases()
        {
            // NaN in either sample, on the pair path and on the Knight path.
            Assert.IsTrue(double.IsNaN(Correlation.KendallsTau(new[] { 1d, 2d, 3d, double.NaN }, new[] { 1d, 2d, 3d, 4d })));
            var nanX = new double[60]; var y60 = new double[60];
            for (int i = 0; i < 60; i++) { nanX[i] = i; y60[i] = (i * 17) % 23; }
            nanX[7] = double.NaN;
            Assert.IsTrue(double.IsNaN(Correlation.KendallsTau(nanX, y60)));

            // Infinite values: scipy gives 1/3 for these four pairs; the same pairs repeated 15 times (Knight path) keep
            // the pattern of ties, so the analytic tau-b follows from the counts.
            Assert.AreEqual(0.3333333333333334, Correlation.KendallsTau(new[] { 1d, 2d, double.PositiveInfinity, 4d }, new[] { 1d, 2d, 3d, 0d }), 1E-15);
            var infX = new double[60]; var infY = new double[60];
            for (int i = 0; i < 60; i++) { infX[i] = i == 59 ? double.PositiveInfinity : i; infY[i] = i == 0 ? double.NegativeInfinity : i; }
            Assert.AreEqual(1d, Correlation.KendallsTau(infX, infY), 1E-15);

            // Differences that underflow (1e-170 squared) or overflow (a range of 2e308) are ranked by comparison.
            var tinyX = new double[60]; var tinyY = new double[60];
            for (int i = 0; i < 60; i++) { tinyX[i] = i * 1E-170; tinyY[i] = i * 1E-170; }
            Assert.AreEqual(1d, Correlation.KendallsTau(tinyX, tinyY), 1E-15);
            Assert.AreEqual(1d, Correlation.KendallsTau(new[] { 1E-170, 2E-170, 3E-170 }, new[] { 1E-170, 2E-170, 3E-170 }), 1E-15);
            Assert.AreEqual(-1d, Correlation.KendallsTau(new[] { -1E308, 0d, 1E308 }, new[] { 3d, 2d, 1d }), 1E-15);

            // n = 70,000: x_i = (7919 i mod 70001) / 70001, y_i = x_i + (104729 i mod 70003) / 70003; scipy 0.5005251658697164.
            int n = 70000;
            var x = new double[n]; var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = (i * 7919L % 70001) / 70001.0;
                y[i] = x[i] + (i * 104729L % 70003) / 70003.0;
            }
            Assert.AreEqual(0.5005251658697164, Correlation.KendallsTau(x, y), 1E-15);
        }

        /// <summary>
        /// On finite samples of ordinary magnitude the comparison-based pair loop reproduces the product-based pair
        /// loop of the 2.2.0 implementation, kept as a bitwise reference, bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_KendallsTau_MatchesLegacyLoop()
        {
            var rng = new System.Random(2024);
            foreach (int n in new[] { 2, 5, 20, 47, 48, 300 })
            {
                for (int scenario = 0; scenario < 3; scenario++)
                {
                    var x = new double[n];
                    var y = new double[n];
                    for (int i = 0; i < n; i++)
                    {
                        double u = rng.NextDouble(), v = rng.NextDouble();
                        x[i] = scenario == 0 ? u : System.Math.Floor(6 * u);
                        y[i] = scenario == 2 ? System.Math.Floor(4 * v) : u + v;
                    }
                    Assert.AreEqual(System.BitConverter.DoubleToInt64Bits(LegacyKendallsTau(x, y)),
                        System.BitConverter.DoubleToInt64Bits(Correlation.KendallsTau(x, y)), $"n = {n}, scenario = {scenario}");
                }
            }
        }

        /// <summary>
        /// The 2.2.0 pair loop (32-bit counters, product of differences), kept as a bitwise reference.
        /// </summary>
        private static double LegacyKendallsTau(double[] sample1, double[] sample2)
        {
            int n = sample1.Length, i = 0, n1 = 0, n2 = 0;
            for (int j = 0; j <= n - 2; j++)
            {
                for (int k = j + 1; k <= n - 1; k++)
                {
                    double a1 = sample1[j] - sample1[k];
                    double a2 = sample2[j] - sample2[k];
                    double aa = a1 * a2;
                    if (aa != 0.0d)
                    {
                        n1 += 1;
                        n2 += 1;
                        i = aa > 0.0 ? ++i : --i;
                    }
                    else
                    {
                        if (a1 != 0.0d) n1 += 1;
                        if (a2 != 0.0d) n2 += 1;
                    }
                }
            }
            return i / (System.Math.Sqrt(n1) * System.Math.Sqrt(n2));
        }

        /// <summary>
        /// Test the Kendall's Tau ranked correlation coefficient on a larger dataset
        /// </summary>
        [TestMethod]
        public void Test_KendallsTau_Big()
        {
            var XArray = new double[] { 230408, 288010, 345611, 403213, 460815, 518417, 576019, 633612, 691223, 748825, 806427, 864029, 921631, 1036834, 1152038 };
            var YArray = new double[] { 1519.7, 1520.5, 1520.9, 1521.7, 1523.5, 1525.9, 1528.4, 1530.9, 1533.2, 1534.7, 1535.9, 1538, 1541.3, 1547.7, 1552.7 };
            double R = Correlation.KendallsTau(XArray, YArray);
            double trueR = 1;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }
    }
}
