using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Sampling;

namespace Data.Statistics
{
    /// <summary>
    /// Unit testing for HypothesisTest and MultipleGrubbsBeckTest classes. Most of these methods were validated with methods from R's "stats" package and known examples.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing,
    /// Vienna, Austria. ISBN 3-900051-07-0, URL http://www.R-project.org/.
    /// </item>
    /// <item>
    /// "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee and F. Ashkar, Water Resources Publications, 1991.
    /// </item>
    /// </list>
    /// </remarks>
    [TestClass]
    public partial class Test_HypothesisTests
    {
        /// <summary>
        /// Test the OneSampleTtest method against known examples and R's "t.test()" function from the "stats" package.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// <see href="http://www.real-statistics.com/students-t-distribution/one-sample-t-test/"/>
        /// </remarks>
        [TestMethod]
        public void Test_OneSampleTtest()
        {
            var data = new double[] { 8.782932, 10.64199, -1.63955, -6.802458, 9.088312, -27.26934, 9.451478, -4.142762, -4.262396, -13.78983, -1.743717, 27.259681, 5.559418, 7.803247, -11.25798, 12.253498, -13.295363, -4.973664, 16.81069, 4.480855, 11.694329, 21.836776, -9.664926, -23.297061, -23.965643, 27.076463, -7.22471, 9.305697, 9.181852, -2.434665 };
            var p = HypothesisTests.OneSampleTtest(data);
            double true_p = 0.6489;
            Assert.AreEqual(true_p, p, 1E-4);

            p = HypothesisTests.OneSampleTtest(data, 10);
            true_p = 0.001823;
            Assert.AreEqual(true_p, p, 1E-4);

            var t = HypothesisTests.OneSampleTtest(new double[] { 23, 15, -5, 7, 1, -10, 12, -8, 20, 8, -2, -5 });
            Assert.AreEqual(0.087585 * 2, t, 1E-6);
        }

        /// <summary>
        /// Test the equal variance T-test method against R's "t.test()" function from the "stats" package
        /// </summary>
        [TestMethod]
        public void Test_EqualVarianceTtest()
        {
            var data1 = new double[] { 8.782932, 10.64199, -1.63955, -6.802458, 9.088312, -27.26934, 9.451478, -4.142762, -4.262396, -13.78983, -1.743717, 27.259681, 5.559418, 7.803247, -11.25798, 12.253498, -13.295363, -4.973664, 16.81069, 4.480855, 11.694329, 21.836776, -9.664926, -23.297061, -23.965643, 27.076463, -7.22471, 9.305697, 9.181852, -2.434665 };
            var data2 = new double[] { 34.3561954, 75.9050064, 71.4757101, 58.9733692, 17.6281358, 24.7356484, 0.2774026, -39.8615073, 63.0320155, 10.7740315, 43.855325, -61.4107418, -21.8079666, 38.1142162, 35.6335516, 53.8218821, -32.3929633, 27.0220976, 27.4956296, -29.203965, -6.2115822, 68.4307799, 11.6077081, 20.5498852, -10.1292962, 18.3386108, 30.7351382, 34.7138599, 74.3519506, -60.4083194 };

            var p = HypothesisTests.EqualVarianceTtest(data1, data2);
            double true_p = 0.0185;
            Assert.AreEqual(true_p, p, 1E-4);
        }

        /// <summary>
        /// Two-sided t p-values far in the tail keep their relative accuracy (the complement 1 - CDF floored them near
        /// 2E-16), against scipy.stats 1.17.1 ttest_1samp, ttest_ind (pooled and Welch) and ttest_rel, and R 4.4.3
        /// summary(lm(y ~ x)) for the linear trend test.
        /// </summary>
        [TestMethod]
        public void Test_Ttests_FarTail()
        {
            var o = new double[] { 10, 10.001, 9.999, 10.0005, 9.9995, 10.0002, 9.9998, 10.0001 };
            Assert.AreEqual(5.609558661898506e-31, HypothesisTests.OneSampleTtest(o, 0), 5.609558661898506e-31 * 1E-8);

            var a = new double[] { 10, 10.001, 9.999, 10.0005 };
            var b = new double[] { 20, 20.001, 19.999, 20.0005, 19.9995 };
            Assert.AreEqual(3.967421244871032e-28, HypothesisTests.EqualVarianceTtest(a, b), 3.967421244871032e-28 * 1E-8);
            Assert.AreEqual(1.5142847333611348e-25, HypothesisTests.UnequalVarianceTtest(a, b), 1.5142847333611348e-25 * 1E-8);
            var c = new double[] { 15, 15.0011, 14.9992, 15.0008 };
            Assert.AreEqual(4.744656587997517e-15, HypothesisTests.PairedTtest(a, c), 4.744656587997517e-15 * 1E-6);

            var time = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            var y = new double[] { 2.1, 3.9, 6.2, 7.8, 10.1, 12.2, 13.8, 16.1, 18.0, 20.3 };
            Assert.AreEqual(8.7635426585689750e-14, HypothesisTests.LinearTrendTest(time, y), 8.7635426585689750e-14 * 1E-6);

            // The two-sided tail 2 P(T > |t|) the tests use: t = 0 gives 1, an infinite t gives 0, NaN propagates,
            // ν < 1 is refused.
            Assert.AreEqual(1d, 2d * new StudentT(5d).CCDF(0d));
            Assert.AreEqual(0d, 2d * new StudentT(5d).CCDF(double.PositiveInfinity));
            Assert.IsTrue(double.IsNaN(2d * new StudentT(5d).CCDF(double.NaN)));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => new StudentT(0.5d).CCDF(1d));
        }

        /// <summary>
        /// A single-value sample adds nothing to the pooled sum of squares: scipy.stats.ttest_ind([1], [1, 2, 3]) 1.17.1
        /// gives p = 0.47776703213290644, in either argument order.
        /// </summary>
        [TestMethod]
        public void Test_EqualVarianceTtest_SingleValueSample()
        {
            Assert.AreEqual(0.47776703213290644, HypothesisTests.EqualVarianceTtest(new double[] { 1 }, new double[] { 1, 2, 3 }), 1E-12);
            Assert.AreEqual(0.47776703213290644, HypothesisTests.EqualVarianceTtest(new double[] { 1, 2, 3 }, new double[] { 1 }), 1E-12);
        }

        /// <summary>
        /// Test the unequal variance T-test method against R's "t.test()" function from the "stats" package
        /// </summary>
        [TestMethod]
        public void Test_UnequalVarianceTtest()
        {
            var data1 = new double[] { 8.782932, 10.64199, -1.63955, -6.802458, 9.088312, -27.26934, 9.451478, -4.142762, -4.262396, -13.78983, -1.743717, 27.259681, 5.559418, 7.803247, -11.25798, 12.253498, -13.295363, -4.973664, 16.81069, 4.480855, 11.694329, 21.836776, -9.664926, -23.297061, -23.965643, 27.076463, -7.22471, 9.305697, 9.181852, -2.434665 };
            var data2 = new double[] { 34.3561954, 75.9050064, 71.4757101, 58.9733692, 17.6281358, 24.7356484, 0.2774026, -39.8615073, 63.0320155, 10.7740315, 43.855325, -61.4107418, -21.8079666, 38.1142162, 35.6335516, 53.8218821, -32.3929633, 27.0220976, 27.4956296, -29.203965, -6.2115822, 68.4307799, 11.6077081, 20.5498852, -10.1292962, 18.3386108, 30.7351382, 34.7138599, 74.3519506, -60.4083194 };

            var p = HypothesisTests.UnequalVarianceTtest(data1, data2);
            double true_p = 0.02043;
            Assert.AreEqual(true_p, p, 1E-4);
        }

        /// <summary>
        /// Test the paired T-test method against R's "t.test()" function from the "stats" package
        /// </summary>
        [TestMethod]
        public void Test_PairedTtest()
        {
            var data1 = new double[] { 200.1, 190.9, 192.7, 213, 241.4, 196.9, 172.2, 185.5, 205.2, 193.7 };
            var data2 = new double[] { 392.9, 393.2, 345.1, 393, 434, 427.9, 422, 383.9, 392.3, 352.2 };

            var p = HypothesisTests.PairedTtest(data1, data2);
            double true_p = 6.2E-9;
            Assert.AreEqual(true_p, p, 1E-4);
        }

        /// <summary>
        /// Test the F-test method against "var.test()" function from the "stats" package
        /// </summary>
        [TestMethod]
        public void Test_Ftest()
        {
            var data1 = new double[] { 200.1, 190.9, 192.7, 213, 241.4, 196.9, 172.2, 185.5, 205.2, 193.7 };
            var data2 = new double[] { 392.9, 393.2, 345.1, 393, 434, 427.9, 422, 383.9, 392.3, 352.2 };

            var p = HypothesisTests.Ftest(data1, data2);
            double true_p = 0.1825;
            Assert.AreEqual(true_p, p, 1E-4);
        }

        /// <summary>
        /// Test the FtestModels method against a known example.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// <see href="https://online.stat.psu.edu/stat501/lesson/6/6.2"/>
        /// </remarks>
        [TestMethod]
        public void Test_FtestModels()
        {
            // full model has more parameters; reduced model has fewer
            // f stat decided whether or not to reject the smaller reduced model in favor of the larger full model
            double fStat = 0; double pValue = 0;

            HypothesisTests.FtestModels(1224.32, 720.27, 49, 48, out fStat, out pValue);
            double trueFVal = 33.5899;
            double truePVal = 0;

            Assert.AreEqual(trueFVal, fStat, 1E-3);
            Assert.AreEqual(truePVal, pValue, 1E-6);
            // scipy.stats.f.sf(33.5899..., 1, 48) = 5.134738187407378e-07.
            Assert.AreEqual(5.134738187407378e-07, pValue, 5.134738187407378e-07 * 1E-8);
        }

        /// <summary>
        /// The nested-model F test is one-sided with df_R − df_F numerator degrees of freedom: scipy.stats.f.sf(5.25,
        /// 2, 7) = 0.04047715405015526 for SSE 100 and 40 on 9 and 7 degrees of freedom. A statistic of zero gives 1, and
        /// reversed degrees of freedom are refused.
        /// </summary>
        [TestMethod]
        public void Test_FtestModels_NumeratorDegreesOfFreedom()
        {
            HypothesisTests.FtestModels(100, 40, 9, 7, out double fStat, out double pValue);
            Assert.AreEqual(5.25, fStat, 1E-12);
            Assert.AreEqual(0.04047715405015526, pValue, 1E-14);

            HypothesisTests.FtestModels(40, 40, 9, 7, out fStat, out pValue);
            Assert.AreEqual(0d, fStat);
            Assert.AreEqual(1d, pValue);

            Assert.ThrowsExactly<System.ArgumentException>(() => HypothesisTests.FtestModels(40, 100, 7, 9, out _, out _));
        }

        /// <summary>
        /// Test the Jarque-Bera method against R's "jarque.test()" method from the "moments" package and a known example.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// <list type="bullet">
        /// <item><description>
        /// Lukasz Komsta (2005). moments: Moments, Cumulants, Skewness, Kurtosis and Related Tests. R package version 0.14.1, https://cran.r-project.org/web/packages/moments
        /// </description></item>
        /// <item><description>
        /// <see href="https://www.statology.org/jarque-bera-test-excel/"/>
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_JarqueBera()
        {
            var data = new double[] { -17.82175266, -2.33394663, 4.66366786, 6.77181741, 48.09893105, -28.26940615, -9.98265593, -0.87518792, 14.97758789, -0.54200675, 8.80374205, -3.40846222, -3.35109891, -12.98362149, -1.42481547, 18.5800533, 10.86238267, -13.65904345, 1.76995771, 13.91485418, 10.8528196, -11.69442361, -11.6048953, 5.89082943, -13.20258835, 3.93329214, 2.62990935, -4.00680666, 18.4215721, 0.14773234, 10.20778973, -7.41284797, 8.42081407, 35.41116192, 59.37166512, 9.36721778, 22.37395361, 22.9971476, 13.47067667, -18.98066759, -22.84094314, 16.7108515, 18.72618308, 29.97227498, -16.078326, -0.26901107, -0.05773469, 14.44571902, -7.23727541, 18.87940528, -10.55665291, -0.40463948, 1.1599797, -9.47746043, -8.83651712, -0.1277879, -7.43500345, 18.02267959, 10.38996171, 6.73008507, -0.78965999, -6.63662283, -0.26534812, -18.26597299, 10.68284417, 4.14715065, -22.73605154, 0.38107214, 17.99480125, 4.67217999, -7.55979566, -5.02964486, 10.07853161, -10.20580542, -3.83664015, -3.13645528, -4.30412819, 12.06651361, 16.46249676, -0.77303738, 10.72787315, 12.09162065, 8.22959713, 5.86544228, -11.14598952, 9.55434186, -4.24740884, -2.84574008, -7.08625811, 0.80619592, 12.92545548, -3.22668772, -25.39204102, 9.92546076, -3.16982112, 18.60432604, -14.00214643, 1.17374306, -13.04390662, 24.21704845, 3.82716675, -5.17619789, -8.06288031, 4.1033081, -13.36564786, -15.91238602, -25.39452748, 18.80121063, 7.80923857, -6.8516946, 6.54494797, 26.80612853, 4.65921504, 23.73597901, 44.15782916, 2.64243694, -24.27815428, 40.02096079, 6.5730404, 25.71086816, 28.14206721, 55.47192889, -14.9762203, -10.03718017, 38.03084527, 7.20256442 };

            // R moments::jarque.test gives 3.444E-05; scipy.stats.jarque_bera 1.17.1 gives 3.443778022828376e-05.
            var p = HypothesisTests.JarqueBeraTest(data);
            double true_p = 3.443778022828376e-05;
            Assert.AreEqual(true_p, p, true_p * 1E-8);

            // scipy.stats.jarque_bera: statistic 1.0040277490351397, p-value 0.605310412198013. (The earlier 0.592128 came
            // from an example that used the bias-corrected skewness and kurtosis.)
            var JB = HypothesisTests.JarqueBeraTest(new double[] { 4, 5, 5, 6, 9, 12, 13, 14, 14, 19, 22, 24, 25 });
            Assert.AreEqual(0.605310412198013, JB, 1E-12);

            // scipy.stats.jarque_bera: p-value 0.6910990901340825.
            Assert.AreEqual(0.6910990901340825, HypothesisTests.JarqueBeraTest(new double[] { 2.1, 3.5, 1.2, 8.9, 4.4, 5.0, 3.3, 2.8, 7.7, 6.1, 0.9, 4.2 }), 1E-12);
        }


        /// <summary>
        /// Test the Wald-Wolfowitz method against a known example from "The Gamma Family..."
        /// </summary>
        [TestMethod]
        public void Test_WaldWolfowitz()
        {
            // Table 1.2 Maximum annual peak discharge values in cms, observed at the Harricana River at Amos (Quebec, Canada)
            var data = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

            var p = HypothesisTests.WaldWolfowitzTest(data);
            double true_p = (1 - Normal.StandardCDF(1.167)) * 2d; // See page 5 of reference.
            Assert.AreEqual(true_p, p, 1E-3);
        }

        /// <summary>
        /// Test the Ljung-Box method against R's "Box.test()" method from the "stats" package.
        /// </summary>
        [TestMethod]
        public void Test_LjungBox()
        {
            var data = new double[] { -17.82175266, -2.33394663, 4.66366786, 6.77181741, 48.09893105, -28.26940615, -9.98265593, -0.87518792, 14.97758789, -0.54200675, 8.80374205, -3.40846222, -3.35109891, -12.98362149, -1.42481547, 18.5800533, 10.86238267, -13.65904345, 1.76995771, 13.91485418, 10.8528196, -11.69442361, -11.6048953, 5.89082943, -13.20258835, 3.93329214, 2.62990935, -4.00680666, 18.4215721, 0.14773234, 10.20778973, -7.41284797, 8.42081407, 35.41116192, 59.37166512, 9.36721778, 22.37395361, 22.9971476, 13.47067667, -18.98066759, -22.84094314, 16.7108515, 18.72618308, 29.97227498, -16.078326, -0.26901107, -0.05773469, 14.44571902, -7.23727541, 18.87940528, -10.55665291, -0.40463948, 1.1599797, -9.47746043, -8.83651712, -0.1277879, -7.43500345, 18.02267959, 10.38996171, 6.73008507, -0.78965999, -6.63662283, -0.26534812, -18.26597299, 10.68284417, 4.14715065, -22.73605154, 0.38107214, 17.99480125, 4.67217999, -7.55979566, -5.02964486, 10.07853161, -10.20580542, -3.83664015, -3.13645528, -4.30412819, 12.06651361, 16.46249676, -0.77303738, 10.72787315, 12.09162065, 8.22959713, 5.86544228, -11.14598952, 9.55434186, -4.24740884, -2.84574008, -7.08625811, 0.80619592, 12.92545548, -3.22668772, -25.39204102, 9.92546076, -3.16982112, 18.60432604, -14.00214643, 1.17374306, -13.04390662, 24.21704845, 3.82716675, -5.17619789, -8.06288031, 4.1033081, -13.36564786, -15.91238602, -25.39452748, 18.80121063, 7.80923857, -6.8516946, 6.54494797, 26.80612853, 4.65921504, 23.73597901, 44.15782916, 2.64243694, -24.27815428, 40.02096079, 6.5730404, 25.71086816, 28.14206721, 55.47192889, -14.9762203, -10.03718017, 38.03084527, 7.20256442 };

            var p = HypothesisTests.LjungBoxTest(data,5);
            double true_p = 0.2314;
            Assert.AreEqual(true_p, p, 1E-4);

            var p2 = HypothesisTests.LjungBoxTest(data, 30);
            var true_p2 = 0.7548;
            Assert.AreEqual(true_p2, p2, 1E-4);
        }


        /// <summary>
        /// Test the Mann-Whitney method against a known example from "The Gamma Family..."
        /// </summary>
        [TestMethod]
        public void Test_MannWhitney()
        {
            // Table 1.2 Maximum annual peak discharge values in cms, observed at the Harricana River at Amos (Quebec, Canada)
            var data1 = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d };
            var data2 = new double[] { 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

            var p = HypothesisTests.MannWhitneyTest(data2, data1);
            double true_p = (1 - Normal.StandardCDF(0.54)) * 2d; // See page 7 of reference.
            Assert.AreEqual(true_p, p, 1E-2);
        }

        /// <summary>
        /// Test the Mann-Kendall method against R's "mk.test()" function from the "trend" package with data values pulled from a known example from "The Gamma Family..."
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Pohlert T (2023). trend: Non-Parametric Trend Tests and Change-Point Detection. R package version 1.1.6, https://CRAN.R-project.org/package=trend
        /// </remarks>
        [TestMethod]
        public void Test_MannKendall()
        {
            // Table 1.2 Maximum annual peak discharge values in cms, observed at the Harricana River at Amos (Quebec, Canada)
            var data = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

            var p = HypothesisTests.MannKendallTest(data);
            double true_p = 0.7757;
            Assert.AreEqual(true_p, p, 1E-4);
        }

        /// <summary>
        /// Test the linear trend test method against the p value of R's "lm()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_LinearTrendTest()
        {
            var time = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100 };
            var data = new double[100];
            for (int i = 0; i < time.Length; i++)
            {
                data[i] = Math.Cos(time[i]);
            }

            var pVal = HypothesisTests.LinearTrendTest(time, data);
            double truePVal = 0.9092;
            Assert.AreEqual(truePVal, pVal, 1E-4);
        }

        /// <summary>
        /// Test the Grubbs-Beck method against a known example from "The Gamma Family..."
        /// </summary>
        [TestMethod]
        public void Test_GrubbsBeck()
        {
            // Table 1.2 Maximum annual peak discharge values in cms, observed at the Harricana River at Amos (Quebec, Canada)
            var data = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

            double xHi = 0, xLo = 0;
            MultipleGrubbsBeckTest.GrubbsBeckTest(data, out xHi, out xLo);

            double true_xHi = 378.2, true_xLo = 91.2; // page 9 of reference.
            Assert.AreEqual(true_xHi, xHi, 1E-1);
            Assert.AreEqual(true_xLo, xLo, 1E-1);
        }

        /// <summary>
        /// Test the Multiple Grubbs Beck method against the HEC-SSP software
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// put link
        /// </remarks>
        [TestMethod]
        public void Test_MultipleGrubbsBeck()
        {
            var sample = new double[] { 423.823157583651d, 769.119600544858d, 1420.9840065658d, 81.1048635557593d, 279.928123967548d, 75.85549285788d, 545.403510679849d, 2765.04138262183d, 99.2328081106953d, 1151.90527161336d, 35.0188163971524d, 93.2163892297505d, 174.065209255604d, 284.811439281534d, 69.4129231117978d, 1393.70251526941d, 366.909211754559d, 57.3577922448949d, 507.512883027978d, 3408.66910217811d, 994.625641160531d, 99.0457917640901d, 253.32702656322d, 155.691675526921d, 644.834237627036d, 81.0277133561875d, 655.861119414072d, 49.0010266733043d, 216.450613452697d, 625.639165872462d };
            var BLU = new double[] { 12903d, 10108d, 7401.3d, 7233.3d, 7167.3d, 7116.3d, 6930.3d, 6929d, 6870d, 6768d, 6742.7d, 6213.3d, 6166.3d, 6071.7d, 6044.3d, 5857d, 5779d, 5289.7d, 5247.3d, 5243d, 5208.3d, 5050.3d, 4887.7d, 4821.3d, 4801d, 4789.7d, 4713d, 4426.7d, 4308.7d, 4265d, 4263.3d, 4199d, 4138.7d, 4093.3d, 4035d, 4023d, 3964d, 3929.7d, 3826.3d, 3698.3d, 3667.7d, 3654d, 3647.3d, 3644d, 3578d, 3532d, 3511.3d, 3510.7d, 3463.3d, 3374.3d, 3349d, 3312.3d, 3263d, 3226d, 3126.7d, 3119.3d, 3084d, 3081d, 3054.3d, 3045.7d, 2968.7d, 2964d, 2907.3d, 2889.3d, 2850.3d, 2759.7d, 2663d, 2643d, 2631.3d, 2596.3d, 2580.3d, 2433.7d, 2295.3d, 2095.3d, 2024d, 1963d, 1938.7d, 1785.3d, 1266.7d, 1142.3d, 74.7d };
            var CGR = new double[] { 20964.3d, 14314d, 10881.7d, 10714.7d, 10053.3d, 9861.7d, 9406.7d, 8858d, 8716.7d, 8589.7d, 8521.7d, 8038.7d, 8026.7d, 7863d, 7709.7d, 7501d, 7353.3d, 7106.7d, 6559d, 6494.7d, 6469d, 6390d, 6193.7d, 6098.7d, 6014.3d, 5963.3d, 5911.7d, 5845d, 5697d, 5551.3d, 5500d, 5242.7d, 5175.7d, 5063.3d, 5050d, 4961.7d, 4944.7d, 4906.3d, 4856d, 4777d, 4696d, 4581.7d, 4567d, 4550.3d, 4513.3d, 4512.7d, 4440.3d, 4411.3d, 4390.7d, 4386.7d, 4295.7d, 4114.7d, 4089d, 3915d, 3899.7d, 3891.3d, 3880d, 3727.3d, 3714.7d, 3702d, 3693d, 3649.3d, 3573.3d, 3505d, 3484.7d, 3472.3d, 3356.7d, 3220.7d, 3214.7d, 3103.3d, 3095.3d, 3071d, 2693.3d, 2675.3d, 2579d, 2263d, 2203.7d, 2070d, 2050.7d, 1423d, 302d };
            var FCK = new double[] { 11433.3d, 10686.7d, 9512d, 9155d, 8964.7d, 8536.7d, 8320d, 8150.3d, 7971d, 7866.7d, 7713.3d, 7643.3d, 7363.3d, 7250d, 6930d, 6683d, 6565.3d, 6466.7d, 6403.3d, 6171d, 5951d, 5854.3d, 5757.3d, 5600d, 5559d, 5559d, 5553.3d, 5475.7d, 5467.7d, 5390.7d, 5390d, 5307.3d, 5215d, 5164.7d, 5051.7d, 5026.7d, 4946.7d, 4841.3d, 4808d, 4790d, 4643.3d, 4477.3d, 4459.7d, 4423.3d, 4369.3d, 4350d, 4283.3d, 4220d, 4136.7d, 4123.3d, 4050d, 4029.7d, 3986.3d, 3984.7d, 3980.7d, 3957.7d, 3878.7d, 3872.7d, 3821.3d, 3763.3d, 3737d, 3522d, 3444.3d, 3440d, 3426.7d, 3417.7d, 3406.3d, 3304.7d, 3197.7d, 3192.7d, 3122d, 2991.3d, 2640d, 2472.3d, 2213.3d, 2133.3d, 2011d, 1950.7d, 1796.3d, 1610.3d, 70.7d };
            var FOS = new double[] { 59843d, 53545.7d, 37129.3d, 36316d, 36183.7d, 36039d, 35913.7d, 35367.7d, 33473.7d, 33335.7d, 31654.3d, 31066.3d, 30669.7d, 30644.3d, 30130d, 29868d, 29333d, 29263.3d, 29012d, 28973d, 28671.7d, 25682.3d, 25485d, 24577.3d, 24187.7d, 23488d, 22258.3d, 21918.7d, 21845.7d, 21721.7d, 21715d, 20967.7d, 20843d, 20561d, 20554.3d, 20454.3d, 20437.3d, 20267.3d, 20044.3d, 19911.3d, 19735.7d, 19602.3d, 19214d, 18623.3d, 18370.7d, 18360.7d, 18093.7d, 17882.3d, 17839.7d, 17510d, 17439d, 17382.3d, 17194.3d, 16199d, 16168.3d, 16162d, 16133.3d, 15987.7d, 15380d, 15360.7d, 15062d, 15045d, 14833d, 14803.7d, 14505.3d, 14358.7d, 14235.3d, 14031d, 13545d, 13192.7d, 12821.3d, 12741.3d, 12402d, 11895d, 11526.3d, 11110d, 9099.7d, 9017.7d, 6200.7d, 5860d, 382d };
            var GPV = new double[] { 38666.7d, 35789d, 25168.7d, 23866.7d, 23633.3d, 23000d, 22427.7d, 21300d, 21086.3d, 20466.7d, 20300d, 20094d, 19333.3d, 18733.3d, 18553d, 18433.3d, 18033d, 17838.3d, 17723.3d, 17625.7d, 17485d, 17433.3d, 16438d, 15725.7d, 15656.7d, 15266.7d, 15233.3d, 14833.3d, 14751.7d, 14333.3d, 14058.7d, 13824.3d, 13666.7d, 13582d, 13115.7d, 12883.3d, 12831d, 12774d, 12719.3d, 12453d, 12379.3d, 12303.3d, 12165.7d, 11966.3d, 11746.7d, 11666.7d, 11510d, 11383.7d, 11259.3d, 10788d, 10756.7d, 10739.3d, 10433.3d, 10341.7d, 10258.7d, 10230d, 10210d, 9990d, 9819.3d, 9758d, 9658d, 9367.3d, 9333.3d, 9230d, 9068.3d, 8941d, 8696.7d, 8593.3d, 8233.3d, 8136.7d, 8114d, 8109.3d, 7563.3d, 7333.3d, 7333.3d, 7065.7d, 5873.3d, 5636.7d, 4033.3d, 3988.3d, 231d };
            var HCK = new double[] { 23520d, 16220d, 15290d, 14750d, 13980d, 13740d, 13030d, 12930d, 12550d, 12140d, 11470d, 10610d, 10220d, 10060d, 10010d, 9880d, 9090d, 8780d, 8700d, 8420d, 8410d, 8280d, 8140d, 7990d, 7770d, 7540d, 7460d, 7400d, 7400d, 6770d, 6590d, 6500d, 6420d, 6200d, 6190d, 6160d, 6100d, 5900d, 5800d, 5730d, 5690d, 5670d, 5630d, 5620d, 5590d, 5580d, 5560d, 5440d, 5410d, 5340d, 5330d, 5310d, 5150d, 5080d, 4690d, 4680d, 4620d, 4600d, 4580d, 4490d, 4380d, 4110d, 4080d, 4050d, 4020d, 4010d, 3980d, 3960d, 3900d, 3840d, 3770d, 3700d, 3670d, 3620d, 3300d, 3220d, 2820d, 2810d, 2670d, 2440d, 2370d, 2030d, 1750d, 1570d, 920d };
            var LOP = new double[] { 60350d, 44630d, 38800d, 34640d, 34510d, 34050d, 33270d, 32570d, 29750d, 28080d, 27860d, 27210d, 26930d, 26920d, 26580d, 24320d, 24260d, 24060d, 23630d, 23530d, 23100d, 21970d, 21930d, 21660d, 21270d, 21260d, 20520d, 19370d, 19370d, 18800d, 17970d, 17630d, 17420d, 16240d, 16060d, 16010d, 15980d, 15970d, 15850d, 15810d, 15700d, 15640d, 15420d, 15150d, 15100d, 15090d, 14890d, 14820d, 14580d, 14130d, 14120d, 14090d, 13460d, 13350d, 13350d, 13330d, 13180d, 12930d, 12740d, 12610d, 12580d, 12350d, 11900d, 11650d, 11420d, 11390d, 11250d, 11120d, 11050d, 10580d, 10180d, 9890d, 9710d, 9130d, 8780d, 8390d, 8180d, 8110d, 6980d, 6230d, 6150d, 6000d, 4630d, 4560d, 2840d };

            int BLU_True = 3;
            int CGR_True = 2;
            int FCK_True = 9;
            int FOS_True = 3;
            int GPV_True = 3;
            int HCK_True = 1;
            int LOP_True = 1;
            int Samp_True = 0;
            int LO;
            LO = MultipleGrubbsBeckTest.Function(BLU);
            Assert.AreEqual(BLU_True, LO);
            LO = MultipleGrubbsBeckTest.Function(CGR);
            Assert.AreEqual(CGR_True, LO);
            LO = MultipleGrubbsBeckTest.Function(FCK);
            Assert.AreEqual(FCK_True, LO);
            LO = MultipleGrubbsBeckTest.Function(FOS);
            Assert.AreEqual(FOS_True, LO);
            LO = MultipleGrubbsBeckTest.Function(GPV);
            Assert.AreEqual(GPV_True, LO);
            LO = MultipleGrubbsBeckTest.Function(HCK);
            Assert.AreEqual(HCK_True, LO);
            LO = MultipleGrubbsBeckTest.Function(LOP);
            Assert.AreEqual(LOP_True, LO);
            LO = MultipleGrubbsBeckTest.Function(sample);
            Assert.AreEqual(Samp_True, LO);
        }


        /// <summary>
        /// GMM test for unimodality is compared against the R package 'mclust'
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Scrucca, L., Fop, M., Murphy, T. B., &amp; Raftery, A. E. (2016). mclust 5: Clustering, Classification and Density Estimation Using Gaussian Finite Mixture Models. The R Journal, 8(1), 289–317. https://doi.org/10.32614/RJ-2016-021
        /// </remarks>
        [TestMethod]
        public void Test_UnimodalityTest()
        {
            var unimodalData = new double[] { 4.5, 5.2, 5.1, 4.9, 5.0, 5.3, 5.4, 4.8, 4.7, 5.2,
                   5.1, 4.6, 5.0, 5.3, 5.1, 4.9, 5.2, 5.0, 4.8, 5.3,
                   4.9, 5.1, 5.2, 4.8, 5.0, 5.1, 5.2, 4.7, 5.3, 5.0 };

            var pval = HypothesisTests.UnimodalityTest(unimodalData);

            double true_pval = 0.4142441;
            Assert.AreEqual(true_pval, pval, 1E-4);

            var bimodalData = new double[] { 3.8, 3.9, 4.0, 3.9, 4.0, 4.1, 4.0, 3.8, 3.9, 4.0,
                4.1, 4.0, 3.9, 4.0, 4.0, 4.3, 4.4, 4.4, 4.5, 4.4, 4.3,
                4.4, 4.5, 4.4, 4.3, 4.4, 4.4, 4.5, 4.4, 4.3 };

            pval = HypothesisTests.UnimodalityTest(bimodalData);

            true_pval = 2.55425752131444e-05;
            Assert.AreEqual(true_pval, pval, 1E-9);

        }

        /// <summary>
        /// The Mann-Kendall variance correction uses the size of each tie group, so a heavily tied
        /// record matches the reference implementation.
        /// </summary>
        /// <remarks>
        /// Reference p-value computed with Python pymannkendall 1.4.3 original_test and confirmed
        /// by an independent evaluation of the Kendall (1975) variance formula with
        /// scipy 1.17.1: p = 2.687544409241127E-08 for this fixture of twelve distinct values
        /// followed by twelve tied values.
        /// </remarks>
        [TestMethod]
        public void Test_MannKendall_HeavyTies()
        {
            var sample = new double[] { 8.017, 8.219, 8.759, 9.016, 9.091, 9.452, 10.002, 10.12, 10.597, 10.714, 10.98, 12.68,
                20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0 };
            double pValue = HypothesisTests.MannKendallTest(sample);
            Assert.AreEqual(2.687544409241127E-08, pValue, 2.687544409241127E-08 * 1E-6);
        }

        /// <summary>
        /// The Mann-Whitney variance correction uses the size of each tie group, so tied samples
        /// match the scipy asymptotic test without continuity correction.
        /// </summary>
        /// <remarks>
        /// Reference p-value computed with Python scipy 1.17.1
        /// scipy.stats.mannwhitneyu(x, y, use_continuity=False, alternative='two-sided',
        /// method='asymptotic') = 0.03509475123312452.
        /// </remarks>
        [TestMethod]
        public void Test_MannWhitney_HeavyTies()
        {
            var sample1 = new double[] { 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 5.0, 5.0 };
            var sample2 = new double[] { 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 5.0, 5.0, 6.0, 6.0, 7.0, 8.0 };
            double pValue = HypothesisTests.MannWhitneyTest(sample1, sample2);
            Assert.AreEqual(0.03509475123312452, pValue, 0.03509475123312452 * 1E-6);
        }

        /// <summary>
        /// With n1 n2 odd the expected U is a half-integer (5 x 17 samples: E[U] = 42.5).
        /// </summary>
        /// <remarks>
        /// Reference p-value computed with Python scipy 1.17.1
        /// scipy.stats.mannwhitneyu(x, y, use_continuity=False, method='asymptotic') = 0.00687292391731027.
        /// </remarks>
        [TestMethod]
        public void Test_MannWhitney_OddProduct()
        {
            var sample1 = new double[] { 1.1, 2.3, 3.2, 4.8, 5.0 };
            var sample2 = new double[] { 2.0, 3.5, 4.1, 5.5, 6.2, 7.7, 8.1, 9.4, 10.3, 11.0, 12.5, 13.1, 14.2, 15.6, 16.0, 17.3, 18.9 };
            double pValue = HypothesisTests.MannWhitneyTest(sample1, sample2);
            Assert.AreEqual(0.00687292391731027, pValue, 0.00687292391731027 * 1E-6);
        }

        /// <summary>
        /// Test Hartigan's dip statistic against R diptest 0.77.2 <c>dip()</c>, including tied, constant,
        /// single-value and empty samples.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Maechler, M. and Ringach, D. (2025). diptest: Hartigan's Dip Test Statistic for Unimodality - Corrected. R package version 0.77-2.
        /// </remarks>
        [TestMethod]
        public void Test_HartiganDip_MatchesR()
        {
            Assert.AreEqual(0d, HypothesisTests.HartiganDip(new double[0]), 0d);
            Assert.AreEqual(0.5, HypothesisTests.HartiganDip(new double[] { 5 }), 0d);
            Assert.AreEqual(0.25, HypothesisTests.HartiganDip(new double[] { 1, 2 }), 0d);
            Assert.AreEqual(0.125, HypothesisTests.HartiganDip(new double[] { 4, 1, 3, 2 }), 0d);
            Assert.AreEqual(0.1, HypothesisTests.HartiganDip(new double[] { 3, 3, 3, 3, 3 }), 1E-15);
            Assert.AreEqual(0.1875, HypothesisTests.HartiganDip(new double[] { 1, 1, 1, 2, 2, 5, 5, 5, 5, 9 }), 1E-15);
            Assert.AreEqual(0.15555555555555556, HypothesisTests.HartiganDip(new double[] { 0.5, 1.5, 1.5, 2, 2, 2, 2.5, 3, 7, 7.5, 8, 8, 8.5, 9, 9 }), 1E-15);
            Assert.AreEqual(0.125, HypothesisTests.HartiganDip(new double[] { 1, 0.5, -1, 0.5, -0.5, -1, 1, 1, -1, 0.5, 2.5, 1, -0.5, 0.5, -1.5, 0.5, -2, 2.5, -1, 1 }), 1E-15);
            Assert.AreEqual(0.086314958296960093, HypothesisTests.HartiganDip(new double[] { -0.2600038160874425, -0.16285100167438976, -0.90087458650539542, -0.23135597425869309,
                0.86195506130218913, 0.3763467681724586, -0.17715607751019694, -0.085177164177269976, -0.43612827351821809, 1.4237823458859504, 3.416678901465473,
                3.8604758569957012, 3.6775664316148138, 3.0153586596756199, 2.0927864607817059, 1.2240458252098512, 3.6884804177102644, 2.4817555475005166,
                1.5391226578204371, 3.8898272439175168 }), 1E-15);
            Assert.AreEqual(0.10710382513661201, HypothesisTests.HartiganDip(BimodalDipFixture), 1E-15);
        }

        /// <summary>
        /// A 30-value bimodal sample: R diptest gives D = 0.10710382513661201, p = 0.0047111 from the
        /// interpolated table and p = 0.0048 from 100,000 simulated uniform samples.
        /// </summary>
        private static readonly double[] BimodalDipFixture =
        {
            -2.1, -1.7, -1.5, -1.2, -1.1, -0.9, -0.6, -0.4, -0.3, -0.1, 0.2, 0.4, 0.6, 0.9, 1.3,
            4.1, 4.4, 4.6, 4.8, 5.0, 5.1, 5.3, 5.6, 5.9, 6.2, 6.5, 6.8, 7.0, 7.4, 8.0
        };

        /// <summary>
        /// Test that the dip statistic rejects NaN, infinite and null input.
        /// </summary>
        [TestMethod]
        public void Test_HartiganDip_InvalidInput()
        {
            Assert.ThrowsExactly<ArgumentException>(() => HypothesisTests.HartiganDip(new double[] { 1, double.NaN, 3, 4 }));
            Assert.ThrowsExactly<ArgumentException>(() => HypothesisTests.HartiganDip(new double[] { 1, double.PositiveInfinity, 3, 4 }));
            Assert.ThrowsExactly<ArgumentNullException>(() => HypothesisTests.HartiganDip(null!));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => HypothesisTests.HartiganDipTest(BimodalDipFixture, 0));
        }

        /// <summary>
        /// Test the simulated dip test p-value against R diptest 0.77.2 <c>dip.test()</c> within Monte
        /// Carlo error, its reproducibility for a fixed seed, and the small-sample and minimum-dip cases.
        /// </summary>
        [TestMethod]
        public void Test_HartiganDipTest_MatchesR()
        {
            // R: p = 0.0048 (100,000 simulations). Default 2,000 realizations: standard error about 0.0015.
            double p = HypothesisTests.HartiganDipTest(BimodalDipFixture);
            Assert.AreEqual(0.0048, p, 0.006);
            Assert.AreEqual(p, HypothesisTests.HartiganDipTest(BimodalDipFixture), 0d);

            // Normal quantiles have the minimum dip 1/(2n); R gives p = 1
            var quantiles = new double[40];
            for (int i = 0; i < 40; i++)
                quantiles[i] = Normal.StandardZ((i + 0.5) / 40);
            Assert.AreEqual(1d / 80, HypothesisTests.HartiganDip(quantiles), 1E-15);
            Assert.AreEqual(1d, HypothesisTests.HartiganDipTest(quantiles), 0d);

            // Three or fewer values and constant samples cannot reject
            Assert.AreEqual(1d, HypothesisTests.HartiganDipTest(new double[] { 1, 2, 5 }), 0d);
            Assert.AreEqual(1d, HypothesisTests.HartiganDipTest(new double[] { 3, 3, 3, 3, 3, 3 }), 0d);
        }

        /// <summary>
        /// Test the parametric-bootstrap p-value of the Gaussian mixture likelihood-ratio test: it is
        /// reproducible for a fixed seed, large for the unimodal fixture and at its minimum 1/(B + 1)
        /// for the bimodal fixture, and the chi-square overload returns its reference value. This test uses
        /// B = 39; the sweep with B = 199 runs with
        /// <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep Unimodality</c>.
        /// </summary>
        [TestMethod]
        public void Test_UnimodalityTest_Bootstrap()
        {
            UnimodalityTestBootstrap(39);
        }

        /// <summary>
        /// The two-sample Kolmogorov-Smirnov statistic equals scipy.stats.ks_2samp(first, second).statistic
        /// bit for bit, including tied values within and across the samples. The expected values are the
        /// scipy 1.17.1 statistics, which its exact mode rounds from the integer count, each the correctly
        /// rounded rational noted beside it.
        /// </summary>
        [TestMethod]
        public void Test_TwoSampleKolmogorovSmirnov_ScipyOracle()
        {
            Assert.AreEqual(1d / 3d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 1d, 2d, 3d }, new[] { 2d, 3d, 4d }), 0d);
            Assert.AreEqual(1d / 3d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 1d, 1d, 2d, 2d, 3d }, new[] { 1d, 2d, 2d, 2d, 5d, 6d }), 0d);
            Assert.AreEqual(0d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 0.5d }, new[] { 0.5d }), 0d);
            Assert.AreEqual(1d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 1d, 2d, 3d, 4d, 5d, 6d, 7d }, new[] { 10d, 11d }), 0d);
            // 8/35
            Assert.AreEqual(0.22857142857142856d, HypothesisTests.TwoSampleKolmogorovSmirnov(
                new[] { -3.2d, -1.0d, 0.0d, 0.4d, 2.2d, 2.2d, 7.5d }, new[] { -1.0d, 0.4d, 0.4d, 1.9d, 3.3d }), 0d);

            // a = sorted((7i mod 23)/4, i < 40), b = sorted((5i mod 17)/3 + 0.5, i < 31): many ties. 77/620
            var a = new double[40];
            var b = new double[31];
            for (int i = 0; i < a.Length; i++) a[i] = (i * 7 % 23) / 4.0;
            for (int i = 0; i < b.Length; i++) b[i] = (i * 5 % 17) / 3.0 + 0.5;
            Array.Sort(a);
            Array.Sort(b);
            Assert.AreEqual(0.12419354838709677d, HypothesisTests.TwoSampleKolmogorovSmirnov(a, b), 0d);
            // The statistic is symmetric, and a list argument gives the same value as an array.
            Assert.AreEqual(HypothesisTests.TwoSampleKolmogorovSmirnov(a, b), HypothesisTests.TwoSampleKolmogorovSmirnov(b, a), 0d);
            Assert.AreEqual(HypothesisTests.TwoSampleKolmogorovSmirnov(a, b),
                HypothesisTests.TwoSampleKolmogorovSmirnov(new System.Collections.Generic.List<double>(a), new System.Collections.Generic.List<double>(b)), 0d);

            // a = sorted((p1 i mod q1)/4, i < n1), b = sorted((p2 i mod q2)/4 + s/4, i < n2): quarter-step ties.
            AssertModularCase(40, 7, 23, 31, 5, 17, 2, 0.2d);                   // 1/5
            AssertModularCase(13, 3, 11, 29, 7, 13, 1, 0.21750663129973474d);   // 82/377
            AssertModularCase(97, 11, 31, 50, 13, 29, 3, 0.1134020618556701d);  // 11/97
            AssertModularCase(200, 17, 41, 7, 3, 5, 0, 0.875d);                 // 7/8
            AssertModularCase(64, 5, 19, 64, 9, 19, 1, 0.09375d);               // 3/32
            AssertModularCase(1, 1, 1, 9, 2, 7, 0, 0.7777777777777778d);        // 7/9
        }

        /// <summary>
        /// Builds the sorted quarter-step modular samples of <see cref="Test_TwoSampleKolmogorovSmirnov_ScipyOracle"/>
        /// and checks the statistic bit for bit.
        /// </summary>
        /// <param name="n1">The first sample size.</param>
        /// <param name="p1">The first multiplier.</param>
        /// <param name="q1">The first modulus.</param>
        /// <param name="n2">The second sample size.</param>
        /// <param name="p2">The second multiplier.</param>
        /// <param name="q2">The second modulus.</param>
        /// <param name="shift">The second sample shift, in quarters.</param>
        /// <param name="expected">The scipy statistic.</param>
        private static void AssertModularCase(int n1, int p1, int q1, int n2, int p2, int q2, int shift, double expected)
        {
            var a = new double[n1];
            var b = new double[n2];
            for (int i = 0; i < n1; i++) a[i] = (i * p1 % q1) / 4.0;
            for (int i = 0; i < n2; i++) b[i] = (i * p2 % q2) / 4.0 + shift / 4.0;
            Array.Sort(a);
            Array.Sort(b);
            Assert.AreEqual(expected, HypothesisTests.TwoSampleKolmogorovSmirnov(a, b), 0d);
        }

        /// <summary>
        /// The statistic is the correctly rounded rational D/(n1 n2). Tie-heavy samples give the exact
        /// quotient, including 1/3 for [0, 1, 1] against [1, 1, 1] and for [1, 2, 3] against [2, 3, 4]; in the
        /// second case the floating-point difference 1 - 2/3 is one unit in the last place above 1/3.
        /// </summary>
        [TestMethod]
        public void Test_TwoSampleKolmogorovSmirnov_ExactRationals()
        {
            Assert.AreEqual(1d / 3d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 0d, 1d, 1d }, new[] { 1d, 1d, 1d }), 0d);
            Assert.AreEqual(1d / 3d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 1d, 2d, 3d }, new[] { 2d, 3d, 4d }), 0d);
            // Unequal sizes with ties: F1 - F2 is 2/3 - 1/7 = 11/21 at 1, then 1 - 4/7 = 3/7 at 2.
            Assert.AreEqual(11d / 21d, HypothesisTests.TwoSampleKolmogorovSmirnov(
                new[] { 1d, 1d, 2d }, new[] { 1d, 2d, 2d, 2d, 3d, 3d, 3d }), 0d);
            // Disjoint samples, including samples tied at a single value each.
            Assert.AreEqual(1d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 2d, 2d, 2d, 2d }, new[] { 5d, 5d, 5d }), 0d);
            Assert.AreEqual(1d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 4d, 5d, 6d }, new[] { 1d, 2d }), 0d);
            // Identical samples, and samples with the same empirical distribution at different sizes.
            var same = new[] { -1d, 0d, 0d, 3d, 3d, 3d, 8d };
            Assert.AreEqual(0d, HypothesisTests.TwoSampleKolmogorovSmirnov(same, (double[])same.Clone()), 0d);
            Assert.AreEqual(0d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 1d, 2d }, new[] { 1d, 1d, 2d, 2d }), 0d);
            // One value shared by both samples: both distributions jump together and the gap is 0.
            Assert.AreEqual(0d, HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 7d, 7d }, new[] { 7d, 7d, 7d }), 0d);
            // A single value against a sample containing it: 7 of the 10 values lie strictly below it.
            Assert.AreEqual(0.7d, HypothesisTests.TwoSampleKolmogorovSmirnov(
                new[] { 0.5d }, new[] { 0.1d, 0.1d, 0.2d, 0.2d, 0.3d, 0.3d, 0.4d, 0.5d, 0.6d, 0.7d }), 0d);
        }

        /// <summary>
        /// On fixed-seed random samples with heavy ties, the statistic equals the correctly rounded
        /// rational computed independently by counting each sample at every pooled value.
        /// </summary>
        [TestMethod]
        public void Test_TwoSampleKolmogorovSmirnov_RandomizedExactRational()
        {
            var prng = new MersenneTwister(20261107);
            for (int trial = 0; trial < 500; trial++)
            {
                int n1 = 1 + prng.Next(60), n2 = 1 + prng.Next(60), levels = 1 + prng.Next(12);
                var a = new double[n1];
                var b = new double[n2];
                for (int i = 0; i < n1; i++) a[i] = prng.Next(levels + 1) / 4.0;
                for (int i = 0; i < n2; i++) b[i] = prng.Next(levels + 1) / 4.0;
                Array.Sort(a);
                Array.Sort(b);

                long numerator = 0;
                foreach (double value in a)
                    numerator = Math.Max(numerator, CountGap(a, b, value));
                foreach (double value in b)
                    numerator = Math.Max(numerator, CountGap(a, b, value));
                // Both integers are exact doubles, so IEEE division gives the correctly rounded quotient.
                double expected = (double)numerator / ((double)n1 * n2);
                Assert.AreEqual(expected, HypothesisTests.TwoSampleKolmogorovSmirnov(a, b), 0d);
            }
        }

        /// <summary>
        /// Returns |i n2 - j n1|, where i and j count the values of each sample at or below the given value.
        /// </summary>
        /// <param name="a">The first sample.</param>
        /// <param name="b">The second sample.</param>
        /// <param name="value">The evaluation point.</param>
        /// <returns>The scaled gap between the two empirical distribution functions.</returns>
        private static long CountGap(double[] a, double[] b, double value)
        {
            long i = 0, j = 0;
            foreach (double v in a) if (v <= value) i++;
            foreach (double v in b) if (v <= value) j++;
            return Math.Abs(i * b.Length - j * a.Length);
        }

        /// <summary>
        /// Null, empty, non-finite and unsorted samples throw the documented exceptions naming the argument.
        /// </summary>
        [TestMethod]
        public void Test_TwoSampleKolmogorovSmirnov_InvalidInputs()
        {
            var good = new[] { 1d, 2d };
            Assert.AreEqual("first", Assert.Throws<ArgumentNullException>(() => HypothesisTests.TwoSampleKolmogorovSmirnov(null!, good)).ParamName);
            Assert.AreEqual("second", Assert.Throws<ArgumentNullException>(() => HypothesisTests.TwoSampleKolmogorovSmirnov(good, null!)).ParamName);
            Assert.AreEqual("first", Assert.Throws<ArgumentException>(() => HypothesisTests.TwoSampleKolmogorovSmirnov(new double[0], good)).ParamName);
            Assert.AreEqual("second", Assert.Throws<ArgumentException>(() => HypothesisTests.TwoSampleKolmogorovSmirnov(good, new[] { 2d, 1d })).ParamName);
            Assert.AreEqual("first", Assert.Throws<ArgumentOutOfRangeException>(() => HypothesisTests.TwoSampleKolmogorovSmirnov(new[] { 1d, double.NaN }, good)).ParamName);
            Assert.AreEqual("second", Assert.Throws<ArgumentOutOfRangeException>(() => HypothesisTests.TwoSampleKolmogorovSmirnov(good, new[] { 1d, double.PositiveInfinity })).ParamName);
        }
    }

}
