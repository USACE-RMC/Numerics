using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Binomial distribution algorithm.
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
    public class Test_Binomial
    {
        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_BinomialDist()
        {
            double true_mean = 1.92d;
            double true_mode = 2.0d;
            double true_median = 2.0d;
            double true_stdDev = Math.Sqrt(1.6896d);
            double true_skew = 0.5847d;
            double true_kurt = 3.2169d;
            double true_pdf = 0.28218979948821621d;
            double true_cdf = 0.12933699143209909d;
            double true_icdf05 = 0.0d;
            double true_icdf95 = 4.0d;
            var B = new Binomial(0.12d, 16);
            Assert.AreEqual(B.Mean, true_mean, 0.0001d);
            Assert.AreEqual(B.Median, true_median, 0.0001d);
            Assert.AreEqual(B.Mode, true_mode, 0.0001d);
            Assert.AreEqual(B.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(B.Skewness, true_skew, 0.0001d);
            Assert.AreEqual(B.Kurtosis, true_kurt, 0.0001d);
            Assert.AreEqual(B.PDF(1.0d), true_pdf, 0.0001d);
            Assert.AreEqual(B.CDF(0.0d), true_cdf, 0.0001d);
            Assert.AreEqual(B.InverseCDF(0.05d), true_icdf05, 0.0001d);
            Assert.AreEqual(B.InverseCDF(0.95d), true_icdf95, 0.0001d);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. See if binomial is being created.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var B = new Binomial(0, 4);
            Assert.AreEqual(0, B.ProbabilityOfSuccess);
            Assert.AreEqual(4,B.NumberOfTrials);

            var B2 = new Binomial(0.3, 3);
            Assert.AreEqual(0.3, B2.ProbabilityOfSuccess);
            Assert.AreEqual(3,B2.NumberOfTrials);

            var B3 = new Binomial(1, 2);
            Assert.AreEqual(1, B3.ProbabilityOfSuccess);
            Assert.AreEqual(2,B3.NumberOfTrials);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. See what parameters fail.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var b = new Binomial(double.NaN, 1); 
            Assert.IsFalse(b.ParametersValid);

            var b2 = new Binomial(-1, 1);
            Assert.IsFalse(b2.ParametersValid);

            var b3 = new Binomial(2, 1);
            Assert.IsFalse(b3.ParametersValid);

            var b4 = new Binomial(0.3, -1);
            Assert.IsFalse(b4.ParametersValid);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking string output.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
            var b = new Binomial(0.5, 4);
            Assert.AreEqual("Probability of Success (p)", b.ParametersToString[0, 0]);
            Assert.AreEqual("Number of Trials (n)", b.ParametersToString[1, 0]);
            Assert.AreEqual("0.5", b.ParametersToString[0, 1]);
            Assert.AreEqual("4", b.ParametersToString[1, 1]);
        }

        /// <summary>
        ///  Verified using MathNet-Numerics testing. Checking mean of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var b = new Binomial(0.5, 4);
            Assert.AreEqual(2, b.Mean);

            var b2 = new Binomial(0, 4);
            Assert.AreEqual(0, b2.Mean);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking median of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var b = new Binomial(0.5, 4);
            Assert.AreEqual(2, b.Mean);

            var b2 = new Binomial(0, 4);
            Assert.AreEqual(0, b2.Mean);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking mode of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var b = new Binomial(0, 4);
            Assert.AreEqual(0, b.Mode);

            var b2 = new Binomial(0.3, 3);
            Assert.AreEqual(1, b2.Mode);

            var b3 = new Binomial(1, 2);
            Assert.AreEqual(2, b3.Mode);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking standard deviation with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var b = new Binomial(0, 4);
            Assert.AreEqual(0, b.StandardDeviation);

            var b2 = new Binomial(0.3, 3);
            Assert.AreEqual(0.793725, b2.StandardDeviation, 1e-04);

            var b3 = new Binomial(1, 2);
            Assert.AreEqual(0, b3.StandardDeviation);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking skewness of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var b = new Binomial(0.5, 5);
            Assert.AreEqual(0, b.Skewness);

            var b2 = new Binomial(0.3, 3);
            Assert.AreEqual(0.503952, b2.Skewness, 1e-04);
        }

        /// <summary>
        /// Checking Kurtosis algorithm.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var b = new Binomial(0.5, 5);
            Assert.AreEqual(2.6, b.Kurtosis);

            var b2 = new Binomial(0.3, 3);
            Assert.AreEqual(2.587301, b2.Kurtosis, 1e-04);
        }

        /// <summary>
        /// Checking minimum function
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var b = new Binomial(0.3, 3);
            Assert.AreEqual(0, b.Minimum);
        }

        /// <summary>
        /// Checking maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var b = new Binomial(0.3, 3);
            Assert.AreEqual(3, b.Maximum);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing PDF function.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var b = new Binomial(0, 1);
            Assert.AreEqual(1, b.PDF(0));
            Assert.AreEqual(0, b.PDF(1));

            var b2 = new Binomial(0.3, 1);
            Assert.AreEqual(0.69999, b2.PDF(0), 1e-04);
            Assert.AreEqual(0.29999, b2.PDF(1), 1e-04);

            var b3 = new Binomial(0.3, 3);
            Assert.AreEqual(0.34299999, b3.PDF(0), 1e-04);
            Assert.AreEqual(0.440999, b3.PDF(1), 1e-04);

            var b4 = new Binomial(0.3, 10);
            Assert.AreEqual(0.0282475, b4.PDF(0), 1e-04);
            Assert.AreEqual(0.1210608, b4.PDF(1), 1e-04);
        }

        /// <summary>
        /// Testing CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var b = new Binomial(0.4, 8);
            Assert.AreEqual(0.5940864, b.CDF(3),1e-04);
        }

        /// <summary>
        /// Testing inverse CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var b = new Binomial(0.3, 100);
            Assert.AreEqual(32, b.InverseCDF(0.7));
        }

        /// <summary>
        /// Verify Binomial InverseCDF correctly returns NumberOfTrials at high probabilities.
        /// Reference: scipy.stats.binom.ppf(p, n, p_success)
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_Boundary()
        {
            // binom(10, 0.5).ppf(0.95) = 8
            var b1 = new Binomial(0.5, 10);
            Assert.AreEqual(8, b1.InverseCDF(0.95));
            // binom(10, 0.5).ppf(0.99) = 9
            Assert.AreEqual(9, b1.InverseCDF(0.99));
            // binom(10, 0.5).ppf(1.0) = 10 (must reach NumberOfTrials)
            Assert.AreEqual(10, b1.InverseCDF(1.0));

            // binom(5, 0.3).ppf(1.0) = 5
            var b2 = new Binomial(0.3, 5);
            Assert.AreEqual(5, b2.InverseCDF(1.0));
        }

        /// <summary>
        /// Median is the smallest k with CDF(k) of at least one half. Reference values from scipy.stats.binom.median.
        /// </summary>
        [TestMethod()]
        public void Test_Median_SmallestKWithCdfAtLeastHalf()
        {
            Assert.AreEqual(0, new Binomial(0.1, 3).Median);
            Assert.AreEqual(5, new Binomial(0.5, 10).Median);
            Assert.AreEqual(2, new Binomial(0.12, 16).Median);
            Assert.AreEqual(3, new Binomial(0.3, 10).Median);
        }

        /// <summary>
        /// Median over a grid of trials and success probabilities, including p = 1/2 with odd n, where CDF((n - 1)/2) is
        /// exactly 1/2. Reference values are the smallest k with CDF(k) of at least 1/2 computed in exact rational
        /// arithmetic (Python fractions) at each double-precision p; R qbinom(0.5, n, p) agrees on every listed case.
        /// </summary>
        [TestMethod()]
        public void Test_Median_ExactGrid()
        {
            var trials = new[] { 1, 2, 3, 5, 7, 10, 11, 16, 30, 35, 39, 99, 100, 101, 1000, 1001 };
            var probabilities = new[] { 0.01, 0.1, 0.25, 1d / 3d, 0.45, 0.5, 0.55, 0.9, 0.99 };
            var medians = new int[,]
            {
                { 0, 0, 0, 0, 0, 0, 1, 1, 1 },  // n = 1
                { 0, 0, 0, 1, 1, 1, 1, 2, 2 },  // n = 2
                { 0, 0, 1, 1, 1, 1, 2, 3, 3 },  // n = 3
                { 0, 0, 1, 2, 2, 2, 3, 5, 5 },  // n = 5
                { 0, 1, 2, 2, 3, 3, 4, 6, 7 },  // n = 7
                { 0, 1, 2, 3, 4, 5, 6, 9, 10 },  // n = 10
                { 0, 1, 3, 4, 5, 5, 6, 10, 11 },  // n = 11
                { 0, 1, 4, 5, 7, 8, 9, 15, 16 },  // n = 16
                { 0, 3, 7, 10, 13, 15, 17, 27, 30 },  // n = 30
                { 0, 3, 9, 12, 16, 17, 19, 32, 35 },  // n = 35
                { 0, 4, 10, 13, 18, 19, 21, 35, 39 },  // n = 39
                { 1, 10, 25, 33, 45, 49, 54, 89, 98 },  // n = 99
                { 1, 10, 25, 33, 45, 50, 55, 90, 99 },  // n = 100
                { 1, 10, 25, 34, 45, 50, 56, 91, 100 },  // n = 101
                { 10, 100, 250, 333, 450, 500, 550, 900, 990 },  // n = 1000
                { 10, 100, 250, 334, 450, 500, 551, 901, 991 },  // n = 1001
            };
            for (int i = 0; i < trials.Length; i++)
            {
                for (int j = 0; j < probabilities.Length; j++)
                {
                    Assert.AreEqual(medians[i, j], new Binomial(probabilities[j], trials[i]).Median, $"n = {trials[i]}, p = {probabilities[j]}");
                }
            }
        }

        /// <summary>
        /// The mass stays finite for large numbers of trials. Reference values from scipy.stats.binom.pmf and scipy.stats.binom.logpmf.
        /// </summary>
        [TestMethod()]
        public void Test_PDF_LargeNumberOfTrials()
        {
            var b = new Binomial(0.5, 1100);
            Assert.AreEqual(0.02405165776823181, b.PDF(550), 1e-12);
            Assert.AreEqual(-3.72755135473426, b.LogPDF(550), 1e-11);
            Assert.AreEqual(-762.4618986159398, b.LogPDF(0), 1e-9);
            Assert.AreEqual(-762.4618986159398, b.LogPDF(1100), 1e-9);
        }

        /// <summary>
        /// LogPDF agrees with the log of PDF and handles the support boundaries. Reference values from scipy.stats.binom.pmf and scipy.stats.binom.logpmf.
        /// </summary>
        [TestMethod()]
        public void Test_LogPDF()
        {
            var b = new Binomial(0.3, 10);
            Assert.AreEqual(0.2668279319999998, b.PDF(3), 1e-14);
            Assert.AreEqual(-1.321151277766889, b.LogPDF(3), 1e-12);
            Assert.AreEqual(double.NegativeInfinity, b.LogPDF(-1));
            Assert.AreEqual(double.NegativeInfinity, b.LogPDF(11));
            Assert.AreEqual(0, b.PDF(11));
            for (int k = 0; k <= 10; k++)
                Assert.AreEqual(Math.Log(b.PDF(k)), b.LogPDF(k), 1e-12);

            var b0 = new Binomial(0, 4);
            Assert.AreEqual(0, b0.LogPDF(0));
            Assert.AreEqual(double.NegativeInfinity, b0.LogPDF(1));
            var b1 = new Binomial(1, 4);
            Assert.AreEqual(0, b1.LogPDF(4));
            Assert.AreEqual(double.NegativeInfinity, b1.LogPDF(3));
        }

        /// <summary>
        /// InverseCDF rejects a NaN probability.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_NaNProbabilityThrows()
        {
            var b = new Binomial(0.3, 10);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => b.InverseCDF(double.NaN));
        }

        /// <summary>
        /// The PMF is zero away from the integers, the CDF floors its argument, and integer values are bit-identical
        /// to the established log-space expression, which is inlined here.
        /// </summary>
        [TestMethod()]
        public void Test_PDF_NonInteger()
        {
            var b = new Binomial(0.3, 10);
            // scipy.stats.binom.pmf(2.5, 10, 0.3) = 0
            Assert.AreEqual(0d, b.PDF(2.5));
            Assert.AreEqual(double.NegativeInfinity, b.LogPDF(2.5));
            Assert.AreEqual(0d, b.PDF(-0.5));
            Assert.AreEqual(b.CDF(2), b.CDF(2.5), 0d);
            Assert.AreEqual(0d, b.PDF(double.PositiveInfinity));
            for (int k = 0; k <= 10; k++)
            {
                double logChoose = Factorial.LogFactorial(10) - Factorial.LogFactorial(k) - Factorial.LogFactorial(10 - k);
                double log = logChoose + k * Math.Log(0.3) + (10d - k) * Math.Log(1d - 0.3);
                Assert.AreEqual(log, b.LogPDF(k), 0d);
                Assert.AreEqual(Math.Exp(log), b.PDF(k), 0d);
            }
        }

        /// <summary>
        /// The per-element static log mass and CDF against exact sums (mpmath 1.3.0, 50 digits).
        /// </summary>
        /// <remarks>The log mass forms ln C(n, k) as ln n! - ln k! - ln (n - k)!, so its absolute error grows with
        /// ln n!: 5.8E-13 at n = 500, p = 0.01, k = 2 and 1.3E-11 at n = 20000, p = 0.001, k = 15. It is checked to
        /// 1E-15 ln n! absolute (at least 1E-14). The CDF error also grows with n (8E-13 relative at n = 20000) and is
        /// checked to 1E-16 n relative (at least 1E-13). scipy.stats.binom 1.18.1 differs
        /// from these references by up to 1.6E-11 in the log mass and 1.7E-14 relative in the CDF.</remarks>
        [TestMethod()]
        public void Test_StaticLogPDF_CDF_Exact()
        {
            double[] p = [0.3, 0.5, 0.01, 0.9, 0.123, 0.6, 0.6, 1E-3];
            int[] n = [10, 1, 500, 40, 1000, 25, 25, 20000];
            double[] k = [3, 1, 2, 35, 100, 0, 25, 15];
            double[] logPmf = [-1.3211512777668886, -0.6931471805599453, -2.4813406134061355, -1.8035711447058922,
                -5.767262876665994, -22.907268296853875, -12.770640594149768, -2.9635377164075183];
            double[] cdf = [0.6496107184000001, 1.0, 0.12338577435358183, 0.37098230346579864,
                0.013360206641529525, 1.1258999068426256E-10, 1.0, 0.15638396191370524];
            double[] actualLog = Binomial.LogPDF(p, n, k);
            double[] actualCdf = Binomial.CDF(p, n, k);
            for (int i = 0; i < k.Length; i++)
            {
                Assert.AreEqual(logPmf[i], actualLog[i], Math.Max(1E-14, 1E-15 * Factorial.LogFactorial(n[i])), "logpmf " + i);
                Assert.AreEqual(cdf[i], actualCdf[i], Math.Max(1E-13, 1E-16 * n[i]) * cdf[i], "cdf " + i);
            }
        }

        /// <summary>
        /// The static members are bitwise identical to the instance members on a grid of parameters and arguments,
        /// including the degenerate success probabilities and non-integer, out-of-range and infinite arguments.
        /// </summary>
        [TestMethod()]
        public void Test_StaticMembers_BitwiseInstance()
        {
            double[] probabilities = [0, 1E-9, 0.01, 0.25, 0.5, 0.77, 0.999, 1];
            int[] trials = [1, 2, 7, 30, 171, 1000, 25000];
            double[] points = [double.NegativeInfinity, -1, 0, 0.5, 1, 2, 6.5, 7, 29, 30, 170, 500, 24999, 25000, 1E6, double.PositiveInfinity];
            int count = probabilities.Length * trials.Length * points.Length;
            var p = new double[count];
            var n = new int[count];
            var k = new double[count];
            int index = 0;
            foreach (double pi in probabilities)
                foreach (int ni in trials)
                    foreach (double ki in points) { p[index] = pi; n[index] = ni; k[index] = ki; index++; }
            double[] log = Binomial.LogPDF(p, n, k), cdf = Binomial.CDF(p, n, k);
            for (int i = 0; i < count; i++)
            {
                var b = new Binomial(p[i], n[i]);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(b.LogPDF(k[i])), BitConverter.DoubleToInt64Bits(log[i]), "LogPDF " + i);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(b.CDF(k[i])), BitConverter.DoubleToInt64Bits(cdf[i]), "CDF " + i);
            }
        }

        /// <summary>
        /// On arrays long enough to be evaluated in parallel ranges, the static members are bitwise identical to the
        /// instance members.
        /// </summary>
        [TestMethod()]
        public void Test_StaticMembers_ParallelBitwiseInstance()
        {
            const int count = 10_000;
            var rng = new Random(54321);
            var p = new double[count];
            var n = new int[count];
            var k = new double[count];
            for (int i = 0; i < count; i++)
            {
                p[i] = i % 50 == 0 ? (i % 100 == 0 ? 0d : 1d) : rng.NextDouble();
                n[i] = i % 4 == 0 ? 1 + rng.Next(100_000) : 1 + rng.Next(200);
                k[i] = i % 9 == 0 ? rng.Next(-2, n[i] + 3) + 0.5 : rng.Next(-1, n[i] + 2);
            }
            double[] log = Binomial.LogPDF(p, n, k), cdf = Binomial.CDF(p, n, k);
            for (int i = 0; i < count; i++)
            {
                var b = new Binomial(p[i], n[i]);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(b.LogPDF(k[i])), BitConverter.DoubleToInt64Bits(log[i]), "LogPDF " + i);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(b.CDF(k[i])), BitConverter.DoubleToInt64Bits(cdf[i]), "CDF " + i);
            }
        }

        /// <summary>
        /// The static members reject null and mismatched arrays and invalid parameters.
        /// </summary>
        [TestMethod()]
        public void Test_StaticMembers_InvalidArguments()
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => Binomial.LogPDF(null!, [1], [0d]));
            Assert.ThrowsExactly<ArgumentNullException>(() => Binomial.CDF([0.5], null!, [0d]));
            Assert.ThrowsExactly<ArgumentNullException>(() => Binomial.CDF([0.5], [1], null!));
            Assert.ThrowsExactly<ArgumentException>(() => Binomial.LogPDF([0.5, 0.5], [1, 1], [0d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Binomial.LogPDF([1.5], [1], [0d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Binomial.CDF([double.NaN], [1], [0d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Binomial.CDF([0.5], [0], [0d]));
            Assert.IsEmpty(Binomial.CDF([], [], []));
        }
    }
}
