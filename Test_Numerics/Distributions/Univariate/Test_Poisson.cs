using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Poisson distribution algorithm.
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
    public class Test_Poisson
    {

        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_PoissonDist()
        {
            double true_mean = 4.2d;
            double true_mode = 4.0d;
            double true_median = 4.0d;
            double true_stdDev = Math.Sqrt(4.2d);
            double true_skew = 0.488d;
            double true_kurt = 3.2381d;
            double true_pdf = 0.19442365170822165d;
            double true_cdf = 0.58982702131057763d;
            double true_icdf05 = 1.0d;
            double true_icdf95 = 8.0d;
            var P = new Poisson(4.2d);
            Assert.AreEqual(true_mean, P.Mean, 0.0001d);
            Assert.AreEqual(true_median, P.Median, 0.0001d);
            Assert.AreEqual(true_mode, P.Mode, 0.0001d);
            Assert.AreEqual(true_stdDev, P.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_skew, P.Skewness, 0.0001d);
            Assert.AreEqual(true_kurt, P.Kurtosis, 0.0001d);
            Assert.AreEqual(true_pdf, P.PDF(4.0d), 0.0001d);
            Assert.AreEqual(true_cdf, P.CDF(4.0d), 0.0001d);
            Assert.AreEqual(true_icdf05, P.InverseCDF(0.05d), 0.0001d);
            Assert.AreEqual(true_icdf95, P.InverseCDF(0.95d), 0.0001d);
        }

        /// <summary>
        /// Verifying that input parameters can create this distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var P = new Poisson();
            Assert.AreEqual(1, P.Lambda);

            var P2 = new Poisson(10);
            Assert.AreEqual(10,P2.Lambda);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var P = new Poisson(double.NaN);
            Assert.IsFalse(P.ParametersValid);

            var P2 = new Poisson(double.PositiveInfinity);
            Assert.IsFalse(P2.ParametersValid);

            var P3 = new Poisson(0);
            Assert.IsFalse(P3.ParametersValid);
        }

        /// <summary>
        /// Testing ParameterToString().
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var P = new Poisson();
            Assert.AreEqual("Rate (λ)",P.ParametersToString[0, 0]);
            Assert.AreEqual("1",P.ParametersToString[0, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Poisson(10);
            var mom = dist.CentralMoments(1000);
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
            var P = new Poisson();
            Assert.AreEqual(1,P.Mean);

            var P2 = new Poisson(10);
            Assert.AreEqual(10,P2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var P = new Poisson();
            Assert.AreEqual(1, P.Median, 1E-4);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var P = new Poisson();
            Assert.AreEqual(1, P.Mode);

            var P2 = new Poisson(2.4);
            Assert.AreEqual(2, P2.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var P = new Poisson();
            Assert.AreEqual(1, P.StandardDeviation);

            var P2 = new Poisson(4);
            Assert.AreEqual(2, P2.StandardDeviation);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var P = new Poisson();
            Assert.AreEqual(1, P.Skewness);

            var P2 = new Poisson(4);
            Assert.AreEqual(0.5, P2.Skewness);
        }

        /// <summary>
        /// Testing kurtosis of this distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var P = new Poisson();
            Assert.AreEqual(4, P.Kurtosis);

            var P2 = new Poisson(4);
            Assert.AreEqual(3.25, P2.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var P = new Poisson();
            Assert.AreEqual(0, P.Minimum);
            Assert.AreEqual(double.PositiveInfinity,P.Maximum);

            var P2 = new Poisson(4);
            Assert.AreEqual(0, P2.Minimum);
            Assert.AreEqual(double.PositiveInfinity, P2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var P = new Poisson(1.5);
            Assert.AreEqual(0.33469, P.PDF(1),  1e-04);
            Assert.AreEqual(0.00000354, P.PDF(10),  1e-08);

            var P2 = new Poisson(5.4);
            Assert.AreEqual(0.024389, P2.PDF(1),  1e-05);
            Assert.AreEqual(0.02624, P2.PDF(10),  1e-05);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var P = new Poisson(1.5);
            Assert.AreEqual(0.55782, P.CDF(1),  1e-04);
            Assert.AreEqual(0.999999, P.CDF(10), 1e-06);

            var P2 = new Poisson(10.8);
            Assert.AreEqual(0.00024, P2.CDF(1),  1e-05);
            Assert.AreEqual(0.483969, P2.CDF(10),  1e-05);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var P = new Poisson();
            Assert.AreEqual(0, P.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity,P.InverseCDF(1) );
        }

        /// <summary>
        /// Below λ = 100,000 the PMF, CDF and quantile at integer k are bit-identical to the established
        /// expressions, which are inlined here.
        /// </summary>
        [TestMethod()]
        public void Test_OrdinaryRates_BitIdentical()
        {
            double[] rates = { 0.001, 0.1, 0.5, 1, 3.5, 4.2, 10, 30.7, 100, 999.9, 5000, 25000, 99999.99 };
            double[] probabilities = { 1E-12, 1E-6, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999999 };
            foreach (double lambda in rates)
            {
                var P = new Poisson(lambda);
                double kMax = lambda + 10d * Math.Sqrt(lambda) + 20d;
                int step = Math.Max(1, (int)(Math.Sqrt(lambda) / 10d));
                for (int k = 0; k <= kMax; k += step)
                {
                    double pdf = Math.Exp(-lambda + k * Math.Log(lambda) - Factorial.LogFactorial(k));
                    Assert.AreEqual(pdf, P.PDF(k), 0d);
                    Assert.AreEqual(Gamma.UpperIncomplete(k + 1d, lambda), P.CDF(k), 0d);
                }
                foreach (double p in probabilities)
                {
                    int k = (int)Math.Max(0, Math.Floor(lambda));
                    while (k > 0 && Gamma.UpperIncomplete(k, lambda) >= p) k--;
                    while (Gamma.UpperIncomplete(k + 1d, lambda) < p) k++;
                    Assert.AreEqual(k, P.InverseCDF(p), 0d);
                }
            }
        }

        /// <summary>
        /// The PMF is zero away from the integers and at +∞; the CDF is one at +∞.
        /// </summary>
        [TestMethod()]
        public void Test_EdgeArguments()
        {
            var P = new Poisson(3.5);
            // scipy.stats.poisson.pmf(2.5, 3.5) = 0
            Assert.AreEqual(0d, P.PDF(2.5));
            Assert.AreEqual(double.NegativeInfinity, P.LogPDF(2.5));
            Assert.AreEqual(P.CDF(2), P.CDF(2.5), 0d);
            Assert.AreEqual(0d, P.PDF(double.PositiveInfinity));
            Assert.AreEqual(1d, P.CDF(double.PositiveInfinity));
            Assert.AreEqual(0d, P.CCDF(double.PositiveInfinity));
            Assert.AreEqual(0d, P.PDF(-1));
            Assert.AreEqual(0d, P.CDF(-0.5));
            Assert.AreEqual(1d, P.CCDF(-0.5));
            Assert.IsTrue(double.IsNaN(P.PDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(P.LogPDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(P.CDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(P.CCDF(double.NaN)));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => P.InverseCDF(double.NaN));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => P.InverseCDF(-0.1));
        }

        /// <summary>
        /// The log PMF and the survival function stay finite where the PMF and 1 − CDF underflow.
        /// </summary>
        [TestMethod()]
        public void Test_LogPDF_CCDF_LogCCDF()
        {
            var P = new Poisson(3.5);
            // scipy.stats.poisson.logpmf(1000, 3.5)
            Assert.AreEqual(-4662.865209992796, P.LogPDF(1000), 4662.865209992796 * 1E-14);
            Assert.AreEqual(Math.Log(P.PDF(4)), P.LogPDF(4), 1E-15);
            // scipy.stats.poisson.sf(30, 3.5)
            Assert.AreEqual(3.028074579820694E-19, P.CCDF(30), 3.028074579820694E-19 * 1E-12);
            // scipy.stats.poisson.logsf(30, 3.5)
            Assert.AreEqual(-42.641189801555306, P.LogCCDF(30), 42.641189801555306 * 1E-13);
            // mpmath.log(mpmath.gammainc(201, 0, 3.5, regularized=True))
            Assert.AreEqual(-620.2124583665938, P.LogCCDF(200), 620.2124583665938 * 1E-13);
            // scipy.stats.poisson.sf(4, 4.2)
            Assert.AreEqual(0.41017297868942243, new Poisson(4.2).CCDF(4), 0.41017297868942243 * 1E-13);
            // scipy.stats.poisson.sf(0, 0.5)
            Assert.AreEqual(0.3934693402873665, new Poisson(0.5).CCDF(0.7), 0.3934693402873665 * 1E-14);
        }

        /// <summary>
        /// From λ = 100,000 the CDF uses the log gamma tail, which keeps near full precision.
        /// </summary>
        [TestMethod()]
        public void Test_CDF_LargeRate()
        {
            // mpmath.gammainc(k + 1, lambda, mpmath.inf, regularized=True), 60 digits
            (double lambda, double k, double expected)[] cases =
            {
                (1E5, 1E5, 0.5008410430993401),
                (1E6, 1E6, 0.5002659614862837),
                (1E7, 1E7, 0.500084104416326),
                (1E8, 1E8, 0.5000265961519927),
                (1E5, 99000, 0.0007742008294447389),
                (1E5, 101500, 0.9999988978481349),
                (1E6, 1003000, 0.9986464098708987),
                (1E7, 9985000, 1.0456666897991976E-06),
                (3E9, 3E9, 0.5000048557708024),
                (3E9, 3E9 + 2E5, 0.999869620764168),
            };
            foreach (var (lambda, k, expected) in cases)
                Assert.AreEqual(expected, new Poisson(lambda).CDF(k), expected * 1E-12);
        }

        /// <summary>
        /// From λ = 100,000 the mass uses the gamma deviance form and matches high-precision values.
        /// </summary>
        [TestMethod()]
        public void Test_PDF_LargeRate()
        {
            // mpmath (50 digits): exp(-lam + k*log(lam) - loggamma(k + 1))
            var cases = new (double lambda, double k, double expected)[]
            {
                (1E5, 1E5, 0.0012615652097053006),
                (1E5, 99000d, 8.4012719339368129E-6),
                (1E6, 1E6, 0.00039894224715624403),
                (1E7, 10003000d, 8.0432571048805763E-5),
                (1E8, 1E8, 3.9894228006898078E-5),
            };
            foreach (var (lambda, k, expected) in cases)
            {
                var P = new Poisson(lambda);
                Assert.AreEqual(expected, P.PDF(k), expected * 1E-12);
                Assert.AreEqual(Math.Log(expected), P.LogPDF(k), 1E-12);
            }
        }

        /// <summary>
        /// Rates beyond the range of <see cref="int"/> give finite PMF values and terminating quantiles.
        /// </summary>
        [TestMethod()]
        public void Test_HugeRate()
        {
            var P = new Poisson(3E9);
            // mpmath: exp(-lam + k*log(lam) - loggamma(k + 1)), lam = k = 3e9 (R dpois(3e9, 3e9) = 7.2837e-06)
            Assert.AreEqual(7.28365620374487E-06, P.PDF(3E9), 7.28365620374487E-06 * 1E-12);
            // mpmath: same expression at k = 3e9 + 5e4
            Assert.AreEqual(4.801653207235253E-06, P.PDF(3E9 + 5E4), 4.801653207235253E-06 * 1E-12);
            Assert.AreEqual(0d, P.PDF(3E9 + 0.5));
            // scipy.stats.poisson.ppf(p, 3e9)
            Assert.AreEqual(2999651582d, P.InverseCDF(1E-10), 0d);
            Assert.AreEqual(2999909908d, P.InverseCDF(0.05), 0d);
            Assert.AreEqual(3000000000d, P.InverseCDF(0.5), 0d);
            Assert.AreEqual(3000169260d, P.InverseCDF(0.999), 0d);
            Assert.AreEqual(3E9, P.Median, 0d);
        }

        /// <summary>
        /// From λ = 100,000 the quantile is found by bracketing and bisection and matches scipy.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_LargeRate()
        {
            // scipy.stats.poisson.ppf(p, lambda)
            var P5 = new Poisson(1E5);
            Assert.AreEqual(97995d, P5.InverseCDF(1E-10), 0d);
            Assert.AreEqual(99480d, P5.InverseCDF(0.05), 0d);
            Assert.AreEqual(100000d, P5.InverseCDF(0.5), 0d);
            Assert.AreEqual(100979d, P5.InverseCDF(0.999), 0d);
            var P7 = new Poisson(1E7);
            Assert.AreEqual(9979890d, P7.InverseCDF(1E-10), 0d);
            Assert.AreEqual(9994799d, P7.InverseCDF(0.05), 0d);
            Assert.AreEqual(10000000d, P7.InverseCDF(0.5), 0d);
            Assert.AreEqual(10009774d, P7.InverseCDF(0.999), 0d);
            // The result is the smallest k with CDF(k) >= p.
            foreach (double p in new[] { 1E-300, 1E-12, 0.3, 0.7, 1 - 1E-12 })
            {
                double k = P5.InverseCDF(p);
                Assert.IsGreaterThanOrEqualTo(p, P5.CDF(k));
                Assert.IsTrue(k == 0 || P5.CDF(k - 1) < p);
            }
        }

        /// <summary>
        /// Maximum likelihood and the method of moments both give λ = mean(x).
        /// </summary>
        [TestMethod()]
        public void Test_Estimation()
        {
            var data = new double[] { 3, 1, 4, 1, 5, 9, 2, 6, 5, 3 };
            var P = new Poisson();
            P.Estimate(data, ParameterEstimationMethod.MaximumLikelihood);
            // Poisson MLE: lambda = mean(data) = 3.9
            Assert.AreEqual(3.9, P.Lambda, 1E-15);
            P.Estimate(data, ParameterEstimationMethod.MethodOfMoments);
            Assert.AreEqual(3.9, P.Lambda, 1E-15);
            Assert.AreEqual(3.9, P.ParametersFromMoments(new[] { 3.9, 2.0 })[0]);
            var moments = P.MomentsFromParameters(new[] { 4.0 });
            Assert.AreEqual(4.0, moments[0]);
            Assert.AreEqual(2.0, moments[1]);
            Assert.AreEqual(0.5, moments[2]);
            Assert.AreEqual(3.25, moments[3]);

            var constraints = P.GetParameterConstraints(data);
            Assert.AreEqual(3.9, constraints.Item1[0], 1E-15);
            Assert.IsTrue(constraints.Item2[0] < 3.9 && 3.9 < constraints.Item3[0]);

            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => P.Estimate(new double[] { 1, 2.5 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => P.Estimate(new double[] { 1, -1 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => P.Estimate(new double[] { 0, 0 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => P.Estimate(new double[0], ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentNullException>(() => P.Estimate(null!, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<NotImplementedException>(() => P.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments));
        }

        /// <summary>
        /// The parametric bootstrap re-estimates a valid distribution from a seeded sample.
        /// </summary>
        [TestMethod()]
        public void Test_Bootstrap()
        {
            var P = new Poisson(6);
            var b = (Poisson)P.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 500, 12345);
            Assert.IsTrue(b.ParametersValid);
            Assert.AreEqual(6d, b.Lambda, 0.5);
        }

        /// <summary>
        /// The quantile function validates the rate before returning the endpoints.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_ValidatesParametersFirst()
        {
            var P = new Poisson(2);
            P.Lambda = -1;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => P.InverseCDF(0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => P.InverseCDF(1));
        }

        /// <summary>
        /// The per-element static quantile against scipy.stats.poisson.ppf (scipy 1.18.1), which is exact in integers.
        /// </summary>
        [TestMethod()]
        public void Test_StaticInverseCDF_Scipy()
        {
            double[] rates = [0.1, 1, 3.7, 25, 150, 1000, 2E5];
            double[] probabilities = [1E-10, 0.01, 0.3, 0.5, 0.9, 0.999999];
            double[] expected =
            [
                0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 9, 0, 0, 3, 4, 6, 16, 1, 14, 22, 25, 32, 52,
                79, 122, 143, 150, 166, 212, 806, 927, 983, 1000, 1041, 1154, 197162, 198960, 199765, 200000, 200573, 202129
            ];
            var r = new double[expected.Length];
            var p = new double[expected.Length];
            for (int i = 0; i < rates.Length; i++)
                for (int j = 0; j < probabilities.Length; j++)
                {
                    r[i * probabilities.Length + j] = rates[i];
                    p[i * probabilities.Length + j] = probabilities[j];
                }
            double[] actual = Poisson.InverseCDF(r, p);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i], 0d, "rate " + r[i] + ", p " + p[i]);
            CollectionAssert.AreEqual(new[] { 0d, double.PositiveInfinity }, Poisson.InverseCDF([2d, 2d], [0d, 1d]));
        }

        /// <summary>
        /// The per-element static log mass and CDF against exact values (mpmath 1.3.0, 60 digits).
        /// </summary>
        /// <remarks>
        /// Below rate 1E5 the log mass is -λ + k ln λ - ln k!, whose rounding is bounded by the size of the cancelling
        /// terms, so it is checked to 2^-52 (λ + k |ln λ| + ln k!) absolute; the error achieved is 1.1E-12 at rate 1000
        /// and k = 980, where the terms reach 1.4E4. The CDF is checked to 1E-13 relative. scipy.stats.poisson 1.18.1
        /// gives -4.56414668187972 at rate 1000, k = 980 and -7.647703528404236 at rate 2E5, k = 200500.
        /// </remarks>
        [TestMethod()]
        public void Test_StaticLogPDF_CDF_Exact()
        {
            double[] rates = [0.5, 3.7, 3.7, 25, 1000, 2E5];
            double[] k = [0, 2, 7, 30, 980, 200500];
            double[] logPmf = [-0.5, -1.7764815412595878, -3.0668316235141626, -3.091961602784142, -4.564146681878576, -7.647703528430046];
            double[] cdf = [0.6065306597126334, 0.28543311310006825, 0.964759269527861, 0.8633088691526638, 0.26981377056126365, 0.8684423557001304];
            double[] actualLog = Poisson.LogPDF(rates, k);
            double[] actualCdf = Poisson.CDF(rates, k);
            for (int i = 0; i < k.Length; i++)
            {
                double terms = rates[i] + k[i] * Math.Abs(Math.Log(rates[i])) + Factorial.LogFactorial((int)k[i]);
                Assert.AreEqual(logPmf[i], actualLog[i], Math.Pow(2, -52) * terms, "logpmf " + i);
                Assert.AreEqual(cdf[i], actualCdf[i], 1E-13 * cdf[i], "cdf " + i);
            }
        }

        /// <summary>
        /// The static members are bitwise identical to the instance members on a grid of rates and arguments,
        /// including non-integer, negative, infinite and not-a-number arguments and the large-rate branch.
        /// </summary>
        [TestMethod()]
        public void Test_StaticMembers_BitwiseInstance()
        {
            double[] rates = [1E-3, 0.1, 0.7, 1, 2.5, 9.99, 37, 400.5, 9999, 99999.5, 1E5, 3.3E5, 1E9];
            double[] points = [double.NaN, double.NegativeInfinity, -1, 0, 0.5, 1, 2, 3.9, 7, 15, 60, 401, 10050, 99000, 1E5 + 300, 3.3E5, 1E9 + 1E4, 3E9, double.PositiveInfinity];
            double[] probabilities = [0, 9.999999999999999E-301, 1E-300, 1E-12, 0.001, 0.25, 0.5, 0.75, 0.999, 1 - 1E-10,
                0.9999999999000001, 1 - 1E-12, 0.9999999999999999, 1];
            int n = rates.Length;
            var r = new double[n * points.Length];
            var x = new double[r.Length];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < points.Length; j++) { r[i * points.Length + j] = rates[i]; x[i * points.Length + j] = points[j]; }
            double[] log = Poisson.LogPDF(r, x), cdf = Poisson.CDF(r, x);
            for (int i = 0; i < r.Length; i++)
            {
                var d = new Poisson(r[i]);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(d.LogPDF(x[i])), BitConverter.DoubleToInt64Bits(log[i]), "LogPDF " + r[i] + ", " + x[i]);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(d.CDF(x[i])), BitConverter.DoubleToInt64Bits(cdf[i]), "CDF " + r[i] + ", " + x[i]);
            }
            var rq = new double[n * probabilities.Length];
            var pq = new double[rq.Length];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < probabilities.Length; j++) { rq[i * probabilities.Length + j] = rates[i]; pq[i * probabilities.Length + j] = probabilities[j]; }
            double[] q = Poisson.InverseCDF(rq, pq);
            for (int i = 0; i < rq.Length; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(new Poisson(rq[i]).InverseCDF(pq[i])), BitConverter.DoubleToInt64Bits(q[i]), "InverseCDF " + rq[i] + ", " + pq[i]);
        }

        /// <summary>
        /// On arrays long enough to be evaluated in parallel ranges, the static members are bitwise identical to the
        /// instance members, and the quantile searched from its Cornish-Fisher start equals the instance search from
        /// the floor of the rate.
        /// </summary>
        [TestMethod()]
        public void Test_StaticMembers_ParallelBitwiseInstance()
        {
            const int n = 10_000;
            var rng = new Random(12345);
            var r = new double[n];
            var p = new double[n];
            var x = new double[n];
            for (int i = 0; i < n; i++)
            {
                r[i] = i % 4 == 0 ? Math.Exp(Math.Log(1E-3) + Math.Log(2E8) * rng.NextDouble()) : 0.5 + 49.5 * rng.NextDouble();
                p[i] = (i % 10) switch
                {
                    8 => Math.Pow(10, -320 * rng.NextDouble()),
                    9 => 1 - Math.Pow(10, -16 * rng.NextDouble()),
                    _ => rng.NextDouble()
                };
                x[i] = Math.Floor(2 * r[i] * rng.NextDouble());
            }
            double[] q = Poisson.InverseCDF(r, p), log = Poisson.LogPDF(r, x), cdf = Poisson.CDF(r, x);
            for (int i = 0; i < n; i++)
            {
                var d = new Poisson(r[i]);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(d.InverseCDF(p[i])), BitConverter.DoubleToInt64Bits(q[i]), "InverseCDF " + r[i] + ", " + p[i]);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(d.LogPDF(x[i])), BitConverter.DoubleToInt64Bits(log[i]), "LogPDF " + r[i] + ", " + x[i]);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(d.CDF(x[i])), BitConverter.DoubleToInt64Bits(cdf[i]), "CDF " + r[i] + ", " + x[i]);
            }
        }

        /// <summary>
        /// The static members reject null and mismatched arrays, invalid rates and invalid probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_StaticMembers_InvalidArguments()
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => Poisson.LogPDF(null!, [1d]));
            Assert.ThrowsExactly<ArgumentNullException>(() => Poisson.CDF([1d], null!));
            Assert.ThrowsExactly<ArgumentException>(() => Poisson.InverseCDF([1d, 2d], [0.5]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Poisson.LogPDF([0d], [1d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Poisson.CDF([double.NaN], [1d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Poisson.InverseCDF([double.PositiveInfinity], [0.5]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Poisson.InverseCDF([1d], [1.5]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Poisson.InverseCDF([1d], [double.NaN]));
            Assert.IsEmpty(Poisson.InverseCDF([], []));
        }
    }
}
