using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Inverse Chi-Squared probability distribution algorithm.
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
    public class Test_InverseChiSquared
    {
        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_InverseChiSquaredDist()
        {
            // Reference: scipy.special.gammaincc(v/2, v*sigma/(2*x)) for CDF
            // InverseCDF: v*sigma / (2 * gammainccinv(v/2, p))
            double true_mean = 0.2;
            double true_median = 0.15758426609d;
            double true_pdf = 0.0000063457380298844403d;
            double true_cdf = 0.9999884277d;
            double true_icdf = 6.27d;
            var IX = new InverseChiSquared(7, (1d / 7d));
            double pdf = IX.PDF(6.27d);
            double cdf = IX.CDF(6.27d);
            double icdf = IX.InverseCDF(cdf);
            Assert.AreEqual(IX.Mean, true_mean, 0.0001d);
            Assert.AreEqual(IX.Median, true_median, 0.0001d);
            Assert.AreEqual(IX.PDF(6.27d), true_pdf, 0.0001d);
            Assert.AreEqual(IX.CDF(6.27d), true_cdf, 0.0001d);
            Assert.AreEqual(IX.InverseCDF(IX.CDF(6.27d)), true_icdf, 0.001d);
        }

        /// <summary>
        /// Checking Inverse Chi-Squared can be created with inputs.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(10, IX.DegreesOfFreedom);
            Assert.AreEqual(1, IX.Sigma);

            var IX2 = new InverseChiSquared(2, 1);
            Assert.AreEqual(2, IX2.DegreesOfFreedom);
            Assert.AreEqual(1, IX2.Sigma);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var IX = new InverseChiSquared(0, 0);
            Assert.IsFalse(IX.ParametersValid);

            var IX2 = new InverseChiSquared(1, double.NaN);
            Assert.IsFalse(IX2.ParametersValid);

            var IX3 = new InverseChiSquared(1, double.PositiveInfinity);
            Assert.IsFalse(IX3.ParametersValid);
        }

        /// <summary>
        /// Testing Parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual("Degrees of Freedom (ν)", IX.ParametersToString[0, 0]);
            Assert.AreEqual("Scale (σ)", IX.ParametersToString[1, 0]);
            Assert.AreEqual("10", IX.ParametersToString[0, 1]);
            Assert.AreEqual("1", IX.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(1.25, IX.Mean);

            var IX2 = new InverseChiSquared(2, 2);
            Assert.AreEqual(double.NaN, IX2.Mean);
        }

        /// <summary>
        /// Testing median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            // Reference: v*sigma / (2 * gammainccinv(v/2, 0.5))
            var IX = new InverseChiSquared();
            Assert.AreEqual(1.07046, IX.Median, 1e-04);

            var IX2 = new InverseChiSquared(7, 1);
            Assert.AreEqual(1.10309, IX2.Median, 1e-04);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(0.8333, IX.Mode,  1e-04);

            var IX2 = new InverseChiSquared(2, 2);
            Assert.AreEqual(1, IX2.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(0.72168, IX.StandardDeviation,  1e-04);

            var IX2 = new InverseChiSquared(2, 2);
            Assert.AreEqual(double.NaN, IX2.StandardDeviation);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(3.46410, IX.Skewness,  1e-04);

            var IX2 = new InverseChiSquared(2, 2);
            Assert.AreEqual(double.NaN, IX2.Skewness);
        }

        /// <summary>
        /// Testing kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(45, IX.Kurtosis);

            var IX2 = new InverseChiSquared(2,2);
            Assert.AreEqual(double.NaN, IX2.Kurtosis);
        }

        /// <summary>
        /// Testing Minimum and Maximum functions are 0 and positive infinity respectively.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(0, IX.Minimum);
            Assert.AreEqual(double.PositiveInfinity, IX.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var IX = new InverseChiSquared(1, 1);
            Assert.AreEqual(0.2419, IX.PDF(1), 1e-04);

            var IX2 = new InverseChiSquared(2, 1);
            Assert.AreEqual(0.15163, IX2.PDF(2),  1e-04);

            var IX3 = new InverseChiSquared();
            Assert.AreEqual(0.16700, IX3.PDF(2),  1e-04);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            // Reference: scipy.special.gammaincc(3.5, 3.5/5) = 0.985571264449
            var IX = new InverseChiSquared(7,1);
            Assert.AreEqual(0.985571264449d, IX.CDF(5), 1e-06);
        }

        /// <summary>
        /// Testing the Inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(0, IX.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity,IX.InverseCDF(1));
            // Reference: scipy.stats.chi2.isf => v*sigma / chi2.isf(p, v)
            Assert.AreEqual(0.84884, IX.InverseCDF(0.3), 1e-04);
        }

        /// <summary>
        /// Verify CDF and InverseCDF are correct and mutually inverse.
        /// Reference: scipy.stats.chi2.sf(v*sigma/x, v) for CDF; v*sigma/chi2.isf(p, v) for InverseCDF.
        /// </summary>
        [TestMethod]
        public void Test_CDF_InverseCDF_Roundtrip()
        {
            // v=10, sigma=1: CDF values from scipy.stats.chi2.sf(v*sigma/x, v)
            var dist = new InverseChiSquared(10, 1);
            Assert.AreEqual(0.0292526881, dist.CDF(0.5), 1E-6);
            Assert.AreEqual(0.4404932851, dist.CDF(1.0), 1E-6);
            Assert.AreEqual(0.8911780189, dist.CDF(2.0), 1E-6);

            // CDF must be non-decreasing
            Assert.IsLessThan(dist.CDF(1.0), dist.CDF(0.5));
            Assert.IsLessThan(dist.CDF(2.0), dist.CDF(1.0));

            // InverseCDF values from v*sigma / scipy.stats.chi2.isf(p, v)
            Assert.AreEqual(0.6255012152, dist.InverseCDF(0.1), 1E-4);
            Assert.AreEqual(0.8488443635, dist.InverseCDF(0.3), 1E-4);
            Assert.AreEqual(1.0704554778, dist.InverseCDF(0.5), 1E-4);
            Assert.AreEqual(1.3760423551, dist.InverseCDF(0.7), 1E-4);
            Assert.AreEqual(2.0554215430, dist.InverseCDF(0.9), 1E-4);

            // Roundtrip: InverseCDF(CDF(x)) ≈ x
            foreach (double x in new[] { 0.5, 1.0, 2.0 })
            {
                Assert.AreEqual(x, dist.InverseCDF(dist.CDF(x)), 1E-4, $"Roundtrip failed for x={x}");
            }

            // v=5, sigma=2: CDF values from scipy.stats.chi2.sf(v*sigma/x, v)
            var dist2 = new InverseChiSquared(5, 2);
            Assert.AreEqual(0.0012497306, dist2.CDF(0.5), 1E-6);
            Assert.AreEqual(0.4158801870, dist2.CDF(2.0), 1E-6);
            Assert.AreEqual(0.8491450361, dist2.CDF(5.0), 1E-6);

            // Roundtrip for v=5, sigma=2
            foreach (double x in new[] { 0.5, 2.0, 5.0 })
            {
                Assert.AreEqual(x, dist2.InverseCDF(dist2.CDF(x)), 1E-3, $"Roundtrip failed for x={x}");
            }
        }
    }
}
