using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Pareto distribution algorithm.
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
    public class Test_Pareto
    {

        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_ParetoDist()
        {
            double true_mean = 0.63d;
            double true_mode = 0.42d;
            double true_median = 0.52916684095584676d;
            double true_stdDev = Math.Sqrt(0.13229999999999997d);
            double true_pdf = 0.057857142857142857d;
            double true_cdf = 0.973d;
            double true_icdf = 1.4d;
            double true_icdf05 = 0.4272d;
            double true_icdf95 = 1.1401d;
            var PA = new Pareto(0.42d, 3d);
            Assert.AreEqual(true_mean, PA.Mean, 0.0001d);
            Assert.AreEqual(true_median, PA.Median, 0.0001d);
            Assert.AreEqual(true_mode, PA.Mode, 0.0001d);
            Assert.AreEqual(true_stdDev, PA.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_pdf, PA.PDF(1.4d), 0.0001d);
            Assert.AreEqual(true_cdf, PA.CDF(1.4d), 0.0001d);
            Assert.AreEqual(true_icdf, PA.InverseCDF(PA.CDF(1.4d)), 0.0001d);
            Assert.AreEqual(true_icdf05, PA.InverseCDF(0.05d), 0.0001d);
            Assert.AreEqual(true_icdf95, PA.InverseCDF(0.95d), 0.0001d);
            PA.SetParameters(new[] { 1d, 10d });
            double true_skew = 2.8111d;
            double true_kurt = 17.8286d;
            Assert.AreEqual(true_skew, PA.Skewness, 0.0001d);
            Assert.AreEqual(true_kurt, PA.Kurtosis, 0.0001d);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var p = new Pareto();
            Assert.AreEqual(1, p.Xm);
            Assert.AreEqual(10, p.Alpha);

            var p2 = new Pareto(1, 1);
            Assert.AreEqual(1, p2.Xm);
            Assert.AreEqual(1, p2.Alpha);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var p = new Pareto(double.NaN,double.PositiveInfinity);
            Assert.IsFalse(p.ParametersValid);

            var p2 = new Pareto(double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(p2.ParametersValid);

            var p3 = new Pareto(0, 0);
            Assert.IsFalse(p3.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var p = new Pareto();
            Assert.AreEqual("Scale (Xm)", p.ParametersToString[0, 0]);
            Assert.AreEqual("Shape (α)", p.ParametersToString[1, 0]);
            Assert.AreEqual("1", p.ParametersToString[0, 1]);
            Assert.AreEqual("10", p.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Pareto();
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
            var p = new Pareto();
            Assert.AreEqual(1.1111, p.Mean,  1e-04);

            var p2 = new Pareto(1,1);
            Assert.AreEqual(double.PositiveInfinity, p2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var p = new Pareto();
            Assert.AreEqual(1.07177, p.Median,  1e-04);

            var p2 = new Pareto(1, 1);
            Assert.AreEqual(2, p2.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var p = new Pareto();
            Assert.AreEqual(1, p.Mode);

            var p2 = new Pareto(2, 1);
            Assert.AreEqual(2, p2.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var p = new Pareto();
            Assert.AreEqual(0.12422, p.StandardDeviation,  1e-04);

            var p2 = new Pareto(1, 2);
            Assert.AreEqual(double.PositiveInfinity, p2.StandardDeviation);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var p = new Pareto();
            Assert.AreEqual(2.81105, p.Skewness, 1e-04);

            var p2 = new Pareto(1, 3);
            Assert.AreEqual(double.NaN, p2.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var p = new Pareto();
            Assert.AreEqual(17.82857, p.Kurtosis, 1e-04);

            var p2 = new Pareto(1, 4);
            Assert.AreEqual(double.NaN, p2.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var p = new Pareto();
            Assert.AreEqual(1, p.Minimum);
            Assert.AreEqual(double.PositiveInfinity, p.Maximum);

            var p2 = new Pareto(2,3);
            Assert.AreEqual(2, p2.Minimum);
            Assert.AreEqual(double.PositiveInfinity, p2.Maximum);
        }

        /// <summary>
        /// Testing PDF
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var p = new Pareto(1,1);
            Assert.AreEqual(1, p.PDF(1));
            Assert.AreEqual(4d / 9d,p.PDF(1.5) );

            var p2 = new Pareto(3, 2);
            Assert.AreEqual(2d / 3d,p2.PDF(3));
            Assert.AreEqual(18d / 125d,p2.PDF(5) );
        }

        /// <summary>
        /// Testing CDF.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var p = new Pareto();
            Assert.AreEqual(0, p.CDF(0));
            Assert.AreEqual(0.9990, p.CDF(2),  1e-04);
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var p = new Pareto();
            Assert.AreEqual(1, p.InverseCDF(0));
            Assert.AreEqual(1.0363, p.InverseCDF(0.3), 1e-04);
        }

        /// <summary>
        /// The density stays finite where Xm^α or x^(α+1) alone overflows.
        /// </summary>
        [TestMethod()]
        public void Test_PDF_LargeShapeAndScale()
        {
            // scipy.stats.pareto.pdf(11, 400, scale=10)
            Assert.AreEqual(1.0083080695310864E-15, new Pareto(10, 400).PDF(11), 1.0083080695310864E-15 * 1E-13);
            // scipy.stats.pareto.pdf(10, 400, scale=10)
            Assert.AreEqual(40d, new Pareto(10, 400).PDF(10), 40d * 1E-15);
            // scipy.stats.pareto.pdf(2e200, 2, scale=1e200)
            Assert.AreEqual(2.5E-201, new Pareto(1E200, 2).PDF(2E200), 2.5E-201 * 1E-15);
            // scipy.stats.pareto.pdf(1.4, 3, scale=0.42)
            Assert.AreEqual(0.05785714285714288, new Pareto(0.42, 3).PDF(1.4), 0.05785714285714288 * 1E-14);
            var p = new Pareto(0.42, 3);
            Assert.AreEqual(0d, p.PDF(0.41));
            Assert.AreEqual(0d, p.PDF(double.PositiveInfinity));
            Assert.IsTrue(double.IsNaN(p.PDF(double.NaN)));
        }

        /// <summary>
        /// The log density, survival function and log survival function are evaluated directly.
        /// </summary>
        [TestMethod()]
        public void Test_LogPDF_CCDF_LogCCDF()
        {
            var p = new Pareto(0.42, 3);
            // scipy.stats.pareto.logpdf(1.4, 3, scale=0.42)
            Assert.AreEqual(-2.849778360930911, p.LogPDF(1.4), 2.849778360930911 * 1E-14);
            // ln(2) - 3 ln(1e200) - ln(1e200) for Pareto(1, 2) at 1e200; scipy.stats.pareto.logpdf(1e200, 2) underflows to -inf
            double expected = Math.Log(2d) - 3d * Math.Log(1E200);
            Assert.AreEqual(expected, new Pareto(1, 2).LogPDF(1E200), Math.Abs(expected) * 1E-14);
            Assert.AreEqual(double.NegativeInfinity, p.LogPDF(0.4));

            var q = new Pareto(1, 2);
            // scipy.stats.pareto.sf(1e10, 2)
            Assert.AreEqual(1E-20, q.CCDF(1E10), 1E-20 * 1E-14);
            // scipy.stats.pareto.logsf(1e10, 2)
            Assert.AreEqual(-46.051701859880914, q.LogCCDF(1E10), 46.051701859880914 * 1E-15);
            // scipy.stats.pareto.sf(1.4, 3, scale=0.42) = 1 - 0.973
            Assert.AreEqual(0.027, p.CCDF(1.4), 0.027 * 1E-13);
            Assert.AreEqual(1d, p.CCDF(0.42));
            Assert.AreEqual(0d, p.LogCCDF(0.3));
            Assert.AreEqual(0d, p.CCDF(double.PositiveInfinity));
            Assert.AreEqual(double.NegativeInfinity, p.LogCCDF(double.PositiveInfinity));
        }

        /// <summary>
        /// The CDF keeps full relative precision just above the scale parameter.
        /// </summary>
        [TestMethod()]
        public void Test_CDF_NearScale()
        {
            // mpmath: 1 - (1/x)**2 at x = 1 + 1e-12 (the double nearest), 50 digits
            Assert.AreEqual(2.0001778011616815E-12, new Pareto(1, 2).CDF(1 + 1E-12), 2.0001778011616815E-12 * 1E-14);
            // mpmath: 1 - (1e5/x)**0.5 at x = 1e5 + 1e-7 (the double nearest), 50 digits
            Assert.AreEqual(5.0000380724631E-13, new Pareto(1E5, 0.5).CDF(1E5 + 1E-7), 5.0000380724631E-13 * 1E-12);
            // scipy.stats.pareto.cdf(1.4, 3, scale=0.42)
            Assert.AreEqual(0.973, new Pareto(0.42, 3).CDF(1.4), 0.973 * 1E-15);
            Assert.AreEqual(1d, new Pareto(0.42, 3).CDF(double.PositiveInfinity));
        }

        /// <summary>
        /// Maximum likelihood has the closed form Xm = min(x), α = n / Σ ln(x/Xm).
        /// </summary>
        [TestMethod()]
        public void Test_MLE()
        {
            var data = new[] { 1.2, 3.4, 1.05, 2.2, 7.9, 1.6, 2.8, 1.33 };
            var p = new Pareto();
            p.Estimate(data, ParameterEstimationMethod.MaximumLikelihood);
            // scipy.stats.pareto.fit(data, floc=0) = (b=1.402355408973135, loc=0, scale=1.05)
            Assert.AreEqual(1.05, p.Xm, 0d);
            Assert.AreEqual(1.402355408973135, p.Alpha, 1.402355408973135 * 1E-14);

            var constraints = p.GetParameterConstraints(data);
            Assert.AreEqual(1.05, constraints.Item1[0], 0d);
            Assert.AreEqual(1.05, constraints.Item3[0], 0d);
            Assert.IsTrue(constraints.Item2[0] > 0 && constraints.Item2[0] < 1.05);
            Assert.IsTrue(constraints.Item2[1] < p.Alpha && p.Alpha < constraints.Item3[1]);

            Assert.ThrowsExactly<NotImplementedException>(() => p.Estimate(data, ParameterEstimationMethod.MethodOfMoments));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => p.Estimate(new[] { 1.0, -2.0, 3.0 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => p.Estimate(new[] { 2.0, 2.0, 2.0 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => p.Estimate(new[] { 2.0 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentNullException>(() => p.Estimate(null!, ParameterEstimationMethod.MaximumLikelihood));
        }

        /// <summary>
        /// The parametric bootstrap re-estimates a valid distribution from a seeded sample.
        /// </summary>
        [TestMethod()]
        public void Test_Bootstrap()
        {
            var p = new Pareto(2, 3);
            var b = (Pareto)p.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 200, 12345);
            Assert.IsTrue(b.ParametersValid);
            Assert.IsTrue(b.Xm >= 2d && b.Xm < 2.1);
            Assert.AreEqual(3d, b.Alpha, 1d);
        }

        /// <summary>
        /// The quantile function validates the parameters before returning the endpoints.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_ValidatesParametersFirst()
        {
            var p = new Pareto(1, 2);
            p.Xm = -1;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => p.InverseCDF(0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => p.InverseCDF(1));
        }
    }
}
