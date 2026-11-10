using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Rayleigh distribution algorithm.
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
    public class Test_Rayleigh
    {

        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_RayleighDist()
        {
            double true_mean = 0.52639193767251d;
            double true_median = 0.49451220943852386d;
            double true_stdDev = Math.Sqrt(0.075711527953380237d);
            double true_pdf = 0.030681905868831811d;
            double true_cdf = 0.99613407986052716d;
            double true_icdf = 1.4d;
            var R = new Rayleigh(0.42);
            Assert.AreEqual(true_mean, R.Mean, 0.0001d);
            Assert.AreEqual(true_median, R.Median, 0.0001d);
            Assert.AreEqual(true_stdDev, R.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_pdf, R.PDF(1.4d), 0.0001d);
            Assert.AreEqual(true_cdf, R.CDF(1.4d), 0.0001d);
            Assert.AreEqual(true_icdf, R.InverseCDF(true_cdf), 0.0001d);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var R = new Rayleigh();
            Assert.AreEqual(10, R.Sigma);

            var R2 = new Rayleigh(2);
            Assert.AreEqual(2, R2.Sigma);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var R = new Rayleigh(double.NaN);
            Assert.IsFalse(R.ParametersValid);

            var R2 = new Rayleigh(double.PositiveInfinity);
            Assert.IsFalse(R2.ParametersValid);

            var R3 = new Rayleigh(0);
            Assert.IsFalse(R3.ParametersValid);
        }

        /// <summary>
        /// Testing Parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var R = new Rayleigh();
            Assert.AreEqual("Scale (σ)", R.ParametersToString[0, 0]);
            Assert.AreEqual("10", R.ParametersToString[0, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Rayleigh(1);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var R = new Rayleigh();
            Assert.AreEqual(12.53314, R.Mean, 1e-04);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(1.25331, R2.Mean, 1e-04);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var R = new Rayleigh();
            Assert.AreEqual(11.7741, R.Median, 1e-04);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(1.1774, R2.Median, 1e-04);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var R = new Rayleigh();
            Assert.AreEqual(10, R.Mode);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(1, R2.Mode);
        }
        
        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var R = new Rayleigh();
            Assert.AreEqual(6.55136, R.StandardDeviation, 1e-05);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(0.65513, R2.StandardDeviation, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var R = new Rayleigh();
            Assert.AreEqual(0.63111, R.Skewness, 1e-04);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(0.63111, R2.Skewness, 1e-04);
        }

        /// <summary>
        /// Testing kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var R = new Rayleigh();
            Assert.AreEqual(3.24508, R.Kurtosis, 1e-05);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(3.24508, R2.Kurtosis, 1e-05);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var R = new Rayleigh();
            Assert.AreEqual(0, R.Minimum);
            Assert.AreEqual(double.PositiveInfinity, R.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var R = new Rayleigh();
            Assert.AreEqual(0, R.PDF(-1));
            Assert.AreEqual(9.9501e-03, R.PDF(1), 1e-06);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(0.019603, R.PDF(2), 1e-05);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var R = new Rayleigh();
            Assert.AreEqual(0, R.CDF(-1));
            Assert.AreEqual(4.9875e-03, R.CDF(1), 1e-04);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var R = new Rayleigh();
            Assert.AreEqual(0, R.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity, R.InverseCDF(1));
            Assert.AreEqual(10.1076, R.InverseCDF(0.4), 1e-04);
        }

        /// <summary>
        /// Verify Rayleigh MoM estimation recovers correct sigma from sample mean.
        /// Reference: scipy.stats.rayleigh(scale=sigma).mean() = sigma * sqrt(pi/2).
        /// Therefore sigma = mean / sqrt(pi/2), NOT sigma = stddev.
        /// </summary>
        [TestMethod()]
        public void Test_MoM_Estimation()
        {
            // Generate a Rayleigh sample with known sigma=3.0 using a fixed seed
            var sigma = 3.0;
            var R = new Rayleigh(sigma);
            // rayleigh(scale=3.0).mean() = 3.7599424119
            Assert.AreEqual(3.7599424119, R.Mean, 1E-4);
            // rayleigh(scale=3.0).std() = 1.9652...
            // MoM: sigma = mean / sqrt(pi/2) = 3.0
            Assert.AreEqual(sigma, R.Mean / Math.Sqrt(Math.PI / 2.0), 1E-10);
        }

        /// <summary>
        /// Ordinary values match scipy.stats.rayleigh to near machine precision.
        /// </summary>
        [TestMethod()]
        public void Test_OrdinaryValues_Scipy()
        {
            var R = new Rayleigh(2);
            // scipy.stats.rayleigh.pdf/cdf/logpdf(x, scale=2)
            Assert.AreEqual(0.12115415430954302, R.PDF(0.5), 0.12115415430954302 * 1E-14);
            Assert.AreEqual(0.030766765523655912, R.CDF(0.5), 0.030766765523655912 * 1E-14);
            Assert.AreEqual(-2.1106915416798357, R.LogPDF(0.5), 2.1106915416798357 * 1E-14);
            Assert.AreEqual(0.3032653298563167, R.PDF(2.0), 0.3032653298563167 * 1E-14);
            Assert.AreEqual(0.3934693402873666, R.CDF(2.0), 0.3934693402873666 * 1E-14);
            Assert.AreEqual(-1.1931471805599454, R.LogPDF(2.0), 1.1931471805599454 * 1E-14);
            Assert.AreEqual(0.00233517133666736, R.PDF(7.3), 0.00233517133666736 * 1E-14);
            Assert.AreEqual(0.9987204540621001, R.CDF(7.3), 0.9987204540621001 * 1E-14);
            Assert.AreEqual(-6.059670012965546, R.LogPDF(7.3), 6.059670012965546 * 1E-14);
            // scipy.stats.rayleigh.ppf(p, scale=2)
            Assert.AreEqual(0.2835536753914707, R.InverseCDF(0.01), 0.2835536753914707 * 1E-14);
            Assert.AreEqual(2.0215353051895795, R.InverseCDF(0.4), 2.0215353051895795 * 1E-14);
            Assert.AreEqual(6.069708517540585, R.InverseCDF(0.99), 6.069708517540585 * 1E-14);
        }

        /// <summary>
        /// Small arguments, small probabilities and far tails keep full precision.
        /// </summary>
        [TestMethod()]
        public void Test_Tails()
        {
            var R = new Rayleigh(1);
            // scipy.stats.rayleigh.cdf(1e-10)
            Assert.AreEqual(5.0000000000000005E-21, R.CDF(1E-10), 5E-21 * 1E-14);
            // scipy.stats.rayleigh.ppf(1e-20)
            Assert.AreEqual(1.414213562373095E-10, R.InverseCDF(1E-20), 1.414213562373095E-10 * 1E-14);
            // scipy.stats.rayleigh.logpdf(50)
            Assert.AreEqual(-1246.0879769945718, R.LogPDF(50), 1246.0879769945718 * 1E-15);
            // scipy.stats.rayleigh.logsf(50)
            Assert.AreEqual(-1250d, R.LogCCDF(50), 0d);
            // scipy.stats.rayleigh.sf(30)
            Assert.AreEqual(3.693883068487256E-196, R.CCDF(30), 3.693883068487256E-196 * 1E-13);
            Assert.AreEqual(1d, R.CCDF(0));
            Assert.AreEqual(0d, R.LogCCDF(-1));
            Assert.AreEqual(double.NegativeInfinity, R.LogPDF(0));
            Assert.AreEqual(double.NegativeInfinity, R.LogPDF(-1));
            Assert.AreEqual(0d, R.PDF(double.PositiveInfinity));
            Assert.AreEqual(double.NegativeInfinity, R.LogPDF(double.PositiveInfinity));
            Assert.AreEqual(1d, R.CDF(double.PositiveInfinity));
            Assert.AreEqual(0d, R.CCDF(double.PositiveInfinity));
            Assert.IsTrue(double.IsNaN(R.PDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(R.CDF(double.NaN)));
        }

        /// <summary>
        /// A huge scale does not overflow through σ².
        /// </summary>
        [TestMethod()]
        public void Test_HugeScale()
        {
            var R = new Rayleigh(1E200);
            // scipy.stats.rayleigh.pdf(1e200, scale=1e200)
            Assert.AreEqual(6.0653065971263344E-201, R.PDF(1E200), 6.0653065971263344E-201 * 1E-14);
            // scipy.stats.rayleigh.cdf(1e200, scale=1e200)
            Assert.AreEqual(0.3934693402873666, R.CDF(1E200), 0.3934693402873666 * 1E-14);
        }

        /// <summary>
        /// Maximum likelihood applies the Γ(n)√n/Γ(n + 1/2) bias correction to the closed-form estimator.
        /// </summary>
        [TestMethod()]
        public void Test_MLE_BiasCorrected()
        {
            var data = new[] { 1.2, 0.4, 2.5, 3.1, 0.9, 1.7, 2.2 };
            var R = new Rayleigh();
            R.Estimate(data, ParameterEstimationMethod.MaximumLikelihood);
            // sqrt(sum(x^2)/(2n)) = scipy.stats.rayleigh.fit(data, floc=0)[1] = 1.3627702877384937, times
            // exp(gammaln(7) + 0.5*log(7) - gammaln(7.5)) from scipy.special
            Assert.AreEqual(1.3873031101602713, R.Sigma, 1.3873031101602713 * 1E-13);

            R.Estimate(data, ParameterEstimationMethod.MethodOfMoments);
            // numpy.mean(data) / sqrt(pi/2)
            Assert.AreEqual(1.3678021042334836, R.Sigma, 1.3678021042334836 * 1E-14);
        }

        /// <summary>
        /// Invalid samples throw instead of producing a silent or negative scale.
        /// </summary>
        [TestMethod()]
        public void Test_Estimate_InvalidSample()
        {
            var R = new Rayleigh();
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => R.Estimate(new double[0], ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => R.Estimate(new[] { -1.0, -2.0, -1.5 }, ParameterEstimationMethod.MethodOfMoments));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => R.Estimate(new[] { 0.0, 0.0 }, ParameterEstimationMethod.MethodOfMoments));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => R.Estimate(new[] { 1.0, double.NaN }, ParameterEstimationMethod.MethodOfMoments));
            Assert.ThrowsExactly<ArgumentNullException>(() => R.Estimate(null!, ParameterEstimationMethod.MethodOfMoments));
            // A single positive observation determines the scale.
            R.Estimate(new[] { 2.0 }, ParameterEstimationMethod.MethodOfMoments);
            Assert.AreEqual(2.0 / Math.Sqrt(Math.PI / 2.0), R.Sigma, 1E-15);
        }

        /// <summary>
        /// The parametric bootstrap re-estimates a valid distribution from a seeded sample.
        /// </summary>
        [TestMethod()]
        public void Test_Bootstrap()
        {
            var R = new Rayleigh(3);
            var b = (Rayleigh)R.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 500, 12345);
            Assert.IsTrue(b.ParametersValid);
            Assert.AreEqual(3d, b.Sigma, 0.3);
        }
    }
}
