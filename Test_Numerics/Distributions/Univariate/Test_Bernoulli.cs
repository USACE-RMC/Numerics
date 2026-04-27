using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// A class testing the Bernoulli Distribution algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> <see href="https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests/Discrete/BernoulliTests.cs"/></item>
    /// <item> <see href="https://en.wikipedia.org/wiki/Kurtosis"/></item>
    /// <item> <see href="https://blogs.sas.com/content/iml/2015/01/28/skewness-and-kurtosis.html"/></item>
    /// </list>
    /// </para>
    /// </remarks>

    [TestClass]
    public class Test_Bernoulli
    {
        /// <summary>
        /// Verified using Palisade's @Risk
        /// </summary>
        [TestMethod()]
        public void Test_Bernoulli_AtRisk()
        {
            double true_mean = 0.7d;
            int true_mode = 1;
            int true_median = 1;
            double true_stdDev = 0.4583d;
            double true_skew = -0.8729d;
            double true_kurt = 1.7619d;
            double true_pdf = 0.3d;
            double true_cdf = 0.3d;
            double true_icdf05 = 0.0d;
            double true_icdf95 = 1.0d;
            var B = new Bernoulli(0.7d);
            Assert.AreEqual(true_mean, B.Mean, 0.0001d);
            Assert.AreEqual(true_median, B.Median, 0.0001d);
            Assert.AreEqual(true_mode, B.Mode, 0.0001d);
            Assert.AreEqual(true_stdDev, B.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_skew, B.Skewness, 0.0001d);
            Assert.AreEqual(true_kurt, B.Kurtosis, 0.0001d);
            Assert.AreEqual(true_pdf, B.PDF(0.0d), 0.0001d);
            Assert.AreEqual(true_cdf, B.CDF(0.5d), 0.0001d);
            Assert.AreEqual(true_icdf05, B.InverseCDF(0.05d), 0.0001d);
            Assert.AreEqual(true_icdf95, B.InverseCDF(0.95d), 0.0001d);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. See if Bernoulli is being created.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var bernoulli = new Bernoulli(0);
            Assert.AreEqual(0, bernoulli.Probability);

            var bernoulli2 = new Bernoulli(0.3);
            Assert.AreEqual(0.3, bernoulli2.Probability);

            var bernoulli3 = new Bernoulli(1);
            Assert.AreEqual(1, bernoulli3.Probability);
        }
        
        /// <summary>
        /// Verified using MathNet-Numerics testing. See what probabilities fail.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var b = new Bernoulli(double.NaN);
            Assert.IsFalse(b.ParametersValid);

            var b2 = new Bernoulli(-1);
            Assert.IsFalse(b2.ParametersValid);

            var b3 = new Bernoulli(2);
            Assert.IsFalse(b3.ParametersValid);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking string output.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var b = new Bernoulli(0.3);
            Assert.AreEqual("Probability (p)", b.ParametersToString[0,0]);
            Assert.AreEqual("0.3",b.ParametersToString[0,1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Bernoulli(0.3);
            var mom = dist.CentralMoments(200);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking mean of distribution with different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var b = new Bernoulli(0);
            Assert.AreEqual(0, b.Mean);

            var b2 = new Bernoulli(0.3);
            Assert.AreEqual(0.3, b2.Mean);

            var b3 = new Bernoulli(1);
            Assert.AreEqual(1, b3.Mean);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking median of distribution with different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var b = new Bernoulli(0);
            Assert.AreEqual(0, b.Median);

            var b2 = new Bernoulli(0.4);
            Assert.AreEqual(0, b2.Median);

            var b3 = new Bernoulli(0.5);
            Assert.AreEqual(0.5, b3.Median);

            var b4 = new Bernoulli(0.6);
            Assert.AreEqual(1, b4.Median);

            var b5 = new Bernoulli(1);
            Assert.AreEqual(1, b5.Median);

        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking mode of distribution with different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var b = new Bernoulli(0);
            Assert.AreEqual(0, b.Mode);

            var b2 = new Bernoulli(0.3);
            Assert.AreEqual(0, b2.Mode);

            var b3 = new Bernoulli(1);
            Assert.AreEqual(1, b3.Mode);

            var b4 = new Bernoulli(0.5);
            Assert.AreEqual(0, b4.Mode);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking standard deviation with different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var b = new Bernoulli(0);
            Assert.AreEqual(0, b.StandardDeviation);

            var b2 = new Bernoulli(0.3);
            Assert.AreEqual(0.458257, b2.StandardDeviation, 1e-04);

            var b3 = new Bernoulli(1);
            Assert.AreEqual(0, b3.StandardDeviation);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var b = new Bernoulli(0.3);
            Assert.AreEqual(0, b.Minimum);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var b = new Bernoulli(0.3);
            Assert.AreEqual(1, b.Maximum);
        }

        /// <summary>
        /// Checking Kurtosis of distribution with different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            // scipy.stats.bernoulli(0).stats(moments='k') => nan
            var b = new Bernoulli(0);
            Assert.AreEqual(double.NaN, b.Kurtosis);

            var b2 = new Bernoulli(0.3d);
            Assert.AreEqual(1.761904762, b2.Kurtosis, 1e-04);

            // scipy.stats.bernoulli(1).stats(moments='k') => nan
            var b3 = new Bernoulli(1);
            Assert.AreEqual(double.NaN, b3.Kurtosis);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking skewness of distribution with different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            // scipy.stats.bernoulli(0).stats(moments='s') => nan
            var b = new Bernoulli(0d);
            Assert.AreEqual(double.NaN, b.Skewness);

            var b2 = new Bernoulli(0.3);
            Assert.AreEqual(0.8728715, b2.Skewness, 1e-04);

            // scipy.stats.bernoulli(1).stats(moments='s') => nan
            var b3 = new Bernoulli(1);
            Assert.AreEqual(double.NaN, b3.Skewness);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing PDF function.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var b = new Bernoulli(0);
            Assert.AreEqual(0, b.PDF(-1));
            Assert.AreEqual(1, b.PDF(0));
            Assert.AreEqual(0, b.PDF(1));
            Assert.AreEqual(0, b.PDF(2));

            var b2 = new Bernoulli(0.3);
            Assert.AreEqual(0, b2.PDF(-1));
            Assert.AreEqual(0.7, b2.PDF(0));
            Assert.AreEqual(0.3, b2.PDF(1));
            Assert.AreEqual(0,b2.PDF(2));

            var b3 = new Bernoulli(1);
            Assert.AreEqual(0, b3.PDF(-1));
            Assert.AreEqual(0, b3.PDF(0));
            Assert.AreEqual(1,b3.PDF(1));
            Assert.AreEqual(0,b3.PDF(-2));
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var b = new Bernoulli(0);
            Assert.AreEqual(0, b.CDF(-1));
            Assert.AreEqual(1, b.CDF(0));
            Assert.AreEqual(1, b.CDF(0.5));
            Assert.AreEqual(1, b.CDF(1));
            Assert.AreEqual(1, b.CDF(2));

            var b2 = new Bernoulli(0.3);
            Assert.AreEqual(0,b2.CDF(-1));
            Assert.AreEqual(0.7, b2.CDF(0));
            Assert.AreEqual(0.7, b2.CDF(0.5));
            Assert.AreEqual(1, b2.CDF(1));
            Assert.AreEqual(1, b2.CDF(2));

            var b3 = new Bernoulli(1);
            Assert.AreEqual(0,b3.CDF(-1));
            Assert.AreEqual(0, b3.CDF(0));
            Assert.AreEqual(0, b3.CDF(0.5));
            Assert.AreEqual(1, b3.CDF(1));
            Assert.AreEqual(1, b3.CDF(2));
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing InverseCDF function.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var b = new Bernoulli(0);
            Assert.AreEqual(0, b.InverseCDF(0));

            var b2 = new Bernoulli(0.3);
            Assert.AreEqual(0, b2.InverseCDF(0.3));

            var b3 = new Bernoulli(1);
            Assert.AreEqual(1, b3.InverseCDF(1));

            var b4 = new Bernoulli(0.7);
            Assert.AreEqual(1, b4.InverseCDF(0.7));

        }

        /// <summary>
        /// Verify Bernoulli returns NaN for skewness/kurtosis at boundary probabilities p=0 and p=1.
        /// Reference: scipy.stats.bernoulli(p).stats(moments='sk') returns nan at p=0,1.
        /// </summary>
        [TestMethod()]
        public void Test_SkewnessKurtosis_Boundary()
        {
            // p=0: skewness and kurtosis are undefined (p*q = 0)
            var b0 = new Bernoulli(0);
            Assert.AreEqual(double.NaN, b0.Skewness);
            Assert.AreEqual(double.NaN, b0.Kurtosis);

            // p=1: same
            var b1 = new Bernoulli(1);
            Assert.AreEqual(double.NaN, b1.Skewness);
            Assert.AreEqual(double.NaN, b1.Kurtosis);

            // p=0.5: well-defined
            // scipy.stats.bernoulli(0.5): skew=0.0, excess_kurt=-2.0, raw_kurt=1.0
            // Note: this library uses raw kurtosis (= excess + 3)
            var b5 = new Bernoulli(0.5);
            Assert.AreEqual(0.0, b5.Skewness, 1E-10);
            Assert.AreEqual(1.0, b5.Kurtosis, 1E-10);
        }
    }
}
