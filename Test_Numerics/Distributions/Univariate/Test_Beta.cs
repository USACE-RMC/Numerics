using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Beta distribution algorithm.
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
    /// <list type="bullet">
    /// <item> <see href = "https://en.wikipedia.org/wiki/Beta_distribution" /></item>
    /// <item> <see href="https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests"/></item>
    /// </list>
    /// </remarks>
    [TestClass]
    public class Test_Beta
    {
        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_BetaDist()
        {
            double true_mean = 0.21105527638190955d;
            double true_mode = 0.0d;
            double true_median = 0.11577711097114812d;
            double true_stdDev = Math.Sqrt(0.055689279830523512d);
            double true_skew = 1.2275d;
            double true_kurt = 3.6048d;
            double true_pdf = 0.94644031936694828d;
            double true_cdf = 0.69358638272337991d;
            double true_icdf = 0.27d;
            double true_icdf05 = 0.000459d;
            double true_icdf95 = 0.7238d;
            var B = new BetaDistribution(0.42d, 1.57d);
            Assert.AreEqual(true_mean, B.Mean, 0.0001d);
            Assert.AreEqual(true_median, B.Median, 0.0001d);
            Assert.AreEqual(true_mode, B.Mode, 0.0001d);
            Assert.AreEqual(true_stdDev, B.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_skew, B.Skewness, 0.0001d);
            Assert.AreEqual(true_kurt, B.Kurtosis, 0.0001d);
            Assert.AreEqual(true_pdf, B.PDF(0.27d), 0.0001d);
            Assert.AreEqual(true_cdf, B.CDF(0.27d), 0.0001d);
            Assert.AreEqual(true_icdf, B.InverseCDF(B.CDF(0.27d)), 0.0001d);
            Assert.AreEqual(true_icdf05, B.InverseCDF(0.05d), 0.0001d);
            Assert.AreEqual(true_icdf95, B.InverseCDF(0.95d), 0.0001d);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. See if beta is being created.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var b = new BetaDistribution(0, 0);
            Assert.AreEqual(0, b.Alpha);
            Assert.AreEqual(0, b.Beta);

            var b2 = new BetaDistribution(0, 1);
            Assert.AreEqual(0, b2.Alpha);
            Assert.AreEqual(1, b2.Beta);

            var b3 = new BetaDistribution(1, 0);
            Assert.AreEqual(1, b3.Alpha);
            Assert.AreEqual(0, b3.Beta);

            var b4 = new BetaDistribution(1, 1);
            Assert.AreEqual(1, b4.Alpha);
            Assert.AreEqual(1, b4.Beta);

            var b5 = new BetaDistribution(9, 1);
            Assert.AreEqual(9, b5.Alpha);
            Assert.AreEqual(1, b5.Beta);

            var b6 = new BetaDistribution(5, 100);
            Assert.AreEqual(5, b6.Alpha);
            Assert.AreEqual(100, b6.Beta);

            var b7 = new BetaDistribution(1,double.PositiveInfinity);
            Assert.AreEqual(1, b7.Alpha);
            Assert.AreEqual(double.PositiveInfinity, b7.Beta);

            var b8 = new BetaDistribution(double.PositiveInfinity,1);
            Assert.AreEqual(double.PositiveInfinity, b8.Alpha);
            Assert.AreEqual(1, b8.Beta);

            var b9 = new BetaDistribution(0,double.PositiveInfinity);
            Assert.AreEqual(0, b9.Alpha);
            Assert.AreEqual(double.PositiveInfinity, b9.Beta);

            var b10 = new BetaDistribution(double.PositiveInfinity, 0);
            Assert.AreEqual(double.PositiveInfinity, b10.Alpha);
            Assert.AreEqual(0,b10.Beta);
        }

        /// <summary>
        /// Checking failure with valid parameter variable.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var b = new BetaDistribution(double.NaN, 0);
            Assert.IsFalse(b.ParametersValid);

            var b2 = new BetaDistribution(-1, 1);
            Assert.IsFalse(b2.ParametersValid);

            var b3 = new BetaDistribution(double.PositiveInfinity, 0);
            Assert.IsFalse(b3.ParametersValid);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking string output.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var b = new BetaDistribution(1d, 1d);
            Assert.AreEqual("Shape (α)", b.ParametersToString[0, 0]);
            Assert.AreEqual("Shape (β)", b.ParametersToString[1, 0]);
            Assert.AreEqual("1", b.ParametersToString[0, 1]);
            Assert.AreEqual("1",b.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new BetaDistribution(9, 1);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        ///  Verified using MathNet-Numerics testing. Checking mean of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var b = new BetaDistribution(1,2);
            Assert.AreEqual((double)1/3,b.Mean);

            var b2 = new BetaDistribution(9,1);
            Assert.AreEqual(0.9, b2.Mean);

            var b3 = new BetaDistribution(1,double.PositiveInfinity);
            Assert.AreEqual(0, b3.Mean);

            var b4 = new BetaDistribution(0, 0);
            Assert.AreEqual(double.NaN, b4.Mean);

            var b5 = new BetaDistribution(double.PositiveInfinity, 1);
            Assert.AreEqual(double.NaN, b5.Mean);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking median of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var b = new BetaDistribution(2, 2);
            Assert.AreEqual(0.5, b.Median);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking mode of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var b = new BetaDistribution(1, 1);
            Assert.AreEqual(0.5, b.Mode);

            var b2 = new BetaDistribution(0, 2);
            Assert.AreEqual(0, b2.Mode);

            var b3 = new BetaDistribution(2, 0);
            Assert.AreEqual(1, b3.Mode);

            var b4 = new BetaDistribution(9, 1);
            Assert.AreEqual(1, b4.Mode);

            var lowerBoundary = new BetaDistribution(0.5d, 1.0d);
            Assert.AreEqual(0.0d, lowerBoundary.Mode);

            var upperBoundary = new BetaDistribution(1.0d, 0.5d);
            Assert.AreEqual(1.0d, upperBoundary.Mode);

            var uShaped = new BetaDistribution(0.5d, 0.5d);
            Assert.AreEqual(0.5d, uShaped.Mode);

            var interior = new BetaDistribution(2.0d, 5.0d);
            Assert.AreEqual(0.2d, interior.Mode, 1E-14);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking standard deviation with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var b = new BetaDistribution(2, 3);
            Assert.AreEqual(0.2, b.StandardDeviation);

            var b2 = new BetaDistribution(5, 1);
            Assert.AreEqual(0.14085,b2.StandardDeviation,1e-04);

            var b3 = new BetaDistribution(7, 3);
            Assert.AreEqual(0.138169, b3.StandardDeviation, 1e-04);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking skewness of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var b = new BetaDistribution(0, 0);
            Assert.AreEqual(double.NaN, b.Skewness);

            var b4 = new BetaDistribution(1, 1);
            Assert.AreEqual(0, b4.Skewness);

            var b5 = new BetaDistribution(9, 1);
            Assert.AreEqual(-1.474055, b5.Skewness,1e-04);
        }

        /// <summary>
        /// Checking Kurtosis of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var b = new BetaDistribution(2, 2);
            Assert.AreEqual(2.14285,b.Kurtosis,1e-04);

            var b2 = new BetaDistribution(5,2);
            Assert.AreEqual(2.88, b2.Kurtosis);

            var b3 = new BetaDistribution(2, 5);
            Assert.AreEqual(2.88,b3.Kurtosis);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var b = new BetaDistribution(1, 2);
            Assert.AreEqual(0,b.Minimum);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var b = new BetaDistribution(1, 2);
            Assert.AreEqual(1, b.Maximum);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing PDF function.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var b = new BetaDistribution(1, 1);
            Assert.AreEqual(1, b.PDF(0), 1E-4);
            Assert.AreEqual(1, b.PDF(0.5), 1E-4);
            Assert.AreEqual(1, b.PDF(1), 1E-4);

            var b2 = new BetaDistribution(9, 1);
            Assert.AreEqual(0, b2.PDF(0), 1E-4);
            Assert.AreEqual(0.035156, b2.PDF(0.5),1e-04);
            Assert.AreEqual(9.0, b2.PDF(1),1e-04);
            Assert.AreEqual(0, b2.PDF(-1), 1E-4);
            Assert.AreEqual(0, b2.PDF(2), 1E-4);

            var b3 = new BetaDistribution(5, 100);
            Assert.AreEqual(0, b3.PDF(0), 1E-4);
            Assert.AreEqual(0, b3.PDF(1), 1E-4);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var b = new BetaDistribution(1, 1);
            Assert.AreEqual(0, b.CDF(0));
            Assert.AreEqual(0.5, b.CDF(0.5));
            Assert.AreEqual(1, b.CDF(1));

            var b2 = new BetaDistribution(9, 1);
            Assert.AreEqual(0,b2.CDF(0));
            Assert.AreEqual(0.001953125, b2.CDF(0.5));
            Assert.AreEqual(1, b2.CDF(1));

            var b3 = new BetaDistribution(5, 100);
            Assert.AreEqual(0, b3.CDF(0));
            Assert.AreEqual(1, b3.CDF(1));
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing InverseCDF function.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var b = new BetaDistribution(1,1);
            Assert.AreEqual(1, b.InverseCDF(1));

            var b2 = new BetaDistribution(9, 1);
            Assert.AreEqual(0,b2.InverseCDF(0));
            Assert.AreEqual(1, b2.InverseCDF(1));

            var b3 = new BetaDistribution(5, 100);
            Assert.AreEqual(0,b3.InverseCDF(0));
        }

        /// <summary>
        /// Verifies the log-space density at large shapes and at the support endpoints.
        /// Reference values from scipy.stats.beta.pdf and scipy.stats.beta.logpdf, and mpmath at 30 digits.
        /// </summary>
        [TestMethod()]
        public void Test_LogPDF_LargeShapesAndEndpoints()
        {
            var large = new BetaDistribution(1000, 1000);
            Assert.AreEqual(35.678022291708641, large.PDF(0.5), 35.678022291708641 * 1E-12);
            Assert.AreEqual(3.574534877131522, large.LogPDF(0.5), 1E-12);

            var b = new BetaDistribution(2, 3);
            Assert.AreEqual(1.764, b.PDF(0.3), 1E-14);
            Assert.AreEqual(0.5675839575845993, b.LogPDF(0.3), 1E-14);
            Assert.AreEqual(double.NegativeInfinity, b.LogPDF(-0.1));
            Assert.AreEqual(double.NegativeInfinity, b.LogPDF(1.1));

            Assert.AreEqual(0.7351051938957226, new BetaDistribution(0.5, 0.5).PDF(0.25), 1E-14);
            Assert.AreEqual(double.PositiveInfinity, new BetaDistribution(0.5, 0.5).PDF(0));
            Assert.AreEqual(1.0, new BetaDistribution(1, 1).PDF(0), 1E-14);
            Assert.AreEqual(1.0, new BetaDistribution(1, 1).PDF(1), 1E-14);
            Assert.AreEqual(3.0, new BetaDistribution(1, 3).PDF(0), 1E-13);
            Assert.AreEqual(3.0, new BetaDistribution(3, 1).PDF(1), 1E-13);
            Assert.AreEqual(0.0, new BetaDistribution(2, 2).PDF(0));
        }

        /// <summary>
        /// Verifies the upper tail is resolved below machine epsilon.
        /// Reference values from scipy.stats.beta.sf and scipy.stats.beta.logsf.
        /// </summary>
        [TestMethod()]
        public void Test_CCDF_UpperTail()
        {
            var b = new BetaDistribution(2, 3);
            double x = 1d - 1E-12;
            Assert.AreEqual(3.9997345452279445E-36, b.CCDF(x), 3.9997345452279445E-36 * 1E-12);
            Assert.AreEqual(-81.50683535256094, b.LogCCDF(x), 1E-12);
            Assert.AreEqual(0.4752, b.CCDF(0.4), 1E-14);
            Assert.AreEqual(1d, b.CCDF(0));
            Assert.AreEqual(0d, b.CCDF(1));
        }

        /// <summary>
        /// Verifies the lower tail of the CCDF keeps a small x that 1 - x would round away.
        /// Reference values from scipy.stats.beta.sf (SciPy 1.17.1).
        /// </summary>
        [TestMethod]
        public void Test_CCDF_LowerTail()
        {
            var b = new BetaDistribution(0.3, 0.4);
            Assert.AreEqual(0.999999, b.CCDF(4.159739860874255E-20), 1E-15);
            Assert.AreEqual(0.9999000121244137, b.CCDF(1.93E-13), 1E-15);
        }

        /// <summary>
        /// Verifies a not-a-number probability is rejected by the inverse CDF.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_NaNProbability()
        {
            var b = new BetaDistribution(2, 3);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => b.InverseCDF(double.NaN));
        }

        /// <summary>
        /// Verifies the CDF when a shape parameter, or the sum of the shapes, lies one ulp below an integer.
        /// Reference values from scipy.stats.beta.cdf; the first matches the closed form 1 - (1 - x)^1.5.
        /// </summary>
        [TestMethod()]
        public void Test_CDF_ShapesJustBelowIntegers()
        {
            Assert.AreEqual(0.025000000000000012, new BetaDistribution(0.9999999999999998, 1.5).CDF(0.01673689415961626), 1E-15);
            Assert.AreEqual(0.6607459491435457, new BetaDistribution(0.49999999999999994, 1.4999999999999998).CDF(0.3), 1E-14);
            Assert.AreEqual(0.16000000000000003, new BetaDistribution(1.9999999999999998, 0.9999999999999998).CDF(0.4), 1E-15);
        }
    }
}
