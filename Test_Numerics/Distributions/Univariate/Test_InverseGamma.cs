using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Inverse Gamma distribution algorithm.
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
    public class Test_InverseGamma
    {
        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_InverseGammaDist()
        {
            double true_median = 3.1072323347401709d;
            double true_pdf = 0.35679850067181362d;
            double true_cdf = 0.042243552114989695d;
            double true_icdf05 = 0.26999994629410995d;
            var IG = new InverseGamma(0.5d, 0.42d);
            Assert.AreEqual(IG.Median, true_median, 0.0001d);
            Assert.AreEqual(IG.PDF(0.27d), true_pdf, 0.0001d);
            Assert.AreEqual(IG.CDF(0.27d), true_cdf, 0.0001d);
            Assert.AreEqual(IG.InverseCDF(IG.CDF(0.27d)), true_icdf05, 0.0001d);
        }

        /// <summary>
        /// Testing InverseGamma is being created.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(0.5,IG.Beta);
            Assert.AreEqual(2,IG.Alpha);

            var IG2 = new InverseGamma(2, 4);
            Assert.AreEqual(2,IG2.Beta);
            Assert.AreEqual(4, IG2.Alpha);
        }

        /// <summary>
        /// Checking inverse gamma distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var IG = new InverseGamma(double.NaN, double.NaN);
            Assert.IsFalse(IG.ParametersValid);

            var IG2 = new InverseGamma(double.PositiveInfinity, double.PositiveInfinity);
            Assert.IsFalse(IG2.ParametersValid);

            var IG3 = new InverseGamma(0, 0);
            Assert.IsFalse(IG3.ParametersValid);
        }

        /// <summary>
        /// Checking ParametersToString()
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var IG = new InverseGamma();
            Assert.AreEqual("Scale (β)",IG.ParametersToString[0, 0] );
            Assert.AreEqual("Shape (α)", IG.ParametersToString[1, 0]);
            Assert.AreEqual("0.5", IG.ParametersToString[0, 1]);
            Assert.AreEqual("2", IG.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(0.5, IG.Mean);

            var IG2 = new InverseGamma(1, 1);
            Assert.AreEqual(double.NaN,IG2.Mean);
        }

        /// <summary>
        /// Testing median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(0.2979, IG.Median, 1e-04);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(0.1666, IG.Mode,  1e-04);

            var IG2 = new InverseGamma(1, 1);
            Assert.AreEqual(0.5, IG2.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(double.NaN, IG.StandardDeviation);

            var IG2 = new InverseGamma(0.5, 3);
            Assert.AreEqual(0.25, IG2.StandardDeviation);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(double.NaN, IG.Skewness);

            var IG2 = new InverseGamma(0.5, 4);
            Assert.AreEqual(5.65685, IG2.Skewness,  1e-04);
        }

        /// <summary>
        /// Testing kurtosis with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(double.NaN, IG.Kurtosis);

            var IG2 = new InverseGamma(0.5, 5);
            Assert.AreEqual(45, IG2.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions are 0 and positive infinity respectively.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(0, IG.Minimum);
            Assert.AreEqual(double.PositiveInfinity, IG.Maximum);

            var IG2 = new InverseGamma(2, 2);
            Assert.AreEqual(0, IG2.Minimum);
            Assert.AreEqual(double.PositiveInfinity, IG2.Maximum);
        }

        /// <summary>
        /// Testing PDF method at different locations with varying parameters.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var IG = new InverseGamma(2,4);
            Assert.AreEqual(0, IG.PDF(-2));
            Assert.AreEqual(0.00057200, IG.PDF(5),  1e-07);
            Assert.AreEqual(1.74443, IG.PDF(0.42),  1e-04);

            var IG2 = new InverseGamma(0.42,2.4);
            Assert.AreEqual(double.NaN, IG2.PDF(0));
            Assert.AreEqual(1.48386, IG2.PDF(0.3),  1e-05);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(0, IG.CDF(-1));
            Assert.AreEqual(1, IG.CDF(double.PositiveInfinity));

            var IG2 = new InverseGamma(2, 2);
            Assert.AreEqual(0.73575, IG2.CDF(2), 1e-04);
        }

        /// <summary>
        /// Testing InverseCDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(0, IG.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity, IG.InverseCDF(1));

            var IG2 = new InverseGamma(2, 2);
            Assert.AreEqual(0.81993, IG2.InverseCDF(0.3), 1e-04);
        }
    }
}
