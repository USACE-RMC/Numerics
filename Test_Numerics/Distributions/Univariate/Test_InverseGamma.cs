/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
            Assert.AreEqual(IG.Beta, 0.5);
            Assert.AreEqual(IG.Alpha, 2);

            var IG2 = new InverseGamma(2, 4);
            Assert.AreEqual(IG2.Beta, 2);
            Assert.AreEqual(IG2.Alpha, 4);
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
            Assert.AreEqual(IG.ParametersToString[0, 0], "Scale (β)");
            Assert.AreEqual(IG.ParametersToString[1, 0], "Shape (α)");
            Assert.AreEqual(IG.ParametersToString[0, 1], "0.5");
            Assert.AreEqual(IG.ParametersToString[1, 1], "2");
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.Mean, 0.5);

            var IG2 = new InverseGamma(1, 1);
            Assert.AreEqual(IG2.Mean, double.NaN);
        }

        /// <summary>
        /// Testing median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.Median, 0.2979,1e-04);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.Mode, 0.1666, 1e-04);

            var IG2 = new InverseGamma(1, 1);
            Assert.AreEqual(IG2.Mode, 0.5);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.StandardDeviation, double.NaN);

            var IG2 = new InverseGamma(0.5, 3);
            Assert.AreEqual(IG2.StandardDeviation, 0.25);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.Skewness, double.NaN);

            var IG2 = new InverseGamma(0.5, 4);
            Assert.AreEqual(IG2.Skewness, 5.65685, 1e-04);
        }

        /// <summary>
        /// Testing kurtosis with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.Kurtosis, double.NaN);

            var IG2 = new InverseGamma(0.5, 5);
            Assert.AreEqual(IG2.Kurtosis, 45);
        }

        /// <summary>
        /// Testing minimum and maximum functions are 0 and positive infinity respectively.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.Minimum, 0);
            Assert.AreEqual(IG.Maximum, double.PositiveInfinity);

            var IG2 = new InverseGamma(2, 2);
            Assert.AreEqual(IG2.Minimum, 0);
            Assert.AreEqual(IG2.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method at different locations with varying parameters.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var IG = new InverseGamma(2,4);
            Assert.AreEqual(IG.PDF(-2), 0);
            Assert.AreEqual(IG.PDF(5), 0.00057200, 1e-07);
            Assert.AreEqual(IG.PDF(0.42), 1.74443, 1e-04);

            var IG2 = new InverseGamma(0.42,2.4);
            Assert.AreEqual(IG2.PDF(0), double.NaN);
            Assert.AreEqual(IG2.PDF(0.3), 1.48386, 1e-05);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.CDF(-1), 0);
            Assert.AreEqual(IG.CDF(double.PositiveInfinity), 1);

            var IG2 = new InverseGamma(2, 2);
            Assert.AreEqual(IG2.CDF(2), 0.73575,1e-04);
        }

        /// <summary>
        /// Testing InverseCDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var IG = new InverseGamma();
            Assert.AreEqual(IG.InverseCDF(0), 0);
            Assert.AreEqual(IG.InverseCDF(1),double.PositiveInfinity);

            var IG2 = new InverseGamma(2, 2);
            Assert.AreEqual(IG2.InverseCDF(0.3), 0.81993,1e-04);
        }
    }
}
