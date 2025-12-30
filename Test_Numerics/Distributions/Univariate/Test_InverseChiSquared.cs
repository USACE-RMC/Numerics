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
            double true_mean = 0.2;
            double true_median = 6.345811068141737d;
            double true_pdf = 0.0000063457380298844403d;
            double true_cdf = 0.50860033566176044d;
            double true_icdf = 6.27d;
            var IX = new InverseChiSquared(7, (1d / 7d));
            double pdf = IX.PDF(6.27d);
            double cdf = IX.CDF(6.27d);
            double icdf = IX.InverseCDF(cdf);
            Assert.AreEqual(IX.Mean, true_mean, 0.0001d);
            Assert.AreEqual(IX.Median, true_median, 0.0001d);
            Assert.AreEqual(IX.PDF(6.27d), true_pdf, 0.0001d);
            Assert.AreEqual(IX.CDF(6.27d), true_cdf, 0.0001d);
            Assert.AreEqual(IX.InverseCDF(IX.CDF(6.27d)), true_icdf, 0.0001d);
        }

        /// <summary>
        /// Checking Inverse Chi-Squared can be created with inputs.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(10,IX.DegreesOfFreedom);
            Assert.AreEqual(1,IX.Sigma);

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
            var IX = new InverseChiSquared();
            Assert.AreEqual(0.93418, IX.Median,1e-04);

            var IX2 = new InverseChiSquared(7, 1);
            Assert.AreEqual(0.906544, IX2.Median,  1e-04);
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
            Assert.AreEqual(double.NaN,IX2.StandardDeviation );
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
            Assert.AreEqual(double.NaN,IX2.Skewness);
        }

        /// <summary>
        /// Testing kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(45,IX.Kurtosis);

            var IX2 = new InverseChiSquared(2,2);
            Assert.AreEqual(double.NaN,IX2.Kurtosis );
        }

        /// <summary>
        /// Testing Minimum and Maximum functions are 0 and positive infinity respectively.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(0,IX.Minimum);
            Assert.AreEqual(double.PositiveInfinity,IX.Maximum );
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
            var IX = new InverseChiSquared(7,1);
            Assert.AreEqual(1.1184e-05, IX.CDF(5),  1e-09);
        }

        /// <summary>
        /// Testing the Inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var IX = new InverseChiSquared();
            Assert.AreEqual(0,IX.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity,IX.InverseCDF(1));
            Assert.AreEqual(1.17807, IX.InverseCDF(0.3), 1e-04);
        }
    }
}
