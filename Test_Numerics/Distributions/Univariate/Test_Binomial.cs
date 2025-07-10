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

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

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
    }
}
