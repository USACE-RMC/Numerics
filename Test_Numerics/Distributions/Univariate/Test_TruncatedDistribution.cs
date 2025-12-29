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
    /// Testing the general Truncated distribution class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list> 
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_TruncatedDistribution
    {

        /// <summary>
        /// This method was verified against the R package "truncnorm"
        /// </summary>
        [TestMethod]
        public void Test_TruncatedNormalDist()
        {
            var tn = new TruncatedDistribution(new Normal(2, 1), 1.10, 2.11);
            var d = tn.PDF(1.5);
            var p = tn.CDF(1.5);
            var q = tn.InverseCDF(p);

            Assert.AreEqual(0.9786791, d, 1E-5);
            Assert.AreEqual(0.3460251, p, 1E-5);
            Assert.AreEqual(1.5, q, 1E-5);

            tn = new TruncatedDistribution(new Normal(10, 3), 8, 25);
            d = tn.PDF(12.75);
            p = tn.CDF(12.75);
            q = tn.InverseCDF(p);

            Assert.AreEqual(0.1168717, d, 1E-5);
            Assert.AreEqual(0.7596566, p, 1E-5);
            Assert.AreEqual(12.75, q, 1E-5);

            tn = new TruncatedDistribution(new Normal(0, 3), 0, 9);
            d = tn.PDF(4.5);
            p = tn.CDF(4.5);
            q = tn.InverseCDF(p);

            Assert.AreEqual(0.08657881, d, 1E-5);
            Assert.AreEqual(0.868731, p, 1E-5);
            Assert.AreEqual(4.5, q, 1E-5);

        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(((Normal)tn.BaseDistribution).Mu, 0.5);
            Assert.AreEqual(((Normal)tn.BaseDistribution).Sigma, 0.2);
            Assert.AreEqual(tn.Min, 0);
            Assert.AreEqual(tn.Max, 1);

            var tn2 = new TruncatedDistribution(new Normal(1, 1), 1, 2);
            Assert.AreEqual(((Normal)tn2.BaseDistribution).Mu, 1);
            Assert.AreEqual(((Normal)tn2.BaseDistribution).Sigma, 1);
            Assert.AreEqual(tn2.Min, 1);
            Assert.AreEqual(tn2.Max, 2);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(tn.ParametersToString[0, 0], "Mean (µ)");
            Assert.AreEqual(tn.ParametersToString[1, 0], "Std Dev (σ)");
            Assert.AreEqual(tn.ParametersToString[2, 0], "Min");
            Assert.AreEqual(tn.ParametersToString[3, 0], "Max");
            Assert.AreEqual(tn.ParametersToString[0, 1], "0.5");
            Assert.AreEqual(tn.ParametersToString[1, 1], "0.2");
            Assert.AreEqual(tn.ParametersToString[2, 1], "0");
            Assert.AreEqual(tn.ParametersToString[3, 1], "1");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new TruncatedDistribution(new Normal(10, 3), 8, 25);
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
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(tn.Mean, 0.5, 1e-4);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(tn.Median, 0.5, 1e-4);
        }

        /// <summary>
        /// Testing mode
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(tn.Mode, 0.5, 1e-4);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(tn.StandardDeviation, 0.19091, 1e-4);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(tn.Skewness, 0, 1E-4);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(tn.Kurtosis, 2.62422, 1e-04);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var tn = new TruncatedDistribution(new Normal(0.5, 0.2), 0, 1);
            Assert.AreEqual(tn.Minimum, 0);
            Assert.AreEqual(tn.Maximum, 1);
        }
    }
}
