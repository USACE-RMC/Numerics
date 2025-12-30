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
    /// Testing the Truncated Normal distribution algorithm.
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
    public class Test_TruncatedNormal
    {

        /// <summary>
        /// This method was verified against the R package "truncnorm"
        /// </summary>
        [TestMethod]
        public void Test_TruncatedNormalDist()
        {
            var tn = new TruncatedNormal(2, 1, 1.10, 2.11);
            var d = tn.PDF(1.5);
            var p = tn.CDF(1.5);
            var q = tn.InverseCDF(p);

            Assert.AreEqual(0.9786791, d,  1E-5);
            Assert.AreEqual(0.3460251, p,  1E-5);
            Assert.AreEqual(1.5, q,  1E-5);

            tn = new TruncatedNormal(10, 3, 8, 25);
            d = tn.PDF(12.75);
            p = tn.CDF(12.75);
            q = tn.InverseCDF(p);

            Assert.AreEqual(0.1168717, d ,1E-5);
            Assert.AreEqual(0.7596566, p, 1E-5);
            Assert.AreEqual(12.75, q, 1E-5);

            tn = new TruncatedNormal(0, 3, 0, 9);
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
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.5,tn.Mu );
            Assert.AreEqual(0.2,tn.Sigma );
            Assert.AreEqual(0,tn.Min);
            Assert.AreEqual(1,tn.Max);

            var tn2 = new TruncatedNormal(1, 1, 1, 2);
            Assert.AreEqual(1,tn2.Mu);
            Assert.AreEqual(1,tn2.Sigma);
            Assert.AreEqual(1,tn2.Min);
            Assert.AreEqual(2,tn2.Max);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var tn = new TruncatedNormal(double.NaN, double.NaN, double.NaN, double.NaN);
            Assert.IsFalse(tn.ParametersValid);

            var tn2 = new TruncatedNormal(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
            Assert.IsFalse(tn2.ParametersValid);

            var tn3 = new TruncatedNormal(0, -1, -1, 0);
            Assert.IsFalse(tn3.ParametersValid);

            var tn4 = new TruncatedNormal(1, 1, 1, 0);
            Assert.IsFalse(tn4.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual("Mean (µ)",tn.ParametersToString[0, 0]);
            Assert.AreEqual("Std Dev (σ)",tn.ParametersToString[1, 0]);
            Assert.AreEqual("Min",tn.ParametersToString[2, 0] );
            Assert.AreEqual("Max", tn.ParametersToString[3, 0]  );
            Assert.AreEqual("0.5",tn.ParametersToString[0, 1]);
            Assert.AreEqual("0.2", tn.ParametersToString[1, 1]);
            Assert.AreEqual("0", tn.ParametersToString[2, 1]);
            Assert.AreEqual("1", tn.ParametersToString[3, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new TruncatedNormal(10, 3, 8, 25);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.5,tn.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.5, tn.Median);
        }

        /// <summary>
        /// Testing mode
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.5, tn.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0.19091, tn.StandardDeviation, 1e-05);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0, tn.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(2.62422, tn.Kurtosis,  1e-04);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var tn = new TruncatedNormal();
            Assert.AreEqual(0, tn.Minimum);
            Assert.AreEqual(1, tn.Maximum);
        }
    }
}
