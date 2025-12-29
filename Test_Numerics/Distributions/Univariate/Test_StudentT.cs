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
    /// Testing the Student T distribution algorithm.
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
    public class Test_StudentT
    {
        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_StudentT_PDF()
        {
            var t = new StudentT(4.2d);
            double pdf = t.PDF(1.4d);
            double result = 0.138377537135553d;
            Assert.AreEqual(result, pdf, 1E-10);
            t = new StudentT(2.5d, 0.5d, 4.2d);
            pdf = t.PDF(1.4d);
            result = 0.0516476521260042d;
            Assert.AreEqual(result, pdf, 1E-10);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_StudentT_CDF()
        {
            var t = new StudentT(4.2d);
            double cdf = t.CDF(1.4d);
            double result = 0.882949686336585d;
            Assert.AreEqual(result, cdf, 1E-10);
            t = new StudentT(2.5d, 0.5d, 4.2d);
            cdf = t.CDF(1.4d);
            result = 0.0463263350898173d;
            Assert.AreEqual(result, cdf, 1E-10);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_StudentT_InverseCDF()
        {
            var t = new StudentT(4.2d);
            double cdf = t.CDF(1.4d);
            double invcdf = t.InverseCDF(cdf);
            double result = 1.4d;
            Assert.AreEqual(result, invcdf, 1E-2);
            t = new StudentT(2.5d, 0.5d, 4.2d);
            cdf = t.CDF(1.4d);
            invcdf = t.InverseCDF(cdf);
            result = 1.4d;
            Assert.AreEqual(result, invcdf, 1E-2);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var t = new StudentT();
            Assert.AreEqual(0, t.Mu);
            Assert.AreEqual(1, t.Sigma);
            Assert.AreEqual(10, t.DegreesOfFreedom);

            var t2 = new StudentT(10, 10, 10);
            Assert.AreEqual(10, t2.Mu);
            Assert.AreEqual(10, t2.Sigma);
            Assert.AreEqual(10, t2.DegreesOfFreedom);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var t = new StudentT(double.NaN, double.NaN, 1);
            Assert.IsFalse(t.ParametersValid);

            var t2 = new StudentT(double.PositiveInfinity, double.PositiveInfinity, 1);
            Assert.IsFalse(t2.ParametersValid);

            var t3 = new StudentT(1, 1, 0);
            Assert.IsFalse(t3.ParametersValid);
        }
        
        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var t = new StudentT();
            Assert.AreEqual(t.ParametersToString[0, 0], "Location (µ)");
            Assert.AreEqual(t.ParametersToString[1, 0], "Scale (σ)");
            Assert.AreEqual(t.ParametersToString[2, 0], "Degrees of Freedom (ν)");
            Assert.AreEqual(t.ParametersToString[0, 1], "0");
            Assert.AreEqual(t.ParametersToString[1,1],"1");
            Assert.AreEqual(t.ParametersToString[2, 1], "10");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new StudentT(10, 1, 100);
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
            var t = new StudentT();
            Assert.AreEqual(0, t.Mean);

            var t2 = new StudentT(1, 1, 1);
            Assert.AreEqual(double.NaN, t2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var t = new StudentT();
            Assert.AreEqual(0, t.Median);

            var t2 = new StudentT(1, 1, 1);
            Assert.AreEqual(1, t2.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var t = new StudentT();
            Assert.AreEqual(0, t.Mode);

            var t2 = new StudentT(1,1,1);
            Assert.AreEqual(1, t2.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var t = new StudentT();
            Assert.AreEqual(1.11803, t.StandardDeviation, 1e-04);

            var t2 = new StudentT(1, 1, 2);
            Assert.AreEqual(double.PositiveInfinity, t2.StandardDeviation);

            var t3 = new StudentT(1, 1, 1);
            Assert.AreEqual(double.NaN, t3.StandardDeviation);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var t = new StudentT();
            Assert.AreEqual(0, t.Skewness);

            var t2 = new StudentT(1, 1, 1);
            Assert.AreEqual(double.NaN, t2.Skewness);
        }

        /// <summary>
        /// Testing kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var t = new StudentT();
            Assert.AreEqual(4, t.Kurtosis);

            var t2 = new StudentT(1, 1, 4);
            Assert.AreEqual(double.PositiveInfinity, t2.Kurtosis);

            var t3 = new StudentT(1, 1, 2);
            Assert.AreEqual(double.NaN, t3.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var t = new StudentT();
            Assert.AreEqual(double.NegativeInfinity, t.Minimum);
            Assert.AreEqual(double.PositiveInfinity, t.Maximum);
        }
    }
}
