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
using System;
using System.Diagnostics;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Noncentral T distribution algorithm.
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
    public class Test_NoncentralT
    {
        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_PDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double pdf = t.PDF(1.4d);
            double result = 0.23552141805184526d;
            Assert.AreEqual(pdf, result, 1E-6);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_CDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double cdf = t.CDF(1.4d);
            double result = 0.15955740661144721d;
            Assert.AreEqual(cdf, result, 1E-6);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_InverseCDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double cdf = t.CDF(1.4d);
            double invcdf = t.InverseCDF(cdf);
            double result = 1.4d;
            Assert.AreEqual(invcdf, result, 1E-6);
            var table = new[,] { { 3.0d, 0.0d, 1d, 0.89758361765043326d }, { 3.0d, 0.0d, 2d, 0.9522670169d }, { 3.0d, 0.0d, 3d, 0.97116555718878128d }, { 3.0d, 0.5d, 1d, 0.8231218864d }, { 3.0d, 0.5d, 2d, 0.904902151d }, { 3.0d, 0.5d, 3d, 0.9363471834d }, { 3.0d, 1.0d, 1d, 0.7301025986d }, { 3.0d, 1.0d, 2d, 0.8335594263d }, { 3.0d, 1.0d, 3d, 0.8774010255d }, { 3.0d, 2.0d, 1d, 0.5248571617d }, { 3.0d, 2.0d, 2d, 0.6293856597d }, { 3.0d, 2.0d, 3d, 0.6800271741d }, { 3.0d, 4.0d, 1d, 0.20590131975d }, { 3.0d, 4.0d, 2d, 0.2112148916d }, { 3.0d, 4.0d, 3d, 0.2074730718d }, { 15.0d, 7.0d, 15d, 0.9981130072d }, { 15.0d, 7.0d, 20d, 0.999487385d }, { 15.0d, 7.0d, 25d, 0.9998391562d }, { 0.05d, 1.0d, 1d, 0.168610566972d }, { 0.05d, 1.0d, 2d, 0.16967950985d }, { 0.05d, 1.0d, 3d, 0.1701041003d }, { 4.0d, 2.0d, 10d, 0.9247683363d }, { 4.0d, 3.0d, 10d, 0.7483139269d }, { 4.0d, 4.0d, 10d, 0.4659802096d }, { 5.0d, 2.0d, 10d, 0.9761872541d }, { 5.0d, 3.0d, 10d, 0.8979689357d }, { 5.0d, 4.0d, 10d, 0.7181904627d }, { 6.0d, 2.0d, 10d, 0.9923658945d }, { 6.0d, 3.0d, 10d, 0.9610341649d }, { 6.0d, 4.0d, 10d, 0.868800735d } };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double x = table[i, 0];
                double delta = table[i, 1];
                double degF = table[i, 2];
                var target = new NoncentralT(degF, delta);
                double expected = table[i, 3];
                double actual = target.CDF(x);
                Assert.AreEqual(expected, actual, 0.000001d);
                double expectedX = target.InverseCDF(actual);
                Assert.AreEqual(expectedX, x, 0.000001d);
            }

        }

        /// <summary>
        /// Verifying input parameters can create Noncentral T distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.DegreesOfFreedom, 10);
            Assert.AreEqual(t.Noncentrality, 0);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(t2.DegreesOfFreedom, 1);
            Assert.AreEqual(t2.Noncentrality, 1);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var t = new NoncentralT(0, 1);
            Assert.IsFalse(t.ParametersValid);

            var t2 = new NoncentralT(1,double.PositiveInfinity);
            Assert.IsFalse(t2.ParametersValid);

            var t3 = new NoncentralT(1,double.NaN);
            Assert.IsFalse(t3.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.ParametersToString[0, 0], "Degrees of Freedom (ν)");
            Assert.AreEqual(t.ParametersToString[1, 0], "Noncentrality (μ)");
            Assert.AreEqual(t.ParametersToString[0, 1], "10");
            Assert.AreEqual(t.ParametersToString[1, 1], "0");
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.Mean, 0);

            var t2 = new NoncentralT(0, 1);
            Assert.AreEqual(t2.Mean, double.NaN);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.Median, 0,1e-04);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(t2.Median, 1.3202,1e-04);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.Mode, 0, 1E-4);

            var t3 = new NoncentralT(10, 1);
            Assert.AreEqual(t3.Mode, 0.9329, 1e-04);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.StandardDeviation,1.1180,1e-04);

            var t2 = new NoncentralT(1, 0);
            Assert.AreEqual(t2.StandardDeviation, double.NaN);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.Skewness, 0.0, 1E-4);
        }

        /// <summary>
        /// Testing Kurtosis
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.Kurtosis, 4.0, 1E-4);
        }

        /// <summary>
        /// Testing min and max functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.Minimum,double.NegativeInfinity);
            Assert.AreEqual(t.Maximum,double.PositiveInfinity);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(t2.Minimum, double.NegativeInfinity);
            Assert.AreEqual(t2.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.PDF(0), 0.38910,1e-04);
            Assert.AreEqual(t.PDF(1),0.23036,1e-04);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.CDF(1), 0.82955,1e-04);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(t.InverseCDF(0), double.NegativeInfinity);
            Assert.AreEqual(t.InverseCDF(1),double.PositiveInfinity);
            Assert.AreEqual(t.InverseCDF(0.4), -0.26018,1e-04);
        }
    }
}
