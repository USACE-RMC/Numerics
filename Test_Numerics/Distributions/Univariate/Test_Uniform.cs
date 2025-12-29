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
    /// Testing the Uniform distribution algorithm.
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
    public class Test_Uniform
    {

        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_UniformDist()
        {
            double true_mean = 0.76d;
            double true_median = 0.76d;
            double true_stdDev = Math.Sqrt(0.03853333333333335d);
            double true_pdf = 1.4705882352941173d;
            double true_cdf = 0.70588235294117641d;
            double true_icdf = 0.9d;
            var U = new Uniform(0.42d, 1.1d);
            Assert.AreEqual(true_mean, U.Mean, 0.0001d);
            Assert.AreEqual(true_median, U.Median, 0.0001d);
            Assert.AreEqual(true_stdDev, U.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_pdf, U.PDF(0.9d), 0.0001d);
            Assert.AreEqual(true_cdf, U.CDF(0.9d), 0.0001d);
            Assert.AreEqual(true_icdf, U.InverseCDF(true_cdf), 0.0001d);
        }

        /// <summary>
        /// Verified using R 'stats'
        /// </summary>
        [TestMethod]
        public void Test_Uniform_R()
        {
            var x = new double[] { 0.1, 0.25, 0.5, 0.75, 0.9 };
            var p = new double[5];
            var true_p = new double[] { 0.1, 0.25, 0.5, 0.75, 0.9 };
            var u = new Uniform(0, 1);
            for (int i = 0; i < 5; i++)
            {
                p[i] = u.CDF(x[i]);
                Assert.AreEqual(true_p[i], p[i], 1E-7);
            }
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(x[i], u.InverseCDF(p[i]), 1E-7);
            }
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.Min);
            Assert.AreEqual(1, U.Max);

            var U2 = new Uniform(2,10);
            Assert.AreEqual(2, U2.Min);
            Assert.AreEqual(10, U2.Max);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            var U = new Uniform(1, 0);
            Assert.IsFalse(U.ParametersValid);

            var U2 = new Uniform(double.NaN, 0);
            Assert.IsFalse(U2.ParametersValid);

            var U3 = new Uniform(0,double.NaN);
            Assert.IsFalse(U3.ParametersValid);

            var U4 = new Uniform(0,double.PositiveInfinity);
            Assert.IsFalse(U4.ParametersValid);

            var U5 = new Uniform(double.PositiveInfinity, 0);
            Assert.IsFalse(U5.ParametersValid);
        }

        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod]
        public void Test_ParametersToString()
        {
            var U = new Uniform();
            Assert.AreEqual("Min", U.ParametersToString[0, 0]);
            Assert.AreEqual("Max", U.ParametersToString[1, 0]);
            Assert.AreEqual("0", U.ParametersToString[0, 1]);
            Assert.AreEqual("1", U.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Uniform(2, 10);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            var U = new Uniform();
            Assert.AreEqual(0.5, U.Mean);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(6, U2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod]
        public void Test_Median()
        {
            var U = new Uniform();
            Assert.AreEqual(0.5, U.Median);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(6, U2.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod]
        public void Test_Mode()
        {
            var U = new Uniform();
            Assert.AreEqual(double.NaN, U.Mode);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(double.NaN, U2.Mode);
        }

        /// <summary>
        /// Testing Standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_StandardDeviation()
        {
            var U = new Uniform();
            Assert.AreEqual(0.288675, U.StandardDeviation, 1e-05);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(2.3094, U2.StandardDeviation, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod]
        public void Test_Skewness()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.Skewness);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(0, U2.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod]
        public void Test_Kurtosis()
        {
            var U = new Uniform();
            Assert.AreEqual(9d / 5d, U.Kurtosis);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(9d / 5d, U2.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod]
        public void Test_MinMax()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.Minimum);
            Assert.AreEqual(1, U.Maximum);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(2, U2.Minimum);
            Assert.AreEqual(10, U2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.PDF(-1));
            Assert.AreEqual(0, U.PDF(2));
            Assert.AreEqual(1, U.PDF(1));
        }

        /// <summary>
        /// Testing CDF.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.CDF(0));
            Assert.AreEqual(1, U.CDF(1));
            Assert.AreEqual(0.5, U.CDF(0.5));
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.InverseCDF(0));
            Assert.AreEqual(1, U.InverseCDF(1));
            Assert.AreEqual(0.3, U.InverseCDF(0.3));        
        }
    }
}
