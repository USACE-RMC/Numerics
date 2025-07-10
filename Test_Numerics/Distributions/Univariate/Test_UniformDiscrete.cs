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
    /// Testing the Discrete Uniform distribution algorithm.
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
    public class Test_UniformDiscrete
    {

        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_UniformDiscreteDist()
        {
            double true_mean = 4.0d;
            double true_median = 4.0d;
            double true_stdDev = Math.Sqrt(1.3333333333333333d);
            int true_skew = 0;
            double true_kurt = 1.7d;
            double true_pdf = 0.2d;
            double true_cdf = 0.2d;
            double true_icdf05 = 2.0d;
            double true_icdf95 = 6.0d;
            var U = new UniformDiscrete(2d, 6d);
            Assert.AreEqual(U.Mean, true_mean, 0.0001d);
            Assert.AreEqual(U.Median, true_median, 0.0001d);
            Assert.AreEqual(U.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(U.Skewness, true_skew, 0.0001d);
            Assert.AreEqual(U.Kurtosis, true_kurt, 0.0001d);
            Assert.AreEqual(U.PDF(4.0d), true_pdf, 0.0001d);
            Assert.AreEqual(U.CDF(2.0d), true_cdf, 0.0001d);
            Assert.AreEqual(U.InverseCDF(0.17d), true_icdf05, 0.0001d);
            Assert.AreEqual(U.InverseCDF(0.87d), true_icdf95, 0.0001d);
        }
        
        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.Min, 0);
            Assert.AreEqual(U.Max, 1);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(U2.Min, 2);
            Assert.AreEqual(U2.Max, 10);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var U = new UniformDiscrete(1, 0);
            Assert.IsFalse(U.ParametersValid);

            var U2 = new UniformDiscrete(double.NaN, 0);
            Assert.IsFalse(U2.ParametersValid);

            var U3 = new UniformDiscrete(0, double.NaN);
            Assert.IsFalse(U3.ParametersValid);

            var U4 = new UniformDiscrete(0, double.PositiveInfinity);
            Assert.IsFalse(U4.ParametersValid);

            var U5 = new UniformDiscrete(double.PositiveInfinity, 0);
            Assert.IsFalse(U5.ParametersValid);
        }
        
        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.ParametersToString[0, 0], "Min");
            Assert.AreEqual(U.ParametersToString[1, 0], "Max");
            Assert.AreEqual(U.ParametersToString[0, 1], "0");
            Assert.AreEqual(U.ParametersToString[1, 1], "1");
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.Mean, 0.5);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(U2.Mean, 6);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod]
        public void Test_Median()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.Median, 0.5);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(U2.Median, 6);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod]
        public void Test_Mode()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.Mode, double.NaN);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(U2.Mode, double.NaN);
        }

        /// <summary>
        /// Testing Standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_StandardDeviation()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.StandardDeviation, 0.288675, 1e-05);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(U2.StandardDeviation, 2.3094, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod]
        public void Test_Skewness()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.Skewness, 0);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(U2.Skewness, 0);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod]
        public void Test_Kurtosis()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.Kurtosis, 1);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(U2.Kurtosis, 1.77);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod]
        public void Test_MinMax()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.Minimum, 0);
            Assert.AreEqual(U.Maximum, 1);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(U2.Minimum, 2);
            Assert.AreEqual(U2.Maximum, 10);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.PDF(-1), 0);
            Assert.AreEqual(U.PDF(2), 0);
            Assert.AreEqual(U.PDF(1), 0.5);
        }

        /// <summary>
        /// Testing CDF.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.CDF(0), 0.5);
            Assert.AreEqual(U.CDF(1), 1);
            Assert.AreEqual(U.CDF(0.5), 0.75);
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(U.InverseCDF(0), 0);
            Assert.AreEqual(U.InverseCDF(1), 1);
            Assert.AreEqual(U.InverseCDF(0.3), 0);
        }
    }

}
