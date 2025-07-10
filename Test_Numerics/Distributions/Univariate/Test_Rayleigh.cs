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
    /// Testing the Rayleigh distribution algorithm.
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
    public class Test_Rayleigh
    {

        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_RayleighDist()
        {
            double true_mean = 0.52639193767251d;
            double true_median = 0.49451220943852386d;
            double true_stdDev = Math.Sqrt(0.075711527953380237d);
            double true_pdf = 0.030681905868831811d;
            double true_cdf = 0.99613407986052716d;
            double true_icdf = 1.4d;
            var R = new Rayleigh(0.42);
            Assert.AreEqual(R.Mean, true_mean, 0.0001d);
            Assert.AreEqual(R.Median, true_median, 0.0001d);
            Assert.AreEqual(R.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(R.PDF(1.4d), true_pdf, 0.0001d);
            Assert.AreEqual(R.CDF(1.4d), true_cdf, 0.0001d);
            Assert.AreEqual(R.InverseCDF(true_cdf), true_icdf, 0.0001d);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.Sigma, 10);

            var R2 = new Rayleigh(2);
            Assert.AreEqual(R2.Sigma, 2);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var R = new Rayleigh(double.NaN);
            Assert.IsFalse(R.ParametersValid);

            var R2 = new Rayleigh(double.PositiveInfinity);
            Assert.IsFalse(R2.ParametersValid);

            var R3 = new Rayleigh(0);
            Assert.IsFalse(R3.ParametersValid);
        }

        /// <summary>
        /// Testing Parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.ParametersToString[0, 0], "Scale (σ)");
            Assert.AreEqual(R.ParametersToString[0, 1], "10");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Rayleigh(1);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.Mean, 12.53314, 1e-04);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(R2.Mean, 1.25331, 1e-04);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.Median, 11.7741, 1e-04);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(R2.Median, 1.1774, 1e-04);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.Mode, 10);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(R2.Mode, 1);
        }
        
        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.StandardDeviation, 6.55136, 1e-05);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(R2.StandardDeviation, 0.65513, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.Skewness, 0.63111, 1e-04);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(R2.Skewness, 0.63111, 1e-04);
        }

        /// <summary>
        /// Testing kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.Kurtosis, 3.24508,1e-05);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(R2.Kurtosis, 3.24508,1e-05);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.Minimum, 0);
            Assert.AreEqual(R.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.PDF(-1), 0);
            Assert.AreEqual(R.PDF(1), 9.9501e-03, 1e-06);

            var R2 = new Rayleigh(1);
            Assert.AreEqual(R.PDF(2), 0.019603, 1e-05);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.CDF(-1), 0);
            Assert.AreEqual(R.CDF(1), 4.9875e-03,1e-04);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var R = new Rayleigh();
            Assert.AreEqual(R.InverseCDF(0), 0);
            Assert.AreEqual(R.InverseCDF(1), double.PositiveInfinity);
            Assert.AreEqual(R.InverseCDF(0.4), 10.1076, 1e-04);
        }
    }
}
