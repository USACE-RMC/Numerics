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
    /// Testing the Poisson distribution algorithm.
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
    public class Test_Poisson
    {

        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_PoissonDist()
        {
            double true_mean = 4.2d;
            double true_mode = 4.0d;
            double true_median = 4.0d;
            double true_stdDev = Math.Sqrt(4.2d);
            double true_skew = 0.488d;
            double true_kurt = 3.2381d;
            double true_pdf = 0.19442365170822165d;
            double true_cdf = 0.58982702131057763d;
            double true_icdf05 = 1.0d;
            double true_icdf95 = 8.0d;
            var P = new Poisson(4.2d);
            Assert.AreEqual(P.Mean, true_mean, 0.0001d);
            Assert.AreEqual(P.Median, true_median, 0.0001d);
            Assert.AreEqual(P.Mode, true_mode, 0.0001d);
            Assert.AreEqual(P.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(P.Skewness, true_skew, 0.0001d);
            Assert.AreEqual(P.Kurtosis, true_kurt, 0.0001d);
            Assert.AreEqual(P.PDF(4.0d), true_pdf, 0.0001d);
            Assert.AreEqual(P.CDF(4.0d), true_cdf, 0.0001d);
            Assert.AreEqual(P.InverseCDF(0.05d), true_icdf05, 0.0001d);
            Assert.AreEqual(P.InverseCDF(0.95d), true_icdf95, 0.0001d);
        }

        /// <summary>
        /// Verifying that input parameters can create this distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var P = new Poisson();
            Assert.AreEqual(P.Lambda, 1);

            var P2 = new Poisson(10);
            Assert.AreEqual(P2.Lambda, 10);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var P = new Poisson(double.NaN);
            Assert.IsFalse(P.ParametersValid);

            var P2 = new Poisson(double.PositiveInfinity);
            Assert.IsFalse(P2.ParametersValid);

            var P3 = new Poisson(0);
            Assert.IsFalse(P3.ParametersValid);
        }

        /// <summary>
        /// Testing ParameterToString().
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var P = new Poisson();
            Assert.AreEqual(P.ParametersToString[0, 0], "Rate (λ)");
            Assert.AreEqual(P.ParametersToString[0, 1], "1");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Poisson(10);
            var mom = dist.CentralMoments(1000);
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
            var P = new Poisson();
            Assert.AreEqual(P.Mean, 1);

            var P2 = new Poisson(10);
            Assert.AreEqual(P2.Mean, 10);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var P = new Poisson();
            Assert.AreEqual(P.Median, 1, 1E-4);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var P = new Poisson();
            Assert.AreEqual(P.Mode, 1);

            var P2 = new Poisson(2.4);
            Assert.AreEqual(P2.Mode, 2);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var P = new Poisson();
            Assert.AreEqual(P.StandardDeviation, 1);

            var P2 = new Poisson(4);
            Assert.AreEqual(P2.StandardDeviation, 2);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var P = new Poisson();
            Assert.AreEqual(P.Skewness, 1);

            var P2 = new Poisson(4);
            Assert.AreEqual(P2.Skewness, 0.5);
        }

        /// <summary>
        /// Testing kurtosis of this distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var P = new Poisson();
            Assert.AreEqual(P.Kurtosis, 4);

            var P2 = new Poisson(4);
            Assert.AreEqual(P2.Kurtosis, 3.25);
        }

        /// <summary>
        /// Testing minimum and maximum
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var P = new Poisson();
            Assert.AreEqual(P.Minimum, 0);
            Assert.AreEqual(P.Maximum,double.PositiveInfinity);

            var P2 = new Poisson(4);
            Assert.AreEqual(P2.Minimum, 0);
            Assert.AreEqual(P2.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var P = new Poisson(1.5);
            Assert.AreEqual(P.PDF(1), 0.33469, 1e-04);
            Assert.AreEqual(P.PDF(10), 0.00000354, 1e-08);

            var P2 = new Poisson(5.4);
            Assert.AreEqual(P2.PDF(1), 0.024389, 1e-05);
            Assert.AreEqual(P2.PDF(10), 0.02624, 1e-05);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var P = new Poisson(1.5);
            Assert.AreEqual(P.CDF(1), 0.55782, 1e-04);
            Assert.AreEqual(P.CDF(10), 0.999999, 1e-06);

            var P2 = new Poisson(10.8);
            Assert.AreEqual(P2.CDF(1), 0.00024, 1e-05);
            Assert.AreEqual(P2.CDF(10), 0.483969, 1e-05);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var P = new Poisson();
            Assert.AreEqual(P.InverseCDF(0), 0);
            Assert.AreEqual(P.InverseCDF(1), double.PositiveInfinity);
        }
    }
}
