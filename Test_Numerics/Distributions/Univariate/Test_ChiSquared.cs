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
    /// Testing the Chi-Squared distribution algorithm.
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
    /// 
    [TestClass]
    public class Test_ChiSquared
    {
        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_ChiSquaredDist()
        {
            double true_mean = 7d;
            double true_median = 6.345811195595612d;
            double true_stdDev = Math.Sqrt(14d);
            double true_pdf = 0.11388708001184455d;
            double true_cdf = 0.49139966433823956d;
            double true_icdf = 6.27d;
            var CHI = new ChiSquared(7);
            Assert.AreEqual(CHI.Mean, true_mean, 0.0001d);
            Assert.AreEqual(CHI.Median, true_median, 0.0001d);
            Assert.AreEqual(CHI.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(CHI.PDF(6.27d), true_pdf, 0.0001d);
            Assert.AreEqual(CHI.CDF(6.27d), true_cdf, 0.0001d);
            Assert.AreEqual(CHI.InverseCDF(true_cdf), true_icdf, 0.0001d);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. See if chi squared is being created.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(1, x.DegreesOfFreedom);

            var x2 = new ChiSquared(3);
            Assert.AreEqual(3, x2.DegreesOfFreedom);
        }

        /// <summary>
        /// Verifying where Chi Squared fails.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var x = new ChiSquared(0);
            Assert.IsFalse(x.ParametersValid);

            var x2 = new ChiSquared(-1);
            Assert.IsFalse(x2.ParametersValid);

            var x3 = new ChiSquared(-100);
            Assert.IsFalse(x3.ParametersValid);
        }

        /// <summary>
        /// Verifying the parameters are being converted to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual("Degrees of Freedom (ν)", x.ParametersToString[0,0]);
            Assert.AreEqual("1", x.ParametersToString[0,1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new ChiSquared(2);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Checking mean is equal to the degrees of freedom.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(1, x.Mean);

            var x2 = new ChiSquared(2);
            Assert.AreEqual(2, x2.Mean);
        }

        /// <summary>
        /// Validating the median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var x = new ChiSquared(2);
            var approx_median = x.DegreesOfFreedom * Math.Pow(1d - 2d / (9d * x.DegreesOfFreedom), 3d);
            Assert.AreEqual(x.Median, approx_median, 1E-1);
        }

        /// <summary>
        /// Validating the mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(0, x.Mode);

            var x2 = new ChiSquared(5);
            Assert.AreEqual(3, x2.Mode);
        }

        /// <summary>
        /// Validating standard deviation function.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(Math.Sqrt(2), x.StandardDeviation);

            var x2 = new ChiSquared(2);
            Assert.AreEqual(2,x2.StandardDeviation);
        }

        /// <summary>
        /// Validating the skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var x = new ChiSquared(2);
            Assert.AreEqual(2, x.Skewness);

            var x2 = new ChiSquared(8);
            Assert.AreEqual(1, x2.Skewness);

            var x3 = new ChiSquared(32);
            Assert.AreEqual(0.5, x3.Skewness);
        }

        /// <summary>
        /// Validating Kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(15, x.Kurtosis);

            var x2 = new ChiSquared(3);
            Assert.AreEqual(7,x2.Kurtosis);
        }

        /// <summary>
        /// Checking minimum of Chi-Squared Distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(0, x.Minimum);
        }

        /// <summary>
        /// Checking maximum of Chi-Squared Distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(double.PositiveInfinity, x.Maximum);
        }

        /// <summary>
        /// Checking PDF function.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(1.2000389,x.PDF(0.1),1e-04);
            Assert.AreEqual(0.2419707, x.PDF(1), 1e-04);
            Assert.AreEqual(0.0108747, x.PDF(5.5), 1e-04);
            Assert.AreEqual(0, x.PDF(double.PositiveInfinity));

            var x2 = new ChiSquared(2);
            Assert.AreEqual(0.5,x2.PDF(0));
            Assert.AreEqual(0.4756147, x2.PDF(0.1), 1e-04);
            Assert.AreEqual(0.3032653, x2.PDF(1), 1e-04);
            Assert.AreEqual(0.0319639, x2.PDF(5.5), 1e-04);
        }

        /// <summary>
        /// Validating CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(0.2481703, x.CDF(0.1), 1e-04);
            Assert.AreEqual(0.682689, x.CDF(1), 1e-04);
            Assert.AreEqual(0.9809835, x.CDF(5.5), 1e-04);
            Assert.AreEqual(1, x.CDF(110.1));

            var x2 = new ChiSquared(2);
            Assert.AreEqual(0, x2.CDF(0));
            Assert.AreEqual(0.04877057, x2.CDF(0.1), 1e-04);
            Assert.AreEqual(0.3934693, x2.CDF(1), 1e-04);
            Assert.AreEqual(0.9360721, x2.CDF(5.5), 1e-04);
        }

        /// <summary>
        /// Validating Inverse CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var x = new ChiSquared(1);
            Assert.AreEqual(x.InverseCDF(0.24817036595415071751), 0.09999,1e-04);
            Assert.AreEqual(x.InverseCDF(0.68268949213708589717), 1, 1e-04);
            Assert.AreEqual(x.InverseCDF(0.9809835), 5.5, 1e-04);

            var x2 = new ChiSquared(2);
            Assert.AreEqual(x2.InverseCDF(0), 0);
            Assert.AreEqual(x2.InverseCDF(0.04877057), 0.1,1e-04);
            Assert.AreEqual(x2.InverseCDF(0.3934693), 1,1e-04);
            Assert.AreEqual(x2.InverseCDF(0.9360721), 5.5, 1e-04);
        }
    }
}
