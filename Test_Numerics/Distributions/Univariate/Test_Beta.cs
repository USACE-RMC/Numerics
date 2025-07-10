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
    /// Testing the Beta distribution algorithm.
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
    /// <list type="bullet">
    /// <item> <see href = "https://en.wikipedia.org/wiki/Beta_distribution" /></item>
    /// <item> <see href="https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests"/></item>
    /// </list>
    /// </remarks>
    [TestClass]
    public class Test_Beta
    {
        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_BetaDist()
        {
            double true_mean = 0.21105527638190955d;
            double true_mode = 0.0d;
            double true_median = 0.11577711097114812d;
            double true_stdDev = Math.Sqrt(0.055689279830523512d);
            double true_skew = 1.2275d;
            double true_kurt = 3.6048d;
            double true_pdf = 0.94644031936694828d;
            double true_cdf = 0.69358638272337991d;
            double true_icdf = 0.27d;
            double true_icdf05 = 0.000459d;
            double true_icdf95 = 0.7238d;
            var B = new BetaDistribution(0.42d, 1.57d);
            Assert.AreEqual(B.Mean, true_mean, 0.0001d);
            Assert.AreEqual(B.Median, true_median, 0.0001d);
            Assert.AreEqual(B.Mode, true_mode, 0.0001d);
            Assert.AreEqual(B.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(B.Skewness, true_skew, 0.0001d);
            Assert.AreEqual(B.Kurtosis, true_kurt, 0.0001d);
            Assert.AreEqual(B.PDF(0.27d), true_pdf, 0.0001d);
            Assert.AreEqual(B.CDF(0.27d), true_cdf, 0.0001d);
            Assert.AreEqual(B.InverseCDF(B.CDF(0.27d)), true_icdf, 0.0001d);
            Assert.AreEqual(B.InverseCDF(0.05d), true_icdf05, 0.0001d);
            Assert.AreEqual(B.InverseCDF(0.95d), true_icdf95, 0.0001d);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. See if beta is being created.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var b = new BetaDistribution(0, 0);
            Assert.AreEqual(0, b.Alpha);
            Assert.AreEqual(0, b.Beta);

            var b2 = new BetaDistribution(0, 1);
            Assert.AreEqual(0, b2.Alpha);
            Assert.AreEqual(1, b2.Beta);

            var b3 = new BetaDistribution(1, 0);
            Assert.AreEqual(1, b3.Alpha);
            Assert.AreEqual(0, b3.Beta);

            var b4 = new BetaDistribution(1, 1);
            Assert.AreEqual(1, b4.Alpha);
            Assert.AreEqual(1, b4.Beta);

            var b5 = new BetaDistribution(9, 1);
            Assert.AreEqual(9, b5.Alpha);
            Assert.AreEqual(1, b5.Beta);

            var b6 = new BetaDistribution(5, 100);
            Assert.AreEqual(5, b6.Alpha);
            Assert.AreEqual(100, b6.Beta);

            var b7 = new BetaDistribution(1,double.PositiveInfinity);
            Assert.AreEqual(1, b7.Alpha);
            Assert.AreEqual(double.PositiveInfinity, b7.Beta);

            var b8 = new BetaDistribution(double.PositiveInfinity,1);
            Assert.AreEqual(double.PositiveInfinity, b8.Alpha);
            Assert.AreEqual(1, b8.Beta);

            var b9 = new BetaDistribution(0,double.PositiveInfinity);
            Assert.AreEqual(0, b9.Alpha);
            Assert.AreEqual(double.PositiveInfinity, b9.Beta);

            var b10 = new BetaDistribution(double.PositiveInfinity, 0);
            Assert.AreEqual(double.PositiveInfinity, b10.Alpha);
            Assert.AreEqual(0,b10.Beta);
        }

        /// <summary>
        /// Checking failure with valid parameter variable.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var b = new BetaDistribution(double.NaN, 0);
            Assert.IsFalse(b.ParametersValid);

            var b2 = new BetaDistribution(-1, 1);
            Assert.IsFalse(b2.ParametersValid);

            var b3 = new BetaDistribution(double.PositiveInfinity, 0);
            Assert.IsFalse(b3.ParametersValid);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking string output.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var b = new BetaDistribution(1d, 1d);
            Assert.AreEqual("Shape (α)", b.ParametersToString[0, 0]);
            Assert.AreEqual("Shape (β)", b.ParametersToString[1, 0]);
            Assert.AreEqual("1", b.ParametersToString[0, 1]);
            Assert.AreEqual("1",b.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new BetaDistribution(9, 1);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        ///  Verified using MathNet-Numerics testing. Checking mean of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var b = new BetaDistribution(1,2);
            Assert.AreEqual((double)1/3,b.Mean);

            var b2 = new BetaDistribution(9,1);
            Assert.AreEqual(0.9, b2.Mean);

            var b3 = new BetaDistribution(1,double.PositiveInfinity);
            Assert.AreEqual(0, b3.Mean);

            var b4 = new BetaDistribution(0, 0);
            Assert.AreEqual(double.NaN, b4.Mean);

            var b5 = new BetaDistribution(double.PositiveInfinity, 1);
            Assert.AreEqual(double.NaN, b5.Mean);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking median of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var b = new BetaDistribution(2, 2);
            Assert.AreEqual(0.5, b.Median);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking mode of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var b = new BetaDistribution(1, 1);
            Assert.AreEqual(0.5, b.Mode);

            var b2 = new BetaDistribution(0, 2);
            Assert.AreEqual(0, b2.Mode);

            var b3 = new BetaDistribution(2, 0);
            Assert.AreEqual(1, b3.Mode);

            var b4 = new BetaDistribution(9, 1);
            Assert.AreEqual(1, b4.Mode);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking standard deviation with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var b = new BetaDistribution(2, 3);
            Assert.AreEqual(0.2, b.StandardDeviation);

            var b2 = new BetaDistribution(5, 1);
            Assert.AreEqual(0.14085,b2.StandardDeviation,1e-04);

            var b3 = new BetaDistribution(7, 3);
            Assert.AreEqual(0.138169, b3.StandardDeviation, 1e-04);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking skewness of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var b = new BetaDistribution(0, 0);
            Assert.AreEqual(double.NaN, b.Skewness);

            var b4 = new BetaDistribution(1, 1);
            Assert.AreEqual(0, b4.Skewness);

            var b5 = new BetaDistribution(9, 1);
            Assert.AreEqual(-1.474055, b5.Skewness,1e-04);
        }

        /// <summary>
        /// Checking Kurtosis of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var b = new BetaDistribution(2, 2);
            Assert.AreEqual(2.14285,b.Kurtosis,1e-04);

            var b2 = new BetaDistribution(5,2);
            Assert.AreEqual(2.88, b2.Kurtosis);

            var b3 = new BetaDistribution(2, 5);
            Assert.AreEqual(2.88,b3.Kurtosis);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var b = new BetaDistribution(1, 2);
            Assert.AreEqual(0,b.Minimum);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var b = new BetaDistribution(1, 2);
            Assert.AreEqual(1, b.Maximum);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing PDF function.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var b = new BetaDistribution(1, 1);
            Assert.AreEqual(1, b.PDF(0));
            Assert.AreEqual(1, b.PDF(0.5));
            Assert.AreEqual(1, b.PDF(1));

            var b2 = new BetaDistribution(9, 1);
            Assert.AreEqual(0, b2.PDF(0));
            Assert.AreEqual(0.035156, b2.PDF(0.5),1e-04);
            Assert.AreEqual(8.9999, b2.PDF(1),1e-04);
            Assert.AreEqual(0, b2.PDF(-1));
            Assert.AreEqual(0, b2.PDF(2));

            var b3 = new BetaDistribution(5, 100);
            Assert.AreEqual(0, b3.PDF(0));
            Assert.AreEqual(0, b3.PDF(1));
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var b = new BetaDistribution(1, 1);
            Assert.AreEqual(0, b.CDF(0));
            Assert.AreEqual(0.5, b.CDF(0.5));
            Assert.AreEqual(1, b.CDF(1));

            var b2 = new BetaDistribution(9, 1);
            Assert.AreEqual(0,b2.CDF(0));
            Assert.AreEqual(0.001953125, b2.CDF(0.5));
            Assert.AreEqual(1, b2.CDF(1));

            var b3 = new BetaDistribution(5, 100);
            Assert.AreEqual(0, b3.CDF(0));
            Assert.AreEqual(1, b3.CDF(1));
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Testing InverseCDF function.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var b = new BetaDistribution(1,1);
            Assert.AreEqual(1, b.InverseCDF(1));

            var b2 = new BetaDistribution(9, 1);
            Assert.AreEqual(0,b2.InverseCDF(0));
            Assert.AreEqual(1, b2.InverseCDF(1));

            var b3 = new BetaDistribution(5, 100);
            Assert.AreEqual(0,b3.InverseCDF(0));
        }
    }
}
