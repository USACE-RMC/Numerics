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
    /// Testing the Triangular distribution algorithm.
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
    public class Test_Triangular
    {

        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_TriangularDist()
        {
            double true_mean = 3.3333333333333335d;
            double true_median = 3.2613872124741694d;
            double true_mode = 3.0d;
            double true_stdDev = Math.Sqrt(1.0555555555555556d);
            double true_pdf = 0.20000000000000001d;
            double true_cdf = 0.10000000000000001d;
            double true_icdf = 2.0d;
            var T = new Triangular(1, 3, 6);

            Assert.AreEqual(T.Mean, true_mean, 0.0001d);
            Assert.AreEqual(T.Median, true_median, 0.0001d);
            Assert.AreEqual(T.Mode, true_mode, 0.0001d);
            Assert.AreEqual(T.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(T.PDF(2.0d), true_pdf, 0.0001d);
            Assert.AreEqual(T.CDF(2.0d), true_cdf, 0.0001d);
            Assert.AreEqual(T.InverseCDF(true_cdf), true_icdf, 0.0001d);
        }

        /// <summary>
        /// Verified using R 'mc2d'
        /// </summary>
        [TestMethod]
        public void Test_Triangular_R()
        {
            var x = new double[] { 0.1, 0.25, 0.5, 0.75, 0.9 };
            var p = new double[5];
            var true_p = new double[] { 0.0400000, 0.2500000, 0.6666667, 0.9166667, 0.9866667 };
            var tri = new Triangular(0, 0.25, 1);
            for (int i = 0; i < 5; i++)
            {
                p[i] = tri.CDF(x[i]);
                Assert.AreEqual(true_p[i], p[i], 1E-7);
            }
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(x[i], tri.InverseCDF(p[i]), 1E-7);
            }
        }

        /// <summary>
        /// Test estimation with moments.
        /// </summary>
        [TestMethod]
        public void Test_Triangular_MOM()
        {
            var dist = new Triangular(-2, 10, 35);
            var sample = dist.GenerateRandomValues(1000, 12345);
            dist.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            Assert.AreEqual(-2, dist.Min, 1);
            Assert.AreEqual(10, dist.MostLikely, 1);
            Assert.AreEqual(35, dist.Max, 1);
        }

        /// <summary>
        /// Test estimation with likelihood.
        /// </summary>
        [TestMethod]
        public void Test_Triangular_MLE()
        {
            var dist = new Triangular(-2, 10, 35);
            var sample = dist.GenerateRandomValues(1000, 12345);
            dist.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(-2, dist.Min, 1);
            Assert.AreEqual(10, dist.MostLikely, 1);
            Assert.AreEqual(35, dist.Max, 1);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var T = new Triangular();
            Assert.AreEqual(T.Min, 0);
            Assert.AreEqual(T.Mode, 0.5);
            Assert.AreEqual(T.Max, 1);

            var T2 = new Triangular(-1,1,2);
            Assert.AreEqual(T2.Min, -1);
            Assert.AreEqual(T2.Mode, 1);
            Assert.AreEqual(T2.Max, 2);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            var T = new Triangular(double.NaN,double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(T.ParametersValid);

            var T2 = new Triangular(double.PositiveInfinity,double.NaN,double.NaN);
            Assert.IsFalse(T2.ParametersValid);

            var T3 = new Triangular(4, 1, 0);
            Assert.IsFalse(T3.ParametersValid);

            var T4 = new Triangular(1, 0, -1);
            Assert.IsFalse(T4.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod]
        public void Test_ParametersToString()
        {
            var T = new Triangular();
            Assert.AreEqual(T.ParametersToString[0, 0], "Min (a)");
            Assert.AreEqual(T.ParametersToString[1, 0], "Most Likely (c)");
            Assert.AreEqual(T.ParametersToString[2, 0], "Max (b)");
            Assert.AreEqual(T.ParametersToString[0, 1], "0");
            Assert.AreEqual(T.ParametersToString[1, 1], "0.5");
            Assert.AreEqual(T.ParametersToString[2, 1], "1");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Triangular(1, 3, 6);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            var T = new Triangular();
            Assert.AreEqual(T.Mean, 0.5);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(T2.Mean, 3.3333, 1e-04);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod]
        public void Test_Median()
        {
            var T = new Triangular();
            Assert.AreEqual(T.Median, 0.5);

            var T2 = new Triangular(1,3,6);
            Assert.AreEqual(T2.Median, 3.26138, 1e-05);
        }

        /// <summary>
        /// Testing mode
        /// </summary>
        [TestMethod]
        public void Test_Mode()
        {
            var T = new Triangular();
            Assert.AreEqual(T.Mode, 0.5);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(T2.Mode, 3);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_StandardDeviation()
        {
            var T = new Triangular();
            Assert.AreEqual(T.StandardDeviation, 0.20412, 1e-04);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(T2.StandardDeviation, 1.02739, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod]
        public void Test_Skewness()
        {
            var T = new Triangular();
            Assert.AreEqual(T.Skewness, 0);
        }

        /// <summary>
        /// Testing kurtosis.
        /// </summary>
        [TestMethod]
        public void Test_Kurtosis()
        {
            var T = new Triangular();
            Assert.AreEqual(T.Kurtosis, 12d / 5d);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(T2.Kurtosis, 12d / 5d);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod]
        public void Test_MinMax()
        {
            var T = new Triangular();
            Assert.AreEqual(T.Minimum, 0);
            Assert.AreEqual(T.Maximum, 1);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(T2.Minimum, 1);
            Assert.AreEqual(T2.Maximum, 6);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var T = new Triangular();
            Assert.AreEqual(T.PDF(-1), 0);
            Assert.AreEqual(T.PDF(0.4), 1.6);
            Assert.AreEqual(T.PDF(0.6), 1.6);
            Assert.AreEqual(T.PDF(0.5), 2);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(T2.PDF(2), 0.2, 1e-04);
        }

        /// <summary>
        /// Testing CDF.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var T = new Triangular();
            Assert.AreEqual(T.CDF(-1), 0);
            Assert.AreEqual(T.CDF(2), 1);
            Assert.AreEqual(T.CDF(0.4), 0.32,1e-04);
            Assert.AreEqual(T.CDF(0.6), 0.68,1e-04);

            var T2 = new Triangular(1,3, 6);
            Assert.AreEqual(T2.CDF(2), 0.1, 1e-04);
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF()
        {
            var T = new Triangular();
            Assert.AreEqual(T.InverseCDF(0), 0);
            Assert.AreEqual(T.InverseCDF(1), 1);
            Assert.AreEqual(T.InverseCDF(0.2), 0.31622, 1e-04);
            Assert.AreEqual(T.InverseCDF(0.5), 0.5);
        }
    }
}
