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
    /// Testing the Pareto distribution algorithm.
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
    public class Test_Pareto
    {

        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_ParetoDist()
        {
            double true_mean = 0.63d;
            double true_mode = 0.42d;
            double true_median = 0.52916684095584676d;
            double true_stdDev = Math.Sqrt(0.13229999999999997d);
            double true_pdf = 0.057857142857142857d;
            double true_cdf = 0.973d;
            double true_icdf = 1.4d;
            double true_icdf05 = 0.4272d;
            double true_icdf95 = 1.1401d;
            var PA = new Pareto(0.42d, 3d);
            Assert.AreEqual(PA.Mean, true_mean, 0.0001d);
            Assert.AreEqual(PA.Median, true_median, 0.0001d);
            Assert.AreEqual(PA.Mode, true_mode, 0.0001d);
            Assert.AreEqual(PA.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(PA.PDF(1.4d), true_pdf, 0.0001d);
            Assert.AreEqual(PA.CDF(1.4d), true_cdf, 0.0001d);
            Assert.AreEqual(PA.InverseCDF(PA.CDF(1.4d)), true_icdf, 0.0001d);
            Assert.AreEqual(PA.InverseCDF(0.05d), true_icdf05, 0.0001d);
            Assert.AreEqual(PA.InverseCDF(0.95d), true_icdf95, 0.0001d);
            PA.SetParameters(new[] { 1d, 10d });
            double true_skew = 2.8111d;
            double true_kurt = 17.8286d;
            Assert.AreEqual(PA.Skewness, true_skew, 0.0001d);
            Assert.AreEqual(PA.Kurtosis, true_kurt, 0.0001d);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var p = new Pareto();
            Assert.AreEqual(p.Xm, 1);
            Assert.AreEqual(p.Alpha, 10);

            var p2 = new Pareto(1, 1);
            Assert.AreEqual(p2.Xm, 1);
            Assert.AreEqual(p2.Alpha, 1);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var p = new Pareto(double.NaN,double.PositiveInfinity);
            Assert.IsFalse(p.ParametersValid);

            var p2 = new Pareto(double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(p2.ParametersValid);

            var p3 = new Pareto(0, 0);
            Assert.IsFalse(p3.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var p = new Pareto();
            Assert.AreEqual(p.ParametersToString[0, 0], "Scale (Xm)");
            Assert.AreEqual(p.ParametersToString[1, 0], "Shape (α)");
            Assert.AreEqual(p.ParametersToString[0, 1], "1");
            Assert.AreEqual(p.ParametersToString[1, 1], "10");
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Pareto();
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
            var p = new Pareto();
            Assert.AreEqual(p.Mean, 1.1111, 1e-04);

            var p2 = new Pareto(1,1);
            Assert.AreEqual(p2.Mean, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var p = new Pareto();
            Assert.AreEqual(p.Median, 1.07177, 1e-04);

            var p2 = new Pareto(1, 1);
            Assert.AreEqual(p2.Median, 2);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var p = new Pareto();
            Assert.AreEqual(p.Mode, 1);

            var p2 = new Pareto(2, 1);
            Assert.AreEqual(p2.Mode, 2);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var p = new Pareto();
            Assert.AreEqual(p.StandardDeviation, 0.12422, 1e-04);

            var p2 = new Pareto(1, 2);
            Assert.AreEqual(p2.StandardDeviation, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var p = new Pareto();
            Assert.AreEqual(p.Skewness, 2.81105, 1e-04);

            var p2 = new Pareto(1, 3);
            Assert.AreEqual(p2.Skewness, double.NaN);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var p = new Pareto();
            Assert.AreEqual(p.Kurtosis, 17.82857, 1e-04);

            var p2 = new Pareto(1, 4);
            Assert.AreEqual(p2.Kurtosis,double.NaN);
        }

        /// <summary>
        /// Testing minimum and maximum functions
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var p = new Pareto();
            Assert.AreEqual(p.Minimum, 1);
            Assert.AreEqual(p.Maximum,double.PositiveInfinity);

            var p2 = new Pareto(2,3);
            Assert.AreEqual(p2.Minimum, 2);
            Assert.AreEqual(p2.Maximum,double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var p = new Pareto(1,1);
            Assert.AreEqual(p.PDF(1), 1);
            Assert.AreEqual(p.PDF(1.5), 4d / 9d);

            var p2 = new Pareto(3, 2);
            Assert.AreEqual(p2.PDF(3), 2d / 3d);
            Assert.AreEqual(p2.PDF(5), 18d / 125d);
        }

        /// <summary>
        /// Testing CDF.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var p = new Pareto();
            Assert.AreEqual(p.CDF(0), 0);
            Assert.AreEqual(p.CDF(2), 0.9990, 1e-04);
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var p = new Pareto();
            Assert.AreEqual(p.InverseCDF(0), 1);
            Assert.AreEqual(p.InverseCDF(0.3), 1.0363, 1e-04);
        }
    }
}
