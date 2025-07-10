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
    /// Testing the Geometric distribution algorithm.
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
    public class Test_Geometric
    {
        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_GeometricDist()
        {
            double true_mean = 1.3809523809523812d;
            double true_mode = 0.0d;
            double true_median = 1.0d;
            double true_stdDev = Math.Sqrt(3.2879818594104315d);
            double true_skew = 2.0746d;
            double true_kurt = 9.3041d;
            double true_pdf = 0.141288d;
            double true_cdf = 0.80488799999999994d;
            double true_icdf05 = 0.0d;
            double true_icdf95 = 5.0d;
            var G = new Geometric(0.42d);
            Assert.AreEqual(G.Mean, true_mean, 0.0001d);
            Assert.AreEqual(G.Median, true_median, 0.0001d);
            Assert.AreEqual(G.Mode, true_mode, 0.0001d);
            Assert.AreEqual(G.StandardDeviation, true_stdDev, 0.0001d);
            Assert.AreEqual(G.Skewness, true_skew, 0.0001d);
            Assert.AreEqual(G.Kurtosis, true_kurt, 0.0001d);
            Assert.AreEqual(G.PDF(2.0d), true_pdf, 0.0001d);
            Assert.AreEqual(G.CDF(2.0d), true_cdf, 0.0001d);
            Assert.AreEqual(G.InverseCDF(0.05d), true_icdf05, 0.0001d);
            Assert.AreEqual(G.InverseCDF(0.95d), true_icdf95, 0.0001d);
        }

        /// <summary>
        /// See if parameters can create geometric distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var G = new Geometric();
            Assert.AreEqual(G.ProbabilityOfSuccess, 0.5);

            var G2 = new Geometric(0);
            Assert.AreEqual(G2.ProbabilityOfSuccess, 0);

            var G3 = new Geometric(1);
            Assert.AreEqual(G3.ProbabilityOfSuccess, 1);
        }

        /// <summary>
        /// Testing geometric distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var G = new Geometric(-1);
            Assert.IsFalse(G.ParametersValid);

            var G2 = new Geometric(2);
            Assert.IsFalse(G2.ParametersValid);

            var G3 = new Geometric(double.NaN);
            Assert.IsFalse(G3.ParametersValid);

            var G4 = new Geometric(double.PositiveInfinity);
            Assert.IsFalse(G4.ParametersValid);
        }

        /// <summary>
        /// Checking parameters to string function.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var G = new Geometric();
            Assert.AreEqual(G.ParametersToString[0, 0], "Probability (p)");
            Assert.AreEqual(G.ParametersToString[0, 1], "0.5");
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var G = new Geometric();
            Assert.AreEqual(G.Mean, 1);

            var G2 = new Geometric(0.3);
            Assert.AreEqual(G2.Mean, 2.3333, 1e-04);
        }

        /// <summary>
        /// Testing median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var G = new Geometric(0.0001);
            Assert.AreEqual(G.Median, 6931);

            var G2 = new Geometric(0.1);
            Assert.AreEqual(G2.Median, 6);

            var G3 = new Geometric(0.9);
            Assert.AreEqual(G3.Median, 0);
        }

        /// <summary>
        /// Testing mode function outputs 0.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var G = new Geometric();
            Assert.AreEqual(G.Mode,0);

            var G2 = new Geometric(0);
            Assert.AreEqual(G2.Mode,0);
        }

        /// <summary>
        /// Testing Standard deviation with different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var G = new Geometric();
            Assert.AreEqual(G.StandardDeviation, 1.41421,1e-04);

            var G2 = new Geometric(0.3);
            Assert.AreEqual(G2.StandardDeviation, 2.78886, 1e-04);
        }

        /// <summary>
        /// Testing skewness function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var G = new Geometric();
            Assert.AreEqual(G.Skewness, 2.12132, 1e-04);

            var G2 = new Geometric(0.3);
            Assert.AreEqual(G2.Skewness, 2.03188, 1e-04);
        }

        /// <summary>
        /// Testing Kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var G = new Geometric();
            Assert.AreEqual(G.Kurtosis, 9.5);

            var G2 = new Geometric(0.3);
            Assert.AreEqual(G2.Kurtosis, 9.12857, 1e-04);
        }

        /// <summary>
        /// Testing minimum is 0 and maximum is positive infinity.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var G = new Geometric();
            Assert.AreEqual(G.Minimum, 0);
            Assert.AreEqual(G.Maximum,double.PositiveInfinity);

            var G2 = new Geometric(0.3);
            Assert.AreEqual(G2.Minimum, 0);
            Assert.AreEqual(G2.Maximum, double.PositiveInfinity);
        }

        /// <summary>
        /// Testing PDF method with different parameters and locations.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var G = new Geometric();
            Assert.AreEqual(G.PDF(0), 0.5);
            Assert.AreEqual(G.PDF(2), 0.125);
            Assert.AreEqual(G.PDF(-1), 0);

            var G2 = new Geometric(0.3);
            Assert.AreEqual(G2.PDF(0), 0.3);
            Assert.AreEqual(G2.PDF(2.5), 0.122989, 1e-05);
        }

        /// <summary>
        /// Testing CDF method with different parameters and locations.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var G = new Geometric();
            Assert.AreEqual(G.CDF(0), 0.5);
            Assert.AreEqual(G.CDF(2), 0.875);
            Assert.AreEqual(G.CDF(-1), 0);
            Assert.AreEqual(G.CDF(double.PositiveInfinity), 1);

            var G2 = new Geometric(0.3);
            Assert.AreEqual(G2.CDF(2), 0.657);
            Assert.AreEqual(G2.CDF(100), 1,1e-04);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var G = new Geometric();
            Assert.AreEqual(G.InverseCDF(0.3), 0);
            Assert.AreEqual(G.InverseCDF(0.7), 1, 1e-04);

            var G2 = new Geometric(0.3);
            Assert.AreEqual(G2.InverseCDF(0.5), 1, 1e-04);
            Assert.AreEqual(G2.InverseCDF(0.9), 6);
        }
    }
}
