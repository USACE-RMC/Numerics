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

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Pert distribution algorithm.
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
    public class Test_Pert
    {
        /// <summary>
        /// Verified against Palisade's @Risk.
        /// </summary>
        [TestMethod]
        public void Test_PertDist()
        {
            var P = new Pert(1d, 2d, 3d);
            var GB = GeneralizedBeta.PERT(1d, 2d, 3d);

            double true_mean = 2d;
            double true_median = 2d;
            double true_mode = 2d;
            double true_pdf = GB.PDF(1.27);
            double true_cdf = GB.CDF(1.27);
            double true_icdf = 1.27d;

            Assert.AreEqual(P.Mean, true_mean, 0.0001d);
            Assert.AreEqual(P.Median, true_median, 0.0001d);
            Assert.AreEqual(P.Mode, true_mode, 0.0001d);
            Assert.AreEqual(P.PDF(1.27d), true_pdf, 0.0001d);
            Assert.AreEqual(P.CDF(1.27d), true_cdf, 0.0001d);
            Assert.AreEqual(P.InverseCDF(P.CDF(1.27d)), true_icdf, 0.0001d);

        }

        /// <summary>
        /// Verified against R "mc2d" package
        /// </summary>
        [TestMethod]
        public void Test_PertDist2()
        {
            var P = new Pert(-1, 0, 1);
            var d = P.PDF(-0.5);
            var p = P.CDF(-0.5);
            var q = P.InverseCDF(p);

            Assert.AreEqual(0.5273438, d,  1E-5);
            Assert.AreEqual(0.1035156, p,  1E-5);
            Assert.AreEqual(-0.5, q,  1E-5);

        }

        /// <summary>
        /// Verified using R 'mc2d'
        /// </summary>
        [TestMethod]
        public void Test_Pert_R()
        {
            var x = new double[] { 0.1, 0.25, 0.5, 0.75, 0.9 };
            var p = new double[5];
            var true_p = new double[] { 0.0814600, 0.3671875, 0.8125000, 0.9843750, 0.9995400 };
            var pert = new Pert(0, 0.25, 1);
            for (int i = 0; i < 5; i++)
            {
                p[i] = pert.CDF(x[i]);
                Assert.AreEqual(true_p[i], p[i], 1E-7);
            }
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(x[i], pert.InverseCDF(p[i]), 1E-7);
            }
        }

        /// <summary>
        /// Testing fit with moments.
        /// </summary>
        [TestMethod]
        public void Test_Pert_MOM()
        {
            var dist = new Pert(-2, 10, 35);
            var sample = dist.GenerateRandomValues(1000, 12345);
            dist.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            Assert.AreEqual(-2, dist.Min, 1);
            Assert.AreEqual(10, dist.MostLikely, 1);
            Assert.AreEqual(35, dist.Max, 5);
        }

        /// <summary>
        /// Testing fit with likelihood.
        /// </summary>
        [TestMethod]
        public void Test_Pert_MLE()
        {
            var dist = new Pert(-2, 10, 35);
            var sample = dist.GenerateRandomValues(1000, 12345);
            dist.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(-2, dist.Min, 1);
            Assert.AreEqual(10, dist.MostLikely, 1);
            Assert.AreEqual(35, dist.Max, 1);
        }

        /// <summary>
        /// Verifying that input parameters create distribution.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var p = new Pert();
            Assert.AreEqual(0,p.Min);
            Assert.AreEqual(1, p.Max);
            Assert.AreEqual(0.5, p.Mode);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            var p = new Pert(double.NaN, double.NaN,double.NaN);
            Assert.IsFalse(p.ParametersValid);
            
            var p2 = new Pert(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
            Assert.IsFalse(p2.ParametersValid);

            var p3 = new Pert(3, 2, 1);
            Assert.IsFalse(p3.ParametersValid);

            var p4 = new Pert(1, 2, 0.5);
            Assert.IsFalse(p4.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString()
        /// </summary>
        [TestMethod]
        public void Test_ParametersToString()
        {
            var p = new Pert();
            Assert.AreEqual("Min (a)",p.ParametersToString[0, 0] );
            Assert.AreEqual("Most Likely (c)", p.ParametersToString[1, 0]);
            Assert.AreEqual("Max (b)", p.ParametersToString[2, 0]);
            Assert.AreEqual("0", p.ParametersToString[0, 1]);
            Assert.AreEqual("0.5", p.ParametersToString[1, 1]);
            Assert.AreEqual("1", p.ParametersToString[2, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Pert(-2, 10, 35);
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
            var p = new Pert();
            Assert.AreEqual(0.5, p.Mean);

            var p2 = new Pert(0, 0, 0);
            Assert.AreEqual(0, p2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod]
        public void Test_Median()
        {
            var p = new Pert();
            Assert.AreEqual(0.5, p.Median);

            var p2 = new Pert(0,0,0);
            Assert.AreEqual(0, p2.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var p = new Pert();
            Assert.AreEqual(0.5, p.Mode);        
        }

        /// <summary>
        /// Testing Standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var p = new Pert();
            Assert.AreEqual(0.1889, p.StandardDeviation, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var p = new Pert();
            Assert.AreEqual(0,p.Skewness);
        }

        /// <summary>
        /// Testing kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var p = new Pert();
            Assert.AreEqual(2.3333, p.Kurtosis, 1e-04);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var p = new Pert();
            Assert.AreEqual(0,p.Minimum);
            Assert.AreEqual(1, p.Maximum);

            var p2 = new Pert(1, 1.5, 2);
            Assert.AreEqual(1, p2.Minimum);
            Assert.AreEqual(2, p2.Maximum);
        }
    }
}
