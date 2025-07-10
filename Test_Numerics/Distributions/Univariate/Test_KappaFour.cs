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
using Numerics.Mathematics;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Kappa-4 distribution algorithm.
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
    public class Test_KappaFour
    {

        [TestMethod]
        public void Test_K4_LMOM()
        {
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var kappa4 = new KappaFour();
            kappa4.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
            double xi = kappa4.Xi;
            double a = kappa4.Alpha;
            double k = kappa4.Kappa;
            double h= kappa4.Hondo;
            double true_xi = 8.68360234;
            double true_a = 3.10384972;
            double true_k = 0.14470737;
            double true_h = -0.07348014;

            Assert.AreEqual(xi, true_xi, 0.0001d);
            Assert.AreEqual(a, true_a, 0.0001d);
            Assert.AreEqual(k, true_k, 0.0001d);
            Assert.AreEqual(h, true_h, 0.0001d);

            var lmom = kappa4.LinearMomentsFromParameters(kappa4.GetParameters);
            Assert.AreEqual(lmom[0], 9.95751634d, 0.0001d);
            Assert.AreEqual(lmom[1], 1.98224114d, 0.0001d);
            Assert.AreEqual(lmom[2], 0.06380885d, 0.0001d);
            Assert.AreEqual(lmom[3], 0.12442297d, 0.0001d);
        }

        /// <summary>
        /// Verification using R lmom package. 
        /// </summary>
        [TestMethod]
        public void Test_K4_CDF()
        {
            var x = new double[] { 5, 10, 12, 15, 18 };
            var p = new double[5];
            var true_p = new double[] { 0.07168831, 0.53317660, 0.73279234, 0.91293987, 0.97980084 };
            var k4 = new KappaFour(8.7, 3.1, 0.14, -0.1);
            for (int i = 0; i < 5; i++)
            {
                p[i] = k4.CDF(x[i]);
                Assert.AreEqual(true_p[i], p[i], 1E-7);
            }
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(x[i], k4.InverseCDF(p[i]), 1E-7);
            }

        }


        [TestMethod]
        public void Test_K4_Dist()
        {
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var kappa4 = new KappaFour();
            kappa4.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);

            double pdf = kappa4.PDF(14.9);
            double cdf = kappa4.CDF(14.9);
            double invcdf = kappa4.InverseCDF(cdf);
            double true_pdf = 0.0385446;
            double true_cdf = 0.9106253;
            double true_invcdf = 14.9;

            Assert.AreEqual(pdf, true_pdf, 0.0001d);
            Assert.AreEqual(cdf, true_cdf, 0.0001d);
            Assert.AreEqual(invcdf, true_invcdf, 0.0001d);

        }

        [TestMethod]
        public void Test_K4_PartialDerivatives()
        {
           
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var kappa4 = new KappaFour();
            kappa4.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
            double p = 0.999;

            var pd1 = kappa4.QuantileGradient(p);
            var pd2 = NumericalDerivative.Gradient(x => 
            {
                var K4 = new KappaFour();
                K4.SetParameters(x);
                return K4.InverseCDF(p);
            },kappa4.GetParameters);

            for (int i = 0; i < pd1.Length; i++)
            {
                Assert.AreEqual(pd1[i], pd2[i], 1E-2);
            }
        }

        /// <summary>
        /// Checking that Kappa-4 is being created with inputs.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var k4 = new KappaFour();
            Assert.AreEqual(k4.Xi, 100);
            Assert.AreEqual(k4.Alpha, 10);
            Assert.AreEqual(k4.Kappa, 0);
            Assert.AreEqual(k4.Hondo, 0);

            var k4ii = new KappaFour(100, 10, 1, 1);
            Assert.AreEqual(k4ii.Xi, 100);
            Assert.AreEqual(k4ii.Alpha, 10);
            Assert.AreEqual(k4ii.Kappa, 1);
            Assert.AreEqual(k4ii.Hondo, 1);
        }

        /// <summary>
        /// Testting Kappa-4 with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var k4 = new KappaFour(double.NaN,double.NaN,double.NaN, double.NaN);
            Assert.IsFalse(k4.ParametersValid);

            var k4ii = new KappaFour(double.PositiveInfinity,double.PositiveInfinity,double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(k4ii.ParametersValid);

            var k4iii = new KappaFour(100, 0, 1, 1);
            Assert.IsFalse(k4iii.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var k4 = new KappaFour();
            Assert.AreEqual(k4.ParametersToString[0, 0], "Location (ξ)");
            Assert.AreEqual(k4.ParametersToString[1, 0], "Scale (α)");
            Assert.AreEqual(k4.ParametersToString[2, 0], "Shape (κ)");
            Assert.AreEqual(k4.ParametersToString[3, 0], "Shape (h)");
            Assert.AreEqual(k4.ParametersToString[0, 1], "100");
            Assert.AreEqual(k4.ParametersToString[1, 1], "10");
            Assert.AreEqual(k4.ParametersToString[2, 1], "0");
            Assert.AreEqual(k4.ParametersToString[3, 1], "0");
        }
     }
}
