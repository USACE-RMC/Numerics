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
    /// Testing the Cauchy distribution algorithm.
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
    public class Test_Cauchy
    {
        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_CauchyDist()
        {
            double true_mean = double.NaN;
            double true_mode = 0.42d;
            double true_median = 0.42d;
            double true_stdDev = double.NaN;
            double true_skew = double.NaN;
            double true_kurt = double.NaN;
            double true_pdf = 0.2009112009763413d;
            double true_cdf = 0.46968025841608563d;
            double true_icdf = 1.5130304686978195d;
            var C = new Cauchy(0.42d, 1.57d);
            Assert.AreEqual(C.Mean, true_mean);
            Assert.AreEqual(C.Median, true_median, 0.0001d);
            Assert.AreEqual(C.Mode, true_mode, 0.0001d);
            Assert.AreEqual(C.StandardDeviation, true_stdDev);
            Assert.AreEqual(C.Skewness, true_skew);
            Assert.AreEqual(C.Kurtosis, true_kurt);
            Assert.AreEqual(C.PDF(0.27d), true_pdf, 0.0001d);
            Assert.AreEqual(C.CDF(0.27d), true_cdf, 0.0001d);
            Assert.AreEqual(C.InverseCDF(0.69358638272337991d), true_icdf, 0.0001d);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. See if Cauchy is being created.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var c = new Cauchy(0, 0.1);
            Assert.AreEqual(0, c.X0);
            Assert.AreEqual(0.1, c.Gamma);

            var c2 = new Cauchy(0, 1);
            Assert.AreEqual(0, c2.X0);
            Assert.AreEqual(1, c2.Gamma);

            var c3 = new Cauchy(0, 10);
            Assert.AreEqual(0, c3.X0);
            Assert.AreEqual(10, c3.Gamma);

            var c4 = new Cauchy(10, 11);
            Assert.AreEqual(10, c4.X0);
            Assert.AreEqual(11, c4.Gamma);

            var c5 = new Cauchy(-5, 100);
            Assert.AreEqual(-5, c5.X0);
            Assert.AreEqual(100, c5.Gamma);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. See what parameters fail.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var c = new Cauchy(double.NaN, 1);
            Assert.IsFalse(c.ParametersValid);

            var c2 = new Cauchy(1,double.NaN);
            Assert.IsFalse(c2.ParametersValid);
            
            var c3 = new Cauchy(double.NaN, double.NaN);
            Assert.IsFalse(c3.ParametersValid);

            var c4 = new Cauchy(1, 0);
            Assert.IsFalse(c4.ParametersValid);
        }

        /// <summary>
        /// Checking parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var c = new Cauchy(0, 1);
            Assert.AreEqual("Location (X0)", c.ParametersToString[0,0]);
            Assert.AreEqual("Scale (γ)", c.ParametersToString[1,0]);
            Assert.AreEqual("0", c.ParametersToString[0,1]);
            Assert.AreEqual("1",c.ParametersToString[1,1]);
        }

        /// <summary>
        /// Checking mean is NaN.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var c = new Cauchy(0, 1);
            Assert.AreEqual(double.NaN, c.Mean);
        }

        /// <summary>
        /// Checking median matches location parameter.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var c = new Cauchy(0, 1);
            Assert.AreEqual(0, c.Median);
        }

        /// <summary>
        /// Checking mode matches location parameter. 
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var c = new Cauchy(0, 1);
            Assert.AreEqual(0, c.Mode);
        }

        /// <summary>
        /// Checking that standard deviation is NaN.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var c = new Cauchy(0, 1);
            Assert.AreEqual(double.NaN,c.StandardDeviation);
        }

        /// <summary>
        /// Checking that skew is NaN.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var c = new Cauchy(0, 1);
            Assert.AreEqual(double.NaN, c.Skewness);
        }

        /// <summary>
        /// Checking that Kurtosis is NaN.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var c = new Cauchy(0,1);
            Assert.AreEqual(double.NaN , c.Kurtosis);
        }

        /// <summary>
        /// Checking minimum.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var c = new Cauchy(0, 1);
            Assert.AreEqual(double.NegativeInfinity, c.Minimum);
        }

        /// <summary>
        /// Checking maximum.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var c = new Cauchy(0, 1);
            Assert.AreEqual(double.PositiveInfinity, c.Maximum);
        }

        /// <summary>
        /// Checking PDF function with different inputs.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var c = new Cauchy(0, 0.1);
            Assert.AreEqual(0.0012727, c.PDF(-5), 1e-04);
            Assert.AreEqual(3.18309, c.PDF(0), 1e-04);
            Assert.AreEqual(0.031515, c.PDF(1), 1e-04);

            var c2 = new Cauchy(-5, 100);
            Assert.AreEqual(0.00318309, c2.PDF(-5), 1e-04);
            Assert.AreEqual(0.00317516, c2.PDF(0), 1e-04);
            Assert.AreEqual(0.00317168, c2.PDF(1), 1e-04);
        }
        
        /// <summary>
        /// Checking CDF function with different inputs.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var c = new Cauchy(0, 0.1);
            Assert.AreEqual(0.0063653, c.CDF(-5),1e-04);
            Assert.AreEqual(0.5, c.CDF(0), 1e-04);
            Assert.AreEqual(0.968274, c.CDF(1), 1e-04);

            var c2 = new Cauchy(-5, 100);
            Assert.AreEqual(0.5, c2.CDF(-5), 1e-04);
            Assert.AreEqual(0.51590, c2.CDF(0), 1e-04);
            Assert.AreEqual(0.51907, c2.CDF(1), 1e-04);
        }

        /// <summary>
        /// Checking inverse CDF function with different probabilities and inputs.
        /// InverseCDF will be undefined if probability is either 0 or 1.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var c = new Cauchy(0, 0.1);
            Assert.AreEqual(-0.307768, c.InverseCDF(0.1), 1e-04);
            Assert.AreEqual(-0.072654, c.InverseCDF(0.3), 1e-04);
            Assert.AreEqual(0.307768, c.InverseCDF(0.9), 1e-04);

            var c2 = new Cauchy(-5, 100);
            Assert.AreEqual(-312.76835, c2.InverseCDF(0.1), 1e-04);
            Assert.AreEqual(-77.65425, c2.InverseCDF(0.3), 1e-04);
            Assert.AreEqual(302.768353, c2.InverseCDF(0.9), 1e-04);
        }
    }
}
