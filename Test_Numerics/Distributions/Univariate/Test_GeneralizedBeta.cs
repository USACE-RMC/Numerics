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
    /// Testing the Generalized Beta distribution algorithm.
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
    public class Test_GeneralizedBeta
    {
        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod]
        public void Test_GenBeta()
        {

            double alpha = 0.42;
            double beta = 1.57;
            var B = new GeneralizedBeta(alpha, beta);

            double true_mean = 0.21105527638190955d;
            double true_median = 0.11577706212908731d;
            double true_mode = 57.999999999999957d;
            double true_var = 0.055689279830523512d;
            double true_pdf = 0.94644031936694828d;
            double true_cdf = 0.69358638272337991d;
            double true_icdf = 0.27d;
            
            Assert.AreEqual(B.Mean, true_mean, 0.0001d);
            Assert.AreEqual(B.Median, true_median, 0.0001d);
            Assert.AreEqual(B.Mode, true_mode, 0.0001d);
            Assert.AreEqual(B.Variance, true_var, 0.0001d);
            Assert.AreEqual(B.PDF(0.27d), true_pdf, 0.0001d);
            Assert.AreEqual(B.CDF(0.27d), true_cdf, 0.0001d);
            Assert.AreEqual(B.InverseCDF(B.CDF(0.27d)), true_icdf, 0.0001d);


        }

        /// <summary>
        /// Verified using R 'mc2d'
        /// </summary>
        [TestMethod]
        public void Test_GenBeta_R()
        {
            var x = new double[] { 0.1, 0.25, 0.5, 0.75, 0.9 };
            var p = new double[5];
            var true_p = new double[] { 0.271000, 0.578125, 0.875000, 0.984375, 0.999000 };
            var B = new GeneralizedBeta(1, 3);
            for (int i = 0; i < 5; i++)
            {
                p[i] = B.CDF(x[i]);
                Assert.AreEqual(true_p[i], p[i], 1E-7);
            }
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(x[i], B.InverseCDF(p[i]), 1E-7);
            }
        }

        /// <summary>
        /// See if Generalized beta is being created.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var b = new GeneralizedBeta(0, 0,0,1);
            Assert.AreEqual(0, b.Alpha);
            Assert.AreEqual(0, b.Beta);
            Assert.AreEqual(0, b.Min);
            Assert.AreEqual(1, b.Max);

            var b2 = new GeneralizedBeta(0, 1,0,1);
            Assert.AreEqual(0, b2.Alpha);
            Assert.AreEqual(1, b2.Beta);
            Assert.AreEqual(0, b2.Min);
            Assert.AreEqual(1, b2.Max);

            var b3 = new GeneralizedBeta(1, 0,0,1);
            Assert.AreEqual(1, b3.Alpha);
            Assert.AreEqual(0, b3.Beta);
            Assert.AreEqual(0, b3.Min);
            Assert.AreEqual(1, b3.Max);

            var b4 = new GeneralizedBeta(1, 1,0,1);
            Assert.AreEqual(1, b4.Alpha);
            Assert.AreEqual(1, b4.Beta);
            Assert.AreEqual(0, b4.Min);
            Assert.AreEqual(1, b4.Max);

            var b5 = new GeneralizedBeta(9, 1, 0, 1);
            Assert.AreEqual(9, b5.Alpha);
            Assert.AreEqual(1, b5.Beta);
            Assert.AreEqual(0, b5.Min);
            Assert.AreEqual(1, b5.Max);
        }

        /// <summary>
        /// Check Generalized beta function with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var b = new GeneralizedBeta(double.NaN, 0);
            Assert.IsFalse(b.ParametersValid);

            var b2 = new GeneralizedBeta(-1, 1);
            Assert.IsFalse(b2.ParametersValid);

            var b3 = new GeneralizedBeta(double.PositiveInfinity, 0);
            Assert.IsFalse(b3.ParametersValid);

            var b4 = new GeneralizedBeta(1, 1, 1, 0);
            Assert.IsFalse(b4.ParametersValid);
        }

        /// <summary>
        /// Testing ParameterToString function.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var b = new GeneralizedBeta(1d, 1d,0,1);
            Assert.AreEqual("Shape (α)", b.ParametersToString[0, 0]);
            Assert.AreEqual("Shape (β)", b.ParametersToString[1, 0]);
            Assert.AreEqual("Min",b.ParametersToString[2, 0]);
            Assert.AreEqual("Max", b.ParametersToString[3,0]);

            Assert.AreEqual("1", b.ParametersToString[0, 1]);
            Assert.AreEqual("1", b.ParametersToString[1, 1]);
            Assert.AreEqual("0", b.ParametersToString[2, 1]);
            Assert.AreEqual("1", b.ParametersToString[3, 1]);    
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new GeneralizedBeta(2, 2, 0, 1);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(mom[0], dist.Mean, 1E-2);
            Assert.AreEqual(mom[1], dist.StandardDeviation, 1E-2);
            Assert.AreEqual(mom[2], dist.Skewness, 1E-2);
            Assert.AreEqual(mom[3], dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Validating the mean of this distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var b = new GeneralizedBeta(2, 2, 0, 1);
            Assert.AreEqual(b.Mean, 0.5);

            var b2 = new GeneralizedBeta(2, 2, -10, 10);
            Assert.AreEqual(b2.Mean, 0);
        }

        /// <summary>
        /// Verified using MathNet-Numerics testing. Checking median of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var b = new GeneralizedBeta(2, 2);
            Assert.AreEqual(0.5, b.Median);
        }

        /// <summary>
        /// Verifies the mode of the distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var b = new GeneralizedBeta();
            Assert.AreEqual(b.Mode, 0.5);

            var b2 = new GeneralizedBeta(2, 2, -10, 10);
            Assert.AreEqual(b2.Mode, 0);
        }

        /// <summary>
        /// Testing Standard Deviation with different Max and min values
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var b = new GeneralizedBeta();
            Assert.AreEqual(b.StandardDeviation, 0.223606, 1e-04);

            var b2 = new GeneralizedBeta(2, 2, -10, 10);
            Assert.AreEqual(b2.StandardDeviation, 4.47213,1e-04);
        }

        /// <summary>
        /// Verifying skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var b = new GeneralizedBeta();
            Assert.AreEqual(b.Skewness, 0);

            var b2 = new GeneralizedBeta(2, 10);
            Assert.AreEqual(b2.Skewness, 0.92140088,1e-04);
        }

        /// <summary>
        /// Checking Kurtosis of distribution with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var b = new GeneralizedBeta(2, 2);
            Assert.AreEqual(2.14285, b.Kurtosis, 1e-04);

            var b2 = new GeneralizedBeta(5, 2);
            Assert.AreEqual(2.88, b2.Kurtosis);

            var b3 = new GeneralizedBeta(2, 5);
            Assert.AreEqual(2.88, b3.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions of this distribution.
        /// </summary>
        [TestMethod()]
        public void Test_MinimumMaximum()
        {
            var b = new GeneralizedBeta();
            Assert.AreEqual(b.Minimum, 0);
            Assert.AreEqual(b.Maximum, 1);

            var b2 = new GeneralizedBeta(2, 2, -10, 10);
            Assert.AreEqual(b2.Minimum, -10);
            Assert.AreEqual(b2.Maximum, 10);
        }

        /// <summary>
        /// Verifying the PDF for Generalized Beta with known inputs from Test_Beta.cs.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var b = new GeneralizedBeta(1, 1,0,2);
            Assert.AreEqual(0.5, b.PDF(0));
            Assert.AreEqual(0.5, b.PDF(0.5));
            Assert.AreEqual(0.5, b.PDF(1));

            var b2 = new GeneralizedBeta(9, 1);
            Assert.AreEqual(0, b2.PDF(0));
            Assert.AreEqual(0.035156, b2.PDF(0.5), 1e-04);
            Assert.AreEqual(8.9999, b2.PDF(1), 1e-04);
            Assert.AreEqual(0, b2.PDF(-1));
            Assert.AreEqual(0, b2.PDF(2));
        }

        /// <summary>
        /// Verifying CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var b = new GeneralizedBeta(2,2,-10,10);
            Assert.AreEqual(b.CDF(-11), 0);
            Assert.AreEqual(b.CDF(11),1);

            var b2 = new GeneralizedBeta(9, 1);
            Assert.AreEqual(0, b2.CDF(0));
            Assert.AreEqual(0.001953125, b2.CDF(0.5));
            Assert.AreEqual(1, b2.CDF(1));
        }

        /// <summary>
        /// Verifying inverse CDF function.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var b = new GeneralizedBeta(1, 1);
            Assert.AreEqual(1, b.InverseCDF(1));

            var b2 = new GeneralizedBeta(9, 1,-10,10);
            Assert.AreEqual(-10, b2.InverseCDF(0));
            Assert.AreEqual(10, b2.InverseCDF(1));

            var b3 = new GeneralizedBeta(5, 100,0,10);
            Assert.AreEqual(0, b3.InverseCDF(0));
        }
    }
}
