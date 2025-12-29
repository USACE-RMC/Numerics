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
    /// Testing the Logistic distribution algorithm.
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
    public class Test_Logistic
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 5.1.1 Tippecanoe River Near Delphi, Indiana (Station 43) Data
        private double[] sample = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

        /// <summary>
        /// Verification of Logistic Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 9.1.1 page 295.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Logistic_MOM_Fit()
        {
            var LO = new Logistic();
            LO.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double x = LO.Xi;
            double a = LO.Alpha;
            double true_x = 12665d;
            double true_a = 2596.62d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x );
            Assert.IsLessThan(0.01d, (a - true_a) / true_a );
        }

        /// <summary>
        /// Verification of Logistic Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 9.1.1 page 295.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Logistic_MLE_Fit()
        {
            var LO = new Logistic();
            LO.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double x = LO.Xi;
            double a = LO.Alpha;
            double true_x = 12628.59d;
            double true_a = 2708.64d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x);
            Assert.IsLessThan(0.01d, (a - true_a) / true_a);
        }

        /// <summary>
        /// Test the quantile function for the Logistic Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 9.1.2 page 297.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Logistic_Quantile()
        {
            var LO = new Logistic(12665d, 2596.62d);
            double q100 = LO.InverseCDF(0.99d);
            double true_100 = 24597d;
            Assert.IsLessThan(0.01d, (q100 - true_100) / true_100);
            double p = LO.CDF(q100);
            double true_p = 0.99d;
            Assert.AreEqual(true_p, p);
        }

        /// <summary>
        /// Test the standard error for the Logistic Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 9.1.3 page 300.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Logistic_StandardError()
        {

            // Method of Moments
            var LO = new Logistic(12665d, 2596.62d);
            double se100 = Math.Sqrt(LO.QuantileVariance(0.99d, 48, ParameterEstimationMethod.MethodOfMoments));
            double true_se100 = 1684d;
            Assert.IsLessThan(0.01d, (se100 - true_se100) / true_se100);

            // Maximum Likelihood
            LO = new Logistic(12628.59d, 2708.64d);
            se100 = Math.Sqrt(LO.QuantileVariance(0.99d, 48, ParameterEstimationMethod.MaximumLikelihood));
            true_se100 = 1648d;
            Assert.IsLessThan(0.01d, (se100 - true_se100) / true_se100);
        }

        /// <summary>
        /// Verifying distribution is being created with inputs.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var LO = new Logistic();
            Assert.AreEqual(0, LO.Xi);
            Assert.AreEqual(0.1, LO.Alpha);

            var LO2 = new Logistic(1, 1);
            Assert.AreEqual(1,LO2.Xi);
            Assert.AreEqual(1, LO2.Alpha);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var LO = new Logistic(double.NaN,double.NaN);
            Assert.IsFalse(LO.ParametersValid);

            var LO2 = new Logistic(double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(LO2.ParametersValid);

            var LO3 = new Logistic(1, 0);
            Assert.IsFalse(LO3.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString().
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var LO = new Logistic();
            Assert.AreEqual("Location (ξ)", LO.ParametersToString[0, 0]);
            Assert.AreEqual("Scale (α)", LO.ParametersToString[1, 0]);
            Assert.AreEqual("0", LO.ParametersToString[0, 1]);
            Assert.AreEqual("0.1", LO.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Logistic();
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var LO = new Logistic();
            Assert.AreEqual(0, LO.Mean);

            var LO2 = new Logistic(1,1);
            Assert.AreEqual(1, LO2.Mean);
        }

        /// <summary>
        /// Testing median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var LO = new Logistic();
            Assert.AreEqual(0, LO.Median);

            var LO2 = new Logistic(1, 1);
            Assert.AreEqual(1, LO2.Median);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var LO = new Logistic();
            Assert.AreEqual(0, LO.Mode);

            var LO2 = new Logistic(1, 1);
            Assert.AreEqual(1, LO2.Mode);
        }

        /// <summary>
        /// Testing Standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var LO = new Logistic();
            Assert.AreEqual(0.18137, LO.StandardDeviation,  1E-04);

            var LO2 = new Logistic(1, 1);
            Assert.AreEqual(1.81379, LO2.StandardDeviation,  1e-04);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var LO = new Logistic();
            Assert.AreEqual(0, LO.Skewness);

            var LO2 = new Logistic(1, 1);
            Assert.AreEqual(0, LO2.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var LO = new Logistic();
            Assert.AreEqual(4.2, LO.Kurtosis);

            var LO2 = new Logistic(1, 1);
            Assert.AreEqual(4.2, LO2.Kurtosis);
        }

        /// <summary>
        /// Testing range for distribution is negative infinity to positive infinity.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var LO = new Logistic();
            Assert.AreEqual(double.NegativeInfinity,LO.Minimum);
            Assert.AreEqual(double.PositiveInfinity, LO.Maximum);

            var LO2 = new Logistic(1, 1);
            Assert.AreEqual(double.NegativeInfinity, LO2.Minimum);
            Assert.AreEqual(double.PositiveInfinity, LO2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var LO = new Logistic(5,2);
            Assert.AreEqual(0.00332, LO.PDF(-5),  1e-04);
            Assert.AreEqual(0.03505, LO.PDF(0),  1e-04);
            Assert.AreEqual(0.125, LO.PDF(5));
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var LO = new Logistic(5,2);
            Assert.AreEqual(0.00669, LO.CDF(-5),  1e-05);
            Assert.AreEqual(0.07585, LO.CDF(0),  1e-04);
            Assert.AreEqual(0.5, LO.CDF(5));
        }

        /// <summary>
        /// Testing inverse cdf method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var LO = new Logistic(5, 2);
            Assert.AreEqual(double.NegativeInfinity, LO.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity, LO.InverseCDF(1));
            Assert.AreEqual(3.3054, LO.InverseCDF(0.3),  1e-04);
        }
    }
}
