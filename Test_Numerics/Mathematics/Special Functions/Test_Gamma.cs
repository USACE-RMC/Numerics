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
using Numerics.Mathematics.SpecialFunctions;
using System;

namespace Mathematics.SpecialFunctions
{
    /// <summary>
    /// Unit tests for the Gamma class
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    ///     Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <b> Description: </b>
    /// Some of these special functions are tested directly with MS Excel and with values provided by:
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// Others must be tested indirectly through the continuous distribution classes.
    /// </remarks>
    [TestClass]
    public class Test_Gamma
    {
        /// <summary>
        /// Test the Gamma Function (which uses the Stirling approximation)
        /// </summary>
        [TestMethod]
        public void Test_Function()
        {
            var testX = new double[] { 1d, 1.2d, 1.4d, 1.6d, 1.8d, 2d, 0.2d, 0.4d, 0.6d, 0.8d, -0.2d, -0.4d, -0.6d, -0.8d, 10d, 20d, 30d, 40d, 50d, 60d };
            var testValid = new double[] { 1d, 0.918168742399761d, 0.887263817503075d, 0.89351534928769d, 0.931383770980243d, 1d, 4.5908437119988d, 2.21815954375769d, 1.48919224881282d, 1.1642297137253d, -5.82114856862652d, -3.72298062203204d, -3.69693257292948d, -5.73855463999851d, 362880d, 121645100408832000d, 8.8417619937397E+30d, 2.03978820811974E+46d, 6.08281864034268E+62d, 1.3868311854569E+80d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.Function(testX[i]);
                Assert.AreEqual(Math.Abs(testValid[i] - testResults[i]) / testValid[i] < 0.01, true);
            }
        }

        /// <summary>
        /// Test the Lanczos approximation to the Gamma Function
        /// </summary>
        [TestMethod]
        public void Test_Lanczos()
        {
            var testX = new double[] { 1d, 1.2d, 1.4d, 1.6d, 1.8d, 2d, 0.2d, 0.4d, 0.6d, 0.8d, -0.2d, -0.4d, -0.6d, -0.8d, 10d, 20d, 30d, 40d, 50d, 60d };
            var testValid = new double[] { 1d, 0.918168742399761d, 0.887263817503075d, 0.89351534928769d, 0.931383770980243d, 1d, 4.5908437119988d, 2.21815954375769d, 1.48919224881282d, 1.1642297137253d, -5.82114856862652d, -3.72298062203204d, -3.69693257292948d, -5.73855463999851d, 362880d, 121645100408832000d, 8.8417619937397E+30d, 2.03978820811974E+46d, 6.08281864034268E+62d, 1.3868311854569E+80d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.Lanczos(testX[i]);
                Assert.AreEqual(Math.Abs(testValid[i] - testResults[i]) / testValid[i] < 0.01, true);
            }
        }

        /// <summary>
        /// Test the digamma function
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_Digamma()
        {
            double[] x = { 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 };

            double[] expected =
            {
                -0.5772156649015329,
                -0.4237549404110768,
                -0.2890398965921883,
                -0.1691908888667997,
                -0.06138454458511615,
                 0.03648997397857652,
                 0.1260474527734763,
                 0.208547874873494,
                 0.2849914332938615,
                 0.3561841611640597,
                 0.4227843350984671,
            };

            for (int i = 0; i < x.Length; i++)
            {
                double xi = x[i];
                double expectedi = expected[i];
                double actual = Gamma.Digamma(xi);

                Assert.AreEqual(expectedi, actual, 1e-10);
            }
        }

        /// <summary>
        /// Test the trigamma function
        /// </summary>
        /// <remarks>
        /// The trigamma function was tested against values calculated from WolframAlpha's "trigamma" function
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// <see href = "https://www.wolframalpha.com/input?i=trigamma%28x%29" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_Trigamma()
        {
            var testX = new double[] {1.2d, 1.4d, 1.6d, 1.8d, 0.2d, 0.4d, 0.6d, 0.8d };
            var testValid = new double[] {1.26738d, 1.02536d, 0.858432d, 0.736974d, 26.2674d, 7.27536d, 3.63621d, 2.29947d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.Trigamma(testX[i]);
                Assert.AreEqual(Math.Abs(testValid[i] - testResults[i]) / testValid[i] < 0.01, true);
            }
        }

        /// <summary>
        /// Test the function that returns the natural logarithm of the gamma function
        /// </summary>
        /// <remarks>
        /// This function was tested against values calculated from MS Excel's "GAMMALN( )" function
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// <see href = "https://learn.microsoft.com/en-us/office/vba/api/excel.worksheetfunction.gammaln" />
        ///  </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_LogGamma()
        {
            var testX = new double[] {1.2d, 1.4d, 1.6d, 1.8d, 0.2d, 0.4d, 0.6d, 0.8d, 10d, 20d, 30d, 40d, 50d, 60d};
            var testValid = new double[] { -0.08537409d, -0.119612914d, -0.112591766d, -0.071083873d,
                1.524063822d, 0.796677818d, 0.398233858d, 0.152059678d, 12.80182748d, 39.33988419d, 71.25703897d,
                106.6317603d, 144.5657439d, 184.5338289d};
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.LogGamma(testX[i]);
                Assert.AreEqual(Math.Abs(testValid[i] - testResults[i]) / testValid[i] < 0.01, true);
            }
        }

        /// <summary>
        /// Test the incomplete gamma integral
        /// </summary>
        /// <remarks>
        /// Test values of x and a, and exact values of the gamma function copied from the testing of the
        /// ACM TOMS algorithm 435 of the modified incomplete Gamma function originally written in FORTRAN77
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>        
        /// Wayne Fullerton,
        /// Algorithm 435: Modified Incomplete Gamma Function,
        /// Communications of the ACM,
        /// Volume 15, Number 11, November 1972, pages 993-995
        /// </description></item>
        /// <item><description>   
        /// <see href = "https://people.math.sc.edu/Burkardt/f77_src/toms435/toms435_prb_output.txt"/>
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_Incomplete()
        {
            var testX = new double[] { 0.03d, 0.3d ,1.5d, 0.075d, 0.75d, 3.5d, 0.1d, 1d, 5d, 0.1d, 1d, 5d,
                0.15d, 1.5d, 7d, 2.5d, 12d, 16d, 25d, 45d};
            var testA = new double[] { 0.1d, 0.1d, 0.1d, 0.5d, 0.5d, 0.5d, 1d, 1d, 1d, 1.1d, 1.1d, 1.1d,
                2d, 2d, 2d, 6d, 6d, 11d, 26d};
            var testValid = new double[] { 0.73823506d, 0.90835798d, 0.98865598d, 0.30146465d, 0.77932864d,
                0.99184901d, 0.95162585E-01d, 0.63212055d, 0.99326205d, 0.72059743E-01d, 0.58918095d,0.99153680d,
                0.10185828E-01d, 0.44217461d, 0.99270493d, 0.42021040E-01d, 0.97965896d, 0.92260396d, 0.44707859d};
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.Incomplete(testX[i], testA[i]);
                Assert.AreEqual(Math.Abs(testValid[i] - testResults[i]) / testValid[i] < 0.01, true);
            }

        }

        /// <summary>
        /// Test the upper incomplete regularized Gamma function
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_UpperIncomplete()
        {
            double expected, actual;

            actual = Gamma.UpperIncomplete(0.000000, 2);
            expected = 1.000000;
            Assert.AreEqual(expected, actual);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(0.250000, 2);
            expected = 0.017286;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(0.500000, 2);
            expected = 0.045500;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(0.750000, 2);
            expected = 0.085056;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(1.000000, 2);
            expected = 0.135335;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(1.250000, 2);
            expected = 0.194847;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(1.500000, 2);
            expected = 0.261464;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(1.750000, 2);
            expected = 0.332706;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(2.000000, 2);
            expected = 0.406006;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(2.250000, 2);
            expected = 0.478944;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(2.500000, 2);
            expected = 0.549416;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(2.750000, 2);
            expected = 0.615734;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));
        }

        /// <summary>
        /// Test the lower incomplete regularized Gamma function
        /// </summary>
        /// <remarks>
        /// Test values of x and a, and exact values of the lower incomplete gamma function copied from the testing of the
        /// "gammaLower_reg(  )" function from the Cephes Math Library
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Cephes Math Library Release 2.8: June, 2000
        /// </description></item>
        /// <item><description>
        /// <see href = "https://www.codecogs.com/library/maths/special/gamma/gamma_lower_reg.php"/>
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_LowerIncomplete()
        {
            double expected, actual;

            actual = Gamma.LowerIncomplete(0.000000, 2);
            expected = 1.000000;
            Assert.AreEqual(expected, actual);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(0.250000, 2);
            expected = 0.982714;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(0.500000, 2);
            expected = 0.954500;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(0.750000, 2);
            expected = 0.914944;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(1.000000, 2);
            expected = 0.864665;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(1.250000, 2);
            expected = 0.805153;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(1.500000, 2);
            expected = 0.738536;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(1.750000, 2);
            expected = 0.667294;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(2.000000, 2);
            expected = 0.593994;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(2.250000, 2);
            expected = 0.521056;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(2.500000, 2);
            expected = 0.450584;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(2.750000, 2);
            expected = 0.384266;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

        }

        /// <summary>
        /// Test the inverse of the lower regularized Gamma function
        /// </summary>
        /// <remarks>
        /// The inverse lower incomplete gamma function was test against values calculated from MATLAB's "gammaincinv( )" function
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// The MathWorks, Inc. (2022). MATLAB version: 9.13.0 (R2022b).
        /// </description></item>
        /// <item><description>
        /// <see href = "https://www.mathworks.com/help/matlab/ref/gammaincinv.html" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_InverseLowerIncomplete()
        {
            var testY = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d, 0.6d, 0.7d, 0.8d, 0.9d };
            var testValid = new double[] { 0.1054d, 0.8244d, 1.9138d, 3.2113d, 4.6709d, 6.2919d, 8.1110d, 10.2325d, 12.9947d};
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.InverseLowerIncomplete(i+1, testY[i]); ;
                Assert.AreEqual(testValid[i], testResults[i], 1e-4);
            }

        }

        /// <summary>
        /// Test the inverse of the upper incomplete regularized Gamma function
        /// </summary>
        /// <remarks>
        /// Part of this code was copied from the Accord Math Library
        /// <para> 
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_InverseUpperIncomplete()
        {
            var vec = new double[100];
            vec[0] = 0.1;
            for (int i = 1; i < 100; i++)
            {
                vec[i] = vec[i - 1] + 0.1;
            }

            foreach (double lambda in vec)
            {
                for (int i = 1; i < 100; i++)
                {
                    double x = Gamma.UpperIncomplete(lambda, i);
                    double j = Gamma.InverseUpperIncomplete(lambda, x);

                    Assert.IsTrue(Math.Abs(i - j) < 1e-2 * Math.Abs(j));
                }
            }
        }
    }
}
