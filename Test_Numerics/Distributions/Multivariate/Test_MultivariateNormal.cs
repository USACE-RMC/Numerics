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
using Numerics.Sampling;

namespace Distributions.Multivariate
{
    /// <summary>
    /// Unit tests for the Multivariate Normal distribution. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_MultivariateNormal
    {

        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_MultivariateNormalDist()
        {

            var true_mean = new[] { 4d, 2d };
            var true_mode = new[] { 4d, 2d };
            var true_median = new[] { 4d, 2d };
            var true_variance = new[] { 0.3d, 0.7d };
            var true_stdDev = new[] { Math.Sqrt(0.3d), Math.Sqrt(0.7d) };
            double true_pdf1 = 0.000000018917884164743237d;
            double true_pdf2 = 0.35588127170858852d;
            double true_pdf3 = 0.000000000036520107734505265d;
            double true_cdf = 0.033944035782101548d;
            var MultiN = new MultivariateNormal(new[] { 4d, 2d }, new[,] { { 0.3d, 0.1d }, { 0.1d, 0.7d } });
            var mean = MultiN.Mean;
            var mode = MultiN.Mode;
            var median = MultiN.Median;
            var variance = MultiN.Variance;
            var stdev = MultiN.StandardDeviation;
            double pdf1 = MultiN.PDF(new[] { 2d, 5d });
            double pdf2 = MultiN.PDF(new[] { 4d, 2d });
            double pdf3 = MultiN.PDF(new[] { 3d, 7d });
            double cdf = MultiN.CDF(new[] { 3d, 5d });
            for (int i = 0; i <= MultiN.Dimension - 1; i++)
            {
                Assert.AreEqual(mean[i], true_mean[i], 0.0001d);
                Assert.AreEqual(median[i], true_median[i], 0.0001d);
                Assert.AreEqual(mode[i], true_mode[i], 0.0001d);
                Assert.AreEqual(stdev[i], true_stdDev[i], 0.0001d);
            }

            Assert.AreEqual(pdf1, true_pdf1, 0.0001d);
            Assert.AreEqual(pdf2, true_pdf2, 0.0001d);
            Assert.AreEqual(pdf3, true_pdf3, 0.0001d);
            Assert.AreEqual(cdf, true_cdf, 0.0001d);

        }

        /// <summary>
        /// Verification against Fortran routine
        /// http://www.math.wsu.edu/faculty/genz/software/fort77/mvndstpack.f
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_Fortran()
        {
            var tester = new MultivariateNormal(5);
            tester.MVNUNI = new MersenneTwister(12345);
            int N = 5;
            int MAXPTS = 50;
            double ABSEPS = 0.00005;
            double RELEPS = 0;
            double[] Lower = new double[] { 0, 0, 1.7817, 1.4755, 1.5949 };
            double[] Upper = new double[] { 0, 1.5198, 1.7817, 1.4755, 1.5949 };
            int[] INFIN = new int[] { 1, 2, 1, 1, 0 };
            double[] CORREL = new double[] { -0.707107, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, .5, .5 };

            double[] expectedErrors = new double[] { 0.00000811074362075292, 0.00000480583149442263, 0.00000660161142196953 };
            double[] expectedValues = new double[] { 0.00286779150981026, 0.00111850297940743, 0.00397918930026649 };
            double[] expectedInforms = new double[] { 0, 0, 0 };

            for (int i = 0; i < 3; i++)
            {
                double ERROR = 0;
                double VAL = 0;
                int INFORM = 0;

                tester.MVNDST(N, Lower, Upper, INFIN, CORREL, MAXPTS, ABSEPS, RELEPS, ref ERROR, ref VAL, ref INFORM);

                Assert.IsTrue(Math.Abs(ERROR - expectedErrors[i]) < 1E-5 && Math.Abs(VAL - expectedValues[i]) < 1E-5 && Math.Abs(INFORM - expectedInforms[i]) == 0);
                INFIN[0] = INFIN[0] - 1;
            }

        }

        /// <summary>
        /// Verified against R "mvtnorm" package
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_R()
        {

            double r = -0.33;
            var mean = new double[] { 0, 0, 0, 0 };
            var covar = new double[,]
                            {{ 1, r, r, r },
                             { r, 1, r, r },
                             { r, r, 1, r },
                             { r, r, r, 1 }};
            var mvn = new MultivariateNormal(mean, covar) { MVNUNI = new MersenneTwister(12345) };

            // AB
            var p = mvn.CDF(new[] { Normal.StandardZ(0.25), Normal.StandardZ(0.35), double.PositiveInfinity, double.PositiveInfinity });
            Assert.AreEqual(p, 0.05011069, 1E-4);

            // AC
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), double.PositiveInfinity, Normal.StandardZ(0.5), double.PositiveInfinity });
            Assert.AreEqual(p, 0.0827451, 1E-4);

            // AD
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), double.PositiveInfinity, double.PositiveInfinity, Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.0827451, 1E-4);

            // BC
            p = mvn.CDF(new[] { double.PositiveInfinity, Normal.StandardZ(0.35), Normal.StandardZ(0.5), double.PositiveInfinity });
            Assert.AreEqual(p, 0.1254504, 1E-4);

            // BD
            p = mvn.CDF(new[] { double.PositiveInfinity, Normal.StandardZ(0.35), double.PositiveInfinity, Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.1254504, 1E-4);

            // CD
            p = mvn.CDF(new[] { double.PositiveInfinity, double.PositiveInfinity, Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.1964756, 1E-4);

            // ABC
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), Normal.StandardZ(0.35), Normal.StandardZ(0.5), double.PositiveInfinity });
            Assert.AreEqual(p, 0.005960125, 1E-4);

            // ABD
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), Normal.StandardZ(0.35), double.PositiveInfinity, Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.005964513, 1E-4);

            // ACD
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), double.PositiveInfinity, Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.0128066, 1E-4);

            // BCD
            p = mvn.CDF(new[] { double.PositiveInfinity, Normal.StandardZ(0.35), Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.02324389, 1E-4);

            // ABCD
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), Normal.StandardZ(0.35), Normal.StandardZ(0.5), Normal.StandardZ(0.5)});
            Assert.AreEqual(p, 3.593582e-13, 1E-4);
        }

        /// <summary>
        /// Verified against R "mvtnorm" package.
        /// Test 1: Perfectly Negative correlation. 
        /// The matrix must be positive semi-definite. 
        /// So, the smallest allowable negative value is -1/(D-1) + an offset for machine double precision. 
        /// For simplicity, I offset by 0.01. 
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_R_PerfectNegative()
        {
            var mean = new double[] { 0, 0, 0 };
            var covar = new double[,]
            { { 1, -0.49,-0.49 },
              { -0.49, 1,-0.49 },
              { -0.49, -0.49, 1 }};
            var mvn = new MultivariateNormal(mean, covar) { MVNUNI = new MersenneTwister(12345) };

            var p = mvn.CDF(new[] { Normal.StandardZ(0.5), Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.002740932, 1E-4);
        }

        /// <summary>
        /// Verified against R "mvtnorm" package.
        /// Test 2: Perfectly positive correlation. 
        /// Again, I offset my 0.01 to keep it positive definite.
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_R_PerfectPositive()
        {

            var mean = new double[] { 0, 0, 0 };
            var covar = new double[,]
            { { 1, 0.99, 0.99 },
              { 0.99, 1, 0.99 },
              { 0.99, 0.99, 1 }};
            var mvn = new MultivariateNormal(mean, covar) { MVNUNI = new MersenneTwister(12345) };

            var p = mvn.CDF(new[] { Normal.StandardZ(0.5), Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.4661416, 1E-4);

        }

        /// <summary>
        /// Verified against R "mvtnorm" package.
        /// Test 3: Independent correlation. 
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_R_Independent()
        {

            var mean = new double[] { 0, 0, 0 };
            var covar = new double[,]
            { { 1, 0, 0 },
              { 0, 1, 0 },
              { 0, 0, 1 }};
            var mvn = new MultivariateNormal(mean, covar) { MVNUNI = new MersenneTwister(12345) };

            var p = mvn.CDF(new[] { Normal.StandardZ(0.5), Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(p, 0.125, 1E-4);

        }
         
    }
}
