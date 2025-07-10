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
using Numerics.Data;
using Numerics.Distributions;
using Numerics.Functions;

namespace Functions
{
    /// <summary>
    /// Unit tests for the function classes
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
    /// </remarks>
    [TestClass]
    public class Test_Functions
    {
        /// <summary>
        /// Test the linear function
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Jiang L, Linton O, Tang H, Zhang Y (2023). drcarlate: Improving 
        /// Estimation Efficiency in CAR with Imperfect Compliance. R package version 1.2.0,
        /// https://CRAN.R-project.org/package=drcarlate
        /// </remarks>
        [TestMethod]
        public void Test_Linear_Function()
        {
            // Default constructor with alpha = 0 and beta = 1
            var func0 = new LinearFunction();
            double y0 = func0.Function(6);
            double valid0 = 6;
            Assert.AreEqual(valid0, y0, 1E-6);

            double alpha = -2;
            double beta = 5;
            double sigma = 3;

            var func1 = new LinearFunction(alpha, beta);
            double y1 = func1.Function(6);
            double valid1 = (5 * 6) + -2;
            Assert.AreEqual(valid1, y1, 1E-6);

            var func2 = new LinearFunction(alpha, beta, sigma);
            double y2 = func2.Function(6);
            Assert.AreEqual(valid1, y2, 1E-6);

            func2.ConfidenceLevel = 0.75;
            double y3 = func2.Function(6);
            // Found using R's norminv() function for ϵ from the "drcarlate" package
            double valid3 = 30.0234692505882;
            Assert.AreEqual(valid3, y3, 1E-6);
        }

        /// <summary>
        /// Test that the inverse of the linear function against the linear function itself
        /// </summary>
        [TestMethod]
        public void Test_Linear_Function_Inverse()
        {
            var func = new LinearFunction(10, 0.5, 20);
            double y = func.Function(400);
            double x = func.InverseFunction(y);
            Assert.AreEqual(400, x, 1E-6);

            func.ConfidenceLevel = 0.75;
            double yy = func.Function(400);
            double xx = func.InverseFunction(yy);
            Assert.AreEqual(400, xx, 1E-6);
        }

        /// <summary>
        /// Test the power function
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Jiang L, Linton O, Tang H, Zhang Y (2023). drcarlate: Improving 
        /// Estimation Efficiency in CAR with Imperfect Compliance. R package version 1.2.0,
        /// https://CRAN.R-project.org/package=drcarlate
        /// </remarks>
        [TestMethod]
        public void Test_Power_Function()
        {
            // Default constructor with alpha = 1, beta = 1.5, and xi = 0
            var func0 = new PowerFunction();
            var y0 = func0.Function(6);
            var valid0 = 1 * Math.Pow(6 - 0, 1.5);
            Assert.AreEqual(valid0, y0, 1E-6);

            double alpha = 5;
            double beta = 2;
            double sigma = 3;
            double xi = 0;

            var func1 = new PowerFunction(alpha, beta, xi);
            double y1 = func1.Function(6);
            double valid1 = alpha * Math.Pow(6 - xi, beta);
            Assert.AreEqual(valid1, y1, 1E-6);

            var func2 = new PowerFunction(alpha, beta, xi, sigma);
            double y2 = func2.Function(6);
            Assert.AreEqual(valid1, y2, 1E-6);

            func2.ConfidenceLevel = 0.75;
            double y3 = func2.Function(6);
            // Found using R's norminv() function for ϵ from the "drcarlate" package
            double valid3 = 1361.61408399941;
            Assert.AreEqual(valid3, y3, 1E-6);
        }

        /// <summary>
        /// Test the inverse of the power function against the power function itself
        /// </summary>
        [TestMethod]
        public void Test_Power_Function_Inverse()
        {
            var func = new PowerFunction(10, 2, 0,0.1);
            double y = func.Function(400);
            double x = func.InverseFunction(y);
            Assert.AreEqual(400, x, 1E-6);

            func.ConfidenceLevel = 0.75;
            double yy = func.Function(400);
            double xx = func.InverseFunction(yy);
            Assert.AreEqual(400, xx, 1E-6);
        }

        /// <summary>
        /// Test the inverse power function
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Jiang L, Linton O, Tang H, Zhang Y (2023). drcarlate: Improving 
        /// Estimation Efficiency in CAR with Imperfect Compliance. R package version 1.2.0,
        /// https://CRAN.R-project.org/package=drcarlate
        /// </remarks>
        [TestMethod]
        public void Test_InversePower_Function()
        {
            double alpha = 5;
            double beta = 2;
            double sigma = 3;
            double xi = 0;

            var func = new PowerFunction(alpha, beta, xi);
            func.IsInverse = true;
            double y = func.Function(6);
            double valid = Math.Sqrt(6 / alpha) + xi;
            Assert.AreEqual(valid, y, 1E-6);

            var func2 = new PowerFunction(alpha, beta, xi, sigma);
            func2.IsInverse = true;
            double y2 = func2.Function(6);
            Assert.AreEqual(valid, y2, 1E-6);

            func2.ConfidenceLevel = 0.75;
            double y3 = func2.Function(6);
            // Found using R's norminv() function for ϵ from the "drcarlate" package
            double valid3 = 0.398290417772997;
            Assert.AreEqual(valid3, y3, 1E-6);
        }

        /// <summary>
        /// Test the inverse of the inverse of the power function against the inverse power function itself
        /// </summary>
        [TestMethod]
        public void Test_InversePower_Function_Inverse()
        {
            var func = new PowerFunction(10, 2, 0, 0.1);
            func.IsInverse = true;
            double y = func.Function(6);
            double x = func.InverseFunction(y);
            Assert.AreEqual(6, x, 1E-6);

            func.ConfidenceLevel = 0.75;
            double yy = func.Function(6);
            double xx = func.InverseFunction(yy);
            Assert.AreEqual(6, xx, 1E-6);
        }

        /// <summary>
        /// Test the tabular function
        /// </summary>
        [TestMethod]
        public void Test_Tabular_Function()
        {
            var XArray = new double[] { 50, 100, 150, 200, 250 };
            var YArray = new UnivariateDistributionBase[] { new Deterministic(100), new Deterministic(200), new Deterministic(300), new Deterministic(400), new Deterministic(500) };
            
            var opd = new UncertainOrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending, UnivariateDistributionType.Deterministic);
            var func = new TabularFunction(opd) { XTransform = Transform.Logarithmic, YTransform = Transform.None };

            // Given X
            double X = 50.0;
            double Y = func.Function(X);
            Assert.AreEqual(Y, 100.0);

            // Given Y
            double Y2 = 100d;
            double X2 = func.InverseFunction(Y2);
            Assert.AreEqual(X, 50.0);

            // Given X - Interpolation
            double X3 = 75.0d;
            double Y3 = func.Function(X3);
            // Basic interpolation formula: y = y1 + (y2 - y1)/(x2 - x1) * (x-x1)
            // Here the x values are logarithmic transformed
            double validY = 100 + (200 - 100) / (Math.Log(100) - Math.Log(50)) * (Math.Log(75) - Math.Log(50));
            Assert.AreEqual(validY, Y3, 1E-6);

            // Given Y - Interpolation
            double Y4 = 100 + (200 - 100) / (Math.Log(100) - Math.Log(50)) * (Math.Log(75) - Math.Log(50));
            double X4 = func.InverseFunction(Y4);
            Assert.AreEqual(75, X4, 1E-6);
        }

    }
}
