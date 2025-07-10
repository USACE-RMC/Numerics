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
using Numerics.Mathematics;

namespace Mathematics.Differentiation
{
    /// <summary>
    /// Unit tests for the methods of numerical differentiation
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    /// Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// </remarks>
    [TestClass]
    public class Test_Differentiation
    {
        /// <summary>
        /// Test function of x^3
        /// </summary>
        public double FX(double x)
        {
            return Math.Pow(x, 3d);
        }

        /// <summary>
        /// Test function of e^x
        /// </summary>
        public double EX(double x)
        {
            return Math.Exp(x);
        }

        /// <summary>
        /// Test function of ln(x)
        /// </summary>
        public double LN(double x)
        {
            return Math.Log(x);
        }

        /// <summary>
        /// Test 2 parameter function
        /// </summary>
        public double FXY(double[] points)
        {
            double x = points[0];
            double y = points[1];
            return Math.Pow(x, 2) * Math.Pow(y, 3);
        }

        /// <summary>
        /// Test multi-parameter function
        /// </summary>
        public double FXYZ(double[] points)
        {
            double x = points[0];
            double y = points[1];
            double z = points[2];
            return Math.Pow(x, 3d) + Math.Pow(y, 4d) + Math.Pow(z, 5d);
        }

        /// <summary>
        /// Test multi-parameter function for Hessian
        /// </summary>
        public double FH(double[] points)
        {
            double x = points[0];
            double y = points[1];
            return Math.Pow(x, 3d) - 2 * x * y - Math.Pow(y, 6);
        }

        /// <summary>
        /// Test multi-parameter function for Jacobian
        /// </summary>
        public double FTXYZ(double t, double[] points)
        {
            double x = points[0];
            double y = points[1];
            double z = points[2];
            return t * Math.Pow(x, 3d) + Math.Pow(y, 4d) + Math.Pow(z, 5d);
        }

        /// <summary>
        /// Test derivative.
        /// </summary>
        [TestMethod()]
        public void Test_Derivative_FX()
        {
            double derivFX = NumericalDerivative.Derivative(FX, 2d);
            double ansFX = 12d;
            Assert.AreEqual(derivFX, ansFX, 0.000001d);
        
            double derivEX = NumericalDerivative.Derivative(EX, 4d);
            double ansEX = Math.Exp(4);
            Assert.AreEqual(derivEX, ansEX, 1E-4);
            
            double derivLN = NumericalDerivative.Derivative(LN, 2d);
            double ansLN = 0.5d;
            Assert.AreEqual(derivLN, ansLN, 1E-4);
        }

        /// <summary>
        /// Test gradient.
        /// </summary>
        [TestMethod()]
        public void Test_Gradient()
        {
            var deriv = NumericalDerivative.Gradient(FXY, new[] { 2d, 2d });
            double x = 2 * 2 * Math.Pow(2, 3);
            double y = 3 * Math.Pow(2, 2) * Math.Pow(2, 2);
            Assert.AreEqual(deriv[0], x, 0.000001d);
            Assert.AreEqual(deriv[1], y, 0.000001d);

            deriv = NumericalDerivative.Gradient(FXYZ, new[] { 2d, 2d, 2d });
            double xx = 3d * Math.Pow(2d, 2d);
            double yy = 4d * Math.Pow(2d, 3d);
            double zz = 5d * Math.Pow(2d, 4d);
            Assert.AreEqual(deriv[0], xx, 0.000001d);
            Assert.AreEqual(deriv[1], yy, 0.000001d);
            Assert.AreEqual(deriv[2], zz, 0.000001d);
        }


        /// <summary>
        /// Test Ridders' method.
        /// </summary>
        [TestMethod()]
        public void Test_RiddersMethod_FX()
        {
            double derivFX = NumericalDerivative.RiddersMethod(FX, 2d, out var err);
            double ansFX = 12d;
            Assert.AreEqual(derivFX, ansFX, 0.000001d);
        
            double derivEX = NumericalDerivative.RiddersMethod(EX, 4d, out err);
            double ansEX = Math.Exp(4);
            Assert.AreEqual(derivEX, ansEX, 1E-4);
        
            double derivLN = NumericalDerivative.RiddersMethod(LN, 2d, out err);
            double ansLN = 0.5d;
            Assert.AreEqual(derivLN, ansLN, 1E-4);
        }

        /// <summary>
        /// Test Jacobian.
        /// </summary>
        [TestMethod()]
        public void Test_Jacobian()
        {
            var jac = NumericalDerivative.Jacobian(FTXYZ, new[] { 4d, 6d, 8d}, new[] { 2d, 2d, 2d });
            var true_jac = new double[,] { { 12d * Math.Pow(2d, 2d), 4d * Math.Pow(2d, 3d), 5d * Math.Pow(2d, 4d) }, 
                                           { 18d * Math.Pow(2d, 2d), 4d * Math.Pow(2d, 3d), 5d * Math.Pow(2d, 4d) }, 
                                           { 24d * Math.Pow(2d, 2d), 4d * Math.Pow(2d, 3d), 5d * Math.Pow(2d, 4d) } };
            for (int i = 0; i < jac.GetLength(0); i++)
                for (int j = 0; j < jac.GetLength(1); j++)
                    Assert.AreEqual(jac[i, j], true_jac[i, j], 1E-4);

        }

        /// <summary>
        /// Test Hessian.
        /// </summary>
        [TestMethod()]
        public void Test_Hessian()
        {
            var hess = NumericalDerivative.Hessian(FH, new[] { 1d, 2d });
            var true_hess = new double[,] { { 6, -2 }, { -2, -480 } };
            Assert.AreEqual(hess[0,0], true_hess[0, 0], 1E-3);
            Assert.AreEqual(hess[0,1], true_hess[0, 1], 1E-3);
            Assert.AreEqual(hess[1,0], true_hess[1, 0], 1E-3);
            Assert.AreEqual(hess[1,1], true_hess[1, 1], 1E-3);

            hess = NumericalDerivative.Hessian(FXYZ, new[] { 2d, 2d, 2d });
            true_hess = new double[,] { { 12d, 0, 0 }, { 0, 48d, 0}, { 0, 0, 160d } };
            for (int i = 0; i < hess.GetLength(0); i++)
                for (int j = 0; j < hess.GetLength(1); j++)
                    Assert.AreEqual(hess[i, j], true_hess[i, j], 1E-3);

        }
    }
}
