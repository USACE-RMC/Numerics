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
using Numerics.Mathematics.RootFinding;
using System;

namespace Mathematics.RootFinding
{
    /// <summary>
    /// A class of various functions unit testing the Newton-Raphson Method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_NewtonRaphson
    {
        /// <summary>
        /// Testing with a quadratic function.
        /// </summary>
        [TestMethod()]
        public void Test_Quadratic()
        {
            double initial = 1.0;
            double X = NewtonRaphson.Solve(TestFunctions.Quadratic, TestFunctions.Quadratic_Deriv, initial);
            double trueX = Math.Sqrt(2);
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing robust Newton-Raphson with a quadratic function.
        /// </summary>
        [TestMethod()]
        public void Test_Robust_Quadratic()
        {
            double initial = 1.0;
            double lower = 0;
            double upper = 4;
            double X = NewtonRaphson.RobustSolve(TestFunctions.Quadratic, TestFunctions.Quadratic_Deriv, initial, lower, upper);
            double trueX = Math.Sqrt(2);
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a cubic function.
        /// </summary>
        [TestMethod()]
        public void Test_Cubic()
        {
            double initial = 1d;
            double X = NewtonRaphson.Solve(TestFunctions.Cubic, TestFunctions.Cubic_Deriv, initial);
            double trueX = 1.32472;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing robust Newton-Raphson with a cubic function.
        /// </summary>
        [TestMethod()]
        public void Test_Robust_Cubic()
        {
            double initial = 1.0;
            double lower = -1;
            double upper = 5;
            double X = NewtonRaphson.RobustSolve(TestFunctions.Cubic, TestFunctions.Cubic_Deriv, initial, lower, upper);
            double trueX = 1.32472d;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a trigonometric function.
        /// </summary>
        [TestMethod()]
        public void Test_Trigonometric()
        {
            double initial = 0.5;
            double X = NewtonRaphson.Solve(TestFunctions.Trigonometric, TestFunctions.Trigonometric_Deriv, initial);
            double trueX = 1.12191713d;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing robust Newton-Raphson with a trigonometric function.
        /// </summary>
        [TestMethod()]
        public void Test_Robust_Trigonometric()
        {
            double initial = 0.5;
            double lower = 0;
            double upper = Math.PI;
            double X = NewtonRaphson.RobustSolve(TestFunctions.Trigonometric, TestFunctions.Trigonometric_Deriv, initial, lower, upper);
            double trueX = 1.12191713d;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with an exponential function.
        /// </summary>
        [TestMethod()]
        public void Test_Exponential()
        {
            double initial = 1d;
            double X = NewtonRaphson.Solve(TestFunctions.Exponential, TestFunctions.Exponential_Deriv, initial);
            double trueX = 0.567143290d;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing robust Newton-Raphson with a exponential function.
        /// </summary>
        [TestMethod()]
        public void Test_Robust_Exponential()
        {
            double initial = 1.0;
            double lower = -2;
            double upper = 2;
            double X = NewtonRaphson.RobustSolve(TestFunctions.Exponential, TestFunctions.Exponential_Deriv, initial, lower, upper);
            double trueX = 0.567143290d;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a power function.
        /// </summary>
        [TestMethod()]
        public void Test_Power()
        {
            double initial = 0.2;
            double X = NewtonRaphson.Solve(TestFunctions.Power, TestFunctions.Power_Deriv, initial);
            double trueX = 1.0;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing robust Newton-Raphson with a power function.
        /// </summary>
        [TestMethod()]
        public void Test_Robust_Power()
        {
            double initial = 0.2;
            double lower = 0;
            double upper = 2;
            double X = NewtonRaphson.RobustSolve(TestFunctions.Power, TestFunctions.Power_Deriv, initial, lower, upper);
            double trueX = 1.0;
            Assert.AreEqual(X, trueX, 1E-5);
        }

    }
}
