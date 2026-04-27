using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.RootFinding;
using System;

namespace Mathematics.RootFinding
{
    /// <summary>
    /// A class of various functions unit testing the Secant Method.
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
    public class Test_Secant
    {

        /// <summary>
        /// Testing with a quadratic function.
        /// </summary>
        [TestMethod()]
        public void Test_Quadratic()
        {
            double lower = 0;
            double upper = 4;
            double X = Secant.Solve(TestFunctions.Quadratic, lower, upper);
            double trueX = Math.Sqrt(2);
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a cubic function.
        /// </summary>
        [TestMethod()]
        public void Test_Cubic()
        {
            double lower = -1;
            double upper = 5;
            double X = Secant.Solve(TestFunctions.Cubic, lower, upper);
            double trueX = 1.32472;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with an exponential function.
        /// </summary>
        [TestMethod()]
        public void Test_Exponential()
        {
            double lower = -2;
            double upper = 2;
            double X = Secant.Solve(TestFunctions.Exponential, lower, upper);
            double trueX = 0.567143290;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing edge case where the function is discontinuous within the interval        
        /// </summary>
        public double Undefined(double x)
        {
            double F = 1 / x;
            return F;
        }

        [TestMethod()]
        public void Test_Edge()
        {
            double lower = -6d;
            double upper = 5d;
            var ex = Assert.Throws<Exception>(() =>
            {
                double X = Secant.Solve(Undefined, lower, upper);
                double trueX = Math.Sqrt(2);
                Assert.AreEqual(X, trueX, 1E-4);
            });
        }

    }
}
