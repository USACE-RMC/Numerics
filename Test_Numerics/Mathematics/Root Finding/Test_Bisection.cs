using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.RootFinding;
using System;

namespace Mathematics.RootFinding
{
    /// <summary>
    /// A class of various functions unit testing the Bisection Method.
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
    public class Test_Bisection
    {

        /// <summary>
        /// Testing with a quadratic function.
        /// </summary>
        [TestMethod()]
        public void Test_Quadratic()
        {
            double initial = 1.0;
            double lower = 0;
            double upper = 4;
            double X = Bisection.Solve(TestFunctions.Quadratic, initial, lower, upper);
            double trueX = Math.Sqrt(2);
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a cubic function.
        /// </summary>
        [TestMethod()]
        public void Test_Cubic()
        {
            double initial = 1.0;
            double lower = -1;
            double upper = 5;
            double X = Bisection.Solve(TestFunctions.Cubic, initial, lower, upper);
            double trueX = 1.32472;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a trigonometric function.
        /// </summary>
        [TestMethod()]
        public void Test_Trigonometric()
        {
            double initial = 0.5;
            double lower = 0;
            double upper = Math.PI;
            double X = Bisection.Solve(TestFunctions.Trigonometric, initial, lower, upper);
            double trueX = 1.12191713d;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with an exponential function.
        /// </summary>
        [TestMethod()]
        public void Test_Exponential()
        {
            double initial = 1.0;
            double lower = -2;
            double upper = 2;
            double X = Bisection.Solve(TestFunctions.Exponential, initial, lower, upper);
            double trueX = 0.567143290;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a power function.
        /// </summary>
        [TestMethod()]
        public void Test_Power()
        {
            double initial = 0.2;
            double lower = 0;
            double upper = 2;
            double X = Bisection.Solve(TestFunctions.Power, initial, lower, upper);
            double trueX = 1.0;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        [TestMethod()]
        public void Test_BisectionEdge1()
        {
            double initial = -1d;
            double lower = 1d;
            double upper = 5d;
            var ex = Assert.Throws<Exception>(() =>
            {
                double X = Bisection.Solve(TestFunctions.Cubic, initial, lower, upper);
            });
        }

        [TestMethod()]
        public void Test_BisectionEdge2()
        {
            double initial = 1d;
            double lower = 5;
            double upper = 0d;
            var ex = Assert.Throws<Exception>(() =>
            {
                double X = Bisection.Solve(TestFunctions.Cubic, initial, lower, upper);
            });
        }

        [TestMethod()]
        public void Test_BisectionEdge3()
        {
            double initial = 3d;
            double lower = 2;
            double upper = 5d;
            var ex = Assert.Throws<Exception>(() =>
            {
                double X = Bisection.Solve(TestFunctions.Cubic, initial, lower, upper);
            });
        }

    }
}
