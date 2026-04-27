using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.RootFinding;
using System;

namespace Mathematics.RootFinding
{
    /// <summary>
    /// A class of various functions unit testing the Bracketing Method.
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
    [TestClass()]
    public class Test_Bracket
    {

        /// <summary>
        /// Test with a quadratic function. 
        /// </summary>
        [TestMethod()]
        public void Test_Quadratic()
        {
            double lower = 0;
            double upper = 4;
            bool X = Brent.Bracket(TestFunctions.Quadratic, ref lower, ref upper, out double f1, out double f2);
            Assert.IsTrue(X);
        }

        /// <summary>
        /// Test with a cubic function. 
        /// </summary>
        [TestMethod()]
        public void Test_Cubic()
        {
            double lower = -1;
            double upper = 5;
            bool X = Brent.Bracket(TestFunctions.Cubic, ref lower, ref upper, out double f1, out double f2);
            Assert.IsTrue(X);
        }

        /// <summary>
        /// Test with a trigonometric function. 
        /// </summary>
        [TestMethod()]
        public void Test_Trigonometric()
        {
            double lower = 0;
            double upper = Math.PI;
            bool X = Brent.Bracket(TestFunctions.Trigonometric, ref lower, ref upper, out double f1, out double f2);
            Assert.IsTrue(X);
        }

        /// <summary>
        /// Test with an exponential function.
        /// </summary>
        [TestMethod()]
        public void Test_Exponential()
        {
            double lower = -2;
            double upper = 2;
            bool X = Brent.Bracket(TestFunctions.Exponential, ref lower, ref upper, out double f1, out double f2);
            Assert.IsTrue(X);
        }

        /// <summary>
        /// Test with a power function.
        /// </summary>
        [TestMethod()]
        public void Test_Power()
        {
            double lower = 0;
            double upper = 2;
            bool X = Brent.Bracket(TestFunctions.Power, ref lower, ref upper, out double f1, out double f2);
            Assert.IsTrue(X);
        }

        /// <summary>
        /// Test bad bracket.
        /// </summary>
        [TestMethod()]
        public void Test_BracketEdge()
        {
            double lower = 1;
            double upper = 1;
            var ex = Assert.Throws<Exception>(() =>
            {
                bool x = Brent.Bracket(
                    TestFunctions.Quadratic,
                    ref lower,
                    ref upper,
                    out double f1,
                    out double f2);
            });
        }
    }
}