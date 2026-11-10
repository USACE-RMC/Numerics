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

        /// <summary>
        /// Tests the first bisection edge case.
        /// </summary>
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

        /// <summary>
        /// Tests the second bisection edge case.
        /// </summary>
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

        /// <summary>
        /// Tests the third bisection edge case.
        /// </summary>
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

        /// <summary>
        /// The search covers the whole bracket regardless of the first guess. Halving that started at
        /// the first guess would return 2.5 for a guess of 2.5 for x^2 - 4 on [0, 3].
        /// </summary>
        [TestMethod()]
        public void Test_FirstGuessDoesNotLimitTheSearch()
        {
            Func<double, double> f = x => x * x - 4d;
            Assert.AreEqual(2d, Bisection.Solve(f, 2.5, 0, 3), 1E-8);
            Assert.AreEqual(2d, Bisection.Solve(f, 0.5, 0, 3), 1E-8);
            Assert.AreEqual(0.1, Bisection.Solve(x => x - 0.1, 0.5, 0, 1), 1E-8);
            // Decreasing function: f(lowerBound) > 0.
            Assert.AreEqual(2d, Bisection.Solve(x => 4d - x * x, 0.1, 0, 3), 1E-8);
        }

        /// <summary>
        /// NaN function values are rejected instead of being treated as a sign.
        /// </summary>
        [TestMethod()]
        public void Test_NaNFunctionValues()
        {
            Assert.ThrowsExactly<ArgumentException>(() => Bisection.Solve(x => double.NaN, 0.5, 0, 1));
            // The bounds bracket a sign change, but the function is NaN on (0.25, 0.75).
            Func<double, double> f = x => x > 0.25 && x < 0.75 ? double.NaN : x - 0.6;
            Assert.ThrowsExactly<ArgumentException>(() => Bisection.Solve(f, 0.5, 0, 1));
            Assert.IsTrue(double.IsNaN(Bisection.Solve(f, 0.5, 0, 1, reportFailure: false)));
        }
    }
}
