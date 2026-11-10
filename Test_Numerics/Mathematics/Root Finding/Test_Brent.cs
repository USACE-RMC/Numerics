using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.RootFinding;
using System;

namespace Mathematics.RootFinding
{
    /// <summary>
    /// A class of various functions unit testing the Brent Method.
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
    public class Test_Brent
    {
        /// <summary>
        /// Testing with a quadratic function.
        /// </summary>
        [TestMethod()]
        public void Test_Quadratic()
        {
            double lower = 0;
            double upper = 4;
            double X = Brent.Solve(TestFunctions.Quadratic, lower, upper);
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
            double X = Brent.Solve(TestFunctions.Cubic, lower, upper);
            double trueX = 1.32472;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a trigonometric function.
        /// </summary>
        [TestMethod()]
        public void Test_Trigonometric()
        {
            double lower = 0;
            double upper = Math.PI;
            double X = Brent.Solve(TestFunctions.Trigonometric, lower, upper);
            double trueX = 1.12191713d;
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
            double X = Brent.Solve(TestFunctions.Exponential, lower, upper);
            double trueX = 0.567143290;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// Testing with a power function.
        /// </summary>
        [TestMethod()]
        public void Test_Power()
        {
            double lower = 0;
            double upper = 2;
            double X = Brent.Solve(TestFunctions.Power, lower, upper);
            double trueX = 1.0;
            Assert.AreEqual(X, trueX, 1E-5);
        }

        /// <summary>
        /// A NaN at one bound, as for a function undefined at the edge of its domain, is tolerated: the other
        /// bound decides the orientation. (Brent.Solve keeps its pre-2.3 NaN handling.)
        /// </summary>
        [TestMethod()]
        public void Test_NaNAtOneBoundIsTolerated()
        {
            Func<double, double> lowerNaN = x => x == 0d ? double.NaN : x - 0.3;
            Func<double, double> upperNaN = x => x == 1d ? double.NaN : x - 0.3;
            Assert.AreEqual(0.3, Bisection.Solve(lowerNaN, 0.5, 0, 1), 1E-8);
            Assert.AreEqual(0.3, Bisection.Solve(upperNaN, 0.5, 0, 1), 1E-8);
            Assert.AreEqual(0.3, NewtonRaphson.RobustSolve(lowerNaN, x => 1d, 0.5, 0, 1), 1E-8);
        }

        /// <summary>The scipy.optimize.brentq oracle cases: the functions, brackets and roots below.</summary>
        private static readonly Func<double, double>[] BrentqFunctions =
        {
                x => x * x - 2.0,
                x => x * x * x - x - 1.0,
                x => x * x * x - 2.0 * x - 5.0,
                x => (x - 1.0) * (x - 1.0) * (x - 1.0),
                x => (x * x - 3.0) / (x + 1.0),
                x => Math.Exp(-x) - x,
                x => x * x * x * x * x * x * x * x * x - 1E-6
        };

        private static readonly double[] BrentqLower = { 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0 };
        private static readonly double[] BrentqUpper = { 2.0, 2.0, 3.0, 3.0, 5.0, 1.0, 1.0 };

        /// <summary>
        /// The root is bitwise equal to scipy.optimize.brentq (SciPy 1.18.1, maxiter=1000) at its default tolerances
        /// (xtol=2e-12, rtol=4 * 2^-52), with xtol=1e-8, with xtol=1e-12 and rtol=1e-10, and with xtol=1e-300 and rtol=1e-14.
        /// </summary>
        [TestMethod()]
        public void Test_BitwiseScipyBrentq()
        {
            double[] defaultTolerance = { 1.4142135623731364, 1.324717957244746, 2.094551481542327, 0.999999999999876, 1.7320508075688776, 0.5671432904097838, 0.21544346900320108 };
            double[] looseTolerance = { 1.4142135623731364, 1.3247179571960475, 2.094551481542327, 1.0000000004382756, 1.7320508083674189, 0.5671432904178182, 0.21544346900320108 };
            double[] relative1E10 = { 1.4142135623731364, 1.3247179572627834, 2.094551481542327, 0.9999999999981231, 1.7320508075688776, 0.5671432904178182, 0.21544346900320108 };
            double[] relative1E14 = { 1.414213562373095, 1.324717957244746, 2.094551481542327, 0.9999999999999978, 1.7320508075688776, 0.5671432904097838, 0.21544346900318836 };
            for (int i = 0; i < BrentqFunctions.Length; i++)
            {
                var f = BrentqFunctions[i];
                double a = BrentqLower[i], b = BrentqUpper[i];
                AssertBitwise(defaultTolerance[i], Brent.Solve(f, a, b, 2E-12), i);
                AssertBitwise(defaultTolerance[i], Brent.Solve(f, a, b, 2E-12, Brent.MinimumRelativeTolerance), i);
                AssertBitwise(looseTolerance[i], Brent.Solve(f, a, b, 1E-8), i);
                AssertBitwise(relative1E10[i], Brent.Solve(f, a, b, 1E-12, 1E-10), i);
                AssertBitwise(relative1E14[i], Brent.Solve(f, a, b, 1E-300, 1E-14), i);
            }
        }

        /// <summary>The minimum relative tolerance is 4 * 2^-52, the default rtol of scipy.optimize.brentq.</summary>
        [TestMethod()]
        public void Test_MinimumRelativeTolerance()
        {
            double expected = 4d * Math.Pow(2d, -52);
            double actual = Brent.MinimumRelativeTolerance;
            Assert.AreEqual(expected, actual, 0d);
        }

        /// <summary>A zero at a bound is returned before the sign test, as brentq.c does.</summary>
        [TestMethod()]
        public void Test_ZeroAtBound()
        {
            Assert.AreEqual(0d, Brent.Solve(x => x * (x + 1d), 0d, 3d), 0d);
            Assert.AreEqual(3d, Brent.Solve(x => x - 3d, 1d, 3d), 0d);
        }

        /// <summary>Invalid bounds, tolerances and brackets throw, and a failed search throws or returns the last estimate.</summary>
        [TestMethod()]
        public void Test_InvalidInputsAndFailure()
        {
            Func<double, double> f = x => x * x - 2d;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Brent.Solve(f, 2d, 0d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Brent.Solve(f, 0d, 2d, -1E-8));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Brent.Solve(f, 0d, 2d, double.NaN));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Brent.Solve(f, 0d, 2d, 1E-8, 4E-16));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Brent.Solve(f, 0d, 2d, 1E-8, double.NaN));
            Assert.ThrowsExactly<ArgumentException>(() => Brent.Solve(f, 2d, 3d));
            Assert.ThrowsExactly<ArgumentException>(() => Brent.Solve(f, 0d, 2d, 1E-300, 1E-14, 3));
            double last = Brent.Solve(f, 0d, 2d, 1E-300, 1E-14, 3, false);
            Assert.IsTrue(last > 0d && last < 2d);
        }

        /// <summary>
        /// Inputs SciPy rejects stay accepted: a zero absolute tolerance converges on the relative term alone (to the
        /// rtol=4 * 2^-52 root of scipy.optimize.brentq with xtol=1e-300, 1.414213562373095), and a negative iteration count runs
        /// no iterations, so the call throws or returns the upper bound.
        /// </summary>
        [TestMethod()]
        public void Test_ZeroToleranceAndNegativeIterationsAreAccepted()
        {
            Func<double, double> f = x => x * x - 2d;
            AssertBitwise(1.414213562373095, Brent.Solve(f, 0d, 2d, 0d), 0);
            Assert.ThrowsExactly<ArgumentException>(() => Brent.Solve(f, 0d, 2d, 1E-8, -1));
            Assert.AreEqual(2d, Brent.Solve(f, 0d, 2d, 1E-8, -1, false), 0d);
        }

        private static void AssertBitwise(double expected, double actual, int index)
        {
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(actual), $"case {index}: expected {expected:R}, actual {actual:R}");
        }
    }
}
