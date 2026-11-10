using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for Simpson's rule of integration
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_SimpsonsRule
    {
        /// <summary>
        /// Test Simpson's rule with a one dimensional function
        /// </summary>
        [TestMethod]
        public void Test_FX3()
        {
            var sr = new SimpsonsRule(Integrands.FX3, 0d, 1d);
            sr.Integrate();
            double result = sr.Result;
            double trueResult = 0.25d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Simpson's rule with the cos(x) function
        /// </summary>
        [TestMethod]
        public void Test_Cosine()
        {
            var sr = new SimpsonsRule(Integrands.Cosine, -1, 1);
            sr.Integrate();
            double result = sr.Result;
            double trueResult = 1.6829419d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Simpson's rule with the sin(x) function
        /// </summary>
        [TestMethod]
        public void Test_Sine()
        {
            var sr = new SimpsonsRule(Integrands.Sine, 0, 1);
            sr.Integrate();
            double result = sr.Result;
            double trueResult = 0.459697694131d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Simpson's rule with a 2nd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXX()
        {
            var sr = new SimpsonsRule(Integrands.FXX, 0, 2);
            sr.Integrate();
            double result = sr.Result;
            double trueResult = 57;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Simpson's rule with a 3rd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXXX()
        {
            var sr = new SimpsonsRule(Integrands.FXXX, 0, 2);
            sr.Integrate();
            double result = sr.Result;
            double trueResult = 89;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Simpson's rule with the Gamma function
        /// </summary>
        [TestMethod]
        public void Test_Gamma()
        {
            var gamma = new GammaDistribution(10, 5);
            var sr = new SimpsonsRule(x => x * gamma.PDF(x), gamma.InverseCDF(1E-16), gamma.InverseCDF(1 - 1E-16));
            sr.Integrate();
            double result = sr.Result;
            double trueResult = 50;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Simpson's rule with the CVaR function
        /// </summary>
        [TestMethod]
        public void Test_CVaR()
        {
            var ln = new LnNormal(10, 2);
            double alpha = 0.99;
            var sr = new SimpsonsRule(x => x * ln.PDF(x), ln.InverseCDF(alpha), ln.InverseCDF(1 - 1E-16));
            sr.Integrate();
            double result = sr.Result / (1 - alpha);
            double trueResult = Math.Exp(ln.Mu + 0.5 * ln.Sigma * ln.Sigma) / (1d - alpha) * (1 - Normal.StandardCDF(Normal.StandardZ(alpha) - ln.Sigma));
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Reversed or equal limits throw ArgumentOutOfRangeException.
        /// </summary>
        [TestMethod]
        public void Test_Constructor_Guards()
        {
            Func<double, double> f = x => x;
            Assert.Throws<ArgumentNullException>(() => new SimpsonsRule(null!, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new SimpsonsRule(f, 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new SimpsonsRule(f, 1, 1));
        }

        /// <summary>
        /// Non-success exits report the latest Simpson estimate (4 T(2n) - T(n)) / 3 instead of NaN.
        /// </summary>
        [TestMethod]
        public void Test_NonSuccessExit_ReportsBestEstimate()
        {
            var iterations = new SimpsonsRule(F, 0, 2) { MaxIterations = 4 };
            iterations.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, iterations.Status);
            Assert.AreEqual((4 * Trapezoid(8) - Trapezoid(4)) / 3, iterations.Result, 1E-15);

            var budget = new SimpsonsRule(F, 0, 2) { MaxFunctionEvaluations = 50 };
            budget.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumFunctionEvaluationsReached, budget.Status);
            Assert.AreEqual((4 * Trapezoid(64) - Trapezoid(32)) / 3, budget.Result, 1E-15);

            // After a single pass there is nothing to extrapolate, so the trapezoid is reported.
            var single = new SimpsonsRule(F, 0, 2) { MaxIterations = 1 };
            single.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, single.Status);
            Assert.AreEqual(Trapezoid(1), single.Result, 1E-15);
        }

        /// <summary>
        /// Composite trapezoid with n panels for exp(-x) cos(3x) on [0, 2], the independent reference.
        /// </summary>
        private static double Trapezoid(int n)
        {
            double h = 2.0 / n, sum = 0.5 * (F(0) + F(2));
            for (int i = 1; i < n; i++) sum += F(i * h);
            return h * sum;
        }

        /// <summary>
        /// The smooth test integrand exp(-x) cos(3x).
        /// </summary>
        private static double F(double x) => Math.Exp(-x) * Math.Cos(3 * x);
    }
}
