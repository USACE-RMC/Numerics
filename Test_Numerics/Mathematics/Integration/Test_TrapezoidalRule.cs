using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit test for the Trapezoidal rule of integration
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_TrapezoidalRule
    {
        /// <summary>
        /// Test Trapezoidal rule with a one dimensional function
        /// </summary>
        [TestMethod]
        public void Test_FX3()
        {
            var tr = new TrapezoidalRule(Integrands.FX3, 0d, 1d);
            tr.Integrate();
            double result = tr.Result;
            double trueResult = 0.25d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Trapezoidal rule with the cos(x) function
        /// </summary>
        [TestMethod]
        public void Test_Cosine()
        {
            var tr = new TrapezoidalRule(Integrands.Cosine, -1, 1);
            tr.Integrate();
            double result = tr.Result;
            double trueResult = 1.6829419d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Trapezoidal rule with the sin(x) function
        /// </summary>
        [TestMethod]
        public void Test_Sine()
        {
            var tr = new TrapezoidalRule(Integrands.Sine, 0, 1);
            tr.Integrate();
            double result = tr.Result;
            double trueResult = 0.459697694131d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Trapezoidal rule with a 2nd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXX()
        {
            var tr = new TrapezoidalRule(Integrands.FXX, 0, 2);
            tr.Integrate();
            double result = tr.Result;
            double trueResult = 57;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Trapezoidal rule with a 3rd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXXX()
        {
            var tr = new TrapezoidalRule(Integrands.FXXX, 0, 2);
            tr.Integrate();
            double result = tr.Result;
            double trueResult = 89;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Trapezoidal rule with the Gamma function
        /// </summary>
        [TestMethod]
        public void Test_Gamma()
        {
            var gamma = new GammaDistribution(10, 5);
            var tr = new TrapezoidalRule(x => x * gamma.PDF(x), gamma.InverseCDF(1E-16), gamma.InverseCDF(1 - 1E-16));
            tr.Integrate();
            double result = tr.Result;
            double trueResult = 50;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test Trapezoidal rule with the CVaR function
        /// </summary>
        [TestMethod]
        public void Test_CVaR()
        {
            var ln = new LnNormal(10, 2);
            double alpha = 0.99;
            var tr = new TrapezoidalRule(x => x * ln.PDF(x), ln.InverseCDF(alpha), ln.InverseCDF(1 - 1E-16));
            tr.Integrate();
            double result = tr.Result / (1 - alpha);
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
            Assert.Throws<ArgumentNullException>(() => new TrapezoidalRule(null!, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new TrapezoidalRule(f, 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new TrapezoidalRule(f, 1, 1));
        }

        /// <summary>
        /// Non-success exits report the latest refinement instead of NaN: four iterations stop at the
        /// 8-panel trapezoid, and a 50-evaluation budget stops at the 64-panel trapezoid.
        /// </summary>
        [TestMethod]
        public void Test_NonSuccessExit_ReportsBestEstimate()
        {
            var iterations = new TrapezoidalRule(F, 0, 2) { MaxIterations = 4 };
            iterations.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, iterations.Status);
            Assert.AreEqual(Trapezoid(8), iterations.Result, 1E-15);

            var budget = new TrapezoidalRule(F, 0, 2) { MaxFunctionEvaluations = 50 };
            budget.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumFunctionEvaluationsReached, budget.Status);
            Assert.AreEqual(65, budget.FunctionEvaluations);
            Assert.AreEqual(Trapezoid(64), budget.Result, 1E-15);

            // A throwing integrand reports the last completed refinement.
            int calls = 0;
            var failing = new TrapezoidalRule(x => ++calls > 3 ? throw new InvalidOperationException("Integrand failure.") : F(x), 0, 2) { ReportFailure = false };
            failing.Integrate();
            Assert.AreEqual(IntegrationStatus.Failure, failing.Status);
            Assert.AreEqual(Trapezoid(2), failing.Result, 1E-15);
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
