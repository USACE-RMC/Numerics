using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the Adaptive Simpson's rule 
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_AdaptiveSimpsonsRule
    {
        /// <summary>
        /// Test the Adaptive Simpson's rule with a one dimensional function
        /// </summary>
        [TestMethod]
        public void Test_FX3()
        {
            var adapt = new AdaptiveSimpsonsRule(Integrands.FX3, 0d, 1d);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 0.25d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Simpson's rule with the cos(x) function
        /// </summary>
        [TestMethod]
        public void Test_Cosine()
        {
            var adapt = new AdaptiveSimpsonsRule(Integrands.Cosine, -1, 1);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 1.6829419d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Simpson's rule with the sin(x) function
        /// </summary>
        [TestMethod]
        public void Test_Sine()
        {
            var adapt = new AdaptiveSimpsonsRule(Integrands.Sine, 0, 1);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 0.459697694131d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Simpson's rule with a 2nd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXX()
        {
            var adapt = new AdaptiveSimpsonsRule(Integrands.FXX, 0, 2);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 57;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Simpson's rule with a 3rd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXXX()
        {
            var adapt = new AdaptiveSimpsonsRule(Integrands.FXXX, 0, 2);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 89;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Simpson's rule with the Gamma function
        /// </summary>
        [TestMethod]
        public void Test_Gamma()
        {
            var gamma = new GammaDistribution(10, 5);
            var adapt = new AdaptiveSimpsonsRule(x => x * gamma.PDF(x), gamma.InverseCDF(1E-16), gamma.InverseCDF(1 - 1E-16));
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 50;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Simpson's rule with the CVaR function
        /// </summary>
        [TestMethod]
        public void Test_CVaR()
        {
            var ln = new LnNormal(10, 2);
            double alpha = 0.99;
            var adapt = new AdaptiveSimpsonsRule(x => x * ln.PDF(x), ln.InverseCDF(alpha), ln.InverseCDF(1 - 1E-16));
            adapt.Integrate();
            double result = adapt.Result / (1 - alpha);
            double trueResult = Math.Exp(ln.Mu + 0.5 * ln.Sigma * ln.Sigma) / (1d - alpha) * (1 - Normal.StandardCDF(Normal.StandardZ(alpha) - ln.Sigma));
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// A smooth integrand still reports Success with the result bits produced before forced
        /// acceptances were counted; the analytic value is (1 + exp(-2) (3 sin 6 - cos 6)) / 10.
        /// </summary>
        [TestMethod]
        public void Test_SmoothIntegrand_SuccessAndBitsUnchanged()
        {
            var asr = new AdaptiveSimpsonsRule(x => Math.Exp(-x) * Math.Cos(3 * x), 0, 2);
            asr.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, asr.Status);
            Assert.AreEqual(0x3FB35E8635146FC1, BitConverter.DoubleToInt64Bits(asr.Result));
            Assert.AreEqual(169, asr.FunctionEvaluations);
            Assert.AreEqual((1 + Math.Exp(-2) * (3 * Math.Sin(6) - Math.Cos(6))) / 10, asr.Result, 1E-11);
        }

        /// <summary>
        /// Intervals accepted at the depth limit without meeting the tolerance report
        /// MaximumIterationsReached with an unchanged result.
        /// </summary>
        [TestMethod]
        public void Test_ForcedAcceptance_ReportsMaximumIterations()
        {
            // A unit step at 0.3 integrates to 0.7.
            var step = new AdaptiveSimpsonsRule(x => x < 0.3 ? 0.0 : 1.0, 0, 1) { MaxDepth = 3 };
            step.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, step.Status);
            Assert.AreEqual(0x3FE6444444444444, BitConverter.DoubleToInt64Bits(step.Result));
            Assert.IsGreaterThan(1E-8, Math.Abs(step.Result - 0.7));

            // The integral of 1/x on [0, 1] diverges, and f(0) is infinite, so the accepted estimate is not finite
            // and integration stops with Failure.
            var divergent = new AdaptiveSimpsonsRule(x => 1.0 / x, 0, 1) { ReportFailure = false };
            divergent.Integrate();
            Assert.AreEqual(IntegrationStatus.Failure, divergent.Status);

            var reported = new AdaptiveSimpsonsRule(x => 1.0 / x, 0, 1);
            Assert.ThrowsExactly<ArithmeticException>(() => reported.Integrate());
            Assert.AreEqual(IntegrationStatus.Failure, reported.Status);
        }

        /// <summary>
        /// Reversed or equal limits throw ArgumentOutOfRangeException.
        /// </summary>
        [TestMethod]
        public void Test_Constructor_Guards()
        {
            Func<double, double> f = x => x;
            Assert.Throws<ArgumentNullException>(() => new AdaptiveSimpsonsRule(null!, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveSimpsonsRule(f, 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveSimpsonsRule(f, 1, 1));
        }
    }
}
