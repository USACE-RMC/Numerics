using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the Adaptive Gauss-Lobatto integration method
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_AdaptiveGaussLobatto
    {
        /// <summary>
        /// Test the Adaptive Gauss Lobatto method with a one dimensional function
        /// </summary>
        [TestMethod]
        public void Test_FX3()
        {
            var adapt = new AdaptiveGaussLobatto(Integrands.FX3, 0d, 1d);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 0.25d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss Lobatto method with the cos(x) function
        /// </summary>
        [TestMethod]
        public void Test_Cosine()
        {
            var adapt = new AdaptiveGaussLobatto(Integrands.Cosine, -1, 1);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 1.6829419d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss Lobatto method with the sin(x) function
        /// </summary>
        [TestMethod]
        public void Test_Sine()
        {
            var adapt = new AdaptiveGaussLobatto(Integrands.Sine, 0, 1);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 0.459697694131d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss Lobatto method with a 2nd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXX()
        {
            var adapt = new AdaptiveGaussLobatto(Integrands.FXX, 0, 2);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 57;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss Lobatto method with a 3rd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXXX()
        {
            var adapt = new AdaptiveGaussLobatto(Integrands.FXXX, 0, 2);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 89;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss Lobatto method with the Gamma function
        /// </summary>
        [TestMethod]
        public void Test_Gamma()
        {
            var gamma = new GammaDistribution(10, 5);
            var adapt = new AdaptiveGaussLobatto(x => x * gamma.PDF(x), gamma.InverseCDF(1E-16), gamma.InverseCDF(1 - 1E-16));
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 50;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss Lobatto method with the CVaR function
        /// </summary>
        [TestMethod]
        public void Test_CVaR()
        {
            var ln = new LnNormal(10, 2);
            double alpha = 0.99;
            var adapt = new AdaptiveGaussLobatto(x => x * ln.PDF(x), ln.InverseCDF(alpha), ln.InverseCDF(1 - 1E-16));
            adapt.Integrate();
            double result = adapt.Result / (1 - alpha);
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
            Assert.Throws<ArgumentNullException>(() => new AdaptiveGaussLobatto(null!, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussLobatto(f, 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussLobatto(f, 1, 1));
        }

        /// <summary>
        /// A NaN integrand stops with Failure once a NaN interval is accepted instead of spending
        /// the evaluation budget (10,007,748 evaluations before), and throws when ReportFailure is true.
        /// A smooth integrand keeps its result bits.
        /// </summary>
        [TestMethod]
        public void Test_NaNIntegrand_FailsFast()
        {
            var silent = new AdaptiveGaussLobatto(x => double.NaN, 0, 1) { ReportFailure = false };
            silent.Integrate();
            Assert.AreEqual(IntegrationStatus.Failure, silent.Status);
            Assert.IsTrue(double.IsNaN(silent.Result));
            Assert.AreEqual(1573, silent.FunctionEvaluations);

            var loud = new AdaptiveGaussLobatto(x => double.NaN, 0, 1);
            Assert.Throws<ArithmeticException>(() => loud.Integrate());

            var smooth = new AdaptiveGaussLobatto(x => Math.Exp(-x) * Math.Cos(3 * x), 0, 2);
            smooth.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, smooth.Status);
            Assert.AreEqual(0x3FB35E863514D493, BitConverter.DoubleToInt64Bits(smooth.Result));
        }
    }
}
