using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.Integration;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the Adaptive Gauss-Kronrod integration method
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_AdaptiveGaussKronrod
    {
        /// <summary>
        /// Test the Adaptive Gauss-Kronrod method with a one dimensional function
        /// </summary>
        [TestMethod]
        public void Test_FX3()
        {
            var adapt = new AdaptiveGaussKronrod(Integrands.FX3, 0d, 1d);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 0.25d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss-Kronrod method with the cos(x) function
        /// </summary>
        [TestMethod]
        public void Test_Cosine()
        {
            var adapt = new AdaptiveGaussKronrod(Integrands.Cosine, -1, 1);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 1.6829419d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss-Kronrod method with the sin(x) function
        /// </summary>
        [TestMethod]
        public void Test_Sine()
        {
            var adapt = new AdaptiveGaussKronrod(Integrands.Sine, 0, 1);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 0.459697694131d;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss-Kronrod method with a 2nd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXX()
        {
            var adapt = new AdaptiveGaussKronrod(Integrands.FXX, 0, 2);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 57;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss-Kronrod method with a 3rd order polynomial
        /// </summary>
        [TestMethod]
        public void Test_FXXX()
        {
            var adapt = new AdaptiveGaussKronrod(Integrands.FXXX, 0, 2);
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 89;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss-Kronrod method with the Gamma function
        /// </summary>
        [TestMethod]
        public void Test_Gamma()
        {
            var gamma = new GammaDistribution(10, 5);
            var adapt = new AdaptiveGaussKronrod(x => x * gamma.PDF(x), gamma.InverseCDF(1E-16), gamma.InverseCDF(1 - 1E-16));
            adapt.Integrate();
            double result = adapt.Result;
            double trueResult = 50;
            Assert.AreEqual(result, trueResult, 1E-3);
        }

        /// <summary>
        /// Test the Adaptive Gauss-Kronrod method with the CVaR function
        /// </summary>
        [TestMethod]
        public void Test_CVaR()
        {
            var ln = new LnNormal(10, 2);
            double alpha = 0.99;
            var adapt = new AdaptiveGaussKronrod(x => x * ln.PDF(x), ln.InverseCDF(alpha), ln.InverseCDF(1 - 1E-16));
            adapt.Integrate();
            double result = adapt.Result / (1 - alpha);
            double trueResult = Math.Exp(ln.Mu + 0.5 * ln.Sigma * ln.Sigma) / (1d - alpha) * (1 - Normal.StandardCDF(Normal.StandardZ(alpha) - ln.Sigma));
            Assert.AreEqual(result, trueResult, 1E-3);
        }

    }
}
