using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics;
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

        /// <summary>
        /// A smooth integrand still reports Success, and its result keeps the bit pattern produced
        /// before forced acceptances were counted. The analytic value is
        /// (1 + exp(-2) (3 sin 6 - cos 6)) / 10.
        /// </summary>
        [TestMethod]
        public void Test_SmoothIntegrand_SuccessAndBitsUnchanged()
        {
            var adapt = new AdaptiveGaussKronrod(x => Math.Exp(-x) * Math.Cos(3 * x), 0, 2);
            adapt.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
            Assert.AreEqual(0x3FB35E863514D492, BitConverter.DoubleToInt64Bits(adapt.Result));
            Assert.AreEqual(21, adapt.FunctionEvaluations);
            Assert.AreEqual((1 + Math.Exp(-2) * (3 * Math.Sin(6) - Math.Cos(6))) / 10, adapt.Result, 1E-15);
        }

        /// <summary>
        /// Intervals accepted at the width or depth limit without meeting the tolerance report
        /// MaximumIterationsReached; the result value is unchanged by the limit status.
        /// </summary>
        [TestMethod]
        public void Test_ForcedAcceptance_ReportsMaximumIterations()
        {
            // The integral of 1/x on [0, 1] diverges; refinement stops at the width limit.
            var divergent = new AdaptiveGaussKronrod(x => 1.0 / x, 0, 1);
            divergent.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, divergent.Status);
            Assert.AreEqual(0x40463986E70A561C, BitConverter.DoubleToInt64Bits(divergent.Result));

            // A unit step at 0.3 integrates to 0.7, but three levels cannot resolve the jump.
            var step = new AdaptiveGaussKronrod(x => x < 0.3 ? 0.0 : 1.0, 0, 1) { MaxDepth = 3 };
            step.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, step.Status);
            Assert.AreEqual(0x3FE671E68616C21C, BitConverter.DoubleToInt64Bits(step.Result));
            Assert.IsGreaterThan(1E-8, Math.Abs(step.Result - 0.7));
        }

        /// <summary>
        /// Non-finite and reversed or equal limits throw ArgumentOutOfRangeException.
        /// </summary>
        [TestMethod]
        public void Test_Constructor_Guards()
        {
            Func<double, double> f = x => x;
            Assert.Throws<ArgumentNullException>(() => new AdaptiveGaussKronrod(null!, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(f, 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(f, 1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(f, double.NaN, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(f, 0, double.PositiveInfinity));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(f, double.NegativeInfinity, 0));
        }

        /// <summary>
        /// A NaN integrand stops with Failure as soon as a NaN interval is accepted instead of
        /// spending the 10,000,000-evaluation budget, and throws when ReportFailure is true.
        /// </summary>
        [TestMethod]
        public void Test_NaNIntegrand_FailsFast()
        {
            var silent = new AdaptiveGaussKronrod(x => double.NaN, 0, 1) { ReportFailure = false };
            silent.Integrate();
            Assert.AreEqual(IntegrationStatus.Failure, silent.Status);
            Assert.IsTrue(double.IsNaN(silent.Result));
            Assert.AreEqual(2247, silent.FunctionEvaluations);

            var loud = new AdaptiveGaussKronrod(x => double.NaN, 0, 1);
            Assert.Throws<ArithmeticException>(() => loud.Integrate());
            Assert.AreEqual(IntegrationStatus.Failure, loud.Status);

            // An infinite value at the center node is a Failure, not an infinite Success.
            var infinite = new AdaptiveGaussKronrod(x => 1.0 / Math.Sqrt(Math.Abs(x)), -1, 1) { ReportFailure = false };
            infinite.Integrate();
            Assert.AreEqual(IntegrationStatus.Failure, infinite.Status);
            Assert.AreEqual(21, infinite.FunctionEvaluations);
        }

        /// <summary>
        /// A NaN on an interval that is subdivided further does not stop integration: sin(x - 0.5)/(x - 0.5)
        /// is 0/0 at the center node of [0, 1], yet the children are finite. The result keeps its bit
        /// pattern and matches 2 Si(0.5) from scipy.special.sici.
        /// </summary>
        [TestMethod]
        public void Test_RemovableNaN_ResultUnchanged()
        {
            var adapt = new AdaptiveGaussKronrod(x => Math.Sin(x - 0.5) / (x - 0.5), 0, 1);
            adapt.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
            Assert.AreEqual(0x3FEF8F126A7A3CFB, BitConverter.DoubleToInt64Bits(adapt.Result));
            Assert.AreEqual(63, adapt.FunctionEvaluations);
            Assert.AreEqual(0.9862148360861335, adapt.Result, 1E-15);
        }

        /// <summary>
        /// The default method is local bisection: results, evaluation counts, statuses and standard errors are
        /// pinned bit for bit (exact double literals in round-trip form).
        /// </summary>
        [TestMethod]
        public void Test_DefaultMethod_IsBitwiseUnchanged()
        {
            Assert.AreEqual(AdaptiveGaussKronrodMethod.LocalBisection, new AdaptiveGaussKronrod(x => x, 0, 1).Method);
            var cases = new (Func<double, double> F, double Rel, double Result, int Evaluations, IntegrationStatus Status, double Error)[]
            {
                (x => 1 / Math.Sqrt(x), 1E-8, 1.9999999945252451, 1911, IntegrationStatus.Success, 8.507218787094601E-09),
                (Math.Log, 1E-8, -0.999999998372131, 819, IntegrationStatus.Success, 9.31325069129541E-09),
                (x => Math.Pow(x, -0.9), 1E-8, 9.88353611413269, 2247, IntegrationStatus.MaximumIterationsReached, 0.023782977988624873),
                (x => 1 / ((x - 0.3) * (x - 0.3) + 1E-4), 1E-8, 309.3986915124149, 357, IntegrationStatus.Success, 1.181691002672109E-08),
                (x => Math.Abs(x - 1.0 / 3), 1E-8, 0.27777777834913675, 399, IntegrationStatus.Success, 4.611653090953075E-09),
                (x => { double t = (x - 0.3) * 60d; return Math.Exp(-t * t) + 0.1 * x; }, 1E-10, 0.07954089751509194, 231, IntegrationStatus.Success, 8.116236945723367E-10),
            };
            foreach (var c in cases)
            {
                var adapt = new AdaptiveGaussKronrod(c.F, 0, 1) { RelativeTolerance = c.Rel };
                adapt.Integrate();
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(c.Result), BitConverter.DoubleToInt64Bits(adapt.Result));
                Assert.AreEqual(c.Evaluations, adapt.FunctionEvaluations);
                Assert.AreEqual(c.Status, adapt.Status);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(c.Error), BitConverter.DoubleToInt64Bits(adapt.StandardError));
                Assert.AreEqual(0, adapt.Iterations);
            }

            var minDepth = new AdaptiveGaussKronrod(x => x * x, 0, 1) { MinDepth = 2 };
            minDepth.Integrate();
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(0.33333333333333337), BitConverter.DoubleToInt64Bits(minDepth.Result));
            Assert.AreEqual(147, minDepth.FunctionEvaluations);
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(1.5540063044689707E-17), BitConverter.DoubleToInt64Bits(minDepth.StandardError));

            // The default keeps the base tolerance range: a zero relative tolerance is rejected.
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(Math.Exp, 0, 1) { RelativeTolerance = 0 }.Integrate());
        }
    }
}
