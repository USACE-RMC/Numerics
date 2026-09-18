using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the Miser algorithm
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_AdaptiveSimpsonsRule2D
    {

        /// <summary>
        /// Test the 2D Adaptive Simpson's Rule algorithm with the Pi function
        /// </summary>
        [TestMethod()]
        public void Test_PI()
        {
            var asr2D = new AdaptiveSimpsonsRule2D(Integrands.PI2D, -1, 1, -1, 1 );
            asr2D.Integrate();
            var result = asr2D.Result;
            double trueResult = 3.14;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the 2D Adaptive Simpson's Rule algorithm with the sum of two normal distributions.
        /// </summary>
        [TestMethod()]
        public void Test_SumOfTwoNormals()
        {
            var asr2D = new AdaptiveSimpsonsRule2D(Integrands.SumOfNormals2D, 1E-15, 1- 1E-15, 1E-15, 1- 1E-15);
            asr2D.Integrate();
            var result = asr2D.Result;
            double trueResult = 40;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the 2D Adaptive Simpson's Rule algorithm with X + Y integrand.
        /// </summary>
        [TestMethod()]
        public void Test_XPlusY()
        {
            // Define the function f(x, y) = x + y
            Func<double, double, double> func = (x, y) => x + y;

            // Set up the 2D Adaptive Simpson's Rule with the bounds [-1, 1] for both x and y
            var asr2D = new AdaptiveSimpsonsRule2D(func, -1, 1, -1, 1);

            // Run the integration process
            asr2D.Integrate();

            // Get the result of the integration
            var result = asr2D.Result;

            // The exact result is 0
            double trueResult = 0;

            // Assert that the result is within a small error margin of 0
            Assert.AreEqual(trueResult, result, 1E-5);
        }

        /// <summary>
        /// Test the 2D Adaptive Simpson's Rule algorithm with X^2 + Y^2 integrand.
        /// </summary>
        [TestMethod()]
        public void Test_XSquaredPlusYSquared()
        {
            // Define the function f(x, y) = x^2 + y^2
            Func<double, double, double> func = (x, y) => x * x + y * y;

            // Set up the 2D Adaptive Simpson's Rule with the bounds [-1, 1] for both x and y
            var asr2D = new AdaptiveSimpsonsRule2D(func, -1, 1, -1, 1);

            // Run the integration process
            asr2D.Integrate();

            // Get the result of the integration
            var result = asr2D.Result;

            // The exact result is 2.666...
            double trueResult = 2.6666667;

            // Assert that the result is within a small error margin of 2
            Assert.AreEqual(trueResult, result, 1E-6);
        }

        /// <summary>
        /// Polynomials at or below Simpson's cubic degree integrate exactly: the integral of
        /// x^3*y^3 over [0,2] x [0,1] is 4 * (1/4) = 1.
        /// </summary>
        [TestMethod]
        public void Test_PolynomialExactness()
        {
            Func<double, double, double> func = (x, y) => x * x * x * y * y * y;
            var asr2D = new AdaptiveSimpsonsRule2D(func, 0, 2, 0, 1);
            asr2D.Integrate();
            Assert.AreEqual(1d, asr2D.Result, 1E-12);
            Assert.AreEqual(IntegrationStatus.Success, asr2D.Status);
        }

        /// <summary>
        /// A separable product f(x)*g(y) must match the product of two one-dimensional adaptive
        /// Gauss-Kronrod reference passes, and the closed form (e^1.5 - 1)*sin(1).
        /// </summary>
        [TestMethod]
        public void Test_SeparableProduct_MatchesTwo1DPasses()
        {
            var asr2D = new AdaptiveSimpsonsRule2D((x, y) => Math.Exp(x) * Math.Cos(y), 0, 1.5, 0, 1);
            asr2D.Integrate();

            var agkX = new AdaptiveGaussKronrod(Math.Exp, 0, 1.5);
            agkX.Integrate();
            var agkY = new AdaptiveGaussKronrod(Math.Cos, 0, 1);
            agkY.Integrate();

            double product = agkX.Result * agkY.Result;
            double exact = (Math.Exp(1.5) - 1d) * Math.Sin(1d);
            Assert.AreEqual(product, asr2D.Result, 1E-6 * Math.Abs(product));
            Assert.AreEqual(exact, asr2D.Result, 1E-6 * Math.Abs(exact));
        }

        /// <summary>
        /// Genz oscillatory family, c = (5, 5), w1 = 0.3, against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzOscillatory()
        {
            var asr2D = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzOscillatory(x, y, 5, 5, 0.3), 0, 1, 0, 1);
            asr2D.Integrate();
            double exact = Integrands.GenzOscillatoryExact(5, 5, 0.3);
            Assert.AreEqual(exact, asr2D.Result, 1E-5 * Math.Abs(exact));
        }

        /// <summary>
        /// Genz product-peak family, c = (10, 10), w = (0.4, 0.6), against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzProductPeak()
        {
            var asr2D = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzProductPeak(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            asr2D.Integrate();
            double exact = Integrands.GenzProductPeakExact(10, 10, 0.4, 0.6);
            Assert.AreEqual(exact, asr2D.Result, 1E-5 * Math.Abs(exact));
        }

        /// <summary>
        /// Genz corner-peak family, c = (5, 5), against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzCornerPeak()
        {
            var asr2D = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzCornerPeak(x, y, 5, 5), 0, 1, 0, 1);
            asr2D.Integrate();
            double exact = Integrands.GenzCornerPeakExact(5, 5);
            Assert.AreEqual(exact, asr2D.Result, 1E-5 * Math.Abs(exact));
        }

        /// <summary>
        /// Genz Gaussian family, c = (10, 10), w = (0.4, 0.6), against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzGaussian()
        {
            var asr2D = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzGaussian(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            asr2D.Integrate();
            double exact = Integrands.GenzGaussianExact(10, 10, 0.4, 0.6);
            Assert.AreEqual(exact, asr2D.Result, 1E-5 * Math.Abs(exact));
        }

        /// <summary>
        /// Genz continuous (C0) family, c = (10, 10), w = (0.4, 0.6), against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzContinuous()
        {
            var asr2D = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzContinuous(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            asr2D.Integrate();
            double exact = Integrands.GenzContinuousExact(10, 10, 0.4, 0.6);
            Assert.AreEqual(exact, asr2D.Result, 1E-4 * Math.Abs(exact));
        }

        /// <summary>
        /// Genz discontinuous family, c = (5, 5), w = (0.4, 0.6), against its closed form. The jump
        /// lines are the hardest case for a fixed low-order rule; the tolerance pins the measured
        /// accuracy of the existing algorithm.
        /// </summary>
        [TestMethod]
        public void Test_GenzDiscontinuous()
        {
            var asr2D = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzDiscontinuous(x, y, 5, 5, 0.4, 0.6), 0, 1, 0, 1);
            asr2D.Integrate();
            double exact = Integrands.GenzDiscontinuousExact(5, 5, 0.4, 0.6);
            Assert.AreEqual(exact, asr2D.Result, 1E-2 * Math.Abs(exact));
        }

        /// <summary>
        /// The reported standard error must bound the actual error at convergence on smooth
        /// integrands; the honesty factor k pins the measured behavior of the existing estimator
        /// (the root-sum-square of the accepted Richardson corrections).
        /// </summary>
        [TestMethod]
        public void Test_ErrorEstimateHonesty()
        {
            var gauss = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzGaussian(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            gauss.Integrate();
            double exactGauss = Integrands.GenzGaussianExact(10, 10, 0.4, 0.6);
            Assert.IsLessThanOrEqualTo(Math.Max(gauss.StandardError, 1E-13 * Math.Abs(exactGauss)), Math.Abs(gauss.Result - exactGauss));

            var peak = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzProductPeak(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            peak.Integrate();
            double exactPeak = Integrands.GenzProductPeakExact(10, 10, 0.4, 0.6);
            Assert.IsLessThanOrEqualTo(Math.Max(peak.StandardError, 1E-13 * Math.Abs(exactPeak)), Math.Abs(peak.Result - exactPeak));
        }

        /// <summary>
        /// MinDepth forces refinement before acceptance: a trivially smooth integrand uses more
        /// function evaluations under a raised MinDepth.
        /// </summary>
        [TestMethod]
        public void Test_MinDepth_ForcesRefinement()
        {
            Func<double, double, double> func = (x, y) => x + y;
            var baseline = new AdaptiveSimpsonsRule2D(func, 0, 1, 0, 1);
            baseline.Integrate();
            var forced = new AdaptiveSimpsonsRule2D(func, 0, 1, 0, 1) { MinDepth = 3 };
            forced.Integrate();
            Assert.IsGreaterThan(baseline.FunctionEvaluations, forced.FunctionEvaluations);
            Assert.AreEqual(1d, forced.Result, 1E-10);
        }

        /// <summary>
        /// Exhausting the evaluation budget is reported through the status: a sharp peak with a
        /// tiny budget stops with MaximumFunctionEvaluationsReached and a finite estimate.
        /// </summary>
        [TestMethod]
        public void Test_MaxFunctionEvaluations_StatusReported()
        {
            var asr2D = new AdaptiveSimpsonsRule2D((x, y) => Integrands.GenzProductPeak(x, y, 50, 50, 0.4, 0.6), 0, 1, 0, 1)
            {
                MaxFunctionEvaluations = 500
            };
            asr2D.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumFunctionEvaluationsReached, asr2D.Status);
            Assert.IsTrue(Tools.IsFinite(asr2D.Result));
        }

        /// <summary>
        /// Constructor guards, pinning the existing exception types: a null function throws
        /// ArgumentNullException, and reversed bounds throw the ArgumentNullException the current
        /// constructor uses (a historical quirk, pinned as-is).
        /// </summary>
        [TestMethod]
        public void Test_Constructor_Guards()
        {
            Func<double, double, double> func = (x, y) => x + y;
            Assert.Throws<ArgumentNullException>(() => new AdaptiveSimpsonsRule2D(null!, 0, 1, 0, 1));
            Assert.Throws<ArgumentNullException>(() => new AdaptiveSimpsonsRule2D(func, 1, 0, 0, 1));
            Assert.Throws<ArgumentNullException>(() => new AdaptiveSimpsonsRule2D(func, 1, 1, 0, 1));
            Assert.Throws<ArgumentNullException>(() => new AdaptiveSimpsonsRule2D(func, 0, 1, 1, 0));
            Assert.Throws<ArgumentNullException>(() => new AdaptiveSimpsonsRule2D(func, 0, 1, 1, 1));
        }

        /// <summary>
        /// A machine-epsilon-width domain completes successfully with an essentially zero integral.
        /// </summary>
        [TestMethod]
        public void Test_DegenerateWidth()
        {
            var thin = new AdaptiveSimpsonsRule2D((x, y) => 1d, 0, 1E-16, 0, 1);
            thin.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, thin.Status);
            Assert.AreEqual(1E-16, thin.Result, 1E-17);
        }

    }
}
