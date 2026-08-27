using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the two-dimensional globally adaptive Gauss-Kronrod integration method.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_AdaptiveGaussKronrod2D
    {

        /// <summary>
        /// Test the 2D Adaptive Gauss-Kronrod algorithm with the Pi function. The disc indicator is
        /// discontinuous along the circle, so the run is expected to stop on the evaluation budget
        /// rather than the tolerance; the estimate must still be accurate.
        /// </summary>
        [TestMethod]
        public void Test_PI()
        {
            var agk2D = new AdaptiveGaussKronrod2D(Integrands.PI2D, -1, 1, -1, 1);
            agk2D.Integrate();
            var result = agk2D.Result;
            double trueResult = 3.14;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the 2D Adaptive Gauss-Kronrod algorithm with the sum of two normal distributions.
        /// </summary>
        [TestMethod]
        public void Test_SumOfTwoNormals()
        {
            var agk2D = new AdaptiveGaussKronrod2D(Integrands.SumOfNormals2D, 1E-15, 1 - 1E-15, 1E-15, 1 - 1E-15);
            agk2D.Integrate();
            var result = agk2D.Result;
            double trueResult = 40;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the 2D Adaptive Gauss-Kronrod algorithm with the X + Y integrand.
        /// </summary>
        [TestMethod]
        public void Test_XPlusY()
        {
            Func<double, double, double> func = (x, y) => x + y;
            var agk2D = new AdaptiveGaussKronrod2D(func, -1, 1, -1, 1);
            agk2D.Integrate();
            var result = agk2D.Result;
            double trueResult = 0;
            Assert.AreEqual(trueResult, result, 1E-5);
        }

        /// <summary>
        /// Test the 2D Adaptive Gauss-Kronrod algorithm with the X^2 + Y^2 integrand.
        /// </summary>
        [TestMethod]
        public void Test_XSquaredPlusYSquared()
        {
            Func<double, double, double> func = (x, y) => x * x + y * y;
            var agk2D = new AdaptiveGaussKronrod2D(func, -1, 1, -1, 1);
            agk2D.Integrate();
            var result = agk2D.Result;
            double trueResult = 2.6666667;
            Assert.AreEqual(trueResult, result, 1E-6);
        }

        /// <summary>
        /// Polynomials below the rule degree integrate exactly from the root evaluation alone.
        /// The G10K21 tensor product is exact for x^3*y^3 + x^5 + y^5, whose integral over
        /// [0,2] x [0,1] is 4*(1/4) + (32/3)*1 + 2*(1/6) = 12.
        /// </summary>
        [TestMethod]
        public void Test_PolynomialExactness()
        {
            Func<double, double, double> func = (x, y) => x * x * x * y * y * y + Math.Pow(x, 5) + Math.Pow(y, 5);
            var agk2D = new AdaptiveGaussKronrod2D(func, 0, 2, 0, 1);
            agk2D.Integrate();
            Assert.AreEqual(12d, agk2D.Result, 1E-11);
            Assert.AreEqual(IntegrationStatus.Success, agk2D.Status);
            // Both the Kronrod and embedded Gauss products are exact, so the root converges alone.
            Assert.AreEqual(441, agk2D.FunctionEvaluations);
        }

        /// <summary>
        /// A separable product f(x)*g(y) must match the product of two one-dimensional adaptive
        /// Gauss-Kronrod passes, and the closed form (e^1.5 - 1)*sin(1).
        /// </summary>
        [TestMethod]
        public void Test_SeparableProduct_MatchesTwo1DPasses()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Math.Exp(x) * Math.Cos(y), 0, 1.5, 0, 1);
            agk2D.Integrate();

            var agkX = new AdaptiveGaussKronrod(Math.Exp, 0, 1.5);
            agkX.Integrate();
            var agkY = new AdaptiveGaussKronrod(Math.Cos, 0, 1);
            agkY.Integrate();

            double product = agkX.Result * agkY.Result;
            double exact = (Math.Exp(1.5) - 1d) * Math.Sin(1d);
            Assert.AreEqual(product, agk2D.Result, 5E-8 * Math.Abs(product));
            Assert.AreEqual(exact, agk2D.Result, 1E-7 * Math.Abs(exact));
        }

        /// <summary>
        /// Genz oscillatory family, c = (5, 5), w1 = 0.3, against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzOscillatory()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzOscillatory(x, y, 5, 5, 0.3), 0, 1, 0, 1);
            agk2D.Integrate();
            double exact = Integrands.GenzOscillatoryExact(5, 5, 0.3);
            Assert.AreEqual(exact, agk2D.Result, 1E-6 * Math.Abs(exact));
            Assert.AreEqual(IntegrationStatus.Success, agk2D.Status);
        }

        /// <summary>
        /// Genz product-peak family, c = (10, 10), w = (0.4, 0.6), against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzProductPeak()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzProductPeak(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            agk2D.Integrate();
            double exact = Integrands.GenzProductPeakExact(10, 10, 0.4, 0.6);
            Assert.AreEqual(exact, agk2D.Result, 1E-6 * Math.Abs(exact));
            Assert.AreEqual(IntegrationStatus.Success, agk2D.Status);
        }

        /// <summary>
        /// Genz corner-peak family, c = (5, 5), against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzCornerPeak()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzCornerPeak(x, y, 5, 5), 0, 1, 0, 1);
            agk2D.Integrate();
            double exact = Integrands.GenzCornerPeakExact(5, 5);
            Assert.AreEqual(exact, agk2D.Result, 1E-6 * Math.Abs(exact));
            Assert.AreEqual(IntegrationStatus.Success, agk2D.Status);
        }

        /// <summary>
        /// Genz Gaussian family, c = (10, 10), w = (0.4, 0.6), against its closed form.
        /// </summary>
        [TestMethod]
        public void Test_GenzGaussian()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzGaussian(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            agk2D.Integrate();
            double exact = Integrands.GenzGaussianExact(10, 10, 0.4, 0.6);
            Assert.AreEqual(exact, agk2D.Result, 1E-6 * Math.Abs(exact));
            Assert.AreEqual(IntegrationStatus.Success, agk2D.Status);
        }

        /// <summary>
        /// Genz continuous (C0) family, c = (10, 10), w = (0.4, 0.6), against its closed form.
        /// The gradient discontinuities along the two ridge lines slow convergence relative to the
        /// smooth families but the tolerance is still met.
        /// </summary>
        [TestMethod]
        public void Test_GenzContinuous()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzContinuous(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            agk2D.Integrate();
            double exact = Integrands.GenzContinuousExact(10, 10, 0.4, 0.6);
            Assert.AreEqual(exact, agk2D.Result, 1E-5 * Math.Abs(exact));
        }

        /// <summary>
        /// Genz discontinuous family, c = (5, 5), w = (0.4, 0.6), against its closed form. As with
        /// the disc indicator, the jump lines defeat tolerance-level convergence, so the run stops
        /// on the budget with an accurate estimate rather than reporting silent success.
        /// </summary>
        [TestMethod]
        public void Test_GenzDiscontinuous()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzDiscontinuous(x, y, 5, 5, 0.4, 0.6), 0, 1, 0, 1);
            agk2D.Integrate();
            double exact = Integrands.GenzDiscontinuousExact(5, 5, 0.4, 0.6);
            Assert.AreEqual(exact, agk2D.Result, 2E-3 * Math.Abs(exact));
        }

        /// <summary>
        /// The reported standard error must bound the actual error at convergence. StandardError is
        /// the sum of per-region |Kronrod - Gauss| deviations - a conservative bound - so the honesty
        /// factor k = 1 with a small floating-point floor.
        /// </summary>
        [TestMethod]
        public void Test_ErrorEstimateHonesty()
        {
            var gauss = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzGaussian(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            gauss.Integrate();
            double exactGauss = Integrands.GenzGaussianExact(10, 10, 0.4, 0.6);
            Assert.IsLessThanOrEqualTo(Math.Max(gauss.StandardError, 1E-13 * Math.Abs(exactGauss)), Math.Abs(gauss.Result - exactGauss));

            var peak = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzProductPeak(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1);
            peak.Integrate();
            double exactPeak = Integrands.GenzProductPeakExact(10, 10, 0.4, 0.6);
            Assert.IsLessThanOrEqualTo(Math.Max(peak.StandardError, 1E-13 * Math.Abs(exactPeak)), Math.Abs(peak.Result - exactPeak));
        }

        /// <summary>
        /// MinDepth forces refinement before convergence may be declared: a trivially smooth
        /// integrand converges at the root with MinDepth = 0 but must be pre-split at MinDepth = 2.
        /// </summary>
        [TestMethod]
        public void Test_MinDepth_ForcesRefinement()
        {
            Func<double, double, double> func = (x, y) => x + y;
            var baseline = new AdaptiveGaussKronrod2D(func, 0, 1, 0, 1);
            baseline.Integrate();
            Assert.AreEqual(441, baseline.FunctionEvaluations);

            var forced = new AdaptiveGaussKronrod2D(func, 0, 1, 0, 1) { MinDepth = 2 };
            forced.Integrate();
            // Root + two depth-1 children + four depth-2 children = 7 evaluated regions.
            Assert.AreEqual(441 * 7, forced.FunctionEvaluations);
            // The integral of x + y over the unit square is 1, from either evaluation count.
            Assert.AreEqual(1d, baseline.Result, 1E-10);
            Assert.AreEqual(1d, forced.Result, 1E-10);
        }

        /// <summary>
        /// Exhausting the evaluation budget is reported through the status rather than silently:
        /// a sharp peak with a tiny budget stops with MaximumFunctionEvaluationsReached and a finite
        /// estimate.
        /// </summary>
        [TestMethod]
        public void Test_MaxFunctionEvaluations_StatusReported()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzProductPeak(x, y, 50, 50, 0.4, 0.6), 0, 1, 0, 1)
            {
                MaxFunctionEvaluations = 2000
            };
            agk2D.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumFunctionEvaluationsReached, agk2D.Status);
            Assert.IsTrue(Tools.IsFinite(agk2D.Result));
        }

        /// <summary>
        /// Exhausting MaxDepth without meeting the tolerance is reported as
        /// MaximumIterationsReached rather than silent success: with MaxDepth = 0 the unconverged
        /// root cannot be refined.
        /// </summary>
        [TestMethod]
        public void Test_MaxDepthExhaustion_StatusReported()
        {
            var agk2D = new AdaptiveGaussKronrod2D((x, y) => Integrands.GenzGaussian(x, y, 10, 10, 0.4, 0.6), 0, 1, 0, 1)
            {
                MaxDepth = 0
            };
            agk2D.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, agk2D.Status);
            Assert.IsTrue(Tools.IsFinite(agk2D.Result));
            Assert.AreEqual(441, agk2D.FunctionEvaluations);
        }

        /// <summary>
        /// A throwing integrand surfaces as Failure without throwing when ReportFailure is false,
        /// and rethrows when ReportFailure is true.
        /// </summary>
        [TestMethod]
        public void Test_ReportFailure()
        {
            Func<double, double, double> throwing = (x, y) => throw new InvalidOperationException("Integrand failure.");
            var silent = new AdaptiveGaussKronrod2D(throwing, 0, 1, 0, 1) { ReportFailure = false };
            silent.Integrate();
            Assert.AreEqual(IntegrationStatus.Failure, silent.Status);

            var loud = new AdaptiveGaussKronrod2D(throwing, 0, 1, 0, 1) { ReportFailure = true };
            Assert.Throws<InvalidOperationException>(() => loud.Integrate());
        }

        /// <summary>
        /// Constructor guards: a null function and invalid or non-finite bounds throw.
        /// </summary>
        [TestMethod]
        public void Test_Constructor_Guards()
        {
            Func<double, double, double> func = (x, y) => x + y;
            Assert.Throws<ArgumentNullException>(() => new AdaptiveGaussKronrod2D(null!, 0, 1, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod2D(func, 1, 1, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod2D(func, 1, 0, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod2D(func, 0, 1, 1, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod2D(func, 0, 1, 1, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod2D(func, double.NaN, 1, 0, 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod2D(func, 0, double.PositiveInfinity, 0, 1));
        }

        /// <summary>
        /// Depth-setting guards: a negative MinDepth or MaxDepth below MinDepth throws at
        /// integration time.
        /// </summary>
        [TestMethod]
        public void Test_DepthGuards()
        {
            Func<double, double, double> func = (x, y) => x + y;
            var negative = new AdaptiveGaussKronrod2D(func, 0, 1, 0, 1) { MinDepth = -1 };
            Assert.Throws<ArgumentOutOfRangeException>(() => negative.Integrate());
            var inverted = new AdaptiveGaussKronrod2D(func, 0, 1, 0, 1) { MinDepth = 3, MaxDepth = 2 };
            Assert.Throws<ArgumentOutOfRangeException>(() => inverted.Integrate());
        }

        /// <summary>
        /// A machine-epsilon-width domain completes successfully with an essentially zero integral,
        /// including under a forced minimum depth where only the wide axis can split.
        /// </summary>
        [TestMethod]
        public void Test_DegenerateWidth()
        {
            var thin = new AdaptiveGaussKronrod2D((x, y) => 1d, 0, 1E-16, 0, 1);
            thin.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, thin.Status);
            Assert.AreEqual(1E-16, thin.Result, 1E-17);

            var forced = new AdaptiveGaussKronrod2D((x, y) => 1d, 0, 1E-16, 0, 1) { MinDepth = 2 };
            forced.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, forced.Status);
            Assert.AreEqual(1E-16, forced.Result, 1E-17);
        }

        /// <summary>
        /// Far from the origin an axis can be wider than the absolute machine-epsilon floor while
        /// its width sits at one unit in the last place, where the midpoint rounds onto an endpoint
        /// and bisection reproduces the region. Every node of such a region rounds onto a single
        /// representable abscissa, so its error estimate is exactly zero and error-driven refinement
        /// never asks to split it; the forced minimum-depth pass splits it regardless, and unguarded
        /// it re-split a bit-identical child once per depth level, burning 882 evaluations per
        /// wasted split. The region must freeze instead, and a splittable domain at the same offset
        /// must still converge.
        /// </summary>
        [TestMethod]
        public void Test_LargeOffsetDomain()
        {
            // The representable spacing at 1e16 is 2, so this x-domain is one unit in the last place
            // wide and cannot be bisected. Unguarded, MinDepth = 30 cost about 27,000 evaluations.
            var ulpWide = new AdaptiveGaussKronrod2D((x, y) => 1d, 1E16, 1E16 + 2d, 0, 1) { MinDepth = 30 };
            ulpWide.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, ulpWide.Status);
            Assert.AreEqual(2d, ulpWide.Result, 1E-8);
            Assert.IsLessThan(2000, ulpWide.FunctionEvaluations);

            // A two-ulp domain bisects once into representable halves and converges to the exact
            // area at the default settings.
            var flat = new AdaptiveGaussKronrod2D((x, y) => 1d, 1E16, 1E16 + 4d, 0, 1);
            flat.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, flat.Status);
            Assert.AreEqual(4d, flat.Result, 1E-8);
        }

        /// <summary>
        /// The recorder reports the final composite rule: weights sum to the domain area, weighted
        /// function values reproduce the result, the node count is a whole number of regions, and
        /// attaching a recorder does not change the computed result.
        /// </summary>
        [TestMethod]
        public void Test_Recorder_MassAndResultIdentities()
        {
            Func<double, double, double> func = (x, y) => Integrands.GenzGaussian(x, y, 10, 10, 0.4, 0.6);
            var plain = new AdaptiveGaussKronrod2D(func, 0, 1, 0, 1);
            plain.Integrate();

            double sumW = 0, sumWF = 0;
            int count = 0;
            var recorded = new AdaptiveGaussKronrod2D(func, 0, 1, 0, 1)
            {
                Recorder = (x, y, w, f) => { sumW += w; sumWF += w * f; count++; }
            };
            recorded.Integrate();

            Assert.AreEqual(BitConverter.DoubleToInt64Bits(plain.Result), BitConverter.DoubleToInt64Bits(recorded.Result));
            Assert.AreEqual(1d, sumW, 1E-12);
            Assert.AreEqual(recorded.Result, sumWF, 1E-12 * Math.Abs(recorded.Result));
            Assert.AreEqual(0, count % 441);
            Assert.IsGreaterThanOrEqualTo(441, count);
        }

        /// <summary>
        /// Identical inputs produce bit-identical results and error bounds on repeated runs.
        /// </summary>
        [TestMethod]
        public void Test_Determinism_RepeatedRunsBitEqual()
        {
            Func<double, double, double> func = (x, y) => Integrands.GenzProductPeak(x, y, 10, 10, 0.4, 0.6);
            var first = new AdaptiveGaussKronrod2D(func, 0, 1, 0, 1);
            first.Integrate();
            var second = new AdaptiveGaussKronrod2D(func, 0, 1, 0, 1);
            second.Integrate();
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(first.Result), BitConverter.DoubleToInt64Bits(second.Result));
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(first.StandardError), BitConverter.DoubleToInt64Bits(second.StandardError));
            Assert.AreEqual(first.FunctionEvaluations, second.FunctionEvaluations);
        }

    }
}
