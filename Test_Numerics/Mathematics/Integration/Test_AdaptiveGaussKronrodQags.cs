using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using System;
using System.Collections.Generic;
using Numerics.Sampling;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the QUADPACK QAGS method of the Adaptive Gauss-Kronrod integrator. The evaluation counts
    /// are those of scipy.integrate.quad 1.18.1 (QUADPACK DQAGSE); the exact values are closed forms.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_AdaptiveGaussKronrodQags
    {
        /// <summary>
        /// Test the Adaptive Gauss-Kronrod method with a one dimensional function
        /// </summary>
        [TestMethod]
        public void Test_FX3()
        {
            var adapt = new AdaptiveGaussKronrod(Integrands.FX3, 0d, 1d) { Method = AdaptiveGaussKronrodMethod.Qags };
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
            var adapt = new AdaptiveGaussKronrod(Integrands.Cosine, -1, 1) { Method = AdaptiveGaussKronrodMethod.Qags };
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
            var adapt = new AdaptiveGaussKronrod(Integrands.Sine, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags };
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
            var adapt = new AdaptiveGaussKronrod(Integrands.FXX, 0, 2) { Method = AdaptiveGaussKronrodMethod.Qags };
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
            var adapt = new AdaptiveGaussKronrod(Integrands.FXXX, 0, 2) { Method = AdaptiveGaussKronrodMethod.Qags };
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
        /// A smooth integrand reports Success, and its result has the bit pattern of scipy.integrate.quad
        /// 1.18.1, whose DQK21 accumulates the sums in the same order. The analytic value is
        /// (1 + exp(-2) (3 sin 6 - cos 6)) / 10.
        /// </summary>
        [TestMethod]
        public void Test_SmoothIntegrand_SuccessAndBitsUnchanged()
        {
            var adapt = new AdaptiveGaussKronrod(x => Math.Exp(-x) * Math.Cos(3 * x), 0, 2) { Method = AdaptiveGaussKronrodMethod.Qags };
            adapt.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
            Assert.AreEqual(0x3FB35E863514D48E, BitConverter.DoubleToInt64Bits(adapt.Result));
            Assert.AreEqual(21, adapt.FunctionEvaluations);
            Assert.AreEqual((1 + Math.Exp(-2) * (3 * Math.Sin(6) - Math.Cos(6))) / 10, adapt.Result, 1E-15);
        }

        /// <summary>
        /// Intervals retired at the depth limit without meeting the tolerance report
        /// MaximumIterationsReached; the result value is unchanged by the limit status.
        /// </summary>
        [TestMethod]
        public void Test_ForcedAcceptance_ReportsMaximumIterations()
        {
            // The integral of 1/x on [0, 1] diverges; refinement stops when [0, 2^-100] reaches the depth
            // limit. The error estimate 9.35056 matches scipy.integrate.quad 1.18.1, which reports the same
            // abserr before its own width limit.
            var divergent = new AdaptiveGaussKronrod(x => 1.0 / x, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags };
            divergent.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, divergent.Status);
            Assert.AreEqual(0x405341C00D250C5A, BitConverter.DoubleToInt64Bits(divergent.Result));
            Assert.AreEqual(100, divergent.Iterations);
            Assert.AreEqual(4221, divergent.FunctionEvaluations);
            Assert.AreEqual(9.35056, divergent.StandardError, 1E-5);

            // A unit step at 0.3 integrates to 0.7, but three levels cannot resolve the jump.
            var step = new AdaptiveGaussKronrod(x => x < 0.3 ? 0.0 : 1.0, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, MaxDepth = 3 };
            step.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumIterationsReached, step.Status);
            Assert.AreEqual(0x3FE671E68616C21D, BitConverter.DoubleToInt64Bits(step.Result));
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
        /// A NaN integrand stops with Failure when a NaN interval reaches the depth limit (bisected depth
        /// first, 21 + 42 * 100 evaluations) instead of spending the 10,000,000-evaluation budget, and
        /// throws when ReportFailure is true.
        /// </summary>
        [TestMethod]
        public void Test_NaNIntegrand_FailsFast()
        {
            var silent = new AdaptiveGaussKronrod(x => double.NaN, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, ReportFailure = false };
            silent.Integrate();
            Assert.AreEqual(IntegrationStatus.Failure, silent.Status);
            Assert.IsTrue(double.IsNaN(silent.Result));
            Assert.AreEqual(4221, silent.FunctionEvaluations);

            var loud = new AdaptiveGaussKronrod(x => double.NaN, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags };
            Assert.Throws<ArithmeticException>(() => loud.Integrate());
            Assert.AreEqual(IntegrationStatus.Failure, loud.Status);
        }

        /// <summary>
        /// An infinite value at the center node is bisected away: 1/sqrt|x| on [-1, 1] integrates to 4.
        /// With the singular point removed, the run is the one scipy.integrate.quad 1.18.1 makes for the
        /// integrand with f(0) = 0: 483 evaluations, result 3.9999999999999813 and error estimate
        /// 5.684341886080802E-14, bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_InfiniteCenterNode_IsBisected()
        {
            var infinite = new AdaptiveGaussKronrod(x => 1.0 / Math.Sqrt(Math.Abs(x)), -1, 1) { Method = AdaptiveGaussKronrodMethod.Qags };
            infinite.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, infinite.Status);
            Assert.AreEqual(3.9999999999999813, infinite.Result, 0d);
            Assert.AreEqual(5.684341886080802E-14, infinite.StandardError, 0d);
            Assert.AreEqual(4d, infinite.Result, infinite.StandardError);
            Assert.AreEqual(483, infinite.FunctionEvaluations);
        }

        /// <summary>
        /// Endpoint singularities are resolved by the epsilon-algorithm extrapolation. Exact values; the
        /// evaluation counts are those of scipy.integrate.quad 1.18.1 (QUADPACK DQAGSE) at
        /// epsabs = epsrel = 1E-8.
        /// </summary>
        [TestMethod]
        public void Test_EndpointSingularities_MatchQuadpack()
        {
            var cases = new (Func<double, double> F, double Exact, int Evaluations, double Tolerance)[]
            {
                (x => 1 / Math.Sqrt(x), 2d, 231, 1E-14),
                (Math.Log, -1d, 231, 1E-14),
                (x => Math.Pow(x, -0.9), 10d, 231, 1E-12),
                (x => Math.Log(x) / Math.Sqrt(x), -4d, 315, 1E-12),
            };
            foreach (var c in cases)
            {
                var adapt = new AdaptiveGaussKronrod(c.F, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags };
                adapt.Integrate();
                Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
                Assert.AreEqual(c.Exact, adapt.Result, c.Tolerance);
                Assert.AreEqual(c.Evaluations, adapt.FunctionEvaluations);
                Assert.AreEqual((c.Evaluations / 21 - 1) / 2, adapt.Iterations);
                Assert.IsLessThanOrEqualTo(1E-8 * Math.Abs(c.Exact), adapt.StandardError);
            }
        }

        /// <summary>
        /// Peaked and oscillatory integrands match the QUADPACK evaluation counts of scipy.integrate.quad
        /// 1.18.1. Exact values: 100 (atan 70 + atan 30), sin(100) / 100, -2 pi / 50 and 5/18.
        /// </summary>
        [TestMethod]
        public void Test_PeakAndOscillatory_MatchQuadpack()
        {
            var cases = new (Func<double, double> F, double B, double Exact, int Evaluations)[]
            {
                (x => 1 / ((x - 0.3) * (x - 0.3) + 1E-4), 1, 100 * (Math.Atan(70) + Math.Atan(30)), 315),
                (x => Math.Cos(100 * x), 1, Math.Sin(100) / 100, 315),
                (x => x * Math.Sin(50 * x), 2 * Math.PI, -2 * Math.PI / 50, 1323),
                (x => Math.Abs(x - 1.0 / 3), 1, 5.0 / 18, 189),
            };
            foreach (var c in cases)
            {
                var adapt = new AdaptiveGaussKronrod(c.F, 0, c.B) { Method = AdaptiveGaussKronrodMethod.Qags };
                adapt.Integrate();
                Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
                Assert.AreEqual(c.Exact, adapt.Result, 1E-8 * Math.Max(1, Math.Abs(c.Exact)));
                Assert.AreEqual(c.Evaluations, adapt.FunctionEvaluations);
            }
        }

        /// <summary>
        /// The QUADPACK tolerance rule: an absolute-only request (relative tolerance 0) is legal; with a zero
        /// absolute tolerance the relative tolerance must be at least max(50 eps, 5E-29) = 1.11E-14.
        /// </summary>
        [TestMethod]
        public void Test_ToleranceRule()
        {
            var absOnly = new AdaptiveGaussKronrod(x => 1 / Math.Sqrt(x), 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, RelativeTolerance = 0, AbsoluteTolerance = 1E-12 };
            absOnly.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, absOnly.Status);
            Assert.AreEqual(2d, absOnly.Result, 1E-12);

            var relOnly = new AdaptiveGaussKronrod(Math.Exp, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, RelativeTolerance = 1.2E-14, AbsoluteTolerance = 0 };
            relOnly.Integrate();
            Assert.AreEqual(Math.E - 1, relOnly.Result, 1E-15);

            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(Math.Exp, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, RelativeTolerance = 1E-14, AbsoluteTolerance = 0 }.Integrate());
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(Math.Exp, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, RelativeTolerance = double.NaN }.Integrate());
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(Math.Exp, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, AbsoluteTolerance = -1E-10 }.Integrate());
            Assert.Throws<ArgumentOutOfRangeException>(() => new AdaptiveGaussKronrod(Math.Exp, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, RelativeTolerance = 2 }.Integrate());
        }

        /// <summary>
        /// The evaluation budget sets the subinterval limit floor((floor(budget / 21) + 1) / 2): a budget of
        /// 100 allows two intervals, so the run stops after one bisection (63 evaluations) and reports
        /// MaximumFunctionEvaluationsReached.
        /// </summary>
        [TestMethod]
        public void Test_EvaluationBudget_IsNeverExceeded()
        {
            var adapt = new AdaptiveGaussKronrod(x => 1 / Math.Sqrt(x), 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, MaxFunctionEvaluations = 100 };
            adapt.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumFunctionEvaluationsReached, adapt.Status);
            Assert.AreEqual(63, adapt.FunctionEvaluations);
            Assert.AreEqual(1, adapt.Iterations);

            var peak = new AdaptiveGaussKronrod(x => 1 / ((x - 0.3) * (x - 0.3) + 1E-6), 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, MaxFunctionEvaluations = 1000 };
            peak.Integrate();
            Assert.IsLessThanOrEqualTo(1000, peak.FunctionEvaluations);
        }

        /// <summary>
        /// MinDepth starts from the uniform partition into 2^MinDepth intervals without evaluating their
        /// parents: x^2 on [0, 1] at MinDepth 2 takes 4 * 21 evaluations.
        /// </summary>
        [TestMethod]
        public void Test_MinDepth_StartsFromUniformPartition()
        {
            var adapt = new AdaptiveGaussKronrod(x => x * x, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, MinDepth = 2 };
            adapt.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
            Assert.AreEqual(84, adapt.FunctionEvaluations);
            Assert.AreEqual(1d / 3, adapt.Result, 1E-15);
        }

        /// <summary>
        /// A NaN on an interval that is subdivided further does not stop integration: sin(x - 0.5)/(x - 0.5)
        /// is 0/0 at the center node of [0, 1], yet the children are finite. The result keeps its bit
        /// pattern and matches 2 Si(0.5) from scipy.special.sici.
        /// </summary>
        [TestMethod]
        public void Test_RemovableNaN_ResultUnchanged()
        {
            var adapt = new AdaptiveGaussKronrod(x => Math.Sin(x - 0.5) / (x - 0.5), 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags };
            adapt.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
            Assert.AreEqual(0x3FEF8F126A7A3CFA, BitConverter.DoubleToInt64Bits(adapt.Result));
            Assert.AreEqual(63, adapt.FunctionEvaluations);
            Assert.AreEqual(0.9862148360861335, adapt.Result, 1E-15);
        }

        /// <summary>A sharp peak that forces deep adaptive subdivision on [0, 1].</summary>
        private static double SharpPeak(double x)
        {
            double t = (x - 0.3) * 60d;
            return Math.Exp(-t * t) + 0.1 * x;
        }

        /// <summary>
        /// When the extrapolated value is returned, the records still partition the domain and reproduce the
        /// composite Kronrod sum of the final partition, which is the less accurate estimate: for x^-0.5 on
        /// [0, 1] the result is 2 while the six-interval composite sum is 1.99426.
        /// </summary>
        [TestMethod]
        public void Test_Recorder_ExtrapolatedResult()
        {
            double weightSum = 0d, composite = 0d;
            int count = 0;
            var gk = new AdaptiveGaussKronrod(x => 1 / Math.Sqrt(x), 0d, 1d)
            {
                Method = AdaptiveGaussKronrodMethod.Qags,
                Recorder = (x, w, fx) => { weightSum += w; composite += w * fx; count++; },
            };
            gk.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, gk.Status);
            Assert.AreEqual(2d, gk.Result, 1E-14);
            Assert.AreEqual(6 * 21, count);
            Assert.AreEqual(1d, weightSum, 1E-12);
            Assert.AreEqual(1.994259303194289, composite, 1E-12);
        }

        /// <summary>
        /// A zero-width stratification bin is evaluated once at its single point with zero weight, as the
        /// TotalRisk ledger expects for a saturated hazard curve.
        /// </summary>
        [TestMethod]
        public void Test_Recorder_ZeroWidthBin()
        {
            var bins = new List<StratificationBin> { new StratificationBin(0d, 0.5d), new StratificationBin(0.5d, 0.5d), new StratificationBin(0.5d, 1d) };
            var records = new List<(double X, double Weight, double Value)>();
            var gk = new AdaptiveGaussKronrod(SharpPeak, 0d, 1d)
            {
                Method = AdaptiveGaussKronrodMethod.Qags,
                RelativeTolerance = 1E-10,
                MinDepth = 2,
                Recorder = (x, w, fx) => records.Add((x, w, fx)),
            };
            gk.Integrate(bins);
            Assert.AreEqual(IntegrationStatus.Success, gk.Status);
            int atHalf = 0;
            double weightSum = 0d;
            foreach (var r in records)
            {
                weightSum += r.Weight;
                if (r.X == 0.5 && r.Weight == 0d) atHalf++;
            }
            Assert.AreEqual(21, atHalf);
            Assert.AreEqual(1d, weightSum, 1E-12);
        }

        /// <summary>
        /// Repeated integration on one instance gives identical bits: the extrapolation history is per run. The
        /// interior singularity and the MinDepth partition both start the epsilon table from a partition.
        /// </summary>
        [TestMethod]
        public void Test_RepeatedIntegration_IsBitwiseIdentical()
        {
            var cases = new[]
            {
                new AdaptiveGaussKronrod(x => 1.0 / Math.Sqrt(Math.Abs(x)), -1, 1) { Method = AdaptiveGaussKronrodMethod.Qags },
                new AdaptiveGaussKronrod(x => 1 / Math.Sqrt(x), 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, MinDepth = 2 },
                new AdaptiveGaussKronrod(Math.Log, 0, 1) { Method = AdaptiveGaussKronrodMethod.Qags, MinDepth = 3 },
            };
            foreach (var adapt in cases)
            {
                adapt.Integrate();
                long bits = BitConverter.DoubleToInt64Bits(adapt.Result);
                int evaluations = adapt.FunctionEvaluations;
                double error = adapt.StandardError;
                for (int run = 0; run < 3; run++)
                {
                    adapt.Integrate();
                    Assert.AreEqual(bits, BitConverter.DoubleToInt64Bits(adapt.Result));
                    Assert.AreEqual(evaluations, adapt.FunctionEvaluations);
                    Assert.AreEqual(error, adapt.StandardError, 0d);
                }
                Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
            }
            Assert.AreEqual(4d, cases[0].Result, 1E-13);
            Assert.AreEqual(2d, cases[1].Result, 1E-14);
            Assert.AreEqual(-1d, cases[2].Result, 1E-14);
        }

        /// <summary>
        /// When the subinterval limit is reached while a non-finite interval is still being bisected, the run
        /// stops with MaximumFunctionEvaluationsReached instead of throwing, even with ReportFailure on.
        /// </summary>
        [TestMethod]
        public void Test_NonFiniteBudgetExhaustion_ReportsMaximumFunctionEvaluations()
        {
            var adapt = new AdaptiveGaussKronrod(x => x < 0.5 ? double.NaN : 1.0, 0, 1)
            {
                Method = AdaptiveGaussKronrodMethod.Qags,
                MaxFunctionEvaluations = 100,
            };
            adapt.Integrate();
            Assert.AreEqual(IntegrationStatus.MaximumFunctionEvaluationsReached, adapt.Status);
            Assert.AreEqual(63, adapt.FunctionEvaluations);
            Assert.IsTrue(double.IsNaN(adapt.Result));

            int records = 0;
            var recorded = new AdaptiveGaussKronrod(x => x < 0.5 ? double.NaN : 1.0, 0, 1)
            {
                Method = AdaptiveGaussKronrodMethod.Qags,
                MaxFunctionEvaluations = 100,
                Recorder = (x, w, fx) => records++,
            };
            recorded.Integrate();
            Assert.AreEqual(2 * 21, records);
        }

        /// <summary>
        /// Refinement forced by MinFunctionEvaluations past convergence is not counted as roundoff, so a smooth
        /// integrand still reports Success after at least the minimum number of evaluations.
        /// </summary>
        [TestMethod]
        public void Test_MinFunctionEvaluations_ReportsSuccess()
        {
            var adapt = new AdaptiveGaussKronrod(Math.Exp, 0, 1)
            {
                Method = AdaptiveGaussKronrodMethod.Qags,
                MinFunctionEvaluations = 1000,
            };
            adapt.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, adapt.Status);
            Assert.IsGreaterThanOrEqualTo(1000, adapt.FunctionEvaluations);
            Assert.AreEqual(Math.E - 1, adapt.Result, 1E-14);
        }
    }
}
