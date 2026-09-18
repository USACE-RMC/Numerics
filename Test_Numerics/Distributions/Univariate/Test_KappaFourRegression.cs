using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using Numerics;
using Numerics.Distributions;
using Numerics.Mathematics.Integration;

namespace Distributions.Univariate
{
    /// <summary>Independent regression contracts for the Kappa Four numerical repair.</summary>
    [TestClass]
    public class Test_KappaFourRegression
    {
        // Exact 38 observations supplied with the September 2026 fitting regression.
        private static readonly double[] FittingSample =
        [
            1.354784607887268, 0.41693252325057983, 0.8899999856948853, 0.8853314518928528,
            2.170344829559326, 1.334205150604248, 0.5150537490844727, 0.8244029879570007,
            0.5099999904632568, 0.44740423560142517, 0.739365816116333, 1.8200000524520874,
            1.9778419733047485, 0.6553186774253845, 0.8488552570343018, 2.2905187606811523,
            0.9865325689315796, 0.7076110243797302, 0.2929774820804596, 1.5241589546203613,
            0.5311949253082275, 0.5221586227416992, 0.907939612865448, 0.2859921157360077,
            0.6391091346740723, 0.6330636739730835, 0.520042359828949, 2.0499587059020996,
            1.83543860912323, 2.450000047683716, 0.6613662838935852, 1.1201441287994385,
            1.020573377609253, 0.5620101094245911, 0.5419405102729797, 2.2691876888275146,
            1.5167182683944702, 3.119999885559082
        ];

        /// <summary>Frozen 90-digit formulas exercise shape signs, zeros, near zeros, and both tails.</summary>
        [TestMethod]
        public void ProbabilityFunctionsMatchIndependentDecimalOracle()
        {
            var assembly = typeof(Test_KappaFourRegression).Assembly;
            using var reader = new StreamReader(assembly.GetManifestResourceStream(assembly.GetManifestResourceNames().Single(n => n.EndsWith("kappa-four-probabilities.csv"))));
            reader.ReadLine();
            var failures = new List<string>();
            int count = 0;
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                double[] r = line.Split(',').Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
                var d = new KappaFour(0, 1, r[0], r[1]);
                string context = $"k={r[0]:R}, h={r[1]:R}, p={r[2]:R}";
                double actual = d.InverseCDF(r[2]);
                if (!Tools.IsFinite(actual) || Math.Abs(actual - r[4]) > 2E-12 * Math.Max(1, Math.Abs(r[4])))
                    failures.Add($"quantile {context}: {actual:R} expected {r[4]:R}");
                double[] gradient = d.QuantileGradient(r[2]);
                for (int component = 2; component < 4; component++)
                    if (!Tools.IsFinite(gradient[component]) || Math.Abs(gradient[component] - r[component + 5]) > 2E-8 * Math.Max(1E-20, Math.Abs(r[component + 5])))
                        failures.Add($"gradient {component} {context}: {gradient[component]:R} expected {r[component + 5]:R}");
                // At a represented support endpoint, the public API returns the declared one-sided
                // limit. Decimal may place the same rounded x just inside the exact real support.
                if (r[3] > d.Minimum && r[3] < d.Maximum)
                {
                    actual = d.CDF(r[3]);
                    if (!Tools.IsFinite(actual) || Math.Abs(actual - r[5]) > 2E-10)
                        failures.Add($"CDF {context}: {actual:R} expected {r[5]:R}");
                    actual = d.PDF(r[3]);
                    if (!Tools.IsFinite(actual) || Math.Abs(actual - r[6]) > 2E-8 * Math.Max(1E-300, r[6]))
                        failures.Add($"PDF {context}: {actual:R} expected {r[6]:R}");
                }
                count++;
            }
            Assert.AreEqual(1546, count);
            Assert.IsEmpty(failures, string.Join(Environment.NewLine, failures.Take(12)));
        }

        /// <summary>The reported sample receives a supported, finite-likelihood MLE initializer.</summary>
        [TestMethod]
        public void FittingInitializerSupportsEveryObservationWithinOriginalBounds()
        {
            var d = new KappaFour();
            var constraints = d.GetParameterConstraints(FittingSample);
            AssertVector(new[] { -10d, Tools.DoubleMachineEpsilon, -10d, -2d }, constraints.Item2, 0);
            AssertVector(new[] { 10d, 100d, 10d, 2d }, constraints.Item3, 0);
            d.SetParameters(constraints.Item1);
            Assert.IsTrue(Tools.IsFinite(d.LogLikelihood(FittingSample)));
            foreach (double observation in FittingSample)
                Assert.IsTrue(observation > d.Minimum && observation < d.Maximum);
        }

        /// <summary>MLE cannot install a failed or nonfinite optimizer result.</summary>
        [TestMethod]
        public void MleDoesNotInstallAnUnsuccessfulEstimate()
        {
            var d = new KappaFour();
            double[] before = d.GetParameters;
            // Nonfinite data produce a deterministic fitting failure independent of the optimizer trajectory.
            double[] invalidSample = (double[])FittingSample.Clone();
            invalidSample[17] = double.NaN;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => d.Estimate(invalidSample, ParameterEstimationMethod.MaximumLikelihood));
            AssertVector(before, d.GetParameters, 0);
        }

        /// <summary>A failed seeded bootstrap propagates failure instead of returning invalid parameters.</summary>
        [TestMethod]
        public void BootstrapRejectsAnInsufficientFittingSample()
        {
            var d = new KappaFour(0, 1, 0, 0);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => d.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 1, 12345));
            AssertVector(new[] { 0d, 1d, 0d, 0d }, d.GetParameters, 0);
        }

        /// <summary>Seeded resampling reaches estimation and propagates failure for values that round to a constant.</summary>
        [TestMethod]
        public void BootstrapPropagatesEstimationFailure()
        {
            var d = new KappaFour(1, double.Epsilon, 0, 0);
            AssertVector(new[] { 1d, 1d, 1d, 1d }, d.GenerateRandomValues(4, 12345), 0d);
            Assert.ThrowsExactly<InvalidOperationException>(() => d.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 4, 12345));
            AssertVector(new[] { 1d, double.Epsilon, 0d, 0d }, d.GetParameters, 0d);
        }

        /// <summary>The endpoint likelihood singularity is preserved, rather than hidden by a fitting penalty.</summary>
        [TestMethod]
        public void UnrestrictedLikelihoodIncreasesTowardTheSampleMinimum()
        {
            double minimum = 0.2859921157360077;
            double previous = double.NegativeInfinity;
            foreach (double epsilon in new[] { 1E-2, 1E-4, 1E-6, 1E-8, 1E-10 })
            {
                var d = new KappaFour(minimum - Math.Log(1.5) - epsilon, 1, 0, 1.5);
                double likelihood = d.LogLikelihood(FittingSample);
                Assert.IsTrue(Tools.IsFinite(likelihood) && likelihood > previous);
                previous = likelihood;
            }
        }

        /// <summary>Near-zero shapes retain their continuous Gumbel limits.</summary>
        [TestMethod]
        public void NearZeroShapesPreserveQuantilesAndProbabilities()
        {
            foreach (double k in new[] { -1E-16, 0d, 1E-16 })
            foreach (double h in new[] { -1E-16, 0d, 1E-16 })
            {
                var d = new KappaFour(0, 1, k, h);
                Assert.AreEqual(0.366512920581664327d, d.InverseCDF(0.5), 2E-14);
                Assert.AreEqual(0.5, d.CDF(0.366512920581664327d), 2E-14);
                Assert.AreEqual(0.346573590279972655d, d.PDF(0.366512920581664327d), 2E-14);
                Assert.AreEqual(27.631043237892857d, d.InverseCDF(0.999999999999), 2E-12);
            }
            Assert.AreEqual(Math.Log(0.2), new KappaFour(0, 1, 1E-16, 0.2).Minimum, 2E-14);
        }

        /// <summary>Log probabilities remain finite beyond the range of ordinary probabilities.</summary>
        [TestMethod]
        public void LogProbabilitiesDoNotRoundThroughDensityOrCdf()
        {
            var gumbel = new KappaFour(0, 1, 0, 0);
            Assert.AreEqual(-1089.6331584284585, gumbel.LogPDF(-7), 2E-12);
            Assert.AreEqual(-40, gumbel.LogCCDF(40), 2E-14);
            Assert.AreEqual(4.248354255291589E-18, gumbel.CCDF(40), 1E-31);
            Assert.AreEqual(0, gumbel.PDF(double.NegativeInfinity));
            Assert.AreEqual(0, gumbel.PDF(-1000));
            var logistic = new KappaFour(0, 1, 0, -1);
            Assert.AreEqual(-1000, logistic.LogPDF(-1000), 1E-12);
            Assert.AreEqual(-1000, logistic.LogCDF(-1000), 1E-12);
        }

        /// <summary>Finite extreme inputs do not overflow intermediate affine or latent transforms.</summary>
        [TestMethod]
        public void ExtremeFiniteParametersRetainRepresentableResults()
        {
            var bounded = new KappaFour(-1E308, 1E308, 0.5, 0);
            Assert.AreEqual(1E308, bounded.Maximum, 1E293);
            Assert.AreEqual(1d, bounded.CDF(1.5E308));
            Assert.AreEqual(0d, bounded.PDF(1.5E308));
            var gumbel = new KappaFour(-1E308, 1E308, 0, 0);
            Assert.AreEqual(1E308, gumbel.InverseCDF(Math.Exp(-Math.Exp(-2))), 1E294);
            Assert.AreEqual(Math.Log(1E308), new KappaFour(0, 1, 0, 1E308).InverseCDF(1E-300), 2E-12);
            Assert.AreEqual(-1d, new KappaFour(0, 1, -1, -1E200).QuantileGradient(0.5)[2]);
            Assert.AreEqual(double.NegativeInfinity, new KappaFour(0, 1, 0, -1E-308).LogPDF(-1000));
            Assert.AreEqual(0d, new KappaFour(0, 1, 0, -1E-308).PDF(-1000));
        }

        /// <summary>Compensated support residuals retain probabilities within a few floating-point steps of either endpoint.</summary>
        [TestMethod]
        public void BoundaryLogProbabilitiesMatchIndependentDecimalOracle()
        {
            var assembly = typeof(Test_KappaFourRegression).Assembly;
            using var reader = new StreamReader(assembly.GetManifestResourceStream(assembly.GetManifestResourceNames().Single(n => n.EndsWith("kappa-four-boundaries.csv"))));
            reader.ReadLine();
            int count = 0;
            string line;
            while ((line = reader.ReadLine()) != null)
            {
                double[] r = line.Split(',').Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
                var d = new KappaFour(0, 1, r[0], r[1]);
                string context = $"k={r[0]:R}, h={r[1]:R}, x={r[2]:R}";
                Assert.AreEqual(r[3], d.LogCDF(r[2]), 5E-10, context);
                Assert.AreEqual(r[4], d.LogPDF(r[2]), 5E-10, context);
                count++;
            }
            Assert.AreEqual(45, count);
            Assert.AreEqual(-74.85989550047409, new KappaFour(0, 1, 2, 0.5).LogCDF(-1.4999999999999998), 5E-12);
        }

        /// <summary>Support endpoint densities equal their one-sided mathematical limits.</summary>
        [TestMethod]
        public void EndpointDensitiesHaveDefinedLimits()
        {
            foreach (var pair in new[] { (-0.5, -1d, 0d), (-1d, -1d, 1d), (-2d, -1d, double.PositiveInfinity), (-1d, 0d, 0d) })
            {
                var d = new KappaFour(0, 1, pair.Item1, pair.Item2);
                Assert.AreEqual(pair.Item3, d.PDF(d.Minimum));
            }
            var singular = new KappaFour(0, 1, 0, 1.5);
            Assert.AreEqual(double.PositiveInfinity, singular.PDF(singular.Minimum));
            Assert.AreEqual(double.PositiveInfinity, singular.LogPDF(singular.Minimum));
            Assert.AreEqual(1d, new KappaFour(0, 1, 1, 0).PDF(1));
        }

        /// <summary>Densities normalize across signs and zero-shape reductions.</summary>
        [TestMethod]
        public void DensityIntegratesToCdfMass()
        {
            foreach (var shapes in new[] { (-0.2, -0.5), (0.2, -0.5), (-0.2, 0.5), (0.2, 0.5), (0d, -1d), (0d, 0d), (0d, 1d), (1d, 1d) })
            {
                var d = new KappaFour(2.5, 1.75, shapes.Item1, shapes.Item2);
                var integral = new AdaptiveGaussKronrod(d.PDF, d.InverseCDF(1E-6), d.InverseCDF(1 - 1E-6));
                integral.Integrate();
                Assert.AreEqual(0.999998, integral.Result, 1E-8);
            }
        }

        /// <summary>Gumbel and exponential L-moments have their analytical values.</summary>
        [TestMethod]
        public void LinearMomentsIncludeExactZeroShapes()
        {
            var d = new KappaFour();
            AssertVector(new[] { 0.5772156649015329, 0.6931471805599453, 0.16992500144231236, 0.15037499278843736 }, d.LinearMomentsFromParameters(new[] { 0d, 1d, 0d, 0d }), 1E-9);
            AssertVector(new[] { 1d, 0.5, 1d / 3d, 1d / 6d }, d.LinearMomentsFromParameters(new[] { 0d, 1d, 0d, 1d }), 1E-9);
            foreach (double k in new[] { -1E-12, 1E-12 })
                AssertVector(new[] { 0.5772156649015329, 0.6931471805599453, 0.16992500144231236, 0.15037499278843736 }, d.LinearMomentsFromParameters(new[] { 0d, 1d, k, 1E-12 }), 1E-8);
            var smallH = d.LinearMomentsFromParameters(new[] { 0d, 1d, 0.2, 0.01 });
            Assert.IsTrue(Array.TrueForAll(smallH, Tools.IsFinite));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => d.LinearMomentsFromParameters(new[] { 0d, 1d, -1d, 0d }));
        }

        /// <summary>Inverse L-moments solve exact limiting cases without singular Newton expressions.</summary>
        [TestMethod]
        public void LinearMomentFitsHandleGumbelAndExponential()
        {
            var d = new KappaFour();
            foreach (double[] moments in new[] { new[] { 1d, 0.5, 1d / 3d, 1d / 6d }, new[] { 0.5772156649015329, 0.6931471805599453, 0.16992500144231236, 0.15037499278843736 } })
                AssertVector(moments, d.LinearMomentsFromParameters(d.ParametersFromLinearMoments(moments)), 1E-6);
        }

        /// <summary>Iteration exhaustion cannot return an unconverged L-moment fit.</summary>
        [TestMethod]
        public void LinearMomentFitRejectsExhaustedIteration()
        {
            Assert.ThrowsExactly<InvalidOperationException>(() => new KappaFour().ParametersFromLinearMoments(new[] { 0d, 1d, -0.9, 0.83375 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new KappaFour().ParametersFromLinearMoments(new[] { double.NaN, 1d, 0d, 0.1 }));
        }

        /// <summary>Finite moments are accurate and divergent moments are reported as undefined.</summary>
        [TestMethod]
        public void MomentsRespectExistenceAndKnownValues()
        {
            var heavy = new KappaFour(0, 1, -0.5, 0);
            Assert.AreEqual(1.544907701811032, heavy.Mean, 1E-8);
            Assert.IsTrue(double.IsNaN(heavy.StandardDeviation));
            Assert.IsTrue(double.IsNaN(heavy.Skewness));
            Assert.IsTrue(double.IsNaN(heavy.Kurtosis));
            var exponential = new KappaFour(0, 1, 0, 1);
            AssertVector(new[] { 1d, 1d, 2d, 9d }, new[] { exponential.Mean, exponential.StandardDeviation, exponential.Skewness, exponential.Kurtosis }, 1E-8);
            exponential.Xi = 3;
            exponential.Alpha = 2;
            Assert.AreEqual(5d, exponential.Mean, 1E-9);
            Assert.AreEqual(2d, exponential.StandardDeviation, 1E-9);
            Assert.IsTrue(double.IsNaN(new KappaFour(0, 1, 1, -1).Mean));
        }

        /// <summary>Modes include boundaries, interior stationary points, and nonunique cases.</summary>
        [TestMethod]
        public void ModesIncludeSupportEndpoints()
        {
            Assert.AreEqual(0d, new KappaFour(0, 1, 0, 1).Mode);
            Assert.AreEqual(1d, new KappaFour(0, 1, 1, 0).Mode);
            Assert.AreEqual(2.5, new KappaFour(2.5, 1.75, 0, 0).Mode, 1E-14);
            Assert.IsTrue(double.IsNaN(new KappaFour(0, 1, 1, 1).Mode));
            Assert.IsTrue(double.IsNaN(new KappaFour(0, 1, 2, 2).Mode));
            Assert.AreEqual(1E-20, new KappaFour(0, 1, 1E-20, 0).Mode, 1E-35);
        }

        /// <summary>Analytical generalized Pareto and reverse exponential L-moments verify gamma ratios independently.</summary>
        [TestMethod]
        public void LinearMomentsMatchAnalyticalFamilies()
        {
            var d = new KappaFour();
            foreach (double h in new[] { -0.5, -0.01, 0d, 1E-12, 0.001, 0.01, 0.5, 1d, 2d })
                AssertVector(new[] { h / (h + 1), 1 / ((h + 1) * (h + 2)), (h - 1) / (h + 3),
                    (h - 1) * (h - 2) / ((h + 3) * (h + 4)) }, d.LinearMomentsFromParameters(new[] { 0d, 1d, 1d, h }), 2E-9);
            foreach (double k in new[] { -0.5, -1E-12, 0d, 1E-12, 0.2, 1d, 2d })
                AssertVector(new[] { 1 / (1 + k), 1 / ((1 + k) * (2 + k)), (1 - k) / (3 + k),
                    (1 - k) * (2 - k) / ((3 + k) * (4 + k)) }, d.LinearMomentsFromParameters(new[] { 0d, 1d, k, 1d }), 2E-9);
        }

        /// <summary>Nonfinite and degenerate data propagate failure while retaining the original distribution.</summary>
        [TestMethod]
        public void InvalidFitsPreserveExistingParameters()
        {
            var d = new KappaFour(2, 3, 0.2, 0.5);
            double[] original = d.GetParameters;
            foreach (var method in new[] { ParameterEstimationMethod.MethodOfLinearMoments, ParameterEstimationMethod.MaximumLikelihood })
            {
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => d.Estimate(new[] { 1d, 2d, 3d, double.NaN }, method));
                AssertVector(original, d.GetParameters, 0d);
            }
            Assert.ThrowsExactly<InvalidOperationException>(() => d.Estimate(new[] { 1d, 1d, 1d, 1d }, ParameterEstimationMethod.MaximumLikelihood));
            AssertVector(original, d.GetParameters, 0d);
        }

        /// <summary>Quantile derivatives have finite, accurate zero-shape limits.</summary>
        [TestMethod]
        public void QuantileGradientsIncludeZeroShapeLimits()
        {
            // Q=xi-alpha*log(-log(p)); shape derivatives follow the Taylor limits.
            double q = 0.366512920581664327d;
            double[] expected = { 1, q, -0.5 * q * q, 0.346573590279972655d };
            foreach (double k in new[] { -1E-12, 0d, 1E-12 })
            foreach (double h in new[] { -1E-12, 0d, 1E-12 })
                AssertVector(expected, new KappaFour(0, 1, k, h).QuantileGradient(0.5), 1E-10);
            var d = new KappaFour(0, 1, 0, 0);
            var j = d.QuantileJacobian(new[] { 0.1, 0.3, 0.6, 0.9 }, out double determinant);
            Assert.IsTrue(Tools.IsFinite(determinant) && determinant != 0);
            Assert.AreEqual(1d, j[0, 0]);
        }

        private static void AssertVector(double[] expected, double[] actual, double tolerance)
        {
            Assert.HasCount(expected.Length, actual);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i], tolerance, $"component {i}");
        }
    }
}
