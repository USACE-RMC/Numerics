using System;
using System.Collections.Generic;
using System.Reflection;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit tests for the NUTS diagonal mass-matrix adaptation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The oracle for these tests is analytic. The target is a zero-mean diagonal Gaussian whose
    /// per-coordinate standard deviations are chosen in advance, so the exact posterior standard
    /// deviation and variance of every coordinate are known in closed form and no reference from
    /// another library is involved. The adapted metric is read back through reflection because the
    /// sampler does not expose it; that is deliberate, so the assertions do not require widening the
    /// public surface of <see cref="NUTS"/>.
    /// </para>
    /// <para>
    /// The scale range is four orders of magnitude in standard deviation, eight in variance. An
    /// identity metric cannot cover that range with one step size, so a working adaptation is
    /// visible both in the recovered spread and in the per-transition leapfrog cost.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_NUTS_MassMatrix
    {
        /// <summary>
        /// Zero-mean diagonal Gaussian with log-spaced per-coordinate standard deviations.
        /// </summary>
        private sealed class DiagonalGaussian
        {
            /// <summary>The exact posterior standard deviation of each coordinate.</summary>
            public readonly double[] StandardDeviations;

            /// <summary>
            /// Creates the target.
            /// </summary>
            /// <param name="dimension">The number of coordinates.</param>
            /// <param name="smallestSd">The standard deviation of the first coordinate.</param>
            /// <param name="largestSd">The standard deviation of the last coordinate.</param>
            public DiagonalGaussian(int dimension, double smallestSd, double largestSd)
            {
                StandardDeviations = new double[dimension];
                double lo = Math.Log10(smallestSd);
                double hi = Math.Log10(largestSd);
                for (int j = 0; j < dimension; j++)
                    StandardDeviations[j] = Math.Pow(10.0, lo + (hi - lo) * j / (dimension - 1.0));
            }

            /// <summary>The log-likelihood, up to an additive constant.</summary>
            /// <param name="x">The parameter vector.</param>
            /// <returns>The log-density.</returns>
            public double LogLikelihood(double[] x)
            {
                double sum = 0d;
                for (int j = 0; j < StandardDeviations.Length; j++)
                {
                    double z = x[j] / StandardDeviations[j];
                    sum += -0.5 * z * z;
                }
                return sum;
            }

            /// <summary>The exact gradient of the log-likelihood.</summary>
            /// <param name="x">The parameter vector.</param>
            /// <returns>The gradient.</returns>
            public Vector Gradient(IList<double> x)
            {
                var g = new Vector(StandardDeviations.Length);
                for (int j = 0; j < StandardDeviations.Length; j++)
                    g[j] = -x[j] / (StandardDeviations[j] * StandardDeviations[j]);
                return g;
            }
        }

        /// <summary>
        /// Builds a single-chain NUTS sampler on the diagonal Gaussian with a deterministic start.
        /// </summary>
        /// <param name="target">The target.</param>
        /// <param name="halfWidth">The half-width of the uniform prior on every coordinate.</param>
        /// <param name="adapt">Whether to adapt the mass matrix.</param>
        /// <param name="warmup">The number of warmup iterations.</param>
        /// <param name="iterations">The number of iterations.</param>
        /// <param name="maxTreeDepth">The maximum tree depth.</param>
        /// <returns>An unsampled sampler.</returns>
        private static NUTS BuildSampler(DiagonalGaussian target, double halfWidth, bool adapt,
                                         int warmup, int iterations, int maxTreeDepth)
        {
            var priors = new List<IUnivariateDistribution>();
            for (int j = 0; j < target.StandardDeviations.Length; j++)
                priors.Add(new Uniform(-halfWidth, halfWidth));

            return new NUTS(priors, target.LogLikelihood, maxTreeDepth: maxTreeDepth,
                            gradientFunction: target.Gradient)
            {
                NumberOfChains = 1,
                ParallelizeChains = false,
                // With one chain and one initial iteration the chain starts at the prior mean,
                // which is the mode of this target, so the run is fully deterministic.
                InitialIterations = 1,
                ThinningInterval = 1,
                WarmupIterations = warmup,
                Iterations = iterations,
                OutputLength = 2000,
                PRNGSeed = 12345,
                AdaptMassMatrix = adapt
            };
        }

        /// <summary>
        /// Reads one of the sampler's private diagonal metric arrays for the given chain.
        /// </summary>
        /// <param name="sampler">The sampler.</param>
        /// <param name="fieldName">The private field name.</param>
        /// <returns>The metric diagonal of chain zero.</returns>
        private static double[] ReadMetric(NUTS sampler, string fieldName)
        {
            var field = typeof(NUTS).GetField(fieldName, BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.IsNotNull(field, "The private field " + fieldName + " was not found on NUTS.");
            var metric = field!.GetValue(sampler) as double[][];
            Assert.IsNotNull(metric, "The private field " + fieldName + " was not a double[][].");
            return metric![0];
        }

        /// <summary>
        /// Computes the sample standard deviation of one coordinate of the recorded output.
        /// </summary>
        /// <param name="sampler">A sampled sampler.</param>
        /// <param name="coordinate">The coordinate index.</param>
        /// <returns>The sample standard deviation.</returns>
        private static double RecoveredStandardDeviation(NUTS sampler, int coordinate)
        {
            var draws = sampler.Output[0];
            double mean = 0d;
            for (int i = 0; i < draws.Count; i++) mean += draws[i].Values[coordinate];
            mean /= draws.Count;
            double m2 = 0d;
            for (int i = 0; i < draws.Count; i++)
            {
                double e = draws[i].Values[coordinate] - mean;
                m2 += e * e;
            }
            return Math.Sqrt(m2 / (draws.Count - 1));
        }

        /// <summary>
        /// With mass-matrix adaptation enabled, the sampler must recover the analytically known
        /// posterior standard deviation of every coordinate across four orders of magnitude.
        /// </summary>
        /// <remarks>
        /// The stated band is recovered sd / true sd in [0.75, 1.30] on all twenty coordinates.
        /// It is wide enough to absorb single-chain Monte Carlo error on 2,000 draws and any
        /// cross-framework divergence of this chaotic trajectory, and narrow enough that the
        /// pre-fix sampler fails it: with the additive prior-scaled floor in place the widest
        /// coordinate returns about 0.33 to 0.37 of its true width.
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_AdaptedMetric_RecoversKnownScales()
        {
            var target = new DiagonalGaussian(20, 1e-3, 1e1);
            var sampler = BuildSampler(target, 1000d, adapt: true, warmup: 600, iterations: 1500, maxTreeDepth: 10);
            sampler.Sample();

            for (int j = 0; j < target.StandardDeviations.Length; j++)
            {
                double ratio = RecoveredStandardDeviation(sampler, j) / target.StandardDeviations[j];
                Assert.IsTrue(ratio >= 0.75 && ratio <= 1.30,
                    $"Coordinate {j} (true sd {target.StandardDeviations[j]:E3}) recovered at {ratio:F4} of its true width; expected [0.75, 1.30].");
            }

            Assert.AreEqual(0, sampler.DivergenceCounts[0], "The adapted run must not diverge on a Gaussian target.");
        }

        /// <summary>
        /// The adapted metric must track the true per-coordinate variances rather than collapsing
        /// toward isotropy, and it must do so in the sense the sampler consumes it.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The estimated posterior variance belongs in the inverse mass: momentum is drawn with
        /// standard deviation sqrt(M), so a coordinate of posterior standard deviation s must carry
        /// M = 1 / s^2 for its leapfrog step to scale like s. This test asserts that directly, within
        /// a factor of four on the variance, which is a factor of two on the standard deviation.
        /// </para>
        /// <para>
        /// It fails on all three of the broken arrangements this fix replaces: an additive
        /// prior-scaled floor drives the whole diagonal to one value, storing the variance in the
        /// mass instead of the inverse mass reverses the ordering, and doing both leaves a metric
        /// anti-correlated with the truth.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_AdaptedMetric_TracksTrueVariances()
        {
            var target = new DiagonalGaussian(20, 1e-3, 1e1);
            var sampler = BuildSampler(target, 1000d, adapt: true, warmup: 600, iterations: 1500, maxTreeDepth: 10);
            sampler.Sample();

            var inverseMass = ReadMetric(sampler, "_inverseMassMatrix");
            var mass = ReadMetric(sampler, "_massMatrix");

            double smallest = double.MaxValue, largest = 0d;
            for (int j = 0; j < target.StandardDeviations.Length; j++)
            {
                double trueVariance = target.StandardDeviations[j] * target.StandardDeviations[j];
                Assert.IsTrue(inverseMass[j] >= 0.25 * trueVariance && inverseMass[j] <= 4.0 * trueVariance,
                    $"Coordinate {j}: inverse mass {inverseMass[j]:E4} is not within a factor of four of the true variance {trueVariance:E4}.");
                Assert.AreEqual(1.0 / inverseMass[j], mass[j], 1e-12 * Math.Abs(1.0 / inverseMass[j]),
                    $"Coordinate {j}: the mass and inverse mass are not reciprocals.");

                if (inverseMass[j] < smallest) smallest = inverseMass[j];
                if (inverseMass[j] > largest) largest = inverseMass[j];
            }

            // The true variances span 1e-8 of their own range. A metric flattened toward isotropy
            // by an additive floor spans less than 1.02.
            Assert.IsGreaterThan(1e4, largest / smallest,
                $"The adapted metric spans only a factor of {largest / smallest:E3}; it has collapsed toward isotropy.");
        }

        /// <summary>
        /// Adaptation must make the ill-conditioned target dramatically cheaper per transition than
        /// the identity metric does.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_AdaptedMetric_ReducesLeapfrogCost()
        {
            var target = new DiagonalGaussian(20, 1e-3, 1e1);

            var adapted = BuildSampler(target, 1000d, adapt: true, warmup: 600, iterations: 1500, maxTreeDepth: 10);
            adapted.Sample();

            // With an identity metric this target saturates the tree-depth cap on essentially every
            // transition, which is 2^10 - 1 = 1023 leapfrog steps.
            Assert.IsLessThan(100d, adapted.MeanLeapfrogSteps[0],
                $"The adapted sampler used {adapted.MeanLeapfrogSteps[0]:F2} leapfrog steps per transition; an identity metric costs on the order of 1,000 here.");
            Assert.IsLessThan(adapted.DiagnosticSampleCounts[0], adapted.MaxTreeDepthHitCounts[0] * 100,
                $"The adapted sampler exhausted the tree-depth cap on {adapted.MaxTreeDepthHitCounts[0]} of {adapted.DiagnosticSampleCounts[0]} transitions.");
        }

        /// <summary>
        /// A window holding fewer than the minimum number of draws cannot estimate a variance, so
        /// the prior-scaled fallback must be used for every coordinate.
        /// </summary>
        /// <remarks>
        /// With 12 warmup iterations the initial and terminal buffers are one iteration each and the
        /// single adaptation window closes at iteration 10, having accumulated nine draws. That is
        /// below the ten-draw minimum, so the metric must be exactly the fallback variance
        /// (priorRange / 6)^2 on every coordinate.
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_ShortAdaptationWindow_UsesFallbackVariance()
        {
            var target = new DiagonalGaussian(4, 1e-2, 1e0);
            var sampler = BuildSampler(target, 1000d, adapt: true, warmup: 12, iterations: 100, maxTreeDepth: 6);
            sampler.Sample();

            double priorRange = 2000d;
            double expected = (priorRange * priorRange) / 36.0;
            var inverseMass = ReadMetric(sampler, "_inverseMassMatrix");
            for (int j = 0; j < inverseMass.Length; j++)
            {
                Assert.AreEqual(expected, inverseMass[j], 1e-9 * expected,
                    $"Coordinate {j} did not fall back to the prior-scaled variance on a nine-draw window.");
            }
        }

        /// <summary>
        /// A window holding at least the minimum number of draws must use its own variance rather
        /// than the prior-scaled fallback.
        /// </summary>
        /// <remarks>
        /// With 14 warmup iterations the single adaptation window closes at iteration 12 with ten
        /// accumulated draws, which is exactly the minimum. The posterior of this target is many
        /// orders of magnitude narrower than the prior, so any use of the window's own variance is
        /// unmistakable.
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_SufficientAdaptationWindow_UsesWindowVariance()
        {
            var target = new DiagonalGaussian(4, 1e-2, 1e0);
            var sampler = BuildSampler(target, 1000d, adapt: true, warmup: 14, iterations: 100, maxTreeDepth: 6);
            sampler.Sample();

            double fallback = (2000d * 2000d) / 36.0;
            var inverseMass = ReadMetric(sampler, "_inverseMassMatrix");
            for (int j = 0; j < inverseMass.Length; j++)
            {
                Assert.IsTrue(inverseMass[j] > 0d && inverseMass[j] < 1d,
                    $"Coordinate {j} returned {inverseMass[j]:E4}; a ten-draw window on this target cannot produce a variance of order the prior's {fallback:E4}.");
            }
        }

        /// <summary>
        /// With <see cref="NUTS.AdaptMassMatrix"/> disabled the sampler must reproduce the exact
        /// draws it produced before the mass-matrix adaptation was corrected.
        /// </summary>
        /// <remarks>
        /// The reference values were captured from the unmodified sampler and are asserted to a
        /// relative tolerance of 1e-12, which is far tighter than any behavioural change could be
        /// and loose enough to absorb the last-place differences the transcendental functions show
        /// between .NET Framework and .NET. The disabled path never enters the adaptation code, so
        /// this is a guard that the change stayed inside it.
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_AdaptationDisabled_ReproducesCapturedDraws()
        {
            var target = new DiagonalGaussian(6, 1e-2, 1e1);
            var priors = new List<IUnivariateDistribution>();
            for (int j = 0; j < 6; j++) priors.Add(new Uniform(-1000d, 1000d));
            var sampler = new NUTS(priors, target.LogLikelihood, maxTreeDepth: 6, gradientFunction: target.Gradient)
            {
                NumberOfChains = 1,
                ParallelizeChains = false,
                InitialIterations = 1,
                ThinningInterval = 1,
                WarmupIterations = 60,
                Iterations = 150,
                OutputLength = 100,
                PRNGSeed = 12345,
                AdaptMassMatrix = false
            };
            sampler.Sample();

            AssertClose(0.010797384971651152, sampler.StepSizes[0], "step size");
            AssertClose(54.85263157894737, sampler.MeanLeapfrogSteps[0], "mean leapfrog steps");

            double[] firstDraw = { 0.0011958505341814842, -0.07526409822522193, -0.21184584905057402,
                                   -0.6582920358936465, 2.1468119667523973, -0.6069130233917713 };
            double[] lastDraw = { 0.004230279378397197, -0.009370317476825615, -0.15308785669000735,
                                  -0.726029878587165, 3.652148817280656, -4.4283508900986215 };
            var last = sampler.Output[0][sampler.Output[0].Count - 1];
            for (int j = 0; j < 6; j++)
            {
                AssertClose(firstDraw[j], sampler.Output[0][0].Values[j], $"first draw coordinate {j}");
                AssertClose(lastDraw[j], last.Values[j], $"last draw coordinate {j}");
            }

            // The metric must still be the identity that the constructor installed.
            var inverseMass = ReadMetric(sampler, "_inverseMassMatrix");
            var mass = ReadMetric(sampler, "_massMatrix");
            for (int j = 0; j < 6; j++)
            {
                Assert.AreEqual(1d, inverseMass[j], 0d, $"Coordinate {j} inverse mass was modified with adaptation disabled.");
                Assert.AreEqual(1d, mass[j], 0d, $"Coordinate {j} mass was modified with adaptation disabled.");
            }
        }

        /// <summary>
        /// Asserts a value matches a captured reference to a relative tolerance of 1e-12.
        /// </summary>
        /// <param name="expected">The captured reference value.</param>
        /// <param name="actual">The value produced by this run.</param>
        /// <param name="what">A description used in the failure message.</param>
        private static void AssertClose(double expected, double actual, string what)
        {
            Assert.AreEqual(expected, actual, 1e-12 * Math.Abs(expected),
                $"The disabled-adaptation path changed: {what} was {actual:R}, expected {expected:R}.");
        }
    }
}
