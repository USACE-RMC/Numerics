using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Functions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;

namespace Correctness
{
    /// <summary>
    /// Regression tests for function state, validation, and serialization behavior.
    /// </summary>
    [TestClass]
    public class Test_FunctionCorrectnessRepairs
    {
        /// <summary>Verifies deterministic restoration and atomic validation for segmented power functions.</summary>
        [TestMethod]
        public void SegmentedPower_DeterministicXmlAndAtomicValidation()
        {
            var function = new SegmentedPowerFunction(1) { IsDeterministic = true, Maximum = 20d };
            function.SetParameters(new[] { 1d, 0.5d, 2d, 0d });
            Assert.IsTrue(function.ParametersValid);

            var restored = (SegmentedPowerFunction)UnivariateFunctionFactory.CreateFromXElement(function.ToXElement());
            Assert.IsTrue(restored.IsDeterministic);
            Assert.AreEqual(20d, restored.Maximum, 0d);
            Assert.AreEqual(function.Function(4d), restored.Function(4d), 0d);

            double originalBeta = function.GetBeta(1);
            function.SetParameters(new[] { 1d, 0.5d, 0d, 0d });
            Assert.IsTrue(function.ParametersValid, "A rejected update must leave the prior state valid.");
            Assert.AreEqual(originalBeta, function.GetBeta(1), 0d);

            Assert.Throws<ArgumentOutOfRangeException>(() => function.Maximum = function.Minimum);
            Assert.AreEqual(20d, function.Maximum, 0d);
            Assert.Throws<ArgumentOutOfRangeException>(() => new SegmentedPowerFunction(new[] { 1d, 0.5d, 0d, 0.1d }));

            XElement malformed = function.ToXElement();
            malformed.SetAttributeValue(nameof(SegmentedPowerFunction.Maximum), "0");
            Assert.Throws<ArgumentOutOfRangeException>(() => SegmentedPowerFunction.FromXElement(malformed));
        }

        /// <summary>Verifies that composite evaluation restores child confidence state after an exception.</summary>
        [TestMethod]
        public void Composite_RestoresChildStateAfterExceptions()
        {
            var child = new ConfidenceProbeFunction { ConfidenceLevel = 0.35d, ThrowOnEvaluation = true };
            var composite = new CompositeFunction(new IUnivariateFunction[] { child }) { ConfidenceLevel = 0.8d };

            Assert.Throws<InvalidOperationException>(() => composite.Function(1d));
            Assert.AreEqual(0.35d, child.ConfidenceLevel, 0d);
            Assert.Throws<InvalidOperationException>(() => composite.InverseFunction(1d));
            Assert.AreEqual(0.35d, child.ConfidenceLevel, 0d);
        }

        /// <summary>Verifies concurrent stateful evaluation and the lock-free deterministic path.</summary>
        [TestMethod]
        public void Composite_SynchronizesOnlyStatefulChildren()
        {
            var child = new ConfidenceProbeFunction { ConfidenceLevel = 0.5d };
            var lower = new CompositeFunction(new IUnivariateFunction[] { child }) { ConfidenceLevel = 0.2d };
            var upper = new CompositeFunction(new IUnivariateFunction[] { child }) { ConfidenceLevel = 0.8d };
            int failures = 0;

            Parallel.For(0, 400, i =>
            {
                double expected = i % 2 == 0 ? 0.2d : 0.8d;
                double actual = i % 2 == 0 ? lower.Function(0d) : upper.Function(0d);
                if (actual != expected) Interlocked.Increment(ref failures);
            });

            Assert.AreEqual(0, failures);
            Assert.AreEqual(0.5d, child.ConfidenceLevel, 0d);

            child.IsDeterministic = true;
            child.ThrowOnConfidenceAssignment = true;
            Assert.AreEqual(0.5d, lower.Function(0d), 0d,
                "Deterministic child evaluation must not mutate confidence state.");
        }

        /// <summary>Verifies atomic weight rejection and composite mode validation.</summary>
        [TestMethod]
        public void Composite_RejectsInvalidStateAtomically()
        {
            var composite = new CompositeFunction(
                new IUnivariateFunction[] { new LinearFunction(), new LinearFunction(1d, 2d) },
                new[] { 0.25d, 0.75d });
            composite.SetParameters(new[] { double.NaN, 0.75d });
            Assert.IsTrue(composite.ParametersValid);
            Assert.AreEqual(0.25d, composite.Weights[0], 0d);
            Assert.Throws<ArgumentOutOfRangeException>(() => composite.Mode = (CompositeFunctionMode)999);

            XElement malformed = composite.ToXElement();
            malformed.SetAttributeValue(nameof(CompositeFunction.Mode), "999");
            Assert.Throws<ArgumentException>(() => CompositeFunction.FromXElement(malformed));
        }

        /// <summary>Verifies ensemble parameter ownership and serialized-set validation.</summary>
        [TestMethod]
        public void Ensemble_OwnsDeepCopiesAndValidatesXmlSets()
        {
            var values = new[] { 1d, 0.5d, 2d, 0.1d };
            var ensemble = new EnsembleFunction(
                new SegmentedPowerFunction(values),
                new[] { new ParameterSet(values, 1d, 0.5d) });

            values[0] = 99d;
            Assert.AreEqual(1d, ((SegmentedPowerFunction)ensemble.Sample(0)).GetBreakpoint(1), 0d);

            ParameterSet exposed = ensemble.ParameterSets[0];
            exposed.Values[0] = 88d;
            Assert.AreEqual(1d, ((SegmentedPowerFunction)ensemble.Sample(0)).GetBreakpoint(1), 0d);
            Assert.Throws<ArgumentOutOfRangeException>(() => ensemble.Sample(double.NaN));

            XElement invalidValues = ensemble.ToXElement();
            invalidValues.Element("ParameterSets")!.Element(nameof(ParameterSet))!
                .SetAttributeValue(nameof(ParameterSet.Values), "1|0.5|0|0.1");
            Assert.Throws<ArgumentException>(() => EnsembleFunction.FromXElement(invalidValues));

            XElement invalidFitness = ensemble.ToXElement();
            invalidFitness.Element("ParameterSets")!.Element(nameof(ParameterSet))!
                .SetAttributeValue(nameof(ParameterSet.Fitness), "NaN");
            Assert.Throws<ArgumentException>(() => EnsembleFunction.FromXElement(invalidFitness));
        }

        private sealed class ConfidenceProbeFunction : IUnivariateFunction
        {
            private double _confidenceLevel = -1d;

            public bool ThrowOnEvaluation { get; set; }
            public bool ThrowOnConfidenceAssignment { get; set; }
            public int NumberOfParameters => 0;
            public bool ParametersValid => true;
            public double Minimum { get; set; } = double.MinValue;
            public double Maximum { get; set; } = double.MaxValue;
            public double[] MinimumOfParameters => Array.Empty<double>();
            public double[] MaximumOfParameters => Array.Empty<double>();
            public bool IsDeterministic { get; set; }
            public double ConfidenceLevel
            {
                get { return _confidenceLevel; }
                set
                {
                    if (ThrowOnConfidenceAssignment) throw new InvalidOperationException("Confidence assignment was not expected.");
                    _confidenceLevel = value;
                }
            }

            public void SetParameters(IList<double> parameters)
            {
                if (parameters == null) throw new ArgumentNullException(nameof(parameters));
                if (parameters.Count != 0) throw new ArgumentException("This test function has no parameters.", nameof(parameters));
            }

            public ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
            {
                return null;
            }

            public double Function(double x)
            {
                if (ThrowOnEvaluation) throw new InvalidOperationException("Test evaluation failure.");
                double first = ConfidenceLevel;
                Thread.SpinWait(10000);
                if (first != ConfidenceLevel) throw new InvalidOperationException("Confidence state changed during evaluation.");
                return ConfidenceLevel;
            }

            public double InverseFunction(double y)
            {
                if (ThrowOnEvaluation) throw new InvalidOperationException("Test inverse failure.");
                return ConfidenceLevel;
            }
        }
    }

    /// <summary>
    /// Regression tests for probability and jackknife edge cases.
    /// </summary>
    [TestClass]
    public class Test_StatisticsCorrectnessRepairs
    {
        /// <summary>Verifies residual mass when capped enumeration omits the no-event row.</summary>
        [TestMethod]
        public void LazyExclusive_CappedWithoutNoEventRowPreservesUnionMass()
        {
            double[] probabilities = { 0.4d, 0.35d, 0.3d, 0.25d, 0.2d, 0.15d };
            var output = new List<double>();
            var indicators = new List<int[]>();

            var status = Probability.IndependentExclusiveLazy(probabilities, output, indicators,
                includeNoEventRow: false, maxEmittedCombinations: 10,
                absoluteTolerance: 0d, relativeTolerance: 0d);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Capped, status);
            double noEventMass = probabilities.Aggregate(1d, (mass, probability) => mass * (1d - probability));
            Assert.AreEqual(1d - noEventMass, output.Sum(), 1E-12);
            Assert.HasCount(11, output);
        }

        /// <summary>Verifies structural and tolerance validation for pooled enumeration.</summary>
        [TestMethod]
        public void PooledExclusive_RejectsMalformedMetadata()
        {
            double[] probabilities = { 0.2d, 0.3d };
            var output = new List<double>();
            var indicators = new List<int[]>();
            int[,] rows = { { 1, 0 }, { 0, 1 }, { 1, 1 } };

            Assert.Throws<ArgumentException>(() => Probability.IndependentExclusive(
                probabilities, new[] { 2 }, rows, output, indicators));
            Assert.Throws<ArgumentException>(() => Probability.IndependentExclusive(
                probabilities, new[] { 1, 2 }, rows, output, indicators));
            Assert.Throws<ArgumentOutOfRangeException>(() => Probability.IndependentExclusive(
                probabilities, new[] { 2, 1 }, rows, output, indicators, double.NaN));
        }

        /// <summary>Verifies public combination tuple validation.</summary>
        [TestMethod]
        public void NextCombination_RejectsInvalidTuples()
        {
            Assert.Throws<ArgumentNullException>(() => Factorial.NextCombination(null!, 3));
            Assert.Throws<ArgumentException>(() => Factorial.NextCombination(Array.Empty<int>(), 3));
            Assert.Throws<ArgumentException>(() => Factorial.NextCombination(new[] { 0, 0 }, 3));
            Assert.Throws<ArgumentException>(() => Factorial.NextCombination(new[] { 0, 3 }, 3));
            Assert.Throws<ArgumentOutOfRangeException>(() => Factorial.NextCombination(new[] { 0 }, -1));
        }

        /// <summary>Verifies single-element behavior and callback sample isolation.</summary>
        [TestMethod]
        public void Jackknife_PreservesSingleElementAndIsolatesCallbacks()
        {
            int calls = 0;
            double standardError = Statistics.JackKnifeStandardError(new[] { 5d }, sample =>
            {
                Interlocked.Increment(ref calls);
                return sample.Count;
            });
            Assert.AreEqual(0d, standardError, 0d);
            Assert.AreEqual(0, calls, "The single-element standard error does not evaluate an empty sample.");

            double[] single = Statistics.JackKnifeSample(new[] { 5d }, sample =>
            {
                Interlocked.Increment(ref calls);
                return sample.Count;
            });
            Assert.IsNotNull(single);
            Assert.AreEqual(0d, single[0], 0d);
            Assert.AreEqual(1, calls);

            double[] original = { 1d, 2d, 3d, 4d };
            var callbackSamples = new List<IList<double>>();
            object sync = new object();
            Statistics.JackKnifeSample(original, sample =>
            {
                lock (sync) callbackSamples.Add(sample);
                if (sample.Count > 0) sample[0] = -100d;
                return sample.Count;
            });
            CollectionAssert.AreEqual(new[] { 1d, 2d, 3d, 4d }, original);
            Assert.HasCount(original.Length, callbackSamples);
            Assert.AreEqual(original.Length, callbackSamples.Distinct().Count());
        }
    }

    /// <summary>
    /// Regression tests for bootstrap fitting and interpolation behavior.
    /// </summary>
    [TestClass]
    public class Test_BootstrapCorrectnessRepairs
    {
        /// <summary>Verifies successful-fit denominators and explicit all-failure errors.</summary>
        [TestMethod]
        public void Bootstrap_UsesOnlySuccessfulFitsAndRejectsAllFailures()
        {
            var parent = new Normal(0d, 1d);
            var analysis = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 1234);
            IUnivariateDistribution[] mixed = { new Normal(0d, 1d), null!, new Normal(1d, 1d) };

            double[] mean = analysis.ExpectedProbabilities(new[] { 0d }, mixed);
            double expected = 0.5d * (new Normal(0d, 1d).CDF(0d) + new Normal(1d, 1d).CDF(0d));
            Assert.AreEqual(expected, mean[0], 1E-14);
            Assert.Throws<InvalidOperationException>(() =>
                analysis.ExpectedProbabilities(new[] { 0d }, new IUnivariateDistribution[] { null!, null! }));

            var aggregate = Assert.Throws<AggregateException>(() => analysis.Distributions(new[]
            {
                new ParameterSet(new[] { 0d, -1d }, 0d),
                new ParameterSet(new[] { double.NaN, 1d }, 0d),
            }));
            Assert.HasCount(2, aggregate.InnerExceptions);
        }

        /// <summary>Verifies sign-preserving transforms for negative quantiles.</summary>
        [TestMethod]
        public void Bootstrap_NormalIntervalsPreserveNegativeQuantiles()
        {
            var parent = new Normal(-10d, 1d);
            var analysis = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 1234);
            IUnivariateDistribution[] fits =
            {
                new Normal(-9.5d, 1d),
                new Normal(-10d, 1.1d),
                null!,
                new Normal(-10.5d, 0.9d),
            };

            double[,] interval = analysis.NormalQuantileCI(new[] { 0.5d }, 0.1d, fits);
            Assert.IsFalse(double.IsNaN(interval[0, 0]) || double.IsInfinity(interval[0, 0]));
            Assert.IsFalse(double.IsNaN(interval[0, 1]) || double.IsInfinity(interval[0, 1]));
            Assert.IsLessThan(0d, interval[0, 0]);
            Assert.IsLessThan(0d, interval[0, 1]);
        }

        /// <summary>Verifies interpolation equivalence for sorted and unsorted ordinates.</summary>
        [TestMethod]
        public void Bootstrap_InterpolationKeepsSortedPairsTogether()
        {
            var parent = new Normal(0d, 1d);
            var analysis = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 1234);
            IUnivariateDistribution[] fits = { new Normal(0d, 1d), new Normal(1d, 2d) };
            double[] sorted = { -3d, -1d, 0d, 2d, 5d };
            double[] unsorted = { 2d, -3d, 5d, 0d, -1d };
            double[] probabilities = { 0.1d, 0.5d, 0.9d };

            double[] expected = analysis.ExpectedProbabilities(sorted, probabilities, fits);
            double[] actual = analysis.ExpectedProbabilities(unsorted, probabilities, fits);
            CollectionAssert.AreEqual(expected, actual);
        }
    }

    /// <summary>
    /// Regression tests for distribution serialization and numerical edge cases.
    /// </summary>
    [TestClass]
    public class Test_DistributionCorrectnessRepairs
    {
        /// <summary>Verifies rejection of malformed serialized distribution data.</summary>
        [TestMethod]
        public void DistributionXml_RejectsMalformedTablesEnumsAndWeights()
        {
            var empirical = new EmpiricalDistribution(new[] { 1d, 2d }, new[] { 0d, 1d });
            XElement invalidOrder = empirical.ToXElement();
            invalidOrder.SetAttributeValue("ProbabilityOrder", "999");
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(invalidOrder));

            XElement invalidTransform = empirical.ToXElement();
            invalidTransform.SetAttributeValue(nameof(EmpiricalDistribution.XTransform), "999");
            Assert.Throws<ArgumentException>(() => EmpiricalDistribution.FromXElement(invalidTransform));

            var kernel = new KernelDensity(new[] { 1d, 2d, 3d }, new[] { 1d, 1d, 1d }, KernelDensity.KernelType.Gaussian, 0.5d);
            XElement zeroWeights = kernel.ToXElement();
            zeroWeights.SetAttributeValue("Weights", "0|0|0");
            Assert.Throws<ArgumentException>(() => KernelDensity.FromXElement(zeroWeights));

            XElement invalidBandwidth = kernel.ToXElement();
            invalidBandwidth.SetAttributeValue(nameof(KernelDensity.Bandwidth), "NaN");
            Assert.Throws<ArgumentException>(() => KernelDensity.FromXElement(invalidBandwidth));

            XElement missingType = new Normal().ToXElement();
            missingType.Attribute(nameof(UnivariateDistributionBase.Type))!.Remove();
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(missingType));
        }

        /// <summary>Verifies occupied support, total mass, and mean for lattice convolution.</summary>
        [TestMethod]
        public void DiscreteConvolution_UsesOccupiedSupportAndPositiveMass()
        {
            EmpiricalDistribution.ConvolveDiscrete(
                new[] { -100d, 2d, 5d }, new[] { 0d, 0.25d, 0.75d },
                new[] { 3d, 7d, 100d }, new[] { 0.5d, 0.5d, 0d },
                256, out double[] values, out double[] masses);

            Assert.AreEqual(5d, values[0], 1E-12);
            double step = values[1] - values[0];
            Assert.IsGreaterThanOrEqualTo(12d, values[values.Length - 1]);
            Assert.IsLessThanOrEqualTo(step + 1E-12, values[values.Length - 1] - 12d,
                "The occupied lattice may extend at most one node beyond the exact support.");
            Assert.AreEqual(1d, masses.Sum(), 1E-12);
            double mean = values.Zip(masses, (value, mass) => value * mass).Sum();
            Assert.AreEqual(9.25d, mean, 1E-10);

            Assert.Throws<ArgumentException>(() => EmpiricalDistribution.ConvolveDiscrete(
                new[] { 1d, 2d }, new[] { 0d, 0d }, new[] { 1d }, new[] { 1d },
                256, out _, out _));
        }

        /// <summary>Verifies early rejection of degenerate logarithmic output support.</summary>
        [TestMethod]
        public void LogConvolution_RejectsDegenerateSupportBeforeFft()
        {
            var point1 = new EmpiricalDistribution(new[] { 2d, 2d }, new[] { 0d, 1d });
            var point2 = new EmpiricalDistribution(new[] { 3d, 3d }, new[] { 0d, 1d });
            Assert.Throws<ArgumentException>(() =>
                EmpiricalDistribution.Convolve(point1, point2, 128, logSpacedOutput: true));
        }

        /// <summary>Verifies competing-risk seed invalidation, cloning, and serialization.</summary>
        [TestMethod]
        public void CompetingRisks_SeedInvalidatesCachesAndRoundTrips()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Normal(0d, 1d),
                new Exponential(2d),
            }) { PRNGSeed = 2468 };

            FieldInfo mvnCreated = typeof(CompetingRisks).GetField("_mvnCreated", BindingFlags.Instance | BindingFlags.NonPublic)!;
            FieldInfo empiricalCreated = typeof(CompetingRisks).GetField("_empiricalCDFCreated", BindingFlags.Instance | BindingFlags.NonPublic)!;
            mvnCreated.SetValue(distribution, true);
            empiricalCreated.SetValue(distribution, true);
            distribution.PRNGSeed = 1357;
            Assert.IsFalse((bool)mvnCreated.GetValue(distribution)!);
            Assert.IsFalse((bool)empiricalCreated.GetValue(distribution)!);

            var clone = (CompetingRisks)distribution.Clone();
            Assert.AreEqual(1357, clone.PRNGSeed);
            var restored = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(distribution.ToXElement());
            Assert.AreEqual(1357, restored.PRNGSeed);

            XElement malformed = distribution.ToXElement();
            malformed.SetAttributeValue(nameof(CompetingRisks.Dependency), "999");
            Assert.Throws<ArgumentException>(() => CompetingRisks.FromXElement(malformed));

            var multivariate = new MultivariateNormal(2);
            Assert.Throws<ArgumentNullException>(() => multivariate.MVNUNI = null!);
        }
    }

    /// <summary>
    /// Regression tests for integration state and tail-focus validation.
    /// </summary>
    [TestClass]
    public class Test_IntegrationCorrectnessRepairs
    {
        /// <summary>Verifies that recorder mutation does not alter an active integration.</summary>
        [TestMethod]
        public void AdaptiveRecorder_IsSnapshottedForAnIntegration()
        {
            int calls = 0;
            AdaptiveGaussKronrod integration = null;
            integration = new AdaptiveGaussKronrod(x => x * x, 0d, 1d)
            {
                MinDepth = 2,
                Recorder = (x, weight, value) =>
                {
                    Interlocked.Increment(ref calls);
                    integration!.Recorder = null;
                },
            };

            integration.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, integration.Status);
            Assert.IsGreaterThan(1, calls, "The recorder snapshot remains active until the integration completes.");
        }

        /// <summary>Verifies finite positive tail-focus validation.</summary>
        [TestMethod]
        public void Vegas_RejectsInvalidTailFocusParametersAtomically()
        {
            var vegas = new Vegas((point, weight) => point[0] * weight, 1, new[] { 0d }, new[] { 1d });
            Assert.Throws<ArgumentOutOfRangeException>(() => vegas.TailFocusParameter = 0d);
            Assert.Throws<ArgumentOutOfRangeException>(() => vegas.TailFocusParameter = double.NaN);
            Assert.Throws<ArgumentOutOfRangeException>(() => vegas.TailFocusParameter = double.PositiveInfinity);
            Assert.AreEqual(1d, vegas.TailFocusParameter, 0d);
            Assert.Throws<ArgumentOutOfRangeException>(() => vegas.ConfigureForRareEvents(0d));
            Assert.Throws<ArgumentOutOfRangeException>(() => vegas.ConfigureForRareEvents(1d));
            Assert.Throws<ArgumentOutOfRangeException>(() => vegas.ConfigureForRareEvents(double.NaN));
        }
    }
}