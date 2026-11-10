using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Distributions;
using Numerics.Functions;
using Numerics.Mathematics.Optimization;

namespace Functions
{
    /// <summary>
    /// Unit tests for <see cref="EnsembleFunction"/>: pure index and percentile sampling of
    /// configured clones, clone independence from the template and each other, parallel
    /// no-shared-mutation behavior, the guards, and the serialization round-trip.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_EnsembleFunction
    {
        /// <summary>Builds a three-draw posterior over a segmented power template.</summary>
        private static EnsembleFunction BuildEnsemble()
        {
            var template = new SegmentedPowerFunction(new[] { 1d, 1.5d, 2d, 0.1d });
            var sets = new[]
            {
                new ParameterSet(new[] { 1.0d, 1.5d, 2.0d, 0.10d }, 0),
                new ParameterSet(new[] { 0.9d, 1.6d, 1.9d, 0.12d }, 0),
                new ParameterSet(new[] { 1.1d, 1.4d, 2.1d, 0.08d }, 0),
            };
            return new EnsembleFunction(template, sets);
        }

        /// <summary>
        /// Test index sampling: each draw returns a fresh clone configured with the indexed
        /// parameter set, range violations throw, and clones are independent of each other and
        /// of later samples.
        /// </summary>
        [TestMethod]
        public void Test_Sample_ByIndex()
        {
            var ensemble = BuildEnsemble();
            Assert.AreEqual(3, ensemble.Count);

            var first = (SegmentedPowerFunction)ensemble.SampleAt(0);
            var second = (SegmentedPowerFunction)ensemble.SampleAt(1);
            Assert.AreEqual(1.0d, first.GetBreakpoint(1), 0);
            Assert.AreEqual(0.9d, second.GetBreakpoint(1), 0);
            Assert.AreEqual(0.12d, second.Sigma, 0);
            Assert.AreNotSame(first, second);

            // Mutating one clone never touches another draw of the same index.
            second.SetParameters(new[] { 5d, 5d, 5d, 5d });
            var secondAgain = (SegmentedPowerFunction)ensemble.SampleAt(1);
            Assert.AreEqual(0.9d, secondAgain.GetBreakpoint(1), 0, "Clones must be independent.");

            Assert.Throws<ArgumentOutOfRangeException>(() => ensemble.SampleAt(-1));
            Assert.Throws<ArgumentOutOfRangeException>(() => ensemble.SampleAt(3));
        }

        /// <summary>
        /// Test percentile sampling: u maps onto the index ladder as min(⌊u·N⌋, N − 1), and
        /// out-of-range percentiles throw.
        /// </summary>
        [TestMethod]
        public void Test_Sample_ByPercentile()
        {
            var ensemble = BuildEnsemble();
            Assert.AreEqual(1.0d, ((SegmentedPowerFunction)ensemble.Sample(0.0)).GetBreakpoint(1), 0);
            Assert.AreEqual(0.9d, ((SegmentedPowerFunction)ensemble.Sample(0.5)).GetBreakpoint(1), 0, "u = 0.5 with N = 3 selects index 1.");
            Assert.AreEqual(1.1d, ((SegmentedPowerFunction)ensemble.Sample(1.0)).GetBreakpoint(1), 0, "u = 1 clamps to the last index.");
            Assert.Throws<ArgumentOutOfRangeException>(() => ensemble.Sample(-0.1));
            Assert.Throws<ArgumentOutOfRangeException>(() => ensemble.Sample(1.1));
        }

        /// <summary>
        /// Verifies that an integer percentile literal uses percentile sampling instead of
        /// silently binding to an index overload.
        /// </summary>
        [TestMethod]
        public void Test_Sample_IntegerPercentileLiteral_SelectsLastDraw()
        {
            var ensemble = BuildEnsemble();

            var last = (SegmentedPowerFunction)ensemble.Sample(1);

            Assert.AreEqual(1.1d, last.GetBreakpoint(1), 0d);
        }

        /// <summary>
        /// Test the thread-safety contract: concurrent sampling shares no mutable state, so
        /// every parallel draw evaluates exactly its own parameter set.
        /// </summary>
        [TestMethod]
        public void Test_Parallel_NoSharedMutation()
        {
            var ensemble = BuildEnsemble();
            double[] expected = { 1.0d, 0.9d, 1.1d };
            var failures = 0;
            Parallel.For(0, 3000, i =>
            {
                int index = i % 3;
                var clone = (SegmentedPowerFunction)ensemble.SampleAt(index);
                if (Math.Abs(clone.GetBreakpoint(1) - expected[index]) > 0d)
                    System.Threading.Interlocked.Increment(ref failures);
            });
            Assert.AreEqual(0, failures, "Every parallel draw must carry exactly its own parameter set.");
        }

        /// <summary>
        /// Test the construction guards: null arguments, empty posteriors, length-mismatched
        /// parameter sets, and non-library templates.
        /// </summary>
        [TestMethod]
        public void Test_Construction_Guards()
        {
            var template = new LinearFunction(0, 1, 1);
            var goodSet = new[] { new ParameterSet(new[] { 0d, 1d, 1d }, 0) };
            Assert.Throws<ArgumentNullException>(() => new EnsembleFunction(null, goodSet));
            Assert.Throws<ArgumentNullException>(() => new EnsembleFunction(template, null));
            Assert.Throws<ArgumentException>(() => new EnsembleFunction(template, Array.Empty<ParameterSet>()));
            Assert.Throws<ArgumentException>(() => new EnsembleFunction(template, new[] { new ParameterSet(new[] { 0d, 1d }, 0) }));
        }

        /// <summary>
        /// Verifies that a tabular template is rejected at construction because its parameter
        /// vector cannot be applied to sampled clones.
        /// </summary>
        [TestMethod]
        public void Test_Construction_RejectsTabularTemplate()
        {
            var data = new UncertainOrderedPairedData(
                new[]
                {
                    new UncertainOrdinate(0d, new Deterministic(0d)),
                    new UncertainOrdinate(1d, new Deterministic(1d)),
                },
                true, SortOrder.Ascending, true, SortOrder.Ascending,
                UnivariateDistributionType.Deterministic);
            var template = new TabularFunction(data);
            var sets = new[] { new ParameterSet(new[] { 0d }, 0d) };

            Assert.Throws<NotSupportedException>(() => new EnsembleFunction(template, sets));
        }

        /// <summary>
        /// Test the XElement round-trip: the template and every posterior draw restore, and
        /// restored samples evaluate identically.
        /// </summary>
        [TestMethod]
        public void Test_Serialization_RoundTrip()
        {
            var original = BuildEnsemble();
            var restored = EnsembleFunction.FromXElement(original.ToXElement());

            Assert.AreEqual(original.Count, restored.Count);
            for (int i = 0; i < original.Count; i++)
            {
                var a = original.SampleAt(i);
                var b = restored.SampleAt(i);
                a.ConfidenceLevel = 0.75;
                b.ConfidenceLevel = 0.75;
                Assert.AreEqual(a.Function(5d), b.Function(5d), 1E-12, $"Draw {i} must evaluate identically after the round-trip.");
            }

            Assert.Throws<ArgumentException>(() => EnsembleFunction.FromXElement(new XElement(nameof(EnsembleFunction))));
        }

        /// <summary>
        /// Test that construction deep-copies the template and every parameter set, so
        /// caller-owned arrays and exposed sets cannot mutate later draws, and that serialized
        /// sets with malformed values or fitness are rejected.
        /// </summary>
        [TestMethod]
        public void Test_OwnsDeepCopies_AndValidatesXmlSets()
        {
            var values = new[] { 1d, 0.5d, 2d, 0.1d };
            var ensemble = new EnsembleFunction(
                new SegmentedPowerFunction(values),
                new[] { new ParameterSet(values, 1d, 0.5d) });

            values[0] = 99d;
            Assert.AreEqual(1d, ((SegmentedPowerFunction)ensemble.SampleAt(0)).GetBreakpoint(1), 0d);

            ParameterSet exposed = ensemble.ParameterSets[0];
            exposed.Values[0] = 88d;
            Assert.AreEqual(1d, ((SegmentedPowerFunction)ensemble.SampleAt(0)).GetBreakpoint(1), 0d);
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

        /// <summary>Builds the composite-template ensemble used by the template-pool regression tests.</summary>
        private static EnsembleFunction BuildPoolFixture()
        {
            var children = new List<IUnivariateFunction>
            {
                new LinearFunction(1.0, 0.5, 0.2) { IsDeterministic = false },
                new PowerFunction(2.0, 1.5, 0.0, 0.1) { IsDeterministic = false },
                new LinearFunction(-0.5, 1.25, 0.3) { IsDeterministic = false }
            };
            var template = new CompositeFunction(children, new[] { 0.2, 0.5, 0.3 });
            var sets = new[]
            {
                new ParameterSet(new[] { 0.2, 0.5, 0.3 }, 0),
                new ParameterSet(new[] { 0.6, 0.1, 0.3 }, 0),
                new ParameterSet(new[] { 0.0, 0.0, 1.0 }, 0),
            };
            return new EnsembleFunction(template, sets);
        }

        /// <summary>Reference values of each fixture draw at confidence level 0.3.</summary>
        private static readonly double[][] PoolFixtureExpected =
        {
            new[] { 0.5548187073790324, 1.4057391387195495, 4.920228816135511, 9.473117576320714 },
            new[] { 0.744474047137527, 1.2046581334056303, 2.7775560688888223, 4.558133820925863 },
            new[] { -0.03232015381241218, 0.5926798461875878, 2.467679846187588, 4.342679846187588 },
        };

        /// <summary>
        /// Test that repeated draws from the pooled parsed template reproduce fixed reference
        /// values exactly.
        /// </summary>
        [TestMethod]
        public void Test_TemplatePool_MatchesReferenceValuesExactly()
        {
            var ensemble = BuildPoolFixture();
            var xs = new[] { 0.5, 1.0, 2.5, 4.0 };
            for (int repeat = 0; repeat < 3; repeat++)
            {
                for (int i = 0; i < ensemble.Count; i++)
                {
                    var draw = ensemble.SampleAt(i);
                    draw.ConfidenceLevel = 0.3;
                    for (int j = 0; j < xs.Length; j++)
                        Assert.AreEqual(PoolFixtureExpected[i][j], draw.Function(xs[j]), 0d);
                }
            }
        }

        /// <summary>
        /// Test that concurrent draws through the shared template pool return independent clones
        /// with exactly the sequential values.
        /// </summary>
        [TestMethod]
        public void Test_TemplatePool_ParallelDrawsMatchSequential()
        {
            var ensemble = BuildPoolFixture();
            var xs = new[] { 0.5, 1.0, 2.5, 4.0 };
            int mismatches = 0;
            Parallel.For(0, 3000, k =>
            {
                int i = k % ensemble.Count;
                var draw = ensemble.Sample((i + 0.5) / ensemble.Count);
                draw.ConfidenceLevel = 0.3;
                for (int j = 0; j < xs.Length; j++)
                {
                    if (draw.Function(xs[j]) != PoolFixtureExpected[i][j])
                        Interlocked.Increment(ref mismatches);
                }
            });
            Assert.AreEqual(0, mismatches);
        }

    }
}
