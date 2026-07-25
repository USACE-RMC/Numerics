using System;
using System.Threading.Tasks;
using Microsoft.VisualStudio.TestTools.UnitTesting;
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

            var first = (SegmentedPowerFunction)ensemble.Sample(0);
            var second = (SegmentedPowerFunction)ensemble.Sample(1);
            Assert.AreEqual(1.0d, first.GetBreakpoint(1), 0);
            Assert.AreEqual(0.9d, second.GetBreakpoint(1), 0);
            Assert.AreEqual(0.12d, second.Sigma, 0);
            Assert.AreNotSame(first, second);

            // Mutating one clone never touches another draw of the same index.
            second.SetParameters(new[] { 5d, 5d, 5d, 5d });
            var secondAgain = (SegmentedPowerFunction)ensemble.Sample(1);
            Assert.AreEqual(0.9d, secondAgain.GetBreakpoint(1), 0, "Clones must be independent.");

            Assert.Throws<ArgumentOutOfRangeException>(() => ensemble.Sample(-1));
            Assert.Throws<ArgumentOutOfRangeException>(() => ensemble.Sample(3));
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
                var clone = (SegmentedPowerFunction)ensemble.Sample(index);
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
                var a = original.Sample(i);
                var b = restored.Sample(i);
                a.ConfidenceLevel = 0.75;
                b.ConfidenceLevel = 0.75;
                Assert.AreEqual(a.Function(5d), b.Function(5d), 1E-12, $"Draw {i} must evaluate identically after the round-trip.");
            }

            Assert.Throws<ArgumentException>(() => EnsembleFunction.FromXElement(new System.Xml.Linq.XElement(nameof(EnsembleFunction))));
        }
    }
}
