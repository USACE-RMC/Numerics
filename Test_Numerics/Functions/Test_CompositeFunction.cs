using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Functions;

namespace Functions
{
    /// <summary>
    /// Unit tests for <see cref="CompositeFunction"/>: the weighted-average and mixture modes,
    /// the single-uniform mixture composition, the deterministic semantics, weight validation,
    /// the numeric inverse, exception-safe child state restoration, stateful-child
    /// synchronization, and the serialization round-trip (including nesting).
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_CompositeFunction
    {
        /// <summary>
        /// Test the weighted-average mode over deterministic children:
        /// 0.25·(2x) + 0.75·(4x + 10) = 3.5x + 7.5 exactly.
        /// </summary>
        [TestMethod]
        public void Test_WeightedAverage_Deterministic()
        {
            var composite = new CompositeFunction(
                new IUnivariateFunction[] { new LinearFunction(0, 2), new LinearFunction(10, 4) },
                new[] { 0.25d, 0.75d });

            foreach (double x in new[] { -3d, 0d, 5d, 42d })
            {
                Assert.AreEqual(3.5d * x + 7.5d, composite.Function(x), 1E-12);
            }

            // The numeric inverse round-trips the monotone composed curve.
            Assert.AreEqual(5d, composite.InverseFunction(composite.Function(5d)), 1E-8);
        }

        /// <summary>
        /// Test that the parent mean sentinel combines independently configured child
        /// realizations without replacing their confidence levels.
        /// </summary>
        [TestMethod]
        public void Test_WeightedAverage_MeanSentinelUsesConfiguredChildLevels()
        {
            var first = new LinearFunction(0d, 1d, 2d) { ConfidenceLevel = 0.2d };
            var second = new LinearFunction(10d, 1d, 3d) { ConfidenceLevel = 0.8d };
            var composite = new CompositeFunction(
                new IUnivariateFunction[] { first, second }, new[] { 0.25d, 0.75d });
            double expected = 0.25d * first.Function(5d) + 0.75d * second.Function(5d);

            Assert.AreEqual(expected, composite.Function(5d), 1E-12);
            Assert.AreEqual(5d, composite.InverseFunction(expected), 1E-8);
            Assert.AreEqual(0.2d, first.ConfidenceLevel, 0d);
            Assert.AreEqual(0.8d, second.ConfidenceLevel, 0d);
        }

        /// <summary>
        /// Test the mixture composition: one uniform selects the child by cumulative weight and
        /// re-scales the remainder as the child's own draw; the mean convention (a confidence
        /// level outside [0, 1]) is the weighted average of the child means.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_SingleUniformComposition()
        {
            var composite = new CompositeFunction(
                new IUnivariateFunction[] { new LinearFunction(0, 1), new LinearFunction(100, 1) },
                new[] { 0.4d, 0.6d })
            {
                Mode = CompositeFunctionMode.Mixture,
            };

            // Below the first cumulative weight the first child is selected; above it, the second.
            composite.ConfidenceLevel = 0.2;
            Assert.AreEqual(5d, composite.Function(5d), 1E-12);
            composite.ConfidenceLevel = 0.4;
            Assert.AreEqual(5d, composite.Function(5d), 1E-12, "The bracket boundary belongs to the lower child.");
            composite.ConfidenceLevel = 0.8;
            Assert.AreEqual(105d, composite.Function(5d), 1E-12);

            // The mixture inverse inverts the selected child.
            Assert.AreEqual(5d, composite.InverseFunction(105d), 1E-12);

            // The re-scaled remainder drives the selected child's own draw: with equal-weight
            // uncertain children, u = 0.75 selects child 2 at remainder 0.5 — its median.
            var uncertain = new CompositeFunction(
                new IUnivariateFunction[] { new LinearFunction(0, 1, 2), new LinearFunction(0, 1, 2) },
                new[] { 0.5d, 0.5d })
            {
                Mode = CompositeFunctionMode.Mixture,
                ConfidenceLevel = 0.75,
            };
            Assert.AreEqual(5d, uncertain.Function(5d), 1E-12, "Remainder 0.5 evaluates the selected child at its median.");

            // The mean convention: the weighted average of the child means.
            composite.ConfidenceLevel = -1;
            Assert.AreEqual(65d, composite.Function(5d), 1E-12);
        }

        /// <summary>
        /// Test the deterministic semantics: averages of deterministic children are
        /// deterministic; a multi-branch mixture is not, even over deterministic children; a
        /// single-reachable-branch mixture is.
        /// </summary>
        [TestMethod]
        public void Test_IsDeterministic_Semantics()
        {
            var children = new IUnivariateFunction[] { new LinearFunction(0, 1), new LinearFunction(100, 1) };
            var average = new CompositeFunction(children, new[] { 0.5d, 0.5d });
            Assert.IsTrue(average.IsDeterministic);

            var mixture = new CompositeFunction(children, new[] { 0.5d, 0.5d }) { Mode = CompositeFunctionMode.Mixture };
            Assert.IsFalse(mixture.IsDeterministic, "A multi-branch mixture varies with the draw.");

            var degenerate = new CompositeFunction(children, new[] { 1d, 0d }) { Mode = CompositeFunctionMode.Mixture };
            Assert.IsTrue(degenerate.IsDeterministic, "A single reachable branch cannot vary.");

            var stochasticChild = new CompositeFunction(
                new IUnivariateFunction[] { new LinearFunction(0, 1, 2) }, new[] { 1d });
            Assert.IsFalse(stochasticChild.IsDeterministic);
        }

        /// <summary>
        /// Test the weight validation guards: count mismatch, negatives, and sums away from one.
        /// </summary>
        [TestMethod]
        public void Test_Weights_Validation()
        {
            var children = new IUnivariateFunction[] { new LinearFunction(0, 1), new LinearFunction(100, 1) };
            var composite = new CompositeFunction(children);
            Assert.IsTrue(composite.ParametersValid, "Equal default weights are valid.");

            Assert.IsNotNull(composite.ValidateParameters(new[] { 0.5d }, false));
            Assert.IsNotNull(composite.ValidateParameters(new[] { -0.1d, 1.1d }, false));
            Assert.IsNotNull(composite.ValidateParameters(new[] { 0.3d, 0.3d }, false));
            Assert.IsNull(composite.ValidateParameters(new[] { 0.3d, 0.7d }, false));

            Assert.Throws<ArgumentException>(() => new CompositeFunction(children, new[] { 0.5d }));
            Assert.Throws<ArgumentException>(() => new CompositeFunction(Array.Empty<IUnivariateFunction>()));
        }

        /// <summary>
        /// Test the XElement round-trip through the factory, including a nested composite child
        /// and the factory dispatch surfaces.
        /// </summary>
        [TestMethod]
        public void Test_Serialization_RoundTrip()
        {
            var inner = new CompositeFunction(
                new IUnivariateFunction[] { new LinearFunction(0, 2), new PowerFunction(5, 2, 0) },
                new[] { 0.3d, 0.7d });
            var original = new CompositeFunction(
                new IUnivariateFunction[] { inner, new LinearFunction(10, 4) },
                new[] { 0.6d, 0.4d })
            {
                Mode = CompositeFunctionMode.Mixture,
            };

            var element = original.ToXElement();
            Assert.IsNull(element.Attribute("ConfidenceLevel"), "ConfidenceLevel is runtime state and must not serialize.");

            var restored = (CompositeFunction)UnivariateFunctionFactory.CreateFromXElement(element);
            Assert.AreEqual(original.Mode, restored.Mode);
            Assert.HasCount(original.Functions.Count, restored.Functions);
            Assert.AreEqual(original.Weights[0], restored.Weights[0], 0d);
            Assert.IsInstanceOfType(restored.Functions[0], typeof(CompositeFunction));
            Assert.IsInstanceOfType(restored.Functions[1], typeof(LinearFunction));

            original.ConfidenceLevel = 0.3;
            restored.ConfidenceLevel = 0.3;
            Assert.AreEqual(original.Function(4d), restored.Function(4d), 1E-12);

            Assert.AreEqual(UnivariateFunctionType.Composite, UnivariateFunctionFactory.GetFunctionType(original));
            Assert.Throws<NotSupportedException>(() => UnivariateFunctionFactory.CreateFunction(UnivariateFunctionType.Composite));
        }

        /// <summary>
        /// Test that evaluation restores a child's configured confidence level when the child
        /// throws, on both the forward and inverse paths.
        /// </summary>
        [TestMethod]
        public void Test_ChildStateRestored_AfterExceptions()
        {
            var child = new ConfidenceProbeFunction { ConfidenceLevel = 0.35d, ThrowOnEvaluation = true };
            var composite = new CompositeFunction(new IUnivariateFunction[] { child }) { ConfidenceLevel = 0.8d };

            Assert.Throws<InvalidOperationException>(() => composite.Function(1d));
            Assert.AreEqual(0.35d, child.ConfidenceLevel, 0d);
            Assert.Throws<InvalidOperationException>(() => composite.InverseFunction(1d));
            Assert.AreEqual(0.35d, child.ConfidenceLevel, 0d);
        }

        /// <summary>
        /// Test that concurrent evaluation over a shared stateful child never interleaves
        /// confidence assignments, and that a deterministic child evaluates through the
        /// lock-free path without any confidence mutation.
        /// </summary>
        [TestMethod]
        public void Test_SynchronizesOnlyStatefulChildren()
        {
            var child = new ConfidenceProbeFunction { ConfidenceLevel = 0.5d };
            var lower = new CompositeFunction(new IUnivariateFunction[] { child }) { ConfidenceLevel = 0.2d };
            var upper = new CompositeFunction(new IUnivariateFunction[] { child }) { ConfidenceLevel = 0.8d };
            var configured = new CompositeFunction(new IUnivariateFunction[] { child });
            int failures = 0;

            Parallel.For(0, 600, i =>
            {
                int branch = i % 3;
                double expected = branch == 0 ? 0.2d : branch == 1 ? 0.8d : 0.5d;
                double actual = branch == 0 ? lower.Function(0d)
                    : branch == 1 ? upper.Function(0d)
                    : configured.Function(0d);
                if (actual != expected) Interlocked.Increment(ref failures);
            });

            Assert.AreEqual(0, failures);
            Assert.AreEqual(0.5d, child.ConfidenceLevel, 0d);

            child.IsDeterministic = true;
            child.ThrowOnConfidenceAssignment = true;
            Assert.AreEqual(0.5d, lower.Function(0d), 0d,
                "Deterministic child evaluation must not mutate confidence state.");
        }

        /// <summary>
        /// Test that a rejected weight update leaves the prior weights intact, and that
        /// undefined mode values are rejected at the setter and on deserialization.
        /// </summary>
        [TestMethod]
        public void Test_RejectsInvalidStateAtomically()
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

        /// <summary>
        /// Test function that exposes confidence-state changes and controlled evaluation failures.
        /// </summary>
        private sealed class ConfidenceProbeFunction : IUnivariateFunction
        {
            /// <summary>The configured confidence level.</summary>
            private double _confidenceLevel = -1d;

            /// <summary>Gets or sets whether evaluations throw a controlled exception.</summary>
            public bool ThrowOnEvaluation { get; set; }

            /// <summary>Gets or sets whether confidence-level assignments throw a controlled exception.</summary>
            public bool ThrowOnConfidenceAssignment { get; set; }

            /// <inheritdoc/>
            public int NumberOfParameters => 0;

            /// <inheritdoc/>
            public bool ParametersValid => true;

            /// <inheritdoc/>
            public double Minimum { get; set; } = double.MinValue;

            /// <inheritdoc/>
            public double Maximum { get; set; } = double.MaxValue;

            /// <inheritdoc/>
            public double[] MinimumOfParameters => Array.Empty<double>();

            /// <inheritdoc/>
            public double[] MaximumOfParameters => Array.Empty<double>();

            /// <inheritdoc/>
            public bool IsDeterministic { get; set; }

            /// <inheritdoc/>
            public double ConfidenceLevel
            {
                get { return _confidenceLevel; }
                set
                {
                    if (ThrowOnConfidenceAssignment) throw new InvalidOperationException("Confidence assignment was not expected.");
                    _confidenceLevel = value;
                }
            }

            /// <summary>
            /// Validates that no parameters are supplied to this parameterless test function.
            /// </summary>
            /// <param name="parameters">The parameter collection, which must be empty.</param>
            /// <exception cref="ArgumentNullException">Thrown when <paramref name="parameters"/> is null.</exception>
            /// <exception cref="ArgumentException">Thrown when <paramref name="parameters"/> is not empty.</exception>
            public void SetParameters(IList<double> parameters)
            {
                if (parameters == null) throw new ArgumentNullException(nameof(parameters));
                if (parameters.Count != 0) throw new ArgumentException("This test function has no parameters.", nameof(parameters));
            }

            /// <summary>
            /// Reports that the parameterless test function has no range-validation error.
            /// </summary>
            /// <param name="parameters">The parameter collection.</param>
            /// <param name="throwException">Ignored because this test function has no parameters.</param>
            /// <returns><see langword="null"/>.</returns>
            public ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
            {
                return null;
            }

            /// <summary>
            /// Returns the configured confidence level after checking for concurrent state changes.
            /// </summary>
            /// <param name="x">The unused evaluation point.</param>
            /// <returns>The configured confidence level.</returns>
            /// <exception cref="InvalidOperationException">Thrown when controlled failure is enabled or confidence state changes during evaluation.</exception>
            public double Function(double x)
            {
                if (ThrowOnEvaluation) throw new InvalidOperationException("Test evaluation failure.");
                double first = ConfidenceLevel;
                Thread.SpinWait(10000);
                if (first != ConfidenceLevel) throw new InvalidOperationException("Confidence state changed during evaluation.");
                return ConfidenceLevel;
            }

            /// <summary>
            /// Returns the configured confidence level as the inverse result.
            /// </summary>
            /// <param name="y">The unused value to invert.</param>
            /// <returns>The configured confidence level.</returns>
            /// <exception cref="InvalidOperationException">Thrown when controlled failure is enabled.</exception>
            public double InverseFunction(double y)
            {
                if (ThrowOnEvaluation) throw new InvalidOperationException("Test inverse failure.");
                return ConfidenceLevel;
            }
        }
    }
}
