using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Functions;

namespace Functions
{
    /// <summary>
    /// Unit tests for <see cref="CompositeFunction"/>: the weighted-average and mixture modes,
    /// the single-uniform mixture composition, the deterministic semantics, weight validation,
    /// the numeric inverse, and the serialization round-trip (including nesting).
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
    }
}
