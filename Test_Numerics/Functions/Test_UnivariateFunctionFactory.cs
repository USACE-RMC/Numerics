using System;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Distributions;
using Numerics.Functions;

namespace Functions
{
    /// <summary>
    /// Unit tests for the univariate function serialization surface: the
    /// <see cref="UnivariateFunctionType"/> enum, the <see cref="UnivariateFunctionFactory"/>
    /// dispatch, and the per-class ToXElement/FromXElement round-trips.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_UnivariateFunctionFactory
    {
        /// <summary>
        /// Test that the factory creates default instances by type, refuses the parameterless
        /// tabular creation, and rejects undefined types.
        /// </summary>
        [TestMethod]
        public void Test_CreateFunction_ByType()
        {
            var linear = UnivariateFunctionFactory.CreateFunction(UnivariateFunctionType.Linear);
            Assert.IsInstanceOfType(linear, typeof(LinearFunction));
            var power = UnivariateFunctionFactory.CreateFunction(UnivariateFunctionType.Power);
            Assert.IsInstanceOfType(power, typeof(PowerFunction));

            Assert.Throws<NotSupportedException>(() => UnivariateFunctionFactory.CreateFunction(UnivariateFunctionType.Tabular));
            Assert.Throws<ArgumentOutOfRangeException>(() => UnivariateFunctionFactory.CreateFunction((UnivariateFunctionType)999));
        }

        /// <summary>
        /// Test that the factory maps live instances back to their enum type and rejects
        /// non-library implementations of the interface.
        /// </summary>
        [TestMethod]
        public void Test_GetFunctionType()
        {
            Assert.AreEqual(UnivariateFunctionType.Linear, UnivariateFunctionFactory.GetFunctionType(new LinearFunction()));
            Assert.AreEqual(UnivariateFunctionType.Power, UnivariateFunctionFactory.GetFunctionType(new PowerFunction()));
            var table = new UncertainOrderedPairedData(
                new[] { new UncertainOrdinate(0d, new Deterministic(0d)), new UncertainOrdinate(1d, new Deterministic(1d)) },
                true, SortOrder.Ascending, false, SortOrder.None, UnivariateDistributionType.Deterministic);
            Assert.AreEqual(UnivariateFunctionType.Tabular, UnivariateFunctionFactory.GetFunctionType(new TabularFunction(table)));
            Assert.Throws<ArgumentNullException>(() => UnivariateFunctionFactory.GetFunctionType(null));
        }

        /// <summary>
        /// Test the linear function XElement round-trip: every serialized member restores, the
        /// runtime confidence level is not serialized, and evaluation matches at a fixed
        /// percentile.
        /// </summary>
        [TestMethod]
        public void Test_LinearFunction_RoundTrip()
        {
            var original = new LinearFunction(-2, 5, 3) { Minimum = -100, Maximum = 250 };
            var element = original.ToXElement();
            Assert.IsNull(element.Attribute("ConfidenceLevel"), "ConfidenceLevel is runtime state and must not serialize.");

            var restored = (LinearFunction)UnivariateFunctionFactory.CreateFromXElement(element);
            Assert.AreEqual(original.Alpha, restored.Alpha, 0);
            Assert.AreEqual(original.Beta, restored.Beta, 0);
            Assert.AreEqual(original.Sigma, restored.Sigma, 0);
            Assert.AreEqual(original.IsDeterministic, restored.IsDeterministic);
            Assert.AreEqual(original.Minimum, restored.Minimum, 0);
            Assert.AreEqual(original.Maximum, restored.Maximum, 0);

            original.ConfidenceLevel = 0.75;
            restored.ConfidenceLevel = 0.75;
            Assert.AreEqual(original.Function(6), restored.Function(6), 1E-12);
        }

        /// <summary>
        /// Test the power function XElement round-trip, including the inverse flag and the
        /// derived minimum (ξ is serialized; Minimum is not).
        /// </summary>
        [TestMethod]
        public void Test_PowerFunction_RoundTrip()
        {
            var original = new PowerFunction(5, 2, 1, 0.1) { IsInverse = true, Maximum = 5000 };
            var element = original.ToXElement();
            Assert.IsNull(element.Attribute("ConfidenceLevel"), "ConfidenceLevel is runtime state and must not serialize.");
            Assert.IsNull(element.Attribute("Minimum"), "Minimum derives from Xi and must not serialize.");

            var restored = (PowerFunction)UnivariateFunctionFactory.CreateFromXElement(element);
            Assert.AreEqual(original.Alpha, restored.Alpha, 0);
            Assert.AreEqual(original.Beta, restored.Beta, 0);
            Assert.AreEqual(original.Xi, restored.Xi, 0);
            Assert.AreEqual(original.Sigma, restored.Sigma, 0);
            Assert.AreEqual(original.IsDeterministic, restored.IsDeterministic);
            Assert.AreEqual(original.IsInverse, restored.IsInverse);
            Assert.AreEqual(original.Minimum, restored.Minimum, 0);
            Assert.AreEqual(original.Maximum, restored.Maximum, 0);

            original.ConfidenceLevel = 0.75;
            restored.ConfidenceLevel = 0.75;
            Assert.AreEqual(original.Function(6), restored.Function(6), 1E-12);
        }

        /// <summary>
        /// Test the tabular function XElement round-trip: the embedded uncertain paired data,
        /// the axis transforms, and the negative-Y policy all restore, and a missing table
        /// throws.
        /// </summary>
        [TestMethod]
        public void Test_TabularFunction_RoundTrip()
        {
            var table = new UncertainOrderedPairedData(
                new[] { new UncertainOrdinate(1d, new Normal(10, 2)), new UncertainOrdinate(100d, new Normal(20, 2)) },
                true, SortOrder.Ascending, false, SortOrder.None, UnivariateDistributionType.Normal);
            var original = new TabularFunction(table)
            {
                // Both axes carry positive hazard-scale values in this fixture, so both use the
                // logarithmic transform (Normal-Z is reserved for probability-valued axes).
                XTransform = Transform.Logarithmic,
                YTransform = Transform.Logarithmic,
                AllowNegativeYValues = false,
            };

            var restored = (TabularFunction)UnivariateFunctionFactory.CreateFromXElement(original.ToXElement());
            Assert.AreEqual(original.XTransform, restored.XTransform);
            Assert.AreEqual(original.YTransform, restored.YTransform);
            Assert.AreEqual(original.AllowNegativeYValues, restored.AllowNegativeYValues);
            Assert.AreEqual(original.PairedData.Count, restored.PairedData.Count);
            Assert.AreEqual(original.PairedData[1].X, restored.PairedData[1].X, 0);
            Assert.AreEqual(original.PairedData[1].Y.Mean, restored.PairedData[1].Y.Mean, 0);

            original.ConfidenceLevel = 0.9;
            restored.ConfidenceLevel = 0.9;
            Assert.AreEqual(original.Function(50), restored.Function(50), 1E-12);

            Assert.Throws<ArgumentException>(() => TabularFunction.FromXElement(new XElement(nameof(TabularFunction))));
        }

        /// <summary>
        /// The extrapolation attribute is written only when non-default, so every pre-existing
        /// serialized form stays byte-identical; a non-default value round-trips, and an absent
        /// attribute reads as None.
        /// </summary>
        [TestMethod]
        public void Test_TabularFunction_Extrapolation_ConditionalPresence()
        {
            var table = new UncertainOrderedPairedData(
                new[] { new UncertainOrdinate(1d, new Normal(10, 2)), new UncertainOrdinate(100d, new Normal(20, 2)) },
                true, SortOrder.Ascending, false, SortOrder.None, UnivariateDistributionType.Normal);

            // Default: no attribute, and the element is byte-identical to a pre-property form.
            var original = new TabularFunction(table);
            var defaultXml = original.ToXElement();
            Assert.IsNull(defaultXml.Attribute(nameof(TabularFunction.Extrapolation)));
            string baseline = defaultXml.ToString();
            original.Extrapolation = ExtrapolationSides.Both;
            original.Extrapolation = ExtrapolationSides.None;
            Assert.AreEqual(baseline, original.ToXElement().ToString());

            // Non-default: attribute present by enum name and restored through the factory.
            original.Extrapolation = ExtrapolationSides.Above;
            var xml = original.ToXElement();
            Assert.AreEqual(nameof(ExtrapolationSides.Above), xml.Attribute(nameof(TabularFunction.Extrapolation))?.Value);
            var restored = (TabularFunction)UnivariateFunctionFactory.CreateFromXElement(xml);
            Assert.AreEqual(ExtrapolationSides.Above, restored.Extrapolation);

            // Absent attribute reads as the default.
            var legacy = TabularFunction.FromXElement(defaultXml);
            Assert.AreEqual(ExtrapolationSides.None, legacy.Extrapolation);
        }

        /// <summary>
        /// Test that the factory rejects null and unknown serialized forms.
        /// </summary>
        [TestMethod]
        public void Test_CreateFromXElement_Guards()
        {
            Assert.Throws<ArgumentNullException>(() => UnivariateFunctionFactory.CreateFromXElement(null));
            Assert.Throws<NotSupportedException>(() => UnivariateFunctionFactory.CreateFromXElement(new XElement("BogusFunction")));
        }
    }
}
