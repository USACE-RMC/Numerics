using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Functions;

namespace Functions
{
    /// <summary>
    /// Unit tests for the BaRatin addition-mode <see cref="SegmentedPowerFunction"/>: parity
    /// with the RMC-BestFit rating-curve prediction, the single-segment degeneracy to
    /// <see cref="PowerFunction"/>, the numeric inverse, the parameter-vector contract, and the
    /// serialization round-trip.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_SegmentedPowerFunction
    {
        /// <summary>
        /// Test parity with the BestFit rating-curve prediction. The reference constants were
        /// evaluated independently from the BaRatin addition-mode closed form
        /// Q(h) = Σₖ 10^(log₁₀αₖ)(h − hₖ)^βₖ·𝟙{h &gt; hₖ}, the exact body of
        /// <c>RMC.BestFit RatingCurve.Predict</c> (Python evaluation, 2026-07-25).
        /// </summary>
        [TestMethod]
        public void Test_BestFitParity()
        {
            // 1 segment: [h₁, log₁₀α₁, β₁, σ] = [1, 1.5, 2, 0.1]
            var oneSegment = new SegmentedPowerFunction(new[] { 1d, 1.5d, 2d, 0.1d }) { IsDeterministic = true };
            Assert.AreEqual(0d, oneSegment.Function(0.5), 0d, "Discharge is zero at and below the cease-to-flow stage h₁.");
            Assert.AreEqual(505.9644256269407, oneSegment.Function(5), 1E-10);

            // 2 segments: [1, 1.5, 2, 3, 1.2, 1.5, 0.1] — below and above the second activation.
            var twoSegments = new SegmentedPowerFunction(new[] { 1d, 1.5d, 2d, 3d, 1.2d, 1.5d, 0.1d }) { IsDeterministic = true };
            Assert.AreEqual(71.15124735378853, twoSegments.Function(2.5), 1E-10, "Only the first control is active below h₂.");
            Assert.AreEqual(550.7919745807667, twoSegments.Function(5), 1E-10, "Both controls add above h₂.");

            // 3 segments: [1, 1.5, 2, 3, 1.2, 1.5, 6, 0.8, 1.1, 0.15]
            var threeSegments = new SegmentedPowerFunction(new[] { 1d, 1.5d, 2d, 3d, 1.2d, 1.5d, 6d, 0.8d, 1.1d, 0.15d }) { IsDeterministic = true };
            Assert.AreEqual(300.45392133976526, threeSegments.Function(4), 1E-10);
            Assert.AreEqual(2277.915263028223, threeSegments.Function(9), 1E-9);

            // The log₁₀-space residual at a fixed confidence level multiplies the curve by
            // 10^(z·σ): 550.79197… · 10^(0.674489750196…·0.1) = 643.3341101314571.
            twoSegments.IsDeterministic = false;
            twoSegments.ConfidenceLevel = 0.75;
            Assert.AreEqual(643.3341101314571, twoSegments.Function(5), 1E-9);
        }

        /// <summary>
        /// Test that one segment degenerates to the plain power law: α = 10^(log₁₀α₁),
        /// β = β₁, ξ = h₁.
        /// </summary>
        [TestMethod]
        public void Test_OneSegment_DegeneratesToPowerFunction()
        {
            var segmented = new SegmentedPowerFunction(new[] { 2d, 0.7d, 1.8d, 0.1d }) { IsDeterministic = true };
            var power = new PowerFunction(Math.Pow(10d, 0.7d), 1.8d, 2d) { IsDeterministic = true };
            foreach (double stage in new[] { 2.5, 4.0, 10.0, 55.0 })
            {
                Assert.AreEqual(power.Function(stage), segmented.Function(stage), 1E-10 * power.Function(stage));
            }
        }

        /// <summary>
        /// Test the numeric inverse: round-trips within and across segments, at a stochastic
        /// confidence level, and the zero/edge policies.
        /// </summary>
        [TestMethod]
        public void Test_InverseFunction_RoundTrip()
        {
            var function = new SegmentedPowerFunction(new[] { 1d, 1.5d, 2d, 3d, 1.2d, 1.5d, 0.1d }) { IsDeterministic = true };
            foreach (double stage in new[] { 1.5, 2.5, 3.5, 5.0, 12.0 })
            {
                double q = function.Function(stage);
                Assert.AreEqual(stage, function.InverseFunction(q), 1E-6, $"Deterministic inverse at stage {stage}.");
            }

            // A stochastic confidence level folds out of the inverse exactly.
            function.IsDeterministic = false;
            function.ConfidenceLevel = 0.9;
            double stochastic = function.Function(5.0);
            Assert.AreEqual(5.0, function.InverseFunction(stochastic), 1E-6);

            // Non-positive discharge maps to the cease-to-flow stage.
            Assert.AreEqual(1d, function.InverseFunction(0d), 0d);
            Assert.AreEqual(1d, function.InverseFunction(-5d), 0d);
        }

        /// <summary>
        /// Test the parameter-vector contract: SetParameters applies a BestFit-layout posterior
        /// draw directly, and the validation guards reject malformed vectors.
        /// </summary>
        [TestMethod]
        public void Test_ParameterVector_Contract()
        {
            var function = new SegmentedPowerFunction(new[] { 1d, 1.5d, 2d, 3d, 1.2d, 1.5d, 0.1d });
            Assert.AreEqual(2, function.NumberOfSegments);
            Assert.AreEqual(7, function.NumberOfParameters);
            Assert.AreEqual(1d, function.GetBreakpoint(1), 0d);
            Assert.AreEqual(1.2d, function.GetLog10Alpha(2), 0d);
            Assert.AreEqual(1.5d, function.GetBeta(2), 0d);
            Assert.AreEqual(0.1d, function.Sigma, 0d);
            Assert.AreEqual(1d, function.Minimum, 0d, "Minimum derives from h₁.");

            // A posterior draw applies in place (same layout, same length).
            function.SetParameters(new[] { 0.8d, 1.6d, 1.9d, 2.9d, 1.1d, 1.4d, 0.12d });
            Assert.IsTrue(function.ParametersValid);
            Assert.AreEqual(0.8d, function.GetBreakpoint(1), 0d);
            Assert.AreEqual(0.12d, function.Sigma, 0d);

            // Guards: wrong length, unordered breakpoints, negative exponents, non-positive σ.
            Assert.IsNotNull(function.ValidateParameters(new[] { 1d, 1.5d, 2d }, false));
            Assert.IsNotNull(function.ValidateParameters(new[] { 3d, 1.5d, 2d, 1d, 1.2d, 1.5d, 0.1d }, false));
            Assert.IsNotNull(function.ValidateParameters(new[] { 1d, 1.5d, -2d, 3d, 1.2d, 1.5d, 0.1d }, false));
            Assert.IsNotNull(function.ValidateParameters(new[] { 1d, 1.5d, 2d, 3d, 1.2d, 1.5d, 0d }, false));
            Assert.Throws<ArgumentOutOfRangeException>(() => new SegmentedPowerFunction(0));
            Assert.Throws<ArgumentException>(() => new SegmentedPowerFunction(new[] { 1d, 2d, 3d, 4d, 5d }));
        }

        /// <summary>
        /// Test the XElement round-trip through the factory, and the factory dispatch surfaces.
        /// </summary>
        [TestMethod]
        public void Test_Serialization_RoundTrip()
        {
            var original = new SegmentedPowerFunction(new[] { 1d, 1.5d, 2d, 3d, 1.2d, 1.5d, 0.1d }) { Maximum = 500d };
            var element = original.ToXElement();
            Assert.IsNull(element.Attribute("ConfidenceLevel"), "ConfidenceLevel is runtime state and must not serialize.");
            Assert.IsNull(element.Attribute("Minimum"), "Minimum derives from h₁ and must not serialize.");

            var restored = (SegmentedPowerFunction)UnivariateFunctionFactory.CreateFromXElement(element);
            Assert.AreEqual(original.NumberOfSegments, restored.NumberOfSegments);
            Assert.AreEqual(original.IsDeterministic, restored.IsDeterministic);
            Assert.AreEqual(original.Maximum, restored.Maximum, 0d);
            for (int k = 1; k <= original.NumberOfSegments; k++)
            {
                Assert.AreEqual(original.GetBreakpoint(k), restored.GetBreakpoint(k), 0d);
                Assert.AreEqual(original.GetLog10Alpha(k), restored.GetLog10Alpha(k), 0d);
                Assert.AreEqual(original.GetBeta(k), restored.GetBeta(k), 0d);
            }
            Assert.AreEqual(original.Sigma, restored.Sigma, 0d);

            original.ConfidenceLevel = 0.75;
            restored.ConfidenceLevel = 0.75;
            Assert.AreEqual(original.Function(5), restored.Function(5), 1E-12);

            Assert.AreEqual(UnivariateFunctionType.SegmentedPower, UnivariateFunctionFactory.GetFunctionType(original));
            Assert.IsInstanceOfType(UnivariateFunctionFactory.CreateFunction(UnivariateFunctionType.SegmentedPower), typeof(SegmentedPowerFunction));
        }
    }
}
