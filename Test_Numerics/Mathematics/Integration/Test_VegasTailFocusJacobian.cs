using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the VEGAS probability-space power transform and its Jacobian correction.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_VegasTailFocusJacobian
    {
        /// <summary>
        /// Runs Vegas over the heavy-upper-tail integrand f(p) = 21·(1 − p)^20 on [0, 1]
        /// (∫ = 1 exactly) at a given tail-focus parameter, accumulating the weights handed to
        /// the integrand.
        /// </summary>
        /// <param name="gamma">The power-transform tail-focus parameter γ.</param>
        /// <returns>The integral estimate, the total weight handed out, and the call count.</returns>
        private static (double Result, double WeightSum, long Calls) Run(double gamma)
        {
            double weightSum = 0d;
            long calls = 0;
            var vegas = new Vegas((x, w) =>
            {
                weightSum += w;
                calls++;
                return 21d * Math.Pow(1d - x[0], 20d);
            }, 1, new[] { 0d }, new[] { 1d })
            {
                UseSobolSequence = false,
                Random = new MersenneTwister(12345),
                TailFocusParameter = gamma,
                IndependentEvaluations = 10,
                FunctionCalls = 10000,
            };
            vegas.Integrate();
            Assert.AreEqual(IntegrationStatus.Success, vegas.Status, $"γ = {gamma}: the integration must succeed.");
            return (vegas.Result, weightSum, calls);
        }

        /// <summary>
        /// Test that the heavy-tail integral is unbiased at every tail-focus setting: the
        /// power transform changes where samples land, and the folded Jacobian must exactly
        /// compensate, so γ ∈ {1, 4, 10} all reproduce ∫ 21(1 − p)^20 dp = 1.
        /// </summary>
        [TestMethod]
        public void Test_HeavyTailIntegral_SameValueAtEveryGamma()
        {
            foreach (double gamma in new[] { 1d, 4d, 10d })
            {
                var run = Run(gamma);
                Assert.AreEqual(1d, run.Result, 0.03d, $"γ = {gamma}: the tail-focused integral must stay unbiased.");
            }
        }

        /// <summary>
        /// Test that the weight handed to the integrand carries the domain measure at every
        /// tail-focus setting: the average weight per evaluation batch must reproduce the
        /// domain volume (a missing Jacobian would distort this by the measure change, an
        /// order-of-magnitude error at γ = 10).
        /// </summary>
        [TestMethod]
        public void Test_WeightSum_EqualsDomainVolumeAtEveryGamma()
        {
            foreach (double gamma in new[] { 1d, 4d, 10d })
            {
                var run = Run(gamma);
                double volumePerBatch = run.WeightSum / (run.Calls / 10000d);
                Assert.AreEqual(1d, volumePerBatch, 0.1d,
                    $"γ = {gamma}: the weights must sum to the domain volume per evaluation batch (Jacobian folded into the weight).");
            }
        }

        /// <summary>
        /// Test that the tail-focus parameter and the rare-event configuration reject
        /// non-finite or non-positive values, leaving the configured state unchanged.
        /// </summary>
        [TestMethod]
        public void Test_TailFocus_RejectsInvalidParametersAtomically()
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
