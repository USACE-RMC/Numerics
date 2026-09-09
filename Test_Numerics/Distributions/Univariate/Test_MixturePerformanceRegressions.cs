using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Regression coverage for allocation-conscious validation of mutable mixtures.</summary>
    [TestClass]
    public class Test_MixturePerformanceRegressions
    {
        /// <summary>Receives diagnostic allocation output from the focused tests.</summary>
        public TestContext TestContext { get; set; }

#if NET8_0_OR_GREATER
        /// <summary>Repeated density evaluation validates live state without rebuilding flattened candidates.</summary>
        [TestMethod]
        public void LogPDF_RepeatedLiveValidationAvoidsFlattenedCandidateAllocations()
        {
            var mixtures = new[]
            {
                new Mixture(new[] { 0.4d, 0.6d }, new UnivariateDistributionBase[]
                    { new Normal(2d, 0.8d), new Normal(6d, 1.4d) }),
                new Mixture(new[] { 0.25d, 0.45d, 0.3d }, new UnivariateDistributionBase[]
                    { new Normal(2d, 0.8d), new Normal(5d, 1.2d), new Normal(9d, 1.5d) }),
                new Mixture(new[] { 0.4d, 0.6d }, new UnivariateDistributionBase[]
                    { new Normal(2d, 0.8d), new Normal(6d, 1.4d) }) { IsZeroInflated = true, ZeroWeight = 0.2d }
            };
            double checksum = 0d;
            foreach (Mixture mixture in mixtures) checksum += mixture.LogPDF(3d);

            long before = GC.GetAllocatedBytesForCurrentThread();
            const int repetitions = 100;
            for (int repetition = 0; repetition < repetitions; repetition++)
                foreach (Mixture mixture in mixtures) checksum += mixture.LogPDF(3d);
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

            Assert.IsTrue(double.IsFinite(checksum));
            TestContext.WriteLine($"Measured mixture LogPDF allocation: {allocated} bytes for {repetitions * mixtures.Length} calls.");
            Assert.IsLessThanOrEqualTo(256L * repetitions * mixtures.Length, allocated,
                $"Repeated mixture LogPDF evaluation allocated {allocated} bytes.");
        }
#endif

        /// <summary>Evaluation observes public weight, component and nested-parameter mutations on every call.</summary>
        [TestMethod]
        public void LogPDF_DirectValidationDetectsEveryLiveMutation()
        {
            var mixture = CreateTwoNormalMixture(new[] { 0.4d, 0.6d });
            mixture.Weights[0] = double.NaN;
            ArgumentOutOfRangeException weightError =
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual("Weights", weightError.ParamName);

            mixture = CreateTwoNormalMixture(new[] { 0.4d, 0.6d });
            mixture.Distributions[0].SetParameters(new[] { 2d, -0.8d });
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));

            mixture = CreateTwoNormalMixture(new[] { 0.4d, 0.6d });
            mixture.Distributions[0] = new Normal(double.NaN, 0.8d);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
        }

        /// <summary>Zero-weight components remain parameter-valid while positive-mass checks remain active-only.</summary>
        [TestMethod]
        public void LogPDF_DirectValidationPreservesInactiveAndZeroInflatedRules()
        {
            var invalidInactive = CreateTwoNormalMixture(new[] { 1d, 0d });
            invalidInactive.Distributions[1].SetParameters(new[] { 6d, -1.4d });
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => invalidInactive.LogPDF(3d));

            var validInactive = new Mixture(new[] { 1d, 0d }, new UnivariateDistributionBase[]
                { new Normal(2d, 0.8d), new Deterministic(0d) })
            {
                IsZeroInflated = true,
                ZeroWeight = 0.2d
            };
            Assert.IsTrue(Numerics.Tools.IsFinite(validInactive.LogPDF(3d)));

            validInactive.Weights[0] = 0d;
            validInactive.Weights[1] = 0.8d;
            ArgumentOutOfRangeException positiveMassError =
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => validInactive.LogPDF(3d));
            StringAssert.Contains(positiveMassError.Message, "positive probability above zero");
        }

        /// <summary>The public validator continues to evaluate supplied hypothetical flattened candidates.</summary>
        [TestMethod]
        public void ValidateParameters_PreservesHypotheticalCandidateContract()
        {
            var mixture = CreateTwoNormalMixture(new[] { 0.4d, 0.6d });
            double[] valid = mixture.GetParameters;
            Assert.IsNull(mixture.ValidateParameters(valid, false));

            double[] invalidWeight = (double[])valid.Clone();
            invalidWeight[0] = double.NaN;
            ArgumentOutOfRangeException weightError =
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.ValidateParameters(invalidWeight, true));
            Assert.AreEqual("Weights", weightError.ParamName);

            double[] invalidComponent = (double[])valid.Clone();
            invalidComponent[3] = -0.8d;
            Assert.IsNotNull(mixture.ValidateParameters(invalidComponent, false));
        }

        private static Mixture CreateTwoNormalMixture(double[] weights) => new Mixture(weights,
            new UnivariateDistributionBase[] { new Normal(2d, 0.8d), new Normal(6d, 1.4d) });
    }
}
