using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Guards allocation-free live Normal validation without weakening mixture validation contracts.</summary>
    [TestClass]
    public class Test_NormalMixtureValidationAllocations
    {
#if NET8_0_OR_GREATER
        /// <summary>Repeated Normal mixture densities validate live parameters without allocating parameter arrays.</summary>
        /// <param name="zeroInflated">Whether to exercise the positive-conditional mixture path.</param>
        [TestMethod]
        [DataRow(false)]
        [DataRow(true)]
        public void LogPDF_LiveNormalValidationAllocatesZeroBytes(bool zeroInflated)
        {
            Mixture mixture = Create();
            if (zeroInflated) { mixture.IsZeroInflated = true; mixture.ZeroWeight = 0.2d; }
            double checksum = 0d;
            for (int i = 0; i < 1000; i++) checksum += mixture.LogPDF(3d);

            long before = GC.GetAllocatedBytesForCurrentThread();
            for (int i = 0; i < 1000; i++) checksum += mixture.LogPDF(3d);
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

            Assert.IsTrue(double.IsFinite(checksum));
            Assert.AreEqual(0L, allocated, $"1000 evaluations allocated {allocated} bytes ({allocated / 1000d} bytes/call).");
        }
#endif

        /// <summary>Live scalar validation preserves weight-first and component-order errors and observes repaired values.</summary>
        [TestMethod]
        public void LogPDF_PreservesErrorOrderAndLiveRecovery()
        {
            Mixture mixture = Create();
            var first = (Normal)mixture.Distributions[0];
            var second = (Normal)mixture.Distributions[1];
            first.Mu = double.NaN;
            first.Sigma = -1d;
            second.Sigma = -2d;
            mixture.Weights[0] = double.NaN;
            var weightError = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual("Weights", weightError.ParamName);

            mixture.Weights[0] = 0.4d;
            AssertScalarError(mixture, first);
            first.Mu = 2d;
            AssertScalarError(mixture, first);
            first.Sigma = 0.8d;
            AssertScalarError(mixture, second);
            second.Sigma = 1.4d;
            Assert.IsFalse(double.IsNaN(mixture.LogPDF(3d)));
            mixture.Distributions[0] = new Normal(double.PositiveInfinity, 1d);
            AssertScalarError(mixture, (Normal)mixture.Distributions[0]);
        }

        /// <summary>Zero-weight Normal components remain validated even when their density will not contribute.</summary>
        [TestMethod]
        public void LogPDF_ValidatesInactiveNormalComponents()
        {
            Mixture mixture = Create();
            mixture.Weights[0] = 1d;
            mixture.Weights[1] = 0d;
            var inactive = (Normal)mixture.Distributions[1];
            inactive.Mu = double.NaN;
            AssertScalarError(mixture, inactive);
        }

        /// <summary>The public candidate validator uses supplied parameters independently of current Normal state.</summary>
        [TestMethod]
        public void ValidateParameters_PreservesCandidateIndependence()
        {
            Mixture mixture = Create();
            double[] validCandidate = mixture.GetParameters;
            double[] invalidCandidate = (double[])validCandidate.Clone();
            invalidCandidate[3] = -0.8d;
            var candidateError = mixture.ValidateParameters(invalidCandidate, false);
            Assert.IsNotNull(candidateError);
            Assert.AreEqual("Sigma", candidateError.ParamName);
            Assert.IsNull(mixture.ValidateParameters(validCandidate, false));

            ((Normal)mixture.Distributions[0]).Mu = double.NaN;
            Assert.IsNull(mixture.ValidateParameters(validCandidate, false));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
        }

        /// <summary>Non-Normal component delegation still detects a directly mutated nested Normal.</summary>
        [TestMethod]
        public void LogPDF_PreservesNestedMutationValidation()
        {
            Mixture inner = Create();
            var outer = new Mixture(new[] { 0.5d, 0.5d }, new UnivariateDistributionBase[] { inner, new Normal(10d, 2d) });
            var nested = (Mixture)outer.Distributions[0];
            ((Normal)nested.Distributions[1]).Sigma = -1d;
            var error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => outer.LogPDF(3d));
            Assert.AreEqual("Sigma", error.ParamName);
        }

        private static void AssertScalarError(Mixture mixture, Normal normal)
        {
            var expected = normal.ValidateParameters(normal.Mu, normal.Sigma, false);
            Assert.IsNotNull(expected);
            var actual = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual(expected.ParamName, actual.ParamName);
            Assert.AreEqual(expected.Message, actual.Message);
            Assert.AreEqual(expected.ActualValue, actual.ActualValue);
        }

        private static Mixture Create() => new Mixture(new[] { 0.4d, 0.6d },
            new UnivariateDistributionBase[] { new Normal(2d, 0.8d), new Normal(6d, 1.4d) });
    }
}
