using System;
using System.IO;
#if NET8_0_OR_GREATER
using System.Runtime.Serialization;
#else
using System.Runtime.Serialization.Formatters.Binary;
#endif
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Guards density and positive-mass normalization when parameter-dependent caches must be refreshed or restored.</summary>
    [TestClass]
    public class Test_NormalLogDensityCache
    {
        /// <summary>Direct and bulk parameter changes refresh the density's scale-dependent value.</summary>
        [TestMethod]
        public void LogPDF_ObservesDirectAndBulkScaleChanges()
        {
            var normal = new Normal(2d, 0.8d);
            Assert.AreEqual(-1.477044981890463d, normal.LogPDF(3d), 2E-15);
            normal.Sigma = 1.4d;
            Assert.AreEqual(-1.5105128106422123d, normal.LogPDF(3d), 2E-15);
            normal.SetParameters(-4d, 2.5d);
            Assert.AreEqual(-5.755229265078827d, normal.LogPDF(3d), 8E-15);
            normal.SetParameters(new[] { 1d, 4d });
            Assert.AreEqual(-2.4302328943245635d, normal.LogPDF(3d), 4E-15);
        }

        /// <summary>Invalid scales retain validation behavior and can subsequently be replaced by valid scales.</summary>
        [TestMethod]
        public void LogPDF_RecoversAfterInvalidScaleAndLocationChanges()
        {
            var normal = new Normal(2d, 0.8d);
            foreach (double invalid in new[] { -1d, double.NaN, double.NegativeInfinity, double.PositiveInfinity })
            {
                normal.Sigma = invalid;
                Assert.IsFalse(normal.ParametersValid);
                var error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => normal.LogPDF(3d));
                Assert.AreEqual(nameof(Normal.Sigma), error.ParamName);
                normal.Sigma = 1.4d;
                Assert.IsTrue(normal.ParametersValid);
                Assert.AreEqual(-1.5105128106422123d, normal.LogPDF(3d), 2E-15);
            }
            normal.Mu = double.NaN;
            normal.Sigma = 0.8d;
            var locationError = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => normal.LogPDF(3d));
            Assert.AreEqual(nameof(Normal.Mu), locationError.ParamName);
            normal.Mu = 2d;
            Assert.AreEqual(-1.477044981890463d, normal.LogPDF(3d), 2E-15);
        }

        /// <summary>The cached value uses the existing normalized scale, including signed zero and tiny positive inputs.</summary>
        [TestMethod]
        public void LogPDF_UsesNormalizedScaleAndRetainsExtremeDensity()
        {
            var normal = new Normal(2d, 0.8d);
            foreach (double scale in new[] { 0d, -0d, double.Epsilon, 1E-20, 1E-16 })
            {
                normal.Sigma = scale;
                Assert.AreEqual(1E-16, normal.Sigma);
                Assert.AreEqual(35.92242295470006d, normal.LogPDF(2d), 2E-14);
            }
            normal.SetParameters(0d, 1E308);
            Assert.AreEqual(-710.1151471753708d, normal.LogPDF(0d), 2E-12);
        }

        /// <summary>Clones initialize their own scale-dependent density state and remain independent after mutation.</summary>
        [TestMethod]
        public void LogPDF_CloneRetainsIndependentScaleState()
        {
            var original = new Normal(2d, 0.8d);
            _ = original.LogPDF(3d);
            var clone = (Normal)original.Clone();
            Assert.AreEqual(-1.477044981890463d, clone.LogPDF(3d), 2E-15);
            clone.Sigma = 1.4d;
            Assert.AreEqual(-1.5105128106422123d, clone.LogPDF(3d), 2E-15);
            Assert.AreEqual(-1.477044981890463d, original.LogPDF(3d), 2E-15);
        }

        /// <summary>XML round trips preserve parameters and density before and after a subsequent scale change.</summary>
        [TestMethod]
        public void LogPDF_XmlRoundTripRetainsMutableScaleState()
        {
            var original = new Normal(2d, 0.8d);
            _ = original.LogPDF(3d);
            var element = original.ToXElement();
            var restored = (Normal)UnivariateDistributionFactory.CreateDistribution(element);
            Assert.AreEqual(element.ToString(), restored.ToXElement().ToString());
            Assert.AreEqual(-1.477044981890463d, restored.LogPDF(3d), 2E-15);
            restored.Sigma = 1.4d;
            Assert.AreEqual(-1.5105128106422123d, restored.LogPDF(3d), 2E-15);
        }

        /// <summary>Field-based deserialization restores density even when constructors and scale setters do not run.</summary>
        [TestMethod]
        public void LogPDF_SerializableRoundTripInitializesTransientScaleState()
        {
            var original = new Normal(2d, 0.8d);
            _ = original.LogPDF(3d);
            var restored = RoundTrip(original);
            Assert.AreEqual(-1.477044981890463d, restored.LogPDF(3d), 2E-15);
            Assert.AreEqual(-1.477044981890463d, restored.LogPDF(3d), 2E-15);
            restored.Sigma = 1.4d;
            Assert.AreEqual(-1.5105128106422123d, restored.LogPDF(3d), 2E-15);

            original.Sigma = -1d;
            restored = RoundTrip(original);
            var error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => restored.LogPDF(3d));
            Assert.AreEqual(nameof(Normal.Sigma), error.ParamName);
            restored.Sigma = 1.4d;
            Assert.AreEqual(-1.5105128106422123d, restored.LogPDF(3d), 2E-15);
        }

        /// <summary>Repeated zero-survival normalization matches the uncached public method after either parameter changes.</summary>
        [TestMethod]
        public void LogCCDFAtZero_ObservesDirectAndBulkParameterChangesExactly()
        {
            var normal = new Normal(2d, 0.8d);
            AssertZeroSurvivalMatchesUncached(normal);
            normal.Mu = -2d;
            AssertZeroSurvivalMatchesUncached(normal);
            normal.Sigma = 1.4d;
            AssertZeroSurvivalMatchesUncached(normal);
            normal.SetParameters(-4d, 2.5d);
            AssertZeroSurvivalMatchesUncached(normal);
            normal.SetParameters(new[] { 1d, 4d });
            AssertZeroSurvivalMatchesUncached(normal);
            normal.SetParameters(0d, 0d);
            Assert.AreEqual(1E-16, normal.Sigma);
            AssertZeroSurvivalMatchesUncached(normal);
            normal.SetParameters(-40d, 1d);
            Assert.AreEqual(-804.6084420137538d, normal.LogCCDFAtZero(), 2E-12);
            AssertZeroSurvivalMatchesUncached(normal);
        }

        /// <summary>Invalid parameters are rejected before a cached value can be returned and repaired values are observed.</summary>
        [TestMethod]
        public void LogCCDFAtZero_PreservesValidationOrderAndRecovery()
        {
            var normal = new Normal(2d, 0.8d);
            AssertZeroSurvivalMatchesUncached(normal);
            normal.Mu = double.NaN;
            normal.Sigma = -1d;
            var error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => normal.LogCCDFAtZero());
            Assert.AreEqual(nameof(Normal.Mu), error.ParamName);
            normal.Mu = -2d;
            error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => normal.LogCCDFAtZero());
            Assert.AreEqual(nameof(Normal.Sigma), error.ParamName);
            normal.Sigma = 1.4d;
            AssertZeroSurvivalMatchesUncached(normal);
            normal.Sigma = double.PositiveInfinity;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => normal.LogCCDFAtZero());
            normal.SetParameters(1d, 4d);
            AssertZeroSurvivalMatchesUncached(normal);
        }

        /// <summary>Clone, XML, and field-based serialization retain zero-survival values and subsequent parameter mutation.</summary>
        [TestMethod]
        public void LogCCDFAtZero_CloneAndSerializationRestoreTransientState()
        {
            var original = new Normal(-2d, 1.4d);
            AssertZeroSurvivalMatchesUncached(original);
            var copies = new[]
            {
                (Normal)original.Clone(),
                (Normal)UnivariateDistributionFactory.CreateDistribution(original.ToXElement()),
                RoundTrip(original)
            };
            foreach (Normal copy in copies)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(original.LogCCDF(0d)),
                    BitConverter.DoubleToInt64Bits(copy.LogCCDFAtZero()));
                AssertZeroSurvivalMatchesUncached(copy);
                copy.Mu = 2d;
                AssertZeroSurvivalMatchesUncached(copy);
                copy.Sigma = 0.8d;
                AssertZeroSurvivalMatchesUncached(copy);
            }
            original.Mu = double.NaN;
            Normal invalid = RoundTrip(original);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => invalid.LogCCDFAtZero());
            invalid.Mu = 2d;
            AssertZeroSurvivalMatchesUncached(invalid);
        }

        /// <summary>Zero-inflated density retains exact uncached normalization after live Normal parameter changes.</summary>
        [TestMethod]
        public void ZeroInflatedLogPDF_ObservesNormalParameterChangesExactly()
        {
            var mixture = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { new Normal(2d, 0.8d) })
            { IsZeroInflated = true, ZeroWeight = 0.2d };
            var normal = (Normal)mixture.Distributions[0];
            AssertSingleComponentNormalization(mixture);
            normal.Mu = -2d;
            AssertSingleComponentNormalization(mixture);
            normal.Sigma = 1.4d;
            AssertSingleComponentNormalization(mixture);
            normal.SetParameters(-4d, 2.5d);
            AssertSingleComponentNormalization(mixture);
            normal.SetParameters(new[] { 1d, 4d });
            AssertSingleComponentNormalization(mixture);
            normal.Mu = double.NaN;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(0d));
            normal.Mu = 2d;
            AssertSingleComponentNormalization(mixture);
        }

        /// <summary>Cached positive mass cannot bypass weight, inactive-parameter, or component-order validation.</summary>
        [TestMethod]
        public void ZeroInflatedLogPDF_PreservesLiveValidationPrecedence()
        {
            var mixture = new Mixture(new[] { 0.4d, 0.6d },
                new UnivariateDistributionBase[] { new Normal(2d, 0.8d), new Normal(6d, 1.4d) })
            { IsZeroInflated = true, ZeroWeight = 0.2d };
            var first = (Normal)mixture.Distributions[0];
            var second = (Normal)mixture.Distributions[1];
            _ = mixture.LogPDF(3d);
            first.SetParameters(-1E308, 1E-16);
            second.Sigma = -1d;
            mixture.Weights[0] = double.NaN;
            var error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual(nameof(Mixture.Weights), error.ParamName);
            mixture.Weights[0] = 0.32d;
            error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual(nameof(Mixture.Distributions), error.ParamName);
            StringAssert.Contains(error.Message, "positive probability above zero");
            mixture.Weights[0] = 0d;
            mixture.Weights[1] = 0.8d;
            first.Mu = double.NaN;
            error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual(nameof(Normal.Mu), error.ParamName);
            first.Mu = -1E308;
            error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual(nameof(Normal.Sigma), error.ParamName);
            second.Sigma = 1.4d;
            Assert.IsTrue(Numerics.Tools.IsFinite(mixture.LogPDF(3d)));
        }

        /// <summary>Non-Normal mixture components continue to evaluate their live zero-survival call after nested mutation.</summary>
        [TestMethod]
        public void ZeroInflatedLogPDF_PreservesNestedComponentNormalization()
        {
            var inner = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { new Normal(-2d, 1.4d) });
            var outer = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { inner })
            { IsZeroInflated = true, ZeroWeight = 0.2d };
            var nested = (Mixture)outer.Distributions[0];
            var normal = (Normal)nested.Distributions[0];
            AssertSingleComponentNormalization(outer);
            normal.Mu = 2d;
            AssertSingleComponentNormalization(outer);
            normal.Sigma = 0.8d;
            AssertSingleComponentNormalization(outer);
            normal.Sigma = -1d;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => outer.LogPDF(3d));
            normal.Sigma = 1.4d;
            AssertSingleComponentNormalization(outer);
        }

        private static void AssertZeroSurvivalMatchesUncached(Normal normal)
        {
            long expected = BitConverter.DoubleToInt64Bits(normal.LogCCDF(0d));
            for (int i = 0; i < 3; i++)
                Assert.AreEqual(expected, BitConverter.DoubleToInt64Bits(normal.LogCCDFAtZero()));
        }

        private static void AssertSingleComponentNormalization(Mixture mixture)
        {
            var component = mixture.Distributions[0];
            double expected = Math.Log(mixture.Weights[0]) + (component.LogPDF(3d) - component.LogCCDF(0d));
            for (int i = 0; i < 3; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(mixture.LogPDF(3d)));
        }

        private static Normal RoundTrip(Normal original)
        {
            using var stream = new MemoryStream();
#if NET8_0_OR_GREATER
            // Honors Serializable/NonSerialized field contracts without enabling the obsolete BinaryFormatter.
            var serializer = new DataContractSerializer(typeof(Normal));
            serializer.WriteObject(stream, original);
            stream.Position = 0;
            return (Normal)serializer.ReadObject(stream);
#else
            // The legacy target still supports the actual existing binary serialization contract.
            var serializer = new BinaryFormatter();
            serializer.Serialize(stream, original);
            stream.Position = 0;
            return (Normal)serializer.Deserialize(stream);
#endif
        }
    }
}
