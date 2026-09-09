using System;
using System.Collections.Generic;
using System.IO;
#if NET8_0_OR_GREATER
using System.Runtime.Serialization;
#else
using System.Runtime.Serialization.Formatters.Binary;
#endif
using System.Threading.Tasks;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Guards live weight reads, callbacks, and transient caches in mixture density evaluation.</summary>
    [TestClass]
    public class Test_MixtureWeightLogCache
    {
#if NET8_0_OR_GREATER
        /// <summary>Repeated evaluations with unchanged weights do not allocate after initial use.</summary>
        /// <param name="count">The number of Normal components.</param>
        /// <param name="zeroInflated">Whether to include the hurdle normalization path.</param>
        [TestMethod]
        [DataRow(2, false)]
        [DataRow(3, false)]
        [DataRow(2, true)]
        public void LogPDF_RepeatedStableWeightsAllocateZeroBytes(int count, bool zeroInflated)
        {
            var mixture = count == 2 ? Create() : new Mixture(new[] { 0.25d, 0.45d, 0.3d },
                new UnivariateDistributionBase[] { new Normal(2d, 0.8d), new Normal(5d, 1.2d), new Normal(9d, 1.5d) });
            if (zeroInflated) { mixture.IsZeroInflated = true; mixture.ZeroWeight = 0.2d; }
            double checksum = 0d;
            for (int i = 0; i < 1000; i++) checksum += mixture.LogPDF(3d);
            long before = GC.GetAllocatedBytesForCurrentThread();
            for (int i = 0; i < 1000; i++) checksum += mixture.LogPDF(3d);
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
            Assert.IsTrue(Numerics.Tools.IsFinite(checksum));
            Assert.AreEqual(0L, allocated);
        }
#endif

        /// <summary>In-place, bulk, replacement, and zero-inflation weight changes match a freshly constructed mixture exactly.</summary>
        [TestMethod]
        public void LogPDF_ObservesEveryWeightUpdateExactly()
        {
            Mixture mixture = Create();
            AssertMatchesFresh(mixture);
            mixture.Weights[0] = 0.7d;
            mixture.Weights[1] = 0.3d;
            AssertMatchesFresh(mixture);
            double[] parameters = mixture.GetParameters;
            parameters[0] = 0.25d;
            parameters[1] = 0.75d;
            mixture.SetParameters(parameters);
            AssertMatchesFresh(mixture);
            mixture.SetParameters(new[] { 0.2d, 0.3d, 0.5d },
                new UnivariateDistributionBase[] { new Normal(1d, 0.8d), new Normal(5d, 1.2d), new Normal(9d, 1.5d) });
            AssertMatchesFresh(mixture);
            mixture.SetParameters(new[] { 0.6d, 0.4d },
                new UnivariateDistributionBase[] { new Normal(2d, 0.8d), new Normal(6d, 1.4d) });
            AssertMatchesFresh(mixture);
            mixture.IsZeroInflated = true;
            mixture.ZeroWeight = 0.2d;
            AssertMatchesFresh(mixture);
            mixture.ZeroWeight = 0.35d;
            AssertMatchesFresh(mixture);
        }

        /// <summary>Warmed caches do not bypass invalid live weights or change weight-first validation and recovery.</summary>
        [TestMethod]
        public void LogPDF_PreservesInvalidWeightPrecedenceAndRecovery()
        {
            Mixture mixture = Create();
            AssertMatchesFresh(mixture);
            var normal = (Normal)mixture.Distributions[0];
            normal.Mu = double.NaN;
            mixture.Weights[0] = double.NaN;
            var error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual(nameof(Mixture.Weights), error.ParamName);
            mixture.Weights[0] = 0.4d;
            error = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => mixture.LogPDF(3d));
            Assert.AreEqual(nameof(Normal.Mu), error.ParamName);
            normal.Mu = 2d;
            AssertMatchesFresh(mixture);
            mixture.Weights[0] = 0d;
            mixture.Weights[1] = 1d;
            AssertMatchesFresh(mixture);
            mixture.Weights[0] = 0.4d;
            mixture.Weights[1] = 0.6d;
            AssertMatchesFresh(mixture);
        }

        /// <summary>Each weight is read after its component callback, without anticipating a later component's mutation.</summary>
        [TestMethod]
        public void LogPDF_ReadsWeightAfterEachComponentCallback()
        {
            var order = new List<int>();
            Mixture mixture = null;
            var first = new CallbackDistribution(() =>
            { order.Add(1); mixture.Weights[0] = 0.25d; mixture.Weights[1] = 0.75d; });
            var second = new CallbackDistribution(() =>
            { order.Add(2); mixture.Weights[0] = 0.5d; mixture.Weights[1] = 0.5d; });
            mixture = new Mixture(new[] { 0.4d, 0.6d }, new UnivariateDistributionBase[] { first, second });
            for (int i = 0; i < 3; i++)
            {
                order.Clear();
                Assert.AreEqual(Math.Log(0.75d), mixture.LogPDF(3d), 2E-15);
                CollectionAssert.AreEqual(new[] { 1, 2 }, order);
            }
        }

        /// <summary>Weights changed by a callback retain the original logarithm's signed-zero and nonfinite semantics.</summary>
        [TestMethod]
        public void LogPDF_CallbackWeightChangesRetainOriginalLogSemantics()
        {
            double target = 1d;
            Mixture mixture = null;
            mixture = new Mixture(new[] { 1d }, new UnivariateDistributionBase[]
                { new CallbackDistribution(() => mixture.Weights[0] = target) });
            foreach (double value in new[] { 0.25d, 0d, -0d, -0.5d, double.NaN,
                BitConverter.Int64BitsToDouble(0x7FF8000000000123L), double.PositiveInfinity, double.NegativeInfinity, 1d })
            {
                target = value;
                mixture.Weights[0] = 1d;
                double expected = DistributionNumerics.LogSum(double.NegativeInfinity, Math.Log(value) + 0d);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(mixture.LogPDF(3d)));
            }
        }

        /// <summary>A prior callback can activate or deactivate a later component before that component's existing weight check.</summary>
        [TestMethod]
        public void LogPDF_PreservesLiveActivationAndCallbackOrder()
        {
            var order = new List<int>();
            bool activateSecond = false;
            Mixture mixture = null;
            var first = new CallbackDistribution(() =>
            {
                order.Add(1);
                mixture.Weights[0] = activateSecond ? 0.5d : 1d;
                mixture.Weights[1] = activateSecond ? 0.5d : 0d;
            });
            var second = new CallbackDistribution(() => order.Add(2));
            mixture = new Mixture(new[] { 0.5d, 0.5d }, new UnivariateDistributionBase[] { first, second });
            Assert.AreEqual(0d, mixture.LogPDF(3d), 0d);
            CollectionAssert.AreEqual(new[] { 1 }, order);
            activateSecond = true;
            order.Clear();
            Assert.AreEqual(0d, mixture.LogPDF(3d), 2E-15);
            CollectionAssert.AreEqual(new[] { 1, 2 }, order);
        }

        /// <summary>Concurrent first use on an unchanged mixture publishes complete weight/log pairs.</summary>
        [TestMethod]
        public void LogPDF_ConcurrentReadOnlyInitializationPreservesExactValues()
        {
            Mixture mixture = Create();
            long expected = BitConverter.DoubleToInt64Bits(mixture.Clone().LogPDF(3d));
            var results = new double[32];
            Parallel.For(0, results.Length, i => results[i] = mixture.LogPDF(3d));
            foreach (double result in results) Assert.AreEqual(expected, BitConverter.DoubleToInt64Bits(result));
        }

        /// <summary>Clone, XML, and field-based deserialization restore transient weight caches and retain later mutations.</summary>
        [TestMethod]
        public void LogPDF_CloneAndSerializationRetainMutableWeightState()
        {
            Mixture mixture = Create();
            AssertMatchesFresh(mixture);
            long expected = BitConverter.DoubleToInt64Bits(mixture.LogPDF(3d));
            var copies = new[] { (Mixture)mixture.Clone(),
                (Mixture)UnivariateDistributionFactory.CreateDistribution(mixture.ToXElement()), RoundTrip(mixture) };
            foreach (Mixture copy in copies)
            {
                Assert.AreEqual(expected, BitConverter.DoubleToInt64Bits(copy.LogPDF(3d)));
                copy.Weights[0] = 0.7d;
                copy.Weights[1] = 0.3d;
                AssertMatchesFresh(copy);
            }
        }

        private static void AssertMatchesFresh(Mixture mixture)
        {
            long expected = BitConverter.DoubleToInt64Bits(mixture.Clone().LogPDF(3d));
            for (int i = 0; i < 3; i++) Assert.AreEqual(expected, BitConverter.DoubleToInt64Bits(mixture.LogPDF(3d)));
        }

        private static Mixture Create() => new Mixture(new[] { 0.4d, 0.6d },
            new UnivariateDistributionBase[] { new Normal(2d, 0.8d), new Normal(6d, 1.4d) });

        private static Mixture RoundTrip(Mixture mixture)
        {
            using var stream = new MemoryStream();
#if NET8_0_OR_GREATER
            var serializer = new DataContractSerializer(typeof(Mixture), new[] { typeof(Normal) });
            serializer.WriteObject(stream, mixture);
            stream.Position = 0;
            return (Mixture)serializer.ReadObject(stream);
#else
            var serializer = new BinaryFormatter();
            serializer.Serialize(stream, mixture);
            stream.Position = 0;
            return (Mixture)serializer.Deserialize(stream);
#endif
        }

        private sealed class CallbackDistribution : Cauchy
        {
            private readonly Action _callback;

            internal CallbackDistribution(Action callback) => _callback = callback;

            /// <inheritdoc/>
            public override double LogPDF(double x)
            {
                _callback();
                return 0d;
            }

            /// <inheritdoc/>
            public override UnivariateDistributionBase Clone() => new CallbackDistribution(_callback);
        }
    }
}
