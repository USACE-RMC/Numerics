using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
#if NETFRAMEWORK
using System.Runtime.Serialization.Formatters.Binary;
#endif
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Guards canonical competing-risk cache invalidation while removing repeated built-in Weibull XML work.</summary>
    [TestClass]
    public class Test_CompetingRisksConfigurationCache
    {
#if NET8_0_OR_GREATER
        /// <summary>Unchanged below-support CDF calls allocate only the same validation/support work as LogPDF.</summary>
        [TestMethod]
        public void CDF_UnchangedWeibullsAvoidConfigurationSerializationAllocations()
        {
            CompetingRisks risks = Create();
            for (int i = 0; i < 100; i++) { _ = risks.CDF(-1d); _ = risks.LogPDF(-1d); }
            long before = GC.GetAllocatedBytesForCurrentThread();
            for (int i = 0; i < 100; i++) _ = risks.LogPDF(-1d);
            long supportWork = GC.GetAllocatedBytesForCurrentThread() - before;
            before = GC.GetAllocatedBytesForCurrentThread();
            for (int i = 0; i < 100; i++) _ = risks.CDF(-1d);
            long cdfWork = GC.GetAllocatedBytesForCurrentThread() - before;
            Assert.AreEqual(supportWork, cdfWork,
                $"Repeated CDF allocated {cdfWork} bytes; matching validation/support work allocated {supportWork} bytes.");
        }
#endif

        /// <summary>Ordered scalar mutations, replacement, and count changes retain the original canonical invalidation decision.</summary>
        [TestMethod]
        public void Configuration_ObservesWeibullMutationsAndReplacement()
        {
            CompetingRisks risks = Create();
            AssertCanonicalInvalidation(risks, () => ((Weibull)risks.Distributions[0]).Lambda = 7d);
            AssertCanonicalInvalidation(risks, () => ((Weibull)risks.Distributions[1]).Kappa = 3d);
            AssertCanonicalInvalidation(risks, () => risks.SetParameters(new[] { 8d, 1.5d, 11d, 2.5d }));
            AssertCanonicalInvalidation(risks, () => risks.SetParameters(new UnivariateDistributionBase[]
                { new Weibull(8d, 1.5d), new Weibull(11d, 2.5d) }));
            AssertCanonicalInvalidation(risks, () => risks.SetParameters(new UnivariateDistributionBase[]
                { new Weibull(11d, 2.5d), new Weibull(8d, 1.5d) }));
            AssertCanonicalInvalidation(risks, () => risks.SetParameters(new UnivariateDistributionBase[]
                { new Weibull(11d, 2.5d) }));
        }

        /// <summary>Flags, seed, matrix shape, live values, and signed zeros follow canonical string equality on each runtime.</summary>
        [TestMethod]
        public void Configuration_ObservesFlagsSeedAndEveryMatrixValue()
        {
            CompetingRisks risks = Create();
            AssertCanonicalInvalidation(risks, () => risks.MinimumOfRandomVariables = false);
            AssertCanonicalInvalidation(risks, () => risks.Dependency = Probability.DependencyType.Independent);
            AssertCanonicalInvalidation(risks, () => risks.PRNGSeed = 9876);
            AssertCanonicalInvalidation(risks, () => risks.XTransform = Transform.Logarithmic);
            AssertCanonicalInvalidation(risks, () => risks.ProbabilityTransform = Transform.None);
            AssertCanonicalInvalidation(risks, () => risks.CorrelationMatrix = new[,] { { 1d, 0d }, { 0d, 1d } });
            AssertCanonicalInvalidation(risks, () => risks.CorrelationMatrix[0, 1] = -0d);
            AssertCanonicalInvalidation(risks, () => risks.CorrelationMatrix[1, 0] = 0.3d);
            AssertCanonicalInvalidation(risks, () => risks.CorrelationMatrix = (double[,])risks.CorrelationMatrix.Clone());
            AssertCanonicalInvalidation(risks, () => risks.CorrelationMatrix = new double[1, 3]);
            AssertCanonicalInvalidation(risks, () => risks.CorrelationMatrix = null);
        }

        /// <summary>Matrix lower bounds do not add an invalidation distinction absent from the original canonical state.</summary>
        [TestMethod]
        public void Configuration_RetainsMatrixEnumerationSemantics()
        {
            CompetingRisks risks = Create();
            risks.CorrelationMatrix = new[,] { { 1d, 0.2d }, { 0.2d, 1d } };
            AssertCanonicalInvalidation(risks, () =>
            {
                var matrix = (double[,])Array.CreateInstance(typeof(double), new[] { 2, 2 }, new[] { 2, 3 });
                matrix[2, 3] = 1d; matrix[2, 4] = 0.2d;
                matrix[3, 3] = 0.2d; matrix[3, 4] = 1d;
                risks.CorrelationMatrix = matrix;
            });
        }

        /// <summary>Derived Weibull and nested/custom children still execute their original XML callbacks on every refresh.</summary>
        [TestMethod]
        public void Configuration_PreservesDerivedAndNestedFallbackCallbacks()
        {
            var child = new CountingWeibull(5d, 1.4d);
            var risks = new CompetingRisks(new UnivariateDistributionBase[] { child });
            _ = risks.InverseCDF(0d);
            int before = child.SerializationCalls;
            _ = risks.InverseCDF(0d);
            _ = risks.InverseCDF(0d);
            Assert.AreEqual(before + 2, child.SerializationCalls);
            child.Token = "changed";
            _ = risks.InverseCDF(0d);
            StringAssert.Contains(CachedConfiguration(risks), "changed");
            var nested = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { child });
            var nestedChild = (CountingWeibull)nested.Distributions[0];
            risks.SetParameters(new UnivariateDistributionBase[] { nested });
            _ = risks.InverseCDF(0d);
            before = nestedChild.SerializationCalls;
            _ = risks.InverseCDF(0d);
            Assert.AreEqual(before + 1, nestedChild.SerializationCalls);
            risks.SetParameters(new UnivariateDistributionBase[] { new Weibull(5d, 1.4d) });
            AssertCanonicalInvalidation(risks, () => ((Weibull)risks.Distributions[0]).Lambda = 6d);
        }

        /// <summary>Correlated CDF evaluations retain the exact sequence of the unchanged generic configuration path.</summary>
        [TestMethod]
        public void CDF_CorrelatedWeibullsMatchGenericRefreshExactly()
        {
            CompetingRisks fast = Create();
            var generic = new CompetingRisks(new UnivariateDistributionBase[]
                { new CountingWeibull(5d, 1.4d), new CountingWeibull(12d, 2.3d) });
            foreach (CompetingRisks risks in new[] { fast, generic })
            {
                risks.Dependency = Probability.DependencyType.CorrelationMatrix;
                risks.CorrelationMatrix = new[,] { { 1d, 0.35d }, { 0.35d, 1d } };
                risks.PRNGSeed = 417;
            }
            for (int round = 0; round < 2; round++)
            {
                foreach (double x in new[] { 0.5d, 1d, 3d, 8d, 15d })
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(generic.CDF(x)), BitConverter.DoubleToInt64Bits(fast.CDF(x)));
                ((Weibull)fast.Distributions[0]).Lambda = 7d;
                ((Weibull)generic.Distributions[0]).Lambda = 7d;
            }
        }

        /// <summary>Clone, XML, and supported legacy binary serialization restore cold snapshots without changing canonical state.</summary>
        [TestMethod]
        public void Configuration_SerializationRestoresColdSnapshot()
        {
            CompetingRisks original = Create();
            _ = original.InverseCDF(0d);
            string expected = CachedConfiguration(original);
            var copies = new List<CompetingRisks> { (CompetingRisks)original.Clone(),
                (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(original.ToXElement()) };
#if NETFRAMEWORK
            copies.Add(RoundTrip(original));
#endif
            foreach (CompetingRisks copy in copies)
            {
                _ = copy.InverseCDF(0d);
                Assert.AreEqual(expected, CachedConfiguration(copy));
                AssertCanonicalInvalidation(copy, () => ((Weibull)copy.Distributions[0]).Lambda = 7d);
            }
        }

        private static readonly FieldInfo Configuration = typeof(CompetingRisks).GetField("_cachedConfiguration", BindingFlags.Instance | BindingFlags.NonPublic);
        private static readonly FieldInfo[] DerivedFlags =
        {
            typeof(CompetingRisks).GetField("_momentsComputed", BindingFlags.Instance | BindingFlags.NonPublic),
            typeof(CompetingRisks).GetField("_empiricalCDFCreated", BindingFlags.Instance | BindingFlags.NonPublic),
            typeof(CompetingRisks).GetField("_mvnCreated", BindingFlags.Instance | BindingFlags.NonPublic)
        };

        private static string CachedConfiguration(CompetingRisks risks) => (string)Configuration.GetValue(risks);

        private static void AssertCanonicalInvalidation(CompetingRisks risks, Action mutate)
        {
            _ = risks.InverseCDF(0d);
            string before = CachedConfiguration(risks);
            mutate();
            string expected = DistributionNumerics.ConfigurationState(risks);
            // Set sentinels after setters so this checks refresh's canonical decision independently of setter invalidation.
            foreach (FieldInfo flag in DerivedFlags) flag.SetValue(risks, true);
            _ = risks.InverseCDF(0d);
            Assert.AreEqual(expected, CachedConfiguration(risks));
            foreach (FieldInfo flag in DerivedFlags)
                Assert.AreEqual(expected == before, (bool)flag.GetValue(risks), flag.Name);
        }

        private static CompetingRisks Create() => new CompetingRisks(new UnivariateDistributionBase[]
            { new Weibull(5d, 1.4d), new Weibull(12d, 2.3d) }) { Dependency = Probability.DependencyType.PerfectlyPositive };

#if NETFRAMEWORK
        private static CompetingRisks RoundTrip(CompetingRisks risks)
        {
            using var stream = new MemoryStream();
            var serializer = new BinaryFormatter();
            serializer.Serialize(stream, risks);
            stream.Position = 0;
            return (CompetingRisks)serializer.Deserialize(stream);
        }
#endif

        private sealed class CountingWeibull : Weibull
        {
            internal int SerializationCalls;
            internal string Token = "original";
            internal CountingWeibull(double scale, double shape) : base(scale, shape) { }
            /// <inheritdoc/>
            public override XElement ToXElement()
            {
                SerializationCalls++;
                XElement element = base.ToXElement();
                element.SetAttributeValue("CustomToken", Token);
                return element;
            }
            /// <inheritdoc/>
            public override UnivariateDistributionBase Clone() => new CountingWeibull(Lambda, Kappa) { Token = Token };
        }
    }
}
