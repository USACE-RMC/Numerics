using System;
using System.Collections.Generic;
using System.IO;
#if NETFRAMEWORK
using System.Runtime.Serialization.Formatters.Binary;
#endif
using System.Threading.Tasks;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Guards reuse of the identical built-in Weibull derivative step without changing dependent densities.</summary>
    [TestClass]
    public class Test_CompetingRisksDensityStepCache
    {
        /// <summary>Live parameter, component, and configuration changes retain the generic derivative's exact results.</summary>
        [TestMethod]
        public void LogPDF_WeibullStepObservesLiveMutations()
        {
            var children = new UnivariateDistributionBase[] { new Weibull(5d, 1.4d), new Weibull(12d, 2.3d) };
            var risks = Create(children);
            AssertGenericParity(risks);
            ((Weibull)children[0]).Lambda = 7d;
            AssertGenericParity(risks);
            ((Weibull)children[1]).Kappa = 3d;
            AssertGenericParity(risks);
            risks.SetParameters(new[] { 8d, 1.5d, 11d, 2.5d });
            AssertGenericParity(risks);
            children[0] = new Weibull(4d, 1d);
            AssertGenericParity(risks);
            risks.SetParameters(new UnivariateDistributionBase[] { new Weibull(11d, 2.5d), new Weibull(4d, 1d) });
            AssertGenericParity(risks);
            risks.MinimumOfRandomVariables = false;
            risks.PRNGSeed = 9876;
            risks.XTransform = Transform.Logarithmic;
            risks.ProbabilityTransform = Transform.None;
            risks.CorrelationMatrix[0, 1] = risks.CorrelationMatrix[1, 0] = .6d;
            AssertGenericParity(risks);
            risks.Dependency = Probability.DependencyType.PerfectlyPositive;
            AssertGenericParity(risks);
        }

        /// <summary>Cached steps never bypass live validity checks and recover after invalid component parameters are repaired.</summary>
        [TestMethod]
        public void LogPDF_WeibullStepPreservesInvalidRecovery()
        {
            var risks = Create(new UnivariateDistributionBase[] { new Weibull(5d, 1.4d), new Weibull(12d, 2.3d) });
            AssertGenericParity(risks);
            var child = (Weibull)risks.Distributions[0];
            child.Lambda = double.NaN;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => risks.LogPDF(3d));
            child.Lambda = 7d;
            AssertGenericParity(risks);
            child.Kappa = 0d;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => risks.LogPDF(3d));
            child.Kappa = 2d;
            AssertGenericParity(risks);
        }

        /// <summary>Unusable component widths keep the observation-dependent fallback, and unresolvable steps still throw.</summary>
        [TestMethod]
        public void LogPDF_WeibullStepPreservesFallbackAndExceptions()
        {
            var risks = Create(new UnivariateDistributionBase[] { new Weibull(1d, .0001d), new Weibull(2d, .0001d) });
            risks.Dependency = Probability.DependencyType.PerfectlyPositive;
            AssertGenericParity(risks, new[] { 0d, .25d, 1d, 2d, 100d, 1E100d });
            // Both IQRs are finite, but their component-derived step underflows to zero.
            risks.SetParameters(new[] { double.Epsilon, 1d, double.Epsilon, 1d });
            _ = risks.CDF(-1d);
            for (int i = 0; i < 3; i++)
                Assert.ThrowsExactly<InvalidOperationException>(() => risks.LogPDF(0d));
            risks.SetParameters(new[] { 5d, 1.4d, 12d, 2.3d });
            AssertGenericParity(risks);
        }

        /// <summary>Derived Weibull callbacks retain their original order, repeated evaluation, and exception behavior.</summary>
        [TestMethod]
        public void LogPDF_DerivedWeibullsRetainQuartileCallbacks()
        {
            var calls = new List<string>();
            var first = new ObservedWeibull(5d, 1.4d) { Name = "first", Calls = calls };
            var second = new ObservedWeibull(12d, 2.3d) { Name = "second", Calls = calls };
            var risks = Create(new UnivariateDistributionBase[] { new Weibull(5d, 1.4d), new Weibull(12d, 2.3d) });
            AssertGenericParity(risks);
            risks.SetParameters(new UnivariateDistributionBase[] { first, second });
            for (int i = 0; i < 3; i++)
            {
                calls.Clear();
                _ = risks.LogPDF(3d);
                CollectionAssert.AreEqual(new[] { "first:upper", "first:lower", "second:upper", "second:lower" }, calls);
            }
            first.ThrowOnUpper = true;
            calls.Clear();
            var error = Assert.ThrowsExactly<InvalidOperationException>(() => risks.LogPDF(3d));
            Assert.AreEqual("quartile callback", error.Message);
            CollectionAssert.AreEqual(new[] { "first:upper" }, calls);
            first.ThrowOnUpper = false;
            _ = risks.LogPDF(3d);
        }

        /// <summary>Concurrent readers can lazily publish a step for an already initialized unchanged configuration.</summary>
        [TestMethod]
        public void LogPDF_WeibullStepInitializesForConcurrentReaders()
        {
            var risks = Create(new UnivariateDistributionBase[] { new Weibull(5d, 1.4d), new Weibull(12d, 2.3d) });
            risks.Dependency = Probability.DependencyType.PerfectlyPositive;
            _ = risks.CDF(-1d);
            var generic = GenericCopy(risks);
            long expected = BitConverter.DoubleToInt64Bits(generic.LogPDF(3d));
            Parallel.For(0, 32, i => Assert.AreEqual(expected, BitConverter.DoubleToInt64Bits(risks.LogPDF(3d))));
        }

        /// <summary>Clone, XML, and actual legacy binary round trips preserve cold-cache evaluation and subsequent mutations.</summary>
        [TestMethod]
        public void LogPDF_WeibullStepSurvivesSerialization()
        {
            var original = Create(new UnivariateDistributionBase[] { new Weibull(5d, 1.4d), new Weibull(12d, 2.3d) });
            AssertGenericParity(original);
            var copies = new List<CompetingRisks> { (CompetingRisks)original.Clone(),
                (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(original.ToXElement()) };
#if NETFRAMEWORK
            using var stream = new MemoryStream();
            var serializer = new BinaryFormatter();
            serializer.Serialize(stream, original);
            stream.Position = 0;
            copies.Add((CompetingRisks)serializer.Deserialize(stream));
#endif
            foreach (CompetingRisks copy in copies)
            {
                AssertGenericParity(copy);
                ((Weibull)copy.Distributions[0]).Lambda = 7d;
                AssertGenericParity(copy);
            }
        }

        private static CompetingRisks Create(UnivariateDistributionBase[] children) => new CompetingRisks(children)
        {
            Dependency = Probability.DependencyType.CorrelationMatrix,
            CorrelationMatrix = new[,] { { 1d, .35d }, { .35d, 1d } },
            PRNGSeed = 417
        };

        private static CompetingRisks GenericCopy(CompetingRisks risks)
        {
            var children = new UnivariateDistributionBase[risks.Distributions.Count];
            for (int i = 0; i < children.Length; i++)
            {
                var child = (Weibull)risks.Distributions[i];
                children[i] = new ObservedWeibull(child.Lambda, child.Kappa);
            }
            return new CompetingRisks(children)
            {
                MinimumOfRandomVariables = risks.MinimumOfRandomVariables,
                Dependency = risks.Dependency,
                CorrelationMatrix = (double[,])risks.CorrelationMatrix.Clone(),
                PRNGSeed = risks.PRNGSeed,
                XTransform = risks.XTransform,
                ProbabilityTransform = risks.ProbabilityTransform
            };
        }

        private static void AssertGenericParity(CompetingRisks risks, double[] observations = null)
        {
            var generic = GenericCopy(risks);
            observations ??= new[] { 0d, .5d, 1d, 3d, 8d, 15d };
            for (int repeat = 0; repeat < 3; repeat++)
                foreach (double x in observations)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(generic.LogPDF(x)), BitConverter.DoubleToInt64Bits(risks.LogPDF(x)), $"x={x:G17}");
        }

        private sealed class ObservedWeibull : Weibull
        {
            internal string Name;
            internal List<string> Calls;
            internal bool ThrowOnUpper;
            internal ObservedWeibull(double scale, double shape) : base(scale, shape) { }
            /// <inheritdoc/>
            public override double InverseCDF(double probability)
            {
                Calls?.Add(Name + (probability == .75d ? ":upper" : ":lower"));
                if (ThrowOnUpper && probability == .75d) throw new InvalidOperationException("quartile callback");
                return base.InverseCDF(probability);
            }
        }
    }
}
