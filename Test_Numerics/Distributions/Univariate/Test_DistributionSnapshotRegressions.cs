using System;
using System.Collections.Generic;
using System.Reflection;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Regression tests for the bitwise configuration snapshot behind the composite
    /// configuration and validation caches.
    /// </summary>
    /// <remarks>
    /// The snapshot has two switches that must stay in lockstep: the capture walk and the
    /// inline compare walk. The round-trip test guards them against drifting apart, and the
    /// reflection sweep guards the per-family scalar lists against new mutable state being
    /// added to a supported family without extending its snapshot arm.
    /// </remarks>
    [TestClass]
    public class Test_DistributionSnapshotRegressions
    {
        /// <summary>One representative, validly parameterized instance per supported family.</summary>
        /// <returns>The instances the snapshot must capture.</returns>
        private static UnivariateDistributionBase[] SupportedInstances()
        {
            return
            [
                new Normal(10, 2),
                new LogNormal() { Mu = 3, Sigma = 0.5 },
                new LnNormal() { Mu = 3, Sigma = 0.5 },
                new LogPearsonTypeIII(3, 0.4, 0.2),
                new PearsonTypeIII(100, 20, 0.4),
                new GammaDistribution(10, 4),
                new Weibull(8, 2.5),
                new Gumbel(100, 12),
                new GeneralizedExtremeValue(100, 12, 0.1),
                new GeneralizedLogistic(100, 12, 0.1),
                new GeneralizedNormal(100, 12, 0.1),
                new GeneralizedPareto(100, 12, 0.1),
                new Exponential(50, 10),
                new KappaFour(100, 12, 0.1, 0.2),
                new Uniform(2, 12),
                new Triangular(2, 5, 12),
                new Pert(2, 5, 12),
                new Deterministic(7),
                new Logistic(100, 12),
                new Cauchy(100, 12),
            ];
        }

        /// <summary>
        /// Capture succeeds for every supported family and the compare walk accepts the
        /// unchanged instance - the round trip that pins the two switches together.
        /// </summary>
        [TestMethod]
        public void CaptureAndCompare_RoundTripEveryFamily()
        {
            foreach (var distribution in SupportedInstances())
            {
                var snapshot = DistributionSnapshot.TryCapture(distribution);
                Assert.IsNotNull(snapshot, distribution.GetType().Name + " must be capturable.");
                Assert.IsTrue(snapshot.Matches(distribution), distribution.GetType().Name + " must match its own capture.");
            }
        }

        /// <summary>
        /// Mutating any public settable double, int, bool, or enum property on a supported
        /// family flips the compare to a mismatch - the completeness guard for the per-family
        /// scalar lists.
        /// </summary>
        [TestMethod]
        public void Capture_EveryPublicSettableScalarFlipsTheMatch()
        {
            foreach (var distribution in SupportedInstances())
            {
                foreach (var property in MutableScalarProperties(distribution.GetType()))
                {
                    var snapshot = DistributionSnapshot.TryCapture(distribution);
                    Assert.IsNotNull(snapshot, distribution.GetType().Name + " must be capturable.");
                    object original = property.GetValue(distribution)!;
                    object mutated = Mutate(property.PropertyType, original);
                    property.SetValue(distribution, mutated);
                    try
                    {
                        Assert.IsFalse(snapshot.Matches(distribution),
                            distribution.GetType().Name + "." + property.Name + " must be part of the snapshot.");
                    }
                    finally
                    {
                        property.SetValue(distribution, original);
                    }
                    Assert.IsTrue(snapshot.Matches(distribution),
                        distribution.GetType().Name + "." + property.Name + " must match again after restoration.");
                }
            }
        }

        /// <summary>
        /// Composite trees capture recursively: mutating a nested grandchild parameter flips
        /// the root compare, and restoring it restores the match.
        /// </summary>
        [TestMethod]
        public void Capture_NestedCompositeGrandchildMutationFlipsTheRoot()
        {
            var grandchild = new Weibull(8, 2.5);
            var nested = new Mixture([0.4, 0.6], new UnivariateDistributionBase[] { grandchild, new LnNormal() { Mu = 3, Sigma = 0.5 } });
            var root = new CompetingRisks(new UnivariateDistributionBase[] { nested, new Weibull(6, 3) });

            var snapshot = DistributionSnapshot.TryCapture(root);
            Assert.IsNotNull(snapshot);
            Assert.IsTrue(snapshot.Matches(root));

            double original = grandchild.Kappa;
            grandchild.Kappa = original + 0.25;
            Assert.IsFalse(snapshot.Matches(root), "A nested grandchild mutation must flip the root compare.");
            grandchild.Kappa = original;
            Assert.IsTrue(snapshot.Matches(root), "Restoring the grandchild must restore the root match.");
        }

        /// <summary>
        /// A derived component defeats capture entirely, keeping the generic canonical-string
        /// path for owners whose components may carry live callbacks.
        /// </summary>
        [TestMethod]
        public void Capture_DerivedComponentDefeatsCapture()
        {
            var root = new Mixture([0.5, 0.5], new UnivariateDistributionBase[] { new DerivedWeibull(8, 2.5), new Weibull(6, 3) });
            Assert.IsNull(DistributionSnapshot.TryCapture(root), "A derived component must defeat capture.");
        }

        /// <summary>
        /// The competing-risks correlation matrix participates in the compare: replacing the
        /// matrix with a bitwise-identical clone still matches, while changing one entry flips it.
        /// </summary>
        [TestMethod]
        public void Capture_CorrelationMatrixEntriesParticipate()
        {
            var root = new CompetingRisks(new UnivariateDistributionBase[] { new Weibull(8, 2.5), new Weibull(6, 3) })
            {
                Dependency = Probability.DependencyType.CorrelationMatrix,
                CorrelationMatrix = new[,] { { 1d, 0.5d }, { 0.5d, 1d } },
            };
            var snapshot = DistributionSnapshot.TryCapture(root);
            Assert.IsNotNull(snapshot);
            Assert.IsTrue(snapshot.Matches(root));

            root.CorrelationMatrix = new[,] { { 1d, 0.5d }, { 0.5d, 1d } };
            Assert.IsTrue(snapshot.Matches(root), "A bitwise-identical matrix clone must still match.");
            root.CorrelationMatrix = new[,] { { 1d, 0.6d }, { 0.6d, 1d } };
            Assert.IsFalse(snapshot.Matches(root), "A changed correlation entry must flip the compare.");
        }

        /// <summary>A minimal Weibull subclass for the derived-component refusal.</summary>
        private sealed class DerivedWeibull : Weibull
        {
            /// <summary>Initializes the derived test distribution.</summary>
            /// <param name="lambda">The scale.</param>
            /// <param name="kappa">The shape.</param>
            internal DerivedWeibull(double lambda, double kappa) : base(lambda, kappa) { }
        }

        /// <summary>Enumerates the public settable scalar-like properties of a family.</summary>
        /// <param name="type">The family type.</param>
        /// <returns>The properties whose mutation must flip the snapshot compare.</returns>
        private static IEnumerable<PropertyInfo> MutableScalarProperties(Type type)
        {
            foreach (var property in type.GetProperties(BindingFlags.Public | BindingFlags.Instance))
            {
                if (property.SetMethod is null || !property.SetMethod.IsPublic) continue;
                var propertyType = property.PropertyType;
                if (propertyType != typeof(double) && propertyType != typeof(int)
                    && propertyType != typeof(bool) && !propertyType.IsEnum) continue;
                yield return property;
            }
        }

        /// <summary>Produces a valid, distinct replacement value for a scalar property.</summary>
        /// <param name="type">The property type.</param>
        /// <param name="original">The original value.</param>
        /// <returns>A value guaranteed to differ from the original.</returns>
        private static object Mutate(Type type, object original)
        {
            if (type == typeof(double))
            {
                double value = (double)original;
                return value * 1.25 + 0.375;
            }
            if (type == typeof(int)) return (int)original + 1;
            if (type == typeof(bool)) return !(bool)original;
            var values = Enum.GetValues(type);
            foreach (var candidate in values)
                if (!candidate.Equals(original)) return candidate;
            return original;
        }
    }
}
