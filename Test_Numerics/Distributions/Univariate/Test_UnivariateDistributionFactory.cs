using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Regression tests for the univariate distribution factory.
    /// </summary>
    [TestClass]
    public class Test_UnivariateDistributionFactory
    {
        /// <summary>
        /// Verifies that every defined enum member is either constructed exactly or rejected explicitly.
        /// </summary>
        [TestMethod]
        public void EveryDefinedDistributionTypeIsHandledExplicitly()
        {
            foreach (UnivariateDistributionType type in Enum.GetValues(typeof(UnivariateDistributionType)))
            {
                if (type == UnivariateDistributionType.CompetingRisks ||
                    type == UnivariateDistributionType.Mixture ||
                    type == UnivariateDistributionType.UserDefined)
                {
                    AssertThrows<NotSupportedException>(
                        () => UnivariateDistributionFactory.CreateDistribution(type));
                    Assert.IsFalse(
                        UnivariateDistributionFactory.TryCreateDistribution(type, out var unsupportedDistribution));
                    Assert.IsNull(unsupportedDistribution);
                }
                else
                {
                    var distribution = UnivariateDistributionFactory.CreateDistribution(type);
                    Assert.AreEqual(type, distribution.Type);
                    Assert.IsTrue(
                        UnivariateDistributionFactory.TryCreateDistribution(type, out var createdDistribution));
                    Assert.IsNotNull(createdDistribution);
                    Assert.AreEqual(type, createdDistribution.Type);
                }
            }
        }

        /// <summary>
        /// Verifies that undefined enum values fail fast rather than silently creating a deterministic distribution.
        /// </summary>
        [TestMethod]
        public void UndefinedDistributionTypeThrowsArgumentOutOfRange()
        {
            AssertThrows<ArgumentOutOfRangeException>(
                () => UnivariateDistributionFactory.CreateDistribution((UnivariateDistributionType)int.MaxValue));
            Assert.IsFalse(UnivariateDistributionFactory.TryCreateDistribution(
                (UnivariateDistributionType)int.MaxValue, out var distribution));
            Assert.IsNull(distribution);
        }

        /// <summary>
        /// Verifies that an action throws the requested exception type on every target framework.
        /// </summary>
        /// <typeparam name="TException">The exception type expected from <paramref name="action"/>.</typeparam>
        /// <param name="action">The action expected to throw.</param>
        private static void AssertThrows<TException>(Action action)
            where TException : Exception
        {
            try
            {
                action();
            }
            catch (TException)
            {
                return;
            }
            catch (Exception exception)
            {
                Assert.Fail("Expected " + typeof(TException).Name + " but received " + exception.GetType().Name + ".");
                return;
            }

            Assert.Fail("Expected " + typeof(TException).Name + ".");
        }
    }
}
