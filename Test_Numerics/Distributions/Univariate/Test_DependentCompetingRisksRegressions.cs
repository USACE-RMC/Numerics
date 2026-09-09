using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics;

namespace Distributions.Univariate
{
    /// <summary>Protects dependent candidate rejection and finite-support density behavior.</summary>
    [TestClass]
    public class Test_DependentCompetingRisksRegressions
    {
        /// <summary>A captured correlated-Weibull candidate remains rejectable without aborting estimation.</summary>
        [TestMethod]
        public void CorrelatedMinimum_UnresolvedCandidateReturnsNegativeInfinity()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Weibull(54.676968161948025d, 43.16940511169378d),
                new Weibull(55.23342531523667d, 47.44412398606073d)
            })
            {
                MinimumOfRandomVariables = true,
                Dependency = Probability.DependencyType.CorrelationMatrix,
                CorrelationMatrix = new[,] { { 1d, .6d }, { .6d, 1d } }
            };

            // Captured from the unchanged BestFit correlated-minimum recovery fixture.
            // The CDF derivative is negative at this rounded upper-tail plateau.
            const double observation = 58.99784456776834d;
            Assert.AreEqual(double.NegativeInfinity, distribution.LogPDF(observation));
            Assert.AreEqual(0d, distribution.PDF(observation));
        }

        /// <summary>A perfectly dependent unit-uniform composite retains unit density at and near its endpoints.</summary>
        [TestMethod]
        public void FiniteSupport_EndpointsAndAdjacentInteriorRetainOneSidedDensity()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Uniform(0d, 1d),
                new Uniform(0d, 1d)
            })
            {
                Dependency = Probability.DependencyType.PerfectlyPositive
            };

            double lowerAdjacent = .25d * NumericalDerivative.CalculateStepSize(0d);
            double upperAdjacent = 1d - .25d * NumericalDerivative.CalculateStepSize(1d);
            foreach (double observation in new[] { 0d, lowerAdjacent, upperAdjacent, 1d })
            {
                Assert.AreEqual(1d, distribution.PDF(observation), 0d, $"x={observation:G17}");
                Assert.AreEqual(0d, distribution.LogPDF(observation), 0d, $"x={observation:G17}");
            }
        }
    }
}
