using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Probability, support, moments and mutable-dependency regressions for composites.</summary>
    [TestClass]
    public class Test_CompositeRobustness
    {
        /// <summary>The minimum of independent unit exponentials is exactly exponential with rate two.</summary>
        [TestMethod]
        public void CompetingExponentialLogDensityHasNoFabricatedFloor()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[] { new Exponential(0, 1), new Exponential(0, 1) });
            Assert.AreEqual(double.NegativeInfinity, distribution.LogPDF(-1));
            Assert.AreEqual(0, distribution.PDF(-1));
            Assert.AreEqual(Math.Log(2) - 2000, distribution.LogPDF(1000), 3E-12);
            Assert.AreEqual(-2000, distribution.LogCCDF(1000), 3E-12);
        }

        /// <summary>Minima and maxima use the appropriate intersection of finite endpoint bounds.</summary>
        [TestMethod]
        public void CompetingSupportUsesMinMaxOperation()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[] { new Uniform(0, 1), new Uniform(2, 3) });
            Assert.AreEqual(1, distribution.InverseCDF(1));
            Assert.AreEqual(0, distribution.PDF(1.5));
            distribution.MinimumOfRandomVariables = false;
            Assert.AreEqual(2, distribution.InverseCDF(0));
            Assert.AreEqual(0, distribution.PDF(1.5));
        }

        /// <summary>Zero-weight singular components cannot contaminate density, tails, or support.</summary>
        [TestMethod]
        public void MixtureIgnoresInactiveSingularComponents()
        {
            var distribution = new Mixture(new[] { 0d, 1d }, new UnivariateDistributionBase[] { new GammaDistribution(1, .5), new Normal(0, 1) });
            Assert.AreEqual(1 / Math.Sqrt(2 * Math.PI), distribution.PDF(0), 2E-15);
            Assert.AreEqual(-.5 * Math.Log(2 * Math.PI), distribution.LogPDF(0), 2E-15);
            var bounded = new Mixture(new[] { 1d, 0d }, new UnivariateDistributionBase[] { new Uniform(0, 1), new Uniform(2, 3) });
            Assert.AreEqual(1, bounded.Maximum);
            Assert.AreEqual(1, bounded.InverseCDF(1));
        }

        /// <summary>Component central moments preserve translation and respond to mutable weights/parameters.</summary>
        [TestMethod]
        public void MixtureMomentsCombineCentrallyAndRefresh()
        {
            var distribution = new Mixture(new[] { .5, .5 }, new UnivariateDistributionBase[] { new Normal(1E10, 1), new Normal(1E10, 1) });
            Assert.AreEqual(1E10, distribution.Mean);
            Assert.AreEqual(1, distribution.StandardDeviation, 1E-14);
            Assert.AreEqual(0, distribution.Skewness, 1E-14);
            Assert.AreEqual(3, distribution.Kurtosis, 1E-14);
            distribution.Distributions[0].SetParameters(new[] { 0d, 1d });
            distribution.Distributions[1].SetParameters(new[] { 4d, 1d });
            distribution.Weights[0] = .25;
            distribution.Weights[1] = .75;
            Assert.AreEqual(3, distribution.Mean, 1E-14);
            Assert.AreEqual(2, distribution.StandardDeviation, 1E-14);
            var separated = new Mixture(new[] { .9, .1 }, new UnivariateDistributionBase[]
            { new Normal(-1E308, 1E200), new Normal(1E308, 1E200) });
            Assert.AreEqual(-8E307, separated.Mean, 2E292);
            Assert.AreEqual(6E307, separated.StandardDeviation, 3E292);
            Assert.AreEqual(8d / 3, separated.Skewness, 3E-14);
            Assert.AreEqual(73d / 9, separated.Kurtosis, 5E-14);
        }

        /// <summary>Hurdle probability tails remain logarithmic and interval endpoint atoms are included correctly.</summary>
        [TestMethod]
        public void HurdleLogTailsAndIntervalAtomsRemainCorrect()
        {
            var distribution = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { new Exponential(0, 1) }) { IsZeroInflated = true, ZeroWeight = .2 };
            Assert.AreEqual(Math.Log(.8) - 1000, distribution.LogCCDF(1000), 2E-12);
            Assert.AreEqual(Math.Log(.2), distribution.LogLikelihood_Intervals(-1, 0), 2E-14);
            Assert.AreEqual(Math.Log(.8) + Math.Log(1 - Math.Exp(-1)), distribution.LogLikelihood_Intervals(0, 1), 2E-14);
            var remote = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { new Normal(-40, 1) }) { IsZeroInflated = true, ZeroWeight = .2 };
            Assert.AreEqual(Math.Log(.8) - .5 * 40.1 * 40.1 - .5 * Math.Log(2 * Math.PI) + 804.6084420137538, remote.LogPDF(.1), 3E-12);
        }

        /// <summary>Candidate validation checks the supplied vector independently of the current valid object.</summary>
        [TestMethod]
        public void CompositeValidationUsesCandidateParameters()
        {
            var competing = new CompetingRisks(new UnivariateDistributionBase[] { new Normal(), new Normal() });
            Assert.IsNotNull(competing.ValidateParameters(new[] { 0d, -1d, 0d, 1d }, false));
            Assert.IsNotNull(competing.ValidateParameters(new[] { 0d }, false));
            var mixture = new Mixture(new[] { .5, .5 }, new UnivariateDistributionBase[] { new Normal(), new Normal() });
            Assert.IsNotNull(mixture.ValidateParameters(new[] { -1d, 2d, 0d, 1d, 0d, 1d }, false));
            Assert.IsNotNull(mixture.ValidateParameters(new[] { .5, .5, 0d, -1d, 0d, 1d }, false));
        }

        /// <summary>Endpoint density uses the joint one-sided limit even when component factors are infinity times zero.</summary>
        [TestMethod]
        public void CompetingEndpointDensityCombinesPowersBeforeTakingLimit()
        {
            // F_Gamma(1/2,1)(x) ~ 2 sqrt(x/pi), hence d(F squared)/dx -> 4/pi.
            var maximum = new CompetingRisks(new UnivariateDistributionBase[] { new GammaDistribution(1, .5), new GammaDistribution(1, .5) })
            { MinimumOfRandomVariables = false };
            Assert.AreEqual(4 / Math.PI, maximum.PDF(0), 2E-14);
            maximum.Distributions[1].SetParameters(new[] { 1d, .25 });
            Assert.AreEqual(double.PositiveInfinity, maximum.LogPDF(0));
            maximum.Distributions[1].SetParameters(new[] { 1d, .75 });
            Assert.AreEqual(double.NegativeInfinity, maximum.LogPDF(0));
            // The minimum of two GPA(0,1,2) variables has constant density two up to .5.
            var minimum = new CompetingRisks(new UnivariateDistributionBase[] { new GeneralizedPareto(0, 1, 2), new GeneralizedPareto(0, 1, 2) });
            Assert.AreEqual(2, minimum.PDF(.5), 2E-14);
            foreach (var component in new UnivariateDistributionBase[]
            { new GeneralizedLogistic(0, 1, 2), new GeneralizedExtremeValue(0, 1, 2), new KappaFour(0, 1, 2, -1) })
            {
                minimum = new CompetingRisks(new[] { component, component.Clone() });
                Assert.AreEqual(2, minimum.PDF(.5), 2E-14, component.DisplayName);
            }
            var kappa = new KappaFour(0, 1, 1, 2);
            maximum = new CompetingRisks(new UnivariateDistributionBase[] { kappa, kappa.Clone() }) { MinimumOfRandomVariables = false };
            Assert.AreEqual(2, maximum.PDF(.5), 2E-14);
            foreach (var component in new UnivariateDistributionBase[]
            { new PearsonTypeIII(.25, .5, 4), new LogPearsonTypeIII(.25, .5, 4) { Base = Math.E } })
            {
                maximum = new CompetingRisks(new[] { component, component.Clone(), component.Clone(), component.Clone() }) { MinimumOfRandomVariables = false };
                // R 4.4.3: exp(-4*lgamma(1.25)). Both transformed endpoints have unit local scale.
                Assert.AreEqual(1.4815477904878236, maximum.PDF(component.Minimum), 3E-14, component.DisplayName);
            }
            var reflected = new LogPearsonTypeIII(0, 2, -2) { Base = Math.E };
            maximum = new CompetingRisks(new UnivariateDistributionBase[] { reflected, reflected.Clone() }) { MinimumOfRandomVariables = false };
            Assert.AreEqual(Math.Exp(-2), maximum.PDF(0), 2E-14);
        }

        /// <summary>Checked integration covers complete support and refreshes after a min/max change.</summary>
        [TestMethod]
        public void CompetingMomentsUseFullSupportAndRefreshConfiguration()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[] { new Exponential(0, 1), new Exponential(0, 1) });
            Assert.AreEqual(.5, distribution.Mean, 2E-9);
            Assert.AreEqual(.5, distribution.StandardDeviation, 2E-9);
            Assert.AreEqual(2, distribution.Skewness, 2E-8);
            Assert.AreEqual(9, distribution.Kurtosis, 2E-7);
            distribution.MinimumOfRandomVariables = false;
            Assert.AreEqual(1.5, distribution.Mean, 2E-9);
            Assert.AreEqual(Math.Sqrt(1.25), distribution.StandardDeviation, 2E-9);
        }

        /// <summary>Positive conditional quantiles remain accurate for remote tails and tiny physical scales.</summary>
        [TestMethod]
        public void HurdleQuantilesAndCentralMomentsAreScaleStable()
        {
            var half = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { new Normal() }) { IsZeroInflated = true, ZeroWeight = .2 };
            double mean = .8 * Math.Sqrt(2 / Math.PI);
            double variance = .8 - mean * mean;
            double central3 = 1.6 * Math.Sqrt(2 / Math.PI) - 3 * mean * .8 + 2 * mean * mean * mean;
            double central4 = 2.4 - 4 * mean * 1.6 * Math.Sqrt(2 / Math.PI) + 6 * mean * mean * .8 - 3 * Math.Pow(mean, 4);
            Assert.AreEqual(mean, half.Mean, 1E-9);
            Assert.AreEqual(Math.Sqrt(variance), half.StandardDeviation, 1E-9);
            Assert.AreEqual(central3 / Math.Pow(variance, 1.5), half.Skewness, 2E-8);
            Assert.AreEqual(central4 / (variance * variance), half.Kurtosis, 2E-8);
            var remote = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { new Normal(-40, 1) }) { IsZeroInflated = true };
            double q = remote.InverseCDF(.5);
            Assert.AreEqual(Math.Log(.5), remote.LogCCDF(q), 2E-6);
            // GNO has no legacy Normal.Sigma minimum-scale clamp.
            var tiny = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { new GeneralizedNormal(0, 1E-100, 0) }) { IsZeroInflated = true };
            Assert.AreEqual(.6744897501960817, tiny.InverseCDF(.5) / 1E-100, 2E-6);
        }

        /// <summary>Ordinary mixtures solve quantiles in scaled coordinates rather than treating nearby physical endpoints as equal.</summary>
        [TestMethod]
        public void MixtureQuantileRespectsTinyScaleAndNonparametricChildren()
        {
            var tiny = new Mixture(new[] { .5, .5 }, new UnivariateDistributionBase[]
            { new GeneralizedNormal(0, 1E-100, 0), new GeneralizedNormal(1E-100, 1E-100, 0) });
            Assert.AreEqual(.5, tiny.InverseCDF(.5) / 1E-100, 2E-6);
            var empirical = new EmpiricalDistribution(new[] { 1d, 2d, 3d }, new[] { .1, .5, .9 });
            var single = new Mixture(new[] { 1d }, new UnivariateDistributionBase[] { empirical });
            Assert.AreEqual(empirical.InverseCDF(.5), single.InverseCDF(.5));
        }

        /// <summary>The generic fixed-bin moment approximation must not subtract two large raw moments.</summary>
        [TestMethod]
        public void NumericalBinMomentsPreserveTranslation()
        {
            var unit = new Normal(0, 1).CentralMoments(300);
            var translated = new Normal(1E10, 1).CentralMoments(300);
            Assert.AreEqual(unit[1], translated[1], 3E-6);
            Assert.AreEqual(unit[2], translated[2], 3E-6);
            Assert.AreEqual(unit[3], translated[3], 3E-5);
        }

        /// <summary>Independent minimum quantiles use the same dimensionless root accuracy at tiny scales.</summary>
        [TestMethod]
        public void CompetingQuantilePreservesPhysicalScale()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[]
            { new GeneralizedNormal(0, 1E-100, 0), new GeneralizedNormal(1E-100, 1E-100, 0) });
            Assert.AreEqual(-.1725495296281549, distribution.InverseCDF(.5) / 1E-100, 2E-6);
        }

        /// <summary>Identical component laws retain their relative weights even at enormous negative log densities.</summary>
        [TestMethod]
        public void LogResponsibilityNormalizationPreservesWeights()
        {
            var responsibilities = new double[3];
            double logRow = MixtureLogWeights.Normalize(new[] { -5E199, -5E199, double.PositiveInfinity }, new[] { .4, .6, 0d }, responsibilities);
            Assert.AreEqual(-5E199, logRow);
            Assert.AreEqual(.4, responsibilities[0], 2E-15);
            Assert.AreEqual(.6, responsibilities[1], 2E-15);
            Assert.AreEqual(0, responsibilities[2]);
        }
    }
}
