using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Independent tail, moment-coordinate and uncertainty regressions for the normal and Pearson families.
    /// </summary>
    [TestClass]
    public class Test_NormalPearsonRobustness
    {
        /// <summary>Scalar uncertainty remains meaningful when physical covariance entries or intermediate squares overflow.</summary>
        [TestMethod]
        public void UncertaintyAtExtremeScaleCombinesFactorsBeforeSquaring()
        {
            const ParameterEstimationMethod method = ParameterEstimationMethod.MaximumLikelihood;
            Assert.AreEqual(1E308, new Normal(0, 1E155).ParameterCovariance(100, method)[0, 0], 3E292);
            Assert.AreEqual(3E307, new Logistic(0, 1E155).ParameterCovariance(1000, method)[0, 0], 8E291);
            foreach (IStandardError distribution in new IStandardError[]
            { new Normal(0, 1E200), new Logistic(0, 1E200), new PearsonTypeIII(0, 1E200, .5) })
                Assert.AreEqual(double.PositiveInfinity, distribution.QuantileVariance(.5, 100, method));
            var natural = new LnNormal { Mu = 200, Sigma = 20 };
            Assert.AreEqual(4 * Math.Exp(400), natural.QuantileVariance(.5, 100, method), 8E161);
            var logged = new LogPearsonTypeIII(370, 1E-10, 0) { Base = Math.E };
            // Full-family zero skew includes variance of estimated skew: median unit variance is 7/(6*n).
            Assert.AreEqual((7d / 6) * 2.38735282838454E299, logged.QuantileVariance(.5, 100, method), 2E287);
        }

        /// <summary>Normal log tails remain finite beyond ordinary probability underflow.</summary>
        [TestMethod]
        public void Normal_LogTailsAndLargeScale_MatchDefiningDensity()
        {
            var normal = new Normal();
            Assert.AreEqual(-804.6084420137538, normal.LogCDF(-40), 2E-12);
            Assert.AreEqual(-804.6084420137538, normal.LogCCDF(40), 2E-12);
            Assert.AreEqual(-710.11514717537079, new Normal(0, 1E308).LogPDF(0), 2E-12);
            Assert.AreEqual(new Normal().LogPDF(2) - Math.Log(1E308), new Normal(-1E308, 1E308).LogPDF(1E308), 2E-12);
            Assert.AreEqual(0d, normal.CCDF(double.PositiveInfinity));
            Assert.AreEqual(1d, normal.CCDF(double.NegativeInfinity));
        }

        /// <summary>Logistic evaluation uses direct log tails and exact infinite endpoints.</summary>
        [TestMethod]
        public void Logistic_LogTailsAndEndpoints_MatchDefiningLogit()
        {
            var logistic = new Logistic(0, 1);
            Assert.AreEqual(-1000d, logistic.LogPDF(-1000), 1E-12);
            Assert.AreEqual(-1000d, logistic.LogPDF(1000), 1E-12);
            Assert.AreEqual(-1000d, logistic.LogCDF(-1000), 1E-12);
            Assert.AreEqual(-1000d, logistic.LogCCDF(1000), 1E-12);
            Assert.AreEqual(0d, logistic.PDF(double.NegativeInfinity));
            Assert.AreEqual(0d, logistic.PDF(double.PositiveInfinity));
            Assert.AreEqual(Math.Exp(-40), logistic.CCDF(40), 1E-32);
        }

        /// <summary>Signed and exactly centered samples have finite feasible optimizer initialization.</summary>
        [TestMethod]
        public void LocationScaleConstraints_SignedZeroAndLargeSamples_AreFiniteAndFeasible()
        {
            foreach (double[] sample in new[] { new[] { -5d, -4d, -2d, -1d }, new[] { -2d, -1d, 1d, 2d } })
            {
                CheckConstraints(new Normal().GetParameterConstraints(sample));
                CheckConstraints(new Logistic().GetParameterConstraints(sample));
                CheckConstraints(new PearsonTypeIII().GetParameterConstraints(sample));
            }
            CheckConstraints(new Normal().GetParameterConstraints(new[] { -1E200, -5E199, 5E199, 1E200 }));
            CheckConstraints(new LogNormal().GetParameterConstraints(new[] { 0.01d, 0.1d, 10d, 100d }));
            CheckConstraints(new LogPearsonTypeIII().GetParameterConstraints(new[] { 0.01d, 0.1d, 10d, 100d }));
        }

        /// <summary>Malformed samples are rejected before fitting or replacing nonpositive observations.</summary>
        [TestMethod]
        public void SampleValidation_RejectsNonfiniteDegenerateAndNonpositiveSamples()
        {
            IMaximumLikelihoodEstimation[] distributions = { new Normal(), new Logistic(), new LnNormal(), new LogNormal(), new PearsonTypeIII(), new LogPearsonTypeIII() };
            foreach (var distribution in distributions)
            {
                Assert.Throws<ArgumentException>(() => distribution.GetParameterConstraints(new[] { 1d, double.NaN, 2d, 3d }));
                Assert.Throws<ArgumentException>(() => distribution.GetParameterConstraints(new[] { 1d, double.PositiveInfinity, 2d, 3d }));
                Assert.Throws<ArgumentException>(() => distribution.GetParameterConstraints(new[] { 2d, 2d, 2d, 2d }));
                Assert.Throws<ArgumentException>(() => distribution.GetParameterConstraints(new[] { 2d }));
                CheckConstraints(distribution.GetParameterConstraints(new[] { 1d, 2d, 3d, 4d }));
            }
            Assert.Throws<ArgumentException>(() => LnNormal.IndirectMethodOfMoments(new[] { -1d, 1d, 2d, 3d }));
            Assert.Throws<ArgumentException>(() => new LogNormal().IndirectMethodOfMoments(new[] { 0d, 1d, 2d, 3d }));
            Assert.Throws<ArgumentException>(() => new LogPearsonTypeIII().IndirectMethodOfMoments(new[] { 0d, 1d, 2d, 3d }));
        }

        /// <summary>LnNormal public conversion coordinates are physical mean and standard deviation.</summary>
        [TestMethod]
        public void LnNormal_PublicMomentCoordinates_AndTinyVarianceRoundTrip()
        {
            var distribution = new LnNormal(10, 2);
            double[] parameters = distribution.ParametersFromMoments(new[] { 10d, 2d });
            Assert.AreEqual(10d, parameters[0], 1E-13);
            Assert.AreEqual(2d, parameters[1], 1E-13);
            Assert.IsFalse(new LnNormal(-10, 2).ParametersValid);
            var tiny = new LnNormal(1, 1E-10);
            Assert.AreEqual(1E-10, tiny.StandardDeviation, 1E-23);
            Assert.AreEqual(1d, tiny.Mean, 1E-14);
            // d[m^2/sqrt(m^2+s^2)]/dm at (m,s)=(10,2) is exactly 135/(26*sqrt(26)).
            Assert.AreEqual(135d / (26d * Math.Sqrt(26d)), distribution.QuantileGradient(0.5)[0], 2E-15);
            CheckGradient(distribution, 0.83, 3E-7);
        }

        /// <summary>The indirect moment estimator transforms log-sample uncertainty into physical coordinates.</summary>
        [TestMethod]
        public void LnNormal_IndirectMomentCovarianceAndMedianVariance_MatchLogEstimator()
        {
            var distribution = new LnNormal(1, 1);
            double[,] covariance = distribution.ParameterCovariance(100, ParameterEstimationMethod.MethodOfMoments);
            double logVariance = Math.Log(2d);
            Assert.AreEqual((logVariance + logVariance * logVariance / 2d) / 100d, covariance[0, 0], 2E-15);
            Assert.AreEqual((logVariance + 1.5d * logVariance * logVariance) / 100d, covariance[0, 1], 2E-15);
            Assert.AreEqual((logVariance + 4.5d * logVariance * logVariance) / 100d, covariance[1, 1], 2E-14);
            Assert.AreEqual(0.0034657359027997265d, distribution.QuantileVariance(0.5, 100, ParameterEstimationMethod.MethodOfMoments), 3E-15);
            var fitted = new LnNormal();
            fitted.Estimate(new[] { Math.Exp(-Math.Sqrt(1.5d)), 1d, 1d, Math.Exp(Math.Sqrt(1.5d)) }, ParameterEstimationMethod.MethodOfMoments);
            Assert.AreEqual(0d, fitted.Mu, 2E-15);
            Assert.AreEqual(1d, fitted.Sigma, 2E-15);
            Assert.AreEqual(Math.Exp(0.5), fitted.GetParameters[0], 2E-14);
            CheckCovarianceAndVariance(distribution, ParameterEstimationMethod.MethodOfMoments, 0.5);
        }

        /// <summary>The MLE delta transformation includes both public coordinates.</summary>
        [TestMethod]
        public void LnNormal_MleCovariance_MatchesTransformedNormalInformation()
        {
            var distribution = new LnNormal(10, 2);
            double variance = distribution.QuantileVariance(0.8, 100, ParameterEstimationMethod.MaximumLikelihood);
            double z = Normal.StandardZ(0.8);
            double expected = Math.Pow(distribution.InverseCDF(0.8), 2) * Math.Log(1.04) / 100 * (1 + z * z / 2);
            Assert.AreEqual(expected, variance, 2E-13);
            CheckCovarianceAndVariance(distribution, ParameterEstimationMethod.MaximumLikelihood, 0.8);
        }

        /// <summary>LogNormal conversion and mode account for the logarithm base.</summary>
        [TestMethod]
        public void LogNormal_MomentConversionModeAndBaseRoundTrip_MatchLognormal()
        {
            foreach (double logarithmBase in new[] { Math.E, 2d, 10d })
            {
                var distribution = new LogNormal { Base = logarithmBase };
                double[] parameters = distribution.ParametersFromMoments(new[] { 10d, 2d });
                distribution.SetParameters(parameters);
                Assert.AreEqual(10d, distribution.Mean, 1E-12);
                Assert.AreEqual(2d, distribution.StandardDeviation, 1E-12);
                Assert.AreEqual(10d / Math.Pow(1.04, 1.5), distribution.Mode, 1E-12);
                double[] moments = distribution.MomentsFromParameters(parameters);
                Assert.AreEqual(10d, moments[0], 1E-12);
                Assert.AreEqual(2d, moments[1], 1E-12);
                CheckGradient(distribution, 0.83, 5E-7);
                CheckCovarianceAndVariance(distribution, ParameterEstimationMethod.MaximumLikelihood, 0.83);
            }
            Assert.AreEqual(0.0860086348330568, new LogNormal().ParametersFromMoments(new[] { 10d, 2d })[1], 2E-16);
        }

        /// <summary>Standardized lognormal moments depend on shape, even when raw moments overflow.</summary>
        [TestMethod]
        public void LogNormal_ShapeMoments_AreScaleIndependent()
        {
            var ordinary = new LogNormal(0, 0.2);
            var huge = new LogNormal(200, 0.2);
            Assert.AreEqual(ordinary.Skewness, huge.Skewness, 2E-14);
            Assert.AreEqual(ordinary.Kurtosis, huge.Kurtosis, 2E-14);
            Assert.IsTrue(IsFinite(huge.StandardDeviation));
            var tiny = new LogNormal(0, 1E-10) { Base = Math.E };
            Assert.AreEqual(1E-10, tiny.StandardDeviation, 1E-23);
        }

        /// <summary>Base one, nonfinite bases and bases below one are invalid.</summary>
        [TestMethod]
        public void LogarithmBases_RequireFiniteValuesGreaterThanOne()
        {
            foreach (double value in new[] { 0d, 1d, double.NaN, double.PositiveInfinity })
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => new LogNormal { Base = value });
                Assert.Throws<ArgumentOutOfRangeException>(() => new LogPearsonTypeIII { Base = value });
            }
        }

        /// <summary>Negative skew inverts the gamma upper tail directly.</summary>
        [TestMethod]
        public void Pearson_NegativeSkewDeepTail_MatchesReflectedExponential()
        {
            var distribution = new PearsonTypeIII(0, 1, -2);
            double quantile = distribution.InverseCDF(1E-20);
            Assert.AreEqual(-45.051701859880914, quantile, 2E-12);
            Assert.AreEqual(Math.Log(1E-20), distribution.LogCDF(quantile), 2E-12);
            Assert.AreEqual(1E-20, distribution.CDF(quantile), 2E-32);
            var reflected = new PearsonTypeIII(0, 1, 2);
            Assert.AreEqual(distribution.LogCDF(quantile), reflected.LogCCDF(-quantile), 2E-12);
        }

        /// <summary>Gamma density endpoint limits cover shape one and shape below one.</summary>
        [TestMethod]
        public void Pearson_DensityEndpointAndModes_MatchGammaLimits()
        {
            Assert.AreEqual(0d, new PearsonTypeIII(1, 1, 2).LogPDF(0), 0d);
            Assert.AreEqual(1d, new PearsonTypeIII(1, 1, 2).PDF(0), 0d);
            foreach (double skew in new[] { -3d, 3d })
            {
                var distribution = new PearsonTypeIII(1, 2, skew);
                Assert.AreEqual(distribution.Xi, distribution.Mode, 0d);
                Assert.AreEqual(double.PositiveInfinity, distribution.LogPDF(distribution.Xi));
            }
        }

        /// <summary>Small nonzero skew retains its signed first-order quantile correction.</summary>
        [TestMethod]
        public void Pearson_SmallSkew_PreservesShapeAndNormalLimit()
        {
            foreach (double skew in new[] { -1E-5, 1E-5 })
            {
                var distribution = new PearsonTypeIII(0, 1, skew);
                Assert.AreEqual(-skew / 6d, distribution.InverseCDF(0.5), 3E-11);
                Assert.AreEqual(1d / 6d, (0.5 - distribution.CDF(0)) / (-skew * Math.Exp(-0.5 * Math.Log(2 * Math.PI))), 3E-5);
                Assert.AreEqual(skew > 0 ? -2d / skew : double.NegativeInfinity, distribution.Minimum);
            }
        }

        /// <summary>Pearson uncertainty is expressed in the same public coordinates as quantiles.</summary>
        [TestMethod]
        public void Pearson_PublicGradientCovarianceAndNormalLimit_AreConsistent()
        {
            foreach (double skew in new[] { -1.2d, -0.2d, 0d, 0.2d, 1.2d })
            {
                var distribution = new PearsonTypeIII(2, 3, skew);
                CheckGradient(distribution, 0.83, 2E-5);
                CheckCovarianceAndVariance(distribution, ParameterEstimationMethod.MaximumLikelihood, 0.83);
                CheckCovarianceAndVariance(distribution, ParameterEstimationMethod.MethodOfMoments, 0.83);
                double[] publicGradient = distribution.QuantileGradient(0.83);
                double[] momentGradient = distribution.QuantileGradientForMoments(0.83);
                for (int i = 0; i < 3; i++) Assert.AreEqual(publicGradient[i], momentGradient[i], 0d);
            }
            double[,] covariance = new PearsonTypeIII(2, 3, 0).ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(0.09, covariance[0, 0], 1E-15);
            Assert.AreEqual(0.045, covariance[1, 1], 1E-15);
            Assert.AreEqual(0.06, covariance[2, 2], 1E-15);
            foreach (double skew in new[] { -2d, -Math.Sqrt(2), Math.Sqrt(2), 2d })
                Assert.Throws<ArgumentOutOfRangeException>(() => new PearsonTypeIII(0, 1, skew).ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood));
        }

        /// <summary>The negative-skew exponential transformed with base e is uniform on (0,e).</summary>
        [TestMethod]
        public void LogPearson_UniformCase_AndTransformedEndpoints()
        {
            var distribution = new LogPearsonTypeIII(0, 1, -2) { Base = Math.E };
            Assert.AreEqual(-1d, distribution.LogPDF(0), 1E-14);
            Assert.AreEqual(-1d, distribution.LogPDF(1), 1E-14);
            Assert.AreEqual(-1d, distribution.LogPDF(Math.E), 1E-14);
            Assert.AreEqual(Math.E / 2, distribution.Mean, 2E-14);
            Assert.AreEqual(Math.E / Math.Sqrt(12), distribution.StandardDeviation, 2E-14);
            Assert.AreEqual(0d, distribution.Skewness, 2E-13);
            Assert.AreEqual(1.8d, distribution.Kurtosis, 2E-12);
            Assert.AreEqual(Math.E * 1E-20, distribution.InverseCDF(1E-20), 2E-31);
        }

        /// <summary>Each LP3 moment has its own existence threshold.</summary>
        [TestMethod]
        public void LogPearson_MomentExistence_IsCheckedForEachOrder()
        {
            var distribution = new LogPearsonTypeIII(0, 0.8, 1) { Base = Math.E };
            Assert.AreEqual(1.55784d, distribution.Mean, 1E-5);
            Assert.AreEqual(4.80099d, distribution.StandardDeviation, 1E-5);
            Assert.AreEqual(double.PositiveInfinity, distribution.Skewness);
            Assert.AreEqual(double.PositiveInfinity, distribution.Kurtosis);
            var normalLimit = new LogPearsonTypeIII(0.3, 0.2, 0) { Base = Math.E };
            Assert.AreEqual(Math.Exp(0.26), normalLimit.Mode, 2E-14);
            var skewed = new LogPearsonTypeIII(0.3, 0.2, 0.8) { Base = Math.E };
            Assert.AreEqual(Math.Exp(skewed.Xi + (skewed.Alpha - 1) * skewed.Beta / (1 + skewed.Beta)), skewed.Mode, 2E-14);
        }

        /// <summary>LP3 quantile derivatives transform actual quantiles with one Jacobian factor.</summary>
        [TestMethod]
        public void LogPearson_PublicCoordinatesBaseRoundTripAndVariance_AreConsistent()
        {
            foreach (double skew in new[] { -0.8d, 0d, 0.8d })
            {
                var distribution = new LogPearsonTypeIII(0.3, 0.2, skew) { Base = Math.E };
                CheckGradient(distribution, 0.83, 3E-5);
                CheckCovarianceAndVariance(distribution, ParameterEstimationMethod.MaximumLikelihood, 0.83);
                CheckCovarianceAndVariance(distribution, ParameterEstimationMethod.MethodOfMoments, 0.83);
                double[] moments = distribution.MomentsFromParameters(distribution.GetParameters);
                Assert.AreEqual(distribution.Mean, moments[0], 1E-13);
                Assert.AreEqual(distribution.StandardDeviation, moments[1], 1E-13);
            }
        }

        /// <summary>Uncertainty methods reject endpoint probabilities and invalid sample sizes.</summary>
        [TestMethod]
        public void Uncertainty_RejectsNoninteriorProbabilitiesAndInvalidSampleSizes()
        {
            IStandardError[] distributions = { new Normal(), new Logistic(), new LnNormal(), new LogNormal(), new PearsonTypeIII(), new LogPearsonTypeIII() };
            foreach (IStandardError distribution in distributions)
            {
                foreach (double probability in new[] { 0d, 1d, double.NaN, double.PositiveInfinity, -0.1d })
                {
                    Assert.Throws<ArgumentOutOfRangeException>(() => distribution.QuantileGradient(probability));
                    Assert.Throws<ArgumentOutOfRangeException>(() => distribution.QuantileVariance(probability, 100, ParameterEstimationMethod.MaximumLikelihood));
                }
                foreach (int sampleSize in new[] { -1, 0 })
                    Assert.Throws<ArgumentOutOfRangeException>(() => distribution.ParameterCovariance(sampleSize, ParameterEstimationMethod.MaximumLikelihood));
            }
        }

        /// <summary>Frozen R stats quantiles, numerical derivatives and independently inverted Fisher matrices agree.</summary>
        [TestMethod]
        public void FrozenRStatsOracle_MatchesPublicQuantilesGradientsAndCovariances()
        {
            int count = 0;
            foreach (Dictionary<string, string> row in DistributionOracle.Read("normal-pearson.csv"))
            {
                double Number(string key) => DistributionOracle.Number(row[key]);
                UnivariateDistributionBase distribution = row["family"] switch
                {
                    "PearsonTypeIII" => new PearsonTypeIII(Number("mu"), Number("sigma"), Number("gamma")),
                    "LogPearsonTypeIII" => new LogPearsonTypeIII(Number("mu"), Number("sigma"), Number("gamma")) { Base = Number("base") },
                    "LogNormal" => new LogNormal(Number("mu"), Number("sigma")) { Base = Number("base") },
                    "LnNormal" => new LnNormal(Number("mu"), Number("sigma")),
                    _ => throw new InvalidOperationException("Unknown oracle distribution.")
                };
                var uncertainty = (IStandardError)distribution;
                double probability = Number("probability"), expectedQuantile = Number("quantile");
                Assert.AreEqual(expectedQuantile, distribution.InverseCDF(probability), 3E-11 * Math.Max(1, Math.Abs(expectedQuantile)));
                double[] gradient = uncertainty.QuantileGradient(probability);
                string[] names = { "gradient_mu", "gradient_sigma", "gradient_gamma" };
                for (int i = 0; i < gradient.Length; i++)
                    Assert.AreEqual(Number(names[i]), gradient[i], 3E-6 * Math.Max(1, Math.Abs(gradient[i])), $"{row["family"]} g={row["gamma"]} p={probability}, gradient {i}");
                double[,] covariance = uncertainty.ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood);
                for (int i = 0; i < gradient.Length; i++)
                    for (int j = i; j < gradient.Length; j++)
                        Assert.AreEqual(Number($"cov{i + 1}{j + 1}"), covariance[i, j], 2E-10 * Math.Max(1, Math.Abs(covariance[i, j])));
                double expectedVariance = Number("variance_mle");
                Assert.AreEqual(expectedVariance, uncertainty.QuantileVariance(probability, 100, ParameterEstimationMethod.MaximumLikelihood), 3E-6 * Math.Max(1, expectedVariance));
                count++;
            }
            Assert.AreEqual(50, count);
        }

        /// <summary>Confidence interval input validation runs before simulation or quantile calculation.</summary>
        [TestMethod]
        public void ConfidenceIntervals_RejectInvalidProbabilitiesAndSampleSizes()
        {
            var normal = new Normal();
            foreach (double probability in new[] { 0d, 1d, double.NaN })
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => normal.NormalConfidenceIntervals(10, new[] { probability }, new[] { 0.5d }));
                Assert.Throws<ArgumentOutOfRangeException>(() => normal.NoncentralTConfidenceIntervals(10, new[] { 0.5d }, new[] { probability }));
                Assert.Throws<ArgumentOutOfRangeException>(() => normal.MonteCarloConfidenceIntervals(10, 1, new[] { probability }, new[] { 0.5d }));
                Assert.Throws<ArgumentOutOfRangeException>(() => new LogNormal().MonteCarloConfidenceIntervals(10, 1, new[] { 0.5d }, new[] { probability }));
                Assert.Throws<ArgumentOutOfRangeException>(() => normal.ExpectedProbability(10, probability));
            }
            Assert.Throws<ArgumentOutOfRangeException>(() => normal.NormalConfidenceIntervals(0, new[] { 0.5d }, new[] { 0.5d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => normal.NoncentralTConfidenceIntervals(1, new[] { 0.5d }, new[] { 0.5d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => normal.MonteCarloConfidenceIntervals(10, 0, new[] { 0.5d }, new[] { 0.5d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new LogNormal().MonteCarloConfidenceIntervals(1, 1, new[] { 0.5d }, new[] { 0.5d }));
        }

        /// <summary>Finite covariance scales are retained when squaring the original scale would overflow.</summary>
        [TestMethod]
        public void LargeScaleCovariance_DividesBeforeSquaring()
        {
            Assert.AreEqual(1E308, new Normal(0, 1E155).ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood)[0, 0], 3E292);
            Assert.AreEqual(3E306, new Logistic(0, 1E155).ParameterCovariance(10000, ParameterEstimationMethod.MaximumLikelihood)[0, 0], 1E291);
        }

        /// <summary>Large LP3 standardized moments retain finite ratios and classify genuine overflow.</summary>
        [TestMethod]
        public void LogPearson_LargeShapeMoments_NormalizeBeforeExponentiating()
        {
            var wide = new LogPearsonTypeIII(0, 16, 0) { Base = Math.E };
            Assert.AreEqual(5.8759900382892355E166, wide.Skewness, 2E153);
            Assert.AreEqual(double.PositiveInfinity, wide.Kurtosis);
            var moderate = new LogPearsonTypeIII(0, 12, 0) { Base = Math.E };
            Assert.AreEqual(1.4243659274306933E250, moderate.Kurtosis, 5E236);
            // Frozen gamma-MGF calculations in normal-pearson.R cover both nonzero skew signs.
            Assert.AreEqual(2.6371094491313571E162, new LogPearsonTypeIII(0, 16, -0.001) { Base = Math.E }.Skewness, 6E151);
            Assert.AreEqual(2.0889174692524114E171, new LogPearsonTypeIII(0, 16, 0.001) { Base = Math.E }.Skewness, 5E160);
            Assert.AreEqual(1.8770290525950149E244, new LogPearsonTypeIII(0, 12, -0.001) { Base = Math.E }.Kurtosis, 4E233);
            Assert.AreEqual(1.9321207191996086E256, new LogPearsonTypeIII(0, 12, 0.001) { Base = Math.E }.Kurtosis, 4E245);
        }

        /// <summary>Lognormal covariance and quantile variance do not square an unscaled large factor.</summary>
        [TestMethod]
        public void LogNormal_Uncertainty_PreservesFiniteScaledSquares()
        {
            foreach (ParameterEstimationMethod method in new[] { ParameterEstimationMethod.MethodOfMoments, ParameterEstimationMethod.MaximumLikelihood })
            {
                Assert.AreEqual(1E308, new LogNormal(0, 1E155).ParameterCovariance(100, method)[0, 0], 3E292);
                var distribution = new LogNormal(370, 1E-10) { Base = Math.E };
                Assert.AreEqual(2.38735282838454E299, distribution.QuantileVariance(0.5, 100, method), 2E286);
            }
        }

        private static void CheckConstraints(Tuple<double[], double[], double[]> constraints)
        {
            for (int i = 0; i < constraints.Item1.Length; i++)
            {
                Assert.IsTrue(IsFinite(constraints.Item1[i]) && IsFinite(constraints.Item2[i]) && IsFinite(constraints.Item3[i]));
                Assert.IsGreaterThan(constraints.Item2[i], constraints.Item1[i]);
                Assert.IsLessThan(constraints.Item3[i], constraints.Item1[i]);
            }
        }

        private static bool IsFinite(double value) => !double.IsNaN(value) && !double.IsInfinity(value);

        private static void CheckGradient(UnivariateDistributionBase distribution, double probability, double tolerance)
        {
            var uncertainty = (IStandardError)distribution;
            double[] gradient = uncertainty.QuantileGradient(probability);
            double[] parameters = distribution.GetParameters;
            for (int i = 0; i < parameters.Length; i++)
            {
                double step = Math.Max(1, Math.Abs(parameters[i])) * 1E-5;
                var plus = distribution.Clone();
                var minus = distribution.Clone();
                double[] upper = (double[])parameters.Clone();
                double[] lower = (double[])parameters.Clone();
                upper[i] += step;
                lower[i] -= step;
                plus.SetParameters(upper);
                minus.SetParameters(lower);
                double difference = (plus.InverseCDF(probability) - minus.InverseCDF(probability)) / (2 * step);
                Assert.AreEqual(difference, gradient[i], tolerance * Math.Max(1, Math.Abs(difference)), $"{distribution.DisplayName}, parameter {i}");
            }
        }

        private static void CheckCovarianceAndVariance(UnivariateDistributionBase distribution, ParameterEstimationMethod method, double probability)
        {
            var uncertainty = (IStandardError)distribution;
            double[] gradient = uncertainty.QuantileGradient(probability);
            double[,] covariance = uncertainty.ParameterCovariance(100, method);
            double expected = 0;
            for (int i = 0; i < gradient.Length; i++)
            {
                Assert.IsTrue(IsFinite(covariance[i, i]) && covariance[i, i] > 0);
                for (int j = 0; j < gradient.Length; j++)
                {
                    Assert.AreEqual(covariance[i, j], covariance[j, i], 1E-14);
                    expected += gradient[i] * covariance[i, j] * gradient[j];
                }
            }
            Assert.IsGreaterThan(0d, expected);
            Assert.AreEqual(expected, uncertainty.QuantileVariance(probability, 100, method), 2E-12 * Math.Max(1, expected));
            double[] probabilities = gradient.Length == 2 ? new[] { 0.2, 0.8 } : new[] { 0.2, 0.5, 0.8 };
            double[,] jacobian = uncertainty.QuantileJacobian(probabilities, out double determinant);
            Assert.IsTrue(IsFinite(determinant));
            for (int row = 0; row < probabilities.Length; row++)
            {
                double[] rowGradient = uncertainty.QuantileGradient(probabilities[row]);
                for (int column = 0; column < gradient.Length; column++)
                    Assert.AreEqual(rowGradient[column], jacobian[row, column], 2E-12 * Math.Max(1, Math.Abs(rowGradient[column])));
            }
        }
    }
}
