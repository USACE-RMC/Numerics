using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Independent probability, moment and local-MLE uncertainty contracts for Hosking's generalized families.</summary>
    [TestClass]
    public class Test_GeneralizedRobustness
    {
        /// <summary>The frozen R values distinguish exact nonzero shapes from the limiting distribution.</summary>
        [TestMethod]
        public void GeneralizedMoments_MatchIndependentRValues()
        {
            int count = 0;
            foreach (var row in DistributionOracle.Read("generalized-fisher.csv").Where(r => r["quantity"] == "moment"))
            {
                var distribution = Create(row);
                double[] moments = { distribution.Mean, distribution.StandardDeviation, distribution.Skewness, distribution.Kurtosis };
                int index = int.Parse(row["row"]) - 1;
                double expected = Number(row, "value");
                Assert.AreEqual(expected, moments[index], Number(row, "absolute_error_estimate"), Context(row));
                count++;
            }
            Assert.AreEqual(88, count);
        }

        /// <summary>Large-shape lognormal moments preserve finite products with very small scales.</summary>
        [TestMethod]
        public void GeneralizedNormal_MomentsAndMode_PreserveRepresentableProducts()
        {
            var distribution = new GeneralizedNormal(0, 1E-200, 30);
            Assert.AreEqual(-Math.Exp(Math.Log(1E-200) + 450 - Math.Log(30)), distribution.Mean, 2E-19);
            double expectedSd = Math.Exp(Math.Log(1E-200) + 900 - Math.Log(30));
            Assert.AreEqual(expectedSd, distribution.StandardDeviation, 2E-13 * expectedSd);
            Assert.AreEqual(1E-200 / 30, distribution.Mode, 1E-215);
            Assert.IsTrue(double.IsNegativeInfinity(distribution.Skewness));
            Assert.IsTrue(double.IsPositiveInfinity(distribution.Kurtosis));
        }

        /// <summary>Direct logarithms and survival probabilities remain informative after rounding and underflow.</summary>
        [TestMethod]
        public void GeneralizedProbabilities_ZeroShape_LogTailsAndInfiniteEndpoints()
        {
            var normal = new GeneralizedNormal(0, 1, 0);
            Assert.AreEqual(-804.6084420137538, normal.LogCDF(-40), 2E-12);
            Assert.AreEqual(-804.6084420137538, normal.LogCCDF(40), 2E-12);
            Assert.AreEqual(-800.9189385332047, normal.LogPDF(40), 2E-12);
            var logistic = new GeneralizedLogistic(0, 1, 0);
            Assert.AreEqual(-1000d, logistic.LogCDF(-1000), 1E-12);
            Assert.AreEqual(-1000d, logistic.LogCCDF(1000), 1E-12);
            Assert.AreEqual(-1000d, logistic.LogPDF(-1000), 1E-12);
            Assert.AreEqual(Math.Exp(-40), logistic.CCDF(40), 1E-32);
            foreach (var distribution in new UnivariateDistributionBase[] { normal, logistic })
            {
                Assert.AreEqual(0d, distribution.PDF(double.NegativeInfinity));
                Assert.AreEqual(0d, distribution.PDF(double.PositiveInfinity));
                Assert.AreEqual(0d, distribution.CDF(double.NegativeInfinity));
                Assert.AreEqual(1d, distribution.CDF(double.PositiveInfinity));
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.InverseCDF(double.NaN));
            }
        }

        /// <summary>The latent shape transform remains finite when physical standardization overflows.</summary>
        [TestMethod]
        [DataRow("GNO", -1, -1072.0411870643629, -1531.7204579917529, 0d)]
        [DataRow("GNO", 1, -1072.0411870643629, -1531.7204579917529, 0d)]
        [DataRow("GLO", -1, -46.201488473558619, -509.71423934592178, 4.3044582966584941E-222)]
        [DataRow("GLO", 1, -46.201488473558619, -509.71423934592178, 4.3044582966584941E-222)]
        public void GeneralizedProbabilities_AffineOverflow_PreservesFiniteTransformedTail(
            string family, int direction, double logTail, double logDensity, double density)
        {
            // Base R 4.4.3 pnorm/dnorm and defining logistic formulas; see generalized-affine-overflow.R/.md.
            // The exact transform has |z|=46.201488473558619 although |x/alpha| exceeds binary64 range.
            UnivariateDistributionBase distribution = family == "GNO"
                ? new GeneralizedNormal(0, 1E-200, -20 * direction)
                : new GeneralizedLogistic(0, 1E-200, -20 * direction);
            double x = direction * 1E200;
            double actualLogTail = direction < 0 ? distribution.LogCDF(x) : distribution.LogCCDF(x);
            double actualTail = direction < 0 ? distribution.CDF(x) : distribution.CCDF(x);
            Assert.AreEqual(logTail, actualLogTail, 5E-12);
            Assert.AreEqual(logDensity, distribution.LogPDF(x), 5E-12);
            Assert.AreEqual(Math.Exp(logTail), actualTail, Math.Exp(logTail) * 3E-12);
            Assert.AreEqual(density, distribution.PDF(x), density * 3E-12);
        }

        /// <summary>Nonzero tiny shapes retain their finite support and first-order quantile correction.</summary>
        [TestMethod]
        public void GeneralizedProbabilities_NearZeroShape_IsNotFlattened()
        {
            foreach (double k in new[] { -1E-6, 1E-6 })
            {
                var normal = new GeneralizedNormal(0, 1, k);
                var logistic = new GeneralizedLogistic(0, 1, k);
                foreach (var distribution in new UnivariateDistributionBase[] { normal, logistic })
                {
                    Assert.AreEqual(1 / k, k < 0 ? distribution.Minimum : distribution.Maximum);
                    double z = distribution is GeneralizedNormal ? 1 : Math.Log(9d);
                    double p = distribution is GeneralizedNormal ? 0.8413447460685429 : 0.9;
                    double expected = z - k * z * z / 2 + k * k * z * z * z / 6;
                    Assert.AreEqual(expected, distribution.InverseCDF(p), 3E-14);
                    Assert.AreEqual(p, distribution.CDF(expected), 3E-14);
                }
            }
        }

        /// <summary>The GLO stationary density point includes its exact boundary transition.</summary>
        [TestMethod]
        public void GeneralizedLogistic_ModeAndEndpointDensity_MatchAnalyticalLimits()
        {
            foreach (double k in new[] { -2d, -1d, -0.2, 0d, 0.2, 1d, 2d })
            {
                var distribution = new GeneralizedLogistic(2, 3, k);
                double expectedMode = Math.Abs(k) >= 1 ? 2 + 3 / k : k == 0 ? 2
                    : 2 - 3 * (Math.Exp(-k * Math.Log((1 + k) / (1 - k))) - 1) / k;
                Assert.AreEqual(expectedMode, distribution.Mode, 3E-14);
                if (k == 0) continue;
                double endpoint = k < 0 ? distribution.Minimum : distribution.Maximum;
                double expectedDensity = Math.Abs(k) < 1 ? 0 : Math.Abs(k) == 1 ? 1d / 3 : double.PositiveInfinity;
                Assert.AreEqual(expectedDensity, distribution.PDF(endpoint));
            }
        }

        /// <summary>Moment divergence is determined by moment order, independently of MLE-information regularity.</summary>
        [TestMethod]
        public void GeneralizedLogistic_Moments_RespectEachExistenceBoundary()
        {
            foreach (double sign in new[] { -1d, 1d })
            {
                Assert.IsTrue(double.IsNaN(new GeneralizedLogistic(0, 1, sign).Mean));
                Assert.IsTrue(double.IsNaN(new GeneralizedLogistic(0, 1, sign * 0.5).StandardDeviation));
                Assert.IsTrue(double.IsNaN(new GeneralizedLogistic(0, 1, sign / 3).Skewness));
                Assert.IsTrue(double.IsNaN(new GeneralizedLogistic(0, 1, sign * 0.25).Kurtosis));
                Assert.IsFalse(double.IsNaN(new GeneralizedLogistic(0, 1, sign * 0.49).Mean));
            }
        }

        /// <summary>Symmetric L-moments produce a valid normal fit in the public xi/alpha/kappa coordinates.</summary>
        [TestMethod]
        public void GeneralizedNormal_LinearMoments_PreserveNormalLimit()
        {
            var distribution = new GeneralizedNormal();
            double[] parameters = distribution.ParametersFromLinearMoments(new[] { 2d, 3d, 0d, 0.122601719540891 });
            Assert.AreEqual(2d, parameters[0]);
            Assert.AreEqual(3 * Math.Sqrt(Math.PI), parameters[1], 1E-14);
            Assert.AreEqual(0d, parameters[2]);
            double[] moments = distribution.LinearMomentsFromParameters(new[] { 2d, 3d, 0d });
            Assert.AreEqual(2d, moments[0]);
            Assert.AreEqual(3 / Math.Sqrt(Math.PI), moments[1], 1E-14);
            Assert.AreEqual(0d, moments[2]);
            Assert.AreEqual(0.12260171954089095, moments[3], 2E-15); // 30*asin(1/3)/pi-9.
            foreach (double k in new[] { -1E-6, 1E-6 })
            {
                moments = distribution.LinearMomentsFromParameters(new[] { 2d, 3d, k });
                parameters = distribution.ParametersFromLinearMoments(moments);
                Assert.AreEqual(2d, parameters[0], 1E-12);
                Assert.AreEqual(3d, parameters[1], 3E-12);
                Assert.AreEqual(k, parameters[2], 2E-13);
            }
        }

        /// <summary>Initialization rejects invalid samples before producing bounds or mutating fitted parameters.</summary>
        [TestMethod]
        public void GeneralizedInitialization_RejectsInvalidSamplesAndKeepsCenteredBounds()
        {
            foreach (var distribution in new IMaximumLikelihoodEstimation[] { new GeneralizedNormal(), new GeneralizedLogistic() })
            {
                foreach (var sample in new[] { new[] { 1d, 2d, 3d }, new[] { 1d, 1d, 1d, 1d }, new[] { 0d, 1d, 2d, double.NaN }, new[] { 0d, 1d, 2d, double.PositiveInfinity } })
                    Assert.Throws<ArgumentException>(() => distribution.GetParameterConstraints(sample));
                var constraints = distribution.GetParameterConstraints(new[] { -2d, -1d, 1d, 2d });
                for (int i = 0; i < 3; i++)
                {
                    Assert.IsTrue(!double.IsInfinity(constraints.Item2[i]) && !double.IsInfinity(constraints.Item3[i]));
                    Assert.IsLessThan(constraints.Item3[i], constraints.Item2[i]);
                    Assert.IsTrue(constraints.Item1[i] >= constraints.Item2[i] && constraints.Item1[i] <= constraints.Item3[i]);
                }
                var candidate = (UnivariateDistributionBase)distribution;
                candidate.SetParameters(constraints.Item1);
                Assert.IsFalse(double.IsInfinity(candidate.LogLikelihood(new[] { -2d, -1d, 1d, 2d })));
            }
        }

        /// <summary>Closed-form normal information agrees with independent complete-tail integration.</summary>
        [TestMethod]
        public void GeneralizedNormal_MleCovariance_MatchesIndependentRInformation() => CheckOracleCovariance("GNO");

        /// <summary>Fixing hondo at minus one requires inversion of the GLO information block.</summary>
        [TestMethod]
        public void GeneralizedLogistic_MleCovariance_MatchesIndependentRInformation() => CheckOracleCovariance("GLO");

        /// <summary>Kappa covariance uses all four score coordinates and includes heavy-tailed regular cases.</summary>
        [TestMethod]
        public void KappaFour_MleCovariance_MatchesIndependentRInformation() => CheckOracleCovariance("K4");

        /// <summary>The zero-shape GNO covariance still estimates shape; it is not the two-parameter normal matrix.</summary>
        [TestMethod]
        public void GeneralizedNormal_MleCovariance_ExactZeroShapeAndLargeShape()
        {
            double[,] covariance = new GeneralizedNormal(0, 1, 0).ParameterCovariance(1, ParameterEstimationMethod.MaximumLikelihood);
            double[,] expected = { { 7d / 6, 0, 1d / 3 }, { 0, 0.5, 0 }, { 1d / 3, 0, 2d / 3 } };
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) Assert.AreEqual(expected[i, j], covariance[i, j], 2E-15);
            covariance = new GeneralizedNormal(0, 1, 30).ParameterCovariance(1, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(1d, covariance[0, 0]);
            Assert.AreEqual(-30d, covariance[0, 1]);
            Assert.AreEqual(900.5, covariance[1, 1]);
            Assert.AreEqual(450d, covariance[2, 2]);
        }

        /// <summary>Scale factors are combined before intermediate exponential or squared-shape overflow.</summary>
        [TestMethod]
        public void GeneralizedExtremeScales_PreserveRepresentableQuantilesAndCovarianceEntries()
        {
            double normalLatent = -37.047096299361199; // R qnorm(1e-300), independently rounded binary64 reference.
            double expectedNormal = -Math.Exp(Math.Log(1E-200) - 30 * normalLatent - Math.Log(30));
            Assert.AreEqual(expectedNormal, new GeneralizedNormal(0, 1E-200, 30).InverseCDF(1E-300), Math.Abs(expectedNormal) * 3E-12);
            double expectedLogistic = -Math.Exp(Math.Log(1E-308) - 2 * Math.Log(1E-300) - Math.Log(2));
            Assert.AreEqual(expectedLogistic, new GeneralizedLogistic(0, 1E-308, 2).InverseCDF(1E-300), Math.Abs(expectedLogistic) * 3E-12);
            double[,] covariance = new GeneralizedNormal(0, 1E200, 30).ParameterCovariance(1, ParameterEstimationMethod.MaximumLikelihood);
            double expectedCross = Math.Exp(Math.Log(1E200) - 900 + Math.Log(900d / 901));
            Assert.AreEqual(expectedCross, covariance[0, 2], expectedCross * 3E-12);
            covariance = new GeneralizedNormal(0, 1E-200, 1E200).ParameterCovariance(1, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(1d, covariance[1, 1], 3E-13);
            covariance = new GeneralizedNormal(0, 1E-154, 1.4E154).ParameterCovariance(1, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(9.8E307, covariance[2, 2], 3E294);
            foreach (var distribution in new UnivariateDistributionBase[] { new GeneralizedNormal(-1E308, 1E308, .5), new GeneralizedLogistic(-1E308, 1E308, .5) })
            {
                Assert.AreEqual(1E308, distribution.Maximum, 1E293);
                Assert.AreEqual(1d, distribution.CDF(1.5E308));
                Assert.AreEqual(double.NegativeInfinity, distribution.LogPDF(1.5E308));
            }
        }

        /// <summary>Covariances transform in xi/alpha coordinates and quantile variance is their delta contraction.</summary>
        [TestMethod]
        public void GeneralizedUncertainty_AffineScaleSampleSizeAndDeltaMethod()
        {
            foreach (var pair in new[]
            {
                (new GeneralizedNormal(0, 1, -0.2) as IStandardError, new GeneralizedNormal(7, 3, -0.2) as IStandardError),
                (new GeneralizedLogistic(0, 1, 0.2) as IStandardError, new GeneralizedLogistic(7, 3, 0.2) as IStandardError),
                (new KappaFour(0, 1, -1, -0.2) as IStandardError, new KappaFour(7, 3, -1, -0.2) as IStandardError)
            })
            {
                double[,] unit = pair.Item1.ParameterCovariance(1, ParameterEstimationMethod.MaximumLikelihood);
                double[,] scaled = pair.Item2.ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood);
                double[] gradient = pair.Item2.QuantileGradient(0.9);
                double expected = 0;
                for (int i = 0; i < gradient.Length; i++)
                for (int j = 0; j < gradient.Length; j++)
                {
                    Assert.AreEqual(unit[i, j] * (i < 2 ? 3 : 1) * (j < 2 ? 3 : 1) / 100, scaled[i, j], 2E-11);
                    expected += gradient[i] * scaled[i, j] * gradient[j];
                }
                Assert.IsGreaterThan(0d, expected);
                Assert.AreEqual(expected, pair.Item2.QuantileVariance(0.9, 100, ParameterEstimationMethod.MaximumLikelihood), 2E-11 * expected);
            }
        }

        /// <summary>Scalar MLE variance preserves affine scale and reports true overflow independently of matrix range.</summary>
        [TestMethod]
        [DataRow("GNO", 0.011666666666666667)]
        [DataRow("GLO", 0.03113972242232403)]
        [DataRow("K4", 0.01550447368065933)]
        public void GeneralizedQuantileVariance_ScaleAndOverflowClassification(string family, double unitVariance)
        {
            // GNO is exact; GLO/K4 contract the frozen R zero-shape matrices with analytical median gradients.
            IStandardError Distribution(double alpha) => family == "GNO" ? new GeneralizedNormal(0, alpha, 0)
                : family == "GLO" ? new GeneralizedLogistic(0, alpha, 0) : new KappaFour(0, alpha, 0, 0);
            foreach (double alpha in new[] { 1E-150, 1d, 1E150 })
            {
                double expected = unitVariance * alpha * alpha;
                Assert.AreEqual(expected, Distribution(alpha).QuantileVariance(.5, 100, ParameterEstimationMethod.MaximumLikelihood), expected * 5E-8, family);
            }
            foreach (double probability in new[] { .1, .5, .9 })
                Assert.IsTrue(double.IsPositiveInfinity(Distribution(1E200).QuantileVariance(probability, 100, ParameterEstimationMethod.MaximumLikelihood)), family);
        }

        /// <summary>A finite scalar can survive underflow or overflow in individual physical covariance entries.</summary>
        [TestMethod]
        [DataRow("GNO", 1.002832088139248E308)]
        [DataRow("GLO", 3.032823215830043E-216)]
        [DataRow("K4", 2.3597144586448996E-307)]
        public void GeneralizedQuantileVariance_RecoversFiniteScalarAfterCovarianceRangeLoss(string family, double expected)
        {
            // Frozen R covariance entries contracted with analytical gradients using 400-digit Decimal arithmetic.
            // The GNO case instead uses its exact closed-form covariance at k=2.
            IStandardError distribution = family == "GNO" ? new GeneralizedNormal(0, 1E155, 2)
                : family == "GLO" ? new GeneralizedLogistic(0, 1E-170, .2) : new KappaFour(0, 1E-170, -1, -.2);
            double probability = family == "GNO" ? .5 : family == "GLO" ? 1E-300 : 1 - 1E-16;
            Assert.AreEqual(expected, distribution.QuantileVariance(probability, 100, ParameterEstimationMethod.MaximumLikelihood), expected * 5E-8, family);
        }

        /// <summary>Exact positive quadratic terms retain variance lost by overflowing gradients or endpoint cancellation.</summary>
        [TestMethod]
        [DataRow(20d, 1E-200, 1E-300, 2.5768579804728683E244)]
        [DataRow(30d, 1E200, 0.9999999999999999, 4.158459661237856E185)]
        public void GeneralizedNormal_QuantileVariance_PreservesExponentialScaleAndCancellation(double kappa, double alpha, double probability, double expected)
        {
            // 450-digit Decimal evaluation of the exact positive quadratic, with R qnorm references
            // z=-37.047096299361199 and z=8.209536151601387 for these binary64 probabilities.
            var distribution = new GeneralizedNormal(0, alpha, kappa);
            Assert.AreEqual(expected, distribution.QuantileVariance(probability, 100, ParameterEstimationMethod.MaximumLikelihood), expected * 2E-10);
        }

        /// <summary>Uncertainty validates its own regularity domain without narrowing distribution validity.</summary>
        [TestMethod]
        public void GeneralizedUncertainty_RejectsInvalidArgumentsAndNonregularShapes()
        {
            foreach (var distribution in new IStandardError[] { new GeneralizedNormal(), new GeneralizedLogistic(), new KappaFour() })
            {
                foreach (double p in new[] { 0d, 1d, -0.1, 1.1, double.NaN, double.PositiveInfinity })
                {
                    Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.QuantileGradient(p));
                    Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.QuantileVariance(p, 10, ParameterEstimationMethod.MaximumLikelihood));
                }
                foreach (int n in new[] { 0, -1 })
                    Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.ParameterCovariance(n, ParameterEstimationMethod.MaximumLikelihood));
            }
            foreach (var distribution in new UnivariateDistributionBase[] { new GeneralizedLogistic(0, 1, -0.5), new GeneralizedLogistic(0, 1, 0.5), new KappaFour(0, 1, 0.5, 0), new KappaFour(0, 1, 0, 0.5), new KappaFour(0, 1, -1, -0.5) })
            {
                Assert.IsTrue(distribution.ParametersValid);
                Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => ((IStandardError)distribution).ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood));
            }
        }

        /// <summary>MLE covariance is never substituted for a different estimator's uncertainty.</summary>
        [TestMethod]
        public void GeneralizedUncertainty_OtherEstimatorsRemainExplicitlyUnsupported()
        {
            foreach (var distribution in new IStandardError[] { new GeneralizedNormal(), new GeneralizedLogistic(), new KappaFour() })
            foreach (var method in new[] { ParameterEstimationMethod.MethodOfMoments, ParameterEstimationMethod.MethodOfLinearMoments })
            {
                Assert.ThrowsExactly<NotImplementedException>(() => distribution.ParameterCovariance(100, method));
                Assert.ThrowsExactly<NotImplementedException>(() => distribution.QuantileVariance(0.9, 100, method));
            }
        }

        /// <summary>All score components integrate to zero and the principal information block matches the frozen R matrix.</summary>
        [TestMethod]
        public void KappaExpectedInformation_MeansErrorsAndMatrix_MatchIndependentR()
        {
            foreach (var group in DistributionOracle.Read("generalized-fisher.csv")
                .Where(r => r["family"] != "GNO" && r["quantity"] == "information")
                .GroupBy(r => r["family"] + ":" + r["kappa"] + ":" + r["hondo"]))
            {
                var first = group.First();
                int count = first["family"] == "K4" ? 4 : 3;
                double[,] information = KappaExpectedInformation.ExpectedInformation(Number(first, "kappa"), Number(first, "hondo"), count,
                    out double[] means, out double[] meanErrors, out double[,] informationErrors);
                for (int i = 0; i < count; i++)
                    Assert.IsLessThanOrEqualTo(4 * meanErrors[i] + 8E-15, Math.Abs(means[i]), Context(first));
                foreach (var row in group)
                {
                    int i = int.Parse(row["row"]) - 1, j = int.Parse(row["column"]) - 1;
                    double expected = Number(row, "value");
                    Assert.AreEqual(expected, information[i, j], Math.Max(5E-9 * Math.Max(1, Math.Abs(expected)),
                        8 * (Number(row, "absolute_error_estimate") + informationErrors[i, j])), Context(row));
                    Assert.IsLessThanOrEqualTo(1E-12 + 1E-10 * Math.Abs(information[i, j]), informationErrors[i, j]);
                }
            }
        }

        /// <summary>A mathematically regular but unresolved boundary approach must not produce a finite-looking covariance.</summary>
        [TestMethod]
        public void KappaExpectedInformation_UnresolvedBoundaryApproach_ThrowsNumericalFailure()
        {
            var distribution = new KappaFour(0, 1, 0.499999999999, 0);
            Assert.IsTrue(distribution.ParametersValid);
            Assert.ThrowsExactly<InvalidOperationException>(() => distribution.ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood));
            // Full Kappa information has a very narrow upper-tail feature at these shapes.
            // Its unresolvable finite-shape calculation must fail explicitly, rather than
            // mistaking log(t)=log(survival) for a uniform approximation in hondo.
            var extremeHondo = new KappaFour(0, 1, 0, -1E200);
            Assert.IsTrue(extremeHondo.ParametersValid);
            Assert.ThrowsExactly<InvalidOperationException>(() => extremeHondo.ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood));
        }

        /// <summary>The first four interior floating-point arguments are checked separately at every finite Kappa endpoint.</summary>
        [TestMethod]
        public void KappaFour_FirstFourAdjacentEndpointValues_MatchIndependentDecimal()
        {
            int count = 0;
            foreach (var row in DistributionOracle.Read("generalized-adjacent-boundaries.csv"))
            {
                var distribution = new KappaFour(0, 1, Number(row, "kappa"), Number(row, "hondo"));
                double x = Number(row, "x");
                string context = $"k={row["kappa"]}, h={row["hondo"]}, {row["endpoint"]} step={row["step"]}, x={x:R}";
                Assert.AreEqual(Number(row, "logcdf"), distribution.LogCDF(x), 5E-10, context);
                Assert.AreEqual(Number(row, "logpdf"), distribution.LogPDF(x), 5E-10, context);
                count++;
            }
            Assert.AreEqual(180, count);
            var changing = new KappaFour(0, 1, -2, 10);
            Assert.AreEqual(49.5, changing.Minimum);
            changing.Xi = 2;
            Assert.AreEqual(51.5, changing.Minimum);
            changing.Alpha = 3;
            Assert.AreEqual(150.5, changing.Minimum);
            changing.Kappa = -1;
            Assert.AreEqual(29d, changing.Minimum);
            changing.Hondo = 2;
            Assert.AreEqual(5d, changing.Minimum);
        }

        /// <summary>Analytical shape derivatives retain their nonzero normal/logistic limit.</summary>
        [TestMethod]
        public void GeneralizedQuantileGradients_ZeroShapeAndSingularJacobians()
        {
            foreach (var distribution in new IStandardError[] { new GeneralizedNormal(2, 3, 0), new GeneralizedLogistic(2, 3, 0) })
            {
                double z = distribution is GeneralizedNormal ? 1 : Math.Log(9d);
                double p = distribution is GeneralizedNormal ? 0.8413447460685429 : 0.9;
                double[] gradient = distribution.QuantileGradient(p);
                Assert.AreEqual(1d, gradient[0]);
                Assert.AreEqual(z, gradient[1], 3E-14);
                Assert.AreEqual(-3 * z * z / 2, gradient[2], 3E-14);
                distribution.QuantileJacobian(new[] { 0.2, 0.2, 0.9 }, out double determinant);
                Assert.AreEqual(0d, determinant);
            }
            new KappaFour(0, 1, 0, 0).QuantileJacobian(new[] { 0.2, 0.2, 0.6, 0.9 }, out double kappaDeterminant);
            Assert.AreEqual(0d, kappaDeterminant);
        }

        private static void CheckOracleCovariance(string family)
        {
            int count = 0;
            foreach (var group in DistributionOracle.Read("generalized-fisher.csv").Where(r => r["family"] == family && r["quantity"] == "covariance").GroupBy(r => r["kappa"] + ":" + r["hondo"]))
            {
                var distribution = (IStandardError)Create(group.First());
                double[,] covariance = distribution.ParameterCovariance(1, ParameterEstimationMethod.MaximumLikelihood);
                foreach (var row in group)
                {
                    int i = int.Parse(row["row"]) - 1;
                    int j = int.Parse(row["column"]) - 1;
                    double expected = Number(row, "value");
                    double tolerance = Math.Max(5E-9 * Math.Max(1, Math.Abs(expected)), 8 * Number(row, "absolute_error_estimate"));
                    Assert.AreEqual(expected, covariance[i, j], tolerance, Context(row));
                    Assert.AreEqual(covariance[i, j], covariance[j, i], 1E-12);
                    if (i == j) Assert.IsGreaterThan(0d, covariance[i, i]);
                }
                count++;
            }
            Assert.AreEqual(family == "GNO" ? 9 : family == "GLO" ? 11 : 18, count);
        }

        private static UnivariateDistributionBase Create(Dictionary<string, string> row)
        {
            double k = Number(row, "kappa");
            if (row["family"] == "GNO") return new GeneralizedNormal(0, 1, k);
            if (row["family"] == "GLO") return new GeneralizedLogistic(0, 1, k);
            return new KappaFour(0, 1, k, Number(row, "hondo"));
        }

        private static double Number(Dictionary<string, string> row, string column) => DistributionOracle.Number(row[column]);
        private static string Context(Dictionary<string, string> row) => $"{row["family"]} k={row["kappa"]} h={row["hondo"]} {row["quantity"]} [{row["row"]},{row["column"]}]";
    }
}
