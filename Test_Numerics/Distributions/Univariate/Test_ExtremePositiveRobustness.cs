using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Independent probability, moment, and uncertainty regressions for six extreme and positive families.</summary>
    [TestClass]
    public class Test_ExtremePositiveRobustness
    {
        /// <summary>Frozen R and defining-formula oracles catch cancelled tails, artificial shape plateaus and endpoint errors.</summary>
        [TestMethod]
        [DataRow("Exponential")]
        [DataRow("Gamma")]
        [DataRow("GEV")]
        [DataRow("GPA")]
        [DataRow("Gumbel")]
        [DataRow("Weibull")]
        public void ProbabilityFunctionsMatchIndependentOracle(string family)
        {
            foreach (var row in Rows(family).Where(r => r["quantity"].StartsWith("Log", StringComparison.Ordinal) || r["quantity"] == "InverseCDF"))
            {
                var distribution = Create(row);
                double x = Number(row, "x"), p = Number(row, "p");
                string quantity = row["quantity"];
                double actual = quantity == "LogPDF" ? distribution.LogPDF(x)
                    : quantity == "LogCDF" ? distribution.LogCDF(x)
                    : quantity == "LogCCDF" ? distribution.LogCCDF(x) : distribution.InverseCDF(p);
                OracleAssert(row, actual);
                if (quantity == "LogPDF") ScalarAssert(Math.Exp(Number(row, "value")), distribution.PDF(x), row, 2E-10);
                if (quantity == "LogCDF") ScalarAssert(Math.Exp(Number(row, "value")), distribution.CDF(x), row, 2E-10);
                if (quantity == "LogCCDF") ScalarAssert(Math.Exp(Number(row, "value")), distribution.CCDF(x), row, 2E-10);
            }
        }

        /// <summary>Independent complete-tail moments distinguish bounded positive shapes from nonexistent heavy-tail moments.</summary>
        [TestMethod]
        [DataRow("Exponential")]
        [DataRow("Gamma")]
        [DataRow("GEV")]
        [DataRow("GPA")]
        [DataRow("Gumbel")]
        [DataRow("Weibull")]
        public void MomentsAndSupportMatchIndependentOracle(string family)
        {
            foreach (var row in Rows(family).Where(r => new[] { "Mean", "StandardDeviation", "Skewness", "Kurtosis", "Minimum", "Maximum" }.Contains(r["quantity"])))
            {
                var distribution = Create(row);
                double actual = row["quantity"] == "Mean" ? distribution.Mean
                    : row["quantity"] == "StandardDeviation" ? distribution.StandardDeviation
                    : row["quantity"] == "Skewness" ? distribution.Skewness
                    : row["quantity"] == "Kurtosis" ? distribution.Kurtosis
                    : row["quantity"] == "Minimum" ? distribution.Minimum : distribution.Maximum;
                // Gumbel's existing rounded published skewness remains its documented precision.
                OracleAssert(row, actual, family == "Gumbel" && row["quantity"] == "Skewness" ? 5E-5 : 0);
            }
        }

        /// <summary>Actual-quantile derivatives use independent R qgamma derivatives and analytic other-family gradients.</summary>
        [TestMethod]
        [DataRow("Exponential")]
        [DataRow("Gamma")]
        [DataRow("GEV")]
        [DataRow("GPA")]
        [DataRow("Gumbel")]
        [DataRow("Weibull")]
        public void QuantileGradientsMatchIndependentOracle(string family)
        {
            foreach (var row in Rows(family).Where(r => r["quantity"].StartsWith("Gradient", StringComparison.Ordinal)))
            {
                if (row["status"] == "unresolved-quantile-underflow") continue;
                int index = row["quantity"] == "GradientScale" ? 0
                    : row["quantity"] == "GradientShape" ? 1 : int.Parse(row["quantity"].Substring(8)) - 1;
                OracleAssert(row, ((IStandardError)Create(row)).QuantileGradient(Number(row, "p"))[index]);
            }
        }

        /// <summary>Covariance matches independent order-statistic, Fisher and moment-Jacobian references in public coordinates.</summary>
        [TestMethod]
        [DataRow("Exponential")]
        [DataRow("Gamma")]
        [DataRow("GEV")]
        [DataRow("GPA")]
        [DataRow("Gumbel")]
        [DataRow("Weibull")]
        public void CovarianceMatchesIndependentOracle(string family)
        {
            var matrices = new Dictionary<string, double[,]>();
            foreach (var row in Rows(family).Where(r => r["quantity"].StartsWith("Covariance", StringComparison.Ordinal)))
            {
                var method = row["case"].StartsWith("MoM", StringComparison.Ordinal)
                    ? ParameterEstimationMethod.MethodOfMoments : ParameterEstimationMethod.MaximumLikelihood;
                var distribution = (IStandardError)Create(row);
                int n = (int)Number(row, "n");
                if (row["quantity"] == "CovarianceDefined")
                {
                    Assert.Throws<ArgumentOutOfRangeException>(() => distribution.ParameterCovariance(n, method), Description(row));
                    continue;
                }
                string key = row["case"] + "/" + row["scale"] + "/" + row["shape"] + "/" + row["n"];
                if (!matrices.TryGetValue(key, out var covariance))
                    matrices.Add(key, covariance = distribution.ParameterCovariance(n, method));
                int i = row["quantity"][10] - '1', j = row["quantity"][11] - '1';
                // Preserve the existing published rounded Gumbel and Weibull Fisher constants.
                double publishedPrecision = family == "Gumbel" ? 9E-5 : family == "Weibull" ? 2E-6 : 0;
                OracleAssert(row, covariance[i, j], publishedPrecision);
                Assert.AreEqual(covariance[i, j], covariance[j, i], Description(row));
            }
        }

        /// <summary>Quantile uncertainty rejects endpoints, nonfinite probabilities and invalid sample sizes for every family.</summary>
        [TestMethod]
        [DataRow("Exponential")]
        [DataRow("Gamma")]
        [DataRow("GEV")]
        [DataRow("GPA")]
        [DataRow("Gumbel")]
        [DataRow("Weibull")]
        public void UncertaintyRejectsInvalidInputs(string family)
        {
            var distribution = (IStandardError)Create(family, 0, 2, family == "Gamma" || family == "Weibull" ? 2 : 0);
            foreach (double p in new[] { -1d, 0d, 1d, 2d, double.NaN, double.PositiveInfinity, double.NegativeInfinity })
            {
                Assert.Throws<ArgumentOutOfRangeException>(() => distribution.QuantileGradient(p), family + "/gradient");
                Assert.Throws<ArgumentOutOfRangeException>(() => distribution.QuantileVariance(p, 100, ParameterEstimationMethod.MaximumLikelihood), family + "/variance");
            }
            foreach (int n in new[] { 0, -1 })
                Assert.Throws<ArgumentOutOfRangeException>(() => distribution.ParameterCovariance(n, ParameterEstimationMethod.MaximumLikelihood), family);
            Assert.Throws<NotImplementedException>(() => distribution.ParameterCovariance(100, ParameterEstimationMethod.MethodOfLinearMoments), family);
            var probabilities = Enumerable.Repeat(.5, ((UnivariateDistributionBase)distribution).NumberOfParameters).ToArray();
            distribution.QuantileJacobian(probabilities, out double determinant);
            Assert.AreEqual(0d, determinant, family + "/duplicate quantile determinant");
        }

        /// <summary>Validation inspects the candidate vector, including wrong lengths and invalid shape entries.</summary>
        [TestMethod]
        [DataRow("Exponential")]
        [DataRow("Gamma")]
        [DataRow("GEV")]
        [DataRow("GPA")]
        [DataRow("Gumbel")]
        [DataRow("Weibull")]
        public void CandidateValidationRejectsMalformedVectors(string family)
        {
            var distribution = Create(family, 0, 2, family == "Gamma" || family == "Weibull" ? 2 : 0);
            foreach (var parameters in new[] { Array.Empty<double>(), new[] { 1d }, new[] { 1d, 2d, 3d, 4d } })
            {
                Assert.IsNotNull(distribution.ValidateParameters(parameters, false), family);
                Assert.Throws<ArgumentOutOfRangeException>(() => distribution.ValidateParameters(parameters, true), family);
            }
            var candidate = distribution.GetParameters;
            candidate[candidate.Length - 1] = double.NaN;
            Assert.IsNotNull(distribution.ValidateParameters(candidate, false), family);
            Assert.Throws<ArgumentOutOfRangeException>(() => distribution.InverseCDF(double.NaN), family);
        }

        /// <summary>Initialization keeps tiny scales feasible, handles signed centers, and rejects unusable samples.</summary>
        [TestMethod]
        [DataRow("Exponential")]
        [DataRow("Gamma")]
        [DataRow("GEV")]
        [DataRow("GPA")]
        [DataRow("Gumbel")]
        [DataRow("Weibull")]
        public void InitializationHasFiniteOrderedFeasibleBounds(string family)
        {
            var distribution = (IMaximumLikelihoodEstimation)Create(family, 0, 2, family == "Gamma" || family == "Weibull" ? 2 : 0);
            foreach (var sample in new[] { new[] { 1d, 1d, 1d, 1d }, new[] { 1d, 2d, 3d, double.NaN }, new[] { 1d, 2d, 3d, double.PositiveInfinity }, new[] { 1d } })
                Assert.Throws<ArgumentOutOfRangeException>(() => distribution.GetParameterConstraints(sample), family);
            if (family == "Gamma" || family == "Weibull")
                Assert.Throws<ArgumentOutOfRangeException>(() => distribution.GetParameterConstraints(new[] { 0d, 1d, 2d, 3d }), family);
            foreach (double scale in new[] { 1E-200, 1d, 1E200 })
            {
                var sample = new[] { 1d, 2d, 4d, 7d, 8d, 10d }.Select(x => x * scale).ToArray();
                CheckConstraints(distribution.GetParameterConstraints(sample), family);
                if (family != "Gamma" && family != "Weibull")
                    CheckConstraints(distribution.GetParameterConstraints(sample.Select(x => x - 5 * scale).ToArray()), family);
            }
        }

        /// <summary>The Weibull median gradient and its variance correspond to the actual inverse CDF.</summary>
        [TestMethod]
        public void WeibullMedianGradientAndVarianceAreCorrect()
        {
            var distribution = new Weibull(2, 2);
            var gradient = distribution.QuantileGradient(.5);
            Assert.AreEqual(.832554611157698, gradient[0], 1E-14);
            Assert.AreEqual(.152571011039570, gradient[1], 1E-14);
            Assert.AreEqual(.00955664647618049, distribution.QuantileVariance(.5, 100, ParameterEstimationMethod.MaximumLikelihood), 1E-15);
        }

        /// <summary>Shape one is a reverse exponential for GEV and a uniform distribution for GPA.</summary>
        [TestMethod]
        public void ShapeOneHasCorrectMedianModeAndEndpointDensity()
        {
            var gev = new GeneralizedExtremeValue(0, 1, 1);
            Assert.AreEqual(.3068528194400547, gev.Median, 1E-15);
            Assert.AreEqual(1d, gev.Mode);
            Assert.AreEqual(1d, gev.PDF(1));
            var gpa = new GeneralizedPareto(0, 1, 1);
            Assert.AreEqual(.5, gpa.Median);
            Assert.IsTrue(double.IsNaN(gpa.Mode), "The uniform mode is nonunique.");
            Assert.AreEqual(1d, gpa.PDF(1));
            Assert.AreEqual(1d, new GeneralizedPareto(0, 2, 2).Mode);
        }

        /// <summary>Finite affine standardized coordinates remain usable when the unscaled subtraction overflows.</summary>
        [TestMethod]
        public void AffineStandardizationAvoidsSpuriousOverflow()
        {
            foreach (string family in new[] { "Exponential", "Gumbel", "GEV", "GPA" })
            {
                var unit = Create(family, 0, 1, 0);
                var shifted = Create(family, -1E308, 1E308, 0);
                Assert.AreEqual(unit.LogCDF(2), shifted.LogCDF(1E308), 1E-13, family);
                Assert.AreEqual(unit.LogCCDF(2), shifted.LogCCDF(1E308), 1E-13, family);
                Assert.AreEqual(unit.LogPDF(2) - Math.Log(1E308), shifted.LogPDF(1E308), 1E-12, family);
            }
        }

        /// <summary>A capped named Gamma approximation is constant beyond the cap and reflects negative skew consistently.</summary>
        [TestMethod]
        public void GammaNamedApproximationClippingAndReflectionAreConsistent()
        {
            foreach (double p in new[] { .01, .25, .5, .75, .99 })
            {
                Assert.AreEqual(GammaDistribution.FrequencyFactorKp(9.75, p), GammaDistribution.FrequencyFactorKp(100, p), 1E-14);
                Assert.AreEqual(-GammaDistribution.FrequencyFactorKp(9.75, 1 - p), GammaDistribution.FrequencyFactorKp(-100, p), 1E-13);
                Assert.AreEqual(-GammaDistribution.FrequencyFactorKp(3, 1 - p), GammaDistribution.FrequencyFactorKp(-3, p), 1E-13);
            }
        }

        /// <summary>Independent unit-scale quadratics retain finite scaling and the correct overflow classification.</summary>
        [TestMethod]
        [DataRow("Exponential", .0048564848377901943)]
        [DataRow("Gamma", .015723154225565583)]
        [DataRow("GEV", .015433980376270423)]
        [DataRow("GPA", .0068559015049324415)]
        [DataRow("Gumbel", .013787478943464875)]
        [DataRow("Weibull", .0023891616190451241)]
        public void QuantileVariancePreservesScaleAndOverflowClassification(string family, double unitVariance)
        {
            double shape = family == "Gamma" || family == "Weibull" ? 2 : 0;
            foreach (double scale in new[] { 1E-150, 1d, 1E150 })
            {
                var distribution = (IStandardError)Create(family, 0, scale, shape);
                double expected = unitVariance * scale * scale;
                Assert.AreEqual(expected, distribution.QuantileVariance(.5, 100, ParameterEstimationMethod.MaximumLikelihood), expected * 5E-8, family);
            }
            var large = (IStandardError)Create(family, 0, 1E200, shape);
            foreach (double probability in new[] { .1, .5, .9 })
                Assert.IsTrue(double.IsPositiveInfinity(large.QuantileVariance(probability, 100, ParameterEstimationMethod.MaximumLikelihood)), family + "/overflowed positive variance");
        }

        /// <summary>A finite affine inverse survives an overflowing displacement followed by an opposite location.</summary>
        [TestMethod]
        [DataRow("Exponential", .8646647167633873)]
        [DataRow("Gumbel", .8734230184931167)]
        [DataRow("GEV", .8734230184931167)]
        [DataRow("GPA", .8646647167633873)]
        public void AffineInverseAvoidsSpuriousOverflow(string family, double probability)
        {
            var distribution = Create(family, -1E308, 1E308, 0);
            Assert.AreEqual(1E308, distribution.InverseCDF(probability), 2E293, family);
        }

        /// <summary>Tiny representable probabilities and adjacent support doubles do not acquire artificial zeros.</summary>
        [TestMethod]
        public void TinyProbabilitiesAndAdjacentSupportPointsRemainDistinct()
        {
            foreach (var distribution in new UnivariateDistributionBase[] { new Exponential(0, 1), new GammaDistribution(1, 1), new Weibull(1, 1), new GeneralizedPareto(0, 1, 0) })
            {
                Assert.AreEqual(double.Epsilon, distribution.InverseCDF(double.Epsilon), distribution.DisplayName);
                Assert.AreEqual(Math.Log(double.Epsilon), distribution.LogCDF(double.Epsilon), 1E-12, distribution.DisplayName);
            }
            double below = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(2d) - 1);
            double above = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(2d) + 1);
            foreach (var distribution in new UnivariateDistributionBase[] { new GeneralizedExtremeValue(0, 1, .5), new GeneralizedPareto(0, 1, .5) })
            {
                Assert.IsGreaterThan(0d, distribution.CCDF(below), distribution.DisplayName);
                Assert.IsFalse(double.IsInfinity(distribution.LogCCDF(below)), distribution.DisplayName);
                Assert.AreEqual(double.NegativeInfinity, distribution.LogCCDF(2));
                Assert.AreEqual(0d, distribution.PDF(above));
            }
            Assert.AreEqual(-1000d, new Gumbel(0, 1).LogCCDF(1000));
            Assert.AreEqual(-1000d, new GeneralizedExtremeValue(0, 1, 0).LogCCDF(1000));
        }

        /// <summary>The zero-shape derivative combines its dimensionless factor before the large scale.</summary>
        [TestMethod]
        public void RepresentableExtremeShapeGradientAvoidsIntermediateOverflow()
        {
            double gevProbability = Math.Exp(-Math.Exp(-1.4));
            double gpaProbability = 1 - Math.Exp(-1.4);
            Assert.AreEqual(-9.8E307, new GeneralizedExtremeValue(0, 1E308, 0).QuantileGradient(gevProbability)[2], 1E294);
            Assert.AreEqual(-9.8E307, new GeneralizedPareto(0, 1E308, 0).QuantileGradient(gpaProbability)[2], 1E294);
        }

        /// <summary>Physical GPA tail values remain finite when their unit-scale counterparts exceed the floating-point range.</summary>
        [TestMethod]
        public void GpaScaleRescuesOverflowedUnitQuantileAndVariance()
        {
            var distribution = new GeneralizedPareto(0, 1E-200, -20);
            double probability = 1 - 1E-16;
            Assert.AreEqual(6.1768265779813681E117, distribution.InverseCDF(probability), 2E105);
            Assert.AreEqual(-2.2660800481989075E119, distribution.QuantileGradient(probability)[2], 1E107);
            Assert.AreEqual(2.2588688104431192E239,
                distribution.QuantileVariance(probability, 100, ParameterEstimationMethod.MaximumLikelihood), 2E227);
        }

        /// <summary>Gradient normalization retains a Gamma variance rescued from a tiny unit quantile by its physical scale.</summary>
        [TestMethod]
        public void GammaScaleRescuesUnderflowedUnitQuantileVariance()
        {
            var distribution = new GammaDistribution(1E200, .001);
            // R contraction of the frozen actual-qgamma gradient and trigamma covariance, after physical scaling.
            Assert.AreEqual(1.3215876806248227E-199,
                distribution.QuantileVariance(.5, 100, ParameterEstimationMethod.MaximumLikelihood), 2E-207);
        }

        /// <summary>Negative-skew approximate frequency factors preserve a tiny lower-tail probability during reflection.</summary>
        [TestMethod]
        public void GammaNamedApproximationReflectsTinyProbabilityWithoutRoundingItAway()
        {
            // Standalone R qnorm and the published magnitude-three Kirby polynomial.
            Assert.AreEqual(-95.085021656541102, GammaDistribution.FrequencyFactorKp(-3, 1E-20), 2E-10);
        }

        /// <summary>Physical tail logarithms remain finite when subtraction, standardization, or the shape product overflows.</summary>
        [TestMethod]
        [DataRow("GEV")]
        [DataRow("GPA")]
        public void TailTransformPreservesLogarithmsBeyondStandardizedRange(string family)
        {
            var rows = Rows(family).Where(row => row["case"] == "overflowing-tail-transform").ToList();
            Assert.IsNotEmpty(rows);
            foreach (var row in rows)
            {
                var distribution = Create(row);
                double x = Number(row, "x");
                double actual = row["quantity"] == "LogPDF" ? distribution.LogPDF(x)
                    : row["quantity"] == "LogCDF" ? distribution.LogCDF(x) : distribution.LogCCDF(x);
                OracleAssert(row, actual);
            }
        }

        /// <summary>Finite covariance-gradient contributions survive mismatched coordinate magnitudes.</summary>
        [TestMethod]
        public void WeibullVariancePreservesMismatchedCoordinateRanges()
        {
            var row = Rows("Weibull").Single(item => item["case"] == "MLE-mismatched-range-variance");
            var distribution = (Weibull)Create(row);
            OracleAssert(row, distribution.QuantileVariance(Number(row, "p"), (int)Number(row, "n"), ParameterEstimationMethod.MaximumLikelihood));
        }

        /// <summary>The Fisher residual stays positive when its leading trigamma term rounds to the reciprocal shape.</summary>
        [TestMethod]
        public void GammaFisherCovarianceRetainsLargeShapeResidual()
        {
            var rows = Rows("Gamma").Where(row => row["case"] == "MLE-large-shape-Fisher-residual").ToList();
            Assert.IsNotEmpty(rows);
            foreach (var row in rows)
            {
                var covariance = ((GammaDistribution)Create(row)).ParameterCovariance((int)Number(row, "n"), ParameterEstimationMethod.MaximumLikelihood);
                OracleAssert(row, covariance[row["quantity"][10] - '1', row["quantity"][11] - '1']);
            }
        }

        /// <summary>The positive mean-direction variance survives covariance cancellation at concentrated Gamma shapes.</summary>
        [TestMethod]
        [DataRow("MLE")]
        [DataRow("MoM")]
        public void GammaMedianVarianceRetainsLargeShapeMeanDirection(string method)
        {
            var row = Rows("Gamma").Single(item => item["case"] == method + "-large-shape-median-variance");
            var distribution = (GammaDistribution)Create(row);
            var estimationMethod = method == "MLE" ? ParameterEstimationMethod.MaximumLikelihood : ParameterEstimationMethod.MethodOfMoments;
            OracleAssert(row, distribution.QuantileVariance(Number(row, "p"), (int)Number(row, "n"), estimationMethod));
        }

        /// <summary>The median property restores physical scale before an unrepresentable unit quantile is rounded to zero.</summary>
        [TestMethod]
        public void WeibullMedianRetainsPhysicallyScaledTinyQuantile()
        {
            var row = Rows("Weibull").Single(item => item["case"] == "scale-rescued-median");
            var distribution = (Weibull)Create(row);
            OracleAssert(row, distribution.Median);
            OracleAssert(row, distribution.InverseCDF(.5));
        }

        /// <summary>GPA skewness and kurtosis remain finite when raw shape powers would overflow.</summary>
        [TestMethod]
        public void GpaHigherMomentsAvoidLargeShapePowerOverflow()
        {
            var rows = Rows("GPA").Where(row => row["case"] == "large-positive-shape-scaled-higher-moments").ToList();
            Assert.IsNotEmpty(rows);
            foreach (var row in rows)
            {
                var distribution = (GeneralizedPareto)Create(row);
                OracleAssert(row, row["quantity"] == "Skewness" ? distribution.Skewness : distribution.Kurtosis);
            }
        }

        /// <summary>Reads a family subset without a runtime R or network dependency.</summary>
        private static IEnumerable<Dictionary<string, string>> Rows(string family) => DistributionOracle.Read("extreme-positive.csv").Where(r => r["family"] == family);

        /// <summary>Creates the public-coordinate distribution used by an independent reference row.</summary>
        private static UnivariateDistributionBase Create(Dictionary<string, string> row) => Create(row["family"], Number(row, "xi"), Number(row, "scale"), Number(row, "shape"));

        /// <summary>Maps fixture families to their actual public constructors.</summary>
        private static UnivariateDistributionBase Create(string family, double xi, double scale, double shape) => family switch
        {
            "Exponential" => new Exponential(xi, scale), "Gamma" => new GammaDistribution(scale, shape),
            "GEV" => new GeneralizedExtremeValue(xi, scale, shape), "GPA" => new GeneralizedPareto(xi, scale, shape),
            "Gumbel" => new Gumbel(xi, scale), "Weibull" => new Weibull(scale, shape),
            _ => throw new ArgumentOutOfRangeException(nameof(family))
        };

        /// <summary>Parses a scalar reference field using the shared invariant reader.</summary>
        private static double Number(Dictionary<string, string> row, string field) => DistributionOracle.Number(row[field]);

        /// <summary>Identifies the independent failing coordinate in a test result.</summary>
        private static string Description(Dictionary<string, string> row) => string.Join("/", new[] { row["family"], row["case"], row["quantity"], row["scale"], row["shape"], row["x"], row["p"], row["n"] });

        /// <summary>Checks classification separately and uses the frozen absolute plus relative tolerance for finite values.</summary>
        private static void OracleAssert(Dictionary<string, string> row, double actual, double minimumRelativeTolerance = 0)
        {
            double expected = Number(row, "value");
            double tolerance = Number(row, "absolute_tolerance") + Math.Max(minimumRelativeTolerance, Number(row, "relative_tolerance")) * Math.Abs(expected);
            if (double.IsNaN(expected)) Assert.IsTrue(double.IsNaN(actual), Description(row));
            else if (double.IsInfinity(expected)) Assert.AreEqual(expected, actual, Description(row));
            else
            {
                Assert.IsTrue(!double.IsNaN(actual) && !double.IsInfinity(actual), Description(row) + " must be finite.");
                Assert.AreEqual(expected, actual, tolerance, Description(row));
            }
        }

        /// <summary>Protects nonzero representable linear probabilities and density classifications.</summary>
        private static void ScalarAssert(double expected, double actual, Dictionary<string, string> row, double relativeTolerance)
        {
            if (double.IsNaN(expected)) Assert.IsTrue(double.IsNaN(actual), Description(row));
            else if (double.IsInfinity(expected) || expected == 0) Assert.AreEqual(expected, actual, Description(row));
            else Assert.AreEqual(expected, actual, Math.Max(double.Epsilon, relativeTolerance * Math.Abs(expected)), Description(row));
        }

        /// <summary>Checks bounds against their returned start without manufacturing expected numerical estimates.</summary>
        private static void CheckConstraints(Tuple<double[], double[], double[]> constraints, string family)
        {
            for (int i = 0; i < constraints.Item1.Length; i++)
            {
                double initial = constraints.Item1[i], lower = constraints.Item2[i], upper = constraints.Item3[i];
                Assert.IsTrue(!double.IsNaN(initial) && !double.IsInfinity(initial), family + "/finite start");
                Assert.IsTrue(!double.IsNaN(lower) && !double.IsInfinity(lower) && !double.IsNaN(upper) && !double.IsInfinity(upper), family + "/finite bounds");
                Assert.IsTrue(lower < upper && lower <= initial && initial <= upper, family + "/ordered feasible bounds");
            }
        }
    }
}
