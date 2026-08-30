using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics;
using Numerics.Distributions;
using Numerics.Data.Statistics;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;
using System;
using System.Reflection;
using System.Xml.Linq;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Competing Risks distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://reliability.readthedocs.io/en/latest/Competing%20risk%20models.html" />
    /// </para>
    /// </remarks>

    [TestClass]
    public class Test_CompetingRisks
    {

        /// <summary>
        /// Test the moments of the Competing Risk model against the Python package 'reliability'.
        /// </summary>
        [TestMethod]
        public void Test_CR_Moments()
        {
            var d1 = new LogNormal(4, 0.1) { Base = Math.E };
            var d2 = new Weibull(50, 2);
            var d3 = new GammaDistribution(30, 1.5);
            var dists = new IUnivariateDistribution[] { d1, d2, d3 };
            var cr = new CompetingRisks(dists);

            // Numerics computes the moments using numerical integration
            Assert.AreEqual(27.0445, cr.Mean, 1E-2);
            Assert.AreEqual(25.0845, cr.Median, 1E-2);
            Assert.AreEqual(16.6581, cr.Mode, 1E-2);
            Assert.AreEqual(15.60225, cr.StandardDeviation, 1E-2);
            Assert.AreEqual(0.3371, cr.Skewness, 1E-2);
            Assert.AreEqual(-0.8719, cr.Kurtosis - 3, 1E-2); // Package reports excess kurtosis.

        }

        /// <summary>
        /// Test the CDF and Inverse CDF of the Competing Risk model against the Python package 'reliability'.
        /// </summary>
        [TestMethod]
        public void Test_CR_CDF()
        {
            var d1 = new LogNormal(4, 0.1) { Base = Math.E };
            var d2 = new Weibull(50, 2);
            var d3 = new GammaDistribution(30, 1.5);
            var dists = new IUnivariateDistribution[] { d1, d2, d3 };
            var cr = new CompetingRisks(dists);

            Assert.AreEqual(4.6431, cr.InverseCDF(0.05), 1E-3);
            Assert.AreEqual(25.0845, cr.InverseCDF(0.50), 1E-3);
            Assert.AreEqual(54.2056, cr.InverseCDF(0.95), 1E-3);

        }

        /// <summary>
        /// Verifies that the PDF does not return Infinity or NaN when evaluated at points
        /// in the left tail where CDF values approach zero under the maximum rule.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     For the maximum of independent random variables, the PDF formula involves
        ///     the ratio f_i(x) / F_i(x). When x is in the left tail, F_i(x) approaches zero,
        ///     causing division by zero if not handled properly.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     Create a competing risks model with two Normal distributions and evaluate
        ///     the PDF at a point approximately 5 standard deviations below the mean,
        ///     where CDF values will be on the order of 1E-7.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The PDF should return a small but finite non-negative value, not Infinity or NaN.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_PDF_MaxRule_SmallX_NoInfinity()
        {
            var dist1 = new Normal(100, 10);
            var dist2 = new Normal(110, 15);
            var cr = new CompetingRisks(new[] { dist1, dist2 });
            cr.MinimumOfRandomVariables = false; // Maximum rule

            // Test at a point where CDFs are very small
            double x = 50; // Far in left tail
            double pdf = cr.PDF(x);

            Assert.IsFalse(double.IsInfinity(pdf), "PDF should not be infinity");
            Assert.IsFalse(double.IsNaN(pdf), "PDF should not be NaN");
            Assert.IsGreaterThanOrEqualTo(0, pdf, "PDF should be non-negative");
        }

        /// <summary>
        /// Verifies that the LogPDF method returns finite values for all points in a
        /// randomly generated sample, ensuring numerical stability for log-likelihood calculations.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     MLE and Bayesian estimation methods work with log-likelihoods rather than
        ///     likelihoods to prevent numerical underflow. The LogPDF method must return
        ///     finite values (not NaN or +Infinity) for all valid inputs to ensure that
        ///     optimization algorithms can compute gradients and evaluate objective functions.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     <list type="number">
        ///         <item>Create a competing risks model with two Exponential distributions</item>
        ///         <item>Generate a random sample of 100 values</item>
        ///         <item>Verify LogPDF is finite for each sample point</item>
        ///         <item>Verify the sum (log-likelihood) is also finite</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     All LogPDF values should be finite negative numbers (since PDF ≤ 1 for
        ///     continuous distributions with unbounded support), and their sum should
        ///     be a finite negative number representing the log-likelihood.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_LogPDF_StableForMLE()
        {
            var dist1 = new Exponential(0.1);
            var dist2 = new Exponential(0.2);
            var cr = new CompetingRisks(new[] { dist1, dist2 });

            // Generate sample
            var sample = cr.GenerateRandomValues(100, 12345);

            // Verify log-likelihood is finite
            double logLik = 0;
            foreach (var x in sample)
            {
                double logPdf = cr.LogPDF(x);
                Assert.IsFalse(double.IsNaN(logPdf), $"LogPDF should not be NaN at x={x}");
                Assert.IsFalse(double.IsPositiveInfinity(logPdf), $"LogPDF should not be +Inf at x={x}");
                logLik += logPdf;
            }

            Assert.IsFalse(double.IsNaN(logLik), "Log-likelihood should not be NaN");
        }

        /// <summary>
        /// Verifies that the PDF integrates to approximately 1 over the support of the distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     A valid probability density function must integrate to 1 over its support.
        ///     This test verifies that the PDF implementation satisfies this fundamental
        ///     requirement for both minimum and maximum rules.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     Use numerical integration (e.g., adaptive quadrature) to compute the integral
        ///     of the PDF from the 1E-10 quantile to the 1-1E-10 quantile, which should
        ///     capture essentially all of the probability mass.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The integral should be within 0.001 of 1.0.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_PDF_IntegratesToOne()
        {
            var dist1 = new Normal(100, 10);
            var dist2 = new Normal(110, 15);

            // Test minimum rule
            var crMin = new CompetingRisks(new[] { dist1.Clone(), dist2.Clone() });
            crMin.MinimumOfRandomVariables = true;

            double lowerMin = crMin.InverseCDF(1E-10);
            double upperMin = crMin.InverseCDF(1 - 1E-10);
            var agkMin = new AdaptiveGaussKronrod(crMin.PDF, lowerMin, upperMin);
            agkMin.Integrate();
            double integralMin = agkMin.Result;
            Assert.AreEqual(1.0, integralMin, 0.001, "PDF (min rule) should integrate to 1");

            // Test maximum rule
            var crMax = new CompetingRisks(new[] { dist1.Clone(), dist2.Clone() });
            crMax.MinimumOfRandomVariables = false;

            double lowerMax = crMax.InverseCDF(1E-10);
            double upperMax = crMax.InverseCDF(1 - 1E-10);
            var agkMax = new AdaptiveGaussKronrod(crMax.PDF, lowerMax, upperMax);
            agkMax.Integrate();
            double integralMax = agkMax.Result;
            Assert.AreEqual(1.0, integralMax, 0.001, "PDF (max rule) should integrate to 1");
        }

        /// <summary>
        /// Verifies consistency between PDF and CDF by checking that the numerical
        /// derivative of the CDF equals the PDF at multiple points.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     By definition, f(x) = dF(x)/dx. This test verifies internal consistency
        ///     between the PDF and CDF implementations, which is critical for MLE where
        ///     both functions may be used.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     Evaluate both the PDF and the numerical derivative of the CDF at several
        ///     quantile points (0.1, 0.25, 0.5, 0.75, 0.9) and verify they match within
        ///     numerical tolerance.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The relative difference between PDF(x) and CDF'(x) should be less than 1E-4
        ///     at all test points.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_PDF_CDF_Consistency()
        {
            var dist1 = new Exponential(0.1);
            var dist2 = new Exponential(0.2);
            var cr = new CompetingRisks(new[] { dist1, dist2 });

            double[] quantiles = { 0.1, 0.25, 0.5, 0.75, 0.9 };

            foreach (double q in quantiles)
            {
                double x = cr.InverseCDF(q);
                double pdf = cr.PDF(x);
                double cdfDerivative = NumericalDerivative.Derivative(cr.CDF, x);

                double relError = Math.Abs(pdf - cdfDerivative) / Math.Max(pdf, 1E-10);
                Assert.IsLessThan(1E-4, relError, $"PDF and CDF derivative should match at quantile {q}. " +
                    $"PDF={pdf}, CDF'={cdfDerivative}, RelError={relError}");
            }
        }

        /// <summary>
        /// Verifies that LogPDF and log(PDF) return consistent values where both are numerically stable.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     The LogPDF method is implemented using log-space arithmetic for numerical stability.
        ///     This test verifies that LogPDF produces results consistent with log(PDF) in
        ///     regions where the standard PDF calculation is stable.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     Compare LogPDF(x) with Math.Log(PDF(x)) at the median (where both should be stable)
        ///     for both minimum and maximum rules.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The values should match within 1E-10 absolute tolerance at stable evaluation points.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_LogPDF_ConsistentWithPDF()
        {
            var dist1 = new Normal(100, 10);
            var dist2 = new Normal(110, 15);

            // Test minimum rule
            var crMin = new CompetingRisks(new[] { dist1.Clone(), dist2.Clone() });
            crMin.MinimumOfRandomVariables = true;
            double xMin = crMin.Median;
            double logPdfMin = crMin.LogPDF(xMin);
            double logOfPdfMin = Math.Log(crMin.PDF(xMin));
            Assert.AreEqual(logOfPdfMin, logPdfMin, 1E-10,
                "LogPDF should equal log(PDF) at median for min rule");

            // Test maximum rule  
            var crMax = new CompetingRisks(new[] { dist1.Clone(), dist2.Clone() });
            crMax.MinimumOfRandomVariables = false;
            double xMax = crMax.Median;
            double logPdfMax = crMax.LogPDF(xMax);
            double logOfPdfMax = Math.Log(crMax.PDF(xMax));
            Assert.AreEqual(logOfPdfMax, logPdfMax, 1E-10,
                "LogPDF should equal log(PDF) at median for max rule");
        }

        /// <summary>
        /// Verifies that dependent numerical-differentiation log-density returns log(0)
        /// when the density is outside the component support.
        /// </summary>
        [TestMethod]
        public void Test_LogPDF_DependentInvalidSupport_ReturnsNegativeInfinity()
        {
            var cr = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Exponential(1.0),
                new Exponential(2.0)
            })
            {
                Dependency = Numerics.Data.Statistics.Probability.DependencyType.PerfectlyPositive
            };

            Assert.AreEqual(double.NegativeInfinity, cr.LogPDF(-1.0));
        }


        #region Minimum Rule - 2 Distributions

        /// <summary>
        /// Verifies that dependency changes rebuild the cached Gaussian copula without mutating user correlation input.
        /// </summary>
        [TestMethod]
        public void Test_DependencyChangeInvalidatesMvnWithoutMutatingCorrelation()
        {
            var correlation = new[,] { { 1d, 0.5d }, { 0.5d, 1d } };
            var risks = new CompetingRisks(new IUnivariateDistribution[]
            {
                new Normal(),
                new Normal()
            })
            {
                CorrelationMatrix = (double[,])correlation.Clone(),
                Dependency = Probability.DependencyType.PerfectlyNegative,
                MinimumOfRandomVariables = true
            };

            double perfectlyNegative = risks.CDF(0d);

            for (int row = 0; row < 2; row++)
            {
                for (int column = 0; column < 2; column++)
                {
                    Assert.AreEqual(correlation[row, column], risks.CorrelationMatrix[row, column], 0d);
                }
            }

            risks.Dependency = Probability.DependencyType.CorrelationMatrix;
            double correlated = risks.CDF(0d);
            double expected = 0.75d - Math.Asin(0.5d) / (2d * Math.PI);

            Assert.AreEqual(1d, perfectlyNegative, 5E-5);
            Assert.AreEqual(expected, correlated, 1E-7);
            Assert.IsGreaterThan(0.1d, perfectlyNegative - correlated);
        }

        /// <summary>
        /// Verifies that the independent simulation path preserves its established seeded sequence.
        /// </summary>
        [TestMethod]
        public void Test_GenerateRandomValues_IndependentPreservesSeededSequence()
        {
            var risks = new CompetingRisks(new IUnivariateDistribution[]
            {
                new Normal(10d, 2d),
                new Normal(20d, 3d)
            })
            {
                Dependency = Probability.DependencyType.Independent,
                MinimumOfRandomVariables = false
            };
            double[] expected =
            {
                23.68205396527114d,
                16.630838190679725d,
                14.73954851356168d,
                22.820525389355137d,
                20.241482757909623d,
                25.128148621977083d,
                19.713254203467052d,
                24.292021772366837d,
                19.035542927284475d,
                16.953998432897098d
            };

            CollectionAssert.AreEqual(expected, risks.GenerateRandomValues(expected.Length, 12345));
        }

        /// <summary>
        /// Verifies that the public override honors every dependency mode by matching the
        /// explicit dependency-aware simulation entry point.
        /// </summary>
        /// <param name="dependency">The dependency mode to exercise.</param>
        [TestMethod]
        [DataRow(Probability.DependencyType.Independent)]
        [DataRow(Probability.DependencyType.PerfectlyPositive)]
        [DataRow(Probability.DependencyType.PerfectlyNegative)]
        [DataRow(Probability.DependencyType.CorrelationMatrix)]
        public void Test_GenerateRandomValues_MatchesDependencyAwarePath(Probability.DependencyType dependency)
        {
            var risks = new CompetingRisks(new IUnivariateDistribution[]
            {
                new Normal(10d, 2d),
                new Normal(20d, 3d)
            })
            {
                CorrelationMatrix = new[,] { { 1d, 0.6d }, { 0.6d, 1d } },
                Dependency = dependency,
                MinimumOfRandomVariables = false
            };

            double[] expected = risks.GenerateRandomValuesWithDependency(128, 24680);
            double[] actual = risks.GenerateRandomValues(128, 24680);

            CollectionAssert.AreEqual(expected, actual);
        }

        /// <summary>
        /// Verifies that correlation-matrix simulation rejects a matrix that is not
        /// positive definite before attempting Gaussian-copula sampling.
        /// </summary>
        [TestMethod]
        public void Test_GenerateRandomValues_CorrelationMatrixRejectsNonPositiveDefiniteMatrix()
        {
            var risks = new CompetingRisks(new IUnivariateDistribution[]
            {
                new Normal(),
                new Normal()
            })
            {
                CorrelationMatrix = new[,] { { 1d, 1d }, { 1d, 1d } },
                Dependency = Probability.DependencyType.CorrelationMatrix
            };

            ArgumentException exception = Assert.ThrowsExactly<ArgumentException>(
                () => risks.GenerateRandomValues(10, 12345));

            StringAssert.Contains(exception.Message, "positive definite");
        }

        private const int RECOVERY_SAMPLE_SIZE = 1000;
        private const int RECOVERY_SEED = 12345;

        /// <summary>
        /// Verifies MLE recovery for the identified independent minimum dog leg formed by
        /// Weibull(50, 1) and Weibull(80, 3).
        /// </summary>
        [TestMethod]
        public void Test_MLE_MinRule_2Dist_Weibull_DogLeg()
        {
            var parent = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Weibull(50d, 1d),
                new Weibull(80d, 3d)
            })
            {
                MinimumOfRandomVariables = true
            };
            var fitted = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Weibull(),
                new Weibull()
            })
            {
                MinimumOfRandomVariables = true
            };

            VerifyIdentifiedMleRecovery(parent, fitted, "independent two-Weibull minimum dog leg");
        }

        /// <summary>
        /// Verifies MLE recovery for the identified independent maximum dog leg formed by
        /// Weibull(100, 3) and Gumbel(80, 20).
        /// </summary>
        [TestMethod]
        public void Test_MLE_MaxRule_2Dist_Weibull_Gumbel_DogLeg()
        {
            var parent = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Weibull(100d, 3d),
                new Gumbel(80d, 20d)
            })
            {
                MinimumOfRandomVariables = false
            };
            var fitted = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Weibull(),
                new Gumbel()
            })
            {
                MinimumOfRandomVariables = false
            };

            VerifyIdentifiedMleRecovery(parent, fitted, "independent Weibull-Gumbel maximum dog leg");
        }

        /// <summary>
        /// Verifies MLE recovery for the identified fixed-correlation minimum dog leg formed by
        /// Weibull(50, 1) and Weibull(80, 3) at latent Gaussian correlation 0.6.
        /// </summary>
        [TestMethod]
        public void Test_MLE_CorrelatedMinRule_2Dist_Weibull_DogLeg()
        {
            double[,] correlation = { { 1d, 0.6d }, { 0.6d, 1d } };
            var parent = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Weibull(50d, 1d),
                new Weibull(80d, 3d)
            })
            {
                CorrelationMatrix = (double[,])correlation.Clone(),
                Dependency = Probability.DependencyType.CorrelationMatrix,
                MinimumOfRandomVariables = true
            };
            var fitted = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Weibull(),
                new Weibull()
            })
            {
                CorrelationMatrix = (double[,])correlation.Clone(),
                Dependency = Probability.DependencyType.CorrelationMatrix,
                MinimumOfRandomVariables = true
            };

            VerifyIdentifiedMleRecovery(parent, fitted, "fixed-correlation two-Weibull minimum dog leg");
        }

        #endregion

        #region Seed and Serialization

        /// <summary>
        /// Test that reseeding invalidates the lazily built multivariate-normal and
        /// empirical-CDF caches, that the seed survives cloning and the XML round-trip, and
        /// that undefined dependency ordinals are rejected on deserialization.
        /// </summary>
        [TestMethod]
        public void Test_PRNGSeed_InvalidatesCachesAndRoundTrips()
        {
            var distribution = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new Normal(0d, 1d),
                new Exponential(2d),
            }) { PRNGSeed = 2468 };

            FieldInfo mvnCreated = typeof(CompetingRisks).GetField("_mvnCreated", BindingFlags.Instance | BindingFlags.NonPublic)!;
            FieldInfo empiricalCreated = typeof(CompetingRisks).GetField("_empiricalCDFCreated", BindingFlags.Instance | BindingFlags.NonPublic)!;
            mvnCreated.SetValue(distribution, true);
            empiricalCreated.SetValue(distribution, true);
            distribution.PRNGSeed = 1357;
            Assert.IsFalse((bool)mvnCreated.GetValue(distribution)!);
            Assert.IsFalse((bool)empiricalCreated.GetValue(distribution)!);

            var clone = (CompetingRisks)distribution.Clone();
            Assert.AreEqual(1357, clone.PRNGSeed);
            var restored = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(distribution.ToXElement());
            Assert.AreEqual(1357, restored.PRNGSeed);

            XElement malformed = distribution.ToXElement();
            malformed.SetAttributeValue(nameof(CompetingRisks.Dependency), "999");
            Assert.Throws<ArgumentException>(() => CompetingRisks.FromXElement(malformed));
        }

        #endregion

        #region Helper Methods

        /// <summary>
        /// Fits one predeclared competing-risk parent and requires central-95-percent coordinate
        /// recovery from the full-likelihood observed-information covariance.
        /// </summary>
        /// <param name="parent">The known generating competing-risk distribution.</param>
        /// <param name="fitted">The fresh competing-risk distribution to estimate.</param>
        /// <param name="label">The scientific fixture label.</param>
        private static void VerifyIdentifiedMleRecovery(
            CompetingRisks parent,
            CompetingRisks fitted,
            string label)
        {
            double[] sample = parent.GenerateRandomValues(RECOVERY_SAMPLE_SIZE, RECOVERY_SEED);
            AssertIdentifiableRecoveryDesign(parent, sample, label);

            double[] rawParameters = fitted.MLE(sample);
            Assert.IsFalse(
                rawParameters.Any(parameter => !Tools.IsFinite(parameter)),
                $"{label}: every fitted coordinate must be finite.");

            Matrix rawCovariance = ComputeObservedInformationCovariance(
                fitted,
                sample,
                rawParameters,
                label);
            (double[] parameters, Matrix covariance) = CanonicalizeRecoveryCoordinates(
                parent,
                rawParameters,
                rawCovariance);
            double[] truth = parent.GetParameters;
            var evidence = new List<string>(truth.Length);
            var standardizedErrors = new double[truth.Length];
            for (int parameterIndex = 0; parameterIndex < truth.Length; parameterIndex++)
            {
                double standardError = Math.Sqrt(covariance[parameterIndex, parameterIndex]);
                standardizedErrors[parameterIndex] =
                    Math.Abs(parameters[parameterIndex] - truth[parameterIndex]) / standardError;
                evidence.Add(
                    $"coordinate {parameterIndex + 1}: fit={parameters[parameterIndex]:G8}, " +
                    $"truth={truth[parameterIndex]:G8}, SE={standardError:G8}, " +
                    $"|z|={standardizedErrors[parameterIndex]:G6}");
            }

            Assert.IsTrue(
                standardizedErrors.All(error => error <= 1.96d),
                $"{label}: one or more identified MLE coordinates missed the central 95% " +
                $"observed-information interval. fitted logL={EvaluateLogLikelihood(fitted, sample, parameters):G12}, " +
                $"parent logL={parent.LogLikelihood(sample):G12}. {string.Join("; ", evidence)}");
        }

        /// <summary>
        /// Computes a symmetric observed-information covariance from the complete competing-risk
        /// sample likelihood without adding a ridge or changing production estimation.
        /// </summary>
        /// <param name="fitted">The fitted distribution defining the model families and rule.</param>
        /// <param name="sample">The generated composite observations.</param>
        /// <param name="parameters">The raw MLE coordinates.</param>
        /// <param name="label">The scientific fixture label.</param>
        /// <returns>The inverse observed-information covariance.</returns>
        private static Matrix ComputeObservedInformationCovariance(
            CompetingRisks fitted,
            IList<double> sample,
            double[] parameters,
            string label)
        {
            Tuple<double[], double[], double[]> constraints = fitted.GetParameterConstraints(sample);
            double LogLikelihood(double[] candidate) => EvaluateLogLikelihood(fitted, sample, candidate);

            Matrix rawInformation = ComputeRecoveryHessian(
                LogLikelihood,
                parameters,
                constraints.Item2,
                constraints.Item3) * -1d;
            var information = new Matrix(parameters.Length, parameters.Length);
            var scaledInformation = new Matrix(parameters.Length, parameters.Length);
            for (int row = 0; row < parameters.Length; row++)
            {
                double rowScale = Math.Max(1d, Math.Abs(parameters[row]));
                for (int column = 0; column < parameters.Length; column++)
                {
                    double value = 0.5d * (rawInformation[row, column] + rawInformation[column, row]);
                    information[row, column] = value;
                    double columnScale = Math.Max(1d, Math.Abs(parameters[column]));
                    scaledInformation[row, column] = value * rowScale * columnScale;
                }
            }

            var scaledDecomposition = new SingularValueDecomposition(scaledInformation);
            Assert.AreEqual(
                parameters.Length,
                scaledDecomposition.Rank(),
                $"{label}: scale-normalized observed information is rank deficient; " +
                $"inverse condition={scaledDecomposition.InverseCondition:G6}.");
            Assert.IsTrue(
                Tools.IsFinite(scaledDecomposition.InverseCondition) &&
                scaledDecomposition.InverseCondition > 0d,
                $"{label}: scale-normalized observed-information condition is unusable.");

            CholeskyDecomposition cholesky;
            try
            {
                cholesky = new CholeskyDecomposition(information);
            }
            catch (Exception exception)
            {
                Assert.Fail($"{label}: observed information is not positive definite: {exception.Message}");
                throw;
            }
            Matrix covariance = cholesky.InverseA();
            for (int parameterIndex = 0; parameterIndex < parameters.Length; parameterIndex++)
            {
                Assert.IsTrue(
                    Tools.IsFinite(covariance[parameterIndex, parameterIndex]) &&
                    covariance[parameterIndex, parameterIndex] > 0d,
                    $"{label}: coordinate {parameterIndex + 1} covariance diagonal is invalid.");
            }
            return covariance;
        }

        /// <summary>
        /// Computes a bounded central-difference Hessian for the test-only uncertainty oracle.
        /// </summary>
        /// <param name="function">The scalar log-likelihood function.</param>
        /// <param name="parameters">The coordinate vector at which to differentiate.</param>
        /// <param name="lowerBounds">The coordinate lower bounds.</param>
        /// <param name="upperBounds">The coordinate upper bounds.</param>
        /// <returns>The symmetric finite-difference Hessian.</returns>
        private static Matrix ComputeRecoveryHessian(
            Func<double[], double> function,
            double[] parameters,
            IReadOnlyList<double> lowerBounds,
            IReadOnlyList<double> upperBounds)
        {
            int parameterCount = parameters.Length;
            var steps = new double[parameterCount];
            for (int parameterIndex = 0; parameterIndex < parameterCount; parameterIndex++)
            {
                double nominal = 1E-4 * (Math.Abs(parameters[parameterIndex]) + 1d);
                double leftRoom = parameters[parameterIndex] - lowerBounds[parameterIndex];
                double rightRoom = upperBounds[parameterIndex] - parameters[parameterIndex];
                steps[parameterIndex] = Math.Min(nominal, 0.25d * Math.Min(leftRoom, rightRoom));
                Assert.IsGreaterThan(0d, steps[parameterIndex],
                    $"Hessian step for coordinate {parameterIndex + 1} must be positive.");
            }

            var hessian = new Matrix(parameterCount, parameterCount);
            double centerValue = function((double[])parameters.Clone());
            for (int row = 0; row < parameterCount; row++)
            {
                double[] forward = (double[])parameters.Clone();
                double[] backward = (double[])parameters.Clone();
                forward[row] += steps[row];
                backward[row] -= steps[row];
                hessian[row, row] =
                    (function(forward) - 2d * centerValue + function(backward)) /
                    (steps[row] * steps[row]);

                for (int column = row + 1; column < parameterCount; column++)
                {
                    double[] plusPlus = (double[])parameters.Clone();
                    double[] plusMinus = (double[])parameters.Clone();
                    double[] minusPlus = (double[])parameters.Clone();
                    double[] minusMinus = (double[])parameters.Clone();
                    plusPlus[row] += steps[row];
                    plusPlus[column] += steps[column];
                    plusMinus[row] += steps[row];
                    plusMinus[column] -= steps[column];
                    minusPlus[row] -= steps[row];
                    minusPlus[column] += steps[column];
                    minusMinus[row] -= steps[row];
                    minusMinus[column] -= steps[column];
                    double mixed =
                        (function(plusPlus) - function(plusMinus) -
                         function(minusPlus) + function(minusMinus)) /
                        (4d * steps[row] * steps[column]);
                    hessian[row, column] = mixed;
                    hessian[column, row] = mixed;
                }
            }
            return hessian;
        }

        /// <summary>
        /// Evaluates the complete sample likelihood on a clone so finite differences cannot leave
        /// the fitted test object in a perturbed coordinate state.
        /// </summary>
        /// <param name="template">The fitted distribution template.</param>
        /// <param name="sample">The observed composite sample.</param>
        /// <param name="parameters">The flattened component coordinates.</param>
        /// <returns>The complete competing-risk log likelihood.</returns>
        private static double EvaluateLogLikelihood(
            CompetingRisks template,
            IList<double> sample,
            double[] parameters)
        {
            var candidate = (CompetingRisks)template.Clone();
            candidate.SetParameters(parameters);
            return candidate.LogLikelihood(sample);
        }

        /// <summary>
        /// Applies the predeclared increasing-Weibull-shape label rule to a point and covariance.
        /// </summary>
        /// <param name="parent">The generating distribution declaring the component families.</param>
        /// <param name="parameters">The raw flattened coordinates.</param>
        /// <param name="covariance">The raw coordinate covariance.</param>
        /// <returns>The canonically ordered coordinates and covariance.</returns>
        private static (double[] Parameters, Matrix Covariance) CanonicalizeRecoveryCoordinates(
            CompetingRisks parent,
            double[] parameters,
            Matrix covariance)
        {
            int[] order = GetRecoveryCoordinateOrder(parent, parameters);
            var orderedParameters = new double[parameters.Length];
            var orderedCovariance = new Matrix(parameters.Length, parameters.Length);
            for (int row = 0; row < parameters.Length; row++)
            {
                orderedParameters[row] = parameters[order[row]];
                for (int column = 0; column < parameters.Length; column++)
                    orderedCovariance[row, column] = covariance[order[row], order[column]];
            }
            return (orderedParameters, orderedCovariance);
        }

        /// <summary>
        /// Gets the scientific component-coordinate order, sorting same-family Weibulls by shape.
        /// </summary>
        /// <param name="parent">The generating distribution declaring the component families.</param>
        /// <param name="parameters">The raw flattened coordinates.</param>
        /// <returns>For each canonical coordinate, the matching raw coordinate index.</returns>
        private static int[] GetRecoveryCoordinateOrder(
            CompetingRisks parent,
            IReadOnlyList<double> parameters)
        {
            if (!parent.Distributions.All(distribution => distribution is Weibull) ||
                parent.Distributions.Count < 2)
            {
                return Enumerable.Range(0, parameters.Count).ToArray();
            }

            const int WeibullParameterCount = 2;
            Assert.HasCount(
                WeibullParameterCount * parent.Distributions.Count,
                parameters);
            return Enumerable.Range(0, parent.Distributions.Count)
                .OrderBy(componentIndex => parameters[WeibullParameterCount * componentIndex + 1])
                .SelectMany(componentIndex => new[]
                {
                    WeibullParameterCount * componentIndex,
                    WeibullParameterCount * componentIndex + 1
                })
                .ToArray();
        }

        /// <summary>
        /// Requires balanced theoretical and realized cause shares plus visible dog-leg crossovers
        /// before the Numerics MLE is called.
        /// </summary>
        /// <param name="parent">The generating competing-risk distribution.</param>
        /// <param name="sample">The production generated sample.</param>
        /// <param name="label">The scientific fixture label.</param>
        private static void AssertIdentifiableRecoveryDesign(
            CompetingRisks parent,
            IReadOnlyList<double> sample,
            string label)
        {
            const int probabilityCount = 1000;
            int componentCount = parent.Distributions.Count;
            (double[] labeledSample, int[] hardWinnerCounts) = GenerateLabeledRecoverySample(parent);
            CollectionAssert.AreEqual(sample.ToArray(), labeledSample,
                $"{label}: verification-only labels must reproduce production generation exactly.");

            var theoreticalShares = new double[componentCount];
            var softCounts = new double[componentCount];
            var dominanceCounts = new int[componentCount];
            var dominantSequence = new int[probabilityCount];
            foreach (double observation in sample)
            {
                double[] responsibilities = ComputeRecoveryResponsibilities(parent, observation, label);
                for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                    softCounts[componentIndex] += responsibilities[componentIndex];
            }

            for (int probabilityIndex = 0; probabilityIndex < probabilityCount; probabilityIndex++)
            {
                double probability = (probabilityIndex + 0.5d) / probabilityCount;
                double[] responsibilities = ComputeRecoveryResponsibilities(
                    parent,
                    parent.InverseCDF(probability),
                    label);
                dominantSequence[probabilityIndex] = IndexOfLargest(responsibilities);
                for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                {
                    theoreticalShares[componentIndex] += responsibilities[componentIndex] / probabilityCount;
                    if (responsibilities[componentIndex] >= 0.5d)
                        dominanceCounts[componentIndex]++;
                }
            }

            double[] crossovers = FindRecoveryCrossovers(dominantSequence, probabilityCount);
            double[] interiorCrossovers = crossovers
                .Where(probability => probability >= 0.10d && probability <= 0.90d)
                .ToArray();
            for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
            {
                double dominanceMass = dominanceCounts[componentIndex] / (double)probabilityCount;
                Assert.IsGreaterThanOrEqualTo(0.15d, theoreticalShares[componentIndex],
                    $"{label}: component {componentIndex + 1} theoretical share " +
                    $"{theoreticalShares[componentIndex]:P2} is below 15%.");
                Assert.IsGreaterThanOrEqualTo(100, hardWinnerCounts[componentIndex],
                    $"{label}: component {componentIndex + 1} has only " +
                    $"{hardWinnerCounts[componentIndex]} hard wins.");
                Assert.IsGreaterThanOrEqualTo(100d, softCounts[componentIndex],
                    $"{label}: component {componentIndex + 1} soft event count " +
                    $"{softCounts[componentIndex]:F1} is below 100.");
                Assert.IsGreaterThanOrEqualTo(0.10d, dominanceMass,
                    $"{label}: component {componentIndex + 1} owns only {dominanceMass:P2} " +
                    "of the composite probability scale.");
            }
            Assert.HasCount(componentCount - 1, interiorCrossovers,
                $"{label}: expected {componentCount - 1} interior dog-leg crossovers, found " +
                $"{interiorCrossovers.Length} inside [0.10, 0.90]; all crossovers are " +
                $"[{string.Join(", ", crossovers.Select(value => value.ToString("F3")))}].");

            Console.WriteLine(
                $"{label}: theoretical shares [{string.Join(", ", theoreticalShares.Select(value => value.ToString("P1")))}], " +
                $"hard wins [{string.Join(", ", hardWinnerCounts)}], soft counts " +
                $"[{string.Join(", ", softCounts.Select(value => value.ToString("F1")))}], " +
                $"crossovers [{string.Join(", ", crossovers.Select(value => value.ToString("F3")))}].");
        }

        /// <summary>
        /// Reproduces production random-number ordering while retaining the latent winning cause.
        /// </summary>
        /// <param name="parent">The generating competing-risk distribution.</param>
        /// <returns>The composite sample and component hard winner counts.</returns>
        private static (double[] Sample, int[] HardWinnerCounts) GenerateLabeledRecoverySample(
            CompetingRisks parent)
        {
            int componentCount = parent.Distributions.Count;
            var sample = new double[RECOVERY_SAMPLE_SIZE];
            var counts = new int[componentCount];
            if (parent.Dependency == Probability.DependencyType.Independent)
            {
                var random = new MersenneTwister(RECOVERY_SEED);
                for (int observationIndex = 0; observationIndex < RECOVERY_SAMPLE_SIZE; observationIndex++)
                {
                    var values = new double[componentCount];
                    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                    {
                        values[componentIndex] = parent.Distributions[componentIndex]
                            .InverseCDF(random.NextDouble());
                    }
                    RecordRecoveryWinner(parent, values, sample, counts, observationIndex);
                }
            }
            else
            {
                Assert.AreEqual(Probability.DependencyType.CorrelationMatrix, parent.Dependency);
                var multivariateNormal = new MultivariateNormal(
                    new double[componentCount],
                    parent.CorrelationMatrix);
                double[,] normals = multivariateNormal.GenerateRandomValues(
                    RECOVERY_SAMPLE_SIZE,
                    RECOVERY_SEED);
                for (int observationIndex = 0; observationIndex < RECOVERY_SAMPLE_SIZE; observationIndex++)
                {
                    var values = new double[componentCount];
                    for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                    {
                        values[componentIndex] = parent.Distributions[componentIndex].InverseCDF(
                            Normal.StandardCDF(normals[observationIndex, componentIndex]));
                    }
                    RecordRecoveryWinner(parent, values, sample, counts, observationIndex);
                }
            }
            return (sample, counts);
        }

        /// <summary>
        /// Records one generated composite observation and hard winning component.
        /// </summary>
        /// <param name="parent">The generating competing-risk distribution.</param>
        /// <param name="values">The latent component values.</param>
        /// <param name="sample">The composite sample under construction.</param>
        /// <param name="counts">The hard winner counts.</param>
        /// <param name="observationIndex">The observation being recorded.</param>
        private static void RecordRecoveryWinner(
            CompetingRisks parent,
            IReadOnlyList<double> values,
            double[] sample,
            int[] counts,
            int observationIndex)
        {
            int winner = 0;
            for (int componentIndex = 1; componentIndex < values.Count; componentIndex++)
            {
                bool replace = parent.MinimumOfRandomVariables
                    ? values[componentIndex] < values[winner]
                    : values[componentIndex] > values[winner];
                if (replace)
                    winner = componentIndex;
            }
            sample[observationIndex] = values[winner];
            counts[winner]++;
        }

        /// <summary>
        /// Computes component cause responsibilities at one observed composite value.
        /// </summary>
        /// <param name="parent">The generating competing-risk distribution.</param>
        /// <param name="location">The observed composite value.</param>
        /// <param name="label">The scientific fixture label.</param>
        /// <returns>The normalized cause responsibility vector.</returns>
        private static double[] ComputeRecoveryResponsibilities(
            CompetingRisks parent,
            double location,
            string label)
        {
            int componentCount = parent.Distributions.Count;
            var contributions = new double[componentCount];
            if (parent.Dependency == Probability.DependencyType.Independent)
            {
                for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                {
                    double contribution = parent.Distributions[componentIndex].PDF(location);
                    for (int otherIndex = 0; otherIndex < componentCount; otherIndex++)
                    {
                        if (componentIndex == otherIndex)
                            continue;
                        contribution *= parent.MinimumOfRandomVariables
                            ? parent.Distributions[otherIndex].CCDF(location)
                            : parent.Distributions[otherIndex].CDF(location);
                    }
                    contributions[componentIndex] = contribution;
                }
            }
            else
            {
                Assert.AreEqual(2, componentCount,
                    $"{label}: correlated responsibilities require two components.");
                Assert.IsTrue(parent.MinimumOfRandomVariables,
                    $"{label}: only the approved correlated minimum is supported.");
                double correlation = parent.CorrelationMatrix[0, 1];
                double conditionalScale = Math.Sqrt(1d - correlation * correlation);
                double[] probabilities = parent.Distributions
                    .Select(distribution => Tools.Clamp(
                        distribution.CDF(location),
                        1E-14,
                        1d - 1E-14))
                    .ToArray();
                double[] normals = probabilities.Select(Normal.StandardZ).ToArray();
                for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                {
                    int otherIndex = 1 - componentIndex;
                    double threshold = (normals[otherIndex] - correlation * normals[componentIndex]) /
                        conditionalScale;
                    contributions[componentIndex] = parent.Distributions[componentIndex].PDF(location) *
                        (1d - Normal.StandardCDF(threshold));
                }
            }

            double total = contributions.Sum();
            Assert.IsTrue(Tools.IsFinite(total) && total > 0d,
                $"{label}: cause contributions are invalid at x={location:G17}.");
            for (int componentIndex = 0; componentIndex < componentCount; componentIndex++)
                contributions[componentIndex] /= total;
            return contributions;
        }

        /// <summary>
        /// Returns the first index containing the largest value.
        /// </summary>
        /// <param name="values">The values to compare.</param>
        /// <returns>The largest-value index.</returns>
        private static int IndexOfLargest(IReadOnlyList<double> values)
        {
            int largest = 0;
            for (int index = 1; index < values.Count; index++)
            {
                if (values[index] > values[largest])
                    largest = index;
            }
            return largest;
        }

        /// <summary>
        /// Finds midpoint-grid probabilities where the dominant cause changes.
        /// </summary>
        /// <param name="dominantSequence">The dominant component over the probability grid.</param>
        /// <param name="probabilityCount">The number of grid midpoints.</param>
        /// <returns>The ordered crossover probabilities.</returns>
        private static double[] FindRecoveryCrossovers(
            IReadOnlyList<int> dominantSequence,
            int probabilityCount)
        {
            var crossovers = new List<double>();
            int previous = dominantSequence[0];
            for (int index = 1; index < dominantSequence.Count; index++)
            {
                if (dominantSequence[index] == previous)
                    continue;
                crossovers.Add(index / (double)probabilityCount);
                previous = dominantSequence[index];
            }
            return crossovers.ToArray();
        }

        /// <summary>
        /// Test that the empirical inverse CDF resolves the quantiles of heavy-tailed, negative-support
        /// and positive-support components to within half a percent of a root-solved inversion for
        /// both x-transforms, because the empirical grid is log-spaced on the offset axis regardless
        /// of the transform.
        /// </summary>
        [TestMethod]
        public void Test_EmpiricalInverseCDF_ResolvesQuantilesRegardlessOfXTransform()
        {
            AssertEmpiricalInverseMatchesRootSolve(
                new CompetingRisks(new UnivariateDistributionBase[] { new GeneralizedExtremeValue(100, 20, -0.2), new GeneralizedExtremeValue(130, 25, -0.15) })
                {
                    MinimumOfRandomVariables = false,
                    XTransform = Transform.None
                },
                new[] { 0.5, 0.9, 0.99, 0.999 },
                "heavy-tailed maxima, no transform");
            AssertEmpiricalInverseMatchesRootSolve(
                new CompetingRisks(new UnivariateDistributionBase[] { new Normal(0, 1), new Normal(5, 2) })
                {
                    MinimumOfRandomVariables = false,
                    XTransform = Transform.None
                },
                new[] { 0.01, 0.5, 0.99 },
                "negative support, no transform");
            AssertEmpiricalInverseMatchesRootSolve(
                new CompetingRisks(new UnivariateDistributionBase[] { new GeneralizedPareto(0, 10, -0.1), new Exponential(0, 8) })
                {
                    MinimumOfRandomVariables = true,
                    XTransform = Transform.Logarithmic
                },
                new[] { 0.1, 0.5, 0.9, 0.99 },
                "positive support, logarithmic transform");
        }

        /// <summary>
        /// Asserts the empirical inverse CDF of a competing-risks model matches the root-solved
        /// inverse of an identical model without an empirical CDF.
        /// </summary>
        /// <param name="competingRisks">The model to evaluate through its empirical CDF.</param>
        /// <param name="probabilities">The non-exceedance probabilities to check.</param>
        /// <param name="context">The assertion context.</param>
        private static void AssertEmpiricalInverseMatchesRootSolve(CompetingRisks competingRisks, double[] probabilities, string context)
        {
            var reference = (CompetingRisks)competingRisks.Clone();
            competingRisks.CreateEmpiricalCDF();
            foreach (double probability in probabilities)
            {
                double expected = reference.InverseCDF(probability);
                double actual = competingRisks.InverseCDF(probability);
                Assert.AreEqual(expected, actual, Math.Abs(expected) * 0.005 + 1E-9, $"{context}: p = {probability}");
            }
        }

        #endregion

        /// <summary>
        /// The empirical machinery under the hood is unchanged by the extrapolation property:
        /// the cumulative incidence functions are built at the default policy, the empirical CDF
        /// keeps its default far-tail endpoint hold, and no extrapolation attribute appears in
        /// the serialized form.
        /// </summary>
        [TestMethod]
        public void Test_CompetingRisks_EmpiricalUnderTheHood_NoRegression()
        {
            var cr = new CompetingRisks(new UnivariateDistributionBase[] { new Normal(10, 2), new Normal(12, 3) });
            var cifs = cr.CumulativeIncidenceFunctions();
            for (int i = 0; i < cifs.Count; i++)
            {
                Assert.AreEqual(ExtrapolationSides.None, cifs[i].Extrapolation);
            }

            cr.CreateEmpiricalCDF();
            // 1 - 1E-17 rounds to exactly 1.0, which takes the container's own support guard —
            // so the far-tail hold is probed at the guard floor and the last double below one.
            double low = cr.InverseCDF(1E-17);
            double high = cr.InverseCDF(1d - 1E-16);
            Assert.IsFalse(double.IsNaN(low) || double.IsInfinity(low));
            Assert.IsFalse(double.IsNaN(high) || double.IsInfinity(high));
            Assert.AreEqual(low, cr.InverseCDF(1E-16), 0d);
            Assert.DoesNotContain("Extrapolation", cr.ToXElement().ToString());
        }
    }
}
