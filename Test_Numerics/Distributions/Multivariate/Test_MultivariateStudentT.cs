using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Sampling;

namespace Distributions.Multivariate
{
    /// <summary>
    /// Unit tests for the Multivariate Student's t-distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// Reference values verified against Python scipy.stats.multivariate_t (v1.11+).
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_MultivariateStudentT
    {

        #region Test Parameters

        /// <summary>
        /// Standard test case: 2D MVT with ν=5, μ=(1,2), Σ=((1,0.5),(0.5,2)).
        /// </summary>
        private static MultivariateStudentT CreateStandard2D()
        {
            return new MultivariateStudentT(
                5.0,
                new[] { 1.0, 2.0 },
                new[,] { { 1.0, 0.5 }, { 0.5, 2.0 } });
        }

        #endregion

        #region Properties Tests

        /// <summary>
        /// Verify basic distribution properties: Mean, Mode, Median, Variance, StandardDeviation, Covariance.
        /// Verified against analytical formulas.
        /// </summary>
        [TestMethod]
        public void Test_Properties()
        {
            var mvt = CreateStandard2D();

            // Mean = location for ν > 1
            Assert.AreEqual(1.0, mvt.Mean[0], 1E-10);
            Assert.AreEqual(2.0, mvt.Mean[1], 1E-10);

            // Mode = location
            Assert.AreEqual(1.0, mvt.Mode[0], 1E-10);
            Assert.AreEqual(2.0, mvt.Mode[1], 1E-10);

            // Median = location
            Assert.AreEqual(1.0, mvt.Median[0], 1E-10);
            Assert.AreEqual(2.0, mvt.Median[1], 1E-10);

            // Variance = ν/(ν-2) · diag(Σ) = 5/3 · [1, 2] = [1.6667, 3.3333]
            Assert.AreEqual(1.6666666666666667, mvt.Variance[0], 1E-10);
            Assert.AreEqual(3.3333333333333335, mvt.Variance[1], 1E-10);

            // StandardDeviation
            Assert.AreEqual(Math.Sqrt(5.0 / 3.0), mvt.StandardDeviation[0], 1E-10);
            Assert.AreEqual(Math.Sqrt(10.0 / 3.0), mvt.StandardDeviation[1], 1E-10);

            // Covariance = ν/(ν-2) · Σ
            var cov = mvt.Covariance;
            Assert.AreEqual(1.6666666666666667, cov[0, 0], 1E-10);
            Assert.AreEqual(0.8333333333333333, cov[0, 1], 1E-10);
            Assert.AreEqual(0.8333333333333333, cov[1, 0], 1E-10);
            Assert.AreEqual(3.3333333333333335, cov[1, 1], 1E-10);

            // Other properties
            Assert.AreEqual(5.0, mvt.DegreesOfFreedom, 1E-10);
            Assert.AreEqual(2, mvt.Dimension);
            Assert.IsTrue(mvt.ParametersValid);
            Assert.IsTrue(mvt.IsPositiveDefinite);
            Assert.AreEqual("Multivariate Student's t", mvt.DisplayName);
            Assert.AreEqual("Multi-T", mvt.ShortDisplayName);
            Assert.AreEqual(MultivariateDistributionType.MultivariateStudentT, mvt.Type);
        }

        /// <summary>
        /// Verify that Mean throws for ν ≤ 1 and Variance throws for ν ≤ 2.
        /// </summary>
        [TestMethod]
        public void Test_UndefinedMoments()
        {
            // ν = 1: Mean undefined
            var mvt1 = new MultivariateStudentT(2, 1.0);
            Assert.Throws<InvalidOperationException>(() => { var _ = mvt1.Mean; });

            // ν = 2: Variance undefined
            var mvt2 = new MultivariateStudentT(2, 2.0);
            Assert.Throws<InvalidOperationException>(() => { var _ = mvt2.Variance; });
            Assert.Throws<InvalidOperationException>(() => { var _ = mvt2.Covariance; });

            // ν = 1: Variance also undefined
            Assert.Throws<InvalidOperationException>(() => { var _ = mvt1.Variance; });

            // ν = 3: Both defined
            var mvt3 = new MultivariateStudentT(2, 3.0);
            Assert.IsNotNull(mvt3.Mean);
            Assert.IsNotNull(mvt3.Variance);
        }

        #endregion

        #region PDF Tests

        /// <summary>
        /// Verify PDF against Python scipy.stats.multivariate_t for 2D case.
        /// MVT(ν=5, μ=(1,2), Σ=((1,0.5),(0.5,2))).
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var mvt = CreateStandard2D();

            // At mode (x = μ)
            double pdf_mode = mvt.PDF(new[] { 1.0, 2.0 });
            Assert.AreEqual(0.12030982838508346, pdf_mode, 1E-12);

            // Away from mode
            double pdf_00 = mvt.PDF(new[] { 0.0, 0.0 });
            Assert.AreEqual(0.032213925218971075, pdf_00, 1E-12);

            double pdf_34 = mvt.PDF(new[] { 3.0, 4.0 });
            Assert.AreEqual(0.012395882561151605, pdf_34, 1E-12);

            double pdf_25 = mvt.PDF(new[] { 2.0, 5.0 });
            Assert.AreEqual(0.012395882561151605, pdf_25, 1E-12);

            double pdf_37 = mvt.PDF(new[] { 3.0, 7.0 });
            Assert.AreEqual(0.0013219841225890843, pdf_37, 1E-12);
        }

        /// <summary>
        /// Verify LogPDF is consistent with log(PDF).
        /// </summary>
        [TestMethod]
        public void Test_LogPDF()
        {
            var mvt = CreateStandard2D();

            var testPoints = new[]
            {
                new[] { 1.0, 2.0 },
                new[] { 0.0, 0.0 },
                new[] { 3.0, 4.0 },
                new[] { -1.0, 5.0 }
            };

            foreach (var x in testPoints)
            {
                double pdf = mvt.PDF(x);
                double logpdf = mvt.LogPDF(x);
                Assert.AreEqual(Math.Log(pdf), logpdf, 1E-10,
                    $"LogPDF inconsistent with log(PDF) at ({x[0]}, {x[1]})");
            }

            // Known values
            Assert.AreEqual(-2.1176849603770576, mvt.LogPDF(new[] { 1.0, 2.0 }), 1E-10);
            Assert.AreEqual(-3.4353564596992499, mvt.LogPDF(new[] { 0.0, 0.0 }), 1E-10);
        }

        /// <summary>
        /// Verify that 1D MVT PDF matches univariate Student's t PDF.
        /// MVT(ν=5, μ=[3], Σ=[[4]]) should equal StudentT(μ=3, σ=2, ν=5).
        /// </summary>
        [TestMethod]
        public void Test_PDF_ReducesToUnivariate()
        {
            double nu = 5.0;
            double mu = 3.0;
            double sigma = 2.0;

            // Note: The Numerics StudentT.PDF returns the standard t density evaluated at
            // Z = (x-μ)/σ, without the 1/σ Jacobian factor. So:
            //   StudentT.PDF(x) = σ · MVT(ν, [μ], [[σ²]]).PDF(x)
            var mvt = new MultivariateStudentT(nu, new[] { mu }, new[,] { { sigma * sigma } });
            var univT = new StudentT(mu, sigma, nu);

            double[] testX = { -5.0, -1.0, 0.0, 3.0, 5.0, 10.0 };
            foreach (double x in testX)
            {
                double mvtPdf = mvt.PDF(new[] { x });
                double univPdf = univT.PDF(x);
                Assert.AreEqual(univPdf / sigma, mvtPdf, 1E-12,
                    $"1D MVT PDF does not match univariate StudentT PDF at x={x}");
            }
        }

        #endregion

        #region Mahalanobis Tests

        /// <summary>
        /// Verify Mahalanobis distance computation against analytical values.
        /// For Σ=((1,0.5),(0.5,2)): Σ⁻¹ = (1/1.75)·((2,-0.5),(-0.5,1)).
        /// </summary>
        [TestMethod]
        public void Test_Mahalanobis()
        {
            var mvt = CreateStandard2D();

            // At mode: Q = 0
            Assert.AreEqual(0.0, mvt.Mahalanobis(new[] { 1.0, 2.0 }), 1E-12);

            // At (0, 0): z = (-1, -2), Q = z'Σ⁻¹z
            // Σ⁻¹ = (1/1.75)·((2,-0.5),(-0.5,1))
            // Q = (1/1.75)[(-1)²·2 + 2·(-1)·(-2)·(-0.5) + (-2)²·1]
            //   = (1/1.75)[2 - 2 + 4] = 4/1.75 ≈ 2.2857
            Assert.AreEqual(2.2857142857142856, mvt.Mahalanobis(new[] { 0.0, 0.0 }), 1E-10);

            // At (3, 4): z = (2, 2), Q = (1/1.75)[4·2 + 2·2·2·(-0.5) + 4·1]
            //   = (1/1.75)[8 - 4 + 4] = 8/1.75 ≈ 4.5714
            Assert.AreEqual(4.5714285714285712, mvt.Mahalanobis(new[] { 3.0, 4.0 }), 1E-10);
        }

        /// <summary>
        /// Verify that Mahalanobis throws for wrong dimension.
        /// </summary>
        [TestMethod]
        public void Test_Mahalanobis_WrongDimension()
        {
            var mvt = CreateStandard2D();
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                mvt.Mahalanobis(new[] { 1.0, 2.0, 3.0 }));
        }

        #endregion

        #region CDF Tests

        /// <summary>
        /// Verify 1D CDF matches univariate Student's t CDF.
        /// </summary>
        [TestMethod]
        public void Test_CDF_1D()
        {
            double nu = 5.0;
            double mu = 3.0;
            double sigma = 2.0;

            var mvt = new MultivariateStudentT(nu, new[] { mu }, new[,] { { sigma * sigma } });
            var univT = new StudentT(mu, sigma, nu);

            // CDF at x=5
            Assert.AreEqual(univT.CDF(5.0), mvt.CDF(new[] { 5.0 }), 1E-10);
            // CDF at x=3 (median)
            Assert.AreEqual(0.5, mvt.CDF(new[] { 3.0 }), 1E-10);
            // CDF at x=0
            Assert.AreEqual(univT.CDF(0.0), mvt.CDF(new[] { 0.0 }), 1E-10);
        }

        /// <summary>
        /// Verify 2D CDF against Monte Carlo reference values from scipy.
        /// These are approximate (MC tolerance ≈ 0.01).
        /// </summary>
        [TestMethod]
        public void Test_CDF_2D()
        {
            var mvt = CreateStandard2D();

            // CDF at (1, 2) ≈ 0.3076 (MC, N=10M)
            double cdf_12 = mvt.CDF(new[] { 1.0, 2.0 });
            Assert.AreEqual(0.3076, cdf_12, 0.01,
                $"CDF at (1,2) = {cdf_12:F4}, expected ≈ 0.3076");

            // CDF at (3, 5) ≈ 0.9166 (MC, N=10M)
            double cdf_35 = mvt.CDF(new[] { 3.0, 5.0 });
            Assert.AreEqual(0.9166, cdf_35, 0.01,
                $"CDF at (3,5) = {cdf_35:F4}, expected ≈ 0.9166");

            // CDF at (0, 0) ≈ 0.0454 (MC, N=10M)
            double cdf_00 = mvt.CDF(new[] { 0.0, 0.0 });
            Assert.AreEqual(0.0454, cdf_00, 0.01,
                $"CDF at (0,0) = {cdf_00:F4}, expected ≈ 0.0454");
        }

        #endregion

        #region Sampling Tests

        /// <summary>
        /// Verify that GenerateRandomValues produces samples with correct mean and covariance.
        /// For ν=5: E[X] = μ, Cov[X] = ν/(ν-2)·Σ = (5/3)·Σ.
        /// </summary>
        [TestMethod]
        public void Test_Sampling_MeanCovariance()
        {
            var mvt = CreateStandard2D();
            int N = 100000;
            var samples = mvt.GenerateRandomValues(N, 12345);

            // Compute sample mean
            double mean0 = 0, mean1 = 0;
            for (int i = 0; i < N; i++)
            {
                mean0 += samples[i, 0];
                mean1 += samples[i, 1];
            }
            mean0 /= N;
            mean1 /= N;

            Assert.AreEqual(1.0, mean0, 0.05, $"Sample mean[0] = {mean0:F4}");
            Assert.AreEqual(2.0, mean1, 0.05, $"Sample mean[1] = {mean1:F4}");

            // Compute sample covariance
            double cov00 = 0, cov01 = 0, cov11 = 0;
            for (int i = 0; i < N; i++)
            {
                double d0 = samples[i, 0] - mean0;
                double d1 = samples[i, 1] - mean1;
                cov00 += d0 * d0;
                cov01 += d0 * d1;
                cov11 += d1 * d1;
            }
            cov00 /= (N - 1);
            cov01 /= (N - 1);
            cov11 /= (N - 1);

            // Theoretical: ν/(ν-2) · Σ = (5/3) · [[1, 0.5], [0.5, 2]]
            Assert.AreEqual(5.0 / 3.0, cov00, 0.1, $"Sample Cov[0,0] = {cov00:F4}");
            Assert.AreEqual(5.0 / 6.0, cov01, 0.1, $"Sample Cov[0,1] = {cov01:F4}");
            Assert.AreEqual(10.0 / 3.0, cov11, 0.2, $"Sample Cov[1,1] = {cov11:F4}");
        }

        /// <summary>
        /// Verify that LatinHypercubeRandomValues produces samples with correct mean and covariance.
        /// </summary>
        [TestMethod]
        public void Test_LHS_MeanCovariance()
        {
            var mvt = CreateStandard2D();
            int N = 100000;
            var samples = mvt.LatinHypercubeRandomValues(N, 12345);

            // Compute sample mean
            double mean0 = 0, mean1 = 0;
            for (int i = 0; i < N; i++)
            {
                mean0 += samples[i, 0];
                mean1 += samples[i, 1];
            }
            mean0 /= N;
            mean1 /= N;

            Assert.AreEqual(1.0, mean0, 0.05, $"LHS mean[0] = {mean0:F4}");
            Assert.AreEqual(2.0, mean1, 0.05, $"LHS mean[1] = {mean1:F4}");

            // Compute sample covariance
            double cov00 = 0, cov01 = 0, cov11 = 0;
            for (int i = 0; i < N; i++)
            {
                double d0 = samples[i, 0] - mean0;
                double d1 = samples[i, 1] - mean1;
                cov00 += d0 * d0;
                cov01 += d0 * d1;
                cov11 += d1 * d1;
            }
            cov00 /= (N - 1);
            cov01 /= (N - 1);
            cov11 /= (N - 1);

            Assert.AreEqual(5.0 / 3.0, cov00, 0.1, $"LHS Cov[0,0] = {cov00:F4}");
            Assert.AreEqual(5.0 / 6.0, cov01, 0.1, $"LHS Cov[0,1] = {cov01:F4}");
            Assert.AreEqual(10.0 / 3.0, cov11, 0.2, $"LHS Cov[1,1] = {cov11:F4}");
        }

        /// <summary>
        /// Verify that StratifiedRandomValues produces valid samples without errors.
        /// </summary>
        [TestMethod]
        public void Test_StratifiedRandomValues()
        {
            var mvt = CreateStandard2D();

            // Create stratification bins from uniform quantiles
            var bins = new System.Collections.Generic.List<StratificationBin>();
            int nBins = 100;
            for (int i = 0; i < nBins; i++)
            {
                double lower = (double)i / nBins;
                double upper = (double)(i + 1) / nBins;
                bins.Add(new StratificationBin(lower, upper, 1));
            }

            var samples = mvt.StratifiedRandomValues(bins, 12345);

            // Basic checks: correct dimensions, finite values
            Assert.AreEqual(nBins, samples.GetLength(0));
            Assert.AreEqual(2, samples.GetLength(1));

            for (int i = 0; i < nBins; i++)
            {
                Assert.IsFalse(double.IsNaN(samples[i, 0]), $"NaN at sample [{i}, 0]");
                Assert.IsFalse(double.IsNaN(samples[i, 1]), $"NaN at sample [{i}, 1]");
                Assert.IsFalse(double.IsInfinity(samples[i, 0]), $"Inf at sample [{i}, 0]");
                Assert.IsFalse(double.IsInfinity(samples[i, 1]), $"Inf at sample [{i}, 1]");
            }
        }

        /// <summary>
        /// Verify that marginal distributions of samples follow univariate Student's t.
        /// Each marginal of MVT(ν, μ, Σ) is StudentT(μ_i, √Σ_ii, ν).
        /// Test via Kolmogorov-Smirnov style check on empirical CDF.
        /// </summary>
        [TestMethod]
        public void Test_Sampling_MarginalDistribution()
        {
            var mvt = CreateStandard2D();
            int N = 50000;
            var samples = mvt.GenerateRandomValues(N, 54321);

            // Extract first marginal and sort
            var x0 = new double[N];
            for (int i = 0; i < N; i++)
                x0[i] = samples[i, 0];
            Array.Sort(x0);

            // Theoretical marginal: StudentT(μ=1, σ=√1=1, ν=5)
            var marginal = new StudentT(1.0, 1.0, 5.0);

            // Compute max |empirical CDF - theoretical CDF|
            double maxDiff = 0;
            for (int i = 0; i < N; i++)
            {
                double empiricalCdf = (i + 1.0) / N;
                double theoreticalCdf = marginal.CDF(x0[i]);
                maxDiff = Math.Max(maxDiff, Math.Abs(empiricalCdf - theoreticalCdf));
            }

            // KS critical value for N=50000 at α=0.01 ≈ 1.63/√N ≈ 0.0073
            Assert.IsLessThan(0.01, maxDiff,
                $"KS statistic = {maxDiff:F5}, marginal does not match StudentT(1, 1, 5)");
        }

        #endregion

        #region Convergence Tests

        /// <summary>
        /// Verify that with very large ν, MVT PDF converges to MVN PDF.
        /// </summary>
        [TestMethod]
        public void Test_HighDF_ConvergesToNormal_PDF()
        {
            double nu = 100000.0;
            var location = new[] { 1.0, 2.0 };
            var sigma = new[,] { { 1.0, 0.5 }, { 0.5, 2.0 } };

            var mvt = new MultivariateStudentT(nu, location, sigma);
            var mvn = new MultivariateNormal(location, sigma);

            var testPoints = new[]
            {
                new[] { 1.0, 2.0 },
                new[] { 0.0, 0.0 },
                new[] { 3.0, 4.0 },
                new[] { -1.0, 5.0 }
            };

            foreach (var x in testPoints)
            {
                double pdfT = mvt.PDF(x);
                double pdfN = mvn.PDF(x);
                double relErr = Math.Abs(pdfT - pdfN) / pdfN;
                Assert.IsLessThan(5E-4, relErr,
                    $"PDF mismatch at ({x[0]},{x[1]}): MVT={pdfT:E6}, MVN={pdfN:E6}, relErr={relErr:E3}");
            }
        }

        /// <summary>
        /// Verify that with very large ν, MVT samples converge to MVN samples (same covariance).
        /// For large ν, Cov → Σ (since ν/(ν-2) → 1).
        /// </summary>
        [TestMethod]
        public void Test_HighDF_ConvergesToNormal_Sampling()
        {
            double nu = 100000.0;
            var location = new[] { 1.0, 2.0 };
            var sigma = new[,] { { 1.0, 0.5 }, { 0.5, 2.0 } };

            var mvt = new MultivariateStudentT(nu, location, sigma);
            int N = 50000;
            var samples = mvt.GenerateRandomValues(N, 12345);

            // Compute sample covariance
            double mean0 = 0, mean1 = 0;
            for (int i = 0; i < N; i++)
            {
                mean0 += samples[i, 0];
                mean1 += samples[i, 1];
            }
            mean0 /= N;
            mean1 /= N;

            double cov00 = 0, cov01 = 0, cov11 = 0;
            for (int i = 0; i < N; i++)
            {
                double d0 = samples[i, 0] - mean0;
                double d1 = samples[i, 1] - mean1;
                cov00 += d0 * d0;
                cov01 += d0 * d1;
                cov11 += d1 * d1;
            }
            cov00 /= (N - 1);
            cov01 /= (N - 1);
            cov11 /= (N - 1);

            // For large ν, Cov → Σ
            Assert.AreEqual(1.0, cov00, 0.05, $"Cov[0,0] = {cov00:F4}");
            Assert.AreEqual(0.5, cov01, 0.05, $"Cov[0,1] = {cov01:F4}");
            Assert.AreEqual(2.0, cov11, 0.1, $"Cov[1,1] = {cov11:F4}");
        }

        #endregion

        #region Validation Tests

        /// <summary>
        /// Verify that invalid parameters are properly rejected.
        /// </summary>
        [TestMethod]
        public void Test_ParameterValidation_NegativeDF()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                new MultivariateStudentT(-1.0, new[] { 0.0 }, new[,] { { 1.0 } }));
        }

        /// <summary>
        /// Verify that zero degrees of freedom is rejected.
        /// </summary>
        [TestMethod]
        public void Test_ParameterValidation_ZeroDF()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                new MultivariateStudentT(0.0, new[] { 0.0 }, new[,] { { 1.0 } }));
        }

        /// <summary>
        /// Verify that null location vector is rejected.
        /// </summary>
        [TestMethod]
        public void Test_ParameterValidation_NullLocation()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                new MultivariateStudentT(5.0, null, new[,] { { 1.0 } }));
        }

        /// <summary>
        /// Verify that null scale matrix is rejected.
        /// </summary>
        [TestMethod]
        public void Test_ParameterValidation_NullScaleMatrix()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                new MultivariateStudentT(5.0, new[] { 0.0 }, null));
        }

        /// <summary>
        /// Verify that non-square scale matrix is rejected.
        /// </summary>
        [TestMethod]
        public void Test_ParameterValidation_NonSquareMatrix()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                new MultivariateStudentT(5.0, new[] { 0.0, 0.0 }, new[,] { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 } }));
        }

        /// <summary>
        /// Verify that dimension mismatch between location and scale matrix is rejected.
        /// </summary>
        [TestMethod]
        public void Test_ParameterValidation_DimensionMismatch()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                new MultivariateStudentT(5.0, new[] { 0.0 }, new[,] { { 1.0, 0.0 }, { 0.0, 1.0 } }));
        }

        /// <summary>
        /// Verify that non-positive-definite scale matrix is rejected.
        /// </summary>
        [TestMethod]
        public void Test_ParameterValidation_NotPositiveDefinite()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                new MultivariateStudentT(5.0, new[] { 0.0, 0.0 }, new[,] { { 1.0, 2.0 }, { 2.0, 1.0 } }));
        }

        /// <summary>
        /// Verify that ValidateParameters returns exception instead of throwing when throwException=false.
        /// </summary>
        [TestMethod]
        public void Test_ValidateParameters_ReturnException()
        {
            var mvt = CreateStandard2D();
            var ex = mvt.ValidateParameters(-1.0, new[] { 0.0 }, new[,] { { 1.0 } }, false);
            Assert.IsNotNull(ex);

            ex = mvt.ValidateParameters(5.0, new[] { 0.0 }, new[,] { { 1.0 } }, false);
            Assert.IsNull(ex);
        }

        #endregion

        #region Clone Tests

        /// <summary>
        /// Verify that Clone creates an independent deep copy.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var original = CreateStandard2D();
            var clone = (MultivariateStudentT)original.Clone();

            // Verify properties match
            Assert.AreEqual(original.DegreesOfFreedom, clone.DegreesOfFreedom, 1E-10);
            Assert.AreEqual(original.Dimension, clone.Dimension);

            for (int i = 0; i < original.Dimension; i++)
            {
                Assert.AreEqual(original.Location[i], clone.Location[i], 1E-10);
            }

            // Verify PDF matches
            var x = new[] { 0.5, 3.0 };
            Assert.AreEqual(original.PDF(x), clone.PDF(x), 1E-14);

            // Verify independence: modify clone parameters
            clone.SetParameters(10.0, new[] { 0.0, 0.0 }, new[,] { { 2.0, 0.0 }, { 0.0, 2.0 } });
            Assert.AreEqual(5.0, original.DegreesOfFreedom, 1E-10);
            Assert.AreEqual(1.0, original.Location[0], 1E-10);
        }

        #endregion

        #region Constructor Tests

        /// <summary>
        /// Verify the dimension-only constructor creates standard MVT (zero mean, identity scale).
        /// </summary>
        [TestMethod]
        public void Test_Constructor_DimensionOnly()
        {
            var mvt = new MultivariateStudentT(3, 5.0);

            Assert.AreEqual(3, mvt.Dimension);
            Assert.AreEqual(5.0, mvt.DegreesOfFreedom, 1E-10);

            // Location should be zero
            for (int i = 0; i < 3; i++)
                Assert.AreEqual(0.0, mvt.Location[i], 1E-10);

            // Scale matrix should be identity
            var scale = mvt.ScaleMatrix;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(i == j ? 1.0 : 0.0, scale[i, j], 1E-10);
        }

        /// <summary>
        /// Verify the location-only constructor creates MVT with identity scale.
        /// </summary>
        [TestMethod]
        public void Test_Constructor_LocationOnly()
        {
            var mvt = new MultivariateStudentT(10.0, new[] { 1.0, 2.0, 3.0 });

            Assert.AreEqual(3, mvt.Dimension);
            Assert.AreEqual(10.0, mvt.DegreesOfFreedom, 1E-10);
            Assert.AreEqual(1.0, mvt.Location[0], 1E-10);
            Assert.AreEqual(2.0, mvt.Location[1], 1E-10);
            Assert.AreEqual(3.0, mvt.Location[2], 1E-10);
        }

        #endregion

        #region 3D Tests

        /// <summary>
        /// Verify PDF for a 3D MVT with identity scale matrix at the mode.
        /// </summary>
        [TestMethod]
        public void Test_PDF_3D_AtMode()
        {
            double nu = 10.0;
            var location = new[] { 0.0, 0.0, 0.0 };
            var mvt = new MultivariateStudentT(nu, location);

            // PDF at mode: C = Γ((10+3)/2) / [Γ(10/2) · (10π)^(3/2) · |I|^(1/2)]
            //             = Γ(6.5) / [Γ(5) · (10π)^1.5]
            // Γ(6.5) = 5.5·4.5·3.5·2.5·1.5·0.5·√π = 287.885 (approx)
            // Γ(5) = 4! = 24
            double expected = Math.Exp(
                Numerics.Mathematics.SpecialFunctions.Gamma.LogGamma(6.5) -
                Numerics.Mathematics.SpecialFunctions.Gamma.LogGamma(5.0) -
                1.5 * Math.Log(10.0 * Math.PI));

            double pdf = mvt.PDF(location);
            Assert.AreEqual(expected, pdf, 1E-10);
        }

        #endregion

        #region Non-Integer DF Tests

        /// <summary>
        /// Verify that non-integer degrees of freedom (e.g., ν=4.5) works correctly.
        /// This is important for Cohn's effective DF which can be non-integer.
        /// </summary>
        [TestMethod]
        public void Test_NonIntegerDF()
        {
            double nu = 4.5;
            var mvt = new MultivariateStudentT(nu, new[] { 0.0, 0.0 }, new[,] { { 1.0, 0.0 }, { 0.0, 1.0 } });

            // PDF at mode should be valid
            double pdfMode = mvt.PDF(new[] { 0.0, 0.0 });
            Assert.IsGreaterThan(0, pdfMode, "PDF at mode should be positive");
            Assert.IsFalse(double.IsNaN(pdfMode), "PDF at mode should not be NaN");

            // Sampling should work
            var samples = mvt.GenerateRandomValues(1000, 12345);
            Assert.AreEqual(1000, samples.GetLength(0));
            Assert.AreEqual(2, samples.GetLength(1));

            // All samples should be finite
            for (int i = 0; i < 1000; i++)
            {
                Assert.IsFalse(double.IsNaN(samples[i, 0]), $"NaN at sample [{i}, 0]");
                Assert.IsFalse(double.IsNaN(samples[i, 1]), $"NaN at sample [{i}, 1]");
                Assert.IsFalse(double.IsInfinity(samples[i, 0]), $"Inf at sample [{i}, 0]");
                Assert.IsFalse(double.IsInfinity(samples[i, 1]), $"Inf at sample [{i}, 1]");
            }

            // Variance should be ν/(ν-2) = 4.5/2.5 = 1.8
            Assert.AreEqual(1.8, mvt.Variance[0], 1E-10);
            Assert.AreEqual(1.8, mvt.Variance[1], 1E-10);
        }

        /// <summary>
        /// Verify that LHS sampling with non-integer DF produces correct statistics.
        /// </summary>
        [TestMethod]
        public void Test_NonIntegerDF_LHS_Statistics()
        {
            double nu = 7.5;
            var mvt = new MultivariateStudentT(nu, new[] { 2.0 }, new[,] { { 3.0 } });
            int N = 50000;
            var samples = mvt.LatinHypercubeRandomValues(N, 99999);

            double mean = 0;
            for (int i = 0; i < N; i++)
                mean += samples[i, 0];
            mean /= N;

            double var_ = 0;
            for (int i = 0; i < N; i++)
                var_ += (samples[i, 0] - mean) * (samples[i, 0] - mean);
            var_ /= (N - 1);

            // E[X] = 2.0
            Assert.AreEqual(2.0, mean, 0.05, $"Sample mean = {mean:F4}");

            // Var[X] = ν/(ν-2) · σ² = 7.5/5.5 · 3.0 = 4.0909
            double expectedVar = (nu / (nu - 2.0)) * 3.0;
            Assert.AreEqual(expectedVar, var_, 0.2, $"Sample var = {var_:F4}, expected = {expectedVar:F4}");
        }

        #endregion

        #region InverseCDF Tests

        /// <summary>
        /// Verify that InverseCDF produces a deterministic, reproducible sample point.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_Deterministic()
        {
            var mvt = CreateStandard2D();

            // Fixed probabilities: 2 for normal dims + 1 for χ²
            var probs = new[] { 0.5, 0.5, 0.5 };

            double[] x1 = mvt.InverseCDF(probs);
            double[] x2 = mvt.InverseCDF(probs);

            // Same input → same output
            Assert.AreEqual(x1[0], x2[0], 1E-14);
            Assert.AreEqual(x1[1], x2[1], 1E-14);
        }

        /// <summary>
        /// Verify that InverseCDF at median probabilities (all 0.5) returns the location vector.
        /// When all normal probabilities are 0.5, z = 0 for each dimension, so x = μ regardless of χ².
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_AtMedian()
        {
            var mvt = CreateStandard2D();

            // All normal probs = 0.5 → z = (0, 0). The χ² prob doesn't matter since L·0 = 0.
            var probs1 = new[] { 0.5, 0.5, 0.3 };
            var probs2 = new[] { 0.5, 0.5, 0.9 };

            double[] x1 = mvt.InverseCDF(probs1);
            double[] x2 = mvt.InverseCDF(probs2);

            // Both should return exactly μ = (1, 2)
            Assert.AreEqual(1.0, x1[0], 1E-10);
            Assert.AreEqual(2.0, x1[1], 1E-10);
            Assert.AreEqual(1.0, x2[0], 1E-10);
            Assert.AreEqual(2.0, x2[1], 1E-10);
        }

        /// <summary>
        /// Verify that InverseCDF throws for wrong-length probability array.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_WrongLength()
        {
            var mvt = CreateStandard2D();

            // Dimension = 2, so needs 3 probabilities
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                mvt.InverseCDF(new[] { 0.5, 0.5 }));  // too few

            Assert.Throws<ArgumentOutOfRangeException>(() =>
                mvt.InverseCDF(new[] { 0.5, 0.5, 0.5, 0.5 }));  // too many
        }

        /// <summary>
        /// Verify that small χ² probability produces more extreme samples (heavier tails).
        /// When the χ² probability is near 0, W is small, so √(ν/W) is large, amplifying z.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_TailBehavior()
        {
            var mvt = new MultivariateStudentT(5.0, new[] { 0.0 }, new[,] { { 1.0 } });

            // Fixed normal prob = 0.975 (z ≈ 1.96)
            double normalProb = 0.975;

            // Small χ² prob → large scaling → more extreme sample
            double[] xSmallChi = mvt.InverseCDF(new[] { normalProb, 0.01 });
            // Large χ² prob → small scaling → closer to normal
            double[] xLargeChi = mvt.InverseCDF(new[] { normalProb, 0.99 });

            Assert.IsGreaterThan(Math.Abs(xLargeChi[0]), Math.Abs(xSmallChi[0]),
                $"|x_small_chi| = {Math.Abs(xSmallChi[0]):F4} should be > |x_large_chi| = {Math.Abs(xLargeChi[0]):F4}");
        }

        /// <summary>
        /// Verify that InverseCDF-generated samples produce correct marginal statistics.
        /// Feed uniform random probabilities through InverseCDF and check mean and variance.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_MarginalStatistics()
        {
            double nu = 5.0;
            var mvt = new MultivariateStudentT(nu, new[] { 1.0, 2.0 }, new[,] { { 1.0, 0.5 }, { 0.5, 2.0 } });

            int N = 50000;
            var rnd = new MersenneTwister(42);
            double mean0 = 0, mean1 = 0;
            var samples0 = new double[N];
            var samples1 = new double[N];

            for (int i = 0; i < N; i++)
            {
                var probs = new[] { rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() };
                var x = mvt.InverseCDF(probs);
                samples0[i] = x[0];
                samples1[i] = x[1];
                mean0 += x[0];
                mean1 += x[1];
            }
            mean0 /= N;
            mean1 /= N;

            // E[X] = μ
            Assert.AreEqual(1.0, mean0, 0.05, $"InverseCDF mean[0] = {mean0:F4}");
            Assert.AreEqual(2.0, mean1, 0.05, $"InverseCDF mean[1] = {mean1:F4}");

            // Var[X] = ν/(ν-2) · diag(Σ)
            double var0 = 0, var1 = 0;
            for (int i = 0; i < N; i++)
            {
                var0 += (samples0[i] - mean0) * (samples0[i] - mean0);
                var1 += (samples1[i] - mean1) * (samples1[i] - mean1);
            }
            var0 /= (N - 1);
            var1 /= (N - 1);

            Assert.AreEqual(5.0 / 3.0, var0, 0.15, $"InverseCDF Var[0] = {var0:F4}");
            Assert.AreEqual(10.0 / 3.0, var1, 0.3, $"InverseCDF Var[1] = {var1:F4}");
        }

        /// <summary>
        /// Verify that InverseCDF matches GenerateRandomValues for the same uniform inputs.
        /// Both methods should produce identical results given the same underlying uniforms.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_ConsistentWithSampling()
        {
            var mvt = CreateStandard2D();

            // Generate a sample via GenerateRandomValues with known seed
            var rnd = new MersenneTwister(777);
            double u0 = rnd.NextDouble();
            double u1 = rnd.NextDouble();
            double u2 = rnd.NextDouble();

            // InverseCDF with the same uniforms
            double[] xInv = mvt.InverseCDF(new[] { u0, u1, u2 });

            // GenerateRandomValues uses the same sequence internally
            var samples = mvt.GenerateRandomValues(1, 777);

            Assert.AreEqual(samples[0, 0], xInv[0], 1E-10,
                $"InverseCDF[0] = {xInv[0]:F6}, GenerateRandomValues[0] = {samples[0, 0]:F6}");
            Assert.AreEqual(samples[0, 1], xInv[1], 1E-10,
                $"InverseCDF[1] = {xInv[1]:F6}, GenerateRandomValues[1] = {samples[0, 1]:F6}");
        }

        #endregion

    }
}