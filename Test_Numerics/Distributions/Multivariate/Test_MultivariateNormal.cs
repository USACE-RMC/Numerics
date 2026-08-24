using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;

namespace Distributions.Multivariate
{
    /// <summary>
    /// Unit tests for the Multivariate Normal distribution. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_MultivariateNormal
    {

        /// <summary>
        /// Asserts the action throws an ArgumentOutOfRangeException (net481-compatible).
        /// </summary>
        /// <param name="action">The action expected to throw.</param>
        private static void AssertThrowsOutOfRange(Action action)
        {
            try
            {
                action();
            }
            catch (ArgumentOutOfRangeException)
            {
                return;
            }
            Assert.Fail("Expected an ArgumentOutOfRangeException.");
        }


        /// <summary>
        /// Verifies the non-throwing mutable-covariance path: a positive-definite swap
        /// evaluates normally, a non-positive-definite swap marks the density invalid
        /// (LogPDF = -inf, PDF = 0) without throwing, and a subsequent valid swap
        /// restores evaluation.
        /// </summary>
        [TestMethod()]
        public void Test_TrySetCovariance_NonThrowingInvalidState()
        {
            var mvn = new MultivariateNormal(new[] { 1d, 2d }, new[,] { { 1d, 0.3d }, { 0.3d, 2d } });
            double validLogPdf = mvn.LogPDF(new[] { 1.2d, 1.8d });
            Assert.IsTrue(mvn.IsDensityValid);
            Assert.IsTrue(!double.IsNaN(validLogPdf) && !double.IsInfinity(validLogPdf));

            // A non-positive-definite covariance: correlation beyond one.
            Assert.IsFalse(mvn.TrySetCovariance(new[,] { { 1d, 1.5d }, { 1.5d, 1d } }));
            Assert.IsFalse(mvn.IsDensityValid);
            Assert.IsTrue(double.IsNegativeInfinity(mvn.LogPDF(new[] { 1.2d, 1.8d })));
            Assert.AreEqual(0d, mvn.PDF(new[] { 1.2d, 1.8d }), 0d);

            // Restore with the original covariance: the density comes back exactly.
            Assert.IsTrue(mvn.TrySetCovariance(new[,] { { 1d, 0.3d }, { 0.3d, 2d } }));
            Assert.IsTrue(mvn.IsDensityValid);
            Assert.AreEqual(validLogPdf, mvn.LogPDF(new[] { 1.2d, 1.8d }), 1E-12);

            // TrySetParameters also moves the mean.
            Assert.IsTrue(mvn.TrySetParameters(new[] { 0d, 0d }, new[,] { { 1d, 0d }, { 0d, 1d } }));
            Assert.AreEqual(-Math.Log(2d * Math.PI), mvn.LogPDF(new[] { 0d, 0d }), 1E-12);
        }

        /// <summary>
        /// Verifies the marginal utility: the sub-mean and sub-covariance at the kept
        /// indices, including a reordered subset, and the index validation matrix.
        /// </summary>
        [TestMethod()]
        public void Test_Marginal_SubsetAndValidation()
        {
            var mean = new[] { 1d, 2d, 3d };
            var covariance = new[,]
            {
                { 4.0d, 1.2d, 0.5d },
                { 1.2d, 9.0d, 2.1d },
                { 0.5d, 2.1d, 16.0d }
            };
            var mvn = new MultivariateNormal(mean, covariance);

            var marginal = mvn.Marginal(2, 0);
            Assert.AreEqual(2, marginal.Dimension);
            Assert.AreEqual(3d, marginal.Mean[0], 0d);
            Assert.AreEqual(1d, marginal.Mean[1], 0d);
            Assert.AreEqual(16d, marginal.Covariance[0, 0], 0d);
            Assert.AreEqual(0.5d, marginal.Covariance[0, 1], 0d);
            Assert.AreEqual(4d, marginal.Covariance[1, 1], 0d);

            // The marginal of a joint Gaussian is its sub-Gaussian: densities agree
            // with a directly constructed distribution.
            var direct = new MultivariateNormal(new[] { 3d, 1d },
                new[,] { { 16d, 0.5d }, { 0.5d, 4d } });
            Assert.AreEqual(direct.LogPDF(new[] { 2.5d, 1.4d }), marginal.LogPDF(new[] { 2.5d, 1.4d }), 1E-12);

            AssertThrowsOutOfRange(() => mvn.Marginal());
            AssertThrowsOutOfRange(() => mvn.Marginal(0, 0));
            AssertThrowsOutOfRange(() => mvn.Marginal(3));
        }

        /// <summary>
        /// Verifies the conditional utility against hand-computed Gaussian
        /// conditioning: for the bivariate case, X1 | X2 = x2 has mean
        /// mu1 + rho*(s1/s2)*(x2 - mu2) and variance s1^2*(1 - rho^2); a trivariate
        /// case checks the Schur complement, and observing everything throws.
        /// </summary>
        [TestMethod()]
        public void Test_Conditional_ClosedFormsAndValidation()
        {
            // Bivariate closed form.
            double mu1 = 4d, mu2 = 2d, s1 = 2d, s2 = 3d, rho = 0.6d;
            var bivariate = MultivariateNormal.Bivariate(mu1, mu2, s1, s2, rho);
            var conditional = bivariate.Conditional(new[] { 1 }, new[] { 3.5d });
            Assert.AreEqual(1, conditional.Dimension);
            Assert.AreEqual(mu1 + rho * (s1 / s2) * (3.5d - mu2), conditional.Mean[0], 1E-12);
            Assert.AreEqual(s1 * s1 * (1d - rho * rho), conditional.Covariance[0, 0], 1E-12);

            // Trivariate Schur complement, observing the middle dimension.
            var mean = new[] { 1d, 2d, 3d };
            var covariance = new[,]
            {
                { 4.0d, 1.2d, 0.5d },
                { 1.2d, 9.0d, 2.1d },
                { 0.5d, 2.1d, 16.0d }
            };
            var trivariate = new MultivariateNormal(mean, covariance);
            var given = trivariate.Conditional(new[] { 1 }, new[] { 4.0d });
            // mu_c = mu_{0,2} + Sigma_{[0,2],1} * (4 - 2) / 9
            Assert.AreEqual(1d + 1.2d * 2d / 9d, given.Mean[0], 1E-12);
            Assert.AreEqual(3d + 2.1d * 2d / 9d, given.Mean[1], 1E-12);
            // Sigma_c = Sigma_{[0,2],[0,2]} - outer(Sigma_{[0,2],1}) / 9
            Assert.AreEqual(4d - 1.2d * 1.2d / 9d, given.Covariance[0, 0], 1E-12);
            Assert.AreEqual(0.5d - 1.2d * 2.1d / 9d, given.Covariance[0, 1], 1E-12);
            Assert.AreEqual(16d - 2.1d * 2.1d / 9d, given.Covariance[1, 1], 1E-12);

            AssertThrowsOutOfRange(() => trivariate.Conditional(new[] { 0, 1, 2 }, new[] { 1d, 2d, 3d }));
            AssertThrowsOutOfRange(() => trivariate.Conditional(new[] { 1 }, new[] { 1d, 2d }));
        }


        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_MultivariateNormalDist()
        {

            var true_mean = new[] { 4d, 2d };
            var true_mode = new[] { 4d, 2d };
            var true_median = new[] { 4d, 2d };
            var true_variance = new[] { 0.3d, 0.7d };
            var true_stdDev = new[] { Math.Sqrt(0.3d), Math.Sqrt(0.7d) };
            double true_pdf1 = 0.000000018917884164743237d;
            double true_pdf2 = 0.35588127170858852d;
            double true_pdf3 = 0.000000000036520107734505265d;
            double true_cdf = 0.033944035782101548d;
            var MultiN = new MultivariateNormal(new[] { 4d, 2d }, new[,] { { 0.3d, 0.1d }, { 0.1d, 0.7d } });
            var mean = MultiN.Mean;
            var mode = MultiN.Mode;
            var median = MultiN.Median;
            var variance = MultiN.Variance;
            var stdev = MultiN.StandardDeviation;
            double pdf1 = MultiN.PDF(new[] { 2d, 5d });
            double pdf2 = MultiN.PDF(new[] { 4d, 2d });
            double pdf3 = MultiN.PDF(new[] { 3d, 7d });
            double cdf = MultiN.CDF(new[] { 3d, 5d });
            for (int i = 0; i <= MultiN.Dimension - 1; i++)
            {
                Assert.AreEqual(true_mean[i], mean[i], 0.0001d);
                Assert.AreEqual(true_median[i], median[i], 0.0001d);
                Assert.AreEqual(true_mode[i], mode[i], 0.0001d);
                Assert.AreEqual(true_stdDev[i], stdev[i], 0.0001d);
            }

            Assert.AreEqual(true_pdf1, pdf1, 0.0001d);
            Assert.AreEqual(true_pdf2, pdf2, 0.0001d);
            Assert.AreEqual(true_pdf3, pdf3, 0.0001d);
            Assert.AreEqual(true_cdf, cdf, 0.0001d);

        }

        /// <summary>
        /// Test that non-finite log-density evaluations return the canonical log(0) sentinel.
        /// </summary>
        [TestMethod]
        public void Test_LogPDF_NonFiniteInput_ReturnsNegativeInfinity()
        {
            var mvn = new MultivariateNormal(new[] { 0d, 0d }, new[,] { { 1d, 0d }, { 0d, 1d } });

            Assert.AreEqual(double.NegativeInfinity, mvn.LogPDF(new[] { double.PositiveInfinity, 0d }));
        }

        /// <summary>
        /// Verification against Fortran routine
        /// http://www.math.wsu.edu/faculty/genz/software/fort77/mvndstpack.f
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_Fortran()
        {
            var tester = new MultivariateNormal(5);
            tester.MVNUNI = new MersenneTwister(12345);
            int N = 5;
            int MAXPTS = 50;
            double ABSEPS = 0.00005;
            double RELEPS = 0;
            double[] Lower = new double[] { 0, 0, 1.7817, 1.4755, 1.5949 };
            double[] Upper = new double[] { 0, 1.5198, 1.7817, 1.4755, 1.5949 };
            int[] INFIN = new int[] { 1, 2, 1, 1, 0 };
            double[] CORREL = new double[] { -0.707107, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, .5, .5 };

            double[] expectedErrors = new double[] { 0.00000811074362075292, 0.00000480583149442263, 0.00000660161142196953 };
            double[] expectedValues = new double[] { 0.00286779150981026, 0.00111850297940743, 0.00397918930026649 };
            double[] expectedInforms = new double[] { 0, 0, 0 };

            for (int i = 0; i < 3; i++)
            {
                double ERROR = 0;
                double VAL = 0;
                int INFORM = 0;

                tester.MVNDST(N, Lower, Upper, INFIN, CORREL, MAXPTS, ABSEPS, RELEPS, ref ERROR, ref VAL, ref INFORM);

                Assert.IsTrue(Math.Abs(ERROR - expectedErrors[i]) < 1E-5 && Math.Abs(VAL - expectedValues[i]) < 1E-5 && Math.Abs(INFORM - expectedInforms[i]) == 0);
                INFIN[0] = INFIN[0] - 1;
            }

        }

        /// <summary>
        /// Verifies that two identically constructed instances left at the default generator produce
        /// bit-identical CDF values above two dimensions, and that the default equals explicit
        /// seeding with the published default seed. Three dimensions are required: two dimensions
        /// use a closed bivariate form that cannot detect a seeding regression.
        /// </summary>
        [TestMethod]
        public void Test_CDF_DefaultSeed_IsBitReproducibleAcrossInstances()
        {
            Assert.AreEqual(12345, MultivariateNormal.DefaultMVNUNISeed);

            var mean = new[] { 0d, 0d, 0d };
            var covariance = new[,]
            {
                { 1.0, 0.5, 0.3 },
                { 0.5, 1.0, 0.4 },
                { 0.3, 0.4, 1.0 }
            };
            var point = new[] { 0.5, 0.2, -0.3 };

            // MVNDST advances the generator on every evaluation, so each instance is evaluated
            // exactly once and the captured values are compared.
            double firstValue = new MultivariateNormal(mean, covariance).CDF(point);
            double secondValue = new MultivariateNormal(mean, covariance).CDF(point);
            Assert.AreEqual(
                BitConverter.DoubleToInt64Bits(firstValue),
                BitConverter.DoubleToInt64Bits(secondValue));

            var explicitlySeeded = new MultivariateNormal(mean, covariance)
            {
                MVNUNI = new MersenneTwister(MultivariateNormal.DefaultMVNUNISeed)
            };
            Assert.AreEqual(
                BitConverter.DoubleToInt64Bits(firstValue),
                BitConverter.DoubleToInt64Bits(explicitlySeeded.CDF(point)));
        }

        /// <summary>
        /// Verifies the documented invalid-dimension termination status without entering MVNDNT initialization.
        /// </summary>
        [TestMethod]
        public void Test_MVNDST_InvalidDimensionReturnsStatusTwo()
        {
            var distribution = new MultivariateNormal(1);
            double error = 0d;
            double value = 0d;
            int inform = 0;

            distribution.MVNDST(0, new double[0], new double[0], new int[0], new double[0], 1, 0d, 0d, ref error, ref value, ref inform);

            Assert.AreEqual(2, inform);
            Assert.AreEqual(0d, value, 0d);
            Assert.AreEqual(1d, error, 0d);
        }

        /// <summary>
        /// Verifies that dimensions flagged as completely unbounded collapse to probability one with successful status.
        /// </summary>
        [TestMethod]
        public void Test_MVNDST_AllUnboundedReturnsExactProbability()
        {
            var distribution = new MultivariateNormal(2);
            double error = -1d;
            double value = -1d;
            int inform = -1;

            distribution.MVNDST(
                2,
                new[] { 0d, 0d },
                new[] { 0d, 0d },
                new[] { -1, -1 },
                new[] { 0d },
                1,
                0d,
                0d,
                ref error,
                ref value,
                ref inform);

            Assert.AreEqual(0, inform);
            Assert.AreEqual(1d, value, 0d);
            Assert.AreEqual(0d, error, 0d);
        }

        /// <summary>
        /// Verifies the exact one-active-dimension Normal limit and its successful initialization status.
        /// </summary>
        [TestMethod]
        public void Test_MVNDST_OneActiveDimensionMatchesNormalProbability()
        {
            var distribution = new MultivariateNormal(2);
            double error = 0d;
            double value = 0d;
            int inform = -1;

            distribution.MVNDST(
                2,
                new[] { 0d, 0d },
                new[] { 0d, 0d },
                new[] { 0, -1 },
                new[] { 0.7d },
                1,
                0d,
                0d,
                ref error,
                ref value,
                ref inform);

            Assert.AreEqual(0, inform);
            Assert.AreEqual(0.5d, value, 1E-15);
            Assert.AreEqual(2E-16, error, 0d);
        }

        /// <summary>
        /// Verifies the analytical bivariate collapse against the zero-threshold quadrant probability for correlation one-half.
        /// </summary>
        [TestMethod]
        public void Test_MVNDST_BivariateCollapseMatchesAnalyticalProbability()
        {
            var distribution = new MultivariateNormal(2);
            double error = 0d;
            double value = 0d;
            int inform = -1;

            distribution.MVNDST(
                2,
                new[] { 0d, 0d },
                new[] { 0d, 0d },
                new[] { 0, 0 },
                new[] { 0.5d },
                1,
                0d,
                0d,
                ref error,
                ref value,
                ref inform);

            Assert.AreEqual(0, inform);
            Assert.AreEqual(1d / 3d, value, 1E-12);
            Assert.AreEqual(2E-16, error, 0d);
        }

        /// <summary>
        /// Verifies that only the lattice integration stage changes successful initialization into a budget-exhaustion status.
        /// </summary>
        [TestMethod]
        public void Test_MVNDST_InsufficientBudgetReturnsStatusOne()
        {
            var distribution = new MultivariateNormal(3) { MVNUNI = new MersenneTwister(12345) };
            double error = 0d;
            double value = 0d;
            int inform = -1;

            distribution.MVNDST(
                3,
                new[] { 0d, 0d, 0d },
                new[] { 0.2d, 0.5d, 1d },
                new[] { 0, 0, 0 },
                new[] { 0.25d, 0.1d, 0.3d },
                1,
                0d,
                0d,
                ref error,
                ref value,
                ref inform);

            Assert.AreEqual(1, inform);
            Assert.IsTrue(value >= 0d && value <= 1d);
            Assert.IsGreaterThanOrEqualTo(0d, error);
            Assert.IsFalse(double.IsNaN(error));
            Assert.IsFalse(double.IsInfinity(error));
        }

        /// <summary>
        /// Verified against R "mvtnorm" package
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_R()
        {

            double r = -0.33;
            var mean = new double[] { 0, 0, 0, 0 };
            var covar = new double[,]
                            {{ 1, r, r, r },
                             { r, 1, r, r },
                             { r, r, 1, r },
                             { r, r, r, 1 }};
            var mvn = new MultivariateNormal(mean, covar) { MVNUNI = new MersenneTwister(12345) };

            // AB
            var p = mvn.CDF(new[] { Normal.StandardZ(0.25), Normal.StandardZ(0.35), double.PositiveInfinity, double.PositiveInfinity });
            Assert.AreEqual(0.05011069, p, 1E-4);

            // AC
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), double.PositiveInfinity, Normal.StandardZ(0.5), double.PositiveInfinity });
            Assert.AreEqual(0.0827451, p, 1E-4);

            // AD
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), double.PositiveInfinity, double.PositiveInfinity, Normal.StandardZ(0.5) });
            Assert.AreEqual(0.0827451, p, 1E-4);

            // BC
            p = mvn.CDF(new[] { double.PositiveInfinity, Normal.StandardZ(0.35), Normal.StandardZ(0.5), double.PositiveInfinity });
            Assert.AreEqual(0.1254504, p, 1E-4);

            // BD
            p = mvn.CDF(new[] { double.PositiveInfinity, Normal.StandardZ(0.35), double.PositiveInfinity, Normal.StandardZ(0.5) });
            Assert.AreEqual(0.1254504, p, 1E-4);

            // CD
            p = mvn.CDF(new[] { double.PositiveInfinity, double.PositiveInfinity, Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(0.1964756, p, 1E-4);

            // ABC
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), Normal.StandardZ(0.35), Normal.StandardZ(0.5), double.PositiveInfinity });
            Assert.AreEqual(0.005960125, p, 1E-4);

            // ABD
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), Normal.StandardZ(0.35), double.PositiveInfinity, Normal.StandardZ(0.5) });
            Assert.AreEqual(0.005964513, p, 1E-4);

            // ACD
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), double.PositiveInfinity, Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(0.0128066, p, 1E-4);

            // BCD
            p = mvn.CDF(new[] { double.PositiveInfinity, Normal.StandardZ(0.35), Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(0.02324389, p, 1E-4);

            // ABCD
            p = mvn.CDF(new[] { Normal.StandardZ(0.25), Normal.StandardZ(0.35), Normal.StandardZ(0.5), Normal.StandardZ(0.5)});
            Assert.AreEqual(3.593582e-13, p, 1E-4);
        }

        /// <summary>
        /// Verified against R "mvtnorm" package.
        /// Test 1: Perfectly Negative correlation. 
        /// The matrix must be positive semi-definite. 
        /// So, the smallest allowable negative value is -1/(D-1) + an offset for machine double precision. 
        /// For simplicity, I offset by 0.01. 
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_R_PerfectNegative()
        {
            var mean = new double[] { 0, 0, 0 };
            var covar = new double[,]
            { { 1, -0.49,-0.49 },
              { -0.49, 1,-0.49 },
              { -0.49, -0.49, 1 }};
            var mvn = new MultivariateNormal(mean, covar) { MVNUNI = new MersenneTwister(12345) };

            var p = mvn.CDF(new[] { Normal.StandardZ(0.5), Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(0.002740932, p, 1E-4);
        }

        /// <summary>
        /// Verified against R "mvtnorm" package.
        /// Test 2: Perfectly positive correlation. 
        /// Again, I offset my 0.01 to keep it positive definite.
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_R_PerfectPositive()
        {

            var mean = new double[] { 0, 0, 0 };
            var covar = new double[,]
            { { 1, 0.99, 0.99 },
              { 0.99, 1, 0.99 },
              { 0.99, 0.99, 1 }};
            var mvn = new MultivariateNormal(mean, covar) { MVNUNI = new MersenneTwister(12345) };

            var p = mvn.CDF(new[] { Normal.StandardZ(0.5), Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(0.4661416, p, 1E-4);

        }

        /// <summary>
        /// Verified against R "mvtnorm" package.
        /// Test 3: Independent correlation. 
        /// </summary>
        [TestMethod]
        public void Test_MultivariateNormalCDF_R_Independent()
        {

            var mean = new double[] { 0, 0, 0 };
            var covar = new double[,]
            { { 1, 0, 0 },
              { 0, 1, 0 },
              { 0, 0, 1 }};
            var mvn = new MultivariateNormal(mean, covar) { MVNUNI = new MersenneTwister(12345) };

            var p = mvn.CDF(new[] { Normal.StandardZ(0.5), Normal.StandardZ(0.5), Normal.StandardZ(0.5) });
            Assert.AreEqual(0.125, p, 1E-4);

        }
         
        /// <summary>
        /// Verifies that perfectly correlated variables collapse to their tightest univariate bound.
        /// </summary>
        [TestMethod]
        public void Test_CDF_PerfectCorrelationCollapsesToNormal()
        {
            double probability = EvaluateStandardCdf(new[] { 0.8d, -0.3d, 0.2d }, new[] { 1d, 1d, 1d });
            Assert.AreEqual(Normal.StandardCDF(-0.3d), probability, 1E-10);
        }

        /// <summary>
        /// Verifies that perfectly anticorrelated variables collapse to a bounded normal interval.
        /// </summary>
        [TestMethod]
        public void Test_CDF_PerfectAnticorrelationCollapsesToNormalInterval()
        {
            double probability = EvaluateStandardCdf(new[] { 0.7d, 0.2d }, new[] { -1d });
            double expected = Normal.StandardCDF(0.7d) - Normal.StandardCDF(-0.2d);

            Assert.AreEqual(expected, probability, 1E-10);
        }

        /// <summary>
        /// Verifies equivalent results for permutations of a rank-deficient covariance matrix.
        /// </summary>
        [TestMethod]
        public void Test_CDF_PermutedRankDeficientMatricesCollapseAnalytically()
        {
            double expected = Normal.StandardCDF(0.2d) * Normal.StandardCDF(-0.4d);
            Assert.AreEqual(
                expected,
                EvaluateStandardCdf(new[] { 0.2d, -0.4d, 0.7d }, new[] { 0d, 1d, 0d }),
                1E-10);
            Assert.AreEqual(
                expected,
                EvaluateStandardCdf(new[] { 0.7d, 0.2d, -0.4d }, new[] { 1d, 0d, 0d }),
                1E-10);
        }

        /// <summary>
        /// Verifies the unchanged nonsingular path immediately above perfect correlation.
        /// </summary>
        [TestMethod]
        public void Test_CDF_NearSingularCorrelationMatchesBivariateFormula()
        {
            const double correlation = 1d - 1E-12;
            var distribution = MultivariateNormal.Bivariate(0d, 0d, 1d, 1d, correlation);
            distribution.MVNUNI = new MersenneTwister(12345);
            double expected = 0.25d + Math.Asin(correlation) / (2d * Math.PI);

            Assert.AreEqual(expected, distribution.CDF(new[] { 0d, 0d }), 1E-10);
        }

        /// <summary>
        /// Evaluates a standard multivariate-normal CDF through the public Genz integration entry point.
        /// </summary>
        /// <param name="upper">The upper integration limits.</param>
        /// <param name="correlations">The packed strict-lower-triangle correlation coefficients.</param>
        /// <returns>The evaluated multivariate-normal probability.</returns>
        private static double EvaluateStandardCdf(double[] upper, double[] correlations)
        {
            int dimensions = upper.Length;
            var distribution = new MultivariateNormal(dimensions)
            {
                MVNUNI = new MersenneTwister(12345)
            };
            var lower = new double[dimensions];
            var infinities = new int[dimensions];
            double error = 0d;
            double value = 0d;
            int inform = 0;

            distribution.MVNDST(
                dimensions,
                lower,
                upper,
                infinities,
                correlations,
                25000,
                1E-10,
                0d,
                ref error,
                ref value,
                ref inform);
            return value;
        }

        /// <summary>
        /// Test that the lattice-rule uniform generator cannot be assigned null.
        /// </summary>
        [TestMethod]
        public void Test_MVNUNI_RejectsNull()
        {
            var multivariate = new MultivariateNormal(2);
            Assert.Throws<ArgumentNullException>(() => multivariate.MVNUNI = null!);
        }

        #region Decomposition method selector (issue #145)

        /// <summary>
        /// Case 1 of the reference set: the mean vector paired with <see cref="Case1Covariance"/>.
        /// </summary>
        private static readonly double[] Case1Mean = { 1d, 2d, 3d };

        /// <summary>
        /// Case 1 of the reference set: a non-singular 3-D covariance.
        /// </summary>
        private static readonly double[,] Case1Covariance =
        {
            { 4d, 1d, 0.5d },
            { 1d, 3d, 0.25d },
            { 0.5d, 0.25d, 2d }
        };

        /// <summary>
        /// Case 2 of the reference set: a singular 3-D covariance of rank two, where the third component
        /// is exactly the first. This is the gridded-data structure described in issue #145.
        /// </summary>
        private static readonly double[,] Case2Covariance =
        {
            { 2d, 0.5d, 2d },
            { 0.5d, 1d, 0.5d },
            { 2d, 0.5d, 2d }
        };

        /// <summary>
        /// Case 3 of the reference set: a singular 2-D covariance of rank one — two perfectly correlated
        /// components.
        /// </summary>
        private static readonly double[,] Case3Covariance =
        {
            { 1d, 1d },
            { 1d, 1d }
        };

        /// <summary>
        /// The tolerance applied to log densities compared against the scipy oracle. The measured
        /// disagreement across every reference point is at most 3.6E-15, so this is that measurement
        /// rounded out.
        /// </summary>
        private const double OracleTolerance = 1E-14;

        /// <summary>
        /// Asserts the action throws, and returns the exception, without requiring a specific type
        /// (net481-compatible).
        /// </summary>
        /// <param name="action">The action expected to throw.</param>
        /// <returns>The exception that was thrown.</returns>
        private static Exception AssertThrowsAny(Action action)
        {
            try
            {
                action();
            }
            catch (Exception ex)
            {
                return ex;
            }
            Assert.Fail("Expected an exception.");
            throw new InvalidOperationException("Unreachable.");
        }

        /// <summary>
        /// Builds the singular value decomposition sampling factor A = U*sqrt(W) for a covariance matrix,
        /// which is the factor the distribution uses under
        /// <see cref="DecompositionMethod.SingularValue"/>.
        /// </summary>
        /// <param name="covariance">The covariance matrix.</param>
        /// <returns>The factor A, which satisfies A*Aᵀ = covariance.</returns>
        private static double[,] SingularValueFactor(double[,] covariance)
        {
            int n = covariance.GetLength(0);
            var svd = new SingularValueDecomposition(new Matrix(covariance));
            double threshold = svd.Threshold;
            var factor = new double[n, n];
            for (int j = 0; j < n; j++)
            {
                double scale = svd.W[j] > threshold ? Math.Sqrt(svd.W[j]) : 0d;
                for (int i = 0; i < n; i++)
                    factor[i, j] = svd.U[i, j] * scale;
            }
            return factor;
        }

        /// <summary>
        /// Verifies that the decomposition selector defaults to Cholesky on every pre-existing constructor
        /// and is carried by the new overloads.
        /// </summary>
        [TestMethod]
        public void Test_Decomposition_DefaultsToCholesky()
        {
            Assert.AreEqual(DecompositionMethod.Cholesky, new MultivariateNormal(2).Decomposition);
            Assert.AreEqual(DecompositionMethod.Cholesky, new MultivariateNormal(new[] { 0d, 0d }).Decomposition);
            Assert.AreEqual(DecompositionMethod.Cholesky, new MultivariateNormal(Case1Mean, Case1Covariance).Decomposition);

            Assert.AreEqual(DecompositionMethod.SingularValue, new MultivariateNormal(2, DecompositionMethod.SingularValue).Decomposition);
            Assert.AreEqual(DecompositionMethod.SingularValue, new MultivariateNormal(new[] { 0d, 0d }, DecompositionMethod.SingularValue).Decomposition);
            Assert.AreEqual(DecompositionMethod.SingularValue, new MultivariateNormal(Case1Mean, Case1Covariance, DecompositionMethod.SingularValue).Decomposition);
        }

        /// <summary>
        /// Case 1: a non-singular covariance, where the Cholesky and singular value paths are
        /// mathematically identical. Both must reproduce the scipy oracle and must agree with each other
        /// to near machine precision.
        /// </summary>
        /// <remarks>
        /// Oracle: <c>scipy.stats.multivariate_normal(mean, cov).logpdf(x)</c>, scipy 1.17.1. This case is
        /// the backward-compatibility anchor for issue #145: the singular value path must not change any
        /// answer that the Cholesky path already gets right.
        /// </remarks>
        [TestMethod]
        public void Test_SVD_Case1_NonSingularMatchesScipyAndCholesky()
        {
            var points = new[]
            {
                new[] { 1d, 2d, 3d },
                new[] { 0d, 0d, 0d },
                new[] { 2.5d, 1d, 4d },
                new[] { -1d, 5d, 2d }
            };
            var oracle = new[]
            {
                -4.2849940472992305d,
                -6.989405812005113d,
                -5.108155812005113d,
                -7.226170517887469d
            };

            var cholesky = new MultivariateNormal(Case1Mean, Case1Covariance);
            var singular = new MultivariateNormal(Case1Mean, Case1Covariance, DecompositionMethod.SingularValue);

            Assert.IsTrue(cholesky.IsPositiveDefinite);
            Assert.IsTrue(singular.IsPositiveDefinite);

            for (int i = 0; i < points.Length; i++)
            {
                Assert.AreEqual(oracle[i], cholesky.LogPDF(points[i]), OracleTolerance);
                Assert.AreEqual(oracle[i], singular.LogPDF(points[i]), OracleTolerance);
                Assert.AreEqual(cholesky.LogPDF(points[i]), singular.LogPDF(points[i]), OracleTolerance);
                Assert.AreEqual(Math.Exp(oracle[i]), singular.PDF(points[i]), 1E-15d);
            }
        }

        /// <summary>
        /// Case 2: a singular covariance of rank two — the gridded-data structure of issue #145, where the
        /// third component repeats the first. The singular value path must reproduce the scipy oracle,
        /// report the covariance as not positive-definite, and return zero density off the support.
        /// </summary>
        /// <remarks>
        /// Oracle: <c>scipy.stats.multivariate_normal(mean, cov, allow_singular=True).logpdf(x)</c>,
        /// scipy 1.17.1. The log pseudo-determinant is 1.2527629684953678.
        /// </remarks>
        [TestMethod]
        public void Test_SVD_Case2_SingularRankTwoMatchesScipy()
        {
            var mean = new[] { 0d, 0d, 0d };
            var singular = new MultivariateNormal(mean, Case2Covariance, DecompositionMethod.SingularValue);

            Assert.IsFalse(singular.IsPositiveDefinite);
            Assert.AreEqual(-2.4642585506570285d, singular.LogPDF(new[] { 0d, 0d, 0d }), OracleTolerance);
            Assert.AreEqual(-2.7499728363713145d, singular.LogPDF(new[] { 1d, 0.5d, 1d }), OracleTolerance);
            Assert.AreEqual(-3.607115693514173d, singular.LogPDF(new[] { -1d, 1d, -1d }), OracleTolerance);
            Assert.AreEqual(Math.Exp(-2.4642585506570285d), singular.PDF(new[] { 0d, 0d, 0d }), 1E-15d);

            // The support is the plane x3 == x1; a point off it has density exactly zero.
            Assert.AreEqual(double.NegativeInfinity, singular.LogPDF(new[] { 1d, 0.5d, 1.5d }));
            Assert.AreEqual(0d, singular.PDF(new[] { 1d, 0.5d, 1.5d }), 0d);
        }

        /// <summary>
        /// Case 3: a singular covariance of rank one — two perfectly correlated components. The singular
        /// value path must reproduce the scipy oracle on the support and return negative infinity off it.
        /// </summary>
        /// <remarks>
        /// Oracle: <c>scipy.stats.multivariate_normal([0, 0], [[1, 1], [1, 1]], allow_singular=True)</c>,
        /// scipy 1.17.1. The log pseudo-determinant is 0.6931471805599453 = log(2).
        /// </remarks>
        [TestMethod]
        public void Test_SVD_Case3_SingularRankOneMatchesScipy()
        {
            var mean = new[] { 0d, 0d };
            var singular = new MultivariateNormal(mean, Case3Covariance, DecompositionMethod.SingularValue);

            Assert.IsFalse(singular.IsPositiveDefinite);
            Assert.AreEqual(-1.2655121234846454d, singular.LogPDF(new[] { 0d, 0d }), OracleTolerance);
            Assert.AreEqual(-1.7655121234846454d, singular.LogPDF(new[] { 1d, 1d }), OracleTolerance);
            Assert.AreEqual(-1.3905121234846454d, singular.LogPDF(new[] { -0.5d, -0.5d }), OracleTolerance);
            Assert.AreEqual(Math.Exp(-1.2655121234846454d), singular.PDF(new[] { 0d, 0d }), 1E-15d);

            foreach (var offSupport in new[] { new[] { 1d, -1d }, new[] { 0d, 1d }, new[] { 2d, 1d } })
            {
                Assert.AreEqual(double.NegativeInfinity, singular.LogPDF(offSupport));
                Assert.AreEqual(0d, singular.PDF(offSupport), 0d);
            }
        }

        /// <summary>
        /// Guards the single most important line of the singular value path: the normalizing constant must
        /// use the rank of the covariance matrix, not its dimension.
        /// </summary>
        /// <remarks>
        /// An SVD fallback carried by this class in the repository's first commit (46d5487, removed in
        /// ad0b293) computed the constant as
        /// <c>-(Math.Log(2d * Math.PI) * _mean.Length + lndet) * 0.5d</c> — the dimension where the rank
        /// belongs. Using the dimension shifts the log density by exactly 0.5 * (d - rank) * log(2π), which
        /// is 0.9189385332046727 per deficient dimension, a factor of 2.5066282746310002 on the density.
        /// Cases 2 and 3 are each rank-deficient by exactly one, so that bug would put every one of their
        /// oracle values off by that fixed amount. This test asserts both that the oracle value is
        /// reproduced and that the value the old constant would have produced is exactly that far away, so
        /// it fails loudly if the constant ever reverts.
        /// </remarks>
        [TestMethod]
        public void Test_SVD_NormalizingConstantUsesRankNotDimension()
        {
            double shiftPerDeficientDimension = 0.5d * Math.Log(2d * Math.PI);
            Assert.AreEqual(0.9189385332046727d, shiftPerDeficientDimension, 1E-15d);

            // Case 2: dimension 3, rank 2 — deficient by one.
            var case2 = new MultivariateNormal(new[] { 0d, 0d, 0d }, Case2Covariance, DecompositionMethod.SingularValue);
            double case2Oracle = -2.4642585506570285d;
            double case2Actual = case2.LogPDF(new[] { 0d, 0d, 0d });
            Assert.AreEqual(case2Oracle, case2Actual, OracleTolerance);
            Assert.IsGreaterThan(0.9d, Math.Abs(case2Actual - (case2Oracle - shiftPerDeficientDimension)),
                "The log density must not carry the dimension-based constant.");

            // Case 3: dimension 2, rank 1 — deficient by one.
            var case3 = new MultivariateNormal(new[] { 0d, 0d }, Case3Covariance, DecompositionMethod.SingularValue);
            double case3Oracle = -1.2655121234846454d;
            double case3Actual = case3.LogPDF(new[] { 0d, 0d });
            Assert.AreEqual(case3Oracle, case3Actual, OracleTolerance);
            Assert.IsGreaterThan(0.9d, Math.Abs(case3Actual - (case3Oracle - shiftPerDeficientDimension)),
                "The log density must not carry the dimension-based constant.");

            // Case 1 is full rank, so rank and dimension agree and the constant is unchanged there.
            var case1 = new MultivariateNormal(Case1Mean, Case1Covariance, DecompositionMethod.SingularValue);
            Assert.AreEqual(-4.2849940472992305d, case1.LogPDF(new[] { 1d, 2d, 3d }), OracleTolerance);
        }

        /// <summary>
        /// Verifies the sampling factor built by the singular value path reproduces the covariance matrix,
        /// A*Aᵀ = Σ, for all three reference cases, and that the distribution actually samples with it.
        /// </summary>
        /// <remarks>
        /// A = U*sqrt(W) is the factor NumPy builds for <c>multivariate_normal(..., method='svd')</c>.
        /// Because Σ is symmetric positive semi-definite the left and right singular vectors coincide for
        /// every singular value above the threshold, so U*W*Uᵀ = Σ. The measured reconstruction error is at
        /// most 2.3E-15 across the three cases. The second half of the test recovers the factor the
        /// distribution itself holds by driving <see cref="MultivariateNormal.InverseCDF"/> with a known
        /// vector of standard normal variates.
        /// </remarks>
        [TestMethod]
        public void Test_SVD_SamplingFactorReproducesCovariance()
        {
            var covariances = new[] { Case1Covariance, Case2Covariance, Case3Covariance };
            foreach (var covariance in covariances)
            {
                int n = covariance.GetLength(0);
                var factor = SingularValueFactor(covariance);
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        double sum = 0d;
                        for (int k = 0; k < n; k++)
                            sum += factor[i, k] * factor[j, k];
                        Assert.AreEqual(covariance[i, j], sum, 1E-14d);
                    }
                }
            }

            // The distribution's own factor: InverseCDF maps probabilities to x = A*z + mu.
            var probabilities = new[] { 0.15d, 0.62d, 0.93d };
            var z = new double[3];
            for (int j = 0; j < 3; j++)
                z[j] = Normal.StandardZ(probabilities[j]);

            var expectedFactor = SingularValueFactor(Case1Covariance);
            var singular = new MultivariateNormal(Case1Mean, Case1Covariance, DecompositionMethod.SingularValue);
            var actual = singular.InverseCDF(probabilities);
            for (int i = 0; i < 3; i++)
            {
                double expected = Case1Mean[i];
                for (int j = 0; j < 3; j++)
                    expected += expectedFactor[i, j] * z[j];
                Assert.AreEqual(expected, actual[i], 1E-12d);
            }
        }

        /// <summary>
        /// Verifies that seeded draws from a rank-one covariance under the singular value path satisfy the
        /// collinearity exactly, which is the property the reporter of issue #145 currently obtains by
        /// performing the decomposition by hand.
        /// </summary>
        /// <remarks>
        /// The null direction of Σ receives a zero column in A = U*sqrt(W), so it carries no noise and every
        /// draw lands on the support x1 = x2. The measured departure across the seeded sample is at most
        /// 8.9E-16, which is roundoff in the matrix-vector product rather than sampling noise.
        /// </remarks>
        [TestMethod]
        public void Test_SVD_SeededDrawsOnRankOneCovarianceAreCollinear()
        {
            var singular = new MultivariateNormal(new[] { 0d, 0d }, Case3Covariance, DecompositionMethod.SingularValue);

            var sample = singular.GenerateRandomValues(200, 4321);
            for (int i = 0; i < 200; i++)
            {
                Assert.AreEqual(sample[i, 0], sample[i, 1], 1E-14d);
                // Every draw is on the support, so the density there is finite.
                Assert.IsGreaterThan(double.NegativeInfinity, singular.LogPDF(new[] { sample[i, 0], sample[i, 1] }));
            }

            var latin = singular.LatinHypercubeRandomValues(50, 777);
            for (int i = 0; i < 50; i++)
                Assert.AreEqual(latin[i, 0], latin[i, 1], 1E-14d);

            var inverse = singular.InverseCDF(new[] { 0.2d, 0.8d });
            Assert.AreEqual(inverse[0], inverse[1], 1E-14d);
        }

        /// <summary>
        /// Pins the seeded output of the default Cholesky path so that any change to the decomposition
        /// plumbing that perturbed it would fail here.
        /// </summary>
        /// <remarks>
        /// The expected values were captured from the library immediately before the decomposition selector
        /// of issue #145 was added, and are asserted exactly — not within a tolerance — so the Cholesky path
        /// is held bit-identical.
        /// </remarks>
        [TestMethod]
        public void Test_Cholesky_SeededOutputIsUnchanged()
        {
            var cholesky = new MultivariateNormal(Case1Mean, Case1Covariance);

            var expected = new[,]
            {
                { 3.9458754639257694d, 4.771800776088886d, 2.7965750395390128d },
                { -1.246107872880183d, -0.05488903400275058d, 0.21419166312288151d },
                { -0.6508844921537047d, 3.146382948321917d, 3.101603919764241d },
                { 1.1609885052730824d, 2.441279589126226d, 5.414221213964924d },
                { 4.611326489815959d, 2.74432691956876d, 3.9917091805407323d }
            };
            var sample = cholesky.GenerateRandomValues(5, 12345);
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(expected[i, j], sample[i, j], 0d);
            }

            var latin = cholesky.LatinHypercubeRandomValues(3, 987);
            var expectedLatin = new[,]
            {
                { 2.742364437119263d, 1.795955104586471d, 2.0437569460738882d },
                { -3.2097827841193878d, 2.29997314661673d, 2.4025926343175597d },
                { 1.359253584263394d, 0.1502253290672808d, 4.165190525051644d }
            };
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(expectedLatin[i, j], latin[i, j], 0d);
            }

            var inverse = cholesky.InverseCDF(new[] { 0.1d, 0.5d, 0.9d });
            Assert.AreEqual(-1.5631031310892016d, inverse[0], 0d);
            Assert.AreEqual(1.3592242172276996d, inverse[1], 0d);
            Assert.AreEqual(4.460838864683783d, inverse[2], 0d);

            Assert.AreEqual(-4.2849940472992305d, cholesky.LogPDF(new[] { 1d, 2d, 3d }), 0d);
            Assert.AreEqual(-6.9894058120051135d, cholesky.LogPDF(new[] { 0d, 0d, 0d }), 0d);
            Assert.AreEqual(-5.108155812005113d, cholesky.LogPDF(new[] { 2.5d, 1d, 4d }), 0d);
            Assert.AreEqual(-7.226170517887466d, cholesky.LogPDF(new[] { -1d, 5d, 2d }), 0d);
        }

        /// <summary>
        /// Verifies the validation contract of the singular value path: a positive semi-definite covariance
        /// is accepted, while an indefinite, an asymmetric, or a non-finite one is still rejected.
        /// </summary>
        /// <remarks>
        /// Singular values are unsigned, so a negative eigenvalue is detected by comparing the matrix
        /// against its own symmetric reconstruction U*W*Uᵀ rather than by inspecting the spectrum. The
        /// matrix <c>{{0, 2}, {2, 0}}</c> is included because its eigenvalues are +2 and -2: a sign test on
        /// the singular vectors is ambiguous there, while the reconstruction test is not.
        /// </remarks>
        [TestMethod]
        public void Test_SVD_ValidationAcceptsSemiDefiniteAndRejectsIndefinite()
        {
            // Accepted: singular but positive semi-definite.
            var accepted = new MultivariateNormal(new[] { 0d, 0d }, Case3Covariance, DecompositionMethod.SingularValue);
            Assert.IsNull(accepted.ValidateParameters(new[] { 0d, 0d }, Case3Covariance, false));

            // Rejected: a negative eigenvalue.
            AssertThrowsOutOfRange(() => new MultivariateNormal(new[] { 0d, 0d }, new[,] { { 1d, 2d }, { 2d, 1d } }, DecompositionMethod.SingularValue));
            AssertThrowsOutOfRange(() => new MultivariateNormal(new[] { 0d, 0d }, new[,] { { -1d, 0d }, { 0d, 2d } }, DecompositionMethod.SingularValue));
            AssertThrowsOutOfRange(() => new MultivariateNormal(new[] { 0d, 0d }, new[,] { { 0d, 2d }, { 2d, 0d } }, DecompositionMethod.SingularValue));

            // Rejected: asymmetric.
            AssertThrowsOutOfRange(() => new MultivariateNormal(new[] { 0d, 0d }, new[,] { { 1d, 0.5d }, { 0.2d, 1d } }, DecompositionMethod.SingularValue));

            // Rejected: non-finite entries.
            AssertThrowsOutOfRange(() => new MultivariateNormal(new[] { 0d, 0d }, new[,] { { 1d, double.NaN }, { double.NaN, 1d } }, DecompositionMethod.SingularValue));
            AssertThrowsOutOfRange(() => new MultivariateNormal(new[] { 0d, 0d }, new[,] { { double.PositiveInfinity, 0d }, { 0d, 1d } }, DecompositionMethod.SingularValue));

            // Accepted: a well-conditioned covariance at a large scale, where the reconstruction residual
            // grows with the magnitude of the entries and the tolerance must scale with it.
            var large = new MultivariateNormal(new[] { 0d, 0d }, new[,] { { 4e8d, 1e8d }, { 1e8d, 3e8d } }, DecompositionMethod.SingularValue);
            Assert.IsTrue(large.IsPositiveDefinite);

            // The non-throwing path reports the same decision without raising.
            var mutable = new MultivariateNormal(new[] { 0d, 0d }, new[,] { { 1d, 0d }, { 0d, 1d } }, DecompositionMethod.SingularValue);
            Assert.IsTrue(mutable.TrySetCovariance(Case3Covariance));
            Assert.IsTrue(mutable.IsDensityValid);
            Assert.AreEqual(-1.7655121234846454d, mutable.LogPDF(new[] { 1d, 1d }), OracleTolerance);
            Assert.IsFalse(mutable.TrySetCovariance(new[,] { { 1d, 2d }, { 2d, 1d } }));
            Assert.IsFalse(mutable.IsDensityValid);
        }

        /// <summary>
        /// Documents what the Cholesky path does with the two singular reference covariances, which is the
        /// behaviour issue #145 was filed about.
        /// </summary>
        /// <remarks>
        /// Case 3, the rank-one covariance, fails cleanly: the Cholesky factorization reaches a
        /// non-positive pivot and throws. Case 2, the rank-two covariance, does <b>not</b> throw — its final
        /// pivot evaluates to about 4.4E-16 rather than exactly zero, so the factorization completes with a
        /// pivot of order 1E-8 and the resulting log density is wrong by about 17 nats. This test
        /// pins both behaviours as they stand today and shows the singular value path getting the right
        /// answer where the Cholesky path does not. Neither Cholesky behaviour is changed here.
        /// </remarks>
        [TestMethod]
        public void Test_Cholesky_BehaviourOnSingularCovariancesIsUnchanged()
        {
            // Case 3: a clean failure.
            var exception = AssertThrowsAny(() => new MultivariateNormal(new[] { 0d, 0d }, Case3Covariance));
            Assert.IsGreaterThanOrEqualTo(0, exception.Message.IndexOf("positive-definite", StringComparison.OrdinalIgnoreCase));

            // Case 2: no failure, but an unusable density — the silent failure mode the issue describes.
            var cholesky = new MultivariateNormal(new[] { 0d, 0d, 0d }, Case2Covariance);
            double choleskyLogPdf = cholesky.LogPDF(new[] { 0d, 0d, 0d });
            double oracle = -2.4642585506570285d;
            Assert.IsGreaterThan(1d, Math.Abs(choleskyLogPdf - oracle),
                "The Cholesky path is expected to be far from the correct density on a rank-deficient covariance.");

            var singular = new MultivariateNormal(new[] { 0d, 0d, 0d }, Case2Covariance, DecompositionMethod.SingularValue);
            Assert.AreEqual(oracle, singular.LogPDF(new[] { 0d, 0d, 0d }), OracleTolerance);
        }

        /// <summary>
        /// Verifies that <see cref="MultivariateNormal.Clone"/> carries the decomposition selector and the
        /// factorization that goes with it.
        /// </summary>
        [TestMethod]
        public void Test_Decomposition_ClonePreservesTheSelector()
        {
            var singular = new MultivariateNormal(new[] { 0d, 0d }, Case3Covariance, DecompositionMethod.SingularValue);
            var singularClone = (MultivariateNormal)singular.Clone();
            Assert.AreEqual(DecompositionMethod.SingularValue, singularClone.Decomposition);
            Assert.IsFalse(singularClone.IsPositiveDefinite);
            Assert.AreEqual(-1.7655121234846454d, singularClone.LogPDF(new[] { 1d, 1d }), OracleTolerance);
            Assert.AreEqual(double.NegativeInfinity, singularClone.LogPDF(new[] { 1d, -1d }));
            var clonedSample = singularClone.GenerateRandomValues(10, 4321);
            for (int i = 0; i < 10; i++)
                Assert.AreEqual(clonedSample[i, 0], clonedSample[i, 1], 1E-14d);

            var cholesky = new MultivariateNormal(Case1Mean, Case1Covariance);
            var choleskyClone = (MultivariateNormal)cholesky.Clone();
            Assert.AreEqual(DecompositionMethod.Cholesky, choleskyClone.Decomposition);
            Assert.IsTrue(choleskyClone.IsPositiveDefinite);
            var original = cholesky.GenerateRandomValues(5, 12345);
            var copied = choleskyClone.GenerateRandomValues(5, 12345);
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(original[i, j], copied[i, j], 0d);
            }
        }

        /// <summary>
        /// Verifies that the paths the selector deliberately does not govern behave as they do today: the
        /// Genz MVNDST integrator behind <see cref="MultivariateNormal.CDF"/> factorizes the correlation
        /// matrix internally, and <see cref="MultivariateNormal.Conditional"/> and
        /// <see cref="MultivariateNormal.Marginal"/> factorize the sub-covariance with their own Cholesky
        /// decomposition and return Cholesky-based distributions.
        /// </summary>
        [TestMethod]
        public void Test_Decomposition_CdfAndConditionalHelpersAreUnaffected()
        {
            var cholesky = new MultivariateNormal(Case1Mean, Case1Covariance);
            var singular = new MultivariateNormal(Case1Mean, Case1Covariance, DecompositionMethod.SingularValue);

            // MVNDST is a randomized lattice rule that advances its generator, so compare first
            // evaluations on freshly constructed instances.
            Assert.AreEqual(cholesky.CDF(new[] { 2d, 3d, 4d }), singular.CDF(new[] { 2d, 3d, 4d }), 0d);

            // The conditional and marginal helpers return distributions on the default selector.
            var conditional = singular.Conditional(new[] { 2 }, new[] { 3d });
            Assert.AreEqual(DecompositionMethod.Cholesky, conditional.Decomposition);
            var marginal = singular.Marginal(0, 1);
            Assert.AreEqual(DecompositionMethod.Cholesky, marginal.Decomposition);

            // They agree with the Cholesky parent, which is the behaviour they have today.
            var choleskyConditional = cholesky.Conditional(new[] { 2 }, new[] { 3d });
            Assert.AreEqual(choleskyConditional.LogPDF(new[] { 1d, 2d }), conditional.LogPDF(new[] { 1d, 2d }), 0d);
            var choleskyMarginal = cholesky.Marginal(0, 1);
            Assert.AreEqual(choleskyMarginal.LogPDF(new[] { 1d, 2d }), marginal.LogPDF(new[] { 1d, 2d }), 0d);
        }

        #endregion
    }
}
