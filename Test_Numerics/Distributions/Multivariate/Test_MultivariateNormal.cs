using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
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
         
    }
}
