using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using System;

namespace Distributions.Univariate
{
    /// <summary>
    /// Unit tests for the closed-form log density overrides of ten univariate distributions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// Each distribution is checked at an ordinary point and at a far-tail point; where the
    /// far-tail density underflows double precision, the log density must stay finite where the
    /// density itself cannot. Reference values computed with Python scipy 1.17.1 (scipy.stats.&lt;dist&gt;.logpdf);
    /// the base-10 logarithmic distributions, which scipy does not provide, are pinned to 50-digit
    /// mpmath 1.4.1 evaluations of the same density and to scipy.stats.pearson3 evaluated on the
    /// base-10 logarithm with the change-of-variables term.
    /// </para>
    /// <para>
    /// Tolerances are measured: the closed-form log densities agree with the references to within
    /// 1E-12 relative error, and the assertions allow 1E-11 relative plus 1E-10 absolute slack.
    /// Every override must also reproduce its own density: exp(LogPDF(x)) is compared to PDF(x)
    /// at the ordinary points to 1E-12 relative error.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_LogPDFOverrides
    {
        /// <summary>
        /// Asserts a log density against its reference value with mixed relative and absolute slack.
        /// </summary>
        /// <param name="expected">The reference log density.</param>
        /// <param name="actual">The evaluated log density.</param>
        private static void AssertLogDensity(double expected, double actual)
        {
            Assert.AreEqual(expected, actual, Math.Abs(expected) * 1E-11 + 1E-10);
        }

        /// <summary>
        /// Asserts that the exponential of the log density reproduces the density at a point where
        /// the density is representable.
        /// </summary>
        /// <param name="dist">The distribution under test.</param>
        /// <param name="x">The evaluation point.</param>
        private static void AssertConsistentWithPDF(UnivariateDistributionBase dist, double x)
        {
            double pdf = dist.PDF(x);
            Assert.AreEqual(pdf, Math.Exp(dist.LogPDF(x)), pdf * 1E-12);
        }

        /// <summary>
        /// The normal log density matches scipy at an ordinary point and stays finite forty
        /// standard deviations out.
        /// </summary>
        [TestMethod]
        public void Test_Normal_LogPDF()
        {
            var d = new Normal(10, 2);
            AssertLogDensity(-2.112085713764618, d.LogPDF(12.0));
            AssertLogDensity(-801.6120857137646, d.LogPDF(90.0));
            AssertConsistentWithPDF(d, 12.0);
        }

        /// <summary>
        /// The natural-log normal log density matches scipy lognorm(s = sigma, scale = exp(mu)) at
        /// an ordinary point and stays finite
        /// where the density underflows, and is negative infinity at and below the support bound.
        /// </summary>
        [TestMethod]
        public void Test_LnNormal_LogPDF()
        {
            var d = new LnNormal() { Mu = 2, Sigma = 0.5 };
            AssertLogDensity(-2.317854811413502, d.LogPDF(8.0));
            AssertLogDensity(-907.4244569387412, d.LogPDF(1E10));
            AssertConsistentWithPDF(d, 8.0);
            Assert.AreEqual(double.NegativeInfinity, d.LogPDF(0.0));
            Assert.AreEqual(double.NegativeInfinity, d.LogPDF(-1.0));
        }

        /// <summary>
        /// The base-10 log-normal log density matches the 50-digit mpmath evaluation at an
        /// ordinary point and deep in the upper tail.
        /// </summary>
        [TestMethod]
        public void Test_LogNormal_LogPDF()
        {
            var d = new LogNormal() { Mu = 1.2, Sigma = 0.25 };
            AssertLogDensity(-3.4440653710776417, d.LogPDF(20.0));
            AssertLogDensity(-388.7073573612851, d.LogPDF(1E8));
            AssertConsistentWithPDF(d, 20.0);
            Assert.AreEqual(double.NegativeInfinity, d.LogPDF(0.0));
        }

        /// <summary>
        /// The Pearson Type III log density matches scipy pearson3 on both skew signs and stays
        /// finite thirty standard deviations into the upper tail.
        /// </summary>
        [TestMethod]
        public void Test_PearsonTypeIII_LogPDF()
        {
            var d = new PearsonTypeIII(10, 2, 0.8);
            AssertLogDensity(-2.6217159334295284, d.LogPDF(12.5));
            AssertLogDensity(-63.159423624324326, d.LogPDF(70.0));
            AssertConsistentWithPDF(d, 12.5);

            var dn = new PearsonTypeIII(10, 2, -0.8);
            AssertLogDensity(-2.1394304489371043, dn.LogPDF(12.5));
            AssertConsistentWithPDF(dn, 12.5);
        }

        /// <summary>
        /// The log-Pearson Type III log density matches scipy pearson3 evaluated on the base-10
        /// logarithm with the change-of-variables term.
        /// </summary>
        [TestMethod]
        public void Test_LogPearsonTypeIII_LogPDF()
        {
            var d = new LogPearsonTypeIII(1.5, 0.3, 0.6);
            AssertLogDensity(-4.844505860354813, d.LogPDF(50.0));
            AssertLogDensity(-58.859703127182755, d.LogPDF(1E7));
            AssertConsistentWithPDF(d, 50.0);
            Assert.AreEqual(double.NegativeInfinity, d.LogPDF(0.0));
        }

        /// <summary>
        /// The generalized extreme value log density matches scipy genextreme and stays finite deep
        /// in its lower tail.
        /// </summary>
        [TestMethod]
        public void Test_GeneralizedExtremeValue_LogPDF()
        {
            var d = new GeneralizedExtremeValue(100, 25, 0.15);
            AssertLogDensity(-4.60976457450874, d.LogPDF(130.0));
            AssertLogDensity(-6133.58747825334, d.LogPDF(-350.0));
            AssertConsistentWithPDF(d, 130.0);
        }

        /// <summary>
        /// The Gumbel log density matches scipy gumbel_r and stays finite where the double
        /// exponential collapses the density.
        /// </summary>
        [TestMethod]
        public void Test_Gumbel_LogPDF()
        {
            var d = new Gumbel(100, 25);
            AssertLogDensity(-4.7200700367804025, d.LogPDF(130.0));
            AssertLogDensity(-78962960182651.9, d.LogPDF(-700.0));
            AssertConsistentWithPDF(d, 130.0);
        }

        /// <summary>
        /// The Weibull log density matches scipy weibull_min at an ordinary point and deep in the
        /// upper tail, and is negative infinity below the support.
        /// </summary>
        [TestMethod]
        public void Test_Weibull_LogPDF()
        {
            var d = new Weibull(50, 2.5);
            AssertLogDensity(-3.902881002765252, d.LogPDF(40.0));
            AssertLogDensity(-732.4019940868236, d.LogPDF(700.0));
            AssertConsistentWithPDF(d, 40.0);
            Assert.AreEqual(double.NegativeInfinity, d.LogPDF(-1.0));
        }

        /// <summary>
        /// The exponential log density matches scipy expon and stays finite where the density
        /// underflows, and is negative infinity below the location bound.
        /// </summary>
        [TestMethod]
        public void Test_Exponential_LogPDF()
        {
            var d = new Exponential(5, 20);
            AssertLogDensity(-4.24573227355399, d.LogPDF(30.0));
            AssertLogDensity(-752.745732273554, d.LogPDF(15000.0));
            AssertConsistentWithPDF(d, 30.0);
            Assert.AreEqual(double.NegativeInfinity, d.LogPDF(4.0));
        }

        /// <summary>
        /// The gamma log density matches scipy and stays finite deep in the upper tail.
        /// </summary>
        [TestMethod]
        public void Test_Gamma_LogPDF()
        {
            var d = new GammaDistribution(15, 3);
            AssertLogDensity(-4.014903020542265, d.LogPDF(30.0));
            AssertLogDensity(-1322.3436560121274, d.LogPDF(20000.0));
            AssertConsistentWithPDF(d, 30.0);
        }

        /// <summary>
        /// A single far-tail observation keeps the log likelihood of a sample finite and equal to
        /// the sum of the log densities.
        /// </summary>
        [TestMethod]
        public void Test_LogLikelihood_FarTailObservation()
        {
            var d = new Normal(0, 1);
            double logLH = d.LogLikelihood(new double[] { 0.0, 50.0 });
            AssertLogDensity(-0.9189385332046727 - 1250.9189385332047, logLH);
        }
    }
}
