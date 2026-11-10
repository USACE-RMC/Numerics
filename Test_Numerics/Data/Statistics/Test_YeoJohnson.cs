using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Data.Statistics;
using System.Collections.Generic;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for Yeo-Johnson transformation fitting and likelihood behavior.
    /// </summary>
    [TestClass]
    public class Test_YeoJohnson
    {
        /// <summary>
        /// FitLambda returns a finite lambda for a non-degenerate finite sample.
        /// </summary>
        [TestMethod]
        public void Test_FitLambda_ValidSample_ReturnsFiniteLambda()
        {
            var sample = new[] { -2d, -1d, -0.25d, 0d, 0.5d, 1d, 3d };

            YeoJohnson.FitLambda(sample, out double lambda);

            Assert.IsTrue(Tools.IsFinite(lambda));
            Assert.IsTrue(lambda >= -5d && lambda <= 5d);
        }

        /// <summary>
        /// FitLambda reports invalid or degenerate samples with NaN instead of throwing through BrentSearch.
        /// </summary>
        [TestMethod]
        public void Test_FitLambda_InvalidSamples_ReturnsNaN()
        {
            IList<double>[] samples =
            {
                null!,
                new[] { 1d },
                new[] { 1d, 1d, 1d },
                new[] { 1d, double.NaN, 2d },
                new[] { 1d, double.PositiveInfinity, 2d },
                new[] { -double.MaxValue, -double.MaxValue / 2d, -double.MaxValue / 4d }
            };

            foreach (IList<double> sample in samples)
            {
                YeoJohnson.FitLambda(sample, out double lambda);
                Assert.IsTrue(double.IsNaN(lambda));
            }
        }

        /// <summary>
        /// LogLikelihood returns negative infinity for unsupported samples or lambda values.
        /// </summary>
        [TestMethod]
        public void Test_LogLikelihood_InvalidSamples_ReturnsNegativeInfinity()
        {
            Assert.AreEqual(double.NegativeInfinity, YeoJohnson.LogLikelihood(new[] { 1d, 1d }, 1d));
            Assert.AreEqual(double.NegativeInfinity, YeoJohnson.LogLikelihood(new[] { 1d, double.NaN }, 1d));
            Assert.AreEqual(double.NegativeInfinity, YeoJohnson.LogLikelihood(new[] { 1d, 2d }, double.NaN));
            Assert.AreEqual(
                double.NegativeInfinity,
                YeoJohnson.LogLikelihood(new[] { -double.MaxValue, -double.MaxValue / 2d, -double.MaxValue / 4d }, 1d));
        }

        /// <summary>
        /// The log-likelihood is the full normal log-likelihood of the transformed values at the maximum likelihood
        /// mean and variance plus the log-Jacobian: scipy 1.17.1 sum(norm.logpdf(yeojohnson(x, 0.7), mu, sigma)) +
        /// (0.7 - 1) sum(sign(x) log1p(|x|)) = -16.78847767213403.
        /// </summary>
        [TestMethod]
        public void Test_LogLikelihood_FullNormalConstant()
        {
            var x = new double[] { -1.5, -0.3, 0.2, 0.8, 1.4, 2.2, 3.9, 5.1 };
            Assert.AreEqual(-16.78847767213403, YeoJohnson.LogLikelihood(x, 0.7), 1E-10);
        }
    }
}
