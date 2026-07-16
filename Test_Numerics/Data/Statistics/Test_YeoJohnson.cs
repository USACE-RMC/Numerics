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
    }
}
