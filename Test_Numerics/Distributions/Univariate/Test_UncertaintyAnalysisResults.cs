using System;
using System.Linq;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Tests for <see cref="UncertaintyAnalysisResults"/>: the mean curve agrees with
    /// <see cref="BootstrapAnalysis.Estimate"/> for continuous and discrete ensembles, and the goodness-of-fit
    /// metrics default to NaN and survive XML and JSON round trips.
    /// </summary>
    [TestClass]
    public class Test_UncertaintyAnalysisResults
    {
        /// <summary>The non-exceedance probabilities evaluated in every case.</summary>
        private static readonly double[] Probabilities = { 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999 };

        /// <summary>A GEV ensemble of 200 members with distinct, deterministic parameters.</summary>
        private static UnivariateDistributionBase[] GevEnsemble()
        {
            var ensemble = new UnivariateDistributionBase[200];
            for (int i = 0; i < ensemble.Length; i++)
                ensemble[i] = new GeneralizedExtremeValue(100d + 0.25 * i, 20d + (i % 7), -0.1 + 0.002 * (i % 11));
            return ensemble;
        }

        /// <summary>A Poisson ensemble of 200 members with rates spread evenly over [3, 7].</summary>
        private static UnivariateDistributionBase[] PoissonEnsemble()
        {
            var ensemble = new UnivariateDistributionBase[200];
            for (int i = 0; i < ensemble.Length; i++)
                ensemble[i] = new Poisson(3d + 4d * i / 199d);
            return ensemble;
        }

        /// <summary>The mean curve computed by <see cref="BootstrapAnalysis.Estimate"/> for the same ensemble.</summary>
        /// <param name="ensemble">The ensemble passed as precomputed bootstrap distributions.</param>
        private static double[] BootstrapMeanCurve(UnivariateDistributionBase[] ensemble)
        {
            var boot = new BootstrapAnalysis(new Normal(100d, 20d), ParameterEstimationMethod.MethodOfMoments, 50, 100);
            return boot.Estimate(Probabilities, 0.1, ensemble.Cast<IUnivariateDistribution>().ToArray(), false).MeanCurve!;
        }

        /// <summary>Asserts two curves are bitwise identical.</summary>
        /// <param name="expected">The reference curve.</param>
        /// <param name="actual">The curve under test.</param>
        private static void AssertBitwise(double[] expected, double[] actual)
        {
            Assert.HasCount(expected.Length, actual);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i]), BitConverter.DoubleToInt64Bits(actual[i]), $"Index {i}.");
        }

        /// <summary>
        /// For a continuous ensemble the mean curve equals the bootstrap mean curve bit for bit.
        /// </summary>
        /// <remarks>
        /// The expected values are the mean curve of this ensemble from the implementation that excluded the first
        /// grid point from the log-transform decision; for a continuous ensemble both decisions agree.
        /// </remarks>
        [TestMethod]
        public void MeanCurve_ContinuousEnsemble_BitIdentical()
        {
            var ensemble = GevEnsemble();
            var results = new UncertaintyAnalysisResults(ensemble[0], ensemble, Probabilities);
            var expected = new double[] { 67.76741657334799, 79.37171905272966, 92.21957589160193, 100.28024357759736, 111.17585685628735, 134.94639370381734, 165.01121449164725, 185.73194847192553, 206.81303931152775, 259.98118501394435, 350.6624427152344 };
            AssertBitwise(expected, results.MeanCurve!);
            AssertBitwise(BootstrapMeanCurve(ensemble), results.MeanCurve!);
        }

        /// <summary>
        /// For a discrete ensemble whose grid starts at zero, the mean curve equals the bootstrap mean curve.
        /// </summary>
        /// <remarks>
        /// The Poisson quantile at probability 0.001 is zero, so the shifted logarithmic grid starts at zero and
        /// the curve must be interpolated on a linear scale. Excluding the first knot from that decision would select
        /// the log scale and give 0.00716 at p = 0.05 instead of 0.867. The mixture has P(X = 0) = 0.0123 and
        /// P(X ≤ 1) = 0.0604 (numpy mean of exp(-λ) and scipy.stats.poisson.cdf(1, λ) over the 200 rates), so its
        /// 0.05 quantile lies between the knots at 0 and 1.
        /// </remarks>
        [TestMethod]
        public void MeanCurve_DiscreteEnsemble_MatchesBootstrapAnalysis()
        {
            var ensemble = PoissonEnsemble();
            var results = new UncertaintyAnalysisResults(new Poisson(5d), ensemble, Probabilities);
            AssertBitwise(BootstrapMeanCurve(ensemble), results.MeanCurve!);

            Assert.AreEqual(0.8669780960321666, results.MeanCurve![2], 0d);
            Assert.IsGreaterThan(0d, results.MeanCurve[2]);
            Assert.IsLessThanOrEqualTo(1d, results.MeanCurve[2]);
            for (int i = 1; i < Probabilities.Length; i++)
                Assert.IsGreaterThanOrEqualTo(results.MeanCurve[i - 1], results.MeanCurve[i]);
        }

        /// <summary>
        /// The default constructor and <see cref="BootstrapAnalysis.Estimate"/> leave the metrics as NaN.
        /// </summary>
        [TestMethod]
        public void DefaultConstructor_MetricsAreNaN()
        {
            var results = new UncertaintyAnalysisResults();
            Assert.IsTrue(double.IsNaN(results.AIC));
            Assert.IsTrue(double.IsNaN(results.BIC));
            Assert.IsTrue(double.IsNaN(results.DIC));
            Assert.IsTrue(double.IsNaN(results.RMSE));
            Assert.IsTrue(double.IsNaN(results.ERL));

            var ensemble = GevEnsemble();
            var boot = new BootstrapAnalysis(new Normal(100d, 20d), ParameterEstimationMethod.MethodOfMoments, 50, 100);
            var estimate = boot.Estimate(Probabilities, 0.1, ensemble.Cast<IUnivariateDistribution>().ToArray(), false);
            Assert.IsTrue(double.IsNaN(estimate.AIC));
            Assert.IsTrue(double.IsNaN(estimate.BIC));
            Assert.IsTrue(double.IsNaN(estimate.DIC));
            Assert.IsTrue(double.IsNaN(estimate.RMSE));
            Assert.IsTrue(double.IsNaN(estimate.ERL));
        }

        /// <summary>
        /// NaN metrics survive the XML and JSON round trips, and assigned values are preserved.
        /// </summary>
        [TestMethod]
        public void RoundTrips_PreserveNaNAndValues()
        {
            var results = new UncertaintyAnalysisResults
            {
                ParentDistribution = new Normal(10d, 2d),
                ModeCurve = new[] { 1d, 2d },
                MeanCurve = new[] { 1.5, 2.5 },
                ConfidenceIntervals = new double[,] { { 0.5, 1.5 }, { 1d, 3d } },
                RMSE = 0.125
            };

            var fromXml = UncertaintyAnalysisResults.FromXElement(results.ToXElement());
            Assert.IsTrue(double.IsNaN(fromXml.AIC));
            Assert.IsTrue(double.IsNaN(fromXml.BIC));
            Assert.IsTrue(double.IsNaN(fromXml.DIC));
            Assert.IsTrue(double.IsNaN(fromXml.ERL));
            Assert.AreEqual(0.125, fromXml.RMSE, 0d);
            CollectionAssert.AreEqual(results.MeanCurve, fromXml.MeanCurve);

            var fromJson = UncertaintyAnalysisResults.FromByteArray(UncertaintyAnalysisResults.ToByteArray(results));
            Assert.IsNotNull(fromJson);
            Assert.IsTrue(double.IsNaN(fromJson.AIC));
            Assert.IsTrue(double.IsNaN(fromJson.BIC));
            Assert.IsTrue(double.IsNaN(fromJson.DIC));
            Assert.IsTrue(double.IsNaN(fromJson.ERL));
            Assert.AreEqual(0.125, fromJson.RMSE, 0d);
            CollectionAssert.AreEqual(results.MeanCurve, fromJson.MeanCurve);
        }

        /// <summary>
        /// XML written with zero-valued metric attributes still reads back as zero.
        /// </summary>
        [TestMethod]
        public void FromXElement_ZeroAttributes_ReadAsZero()
        {
            var xml = XElement.Parse("<UncertaintyAnalysisResults AIC=\"0\" BIC=\"0\" DIC=\"0\" RMSE=\"0\" ERL=\"0\" ModeCurve=\"1|2\" MeanCurve=\"1.5|2.5\" ConfidenceIntervals=\"0.5,1.5|1,3\" />");
            var results = UncertaintyAnalysisResults.FromXElement(xml);
            Assert.AreEqual(0d, results.AIC, 0d);
            Assert.AreEqual(0d, results.BIC, 0d);
            Assert.AreEqual(0d, results.DIC, 0d);
            Assert.AreEqual(0d, results.RMSE, 0d);
            Assert.AreEqual(0d, results.ERL, 0d);
            CollectionAssert.AreEqual(new[] { 1.5, 2.5 }, results.MeanCurve);
            Assert.AreEqual(3d, results.ConfidenceIntervals![1, 1], 0d);
        }
    }
}
