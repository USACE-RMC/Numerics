using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Distributions;
using Numerics.Distributions.Regional;

namespace Distributions.Regional
{
    /// <summary>
    /// Tests of the Hosking-Wallis regional frequency analysis: regional data, discordancy, heterogeneity,
    /// goodness of fit and the regional growth curve.
    /// </summary>
    /// <remarks>
    /// The oracle is the North Cascades example of Hosking's LMOMENTS Fortran package, version 3.04
    /// (https://lib.stat.cmu.edu/general/lmoments): CASCADES.DAT, the L-moments of annual precipitation at 19 sites,
    /// and CASCADES.OUT, the output of PROGRAM XTEST (REGTST with 500 simulations) on those data. The same data are
    /// lmomRFA's Cascades data set, and Hosking and Wallis (1997) discuss the example in their Table 3.4. Fixed values
    /// are checked to the printed precision; H and Z are checked within Monte Carlo tolerance, because the random
    /// number generators differ.
    /// </remarks>
    [TestClass]
    public class Test_RegionalFrequencyAnalysis
    {
        /// <summary>The CASCADES.DAT site identifiers.</summary>
        private static readonly string[] Names =
        [
            "350304", "351433", "351862", "351897", "352997", "353445", "353770", "356907", "357169", "357331",
            "357354", "358466", "450945", "451233", "453284", "454764", "454769", "457773", "458773",
        ];

        /// <summary>The CASCADES.DAT record lengths.</summary>
        private static readonly int[] N = [98, 59, 90, 61, 65, 86, 78, 72, 67, 99, 49, 61, 69, 73, 70, 66, 59, 74, 82];

        /// <summary>The CASCADES.DAT site means.</summary>
        private static readonly double[] Mean =
        [
            19.685, 62.580, 40.852, 46.045, 45.021, 31.042, 80.143, 41.305, 30.585, 32.932,
            17.560, 69.518, 47.653, 102.501, 52.413, 79.696, 44.643, 58.655, 39.024,
        ];

        /// <summary>The CASCADES.DAT L-CV.</summary>
        private static readonly double[] LCV =
        [
            0.1209, 0.0915, 0.1124, 0.1032, 0.0967, 0.1328, 0.1008, 0.1143, 0.1107, 0.1179,
            0.1308, 0.1119, 0.1018, 0.1025, 0.1054, 0.1174, 0.1115, 0.1003, 0.1046,
        ];

        /// <summary>The CASCADES.DAT L-skewness.</summary>
        private static readonly double[] T3 =
        [
            0.0488, 0.0105, 0.0614, 0.0417, -0.0134, -0.0176, 0.0943, 0.0555, 0.0478, 0.0492,
            0.0940, -0.0429, 0.0435, 0.0182, -0.0224, 0.0124, -0.0346, 0.0446, 0.0128,
        ];

        /// <summary>The CASCADES.DAT L-kurtosis.</summary>
        private static readonly double[] T4 =
        [
            0.1433, 0.1569, 0.1541, 0.1429, 0.1568, 0.1206, 0.1967, 0.1210, 0.1371, 0.0900,
            0.1273, 0.0927, 0.1446, 0.1047, 0.1664, 0.1317, 0.1032, 0.1450, 0.1583,
        ];

        /// <summary>The CASCADES.DAT t5.</summary>
        private static readonly double[] T5 =
        [
            -0.0004, 0.0020, -0.0058, -0.0022, 0.0173, 0.0235, 0.0856, 0.0487, 0.0316, 0.0225,
            0.0352, -0.0061, -0.0056, -0.0221, 0.0035, -0.0176, 0.0083, -0.0379, 0.0443,
        ];

        /// <summary>The quantile probabilities of CASCADES.OUT.</summary>
        private static readonly double[] Probabilities = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.9, 0.95, 0.99, 0.999];

        /// <summary>The shared test of the Cascades region with 500 simulations.</summary>
        private static readonly Lazy<RegionalHomogeneityTest> _cascadesTest = new(() => new RegionalHomogeneityTest(Cascades(), 500, 12345));

        /// <summary>Returns the Cascades regional data.</summary>
        private static RegionalData Cascades() => new RegionalData(Names, N, Mean, LCV, T3, T4, T5);

        /// <summary>Returns the Cascades test with 500 simulations, computed once.</summary>
        private static RegionalHomogeneityTest CascadesTest() => _cascadesTest.Value;

        /// <summary>
        /// The discordancy measures match CASCADES.OUT to two decimals, and their mean is 1.
        /// </summary>
        [TestMethod]
        public void Discordancy_MatchesCascades()
        {
            double[] expected = [0.60, 1.02, 0.38, 0.23, 0.93, 2.63, 2.12, 0.45, 0.11, 1.61, 2.08, 1.52, 0.31, 1.30, 1.58, 0.29, 1.04, 0.43, 0.38];
            var d = Cascades().Discordancy();
            double sum = 0d;
            for (int i = 0; i < expected.Length; i++)
            {
                Assert.AreEqual(expected[i], d[i], 0.005, $"site {Names[i]}");
                sum += d[i];
            }
            Assert.AreEqual(1d, sum / d.Length, 1E-12);
        }

        /// <summary>
        /// The record-length weighted regional averages match CASCADES.OUT (WEIGHTED MEANS 0.1103 0.0279 0.1366).
        /// </summary>
        [TestMethod]
        public void RegionalAverage_MatchesCascades()
        {
            var r = Cascades().RegionalAverage();
            Assert.HasCount(5, r);
            Assert.AreEqual(1d, r[0]);
            Assert.AreEqual(0.1103, r[1], 0.00005);
            Assert.AreEqual(0.0279, r[2], 0.00005);
            Assert.AreEqual(0.1366, r[3], 0.00005);
            Assert.IsLessThan(0.05, Math.Abs(r[4]));
        }

        /// <summary>
        /// Equal weights give the unweighted means, and a t5 of NaN leaves the site out of t5^R.
        /// </summary>
        [TestMethod]
        public void RegionalAverage_WeightsAndMissingTau5()
        {
            var data = new RegionalData(["a", "b"], [10, 30], [1, 2], [0.1, 0.2], [0.0, 0.1], [0.1, 0.2], [double.NaN, 0.05]);
            var r = data.RegionalAverage();
            Assert.AreEqual(0.175, r[1], 1E-15);
            Assert.AreEqual(0.075, r[2], 1E-15);
            Assert.AreEqual(0.175, r[3], 1E-15);
            Assert.AreEqual(0.05, r[4], 1E-15);
            var equal = data.RegionalAverage([1d, 1d]);
            Assert.AreEqual(0.15, equal[1], 1E-15);
            Assert.ThrowsExactly<ArgumentException>(() => data.RegionalAverage([0d, 0d]));
            Assert.ThrowsExactly<ArgumentException>(() => data.RegionalAverage([1d]));
            Assert.ThrowsExactly<ArgumentException>(() => data.RegionalAverage([-1d, 2d]));
        }

        /// <summary>
        /// The critical values follow Table 3.1 of Hosking and Wallis (1997) and lmomRFA's 5% table.
        /// </summary>
        [TestMethod]
        public void DiscordancyCriticalValue_Table()
        {
            Assert.AreEqual(3d, RegionalData.DiscordancyCriticalValue(3));
            Assert.AreEqual(1.333, RegionalData.DiscordancyCriticalValue(5), 1E-12);
            Assert.AreEqual(2.1401, RegionalData.DiscordancyCriticalValue(8), 1E-12);
            Assert.AreEqual(2.9709, RegionalData.DiscordancyCriticalValue(14), 1E-12);
            Assert.AreEqual(3d, RegionalData.DiscordancyCriticalValue(15));
            Assert.AreEqual(3d, RegionalData.DiscordancyCriticalValue(500));
            Assert.AreEqual(1.3333, RegionalData.DiscordancyCriticalValue(5, 0.05), 1E-12);
            Assert.AreEqual(3.9816, RegionalData.DiscordancyCriticalValue(18, 0.05), 1E-12);
            Assert.AreEqual(4d, RegionalData.DiscordancyCriticalValue(19, 0.05));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => RegionalData.DiscordancyCriticalValue(0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => RegionalData.DiscordancyCriticalValue(10, 0.01));
        }

        /// <summary>
        /// Small regions give D_i = 1, as REGTST does. Coplanar ratios make the cross-product matrix singular and give
        /// NaN, a deliberate difference from REGTST, which leaves D_i = 0; lmomRFA also reports NA.
        /// </summary>
        [TestMethod]
        public void Discordancy_EdgeCases()
        {
            var three = new RegionalData(["a", "b", "c"], [10, 10, 10], [1, 1, 1], [0.1, 0.2, 0.3], [0.0, 0.1, 0.0], [0.1, 0.2, 0.1]);
            CollectionAssert.AreEqual(new[] { 1d, 1d, 1d }, three.Discordancy());
            var flat = new RegionalData(["a", "b", "c", "d", "e"], [10, 10, 10, 10, 10], [1, 1, 1, 1, 1],
                [0.1, 0.2, 0.3, 0.4, 0.5], [0.0, 0.1, 0.2, 0.3, 0.4], [0.1, 0.1, 0.1, 0.1, 0.1]);
            foreach (double d in flat.Discordancy()) Assert.IsTrue(double.IsNaN(d));
        }

        /// <summary>
        /// Raw samples give the exact unbiased sample L-moments, including t5, and missing values are left out.
        /// </summary>
        /// <remarks>
        /// Oracle: exact rational arithmetic of the unbiased probability weighted moments (Python fractions, printed
        /// with mpmath 1.3.0 to 20 digits).
        /// </remarks>
        [TestMethod]
        public void FromSamples_MatchesExactLinearMoments()
        {
            var a = new List<double> { 12.1, 15.3, 9.8, 22.4, 18.0, 13.7, 11.2, 25.9, 16.4, 14.8 };
            var b = new List<double> { 3.2, double.NaN, 4.1, 2.7, 5.6, 3.9, 4.4, 6.8, 3.1, 4.9 };
            var data = RegionalData.FromSamples(new List<IList<double>> { a, b }, ["A", "B"]);
            Assert.AreEqual(10, data.RecordLengths[0]);
            Assert.AreEqual(9, data.RecordLengths[1]);
            Assert.AreEqual(15.96, data.Mean[0], 1E-13);
            Assert.AreEqual(2.9044444444444444444 / 15.96, data.LCV[0], 1E-14);
            Assert.AreEqual(0.23297628156082631982, data.Tau3[0], 1E-13);
            Assert.AreEqual(0.18297081648267570226, data.Tau4[0], 1E-13);
            Assert.AreEqual(0.045524100994644223412, data.Tau5[0], 1E-13);
            Assert.AreEqual(4.3, data.Mean[1], 1E-13);
            Assert.AreEqual(0.77222222222222222222 / 4.3, data.LCV[1], 1E-14);
            Assert.AreEqual(0.20041109969167523124, data.Tau3[1], 1E-13);
            Assert.AreEqual(0.13669064748201438849, data.Tau4[1], 1E-13);
            Assert.AreEqual(0.099691675231243576567, data.Tau5[1], 1E-13);
            var site = data.SiteLinearMoments(0);
            Assert.AreEqual(2.9044444444444444444, site[1], 1E-13);
            Assert.AreEqual("B", data.Names[1]);

            var unnamed = RegionalData.FromSamples(new List<IList<double>> { a, new List<double> { 1, 2, 3, 4 } });
            Assert.AreEqual("2", unnamed.Names[1]);
            Assert.IsTrue(double.IsNaN(unnamed.Tau5[1]));
        }

        /// <summary>
        /// Invalid regional data are refused.
        /// </summary>
        [TestMethod]
        public void RegionalData_InvalidInputs()
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => RegionalData.FromSamples(null!));
            Assert.ThrowsExactly<ArgumentException>(() => RegionalData.FromSamples(new List<IList<double>>()));
            Assert.ThrowsExactly<ArgumentException>(() => RegionalData.FromSamples(new List<IList<double>> { new List<double> { 1, 2, 3 } }));
            Assert.ThrowsExactly<ArgumentException>(() => RegionalData.FromSamples(new List<IList<double>> { new List<double> { 1, 2, 3, double.PositiveInfinity } }));
            Assert.ThrowsExactly<ArgumentException>(() => RegionalData.FromSamples(new List<IList<double>> { new List<double> { -1, -2, -3, -4 } }));
            Assert.ThrowsExactly<ArgumentException>(() => new RegionalData(["a"], [10, 20], [1], [0.1], [0], [0.1]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalData(["a"], [0], [1], [0.1], [0], [0.1]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalData(["a"], [10], [1], [double.NaN], [0], [0.1]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Cascades().SiteLinearMoments(19));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalData(["a"], [10], [1], [0.0], [0], [0.1]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalData(["a"], [10], [1], [-0.1], [0], [0.1]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalData(["a"], [10], [1], [0.1], [0], [0.1], [double.PositiveInfinity]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => RegionalData.FromSamples(new List<IList<double>> { new List<double> { 2, 2, 2, 2, 2 } }));
            var missingTau5 = new RegionalData(["a"], [10], [1], [0.1], [0], [0.1], [double.NaN]);
            Assert.IsTrue(double.IsNaN(missingTau5.Tau5[0]));
        }

        /// <summary>
        /// The kappa parent matches CASCADES.OUT (0.9542 0.1533 0.1236 -0.2955), and the observed V1 to V3 match the
        /// observed dispersions of CASCADES.OUT (0.0104, 0.0339, 0.0405).
        /// </summary>
        [TestMethod]
        public void HomogeneityTest_FixedValuesMatchCascades()
        {
            var test = CascadesTest();
            Assert.IsFalse(test.IsGeneralizedLogisticParent);
            Assert.AreEqual(0.9542, test.KappaParameters[0], 0.00005);
            Assert.AreEqual(0.1533, test.KappaParameters[1], 0.00005);
            Assert.AreEqual(0.1236, test.KappaParameters[2], 0.00005);
            Assert.AreEqual(-0.2955, test.KappaParameters[3], 0.00005);
            Assert.AreEqual(0.0104, test.ObservedDispersion[0], 0.00005);
            Assert.AreEqual(0.0339, test.ObservedDispersion[1], 0.00005);
            Assert.AreEqual(0.0405, test.ObservedDispersion[2], 0.00005);
            Assert.AreEqual(3d, test.DiscordancyCriticalValue);

            // The L-kurtosis of the fitted candidates, printed to three decimals.
            double[] tau4 = [0.167, 0.111, 0.123, 0.123, 0.006];
            for (int d = 0; d < tau4.Length; d++)
                Assert.AreEqual(tau4[d], test.DistributionTau4[d], 0.0005, test.Distributions[d].ToString());
        }

        /// <summary>
        /// H and Z lie within Monte Carlo tolerance of CASCADES.OUT: H = 0.62, -1.49, -2.37 and Z = 3.46, -2.94,
        /// -1.51, -1.60, -14.75; the simulated means and standard deviations of V are 0.0095, 0.0447, 0.0578 and
        /// 0.0016, 0.0072, 0.0073.
        /// </summary>
        /// <remarks>
        /// With 500 simulations the standard error of a simulated mean is about σ/√500 (0.045 in H units), and of a
        /// simulated standard deviation about 3% (so about 0.03|H|); the published values carry the same error. The
        /// tolerances 0.3 + 0.1|H| and 0.3 + 0.1|Z| are about three combined standard errors.
        /// </remarks>
        [TestMethod]
        public void HomogeneityTest_MonteCarloValuesMatchCascades()
        {
            var test = CascadesTest();
            double[] h = [0.62, -1.49, -2.37];
            for (int j = 0; j < 3; j++)
                Assert.AreEqual(h[j], test.Heterogeneity[j], 0.3 + 0.1 * Math.Abs(h[j]), $"H{j + 1}");
            double[] vbar = [0.0095, 0.0447, 0.0578];
            double[] vsd = [0.0016, 0.0072, 0.0073];
            for (int j = 0; j < 3; j++)
            {
                Assert.AreEqual(vbar[j], test.SimulatedDispersionMean[j], 0.1 * vbar[j], $"mean V{j + 1}");
                Assert.AreEqual(vsd[j], test.SimulatedDispersionStandardDeviation[j], 0.2 * vsd[j], $"sd V{j + 1}");
            }
            double[] z = [3.46, -2.94, -1.51, -1.60, -14.75];
            for (int d = 0; d < z.Length; d++)
                Assert.AreEqual(z[d], test.Z[d], 0.3 + 0.1 * Math.Abs(z[d]), test.Distributions[d].ToString());
            Assert.AreEqual(test.Z[2], test.GetZ(UnivariateDistributionType.GeneralizedNormal));
            var accepted = test.AcceptableDistributions();
            CollectionAssert.AreEquivalent(new[] { UnivariateDistributionType.GeneralizedNormal, UnivariateDistributionType.PearsonTypeIII }, accepted);
        }

        /// <summary>
        /// A seed gives the same results whatever the number of threads, and another seed gives different results.
        /// </summary>
        [TestMethod]
        public void HomogeneityTest_DeterministicAcrossThreadCounts()
        {
            var data = Cascades();
            var serial = new RegionalHomogeneityTest(data, 200, 777, 1);
            var parallel = new RegionalHomogeneityTest(data, 200, 777, 8);
            var automatic = new RegionalHomogeneityTest(data, 200, 777);
            CollectionAssert.AreEqual(serial.Heterogeneity, parallel.Heterogeneity);
            CollectionAssert.AreEqual(serial.Heterogeneity, automatic.Heterogeneity);
            CollectionAssert.AreEqual(serial.Z, parallel.Z);
            Assert.AreEqual(serial.Tau4Bias, automatic.Tau4Bias);
            var other = new RegionalHomogeneityTest(data, 200, 778);
            Assert.AreNotEqual(serial.Heterogeneity[0], other.Heterogeneity[0]);
        }

        /// <summary>
        /// Regional L-moments above the generalized logistic line fall back to the generalized logistic parent.
        /// </summary>
        [TestMethod]
        public void HomogeneityTest_GeneralizedLogisticFallback()
        {
            // t4 = 0.30 lies above the GLO line (1 + 5 t3^2) / 6 = 0.1767 at t3 = 0.1.
            var data = new RegionalData(["a", "b", "c", "d", "e"], [40, 50, 60, 45, 55], [10, 12, 9, 11, 10],
                [0.20, 0.22, 0.18, 0.21, 0.19], [0.08, 0.12, 0.10, 0.09, 0.11], [0.29, 0.31, 0.30, 0.32, 0.28]);
            var test = new RegionalHomogeneityTest(data, 100, 1);
            Assert.IsTrue(test.IsGeneralizedLogisticParent);
            Assert.AreEqual(-1d, test.KappaParameters[3]);
            // The generalized logistic shape is -t3^R.
            Assert.AreEqual(-test.RegionalLinearMoments[2], test.KappaParameters[2], 1E-15);
            foreach (double h in test.Heterogeneity) Assert.IsTrue(Tools.IsFinite(h));
        }

        /// <summary>
        /// The result collections are read-only, so callers cannot change the statistics the test reports.
        /// </summary>
        [TestMethod]
        public void HomogeneityTest_ResultsAreReadOnly()
        {
            var test = new RegionalHomogeneityTest(Cascades(), 50, 3);
            var accepted = test.AcceptableDistributions(2.0);
            Assert.ThrowsExactly<NotSupportedException>(() => ((IList<double>)test.Z)[2] = 100d);
            Assert.ThrowsExactly<NotSupportedException>(() => ((IList<double>)test.Heterogeneity)[0] = 0d);
            Assert.ThrowsExactly<NotSupportedException>(() => ((IList<UnivariateDistributionType>)test.Distributions)[0] = UnivariateDistributionType.Normal);
            CollectionAssert.AreEqual(accepted, test.AcceptableDistributions(2.0));
            Assert.AreEqual(0, test.NonpositiveMeanSimulations);
        }

        /// <summary>
        /// A parent with a large L-CV and short records produces simulated sites with a nonpositive mean, which the
        /// test counts; the algorithm, as REGTST's, keeps those regions.
        /// </summary>
        [TestMethod]
        public void HomogeneityTest_CountsNonpositiveSimulatedMeans()
        {
            var data = new RegionalData(["a", "b", "c", "d", "e"], [4, 5, 4, 5, 4], [1, 1, 1, 1, 1],
                [0.70, 0.72, 0.68, 0.71, 0.69], [0.00, 0.01, -0.01, 0.02, -0.02], [0.12, 0.13, 0.11, 0.12, 0.12]);
            var test = new RegionalHomogeneityTest(data, 200, 5);
            Assert.IsGreaterThan(0, test.NonpositiveMeanSimulations);
            Assert.IsLessThanOrEqualTo(200, test.NonpositiveMeanSimulations);
        }

        /// <summary>
        /// Invalid test inputs are refused.
        /// </summary>
        [TestMethod]
        public void HomogeneityTest_InvalidInputs()
        {
            var one = new RegionalData(["a"], [10], [1], [0.1], [0], [0.1]);
            Assert.ThrowsExactly<ArgumentNullException>(() => new RegionalHomogeneityTest(null!));
            Assert.ThrowsExactly<ArgumentException>(() => new RegionalHomogeneityTest(one));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalHomogeneityTest(Cascades(), 1));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalHomogeneityTest(Cascades(), 10, 1, 0));
            var shortRecord = new RegionalData(["a", "b"], [3, 10], [1, 1], [0.1, 0.1], [0, 0], [0.1, 0.1]);
            Assert.ThrowsExactly<ArgumentException>(() => new RegionalHomogeneityTest(shortRecord));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => CascadesTest().GetZ(UnivariateDistributionType.Normal));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => CascadesTest().AcceptableDistributions(0));
        }

        /// <summary>
        /// The generalized normal and Pearson type III growth curves match the parameters and quantile table of
        /// CASCADES.OUT, printed to three decimals.
        /// </summary>
        [TestMethod]
        public void GrowthCurve_MatchesCascades()
        {
            var data = Cascades();
            var gno = new RegionalGrowthCurve(UnivariateDistributionType.GeneralizedNormal, data);
            AssertClose([0.994, 0.195, -0.057], gno.Parameters, 0.0005);
            AssertClose([0.569, 0.616, 0.688, 0.753, 0.834, 0.994, 1.254, 1.331, 1.480, 1.654], gno.GrowthFactors(Probabilities), 0.0005);

            var pe3 = new RegionalGrowthCurve(UnivariateDistributionType.PearsonTypeIII, data);
            AssertClose([1.000, 0.196, 0.171], pe3.Parameters, 0.0005);
            AssertClose([0.570, 0.616, 0.688, 0.753, 0.834, 0.994, 1.254, 1.331, 1.480, 1.653], pe3.GrowthFactors(Probabilities), 0.0005);

            Assert.AreEqual(1d, gno.GetDistribution().Mean, 1E-9);
            Assert.AreSame(data, gno.Data);
            Assert.AreEqual(gno.GrowthFactor(0.99), gno.GrowthFactors([0.99])[0]);
        }

        /// <summary>
        /// Site quantiles are the site means times the growth factors, and every supported family fits.
        /// </summary>
        [TestMethod]
        public void GrowthCurve_SiteQuantilesAndFamilies()
        {
            var data = Cascades();
            var gev = new RegionalGrowthCurve(UnivariateDistributionType.GeneralizedExtremeValue, data);
            var q = gev.GrowthFactors(Probabilities);
            var all = gev.SiteQuantiles(data, Probabilities);
            var one = gev.SiteQuantiles(Mean[13], Probabilities);
            for (int j = 0; j < Probabilities.Length; j++)
            {
                Assert.AreEqual(Mean[13] * q[j], all[13, j], 1E-12);
                Assert.AreEqual(all[13, j], one[j], 1E-12);
            }
            var r = data.RegionalAverage();
            foreach (var type in new[] { UnivariateDistributionType.GeneralizedLogistic, UnivariateDistributionType.GeneralizedExtremeValue,
                UnivariateDistributionType.GeneralizedNormal, UnivariateDistributionType.PearsonTypeIII,
                UnivariateDistributionType.GeneralizedPareto, UnivariateDistributionType.KappaFour })
            {
                var curve = new RegionalGrowthCurve(type, r[1], r[2], r[3]);
                Assert.AreEqual(1d, curve.GetDistribution().Mean, 1E-6, type.ToString());
                Assert.IsGreaterThan(curve.GrowthFactor(0.5), curve.GrowthFactor(0.99), type.ToString());
            }
        }

        /// <summary>
        /// Invalid growth curve inputs are refused.
        /// </summary>
        [TestMethod]
        public void GrowthCurve_InvalidInputs()
        {
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalGrowthCurve(UnivariateDistributionType.Normal, 0.1, 0.0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalGrowthCurve(UnivariateDistributionType.GeneralizedNormal, -0.1, 0.0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new RegionalGrowthCurve(UnivariateDistributionType.KappaFour, 0.1, 0.0));
            Assert.ThrowsExactly<ArgumentNullException>(() => new RegionalGrowthCurve(UnivariateDistributionType.GeneralizedNormal, (RegionalData)null!));
            var curve = new RegionalGrowthCurve(UnivariateDistributionType.GeneralizedNormal, 0.1, 0.0);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => curve.SiteQuantiles(0d, Probabilities));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => curve.GrowthFactor(1.5));
            Assert.ThrowsExactly<ArgumentNullException>(() => curve.GrowthFactors(null!));
        }

        /// <summary>Asserts that two arrays agree element by element within a tolerance.</summary>
        private static void AssertClose(double[] expected, double[] actual, double tolerance)
        {
            Assert.HasCount(expected.Length, actual);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i], tolerance, $"index {i}");
        }
    }
}
