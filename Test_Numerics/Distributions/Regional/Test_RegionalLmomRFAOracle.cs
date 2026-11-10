using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Distributions.Regional;

namespace Distributions.Regional
{
    /// <summary>
    /// The regional frequency analysis of the North Cascades data against the R package lmomRFA at full precision.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Oracle: R 4.4.3 (2025-02-28 ucrt), lmom 3.3 and lmomRFA 3.8, on lmomRFA's Cascades data set (the CASCADES.DAT of
    /// Hosking's LMOMENTS package): regavlmom(Cascades), regtst(Cascades, nsim = 500) after set.seed(1), and
    /// regquant(f, regfit(Cascades, dist)) for each distribution, printed with 17 significant digits.
    /// </para>
    /// <para>
    /// The deterministic values agree to rounding, except where lmom itself is less precise:
    /// <list type="bullet">
    /// <item><description>The kappa parameters come from Newton iterations that stop at a 1E-6 L-moment residual in
    /// both codes, so the two iterates agree to about 3E-12, not to rounding.</description></item>
    /// <item><description>lmom's PELGEV takes the GEV shape from Hosking's rational approximation without iteration for
    /// t3 above -0.8; it is 1.4E-7 from the exact shape here, while Numerics agrees with the exact shape (mpmath, 40
    /// digits) to 4E-15. The GEV parameters, quantiles and τ4 are therefore checked against R within that
    /// approximation error and against the exact values within 1E-13.</description></item>
    /// <item><description>lmom computes the GNO τ4 by numerical integration, exact here to 4E-16, while Numerics'
    /// GeneralizedNormal.LinearMomentsFromParameters uses Hosking's rational approximation, 1.5E-11 from it.</description></item>
    /// <item><description>lmom's PELPE3 uses the eight-digit constants 3π = 9.4247780 and √π = 1.7724539, which move
    /// its PE3 parameters and quantiles by up to 2E-8.</description></item>
    /// </list>
    /// </para>
    /// <para>
    /// H and Z depend on the random number streams, which differ, so they are checked within Monte Carlo tolerance.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_RegionalLmomRFAOracle
    {
        private static readonly int[] N = [98, 59, 90, 61, 65, 86, 78, 72, 67, 99, 49, 61, 69, 73, 70, 66, 59, 74, 82];

        private static readonly double[] Mean =
        [
            19.685, 62.580, 40.852, 46.045, 45.021, 31.042, 80.143, 41.305, 30.585, 32.932,
            17.560, 69.518, 47.653, 102.501, 52.413, 79.696, 44.643, 58.655, 39.024,
        ];

        private static readonly double[] LCV =
        [
            0.1209, 0.0915, 0.1124, 0.1032, 0.0967, 0.1328, 0.1008, 0.1143, 0.1107, 0.1179,
            0.1308, 0.1119, 0.1018, 0.1025, 0.1054, 0.1174, 0.1115, 0.1003, 0.1046,
        ];

        private static readonly double[] T3 =
        [
            0.0488, 0.0105, 0.0614, 0.0417, -0.0134, -0.0176, 0.0943, 0.0555, 0.0478, 0.0492,
            0.0940, -0.0429, 0.0435, 0.0182, -0.0224, 0.0124, -0.0346, 0.0446, 0.0128,
        ];

        private static readonly double[] T4 =
        [
            0.1433, 0.1569, 0.1541, 0.1429, 0.1568, 0.1206, 0.1967, 0.1210, 0.1371, 0.0900,
            0.1273, 0.0927, 0.1446, 0.1047, 0.1664, 0.1317, 0.1032, 0.1450, 0.1583,
        ];

        private static readonly double[] T5 =
        [
            -0.0004, 0.0020, -0.0058, -0.0022, 0.0173, 0.0235, 0.0856, 0.0487, 0.0316, 0.0225,
            0.0352, -0.0061, -0.0056, -0.0221, 0.0035, -0.0176, 0.0083, -0.0379, 0.0443,
        ];

        private static readonly double[] F = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.9, 0.95, 0.99, 0.999];

        private static readonly Lazy<RegionalHomogeneityTest> _test = new(() => new RegionalHomogeneityTest(Cascades(), 500, 12345));

        /// <summary>Returns the Cascades regional data.</summary>
        private static RegionalData Cascades()
        {
            var names = new string[N.Length];
            for (int i = 0; i < names.Length; i++) names[i] = (i + 1).ToString();
            return new RegionalData(names, N, Mean, LCV, T3, T4, T5);
        }

        /// <summary>
        /// D_i, the regional average L-moments and the observed V statistics match lmomRFA to rounding.
        /// </summary>
        [TestMethod]
        public void DiscordancyAveragesAndDispersion_MatchLmomRFA()
        {
            double[] d =
            [
                0.59749915864725933, 1.01792414800512110, 0.37895603475368045, 0.22849269297693847, 0.93079346376093008,
                2.63353718463120812, 2.12016691716794448, 0.45067562007606365, 0.11105346937902152, 1.61498444279377740,
                2.07758982464012609, 1.52107887661495744, 0.31444153391228996, 1.29739409494191915, 1.57707124936012222,
                0.28545596939310147, 1.03911449809153256, 0.42795699746744864, 0.37581382338655661,
            ];
            AssertClose(d, Cascades().Discordancy(), 1E-13);
            AssertClose([1d, 0.110298476052249628, 0.027859216255442670, 0.136613062409288832, 0.012227939042089986],
                Cascades().RegionalAverage(), 1E-15);
            var test = _test.Value;
            AssertClose([0.010438443057687518, 0.033922991896402549, 0.040468294227217538], [.. test.ObservedDispersion], 1E-15);
        }

        /// <summary>
        /// The kappa parent matches lmomRFA's rpara to the agreement of the two Newton iterates.
        /// </summary>
        [TestMethod]
        public void KappaParameters_MatchLmomRFA()
        {
            AssertClose([0.95416196659736985, 0.15327105718357470, 0.12359467971070576, -0.29549146940618615],
                [.. _test.Value.KappaParameters], 1E-11);
        }

        /// <summary>
        /// The candidates' τ4 match lmomRFA's t4fit; the GEV within lmom's shape approximation and the exact value.
        /// </summary>
        [TestMethod]
        public void CandidateTau4_MatchLmomRFA()
        {
            double[] r = [0.1673134466086395877, 0.1109982218873324378, 0.1232115397472264701, 0.1228309348996898392, 0.0063128052202939267];
            var tau4 = _test.Value.DistributionTau4;
            Assert.AreEqual(r[0], tau4[0], 1E-15, "GLO");
            Assert.AreEqual(r[1], tau4[1], 2E-8, "GEV against lmom");
            Assert.AreEqual(0.1109982343661195, tau4[1], 1E-13, "GEV against the exact shape (mpmath)");
            // lmom integrates the GNO L-moments numerically and agrees with mpmath quadrature (0.12321153974722603) to
            // 4E-16; Numerics' GeneralizedNormal.LinearMomentsFromParameters uses Hosking's rational approximation of
            // τ4, which is 1.5E-11 from it here.
            Assert.AreEqual(r[2], tau4[2], 2E-11, "GNO");
            Assert.AreEqual(r[3], tau4[3], 1E-11, "PE3");
            Assert.AreEqual(r[4], tau4[4], 1E-15, "GPA");
        }

        /// <summary>
        /// The growth curves match regfit and regquant.
        /// </summary>
        [TestMethod]
        public void GrowthCurves_MatchLmomRFA()
        {
            var data = Cascades();
            Check(UnivariateDistributionType.GeneralizedNormal, data, [0.994428518099721548, 0.195234211244999262, -0.057028459675727901],
                [0.56908430067329852, 0.61605675454949282, 0.68789834690449880, 0.75315010246273062, 0.83399614278899281,
                 0.99442851809972155, 1.25400111753404175, 1.33110405131376042, 1.48011669195602691, 1.65417491123852312], 1E-14);
            Check(UnivariateDistributionType.GeneralizedLogistic, data, [0.99494733433507043, 0.11015771274682390, -0.02785921625544267],
                [0.51982018103957073, 0.58865687628584218, 0.68354154082723417, 0.76016525938926971, 0.84514765017425775,
                 0.99494733433507043, 1.24455009527781346, 1.33297462242296283, 1.53496340358443484, 1.83390778341806593], 1E-14);
            Check(UnivariateDistributionType.GeneralizedPareto, data, [0.68106274295613911, 0.60329646341548082, 0.89158353278404956],
                [0.68709899010102815, 0.69314185137525364, 0.71131087668915594, 0.74173218816171504, 0.80313850668522790,
                 0.99298674742566884, 1.27086691838290600, 1.31090426287912765, 1.34657187413656976, 1.35628906598506282], 1E-14);

            // lmom's eight-digit constants in PELPE3.
            Check(UnivariateDistributionType.PearsonTypeIII, data, [1.0, 0.19567766921853608, 0.17099074266633735],
                [0.56951806585921938, 0.61628821001212408, 0.68792659045791149, 0.75308342841179599, 0.83389887651363814,
                 0.99442591026948657, 1.25407989914629736, 1.33108783514689488, 1.47965468202921890, 1.65263235514281082], 2E-8);

            // lmom's rational approximation of the GEV shape, then the exact shape (mpmath, 40 digits).
            Check(UnivariateDistributionType.GeneralizedExtremeValue, data, [0.92703876194538626, 0.18950153059337960, 0.23436548603407070],
                [0.57907063326659935, 0.62245166191219403, 0.68994078595869635, 0.75248377099772756, 0.83163654704531309,
                 0.99359411810614284, 1.25844654258915822, 1.33252166489173152, 1.46050602855695866, 1.57540738683961035], 2E-7);
            var gev = new RegionalGrowthCurve(UnivariateDistributionType.GeneralizedExtremeValue, data);
            AssertClose([0.92703874998556148, 0.1895015154841016, 0.23436534939759429], gev.Parameters, 1E-13);

            // sitequant(c(0.5, 0.99), regfit(Cascades, "gno")) for the first site.
            var gno = new RegionalGrowthCurve(UnivariateDistributionType.GeneralizedNormal, data);
            AssertClose([19.575325378793018, 29.136097081154389], gno.SiteQuantiles(Mean[0], [0.5, 0.99]), 1E-12);
        }

        /// <summary>
        /// H, Z and the simulated V statistics lie within Monte Carlo tolerance of lmomRFA's (set.seed(1), nsim = 500):
        /// H = 0.677, -1.375, -2.231 and Z = 3.609, -2.922, -1.506, -1.550, -15.063.
        /// </summary>
        [TestMethod]
        public void SimulatedStatistics_WithinMonteCarloToleranceOfLmomRFA()
        {
            var test = _test.Value;
            double[] h = [0.67740838504917988, -1.37476201503540807, -2.23083059771975600];
            double[] z = [3.6090552726796936, -2.9220214014689359, -1.5055995407116798, -1.5497396370020895, -15.0627610735764890];
            double[] vbar = [0.0093545868755127767, 0.0446053212122454260, 0.0574076934453318505];
            double[] vsd = [0.0016000040833507746, 0.0077703116605005599, 0.0075933149004809792];
            for (int j = 0; j < 3; j++)
            {
                Assert.AreEqual(h[j], test.Heterogeneity[j], 0.3 + 0.1 * Math.Abs(h[j]), $"H{j + 1}");
                Assert.AreEqual(vbar[j], test.SimulatedDispersionMean[j], 0.1 * vbar[j], $"mean V{j + 1}");
                Assert.AreEqual(vsd[j], test.SimulatedDispersionStandardDeviation[j], 0.2 * vsd[j], $"sd V{j + 1}");
            }
            for (int d = 0; d < z.Length; d++)
                Assert.AreEqual(z[d], test.Z[d], 0.3 + 0.1 * Math.Abs(z[d]), test.Distributions[d].ToString());
        }

        /// <summary>Fits a growth curve and checks its parameters and quantiles.</summary>
        private static void Check(UnivariateDistributionType type, RegionalData data, double[] parameters, double[] quantiles, double tolerance)
        {
            var curve = new RegionalGrowthCurve(type, data);
            AssertClose(parameters, curve.Parameters, tolerance);
            AssertClose(quantiles, curve.GrowthFactors(F), tolerance);
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
