using System;
using System.Globalization;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Integration;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Pert Percentile distribution algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil</item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests" />
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_PertPercentileDists
    {

        /// <summary>
        /// Verified using Palisade's @Risk
        /// </summary>
        [TestMethod()]
        public void Test_PertPercentile_AtRisk()
        {
            double true_min = -2.6453d;
            double true_max = 2.6453d;
            int true_mean = 0;
            int true_mode = 0;
            int true_median = 0;
            double true_stdDev = 0.9998d;
            int true_skew = 0;
            double true_kurt = 2.3333d;
            double true_icdf05 = -1.644;
            double true_icdf50 = 0.0;
            double true_icdf95 = 1.644;
            var pert = new PertPercentile(true_icdf05, true_icdf50, true_icdf95);
            pert.SolveParameters();
            Assert.AreEqual(pert.Minimum, true_min, 1E-3);
            Assert.AreEqual(pert.Maximum, true_max, 1E-3);
            Assert.AreEqual(pert.Mean, true_mean, 1E-3);
            Assert.AreEqual(pert.Median, true_median, 1E-3);
            Assert.AreEqual(pert.Mode, true_mode, 1E-3);
            Assert.AreEqual(pert.StandardDeviation, true_stdDev, 1E-3);
            Assert.AreEqual(pert.Skewness, true_skew, 1E-3);
            Assert.AreEqual(pert.Kurtosis, true_kurt, 1E-3);
            Assert.AreEqual(0.05, pert.CDF(true_icdf05),  1E-3);
            Assert.AreEqual(0.5, pert.CDF(true_icdf50), 1E-3);
            Assert.AreEqual(0.95, pert.CDF(true_icdf95), 1E-3);
            Assert.AreEqual(pert.InverseCDF(0.05d), true_icdf05, 1E-3);
            Assert.AreEqual(pert.InverseCDF(0.5d), true_icdf50, 1E-3);
            Assert.AreEqual(pert.InverseCDF(0.95d), true_icdf95, 1E-3);
        }

        /// <summary>
        /// Verified using Palisade's @Risk
        /// </summary>
        [TestMethod()]
        public void Test_PertPercentileZ_AtRisk()
        {
            double true_min = Normal.StandardCDF(-2.6453d);
            double true_max = Normal.StandardCDF(2.6453d);
            double true_mean = Normal.StandardCDF(0);
            // The density of X = Phi(Z) is g(z) / phi(z), which for this wide symmetric PERT is bimodal with nearly
            // equal peaks at z = -1.7312987 and z = 1.7312986. Global maximizer of
            // (alpha-1)*log(z-a) + (beta-1)*log(b-z) + z^2/2 by mpmath.findroot on its derivative (40 digits) for the
            // solved z-space PERT alpha = 3.0000000051165947, beta = 2.9999999829341211, a = -2.6452589821042789,
            // b = 2.6452589561908502; mode = mpmath.ncdf(-1.73129866479) = 0.041699254738463711.
            double true_mode = 0.041699254738463711;
            double true_median = Normal.StandardCDF(0);
            double true_icdf05 = Normal.StandardCDF(-1.644);
            double true_icdf50 = Normal.StandardCDF(0.0);
            double true_icdf95 = Normal.StandardCDF(1.644);
            var pert = new PertPercentileZ(true_icdf05, true_icdf50, true_icdf95);
            pert.SolveParameters();
            Assert.AreEqual(pert.Minimum, true_min, 1E-2);
            Assert.AreEqual(pert.Maximum, true_max, 1E-2);
            Assert.AreEqual(pert.Mean, true_mean, 1E-2);
            Assert.AreEqual(pert.Median, true_median, 1E-2);
            Assert.AreEqual(pert.Mode, true_mode, 1E-2);
            Assert.AreEqual(0.05, pert.CDF(true_icdf05),  1E-2);
            Assert.AreEqual(0.5, pert.CDF(true_icdf50), 1E-2);
            Assert.AreEqual(0.95, pert.CDF(true_icdf95), 1E-2);
            Assert.AreEqual(pert.InverseCDF(0.05d), true_icdf05, 1E-2);
            Assert.AreEqual(pert.InverseCDF(0.5d), true_icdf50, 1E-2);
            Assert.AreEqual(pert.InverseCDF(0.95d), true_icdf95, 1E-2);
        }

        /// <summary>
        /// Test several Pert Percentile configurations.
        /// </summary>
        [TestMethod()]
        public void Test_PertPercentile()
        {
            double fifth = 12;
            double fiftieth = 21;
            double ninetyfifth = 32;
            var pert = new PertPercentile(fifth, fiftieth, ninetyfifth);
            pert.SolveParameters();
            var p1 = pert.CDF(fifth);
            var p2 = pert.CDF(fiftieth);
            var p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(0.05, pert.CDF(fifth), 1E-1);
            Assert.AreEqual(0.5, pert.CDF(fiftieth), 1E-1);
            Assert.AreEqual(0.95, pert.CDF(ninetyfifth), 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

            fifth = 0;
            fiftieth = 3;
            ninetyfifth = 6;
            pert = new PertPercentile(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(0.05, pert.CDF(fifth), 1E-1);
            Assert.AreEqual(0.5, pert.CDF(fiftieth), 1E-1);
            Assert.AreEqual(0.95, pert.CDF(ninetyfifth), 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

            fifth = 0;
            fiftieth = 2;
            ninetyfifth = 3;
            pert = new PertPercentile(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(0.05, pert.CDF(fifth), 1E-1);
            Assert.AreEqual(0.5, pert.CDF(fiftieth), 1E-1);
            Assert.AreEqual(0.95, pert.CDF(ninetyfifth), 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);


            fifth = 14;
            fiftieth = 29;
            ninetyfifth = 36;
            pert = new PertPercentile(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(0.05, pert.CDF(fifth), 1E-1);
            Assert.AreEqual(0.5, pert.CDF(fiftieth), 1E-1);
            Assert.AreEqual(0.95, pert.CDF(ninetyfifth), 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

        }


        /// <summary>
        /// Test several Pert PercentileZ configurations.
        /// </summary>
        [TestMethod()]
        public void Test_PertPercentileZ()
        {
            double fifth = 0.12;
            double fiftieth = 0.21;
            double ninetyfifth = 0.32;
            var pert = new PertPercentileZ(fifth, fiftieth, ninetyfifth);
            pert.SolveParameters();
            var p1 = pert.CDF(fifth);
            var p2 = pert.CDF(fiftieth);
            var p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(0.05, pert.CDF(fifth), 1E-1);
            Assert.AreEqual(0.5, pert.CDF(fiftieth), 1E-1);
            Assert.AreEqual(0.95, pert.CDF(ninetyfifth), 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

            fifth = 0.01;
            fiftieth = 0.3;
            ninetyfifth = 0.6;
            pert = new PertPercentileZ(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(0.05, pert.CDF(fifth), 1E-1);
            Assert.AreEqual(0.5, pert.CDF(fiftieth), 1E-1);
            Assert.AreEqual(0.95, pert.CDF(ninetyfifth), 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

            fifth = 0.01;
            fiftieth = 0.02;
            ninetyfifth = 0.03;
            pert = new PertPercentileZ(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(0.05, pert.CDF(fifth), 1E-1);
            Assert.AreEqual(0.5, pert.CDF(fiftieth), 1E-1);
            Assert.AreEqual(0.95, pert.CDF(ninetyfifth), 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);


            fifth = 0.014;
            fiftieth = 0.029;
            ninetyfifth = 0.036;
            pert = new PertPercentileZ(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(0.05, pert.CDF(fifth), 1E-1);
            Assert.AreEqual(0.5, pert.CDF(fiftieth), 1E-1);
            Assert.AreEqual(0.95, pert.CDF(ninetyfifth), 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

        }

        /// <summary>
        /// Rebuilds the solved z-space Beta exactly from the G17 serialized parameters.
        /// </summary>
        private static GeneralizedBeta SolvedBeta(UnivariateDistributionBase distribution)
        {
            var parts = distribution.ToXElement().Attribute("BetaParameters")!.Value.Split('|');
            var values = new double[parts.Length];
            for (int i = 0; i < parts.Length; i++)
                values[i] = double.Parse(parts[i], CultureInfo.InvariantCulture);
            return new GeneralizedBeta(values[0], values[1], values[2], values[3]);
        }

        /// <summary>
        /// The density of X = Phi(Z) includes the Jacobian 1 / phi(z) and integrates to one.
        /// </summary>
        [TestMethod]
        public void Test_PertPercentileZ_PDF_Jacobian()
        {
            var dist = new PertPercentileZ(0.05, 0.2, 0.6);
            // mpmath (40 digits): g(z) / npdf(z) at z = sqrt(2) * erfinv(2 * 0.3 - 1) for the solved PERT
            // alpha = 2.1453766604955331, beta = 3.8546233395044638, a = -1.9342278025553861, b = 1.2714160325182078
            Assert.AreEqual(1.4974442027737657, dist.PDF(0.3), 1E-12 * 1.4974442027737657);

            var beta = SolvedBeta(dist);
            double x = 0.3, z = Normal.StandardZ(x);
            Assert.AreEqual(beta.PDF(z) / Normal.StandardPDF(z), dist.PDF(x), 0d);

            // Central difference of the CDF agrees with the density
            double h = 1E-6;
            Assert.AreEqual((dist.CDF(x + h) - dist.CDF(x - h)) / (2 * h), dist.PDF(x), 1E-6);

            // mpmath.quad of the density over the support: 1
            var agk = new AdaptiveGaussKronrod(dist.PDF, dist.Minimum, dist.Maximum) { RelativeTolerance = 1E-10 };
            agk.Integrate();
            Assert.AreEqual(1d, agk.Result, 1E-8);
        }

        /// <summary>
        /// Moments, mode and median of PertPercentileZ are those of X = Phi(Z).
        /// </summary>
        [TestMethod]
        public void Test_PertPercentileZ_Moments()
        {
            // mpmath.quad (40 digits) of h(ncdf(z)) * g(z) over [a, b] for the solved PERT
            // alpha = 2.1453766604955331, beta = 3.8546233395044638, a = -1.9342278025553861, b = 1.2714160325182078
            var dist = new PertPercentileZ(0.05, 0.2, 0.6);
            Assert.AreEqual(0.2471240592128751, dist.Mean, 1E-10 * 0.2471240592128751);
            Assert.AreEqual(0.17330030019022896, dist.StandardDeviation, 1E-10 * 0.17330030019022896);
            Assert.AreEqual(0.96965417127997913, dist.Skewness, 1E-9);
            Assert.AreEqual(3.2413257093916389, dist.Kurtosis, 1E-9);
            // mpmath.ncdf of the root of (alpha-1)/(z-a) - (beta-1)/(b-z) + z = 0
            Assert.AreEqual(0.06954738548994736, dist.Mode, 1E-10);
            Assert.AreEqual(dist.PDF(dist.Mode), Math.Max(dist.PDF(dist.Mode - 1E-4), Math.Max(dist.PDF(dist.Mode), dist.PDF(dist.Mode + 1E-4))), 0d);

            // The median is Phi of the z-space median
            var beta = SolvedBeta(dist);
            Assert.AreEqual(Normal.StandardCDF(beta.Median), dist.Median, 0d);

            // Narrow left-tail case: alpha = 4.1907820674280813, beta = 1.8092179325719133,
            // a = -2.6245808826022099, b = -1.8335358126092394 (mpmath, 40 digits)
            var narrow = new PertPercentileZ(0.01, 0.02, 0.03);
            Assert.AreEqual(0.02000663990831981, narrow.Mean, 1E-10 * 0.02);
            Assert.AreEqual(0.0061338232454936886, narrow.StandardDeviation, 1E-9 * 0.0061338232454936886);
            Assert.AreEqual(-0.016669330257079791, narrow.Skewness, 1E-9);
            Assert.AreEqual(2.1669057631789149, narrow.Kurtosis, 1E-9);
            Assert.AreEqual(0.019734574512671628, narrow.Mode, 1E-12);

            // Base-class adaptive integration of x f(x) agrees
            var moments = dist.CentralMoments(1E-10);
            Assert.AreEqual(dist.Mean, moments[0], 1E-8);
            Assert.AreEqual(dist.StandardDeviation, moments[1], 1E-8);
        }

        /// <summary>
        /// CDF and InverseCDF of PertPercentileZ inside (0, 1) are Phi-transforms of the z-space PERT.
        /// </summary>
        [TestMethod]
        public void Test_PertPercentileZ_CDF_InverseCDF_MatchBetaExpressions()
        {
            var dist = new PertPercentileZ(0.05, 0.2, 0.6);
            dist.SolveParameters();
            var beta = SolvedBeta(dist);
            foreach (double x in new[] { 1E-12, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.8, 0.9, 0.99, 1 - 1E-12 })
                Assert.AreEqual(beta.CDF(Normal.StandardZ(x)), dist.CDF(x), 0d);
            foreach (double p in new[] { 0d, 1E-10, 0.05, 0.25, 0.5, 0.75, 0.95, 1 - 1E-10, 1d })
                Assert.AreEqual(Normal.StandardCDF(beta.InverseCDF(p)), dist.InverseCDF(p), 0d);
        }

        /// <summary>
        /// PertPercentileZ returns zero density and a 0/1 CDF outside its support, and requires percentiles in (0, 1).
        /// </summary>
        [TestMethod]
        public void Test_PertPercentileZ_SupportAndValidation()
        {
            var dist = new PertPercentileZ(0.05, 0.2, 0.6);
            foreach (double x in new[] { double.NegativeInfinity, -0.5, 0d, 1d, 1.5, double.PositiveInfinity })
                Assert.AreEqual(0d, dist.PDF(x));
            Assert.AreEqual(0d, dist.CDF(-0.5));
            Assert.AreEqual(0d, dist.CDF(0d));
            Assert.AreEqual(1d, dist.CDF(1d));
            Assert.AreEqual(1d, dist.CDF(1.5));
            Assert.AreEqual(0d, dist.PDF(dist.Minimum / 2));
            Assert.AreEqual(0d, dist.CDF(dist.Minimum / 2));

            Assert.IsFalse(new PertPercentileZ(0d, 0.5, 0.9).ParametersValid);
            Assert.IsFalse(new PertPercentileZ(0.1, 0.5, 1d).ParametersValid);
            Assert.IsFalse(new PertPercentileZ(0d, 0d, 0d).ParametersValid);
            Assert.IsFalse(new PertPercentileZ(1d, 1d, 1d).ParametersValid);
            Assert.IsTrue(new PertPercentileZ(1E-10, 0.5, 1 - 1E-10).ParametersValid);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new PertPercentileZ(0d, 0.5, 0.9).InverseCDF(0.5));
        }

        /// <summary>
        /// Summary properties solve the parameters on first use and match the values after an explicit evaluation.
        /// </summary>
        [TestMethod]
        public void Test_PertPercentile_PropertiesSolveParameters()
        {
            var unsolved = new PertPercentile(1, 5, 10);
            var solved = new PertPercentile(1, 5, 10);
            solved.InverseCDF(0.5);
            // Regression: solved PERT mean, (a + 4c + b) / 6 for the Nelder-Mead fit of the three percentiles
            Assert.AreEqual(5.183754401166906, unsolved.Mean, 1E-12);
            Assert.AreEqual(solved.Mean, unsolved.Mean, 0d);
            Assert.AreEqual(solved.Median, new PertPercentile(1, 5, 10).Median, 0d);
            Assert.AreEqual(solved.Mode, new PertPercentile(1, 5, 10).Mode, 0d);
            Assert.AreEqual(solved.StandardDeviation, new PertPercentile(1, 5, 10).StandardDeviation, 0d);
            Assert.AreEqual(solved.Skewness, new PertPercentile(1, 5, 10).Skewness, 0d);
            Assert.AreEqual(solved.Kurtosis, new PertPercentile(1, 5, 10).Kurtosis, 0d);
            Assert.AreEqual(solved.Minimum, new PertPercentile(1, 5, 10).Minimum, 0d);
            Assert.AreEqual(solved.Maximum, new PertPercentile(1, 5, 10).Maximum, 0d);
            var pert = new PertPercentile(1, 5, 10).ToPert();
            var solvedPert = solved.ToPert();
            Assert.AreEqual(solvedPert.Min, pert.Min, 0d);
            Assert.AreEqual(solvedPert.MostLikely, pert.MostLikely, 0d);
            Assert.AreEqual(solvedPert.Max, pert.Max, 0d);
            Assert.AreNotEqual(0d, pert.Min);

            var z = new PertPercentileZ(0.1, 0.1, 0.1);
            Assert.AreEqual(0.1, z.Mean, 1E-15);
            var zSolved = new PertPercentileZ(0.05, 0.2, 0.6);
            zSolved.CDF(0.3);
            var zUnsolved = new PertPercentileZ(0.05, 0.2, 0.6);
            Assert.AreEqual(zSolved.Median, zUnsolved.Median, 0d);
            Assert.AreEqual(zSolved.Minimum, new PertPercentileZ(0.05, 0.2, 0.6).Minimum, 0d);
            Assert.AreEqual(zSolved.Maximum, new PertPercentileZ(0.05, 0.2, 0.6).Maximum, 0d);
            Assert.AreEqual(zSolved.Mean, new PertPercentileZ(0.05, 0.2, 0.6).Mean, 0d);
        }

        /// <summary>
        /// Equal percentiles give the fixed-value PERT.
        /// </summary>
        [TestMethod]
        public void Test_PertPercentile_EqualPercentiles_ToPert()
        {
            var pert = new PertPercentile(3, 3, 3).ToPert();
            Assert.IsTrue(pert.ParametersValid);
            Assert.AreEqual(3d, pert.Min);
            Assert.AreEqual(3d, pert.MostLikely);
            Assert.AreEqual(3d, pert.Max);

            var zPert = new PertPercentileZ(0.1, 0.1, 0.1).ToPert();
            double z = Normal.StandardZ(0.1);
            Assert.IsTrue(zPert.ParametersValid);
            Assert.AreEqual(z, zPert.Min);
            Assert.AreEqual(z, zPert.MostLikely);
            Assert.AreEqual(z, zPert.Max);
        }

        /// <summary>
        /// Allowable bounds censor the PERT: the CDF is zero below the minimum allowable value and carries the point mass
        /// at it, the density is zero outside the bounds, and the moments are those of the censored variable.
        /// </summary>
        [TestMethod]
        public void Test_AllowableBounds_Censoring()
        {
            var dist = new PertPercentile(1d, 5d, 10d) { MinAllowableValue = 3d };
            // Solved Beta: alpha 2.3697534643459002, beta 3.6302465356541025 on [-0.68495321122778774, 14.174080103483483].
            // mpmath (40 digits): F(3) = 0.24192706859137055; moments of max(X, 3) by mpmath.quad plus the point mass at 3.
            Assert.AreEqual(0d, dist.CDF(2.999999));
            Assert.AreEqual(0d, dist.CDF(-100d));
            Assert.AreEqual(0.24192706859137055, dist.CDF(3d), 1E-14);
            Assert.AreEqual(0d, dist.PDF(2.9d));
            Assert.AreEqual(0d, dist.PDF(-100d));
            Assert.AreEqual(3d, dist.InverseCDF(0.01));
            Assert.AreEqual(3d, dist.Minimum);
            Assert.AreEqual(5.4829387426215514, dist.Mean, 1E-12);
            Assert.AreEqual(2.3660486241767284, dist.StandardDeviation, 1E-12);
            Assert.AreEqual(0.74304678881611099, dist.Skewness, 1E-11);
            Assert.AreEqual(2.6427506319325361, dist.Kurtosis, 1E-11);

            var unit = new PertPercentile(0.05, 0.5, 0.95) { MinAllowableValue = 0d, MaxAllowableValue = 1d };
            // Solved Beta: alpha 3.0000000692627484, beta 2.9999999665728749 on [-0.22406724551791274, 1.2240672335723484].
            Assert.AreEqual(0.028977858026146096, unit.CDF(0d), 1E-14);
            Assert.AreEqual(1d, unit.CDF(1d));
            Assert.AreEqual(0.50000000634771135, unit.Mean, 1E-12);
            Assert.AreEqual(0.26675828116536562, unit.StandardDeviation, 1E-12);
            Assert.AreEqual(2.137979762618169, unit.Kurtosis, 1E-11);
        }

        /// <summary>
        /// Bounds outside the solved PERT support do not censor, so every value equals the PERT value bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_AllowableBounds_NonCensoringMatchPert()
        {
            var bounded = new PertPercentile(1d, 5d, 10d) { MinAllowableValue = -100d, MaxAllowableValue = 100d };
            var plain = new PertPercentile(1d, 5d, 10d);
            var pert = plain.ToPert();
            Assert.AreEqual(plain.Mean, bounded.Mean, 0d);
            Assert.AreEqual(plain.StandardDeviation, bounded.StandardDeviation, 0d);
            Assert.AreEqual(plain.Skewness, bounded.Skewness, 0d);
            Assert.AreEqual(plain.Kurtosis, bounded.Kurtosis, 0d);
            for (int i = -2; i <= 16; i++)
            {
                Assert.AreEqual(plain.PDF(i), bounded.PDF(i), 0d);
                Assert.AreEqual(plain.CDF(i), bounded.CDF(i), 0d);
            }
            for (int i = 1; i < 20; i++)
                Assert.AreEqual(plain.InverseCDF(i / 20d), bounded.InverseCDF(i / 20d), 0d);
            Assert.AreEqual(pert.Min, plain.Minimum, 0d);
        }

        /// <summary>
        /// Finite allowable bounds survive an XML round trip; infinite bounds add no attributes, and invalid percentiles
        /// read from XML are reported as invalid.
        /// </summary>
        [TestMethod]
        public void Test_AllowableBounds_XElementRoundTrip()
        {
            var dist = new PertPercentile(0.01, 0.1, 0.3) { MinAllowableValue = 0d, MaxAllowableValue = 1d };
            double before = dist.InverseCDF(0.999);
            var copy = PertPercentile.FromXElement(dist.ToXElement())!;
            Assert.AreEqual(0d, copy.MinAllowableValue);
            Assert.AreEqual(1d, copy.MaxAllowableValue);
            Assert.AreEqual(before, copy.InverseCDF(0.999), 0d);
            Assert.AreEqual(dist.Mean, copy.Mean, 0d);

            var unbounded = new PertPercentile(0.01, 0.1, 0.3).ToXElement();
            Assert.IsNull(unbounded.Attribute(nameof(PertPercentile.MinAllowableValue)));
            Assert.IsNull(unbounded.Attribute(nameof(PertPercentile.MaxAllowableValue)));
            var unboundedCopy = PertPercentile.FromXElement(unbounded)!;
            Assert.AreEqual(double.NegativeInfinity, unboundedCopy.MinAllowableValue);
            Assert.AreEqual(double.PositiveInfinity, unboundedCopy.MaxAllowableValue);

            var reversed = new System.Xml.Linq.XElement("Distribution",
                new System.Xml.Linq.XAttribute("Type", "PertPercentile"),
                new System.Xml.Linq.XAttribute("ParametersSolved", "False"),
                new System.Xml.Linq.XAttribute("Parameters", "10|5|1"),
                new System.Xml.Linq.XAttribute("BetaParameters", "2|2|0|1"));
            Assert.IsFalse(PertPercentile.FromXElement(reversed)!.ParametersValid);
        }
    }
}
