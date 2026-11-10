using System;
using System.Globalization;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Tests the Gamma distribution tails against the mpmath oracle of the incomplete gamma functions and R's pgamma.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// The oracle rows are docs/mathematics/oracles/special-functions/special-functions-mpmath.csv (gammap and gammaq,
    /// 50-digit references). The tolerance is that of the incomplete gamma oracle test: 5E-14 relative plus 16 |ln v| units
    /// of roundoff, the conditioning of the factor x^a e^(-x) / Γ(a).
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_GammaDistributionTails
    {
        /// <summary>The oracle rows: function, shape, argument and reference value.</summary>
        private static readonly (string Function, double A, double X, double Value)[] Rows =
            DistributionOracle.Read("special-functions-mpmath.csv")
                .Where(r => r["function"] == "gammap" || r["function"] == "gammaq")
                .Select(r => (r["function"], Parse(r["a"]), Parse(r["x"]), Parse(r["value"])))
                .ToArray();

        /// <summary>Parses an oracle number written by Python repr.</summary>
        private static double Parse(string text) => double.Parse(text, CultureInfo.InvariantCulture);

        /// <summary>The incomplete gamma tolerance: 5E-14 relative plus 16 |ln v| units of roundoff, and two subnormal steps.</summary>
        private static double Tolerance(double v) => v == 0 ? 0 :
            Math.Max((5E-14 + 16 * 1.1102230246251565E-16 * Math.Abs(Math.Log(v))) * v, 2 * double.Epsilon);

        /// <summary>
        /// CDF and CCDF match mpmath for shapes from 1E-3 to 1E6, with the scale 1 and the exact scale 0.5.
        /// </summary>
        [TestMethod]
        public void Test_CDF_CCDF_MatchMpmath()
        {
            int count = 0;
            foreach (var row in Rows)
            {
                foreach (double scale in new[] { 1.0, 0.5 })
                {
                    var gamma = new GammaDistribution(scale, row.A);
                    double actual = row.Function == "gammap" ? gamma.CDF(row.X * scale) : gamma.CCDF(row.X * scale);
                    Assert.AreEqual(row.Value, actual, Tolerance(row.Value), $"{row.Function}(a={row.A:R}, x={row.X:R}), scale {scale}");
                }
                count++;
            }
            Assert.IsGreaterThan(1400, count);
        }

        /// <summary>
        /// LogCDF and LogCCDF match the logarithm of the mpmath value to 2E-15 absolute plus 4E-15 relative and the
        /// conditioning allowance, which covers rows whose probability is far below 1E-300.
        /// </summary>
        [TestMethod]
        public void Test_LogCDF_LogCCDF_MatchMpmath()
        {
            foreach (var row in Rows.Where(r => r.Value > 0))
            {
                var gamma = new GammaDistribution(1.0, row.A);
                double expected = Math.Log(row.Value);
                double actual = row.Function == "gammap" ? gamma.LogCDF(row.X) : gamma.LogCCDF(row.X);
                Assert.AreEqual(expected, actual, 2E-15 + 4E-15 * Math.Abs(expected) + 2 * Tolerance(row.Value) / row.Value,
                    $"{row.Function}(a={row.A:R}, x={row.X:R})");
            }
        }

        /// <summary>
        /// At a = 21544 and x = 1.25 a the upper tail, 5.64E-254, is within 2E-13 relative of mpmath.
        /// </summary>
        [TestMethod]
        public void Test_LargeShapeUpperTail()
        {
            var row = Rows.Single(r => r.Function == "gammaq" && Math.Abs(r.A - 21544.346900318822) < 1E-6 && Math.Abs(r.X / r.A - 1.25) < 1E-12);
            var gamma = new GammaDistribution(1.0, row.A);
            Assert.AreEqual(row.Value, gamma.CCDF(row.X), 2E-13 * row.Value);
            Assert.AreEqual(Math.Log(row.Value), gamma.LogCCDF(row.X), 2E-13);
        }

        /// <summary>
        /// Shapes from 1E8 to 4E12, where the power series would need of order sqrt(a) terms, match R 4.4.3 pgamma,
        /// including log tails far below the double range.
        /// </summary>
        [TestMethod]
        public void Test_HugeShapes_MatchR()
        {
            // R: pgamma(x, a), pgamma(x, a, log.p = TRUE) and pgamma(x, a, lower.tail = FALSE, log.p = TRUE).
            var cases = new (double A, double X, double P, double LogP, double LogQ)[]
            {
                (1E8, 99970000, 0.0013487164491615515, -6.6086019172071317, -0.0013496267858077394),
                (1E8, 95000000, 0, -129336.57230548575, 0),
                (1E8, 105000000, 1, 0, -120990.7166076321),
                (1E10, 9999700000, 0.0013497798514433145, -6.6078137728345991, -0.0013506916248215723),
                (1E10, 10000300000, 0.99864998378330805, -0.0013509283095704428, -6.6076386742395705),
                (4E12, 3999994000000, 0.0013498921225049716, -6.6077305989802904, -0.0013508040476355896),
                (4E12, 3800000000000, 0, -5173177562.633997, 0),
                (4E12, 4200000000000, 1, 0, -4839343334.7038517),
            };
            foreach (var c in cases)
            {
                var gamma = new GammaDistribution(1.0, c.A);
                Assert.AreEqual(c.P, gamma.CDF(c.X), 1E-14 * c.P, $"CDF a={c.A}");
                Assert.AreEqual(c.P, Gamma.LowerIncomplete(c.A, c.X), 1E-14 * c.P, $"P a={c.A}");
                // The logarithm is conditioned by |ln P| units of roundoff in the argument.
                Assert.AreEqual(c.LogP, gamma.LogCDF(c.X), 1E-15 + 1E-14 * Math.Abs(c.LogP), $"LogCDF a={c.A}");
                Assert.AreEqual(c.LogQ, gamma.LogCCDF(c.X), 1E-15 + 1E-14 * Math.Abs(c.LogQ), $"LogCCDF a={c.A}");
            }
        }

        /// <summary>
        /// The edges: zero below the support, one at +∞, NaN for NaN, and the log path when x / θ underflows.
        /// </summary>
        [TestMethod]
        public void Test_Edges()
        {
            var gamma = new GammaDistribution(2.0, 3.5);
            Assert.AreEqual(0.0, gamma.CDF(-1.0), 0.0);
            Assert.AreEqual(0.0, gamma.CDF(0.0), 0.0);
            Assert.AreEqual(1.0, gamma.CCDF(0.0), 0.0);
            Assert.AreEqual(1.0, gamma.CDF(double.PositiveInfinity), 0.0);
            Assert.AreEqual(0.0, gamma.CCDF(double.PositiveInfinity), 0.0);
            Assert.IsTrue(double.IsNaN(gamma.CDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(gamma.CCDF(double.NaN)));
            Assert.AreEqual(double.NegativeInfinity, gamma.LogCDF(0.0));
            Assert.AreEqual(0.0, gamma.LogCCDF(0.0), 0.0);
            // x / θ underflows to zero: the CDF is the exponential of the log form, (x/θ)^κ / Γ(κ + 1).
            var wide = new GammaDistribution(1E300, 0.5);
            double x = 1E-30;
            double expected = Math.Exp(0.5 * (Math.Log(x) - Math.Log(1E300)) - Gamma.LogGamma(1.5));
            Assert.AreEqual(expected, wide.CDF(x), 1E-14 * expected);
            Assert.AreEqual(1.0, wide.CCDF(x), 0.0);
            // CDF and CCDF are the incomplete gamma functions of x / θ.
            Assert.AreEqual(Gamma.LowerIncomplete(3.5, 2.5), gamma.CDF(5.0), 0.0);
            Assert.AreEqual(Gamma.UpperIncomplete(3.5, 2.5), gamma.CCDF(5.0), 0.0);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new GammaDistribution(2.0, -1.0).CDF(1.0));
        }
    }
}
