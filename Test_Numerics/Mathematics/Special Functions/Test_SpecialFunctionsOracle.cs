using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using Distributions.Univariate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.SpecialFunctions;

namespace Mathematics.SpecialFunctions
{
    /// <summary>
    /// Compares the Bessel, Debye, incomplete gamma and complementary error functions and the standard Normal CDF with the mpmath oracle
    /// docs/mathematics/oracles/special-functions/special-functions-mpmath.csv (50-digit references).
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// The tolerances are relative except for the oscillatory J and Y functions above 8, where the error is
    /// measured against the larger of |f| and the envelope sqrt(2/(πx)), and for the incomplete gamma tails, where
    /// the allowance grows with |ln P|, the conditioning of the factor x^a e^(-x) / Γ(a). Each tolerance is about
    /// twice the worst error measured on the oracle grid.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_SpecialFunctionsOracle
    {
        /// <summary>The oracle rows: function, integer order, shape, argument and reference value.</summary>
        private static readonly List<(string Function, int N, double A, double X, double Value)> Rows =
            DistributionOracle.Read("special-functions-mpmath.csv")
                .Select(r => (r["function"], int.Parse(r["n"], CultureInfo.InvariantCulture), Parse(r["a"]), Parse(r["x"]), Parse(r["value"])))
                .ToList();

        /// <summary>Parses an oracle number written by Python repr.</summary>
        private static double Parse(string text) => double.Parse(text, CultureInfo.InvariantCulture);

        /// <summary>Asserts every oracle row of a function and returns the number of rows checked.</summary>
        private static int Check(string function, Func<int, double, double, double> evaluate, Func<double, double, double, double, double> tolerance)
        {
            int count = 0;
            foreach (var row in Rows.Where(r => r.Function == function))
            {
                double actual = evaluate(row.N, row.A, row.X);
                double allowed = tolerance(row.A, row.X, row.Value, actual);
                Assert.IsLessThanOrEqualTo(allowed, Math.Abs(actual - row.Value),
                    $"{function}(n={row.N}, a={row.A:R}, x={row.X:R}) = {actual:R}, expected {row.Value:R}");
                count++;
            }
            return count;
        }

        /// <summary>A relative tolerance.</summary>
        private static Func<double, double, double, double, double> Relative(double rel) => (a, x, v, f) => rel * Math.Abs(v);

        /// <summary>A tolerance relative to the larger of |f| and the envelope sqrt(2/(πx)) of J and Y.</summary>
        private static Func<double, double, double, double, double> Envelope(double rel) =>
            (a, x, v, f) => rel * Math.Max(Math.Abs(v), Math.Abs(x) > 1 ? Math.Sqrt(2 / (Math.PI * Math.Abs(x))) : 1.0);

        /// <summary>Relative for |x| &lt;= 8, where the zeros are factored out, and envelope-relative above.</summary>
        private static Func<double, double, double, double, double> Oscillatory(double rel, double envelope) =>
            (a, x, v, f) => Math.Abs(x) <= 8 ? rel * Math.Abs(v) : Envelope(envelope)(a, x, v, f);

        /// <summary>The incomplete gamma tolerance: 5E-14 relative plus 16 |ln v| units of roundoff, and two subnormal steps.</summary>
        private static double GammaTolerance(double v) => v == 0 ? 0 :
            Math.Max((5E-14 + 16 * 1.1102230246251565E-16 * Math.Abs(Math.Log(v))) * v, 2 * double.Epsilon);

        /// <summary>
        /// I0, I1, In, K0, K1 and Kn match mpmath over |x| from 1E-8 to 700 (orders 2 to 50 from 1E-3 to 600).
        /// </summary>
        [TestMethod]
        public void Test_ModifiedBessel_MatchMpmath()
        {
            Assert.IsGreaterThan(200, Check("i0", (n, a, x) => Bessel.I0(x), Relative(2E-15)));
            Assert.IsGreaterThan(200, Check("i1", (n, a, x) => Bessel.I1(x), Relative(2E-15)));
            Assert.IsGreaterThan(300, Check("in", (n, a, x) => Bessel.In(n, x), Relative(1E-14)));
            Assert.IsGreaterThan(200, Check("k0", (n, a, x) => Bessel.K0(x), Relative(2E-15)));
            Assert.IsGreaterThan(200, Check("k1", (n, a, x) => Bessel.K1(x), Relative(2E-15)));
            Assert.IsGreaterThan(300, Check("kn", (n, a, x) => Bessel.Kn(n, x), Relative(1E-14)));
        }

        /// <summary>
        /// J0, J1, Y0 and Y1 match mpmath to 3E-15 relative for x &lt;= 8, including within 1E-6 of the first zeros,
        /// and to 2E-15 of the envelope above 8, from 1E-8 to 1E4.
        /// </summary>
        [TestMethod]
        public void Test_Bessel_J_Y_MatchMpmath()
        {
            Assert.IsGreaterThan(600, Check("j0", (n, a, x) => Bessel.J0(x), Oscillatory(3E-15, 2E-15)));
            Assert.IsGreaterThan(600, Check("j1", (n, a, x) => Bessel.J1(x), Oscillatory(3E-15, 2E-15)));
            Assert.IsGreaterThan(600, Check("y0", (n, a, x) => Bessel.Y0(x), Oscillatory(3E-15, 2E-15)));
            Assert.IsGreaterThan(600, Check("y1", (n, a, x) => Bessel.Y1(x), Oscillatory(3E-15, 2E-15)));
        }

        /// <summary>
        /// Jn and Yn of orders 2 to 50 match mpmath to 1E-14 of the larger of |f| and the envelope, for x from 0.01 to 1000.
        /// </summary>
        [TestMethod]
        public void Test_Bessel_Jn_Yn_MatchMpmath()
        {
            Assert.IsGreaterThan(800, Check("jn", (n, a, x) => Bessel.Jn(n, x), Envelope(1E-14)));
            Assert.IsGreaterThan(800, Check("yn", (n, a, x) => Bessel.Yn(n, x), Envelope(1E-14)));
        }

        /// <summary>
        /// Near each of the first zeros of J0, J1, Y0 and Y1 the relative error stays small, because the zeros are factored out.
        /// </summary>
        [TestMethod]
        public void Test_Bessel_NearZeros_Relative()
        {
            // Y0(y₀,₁ + 1E-6), J0(j₀,₂ - 1E-6) and Y1(y₁,₂ + 1E-6) from mpmath at 50 digits.
            foreach (var row in Rows.Where(r => (r.Function == "j0" || r.Function == "y0" || r.Function == "j1" || r.Function == "y1")
                && r.X > 0 && r.X <= 8 && Math.Abs(r.Value) < 1E-5))
            {
                double actual = row.Function == "j0" ? Bessel.J0(row.X) : row.Function == "j1" ? Bessel.J1(row.X)
                    : row.Function == "y0" ? Bessel.Y0(row.X) : Bessel.Y1(row.X);
                Assert.AreEqual(row.Value, actual, 3E-15 * Math.Abs(row.Value), $"{row.Function}({row.X:R})");
            }
        }

        /// <summary>
        /// The order-3 Debye function matches mpmath to 1E-15 relative over x from 1E-6 to 1000, with its limits.
        /// </summary>
        [TestMethod]
        public void Test_Debye3_MatchMpmath()
        {
            Assert.IsGreaterThan(250, Check("debye3", (n, a, x) => Debye.Function(x), Relative(1E-15)));
            Assert.AreEqual(1.0, Debye.Function(0.0), 0.0);
            Assert.AreEqual(0.0, Debye.Function(double.PositiveInfinity), 0.0);
            Assert.IsTrue(double.IsNaN(Debye.Function(double.NaN)));
            // Far above the series range D3(x) = (π⁴/5) / x³ exactly in double: 3 * π⁴/15 / 1E3 / 1E3 / 1E3.
            Assert.AreEqual(3.0 * 6.493939402266829 / 1E3 / 1E3 / 1E3, Debye.Function(1000.0), 0.0);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Debye.Function(-0.5));
        }

        /// <summary>
        /// P and Q match mpmath over a from 1E-3 to 1E6 and x from 1E-4 a to 30 a, and AS 239 (Gamma.Incomplete)
        /// returns the same P.
        /// </summary>
        [TestMethod]
        public void Test_IncompleteGamma_MatchMpmath()
        {
            Assert.IsGreaterThan(700, Check("gammap", (n, a, x) => Gamma.LowerIncomplete(a, x), (a, x, v, f) => GammaTolerance(v)));
            Assert.IsGreaterThan(700, Check("gammaq", (n, a, x) => Gamma.UpperIncomplete(a, x), (a, x, v, f) => GammaTolerance(v)));
            Assert.IsGreaterThan(700, Check("gammap", (n, a, x) => Gamma.Incomplete(x, a), (a, x, v, f) => GammaTolerance(v)));
        }

        /// <summary>
        /// The incomplete gamma integral underflows to zero only where the value is below the double range, far below
        /// the exp(-180) floor of algorithm AS 239.
        /// </summary>
        [TestMethod]
        public void Test_IncompleteGamma_NoPrematureUnderflow()
        {
            // mpmath: P(46.41588833612782, 0.004641588833612782) = 1.8173331985826257E-167.
            Assert.AreEqual(1.8173331985826257E-167, Gamma.Incomplete(0.004641588833612782, 46.41588833612782), 1E-13 * 1.8173331985826257E-167);
            // mpmath: P(2, 1E-160) = 5E-321, a subnormal.
            Assert.AreEqual(5E-321, Gamma.LowerIncomplete(2.0, 1E-160), 2 * double.Epsilon);
            // mpmath: Q(1.5, 700) = 2.945619361016309E-303.
            Assert.AreEqual(2.945619361016309E-303, Gamma.UpperIncomplete(1.5, 700.0), 1E-13 * 2.945619361016309E-303);
        }

        /// <summary>
        /// The inverse lower and upper incomplete gamma functions return the mpmath root of P(a, x) = p or Q(a, x) = q
        /// for a from 0.01 to 1E5 and probabilities from 1E-300 to 1 - 1E-10. On these 375 roots the largest relative error
        /// is 1.8E-15 for a ≥ 1 and 2.8E-14 for 0.01 ≤ a &lt; 1, where the condition number near 1/a of a small quantile
        /// amplifies the error of P.
        /// </summary>
        [TestMethod]
        public void Test_InverseIncompleteGamma_MatchMpmath()
        {
            Func<double, double, double, double, double> tolerance = (a, p, v, f) => (a >= 1d ? 2.5E-15 : 3.5E-14) * v;
            Assert.IsGreaterThan(150, Check("gammapinv", (n, a, p) => Gamma.InverseLowerIncomplete(a, p), tolerance));
            Assert.IsGreaterThan(150, Check("gammaqinv", (n, a, q) => Gamma.InverseUpperIncomplete(a, q), tolerance));
            // An upper probability near one is inverted through the exact complement 1 - q; inverting Q directly errs by 2.6E-5.
            Assert.AreEqual(6.073053387272577E-101, Gamma.InverseUpperIncomplete(0.1, 0.9999999999), 1E-13 * 6.073053387272577E-101);
        }

        /// <summary>A relative tolerance for normal results, and two subnormal steps once the result is subnormal.</summary>
        private static Func<double, double, double, double, double> RelativeOrSubnormal(double rel) =>
            (a, x, v, f) => Math.Max(rel * Math.Abs(v), Math.Abs(v) < 2.2250738585072014E-308 ? 2 * double.Epsilon : 0);

        /// <summary>
        /// Erfc matches mpmath to 1E-15 relative (the worst measured error is about 5E-16) for x from -6 to 26.5 and
        /// for tiny arguments, and to two subnormal steps in the underflow tail.
        /// </summary>
        [TestMethod]
        public void Test_Erfc_MatchMpmath()
        {
            Assert.IsGreaterThan(690, Check("erfc", (n, a, x) => Erf.Erfc(x), RelativeOrSubnormal(1E-15)));
        }

        /// <summary>
        /// The standard Normal CDF matches mpmath to 1E-15 relative (the worst measured error is about 6E-16) for z from
        /// -37.5 to 8.5, and to two subnormal steps below -37.5.
        /// </summary>
        [TestMethod]
        public void Test_StandardNormalCDF_MatchMpmath()
        {
            Assert.IsGreaterThan(920, Check("ndtr", (n, a, z) => Numerics.Distributions.Normal.StandardCDF(z), RelativeOrSubnormal(1E-15)));
        }

        /// <summary>
        /// Erfc and the standard Normal CDF keep their limits, symmetry and underflow edges.
        /// </summary>
        [TestMethod]
        public void Test_Erfc_StandardNormalCDF_Edges()
        {
            Assert.IsTrue(double.IsNaN(Erf.Erfc(double.NaN)));
            Assert.AreEqual(0.0, Erf.Erfc(double.PositiveInfinity), 0.0);
            Assert.AreEqual(2.0, Erf.Erfc(double.NegativeInfinity), 0.0);
            Assert.AreEqual(1.0, Erf.Erfc(0.0), 0.0);
            Assert.AreEqual(1.0, Erf.Erfc(-0.0), 0.0);
            Assert.AreEqual(0.0, Erf.Erfc(27.3), 0.0);
            Assert.IsGreaterThan(0.0, Erf.Erfc(27.2));
            Assert.IsTrue(double.IsNaN(Numerics.Distributions.Normal.StandardCDF(double.NaN)));
            Assert.AreEqual(0.0, Numerics.Distributions.Normal.StandardCDF(double.NegativeInfinity), 0.0);
            Assert.AreEqual(1.0, Numerics.Distributions.Normal.StandardCDF(double.PositiveInfinity), 0.0);
            Assert.AreEqual(0.5, Numerics.Distributions.Normal.StandardCDF(0.0), 1E-16);
            Assert.AreEqual(0.0, Numerics.Distributions.Normal.StandardCDF(-38.7), 0.0);
            Assert.IsGreaterThan(0.0, Numerics.Distributions.Normal.StandardCDF(-38.4));
            // Both sides of the switch to the asymptotic polynomial at x = 3: mpmath erfc(3) = 2.209049699858544E-05 and
            // erfc(2.9999999999999996) = 2.2090496998585502E-05.
            Assert.AreEqual(2.209049699858544E-05, Erf.Erfc(3.0), 1E-15 * 2.209049699858544E-05);
            Assert.AreEqual(2.2090496998585502E-05, Erf.Erfc(2.9999999999999996), 1E-15 * 2.209049699858544E-05);
            for (double z = -9; z <= 9; z += 0.37)
                Assert.AreEqual(1.0, Numerics.Distributions.Normal.StandardCDF(z) + Numerics.Distributions.Normal.StandardCDF(-z), 2E-16);
        }

        /// <summary>
        /// An infinite shape gives the limits P = 0 and Q = 1 for finite x; signed zero keeps the odd symmetry of J1.
        /// </summary>
        [TestMethod]
        public void Test_EdgeValues()
        {
            Assert.AreEqual(0.0, Gamma.LowerIncomplete(double.PositiveInfinity, 3.0), 0.0);
            Assert.AreEqual(1.0, Gamma.UpperIncomplete(double.PositiveInfinity, 3.0), 0.0);
            Assert.IsTrue(double.IsNaN(Gamma.LowerIncomplete(double.PositiveInfinity, double.PositiveInfinity)));
            Assert.IsTrue(double.IsNaN(Gamma.UpperIncomplete(2.0, double.NaN)));
            Assert.AreEqual(1.0, Gamma.LowerIncomplete(2.0, double.PositiveInfinity), 0.0);
            Assert.AreEqual(0.0, Gamma.UpperIncomplete(2.0, double.PositiveInfinity), 0.0);
            Assert.AreEqual(long.MinValue, BitConverter.DoubleToInt64Bits(Bessel.J1(-0.0)));
            Assert.IsTrue(double.IsNaN(Bessel.J0(double.NaN)));
            Assert.IsTrue(double.IsNaN(Bessel.Y1(double.NaN)));
            Assert.AreEqual(0.0, Bessel.K0(double.PositiveInfinity), 0.0);
            Assert.AreEqual(0.0, Bessel.K1(800.0), 0.0);
            Assert.IsGreaterThan(0.0, Bessel.K0(740.0));
        }
    }
}
