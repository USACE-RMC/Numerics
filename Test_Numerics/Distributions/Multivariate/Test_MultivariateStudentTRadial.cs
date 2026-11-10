using System;
using System.Globalization;
using System.Linq;
using Distributions.Univariate;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;

namespace Distributions.Multivariate
{
    /// <summary>
    /// Tests the inverse incomplete gamma ratio on the radial coordinate of the multivariate t integrand.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// The quantile oracle is the gammapinv rows of docs/mathematics/oracles/special-functions/special-functions-mpmath.csv,
    /// roots of P(a, x) = p solved in mpmath at 50 digits. The three-dimensional CDF reference is R 4.4.3 mvtnorm pmvt.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_MultivariateStudentTRadial
    {
        /// <summary>Parses an oracle number written by Python repr.</summary>
        private static double Parse(string text) => double.Parse(text, CultureInfo.InvariantCulture);

        /// <summary>
        /// The radial quantile matches the mpmath root of P(a, x) = p for a from 0.01 to 1E5 and p from 1E-300 to 0.999,
        /// to the tolerance of the general inverse incomplete gamma test.
        /// </summary>
        [TestMethod]
        public void Test_RadialQuantile_MatchesMpmath()
        {
            var rows = DistributionOracle.Read("special-functions-mpmath.csv").Where(r => r["function"] == "gammapinv")
                .Select(r => (A: Parse(r["a"]), P: Parse(r["x"]), X: Parse(r["value"]))).ToList();
            int count = 0;
            foreach (var (a, p, expected) in rows)
            {
                if (!(p < 1d)) continue;
                double actual = Gamma.InverseIncompleteRatio(a, p, false, Gamma.LogGamma(a));
                Assert.AreEqual(expected, actual, (5E-14 + 2E-16 / a) * expected, $"a={a:R}, p={p:R}");
                count++;
            }
            Assert.IsGreaterThan(150, count);
        }

        /// <summary>
        /// Uniform coordinates near one keep their relative precision through the exact complement 1 - u, and the
        /// residual of the root is at rounding level on both sides of u = 1/2.
        /// </summary>
        [TestMethod]
        public void Test_RadialQuantile_UpperTailAndResidual()
        {
            foreach (double a in new[] { 0.5, 1.25, 2.0, 15.0, 50.0 })
            {
                double logGamma = Gamma.LogGamma(a);
                foreach (double u in new[] { 1E-12, 0.01, 0.3, 0.4999, 0.5, 0.7, 0.99, 1 - 1E-10, 1d - 1.1102230246251565E-16 })
                {
                    double x = Gamma.InverseIncompleteRatio(a, u, false, logGamma);
                    double residual = u < 0.5 ? Gamma.LowerIncomplete(a, x) / u - 1 : Gamma.UpperIncomplete(a, x) / (1 - u) - 1;
                    Assert.IsLessThan(2E-14, Math.Abs(residual), $"a={a}, u={u:R}");
                }
            }
        }

        /// <summary>
        /// A zero coordinate gives zero, negative and NaN coordinates give NaN, and the radial factor is sqrt(x / a).
        /// </summary>
        [TestMethod]
        public void Test_RadialQuantile_Edges()
        {
            double logGamma = Gamma.LogGamma(2.0);
            Assert.AreEqual(0.0, MultivariateStudentT.RadialFactor(2.0, logGamma, 0.0), 0.0);
            Assert.IsTrue(double.IsNaN(MultivariateStudentT.RadialFactor(2.0, logGamma, -0.5)));
            Assert.IsTrue(double.IsNaN(MultivariateStudentT.RadialFactor(2.0, logGamma, double.NaN)));
            double x = Gamma.InverseIncompleteRatio(2.0, 0.3, false, logGamma);
            Assert.AreEqual(Math.Sqrt(x / 2.0), MultivariateStudentT.RadialFactor(2.0, logGamma, 0.3), 0.0);
            // The public inverse is the same method.
            Assert.AreEqual(Gamma.InverseLowerIncomplete(2.0, 0.3), x, 0.0);
        }

        /// <summary>
        /// The three-dimensional CDF with correlated scale and ν = 4 matches R pmvt (GenzBretz(maxpts = 5E7, abseps = 1E-11),
        /// 0.420705484459314 with error 5.6E-9) within its default error control, and to 1E-6 with a tight budget.
        /// </summary>
        [TestMethod]
        public void Test_ThreeDimensionalCDF_MatchesR()
        {
            var scale = new double[,] { { 1, 0.5, 0.3 }, { 0.5, 1, 0.4 }, { 0.3, 0.4, 1 } };
            var point = new[] { 0.5, 0.2, 1.0 };
            const double reference = 0.420705484459314;
            var t = new MultivariateStudentT(4.0, new double[3], scale);
            Assert.AreEqual(reference, t.CDF(point), 2E-5);
            var tight = new MultivariateStudentT(4.0, new double[3], scale) { AbsoluteError = 1E-7, MaxEvaluations = 2_000_000 };
            tight.MVNUNI = new MersenneTwister(11);
            Assert.AreEqual(reference, tight.CDF(point), 1E-6);
        }
    }
}
