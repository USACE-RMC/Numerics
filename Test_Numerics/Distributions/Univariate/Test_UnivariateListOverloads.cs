using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Tests for the list forms <see cref="UnivariateDistributionBase.CCDF(IList{double})"/>,
    /// <see cref="UnivariateDistributionBase.LogPDF(IList{double})"/>, <see cref="UnivariateDistributionBase.LogCDF(IList{double})"/>
    /// and <see cref="UnivariateDistributionBase.LogCCDF(IList{double})"/>.
    /// </summary>
    /// <remarks>
    /// Every value of a list form must be bitwise the scalar member at the same point, for every family the factory
    /// creates. The SciPy 1.18.1 reference values (<c>scipy.stats.genextreme.sf</c>, <c>logsf</c>, <c>logcdf</c>,
    /// <c>scipy.stats.genpareto.sf</c> and <c>scipy.stats.norm.logpdf</c>; SciPy's shape <c>c</c> is Hosking's kappa
    /// for the GEV and its negative for the generalized Pareto) are independent oracles; they are compared to 1E-12
    /// relative error because the two libraries evaluate the same closed forms in different operation orders.
    /// </remarks>
    [TestClass]
    public class Test_UnivariateListOverloads
    {
        /// <summary>Points spanning the supports, both zeros, the far tails, the infinities and NaN.</summary>
        private static readonly double[] Points =
        {
            double.NegativeInfinity, -1E300, -50.0, -1.0, -0.0, 0.0, 1E-300, 0.25, 0.5, 1.0, 2.0, 5.0, 17.5, 50.0, 99.0,
            100.0, 150.0, 300.0, 1E6, 1E300, double.PositiveInfinity, double.NaN,
        };

        /// <summary>
        /// Asserts that each list value is bitwise the scalar value; where the scalar member throws at some point,
        /// the list form must throw the same exception type (the first point that throws stops the list).
        /// </summary>
        private static void AssertListMatchesScalar(string label, Func<IList<double>, double[]> list, Func<double, double> scalar, IList<double> points)
        {
            var expected = new double[points.Count];
            Exception scalarError = null!;
            for (int i = 0; i < points.Count && scalarError == null; i++)
            {
                try { expected[i] = scalar(points[i]); }
                catch (Exception e) { scalarError = e; }
            }

            double[] values = new double[0];
            Exception listError = null!;
            try { values = list(points); }
            catch (Exception e) { listError = e; }

            if (scalarError != null)
            {
                Assert.IsNotNull(listError, $"{label}: the scalar member threw {scalarError.GetType().Name} but the list form did not.");
                Assert.AreEqual(scalarError.GetType(), listError.GetType(), label);
                return;
            }
            Assert.IsNull(listError, $"{label}: the list form threw but the scalar member did not.");
            for (int i = 0; i < points.Count; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i]), BitConverter.DoubleToInt64Bits(values[i]), $"{label} at x = {points[i]}");
        }

        /// <summary>Every family the factory creates: the four list forms are bitwise their scalar members.</summary>
        [TestMethod]
        public void ListFormsMatchScalarMembersForEveryFamily()
        {
            int families = 0;
            foreach (UnivariateDistributionType type in Enum.GetValues(typeof(UnivariateDistributionType)))
            {
                if (!UnivariateDistributionFactory.TryCreateDistribution(type, out var distribution) || distribution == null || !distribution.ParametersValid)
                    continue;
                families++;
                var points = new List<double>(Points);
                // Points scaled to the family's own spread, so that every family sees interior values.
                double location = SafeQuantile(distribution, 0.5), spread = Math.Abs(SafeQuantile(distribution, 0.9) - SafeQuantile(distribution, 0.1));
                if (!double.IsNaN(location) && !double.IsNaN(spread))
                    for (int j = -8; j <= 8; j++) points.Add(location + 0.4 * j * spread);

                string name = type.ToString();
                AssertListMatchesScalar(name + ".CCDF", distribution.CCDF, distribution.CCDF, points);
                AssertListMatchesScalar(name + ".LogPDF", distribution.LogPDF, distribution.LogPDF, points);
                AssertListMatchesScalar(name + ".LogCDF", distribution.LogCDF, distribution.LogCDF, points);
                AssertListMatchesScalar(name + ".LogCCDF", distribution.LogCCDF, distribution.LogCCDF, points);
            }
            Assert.IsGreaterThan(30, families, $"Only {families} families were created.");
        }

        /// <summary>Returns a quantile, or NaN when the family cannot evaluate it.</summary>
        private static double SafeQuantile(UnivariateDistributionBase distribution, double probability)
        {
            try { return distribution.InverseCDF(probability); }
            catch (ArithmeticException) { return double.NaN; }
            catch (ArgumentException) { return double.NaN; }
        }

        /// <summary>SciPy oracle values for the generalized extreme value, generalized Pareto and Normal list forms.</summary>
        [TestMethod]
        public void ListFormsMatchSciPy()
        {
            var gev = new GeneralizedExtremeValue(100, 20, 0.1);
            var ccdf = gev.CCDF(new[] { -5.0, 0.0, 10.0, 50.0, 200.0 });
            var expectedCcdf = new[] { 1.0, 1.0, 1.0, 0.9999097769613546, 0.0009760858180243385 };
            for (int i = 0; i < ccdf.Length; i++) Assert.AreEqual(expectedCcdf[i], ccdf[i], 1E-12 * expectedCcdf[i]);

            var logCcdf = gev.LogCCDF(new[] { 50.0, 100.0, 150.0, 300.0 });
            var expectedLogCcdf = new[] { -9.022710898854306E-05, -0.45867514538708193, -2.9048453515334502 };
            for (int i = 0; i < 3; i++) Assert.AreEqual(expectedLogCcdf[i], logCcdf[i], 1E-12 * Math.Abs(expectedLogCcdf[i]));
            Assert.IsTrue(double.IsNegativeInfinity(logCcdf[3]));

            var logCdf = gev.LogCDF(new[] { 50.0, 100.0, 150.0, 300.0 });
            var expectedLogCdf = new[] { -9.313225746154785, -1.0, -0.056313514709472684, 0.0 };
            for (int i = 0; i < 4; i++) Assert.AreEqual(expectedLogCdf[i], logCdf[i], 1E-12 * Math.Max(1.0, Math.Abs(expectedLogCdf[i])));

            var gpa = new GeneralizedPareto(0, 10, -0.1);
            var gpaCcdf = gpa.CCDF(new[] { 0.0, 5.0, 20.0, 60.0 });
            var expectedGpa = new[] { 1.0, 0.6139132535407594, 0.1615055828898457, 0.009094947017729279 };
            for (int i = 0; i < 4; i++) Assert.AreEqual(expectedGpa[i], gpaCcdf[i], 1E-12 * expectedGpa[i]);

            var normal = new Normal(0, 1);
            var logPdf = normal.LogPDF(new[] { -30.0, 0.0, 1.5, 40.0 });
            var expectedLogPdf = new[] { -450.91893853320465, -0.9189385332046727, -2.0439385332046727, -800.9189385332047 };
            for (int i = 0; i < 4; i++) Assert.AreEqual(expectedLogPdf[i], logPdf[i], 1E-12 * Math.Abs(expectedLogPdf[i]));
        }

        /// <summary>An empty list gives an empty array; a null list is refused with <see cref="ArgumentNullException"/>.</summary>
        [TestMethod]
        public void EmptyAndNullLists()
        {
            var normal = new Normal(10, 2);
            Assert.IsEmpty(normal.CCDF(new double[0]));
            Assert.IsEmpty(normal.LogPDF(new List<double>()));
            Assert.IsEmpty(normal.LogCDF(new double[0]));
            Assert.IsEmpty(normal.LogCCDF(new double[0]));
            Assert.ThrowsExactly<ArgumentNullException>(() => normal.CCDF((IList<double>)null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => normal.LogPDF((IList<double>)null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => normal.LogCDF((IList<double>)null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => normal.LogCCDF((IList<double>)null!));
        }

        /// <summary>Invalid parameters raise the same exception as the scalar member.</summary>
        [TestMethod]
        public void InvalidParametersThrowAsTheScalarMember()
        {
            var gev = new GeneralizedExtremeValue(100, -1, 0.1);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.CCDF(150.0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.CCDF(new[] { 150.0 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.LogCCDF(new[] { 150.0 }));
        }
    }
}
