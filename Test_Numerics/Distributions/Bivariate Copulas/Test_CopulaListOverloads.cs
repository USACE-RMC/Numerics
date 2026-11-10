using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions.Copulas;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Tests for the paired-list forms of <see cref="BivariateCopula.PDF(IList{double}, IList{double})"/>,
    /// <see cref="BivariateCopula.CDF(IList{double}, IList{double})"/>,
    /// <see cref="BivariateCopula.ConditionalCDF(IList{double}, IList{double})"/> and
    /// <see cref="BivariateCopula.InverseConditionalCDF(IList{double}, IList{double})"/>.
    /// </summary>
    /// <remarks>
    /// Each list value must be bitwise the scalar member for every shipped copula. The Normal copula's conditional
    /// distribution has the closed form h(v|u) = Φ((Φ⁻¹(v) − ρΦ⁻¹(u)) / √(1 − ρ²)) and its inverse
    /// Φ(ρΦ⁻¹(u) + √(1 − ρ²)Φ⁻¹(t)); the reference values were computed with SciPy 1.18.1 (<c>scipy.special.ndtr</c>
    /// and <c>ndtri</c>) and are compared to 1E-12, the agreement of two independent normal-integral implementations.
    /// </remarks>
    [TestClass]
    public class Test_CopulaListOverloads
    {
        /// <summary>One copula of every shipped family.</summary>
        private static IEnumerable<BivariateCopula> Copulas()
        {
            yield return new AMHCopula(0.5);
            yield return new ClaytonCopula(2.0);
            yield return new FrankCopula(5.0);
            yield return new GumbelCopula(2.0);
            yield return new IndependenceCopula();
            yield return new JoeCopula(2.0);
            yield return new NormalCopula(0.6);
            yield return new StudentTCopula(0.6, 5.0);
        }

        /// <summary>The four list forms are bitwise their scalar members for every family, the unit-square edges included.</summary>
        [TestMethod]
        public void ListFormsMatchScalarMembers()
        {
            var random = new Random(20261023);
            var u = new List<double> { 0.0, 1.0, 0.5, 1E-12, 1 - 1E-12, 0.0, 1.0 };
            var v = new List<double> { 0.5, 0.5, 0.0, 0.3, 0.7, 1.0, 1.0 };
            for (int i = 0; i < 200; i++)
            {
                u.Add(random.NextDouble());
                v.Add(random.NextDouble());
            }
            foreach (var copula in Copulas())
            {
                string name = copula.GetType().Name;
                var pdf = copula.PDF(u, v);
                var cdf = copula.CDF(u, v);
                var conditional = copula.ConditionalCDF(u, v);
                var inverse = copula.InverseConditionalCDF(u, v);
                for (int i = 0; i < u.Count; i++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(copula.PDF(u[i], v[i])), BitConverter.DoubleToInt64Bits(pdf[i]), $"{name}.PDF[{i}]");
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(copula.CDF(u[i], v[i])), BitConverter.DoubleToInt64Bits(cdf[i]), $"{name}.CDF[{i}]");
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(copula.ConditionalCDF(u[i], v[i])), BitConverter.DoubleToInt64Bits(conditional[i]), $"{name}.ConditionalCDF[{i}]");
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(copula.InverseConditionalCDF(u[i], v[i])), BitConverter.DoubleToInt64Bits(inverse[i]), $"{name}.InverseConditionalCDF[{i}]");
                }
            }
        }

        /// <summary>The conditional inverse list equals the second coordinate of <see cref="BivariateCopula.InverseCDF(double, double)"/>, bit for bit.</summary>
        [TestMethod]
        public void InverseConditionalListEqualsInverseCdfSecondCoordinate()
        {
            var random = new Random(5);
            var u = new double[100];
            var t = new double[100];
            for (int i = 0; i < u.Length; i++)
            {
                u[i] = random.NextDouble();
                t[i] = random.NextDouble();
            }
            foreach (var copula in Copulas())
            {
                var inverse = copula.InverseConditionalCDF(u, t);
                for (int i = 0; i < u.Length; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(copula.InverseCDF(u[i], t[i])[1]), BitConverter.DoubleToInt64Bits(inverse[i]), copula.GetType().Name);
            }
        }

        /// <summary>SciPy oracle values of the Normal copula's conditional distribution and its inverse.</summary>
        [TestMethod]
        public void NormalCopulaMatchesSciPy()
        {
            var copula = new NormalCopula(0.6);
            var u = new[] { 0.1, 0.5, 0.9, 0.3 };
            var t = new[] { 0.2, 0.5, 0.95, 0.7 };
            var inverse = copula.InverseConditionalCDF(u, t);
            var expectedInverse = new[] { 0.07461904244765404, 0.5, 0.9814568990712572, 0.5417645261243141 };
            var conditional = copula.ConditionalCDF(u, t);
            var expectedConditional = new[] { 0.4638007776485923, 0.5, 0.8632205114546854, 0.8528651472793856 };
            for (int i = 0; i < u.Length; i++)
            {
                Assert.AreEqual(expectedInverse[i], inverse[i], 1E-12);
                Assert.AreEqual(expectedConditional[i], conditional[i], 1E-12);
            }
        }

        /// <summary>Null lists and lists of different lengths are refused; empty lists give empty arrays.</summary>
        [TestMethod]
        public void ArgumentErrors()
        {
            var copula = new NormalCopula(0.6);
            var one = new[] { 0.5 };
            var two = new[] { 0.5, 0.5 };
            Assert.ThrowsExactly<ArgumentNullException>(() => copula.PDF(null!, one));
            Assert.ThrowsExactly<ArgumentNullException>(() => copula.CDF(one, (IList<double>)null!));
            Assert.ThrowsExactly<ArgumentException>(() => copula.ConditionalCDF(one, two));
            Assert.ThrowsExactly<ArgumentException>(() => copula.InverseConditionalCDF(two, one));
            Assert.IsEmpty(copula.PDF(new double[0], new double[0]));
            Assert.IsEmpty(copula.InverseConditionalCDF(new List<double>(), new List<double>()));
        }
    }
}
