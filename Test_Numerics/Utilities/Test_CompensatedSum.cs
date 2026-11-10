using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;

namespace Utilities
{
    /// <summary>
    /// Tests for the compensated (Neumaier) sums <see cref="Tools.CompensatedSum(IList{double})"/>,
    /// <see cref="Tools.CompensatedSum(IList{double}, IList{int}, bool)"/>,
    /// <see cref="Tools.CompensatedSum(IList{double}, IList{int}, int, int)"/> and
    /// <see cref="Tools.CompensatedCumulativeSum(IList{double})"/>.
    /// </summary>
    /// <remarks>
    /// The expected values are Python 3.14 <c>math.fsum</c> results (the correctly rounded sum, Shewchuk's algorithm),
    /// checked against mpmath at 400 bits.
    /// </remarks>
    [TestClass]
    public class Test_CompensatedSum
    {
        /// <summary>A 24-term fixture spanning twenty decades with both signs (Python random.Random(20261107)).</summary>
        private static readonly double[] Mixed =
        {
            -4.415159360150449e-08, 7.317112953068489, 9.868669377084806e-07, -0.0006463016681697819, -7.502006205048129e-05,
            -0.01673876386175016, -0.0006856973900070548, -825842102090.3063, 75526526787.25735, 2.3726148037552075e-08,
            -9.761977681447285e-08, -51658895725.59849, 1.3059917476939964e-07, -65237681874.17293, 994405.7906876896,
            -2.0600436727650427e-07, -2.2111435313817664e-09, -349.5172509335271, -9.36662463129794e-09, -3.733797169087079e-09,
            0.024886648499566268, 45263.7648646135, 0.8450434057865094, 4.4580288798240146e-09,
        };

        /// <summary>
        /// Ill-conditioned sums equal math.fsum: [1e100, 1, -1e100] is 1 where the plain sum is 0, a term larger than the
        /// running sum ([1, 1e100, 1, -1e100] is 2, where Kahan summation gives 0), ten 0.1 terms are 1 where the plain sum is
        /// 0.9999999999999999, and the mixed fixture is the correctly rounded sum.
        /// </summary>
        [TestMethod]
        public void CompensatedSum_MatchesFsum()
        {
            Assert.AreEqual(1d, Tools.CompensatedSum(new[] { 1E100, 1d, -1E100 }), 0d);
            Assert.AreEqual(0d, Tools.Sum(new[] { 1E100, 1d, -1E100 }), 0d);
            Assert.AreEqual(2d, Tools.CompensatedSum(new[] { 1d, 1E100, 1d, -1E100 }), 0d);
            Assert.AreEqual(4d, Tools.CompensatedSum(new[] { 1E16, 1d, -1E16, 1E-16, 3d, -Math.Pow(2, -60) }), 0d);
            var tenths = new double[10];
            for (int i = 0; i < tenths.Length; i++) tenths[i] = 0.1;
            Assert.AreEqual(1d, Tools.CompensatedSum(tenths), 0d);
            Assert.AreEqual(0.9999999999999999, Tools.Sum(tenths), 0d);
            Assert.AreEqual(unchecked((long)0xC2693D39EE8CD39FUL), BitConverter.DoubleToInt64Bits(Tools.CompensatedSum(Mixed)));
            Assert.AreEqual(5d, Tools.CompensatedSum(new[] { 5d }), 0d);
        }

        /// <summary>
        /// Non-finite input gives the IEEE sum: NaN propagates and a single infinity is returned with its sign, as math.fsum
        /// gives, while +∞ and −∞ together give NaN and an intermediate overflow of finite terms returns +∞, where math.fsum
        /// raises.
        /// The empty list returns NaN, as <see cref="Tools.Sum(IList{double})"/> does (math.fsum returns 0).
        /// </summary>
        [TestMethod]
        public void CompensatedSum_NonFiniteAndEmpty()
        {
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(new[] { 1d, double.NaN, 2d })));
            Assert.AreEqual(double.PositiveInfinity, Tools.CompensatedSum(new[] { 1d, double.PositiveInfinity, -1E308 }));
            Assert.AreEqual(double.NegativeInfinity, Tools.CompensatedSum(new[] { double.NegativeInfinity, 1d, 1E100 }));
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(new[] { double.PositiveInfinity, 1d, double.NegativeInfinity })));
            Assert.AreEqual(double.PositiveInfinity, Tools.CompensatedSum(new[] { 1E308, 1E308, -1E308 }));
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(new double[0])));
            Assert.AreEqual("values", Assert.Throws<ArgumentNullException>(() => Tools.CompensatedSum(null!)).ParamName);
        }

        /// <summary>
        /// The indicator overload sums the values marked 1, or those marked 0 with the complement flag, as
        /// <see cref="Tools.Sum(IList{double}, IList{int}, bool)"/> selects them, with the compensated result; no selected
        /// value gives 0, and an empty list or a length mismatch gives NaN.
        /// </summary>
        [TestMethod]
        public void CompensatedSum_Indicators()
        {
            var values = new[] { 1E100, 1d, -1E100, 7d, double.NaN };
            var indicators = new[] { 1, 1, 1, 0, 2 };
            Assert.AreEqual(1d, Tools.CompensatedSum(values, indicators), 0d);
            Assert.AreEqual(7d, Tools.CompensatedSum(values, indicators, true), 0d);
            Assert.AreEqual(0d, Tools.CompensatedSum(new[] { 1d, 2d }, new[] { 0, 0 }), 0d);
            Assert.AreEqual(Tools.Sum(new[] { 0.5, 0.25, 2d }, new[] { 1, 0, 1 }), Tools.CompensatedSum(new[] { 0.5, 0.25, 2d }, new[] { 1, 0, 1 }), 0d);
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(new double[0], new int[0])));
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(new[] { 1d, 2d }, new[] { 1 })));
            Assert.AreEqual("values", Assert.Throws<ArgumentNullException>(() => Tools.CompensatedSum(null!, indicators)).ParamName);
            Assert.AreEqual("indicators", Assert.Throws<ArgumentNullException>(() => Tools.CompensatedSum(values, (IList<int>)null!)).ParamName);
        }

        /// <summary>
        /// The gather overload sums values[indices[i]] over the requested range in index-list order, equal to math.fsum of
        /// the gathered values, accepts repeated positions, returns NaN for an empty range, and refuses a range or a gathered
        /// index out of bounds.
        /// </summary>
        [TestMethod]
        public void CompensatedSum_Gather()
        {
            var values = new[] { 1E100, 1d, -1E100, 3d };
            var indices = new[] { 3, 0, 1, 2, 1, 1 };
            Assert.AreEqual(1d, Tools.CompensatedSum(values, indices, 1, 3), 0d);
            Assert.AreEqual(2d, Tools.CompensatedSum(values, indices, 4, 2), 0d);
            Assert.AreEqual(4d, Tools.CompensatedSum(values, indices, 0, 4), 0d);
            var all = new int[Mixed.Length];
            for (int i = 0; i < all.Length; i++) all[i] = i;
            Assert.AreEqual(Tools.CompensatedSum(Mixed), Tools.CompensatedSum(Mixed, all, 0, all.Length), 0d);
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(values, indices, 2, 0)));
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(values, indices, 6, 0)));
            Assert.AreEqual("start", Assert.Throws<ArgumentOutOfRangeException>(() => Tools.CompensatedSum(values, indices, -1, 2)).ParamName);
            Assert.AreEqual("count", Assert.Throws<ArgumentOutOfRangeException>(() => Tools.CompensatedSum(values, indices, 0, -1)).ParamName);
            Assert.AreEqual("count", Assert.Throws<ArgumentOutOfRangeException>(() => Tools.CompensatedSum(values, indices, 4, 3)).ParamName);
            Assert.AreEqual("indices", Assert.Throws<ArgumentOutOfRangeException>(() => Tools.CompensatedSum(values, new[] { 0, 4 }, 0, 2)).ParamName);
            Assert.AreEqual("indices", Assert.Throws<ArgumentOutOfRangeException>(() => Tools.CompensatedSum(values, new[] { -1 }, 0, 1)).ParamName);
            Assert.AreEqual("values", Assert.Throws<ArgumentNullException>(() => Tools.CompensatedSum(null!, indices, 0, 1)).ParamName);
            Assert.AreEqual("indices", Assert.Throws<ArgumentNullException>(() => Tools.CompensatedSum(values, null!, 0, 1)).ParamName);
        }

        /// <summary>
        /// The cumulative sums equal math.fsum of each prefix, and each entry is bitwise the compensated sum of the prefix,
        /// including after the running sum becomes infinite or NaN; the empty list gives an empty array.
        /// </summary>
        [TestMethod]
        public void CompensatedCumulativeSum_MatchesPrefixSums()
        {
            var values = new[] { 1E100, 1d, -1E100, 1E-100, 0.1, 0.2 };
            CollectionAssert.AreEqual(new[] { 1E100, 1E100, 1d, 1d, 1.1, 1.3 }, Tools.CompensatedCumulativeSum(values));

            var withInfinities = new List<double>(Mixed) { double.PositiveInfinity, 1d, double.NegativeInfinity, 2d };
            var cumulative = Tools.CompensatedCumulativeSum(withInfinities);
            Assert.HasCount(withInfinities.Count, cumulative);
            for (int i = 0; i < withInfinities.Count; i++)
            {
                double prefix = Tools.CompensatedSum(withInfinities.GetRange(0, i + 1));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(prefix), BitConverter.DoubleToInt64Bits(cumulative[i]), $"entry {i}");
            }
            Assert.AreEqual(double.PositiveInfinity, cumulative[Mixed.Length]);
            Assert.IsTrue(double.IsNaN(cumulative[cumulative.Length - 1]));
            Assert.IsEmpty(Tools.CompensatedCumulativeSum(new double[0]));
            Assert.AreEqual("values", Assert.Throws<ArgumentNullException>(() => Tools.CompensatedCumulativeSum(null!)).ParamName);
        }

        /// <summary>
        /// The Kahan option reproduces the classic loop a = v − c, f = t + a, c = (f − t) − a, t = f bit for bit, for the
        /// plain, indicator, gather and cumulative forms, on arrays and on other lists; without it every form keeps its
        /// Neumaier bits.
        /// </summary>
        [TestMethod]
        public void KahanOptionMatchesClassicLoopBitwise()
        {
            var random = new Random(4);
            foreach (int n in new[] { 1, 2, 7, 50, 5000 })
            {
                var values = new double[n];
                for (int i = 0; i < n; i++) values[i] = random.NextDouble() * (random.NextDouble() < 0.1 ? 1e6 : 1d) * (i % 3 == 0 ? -1 : 1);
                var indicators = new int[n];
                var indices = new int[n];
                for (int i = 0; i < n; i++) { indicators[i] = random.Next(2); indices[i] = random.Next(n); }

                double t = 0, c = 0;
                var running = new double[n];
                for (int i = 0; i < n; i++) { double a = values[i] - c; double f = t + a; c = (f - t) - a; t = f; running[i] = t; }
                double ti = 0, ci = 0, tg = 0, cg = 0;
                for (int i = 0; i < n; i++)
                {
                    if (indicators[i] == 1) { double a = values[i] - ci; double f = ti + a; ci = (f - ti) - a; ti = f; }
                    double b = values[indices[i]] - cg; double g = tg + b; cg = (g - tg) - b; tg = g;
                }

                var list = new List<double>(values);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(t), BitConverter.DoubleToInt64Bits(Tools.CompensatedSum(values, true)));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(t), BitConverter.DoubleToInt64Bits(Tools.CompensatedSum(list, true)));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(ti), BitConverter.DoubleToInt64Bits(Tools.CompensatedSum(values, indicators, false, true)));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(tg), BitConverter.DoubleToInt64Bits(Tools.CompensatedSum(values, indices, 0, n, true)));
                CollectionAssert.AreEqual(running, Tools.CompensatedCumulativeSum(values, true));
                CollectionAssert.AreEqual(running, Tools.CompensatedCumulativeSum(list, true));

                Assert.AreEqual(Tools.CompensatedSum(list), Tools.CompensatedSum(values, false), 0d);
                Assert.AreEqual(Tools.CompensatedSum(list), Tools.CompensatedSum(values), 0d);
                Assert.AreEqual(Tools.CompensatedSum(values, indicators), Tools.CompensatedSum(values, indicators, false, false), 0d);
                Assert.AreEqual(Tools.CompensatedSum(values, indices, 0, n), Tools.CompensatedSum(values, indices, 0, n, false), 0d);
                CollectionAssert.AreEqual(Tools.CompensatedCumulativeSum(list), Tools.CompensatedCumulativeSum(values, false));
            }

            // Kahan loses the large term that Neumaier keeps; once the sum is infinite the Kahan correction is NaN.
            var hard = new[] { 1d, 1e100, 1d, -1e100 };
            Assert.AreEqual(2d, Tools.CompensatedSum(hard), 0d);
            Assert.AreEqual(0d, Tools.CompensatedSum(hard, true), 0d);
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(new[] { 1d, double.PositiveInfinity, 1d }, true)));
            Assert.AreEqual(double.PositiveInfinity, Tools.CompensatedSum(new[] { 1d, double.PositiveInfinity }, true));
            Assert.AreEqual(double.PositiveInfinity, Tools.CompensatedSum(new[] { 1d, double.PositiveInfinity, 1d }));
            Assert.IsTrue(double.IsNaN(Tools.CompensatedSum(new double[0], true)));
            Assert.IsEmpty(Tools.CompensatedCumulativeSum(new double[0], true));
            Assert.ThrowsExactly<ArgumentNullException>(() => Tools.CompensatedSum(null!, true));
        }
    }
}
