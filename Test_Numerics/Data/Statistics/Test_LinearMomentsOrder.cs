using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Stats = Numerics.Data.Statistics.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Tests for <see cref="Numerics.Data.Statistics.Statistics.LinearMoments(System.Collections.Generic.IList{double}, int)"/>.
    /// </summary>
    /// <remarks>
    /// The expected values are the exact rational evaluation of the unbiased estimator of Hosking (1990) on the double
    /// data, rounded once to double: Python 3.13 <c>fractions.Fraction</c> sums of x(i) C(i − 1, r) / (n C(n − 1, r)),
    /// then λ_(r+1) = Σ_k (−1)^(r−k) (r + k)! / ((k!)² (r − k)!) b_k and τr = λr / λ2.
    /// </remarks>
    [TestClass]
    public class Test_LinearMomentsOrder
    {
        /// <summary>Ten values with mean 13.3.</summary>
        private static readonly double[] Small = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

        /// <summary>Thirty skewed values (seeded lognormal, rounded to 0.1).</summary>
        private static readonly double[] Skewed =
        {
            4865.5, 1203.2, 2857.6, 3313.4, 4127.6, 4252.0, 3423.6, 2036.1, 3378.9, 5313.4, 2511.6, 1780.8, 2340.7, 5604.9, 2068.4,
            3727.2, 2150.0, 829.8, 2325.6, 3338.6, 1574.5, 12000.0, 3189.4, 2370.2, 2625.4, 5201.9, 2754.7, 1153.2, 4450.0, 5614.7
        };

        /// <summary>Forty values in [1000, 1015], whose mean is 429 times λ2.</summary>
        private static readonly double[] Offset =
        {
            1001.787, 1007.538, 1007.677, 1012.9, 1001.54, 1003.349, 1009.015, 1008.348, 1011.751, 1008.217, 1010.958, 1011.522,
            1011.265, 1008.798, 1003.599, 1009.213, 1001.661, 1012.251, 1006.745, 1012.221, 1010.278, 1010.191, 1003.153, 1003.773,
            1014.699, 1013.94, 1012.077, 1014.995, 1007.704, 1001.145, 1005.106, 1008.362, 1004.382, 1006.862, 1013.147, 1010.56,
            1008.661, 1000.371, 1010.663, 1005.453
        };

        /// <summary>
        /// Known values of λ1, λ2 and τ3 to τ6 against the exact rational oracle.
        /// </summary>
        [TestMethod]
        public void KnownValuesMatchExactRationalOracle()
        {
            double[] small = { 13.3, 1.16, 0.06752873563218392, 0.09995894909688034, 0.10296934865900396, 0.1307471264367814 };
            var actual = Stats.LinearMoments(Small, 6);
            Assert.HasCount(6, actual);
            double[] smallTolerance = { 1E-14, 1E-14, 1E-14, 1E-14, 1E-13, 1E-12 };
            for (int j = 0; j < 6; j++) Assert.AreEqual(small[j], actual[j], smallTolerance[j], $"small, j = {j}");

            double[] skewed = { 3412.7633333333333, 1032.017816091954, 0.2762114202754707, 0.24966756546615443, 0.15736743996068578 };
            actual = Stats.LinearMoments(Skewed, 5);
            Assert.AreEqual(skewed[0], actual[0], 1E-15 * skewed[0]);
            Assert.AreEqual(skewed[1], actual[1], 1E-15 * skewed[1]);
            for (int j = 2; j < 5; j++) Assert.AreEqual(skewed[j], actual[j], 1E-13, $"skewed, j = {j}");

            double[] offset = { 1008.146925, 2.3490788461538474, -0.0813330696988244, 0.05123784048803355, 0.03459415861705856 };
            actual = Stats.LinearMoments(Offset, 5);
            Assert.AreEqual(offset[0], actual[0], 1E-15 * offset[0]);
            Assert.AreEqual(offset[1], actual[1], 1E-13 * offset[1]);
            for (int j = 2; j < 4; j++) Assert.AreEqual(offset[j], actual[j], 1E-12, $"offset, j = {j}");
            Assert.AreEqual(offset[4], actual[4], 1E-11);
        }

        /// <summary>
        /// The first four values are bitwise the four-value member, for every order and for small and large samples, so
        /// the λ2 cancellation of large means is no worse than in that member; lower orders are prefixes of higher ones.
        /// </summary>
        [TestMethod]
        public void LowerOrdersMatchTheFourValueMemberBitwise()
        {
            var random = new Random(31);
            foreach (int n in new[] { 4, 5, 10, 57, 2000, 50000 })
            {
                var x = new double[n];
                for (int i = 0; i < n; i++) x[i] = 1E3 + Math.Exp(random.NextDouble() * 3d);
                var four = Stats.LinearMoments(x);
                for (int nmom = 1; nmom <= Math.Min(n, 7); nmom++)
                {
                    var moments = Stats.LinearMoments(x, nmom);
                    Assert.HasCount(nmom, moments);
                    for (int j = 0; j < Math.Min(nmom, 4); j++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(four[j]), BitConverter.DoubleToInt64Bits(moments[j]), $"n = {n}, nmom = {nmom}, j = {j}");
                    var seven = Stats.LinearMoments(x, Math.Min(n, 7));
                    for (int j = 0; j < nmom; j++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(seven[j]), BitConverter.DoubleToInt64Bits(moments[j]), $"prefix n = {n}, nmom = {nmom}, j = {j}");
                }
            }
        }

        /// <summary>
        /// The four-value member is bitwise the reference accumulation below, the sequential probability weighted moment
        /// sums with each coefficient divided in its term, on samples up to 50000 values.
        /// </summary>
        [TestMethod]
        public void FourValueMemberMatchesSequentialReferenceBitwise()
        {
            var random = new Random(47);
            foreach (int n in new[] { 4, 9, 100, 1500, 50000 })
            {
                var x = new double[n];
                for (int i = 0; i < n; i++) x[i] = 50d + 10d * Math.Exp(random.NextDouble());
                var expected = ReferenceFourMoments(x);
                var actual = Stats.LinearMoments(x);
                for (int j = 0; j < 4; j++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[j]), BitConverter.DoubleToInt64Bits(actual[j]), $"n = {n}, j = {j}");
            }
        }

        /// <summary>
        /// The sequential reference: b0 to b3 with the coefficients (i − 1)/(N − 1), (i − 2)(i − 1)/((N − 2)(N − 1)) and
        /// (i − 3)(i − 2)(i − 1)/((N − 3)(N − 2)(N − 1)), then the closed forms of λ2, τ3 and τ4.
        /// </summary>
        private static double[] ReferenceFourMoments(double[] data)
        {
            var sorted = (double[])data.Clone();
            Array.Sort(sorted);
            double N = sorted.Length;
            double B0 = 0, B1 = 0, B2 = 0, B3 = 0;
            for (int i = 1; i <= N; i++)
            {
                double di = i;
                B0 += sorted[i - 1];
                if (i > 1) B1 += (i - 1) / (N - 1) * sorted[i - 1];
                if (i > 2) B2 += (di - 2) * (di - 1) / ((N - 2) * (N - 1)) * sorted[i - 1];
                if (i > 3) B3 += (di - 3) * (di - 2) * (di - 1) / ((N - 3) * (N - 2) * (N - 1)) * sorted[i - 1];
            }
            B0 /= N; B1 /= N; B2 /= N; B3 /= N;
            return new[] { B0, 2 * B1 - B0, 2 * (3 * B2 - B0) / (2 * B1 - B0) - 3, 5 * (2 * (2 * B3 - 3 * B2) + B0) / (2 * B1 - B0) + 6 };
        }

        /// <summary>
        /// Exact L-moments of linear data: for x(i) = i, λ1 = (n + 1)/2, λ2 = (n + 1)/6 and every higher λr is zero, so the
        /// ratios vanish; the orders also stay finite for a large sample.
        /// </summary>
        [TestMethod]
        public void LinearDataHasVanishingHigherRatios()
        {
            int n = 101;
            var x = new double[n];
            for (int i = 0; i < n; i++) x[i] = n - i;
            var moments = Stats.LinearMoments(x, 8);
            Assert.AreEqual(51d, moments[0], 1E-13);
            Assert.AreEqual(17d, moments[1], 1E-13);
            for (int j = 2; j < 8; j++) Assert.AreEqual(0d, moments[j], 1E-11, $"j = {j}");
        }

        /// <summary>
        /// Too few values give NaN of the requested length, as does a NaN entry; null data and an order below one are refused.
        /// </summary>
        [TestMethod]
        public void EdgeCasesAndArgumentErrors()
        {
            var shortSample = Stats.LinearMoments(new[] { 1d, 2d, 3d, 4d }, 5);
            Assert.HasCount(5, shortSample);
            foreach (double value in shortSample) Assert.IsTrue(double.IsNaN(value));

            Assert.AreEqual(7d, Stats.LinearMoments(new[] { 7d }, 1)[0]);
            var pair = Stats.LinearMoments(new[] { 1d, 3d }, 2);
            Assert.AreEqual(2d, pair[0]);
            Assert.AreEqual(1d, pair[1]);

            foreach (double value in Stats.LinearMoments(new[] { 1d, double.NaN, 3d, 4d, 5d, 6d }, 5)) Assert.IsTrue(double.IsNaN(value));

            var constant = Stats.LinearMoments(new[] { 3.5, 3.5, 3.5, 3.5, 3.5, 3.5 }, 5);
            Assert.AreEqual(3.5, constant[0]);
            Assert.AreEqual(0d, constant[1]);
            for (int j = 2; j < 5; j++) Assert.IsTrue(double.IsNaN(constant[j]) || double.IsInfinity(constant[j]), $"j = {j}");

            Assert.ThrowsExactly<ArgumentNullException>(() => Stats.LinearMoments(null!, 5));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Stats.LinearMoments(Small, 0));
        }
    }
}
