#nullable disable
using System;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Stats = Numerics.Data.Statistics.Statistics;
using PPType = Numerics.Data.Statistics.PlottingPositions.PlottingPostionType;

namespace Data.Statistics
{
    /// <summary>
    /// The randomized plotting-position percentile comparisons, shared by the unit tests and the full sweeps in
    /// Benchmarks_Numerics.
    /// </summary>
    public partial class Test_PlottingPositionPercentile
    {
        /// <summary>
        /// For Blom and Hazen the weight formed from k and the bracketing positions, (k − p_r) / (p_(r+1) − p_r) with
        /// p_r = (r − α)/(n + 1 − 2α), gives percentiles at least as close to exact rational evaluation as h − floor(h), in
        /// mean and maximum error, at the return-period probabilities 1 − 1/T over seeded samples with n = 5 to 120
        /// and T = 2 to 100. The exact value is evaluated in integers, as in the Weibull test, with α = a/8 and the rank
        /// scaled by 8.
        /// </summary>
        /// <param name="cases">The number of random cases.</param>
        public static void BlomAndHazenWeightCloserToExactRationalSweep(int cases)
        {
            foreach (var (type, eighths) in new[] { (PPType.Blom, 3), (PPType.Hazen, 4) })
            {
                double alpha = eighths / 8d;
                var random = new Random(12345);
                var scale = BigInteger.Pow(2, 1074);
                var scale8 = 8 * scale;
                double maxBracket = 0, maxRank = 0, sumBracket = 0, sumRank = 0;
                int count = 0;
                for (int c = 0; c < cases; c++)
                {
                    int n = random.Next(5, 121);
                    int T = random.Next(2, 101);
                    double k = 1d - 1d / T;
                    var x = new double[n];
                    for (int i = 0; i < n; i++) x[i] = Math.Exp(3d + 2d * (random.NextDouble() - 0.5));
                    Array.Sort(x);
                    double h = (n + 1 - 2d * alpha) * k + alpha;
                    if (h <= 1d || h >= n) continue;
                    int r = (int)Math.Floor(h);
                    double rankWeightValue = x[r - 1] + (h - r) * (x[r] - x[r - 1]);
                    double bracketWeightValue = Stats.Percentile(x, k, type, true);
                    // 8h in units of 2^-1074, and the exact value scaled by 8 * 2^2148.
                    var rank8 = (8 * (n + 1) - 2 * eighths) * Scaled(k) + eighths * scale;
                    int j = (int)BigInteger.Divide(rank8, scale8);
                    if (j < 1 || j >= n) continue;
                    var lo = Scaled(x[j - 1]);
                    var exact = lo * scale8 + (rank8 - j * scale8) * (Scaled(x[j]) - lo);
                    double ulpBracket = UlpError(bracketWeightValue, exact, scale8), ulpRank = UlpError(rankWeightValue, exact, scale8);
                    maxBracket = Math.Max(maxBracket, ulpBracket); maxRank = Math.Max(maxRank, ulpRank);
                    sumBracket += ulpBracket; sumRank += ulpRank;
                    count++;
                }
                Assert.IsGreaterThan(cases / 2, count);
                Assert.IsLessThan(sumRank / count, sumBracket / count, $"{type} mean {sumBracket / count} vs {sumRank / count}");
                Assert.IsLessThanOrEqualTo(maxRank, maxBracket, $"{type} max {maxBracket} vs {maxRank}");
            }
        }

        /// <summary>
        /// The Blom and Hazen weights themselves, against the exact weight h − floor(h) evaluated in integers: on the
        /// cases of k uniform in (0, 1) with n = 2 to 1000, of k = j/1000 with n = 5 to 200, and of k = 1 − 1/T with
        /// n = 5 to 120, the bracket weight has a lower mean and a maximum no larger than h − floor(h). The library
        /// result is checked to be the interpolation with exactly that weight.
        /// </summary>
        /// <param name="cases">The number of random cases.</param>
        public static void BlomAndHazenWeightCloserOnProbabilityGridsSweep(int cases)
        {
            var scale = BigInteger.Pow(2, 1074);
            var unit = BigInteger.Pow(2, 1074 - 53 + 3);
            var grids = new Func<Random, (int n, double k)>[]
            {
                rng => (rng.Next(2, 1001), rng.NextDouble()),
                rng => (rng.Next(5, 201), rng.Next(1, 1000) / 1000d),
                rng => (rng.Next(5, 121), 1d - 1d / rng.Next(2, 101)),
            };
            foreach (var (type, eighths) in new[] { (PPType.Blom, 3), (PPType.Hazen, 4) })
            {
                double alpha = eighths / 8d;
                for (int g = 0; g < grids.Length; g++)
                {
                    var random = new Random(3 + g);
                    double maxBracket = 0, maxRank = 0, sumBracket = 0, sumRank = 0;
                    int count = 0;
                    for (int c = 0; c < cases; c++)
                    {
                        var (n, k) = grids[g](random);
                        double h = (n + 1 - 2d * alpha) * k + alpha;
                        if (h <= 1d || h >= n) continue;
                        int r = (int)Math.Floor(h);
                        var rank8 = (8 * (n + 1) - 2 * eighths) * Scaled(k) + eighths * scale;
                        if ((int)BigInteger.Divide(rank8, 8 * scale) != r) continue;
                        var exactWeight8 = rank8 - r * 8 * scale;
                        double d = n + 1 - 2d * alpha;
                        double pLower = (r - alpha) / d;
                        double bracket = Math.Min(Math.Max((k - pLower) / ((r + 1 - alpha) / d - pLower), 0d), 1d);
                        double rankWeight = h - r;
                        double errorBracket = (double)(BigInteger.Abs(8 * Scaled(bracket) - exactWeight8) * 1000 / unit) / 1000d;
                        double errorRank = (double)(BigInteger.Abs(8 * Scaled(rankWeight) - exactWeight8) * 1000 / unit) / 1000d;
                        maxBracket = Math.Max(maxBracket, errorBracket); maxRank = Math.Max(maxRank, errorRank);
                        sumBracket += errorBracket; sumRank += errorRank;
                        count++;
                        if (c % 97 == 0)
                        {
                            var x = new double[n];
                            for (int i = 0; i < n; i++) x[i] = i + 1 + 0.5 * random.NextDouble();
                            double expected = x[r - 1] + bracket * (x[r] - x[r - 1]);
                            Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(Stats.Percentile(x, k, type, true)), $"{type} grid {g} case {c}");
                        }
                    }
                    Assert.IsGreaterThan(cases / 2, count);
                    Assert.IsLessThan(sumRank / count, sumBracket / count, $"{type} grid {g} mean {sumBracket / count} vs {sumRank / count}");
                    Assert.IsLessThanOrEqualTo(maxRank, maxBracket, $"{type} grid {g} max {maxBracket} vs {maxRank}");
                }
            }
        }

        /// <summary>
        /// The Weibull weight is unchanged, bit for bit (k − r/d) / ((r + 1)/d − r/d) with d = n + 1, and the Median,
        /// Cunnane and Gringorten formulas keep the weight h − floor(h).
        /// </summary>
        /// <param name="cases">The number of random cases.</param>
        public static void OtherFormulasKeepTheirWeightsBitwiseSweep(int cases)
        {
            var random = new Random(8);
            for (int c = 0; c < cases; c++)
            {
                int n = random.Next(2, 300);
                double k = random.NextDouble();
                var x = new double[n];
                for (int i = 0; i < n; i++) x[i] = Math.Exp(4d * random.NextDouble());
                Array.Sort(x);

                double hw = (n + 1) * k;
                double weibull;
                if (hw <= 1d) weibull = x[0];
                else if (hw >= n) weibull = x[n - 1];
                else
                {
                    int r = (int)Math.Floor(hw);
                    double d = n + 1;
                    double pLower = r / d;
                    double w = Math.Min(Math.Max((k - pLower) / ((r + 1) / d - pLower), 0d), 1d);
                    weibull = x[r - 1] + w * (x[r] - x[r - 1]);
                }
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(weibull), BitConverter.DoubleToInt64Bits(Stats.Percentile(x, k, PPType.Weibull, true)), $"Weibull case {c}");

                foreach (var type in new[] { PPType.Median, PPType.Cunnane, PPType.Gringorten })
                {
                    double alpha = PlottingPositions.Alpha(type);
                    double h = (n + 1 - 2d * alpha) * k + alpha;
                    double expected;
                    if (h <= 1d) expected = x[0];
                    else if (h >= n) expected = x[n - 1];
                    else
                    {
                        int r = (int)Math.Floor(h);
                        expected = x[r - 1] + (h - r) * (x[r] - x[r - 1]);
                    }
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(Stats.Percentile(x, k, type, true)), $"{type} case {c}");
                }
            }
        }

        /// <summary>
        /// The weighted plotting-position percentile of a list of probabilities equals the scalar member bit for bit, for
        /// every formula, sorted and unsorted data, ascending and unordered probabilities with repeats and the ends, and
        /// list inputs that are not arrays.
        /// </summary>
        /// <param name="trials">The number of random trials.</param>
        public static void WeightedListOfProbabilitiesMatchesScalarBitwiseSweep(int trials)
        {
            var random = new Random(91);
            var types = (PPType[])Enum.GetValues(typeof(PPType));
            for (int trial = 0; trial < trials; trial++)
            {
                int n = random.Next(1, 80);
                var x = new double[n];
                var w = new double[n];
                for (int i = 0; i < n; i++)
                {
                    x[i] = Math.Round(random.NextDouble() * 20d, random.Next(0, 3));
                    w[i] = random.NextDouble() < 0.1 ? 0d : random.NextDouble() * Math.Pow(10, random.Next(-6, 3));
                }
                if (w[0] == 0d) w[0] = 0.5;
                var order = new int[n];
                for (int i = 0; i < n; i++) order[i] = i;
                Array.Sort(order, (a, c) => { int v = x[a].CompareTo(x[c]); return v != 0 ? v : a.CompareTo(c); });
                var xs = new double[n];
                var ws = new double[n];
                for (int i = 0; i < n; i++) { xs[i] = x[order[i]]; ws[i] = w[order[i]]; }

                int m = random.Next(1, 40);
                var k = new double[m];
                for (int j = 0; j < m; j++)
                {
                    double u = random.NextDouble();
                    k[j] = u < 0.05 ? 0d : u > 0.95 ? 1d : u < 0.15 && j > 0 ? k[j - 1] : random.NextDouble();
                }
                var ascending = (double[])k.Clone();
                Array.Sort(ascending);

                foreach (var type in types)
                {
                    var unsorted = Stats.Percentile(x, k, w, type);
                    var sorted = Stats.Percentile(xs, k, ws, type, true);
                    var sortedAscending = Stats.Percentile(xs, ascending, ws, type, true);
                    var lists = Stats.Percentile(new System.Collections.Generic.List<double>(x), new System.Collections.Generic.List<double>(k), new System.Collections.Generic.List<double>(w), type);
                    for (int j = 0; j < m; j++)
                    {
                        long scalar = BitConverter.DoubleToInt64Bits(Stats.Percentile(x, k[j], w, type));
                        long scalarSorted = BitConverter.DoubleToInt64Bits(Stats.Percentile(xs, k[j], ws, type, true));
                        Assert.AreEqual(scalar, BitConverter.DoubleToInt64Bits(unsorted[j]), $"{type} unsorted trial {trial}, j = {j}");
                        Assert.AreEqual(scalar, BitConverter.DoubleToInt64Bits(lists[j]), $"{type} list trial {trial}, j = {j}");
                        Assert.AreEqual(scalarSorted, BitConverter.DoubleToInt64Bits(sorted[j]), $"{type} sorted trial {trial}, j = {j}");
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(Stats.Percentile(xs, ascending[j], ws, type, true)),
                            BitConverter.DoubleToInt64Bits(sortedAscending[j]), $"{type} ascending trial {trial}, j = {j}");
                    }
                }
            }
        }

        /// <summary>
        /// The Weibull weight formed from k and the bracketing positions, (k − r/d) / ((r + 1)/d − r/d), is closer to
        /// exact rational evaluation than h − floor(h) at the return-period probabilities 1 − 1/T, over seeded
        /// samples with n = 5 to 120 and T = 2 to 100. The exact value is evaluated in integers (every double is an
        /// integer multiple of 2^-1074), and the error is in units in the last place of the computed value.
        /// </summary>
        /// <param name="cases">The number of random cases.</param>
        public static void WeibullWeightCloserToExactRationalSweep(int cases)
        {
            var random = new Random(12345);
            var scale = BigInteger.Pow(2, 1074);
            double maxBracket = 0, maxRank = 0, sumBracket = 0, sumRank = 0;
            int count = 0;
            for (int c = 0; c < cases; c++)
            {
                int n = random.Next(5, 121);
                int T = random.Next(2, 101);
                double k = 1d - 1d / T;
                var x = new double[n];
                for (int i = 0; i < n; i++) x[i] = Math.Exp(3d + 2d * (random.NextDouble() - 0.5));
                Array.Sort(x);
                double h = (n + 1) * k;
                if (h <= 1d || h >= n) continue;
                int r = (int)Math.Floor(h);
                double rankWeightValue = x[r - 1] + (h - r) * (x[r] - x[r - 1]);
                double bracketWeightValue = Stats.Percentile(x, k, PPType.Weibull, true);
                // Exact: x(j) + ((n + 1) k - j)(x(j+1) - x(j)) with j the floor of the exact rank, scaled by 2^2148.
                var rank = (n + 1) * Scaled(k);
                int j = (int)BigInteger.Divide(rank, scale);
                if (j < 1 || j >= n) continue;
                var lo = Scaled(x[j - 1]);
                var exact = lo * scale + (rank - j * scale) * (Scaled(x[j]) - lo);
                double ulpBracket = UlpError(bracketWeightValue, exact, scale), ulpRank = UlpError(rankWeightValue, exact, scale);
                maxBracket = Math.Max(maxBracket, ulpBracket); maxRank = Math.Max(maxRank, ulpRank);
                sumBracket += ulpBracket; sumRank += ulpRank;
                count++;
            }
            Assert.IsGreaterThan(cases / 2, count);
            Assert.IsLessThan(sumRank / count, sumBracket / count, $"mean {sumBracket / count} vs {sumRank / count}");
            Assert.IsLessThanOrEqualTo(maxRank, maxBracket, $"max {maxBracket} vs {maxRank}");
        }

        /// <summary>
        /// The weighted Hazen percentile equals, bit for bit, the sequential reference below on random weighted sets with
        /// ties and zero weights, for sorted input and for unsorted input sorted stably by value; and the other formulas
        /// agree with the positions of
        /// <see cref="PlottingPositions.Function(System.Collections.Generic.IList{double}, PPType)"/>.
        /// </summary>
        /// <param name="cases">The number of random cases.</param>
        public static void WeightedHazenMatchesSequentialReferenceBitwiseSweep(int cases)
        {
            var random = new Random(77);
            for (int trial = 0; trial < cases; trial++)
            {
                int n = random.Next(1, 60);
                var x = new double[n];
                var w = new double[n];
                for (int i = 0; i < n; i++)
                {
                    x[i] = Math.Round(random.NextDouble() * 20d, random.Next(0, 3));
                    w[i] = random.NextDouble() < 0.1 ? 0d : random.NextDouble() * Math.Pow(10, random.Next(-6, 3));
                }
                if (w[0] == 0d) w[0] = 0.5;
                var order = new int[n];
                for (int i = 0; i < n; i++) order[i] = i;
                Array.Sort(order, (a, c) => { int v = x[a].CompareTo(x[c]); return v != 0 ? v : a.CompareTo(c); });
                var xs = new double[n];
                var ws = new double[n];
                for (int i = 0; i < n; i++) { xs[i] = x[order[i]]; ws[i] = w[order[i]]; }
                double k = random.NextDouble();
                if (k == 0d) continue;
                double expected = ReferenceWeightedHazen(xs, ws, k);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(Stats.Percentile(xs, k, ws, PPType.Hazen, true)), $"sorted trial {trial}");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(Stats.Percentile(x, k, w, PPType.Hazen)), $"unsorted trial {trial}");

                foreach (var type in new[] { PPType.Weibull, PPType.Median, PPType.Blom, PPType.Cunnane, PPType.Gringorten, PPType.KaplanMeier })
                {
                    var positions = PlottingPositions.Function(ws, type);
                    double value = Stats.Percentile(xs, k, ws, type, true);
                    Assert.AreEqual(Interpolate(xs, positions, k), value, 1e-9 * Math.Max(1d, Math.Abs(value)), $"{type} trial {trial}");
                }
            }
        }

        /// <summary>
        /// The sequential weighted Hazen reference: positions (c_(i−1) + w_i / 2) / W formed in a running loop, then linear
        /// interpolation between them, held at the end values.
        /// </summary>
        private static double ReferenceWeightedHazen(double[] values, double[] weights, double probability)
        {
            int count = values.Length;
            double totalWeight = 0.0;
            for (int i = 0; i < count; i++) totalWeight += weights[i];
            double cumulativeBefore = 0.0;
            double previousPosition = double.NaN;
            double previousValue = double.NaN;
            for (int i = 0; i < count; i++)
            {
                double position = (cumulativeBefore + 0.5 * weights[i]) / totalWeight;
                if (probability <= position)
                {
                    if (i == 0) return values[0];
                    double span = position - previousPosition;
                    if (span <= 0.0) return values[i];
                    return previousValue + (values[i] - previousValue) * (probability - previousPosition) / span;
                }
                previousPosition = position;
                previousValue = values[i];
                cumulativeBefore += weights[i];
            }
            return values[count - 1];
        }

        /// <summary>Linear interpolation of (position, value) points, held at the ends.</summary>
        private static double Interpolate(double[] values, double[] positions, double k)
        {
            for (int i = 0; i < values.Length; i++)
            {
                if (k <= positions[i])
                {
                    if (i == 0) return values[0];
                    double span = positions[i] - positions[i - 1];
                    if (span <= 0d) return values[i];
                    return values[i - 1] + (values[i] - values[i - 1]) * (k - positions[i - 1]) / span;
                }
            }
            return values[values.Length - 1];
        }

        /// <summary>Returns a finite non-negative double as an integer multiple of 2^-1074.</summary>
        private static BigInteger Scaled(double value)
        {
            long bits = BitConverter.DoubleToInt64Bits(value);
            int exponent = (int)((bits >> 52) & 0x7FF);
            long mantissa = bits & 0xFFFFFFFFFFFFFL;
            if (exponent == 0) return mantissa;
            return new BigInteger(mantissa | (1L << 52)) << (exponent - 1);
        }

        /// <summary>Returns |value - exact| in units in the last place of value, with exact scaled by 2^-2148.</summary>
        private static double UlpError(double value, BigInteger exact, BigInteger scale)
        {
            var diff = BigInteger.Abs(Scaled(value) * scale - exact);
            int exponent = (int)((BitConverter.DoubleToInt64Bits(value) >> 52) & 0x7FF);
            var ulp = (exponent == 0 ? BigInteger.One : BigInteger.Pow(2, exponent - 1)) * scale;
            return (double)(diff * 1000 / ulp) / 1000d;
        }
    }
}
