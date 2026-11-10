#nullable disable
using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;

namespace Utilities
{
    /// <summary>
    /// The reference formulas and the bitwise comparison grids of <see cref="Test_ToolsBitwiseReference"/>, shared by the
    /// unit tests and the full sweep in Benchmarks_Numerics.
    /// </summary>
    public partial class Test_ToolsBitwiseReference
    {
        private const double DoubleMachineEpsilon = 1.11022302462516E-16;

        /// <summary>Reference formula for log(1 + x).</summary>
        private static double ReferenceLog1p(double x)
        {
            double y = 1.0 + x;
            if (!(y > 0d) || double.IsInfinity(y))
                return Math.Log(y);
            return Math.Log(y) - ((y - 1.0) - x) / y;
        }

        /// <summary>Reference formula for exp(x) - 1: the direct difference for |x| &gt; 1/2, otherwise the compensated form.</summary>
        private static double ReferenceExpm1(double x)
        {
            return Math.Abs(x) > 0.5 ? Math.Exp(x) - 1.0 : CompensatedExpm1(x);
        }

        /// <summary>The compensated formula for exp(x) - 1, used for |x| &lt;= 1/2 and for NaN.</summary>
        private static double CompensatedExpm1(double x)
        {
            double u = Math.Exp(x);
            if (u == 1.0) return x;
            if (double.IsPositiveInfinity(u)) return u;
            if (u <= DoubleMachineEpsilon / 2.0) return -1.0;
            double numerator = (u - 1.0) * x;
            if (double.IsInfinity(numerator)) return u - 1.0;
            return numerator / Math.Log(u);
        }

        /// <summary>Reference formula for log(1 - exp(a)).</summary>
        private static double ReferenceLog1mExp(double a)
        {
            if (a > 0 || double.IsNaN(a)) return double.NaN;
            return a < -0.69314718055994530942 ? ReferenceLog1p(-Math.Exp(a)) : Math.Log(-ReferenceExpm1(a));
        }

        /// <summary>Returns the double <paramref name="ulps"/> bit patterns away from <paramref name="x"/>.</summary>
        private static double Step(double x, long ulps) => BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(x) + ulps);

        /// <summary>Builds the deterministic part of the grid, adding walks across the given branch boundaries.</summary>
        /// <param name="perDecade">The number of log-spaced magnitudes per decade.</param>
        /// <param name="walk">The number of bit patterns walked on each side of each boundary.</param>
        /// <param name="boundaries">The branch boundaries.</param>
        private static List<double> Grid(int perDecade, int walk, params double[] boundaries)
        {
            var points = new List<double>();
            // Log-spaced magnitudes from 1e-320 to 1e308, both signs.
            for (int k = -320 * perDecade; k <= 308 * perDecade; k++)
            {
                double m = Math.Pow(10d, k / (double)perDecade);
                if (double.IsInfinity(m)) continue;
                points.Add(m);
                points.Add(-m);
            }
            // Signed zeros, infinities, extremes and subnormals.
            points.AddRange(new[]
            {
                0d, -0d, double.PositiveInfinity, double.NegativeInfinity, double.MaxValue, -double.MaxValue,
                double.Epsilon, -double.Epsilon, 2.2250738585072014E-308, -2.2250738585072014E-308,
            });
            for (long b = 1; b < 0x0010000000000000L; b = b * 3 + 1)
            {
                points.Add(BitConverter.Int64BitsToDouble(b));
                points.Add(-BitConverter.Int64BitsToDouble(b));
            }
            points.Add(BitConverter.Int64BitsToDouble(0x000FFFFFFFFFFFFFL));
            points.Add(-BitConverter.Int64BitsToDouble(0x000FFFFFFFFFFFFFL));
            // NaNs: the default NaN, both signs, quiet and signaling patterns with payloads.
            foreach (ulong bits in new ulong[]
            {
                0x7FF8000000000000UL, 0xFFF8000000000000UL, 0x7FF8000000000001UL, 0xFFF8000000000001UL,
                0x7FFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0x7FF0000000000001UL, 0xFFF0000000000001UL,
                0x7FF4000000000000UL, 0x7FF800000000ABCDUL, 0xFFF0DEADBEEF0001UL, 0x7FF7FFFFFFFFFFFFUL,
            })
            {
                points.Add(BitConverter.Int64BitsToDouble(unchecked((long)bits)));
            }
            points.Add(double.NaN);
            // Dense walks across each boundary: every bit pattern within the walk, then relative offsets out to 1e-2.
            foreach (double c in boundaries)
            {
                for (long k = -walk; k <= walk; k++) points.Add(Step(c, k));
                for (int e = -160; e <= -20; e++)
                {
                    double d = Math.Pow(10d, e / 10d);
                    points.Add(c * (1d + d));
                    points.Add(c * (1d - d));
                }
            }
            return points;
        }

        /// <summary>Compares a library method with its reference over the grid and random inputs, bit for bit.</summary>
        private static void AssertBitwise(Func<double, double> actual, Func<double, double> reference, List<double> grid,
            double low, double high, int seed, int randomCount, string name)
        {
            var random = new Random(seed);
            var buffer = new byte[8];
            int count = 0, diffs = 0;
            string first = null;
            void Check(double x)
            {
                count++;
                long a = BitConverter.DoubleToInt64Bits(actual(x)), r = BitConverter.DoubleToInt64Bits(reference(x));
                if (a == r) return;
                diffs++;
                first ??= $"{name}(0x{BitConverter.DoubleToInt64Bits(x):X16}) = 0x{a:X16}, reference 0x{r:X16}";
            }
            foreach (double x in grid) Check(x);
            for (int i = 0; i < randomCount; i++)
            {
                random.NextBytes(buffer);
                Check(BitConverter.ToDouble(buffer, 0));
                Check(low + (high - low) * random.NextDouble());
            }
            Assert.AreEqual(grid.Count + 2 * randomCount, count);
            Assert.AreEqual(0, diffs, $"{diffs} of {count} results differ; first: {first}");
        }

        /// <summary>Expm1 matches its reference bit for bit.</summary>
        /// <param name="perDecade">The number of log-spaced magnitudes per decade.</param>
        /// <param name="walk">The number of bit patterns walked on each side of each branch boundary.</param>
        /// <param name="randomCount">The number of random bit patterns, and of random values in the working range.</param>
        public static void Expm1MatchesReferenceBitwiseSweep(int perDecade, int walk, int randomCount)
        {
            // u = 1 (|x| near 2^-53), u at 2^-54 (x near -37.43), exp overflow (x near 709.78), the former
            // compensating product's overflow (x near 703.2), the switch at |x| = 1/2, and |x| = 1.
            var grid = Grid(perDecade, walk, 1.1102230246251565E-16, -1.1102230246251565E-16, 5.551115123125783E-17, -5.551115123125783E-17,
                -37.42994775023705, 709.782712893384, 703.2, 703.25, 0.5, -0.5, 1d, -1d);
            AssertBitwise(x => Tools.Expm1(x), ReferenceExpm1, grid, -40d, 40d, 101, randomCount, "Expm1");
        }

        /// <summary>Log1p matches its reference bit for bit.</summary>
        /// <param name="perDecade">The number of log-spaced magnitudes per decade.</param>
        /// <param name="walk">The number of bit patterns walked on each side of each branch boundary.</param>
        /// <param name="randomCount">The number of random bit patterns, and of random values in the working range.</param>
        public static void Log1pMatchesReferenceBitwiseSweep(int perDecade, int walk, int randomCount)
        {
            // 1 + x at 0 (x = -1), at the overflow to infinity, |x| = 1/2, |x| = 1, and where 1 + x rounds to 1.
            var grid = Grid(perDecade, walk, -1d, double.MaxValue, 0.5, -0.5, 1.1102230246251565E-16, -1.1102230246251565E-16, 1d);
            AssertBitwise(x => Tools.Log1p(x), ReferenceLog1p, grid, -1d, 2d, 202, randomCount, "Log1p");
        }

        /// <summary>Log1mExp matches its reference bit for bit.</summary>
        /// <param name="perDecade">The number of log-spaced magnitudes per decade.</param>
        /// <param name="walk">The number of bit patterns walked on each side of each branch boundary.</param>
        /// <param name="randomCount">The number of random bit patterns, and of random values in the working range.</param>
        public static void Log1mExpMatchesReferenceBitwiseSweep(int perDecade, int walk, int randomCount)
        {
            // The -ln 2 switch, a = 0 from below, the -1/2 point, the inner Expm1 thresholds, and where exp(a) underflows.
            var grid = Grid(perDecade, walk, -0.69314718055994530942, -0.5, -1.1102230246251565E-16, -37.42994775023705, -745.1332191019411);
            AssertBitwise(a => LibraryLog1mExp(a), ReferenceLog1mExp, grid, -40d, 0d, 303, randomCount, "Log1mExp");
        }
    }
}
