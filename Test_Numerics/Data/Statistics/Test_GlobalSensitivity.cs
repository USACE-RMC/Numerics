using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Sampling;
using System;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the given-data global sensitivity estimators.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_GlobalSensitivity
    {

        /// <summary>
        /// Generates the fixed-seed Ishigami sample: y = sin(x1) + a*sin(x2)^2 + b*x3^4*sin(x1)
        /// with the inputs uniform on (-pi, pi).
        /// </summary>
        /// <param name="n">The sample size.</param>
        /// <param name="x1">Output. The first input sample.</param>
        /// <param name="x2">Output. The second input sample.</param>
        /// <param name="x3">Output. The third input sample.</param>
        /// <param name="y">Output. The output sample.</param>
        private static void IshigamiSample(int n, out double[] x1, out double[] x2, out double[] x3, out double[] y)
        {
            const double a = 7d, b = 0.1d;
            var prng = new MersenneTwister(12345);
            x1 = new double[n];
            x2 = new double[n];
            x3 = new double[n];
            y = new double[n];
            for (int i = 0; i < n; i++)
            {
                x1[i] = -Math.PI + 2d * Math.PI * prng.NextDouble();
                x2[i] = -Math.PI + 2d * Math.PI * prng.NextDouble();
                x3[i] = -Math.PI + 2d * Math.PI * prng.NextDouble();
                y[i] = Math.Sin(x1[i]) + a * Math.Sin(x2[i]) * Math.Sin(x2[i]) + b * Math.Pow(x3[i], 4) * Math.Sin(x1[i]);
            }
        }

        /// <summary>
        /// The given-data first-order Sobol estimator recovers the analytic Ishigami indices
        /// (a = 7, b = 0.1): S1 = V1/V and S2 = V2/V from the closed-form variance decomposition,
        /// and S3 = 0 (the third input acts only through interaction). The 0.03 tolerance covers
        /// the estimator's binning bias and the fixed-seed Monte Carlo noise at 2^14 samples.
        /// </summary>
        [TestMethod]
        public void Test_FirstOrderSobol_IshigamiOracle()
        {
            const double a = 7d, b = 0.1d;
            double pi4 = Math.Pow(Math.PI, 4);
            double v1 = 0.5d * Math.Pow(1d + b * pi4 / 5d, 2);
            double v2 = a * a / 8d;
            double v = 0.5d + a * a / 8d + b * pi4 / 5d + b * b * Math.Pow(Math.PI, 8) / 18d;

            IshigamiSample(16384, out var x1, out var x2, out var x3, out var y);
            Assert.AreEqual(v1 / v, GlobalSensitivity.FirstOrderSobol(x1, y), 0.03d);
            Assert.AreEqual(v2 / v, GlobalSensitivity.FirstOrderSobol(x2, y), 0.03d);
            Assert.IsLessThan(0.02d, GlobalSensitivity.FirstOrderSobol(x3, y));
        }

        /// <summary>
        /// The given-data first-order Sobol estimator recovers the analytic Sobol-g indices for
        /// a = (0, 1, 4.5, 9): V(i) = 1/(3*(1+a(i))^2) and V = prod(1+V(i)) - 1.
        /// </summary>
        [TestMethod]
        public void Test_FirstOrderSobol_SobolGOracle()
        {
            var a = new double[] { 0d, 1d, 4.5d, 9d };
            int n = 16384;
            var prng = new MersenneTwister(45678);
            var xs = new double[4][];
            for (int j = 0; j < 4; j++) xs[j] = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                double product = 1d;
                for (int j = 0; j < 4; j++)
                {
                    double u = prng.NextDouble();
                    xs[j][i] = u;
                    product *= (Math.Abs(4d * u - 2d) + a[j]) / (1d + a[j]);
                }
                y[i] = product;
            }

            double v = 1d;
            var vi = new double[4];
            for (int j = 0; j < 4; j++)
            {
                vi[j] = 1d / (3d * (1d + a[j]) * (1d + a[j]));
                v *= 1d + vi[j];
            }
            v -= 1d;
            for (int j = 0; j < 4; j++)
            {
                Assert.AreEqual(vi[j] / v, GlobalSensitivity.FirstOrderSobol(xs[j], y), 0.03d, $"input {j}");
            }
        }

        /// <summary>
        /// An independent input reads near zero on all three estimators, at each estimator's own
        /// finite-sample bias scale: the Sobol index at the bin-count-over-sample-size scale, the
        /// PAWN median at the two-sample Kolmogorov-Smirnov noise scale, and the double-histogram
        /// delta at the class-fluctuation scale (the largest of the three by construction).
        /// </summary>
        [TestMethod]
        public void Test_IndependentInput_NearZero()
        {
            int n = 16384;
            var prng = new MersenneTwister(999);
            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = prng.NextDouble();
                y[i] = Normal.StandardZ(prng.NextDouble());
            }
            Assert.IsLessThan(0.02d, GlobalSensitivity.FirstOrderSobol(x, y));
            Assert.IsLessThan(0.06d, GlobalSensitivity.PawnMedian(x, y));
            Assert.IsLessThan(0.1d, GlobalSensitivity.BorgonovoDelta(x, y));
        }

        /// <summary>
        /// A functionally dependent output saturates the estimators: on the exact grid y = x with
        /// the sample size divisible by the bin counts, the Sobol index is essentially one, every
        /// conditional distribution is far from the unconditional one, and the double-histogram
        /// delta hits its exact saturation value 1 - 1/yBins because each input bin maps onto
        /// exactly one output class.
        /// </summary>
        [TestMethod]
        public void Test_FunctionalDependence_Saturates()
        {
            int n = 1000;
            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = i + 0.5d;
                y[i] = x[i];
            }
            Assert.IsGreaterThan(0.99d, GlobalSensitivity.FirstOrderSobol(x, y));
            Assert.IsGreaterThanOrEqualTo(0.45d, GlobalSensitivity.PawnMedian(x, y));
            Assert.AreEqual(1d - 1d / 20d, GlobalSensitivity.BorgonovoDelta(x, y), 1E-12);
        }

        /// <summary>
        /// The Borgonovo delta is exactly invariant to strictly increasing transforms of the
        /// output - the class partition is rank-based, so cubing the Ishigami output changes no
        /// class assignment and the estimate is bit-identical.
        /// </summary>
        [TestMethod]
        public void Test_BorgonovoDelta_MonotoneTransformInvariant()
        {
            IshigamiSample(4096, out var x1, out _, out _, out var y);
            var cubed = new double[y.Length];
            for (int i = 0; i < y.Length; i++)
                cubed[i] = y[i] * y[i] * y[i];
            Assert.AreEqual(GlobalSensitivity.BorgonovoDelta(x1, y), GlobalSensitivity.BorgonovoDelta(x1, cubed), 0d);
        }

        /// <summary>
        /// The estimators are deterministic: repeated calls on the same sample are bit-identical.
        /// </summary>
        [TestMethod]
        public void Test_Determinism_RepeatedCallsBitEqual()
        {
            IshigamiSample(4096, out var x1, out _, out _, out var y);
            Assert.AreEqual(GlobalSensitivity.FirstOrderSobol(x1, y), GlobalSensitivity.FirstOrderSobol(x1, y), 0d);
            Assert.AreEqual(GlobalSensitivity.PawnMedian(x1, y), GlobalSensitivity.PawnMedian(x1, y), 0d);
            Assert.AreEqual(GlobalSensitivity.BorgonovoDelta(x1, y), GlobalSensitivity.BorgonovoDelta(x1, y), 0d);
            var first = GlobalSensitivity.Pawn(x1, y);
            var second = GlobalSensitivity.Pawn(x1, y);
            for (int b = 0; b < first.Length; b++)
                Assert.AreEqual(first[b], second[b], 0d);
        }

        /// <summary>
        /// Regression pin: the PAWN statistics and median of every Ishigami input at 2^12 samples equal, bit
        /// for bit, the recorded values. Each per-bin statistic is the correctly rounded rational D/(nb n)
        /// for its bin size nb.
        /// </summary>
        [TestMethod]
        public void Test_Pawn_BitwiseRegression()
        {
            IshigamiSample(4096, out var x1, out var x2, out var x3, out var y);
            var bits = new System.Text.StringBuilder();
            foreach (var x in new[] { x1, x2, x3 })
            {
                foreach (double value in GlobalSensitivity.Pawn(x, y))
                    bits.Append(BitConverter.DoubleToInt64Bits(value).ToString("X16"));
                bits.Append(BitConverter.DoubleToInt64Bits(GlobalSensitivity.PawnMedian(x, y)).ToString("X16"));
            }
            Assert.AreEqual(PawnRegressionBits, bits.ToString());
        }

        /// <summary>
        /// The concatenated big-endian bit patterns pinned by <see cref="Test_Pawn_BitwiseRegression"/>.
        /// </summary>
        private const string PawnRegressionBits =
            "3FC4A800000000003FCA7800000000003FCCF98F9C18F9C23FD472A2576A25773FD290D1CB8D1CB93FD4CA82828282833FD1C081DF881DF9"
            + "3FCC026B6526B6523FC8F3BF103BF1043FC5084FEC04FEC03FC60000000000003FCB3000000000003FCD7800000000003FCF3027F6027F60"
            + "3FD563562A7562A73FD42FC3C3C3C3C43FD104F9C18F9C193FCBB7C40EFC40F03FC8B000000000003FC56800000000003FCC7DFD809FD80A"
            + "3FDF85B9B9B9B9BA3FD6D2CF4C2CF4C33FCA4C86DE486DE53FD5521B7921B7923FDD2ECF4C2CF4C33FE04F23232323233FDC41B7921B7922"
            + "3FCA44AED44AED453FD5527F6027F6023FE0AAE5C68E5C693FDFCC14141414143FD77B1F3831F3833FC42CAED44AED453FD772436F2436F2"
            + "3FDF8AE3472E34733FE11214141414143FD8626B6526B6523FCA2E167A6167A63FD806077E2077E23FE07721B7921B793FD8343971A3971A"
            + "3FCF3969696969693FBFD103BF103BF13FBCD2E3472E34733FB9FF4C2CF4C2CF3FB40C7CE0C7CE0C3FB18000000000003FBC690DBC90DBC9"
            + "3FBBB31F3831F3833FB8309FD809FD813FB9040EFC40EFC43FBA5373737373733FBAC0EFC40EFC413FBB80EFC40EFC413FB5DA9D58A9D58B"
            + "3FB4945EE845EE843FB84414141414143FBE70DBC90DBC913FB2A063E7063E703FC0A1999999999A3FD10E895DA895DB3FBA8A319BC137DA";

        /// <summary>
        /// Each per-bin PAWN statistic equals the two-sample Kolmogorov-Smirnov statistic of the sorted bin and
        /// the sorted full output sample.
        /// </summary>
        [TestMethod]
        public void Test_Pawn_BinsMatchTwoSampleKolmogorovSmirnov()
        {
            IshigamiSample(1000, out var x1, out _, out _, out var y);
            const int bins = 10;
            var pawn = GlobalSensitivity.Pawn(x1, y, bins);
            var order = new int[y.Length];
            for (int i = 0; i < order.Length; i++) order[i] = i;
            Array.Sort(order, (first, second) => { int c = x1[first].CompareTo(x1[second]); return c != 0 ? c : first.CompareTo(second); });
            var sortedY = (double[])y.Clone();
            Array.Sort(sortedY);
            for (int b = 0; b < bins; b++)
            {
                int start = b * y.Length / bins, end = (b + 1) * y.Length / bins;
                var bin = new double[end - start];
                for (int i = start; i < end; i++) bin[i - start] = y[order[i]];
                Array.Sort(bin);
                Assert.AreEqual(pawn[b], HypothesisTests.TwoSampleKolmogorovSmirnov(bin, sortedY), 0d);
            }
        }

        /// <summary>
        /// Equal rank values are resolved by their original sample indices, so a tied input with
        /// outputs split by original position forms the same two pure bins on every framework.
        /// </summary>
        [TestMethod]
        public void Test_TiedRanks_UseOriginalIndexOrder()
        {
            var x = new double[32];
            var y = new double[32];
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = 1d;
                y[i] = i < x.Length / 2 ? 0d : 1d;
            }

            Assert.AreEqual(1d, GlobalSensitivity.FirstOrderSobol(x, y, 2), 0d);
            CollectionAssert.AreEqual(new double[] { 0.5d, 0.5d }, GlobalSensitivity.Pawn(x, y, 2));
            Assert.AreEqual(0.5d, GlobalSensitivity.BorgonovoDelta(x, y, 2, 2), 0d);
        }

        /// <summary>
        /// Guard matrix: nulls, length mismatches, degenerate bin counts, samples shorter than the
        /// bin counts, and non-finite values all throw the documented exceptions, and a
        /// zero-variance output returns a zero Sobol index.
        /// </summary>
        [TestMethod]
        public void Test_GuardMatrix()
        {
            var x = new double[] { 1d, 2d, 3d, 4d };
            var y = new double[] { 5d, 6d, 7d, 8d };
            Assert.Throws<ArgumentNullException>(() => GlobalSensitivity.FirstOrderSobol(null!, y, 2));
            Assert.Throws<ArgumentNullException>(() => GlobalSensitivity.FirstOrderSobol(x, null!, 2));
            Assert.Throws<ArgumentException>(() => GlobalSensitivity.FirstOrderSobol(x, new double[] { 1d }, 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => GlobalSensitivity.FirstOrderSobol(x, y, 1));
            Assert.Throws<ArgumentException>(() => GlobalSensitivity.FirstOrderSobol(x, y, 5));
            // The finiteness guards name the sample the offending value came from.
            var nonFiniteInput = Assert.Throws<ArgumentOutOfRangeException>(() => GlobalSensitivity.FirstOrderSobol(new double[] { 1d, double.NaN, 3d, 4d }, y, 2));
            Assert.AreEqual("x", nonFiniteInput.ParamName);
            var nonFiniteOutput = Assert.Throws<ArgumentOutOfRangeException>(() => GlobalSensitivity.Pawn(x, new double[] { 1d, 2d, double.PositiveInfinity, 4d }, 2));
            Assert.AreEqual("y", nonFiniteOutput.ParamName);
            Assert.Throws<ArgumentOutOfRangeException>(() => GlobalSensitivity.BorgonovoDelta(x, y, 2, 1));
            Assert.Throws<ArgumentException>(() => GlobalSensitivity.BorgonovoDelta(x, y, 2, 5));
            Assert.AreEqual(0d, GlobalSensitivity.FirstOrderSobol(x, new double[] { 3d, 3d, 3d, 3d }, 2), 0d);
        }

        /// <summary>
        /// Tied outputs share one class. A constant output carries no information about any input (analytic delta
        /// 0; the first-order Sobol index agrees), and y = x mod 2 over x = 0..199 puts five even and five odd values
        /// in every one of the 20 input bins, so each conditional distribution equals the marginal (analytic delta 0).
        /// </summary>
        [TestMethod]
        public void Test_BorgonovoDelta_TiedOutputs()
        {
            var x = new double[200];
            var constant = new double[200];
            var parity = new double[200];
            for (int i = 0; i < 200; i++)
            {
                x[i] = i;
                constant[i] = 3d;
                parity[i] = i % 2;
            }
            Assert.AreEqual(0d, GlobalSensitivity.BorgonovoDelta(x, constant), 0d);
            Assert.AreEqual(0d, GlobalSensitivity.BorgonovoDelta(x, parity), 1E-15);
            Assert.AreEqual(0d, GlobalSensitivity.FirstOrderSobol(x, parity), 1E-12);
        }

        /// <summary>
        /// Untied outputs are classed exactly as without the tie rule: the delta equals, bit for bit, the
        /// equal-frequency rank classes of the 2.2.0 estimator, kept as a bitwise reference.
        /// </summary>
        [TestMethod]
        public void Test_BorgonovoDelta_UntiedUnchanged()
        {
            var rng = new Random(77);
            var x = new double[1000];
            var y = new double[1000];
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = rng.NextDouble();
                y[i] = Math.Sin(3 * x[i]) + 0.3 * rng.NextDouble();
            }
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(LegacyBorgonovoDelta(x, y, 20, 20)),
                BitConverter.DoubleToInt64Bits(GlobalSensitivity.BorgonovoDelta(x, y)));
        }

        /// <summary>
        /// The 2.2.0 estimator (output classes split ties by sample index), a bitwise reference for untied outputs.
        /// </summary>
        private static double LegacyBorgonovoDelta(double[] x, double[] y, int xBins, int yBins)
        {
            int n = x.Length;
            int[] Order(double[] v)
            {
                var o = new int[n];
                for (int i = 0; i < n; i++) o[i] = i;
                Array.Sort(o, (a, b) => { int c = v[a].CompareTo(v[b]); return c != 0 ? c : a.CompareTo(b); });
                return o;
            }
            var xOrder = Order(x);
            var yOrder = Order(y);
            var yClass = new int[n];
            var classCounts = new int[yBins];
            for (int c = 0; c < yBins; c++)
            {
                int start = (int)((long)c * n / yBins), end = (int)((long)(c + 1) * n / yBins);
                classCounts[c] = end - start;
                for (int i = start; i < end; i++) yClass[yOrder[i]] = c;
            }
            double delta = 0;
            var counts = new int[yBins];
            for (int b = 0; b < xBins; b++)
            {
                int start = (int)((long)b * n / xBins), end = (int)((long)(b + 1) * n / xBins), size = end - start;
                Array.Clear(counts, 0, yBins);
                for (int i = start; i < end; i++) counts[yClass[xOrder[i]]]++;
                double tv = 0;
                for (int c = 0; c < yBins; c++) tv += Math.Abs((double)counts[c] / size - (double)classCounts[c] / n);
                delta += (double)size / n * tv;
            }
            return 0.5d * delta;
        }
    }
}
