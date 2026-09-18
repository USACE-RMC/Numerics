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

    }
}
