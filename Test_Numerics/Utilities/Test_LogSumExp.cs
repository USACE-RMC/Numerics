using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;

namespace Utilities
{
    /// <summary>
    /// Tests for <see cref="Tools.LogSumExp(double, double)"/>, <see cref="Tools.LogSumExp(IList{double})"/> and
    /// <see cref="Tools.LogSumExp(double[,], int)"/>.
    /// </summary>
    /// <remarks>
    /// The expected values are <c>scipy.special.logsumexp</c> results from SciPy 1.18.1 with NumPy 2.5.4 (Python 3.14),
    /// printed with <c>repr</c>. Agreement is required to 2 ulps of max(|result|, |largest input|), the scale at which
    /// the final addition of the largest input is rounded.
    /// </remarks>
    [TestClass]
    public class Test_LogSumExp
    {
        /// <summary>Asserts agreement to a number of ulps of max(|expected|, |scale|).</summary>
        private static void AssertUlps(double expected, double actual, double scale, double maxUlps = 2)
        {
            double reference = Math.Max(Math.Abs(expected), Math.Abs(scale));
            double ulp = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(reference) + 1) - reference;
            Assert.IsLessThanOrEqualTo(maxUlps * ulp, Math.Abs(expected - actual),
                $"expected {expected:R}, actual {actual:R}, {Math.Abs(expected - actual) / ulp} ulps");
        }

        /// <summary>The largest value of a list, the scale at which its result is rounded.</summary>
        private static double Largest(IList<double> values)
        {
            double max = double.NegativeInfinity;
            foreach (double v in values) if (v > max) max = v;
            return max;
        }

        /// <summary>
        /// Lists with ties, results near zero, a dominant term and the edges of the exponential range agree with SciPy.
        /// </summary>
        [TestMethod]
        public void List_MatchesScipy()
        {
            var cases = new (double[] Values, double Expected)[]
            {
                (new[] { 1000.01, 1000.02 }, 1000.7081596805078),
                (new[] { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 }, 9.45862974442671),
                (new[] { 0.0, -30.0 }, 9.357622968839737e-14),
                (new[] { -0.5, -1.2, -2.0, -0.5, -0.5 }, 0.8136471771714695),
                (new[] { -1.3862943611198906, -1.3862943611198906, -0.6931471805599453 }, 0.0),
                (new[] { 1e-20, -50.0, -3.0 }, 0.048587351573742055),
                (new[] { -745.0, -700.0, 709.0, 708.5, -1e+300 }, 709.4740769841801),
                (Random64, 56.932863328731095),
            };
            foreach (var (values, expected) in cases)
            {
                AssertUlps(expected, Tools.LogSumExp(values), Largest(values));
                AssertUlps(expected, Tools.LogSumExp(new List<double>(values)), Largest(values));
            }
        }

        /// <summary>
        /// A result near zero keeps its relative accuracy: log(1 + exp(-30)) is 9.357622968839737e-14 to the last bit,
        /// where max + log(exp(0) + exp(-30)) gives 9.348077867343381e-14.
        /// </summary>
        [TestMethod]
        public void NearZeroResult_IsRelativelyAccurate()
        {
            Assert.AreEqual(9.357622968839737e-14, Tools.LogSumExp(0.0, -30.0), 0d);
            Assert.AreEqual(9.357622968839737e-14, Tools.LogSumExp(new[] { 0.0, -30.0 }), 0d);
        }

        /// <summary>Pairs agree with SciPy and are bitwise the list form of the same two values.</summary>
        [TestMethod]
        public void Pair_MatchesScipyAndList()
        {
            var cases = new (double U, double V, double Expected)[]
            {
                (0.0, -30.0, 9.357622968839737e-14),
                (-0.001, -0.002, 0.6916473055599401),
                (3.0, 3.0, 3.6931471805599454),
                (-5.25, 2.5, 2.500430649797639),
                (1000.01, 1000.02, 1000.7081596805078),
            };
            foreach (var (u, v, expected) in cases)
            {
                double pair = Tools.LogSumExp(u, v);
                AssertUlps(expected, pair, Math.Max(u, v));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.LogSumExp(new[] { u, v })), BitConverter.DoubleToInt64Bits(pair));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(pair), BitConverter.DoubleToInt64Bits(Tools.LogSumExp(v, u)));
            }
        }

        /// <summary>
        /// The special values follow SciPy: +∞ gives +∞, NaN propagates (also beside +∞), every value −∞ gives −∞, an
        /// empty list gives −∞ (the log of an empty sum), and −∞ terms add nothing.
        /// </summary>
        [TestMethod]
        public void SpecialValues_FollowScipy()
        {
            double inf = double.PositiveInfinity, ninf = double.NegativeInfinity, nan = double.NaN;
            Assert.AreEqual(inf, Tools.LogSumExp(inf, 5d));
            Assert.AreEqual(inf, Tools.LogSumExp(5d, inf));
            Assert.AreEqual(inf, Tools.LogSumExp(inf, ninf));
            Assert.AreEqual(inf, Tools.LogSumExp(inf, inf));
            Assert.IsTrue(double.IsNaN(Tools.LogSumExp(nan, inf)));
            Assert.IsTrue(double.IsNaN(Tools.LogSumExp(inf, nan)));
            Assert.AreEqual(ninf, Tools.LogSumExp(ninf, ninf));
            Assert.AreEqual(5d, Tools.LogSumExp(ninf, 5d));

            Assert.AreEqual(inf, Tools.LogSumExp(new[] { 1d, inf, -2d }));
            Assert.AreEqual(inf, Tools.LogSumExp(new[] { ninf, inf }));
            Assert.IsTrue(double.IsNaN(Tools.LogSumExp(new[] { 1d, inf, nan })));
            Assert.IsTrue(double.IsNaN(Tools.LogSumExp(new[] { nan, 1d })));
            Assert.AreEqual(ninf, Tools.LogSumExp(new[] { ninf, ninf, ninf }));
            Assert.AreEqual(ninf, Tools.LogSumExp(Array.Empty<double>()));
            Assert.AreEqual(ninf, Tools.LogSumExp(new List<double>()));
            Assert.AreEqual(2d, Tools.LogSumExp(new[] { ninf, 2d, ninf }), 0d);
            Assert.ThrowsExactly<ArgumentNullException>(() => Tools.LogSumExp((IList<double>)null!));
        }

        /// <summary>
        /// Both axes of a matrix agree with SciPy's <c>logsumexp(m, axis=...)</c> and are bitwise the list form of each
        /// column or row, including a tie in a column.
        /// </summary>
        [TestMethod]
        public void Matrix_MatchesScipyAndList()
        {
            var m = new double[,]
            {
                { -1.612921, -0.452752, 1.056098, 1.128163 },
                { 0.288056, 1.386014, 1.056098, -2.076967 },
                { 1.78159, 3.463898, -0.622853, -6.090996 },
            };
            var axis0 = new[] { 2.011220111055738, 3.5993921686132095, 1.838431916161694, 1.168619533296508 };
            var axis1 = new[] { 1.9169646288318463, 2.1202416229067627, 3.648557556827851 };
            var columns = Tools.LogSumExp(m, 0);
            var rows = Tools.LogSumExp(m, 1);
            Assert.HasCount(4, columns);
            Assert.HasCount(3, rows);
            for (int c = 0; c < 4; c++)
            {
                var column = new[] { m[0, c], m[1, c], m[2, c] };
                AssertUlps(axis0[c], columns[c], Largest(column));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.LogSumExp(column)), BitConverter.DoubleToInt64Bits(columns[c]));
            }
            for (int r = 0; r < 3; r++)
            {
                var row = new[] { m[r, 0], m[r, 1], m[r, 2], m[r, 3] };
                AssertUlps(axis1[r], rows[r], Largest(row));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.LogSumExp(row)), BitConverter.DoubleToInt64Bits(rows[r]));
            }
        }

        /// <summary>
        /// Each column and row of a random matrix with special values is bitwise the list form, so the special values of
        /// the matrix form are those of <see cref="Tools.LogSumExp(IList{double})"/>, including the NaN bit pattern.
        /// </summary>
        [TestMethod]
        public void Matrix_SpecialValuesBitwiseList()
        {
            var random = new Random(20261108);
            var m = new double[7, 5];
            for (int r = 0; r < 7; r++)
                for (int c = 0; c < 5; c++)
                    m[r, c] = 40 * random.NextDouble() - 20;
            m[2, 0] = BitConverter.Int64BitsToDouble(0x7FF8000000000000); // a NaN with the sign bit clear, as NumPy produces
            m[4, 1] = double.PositiveInfinity;
            m[1, 1] = double.NaN;
            for (int r = 0; r < 7; r++) m[r, 2] = double.NegativeInfinity;
            m[3, 3] = m[5, 3] = 25;
            m[6, 4] = double.PositiveInfinity;
            var columns = Tools.LogSumExp(m, 0);
            var rows = Tools.LogSumExp(m, 1);
            for (int c = 0; c < 5; c++)
            {
                var column = new double[7];
                for (int r = 0; r < 7; r++) column[r] = m[r, c];
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.LogSumExp(column)), BitConverter.DoubleToInt64Bits(columns[c]));
            }
            for (int r = 0; r < 7; r++)
            {
                var row = new double[5];
                for (int c = 0; c < 5; c++) row[c] = m[r, c];
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.LogSumExp(row)), BitConverter.DoubleToInt64Bits(rows[r]));
            }
            Assert.IsTrue(double.IsNaN(columns[0]) && double.IsNaN(columns[1]));
            Assert.AreEqual(double.NegativeInfinity, columns[2]);
            Assert.AreEqual(double.PositiveInfinity, columns[4]);
        }

        /// <summary>
        /// An empty reduced dimension gives −∞ per entry and an empty kept dimension an empty result; a null matrix or an
        /// axis other than 0 or 1 throws.
        /// </summary>
        [TestMethod]
        public void Matrix_EmptyAndInvalid()
        {
            var noRows = new double[0, 3];
            CollectionAssert.AreEqual(new[] { double.NegativeInfinity, double.NegativeInfinity, double.NegativeInfinity }, Tools.LogSumExp(noRows, 0));
            Assert.IsEmpty(Tools.LogSumExp(noRows, 1));
            var noColumns = new double[2, 0];
            CollectionAssert.AreEqual(new[] { double.NegativeInfinity, double.NegativeInfinity }, Tools.LogSumExp(noColumns, 1));
            Assert.IsEmpty(Tools.LogSumExp(noColumns, 0));
            Assert.ThrowsExactly<ArgumentNullException>(() => Tools.LogSumExp((double[,])null!, 0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Tools.LogSumExp(new double[1, 1], 2));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Tools.LogSumExp(new double[1, 1], -1));
        }

        /// <summary>64 normal values with standard deviation 20 (NumPy default_rng(20261108), after a 3 by 4 draw).</summary>
        private static readonly double[] Random64 =
        {
            0.21442459252790846, -12.937574556161044, -9.430370624444988, 21.094252318144385, 8.377650153771869,
            6.048339106072779, 15.345112423954275, -4.339761154491149, 12.81924300866235, 2.571133428130227,
            5.145922594537518, 8.212841764557478, 33.96889867402037, -22.355661926445055, 16.043469414023487,
            1.167600900961216, -3.2204317696412845, -6.362527690528225, -32.41075593274388, 56.9328629021649,
            -5.013014241811421, 42.25885278115213, 25.55787601807968, 5.628678606025154, 17.07445646281326,
            37.18291090048283, 17.635047650593453, 19.641922273963033, 32.05385579124055, 22.56146630908091,
            5.764611754498846, 2.3804765004570956, 15.174936313840634, 15.853645885706417, -15.472461406597269,
            22.32629264479104, 1.0498851243698102, -53.38488202059396, -17.210269799084156, 11.03627345159665,
            15.730789104740593, 6.399908347587733, -26.54421344255624, 6.268638711780149, 5.158294491675797,
            2.655992156391787, -4.765519072277532, 7.97123466425509, -11.348980745484798, 12.200040565153758,
            3.8608220707389256, -13.719742881035852, 11.502197606005515, 9.042158829299096, 2.5119700491044568,
            -17.095359956840824, -17.260161489799827, -11.788443204066992, -0.418143515964807, -0.7850690450111197,
            8.232380647195765, 25.19364404235047, -2.858677795134685, -21.782496575602085,
        };
    }
}
