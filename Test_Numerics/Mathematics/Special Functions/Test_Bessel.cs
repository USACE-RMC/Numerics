using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Mathematics.SpecialFunctions;

namespace Mathematics.SpecialFunctions
{
    /// <summary>
    /// Unit tests for the Bessel functions class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// Reference values computed using Wolfram Alpha and Abramowitz &amp; Stegun tables.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Bessel
    {

        #region Modified Bessel Functions of the First Kind (I)

        /// <summary>
        /// Test I0 at known values.
        /// </summary>
        [TestMethod]
        public void Test_I0()
        {
            // I0(0) = 1
            Assert.AreEqual(1.0, Bessel.I0(0), 1e-7);
            // I0(1) ≈ 1.2660658778
            Assert.AreEqual(1.2660658778, Bessel.I0(1), 1e-5);
            // I0(2) ≈ 2.2795853024
            Assert.AreEqual(2.2795853024, Bessel.I0(2), 1e-5);
            // I0(5) ≈ 27.2398718236
            Assert.AreEqual(27.2398718236, Bessel.I0(5), 1e-2);
            // I0(10) ≈ 2815.7166284
            Assert.AreEqual(2815.7166284, Bessel.I0(10), 1e0);
            // I0 is even: I0(-x) = I0(x)
            Assert.AreEqual(Bessel.I0(3.5), Bessel.I0(-3.5), 1e-10);
        }

        /// <summary>
        /// Test I1 at known values.
        /// </summary>
        [TestMethod]
        public void Test_I1()
        {
            // I1(0) = 0
            Assert.AreEqual(0.0, Bessel.I1(0), 1e-7);
            // I1(1) ≈ 0.5651591040
            Assert.AreEqual(0.5651591040, Bessel.I1(1), 1e-5);
            // I1(2) ≈ 1.5906368546
            Assert.AreEqual(1.5906368546, Bessel.I1(2), 1e-5);
            // I1(5) ≈ 24.3356421088
            Assert.AreEqual(24.3356421088, Bessel.I1(5), 1e-2);
            // I1 is odd: I1(-x) = -I1(x)
            Assert.AreEqual(-Bessel.I1(2.5), Bessel.I1(-2.5), 1e-10);
        }

        /// <summary>
        /// I0e and I1e match scipy.special.i0e/i1e (Cephes) and mpmath over [0, 1E4], including small and
        /// near-breakpoint arguments.
        /// </summary>
        /// <remarks>
        /// Oracle: scipy 1.17.1 scipy.special.i0e(x) and i1e(x) (printed with repr), cross-checked with mpmath 1.4.1
        /// besseli(n, x) * exp(-x) at 50 digits. Cephes states peak relative errors of 5.4E-16 (i0e) and 2.0E-15
        /// (i1e); the largest mpmath difference on this grid is 8.6E-16 at i1e(7.9), which scipy shares.
        /// </remarks>
        [TestMethod]
        public void Test_I0e_I1e_Oracle()
        {
            double[] x = { 1e-300, 1e-10, 1e-3, 0.5, 1.0, 2.5, 3.75, 7.9, 8.0, 8.1, 10.0, 30.0, 100.0, 500.0, 709.0, 800.0, 1000.0, 5000.0, 10000.0 };
            // scipy.special.i0e(x)
            double[] i0e = { 1.0, 0.9999999998999998, 0.9990007495835156, 0.64503527044915, 0.46575960759364043, 0.27004644161220276, 0.2144570512300487, 0.14436986414104191, 0.1434317818568503, 0.14251180948829525, 0.1278333371634286, 0.0731459464822373, 0.03994437929909668, 0.017845706500153165, 0.014985235429889937, 0.014106945005869183, 0.012617240455891257, 0.005642036898744589, 0.0039894726746047314 };
            // scipy.special.i1e(x)
            double[] i1e = { 5.000000000000001e-301, 4.999999999500001e-11, 0.0004995003123542213, 0.15642080318487173, 0.2079104153497085, 0.20658464953126654, 0.1829684209308909, 0.13489649943989365, 0.13414249329269812, 0.13340068832583662, 0.1212626813844555, 0.07191633059864755, 0.03974415302513025, 0.017827851852898053, 0.014974663831277435, 0.014098125406526996, 0.01261093025692863, 0.005641472666838885, 0.003989273195983662 };
            for (int i = 0; i < x.Length; i++)
            {
                Assert.AreEqual(i0e[i], Bessel.I0e(x[i]), 1E-15 * i0e[i], $"I0e({x[i]})");
                Assert.AreEqual(i1e[i], Bessel.I1e(x[i]), 1E-15 * i1e[i], $"I1e({x[i]})");
            }
            // mpmath besseli(1, 7.9) * exp(-7.9) = 0.13489649943989377 (worst case on the grid)
            Assert.AreEqual(0.13489649943989377, Bessel.I1e(7.9), 2E-15 * 0.13489649943989377);
        }

        /// <summary>
        /// I0e and I1e symmetry, limits and special values.
        /// </summary>
        /// <remarks>
        /// Oracle: scipy.special.i0e(0) = 1, i1e(0) = 0, i0e(inf) = 0, i1e(inf) = 0, i1e(-inf) = -0, i0e(nan) = nan,
        /// i0e(-2.5) = 0.27004644161220276, i1e(-2.5) = -0.20658464953126654, i1e(-800) = -0.014098125406526996.
        /// </remarks>
        [TestMethod]
        public void Test_I0e_I1e_EdgeCases()
        {
            Assert.AreEqual(1.0, Bessel.I0e(0.0));
            Assert.AreEqual(0.0, Bessel.I1e(0.0));
            Assert.AreEqual(0.0, Bessel.I0e(double.PositiveInfinity));
            Assert.AreEqual(0.0, Bessel.I0e(double.NegativeInfinity));
            Assert.AreEqual(0.0, Bessel.I1e(double.PositiveInfinity));
            Assert.AreEqual(0.0, Bessel.I1e(double.NegativeInfinity));
            Assert.IsTrue(double.IsNaN(Bessel.I0e(double.NaN)));
            Assert.IsTrue(double.IsNaN(Bessel.I1e(double.NaN)));
            Assert.AreEqual(0.27004644161220276, Bessel.I0e(-2.5), 1E-15 * 0.27004644161220276);
            Assert.AreEqual(-0.20658464953126654, Bessel.I1e(-2.5), 1E-15 * 0.20658464953126654);
            Assert.AreEqual(-0.014098125406526996, Bessel.I1e(-800.0), 1E-15 * 0.014098125406526996);
            Assert.AreEqual(Bessel.I0e(3.3), Bessel.I0e(-3.3), 0d);
            Assert.AreEqual(-Bessel.I1e(12.0), Bessel.I1e(-12.0), 0d);
        }

        /// <summary>
        /// I0 and I1 are the scaled expansions times exp(|x|): rescaled I0e and I1e agree with them to a few ulp,
        /// including across the branch point 8, the split exponential above 709 and negative arguments.
        /// </summary>
        [TestMethod]
        public void Test_I0e_I1e_ConsistentWithUnscaled()
        {
            for (double x = -705.0; x <= 705.0; x += 0.37)
            {
                foreach (double v in new[] { x, 8.0, -8.0, 0.0 })
                {
                    double i0 = Bessel.I0e(v) * Math.Exp(Math.Abs(v));
                    double i1 = Bessel.I1e(v) * Math.Exp(Math.Abs(v));
                    Assert.AreEqual(i0, Bessel.I0(v), 4E-16 * i0, $"I0({v})");
                    Assert.AreEqual(i1, Bessel.I1(v), 4E-16 * Math.Abs(i1), $"I1({v})");
                }
            }
            // Above 709 the exponential is split, so I0 stays finite until I0 itself overflows near 713.9.
            Assert.IsTrue(Tools.IsFinite(Bessel.I0(713.0)));
            Assert.IsTrue(double.IsPositiveInfinity(Bessel.I0(714.0)));
            Assert.AreEqual(-Bessel.I1(713.0), Bessel.I1(-713.0), 0d);
        }

        /// <summary>
        /// Test In for arbitrary integer order.
        /// </summary>
        [TestMethod]
        public void Test_In()
        {
            // In(0, x) = I0(x)
            Assert.AreEqual(Bessel.I0(2), Bessel.In(0, 2), 1e-7);
            // In(1, x) = I1(x)
            Assert.AreEqual(Bessel.I1(2), Bessel.In(1, 2), 1e-7);
            // I2(1) ≈ 0.1357476698
            Assert.AreEqual(0.1357476698, Bessel.In(2, 1), 1e-5);
            // I3(2) ≈ 0.2127399592
            Assert.AreEqual(0.2127399592, Bessel.In(3, 2), 1e-5);
            // I5(3) ≈ 0.0912064773
            Assert.AreEqual(0.0912064773, Bessel.In(5, 3), 1e-5);
            // In(n, 0) = 0 for n > 0
            Assert.AreEqual(0.0, Bessel.In(5, 0), 1e-10);
            // Parity: In(n, -x) = (-1)^n * In(n, x)
            Assert.AreEqual(Bessel.In(2, 3.0), Bessel.In(2, -3.0), 1e-7); // even n
            Assert.AreEqual(-Bessel.In(3, 3.0), Bessel.In(3, -3.0), 1e-7); // odd n
        }

        /// <summary>
        /// Test that In throws for negative order.
        /// </summary>
        [TestMethod]
        public void Test_In_NegativeOrder()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.In(-1, 1.0));
        }

        #endregion

        #region Modified Bessel Functions of the Second Kind (K)

        /// <summary>
        /// Test K0 at known values.
        /// </summary>
        [TestMethod]
        public void Test_K0()
        {
            // K0(1) ≈ 0.4210244382
            Assert.AreEqual(0.4210244382, Bessel.K0(1), 1e-5);
            // K0(2) ≈ 0.1138938727
            Assert.AreEqual(0.1138938727, Bessel.K0(2), 1e-5);
            // K0(5) ≈ 0.003691098334
            Assert.AreEqual(0.003691098334, Bessel.K0(5), 1e-6);
            // K0(0.1) ≈ 2.4270690248
            Assert.AreEqual(2.4270690248, Bessel.K0(0.1), 1e-4);
            // K0 is monotonically decreasing
            Assert.IsGreaterThan(Bessel.K0(2), Bessel.K0(1));
            Assert.IsGreaterThan(Bessel.K0(5), Bessel.K0(2));
        }

        /// <summary>
        /// Test K1 at known values.
        /// </summary>
        [TestMethod]
        public void Test_K1()
        {
            // K1(1) ≈ 0.6019072302
            Assert.AreEqual(0.6019072302, Bessel.K1(1), 1e-5);
            // K1(2) ≈ 0.1398658818
            Assert.AreEqual(0.1398658818, Bessel.K1(2), 1e-5);
            // K1(5) ≈ 0.004044613184
            Assert.AreEqual(0.004044613184, Bessel.K1(5), 1e-6);
            // K1 is monotonically decreasing for x > 0
            Assert.IsGreaterThan(Bessel.K1(2), Bessel.K1(1));
        }

        /// <summary>
        /// Test Kn for arbitrary integer order.
        /// </summary>
        [TestMethod]
        public void Test_Kn()
        {
            // Kn(0, x) = K0(x)
            Assert.AreEqual(Bessel.K0(2), Bessel.Kn(0, 2), 1e-7);
            // Kn(1, x) = K1(x)
            Assert.AreEqual(Bessel.K1(2), Bessel.Kn(1, 2), 1e-7);
            // K2(1) ≈ 1.6248388986
            Assert.AreEqual(1.6248388986, Bessel.Kn(2, 1), 1e-4);
            // K3(2) ≈ 0.6473853909
            Assert.AreEqual(0.6473853909, Bessel.Kn(3, 2), 1e-4);
            // Kn increases with n for fixed x
            Assert.IsGreaterThan(Bessel.Kn(2, 2), Bessel.Kn(3, 2));
        }

        /// <summary>
        /// Test that K functions throw for non-positive argument.
        /// </summary>
        [TestMethod]
        public void Test_K_InvalidArgument()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.K0(0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.K0(-1));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.K1(0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.K1(-1));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Kn(2, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Kn(2, -1));
        }

        #endregion

        #region Bessel Functions of the First Kind (J)

        /// <summary>
        /// Test J0 at known values.
        /// </summary>
        [TestMethod]
        public void Test_J0()
        {
            // J0(0) = 1
            Assert.AreEqual(1.0, Bessel.J0(0), 1e-7);
            // J0(1) ≈ 0.7651976866
            Assert.AreEqual(0.7651976866, Bessel.J0(1), 1e-7);
            // J0(2) ≈ 0.2238907791
            Assert.AreEqual(0.2238907791, Bessel.J0(2), 1e-7);
            // J0(5) ≈ -0.1775967713
            Assert.AreEqual(-0.1775967713, Bessel.J0(5), 1e-7);
            // J0(10) ≈ -0.2459357645
            Assert.AreEqual(-0.2459357645, Bessel.J0(10), 1e-7);
            // J0(20) ≈ 0.1670246643
            Assert.AreEqual(0.1670246643, Bessel.J0(20), 1e-6);
            // J0 is even: J0(-x) = J0(x)
            Assert.AreEqual(Bessel.J0(3.0), Bessel.J0(-3.0), 1e-10);
        }

        /// <summary>
        /// Test J1 at known values.
        /// </summary>
        [TestMethod]
        public void Test_J1()
        {
            // J1(0) = 0
            Assert.AreEqual(0.0, Bessel.J1(0), 1e-7);
            // J1(1) ≈ 0.4400505857
            Assert.AreEqual(0.4400505857, Bessel.J1(1), 1e-7);
            // J1(2) ≈ 0.5767248078
            Assert.AreEqual(0.5767248078, Bessel.J1(2), 1e-7);
            // J1(5) ≈ -0.3275791376
            Assert.AreEqual(-0.3275791376, Bessel.J1(5), 1e-7);
            // J1(10) ≈ 0.0434727462
            Assert.AreEqual(0.0434727462, Bessel.J1(10), 1e-7);
            // J1 is odd: J1(-x) = -J1(x)
            Assert.AreEqual(-Bessel.J1(3.0), Bessel.J1(-3.0), 1e-10);
        }

        /// <summary>
        /// Test Jn for arbitrary integer order.
        /// </summary>
        [TestMethod]
        public void Test_Jn()
        {
            // Jn(0, x) = J0(x)
            Assert.AreEqual(Bessel.J0(2), Bessel.Jn(0, 2), 1e-7);
            // Jn(1, x) = J1(x)
            Assert.AreEqual(Bessel.J1(2), Bessel.Jn(1, 2), 1e-7);
            // J2(1) ≈ 0.1149034849
            Assert.AreEqual(0.1149034849, Bessel.Jn(2, 1), 1e-7);
            // J3(2) ≈ 0.1289432494
            Assert.AreEqual(0.1289432494, Bessel.Jn(3, 2), 1e-7);
            // J5(5) ≈ 0.2611405461
            Assert.AreEqual(0.2611405461, Bessel.Jn(5, 5), 1e-6);
            // Jn(n, 0) = 0 for n > 0
            Assert.AreEqual(0.0, Bessel.Jn(5, 0), 1e-10);
            // J10(10) ≈ 0.2074861066
            Assert.AreEqual(0.2074861066, Bessel.Jn(10, 10), 1e-5);
            // Parity: Jn(n, -x) = (-1)^n * Jn(n, x)
            Assert.AreEqual(Bessel.Jn(2, 3.0), Bessel.Jn(2, -3.0), 1e-7);   // even n
            Assert.AreEqual(-Bessel.Jn(3, 3.0), Bessel.Jn(3, -3.0), 1e-7);  // odd n
        }

        /// <summary>
        /// Test Jn using the forward recurrence path (x > n).
        /// </summary>
        [TestMethod]
        public void Test_Jn_ForwardRecurrence()
        {
            // J2(10) ≈ 0.2546303137: uses forward recurrence since |x| > n
            Assert.AreEqual(0.2546303137, Bessel.Jn(2, 10), 1e-6);
            // J3(10) ≈ 0.0583793794
            Assert.AreEqual(0.0583793794, Bessel.Jn(3, 10), 1e-6);
        }

        /// <summary>
        /// Test Jn using the Miller's downward recurrence path (x &lt;= n).
        /// </summary>
        [TestMethod]
        public void Test_Jn_MillerRecurrence()
        {
            // J10(5) ≈ 0.0014678027: uses Miller's since |x| <= n
            Assert.AreEqual(0.0014678027, Bessel.Jn(10, 5), 1e-7);
            // J20(10) ≈ 1.15134e-5 (very small, deep into Miller's regime)
            Assert.AreEqual(1.1513369e-5, Bessel.Jn(20, 10), 1e-9);
        }

        /// <summary>
        /// Test that Jn throws for negative order.
        /// </summary>
        [TestMethod]
        public void Test_Jn_NegativeOrder()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Jn(-1, 1.0));
        }

        #endregion

        #region Bessel Functions of the Second Kind (Y)

        /// <summary>
        /// Test Y0 at known values.
        /// </summary>
        [TestMethod]
        public void Test_Y0()
        {
            // Y0(1) ≈ 0.0882569642
            Assert.AreEqual(0.0882569642, Bessel.Y0(1), 1e-7);
            // Y0(2) ≈ 0.5103756726
            Assert.AreEqual(0.5103756726, Bessel.Y0(2), 1e-7);
            // Y0(5) ≈ -0.3085176252
            Assert.AreEqual(-0.3085176252, Bessel.Y0(5), 1e-7);
            // Y0(10) ≈ 0.0556711673
            Assert.AreEqual(0.0556711673, Bessel.Y0(10), 1e-7);
            // Y0(20) ≈ 0.0626405967
            Assert.AreEqual(0.0626405967, Bessel.Y0(20), 1e-6);
        }

        /// <summary>
        /// Test Y1 at known values.
        /// </summary>
        [TestMethod]
        public void Test_Y1()
        {
            // Y1(1) ≈ -0.7812128213
            Assert.AreEqual(-0.7812128213, Bessel.Y1(1), 1e-7);
            // Y1(2) ≈ -0.1070324315
            Assert.AreEqual(-0.1070324315, Bessel.Y1(2), 1e-7);
            // Y1(5) ≈ 0.1478631434
            Assert.AreEqual(0.1478631434, Bessel.Y1(5), 1e-7);
            // Y1(10) ≈ 0.2490154242
            Assert.AreEqual(0.2490154242, Bessel.Y1(10), 1e-7);
        }

        /// <summary>
        /// Test Yn for arbitrary integer order.
        /// </summary>
        [TestMethod]
        public void Test_Yn()
        {
            // Yn(0, x) = Y0(x)
            Assert.AreEqual(Bessel.Y0(2), Bessel.Yn(0, 2), 1e-7);
            // Yn(1, x) = Y1(x)
            Assert.AreEqual(Bessel.Y1(2), Bessel.Yn(1, 2), 1e-7);
            // Y2(1) ≈ -1.6506826068
            Assert.AreEqual(-1.6506826068, Bessel.Yn(2, 1), 1e-5);
            // Y3(2) ≈ -1.1277837769
            Assert.AreEqual(-1.1277837769, Bessel.Yn(3, 2), 1e-5);
            // Y5(5) ≈ -0.4536948225
            Assert.AreEqual(-0.4536948225, Bessel.Yn(5, 5), 1e-5);
        }

        /// <summary>
        /// Test that Y functions throw for non-positive argument.
        /// </summary>
        [TestMethod]
        public void Test_Y_InvalidArgument()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Y0(0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Y0(-1));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Y1(0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Y1(-1));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Yn(2, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Bessel.Yn(2, -1));
        }

        #endregion

        #region Cross-Function Identities

        /// <summary>
        /// Verify the Wronskian identity: J0(x)Y1(x) - J1(x)Y0(x) = -2/(pi*x).
        /// </summary>
        [TestMethod]
        public void Test_Wronskian_J0Y1_J1Y0()
        {
            double[] testValues = { 0.5, 1.0, 2.0, 5.0, 10.0 };
            foreach (double x in testValues)
            {
                double wronskian = Bessel.J0(x) * Bessel.Y1(x) - Bessel.J1(x) * Bessel.Y0(x);
                double expected = -2.0 / (Math.PI * x);
                Assert.AreEqual(expected, wronskian, 1e-6, $"Wronskian failed at x={x}");
            }
        }

        /// <summary>
        /// Verify the modified Wronskian identity: I0(x)K1(x) + I1(x)K0(x) = 1/x.
        /// </summary>
        [TestMethod]
        public void Test_Wronskian_I0K1_I1K0()
        {
            double[] testValues = { 0.5, 1.0, 2.0, 5.0, 10.0 };
            foreach (double x in testValues)
            {
                double wronskian = Bessel.I0(x) * Bessel.K1(x) + Bessel.I1(x) * Bessel.K0(x);
                double expected = 1.0 / x;
                Assert.AreEqual(expected, wronskian, 1e-5, $"Modified Wronskian failed at x={x}");
            }
        }

        /// <summary>
        /// Verify that I0'(x) = I1(x) numerically using finite differences.
        /// </summary>
        [TestMethod]
        public void Test_I0_Derivative_Equals_I1()
        {
            double[] testValues = { 0.5, 1.0, 2.0, 5.0 };
            double h = 1e-6;
            foreach (double x in testValues)
            {
                double numericalDerivative = (Bessel.I0(x + h) - Bessel.I0(x - h)) / (2.0 * h);
                Assert.AreEqual(Bessel.I1(x), numericalDerivative, 1e-4, $"I0'(x) != I1(x) at x={x}");
            }
        }

        /// <summary>
        /// Verify the addition formula: J0(x) + 2*J2(x) + 2*J4(x) + 2*J6(x) + ... = 1.
        /// </summary>
        [TestMethod]
        public void Test_J_SumIdentity()
        {
            double x = 3.0;
            double sum = Bessel.J0(x);
            for (int k = 1; k <= 20; k++)
            {
                sum += 2.0 * Bessel.Jn(2 * k, x);
            }
            Assert.AreEqual(1.0, sum, 1e-6);
        }

        #endregion

        #region Large Argument Tests

        /// <summary>
        /// Test J0 for large arguments where the asymptotic expansion is used.
        /// </summary>
        [TestMethod]
        public void Test_J0_LargeArgument()
        {
            // J0(50) ≈ 0.0558123276
            Assert.AreEqual(0.0558123276, Bessel.J0(50), 1e-5);
            // J0(100) ≈ 0.0199858503
            Assert.AreEqual(0.0199858503, Bessel.J0(100), 1e-4);
        }

        /// <summary>
        /// Test I0 for large arguments where the asymptotic expansion is used.
        /// </summary>
        [TestMethod]
        public void Test_I0_LargeArgument()
        {
            // I0(20) ≈ 4.355828256e7
            Assert.AreEqual(4.355828256e7, Bessel.I0(20), 1e3);
        }

        /// <summary>
        /// I0 and I1 stay finite up to the edge of the double range instead of overflowing at 709.
        /// </summary>
        /// <remarks>
        /// Oracle: mpmath besseli at 40 digits, rounded to double. The tolerance is 1E-13 relative; the
        /// results equal the rounded references.
        /// </remarks>
        [TestMethod]
        public void Test_I0_I1_NearOverflow()
        {
            Assert.AreEqual(9.087162727263793E+306, Bessel.I0(711), 1E-13 * 9.087162727263793E+306);
            Assert.AreEqual(9.080770067322847E+306, Bessel.I1(711), 1E-13 * 9.080770067322847E+306);
            Assert.AreEqual(-9.080770067322847E+306, Bessel.I1(-711), 1E-13 * 9.080770067322847E+306);
            Assert.AreEqual(2.029763060075537E+306, Bessel.I0(709.5), 1E-13 * 2.029763060075537E+306);
            Assert.AreEqual(double.PositiveInfinity, Bessel.I0(720));
            Assert.AreEqual(double.PositiveInfinity, Bessel.I0(double.PositiveInfinity));
        }

        /// <summary>
        /// The Miller start index for In grows with x when x &gt; n, so large arguments converge.
        /// </summary>
        /// <remarks>
        /// Oracle: mpmath besseli and scipy 1.17 iv. The tolerance is 1E-13 relative; the measured
        /// deviation is at most 6E-16. The small-argument cases use the base start index.
        /// </remarks>
        [TestMethod]
        public void Test_In_LargeArgument()
        {
            Assert.AreEqual(1.0523843193243106E+42, Bessel.In(2, 100), 1E-13 * 1.0523843193243106E+42);
            Assert.AreEqual(2.442904816107921E+215, Bessel.In(5, 500), 1E-13 * 2.442904816107921E+215);
            Assert.AreEqual(39312785.221040756, Bessel.In(2, 20), 1E-13 * 39312785.221040756);
            Assert.AreEqual(0.33783461833568074, Bessel.In(2, 1.5), 1E-13 * 0.33783461833568074);
            Assert.AreEqual(0.002645111968990286, Bessel.In(3, 0.5), 1E-13 * 0.002645111968990286);
            Assert.AreEqual(0.004580044419176052, Bessel.In(10, 5), 1E-13 * 0.004580044419176052);
        }

        /// <summary>
        /// Test K0 for moderate arguments.
        /// </summary>
        [TestMethod]
        public void Test_K0_ModerateArgument()
        {
            // K0(10) ≈ 1.778006232e-5
            Assert.AreEqual(1.778006232e-5, Bessel.K0(10), 1e-8);
        }

        #endregion

    }
}
