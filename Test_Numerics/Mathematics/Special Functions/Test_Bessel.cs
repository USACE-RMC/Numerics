/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
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
            Assert.IsTrue(Bessel.K0(1) > Bessel.K0(2));
            Assert.IsTrue(Bessel.K0(2) > Bessel.K0(5));
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
            Assert.IsTrue(Bessel.K1(1) > Bessel.K1(2));
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
            Assert.IsTrue(Bessel.Kn(3, 2) > Bessel.Kn(2, 2));
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
            // J2(10) ≈ 0.2546303137 — uses forward recurrence since |x| > n
            Assert.AreEqual(0.2546303137, Bessel.Jn(2, 10), 1e-6);
            // J3(10) ≈ 0.0583793794
            Assert.AreEqual(0.0583793794, Bessel.Jn(3, 10), 1e-6);
        }

        /// <summary>
        /// Test Jn using the Miller's downward recurrence path (x <= n).
        /// </summary>
        [TestMethod]
        public void Test_Jn_MillerRecurrence()
        {
            // J10(5) ≈ 0.0014678027 — uses Miller's since |x| <= n
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
