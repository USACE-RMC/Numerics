using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;

namespace Utilities
{
    /// <summary>
    /// Accuracy tests for <see cref="Tools.Log1p(double)"/>, which applies the compensated evaluation at every finite
    /// argument above −1.
    /// </summary>
    /// <remarks>
    /// The reference values are Python 3.14 <c>math.log1p</c> (the C library's <c>log1p</c>). On Windows (UCRT) the
    /// compensated form agrees with it bit for bit over 301,000 random arguments spanning (−1, 1e300); on other
    /// platforms the C library rounds <c>log</c> and <c>log1p</c> differently, and the values agree within one unit in
    /// the last place, which is the tolerance the assertions allow. The result is within one unit in the last place
    /// of the correctly rounded value, not correctly rounded everywhere.
    /// </remarks>
    [TestClass]
    public class Test_Log1pAccuracy
    {
        /// <summary>Arguments, including the range 1e-4 to 1e-2 where the former form lost accuracy.</summary>
        private static readonly double[] X =
        {
            -0.9999999, -0.5, -0.01, -0.0050001, -0.001, -0.00012, -0.00010001, 0.00010001, 0.00012, 0.000731, 0.001,
            0.0025, 0.0050001, 0.0099, 0.01, 0.05, 0.3, 1.0, 7.5, 1000.0, 1E15, 1E300,
        };

        /// <summary>Python <c>math.log1p</c> at <see cref="X"/>.</summary>
        private static readonly double[] Expected =
        {
            -16.118095651484676, -0.6931471805599453, -0.010050335853501442, -0.005012642326061896, -0.0010005003335835335,
            -0.00012000720057605185, -0.00010001500133350836, 0.00010000499933335833, 0.00011999280057594817,
            0.00073073294963462, 0.0009995003330835331, 0.002496880198587199, 0.004987641013521686, 0.009851316050374183,
            0.009950330853168083, 0.048790164169432, 0.262364264467491, 0.6931471805599453, 2.1400661634962708,
            6.90875477931522, 34.538776394910684, 690.7755278982137,
        };

        /// <summary>Returns the distance in units in the last place between two finite doubles of one sign.</summary>
        private static long Ulps(double a, double b) => Math.Abs(BitConverter.DoubleToInt64Bits(a) - BitConverter.DoubleToInt64Bits(b));

        /// <summary>The values agree with the C library's log1p to one unit in the last place.</summary>
        [TestMethod]
        public void MatchesTheCLibraryLog1p()
        {
            for (int i = 0; i < X.Length; i++)
                Assert.IsLessThanOrEqualTo(1L, Ulps(Expected[i], Tools.Log1p(X[i])), $"Log1p({X[i]})");
        }

        /// <summary>For |x| ≤ 1e-4 the value is bitwise the former compensated branch, which is unchanged.</summary>
        [TestMethod]
        public void SmallArgumentsAreUnchanged()
        {
            var random = new Random(43);
            for (int i = 0; i < 20000; i++)
            {
                double x = (2 * random.NextDouble() - 1) * Math.Pow(10, -4 - 12 * random.NextDouble());
                double y = 1.0 + x;
                double former = Math.Log(y) - ((y - 1.0) - x) / y;
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(former), BitConverter.DoubleToInt64Bits(Tools.Log1p(x)), $"Log1p({x})");
            }
            foreach (double x in new[] { 1E-4, -1E-4, 0.0, -0.0, 1E-300, -1E-300 })
            {
                double y = 1.0 + x;
                double former = Math.Log(y) - ((y - 1.0) - x) / y;
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(former), BitConverter.DoubleToInt64Bits(Tools.Log1p(x)), $"Log1p({x})");
            }
        }

        /// <summary>The special values are those of log(1 + x).</summary>
        [TestMethod]
        public void SpecialValues()
        {
            Assert.IsTrue(double.IsNegativeInfinity(Tools.Log1p(-1.0)));
            Assert.IsTrue(double.IsNaN(Tools.Log1p(-1.5)));
            Assert.IsTrue(double.IsNaN(Tools.Log1p(double.NegativeInfinity)));
            Assert.IsTrue(double.IsNaN(Tools.Log1p(double.NaN)));
            Assert.IsTrue(double.IsPositiveInfinity(Tools.Log1p(double.PositiveInfinity)));
            Assert.AreEqual(Math.Log(double.MaxValue), Tools.Log1p(double.MaxValue));
        }
    }
}
