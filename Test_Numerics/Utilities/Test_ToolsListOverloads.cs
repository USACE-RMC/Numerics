using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;

namespace Utilities
{
    /// <summary>
    /// Tests for <see cref="Tools.Expm1(IList{double})"/> and <see cref="Tools.Log1p(IList{double})"/>.
    /// </summary>
    /// <remarks>
    /// Each list value must be bitwise the scalar member. NumPy 2.5.4 <c>numpy.expm1</c> and <c>numpy.log1p</c> (the C
    /// library functions) are independent oracles; Numerics uses compensated forms built on <c>exp</c> and <c>log</c>,
    /// a different algorithm, so the oracle comparison allows 1E-15 relative error.
    /// </remarks>
    [TestClass]
    public class Test_ToolsListOverloads
    {
        /// <summary>Arguments around zero, at the 1E-4 switch of the compensated forms, far out and at the special values.</summary>
        private static readonly double[] Values =
        {
            double.NegativeInfinity, -800.0, -40.0, -1.0, -0.5, -1.0001E-4, -1E-4, -9.999E-5, -1E-12, -1E-300, -0.0, 0.0,
            1E-300, 1E-12, 9.999E-5, 1E-4, 1.0001E-4, 0.7, 20.0, 709.0, 710.0, 1E300, double.PositiveInfinity, double.NaN,
        };

        /// <summary>The list forms are bitwise the scalar members.</summary>
        [TestMethod]
        public void ListFormsMatchScalarMembers()
        {
            var expm1 = Tools.Expm1(Values);
            var log1p = Tools.Log1p(new List<double>(Values));
            Assert.HasCount(Values.Length, expm1);
            for (int i = 0; i < Values.Length; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.Expm1(Values[i])), BitConverter.DoubleToInt64Bits(expm1[i]), $"Expm1({Values[i]})");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.Log1p(Values[i])), BitConverter.DoubleToInt64Bits(log1p[i]), $"Log1p({Values[i]})");
            }
        }

        /// <summary>NumPy oracle values.</summary>
        [TestMethod]
        public void ListFormsMatchNumPy()
        {
            var x = new[] { -0.5, -1E-5, -1E-12, 0.0, 1E-9, 3E-5, 0.7, 20.0 };
            var expm1 = new[] { -0.39346934028736663, -9.999950000166668E-06, -9.999999999995E-13, 0.0, 1.0000000005000001E-09, 3.000045000450003E-05, 1.0137527074704766, 485165194.4097903 };
            var log1p = new[] { -0.6931471805599453, -1.0000050000333335E-05, -1.0000000000005E-12, 0.0, 9.999999995E-10, 2.99995500089998E-05, 0.5306282510621704, 3.044522437723423 };
            var e = Tools.Expm1(x);
            var l = Tools.Log1p(x);
            for (int i = 0; i < x.Length; i++)
            {
                Assert.AreEqual(expm1[i], e[i], 1E-15 * Math.Abs(expm1[i]), $"Expm1({x[i]})");
                Assert.AreEqual(log1p[i], l[i], 1E-15 * Math.Abs(log1p[i]), $"Log1p({x[i]})");
            }
        }

        /// <summary>An empty list gives an empty array; a null list is refused.</summary>
        [TestMethod]
        public void EmptyAndNullLists()
        {
            Assert.IsEmpty(Tools.Expm1(new double[0]));
            Assert.IsEmpty(Tools.Log1p(new List<double>()));
            Assert.ThrowsExactly<ArgumentNullException>(() => Tools.Expm1((IList<double>)null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => Tools.Log1p((IList<double>)null!));
        }
    }
}
