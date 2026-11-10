using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.SpecialFunctions;

namespace Mathematics.SpecialFunctions
{
    /// <summary>
    /// Tests the IList&lt;double&gt; overloads of the Gamma, Erf, Bessel, Debye and Beta special functions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// Every list value must be bitwise the scalar member at that value, and a value at which the scalar member
    /// throws must make the list call throw the same exception type. The scalar members have their own oracle
    /// tests, so the list forms need no separate oracle.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_SpecialFunctionListOverloads
    {
        /// <summary>Arguments covering signs, zero, the branch points, large values and the special values.</summary>
        private static readonly double[] Values =
        {
            double.NegativeInfinity, -750.0, -30.0, -8.0, -2.5, -1.0, -1E-12, -0.0, 0.0, 1E-300, 1E-12, 0.1, 0.5, 0.75,
            1.0, 2.0, 2.404825557695773, 3.0, 7.9, 8.0, 8.1, 25.0, 26.0, 100.0, 700.0, 750.0, 1E6, double.PositiveInfinity, double.NaN,
        };

        /// <summary>Second arguments for the two-argument functions.</summary>
        private static readonly double[] Second =
        {
            0.5, 3.0, 0.25, 1.0, 10.0, 0.9, 2.5, 30.0, 0.01, 1E-5, 50.0, 1.5, 0.75, 4.0, 0.3, 1E3, 0.5, 2.0, 0.2, 25.0,
            8.0, 0.6, 1.0, 0.05, 2.0, 0.99, 7.0, 0.5, 1.0,
        };

        /// <summary>Third arguments for the three-argument functions.</summary>
        private static readonly double[] Third =
        {
            0.3, 0.5, 0.9, 0.1, 0.999, 1E-6, 0.25, 0.75, 0.5, 0.01, 0.6, 0.4, 0.2, 0.95, 0.05, 0.7, 0.33, 0.8, 0.45, 0.5,
            0.15, 0.85, 0.65, 0.35, 0.55, 0.02, 0.98, 0.5, 0.5,
        };

        /// <summary>Asserts that a one-argument list form is bitwise the scalar member, value by value.</summary>
        private static void Same(string name, Func<double, double> scalar, Func<IList<double>, double[]> list, IList<double> values)
        {
            var good = new List<double>();
            foreach (double v in values)
            {
                Type thrown = Thrown(() => scalar(v));
                if (thrown == null) { good.Add(v); continue; }
                Assert.AreEqual(thrown, Thrown(() => list(new[] { v })), $"{name}({v}) exception");
            }
            var result = list(good);
            Assert.HasCount(good.Count, result);
            for (int i = 0; i < good.Count; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(scalar(good[i])), BitConverter.DoubleToInt64Bits(result[i]), $"{name}({good[i]})");
        }

        /// <summary>Asserts that a list form over parallel lists is bitwise the scalar member, row by row.</summary>
        private static void SameRows(string name, Func<double[], double> scalar, Func<IList<double>[], double[]> list, params IList<double>[] columns)
        {
            int count = columns[0].Count;
            var good = new List<int>();
            for (int i = 0; i < count; i++)
            {
                double[] row = columns.Select(c => c[i]).ToArray();
                Type thrown = Thrown(() => scalar(row));
                if (thrown == null) { good.Add(i); continue; }
                Assert.AreEqual(thrown, Thrown(() => list(columns.Select(c => (IList<double>)new[] { c[i] }).ToArray())), $"{name} row {i} exception");
            }
            var result = list(columns.Select(c => (IList<double>)good.Select(i => c[i]).ToArray()).ToArray());
            Assert.HasCount(good.Count, result);
            for (int k = 0; k < good.Count; k++)
            {
                double[] row = columns.Select(c => c[good[k]]).ToArray();
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(scalar(row)), BitConverter.DoubleToInt64Bits(result[k]), $"{name} row {good[k]}");
            }
        }

        /// <summary>Returns the type of the exception the action throws, or null.</summary>
        private static Type Thrown(Action action)
        {
            try { action(); return null; }
            catch (Exception e) { return e.GetType(); }
        }

        /// <summary>The non-negative values, for the functions defined on [0, ∞) or (0, ∞).</summary>
        private static double[] Positive => Values.Select(Math.Abs).ToArray();

        /// <summary>Probabilities in [0, 1] and just outside, for the inverse functions.</summary>
        private static readonly double[] Probabilities = { -0.1, 0.0, 1E-300, 1E-10, 0.01, 0.3, 0.5, 0.7, 0.99, 1.0 - 1E-10, 1.0, 1.1, double.NaN };

        /// <summary>The Gamma list forms are bitwise the scalar members.</summary>
        [TestMethod]
        public void Gamma_ListFormsMatchScalar()
        {
            Same("Stirling", Gamma.Stirling, Gamma.Stirling, Positive);
            Same("Lanczos", Gamma.Lanczos, Gamma.Lanczos, Values);
            Same("Function", Gamma.Function, Gamma.Function, Values);
            Same("Digamma", Gamma.Digamma, Gamma.Digamma, Values);
            Same("Trigamma", Gamma.Trigamma, Gamma.Trigamma, Values);
            Same("LogGamma", Gamma.LogGamma, Gamma.LogGamma, Values);
            SameRows("Incomplete", r => Gamma.Incomplete(r[0], r[1]), c => Gamma.Incomplete(c[0], c[1]), Values, Second);
            SameRows("LowerIncomplete", r => Gamma.LowerIncomplete(r[0], r[1]), c => Gamma.LowerIncomplete(c[0], c[1]), Second, Values);
            SameRows("UpperIncomplete", r => Gamma.UpperIncomplete(r[0], r[1]), c => Gamma.UpperIncomplete(c[0], c[1]), Second, Values);
            var shapes = Probabilities.Select((p, i) => Second[i]).ToArray();
            SameRows("InverseLowerIncomplete", r => Gamma.InverseLowerIncomplete(r[0], r[1]), c => Gamma.InverseLowerIncomplete(c[0], c[1]), shapes, Probabilities);
            SameRows("InverseUpperIncomplete", r => Gamma.InverseUpperIncomplete(r[0], r[1]), c => Gamma.InverseUpperIncomplete(c[0], c[1]), shapes, Probabilities);
        }

        /// <summary>The Erf list forms are bitwise the scalar members.</summary>
        [TestMethod]
        public void Erf_ListFormsMatchScalar()
        {
            Same("Function", Erf.Function, Erf.Function, Values);
            Same("Erfc", Erf.Erfc, Erf.Erfc, Values);
            var signed = Probabilities.Concat(Probabilities.Select(p => -p)).ToArray();
            Same("InverseErf", Erf.InverseErf, Erf.InverseErf, signed);
            Same("InverseErfc", Erf.InverseErfc, Erf.InverseErfc, Probabilities.Concat(new[] { 1.5, 2.0, 2.5 }).ToArray());
        }

        /// <summary>The Bessel list forms are bitwise the scalar members, including the order-n forms.</summary>
        [TestMethod]
        public void Bessel_ListFormsMatchScalar()
        {
            Same("I0", Bessel.I0, Bessel.I0, Values);
            Same("I1", Bessel.I1, Bessel.I1, Values);
            Same("I0e", Bessel.I0e, Bessel.I0e, Values);
            Same("I1e", Bessel.I1e, Bessel.I1e, Values);
            Same("K0", Bessel.K0, Bessel.K0, Values);
            Same("K1", Bessel.K1, Bessel.K1, Values);
            Same("J0", Bessel.J0, Bessel.J0, Values);
            Same("J1", Bessel.J1, Bessel.J1, Values);
            Same("Y0", Bessel.Y0, Bessel.Y0, Values);
            Same("Y1", Bessel.Y1, Bessel.Y1, Values);
            var finite = Values.Where(v => Math.Abs(v) < 1E3).ToArray();
            foreach (int n in new[] { 0, 1, 2, 5, 20, -1 })
            {
                Same($"In({n})", x => Bessel.In(n, x), x => Bessel.In(n, x), finite);
                Same($"Kn({n})", x => Bessel.Kn(n, x), x => Bessel.Kn(n, x), finite);
                Same($"Jn({n})", x => Bessel.Jn(n, x), x => Bessel.Jn(n, x), finite);
                Same($"Yn({n})", x => Bessel.Yn(n, x), x => Bessel.Yn(n, x), finite);
            }
        }

        /// <summary>The Debye list forms are bitwise the scalar members.</summary>
        [TestMethod]
        public void Debye_ListFormsMatchScalar()
        {
            Same("Function", Debye.Function, Debye.Function, Values);
            Same("FunctionOrderOne", Debye.FunctionOrderOne, Debye.FunctionOrderOne, Values);
        }

        /// <summary>The Beta list forms are bitwise the scalar members.</summary>
        [TestMethod]
        public void Beta_ListFormsMatchScalar()
        {
            var a = Second;
            var b = Enumerable.Reverse(Second).ToArray();
            SameRows("Function", r => Beta.Function(r[0], r[1]), c => Beta.Function(c[0], c[1]), a, b);
            SameRows("LogFunction", r => Beta.LogFunction(r[0], r[1]), c => Beta.LogFunction(c[0], c[1]), a, b);
            SameRows("Incomplete", r => Beta.Incomplete(r[0], r[1], r[2]), c => Beta.Incomplete(c[0], c[1], c[2]), a, b, Third);
            SameRows("Incbcf", r => Beta.Incbcf(r[0], r[1], r[2]), c => Beta.Incbcf(c[0], c[1], c[2]), a, b, Third);
            SameRows("Incbd", r => Beta.Incbd(r[0], r[1], r[2]), c => Beta.Incbd(c[0], c[1], c[2]), a, b, Third);
            SameRows("IncompleteInverse", r => Beta.IncompleteInverse(r[0], r[1], r[2]), c => Beta.IncompleteInverse(c[0], c[1], c[2]), a, b, Third);
            SameRows("PowerSeries", r => Beta.PowerSeries(r[0], r[1], r[2]), c => Beta.PowerSeries(c[0], c[1], c[2]), a, b, Third);
            var logBeta = a.Select((v, i) => Math.Log(Beta.Function(v, b[i]))).ToArray();
            SameRows("IncompleteRatio", r => Beta.IncompleteRatio(r[0], r[1], r[2], r[3]), c => Beta.IncompleteRatio(c[0], c[1], c[2], c[3]), Third, a, b, logBeta);
        }

        /// <summary>Null lists and lists of different lengths are refused; empty lists give empty arrays.</summary>
        [TestMethod]
        public void NullMismatchedAndEmptyLists()
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => Bessel.J0((IList<double>)null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => Bessel.Jn(3, (IList<double>)null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => Gamma.LowerIncomplete(new[] { 1.0 }, null!));
            Assert.ThrowsExactly<ArgumentNullException>(() => Beta.Incomplete(new[] { 1.0 }, null!, new[] { 0.5 }));
            Assert.ThrowsExactly<ArgumentException>(() => Gamma.LowerIncomplete(new[] { 1.0, 2.0 }, new[] { 1.0 }));
            Assert.ThrowsExactly<ArgumentException>(() => Beta.IncompleteRatio(new[] { 0.5 }, new[] { 1.0 }, new[] { 1.0, 2.0 }, new[] { 0.0 }));
            Assert.IsEmpty(Debye.Function(new double[0]));
            Assert.IsEmpty(Gamma.UpperIncomplete(new double[0], new double[0]));
            Assert.IsEmpty(Bessel.Kn(4, new List<double>()));
        }
    }
}
