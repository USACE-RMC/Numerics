using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Tests for <see cref="PlottingPositions.Alpha(PlottingPositions.PlottingPostionType)"/>, the single table of the
    /// alpha coefficients, and for the typed weighted overload
    /// <see cref="PlottingPositions.Function(IList{double}, PlottingPositions.PlottingPostionType)"/>.
    /// </summary>
    /// <remarks>
    /// The named formulas read their coefficient from the table; the regression cases pin their positions bitwise
    /// to the formula evaluated at the literal coefficients. The alpha values are those of the published
    /// formulas (Weibull 0, Median 0.3175, Blom 0.375, Cunnane 0.40, Gringorten 0.44, Hazen 0.50).
    /// </remarks>
    [TestClass]
    public class Test_PlottingPositionAlpha
    {
        /// <summary>The alpha formulas and the literal coefficient each named method used before.</summary>
        private static readonly (PlottingPositions.PlottingPostionType Kind, double Alpha, Func<int, double[]> Named)[] Formulas =
        {
            (PlottingPositions.PlottingPostionType.Weibull, 0d, PlottingPositions.Weibull),
            (PlottingPositions.PlottingPostionType.Median, 0.3175d, PlottingPositions.Median),
            (PlottingPositions.PlottingPostionType.Blom, 0.375d, PlottingPositions.Blom),
            (PlottingPositions.PlottingPostionType.Cunnane, 0.4d, PlottingPositions.Cunnane),
            (PlottingPositions.PlottingPostionType.Gringorten, 0.44d, PlottingPositions.Gringorten),
            (PlottingPositions.PlottingPostionType.Hazen, 0.5d, PlottingPositions.Hazen),
        };

        /// <summary>The table holds the published coefficients; Kaplan-Meier and undefined values are refused.</summary>
        [TestMethod]
        public void AlphaTableAndRefusals()
        {
            foreach (var f in Formulas)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(f.Alpha), BitConverter.DoubleToInt64Bits(PlottingPositions.Alpha(f.Kind)), f.Kind.ToString());
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => PlottingPositions.Alpha(PlottingPositions.PlottingPostionType.KaplanMeier));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => PlottingPositions.Alpha((PlottingPositions.PlottingPostionType)42));
        }

        /// <summary>The named methods and the typed unweighted overload keep their positions bit for bit; an undefined type still returns null.</summary>
        [TestMethod]
        public void NamedFormulasAreUnchanged()
        {
            foreach (int n in new[] { 3, 10, 101 })
            {
                foreach (var f in Formulas)
                {
                    var literal = PlottingPositions.Function(n, f.Alpha);
                    var named = f.Named(n);
                    var typed = PlottingPositions.Function(n, f.Kind)!;
                    for (int i = 0; i < n; i++)
                    {
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(literal[i]), BitConverter.DoubleToInt64Bits(named[i]), $"{f.Kind} n={n}");
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(literal[i]), BitConverter.DoubleToInt64Bits(typed[i]), $"{f.Kind} n={n}");
                    }
                }
            }
            Assert.IsNull(PlottingPositions.Function(10, (PlottingPositions.PlottingPostionType)42));
        }

        /// <summary>Unit weights give bitwise the unweighted positions of every formula, Kaplan-Meier included.</summary>
        [TestMethod]
        public void UnitWeightsGiveTheUnweightedPositions()
        {
            foreach (int n in new[] { 3, 7, 50, 1000 })
            {
                var ones = new double[n];
                for (int i = 0; i < n; i++) ones[i] = 1d;
                foreach (PlottingPositions.PlottingPostionType kind in Enum.GetValues(typeof(PlottingPositions.PlottingPostionType)))
                {
                    var weighted = PlottingPositions.Function(ones, kind);
                    var plain = PlottingPositions.Function(n, kind)!;
                    for (int i = 0; i < n; i++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(plain[i]), BitConverter.DoubleToInt64Bits(weighted[i]), $"{kind} n={n} i={i}");
                }
            }
        }

        /// <summary>The alpha formulas equal the alpha overload bitwise; Kaplan-Meier is the cumulative weight share (hand values).</summary>
        [TestMethod]
        public void WeightedPositions()
        {
            var random = new Random(3);
            var weights = new double[60];
            for (int i = 0; i < weights.Length; i++) weights[i] = random.NextDouble() * 3;
            weights[4] = 0d;
            foreach (var f in Formulas)
            {
                var typed = PlottingPositions.Function(weights, f.Kind);
                var direct = PlottingPositions.Function(weights, f.Alpha);
                for (int i = 0; i < weights.Length; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(direct[i]), BitConverter.DoubleToInt64Bits(typed[i]), f.Kind.ToString());
            }
            var km = PlottingPositions.Function(new[] { 1d, 3d, 0d, 4d }, PlottingPositions.PlottingPostionType.KaplanMeier);
            CollectionAssert.AreEqual(new[] { 0.125, 0.5, 0.5, 1.0 }, km);
            Assert.AreEqual(1d, PlottingPositions.Function(new[] { 2.5 }, PlottingPositions.PlottingPostionType.KaplanMeier)[0]);
        }

        /// <summary>Refusals of the typed weighted overload.</summary>
        [TestMethod]
        public void WeightedRefusals()
        {
            var km = PlottingPositions.PlottingPostionType.KaplanMeier;
            var hazen = PlottingPositions.PlottingPostionType.Hazen;
            Assert.ThrowsExactly<ArgumentNullException>(() => PlottingPositions.Function((IList<double>)null!, km));
            Assert.ThrowsExactly<ArgumentNullException>(() => PlottingPositions.Function((IList<double>)null!, hazen));
            Assert.ThrowsExactly<ArgumentException>(() => PlottingPositions.Function(new double[0], km));
            Assert.ThrowsExactly<ArgumentException>(() => PlottingPositions.Function(new[] { 0d, 0d }, km));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1d, -1d }, km));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1d, double.NaN }, km));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { double.MaxValue, double.MaxValue }, km));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1d, 2d }, (PlottingPositions.PlottingPostionType)42));
        }
    }
}
