using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Distributions;

namespace Utilities
{
    /// <summary>
    /// Tests for the public compensated steps <see cref="Tools.NeumaierAdd(ref double, ref double, double)"/>,
    /// <see cref="Tools.NeumaierResult(double, double)"/> and <see cref="Tools.KahanAdd(ref double, ref double, double)"/>.
    /// </summary>
    /// <remarks>
    /// The steps must reproduce the compensated sums bit for bit, and the quantile variance contraction that uses them must
    /// reproduce its earlier inline Neumaier loop bit for bit.
    /// </remarks>
    [TestClass]
    public class Test_CompensatedStep
    {
        /// <summary>Random lists spanning many decades with both signs, plus non-finite cases.</summary>
        private static IEnumerable<double[]> Lists()
        {
            var random = new Random(20261108);
            for (int trial = 0; trial < 200; trial++)
            {
                var values = new double[random.Next(1, 300)];
                for (int i = 0; i < values.Length; i++)
                    values[i] = (random.NextDouble() - 0.5) * Math.Pow(10, 40 * random.NextDouble() - 20);
                yield return values;
            }
            yield return new[] { 1E100, 1d, -1E100 };
            yield return new[] { 1d, 1E100, 1d, -1E100 };
            yield return new[] { 1d, double.PositiveInfinity, -1E308 };
            yield return new[] { 1d, double.NaN, 2d };
            yield return new[] { double.PositiveInfinity, 1d, double.NegativeInfinity };
            yield return new[] { 1E308, 1E308, -1E308 };
        }

        /// <summary>The Neumaier step and result reproduce <see cref="Tools.CompensatedSum(IList{double})"/> and each entry
        /// of <see cref="Tools.CompensatedCumulativeSum(IList{double})"/> bit for bit.</summary>
        [TestMethod]
        public void NeumaierStep_BitwiseCompensatedSum()
        {
            foreach (var values in Lists())
            {
                double sum = 0d, compensation = 0d;
                var cumulative = Tools.CompensatedCumulativeSum(values);
                for (int i = 0; i < values.Length; i++)
                {
                    Tools.NeumaierAdd(ref sum, ref compensation, values[i]);
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(cumulative[i]), BitConverter.DoubleToInt64Bits(Tools.NeumaierResult(sum, compensation)));
                }
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.CompensatedSum(values)), BitConverter.DoubleToInt64Bits(Tools.NeumaierResult(sum, compensation)));
            }
        }

        /// <summary>The Kahan step reproduces <see cref="Tools.CompensatedSum(IList{double}, bool)"/> with kahan true and
        /// the classic loop bit for bit.</summary>
        [TestMethod]
        public void KahanStep_BitwiseCompensatedSum()
        {
            foreach (var values in Lists())
            {
                double sum = 0d, compensation = 0d, t = 0d, c = 0d;
                var cumulative = Tools.CompensatedCumulativeSum(values, true);
                for (int i = 0; i < values.Length; i++)
                {
                    Tools.KahanAdd(ref sum, ref compensation, values[i]);
                    double a = values[i] - c, f = t + a;
                    c = (f - t) - a;
                    t = f;
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(t), BitConverter.DoubleToInt64Bits(sum));
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(c), BitConverter.DoubleToInt64Bits(compensation));
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(cumulative[i]), BitConverter.DoubleToInt64Bits(sum));
                }
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Tools.CompensatedSum(values, true)), BitConverter.DoubleToInt64Bits(sum));
            }
        }

        /// <summary>Known values: Neumaier recovers 2 for [1, 1e100, 1, -1e100] where Kahan gives 0; a non-finite sum is
        /// returned as the IEEE sum.</summary>
        [TestMethod]
        public void Steps_KnownValues()
        {
            double sum = 0d, compensation = 0d;
            foreach (double v in new[] { 1d, 1E100, 1d, -1E100 }) Tools.NeumaierAdd(ref sum, ref compensation, v);
            Assert.AreEqual(0d, sum, 0d);
            Assert.AreEqual(2d, Tools.NeumaierResult(sum, compensation), 0d);
            double ks = 0d, kc = 0d;
            foreach (double v in new[] { 1d, 1E100, 1d, -1E100 }) Tools.KahanAdd(ref ks, ref kc, v);
            Assert.AreEqual(0d, ks, 0d);
            Assert.AreEqual(double.PositiveInfinity, Tools.NeumaierResult(double.PositiveInfinity, double.NaN));
            Assert.IsTrue(double.IsNaN(Tools.NeumaierResult(double.NaN, 1d)));
            Assert.AreEqual(3d, Tools.NeumaierResult(2d, 1d), 0d);
        }

        /// <summary>
        /// The quantile variance contraction is bitwise the earlier inline Neumaier loop over random covariance matrices
        /// and gradients with mixed signs, magnitudes and zero entries.
        /// </summary>
        [TestMethod]
        public void ScaledQuantileVariance_BitwiseInlineLoop()
        {
            var random = new Random(20261108);
            for (int trial = 0; trial < 300; trial++)
            {
                int n = random.Next(1, 6);
                var a = new double[n, n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        a[i, j] = (random.NextDouble() - 0.5) * Math.Pow(10, 6 * random.NextDouble() - 3);
                var covariance = new double[n, n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        for (int k = 0; k < n; k++)
                            covariance[i, j] += a[i, k] * a[j, k];
                var gradient = new double[n];
                for (int i = 0; i < n; i++)
                    gradient[i] = random.Next(5) == 0 ? 0d : (random.NextDouble() - 0.5) * Math.Pow(10, 8 * random.NextDouble() - 4);
                double scale = Math.Pow(10, 4 * random.NextDouble() - 2);
                double expected = InlineReference(covariance, gradient, scale);
                double actual = DistributionNumerics.ScaledQuantileVariance(covariance, gradient, scale);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(actual));
            }
        }

        /// <summary>The contraction as written before the public step, kept verbatim as the bitwise reference.</summary>
        private static double InlineReference(double[,] covariance, double[] gradient, double scale)
        {
            var gradientLogs = new double[gradient.Length];
            for (int i = 0; i < gradient.Length; i++) gradientLogs[i] = Math.Log(Math.Abs(gradient[i]));
            double largestLogTerm = double.NegativeInfinity;
            for (int i = 0; i < gradient.Length; i++)
                for (int j = 0; j < gradient.Length; j++)
                {
                    if (covariance[i, j] == 0 || gradient[i] == 0 || gradient[j] == 0) continue;
                    largestLogTerm = Math.Max(largestLogTerm, Math.Log(Math.Abs(covariance[i, j])) + gradientLogs[i] + gradientLogs[j]);
                }
            if (double.IsNegativeInfinity(largestLogTerm)) return 0;
            double quadratic = 0, correction = 0;
            for (int i = 0; i < gradient.Length; i++)
                for (int j = 0; j < gradient.Length; j++)
                {
                    if (covariance[i, j] == 0 || gradient[i] == 0 || gradient[j] == 0) continue;
                    double logTerm = Math.Log(Math.Abs(covariance[i, j])) + gradientLogs[i] + gradientLogs[j];
                    double term = Math.Sign(covariance[i, j]) * Math.Sign(gradient[i]) * Math.Sign(gradient[j])
                        * Math.Exp(logTerm - largestLogTerm);
                    double next = quadratic + term;
                    correction += Math.Abs(quadratic) >= Math.Abs(term) ? (quadratic - next) + term : (term - next) + quadratic;
                    quadratic = next;
                }
            quadratic += correction;
            if (!Tools.IsFinite(quadratic) || quadratic < 0) throw new InvalidOperationException();
            if (quadratic == 0) return 0;
            return Math.Exp(2 * Math.Log(scale) + largestLogTerm + Math.Log(quadratic));
        }
    }
}
