using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;

namespace Data.Statistics
{
    /// <summary>
    /// Bit-equality tests for the element-wise running covariance update and the vectorized correlation matrix.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// Each test compares the library result with a reference that performs the same operations in the same order,
    /// on data containing −0.0, NaN and infinities. Results must agree bit for bit, with any NaN matching any NaN.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_VectorizedStatistics
    {
        private static void AssertSame(double expected, double actual, string where)
        {
            if (double.IsNaN(expected) && double.IsNaN(actual)) return;
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(actual), $"{where}: expected {expected:R}, actual {actual:R}");
        }

        /// <summary>
        /// The running covariance update matches the matrix-expression Welford update bit for bit, including −0.0
        /// inputs, a NaN and an infinity, and earlier references to the mean and co-moment matrices keep their values.
        /// </summary>
        [TestMethod]
        public void Test_RunningCovarianceMatrix_BitIdentical()
        {
            foreach (int size in new[] { 1, 3, 9 })
            {
                var rng = new Random(size);
                var running = new RunningCovarianceMatrix(size);
                var mean = new Matrix(size, 1);
                var covariance = Matrix.Identity(size);
                for (int step = 1; step <= 40; step++)
                {
                    var values = new double[size];
                    for (int i = 0; i < size; i++)
                    {
                        values[i] = (rng.NextDouble() - 0.5) * 10;
                        if (step == 1) values[i] = -0.0;
                        if (step == 7 && i == 0) values[i] = 0.0;
                    }
                    if (step == 30) values[size - 1] = double.NaN;
                    if (step == 35) values[0] = double.PositiveInfinity;

                    // The matrix-expression update.
                    var x = new Matrix(values);
                    var oldMean = mean;
                    mean += step == 1 ? x : (x - mean) * (1d / step);
                    covariance += (x - oldMean) * Matrix.Transpose(x - mean);

                    var previousMean = running.Mean;
                    var previousCovariance = running.Covariance;
                    var copyMean = previousMean.Clone();
                    running.Push(values);
                    Assert.AreEqual(step, running.N);
                    for (int i = 0; i < size; i++)
                    {
                        AssertSame(mean[i, 0], running.Mean[i, 0], $"mean {size}/{step}/{i}");
                        AssertSame(copyMean[i, 0], previousMean[i, 0], "earlier mean unchanged");
                        for (int j = 0; j < size; j++)
                            AssertSame(covariance[i, j], running.Covariance[i, j], $"cov {size}/{step}/{i},{j}");
                    }
                    Assert.AreNotSame(previousCovariance, running.Covariance);
                }
            }
            var r = new RunningCovarianceMatrix(3);
            Assert.ThrowsExactly<ArgumentException>(() => r.Push(new double[2]));
        }

        /// <summary>
        /// The scalar correlation-matrix loops of <see cref="Correlation.Pearson(double[,])"/>.
        /// </summary>
        private static double[,] ReferencePearson(double[,] samples)
        {
            int n = samples.GetLength(0), p = samples.GetLength(1);
            var means = new double[p];
            for (int j = 0; j < p; j++) for (int i = 0; i < n; i++) means[j] += samples[i, j];
            for (int j = 0; j < p; j++) means[j] /= n;
            var ss = new double[p];
            var cov = new double[p, p];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < p; j++)
                {
                    double dx = samples[i, j] - means[j];
                    ss[j] += dx * dx;
                    for (int k = j; k < p; k++) cov[j, k] += dx * (samples[i, k] - means[k]);
                }
            for (int j = 0; j < p; j++) for (int k = j + 1; k < p; k++) cov[k, j] = cov[j, k];
            var corr = new double[p, p];
            for (int j = 0; j < p; j++) for (int k = 0; k < p; k++) corr[j, k] = cov[j, k] / Math.Sqrt(ss[j] * ss[k]);
            return corr;
        }

        /// <summary>
        /// The Pearson and Spearman correlation matrices match the scalar loops bit for bit below and above the
        /// vectorization width, including −0.0, a constant column (0/0) and a NaN.
        /// </summary>
        [TestMethod]
        public void Test_PearsonMatrix_BitIdentical()
        {
            foreach (var (n, p) in new[] { (30, 3), (50, 7), (64, 8), (101, 23) })
            {
                var rng = new Random(n + p);
                var s = new double[n, p];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < p; j++) s[i, j] = (rng.NextDouble() - 0.5) * Math.Pow(10, j % 5) + j;
                s[3, 0] = -0.0;
                for (int i = 0; i < n; i++) s[i, p - 1] = 2.5;
                var expected = ReferencePearson(s);
                var actual = Correlation.Pearson(s);
                for (int j = 0; j < p; j++)
                    for (int k = 0; k < p; k++) AssertSame(expected[j, k], actual[j, k], $"Pearson {n}x{p} [{j},{k}]");

                s[5, 1] = double.NaN;
                expected = ReferencePearson(s);
                actual = Correlation.Pearson(s);
                for (int j = 0; j < p; j++)
                    for (int k = 0; k < p; k++) AssertSame(expected[j, k], actual[j, k], $"Pearson NaN {n}x{p} [{j},{k}]");
            }
        }
    }
}
