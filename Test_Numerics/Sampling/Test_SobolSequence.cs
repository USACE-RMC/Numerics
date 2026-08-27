using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Sampling;
using System;

namespace Sampling
{
    /// <summary>
    /// Unit test for the Sobol Sequence class. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// R Core Team (2024). _R: A Language and Environment for Statistical Computing_.R Foundation for Statistical Computing, Vienna,
    /// Austria. <see href="https://www.R-project.org/"/>
    /// </remarks>
    [TestClass]
    public class Test_SobolSequence
    {

        /// <summary>
        /// Tested against the 'sobol' method in the 'randtoolbox' R package.
        /// </summary>
        [TestMethod]
        public void Test_Sobol()
        {
            // the results from R
            var trueResult = new double[,]
            { 
                { 0.5000, 0.5000 },
                { 0.7500, 0.2500 },
                { 0.2500, 0.7500 },
                { 0.3750, 0.3750 },
                { 0.8750, 0.8750 },
                { 0.6250, 0.1250 },
                { 0.1250, 0.6250 },
                { 0.1875, 0.3125 },
                { 0.6875, 0.8125 },
                { 0.9375, 0.0625 }
            };

            var sobol = new SobolSequence(2);
            for (int i = 0; i < 10; i++)
            {
                var rnd = sobol.NextDouble();
                for (int j = 0; j < 2; j++)
                {
                    Assert.AreEqual(trueResult[i, j], rnd[j]);
                }
            }

        }

        /// <summary>
        /// The same seed reproduces the scrambled sequence bit-for-bit, and the Seed property
        /// reports the scrambling seed (null for the unrandomized sequence).
        /// </summary>
        [TestMethod]
        public void Test_ScrambledSobol_SameSeed_BitIdentical()
        {
            var first = new SobolSequence(3, 12345);
            var second = new SobolSequence(3, 12345);
            for (int i = 0; i < 100; i++)
            {
                var a = first.NextDouble();
                var b = second.NextDouble();
                for (int j = 0; j < 3; j++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(a[j]), BitConverter.DoubleToInt64Bits(b[j]));
                }
            }
            Assert.AreEqual(12345, first.Seed);
            Assert.IsNull(new SobolSequence(3).Seed);
        }

        /// <summary>
        /// Different seeds produce different scrambled sequences.
        /// </summary>
        [TestMethod]
        public void Test_ScrambledSobol_DifferentSeeds_Differ()
        {
            var first = new SobolSequence(2, 1);
            var second = new SobolSequence(2, 2);
            bool anyDifferent = false;
            for (int i = 0; i < 10 && !anyDifferent; i++)
            {
                var a = first.NextDouble();
                var b = second.NextDouble();
                for (int j = 0; j < 2; j++)
                {
                    if (a[j] != b[j]) anyDifferent = true;
                }
            }
            Assert.IsTrue(anyDifferent);
        }

        /// <summary>
        /// SkipTo on a scrambled sequence reproduces the sequential emissions under the class's
        /// historical indexing, pinned as-is: the index positions the next emission, so SkipTo(0)
        /// returns the first emitted point and SkipTo(i) for i &gt;= 1 returns the i-th emitted
        /// point (SkipTo(0) and SkipTo(1) coincide). The unrandomized sequence shares the contract.
        /// </summary>
        [TestMethod]
        public void Test_ScrambledSobol_SkipTo_MatchesSequential()
        {
            foreach (int? seed in new int?[] { null, 42 })
            {
                var sequential = seed.HasValue ? new SobolSequence(3, seed.Value) : new SobolSequence(3);
                var points = new double[8][];
                for (int i = 0; i < 8; i++)
                    points[i] = sequential.NextDouble();

                for (int i = 0; i < 8; i++)
                {
                    var skipped = seed.HasValue ? new SobolSequence(3, seed.Value) : new SobolSequence(3);
                    var point = skipped.SkipTo(i);
                    var expected = points[i == 0 ? 0 : i - 1];
                    for (int j = 0; j < 3; j++)
                    {
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[j]), BitConverter.DoubleToInt64Bits(point[j]));
                    }
                }
            }
        }

        /// <summary>
        /// Linear matrix scrambling with a digital shift preserves the base-2 (0,1)-sequence
        /// property of every dimension: each full dyadic block of 2^k consecutive points (standard
        /// indices 2^k through 2^(k+1)-1; this generator drops the origin, so those are emitted
        /// positions 2^k through 2^(k+1)-1) places exactly one point in each of the 2^k equal bins.
        /// </summary>
        [TestMethod]
        public void Test_ScrambledSobol_PreservesDyadicStratification()
        {
            foreach (int? seed in new int?[] { null, 11, 12345 })
            {
                for (int k = 1; k <= 7; k++)
                {
                    var sobol = seed.HasValue ? new SobolSequence(3, seed.Value) : new SobolSequence(3);
                    int block = 1 << k;
                    // Discard emitted positions 1 .. 2^k - 1.
                    for (int i = 1; i < block; i++)
                        sobol.NextDouble();
                    var counts = new int[3, block];
                    for (int i = 0; i < block; i++)
                    {
                        var point = sobol.NextDouble();
                        for (int j = 0; j < 3; j++)
                        {
                            int bin = (int)(point[j] * block);
                            counts[j, bin]++;
                        }
                    }
                    for (int j = 0; j < 3; j++)
                    {
                        for (int b = 0; b < block; b++)
                        {
                            Assert.AreEqual(1, counts[j, b], $"seed {seed}, k {k}, dimension {j}, bin {b}");
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Scrambled marginals are uniform: over 4096 points each dimension's mean is near 1/2 and
        /// variance near 1/12, far inside Monte Carlo tolerances at this sample size.
        /// </summary>
        [TestMethod]
        public void Test_ScrambledSobol_MarginalUniformity()
        {
            var sobol = new SobolSequence(2, 777);
            int n = 4096;
            var sums = new double[2];
            var sumSquares = new double[2];
            for (int i = 0; i < n; i++)
            {
                var point = sobol.NextDouble();
                for (int j = 0; j < 2; j++)
                {
                    sums[j] += point[j];
                    sumSquares[j] += point[j] * point[j];
                }
            }
            for (int j = 0; j < 2; j++)
            {
                double mean = sums[j] / n;
                double variance = sumSquares[j] / n - mean * mean;
                Assert.AreEqual(0.5d, mean, 0.01d);
                Assert.AreEqual(1d / 12d, variance, 0.005d);
            }
        }

        /// <summary>
        /// Scrambled quasi-Monte Carlo beats plain Monte Carlo on a smooth integrand: across twenty
        /// seeds at 1,024 points, the root-mean-square error of the scrambled-Sobol estimate of
        /// the integral of x^2*e^y over the unit square must be under half the Monte Carlo RMSE
        /// (the measured advantage is roughly an order of magnitude; one half is the guarded bound).
        /// </summary>
        [TestMethod]
        public void Test_ScrambledSobol_VarianceReductionVsMonteCarlo()
        {
            double exact = (Math.E - 1d) / 3d;
            int n = 1024;
            double sumSqQmc = 0, sumSqMc = 0;
            for (int seed = 101; seed <= 120; seed++)
            {
                var sobol = new SobolSequence(2, seed);
                double qmc = 0;
                for (int i = 0; i < n; i++)
                {
                    var p = sobol.NextDouble();
                    qmc += p[0] * p[0] * Math.Exp(p[1]);
                }
                qmc /= n;
                sumSqQmc += (qmc - exact) * (qmc - exact);

                var prng = new MersenneTwister(seed);
                double mc = 0;
                for (int i = 0; i < n; i++)
                {
                    double u1 = prng.NextDouble();
                    double u2 = prng.NextDouble();
                    mc += u1 * u1 * Math.Exp(u2);
                }
                mc /= n;
                sumSqMc += (mc - exact) * (mc - exact);
            }
            double rmseQmc = Math.Sqrt(sumSqQmc / 20d);
            double rmseMc = Math.Sqrt(sumSqMc / 20d);
            Assert.IsLessThan(0.5d * rmseMc, rmseQmc);
        }

    }
}
