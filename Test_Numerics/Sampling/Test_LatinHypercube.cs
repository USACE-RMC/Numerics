using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Sampling;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using System.Diagnostics;

namespace Sampling
{
    /// <summary>
    /// Unit tests for the Latin Hypercube class. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_LatinHypercube
    {

        /// <summary>
        /// Test for Latin Hypercube Sampling
        /// </summary>
        [TestMethod]
        public void Test_LHS()
        {
            // There is not a great way to test LHS against a known theoretical solution,
            // or other code solutions. So this test only checks that it returns the correct
            // sampling statistics.

            int N = 10000;
            var norm1 = new Normal(100, 15);
            var norm2 = new Normal(200, 30);
            var x1 = new double[N];
            var x2 = new double[N];

            // LHS
            var lhs = LatinHypercube.Random(N, 2, 45678);
            for (int i = 0; i < N; i++)
            {
                x1[i] = norm1.InverseCDF(lhs[i, 0]);
                x2[i] = norm2.InverseCDF(lhs[i, 1]);
            }

            // Test mean
            Assert.AreEqual(norm1.Mean, Statistics.Mean(x1), 1E-2);
            Assert.AreEqual(norm2.Mean, Statistics.Mean(x2), 1E-2);

            // Test standard deviation
            Assert.AreEqual(norm1.StandardDeviation, Statistics.StandardDeviation(x1), 1E-2);
            Assert.AreEqual(norm2.StandardDeviation, Statistics.StandardDeviation(x2), 1E-2);

            // Test correlation. Should be close to zero.
            Assert.AreEqual(0.0, Correlation.Pearson(x1, x2), 1E-2);

        }

        /// <summary>
        /// Test for Latin Hypercube Sampling using the median value per bin.
        /// </summary>
        [TestMethod]
        public void Test_LHS_Median()
        {
            // There is not a great way to test LHS against a known theoretical solution,
            // or other code solutions. So this test only checks that it returns the correct
            // sampling statistics.

            int N = 10000;
            var norm1 = new Normal(100, 15);
            var norm2 = new Normal(200, 30);
            var x1 = new double[N];
            var x2 = new double[N];

            // LHS
            var lhs = LatinHypercube.Median(N, 2, 45678);
            for (int i = 0; i < N; i++)
            {
                x1[i] = norm1.InverseCDF(lhs[i, 0]);
                x2[i] = norm2.InverseCDF(lhs[i, 1]);
            }

            // Test mean
            Assert.AreEqual(norm1.Mean, Statistics.Mean(x1), 1E-2);
            Assert.AreEqual(norm2.Mean, Statistics.Mean(x2), 1E-2);

            // Test standard deviation
            Assert.AreEqual(norm1.StandardDeviation, Statistics.StandardDeviation(x1), 1E-2);
            Assert.AreEqual(norm2.StandardDeviation, Statistics.StandardDeviation(x2), 1E-2);

            // Test correlation. Should be close to zero.
            // Correlation is higher with the median option
            Assert.AreEqual(0.0, Correlation.Pearson(x1, x2), 0.05);

        }


        /// <summary>
        /// The generator overloads reproduce the seeded overloads exactly when given the same seed.
        /// </summary>
        [TestMethod]
        public void Test_LHS_GeneratorOverload_MatchesSeedOverload()
        {
            var expectedRandom = LatinHypercube.Random(50, 3, 45678);
            var actualRandom = LatinHypercube.Random(50, 3, new MersenneTwister(45678));
            CollectionAssert.AreEqual(expectedRandom, actualRandom);

            var expectedMedian = LatinHypercube.Median(50, 3, 45678);
            var actualMedian = LatinHypercube.Median(50, 3, new MersenneTwister(45678));
            CollectionAssert.AreEqual(expectedMedian, actualMedian);
        }

        /// <summary>
        /// A shared generator advances one value per column, so consecutive calls are reproducible and differ.
        /// </summary>
        [TestMethod]
        public void Test_LHS_GeneratorOverload_AdvancesSharedGenerator()
        {
            var first = new MersenneTwister(2024);
            var a1 = LatinHypercube.Random(20, 2, first);
            var a2 = LatinHypercube.Random(20, 2, first);

            var second = new MersenneTwister(2024);
            var b1 = LatinHypercube.Random(20, 2, second);
            var b2 = LatinHypercube.Random(20, 2, second);

            CollectionAssert.AreEqual(a1, b1);
            CollectionAssert.AreEqual(a2, b2);
            CollectionAssert.AreNotEqual(a1, a2);

            // Two columns consume two values, so a third value is the next one in the stream.
            var reference = new MersenneTwister(2024);
            for (int i = 0; i < 4; i++) reference.Next();
            Assert.AreEqual(reference.Next(), first.Next());
        }

        /// <summary>
        /// The generator overloads keep exactly one point in each of the n strata of every column.
        /// </summary>
        [TestMethod]
        public void Test_LHS_GeneratorOverload_IsStratified()
        {
            int n = 40;
            foreach (var lhs in new[] { LatinHypercube.Random(n, 3, new MersenneTwister(7)), LatinHypercube.Median(n, 3, new MersenneTwister(7)) })
            {
                for (int col = 0; col < 3; col++)
                {
                    var used = new bool[n];
                    for (int row = 0; row < n; row++)
                    {
                        int bin = (int)(lhs[row, col] * n);
                        Assert.IsFalse(used[bin], $"Column {col} has two points in bin {bin}.");
                        used[bin] = true;
                    }
                }
            }
        }

        /// <summary>
        /// The generator overloads reject a null generator.
        /// </summary>
        [TestMethod]
        public void Test_LHS_GeneratorOverload_NullGenerator_Throws()
        {
            Assert.Throws<System.ArgumentNullException>(() => LatinHypercube.Random(10, 2, (System.Random)null));
            Assert.Throws<System.ArgumentNullException>(() => LatinHypercube.Median(10, 2, (System.Random)null));
        }

    }
}
