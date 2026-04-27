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
            // or other code solutions. So here I am just testing that I get back the correct
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
            // or other code solutions. So here I am just testing that I get back the correct
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

    }
}
