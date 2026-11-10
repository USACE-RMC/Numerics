using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the RunningStatistics class. Most of these methods were tested against the same validation methods used in testing the Statistics
    /// class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_RunningStatistics
    {
        /// <summary>
        /// Test the construction of the running covariance matrix.
        /// </summary>
        [TestMethod]
        public void Test_RunningMatrix()
        {
            var matrix = new RunningCovarianceMatrix(5);
            var mean = matrix.Mean;
            var covariance = matrix.Covariance;

            var trueMean = new Matrix(new double[,] { { 0 }, { 0 }, { 0 }, { 0 }, { 0 } });
            var trueCovariance = new Matrix(new double[,] { { 1, 0, 0, 0, 0 },
                                                            { 0, 1, 0, 0, 0 },
                                                            { 0, 0, 1, 0, 0 },
                                                            { 0, 0, 0, 1, 0 },
                                                            { 0, 0, 0, 0, 1 } });

            for(int i = 0; i < mean.NumberOfRows; i++)
            {
                CollectionAssert.AreEqual(trueMean.Row(i), mean.Row(i));
                CollectionAssert.AreEqual(trueCovariance.Row(i), covariance.Row(i));
            }
        }

        /// <summary>
        /// Test the push method of the running covariance matrix.
        /// </summary>
        [TestMethod]
        public void Test_RunningMatrix_Push()
        {
            var matrix = new RunningCovarianceMatrix(5);
            var values = new double[] { 1, 2, 3, 4, 5 };

            matrix.Push(values);
            var mean = matrix.Mean;
            var covariance = matrix.Covariance;

            var trueMean = new Matrix(new double[,] { { 1 }, { 2 }, { 3 }, { 4 }, { 5 } });
            var trueCovariance = new Matrix(new double[,] { { 1, 0, 0, 0, 0 },
                                                            { 0, 1, 0, 0, 0 },
                                                            { 0, 0, 1, 0, 0 },
                                                            { 0, 0, 0, 1, 0 },
                                                            { 0, 0, 0, 0, 1 } });

            for (int i = 0; i < mean.NumberOfRows; i++)
            {
                CollectionAssert.AreEqual(trueMean.Row(i), mean.Row(i));
                CollectionAssert.AreEqual(trueCovariance.Row(i), covariance.Row(i));
            }
        }

        /// <summary>
        /// The documented regularization: SampleCovariance is (M2 + I) / (N - 1), so removing I / (N - 1) leaves
        /// numpy.cov of (1, 2), (2, 1), (4, 5), (3, 3): [[5/3, 5.5/3], [5.5/3, 8.75/3]].
        /// </summary>
        [TestMethod]
        public void Test_RunningMatrix_IdentityPrior()
        {
            var matrix = new RunningCovarianceMatrix(2);
            matrix.Push(new double[] { 1, 2 });
            matrix.Push(new double[] { 2, 1 });
            matrix.Push(new double[] { 4, 5 });
            matrix.Push(new double[] { 3, 3 });
            var sample = matrix.SampleCovariance;
            var numpy = new double[,] { { 5d / 3d, 5.5 / 3d }, { 5.5 / 3d, 8.75 / 3d } };
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    Assert.AreEqual(numpy[i, j] + (i == j ? 1d / 3d : 0d), sample[i, j], 1E-12);
        }

        /// <summary>
        /// Test the push method of the running covariance matrix with mulitple pushes.
        /// </summary>
        [TestMethod]
        public void Test_RunningMatrix_MultiplePush()
        {
            var matrix = new RunningCovarianceMatrix(5);
            matrix.Push(new double[] { 1, 1, 1, 1, 1 });
            matrix.Push(new double[] { 2, 2, 2, 2, 2 });
            matrix.Push(new double[] { 3, 3, 3, 3, 3 });
            matrix.Push(new double[] { 4, 4, 4, 4, 4 });
            matrix.Push(new double[] { 5, 5, 5, 5, 5 });

            var mean = matrix.Mean;
            var covariance = matrix.Covariance;

            var trueMean = new Matrix(new double[,] { { 3 }, { 3 }, { 3 }, { 3 }, { 3 } });
            var trueCovariance = new Matrix(new double[,] { { 11, 10, 10, 10, 10 },
                                                            { 10, 11, 10, 10, 10 },
                                                            { 10, 10, 11, 10, 10 },
                                                            { 10, 10, 10, 11, 10 },
                                                            { 10, 10, 10, 10, 11 } });

            for (int i = 0; i < mean.NumberOfRows; i++)
            {
                CollectionAssert.AreEqual(trueMean.Row(i), mean.Row(i));
                CollectionAssert.AreEqual(trueCovariance.Row(i), covariance.Row(i));
            }
        }

        /// <summary>
        /// Test the construction of the running statistic class.
        /// </summary>
        [TestMethod]
        public void Test_RunningStat_Construction()
        {
            var runningstat1 = new RunningStatistics();

            var values = new double[] { 1, 2, 3, 4, 5 };
            var runningstat2 = new RunningStatistics(values);

            var test1 = runningstat1.Count;
            Assert.AreEqual(0, test1);

            var test2 = runningstat2.Count;
            Assert.AreEqual(5, test2);
        }

        /// <summary>
        /// Test the sample methods of the running statistics class. These methods were validated against the values found in the testing of the
        /// Statistics class.
        /// </summary>
        [TestMethod]
        public void Test_Stats()
        {
            var sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166, 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var runningstat = new RunningStatistics(sample1);

            var min = runningstat.Minimum;
            Assert.AreEqual(98.8, min, 1E-10);

            var max = runningstat.Maximum;
            Assert.AreEqual(337.0, max, 1E-10);

            var mean = runningstat.Mean;
            Assert.AreEqual(191.317391304348, mean, 1E-10);

            var variance = runningstat.Variance;
            Assert.AreEqual(2300.31616368286, variance, 1E-10);

            var sd = runningstat.StandardDeviation;
            Assert.AreEqual(47.9616113541118, sd, 1E-10);

            var coeffOfVariation = runningstat.CoefficientOfVariation;
            Assert.AreEqual(0.250691330396694, coeffOfVariation, 1E-10);

            var skewness = runningstat.Skewness;
            Assert.AreEqual(0.8605451107461, skewness, 1E-10);

            var kurtosis = runningstat.Kurtosis;
            Assert.AreEqual(1.3434868130194, kurtosis, 1E-10);
        }

        /// <summary>
        /// Test the population method of the running statistics class. These methods were validated against the values found in the testing of the
        /// Statistics class, or through the direct equation for the value calculated using R.
        /// </summary>
        [TestMethod]
        public void Test_PopStats()
        {
            var sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166, 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var runningstat = new RunningStatistics(sample1);

            var popVariance = runningstat.PopulationVariance;
            Assert.AreEqual(2266.97824826717, popVariance, 1E-10);

            var popSD = runningstat.PopulationStandardDeviation;
            Assert.AreEqual(47.6127950058298, popSD, 1E-10);

            // sqrt(n)* sum((sample1-mean(sample1))**3) / sum((sample1-mean(sample1))**2)**(3/2)  -- R eqm
            var popSkewness = runningstat.PopulationSkewness;
            Assert.AreEqual(0.84172348076203, popSkewness, 1E-10);

            // (sum((sample1-mean(sample1))**4)/n) / (sum((sample1-mean(sample1))**2)/n)**2 - 3  -- R eqn
            var popKurtosis = runningstat.PopulationKurtosis;
            Assert.AreEqual(1.16237367377557, popKurtosis, 1E-10);
        }

        /// <summary>
        /// Test the push methods of the running statistics class.
        /// </summary>
        [TestMethod]
        public void Test_Push()
        {
            var values = new double[] { 1, 2, 3, 4, 5 };
            var runningstat = new RunningStatistics(values);

            // push new min
            runningstat.Push(-1);
            var newMin = runningstat.Minimum;
            Assert.AreEqual(-1, newMin);

            // push new max
            runningstat.Push(6);
            var newMax = runningstat.Maximum;
            Assert.AreEqual(6, newMax);

            // push new max and min in a list
            runningstat.Push(new double[] { -5, 2.5, 3.33, 12, 20 });
            var newMinList = runningstat.Minimum;
            Assert.AreEqual(-5, newMinList);
            var newMaxList = runningstat.Maximum;
            Assert.AreEqual(20, newMaxList);
        }

        /// <summary>
        /// Test the combine method of the running statistics class.
        /// </summary>
        [TestMethod]
        public void Test_Combine()
        {
            var sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166 };
            var runningstat1 = new RunningStatistics(sample1);

            var sample2 = new double[] { 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var runningstat2 = new RunningStatistics(sample2);

            var complete = RunningStatistics.Combine(runningstat1, runningstat2);

            // values from previous tests
            var min = complete.Minimum;
            Assert.AreEqual(98.8, min, 1E-10);

            var max = complete.Maximum;
            Assert.AreEqual(337.0, max, 1E-10);

            var mean = complete.Mean;
            Assert.AreEqual(191.317391304348, mean, 1E-10);

            var variance = complete.Variance;
            Assert.AreEqual(2300.31616368286, variance, 1E-10);

            var sd = complete.StandardDeviation;
            Assert.AreEqual(47.9616113541118, sd, 1E-10);

            var coeffOfVariation = complete.CoefficientOfVariation;
            Assert.AreEqual(0.250691330396694, coeffOfVariation, 1E-10);

            var skewness = complete.Skewness;
            Assert.AreEqual(0.8605451107461, skewness, 1E-10);

            var kurtosis = complete.Kurtosis;
            Assert.AreEqual(1.3434868130194, kurtosis, 1E-10);
        }

        /// <summary>
        /// Test the overload of the + operator for the running statistics class.
        /// </summary>
        [TestMethod]
        public void Test_Add()
        {
            var sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166 };
            var runningstat1 = new RunningStatistics(sample1);

            var sample2 = new double[] { 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var runningstat2 = new RunningStatistics(sample2);

            var complete = runningstat1 + runningstat2;

            // values from previous tests
            var min = complete.Minimum;
            Assert.AreEqual(98.8, min, 1E-10);

            var max = complete.Maximum;
            Assert.AreEqual(337.0, max, 1E-10);

            var mean = complete.Mean;
            Assert.AreEqual(191.317391304348, mean, 1E-10);

            var variance = complete.Variance;
            Assert.AreEqual(2300.31616368286, variance, 1E-10);

            var sd = complete.StandardDeviation;
            Assert.AreEqual(47.9616113541118, sd, 1E-10);

            var coeffOfVariation = complete.CoefficientOfVariation;
            Assert.AreEqual(0.250691330396694, coeffOfVariation, 1E-10);

            var skewness = complete.Skewness;
            Assert.AreEqual(0.8605451107461, skewness, 1E-10);

            var kurtosis = complete.Kurtosis;
            Assert.AreEqual(1.3434868130194, kurtosis, 1E-10);
        }

        /// <summary>
        /// Test the clone method of the running statistics class. The clone should match the original statistics
        /// and remain independent of the original after either instance is updated.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var values = new double[] { 1, 2, 3, 4, 5 };
            var runningstat = new RunningStatistics(values);
            var clone = runningstat.Clone();

            // The clone should be a new instance that matches the original statistics.
            Assert.AreNotSame(runningstat, clone);
            Assert.AreEqual(5, clone.Count);
            Assert.AreEqual(1, clone.Minimum);
            Assert.AreEqual(5, clone.Maximum);
            Assert.AreEqual(3, clone.Mean, 1E-10);
            Assert.AreEqual(2.5, clone.Variance, 1E-10);

            // Pushing to the original should not change the clone.
            runningstat.Push(new double[] { -10, 30 });
            Assert.AreEqual(5, clone.Count);
            Assert.AreEqual(1, clone.Minimum);
            Assert.AreEqual(5, clone.Maximum);
            Assert.AreEqual(3, clone.Mean, 1E-10);
            Assert.AreEqual(2.5, clone.Variance, 1E-10);

            // Pushing to the clone should not change the original.
            var originalCount = runningstat.Count;
            var originalMean = runningstat.Mean;
            clone.Push(100);
            Assert.AreEqual(originalCount, runningstat.Count);
            Assert.AreEqual(originalMean, runningstat.Mean, 1E-10);
        }

        /// <summary>
        /// Test that combining with an empty running statistics returns a new instance
        /// rather than aliasing the non-empty operand.
        /// </summary>
        [TestMethod]
        public void Test_Combine_WithEmpty()
        {
            var values = new double[] { 1, 2, 3, 4, 5 };
            var stats = new RunningStatistics(values);

            // Empty first operand. The result should match the non-empty operand.
            var combined = RunningStatistics.Combine(new RunningStatistics(), stats);
            Assert.AreNotSame(stats, combined);
            Assert.AreEqual(5, combined.Count);
            Assert.AreEqual(1, combined.Minimum);
            Assert.AreEqual(5, combined.Maximum);
            Assert.AreEqual(3, combined.Mean, 1E-10);
            Assert.AreEqual(2.5, combined.Variance, 1E-10);

            // Mutating the result should not mutate the operand.
            combined.Push(1000);
            Assert.AreEqual(5, stats.Count);
            Assert.AreEqual(3, stats.Mean, 1E-10);

            // Empty second operand.
            var combined2 = RunningStatistics.Combine(stats, new RunningStatistics());
            Assert.AreNotSame(stats, combined2);
            combined2.Push(1000);
            Assert.AreEqual(5, stats.Count);
            Assert.AreEqual(3, stats.Mean, 1E-10);

            // The + operator should behave the same way.
            var combined3 = new RunningStatistics() + stats;
            Assert.AreNotSame(stats, combined3);
            combined3.Push(1000);
            Assert.AreEqual(5, stats.Count);
            Assert.AreEqual(3, stats.Mean, 1E-10);

            // Combining two empty instances should return a new empty instance.
            var empty1 = new RunningStatistics();
            var empty2 = new RunningStatistics();
            var combinedEmpty = RunningStatistics.Combine(empty1, empty2);
            Assert.AreNotSame(empty1, combinedEmpty);
            Assert.AreNotSame(empty2, combinedEmpty);
            Assert.AreEqual(0, combinedEmpty.Count);
        }

        /// <summary>
        /// A single value has population variance and standard deviation 0, as Statistics.PopulationVariance gives;
        /// an empty instance gives NaN and the sample variance of one value stays NaN.
        /// </summary>
        [TestMethod]
        public void Test_PopulationVariance_SingleValue()
        {
            var one = new RunningStatistics(new[] { 5.0 });
            Assert.AreEqual(0d, one.PopulationVariance);
            Assert.AreEqual(0d, one.PopulationStandardDeviation);
            Assert.AreEqual(Numerics.Data.Statistics.Statistics.PopulationVariance(new[] { 5.0 }), one.PopulationVariance);
            Assert.IsTrue(double.IsNaN(one.Variance));
            var empty = new RunningStatistics();
            Assert.IsTrue(double.IsNaN(empty.PopulationVariance));
            Assert.IsTrue(double.IsNaN(empty.PopulationStandardDeviation));
        }

        /// <summary>
        /// Combining two accumulators whose total exceeds 2,097,151 values (where n³ no longer fits a long)
        /// against the analytic moments of an equal two-point sample: 1.5M zeros and 1.5M threes have population
        /// excess kurtosis exactly -2 and population skewness 0.
        /// </summary>
        [TestMethod]
        public void Test_Combine_LargeCounts()
        {
            var a = new RunningStatistics();
            var b = new RunningStatistics();
            for (int i = 0; i < 1_500_000; i++)
            {
                a.Push(0d);
                b.Push(3d);
            }
            var c = RunningStatistics.Combine(a, b);
            double n = 3_000_000d;
            Assert.AreEqual(3_000_000L, c.Count);
            Assert.AreEqual(1.5, c.Mean, 1E-12);
            Assert.AreEqual(2.25, c.PopulationVariance, 1E-12);
            Assert.AreEqual(0d, c.PopulationSkewness, 1E-12);
            Assert.AreEqual(-2d, c.PopulationKurtosis, 1E-9);
            // Sample excess kurtosis from n m4 / m2² = 1: (n² - 1) / ((n - 2)(n - 3)) (1 - 3 + 6 / (n + 1)).
            double g2 = (n * n - 1d) / ((n - 2d) * (n - 3d)) * (-2d + 6d / (n + 1d));
            Assert.AreEqual(g2, c.Kurtosis, 1E-9);
        }

        /// <summary>
        /// Clear returns the instance to the empty state, and pushing values afterwards reproduces a fresh
        /// instance bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_Clear()
        {
            var stats = new RunningStatistics(new double[] { 100, -3, 7.5, 42, 1E6 });
            stats.Clear();
            Assert.AreEqual(0L, stats.Count);
            Assert.IsTrue(double.IsNaN(stats.Mean));
            Assert.IsTrue(double.IsNaN(stats.Minimum));
            Assert.IsTrue(double.IsNaN(stats.Maximum));
            Assert.IsTrue(double.IsNaN(stats.Variance));

            var values = new double[] { 1.25, 3.5, -2.75, 8.125, 4.0, 0.3 };
            stats.Push(values);
            var fresh = new RunningStatistics(values);
            Assert.AreEqual(fresh.Count, stats.Count);
            Assert.AreEqual(fresh.Minimum, stats.Minimum);
            Assert.AreEqual(fresh.Maximum, stats.Maximum);
            Assert.AreEqual(fresh.Mean, stats.Mean);
            Assert.AreEqual(fresh.Variance, stats.Variance);
            Assert.AreEqual(fresh.Skewness, stats.Skewness);
            Assert.AreEqual(fresh.Kurtosis, stats.Kurtosis);

            // A cleared instance combines like an empty one.
            stats.Clear();
            var combined = RunningStatistics.Combine(stats, fresh);
            Assert.AreEqual(fresh.Count, combined.Count);
            Assert.AreEqual(fresh.Kurtosis, combined.Kurtosis);
        }
    }
}
