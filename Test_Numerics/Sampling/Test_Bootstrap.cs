using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace Sampling
{

    /// <summary>
    /// Unit tests for the general-purpose Bootstrap class.
    /// Tests compare against the 'boot' R package reference values and Monte Carlo CIs,
    /// following the same validation approach as Test_BootstrapAnalysis.
    /// </summary>
    [TestClass]
    public class Test_Bootstrap
    {

        private double _mu = 3.122599;
        private double _sigma = 0.5573654;
        private int _sampleSize = 100;
        private double[] _probabilities = new double[] { 0.999, 0.99, 0.95, 0.9, 0.5, 0.1, 0.05, 0.01 };

        /// <summary>
        /// Creates a generic bootstrap configured to replicate parametric bootstrap of a Normal distribution
        /// using the method of moments. This mirrors the BootstrapAnalysis setup used in Test_BootstrapAnalysis.
        /// </summary>
        private Bootstrap<double[]> CreateNormalBootstrap()
        {
            // The "original data" is not used directly for parametric bootstrap;
            // resampling is driven by the parameters. We pass null.
            var parms = new ParameterSet(new double[] { _mu, _sigma }, double.NaN);
            var boot = new Bootstrap<double[]>(null, parms);
            boot.Replicates = 10000;
            boot.PRNGSeed = 12345;

            // Resample: generate random values from the Normal distribution defined by the parameter set
            boot.ResampleFunction = (data, ps, rng) =>
            {
                var d = new Normal(ps.Values[0], ps.Values[1]);
                return d.GenerateRandomValues(_sampleSize, rng.Next());
            };

            // Fit: estimate Normal parameters from sample using method of moments
            boot.FitFunction = (sample) =>
            {
                var d = new Normal();
                ((IEstimation)d).Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
                if (!d.ParametersValid)
                    throw new Exception("Invalid parameters.");
                return new ParameterSet(d.GetParameters, double.NaN);
            };

            // Statistic: compute quantiles at the specified probabilities
            boot.StatisticFunction = (ps) =>
            {
                var d = new Normal(ps.Values[0], ps.Values[1]);
                var result = new double[_probabilities.Length];
                for (int i = 0; i < _probabilities.Length; i++)
                    result[i] = d.InverseCDF(_probabilities[i]);
                return result;
            };

            return boot;
        }

        /// <summary>
        /// Test that the percentile method produces results consistent with the 'boot' R package.
        /// Reference values from Test_BootstrapAnalysis.Test_PercentileCI().
        /// </summary>
        [TestMethod]
        public void Test_PercentileCI()
        {
            var boot = CreateNormalBootstrap();
            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile);

            // Verify we got results for all probabilities
            Assert.HasCount(_probabilities.Length, results.StatisticResults);

            // Verify valid count
            for (int i = 0; i < _probabilities.Length; i++)
            {
                Assert.IsGreaterThan(0.95 * boot.Replicates, results.StatisticResults[i].ValidCount,
                    $"Too many failures at probability {_probabilities[i]}");
            }

            // Compare against BootstrapAnalysis
            var dist = new Normal(_mu, _sigma);
            var ba = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, _sampleSize);
            var baCIs = ba.PercentileQuantileCI(_probabilities);

            for (int i = 0; i < _probabilities.Length; i++)
            {
                Assert.AreEqual(baCIs[i, 0], results.StatisticResults[i].LowerCI, 0.01 * Math.Abs(baCIs[i, 0]),
                    $"Lower CI mismatch at p={_probabilities[i]}");
                Assert.AreEqual(baCIs[i, 1], results.StatisticResults[i].UpperCI, 0.01 * Math.Abs(baCIs[i, 1]),
                    $"Upper CI mismatch at p={_probabilities[i]}");
            }
        }

        /// <summary>
        /// Test that the Normal (standard) method produces results consistent with the 'boot' R package.
        /// Reference values from Test_BootstrapAnalysis.Test_NormalCI().
        /// </summary>
        [TestMethod]
        public void Test_NormalCI()
        {
            var boot = CreateNormalBootstrap();
            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.Normal);

            // Compare against BootstrapAnalysis
            var dist = new Normal(_mu, _sigma);
            var ba = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, _sampleSize);
            var baCIs = ba.NormalQuantileCI(_probabilities);

            for (int i = 0; i < _probabilities.Length; i++)
            {
                Assert.AreEqual(baCIs[i, 0], results.StatisticResults[i].LowerCI, 0.01 * Math.Abs(baCIs[i, 0]),
                    $"Lower CI mismatch at p={_probabilities[i]}");
                Assert.AreEqual(baCIs[i, 1], results.StatisticResults[i].UpperCI, 0.01 * Math.Abs(baCIs[i, 1]),
                    $"Upper CI mismatch at p={_probabilities[i]}");
            }
        }

        /// <summary>
        /// Test that the bias-corrected (BC) method produces results consistent with BootstrapAnalysis.
        /// </summary>
        [TestMethod]
        public void Test_BiasCorrectedCI()
        {
            var boot = CreateNormalBootstrap();
            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.BiasCorrected);

            // Compare against BootstrapAnalysis
            var dist = new Normal(_mu, _sigma);
            var ba = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, _sampleSize);
            var baCIs = ba.BiasCorrectedQuantileCI(_probabilities);

            for (int i = 0; i < _probabilities.Length; i++)
            {
                Assert.AreEqual(baCIs[i, 0], results.StatisticResults[i].LowerCI, 0.01 * Math.Abs(baCIs[i, 0]),
                    $"Lower CI mismatch at p={_probabilities[i]}");
                Assert.AreEqual(baCIs[i, 1], results.StatisticResults[i].UpperCI, 0.01 * Math.Abs(baCIs[i, 1]),
                    $"Upper CI mismatch at p={_probabilities[i]}");
            }
        }

        /// <summary>
        /// Test that the Bootstrap-t (studentized) method produces results consistent with the "true"
        /// Monte Carlo confidence intervals for the Normal distribution.
        /// </summary>
        [TestMethod]
        public void Test_BootstrapTCI()
        {
            var boot = CreateNormalBootstrap();
            boot.RunWithStudentizedBootstrap();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.BootstrapT);

            // Compare against the true Monte Carlo CIs (same approach as Test_BootstrapAnalysis.Test_BootstrapTCI)
            var dist = new Normal(_mu, _sigma);
            var trueCIs = dist.MonteCarloConfidenceIntervals(_sampleSize, 10000, _probabilities, new double[] { 0.05, 0.95 });

            for (int i = 0; i < _probabilities.Length; i++)
            {
                Assert.AreEqual(trueCIs[i, 0], results.StatisticResults[i].LowerCI, 0.01 * Math.Abs(trueCIs[i, 0]),
                    $"Lower CI mismatch at p={_probabilities[i]}");
                Assert.AreEqual(trueCIs[i, 1], results.StatisticResults[i].UpperCI, 0.01 * Math.Abs(trueCIs[i, 1]),
                    $"Upper CI mismatch at p={_probabilities[i]}");
            }
        }

        /// <summary>
        /// Test that the BCa method produces results consistent with the "true"
        /// Monte Carlo confidence intervals for the Normal distribution.
        /// </summary>
        [TestMethod]
        public void Test_BCaCI()
        {
            // Use the same sample data as Test_BootstrapAnalysis.Test_BCaCI
            var sampleData = new double[] { 3.292764, 3.354733, 2.945348, 2.773251, 3.302944, 2.091022, 3.315049, 2.861908, 2.85792, 2.540339, 2.941876, 3.908656, 3.185314, 3.260108, 2.624734, 3.40845, 2.556821, 2.834211, 3.560356, 3.149362, 3.389811, 3.727893, 2.677836, 2.223431, 2.201145, 3.902549, 2.759176, 3.31019, 3.306062, 2.918845, 3.405937, 4.098417, 4.024595, 3.816223, 3.127136, 3.245594, 2.837957, 2.168975, 3.883867, 3.012901, 3.564255, 1.809821, 2.469867, 3.46857, 3.427226, 3.730365, 2.293451, 3.283702, 3.291594, 2.346601, 2.729807, 3.973846, 3.026795, 3.175831, 2.664512, 3.138977, 3.345586, 3.411898, 4.072533, 1.826528, 3.074796, 2.328734, 3.276652, 3.794981, 2.70656, 2.083811, 3.44407, 3.796744, 3.258427, 2.352164, 3.027308, 2.607675, 2.475324, 4.165256, 3.701353, 3.4713, 3.413129, 2.59423, 3.238124, 3.510629, 3.322692, 3.521572, 2.847815, 4.238555, 3.48561, 3.93355, 3.336021, 2.846023, 3.268262, 3.412435, 2.518049, 2.572459, 3.943473, 2.80409, 2.509684, 3.343666, 2.747478, 4.07886, 2.700101, 2.652727 };

            // Fit the distribution from sample data (matching BCaQuantileCI which re-estimates)
            var dist = new Normal();
            ((IEstimation)dist).Estimate(sampleData, ParameterEstimationMethod.MethodOfMoments);

            var parms = new ParameterSet(dist.GetParameters, double.NaN);
            var boot = new Bootstrap<double[]>(sampleData, parms);
            boot.Replicates = 10000;
            boot.PRNGSeed = 12345;

            boot.ResampleFunction = (data, ps, rng) =>
            {
                var d = new Normal(ps.Values[0], ps.Values[1]);
                return d.GenerateRandomValues(sampleData.Length, rng.Next());
            };

            boot.FitFunction = (sample) =>
            {
                var d = new Normal();
                ((IEstimation)d).Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
                if (!d.ParametersValid) throw new Exception("Invalid parameters.");
                return new ParameterSet(d.GetParameters, double.NaN);
            };

            boot.StatisticFunction = (ps) =>
            {
                var d = new Normal(ps.Values[0], ps.Values[1]);
                var result = new double[_probabilities.Length];
                for (int i = 0; i < _probabilities.Length; i++)
                    result[i] = d.InverseCDF(_probabilities[i]);
                return result;
            };

            // Set up jackknife delegates for BCa
            boot.JackknifeFunction = (data, idx) =>
            {
                var list = new List<double>(data);
                list.RemoveAt(idx);
                return list.ToArray();
            };
            boot.SampleSizeFunction = (data) => data.Length;

            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.BCa);

            // Compare against the true Monte Carlo CIs
            var trueCIs = dist.MonteCarloConfidenceIntervals(sampleData.Length, 10000, _probabilities, new double[] { 0.05, 0.95 });
            for (int i = 0; i < _probabilities.Length; i++)
            {
                Assert.AreEqual(trueCIs[i, 0], results.StatisticResults[i].LowerCI, 0.01 * Math.Abs(trueCIs[i, 0]),
                    $"Lower CI mismatch at p={_probabilities[i]}");
                Assert.AreEqual(trueCIs[i, 1], results.StatisticResults[i].UpperCI, 0.01 * Math.Abs(trueCIs[i, 1]),
                    $"Upper CI mismatch at p={_probabilities[i]}");
            }
        }

        /// <summary>
        /// Test that parameter-level confidence intervals are computed correctly.
        /// </summary>
        [TestMethod]
        public void Test_ParameterCIs()
        {
            var boot = CreateNormalBootstrap();
            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile);

            // Should have 2 parameter results (mu, sigma)
            Assert.HasCount(2, results.ParameterResults);

            // Population estimates should match
            Assert.AreEqual(_mu, results.ParameterResults[0].PopulationEstimate, 1e-10);
            Assert.AreEqual(_sigma, results.ParameterResults[1].PopulationEstimate, 1e-10);

            // CIs should bracket the population values
            Assert.IsLessThan(_mu, results.ParameterResults[0].LowerCI);
            Assert.IsGreaterThan(_mu, results.ParameterResults[0].UpperCI);
            Assert.IsLessThan(_sigma, results.ParameterResults[1].LowerCI);
            Assert.IsGreaterThan(_sigma, results.ParameterResults[1].UpperCI);

            // Valid count should be high
            Assert.IsGreaterThan(0.95 * boot.Replicates, results.ParameterResults[0].ValidCount);
        }

        /// <summary>
        /// Test that the bootstrap handles some failed replicates and tracks valid count correctly.
        /// </summary>
        [TestMethod]
        public void Test_ErrorHandlingAndValidCount()
        {
            var parms = new ParameterSet(new double[] { _mu, _sigma }, double.NaN);
            var boot = new Bootstrap<double[]>(null, parms);
            boot.Replicates = 1000;
            boot.PRNGSeed = 12345;
            boot.MaxRetries = 2;

            int callCount = 0;

            boot.ResampleFunction = (data, ps, rng) =>
            {
                var d = new Normal(ps.Values[0], ps.Values[1]);
                return d.GenerateRandomValues(_sampleSize, rng.Next());
            };

            // Fit function that fails ~10% of the time
            boot.FitFunction = (sample) =>
            {
                int count = Interlocked.Increment(ref callCount);
                if (count % 10 == 0) throw new Exception("Simulated failure");
                var d = new Normal();
                ((IEstimation)d).Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
                return new ParameterSet(d.GetParameters, double.NaN);
            };

            boot.StatisticFunction = (ps) =>
            {
                var d = new Normal(ps.Values[0], ps.Values[1]);
                return new double[] { d.InverseCDF(0.99) };
            };

            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile);

            // With retries, most replicates should succeed
            Assert.IsGreaterThan(0.8 * boot.Replicates, results.StatisticResults[0].ValidCount,
                $"Valid count {results.StatisticResults[0].ValidCount} too low");

            // CIs should still be reasonable (non-NaN)
            Assert.IsFalse(double.IsNaN(results.StatisticResults[0].LowerCI));
            Assert.IsFalse(double.IsNaN(results.StatisticResults[0].UpperCI));
        }

        /// <summary>
        /// Test that requesting BCa without JackknifeFunction throws.
        /// </summary>
        [TestMethod]
        public void Test_BCa_RequiresJackknifeFunction()
        {
            var boot = CreateNormalBootstrap();
            boot.Run();
            Assert.Throws<InvalidOperationException>(() => boot.GetConfidenceIntervals(BootstrapCIMethod.BCa));
        }

        /// <summary>
        /// Test that requesting Bootstrap-t without RunWithStudentizedBootstrap throws.
        /// </summary>
        [TestMethod]
        public void Test_BootstrapT_RequiresStudentizedRun()
        {
            var boot = CreateNormalBootstrap();
            boot.Run();
            Assert.Throws<InvalidOperationException>(() => boot.GetConfidenceIntervals(BootstrapCIMethod.BootstrapT));
        }

        /// <summary>
        /// Test that calling GetConfidenceIntervals before Run() throws.
        /// </summary>
        [TestMethod]
        public void Test_GetCI_RequiresRun()
        {
            var boot = CreateNormalBootstrap();
            Assert.Throws<InvalidOperationException>(() => boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile));
        }

        /// <summary>
        /// Test that the double bootstrap method runs and produces reasonable results.
        /// </summary>
        [TestMethod]
        public void Test_DoubleBootstrap()
        {
            var boot = CreateNormalBootstrap();
            boot.Replicates = 1000;
            boot.RunDoubleBootstrap(100);
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile);

            // CIs should bracket the population quantiles
            var dist = new Normal(_mu, _sigma);
            for (int i = 0; i < _probabilities.Length; i++)
            {
                double q = dist.InverseCDF(_probabilities[i]);
                Assert.IsLessThan(q, results.StatisticResults[i].LowerCI,
                    $"Lower CI {results.StatisticResults[i].LowerCI} not below quantile {q} at p={_probabilities[i]}");
                Assert.IsGreaterThan(q, results.StatisticResults[i].UpperCI,
                    $"Upper CI {results.StatisticResults[i].UpperCI} not above quantile {q} at p={_probabilities[i]}");
            }
        }

        /// <summary>
        /// Test that standard error and mean are computed correctly for percentile method.
        /// </summary>
        [TestMethod]
        public void Test_StandardErrorAndMean()
        {
            var boot = CreateNormalBootstrap();
            boot.Replicates = 5000;
            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile);

            // For the median (p=0.5), the SE should be approximately sigma/sqrt(n) ~ 0.0557
            int medianIdx = Array.IndexOf(_probabilities, 0.5);
            Assert.IsGreaterThan(0.0, results.StatisticResults[medianIdx].StandardError);
            Assert.IsLessThan(0.15, results.StatisticResults[medianIdx].StandardError,
                $"SE for median = {results.StatisticResults[medianIdx].StandardError} seems too large");

            // Mean should be close to the population median
            Assert.AreEqual(_mu, results.StatisticResults[medianIdx].Mean, 0.05,
                $"Mean for median = {results.StatisticResults[medianIdx].Mean} too far from {_mu}");
        }

        /// <summary>
        /// Test that FailedReplicates count is zero when no failures occur.
        /// </summary>
        [TestMethod]
        public void Test_NoFailures()
        {
            var boot = CreateNormalBootstrap();
            boot.Replicates = 500;
            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile);

            Assert.AreEqual(0, results.FailedReplicates);
            for (int i = 0; i < _probabilities.Length; i++)
            {
                Assert.AreEqual(500, results.StatisticResults[i].ValidCount);
                Assert.AreEqual(500, results.StatisticResults[i].TotalCount);
            }
        }

    }
}
