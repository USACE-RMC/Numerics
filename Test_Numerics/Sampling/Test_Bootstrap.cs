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
            // resampling is driven by the parameters, so null is passed.
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

            // Verify there are results for all probabilities
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
        /// FailedAttempts counts every attempt that threw: with a fit that throws on every tenth call, it
        /// equals the number of calls divided by ten, and a run without failures reports zero.
        /// </summary>
        [TestMethod]
        public void Test_FailedAttempts_CountsEveryThrownAttempt()
        {
            int callCount = 0;
            var boot = CreateNormalBootstrap();
            boot.Replicates = 500;
            boot.MaxRetries = 3;
            var fit = boot.FitFunction!;
            boot.FitFunction = sample =>
            {
                if (Interlocked.Increment(ref callCount) % 10 == 0) throw new InvalidOperationException("Simulated failure");
                return fit(sample);
            };
            boot.Run();

            // Every call either throws or yields a replicate, so calls = successful replicates + thrown attempts.
            Assert.AreEqual(callCount / 10, boot.FailedAttempts);
            Assert.AreEqual(boot.Replicates - boot.FailedReplicates + boot.FailedAttempts, callCount);

            boot.FitFunction = fit;
            boot.Run();
            Assert.AreEqual(0, boot.FailedAttempts);
        }

        /// <summary>
        /// When every fit throws, each replicate contributes one failed attempt per retry.
        /// </summary>
        [TestMethod]
        public void Test_FailedAttempts_AllFail_EqualsReplicatesTimesRetries()
        {
            var boot = CreateNormalBootstrap();
            boot.Replicates = 100;
            boot.MaxRetries = 4;
            boot.FitFunction = sample => throw new InvalidOperationException("Simulated failure");
            boot.Run();

            Assert.AreEqual(100, boot.FailedReplicates);
            Assert.AreEqual(400, boot.FailedAttempts);
        }

        /// <summary>
        /// A retry continues the replicate's own generator: with one replicate whose first two fits fail, the
        /// three attempts together consume exactly the first draws of the replicate's MersenneTwister, as an
        /// independent reconstruction from PRNGSeed gives.
        /// </summary>
        [TestMethod]
        public void Test_Retries_ContinueTheReplicateStream()
        {
            var draws = new List<double>();
            int fits = 0;
            var boot = new Bootstrap<double[]>(null, new ParameterSet(new double[] { 0d }, double.NaN))
            {
                Replicates = 1,
                PRNGSeed = 2468,
                MaxRetries = 5,
                ResampleFunction = (data, ps, rng) =>
                {
                    var sample = new[] { rng.NextDouble(), rng.NextDouble(), rng.NextDouble() };
                    draws.AddRange(sample);
                    return sample;
                },
                FitFunction = sample =>
                {
                    if (++fits <= 2) throw new InvalidOperationException("Simulated failure");
                    return new ParameterSet(new[] { sample.Average() }, double.NaN);
                },
                StatisticFunction = ps => new[] { ps.Values[0] }
            };
            boot.Run();

            var reference = new MersenneTwister(new MersenneTwister(2468).NextIntegers(1)[0]);
            Assert.HasCount(9, draws);
            for (int i = 0; i < draws.Count; i++)
                Assert.AreEqual(reference.NextDouble(), draws[i], 0d);
            Assert.AreEqual(2, boot.FailedAttempts);
            Assert.AreEqual((draws[6] + draws[7] + draws[8]) / 3d, boot.BootstrapStatistics[0, 0], 1E-15);
        }

        /// <summary>
        /// Records every reported progress value synchronously, so a test can inspect them after a run.
        /// </summary>
        private sealed class RecordingProgress : IProgress<int>
        {
            private int _calls;
            private int _maximum;

            /// <summary>Gets the number of reports.</summary>
            public int Calls => _calls;

            /// <summary>Gets the largest reported value.</summary>
            public int Maximum => _maximum;

            /// <inheritdoc/>
            public void Report(int value)
            {
                Interlocked.Increment(ref _calls);
                int current;
                while (value > (current = Volatile.Read(ref _maximum)))
                    Interlocked.CompareExchange(ref _maximum, value, current);
            }
        }

        /// <summary>
        /// Every outer replicate is reported once, the last count is the replicate total, and the hooks
        /// change no bootstrap statistic.
        /// </summary>
        [TestMethod]
        public void Test_ProgressAndToken_ReportEveryReplicate_AndDoNotChangeResults()
        {
            var plain = CreateNormalBootstrap();
            plain.Replicates = 500;
            plain.Run();

            var progress = new RecordingProgress();
            using var source = new CancellationTokenSource();
            var hooked = CreateNormalBootstrap();
            hooked.Replicates = 500;
            hooked.Progress = progress;
            hooked.CancellationToken = source.Token;
            hooked.Run();

            Assert.AreEqual(500, progress.Calls);
            Assert.AreEqual(500, progress.Maximum);
            CollectionAssert.AreEqual(plain.BootstrapStatistics, hooked.BootstrapStatistics);
        }

        /// <summary>
        /// The double bootstrap reports outer replicates only.
        /// </summary>
        [TestMethod]
        public void Test_Progress_DoubleBootstrap_ReportsOuterReplicates()
        {
            var progress = new RecordingProgress();
            var boot = CreateNormalBootstrap();
            boot.Replicates = 50;
            boot.Progress = progress;
            boot.RunDoubleBootstrap(20);

            Assert.AreEqual(50, progress.Calls);
            Assert.AreEqual(50, progress.Maximum);
        }

        /// <summary>
        /// A cancelled token stops the run with <see cref="OperationCanceledException"/> and leaves no
        /// results to request intervals from.
        /// </summary>
        [TestMethod]
        public void Test_CancellationToken_Cancelled_ThrowsAndLeavesNoResults()
        {
            using var source = new CancellationTokenSource();
            source.Cancel();
            var boot = CreateNormalBootstrap();
            boot.Replicates = 500;
            boot.CancellationToken = source.Token;

            Assert.Throws<OperationCanceledException>(() => boot.Run());
            Assert.Throws<InvalidOperationException>(() => boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile));
        }

        /// <summary>
        /// A token cancelled part way through stops the run before every replicate is attempted.
        /// </summary>
        [TestMethod]
        public void Test_CancellationToken_CancelledDuringRun_StopsEarly()
        {
            using var source = new CancellationTokenSource();
            var progress = new CancellingProgress(source, 50);
            var boot = CreateNormalBootstrap();
            boot.Replicates = 5000;
            boot.Progress = progress;
            boot.CancellationToken = source.Token;

            Assert.Throws<OperationCanceledException>(() => boot.Run());
            Assert.IsLessThan(5000, progress.Calls);
        }

        /// <summary>
        /// A studentized run cancelled after a completed regular run leaves no results to request intervals from,
        /// rather than the partly overwritten results of the earlier run.
        /// </summary>
        [TestMethod]
        public void Test_CancellationToken_CancelledStudentizedAfterRun_LeavesNoResults()
        {
            var boot = CreateNormalBootstrap();
            boot.Replicates = 500;
            boot.InnerReplicates = 20;
            boot.UseParallel = false;
            boot.Run();
            Assert.IsNotNull(boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile));

            using var source = new CancellationTokenSource();
            boot.Progress = new CancellingProgress(source, 50);
            boot.CancellationToken = source.Token;

            Assert.Throws<OperationCanceledException>(() => boot.RunWithStudentizedBootstrap());
            Assert.ThrowsExactly<InvalidOperationException>(() => boot.GetConfidenceIntervals(BootstrapCIMethod.Percentile));
        }

        /// <summary>
        /// Cancels a token once a set number of replicates have been reported.
        /// </summary>
        private sealed class CancellingProgress : IProgress<int>
        {
            private readonly CancellationTokenSource _source;
            private readonly int _threshold;
            private int _calls;

            /// <summary>
            /// Creates the progress receiver.
            /// </summary>
            /// <param name="source">The source to cancel.</param>
            /// <param name="threshold">The report count at which to cancel.</param>
            public CancellingProgress(CancellationTokenSource source, int threshold)
            {
                _source = source;
                _threshold = threshold;
            }

            /// <summary>Gets the number of reports.</summary>
            public int Calls => _calls;

            /// <inheritdoc/>
            public void Report(int value)
            {
                if (Interlocked.Increment(ref _calls) == _threshold) _source.Cancel();
            }
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

        /// <summary>
        /// Creates a deterministic sample used by the BCa jackknife tests.
        /// </summary>
        /// <returns>A fixed Normal sample.</returns>
        private static double[] CreateJackknifeSample()
        {
            return new Normal(3.1d, 0.55d).GenerateRandomValues(40, 8675309);
        }

        /// <summary>
        /// Creates a small BCa-ready bootstrap over the supplied sample. The jackknife delegate is left
        /// unset so each test can supply its own leave-one-out behavior.
        /// </summary>
        /// <param name="sampleData">The observed sample.</param>
        /// <param name="replicates">The number of bootstrap replicates.</param>
        /// <returns>The configured bootstrap.</returns>
        private static Bootstrap<double[]> CreateJackknifeBootstrap(double[] sampleData, int replicates)
        {
            var dist = new Normal();
            ((IEstimation)dist).Estimate(sampleData, ParameterEstimationMethod.MethodOfMoments);

            var boot = new Bootstrap<double[]>(sampleData, new ParameterSet(dist.GetParameters, double.NaN))
            {
                Replicates = replicates,
                PRNGSeed = 12345
            };

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
                return new double[] { d.InverseCDF(0.5d), d.InverseCDF(0.99d) };
            };

            boot.SampleSizeFunction = (data) => data.Length;
            return boot;
        }

        /// <summary>
        /// Reproduces the BCa limit that a zero acceleration constant would produce for one statistic,
        /// from the same bootstrap ensemble the analysis used.
        /// </summary>
        /// <param name="values">The bootstrap values for one statistic.</param>
        /// <param name="populationEstimate">The original statistic estimate.</param>
        /// <param name="alpha">The two-sided alpha level.</param>
        /// <param name="upper">True for the upper limit, false for the lower limit.</param>
        /// <returns>The limit implied by an acceleration constant of exactly zero.</returns>
        /// <remarks>
        /// With a = 0 the BCa denominator 1 - a(Z₀ + z) is exactly one, so the division leaves its
        /// numerator bit-unchanged and the adjusted probability is Φ(Z₀ + (Z₀ + z)). Every step below
        /// therefore reproduces the library's arithmetic exactly, which lets the caller compare bit
        /// patterns instead of choosing a tolerance.
        /// </remarks>
        private static double ZeroAccelerationBCaLimit(double[] values, double populationEstimate, double alpha, bool upper)
        {
            var validValues = values.Where(Tools.IsFinite).ToArray();
            int countLeq = 0;
            for (int i = 0; i < validValues.Length; i++)
                if (validValues[i] <= populationEstimate) countLeq++;

            double p0 = (double)(countLeq + 1) / (validValues.Length + 1);
            Array.Sort(validValues);

            double z0 = Normal.StandardZ(p0);
            double z = Normal.StandardZ(upper ? 1d - alpha / 2d : alpha / 2d);
            double adjusted = Normal.StandardCDF(z0 + (z0 + z) / 1d);
            return Statistics.Percentile(validValues, adjusted, true);
        }

        /// <summary>
        /// Test that a wholly failed jackknife reports the cause instead of producing NaN acceleration.
        /// </summary>
        /// <remarks>
        /// A NaN acceleration constant reaches Statistics.Percentile as the percentile argument, where it
        /// silently returned NaN bounds on .NET Core and raised an index error on .NET Framework.
        /// </remarks>
        [TestMethod]
        public void Test_BCa_AllJackknifeReplicatesFail()
        {
            var sampleData = CreateJackknifeSample();
            var boot = CreateJackknifeBootstrap(sampleData, 500);
            boot.JackknifeFunction = (data, idx) => throw new InvalidOperationException("Jackknife replicate rejected.");

            boot.Run();

            var exception = Assert.Throws<InvalidOperationException>(() => boot.GetConfidenceIntervals(BootstrapCIMethod.BCa));
            StringAssert.Contains(exception.Message, nameof(Bootstrap<double[]>.JackknifeFunction));
            StringAssert.Contains(exception.Message, "Jackknife replicate rejected.");
            Assert.AreEqual(sampleData.Length, boot.FailedJackknifeReplicates);
        }

        /// <summary>
        /// Test that a partially failing jackknife still produces finite BCa bounds and reports the failures.
        /// </summary>
        [TestMethod]
        public void Test_BCa_PartialJackknifeFailuresAreReported()
        {
            var sampleData = CreateJackknifeSample();
            var boot = CreateJackknifeBootstrap(sampleData, 500);
            boot.JackknifeFunction = (data, idx) =>
            {
                if (idx % 5 == 0) throw new InvalidOperationException("Jackknife replicate rejected.");
                var list = new List<double>(data);
                list.RemoveAt(idx);
                return list.ToArray();
            };

            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.BCa);

            Assert.AreEqual(sampleData.Length / 5, boot.FailedJackknifeReplicates);
            for (int i = 0; i < results.StatisticResults.Length; i++)
            {
                Assert.IsTrue(Tools.IsFinite(results.StatisticResults[i].LowerCI), $"Lower CI is not finite for statistic {i}.");
                Assert.IsTrue(Tools.IsFinite(results.StatisticResults[i].UpperCI), $"Upper CI is not finite for statistic {i}.");
                Assert.IsLessThan(results.StatisticResults[i].UpperCI, results.StatisticResults[i].LowerCI);
            }
        }

        /// <summary>
        /// Test that a statistic with no jackknife variation yields a zero acceleration constant, so BCa
        /// degenerates to the bias-corrected interval instead of dividing zero by zero.
        /// </summary>
        [TestMethod]
        public void Test_BCa_DegenerateJackknifeYieldsBiasCorrectedInterval()
        {
            var sampleData = CreateJackknifeSample();
            var boot = CreateJackknifeBootstrap(sampleData, 2000);

            // Returning the whole sample for every index makes each jackknife statistic identical to the
            // population estimate, so the jackknife second moment is exactly zero.
            boot.JackknifeFunction = (data, idx) => (double[])data.Clone();

            boot.Run();
            var bca = boot.GetConfidenceIntervals(BootstrapCIMethod.BCa);
            var biasCorrected = boot.GetConfidenceIntervals(BootstrapCIMethod.BiasCorrected);

            Assert.AreEqual(0, boot.FailedJackknifeReplicates);
            for (int i = 0; i < bca.StatisticResults.Length; i++)
            {
                Assert.IsTrue(Tools.IsFinite(bca.StatisticResults[i].LowerCI), $"Lower CI is not finite for statistic {i}.");
                Assert.IsTrue(Tools.IsFinite(bca.StatisticResults[i].UpperCI), $"Upper CI is not finite for statistic {i}.");

                // Pin the acceleration to exactly zero. Reproducing the a = 0 limit from the bootstrap
                // ensemble is bit-exact, so any non-zero acceleration moves the limit and fails here.
                var values = boot.BootstrapStatistics.GetColumn(i);
                double populationEstimate = bca.StatisticResults[i].PopulationEstimate;
                Assert.AreEqual(
                    BitConverter.DoubleToInt64Bits(ZeroAccelerationBCaLimit(values, populationEstimate, 0.1d, false)),
                    BitConverter.DoubleToInt64Bits(bca.StatisticResults[i].LowerCI),
                    $"The lower CI does not match a zero acceleration for statistic {i}.");
                Assert.AreEqual(
                    BitConverter.DoubleToInt64Bits(ZeroAccelerationBCaLimit(values, populationEstimate, 0.1d, true)),
                    BitConverter.DoubleToInt64Bits(bca.StatisticResults[i].UpperCI),
                    $"The upper CI does not match a zero acceleration for statistic {i}.");

                // With a zero acceleration the BCa limits also reduce to the bias-corrected limits. The
                // residual gap is not zero because the two methods use different plotting positions for
                // the bias proportion (countLeq / (B + 1) versus (countLeq + 1) / (B + 1)), an offset of
                // 1/(B+1) in P0 that propagates to roughly 2E-4 in the limit here. The tolerance bounds
                // that offset and nothing wider.
                Assert.AreEqual(biasCorrected.StatisticResults[i].LowerCI, bca.StatisticResults[i].LowerCI,
                    1E-3 * Math.Abs(biasCorrected.StatisticResults[i].LowerCI));
                Assert.AreEqual(biasCorrected.StatisticResults[i].UpperCI, bca.StatisticResults[i].UpperCI,
                    1E-3 * Math.Abs(biasCorrected.StatisticResults[i].UpperCI));
            }
        }

        /// <summary>
        /// Test that a leave-one-out replicate whose statistic is non-finite is counted as a failure
        /// rather than silently poisoning the jackknife moments.
        /// </summary>
        /// <remarks>
        /// A non-finite statistic drives the second moment to NaN. Because NaN fails every ordered
        /// comparison, it would otherwise fall through the zero-acceleration fallback and quietly
        /// degenerate BCa to the bias-corrected interval while reporting no failures at all.
        /// </remarks>
        [TestMethod]
        public void Test_BCa_NonFiniteJackknifeStatisticIsCountedAsFailure()
        {
            var sampleData = CreateJackknifeSample();
            var boot = CreateJackknifeBootstrap(sampleData, 1000);

            // One leave-one-out index returns a sample two elements short. Bootstrap resamples have the
            // full length and every other leave-one-out sample has length n-1, so carrying the length
            // through the parameter set's fitness field marks exactly one replicate.
            int markerLength = sampleData.Length - 2;
            boot.JackknifeFunction = (data, idx) =>
            {
                var list = new List<double>(data);
                list.RemoveAt(idx);
                if (idx == 7) list.RemoveAt(0);
                return list.ToArray();
            };

            boot.FitFunction = (sample) =>
            {
                var d = new Normal();
                ((IEstimation)d).Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
                if (!d.ParametersValid) throw new Exception("Invalid parameters.");
                return new ParameterSet(d.GetParameters, sample.Length);
            };

            boot.StatisticFunction = (ps) =>
            {
                if (ps.Fitness == markerLength) return new double[] { double.NaN, double.NaN };
                var d = new Normal(ps.Values[0], ps.Values[1]);
                return new double[] { d.InverseCDF(0.5d), d.InverseCDF(0.99d) };
            };

            boot.Run();
            var results = boot.GetConfidenceIntervals(BootstrapCIMethod.BCa);

            Assert.AreEqual(1, boot.FailedJackknifeReplicates);
            for (int i = 0; i < results.StatisticResults.Length; i++)
            {
                Assert.IsTrue(Tools.IsFinite(results.StatisticResults[i].LowerCI), $"Lower CI is not finite for statistic {i}.");
                Assert.IsTrue(Tools.IsFinite(results.StatisticResults[i].UpperCI), $"Upper CI is not finite for statistic {i}.");

                // The surviving 39 replicates still carry jackknife variation, so the acceleration must
                // not have collapsed to the zero fallback.
                var values = boot.BootstrapStatistics.GetColumn(i);
                double populationEstimate = results.StatisticResults[i].PopulationEstimate;
                Assert.AreNotEqual(
                    BitConverter.DoubleToInt64Bits(ZeroAccelerationBCaLimit(values, populationEstimate, 0.1d, false)),
                    BitConverter.DoubleToInt64Bits(results.StatisticResults[i].LowerCI),
                    $"The lower CI collapsed to a zero acceleration for statistic {i}.");
            }
        }

        /// <summary>
        /// Test that two BCa runs over the same fixture produce bit-identical bounds. The jackknife
        /// reduction must not depend on the thread scheduler.
        /// </summary>
        [TestMethod]
        public void Test_BCa_BoundsAreBitReproducible()
        {
            var sampleData = CreateJackknifeSample();

            var first = RunReproducibleBCa(sampleData);
            var second = RunReproducibleBCa(sampleData);

            Assert.HasCount(first.Length, second);
            for (int i = 0; i < first.Length; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(first[i].LowerCI), BitConverter.DoubleToInt64Bits(second[i].LowerCI),
                    $"Lower CI is not bit-identical for statistic {i}.");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(first[i].UpperCI), BitConverter.DoubleToInt64Bits(second[i].UpperCI),
                    $"Upper CI is not bit-identical for statistic {i}.");
            }
        }

        /// <summary>
        /// Runs a complete BCa analysis over a fresh bootstrap instance.
        /// </summary>
        /// <param name="sampleData">The observed sample.</param>
        /// <returns>The BCa statistic results.</returns>
        private static BootstrapStatisticResult[] RunReproducibleBCa(double[] sampleData)
        {
            var boot = CreateJackknifeBootstrap(sampleData, 1000);
            boot.JackknifeFunction = (data, idx) =>
            {
                var list = new List<double>(data);
                list.RemoveAt(idx);
                return list.ToArray();
            };

            boot.Run();
            return boot.GetConfidenceIntervals(BootstrapCIMethod.BCa).StatisticResults;
        }

        /// <summary>
        /// Asserts two bootstrap runs stored bit-identical statistics, parameter sets and counts.
        /// </summary>
        private static void AssertSameRun(Bootstrap<double[]> expected, Bootstrap<double[]> actual)
        {
            Assert.AreEqual(expected.FailedReplicates, actual.FailedReplicates);
            Assert.AreEqual(expected.FailedAttempts, actual.FailedAttempts);
            var a = expected.BootstrapStatistics;
            var b = actual.BootstrapStatistics;
            Assert.AreEqual(a.GetLength(0), b.GetLength(0));
            for (int i = 0; i < a.GetLength(0); i++)
            {
                for (int j = 0; j < a.GetLength(1); j++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(a[i, j]), BitConverter.DoubleToInt64Bits(b[i, j]));
                var pa = expected.BootstrapParameterSets[i].Values;
                var pb = actual.BootstrapParameterSets[i].Values;
                for (int j = 0; j < pa.Length; j++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(pa[j]), BitConverter.DoubleToInt64Bits(pb[j]));
            }
        }

        /// <summary>
        /// UseParallel defaults to true; serial and parallel runs of every run method give bit-identical results.
        /// </summary>
        [TestMethod]
        public void Test_UseParallel_SerialMatchesParallelBits()
        {
            Assert.IsTrue(CreateNormalBootstrap().UseParallel);

            var parallel = CreateNormalBootstrap();
            parallel.Replicates = 400;
            parallel.Run();
            var serial = CreateNormalBootstrap();
            serial.Replicates = 400;
            serial.UseParallel = false;
            serial.Run();
            AssertSameRun(parallel, serial);
            var pci = parallel.GetConfidenceIntervals(BootstrapCIMethod.Percentile).StatisticResults;
            var sci = serial.GetConfidenceIntervals(BootstrapCIMethod.Percentile).StatisticResults;
            for (int i = 0; i < pci.Length; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(pci[i].LowerCI), BitConverter.DoubleToInt64Bits(sci[i].LowerCI));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(pci[i].UpperCI), BitConverter.DoubleToInt64Bits(sci[i].UpperCI));
            }

            parallel = CreateNormalBootstrap();
            parallel.Replicates = 60;
            parallel.RunDoubleBootstrap(20);
            serial = CreateNormalBootstrap();
            serial.Replicates = 60;
            serial.UseParallel = false;
            serial.RunDoubleBootstrap(20);
            AssertSameRun(parallel, serial);

            parallel = CreateNormalBootstrap();
            parallel.Replicates = 60;
            parallel.InnerReplicates = 20;
            parallel.RunWithStudentizedBootstrap();
            serial = CreateNormalBootstrap();
            serial.Replicates = 60;
            serial.InnerReplicates = 20;
            serial.UseParallel = false;
            serial.RunWithStudentizedBootstrap();
            AssertSameRun(parallel, serial);
            var pt = parallel.GetConfidenceIntervals(BootstrapCIMethod.BootstrapT).StatisticResults;
            var st = serial.GetConfidenceIntervals(BootstrapCIMethod.BootstrapT).StatisticResults;
            for (int i = 0; i < pt.Length; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(pt[i].LowerCI), BitConverter.DoubleToInt64Bits(st[i].LowerCI));
        }

        /// <summary>
        /// A serial run calls the delegates on the calling thread, in replicate order, and reports progress in order.
        /// </summary>
        [TestMethod]
        public void Test_UseParallel_False_RunsInOrderOnCallingThread()
        {
            int caller = Environment.CurrentManagedThreadId;
            var threads = new List<int>();
            var boot = CreateNormalBootstrap();
            boot.Replicates = 50;
            boot.UseParallel = false;
            var fit = boot.FitFunction!;
            boot.FitFunction = sample => { threads.Add(Environment.CurrentManagedThreadId); return fit(sample); };
            var progress = new OrderedProgress();
            boot.Progress = progress;
            boot.Run();

            Assert.HasCount(50, threads);
            Assert.IsTrue(threads.All(id => id == caller));
            CollectionAssert.AreEqual(Enumerable.Range(1, 50).ToArray(), progress.Values.ToArray());
        }

        /// <summary>
        /// A serial run honours cancellation and wraps an escaping exception as the parallel run does.
        /// </summary>
        [TestMethod]
        public void Test_UseParallel_False_CancellationAndExceptions()
        {
            using var cancelled = new CancellationTokenSource();
            cancelled.Cancel();
            var boot = CreateNormalBootstrap();
            boot.Replicates = 20;
            boot.UseParallel = false;
            boot.CancellationToken = cancelled.Token;
            Assert.Throws<OperationCanceledException>(() => boot.Run());

            using var source = new CancellationTokenSource();
            var stopping = new CancellingProgress(source, 5);
            boot = CreateNormalBootstrap();
            boot.Replicates = 1000;
            boot.UseParallel = false;
            boot.Progress = stopping;
            boot.CancellationToken = source.Token;
            Assert.Throws<OperationCanceledException>(() => boot.Run());
            Assert.AreEqual(5, stopping.Calls);

            foreach (bool useParallel in new[] { true, false })
            {
                boot = CreateNormalBootstrap();
                boot.Replicates = 10;
                boot.UseParallel = useParallel;
                boot.Progress = new ThrowingProgress();
                var error = Assert.Throws<AggregateException>(() => boot.Run());
                Assert.IsInstanceOfType<InvalidOperationException>(error.InnerExceptions[0]);
            }
        }

        /// <summary>
        /// With UseParallel false the BCa leave-one-out jackknife runs in index order on the calling thread, and its
        /// intervals and failure count are bit-identical to the parallel jackknife.
        /// </summary>
        [TestMethod]
        public void Test_UseParallel_False_SerialisesBCaJackknife()
        {
            // 200 points spread over the 64 accumulation chunks, with some failing replicates.
            var sampleData = new Normal(3.1d, 0.55d).GenerateRandomValues(200, 8675309);
            BootstrapResults parallelResults = null;
            int parallelFailures = -1;
            foreach (bool useParallel in new[] { true, false })
            {
                int caller = Environment.CurrentManagedThreadId;
                var calls = new List<(int Index, int Thread)>();
                var boot = CreateJackknifeBootstrap(sampleData, 300);
                boot.UseParallel = useParallel;
                boot.JackknifeFunction = (data, idx) =>
                {
                    lock (calls) calls.Add((idx, Environment.CurrentManagedThreadId));
                    if (idx % 17 == 3) throw new InvalidOperationException("Jackknife replicate rejected.");
                    var list = new List<double>(data);
                    list.RemoveAt(idx);
                    return list.ToArray();
                };
                boot.Run();
                var results = boot.GetConfidenceIntervals(BootstrapCIMethod.BCa);
                Assert.AreEqual(12, boot.FailedJackknifeReplicates);
                Assert.HasCount(sampleData.Length, calls);

                if (useParallel)
                {
                    parallelResults = results;
                    parallelFailures = boot.FailedJackknifeReplicates;
                    continue;
                }

                CollectionAssert.AreEqual(Enumerable.Range(0, sampleData.Length).ToArray(), calls.Select(c => c.Index).ToArray());
                Assert.IsTrue(calls.All(c => c.Thread == caller));
                Assert.AreEqual(parallelFailures, boot.FailedJackknifeReplicates);
                for (int i = 0; i < results.StatisticResults.Length; i++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(parallelResults.StatisticResults[i].LowerCI), BitConverter.DoubleToInt64Bits(results.StatisticResults[i].LowerCI));
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(parallelResults.StatisticResults[i].UpperCI), BitConverter.DoubleToInt64Bits(results.StatisticResults[i].UpperCI));
                }
            }
        }

        /// <summary>
        /// Records the reported values in arrival order; for serial runs only.
        /// </summary>
        private sealed class OrderedProgress : IProgress<int>
        {
            /// <summary>Gets the reported values in order.</summary>
            public List<int> Values { get; } = new List<int>();

            /// <inheritdoc/>
            public void Report(int value) => Values.Add(value);
        }

        /// <summary>
        /// A progress receiver that throws on its first report.
        /// </summary>
        private sealed class ThrowingProgress : IProgress<int>
        {
            /// <inheritdoc/>
            public void Report(int value) => throw new InvalidOperationException("Progress failed.");
        }

    }
}
