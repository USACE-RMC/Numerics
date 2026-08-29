using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading;

namespace Distributions.Univariate
{

    /// <summary>
    /// Unit test for bootstrap analysis of univariate distributions. 
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
    public class Test_BootstrapAnalysis
    {

        /// <summary>
        /// This test compares the quantile confidence intervals obtained by the 'Normal' method for the Normal distribution with those from the 'boot' package. 
        /// </summary>
        [TestMethod]
        public void Test_NormalCI()
        {
            var probabilities = new double[] { 0.999999, 0.999998, 0.999995, 0.99999, 0.99998, 0.99995, 0.9999, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);
            var CIs = boot.NormalQuantileCI(probabilities);

            /* Below are the results from 'boot' using comparable settings:
            *  Since bootstrap methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 1% of 'boot' results. 
            */

            var true05 = new double[] { 5.451061, 5.3807, 5.284459, 5.208962, 5.130887, 5.023248, 4.938038, 4.849109, 4.724954, 4.625182, 4.519401, 4.368283, 4.243235, 4.106139, 3.899265, 3.713675, 3.485485, 3.317521, 3.031074, 2.732176, 2.546064, 2.283231, 2.063337, 1.813848, 1.646703 };
            var true95 = new double[] { 6.098854, 6.010689, 5.890184, 5.795726, 5.698123, 5.563704, 5.457429, 5.34666, 5.192298, 5.068532, 4.937633, 4.751333, 4.597948, 4.430811, 4.181338, 3.961466, 3.698674, 3.512604, 3.213796, 2.927438, 2.759516, 2.531367, 2.345799, 2.138941, 2.001853 };

            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(true05[i], CIs[i, 0], 0.01 * true05[i]);
                Assert.AreEqual(true95[i], CIs[i, 1], 0.01 * true95[i]);
            }
        }

        /// <summary>
        /// This test compares the quantile confidence intervals obtained by the 'Percentile' method for the Normal distribution with those from the 'boot' package. 
        /// </summary>
        [TestMethod]
        public void Test_PercentileCI()
        {
            var probabilities = new double[] { 0.999999, 0.999998, 0.999995, 0.99999, 0.99998, 0.99995, 0.9999, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);
            var CIs = boot.PercentileQuantileCI(probabilities);

            /* Below are the results from 'boot' using comparable settings:
            *  Since bootstrap methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 1% of 'boot' results. 
            */

            var true05 = new double[] { 5.448065, 5.378149, 5.280887, 5.205002, 5.127415, 5.020843, 4.935679, 4.847712, 4.723488, 4.623268, 4.517452, 4.367232, 4.24338, 4.106147, 3.900035, 3.714043, 3.486811, 3.318758, 3.032698, 2.734326, 2.548027, 2.284871, 2.06533, 1.816327, 1.648123 };
            var true95 = new double[] { 6.095642, 6.00775, 5.887441, 5.792823, 5.696038, 5.562841, 5.457508, 5.346516, 5.191728, 5.067822, 4.9371, 4.750642, 4.596812, 4.428688, 4.178886, 3.95975, 3.69779, 3.51265, 3.215176, 2.930193, 2.761301, 2.533379, 2.3472, 2.139633, 2.003285 };

            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(true05[i], CIs[i, 0], 0.01 * true05[i]);
                Assert.AreEqual(true95[i], CIs[i, 1], 0.01 * true95[i]);
            }

        }

        /// <summary>
        /// This test compares the quantile confidence intervals obtained by the 'Bootstrap-t' or 'Student' method for the Normal distribution with the "true" Noncentral-t intervals. 
        /// </summary>
        [TestMethod]
        public void Test_BootstrapTCI()
        {
            var probabilities = new double[] { 0.999999, 0.999998, 0.999995, 0.99999, 0.99998, 0.99995, 0.9999, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);
            var CIs = boot.BootstrapTQuantileCI(probabilities);

            /* Since bootstrap methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 1% of 'boot' results. 
            */

            var trueCIs = dist.MonteCarloConfidenceIntervals(100, 10000, probabilities, new double[] { 0.05, 0.95 });
            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(trueCIs[i, 0], CIs[i, 0], 0.01 * trueCIs[i, 0]);
                Assert.AreEqual(trueCIs[i, 1], CIs[i, 1], 0.01 * trueCIs[i, 1]);
            }

        }


        /// <summary>
        /// This test compares the quantile confidence intervals obtained by the 'BCa' method for the Normal distribution with the "true" Noncentral-t intervals. 
        /// </summary>
        [TestMethod]
        public void Test_BCaCI()
        {
            var sampleData = new double[] { 3.292764, 3.354733, 2.945348, 2.773251, 3.302944, 2.091022, 3.315049, 2.861908, 2.85792, 2.540339, 2.941876, 3.908656, 3.185314, 3.260108, 2.624734, 3.40845, 2.556821, 2.834211, 3.560356, 3.149362, 3.389811, 3.727893, 2.677836, 2.223431, 2.201145, 3.902549, 2.759176, 3.31019, 3.306062, 2.918845, 3.405937, 4.098417, 4.024595, 3.816223, 3.127136, 3.245594, 2.837957, 2.168975, 3.883867, 3.012901, 3.564255, 1.809821, 2.469867, 3.46857, 3.427226, 3.730365, 2.293451, 3.283702, 3.291594, 2.346601, 2.729807, 3.973846, 3.026795, 3.175831, 2.664512, 3.138977, 3.345586, 3.411898, 4.072533, 1.826528, 3.074796, 2.328734, 3.276652, 3.794981, 2.70656, 2.083811, 3.44407, 3.796744, 3.258427, 2.352164, 3.027308, 2.607675, 2.475324, 4.165256, 3.701353, 3.4713, 3.413129, 2.59423, 3.238124, 3.510629, 3.322692, 3.521572, 2.847815, 4.238555, 3.48561, 3.93355, 3.336021, 2.846023, 3.268262, 3.412435, 2.518049, 2.572459, 3.943473, 2.80409, 2.509684, 3.343666, 2.747478, 4.07886, 2.700101, 2.652727 };
            var probabilities = new double[] { 0.999999, 0.999998, 0.999995, 0.99999, 0.99998, 0.99995, 0.9999, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);
            var CIs = boot.BCaQuantileCI(sampleData, probabilities);

            /* Since bootstrap methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 1% of 'boot' results. 
            */

            var trueCIs = dist.MonteCarloConfidenceIntervals(100, 10000, probabilities, new double[] { 0.05, 0.95 });
            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(trueCIs[i, 0], CIs[i, 0], 0.01 * trueCIs[i, 0]);
                Assert.AreEqual(trueCIs[i, 1], CIs[i, 1], 0.01 * trueCIs[i, 1]);
            }

        }

        /// <summary>
        /// This test verifies that UncertaintyAnalysisResults produces equivalent results to BootstrapAnalysis.Estimate.
        /// </summary>
        [TestMethod]
        public void Test_BootstrapAnalysis_UncertaintyAnalysisResults_Equivalence()
        {
            var probabilities = new double[] { 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            double alpha = 0.1;
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);

            // Generate bootstrap distributions once and share between both methods
            var bootDistributions = boot.Distributions();

            // Reference result from BootstrapAnalysis.Estimate
            var reference = boot.Estimate(probabilities, alpha, bootDistributions, recordParameterSets: false);

            // Result from UncertaintyAnalysisResults constructor
            var sampledDists = bootDistributions.Cast<UnivariateDistributionBase>().ToArray();
            var result = new UncertaintyAnalysisResults(dist, sampledDists, probabilities, alpha);

            // Compare ModeCurve
            Assert.HasCount(reference.ModeCurve.Length, result.ModeCurve);
            for (int i = 0; i < reference.ModeCurve.Length; i++)
            {
                Assert.AreEqual(reference.ModeCurve[i], result.ModeCurve[i], 1E-8,
                    $"ModeCurve mismatch at index {i}");
            }

            // Compare ConfidenceIntervals
            Assert.AreEqual(reference.ConfidenceIntervals.GetLength(0), result.ConfidenceIntervals.GetLength(0));
            Assert.AreEqual(reference.ConfidenceIntervals.GetLength(1), result.ConfidenceIntervals.GetLength(1));
            for (int i = 0; i < reference.ConfidenceIntervals.GetLength(0); i++)
            {
                for (int j = 0; j < reference.ConfidenceIntervals.GetLength(1); j++)
                {
                    Assert.AreEqual(reference.ConfidenceIntervals[i, j], result.ConfidenceIntervals[i, j], 1E-8,
                        $"ConfidenceIntervals mismatch at [{i},{j}]");
                }
            }

            // Compare MeanCurve
            Assert.HasCount(reference.MeanCurve.Length, result.MeanCurve);
            for (int i = 0; i < reference.MeanCurve.Length; i++)
            {
                Assert.AreEqual(reference.MeanCurve[i], result.MeanCurve[i], 1E-8,
                    $"MeanCurve mismatch at index {i}");
            }
        }

        /// <summary>
        /// Verifies Estimate() is bit-reproducible across calls at the same seed. Compares raw
        /// bits: a tolerance assert cannot detect a reduction-order difference.
        /// </summary>
        [TestMethod]
        public void Test_Estimate_IsBitReproducible()
        {
            var probabilities = new double[] { 0.999, 0.99, 0.9, 0.5, 0.1, 0.01, 0.001 };
            var dist = new Normal(3.122599, 0.5573654);

            var first = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100, 1000).Estimate(probabilities);
            var second = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100, 1000).Estimate(probabilities);

            Assert.HasCount(first.MeanCurve.Length, second.MeanCurve);
            for (int i = 0; i < first.MeanCurve.Length; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(first.MeanCurve[i]),
                    BitConverter.DoubleToInt64Bits(second.MeanCurve[i]), $"MeanCurve differs at index {i}.");
            }
            for (int i = 0; i < first.ConfidenceIntervals.GetLength(0); i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(first.ConfidenceIntervals[i, j]),
                        BitConverter.DoubleToInt64Bits(second.ConfidenceIntervals[i, j]), $"ConfidenceIntervals differ at [{i},{j}].");
                }
            }
        }

        /// <summary>
        /// Verifies the mean curve does not depend on how many threads compute it. Constraining
        /// one arm to a single worker stands in for running on a machine with a different core
        /// count.
        /// </summary>
        [TestMethod]
        public void Test_ExpectedProbabilities_IsThreadCountIndependent()
        {
            var probabilities = new double[] { 0.99, 0.9, 0.5, 0.1, 0.01 };
            var quantiles = new double[] { 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100, 1000);
            var distributions = boot.Distributions();
            Assert.AreEqual(0, boot.FailedReplications);

            var parallelResult = boot.ExpectedProbabilities(quantiles, probabilities, distributions);

            double[] serialResult;
            ThreadPool.GetMinThreads(out int minWorker, out int minIO);
            try
            {
                ThreadPool.SetMinThreads(1, minIO);
                serialResult = boot.ExpectedProbabilities(quantiles, probabilities, distributions);
            }
            finally
            {
                ThreadPool.SetMinThreads(minWorker, minIO);
            }

            for (int i = 0; i < parallelResult.Length; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(parallelResult[i]),
                    BitConverter.DoubleToInt64Bits(serialResult[i]), $"Expected probabilities differ at index {i}.");
            }
        }

        /// <summary>
        /// Pins the bias-correction proportion count(θ* ≤ θ̂) / (B + 1). With the estimate at the
        /// bootstrap median of nine replicates the proportion is 5/10, the bias correction is zero,
        /// and the bias-corrected limits equal the unadjusted percentile limits — any shift of the
        /// count numerator breaks this identity.
        /// </summary>
        [TestMethod]
        public void Test_BiasCorrected_ZeroBias_MatchesPercentileLimits()
        {
            var parent = new Normal(5d, 1d);
            var boot = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 100);
            var values = new double[] { 1d, 2d, 3d, 4d, 4.5, 6d, 7d, 8d, 9d };
            var replicates = new IUnivariateDistribution[values.Length];
            for (int i = 0; i < values.Length; i++)
                replicates[i] = new Deterministic(values[i]);

            double alpha = 0.1;
            var ci = boot.BiasCorrectedQuantileCI(new double[] { 0.5 }, alpha, replicates);

            Assert.AreEqual(Statistics.Percentile(values, alpha / 2d, true), ci[0, 0], 1E-12);
            Assert.AreEqual(Statistics.Percentile(values, 1d - alpha / 2d, true), ci[0, 1], 1E-12);
        }

        /// <summary>
        /// Pins the full bias-corrected limit formula at an asymmetric proportion: two of nine
        /// replicates at or below the estimate give proportion 2/10, and the limits follow
        /// Φ(2·Φ⁻¹(0.2) + z) exactly. A shifted numerator (3/10) or a B denominator (2/9) breaks it.
        /// </summary>
        [TestMethod]
        public void Test_BiasCorrected_AsymmetricProportion_MatchesFormula()
        {
            var parent = new Normal(2.5, 1d);
            var boot = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 100);
            var values = new double[] { 1d, 2d, 3d, 4d, 4.5, 6d, 7d, 8d, 9d };
            var replicates = new IUnivariateDistribution[values.Length];
            for (int i = 0; i < values.Length; i++)
                replicates[i] = new Deterministic(values[i]);

            double alpha = 0.1;
            var ci = boot.BiasCorrectedQuantileCI(new double[] { 0.5 }, alpha, replicates);

            double bias = Normal.StandardZ(2d / 10d);
            double lower = Statistics.Percentile(values, Normal.StandardCDF(2d * bias + Normal.StandardZ(alpha / 2d)), true);
            double upper = Statistics.Percentile(values, Normal.StandardCDF(2d * bias + Normal.StandardZ(1d - alpha / 2d)), true);
            Assert.AreEqual(lower, ci[0, 0], 1E-12);
            Assert.AreEqual(upper, ci[0, 1], 1E-12);
        }

        /// <summary>
        /// A sampled distribution that failed to fit contributes a NaN parameter set rather than a
        /// null entry, matching BootstrapAnalysis.ParameterSets. Consumers index the array directly.
        /// </summary>
        [TestMethod]
        public void Test_ProcessParameterSets_FillsFailuresWithNaN()
        {
            var probabilities = new double[] { 0.99, 0.9, 0.5, 0.1, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100, 200);
            var sampled = boot.Distributions().Cast<UnivariateDistributionBase>().ToArray();
            sampled[7] = null;

            var results = new UncertaintyAnalysisResults(dist, sampled, probabilities, recordParameterSets: true);
            var reference = boot.ParameterSets(sampled.Cast<IUnivariateDistribution>().ToArray());

            Assert.HasCount(sampled.Length, results.ParameterSets);
            for (int i = 0; i < sampled.Length; i++)
            {
                Assert.IsNotNull(results.ParameterSets[i], $"Parameter set {i} is null.");
                Assert.HasCount(dist.NumberOfParameters, results.ParameterSets[i].Values);
                for (int j = 0; j < dist.NumberOfParameters; j++)
                {
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference[i].Values[j]),
                        BitConverter.DoubleToInt64Bits(results.ParameterSets[i].Values[j]),
                        $"Parameter set {i} value {j} differs from the BootstrapAnalysis form.");
                }
            }
            Assert.IsTrue(double.IsNaN(results.ParameterSets[7].Values[0]));
        }

        /// <summary>
        /// Test that summary probabilities average only the successful fits, and that a
        /// replication set with no successful fit is rejected loudly — as a single error for
        /// the all-null expected-probability path and as an aggregate of the per-set failures
        /// for the distribution builder.
        /// </summary>
        [TestMethod]
        public void Test_UsesOnlySuccessfulFits_AndRejectsAllFailures()
        {
            var parent = new Normal(0d, 1d);
            var analysis = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 1234);
            IUnivariateDistribution[] mixed = { new Normal(0d, 1d), null!, new Normal(1d, 1d) };

            double[] mean = analysis.ExpectedProbabilities(new[] { 0d }, mixed);
            double expected = 0.5d * (new Normal(0d, 1d).CDF(0d) + new Normal(1d, 1d).CDF(0d));
            Assert.AreEqual(expected, mean[0], 1E-14);
            Assert.Throws<InvalidOperationException>(() =>
                analysis.ExpectedProbabilities(new[] { 0d }, new IUnivariateDistribution[] { null!, null! }));

            var aggregate = Assert.Throws<AggregateException>(() => analysis.Distributions(new[]
            {
                new ParameterSet(new[] { 0d, -1d }, 0d),
                new ParameterSet(new[] { double.NaN, 1d }, 0d),
            }));
            Assert.HasCount(2, aggregate.InnerExceptions);
        }

        /// <summary>
        /// Verifies that the released one-argument bootstrap summary method tokens remain
        /// available to already-compiled consumers.
        /// </summary>
        [TestMethod]
        public void Test_OneArgumentSummaryMethods_RetainBinarySignatures()
        {
            Type listType = typeof(IList<double>);

            Assert.IsNotNull(typeof(BootstrapAnalysis).GetMethod(
                nameof(BootstrapAnalysis.Quantiles), new[] { listType }));
            Assert.IsNotNull(typeof(BootstrapAnalysis).GetMethod(
                nameof(BootstrapAnalysis.Probabilities), new[] { listType }));
        }

        /// <summary>
        /// Test that the normal-approximation quantile interval preserves the sign of negative
        /// quantiles: the transform applied around the point estimate must remain finite and
        /// keep both interval endpoints on the data's side of zero.
        /// </summary>
        [TestMethod]
        public void Test_NormalQuantileCI_PreservesNegativeQuantiles()
        {
            var parent = new Normal(-10d, 1d);
            var analysis = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 1234);
            IUnivariateDistribution[] fits =
            {
                new Normal(-9.5d, 1d),
                new Normal(-10d, 1.1d),
                null!,
                new Normal(-10.5d, 0.9d),
            };

            double[,] interval = analysis.NormalQuantileCI(new[] { 0.5d }, 0.1d, fits);
            Assert.IsTrue(Tools.IsFinite(interval[0, 0]));
            Assert.IsTrue(Tools.IsFinite(interval[0, 1]));
            Assert.IsLessThan(0d, interval[0, 0]);
            Assert.IsLessThan(0d, interval[0, 1]);
        }

        /// <summary>
        /// Test that expected probabilities pair each quantile with its own probability after
        /// the internal sort: unsorted input ordinates produce exactly the same value array as
        /// the pre-sorted equivalent.
        /// </summary>
        [TestMethod]
        public void Test_ExpectedProbabilities_InterpolationKeepsPairsSorted()
        {
            var parent = new Normal(0d, 1d);
            var analysis = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 1234);
            IUnivariateDistribution[] fits = { new Normal(0d, 1d), new Normal(1d, 2d) };
            double[] sorted = { -3d, -1d, 0d, 2d, 5d };
            double[] unsorted = { 2d, -3d, 5d, 0d, -1d };
            double[] probabilities = { 0.1d, 0.5d, 0.9d };

            double[] expected = analysis.ExpectedProbabilities(sorted, probabilities, fits);
            double[] actual = analysis.ExpectedProbabilities(unsorted, probabilities, fits);
            CollectionAssert.AreEqual(expected, actual);
        }

        /// <summary>
        /// Test that Quantiles dimensions its output by the supplied distributions array — not by
        /// the replication count — with orientation [replication, ordinate], writes a NaN row for
        /// each null entry, and fills every non-null entry's row with that distribution's own
        /// InverseCDF at the requested probabilities.
        /// </summary>
        [TestMethod]
        public void Test_Quantiles_SuppliedDistributions_DimensionsNaNRowsAndValues()
        {
            var probabilities = new double[] { 0.1d, 0.5d, 0.9d };
            var parent = new Normal(0d, 1d);
            var boot = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 200);
            IUnivariateDistribution[] fits = { new Normal(0d, 1d), null!, new Normal(1d, 2d) };

            double[,] result = boot.Quantiles(probabilities, fits);

            Assert.AreEqual(3, result.GetLength(0));
            Assert.AreEqual(probabilities.Length, result.GetLength(1));
            for (int j = 0; j < probabilities.Length; j++)
            {
                Assert.IsTrue(double.IsNaN(result[1, j]), $"The null entry's row must be NaN at column {j}.");
                Assert.IsTrue(Tools.IsFinite(result[0, j]));
                Assert.IsTrue(Tools.IsFinite(result[2, j]));
                Assert.AreEqual(fits[0].InverseCDF(probabilities[j]), result[0, j], 0d);
                Assert.AreEqual(fits[2].InverseCDF(probabilities[j]), result[2, j], 0d);
            }
        }

        /// <summary>
        /// Test that Probabilities dimensions its output by the supplied distributions array — not
        /// by the replication count — with orientation [replication, ordinate], writes a NaN row
        /// for each null entry, and fills every non-null entry's row with that distribution's own
        /// CDF at the requested quantiles.
        /// </summary>
        [TestMethod]
        public void Test_Probabilities_SuppliedDistributions_DimensionsNaNRowsAndValues()
        {
            var quantiles = new double[] { -1d, 0.5d, 2d };
            var parent = new Normal(0d, 1d);
            var boot = new BootstrapAnalysis(parent, ParameterEstimationMethod.MethodOfMoments, 10, 200);
            IUnivariateDistribution[] fits = { new Normal(0d, 1d), null!, new Normal(1d, 2d) };

            double[,] result = boot.Probabilities(quantiles, fits);

            Assert.AreEqual(3, result.GetLength(0));
            Assert.AreEqual(quantiles.Length, result.GetLength(1));
            for (int j = 0; j < quantiles.Length; j++)
            {
                Assert.IsTrue(double.IsNaN(result[1, j]), $"The null entry's row must be NaN at column {j}.");
                Assert.IsTrue(Tools.IsFinite(result[0, j]));
                Assert.IsTrue(Tools.IsFinite(result[2, j]));
                Assert.AreEqual(fits[0].CDF(quantiles[j]), result[0, j], 0d);
                Assert.AreEqual(fits[2].CDF(quantiles[j]), result[2, j], 0d);
            }
        }

    }
}
