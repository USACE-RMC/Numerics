/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// A class for performing the bootstrap uncertainty analysis.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Bootstrapping_(statistics)" />
    /// </para>
    /// </remarks>
    public class BootstrapAnalysis
    {

        #region Construction

        /// <summary>
        /// Construct a new Bootstrap Analysis.
        /// </summary>
        /// <param name="distribution">The univariate distribution to bootstrap.</param>
        /// <param name="estimationMethod">The parameter estimation method.</param>
        /// <param name="sampleSize">Size of the bootstrap sample to generate.</param>
        /// <param name="replications">The number of bootstrap replications to be sampled.</param>
        /// <param name="seed">Optional. Seed for random number generator. Default = 12345.</param>
        public BootstrapAnalysis(IUnivariateDistribution distribution, ParameterEstimationMethod estimationMethod, int sampleSize, int replications = 10000, int seed = 12345)
        {
            if (distribution as IBootstrappable == null) throw new ArgumentException("The distribution must implement IBootstrappable.", nameof(Distribution));
            if (sampleSize < 10) throw new ArgumentOutOfRangeException(nameof(SampleSize), "The sample size must at least 10.");
            if (replications < 100) throw new ArgumentOutOfRangeException(nameof(Replications), "The number of bootstrap replications must be at least 100.");

            Distribution = (IBootstrappable)distribution;
            EstimationMethod = estimationMethod;
            SampleSize = sampleSize;
            Replications = replications;
            PRNGSeed = seed;
        }

        #endregion

        #region Members

        private int _retries = 20;

        /// <summary>
        /// The distribution parameter estimation method.
        /// </summary>
        public ParameterEstimationMethod EstimationMethod { get; private set; }

        /// <summary>
        /// The univariate distribution to bootstrap.
        /// </summary>
        public IBootstrappable Distribution { get; private set; }

        /// <summary>
        /// Size of the bootstrap sample to generate.
        /// </summary>
        public int SampleSize { get; private set; }

        /// <summary>
        /// The number of bootstrap replications to be sampled.
        /// </summary>
        public int Replications { get; private set; }

        /// <summary>
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; private set; }

        #endregion

        #region Methods

        /// <summary>
        /// Bootstrap a list of fitted distributions.
        /// </summary>
        public IUnivariateDistribution[] Distributions()
        {
            var bootDistributions = new IUnivariateDistribution[Replications];
            var r = new MersenneTwister(PRNGSeed);
            var seeds = r.NextIntegers(Replications);
            Parallel.For(0, Replications, idx =>
            {
                bool failed = false;
                for (int m = 0; m < _retries; m++)
                {
                    try
                    {
                        bootDistributions[idx] = Distribution.Bootstrap(EstimationMethod, SampleSize, seeds[idx] + 10 * m);
                        failed = false;
                    }
                    catch (Exception)
                    {
                        failed = true;
                    };

                    if (failed == false) break;
                }

                // MLE and certain L-moments methods can fail to find a solution
                // On fail, set to null
                if (failed == true) bootDistributions[idx] = null;

            });
            return bootDistributions;
        }

        /// <summary>
        /// Return a list of distributions given an array of parameter sets.
        /// </summary>
        /// <param name="parameterSets">An array of parameter sets.</param>
        public IUnivariateDistribution[] Distributions(ParameterSet[] parameterSets)
        {
            var bootDistributions = new IUnivariateDistribution[parameterSets.Length];
            Parallel.For(0, parameterSets.Length, idx =>
            {
                bool failed = false;

                try
                {
                    var dist = ((UnivariateDistributionBase)Distribution).Clone();
                    dist.SetParameters(parameterSets[idx].Values);
                    bootDistributions[idx] = dist;
                    failed = false;
                }
                catch (Exception)
                {
                    failed = true;
                };

                // On fail, set to null
                if (failed == true) bootDistributions[idx] = null;

            });
            return bootDistributions;
        }


        /// <summary>
        /// Bootstrap an array of distribution parameters.
        /// </summary>
        public double[,] Parameters(IUnivariateDistribution[] distributions = null)
        {
            var bootDistributions = distributions != null ? distributions : Distributions();
            var bootParameters = new double[bootDistributions.Count(), Distribution.NumberOfParameters];
            Parallel.For(0, bootDistributions.Count(), idx =>
            {
                if (bootDistributions[idx] != null)
                {
                    var parameters = bootDistributions[idx].GetParameters;
                    for (int i = 0; i < Distribution.NumberOfParameters; i++)
                        bootParameters[idx, i] = parameters[i];
                }
                else
                {
                    for (int i = 0; i < Distribution.NumberOfParameters; i++)
                        bootParameters[idx, i] = double.NaN;
                }

            });
            return bootParameters;
        }

        /// <summary>
        /// Bootstrap an array of distribution parameter sets.
        /// </summary>
        public ParameterSet[] ParameterSets(IUnivariateDistribution[] distributions = null)
        {
            var bootDistributions = distributions != null ? distributions : Distributions();
            var bootParameters = new ParameterSet[bootDistributions.Count()];
            Parallel.For(0, bootDistributions.Count(), idx =>
            {
                if (bootDistributions[idx] != null)
                {
                    bootParameters[idx] = new ParameterSet(bootDistributions[idx].GetParameters, double.NaN);
                }
                else
                {
                    var parameters = new double[Distribution.NumberOfParameters];
                    for (int i = 0; i < Distribution.NumberOfParameters; i++)
                        parameters[i] = double.NaN;
                    bootParameters[idx] = new ParameterSet(parameters, double.NaN);
                }

            });
            return bootParameters;
        }

        /// <summary>
        /// Bootstrap a list of product moments for each bootstrapped sample.
        /// </summary>
        public double[,] ProductMoments()
        {
            var bootMoments = new double[Replications, 4];
            var r = new MersenneTwister(PRNGSeed);
            var seeds = r.NextIntegers(Replications);
            Parallel.For(0, Replications, idx =>
            {
                var moments = Statistics.ProductMoments(Distribution.GenerateRandomValues(SampleSize, seeds[idx]));
                for (int i = 0; i < moments.Length; i++) bootMoments[idx, i] = moments[i];
            });
            return bootMoments;
        }

        /// <summary>
        /// Bootstrap a list of linear moments for each bootstrapped sample.
        /// </summary>
        public double[,] LinearMoments()
        {
            var bootMoments = new double[Replications, 4];
            var r = new MersenneTwister(PRNGSeed);
            var seeds = r.NextIntegers(Replications);
            Parallel.For(0, Replications, idx =>
            {
                var moments = Statistics.LinearMoments(Distribution.GenerateRandomValues(SampleSize, seeds[idx]));
                for (int i = 0; i < moments.Length; i++) bootMoments[idx, i] = moments[i];
            });
            return bootMoments;
        }

        /// <summary>
        /// Bootstrap a list of quantiles given the input non-exceedance probabilities.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        public double[,] Quantiles(IList<double> probabilities)
        {
            var Output = new double[Replications, probabilities.Count];
            var bootDistributions = Distributions();
            for (int i = 0; i < probabilities.Count; i++)
                Parallel.For(0, Replications, idx => { Output[idx, i] = bootDistributions[idx] != null ? bootDistributions[idx].InverseCDF(probabilities[i]) : double.NaN; });
            return Output;
        }

        /// <summary>
        /// Bootstrap a list of non-exceedance probabilities given the input quantile values.
        /// </summary>
        /// <param name="quantiles">List quantile values.</param>
        public double[,] Probabilities(IList<double> quantiles)
        {
            var Output = new double[Replications, quantiles.Count];
            var bootDistributions = Distributions();
            for (int i = 0; i < quantiles.Count; i++)
                Parallel.For(0, Replications, idx => { Output[idx, i] = bootDistributions[idx] != null ? bootDistributions[idx].CDF(quantiles[i]) : double.NaN; });
            return Output;
        }

        /// <summary>
        /// Bootstrap full uncertainty analysis results using the percentile method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        /// <param name="recordParameterSets">Optional. Determines whether to record parameter sets. Default = true.</param>
        public UncertaintyAnalysisResults Estimate(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[] distributions = null, bool recordParameterSets = true)
        {
            var results = new UncertaintyAnalysisResults();
            results.ParentDistribution = (UnivariateDistributionBase)Distribution;

            // get mode curve
            results.ModeCurve = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                results.ModeCurve[i] = Distribution.InverseCDF(probabilities[i]);

            // get bootstrapped list of distributions
            var bootDistributions = distributions != null ? distributions : Distributions();

            // get parameter sets
            if (recordParameterSets == true)
                results.ParameterSets = ParameterSets(bootDistributions);

            // get confidence intervals
            results.ConfidenceIntervals = PercentileQuantileCI(probabilities, alpha, bootDistributions);
     
            // create list of quantiles
            var minMax = ComputeMinMaxQuantiles(0.001, 1 - 1E-9, bootDistributions);
            int bins = 200;
            List<double> quantiles = new List<double>();
            double shift = 0;
            if (minMax[0] <= 0) shift = Math.Abs(minMax[0]) + 1d;
            double min = minMax[0] + shift;
            double max = minMax[1] + shift;
            int order = (int)Math.Floor(Math.Log10(max) - Math.Log10(min));
            bins = Math.Max(200, Math.Min(1000, 100 * order));
            double delta = (Math.Log10(max) - Math.Log10(min)) / (bins - 1);
            double x = Math.Log10(min);
            quantiles.Add(Math.Pow(10, x) - shift);
            for (int i = 1; i <= bins - 1; i++)
            {
                x = Math.Log10(quantiles[i - 1] + shift) + delta;
                quantiles.Add(Math.Pow(10, x) - shift);
            }

            // get mean curve
            results.MeanCurve = ExpectedProbabilities(quantiles, probabilities, bootDistributions);

            return results;
        }

        /// <summary>
        /// Bootstrap the expected non-exceedance probabilities given the input quantile values. Returns the x-values interpolated from the list of desired non-exceedance probabilities.
        /// </summary>
        /// <param name="quantiles">List quantile values.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[] ExpectedProbabilities(IList<double> quantiles, IList<double> probabilities, IUnivariateDistribution[] distributions = null)
        {
            var quants = quantiles.ToArray();
            var probs = probabilities.ToArray();
            Array.Sort(quants);
            var expected = new double[quantiles.Count];
            var output = new double[probabilities.Count];
            var bootDistributions = distributions != null ? distributions : Distributions();
            for (int i = 0; i < quantiles.Count; i++)
            {
                double total = 0d;
                Parallel.For(0, bootDistributions.Count(), () => 0d, (j, loop, sum) =>
                {
                    if (bootDistributions[j] != null)
                    {
                        sum += bootDistributions[j].CDF(quants[i]);
                    }
                    return sum;
                }, z => Tools.ParallelAdd(ref total, z));
                expected[i] = total / bootDistributions.Count();
            }

            double minY = double.MaxValue;
            double maxY = double.MinValue;
            var yVals = new List<double>();
            var xVals = new List<double>();
            yVals.Add(quantiles[0]);
            xVals.Add(expected[0]);
            for (int i = 1; i < quantiles.Count; i++)
            {
                if (expected[i] > xVals.Last())
                {
                    minY = Math.Min(minY, quantiles[i]);
                    maxY = Math.Max(maxY, quantiles[i]);
                    yVals.Add(quantiles[i]);
                    xVals.Add(expected[i]);
                }
            }
            bool useLogTransform = false;
            if (minY > 0 && (Math.Log10(maxY) - Math.Log10(minY)) > 1)
                useLogTransform = true;

            Linear linint = new Linear(xVals, yVals) { XTransform = Transform.NormalZ, YTransform = useLogTransform ?  Transform.Logarithmic : Transform.None };
            output = linint.Interpolate(probs);
            return output;
        }

        /// <summary>
        /// Bootstrap the expected non-exceedance probabilities given the input quantile values.
        /// </summary>
        /// <param name="quantiles">List quantile values.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[] ExpectedProbabilities(IList<double> quantiles, IUnivariateDistribution[] distributions = null)
        {
            var quants = quantiles.ToArray();
            Array.Sort(quants);
            var expected = new double[quantiles.Count];
            var bootDistributions = distributions != null ? distributions : Distributions();
            for (int i = 0; i < quantiles.Count; i++)
            {
                double total = 0d;
                Parallel.For(0, bootDistributions.Count(), () => 0d, (j, loop, sum) =>
                {
                    if (bootDistributions[j] != null)
                    {
                        sum += bootDistributions[j].CDF(quants[i]);
                    }
                    return sum;
                }, z => Tools.ParallelAdd(ref total, z));
                expected[i] = total / bootDistributions.Count();
            }
            return expected;
        }

        /// <summary>
        /// Returns the min and max quantiles from a bootstrap analysis.
        /// </summary>
        /// <param name="minProbability">The minimum probability to compute quantiles.</param>
        /// <param name="maxProbability">The maximum probability to compute quantiles.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[] ComputeMinMaxQuantiles(double minProbability, double maxProbability, IUnivariateDistribution[] distributions)
        {
            var output = new double[] { double.MaxValue, double.MinValue };
            Parallel.For(0, distributions.Count(), j =>
            {
                if (distributions[j] != null)
                {
                    var minX = distributions[j].InverseCDF(minProbability);
                    if (minX < output[0]) output[0] = minX;

                    var maxX = distributions[j].InverseCDF(maxProbability);
                    if (maxX > output[1]) output[1] = maxX;
                }
            });
            return output;
        }

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the percentile method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[,] PercentileQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[] distributions = null)
        {
            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];
            var bootDistributions = distributions != null ? distributions : Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var XValues = new double[bootDistributions.Count()];
                Parallel.For(0, bootDistributions.Count(), idx => { XValues[idx] = bootDistributions[idx] != null ? bootDistributions[idx].InverseCDF(probabilities[i]) : double.NaN; });

                // sort X values
                var validValues = XValues.Where(x => !double.IsNaN(x)).ToArray();
                Array.Sort(validValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                    Output[i, j] = Statistics.Percentile(validValues, CIs[j], true);
            }
            return Output;
        }

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the bias-corrected percentile method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[,] BiasCorrectedQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[] distributions = null)
        {
            // Create list of original X values given probability values
            var populationXValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationXValues[i] = Distribution.InverseCDF(probabilities[i]);

            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];
            var bootDistributions = distributions != null ? distributions : Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                double P0 = 0d; // proportions of values less than population
                var XValues = new double[bootDistributions.Count()];
                Parallel.For(0, bootDistributions.Count(), () => 0d, (idx, loop, subP0) =>
                {
                    XValues[idx] = bootDistributions[idx] != null ? bootDistributions[idx].InverseCDF(probabilities[i]) : double.NaN;
                    if (XValues[idx] != double.NaN && XValues[idx] <= populationXValues[i]) subP0 += 1d;
                    return subP0;
                }, z => Tools.ParallelAdd(ref P0, z));

                // get proportion
                P0 = P0 / (bootDistributions.Count() + 1);

                // sort X values
                var validValues = XValues.Where(x => !double.IsNaN(x)).ToArray();
                Array.Sort(validValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                {
                    double Z0 = Normal.StandardZ(P0);
                    double Z = Normal.StandardZ(CIs[j]);
                    double BC = Normal.StandardCDF(2d * Z0 + Z);
                    Output[i, j] = Statistics.Percentile(validValues, BC, true);
                }
            }
            return Output;
        }

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the Normal, or standard method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        /// <param name="distributions">Optional. Pass in an array of bootstrapped distributions. Default = null.</param>
        public double[,] NormalQuantileCI(IList<double> probabilities, double alpha = 0.1, IUnivariateDistribution[] distributions = null)
        {

            // Create list of original X values given probability values
            // Use a cube-root transform to make results transformation invariant
            var populationXValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationXValues[i] = Math.Pow(Distribution.InverseCDF(probabilities[i]), 1d / 3d);

            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];
            var bootDistributions = distributions != null ? distributions : Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                var XValues = new double[bootDistributions.Count()];
                Parallel.For(0, bootDistributions.Count(), idx => { XValues[idx] = bootDistributions[idx] != null ? Math.Pow(bootDistributions[idx].InverseCDF(probabilities[i]), 1d / 3d) : double.NaN; });

                // Get Standard error
                var validValues = XValues.Where(x => !double.IsNaN(x)).ToArray();
                double SE = Statistics.StandardDeviation(validValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                {
                    double Z = Normal.StandardZ(CIs[j]);
                    Output[i, j] = Math.Pow(populationXValues[i] + SE * Z, 3d);
                }
            }
            return Output;
        }

        #region Bias-Corrected and Accelerated

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the bias-corrected and accelerated (BCa) percentile method.
        /// </summary>
        /// <param name="sampleData">Sample of data.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,] BCaQuantileCI(IList<double> sampleData, IList<double> probabilities, double alpha = 0.1)
        {
            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];

            // Estimate distribution
            SampleSize = sampleData.Count;
            ((IEstimation)Distribution).Estimate(sampleData, EstimationMethod);

            // Create list of original X values given probability values
            var populationXValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationXValues[i] = Distribution.InverseCDF(probabilities[i]);

            // Get acceleration constants
            var a = AccelerationConstants(sampleData, probabilities, populationXValues);

            // Get bootstrapped distributions
            var bootDistributions = Distributions();
            for (int i = 0; i < probabilities.Count; i++)
            {
                double P0 = 0d; // proportions of values less than population
                var XValues = new double[Replications];
                Parallel.For(0, Replications, () => 0d, (idx, loop, subP0) =>
                {
                    XValues[idx] = bootDistributions[idx] != null ? bootDistributions[idx].InverseCDF(probabilities[i]) : double.NaN;
                    if (XValues[idx] != double.NaN && XValues[idx] <= populationXValues[i]) subP0 += 1d;
                    return subP0;
                }, z => Tools.ParallelAdd(ref P0, z));

                // get proportion
                P0 = P0 / (Replications + 1);

                // sort X values
                var validValues = XValues.Where(x => !double.IsNaN(x)).ToArray();
                Array.Sort(validValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                {
                    double Z0 = Normal.StandardZ(P0);
                    double Z = Normal.StandardZ(CIs[j]);
                    double num = Z0 + Z;
                    double den = 1 - a[i] * (Z0 + Z);
                    double BC = Normal.StandardCDF(Z0 + num / den);
                    Output[i, j] = Statistics.Percentile(validValues, BC, true);
                }
            }
            return Output;
        }

        /// <summary>
        /// Estimates the acceleration constants for each probability.
        /// </summary>
        /// <param name="sampleData">Sample of data.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="thetaHats">The list of best-estimate quantiles.</param>
        private double[] AccelerationConstants(IList<double> sampleData, IList<double> probabilities, IList<double> thetaHats)
        {
            var N = sampleData.Count;
            var I2 = new double[probabilities.Count];
            var I3 = new double[probabilities.Count];
            var a = new double[probabilities.Count];

            // Perform Jackknife
            Parallel.For(0, N, idx =>
            {
                // Remove data point
                var jackSample = new List<double>(sampleData);
                jackSample.RemoveAt(idx);

                // Estimate distribution
                var newDistribution = ((UnivariateDistributionBase)Distribution).Clone();

                try
                {
                    ((IEstimation)newDistribution).Estimate(jackSample, EstimationMethod);
                    // Get quantiles from new distribution
                    var thetaJack = new double[probabilities.Count];
                    for (int i = 0; i < probabilities.Count; i++)
                    {
                        thetaJack[i] = newDistribution.InverseCDF(probabilities[i]);
                        I2[i] += Math.Pow((N - 1) * (thetaHats[i] - thetaJack[i]), 2);
                        I3[i] += Math.Pow((N - 1) * (thetaHats[i] - thetaJack[i]), 3);
                    }
                }
                catch (Exception)
                {
                    // MLE and certain L-moments methods can fail to find a solution
                };

            });
            // Get acceleration constant
            for (int i = 0; i < probabilities.Count; i++)
                a[i] = I3[i] / (Math.Pow(I2[i], 1.5) * 6);

            return a;
        }

        #endregion

        #region Bootstrap-t (aka Student-t Bootstrap)

        /// <summary>
        /// Bootstrap confidence intervals for a list of quantiles using the Bootstrap-t method.
        /// </summary>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,] BootstrapTQuantileCI(IList<double> probabilities, double alpha = 0.1)
        {
            // Create list of original X values given probability values
            // Use a cube-root transform to make results transformation invariant
            var populationXValues = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                populationXValues[i] = Math.Pow(Distribution.InverseCDF(probabilities[i]), 1d / 3d);

            var xValues = new double[Replications, probabilities.Count];
            var studentT = new double[Replications, probabilities.Count];
            var CIs = new double[] { alpha / 2d, 1d - alpha / 2d };
            var Output = new double[probabilities.Count, 2];

            // First create list of bootstrap distributions, 
            // and estimate standard error for each quantiles          
            var bootDistributions = new IUnivariateDistribution[Replications];
            var r = new MersenneTwister(PRNGSeed);
            var seeds = r.NextIntegers(Replications);
            Parallel.For(0, Replications, i =>
            {
                try
                {
                    var newDistribution = ((UnivariateDistributionBase)Distribution).Clone();
                    var sample = newDistribution.GenerateRandomValues(SampleSize, seeds[i]);
                    ((IEstimation)newDistribution).Estimate(sample, EstimationMethod);
                    bootDistributions[i] = newDistribution;

                    // Record inner boot thetas
                    var bootXValues = new double[probabilities.Count];
                    for (int j = 0; j < probabilities.Count; j++)
                        bootXValues[j] = Math.Pow(bootDistributions[i].InverseCDF(probabilities[j]), 1d / 3d);

                    // Now estimate the standard error at each quantile using the jackknife method
                    var bootSE = StandardError(sample, probabilities, bootXValues);
                    for (int j = 0; j < probabilities.Count; j++)
                    {
                        xValues[i, j] = bootXValues[j];
                        studentT[i, j] = (populationXValues[j] - bootXValues[j]) / bootSE[j];
                    }

                }
                catch (Exception)
                {
                    // MLE and certain L-moments methods can fail to find a solution
                    // On fail, set to null
                    bootDistributions[i] = null;
                    for (int j = 0; j < probabilities.Count; j++)
                    {
                        xValues[i, j] = double.NaN;
                        studentT[i, j] = double.NaN;
                    }
                };

            });


            for (int i = 0; i < probabilities.Count; i++)
            {
                var XValues = xValues.GetColumn(i).Where(x => !double.IsNaN(x)).ToArray();
                var TValues = studentT.GetColumn(i).Where(x => !double.IsNaN(x)).ToArray();

                // Get Standard error
                double SE = Statistics.StandardDeviation(XValues);
                Array.Sort(TValues);

                // Record percentiles for CIs
                for (int j = 0; j < 2; j++)
                {
                    double T = Statistics.Percentile(TValues, CIs[j], true);
                    Output[i, j] = Math.Pow(populationXValues[i] + SE * T, 3d);
                }
            }

            return Output;
        }

        /// <summary>
        /// Estimates the standard error for each probability.
        /// </summary>
        /// <param name="sampleData">Sample of data.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="thetaHats">The list of best-estimate quantiles.</param>
        private double[] StandardError(IList<double> sampleData, IList<double> probabilities, IList<double> thetaHats)
        {
            var N = sampleData.Count;
            var I2 = new double[probabilities.Count];
            var se = new double[probabilities.Count];

            // Perform Jackknife
            Parallel.For(0, N, idx =>
            {
                // Remove data point
                var jackSample = new List<double>(sampleData);
                jackSample.RemoveAt(idx);
                // Estimate distribution
                var newDistribution = ((UnivariateDistributionBase)Distribution).Clone();

                try
                {
                    ((IEstimation)newDistribution).Estimate(jackSample, EstimationMethod);
                    // Get quantiles from new distribution
                    var thetaJack = new double[probabilities.Count];
                    for (int i = 0; i < probabilities.Count; i++)
                    {
                        thetaJack[i] = Math.Pow(newDistribution.InverseCDF(probabilities[i]), 1d / 3d);
                        I2[i] += Math.Pow(thetaHats[i] - thetaJack[i], 2);
                    }
                }
                catch (Exception)
                {
                    // MLE and certain L-moments methods can fail to find a solution
                    // On fail, set to null
                    for (int i = 0; i < probabilities.Count; i++)
                        I2[i] += 0;
                };

            });
            // Get standard error
            for (int i = 0; i < probabilities.Count; i++)
                se[i] = Math.Sqrt((N - 1) / (double)N * I2[i]);

            return se;
        }

        #endregion

        #endregion

    }
}