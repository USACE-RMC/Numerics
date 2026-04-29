using Numerics.Mathematics.Optimization;
using Numerics.Utilities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace Numerics.Sampling.MCMC
{
    /// <summary>
    /// A class for post-processing and saving Bayesian MCMC results.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
   [Serializable] 
    public class MCMCResults
    {
        /// <summary>
        /// Constructs an empty MCMC results.
        /// </summary>
        public MCMCResults() { }

        /// <summary>
        /// Constructs and post-processes MCMC results. 
        /// </summary>
        /// <param name="sampler">The MCMC sampler to post-process.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param> 
        public MCMCResults(MCMCSampler sampler, double alpha = 0.1)
        {
            // Clone the Markov Chains and Output
            MarkovChains = new List<ParameterSet>[sampler.NumberOfChains];
            Output = new List<ParameterSet>();
            for (int i = 0; i < sampler.NumberOfChains; i++)
            {
                MarkovChains[i] = sampler.MarkovChains[i].ToList();
                Output.AddRange(sampler.Output[i].ToList());
            }
            AcceptanceRates = sampler.AcceptanceRates.ToArray();
            MeanLogLikelihood = sampler.MeanLogLikelihood.ToList();
            MAP = sampler.MAP.Clone();
            ProcessParameterResults(sampler, alpha);
        }

        /// <summary>
        /// Constructs and post-processes MCMC results. 
        /// </summary>
        /// <param name="map">The output parameter set that produced the maximum likelihood.</param>
        /// <param name="parameterSets">The list of parameter sets to process.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param> 
        public MCMCResults(ParameterSet map, IList<ParameterSet> parameterSets, double alpha = 0.1)
        {
            
            MAP = map.Clone();
            Output = parameterSets.ToList();
            ProcessParameterResults(alpha);
        }

        /// <summary>
        /// The list of sampled Markov Chains.
        /// </summary>
        [JsonInclude]
        public List<ParameterSet>[]? MarkovChains { get; private set; }

        /// <summary>
        /// Output posterior parameter sets.
        /// </summary>
        [JsonInclude]
        public List<ParameterSet> Output { get; private set; } = new List<ParameterSet>();

        /// <summary>
        /// The average log-likelihood across each chain for each iteration.
        /// </summary>
        [JsonInclude]
        public List<double>? MeanLogLikelihood { get; private set; }

        /// <summary>
        /// The acceptance rate for each chain.
        /// </summary>
        [JsonInclude]
        public double[] AcceptanceRates { get; private set; } = null!;

        /// <summary>
        /// Parameter results using the output posterior parameter sets.
        /// </summary>
        [JsonInclude]
        public ParameterResults[] ParameterResults { get; private set; } = null!;

        /// <summary>
        /// The output parameter set that produced the maximum likelihood.
        /// This is referred to as the maximum a posteriori (MAP).
        /// </summary>
        [JsonInclude]
        public ParameterSet MAP { get; private set; } = new ParameterSet();

        /// <summary>
        /// The mean of the posterior distribution of each parameter.
        /// </summary>
        [JsonInclude]
        public ParameterSet PosteriorMean { get; private set; } = new ParameterSet();

        /// <summary>
        /// Process the parameter results.
        /// </summary>
        /// <param name="sampler">The MCMC sampler to post-process.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param> 
        private void ProcessParameterResults(MCMCSampler sampler, double alpha = 0.1)
        {
            // Compute the Gelman-Rubin diagnostic using the post-warm up period
            var GR = MCMCDiagnostics.GelmanRubin(sampler.MarkovChains, sampler.WarmupIterations);
            // Compute the effective sample size using the output
            var ESS = MCMCDiagnostics.EffectiveSampleSize(sampler.Output, out var averageACF);
            // Compute parameter summary statistics
            var postMean = new double[sampler.NumberOfParameters];
            ParameterResults = new ParameterResults[sampler.NumberOfParameters];
            for (int i = 0; i < sampler.NumberOfParameters; i++)
            {
                var x = Output.Select(set => set.Values[i]).ToArray();
                ParameterResults[i] = new ParameterResults(x, alpha);
                ParameterResults[i].SummaryStatistics.Rhat = GR[i];
                ParameterResults[i].SummaryStatistics.ESS = ESS[i];
                ParameterResults[i].Autocorrelation = averageACF[i];
                postMean[i] = ParameterResults[i].SummaryStatistics.Mean;
            }
            // Set the posterior mean parameter set. 
            var postMeanLogLH = sampler.LogLikelihoodFunction(postMean);
            PosteriorMean = new ParameterSet(postMean, postMeanLogLH);
        }

        /// <summary>
        /// Process a parameter results using the output list.
        /// </summary>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        private void ProcessParameterResults(double alpha = 0.1)
        {
            int p = Output.First().Values.Length;
            var postMean = new double[p];
            ParameterResults = new ParameterResults[p];
            for (int i = 0; i < p; i++)
            {
                var x = Output.Select(set => set.Values[i]).ToArray();
                ParameterResults[i] = new ParameterResults(x, alpha);
                postMean[i] = ParameterResults[i].SummaryStatistics.Mean;
            }
            // Set the posterior mean parameter set.
            PosteriorMean = new ParameterSet(postMean, double.NaN);
        }

        /// <summary>
        /// Recompute parameter summary statistics at a new credible-interval level
        /// (alpha) without rerunning the chain. Preserves Rhat, ESS, autocorrelation,
        /// MarkovChains, AcceptanceRates, MeanLogLikelihood, MAP, and Output.
        /// </summary>
        /// <param name="alpha">
        /// The new significance level (e.g., 0.05 for 95% credible intervals,
        /// 0.10 for 90% credible intervals).
        /// </param>
        /// <remarks>
        /// Used when the user changes the credible-interval width on a Bayesian
        /// analysis whose MCMC chain has already converged. The chain output is
        /// independent of alpha; only the posterior summary percentiles
        /// (LowerCI/UpperCI) need to be recomputed. Snapshots the per-parameter
        /// Rhat, ESS, and autocorrelation diagnostics (which depend on the chain,
        /// not on alpha) before delegating to the existing private reprocess, then
        /// restores them on the freshly-built ParameterResults instances.
        /// </remarks>
        public void RecomputeParameterResults(double alpha)
        {
            if (ParameterResults == null || Output == null || Output.Count == 0) return;

            // Snapshot diagnostics that are independent of alpha. Length matches
            // the parameter count of the existing ParameterResults array.
            int n = ParameterResults.Length;
            var rhats = new double[n];
            var esss = new double[n];
            var acfs = new double[n][,];
            for (int i = 0; i < n; i++)
            {
                rhats[i] = ParameterResults[i].SummaryStatistics.Rhat;
                esss[i] = ParameterResults[i].SummaryStatistics.ESS;
                acfs[i] = ParameterResults[i].Autocorrelation;
            }

            // Recompute (replaces ParameterResults[] and PosteriorMean — Output untouched)
            ProcessParameterResults(alpha);

            // Restore preserved diagnostics
            for (int i = 0; i < ParameterResults.Length; i++)
            {
                ParameterResults[i].SummaryStatistics.Rhat = rhats[i];
                ParameterResults[i].SummaryStatistics.ESS = esss[i];
                ParameterResults[i].Autocorrelation = acfs[i];
            }
        }

        #region Serialization

        /// <summary>
        /// Converts the MCMC Results to a byte array.
        /// </summary>
        /// <param name="mcmcResults">The MCMC Results.</param>
        public static byte[] ToByteArray(MCMCResults mcmcResults)
        {
            var options = new JsonSerializerOptions
            {
                WriteIndented = false,
                DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
                IncludeFields = true,
                NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals
            };
            options.Converters.Add(new Double2DArrayConverter());
            options.Converters.Add(new HistogramConverter());
            return JsonSerializer.SerializeToUtf8Bytes(mcmcResults, options);
        }

        /// <summary>
        /// Creates MCMC Results from a byte array.
        /// </summary>
        /// <param name="bytes">Byte array.</param>
        public static MCMCResults? FromByteArray(byte[] bytes)
        {
            var options = new JsonSerializerOptions
            {
                DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
                IncludeFields = true,
                NumberHandling = JsonNumberHandling.AllowNamedFloatingPointLiterals
            };
            options.Converters.Add(new Double2DArrayConverter());
            options.Converters.Add(new HistogramConverter());
            return JsonSerializer.Deserialize<MCMCResults>(bytes, options);
        }

        #endregion

    }
}
