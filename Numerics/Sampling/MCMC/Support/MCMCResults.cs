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

using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.IO;
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
        public List<ParameterSet>[] MarkovChains { get; private set; }

        /// <summary>
        /// Output posterior parameter sets.
        /// </summary>
        [JsonInclude]
        public List<ParameterSet> Output { get; private set; }

        /// <summary>
        /// The average log-likelihood across each chain for each iteration.
        /// </summary>
        [JsonInclude]
        public List<double> MeanLogLikelihood { get; private set; }

        /// <summary>
        /// The acceptance rate for each chain.
        /// </summary>
        [JsonInclude]
        public double[] AcceptanceRates { get; private set; }

        /// <summary>
        /// Parameter results using the output posterior parameter sets.
        /// </summary>
        [JsonInclude]
        public ParameterResults[] ParameterResults { get; private set; }

        /// <summary>
        /// The output parameter set that produced the maximum likelihood.
        /// This is referred to as the maximum a posteriori (MAP).
        /// </summary>
        [JsonInclude]
        public ParameterSet MAP { get; private set; }

        /// <summary>
        /// The mean of the posterior distribution of each parameter.
        /// </summary>
        [JsonInclude]
        public ParameterSet PosteriorMean { get; private set; }

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
                IncludeFields = true
            };
            return JsonSerializer.SerializeToUtf8Bytes(mcmcResults, options);
        }

        /// <summary>
        /// Creates MCMC Results from a byte array.
        /// </summary>
        /// <param name="bytes">Byte array.</param>
        public static MCMCResults FromByteArray(byte[] bytes)
        {
            var options = new JsonSerializerOptions
            {
                DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
                IncludeFields = true
            };
            try
            {
                return JsonSerializer.Deserialize<MCMCResults>(bytes, options);
            }
            catch
            {
                ///Previous serialization used Binary Formatter, which won't deserialize cleanly as JSON. 
                /// If this fails, then it's probably the bf bytes. fall back to legacy.
                return FromByteArrayLegacy(bytes);
            }
        }

        /// <summary>
        /// Creates MCMC Results from a byte array.
        /// </summary>
        /// <param name="bytes">Byte array.</param>
        private static MCMCResults FromByteArrayLegacy(byte[] bytes)
        {
            using var ms = new MemoryStream();
            #pragma warning disable SYSLIB0011 // Suppress obsolete BinaryFormatter warning for legacy support
            var bf = new System.Runtime.Serialization.Formatters.Binary.BinaryFormatter();
            ms.Write(bytes, 0, bytes.Length);
            ms.Seek(0L, SeekOrigin.Begin);
            var obj = bf.Deserialize(ms);
            #pragma warning disable SYSLIB0011 // Suppress obsolete BinaryFormatter warning for legacy support
            return (MCMCResults)obj;
        }

        #endregion

    }
}
