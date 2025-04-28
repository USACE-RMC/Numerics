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

using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace Numerics.Sampling.MCMC
{
    /// <summary>
    /// This class performs Bayesian inference using the self-normalizing importance sampling (SNIS) method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class SNIS : MCMCSampler
    {
        /// <summary>
        /// Constructs a self-normalizing importance sampler.
        /// </summary>
        /// <param name="priorDistributions">The list of prior distributions for the model parameters.</param>
        /// <param name="logLikelihoodFunction">The Log-Likelihood function to evaluate.</param>
        /// <param name="multivariateNormal">Optional. The multivariate Normal distribution is used for importance sampling. If null, naive Monte Carlo is performed.</param>
        public SNIS(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction, MultivariateNormal multivariateNormal = null) : base(priorDistributions, logLikelihoodFunction)
        {
            mvn = multivariateNormal;
            useImportanceSampling = multivariateNormal != null ? true : false;

            // Create default settings
            _numberOfChains = 1;
            _warmupIterations = 0;
            _thinningInterval = 1;
            _initialIterations = 1;
            _iterations = 100000;
            OutputLength = 10000;
            Reset();
        }

        private bool useImportanceSampling = false;
        private MultivariateNormal mvn = null;

        /// <inheritdoc/>
        protected override void InitializeCustomSettings()
        {
            if (mvn == null && useImportanceSampling == false && Initialize == InitializationType.MAP && _mapSuccessful)
            {
                mvn = (MultivariateNormal)_MVN.Clone();
                var covar = new Matrix(mvn.Covariance);
                covar = covar * 2.5; // Inflate to cover the posterior
                mvn.SetParameters(mvn.Mean, covar.ToArray());
                useImportanceSampling = true;
            }
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {
            return new ParameterSet();
        }

        /// <inheritdoc/>
        protected override void ValidateSettings()
        {
            if (NumberOfChains != 1) throw new ArgumentException(nameof(InitialIterations), "There can only be 1 chain with this method.");
            if (OutputLength < 100) throw new ArgumentException(nameof(OutputLength), "The output length must be at least 100.");
            if (Iterations < OutputLength) throw new ArgumentException(nameof(Iterations), "The number of iterations cannot be less than the output length.");
            if (WarmupIterations != 0) throw new ArgumentException(nameof(WarmupIterations), "There are no warmup iterations with this method.");
            if (ThinningInterval != 1) throw new ArgumentException(nameof(ThinningInterval), "The thinning interval must be 1 for this method.");
            if (InitialIterations != 1) throw new ArgumentException(nameof(InitialIterations), "The initial population must be 1 for this method.");
            if (mvn != null && mvn.ParametersValid == false)
                throw new ArgumentException(nameof(MultivariateNormal), "The multivariate Normal importance distribution is invalid.");
        }

        /// <summary>
        /// Perform importance sampling.
        /// </summary>
        public override void Sample()
        {
            InitializeChains();
            InitializeCustomSettings();
            ValidateSettings();
            ParallelizeChains = true;
            ResumeSimulation = false;
            CancellationTokenSource = new CancellationTokenSource();

            // Create inputs 
            _masterPRNG = new MersenneTwister(PRNGSeed);
            MarkovChains = new List<ParameterSet>[NumberOfChains];
            var sets = new ParameterSet[Iterations];
            MarkovChains[0] = sets.ToList();
            Output = new List<ParameterSet>[NumberOfChains];
            Output[0] = new List<ParameterSet>();

            // Create sample & accept counts
            AcceptCount = new int[NumberOfChains];
            SampleCount = new int[NumberOfChains];

            // Create mean log-likelihood list
            MeanLogLikelihood = new List<double>();

            // Keeps track of best parameter set
            MAP = new ParameterSet([], double.MinValue, double.MinValue);

            // Create parameter random values
            var rnds = _masterPRNG.NextDoubles(Iterations, NumberOfParameters);

            // Perform sampling
            Parallel.For(0, Iterations, (idx) =>
            {
                // Sample parameters
                var parameters = new double[NumberOfParameters];
                double logLH = 0;
                double weight = 0;
                if (useImportanceSampling == true)
                {
                    parameters = mvn.InverseCDF(rnds.GetRow(idx));
                    logLH = LogLikelihoodFunction(parameters);
                    weight = logLH - mvn.LogPDF(parameters);
                }
                else
                {
                    for (int i = 0; i < NumberOfParameters; i++)
                        parameters[i] = PriorDistributions[i].InverseCDF(rnds[idx, i]);
                    logLH = LogLikelihoodFunction(parameters);
                    weight = logLH;
                }
                // Record sample
                MarkovChains[0][idx] = new ParameterSet(parameters, logLH, weight);
            });
       
            // Get the maximum a posteriori
            for (int i = 0; i < Iterations; i++)
                if (MarkovChains[0][i].Weight > MAP.Weight)
                    MAP = MarkovChains[0][i].Clone();

            // Get the normalization factor
            double max = MAP.Weight;
            double sum = 0;
            for (int i = 0; i < Iterations; i++)
            {
                if (MarkovChains[0][i].Weight != double.MinValue)
                {
                    sum += Math.Exp(MarkovChains[0][i].Weight - max);
                }            
            }               
            double normalization = max + Math.Log(sum);

            // Compute the posterior weights
            Parallel.For(0, Iterations, (idx) => 
            {
                double w = MarkovChains[0][idx].Weight != double.MinValue ? Math.Exp(MarkovChains[0][idx].Weight - normalization) : 0d;
                MarkovChains[0][idx] = new ParameterSet(MarkovChains[0][idx].Values, MarkovChains[0][idx].Fitness, w);
            });

            // Sort list in ascending order of posterior weights
            MarkovChains[0].Sort((x, y) => x.Fitness.CompareTo(y.Fitness));
            var cdf = new double[Iterations];
            cdf[0] = MarkovChains[0][0].Weight;

            // Create CDF
            for (int i = 1; i < Iterations; i++)
                cdf[i] = cdf[i-1] + MarkovChains[0][i].Weight;

            // Record output
            var rndOut  = PlottingPositions.Weibull(OutputLength);
            int idx = 0;
            for (int i = 0; i < OutputLength; i++)
            {
                idx = Search.Sequential(rndOut[i], cdf, idx);
                Output[0].Add(MarkovChains[0][idx].Clone());
            }
        }

    }
}
