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

using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Data.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Sampling.MCMC
{
    /// <summary>
    /// This class performs Bayesian MCMC using the adaptive random walk Metropolis-Hastings (RWMH) method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    ///    <see href="https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public class ARWMH : MCMCSampler
    {

        /// <summary>
        /// Constructs a new ARWMH sampler.
        /// </summary>
        /// <param name="priorDistributions">The list of prior distributions for the model parameters.</param>
        /// <param name="logLikelihoodFunction">The Log-Likelihood function to evaluate.</param>      
        public ARWMH(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction) : base(priorDistributions, logLikelihoodFunction)
        {
            // Set initial
            InitialIterations = 100 * NumberOfParameters;

            // The optimal scale & adaptive covariance matrix
            Scale = 2.38 * 2.38 / NumberOfParameters;
            // The initial scale & identity covariance matrix
            sigmaIdentity = Matrix.Identity(NumberOfParameters) * (0.1 * 0.1 / NumberOfParameters);
            Beta = 0.05;
        }


        private Matrix sigmaIdentity;
        private RunningCovarianceMatrix[] sigma;
        private MultivariateNormal[] mvn;

        /// <summary>
        /// The scaling parameter used to scale the adaptive covariance matrix.
        /// </summary>
        public double Scale { get; set; }

        /// <summary>
        /// Determines how often to sample from the small identity covariance matrix; e.g., 0.05 will result in sampling 5% of the time.
        /// </summary>
        public double Beta { get; set; }

        /// <summary>
        /// The covariance matrix Σ (sigma) for the proposal distribution.
        /// </summary>
        public Matrix[] ProposalSigma
        {
            get
            {
                var sigmas = new Matrix[NumberOfChains];
                for (int i = 0; i < NumberOfChains; i++)
                {
                    if (sigma[i].N < 2)
                        sigmas[i] = sigmaIdentity;
                    else
                        sigmas[i] = sigma[i].Covariance * (1d / (sigma[i].N - 1));
                }
                return sigmas;
            }
        }

        /// <inheritdoc/>
        protected override void ValidateCustomSettings()
        {
            if (Scale <= 0) throw new ArgumentException(nameof(Scale), "The scale parameter must greater than 0.");
            if (Beta < 0 || Beta > 1) throw new ArgumentException(nameof(Beta), "Beta must be between 0 and 1.");
        }

        /// <inheritdoc/>
        protected override void InitializeCustomSettings()
        {
            // Set up multivariate Normal distributions and 
            // adaptive covariance matrix for each chain
            mvn = new MultivariateNormal[NumberOfChains];
            sigma = new RunningCovarianceMatrix[NumberOfChains];
            for (int i = 0; i < NumberOfChains; i++)
            {
                mvn[i] = new MultivariateNormal(NumberOfParameters);
                sigma[i] = new RunningCovarianceMatrix(NumberOfParameters);

                if (Initialize == InitializationType.MAP && _mapSuccessful)
                {
                    // Hot start the covariance matrix
                    for (int j = 0; j < NumberOfParameters; j++)
                    {
                        for (int k = 0; k < NumberOfParameters; k++)
                        {
                            sigma[i].Covariance[j, k] = _MVN.Covariance[j, k];
                        }
                    }
                }
            }
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {

            // Update the sample count
            SampleCount[index] += 1;

            if (_chainPRNGs[index].NextDouble() <= Beta || SampleCount[index] <= 100 * NumberOfParameters)
            {
                // Use the identity matrix the first 100*D samples
                mvn[index].SetParameters(state.Values, sigmaIdentity.Array);
            }
            else
            {
                // Use the adaptive covariance matrix
                var proposalCov = sigma[index].Covariance * (Scale / (sigma[index].N - 1));
                mvn[index].SetParameters(state.Values, proposalCov.Array);
            }

            // Get proposal vector
            var xp = mvn[index].InverseCDF(_chainPRNGs[index].NextDoubles(NumberOfParameters));

            // Check if the parameter is feasible (within the constraints)
            for (int i = 0; i < NumberOfParameters; i++)
            {
                if (xp[i] < PriorDistributions[i].Minimum || xp[i] > PriorDistributions[i].Maximum)
                {
                    // The proposed parameter vector was infeasible, so leave xi unchanged.
                    // Adapt Covariance Matrix after warmup
                    if (SampleCount[index] > ThinningInterval * WarmupIterations)
                        sigma[index].Push(state.Values);
                    return state;
                }
            }

            // Evaluate fitness
            var logLHp = LogLikelihoodFunction(xp);
            var logLHi = state.Fitness;

            // Calculate the Metropolis ratio
            var logRatio = logLHp - logLHi;

            // Accept the proposal with probability min(1,r)
            // otherwise leave xi unchanged
            var logU = Math.Log(_chainPRNGs[index].NextDouble());
            if (logU <= logRatio)
            {
                // The proposal is accepted
                AcceptCount[index] += 1;
                // Adapt Covariance Matrix
                sigma[index].Push(xp);
                return new ParameterSet(xp, logLHp);
            }
            else
            {
                // Adapt Covariance Matrix after warmup
                if (SampleCount[index] > ThinningInterval * WarmupIterations)
                    sigma[index].Push(state.Values);
                return state;
            }
        }

    }
}
