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
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Sampling.MCMC
{

    /// <summary>
    /// This class performs Bayesian MCMC using the random walk Metropolis-Hastings (RWMH) method.
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
    public class RWMH : MCMCSampler
    {

        /// <summary>
        /// Constructs a new RWMH sampler.
        /// </summary>
        /// <param name="priorDistributions">The list of prior distributions for the model parameters.</param>
        /// <param name="logLikelihoodFunction">The Log-Likelihood function to evaluate.</param>
        /// <param name="proposalSigma">The covariance matrix Σ (sigma) for the proposal distribution.</param>
        public RWMH(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction, Matrix proposalSigma) : base(priorDistributions, logLikelihoodFunction)
        {
            InitialIterations = 100 * NumberOfParameters;
            ProposalSigma = proposalSigma;
        }

        private MultivariateNormal[] mvn = null!;

        /// <summary>
        /// The covariance matrix Σ (sigma) for the proposal distribution.
        /// </summary>
        public Matrix ProposalSigma { get; private set; }

        /// <inheritdoc/>
        protected override void ValidateCustomSettings()
        {
            if (ProposalSigma == null) throw new ArgumentException(nameof(ProposalSigma), "The proposal covariance matrix cannot be null.");
            if (ProposalSigma.NumberOfRows != ProposalSigma.NumberOfColumns) throw new ArgumentException(nameof(ProposalSigma), "The proposal covariance matrix must be square.");
            if (ProposalSigma.NumberOfRows != NumberOfParameters) throw new ArgumentException(nameof(ProposalSigma), "The proposal covariance matrix must have the same number of rows and columns as the number of parameters.");
        }

        /// <inheritdoc/>
        protected override void InitializeCustomSettings()
        {
            // Set up multivariate Normal distributions for each chain
            mvn = new MultivariateNormal[NumberOfChains];
            for (int i = 0; i < NumberOfChains; i++)
            {
                mvn[i] = new MultivariateNormal(NumberOfParameters);
            }
            // Set up proposal matrix
            if (Initialize == InitializationType.MAP && _mapSuccessful)
            {
                ProposalSigma = new Matrix(base._MVN.Covariance);
            }
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {
            // Update the sample count
            SampleCount[index] += 1;

            // Get proposal vector
            mvn[index].SetParameters(state.Values, ProposalSigma.Array);
            var xp = mvn[index].InverseCDF(_chainPRNGs[index].NextDoubles(NumberOfParameters));

            for (int i = 0; i < NumberOfParameters; i++)
            {
                // Check if the parameter is feasible (within the constraints)
                if (xp[i] < PriorDistributions[i].Minimum || xp[i] > PriorDistributions[i].Maximum)
                {
                    // The proposed parameter vector was infeasible, 
                    // so leave xi unchanged.
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
                return new ParameterSet(xp, logLHp);
            }
            else
            {
                return state;
            }
        }


    }
}
