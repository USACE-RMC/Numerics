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
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;

namespace Numerics.Sampling.MCMC
{
    /// <summary>
    /// This class performs Bayesian MCMC using the Gibbs method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b>Authors:</b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References:</b>
    /// </para>
    /// <para> 
    ///    <see href="https://en.wikipedia.org/wiki/Gibbs_sampling"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public class Gibbs : MCMCSampler
    {

        /// <summary>
        /// The proposal function for creating a proposal vector of parameters to evaluate.
        /// </summary>
        /// <param name="parameters">The list of parameters to evaluate.</param>
        /// <param name="prng">Random number generator.</param>
        /// <returns>Returns a proposal vector.</returns>
        public delegate double[] Proposal(double[] parameters, Random prng);

        /// <summary>
        /// Constructs a new Gibbs sampler.
        /// </summary>
        /// <param name="priorDistributions">The list of prior distributions for the model parameters.</param>
        /// <param name="logLikelihoodFunction">The Log-Likelihood function to evaluate.</param>
        /// <param name="proposalFunction">The conditional proposal function.</param>
        public Gibbs(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction, Proposal proposalFunction) : base(priorDistributions, logLikelihoodFunction)
        {

            ProposalFunction = proposalFunction;
            // Create default settings
            _numberOfChains = 1;
            _warmupIterations = 1;
            _thinningInterval = 1;
            _initialIterations = 1;
            _iterations = 100000;
            OutputLength = 10000;
            Reset();
        }

        /// <summary>
        /// The proposal function for creating a proposal vector of parameters to evaluate.
        /// </summary>
        public Proposal ProposalFunction { get; }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {
            var xp = ProposalFunction(state.Values, _chainPRNGs[index]);
            return new ParameterSet(xp, LogLikelihoodFunction(xp));
        }
    }
}
