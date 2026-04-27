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
