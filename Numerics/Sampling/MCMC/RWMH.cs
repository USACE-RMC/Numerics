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
            if (Initialize == InitializationType.MAP && _mapSuccessful && _MVN != null)
            {
                ProposalSigma = new Matrix(_MVN.Covariance);
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
