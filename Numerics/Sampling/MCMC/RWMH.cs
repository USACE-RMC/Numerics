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
            if (ProposalSigma == null) throw new ArgumentException("The proposal covariance matrix cannot be null.", nameof(ProposalSigma));
            if (ProposalSigma.NumberOfRows != ProposalSigma.NumberOfColumns) throw new ArgumentException("The proposal covariance matrix must be square.", nameof(ProposalSigma));
            if (ProposalSigma.NumberOfRows != NumberOfParameters) throw new ArgumentException("The proposal covariance matrix must have the same number of rows and columns as the number of parameters.", nameof(ProposalSigma));
        }

        /// <inheritdoc/>
        protected override void InitializeCustomSettings()
        {
            // Set up proposal matrix
            if (Initialize == InitializationType.MAP && _mapSuccessful && _MVN != null)
            {
                ProposalSigma = new Matrix(_MVN.Covariance);
            }
            // Set up multivariate Normal distributions for each chain. The proposal covariance is
            // fixed for the whole run, so each chain's proposal is factorized exactly once here and
            // ChainIteration then translates only the mean, keeping the factorization. A covariance
            // that fails the factorization therefore throws here, before sampling starts, instead of
            // from inside the first chain iteration.
            mvn = new MultivariateNormal[NumberOfChains];
            for (int i = 0; i < NumberOfChains; i++)
            {
                mvn[i] = new MultivariateNormal(NumberOfParameters);
                mvn[i].SetParameters(new double[NumberOfParameters], ProposalSigma.Array);
            }
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {
            // Update the sample count
            SampleCount[index] += 1;

            // Get proposal vector. The proposal covariance was factorized once at initialization, so
            // only the mean moves with the chain state — a translation changes nothing the
            // factorization derives from the covariance.
            mvn[index].SetMean(state.Values);
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
