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
    /// <para>
    /// <see cref="MCMCSampler.LogLikelihoodFunction"/> returns the full log-posterior, as for the MCMC
    /// samplers, in both modes. Each draw is weighted by the delegate value minus the log-density of the
    /// distribution it was drawn from: the joint prior, the sum of the
    /// <see cref="MCMCSampler.PriorDistributions"/> log-densities, when no importance distribution is
    /// given, or the importance distribution q, supplied to the constructor or built from the MAP estimate
    /// when <see cref="MCMCSampler.Initialize"/> is <see cref="MCMCSampler.InitializationType.MAP"/>.
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
        public SNIS(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction, MultivariateNormal multivariateNormal = null!) : base(priorDistributions, logLikelihoodFunction)
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
        private MultivariateNormal mvn = null!;

        /// <inheritdoc/>
        protected override void InitializeCustomSettings()
        {
            if (mvn == null && useImportanceSampling == false && Initialize == InitializationType.MAP && _mapSuccessful && _MVN != null)
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
            if (NumberOfChains != 1) throw new ArgumentException("There can only be 1 chain with this method.", nameof(NumberOfChains));
            if (OutputLength < 100) throw new ArgumentException("The output length must be at least 100.", nameof(OutputLength));
            if (Iterations < OutputLength) throw new ArgumentException("The number of iterations cannot be less than the output length.", nameof(Iterations));
            if (WarmupIterations != 0) throw new ArgumentException("There are no warmup iterations with this method.", nameof(WarmupIterations));
            if (ThinningInterval != 1) throw new ArgumentException("The thinning interval must be 1 for this method.", nameof(ThinningInterval));
            if (InitialIterations != 1) throw new ArgumentException("The initial population must be 1 for this method.", nameof(InitialIterations));
            if (mvn != null && mvn.ParametersValid == false)
                throw new ArgumentException("The multivariate Normal importance distribution is invalid.", nameof(MultivariateNormal));
        }

        /// <summary>
        /// Perform importance sampling.
        /// </summary>
        /// <remarks>
        /// <para>
        /// After the normalized posterior weights are computed, the sample list is sorted ascending on
        /// <c>ParameterSet.Weight</c>, the quantity the resampling CDF then accumulates, and the output is
        /// drawn from that CDF at the Weibull plotting positions k / (OutputLength + 1): each position
        /// selects the first draw whose cumulative weight reaches it, so every draw is output in proportion
        /// to its own weight. <see cref="MCMCSampler.MAP"/> is the draw
        /// with the largest finite <c>Fitness</c>, the largest log-posterior.
        /// </para>
        /// <para>
        /// The sort is stable, so draws with tied weights (commonly many zero weights from -Infinity
        /// log-posteriors under wide priors) keep their original draw order and a seeded run resamples the
        /// same output on every run and platform.
        /// </para>
        /// <para>
        /// The draws are evaluated with <see cref="Parallel.For(int, int, Action{int})"/> when
        /// <see cref="MCMCSampler.ParallelizeChains"/> is true (the default), so the log-likelihood is then called
        /// concurrently and must be thread safe; when false they are evaluated in order on the calling thread. Every
        /// draw comes from a matrix generated up front from <see cref="MCMCSampler.PRNGSeed"/> and is stored at its own
        /// index, so both settings give bit-identical output.
        /// </para>
        /// </remarks>
        public override void Sample()
        {
            // Validate before initializing, because MAP initialization runs a full optimization.
            ValidateSettings();
            InitializeChains();
            InitializeCustomSettings();
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
            MAP = new ParameterSet([], double.NegativeInfinity, double.NegativeInfinity);

            // Create parameter random values
            var rnds = _masterPRNG.NextDoubles(Iterations, NumberOfParameters);

            // Perform sampling
            ForEachDraw(idx =>
            {
                // Sample parameters
                var parameters = new double[NumberOfParameters];
                double logLH = 0;
                double weight = 0;
                if (useImportanceSampling == true && mvn != null)
                {
                    parameters = mvn.InverseCDF(rnds.GetRow(idx));
                    logLH = LogLikelihoodFunction(parameters);
                    weight = logLH - mvn.LogPDF(parameters);
                }
                else
                {
                    // The draw comes from the joint prior, so its log-density is the importance correction.
                    double logPrior = 0;
                    for (int i = 0; i < NumberOfParameters; i++)
                    {
                        parameters[i] = PriorDistributions[i].InverseCDF(rnds[idx, i]);
                        logPrior += PriorDistributions[i].LogPDF(parameters[i]);
                    }
                    logLH = LogLikelihoodFunction(parameters);
                    weight = logLH - logPrior;
                }
                // Record sample
                MarkovChains[0][idx] = new ParameterSet(parameters, logLH, weight);
            });
       
            // Get the maximum a posteriori, the draw with the largest log-posterior
            for (int i = 0; i < Iterations; i++)
                if (Tools.IsFinite(MarkovChains[0][i].Fitness) && MarkovChains[0][i].Fitness > MAP.Fitness)
                    MAP = MarkovChains[0][i].Clone();

            // Get the normalization factor from the largest finite log-weight
            double max = double.NegativeInfinity;
            for (int i = 0; i < Iterations; i++)
                if (Tools.IsFinite(MarkovChains[0][i].Weight) && MarkovChains[0][i].Weight > max)
                    max = MarkovChains[0][i].Weight;
            if (!Tools.IsFinite(max))
                throw new InvalidOperationException("SNIS failed because all importance weights are non-finite.");

            double sum = 0;
            for (int i = 0; i < Iterations; i++)
            {
                if (Tools.IsFinite(MarkovChains[0][i].Weight))
                {
                    sum += Math.Exp(MarkovChains[0][i].Weight - max);
                }            
            }
            if (!Tools.IsFinite(sum) || sum <= 0d)
                throw new InvalidOperationException("SNIS failed because the finite importance weights could not be normalized.");

            double normalization = max + Math.Log(sum);
            if (!Tools.IsFinite(normalization))
                throw new InvalidOperationException("SNIS failed because the importance-weight normalization is non-finite.");

            // Compute the posterior weights
            ForEachDraw(idx =>
            {
                double w = Tools.IsFinite(MarkovChains[0][idx].Weight) ? Math.Exp(MarkovChains[0][idx].Weight - normalization) : 0d;
                MarkovChains[0][idx] = new ParameterSet(MarkovChains[0][idx].Values, MarkovChains[0][idx].Fitness, w);
            });

            // Sort ascending on the normalized weight the CDF below accumulates. OrderBy is a stable sort,
            // so tied weights (commonly many zero weights from -Infinity draws under wide priors) keep their
            // draw order; List<T>.Sort's tie order is implementation-defined and differs between target
            // frameworks, so a seeded run would resample differently per runtime.
            MarkovChains[0] = MarkovChains[0].OrderBy(x => x.Weight).ToList();
            var cdf = new double[Iterations];
            cdf[0] = Math.Max(0.0, MarkovChains[0][0].Weight);

            // Create CDF
            for (int i = 1; i < Iterations; i++)
                cdf[i] = cdf[i-1] + Math.Max(0.0, MarkovChains[0][i].Weight);

            // Record output. Each plotting position p selects the first draw whose cumulative weight
            // reaches p, the inverse of the weight CDF, so draw j is selected for p in (cdf[j-1], cdf[j]]
            // and is output in proportion to its own weight. The positions ascend, so the scan only moves
            // forward; it stops at the last draw in case rounding leaves cdf[^1] just below a position.
            var rndOut  = PlottingPositions.Weibull(OutputLength);
            int idx = 0;
            for (int i = 0; i < OutputLength; i++)
            {
                while (idx < Iterations - 1 && cdf[idx] < rndOut[i])
                    idx++;
                Output[0].Add(MarkovChains[0][idx].Clone());
            }
        }

        /// <summary>
        /// Runs a loop over the draws, in parallel or in index order as <see cref="MCMCSampler.ParallelizeChains"/> selects.
        /// </summary>
        /// <param name="body">The body, called once for each draw index from 0 to <see cref="MCMCSampler.Iterations"/> - 1.</param>
        /// <remarks>
        /// The loop does not use the machine-learning range helper, which chooses serial or parallel execution from an
        /// arithmetic work estimate: a draw's cost is that of the caller's log-likelihood, and
        /// <see cref="MCMCSampler.ParallelizeChains"/> selects the mode explicitly. Each draw writes only its own index,
        /// so either schedule gives the same results.
        /// </remarks>
        private void ForEachDraw(Action<int> body)
        {
            if (ParallelizeChains)
            {
                Parallel.For(0, Iterations, body);
                return;
            }
            for (int idx = 0; idx < Iterations; idx++)
                body(idx);
        }

    }
}
