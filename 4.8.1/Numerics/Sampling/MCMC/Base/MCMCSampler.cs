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
using System.Threading;
using System.Threading.Tasks;

namespace Numerics.Sampling.MCMC
{

    /// <summary>
    /// The log-likelihood function to evaluate.
    /// </summary>
    /// <param name="parameters">The list of parameters to evaluate.</param>
    /// <returns>The log-Likelihood given the parameter set.</returns>
    /// <remarks>
    /// This function should account for the data likelihood 
    /// as well as the prior likelihood of the model parameters.
    /// </remarks>
    [Serializable]
    public delegate double LogLikelihood(double[] parameters);

    /// <summary>
    /// A base class for all Markov Chain Monte Carlo (MCMC) samplers.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class MCMCSampler
    {

        /// <summary>
        /// Constructs a new MCMC sampler.
        /// </summary>
        /// <param name="priorDistributions">The list of prior distributions for the model parameters.</param>
        /// <param name="logLikelihoodFunction">The Log-Likelihood function to evaluate.</param>    
        public MCMCSampler(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction)
        {
            PriorDistributions = priorDistributions;
            LogLikelihoodFunction = logLikelihoodFunction;

            // Initialize output arrays
            Reset();
        }


        #region Inputs

        protected int _prngSeed = 12345;
        protected int _initialIterations = 10;
        protected int _warmupIterations = 1750;
        protected int _iterations = 3500;
        protected int _numberOfChains = 4;
        protected int _thinningInterval = 20;

        /// <summary>
        /// Gets and sets the pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed
        {
            get { return _prngSeed; }
            set
            {
                _prngSeed = value;
                Reset();
            }
        }

        /// <summary>
        /// Determines the number of iterations used to initialize the chains. It is recommended that the initial iterations be at least 10 x number of parameters in length.
        /// </summary>

        public int InitialIterations
        {
            get { return _initialIterations; }
            set
            {
                _initialIterations = value;
                Reset();
            }
        }

        /// <summary>
        /// Gets and sets the number of warm up MCMC iterations to discard at the beginning of the simulation.
        /// </summary>
        public int WarmupIterations
        {
            get { return _warmupIterations; }
            set
            {
                _warmupIterations = value;
                Reset();
            }
        }

        /// <summary>
        /// Gets and sets the number of MCMC iterations to simulate.
        /// </summary>
        public int Iterations
        {
            get { return _iterations; }
            set
            {
                _iterations = value;
                Reset();
            }
        }

        /// <summary>
        /// Gets and sets the number of Markov Chains.
        /// </summary>
        public int NumberOfChains
        {
            get { return _numberOfChains; }
            set
            {
                _numberOfChains = value;
                Reset();
            }
        }

        /// <summary>
        /// Gets and sets the thinning interval. This determines how often the MCMC iterations will be recorded and evaluated.
        /// </summary>
        public int ThinningInterval
        {
            get { return _thinningInterval; }
            set
            {
                _thinningInterval = value;
                Reset();
            }
        }

        /// <summary>
        /// The number of simulations that have been run with this instance of the sampler. 
        /// </summary>
        protected int _simulations = 0;

        /// <summary>
        /// The master pseudo random number generator (PRNG).
        /// </summary>
        protected Random _masterPRNG;

        /// <summary>
        /// The PRNG for each Markov Chain.
        /// </summary>
        protected Random[] _chainPRNGs;

        /// <summary>
        /// The current states of each chain. 
        /// </summary>
        protected ParameterSet[] _chainStates;

        /// <summary>
        /// The Log-Likelihood function to evaluate. 
        /// </summary>
        public LogLikelihood LogLikelihoodFunction { get; protected set; }

        /// <summary>
        /// Gets and sets the list of prior distributions for the model parameters.
        /// </summary>
        public List<IUnivariateDistribution> PriorDistributions { get; protected set; }

        /// <summary>
        /// Gets the number of parameters to evaluate.
        /// </summary>
        public int NumberOfParameters => PriorDistributions.Count;

        /// <summary>
        /// Determines whether to update the population matrix when the chain states are recorded.
        /// </summary>
        public bool IsPopulationSampler { get; protected set; } = false;

        /// <summary>
        /// Determines if the chains should be sampled in parallel. Default = true.
        /// </summary>
        public bool ParallelizeChains { get; set; } = true;

        /// <summary>
        /// Determines if the MCMC simulation should be resumed. Default = false.
        /// </summary>
        public bool ResumeSimulation { get; set; } = true;

        /// <summary>
        /// Enumerates the initialization types.
        /// </summary>
        public enum InitializationType
        {
            /// <summary>
            /// Initialize the chains using the Maximum a Posteriori (MAP) estimate and covariance matrix.
            /// If the MAP optimization fails, chains will be automatically initialization with random samples from the priors.
            /// </summary>
            MAP,
            /// <summary>
            /// Automatically initialize the chains with random samples from the priors. This is the default.
            /// </summary>
            Randomize,
            /// <summary>
            /// Initialize the chains from user-defined points. 
            /// </summary>
            UserDefined,
        }

        /// <summary>
        /// Determines whether to initialize the chains using the Maximum a Posteriori (MAP) estimate and covariance matrix.
        /// </summary>
        public InitializationType Initialize { get; set; } = InitializationType.Randomize;

        /// <summary>
        /// Determines if the Maximum a Posteriori (MAP) estimate was successful.
        /// </summary>
        protected bool _mapSuccessful = false;

        /// <summary>
        /// The Multivariate Normal proposal distribution set from the MAP estimate.
        /// </summary>
        protected MultivariateNormal _MVN;

        /// <summary>
        /// Event is raised when the simulation progress changes.
        /// </summary>
        public event ProgressChangedEventHandler ProgressChanged;

        /// <summary>
        /// Event is raised when the simulation progress changes.
        /// </summary>
        /// <param name="percentComplete">The percent complete as decimal between 0 and 1; e.g. 10% complete is passed through the event as 0.1.</param>
        /// <param name="progressText"></param>
        public delegate void ProgressChangedEventHandler(double percentComplete, string progressText);

        /// <summary>
        /// Get and set the progress changed rate. The default is to update progress with every 1% (0.01) change in progress. 
        /// </summary>
        public double ProgressChangedRate { get; set; } = 0.01;

        /// <summary>
        /// Cancellation token source.
        /// </summary>
        public CancellationTokenSource CancellationTokenSource { get; set; }

        #endregion

        #region Outputs

        /// <summary>
        /// Gets the population matrix used for population-based samplers.
        /// </summary>
        public List<ParameterSet> PopulationMatrix { get; protected set; }

        /// <summary>
        /// Gets the list of sampled Markov Chains.
        /// </summary>
        public List<ParameterSet>[] MarkovChains { get; protected set; }

        /// <summary>
        /// Keeps track of the number of accepted samples per chain.
        /// </summary>
        public int[] AcceptCount { get; protected set; }

        /// <summary>
        /// Keeps track of the number of calls to the proposal sampler per chain.
        /// </summary>
        public int[] SampleCount { get; protected set; }

        /// <summary>
        /// The acceptance rate per chain.
        /// </summary>
        public double[] AcceptanceRates
        {
            get
            {
                var ar = new double[NumberOfChains];
                for (int i = 0; i < NumberOfChains; i++)
                    ar[i] = (double)AcceptCount[i] / (double)SampleCount[i];
                return ar;
            }
        }

        /// <summary>
        /// The average log-likelihood across each chain for each iteration.
        /// </summary>
        public List<double>  MeanLogLikelihood { get; protected set; }

        /// <summary>
        /// Gets and sets the number of posterior parameter sets to output.
        /// </summary>
        public int OutputLength { get; set; } = 10000;

        /// <summary>
        /// Output posterior parameter sets. These are recorded after the iterations have been completed. 
        /// </summary>
        public List<ParameterSet>[] Output { get; protected set; }

        /// <summary>
        /// The output parameter set that produced the maximum likelihood. 
        /// This is referred to as the maximum a posteriori (MAP). 
        /// </summary>
        public ParameterSet MAP { get; protected set; }

        #endregion

        #region Simulation Methods

        /// <summary>
        /// Validate the sampler settings.
        /// </summary>
        protected virtual void ValidateSettings()
        {
            if (NumberOfChains < 1) throw new ArgumentException(nameof(NumberOfChains), "There must be at least 1 chain.");
            if (Iterations < 100) throw new ArgumentException(nameof(Iterations), "The number of iterations cannot be less than 100.");
            if (WarmupIterations < 1) throw new ArgumentException(nameof(WarmupIterations), "The number of warm up iterations cannot be less than 1.");
            if (WarmupIterations > (int)(0.5 * Iterations)) throw new ArgumentException(nameof(WarmupIterations), "The number of warm up iterations cannot be greater than half the number of iterations.");
            if (ThinningInterval < 1) throw new ArgumentException(nameof(ThinningInterval), "The thinning interval cannot be less than 1.");
            if (InitialIterations < NumberOfChains) throw new ArgumentException(nameof(InitialIterations), "The initial population cannot be less than the number of chains.");
            if (OutputLength < 100) throw new ArgumentException(nameof(OutputLength), "The output length must be at least 100.");
            ValidateCustomSettings();
        }

        /// <summary>
        /// Validate any custom MCMC sampler settings. 
        /// </summary>
        protected virtual void ValidateCustomSettings() { }

        /// <summary>
        /// Initialize any custom MCMC sampler settings.
        /// </summary>
        protected virtual void InitializeCustomSettings() { }

        /// <summary>
        /// Initialize the Markov Chains.
        /// </summary>
        protected virtual ParameterSet[] InitializeChains()
        {
            if (Initialize == InitializationType.UserDefined)
            {
                // If user-defined, return the last states of the chains.
                var chainStates = new ParameterSet[NumberOfChains];
                for (int i = 0; i < NumberOfChains; i++)
                {
                    chainStates[i] = MarkovChains[i].Last();
                }
                return chainStates;
            }

            var prng = new Random(PRNGSeed);
            var rnds = LatinHypercube.Random(InitialIterations, NumberOfParameters, prng.Next());
            var parameters = new double[NumberOfParameters];
            var tempPopulation = new List<ParameterSet>();       
            var initials = new ParameterSet[NumberOfChains];
            double logLH = 0;

            if (Initialize == InitializationType.MAP)
            {
                // Use differential evolution to find a global optimum
                var lowerBounds = PriorDistributions.Select(x => x.Minimum).ToArray();
                var upperBounds = PriorDistributions.Select(x => x.Maximum).ToArray();
                var inititals = lowerBounds.Add(upperBounds).Divide(2d);
                var DE = new DifferentialEvolution((x) => { return LogLikelihoodFunction(x); }, NumberOfParameters, lowerBounds, upperBounds);
                DE.ReportFailure = false;
                DE.Maximize();
                if (DE.Status == OptimizationStatus.Success)
                {
                    try
                    {
                        _mapSuccessful = true;
                        // Get MAP
                        MAP = DE.BestParameterSet.Clone();
                        // Get Fisher Information Matrix
                        var fisher = DE.Hessian * -1d;
                        // Invert it to get the covariance matrix, and scale to give wider coverage
                        var covar = fisher.Inverse() * 2.5;
                        

                        // Set up proposal distribution
                        _MVN = new MultivariateNormal(MAP.Values, covar.ToArray());
                        // Then randomly sample from the proposal
                        for (int i = 0; i < InitialIterations; i++)
                        {
                            parameters = _MVN.InverseCDF(rnds.GetRow(i));
                            logLH = LogLikelihoodFunction(parameters);
                            if (IsPopulationSampler) PopulationMatrix.Add(new ParameterSet(parameters, logLH));
                            tempPopulation.Add(new ParameterSet(parameters, logLH));
                        }

                        // Set the initial vectors randomly from the MVN proposal
                        for (int i = 0; i < NumberOfChains; i++)
                            initials[i] = tempPopulation[i].Clone();

                        return initials;

                    }
                    catch (Exception) 
                    {
                        // If this fails go to naive initialization below
                        Initialize = InitializationType.Randomize;
                    }
                }
            }


            // *** If not using MAP or if MAP fails, then use random initialization *** //

            // First add the mean of the priors
            for (int j = 0; j < NumberOfParameters; j++)
                parameters[j] = PriorDistributions[j].Mean;
            logLH = LogLikelihoodFunction(parameters);
            if (IsPopulationSampler) PopulationMatrix.Add(new ParameterSet(parameters, logLH));
            tempPopulation.Add(new ParameterSet(parameters, logLH));

            // If the initial population and the number of chains is 1, 
            // then just take the mean of the priors
            if (InitialIterations == 1 && NumberOfChains == 1)
            {
                initials[0] = tempPopulation.First().Clone();
                return initials;
            }

            // Then randomly sample from the priors
            for (int i = 1; i < InitialIterations; i++)
            {
                for (int j = 0; j < NumberOfParameters; j++)
                    parameters[j] = PriorDistributions[j].InverseCDF(rnds[i, j]);
                logLH = LogLikelihoodFunction(parameters);
                if (IsPopulationSampler) PopulationMatrix.Add(new ParameterSet(parameters, logLH));
                tempPopulation.Add(new ParameterSet(parameters, logLH));
            }
            
            // Sort temp population by log-likelihood in descending order
            tempPopulation.Sort((x, y) => -1 * x.Fitness.CompareTo(y.Fitness));

            // Set the initial vectors to the best performing parameter sets
            for (int i = 0; i < NumberOfChains; i++)
                initials[i] = tempPopulation[i].Clone();

            return initials;
        }

        /// <summary>
        /// Sample the Markov Chain. 
        /// </summary>
        /// <param name="index">The Markov Chain zero-based index</param>
        /// <param name="state">The initial state.</param>
        protected virtual ParameterSet SampleChain(int index, ParameterSet state)
        {
            // Sample until thinning interval, then return the last state.
            for (int j = 1; j <= ThinningInterval; j++)
            {
                state = ChainIteration(index, state).Clone();
            }
            return state;
        }

        /// <summary>
        /// Returns a proposed MCMC parameter set and its fitness. 
        /// </summary>
        /// <param name="index">The Markov Chain zero-based index.</param>
        /// <param name="state">The current chain state to compare against.</param>
        protected abstract ParameterSet ChainIteration(int index, ParameterSet state);

        /// <summary>
        /// Sample the Markov Chains.
        /// </summary>
        public virtual void Sample()
        {
            // Validate the input settings
            ValidateSettings();

            CancellationTokenSource = new CancellationTokenSource();

            // Setup the sampler
            if (ResumeSimulation = false || _simulations < 1)
            {
                // Initialize the chains
                _chainStates = InitializeChains();

                // Initialize custom settings
                InitializeCustomSettings();
            }

            // Output settings
            int outputIterations = (int)Math.Ceiling(OutputLength / (double)NumberOfChains);
            int totalIterations = Iterations + outputIterations;
            int outputCount = 0;
            Output = new List<ParameterSet>[NumberOfChains];
            for (int i = 0; i < NumberOfChains; i++)
                Output[i] = new List<ParameterSet>();

            // progress counter
            int progress = 0;

            // Sample chains
            for (int i = 1; i <= totalIterations; i++)
            {

                if (ParallelizeChains)
                {
                    Parallel.For(0, NumberOfChains, (j) => { _chainStates[j] = SampleChain(j, _chainStates[j]); });
                }
                else
                {
                    for (int j = 0; j < NumberOfChains; j++)
                        _chainStates[j] = SampleChain(j, _chainStates[j]);
                }

                // Save chain states
                if (i <= Iterations)
                    MeanLogLikelihood.Add(0);

                // Record output
                for (int j = 0; j < NumberOfChains; j++)
                {
                    // Update population
                    if (IsPopulationSampler == true)
                        PopulationMatrix.Add(_chainStates[j].Clone());

                    if (i <= Iterations)
                    {
                        // Record mean log-likelihood
                        MeanLogLikelihood[i - 1] += _chainStates[j].Fitness / NumberOfChains;

                        // Save chain state
                        MarkovChains[j].Add(_chainStates[j].Clone());
                    }
                    else if (i > Iterations && outputCount < OutputLength)
                    {
                        // Record the output and keep track of MAP
                        Output[j].Add(_chainStates[j].Clone());
                        outputCount++;
                        if (_chainStates[j].Fitness > MAP.Fitness)
                            MAP = _chainStates[j].Clone();
                    }
                }

                // Check for cancellation
                if (CancellationTokenSource.Token.IsCancellationRequested)
                    break;

                // Update progress
                progress += 1;
                if (progress % (int)(totalIterations * ProgressChangedRate) == 0)
                {
                    ReportProgress((double)progress / totalIterations);
                }

            }

            _simulations += 1;
        }

        /// <summary>
        /// Cancel the MCMC simulation.
        /// </summary>
        public void CancelSimulation()
        {
            if (CancellationTokenSource != null && CancellationTokenSource.Token != null  && CancellationTokenSource.Token.CanBeCanceled == true)
            {
                CancellationTokenSource.Cancel();
            }
        }

        /// <summary>
        /// Report the simulation progress.
        /// </summary>
        /// <param name="percentComplete">The percent complete as decimal between 0 and 1; e.g. 10% complete is passed through the event as 0.1.</param>
        public void ReportProgress(double percentComplete)
        {
            ProgressChanged?.Invoke(percentComplete, (percentComplete * 100) + "%");
        }

        /// <summary>
        /// Reset simulation results.
        /// </summary>
        public void Reset()
        {
            _simulations = 0;
            // Clear old memory and re-instantiate the result storage
            _masterPRNG = new Random(PRNGSeed);
            _chainPRNGs = new Random[NumberOfChains];
            PopulationMatrix = new List<ParameterSet>();
            MarkovChains = new List<ParameterSet>[NumberOfChains];
            for (int i = 0; i < NumberOfChains; i++)
            {
                _chainPRNGs[i] = new Random(_masterPRNG.Next());
                MarkovChains[i] = new List<ParameterSet>();
            }
            AcceptCount = new int[NumberOfChains];
            SampleCount = new int[NumberOfChains];
            MeanLogLikelihood = new List<double>();
            MAP = new ParameterSet([], double.MinValue);
            Output = new List<ParameterSet>[NumberOfChains];
        }

        #endregion

    }
}
