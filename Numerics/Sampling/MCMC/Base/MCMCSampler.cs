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

        /// <summary>
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        protected int _prngSeed = 12345;
        /// <summary>
        /// The number of initial iterations before adaptation begins.
        /// </summary>
        protected int _initialIterations = 10;
        /// <summary>
        /// The number of warmup (burn-in) iterations to discard.
        /// </summary>
        protected int _warmupIterations = 1750;
        /// <summary>
        /// The number of sampling iterations after warmup.
        /// </summary>
        protected int _iterations = 3500;
        /// <summary>
        /// The number of parallel Markov chains to run.
        /// </summary>
        protected int _numberOfChains = 4;
        /// <summary>
        /// The thinning interval for reducing autocorrelation.
        /// </summary>
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
        /// The number of chain transitions a single Markov chain performs during a call to <see cref="Sample"/>.
        /// </summary>
        /// <returns>
        /// (<see cref="Iterations"/> + ceil(<see cref="OutputLength"/> / <see cref="NumberOfChains"/>))
        /// × <see cref="ThinningInterval"/>.
        /// </returns>
        /// <remarks>
        /// <para>
        /// <b>Why this is not <see cref="Iterations"/>.</b> The configured iteration count does not describe
        /// the work a run actually does, and it under-reports it substantially at the defaults. Two
        /// multipliers sit between the two numbers. First, <see cref="Sample"/> runs
        /// ceil(<see cref="OutputLength"/> / <see cref="NumberOfChains"/>) recorded iterations beyond
        /// <see cref="Iterations"/> in order to collect the posterior output. Second, every recorded
        /// iteration advances the chain <see cref="ThinningInterval"/> times, because thinning is applied by
        /// discarding intermediate transitions rather than by discarding recorded draws. Each of those
        /// transitions costs at least one evaluation of <see cref="LogLikelihoodFunction"/>, so this
        /// property, not <see cref="Iterations"/>, is what a runtime estimate should be built on.
        /// </para>
        /// <para>
        /// <b>Warmup is already included.</b> <see cref="WarmupIterations"/> is a subset of
        /// <see cref="Iterations"/>, not an addition to it — <c>ValidateSettings</c> rejects a warmup longer
        /// than half of <see cref="Iterations"/> — so it must not be added again when reasoning about total
        /// cost.
        /// </para>
        /// <para>
        /// <b>At the defaults</b> (<see cref="Iterations"/> = 3,500, <see cref="OutputLength"/> = 10,000,
        /// <see cref="NumberOfChains"/> = 4, <see cref="ThinningInterval"/> = 20) this is
        /// (3,500 + 2,500) × 20 = 120,000 transitions per chain, against a user-facing
        /// <see cref="Iterations"/> that reads 3,500 — a factor of roughly 34.
        /// </para>
        /// <para>
        /// The type is <see cref="long"/> because the product overflows <see cref="int"/> well inside the
        /// range of settings the sampler accepts: nothing bounds <see cref="Iterations"/> from above, and at
        /// the default output length, chain count and thinning interval the per-chain product passes
        /// <see cref="int.MaxValue"/> at about 1.074E8 iterations. <see cref="TotalTransitionCount"/> passes
        /// it four times sooner, at about 2.68E7 iterations, because of the multiplication by the default
        /// four chains.
        /// </para>
        /// <para>
        /// This reports the settings as they are currently configured and does not validate them; the
        /// settings are checked by <c>ValidateSettings</c> when <see cref="Sample"/> is called. The one
        /// exception is <see cref="NumberOfChains"/>: a value below one is not a samplable configuration and
        /// would otherwise divide by zero here, so it reports 0 rather than a framework-dependent number.
        /// </para>
        /// <para>
        /// The count describes the base <see cref="Sample"/> loop and the base <c>SampleChain</c>, neither of
        /// which any chain sampler in this library overrides. It does <b>not</b> describe
        /// <see cref="SNIS"/>, which replaces <see cref="Sample"/> with a single non-Markovian importance
        /// sampling pass and never advances a chain.
        /// </para>
        /// </remarks>
        public long TransitionCount
        {
            get
            {
                // NumberOfChains is only validated by ValidateSettings when Sample() runs, but this property
                // is readable at any time. At zero chains the division below would be Infinity, and casting
                // Infinity to int saturates to int.MaxValue on .NET Core while being unspecified on .NET
                // Framework, so the answer would differ by target framework. Report 0 for a configuration
                // that cannot be sampled instead.
                if (NumberOfChains < 1) return 0L;

                // OutputIterations is the same member Sample() uses, so the two cannot drift apart. Each
                // recorded iteration advances the chain ThinningInterval times. The sum is widened to long
                // before the multiplication so that large settings do not overflow.
                return ((long)Iterations + OutputIterations) * ThinningInterval;
            }
        }

        /// <summary>
        /// The number of recorded iterations that <see cref="Sample"/> runs beyond <see cref="Iterations"/>
        /// in order to collect the posterior output, ceil(<see cref="OutputLength"/> / <see cref="NumberOfChains"/>).
        /// </summary>
        /// <remarks>
        /// Shared by <see cref="Sample"/> and <see cref="TransitionCount"/> so that the reported work and the
        /// work actually performed are computed from a single expression rather than two copies of it.
        /// Callers must ensure <see cref="NumberOfChains"/> is at least one; <see cref="Sample"/> does so via
        /// <c>ValidateSettings</c> and <see cref="TransitionCount"/> guards it directly.
        /// </remarks>
        private int OutputIterations => (int)Math.Ceiling(OutputLength / (double)NumberOfChains);

        /// <summary>
        /// The number of chain transitions performed across all chains during a call to <see cref="Sample"/>.
        /// </summary>
        /// <returns><see cref="TransitionCount"/> × <see cref="NumberOfChains"/>.</returns>
        /// <remarks>
        /// At the defaults this is 120,000 × 4 = 480,000 transitions. Treat it as a <b>lower bound</b> on the
        /// evaluation count, not as a budget: every transition costs at least one evaluation of
        /// <see cref="LogLikelihoodFunction"/>, but a gradient-based sampler such as HMC or NUTS spends many
        /// likelihood and gradient evaluations per transition, and chain initialization adds further
        /// evaluations on top of all of these. See
        /// <see cref="TransitionCount"/> for why this differs so widely from <see cref="Iterations"/>. When
        /// <see cref="ParallelizeChains"/> is true these are distributed across worker threads, so this is
        /// the total work rather than the critical path.
        /// </remarks>
        public long TotalTransitionCount => TransitionCount * NumberOfChains;

        /// <summary>
        /// The number of simulations that have been run with this instance of the sampler.
        /// </summary>
        protected int _simulations = 0;

        /// <summary>
        /// The master pseudo random number generator (PRNG).
        /// </summary>
        protected Random _masterPRNG = null!;

        /// <summary>
        /// The PRNG for each Markov Chain.
        /// </summary>
        protected Random[] _chainPRNGs = null!;

        /// <summary>
        /// The current states of each chain.
        /// </summary>
        protected ParameterSet[] _chainStates = null!;

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
        /// <remarks>
        /// <para>
        /// <b>Thread-safety requirement.</b> When this is true — which is the default —
        /// <see cref="Sample"/> advances all <see cref="NumberOfChains"/> chains inside a
        /// <see cref="System.Threading.Tasks.Parallel.For(int, int, Action{int})"/>, and every chain calls
        /// the <i>same</i> <see cref="LogLikelihoodFunction"/> delegate instance. The log-likelihood is
        /// therefore invoked concurrently from multiple threads, as is the gradient delegate for
        /// gradient-based samplers. That delegate must be thread-safe or stateless. A likelihood that closes
        /// over mutable state — an automatic-differentiation tape, a reused workspace or buffer, a native
        /// solver handle, a cached factorization, a non-thread-safe PRNG — races by default, and the
        /// resulting corruption is silent: it surfaces as an implausible posterior rather than as an
        /// exception.
        /// </para>
        /// <para>
        /// <b>If your likelihood is not thread-safe, set this to false.</b> That is the supported answer, and
        /// it is the only one. Per-chain resources cannot be selected from inside the callback, because
        /// neither <see cref="LogLikelihood"/> nor the gradient delegate receives a chain index — they take
        /// only the parameter vector, so a callback has no way to learn which chain is calling it. Those
        /// signatures are part of the public API and are not going to change to add one.
        /// </para>
        /// </remarks>
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
        /// Indicates whether MAP initialization was attempted but failed, resulting in a fallback to random initialization.
        /// </summary>
        public bool MAPInitializationFailed { get; protected set; } = false;

        /// <summary>
        /// The Multivariate Normal proposal distribution set from the MAP estimate.
        /// </summary>
        protected MultivariateNormal? _MVN;

        /// <summary>
        /// Event is raised when the simulation progress changes.
        /// </summary>
        public event ProgressChangedEventHandler? ProgressChanged;

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
        public CancellationTokenSource? CancellationTokenSource { get; set; }

        #endregion

        #region Outputs

        /// <summary>
        /// Gets the population matrix used for population-based samplers.
        /// </summary>
        public List<ParameterSet> PopulationMatrix { get; protected set; } = null!;

        /// <summary>
        /// Gets the list of sampled Markov Chains.
        /// </summary>
        public List<ParameterSet>[] MarkovChains { get; protected set; } = null!;

        /// <summary>
        /// Keeps track of the number of accepted samples per chain.
        /// </summary>
        public int[] AcceptCount { get; protected set; } = null!;

        /// <summary>
        /// Keeps track of the number of calls to the proposal sampler per chain.
        /// </summary>
        public int[] SampleCount { get; protected set; } = null!;

        /// <summary>
        /// The acceptance rate per chain.
        /// </summary>
        public double[] AcceptanceRates
        {
            get
            {
                var acceptanceRates = new double[NumberOfChains];
                for (int i = 0; i < NumberOfChains; i++)
                {
                    acceptanceRates[i] = SampleCount[i] > 0
                        ? (double)AcceptCount[i] / SampleCount[i]
                        : 0d;
                }
                return acceptanceRates;
            }
        }

        /// <summary>
        /// The average log-likelihood across each chain for each iteration.
        /// </summary>
        public List<double> MeanLogLikelihood { get; protected set; } = null!;

        /// <summary>
        /// Gets and sets the number of posterior parameter sets to output.
        /// </summary>
        public int OutputLength { get; set; } = 10000;

        /// <summary>
        /// Output posterior parameter sets. These are recorded after the iterations have been completed.
        /// </summary>
        public List<ParameterSet>[] Output { get; protected set; } = null!;

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
            if (NumberOfChains < 1) throw new ArgumentException("There must be at least 1 chain.", nameof(NumberOfChains));
            if (Iterations < 100) throw new ArgumentException("The number of iterations cannot be less than 100.", nameof(Iterations));
            if (WarmupIterations < 1) throw new ArgumentException("The number of warm up iterations cannot be less than 1.", nameof(WarmupIterations));
            if (WarmupIterations > (int)(0.5 * Iterations)) throw new ArgumentException("The number of warm up iterations cannot be greater than half the number of iterations.", nameof(WarmupIterations));
            if (ThinningInterval < 1) throw new ArgumentException("The thinning interval cannot be less than 1.", nameof(ThinningInterval));
            if (InitialIterations < NumberOfChains) throw new ArgumentException("The initial population cannot be less than the number of chains.", nameof(InitialIterations));
            if (OutputLength < 100) throw new ArgumentException("The output length must be at least 100.", nameof(OutputLength));
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
        /// <returns>The initialized parameter state for each Markov chain.</returns>
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

            var prng = new MersenneTwister(PRNGSeed);
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
                var DE = new DifferentialEvolution((x) => { return LogLikelihoodFunction(x); }, NumberOfParameters, lowerBounds, upperBounds);
                DE.ReportFailure = false;
                DE.Maximize();
                if (DE.Status == OptimizationStatus.Success)
                {
                    try
                    {
                        _mapSuccessful = true;
                        // Get MAP
                        MAP = new ParameterSet((double[])DE.BestParameterSet.Values.Clone(), -DE.BestParameterSet.Fitness);
                        // Get Fisher Information Matrix
                        if (DE.Hessian == null) throw new InvalidOperationException("Hessian matrix is not available.");
                        var fisher = DE.Hessian * -1d;
                        // Invert it to get the covariance matrix, and scale to give wider coverage
                        var covar = fisher.Inverse() * 2;
                        
                        // Set up proposal distribution
                        _MVN = new MultivariateNormal(MAP.Values, covar.Array);
                        // Then randomly sample from the proposal
                        for (int i = 0; i < InitialIterations; i++)
                        {
                            parameters = _MVN.InverseCDF(rnds.GetRow(i));
                            logLH = LogLikelihoodFunction(parameters);
                            if (IsPopulationSampler) PopulationMatrix.Add(new ParameterSet((double[])parameters.Clone(), logLH));
                            tempPopulation.Add(new ParameterSet((double[])parameters.Clone(), logLH));
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
                        MAPInitializationFailed = true;
                    }
                }
            }


            // *** If not using MAP or if MAP fails, then use random initialization *** //

            // First add the mean of the priors
            for (int j = 0; j < NumberOfParameters; j++)
                parameters[j] = PriorDistributions[j].Mean;
            logLH = LogLikelihoodFunction(parameters);
            if (IsPopulationSampler) PopulationMatrix.Add(new ParameterSet((double[])parameters.Clone(), logLH));
            tempPopulation.Add(new ParameterSet((double[])parameters.Clone(), logLH));

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
                if (IsPopulationSampler) PopulationMatrix.Add(new ParameterSet((double[])parameters.Clone(), logLH));
                tempPopulation.Add(new ParameterSet((double[])parameters.Clone(), logLH));
            }
            
            // Sort temp population by log-likelihood in descending order.
            // List.Sort is an unstable introspective sort, and the chain starting states are taken from
            // the front of this list, so ties would decide the starting states and with them the entire
            // trajectory of every chain. A wide prior can leave many draws at exactly negative infinity.
            // OrderByDescending is a stable sort using the same default double comparison, so the order
            // of distinct log-likelihoods is unchanged and ties keep their draw order. Do not replace it
            // with Sort.
            tempPopulation = tempPopulation.OrderByDescending(x => x.Fitness).ToList();

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
                state = ChainIteration(index, state).Clone(false);
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
            if (!ResumeSimulation || _simulations < 1)
            {
                // Initialize the chains
                _chainStates = InitializeChains();

                // Initialize custom settings
                InitializeCustomSettings();
            }

            // Output settings. OutputIterations is shared with TransitionCount so that the advertised work
            // and the work performed here cannot drift apart.
            int outputIterations = OutputIterations;
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

                // Record output
                for (int j = 0; j < NumberOfChains; j++)
                {
                    // Update population
                    if (IsPopulationSampler == true)
                        PopulationMatrix.Add(_chainStates[j].Clone(false));

                    if (i <= Iterations)
                    {
                        // Record mean log-likelihood
                        if (MeanLogLikelihood.Count < i)
                            MeanLogLikelihood.Add(0);
                        MeanLogLikelihood[i - 1] += _chainStates[j].Fitness / NumberOfChains;

                        // Save chain state
                        MarkovChains[j].Add(_chainStates[j].Clone(false));
                    }
                    else if (i > Iterations && outputCount < OutputLength)
                    {
                        // Record the output and keep track of MAP
                        Output[j].Add(_chainStates[j].Clone(false));
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
                if (progress % Math.Max(1, (int)(totalIterations * ProgressChangedRate)) == 0)
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
            if (CancellationTokenSource != null && CancellationTokenSource.Token.CanBeCanceled == true)
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
            MAPInitializationFailed = false;
            // Clear old memory and re-instantiate the result storage
            _masterPRNG = new MersenneTwister(PRNGSeed);
            _chainPRNGs = new MersenneTwister[NumberOfChains];
            PopulationMatrix = new List<ParameterSet>();
            MarkovChains = new List<ParameterSet>[NumberOfChains];
            Output = new List<ParameterSet>[NumberOfChains];
            for (int i = 0; i < NumberOfChains; i++)
            {
                _chainPRNGs[i] = new MersenneTwister(_masterPRNG.Next());
                MarkovChains[i] = new List<ParameterSet>();
                Output[i] = new List<ParameterSet>();
            }
            AcceptCount = new int[NumberOfChains];
            SampleCount = new int[NumberOfChains];
            MeanLogLikelihood = new List<double>();
            MAP = new ParameterSet([], double.NegativeInfinity);
        }

        #endregion

    }
}
