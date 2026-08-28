using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using System;
using System.Collections.Generic;

namespace Numerics.Sampling.MCMC
{

    /// <summary>
    /// The No-U-Turn Sampler (NUTS), an adaptive extension of Hamiltonian Monte Carlo that
    /// automatically tunes the trajectory length.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// NUTS eliminates the need to manually tune the number of leapfrog steps required by standard HMC.
    /// It builds a balanced binary tree of leapfrog states by repeatedly doubling the trajectory in a
    /// randomly chosen direction. The trajectory stops when a U-turn is detected (the trajectory starts
    /// heading back towards the starting point) or when a maximum tree depth is reached.
    /// A candidate state is selected via multinomial sampling weighted by the exponential of the
    /// negative Hamiltonian.
    /// </para>
    /// <para>
    /// During warmup, the leapfrog step size is automatically adapted using dual averaging
    /// to achieve a target Metropolis acceptance probability of approximately 80%. A diagonal mass matrix
    /// is estimated during warmup using Stan-style windowed adaptation (Welford's online algorithm) to
    /// precondition the Hamiltonian dynamics for multi-scale posteriors.
    /// </para>
    /// <para>
    /// NUTS is suitable for Bayesian parameter estimation in models where manual tuning of HMC
    /// trajectory lengths is impractical, including hierarchical and mixture models.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Hoffman, M.D. and Gelman, A. (2014). "The No-U-Turn Sampler: Adaptively Setting Path Lengths
    /// in Hamiltonian Monte Carlo." Journal of Machine Learning Research, 15, 1593-1623.
    /// </description></item>
    /// <item><description>
    /// Betancourt, M. (2017). "A Conceptual Introduction to Hamiltonian Monte Carlo." arXiv:1701.02434.
    /// </description></item>
    /// <item><description>
    /// Stan Development Team (2024). Stan Reference Manual, Section 15.2: HMC Algorithm Parameters.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Hamiltonian_Monte_Carlo#No_U-Turn_Sampler"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class NUTS : MCMCSampler
    {

        /// <summary>
        /// Constructs a new NUTS sampler.
        /// </summary>
        /// <param name="priorDistributions">The list of prior distributions for the model parameters.</param>
        /// <param name="logLikelihoodFunction">The log-likelihood function to evaluate.</param>
        /// <param name="mass">Optional. The initial mass vector for the momentum distribution. Default = Identity. Will be adapted during warmup.</param>
        /// <param name="stepSize">Optional. The initial leapfrog step size. Will be adapted during warmup. Default = 0.1.</param>
        /// <param name="maxTreeDepth">Optional. The maximum binary tree depth. Default = 10.</param>
        /// <param name="gradientFunction">Optional. The function for evaluating the gradient of the log-likelihood. Numerical finite difference will be used by default.</param>
        public NUTS(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction, Vector? mass = null, double stepSize = 0.1, int maxTreeDepth = 10, HMC.Gradient? gradientFunction = null) : base(priorDistributions, logLikelihoodFunction)
        {
            InitialIterations = 100 * NumberOfParameters;

            // Set the initial mass vector (will be overridden by adaptation during warmup)
            if (mass == null)
            {
                Mass = new Vector(NumberOfParameters, 1d);
            }
            else
            {
                Mass = mass;
            }

            // Set the inverse mass vector
            _inverseMass = new Vector(NumberOfParameters);
            for (int i = 0; i < NumberOfParameters; i++)
                _inverseMass[i] = 1d / Mass[i];

            // Set defaults
            _initialStepSize = stepSize;
            MaxTreeDepth = maxTreeDepth;

            // Cache prior distribution bounds for the gradient function
            _lowerBounds = new double[NumberOfParameters];
            _upperBounds = new double[NumberOfParameters];
            for (int i = 0; i < NumberOfParameters; i++)
            {
                _lowerBounds[i] = priorDistributions[i].Minimum;
                _upperBounds[i] = priorDistributions[i].Maximum;
            }

            // Set the gradient function with prior bounds so finite-difference probes stay in valid region
            if (gradientFunction == null)
            {
                GradientFunction = (x) => new Vector(NumericalDerivative.Gradient((y) => SafeLogLikelihood(y), x.ToArray(), _lowerBounds, _upperBounds));
            }
            else
            {
                GradientFunction = gradientFunction;
            }
        }

        private Vector _inverseMass;
        private double _initialStepSize;
        private double[] _lowerBounds;
        private double[] _upperBounds;

        // Per-chain dual averaging state
        private double[] _chainStepSizes = null!;
        private double[] _chainLogEpsBar = null!;
        private double[] _chainHBar = null!;
        private double[] _chainMu = null!;
        private int[] _chainAdaptStep = null!;

        // Per-chain post-warmup diagnostics. These are streaming accumulators so
        // diagnostic collection adds no target/gradient evaluations or draw storage.
        private double[] _hamiltonianAcceptanceSums = Array.Empty<double>();
        private int[] _diagnosticSampleCounts = Array.Empty<int>();
        private int[] _divergenceCounts = Array.Empty<int>();
        private int[] _maxTreeDepthHitCounts = Array.Empty<int>();
        private double[] _treeDepthSums = Array.Empty<double>();
        private double[] _leapfrogStepSums = Array.Empty<double>();
        private double[] _energyMeans = Array.Empty<double>();
        private double[] _energyM2 = Array.Empty<double>();
        private double[] _energySquaredDifferenceSums = Array.Empty<double>();
        private double[] _previousEnergy = Array.Empty<double>();
        private bool[] _hasPreviousEnergy = Array.Empty<bool>();

        // Dual averaging hyperparameters (Hoffman & Gelman 2014, Section 3.2).
        // The adaptation reads TargetAcceptanceRate; DELTA_TARGET is that property's default.
        private const double DELTA_TARGET = 0.80;
        private const double GAMMA = 0.05;
        private const double T0 = 10.0;
        private const double KAPPA = 0.75;

        // Divergence threshold: if H - H0 exceeds this, the trajectory is considered divergent
        private const double MAX_DELTA_H = 1000.0;

        // Per-chain diagonal mass matrix adaptation (Welford's online algorithm)
        private double[][] _welfordMean = null!;
        private double[][] _welfordM2 = null!;
        private int[] _welfordCount = null!;
        private double[][] _massMatrix = null!;
        private double[][] _inverseMassMatrix = null!;

        // Adaptation window boundaries (computed per chain in InitializeCustomSettings)
        private int _initBuffer;
        private int _termBuffer;
        private int[] _adaptWindowEnds = null!;

        /// <summary>
        /// The number of recently evaluated positions each chain's gradient memo retains.
        /// </summary>
        /// <remarks>
        /// A doubling extends the trajectory in one direction, so the join a memo has to cover is
        /// between two consecutive leaves and only one entry is strictly needed. A size of four leaves
        /// room for the step-size heuristic, which re-enters the same starting position on every trial
        /// step, and for the forward and backward endpoints to be resident at the same time, at a cost
        /// of 8 × <see cref="MCMCSampler.NumberOfParameters"/> doubles per chain.
        /// </remarks>
        private const int GRADIENT_CACHE_SIZE = 4;

        // Per-chain gradient memo: positions, the gradient at each, and which slots hold a value.
        // Keyed by chain because Sample() runs chains under Parallel.For. See EvaluateGradient.
        private double[][][] _gradientCachePositions = null!;
        private double[][][] _gradientCacheValues = null!;
        private bool[][] _gradientCacheOccupied = null!;
        private int[] _gradientCacheNextSlot = null!;

        /// <summary>
        /// The mass vector for the momentum distribution.
        /// </summary>
        /// <remarks>
        /// When <see cref="AdaptMassMatrix"/> is <see langword="true"/> (the default), the supplied
        /// mass seeds the metric and is replaced at the end of each adaptation window; set
        /// <see cref="AdaptMassMatrix"/> to <see langword="false"/> to sample with the fixed
        /// supplied metric.
        /// </remarks>
        public Vector Mass { get; }

        /// <summary>
        /// The maximum binary tree depth. Default = 10. This caps the trajectory length at 2^MaxTreeDepth leapfrog steps.
        /// </summary>
        public int MaxTreeDepth { get; set; }

        /// <summary>
        /// The function for evaluating the gradient of the log-likelihood.
        /// </summary>
        public HMC.Gradient GradientFunction { get; }

        /// <summary>
        /// Gets or sets the target Metropolis acceptance probability that dual averaging adapts the
        /// leapfrog step size toward during warmup. Must be strictly between 0 and 1. Default = 0.80.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The default of 0.80 is the value recommended by Hoffman and Gelman (2014) and is what Stan's
        /// <c>adapt_delta</c> defaults to.
        /// </para>
        /// <para>
        /// Raising it toward 0.90 or 0.95 makes dual averaging settle on a shorter step size, which
        /// integrates the Hamiltonian more accurately and is the standard response to a run that reports
        /// divergent transitions through <see cref="DivergenceCounts"/>. The cost is proportional: a
        /// shorter step needs more leapfrog steps, and therefore more gradient evaluations, to cover the
        /// same trajectory length. A target close to 1 can drive the step size to the lower clamp and
        /// exhaust <see cref="MaxTreeDepth"/> on every transition, so raise it in steps and watch
        /// <see cref="MeanLeapfrogSteps"/> alongside <see cref="DivergenceCounts"/>.
        /// </para>
        /// <para>
        /// Persistent divergences that survive a raised target usually indicate a posterior geometry the
        /// step size cannot fix on its own, such as a funnel, rather than an adaptation failure.
        /// </para>
        /// <para>
        /// The value is validated when sampling starts, not at assignment.
        /// </para>
        /// </remarks>
        public double TargetAcceptanceRate { get; set; } = DELTA_TARGET;

        /// <summary>
        /// Gets the mean post-warmup Hamiltonian acceptance probability for each chain.
        /// </summary>
        /// <remarks>
        /// This is the mean tree acceptance statistic used by dual averaging, not the
        /// fraction of NUTS iterations that returned a retained state.
        /// </remarks>
        public double[] HamiltonianAcceptanceRates
        {
            get
            {
                if (_diagnosticSampleCounts.Length != NumberOfChains ||
                    _hamiltonianAcceptanceSums.Length != NumberOfChains)
                    return Array.Empty<double>();

                var values = new double[NumberOfChains];
                for (int i = 0; i < NumberOfChains; i++)
                {
                    values[i] = _diagnosticSampleCounts[i] == 0
                        ? 0d
                        : _hamiltonianAcceptanceSums[i] / _diagnosticSampleCounts[i];
                }
                return values;
            }
        }

        /// <summary>
        /// Gets the number of post-warmup transitions contributing diagnostics per chain.
        /// </summary>
        public int[] DiagnosticSampleCounts => _diagnosticSampleCounts.Length != NumberOfChains
            ? Array.Empty<int>()
            : (int[])_diagnosticSampleCounts.Clone();

        /// <summary>
        /// Gets the number of divergent post-warmup transitions per chain.
        /// </summary>
        public int[] DivergenceCounts => _divergenceCounts.Length != NumberOfChains
            ? Array.Empty<int>()
            : (int[])_divergenceCounts.Clone();

        /// <summary>
        /// Gets the number of post-warmup transitions that exhausted <see cref="MaxTreeDepth"/> per chain.
        /// </summary>
        public int[] MaxTreeDepthHitCounts => _maxTreeDepthHitCounts.Length != NumberOfChains
            ? Array.Empty<int>()
            : (int[])_maxTreeDepthHitCounts.Clone();

        /// <summary>
        /// Gets the mean post-warmup tree depth for each chain.
        /// </summary>
        public double[] MeanTreeDepths => ComputeDiagnosticMeans(_treeDepthSums);

        /// <summary>
        /// Gets the mean number of post-warmup leapfrog steps per transition for each chain.
        /// </summary>
        public double[] MeanLeapfrogSteps => ComputeDiagnosticMeans(_leapfrogStepSums);

        /// <summary>
        /// Gets the current adapted leapfrog step size for each chain.
        /// </summary>
        public double[] StepSizes => _chainStepSizes == null
            ? Array.Empty<double>()
            : (double[])_chainStepSizes.Clone();

        /// <summary>
        /// Gets the post-warmup energy Bayesian fraction of missing information for each chain.
        /// </summary>
        /// <remarks>
        /// E-BFMI is the mean squared successive Hamiltonian difference divided by
        /// sample variance of the Hamiltonian, matching the Stan diagnostic definition.
        /// </remarks>
        public double[] EnergyBayesianFractionOfMissingInformation
        {
            get
            {
                if (_diagnosticSampleCounts.Length != NumberOfChains ||
                    _energyM2.Length != NumberOfChains ||
                    _energySquaredDifferenceSums.Length != NumberOfChains)
                    return Array.Empty<double>();

                var values = new double[NumberOfChains];
                for (int i = 0; i < NumberOfChains; i++)
                {
                    values[i] = ComputeEnergyBayesianFractionOfMissingInformation(
                        _diagnosticSampleCounts[i],
                        _energyM2[i],
                        _energySquaredDifferenceSums[i]);
                }
                return values;
            }
        }


        /// <summary>
        /// Computes per-chain diagnostic means from streaming sums.
        /// </summary>
        /// <param name="sums">Per-chain diagnostic sums.</param>
        /// <returns>The corresponding per-chain means.</returns>
        private double[] ComputeDiagnosticMeans(double[] sums)
        {
            if (sums.Length != NumberOfChains || _diagnosticSampleCounts.Length != NumberOfChains)
                return Array.Empty<double>();

            var values = new double[NumberOfChains];
            for (int i = 0; i < NumberOfChains; i++)
            {
                if (_diagnosticSampleCounts[i] > 0)
                    values[i] = sums[i] / _diagnosticSampleCounts[i];
            }
            return values;
        }

        /// <summary>
        /// Gets or sets whether to adapt the diagonal mass matrix during warmup. Default = true.
        /// </summary>
        /// <remarks>
        /// <para>
        /// When enabled, the sampler uses Stan-style windowed adaptation with Welford's online algorithm
        /// to estimate the posterior variance of each parameter during warmup, and takes that variance as
        /// the coordinate's inverse mass. Each leapfrog step then scales with the coordinate's own posterior
        /// width, which is what lets one step size serve a posterior whose parameters differ in scale.
        /// </para>
        /// <para>
        /// With an identity metric the step size must track the narrowest posterior direction while the
        /// trajectory has to span the widest, so on an ill-conditioned posterior NUTS saturates
        /// <see cref="MaxTreeDepth"/> on nearly every transition; adaptation removes that failure mode.
        /// On small, well-conditioned fits the metric has little to correct and the adaptation can cost
        /// up to about 38% more leapfrog steps per transition.
        /// </para>
        /// <para>
        /// Set this to <see langword="false"/> to sample with the fixed metric supplied through
        /// <see cref="Mass"/>.
        /// </para>
        /// </remarks>
        public bool AdaptMassMatrix { get; set; } = true;

        /// <inheritdoc/>
        protected override void ValidateCustomSettings()
        {
            if (Mass.Length != NumberOfParameters) throw new ArgumentException("The mass vector must be the same length as the number of parameters.", nameof(Mass));
            if (_initialStepSize <= 0) throw new ArgumentException("The leapfrog step size must be positive.", "stepSize");
            if (MaxTreeDepth < 1) throw new ArgumentException("The maximum tree depth must be at least 1.", nameof(MaxTreeDepth));
            if (!Tools.IsFinite(TargetAcceptanceRate) || TargetAcceptanceRate <= 0d || TargetAcceptanceRate >= 1d) throw new ArgumentException("The target acceptance rate must be greater than 0 and less than 1.", nameof(TargetAcceptanceRate));
        }

        /// <inheritdoc/>
        protected override void InitializeCustomSettings()
        {
            int D = NumberOfParameters;
            int N = NumberOfChains;

            // Initialize dual averaging state
            _chainStepSizes = new double[N];
            _chainLogEpsBar = new double[N];
            _chainHBar = new double[N];
            _chainMu = new double[N];
            _chainAdaptStep = new int[N];

            // Initialize post-warmup diagnostic accumulators.
            _hamiltonianAcceptanceSums = new double[N];
            _diagnosticSampleCounts = new int[N];
            _divergenceCounts = new int[N];
            _maxTreeDepthHitCounts = new int[N];
            _treeDepthSums = new double[N];
            _leapfrogStepSums = new double[N];
            _energyMeans = new double[N];
            _energyM2 = new double[N];
            _energySquaredDifferenceSums = new double[N];
            _previousEnergy = new double[N];
            _hasPreviousEnergy = new bool[N];

            // Initialize diagonal mass matrix and Welford accumulators
            _welfordMean = new double[N][];
            _welfordM2 = new double[N][];
            _welfordCount = new int[N];
            _massMatrix = new double[N][];
            _inverseMassMatrix = new double[N][];

            // Allocate the gradient memo empty. Allocating it here is also what resets it, so a
            // re-run cannot serve an entry from the previous run.
            _gradientCachePositions = new double[N][][];
            _gradientCacheValues = new double[N][][];
            _gradientCacheOccupied = new bool[N][];
            _gradientCacheNextSlot = new int[N];

            for (int i = 0; i < N; i++)
            {
                _gradientCachePositions[i] = new double[GRADIENT_CACHE_SIZE][];
                _gradientCacheValues[i] = new double[GRADIENT_CACHE_SIZE][];
                _gradientCacheOccupied[i] = new bool[GRADIENT_CACHE_SIZE];
                for (int s = 0; s < GRADIENT_CACHE_SIZE; s++)
                {
                    _gradientCachePositions[i][s] = new double[D];
                    _gradientCacheValues[i][s] = new double[D];
                }

                // Start with identity mass matrix (or user-provided mass)
                _massMatrix[i] = new double[D];
                _inverseMassMatrix[i] = new double[D];
                for (int j = 0; j < D; j++)
                {
                    _massMatrix[i][j] = Mass[j];
                    _inverseMassMatrix[i][j] = 1.0 / Mass[j];
                }

                _welfordMean[i] = new double[D];
                _welfordM2[i] = new double[D];
                _welfordCount[i] = 0;
            }

            // Use Stan's initial buffer, doubling mass-matrix windows, and terminal buffer.
            // Default lengths are init_buffer=75, base_window=25, and term_buffer=50.
            int totalWarmup = WarmupIterations * ThinningInterval;

            // Stan default window sizes, scaled if warmup is too short
            _initBuffer = 75;
            _termBuffer = 50;
            int baseWindow = 25;
            if (_initBuffer + baseWindow + _termBuffer > totalWarmup)
            {
                // Fallback: redistribute proportionally
                _initBuffer = Math.Max(1, (int)(0.15 * totalWarmup));
                _termBuffer = Math.Max(1, (int)(0.10 * totalWarmup));
                baseWindow = Math.Max(1, totalWarmup - _initBuffer - _termBuffer);
            }

            // Build doubling window boundaries with Stan's look-ahead merging.
            // If the next-next window would overshoot, stretch the current window.
            var windowEnds = new List<int>();
            int adaptEnd = totalWarmup - _termBuffer;
            if (adaptEnd > _initBuffer)
            {
                int windowSize = baseWindow;
                int nextWindowEnd = _initBuffer + windowSize - 1;

                while (nextWindowEnd < adaptEnd)
                {
                    // Look ahead: if (current end + 2 * next window size) > adaptEnd,
                    // stretch this window to fill the remaining space
                    int nextNextEnd = nextWindowEnd + 2 * windowSize;
                    if (nextNextEnd >= adaptEnd)
                    {
                        nextWindowEnd = adaptEnd - 1;
                    }
                    windowEnds.Add(nextWindowEnd);

                    windowSize *= 2;
                    nextWindowEnd += windowSize;
                }
                // Ensure we always have the final boundary
                if (windowEnds.Count == 0 || windowEnds[windowEnds.Count - 1] != adaptEnd - 1)
                {
                    windowEnds.Add(adaptEnd - 1);
                }
            }
            _adaptWindowEnds = windowEnds.ToArray();

            // Find reasonable initial step size per chain (Hoffman & Gelman Algorithm 4)
            // and initialize dual averaging state.
            // Note: _chainStates is populated by InitializeChains() before this method is called.
            for (int i = 0; i < N; i++)
            {
                double[] theta0 = _chainStates[i].Values;
                double logLH0 = _chainStates[i].Fitness;

                double eps0;
                try
                {
                    eps0 = FindReasonableEpsilon(theta0, logLH0, i);
                }
                catch
                {
                    eps0 = _initialStepSize;
                }

                _chainStepSizes[i] = eps0;
                _chainLogEpsBar[i] = Math.Log(eps0);
                _chainHBar[i] = 0.0;
                _chainMu[i] = Math.Log(10.0 * eps0);
                _chainAdaptStep[i] = 0;
            }
        }

        /// <summary>
        /// Finds a reasonable initial step size using the heuristic from Hoffman and Gelman (2014), Algorithm 4.
        /// Searches for a step size that gives roughly 50% acceptance probability for a single leapfrog step.
        /// </summary>
        /// <param name="theta0">The initial position.</param>
        /// <param name="logLH0">The log-likelihood at the initial position.</param>
        /// <param name="chainIndex">The chain index for RNG and mass matrix access.</param>
        /// <returns>A reasonable initial step size.</returns>
        private double FindReasonableEpsilon(double[] theta0, double logLH0, int chainIndex)
        {
            int D = NumberOfParameters;
            double epsilon = 1.0;

            // Sample momentum from current mass matrix
            var r0 = new double[D];
            for (int j = 0; j < D; j++)
                r0[j] = Math.Sqrt(_massMatrix[chainIndex][j]) * Normal.StandardZ(_chainPRNGs[chainIndex].NextDouble());

            // Compute initial Hamiltonian
            double H0 = -logLH0 + 0.5 * DiagonalQuadraticForm(r0, _inverseMassMatrix[chainIndex]);

            var thetaPrime = (double[])theta0.Clone();
            var rPrime = (double[])r0.Clone();
            double logAlpha = TrySingleStepLogAcceptance(theta0, r0, thetaPrime, rPrime, epsilon, H0, chainIndex);

            // Determine direction: double epsilon if acceptance too high, halve if too low
            double a = logAlpha > Math.Log(0.5) ? 1.0 : -1.0;

            for (int iter = 0; iter < 25; iter++)
            {
                // Check if we've crossed the 0.5 threshold
                if (a * logAlpha <= -a * Math.Log(2.0))
                    break;

                epsilon *= Math.Pow(2.0, a);

                // Safety bounds: don't let epsilon get absurdly small or large
                if (epsilon < 1e-8 || epsilon > 1e6)
                    break;

                Array.Copy(theta0, thetaPrime, D);
                Array.Copy(r0, rPrime, D);
                logAlpha = TrySingleStepLogAcceptance(theta0, r0, thetaPrime, rPrime, epsilon, H0, chainIndex);
            }

            return Math.Max(1e-8, Math.Min(epsilon, 1e6));
        }

        /// <summary>
        /// Attempts one leapfrog step and converts invalid Hamiltonian states to a low log-acceptance value.
        /// </summary>
        /// <param name="theta0">The starting position.</param>
        /// <param name="r0">The starting momentum.</param>
        /// <param name="thetaPrime">Working buffer for the proposed position.</param>
        /// <param name="rPrime">Working buffer for the proposed momentum.</param>
        /// <param name="epsilon">The leapfrog step size.</param>
        /// <param name="initialHamiltonian">The initial Hamiltonian.</param>
        /// <param name="chainIndex">The chain index for mass matrix access.</param>
        /// <returns>The log acceptance statistic for the trial step, or a low value when the trial is invalid.</returns>
        private double TrySingleStepLogAcceptance(double[] theta0, double[] r0, double[] thetaPrime, double[] rPrime, double epsilon, double initialHamiltonian, int chainIndex)
        {
            try
            {
                LeapfrogInPlace(thetaPrime, rPrime, epsilon, chainIndex);

                double logLH = SafeLogLikelihood(thetaPrime);
                double hamiltonian = -logLH + 0.5 * DiagonalQuadraticForm(rPrime, _inverseMassMatrix[chainIndex]);
                double logAlpha = initialHamiltonian - hamiltonian;
                if (!Tools.IsFinite(logAlpha))
                    return -1000.0;

                return logAlpha;
            }
            catch (ArithmeticException)
            {
                Array.Copy(theta0, thetaPrime, theta0.Length);
                Array.Copy(r0, rPrime, r0.Length);
                return -1000.0;
            }
        }

        /// <summary>
        /// Evaluates the gradient of the log-likelihood at a position, reusing a value this chain
        /// computed at the identical position within the last <see cref="GRADIENT_CACHE_SIZE"/>
        /// evaluations.
        /// </summary>
        /// <param name="position">The position to evaluate at. Not retained; a copy is stored.</param>
        /// <param name="chainIndex">The chain index whose memo is consulted.</param>
        /// <returns>
        /// The gradient. The array is owned by the memo and must be treated as read-only by the caller.
        /// It is valid only until the next miss on this chain: a hit does not advance the ring, so a hit
        /// on the slot the ring is currently pointing at is overwritten by the very next miss. Both
        /// callers consume the array before evaluating again.
        /// </returns>
        /// <remarks>
        /// <para>
        /// Consecutive leaves of a doubling chain end to end, so each leapfrog step opens from the
        /// position the previous step landed on and would otherwise recompute its gradient, and the
        /// step-size heuristic re-enters the same starting position on every trial step. With the
        /// default finite-difference gradient each avoided evaluation saves on the order of
        /// 2 × <see cref="MCMCSampler.NumberOfParameters"/> log-likelihood evaluations.
        /// </para>
        /// <para>
        /// Positions are compared bitwise through <see cref="BitConverter.DoubleToInt64Bits(double)"/>,
        /// never with <c>==</c> and never within a tolerance, so a hit is possible only at the exact
        /// point the delegate was called on and returns exactly the value a recomputation would produce,
        /// with <c>+0</c> and <c>-0</c> kept distinct. The memo requires the gradient delegate to be a
        /// deterministic function of its argument, which a seeded, reproducible run already requires.
        /// </para>
        /// <para>
        /// The metric is not part of the key: <see cref="GradientFunction"/> takes a position and
        /// nothing else, so the gradient does not depend on the mass matrix, the step size, or the
        /// adaptation window, and an entry recorded under one metric remains valid after
        /// <see cref="UpdateMassMatrix"/> replaces the metric. The chain index is part of the key for
        /// thread safety: every array is indexed by chain first, and <see cref="MCMCSampler.Sample"/>
        /// gives each chain index to exactly one <see cref="System.Threading.Tasks.Parallel"/> iteration
        /// at a time.
        /// </para>
        /// <para>
        /// Both the position and the gradient are copied into the memo: the caller's position array is
        /// mutated in place by <see cref="LeapfrogInPlace"/>, and a caller-supplied gradient delegate is
        /// free to return the same <see cref="Vector"/> every call, so neither may be retained by
        /// reference. The position is copied before the delegate runs, so a delegate that writes through
        /// its argument cannot pair a mutated position with the gradient of a different point. A
        /// gradient whose length does not match the parameter count is returned without being stored.
        /// </para>
        /// </remarks>
        private double[] EvaluateGradient(double[] position, int chainIndex)
        {
            int D = NumberOfParameters;
            var positions = _gradientCachePositions[chainIndex];
            var values = _gradientCacheValues[chainIndex];
            var occupied = _gradientCacheOccupied[chainIndex];

            for (int slot = 0; slot < GRADIENT_CACHE_SIZE; slot++)
            {
                if (!occupied[slot]) continue;

                var stored = positions[slot];
                bool identical = true;
                for (int j = 0; j < D; j++)
                {
                    if (BitConverter.DoubleToInt64Bits(stored[j]) != BitConverter.DoubleToInt64Bits(position[j]))
                    {
                        identical = false;
                        break;
                    }
                }

                if (identical) return values[slot];
            }

            // A miss claims a slot and records the position before the delegate runs, and the slot stays
            // marked empty for the whole window, so a throw, a malformed length, or a delegate that
            // writes through its argument leaves nothing that could be matched later.
            int next = _gradientCacheNextSlot[chainIndex];
            var slotPosition = positions[next];
            var slotValue = values[next];

            occupied[next] = false;
            for (int j = 0; j < D; j++)
                slotPosition[j] = position[j];

            double[] gradient = GradientFunction(position).Array;
            if (gradient.Length != D) return gradient;

            for (int j = 0; j < D; j++)
                slotValue[j] = gradient[j];

            occupied[next] = true;
            _gradientCacheNextSlot[chainIndex] = (next + 1) % GRADIENT_CACHE_SIZE;

            return slotValue;
        }

        /// <summary>
        /// Performs a single leapfrog step in-place on raw arrays, using the per-chain mass matrix.
        /// Used by FindReasonableEpsilon to avoid Vector allocations.
        /// </summary>
        /// <param name="theta">Position array (modified in-place).</param>
        /// <param name="momentum">Momentum array (modified in-place).</param>
        /// <param name="epsilon">The step size.</param>
        /// <param name="chainIndex">Chain index for mass matrix access.</param>
        private void LeapfrogInPlace(double[] theta, double[] momentum, double epsilon, int chainIndex)
        {
            int D = NumberOfParameters;
            double halfEps = epsilon * 0.5;
            double[] invMass = _inverseMassMatrix[chainIndex];

            // Half-step momentum update. grad is the memo's own buffer and must be read only; write
            // through it and this chain's entry is silently corrupted.
            double[] grad = EvaluateGradient(theta, chainIndex);
            for (int j = 0; j < D; j++)
                momentum[j] += grad[j] * halfEps;

            // Full-step position update
            for (int j = 0; j < D; j++)
            {
                theta[j] += invMass[j] * momentum[j] * epsilon;
                if (theta[j] < _lowerBounds[j])
                    theta[j] = _lowerBounds[j] + Tools.DoubleMachineEpsilon;
                if (theta[j] > _upperBounds[j])
                    theta[j] = _upperBounds[j] - Tools.DoubleMachineEpsilon;
            }

            // Half-step momentum update. As above, grad is memo-owned and must be read only.
            grad = EvaluateGradient(theta, chainIndex);
            for (int j = 0; j < D; j++)
                momentum[j] += grad[j] * halfEps;
        }

        /// <summary>
        /// Computes the diagonal quadratic form φᵀ M⁻¹ φ using raw arrays.
        /// </summary>
        /// <param name="momentum">The momentum vector.</param>
        /// <param name="inverseMass">The diagonal inverse mass matrix.</param>
        /// <returns>The scalar result Σᵢ momentum[i]² × inverseMass[i].</returns>
        private static double DiagonalQuadraticForm(double[] momentum, double[] inverseMass)
        {
            double sum = 0;
            for (int j = 0; j < momentum.Length; j++)
                sum += momentum[j] * momentum[j] * inverseMass[j];
            return sum;
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {
            // Update the sample count
            SampleCount[index] += 1;
            int sampleNum = SampleCount[index];
            int warmupSteps = WarmupIterations * ThinningInterval;

            double eps = _chainStepSizes[index];
            int D = NumberOfParameters;

            // Sample momentum from N(0, M) using the per-chain mass matrix.
            var phi = new Vector(D);
            for (int i = 0; i < D; i++)
                phi[i] = Math.Sqrt(_massMatrix[index][i]) * Normal.StandardZ(_chainPRNGs[index].NextDouble());

            // Compute initial Hamiltonian using per-chain inverse mass matrix
            double H0 = -state.Fitness + 0.5 * DiagonalQuadraticFormVec(phi, _inverseMassMatrix[index]);

            // Initialize the trajectory tree.
            var theta = new Vector(state.Values);
            var thetaMinus = theta.Clone();
            var thetaPlus = theta.Clone();
            var rMinus = phi.Clone();
            var rPlus = phi.Clone();

            var candidate = theta.Clone();
            double candidateLogLH = state.Fitness;
            double logSumWeight = -H0;

            int depth = 0;
            double sumAlpha = 0;
            int numAlpha = 0;
            int leapfrogSteps = 0;
            int trajectoryDepth = 0;
            bool trajectoryDivergent = false;

            // Double the trajectory tree until a U-turn or the maximum depth.
            while (depth < MaxTreeDepth)
            {
                // Choose a random direction
                int v = _chainPRNGs[index].NextDouble() < 0.5 ? -1 : 1;

                TreeState subtree;
                if (v == -1)
                {
                    subtree = BuildTree(thetaMinus, rMinus, -eps, depth, H0, index);
                    thetaMinus = subtree.ThetaMinus;
                    rMinus = subtree.MomentumMinus;
                }
                else
                {
                    subtree = BuildTree(thetaPlus, rPlus, eps, depth, H0, index);
                    thetaPlus = subtree.ThetaPlus;
                    rPlus = subtree.MomentumPlus;
                }

                // If the subtree is valid, consider accepting its candidate
                if (subtree.Valid)
                {
                    double logSumWeightNew = LogSumExp(logSumWeight, subtree.LogSumWeight);
                    double acceptProb = Math.Exp(subtree.LogSumWeight - logSumWeightNew);
                    if (_chainPRNGs[index].NextDouble() < acceptProb)
                    {
                        candidate = subtree.ThetaPrime;
                        candidateLogLH = subtree.LogLikelihoodPrime;
                    }
                    logSumWeight = logSumWeightNew;
                }

                leapfrogSteps += subtree.LeafCount;
                trajectoryDepth = depth + 1;
                trajectoryDivergent |= subtree.Divergent;
                // Accumulate adaptation statistics
                sumAlpha += subtree.SumAlpha;
                numAlpha += subtree.NumAlpha;

                // Check stopping criterion: divergence or U-turn at the top level
                if (!subtree.Valid)
                    break;

                var dTheta = thetaPlus - thetaMinus;
                if (Vector.DotProduct(dTheta, rMinus) < 0 || Vector.DotProduct(dTheta, rPlus) < 0)
                    break;

                depth++;
            }

            double averageAcceptanceProbability = numAlpha > 0 ? sumAlpha / numAlpha : 0d;
            // Adapt the step size and mass matrix during warmup.
            if (sampleNum <= warmupSteps)
            {
                // Always do dual averaging step size adaptation during warmup
                // Preserve the original neutral fallback when no subtree contributed.
                double adaptationAcceptanceProbability = numAlpha > 0
                    ? averageAcceptanceProbability
                    : TargetAcceptanceRate;
                DualAveragingUpdate(index, adaptationAcceptanceProbability);

                // Accumulate Welford statistics during mass-matrix adaptation windows.
                if (AdaptMassMatrix && sampleNum > _initBuffer && sampleNum <= warmupSteps - _termBuffer)
                {
                    AccumulateWelfordStatistics(index, candidate.Array);

                    // Check if we're at the end of an adaptation window
                    if (IsEndOfAdaptationWindow(sampleNum))
                    {
                        var currentState = new ParameterSet(candidate.Array, candidateLogLH);
                        UpdateMassMatrix(index, currentState);
                    }
                }
            }
            else if (sampleNum == warmupSteps + 1)
            {
                // After warmup, fix step size to the smoothed value
                _chainStepSizes[index] = Math.Exp(_chainLogEpsBar[index]);
            }

            if (sampleNum > warmupSteps)
            {
                RecordDiagnostics(
                    index,
                    averageAcceptanceProbability,
                    trajectoryDivergent,
                    depth >= MaxTreeDepth,
                    trajectoryDepth,
                    leapfrogSteps,
                    H0);
            }

            // NUTS always accepts
            AcceptCount[index] += 1;
            return new ParameterSet(candidate.Array, candidateLogLH);
        }

        /// <summary>
        /// Records one post-warmup NUTS transition using constant-memory accumulators.
        /// </summary>
        /// <param name="chainIndex">Zero-based chain index.</param>
        /// <param name="acceptanceProbability">Mean tree acceptance probability.</param>
        /// <param name="divergent">Whether the trajectory contained a divergence.</param>
        /// <param name="hitMaximumDepth">Whether tree construction exhausted the configured maximum depth.</param>
        /// <param name="treeDepth">Tree depth attempted by the transition.</param>
        /// <param name="leapfrogSteps">Number of leapfrog steps built by the transition.</param>
        /// <param name="energy">Initial Hamiltonian after momentum resampling.</param>
        private void RecordDiagnostics(int chainIndex, double acceptanceProbability, bool divergent,
            bool hitMaximumDepth, int treeDepth, int leapfrogSteps, double energy)
        {
            int newCount = ++_diagnosticSampleCounts[chainIndex];
            _hamiltonianAcceptanceSums[chainIndex] += acceptanceProbability;
            if (divergent)
                _divergenceCounts[chainIndex]++;
            if (hitMaximumDepth)
                _maxTreeDepthHitCounts[chainIndex]++;
            _treeDepthSums[chainIndex] += treeDepth;
            _leapfrogStepSums[chainIndex] += leapfrogSteps;

            if (_hasPreviousEnergy[chainIndex])
            {
                double energyDifference = energy - _previousEnergy[chainIndex];
                _energySquaredDifferenceSums[chainIndex] += energyDifference * energyDifference;
            }
            _previousEnergy[chainIndex] = energy;
            _hasPreviousEnergy[chainIndex] = true;

            double delta = energy - _energyMeans[chainIndex];
            _energyMeans[chainIndex] += delta / newCount;
            double centeredEnergy = energy - _energyMeans[chainIndex];
            _energyM2[chainIndex] += delta * centeredEnergy;
        }

        /// <summary>
        /// Computes E-BFMI from streaming energy accumulators.
        /// </summary>
        /// <param name="sampleCount">Number of energy observations.</param>
        /// <param name="energyM2">Sum of squared deviations from the running mean.</param>
        /// <param name="squaredDifferenceSum">Sum of squared successive energy differences.</param>
        /// <returns>E-BFMI, or <see cref="double.NaN"/> when fewer than two energies or zero energy variance are available.</returns>
        /// <remarks>
        /// The formula is <c>mean(diff(E)^2) / var(E)</c>, with the same sample-variance
        /// convention used by Stan interfaces.
        /// </remarks>
        internal static double ComputeEnergyBayesianFractionOfMissingInformation(
            int sampleCount,
            double energyM2,
            double squaredDifferenceSum)
        {
            if (sampleCount < 2 || !(energyM2 > 0d))
                return double.NaN;

            double meanSquaredDifference = squaredDifferenceSum / sampleCount;
            double sampleVariance = energyM2 / (sampleCount - 1d);
            return meanSquaredDifference / sampleVariance;
        }

        /// <summary>
        /// Accumulates sample statistics using Welford's online algorithm for computing
        /// the diagonal mass matrix during warmup.
        /// </summary>
        /// <param name="chainIndex">The chain index.</param>
        /// <param name="sample">The parameter values from the current iteration.</param>
        private void AccumulateWelfordStatistics(int chainIndex, double[] sample)
        {
            _welfordCount[chainIndex]++;
            int n = _welfordCount[chainIndex];
            for (int j = 0; j < NumberOfParameters; j++)
            {
                double delta = sample[j] - _welfordMean[chainIndex][j];
                _welfordMean[chainIndex][j] += delta / n;
                double delta2 = sample[j] - _welfordMean[chainIndex][j];
                _welfordM2[chainIndex][j] += delta * delta2;
            }
        }

        /// <summary>
        /// The smallest number of draws an adaptation window must contain before its Welford variance
        /// is trusted. Shorter windows use the prior-scaled fallback variance instead.
        /// </summary>
        private const int MIN_ADAPT_WINDOW_COUNT = 10;

        /// <summary>
        /// The smallest per-coordinate variance retained in an adapted metric, as a fraction of the
        /// largest <i>measured</i> variance in the same window; fallback values never set that scale.
        /// This bounds the ratio of the largest measured window variance to any retained variance at
        /// 1e12; fallback values are outside that ratio.
        /// </summary>
        private const double RELATIVE_VARIANCE_FLOOR = 1e-12;

        /// <summary>
        /// Updates the diagonal mass matrix from the accumulated Welford statistics at the end
        /// of an adaptation window. Resets Welford accumulators and dual averaging state so the
        /// step size can re-adapt to the new mass matrix.
        /// </summary>
        /// <param name="chainIndex">The chain index.</param>
        /// <param name="currentState">The current parameter state, used to find a new reasonable step size after the metric change.</param>
        /// <remarks>
        /// <para>
        /// The estimated posterior variance is the <i>inverse</i> mass, not the mass. This class draws
        /// momentum with standard deviation <c>sqrt(M)</c>, evaluates kinetic energy as
        /// <c>0.5 * r' * M^-1 * r</c>, and flows position as <c>q += M^-1 * r * epsilon</c>. A coordinate
        /// with posterior standard deviation <c>s</c> therefore needs <c>M = 1 / s^2</c> for its leapfrog
        /// step to scale like <c>s</c>. The window variance is stored in <see cref="_inverseMassMatrix"/>
        /// and its reciprocal in <see cref="_massMatrix"/>, the same correspondence Stan uses when it
        /// keeps the estimated variance in <c>inv_e_metric</c>.
        /// </para>
        /// <para>
        /// Stan shrinks the window variance toward a small absolute constant because it works in an
        /// unconstrained space where the variance is O(1). This sampler works on the natural scale, where
        /// no absolute constant is meaningful, so the window variance is used directly whenever it is
        /// usable and the prior-scaled fallback <c>(priorRange / 6)^2</c> is engaged only when it is not.
        /// An unconditional blend would impose a floor tied to the prior width rather than the posterior,
        /// which on a diffuse prior flattens the metric to isotropy and erases the adaptation.
        /// </para>
        /// <para>
        /// Two conditions decide that a window cannot produce a usable variance. First, the window must
        /// contain at least <see cref="MIN_ADAPT_WINDOW_COUNT"/> draws: the relative standard error of a
        /// sample variance is approximately <c>sqrt(2 / (n - 1))</c>, still about 47% at ten draws, so a
        /// shorter window is noise rather than an estimate. Second, the variance must be finite and
        /// strictly positive. Windows failing either condition take the fallback. With the shipped buffer
        /// sizes the first condition is only reachable at a total warmup of twelve transitions or fewer.
        /// </para>
        /// <para>
        /// Retained variances are then floored at <see cref="RELATIVE_VARIANCE_FLOOR"/> times the largest
        /// <i>measured</i> variance in the same window; fallback values never set that scale, because
        /// letting them do so would put the prior range back into the floor for every other coordinate.
        /// The floor bounds the ratio of the largest <i>measured</i> window variance to any retained
        /// variance at 1e12 — fallback values are outside that ratio — and prevents a coordinate that
        /// is numerically degenerate over the window from producing an unbounded mass. It does not
        /// correct a coordinate that merely under-explored: a variance that comes back at 1e-4 to 1e-6 of
        /// the truth is far above the floor and passes through, yielding a mass that is too large and a
        /// step that is too small until the next window re-estimates it. A floor tight enough to catch
        /// that case would have to encode an expectation about how well the window mixed, which is the
        /// prior-scaled assumption this method avoids.
        /// </para>
        /// </remarks>
        private void UpdateMassMatrix(int chainIndex, ParameterSet currentState)
        {
            int n = _welfordCount[chainIndex];
            if (n < 2) return;

            // Estimate the posterior variance of every coordinate from this window, falling back to the
            // prior-scaled variance only for coordinates whose window estimate is unusable.
            var estimatedVariance = new double[NumberOfParameters];
            double largestWindowVariance = 0d;
            for (int j = 0; j < NumberOfParameters; j++)
            {
                double variance = _welfordM2[chainIndex][j] / (n - 1);

                // Fallback: (prior_range / 6)^2 as a conservative variance estimate.
                double priorRange = _upperBounds[j] - _lowerBounds[j];
                double fallbackVariance = (priorRange * priorRange) / 36.0;
                if (!Tools.IsFinite(fallbackVariance) || fallbackVariance <= 0)
                    fallbackVariance = 1.0;

                bool windowIsUsable = n >= MIN_ADAPT_WINDOW_COUNT && Tools.IsFinite(variance) && variance > 0;
                estimatedVariance[j] = windowIsUsable ? variance : fallbackVariance;

                // Only measured variances set the window's scale. Letting a fallback set it would put
                // the prior range back into the floor for every other coordinate.
                if (windowIsUsable && variance > largestWindowVariance)
                    largestWindowVariance = variance;
            }

            // Bound the metric's condition number against the window's own scale. When no coordinate
            // produced a usable variance there is no such scale, and the floor is inert.
            double varianceFloor = largestWindowVariance * RELATIVE_VARIANCE_FLOOR;
            for (int j = 0; j < NumberOfParameters; j++)
            {
                double inverseMass = Math.Max(estimatedVariance[j], varianceFloor);
                _inverseMassMatrix[chainIndex][j] = inverseMass;
                _massMatrix[chainIndex][j] = 1.0 / inverseMass;
            }

            // Reset Welford accumulators for next window
            Array.Clear(_welfordMean[chainIndex], 0, NumberOfParameters);
            Array.Clear(_welfordM2[chainIndex], 0, NumberOfParameters);
            _welfordCount[chainIndex] = 0;

            // Critical: after metric change, find a new reasonable step size (Stan's init_stepsize)
            // then reset dual averaging to re-adapt from the new starting point
            double eps0;
            try
            {
                eps0 = FindReasonableEpsilon(currentState.Values, currentState.Fitness, chainIndex);
            }
            catch
            {
                eps0 = _chainStepSizes[chainIndex];
            }
            _chainStepSizes[chainIndex] = eps0;
            _chainHBar[chainIndex] = 0.0;
            _chainAdaptStep[chainIndex] = 0;
            _chainMu[chainIndex] = Math.Log(10.0 * eps0);
            _chainLogEpsBar[chainIndex] = Math.Log(eps0);
        }

        /// <summary>
        /// Checks whether the given sample number falls at the end of a mass matrix adaptation window.
        /// </summary>
        /// <param name="sampleNum">The current sample number within warmup.</param>
        /// <returns>True if this sample marks the end of an adaptation window.</returns>
        private bool IsEndOfAdaptationWindow(int sampleNum)
        {
            for (int i = 0; i < _adaptWindowEnds.Length; i++)
            {
                if (sampleNum == _adaptWindowEnds[i])
                    return true;
            }
            return false;
        }

        /// <summary>
        /// Recursively builds a balanced binary tree of leapfrog states.
        /// </summary>
        /// <param name="theta">Starting position.</param>
        /// <param name="momentum">Starting momentum.</param>
        /// <param name="epsilon">Signed step size (negative = backward direction).</param>
        /// <param name="depth">Current tree depth (0 = single leapfrog step).</param>
        /// <param name="H0">Initial Hamiltonian for the trajectory.</param>
        /// <param name="chainIndex">The chain index for RNG and mass matrix access.</param>
        /// <returns>The tree state containing endpoints, candidate, and diagnostics.</returns>
        private TreeState BuildTree(Vector theta, Vector momentum, double epsilon, int depth, double H0, int chainIndex)
        {
            if (depth == 0)
            {
                // Base case: take one leapfrog step
                Vector thetaPrime;
                Vector momentumPrime;
                try
                {
                    (thetaPrime, momentumPrime) = Leapfrog(theta, momentum, epsilon, chainIndex);
                }
                catch (ArithmeticException)
                {
                    return InvalidTreeState(theta, momentum);
                }

                double logLH = SafeLogLikelihood(thetaPrime.Array);
                double H = -logLH + 0.5 * DiagonalQuadraticFormVec(momentumPrime, _inverseMassMatrix[chainIndex]);
                double logWeight = -H;
                if (!Tools.IsFinite(logLH) || !Tools.IsFinite(H) || !Tools.IsFinite(logWeight))
                    return InvalidTreeState(thetaPrime, momentumPrime);

                bool divergent = (H - H0) > MAX_DELTA_H;
                double alpha = Math.Min(1.0, Math.Exp(H0 - H));
                if (double.IsNaN(alpha)) alpha = 0;

                return new TreeState
                {
                    ThetaMinus = thetaPrime,
                    MomentumMinus = momentumPrime,
                    ThetaPlus = thetaPrime,
                    MomentumPlus = momentumPrime,
                    ThetaPrime = thetaPrime,
                    LogSumWeight = logWeight,
                    LogLikelihoodPrime = logLH,
                    LeafCount = 1,
                    Valid = !divergent,
                    Divergent = divergent,
                    SumAlpha = alpha,
                    NumAlpha = 1
                };
            }

            // Recursive case: build first half-tree
            var tree = BuildTree(theta, momentum, epsilon, depth - 1, H0, chainIndex);

            if (tree.Valid)
            {
                // Build second half-tree
                TreeState tree2;
                if (epsilon > 0)
                {
                    tree2 = BuildTree(tree.ThetaPlus, tree.MomentumPlus, epsilon, depth - 1, H0, chainIndex);
                    tree.ThetaPlus = tree2.ThetaPlus;
                    tree.MomentumPlus = tree2.MomentumPlus;
                }
                else
                {
                    tree2 = BuildTree(tree.ThetaMinus, tree.MomentumMinus, epsilon, depth - 1, H0, chainIndex);
                    tree.ThetaMinus = tree2.ThetaMinus;
                    tree.MomentumMinus = tree2.MomentumMinus;
                }

                // Multinomial sampling: accept candidate from tree2 with appropriate probability
                double logSumWeightNew = LogSumExp(tree.LogSumWeight, tree2.LogSumWeight);
                double acceptTree2Prob = Math.Exp(tree2.LogSumWeight - logSumWeightNew);
                if (_chainPRNGs[chainIndex].NextDouble() < acceptTree2Prob)
                {
                    tree.ThetaPrime = tree2.ThetaPrime;
                    tree.LogLikelihoodPrime = tree2.LogLikelihoodPrime;
                }

                tree.LogSumWeight = logSumWeightNew;
                tree.LeafCount += tree2.LeafCount;
                tree.Divergent |= tree2.Divergent;
                tree.SumAlpha += tree2.SumAlpha;
                tree.NumAlpha += tree2.NumAlpha;

                // Check U-turn criterion on the combined tree
                var dTheta = tree.ThetaPlus - tree.ThetaMinus;
                bool uturn = Vector.DotProduct(dTheta, tree.MomentumMinus) < 0 ||
                             Vector.DotProduct(dTheta, tree.MomentumPlus) < 0;
                tree.Valid = tree2.Valid && !uturn;
            }

            return tree;
        }

        /// <summary>
        /// Creates an invalid tree state for a trajectory that entered a non-finite log-density region.
        /// </summary>
        /// <param name="theta">The position to preserve as the tree endpoint.</param>
        /// <param name="momentum">The momentum to preserve as the tree endpoint.</param>
        /// <returns>An invalid tree state with zero acceptance contribution.</returns>
        private static TreeState InvalidTreeState(Vector theta, Vector momentum)
        {
            return new TreeState
            {
                ThetaMinus = theta.Clone(),
                MomentumMinus = momentum.Clone(),
                ThetaPlus = theta.Clone(),
                MomentumPlus = momentum.Clone(),
                ThetaPrime = theta.Clone(),
                LogSumWeight = double.NegativeInfinity,
                LogLikelihoodPrime = double.NegativeInfinity,
                LeafCount = 1,
                Valid = false,
                Divergent = true,
                SumAlpha = 0d,
                NumAlpha = 1
            };
        }

        /// <summary>
        /// Performs a single leapfrog integration step with boundary enforcement,
        /// using the per-chain diagonal mass matrix.
        /// </summary>
        /// <param name="theta">Current position.</param>
        /// <param name="momentum">Current momentum.</param>
        /// <param name="epsilon">The step size (signed: positive = forward, negative = backward).</param>
        /// <param name="chainIndex">The chain index for mass matrix access.</param>
        /// <returns>The updated (position, momentum) after one leapfrog step.</returns>
        private (Vector theta, Vector momentum) Leapfrog(Vector theta, Vector momentum, double epsilon, int chainIndex)
        {
            int D = NumberOfParameters;
            double[] invMass = _inverseMassMatrix[chainIndex];

            // Half-step momentum update. new Vector(double[]) wraps without copying, so grad is a live
            // view onto the memo's own buffer and must be read only. Vector's arithmetic operators all
            // allocate their result, so the expression below cannot write through it, but Vector's
            // indexer setter is public: do not accumulate, clamp, or negate in place through grad.
            var grad = new Vector(EvaluateGradient(theta.Array, chainIndex));
            var r = momentum + grad * (epsilon * 0.5);

            // Full-step position update using per-chain inverse mass matrix
            var q = new Vector(D);
            for (int j = 0; j < D; j++)
                q[j] = theta[j] + invMass[j] * r[j] * epsilon;

            // Enforce parameter bounds
            for (int j = 0; j < D; j++)
            {
                if (q[j] < PriorDistributions[j].Minimum)
                    q[j] = PriorDistributions[j].Minimum + Tools.DoubleMachineEpsilon;
                if (q[j] > PriorDistributions[j].Maximum)
                    q[j] = PriorDistributions[j].Maximum - Tools.DoubleMachineEpsilon;
            }

            // Half-step momentum update. As above, grad wraps the memo's buffer and must be read only.
            grad = new Vector(EvaluateGradient(q.Array, chainIndex));
            r = r + grad * (epsilon * 0.5);

            return (q, r);
        }

        /// <summary>
        /// Computes the diagonal quadratic form φᵀ M⁻¹ φ using a Vector and raw array.
        /// </summary>
        /// <param name="phi">The momentum Vector.</param>
        /// <param name="inverseMass">The diagonal inverse mass matrix as a raw array.</param>
        /// <returns>The scalar result Σᵢ phi[i]² × inverseMass[i].</returns>
        private static double DiagonalQuadraticFormVec(Vector phi, double[] inverseMass)
        {
            double sum = 0;
            for (int j = 0; j < phi.Length; j++)
                sum += phi[j] * phi[j] * inverseMass[j];
            return sum;
        }

        /// <summary>
        /// Updates the step size using the dual averaging scheme from Hoffman and Gelman (2014), Algorithm 5.
        /// </summary>
        /// <param name="chainIndex">The chain index.</param>
        /// <param name="avgAcceptProb">The average Metropolis acceptance probability from the current tree.</param>
        private void DualAveragingUpdate(int chainIndex, double avgAcceptProb)
        {
            _chainAdaptStep[chainIndex]++;
            int m = _chainAdaptStep[chainIndex];

            // Update running average of the acceptance statistic
            _chainHBar[chainIndex] = (1.0 - 1.0 / (m + T0)) * _chainHBar[chainIndex]
                                   + (TargetAcceptanceRate - avgAcceptProb) / (m + T0);

            // Compute new log step size
            double logEps = _chainMu[chainIndex] - Math.Sqrt(m) / GAMMA * _chainHBar[chainIndex];

            // Update smoothed log step size (exponential moving average)
            double mPow = Math.Pow(m, -KAPPA);
            _chainLogEpsBar[chainIndex] = mPow * logEps + (1.0 - mPow) * _chainLogEpsBar[chainIndex];

            // Set current step size (during adaptation, use the un-smoothed value)
            _chainStepSizes[chainIndex] = Math.Exp(logEps);

            // Clamp step size to prevent extreme values
            if (_chainStepSizes[chainIndex] < 1e-10)
                _chainStepSizes[chainIndex] = 1e-10;
            if (_chainStepSizes[chainIndex] > 1e5)
                _chainStepSizes[chainIndex] = 1e5;
        }

        /// <summary>
        /// Computes log(exp(a) + exp(b)) in a numerically stable way.
        /// </summary>
        private static double LogSumExp(double a, double b)
        {
            double max = Math.Max(a, b);
            if (double.IsNegativeInfinity(max)) return double.NegativeInfinity;
            return max + Math.Log(Math.Exp(a - max) + Math.Exp(b - max));
        }

        /// <summary>
        /// Evaluates the log-likelihood, returning negative infinity if the parameters are out of range.
        /// This prevents ArgumentOutOfRangeException from propagating during leapfrog integration
        /// when the sampler explores parameter values that violate distribution constraints.
        /// </summary>
        private double SafeLogLikelihood(double[] parameters)
        {
            try
            {
                return LogLikelihoodFunction(parameters);
            }
            catch (ArgumentOutOfRangeException)
            {
                return double.NegativeInfinity;
            }
        }

        /// <summary>
        /// Internal state of a binary tree node used during the NUTS tree-building recursion.
        /// </summary>
        private struct TreeState
        {
            /// <summary>Leftmost position in the subtree.</summary>
            public Vector ThetaMinus;
            /// <summary>Leftmost momentum in the subtree.</summary>
            public Vector MomentumMinus;
            /// <summary>Rightmost position in the subtree.</summary>
            public Vector ThetaPlus;
            /// <summary>Rightmost momentum in the subtree.</summary>
            public Vector MomentumPlus;
            /// <summary>Candidate position selected by multinomial sampling.</summary>
            public Vector ThetaPrime;
            /// <summary>Log of the sum of weights for multinomial sampling.</summary>
            public double LogSumWeight;
            /// <summary>Log-likelihood of the candidate position.</summary>
            public double LogLikelihoodPrime;
            /// <summary>Number of leaf nodes in the subtree.</summary>
            public int LeafCount;
            /// <summary>Whether the subtree is valid (no divergence, no U-turn).</summary>
            public bool Valid;
            /// <summary>Whether the subtree contains a divergent or non-finite trajectory.</summary>
            public bool Divergent;
            /// <summary>Sum of per-leaf Metropolis acceptance probabilities (for dual averaging).</summary>
            public double SumAlpha;
            /// <summary>Number of leaves contributing to SumAlpha.</summary>
            public int NumAlpha;
        }

    }
}
