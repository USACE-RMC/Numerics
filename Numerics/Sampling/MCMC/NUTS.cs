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
    /// During the warmup phase, the leapfrog step size is automatically adapted using dual averaging
    /// to achieve a target Metropolis acceptance probability of approximately 80%. A diagonal mass matrix
    /// is estimated during warmup using Stan-style windowed adaptation (Welford's online algorithm) to
    /// precondition the Hamiltonian dynamics for multi-scale posteriors.
    /// </para>
    /// <para>
    /// Key applications include:
    /// <list type="bullet">
    /// <item><description>Bayesian flood frequency analysis in RMC-BestFit where manual HMC tuning is impractical.</description></item>
    /// <item><description>Bayesian parameter estimation for complex hydrologic models in RFA.</description></item>
    /// <item><description>Posterior inference for mixture models and hierarchical models in TotalRisk.</description></item>
    /// </list>
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

        // Dual averaging hyperparameters (Hoffman & Gelman 2014, Section 3.2)
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
        /// The mass vector for the momentum distribution.
        /// </summary>
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
        /// The target Metropolis acceptance probability for dual averaging adaptation. Default = 0.80.
        /// </summary>
        public double TargetAcceptanceRate => DELTA_TARGET;

        /// <summary>
        /// Gets or sets whether to adapt the diagonal mass matrix during warmup.
        /// When enabled, uses Stan-style windowed adaptation with Welford's online algorithm
        /// to estimate the posterior variance per parameter and precondition the Hamiltonian dynamics.
        /// Default = false.
        /// </summary>
        public bool AdaptMassMatrix { get; set; } = false;

        /// <inheritdoc/>
        protected override void ValidateCustomSettings()
        {
            if (Mass.Length != NumberOfParameters) throw new ArgumentException(nameof(Mass), "The mass vector must be the same length as the number of parameters.");
            if (_initialStepSize <= 0) throw new ArgumentException("stepSize", "The leapfrog step size must be positive.");
            if (MaxTreeDepth < 1) throw new ArgumentException(nameof(MaxTreeDepth), "The maximum tree depth must be at least 1.");
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

            // Initialize diagonal mass matrix and Welford accumulators
            _welfordMean = new double[N][];
            _welfordM2 = new double[N][];
            _welfordCount = new int[N];
            _massMatrix = new double[N][];
            _inverseMassMatrix = new double[N][];

            for (int i = 0; i < N; i++)
            {
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

            // Compute adaptation window boundaries following Stan exactly.
            // Stan defaults: init_buffer=75, term_buffer=50, base_window=25
            // Phase 1 (init_buffer): step size adaptation only, identity mass matrix
            // Phase 2 (slow adaptation): mass matrix + step size in doubling windows
            // Phase 3 (term_buffer): final step size tuning with fixed mass matrix
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

            // Take one leapfrog step
            var thetaPrime = (double[])theta0.Clone();
            var rPrime = (double[])r0.Clone();
            LeapfrogInPlace(thetaPrime, rPrime, epsilon, chainIndex);

            double logLH1 = SafeLogLikelihood(thetaPrime);
            double H1 = -logLH1 + 0.5 * DiagonalQuadraticForm(rPrime, _inverseMassMatrix[chainIndex]);
            double logAlpha = H0 - H1;
            if (double.IsNaN(logAlpha) || double.IsNegativeInfinity(logAlpha))
                logAlpha = -1000.0;

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

                // Recompute with new epsilon
                Array.Copy(theta0, thetaPrime, D);
                Array.Copy(r0, rPrime, D);
                LeapfrogInPlace(thetaPrime, rPrime, epsilon, chainIndex);

                logLH1 = SafeLogLikelihood(thetaPrime);
                H1 = -logLH1 + 0.5 * DiagonalQuadraticForm(rPrime, _inverseMassMatrix[chainIndex]);
                logAlpha = H0 - H1;
                if (double.IsNaN(logAlpha) || double.IsNegativeInfinity(logAlpha))
                    logAlpha = -1000.0;
            }

            return Math.Max(1e-8, Math.Min(epsilon, 1e6));
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

            // Half-step momentum update
            double[] grad = NumericalDerivative.Gradient((y) => SafeLogLikelihood(y), theta, _lowerBounds, _upperBounds);
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

            // Half-step momentum update
            grad = NumericalDerivative.Gradient((y) => SafeLogLikelihood(y), theta, _lowerBounds, _upperBounds);
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

            double eps = _chainStepSizes[index];
            int D = NumberOfParameters;

            // Step 1: Sample momentum from N(0, M) using per-chain mass matrix
            var phi = new Vector(D);
            for (int i = 0; i < D; i++)
                phi[i] = Math.Sqrt(_massMatrix[index][i]) * Normal.StandardZ(_chainPRNGs[index].NextDouble());

            // Compute initial Hamiltonian using per-chain inverse mass matrix
            double H0 = -state.Fitness + 0.5 * DiagonalQuadraticFormVec(phi, _inverseMassMatrix[index]);

            // Step 2: Initialize tree
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

            // Step 3: Build tree by doubling until U-turn or max depth
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

            // Step 4: Warmup adaptation (step size + mass matrix)
            int warmupSteps = WarmupIterations * ThinningInterval;
            int sampleNum = SampleCount[index];

            if (sampleNum <= warmupSteps)
            {
                // Always do dual averaging step size adaptation during warmup
                double avgAlpha = numAlpha > 0 ? sumAlpha / numAlpha : DELTA_TARGET;
                DualAveragingUpdate(index, avgAlpha);

                // Accumulate Welford statistics during mass matrix adaptation windows (Phase 2)
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

            // NUTS always accepts
            AcceptCount[index] += 1;
            return new ParameterSet(candidate.Array, candidateLogLH);
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
        /// Updates the diagonal mass matrix from the accumulated Welford statistics at the end
        /// of an adaptation window. Resets Welford accumulators and dual averaging state so the
        /// step size can re-adapt to the new mass matrix.
        /// </summary>
        /// <param name="chainIndex">The chain index.</param>
        /// <param name="currentState">The current parameter state, used to find a new reasonable step size after the metric change.</param>
        private void UpdateMassMatrix(int chainIndex, ParameterSet currentState)
        {
            int n = _welfordCount[chainIndex];
            if (n < 2) return;

            for (int j = 0; j < NumberOfParameters; j++)
            {
                double variance = _welfordM2[chainIndex][j] / (n - 1);
                // Stan regularization: (n/(n+5)) * var + 1e-3 * (5/(n+5))
                // Stan operates in unconstrained space where variance ~ O(1), so 1e-3 is fine.
                // We operate in natural scale, so use a scale-aware fallback instead.
                // Fallback: (prior_range / 6)^2 as a conservative variance estimate.
                double priorRange = _upperBounds[j] - _lowerBounds[j];
                double fallbackVariance = (priorRange * priorRange) / 36.0;
                if (!Tools.IsFinite(fallbackVariance) || fallbackVariance <= 0)
                    fallbackVariance = 1.0;
                double shrinkage = 5.0;
                double regularized = (n / (n + shrinkage)) * variance + (shrinkage / (n + shrinkage)) * fallbackVariance;
                _massMatrix[chainIndex][j] = regularized;
                _inverseMassMatrix[chainIndex][j] = 1.0 / regularized;
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
                var (thetaPrime, momentumPrime) = Leapfrog(theta, momentum, epsilon, chainIndex);
                double logLH = SafeLogLikelihood(thetaPrime.Array);
                double H = -logLH + 0.5 * DiagonalQuadraticFormVec(momentumPrime, _inverseMassMatrix[chainIndex]);
                double logWeight = -H;
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

            // Half-step momentum update
            var grad = GradientFunction(theta.Array);
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

            // Half-step momentum update
            grad = GradientFunction(q.Array);
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
                                   + (DELTA_TARGET - avgAcceptProb) / (m + T0);

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
            /// <summary>Sum of per-leaf Metropolis acceptance probabilities (for dual averaging).</summary>
            public double SumAlpha;
            /// <summary>Number of leaves contributing to SumAlpha.</summary>
            public int NumAlpha;
        }

    }
}
