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
    /// to achieve a target Metropolis acceptance probability of approximately 80%.
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
        /// <param name="mass">Optional. The mass vector for the momentum distribution. Default = Identity.</param>
        /// <param name="stepSize">Optional. The initial leapfrog step size. Will be adapted during warmup. Default = 0.1.</param>
        /// <param name="maxTreeDepth">Optional. The maximum binary tree depth. Default = 10.</param>
        /// <param name="gradientFunction">Optional. The function for evaluating the gradient of the log-likelihood. Numerical finite difference will be used by default.</param>
        public NUTS(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction, Vector? mass = null, double stepSize = 0.1, int maxTreeDepth = 10, HMC.Gradient? gradientFunction = null) : base(priorDistributions, logLikelihoodFunction)
        {
            InitialIterations = 100 * NumberOfParameters;

            // Set the mass vector
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

            // Set the gradient function
            if (gradientFunction == null)
            {
                GradientFunction = (x) => new Vector(NumericalDerivative.Gradient((y) => LogLikelihoodFunction(y), x.ToArray()));
            }
            else
            {
                GradientFunction = gradientFunction;
            }
        }

        private Vector _inverseMass;
        private double _initialStepSize;

        // Per-chain dual averaging state
        private double[] _chainStepSizes;
        private double[] _chainLogEpsBar;
        private double[] _chainHBar;
        private double[] _chainMu;
        private int[] _chainAdaptStep;

        // Dual averaging hyperparameters (Hoffman & Gelman 2014, Section 3.2)
        private const double DELTA_TARGET = 0.80;
        private const double GAMMA = 0.05;
        private const double T0 = 10.0;
        private const double KAPPA = 0.75;

        // Divergence threshold: if H - H0 exceeds this, the trajectory is considered divergent
        private const double MAX_DELTA_H = 1000.0;

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
            _chainStepSizes = new double[NumberOfChains];
            _chainLogEpsBar = new double[NumberOfChains];
            _chainHBar = new double[NumberOfChains];
            _chainMu = new double[NumberOfChains];
            _chainAdaptStep = new int[NumberOfChains];

            for (int i = 0; i < NumberOfChains; i++)
            {
                _chainStepSizes[i] = _initialStepSize;
                _chainLogEpsBar[i] = Math.Log(_initialStepSize);
                _chainHBar[i] = 0.0;
                _chainMu[i] = Math.Log(10.0 * _initialStepSize);
                _chainAdaptStep[i] = 0;
            }
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {
            // Update the sample count
            SampleCount[index] += 1;

            double eps = _chainStepSizes[index];

            // Step 1: Sample momentum from N(0, M)
            var phi = new Vector(NumberOfParameters);
            for (int i = 0; i < NumberOfParameters; i++)
                phi[i] = Math.Sqrt(Mass[i]) * Normal.StandardZ(_chainPRNGs[index].NextDouble());

            // Compute initial Hamiltonian: H = -log p(theta) + 0.5 * phi^T M^{-1} phi
            double H0 = -state.Fitness + 0.5 * HMC.QuadraticForm(phi, _inverseMass);

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

            // Step 4: Dual averaging step size adaptation during warmup
            int warmupSteps = WarmupIterations * ThinningInterval;
            if (SampleCount[index] <= warmupSteps)
            {
                double avgAlpha = numAlpha > 0 ? sumAlpha / numAlpha : DELTA_TARGET;
                DualAveragingUpdate(index, avgAlpha);
            }
            else if (SampleCount[index] == warmupSteps + 1)
            {
                // After warmup, fix step size to the smoothed value
                _chainStepSizes[index] = Math.Exp(_chainLogEpsBar[index]);
            }

            // NUTS always accepts
            AcceptCount[index] += 1;
            return new ParameterSet(candidate.Array, candidateLogLH);
        }

        /// <summary>
        /// Recursively builds a balanced binary tree of leapfrog states.
        /// </summary>
        /// <param name="theta">Starting position.</param>
        /// <param name="momentum">Starting momentum.</param>
        /// <param name="epsilon">Signed step size (negative = backward direction).</param>
        /// <param name="depth">Current tree depth (0 = single leapfrog step).</param>
        /// <param name="H0">Initial Hamiltonian for the trajectory.</param>
        /// <param name="chainIndex">The chain index for RNG access.</param>
        /// <returns>The tree state containing endpoints, candidate, and diagnostics.</returns>
        private TreeState BuildTree(Vector theta, Vector momentum, double epsilon, int depth, double H0, int chainIndex)
        {
            if (depth == 0)
            {
                // Base case: take one leapfrog step
                var (thetaPrime, momentumPrime) = Leapfrog(theta, momentum, epsilon);
                double logLH = LogLikelihoodFunction(thetaPrime.Array);
                double H = -logLH + 0.5 * HMC.QuadraticForm(momentumPrime, _inverseMass);
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
        /// Performs a single leapfrog integration step with boundary enforcement.
        /// </summary>
        /// <param name="theta">Current position.</param>
        /// <param name="momentum">Current momentum.</param>
        /// <param name="epsilon">The step size (signed: positive = forward, negative = backward).</param>
        /// <returns>The updated (position, momentum) after one leapfrog step.</returns>
        private (Vector theta, Vector momentum) Leapfrog(Vector theta, Vector momentum, double epsilon)
        {
            // Half-step momentum update
            var grad = GradientFunction(theta.Array);
            var r = momentum + grad * (epsilon * 0.5);

            // Full-step position update
            var q = theta + _inverseMass * r * epsilon;

            // Enforce parameter bounds
            for (int j = 0; j < NumberOfParameters; j++)
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
