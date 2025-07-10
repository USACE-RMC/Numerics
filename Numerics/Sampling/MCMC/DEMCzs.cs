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

namespace Numerics.Sampling.MCMC
{
    /// <summary>
    /// This class performs Bayesian MCMC using the adaptive Differential Evolution Markov Chain with snooker update (DE-MCzs) method.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    ///    <b>  References: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     This class is based on an algorithm described in:
    ///     Braak and Vrugt "Differential Evolution Markov Chain with snooker update
    ///     and fewer chains" (2008) Statistics and Computing.
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class DEMCzs : MCMCSampler
    {

        /// <summary>
        /// Constructs a new DEMCzs sampler.
        /// </summary>
        /// <param name="priorDistributions">The list of prior distributions for the model parameters.</param>
        /// <param name="logLikelihoodFunction">The Log-Likelihood function to evaluate.</param>
        public DEMCzs(List<IUnivariateDistribution> priorDistributions, LogLikelihood logLikelihoodFunction) : base(priorDistributions, logLikelihoodFunction)
        {
            InitialIterations = 100 * NumberOfParameters;
            // DE-MCz options
            IsPopulationSampler = true;
            // Jump parameter. Default = 2.38/SQRT(2*D)
            Jump = 2.38d / Math.Sqrt(2.0d * NumberOfParameters);
            // Adaptation threshold. Default = 0.1 or 10% of the time. 
            JumpThreshold = 0.1d;
            // Snooker update. Default = 0.1 or 10% of the time.
            SnookerThreshold = 0.1d;
            _b = new Normal(0, _noise);
            _g = new Uniform(1.2, 2.2);
        }

        private double _noise = 1E-3;
        private Normal _b;
        private Uniform _g;

        /// <summary>
        /// The jumping parameter used to jump from one mode region to another in the target distribution.
        /// </summary>
        public double Jump { get; set; }

        /// <summary>
        /// Determines how often the jump parameter switches to 1.0; e.g., 0.10 will result in a large jump 10% of the time.
        /// </summary>
        public double JumpThreshold { get; set; }

        /// <summary>
        /// Determines how often to perform the snooker update; e.g., 0.10 will result in an update 10% of the time. 
        /// </summary>
        public double SnookerThreshold { get; set; }

        /// <summary>
        /// The noise parameter (b).
        /// </summary>
        public double Noise
        {
            get { return _noise; }
            set
            {
                _noise = value;
                _b = new Normal(0, _noise);
            }
        }

        /// <inheritdoc/>
        protected override void ValidateCustomSettings()
        {
            if (NumberOfChains < 3) throw new ArgumentException(nameof(NumberOfChains), "There must be at least 3 chains.");
            if (Jump <= 0 || Jump >= 2) throw new ArgumentException(nameof(Jump), "The jump parameter must be between 0 and 2.");
            if (JumpThreshold < 0 || JumpThreshold > 1) throw new ArgumentException(nameof(JumpThreshold), "The jump threshold must be between 0 and 1.");
            if (SnookerThreshold < 0 || SnookerThreshold > 0.5) throw new ArgumentException(nameof(SnookerThreshold), "The snooker threshold must be between 0 and 0.5.");
            if (Noise < 0) throw new ArgumentException(nameof(Noise), "The noise parameter must be greater than 0.");
        }

        /// <inheritdoc/>
        protected override ParameterSet ChainIteration(int index, ParameterSet state)
        {
            // Update the sample count
            SampleCount[index] += 1;

            // 10% snooker updates and 90% parallel direction updates.
            if (_chainPRNGs[index].NextDouble() <= SnookerThreshold && SampleCount[index] > 5 * ThinningInterval)
            {
                return SnookerUpdate(index, state);
            }
            else
            {

                // The adaptation for the algorithm to allow for 
                // jumps from one mode region to another in the target
                // distribution.         
                var G = _chainPRNGs[index].NextDouble() <= JumpThreshold ? 1.0d : Jump;

                // Sample uniformly at random without replacement two numbers R1 and R2
                // from the numbers 1, 2, ..., M. 
                int r1, r2, M = PopulationMatrix.Count;
                r1 = _chainPRNGs[index].Next(0, M); 
                do r2 = _chainPRNGs[index].Next(0, M); while (r2 == r1);

                // Calculate the proposal vector
                // x* ← xi + γ(zR1 − zR2) + e
                // where zR1 and zR2 are rows R1 and R2 of the population matrix.
                var xp = new double[NumberOfParameters];
                for (int i = 0; i < NumberOfParameters; i++)
                {
                    var xi = state.Values[i];
                    var zr1 = PopulationMatrix[r1].Values[i];
                    var zr2 = PopulationMatrix[r2].Values[i];
                    var e = _b.InverseCDF(_chainPRNGs[index].NextDouble());
                    xp[i] = xi + G * (zr1 - zr2) + e;

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

        /// <summary>
        /// Returns a proposed MCMC iteration based on the Snooker Update method. 
        /// </summary>
        /// <param name="index">The Markov Chain zero-based index.</param>
        /// <param name="state">The current chain state to compare against.</param>
        private ParameterSet SnookerUpdate(int index, ParameterSet state)
        {
            // Do snooker update
            // Get Jump -- uniform random number between 1.2 and 2.2         
            double G = _g.InverseCDF(_chainPRNGs[index].NextDouble());
            //double G = 1.7;

            // Select another chain, which is in state z
            int c = index;
            do c = _chainPRNGs[index].Next(0, NumberOfChains); while (c == index);
            // Select two other random chains, zR1 and zR2
            int c1 = c, c2 = c;
            do c1 = _chainPRNGs[index].Next(0, NumberOfChains); while (c1 == c);
            do c2 = _chainPRNGs[index].Next(0, NumberOfChains); while (c2 == c1 || c2 == c);

            // Define z
            int n = MarkovChains[c].Count();
            var z = new Vector(MarkovChains[c][n - 1].Values);
            var xi = new Vector(state.Values);
            // Define line xi - z
            var line = xi - z;
            // Orthogonally project zR1 and zR2 onto the line xi - z
            var zr1 = new Vector(MarkovChains[c1][n - 1].Values);
            var zr2 = new Vector(MarkovChains[c2][n - 1].Values);
            var zp1 = Vector.Project(zr1, line);
            var zp2 = Vector.Project(zr2, line);

            // Calculate the proposal vector
            // x* ← xi + γ(zP1 − zP2)
            var xp = new Vector(NumberOfParameters);
            for (int i = 0; i < NumberOfParameters; i++)
            {
                xp[i] = xi[i] + G * (zp1[i] - zp2[i]);

                // Check if the parameter is feasible (within the constraints)
                if (xp[i] < PriorDistributions[i].Minimum || xp[i] > PriorDistributions[i].Maximum)
                {
                    // The proposed parameter vector was infeasible, 
                    // so leave xi unchanged.
                    return state;
                }
            }

            // Evaluate fitness
            var logLHp = LogLikelihoodFunction(xp.Array);
            var logLHi = state.Fitness;

            // Euclidean distance terms for Jacobian correction
            double distP = Vector.Distance(xp, z);
            double distI = Vector.Distance(xi, z);

            // Avoid underflow in log terms
            if (distP < 1e-12 || distI < 1e-12)
                return state;

            double logEDp = (NumberOfParameters - 1) * Math.Log(distP);
            double logEDi = (NumberOfParameters - 1) * Math.Log(distI);

            // Calculate the Metropolis ratio
            var logRatio = logLHp + logEDp - logLHi - logEDi;
            // Accept the proposal with probability min(1,r)
            // otherwise leave xi unchanged
            var logU = Math.Log(_chainPRNGs[index].NextDouble());
            if (logU <= logRatio)
            {
                // The proposal is accepted
                AcceptCount[index] += 1;
                return new ParameterSet(xp.Array, logLHp);
            }
            else
            {
                return state;
            }
        }

    }
}
