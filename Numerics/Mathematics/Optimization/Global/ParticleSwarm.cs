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

using Numerics.Data.Statistics;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// The Particle Swarm optimization algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This method optimizes a problem by iteratively trying to improve a candidate solution with a given measure of quality.
    /// Using a population of candidate solutions, called particles, these particles move around in the search space. Each 
    /// particles's movement is influenced by its local best known position. These positions are updated as better positions are 
    /// found by other particles, and the swarm is expected to move toward the best solutions.
    /// </para>
    /// <para>
    ///     <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Implements routine described by Kockenderfer and Wheeler. "Algorithms for Optimization" (2019).
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Particle_swarm_optimization"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class ParticleSwarm : Optimizer
    {
        /// <summary>
        /// Construct a new particle swarm optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public ParticleSwarm(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
        {
            // Check if the length of the lower and upper bounds equal the number of parameters
            if (lowerBounds.Count != numberOfParameters || upperBounds.Count != numberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(lowerBounds), "The lower and upper bounds must be the same length as the number of parameters.");
            }
            // Check if the lower bounds are less than the upper bounds
            for (int i = 0; i < lowerBounds.Count; i++)
            {
                if (upperBounds[i] <= lowerBounds[i])
                {
                    throw new ArgumentOutOfRangeException(nameof(upperBounds), "The upper bound cannot be less than or equal to the lower bound.");
                }
            }
            LowerBounds = lowerBounds.ToArray();
            UpperBounds = upperBounds.ToArray();
        }

        /// <summary>
        /// An array of lower bounds (inclusive) of the interval containing the optimal point. 
        /// </summary>
        public double[] LowerBounds { get; private set; }

        /// <summary>
        /// An array of upper bounds (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double[] UpperBounds { get; private set; }

        /// <summary>
        /// The total population size. Default = 30.
        /// </summary>
        public int PopulationSize { get; set; } = 30;

        /// <summary>
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <inheritdoc/>
        /// <remarks>
        /// <b> References: </b>
        /// Alam "Particle Swarm Optimization: Algorithm and its Codes in MATLAB" (2016)
        /// </remarks>
        protected override void Optimize()
        {
            if (PopulationSize < 1) throw new ArgumentOutOfRangeException(nameof(PopulationSize), "The population size must be greater than 0.");

            int i, j, D = NumberOfParameters;
            bool cancel = false;

            // These parameters are recommended by Alam "Particle Swarm Optimization: Algorithm and its Codes in MATLAB" (2016)
            double wmin = 0.4;
            double wmax = 0.9;
            double c1 = 2.05;
            double c2 = 2.05;

            // Initialize the population of points
            var prng = new MersenneTwister(PRNGSeed);
            var rnd = LatinHypercube.Random(PopulationSize, D, PRNGSeed);
            var Xp = new List<Particle>();
            for (i = 0; i < PopulationSize; i++)
            {
                var values = new double[D];
                var velocity = new double[D];
                for (j = 0; j < D; j++)
                {
                    values[j] = LowerBounds[j] + rnd[i, j] * (UpperBounds[j] - LowerBounds[j]);
                    velocity[j] = 0.1 * (UpperBounds[j] - LowerBounds[j]) * (prng.NextDouble() - 0.5);

                }

                var fitness = Evaluate(values, ref cancel);
                if (cancel) return;

                Xp.Add(new Particle { ParameterSet = new ParameterSet(values, fitness),
                                      BestParameterSet = new ParameterSet(values, fitness),
                                      Velocity = velocity});
            }

            Iterations += 1;

            // Perform Particle Swarm
            while (Iterations < MaxIterations)
            {
                // Keep track of population statistics to assess convergence
                var statistics = new RunningStatistics();

                // Update momentum for each particle
                for (i = 0; i < PopulationSize; i++)
                {
                    // Generate trial vector
                    var u = new double[D];
                    var v = new double[D];
                    for (j = 0; j < D; j++)
                    {
                        // New velocity
                        var w = wmax - (wmax - wmin) * Iterations / MaxIterations;
                        v[j] = w * Xp[i].Velocity[j] + c1 * prng.NextDouble() * (Xp[i].BestParameterSet.Values[j] - Xp[i].ParameterSet.Values[j]) 
                                                     + c2 * prng.NextDouble() * (BestParameterSet.Values[j] - Xp[i].ParameterSet.Values[j]);                      
                        // New position
                        u[j] = Xp[i].ParameterSet.Values[j] + v[j];
                        u[j] = RepairParameter(u[j], LowerBounds[j], UpperBounds[j]);
                    }

                    // Evaluate fitness
                    var fitness = Evaluate(u, ref cancel);
                    if (cancel) return;

                    Xp[i].ParameterSet = new ParameterSet(u, fitness);
                    Xp[i].Velocity = v.ToArray();

                    // Update population
                    if (fitness <= Xp[i].BestParameterSet.Fitness)
                    {
                        Xp[i].BestParameterSet = new ParameterSet(u, fitness);
                    }
             
                    // Keep running stats of population
                    statistics.Push(Xp[i].BestParameterSet.Fitness);
                }

                // Evaluate convergence 
                if (Iterations >= 10 && statistics.StandardDeviation < AbsoluteTolerance + RelativeTolerance * Math.Abs(statistics.Mean))
                {
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }

                Iterations += 1;

            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

        /// <summary>
        /// Class for storing particles. 
        /// </summary>
        private class Particle
        {
            /// <summary>
            /// Gets or sets the current parameter set (position) of the particle in the search space.
            /// </summary>
            public ParameterSet ParameterSet { get; set; }

            /// <summary>
            /// Gets or sets the best parameter set (position) this particle has found during the search.
            /// </summary>
            /// <remarks>
            /// This represents the particle's "personal best" or "pbest" in PSO terminology.
            /// </remarks>
            public ParameterSet BestParameterSet { get; set; }

            /// <summary>
            /// Gets or sets the velocity vector of the particle in parameter space.
            /// </summary>
            /// <remarks>
            /// The velocity determines how the particle moves through the search space.
            /// It is updated based on the particle's personal best and the global best.
            /// </remarks>
            public double[] Velocity { get; set; }
        }
    }
}
