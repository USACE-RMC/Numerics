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
    /// The Differential Evolution (DE) algorithm, which finds a global minima when no gradient is available.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This method optimizes a problem by trying to improve a candidate solution with regard to a given measure
    /// of quality. DE maintains a population of candidate solutions, creating new ones by combining existing 
    /// ones, and keeping whichever solution has the best score or fitness on the optimization problem. Because 
    /// the problem is treated as merely a measure of quality for the candidate solution, the gradient is not 
    /// needed, and the problem does not need to be differentiable.
    /// </para>
    /// <para>
    ///     <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Implements routine described by Price et al. "Differential Evolution: A Practical Approach to Global Optimization" (1998).
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Differential_evolution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class DifferentialEvolution : Optimizer
    {

        /// <summary>
        /// Construct a new differential evolution optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public DifferentialEvolution(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
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

            PopulationSize = 10 * NumberOfParameters;
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
        /// The total population size. Default = 10 * D (Storn &amp; Price, 1997).
        /// </summary>
        public int PopulationSize { get; set; } = 30;

        /// <summary>
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// The mutation constant or differential weight, in the range [0, 2]. Increasing the mutation constant increases the search radius, 
        /// but will slow down convergence. The default is 0.75. 
        /// </summary>
        public double Mutation { get; set; } = 0.75;

        /// <summary>
        /// Determines how often the mutation constant dithers between 0.5 and 1.0; e.g., 0.90 will result in dithering 90% of the time.
        /// </summary>
        public double DitherRate { get; set; } = 0.9;

        /// <summary>
        /// The crossover probability or recombination constant, in the range [0, 1]. Increasing this value allows a larger number of mutants
        /// to progress into the next generation, but at the risk of population stability. 
        /// </summary>
        public double CrossoverProbability { get; set; } = 0.9;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            if (PopulationSize < 1) throw new ArgumentOutOfRangeException(nameof(PopulationSize), "The population size must be greater than 0.");
            if (Mutation < 0 || Mutation > 2) throw new ArgumentOutOfRangeException(nameof(Mutation), "The mutation parameter must be between 0 and 2.");
            if (DitherRate < 0 || DitherRate > 1) throw new ArgumentOutOfRangeException(nameof(DitherRate), "The dithering rate must be between 0 and 1.");
            if (CrossoverProbability < 0 || CrossoverProbability > 1) throw new ArgumentOutOfRangeException(nameof(CrossoverProbability), "The crossover probability must be between 0 and 1.");

            int i, j, D = NumberOfParameters;
            bool cancel = false;

            // Initialize the population of points
            var prng = new MersenneTwister(PRNGSeed);
            var rnd = LatinHypercube.Random(PopulationSize, D, PRNGSeed);
            var Xp = new List<ParameterSet>();
            for (i = 0; i < PopulationSize; i++)
            {
                var values = new double[D];
                for (j = 0; j < D; j++)
                    values[j] = LowerBounds[j] + rnd[i, j] * (UpperBounds[j] - LowerBounds[j]);
                var fitness = Evaluate(values, ref cancel);
                Xp.Add(new ParameterSet(values, fitness));
                if (cancel == true) return;
            }

            Iterations += 1;

            // Perform Differential Evolution
            while (Iterations < MaxIterations)
            {
                // Keep track of population statistics to assess convergence
                var statistics = new RunningStatistics();

                // Mutate and recombine population
                for (i = 0; i < PopulationSize; i++)
                {

                    // Randomly select three vectors indexes
                    var indices = Enumerable.Range(0, PopulationSize).Where(x => x != i).OrderBy(_ => prng.NextDouble()).Take(3).ToArray();
                    int r0 = indices[0], r1 = indices[1], r2 = indices[2];

                    // Generate trial vector
                    double G = prng.NextDouble() <= DitherRate ? 0.5 + prng.NextDouble() * 0.5 : Mutation; // Scale factor
                    double jRand = prng.Next(0, D);
                    var u = new double[D];
                    for (j = 0; j < D; j++)
                    {
                        var rr = prng.NextDouble();
                        if (rr <= CrossoverProbability || j == jRand)
                        {
                            u[j] = Xp[r0].Values[j] + G * (Xp[r1].Values[j] - Xp[r2].Values[j]);
                            u[j] = RepairParameter(u[j], LowerBounds[j], UpperBounds[j]);
                        }
                        else
                        {
                            u[j] = Xp[i].Values[j];
                        }
                    }

                    // Evaluate fitness
                    var fitness = Evaluate(u, ref cancel);
                    if (cancel) return;

                    // Update population
                    if (fitness <= Xp[i].Fitness)
                    {
                        Xp[i] = new ParameterSet(u, fitness);
                    }
              
                    // Keep running stats of population
                    statistics.Push(Xp[i].Fitness);

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
        
    }
}
