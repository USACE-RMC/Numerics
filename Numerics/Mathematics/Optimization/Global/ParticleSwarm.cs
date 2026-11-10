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
        public double[] LowerBounds { get; private set; } = null!;

        /// <summary>
        /// An array of upper bounds (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double[] UpperBounds { get; private set; } = null!;

        /// <inheritdoc />
        protected override double[]? ParameterLowerBounds => LowerBounds;

        /// <inheritdoc />
        protected override double[]? ParameterUpperBounds => UpperBounds;

        /// <summary>
        /// The total population size. Default = 30.
        /// </summary>
        public int PopulationSize { get; set; } = 30;

        /// <summary>
        /// Determines whether the best point is refined by a local search when the global search ends. Default = false.
        /// </summary>
        /// <remarks>
        /// The polish runs after convergence or after <see cref="Optimizer.MaxIterations"/>, never after the
        /// evaluation budget is spent. Its evaluations count toward <see cref="Optimizer.FunctionEvaluations"/>
        /// and the reported status is the global search's.
        /// </remarks>
        public bool Polish { get; set; }

        /// <summary>
        /// The local method used to polish: BFGS, Nelder-Mead, or Powell. Default = Nelder-Mead.
        /// </summary>
        public LocalMethod PolishMethod { get; set; } = LocalMethod.NelderMead;

        /// <summary>
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// Determines whether each iteration's particles are evaluated in parallel. Default = false.
        /// </summary>
        /// <remarks>
        /// When true, the swarm moves synchronously: every particle is moved toward the global best of the
        /// previous iteration, the particles are evaluated concurrently, and the personal and global bests are
        /// then updated in index order. The random draws are made in the same order on one thread, so the
        /// result is reproducible and does not depend on the number of threads. It differs from the default
        /// asynchronous update, where a particle sees a global best improved earlier in the same iteration.
        /// The objective function must be thread safe.
        /// </remarks>
        public bool UseParallel { get; set; }

        /// <summary>
        /// The maximum number of concurrent evaluations when <see cref="UseParallel"/> is true, or -1 for no
        /// limit. Default = -1. The result does not depend on it.
        /// </summary>
        public int MaxDegreeOfParallelism { get; set; } = -1;

        /// <inheritdoc/>
        /// <remarks>
        /// <b> References: </b>
        /// Alam "Particle Swarm Optimization: Algorithm and its Codes in MATLAB" (2016)
        /// </remarks>
        protected override void Optimize()
        {
            if (PopulationSize < 1) throw new ArgumentOutOfRangeException(nameof(PopulationSize), "The population size must be greater than 0.");

            if (Polish) ValidatePolishMethod(PolishMethod);
            if (MaxDegreeOfParallelism == 0 || MaxDegreeOfParallelism < -1) throw new ArgumentOutOfRangeException(nameof(MaxDegreeOfParallelism), "The maximum degree of parallelism must be positive or -1.");
            if (UseParallel)
            {
                OptimizeSynchronous();
                return;
            }
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
                    Xp[i].Velocity = v;

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
                    if (Polish && !PolishBestParameterSet(PolishMethod)) return;
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }

                Iterations += 1;

            }

            // Reaching this point means the maximum iterations were reached.
            if (Polish && !PolishBestParameterSet(PolishMethod)) return;
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

        /// <summary>
        /// The synchronous particle swarm with parallel evaluation of each iteration. See <see cref="UseParallel"/>.
        /// </summary>
        private void OptimizeSynchronous()
        {
            int i, j, D = NumberOfParameters, N = PopulationSize;
            bool cancel = false;
            double wmin = 0.4, wmax = 0.9, c1 = 2.05, c2 = 2.05;

            var prng = new MersenneTwister(PRNGSeed);
            var rnd = LatinHypercube.Random(N, D, PRNGSeed);
            var positions = new double[N][];
            var velocities = new double[N][];
            var fitness = new double[N];
            for (i = 0; i < N; i++)
            {
                positions[i] = new double[D];
                velocities[i] = new double[D];
                for (j = 0; j < D; j++)
                {
                    positions[i][j] = LowerBounds[j] + rnd[i, j] * (UpperBounds[j] - LowerBounds[j]);
                    velocities[i][j] = 0.1 * (UpperBounds[j] - LowerBounds[j]) * (prng.NextDouble() - 0.5);
                }
            }
            int done = EvaluateBatch(positions, N, fitness, MaxDegreeOfParallelism, ref cancel);
            if (cancel) return;
            var personal = new ParameterSet[N];
            for (i = 0; i < N; i++)
                personal[i] = new ParameterSet((double[])positions[i].Clone(), fitness[i]);

            Iterations += 1;
            var global = new double[D];

            while (Iterations < MaxIterations)
            {
                Array.Copy(BestParameterSet.Values, global, D);
                var w = wmax - (wmax - wmin) * Iterations / MaxIterations;
                for (i = 0; i < N; i++)
                {
                    for (j = 0; j < D; j++)
                    {
                        velocities[i][j] = w * velocities[i][j] + c1 * prng.NextDouble() * (personal[i].Values[j] - positions[i][j])
                                                                + c2 * prng.NextDouble() * (global[j] - positions[i][j]);
                        positions[i][j] = RepairParameter(positions[i][j] + velocities[i][j], LowerBounds[j], UpperBounds[j]);
                    }
                }

                done = EvaluateBatch(positions, N, fitness, MaxDegreeOfParallelism, ref cancel);
                var statistics = new RunningStatistics();
                for (i = 0; i < done; i++)
                {
                    if (fitness[i] <= personal[i].Fitness)
                        personal[i] = new ParameterSet((double[])positions[i].Clone(), fitness[i]);
                    statistics.Push(personal[i].Fitness);
                }
                if (cancel) return;

                if (Iterations >= 10 && statistics.StandardDeviation < AbsoluteTolerance + RelativeTolerance * Math.Abs(statistics.Mean))
                {
                    if (Polish && !PolishBestParameterSet(PolishMethod)) return;
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }

                Iterations += 1;
            }

            if (Polish && !PolishBestParameterSet(PolishMethod)) return;
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
            public double[] Velocity { get; set; } = null!;
        }
    }
}
