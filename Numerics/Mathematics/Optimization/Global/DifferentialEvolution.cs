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

        /// <inheritdoc />
        protected override double[]? ParameterLowerBounds => LowerBounds;

        /// <inheritdoc />
        protected override double[]? ParameterUpperBounds => UpperBounds;

        /// <summary>
        /// The total population size. Default = 10 * D (Storn &amp; Price, 1997).
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

        /// <summary>
        /// Optional starting members, for example a previous optimum or a baseline design. Default = null.
        /// </summary>
        /// <remarks>
        /// Row k replaces member k of the Latin hypercube initial population; the other members and the
        /// random stream are unchanged. Each row must have one value per parameter, inside the bounds, and
        /// there can be at most <see cref="PopulationSize"/> rows. With null the search is unchanged.
        /// </remarks>
        public IList<double[]>? InitialPopulation { get; set; }

        /// <summary>
        /// Optional flags marking integer-valued parameters, one per parameter. Default = null (all continuous).
        /// </summary>
        /// <remarks>
        /// An integer parameter is rounded to the nearest integer (ties to even) and kept within
        /// [ceiling(lower bound), floor(upper bound)] before every evaluation, as in the integrality option of
        /// SciPy's <c>differential_evolution</c>. Its bounds must contain an integer. The polish step, when
        /// enabled, is skipped for a problem with integer parameters, and a Hessian is not meaningful for them:
        /// set <see cref="Optimizer.ComputeHessian"/> to false.
        /// </remarks>
        public IList<bool>? IntegerParameters { get; set; }

        /// <summary>
        /// Determines whether each generation is evaluated in parallel. Default = false.
        /// </summary>
        /// <remarks>
        /// <para>
        /// When true, the search uses deferred updating: every trial vector of a generation is built from the
        /// population at the start of the generation, the trials are evaluated concurrently, and selection then
        /// runs in index order (as SciPy's <c>updating='deferred'</c>). The random draws are made in the same
        /// order on one thread, so the result is reproducible and does not depend on the number of threads.
        /// It differs from the default immediate updating, where an accepted trial can serve as a donor later
        /// in the same generation, so the two modes give different, equally valid, searches.
        /// </para>
        /// <para>
        /// The objective function must be thread safe. The speedup approaches the number of cores when an
        /// evaluation costs much more than the bookkeeping, as for a simulation model.
        /// </para>
        /// </remarks>
        public bool UseParallel { get; set; }

        /// <summary>
        /// The maximum number of concurrent evaluations when <see cref="UseParallel"/> is true, or -1 for no
        /// limit. Default = -1. The result does not depend on it.
        /// </summary>
        public int MaxDegreeOfParallelism { get; set; } = -1;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            if (PopulationSize < 1) throw new ArgumentOutOfRangeException(nameof(PopulationSize), "The population size must be greater than 0.");
            if (Mutation < 0 || Mutation > 2) throw new ArgumentOutOfRangeException(nameof(Mutation), "The mutation parameter must be between 0 and 2.");
            if (DitherRate < 0 || DitherRate > 1) throw new ArgumentOutOfRangeException(nameof(DitherRate), "The dithering rate must be between 0 and 1.");
            if (CrossoverProbability < 0 || CrossoverProbability > 1) throw new ArgumentOutOfRangeException(nameof(CrossoverProbability), "The crossover probability must be between 0 and 1.");

            if (Polish) ValidatePolishMethod(PolishMethod);
            int i, j, D = NumberOfParameters;
            bool cancel = false;
            bool[]? isInteger = ValidateIntegerParameters();
            ValidateInitialPopulation();
            if (MaxDegreeOfParallelism == 0 || MaxDegreeOfParallelism < -1) throw new ArgumentOutOfRangeException(nameof(MaxDegreeOfParallelism), "The maximum degree of parallelism must be positive or -1.");
            if (UseParallel)
            {
                OptimizeDeferred(isInteger);
                return;
            }

            // Initialize the population of points
            var prng = new MersenneTwister(PRNGSeed);
            var rnd = LatinHypercube.Random(PopulationSize, D, PRNGSeed);
            var Xp = new List<ParameterSet>();
            for (i = 0; i < PopulationSize; i++)
            {
                var values = new double[D];
                if (InitialPopulation != null && i < InitialPopulation.Count)
                {
                    Array.Copy(InitialPopulation[i], values, D);
                }
                else
                {
                    for (j = 0; j < D; j++)
                        values[j] = LowerBounds[j] + rnd[i, j] * (UpperBounds[j] - LowerBounds[j]);
                }
                if (isInteger != null) RoundIntegers(values, isInteger);
                var fitness = Evaluate(values, ref cancel);
                Xp.Add(new ParameterSet(values, fitness));
                if (cancel == true) return;
            }

            Iterations += 1;

            // Trial vector buffer. A rejected trial leaves it free for the next one; an accepted trial
            // moves into the population and a new buffer is allocated.
            var u = new double[D];

            // Perform Differential Evolution
            while (Iterations < MaxIterations)
            {
                // Keep track of population statistics to assess convergence
                var statistics = new RunningStatistics();

                // Mutate and recombine population
                for (i = 0; i < PopulationSize; i++)
                {

                    // Randomly select three vectors indexes
                    int r0, r1, r2;
                    if (PopulationSize >= 4)
                    {
                        SelectDonors(prng, PopulationSize, i, out r0, out r1, out r2);
                    }
                    else
                    {
                        // With fewer than four members there are not three distinct donors, and indexing fails.
                        var indices = Enumerable.Range(0, PopulationSize).Where(x => x != i).OrderBy(_ => prng.NextDouble()).Take(3).ToArray();
                        r0 = indices[0]; r1 = indices[1]; r2 = indices[2];
                    }

                    // Generate trial vector
                    double G = prng.NextDouble() <= DitherRate ? 0.5 + prng.NextDouble() * 0.5 : Mutation; // Scale factor
                    double jRand = prng.Next(0, D);
                    for (j = 0; j < D; j++)
                    {
                        var rr = prng.NextDouble();
                        if (rr <= CrossoverProbability || j == jRand)
                        {
                            u[j] = Xp[r0].Values[j] + G * (Xp[r1].Values[j] - Xp[r2].Values[j]);
                            // Repair halfway toward the target so infeasible trials neither accumulate exactly
                            // on a bound nor consume an extra random draw that shifts the seeded DE trajectory.
                            if (u[j] < LowerBounds[j])
                                u[j] = 0.5d * Xp[i].Values[j] + 0.5d * LowerBounds[j];
                            else if (u[j] > UpperBounds[j])
                                u[j] = 0.5d * Xp[i].Values[j] + 0.5d * UpperBounds[j];
                        }
                        else
                        {
                            u[j] = Xp[i].Values[j];
                        }
                    }

                    // Evaluate fitness
                    if (isInteger != null) RoundIntegers(u, isInteger);
                    var fitness = Evaluate(u, ref cancel);
                    if (cancel) return;

                    // Update population
                    if (fitness <= Xp[i].Fitness)
                    {
                        Xp[i] = new ParameterSet(u, fitness);
                        u = new double[D];
                    }
              
                    // Keep running stats of population
                    statistics.Push(Xp[i].Fitness);

                }

                // Evaluate convergence 
                if (Iterations >= 10 && statistics.StandardDeviation < AbsoluteTolerance + RelativeTolerance * Math.Abs(statistics.Mean))
                {
                    if (Polish && isInteger == null && !PolishBestParameterSet(PolishMethod)) return;
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }

                Iterations += 1;

            }

            // Reaching this point means the maximum iterations were reached.
            if (Polish && isInteger == null && !PolishBestParameterSet(PolishMethod)) return;
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

        /// <summary>
        /// Differential evolution with deferred updating and parallel evaluation of each generation.
        /// </summary>
        /// <param name="isInteger">The integer flags, or null.</param>
        /// <remarks>
        /// The mutation, crossover and repair rules and the random draws per trial are those of the immediate
        /// search; only the donors' generation and the evaluation order differ. See <see cref="UseParallel"/>.
        /// </remarks>
        private void OptimizeDeferred(bool[]? isInteger)
        {
            int i, j, D = NumberOfParameters, N = PopulationSize;
            bool cancel = false;
            var prng = new MersenneTwister(PRNGSeed);
            var rnd = LatinHypercube.Random(N, D, PRNGSeed);
            var trials = new double[N][];
            var fitness = new double[N];
            for (i = 0; i < N; i++)
            {
                trials[i] = new double[D];
                if (InitialPopulation != null && i < InitialPopulation.Count)
                    Array.Copy(InitialPopulation[i], trials[i], D);
                else
                    for (j = 0; j < D; j++)
                        trials[i][j] = LowerBounds[j] + rnd[i, j] * (UpperBounds[j] - LowerBounds[j]);
                if (isInteger != null) RoundIntegers(trials[i], isInteger);
            }
            int done = EvaluateBatch(trials, N, fitness, MaxDegreeOfParallelism, ref cancel);
            var Xp = new ParameterSet[N];
            for (i = 0; i < done; i++)
            {
                Xp[i] = new ParameterSet(trials[i], fitness[i]);
                trials[i] = new double[D];
            }
            if (cancel) return;

            Iterations += 1;

            while (Iterations < MaxIterations)
            {
                // Build every trial from the population at the start of the generation.
                for (i = 0; i < N; i++)
                {
                    int r0, r1, r2;
                    if (N >= 4)
                    {
                        SelectDonors(prng, N, i, out r0, out r1, out r2);
                    }
                    else
                    {
                        var indices = Enumerable.Range(0, N).Where(x => x != i).OrderBy(_ => prng.NextDouble()).Take(3).ToArray();
                        r0 = indices[0]; r1 = indices[1]; r2 = indices[2];
                    }
                    double G = prng.NextDouble() <= DitherRate ? 0.5 + prng.NextDouble() * 0.5 : Mutation;
                    double jRand = prng.Next(0, D);
                    var u = trials[i];
                    for (j = 0; j < D; j++)
                    {
                        var rr = prng.NextDouble();
                        if (rr <= CrossoverProbability || j == jRand)
                        {
                            u[j] = Xp[r0].Values[j] + G * (Xp[r1].Values[j] - Xp[r2].Values[j]);
                            if (u[j] < LowerBounds[j])
                                u[j] = 0.5d * Xp[i].Values[j] + 0.5d * LowerBounds[j];
                            else if (u[j] > UpperBounds[j])
                                u[j] = 0.5d * Xp[i].Values[j] + 0.5d * UpperBounds[j];
                        }
                        else
                        {
                            u[j] = Xp[i].Values[j];
                        }
                    }
                    if (isInteger != null) RoundIntegers(u, isInteger);
                }

                // Evaluate the generation concurrently, then select in index order.
                done = EvaluateBatch(trials, N, fitness, MaxDegreeOfParallelism, ref cancel);
                var statistics = new RunningStatistics();
                for (i = 0; i < done; i++)
                {
                    if (fitness[i] <= Xp[i].Fitness)
                    {
                        Xp[i] = new ParameterSet(trials[i], fitness[i]);
                        trials[i] = new double[D];
                    }
                    statistics.Push(Xp[i].Fitness);
                }
                if (cancel) return;

                if (Iterations >= 10 && statistics.StandardDeviation < AbsoluteTolerance + RelativeTolerance * Math.Abs(statistics.Mean))
                {
                    if (Polish && isInteger == null && !PolishBestParameterSet(PolishMethod)) return;
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }

                Iterations += 1;
            }

            if (Polish && isInteger == null && !PolishBestParameterSet(PolishMethod)) return;
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

        /// <summary>
        /// Validates <see cref="IntegerParameters"/> and returns the flags, or null when no parameter is an integer.
        /// </summary>
        /// <returns>The integer flags, or null.</returns>
        private bool[]? ValidateIntegerParameters()
        {
            if (IntegerParameters == null) return null;
            if (IntegerParameters.Count != NumberOfParameters)
                throw new ArgumentOutOfRangeException(nameof(IntegerParameters), "There must be one integer flag per parameter.");
            bool any = false;
            for (int j = 0; j < NumberOfParameters; j++)
            {
                if (!IntegerParameters[j]) continue;
                any = true;
                if (Math.Ceiling(LowerBounds[j]) > Math.Floor(UpperBounds[j]))
                    throw new ArgumentOutOfRangeException(nameof(IntegerParameters), "The bounds of an integer parameter must contain an integer.");
            }
            return any ? IntegerParameters.ToArray() : null;
        }

        /// <summary>
        /// Validates <see cref="InitialPopulation"/>.
        /// </summary>
        private void ValidateInitialPopulation()
        {
            if (InitialPopulation == null) return;
            if (InitialPopulation.Count > PopulationSize)
                throw new ArgumentOutOfRangeException(nameof(InitialPopulation), "There cannot be more initial members than the population size.");
            foreach (var row in InitialPopulation)
            {
                if (row == null || row.Length != NumberOfParameters)
                    throw new ArgumentOutOfRangeException(nameof(InitialPopulation), "Each initial member must have one value per parameter.");
                for (int j = 0; j < NumberOfParameters; j++)
                    if (!(row[j] >= LowerBounds[j] && row[j] <= UpperBounds[j]))
                        throw new ArgumentOutOfRangeException(nameof(InitialPopulation), "The initial members must be within the lower and upper bounds.");
            }
        }

        /// <summary>
        /// Rounds the integer parameters to the nearest integer within their bounds.
        /// </summary>
        /// <param name="values">The point, modified in place.</param>
        /// <param name="isInteger">The integer flags.</param>
        private void RoundIntegers(double[] values, bool[] isInteger)
        {
            for (int j = 0; j < values.Length; j++)
            {
                if (!isInteger[j]) continue;
                double r = Math.Round(values[j], MidpointRounding.ToEven);
                double lo = Math.Ceiling(LowerBounds[j]), hi = Math.Floor(UpperBounds[j]);
                values[j] = Tools.Clamp(r, lo, hi);
            }
        }

        /// <summary>
        /// Selects three distinct donor indexes other than the target, in one pass with no allocation.
        /// </summary>
        /// <param name="prng">The generator. One uniform draw is consumed per population member other than the target.</param>
        /// <param name="populationSize">The population size, at least four.</param>
        /// <param name="target">The index of the target vector, which is never selected.</param>
        /// <param name="r0">The donor with the smallest draw.</param>
        /// <param name="r1">The donor with the second smallest draw.</param>
        /// <param name="r2">The donor with the third smallest draw.</param>
        /// <remarks>
        /// This is a stable random-key shuffle, equal to
        /// <c>Enumerable.Range(0, populationSize).Where(x => x != target).OrderBy(_ => prng.NextDouble()).Take(3)</c>:
        /// the keys are drawn in index order and the three smallest are kept, ties going to the lower index.
        /// It costs O(N) with no allocation.
        /// </remarks>
        internal static void SelectDonors(MersenneTwister prng, int populationSize, int target, out int r0, out int r1, out int r2)
        {
            r0 = r1 = r2 = -1;
            double k0 = double.PositiveInfinity, k1 = double.PositiveInfinity, k2 = double.PositiveInfinity;
            for (int x = 0; x < populationSize; x++)
            {
                if (x == target) continue;
                double key = prng.NextDouble();
                // Strict comparisons keep an earlier index ahead of a later one with an equal key.
                if (r0 < 0 || key < k0)
                {
                    r2 = r1; k2 = k1;
                    r1 = r0; k1 = k0;
                    r0 = x; k0 = key;
                }
                else if (r1 < 0 || key < k1)
                {
                    r2 = r1; k2 = k1;
                    r1 = x; k1 = key;
                }
                else if (r2 < 0 || key < k2)
                {
                    r2 = x; k2 = key;
                }
            }
        }

    }
}
