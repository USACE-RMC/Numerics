using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The adaptive simulated annealing algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This algorithm is interpreted as a slow decrease in the probability of accepting worse solutions as the solution
    /// space is explored. Accepting the worse solutions allows for a more extensive search for the global solution.
    /// At each step, the algorithm randomly selects a solution close to the current one, measures its quality, and 
    /// probabilistically decides moving between the two states.
    /// </para>
    /// <para>
    ///     <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Implements routine described by Corana et al. "Minimizing Multimodal Functions of Continuous Variables with 'Simulated Annealing' Algorithm" (1987).
    /// </description></item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Simulated_annealing"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class SimulatedAnnealing : Optimizer
    {

        /// <summary>
        /// Construct a new simulated annealing optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public SimulatedAnnealing(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
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

        /// <inheritdoc />
        protected override double[]? ParameterLowerBounds => LowerBounds;

        /// <inheritdoc />
        protected override double[]? ParameterUpperBounds => UpperBounds;

        /// <summary>
        /// Determines whether the best point is refined by a local search when the global search ends. Default = true.
        /// </summary>
        /// <remarks>
        /// The polish runs after convergence or after <see cref="Optimizer.MaxIterations"/>, never after the
        /// evaluation budget is spent. Its evaluations count toward <see cref="Optimizer.FunctionEvaluations"/>
        /// and the reported status is the global search's.
        /// </remarks>
        public bool Polish { get; set; } = true;

        /// <summary>
        /// The local method used to polish: BFGS, Nelder-Mead, or Powell. Default = Nelder-Mead.
        /// </summary>
        public LocalMethod PolishMethod { get; set; } = LocalMethod.NelderMead;

        /// <summary>
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// The initial temperature at the start of the algorithm. Default = 10. 
        /// </summary>
        public double InitialTemperature { get; set; } = 10;
        
        /// <summary>
        /// The minimum temperature allowable. The temperature will be held constant when it reaches this point. 
        /// </summary>
        public double MinTemperature { get; set; } = 0.1;

        /// <summary>
        /// The cooling rate for the annealing schedule. Default = 0.95.
        /// </summary>
        public double CoolingRate { get; set; } = 0.95;

        /// <summary>
        /// The number of cycles before updating the step size. Default = 10.
        /// </summary>
        public int UpdateCycles { get; set; } = 4;

        /// <summary>
        /// The number of cycles before reducing the temperature. 10. 
        /// </summary>
        public int TemperatureCycles { get; set; } = 10;

        /// <summary>
        /// The number of successive temperature reductions to test for termination. Default = 20.
        /// </summary>
        /// <remarks>
        /// The search stops with <see cref="OptimizationStatus.Success"/> when the best fitness has changed by
        /// no more than the tolerances (<see cref="Optimizer.CheckConvergence"/>) over this many successive
        /// temperature reductions. The best point is then polished when <see cref="Polish"/> is true.
        /// </remarks>
        public int ToleranceSteps { get; set; } = 20;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            if (Polish) ValidatePolishMethod(PolishMethod);
            int i, j, k, D = NumberOfParameters;

            // Validate
            if (InitialTemperature < 1) throw new ArgumentOutOfRangeException(nameof(InitialTemperature), "The initial temperature must be greater than 0.");
            if (UpdateCycles < 4) throw new ArgumentOutOfRangeException(nameof(UpdateCycles), "The number of cycles before updating the step size must be at least 4.");
            if (TemperatureCycles < 4) throw new ArgumentOutOfRangeException(nameof(TemperatureCycles), "The number of cycles before reducing the temperature must be at least 4.");
            if (ToleranceSteps < 4) throw new ArgumentOutOfRangeException(nameof(ToleranceSteps), "The number of tolerance steps must be greater than 3.");
 
            // Setup variables
            double T = InitialTemperature;  
            bool cancel = false;
            var x = new Vector(D);
            var xp = new Vector(D);
            var trial = new double[D];
            var fitness = new List<double>();
            var prng = new MersenneTwister(PRNGSeed);

            // Initial evaluation
            for (j = 0; j < D; j++) x[j] = 0.5 * (LowerBounds[j] + UpperBounds[j]);
            double fx = Evaluate(x.ToArray(), ref cancel);
            double fxp, df; 
            
            // variables for corona update
            var acceptances = new int[D];
            var v = new double[D];
            var c = new double[D];
            for (j = 0; j < D; ++j)
            {
                v[j] = 1.0 / InitialTemperature;
                c[j] = 2.0;
            }

            // Perform adaptive simulated annealing
            while (Iterations < MaxIterations)
            {
                // Reduce temperature
                T *= CoolingRate;
                T = Math.Max(MinTemperature, T);

                for (i = 1; i <= TemperatureCycles; i++)
                {
                    for (j = 1; j <= UpdateCycles; j++)
                    {
                        // Perform a cycle of random moves, each along one coordinate direction. 
                        // Accept or reject each point according to the Metropolis criterion. 
                        // Record the optimum point reached so far.             
                        for (k = 0; k < D; k++)
                        {
                            // Create the proposal x + e_k * v[k] * u in place, element by element.
                            double step = prng.NextDouble() * 2d - 1d;
                            for (int m = 0; m < D; m++)
                                xp[m] = x[m] + (m == k ? 1d : 0d) * v[k] * step;

                            // If XP is out of bounds, continue
                            if (xp[k] < LowerBounds[k] || xp[k] > UpperBounds[k]) continue;
                                
                            // Evaluate proposal vector
                            // The objective gets its own copy, so it cannot alter the proposal.
                            Array.Copy(xp.Array, trial, D);
                            fxp = Evaluate(trial, ref cancel);
                            if (cancel) return;
                            df = fxp - fx;

                            // Metropolis rule
                            if (df < 0 || prng.NextDouble() < Math.Exp(-df / T))
                            {
                                // Accept proposal
                                acceptances[k] += 1;
                                Array.Copy(xp.Array, x.Array, D);
                                fx = fxp;

                            }
                        }
                    }

                    // Do Corana et al. update
                    // Goal is to keep the acceptance rate near 50%
                    for (j = 0; j < D; j++)
                    {
                        // Get acceptance rate
                        double rate = (double)acceptances[j] / UpdateCycles;
                        if (rate > 0.6)
                        {
                            v[j] *= 1d + c[j] * (rate - 0.6) / 0.4;
                        }
                        else if (rate < 0.4)
                        {
                            v[j] /= 1d + c[j] * (0.4 - rate) / 0.4;
                        }
                    }
                    acceptances.Fill(0);
                }

                Iterations += 1;

                // Stop when the best value has not changed, within the tolerances, over the last
                // ToleranceSteps temperature reductions, then polish the best point.
                fitness.Add(BestParameterSet.Fitness);
                if (fitness.Count > ToleranceSteps)
                {
                    bool stalled = true;
                    for (int u = 1; u <= ToleranceSteps && stalled; u++)
                        stalled = CheckConvergence(fitness[fitness.Count - 1 - u], BestParameterSet.Fitness);
                    if (stalled)
                    {
                        if (Polish && !PolishBestParameterSet(PolishMethod)) return;
                        UpdateStatus(OptimizationStatus.Success);
                        return;
                    }
                }

                // Set current point to the optimum.           
                // Copied, because x is updated in place.
                Array.Copy(BestParameterSet.Values, x.Array, D);
                fx = BestParameterSet.Fitness;

            }

            // Reaching this point means the maximum iterations were reached.
            if (Polish && !PolishBestParameterSet(PolishMethod)) return;
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

    }
}
