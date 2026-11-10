using System;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// A snapshot of an optimization run, reported through <see cref="Optimizer.Progress"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class OptimizationProgress
    {
        /// <summary>
        /// Constructs a progress snapshot.
        /// </summary>
        /// <param name="iterations">The iterations completed so far.</param>
        /// <param name="functionEvaluations">The function evaluations made so far.</param>
        /// <param name="bestFitness">The best fitness so far, as minimized by the solver.</param>
        /// <param name="bestValues">A copy of the best parameter values so far.</param>
        public OptimizationProgress(int iterations, int functionEvaluations, double bestFitness, double[] bestValues)
        {
            Iterations = iterations;
            FunctionEvaluations = functionEvaluations;
            BestFitness = bestFitness;
            BestValues = bestValues;
        }

        /// <summary>
        /// The iterations completed so far.
        /// </summary>
        public int Iterations { get; }

        /// <summary>
        /// The function evaluations made so far.
        /// </summary>
        public int FunctionEvaluations { get; }

        /// <summary>
        /// The best fitness so far, as minimized by the solver. After <see cref="Optimizer.Maximize"/> this is
        /// the negated objective value, as in <see cref="Optimizer.BestParameterSet"/>.
        /// </summary>
        public double BestFitness { get; }

        /// <summary>
        /// A copy of the best parameter values so far.
        /// </summary>
        public double[] BestValues { get; }
    }
}
