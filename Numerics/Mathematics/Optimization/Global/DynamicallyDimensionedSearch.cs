using System;
using System.Collections.Generic;
using System.Linq;
using Numerics.Distributions;
using Numerics.Sampling;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The Dynamically Dimensioned Search (DDS) algorithm, which finds a good global solution within a fixed
    /// budget of objective function evaluations when no gradient is available.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// DDS is a greedy, single-solution search designed for calibrating models whose evaluations are expensive.
    /// The search scales itself to the evaluation budget m: it is global at the start and becomes local as the
    /// budget is spent, by perturbing fewer decision variables of the best solution at each iteration. The steps,
    /// as given in Figure 1 of Tolson and Shoemaker (2007), are as follows.
    /// <list type="number">
    /// <item><description>
    /// Define the perturbation size r (default 0.2), the budget m, the bounds, and an initial solution. Without
    /// an initial solution, the search starts from the best of M = max(5, floor(0.005 m)) uniform random
    /// solutions (paragraph 18), and these evaluations count toward m.
    /// </description></item>
    /// <item><description>
    /// Set the counter i to the number of evaluations made so far and take the best solution as current.
    /// </description></item>
    /// <item><description>
    /// Include each of the D decision variables in the neighborhood with probability P(i) = 1 - ln(i) / ln(m).
    /// If none is included, include one at random.
    /// </description></item>
    /// <item><description>
    /// Perturb each selected variable as x_new = x_best + σ N(0, 1), where σ = r (x_max - x_min). A value below
    /// x_min is reflected to x_min + (x_min - x_new), and if that exceeds x_max it is set to x_min. A value above
    /// x_max is reflected to x_max - (x_new - x_max), and if that is below x_min it is set to x_max.
    /// </description></item>
    /// <item><description>
    /// Evaluate the candidate and accept it as the new best when F(x_new) &lt;= F_best.
    /// </description></item>
    /// <item><description>
    /// Set i = i + 1 and stop when i = m; otherwise return to Step 3.
    /// </description></item>
    /// </list>
    /// </para>
    /// <para>
    /// The paper defines the algorithm in its Figure 1 and paragraphs 16 to 18 rather than in numbered
    /// equations. Boundaries are fully reflective as published; later variants that absorb at a bound with 50%
    /// probability (as in spotpy and OSTRICH) are not used.
    /// </para>
    /// <para>
    /// Spending the whole budget <see cref="EvaluationBudget"/> is the natural end of the search and is reported
    /// as <see cref="OptimizationStatus.Success"/>. <see cref="Optimizer.MaxFunctionEvaluations"/> and
    /// <see cref="Optimizer.MaxIterations"/> remain hard caps: if either is reached first, the run stops early
    /// with the corresponding status, as with the other solvers. <see cref="Optimizer.Iterations"/> counts the
    /// DDS iterations (Steps 3 to 6), which is the budget less the initial evaluations.
    /// </para>
    /// <para>
    ///    <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Tolson, B. A., and C. A. Shoemaker (2007), Dynamically dimensioned search algorithm for computationally
    /// efficient watershed model calibration, Water Resources Research, 43, W01413, doi:10.1029/2005WR004723.
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class DynamicallyDimensionedSearch : Optimizer
    {

        /// <summary>
        /// Construct a new dynamically dimensioned search (DDS) optimization method that starts from the best of a
        /// small uniform random sample.
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public DynamicallyDimensionedSearch(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
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
        /// Construct a new dynamically dimensioned search (DDS) optimization method that starts from a given point.
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values to evaluate, within the bounds.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public DynamicallyDimensionedSearch(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds) : this(objectiveFunction, numberOfParameters, lowerBounds, upperBounds)
        {
            if (initialValues.Count != numberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(initialValues), "The initial values must be the same length as the number of parameters.");
            }
            for (int i = 0; i < initialValues.Count; i++)
            {
                if (double.IsNaN(initialValues[i]) || initialValues[i] < LowerBounds[i] || initialValues[i] > UpperBounds[i])
                {
                    throw new ArgumentOutOfRangeException(nameof(initialValues), "The initial values must be between the upper and lower bounds.");
                }
            }
            InitialValues = initialValues.ToArray();
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
        /// The initial solution, or null when the search starts from the best of a uniform random sample.
        /// </summary>
        public double[]? InitialValues { get; private set; }

        /// <inheritdoc />
        protected override double[]? ParameterLowerBounds => LowerBounds;

        /// <inheritdoc />
        protected override double[]? ParameterUpperBounds => UpperBounds;

        /// <summary>
        /// Determines whether the best point is refined by a local search when the global search ends. Default = false.
        /// </summary>
        /// <remarks>
        /// The polish runs after the planned budget is spent or after <see cref="Optimizer.MaxIterations"/>,
        /// never after <see cref="Optimizer.MaxFunctionEvaluations"/> is reached. Its evaluations are counted in
        /// <see cref="Optimizer.FunctionEvaluations"/>, beyond <see cref="EvaluationBudget"/>, and the reported
        /// status is the global search's.
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
        /// The planned number of objective function evaluations m, including the initial evaluations. Must be at
        /// least 10. Default = 1,000.
        /// </summary>
        /// <remarks>
        /// The search is scaled to this budget through P(i), so it is an input to the algorithm rather than a
        /// safety cap. Spending it ends the run with <see cref="OptimizationStatus.Success"/>. Keep
        /// <see cref="Optimizer.MaxFunctionEvaluations"/> above it, since reaching that cap stops the run early.
        /// </remarks>
        public int EvaluationBudget { get; set; } = 1000;

        /// <summary>
        /// The neighborhood perturbation size r, the standard deviation of a perturbation as a fraction of the
        /// decision variable range, in the range (0, 1]. Default = 0.2 (Tolson and Shoemaker, 2007).
        /// </summary>
        public double PerturbationSize { get; set; } = 0.2;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            if (Polish) ValidatePolishMethod(PolishMethod);
            if (EvaluationBudget < 10) throw new ArgumentOutOfRangeException(nameof(EvaluationBudget), "The evaluation budget must be at least 10.");
            if (double.IsNaN(PerturbationSize) || PerturbationSize <= 0 || PerturbationSize > 1) throw new ArgumentOutOfRangeException(nameof(PerturbationSize), "The perturbation size must be greater than 0 and less than or equal to 1.");

            int i, j, D = NumberOfParameters, m = EvaluationBudget;
            bool cancel = false;
            var prng = new MersenneTwister(PRNGSeed);
            var xBest = new double[D];
            var xNew = new double[D];
            double fBest = double.PositiveInfinity, fitness;
            int count = 0;

            // Steps 1 and 2. Start from the initial solution, or the best of M uniform random solutions.
            if (InitialValues != null)
            {
                Array.Copy(InitialValues, xBest, D);
                fBest = Evaluate(xBest, ref cancel);
                count = 1;
                if (cancel) return;
            }
            else
            {
                int M = Math.Max(5, (int)Math.Floor(0.005 * m));
                for (i = 0; i < M; i++)
                {
                    for (j = 0; j < D; j++)
                        xNew[j] = LowerBounds[j] + (UpperBounds[j] - LowerBounds[j]) * prng.NextDouble();
                    fitness = Evaluate(xNew, ref cancel);
                    count += 1;
                    if (cancel) return;
                    if (i == 0 || fitness <= fBest)
                    {
                        fBest = fitness;
                        Array.Copy(xNew, xBest, D);
                    }
                }
            }

            // Each variable's perturbation standard deviation, and the neighborhood selection flags.
            var sigma = new double[D];
            for (j = 0; j < D; j++)
                sigma[j] = PerturbationSize * (UpperBounds[j] - LowerBounds[j]);
            var selected = new bool[D];
            double logM = Math.Log(m);

            while (count < m)
            {
                if (Iterations >= MaxIterations)
                {
                    if (Polish && !PolishBestParameterSet(PolishMethod)) return;
                    UpdateStatus(OptimizationStatus.MaximumIterationsReached);
                    return;
                }

                // Step 3. Select the neighborhood with probability P(i) = 1 - ln(i) / ln(m).
                double P = 1d - Math.Log(count) / logM;
                int J = 0;
                for (j = 0; j < D; j++)
                {
                    selected[j] = prng.NextDouble() < P;
                    if (selected[j]) J += 1;
                }
                if (J == 0) selected[prng.Next(0, D)] = true;

                // Step 4. Perturb the selected variables of the best solution, reflecting at the bounds.
                for (j = 0; j < D; j++)
                {
                    if (!selected[j])
                    {
                        xNew[j] = xBest[j];
                        continue;
                    }
                    xNew[j] = Reflect(xBest[j] + sigma[j] * Normal.StandardZ(prng.NextDouble()), LowerBounds[j], UpperBounds[j]);
                }

                // Step 5. Greedy acceptance.
                fitness = Evaluate(xNew, ref cancel);
                count += 1;
                Iterations += 1;
                if (cancel) return;
                if (fitness <= fBest)
                {
                    fBest = fitness;
                    Array.Copy(xNew, xBest, D);
                }
            }

            // Step 6. The planned budget is spent, which is the natural end of the search.
            if (Polish && !PolishBestParameterSet(PolishMethod)) return;
            UpdateStatus(OptimizationStatus.Success);
        }

        /// <summary>
        /// Applies the reflecting boundary rule of Step 4 in Figure 1 of Tolson and Shoemaker (2007).
        /// </summary>
        /// <param name="value">The perturbed value.</param>
        /// <param name="lowerBound">The lower bound.</param>
        /// <param name="upperBound">The upper bound.</param>
        /// <returns>The value reflected into the bounds, or the violated bound when one reflection overshoots the other bound.</returns>
        internal static double Reflect(double value, double lowerBound, double upperBound)
        {
            if (value < lowerBound)
            {
                value = lowerBound + (lowerBound - value);
                if (value > upperBound) value = lowerBound;
            }
            else if (value > upperBound)
            {
                value = upperBound - (value - upperBound);
                if (value < lowerBound) value = upperBound;
            }
            return value;
        }

    }
}
