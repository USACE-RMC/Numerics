using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// The Gradient Descent optimization algorithm. The objective function must be differentiable and convex.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Gradient descent is a method for unconstrained mathematical optimization. It is a first-order 
    /// iterative algorithm for finding a local minimum of a differentiable multivariate function.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Gradient_descent"/>
    /// </description></item>
    /// <item><description>
    /// <see href="https://towardsdatascience.com/gradient-descent-algorithm-a-deep-dive-cf04e8115f21"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class GradientDescent : Optimizer
    {

        /// <summary>
        /// Construct a new Gradient Descent optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values to evaluate.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="alpha">Optional. The step size or learning rate. Default = 0.001.</param>
        /// <param name="gradient">Optional. Function to evaluate the gradient. Default uses finite difference.</param>
        public GradientDescent(Func<double[], double> objectiveFunction, int numberOfParameters,
                    IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds, double alpha = 0.001,
                    Func<double[], double[]>? gradient = null) : base(objectiveFunction, numberOfParameters)
        {
            // Check if the length of the initial, lower and upper bounds equal the number of parameters
            if (initialValues.Count != numberOfParameters || lowerBounds.Count != numberOfParameters || upperBounds.Count != numberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(lowerBounds), "The initial values and lower and upper bounds must be the same length as the number of parameters.");
            }
            // Check if the initial values are between the lower and upper values
            for (int j = 0; j < initialValues.Count; j++)
            {
                if (upperBounds[j] < lowerBounds[j])
                {
                    throw new ArgumentOutOfRangeException(nameof(upperBounds), "The upper bound cannot be less than the lower bound.");
                }
                if (initialValues[j] < lowerBounds[j] || initialValues[j] > upperBounds[j])
                {
                    throw new ArgumentOutOfRangeException(nameof(initialValues), "The initial values must be between the upper and lower bounds.");
                }
            }
            InitialValues = initialValues.ToArray();
            LowerBounds = lowerBounds.ToArray();
            UpperBounds = upperBounds.ToArray();
            Alpha = alpha;
            Gradient = gradient;
        }

        /// <summary>
        /// An array of initial values to evaluate. 
        /// </summary>
        public double[] InitialValues { get; private set; }

        /// <summary>
        /// An array of lower bounds (inclusive) of the interval containing the optimal point. 
        /// </summary>
        public double[] LowerBounds { get; private set; }

        /// <summary>
        /// An array of upper bounds (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double[] UpperBounds { get; private set; }

        /// <summary>
        /// Gets and sets the step size, or learning rate. Default = 0.001.
        /// </summary>
        public double Alpha { get; set; } = 0.001;

        /// <summary>
        /// The function for evaluating the gradient of the objective function.
        /// </summary>
        public Func<double[], double[]>? Gradient;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            int D = NumberOfParameters;
            bool cancel = false;
            var p0 = InitialValues.ToArray();
            var g = new double[D];
            double f0 = Evaluate(p0, ref cancel);
            double f1;

            while (Iterations < MaxIterations)
            {
                // Get gradient with respect to objective function
                g = Gradient != null ? Gradient(p0) : NumericalDerivative.Gradient((x) => Evaluate(x, ref cancel), p0);
                if (cancel) return;

                // Update parameters
                for (int i = 0; i < D; i++)
                {
                    // Update parameter
                    p0[i] = p0[i] - Alpha * g[i];
                    // Make sure the parameter is within the bounds.
                    p0[i] = RepairParameter(p0[i], LowerBounds[i], UpperBounds[i]);
                }

                // Evaluate the objective function at the new point
                f1 = Evaluate(p0, ref cancel);
                if (cancel) return;

                // Check convergence.
                if (CheckConvergence(f0, f1))
                {
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }

                f0 = f1;
                Iterations += 1;
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }
    }
}
