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

using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// The Adaptive Movement (Adam) optimization algorithm. The objective function must be differentiable and convex.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Adaptive Movement Estimation algorithm, or Adam for short, is an extension to the gradient descent optimization algorithm.
    /// Adam is designed to accelerate the optimization process, e.g. decrease the number of function evaluations required to reach the optima, 
    /// or to improve the capability of the optimization algorithm, e.g. result in a better final result.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Kingma and Ba (2014) Adam: A method for Stochastic Optimization
    /// <see href="https://arxiv.org/abs/1412.6980"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class ADAM : Optimizer
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
        public ADAM(Func<double[], double> objectiveFunction, int numberOfParameters,
                    IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds, double alpha = 0.001,
                    Func<double[], double[]> gradient = null) : base(objectiveFunction, numberOfParameters)
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
        /// Gets and sets the decay factor for the first momentum. Default = 0.9.
        /// </summary>
        public double Beta1 { get; set; } = 0.9;

        /// <summary>
        /// Gets and sets the decay factor for infinity norm. Default = 0.999.
        /// </summary>
        public double Beta2 { get; set; } = 0.999;

        /// <summary>
        /// The function for evaluating the gradient of the objective function.
        /// </summary>
        public Func<double[], double[]> Gradient;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            int D = NumberOfParameters;
            bool cancel = false;
            var p = InitialValues.ToArray();
            var g = new double[D];
            var m = new double[D];
            var v = new double[D];
            double mhat, vhat;
            double f0 = Evaluate(p, ref cancel);
            double f1;

            while (Iterations < MaxIterations)
            {
                // Get gradient with respect to objective function
                g = Gradient != null ? Gradient(p) : NumericalDerivative.Gradient((x) => Evaluate(x, ref cancel), p);
                if (cancel) return;

                // Update parameters
                for (int i = 0; i < D; i++)
                {
                    // Update biased first and second moment estimates
                    m[i] = Beta1 * m[i] + (1d - Beta1) * g[i];
                    v[i] = Beta2 * v[i] + (1d - Beta2) * g[i] * g[i];
                    // Compute biased corrected moment estimates
                    mhat = m[i] / (1d - Math.Pow(Beta1, Iterations + 1));
                    vhat = v[i] / (1d - Math.Pow(Beta2, Iterations + 1));
                    // Update parameter
                    p[i] = p[i] - Alpha * mhat / (Math.Sqrt(vhat) + Tools.DoubleMachineEpsilon);
                    // Make sure the parameter is within the bounds.
                    p[i] = RepairParameter(p[i], LowerBounds[i], UpperBounds[i]);
                }

                // Evaluate the objective function at the new point
                f1 = Evaluate(p, ref cancel);
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
