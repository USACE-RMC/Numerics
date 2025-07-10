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
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The Multi-Start (MS) optimization method. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Brian Skahill, USACE Engineer Research and Development Center Coastal and Hydraulics Laboratory, brian.e.skahill@usace.army.mil
    /// </para>
    /// 
    /// <para>
    /// <b> Description: </b>
    ///     With the Multi-Start (MS) approach, a local search procedure P is applied to each point in the random sample; 
    ///     the best local minimum found in this way is our candidate for the global minimum
    /// <list type="number">
    /// <item><description>
    ///     Draw a point from the uniform distribution over S.
    ///     </description></item>
    /// <item><description>
    ///     Apply P to the new sample point.
    ///     </description></item>
    /// <item><description>
    ///     The local minimum x* identified with the lowest function value is the
    ///     candidate value for x^,. Return to Step 1, unless a stopping criterion
    ///     is satisfied. 
    /// </description></item>
    /// </list>
    /// </para>
    /// <para>
    ///     While optimal Bayesian stopping rules can be specified for Multi-Start, the only stopping criterion implemented for this 
    ///     simple ("folklore") global optimization method is the maximum number of local searches to perform.
    /// </para>
    /// 
    /// <para>
    ///     <b> References: </b>
    /// </para>
    /// <para>
    ///     Implements routine described by 
    ///     Kan A.H.G.R., Boender C.G.E., Timmer G.T. (1985) A Stochastic Approach to Global Optimization. 
    ///     In: Schittkowski K. (eds) Computational Mathematical Programming. NATO ASI Series (Series F: Computer and Systems Sciences), 
    ///     vol 15. Springer, Berlin, Heidelberg. <see href="https://doi.org/10.1007/978-3-642-82450-0_10."/>
    /// </para>
    /// </remarks>
    [Serializable]
    public class MultiStart : Optimizer
    {

        /// <summary>
        /// Construct a new multi-start optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values to evaluate.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="method">The local optimization method to use. Default = BFGS.</param>
        public MultiStart(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds, LocalMethod method = LocalMethod.BFGS) : base(objectiveFunction, numberOfParameters)
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
            Method = method;
            MaxIterations = 100;
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
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// The local search method to use. Default = BFGS.
        /// </summary>
        public LocalMethod Method { get; set; } = LocalMethod.BFGS;

        /// <summary>
        /// The desired absolute tolerance for the local solution. Default = ~Sqrt(1E-16), or 1E-8.
        /// </summary>
        public double LocalAbsoluteTolerance { get; set; } = 1E-8;

        /// <summary>
        /// The desired relative tolerance for the local solution. Default = ~Sqrt(1E-16), or 1E-8.
        /// </summary>
        public double LocalRelativeTolerance { get; set; } = 1E-8;

        /// <summary>
        /// If true (default), then a final local search is used to polish the best population member at the end, which can improve the optimization slightly. 
        /// </summary>
        public bool Polish { get; set; } = true;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            int i, j, D = NumberOfParameters;
            bool cancel = false;
            Optimizer solver = null;

            // Set lower and upper bounds and
            // create uniform distributions for each parameter
            var uniformDists = new List<Uniform>();
            for (i = 0; i < D; i++)
                uniformDists.Add(new Uniform(LowerBounds[i], UpperBounds[i]));
            
            // Set solver parameters
            var prng = new MersenneTwister(PRNGSeed);
            var values = new double[D];

            while (Iterations < MaxIterations)
            {
                if (Iterations == 0)
                {
                    values = InitialValues;
                }
                else
                {
                    // Step 1. Draw a point from the uniform distribution over S.              
                    for (j = 0; j < D; j++)
                        values[j] = uniformDists[j].InverseCDF(prng.NextDouble());                
                }

                // Step 2. Apply P to the new sample point.
                solver = GetLocalOptimizer(values, LocalRelativeTolerance, LocalAbsoluteTolerance, ref cancel);
                solver.Minimize();
                if (cancel) return;

                Iterations += 1;
            }

            // Polish the final result
            if (Polish)
            {
                solver = GetLocalOptimizer(BestParameterSet.Values, RelativeTolerance, AbsoluteTolerance, ref cancel);
                solver.Minimize();
                if (cancel) return;
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);

        }

        /// <summary>
        /// Returns an optimizer for the local search. 
        /// </summary>
        /// <param name="initialValues"> An array of initial values to evaluate. </param>
        /// <param name="relativeTolerance">The desired relative tolerance for the solution.</param>
        /// <param name="absoluteTolerance">The desired absolute tolerance for the solution.</param>
        /// <param name="cancel">By ref. Determines if the solver should be canceled.</param>
        private Optimizer GetLocalOptimizer(IList<double> initialValues, double relativeTolerance, double absoluteTolerance, ref bool cancel)
        {
            bool localCancel = false;
            Optimizer solver = null;

            // Make sure the parameters are within the bounds.
            for (int i = 0; i < NumberOfParameters; i++)
                initialValues[i] = RepairParameter(initialValues[i], LowerBounds[i], UpperBounds[i]);

            if (Method == LocalMethod.BFGS)
            {
                solver = new BFGS((x) => Evaluate(x, ref localCancel), NumberOfParameters, initialValues, LowerBounds, UpperBounds) { RelativeTolerance = relativeTolerance, AbsoluteTolerance = absoluteTolerance, MaxFunctionEvaluations = MaxFunctionEvaluations - FunctionEvaluations };
            }
            else if (Method == LocalMethod.NelderMead)
            {
                solver = new NelderMead((x) => Evaluate(x, ref localCancel), NumberOfParameters, initialValues, LowerBounds, UpperBounds) { RelativeTolerance = relativeTolerance, AbsoluteTolerance = absoluteTolerance, MaxFunctionEvaluations = MaxFunctionEvaluations - FunctionEvaluations };
            }
            else if (Method == LocalMethod.Powell)
            {
                solver = new Powell((x) => Evaluate(x, ref localCancel), NumberOfParameters, initialValues, LowerBounds, UpperBounds) { RelativeTolerance = relativeTolerance, AbsoluteTolerance = absoluteTolerance, MaxFunctionEvaluations = MaxFunctionEvaluations - FunctionEvaluations };
            }
            else
            {
                throw new NotSupportedException("Unsupported local method");
            }
            cancel = cancel || localCancel;
            return solver;
        }


    }
}
