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
    /// The Multi-Level Single Linkage (MLSL) optimization method. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// <item><description>
    /// </description></item>
    ///     Brian Skahill, USACE Engineer Research and Development Center Coastal and Hydraulics Laboratory, brian.e.skahill@usace.army.mil
    /// </description></item>
    /// </list>
    /// </para>
    /// <b> Description: </b>
    /// 
    ///     The key optimization steps are as follows:
    /// <list type="number">
    /// <item><description>
    ///     Generate sample points and function values. Add N points, drawn from a uniform distribution over S, to the (initially empty) set of sample points, 
    ///     and evaluate f(x) at each new sample point.
    /// </description></item>
    /// <item><description>
    ///     Reduce the sample points. Sort the entire sample of kN points in order of increasing object function values.
    ///     Select the γkN points with the lowest objective function values. This resultant set, Rk, is called the reduced sample.
    /// </description></item>
    /// <item><description>
    ///     Select start points for local searches. Determine a (possibly empty) subset of the sample points from which to start local searches.
    /// </description></item>
    /// <item><description>
    ///     Decide whether to stop. If stopping rule is satisfied, regard the lowest local minimizer as the global minimizer, otherwise go to step 1.
    /// </description></item>
    /// </list>
    /// <para>
    ///     <b> References: </b>
    /// </para>
    /// <para> 
    ///     Implements routine described by 
    ///     Kan A.H.G.R., Boender C.G.E., Timmer G.T. (1985) A Stochastic Approach to Global Optimization. 
    ///     In: Schittkowski K. (eds) Computational Mathematical Programming. NATO ASI Series (Series F: Computer and Systems Sciences), 
    ///     vol 15. Springer, Berlin, Heidelberg. <see href="https://doi.org/10.1007/978-3-642-82450-0_10"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public class MLSL : Optimizer
    {
        /// <summary>
        /// Constructs a new multi-level single linkage (MLSL) optimization method.
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values to evaluate.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="method">The local search method to use. Default = BFGS.</param>
        public MLSL(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds, LocalMethod method = LocalMethod.BFGS) : base(objectiveFunction, numberOfParameters)
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

        /// <summary>
        /// The number of random samples to evaluate each iteration. Default = 30
        /// </summary>
        public int SampleSize { get; set; } = 30;

        /// <summary>
        /// Determines the reduced sample size. Must be between 0 and 1. Default = 0.05.
        /// </summary>
        public double Gamma { get; set; } = 0.05;

        /// <summary>
        /// The scale parameter for determining the critical distance. Default = 2.0. 
        /// </summary>
        public double Sigma { get; set; } = 2.0;

        /// <summary>
        /// The list of all sampled points.
        /// </summary>
        public List<SamplePoint> SampledPoints { get; private set; }

        /// <summary>
        /// The list of all local optimums.
        /// </summary>
        public List<SamplePoint> LocalMinimums { get; private set; }

        /// <summary>
        /// The minimum number of iterations to carry out with no improvement. Default = 5.
        /// </summary>
        public int MinNoImprovement { get; set; } = 5;

        /// <summary>
        /// The maximum number of iterations to carry out with no improvement. Default = 10.
        /// </summary>
        public int MaxNoImprovement { get; set; } = 10;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            if (SampleSize < 4) throw new ArgumentOutOfRangeException(nameof(SampleSize), "The sample size must be greater than or equal to 4.");
            if (Gamma <= 0 || Gamma >= 1) throw new ArgumentOutOfRangeException(nameof(Gamma), "The reduction parameter must be between 0 and 1.");
            if (Sigma <= 0) throw new ArgumentOutOfRangeException(nameof(Sigma), "The scale parameter must be positive.");
            if (MinNoImprovement < 1) throw new ArgumentOutOfRangeException(nameof(MinNoImprovement), "The minimum number of iterations to carry out with no improvement must be greater than 0.");
            if (MaxNoImprovement < 1) throw new ArgumentOutOfRangeException(nameof(MaxNoImprovement), "The maximum number of iterations to carry out with no improvement must be greater than 0.");
            if (MaxNoImprovement < MinNoImprovement) throw new ArgumentOutOfRangeException(nameof(MaxNoImprovement), "The maximum number of iterations must be greater than or equal to the minimum number.");

            int D = NumberOfParameters;
            int N = SampleSize;
            double oldFit = double.MaxValue;
            int noImprovement = 0;
            bool cancel = false;
            Optimizer solver = null;
            var prng = new MersenneTwister(PRNGSeed);

            // Set lower and upper bounds and
            // create uniform distributions for each parameter
            var uniformDists = new List<Uniform>();
            for (int i = 0; i < D; i++)
                uniformDists.Add(new Uniform(LowerBounds[i], UpperBounds[i]));

            // Critical distance
            double rkfactor = Math.Sqrt(Math.PI) * Math.Pow(SpecialFunctions.Gamma.Function(D) * Sigma, 1.0 / D);
            for (int i = 0; i < D; i++)
                rkfactor *= Math.Pow(UpperBounds[i] - LowerBounds[i], 1.0 / D);
     
            // Create lists for sample points
            SampledPoints = new List<SamplePoint>();
            LocalMinimums = new List<SamplePoint>();

            while (Iterations < MaxIterations)
            {
                // Step 1. Generate sample points and function values
                // Add N points, drawn from a uniform distribution over S, to the (initially empty) set of sample points,
                // and evaluate f(x) at each new sample point.

                if (Iterations == 0)
                {
                    // On the first iteration, add the user-defined initial starting points
                    // This can often be very close to the true minimum
                    SampledPoints.Add(new SamplePoint { ParameterSet = new ParameterSet(InitialValues, Evaluate(InitialValues, ref cancel)), Minimized = true });

                    // Perform local minimizations from initial values
                    solver = GetLocalOptimizer(InitialValues, LocalRelativeTolerance, LocalAbsoluteTolerance, ref cancel);
                    solver.Minimize();
                    if (cancel) return;

                    LocalMinimums.Add(new SamplePoint { ParameterSet = solver.BestParameterSet.Clone(), Minimized = true });
                }

                // Sample a point from the uniform distribution over S.  
                for (int i = 0; i < (Iterations == 0 ? N - 1 : N); i++)
                {
                    var values = new double[(int)D];
                    for (int j = 0; j < D; j++)
                        values[j] = uniformDists[j].InverseCDF(prng.NextDouble());
                   
                    var fitness = Evaluate(values, ref cancel);
                    if (cancel) return;
                    SampledPoints.Add(new SamplePoint { ParameterSet = new ParameterSet(values, fitness) });

                }

                // Step 2. Reduce the sample points 
                // Sort the entire sample of kN points in order of increasing object function values.
                // Select the γkN points with the lowest objective function values. 
                // This resultant set, Rk, is called the reduced sample.

                SampledPoints.Sort((x, y) => x.ParameterSet.Fitness.CompareTo(y.ParameterSet.Fitness));
                int gkN = (int)Math.Ceiling(Gamma * (Iterations + 1) * N);
                var Rk = SampledPoints.Take(gkN).ToList();

                // Step 3. Select start points for local searches
                // Determine a (possibly empty) subset of the sample points from which to start local searches.
                //
                // At iteration k, each sample point x is selected as a start point for a local minimization if it has:
                // -> not been used as a start point at a previous iteration
                // -> there is no sample point y within the critical distance rk of x with a lower function value:
                //    ||x - y|| <= rk and also f(x) < f(y)

                double rk = rkfactor * Math.Pow(Math.Log((Iterations + 1) * N) / ((Iterations + 1) * N), 1.0 / D);

                foreach (var point in Rk)
                {

                    // Check if local search has already been performed
                    if (point.Minimized) continue;

                    bool canMinimize = true;

                    // Check if there is any other point y in Rk such that ||xi - y|| <= rk and also f(x) < f(y)
                    foreach (var other in Rk)
                    {
                        if (point == other) continue;
                        if (Tools.NormalizedDistance(point.ParameterSet.Values, other.ParameterSet.Values, LowerBounds, UpperBounds) <= rk &&
                            other.ParameterSet.Fitness < point.ParameterSet.Fitness)
                        {
                            canMinimize = false;
                            break;
                        }
                    }

                    // Check if there is any other local minimum y such that ||xi - y|| <= rk and also f(x) < f(y)
                    if (canMinimize)
                    {
                        foreach (var localMin in LocalMinimums)
                        {
                            if (Tools.NormalizedDistance(point.ParameterSet.Values, localMin.ParameterSet.Values, LowerBounds, UpperBounds) <= rk &&
                                localMin.ParameterSet.Fitness < point.ParameterSet.Fitness)
                            {
                                canMinimize = false;
                                break;
                            }
                        }
                    }

                    // Check if point is a potential minimizer
                    if (canMinimize)
                    {
                        // Perform local minimizations from all start points
                        solver = GetLocalOptimizer(point.ParameterSet.Values, LocalRelativeTolerance, LocalAbsoluteTolerance, ref cancel);
                        solver.Minimize();
                        if (cancel) return;

                        point.Minimized = true;
                        LocalMinimums.Add(new SamplePoint { ParameterSet = solver.BestParameterSet.Clone(), Minimized = true });
                    }
                    
                }

                // Step 4. Decide whether to stop
                // If stopping rule is satisfied, regard the lowest local minimizer as the global minimizer, otherwise go to step 1.
                double w = LocalMinimums.Count;
                double s = SampledPoints.Count;
                double B1 = (w * (s - 1d)) / (s - w - 2d);
                double B2 = (s - w - 1d) * (s + w) / (s * (s - 1d));

                if (Iterations >= 1 && w >= 1 && oldFit == BestParameterSet.Fitness)
                {
                    noImprovement++;
                }
                else
                {
                    noImprovement = 0;
                    oldFit = BestParameterSet.Fitness;
                }

                // And the Bayesian stopping rule is satisfied
                if (noImprovement >= MaxNoImprovement || (noImprovement >= MinNoImprovement && w >= 1 && B1 - w < 0.5 && B2 >= 0.995))
                {
                    // Polish the final result
                    if (Polish)
                    {
                        solver = GetLocalOptimizer(BestParameterSet.Values, RelativeTolerance, AbsoluteTolerance, ref cancel);
                        solver.Minimize();
                        UpdateStatus(solver.Status);
                        return;
                    }

                    // Stop and return
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }


                Iterations += 1;
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

        /// <summary>
        /// Class for storing sampled points
        /// </summary>
        public class SamplePoint
        {
            /// <summary>
            /// The sample point parameter set.
            /// </summary>
            public ParameterSet ParameterSet { get; set; }

            /// <summary>
            /// Determines if the sample point has already been minimized
            /// </summary>
            public bool Minimized { get; set; } = false;

            /// <summary>
            /// Create new sample point
            /// </summary>
            public SamplePoint()
            {
                ParameterSet = new ParameterSet();
            }

        }

    }
}
