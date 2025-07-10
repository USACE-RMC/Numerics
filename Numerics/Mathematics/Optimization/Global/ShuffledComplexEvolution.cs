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
using Numerics.Distributions;
using Numerics.Sampling;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The Shuffled Complex Evolution (SCE-UA) algorithm, which finds a global minima when no gradient is available.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The main steps of the algorithm are as follows
    /// <list type="number">
    /// <item><description>
    /// Randomly initialize s sample points in the feasible region, and calculate their objective function values.
    /// </description></item>
    /// <item><description>
    /// Rank the sample points in ascending order by objective function values.
    /// </description></item>
    /// <item><description>
    /// Divide the s sorted points into p complexes, and each complex contains m points.
    /// </description></item>
    /// <item><description>
    /// Evolve each complex β times with the competitive complex evolution (CCE) algorithm. In each evolution, 
    /// q points in a complex are selected to form a subcomplex, and the subcomplex is adjusted α times using reflection, 
    /// contraction, and mutation steps.
    /// </description></item>
    /// <item><description>
    /// Mix all points in the evolved complexes and sort them in ascending order of objective function value.
    /// </description></item>
    /// <item><description>
    /// Determine whether the termination criteria are met, if not, return to Step 3.
    /// </description></item>
    /// </list>
    /// </para>
    /// <para>
    ///    <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Implements Duan et al. (1992) SCE-UA probabilistic search algorithm.
    /// </description></item>
    /// <item><description>
    /// <see href="https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2022.1037173/full"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class ShuffledComplexEvolution : Optimizer
    {

        /// <summary>
        /// Construct a new shuffled complex evolution (SCE-UA) optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public ShuffledComplexEvolution(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
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
            CCEIterations = 2 * NumberOfParameters + 1;
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
        /// The pseudo random number generator (PRNG) seed.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// The number of complexes. Default = 5.
        /// </summary>
        public int Complexes { get; set; } = 5;

        /// <summary>
        /// The number of iterations in the inner loop (CCE algorithm).  Defaults = 2 * NumberOfParameters + 1, as recommended by Duan et al (1994). 
        /// </summary>
        public int CCEIterations { get; set; }

        /// <summary>
        /// The number of iterations where the improvement is within the relative tolerance required to confirm convergence. Default = 20.
        /// </summary>
        public int ToleranceSteps { get; set; } = 20;

        /// <inheritdoc/>
        protected override void Optimize()
        {
            if (Complexes < 1) throw new ArgumentOutOfRangeException(nameof(Complexes), "The number of complexes must be greater than 0.");
            if (CCEIterations < 1) throw new ArgumentOutOfRangeException(nameof(CCEIterations), "The number of CCE iterations must be greater than 0.");
            if (ToleranceSteps < 1) throw new ArgumentOutOfRangeException(nameof(ToleranceSteps), "The number of tolerance steps must be greater than 0.");

            int i, j, D = NumberOfParameters;
            bool cancel = false;
            int populationSize = Complexes * CCEIterations; 
            int converge = 0, alpha = 1;
            double oldFit, newFit;

            // Initialize the population of points
            var prng = new MersenneTwister(PRNGSeed);
            var rnd = LatinHypercube.Random(populationSize, D, PRNGSeed);
            var Dpoints = new List<PointFitness>();
            for (i = 0; i < populationSize; i++)
            {
                var values = new double[D];
                for (j = 0; j < D; j++)
                    values[j] = LowerBounds[j] + (UpperBounds[j] - LowerBounds[j]) * rnd[i, j];

                Dpoints.Add(new PointFitness { ParameterSet = new ParameterSet(values, Evaluate(values, ref cancel)), Index = i });
                if (cancel) return;
            }
            oldFit = Dpoints[0].ParameterSet.Fitness;
            Iterations += 1;

            // Initialize shuffled complex evolution (SCE) algorithm
            // Create list of complex prngs
            var complexPRNGs = new List<Random>();
            for (i = 0; i < Complexes; i++)
                complexPRNGs.Add(new Random(prng.Next()));

            // Create trapezoidal cumulative probability for points in complex
            var cdf = Trapezoidal(CCEIterations);

            // Perform the Shuffled complex evolution loop
            while (Iterations < MaxIterations)
            {
                // Partition populations into p complexes distributing points evenly between complexes
                var Acomplex = new List<List<PointFitness>>(); // The A complex
                for (i = 0; i < Complexes; i++)
                {
                    Acomplex.Add(new List<PointFitness>());
                    for (j = 0; j < CCEIterations; j++)
                        Acomplex[i].Add(Dpoints[i + Complexes * j].Clone());
                }

                // Evolve each complex according to competitive complex evolution (CCE) algorithm
                for (i = 0; i < Complexes; i++)
                {
                    var argAcomplex = Acomplex[i];
                    EvolveComplex(alpha, ref argAcomplex, cdf, complexPRNGs[i], ref cancel);
                    Acomplex[i] = argAcomplex;
                    if (cancel) return;
                }

                // Replace A into D
                for (i = 0; i < Complexes; i++)
                {
                    for (j = 0; j < CCEIterations; j++)
                        Dpoints[i + Complexes * j] = Acomplex[i][j].Clone();
                }

                // Rank the points in terms of fitness (i.e., order of increasing function value).
                Dpoints.Sort((x, y) => x.ParameterSet.Fitness.CompareTo(y.ParameterSet.Fitness));

                // Next reset the indexes, so that i = 0 represents the point With the smallest function value.
                for (i = 0; i <= Dpoints.Count - 1; i++)
                    Dpoints[i].Index = i;

                newFit = Dpoints.First().ParameterSet.Fitness;

                // Check convergence
                if (CheckConvergence(oldFit, newFit))
                {
                    // If no improvement increase alpha to increase chance of sub-complexes evolving
                    converge += 1;
                    alpha = Math.Min(3, converge + 1);
                    oldFit = newFit;

                    if (converge >= ToleranceSteps)
                    {
                        UpdateStatus(OptimizationStatus.Success);
                        return;
                    }

                }
                else
                {
                    oldFit = newFit;
                    converge = 0;
                    alpha = 1;
                }

                Iterations += 1;
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }


        /// <summary>
        /// Evolve complex according to the competitive complex evolution (CCE) algorithm.
        /// </summary>
        /// <param name="alpha">The worst point in the complex is reflected or contracted to seek an improvement alpha times.</param>
        /// <param name="Acomplex">The complex to evolve.</param>
        /// <param name="cdf">The trapezoidal cumulative probability for points in complex.</param>
        /// <param name="prng">The prng.</param>
        /// <param name="cancel">By ref. Determines if the solver should be canceled.</param>
        private void EvolveComplex(int alpha, ref List<PointFitness> Acomplex, double[] cdf, Random prng, ref bool cancel)
        {
            int beta = Acomplex.Count;
            int q = NumberOfParameters + 1;
            int i, j, k;
            var g = new double[NumberOfParameters]; // centroid g
            var p = new double[NumberOfParameters];
            double fitness;

            // The Beta loop. Allow sub-complex to evolve beta times.
            for (i = 0; i < beta; i++)
            {

                // Select parents by randomly choosing q distinct points from A complex
                var B = new List<PointFitness>();
                for (j = 1; j <= q; j++)
                {
                    double rnd = prng.NextDouble();
                    var rndi = Array.FindIndex(cdf, p => rnd <= p);
                    if (rndi == -1) rndi = cdf.Length - 1;

                    B.Add(Acomplex[rndi].Clone());
                    B.Last().Index = rndi;
                }

                // Alpha loop. The worst point in the complex is reflected or contracted to seek an improvement alpha times.
                for (j = 0; j < alpha; j++)
                {

                    // Rank the B points in terms of fitness (i.e., order of increasing function value.
                    B.Sort((x, y) => x.ParameterSet.Fitness.CompareTo(y.ParameterSet.Fitness));

                    // Find centroid g excluding worst point and compute reflection of worst point about centroid r = 2g - u(worst)
                    Reflection(ref B, ref g, ref p);

                    // Check if r is feasible
                    if (IsFeasible(p) == true)
                    {
                        // The parameter set is feasible, so evaluate it.
                        fitness = Evaluate(p, ref cancel);
                        if (cancel) return;
                    }
                    else
                    {
                        // The parameters were infeasible, so perform mutation
                        // Compute the smallest hypercube enclosing A complex and randomly sample a point within it. 
                        SmallestHypercube(ref Acomplex, ref prng, ref p);
                        fitness = Evaluate(p, ref cancel);
                        if (cancel) return;
                    }

                    // Now, either replace worst point with better point
                    if (fitness < B.Last().ParameterSet.Fitness)
                    {
                        B.Last().ParameterSet = new ParameterSet(p, fitness);
                    }
                    else
                    {
                        // Or contract to midpoint between centroid and worst point
                        Contraction(ref B, ref g, ref p);
                        fitness = Evaluate(p, ref cancel);
                        if (cancel) return;

                        // If better than worst point replace worst point
                        if (fitness < B.Last().ParameterSet.Fitness)
                        {
                            B.Last().ParameterSet = new ParameterSet(p, fitness);
                        }
                        else
                        {
                            // Otherwise perform mutation
                            // Compute the smallest hypercube enclosing A complex and randomly sample a point within it. 
                            SmallestHypercube(ref Acomplex, ref prng, ref p);
                            fitness = Evaluate(p, ref cancel);
                            if (cancel) return;
                            // Replace worst point with new point regardless of its value
                            B.Last().ParameterSet = new ParameterSet(p, fitness);
                        }
                    }
                }

                // Replace B into A according to L And sort A in order of increasing function value
                for (j = 0; j < Acomplex.Count; j++)
                {
                    for (k = 0; k < q; k++)
                    {
                        if (B[k].Index == j && B[k].ParameterSet.Fitness < Acomplex[j].ParameterSet.Fitness)
                        {
                            Acomplex[j].ParameterSet = B[k].ParameterSet.Clone();
                        }
                    }
                }

                Acomplex.Sort((x, y) => x.ParameterSet.Fitness.CompareTo(y.ParameterSet.Fitness));
            }
        }

        /// <summary>
        /// Determines if the point is within the feasible parameter space.
        /// </summary>
        /// <param name="point">The point to evaluate.</param>
        private bool IsFeasible(double[] point)
        {
            for (int i = 0; i < point.Length; i++)
                if (point[i] < LowerBounds[i] || point[i] > UpperBounds[i])
                    return false;
            return true;
        }

        /// <summary>
        /// Compute reflection of worst point about centroid r = 2g - u(worst)
        /// </summary>
        /// <param name="U">The sub-complex.</param>
        /// <param name="g">The centroid.</param>
        /// <param name="r">The reflection point.</param>
        private void Reflection(ref List<PointFitness> U, ref double[] g, ref double[] r)
        {
            for (int i = 0; i < NumberOfParameters; i++)
            {
                g[i] = 0.0d;
                for (int j = 0; j < U.Count - 1; j++)
                    g[i] += U[j].ParameterSet.Values[i];
                g[i] /= U.Count - 1;
                r[i] = 2d * g[i] - U.Last().ParameterSet.Values[i]; // u(worst)
            }
        }

        /// <summary>
        /// Contract to midpoint between centroid and worst point
        /// </summary>
        /// <param name="U">The sub-complex.</param>
        /// <param name="g">The centroid.</param>
        /// <param name="c">The contraction point.</param>
        private void Contraction(ref List<PointFitness> U, ref double[] g, ref double[] c)
        {
            for (int i = 0; i < NumberOfParameters; i++)
                c[i] = 0.5d * (g[i] + U.Last().ParameterSet.Values[i]); // u(worst)
        }

        /// <summary>
        /// Compute the smallest hypercube enclosing A complex and randomly sample a point within it.
        /// </summary>
        /// <param name="Acomplex">The complex.</param>
        /// <param name="prng">The prng.</param>
        /// <param name="z">The randomly generated point.</param>
        private void SmallestHypercube(ref List<PointFitness> Acomplex, ref Random prng, ref double[] z)
        {
            for (int i = 0; i < NumberOfParameters; i++)
            {
                double low = Acomplex[0].ParameterSet.Values[i];
                double high = low;
                for (int j = 1; j < Acomplex.Count; j++)
                {
                    low = Math.Min(low, Acomplex[j].ParameterSet.Values[i]);
                    high = Math.Max(high, Acomplex[j].ParameterSet.Values[i]);
                }
                z[i] = low + (high - low) * prng.NextDouble();
            }
        }

        /// <summary>
        /// The Trapezoidal cumulative probability distribution.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>The CDF.</returns>
        private double[] Trapezoidal(int N)
        {
            var PP = new double[N];
            for (int i = 1; i <= N; i++)
                PP[i - 1] = 2 * (N + 1 - i) / (N * (N + 1.0));
            for (int i = 1; i < N; i++)
                PP[i] += PP[i - 1]; // cumulative sum
            return PP;
        }

        /// <summary>
        /// Class for keeping track of parameter sets with indexes. 
        /// </summary>
        private class PointFitness
        {
            public ParameterSet ParameterSet;
            public int Index;
            public PointFitness Clone()
            {
                return new PointFitness { ParameterSet = ParameterSet.Clone(), Index = Index };
            }
        }

    }
}