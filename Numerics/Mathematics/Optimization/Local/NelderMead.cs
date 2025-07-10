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
    /// Contains the Nelder-Mead downhill simplex algorithm, which finds a minima when no gradient is available.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// A simplex is the convex hull of n+1 vertices. In 2D space, a simplex is a triangle, and in 3D space, a simplex is a
    /// tetragedron. The Nelder-Mead algoirthm starts with a set of n+1 points that are the vertices of a working simplex, S,
    /// and the corresponding set of function values. Through a sequence of transformations of the working simplex S, the 
    /// function values at the vertices are decreased. At each step, the transformation is determined by computing one or more 
    /// test points, and their function values, and by comparison of these function values with those vertices. The process
    /// is terminated with the S becomes sufficiently small or the function values are close enough.
    /// </para>
    /// <para>
    /// <b>References: </b>
    /// <list type="bullet">
    /// <item><description> 
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </description></item>
    /// <item><description> 
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017
    /// </description></item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method"/>
    /// </description></item>
    /// <item><description> 
    /// <see href="http://www.scholarpedia.org/article/Nelder-Mead_algorithm"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class NelderMead : Optimizer
    {

        /// <summary>
        /// Construct a new Nelder-Mead optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values to evaluate.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public NelderMead(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
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
        /// The relative size of the initial simplex. Typical values range from 0.01 to 0.1.
        /// </summary>
        public double SimplexScale { get; set; } = 0.05;


        /// <inheritdoc/>
        protected override void Optimize()
        {

            // Define all variables
            // Three parameters which define the expansion and contractions
            double alpha = 1.0d;
            double beta = 0.5d;
            double gamma = 2.0d;

            int D = NumberOfParameters;
            bool cancel = false;
            int mpts = D + 1;
            var p = new double[mpts, D];
            var f = new double[mpts];
            var pr = new double[D];
            var prr = new double[D];
            var pbar = new double[D];
            int ilo, ihi, inhi;
            double fpr, fprr;

            // Define D + 1 vertices (+1 by row).
            // Scaled Simplex Initialization (per Nelder & Mead 1965 and Numerical Recipes)
            for (int i = 0; i < mpts; i++)
            {
                for (int j = 0; j < D; j++)
                    p[i, j] = InitialValues[j];
                //if (i != 0)
                //    p[i, i - 1] += (p[i, i - 1] == 0d) ? 0.00025 : p[i, i - 1] * 0.05;

                if (i > 0)
                {
                    double delta = SimplexScale * (InitialValues[i - 1] != 0.0 ? InitialValues[i - 1] : 1.0);
                    p[i, i - 1] += delta;
                    p[i, i - 1] = RepairParameter(p[i, i - 1], LowerBounds[i - 1], UpperBounds[i - 1]);
                }
            }

            // Initialize f to the values of the objective function evaluated 
            // at the D + 1 vertices (rows) of p.
            var PT = new double[D];
            for (int i = 0; i < mpts; i++)
            {
                for (int j = 0; j < D; j++)
                    PT[j] = p[i, j];
                f[i] = Evaluate(PT, ref cancel);
            }

            // Do Nelder-Mead Loop
            while (Iterations < MaxIterations)
            {
                // First we determine which point is the highest (worst),
                // next-highest, and lowest (best),
                ilo = 1;
                if (f[0] > f[1])
                {
                    ihi = 0;
                    inhi = 1;
                }
                else
                {
                    ihi = 1;
                    inhi = 0;
                }

                // by looping over the points in the simplex.
                for (int i = 0; i < mpts; i++)
                {
                    if (f[i] < f[ilo])
                        ilo = i;
                    if (f[i] > f[ihi])
                    {
                        inhi = ihi;
                        ihi = i;
                    }
                    else if (f[i] > f[inhi])
                    {
                        if (i != ihi)
                            inhi = i;
                    }
                }

                // Check convergence.
                if (CheckConvergence(f[ihi], f[ilo]))
                {
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }

                Iterations += 1;
                for (int j = 0; j < D; j++)
                    pbar[j] = 0d;

                // Begin a new iteration. Compute the vector average of all
                // points except the highest, i.e. the center of the "face" of the
                // simplex across from the high point. We will subsequently explore
                // along the ray from the high point through that center. 
                for (int i = 0; i < mpts; i++)
                {
                    if (i == ihi) continue;
                    for (int j = 0; j < D; j++)
                        pbar[j] += p[i, j];
                }

                // Extrapolate by a factor ALPHA through the face. i.e. reflect the
                // simplex from the high point.
                for (int j = 0; j < D; j++)
                {
                    pbar[j] = pbar[j] / D;
                    pr[j] = (1.0d + alpha) * pbar[j] - alpha * p[ihi, j];
                    // Make sure the parameter is within bounds
                    pr[j] = RepairParameter(pr[j], LowerBounds[j], UpperBounds[j]);
                }

                // Evaluate the function at the reflection point.
                fpr = Evaluate(pr, ref cancel);
                if (cancel) return;

                if (fpr <= f[ilo])
                {
                    // Gives a result better than the best point, so 
                    // try an additional extrapolation by a factor of GAMMA,
                    for (int j = 0; j < D; j++)
                    {
                        prr[j] = gamma * pr[j] + (1.0d - gamma) * pbar[j];
                        // Make sure the parameter is within bounds
                        prr[j] = RepairParameter(prr[j], LowerBounds[j], UpperBounds[j]);
                    }

                    // and check out the function there.
                    fprr = Evaluate(prr, ref cancel);
                    if (cancel) return; 

                    if (fprr < f[ilo])
                    {
                        // The additional extrapolation succeeded, and replaces
                        // the high point.
                        for (int j = 0; j < D; j++)
                            p[ihi, j] = prr[j];
                        f[ihi] = fprr;
                    }
                    else
                    {
                        // The additional extrapolation failed, but we 
                        // can still use the reflected point.
                        for (int j = 0; j < D; j++)
                            p[ihi, j] = pr[j];
                        f[ihi] = fpr;
                    }
                }
                else if (fpr >= f[inhi])
                {
                    // The reflected point is worse than the second-highest.
                    if (fpr < f[ihi])
                    {
                        // If it's better than the highest, the replace the highest,
                        for (int j = 0; j < D; j++)
                            p[ihi, j] = pr[j];
                        f[ihi] = fpr;
                    }
                    // but look for an intermediate lower point,
                    for (int j = 0; j < D; j++)
                    {
                        prr[j] = beta * p[ihi, j] + (1.0d - beta) * pbar[j];
                        // Make sure the parameter is within bounds
                        prr[j] = RepairParameter(prr[j], LowerBounds[j], UpperBounds[j]);
                    }
                    // in other words, perform a contraction of the simplex along
                    // one dimension. Then evaluate the function.
                    fprr = Evaluate(prr, ref cancel);
                    if (cancel) return;

                    if (fprr < f[ihi])
                    {
                        // Contraction gives an improvement, so accept it. 
                        for (int j = 0; j < D; j++)
                            p[ihi, j] = prr[j];
                        f[ihi] = fprr;
                    }
                    else
                    {
                        // Can't seem to get rid of the high point. Better contract
                        // around the lowest (best) point.
                        for (int i = 0; i < mpts; i++)
                        {
                            if (i == ilo) continue;
                            for (int j = 0; j < D; j++)
                            {
                                pr[j] = 0.5d * (p[i, j] + p[ilo, j]);
                                // Make sure the parameter is within bounds
                                pr[j] = RepairParameter(pr[j], LowerBounds[j], UpperBounds[j]);
                                p[i, j] = pr[j];
                            }

                            f[i] = Evaluate(pr, ref cancel);
                            if (cancel) return;
                         
                        }
                    }
                }
                else
                {
                    // We arrive here if the original reflection gives a middling point.
                    // Replace the old high point and continue
                    for (int j = 0; j < D; j++)
                        p[ihi, j] = pr[j];
                    f[ihi] = fpr;
                }
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);

        }

    }
}