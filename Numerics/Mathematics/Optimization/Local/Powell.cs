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
    /// Contains the Powell optimization algorithm. The function need not be differentiable, and no derivatives are taken.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors:</b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This method minimizes the function by bi-directionally searching along search vectors via Brent's method. The minima of
    /// these search vectors are recorded and used to create new search vectors, as the current ones are deleted after use. 
    /// The algorithm iterates until no significant improvement is made.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description> 
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </description></item>
    /// <item><description> 
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017
    /// </description></item>
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Powell%27s_method"/>
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class Powell : Optimizer
    {
        /// <summary>
        /// Construct a new Powell optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="numberOfParameters">The number of parameters in the objective function.</param>
        /// <param name="initialValues">An array of initial values to evaluate.</param>
        /// <param name="lowerBounds">An array of lower bounds (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBounds">An array of upper bounds (inclusive) of the interval containing the optimal point.</param>
        public Powell(Func<double[], double> objectiveFunction, int numberOfParameters, IList<double> initialValues, IList<double> lowerBounds, IList<double> upperBounds) : base(objectiveFunction, numberOfParameters)
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

        /// <inheritdoc/>
        protected override void Optimize()
        {
            // Set variables
            int i, j, D = NumberOfParameters, ibig;
            bool cancel = false;
            double t, fret, fp, fptt, delta;
            var p = InitialValues.ToArray();
            var pt = new double[D];
            var ptt = new double[D];
            var xi = new double[D]; // Direction vector
            // Set the initial matrix for directions
            // and save the initial point
            var ximat = new double[D, D];
            for (i = 0; i < D; i++)
            {
                ximat[i, i] = 1d;
                pt[i] = p[i];
            }
            // initial function evaluation
            fret = Evaluate(p, ref cancel);
            while (Iterations < MaxIterations)
            {
                fp = fret;
                ibig = 0;
                delta = 0.0; // Will be the biggest function decrease.
                // In each iteration, loop over all directions in the set.
                for (i = 0; i < D; i++)
                {
                    // Copy the direction
                    for (j = 0; j < D; j++) xi[j] = ximat[j, i];
                    fptt = fret;
                    fret = LineMinimization(p, xi, ref cancel);
                    if (cancel == true) return;
                    // And record it if it is the larges decrease so far.
                    if (fptt - fret > delta)
                    {
                        delta = fptt - fret;
                        ibig = i + 1;
                    }
                }
                // Check convergence
                if (CheckConvergence(fp, fret))
                {
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }
                // Construct the extrapolated point and save the average direction moved.
                // Save the old starting point.
                for (j = 0; j < D; j++)
                {
                    ptt[j] = 2.0 * p[j] - pt[j];
                    xi[j] = p[j] - pt[j];
                    pt[j] = p[j];
                }
                // Function evaluated at the extrapolated point
                fptt = Evaluate(ptt, ref cancel);
                if (cancel == true) return;
                if (fptt < fp)
                {
                    t = 2.0 * (fp - 2.0 * fret + fptt) * Tools.Sqr(fp - fret - delta) - delta * Tools.Sqr(fp - fptt);
                    if (t < 0.0)
                    {
                        // Move to the minimum of the new direction and save the new direction
                        fret = LineMinimization(p, xi, ref cancel);
                        if (cancel == true) return;
                        for (j = 0; j < D; j++)
                        {
                            ximat[j, ibig - 1] = ximat[j, D - 1];
                            ximat[j, D - 1] = xi[j];
                        }
                    }
                }

                Iterations += 1;
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);

        }

        /// <summary>
        /// Auxiliary line minimization routine. 
        /// </summary>
        /// <param name="startPoint">The initial point.</param>
        /// <param name="direction">The initial direction.</param>
        /// <param name="cancel">Determines if the solver should be canceled.</param>
        private double LineMinimization(double[] startPoint, double[] direction, ref bool cancel)
        {
            // Line-minimization routine, Given an n-dimensional point p[0..n-1] and an n-dimension 
            // direction xi[0..n-1], moves and resets p to where the function of functor func(p) takes on
            // a minimum along the direction xi from p, and replaces xi by the actual vector displacement
            // that p was moved. Also returns the value of func at the return location p. This is actually
            // all accomplished by calling the Brent minimize routine. 
            int D = NumberOfParameters;
            bool c = cancel;
            double func(double alpha)
            {
                var x = new double[D];
                for (int i = 0; i < D; i++)
                    x[i] = startPoint[i] + alpha * direction[i];
                return Evaluate(x, ref c);
            }
            var brent = new BrentSearch(func, 0d, 1d) { RelativeTolerance = RelativeTolerance, AbsoluteTolerance = AbsoluteTolerance };
            brent.Bracket(0.1);
            brent.Minimize();
            cancel = c;
            if (cancel) return double.NaN;
            double xmin = brent.BestParameterSet.Values[0];
            for (int j = 0; j < NumberOfParameters; j++)
            {
                direction[j] *= xmin;
                startPoint[j] += direction[j];
                // Make sure the parameter is within bounds
                startPoint[j] = RepairParameter(startPoint[j], LowerBounds[j], UpperBounds[j]);
            }
            return brent.BestParameterSet.Fitness;
        }

    }
}
