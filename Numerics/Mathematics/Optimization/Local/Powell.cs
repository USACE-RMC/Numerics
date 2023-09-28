﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{
    /// <summary>
    /// Contains the Powell optimization algorithm. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// References:
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017.
    /// <see href="https://en.wikipedia.org/wiki/Powell%27s_method"/>
    /// </remarks>
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

        /// <summary>
        /// Implements the actual optimization algorithm. This method should minimize the objective function. 
        /// </summary>
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
