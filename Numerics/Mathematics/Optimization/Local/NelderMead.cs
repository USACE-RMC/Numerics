﻿using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// Contains the Nelder-Mead downhill simplex algorithm, which finds a minima when no gradient is available.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// References:
    /// </para>
    /// <para>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017.
    /// <see href="https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method"/>
    /// </para>
    /// </remarks>
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
        /// Implements the actual optimization algorithm. This method should minimize the objective function. 
        /// </summary>
        protected override void Optimize()
        {

            // Define all variables
            // Three parameters which define the expansion and contractions
            double alpha = 1.0d;
            double beta = 0.5d;
            double gamma = 2.0d;

            int i, j, D = NumberOfParameters;
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
            // This provides a heuristic way to initialize the simplex.
            // http://www.scholarpedia.org/article/Nelder-Mead_algorithm
            for (i = 0; i < mpts; i++)
            {
                for (j = 0; j < D; j++)
                    p[i, j] = InitialValues[j];
                if (i != 0)
                {
                    if (p[i, i - 1] == 0d)
                    {
                        p[i, i - 1] += 0.00025d;
                    }
                    else
                    {
                        p[i, i - 1] += p[i, i - 1] * 0.05d;
                    }
                }
            }
            // Initialize f to the values of the objective function evaluated 
            // at the D + 1 vertices (rows) of p.
            var PT = new double[D];
            for (i = 0; i < mpts; i++)
            {
                for (j = 0; j < D; j++)
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
                for (i = 0; i < mpts; i++)
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
                for (j = 0; j < D; j++)
                    pbar[j] = 0d;

                // Begin a new iteration. Compute the vector average of all
                // points except the highest, i.e. the center of the "face" of the
                // simplex across from the high point. We will subsequently explore
                // along the ray from the high point through that center. 
                for (i = 0; i < mpts; i++)
                {
                    if (i != ihi)
                    {
                        for (j = 0; j < D; j++)
                            pbar[j] += p[i, j];
                    }
                }

                // Extrapolate by a factor ALPHA through the face. i.e. reflect the
                // simplex from the high point.
                for (j = 0; j < D; j++)
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
                    for (j = 0; j < D; j++)
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
                        for (j = 0; j < D; j++)
                            p[ihi, j] = prr[j];
                        f[ihi] = fprr;
                    }
                    else
                    {
                        // The additional extrapolation failed, but we 
                        // can still use the reflected point.
                        for (j = 0; j < D; j++)
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
                        for (j = 0; j < D; j++)
                            p[ihi, j] = pr[j];
                        f[ihi] = fpr;
                    }
                    // but look for an intermediate lower point,
                    for (j = 0; j < D; j++)
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
                        for (j = 0; j < D; j++)
                            p[ihi, j] = prr[j];
                        f[ihi] = fprr;
                    }
                    else
                    {
                        // Can't seem to get rid of the high point. Better contract
                        // around the lowest (best) point.
                        for (i = 0; i < mpts; i++)
                        {
                            if (i != ilo)
                            {
                                for (j = 0; j < D; j++)
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
                }
                else
                {
                    // We arrive here if the original reflection gives a middling point.
                    // Replace the old high point and continue
                    for (j = 0; j < D; j++)
                        p[ihi, j] = pr[j];
                    f[ihi] = fpr;
                }
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);

        }

    }
}