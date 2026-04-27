using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// An interface for Maximum Likelihood Estimation.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface IMaximumLikelihoodEstimation
    {

        /// <summary>
        /// Get the initial, lower, and upper values for the distribution parameters for constrained optimization.
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        /// <returns>Returns a Tuple of initial, lower, and upper values.</returns>
        Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample);

        /// <summary>
        /// Estimates the distribution parameters using the method of Maximum Likelihood Estimation (MLE).
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        /// <returns>Returns the MLE parameter set.</returns>
        double[] MLE(IList<double> sample);
    }
}