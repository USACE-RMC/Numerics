using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// An interface for calculating the standard error for a probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface IStandardError
    {

        /// <summary>
        /// Returns a matrix containing the covariances of the parameters given the sample size.
        /// </summary>
        /// <param name="sampleSize">The sample size.</param>
        /// <param name="estimationMethod">The distribution parameter estimation method.</param>
        double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod);

        /// <summary>
        /// Returns the quantile variance given probability and sample size.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <param name="sampleSize">The sample size.</param>
        /// <param name="estimationMethod">The distribution parameter estimation method.</param>
        double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod);

        /// <summary>
        /// Returns a list of partial derivatives of X given probability with respect to each parameter.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        double[] QuantileGradient(double probability);

        /// <summary>
        /// Returns the Jacobian matrix of the quantile function with respect to each parameter.
        /// </summary>
        /// <param name="probabilities">List of probabilities, must be the same length as the number of distribution parameters.</param>
        /// <param name="determinant">Output. The determinant of the Jacobian matrix.</param>
        double[,] QuantileJacobian(IList<double> probabilities, out double determinant);

    }
}