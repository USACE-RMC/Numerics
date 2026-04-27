using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// Interface for Univariate Probability Distributions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface IUnivariateDistribution : IDistribution
    {

        /// <summary>
        /// Returns the univariate distribution type.
        /// </summary>
        UnivariateDistributionType Type { get; }

        /// <summary>
        /// Returns the number of distribution parameters.
        /// </summary>
        int NumberOfParameters { get; }

        /// <summary>
        /// Returns the distribution parameters in 2-column array of string.
        /// </summary>
        string[,] ParametersToString { get; }

        /// <summary>
        /// Returns the distribution parameter names in short form (e.g. µ, σ) in an array of string.
        /// </summary>
        string[] ParameterNamesShortForm { get; }

        /// <summary>
        /// Returns a short label of the distribution and parameters as a string.
        /// </summary>
        string DisplayLabel { get; }

        /// <summary>
        /// Returns the distribution parameters in an array of double.
        /// </summary>
        double[] GetParameters { get; }

        /// <summary>
        /// Returns the distribution parameter property names in an array of string.
        /// </summary>
        string[] GetParameterPropertyNames { get; }

        /// <summary>
        /// Returns a boolean value describing if the current parameters are valid or not.
        /// If not, an ArgumentOutOfRange exception will be thrown when trying to use statistical functions (e.g. CDF())
        /// </summary>
        bool ParametersValid { get; }

        /// <summary>
        /// Gets the mean of the distribution.
        /// </summary>
        double Mean { get; }

        /// <summary>
        /// Gets the median of the distribution.
        /// </summary>
        double Median { get; }

        /// <summary>
        /// Gets the mode of the distribution.
        /// </summary>
        double Mode { get; }

        /// <summary>
        /// Gets the variance of the distribution.
        /// </summary>
        double Variance { get; }

        /// <summary>
        /// Gets the standard deviation of the distribution.
        /// </summary>
        double StandardDeviation { get; }

        /// <summary>
        /// Gets the coefficient of variation of the distribution.
        /// </summary>
        double CoefficientOfVariation { get; }

        /// <summary>
        /// Gets the skewness of the distribution.
        /// </summary>
        double Skewness { get; }

        /// <summary>
        /// Gets the kurtosis of the distribution.
        /// </summary>
        double Kurtosis { get; }

        /// <summary>
        /// Gets the minimum of the distribution.
        /// </summary>
        double Minimum { get; }

        /// <summary>
        /// Gets the maximum of the distribution.
        /// </summary>
        double Maximum { get; }

        /// <summary>
        /// Gets the minimum values allowable for each parameter.
        /// </summary>
        double[] MinimumOfParameters { get; }

        /// <summary>
        /// Gets the maximum values allowable for each parameter.
        /// </summary>
        double[] MaximumOfParameters { get; }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="parameters">Array of parameters.</param>
        void SetParameters(IList<double> parameters);

        /// <summary>
        /// Validate the distribution parameters.
        /// </summary>
        /// <param name="parameters">Array of parameters.</param>
        /// <param name="throwException">Boolean indicating whether to throw the exception or not.</param>
        /// <returns>Nothing if the parameters are valid and the exception if invalid parameters were found.</returns>
        ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException);

        /// <summary>
        /// The Probability Density Function (PDF) of the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        /// <returns>The probability of X occurring.</returns>
        /// <remarks>
        /// The PDF describes the probability that X will occur. Returns the Probability Mass Function (PMF) for discrete distributions.
        /// </remarks>
        double PDF(double x);

        /// <summary>
        /// Returns the natural log of the PDF.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        double LogPDF(double x);

        /// <summary>
        /// Returns the hazard function (HF) of the distribution evaluated at a point x.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        double HF(double x);

        /// <summary>
        /// The Cumulative Distribution Function (CDF) for the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        /// <returns>The non-exceedance probability given a point X.</returns>
        /// <remarks>
        /// The CDF describes the cumulative probability that a given value or any value smaller than it will occur.
        /// </remarks>
        double CDF(double x);

        /// <summary>
        /// Returns the natural log of the CDF.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        double LogCDF(double x);

        /// <summary>
        /// The complement of the CDF. This function is also known as the survival function.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        double CCDF(double x);

        /// <summary>
        /// Returns the natural log of the CCDF.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        double LogCCDF(double x);

        /// <summary>
        /// The Inverse Cumulative Distribution Function (InvCFD) of the distribution evaluated at a probability.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <returns>
        /// Returns for a given probability in the probability distribution of a random variable,
        /// the value at which the probability of the random variable is less than or equal to the
        /// given probability.
        /// </returns>
        /// <remarks>
        /// This function is also know as the Quantile Function.
        /// </remarks>
        double InverseCDF(double probability);

        /// <summary>
        /// Generate random values of a distribution given a sample size.
        /// </summary>
        /// <param name="sampleSize"> Size of random sample to generate. </param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        double[] GenerateRandomValues(int sampleSize, int seed = -1);

    }
}