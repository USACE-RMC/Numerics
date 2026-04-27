using System;

namespace Numerics.Distributions
{

    /// <summary>
    /// Interface for Multivariate Probability Distributions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface IMultivariateDistribution : IDistribution
    {

        /// <summary>
        /// Gets the number of variables for the distribution.
        /// </summary>
        int Dimension { get; }

        /// <summary>
        /// Returns the multivariate distribution type.
        /// </summary>
        MultivariateDistributionType Type { get; }

        /// <summary>
        /// Returns a boolean value describing if the current parameters are valid or not.
        /// If not, an ArgumentOutOfRange exception will be thrown when trying to use statistical functions (e.g. CDF())
        /// </summary>
        bool ParametersValid { get; }

        /// <summary>
        /// The Probability Density Function (PDF) of the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        double PDF(double[] x);

        /// <summary>
        /// Returns the natural log of the PDF.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        double LogPDF(double[] x);

        /// <summary>
        /// The Cumulative Distribution Function (CDF) for the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        double CDF(double[] x);

        /// <summary>
        /// Returns the natural log of the CDF.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        double LogCDF(double[] x);

        /// <summary>
        /// The complement of the CDF. This function is also known as the survival function.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        double CCDF(double[] x);

        /// <summary>
        /// Returns the natural log of the CCDF.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        double LogCCDF(double[] x);

    }
}