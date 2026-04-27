using System;

namespace Numerics.Distributions
{

    /// <summary>
    /// Declares common functionality for Multivariate Probability Distributions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class MultivariateDistribution : IMultivariateDistribution
    {

        /// <summary>
        /// Gets the number of variables for the distribution.
        /// </summary>
        public abstract int Dimension { get; }

        /// <summary>
        /// Returns the multivariate distribution type.
        /// </summary>
        public abstract MultivariateDistributionType Type { get; }

        /// <summary>
        /// Returns the display name of the distribution type as a string.
        /// </summary>
        public abstract string DisplayName { get; }

        /// <summary>
        /// Returns the short display name of the distribution as a string.
        /// </summary>
        public abstract string ShortDisplayName { get; }

        /// <summary>
        /// Returns a boolean value describing if the current parameters are valid or not.
        /// If not, an ArgumentOutOfRange exception will be thrown when trying to use statistical functions (e.g. CDF())
        /// </summary>
        public abstract bool ParametersValid { get; }
       
        /// <summary>
        /// The Probability Density Function (PDF) of the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        public abstract double PDF(double[] x);

        /// <summary>
        /// Returns the natural log of the PDF.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        public virtual double LogPDF(double[] x)
        {
            double f = PDF(x);
            // If the PDF returns an invalid probability, then return the worst log-probability.
            if (double.IsNaN(f) || double.IsInfinity(f) || f <= 0d) return double.MinValue;
            return Math.Log(f);
        }

        /// <summary>
        /// The Cumulative Distribution Function (CDF) for the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        public abstract double CDF(double[] x);

        /// <summary>
        /// Returns the natural log of the CDF.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        public virtual double LogCDF(double[] x)
        {
            double F = CDF(x);
            // If the CDF returns an invalid probability, then return the worst log-probability.
            if (double.IsNaN(F) || double.IsInfinity(F) || F <= 0d) return double.MinValue;
            return Math.Log(F);
        }

        /// <summary>
        /// The complement of the CDF. This function is also known as the survival function.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        public virtual double CCDF(double[] x)
        {
            return 1.0d - CDF(x);
        }

        /// <summary>
        /// Returns the natural log of the CCDF.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        public virtual double LogCCDF(double[] x)
        {
            double cF = CCDF(x);
            // If the CCDF returns an invalid probability, then return the worst log-probability.
            if (double.IsNaN(cF) || cF <= 0d)
                return int.MinValue;
            return Math.Log(cF);
        }

        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public abstract MultivariateDistribution Clone();

    }
}