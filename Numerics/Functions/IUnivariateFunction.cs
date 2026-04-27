using System;
using System.Collections.Generic;

namespace Numerics.Functions
{
    /// <summary>
    /// Interface for Univariate Functions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors:  </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public interface IUnivariateFunction
    {
        /// <summary>
        /// Returns the number of function parameters.
        /// </summary>
        int NumberOfParameters { get; }

        /// <summary>
        /// Returns a boolean value describing if the current parameters are valid or not.
        /// If not, an ArgumentOutOfRange exception will be thrown when trying to use function.
        /// </summary>
        bool ParametersValid { get; }

        /// <summary>
        /// Gets and sets the minimum X value supported by the function. Default = double.MinValue.
        /// </summary>
        double Minimum { get; set;  }

        /// <summary>
        /// Gets and sets the maximum X value supported by the function. Default = double.MaxValue.
        /// </summary>
        double Maximum { get; set;  }

        /// <summary>
        /// Gets the minimum values allowable for each parameter.
        /// </summary>
        double[] MinimumOfParameters { get; }

        /// <summary>
        /// Gets the maximum values allowable for each parameter.
        /// </summary>
        double[] MaximumOfParameters { get; }

        /// <summary>
        /// Determines if the function is deterministic or if it has uncertainty. 
        /// </summary>
        bool IsDeterministic { get; set; }

        /// <summary>
        /// The confidence level to estimate when the function has uncertainty. 
        /// </summary>
        double ConfidenceLevel { get; set; }

        /// <summary>
        /// Set the function parameters.
        /// </summary>
        /// <param name="parameters">Array of parameters.</param>
        void SetParameters(IList<double> parameters);

        /// <summary>
        /// Test to see if function parameters are valid.
        /// </summary>
        /// <param name="parameters">Array of parameters.</param>
        /// <param name="throwException">Boolean indicating whether to throw the exception or not.</param>
        /// <returns>Nothing if the parameters are valid and the exception if invalid parameters were found.</returns>
        ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException);

        /// <summary>
        /// Returns the function evaluated at a point x. If function is uncertain, the function is computed at the set confidence level. 
        /// </summary>
        /// <param name="x">The x-value in the function to evaluate.</param>
        double Function(double x);

        /// <summary>
        /// Returns the inverse function evaluated at a point y. If function is uncertain, the function is computed at the set confidence level. 
        /// </summary>
        /// <param name="y">The y-value in the inverse function to evaluate.</param>
        double InverseFunction(double y);

    }
}
