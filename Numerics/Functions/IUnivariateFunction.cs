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
        ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException);

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
