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