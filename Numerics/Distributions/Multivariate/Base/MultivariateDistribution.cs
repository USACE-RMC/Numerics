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