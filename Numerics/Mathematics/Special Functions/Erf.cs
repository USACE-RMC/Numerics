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

using Numerics.Distributions;

namespace Numerics.Mathematics.SpecialFunctions
{
    /// <summary>
    /// The error function.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The error function, also known as the Gauss error function, is defined 
    /// as:
    /// </para>
    /// <code>
    ///                        X
    ///     erf(X) = 2/sqrt(𝜋) ∫  e^(-t^2) dt
    ///                        0           
    /// </code>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Error_function#cite_note-2"/>
    ///  </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public partial class Erf
    {
        /// <summary>
        /// Computes the error function
        /// </summary>
        /// <param name="X"> The upper bound </param>>
        /// <returns>
        /// The error function evaluated with the given upper bound
        /// </returns>
        public static double Function(double X)
        {
            if (X < 0d)
            {
                return -Gamma.LowerIncomplete(0.5d, X * X);
            }
            else
            {
                return Gamma.LowerIncomplete(0.5d, X * X);
            }
        }

        /// <summary>
        /// Computes the complement of the error function
        /// </summary>
        /// <param name="X"> The upper bound </param>>
        /// <returns>
        /// The complement of the error function evaluated with the given upper bound
        /// </returns>
        public static double Erfc(double X)
        {
            return 1d - Function(X);
        }

        /// <summary>
        /// Computes the inverse error function
        /// </summary>
        ///<param name="y"> The value to be evaluated (such that y = erf(erf^-1(y)) ) </param>
        /// <returns>
        /// The inverse error function evaluated at the given y
        /// </returns>
        public static double InverseErf(double y)
        {
            double s = Normal.StandardZ(0.5d * y + 0.5d);
            double r = s * Tools.Sqrt2 / 2.0d;
            return r;
        }

        /// <summary>
        /// Computes the inverse of the complement of the error function
        /// </summary>
        ///<param name="y"> The value to be evaluated (such that y = erf(erf^-1(y)) ) </param>
        /// <returns>
        /// The inverse of the complement of the error function evaluated at the given y
        /// </returns>
        public static double InverseErfc(double y)
        {
            double s = Normal.StandardZ(-0.5d * y + 1d);
            double r = s * Tools.Sqrt2 / 2.0d;
            return r;
        }
    }
}