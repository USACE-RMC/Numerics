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

namespace Numerics.Sampling.MCMC
{

    /// <summary>
    /// A class for storing parameter statistics.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class ParameterStatistics
    {

        /// <summary>
        /// Constructs an empty parameter summary statistics class.
        /// </summary>
        public ParameterStatistics() { }

        /// <summary>
        /// The Gelman-Rubin diagnostic.
        /// </summary>
        public double Rhat { get; set; }

        /// <summary>
        /// The effective sample size.
        /// </summary>
        public double ESS { get; set; }

        /// <summary>
        /// The total sample size.
        /// </summary>
        public int N { get; set; }

        /// <summary>
        /// The parameter mean.
        /// </summary>
        public double Mean { get; set; }

        /// <summary>
        /// The parameter median.
        /// </summary>
        public double Median { get; set; }

        /// <summary>
        /// The parameter standard deviation.
        /// </summary>
        public double StandardDeviation { get; set; }

        /// <summary>
        /// The lower confidence interval.
        /// </summary>
        public double LowerCI { get; set; }

        /// <summary>
        /// The upper confidence interval.
        /// </summary>
        public double UpperCI { get; set; }

    }
}
