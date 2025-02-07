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

using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Runtime.Remoting.Messaging;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// A class for keeping track of a running covariance matrix.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class RunningCovarianceMatrix
    {
        /// <summary>
        /// Construct a covariance matrix.
        /// </summary>
        /// <param name="size">The number of rows and columns.</param>
        public RunningCovarianceMatrix(int size)
        {
            // The mean vector set at zero.
            Mean = new Matrix(size, 1);
            // The initial covariance  is the identity matrix
            Covariance = Matrix.Identity(size);
        }

        /// <summary>
        /// The sample size N.
        /// </summary>
        public int N { get; private set; }

        /// <summary>
        /// The mean vector
        /// </summary>
        public Matrix Mean { get; private set; }

        /// <summary>
        /// The covariance matrix. This is unadjusted by the sample size.
        /// </summary>
        public Matrix Covariance { get; private set; }

        /// <summary>
        /// The sample covariance matrix corrected by the sample size with the degrees of freedom adjustment.
        /// </summary>
        public Matrix SampleCovariance => (1d / (N - 1)) * Covariance;

        /// <summary>
        /// The population covariance matrix corrected by the total sample size.
        /// </summary>
        public Matrix PopulationCovariance => (1d / N) * Covariance;

        /// <summary>
        /// The sample correlation matrix.
        /// </summary>
        public Matrix SampleCorrelation
        { 
            get
            {
                var D = Matrix.Diagonal(SampleCovariance);
                D.Sqrt();
                var invSqrtD = !D;
                return (invSqrtD * SampleCovariance) * invSqrtD;
            }       
        }

        /// <summary>
        /// The population correlation matrix.
        /// </summary>
        public Matrix PopulationCorrelation
        {
            get
            {
                var D = Matrix.Diagonal(PopulationCovariance);
                D.Sqrt();
                var invSqrtD = !D;
                return (invSqrtD * PopulationCovariance) * invSqrtD;
            }
        }

        /// <summary>
        /// Add a new vector to the running statistics. 
        /// </summary>
        /// <param name="values">Vector of data values. The length of the vector but be the same as 
        /// the number of rows</param>
        public void Push(IList<double> values)
        {
            N += 1;
            var x = new Matrix(values.ToArray());
            // Update mean
            var oldMean = Mean;
            Mean += N == 1 ? x : (x - Mean) * (1d / N);
            // Update covariance                 
            Covariance += (x - oldMean) * Matrix.Transpose(x - Mean);
        }

    }
}
