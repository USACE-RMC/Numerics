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

namespace Numerics.MachineLearning
{
    /// <summary>
    /// Supporting class for a Jenks natural breaks cluster.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The Jenks optimization method, also called the Jenks natural breaks classification method, 
    /// is a data clustering method designed to determine the best arrangement of values into different classes.
    /// </para>
    /// </remarks>
    public class JenksCluster
    {
        /// <summary>
        /// Creates a new Jenks cluster. 
        /// </summary>
        /// <param name="data">The sorted input data array.</param>
        /// <param name="startIndex">The starting index of the cluster.</param>
        /// <param name="endIndex">The ending index of the cluster.</param>
        public JenksCluster(double[] data, int startIndex,  int endIndex)
        {
            StartIndex = startIndex;
            EndIndex = endIndex;
            MinValue = data[startIndex];
            MaxValue = data[endIndex];

            // Compute summary statistics
            double X = 0;     // sum
            double X2 = 0;    // sum of X^2
            for (int i = startIndex; i <= endIndex; i++) 
            {
                X += data[i];
                X2 += Math.Pow(data[i], 2d);
            }

            double U1 = X / Count;
            double U2 = X2 / Count;
            Sum = X;
            Average = Count == 1 ? X : U1;
            Variance = Count == 1 ? 0 : (U2 - Math.Pow(U1, 2d)) * (Count / (double)(Count - 1));
            SumOfSquaredDeviations = Variance * (double)(Count - 1);
        }

        /// <summary>
        /// The starting index of the cluster.
        /// </summary>
        public int StartIndex { get; private set; }

        /// <summary>
        /// The ending index of the cluster.
        /// </summary>
        public int EndIndex { get; private set; }

        /// <summary>
        /// The number of data points in the cluster.
        /// </summary>
        public int Count { get { return EndIndex - StartIndex + 1; } }

        /// <summary>
        /// The minimum value of the cluster.
        /// </summary>
        public double MinValue { get; private set; }

        /// <summary>
        /// The maximum value of the cluster.
        /// </summary>
        public double MaxValue { get; private set; }

        /// <summary>
        /// The sum of the values in the cluster.
        /// </summary>
        public double Sum { get; private set; }

        /// <summary>
        /// The average value of the cluster. 
        /// </summary>
        public double Average { get; private set; }

        /// <summary>
        /// The variance of the cluster.
        /// </summary>
        public double Variance { get; private set; }

        /// <summary>
        /// The sum of squared deviations.
        /// </summary>
        public double SumOfSquaredDeviations { get; private set; }

    }
}
