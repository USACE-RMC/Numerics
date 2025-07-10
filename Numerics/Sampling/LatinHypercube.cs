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

namespace Numerics.Sampling
{
    /// <summary>
    /// A class to perform Latin hypercube sampling (LHS).
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
    /// Stratified sampling is used to generate uniform random numbers by dividing the interval [0,1) into n bins 
    /// of equal probability, where n is the total number of samples required. 
    /// <para>
    /// In each iteration:
    /// </para>
    /// <list type="number">
    /// <item>
    /// A random number is generated to select one of the remaining bins.
    /// </item>
    /// <item>
    /// A second random number is generated to select a value within the chosen bin.
    /// </item>
    /// <item>
    /// The selected bin is marked as used and will not be selected in subsequent iterations.
    /// </item>
    /// </list>
    /// <para>
    /// This process is repeated until all n bins have been sampled, ensuring that each bin is selected exactly once.
    /// </para>
    /// </para>
    /// </remarks>
    public class LatinHypercube
    {

        /// <summary>
        /// Generate Latin Hypercube samples with uniform stratified sampling.
        /// </summary>
        /// <param name="sampleSize">Number of samples (rows).</param>
        /// <param name="dimension">Number of dimensions (columns).</param>
        /// <param name="seed">Random seed (optional).</param>
        public static double[,] Random(int sampleSize, int dimension = 1, int seed = -1)
        {
            double[,] lhs = new double[sampleSize, dimension];
            var rndMaster = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            for (int col = 0; col < dimension; col++)
            {
                // Create bins
                double[] bins = new double[sampleSize];
                var rnd = new MersenneTwister(rndMaster.Next());

                for (int i = 0; i < sampleSize; i++)
                {
                    double delta = rnd.NextDouble();
                    bins[i] = (i + delta) / sampleSize;
                }

                // Shuffle bins
                Shuffle(bins, rnd);

                // Assign to LHS matrix
                for (int row = 0; row < sampleSize; row++)
                    lhs[row, col] = bins[row];
            }

            return lhs;
        }

        /// <summary>
        /// Generate Latin Hypercube samples using median bin locations.
        /// </summary>
        /// <param name="sampleSize">Number of samples (rows).</param>
        /// <param name="dimension">Number of dimensions (columns).</param>
        /// <param name="seed">Random seed (optional).</param>
        public static double[,] Median(int sampleSize, int dimension = 1, int seed = -1)
        {
            double[,] lhs = new double[sampleSize, dimension];
            var rndMaster = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            for (int col = 0; col < dimension; col++)
            {
                // Create median-centered bins
                double[] bins = new double[sampleSize];
                for (int i = 0; i < sampleSize; i++)
                    bins[i] = (i + 0.5) / sampleSize;

                // Shuffle bins
                Shuffle(bins, new MersenneTwister(rndMaster.Next()));

                // Assign to LHS matrix
                for (int row = 0; row < sampleSize; row++)
                    lhs[row, col] = bins[row];
            }

            return lhs;
        }

        /// <summary>
        /// Performs in-place Fisher–Yates shuffle.
        /// </summary>
        private static void Shuffle(double[] array, Random rnd)
        {
            for (int i = array.Length - 1; i > 0; i--)
            {
                int j = rnd.Next(i + 1);
                (array[i], array[j]) = (array[j], array[i]);
            }
        }
    }
}
