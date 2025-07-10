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
using System.Linq;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// Jenks natural breaks optimization.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors:</b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The Jenks optimization method, also called the Jenks natural breaks classification method, 
    /// is a data clustering method designed to determine the best arrangement of values into different classes.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> <see href = "http://en.wikipedia.org/wiki/Jenks_natural_breaks_optimization" /> </item>
    /// <item> <see href = "https://github.com/simple-statistics/simple-statistics/blob/main/src/jenks.js" /> </item>
    /// <item> <see href="https://medium.com/analytics-vidhya/jenks-natural-breaks-best-range-finder-algorithm-8d1907192051"/></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class JenksNaturalBreaks
    {
        /// <summary>
        /// Creates a new Jenks natural breaks optimization.
        /// </summary>
        /// <param name="data">The input data array to be classified.</param>
        /// <param name="numberOfClusters">The number of desired clusters (or classes).</param>
        /// <param name="isDataSorted">Determines if the data array is sorted. Default = false.</param>
        public JenksNaturalBreaks(IList<double> data, int numberOfClusters, bool isDataSorted = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data), "The data array is null.");
            if (data.Count == 0) throw new ArgumentException("The data array is empty.", nameof(data));
            if (numberOfClusters <= 0) throw new ArgumentNullException(nameof(numberOfClusters), "The number of clusters must be greater than zero.");
            if (numberOfClusters > data.Count) throw new ArgumentException("The number of clusters cannot be greater than the length of the data array.", nameof(numberOfClusters));

            // Sort the data in numerical order
            double[] sortedData = data.ToArray();
            if (isDataSorted == false)
            {
                Array.Sort(sortedData);
            }

            SortedData = sortedData;
            NumberOfClusters = numberOfClusters;
            Estimate();
        }

        /// <summary>
        /// Creates a new Jenks natural breaks optimization.
        /// </summary>
        /// <param name="data">The input data array to be classified.</param>
        /// <param name="numberOfClusters">The number of desired clusters (or classes).</param>
        /// <param name="isDataSorted">Determines if the data array is sorted. Default = false.</param>
        public JenksNaturalBreaks(IList<float> data, int numberOfClusters, bool isDataSorted = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data), "The data array is null.");
            if (data.Count == 0) throw new ArgumentException("The data array is empty.", nameof(data));
            if (numberOfClusters <= 0) throw new ArgumentNullException(nameof(numberOfClusters), "The number of clusters must be greater than zero.");
            if (numberOfClusters > data.Count) throw new ArgumentException("The number of clusters cannot be greater than the length of the data array.", nameof(numberOfClusters));

            // Sort the data in numerical order
            float[] sortedData = data.ToArray();
            if (isDataSorted == false)
            {
                Array.Sort(sortedData);
            }

            SortedData = new double[sortedData.Length];
            for (int i = 0; i < sortedData.Length; i++)
            {
                SortedData[i] = sortedData[i];
            }

            NumberOfClusters = numberOfClusters;
            Estimate();
        }

        /// <summary>
        /// The array of sorted input data.
        /// </summary>
        public double[] SortedData { get; private set; }

        /// <summary>
        /// The number of clusters.
        /// </summary>
        public int NumberOfClusters { get; private set; }

        /// <summary>
        /// Gets the array of estimated clusters.
        /// </summary>
        public JenksCluster[] Clusters { get; private set; }

        /// <summary>
        /// The array of break points.
        /// </summary>
        public double[] Breaks { get; private set; }

        /// <summary>
        /// The goodness of fit measure. The closer to 1, the better the fit.
        /// </summary>
        public double GoodnessOfVarianceFit 
        {
            get 
            {
                double gvf = 0;
                double mean = Tools.Mean(SortedData);
                // Sum of squared deviations from data
                double sdam = 0;
                for (int i = 0; i < SortedData.Length; i++)
                {
                    sdam += Tools.Sqr(SortedData[i] - mean);
                }
                // Sum of squared deviations from clusters
                double sdcm = 0;
                for (int i = 0; i < Clusters.Length; i++)
                {
                    sdcm += Clusters[i].SumOfSquaredDeviations;
                }
                gvf = (sdam - sdcm) / sdam;
                return gvf;

            }
        }

        /// <summary>
        /// Estimate the Jenks natural breaks. 
        /// </summary>
        private void Estimate()
        {
            // Initialize the matrices
            var lowerClassLimits = new int[SortedData.Length + 1, NumberOfClusters + 1];
            var varianceCombinations = new double[SortedData.Length + 1, NumberOfClusters + 1];
            for (int i = 1; i <= NumberOfClusters; i++)
            {
                lowerClassLimits[1, i] = 1;
                varianceCombinations[1, i] = 0;
                for (int j = 2; j <= SortedData.Length; j++)
                    varianceCombinations[j, i] = double.MaxValue;
            }

            double v = 0;
            for (int l = 2; l <= SortedData.Length; l++)
            {
                double s1 = 0;
                double s2 = 0;
                double w = 0;
                for (int m = 1; m <= l; m++)
                {
                    int i3 = l - m + 1;
                    double val = SortedData[i3 - 1];

                    s2 += val * val;
                    s1 += val;

                    w++;
                    v = s2 - (s1 * s1) / w;
                    int i4 = i3 - 1;
                    if (i4 != 0)
                    {
                        for (int j = 2; j <= NumberOfClusters; j++)
                        {
                            if (varianceCombinations[l, j] >= (v + varianceCombinations[i4, j - 1]))
                            {
                                lowerClassLimits[l, j] = i3;
                                varianceCombinations[l, j] = v + varianceCombinations[i4, j - 1];

                            };
                        };
                    };
                };
                lowerClassLimits[l, 1] = 1;
                varianceCombinations[l, 1] = v;
            };

            // Create clusters
            int k = SortedData.Length;
            int[] kclass = new int[NumberOfClusters];
            kclass[NumberOfClusters - 1] = SortedData.Length - 1;
            for (int j = NumberOfClusters; j >= 2; j--)
            {
                int id = lowerClassLimits[k, j] - 2;
                kclass[j - 2] = id;
                k = lowerClassLimits[k, j] - 1;
            };

            Clusters = new JenksCluster[NumberOfClusters];
            Breaks = new double[NumberOfClusters];
            Clusters[0] = new JenksCluster(SortedData, 0, kclass[0]);
            Breaks[0] = Clusters[0].MaxValue;
            for (int i = 1; i < NumberOfClusters; i++)
            {
                Clusters[i] = new JenksCluster(SortedData, kclass[i - 1] + 1, kclass[i]);
                Breaks[i] = Clusters[i].MaxValue;
            }
        }

    }
}
