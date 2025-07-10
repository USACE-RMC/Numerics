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
using Numerics.Sampling;
using System;
using System.Linq;
using System.Threading.Tasks;

namespace Numerics.MachineLearning
{

    /// <summary>
    /// k-Means clustering.
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
    /// k-Means clustering is a method of vector quantization, originally from signal processing, 
    /// that aims to partition n observations into k clusters in which each observation belongs to 
    /// the cluster with the nearest mean (cluster centers or cluster centroid), serving as a prototype of the cluster.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/K-means_clustering" />
    /// </para>
    /// </remarks>
    public class KMeans
    {
        /// <summary>
        /// Creates a new k-Means clustering analysis.
        /// </summary>
        /// <param name="X">The 1D array of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        public KMeans(float[] X, int k)
        {
            this.K = k;
            this.X = new Matrix(X);
            Dimension = this.X.NumberOfColumns;
            Means = new double[K, Dimension];
            Labels = new int[this.X.NumberOfRows];
        }

        /// <summary>
        /// Creates a new k-Means clustering analysis.
        /// </summary>
        /// <param name="X">The 1D array of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        public KMeans(double[] X, int k) 
        { 
            this.K = k;
            this.X = new Matrix(X);
            Dimension = this.X.NumberOfColumns;
            Means = new double[K, Dimension];
            Labels = new int[this.X.NumberOfRows];
        }

        /// <summary>
        /// Creates a new k-Means clustering analysis.
        /// </summary>
        /// <param name="X">The 2D array of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        public KMeans(double[,] X, int k)
        {
            this.K = k;
            this.X = new Matrix(X);
            Dimension = this.X.NumberOfColumns;
            Means = new double[K, Dimension];
            Labels = new int[this.X.NumberOfRows];
        }

        /// <summary>
        /// Creates a new k-Means clustering analysis.
        /// </summary>
        /// <param name="X">The matrix of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        public KMeans(Matrix X, int k)
        {
            this.K = k;
            this.X = X;
            Dimension = this.X.NumberOfColumns;
            Means = new double[K, Dimension];
            Labels = new int[this.X.NumberOfRows];
        }

        /// <summary>
        /// The number of clusters.
        /// </summary>
        public int K { get; private set; }

        /// <summary>
        /// The matrix of predictor values. 
        /// </summary>
        public Matrix X { get; private set; }

        /// <summary>
        /// The dimensionality (or number of features) of the data space.
        /// </summary>
        public int Dimension { get; private set; }

        /// <summary>
        /// The cluster means.
        /// </summary>
        public double[,] Means { get; private set; }

        /// <summary>
        /// The array of cluster labels assigned to each of the data points.
        /// </summary>
        public int[] Labels { get; private set; }

        /// <summary>
        /// The maximum iterations in the clustering algorithm. Default = 1,000. 
        /// </summary>
        public int MaxIterations { get; set; } = 1000;

        /// <summary>
        /// The total number of iterations required to find the clusters.
        /// </summary>
        public int Iterations { get; private set; }

        /// <summary>
        /// Estimate the k-Means clusters.
        /// </summary>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="kMeansPlusPlus">Determines whether to use random initialization or to use the k-Means++ method. Default is to use k-Means++.</param>
        public void Train(int seed = -1, bool kMeansPlusPlus= true)
        {

            // 1. Initialize cluster centers 
            Means = Initialize(X, K, seed, kMeansPlusPlus);

            // 2. Optimize clusters
            for (Iterations = 1; Iterations <= MaxIterations; Iterations++)
            {
                // Perform E-step
                // Assign samples to closest centroids (create clusters)
                var oldLabels = Labels.ToArray();
                Labels = GetLabels(Means);

                // Check if the labels changed
                bool labelsChanged = false;
                for (int i = 0;  i < Labels.Length; i++)
                {
                    if (oldLabels[i] != Labels[i])
                    {
                        labelsChanged = true;
                        break;
                    }

                }
                // Stop when the E-step doesn't change the assignment of any data point
                if (labelsChanged == false)
                    break;

                // Perform M-step
                // Calculate new centroids from the clusters
                Means = GetCentroids(Labels);
 
            }

        }

        /// <summary>
        /// Initializes the centroids of the k-Means clusters.
        /// </summary>
        /// <param name="X">The matrix of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="kMeansPlusPlus">Determines whether to use random initialization or to use the k-Means++ method. Default is to use k-Means++.</param>
        public static double[,] Initialize(Matrix X, int k, int seed = -1, bool kMeansPlusPlus = true)
        {
            var centroids = new double[k, X.NumberOfColumns];
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            if (kMeansPlusPlus == false)
            {
                
                var rndIdxs = rnd.NextIntegers(0, X.NumberOfRows, k, false);
                Array.Sort(rndIdxs);
                for (int i = 0; i < k; i++)
                    for (int j = 0; j < X.NumberOfColumns; j++)
                        centroids[i, j] = X[rndIdxs[i], j];
            }
            else
            {
                // Initialize using K-Means++
                // http://en.wikipedia.org/wiki/K-means%2B%2B


                // 1. Choose one center uniformly at random from among the data points.
                int idx = rnd.Next(0, X.NumberOfRows);
                for (int j = 0; j < X.NumberOfColumns; j++)
                    centroids[0, j] = X[idx, j];

                for (int c = 1; c < k; c++)
                {
                    // 2. For each data point x not chosen yet, compute D(x),
                    // the distance between x and the nearest center that has already been chosen.

                    double sum = 0;
                    var D = new double[X.NumberOfRows];
                    for (int i = 0; i < X.NumberOfRows; i++)
                    {
                        var x = X.Row(i);

                        double min = Tools.Distance(x, centroids.GetRow(0));
                        for (int j = 1; j < c; j++)
                        {
                            double d = Tools.Distance(x, centroids.GetRow(0));

                            if (d < min)
                                min = d;
                        }

                        D[i] = min;
                        sum += min;
                    }

                    // Following Accord.Net checks:
                    // https://github.com/accord-net/framework/blob/development/Sources/Accord.MachineLearning/Clustering/KMeans/KMeans.cs
                    
                    // Note: the following checks could have been avoided if we added
                    // a small value to each distance, but is kept as this to avoid 022
                    // breaking the random pattern in existing code.

                    if (sum == 0)
                    {
                        // Degenerate case: all points are the same, chose any of them
                        idx = rnd.Next(0, X.NumberOfRows);
                    }
                    else
                    {
                        // 3. Choose one new data point at random as a new center, using a weighted
                        //    probability distribution where a point x is chosen with probability 
                        //    proportional to D(x)^2.
                        var u = rnd.NextDouble();
                        var cdf = new double[X.NumberOfRows];
                        for (int i = 0; i < D.Length; i++)
                        {
                            D[i] /= sum;
                            cdf[i] = i == 0 ? D[i] : cdf[i - 1] + D[i];
                            if (u <= cdf[i])
                            {
                                idx = i;
                                break;
                            }
                        }
                            
                    }
                    for (int j = 0; j < X.NumberOfColumns; j++) 
                        centroids[c, j] = X[idx, j];
                }
            }

            return centroids;
        }

        /// <summary>
        /// Gets the array of cluster labels given the list of centroids.
        /// </summary>
        /// <param name="centroids">The list of centroids.</param>
        private int[] GetLabels(double[,] centroids)
        {
            // Assign samples to the closest centroids
            var labels = new int[X.NumberOfRows];       
            Parallel.For(0, X.NumberOfRows, idx =>  {  labels[idx] = GetClosestCentroid(X.Row(idx), centroids); });
            return labels;
        }

        /// <summary>
        /// Gets the centroids given the specified cluster labels.
        /// </summary>
        /// <param name="labels">The array of cluster labels assigned to each of the data points.</param>
        private double[,] GetCentroids(int[] labels)
        {
            var count = new double[K];
            var centroids = new double[K, Dimension];

            // Get sums and counts
            for (int i = 0; i < X.NumberOfRows; i++)
            {
                count[labels[i]]++;
                for (int j = 0; j < Dimension; j++)
                    centroids[labels[i], j] += X[i, j];
            };

            // Get mean of clusters
            for (int k = 0; k < K; k++)
                for (int j = 0; j < Dimension; j++)
                    centroids[k, j] /= count[k] > 0 ? count[k] : 1;
            
            return centroids;
        }

        /// <summary>
        /// Returns the index of the centroid closest to the sample vector.
        /// </summary>
        /// <param name="sample">The sample vector.</param>
        /// <param name="centroids">The list of centroids.</param>
        private int GetClosestCentroid(double[] sample, double[,] centroids)
        {
            double min = double.MaxValue;
            int minIdx = 0;
            for (int k = 0; k < K; k++)
            {
                var dist = Tools.Distance(sample, centroids.GetRow(k));
                if (dist < min)
                {
                    min = dist;
                    minIdx = k;
                }
            }
            return minIdx;
        }

    }

}
