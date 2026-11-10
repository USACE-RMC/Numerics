using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;
using System;
using System.Linq;

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
            if (k < 1 || k > this.X.NumberOfRows)
                throw new ArgumentOutOfRangeException(nameof(k), "The number of clusters must be between 1 and the number of data rows.");
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
            if (k < 1 || k > this.X.NumberOfRows)
                throw new ArgumentOutOfRangeException(nameof(k), "The number of clusters must be between 1 and the number of data rows.");
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
            if (k < 1 || k > this.X.NumberOfRows)
                throw new ArgumentOutOfRangeException(nameof(k), "The number of clusters must be between 1 and the number of data rows.");
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
            if (k < 1 || k > this.X.NumberOfRows)
                throw new ArgumentOutOfRangeException(nameof(k), "The number of clusters must be between 1 and the number of data rows.");
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
        /// <remarks>
        /// The k-Means++ seeding chooses each new center with probability proportional to the distance
        /// D(x) to the nearest chosen center, a variant of Arthur and Vassilvitskii (2007), who use D(x)^2.
        /// </remarks>
        public void Train(int seed = -1, bool kMeansPlusPlus= true)
        {

            // 1. Initialize cluster centers
            var rows = new double[X.NumberOfRows][];
            for (int i = 0; i < rows.Length; i++)
                rows[i] = X.Row(i);
            Means = Initialize(rows, Dimension, K, seed, kMeansPlusPlus);

            // 2. Optimize clusters
            for (Iterations = 1; Iterations <= MaxIterations; Iterations++)
            {
                // Perform E-step
                // Assign samples to closest centroids (create clusters)
                var oldLabels = Labels;
                Labels = GetLabels(rows, Means);

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
                {
                    // With a single cluster, and for any fit whose first assignment matches the
                    // zero-initialized labels, the loop exits before a centroid has ever been
                    // computed, leaving the initializer's seed point as the reported mean. Run the
                    // M-step once so the reported means are the centroids of the final labels; on a
                    // later iteration convergence means the labels did not change, so the centroids
                    // already reflect them and the reported means are untouched.
                    if (Iterations == 1)
                        Means = GetCentroids(rows, Labels);
                    break;
                }

                // Perform M-step
                // Calculate new centroids from the clusters
                Means = GetCentroids(rows, Labels);

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
            var rows = new double[X.NumberOfRows][];
            for (int i = 0; i < rows.Length; i++)
                rows[i] = X.Row(i);
            return Initialize(rows, X.NumberOfColumns, k, seed, kMeansPlusPlus);
        }

        /// <summary>
        /// Initializes the centroids of the k-Means clusters from the data rows.
        /// </summary>
        /// <param name="rows">The rows of predictor values.</param>
        /// <param name="dimension">The number of columns.</param>
        /// <param name="k">The number of clusters.</param>
        /// <param name="seed">The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="kMeansPlusPlus">Determines whether to use random initialization or the k-Means++ method.</param>
        /// <returns>The initial centroids.</returns>
        /// <remarks>
        /// For k-Means++, each point keeps its distance to the nearest centroid chosen so far and
        /// compares it only with the newest centroid. These are the comparisons, in order, of a scan
        /// over every chosen centroid, so the seeding costs O(n k) and the draws match that scan.
        /// </remarks>
        private static double[,] Initialize(double[][] rows, int dimension, int k, int seed, bool kMeansPlusPlus)
        {
            int n = rows.Length;
            var centroids = new double[k, dimension];
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            if (kMeansPlusPlus == false)
            {

                var rndIdxs = rnd.NextIntegers(0, n, k, false);
                Array.Sort(rndIdxs);
                for (int i = 0; i < k; i++)
                    for (int j = 0; j < dimension; j++)
                        centroids[i, j] = rows[rndIdxs[i]][j];
            }
            else
            {
                // Initialize using K-Means++
                // http://en.wikipedia.org/wiki/K-means%2B%2B


                // 1. Choose one center uniformly at random from among the data points.
                int idx = rnd.Next(0, n);
                for (int j = 0; j < dimension; j++)
                    centroids[0, j] = rows[idx][j];

                var minDistance = new double[n];
                var D = new double[n];
                var cdf = new double[n];
                for (int c = 1; c < k; c++)
                {
                    // 2. For each data point x not chosen yet, compute D(x),
                    // the distance between x and the nearest center that has already been chosen.
                    int newest = c - 1;
                    var center = new double[dimension];
                    for (int j = 0; j < dimension; j++)
                        center[j] = centroids[newest, j];
                    RangeParallel.For(n, (long)n * dimension, i =>
                    {
                        double d = Tools.Distance(rows[i], center);
                        if (newest == 0 || d < minDistance[i])
                            minDistance[i] = d;
                    });

                    double sum = 0;
                    for (int i = 0; i < n; i++)
                    {
                        D[i] = minDistance[i];
                        sum += D[i];
                    }

                    // Following Accord.Net checks:
                    // https://github.com/accord-net/framework/blob/development/Sources/Accord.MachineLearning/Clustering/KMeans/KMeans.cs

                    // Note: adding a small value to each distance would avoid the
                    // following checks, but it would change the random sequence
                    // drawn for a given seed.

                    if (sum == 0)
                    {
                        // Degenerate case: all points are the same, chose any of them
                        idx = rnd.Next(0, n);
                    }
                    else
                    {
                        // 3. Choose one new data point at random as a new center, using a weighted
                        //    probability distribution where a point x is chosen with probability
                        //    proportional to D(x). Arthur and Vassilvitskii (2007) sample proportional
                        //    to D(x)^2; this variant samples proportional to the distance itself.
                        var u = rnd.NextDouble();
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
                    for (int j = 0; j < dimension; j++)
                        centroids[c, j] = rows[idx][j];
                }
            }

            return centroids;
        }

        /// <summary>
        /// Gets the array of cluster labels given the list of centroids.
        /// </summary>
        /// <param name="rows">The rows of predictor values.</param>
        /// <param name="centroids">The list of centroids.</param>
        private int[] GetLabels(double[][] rows, double[,] centroids)
        {
            // Assign samples to the closest centroids
            var labels = new int[X.NumberOfRows];
            var c = new double[K][];
            for (int k = 0; k < K; k++)
            {
                c[k] = new double[Dimension];
                for (int j = 0; j < Dimension; j++)
                    c[k][j] = centroids[k, j];
            }
            RangeParallel.For(X.NumberOfRows, (long)X.NumberOfRows * K * Dimension, idx => { labels[idx] = GetClosestCentroid(rows[idx], c); });
            return labels;
        }

        /// <summary>
        /// Gets the centroids given the specified cluster labels.
        /// </summary>
        /// <param name="rows">The rows of predictor values.</param>
        /// <param name="labels">The array of cluster labels assigned to each of the data points.</param>
        private double[,] GetCentroids(double[][] rows, int[] labels)
        {
            var count = new double[K];
            var centroids = new double[K, Dimension];

            // Get sums and counts
            for (int i = 0; i < X.NumberOfRows; i++)
            {
                count[labels[i]]++;
                for (int j = 0; j < Dimension; j++)
                    centroids[labels[i], j] += rows[i][j];
            };

            // Get mean of clusters
            for (int k = 0; k < K; k++)
                for (int j = 0; j < Dimension; j++)
                    centroids[k, j] /= count[k] > 0 ? count[k] : 1;

            // Relocate each empty cluster to the point farthest from its assigned centroid, one
            // distinct point per empty cluster, matching the scikit-learn treatment. An empty
            // cluster otherwise sits at the origin, which is not a property of the data.
            var consumed = new bool[X.NumberOfRows];
            for (int k = 0; k < K; k++)
            {
                if (count[k] > 0)
                    continue;
                int farthest = -1;
                double maxDistance = -1;
                for (int i = 0; i < X.NumberOfRows; i++)
                {
                    if (consumed[i])
                        continue;
                    double distance = 0;
                    for (int j = 0; j < Dimension; j++)
                        distance += Tools.Sqr(rows[i][j] - centroids[labels[i], j]);
                    if (distance > maxDistance)
                    {
                        maxDistance = distance;
                        farthest = i;
                    }
                }
                if (farthest < 0)
                    continue;
                consumed[farthest] = true;
                for (int j = 0; j < Dimension; j++)
                    centroids[k, j] = rows[farthest][j];
            }

            return centroids;
        }

        /// <summary>
        /// Returns the index of the centroid closest to the sample vector.
        /// </summary>
        /// <param name="row">The sample vector.</param>
        /// <param name="centroids">The centroids.</param>
        private int GetClosestCentroid(double[] row, double[][] centroids)
        {
            double min = double.MaxValue;
            int minIdx = 0;
            for (int k = 0; k < K; k++)
            {
                var dist = Tools.Distance(row, centroids[k]);
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
