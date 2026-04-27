using Numerics.Mathematics.LinearAlgebra;
using System;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// Gaussian Mixture Model (GMM) for classification.
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
    /// A Gaussian mixture model is a probabilistic model that assumes all the data points are 
    /// generated from a mixture of a finite number of Gaussian distributions with unknown parameters. 
    /// One can think of mixture models as generalizing k-means clustering to incorporate information 
    /// about the covariance structure of the data as well as the centers of the latent Gaussian distributions.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/EM_algorithm_and_GMM_model" />
    /// </para>
    /// </remarks>
    public class GaussianMixtureModel
    {
        /// <summary>
        /// Creates a new Gaussian Mixture Model.
        /// </summary>
        /// <param name="X">The 1D array of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        public GaussianMixtureModel(float[] X, int k)
        {
            this.K = k;
            this.X = new Matrix(X);
            Dimension = this.X.NumberOfColumns;
            Means = new double[K, Dimension];
            Sigmas = new Matrix[K];
            Labels = new int[this.X.NumberOfRows];
        }

        /// <summary>
        /// Creates a new Gaussian Mixture Model.
        /// </summary>
        /// <param name="X">The 1D array of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        public GaussianMixtureModel(double[] X, int k)
        {
            this.K = k;
            this.X = new Matrix(X);
            Dimension = this.X.NumberOfColumns;
            Means = new double[K, Dimension];
            Sigmas = new Matrix[K];
            Labels = new int[this.X.NumberOfRows];
        }

        /// <summary>
        /// Creates a new Gaussian Mixture Model.
        /// </summary>
        /// <param name="X">The 2D array of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        public GaussianMixtureModel(double[,] X, int k)
        {
            this.K = k;
            this.X = new Matrix(X);
            Dimension = this.X.NumberOfColumns;
            Means = new double[K, Dimension];
            Sigmas = new Matrix[K];
            Labels = new int[this.X.NumberOfRows];
        }

        /// <summary>
        /// Creates a new Gaussian Mixture Model.
        /// </summary>
        /// <param name="X">The matrix of predictor values.</param>
        /// <param name="k">The number of clusters.</param>
        public GaussianMixtureModel(Matrix X, int k)
        {
            this.K = k;
            this.X = X;
            Dimension = this.X.NumberOfColumns;
            Means = new double[K, Dimension];
            Sigmas = new Matrix[K];
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
        /// The cluster covariance matrices.
        /// </summary>
        public Matrix[] Sigmas { get; private set; }

        /// <summary>
        /// The array of cluster labels assigned to each of the data points.
        /// </summary>
        public int[] Labels { get; private set; }

        /// <summary>
        /// The mixing weights 
        /// </summary>
        public double[] Weights { get; private set; } = null!;

        /// <summary>
        /// The likelihood of each data point (row) and for each cluster (column).
        /// </summary>
        public double[,] LikelihoodMatrix { get; private set; } = new double[0, 0];

        /// <summary>
        /// The total log-likelihood of the fit.
        /// </summary>
        public double LogLikelihood { get; private set; }

        /// <summary>
        /// The maximum iterations in the clustering algorithm. Default = 1,000. 
        /// </summary>
        public int MaxIterations { get; set; } = 1000;

        /// <summary>
        /// The relative tolerance for convergence. Default = 1E-8.
        /// </summary>
        public double Tolerance { get; set; } = 1E-8;

        /// <summary>
        /// The total number of iterations required to find the clusters.
        /// </summary>
        public int Iterations { get; private set; }

        /// <summary>
        /// Estimate the Gaussian Mixture Model.
        /// </summary>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="kMeansPlusPlus">Determines whether to use random initialization or to use the k-Means++ method. Default is to use k-Means++.</param>
        public void Train(int seed = -1, bool kMeansPlusPlus = true)
        {
            // 1. Initialize clusters from k-Means
            var kMeans = new KMeans(X, K);
            kMeans.Train(seed, kMeansPlusPlus);
            Means = kMeans.Means;

            // Give equal weight to each cluster
            // and initialize the covariance matrix from cluster data variance.
            // Using actual data variance (instead of a tiny constant like 1e-10)
            // prevents the first E-step from computing likelihoods from near-degenerate
            // Gaussians, which can cause numerical overflow or underflow.
            Weights = new double[K];
            LikelihoodMatrix = new double[X.NumberOfRows, K];
            Sigmas = new Matrix[K];
            for (int k = 0; k < K; k++)
            {
                Weights[k] = 1d / K;
                Sigmas[k] = new Matrix(Dimension);

                // Compute within-cluster covariance from k-means labels
                int clusterCount = 0;
                for (int i = 0; i < X.NumberOfRows; i++)
                {
                    if (kMeans.Labels[i] == k)
                        clusterCount++;
                }

                if (clusterCount > 1)
                {
                    // Compute covariance from cluster members
                    for (int d = 0; d < Dimension; d++)
                    {
                        for (int j = 0; j < Dimension; j++)
                        {
                            double sum = 0;
                            for (int i = 0; i < X.NumberOfRows; i++)
                            {
                                if (kMeans.Labels[i] == k)
                                    sum += (X[i, d] - Means[k, d]) * (X[i, j] - Means[k, j]);
                            }
                            Sigmas[k][d, j] = sum / clusterCount;
                        }
                    }
                }

                // Ensure positive-definite: floor diagonal at fraction of overall data variance
                for (int d = 0; d < Dimension; d++)
                {
                    double colVar = 0;
                    double colMean = 0;
                    for (int i = 0; i < X.NumberOfRows; i++)
                        colMean += X[i, d];
                    colMean /= X.NumberOfRows;
                    for (int i = 0; i < X.NumberOfRows; i++)
                        colVar += (X[i, d] - colMean) * (X[i, d] - colMean);
                    colVar /= X.NumberOfRows;

                    Sigmas[k][d, d] = Math.Max(Sigmas[k][d, d], 1E-6 * colVar);
                }
            }

            // 2. Optimize clusters
            double oldLogLH = double.MinValue, newLogLH = double.MinValue;
            for (Iterations = 1; Iterations <= MaxIterations; Iterations++)
            {
                // Perform the expectation step
                newLogLH = EStep();

                // Check convergence
                if (Math.Abs((oldLogLH - newLogLH) / oldLogLH) < Tolerance)
                {
                    LogLikelihood = newLogLH;
                    break;
                }

                // Perform the maximization step
                MStep();

                // Update log-likelihood state
                oldLogLH = newLogLH;

            }

        }


        /// <summary>
        /// The expectation step. 
        /// </summary>
        /// <returns>Returns the log-likelihood.</returns>
        private double EStep()
        {
            var logDet = new double[K];

            // Outer loop for computing the likelihoods
            for (int k = 0; k < K; k++)
            {
                // Decompose the covariance in the outer loop
                var cholesky = new CholeskyDecomposition(Sigmas[k]);
                logDet[k] = cholesky.LogDeterminant();
                for (int i = 0; i < X.NumberOfRows; i++)
                {
                    // Inner loop for likelihoods
                    var u = new Vector(Dimension);
                    for (int d = 0; d < Dimension; d++)
                        u[d] = X[i, d] - Means[k, d];
                    // Solve L*v=u
                    var v = cholesky.Forward(u);
                    double sum = 0;
                    for (int d = 0; d < Dimension; d++)
                        sum += Tools.Sqr(v[d]);
                    LikelihoodMatrix[i, k] = -0.5 * (sum + logDet[k]) + Math.Log(Weights[k]);
                }
            }
            // At this point we have unnormalized logs of the likelihoods.
            // We need to normalize using log-sum-exp and compute the log-likelihoods.
            double logLH = 0;
            for (int i = 0; i < X.NumberOfRows; i++)
            {
                // Get max likelihood and index
                double max = double.MinValue;
                int idx = -1;
                for (int k = 0; k < K; k++)
                {
                    if (LikelihoodMatrix[i, k] > max)
                    {
                        max = LikelihoodMatrix[i, k];
                        idx = k;
                    }
                }
                // update labels given the likelihoods
                Labels[i] = idx;

                // log-sum-exp trick begins here
                double sum = 0;
                for (int k = 0; k < K; k++)
                    sum += Math.Exp(LikelihoodMatrix[i, k] - max);
                double tmp = max + Math.Log(sum);
                for (int k = 0; k < K; k++)
                    LikelihoodMatrix[i, k] = Math.Exp(LikelihoodMatrix[i, k] - tmp);
                logLH += tmp;
            }
            return logLH;
        }
    
        /// <summary>
        /// The maximization step. 
        /// </summary>
        private void MStep()
        {
            for (int k = 0; k < K; k++)
            {
                double wgt = 0d;
                for (int i = 0; i < X.NumberOfRows; i++)
                    wgt += LikelihoodMatrix[i, k];
                Weights[k] = wgt / X.NumberOfRows;
                for (int d = 0; d < Dimension; d++)
                {
                    // Compute centroids
                    double sum = 0;
                    for (int i = 0; i < X.NumberOfRows; i++)
                        sum += LikelihoodMatrix[i, k] * X[i, d];
                    Means[k, d] = sum / wgt;
                    // Compute covariance
                    for (int j = 0; j < Dimension; j++)
                    {
                        sum = 0;
                        for (int i = 0; i < X.NumberOfRows; i++)
                        {
                            sum += LikelihoodMatrix[i, k] * (X[i, d] - Means[k, d]) * (X[i, j] - Means[k, j]);
                        }
                        Sigmas[k][d, j] = sum / wgt;
                    }
                }

                // Floor diagonal at a fraction of the overall data variance to prevent
                // component collapse. When a component captures very few points, its
                // covariance can become singular, causing Cholesky decomposition in the
                // E-step to fail. This mirrors sklearn's reg_covar parameter.
                for (int d = 0; d < Dimension; d++)
                {
                    double colVar = 0;
                    double colMean = 0;
                    for (int i = 0; i < X.NumberOfRows; i++)
                        colMean += X[i, d];
                    colMean /= X.NumberOfRows;
                    for (int i = 0; i < X.NumberOfRows; i++)
                        colVar += (X[i, d] - colMean) * (X[i, d] - colMean);
                    colVar /= X.NumberOfRows;

                    Sigmas[k][d, d] = Math.Max(Sigmas[k][d, d], 1E-6 * colVar);
                }

                // Ensure the full covariance matrix remains symmetric positive-definite
                MatrixRegularization.MakeSymmetricPositiveDefinite(Sigmas[k]);
            }
        }

    }
}
