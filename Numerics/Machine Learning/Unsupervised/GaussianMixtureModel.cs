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
        /// The column variances of <see cref="X"/>, computed once per fit for the covariance floors.
        /// </summary>
        private double[]? _columnVariances;

        /// <summary>
        /// Estimate the Gaussian Mixture Model.
        /// </summary>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="kMeansPlusPlus">Determines whether to use random initialization or to use the k-Means++ method. Default is to use k-Means++.</param>
        /// <exception cref="InvalidOperationException">Thrown when a component covariance cannot be factorized.</exception>
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
            var colVar = ColumnVariances();
            _columnVariances = colVar;
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
                    Sigmas[k][d, d] = Math.Max(Sigmas[k][d, d], 1E-6 * colVar[d]);
            }

            // 2. Optimize clusters
            double oldLogLH = double.MinValue, newLogLH = double.MinValue;
            for (Iterations = 1; Iterations <= MaxIterations; Iterations++)
            {
                // Perform the expectation step. The log-likelihood records every iteration, so a
                // run that exhausts its iteration budget reports its final evaluated value.
                newLogLH = EStep();
                LogLikelihood = newLogLH;

                // Check convergence
                if (Math.Abs((oldLogLH - newLogLH) / oldLogLH) < Tolerance)
                {
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
        /// <remarks>
        /// Each row's component likelihoods, label and normalization depend only on that row, so rows
        /// are processed in parallel for large problems; the log-likelihood is then summed in row
        /// order. Each row is whitened by <see cref="CholeskyDecomposition.Forward(Vector)"/> with the
        /// Cholesky factor of its component covariance; the factorizations are only read, so the rows
        /// share them.
        /// </remarks>
        private double EStep()
        {
            int n = X.NumberOfRows;
            var logDet = new double[K];
            var factors = new CholeskyDecomposition[K];
            var logWeights = new double[K];

            // Decompose each covariance once
            for (int k = 0; k < K; k++)
            {
                CholeskyDecomposition cholesky;
                try
                {
                    cholesky = new CholeskyDecomposition(Sigmas[k]);
                }
                catch (Exception exception)
                {
                    throw new InvalidOperationException(
                        $"Gaussian mixture component {k + 1} covariance could not be factorized.",
                        exception);
                }
                logDet[k] = cholesky.LogDeterminant();
                factors[k] = cholesky;
                logWeights[k] = Math.Log(Weights[k]);
            }

            var rowLogLH = new double[n];
            RangeParallel.For(n, (long)n * K * Dimension * (Dimension + 2), i =>
            {
                var u = new double[Dimension];
                var residual = new Vector(u);
                for (int k = 0; k < K; k++)
                {
                    for (int d = 0; d < Dimension; d++)
                        u[d] = X[i, d] - Means[k, d];
                    var v = factors[k].Forward(residual);
                    double sum = 0;
                    for (int d = 0; d < Dimension; d++)
                        sum += Tools.Sqr(v[d]);
                    LikelihoodMatrix[i, k] = -0.5 * (sum + logDet[k]) + logWeights[k];
                }

                // The row holds unnormalized log-likelihoods; normalize with log-sum-exp and record the most likely label.
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
                double total = 0;
                for (int k = 0; k < K; k++)
                    total += Math.Exp(LikelihoodMatrix[i, k] - max);
                double tmp = max + Math.Log(total);
                for (int k = 0; k < K; k++)
                    LikelihoodMatrix[i, k] = Math.Exp(LikelihoodMatrix[i, k] - tmp);
                rowLogLH[i] = tmp;
            });

            double logLH = 0;
            for (int i = 0; i < n; i++)
                logLH += rowLogLH[i];
            // The likelihood matrix omits the constant -D/2 * ln(2*pi) per point because it cancels
            // in the responsibilities; the reported log-likelihood restores it.
            return logLH - 0.5 * Dimension * X.NumberOfRows * Math.Log(2.0 * Math.PI);
        }

        /// <summary>
        /// Returns the population variance of each column of <see cref="X"/>, used to floor the
        /// component covariance diagonals.
        /// </summary>
        /// <returns>The column variances.</returns>
        /// <remarks>
        /// The variances depend only on the data, so they are computed once and shared by every
        /// component.
        /// These are population variances (divided by n), which differ from the sample variances (divided by n - 1) of <see cref="Numerics.Data.Statistics.Statistics.ColumnVariances(double[,])"/>.
        /// </remarks>
        private double[] ColumnVariances()
        {
            var result = new double[Dimension];
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
                result[d] = colVar;
            }
            return result;
        }

        /// <summary>
        /// The maximization step.
        /// </summary>
        /// <remarks>
        /// The means are updated first, and every covariance entry (d, j) is then formed about the updated
        /// means d and j, as in the EM update of Bishop (2006), Section 9.2.2. Each (component, dimension)
        /// task sums over rows in row order, so the tasks run in parallel for large problems.
        /// </remarks>
        private void MStep()
        {
            var colVar = _columnVariances ?? ColumnVariances();
            int n = X.NumberOfRows;
            var wgts = new double[K];
            for (int k = 0; k < K; k++)
            {
                double wgt = 0d;
                for (int i = 0; i < n; i++)
                    wgt += LikelihoodMatrix[i, k];
                wgts[k] = wgt;
                Weights[k] = wgt / X.NumberOfRows;
            }

            // Compute centroids
            RangeParallel.For(K * Dimension, (long)n * K * Dimension, t =>
            {
                int k = t / Dimension, d = t % Dimension;
                double wgt = wgts[k];
                // A component with no responsibility keeps its previous parameters; updating it
                // would divide by zero and spread NaN through the model
                if (wgt <= 0)
                    return;
                double sum = 0;
                for (int i = 0; i < n; i++)
                    sum += LikelihoodMatrix[i, k] * X[i, d];
                Means[k, d] = sum / wgt;
            });

            // Compute covariance
            RangeParallel.For(K * Dimension, (long)n * K * Dimension * Dimension, t =>
            {
                int k = t / Dimension, d = t % Dimension;
                double wgt = wgts[k];
                if (wgt <= 0)
                    return;
                double meanD = Means[k, d];
                for (int j = 0; j < Dimension; j++)
                {
                    double meanJ = Means[k, j];
                    double sum = 0;
                    for (int i = 0; i < n; i++)
                    {
                        sum += LikelihoodMatrix[i, k] * (X[i, d] - meanD) * (X[i, j] - meanJ);
                    }
                    Sigmas[k][d, j] = sum / wgt;
                }
            });

            for (int k = 0; k < K; k++)
            {
                if (wgts[k] <= 0)
                    continue;

                // Floor diagonal at a fraction of the overall data variance to prevent
                // component collapse. When a component captures very few points, its
                // covariance can become singular, causing Cholesky decomposition in the
                // E-step to fail. This mirrors sklearn's reg_covar parameter.
                for (int d = 0; d < Dimension; d++)
                    Sigmas[k][d, d] = Math.Max(Sigmas[k][d, d], 1E-6 * colVar[d]);

                // Ensure the full covariance matrix remains symmetric positive-definite. The helper is
                // pure: it returns the symmetrized covariance unchanged when usable and adds a
                // trace-scaled ridge only when needed. Its result must be assigned so an actual repair
                // reaches the covariance consumed by the next E-step's Cholesky factorization.
                Sigmas[k] = MatrixRegularization.MakeSymmetricPositiveDefinite(Sigmas[k]);
            }
        }

    }
}
