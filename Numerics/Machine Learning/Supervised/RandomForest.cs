using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;
using System;
using System.Threading.Tasks;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// The Random Forest method for regression and classification.
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
    /// <b> References:</b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Random_forest" />
    /// </para>
    /// </remarks>
    public class RandomForest
    {

        #region Construction

        /// <summary>
        /// Create new Random Forest.
        /// </summary>
        /// <param name="x">The training 1D array of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public RandomForest(double[] x, double[] y, int seed = -1)
        {
            // Set inputs
            Y = new Vector(y);
            X = new Matrix(x);
            Dimensions = X.NumberOfColumns;
            Features = Math.Max(1, Dimensions - 1);
            Random = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");

        }

        /// <summary>
        /// Create new Random Forest.
        /// </summary>
        /// <param name="x">The training 2D array of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public RandomForest(double[,] x, double[] y, int seed = -1)
        {
            // Set inputs
            Y = new Vector(y);
            X = new Matrix(x);
            Dimensions = X.NumberOfColumns;
            Features = Math.Max(1, Dimensions - 1);
            Random = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");

        }

        /// <summary>
        /// Create new Random Forest.
        /// </summary>
        /// <param name="x">The training matrix of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public RandomForest(Matrix x, Vector y, int seed = -1)
        {
            // Set inputs
            Y = y;
            X = x;
            Dimensions = X.NumberOfColumns;
            Features = Math.Max(1, Dimensions - 1);
            Random = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");

        }

        #endregion

        #region Members

        /// <summary>
        /// The number of trees to use in the Random Forest. Default = 1000.
        /// </summary>
        public int NumberOfTrees { get; set; } = 1000;

        /// <summary>
        /// The minimum split size of the samples. Default = 2.
        /// </summary>
        public int MinimumSplitSize { get; set; } = 2;

        /// <summary>
        /// The maximum recursion depth. Default = 100.
        /// </summary>
        public int MaxDepth { get; set; } = 100;

        /// <summary>
        /// The dimensionality (or total number of features) of the data space.
        /// </summary>
        public int Dimensions { get; private set; }

        /// <summary>
        /// The number of random sub features to evaluate in the tree recursion.
        /// </summary>
        public int Features { get; set; }

        /// <summary>
        /// The random number generator to be used within the decision tree estimation.
        /// </summary>
        public Random Random { get; set; }

        /// <summary>
        /// The training vector of response values.
        /// </summary>
        public Vector Y { get; private set; }

        /// <summary>
        /// The training matrix of predictor values. 
        /// </summary>
        public Matrix X { get; private set; }

        /// <summary>
        /// The array of decision trees.
        /// </summary>
        public DecisionTree[] DecisionTrees { get; private set; } = null!;

        /// <summary>
        /// Determines whether this is for regression or classification. Default = regression.
        /// </summary>
        public bool IsRegression { get; set; } = true;

        /// <summary>
        /// Determines if the Random Forest has been trained.
        /// </summary>
        public bool IsTrained { get; private set; } = false;

        #endregion

        #region Methods

        /// <summary>
        /// Train the Random Forest.
        /// </summary>
        public void Train()
        {
            IsTrained = false;
            Features = Math.Min(Features, Dimensions);
            DecisionTrees = new DecisionTree[NumberOfTrees];
            var seeds = Random.NextIntegers(NumberOfTrees);

            // Estimate trees in parallel
            Parallel.For(0, NumberOfTrees, idx =>
            {
                DecisionTrees[idx] = BootstrapDecisionTree(seeds[idx]);
                DecisionTrees[idx].Train();
            });

            IsTrained = true;
        }

        /// <summary>
        /// Returns a bootstrapped decision tree.
        /// </summary>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        private DecisionTree BootstrapDecisionTree(int seed = -1)
        {
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var idxs = rnd.NextIntegers(0, X.NumberOfRows, X.NumberOfRows);
            var bootX = new Matrix(X.NumberOfRows, X.NumberOfColumns);
            var bootY = new Vector(Y.Length);
            for (int i = 0; i < X.NumberOfRows; i++)
            {
                for (int j = 0; j < X.NumberOfColumns; j++)
                {
                    bootX[i, j] = X[idxs[i], j];
                }
                bootY[i] = Y[idxs[i]];
            }
            return new DecisionTree(bootX, bootY, seed) { MinimumSplitSize = MinimumSplitSize, MaxDepth = MaxDepth, Features = Features, IsRegression = IsRegression };        
        }

        /// <summary>
        /// Returns the prediction intervals in a 2D array with columns: lower, median, upper, mean. 
        /// </summary>
        /// <param name="X">The 1D array of predictors.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,]? Predict(double[] X, double alpha = 0.1)
        {
            return Predict(new Matrix(X), alpha);
        }

        /// <summary>
        /// Returns the prediction intervals in a 2D array with columns: lower, median, upper, mean. 
        /// </summary>
        /// <param name="X">The 2D array of predictors.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,]? Predict(double[,] X, double alpha = 0.1)
        {
            return Predict(new Matrix(X), alpha);
        }

        /// <summary>
        /// Returns the prediction intervals in a 2D array with columns: lower, median, upper, mean. 
        /// </summary>
        /// <param name="X">The matrix of predictors.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,]? Predict(Matrix X, double alpha = 0.1)
        {
            if (!IsTrained) return null!;

            var percentiles = new double[] { alpha / 2d, 0.5, 1d - alpha / 2d };
            var output = new double[X.NumberOfRows, 4]; // lower, median, upper, mean

            var bootResults = new double[X.NumberOfRows, NumberOfTrees];

            // Bootstrap the predictions
            Parallel.For(0, NumberOfTrees, idx => { bootResults.SetColumn(idx, DecisionTrees[idx].Predict(X)!); });

            // Process results
            Parallel.For(0, X.NumberOfRows, idx =>
            {
                var values = bootResults.GetRow(idx);
                Array.Sort(values);

                if (IsRegression)
                {
                    // Record percentiles for CIs
                    for (int j = 0; j < percentiles.Length; j++)
                        output[idx, j] = Statistics.Percentile(values, percentiles[j], true);
                    output[idx, 3] = Statistics.ParallelMean(values);
                }
                else
                {
                    // Record percentiles for CIs
                    for (int j = 0; j < percentiles.Length; j++)
                        output[idx, j] = Math.Floor(Statistics.Percentile(values, percentiles[j], true));
                    output[idx, 3] = Math.Floor(Statistics.ParallelMean(values));
                }

            });

            return output;
        }



        #endregion

    }
}
