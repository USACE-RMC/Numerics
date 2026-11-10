using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;
using System;
using System.Linq;
using System.Threading.Tasks;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// The k-Nearest Neighbors (k-NN) algorithm.
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
    /// k-NN is a simple method for applying supervised learning. When training a k-NN model,
    /// the model practically memorizes all locations of all the points and their values.
    /// When it is time to use the model for predictions, the model takes the input data and 
    /// calculates the k-nearest points to that input. 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm" />
    /// </para>
    /// </remarks>
    public class KNearestNeighbors
    {

        #region Construction

        /// <summary>
        /// Create new k-NN method.
        /// </summary>
        /// <param name="x">The training 1D array of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="k">The number of nearest neighbors.</param>
        public KNearestNeighbors(double[] x, double[] y, int k)
        {
            K = k;
            Y = new Vector(y);
            X = new Matrix(x);
            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");
        }

        /// <summary>
        /// Create new k-NN method.
        /// </summary>
        /// <param name="x">The training 2D array of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="k">The number of nearest neighbors.</param>
        public KNearestNeighbors(double[,] x, double[] y, int k)
        {
            K = k;
            Y = new Vector(y);
            X = new Matrix(x);
            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");
        }

        /// <summary>
        /// Create new k-NN method.
        /// </summary>
        /// <param name="x">The training matrix of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="k">The number of nearest neighbors.</param>
        public KNearestNeighbors(Matrix x, Vector y, int k)
        {
            K = k;
            Y = y;
            X = x;
            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");
        }

        #endregion

        #region Members

        /// <summary>
        /// The number of nearest neighbors.
        /// </summary>
        public int K { get; private set; }

        /// <summary>
        /// The training vector of response values.
        /// </summary>
        public Vector Y { get; private set; }

        /// <summary>
        /// The training matrix of predictor values. 
        /// </summary>
        public Matrix X { get; private set; }

        /// <summary>
        /// Gets the number of features in the training matrix of predictor values.
        /// </summary>
        public int NumberOfFeatures => X.NumberOfColumns;

        /// <summary>
        /// Determines whether this is for regression or classification. Default = regression.
        /// </summary>
        public bool IsRegression { get; set; } = true;

        #endregion

        #region Methods

        /// <summary>
        /// Returns the indexes of the k-Nearest Neighbors. 
        /// </summary>
        /// <param name="X">The 1D array of predictors.</param>
        public int[]? GetNeighbors(double[] X)
        {
            return kNN(this.X, this.Y, new Matrix(X));
        }

        /// <summary>
        /// Returns the indexes of the k-Nearest Neighbors. 
        /// </summary>
        /// <param name="X">The 2D array of predictors.</param>
        public int[]? GetNeighbors(double[,] X)
        {
            return kNN(this.X, this.Y, new Matrix(X));
        }

        /// <summary>
        /// Returns the indexes of the k-Nearest Neighbors. 
        /// </summary>
        /// <param name="X">The matrix of predictors.</param>
        public int[]? GetNeighbors(Matrix X)
        {
            return kNN(this.X, this.Y, X);
        }

        /// <summary>
        /// Returns the prediction from k-Nearest neighbors.
        /// </summary>
        /// <param name="X">The 1D array of predictors.</param>
        public double[]? Predict(double[] X)
        {
            return kNNPredict(this.X, this.Y, new Matrix(X));
        }

        /// <summary>
        /// Returns the prediction from k-Nearest neighbors.
        /// </summary>
        /// <param name="X">The 2D array of predictors.</param>
        public double[]? Predict(double[,] X)
        {
            return kNNPredict(this.X, this.Y, new Matrix(X));
        }

        /// <summary>
        /// Returns the prediction from k-Nearest neighbors.
        /// </summary>
        /// <param name="X">The matrix of predictors.</param>
        public double[]? Predict(Matrix X)
        {
            return kNNPredict(this.X, this.Y, X);
        }

        /// <summary>
        /// Returns the bootstrapped prediction from k-Nearest neighbors.
        /// </summary>
        /// <param name="X">The 1D array of predictors.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public double[]? BootstrapPredict(double[] X, int seed = -1)
        {
            return kNNBootstrapPredict(this.X, this.Y, new Matrix(X), seed);
        }

        /// <summary>
        /// Returns the bootstrapped prediction from k-Nearest neighbors.
        /// </summary>
        /// <param name="X">The 2D array of predictors.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public double[]? BootstrapPredict(double[,] X, int seed = -1)
        {
            return kNNBootstrapPredict(this.X, this.Y, new Matrix(X), seed);
        }

        /// <summary>
        /// Returns the bootstrapped prediction from k-Nearest neighbors.
        /// </summary>
        /// <param name="X">The matrix of predictors.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public double[]? BootstrapPredict(Matrix X, int seed = -1)
        {
            return kNNBootstrapPredict(this.X, this.Y, X, seed);
        }

        /// <summary>
        /// Returns the bootstrapped prediction intervals in a 2D array with columns: lower, median, upper, mean. 
        /// </summary>
        /// <param name="X">The 1D array of predictors.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="realizations">The number of bootstrap realizations. Default = 1,000.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,]? PredictionIntervals(double[] X, int seed = -1, int realizations = 1000, double alpha = 0.1)
        {
            return kNNPredictionIntervals(this.X, this.Y, new Matrix(X), seed, realizations, alpha);
        }

        /// <summary>
        /// Returns the bootstrapped prediction intervals in a 2D array with columns: lower, median, upper, mean. 
        /// </summary>
        /// <param name="X">The 2D array of predictors.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="realizations">The number of bootstrap realizations. Default = 1,000.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,]? PredictionIntervals(double[,] X, int seed = -1, int realizations = 1000, double alpha = 0.1)
        {
            return kNNPredictionIntervals(this.X, this.Y, new Matrix(X), seed, realizations, alpha);
        }

        /// <summary>
        /// Returns the bootstrapped prediction intervals in a 2D array with columns: lower, median, upper, mean. 
        /// </summary>
        /// <param name="X">The matrix of predictors.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="realizations">The number of bootstrap realizations. Default = 1,000.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        public double[,]? PredictionIntervals(Matrix X, int seed = -1, int realizations = 1000, double alpha = 0.1)
        {
            return kNNPredictionIntervals(this.X, this.Y, X, seed, realizations, alpha);
        }

        /// <summary>
        /// Returns the indexes of the k-Nearest neighbors.
        /// </summary>
        /// <param name="xTrain">The training matrix of predictors.</param>
        /// <param name="yTrain">The training response vector.</param>
        /// <param name="xTest">The test matrix of predictors</param>
        private int[]? kNN(Matrix xTrain, Vector yTrain, Matrix xTest)
        {
            // The guard must compare the query to the training matrix, matching kNNPredict; a query
            // with the wrong column count would otherwise compute partial-dimension distances or
            // index past the end of the training rows.
            if (xTest.NumberOfColumns != xTrain.NumberOfColumns) return null!;
            int R = xTest.NumberOfRows;
            var result = new int[R * K];
            if (K < 1 || K > xTrain.NumberOfRows)
            {
                // Outside 1 <= K <= n the full sort reproduces the existing failure modes.
                for (int i = 0; i < R; i++)
                {
                    var items = SortedItems(xTest.Row(i), xTrain);
                    for (int j = 0; j < K; j++)
                        result[i * K + j] = items[j].Index;
                }
                return result;
            }
            ForEachQuery(xTest, xTrain, (i, nearest) =>
            {
                for (int j = 0; j < K; j++)
                    result[i * K + j] = nearest[j].Index;
            });

            return result;
        }

        /// <summary>
        /// Returns the prediction from k-Nearest neighbors.
        /// </summary>
        /// <param name="xTrain">The training matrix of predictors.</param>
        /// <param name="yTrain">The training response vector.</param>
        /// <param name="xTest">The test matrix of predictors</param>
        private double[]? kNNPredict(Matrix xTrain, Vector yTrain, Matrix xTest)
        {
            if (xTest.NumberOfColumns != xTrain.NumberOfColumns) return null!;
            int R = xTest.NumberOfRows;
            var result = new double[R];
            if (K < 1 || K > xTrain.NumberOfRows)
            {
                // Outside 1 <= K <= n the full sort reproduces the existing failure modes.
                for (int i = 0; i < R; i++)
                    result[i] = Aggregate(SortedItems(xTest.Row(i), xTrain), yTrain);
                return result;
            }
            ForEachQuery(xTest, xTrain, (i, nearest) => result[i] = Aggregate(nearest, yTrain));

            return result;
        }

        /// <summary>
        /// Returns the regression or classification prediction from the nearest items.
        /// </summary>
        /// <param name="items">The items in ascending (distance, index) order; the first <see cref="K"/> are used.</param>
        /// <param name="yTrain">The training response vector.</param>
        /// <returns>The inverse squared distance weighted average (the mean of the exact matches when any neighbor is at zero distance), or the most common class.</returns>
        private double Aggregate(kNNItem[] items, Vector yTrain)
        {
            var knn = new double[K];
            if (IsRegression == true)
            {
                // Neighbors at zero distance take all the weight, averaged equally, as in scikit-learn's
                // distance weighting; an inverse squared distance weight would be infinite.
                int exact = 0;
                double exactSum = 0;
                for (int j = 0; j < K; j++)
                {
                    if (items[j].Distance == 0)
                    {
                        exact++;
                        exactSum += yTrain[items[j].Index];
                    }
                }
                if (exact > 0)
                    return exactSum / exact;

                // get the inverse distance weighted average of the values of k nearest neighbors
                double sum = 0;
                double avg = 0;
                for (int j = 0; j < K; j++)
                {
                    double w = items[j].Distance > 0 ? 1d / Tools.Sqr(items[j].Distance) : 1;
                    knn[j] = yTrain[items[j].Index] * w;
                    sum += w;
                }
                for (int j = 0; j < K; j++)
                    avg += knn[j] / sum;
                return avg;
            }

            // get the most common value
            for (int j = 0; j < K; j++)
                knn[j] = yTrain[items[j].Index];
            return knn.GroupBy(i => i).OrderByDescending(grp => grp.Count()).Select(grp => grp.Key).First();
        }

        /// <summary>
        /// Returns every training item sorted in ascending (distance, index) order.
        /// </summary>
        /// <param name="point">The query point.</param>
        /// <param name="xTrain">The training matrix of predictors.</param>
        /// <returns>The sorted items.</returns>
        private static kNNItem[] SortedItems(double[] point, Matrix xTrain)
        {
            var items = new kNNItem[xTrain.NumberOfRows];
            RangeParallel.For(xTrain.NumberOfRows, (long)xTrain.NumberOfRows * xTrain.NumberOfColumns, idx =>
            {
                items[idx].Index = idx;
                items[idx].Distance = Tools.Distance(point, xTrain.Row(idx));
            });
            Array.Sort(items, Compare);
            return items;
        }

        /// <summary>
        /// Orders items by distance, then by training index.
        /// </summary>
        /// <param name="a">The first item.</param>
        /// <param name="b">The second item.</param>
        /// <returns>A negative value, zero or a positive value as <paramref name="a"/> precedes, equals or follows <paramref name="b"/>.</returns>
        /// <remarks>
        /// Array.Sort is an unstable introspective sort, so the training index is an explicit
        /// secondary key. Duplicated rows and coded features produce exact distance ties, and without
        /// the secondary key the selected neighbors would be implementation-defined.
        /// </remarks>
        private static int Compare(kNNItem a, kNNItem b)
        {
            int c = a.Distance.CompareTo(b.Distance);
            return c != 0 ? c : a.Index.CompareTo(b.Index);
        }

        /// <summary>
        /// Finds the <see cref="K"/> nearest training items of every query row, for
        /// 1 &lt;= <see cref="K"/> &lt;= n, and passes them to <paramref name="consume"/>.
        /// </summary>
        /// <param name="xTest">The test matrix of predictors.</param>
        /// <param name="xTrain">The training matrix of predictors.</param>
        /// <param name="consume">Receives the query row index and its nearest items in ascending (distance, index) order.</param>
        /// <remarks>
        /// <see cref="Compare"/> is a strict total order on the training items, so the K smallest
        /// items kept by insertion are exactly the first K of the fully sorted array. With at least as
        /// many queries as processors the queries run in parallel and each computes its distances
        /// serially; otherwise the queries run in turn and the distances are computed in parallel.
        /// Each query writes only its own outputs, so the results do not depend on the schedule.
        /// The query loop keeps a <see cref="Parallel"/> loop because each worker reuses one
        /// neighbor buffer, which <see cref="RangeParallel"/> has no per-worker state for.
        /// </remarks>
        private void ForEachQuery(Matrix xTest, Matrix xTrain, Action<int, kNNItem[]> consume)
        {
            int R = xTest.NumberOfRows;
            int n = xTrain.NumberOfRows;
            var rows = new double[n][];
            for (int idx = 0; idx < n; idx++)
                rows[idx] = xTrain.Row(idx);
            if (R >= Environment.ProcessorCount)
            {
                Parallel.For(0, R, () => new kNNItem[K], (i, state, nearest) =>
                {
                    var point = xTest.Row(i);
                    int count = 0;
                    for (int idx = 0; idx < n; idx++)
                        count = Insert(nearest, count, new kNNItem { Index = idx, Distance = Tools.Distance(point, rows[idx]) });
                    consume(i, nearest);
                    return nearest;
                }, nearest => { });
            }
            else
            {
                var distances = new double[n];
                var nearest = new kNNItem[K];
                for (int i = 0; i < R; i++)
                {
                    var point = xTest.Row(i);
                    RangeParallel.For(n, (long)n * xTrain.NumberOfColumns, idx => distances[idx] = Tools.Distance(point, rows[idx]));
                    int count = 0;
                    for (int idx = 0; idx < n; idx++)
                        count = Insert(nearest, count, new kNNItem { Index = idx, Distance = distances[idx] });
                    consume(i, nearest);
                }
            }
        }

        /// <summary>
        /// Inserts an item into the sorted buffer of the smallest items seen so far.
        /// </summary>
        /// <param name="nearest">The sorted buffer.</param>
        /// <param name="count">The number of items in the buffer.</param>
        /// <param name="item">The candidate item.</param>
        /// <returns>The new number of items in the buffer.</returns>
        private static int Insert(kNNItem[] nearest, int count, kNNItem item)
        {
            int k = nearest.Length;
            if (count == k && Compare(item, nearest[k - 1]) >= 0)
                return count;
            int pos = count < k ? count : k - 1;
            while (pos > 0 && Compare(item, nearest[pos - 1]) < 0)
            {
                nearest[pos] = nearest[pos - 1];
                pos--;
            }
            nearest[pos] = item;
            return count < k ? count + 1 : count;
        }

        /// <summary>
        /// Returns the bootstrapped prediction from k-Nearest neighbors.
        /// </summary>
        /// <param name="xTrain">The training matrix of predictors.</param>
        /// <param name="yTrain">The training response vector.</param>
        /// <param name="xTest">The test matrix of predictors</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        private double[]? kNNBootstrapPredict(Matrix xTrain, Vector yTrain, Matrix xTest, int seed = -1)
        {
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var idxs = rnd.NextIntegers(0, xTrain.NumberOfRows, xTrain.NumberOfRows);
            var bootX = new Matrix(xTrain.NumberOfRows, xTrain.NumberOfColumns);
            var bootY = new Vector(yTrain.Length);
            for (int i = 0; i < xTrain.NumberOfRows; i++)
            {
                for (int j = 0; j < xTrain.NumberOfColumns; j++)
                {
                    bootX[i, j] = xTrain[idxs[i], j];
                }
                bootY[i] = yTrain[idxs[i]];
            }
            return kNNPredict(bootX, bootY, xTest);
        }

        /// <summary>
        /// Returns the bootstrapped prediction intervals in a 2D array with columns: lower, median, upper, mean. 
        /// </summary>
        /// <param name="xTrain">The training matrix of predictors.</param>
        /// <param name="yTrain">The training response vector.</param>
        /// <param name="xTest">The test matrix of predictors</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="realizations">The number of bootstrap realizations. Default = 1,000.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% confidence intervals.</param>
        private double[,] kNNPredictionIntervals(Matrix xTrain, Vector yTrain, Matrix xTest, int seed = -1, int realizations = 1000, double alpha = 0.1)
        {
            var percentiles = new double[] { alpha / 2d, 0.5, 1d - alpha / 2d };
            var output = new double[xTest.NumberOfRows, 4]; // lower, median, upper, mean

            var bootResults = new double[xTest.NumberOfRows, realizations];
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var seeds = rnd.NextIntegers(realizations);

            // Bootstrap the predictions
            RangeParallel.For(realizations, (long)realizations * xTest.NumberOfRows * xTrain.NumberOfRows * xTrain.NumberOfColumns,
                idx => { bootResults.SetColumn(idx, kNNBootstrapPredict(xTrain, yTrain, xTest, seeds[idx])!); });

            // Process results
            RangeParallel.For(xTest.NumberOfRows, (long)xTest.NumberOfRows * realizations, idx =>
            {
                var values = bootResults.GetRow(idx);
                Array.Sort(values);

                // Record percentiles for CIs
                for (int j = 0; j < percentiles.Length; j++)
                    output[idx, j] = Statistics.Percentile(values, percentiles[j], true);

                // The mean is accumulated sequentially so the reduction is deterministic on every
                // host regardless of processor count.
                output[idx, 3] = Statistics.Mean(values);
            });

            return output;
        }

        /// <summary>
        /// A structure for storing a k-NN item.
        /// </summary>
        private struct kNNItem
        {
            /// <summary>
            /// The index of the item.
            /// </summary>
            public int Index;

            /// <summary>
            /// The distance of the item.
            /// </summary>
            public double Distance;
        }

        #endregion

    }
}
