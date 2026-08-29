using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// The Decision Tree learning algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b> 
    /// </para>
    /// <para> <b> References: </b> </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Decision_tree_learning" />
    /// </para>
    /// </remarks>
    public class DecisionTree
    {

        #region Construction

        /// <summary>
        /// Create new Decision Tree.
        /// </summary>
        /// <param name="x">The training 1D array of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public DecisionTree(double[] x, double[] y, int seed = -1)
        {
            // Set inputs
            Y = new Vector(y);
            X = new Matrix(x);
            Dimensions = X.NumberOfColumns;
            Features = Math.Max(1, Dimensions - 1);
            Root = new DecisionNode();
            Random = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");

        }

        /// <summary>
        /// Create new Decision Tree.
        /// </summary>
        /// <param name="x">The training 2D array of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public DecisionTree(double[,] x, double[] y, int seed = -1)
        {
            // Set inputs
            Y = new Vector(y);
            X = new Matrix(x);
            Dimensions = X.NumberOfColumns;
            Features = Math.Max(1, Dimensions - 1);
            Root = new DecisionNode();
            Random = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");

        }

        /// <summary>
        /// Create new Decision Tree.
        /// </summary>
        /// <param name="x">The training matrix of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        public DecisionTree(Matrix x, Vector y, int seed = -1)
        {
            // Set inputs
            Y = y;
            X = x;
            Dimensions = X.NumberOfColumns;
            Features = Math.Max(1, Dimensions - 1);
            Root = new DecisionNode();
            Random = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");

        }

        /// <summary>
        /// Create new Decision Tree that trains on a caller-specified multiset of shared training rows.
        /// </summary>
        /// <param name="x">The training matrix of predictor values, shared with the caller.</param>
        /// <param name="y">The training response vector, shared with the caller.</param>
        /// <param name="sampleIndices">The rows of <paramref name="x"/> that form the training multiset, in training order.</param>
        /// <param name="seed">The prng seed.</param>
        internal DecisionTree(Matrix x, Vector y, int[] sampleIndices, int seed)
        {
            // Set inputs
            Y = y;
            X = x;
            Dimensions = X.NumberOfColumns;
            Features = Math.Max(1, Dimensions - 1);
            Root = new DecisionNode();
            Random = new MersenneTwister(seed);
            _sampleIndices = sampleIndices;

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");

        }

        #endregion

        #region Members

        // The training rows this tree trains on; null trains on every row.
        private readonly int[]? _sampleIndices;

        // Split-search scratch, allocated per training run and shared down the single-threaded recursion
        private double[] _keyScratch = Array.Empty<double>();
        private double[] _valScratch = Array.Empty<double>();
        private double[] _rightAccScratch = Array.Empty<double>();
        private int[] _partitionScratch = Array.Empty<int>();

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
        /// The root node of the decision tree.
        /// </summary>
        public DecisionNode Root { get; private set; }

        /// <summary>
        /// The training vector of response values.
        /// </summary>
        public Vector Y { get; private set; }

        /// <summary>
        /// The training matrix of predictor values. 
        /// </summary>
        public Matrix X { get; private set; }

        /// <summary>
        /// Determines whether this is for regression or classification. Default = regression.
        /// </summary>
        public bool IsRegression { get; set; } = true;

        /// <summary>
        /// Determines if the decision tree has been trained.
        /// </summary>
        public bool IsTrained { get; private set; } = false;

        #endregion

        #region Methods

        /// <summary>
        /// Train the decision tree.
        /// </summary>
        public void Train()
        {
            IsTrained = false;
            Features = Math.Min(Features, Dimensions);
            int n = X.NumberOfRows;

            // Row indices define the training multiset, so a bootstrap tree shares the parent
            // matrix rather than holding its own copy.
            var indices = new int[_sampleIndices == null ? n : _sampleIndices.Length];
            if (_sampleIndices == null)
            {
                for (int i = 0; i < n; i++) indices[i] = i;
            }
            else
            {
                Array.Copy(_sampleIndices, indices, _sampleIndices.Length);
            }

            // Scratch buffers sized once per training run and shared down the recursion, which is
            // single threaded within a tree.
            int m = indices.Length;
            _keyScratch = new double[m];
            _valScratch = new double[m];
            _rightAccScratch = new double[m + 1];
            _partitionScratch = new int[m];

            Root = GrowTree(indices, 0, m, 0);

            _keyScratch = _valScratch = _rightAccScratch = Array.Empty<double>();
            _partitionScratch = Array.Empty<int>();
            IsTrained = true;
        }

        /// <summary>
        /// Grow the decision tree recursively over a range of training row indices.
        /// </summary>
        /// <param name="indices">The training row indices, partitioned in place as the tree grows.</param>
        /// <param name="lo">The inclusive start of the node's range within <paramref name="indices"/>.</param>
        /// <param name="hi">The exclusive end of the node's range within <paramref name="indices"/>.</param>
        /// <param name="depth">The depth of the recursion.</param>
        /// <returns>The subtree root for the range: a decision node when a split is found, otherwise a leaf.</returns>
        private DecisionNode GrowTree(int[] indices, int lo, int hi, int depth)
        {
            int numberOfSamples = hi - lo;
            // Count distinct responses for BOTH modes so the pure-node guard below can fire —
            // scikit-learn's rule: a node is never split once it is pure. Substituting the sample
            // count makes the guard redundant with the minimum split size, and a zero-gain split of
            // a pure node beats the double.MinValue seed, so a default regression tree recurses to
            // one observation per leaf.
            int numberOfLabels = CountDistinctLabels(indices, lo, hi);

            // The feature subset is drawn for every node, split or leaf, so the generator consumes
            // exactly one draw per node and the draw schedule is independent of the stopping conditions.
            var featureIdxs = Random.NextIntegers(0, Dimensions, Features, false);

            // Leaf conditions that need no split search
            if (depth >= MaxDepth || numberOfLabels <= 1 || numberOfSamples < MinimumSplitSize)
            {
                return CreateLeaf(indices, lo, hi);
            }

            // Find the best split
            BestSplit(indices, lo, hi, featureIdxs, out int bestIndex, out double bestThreshold);
            if (bestIndex == -1)
            {
                return CreateLeaf(indices, lo, hi);
            }

            // Partition the index range in place, keeping the relative row order stable on both
            // sides. Rows with a NaN predictor fail x <= t and fall to the right.
            int split = StablePartition(indices, lo, hi, bestIndex, bestThreshold);

            var left = GrowTree(indices, lo, split, depth + 1);
            var right = GrowTree(indices, split, hi, depth + 1);

            // Return decision node
            return new DecisionNode() { FeatureIndex = bestIndex, Threshold = bestThreshold, Left = left, Right = right };
        }

        /// <summary>
        /// Counts the distinct response values within an index range — class labels for classification,
        /// response values for regression.
        /// </summary>
        /// <param name="indices">The training row indices.</param>
        /// <param name="lo">The inclusive start of the range.</param>
        /// <param name="hi">The exclusive end of the range.</param>
        /// <returns>The number of distinct values, counting NaN responses as one value.</returns>
        private int CountDistinctLabels(int[] indices, int lo, int hi)
        {
            var labels = new HashSet<double>();
            for (int i = lo; i < hi; i++)
                labels.Add(Y[indices[i]]);
            return labels.Count;
        }

        /// <summary>
        /// Creates a leaf node from the responses within an index range.
        /// </summary>
        /// <param name="indices">The training row indices.</param>
        /// <param name="lo">The inclusive start of the range.</param>
        /// <param name="hi">The exclusive end of the range.</param>
        /// <returns>A leaf node holding the mean response for regression or the most common label for classification.</returns>
        /// <remarks>
        /// The classification mode is the first label in sample order to reach the maximum count, so
        /// ties resolve to the earliest observed label.
        /// </remarks>
        private DecisionNode CreateLeaf(int[] indices, int lo, int hi)
        {
            if (IsRegression)
            {
                // The average of Y over the node, accumulated in sample order
                double sum = 0;
                for (int i = lo; i < hi; i++)
                    sum += Y[indices[i]];
                return new DecisionNode() { Value = sum / (hi - lo), IsLeafNode = true };
            }
            else
            {
                // The most common label over the node
                var counts = new Dictionary<double, int>();
                for (int i = lo; i < hi; i++)
                {
                    double v = Y[indices[i]];
                    counts.TryGetValue(v, out int c);
                    counts[v] = c + 1;
                }
                double most = double.NaN;
                int best = -1;
                for (int i = lo; i < hi; i++)
                {
                    double v = Y[indices[i]];
                    if (counts[v] > best)
                    {
                        best = counts[v];
                        most = v;
                    }
                }
                return new DecisionNode() { Value = most, IsLeafNode = true };
            }
        }

        /// <summary>
        /// Returns the best split feature index and threshold for an index range.
        /// </summary>
        /// <param name="indices">The training row indices.</param>
        /// <param name="lo">The inclusive start of the node's range.</param>
        /// <param name="hi">The exclusive end of the node's range.</param>
        /// <param name="featureIdxs">The feature indexes to evaluate.</param>
        /// <param name="bestFeatureIndex">Output. The best feature index, or -1 when no candidate separates the node.</param>
        /// <param name="bestThreshold">Output. The best threshold for splitting the tree.</param>
        /// <remarks>
        /// <para>
        /// Each candidate feature is sorted once and every distinct value is evaluated as a
        /// threshold in a single sweep that maintains running child statistics, so a node costs
        /// O(features × m log m) rather than one full pass per candidate. Rows with a NaN predictor
        /// fail every x &lt;= t comparison and are held on the right side throughout the sweep.
        /// </para>
        /// <para>
        /// Regression splits maximize the variance reduction, with child variances accumulated by
        /// the Youngs and Cramer (1971) updating formula, the same recurrence used by
        /// <see cref="Statistics.PopulationVariance(System.Collections.Generic.IList{double})"/>.
        /// Classification splits maximize the information gain of the per-sample entropy, in which
        /// each label's term is weighted by its own empirical probability.
        /// </para>
        /// <para>
        /// Rows with equal predictor values retain their order in the current node. Candidate
        /// thresholds are evaluated in ascending order, so the smallest threshold wins an exact
        /// within-feature gain tie; the first sampled feature wins an exact cross-feature tie.
        /// </para>
        /// <para><b> References: </b></para>
        /// <list type="bullet">
        /// <item><description>
        /// Youngs, E. A. and Cramer, E. M. (1971). Some results relevant to choice of sum and sum-of-product algorithms. Technometrics, 13(3), 657-665.
        /// </description></item>
        /// </list>
        /// </remarks>
        private void BestSplit(int[] indices, int lo, int hi, int[] featureIdxs, out int bestFeatureIndex, out double bestThreshold)
        {
            double best = double.MinValue;
            bestFeatureIndex = -1;
            bestThreshold = 0;

            var keys = _keyScratch;
            var vals = _valScratch;
            var order = _partitionScratch;
            var sortedVals = _rightAccScratch;

            for (int f = 0; f < featureIdxs.Length; f++)
            {
                int feature = featureIdxs[f];

                // Gather the feature and response pairs, holding NaN predictors aside on the right
                int valid = 0, nanCount = 0;
                for (int i = lo; i < hi; i++)
                {
                    double key = X[indices[i], feature];
                    if (double.IsNaN(key))
                    {
                        nanCount++;
                    }
                    else
                    {
                        keys[valid] = key;
                        vals[valid] = Y[indices[i]];
                        order[valid] = valid;
                        valid++;
                    }
                }
                if (valid == 0)
                    continue;

                Array.Sort(keys, order, 0, valid);
                int runStart = 0;
                while (runStart < valid)
                {
                    int runEnd = runStart + 1;
                    while (runEnd < valid && keys[runEnd] == keys[runStart])
                        runEnd++;
                    Array.Sort(order, runStart, runEnd - runStart);
                    runStart = runEnd;
                }
                for (int i = 0; i < valid; i++)
                    sortedVals[i] = vals[order[i]];
                Array.Copy(sortedVals, vals, valid);

                double performance;
                double threshold;
                if (IsRegression)
                {
                    if (!SweepVariance(indices, lo, hi, feature, valid, nanCount, out performance, out threshold))
                        continue;
                }
                else
                {
                    if (!SweepEntropy(indices, lo, hi, feature, valid, nanCount, out performance, out threshold))
                        continue;
                }

                // Keep track of the best value
                if (performance > best)
                {
                    best = performance;
                    bestFeatureIndex = feature;
                    bestThreshold = threshold;
                }
            }
        }

        /// <summary>
        /// Sweeps the sorted candidate thresholds of one feature and returns the best variance reduction.
        /// </summary>
        /// <param name="indices">The training row indices.</param>
        /// <param name="lo">The inclusive start of the node's range.</param>
        /// <param name="hi">The exclusive end of the node's range.</param>
        /// <param name="feature">The feature under evaluation.</param>
        /// <param name="valid">The number of sorted non-NaN feature values in the scratch buffers.</param>
        /// <param name="nanCount">The number of rows whose feature value is NaN, held on the right side.</param>
        /// <param name="bestGain">Output. The largest variance reduction over the candidate thresholds.</param>
        /// <param name="bestThreshold">Output. The threshold attaining <paramref name="bestGain"/>.</param>
        /// <returns>True when at least one threshold produces two non-empty children.</returns>
        private bool SweepVariance(int[] indices, int lo, int hi, int feature, int valid, int nanCount, out double bestGain, out double bestThreshold)
        {
            int m = hi - lo;
            var keys = _keyScratch;
            var vals = _valScratch;

            // Parent variance over the full node in sample order
            double parentVariance = RangeVariance(indices, lo, hi);

            // Right-side running statistics from the top of the sort down, seeded with the NaN rows
            var rightAcc = _rightAccScratch;
            int rightSeed = 0;
            double accR = 0, sumR = 0;
            if (nanCount > 0)
            {
                for (int i = lo; i < hi; i++)
                {
                    if (!double.IsNaN(X[indices[i], feature]))
                        continue;
                    rightSeed++;
                    double v = Y[indices[i]];
                    if (rightSeed == 1)
                    {
                        sumR = v;
                    }
                    else
                    {
                        sumR += v;
                        double diff = rightSeed * v - sumR;
                        accR += diff * diff / (rightSeed * (rightSeed - 1.0));
                    }
                }
            }
            rightAcc[valid] = accR;
            for (int i = valid - 1; i >= 0; i--)
            {
                int k = rightSeed + (valid - i);
                double v = vals[i];
                if (k == 1)
                {
                    sumR = v;
                }
                else
                {
                    sumR += v;
                    double diff = k * v - sumR;
                    accR += diff * diff / (k * (k - 1.0));
                }
                rightAcc[i] = accR;
            }

            // Forward sweep: grow the left side one sorted value at a time and score every distinct
            // boundary
            bestGain = double.MinValue;
            bestThreshold = 0;
            bool found = false;
            double n = m;
            int countLeft = 0;
            double sumLeft = 0, accLeft = 0;
            for (int i = 0; i < valid; i++)
            {
                countLeft++;
                double v = vals[i];
                if (countLeft == 1)
                {
                    sumLeft = v;
                }
                else
                {
                    sumLeft += v;
                    double diff = countLeft * v - sumLeft;
                    accLeft += diff * diff / (countLeft * (countLeft - 1.0));
                }

                // A boundary exists after the final occurrence of each distinct value
                if (i + 1 < valid && keys[i + 1] == keys[i])
                    continue;
                int countRight = m - countLeft;
                if (countRight == 0)
                    break;

                double varLeft = accLeft / countLeft;
                double varRight = rightAcc[i + 1] / countRight;
                double childVariance = countLeft / n * varLeft + countRight / n * varRight;
                double gain = parentVariance - childVariance;
                if (gain > bestGain)
                {
                    bestGain = gain;
                    bestThreshold = keys[i];
                    found = true;
                }
            }
            return found;
        }

        /// <summary>
        /// Sweeps the sorted candidate thresholds of one feature and returns the best information gain.
        /// </summary>
        /// <param name="indices">The training row indices.</param>
        /// <param name="lo">The inclusive start of the node's range.</param>
        /// <param name="hi">The exclusive end of the node's range.</param>
        /// <param name="feature">The feature under evaluation.</param>
        /// <param name="valid">The number of sorted non-NaN feature values in the scratch buffers.</param>
        /// <param name="nanCount">The number of rows whose feature value is NaN, held on the right side.</param>
        /// <param name="bestGain">Output. The largest information gain over the candidate thresholds.</param>
        /// <param name="bestThreshold">Output. The threshold attaining <paramref name="bestGain"/>.</param>
        /// <returns>True when at least one threshold produces two non-empty children.</returns>
        private bool SweepEntropy(int[] indices, int lo, int hi, int feature, int valid, int nanCount, out double bestGain, out double bestThreshold)
        {
            int m = hi - lo;
            var keys = _keyScratch;
            var vals = _valScratch;

            // Class counts for the full node; NaN labels carry zero probability and are excluded
            // from the counts while remaining in the child sizes
            var rightCounts = new Dictionary<double, int>();
            for (int i = lo; i < hi; i++)
            {
                double v = Y[indices[i]];
                if (!double.IsNaN(v))
                {
                    rightCounts.TryGetValue(v, out int c);
                    rightCounts[v] = c + 1;
                }
            }
            double parentEntropy = CountsEntropy(rightCounts, m);
            var leftCounts = new Dictionary<double, int>();

            bestGain = double.MinValue;
            bestThreshold = 0;
            bool found = false;
            double n = m;
            int countLeft = 0;
            for (int i = 0; i < valid; i++)
            {
                countLeft++;
                double v = vals[i];
                if (!double.IsNaN(v))
                {
                    leftCounts.TryGetValue(v, out int c);
                    leftCounts[v] = c + 1;
                    rightCounts[v] = rightCounts[v] - 1;
                }

                // A boundary exists after the final occurrence of each distinct value
                if (i + 1 < valid && keys[i + 1] == keys[i])
                    continue;
                int countRight = m - countLeft;
                if (countRight == 0)
                    break;

                double el = CountsEntropy(leftCounts, countLeft);
                double er = CountsEntropy(rightCounts, countRight);
                double childrenE = countLeft / n * el + countRight / n * er;
                double gain = parentEntropy - childrenE;
                if (gain > bestGain)
                {
                    bestGain = gain;
                    bestThreshold = keys[i];
                    found = true;
                }
            }
            return found;
        }

        /// <summary>
        /// Computes the per-sample entropy of a group from its label counts.
        /// </summary>
        /// <param name="counts">The label counts of the group.</param>
        /// <param name="size">The group size, including any rows whose label is NaN.</param>
        /// <returns>The entropy of the labels.</returns>
        /// <remarks>
        /// Each label contributes count × p × ln(p) with p equal to its count divided by the group
        /// size, matching a per-sample accumulation in which every sample carries the empirical
        /// probability of its own label.
        /// </remarks>
        private static double CountsEntropy(Dictionary<double, int> counts, int size)
        {
            double sum = 0;
            foreach (var pair in counts)
            {
                if (pair.Value <= 0)
                    continue;
                double p = (double)pair.Value / size;
                sum += pair.Value * p * Math.Log(p);
            }
            return -sum;
        }

        /// <summary>
        /// Computes the population variance of the responses within an index range.
        /// </summary>
        /// <param name="indices">The training row indices.</param>
        /// <param name="lo">The inclusive start of the range.</param>
        /// <param name="hi">The exclusive end of the range.</param>
        /// <returns>The population variance, accumulated in sample order with the same recurrence as
        /// <see cref="Statistics.PopulationVariance(System.Collections.Generic.IList{double})"/>.</returns>
        private double RangeVariance(int[] indices, int lo, int hi)
        {
            int count = 0;
            double sum = 0, acc = 0;
            for (int i = lo; i < hi; i++)
            {
                count++;
                double v = Y[indices[i]];
                if (count == 1)
                {
                    sum = v;
                }
                else
                {
                    sum += v;
                    double diff = count * v - sum;
                    acc += diff * diff / (count * (count - 1.0));
                }
            }
            return count == 0 ? double.NaN : acc / count;
        }

        /// <summary>
        /// Partitions an index range in place around a threshold, preserving the relative order on
        /// both sides.
        /// </summary>
        /// <param name="indices">The training row indices.</param>
        /// <param name="lo">The inclusive start of the range.</param>
        /// <param name="hi">The exclusive end of the range.</param>
        /// <param name="feature">The splitting feature.</param>
        /// <param name="threshold">The splitting threshold; rows with x &lt;= t go left and all others, NaN included, go right.</param>
        /// <returns>The index of the first right-side element.</returns>
        private int StablePartition(int[] indices, int lo, int hi, int feature, double threshold)
        {
            var scratch = _partitionScratch;
            int left = lo, rightCount = 0;
            for (int i = lo; i < hi; i++)
            {
                int row = indices[i];
                if (X[row, feature] <= threshold)
                {
                    indices[left++] = row;
                }
                else
                {
                    scratch[rightCount++] = row;
                }
            }
            Array.Copy(scratch, 0, indices, left, rightCount);
            return left;
        }

        /// <summary>
        /// Traverses the tree and return the leaf node value.
        /// </summary>
        /// <param name="x">The row of x-value predictors.</param>
        /// <param name="node">The decision node to traverse.</param>
        private double TraverseTree(double[] x, DecisionNode node)
        {
            if (node.IsLeafNode == true)
                return node.Value;
            if (x[node.FeatureIndex] <= node.Threshold)
                return node.Left != null ? TraverseTree(x, node.Left) : node.Value;
            return node.Right != null ? TraverseTree(x, node.Right) : node.Value;
        }

        /// <summary>
        /// Returns the prediction for a single row of predictors.
        /// </summary>
        /// <param name="x">The row of x-value predictors.</param>
        /// <returns>The leaf value the row traverses to.</returns>
        internal double PredictRow(double[] x)
        {
            return TraverseTree(x, Root);
        }

        /// <summary>
        /// Returns the prediction from the Decision Tree.
        /// </summary>
        /// <param name="X">The 1D array of predictors.</param>
        public double[]? Predict(double[] X)
        {
            return Predict(new Matrix(X));
        }

        /// <summary>
        /// Returns the prediction from the Decision Tree.
        /// </summary>
        /// <param name="X">The 2D array of predictors.</param>
        public double[]? Predict(double[,] X)
        {
            return Predict(new Matrix(X));
        }

        /// <summary>
        /// Returns the prediction from the Decision Tree.
        /// </summary>
        /// <param name="X">The matrix of predictors.</param>
        public double[]? Predict(Matrix X)
        {
            if (!IsTrained || X.NumberOfColumns != Dimensions) return null!;
            var result = new double[X.NumberOfRows];
            for (int i = 0; i < X.NumberOfRows; i++)
            {
                result[i] = TraverseTree(X.Row(i), Root);
            }
            return result;
        }

        #endregion

    }
}
