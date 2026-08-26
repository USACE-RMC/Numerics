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

        #endregion

        #region Members

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
            Root = GrowTree(X, Y);
            IsTrained = true;
        }

        /// <summary>
        /// Grow the decision tree recursively. 
        /// </summary>
        /// <param name="xTrain">The training matrix of predictor values.</param>
        /// <param name="yTrain">The training vector of response values.</param>
        /// <param name="depth">The depth of the recursion.</param>
        private DecisionNode GrowTree(Matrix xTrain, Vector yTrain, int depth = 0)
        {
            int numberOfSamples = xTrain.NumberOfRows;
            int numberOfLabels = IsRegression ? yTrain.Length : yTrain.ToList().Distinct().Count();

            // The feature subset is drawn for every node, split or leaf, so the generator consumes
            // exactly one draw per node and the draw schedule is independent of the stopping conditions.
            var featureIdxs = Random.NextIntegers(0, Dimensions, Features, false);

            // Leaf conditions that need no split search
            if (depth >= MaxDepth || numberOfLabels <= 1 || numberOfSamples < MinimumSplitSize)
            {
                return CreateLeaf(yTrain);
            }

            // Find the best split
            int bestIndex = 0; double bestThreshold = 0;
            BestSplit(xTrain, yTrain.ToArray(), featureIdxs, out bestIndex, out bestThreshold);
            if (bestIndex == -1)
            {
                return CreateLeaf(yTrain);
            }


            // Create child nodes
            var leftIdxs = new List<int>();
            var rightIdxs = new List<int>();
            Split(xTrain.Column(bestIndex), bestThreshold, out leftIdxs, out rightIdxs);

            // Split to the left
            var xLeft = new Matrix(leftIdxs.Count, xTrain.NumberOfColumns);
            var yLeft = new Vector(leftIdxs.Count);
            for (int i = 0; i < leftIdxs.Count; i++)
            {
                yLeft[i] = yTrain[leftIdxs[i]];
                for (int j = 0; j < xTrain.NumberOfColumns; j++)
                {
                    xLeft[i, j] = xTrain[leftIdxs[i], j];
                }
            }
            var left = GrowTree(xLeft, yLeft, depth + 1);

            // Split to the right
            var xRight = new Matrix(rightIdxs.Count, xTrain.NumberOfColumns);
            var yRight = new Vector(rightIdxs.Count);
            for (int i = 0; i < rightIdxs.Count; i++)
            {
                yRight[i] = yTrain[rightIdxs[i]];
                for (int j = 0; j < xTrain.NumberOfColumns; j++)
                {
                    xRight[i, j] = xTrain[rightIdxs[i], j];
                }
            }
            var right = GrowTree(xRight, yRight, depth + 1);

            // Return decision node
            return new DecisionNode() { FeatureIndex = bestIndex, Threshold = bestThreshold, Left = left, Right = right };
        }

        /// <summary>
        /// Creates a leaf node from the response values.
        /// </summary>
        /// <param name="yTrain">The training vector of response values.</param>
        /// <returns>A leaf node holding the mean response for regression or the most common label for classification.</returns>
        private DecisionNode CreateLeaf(Vector yTrain)
        {
            if (IsRegression)
            {
                // If regression return the average of Y
                var avg = Tools.Mean(yTrain.ToArray());
                return new DecisionNode() { Value = avg, IsLeafNode = true };
            }
            else
            {
                // If classification, return the most common value
                var most = yTrain.ToList().GroupBy(i => i).OrderByDescending(grp => grp.Count()).Select(grp => grp.Key).First();
                return new DecisionNode() { Value = most, IsLeafNode = true };
            }
        }

        /// <summary>
        /// Returns the best split feature index and threshold.
        /// </summary>
        /// <param name="xTrain">The matrix of predictor values.</param>
        /// <param name="yTrain">The array of y-values.</param>
        /// <param name="indices">The feature indexes to evaluate.</param>
        /// <param name="bestFeatureIndex">Output. The best feature index.</param>
        /// <param name="bestThreshold">Output. The best threshold for splitting the tree.</param>
        /// <remarks>
        /// The parent impurity is constant within a node, so it is evaluated once and shared across
        /// every candidate. Duplicate candidate thresholds produce identical splits, so only the
        /// first occurrence of each value is evaluated; because ties break toward the earliest
        /// candidate, the selected split is the same one an exhaustive scan selects.
        /// </remarks>
        private void BestSplit(Matrix xTrain, double[] yTrain, int[] indices, out int bestFeatureIndex, out double bestThreshold)
        {
            double best = double.MinValue;
            bestFeatureIndex = -1;
            bestThreshold = 0;

            double parentImpurity = IsRegression ? Statistics.PopulationVariance(yTrain) : Entropy(yTrain);
            var seen = new HashSet<double>();

            for (int i = 0; i < indices.Length; i++)
            {
                var x = xTrain.Column(indices[i]);
                seen.Clear();
                for (int j = 0; j < x.Length; j++)
                {
                    if (!seen.Add(x[j]))
                    {
                        continue;
                    }
                    // Test if the split variance reduction or information gain
                    double performance = IsRegression ? VarianceReduction(x, yTrain, x[j], parentImpurity) : InformationGain(x, yTrain, x[j], parentImpurity);
                    // Keep track of the best value
                    if (performance > best)
                    {
                        best = performance;
                        bestFeatureIndex = indices[i];
                        bestThreshold = x[j];
                    }
                }
            }

        }

        /// <summary>
        /// Computes the variance reduction for the threshold.
        /// </summary>
        /// <param name="x">The column of x-values.</param>
        /// <param name="y">The column of y-values.</param>
        /// <param name="threshold">The split threshold.</param>
        /// <param name="parentVariance">The population variance of the parent node response values.</param>
        /// <returns>The variance reduction of the split, or double.MinValue when the threshold places every sample in one child.</returns>
        /// <remarks>
        /// The child variances are accumulated in sample order with the same recurrence as
        /// <see cref="Statistics.PopulationVariance(System.Collections.Generic.IList{double})"/>, so
        /// the reduction is identical to evaluating that method on materialized child arrays.
        /// </remarks>
        private double VarianceReduction(double[] x, double[] y, double threshold, double parentVariance)
        {
            // Stream both child variances in one pass over the sample
            int countLeft = 0, countRight = 0;
            double sumLeft = 0, sumRight = 0, accLeft = 0, accRight = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double v = y[i];
                if (x[i] <= threshold)
                {
                    countLeft++;
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
                }
                else
                {
                    countRight++;
                    if (countRight == 1)
                    {
                        sumRight = v;
                    }
                    else
                    {
                        sumRight += v;
                        double diff = countRight * v - sumRight;
                        accRight += diff * diff / (countRight * (countRight - 1.0));
                    }
                }
            }

            if (countLeft == 0 || countRight == 0)
                return double.MinValue;

            // calculate the weighted average variance of the children
            var n = (double)y.Length;
            var nLeft = (double)countLeft;
            var nRight = (double)countRight;
            var varLeft = accLeft / countLeft;
            var varRight = accRight / countRight;
            var childVariance = nLeft / n * varLeft + nRight / n * varRight;

            // Return variance reduction
            return parentVariance - childVariance;
        }

        /// <summary>
        /// Returns the information gain of the split threshold.
        /// </summary>
        /// <param name="x">The column of x-values.</param>
        /// <param name="y">The column of y-values.</param>
        /// <param name="threshold">The split threshold.</param>
        /// <param name="parentEntropy">The entropy of the parent node response values.</param>
        /// <returns>The information gain of the split, or double.MinValue when the threshold places every sample in one child.</returns>
        /// <remarks>
        /// Label probabilities are the per-child label counts divided by the child size, and entropy
        /// terms are accumulated in sample order, so each child entropy is identical to evaluating
        /// <see cref="Entropy(double[])"/> on a materialized child array.
        /// </remarks>
        private double InformationGain(double[] x, double[] y, double threshold, double parentEntropy)
        {
            // Count the labels on each side of the threshold
            var leftCounts = new Dictionary<double, int>();
            var rightCounts = new Dictionary<double, int>();
            int countLeft = 0, countRight = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double v = y[i];
                if (x[i] <= threshold)
                {
                    countLeft++;
                    if (!double.IsNaN(v))
                    {
                        leftCounts.TryGetValue(v, out int c);
                        leftCounts[v] = c + 1;
                    }
                }
                else
                {
                    countRight++;
                    if (!double.IsNaN(v))
                    {
                        rightCounts.TryGetValue(v, out int c);
                        rightCounts[v] = c + 1;
                    }
                }
            }

            if (countLeft == 0 || countRight == 0)
                return double.MinValue;

            // calculate the weighted average entropy of children
            double sumLeft = 0, sumRight = 0;
            for (int i = 0; i < x.Length; i++)
            {
                double v = y[i];
                if (double.IsNaN(v))
                    continue;
                if (x[i] <= threshold)
                {
                    double p = (double)leftCounts[v] / countLeft;
                    sumLeft += p * Math.Log(p);
                }
                else
                {
                    double p = (double)rightCounts[v] / countRight;
                    sumRight += p * Math.Log(p);
                }
            }

            var n = (double)y.Length;
            var nl = (double)countLeft;
            var nr = (double)countRight;
            var childrenE = nl / n * (-sumLeft) + nr / n * (-sumRight);

            return parentEntropy - childrenE;
        }

        /// <summary>
        /// Computes the entropy for a vector of classification labels.
        /// </summary>
        /// <param name="y">The column of y-values.</param>
        /// <returns>The entropy of the labels.</returns>
        /// <remarks>
        /// The empirical probability of each label is its count divided by the sample length, and
        /// terms are accumulated in sample order. Labels of NaN carry an empirical probability of
        /// zero and contribute nothing.
        /// </remarks>
        private double Entropy(double[] y)
        {
            var counts = new Dictionary<double, int>();
            for (int i = 0; i < y.Length; i++)
            {
                double v = y[i];
                if (!double.IsNaN(v))
                {
                    counts.TryGetValue(v, out int c);
                    counts[v] = c + 1;
                }
            }

            double sum = 0;
            for (int i = 0; i < y.Length; i++)
            {
                double v = y[i];
                if (double.IsNaN(v))
                    continue;
                double p = (double)counts[v] / y.Length;
                sum += p * Math.Log(p);
            }
            return -sum;
        }

        /// <summary>
        /// Splits the x-column based on the threshold.
        /// </summary>
        /// <param name="x">The column of x-values.</param>
        /// <param name="threshold">The split threshold.</param>
        /// <param name="leftIndices">Output. A list of left indexes.</param>
        /// <param name="rightIndices">Output. A list of right indexes.</param>
        private void Split(double[] x, double threshold, out List<int> leftIndices, out List<int> rightIndices)
        {
            leftIndices = new List<int>();
            rightIndices = new List<int>();
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] <= threshold)
                {
                    leftIndices.Add(i);
                }
                else
                {
                    rightIndices.Add(i);
                }
            }
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
