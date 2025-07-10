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

using Numerics.Data.Statistics;
using Numerics.Distributions;
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

            // Find the best split
            var featureIdxs = Random.NextIntegers(0, Dimensions, Features, false);
            int bestIndex = 0; double bestThreshold = 0;
            BestSplit(xTrain, yTrain.ToArray(), featureIdxs, out bestIndex, out bestThreshold);

            // Check stopping criteria
            if (bestIndex == -1 || depth >= MaxDepth || numberOfLabels <= 1 || numberOfSamples < MinimumSplitSize)
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
        /// Returns the best split feature index and threshold.
        /// </summary>
        /// <param name="xTrain">The matrix of predictor values.</param>
        /// <param name="yTrain">The array of y-values.</param>
        /// <param name="indices">The feature indexes to evaluate.</param>
        /// <param name="bestFeatureIndex">Output. The best feature index.</param>
        /// <param name="bestThreshold">Output. The best threshold for splitting the tree.</param>
        private void BestSplit(Matrix xTrain, double[] yTrain, int[] indices, out int bestFeatureIndex, out double bestThreshold)
        {
            double best = double.MinValue;
            bestFeatureIndex = -1;
            bestThreshold = 0;

            for (int i = 0; i < indices.Length; i++)
            {
                var x = xTrain.Column(indices[i]);
                var thresholds = IsRegression ? x : x.Distinct().ToArray();
                for (int j = 0; j < thresholds.Count(); j++)
                {
                    // Test if the split variance reduction or information gain
                    double performance = IsRegression ? VarianceReduction(x, yTrain, thresholds[j]) : InformationGain(x, yTrain, thresholds[j]);
                    // Keep track of the best value
                    if (performance > best)
                    {
                        best = performance;
                        bestFeatureIndex = indices[i];
                        bestThreshold = thresholds[j];
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
        private double VarianceReduction(double[] x, double[] y, double threshold)
        {
            // parent entropy
            var parentVariance = Statistics.PopulationVariance(y);

            // create children
            var leftIdxs = new List<int>();
            var rightIdxs = new List<int>();
            Split(x, threshold, out leftIdxs, out rightIdxs);

            if (leftIdxs.Count == 0 || rightIdxs.Count == 0)
                return double.MinValue;

            // calculate the weighted average variance of the children
            var n = (double)y.Length;
            var nLeft = (double)leftIdxs.Count;
            var nRight = (double)rightIdxs.Count;
            var yLeft = new double[leftIdxs.Count];
            for (int i = 0; i < leftIdxs.Count; i++)
                yLeft[i] = y[leftIdxs[i]];
            var yRight = new double[rightIdxs.Count];
            for (int i = 0; i < rightIdxs.Count; i++)
                yRight[i] = y[rightIdxs[i]];
            var varLeft = Statistics.PopulationVariance(yLeft);
            var varRight = Statistics.PopulationVariance(yRight);
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
        private double InformationGain(double[] x, double[] y, double threshold)
        {
            // parent entropy
            var parentE = Entropy(y);
            
            // create children
            var leftIdxs = new List<int>();
            var rightIdxs = new List<int>();
            Split(x, threshold, out leftIdxs, out rightIdxs);

            if (leftIdxs.Count == 0 || rightIdxs.Count == 0)
                return double.MinValue;

            // calculate the weighted average entropy of children
            var n = (double)y.Length;
            var nl = (double)leftIdxs.Count;
            var nr = (double)rightIdxs.Count;
            var yl = new double[leftIdxs.Count];
            for (int i = 0; i < leftIdxs.Count; i++)
                yl[i] = y[leftIdxs[i]];
            var yr = new double[rightIdxs.Count];
            for (int i = 0; i < rightIdxs.Count; i++)
                yr[i] = y[rightIdxs[i]];
            var el = Entropy(yl);
            var er = Entropy(yr);
            var childrenE = nl / n * el + nr / n * er;

            return parentE - childrenE;
        }

        /// <summary>
        /// Computes the entropy for vector of y-values.
        /// </summary>
        /// <param name="y">The column of y-values.</param>
        private double Entropy(double[] y)
        {
            if (IsRegression == true)
            {
                // use kernel density
                var kde = new KernelDensity(y);
                return Statistics.Entropy(y, kde.PDF);
            }
            else
            {
                // Use histogram
                return Statistics.Entropy(y, (x) => 
                {
                    double n = 0;
                    for (int i = 0; i < y.Length; i++)
                    {
                        if (x == y[i]){ n++; }
                    }
                    return n / y.Length;          
                });
            }
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
                return TraverseTree(x, node.Left);
            return TraverseTree(x, node.Right);
        }

        /// <summary>
        /// Returns the prediction from the Decision Tree.
        /// </summary>
        /// <param name="X">The 1D array of predictors.</param>
        public double[] Predict(double[] X)
        {
            return Predict(new Matrix(X));
        }

        /// <summary>
        /// Returns the prediction from the Decision Tree.
        /// </summary>
        /// <param name="X">The 2D array of predictors.</param>
        public double[] Predict(double[,] X)
        {
            return Predict(new Matrix(X));
        }

        /// <summary>
        /// Returns the prediction from the Decision Tree.
        /// </summary>
        /// <param name="X">The matrix of predictors.</param>
        public double[] Predict(Matrix X)
        {
            if (!IsTrained || X.NumberOfColumns != Dimensions) return null;
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
