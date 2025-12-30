# Machine Learning

[← Back to Index](../index.md)

The ***Numerics*** library provides machine learning algorithms for both supervised and unsupervised learning tasks. These implementations are designed for engineering and scientific applications including classification, regression, and clustering.

## Overview

**Supervised Learning:**
- Generalized Linear Models (GLM)
- Decision Trees
- Random Forests
- k-Nearest Neighbors (KNN)
- Naive Bayes

**Unsupervised Learning:**
- k-Means Clustering
- Gaussian Mixture Models (GMM)
- Jenks Natural Breaks

---

## Supervised Learning

### Generalized Linear Models (GLM)

GLMs extend linear regression to non-normal response distributions [[1]](#1):

```cs
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;

// Training data
double[,] X = {
    { 1, 2.5, 1.2 },  // Observation 1: [intercept, feature1, feature2]
    { 1, 3.1, 1.5 },
    { 1, 2.8, 1.1 },
    { 1, 3.5, 1.8 },
    { 1, 2.2, 0.9 }
};

double[] y = { 45.2, 52.3, 47.8, 58.1, 42.5 };  // Response variable

// Create GLM
var glm = new GeneralizedLinearModel(
    x: new Matrix(X),
    y: new Vector(y),
    family: GLMFamily.Normal,        // Distribution family
    linkFunction: LinkFunction.Identity  // Link function
);

// Set optimizer (optional)
glm.SetOptimizer(LocalMethod.NelderMead);

// Train model
glm.Train();

Console.WriteLine("GLM Results:");
Console.WriteLine($"Parameters: [{string.Join(", ", glm.Parameters.Select(p => p.ToString("F4")))}]");
Console.WriteLine($"Standard Errors: [{string.Join(", ", glm.ParameterStandardErrors.Select(se => se.ToString("F4")))}]");
Console.WriteLine($"p-values: [{string.Join(", ", glm.ParameterPValues.Select(p => p.ToString("F4")))}]");

// Model selection criteria
Console.WriteLine($"\nModel Selection:");
Console.WriteLine($"  AIC: {glm.AIC:F2}");
Console.WriteLine($"  BIC: {glm.BIC:F2}");
Console.WriteLine($"  Standard Error: {glm.StandardError:F4}");

// Make predictions
double[,] XNew = {
    { 1, 3.0, 1.4 },
    { 1, 2.6, 1.0 }
};

double[] predictions = glm.Predict(new Matrix(XNew));

Console.WriteLine($"\nPredictions:");
for (int i = 0; i < predictions.Length; i++)
{
    Console.WriteLine($"  X_new[{i}] → {predictions[i]:F2}");
}

// Prediction intervals (alpha = 0.1 for 90% interval)
double[,] intervals = glm.Predict(new Matrix(XNew), alpha: 0.1);

Console.WriteLine($"\n90% Prediction Intervals:");
for (int i = 0; i < XNew.GetLength(0); i++)
{
    Console.WriteLine($"  X_new[{i}]: [{intervals[i, 0]:F2}, {intervals[i, 1]:F2}]");
}
```

**Supported Families:**
- `GLMFamily.Normal` - Gaussian (linear regression)
- `GLMFamily.Binomial` - Binary outcomes (logistic regression)
- `GLMFamily.Poisson` - Count data
- `GLMFamily.Gamma` - Positive continuous data

**Link Functions:**
- `LinkFunction.Identity` - g(μ) = μ
- `LinkFunction.Log` - g(μ) = log(μ)
- `LinkFunction.Logit` - g(μ) = log(μ/(1-μ))
- `LinkFunction.Probit` - g(μ) = Φ⁻¹(μ)

### Decision Trees

Classification and regression trees [[2]](#2):

```cs
using Numerics.MachineLearning;

// Classification example
double[,] X = {
    { 5.1, 3.5, 1.4, 0.2 },  // Iris features
    { 4.9, 3.0, 1.4, 0.2 },
    { 7.0, 3.2, 4.7, 1.4 },
    { 6.4, 3.2, 4.5, 1.5 },
    { 6.3, 3.3, 6.0, 2.5 },
    { 5.8, 2.7, 5.1, 1.9 }
};

double[] y = { 0, 0, 1, 1, 2, 2 };  // Classes: Setosa(0), Versicolor(1), Virginica(2)

// Create decision tree
var tree = new DecisionTree(
    X: X,
    y: y,
    maxDepth: 5,              // Maximum tree depth
    minSamplesSplit: 2,       // Minimum samples to split node
    minSamplesLeaf: 1         // Minimum samples in leaf
);

// Train
tree.Train();

Console.WriteLine($"Decision Tree Trained: {tree.IsTrained}");

// Predict
double[] testSample = { 5.0, 3.0, 1.6, 0.2 };
double[] prediction = tree.Predict(testSample);

Console.WriteLine($"Prediction for test sample: Class {prediction[0]}");

// Predict multiple samples
double[,] testSamples = {
    { 5.0, 3.0, 1.6, 0.2 },
    { 6.0, 3.0, 4.5, 1.5 },
    { 6.5, 3.0, 5.5, 2.0 }
};

double[] predictions = tree.Predict(testSamples);

Console.WriteLine("\nBatch predictions:");
for (int i = 0; i < predictions.Length; i++)
{
    Console.WriteLine($"  Sample {i}: Class {predictions[i]}");
}
```

### Random Forests

Ensemble of decision trees for improved accuracy [[3]](#3):

```cs
using Numerics.MachineLearning;

double[,] X = {
    // Same Iris data as above
    { 5.1, 3.5, 1.4, 0.2 },
    { 4.9, 3.0, 1.4, 0.2 },
    { 7.0, 3.2, 4.7, 1.4 },
    { 6.4, 3.2, 4.5, 1.5 },
    { 6.3, 3.3, 6.0, 2.5 },
    { 5.8, 2.7, 5.1, 1.9 }
};

double[] y = { 0, 0, 1, 1, 2, 2 };

// Create random forest
var forest = new RandomForest(
    X: X,
    y: y,
    nTrees: 100,             // Number of trees
    maxDepth: 5,
    minSamplesSplit: 2,
    minSamplesLeaf: 1,
    maxFeatures: 2,          // Features per split
    bootstrap: true,         // Bootstrap sampling
    seed: 12345
);

// Train
forest.Train();

Console.WriteLine($"Random Forest Trained: {forest.IsTrained}");
Console.WriteLine($"Number of trees: {forest.NTrees}");

// Predict with confidence intervals
double[] testSample = { 5.0, 3.0, 1.6, 0.2 };
double[,] result = forest.Predict(testSample, alpha: 0.1);  // 90% CI

Console.WriteLine($"\nPrediction:");
Console.WriteLine($"  Predicted class: {result[0, 0]:F0}");
Console.WriteLine($"  90% CI: [{result[0, 1]:F2}, {result[0, 2]:F2}]");

// Batch prediction
double[,] testSamples = {
    { 5.0, 3.0, 1.6, 0.2 },
    { 6.0, 3.0, 4.5, 1.5 }
};

double[,] results = forest.Predict(testSamples, alpha: 0.1);

Console.WriteLine($"\nBatch predictions:");
for (int i = 0; i < testSamples.GetLength(0); i++)
{
    Console.WriteLine($"  Sample {i}: Class {results[i, 0]:F0}, " +
                     $"CI [{results[i, 1]:F2}, {results[i, 2]:F2}]");
}
```

**Advantages of Random Forests:**
- Reduces overfitting compared to single tree
- Provides prediction uncertainty
- Handles missing values well
- Works with mixed feature types

### k-Nearest Neighbors (KNN)

Non-parametric classification and regression [[4]](#4):

```cs
using Numerics.MachineLearning;

double[,] X = {
    { 1.0, 2.0 },
    { 1.5, 1.8 },
    { 5.0, 8.0 },
    { 8.0, 8.0 },
    { 1.0, 0.6 },
    { 9.0, 11.0 }
};

double[] y = { 0, 0, 1, 1, 0, 1 };  // Binary classification

// Create KNN classifier
var knn = new KNearestNeighbors(
    X: X,
    y: y,
    k: 3,                    // Number of neighbors
    weights: "uniform"       // "uniform" or "distance"
);

// KNN doesn't require explicit training
// Prediction happens at query time

// Predict
double[] testPoint = { 2.0, 3.0 };
double prediction = knn.Predict(testPoint);

Console.WriteLine($"KNN Prediction for [{testPoint[0]}, {testPoint[1]}]: Class {prediction}");

// Predict with probability estimates
double[,] probs = knn.PredictProba(testPoint);

Console.WriteLine($"Class probabilities:");
for (int i = 0; i < probs.GetLength(0); i++)
{
    Console.WriteLine($"  Class {i}: {probs[i, 0]:P1}");
}
```

**Distance Metrics:**
- Euclidean (default)
- Manhattan
- Minkowski

**Choosing k:**
- Small k: More sensitive to noise
- Large k: Smoother boundaries
- Rule of thumb: k = √n or use cross-validation

### Naive Bayes

Probabilistic classifier based on Bayes' theorem [[5]](#5):

```cs
using Numerics.MachineLearning;

// Text classification example (word counts)
double[,] X = {
    { 2, 1, 0, 1 },  // Document 1: word counts
    { 1, 1, 1, 0 },
    { 0, 3, 2, 1 },
    { 1, 0, 1, 2 }
};

double[] y = { 0, 0, 1, 1 };  // Classes: spam(1), ham(0)

// Create Naive Bayes
var nb = new NaiveBayes(X: X, y: y);

// Train
nb.Train();

Console.WriteLine("Naive Bayes trained");

// Predict
double[] testDoc = { 1, 2, 0, 1 };
double prediction = nb.Predict(testDoc);

Console.WriteLine($"Prediction: Class {prediction}");

// Class probabilities
double[] probabilities = nb.PredictProba(testDoc);

Console.WriteLine($"Class probabilities:");
Console.WriteLine($"  Class 0 (ham): {probabilities[0]:P1}");
Console.WriteLine($"  Class 1 (spam): {probabilities[1]:P1}");
```

**Assumptions:**
- Features are conditionally independent given class
- Works well despite violation of independence
- Fast training and prediction
- Good for text classification

---

## Unsupervised Learning

### k-Means Clustering

Partition data into k clusters [[6]](#6):

```cs
using Numerics.MachineLearning;

// 2D data points
double[,] X = {
    { 1.0, 2.0 },
    { 1.5, 1.8 },
    { 5.0, 8.0 },
    { 8.0, 8.0 },
    { 1.0, 0.6 },
    { 9.0, 11.0 },
    { 8.0, 2.0 },
    { 10.0, 2.0 },
    { 9.0, 3.0 }
};

// Create k-means with 3 clusters
var kmeans = new KMeans(X: X, k: 3);

// Configure
kmeans.MaxIterations = 100;
kmeans.Tolerance = 1e-4;
kmeans.Seed = 12345;

// Fit
kmeans.Fit();

Console.WriteLine($"k-Means Clustering (k={kmeans.K}):");
Console.WriteLine($"Converged: {kmeans.HasConverged}");
Console.WriteLine($"Iterations: {kmeans.Iterations}");
Console.WriteLine($"Inertia: {kmeans.Inertia:F2}");

// Cluster centers
Console.WriteLine($"\nCluster Centers:");
for (int i = 0; i < kmeans.K; i++)
{
    Console.WriteLine($"  Cluster {i}: [{kmeans.Means[i, 0]:F2}, {kmeans.Means[i, 1]:F2}]");
}

// Cluster labels
Console.WriteLine($"\nCluster Assignments:");
for (int i = 0; i < X.GetLength(0); i++)
{
    Console.WriteLine($"  Point [{X[i, 0]:F1}, {X[i, 1]:F1}] → Cluster {kmeans.Labels[i]}");
}

// Predict cluster for new point
double[] newPoint = { 2.0, 3.0 };
int cluster = kmeans.Predict(newPoint);

Console.WriteLine($"\nNew point [{newPoint[0]}, {newPoint[1]}] → Cluster {cluster}");

// Cluster sizes
var clusterSizes = kmeans.Labels.GroupBy(l => l).Select(g => g.Count()).ToArray();
Console.WriteLine($"\nCluster sizes: [{string.Join(", ", clusterSizes)}]");
```

**Choosing k:**
- Elbow method (plot inertia vs. k)
- Silhouette analysis
- Domain knowledge

**Initialization Methods:**
- Random selection
- k-means++ (default, better initialization)

### Gaussian Mixture Models (GMM)

Probabilistic clustering with soft assignments [[7]](#7):

```cs
using Numerics.MachineLearning;

double[,] X = {
    // Same data as k-means example
    { 1.0, 2.0 }, { 1.5, 1.8 }, { 5.0, 8.0 },
    { 8.0, 8.0 }, { 1.0, 0.6 }, { 9.0, 11.0 }
};

// Create GMM with 2 components
var gmm = new GaussianMixtureModel(
    X: X,
    nComponents: 2,
    covarianceType: "full"   // "full", "tied", "diag", "spherical"
);

// Configure
gmm.MaxIterations = 100;
gmm.Tolerance = 1e-3;
gmm.Seed = 12345;

// Fit using EM algorithm
gmm.Fit();

Console.WriteLine($"GMM Clustering ({gmm.NComponents} components):");
Console.WriteLine($"Converged: {gmm.HasConverged}");
Console.WriteLine($"Log-likelihood: {gmm.LogLikelihood:F2}");
Console.WriteLine($"BIC: {gmm.BIC:F2}");
Console.WriteLine($"AIC: {gmm.AIC:F2}");

// Component parameters
Console.WriteLine($"\nComponent Parameters:");
for (int i = 0; i < gmm.NComponents; i++)
{
    Console.WriteLine($"  Component {i}:");
    Console.WriteLine($"    Weight: {gmm.Weights[i]:F3}");
    Console.WriteLine($"    Mean: [{string.Join(", ", gmm.Means[i].Select(m => m.ToString("F2")))}]");
}

// Predict (hard assignment)
double[] newPoint = { 2.0, 3.0 };
int component = gmm.Predict(newPoint);

Console.WriteLine($"\nNew point [{newPoint[0]}, {newPoint[1]}] → Component {component}");

// Predict probabilities (soft assignment)
double[] probabilities = gmm.PredictProba(newPoint);

Console.WriteLine($"Component probabilities:");
for (int i = 0; i < probabilities.Length; i++)
{
    Console.WriteLine($"  Component {i}: {probabilities[i]:P1}");
}
```

**Advantages over k-Means:**
- Soft clustering (probabilistic assignments)
- Flexible cluster shapes (elliptical vs. spherical)
- Provides uncertainty quantification
- Can model overlapping clusters

### Jenks Natural Breaks

Optimal classification for univariate data [[8]](#8):

```cs
using Numerics.MachineLearning;

// Data values (e.g., elevation, rainfall, etc.)
double[] data = { 10, 12, 15, 18, 22, 25, 28, 35, 40, 45, 50, 55, 60, 70, 80 };

// Find natural breaks with 4 classes
int nClasses = 4;
var jenks = new JenksNaturalBreaks(data, nClasses);

jenks.Compute();

Console.WriteLine($"Jenks Natural Breaks ({nClasses} classes):");
Console.WriteLine($"Class breaks: [{string.Join(", ", jenks.Breaks.Select(b => b.ToString("F1")))}]");
Console.WriteLine($"Goodness of variance fit: {jenks.GoodnessOfVarianceFit:F4}");

// Classify data
int[] classes = jenks.Classify(data);

Console.WriteLine($"\nData classification:");
for (int i = 0; i < Math.Min(10, data.Length); i++)
{
    Console.WriteLine($"  Value {data[i]:F1} → Class {classes[i]}");
}

// Class statistics
for (int c = 0; c < nClasses; c++)
{
    var classData = data.Where((v, i) => classes[i] == c).ToArray();
    Console.WriteLine($"\nClass {c}:");
    Console.WriteLine($"  Range: [{classData.Min():F1}, {classData.Max():F1}]");
    Console.WriteLine($"  Count: {classData.Length}");
    Console.WriteLine($"  Mean: {classData.Average():F1}");
}
```

**Applications:**
- Choropleth map classification
- Data binning for visualization
- Natural grouping identification
- Minimizes within-class variance

---

## Practical Examples

### Example 1: Regression with GLM

```cs
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;

// Predict home prices
double[,] features = {
    { 1, 1500, 3, 20 },  // [intercept, sqft, bedrooms, age]
    { 1, 1800, 4, 15 },
    { 1, 1200, 2, 30 },
    { 1, 2000, 4, 10 },
    { 1, 1600, 3, 25 }
};

double[] prices = { 250000, 320000, 190000, 380000, 270000 };  // $

var glm = new GeneralizedLinearModel(
    new Matrix(features),
    new Vector(prices),
    GLMFamily.Normal,
    LinkFunction.Identity
);

glm.Train();

Console.WriteLine("Home Price Prediction Model:");
Console.WriteLine($"Coefficients:");
Console.WriteLine($"  Intercept: ${glm.Parameters[0]:F0}");
Console.WriteLine($"  Per sqft: ${glm.Parameters[1]:F2}");
Console.WriteLine($"  Per bedroom: ${glm.Parameters[2]:F0}");
Console.WriteLine($"  Per year age: ${glm.Parameters[3]:F0}");

// Predict new home
double[,] newHome = { { 1, 1700, 3, 12 } };
double predicted = glm.Predict(new Matrix(newHome))[0];
double[,] interval = glm.Predict(new Matrix(newHome), alpha: 0.1);

Console.WriteLine($"\nPrediction for 1700 sqft, 3BR, 12 years:");
Console.WriteLine($"  Predicted price: ${predicted:F0}");
Console.WriteLine($"  90% Interval: [${interval[0, 0]:F0}, ${interval[0, 1]:F0}]");
```

### Example 2: Classification Pipeline

```cs
// Iris classification
double[,] X_train = LoadIrisFeatures();  // Load training data
double[] y_train = LoadIrisLabels();
double[,] X_test = LoadIrisTestFeatures();
double[] y_test = LoadIrisTestLabels();

// Train random forest
var rf = new RandomForest(X_train, y_train, nTrees: 100, seed: 42);
rf.Train();

// Evaluate
double[,] predictions = rf.Predict(X_test);
int correct = 0;
for (int i = 0; i < y_test.Length; i++)
{
    if (predictions[i, 0] == y_test[i])
        correct++;
}

double accuracy = (double)correct / y_test.Length;

Console.WriteLine($"Random Forest Classification:");
Console.WriteLine($"  Accuracy: {accuracy:P1}");
Console.WriteLine($"  Correct: {correct}/{y_test.Length}");
```

### Example 3: Customer Segmentation

```cs
// Customer data: [annual_spending, visit_frequency, avg_basket_size]
double[,] customers = {
    { 1200, 24, 50 },   // Regular customer
    { 5000, 52, 95 },   // High-value customer
    { 300, 6, 45 },     // Occasional customer
    { 4800, 48, 100 },  // High-value customer
    { 800, 12, 65 },    // Regular customer
    { 250, 4, 55 },     // Occasional customer
    { 6000, 60, 105 }   // VIP customer
};

// Cluster into 3 segments
var kmeans = new KMeans(customers, k: 3);
kmeans.Fit();

Console.WriteLine("Customer Segmentation:");
for (int i = 0; i < 3; i++)
{
    var segment = Enumerable.Range(0, customers.GetLength(0))
                            .Where(j => kmeans.Labels[j] == i)
                            .ToArray();
    
    Console.WriteLine($"\nSegment {i} ({segment.Length} customers):");
    Console.WriteLine($"  Avg spending: ${segment.Average(j => customers[j, 0]):F0}");
    Console.WriteLine($"  Avg visits: {segment.Average(j => customers[j, 1]):F0}/year");
    Console.WriteLine($"  Avg basket: ${segment.Average(j => customers[j, 2]):F0}");
}
```

## Model Selection and Evaluation

### Cross-Validation

```cs
// Simple k-fold cross-validation
int k = 5;
int n = X.GetLength(0);
int foldSize = n / k;

double[] accuracies = new double[k];

for (int fold = 0; fold < k; fold++)
{
    // Split data into train/test
    var trainIndices = Enumerable.Range(0, n)
                                 .Where(i => i < fold * foldSize || i >= (fold + 1) * foldSize)
                                 .ToArray();
    
    var testIndices = Enumerable.Range(fold * foldSize, foldSize).ToArray();
    
    // Train and evaluate
    // ... (extract train/test sets, train model, compute accuracy)
    
    accuracies[fold] = ComputeAccuracy(testIndices);
}

Console.WriteLine($"Cross-Validation Results:");
Console.WriteLine($"  Mean accuracy: {accuracies.Average():P1}");
Console.WriteLine($"  Std dev: {Statistics.StandardDeviation(accuracies):F4}");
```

### Model Comparison

```cs
// Compare models on same dataset
var models = new[] {
    ("Decision Tree", new DecisionTree(X, y)),
    ("Random Forest", new RandomForest(X, y, nTrees: 50)),
    ("KNN (k=3)", new KNearestNeighbors(X, y, k: 3))
};

Console.WriteLine("Model Comparison:");
foreach (var (name, model) in models)
{
    model.Train();
    double accuracy = EvaluateModel(model, X_test, y_test);
    Console.WriteLine($"  {name}: {accuracy:P1}");
}
```

## Best Practices

### Supervised Learning
1. **Split data** - Use train/test split or cross-validation
2. **Normalize features** - Especially for distance-based methods (KNN)
3. **Handle imbalanced classes** - Use stratified sampling or class weights
4. **Tune hyperparameters** - Grid search or random search
5. **Validate assumptions** - Check residuals for GLM
6. **Ensemble methods** - Random Forests often outperform single trees

### Unsupervised Learning
1. **Scale features** - Clustering sensitive to feature scales
2. **Choose k carefully** - Use elbow method or silhouette scores
3. **Multiple runs** - k-Means sensitive to initialization
4. **Validate clusters** - Inspect cluster characteristics
5. **Consider GMM** - When clusters overlap or have different shapes

---

## References

<a id="1">[1]</a> Nelder, J. A., & Wedderburn, R. W. M. (1972). Generalized linear models. *Journal of the Royal Statistical Society: Series A*, 135(3), 370-384.

<a id="2">[2]</a> Breiman, L., Friedman, J., Stone, C. J., & Olshen, R. A. (1984). *Classification and Regression Trees*. CRC Press.

<a id="3">[3]</a> Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.

<a id="4">[4]</a> Cover, T., & Hart, P. (1967). Nearest neighbor pattern classification. *IEEE Transactions on Information Theory*, 13(1), 21-27.

<a id="5">[5]</a> Zhang, H. (2004). The optimality of naive Bayes. *AA*, 1(2), 3.

<a id="6">[6]</a> MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability*, 1(14), 281-297.

<a id="7">[7]</a> Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer.

<a id="8">[8]</a> Jenks, G. F. (1967). The data model concept in statistical mapping. *International Yearbook of Cartography*, 7, 186-190.

---

[← Back to Index](../index.md)
