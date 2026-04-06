# Machine Learning

[← Previous: Multivariate Distributions](../distributions/multivariate.md) | [Back to Index](../index.md) | [Next: Random Generation →](../sampling/random-generation.md)

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

GLMs extend linear regression to non-normal response distributions [[1]](#1). A GLM relates the expected value of the response variable $\mu = E(y)$ to a linear predictor $\eta = X\beta$ through a link function $g$:

```math
g(\mu) = \eta = X\beta = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \cdots + \beta_p x_p
```

The inverse link function maps from the linear predictor back to the mean of the response: $\mu = g^{-1}(\eta)$. Each link function corresponds to a distribution family:

| Link Function | $g(\mu)$ | $g^{-1}(\eta)$ | Distribution Family |
|---|---|---|---|
| Identity | $\mu$ | $\eta$ | Normal |
| Log | $\log(\mu)$ | $\exp(\eta)$ | Poisson |
| Logit | $\log\!\left(\frac{\mu}{1-\mu}\right)$ | $\frac{1}{1+\exp(-\eta)}$ | Binomial |
| Probit | $\Phi^{-1}(\mu)$ | $\Phi(\eta)$ | Binomial |
| CLogLog | $\log(-\log(1-\mu))$ | $1-\exp(-\exp(\eta))$ | Binomial |

**Parameter estimation.** Model parameters $\beta$ are estimated by maximum likelihood estimation (MLE), which maximizes the total log-likelihood over all $n$ observations:

```math
\hat{\beta} = \arg\max_{\beta} \sum_{i=1}^{n} \ell(\mu_i, y_i)
```

where the individual log-likelihood contribution $\ell(\mu_i, y_i)$ depends on the distribution family:

- **Normal family** (Identity link):

```math
\ell(\mu_i, y_i) = -\tfrac{1}{2}(y_i - \mu_i)^2
```

- **Poisson family** (Log link):

```math
\ell(\mu_i, y_i) = y_i \log(\mu_i) - \mu_i
```

- **Binomial family** (Logit, Probit, or CLogLog link):

```math
\ell(\mu_i, y_i) = y_i \log(\mu_i) + (1 - y_i)\log(1 - \mu_i)
```

**Model selection.** The Akaike Information Criterion (AIC) balances goodness of fit against model complexity:

```math
\text{AIC} = 2k - 2\ln(\hat{L})
```

where $k$ is the number of estimated parameters and $\hat{L}$ is the maximized likelihood. The corrected AIC (AICc) adjusts for small sample sizes, and the Bayesian Information Criterion (BIC) applies a stronger penalty for additional parameters.

```cs
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Functions;

// Training data (no intercept column needed — GLM adds it automatically when hasIntercept = true)
double[,] X = {
    { 2.5, 1.2 },  // Observation 1: [feature1, feature2]
    { 3.1, 1.5 },
    { 2.8, 1.1 },
    { 3.5, 1.8 },
    { 2.2, 0.9 }
};

double[] y = { 45.2, 52.3, 47.8, 58.1, 42.5 };  // Response variable

// Create GLM
var glm = new GeneralizedLinearModel(
    x: new Matrix(X),
    y: new Vector(y),
    hasIntercept: true,                        // Adds intercept column automatically
    linkType: LinkFunctionType.Identity        // Link function type
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
    { 3.0, 1.4 },
    { 2.6, 1.0 }
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
    Console.WriteLine($"  X_new[{i}]: [{intervals[i, 0]:F2}, {intervals[i, 2]:F2}]");
}
```

**Link Function Types** (`LinkFunctionType` in `Numerics.Functions`):
- `LinkFunctionType.Identity` - g(μ) = μ (Normal/Gaussian family)
- `LinkFunctionType.Log` - g(μ) = log(μ) (Poisson family)
- `LinkFunctionType.Logit` - g(μ) = log(μ/(1-μ)) (Binomial family)
- `LinkFunctionType.Probit` - g(μ) = Φ⁻¹(μ) (Binomial family, alternative)
- `LinkFunctionType.ComplementaryLogLog` - g(μ) = log(-log(1-μ)) (Asymmetric binary response)

### Decision Trees

Classification and regression trees [[2]](#2) recursively partition the feature space using binary splits of the form $x_j \leq t$ (left child) and $x_j > t$ (right child), where $x_j$ is a feature and $t$ is a threshold. At each internal node, the algorithm selects the feature and threshold that maximize a splitting criterion.

**Classification trees** use information gain based on Shannon entropy. The entropy of a node $S$ with $C$ classes is:

```math
H(S) = -\sum_{c=1}^{C} p_c \log(p_c)
```

where $p_c$ is the proportion of observations belonging to class $c$. The information gain from splitting node $S$ into children $S_L$ and $S_R$ is:

```math
\text{IG}(S, S_L, S_R) = H(S) - \frac{n_L}{n_S} H(S_L) - \frac{n_R}{n_S} H(S_R)
```

where $n_L$, $n_R$, and $n_S$ are the number of observations in the left child, right child, and parent node, respectively.

**Regression trees** use variance reduction. The population variance of a node $S$ is:

```math
\text{Var}(S) = \frac{1}{n_S} \sum_{i \in S} (y_i - \bar{y}_S)^2
```

The variance reduction from a split is:

```math
\text{VR}(S, S_L, S_R) = \text{Var}(S) - \frac{n_L}{n_S} \text{Var}(S_L) - \frac{n_R}{n_S} \text{Var}(S_R)
```

**Leaf node predictions.** For regression, the prediction is the mean of the training observations at the leaf: $\hat{y} = \bar{y}_{\text{leaf}}$. For classification, the prediction is the most frequent class among the leaf observations.

**Stopping criteria.** Tree growth stops when any of the following conditions is met: the maximum depth (`MaxDepth`, default 100) is reached, the number of observations at a node falls below `MinimumSplitSize` (default 2), or all observations at a node belong to the same class.

```cs
using Numerics.MachineLearning;

// Classification example (Iris-like data, requires at least 10 training samples)
double[,] X = {
    { 5.1, 3.5, 1.4, 0.2 },  // Iris features: sepal length, sepal width, petal length, petal width
    { 4.9, 3.0, 1.4, 0.2 },
    { 4.7, 3.2, 1.3, 0.2 },
    { 5.0, 3.6, 1.4, 0.2 },
    { 7.0, 3.2, 4.7, 1.4 },
    { 6.4, 3.2, 4.5, 1.5 },
    { 6.9, 3.1, 4.9, 1.5 },
    { 6.3, 3.3, 6.0, 2.5 },
    { 5.8, 2.7, 5.1, 1.9 },
    { 7.1, 3.0, 5.9, 2.1 },
    { 6.5, 3.0, 5.8, 2.2 }
};

double[] y = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2 };  // Classes: Setosa(0), Versicolor(1), Virginica(2)

// Create and train decision tree
var tree = new DecisionTree(X, y);
tree.MaxDepth = 5;  // Optional: limit tree depth (default: 100)
tree.Train();

Console.WriteLine($"Decision Tree Trained: {tree.IsTrained}");

// Predict single sample (pass as 2D array with 1 row for multi-feature input)
double[,] testSample = { { 5.0, 3.0, 1.6, 0.2 } };
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

Ensemble of decision trees for improved accuracy [[3]](#3) [[9]](#9). Random forests use bootstrap aggregation (bagging) to reduce variance and improve generalization. Given a training set of $n$ observations, the algorithm builds $B$ independent decision trees (default $B = 1000$), each trained on a bootstrap sample drawn with replacement from the original data.

**Bootstrap sampling.** For each tree $b = 1, \ldots, B$, draw $n$ samples with replacement from the original $n$ observations. On average, each bootstrap sample contains approximately $1 - 1/e \approx 63.2\%$ of the unique training observations, leaving the remainder as out-of-bag (OOB) observations that can be used for validation.

**Feature subsampling.** At each split within a tree, a random subset of features is considered as split candidates. The `Features` property controls the number of candidate features (default $d - 1$, where $d$ is the total number of features). This decorrelates the trees and further reduces variance.

**Aggregation.** The ensemble prediction for a new input $\mathbf{x}$ aggregates the predictions of all $B$ trees:

```math
\hat{f}_{\text{bag}}(\mathbf{x}) = \frac{1}{B} \sum_{b=1}^{B} \hat{f}_b(\mathbf{x})
```

For regression, predictions are the direct percentiles across all tree outputs. For classification, predictions are the floor of the percentiles. The `Predict` method returns lower, median, upper, and mean values across the ensemble, providing built-in prediction uncertainty quantification.

```cs
using Numerics.MachineLearning;

double[,] X = {
    // Same Iris-like data as above (requires at least 10 training samples)
    { 5.1, 3.5, 1.4, 0.2 },
    { 4.9, 3.0, 1.4, 0.2 },
    { 4.7, 3.2, 1.3, 0.2 },
    { 5.0, 3.6, 1.4, 0.2 },
    { 7.0, 3.2, 4.7, 1.4 },
    { 6.4, 3.2, 4.5, 1.5 },
    { 6.9, 3.1, 4.9, 1.5 },
    { 6.3, 3.3, 6.0, 2.5 },
    { 5.8, 2.7, 5.1, 1.9 },
    { 7.1, 3.0, 5.9, 2.1 },
    { 6.5, 3.0, 5.8, 2.2 }
};

double[] y = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2 };

// Create and train random forest
var forest = new RandomForest(X, y, seed: 12345);
forest.NumberOfTrees = 100;   // Default: 1000
forest.MaxDepth = 5;          // Default: 100
forest.Train();

Console.WriteLine($"Random Forest Trained: {forest.IsTrained}");
Console.WriteLine($"Number of trees: {forest.NumberOfTrees}");

// Predict with confidence intervals (pass as 2D array for multi-feature input)
// Predict returns double[,] with columns: lower(0), median(1), upper(2), mean(3)
double[,] testSample = { { 5.0, 3.0, 1.6, 0.2 } };
double[,] result = forest.Predict(testSample, alpha: 0.1);  // 90% CI

Console.WriteLine($"\nPrediction:");
Console.WriteLine($"  Predicted class (median): {result[0, 1]:F0}");
Console.WriteLine($"  Mean: {result[0, 3]:F2}");
Console.WriteLine($"  90% CI: [{result[0, 0]:F2}, {result[0, 2]:F2}]");

// Batch prediction
double[,] testSamples = {
    { 5.0, 3.0, 1.6, 0.2 },
    { 6.0, 3.0, 4.5, 1.5 }
};

double[,] results = forest.Predict(testSamples, alpha: 0.1);

Console.WriteLine($"\nBatch predictions:");
for (int i = 0; i < testSamples.GetLength(0); i++)
{
    Console.WriteLine($"  Sample {i}: Class {results[i, 1]:F0}, " +
                     $"CI [{results[i, 0]:F2}, {results[i, 2]:F2}]");
}
```

**Advantages of Random Forests:**
- Reduces overfitting compared to single tree
- Provides prediction uncertainty
- Handles missing values well
- Works with mixed feature types

### k-Nearest Neighbors (KNN)

Non-parametric classification and regression [[4]](#4). KNN is a lazy learning algorithm that stores the training data and defers computation until prediction time. Given a new input $\mathbf{x}$, the algorithm finds the $k$ nearest training observations using the Euclidean distance:

```math
d(\mathbf{x}, \mathbf{x}_i) = \sqrt{\sum_{j=1}^{p} (x_j - x_{ij})^2}
```

where $p$ is the number of features.

**Regression prediction.** For regression, the prediction is an inverse distance weighted average of the $k$ nearest neighbors:

```math
\hat{y} = \frac{\sum_{i=1}^{k} w_i \cdot y_i}{\sum_{i=1}^{k} w_i}, \quad w_i = \frac{1}{d(\mathbf{x}, \mathbf{x}_i)^2}
```

If the distance to a neighbor is exactly zero (i.e., the query point coincides with a training point), that neighbor receives a weight of $w_i = 1$ and all other weights are set accordingly.

**Classification prediction.** For classification, the prediction is determined by majority vote among the $k$ nearest neighbors: the class that appears most frequently is assigned to the new observation.

**Choosing $k$.** The choice of $k$ controls the bias-variance tradeoff. A small $k$ produces a flexible model sensitive to noise, while a large $k$ produces smoother decision boundaries at the cost of increased bias. A common rule of thumb is $k = \sqrt{n}$, where $n$ is the number of training observations, though cross-validation is preferred for rigorous selection.

**Curse of dimensionality.** As the number of features $p$ grows, Euclidean distances become increasingly uniform, degrading the discriminative power of nearest neighbor methods. Feature scaling (normalization) is critical because KNN is sensitive to the relative magnitudes of features.

```cs
using Numerics.MachineLearning;

double[,] X = {
    { 1.0, 2.0 },
    { 1.5, 1.8 },
    { 1.0, 0.6 },
    { 2.0, 1.5 },
    { 1.2, 1.0 },
    { 5.0, 8.0 },
    { 8.0, 8.0 },
    { 9.0, 11.0 },
    { 7.0, 9.0 },
    { 6.5, 7.5 }
};

double[] y = { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 };  // Binary classification

// Create KNN classifier (no explicit training needed — lazy learner)
var knn = new KNearestNeighbors(X, y, k: 3);

// Predict single sample
double[,] testPoint = { { 2.0, 3.0 } };
double[] prediction = knn.Predict(testPoint);
Console.WriteLine($"KNN Prediction for [{testPoint[0, 0]}, {testPoint[0, 1]}]: Class {prediction[0]}");

// Predict multiple samples
double[,] testPoints = { { 2.0, 3.0 }, { 7.0, 8.0 } };
double[] predictions = knn.Predict(testPoints);

Console.WriteLine($"Batch predictions:");
for (int i = 0; i < predictions.Length; i++)
{
    Console.WriteLine($"  Point {i}: Class {predictions[i]}");
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

Probabilistic classifier based on Bayes' theorem [[5]](#5). The Naive Bayes classifier applies Bayes' theorem with the assumption of conditional independence among features. Given a feature vector $\mathbf{x} = (x_1, x_2, \ldots, x_p)$, the posterior probability of class $c$ is:

```math
P(c \mid \mathbf{x}) \propto P(c) \cdot P(\mathbf{x} \mid c)
```

The conditional independence assumption simplifies the class-conditional likelihood to a product over individual features:

```math
P(\mathbf{x} \mid c) = \prod_{j=1}^{p} P(x_j \mid c)
```

This implementation uses Gaussian Naive Bayes, where each feature conditional on the class follows a normal distribution:

```math
P(x_j \mid c) = \mathcal{N}(x_j \mid \mu_{c,j},\, \sigma_{c,j}) = \frac{1}{\sqrt{2\pi}\,\sigma_{c,j}} \exp\!\left(-\frac{(x_j - \mu_{c,j})^2}{2\sigma_{c,j}^2}\right)
```

where $\mu_{c,j}$ and $\sigma_{c,j}$ are the mean and standard deviation of feature $j$ within class $c$, estimated from the training data.

**Prior probabilities.** The class prior $P(c)$ is estimated as the relative frequency of class $c$ in the training data:

```math
P(c) = \frac{n_c}{n}
```

where $n_c$ is the number of training observations in class $c$ and $n$ is the total number of observations.

**Classification rule.** The predicted class is determined by Maximum A Posteriori (MAP) estimation. Working in log-space for numerical stability:

```math
\hat{c} = \arg\max_{c} \left[ \log P(c) + \sum_{j=1}^{p} \log \mathcal{N}(x_j \mid \mu_{c,j},\, \sigma_{c,j}) \right]
```

```cs
using Numerics.MachineLearning;

// Text classification example (word counts, requires at least 10 training samples)
double[,] X = {
    { 2, 1, 0, 1 },  // Document 1: word counts
    { 1, 1, 1, 0 },
    { 3, 2, 0, 1 },
    { 2, 0, 1, 0 },
    { 1, 2, 0, 2 },
    { 0, 3, 2, 1 },
    { 1, 0, 1, 2 },
    { 0, 1, 3, 2 },
    { 1, 0, 2, 3 },
    { 0, 2, 2, 1 }
};

double[] y = { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 };  // Classes: spam(1), ham(0)

// Create and train Naive Bayes
var nb = new NaiveBayes(X, y);
nb.Train();

Console.WriteLine($"Naive Bayes trained: {nb.IsTrained}");

// Predict single sample
double[,] testDoc = { { 1, 2, 0, 1 } };
double[] prediction = nb.Predict(testDoc);
Console.WriteLine($"Prediction: Class {prediction[0]}");

// Predict multiple samples
double[,] testDocs = { { 1, 2, 0, 1 }, { 0, 3, 1, 0 } };
double[] predictions = nb.Predict(testDocs);
for (int i = 0; i < predictions.Length; i++)
    Console.WriteLine($"  Doc {i}: Class {predictions[i]}");
```

**Assumptions:**
- Features are conditionally independent given class
- Works well despite violation of independence
- Fast training and prediction
- Good for text classification

---

## Unsupervised Learning

### k-Means Clustering

Partition data into $k$ clusters [[6]](#6). The $k$-means algorithm seeks to minimize the within-cluster sum of squares (WCSS) objective function:

```math
J = \sum_{k=1}^{K} \sum_{\mathbf{x} \in C_k} \|\mathbf{x} - \boldsymbol{\mu}_k\|^2
```

where $C_k$ is the set of observations assigned to cluster $k$ and $\boldsymbol{\mu}_k$ is the centroid (mean) of cluster $k$.

**Lloyd's algorithm.** The algorithm alternates between two steps until convergence:

1. **E-step (assignment):** Assign each observation $\mathbf{x}_i$ to the nearest centroid using Euclidean distance:

```math
c_i = \arg\min_{k} \|\mathbf{x}_i - \boldsymbol{\mu}_k\|^2
```

2. **M-step (update):** Recompute each centroid as the mean of all observations assigned to its cluster:

```math
\boldsymbol{\mu}_k = \frac{1}{|C_k|} \sum_{\mathbf{x} \in C_k} \mathbf{x}
```

The algorithm converges when the cluster labels do not change between iterations, or when `MaxIterations` (default 1000) is reached.

**$k$-Means++ initialization** [[11]](#11). The default initialization method selects initial centroids that are well-separated. The first centroid is chosen uniformly at random. Each subsequent centroid is selected with probability proportional to $D(\mathbf{x})^2$, where $D(\mathbf{x})$ is the Euclidean distance from $\mathbf{x}$ to the nearest already-chosen centroid. This initialization significantly reduces the chance of poor convergence compared to random initialization.

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
var kmeans = new KMeans(X, k: 3);
kmeans.MaxIterations = 100;

// Train (use seed for reproducibility, k-means++ initialization by default)
kmeans.Train(seed: 12345);

Console.WriteLine($"k-Means Clustering (k={kmeans.K}):");
Console.WriteLine($"Iterations: {kmeans.Iterations}");

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

Probabilistic clustering with soft assignments [[7]](#7) [[10]](#10). A Gaussian Mixture Model represents the data as a mixture of $K$ multivariate Gaussian components. The probability density function of the mixture is:

```math
p(\mathbf{x}) = \sum_{k=1}^{K} \pi_k \cdot \mathcal{N}(\mathbf{x} \mid \boldsymbol{\mu}_k, \boldsymbol{\Sigma}_k)
```

where $\pi_k$ are the mixing weights (with $\sum_k \pi_k = 1$ and $\pi_k \geq 0$), $\boldsymbol{\mu}_k$ is the mean vector, and $\boldsymbol{\Sigma}_k$ is the covariance matrix of component $k$.

**Expectation-Maximization (EM) algorithm.** The parameters are estimated by maximizing the log-likelihood using the EM algorithm, initialized from a $k$-means solution with equal weights and near-identity covariance matrices.

The **E-step** computes the responsibility of each component $k$ for each observation $\mathbf{x}_i$:

```math
r_{ik} = \frac{\pi_k \cdot \mathcal{N}(\mathbf{x}_i \mid \boldsymbol{\mu}_k, \boldsymbol{\Sigma}_k)}{\sum_{j=1}^{K} \pi_j \cdot \mathcal{N}(\mathbf{x}_i \mid \boldsymbol{\mu}_j, \boldsymbol{\Sigma}_j)}
```

The implementation uses the log-sum-exp trick with Cholesky decomposition for numerical stability.

The **M-step** updates the parameters using the responsibilities:

```math
\pi_k = \frac{1}{n}\sum_{i=1}^{n} r_{ik}, \qquad \boldsymbol{\mu}_k = \frac{\sum_{i=1}^{n} r_{ik}\,\mathbf{x}_i}{\sum_{i=1}^{n} r_{ik}}, \qquad \boldsymbol{\Sigma}_k = \frac{\sum_{i=1}^{n} r_{ik}\,(\mathbf{x}_i - \boldsymbol{\mu}_k)(\mathbf{x}_i - \boldsymbol{\mu}_k)^\top}{\sum_{i=1}^{n} r_{ik}}
```

**Convergence.** The algorithm iterates until the relative change in log-likelihood falls below `Tolerance` (default $10^{-8}$) or `MaxIterations` (default 1000) is reached. The total log-likelihood is:

```math
\mathcal{L} = \sum_{i=1}^{n} \log\!\left(\sum_{k=1}^{K} \pi_k \cdot \mathcal{N}(\mathbf{x}_i \mid \boldsymbol{\mu}_k, \boldsymbol{\Sigma}_k)\right)
```

**Connection to $k$-means.** GMM generalizes $k$-means: when all covariance matrices are constrained to $\sigma^2 \mathbf{I}$ and $\sigma^2 \to 0$, the soft responsibilities $r_{ik}$ converge to hard 0/1 assignments, recovering the $k$-means solution.

```cs
using Numerics.MachineLearning;

double[,] X = {
    // Same data as k-means example
    { 1.0, 2.0 }, { 1.5, 1.8 }, { 5.0, 8.0 },
    { 8.0, 8.0 }, { 1.0, 0.6 }, { 9.0, 11.0 }
};

// Create GMM with 2 components
var gmm = new GaussianMixtureModel(X, k: 2);
gmm.MaxIterations = 100;
gmm.Tolerance = 1e-3;

// Train using EM algorithm
gmm.Train(seed: 12345);

Console.WriteLine($"GMM Clustering ({gmm.K} components):");
Console.WriteLine($"Iterations: {gmm.Iterations}");
Console.WriteLine($"Log-likelihood: {gmm.LogLikelihood:F2}");

// Component parameters
Console.WriteLine($"\nComponent Parameters:");
for (int i = 0; i < gmm.K; i++)
{
    Console.WriteLine($"  Component {i}:");
    Console.WriteLine($"    Weight: {gmm.Weights[i]:F3}");
    int dims = X.GetLength(1);
    Console.Write($"    Mean: [");
    for (int d = 0; d < dims; d++)
        Console.Write($"{gmm.Means[i, d]:F2}{(d < dims - 1 ? ", " : "")}");
    Console.WriteLine("]");
}

// Cluster labels
Console.WriteLine($"\nCluster Assignments:");
for (int i = 0; i < X.GetLength(0); i++)
{
    Console.WriteLine($"  Point [{X[i, 0]:F1}, {X[i, 1]:F1}] → Component {gmm.Labels[i]}");
}
```

**Advantages over k-Means:**
- Soft clustering (probabilistic assignments)
- Flexible cluster shapes (elliptical vs. spherical)
- Provides uncertainty quantification
- Can model overlapping clusters

### Jenks Natural Breaks

Optimal classification for univariate data [[8]](#8). The Jenks Natural Breaks algorithm (also known as the Fisher-Jenks algorithm) partitions a sorted univariate dataset into $k$ classes by minimizing the within-class sum of squared deviations while maximizing the between-class separation. The quality of a classification is measured by the Goodness of Variance Fit (GVF):

```math
\text{GVF} = \frac{\text{SDAM} - \text{SDCM}}{\text{SDAM}}
```

where SDAM is the sum of squared deviations from the array (overall) mean:

```math
\text{SDAM} = \sum_{i=1}^{n} (x_i - \bar{x})^2
```

and SDCM is the sum of squared deviations from the class means:

```math
\text{SDCM} = \sum_{c=1}^{k} \sum_{x \in C_c} (x - \bar{x}_c)^2
```

A GVF of 1.0 indicates a perfect fit (zero within-class variance), while values closer to 0 indicate poor classification. The algorithm uses dynamic programming to efficiently find the optimal break points that maximize the GVF, avoiding the exponential cost of exhaustive search over all possible partitions.

```cs
using Numerics.MachineLearning;

// Data values (e.g., elevation, rainfall, etc.)
double[] data = { 10, 12, 15, 18, 22, 25, 28, 35, 40, 45, 50, 55, 60, 70, 80 };

// Find natural breaks with 4 classes (computation happens in constructor)
int nClasses = 4;
var jenks = new JenksNaturalBreaks(data, nClasses);

Console.WriteLine($"Jenks Natural Breaks ({nClasses} classes):");
Console.WriteLine($"Break points: [{string.Join(", ", jenks.Breaks.Select(b => b.ToString("F1")))}]");
Console.WriteLine($"Goodness of variance fit: {jenks.GoodnessOfVarianceFit:F4}");

// Access cluster details
Console.WriteLine($"\nCluster details:");
for (int c = 0; c < jenks.Clusters.Length; c++)
{
    var cluster = jenks.Clusters[c];
    Console.WriteLine($"  Cluster {c}: {cluster.Count} values");
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
using Numerics.Functions;

// Predict home prices (no intercept column — GLM adds it automatically)
double[,] features = {
    { 1500, 3, 20 },  // [sqft, bedrooms, age]
    { 1800, 4, 15 },
    { 1200, 2, 30 },
    { 2000, 4, 10 },
    { 1600, 3, 25 }
};

double[] prices = { 250000, 320000, 190000, 380000, 270000 };  // $

var glm = new GeneralizedLinearModel(
    new Matrix(features),
    new Vector(prices),
    hasIntercept: true,
    linkType: LinkFunctionType.Identity
);

glm.Train();

Console.WriteLine("Home Price Prediction Model:");
Console.WriteLine($"Coefficients:");
Console.WriteLine($"  Intercept: ${glm.Parameters[0]:F0}");
Console.WriteLine($"  Per sqft: ${glm.Parameters[1]:F2}");
Console.WriteLine($"  Per bedroom: ${glm.Parameters[2]:F0}");
Console.WriteLine($"  Per year age: ${glm.Parameters[3]:F0}");

// Predict new home
double[,] newHome = { { 1700, 3, 12 } };
double predicted = glm.Predict(new Matrix(newHome))[0];
// Predict returns columns: lower(0), mean(1), upper(2)
double[,] interval = glm.Predict(new Matrix(newHome), alpha: 0.1);

Console.WriteLine($"\nPrediction for 1700 sqft, 3BR, 12 years:");
Console.WriteLine($"  Predicted price: ${predicted:F0}");
Console.WriteLine($"  90% Interval: [${interval[0, 0]:F0}, ${interval[0, 2]:F0}]");
```

### Example 2: Classification Pipeline

```cs
// Sample binary classification data (requires at least 10 training samples)
double[,] X_train = {
    { 2.5, 3.2 }, { 3.1, 2.8 }, { 2.8, 3.5 }, { 3.3, 2.9 }, { 2.6, 3.1 },  // Class 0
    { 6.2, 5.8 }, { 5.9, 6.1 }, { 6.5, 5.5 }, { 5.8, 6.3 }, { 6.1, 5.7 }   // Class 1
};
double[] y_train = { 0, 0, 0, 0, 0, 1, 1, 1, 1, 1 };

double[,] X_test = {
    { 2.9, 3.0 }, { 6.0, 5.9 }  // One from each class
};
double[] y_test = { 0, 1 };

// Train random forest (no nTrees constructor parameter — set via property)
var rf = new RandomForest(X_train, y_train, seed: 42);
rf.NumberOfTrees = 100;
rf.Train();

// Evaluate — Predict returns double[,] with columns: lower(0), median(1), upper(2), mean(3)
double[,] predictions = rf.Predict(X_test);
int correct = 0;
for (int i = 0; i < y_test.Length; i++)
{
    if (predictions[i, 1] == y_test[i])  // Use median (column 1) as predicted class
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
kmeans.Train(seed: 42);

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

### Example 4: Flood Damage Prediction with Logistic GLM

Predicting whether flood damage occurs based on hydraulic variables using data from [`example-data/flood-damage-glm.csv`](example-data/flood-damage-glm.csv):

```cs
using System.IO;
using System.Linq;
using Numerics.MachineLearning;
using Numerics.Functions;
using Numerics.Mathematics.LinearAlgebra;

// Load CSV data (skip comment lines starting with #)
string[] lines = File.ReadAllLines("example-data/flood-damage-glm.csv");
var dataLines = lines
    .Where(line => !line.StartsWith("#") && !string.IsNullOrWhiteSpace(line))
    .Skip(1) // Skip header
    .ToArray();

int n = dataLines.Length;
double[,] features = new double[n, 3];
double[] damage = new double[n];

for (int i = 0; i < n; i++)
{
    var parts = dataLines[i].Split(',');
    features[i, 0] = double.Parse(parts[0]); // FloodStage_ft
    features[i, 1] = double.Parse(parts[1]); // Duration_hr
    features[i, 2] = double.Parse(parts[2]); // Velocity_fps
    damage[i] = double.Parse(parts[3]);       // DamageOccurred (0/1)
}

// Fit logistic regression (GLM with Logit link)
var glm = new GeneralizedLinearModel(
    new Matrix(features),
    new Vector(damage),
    hasIntercept: true,
    linkType: LinkFunctionType.Logit
);

glm.Train();

// Print R-style summary
Console.WriteLine("Flood Damage Logistic Regression");
Console.WriteLine("=" + new string('=', 50));
foreach (var line in glm.Summary())
{
    Console.WriteLine(line);
}

// Model fit statistics
Console.WriteLine($"\nAIC: {glm.AIC:F2}");
Console.WriteLine($"AICc: {glm.AICc:F2}");
Console.WriteLine($"BIC: {glm.BIC:F2}");

// Predict damage probability for a new flood event
double[,] newEvent = { { 14.0, 10, 3.5 } };  // Stage=14ft, Duration=10hr, Velocity=3.5fps
double[] prob = glm.Predict(new Matrix(newEvent));
Console.WriteLine($"\nPredicted damage probability: {prob[0]:P1}");

// Prediction with 90% confidence interval
double[,] interval = glm.Predict(new Matrix(newEvent), alpha: 0.1);
Console.WriteLine($"90% Prediction interval: [{interval[0, 0]:P1}, {interval[0, 2]:P1}]");
Console.WriteLine($"  Mean prediction: {interval[0, 1]:P1}");

// Compare link functions using AIC
var linkTypes = new[] {
    LinkFunctionType.Logit,
    LinkFunctionType.Probit,
    LinkFunctionType.ComplementaryLogLog
};

Console.WriteLine("\nLink Function Comparison:");
Console.WriteLine("  Link              |     AIC |     BIC");
Console.WriteLine("  ------------------|---------|--------");

foreach (var link in linkTypes)
{
    var model = new GeneralizedLinearModel(
        new Matrix(features),
        new Vector(damage),
        hasIntercept: true,
        linkType: link
    );
    model.Train();
    Console.WriteLine($"  {link,-18} | {model.AIC,7:F2} | {model.BIC,7:F2}");
}
```

## Model Selection and Evaluation

### Cross-Validation

```cs
using Numerics.Data.Statistics;

// Simple k-fold cross-validation
int k = 5;
int n = X.GetLength(0);
int dims = X.GetLength(1);
int foldSize = n / k;

double[] accuracies = new double[k];

for (int fold = 0; fold < k; fold++)
{
    // Split data into train/test indices
    var testIndices = Enumerable.Range(fold * foldSize, foldSize).ToArray();
    var trainIndices = Enumerable.Range(0, n)
                                 .Where(i => !testIndices.Contains(i))
                                 .ToArray();

    // Extract train/test data by copying rows
    double[,] X_trainFold = new double[trainIndices.Length, dims];
    double[] y_trainFold = new double[trainIndices.Length];
    for (int i = 0; i < trainIndices.Length; i++)
    {
        for (int d = 0; d < dims; d++)
            X_trainFold[i, d] = X[trainIndices[i], d];
        y_trainFold[i] = y[trainIndices[i]];
    }

    double[,] X_testFold = new double[testIndices.Length, dims];
    double[] y_testFold = new double[testIndices.Length];
    for (int i = 0; i < testIndices.Length; i++)
    {
        for (int d = 0; d < dims; d++)
            X_testFold[i, d] = X[testIndices[i], d];
        y_testFold[i] = y[testIndices[i]];
    }

    // Train and evaluate
    var model = new DecisionTree(X_trainFold, y_trainFold);
    model.Train();
    var predictions = model.Predict(X_testFold);

    int correct = 0;
    for (int i = 0; i < testIndices.Length; i++)
        if (predictions[i] == y_testFold[i]) correct++;

    accuracies[fold] = (double)correct / testIndices.Length;
}

Console.WriteLine($"Cross-Validation Results:");
Console.WriteLine($"  Mean accuracy: {accuracies.Average():P1}");
Console.WriteLine($"  Std dev: {Statistics.StandardDeviation(accuracies):F4}");
```

### Model Comparison

```cs
// Compare different classifiers on the same dataset
double[,] X = { /* training features (at least 10 rows) */ };
double[] y = { /* training labels */ };
double[,] X_test = { /* test features */ };
double[] y_test = { /* test labels */ };

// Train models (KNN is lazy — no Train() call needed)
var decisionTree = new DecisionTree(X, y);
decisionTree.Train();

var randomForest = new RandomForest(X, y, seed: 42);
randomForest.NumberOfTrees = 50;
randomForest.Train();

var knn = new KNearestNeighbors(X, y, k: 3);

// Evaluate each model
// Note: DecisionTree and KNN return double[], RandomForest returns double[,]
Console.WriteLine("Model Comparison:");

// DecisionTree returns double[]
double[] dtPredictions = decisionTree.Predict(X_test);
int dtCorrect = 0;
for (int i = 0; i < y_test.Length; i++)
    if (dtPredictions[i] == y_test[i]) dtCorrect++;
Console.WriteLine($"  Decision Tree: {(double)dtCorrect / y_test.Length:P1}");

// RandomForest returns double[,] with columns: lower(0), median(1), upper(2), mean(3)
double[,] rfPredictions = randomForest.Predict(X_test);
int rfCorrect = 0;
for (int i = 0; i < y_test.Length; i++)
    if (rfPredictions[i, 1] == y_test[i]) rfCorrect++;  // Use median (column 1)
Console.WriteLine($"  Random Forest: {(double)rfCorrect / y_test.Length:P1}");

// KNN returns double[]
double[] knnPredictions = knn.Predict(X_test);
int knnCorrect = 0;
for (int i = 0; i < y_test.Length; i++)
    if (knnPredictions[i] == y_test[i]) knnCorrect++;
Console.WriteLine($"  KNN (k=3): {(double)knnCorrect / y_test.Length:P1}");
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

<a id="2">[2]</a> Breiman, L., Friedman, J. H., Olshen, R. A., & Stone, C. J. (1984). *Classification and Regression Trees*. Wadsworth.

<a id="3">[3]</a> Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.

<a id="4">[4]</a> Cover, T., & Hart, P. (1967). Nearest neighbor pattern classification. *IEEE Transactions on Information Theory*, 13(1), 21-27.

<a id="5">[5]</a> Zhang, H. (2004). The optimality of naive Bayes. *Proceedings of the Seventeenth International FLAIRS Conference*, 562-567.

<a id="6">[6]</a> MacQueen, J. (1967). Some methods for classification and analysis of multivariate observations. *Proceedings of the Fifth Berkeley Symposium on Mathematical Statistics and Probability*, 1, 281-297.

<a id="7">[7]</a> Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer.

<a id="8">[8]</a> Jenks, G. F. (1967). The data model concept in statistical mapping. *International Yearbook of Cartography*, 7, 186-190.

<a id="9">[9]</a> Hastie, T., Tibshirani, R., & Friedman, J. (2009). *The Elements of Statistical Learning* (2nd ed.). Springer.

<a id="10">[10]</a> Dempster, A. P., Laird, N. M., & Rubin, D. B. (1977). Maximum likelihood from incomplete data via the EM algorithm. *Journal of the Royal Statistical Society: Series B*, 39(1), 1-38.

<a id="11">[11]</a> Arthur, D., & Vassilvitskii, S. (2007). k-means++: The advantages of careful seeding. *Proceedings of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms*, 1027-1035.

---

[← Previous: Multivariate Distributions](../distributions/multivariate.md) | [Back to Index](../index.md) | [Next: Random Generation →](../sampling/random-generation.md)
