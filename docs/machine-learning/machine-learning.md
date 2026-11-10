# Machine Learning

[← Previous: Univariate Functions](../functions/index.md) | [Back to Index](../index.md) | [Next: Random Generation →](../sampling/random-generation.md)

The ***Numerics*** library provides machine learning algorithms for both supervised and unsupervised learning tasks. These implementations are designed for engineering and scientific applications including classification, regression, and clustering.

## Overview

**Supervised Learning:**
- Generalized Linear Models (GLM)
- Decision Trees
- Random Forests
- k-Nearest Neighbors (KNN)
- Naive Bayes
- Gaussian Process Regression (GPR)
- Kriging
- Elastic Net Regression (lasso, ridge and elastic net paths with cross-validation)

**Unsupervised Learning:**
- k-Means Clustering
- Gaussian Mixture Models (GMM)
- Jenks Natural Breaks
- Principal Component Analysis (PCA)

**Parallel execution.** k-Means, Gaussian mixture models, KNN and Naive Bayes split large passes (cluster assignment, the EM steps, query batches and prediction rows) across processors and run small problems serially. Each row or component is computed independently and every sum is accumulated in a fixed order, so results are identical for any processor count.

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

**Fitting.** The estimate is computed by iteratively reweighted least squares (Fisher scoring), the algorithm used by R's `glm`, which also centers the optimizer's starting point and bounds; the optimizer selected with `SetOptimizer` then refines it, and the fit with the higher likelihood is kept. If Fisher scoring does not converge, as under complete separation in a binary model, the optimizer starts from family-specific values inside fixed bounds.

**Standard errors.** The parameter covariance is the inverse Fisher information scaled by the dispersion:

```math
\operatorname{Cov}(\hat{\beta}) = \phi\,(X^\top W X)^{-1}, \qquad W_{ii} = \frac{(d\mu_i/d\eta_i)^2}{V(\mu_i)}
```

with variance function $V(\mu) = 1$, $\mu$ and $\mu(1-\mu)$ for the Normal, Poisson and Binomial families, and dispersion $\phi$ equal to the residual variance for the Normal family and 1 otherwise. Setting `UseRobustSE = true` gives the HC0 sandwich estimator instead.

**Confidence intervals for the mean.** `Predict(x, alpha)` forms $\hat{\eta} \pm z_{1-\alpha/2}\,\text{se}(\hat{\eta})$ on the link scale, with $\text{se}(\hat{\eta})^2 = x^\top \operatorname{Cov}(\hat{\beta})\, x$, and maps both ends through the inverse link, so the interval stays inside the response range. It describes the uncertainty of the mean response, not the spread of new observations.

**Model selection.** The Akaike Information Criterion (AIC) balances goodness of fit against model complexity:

```math
\text{AIC} = 2k - 2\ln(\hat{L})
```

where $k$ is the number of regression coefficients and $\hat{L}$ is the maximized likelihood. For the Normal family the likelihood is evaluated at the maximum likelihood variance $\hat{\sigma}^2 = \text{SSE}/n$, which matches statsmodels GLM; R's `glm` also counts $\sigma^2$ as a parameter, so its Gaussian AIC is larger by 2. The residual standard error `StandardError` keeps the divisor $n - k$. With `hasIntercept = true`, `ParameterNames` is `"Intercept"` followed by the predictor matrix header names (or β1, β2, ... when the matrix has no header). The corrected AIC (AICc) adjusts for small sample sizes, and the Bayesian Information Criterion (BIC) applies a stronger penalty for additional parameters.

```cs
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Functions;

// Training data (no intercept column needed; GLM adds it automatically when hasIntercept = true)
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

// Confidence intervals for the mean response (alpha = 0.1 for 90% intervals)
double[,] intervals = glm.Predict(new Matrix(XNew), alpha: 0.1);

Console.WriteLine($"\n90% Confidence Intervals for the Mean:");
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
tree.IsRegression = false;  // Classification (default: regression)
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
forest.IsRegression = false;  // Classification (default: regression)
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
- Accepts missing (NaN) predictors: at each split such rows go to the right child
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

If one or more of the $k$ neighbors are at exactly zero distance (the query coincides with training points), the prediction is the plain average of those neighbors' responses and the other neighbors receive no weight, the rule scikit-learn applies for distance weighting.

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

// Create KNN classifier (no explicit training needed; KNN is a lazy learner)
var knn = new KNearestNeighbors(X, y, k: 3);
knn.IsRegression = false;  // Classification (default: regression)

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

**Distance metric:** Euclidean. Scale the features first when they have different units.

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

### Gaussian Process Regression (GPR)

Gaussian process regression [[17]](#17) is a nonparametric Bayesian regression method that gives a predictive mean and a predictive variance at every point. It is well suited to surrogate models of expensive simulations, where the data are a few hundred to a few thousand runs and an honest measure of uncertainty is needed. The response is modeled as

```math
y = m(\mathbf{x}) + f(\mathbf{x}) + \varepsilon, \qquad f \sim \mathcal{GP}\left(0, \sigma_f^2 \rho(r)\right), \qquad \varepsilon \sim N(0, \sigma_n^2)
```

where $\sigma_f^2$ is the signal variance, $\sigma_n^2$ is the noise (nugget) variance, and $\rho$ is a stationary correlation of the scaled distance

```math
r = \sqrt{\sum_{j=1}^{d} \frac{(x_j - x'_j)^2}{\ell_j^2}}
```

With automatic relevance determination (ARD, the default) each feature has its own length scale $\ell_j$, so features that do not affect the response get long length scales. Set `UseARD = false` to share one length scale.

**Covariance functions** (`GaussianProcessKernel`):

| Kernel | $\rho(r)$ | Smoothness |
|---|---|---|
| `SquaredExponential` | $\exp(-r^2/2)$ | Infinitely differentiable |
| `Matern52` (default) | $(1 + \sqrt{5} r + 5r^2/3)\exp(-\sqrt{5} r)$ | Twice differentiable |
| `Matern32` | $(1 + \sqrt{3} r)\exp(-\sqrt{3} r)$ | Once differentiable |
| `Exponential` | $\exp(-r)$ | Continuous |

**Mean functions** (`GaussianProcessMean`): `Zero`, `Constant` (default) or `Linear` in the features. The mean coefficients $\boldsymbol\beta$ are estimated by generalized least squares for the current covariance, as in universal kriging, and the predictive variance includes their estimation uncertainty.

**Prediction.** With $K = \sigma_f^2 R + \sigma_n^2 I$ the covariance of the training responses, $F$ the mean basis, $\mathbf{k}_*$ the covariances between a new point and the training points, and $\boldsymbol\alpha = K^{-1}(\mathbf{y} - F\boldsymbol\beta)$:

```math
\mu_* = \mathbf{f}_*^\top \boldsymbol\beta + \mathbf{k}_*^\top \boldsymbol\alpha, \qquad
\sigma_*^2 = \sigma_f^2 - \mathbf{k}_*^\top K^{-1} \mathbf{k}_* + \mathbf{u}^\top (F^\top K^{-1} F)^{-1} \mathbf{u}, \quad \mathbf{u} = \mathbf{f}_* - F^\top K^{-1} \mathbf{k}_*
```

`PredictWithVariance(x)` returns the variance of the latent function $f$; `PredictWithVariance(x, includeNoise: true)` adds $\sigma_n^2$ to give the variance of a new observation.

**Hyperparameter estimation.** `Train()` maximizes the log marginal likelihood

```math
\ln p(\mathbf{y} \mid \boldsymbol\theta) = -\tfrac{1}{2} (\mathbf{y} - F\boldsymbol\beta)^\top K^{-1} (\mathbf{y} - F\boldsymbol\beta) - \tfrac{1}{2} \ln |K| - \tfrac{n}{2} \ln 2\pi
```

with respect to the hyperparameters. The covariance is written $K = \sigma^2(R + gI)$ with the noise-to-signal ratio $g = \sigma_n^2/\sigma_f^2$, and $\sigma^2$ is profiled out at its closed-form maximum $\hat\sigma^2 = (\mathbf{y} - F\boldsymbol\beta)^\top (R + gI)^{-1} (\mathbf{y} - F\boldsymbol\beta)/n$, as in DiceKriging [[19]](#19). `LBFGSB` maximizes the resulting concentrated likelihood over $(\ln \ell_1, \ldots, \ln \ell_d, \ln g)$ with the analytic gradient $\partial \ln p / \partial \theta_k = \tfrac{1}{2}\operatorname{tr}\left((\boldsymbol\alpha\boldsymbol\alpha^\top - K^{-1})\,\partial K/\partial\theta_k\right)$ [[17]](#17), stopping at a projected gradient of `Tolerance` (default $10^{-5}$) or a relative reduction of $2.2 \times 10^{-9}$, the L-BFGS-B defaults that scikit-learn also uses; see the [optimization guide](../mathematics/optimization.md). The likelihood usually has several local maxima that differ in which features are judged irrelevant (their length scales go to the upper bound). Two data-based starts are therefore always run, $\ell_j = r_j/2$ and $\ell_j = r_j$ with $g = 0.01$, where $r_j$ is the range of feature $j$, and the better result is kept; `NumberOfRestarts` adds starts drawn uniformly from a box of plausible values by a Mersenne Twister seeded with `PRNGSeed`. `Train(noiseVariance)` holds the noise fixed (zero gives an interpolating model), and `SetHyperparameters` conditions the model on given values without optimization. `MeanCoefficients` holds the generalized least squares coefficients $\boldsymbol\beta$ of the mean function: empty for a zero mean, $[\beta_0]$ for a constant mean and $[\beta_0, \beta_1, \ldots, \beta_d]$ for a linear mean. `EvaluateLogMarginalLikelihood(theta, gradient)` evaluates the log marginal likelihood, and optionally its gradient, at the log hyperparameters $(\ln\sigma_f^2, \ln\ell_1, \ldots, \ln\ell_q, \ln\sigma_n^2)$ without changing the model, with $q = d$ under ARD and $q = 1$ otherwise. It returns NaN when the covariance matrix is not numerically positive-definite, and it is not available for a model with a custom distance.

```cs
using Numerics.MachineLearning;

// A deterministic two-input test function, sampled at 40 points
int n = 40;
var x = new double[n, 2];
var y = new double[n];
for (int i = 0; i < n; i++)
{
    x[i, 0] = (i * 0.6180339887) % 1.0 * 4.0;
    x[i, 1] = (i * 0.7548776662) % 1.0 * 4.0;
    y[i] = Math.Sin(x[i, 0]) + 0.5 * Math.Cos(2.0 * x[i, 1]);
}

var gp = new GaussianProcessRegression(x, y, GaussianProcessKernel.Matern52, GaussianProcessMean.Constant)
{
    NumberOfRestarts = 5,
    PRNGSeed = 12345
};
gp.Train();

Console.WriteLine($"Signal variance: {gp.SignalVariance:F4}");
Console.WriteLine($"Length scales: {string.Join(", ", gp.LengthScales.Select(l => l.ToString("F3")))}");
Console.WriteLine($"Noise variance: {gp.NoiseVariance:E2}");
Console.WriteLine($"Log marginal likelihood: {gp.LogMarginalLikelihood:F3}");

// Predictive mean and standard deviation at new points
double[,] xNew = { { 1.0, 1.0 }, { 2.5, 3.0 } };
var (mean, variance) = gp.PredictWithVariance(xNew);
for (int i = 0; i < mean.Length; i++)
    Console.WriteLine($"f({xNew[i, 0]}, {xNew[i, 1]}) = {mean[i]:F4} ± {1.96 * Math.Sqrt(variance[i]):F4}");

// Leave-one-out cross-validation without refitting
var (looMean, looVariance) = gp.LeaveOneOut();
double rmse = Math.Sqrt(Enumerable.Range(0, n).Average(i => Math.Pow(y[i] - looMean[i], 2)));
Console.WriteLine($"Leave-one-out RMSE: {rmse:F4}");

// The log marginal likelihood and its gradient at given log hyperparameters
double[] theta = { Math.Log(gp.SignalVariance), Math.Log(gp.LengthScales[0]), Math.Log(gp.LengthScales[1]), Math.Log(gp.NoiseVariance) };
var gradient = new double[theta.Length];
double lnp = gp.EvaluateLogMarginalLikelihood(theta, gradient);   // the fitted values give gp.LogMarginalLikelihood
```

**Leave-one-out cross-validation.** `LeaveOneOut()` returns, for every observation, the mean and variance predicted from the others with the same hyperparameters and re-estimated mean coefficients, in closed form [[18]](#18) at the cost of one extra matrix inversion. The standardized residuals $(y_i - \mu_{-i})/\sigma_{-i}$ should look like standard normal draws when the model is adequate.

**Numerical notes.** The covariance matrix is factored by a blocked Cholesky factorization with 256-bit vector instructions, and the work is split in fixed tiles so the results do not depend on the number of threads. When the matrix is not numerically positive-definite (for example duplicate inputs without noise), the factorization falls back to `CholeskyDecomposition.TryCholesky`, which adds the smallest diagonal jitter that works; it is reported in `Jitter` and acts as extra noise. Each likelihood evaluation costs $O(n^3)$ time and $O(n^2)$ memory, which is practical up to a few thousand observations. Results agree with scikit-learn `GaussianProcessRegressor` and R `DiceKriging::km` [[19]](#19) to about $10^{-15}$ for fixed hyperparameters; the oracle scripts are in [oracles](oracles/).

**Practical tips:**
- Scale the features to comparable ranges; the default bounds for each length scale are $10^{-3}$ to $10^{3}$ times the feature range.
- Use `Matern52` for physical responses; the squared exponential is often too smooth and makes the covariance matrix ill-conditioned.
- Add restarts with `NumberOfRestarts` when $d$ is large or the data are few; each start costs one optimization.
- For a deterministic simulator, `Train(0.0)` or a small fixed noise interpolates the runs.

### Kriging

Kriging [[20]](#20) is the geostatistical form of Gaussian process prediction, parameterized by a variogram instead of a covariance function. The `Kriging` class is a thin layer over `GaussianProcessRegression` with fixed hyperparameters. The semivariogram is

```math
\gamma(h) = c_0 + c\,(1 - \rho(h)) \quad (h > 0), \qquad \gamma(0) = 0
```

with nugget $c_0$, partial sill $c$ and range parameter $a$. The covariance is $C(h) = c\,\rho(h)$ for $h > 0$ and $C(0) = c_0 + c$.

| `VariogramModel` | $\rho(h)$ | Practical range |
|---|---|---|
| `Spherical` | $1 - 1.5(h/a) + 0.5(h/a)^3$ for $h < a$, else 0 | $a$ |
| `Exponential` | $\exp(-h/a)$ | $3a$ |
| `Gaussian` | $\exp(-(h/a)^2)$ | $\sqrt{3}\,a$ |

The range parameter follows the R package gstat [[21]](#21) (`vgm`). Isaaks and Srivastava [[20]](#20) write the exponential and Gaussian models with the practical range, so divide their range by 3 or $\sqrt{3}$.

**Kriging variants** (`KrigingType`):
- `Simple`: known mean $m$, passed as the `mean` argument.
- `Ordinary` (default): unknown constant mean.
- `Universal`: unknown mean linear in the coordinates, $\beta_0 + \sum_j \beta_j x_j$.

`KnownMean` holds the simple kriging mean (NaN for the other variants), and `TrendCoefficients` the estimated trend: empty for simple kriging, $[\beta_0]$ for ordinary kriging and $[\beta_0, \beta_1, \ldots, \beta_d]$ for universal kriging. `Semivariance(h)` and `Covariance(h)` evaluate the fitted variogram model at a separation distance $h \ge 0$; at $h = 0$ they return 0 and the total sill `Sill` $= c_0 + c$.

Ordinary and universal kriging estimate the trend by generalized least squares, which gives the same predictor and kriging variance as the kriging system with Lagrange multipliers. The nugget is part of the covariance at zero distance, so the prediction at an observed location equals the observation with zero kriging variance. Elsewhere the kriging variance $E[(Z(\mathbf{x}_0) - \hat Z(\mathbf{x}_0))^2]$ includes the nugget.

The following is the ordinary kriging example of Isaaks and Srivastava [[20]](#20) (chapter 12): seven samples and the covariance $C(h) = 10\exp(-0.3h)$, which is the exponential model with $c = 10$ and $a = 1/0.3$.

```cs
using Numerics.MachineLearning;

double[,] locations = { { 61, 139 }, { 63, 140 }, { 64, 129 }, { 68, 128 }, { 71, 140 }, { 73, 141 }, { 75, 128 } };
double[] values = { 477, 696, 227, 646, 606, 791, 783 };

var kriging = new Kriging(locations, values, VariogramModel.Exponential,
    nugget: 0.0, partialSill: 10.0, range: 1.0 / 0.3, type: KrigingType.Ordinary);

var (estimate, variance) = kriging.PredictWithVariance(new double[] { 65, 137 });
Console.WriteLine($"Estimate: {estimate:F1}");          // 592.7
Console.WriteLine($"Kriging variance: {variance:F2}");  // 8.96
Console.WriteLine($"Trend: {kriging.TrendCoefficients[0]:F1}");   // 605.1, the GLS estimate of the constant mean
Console.WriteLine($"γ(2) = {kriging.Semivariance(2.0):F4}, C(2) = {kriging.Covariance(2.0):F4}");   // 4.5119, 5.4881

// Leave-one-out cross-validation, equivalent to gstat::krige.cv
var (cvMean, cvVariance) = kriging.CrossValidate();
for (int i = 0; i < values.Length; i++)
    Console.WriteLine($"Sample {i}: observed {values[i]}, cross-validated {cvMean[i]:F1}");
```

**Custom distances.** A distance function can replace the Euclidean distance, for example a great-circle distance for latitude and longitude. The variogram model must give a valid covariance for that metric. The exponential model is valid with the great-circle distance on the sphere, because it is completely monotone [[22]](#22); a model that is valid in the plane need not be.

```cs
using Numerics.MachineLearning;

// Great-circle distance in kilometers between [latitude, longitude] pairs in degrees
Func<double[], double[], double> greatCircle = (a, b) =>
{
    double toRad = Math.PI / 180.0;
    double dLat = (b[0] - a[0]) * toRad, dLon = (b[1] - a[1]) * toRad;
    double h = Math.Pow(Math.Sin(dLat / 2), 2) + Math.Cos(a[0] * toRad) * Math.Cos(b[0] * toRad) * Math.Pow(Math.Sin(dLon / 2), 2);
    return 2.0 * 6371.0 * Math.Asin(Math.Min(1.0, Math.Sqrt(h)));
};

double[,] sites = { { 38.6, -90.2 }, { 39.1, -94.6 }, { 41.3, -96.0 }, { 36.2, -86.8 }, { 35.1, -90.0 } };
double[] errors = { 0.12, -0.05, 0.08, -0.10, 0.02 };

// Simple kriging of zero-mean residuals, as for spatially correlated regression errors
var residualKriging = new Kriging(sites, errors, VariogramModel.Exponential,
    nugget: 0.0, partialSill: 0.01, range: 400.0, distance: greatCircle, type: KrigingType.Simple, mean: 0.0);
var (error, errorVariance) = residualKriging.PredictWithVariance(new double[] { 38.0, -92.0 });
Console.WriteLine($"Predicted error: {error:F4} (variance {errorVariance:F5})");
```

**Replacing a hand-written simple kriging.** Code that builds $K$ from a correlation function, factors it, and returns $\mathbf{k}_*^\top K^{-1}\boldsymbol\varepsilon$ and $\sigma^2 - \mathbf{k}_*^\top K^{-1}\mathbf{k}_*$ for a zero-mean process is `Kriging` with `KrigingType.Simple`, `mean: 0`, `nugget: 0`, `partialSill` equal to $\sigma^2$ and the same distance function. The class reuses one factorization for all predictions, adds jitter instead of failing when the matrix is singular, and also offers the ordinary and universal variants and cross-validation.

**Accuracy.** Predictions, variances and cross-validation agree with gstat `krige` and `krige.cv` on the `sp::meuse` data to about $10^{-13}$ relative ($10^{-11}$ for universal kriging, whose trend uses large coordinates), and the textbook example reproduces the published estimate, weights and variance. The oracle scripts are in [oracles](oracles/).

**Performance.** Fitting and prediction are benchmarked against scikit-learn and DiceKriging with the scripts in [benchmarks](benchmarks/) on the same deterministic data: $d = 10$ features, a Matérn 5/2 kernel with automatic relevance determination and estimated noise (Numerics with its two default starts and a constant mean; scikit-learn and DiceKriging with one start), and the mean and variance at 1000 new points. Median wall times on an Intel Core Ultra 9 185H (scikit-learn 1.8.0 with OpenBLAS 0.3.31.dev; DiceKriging 1.6.1 with the R 4.4.3 reference BLAS):

| $n$ | Fit, Numerics | Fit, scikit-learn | Fit, DiceKriging | Predict, Numerics | Predict, scikit-learn | Predict, DiceKriging | ln p, Numerics | ln p, scikit-learn |
|---|---|---|---|---|---|---|---|---|
| 100 | 0.064 s | 0.26 s | 0.60 s | 0.002 s | 0.004 s | 0.06 s | -34.18 | -34.18 |
| 500 | 0.67 s | 5.3 s | 16.7 s | 0.006 s | 0.021 s | 0.34 s | -25.24 | -45.15 |
| 1000 | 4.7 s | 31 s | 111 s | 0.023 s | 0.041 s | 0.73 s | 211.15 | 153.74 |
| 2000 | 23 s | 136 s | 2086 s | 0.060 s | 0.099 s | 41 s | 904.82 | 904.81 |

Fit times include the optimizer iterations, whose number differs between libraries. With the same start, bounds and parameterization as scikit-learn, `LBFGSB` reaches scikit-learn's log marginal likelihood at every size; the higher values at $n$ = 500 and 1000 come from the second default start, which finds a better local maximum. Timings on a laptop processor with performance and efficiency cores vary by a factor of two or more between runs.

### Elastic Net Regression

Penalized linear regression for variable selection and shrinkage [[12]](#12) [[13]](#13). `ElasticNetRegression` fits a Gaussian linear model along a decreasing sequence of penalties $\lambda$, minimizing at each one

```math
\min_{\beta_0,\,\beta}\; \frac{1}{2n}\sum_{i=1}^{n}\left(y_i - \beta_0 - \mathbf{x}_i^\top\beta\right)^2 + \lambda \sum_{j=1}^{p} v_j\left[\frac{1-\alpha}{2}\beta_j^2 + \alpha\,|\beta_j|\right]
```

where $\alpha \in [0, 1]$ mixes ridge regression ($\alpha = 0$) and the lasso ($\alpha = 1$), and $v_j$ are optional penalty factors (`PenaltyFactors`; $v_j = 0$ leaves a predictor unpenalized). The $L_1$ part sets coefficients exactly to zero, so the set of nonzero coefficients along the path is a ranking of the predictors; the $L_2$ part keeps groups of correlated predictors together instead of picking one at random. Typical uses are covariate selection for regional regression and screening of calibration parameters.

**Algorithm.** The path is computed by cyclic coordinate descent as in glmnet [[12]](#12). The predictors are standardized to unit population variance (`Standardize`, default true) and the response to unit population variance; coefficients are reported on the original scale. On the standardized scale each coordinate update is a soft-thresholding step on the current gradient $g_j$, with $\lambda$ divided by $s_y$:

```math
\beta_j \leftarrow \frac{S\!\left(g_j + \beta_j,\; \lambda\alpha v_j\right)}{1 + \lambda(1-\alpha)v_j}, \qquad S(u, t) = \operatorname{sign}(u)\max(|u| - t, 0)
```

The gradient is kept current from cached inner products of the predictors (covariance updates), so an update costs $O(p)$ rather than $O(n)$. Full sweeps over all predictors alternate with sweeps over the active set, and each solution warm-starts the next penalty. The default sequence has 100 penalties log-spaced from $\lambda_{\max}$, the smallest penalty with all coefficients zero, down to $\epsilon\lambda_{\max}$ with $\epsilon = 10^{-4}$ ($10^{-2}$ when $n < p$). A sweep has converged when the largest weighted squared coefficient change is below `Tolerance` (default $10^{-7}$), and the path stops early when the fraction of deviance explained exceeds 0.999 or changes by less than $10^{-5}$.

**Cross-validation.** `CrossValidate(folds, seed)` assigns observations to folds at random (reproducibly from the seed), or `CrossValidate(foldIds)` uses given fold labels. The labels of the last run are kept in `FoldIds`, so `enet.CrossValidate(enet.FoldIds)` repeats a seeded run with the same folds and the same errors. Each fold is solved exactly at the full-data penalties, and the error is the fold-size weighted mean squared prediction error with standard error

```math
\mathrm{SE}(\lambda) = \sqrt{\frac{\sum_k n_k\,(\mathrm{MSE}_k(\lambda) - \mathrm{CV}(\lambda))^2}{n\,(K-1)}}
```

`LambdaMin` is the largest penalty with the smallest error, and `LambdaOneStandardError` is the largest penalty whose error is within one standard error of that minimum, a sparser model whose error is statistically indistinguishable from the best [[9]](#9).

The example uses the prostate cancer data of Stamey et al. [[15]](#15), the standard lasso example in Hastie, Tibshirani and Friedman [[9]](#9) (`docs/example-data/prostate-cancer.csv`):

```cs
using System;
using System.IO;
using System.Linq;
using Numerics.MachineLearning;

// Prostate cancer training set (docs/example-data/prostate-cancer.csv): 67 cases, 8 predictors
double[][] rows = File.ReadAllLines("prostate-cancer.csv")
    .Where(line => !line.StartsWith("#"))
    .Skip(1)                                     // header
    .Select(line => line.Split(',').Select(double.Parse).ToArray())
    .Where(r => r[9] == 1)                       // training cases
    .ToArray();
string[] names = { "lcavol", "lweight", "age", "lbph", "svi", "lcp", "gleason", "pgg45" };
var X = new double[rows.Length, 8];
var y = new double[rows.Length];
for (int i = 0; i < rows.Length; i++)
{
    for (int j = 0; j < 8; j++) X[i, j] = rows[i][j];
    y[i] = rows[i][8];                           // lpsa
}

// Fit the full penalty path (alpha = 0.5 mixes ridge and lasso equally)
var enet = new ElasticNetRegression(X, y) { Alpha = 0.5 };
enet.Fit();
int last = enet.Lambdas.Length - 1;
Console.WriteLine($"Path: {enet.Lambdas.Length} penalties from {enet.Lambdas[0]:F4} to {enet.Lambdas[last]:F6}");

// 10-fold cross-validation with a seeded fold assignment
enet.CrossValidate(folds: 10, seed: 12345);
Console.WriteLine($"lambda_min = {enet.LambdaMin:F4}, lambda_1se = {enet.LambdaOneStandardError:F4}");

// Coefficients and selected variables at lambda_1se
int k = enet.LambdaOneStandardErrorIndex;
double[] beta = enet.GetCoefficients(k);
Console.WriteLine($"Selected: {string.Join(", ", enet.SelectedVariables(k).Select(j => names[j]))}");
Console.WriteLine($"  intercept {enet.Intercepts[k],9:F4}");
for (int j = 0; j < 8; j++) Console.WriteLine($"  {names[j],-9} {beta[j],9:F4}");

// Predict a new case at lambda_1se
double[] patient = { 1.5, 3.6, 65, 0.1, 0, -0.5, 7, 15 };
Console.WriteLine($"Predicted lpsa: {enet.Predict(patient, k):F3}");

// Keep lcavol in every model (penalty factor 0) and trace the lasso path
var lasso = new ElasticNetRegression(X, y)
{
    Alpha = 1.0,
    PenaltyFactors = new double[] { 0, 1, 1, 1, 1, 1, 1, 1 }
};
lasso.Fit();
for (int m = 0; m < lasso.Lambdas.Length; m += 15)
    Console.WriteLine($"lambda = {lasso.Lambdas[m]:F4}  nonzero = {lasso.NonzeroCount[m]}  R2 = {lasso.DevianceRatio[m]:F3}");
```

Output:

```
Path: 73 penalties from 1.7578 to 0.002167
lambda_min = 0.0322, lambda_1se = 0.3967
Selected: lcavol, lweight, lbph, svi, pgg45
  intercept    0.4846
  lcavol       0.3681
  lweight      0.3802
  age          0.0000
  lbph         0.0171
  svi          0.3148
  lcp          0.0000
  gleason      0.0000
  pgg45        0.0013
Predicted lpsa: 2.427
lambda = 0.2781  nonzero = 1  R2 = 0.538
lambda = 0.0689  nonzero = 5  R2 = 0.648
lambda = 0.0171  nonzero = 7  R2 = 0.690
lambda = 0.0042  nonzero = 8  R2 = 0.694
```

**Key members:**

| Member | Description |
|---|---|
| `Alpha` | Mixing parameter $\alpha$ (default 1, the lasso) |
| `NumberOfLambdas`, `LambdaMinRatio`, `UserLambdas` | Default path length (100) and $\epsilon$, or a user-supplied penalty sequence |
| `Standardize`, `HasIntercept` | Standardize the predictors (default true); fit an unpenalized intercept (default true) |
| `PenaltyFactors` | Relative penalty $v_j$ per predictor, rescaled to sum to $p$ |
| `Tolerance`, `MaxIterations` | Convergence threshold ($10^{-7}$) and maximum total sweeps (100,000) |
| `Lambdas`, `Coefficients`, `Intercepts` | The fitted path: penalties, coefficients `[penalty, predictor]`, intercepts |
| `DevianceRatio`, `NonzeroCount` | Fraction of variance explained and number of nonzero coefficients per penalty |
| `CrossValidationError`, `CrossValidationStandardError` | Cross-validated error and its standard error per penalty |
| `LambdaMinIndex`, `LambdaOneStandardErrorIndex` | Indexes of $\lambda_{\min}$ and $\lambda_{1se}$ |
| `GetCoefficients(k)`, `SelectedVariables(k)`, `Predict(x, k)` | Coefficients, nonzero predictors and predictions at path index `k` |

**Relation to glmnet and scikit-learn.** With the same settings, `ElasticNetRegression` reproduces R glmnet's penalties, coefficients, deviance ratios and path length to about $10^{-13}$, for $\alpha \in [0, 1]$, with or without intercept, standardization and penalty factors. scikit-learn parameterizes the same model differently and does not standardize the response:

| ElasticNetRegression / glmnet | scikit-learn `ElasticNet` |
|---|---|
| `Alpha` = $\alpha$ | `l1_ratio` = $\alpha$ |
| penalty $\lambda$ | `alpha` = $\lambda / s_y$ on standardized $X$ and $y$, where $s_y$ is the population standard deviation of $y$; for the lasso ($\alpha = 1$) the scaling cancels and `alpha` = $\lambda$ |
| `Standardize = true` | divide each column of $X$ by its population standard deviation before fitting, then back-transform $\beta_j = b_j\,s_y / s_{x_j}$ |

Cross-validation solves every fold at the full-data penalties, as scikit-learn's `ElasticNetCV` does. R's `cv.glmnet` fits each fold on its own penalty sequence and interpolates the predictions, so its cross-validated errors can differ slightly; on the prostate data both select the same $\lambda_{\min}$ and $\lambda_{1se}$.

**Performance.** The inner products are computed only for predictors that enter the model, in batches, by a cache-blocked SIMD kernel, and the cross-validation folds share them through per-fold partial sums, so $k$-fold cross-validation costs about one more pass over the data rather than $k$ refits. Fits and cross-validation are reproducible on a given machine; because the SIMD register width sets the order of the partial sums, the last digits can differ between processor architectures.

Median wall time on identical data (Intel Core Ultra 9 185H (16 cores), Windows 11, .NET 10; Python 3.14.3, NumPy 2.4.2, SciPy 1.17.1, scikit-learn 1.8.0 with OpenBLAS 0.3.31 (scipy-openblas, multithreaded); R 4.4.3 with glmnet 4.1.10 and R's reference BLAS and LAPACK; scripts in `docs/machine-learning/benchmarks`), $\alpha = 0.5$, 100 penalties, 10 folds:

| n × p | Path: C# | scikit-learn | glmnet | 10-fold CV: C# | scikit-learn | cv.glmnet |
|---|---|---|---|---|---|---|
| 100 × 200 | 2.5 ms | 10.9 ms | 5.2 ms | 27 ms | 125 ms | 110 ms |
| 1,000 × 50 | 1.3 ms | 10.5 ms | 4.1 ms | 8.3 ms | 102 ms | 165 ms |
| 10,000 × 100 | 15 ms | 33 ms | 45 ms | 46 ms | 271 ms | 890 ms |
| 100,000 × 500 | 0.82 s | 16.6 s | 6.7 s | 1.6 s | 203 s | 48 s |

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

**$k$-Means++ initialization** [[11]](#11). The default initialization method selects initial centroids that are well-separated. The first centroid is chosen uniformly at random. Each subsequent centroid is selected with probability proportional to $D(\mathbf{x})$, where $D(\mathbf{x})$ is the Euclidean distance from $\mathbf{x}$ to the nearest already-chosen centroid. This is a variant of the method of Arthur and Vassilvitskii, which samples proportional to $D(\mathbf{x})^2$; on clustered data the two seedings give fits of similar quality. This initialization significantly reduces the chance of poor convergence compared to random initialization.

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

**Expectation-Maximization (EM) algorithm.** The parameters are estimated by maximizing the log-likelihood using the EM algorithm, initialized from a $k$-means solution with equal weights. Each initial covariance is the within-cluster covariance of the $k$-means cluster, with each diagonal entry floored at $10^{-6}$ times the variance of that column over all the data.

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

### Principal Component Analysis (PCA)

Orthogonal directions of maximum variance [[14]](#14). With the $n \times p$ data matrix centered (and, on the correlation basis, scaled to unit sample standard deviation) to $Z$, the sample covariance matrix has the eigen-decomposition

```math
S = \frac{1}{n-1} Z^\top Z = V \Lambda V^\top, \qquad \lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_p \ge 0
```

The columns of $V$ are the principal axes, $\lambda_i$ is the variance explained by component $i$, and the scores are $T = ZV$. Use the covariance basis (`PrincipalComponentBasis.Covariance`, the default) when the variables share units, and the correlation basis when they do not; on the correlation basis the eigenvalues sum to $p$ and the loadings $L_{ji} = V_{ji}\sqrt{\lambda_i}$ are the correlations between variable $j$ and component $i$.

**Algorithm.** $S$ is accumulated directly from the rows of the data by a SIMD cross-product update (no transposed copy) and decomposed by Householder tridiagonalization and the implicit QL algorithm [[16]](#16), at a cost of $O(np^2 + p^3)$. Each axis is oriented so that its largest-magnitude element is positive, the convention of scikit-learn's `PCA`. R's `prcomp` does not fix the sign, so its rotation and scores agree up to a sign per component.

```cs
using System;
using System.IO;
using System.Linq;
using Numerics.MachineLearning;

// Fisher's iris measurements (docs/example-data/iris-dataset.csv)
double[][] rows = File.ReadAllLines("iris-dataset.csv")
    .Where(line => !line.StartsWith("#"))
    .Skip(1)                                     // header
    .Select(line => line.Split(',').Select(double.Parse).ToArray())
    .ToArray();
string[] names = { "SepalLength", "SepalWidth", "PetalLength", "PetalWidth" };
var X = new double[rows.Length, 4];
for (int i = 0; i < rows.Length; i++)
    for (int j = 0; j < 4; j++) X[i, j] = rows[i][j];

// Principal components of the correlation matrix (variables in different units)
var pca = new PrincipalComponentAnalysis(X, PrincipalComponentBasis.Correlation);
pca.Train();
for (int i = 0; i < pca.NumberOfComponents; i++)
{
    Console.WriteLine($"PC{i + 1}: variance {pca.ExplainedVariance[i]:F4}  " +
                      $"ratio {pca.ExplainedVarianceRatio[i]:P1}  cumulative {pca.CumulativeExplainedVarianceRatio[i]:P1}");
}

// Loadings: correlations between each variable and the first two components
for (int j = 0; j < 4; j++)
    Console.WriteLine($"  {names[j],-12} {pca.Loadings[j, 0],7:F3} {pca.Loadings[j, 1],7:F3}");

// Scores of the first observation
Console.WriteLine($"Scores[0]: {pca.Scores[0, 0]:F4}, {pca.Scores[0, 1]:F4}");

// Keep two components, project a new flower and map it back
var two = new PrincipalComponentAnalysis(X, PrincipalComponentBasis.Correlation, numberOfComponents: 2);
two.Train();
double[] t = two.Transform(new double[] { 6.0, 3.0, 4.5, 1.5 });
double[] back = two.InverseTransform(t);
Console.WriteLine($"Scores: {t[0]:F4}, {t[1]:F4}");
Console.WriteLine($"Reconstruction: {string.Join(", ", back.Select(v => v.ToString("F3")))}");
```

Output (the variances match `prcomp(iris[, 1:4], scale. = TRUE)$sdev^2`):

```
PC1: variance 2.9142  ratio 72.9%  cumulative 72.9%
PC2: variance 0.9172  ratio 22.9%  cumulative 95.8%
PC3: variance 0.1479  ratio 3.7%  cumulative 99.5%
PC4: variance 0.0207  ratio 0.5%  cumulative 100.0%
  SepalLength    0.890   0.360
  SepalWidth    -0.456   0.885
  PetalLength    0.992   0.022
  PetalWidth     0.965   0.063
Scores[0]: -2.2618, 0.4843
Scores: 0.5986, -0.0128
Reconstruction: 6.099, 2.981, 4.372, 1.457
```

**Key members:**

| Member | Description |
|---|---|
| `Basis`, `NumberOfComponents` | Covariance or correlation basis; number of retained components $k$ (default $\min(n, p)$) |
| `Means`, `Scales` | Column means, and column standard deviations on the correlation basis (ones otherwise) |
| `Components` | $k \times p$ principal axes (scikit-learn `components_`; the transpose of `prcomp` `rotation`) |
| `ExplainedVariance`, `ExplainedVarianceRatio`, `CumulativeExplainedVarianceRatio` | $\lambda_i$ (`prcomp` `sdev^2`), $\lambda_i / \sum \lambda$ and its running sum |
| `SingularValues`, `TotalVariance` | $\sqrt{(n-1)\lambda_i}$ and the trace of $S$ |
| `Loadings`, `Scores` | $p \times k$ loadings $V\sqrt{\Lambda}$; $n \times k$ scores of the training data (computed on first access) |
| `Transform`, `InverseTransform` | Project observations onto the components; map scores back to the variables (exact when $k = p$) |

**Numerical notes.** Forming $S$ squares the condition number of $Z$, so eigenvalues much smaller than $\epsilon\lambda_1$ lose relative accuracy and round-off below zero is set to zero; with more variables than observations the last $p - n + 1$ variances are zero. A constant column has zero variance on the covariance basis and is rejected on the correlation basis. Results are reproducible on a given machine; the SIMD register width sets the order of the partial sums, so the last digits can differ between processor architectures.

**Performance.** Median wall time to fit all components on identical data (Intel Core Ultra 9 185H (16 cores), Windows 11, .NET 10; Python 3.14.3, NumPy 2.4.2, SciPy 1.17.1, scikit-learn 1.8.0 with OpenBLAS 0.3.31 (scipy-openblas, multithreaded); R 4.4.3 with glmnet 4.1.10 and R's reference BLAS and LAPACK; scripts in `docs/machine-learning/benchmarks`):

| n × p | C# | scikit-learn `PCA().fit` | R `prcomp` |
|---|---|---|---|
| 1,000 × 20 | 0.24 ms | 0.35 ms | 1.5 ms |
| 10,000 × 50 | 1.1 ms | 2.7 ms | 78 ms |
| 100,000 × 100 | 16 ms | 48 ms | 4.4 s |

---

## Practical Examples

### Example 1: Regression with GLM

```cs
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Functions;

// Predict home prices (no intercept column; GLM adds it automatically)
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
Console.WriteLine($"  90% CI for the mean: [${interval[0, 0]:F0}, ${interval[0, 2]:F0}]");
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

// Train random forest (no nTrees constructor parameter; set it via the property)
var rf = new RandomForest(X_train, y_train, seed: 42);
rf.IsRegression = false;
rf.NumberOfTrees = 100;
rf.Train();

// Evaluate: Predict returns double[,] with columns: lower(0), median(1), upper(2), mean(3)
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

Predicting whether flood damage occurs based on hydraulic variables using data from [`example-data/flood-damage-glm.csv`](../example-data/flood-damage-glm.csv):

These synthetic data are completely separated: a combination of stage and velocity divides the damaged from the undamaged structures exactly, so the maximum likelihood estimates do not exist (R's `glm` reports fitted probabilities of 0 or 1 and coefficients in the thousands). Fisher scoring does not converge, the optimizer stops at its fixed bounds, and the coefficients and standard errors below are not meaningful. The example shows the API; with real data, check for separation before interpreting a binary GLM.

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

// 90% confidence interval for the mean damage probability
double[,] interval = glm.Predict(new Matrix(newEvent), alpha: 0.1);
Console.WriteLine($"90% CI for the mean: [{interval[0, 0]:P1}, {interval[0, 2]:P1}]");
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
    var model = new DecisionTree(X_trainFold, y_trainFold) { IsRegression = false };
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

// Train classifiers (KNN is lazy; no Train() call needed)
var decisionTree = new DecisionTree(X, y) { IsRegression = false };
decisionTree.Train();

var randomForest = new RandomForest(X, y, seed: 42) { IsRegression = false };
randomForest.NumberOfTrees = 50;
randomForest.Train();

var knn = new KNearestNeighbors(X, y, k: 3) { IsRegression = false };

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
7. **Select covariates with a penalty path** - Use `ElasticNetRegression` with cross-validation and prefer $\lambda_{1se}$ for a sparser model

### Unsupervised Learning
1. **Scale features** - Clustering sensitive to feature scales
2. **Choose k carefully** - Use elbow method or silhouette scores
3. **Multiple runs** - k-Means sensitive to initialization
4. **Validate clusters** - Inspect cluster characteristics
5. **Consider GMM** - When clusters overlap or have different shapes
6. **Reduce dimension first** - Use PCA on the correlation basis when variables have different units

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

<a id="12">[12]</a> Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software*, 33(1), 1-22. https://doi.org/10.18637/jss.v033.i01

<a id="13">[13]</a> Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. *Journal of the Royal Statistical Society: Series B*, 67(2), 301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x

<a id="14">[14]</a> Jolliffe, I. T. (2002). *Principal Component Analysis* (2nd ed.). Springer. https://doi.org/10.1007/b98835

<a id="15">[15]</a> Stamey, T. A., Kabalin, J. N., McNeal, J. E., Johnstone, I. M., Freiha, F., Redwine, E. A., & Yang, N. (1989). Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy treated patients. *Journal of Urology*, 141(5), 1076-1083. https://doi.org/10.1016/S0022-5347(17)41175-X

<a id="16">[16]</a> Wilkinson, J. H., & Reinsch, C. (1971). *Handbook for Automatic Computation, Volume II: Linear Algebra*. Springer. https://doi.org/10.1007/978-3-642-86940-2

<a id="17">[17]</a> Rasmussen, C. E., & Williams, C. K. I. (2006). *Gaussian Processes for Machine Learning*. MIT Press.

<a id="18">[18]</a> Dubrule, O. (1983). Cross validation of kriging in a unique neighborhood. *Mathematical Geology*, 15(6), 687-699. https://doi.org/10.1007/BF01033232

<a id="19">[19]</a> Roustant, O., Ginsbourger, D., & Deville, Y. (2012). DiceKriging, DiceOptim: Two R packages for the analysis of computer experiments by kriging-based metamodeling and optimization. *Journal of Statistical Software*, 51(1), 1-55. https://doi.org/10.18637/jss.v051.i01

<a id="20">[20]</a> Isaaks, E. H., & Srivastava, R. M. (1989). *An Introduction to Applied Geostatistics*. Oxford University Press.

<a id="21">[21]</a> Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat package. *Computers & Geosciences*, 30(7), 683-691. https://doi.org/10.1016/j.cageo.2004.03.012

<a id="22">[22]</a> Gneiting, T. (2013). Strictly and non-strictly positive definite functions on spheres. *Bernoulli*, 19(4), 1327-1349. https://doi.org/10.3150/12-BEJSP06

---

[← Previous: Univariate Functions](../functions/index.md) | [Back to Index](../index.md) | [Next: Random Generation →](../sampling/random-generation.md)
