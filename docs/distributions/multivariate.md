# Multivariate Distributions

[← Previous: Copulas](copulas.md) | [Back to Index](../index.md) | [Next: Machine Learning →](../machine-learning/machine-learning.md)

The ***Numerics*** library provides several multivariate distributions for modeling correlated random variables: the **Multivariate Normal**, **Multivariate Student-t**, **Dirichlet**, and **Multinomial** distributions. These are fundamental in multivariate statistics, risk assessment, and uncertainty quantification.

## Multivariate Normal Distribution

The multivariate normal (Gaussian) distribution generalizes the univariate normal to multiple dimensions with a specified covariance structure [[1]](#1).

### Mathematical Definition

A random vector $\mathbf{X} = (X_1, X_2, \ldots, X_k)^T$ follows a $k$-dimensional multivariate normal distribution, written $\mathbf{X} \sim \mathcal{N}_k(\boldsymbol{\mu}, \boldsymbol{\Sigma})$, if its probability density function is:

```math
f(\mathbf{x}) = (2\pi)^{-k/2} \, |\boldsymbol{\Sigma}|^{-1/2} \, \exp\!\left( -\tfrac{1}{2} (\mathbf{x} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{x} - \boldsymbol{\mu}) \right)
```

where:
- $\boldsymbol{\mu} \in \mathbb{R}^k$ is the mean vector,
- $\boldsymbol{\Sigma}$ is the $k \times k$ positive-definite covariance matrix,
- $|\boldsymbol{\Sigma}|$ is the determinant of $\boldsymbol{\Sigma}$.

**Mahalanobis distance.** The quadratic form in the exponent defines the squared Mahalanobis distance between a point $\mathbf{x}$ and the distribution center $\boldsymbol{\mu}$:

```math
D^2(\mathbf{x}) = (\mathbf{x} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{x} - \boldsymbol{\mu})
```

The PDF can be expressed compactly in terms of this distance as $`f(\mathbf{x}) = (2\pi)^{-k/2} |\boldsymbol{\Sigma}|^{-1/2} \exp(-D^2/2)`$. Surfaces of constant density are ellipsoids defined by $D^2 = c$, and the squared Mahalanobis distance follows a chi-squared distribution: $D^2 \sim \chi^2_k$. This property is used for multivariate outlier detection -- a point is flagged as an outlier if $D^2$ exceeds the $\chi^2_k$ critical value at the desired significance level.

**Marginal distributions.** Any subset of variables from a multivariate normal is itself multivariate normal. If the full vector is partitioned as $\mathbf{X} = (\mathbf{X}_a, \mathbf{X}_b)^T$ with corresponding partitioned mean and covariance:

```math
\boldsymbol{\mu} = \begin{pmatrix} \boldsymbol{\mu}_a \\ \boldsymbol{\mu}_b \end{pmatrix}, \quad \boldsymbol{\Sigma} = \begin{pmatrix} \boldsymbol{\Sigma}_{aa} & \boldsymbol{\Sigma}_{ab} \\ \boldsymbol{\Sigma}_{ba} & \boldsymbol{\Sigma}_{bb} \end{pmatrix}
```

then the marginal distribution is $`\mathbf{X}_a \sim \mathcal{N}(\boldsymbol{\mu}_a, \boldsymbol{\Sigma}_{aa})`$.

**Conditional distributions.** The conditional distribution of $\mathbf{X}_a$ given $\mathbf{X}_b = \mathbf{x}_b$ is also multivariate normal:

```math
\mathbf{X}_a \mid \mathbf{X}_b = \mathbf{x}_b \sim \mathcal{N}\!\left( \boldsymbol{\mu}_{a|b}, \, \boldsymbol{\Sigma}_{a|b} \right)
```

with conditional mean and covariance:

```math
\boldsymbol{\mu}_{a|b} = \boldsymbol{\mu}_a + \boldsymbol{\Sigma}_{ab} \boldsymbol{\Sigma}_{bb}^{-1} (\mathbf{x}_b - \boldsymbol{\mu}_b)
```

```math
\boldsymbol{\Sigma}_{a|b} = \boldsymbol{\Sigma}_{aa} - \boldsymbol{\Sigma}_{ab} \boldsymbol{\Sigma}_{bb}^{-1} \boldsymbol{\Sigma}_{ba}
```

Note that the conditional covariance $\boldsymbol{\Sigma}_{a|b}$ does not depend on the observed value $\mathbf{x}_b$. For the bivariate case ($k = 2$), these reduce to:

```math
\mu_{Y|X=x} = \mu_Y + \rho \frac{\sigma_Y}{\sigma_X}(x - \mu_X), \quad \sigma^2_{Y|X} = \sigma^2_Y (1 - \rho^2)
```

**Cholesky sampling.** Random samples are generated via the Cholesky decomposition $\boldsymbol{\Sigma} = \mathbf{L}\mathbf{L}^T$, where $\mathbf{L}$ is a lower-triangular matrix. Given a vector of independent standard normal variates $\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_k)$:

```math
\mathbf{x} = \boldsymbol{\mu} + \mathbf{L}\mathbf{z}
```

This produces a sample $\mathbf{x} \sim \mathcal{N}_k(\boldsymbol{\mu}, \boldsymbol{\Sigma})$. The Cholesky approach is preferred because it is computationally efficient ($O(k^3)$ once, then $O(k^2)$ per sample), and the decomposition itself serves as a check that $\boldsymbol{\Sigma}$ is positive definite.

### Creating Multivariate Normal Distributions

```cs
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;

// Example 1: Standard multivariate normal (dimension 3)
// Mean = [0, 0, 0], Covariance = Identity
var mvn1 = new MultivariateNormal(dimension: 3);

// Example 2: Specified mean, identity covariance
double[] mean = { 100, 50, 75 };
var mvn2 = new MultivariateNormal(mean);

// Example 3: Full specification with covariance matrix
double[] mu = { 100, 50, 75 };
double[,] sigma = {
    { 225,  75,  50 },  // Var(X₁)=225, Cov(X₁,X₂)=75, Cov(X₁,X₃)=50
    {  75, 100,  25 },  // Var(X₂)=100, Cov(X₂,X₃)=25
    {  50,  25, 144 }   // Var(X₃)=144
};

var mvn3 = new MultivariateNormal(mu, sigma);

Console.WriteLine($"Multivariate Normal Distribution:");
Console.WriteLine($"  Dimension: {mvn3.Dimension}");
Console.WriteLine($"  Mean vector: [{string.Join(", ", mvn3.Mean)}]");
```

### Properties

```cs
var mvn = new MultivariateNormal(mu, sigma);

// Basic properties
int dim = mvn.Dimension;               // Number of variables
double[] mean = mvn.Mean;              // Mean vector μ
double[,] cov = mvn.Covariance;        // Covariance matrix Σ

// Individual variances (diagonal of covariance)
double var1 = cov[0, 0];               // Var(X₁) = 225
double var2 = cov[1, 1];               // Var(X₂) = 100
double var3 = cov[2, 2];               // Var(X₃) = 144

// Covariances (off-diagonal)
double cov12 = cov[0, 1];              // Cov(X₁, X₂) = 75
double cov13 = cov[0, 2];              // Cov(X₁, X₃) = 50
double cov23 = cov[1, 2];              // Cov(X₂, X₃) = 25

Console.WriteLine($"Standard deviations: [{Math.Sqrt(var1):F1}, {Math.Sqrt(var2):F1}, {Math.Sqrt(var3):F1}]");
```

### Probability Density Function (PDF)

```cs
double[] x = { 105, 55, 80 };

// Compute PDF at point x
double pdf = mvn.PDF(x);
double logPdf = mvn.LogPDF(x);

Console.WriteLine($"PDF at x = [{string.Join(", ", x)}]:");
Console.WriteLine($"  f(x) = {pdf:E4}");
Console.WriteLine($"  log f(x) = {logPdf:F4}");
```

The PDF formula is given in the [Mathematical Definition](#mathematical-definition) section above. The library computes it as $f(\mathbf{x}) = \exp(-\tfrac{1}{2}D^2 + \ln C)$, where $D^2$ is the Mahalanobis distance and $\ln C = -\tfrac{1}{2}(k \ln 2\pi + \ln|\boldsymbol{\Sigma}|)$ is a precomputed constant.

### Cumulative Distribution Function (CDF)

```cs
double[] x = { 105, 55, 80 };

// Compute CDF at point x
// P(X₁ ≤ 105, X₂ ≤ 55, X₃ ≤ 80)
double cdf = mvn.CDF(x);

Console.WriteLine($"CDF at x = [{string.Join(", ", x)}]:");
Console.WriteLine($"  P(X ≤ x) = {cdf:F6}");

// For dimensions > 2, uses Monte Carlo integration
// Can control accuracy with properties:
mvn.MaxEvaluations = 100000;    // Max function evaluations
mvn.AbsoluteError = 1e-4;       // Absolute error tolerance
mvn.RelativeError = 1e-4;       // Relative error tolerance
```

### Bivariate CDF (Special Case)

For two-dimensional case, exact calculation available:

```cs
// Bivariate standard normal with correlation ρ
double z1 = 1.0;
double z2 = 1.5;
double rho = 0.6;

double bivCDF = MultivariateNormal.BivariateCDF(z1, z2, rho);

Console.WriteLine($"Bivariate CDF:");
Console.WriteLine($"  P(Z₁ ≤ {z1}, Z₂ ≤ {z2} | ρ={rho}) = {bivCDF:F6}");

// Useful for correlation analysis and joint probabilities
```

### Inverse CDF (Quantile Function)

The `InverseCDF(double[] probabilities)` method takes one probability per dimension and returns a single point in k-dimensional space where each marginal is at its respective probability level:

```cs
// Get the point where dim1 is at 5%, dim2 is at 50%, dim3 is at 95%
double[] point = mvn.InverseCDF(new double[] { 0.05, 0.50, 0.95 });

Console.WriteLine("Multivariate quantile point:");
for (int i = 0; i < mvn.Dimension; i++)
{
    Console.WriteLine($"  X{i + 1} at p={new[] { 0.05, 0.50, 0.95 }[i]}: {point[i]:F2}");
}

// For marginal quantiles, use the marginal normal distributions directly
Console.WriteLine("\nMarginal quantile ranges:");
for (int i = 0; i < mvn.Dimension; i++)
{
    double sd = Math.Sqrt(mvn.Covariance[i, i]);
    double q05 = mvn.Mean[i] + sd * new Normal(0, 1).InverseCDF(0.05);
    double q50 = mvn.Mean[i];
    double q95 = mvn.Mean[i] + sd * new Normal(0, 1).InverseCDF(0.95);

    Console.WriteLine($"  X{i + 1}: 5%={q05:F2}, 50%={q50:F2}, 95%={q95:F2}");
}
```

### Random Sample Generation

```cs
// Generate random samples
int n = 1000;
int seed = 12345;

double[,] samples = mvn.GenerateRandomValues(n, seed);

Console.WriteLine($"Generated {n} samples from MVN distribution");

// Compute sample statistics
double[] sampleMeans = new double[mvn.Dimension];
double[,] sampleCov = new double[mvn.Dimension, mvn.Dimension];

// Calculate sample means
for (int j = 0; j < mvn.Dimension; j++)
{
    for (int i = 0; i < n; i++)
    {
        sampleMeans[j] += samples[i, j];
    }
    sampleMeans[j] /= n;
}

Console.WriteLine($"Sample means: [{string.Join(", ", sampleMeans.Select(m => m.ToString("F2")))}]");
Console.WriteLine($"True means:   [{string.Join(", ", mvn.Mean.Select(m => m.ToString("F2")))}]");
```

## Correlation Structure

### Creating Correlation Matrix

```cs
// Convert covariance to correlation
int dim = 3;
double[,] cov = mvn.Covariance;
double[,] corr = new double[dim, dim];

for (int i = 0; i < dim; i++)
{
    for (int j = 0; j < dim; j++)
    {
        corr[i, j] = cov[i, j] / Math.Sqrt(cov[i, i] * cov[j, j]);
    }
}

Console.WriteLine("Correlation matrix:");
Console.WriteLine("     X₁    X₂    X₃");
for (int i = 0; i < dim; i++)
{
    Console.Write($"X{i + 1}: ");
    for (int j = 0; j < dim; j++)
    {
        Console.Write($"{corr[i, j],6:F3}");
    }
    Console.WriteLine();
}
```

### Building Covariance from Correlation

```cs
// Given correlation and standard deviations
double[] stdDevs = { 15, 10, 12 };
double[,] corrMatrix = {
    { 1.0, 0.5, 0.3 },
    { 0.5, 1.0, 0.2 },
    { 0.3, 0.2, 1.0 }
};

// Convert to covariance: Σᵢⱼ = ρᵢⱼ·σᵢ·σⱼ
double[,] covMatrix = new double[3, 3];
for (int i = 0; i < 3; i++)
{
    for (int j = 0; j < 3; j++)
    {
        covMatrix[i, j] = corrMatrix[i, j] * stdDevs[i] * stdDevs[j];
    }
}

var mvn = new MultivariateNormal(new double[] { 100, 50, 75 }, covMatrix);

Console.WriteLine("Created MVN from correlation structure");
```

## Practical Applications

### Example 1: Correlated Risk Factors

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

// Model correlated financial returns
// X₁ = Stock A return, X₂ = Stock B return, X₃ = Market return

double[] meanReturns = { 0.08, 0.10, 0.09 };  // 8%, 10%, 9% annual
double[] volatilities = { 0.20, 0.25, 0.15 };  // 20%, 25%, 15% annual

// Correlation structure
double[,] corr = {
    { 1.0, 0.6, 0.8 },  // Stock A vs B: 0.6, A vs Market: 0.8
    { 0.6, 1.0, 0.7 },  // B vs Market: 0.7
    { 0.8, 0.7, 1.0 }
};

// Build covariance matrix
double[,] cov = new double[3, 3];
for (int i = 0; i < 3; i++)
{
    for (int j = 0; j < 3; j++)
    {
        cov[i, j] = corr[i, j] * volatilities[i] * volatilities[j];
    }
}

var returns = new MultivariateNormal(meanReturns, cov);

// Simulate portfolio returns
int nYears = 1000;
double[,] scenarios = returns.GenerateRandomValues(nYears, seed: 123);

// Portfolio: 50% Stock A, 30% Stock B, 20% Market Index
double[] weights = { 0.5, 0.3, 0.2 };

double[] portfolioReturns = new double[nYears];
for (int i = 0; i < nYears; i++)
{
    portfolioReturns[i] = weights[0] * scenarios[i, 0] +
                          weights[1] * scenarios[i, 1] +
                          weights[2] * scenarios[i, 2];
}

Console.WriteLine($"Portfolio Analysis ({nYears} simulations):");
Console.WriteLine($"  Mean return: {portfolioReturns.Average():P2}");
Console.WriteLine($"  Volatility: {Statistics.StandardDeviation(portfolioReturns):P2}");
Console.WriteLine($"  5th percentile: {Statistics.Percentile(portfolioReturns, 0.05):P2}");
Console.WriteLine($"  95th percentile: {Statistics.Percentile(portfolioReturns, 0.95):P2}");
```

### Example 2: Multivariate Exceedance Probability

```cs
// Joint probability of exceeding thresholds
// Example: High temperature AND low rainfall AND high wind speed

double[] means = { 75, 2.5, 15 };  // Temp (°F), Rain (in), Wind (mph)
double[,] cov = {
    { 100, -5, 10 },    // Temp variance=100, negative corr with rain
    { -5, 1.5, -2 },    // Rain variance=1.5, negative corr with wind
    { 10, -2, 25 }      // Wind variance=25
};

var weather = new MultivariateNormal(means, cov);

// Critical thresholds
double[] thresholds = { 85, 1.0, 20 };  // High temp, low rain, high wind

// P(Temp > 85 AND Rain < 1.0 AND Wind > 20)
// This is complex for MVN, but we can estimate with simulation

int nSims = 100000;
double[,] samples = weather.GenerateRandomValues(nSims, seed: 456);

int exceedances = 0;
for (int i = 0; i < nSims; i++)
{
    if (samples[i, 0] > thresholds[0] &&  // High temp
        samples[i, 1] < thresholds[1] &&  // Low rain
        samples[i, 2] > thresholds[2])    // High wind
    {
        exceedances++;
    }
}

double prob = (double)exceedances / nSims;

Console.WriteLine($"Joint Exceedance Probability:");
Console.WriteLine($"  P(Temp>{thresholds[0]}, Rain<{thresholds[1]}, Wind>{thresholds[2]})");
Console.WriteLine($"  = {prob:P4} ({exceedances} out of {nSims})");
Console.WriteLine($"  Return period: {1.0 / prob:F0} events");
```

### Example 3: Conditional Distributions

```cs
// Conditional distribution of Y given X for bivariate normal

double muX = 100, muY = 50;
double sigmaX = 15, sigmaY = 10;
double rho = 0.6;

double[,] cov = {
    { sigmaX * sigmaX, rho * sigmaX * sigmaY },
    { rho * sigmaX * sigmaY, sigmaY * sigmaY }
};

var bivariate = new MultivariateNormal(new[] { muX, muY }, cov);

// Given X = 110, what is distribution of Y?
double xObs = 110;

// Conditional mean: E[Y|X=x] = μ_Y + ρ(σ_Y/σ_X)(x - μ_X)
double muY_given_X = muY + rho * (sigmaY / sigmaX) * (xObs - muX);

// Conditional variance: Var[Y|X=x] = σ_Y²(1 - ρ²)
double sigmaY_given_X = sigmaY * Math.Sqrt(1 - rho * rho);

Console.WriteLine($"Conditional Distribution of Y given X={xObs}:");
Console.WriteLine($"  E[Y|X={xObs}] = {muY_given_X:F2}");
Console.WriteLine($"  SD[Y|X={xObs}] = {sigmaY_given_X:F2}");

// Create conditional distribution
var conditional = new Normal(muY_given_X, sigmaY_given_X);

// 95% prediction interval
double lower = conditional.InverseCDF(0.025);
double upper = conditional.InverseCDF(0.975);

Console.WriteLine($"  95% prediction interval: [{lower:F2}, {upper:F2}]");
```

### Example 4: Mahalanobis Distance

```cs
// Measure distance from mean accounting for correlations

double[] means = { 100, 50, 75 };
double[,] cov = {
    { 225, 75, 50 },
    { 75, 100, 25 },
    { 50, 25, 144 }
};

var mvn = new MultivariateNormal(means, cov);

// Test point
double[] x = { 120, 60, 90 };

// Mahalanobis distance: D² = (x-μ)ᵀΣ⁻¹(x-μ)
var diff = new Vector(x.Select((xi, i) => xi - means[i]).ToArray());
var covMatrix = new Matrix(cov);
var covInv = covMatrix.Inverse();

double mahalanobis = Math.Sqrt(Vector.DotProduct(diff, covInv.Multiply(diff)));

// Chi-squared test: D² ~ χ²(k) under null hypothesis
double chiSqCritical = 7.815;  // 95th percentile of χ²(3)

Console.WriteLine($"Mahalanobis Distance Analysis:");
Console.WriteLine($"  Point: [{string.Join(", ", x)}]");
Console.WriteLine($"  Distance: {mahalanobis:F3}");
Console.WriteLine($"  D²: {mahalanobis * mahalanobis:F3}");
Console.WriteLine($"  Critical value (95%): {chiSqCritical:F3}");

if (mahalanobis * mahalanobis > chiSqCritical)
    Console.WriteLine("  → Point is an outlier (p < 0.05)");
else
    Console.WriteLine("  → Point is not an outlier");
```

## Multivariate Student-t Distribution

The multivariate Student-t distribution generalizes the univariate Student-t to multiple dimensions, providing heavier tails than the multivariate normal. It is used for robust modeling when data may contain outliers [[2]](#2).

### Mathematical Definition

A random vector $\mathbf{X} = (X_1, \ldots, X_k)^T$ follows a $k$-dimensional multivariate Student-t distribution, written $\mathbf{X} \sim t_k(\nu, \boldsymbol{\mu}, \boldsymbol{\Sigma})$, if its probability density function is:

```math
f(\mathbf{x}) = \frac{\Gamma\!\left(\frac{\nu + k}{2}\right)}{\Gamma\!\left(\frac{\nu}{2}\right) (\nu\pi)^{k/2} |\boldsymbol{\Sigma}|^{1/2}} \left[ 1 + \frac{1}{\nu} (\mathbf{x} - \boldsymbol{\mu})^T \boldsymbol{\Sigma}^{-1} (\mathbf{x} - \boldsymbol{\mu}) \right]^{-(\nu+k)/2}
```

where:
- $\nu > 0$ is the degrees of freedom,
- $\boldsymbol{\mu} \in \mathbb{R}^k$ is the location vector,
- $\boldsymbol{\Sigma}$ is the $k \times k$ positive-definite scale matrix (not the covariance matrix).

**Important:** The scale matrix $\boldsymbol{\Sigma}$ is distinct from the covariance matrix. The moments are:

```math
\text{Mean} = \boldsymbol{\mu} \quad (\text{for } \nu > 1), \qquad \text{Cov}(\mathbf{X}) = \frac{\nu}{\nu - 2}\,\boldsymbol{\Sigma} \quad (\text{for } \nu > 2)
```

The mean is undefined when $\nu \leq 1$, and the covariance is undefined when $\nu \leq 2$.

**Relationship to MVN.** As $\nu \to \infty$, the multivariate Student-t distribution converges to the multivariate normal $\mathcal{N}_k(\boldsymbol{\mu}, \boldsymbol{\Sigma})$. For finite $\nu$, the tails are heavier, controlled by the degrees of freedom parameter. This makes the multivariate Student-t distribution useful for robust statistical inference where normal approximations may underestimate tail probabilities.

**Sampling algorithm.** The multivariate Student-t can be represented as a scale mixture of normals. A sample is generated as:

```math
\mathbf{x} = \boldsymbol{\mu} + \frac{\mathbf{L}\mathbf{z}}{\sqrt{W/\nu}}
```

where $\mathbf{L}$ is the Cholesky factor of $\boldsymbol{\Sigma}$, $\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_k)$, and $W \sim \chi^2(\nu)$. The library generates $W$ via $W \sim \text{Gamma}(\nu/2, 2)$ to support non-integer degrees of freedom.

### Creating Multivariate Student-t Distributions

```cs
using Numerics.Distributions;

// Bivariate Student-t with 5 degrees of freedom
double nu = 5;
double[] mu = { 100, 50 };
double[,] sigma = {
    { 225, 75 },
    { 75, 100 }
};

var mvt = new MultivariateStudentT(nu, mu, sigma);

Console.WriteLine($"Multivariate Student-t Distribution:");
Console.WriteLine($"  Dimension: {mvt.Dimension}");
Console.WriteLine($"  Degrees of freedom: {nu}");
```

### PDF, CDF, and Sampling

```cs
// Evaluate density
double[] x = { 105, 55 };
double pdf = mvt.PDF(x);
double logPdf = mvt.LogPDF(x);

Console.WriteLine($"PDF at x: {pdf:E4}");
Console.WriteLine($"Log-PDF at x: {logPdf:F4}");

// CDF
double cdf = mvt.CDF(x);
Console.WriteLine($"CDF at x: {cdf:F6}");

// Generate random samples
double[,] samples = mvt.GenerateRandomValues(1000, seed: 42);
```

## Dirichlet Distribution

The Dirichlet distribution $\text{Dir}(\alpha_1, \alpha_2, \ldots, \alpha_k)$ is a multivariate generalization of the Beta distribution, defined on the probability simplex (components sum to 1). It is the conjugate prior for the categorical and multinomial distributions in Bayesian statistics [[3]](#3).

### Mathematical Definition

A random vector $\mathbf{X} = (X_1, X_2, \ldots, X_k)$ follows a Dirichlet distribution $\text{Dir}(\alpha_1, \ldots, \alpha_k)$ if its probability density function on the $(k-1)$-dimensional simplex is:

```math
f(x_1, \ldots, x_k) = \frac{\Gamma\!\left(\sum_{i=1}^{k} \alpha_i\right)}{\prod_{i=1}^{k} \Gamma(\alpha_i)} \prod_{i=1}^{k} x_i^{\alpha_i - 1}
```

with constraints $x_i > 0$ for all $i$ and $\sum_{i=1}^{k} x_i = 1$. Each $\alpha_i > 0$ is a concentration parameter. The normalizing constant involves the multivariate Beta function $B(\boldsymbol{\alpha}) = \prod \Gamma(\alpha_i) / \Gamma(\sum \alpha_i)$.

**Moments.** Let $\alpha_0 = \sum_{i=1}^{k} \alpha_i$ denote the total concentration. Then:

```math
\text{E}[X_i] = \frac{\alpha_i}{\alpha_0}, \qquad \text{Var}(X_i) = \frac{\alpha_i (\alpha_0 - \alpha_i)}{\alpha_0^2 (\alpha_0 + 1)}
```

```math
\text{Cov}(X_i, X_j) = \frac{-\alpha_i \alpha_j}{\alpha_0^2 (\alpha_0 + 1)} \quad (i \neq j)
```

The negative covariance reflects that components are constrained to sum to 1 -- if one component increases, others must decrease. The mode exists when all $\alpha_i > 1$ and is given by $\text{Mode}(X_i) = (\alpha_i - 1) / (\alpha_0 - k)$.

**Connection to Beta distribution.** For $k = 2$, the Dirichlet reduces to the Beta distribution: if $(X_1, X_2) \sim \text{Dir}(\alpha_1, \alpha_2)$, then $X_1 \sim \text{Beta}(\alpha_1, \alpha_2)$.

**Sampling algorithm.** Samples are generated using the gamma distribution representation. Draw $k$ independent gamma variates $Y_i \sim \text{Gamma}(\alpha_i, 1)$, then normalize:

```math
X_i = \frac{Y_i}{\sum_{j=1}^{k} Y_j}
```

The resulting vector $(X_1, \ldots, X_k)$ lies on the simplex and follows the Dirichlet distribution.

**Use cases.** The Dirichlet distribution arises naturally when modeling proportions or mixture weights: Bayesian priors for categorical data, topic modeling (document-topic proportions), compositional data analysis, and mixture model weight estimation in RMC-BestFit.

### Creating Dirichlet Distributions

```cs
using Numerics.Distributions;

// Symmetric Dirichlet with K=3 categories, all alpha=2
var dir1 = new Dirichlet(dimension: 3, alpha: 2.0);

// Asymmetric Dirichlet with specified concentration parameters
var dir2 = new Dirichlet(new double[] { 1.0, 2.0, 5.0 });

Console.WriteLine($"Dirichlet Distribution:");
Console.WriteLine($"  Dimension: {dir2.Dimension}");
Console.WriteLine($"  Alpha sum: {dir2.AlphaSum}");
```

### Properties

```cs
var dir = new Dirichlet(new double[] { 2.0, 3.0, 5.0 });

// Mean vector: E[Xi] = alpha_i / sum(alpha)
double[] mean = dir.Mean;
Console.WriteLine($"Mean: [{string.Join(", ", mean.Select(m => m.ToString("F3")))}]");

// Variance vector
double[] variance = dir.Variance;
Console.WriteLine($"Variance: [{string.Join(", ", variance.Select(v => v.ToString("F4")))}]");

// Mode (requires all alpha_i > 1)
double[] mode = dir.Mode;
Console.WriteLine($"Mode: [{string.Join(", ", mode.Select(m => m.ToString("F3")))}]");

// Covariance between components
double cov01 = dir.Covariance(0, 1);  // Negative (components are negatively correlated)
Console.WriteLine($"Cov(X0, X1): {cov01:F4}");
```

### PDF and Sampling

```cs
// Evaluate density at a point on the simplex
double[] x = { 0.2, 0.3, 0.5 };
double pdf = dir.PDF(x);
double logPdf = dir.LogPDF(x);

Console.WriteLine($"PDF at x: {pdf:F6}");
Console.WriteLine($"Log-PDF at x: {logPdf:F4}");

// Generate random samples (each row sums to 1)
double[,] samples = dir.GenerateRandomValues(1000, seed: 42);
```

### Application: Bayesian Mixture Weights

```cs
// Prior for 3-component mixture model weights
var prior = new Dirichlet(dimension: 3, alpha: 1.0);  // Uniform on simplex

// After observing data, update to posterior
// (conjugate update: add counts to alpha)
int[] counts = { 50, 120, 80 };
double[] posteriorAlpha = { 1.0 + counts[0], 1.0 + counts[1], 1.0 + counts[2] };
var posterior = new Dirichlet(posteriorAlpha);

Console.WriteLine($"Posterior mean weights: [{string.Join(", ", posterior.Mean.Select(m => m.ToString("F3")))}]");
```

## Multinomial Distribution

The Multinomial distribution models the number of outcomes in each of $k$ categories over $n$ independent trials, where each trial has a fixed probability vector. It generalizes the Binomial distribution to more than two categories [[4]](#4).

### Mathematical Definition

A random vector $\mathbf{X} = (X_1, X_2, \ldots, X_k)$ follows a multinomial distribution $\text{Mult}(n, \mathbf{p})$ with $n$ trials and probability vector $\mathbf{p} = (p_1, \ldots, p_k)$ if its probability mass function is:

```math
P(X_1 = x_1, \ldots, X_k = x_k) = \frac{n!}{\prod_{i=1}^{k} x_i!} \prod_{i=1}^{k} p_i^{x_i}
```

with constraints $x_i \in \lbrace 0, 1, \ldots, n\rbrace$, $\sum_{i=1}^{k} x_i = n$, $p_i \geq 0$, and $\sum_{i=1}^{k} p_i = 1$.

The multinomial coefficient $n! / \prod x_i!$ counts the number of ways to arrange $n$ trials into $k$ categories with the specified counts.

**Moments:**

```math
\text{E}[X_i] = n p_i, \qquad \text{Var}(X_i) = n p_i (1 - p_i)
```

```math
\text{Cov}(X_i, X_j) = -n p_i p_j \quad (i \neq j)
```

The negative covariance follows from the constraint $\sum X_i = n$: observing more counts in one category necessarily reduces the counts available for other categories.

**Connection to Binomial.** For $k = 2$, the multinomial reduces to the binomial distribution: if $(X_1, X_2) \sim \text{Mult}(n, (p, 1-p))$, then $X_1 \sim \text{Binomial}(n, p)$.

**Sampling algorithm.** The library generates samples via sequential conditional binomial draws. For each category $i = 1, \ldots, k-1$, draw $X_i \sim \text{Binomial}(n_{\text{remaining}}, p_i / p_{\text{remaining}})$, then set the last category to the remainder. This method is exact and efficient.

**Use cases.** The multinomial distribution is used for modeling count data across categories: MCMC trajectory state selection in the NUTS sampler, categorical data modeling in LifeSim, and as the likelihood function for categorical observations in Bayesian inference.

### Creating Multinomial Distributions

```cs
using Numerics.Distributions;

// 10 trials with 3 categories having probabilities 0.2, 0.3, 0.5
var mult = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });

Console.WriteLine($"Multinomial Distribution:");
Console.WriteLine($"  Dimension: {mult.Dimension}");
Console.WriteLine($"  Number of trials: {mult.NumberOfTrials}");
```

### Properties

```cs
// Mean vector: E[Xi] = N * p_i
double[] mean = mult.Mean;
Console.WriteLine($"Mean: [{string.Join(", ", mean.Select(m => m.ToString("F1")))}]");

// Variance vector: Var[Xi] = N * p_i * (1 - p_i)
double[] variance = mult.Variance;
Console.WriteLine($"Variance: [{string.Join(", ", variance.Select(v => v.ToString("F2")))}]");

// Covariance: Cov(Xi, Xj) = -N * p_i * p_j
double cov01 = mult.Covariance(0, 1);
Console.WriteLine($"Cov(X0, X1): {cov01:F2}");
```

### PMF and Sampling

```cs
// Probability mass function
double[] counts = { 2, 3, 5 };
double pmf = mult.PDF(counts);
double logPmf = mult.LogPMF(counts);

Console.WriteLine($"P(X = [2,3,5]): {pmf:F6}");

// Generate random samples (each row sums to N)
double[,] samples = mult.GenerateRandomValues(1000, seed: 42);

// Weighted categorical sampling (static utility)
var rng = new Random(42);
double[] weights = { 1.0, 3.0, 6.0 };
int category = Multinomial.Sample(weights, rng);
Console.WriteLine($"Sampled category: {category}");
```

## Key Concepts

### Positive Definiteness

The covariance matrix must be positive definite:

```cs
// Valid covariance matrix
double[,] validCov = {
    { 4.0, 1.5 },
    { 1.5, 2.0 }
};

// Invalid - not positive definite
double[,] invalidCov = {
    { 1.0, 2.0 },
    { 2.0, 1.0 }  // Correlation would be > 1
};

try
{
    var invalid = new MultivariateNormal(new[] { 0.0, 0.0 }, invalidCov);
}
catch (Exception ex)
{
    Console.WriteLine($"Error: {ex.Message}");
}
```

### Marginal Distributions

Each component follows univariate normal:

```cs
var mvn = new MultivariateNormal(means, cov);

// Marginal distribution of X₁
double mu1 = mvn.Mean[0];
double sigma1 = Math.Sqrt(mvn.Covariance[0, 0]);
var marginal1 = new Normal(mu1, sigma1);

Console.WriteLine($"Marginal X₁ ~ N({mu1:F2}, {sigma1:F2})");
```

### Independence

Variables are independent if covariance matrix is diagonal:

```cs
// Independent variables
double[,] indepCov = {
    { 100, 0, 0 },
    { 0, 64, 0 },
    { 0, 0, 225 }
};

var independent = new MultivariateNormal(new[] { 0.0, 0.0, 0.0 }, indepCov);
Console.WriteLine("Independent MVN → Diagonal covariance matrix");
```

## Best Practices

1. **Validate covariance matrix** - Must be symmetric and positive definite
2. **Check correlations** - Ensure |ρ| ≤ 1 for all pairs
3. **Use Cholesky decomposition** - Efficient for generation and PDF
4. **Consider dimension** - CDF computation expensive for high dimensions
5. **Monte Carlo for CDF** - Set appropriate tolerance for accuracy
6. **Conditional distributions** - Use formula for bivariate case
7. **Outlier detection** - Use Mahalanobis distance

## Sampling Algorithms

This section summarizes the algorithms used by the library to generate random samples from each multivariate distribution.

### MVN via Cholesky Decomposition

The `MultivariateNormal.GenerateRandomValues` method uses the Cholesky decomposition $\boldsymbol{\Sigma} = \mathbf{L}\mathbf{L}^T$ to transform independent standard normal variates into correlated samples:

```math
\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_k), \quad \mathbf{x} = \boldsymbol{\mu} + \mathbf{L}\mathbf{z}
```

The Cholesky approach is preferred for three reasons: (1) it is computationally efficient -- the decomposition is $O(k^3)$ but computed only once, after which each sample requires only $O(k^2)$; (2) the decomposition succeeds if and only if $\boldsymbol{\Sigma}$ is positive definite, providing a built-in validity check; and (3) the triangular structure of $\mathbf{L}$ makes the matrix-vector product efficient. The library also provides `LatinHypercubeRandomValues` for improved space-filling properties via stratified sampling.

### Multivariate Student-t Sampling

The `MultivariateStudentT.GenerateRandomValues` method exploits the scale-mixture representation. A multivariate Student-t variate can be constructed by scaling a multivariate normal variate by an independent chi-squared random variable:

```math
\mathbf{z} \sim \mathcal{N}(\mathbf{0}, \mathbf{I}_k), \quad W \sim \chi^2(\nu), \quad \mathbf{x} = \boldsymbol{\mu} + \mathbf{L}\mathbf{z} \cdot \sqrt{\frac{\nu}{W}}
```

The library generates $W$ via $W \sim \text{Gamma}(\nu/2, 2)$ to support non-integer degrees of freedom.

### Dirichlet Sampling

The `Dirichlet.GenerateRandomValues` method uses the gamma distribution representation. Draw $k$ independent gamma variates and normalize:

```math
Y_i \sim \text{Gamma}(\alpha_i, 1), \quad X_i = \frac{Y_i}{\sum_{j=1}^{k} Y_j}
```

### Multinomial Sampling

The `Multinomial.GenerateRandomValues` method uses sequential conditional binomial draws. For category $i = 1, \ldots, k-1$:

```math
X_i \sim \text{Binomial}\!\left(n_{\text{remaining}},\; \frac{p_i}{p_{\text{remaining}}}\right), \quad X_k = n - \sum_{i=1}^{k-1} X_i
```

This method is exact and avoids the need to enumerate the full multinomial support. The library uses direct simulation for small $n$ and the BTPE algorithm for large $n$.

## Computational Notes

- **PDF**: O(k^3) for Cholesky decomposition (one-time), O(k^2) per evaluation
- **CDF**: O(k) for k<=2, O(n*k) for k>2 (Monte Carlo with n evaluations)
- **Sampling**: O(k^2) per sample using Cholesky (MVN and Student-t)
- **Storage**: O(k^2) for covariance/scale matrix

---

## Bivariate Empirical Distribution

The `BivariateEmpirical` class represents a bivariate empirical CDF defined on a grid of values. It is useful when the joint distribution is known only through tabulated CDF values rather than a parametric form.

```cs
using Numerics.Distributions;
using Numerics.Data;

// Define grid values
double[] x1Values = { 100, 200, 300, 400, 500 };  // e.g., flow (cfs)
double[] x2Values = { 10, 20, 30, 40 };            // e.g., stage (ft)

// Define CDF values on the grid (must be non-decreasing, values in [0, 1])
double[,] pValues = {
    { 0.01, 0.02, 0.03, 0.04 },
    { 0.05, 0.15, 0.20, 0.22 },
    { 0.10, 0.35, 0.50, 0.55 },
    { 0.15, 0.50, 0.75, 0.82 },
    { 0.20, 0.60, 0.90, 1.00 }
};

// Create bivariate empirical distribution
var bivariate = new BivariateEmpirical(x1Values, x2Values, pValues);

// Evaluate CDF at a point
double prob = bivariate.CDF(250, 25);
Console.WriteLine($"P(X1 ≤ 250, X2 ≤ 25) = {prob:F4}");

// Optional: apply transforms for interpolation
var bivariateLog = new BivariateEmpirical(
    x1Values, x2Values, pValues,
    x1Transform: Transform.Logarithmic,
    x2Transform: Transform.None,
    probabilityTransform: Transform.NormalZ
);
```

**Properties:**
- `X1Values`, `X2Values` — grid axis values
- `ProbabilityValues` — CDF values on the grid
- `Dimension` — always 2
- Supports `Transform.None`, `Transform.Logarithmic`, and `Transform.NormalZ` for interpolation

**Note:** The `PDF()` method returns `double.NaN` because this is a purely empirical CDF — no density function is defined.

### Practical Example: Joint Flood Stage-Duration Exceedance

```cs
using Numerics.Distributions;
using Numerics.Data;

// Define flood stage grid (ft)
double[] stages = { 10, 12, 14, 16, 18, 20 };
// Define flood duration grid (hr)
double[] durations = { 6, 12, 24, 48, 72 };

// Joint CDF values P(Stage ≤ s, Duration ≤ d)
// Based on observed flood event records
double[,] jointCDF = {
    { 0.30, 0.35, 0.38, 0.39, 0.40 },
    { 0.40, 0.50, 0.55, 0.57, 0.58 },
    { 0.45, 0.60, 0.70, 0.73, 0.75 },
    { 0.48, 0.65, 0.78, 0.83, 0.85 },
    { 0.49, 0.68, 0.82, 0.90, 0.93 },
    { 0.50, 0.70, 0.85, 0.95, 1.00 }
};

// Create with log transform on stages for better interpolation
var jointDist = new BivariateEmpirical(
    stages, durations, jointCDF,
    x1Transform: Transform.Logarithmic,
    x2Transform: Transform.None,
    probabilityTransform: Transform.NormalZ
);

// Probability that a flood has stage ≤ 15 ft AND duration ≤ 36 hr
double prob = jointDist.CDF(15, 36);
Console.WriteLine($"P(Stage ≤ 15, Duration ≤ 36hr) = {prob:F4}");

// Joint exceedance probability
double jointExceedance = 1 - jointDist.CDF(16, 48);
Console.WriteLine($"P(Stage > 16 AND/OR Duration > 48hr) = {jointExceedance:F4}");
```

---

## References

<a id="1">[1]</a> Anderson, T. W. (2003). *An Introduction to Multivariate Statistical Analysis* (3rd ed.). Wiley.

<a id="2">[2]</a> Kotz, S., Balakrishnan, N., & Johnson, N. L. (2000). *Continuous Multivariate Distributions, Volume 1: Models and Applications* (2nd ed.). Wiley.

<a id="3">[3]</a> Kotz, S., Balakrishnan, N. & Johnson, N. L. (2000). *Continuous Multivariate Distributions, Volume 1: Models and Applications* (2nd ed.). Wiley. Chapter 49 (Dirichlet distribution).

<a id="4">[4]</a> Johnson, N. L., Kotz, S. & Balakrishnan, N. (1997). *Discrete Multivariate Distributions*. Wiley. Chapter 35 (Multinomial distribution).

<a id="5">[5]</a> Tong, Y. L. (2012). *The Multivariate Normal Distribution*. Springer Science & Business Media.

<a id="6">[6]</a> Kotz, S. & Nadarajah, S. (2004). *Multivariate t Distributions and Their Applications*. Cambridge University Press.

---

[← Previous: Copulas](copulas.md) | [Back to Index](../index.md) | [Next: Machine Learning →](../machine-learning/machine-learning.md)
