# Multivariate Distributions

[← Previous: Copulas](copulas.md) | [Back to Index](../index.md)

The ***Numerics*** library provides the **Multivariate Normal** distribution for modeling correlated random variables. This distribution is fundamental in multivariate statistics, risk assessment, and uncertainty quantification.

## Multivariate Normal Distribution

The multivariate normal (Gaussian) distribution generalizes the univariate normal to multiple dimensions with a specified covariance structure [[1]](#1).

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

**Formula:**
```
f(x) = (2π)^(-k/2) |Σ|^(-1/2) exp[-½(x-μ)ᵀΣ⁻¹(x-μ)]

Where:
- k = dimension
- μ = mean vector
- Σ = covariance matrix
- |Σ| = determinant of Σ
```

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

```cs
// Generate quantiles for each marginal
double[] probabilities = { 0.05, 0.50, 0.95 };

double[] quantiles = mvn.InverseCDF(probabilities);

Console.WriteLine("Marginal quantiles:");
for (int i = 0; i < mvn.Dimension; i++)
{
    double q05 = mvn.Mean[i] + Math.Sqrt(mvn.Covariance[i, i]) * 
                 new Normal(0, 1).InverseCDF(0.05);
    double q50 = mvn.Mean[i];
    double q95 = mvn.Mean[i] + Math.Sqrt(mvn.Covariance[i, i]) * 
                 new Normal(0, 1).InverseCDF(0.95);
    
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
Console.WriteLine($"  5th percentile: {Statistics.Percentile(portfolioReturns.OrderBy(r => r).ToArray(), 5):P2}");
Console.WriteLine($"  95th percentile: {Statistics.Percentile(portfolioReturns.OrderBy(r => r).ToArray(), 95):P2}");
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

double mahalanobis = Math.Sqrt(diff.DotProduct(covInv.Multiply(diff)));

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

## Computational Notes

- **PDF**: O(k³) for Cholesky decomposition
- **CDF**: O(k) for k≤2, O(n·k) for k>2 (Monte Carlo with n evaluations)
- **Sampling**: O(k²) per sample using Cholesky
- **Storage**: O(k²) for covariance matrix

---

## References

<a id="1">[1]</a> Tong, Y. L. (2012). *The Multivariate Normal Distribution*. Springer Science & Business Media.

<a id="2">[2]</a> Kotz, S., Balakrishnan, N., & Johnson, N. L. (2004). *Continuous Multivariate Distributions, Volume 1: Models and Applications* (2nd ed.). John Wiley & Sons.

---

[← Previous: Copulas](copulas.md) | [Back to Index](../index.md)
