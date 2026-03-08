# Linear Regression

[← Previous: Interpolation](interpolation.md) | [Back to Index](../index.md) | [Next: Time Series →](time-series.md)

Linear regression estimates the relationship between a scalar response variable $y$ and one or more predictor variables $\mathbf{x}$. The model assumes a linear form:

```math
y_i = \beta_0 + \beta_1 x_{i1} + \beta_2 x_{i2} + \cdots + \beta_p x_{ip} + \varepsilon_i
```

where $\beta_0$ is the intercept, $\beta_1, \ldots, \beta_p$ are the regression coefficients, and $\varepsilon_i \sim N(0, \sigma^2)$ are independent error terms. In matrix notation, this becomes $\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \boldsymbol{\varepsilon}$.

## The Normal Equations

The ordinary least squares (OLS) estimator minimizes the sum of squared residuals $\sum_{i=1}^n (y_i - \hat{y}_i)^2$. Setting the gradient of this objective to zero yields the **normal equations**:

```math
\mathbf{X}^T \mathbf{X} \boldsymbol{\hat{\beta}} = \mathbf{X}^T \mathbf{y}
```

When $\mathbf{X}^T \mathbf{X}$ is invertible, the solution is:

```math
\boldsymbol{\hat{\beta}} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \mathbf{y}
```

However, directly computing $(\mathbf{X}^T \mathbf{X})^{-1}$ is numerically ill-conditioned, especially when predictor variables are correlated or span very different scales. The condition number of $\mathbf{X}^T \mathbf{X}$ is the square of the condition number of $\mathbf{X}$, so even moderate collinearity can cause significant loss of precision.

### Why SVD?

The ***Numerics*** library uses **Singular Value Decomposition** (SVD) to solve the least squares problem instead [[1]](#1). The SVD factorizes the design matrix as $\mathbf{X} = \mathbf{U} \mathbf{W} \mathbf{V}^T$, where $\mathbf{U}$ and $\mathbf{V}$ are orthogonal matrices and $\mathbf{W}$ is a diagonal matrix of singular values. The least squares solution is then:

```math
\boldsymbol{\hat{\beta}} = \mathbf{V} \mathbf{W}^{-1} \mathbf{U}^T \mathbf{y}
```

SVD provides several advantages:
- **Numerical stability**: Works correctly even when $\mathbf{X}^T \mathbf{X}$ is nearly singular
- **Rank detection**: Small singular values (below a threshold of $10^{-12}$) are set to zero, effectively handling rank-deficient problems
- **Covariance computation**: The parameter covariance matrix is computed directly from the SVD components as $\text{Cov}(\hat{\beta}_{i}, \hat{\beta}_{j}) = \sigma^2 \sum_k \frac{V_{ik} V_{jk}}{W_k^2}$

## Creating a Linear Regression Model

The model is fitted automatically in the constructor — no separate `Train()` call is needed:

```cs
using Numerics.Data;
using Numerics.Mathematics.LinearAlgebra;

// Predictor data (no intercept column needed — added automatically)
double[,] X = {
    { 1.0 },
    { 2.0 },
    { 3.0 },
    { 4.0 },
    { 5.0 }
};

double[] y = { 2.1, 4.0, 5.8, 8.1, 9.9 };

// Create and fit the model (fitting happens in the constructor)
var lm = new LinearRegression(new Matrix(X), new Vector(y));

Console.WriteLine($"Intercept: {lm.Parameters[0]:F4}");
Console.WriteLine($"Slope: {lm.Parameters[1]:F4}");
Console.WriteLine($"R²: {lm.RSquared:F4}");
Console.WriteLine($"Adjusted R²: {lm.AdjRSquared:F4}");
Console.WriteLine($"Standard Error: {lm.StandardError:F4}");
```

## Model Properties

After construction, the following properties are available:

| Property | Type | Description |
|----------|------|-------------|
| `Parameters` | `List<double>` | Estimated coefficients (intercept first if `HasIntercept = true`) |
| `ParameterNames` | `List<string>` | Names of each parameter |
| `ParameterStandardErrors` | `List<double>` | Standard errors of each coefficient |
| `ParameterTStats` | `List<double>` | t-statistics for each coefficient |
| `Covariance` | `Matrix` | Parameter covariance matrix |
| `Residuals` | `double[]` | Model residuals (observed − predicted) |
| `StandardError` | `double` | Residual standard error |
| `SampleSize` | `int` | Number of observations |
| `DegreesOfFreedom` | `int` | Residual degrees of freedom ($n - p$) |
| `RSquared` | `double` | Coefficient of determination |
| `AdjRSquared` | `double` | Adjusted $R^2$ |
| `HasIntercept` | `bool` | Whether the model includes an intercept |

### Understanding Key Statistics

**Coefficient of determination** ($R^2$) measures the proportion of variance in $y$ explained by the model:

```math
R^2 = 1 - \frac{\sum (y_i - \hat{y}_i)^2}{\sum (y_i - \bar{y})^2} = 1 - \frac{SS_{res}}{SS_{tot}}
```

**Adjusted $R^2$** penalizes for additional predictors, preventing overfitting:

```math
R^2_{adj} = 1 - \frac{SS_{res} / (n - p)}{SS_{tot} / (n - 1)}
```

where $n$ is the sample size and $p$ is the number of parameters (including intercept).

**Parameter standard errors** quantify uncertainty in each coefficient estimate. The standard error of $\hat{\beta}_j$ is:

```math
SE(\hat{\beta}_j) = \hat{\sigma} \sqrt{C_{jj}}
```

where $\hat{\sigma}$ is the residual standard error and $C_{jj}$ is the $j$-th diagonal element of $(\mathbf{X}^T \mathbf{X})^{-1}$ (computed via SVD in practice).

**t-statistics** test whether each coefficient is significantly different from zero: $t_j = \hat{\beta}_j / SE(\hat{\beta}_j)$. Under the null hypothesis $\beta_j = 0$, this follows a Student's t-distribution with $n - p$ degrees of freedom.

## Summary Output

The `Summary()` method produces an R-style summary table:

```cs
using Numerics.Data;
using Numerics.Mathematics.LinearAlgebra;

double[,] X = {
    { 12.5 }, { 25.0 }, { 48.3 }, { 75.0 }, { 102.0 },
    { 18.7 }, { 55.0 }, { 130.0 }, { 8.2 }, { 200.0 }
};

double[] y = { 1450, 2680, 4520, 6100, 7850, 2050, 5200, 9300, 980, 12500 };

var lm = new LinearRegression(new Matrix(X), new Vector(y));

// Print R-style summary
foreach (var line in lm.Summary())
{
    Console.WriteLine(line);
}
```

The summary includes:
- Coefficient estimates with standard errors, t-statistics, and p-values
- Significance codes (*** < 0.001, ** < 0.01, * < 0.05)
- Residual standard error and degrees of freedom
- $R^2$ and adjusted $R^2$
- F-statistic with p-value
- Five-number summary of residuals

## Prediction

### Point Predictions

```cs
// Predict for new predictor values
double[,] XNew = {
    { 50.0 },
    { 100.0 },
    { 150.0 }
};

double[] predictions = lm.Predict(new Matrix(XNew));

Console.WriteLine("Predictions:");
for (int i = 0; i < predictions.Length; i++)
{
    Console.WriteLine($"  X = {XNew[i, 0]:F1} → ŷ = {predictions[i]:F1}");
}
```

### Prediction Intervals

Prediction intervals account for both the uncertainty in the estimated regression line and the inherent variability of individual observations:

```math
\hat{y}_0 \pm t_{\alpha/2, n-p} \cdot \hat{\sigma} \sqrt{1 + \mathbf{x}_0^T (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{x}_0}
```

The $1$ inside the square root represents the irreducible variance of a new observation; the second term represents uncertainty in the estimated mean.

```cs
// 90% prediction intervals (alpha = 0.1)
double[,] intervals = lm.PredictionIntervals(new Matrix(XNew), alpha: 0.1);

Console.WriteLine("90% Prediction Intervals:");
Console.WriteLine("    X      |   Lower   |   Mean    |   Upper");
Console.WriteLine("-----------|-----------|-----------|----------");

for (int i = 0; i < XNew.GetLength(0); i++)
{
    Console.WriteLine($"  {XNew[i, 0],7:F1} | {intervals[i, 0],9:F1} | {intervals[i, 2],9:F1} | {intervals[i, 1],9:F1}");
}
```

**Note:** `PredictionIntervals` returns a 2D array with columns: lower (index 0), upper (index 1), mean (index 2).

## Multiple Regression

Multiple predictor variables are supported:

```cs
using Numerics.Data;
using Numerics.Mathematics.LinearAlgebra;

// Three predictor variables
double[,] X = {
    { 12.5, 38.2, 820 },
    { 25.0, 40.1, 790 },
    { 48.3, 39.5, 850 },
    { 75.0, 41.0, 780 },
    { 102.0, 42.3, 760 },
    { 130.0, 43.0, 740 },
    { 200.0, 44.2, 720 }
};

double[] y = { 1450, 2680, 4520, 6100, 7850, 9300, 12500 };

var lm = new LinearRegression(new Matrix(X), new Vector(y));

Console.WriteLine("Multiple Regression Results:");
for (int i = 0; i < lm.Parameters.Count; i++)
{
    Console.WriteLine($"  {lm.ParameterNames[i]}: {lm.Parameters[i]:F4} (SE: {lm.ParameterStandardErrors[i]:F4})");
}

Console.WriteLine($"\nR²: {lm.RSquared:F4}");
Console.WriteLine($"Adjusted R²: {lm.AdjRSquared:F4}");
```

### Multicollinearity

When predictor variables are highly correlated, the regression coefficients become unstable — small changes in the data lead to large changes in the estimated coefficients, and standard errors inflate. This is called **multicollinearity**.

Signs of multicollinearity:
- Large standard errors on coefficients that should be significant
- $R^2$ is high but individual t-statistics are low
- Coefficient signs are opposite to what domain knowledge suggests

One diagnostic is the **Variance Inflation Factor** (VIF), defined for each predictor $j$ as:

```math
\text{VIF}_j = \frac{1}{1 - R^2_j}
```

where $R^2_j$ is the $R^2$ from regressing $x_j$ on all other predictors. A VIF above 10 suggests problematic collinearity. While the ***Numerics*** library does not compute VIF directly, you can fit auxiliary regressions to diagnose it.

## Regression Without Intercept

```cs
// Force model through the origin (no intercept)
var lm = new LinearRegression(new Matrix(X), new Vector(y), hasIntercept: false);

Console.WriteLine("No-intercept model:");
Console.WriteLine($"Slope: {lm.Parameters[0]:F4}");
```

**Caution**: Removing the intercept changes the interpretation of $R^2$ and can bias coefficients if the true relationship does not pass through the origin. Only use this when the physics of the problem demands it.

## Regression Assumptions

The validity of OLS inference (standard errors, p-values, prediction intervals) depends on several assumptions about the error terms $\varepsilon_i$:

1. **Linearity**: The relationship between $\mathbf{x}$ and $y$ is linear. Check by plotting residuals vs. fitted values — a systematic pattern indicates nonlinearity.

2. **Independence**: Errors are independent. Violated with time series data or spatial data. The Ljung-Box test (available in ***Numerics*** via `HypothesisTests.LjungBoxTest`) can detect serial correlation.

3. **Homoscedasticity**: Errors have constant variance $\text{Var}(\varepsilon_i) = \sigma^2$. A fan-shaped pattern in the residual plot indicates heteroscedasticity. Consider log-transforming the response variable if variance increases with the mean.

4. **Normality**: Errors are normally distributed. Check with the Jarque-Bera test (available via `HypothesisTests.JarqueBeraTest`). Mild departures from normality have little impact on coefficient estimates but affect confidence intervals and p-values.

## Practical Examples

### Example 1: Regional Streamflow Regression

Predicting annual peak streamflow from watershed characteristics using data from [`example-data/streamflow-regression.csv`](../example-data/streamflow-regression.csv):

```cs
using System.IO;
using System.Linq;
using Numerics.Data;
using Numerics.Mathematics.LinearAlgebra;

// Load CSV data (skip comment lines starting with #)
string[] lines = File.ReadAllLines("example-data/streamflow-regression.csv");
var dataLines = lines
    .Where(line => !line.StartsWith("#") && !string.IsNullOrWhiteSpace(line))
    .Skip(1) // Skip header
    .ToArray();

int n = dataLines.Length;
double[,] features = new double[n, 3];
double[] flow = new double[n];

for (int i = 0; i < n; i++)
{
    var parts = dataLines[i].Split(',');
    features[i, 0] = double.Parse(parts[0]); // DrainageArea_sqmi
    features[i, 1] = double.Parse(parts[1]); // MeanAnnualPrecip_in
    features[i, 2] = double.Parse(parts[2]); // MeanElevation_ft
    flow[i] = double.Parse(parts[3]);         // AnnualPeakFlow_cfs
}

// Fit regression model
var lm = new LinearRegression(new Matrix(features), new Vector(flow));

// Print summary
Console.WriteLine("Regional Streamflow Regression");
Console.WriteLine("=" + new string('=', 50));
foreach (var line in lm.Summary())
{
    Console.WriteLine(line);
}

// Predict for a new watershed
double[,] newSite = { { 100.0, 41.0, 780 } };
double[] predicted = lm.Predict(new Matrix(newSite));
double[,] interval = lm.PredictionIntervals(new Matrix(newSite), alpha: 0.1);

Console.WriteLine($"\nPrediction for 100 sq mi watershed:");
Console.WriteLine($"  Predicted peak flow: {predicted[0]:F0} cfs");
Console.WriteLine($"  90% Prediction interval: [{interval[0, 0]:F0}, {interval[0, 1]:F0}] cfs");
```

### Example 2: Residual Analysis

```cs
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;

double[,] X = {
    { 12.5 }, { 25.0 }, { 48.3 }, { 75.0 }, { 102.0 },
    { 18.7 }, { 55.0 }, { 130.0 }, { 8.2 }, { 200.0 }
};

double[] y = { 1450, 2680, 4520, 6100, 7850, 2050, 5200, 9300, 980, 12500 };

var lm = new LinearRegression(new Matrix(X), new Vector(y));

// Residual diagnostics
Console.WriteLine("Residual Diagnostics:");
Console.WriteLine($"  Mean residual: {lm.Residuals.Average():F2}");
Console.WriteLine($"  Std dev of residuals: {Statistics.StandardDeviation(lm.Residuals):F2}");
Console.WriteLine($"  Min residual: {lm.Residuals.Min():F2}");
Console.WriteLine($"  Max residual: {lm.Residuals.Max():F2}");

// Check normality of residuals
double jbPValue = HypothesisTests.JarqueBeraTest(lm.Residuals);
Console.WriteLine($"\nJarque-Bera test p-value: {jbPValue:F4}");

if (jbPValue > 0.05)
    Console.WriteLine("Residuals are consistent with normality (p > 0.05)");
else
    Console.WriteLine("Residuals may not be normally distributed (p < 0.05)");
```

## Best Practices

1. **Check assumptions** — Plot residuals vs. fitted values to detect nonlinearity and heteroscedasticity. Test residual normality with Jarque-Bera
2. **Examine residuals** — Residuals should appear random with no systematic pattern. Look for outliers (residuals beyond $\pm 3\sigma$)
3. **Watch for multicollinearity** — Highly correlated predictors inflate standard errors and destabilize coefficients. Compute VIF if suspected
4. **Use adjusted $R^2$** — Better for comparing models with different numbers of predictors, as it penalizes model complexity
5. **Validate predictions** — Don't extrapolate far beyond the range of training data. Prediction intervals widen rapidly outside the data range
6. **Consider transformations** — Log-transform skewed variables before fitting. For nonlinear relationships, consider the GLM framework available in the [Machine Learning](../machine-learning.md) documentation

---

## References

<a id="1">[1]</a> W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery, *Numerical Recipes: The Art of Scientific Computing*, 3rd ed., Cambridge, UK: Cambridge University Press, 2007.

<a id="2">[2]</a> G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed., Baltimore: Johns Hopkins University Press, 2013.

<a id="3">[3]</a> R. L. Burden and J. D. Faires, *Numerical Analysis*, 9th ed., Boston: Brooks/Cole, 2010.

---

[← Previous: Interpolation](interpolation.md) | [Back to Index](../index.md) | [Next: Time Series →](time-series.md)
