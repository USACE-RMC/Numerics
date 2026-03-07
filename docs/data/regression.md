# Linear Regression

[← Previous: Interpolation](interpolation.md) | [Back to Index](../index.md) | [Next: Time Series →](time-series.md)

The ***Numerics*** library provides a `LinearRegression` class for estimating the linear relationship between a scalar response variable and one or more predictor variables. The model Y = α + βX + ε is fitted using Singular Value Decomposition (SVD) for numerical stability.

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
| `DegreesOfFreedom` | `int` | Residual degrees of freedom (n − p) |
| `RSquared` | `double` | Coefficient of determination |
| `AdjRSquared` | `double` | Adjusted R² |
| `HasIntercept` | `bool` | Whether the model includes an intercept |

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
- R² and adjusted R²
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

## Regression Without Intercept

```cs
// Force model through the origin (no intercept)
var lm = new LinearRegression(new Matrix(X), new Vector(y), hasIntercept: false);

Console.WriteLine("No-intercept model:");
Console.WriteLine($"Slope: {lm.Parameters[0]:F4}");
```

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

1. **Check assumptions** — Linear regression assumes linearity, independence, normality of errors, and constant variance
2. **Examine residuals** — Use the `Residuals` property to check for patterns
3. **Watch for multicollinearity** — Highly correlated predictors can inflate standard errors
4. **Use adjusted R²** — Better for comparing models with different numbers of predictors
5. **Validate predictions** — Don't extrapolate far beyond the range of training data
6. **Consider transformations** — Log-transform skewed variables before fitting

---

[← Previous: Interpolation](interpolation.md) | [Back to Index](../index.md) | [Next: Time Series →](time-series.md)
