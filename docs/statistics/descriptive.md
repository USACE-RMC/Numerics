# Descriptive Statistics

The `Statistics` class in ***Numerics*** provides comprehensive functions for computing summary statistics, moments, and data transformations.

## Basic Statistics

### Central Tendency

```cs
using Numerics.Data.Statistics;

double[] data = { 10, 15, 12, 18, 14, 16, 13, 17, 11, 19, 20, 8 };

// Measures of central tendency
double mean = Statistics.Mean(data);                 // Arithmetic mean
double median = Statistics.Median(data);             // 50th percentile
double geometricMean = Statistics.GeometricMean(data); // For positive data
double harmonicMean = Statistics.HarmonicMean(data);   // For rates

Console.WriteLine($"Mean:           {mean:F2}");
Console.WriteLine($"Median:         {median:F2}");
Console.WriteLine($"Geometric Mean: {geometricMean:F2}");
Console.WriteLine($"Harmonic Mean:  {harmonicMean:F2}");
```

### Dispersion

```cs
double variance = Statistics.Variance(data);           // Sample variance
double stdDev = Statistics.StandardDeviation(data);    // Sample std dev
double popVar = Statistics.PopulationVariance(data);   // Population variance
double popStdDev = Statistics.PopulationStandardDeviation(data);

double range = Statistics.Maximum(data) - Statistics.Minimum(data);
double iqr = Statistics.InterquartileRange(data);      // Q3 - Q1

Console.WriteLine($"Variance:    {variance:F2}");
Console.WriteLine($"Std Dev:     {stdDev:F2}");
Console.WriteLine($"Range:       {range:F2}");
Console.WriteLine($"IQR:         {iqr:F2}");
```

### Extremes

```cs
double min = Statistics.Minimum(data);
double max = Statistics.Maximum(data);
double sum = Statistics.Sum(data);

Console.WriteLine($"Minimum: {min}");
Console.WriteLine($"Maximum: {max}");
Console.WriteLine($"Sum:     {sum}");
```

---

## Percentiles and Quantiles

### Computing Percentiles

```cs
// Single percentile
double p50 = Statistics.Percentile(data, 0.50);  // Median
double p90 = Statistics.Percentile(data, 0.90);  // 90th percentile
double p95 = Statistics.Percentile(data, 0.95);  // 95th percentile

// Multiple percentiles at once
double[] probs = { 0.10, 0.25, 0.50, 0.75, 0.90 };
double[] percentiles = Statistics.Percentiles(data, probs);

Console.WriteLine("Percentile Summary:");
for (int i = 0; i < probs.Length; i++)
{
    Console.WriteLine($"  {probs[i]*100:F0}%: {percentiles[i]:F2}");
}
```

### Quartiles

```cs
double q1 = Statistics.Percentile(data, 0.25);  // First quartile
double q2 = Statistics.Percentile(data, 0.50);  // Second quartile (median)
double q3 = Statistics.Percentile(data, 0.75);  // Third quartile
double iqr = q3 - q1;                           // Interquartile range

Console.WriteLine($"Q1: {q1:F2}");
Console.WriteLine($"Q2: {q2:F2}");
Console.WriteLine($"Q3: {q3:F2}");
Console.WriteLine($"IQR: {iqr:F2}");
```

### Ranks and Order Statistics

```cs
// Get ranks (1-based)
int[] ranks = Statistics.Ranks(data);

// Get sorted data with original indices
var sorted = data.Select((val, idx) => new { Value = val, Index = idx })
                 .OrderBy(x => x.Value)
                 .ToList();
```

---

## Moments

### Product Moments

Product moments characterize the shape of a distribution [[1]](#ref1).

```cs
// Central moments
double m2 = Statistics.CentralMoment(data, 2);  // Variance (biased)
double m3 = Statistics.CentralMoment(data, 3);  // Third central moment
double m4 = Statistics.CentralMoment(data, 4);  // Fourth central moment

// Standardized moments
double skewness = Statistics.Skewness(data);    // Third standardized moment
double kurtosis = Statistics.Kurtosis(data);    // Fourth standardized moment
double excessKurtosis = Statistics.ExcessKurtosis(data);  // Kurtosis - 3

Console.WriteLine($"Skewness: {skewness:F3}");
Console.WriteLine($"Kurtosis: {kurtosis:F3}");
Console.WriteLine($"Excess Kurtosis: {excessKurtosis:F3}");
```

### Interpreting Skewness

| Skewness | Interpretation |
|----------|----------------|
| γ < -1 | Highly left-skewed |
| -1 ≤ γ < -0.5 | Moderately left-skewed |
| -0.5 ≤ γ ≤ 0.5 | Approximately symmetric |
| 0.5 < γ ≤ 1 | Moderately right-skewed |
| γ > 1 | Highly right-skewed |

### Interpreting Kurtosis

| Excess Kurtosis | Interpretation |
|-----------------|----------------|
| κ < 0 | Platykurtic (lighter tails than Normal) |
| κ ≈ 0 | Mesokurtic (similar to Normal) |
| κ > 0 | Leptokurtic (heavier tails than Normal) |

---

## L-Moments

L-moments are linear combinations of order statistics, more robust than product moments [[2]](#ref2).

### Computing L-Moments

```cs
// First four L-moments
double[] lmom = Statistics.LinearMoments(data);

double L1 = lmom[0];  // L-location (mean)
double L2 = lmom[1];  // L-scale
double T3 = lmom[2];  // L-skewness (τ₃)
double T4 = lmom[3];  // L-kurtosis (τ₄)

Console.WriteLine($"L1 (location): {L1:F3}");
Console.WriteLine($"L2 (scale):    {L2:F3}");
Console.WriteLine($"τ3 (L-skew):   {T3:F4}");
Console.WriteLine($"τ4 (L-kurt):   {T4:F4}");
```

### L-Moment Definitions

```math
\lambda_1 = E[X]
```

```math
\lambda_2 = \frac{1}{2}E[X_{2:2} - X_{1:2}]
```

```math
\lambda_3 = \frac{1}{3}E[X_{3:3} - 2X_{2:3} + X_{1:3}]
```

```math
\lambda_4 = \frac{1}{4}E[X_{4:4} - 3X_{3:4} + 3X_{2:4} - X_{1:4}]
```

### L-Moment Ratios

| Ratio | Symbol | Range | Interpretation |
|-------|--------|-------|----------------|
| L-CV | τ = λ₂/λ₁ | [0, 1] | Coefficient of L-variation |
| L-Skewness | τ₃ = λ₃/λ₂ | [-1, 1] | Asymmetry |
| L-Kurtosis | τ₄ = λ₄/λ₂ | [-0.25, 1] | Tail heaviness |

### L-Moment Ratio Diagram

L-moment ratios help identify distribution families:

```cs
// Sample L-moments
double[] lmom = Statistics.LinearMoments(data);
double tau3 = lmom[2];
double tau4 = lmom[3];

// Compare to theoretical L-moment ratios
// Normal:     τ3 = 0,      τ4 = 0.1226
// Gumbel:     τ3 = 0.1699, τ4 = 0.1504
// Exponential: τ3 = 0.3333, τ4 = 0.1667

Console.WriteLine($"Sample: τ3 = {tau3:F4}, τ4 = {tau4:F4}");
Console.WriteLine($"Closest to: {IdentifyDistribution(tau3, tau4)}");
```

---

## Covariance and Correlation

### Covariance

```cs
double[] x = { 1, 2, 3, 4, 5 };
double[] y = { 2.1, 4.2, 5.8, 8.1, 9.9 };

// Sample covariance
double cov = Statistics.Covariance(x, y);

// Population covariance
double popCov = Statistics.PopulationCovariance(x, y);

Console.WriteLine($"Covariance: {cov:F4}");
```

### Pearson Correlation

Measures linear association [[3]](#ref3):

```math
r = \frac{\sum_{i=1}^{n}(x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum_{i=1}^{n}(x_i - \bar{x})^2}\sqrt{\sum_{i=1}^{n}(y_i - \bar{y})^2}}
```

```cs
double pearsonR = Statistics.Correlation(x, y);
Console.WriteLine($"Pearson r: {pearsonR:F4}");
```

### Spearman Rank Correlation

Non-parametric measure based on ranks:

```cs
double spearmanRho = Statistics.SpearmanCorrelation(x, y);
Console.WriteLine($"Spearman ρ: {spearmanRho:F4}");
```

### Kendall's Tau

Measures concordance:

```cs
double kendallTau = Statistics.KendallTau(x, y);
Console.WriteLine($"Kendall τ: {kendallTau:F4}");
```

### Correlation Interpretation

| |r| | Interpretation |
|-----|----------------|
| 0.0 - 0.2 | Very weak |
| 0.2 - 0.4 | Weak |
| 0.4 - 0.6 | Moderate |
| 0.6 - 0.8 | Strong |
| 0.8 - 1.0 | Very strong |

---

## Autocorrelation

Measures correlation of a time series with itself at different lags [[4]](#ref4).

### Sample Autocorrelation

```cs
double[] timeSeries = { /* ... */ };

// Autocorrelation at specific lag
int lag = 1;
double acf1 = Statistics.Autocorrelation(timeSeries, lag);

// Autocorrelation function for multiple lags
int maxLag = 20;
double[] acf = new double[maxLag + 1];
for (int k = 0; k <= maxLag; k++)
{
    acf[k] = Statistics.Autocorrelation(timeSeries, k);
}

Console.WriteLine("Lag   ACF");
for (int k = 0; k <= 10; k++)
{
    Console.WriteLine($"{k,3}   {acf[k]:F4}");
}
```

### Significance Bounds

For white noise, approximate 95% confidence bounds are:

```math
\pm \frac{1.96}{\sqrt{n}}
```

```cs
int n = timeSeries.Length;
double bound = 1.96 / Math.Sqrt(n);
Console.WriteLine($"95% bounds: ±{bound:F4}");
```

### Lag-1 Autocorrelation Adjustment

For frequency analysis, adjust for serial correlation:

```cs
double rho1 = Statistics.Autocorrelation(data, 1);
int n = data.Length;

// Effective sample size
double nEffective = n * (1 - rho1) / (1 + rho1);

Console.WriteLine($"Actual n: {n}");
Console.WriteLine($"Lag-1 autocorrelation: {rho1:F4}");
Console.WriteLine($"Effective n: {nEffective:F1}");
```

---

## Data Transformations

### Log Transformation

```cs
// Natural log transform
double[] logData = data.Select(x => Math.Log(x)).ToArray();

// Log10 transform
double[] log10Data = data.Select(x => Math.Log10(x)).ToArray();

// Statistics on transformed data
double logMean = Statistics.Mean(logData);
double logStd = Statistics.StandardDeviation(logData);
```

### Box-Cox Transformation

The Box-Cox transformation finds the optimal power transformation [[5]](#ref5):

```math
y^{(\lambda)} = \begin{cases}
\frac{y^\lambda - 1}{\lambda} & \lambda \neq 0 \\
\ln(y) & \lambda = 0
\end{cases}
```

```cs
// Apply Box-Cox transformation
double lambda = 0.5;  // Square root transformation
double[] transformed = data.Select(x => 
    lambda != 0 ? (Math.Pow(x, lambda) - 1) / lambda : Math.Log(x)
).ToArray();

// Find optimal lambda (minimize skewness)
double optimalLambda = FindOptimalBoxCox(data);
```

### Standardization (Z-scores)

```cs
double mean = Statistics.Mean(data);
double std = Statistics.StandardDeviation(data);

double[] zScores = data.Select(x => (x - mean) / std).ToArray();

// Z-scores have mean ≈ 0 and std ≈ 1
Console.WriteLine($"Z-score mean: {Statistics.Mean(zScores):F6}");
Console.WriteLine($"Z-score std:  {Statistics.StandardDeviation(zScores):F6}");
```

### Normalization (Min-Max Scaling)

```cs
double min = Statistics.Minimum(data);
double max = Statistics.Maximum(data);

// Scale to [0, 1]
double[] normalized = data.Select(x => (x - min) / (max - min)).ToArray();

// Scale to [a, b]
double a = 10, b = 100;
double[] scaled = data.Select(x => a + (b - a) * (x - min) / (max - min)).ToArray();
```

---

## Outlier Detection

### Multiple Grubbs-Beck Test

The Multiple Grubbs-Beck (MGB) test identifies low outliers in flood data [[6]](#ref6):

```cs
double[] annualMax = { 12500, 15200, 11800, 500, 18900, 14200, 16500, 800 };

// Detect low outliers
var mgbResult = Statistics.MultipleGrubbsBeck(annualMax, alpha: 0.10);

Console.WriteLine($"Low outlier threshold: {mgbResult.Threshold:F0}");
Console.WriteLine($"Number of outliers: {mgbResult.OutlierCount}");
Console.WriteLine($"Outliers: {string.Join(", ", mgbResult.Outliers)}");
```

### IQR Method

```cs
double q1 = Statistics.Percentile(data, 0.25);
double q3 = Statistics.Percentile(data, 0.75);
double iqr = q3 - q1;

double lowerFence = q1 - 1.5 * iqr;
double upperFence = q3 + 1.5 * iqr;

var outliers = data.Where(x => x < lowerFence || x > upperFence);

Console.WriteLine($"Lower fence: {lowerFence:F2}");
Console.WriteLine($"Upper fence: {upperFence:F2}");
Console.WriteLine($"Outliers: {string.Join(", ", outliers)}");
```

### Z-Score Method

```cs
double mean = Statistics.Mean(data);
double std = Statistics.StandardDeviation(data);
double threshold = 3.0;  // Standard threshold

var outliers = data.Where(x => Math.Abs((x - mean) / std) > threshold);
```

---

## Plotting Positions

Plotting positions assign probabilities to ordered observations for frequency analysis [[7]](#ref7).

### Common Formulas

| Method | Formula | α value | Best For |
|--------|---------|---------|----------|
| Weibull | $\frac{i}{n+1}$ | 0 | General use |
| Hazen | $\frac{i-0.5}{n}$ | 0.5 | Symmetric distributions |
| Cunnane | $\frac{i-0.4}{n+0.2}$ | 0.4 | GEV, Gumbel |
| Gringorten | $\frac{i-0.44}{n+0.12}$ | 0.44 | Gumbel |
| Blom | $\frac{i-0.375}{n+0.25}$ | 0.375 | Normal |

### General Formula

```math
p_i = \frac{i - \alpha}{n + 1 - 2\alpha}
```

```cs
double[] sortedData = data.OrderBy(x => x).ToArray();
int n = sortedData.Length;

// Cunnane plotting positions (α = 0.4)
double alpha = 0.4;
double[] plottingPositions = new double[n];

for (int i = 0; i < n; i++)
{
    plottingPositions[i] = (i + 1 - alpha) / (n + 1 - 2 * alpha);
}

Console.WriteLine("Value      Exceedance Prob    Return Period");
for (int i = 0; i < n; i++)
{
    double exceedProb = 1 - plottingPositions[i];
    double returnPeriod = 1 / exceedProb;
    Console.WriteLine($"{sortedData[i],10:F0}  {exceedProb,12:F4}    {returnPeriod,12:F1}");
}
```

---

## Summary Statistics Table

```cs
public static void PrintSummary(double[] data)
{
    Console.WriteLine("=== Summary Statistics ===");
    Console.WriteLine($"n:        {data.Length}");
    Console.WriteLine($"Mean:     {Statistics.Mean(data):F3}");
    Console.WriteLine($"Std Dev:  {Statistics.StandardDeviation(data):F3}");
    Console.WriteLine($"Variance: {Statistics.Variance(data):F3}");
    Console.WriteLine($"Minimum:  {Statistics.Minimum(data):F3}");
    Console.WriteLine($"Q1:       {Statistics.Percentile(data, 0.25):F3}");
    Console.WriteLine($"Median:   {Statistics.Median(data):F3}");
    Console.WriteLine($"Q3:       {Statistics.Percentile(data, 0.75):F3}");
    Console.WriteLine($"Maximum:  {Statistics.Maximum(data):F3}");
    Console.WriteLine($"IQR:      {Statistics.InterquartileRange(data):F3}");
    Console.WriteLine($"Skewness: {Statistics.Skewness(data):F4}");
    Console.WriteLine($"Kurtosis: {Statistics.Kurtosis(data):F4}");
    
    var lmom = Statistics.LinearMoments(data);
    Console.WriteLine($"L-CV:     {lmom[1]/lmom[0]:F4}");
    Console.WriteLine($"L-Skew:   {lmom[2]:F4}");
    Console.WriteLine($"L-Kurt:   {lmom[3]:F4}");
}
```

---

## References

<a id="ref1">[1]</a> Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society: Series B*, 52(1), 105-124.

<a id="ref2">[2]</a> Hosking, J. R. M., & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press.

<a id="ref3">[3]</a> Wilks, D. S. (2011). *Statistical Methods in the Atmospheric Sciences* (3rd ed.). Academic Press.

<a id="ref4">[4]</a> Box, G. E. P., Jenkins, G. M., Reinsel, G. C., & Ljung, G. M. (2015). *Time Series Analysis: Forecasting and Control* (5th ed.). Wiley.

<a id="ref5">[5]</a> Box, G. E. P., & Cox, D. R. (1964). An analysis of transformations. *Journal of the Royal Statistical Society: Series B*, 26(2), 211-252.

<a id="ref6">[6]</a> Cohn, T. A., England, J. F., Berenbrock, C. E., Mason, R. R., Stedinger, J. R., & Lamontagne, J. R. (2013). A generalized Grubbs-Beck test statistic for detecting multiple potentially influential low outliers in flood series. *Water Resources Research*, 49(8), 5047-5058.

<a id="ref7">[7]</a> Cunnane, C. (1978). Unbiased plotting positions—A review. *Journal of Hydrology*, 37(3-4), 205-222.
