# Descriptive Statistics

[← Previous: Time Series](../data/time-series.md) | [Back to Index](../index.md) | [Next: Goodness-of-Fit →](goodness-of-fit.md)

The ***Numerics*** library provides comprehensive functions for computing descriptive statistics from data samples. The `Statistics` class contains static methods for all common statistical measures, sample moments, percentiles, and specialized analyses.

## Mathematical Foundations

This section defines the mathematical formulas underlying the descriptive statistics computed by the library. All formulas have been verified against the source implementation. The library uses numerically stable algorithms internally, but the results are mathematically equivalent to the definitions below.

### Sample Mean

The arithmetic mean of $n$ observations:

```math
\bar{x} = \frac{1}{n}\sum_{i=1}^{n} x_i
```

### Sample Variance (Bessel's Correction)

The unbiased estimator of the population variance divides by $n-1$ rather than $n$. This is known as Bessel's correction, and it compensates for the bias that arises when using the sample mean in place of the true population mean. The `Variance` method returns this estimator:

```math
s^2 = \frac{1}{n-1}\sum_{i=1}^{n}(x_i - \bar{x})^2
```

The source code implements this using a numerically stable one-pass recurrence rather than the naive two-pass formula shown above. The two approaches are mathematically equivalent, but the one-pass algorithm avoids catastrophic cancellation when the mean is large relative to the variance.

### Population Variance

When the data represents the entire population (not a sample), the `PopulationVariance` method divides by $n$:

```math
\sigma^2 = \frac{1}{n}\sum_{i=1}^{n}(x_i - \bar{x})^2
```

### Skewness (Fisher's Adjusted)

The `Skewness` method computes the bias-corrected Fisher skewness (type 2), which adjusts the raw sample skewness for sample size:

```math
G_1 = \frac{\sqrt{n(n-1)}}{n-2} \cdot g_1
```

where $g_1 = m_3 / m_2^{3/2}$ is the sample skewness computed from the central moments $m_k = \frac{1}{n}\sum_{i=1}^{n}(x_i - \bar{x})^k$.

**Interpretation:**
- $G_1 > 0$: right-skewed (long right tail)
- $G_1 < 0$: left-skewed (long left tail)
- $|G_1| < 0.5$: approximately symmetric

### Excess Kurtosis (Fisher's Adjusted)

The `Kurtosis` method computes bias-corrected excess kurtosis (type 2). Excess kurtosis subtracts 3 from Pearson's kurtosis so that the normal distribution has a value of zero:

```math
G_2 = \frac{n(n+1)}{(n-1)(n-2)(n-3)} \cdot \frac{\sum(x_i - \bar{x})^4}{s^4} - \frac{3(n-1)^2}{(n-2)(n-3)}
```

where $s^2 = \frac{1}{n-1}\sum(x_i - \bar{x})^2$ is the sample variance.

**Interpretation:**
- $G_2 > 0$: leptokurtic (heavier tails than normal)
- $G_2 < 0$: platykurtic (lighter tails than normal)
- $G_2 \approx 0$: mesokurtic (similar to normal distribution)

### Coefficient of Variation

The coefficient of variation expresses the standard deviation as a fraction of the mean, providing a dimensionless measure of relative variability:

```math
\text{CV} = \frac{s}{\bar{x}}
```

## Basic Statistics

### Central Tendency

Measures of the center or typical value of a dataset:

```cs
using Numerics.Data.Statistics;

double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Arithmetic mean
double mean = Statistics.Mean(data);  // 13.3

// Geometric mean (for positive data)
double geomMean = Statistics.GeometricMean(data);  // 13.13

// Harmonic mean
double harmMean = Statistics.HarmonicMean(data);  // 12.96

Console.WriteLine($"Arithmetic mean: {mean:F2}");
Console.WriteLine($"Geometric mean: {geomMean:F2}");
Console.WriteLine($"Harmonic mean: {harmMean:F2}");
```

**Note on Means:**
- **Arithmetic mean**: Best for symmetric data, $\bar{x} = \frac{1}{n}\sum x_i$
- **Geometric mean**: For multiplicative data (growth rates, ratios), $\sqrt[n]{\prod x_i}$
- **Harmonic mean**: For rates and ratios, $\frac{n}{\sum \frac{1}{x_i}}$
- Relationship: $\text{Harmonic} \leq \text{Geometric} \leq \text{Arithmetic}$

### Dispersion

Measures of spread or variability:

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Sample variance (divides by n-1)
double variance = Statistics.Variance(data);  // 4.08

// Population variance (divides by n)
double popVariance = Statistics.PopulationVariance(data);  // 3.67

// Sample standard deviation
double stdDev = Statistics.StandardDeviation(data);  // 2.02

// Population standard deviation
double popStdDev = Statistics.PopulationStandardDeviation(data);  // 1.92

// Coefficient of variation (relative variability)
double cv = Statistics.CoefficientOfVariation(data);  // 0.152 (15.2%)

Console.WriteLine($"Variance: {variance:F2}");
Console.WriteLine($"Std Dev: {stdDev:F2}");
Console.WriteLine($"CV: {cv:P1}");
```

### Efficient Combined Computations

For better performance when needing multiple related statistics:

```cs
// Compute mean and variance together (single pass through data)
var (mean, variance) = Statistics.MeanVariance(data);

// Compute mean and standard deviation together
var (mean2, stdDev) = Statistics.MeanStandardDeviation(data);

Console.WriteLine($"Mean: {mean:F2}, Variance: {variance:F2}");
Console.WriteLine($"Mean: {mean2:F2}, Std Dev: {stdDev:F2}");
```

### Range

```cs
// Minimum and maximum
double min = Statistics.Minimum(data);  // 10.5
double max = Statistics.Maximum(data);  // 16.8
double range = max - min;               // 6.3

// Sum
double sum = Statistics.Sum(data);      // 133.0

Console.WriteLine($"Range: [{min:F1}, {max:F1}] with width {range:F1}");
```

## Shape Statistics

### Skewness

Measures asymmetry of the distribution:

```cs
double[] symmetric = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
double[] rightSkewed = { 1, 1, 2, 2, 3, 5, 8, 13, 21 };
double[] leftSkewed = { 21, 13, 8, 5, 3, 2, 2, 1, 1 };

double skew1 = Statistics.Skewness(symmetric);     // ≈ 0 (symmetric)
double skew2 = Statistics.Skewness(rightSkewed);   // > 0 (right-skewed)
double skew3 = Statistics.Skewness(leftSkewed);    // < 0 (left-skewed)

Console.WriteLine($"Symmetric data skewness: {skew1:F3}");
Console.WriteLine($"Right-skewed data: {skew2:F3}");
Console.WriteLine($"Left-skewed data: {skew3:F3}");

// Interpretation
if (Math.Abs(skew1) < 0.5)
    Console.WriteLine("Approximately symmetric");
else if (skew1 > 0)
    Console.WriteLine("Right-skewed (tail extends right)");
else
    Console.WriteLine("Left-skewed (tail extends left)");
```

### Kurtosis

Measures tail heaviness (peakedness):

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Excess kurtosis (subtract 3 from Pearson's kurtosis)
double kurtosis = Statistics.Kurtosis(data);

Console.WriteLine($"Kurtosis: {kurtosis:F3}");

// Interpretation
if (kurtosis > 0)
    Console.WriteLine("Leptokurtic - heavier tails than normal (excess kurtosis > 0)");
else if (kurtosis < 0)
    Console.WriteLine("Platykurtic - lighter tails than normal (excess kurtosis < 0)");
else
    Console.WriteLine("Mesokurtic - similar to normal distribution (excess kurtosis ≈ 0)");
```

## Moments

### Product Moments

The first four product moments (mean, standard deviation, skewness, kurtosis):

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

double[] moments = Statistics.ProductMoments(data);

Console.WriteLine("Product Moments:");
Console.WriteLine($"  Mean (μ): {moments[0]:F2}");
Console.WriteLine($"  Standard Deviation (σ): {moments[1]:F2}");
Console.WriteLine($"  Skewness (γ₁): {moments[2]:F3}");
Console.WriteLine($"  Kurtosis (γ₂): {moments[3]:F3}");
```

### Linear Moments (L-Moments)

Linear moments (L-moments) are robust alternatives to conventional product moments, introduced by Hosking [[1]](#1). They are defined through probability weighted moments (PWMs), which use order statistics rather than powers of deviations from the mean.

The PWMs are computed as:

```math
\beta_r = \frac{1}{n}\sum_{i=1}^{n}\frac{\binom{i-1}{r}}{\binom{n-1}{r}} x_{i:n}
```

where $x_{i:n}$ are the order statistics (sorted data values). The library computes PWMs $\beta_0$ through $\beta_3$, from which the first four L-moments are derived:

```math
\lambda_1 = \beta_0, \quad \lambda_2 = 2\beta_1 - \beta_0
```

The `LinearMoments` method returns $\lambda_1$ (L-location), $\lambda_2$ (L-scale), and the L-moment ratios $\tau_3$ (L-skewness) and $\tau_4$ (L-kurtosis):

```math
\tau_3 = \frac{\lambda_3}{\lambda_2}, \quad \tau_4 = \frac{\lambda_4}{\lambda_2}
```

where $\lambda_3$ and $\lambda_4$ are the third and fourth L-moments, computed from PWMs $\beta_0$ through $\beta_3$.

**Why use L-moments over product moments?**

- **Robust to outliers**: L-moments are based on order statistics, not powers of deviations, making them far less sensitive to extreme values.
- **Better estimators for small samples**: When $n < 50$, L-moment estimators have lower bias and variance than product moment estimators.
- **Unique characterization**: L-moments uniquely characterize a distribution, similar to moment generating functions.
- **Bounded L-moment ratios**: $-1 \leq \tau_3 \leq 1$ and $\tau_4 \geq (5\tau_3^2 - 1)/4$, which makes them easy to interpret and compare.
- **Standard in hydrology**: L-moments are the preferred approach for flood frequency analysis per Hosking (1990) [[1]](#1).

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

double[] lMoments = Statistics.LinearMoments(data);

Console.WriteLine("L-Moments:");
Console.WriteLine($"  λ₁ (L-location): {lMoments[0]:F2}");
Console.WriteLine($"  λ₂ (L-scale): {lMoments[1]:F2}");
Console.WriteLine($"  τ₃ (L-skewness): {lMoments[2]:F4}");
Console.WriteLine($"  τ₄ (L-kurtosis): {lMoments[3]:F4}");

// L-moments are preferred for:
// - Small samples (n < 50)
// - Data with outliers
// - Hydrological data
// - Extreme value analysis
```

## Percentiles and Quantiles

The library computes percentiles using Type 7 linear interpolation, which is the default method in R and Excel. For a given probability $p \in [0, 1]$:

```math
Q(p) = x_{\lfloor h \rfloor} + (h - \lfloor h \rfloor)(x_{\lceil h \rceil} - x_{\lfloor h \rfloor})
```

where $h = (n-1)p$ and $x_i$ are the sorted data values (0-indexed). When $h$ falls exactly on an integer, the result is the corresponding order statistic; otherwise, it linearly interpolates between the two adjacent order statistics.

### Computing Percentiles

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Single percentile (k in [0, 1])
double median = Statistics.Percentile(data, 0.50);  // 50th percentile
double p90 = Statistics.Percentile(data, 0.90);     // 90th percentile
double p95 = Statistics.Percentile(data, 0.95);     // 95th percentile

Console.WriteLine($"Median (50th percentile): {median:F2}");
Console.WriteLine($"90th percentile: {p90:F2}");
Console.WriteLine($"95th percentile: {p95:F2}");

// Multiple percentiles at once (more efficient)
double[] percentiles = Statistics.Percentile(data, new double[] { 0.25, 0.50, 0.75, 0.90, 0.95 });

Console.WriteLine("\nPercentiles:");
Console.WriteLine($"  25th: {percentiles[0]:F2}");
Console.WriteLine($"  50th: {percentiles[1]:F2}");
Console.WriteLine($"  75th: {percentiles[2]:F2}");
Console.WriteLine($"  90th: {percentiles[3]:F2}");
Console.WriteLine($"  95th: {percentiles[4]:F2}");

// Note: Can specify if data is already sorted for efficiency
bool isSorted = false;
double q25 = Statistics.Percentile(data, 0.25, dataIsSorted: isSorted);
```

### Five-Number Summary

Box plot statistics (Tukey's five-number summary):

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

double[] fiveNum = Statistics.FiveNumberSummary(data);

Console.WriteLine("Five-Number Summary:");
Console.WriteLine($"  Minimum: {fiveNum[0]:F2}");
Console.WriteLine($"  Q1 (25th percentile): {fiveNum[1]:F2}");
Console.WriteLine($"  Median (50th percentile): {fiveNum[2]:F2}");
Console.WriteLine($"  Q3 (75th percentile): {fiveNum[3]:F2}");
Console.WriteLine($"  Maximum: {fiveNum[4]:F2}");

// Interquartile range
double iqr = fiveNum[3] - fiveNum[1];
Console.WriteLine($"  IQR: {iqr:F2}");

// Outlier bounds (Tukey's fences)
double lowerFence = fiveNum[1] - 1.5 * iqr;
double upperFence = fiveNum[3] + 1.5 * iqr;
Console.WriteLine($"  Outlier bounds: [{lowerFence:F2}, {upperFence:F2}]");
```

### Seven-Number Summary

Extended summary including additional percentiles:

```cs
double[] sevenNum = Statistics.SevenNumberSummary(data);

Console.WriteLine("Seven-Number Summary:");
Console.WriteLine($"  Minimum: {sevenNum[0]:F2}");
Console.WriteLine($"  5th percentile: {sevenNum[1]:F2}");
Console.WriteLine($"  Q1 (25th): {sevenNum[2]:F2}");
Console.WriteLine($"  Median (50th): {sevenNum[3]:F2}");
Console.WriteLine($"  Q3 (75th): {sevenNum[4]:F2}");
Console.WriteLine($"  95th percentile: {sevenNum[5]:F2}");
Console.WriteLine($"  Maximum: {sevenNum[6]:F2}");
```

## Covariance and Correlation

### Covariance

Measures how two variables vary together:

```cs
double[] x = { 1, 2, 3, 4, 5 };
double[] y = { 2, 4, 5, 4, 5 };

// Sample covariance
double cov = Statistics.Covariance(x, y);

// Population covariance
double popCov = Statistics.PopulationCovariance(x, y);

Console.WriteLine($"Sample covariance: {cov:F3}");
Console.WriteLine($"Population covariance: {popCov:F3}");

// Interpretation
if (cov > 0)
    Console.WriteLine("Positive association: as x increases, y tends to increase");
else if (cov < 0)
    Console.WriteLine("Negative association: as x increases, y tends to decrease");
else
    Console.WriteLine("No linear association");
```

### Correlation

For correlation coefficients, use the `Correlation` class. Three measures are provided, each capturing a different aspect of association between two variables.

**Pearson correlation** measures the strength and direction of the linear relationship between two variables [[3]](#3):

```math
r = \frac{\sum_{i=1}^{n}(x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum_{i=1}^{n}(x_i - \bar{x})^2 \cdot \sum_{i=1}^{n}(y_i - \bar{y})^2}}
```

**Spearman rank correlation** is the Pearson correlation applied to the ranks of the data rather than the data values themselves. This makes it a nonparametric measure of monotonic association:

```math
\rho_s = r(\text{rank}(x), \text{rank}(y))
```

**Kendall's tau-b** [[4]](#4) measures the strength of ordinal association, adjusted for ties. The library implements tau-b:

```math
\tau_b = \frac{n_c - n_d}{\sqrt{n_1 \cdot n_2}}
```

where $n_c$ is the number of concordant pairs, $n_d$ is the number of discordant pairs, $n_1$ is the number of pairs not tied in $x$, and $n_2$ is the number of pairs not tied in $y$.

**When to use which correlation:**

| Correlation | Measures | Robust to Outliers | Handles Nonlinear |
|------------|----------|-------------------|------------------|
| Pearson | Linear association | No | No |
| Spearman | Monotonic association | Yes | Yes (monotonic) |
| Kendall's $\tau_b$ | Concordance | Yes | Yes (monotonic) |

Rules of thumb:
- **Pearson**: use for normally distributed data with linear relationships.
- **Spearman**: use when data has outliers or the relationship is monotonic but not linear.
- **Kendall's $\tau_b$**: more robust than Spearman for small samples and has better statistical properties. Preferred when there are many tied values.

The `Correlation` class also provides matrix versions (`Pearson(double[,])` and `Spearman(double[,])`) that compute the full $p \times p$ correlation matrix for multivariate data.

```cs
using Numerics.Data.Statistics;

double[] x = { 1, 2, 3, 4, 5 };
double[] y = { 2, 4, 5, 4, 5 };

// Pearson correlation coefficient
double pearson = Correlation.Pearson(x, y);

// Spearman rank correlation
double spearman = Correlation.Spearman(x, y);

// Kendall tau correlation
double kendall = Correlation.KendallsTau(x, y);

Console.WriteLine($"Pearson r: {pearson:F3}");
Console.WriteLine($"Spearman ρ: {spearman:F3}");
Console.WriteLine($"Kendall τ: {kendall:F3}");

// Interpretation of Pearson r
if (Math.Abs(pearson) > 0.7)
    Console.WriteLine("Strong correlation");
else if (Math.Abs(pearson) > 0.3)
    Console.WriteLine("Moderate correlation");
else
    Console.WriteLine("Weak correlation");
```

## Ranking

### Rank Statistics

```cs
using System.Linq;

double[] data = { 5.2, 3.1, 7.8, 3.1, 9.2, 5.2 };

// Compute ranks (in-place, modifies array)
double[] dataCopy = (double[])data.Clone();
double[] ranks = Statistics.RanksInPlace(dataCopy);

Console.WriteLine("Value | Rank");
for (int i = 0; i < data.Length; i++)
{
    Console.WriteLine($"{data[i],5:F1} | {ranks[i],4:F1}");
}

// Ranks with ties reported
double[] dataCopy2 = (double[])data.Clone();
double[] ranks2 = Statistics.RanksInPlace(dataCopy2, out double[] ties);

Console.WriteLine($"\nNumber of tied groups: {ties.Count(t => t > 1)}");
```

## Entropy

Shannon entropy for continuous distributions:

```cs
using Numerics.Distributions;

double[] sample = new Normal(0, 1).GenerateRandomValues(1000);

// Estimate entropy using kernel density
var kde = new KernelDensity(sample, KernelDensity.KernelType.Gaussian, 0.5);
Func<double, double> pdf = x => kde.PDF(x);

double entropy = Statistics.Entropy(sample, pdf);

Console.WriteLine($"Estimated entropy: {entropy:F3} nats");
Console.WriteLine($"In bits: {entropy / Math.Log(2):F3}");

// For Normal(0,1), theoretical entropy ≈ 1.42 nats
```

## Jackknife Resampling

The jackknife is a resampling technique that estimates the variability of a statistic $\hat{\theta}$ by systematically leaving out one observation at a time. It is particularly useful for estimating standard errors of statistics that lack closed-form variance expressions (e.g., the median, L-moment ratios, or custom estimators).

### Jackknife Standard Error

The `JackKnifeStandardError` method computes:

```math
\text{SE}_{\text{jack}} = \sqrt{\frac{n-1}{n}\sum_{i=1}^{n}(\hat{\theta}_{(-i)} - \hat{\theta})^2}
```

where $\hat{\theta}_{(-i)}$ is the statistic computed without the $i$-th observation, and $\hat{\theta}$ is the statistic computed on the full sample.

### Jackknife Bias Estimate

The `JackKnifeSample` method returns the array of leave-one-out estimates $\hat{\theta}_{(-i)}$, from which the jackknife bias estimate can be computed:

```math
\text{bias}_{\text{jack}} = (n-1)(\bar{\theta}_{(\cdot)} - \hat{\theta})
```

where $\bar{\theta}_{(\cdot)} = \frac{1}{n}\sum_{i=1}^{n}\hat{\theta}_{(-i)}$ is the mean of the jackknife replications.

Note that the library parallelizes the jackknife loop internally using `Parallel.For`, so it scales well with multi-core processors.

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Define a statistic function (e.g., median)
Func<IList<double>, double> medianFunc = sample => Statistics.Percentile(sample.ToArray(), 0.50);

// Jackknife standard error
double jackknifeSE = Statistics.JackKnifeStandardError(data, medianFunc);

Console.WriteLine($"Median: {medianFunc(data):F2}");
Console.WriteLine($"Jackknife SE: {jackknifeSE:F3}");

// Get all jackknife samples (returns null if data is empty)
double[]? jackknifeValues = Statistics.JackKnifeSample(data, medianFunc);

if (jackknifeValues != null)
{
    Console.WriteLine($"Jackknife samples: {jackknifeValues.Length}");
    Console.WriteLine($"Mean of jackknife estimates: {jackknifeValues.Average():F2}");
}
```

## Practical Examples

### Example 1: Complete Streamflow Data Summary

This example analyzes the Tippecanoe River annual peak streamflow data. Statistics are validated against R's `base`, `psych`, `EnvStats`, and `lmom` packages.

**Data source:** Rao, A. R. & Hamed, K. H. (2000). *Flood Frequency Analysis*. CRC Press, Table 5.1.1.
See also: [`example-data/tippecanoe-river-streamflow.csv`](../example-data/tippecanoe-river-streamflow.csv)

```cs
using Numerics.Data.Statistics;

// Tippecanoe River near Delphi, IN — 48 years of annual peak streamflow (cfs)
double[] annualPeaks = {
    6290, 2700, 13100, 16900, 14600, 9600, 7740, 8490, 8130, 12000,
    17200, 15000, 12400, 6960, 6500, 5840, 10400, 18800, 21400, 22600,
    14200, 11000, 12800, 15700, 4740, 6950, 11800, 12100, 20600, 14600,
    14600, 8900, 10600, 14200, 14100, 14100, 12500, 7530, 13400, 17600,
    13400, 19200, 16900, 15500, 14500, 21900, 10400, 7460
};

Console.WriteLine("Tippecanoe River Annual Peak Streamflow Analysis (cfs)");
Console.WriteLine("=" + new string('=', 55));

// Central tendency
Console.WriteLine("\nCentral Tendency:");
Console.WriteLine($"  Mean:           {Statistics.Mean(annualPeaks):F1} cfs");
Console.WriteLine($"  Median:         {Statistics.Percentile(annualPeaks, 0.50):F1} cfs");

// Dispersion
Console.WriteLine("\nDispersion:");
Console.WriteLine($"  Range:          {Statistics.Minimum(annualPeaks):F0} - {Statistics.Maximum(annualPeaks):F0} cfs");
Console.WriteLine($"  Std Dev:        {Statistics.StandardDeviation(annualPeaks):F1} cfs");
Console.WriteLine($"  CV:             {Statistics.CoefficientOfVariation(annualPeaks):P1}");

// Shape
Console.WriteLine("\nShape:");
Console.WriteLine($"  Skewness:       {Statistics.Skewness(annualPeaks):F4}");
Console.WriteLine($"  Kurtosis:       {Statistics.Kurtosis(annualPeaks):F4}");

// L-Moments (more robust for hydrological data)
double[] lMoments = Statistics.LinearMoments(annualPeaks);
Console.WriteLine("\nL-Moments:");
Console.WriteLine($"  λ₁ (L-mean):    {lMoments[0]:F1}");
Console.WriteLine($"  λ₂ (L-scale):   {lMoments[1]:F1}");
Console.WriteLine($"  τ₃ (L-skew):    {lMoments[2]:F4}");
Console.WriteLine($"  τ₄ (L-kurtosis):{lMoments[3]:F4}");

// Percentiles (five-number summary)
Console.WriteLine("\nFive-Number Summary:");
Console.WriteLine($"  Min:  {Statistics.Minimum(annualPeaks):F0} cfs");
Console.WriteLine($"  Q1:   {Statistics.Percentile(annualPeaks, 0.25):F0} cfs");
Console.WriteLine($"  Q2:   {Statistics.Percentile(annualPeaks, 0.50):F0} cfs");
Console.WriteLine($"  Q3:   {Statistics.Percentile(annualPeaks, 0.75):F0} cfs");
Console.WriteLine($"  Max:  {Statistics.Maximum(annualPeaks):F0} cfs");
```

### Example 2: Comparing Two Datasets

```cs
double[] before = { 85, 92, 78, 95, 88, 91, 82, 89 };
double[] after = { 88, 95, 81, 98, 91, 94, 85, 92 };

Console.WriteLine("Before vs After Treatment");
Console.WriteLine("=" + new string('=', 50));

Console.WriteLine($"\n{"Statistic",-20} | {"Before",10} | {"After",10} | {"Change",10}");
Console.WriteLine(new string('-', 55));

double meanBefore = Statistics.Mean(before);
double meanAfter = Statistics.Mean(after);
Console.WriteLine($"{"Mean",-20} | {meanBefore,10:F2} | {meanAfter,10:F2} | {meanAfter - meanBefore,10:F2}");

double sdBefore = Statistics.StandardDeviation(before);
double sdAfter = Statistics.StandardDeviation(after);
Console.WriteLine($"{"Std Dev",-20} | {sdBefore,10:F2} | {sdAfter,10:F2} | {sdAfter - sdBefore,10:F2}");

double medBefore = Statistics.Percentile(before, 0.50);
double medAfter = Statistics.Percentile(after, 0.50);
Console.WriteLine($"{"Median",-20} | {medBefore,10:F2} | {medAfter,10:F2} | {medAfter - medBefore,10:F2}");

// Effect size (Cohen's d)
double pooledSD = Math.Sqrt((Statistics.Variance(before) + Statistics.Variance(after)) / 2);
double cohenD = (meanAfter - meanBefore) / pooledSD;
Console.WriteLine($"\nEffect size (Cohen's d): {cohenD:F3}");
```

### Example 3: Time Series Summary

```cs
double[] monthlyFlow = { 125, 135, 180, 220, 250, 280, 260, 230, 190, 150, 130, 120 };
string[] months = { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" };

Console.WriteLine("Monthly Streamflow Summary (cfs)");
Console.WriteLine("=" + new string('=', 50));

// Overall statistics
Console.WriteLine($"\nAnnual Statistics:");
Console.WriteLine($"  Mean: {Statistics.Mean(monthlyFlow):F0} cfs");
Console.WriteLine($"  Std Dev: {Statistics.StandardDeviation(monthlyFlow):F0} cfs");
Console.WriteLine($"  Min: {Statistics.Minimum(monthlyFlow):F0} cfs ({months[Array.IndexOf(monthlyFlow, Statistics.Minimum(monthlyFlow))]})");
Console.WriteLine($"  Max: {Statistics.Maximum(monthlyFlow):F0} cfs ({months[Array.IndexOf(monthlyFlow, Statistics.Maximum(monthlyFlow))]})");

// Seasonal means
double springMean = Statistics.Mean(new[] { monthlyFlow[2], monthlyFlow[3], monthlyFlow[4] });
double summerMean = Statistics.Mean(new[] { monthlyFlow[5], monthlyFlow[6], monthlyFlow[7] });
double fallMean = Statistics.Mean(new[] { monthlyFlow[8], monthlyFlow[9], monthlyFlow[10] });
double winterMean = Statistics.Mean(new[] { monthlyFlow[11], monthlyFlow[0], monthlyFlow[1] });

Console.WriteLine($"\nSeasonal Means:");
Console.WriteLine($"  Spring (MAM): {springMean:F0} cfs");
Console.WriteLine($"  Summer (JJA): {summerMean:F0} cfs");
Console.WriteLine($"  Fall (SON): {fallMean:F0} cfs");
Console.WriteLine($"  Winter (DJF): {winterMean:F0} cfs");
```

## Running Statistics

The `RunningStatistics` class implements Welford's online algorithm [[2]](#2) for numerically stable computation of variance and higher moments in a single pass through the data. This is essential for streaming data or very large datasets where it is impractical to store all values in memory.

### Welford's Algorithm

The classical two-pass approach (compute mean first, then deviations) requires two full scans of the data and storing the entire dataset. Welford's algorithm updates the running mean and sum of squared deviations incrementally with each new observation $x_n$:

```math
\delta = x_n - M_1^{(n-1)}
```

```math
M_1^{(n)} = M_1^{(n-1)} + \frac{\delta}{n}
```

```math
M_2^{(n)} = M_2^{(n-1)} + \delta(x_n - M_1^{(n)})
```

where $M_1$ is the running mean and $M_2$ is the running sum of squared deviations. The sample variance is then $s^2 = M_2 / (n-1)$.

The `RunningStatistics` class extends this to the fourth moment ($M_3$ and $M_4$) for skewness and kurtosis computation. It also supports combining two `RunningStatistics` objects via the `Combine` method (or `+` operator), which is useful for parallel processing -- each thread can accumulate statistics independently and then merge results.

**Available properties:** `Count`, `Minimum`, `Maximum`, `Mean`, `Variance`, `PopulationVariance`, `StandardDeviation`, `PopulationStandardDeviation`, `CoefficientOfVariation`, `Skewness`, `PopulationSkewness`, `Kurtosis`, `PopulationKurtosis`.

```cs
using Numerics.Data.Statistics;

var runningStats = new RunningStatistics();

// Add data points incrementally
double[] newData = { 10.5, 12.3, 11.8, 15.2, 13.7 };

foreach (var value in newData)
{
    runningStats.Push(value);
}

Console.WriteLine($"Count: {runningStats.Count}");
Console.WriteLine($"Mean: {runningStats.Mean:F2}");
Console.WriteLine($"Variance: {runningStats.Variance:F2}");
Console.WriteLine($"Std Dev: {runningStats.StandardDeviation:F2}");
Console.WriteLine($"Skewness: {runningStats.Skewness:F3}");
Console.WriteLine($"Kurtosis: {runningStats.Kurtosis:F3}");
Console.WriteLine($"Min: {runningStats.Minimum:F2}");
Console.WriteLine($"Max: {runningStats.Maximum:F2}");
```

## Best Practices

1. **Check for NaN and Inf**: All methods return NaN if data contains NaN or is empty
2. **Sort when needed**: Some methods require sorted data - check documentation
3. **Use appropriate sample size**: Small samples (n < 30) need careful interpretation
4. **Consider outliers**: L-moments are more robust than product moments
5. **Understand bias**: Sample statistics (n-1 denominator) vs population (n denominator)
6. **Use running statistics**: For streaming data or memory constraints
7. **Vectorize operations**: Use batch methods for multiple percentiles

---

## References

<a id="1">[1]</a> Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society: Series B (Methodological)*, 52(1), 105-124.

<a id="2">[2]</a> Welford, B. P. (1962). Note on a method for calculating corrected sums of squares and products. *Technometrics*, 4(3), 419-420.

<a id="3">[3]</a> Fisher, R. A. (1930). The moments of the distribution for normal samples of measures of departure from normality. *Proceedings of the Royal Society of London. Series A*, 130(812), 16-28.

<a id="4">[4]</a> Kendall, M. G. (1938). A new measure of rank correlation. *Biometrika*, 30(1/2), 81-93.

---

[← Previous: Time Series](../data/time-series.md) | [Back to Index](../index.md) | [Next: Goodness-of-Fit →](goodness-of-fit.md)
