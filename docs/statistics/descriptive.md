# Descriptive Statistics

[← Back to Index](../index.md) | [Next: Goodness-of-Fit →](goodness-of-fit.md)

The ***Numerics*** library provides comprehensive functions for computing descriptive statistics from data samples. The `Statistics` class contains static methods for all common statistical measures, sample moments, percentiles, and specialized analyses.

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

Linear moments are robust alternatives to product moments [[1]](#1):

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

### Computing Percentiles

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Single percentile (k as decimal: 0-100)
double median = Statistics.Percentile(data, 50);  // 50th percentile
double p90 = Statistics.Percentile(data, 90);     // 90th percentile
double p95 = Statistics.Percentile(data, 95);     // 95th percentile

Console.WriteLine($"Median (50th percentile): {median:F2}");
Console.WriteLine($"90th percentile: {p90:F2}");
Console.WriteLine($"95th percentile: {p95:F2}");

// Multiple percentiles at once (more efficient)
double[] percentiles = Statistics.Percentile(data, new double[] { 25, 50, 75, 90, 95 });

Console.WriteLine("\nPercentiles:");
Console.WriteLine($"  25th: {percentiles[0]:F2}");
Console.WriteLine($"  50th: {percentiles[1]:F2}");
Console.WriteLine($"  75th: {percentiles[2]:F2}");
Console.WriteLine($"  90th: {percentiles[3]:F2}");
Console.WriteLine($"  95th: {percentiles[4]:F2}");

// Note: Can specify if data is already sorted for efficiency
bool isSorted = false;
double q25 = Statistics.Percentile(data, 25, dataIsSorted: isSorted);
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
Console.WriteLine($"  10th percentile: {sevenNum[1]:F2}");
Console.WriteLine($"  Q1 (25th): {sevenNum[2]:F2}");
Console.WriteLine($"  Median (50th): {sevenNum[3]:F2}");
Console.WriteLine($"  Q3 (75th): {sevenNum[4]:F2}");
Console.WriteLine($"  90th percentile: {sevenNum[5]:F2}");
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

For correlation coefficients, use the `Correlation` class:

```cs
using Numerics.Data.Statistics;

double[] x = { 1, 2, 3, 4, 5 };
double[] y = { 2, 4, 5, 4, 5 };

// Pearson correlation coefficient
double pearson = Correlation.Pearson(x, y);

// Spearman rank correlation
double spearman = Correlation.Spearman(x, y);

// Kendall tau correlation
double kendall = Correlation.Kendall(x, y);

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
Func<double, double> pdf = x =>
{
    var kde = new KernelDensity(sample, bandwidth: 0.5);
    return kde.PDF(x);
};

double entropy = Statistics.Entropy(sample, pdf);

Console.WriteLine($"Estimated entropy: {entropy:F3} nats");
Console.WriteLine($"In bits: {entropy / Math.Log(2):F3}");

// For Normal(0,1), theoretical entropy ≈ 1.42 nats
```

## Jackknife Resampling

Leave-one-out resampling for standard error estimation:

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Define a statistic function (e.g., median)
Func<IList<double>, double> medianFunc = sample => Statistics.Percentile(sample.ToArray(), 50);

// Jackknife standard error
double jackknifeSE = Statistics.JackKnifeStandardError(data, medianFunc);

Console.WriteLine($"Median: {medianFunc(data):F2}");
Console.WriteLine($"Jackknife SE: {jackknifeSE:F3}");

// Get all jackknife samples
double[] jackknifeValues = Statistics.JackKnifeSample(data, medianFunc);

Console.WriteLine($"Jackknife samples: {jackknifeValues.Length}");
Console.WriteLine($"Mean of jackknife estimates: {jackknifeValues.Average():F2}");
```

## Practical Examples

### Example 1: Complete Data Summary

```cs
using Numerics.Data.Statistics;

double[] annualRainfall = { 850, 920, 780, 1050, 890, 950, 820, 1100, 870, 980 };

Console.WriteLine("Annual Rainfall Analysis (mm)");
Console.WriteLine("=" + new string('=', 50));

// Central tendency
Console.WriteLine("\nCentral Tendency:");
Console.WriteLine($"  Mean: {Statistics.Mean(annualRainfall):F1} mm");
Console.WriteLine($"  Median: {Statistics.Percentile(annualRainfall, 50):F1} mm");

// Dispersion
Console.WriteLine("\nDispersion:");
Console.WriteLine($"  Range: {Statistics.Minimum(annualRainfall):F0} - {Statistics.Maximum(annualRainfall):F0} mm");
Console.WriteLine($"  Std Dev: {Statistics.StandardDeviation(annualRainfall):F1} mm");
Console.WriteLine($"  CV: {Statistics.CoefficientOfVariation(annualRainfall):P1}");

// Shape
Console.WriteLine("\nShape:");
Console.WriteLine($"  Skewness: {Statistics.Skewness(annualRainfall):F3}");
Console.WriteLine($"  Kurtosis: {Statistics.Kurtosis(annualRainfall):F3}");

// Percentiles
var percentiles = Statistics.Percentile(annualRainfall, new double[] { 10, 25, 50, 75, 90 });
Console.WriteLine("\nPercentiles:");
Console.WriteLine($"  10th: {percentiles[0]:F1} mm");
Console.WriteLine($"  25th: {percentiles[1]:F1} mm");
Console.WriteLine($"  50th: {percentiles[2]:F1} mm");
Console.WriteLine($"  75th: {percentiles[3]:F1} mm");
Console.WriteLine($"  90th: {percentiles[4]:F1} mm");
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

double medBefore = Statistics.Percentile(before, 50);
double medAfter = Statistics.Percentile(after, 50);
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

For streaming data or very large datasets:

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

---

[← Back to Index](../index.md) | [Next: Goodness-of-Fit →](goodness-of-fit.md)
