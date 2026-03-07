# Hypothesis Tests

[← Previous: Goodness-of-Fit](goodness-of-fit.md) | [Back to Index](../index.md) | [Next: Convergence Diagnostics →](../sampling/convergence-diagnostics.md)

Statistical hypothesis testing is a method of statistical inference used to decide whether observed data sufficiently support a particular hypothesis. The ***Numerics*** library provides a comprehensive set of hypothesis tests that return **p-values** for making statistical decisions. If the p-value is less than the chosen significance level (typically α = 0.05), the null hypothesis is rejected in favor of the alternative hypothesis.

## Understanding p-values

All hypothesis tests in ***Numerics*** return **p-values**, not test statistics. The p-value represents the probability of obtaining results at least as extreme as the observed data, assuming the null hypothesis is true.

**Decision rule:**
- p < 0.01: Very strong evidence against H₀
- p < 0.05: Strong evidence against H₀
- p < 0.10: Weak evidence against H₀
- p ≥ 0.10: Insufficient evidence to reject H₀

## t-Tests

### One-Sample t-Test

Tests if the sample mean differs from a hypothesized population mean:

```cs
using Numerics.Data.Statistics;

double[] sample = { 12.5, 13.2, 11.8, 14.1, 12.9, 13.5, 12.2, 13.8 };
double mu0 = 12.0;  // Hypothesized mean

// Returns the two-sided p-value
double pValue = HypothesisTests.OneSampleTtest(sample, mu0);

Console.WriteLine($"One-sample t-test:");
Console.WriteLine($"  Sample mean: {Statistics.Mean(sample):F2}");
Console.WriteLine($"  Hypothesized mean: {mu0:F2}");
Console.WriteLine($"  p-value: {pValue:F4}");
Console.WriteLine($"  Degrees of freedom: {sample.Length - 1}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Reject H₀ - mean differs significantly from 12.0");
else
    Console.WriteLine("  Result: Fail to reject H₀ - insufficient evidence of difference");
```

**Hypotheses:**
- H₀: μ = μ₀
- H₁: μ ≠ μ₀ (two-tailed)

### Two-Sample t-Tests

#### Equal Variance (Pooled) t-Test

Tests if means of two independent samples are equal, assuming equal variances:

```cs
double[] sample1 = { 12.5, 13.2, 11.8, 14.1, 12.9 };
double[] sample2 = { 15.3, 14.8, 15.9, 14.5, 15.1 };

// Returns the two-sided p-value
double pValue = HypothesisTests.EqualVarianceTtest(sample1, sample2);

Console.WriteLine($"Equal variance t-test:");
Console.WriteLine($"  Sample 1 mean: {Statistics.Mean(sample1):F2}");
Console.WriteLine($"  Sample 2 mean: {Statistics.Mean(sample2):F2}");
Console.WriteLine($"  p-value: {pValue:F4}");
Console.WriteLine($"  Degrees of freedom: {sample1.Length + sample2.Length - 2}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Significant difference between means");
else
    Console.WriteLine("  Result: No significant difference between means");
```

**Hypotheses:**
- H₀: μ₁ = μ₂
- H₁: μ₁ ≠ μ₂

#### Unequal Variance (Welch's) t-Test

Use when variances may be different between groups:

```cs
// Returns the two-sided p-value (Welch's approximation for df)
double pValue = HypothesisTests.UnequalVarianceTtest(sample1, sample2);

Console.WriteLine($"Welch's t-test (unequal variances):");
Console.WriteLine($"  p-value: {pValue:F4}");
Console.WriteLine("  Use when variances appear different between groups");
```

#### Paired t-Test

For matched pairs or before/after comparisons:

```cs
double[] before = { 120, 135, 118, 142, 128 };
double[] after = { 115, 130, 112, 138, 125 };

// Returns the two-sided p-value
double pValue = HypothesisTests.PairedTtest(before, after);

Console.WriteLine($"Paired t-test:");
Console.WriteLine($"  Mean before: {Statistics.Mean(before):F1}");
Console.WriteLine($"  Mean after: {Statistics.Mean(after):F1}");
Console.WriteLine($"  Mean difference: {Statistics.Mean(before) - Statistics.Mean(after):F1}");
Console.WriteLine($"  p-value: {pValue:F4}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Treatment had a significant effect");
else
    Console.WriteLine("  Result: No significant treatment effect detected");
```

**Hypotheses:**
- H₀: μ_diff = 0
- H₁: μ_diff ≠ 0

## F-Test for Variance Equality

Tests if two populations have equal variances:

```cs
double[] sample1 = { 10, 12, 14, 16, 18 };
double[] sample2 = { 11, 13, 15, 17, 19, 21, 23 };

// Returns the p-value
double pValue = HypothesisTests.Ftest(sample1, sample2);

Console.WriteLine($"F-test for equal variances:");
Console.WriteLine($"  Variance 1: {Statistics.Variance(sample1):F2}");
Console.WriteLine($"  Variance 2: {Statistics.Variance(sample2):F2}");
Console.WriteLine($"  p-value: {pValue:F4}");
Console.WriteLine($"  df1 = {sample1.Length - 1}, df2 = {sample2.Length - 1}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Variances are significantly different");
else
    Console.WriteLine("  Result: No significant difference in variances");
```

**Hypotheses:**
- H₀: σ₁² = σ₂²
- H₁: σ₁² ≠ σ₂²

### F-Test for Nested Models

Compares restricted and full regression models:

```cs
// SSE values from model fitting
double sseRestricted = 150.0;  // SSE of restricted model
double sseFull = 120.0;        // SSE of full model
int dfRestricted = 47;         // n - k_restricted - 1
int dfFull = 45;               // n - k_full - 1

HypothesisTests.FtestModels(sseRestricted, sseFull, dfRestricted, dfFull,
                            out double fStat, out double pValue);

Console.WriteLine($"F-test for model comparison:");
Console.WriteLine($"  F-statistic: {fStat:F3}");
Console.WriteLine($"  p-value: {pValue:F4}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Full model is significantly better");
else
    Console.WriteLine("  Result: Models not significantly different");
```

## Normality Test

### Jarque-Bera Test

Tests if data follows a normal distribution using skewness and kurtosis:

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Returns the p-value (chi-squared with df=2)
double pValue = HypothesisTests.JarqueBeraTest(data);

Console.WriteLine($"Jarque-Bera normality test:");
Console.WriteLine($"  Skewness: {Statistics.Skewness(data):F3}");
Console.WriteLine($"  Kurtosis: {Statistics.Kurtosis(data):F3}");
Console.WriteLine($"  p-value: {pValue:F4}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Data is not normally distributed");
else
    Console.WriteLine("  Result: Cannot reject normality assumption");
```

**Hypotheses:**
- H₀: Data is normally distributed
- H₁: Data is not normally distributed

## Randomness and Independence Tests

### Wald-Wolfowitz Runs Test

Tests if a sequence is random (tests for independence and stationarity):

```cs
double[] sequence = { 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0 };

// Returns the two-sided p-value
double pValue = HypothesisTests.WaldWolfowitzTest(sequence);

Console.WriteLine($"Wald-Wolfowitz runs test:");
Console.WriteLine($"  p-value: {pValue:F4}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Sequence is not random");
else
    Console.WriteLine("  Result: Cannot reject randomness assumption");
```

### Ljung-Box Test

Tests for autocorrelation in time series data:

```cs
double[] timeSeries = { 12.5, 13.2, 11.8, 14.1, 12.9, 13.5, 12.2, 13.8, 14.5, 13.1 };
int lagMax = 5;  // Test lags 1 through 5

// Returns the p-value (chi-squared with df = lagMax)
double pValue = HypothesisTests.LjungBoxTest(timeSeries, lagMax);

Console.WriteLine($"Ljung-Box test for autocorrelation:");
Console.WriteLine($"  Lags tested: {lagMax}");
Console.WriteLine($"  p-value: {pValue:F4}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Significant autocorrelation detected");
else
    Console.WriteLine("  Result: No significant autocorrelation");
```

**Hypotheses:**
- H₀: No autocorrelation up to lag k
- H₁: Some autocorrelation exists

## Non-Parametric Tests

### Mann-Whitney U Test

Non-parametric alternative to two-sample t-test (tests if distributions differ):

```cs
double[] group1 = { 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42 };
double[] group2 = { 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40 };

// Note: First sample must have length ≤ second sample
// Combined samples must have length > 20 (here: 11 + 11 = 22)
// Returns the two-sided p-value
double pValue = HypothesisTests.MannWhitneyTest(group1, group2);

Console.WriteLine($"Mann-Whitney U test:");
Console.WriteLine($"  p-value: {pValue:F4}");
Console.WriteLine("  (Rank-based, no normality assumption required)");

if (pValue < 0.05)
    Console.WriteLine("  Result: Distributions differ significantly");
else
    Console.WriteLine("  Result: No significant difference in distributions");
```

**Hypotheses:**
- H₀: Distributions are equal
- H₁: Distributions differ

## Trend Tests

### Mann-Kendall Trend Test

Detects monotonic trends in time series (non-parametric):

```cs
double[] timeSeries = { 10, 12, 11, 15, 14, 18, 17, 21, 20, 24 };

// Requires sample size ≥ 10
// Returns the two-sided p-value
double pValue = HypothesisTests.MannKendallTest(timeSeries);

Console.WriteLine($"Mann-Kendall trend test:");
Console.WriteLine($"  p-value: {pValue:F4}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Significant monotonic trend detected");
else
    Console.WriteLine("  Result: No significant trend");
```

**Hypotheses:**
- H₀: No monotonic trend
- H₁: Monotonic trend exists

### Linear Trend Test

Tests for significant linear relationship (parametric):

```cs
double[] time = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
double[] values = { 10.5, 11.2, 12.8, 13.1, 14.5, 15.2, 16.1, 16.8, 17.5, 18.2 };

// Returns the two-sided p-value
double pValue = HypothesisTests.LinearTrendTest(time, values);

Console.WriteLine($"Linear trend test:");
Console.WriteLine($"  p-value: {pValue:F4}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Significant linear trend detected");
else
    Console.WriteLine("  Result: No significant linear trend");
```

## Unimodality Test

Tests if distribution has a single peak using Gaussian Mixture Model comparison:

```cs
double[] data = { 10, 11, 12, 13, 14, 15, 14, 13, 12, 11, 10 };

// Requires sample size ≥ 10
// Returns the p-value (likelihood ratio test with chi-squared df=3)
double pValue = HypothesisTests.UnimodalityTest(data);

Console.WriteLine($"Unimodality test:");
Console.WriteLine($"  p-value: {pValue:F4}");

if (pValue < 0.05)
    Console.WriteLine("  Result: Evidence of multimodality (bimodal or more)");
else
    Console.WriteLine("  Result: Cannot reject unimodality");
```

## Hydrological Applications

### Example 1: Testing for Stationarity in Annual Maximum Floods

A fundamental assumption in flood frequency analysis is that the annual maximum flood series is stationary (no trend over time). This example demonstrates how to test this assumption:

```cs
using Numerics.Data.Statistics;
using Numerics.Distributions;

// Annual maximum peak flows (cfs) from 1990-2019
double[] annualMaxFlows = {
    12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200, 10500, 19300,
    14800, 16200, 13100, 18500, 15600, 17800, 12800, 19100, 14300, 16500,
    13800, 18200, 15100, 17400, 12300, 19800, 14600, 16900, 13500, 18900
};

double[] years = Enumerable.Range(1990, annualMaxFlows.Length).Select(y => (double)y).ToArray();

Console.WriteLine("Stationarity Analysis of Annual Maximum Floods");
Console.WriteLine("=" + new string('=', 50));
Console.WriteLine($"Record length: {annualMaxFlows.Length} years (1990-2019)");
Console.WriteLine($"Mean: {Statistics.Mean(annualMaxFlows):F0} cfs");
Console.WriteLine($"Std Dev: {Statistics.StandardDeviation(annualMaxFlows):F0} cfs");

// Test 1: Mann-Kendall test for monotonic trend
double pMK = HypothesisTests.MannKendallTest(annualMaxFlows);
Console.WriteLine($"\nMann-Kendall Trend Test:");
Console.WriteLine($"  p-value: {pMK:F4}");
Console.WriteLine($"  Result: {(pMK < 0.05 ? "Trend detected - stationarity violated" : "No significant trend")}");

// Test 2: Linear trend test
double pLinear = HypothesisTests.LinearTrendTest(years, annualMaxFlows);
Console.WriteLine($"\nLinear Trend Test:");
Console.WriteLine($"  p-value: {pLinear:F4}");
Console.WriteLine($"  Result: {(pLinear < 0.05 ? "Linear trend detected" : "No significant linear trend")}");

// Test 3: Ljung-Box test for serial correlation
double pLB = HypothesisTests.LjungBoxTest(annualMaxFlows, lagMax: 5);
Console.WriteLine($"\nLjung-Box Autocorrelation Test (lag=5):");
Console.WriteLine($"  p-value: {pLB:F4}");
Console.WriteLine($"  Result: {(pLB < 0.05 ? "Autocorrelation detected" : "No significant autocorrelation")}");

// Overall assessment
Console.WriteLine("\n" + new string('-', 50));
bool isStationary = pMK >= 0.05 && pLinear >= 0.05 && pLB >= 0.05;
if (isStationary)
{
    Console.WriteLine("Assessment: Data appears stationary - proceed with frequency analysis");
}
else
{
    Console.WriteLine("Assessment: Potential non-stationarity detected");
    Console.WriteLine("Consider: detrending, using shorter record, or non-stationary methods");
}
```

### Example 2: Comparing Flood Records Between Two Periods

Test whether flood characteristics have changed between historical and recent periods:

```cs
// Split record into two periods
double[] period1 = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200, 10500, 19300 };  // 1990-1999
double[] period2 = { 14800, 16200, 13100, 18500, 15600, 17800, 12800, 19100, 14300, 16500 };  // 2000-2009

Console.WriteLine("Comparison of Flood Characteristics: 1990s vs 2000s");
Console.WriteLine("=" + new string('=', 55));

Console.WriteLine($"\nPeriod 1 (1990-1999):");
Console.WriteLine($"  Mean: {Statistics.Mean(period1):F0} cfs");
Console.WriteLine($"  Std Dev: {Statistics.StandardDeviation(period1):F0} cfs");

Console.WriteLine($"\nPeriod 2 (2000-2009):");
Console.WriteLine($"  Mean: {Statistics.Mean(period2):F0} cfs");
Console.WriteLine($"  Std Dev: {Statistics.StandardDeviation(period2):F0} cfs");

// Test for difference in means
double pMean = HypothesisTests.EqualVarianceTtest(period1, period2);
Console.WriteLine($"\nt-test for difference in means:");
Console.WriteLine($"  p-value: {pMean:F4}");
Console.WriteLine($"  Result: {(pMean < 0.05 ? "Means differ significantly" : "No significant difference in means")}");

// Test for difference in variances
double pVar = HypothesisTests.Ftest(period1, period2);
Console.WriteLine($"\nF-test for difference in variances:");
Console.WriteLine($"  p-value: {pVar:F4}");
Console.WriteLine($"  Result: {(pVar < 0.05 ? "Variances differ significantly" : "No significant difference in variances")}");

// Interpretation for flood frequency analysis
Console.WriteLine("\n" + new string('-', 55));
if (pMean >= 0.05 && pVar >= 0.05)
    Console.WriteLine("Conclusion: Records can be combined for frequency analysis");
else
    Console.WriteLine("Conclusion: Consider analyzing periods separately or investigating causes");
```

### Example 3: Testing Normality of Log-Transformed Flood Data

Log-Pearson Type III analysis assumes the log-transformed data follows a Pearson Type III distribution. Testing normality of log-transformed data can indicate if a simpler Log-Normal model might suffice:

```cs
double[] annualPeaks = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200, 10500, 19300,
                         14800, 16200, 13100, 18500, 15600, 17800, 12800, 19100, 14300, 16500 };

// Log-transform the data
double[] logPeaks = annualPeaks.Select(x => Math.Log10(x)).ToArray();

Console.WriteLine("Normality Test for Log-Transformed Annual Peak Flows");
Console.WriteLine("=" + new string('=', 55));

Console.WriteLine($"\nLog-transformed data statistics:");
Console.WriteLine($"  Mean: {Statistics.Mean(logPeaks):F4}");
Console.WriteLine($"  Std Dev: {Statistics.StandardDeviation(logPeaks):F4}");
Console.WriteLine($"  Skewness: {Statistics.Skewness(logPeaks):F4}");
Console.WriteLine($"  Kurtosis: {Statistics.Kurtosis(logPeaks):F4}");

// Jarque-Bera test for normality
double pJB = HypothesisTests.JarqueBeraTest(logPeaks);

Console.WriteLine($"\nJarque-Bera Normality Test:");
Console.WriteLine($"  p-value: {pJB:F4}");

if (pJB >= 0.05)
{
    Console.WriteLine("  Result: Cannot reject normality of log-transformed data");
    Console.WriteLine("  Implication: Log-Normal distribution may be appropriate");
}
else
{
    Console.WriteLine("  Result: Log-transformed data departs from normality");
    Console.WriteLine("  Implication: Consider Log-Pearson Type III or GEV distribution");
}

// Additional check: skewness significance
double skew = Statistics.Skewness(logPeaks);
if (Math.Abs(skew) < 0.5)
{
    Console.WriteLine($"\n  Note: Skewness ({skew:F3}) is small - Log-Normal may be adequate");
}
else
{
    Console.WriteLine($"\n  Note: Skewness ({skew:F3}) is substantial - LP3 recommended");
}
```

## Test Selection Guide

| Question | Test | Notes |
|----------|------|-------|
| One sample mean vs. hypothesized value | `OneSampleTtest` | Parametric, assumes normality |
| Two independent means (equal variance) | `EqualVarianceTtest` | Use F-test first to check variance equality |
| Two independent means (unequal variance) | `UnequalVarianceTtest` | Welch's t-test, more robust |
| Two paired measurements | `PairedTtest` | Before/after, matched pairs |
| Two variances equal? | `Ftest` | Sensitive to non-normality |
| Model comparison | `FtestModels` | Nested regression models |
| Data normally distributed? | `JarqueBeraTest` | Uses skewness and kurtosis |
| Sequence random? | `WaldWolfowitzTest` | Independence, stationarity |
| Time series autocorrelation? | `LjungBoxTest` | Tests multiple lags |
| Distribution difference? | `MannWhitneyTest` | Non-parametric, rank-based |
| Monotonic trend? | `MannKendallTest` | Non-parametric trend test |
| Linear trend? | `LinearTrendTest` | Parametric trend test |
| Unimodal distribution? | `UnimodalityTest` | GMM-based comparison |

## Best Practices

1. **Choose significance level a priori** - Typically α = 0.05, but consider α = 0.10 for exploratory analysis
2. **Check test assumptions** - Many tests assume normality or independence
3. **Use non-parametric alternatives** - When assumptions are violated (e.g., Mann-Whitney instead of t-test)
4. **Report p-values, not just decisions** - p = 0.049 and p = 0.051 are essentially equivalent
5. **Consider multiple testing corrections** - Bonferroni or FDR when performing many tests
6. **Report effect sizes** - Statistical significance ≠ practical significance
7. **Visualize data** - Plots often reveal more than hypothesis tests

## Important Notes

- All tests return **p-values**, not test statistics
- Two-sided tests are used by default
- Sample size requirements vary by test (see individual test documentation)
- For hydrological applications, consider the impact of outliers and measurement uncertainty

---

## References

<a id="1">[1]</a> Helsel, D. R., Hirsch, R. M., Ryberg, K. R., Archfield, S. A., & Gilroy, E. J. (2020). *Statistical Methods in Water Resources*. U.S. Geological Survey Techniques and Methods, Book 4, Chapter A3.

<a id="2">[2]</a> Hirsch, R. M., Slack, J. R., & Smith, R. A. (1982). Techniques of trend analysis for monthly water quality data. *Water Resources Research*, 18(1), 107-121.

---

[← Previous: Goodness-of-Fit](goodness-of-fit.md) | [Back to Index](../index.md) | [Next: Convergence Diagnostics →](../sampling/convergence-diagnostics.md)
