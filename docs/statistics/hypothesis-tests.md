# Hypothesis Tests

[← Previous: Goodness-of-Fit](goodness-of-fit.md) | [Back to Index](../index.md) | [Next: Convergence Diagnostics →](../sampling/convergence-diagnostics.md)

The ***Numerics*** library provides statistical hypothesis tests for comparing samples, testing distributions, and detecting trends. These tests are essential for data analysis, model validation, and quality control.

## t-Tests

### One-Sample t-Test

Test if sample mean differs from hypothesized value:

```cs
using Numerics.Data.Statistics;

double[] sample = { 12.5, 13.2, 11.8, 14.1, 12.9, 13.5, 12.2, 13.8 };
double mu0 = 12.0;  // Hypothesized mean

// Compute t-statistic
double t = HypothesisTests.OneSampleTtest(sample, mu0);

Console.WriteLine($"One-sample t-test:");
Console.WriteLine($"  t-statistic: {t:F3}");
Console.WriteLine($"  Sample mean: {Statistics.Mean(sample):F2}");
Console.WriteLine($"  Hypothesized mean: {mu0:F2}");

// Compare with critical value
int df = sample.Length - 1;
Console.WriteLine($"  Degrees of freedom: {df}");
Console.WriteLine("  If |t| > t_critical (e.g., 2.365 for α=0.05, df=7), reject H₀");
```

**Hypotheses:**
- H₀: μ = μ₀
- H₁: μ ≠ μ₀ (two-tailed)

### Two-Sample t-Tests

#### Equal Variance (Pooled) t-Test

```cs
double[] sample1 = { 12.5, 13.2, 11.8, 14.1, 12.9 };
double[] sample2 = { 15.3, 14.8, 15.9, 14.5, 15.1 };

// Test if means are equal (assuming equal variances)
double t = HypothesisTests.EqualVarianceTtest(sample1, sample2);

Console.WriteLine($"Equal variance t-test:");
Console.WriteLine($"  t-statistic: {t:F3}");
Console.WriteLine($"  Sample 1 mean: {Statistics.Mean(sample1):F2}");
Console.WriteLine($"  Sample 2 mean: {Statistics.Mean(sample2):F2}");

int df = sample1.Length + sample2.Length - 2;
Console.WriteLine($"  Degrees of freedom: {df}");
```

**Hypotheses:**
- H₀: μ₁ = μ₂
- H₁: μ₁ ≠ μ₂

#### Unequal Variance (Welch's) t-Test

```cs
// Test if means are equal (not assuming equal variances)
double t_welch = HypothesisTests.UnequalVarianceTtest(sample1, sample2);

Console.WriteLine($"\nWelch's t-test:");
Console.WriteLine($"  t-statistic: {t_welch:F3}");
Console.WriteLine("  Use when variances appear different");
```

#### Paired t-Test

For before/after or matched pairs:

```cs
double[] before = { 120, 135, 118, 142, 128 };
double[] after = { 115, 130, 112, 138, 125 };

// Test if treatment had effect
double t_paired = HypothesisTests.PairedTtest(before, after);

Console.WriteLine($"Paired t-test:");
Console.WriteLine($"  t-statistic: {t_paired:F3}");
Console.WriteLine($"  Mean difference: {Statistics.Mean(before) - Statistics.Mean(after):F2}");

// Differences
double[] diffs = new double[before.Length];
for (int i = 0; i < before.Length; i++)
    diffs[i] = before[i] - after[i];

Console.WriteLine($"  SE of difference: {Statistics.StandardDeviation(diffs) / Math.Sqrt(diffs.Length):F2}");
```

**Hypotheses:**
- H₀: μ_diff = 0
- H₁: μ_diff ≠ 0

## F-Test

Test equality of variances:

```cs
double[] sample1 = { 10, 12, 14, 16, 18 };
double[] sample2 = { 11, 13, 15, 17, 19, 21, 23 };

// Test if variances are equal
double f = HypothesisTests.Ftest(sample1, sample2);

Console.WriteLine($"F-test for equal variances:");
Console.WriteLine($"  F-statistic: {f:F3}");
Console.WriteLine($"  Variance 1: {Statistics.Variance(sample1):F2}");
Console.WriteLine($"  Variance 2: {Statistics.Variance(sample2):F2}");
Console.WriteLine($"  df1 = {sample1.Length - 1}, df2 = {sample2.Length - 1}");
```

**Hypotheses:**
- H₀: σ₁² = σ₂²
- H₁: σ₁² ≠ σ₂²

### F-Test for Nested Models

Compare restricted and full models:

```cs
// Example: Testing if additional predictors improve model
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
    Console.WriteLine("  Reject H₀: Full model is significantly better");
else
    Console.WriteLine("  Fail to reject H₀: Models not significantly different");
```

## Normality Tests

### Jarque-Bera Test

Tests if data follows normal distribution using skewness and kurtosis:

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };

// Test for normality
double jb = HypothesisTests.JarqueBeraTest(data);

Console.WriteLine($"Jarque-Bera normality test:");
Console.WriteLine($"  JB statistic: {jb:F3}");
Console.WriteLine($"  Skewness: {Statistics.Skewness(data):F3}");
Console.WriteLine($"  Kurtosis: {Statistics.Kurtosis(data):F3}");
Console.WriteLine("  Critical value (α=0.05): 5.99 (χ² with df=2)");

if (jb < 5.99)
    Console.WriteLine("  Fail to reject H₀: Data appears normally distributed");
else
    Console.WriteLine("  Reject H₀: Data not normally distributed");
```

**Hypotheses:**
- H₀: Data is normally distributed
- H₁: Data is not normal

## Randomness Tests

### Wald-Wolfowitz Runs Test

Tests if sequence is random:

```cs
double[] sequence = { 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0 };

double z = HypothesisTests.WaldWolfowitzTest(sequence);

Console.WriteLine($"Wald-Wolfowitz runs test:");
Console.WriteLine($"  z-statistic: {z:F3}");

if (Math.Abs(z) > 1.96)
    Console.WriteLine("  Reject H₀: Sequence is not random");
else
    Console.WriteLine("  Fail to reject H₀: Sequence appears random");
```

### Ljung-Box Test

Tests for autocorrelation in time series:

```cs
double[] timeSeries = { 12.5, 13.2, 11.8, 14.1, 12.9, 13.5, 12.2, 13.8, 14.5, 13.1 };
int lagMax = 5;  // Test lags 1 through 5

double q = HypothesisTests.LjungBoxTest(timeSeries, lagMax);

Console.WriteLine($"Ljung-Box test for autocorrelation:");
Console.WriteLine($"  Q-statistic: {q:F3}");
Console.WriteLine($"  Lags tested: {lagMax}");
Console.WriteLine($"  Critical value (α=0.05, df={lagMax}): ~11.07");

if (q > 11.07)
    Console.WriteLine("  Reject H₀: Significant autocorrelation present");
else
    Console.WriteLine("  Fail to reject H₀: No significant autocorrelation");
```

**Hypotheses:**
- H₀: No autocorrelation up to lag k
- H₁: Some autocorrelation exists

## Non-Parametric Tests

### Mann-Whitney U Test

Non-parametric alternative to two-sample t-test:

```cs
double[] group1 = { 12, 15, 18, 21, 24 };
double[] group2 = { 10, 13, 16, 19, 22, 25 };

double u = HypothesisTests.MannWhitneyTest(group1, group2);

Console.WriteLine($"Mann-Whitney U test:");
Console.WriteLine($"  U-statistic: {u:F3}");
Console.WriteLine("  Tests if distributions are different (rank-based)");
Console.WriteLine("  No assumption of normality required");
```

**Hypotheses:**
- H₀: Distributions are equal
- H₁: Distributions differ

## Trend Tests

### Mann-Kendall Trend Test

Detects monotonic trends in time series:

```cs
double[] timeSeries = { 10, 12, 11, 15, 14, 18, 17, 21, 20, 24 };

double s = HypothesisTests.MannKendallTest(timeSeries);

Console.WriteLine($"Mann-Kendall trend test:");
Console.WriteLine($"  S-statistic: {s:F3}");

if (s > 1.96)
    Console.WriteLine("  Significant increasing trend detected");
else if (s < -1.96)
    Console.WriteLine("  Significant decreasing trend detected");
else
    Console.WriteLine("  No significant trend");
```

**Hypotheses:**
- H₀: No monotonic trend
- H₁: Monotonic trend exists

### Linear Trend Test

Tests for linear relationship:

```cs
double[] time = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
double[] values = { 10.5, 11.2, 12.8, 13.1, 14.5, 15.2, 16.1, 16.8, 17.5, 18.2 };

double t = HypothesisTests.LinearTrendTest(time, values);

Console.WriteLine($"Linear trend test:");
Console.WriteLine($"  t-statistic: {t:F3}");
Console.WriteLine("  Tests if slope is significantly different from zero");

if (Math.Abs(t) > 2.306)  // Critical value for df=8, α=0.05
    Console.WriteLine("  Significant linear trend detected");
else
    Console.WriteLine("  No significant linear trend");
```

## Unimodality Test

Tests if distribution has single peak:

```cs
double[] data = { 10, 11, 12, 13, 14, 15, 14, 13, 12, 11, 10 };

double u = HypothesisTests.UnimodalityTest(data);

Console.WriteLine($"Unimodality test:");
Console.WriteLine($"  Test statistic: {u:F3}");

if (u < -1.96)
    Console.WriteLine("  Reject H₀: Data is multimodal");
else
    Console.WriteLine("  Fail to reject H₀: Data appears unimodal");
```

## Practical Examples

### Example 1: Comparing Treatment Groups

```cs
using Numerics.Data.Statistics;

// Control and treatment groups
double[] control = { 120, 135, 118, 142, 128, 133, 125, 138 };
double[] treatment = { 115, 125, 110, 130, 120, 128, 118, 130 };

Console.WriteLine("Treatment Comparison Study");
Console.WriteLine("=" + new string('=', 50));

// Descriptive statistics
Console.WriteLine($"\nControl: mean={Statistics.Mean(control):F1}, " +
                 $"SD={Statistics.StandardDeviation(control):F1}");
Console.WriteLine($"Treatment: mean={Statistics.Mean(treatment):F1}, " +
                 $"SD={Statistics.StandardDeviation(treatment):F1}");

// Test for equal variances
double f = HypothesisTests.Ftest(control, treatment);
Console.WriteLine($"\nF-test for equal variances: F={f:F3}");

// Choose appropriate t-test
double t;
if (f < 4.0)  // Approximate F-critical value
{
    t = HypothesisTests.EqualVarianceTtest(control, treatment);
    Console.WriteLine($"Equal variance t-test: t={t:F3}");
}
else
{
    t = HypothesisTests.UnequalVarianceTtest(control, treatment);
    Console.WriteLine($"Unequal variance t-test: t={t:F3}");
}

if (Math.Abs(t) > 2.145)  // Approximate critical value
    Console.WriteLine("Conclusion: Significant difference detected (p < 0.05)");
else
    Console.WriteLine("Conclusion: No significant difference (p ≥ 0.05)");
```

### Example 2: Time Series Analysis

```cs
double[] monthlyData = { 125, 130, 135, 132, 138, 145, 142, 148, 155, 152, 158, 165 };

Console.WriteLine("Time Series Analysis");
Console.WriteLine("=" + new string('=', 50));

// Test for trend
double[] months = Enumerable.Range(1, monthlyData.Length).Select(i => (double)i).ToArray();
double tTrend = HypothesisTests.LinearTrendTest(months, monthlyData);

Console.WriteLine($"\nLinear trend test: t={tTrend:F3}");
if (Math.Abs(tTrend) > 2.228)  // Critical value for df=10
    Console.WriteLine("Significant trend detected");

// Test for autocorrelation
double q = HypothesisTests.LjungBoxTest(monthlyData, lagMax: 3);

Console.WriteLine($"\nLjung-Box test (lag 3): Q={q:F3}");
if (q > 7.815)  // Chi-squared critical value
    Console.WriteLine("Significant autocorrelation detected");
else
    Console.WriteLine("No significant autocorrelation");

// Mann-Kendall for monotonic trend
double s = HypothesisTests.MannKendallTest(monthlyData);

Console.WriteLine($"\nMann-Kendall test: S={s:F3}");
if (s > 1.96)
    Console.WriteLine("Significant increasing trend (non-parametric)");
```

### Example 3: Quality Control

```cs
// Historical process mean
double mu0 = 50.0;

// New sample
double[] newSample = { 51.2, 52.1, 49.8, 51.5, 50.9, 52.3, 51.0, 50.5 };

Console.WriteLine("Quality Control Check");
Console.WriteLine("=" + new string('=', 50));

// One-sample t-test
double t = HypothesisTests.OneSampleTtest(newSample, mu0);

Console.WriteLine($"\nHistorical mean: {mu0:F1}");
Console.WriteLine($"Current sample mean: {Statistics.Mean(newSample):F2}");
Console.WriteLine($"t-statistic: {t:F3}");
Console.WriteLine($"Critical value (two-tailed, α=0.05): ±2.365");

if (Math.Abs(t) > 2.365)
{
    Console.WriteLine("\nProcess has shifted significantly!");
    Console.WriteLine("Action: Investigate and adjust process");
}
else
{
    Console.WriteLine("\nProcess remains in control");
}
```

## Interpreting Results

### p-values
- p < 0.01: Very strong evidence against H₀
- p < 0.05: Strong evidence against H₀
- p < 0.10: Weak evidence against H₀
- p ≥ 0.10: Insufficient evidence to reject H₀

### Effect Size
Statistical significance ≠ practical significance. Consider:
- Cohen's d for t-tests: (mean difference) / pooled SD
- Small: d = 0.2, Medium: d = 0.5, Large: d = 0.8

### Power
- Probability of detecting true effect
- Influenced by sample size, effect size, α level
- Power ≥ 0.80 typically desired

## Best Practices

1. **Check assumptions** before parametric tests (normality, equal variance)
2. **Use non-parametric tests** when assumptions violated
3. **Report effect sizes** along with p-values
4. **Consider multiple testing** corrections if doing many tests
5. **Visualize data** before and after testing
6. **Understand context** - statistical vs practical significance

## Test Selection Guide

| Question | Test |
|----------|------|
| One sample mean vs. value | One-sample t-test |
| Two independent means | Two-sample t-test (or Mann-Whitney) |
| Two paired measurements | Paired t-test |
| Two variances | F-test |
| Normality | Jarque-Bera |
| Trend existence | Mann-Kendall |
| Linear relationship | Linear trend test |
| Autocorrelation | Ljung-Box |
| Randomness | Wald-Wolfowitz |

---

[← Previous: Goodness-of-Fit](goodness-of-fit.md) | [Back to Index](../index.md) | [Next: Convergence Diagnostics →](../sampling/convergence-diagnostics.md)
