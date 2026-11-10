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

The source code uses the corrected two-pass algorithm of Chan, Golub and LeVeque (1983): it sums $(x_i - \bar{x})^2$ about the computed mean and subtracts $\left(\sum (x_i - \bar{x})\right)^2 / n$, which removes the rounding error of the mean to first order. It keeps its accuracy when the mean is large relative to the spread: on 28,317 values $10^6 + (i \bmod 997) \cdot 10^{-6}$ the error is 5.5E-14 relative, where a running sum of squares was 2.5E-7 off. `Covariance`, `PopulationCovariance` and the weighted `Variance` use the same correction.

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

### Median

`Statistics.Median` returns the middle order statistic for an odd count and the midpoint of the two middle order statistics for an even count:

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9, 11.2, 14.5 };
double median = Statistics.Median(data);  // 13.3, the midpoint of 12.9 and 13.7
```

The even-count midpoint is computed as $\tfrac{1}{2}(x_{lo} + x_{hi})$, which is correctly rounded unless the sum overflows; when it overflows, the midpoint is $\tfrac{1}{2}x_{lo} + \tfrac{1}{2}x_{hi}$, so the median of finite values is always finite. Any NaN value makes the median NaN, as in `numpy.median`. `Percentile(data, 0.5)` returns the same value up to its last bit, because Type 7 interpolation computes it as $x_{lo} + \tfrac{1}{2}(x_{hi} - x_{lo})$ with two roundings. Pass `dataIsSorted: true` to skip the sorting copy.

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

### Compensated Summation

`Tools.Sum` adds the values in order, so its rounding error can grow with the number of terms and is large when terms cancel. `Tools.CompensatedSum` uses Neumaier's improvement of Kahan summation [[8]](#8): a running correction collects the exact rounding error of each addition and is added back at the end. Its error is bounded by 2u|s| + O(n u²) Σ|x_i|, so to first order it does not grow with n, and unlike plain Kahan summation it stays accurate when a term is larger than the running sum.

```cs
using Numerics;

double[] x = { 1E100, 1.0, -1E100 };
double plain = Tools.Sum(x);                       // 0
double exact = Tools.CompensatedSum(x);            // 1, as Python math.fsum
double kahanCase = Tools.CompensatedSum(new[] { 1.0, 1E100, 1.0, -1E100 });   // 2 (Kahan gives 0)

// The values marked 1 (or 0 with useComplement), as Tools.Sum selects them
double selected = Tools.CompensatedSum(new[] { 1E100, 1.0, -1E100, 7.0 }, new[] { 1, 1, 1, 0 });   // 1

// A gather sum of values[indices[i]] for i in [start, start + count), without copying
double gathered = Tools.CompensatedSum(new[] { 1E100, 1.0, -1E100, 3.0 }, new[] { 3, 0, 1, 2 }, 1, 3);   // 1

// Running totals: entry i is the compensated sum of the first i + 1 values
double[] running = Tools.CompensatedCumulativeSum(new[] { 1E100, 1.0, -1E100, 0.1 });   // { 1E100, 1E100, 1, 1.1 }

// Classic Kahan summation instead, bit for bit, for reproducing results of the classic Kahan loop
double kahan = Tools.CompensatedSum(new[] { 1.0, 1E100, 1.0, -1E100 }, kahan: true);       // 0
double[] kahanRunning = Tools.CompensatedCumulativeSum(new[] { 0.1, 0.2, 0.3 }, kahan: true);
```

With `kahan: true` an infinite final sum is returned as is, but any value added after the sum becomes infinite gives NaN, as in the classic Kahan loop; the default Neumaier form returns the IEEE sum in both cases. A NaN, or +∞ and −∞ together, gives NaN, and a single infinity is returned with its sign, as in IEEE arithmetic; an intermediate overflow of finite terms returns ±∞. Python's `math.fsum` gives the same NaN or infinity, but raises for +∞ with −∞ and for an intermediate overflow. As for `Tools.Sum`, an empty list (or an empty gather range) returns NaN, and the indicator form returns NaN when the indicators differ in length from the values. The gather form throws `ArgumentOutOfRangeException` for a range outside the indices or an index outside the values.

The step of each method is public, for loops that form their terms on the fly, without collecting them in a list. `Tools.NeumaierAdd(ref sum, ref compensation, value)` adds one value and `Tools.NeumaierResult(sum, compensation)` returns the result; started from zero and called in order, they give `Tools.CompensatedSum` bit for bit. `Tools.KahanAdd` is the classic Kahan step, whose running sum is the `kahan: true` result bit for bit.

```cs
using Numerics;

double sum = 0, compensation = 0;
for (int i = 1; i <= 1000; i++)
    Tools.NeumaierAdd(ref sum, ref compensation, 1.0 / i);
double harmonic = Tools.NeumaierResult(sum, compensation);   // 7.485470860550345, as math.fsum

double t = 0, c = 0;
foreach (double v in new[] { 0.1, 0.2, 0.3 }) Tools.KahanAdd(ref t, ref c, v);   // t is the Kahan sum
```

### Log-Sum-Exp

`Tools.LogSumExp` evaluates log(Σ exp(x_i)) of log values without overflow, as `scipy.special.logsumexp` does. With m the largest value and k the number of values equal to it, the result is log1p(s / k) + log(k) + m, where s is the compensated sum of exp(x_i − m) over the other values. Separating the largest terms and using `Tools.Log1p` keeps a result near zero accurate (Blanchard, Higham and Higham, 2021), where m + log(Σ exp(x_i − m)) loses its relative accuracy. The results agree with SciPy to 2 ulps of the larger of the result and m.

```cs
using System;
using Numerics;

double pair = Tools.LogSumExp(1000.0, 1000.0);                 // 1000 + ln 2, without overflow
double small = Tools.LogSumExp(new[] { 0.0, -30.0 });         // 9.357622968839737e-14 = log1p(exp(-30))
double empty = Tools.LogSumExp(new double[0]);                // -∞, the log of an empty sum

// One value per column (axis 0) or per row (axis 1), each bitwise the list form of that column or row
var logs = new double[,] { { Math.Log(0.25), Math.Log(0.5) }, { Math.Log(0.75), Math.Log(0.5) } };
double[] columns = Tools.LogSumExp(logs, 0);                  // about { 0, 0 }
double[] rows = Tools.LogSumExp(logs, 1);                     // about { log(0.75), log(1.25) }
```

The special values follow SciPy: an empty list or a list of −∞ values returns −∞, any +∞ returns +∞, and NaN propagates (also beside +∞).

### Hypotenuse

`Tools.Hypot(a, b)` returns √(a² + b²) without the overflow or underflow of squaring the arguments. The larger magnitude m is factored out and the result is m·√(1 + r²) with r = min/max ≤ 1, within 2 ulps of the exact value. Special values follow IEEE 754 and C99 `hypot`: an infinite argument gives +∞ even when the other is NaN, otherwise a NaN argument gives NaN, and (0, 0) gives 0.

```cs
using Numerics;

double big = Tools.Hypot(3e200, 4e200);      // about 5e200; Math.Sqrt(a * a + b * b) overflows to ∞
double tiny = Tools.Hypot(3e-200, 4e-200);   // 5e-200; Math.Sqrt(a * a + b * b) underflows to 0
double inf = Tools.Hypot(double.PositiveInfinity, double.NaN);   // ∞
```

### Column Statistics

For a matrix with one row per member (for example one row per year) and one column per variable (for example one column per grid cell), the column members evaluate a statistic for every column in one call, in parallel over the columns. NaN and infinite entries are excluded from each column and counted in the optional `excluded` array, and every column's value is bitwise the one-sample member applied to that column's finite entries:

```cs
double[,] maxima =
{
    { 3.2, 101.5 },
    { 1.7, 99.25 },
    { 4.4, double.NaN },
    { 2.9, 98.5 },
    { 5.1, 102.0 },
};
double[] means = Statistics.ColumnMeans(maxima);
double[] variances = Statistics.ColumnVariances(maxima, out int[] excluded);  // { 1.763, 2.890625 }, excluded { 0, 1 }
double[] deviations = Statistics.ColumnStandardDeviations(maxima);
double[,] lmoments = Statistics.ColumnLinearMoments(maxima);  // 4 rows: λ1, λ2, τ3, τ4
double lcv = lmoments[1, 0] / lmoments[0, 0];  // bitwise Statistics.LinearCoefficientOfVariation of column 0
```

A column with too few finite entries gives NaN: fewer than two for the variance and standard deviation, fewer than four for the L-moments. `ColumnPercentiles` follows the same rules.

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

`LinearMoments(data, nmom)` returns the first `nmom` L-moments $\lambda_1, \lambda_2, \tau_3, \ldots, \tau_{nmom}$ with the same unbiased estimator (the SAMLMR estimator of Hosking's LMOMENTS package), extended through $\lambda_{r+1} = \sum_{k=0}^{r} p^*_{r,k}\, b_k$ with the shifted Legendre coefficients $p^*_{r,k} = (-1)^{r-k}(r+k)!/((k!)^2(r-k)!)$ [[1]](#1). Its first four values are bitwise those of `LinearMoments(data)`, and $\tau_5$ is the extra ratio used in regional frequency analysis. Every value is NaN when the sample has fewer than `nmom` entries. The error of $	au_r$ grows by roughly a factor of 10 to 30 per order (about $10^{-10}$ at $	au_{10}$ and $10^{-3}$ at $	au_{20}$ for lognormal samples), so keep `nmom` at 10 or below for accurate ratios; at very high orders the coefficient products overflow and the ratios become infinite or NaN.

```cs
double[] moments = Statistics.LinearMoments(data, 5);   // λ1, λ2, τ3, τ4, τ5
```

The L-coefficient of variation (L-CV), $\tau = \lambda_2 / \lambda_1$ [[1]](#1), is returned by `LinearCoefficientOfVariation`. It is bitwise `LinearMoments(data)[1] / LinearMoments(data)[0]`, and like `LinearMoments` it is NaN for fewer than four values.

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

// L-CV, τ = λ₂ / λ₁
double lcv = Statistics.LinearCoefficientOfVariation(data);
Console.WriteLine($"  τ (L-CV): {lcv:F4}");

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

With a plotting-position formula, the i-th order statistic sits at its plotting position $p_i = (i - \alpha)/(n + 1 - 2\alpha)$ and the percentile interpolates linearly between them, held at the end values. Weibull is Hyndman and Fan type 6, Hazen type 5 and Kaplan-Meier type 4. For Weibull, Blom and Hazen the interpolation weight is formed from $k$ and the two bracketing positions, $(k - p_j)/(p_{j+1} - p_j)$, which keeps bits of $k$ that the rank $(n + 1 - 2\alpha)k + \alpha$ loses; the other formulas use the fractional part of that rank.

```cs
double p90 = Statistics.Percentile(data, 0.9, PlottingPositions.PlottingPostionType.Weibull);
```

With weights and a plotting-position formula, point i of the sorted data sits at the weighted position $(c_{i-1} + (1-\alpha) w_i) / (W + (1 - 2\alpha)\bar{w})$, where $c_{i-1}$ is the weight before it, $W$ the total and $\bar{w}$ the mean weight; for Hazen this is $(c_{i-1} + w_i/2)/W$. The percentile interpolates linearly between the positions and is held at the end values. Zero weights keep their position. With sorted data no array is allocated:

```cs
double[] flows = { 120, 150, 180, 260 };          // ascending
double[] weights = { 0.4, 0.3, 0.2, 0.1 };
double p50 = Statistics.Percentile(flows, 0.5, weights, PlottingPositions.PlottingPostionType.Hazen, dataIsSorted: true);   // 145.71, between 120 at 0.2 and 150 at 0.55
```

For several probabilities, the list overload validates and sorts the sample once and serves all of them in one walk over the positions, in any order; each value is bitwise the scalar call:

```cs
double[] p = Statistics.Percentile(flows, new[] { 0.9, 0.1, 0.5 }, weights, PlottingPositions.PlottingPostionType.Hazen);   // 233.33, 120, 145.71
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

Ties are exact equality and concordance is decided by comparing values, so infinite values and values of any magnitude are ranked exactly; a NaN in either sample gives NaN, as `scipy.stats.kendalltau` does by default. For samples of 48 or more pairs, `KendallsTau` obtains the counts in $O(n \log n)$ time with Knight's algorithm [[5]](#5): sort the pairs by $x$ (ties by $y$), count the tied pairs, then count the discordant pairs as the inversions removed by a merge sort of $y$. Smaller samples use the $O(n^2)$ pair loop, which gives the same integer counts. The counts are 64-bit, so there is no practical limit on the sample size. At $n = 1{,}344$ the fast path is about 36 times faster than the pair loop, and at $n = 10{,}000$ about 170 times faster.

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

The Pearson sums use the corrected two-pass algorithm of `Statistics.Variance` and `Statistics.Covariance`: each sum of products of deviations is reduced by the product of the deviation sums divided by $n$, which removes the rounding error of the computed means [[9]](#9). The coefficient therefore stays accurate when the means are large relative to the spread, for example values near $10^9$ that vary in the first decimal.

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

// Ranks with ties reported. Each entry of "ties" holds a tie run's length minus one at the
// run's last sorted position (zero elsewhere), so a tied group of any size reports a value
// greater than zero.
double[] dataCopy2 = (double[])data.Clone();
double[] ranks2 = Statistics.RanksInPlace(dataCopy2, out double[] ties);

Console.WriteLine($"\nNumber of tied groups: {ties.Count(t => t > 0)}");
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

// Tippecanoe River near Delphi, IN: 48 years of annual peak streamflow (cfs)
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

## Binomial Proportion Interval

`Statistics.WilsonScoreInterval(successes, trials, confidenceLevel)` returns the Wilson score interval [[10]](#10) for a binomial proportion. With $\hat{p} = s/n$ and $z$ the standard normal quantile at $0.5 + \text{level}/2$:

```math
\frac{\hat{p} + \frac{z^2}{2n}}{1 + \frac{z^2}{n}} \pm \frac{z}{1 + \frac{z^2}{n}} \sqrt{\frac{\hat{p}(1 - \hat{p})}{n} + \frac{z^2}{4n^2}}
```

Unlike the normal-approximation (Wald) interval, it stays inside [0, 1] and has a positive width when there are no successes or no failures. The bounds are clamped to [0, 1], and the counts may be fractional, for example weighted counts. A non-positive or non-finite number of trials, successes outside [0, trials], or a level outside (0, 1) throws `ArgumentOutOfRangeException`.

```cs
using Numerics.Data.Statistics;

var (lower, upper) = Statistics.WilsonScoreInterval(successes: 7, trials: 20, confidenceLevel: 0.9);
Console.WriteLine($"[{lower:F4}, {upper:F4}]");   // [0.2023, 0.5335]

var (lower0, upper0) = Statistics.WilsonScoreInterval(0, 20, 0.9);
Console.WriteLine($"[{lower0:F4}, {upper0:F4}]"); // [0.0000, 0.1192]
```

## Effective Sample Size of Weights

`Statistics.EffectiveSampleSize(weights)` returns the Kish effective sample size [[7]](#7) of a weight vector, the number of equally weighted observations that carry the same information:

```math
n_{\text{eff}} = \frac{\left(\sum_i w_i\right)^2}{\sum_i w_i^2}
```

Equal positive weights give the count, zero weights drop out, and one dominant weight gives a value near one. This is the effective sample size the reliability-weighted overloads (`WeightType.Reliability`) use in their bias corrections.

```cs
using Numerics.Data.Statistics;

double[] weights = { 1, 1, 2 };
double neff = Statistics.EffectiveSampleSize(weights);   // 16 / 6 = 2.6667
```

Both sums run left to right from zero in index order, and the result is W * W / S. The bits therefore match a plain sequential loop in any language, and can differ in the last bits from a pairwise sum such as `numpy.sum` beyond eight weights. The weights are not rescaled, so weights whose squares underflow (for example ten weights of 1e-200) give NaN; divide such weights by their maximum first. Negative or non-finite weights throw `ArgumentOutOfRangeException`, all-zero weights throw `ArgumentException`, and an empty vector returns NaN.

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

**Methods:** `Push(value)` and `Push(values)` add data, `Clone()` returns an independent copy, `Combine(a, b)` (or `a + b`) merges two accumulators, and `Clear()` resets an instance to the empty state without reallocating. After `Clear()`, pushing the same values gives bitwise the same statistics as a new instance.

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

// Reuse the accumulator for the next batch
runningStats.Clear();
runningStats.Push(new double[] { 20.1, 22.4, 19.8 });
Console.WriteLine($"After Clear: Count = {runningStats.Count}, Mean = {runningStats.Mean:F2}");
```

## Box-Cox Transformation

The `BoxCox` class applies the Box-Cox power transformation [[6]](#6) to positive data, to make a skewed sample closer to normal:

```math
y^{(\lambda)} = \begin{cases} \dfrac{x^\lambda - 1}{\lambda} & \lambda \neq 0 \\ \ln x & \lambda = 0 \end{cases}, \qquad \frac{dy^{(\lambda)}}{dx} = x^{\lambda - 1}
```

The library uses the logarithm when $|\lambda| < 10^{-8}$ and limits $\lambda$ to $[-5, 5]$. Values $x \le 0$ or exponents outside that range give NaN.

| Method | Description |
|--------|-------------|
| `FitLambda(values, out lambda)` | Maximum likelihood exponent by Brent search on $[-5, 5]$; NaN when the fit fails |
| `TryFitLambda(values, out lambda)` | The same exponent, returning `false` instead of a silent NaN when the sample cannot support a fit (fewer than two values, a non-positive or non-finite value, all values equal) or the search fails |
| `Transform(value, lambda)` / `Transform(values, lambda)` | The transformation |
| `InverseTransform(value, lambda)` / `InverseTransform(values, lambda)` | The back-transformation $(\lambda y + 1)^{1/\lambda}$, or $e^y$ |
| `Derivative(value, lambda)` | $x^{\lambda - 1}$ (or $1/x$ on the logarithm branch), the Jacobian of the change of variable |
| `LogLikelihood(values, lambda)` | The normal log-likelihood of the transformed values at their maximum likelihood mean and variance, $-\frac{n}{2}\ln(2\pi\hat{\sigma}^2) - \frac{n}{2} + (\lambda - 1)\sum \ln x_i$, maximized by `FitLambda`; comparable with the log-likelihood of an untransformed normal fit |

```cs
using Numerics.Data.Statistics;

double[] flows = { 1, 3, 59, 1, 6, 4, 9, 13, 5, 84, 35, 8, 31, 34, 9, 1, 35, 66, 1, 65, 4, 68, 46 };

if (BoxCox.TryFitLambda(flows, out double lambda))
{
    double y = BoxCox.Transform(59.0, lambda);
    double back = BoxCox.InverseTransform(y, lambda);
    double dydx = BoxCox.Derivative(59.0, lambda);
    Console.WriteLine($"λ = {lambda:F4}, y = {y:F4}, back = {back:F2}, dy/dx = {dydx:F5}");
}
else
{
    Console.WriteLine("The sample cannot support a Box-Cox fit.");
}
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

<a id="5">[5]</a> Knight, W. R. (1966). A computer method for calculating Kendall's tau with ungrouped data. *Journal of the American Statistical Association*, 61(314), 436-439.

<a id="6">[6]</a> Box, G. E. P., & Cox, D. R. (1964). An analysis of transformations. *Journal of the Royal Statistical Society: Series B*, 26(2), 211-252.

<a id="7">[7]</a> Kish, L. (1965). *Survey Sampling*. Wiley, New York.

<a id="8">[8]</a> Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen. *ZAMM - Zeitschrift für Angewandte Mathematik und Mechanik*, 54(1), 39-51. https://doi.org/10.1002/zamm.19740540106

<a id="9">[9]</a> Chan, T. F., Golub, G. H., & LeVeque, R. J. (1983). Algorithms for computing the sample variance: Analysis and recommendations. *The American Statistician*, 37(3), 242-247. https://doi.org/10.1080/00031305.1983.10483115

<a id="10">[10]</a> Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. *Journal of the American Statistical Association*, 22(158), 209-212. https://doi.org/10.1080/01621459.1927.10502953

---

[← Previous: Time Series](../data/time-series.md) | [Back to Index](../index.md) | [Next: Goodness-of-Fit →](goodness-of-fit.md)
