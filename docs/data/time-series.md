# Time Series

[← Previous: Linear Regression](regression.md) | [Back to Index](../index.md) | [Next: Descriptive Statistics →](../statistics/descriptive.md)

The ***Numerics*** library provides a comprehensive `TimeSeries` class for working with time-indexed data. This class supports regular and irregular time intervals, statistical operations, transformations, and analysis methods essential for hydrological and environmental data.

## Creating Time Series

### Empty Time Series

```cs
using Numerics.Data;

// Create empty time series
var ts = new TimeSeries();

// Create with time interval
var dailyData = new TimeSeries(TimeInterval.OneDay);
var monthlyData = new TimeSeries(TimeInterval.OneMonth);
```

### Time Series with Date Range

```cs
// Create time series with start and end dates
DateTime start = new DateTime(2020, 1, 1);
DateTime end = new DateTime(2020, 12, 31);

// Create empty time series (values default to NaN)
var ts1 = new TimeSeries(TimeInterval.OneDay, start, end);

// With fixed value
var ts2 = new TimeSeries(TimeInterval.OneDay, start, end, fixedValue: 0.0);

Console.WriteLine($"Created time series with {ts1.Count} daily values");
```

### Time Series from Data

```cs
// Create from array of values
double[] dailyFlow = { 125.0, 130.0, 135.0, 132.0, 138.0 };
DateTime start = new DateTime(2024, 1, 1);

var ts = new TimeSeries(TimeInterval.OneDay, start, dailyFlow);

Console.WriteLine("Daily Flow Data:");
for (int i = 0; i < ts.Count; i++)
{
    Console.WriteLine($"{ts[i].Index:yyyy-MM-dd}: {ts[i].Value:F1} cfs");
}
```

## Time Intervals

Supported time intervals:

```cs
public enum TimeInterval
{
    OneMinute,
    FiveMinute,
    FifteenMinute,
    ThirtyMinute,
    OneHour,
    SixHour,
    TwelveHour,
    OneDay,
    SevenDay,
    OneMonth,
    OneQuarter,
    OneYear,
    Irregular       // No fixed interval
}

// Example usage
var hourlyData = new TimeSeries(TimeInterval.OneHour);
var yearlyData = new TimeSeries(TimeInterval.OneYear);
```

## Accessing Data

### Indexing

```cs
using System.Linq;
using Numerics.Data;

var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 1, 1),
                        new[] { 10.0, 15.0, 20.0, 25.0, 30.0 });

// Access by index
double value = ts[2].Value;                // 20.0
DateTime date = ts[2].Index;               // 2024-01-03

// Access by date (use LINQ to find ordinate by date)
DateTime queryDate = new DateTime(2024, 1, 3);
var ordinate = ts.First(o => o.Index == queryDate);
Console.WriteLine($"Flow on {queryDate:yyyy-MM-dd}: {ordinate.Value:F1}");

// Properties
int count = ts.Count;
DateTime firstDate = ts.StartDate;
DateTime lastDate = ts.EndDate;
double[] values = ts.ValuesToArray();
```

### Missing Values

```cs
// Check for missing values
bool hasMissing = ts.HasMissingValues;
int numMissing = ts.NumberOfMissingValues();

Console.WriteLine($"Missing values: {numMissing} out of {ts.Count}");

// Identify missing
for (int i = 0; i < ts.Count; i++)
{
    if (double.IsNaN(ts[i].Value))
    {
        Console.WriteLine($"Missing value at {ts[i].Index:yyyy-MM-dd}");
    }
}
```

## Mathematical Operations

### Basic Arithmetic

```cs
var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 1, 1),
                        new[] { 10.0, 15.0, 20.0, 25.0, 30.0 });

// Add constant to all values
ts.Add(5.0);                    // Now: 15, 20, 25, 30, 35

// Subtract constant
ts.Subtract(3.0);               // Now: 12, 17, 22, 27, 32

// Multiply by constant
ts.Multiply(2.0);               // Now: 24, 34, 44, 54, 64

// Divide by constant
ts.Divide(2.0);                 // Back to: 12, 17, 22, 27, 32

Console.WriteLine("Transformed values:");
foreach (var ord in ts)
{
    Console.WriteLine($"{ord.Index:yyyy-MM-dd}: {ord.Value:F1}");
}
```

### Operations on Subsets

```cs
// Apply operations to specific indexes
var indexes = new List<int> { 0, 2, 4 };  // First, third, fifth values

ts.Add(10.0, indexes);          // Add 10 to selected values only
ts.Multiply(1.5, indexes);      // Multiply selected values by 1.5
```

### Transformations

```cs
var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 1, 1),
                        new[] { 10.0, 15.0, 20.0, 25.0, 30.0 });

// Absolute value
ts.AbsoluteValue();

// Exponentiation
ts.Exponentiate(2.0);           // Square all values

// Logarithm
ts.LogTransform(baseValue: 10); // Log₁₀ transform

// Standardize (z-score)
ts.Standardize();               // (x - μ) / σ

// Inverse
ts.Inverse();                   // 1 / x
```

## Time Series Analysis

### Autocorrelation

The **autocorrelation function** (ACF) measures the correlation of a time series with a lagged copy of itself. At lag $k$, the sample autocorrelation is:

```math
\hat{\rho}(k) = \frac{\sum_{t=1}^{n-k}(x_t - \bar{x})(x_{t+k} - \bar{x})}{\sum_{t=1}^{n}(x_t - \bar{x})^2}
```

where $\bar{x}$ is the sample mean and $n$ is the series length. By definition, $\hat{\rho}(0) = 1$.

Autocorrelation is central to time series analysis: significant autocorrelation at lag $k$ indicates that values $k$ time steps apart are linearly related. For an independent series, $\hat{\rho}(k) \approx 0$ for all $k > 0$, and the approximate 95% confidence bounds are $\pm 1.96 / \sqrt{n}$.

The ***Numerics*** library uses the ACF internally in several statistical tests. The **Ljung-Box test** (`SummaryHypothesisTest()`) checks whether a group of autocorrelations are jointly significant:

```math
Q = n(n+2) \sum_{k=1}^{h} \frac{\hat{\rho}(k)^2}{n - k}
```

Under the null hypothesis of independence, $Q \sim \chi^2(h)$.

### Cumulative Sum

```cs
double[] dailyRainfall = { 0.5, 1.2, 0.8, 0.0, 2.1, 1.5 };
var rainfall = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 1, 1), dailyRainfall);

// Compute cumulative rainfall
var cumulative = rainfall.CumulativeSum();

Console.WriteLine("Day | Daily | Cumulative");
for (int i = 0; i < rainfall.Count; i++)
{
    Console.WriteLine($"{i + 1,3} | {rainfall[i].Value,5:F1} | {cumulative[i].Value,10:F1}");
}
```

### Differencing

The **difference operator** $\nabla$ removes trends from a time series. The first difference at lag $d$ is:

```math
\nabla_d x_t = x_t - x_{t-d}
```

First differencing ($d = 1$) removes a linear trend. Applying the operator twice ($\nabla^2 x_t = \nabla(\nabla x_t)$) removes a quadratic trend. Seasonal differencing uses a lag equal to the seasonal period — for example, $d = 12$ for monthly data with an annual cycle removes the seasonal component directly.

```cs
// First difference (change from previous)
var diff1 = ts.Difference(lag: 1, differences: 1);

// Second difference
var diff2 = ts.Difference(lag: 1, differences: 2);

// Seasonal difference (e.g., monthly lag for annual pattern)
var seasonalDiff = ts.Difference(lag: 12, differences: 1);

Console.WriteLine("Original | First Diff | Second Diff");
for (int i = 0; i < Math.Min(5, ts.Count); i++)
{
    double orig = i < ts.Count ? ts[i].Value : double.NaN;
    double d1 = i < diff1.Count ? diff1[i].Value : double.NaN;
    double d2 = i < diff2.Count ? diff2[i].Value : double.NaN;
    
    Console.WriteLine($"{orig,8:F2} | {d1,10:F2} | {d2,11:F2}");
}
```

## Statistical Methods

### Descriptive Statistics

```cs
using System.Linq;
using Numerics.Data;
using Numerics.Data.Statistics;

var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 1, 1),
                        new[] { 125.0, 130.0, 135.0, 132.0, 138.0, 145.0 });

// Compute statistics
double mean = Statistics.Mean(ts.ValuesToArray());
double std = Statistics.StandardDeviation(ts.ValuesToArray());
double min = ts.ValuesToArray().Min();
double max = ts.ValuesToArray().Max();

Console.WriteLine($"Mean: {mean:F1}");
Console.WriteLine($"Std Dev: {std:F1}");
Console.WriteLine($"Range: [{min:F1}, {max:F1}]");

// Percentiles
double[] values = ts.ValuesToArray();
double p25 = Statistics.Percentile(values, 0.25);
double p50 = Statistics.Percentile(values, 0.50);
double p75 = Statistics.Percentile(values, 0.75);

Console.WriteLine($"25th percentile: {p25:F1}");
Console.WriteLine($"Median: {p50:F1}");
Console.WriteLine($"75th percentile: {p75:F1}");
```

### Moving Average

A **simple moving average** (SMA) of period $m$ smooths the series by replacing each value with the average of the preceding $m$ observations:

```math
\text{SMA}_t = \frac{1}{m} \sum_{j=0}^{m-1} x_{t-j}
```

The moving average acts as a low-pass filter: it attenuates fluctuations with period shorter than $m$ while preserving longer-term trends. The output series has $n - m + 1$ values because the first $m - 1$ observations lack a full window.

The ***Numerics*** library provides built-in `MovingAverage` and `MovingSum` methods that use a sliding window for $O(n)$ efficiency:

```cs
var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2024, 1, 1),
                        new[] { 125.0, 130.0, 135.0, 132.0, 138.0, 145.0 });

// Built-in moving average
var movingAvg = ts.MovingAverage(period: 3);

// Built-in moving sum
var movingSum = ts.MovingSum(period: 3);

Console.WriteLine("Original | 3-day MA | 3-day Sum");
for (int i = 0; i < movingAvg.Count; i++)
{
    Console.WriteLine($"{movingAvg[i].Index:yyyy-MM-dd}: {movingAvg[i].Value,8:F1} | {movingSum[i].Value,8:F1}");
}
```

## Block Series

The ***Numerics*** library can aggregate a time series into annual, monthly, or quarterly blocks using a specified function (minimum, maximum, average, or sum). This is essential for extracting annual maxima for flood frequency analysis, computing monthly means, or aggregating sub-daily data.

### Calendar and Water Year Series

```cs
// Annual maximum flow (calendar year: Jan–Dec)
var annualMax = dailyFlow.CalendarYearSeries(BlockFunctionType.Maximum);

// Water year maximum (Oct–Sep, standard in US hydrology)
var waterYearMax = dailyFlow.CustomYearSeries(startMonth: 10, BlockFunctionType.Maximum);

// Custom season: June–August summer average
var summerAvg = dailyFlow.CustomYearSeries(
    startMonth: 6, endMonth: 8, BlockFunctionType.Average);
```

### Monthly and Quarterly Series

```cs
// Monthly average flow
var monthlyAvg = dailyFlow.MonthlySeries(BlockFunctionType.Average);

// Quarterly maximum
var quarterlyMax = dailyFlow.QuarterlySeries(BlockFunctionType.Maximum);
```

### Smoothing Before Aggregation

Block methods support optional smoothing before aggregation. For example, to find the annual maximum 7-day average flow (a common low-flow statistic):

```cs
// Annual maximum of 7-day moving average
var annualMax7Day = dailyFlow.CalendarYearSeries(
    BlockFunctionType.Maximum,
    SmoothingFunctionType.MovingAverage,
    period: 7);
```

### Peaks Over Threshold

Extract independent peaks that exceed a threshold, with a minimum separation between events:

```cs
// Find flood peaks exceeding 500 cfs, at least 7 days apart
var peaks = dailyFlow.PeaksOverThresholdSeries(
    threshold: 500.0, minStepsBetweenEvents: 7);
```

## Time Interval Conversion

Convert between time intervals by aggregation (downsampling) or interpolation (upsampling):

```cs
// Convert daily to monthly (average)
var monthly = dailyFlow.ConvertTimeInterval(TimeInterval.OneMonth, average: true);

// Convert daily to monthly (sum — e.g., for precipitation)
var monthlySum = dailyPrecip.ConvertTimeInterval(TimeInterval.OneMonth, average: false);
```

When downsampling, `average: true` computes the block mean and `average: false` computes the block sum. When upsampling, `average: true` uses linear interpolation and `average: false` disaggregates proportionally.

## Missing Data Interpolation

The `InterpolateMissingData` method fills gaps using linear interpolation in date-space, but only when the gap is smaller than a specified maximum:

```cs
// Fill gaps of up to 3 missing values by linear interpolation
ts.InterpolateMissingData(maxNumberOfMissing: 3);
```

This prevents unreliable interpolation across long gaps while filling short data dropouts.

## Resampling Methods

### Block Bootstrap

The **block bootstrap** preserves temporal dependence by resampling contiguous blocks rather than individual observations. Given a block size $b$, the method randomly selects blocks of $b$ consecutive values (with replacement) and concatenates them:

```cs
// Generate a 1000-step synthetic series preserving 30-day temporal structure
var resampled = ts.ResampleWithBlockBootstrap(
    timeSteps: 1000, blockSize: 30, seed: 42);
```

Unlike the standard bootstrap (which destroys autocorrelation), the block bootstrap retains the short-range dependence structure within each block.

### k-Nearest Neighbors

The **k-nearest neighbors** (k-NN) resampling method generates synthetic time series that preserve the multivariate dependence structure. At each step, it finds the $k$ nearest neighbors of the current state in the historical record (using Euclidean distance on standardized values) and randomly selects one as the next value:

```cs
// Generate synthetic series using 5-nearest neighbors
var synthetic = ts.ResampleWithKNN(
    timeSteps: 500, k: 5, seed: 42);
```

## Sorting and Filtering

### Sorting

```cs
using System.ComponentModel;

// Sort by time (ascending or descending)
ts.SortByTime(ListSortDirection.Ascending);

// Sort by value
ts.SortByValue(ListSortDirection.Descending);  // Largest first
```

### Filtering by Date Range

```cs
DateTime filterStart = new DateTime(2024, 6, 1);
DateTime filterEnd = new DateTime(2024, 8, 31);

var filtered = new TimeSeries(ts.TimeInterval);

foreach (var ordinate in ts)
{
    if (ordinate.Index >= filterStart && ordinate.Index <= filterEnd)
    {
        filtered.Add(ordinate);
    }
}

Console.WriteLine($"Filtered to {filtered.Count} values in date range");
```

## Practical Examples

### Example 1: Annual Peak Flow Analysis

```cs
using System.Linq;
using Numerics.Data;

// Monthly flow data
double[] monthlyFlow = { 125, 135, 180, 220, 250, 280, 260, 230, 190, 150, 130, 120 };
var flowData = new TimeSeries(TimeInterval.OneMonth, new DateTime(2024, 1, 1), monthlyFlow);

Console.WriteLine("Monthly Streamflow Analysis");
Console.WriteLine("=" + new string('=', 50));

// Find annual peak
double peakFlow = flowData.ValuesToArray().Max();
int peakMonth = Array.IndexOf(monthlyFlow, peakFlow) + 1;

Console.WriteLine($"Peak flow: {peakFlow:F0} cfs");
Console.WriteLine($"Peak month: Month {peakMonth}");

// Compute seasonal statistics
var spring = monthlyFlow.Skip(2).Take(3).ToArray();  // MAM
var summer = monthlyFlow.Skip(5).Take(3).ToArray();  // JJA

Console.WriteLine($"\nSeasonal Means:");
Console.WriteLine($"  Spring (MAM): {spring.Average():F0} cfs");
Console.WriteLine($"  Summer (JJA): {summer.Average():F0} cfs");
```

### Example 2: Filling Missing Values

```cs
using System.Linq;
using Numerics.Data;

// Time series with gaps
var dates = new[] { 
    new DateTime(2024, 1, 1),
    new DateTime(2024, 1, 2),
    new DateTime(2024, 1, 3),
    new DateTime(2024, 1, 6),
    new DateTime(2024, 1, 7)
};
var values = new[] { 10.0, 12.0, 11.0, 15.0, 14.0 };

var ts = new TimeSeries(TimeInterval.Irregular);
for (int i = 0; i < dates.Length; i++)
{
    ts.Add(new SeriesOrdinate<DateTime, double>(dates[i], values[i]));
}

Console.WriteLine("Original data has gaps:");
for (int i = 0; i < ts.Count - 1; i++)
{
    TimeSpan gap = ts[i + 1].Index - ts[i].Index;
    if (gap.TotalDays > 1)
    {
        Console.WriteLine($"  Gap of {gap.TotalDays:F0} days after {ts[i].Index:yyyy-MM-dd}");
    }
}

// Fill gaps with linear interpolation
var filled = new TimeSeries(TimeInterval.OneDay, dates.Min(), dates.Max());
foreach (var ord in filled)
{
    // Find surrounding values
    var before = ts.Where(o => o.Index <= ord.Index).OrderBy(o => o.Index).LastOrDefault();
    var after = ts.Where(o => o.Index >= ord.Index).OrderBy(o => o.Index).FirstOrDefault();
    
    if (before != null && after != null && before.Index != after.Index)
    {
        // Linear interpolation
        double t = (ord.Index - before.Index).TotalDays / (after.Index - before.Index).TotalDays;
        ord.Value = before.Value + t * (after.Value - before.Value);
    }
    else if (ts.Any(o => o.Index == ord.Index))
    {
        ord.Value = ts.First(o => o.Index == ord.Index).Value;
    }
}

Console.WriteLine($"\nFilled series: {filled.Count} continuous daily values");
```

### Example 3: Trend Detection

```cs
using System.Linq;
using Numerics.Data;
using Numerics.Data.Statistics;

double[] annualPeaks = { 1200, 1250, 1180, 1300, 1320, 1280, 1350, 1400, 1380, 1450 };
var years = Enumerable.Range(2015, 10).Select(y => new DateTime(y, 1, 1)).ToArray();

var peakSeries = new TimeSeries(TimeInterval.OneYear);
for (int i = 0; i < years.Length; i++)
{
    peakSeries.Add(new SeriesOrdinate<DateTime, double>(years[i], annualPeaks[i]));
}

Console.WriteLine("Annual Peak Flow Trend Analysis");
Console.WriteLine("=" + new string('=', 50));

// Trend analysis using built-in hypothesis tests
double[] indices = Enumerable.Range(0, peakSeries.Count).Select(i => (double)i).ToArray();
double[] y = peakSeries.ValuesToArray();

// Parametric: Linear regression t-test (returns p-value)
double linearPValue = HypothesisTests.LinearTrendTest(indices, y);
Console.WriteLine($"Linear trend test p-value: {linearPValue:F4}");

// Non-parametric: Mann-Kendall test (returns p-value)
double mkPValue = HypothesisTests.MannKendallTest(y);
Console.WriteLine($"Mann-Kendall p-value: {mkPValue:F4}");

if (linearPValue < 0.05 || mkPValue < 0.05)
    Console.WriteLine("Trend is statistically significant (p < 0.05)");
else
    Console.WriteLine("Trend is not statistically significant");
```

### Example 4: Seasonal Decomposition

A time series can be decomposed into trend ($T_t$), seasonal ($S_t$), and residual ($R_t$) components. The two standard models are:

- **Additive**: $x_t = T_t + S_t + R_t$ — when seasonal fluctuations are roughly constant in magnitude
- **Multiplicative**: $x_t = T_t \cdot S_t \cdot R_t$ — when seasonal fluctuations scale with the level

The `SeasonalDecompose` method performs classical additive decomposition using a moving average for the trend and FFT-based extraction for the seasonal component:

```cs
using System.Linq;
using Numerics.Data;

// Monthly temperature data (5 years)
int nYears = 5;
int period = 12;
var random = new Random(123);
var monthlyTemp = new TimeSeries(TimeInterval.OneMonth);

for (int i = 0; i < nYears * period; i++)
{
    double trend = 15.0 + 0.05 * i;   // Slight warming trend
    double seasonal = 10.0 * Math.Sin(2 * Math.PI * i / period);
    double noise = (random.NextDouble() - 0.5) * 2;
    monthlyTemp.Add(new SeriesOrdinate<DateTime, double>(
        new DateTime(2020, 1, 1).AddMonths(i), trend + seasonal + noise));
}

// Decompose into trend, seasonal, and residual components
var (trend, seasonal, residual) = monthlyTemp.SeasonalDecompose(period);

Console.WriteLine("Seasonal Decomposition:");
Console.WriteLine($"  Trend range: {trend.MinValue():F1} to {trend.MaxValue():F1}");
Console.WriteLine($"  Seasonal amplitude: {seasonal.Max() - seasonal.Min():F1}");
Console.WriteLine($"  Residual count: {residual.Count}");
```

### Example 5: Seasonal Analysis

```cs
using System.Linq;
using Numerics.Data;

// Multi-year daily data
int nYears = 3;
int daysPerYear = 365;
var random = new Random(123);

var dailyTemp = new TimeSeries(TimeInterval.OneDay);

// Generate seasonal temperature pattern
for (int day = 0; day < nYears * daysPerYear; day++)
{
    // Sinusoidal pattern + noise
    double seasonalTemp = 15 + 10 * Math.Sin(2 * Math.PI * day / 365.0);
    double noise = (random.NextDouble() - 0.5) * 4;
    dailyTemp.Add(new SeriesOrdinate<DateTime, double>(
        new DateTime(2022, 1, 1).AddDays(day),
        seasonalTemp + noise
    ));
}

Console.WriteLine("Seasonal Temperature Analysis");
Console.WriteLine("=" + new string('=', 50));

// Compute monthly averages using built-in MonthlySeries
var monthlyAvgSeries = dailyTemp.MonthlySeries(BlockFunctionType.Average);

Console.WriteLine("\nDate       | Avg Temp (°C)");
Console.WriteLine("-----------|-------------");
foreach (var ord in monthlyAvgSeries)
{
    Console.WriteLine($"{ord.Index:yyyy-MM} | {ord.Value,13:F1}");
}
```

## Best Practices

1. **Choose appropriate time interval** - Use regular intervals when possible
2. **Handle missing values** - Check for and appropriately handle NaN values
3. **Validate dates** - Ensure dates are in correct order
4. **Consider time zones** - Be aware of daylight saving time issues
5. **Document transformations** - Keep track of applied operations
6. **Use appropriate statistics** - Account for autocorrelation in time series data

## Common Operations Summary

| Operation | Method | Use Case |
|-----------|--------|----------|
| Add constant | `Add(value)` | Baseline adjustment |
| Transform | `LogTransform()` | Variance stabilization |
| Standardize | `Standardize()` | Compare different scales |
| Cumulative | `CumulativeSum()` | Total accumulation |
| Difference | `Difference(lag)` | Remove trends |
| Moving average | `MovingAverage(period)` | Smoothing |
| Moving sum | `MovingSum(period)` | Accumulation over window |
| Annual block | `CalendarYearSeries()` | Annual statistics |
| Monthly block | `MonthlySeries()` | Monthly aggregation |
| Convert interval | `ConvertTimeInterval()` | Up/downsampling |
| Peaks over threshold | `PeaksOverThresholdSeries()` | Event extraction |
| Block bootstrap | `ResampleWithBlockBootstrap()` | Synthetic generation |
| k-NN resampling | `ResampleWithKNN()` | Synthetic generation |
| Sort | `SortByTime()` | Ensure chronological order |

---

## References

<a id="1">[1]</a> Box, G. E. P., Jenkins, G. M., Reinsel, G. C., & Ljung, G. M. (2015). *Time Series Analysis: Forecasting and Control* (5th ed.). Wiley.

<a id="2">[2]</a> Kundzewicz, Z. W. & Robson, A. J. (2004). Change detection in hydrological records — a review of the methodology. *Hydrological Sciences Journal*, 49(1), 7–19.

---

[← Previous: Linear Regression](regression.md) | [Back to Index](../index.md) | [Next: Descriptive Statistics →](../statistics/descriptive.md)
