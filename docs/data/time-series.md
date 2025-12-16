# Time Series

[← Previous: Interpolation](interpolation.md) | [Back to Index](../index.md) | [Next: Random Generation →](../sampling/random-generation.md)

The ***Numerics*** library provides a comprehensive `TimeSeries` class for working with time-indexed data. This class supports regular and irregular time intervals, statistical operations, transformations, and analysis methods essential for hydrological and environmental data.

## Creating Time Series

### Empty Time Series

```cs
using Numerics.Data;

// Create empty time series
var ts = new TimeSeries();

// Create with time interval
var dailyData = new TimeSeries(TimeInterval.Daily);
var monthlyData = new TimeSeries(TimeInterval.Monthly);
```

### Time Series with Date Range

```cs
// Create time series with start and end dates
DateTime start = new DateTime(2020, 1, 1);
DateTime end = new DateTime(2020, 12, 31);

// With NaN values (placeholder)
var ts1 = new TimeSeries(TimeInterval.Daily, start, end);

// With fixed value
var ts2 = new TimeSeries(TimeInterval.Daily, start, end, fixedValue: 0.0);

Console.WriteLine($"Created time series with {ts1.Count} daily values");
```

### Time Series from Data

```cs
// Create from array of values
double[] dailyFlow = { 125.0, 130.0, 135.0, 132.0, 138.0 };
DateTime start = new DateTime(2024, 1, 1);

var ts = new TimeSeries(TimeInterval.Daily, start, dailyFlow);

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
    Irregular,      // No fixed interval
    OneMinute,
    FiveMinute,
    TenMinute,
    FifteenMinute,
    ThirtyMinute,
    OneHour,
    SixHour,
    TwelveHour,
    Daily,
    Weekly,
    Monthly,
    Annual
}

// Example usage
var hourlyData = new TimeSeries(TimeInterval.OneHour);
var yearlyData = new TimeSeries(TimeInterval.Annual);
```

## Accessing Data

### Indexing

```cs
var ts = new TimeSeries(TimeInterval.Daily, new DateTime(2024, 1, 1), 
                        new[] { 10.0, 15.0, 20.0, 25.0, 30.0 });

// Access by index
double value = ts[2].Value;                // 20.0
DateTime date = ts[2].Index;               // 2024-01-03

// Access by date
DateTime queryDate = new DateTime(2024, 1, 3);
var ordinate = ts.GetByIndex(queryDate);
Console.WriteLine($"Flow on {queryDate:yyyy-MM-dd}: {ordinate.Value:F1}");

// Properties
int count = ts.Count;
DateTime firstDate = ts.FirstIndex;
DateTime lastDate = ts.LastIndex;
double[] values = ts.Values.ToArray();
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
var ts = new TimeSeries(TimeInterval.Daily, new DateTime(2024, 1, 1),
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
var ts = new TimeSeries(TimeInterval.Daily, new DateTime(2024, 1, 1),
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

### Cumulative Sum

```cs
double[] dailyRainfall = { 0.5, 1.2, 0.8, 0.0, 2.1, 1.5 };
var rainfall = new TimeSeries(TimeInterval.Daily, new DateTime(2024, 1, 1), dailyRainfall);

// Compute cumulative rainfall
var cumulative = rainfall.CumulativeSum();

Console.WriteLine("Day | Daily | Cumulative");
for (int i = 0; i < rainfall.Count; i++)
{
    Console.WriteLine($"{i + 1,3} | {rainfall[i].Value,5:F1} | {cumulative[i].Value,10:F1}");
}
```

### Differencing

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
using Numerics.Data.Statistics;

var ts = new TimeSeries(TimeInterval.Daily, new DateTime(2024, 1, 1),
                        new[] { 125.0, 130.0, 135.0, 132.0, 138.0, 145.0 });

// Compute statistics
double mean = Statistics.Mean(ts.Values.ToArray());
double std = Statistics.StandardDeviation(ts.Values.ToArray());
double min = ts.Values.Min();
double max = ts.Values.Max();

Console.WriteLine($"Mean: {mean:F1}");
Console.WriteLine($"Std Dev: {std:F1}");
Console.WriteLine($"Range: [{min:F1}, {max:F1}]");

// Percentiles
double[] values = ts.Values.ToArray();
double p25 = Statistics.Percentile(values, 25);
double p50 = Statistics.Percentile(values, 50);
double p75 = Statistics.Percentile(values, 75);

Console.WriteLine($"25th percentile: {p25:F1}");
Console.WriteLine($"Median: {p50:F1}");
Console.WriteLine($"75th percentile: {p75:F1}");
```

### Moving Average

```cs
// Compute moving average
int window = 3;
var movingAvg = new TimeSeries(ts.TimeInterval);

for (int i = window - 1; i < ts.Count; i++)
{
    double sum = 0;
    for (int j = 0; j < window; j++)
    {
        sum += ts[i - j].Value;
    }
    double avg = sum / window;
    movingAvg.Add(new SeriesOrdinate<DateTime, double>(ts[i].Index, avg));
}

Console.WriteLine("Original | 3-day Moving Average");
for (int i = 0; i < ts.Count; i++)
{
    string ma = i < movingAvg.Count ? $"{movingAvg[i].Value:F1}" : "N/A";
    Console.WriteLine($"{ts[i].Value,8:F1} | {ma,22}");
}
```

## Sorting and Filtering

### Sorting

```cs
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
// Monthly flow data
double[] monthlyFlow = { 125, 135, 180, 220, 250, 280, 260, 230, 190, 150, 130, 120 };
var flowData = new TimeSeries(TimeInterval.Monthly, new DateTime(2024, 1, 1), monthlyFlow);

Console.WriteLine("Monthly Streamflow Analysis");
Console.WriteLine("=" + new string('=', 50));

// Find annual peak
double peakFlow = flowData.Values.Max();
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
var filled = new TimeSeries(TimeInterval.Daily, dates.Min(), dates.Max());
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
double[] annualPeaks = { 1200, 1250, 1180, 1300, 1320, 1280, 1350, 1400, 1380, 1450 };
var years = Enumerable.Range(2015, 10).Select(y => new DateTime(y, 1, 1)).ToArray();

var peakSeries = new TimeSeries(TimeInterval.Annual);
for (int i = 0; i < years.Length; i++)
{
    peakSeries.Add(new SeriesOrdinate<DateTime, double>(years[i], annualPeaks[i]));
}

Console.WriteLine("Annual Peak Flow Trend Analysis");
Console.WriteLine("=" + new string('=', 50));

// Linear regression for trend
double[] x = Enumerable.Range(0, peakSeries.Count).Select(i => (double)i).ToArray();
double[] y = peakSeries.Values.ToArray();

double xMean = x.Average();
double yMean = y.Average();

double slope = x.Zip(y, (xi, yi) => (xi - xMean) * (yi - yMean)).Sum() /
               x.Sum(xi => Math.Pow(xi - xMean, 2));
double intercept = yMean - slope * xMean;

Console.WriteLine($"Trend: {slope:F1} cfs/year");
Console.WriteLine($"Direction: {(slope > 0 ? "Increasing" : "Decreasing")}");

// Mann-Kendall test for significance
double mkStat = HypothesisTests.MannKendallTest(y);
Console.WriteLine($"Mann-Kendall statistic: {mkStat:F2}");

if (Math.Abs(mkStat) > 1.96)
    Console.WriteLine("Trend is statistically significant (p < 0.05)");
else
    Console.WriteLine("Trend is not statistically significant");
```

### Example 4: Seasonal Analysis

```cs
// Multi-year daily data
int years = 3;
int daysPerYear = 365;
var random = new Random(123);

var dailyTemp = new TimeSeries(TimeInterval.Daily, new DateTime(2022, 1, 1));

// Generate seasonal temperature pattern
for (int day = 0; day < years * daysPerYear; day++)
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

// Compute monthly averages
var monthlyAvg = new Dictionary<int, List<double>>();
for (int m = 1; m <= 12; m++)
    monthlyAvg[m] = new List<double>();

foreach (var ord in dailyTemp)
{
    monthlyAvg[ord.Index.Month].Add(ord.Value);
}

Console.WriteLine("\nMonth | Avg Temp (°C)");
Console.WriteLine("------|-------------");
for (int m = 1; m <= 12; m++)
{
    double avg = monthlyAvg[m].Average();
    Console.WriteLine($"{m,5} | {avg,13:F1}");
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
| Moving average | Custom loop | Smoothing |
| Sort | `SortByTime()` | Ensure chronological order |

---

[← Previous: Interpolation](interpolation.md) | [Back to Index](../index.md) | [Next: Random Generation →](../sampling/random-generation.md)
