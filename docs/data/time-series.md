# Time Series

The ***Numerics*** library provides classes for working with time series data, including data structures, operations, and downloading from hydrologic data sources.

## Overview

| Class | Description |
|-------|-------------|
| `TimeSeries` | Basic time series container |
| `PairedData` | X-Y paired data (e.g., rating curves) |
| `TimeSeriesDownload` | Download from USGS, Canadian, Australian sources |

---

## TimeSeries Class

### Creating Time Series

```cs
using Numerics.Data.TimeSeries;

// From arrays
DateTime[] times = {
    new DateTime(2020, 1, 1),
    new DateTime(2020, 1, 2),
    new DateTime(2020, 1, 3),
    new DateTime(2020, 1, 4),
    new DateTime(2020, 1, 5)
};
double[] values = { 100, 150, 120, 180, 140 };

var ts = new TimeSeries(times, values);

// With metadata
ts.Name = "Daily Streamflow";
ts.Units = "cfs";
ts.StationID = "12345678";
```

### Accessing Data

```cs
// Get values
DateTime[] allTimes = ts.Times;
double[] allValues = ts.Values;

// Single value
double val = ts[2];  // Third value
DateTime time = ts.Times[2];

// Length
int n = ts.Count;

// Time range
DateTime start = ts.StartTime;
DateTime end = ts.EndTime;
TimeSpan duration = ts.Duration;
```

### Basic Operations

```cs
// Statistics
double mean = ts.Mean();
double min = ts.Minimum();
double max = ts.Maximum();
double sum = ts.Sum();
double stdDev = ts.StandardDeviation();

// Missing data
int missingCount = ts.MissingCount;
double percentMissing = ts.PercentMissing;
bool hasMissing = ts.HasMissingValues;

// Fill missing values
ts.FillMissing(FillMethod.Linear);      // Linear interpolation
ts.FillMissing(FillMethod.Previous);    // Carry forward
ts.FillMissing(FillMethod.Next);        // Carry backward
ts.FillMissing(FillMethod.Mean);        // Replace with mean
ts.FillMissing(0.0);                    // Replace with constant
```

### Time Series Arithmetic

```cs
var ts1 = new TimeSeries(times, values1);
var ts2 = new TimeSeries(times, values2);

// Element-wise operations (requires matching times)
var sum = ts1 + ts2;
var diff = ts1 - ts2;
var product = ts1 * ts2;
var quotient = ts1 / ts2;

// Scalar operations
var scaled = ts1 * 2.0;
var offset = ts1 + 10.0;
```

### Resampling

```cs
// Resample to different frequency
var daily = ts.Resample(TimeSpan.FromDays(1), AggregationMethod.Mean);
var monthly = ts.ResampleMonthly(AggregationMethod.Sum);
var annual = ts.ResampleAnnual(AggregationMethod.Maximum);

// Aggregation methods
// - Mean, Sum, Minimum, Maximum
// - First, Last
// - Count
```

### Subsetting

```cs
// By time range
var subset = ts.GetRange(
    new DateTime(2020, 3, 1),
    new DateTime(2020, 6, 30)
);

// By water year
var wy2020 = ts.GetWaterYear(2020);  // Oct 2019 - Sep 2020

// By month
var januaryData = ts.GetMonth(1);

// By season
var summer = ts.GetSeason(Season.Summer);  // Jun-Aug
```

---

## Annual Maximum Series

Extract annual maxima for flood frequency analysis:

```cs
// Daily streamflow data
var dailyFlow = LoadDailyData();

// Extract annual maxima (water year basis)
var annualMax = dailyFlow.AnnualMaxima(WaterYearStart.October);

Console.WriteLine("Water Year   Annual Max (cfs)");
Console.WriteLine("----------   ----------------");

for (int i = 0; i < annualMax.Count; i++)
{
    int waterYear = annualMax.Times[i].Year;
    if (annualMax.Times[i].Month >= 10) waterYear++;
    
    Console.WriteLine($"    {waterYear}        {annualMax.Values[i]:N0}");
}

// Use for frequency analysis
double[] maxValues = annualMax.Values;
var gev = new GeneralizedExtremeValue();
gev.SetParameters(gev.ParametersFromLinearMoments(maxValues));
```

### Peak Over Threshold

```cs
double threshold = 10000;  // cfs

// Extract peaks above threshold
var peaks = dailyFlow.PeaksOverThreshold(threshold, 
    minimumSeparation: TimeSpan.FromDays(7));

Console.WriteLine($"Found {peaks.Count} peaks above {threshold} cfs");
```

---

## Data Download

### USGS WaterServices

Download streamflow and other parameters from USGS:

```cs
using Numerics.Data.TimeSeries;

// Download daily discharge
var ts = TimeSeriesDownload.FromUSGS(
    stationNumber: "09380000",           // Colorado River at Lees Ferry
    parameterCode: "00060",              // Discharge
    startDate: new DateTime(2020, 1, 1),
    endDate: new DateTime(2020, 12, 31),
    dataType: USGSDataType.DailyValues
);

Console.WriteLine($"Downloaded {ts.Count} values");
Console.WriteLine($"Station: {ts.StationID}");
Console.WriteLine($"Period: {ts.StartTime:d} to {ts.EndTime:d}");
```

**Common USGS Parameter Codes**:

| Code | Parameter |
|------|-----------|
| 00060 | Discharge (cfs) |
| 00065 | Gage height (ft) |
| 00010 | Temperature (°C) |
| 00045 | Precipitation (in) |
| 00400 | pH |

**Data Types**:
- `DailyValues` - Daily mean values
- `InstantaneousValues` - Unit values (15-min, hourly, etc.)
- `MonthlyStatistics` - Monthly statistics
- `AnnualStatistics` - Annual statistics

### Environment Canada

Download hydrometric data from Canadian stations:

```cs
var ts = TimeSeriesDownload.FromCanada(
    stationNumber: "05BB001",            // Bow River at Banff
    parameterType: CanadaParameterType.Discharge,
    startDate: new DateTime(2020, 1, 1),
    endDate: new DateTime(2020, 12, 31)
);
```

### Australian Bureau of Meteorology

Download from Australian KiWIS API:

```cs
var ts = TimeSeriesDownload.FromABOM(
    stationNumber: "410730",             // Cotter River at Gingera
    parameterType: ABOMParameterType.Discharge,
    startDate: new DateTime(2020, 1, 1),
    endDate: new DateTime(2020, 12, 31)
);
```

---

## PairedData Class

For X-Y relationships like rating curves:

```cs
// Stage-discharge rating curve
double[] stage = { 0, 1, 2, 3, 4, 5 };
double[] discharge = { 0, 100, 400, 900, 1600, 2500 };

var rating = new PairedData(stage, discharge);
rating.Name = "Rating Curve";
rating.XLabel = "Stage (ft)";
rating.YLabel = "Discharge (cfs)";

// Interpolate
double q = rating.Interpolate(2.5);
Console.WriteLine($"Q at stage 2.5 ft: {q:F0} cfs");

// Inverse interpolate (find stage for given Q)
double h = rating.InverseInterpolate(500);
Console.WriteLine($"Stage at Q=500 cfs: {h:F2} ft");
```

### Extending Rating Curves

```cs
// Fit power function: Q = a * (h - h0)^b
var (a, b, h0) = rating.FitPowerFunction();

Console.WriteLine($"Q = {a:F2} * (h - {h0:F2})^{b:F2}");

// Extrapolate beyond measured range
double highStage = 7.0;
double extrapolatedQ = a * Math.Pow(highStage - h0, b);
```

---

## Flow Duration Curves

```cs
var dailyFlow = LoadDailyData();

// Compute flow duration curve
var fdc = dailyFlow.FlowDurationCurve();

Console.WriteLine("Exceedance %   Flow (cfs)");
Console.WriteLine("-----------   ----------");

double[] exceedances = { 10, 25, 50, 75, 90, 95, 99 };
foreach (double pct in exceedances)
{
    double flow = fdc.Interpolate(pct / 100.0);
    Console.WriteLine($"    {pct,3:F0}%        {flow:N0}");
}

// Q50 (median flow)
double q50 = fdc.Interpolate(0.50);

// Q90 (low flow - exceeded 90% of time)
double q90 = fdc.Interpolate(0.90);
```

---

## Moving Window Statistics

```cs
var ts = LoadTimeSeries();

// Moving average
int windowSize = 7;  // 7-day moving average
var movingAvg = ts.MovingAverage(windowSize);

// Moving standard deviation
var movingStd = ts.MovingStandardDeviation(windowSize);

// Moving min/max
var movingMin = ts.MovingMinimum(windowSize);
var movingMax = ts.MovingMaximum(windowSize);
```

---

## Hydrologic Calculations

### Baseflow Separation

```cs
var dailyFlow = LoadDailyData();

// Digital filter baseflow separation
var baseflow = dailyFlow.BaseflowSeparation(
    filterParameter: 0.925,
    passes: 3
);

var quickflow = dailyFlow - baseflow;

double baseflowIndex = baseflow.Sum() / dailyFlow.Sum();
Console.WriteLine($"Baseflow Index: {baseflowIndex:P1}");
```

### Cumulative Volume

```cs
// Convert flow to volume
var dailyFlow = LoadDailyData();  // cfs

// Cumulative volume in acre-feet
// 1 cfs-day = 1.9835 acre-feet
var cumVolume = dailyFlow.CumulativeSum() * 1.9835;

Console.WriteLine($"Total annual volume: {cumVolume.Values.Last():N0} acre-ft");
```

### Unit Conversion

```cs
// cfs to cms
var flowCMS = dailyFlow * 0.028317;

// Inches to mm
var precipMM = precipInches * 25.4;
```

---

## Data Quality

### Gap Detection

```cs
var ts = LoadTimeSeries();

var gaps = ts.FindGaps(expectedInterval: TimeSpan.FromDays(1));

Console.WriteLine($"Found {gaps.Count} gaps:");
foreach (var gap in gaps)
{
    Console.WriteLine($"  {gap.Start:d} to {gap.End:d} ({gap.Duration.Days} days)");
}
```

### Outlier Detection

```cs
// Flag values outside ±3 standard deviations
var flags = ts.FlagOutliers(threshold: 3.0);

// Replace outliers with interpolated values
ts.ReplaceOutliers(threshold: 3.0, FillMethod.Linear);
```

### Consistency Checks

```cs
// Check for negative values (invalid for discharge)
bool hasNegative = ts.Values.Any(v => v < 0);

// Check for flat-lining (sensor malfunction)
int consecutiveSame = ts.MaxConsecutiveSameValues();
if (consecutiveSame > 10)
    Console.WriteLine("Warning: Possible sensor malfunction");

// Check for spikes
var spikes = ts.FindSpikes(threshold: 5.0);  // 5x median absolute deviation
```

---

## Export and Import

### CSV Export

```cs
var ts = LoadTimeSeries();

// Export to CSV
ts.ToCsv("output.csv");

// With custom format
ts.ToCsv("output.csv", 
    dateFormat: "yyyy-MM-dd",
    valueFormat: "F2",
    includeHeader: true);
```

### CSV Import

```cs
var ts = TimeSeries.FromCsv("input.csv",
    dateColumn: 0,
    valueColumn: 1,
    dateFormat: "MM/dd/yyyy",
    hasHeader: true);
```

---

## Example: Annual Flood Frequency Analysis

```cs
using Numerics.Data.TimeSeries;
using Numerics.Distributions;
using Numerics.Data.Statistics;

// Download historical daily flows
var daily = TimeSeriesDownload.FromUSGS(
    stationNumber: "09380000",
    parameterCode: "00060",
    startDate: new DateTime(1920, 1, 1),
    endDate: new DateTime(2023, 12, 31),
    dataType: USGSDataType.DailyValues
);

Console.WriteLine($"Downloaded {daily.Count} daily values");
Console.WriteLine($"Period: {daily.StartTime:yyyy} to {daily.EndTime:yyyy}");

// Extract annual maxima
var annualMax = daily.AnnualMaxima(WaterYearStart.October);
double[] maxValues = annualMax.Values;

Console.WriteLine($"\n{annualMax.Count} years of annual maxima");
Console.WriteLine($"Mean: {Statistics.Mean(maxValues):N0} cfs");
Console.WriteLine($"Max:  {Statistics.Maximum(maxValues):N0} cfs");

// Fit Log-Pearson Type III
var lp3 = new LogPearsonTypeIII();
lp3.SetParameters(lp3.ParametersFromLinearMoments(maxValues));

// Compute flood quantiles
Console.WriteLine("\nFlood Frequency Estimates:");
Console.WriteLine("Return Period    Discharge (cfs)");
Console.WriteLine("-------------    ---------------");

double[] returnPeriods = { 2, 5, 10, 25, 50, 100, 200, 500 };
foreach (double T in returnPeriods)
{
    double p = 1 - 1.0 / T;
    double q = lp3.InverseCDF(p);
    Console.WriteLine($"    {T,3:F0}-year        {q,12:N0}");
}
```

---

## References

<a id="ref1">[1]</a> Helsel, D. R., Hirsch, R. M., Ryberg, K. R., Archfield, S. A., & Gilroy, E. J. (2020). *Statistical Methods in Water Resources*. U.S. Geological Survey Techniques and Methods, Book 4, Chapter A3.

<a id="ref2">[2]</a> Eckhardt, K. (2005). How to construct recursive digital filters for baseflow separation. *Hydrological Processes*, 19(2), 507-515.

<a id="ref3">[3]</a> U.S. Geological Survey (2023). *USGS Water Services*. https://waterservices.usgs.gov/
