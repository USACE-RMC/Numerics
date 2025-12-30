# Uncertainty Analysis

[← Previous: Parameter Estimation](parameter-estimation.md) | [Back to Index](../index.md) | [Next: Copulas →](copulas.md)

Uncertainty analysis quantifies the confidence in estimated distribution parameters and derived quantities (like design floods or failure probabilities). The ***Numerics*** library provides comprehensive bootstrap resampling methods for assessing parameter and quantile uncertainty, which are essential for risk-informed decision making.

## Bootstrap Analysis Overview

Bootstrap resampling [[1]](#1) is a nonparametric method for estimating sampling distributions and confidence intervals. It works by:

1. Resampling the original data with replacement
2. Fitting the distribution to each bootstrap sample
3. Computing statistics from the ensemble of fitted distributions
4. Using percentiles of the bootstrap distribution for confidence intervals

The bootstrap is particularly valuable when:
- Analytical confidence intervals are unavailable
- Sample sizes are small to moderate
- Distribution of estimators is unknown or complex
- Dealing with extreme value distributions

## Creating a Bootstrap Analysis

```cs
using Numerics.Distributions;

// Original data
double[] annualPeakFlows = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

// Fit parent distribution
var gev = new GeneralizedExtremeValue();
gev.Estimate(annualPeakFlows, ParameterEstimationMethod.MethodOfLinearMoments);

// Create bootstrap analysis
var bootstrap = new BootstrapAnalysis(
    distribution: gev,
    estimationMethod: ParameterEstimationMethod.MethodOfLinearMoments,
    sampleSize: annualPeakFlows.Length,
    replications: 10000,
    seed: 12345);

Console.WriteLine("Bootstrap analysis configured:");
Console.WriteLine($"  Sample size: {bootstrap.SampleSize}");
Console.WriteLine($"  Replications: {bootstrap.Replications}");
Console.WriteLine($"  Estimation method: {bootstrap.EstimationMethod}");
```

### Constructor Parameters

- **distribution**: The fitted parent distribution (must implement `IBootstrappable`)
- **estimationMethod**: Parameter estimation method for each bootstrap sample
- **sampleSize**: Size of each bootstrap sample (typically original sample size)
- **replications**: Number of bootstrap replications (recommended: 10,000-50,000)
- **seed**: Random seed for reproducibility (default: 12345)

## Complete Uncertainty Analysis

The `Estimate()` method performs a comprehensive uncertainty analysis:

```cs
using Numerics.Distributions;

double[] data = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

// Fit parent distribution
var gev = new GeneralizedExtremeValue();
gev.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);

// Configure bootstrap
var bootstrap = new BootstrapAnalysis(gev, 
    ParameterEstimationMethod.MethodOfLinearMoments,
    sampleSize: data.Length,
    replications: 10000);

// Define probabilities of interest
var probabilities = new double[] { 0.5, 0.9, 0.95, 0.98, 0.99, 0.998 };

// Perform complete uncertainty analysis
var results = bootstrap.Estimate(
    probabilities: probabilities,
    alpha: 0.1,  // 90% confidence intervals (1-α)
    distributions: null,  // Will generate new bootstrap samples
    recordParameterSets: true);

Console.WriteLine("Uncertainty Analysis Complete:");
Console.WriteLine($"  Parent distribution: {results.ParentDistribution.DisplayName}");
Console.WriteLine($"  Confidence level: {(1 - 0.1) * 100:F0}%");
```

### UncertaintyAnalysisResults Structure

The results object contains:

```cs
// Parent (point estimate) distribution
UnivariateDistributionBase ParentDistribution

// Mode curve: quantiles from parent distribution
double[,] ModeCurve  // [probability, quantile]

// Mean curve: expected quantiles from bootstrap ensemble
double[,] MeanCurve  // [probability, quantile]

// Confidence intervals for quantiles
double[,] ConfidenceIntervals  // [probability, lower, upper]

// Bootstrap parameter sets (if recorded)
ParameterSet[] ParameterSets

// Goodness-of-fit metrics
double AIC, BIC, DIC, RMSE, ERL
```

## Accessing Uncertainty Results

### Mode and Mean Curves

```cs
var results = bootstrap.Estimate(probabilities, alpha: 0.1);

Console.WriteLine("Quantile Estimates:");
Console.WriteLine("AEP   | T (years) | Mode     | Mean     | 90% CI");
Console.WriteLine("--------------------------------------------------------");

for (int i = 0; i < probabilities.Length; i++)
{
    double prob = probabilities[i];
    double aep = 1 - prob;
    double T = 1.0 / aep;
    
    double mode = results.ModeCurve[i, 1];     // Point estimate
    double mean = results.MeanCurve[i, 1];     // Expected value
    double lower = results.ConfidenceIntervals[i, 1];  // Lower bound
    double upper = results.ConfidenceIntervals[i, 2];  // Upper bound
    
    Console.WriteLine($"{aep:F3} | {T,9:F1} | {mode,8:F0} | {mean,8:F0} | [{lower,6:F0}, {upper,6:F0}]");
}
```

### Visualizing Uncertainty

```cs
// Create plotting data for frequency curve with uncertainty bounds
var plotProbs = Enumerable.Range(1, 99).Select(i => i / 100.0).ToArray();

var results = bootstrap.Estimate(plotProbs, alpha: 0.1);

// Export for plotting
using (var writer = new System.IO.StreamWriter("frequency_curve.csv"))
{
    writer.WriteLine("Probability,AEP,ReturnPeriod,Mode,Mean,Lower90,Upper90");
    
    for (int i = 0; i < plotProbs.Length; i++)
    {
        double p = plotProbs[i];
        double aep = 1 - p;
        double T = 1.0 / aep;
        
        writer.WriteLine($"{p:F4},{aep:F6},{T:F2}," +
                        $"{results.ModeCurve[i, 1]:F2}," +
                        $"{results.MeanCurve[i, 1]:F2}," +
                        $"{results.ConfidenceIntervals[i, 1]:F2}," +
                        $"{results.ConfidenceIntervals[i, 2]:F2}");
    }
}

Console.WriteLine("Frequency curve data exported to frequency_curve.csv");
```

## Bootstrap Parameter Distributions

Examine the distribution of estimated parameters:

```cs
var bootstrap = new BootstrapAnalysis(gev, 
    ParameterEstimationMethod.MethodOfLinearMoments,
    data.Length, 10000);

// Generate bootstrap distributions
var bootstrapDists = bootstrap.Distributions();

// Extract parameters
var parameterSets = bootstrap.ParameterSets(bootstrapDists);

// Get parameter matrix [replication, parameter]
double[,] params = bootstrap.Parameters(bootstrapDists);

Console.WriteLine("Bootstrap Parameter Statistics:");
Console.WriteLine($"Number of successful replications: {bootstrapDists.Count(d => d != null)}");

// Analyze each parameter
string[] paramNames = gev.ParameterNamesShortForm;
for (int j = 0; j < gev.NumberOfParameters; j++)
{
    var values = new List<double>();
    for (int i = 0; i < parameterSets.Length; i++)
    {
        if (parameterSets[i].Values != null)
            values.Add(parameterSets[i].Values[j]);
    }
    
    Console.WriteLine($"\n{paramNames[j]}:");
    Console.WriteLine($"  Mean: {values.Average():F4}");
    Console.WriteLine($"  Std Dev: {Statistics.StandardDeviation(values.ToArray()):F4}");
    Console.WriteLine($"  5th percentile: {Statistics.Quantile(values.OrderBy(x => x).ToArray(), 0.05):F4}");
    Console.WriteLine($"  95th percentile: {Statistics.Quantile(values.OrderBy(x => x).ToArray(), 0.95):F4}");
}
```

## Confidence Interval Methods

The ***Numerics*** library provides multiple methods for computing bootstrap confidence intervals, each with different properties.

### 1. Percentile Method (Default)

The simplest method - uses percentiles of the bootstrap distribution:

```cs
var bootstrap = new BootstrapAnalysis(gev, 
    ParameterEstimationMethod.MethodOfLinearMoments,
    data.Length, 10000);

var probabilities = new double[] { 0.9, 0.99 };

// Percentile confidence intervals
double[,] percentileCI = bootstrap.PercentileQuantileCI(
    probabilities: probabilities,
    alpha: 0.1);  // 90% CI

Console.WriteLine("Percentile Confidence Intervals:");
for (int i = 0; i < probabilities.Length; i++)
{
    Console.WriteLine($"P={probabilities[i]:F2}: [{percentileCI[i, 0]:F0}, {percentileCI[i, 1]:F0}]");
}
```

### 2. Bias-Corrected (BC) Method

Corrects for bias in the bootstrap distribution:

```cs
// Bias-corrected CI
double[,] bcCI = bootstrap.BiasCorrectedQuantileCI(probabilities, alpha: 0.1);

Console.WriteLine("Bias-Corrected Confidence Intervals:");
for (int i = 0; i < probabilities.Length; i++)
{
    Console.WriteLine($"P={probabilities[i]:F2}: [{bcCI[i, 0]:F0}, {bcCI[i, 1]:F0}]");
}
```

### 3. Normal Approximation Method

Assumes normal distribution for bootstrap statistics:

```cs
// Normal approximation CI
double[,] normalCI = bootstrap.NormalQuantileCI(probabilities, alpha: 0.1);

Console.WriteLine("Normal Approximation Confidence Intervals:");
for (int i = 0; i < probabilities.Length; i++)
{
    Console.WriteLine($"P={probabilities[i]:F2}: [{normalCI[i, 0]:F0}, {normalCI[i, 1]:F0}]");
}
```

### 4. BCa (Bias-Corrected and Accelerated)

The most accurate but computationally intensive method [[2]](#2):

```cs
// BCa method requires original sample data
double[,] bcaCI = bootstrap.BCaQuantileCI(
    sampleData: data,
    probabilities: probabilities,
    alpha: 0.1);

Console.WriteLine("BCa Confidence Intervals:");
for (int i = 0; i < probabilities.Length; i++)
{
    Console.WriteLine($"P={probabilities[i]:F2}: [{bcaCI[i, 0]:F0}, {bcaCI[i, 1]:F0}]");
}
```

### 5. Bootstrap-t Method

Uses studentized bootstrap for improved coverage:

```cs
// Bootstrap-t CI
double[,] btCI = bootstrap.BootstrapTQuantileCI(probabilities, alpha: 0.1);

Console.WriteLine("Bootstrap-t Confidence Intervals:");
for (int i = 0; i < probabilities.Length; i++)
{
    Console.WriteLine($"P={probabilities[i]:F2}: [{btCI[i, 0]:F0}, {btCI[i, 1]:F0}]");
}
```

### Comparing CI Methods

```cs
var methods = new[]
{
    ("Percentile", bootstrap.PercentileQuantileCI(probabilities, 0.1)),
    ("Bias-Corrected", bootstrap.BiasCorrectedQuantileCI(probabilities, 0.1)),
    ("Normal", bootstrap.NormalQuantileCI(probabilities, 0.1)),
    ("BCa", bootstrap.BCaQuantileCI(data, probabilities, 0.1)),
    ("Bootstrap-t", bootstrap.BootstrapTQuantileCI(probabilities, 0.1))
};

Console.WriteLine("Comparison of CI Methods (90% intervals):");
Console.WriteLine("Probability | Percentile      | BC              | Normal          | BCa             | Bootstrap-t");
Console.WriteLine("------------------------------------------------------------------------------------------------------");

for (int i = 0; i < probabilities.Length; i++)
{
    Console.Write($"{probabilities[i],11:F2} | ");
    foreach (var (name, ci) in methods)
    {
        Console.Write($"[{ci[i, 0],6:F0},{ci[i, 1],6:F0}] | ");
    }
    Console.WriteLine();
}
```

## Bootstrap Moments

Compute product moments and L-moments from bootstrap ensemble:

```cs
// Product moments from bootstrap replications
double[,] productMoments = bootstrap.ProductMoments();

Console.WriteLine("Bootstrap Product Moments:");
Console.WriteLine($"Mean of means: {productMoments[0, 0]:F2}");
Console.WriteLine($"Mean of std devs: {productMoments[1, 0]:F2}");
Console.WriteLine($"Mean of skewness: {productMoments[2, 0]:F4}");
Console.WriteLine($"Mean of kurtosis: {productMoments[3, 0]:F4}");

// L-moments from bootstrap replications
double[,] lMoments = bootstrap.LinearMoments();

Console.WriteLine("\nBootstrap L-Moments:");
Console.WriteLine($"Mean λ₁: {lMoments[0, 0]:F2}");
Console.WriteLine($"Mean λ₂: {lMoments[1, 0]:F2}");
Console.WriteLine($"Mean τ₃: {lMoments[2, 0]:F4}");
Console.WriteLine($"Mean τ₄: {lMoments[3, 0]:F4}");
```

## Expected Probability (Rare Events)

For very rare events, compute expected probabilities:

```cs
// Quantiles of interest (e.g., design floods)
var quantiles = new double[] { 15000, 20000, 25000, 30000 };

// Expected probabilities from bootstrap ensemble
double[] expectedProbs = bootstrap.ExpectedProbabilities(quantiles);

Console.WriteLine("Expected Probabilities for Design Floods:");
Console.WriteLine("Discharge | Expected AEP | Expected Return Period");
Console.WriteLine("----------------------------------------------------");

for (int i = 0; i < quantiles.Length; i++)
{
    double aep = 1 - expectedProbs[i];
    double T = 1.0 / aep;
    Console.WriteLine($"{quantiles[i],9:F0} | {aep,12:E4} | {T,22:F1} years");
}
```

## Practical Example: Complete Flood Frequency Analysis with Uncertainty

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

// Step 1: Load annual peak flow data
double[] annualPeaks = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200, 10500, 19300 };

Console.WriteLine($"Flood Frequency Analysis");
Console.WriteLine($"Sample size: {annualPeaks.Length}");
Console.WriteLine($"Sample mean: {annualPeaks.Average():F0} cfs");
Console.WriteLine($"Sample std dev: {Statistics.StandardDeviation(annualPeaks):F0} cfs");

// Step 2: Compute sample L-moments
double[] sampleLM = Statistics.LinearMoments(annualPeaks);
Console.WriteLine($"\nSample L-moments:");
Console.WriteLine($"  λ₁: {sampleLM[0]:F0}");
Console.WriteLine($"  λ₂: {sampleLM[1]:F0}");
Console.WriteLine($"  τ₃: {sampleLM[2]:F4}");

// Step 3: Fit LP3 distribution (USGS standard)
var lp3 = new LogPearsonTypeIII();
lp3.Estimate(annualPeaks, ParameterEstimationMethod.MethodOfLinearMoments);

Console.WriteLine($"\nFitted LP3 Distribution:");
Console.WriteLine($"  μ: {lp3.Mu:F4}");
Console.WriteLine($"  σ: {lp3.Sigma:F4}");
Console.WriteLine($"  γ: {lp3.Gamma:F4}");

// Step 4: Bootstrap uncertainty analysis
Console.WriteLine($"\nPerforming bootstrap analysis...");
var bootstrap = new BootstrapAnalysis(
    distribution: lp3,
    estimationMethod: ParameterEstimationMethod.MethodOfLinearMoments,
    sampleSize: annualPeaks.Length,
    replications: 10000,
    seed: 12345);

// Return periods of interest
var returnPeriods = new int[] { 2, 5, 10, 25, 50, 100, 200, 500 };
var probabilities = returnPeriods.Select(T => 1.0 - 1.0 / T).ToArray();

// Perform uncertainty analysis with 90% confidence intervals
var results = bootstrap.Estimate(probabilities, alpha: 0.1, recordParameterSets: true);

// Step 5: Present results
Console.WriteLine($"\nFlood Frequency Analysis Results (90% Confidence):");
Console.WriteLine("─────────────────────────────────────────────────────────────────────");
Console.WriteLine("Return   Annual         Point    Expected    90% Confidence Interval");
Console.WriteLine("Period   Exceedance     Estimate  Value      Lower         Upper    ");
Console.WriteLine("(years)  Probability    (cfs)     (cfs)      (cfs)         (cfs)    ");
Console.WriteLine("─────────────────────────────────────────────────────────────────────");

for (int i = 0; i < returnPeriods.Length; i++)
{
    int T = returnPeriods[i];
    double aep = 1.0 / T;
    double point = results.ModeCurve[i, 1];
    double mean = results.MeanCurve[i, 1];
    double lower = results.ConfidenceIntervals[i, 1];
    double upper = results.ConfidenceIntervals[i, 2];
    
    Console.WriteLine($"{T,6}   {aep,11:F5}      {point,8:F0}  {mean,8:F0}   {lower,8:F0}      {upper,8:F0}");
}

Console.WriteLine("─────────────────────────────────────────────────────────────────────");

// Step 6: Parameter uncertainty
Console.WriteLine($"\nParameter Uncertainty:");
Console.WriteLine($"Parameter sets recorded: {results.ParameterSets.Length}");

double[] muValues = results.ParameterSets.Select(ps => ps.Values[0]).ToArray();
double[] sigmaValues = results.ParameterSets.Select(ps => ps.Values[1]).ToArray();
double[] gammaValues = results.ParameterSets.Select(ps => ps.Values[2]).ToArray();

Console.WriteLine($"μ:  {muValues.Average():F4} ± {Statistics.StandardDeviation(muValues):F4}");
Console.WriteLine($"σ:  {sigmaValues.Average():F4} ± {Statistics.StandardDeviation(sigmaValues):F4}");
Console.WriteLine($"γ:  {gammaValues.Average():F4} ± {Statistics.StandardDeviation(gammaValues):F4}");

// Step 7: Uncertainty bounds width
Console.WriteLine($"\nUncertainty Analysis Summary:");
for (int i = 0; i < returnPeriods.Length; i++)
{
    double width = results.ConfidenceIntervals[i, 2] - results.ConfidenceIntervals[i, 1];
    double relativeWidth = width / results.ModeCurve[i, 1] * 100;
    Console.WriteLine($"{returnPeriods[i]}-year: ±{relativeWidth:F1}% relative uncertainty");
}
```

## Advanced Bootstrap Techniques

### Reusing Bootstrap Samples

Generate bootstrap distributions once and reuse for multiple analyses:

```cs
// Generate bootstrap distributions once
var bootstrapDists = bootstrap.Distributions();

Console.WriteLine($"Generated {bootstrapDists.Count(d => d != null)} valid bootstrap replications");

// Reuse for different probability sets
var probs1 = new double[] { 0.9, 0.95, 0.99 };
var results1 = bootstrap.Estimate(probs1, alpha: 0.1, distributions: bootstrapDists);

var probs2 = new double[] { 0.5, 0.75, 0.98 };
var results2 = bootstrap.Estimate(probs2, alpha: 0.1, distributions: bootstrapDists);

// Much faster since bootstrap samples are reused
```

### Custom Quantile Computations

```cs
// Compute quantiles from bootstrap ensemble
var probsOfInterest = new double[] { 0.9, 0.95, 0.98, 0.99, 0.998 };

double[,] quantiles = bootstrap.Quantiles(probsOfInterest);

Console.WriteLine("Bootstrap Quantiles:");
Console.WriteLine("Prob   | Mean      | Std Dev   | 5th %ile  | 95th %ile");
Console.WriteLine("------------------------------------------------------------");

for (int i = 0; i < probsOfInterest.Length; i++)
{
    Console.WriteLine($"{probsOfInterest[i]:F3} | {quantiles[i, 0],9:F0} | " +
                     $"{quantiles[i, 1],9:F0} | {quantiles[i, 2],9:F0} | {quantiles[i, 3],9:F0}");
}
```

### Computing Probabilities

Reverse direction - find probabilities for given quantiles:

```cs
var designFlows = new double[] { 15000, 20000, 25000, 30000 };

double[,] probabilities = bootstrap.Probabilities(designFlows);

Console.WriteLine("Probabilities for Design Flows:");
Console.WriteLine("Flow  | Mean Prob | Std Dev   | 5th %ile  | 95th %ile");
Console.WriteLine("--------------------------------------------------------");

for (int i = 0; i < designFlows.Length; i++)
{
    double meanAEP = 1 - probabilities[i, 0];
    double meanT = 1.0 / meanAEP;
    
    Console.WriteLine($"{designFlows[i],5:F0} | {probabilities[i, 0],9:F4} | " +
                     $"{probabilities[i, 1],9:F4} | {probabilities[i, 2],9:F4} | {probabilities[i, 3],9:F4}");
    Console.WriteLine($"      | T={meanT:F1} years");
}
```

## Choosing Number of Replications

```cs
// Assess convergence with different replication counts
var replicationCounts = new int[] { 1000, 5000, 10000, 20000 };
var testProb = 0.99;  // 100-year event

Console.WriteLine("Bootstrap Convergence Analysis:");
Console.WriteLine("Replications | 100-yr Estimate | CI Width | Time (relative)");

foreach (var nRep in replicationCounts)
{
    var boot = new BootstrapAnalysis(gev, 
        ParameterEstimationMethod.MethodOfLinearMoments,
        data.Length, nRep);
    
    var result = boot.Estimate(new[] { testProb }, alpha: 0.1);
    
    double estimate = result.ModeCurve[0, 1];
    double width = result.ConfidenceIntervals[0, 2] - result.ConfidenceIntervals[0, 1];
    
    Console.WriteLine($"{nRep,12} | {estimate,15:F0} | {width,8:F0} | {nRep / 1000.0,4:F1}×");
}

// Recommendation: 10,000-20,000 replications for most applications
// 50,000+ for critical infrastructure or very rare events
```

## Bootstrap vs. Analytical Uncertainty

```cs
// Compare bootstrap CI with analytical (if available)
var normal = new Normal(100, 15);
normal.Estimate(data, ParameterEstimationMethod.MethodOfMoments);

// Bootstrap CI
var bootstrap = new BootstrapAnalysis(normal, 
    ParameterEstimationMethod.MethodOfMoments,
    data.Length, 10000);

var probs = new double[] { 0.5, 0.9, 0.95 };
var bootResults = bootstrap.Estimate(probs, alpha: 0.1);

// Analytical CI (for normal distribution, if implemented)
// This would require separate implementation

Console.WriteLine("Bootstrap vs Analytical Confidence Intervals:");
Console.WriteLine("Both methods should give similar results for Normal distribution");
Console.WriteLine("Bootstrap is more general and works for any distribution");
```

## Best Practices

### 1. Always Use Adequate Replications

```cs
// Too few replications give unstable confidence intervals
int nReplications = Math.Max(10000, 20 * data.Length);
Console.WriteLine($"Recommended replications: {nReplications}");
```

### 2. Check for Failed Replications

```cs
var bootstrapDists = bootstrap.Distributions();
int nFailed = bootstrapDists.Count(d => d == null);

if (nFailed > bootstrapDists.Length * 0.05)
{
    Console.WriteLine($"Warning: {nFailed} replications failed ({100.0 * nFailed / bootstrapDists.Length:F1}%)");
    Console.WriteLine("Consider using more robust estimation method");
}
```

### 3. Report Uncertainty Appropriately

```cs
// Report point estimate ± uncertainty
double point = results.ModeCurve[0, 1];
double lower = results.ConfidenceIntervals[0, 1];
double upper = results.ConfidenceIntervals[0, 2];

Console.WriteLine($"100-year flood: {point:F0} cfs");
Console.WriteLine($"90% CI: [{lower:F0}, {upper:F0}] cfs");
Console.WriteLine($"Or: {point:F0} ± {(upper - point):F0}/-{(point - lower):F0} cfs");
```

### 4. Consider Sample Size Effects

```cs
// Uncertainty increases as sample size decreases
Console.WriteLine("Effect of Sample Size on Uncertainty:");

foreach (var n in new[] { 10, 20, 50, 100 })
{
    // Simulate sampling from fitted distribution
    var sample = gev.GenerateRandomValues(n);
    var testDist = new GeneralizedExtremeValue();
    testDist.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
    
    var testBoot = new BootstrapAnalysis(testDist, 
        ParameterEstimationMethod.MethodOfLinearMoments, n, 1000);
    var testResults = testBoot.Estimate(new[] { 0.99 }, alpha: 0.1);
    
    double width = testResults.ConfidenceIntervals[0, 2] - testResults.ConfidenceIntervals[0, 1];
    Console.WriteLine($"n={n,3}: CI width = {width:F0}");
}
```

## Common Pitfalls

1. **Too few replications** - Use at least 10,000 for final analyses
2. **Wrong estimation method** - Use same method for parent and bootstrap
3. **Not checking convergence** - Monitor failed replications
4. **Ignoring small sample bias** - Bootstrap can't fix fundamental data limitations
5. **Overinterpreting precision** - CI width reflects sampling uncertainty only

---

## References

<a id="1">[1]</a> Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the Bootstrap*. Chapman & Hall/CRC.

<a id="2">[2]</a> Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the American Statistical Association*, 82(397), 171-185.

<a id="3">[3]</a> Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap Methods and Their Application*. Cambridge University Press.

---

[← Previous: Parameter Estimation](parameter-estimation.md) | [Back to Index](../index.md) | [Next: Copulas →](copulas.md)
