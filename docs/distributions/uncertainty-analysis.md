# Uncertainty Analysis

[← Previous: Parameter Estimation](parameter-estimation.md) | [Back to Index](../index.md) | [Next: Copulas →](copulas.md)

Uncertainty analysis quantifies the confidence in estimated distribution parameters and derived quantities (like design floods or failure probabilities). The ***Numerics*** library provides comprehensive bootstrap resampling methods for assessing parameter and quantile uncertainty, which are essential for risk-informed decision making.

## Bootstrap Analysis Overview

Bootstrap resampling [[1]](#1) is a powerful method for estimating sampling distributions and confidence intervals. The core idea is deceptively simple: rather than deriving the sampling distribution of an estimator analytically, approximate it empirically by repeatedly resampling and re-estimating.

The ***Numerics*** library implements a **parametric bootstrap**, which generates new samples from the fitted distribution rather than resampling the original observations with replacement. This is the natural approach when working with fitted probability distributions, and is particularly effective for extreme value analysis where the fitted distribution captures tail behavior that the limited observed data cannot.

### The Parametric Bootstrap Algorithm

Given observed data $x_1, x_2, \ldots, x_n$ and a fitted distribution $\hat{F}$ with estimated parameters $\hat{\theta}$, the parametric bootstrap proceeds as follows:

**Step 1.** Generate $B$ bootstrap samples, each of size $n$, by drawing from the fitted distribution:

```math
x^{*}_{b,1}, x^{*}_{b,2}, \ldots, x^{*}_{b,n} \sim \hat{F}(\hat{\theta}), \quad b = 1, 2, \ldots, B
```

**Step 2.** For each bootstrap sample $b$, re-estimate the distribution parameters using the same estimation method (e.g., L-moments, MLE) to obtain $\hat{\theta}^{*}_b$.

**Step 3.** Compute the statistic of interest from each refitted distribution. For quantile estimation at non-exceedance probability $p$:

```math
\hat{Q}^{*}_b = \hat{F}^{-1}(p \mid \hat{\theta}^{*}_b), \quad b = 1, 2, \ldots, B
```

**Step 4.** The empirical distribution of $\lbrace\hat{Q}^{*}_1, \hat{Q}^{*}_2, \ldots, \hat{Q}^{*}_B\rbrace$ approximates the sampling distribution of the quantile estimator $\hat{Q}$.

From this bootstrap distribution, we can compute several useful summaries. The **bootstrap standard error** is the sample standard deviation of the bootstrap replicates:

```math
\widehat{SE}_{boot} = \sqrt{\frac{1}{B-1}\sum_{b=1}^{B}\left(\hat{Q}^{*}_b - \bar{Q}^{*}\right)^2}
```

where $\bar{Q}^{*} = \frac{1}{B}\sum_{b=1}^{B}\hat{Q}^{*}_b$ is the mean of the bootstrap replicates. The **bootstrap estimate of bias** is:

```math
\widehat{bias} = \bar{Q}^{*} - \hat{Q}
```

where $\hat{Q}$ is the original point estimate from the parent distribution.

### Why Parametric Bootstrap?

The parametric bootstrap is preferred over the nonparametric bootstrap for flood frequency analysis and similar applications because:
- It leverages the assumed distributional form to generate realistic samples, including plausible extreme values
- It produces smoother confidence intervals in the tails where data are sparse
- It naturally handles the small sample sizes typical of annual peak flow records
- It maintains consistency with the fitted distribution model used for design

The bootstrap is particularly valuable when:
- Analytical confidence intervals are unavailable or intractable
- Sample sizes are small to moderate (common in hydrology)
- The sampling distribution of the estimator is unknown or complex
- Dealing with extreme value distributions where tail uncertainty is critical

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
double[] ModeCurve  // ModeCurve[i] = quantile at probabilities[i]

// Mean curve: expected quantiles from bootstrap ensemble
double[] MeanCurve  // MeanCurve[i] = expected quantile at probabilities[i]

// Confidence intervals for quantiles
double[,] ConfidenceIntervals  // [i, 0] = lower, [i, 1] = upper

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
    
    double mode = results.ModeCurve[i];          // Point estimate
    double mean = results.MeanCurve[i];          // Expected value
    double lower = results.ConfidenceIntervals[i, 0];  // Lower bound
    double upper = results.ConfidenceIntervals[i, 1];  // Upper bound
    
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
                        $"{results.ModeCurve[i]:F2}," +
                        $"{results.MeanCurve[i]:F2}," +
                        $"{results.ConfidenceIntervals[i, 0]:F2}," +
                        $"{results.ConfidenceIntervals[i, 1]:F2}");
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
    Console.WriteLine($"  5th percentile: {Statistics.Percentile(values.ToArray(), 0.05):F4}");
    Console.WriteLine($"  95th percentile: {Statistics.Percentile(values.ToArray(), 0.95):F4}");
}
```

## Confidence Interval Methods

The ***Numerics*** library provides five methods for computing bootstrap confidence intervals, each with different statistical properties. They range from simple (Percentile) to sophisticated (BCa, Bootstrap-t), trading off computational cost against accuracy of coverage. All methods construct a $(1-\alpha)\times 100\%$ confidence interval for a quantile $\hat{Q}$ at a given non-exceedance probability.

### 1. Percentile Method (Default)

The Percentile method [[1]](#1) is the simplest bootstrap confidence interval. It uses the quantiles of the bootstrap distribution directly as confidence limits:

```math
CI_{1-\alpha} = \left[\hat{Q}^{*}_{(\alpha/2)},\;\hat{Q}^{*}_{(1-\alpha/2)}\right]
```

where $\hat{Q}^{*}_{(p)}$ denotes the $p$-th percentile of the bootstrap distribution $\lbrace\hat{Q}^{*}_1, \ldots, \hat{Q}^{*}_B\rbrace$. For a 90% confidence interval ($\alpha = 0.1$), this takes the 5th and 95th percentiles of the bootstrap replicates.

The Percentile method is intuitive and easy to implement, but it does **not** correct for bias or skewness in the bootstrap distribution. It works well when the bootstrap distribution is approximately symmetric and the estimator is approximately unbiased. This is the default method used by the `Estimate()` method.

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

The Bias-Corrected (BC) method [[2]](#2) adjusts the percentiles used for the confidence interval to correct for median bias in the bootstrap distribution. If the estimator is biased, the simple Percentile method produces intervals that are shifted; the BC method fixes this.

First, compute the bias-correction factor $z_0$, which measures how far the bootstrap distribution median is from the original estimate:

```math
z_0 = \Phi^{-1}\!\left(\frac{\#\{\hat{Q}^{*}_b \le \hat{Q}\}}{B + 1}\right)
```

where $\Phi^{-1}$ is the inverse of the standard Normal CDF, $\hat{Q}$ is the original point estimate, and the fraction counts the proportion of bootstrap replicates that fall at or below the original estimate. If the estimator is unbiased, approximately half the replicates will be below $\hat{Q}$, giving $z_0 \approx 0$.

The adjusted percentiles for the confidence interval are:

```math
\alpha_1 = \Phi\!\left(2z_0 + z_{\alpha/2}\right), \quad \alpha_2 = \Phi\!\left(2z_0 + z_{1-\alpha/2}\right)
```

where $z_{\alpha/2} = \Phi^{-1}(\alpha/2)$ is the standard Normal quantile. The confidence interval is then:

```math
CI_{1-\alpha} = \left[\hat{Q}^{*}_{(\alpha_1)},\;\hat{Q}^{*}_{(\alpha_2)}\right]
```

When $z_0 = 0$ (no bias), this reduces to the standard Percentile method. The BC method is better than the Percentile method when the estimator has median bias, which is common for quantile estimators of skewed distributions.

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

The Normal method assumes the bootstrap distribution of the statistic is approximately Gaussian. In ***Numerics***, this method applies a cube-root transformation to improve normality before computing the interval, then back-transforms the result:

```math
\tilde{Q} = \hat{Q}^{1/3}, \quad \tilde{Q}^{*}_b = \left(\hat{Q}^{*}_b\right)^{1/3}
```

The bootstrap standard error is computed on the transformed scale:

```math
\widetilde{SE} = \sqrt{\frac{1}{B-1}\sum_{b=1}^{B}\left(\tilde{Q}^{*}_b - \bar{\tilde{Q}}^{*}\right)^2}
```

The confidence interval is constructed in the transformed space and back-transformed:

```math
CI_{1-\alpha} = \left[\left(\tilde{Q} + z_{\alpha/2}\cdot\widetilde{SE}\right)^3,\;\left(\tilde{Q} + z_{1-\alpha/2}\cdot\widetilde{SE}\right)^3\right]
```

where $z_{\alpha/2} = \Phi^{-1}(\alpha/2)$ is the standard Normal quantile. The cube-root transformation is a variance-stabilizing power transformation that makes the method approximately transformation-invariant, improving performance for skewed distributions common in hydrology.

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

The BCa method [[2]](#2) extends the BC method by adding an **acceleration constant** $\hat{a}$ that corrects for skewness in the bootstrap distribution. This makes the BCa interval **second-order accurate** and **transformation-respecting** -- meaning it gives the same answer regardless of what monotone transformation is applied to the data.

The bias-correction factor $z_0$ is computed the same way as in the BC method. The acceleration constant $\hat{a}$ is estimated using the jackknife. For each observation $i = 1, \ldots, n$, the distribution is re-estimated with that observation removed and the statistic is recomputed, yielding jackknife replicates $\hat{Q}_{(i)}$. The acceleration constant is:

```math
\hat{a} = \frac{\sum_{i=1}^{n}\left(\hat{Q} - \hat{Q}_{(i)}\right)^3}{6\left[\sum_{i=1}^{n}\left(\hat{Q} - \hat{Q}_{(i)}\right)^2\right]^{3/2}}
```

where $\hat{Q}$ is the original point estimate and $\hat{Q}_{(i)}$ is the estimate computed from the sample with observation $i$ removed. The acceleration constant measures the rate at which the standard error of $\hat{Q}$ changes as the true parameter value changes.

The adjusted percentiles incorporate both bias correction and acceleration:

```math
\alpha_1 = \Phi\!\left(z_0 + \frac{z_0 + z_{\alpha/2}}{1 - \hat{a}(z_0 + z_{\alpha/2})}\right), \quad \alpha_2 = \Phi\!\left(z_0 + \frac{z_0 + z_{1-\alpha/2}}{1 - \hat{a}(z_0 + z_{1-\alpha/2})}\right)
```

The confidence interval is then $CI_{1-\alpha} = [\hat{Q}^{*}_{(\alpha_1)}, \hat{Q}^{*}_{(\alpha_2)}]$.

When $\hat{a} = 0$, the BCa method reduces to the BC method. When both $z_0 = 0$ and $\hat{a} = 0$, it reduces to the Percentile method. The BCa method is the most accurate of the percentile-based methods, but it requires the original sample data and is computationally expensive due to the jackknife (which requires $n$ additional distribution fits).

> **Note:** The `BCaQuantileCI` method requires the original sample data because it re-estimates the distribution parameters internally and performs the jackknife. This means it calls `Estimate()` on the distribution, which will update the distribution's parameters.

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

The Bootstrap-t (studentized bootstrap) method [[3]](#3) is the bootstrap analog of the Student-t confidence interval. Rather than using percentiles of the bootstrap distribution directly, it constructs a pivotal quantity by standardizing each bootstrap replicate by its own standard error. This approach typically achieves the most accurate coverage for location-type parameters.

Like the Normal method, ***Numerics*** applies a cube-root transformation for variance stabilization. For each bootstrap replicate $b$, the method computes both the transformed quantile estimate $\tilde{Q}^{*}_b = (\hat{Q}^{*}_b)^{1/3}$ and its standard error $\widetilde{SE}^{*}_b$ (estimated via an inner bootstrap of 300 replications). The studentized statistic is:

```math
t^{*}_b = \frac{\tilde{Q} - \tilde{Q}^{*}_b}{\widetilde{SE}^{*}_b}, \quad b = 1, \ldots, B
```

where $\tilde{Q} = \hat{Q}^{1/3}$ is the transformed original estimate. The confidence interval is constructed using the percentiles of the studentized distribution and the overall bootstrap standard error $\widetilde{SE}$:

```math
CI_{1-\alpha} = \left[\left(\tilde{Q} + t^{*}_{(\alpha/2)}\cdot\widetilde{SE}\right)^3,\;\left(\tilde{Q} + t^{*}_{(1-\alpha/2)}\cdot\widetilde{SE}\right)^3\right]
```

where $t^{*}_{(p)}$ is the $p$-th percentile of $\lbrace t^{*}_1, \ldots, t^{*}_B\rbrace$, and $\widetilde{SE}$ is the standard deviation of the transformed bootstrap replicates $\lbrace\tilde{Q}^{*}_1, \ldots, \tilde{Q}^{*}_B\rbrace$.

The Bootstrap-t method is the most computationally expensive method because it requires a **double bootstrap**: each of the $B$ outer replications requires an inner bootstrap (300 replications by default) to estimate the standard error. However, it can provide the most accurate coverage probabilities for location parameters and is second-order accurate.

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

Compute product moments and L-moments from the bootstrap ensemble. Both methods return a `double[Replications, 4]` array where each row is one bootstrap replication and each column is a moment:

```cs
using Numerics.Data.Statistics;

// Product moments: [replication, moment]
//   Column 0 = mean, 1 = std dev, 2 = skewness, 3 = kurtosis
double[,] productMoments = bootstrap.ProductMoments();

int R = productMoments.GetLength(0);  // number of replications

// Summarize across replications for each moment
Console.WriteLine("Bootstrap Product Moments:");
Console.WriteLine($"Mean of means:    {Enumerable.Range(0, R).Average(i => productMoments[i, 0]):F2}");
Console.WriteLine($"Mean of std devs: {Enumerable.Range(0, R).Average(i => productMoments[i, 1]):F2}");
Console.WriteLine($"Mean of skewness: {Enumerable.Range(0, R).Average(i => productMoments[i, 2]):F4}");
Console.WriteLine($"Mean of kurtosis: {Enumerable.Range(0, R).Average(i => productMoments[i, 3]):F4}");

// L-moments: [replication, moment]
//   Column 0 = λ₁, 1 = λ₂, 2 = τ₃, 3 = τ₄
double[,] lMoments = bootstrap.LinearMoments();

Console.WriteLine("\nBootstrap L-Moments:");
Console.WriteLine($"Mean λ₁: {Enumerable.Range(0, R).Average(i => lMoments[i, 0]):F2}");
Console.WriteLine($"Mean λ₂: {Enumerable.Range(0, R).Average(i => lMoments[i, 1]):F2}");
Console.WriteLine($"Mean τ₃: {Enumerable.Range(0, R).Average(i => lMoments[i, 2]):F4}");
Console.WriteLine($"Mean τ₄: {Enumerable.Range(0, R).Average(i => lMoments[i, 3]):F4}");
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
    double point = results.ModeCurve[i];
    double mean = results.MeanCurve[i];
    double lower = results.ConfidenceIntervals[i, 0];
    double upper = results.ConfidenceIntervals[i, 1];
    
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
    double width = results.ConfidenceIntervals[i, 1] - results.ConfidenceIntervals[i, 0];
    double relativeWidth = width / results.ModeCurve[i] * 100;
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

The `Quantiles()` method returns a `double[Replications, probabilities.Count]` array containing the raw quantile value from each bootstrap replication at each probability. You can then compute your own summary statistics:

```cs
using Numerics.Data.Statistics;

// Compute quantiles from bootstrap ensemble
var probsOfInterest = new double[] { 0.9, 0.95, 0.98, 0.99, 0.998 };

// Returns double[Replications, probabilities.Count]
double[,] quantiles = bootstrap.Quantiles(probsOfInterest);

int R = quantiles.GetLength(0);  // number of replications

Console.WriteLine("Bootstrap Quantiles:");
Console.WriteLine("Prob   | Mean      | Std Dev   | 5th %ile  | 95th %ile");
Console.WriteLine("------------------------------------------------------------");

for (int i = 0; i < probsOfInterest.Length; i++)
{
    // Extract column for this probability across all replications
    double[] values = Enumerable.Range(0, R)
        .Select(r => quantiles[r, i])
        .Where(v => !double.IsNaN(v))
        .ToArray();

    double mean = values.Average();
    double sd = Statistics.StandardDeviation(values);
    double p05 = Statistics.Percentile(values, 0.05);
    double p95 = Statistics.Percentile(values, 0.95);

    Console.WriteLine($"{probsOfInterest[i]:F3} | {mean,9:F0} | {sd,9:F0} | {p05,9:F0} | {p95,9:F0}");
}
```

### Computing Probabilities

Reverse direction — find CDF probabilities for given quantile values. The `Probabilities()` method returns a `double[Replications, quantiles.Count]` array of raw CDF values from each replication:

```cs
var designFlows = new double[] { 15000, 20000, 25000, 30000 };

// Returns double[Replications, quantiles.Count]
double[,] probs = bootstrap.Probabilities(designFlows);

int R = probs.GetLength(0);

Console.WriteLine("Probabilities for Design Flows:");
Console.WriteLine("Flow  | Mean Prob | Std Dev   | 5th %ile  | 95th %ile");
Console.WriteLine("--------------------------------------------------------");

for (int i = 0; i < designFlows.Length; i++)
{
    // Extract column for this quantile across all replications
    double[] values = Enumerable.Range(0, R)
        .Select(r => probs[r, i])
        .Where(v => !double.IsNaN(v))
        .ToArray();

    double meanProb = values.Average();
    double sd = Statistics.StandardDeviation(values);
    double p05 = Statistics.Percentile(values, 0.05);
    double p95 = Statistics.Percentile(values, 0.95);

    double meanAEP = 1 - meanProb;
    double meanT = 1.0 / meanAEP;

    Console.WriteLine($"{designFlows[i],5:F0} | {meanProb,9:F4} | {sd,9:F4} | {p05,9:F4} | {p95,9:F4}");
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
    
    double estimate = result.ModeCurve[0];
    double width = result.ConfidenceIntervals[0, 1] - result.ConfidenceIntervals[0, 0];
    
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
double point = results.ModeCurve[0];
double lower = results.ConfidenceIntervals[0, 0];
double upper = results.ConfidenceIntervals[0, 1];

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
    
    double width = testResults.ConfidenceIntervals[0, 1] - testResults.ConfidenceIntervals[0, 0];
    Console.WriteLine($"n={n,3}: CI width = {width:F0}");
}
```

## Common Pitfalls

1. **Too few replications** - Use at least 10,000 for final analyses
2. **Wrong estimation method** - Use same method for parent and bootstrap
3. **Not checking convergence** - Monitor failed replications
4. **Ignoring small sample bias** - Bootstrap can't fix fundamental data limitations
5. **Overinterpreting precision** - CI width reflects sampling uncertainty only

## Rules of Thumb for Number of Replications

The number of bootstrap replications $B$ directly affects the stability and precision of the results. The following guidelines apply:

| Purpose | Minimum $B$ | Recommended $B$ |
|---|---|---|
| Standard error estimation | 200 | 1,000 |
| Confidence intervals | 1,000 | 10,000 |
| Precise CI endpoints | 5,000 | 20,000--50,000 |
| Extreme quantiles (99th percentile and beyond) | 10,000 | 50,000+ |

For life-safety applications, use at least 10,000 replications for confidence intervals and verify convergence by comparing results at different replication counts. The ***Numerics*** library enforces a minimum of 100 replications and a minimum sample size of 10.

## When Bootstrap Fails

The bootstrap is not a universal solution. Practitioners working on safety-critical projects must be aware of these limitations:

- **Very small samples ($n < 15$):** The bootstrap distribution is a poor approximation of the true sampling distribution because the fitted parametric model may itself be unreliable. Confidence intervals will be too narrow (overly optimistic). The ***Numerics*** library enforces $n \ge 10$.

- **Heavy-tailed distributions:** When the underlying distribution has very heavy tails (e.g., GEV with large positive shape parameter), the bootstrap variance estimate can itself be highly variable. More replications are needed, and results should be interpreted cautiously.

- **Dependent data:** The standard parametric bootstrap assumes that observations are independent and identically distributed (i.i.d.). For time series data with serial correlation, the bootstrap underestimates uncertainty. Block bootstrap or other specialized methods are needed for dependent data.

- **Parameters near boundary of parameter space:** If the true parameter lies on or near the boundary of the parameter space (e.g., shape parameter near zero for GEV), the bootstrap distribution may be inconsistent. Some bootstrap replications may fail to converge, which is why the library allows up to 20 retries per replication.

- **Model misspecification:** The parametric bootstrap inherits any bias from the assumed distributional form. If the wrong distribution family is used, the bootstrap confidence intervals may have poor coverage even with large $B$.

## Confidence Interval Method Selection Guide

Choosing the right confidence interval method involves balancing accuracy against computational cost:

| Method | Strengths | Weaknesses | Best For |
|---|---|---|---|
| **Percentile** | Simplest, fastest | No bias/skewness correction | Quick estimates; symmetric, unbiased cases |
| **Normal** | Simple, uses cube-root transform | Assumes approximate normality | Distributions with near-Normal quantile estimators |
| **Bias-Corrected (BC)** | Corrects for median bias | Does not correct for skewness | Moderately biased estimators |
| **BCa** | Second-order accurate; handles bias and skewness | Expensive (jackknife); needs original data | Best general-purpose accuracy |
| **Bootstrap-t** | Most accurate for location parameters; second-order accurate | Most expensive (double bootstrap) | Location parameters; when coverage accuracy is critical |

**Practical recommendations:**

- For **routine analyses**, the Percentile method (default) is adequate and fastest.
- For **design-level analyses**, use BC or BCa for improved accuracy.
- For **critical infrastructure** where confidence interval coverage must be precise, consider BCa or Bootstrap-t with $B \ge 20{,}000$.
- When **computational cost matters**, the Normal method offers a good balance of speed and accuracy through its cube-root variance-stabilizing transform.
- When in doubt, **compare methods**: if all five methods give similar intervals, any method is adequate. If they disagree substantially, the BCa or Bootstrap-t results should be preferred.

---

## References

<a id="1">[1]</a> Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the Bootstrap*. Chapman & Hall/CRC.

<a id="2">[2]</a> Efron, B. (1987). Better bootstrap confidence intervals. *Journal of the American Statistical Association*, 82(397), 171-185.

<a id="3">[3]</a> Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap Methods and Their Application*. Cambridge University Press.

---

[← Previous: Parameter Estimation](parameter-estimation.md) | [Back to Index](../index.md) | [Next: Copulas →](copulas.md)
