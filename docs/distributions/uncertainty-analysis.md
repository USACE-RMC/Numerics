# Uncertainty Analysis

This document describes methods for quantifying uncertainty in distribution parameter estimates and derived quantities such as quantiles.

## Overview

Uncertainty in frequency analysis arises from:
1. **Sampling uncertainty**: Limited sample size
2. **Model uncertainty**: Choice of distribution family
3. **Measurement uncertainty**: Errors in observed data

***Numerics*** provides tools for quantifying sampling uncertainty through analytical standard errors and bootstrap resampling.

## Analytical Standard Errors

For distributions implementing `IStandardError`, analytical formulas provide variance estimates for parameters and quantiles.

### Parameter Covariance Matrix

The parameter covariance matrix captures the variance of and correlation between parameter estimates:

```cs
using Numerics.Distributions;

var normal = new Normal(100, 15);
int sampleSize = 50;

// Get parameter covariance matrix
double[,] covar = normal.ParameterCovariance(
    sampleSize, 
    ParameterEstimationMethod.MaximumLikelihood
);

double varMu = covar[0, 0];    // Variance of mean estimate
double varSigma = covar[1, 1]; // Variance of std dev estimate
double covMuSigma = covar[0, 1]; // Covariance

Console.WriteLine($"SE(μ) = {Math.Sqrt(varMu):F3}");
Console.WriteLine($"SE(σ) = {Math.Sqrt(varSigma):F3}");
```

### Quantile Standard Error

The standard error of a quantile estimate is computed using the delta method [[1]](#ref1):

```math
\text{Var}(\hat{x}_p) = \nabla Q_p^T \cdot \Sigma \cdot \nabla Q_p
```

where $\nabla Q_p$ is the gradient of the quantile with respect to parameters and $\Sigma$ is the parameter covariance matrix.

```cs
var gev = new GeneralizedExtremeValue(1000, 200, -0.1);
int n = 50;

// Standard error of 100-year flood estimate
double se100 = gev.QuantileStandardError(
    0.99,  // Exceedance probability for 100-year event
    n,
    ParameterEstimationMethod.MethodOfLinearMoments
);

double q100 = gev.InverseCDF(0.99);
Console.WriteLine($"Q100 = {q100:F0} ± {1.96 * se100:F0} (95% CI)");
```

### Confidence Intervals

Construct confidence intervals for quantiles:

```cs
double alpha = 0.05;  // 95% confidence
double zAlpha = 1.96;

double q100 = gev.InverseCDF(0.99);
double se = gev.QuantileStandardError(0.99, n, method);

// Normal approximation CI
double lowerCI = q100 - zAlpha * se;
double upperCI = q100 + zAlpha * se;

Console.WriteLine($"95% CI: [{lowerCI:F0}, {upperCI:F0}]");
```

---

## Bootstrap Methods

Bootstrap resampling provides distribution-free uncertainty estimates [[2]](#ref2).

### Parametric Bootstrap

Resample from the fitted distribution:

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

double[] data = { 12500, 15200, 11800, 18900, 14200, 16500 };

// Fit original distribution
var gev = new GeneralizedExtremeValue();
gev.SetParameters(gev.ParametersFromLinearMoments(data));

// Bootstrap settings
int nBootstrap = 1000;
int n = data.Length;
var rng = new Random(12345);

double[] q100Bootstrap = new double[nBootstrap];

for (int b = 0; b < nBootstrap; b++)
{
    // Generate bootstrap sample from fitted distribution
    double[] bootSample = new double[n];
    gev.GenerateRandomValues(bootSample, seed: rng.Next());
    
    // Refit distribution to bootstrap sample
    var bootGEV = new GeneralizedExtremeValue();
    bootGEV.SetParameters(bootGEV.ParametersFromLinearMoments(bootSample));
    
    // Store bootstrap quantile estimate
    q100Bootstrap[b] = bootGEV.InverseCDF(0.99);
}

// Bootstrap confidence interval (percentile method)
Array.Sort(q100Bootstrap);
double lower = Statistics.Percentile(q100Bootstrap, 0.025);
double upper = Statistics.Percentile(q100Bootstrap, 0.975);

Console.WriteLine($"Q100 = {gev.InverseCDF(0.99):F0}");
Console.WriteLine($"Bootstrap 95% CI: [{lower:F0}, {upper:F0}]");
Console.WriteLine($"Bootstrap SE: {Statistics.StandardDeviation(q100Bootstrap):F0}");
```

### Non-Parametric Bootstrap

Resample with replacement from the original data:

```cs
double[] q100Bootstrap = new double[nBootstrap];

for (int b = 0; b < nBootstrap; b++)
{
    // Resample with replacement
    double[] bootSample = new double[n];
    for (int i = 0; i < n; i++)
    {
        int idx = rng.Next(n);
        bootSample[i] = data[idx];
    }
    
    // Refit distribution
    var bootGEV = new GeneralizedExtremeValue();
    bootGEV.SetParameters(bootGEV.ParametersFromLinearMoments(bootSample));
    
    q100Bootstrap[b] = bootGEV.InverseCDF(0.99);
}
```

### Bootstrap Confidence Interval Methods

| Method | Description | Code |
|--------|-------------|------|
| Percentile | Direct percentiles of bootstrap distribution | `Percentile(boots, α/2)` to `Percentile(boots, 1-α/2)` |
| Basic | $2\hat{\theta} - \theta^*_{1-\alpha/2}$ to $2\hat{\theta} - \theta^*_{\alpha/2}$ | Reflects around estimate |
| BCa | Bias-corrected and accelerated | Adjusts for bias and skewness |

```cs
// Percentile method (simplest)
double ciLowerPercentile = Statistics.Percentile(q100Bootstrap, 0.025);
double ciUpperPercentile = Statistics.Percentile(q100Bootstrap, 0.975);

// Basic method
double q100Original = gev.InverseCDF(0.99);
double ciLowerBasic = 2 * q100Original - Statistics.Percentile(q100Bootstrap, 0.975);
double ciUpperBasic = 2 * q100Original - Statistics.Percentile(q100Bootstrap, 0.025);

Console.WriteLine($"Percentile CI: [{ciLowerPercentile:F0}, {ciUpperPercentile:F0}]");
Console.WriteLine($"Basic CI:      [{ciLowerBasic:F0}, {ciUpperBasic:F0}]");
```

---

## Expected Probability Analysis

Account for uncertainty when computing expected (average) exceedance probabilities [[3]](#ref3).

### Theory

The expected probability integrates over parameter uncertainty:

```math
\bar{P}(X > x) = \int P(X > x | \theta) \pi(\theta | \text{data}) d\theta
```

For quantiles, the expected probability curve lies above the fitted curve in the tail.

### Implementation

```cs
// Generate parameter samples via MCMC or bootstrap
double[,] parameterSamples = new double[nSamples, nParams];
// ... fill from MCMC or bootstrap ...

// Compute expected exceedance probability at a given flow
double flow = 25000;
double expectedProb = 0;

for (int i = 0; i < nSamples; i++)
{
    var tempGEV = new GeneralizedExtremeValue(
        parameterSamples[i, 0],
        parameterSamples[i, 1],
        parameterSamples[i, 2]
    );
    expectedProb += tempGEV.CCDF(flow);
}
expectedProb /= nSamples;

Console.WriteLine($"Expected P(X > {flow}) = {expectedProb:F6}");
```

---

## Plotting Position Uncertainty

Uncertainty in empirical plotting positions affects comparison with fitted distributions [[4]](#ref4).

### Plotting Position Variance

For the Weibull plotting position $p_i = i/(n+1)$:

```math
\text{Var}(p_i) = \frac{p_i(1-p_i)}{n+2}
```

```cs
int n = data.Length;
Array.Sort(data);

Console.WriteLine("Rank  Value     Plot Pos   SE(p)");
for (int i = 1; i <= n; i++)
{
    double p = (double)i / (n + 1);  // Weibull plotting position
    double seP = Math.Sqrt(p * (1 - p) / (n + 2));
    
    Console.WriteLine($"{i,4}  {data[i-1],8:F0}  {p,8:F4}   {seP:F4}");
}
```

### Confidence Bands for Frequency Curve

```cs
// Generate confidence bands at each plotting position
double[] pValues = new double[100];
double[] qFitted = new double[100];
double[] qLower = new double[100];
double[] qUpper = new double[100];

for (int i = 0; i < 100; i++)
{
    pValues[i] = 0.01 + 0.98 * i / 99.0;  // 0.01 to 0.99
    
    qFitted[i] = gev.InverseCDF(pValues[i]);
    double se = gev.QuantileStandardError(pValues[i], n, method);
    
    qLower[i] = qFitted[i] - 1.96 * se;
    qUpper[i] = qFitted[i] + 1.96 * se;
}
```

---

## Regional Uncertainty

For regional frequency analysis, account for index flood uncertainty [[5]](#ref5).

### Index Flood Method

```cs
// At-site mean (index flood)
double indexFlood = Statistics.Mean(data);
double seIndex = Statistics.StandardDeviation(data) / Math.Sqrt(data.Length);

// Regional growth curve (from pooled regional data)
var regionalGEV = new GeneralizedExtremeValue(1.0, 0.2, -0.1);

// At-site quantile estimate
double growthFactor100 = regionalGEV.InverseCDF(0.99);
double q100 = indexFlood * growthFactor100;

// Combined uncertainty (approximately)
double seGrowth = /* from regional analysis */;
double seQ100 = q100 * Math.Sqrt(
    Math.Pow(seIndex / indexFlood, 2) + 
    Math.Pow(seGrowth / growthFactor100, 2)
);

Console.WriteLine($"Q100 = {q100:F0} ± {1.96 * seQ100:F0}");
```

---

## Monte Carlo Uncertainty Propagation

For complex risk calculations, propagate uncertainty through Monte Carlo simulation:

```cs
using Numerics.Sampling.MCMC;

// Get posterior parameter samples from MCMC
var sampler = new DEMCz(priors, logLikelihood);
sampler.Sample();
var posteriorSamples = sampler.Output;

// Propagate through risk calculation
int nSim = 10000;
double[] riskEstimates = new double[nSim];

for (int i = 0; i < nSim; i++)
{
    // Sample parameters from posterior
    double[] theta = posteriorSamples.GetSample(i % posteriorSamples.Length);
    
    // Create distribution with sampled parameters
    var dist = new GeneralizedExtremeValue(theta[0], theta[1], theta[2]);
    
    // Compute risk metric
    riskEstimates[i] = ComputeRisk(dist);
}

// Summarize uncertainty in risk
double meanRisk = Statistics.Mean(riskEstimates);
double ciLower = Statistics.Percentile(riskEstimates, 0.05);
double ciUpper = Statistics.Percentile(riskEstimates, 0.95);

Console.WriteLine($"Risk = {meanRisk:F4} (90% CI: [{ciLower:F4}, {ciUpper:F4}])");
```

---

## Summary of Uncertainty Methods

| Method | Best For | Assumptions |
|--------|----------|-------------|
| Analytical SE | Quick estimates, Normal-based CI | Large sample, MLE |
| Parametric Bootstrap | Model-based uncertainty | Correct model specification |
| Non-parametric Bootstrap | Distribution-free | IID samples |
| MCMC | Full posterior, complex models | Proper priors |
| Regional Analysis | Short records | Regional homogeneity |

---

## References

<a id="ref1">[1]</a> Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.

<a id="ref2">[2]</a> Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the Bootstrap*. Chapman & Hall.

<a id="ref3">[3]</a> Stedinger, J. R. (1983). Confidence intervals for design events. *Journal of Hydraulic Engineering*, 109(1), 13-27.

<a id="ref4">[4]</a> Hirsch, R. M., & Stedinger, J. R. (1987). Plotting positions for historical floods and their precision. *Water Resources Research*, 23(4), 715-727.

<a id="ref5">[5]</a> Hosking, J. R. M., & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press.
