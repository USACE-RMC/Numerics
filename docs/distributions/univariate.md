# Univariate Distributions

[← Back to Index](../index.md) | [Next: Parameter Estimation →](parameter-estimation.md)

The ***Numerics*** library provides over 40 univariate probability distributions for statistical analysis, risk assessment, and uncertainty quantification. All distributions implement a common interface with consistent methods for computing probability density functions (PDF), cumulative distribution functions (CDF), quantiles, and statistical moments.

## Available Distributions

### Continuous Distributions

| Distribution | Parameters | Typical Applications |
|--------------|-----------|---------------------|
| **Normal** | μ (mean), σ (std dev) | General purpose, natural phenomena |
| **Log-Normal** | μ, σ | Right-skewed data, multiplicative processes |
| **Uniform** | a (min), b (max) | Maximum entropy, prior distributions |
| **Exponential** | λ (rate) | Time between events, survival analysis |
| **Gamma** | α (shape), β (scale) | Waiting times, rainfall |
| **Beta** | α, β | Probabilities, proportions, [0,1] bounded |
| **Weibull** | α (scale), β (shape) | Failure times, wind speed |
| **Gumbel** | ξ (location), α (scale) | Extreme values (maxima) |
| **Generalized Extreme Value (GEV)** | ξ, α, κ (shape) | Block maxima, floods, earthquakes |
| **Generalized Pareto (GP)** | ξ, α, κ | Exceedances over threshold |
| **Log-Pearson Type III (LP3)** | μ, σ, γ (skew) | USGS flood frequency analysis |
| **Pearson Type III (P3)** | μ, σ, γ | Flood frequency, rainfall |
| **Kappa-Four (K4)** | ξ, α, κ, h | Flexible 4-parameter family |
| **Generalized Logistic (GLO)** | ξ, α, κ | Growth models, extreme values |
| **Generalized Normal (GNO)** | ξ, α, κ | Flexible alternative to GEV |
| **Generalized Beta (GB)** | a, b, α, β | [a,b] bounded with flexibility |
| **Triangular** | a, b, c (mode) | Simple uncertainty modeling |
| **Rayleigh** | σ | Wind speed, wave height |
| **Cauchy** | x₀ (location), γ (scale) | Heavy-tailed phenomena |
| **Logistic** | μ, s | Growth processes, neural networks |
| **Student's t** | ν (degrees of freedom) | Heavy-tailed alternative to Normal |
| **Noncentral t** | ν, δ (noncentrality) | Power analysis, hypothesis testing |
| **Chi-Squared** | k (degrees of freedom) | Variance estimation, goodness-of-fit |
| **Inverse Gamma** | α, β | Bayesian priors for variance |
| **Inverse Chi-Squared** | ν | Bayesian inference |
| **Pareto** | xₘ (scale), α (shape) | Income distributions, city sizes |
| **PERT** | a, b, c | Project management, expert judgment |
| **PERT Percentile** | P₁₀, P₅₀, P₉₀ | Expert percentile elicitation |
| **PERT Percentile Z** | Similar to PERT Percentile | Alternative parametrization |
| **Truncated Normal** | μ, σ, a, b | Bounded normal distributions |
| **Truncated Distribution** | Any distribution + bounds | Bounded versions of distributions |
| **Mixture** | Multiple distributions | Multi-modal data |
| **Empirical** | Sample data | Non-parametric, data-driven |
| **Kernel Density** | Sample data, bandwidth | Smooth non-parametric estimation |
| **Deterministic** | Single value | Point estimates, constants |
| **Competing Risks** | Multiple distributions | Failure analysis with multiple causes |

### Discrete Distributions

| Distribution | Parameters | Typical Applications |
|--------------|-----------|---------------------|
| **Bernoulli** | p (success probability) | Binary outcomes |
| **Binomial** | n (trials), p (success prob) | Number of successes in n trials |
| **Poisson** | λ (rate) | Count data, rare events |
| **Geometric** | p | Number of trials until first success |
| **Uniform Discrete** | a, b | Discrete uniform outcomes |

## Common Interface

All univariate distributions in ***Numerics*** implement the `IUnivariateDistribution` interface, providing:

### Statistical Properties
```cs
double Mean              // E[X]
double Median            // 50th percentile
double Mode              // Most likely value
double Variance          // Var(X)
double StandardDeviation // √Var(X)
double Skewness          // Measure of asymmetry
double Kurtosis          // Measure of tail heaviness
double Minimum           // Support lower bound
double Maximum           // Support upper bound
```

### Probability Functions
```cs
double PDF(double x)         // Probability density (or mass for discrete)
double CDF(double x)         // P(X ≤ x)
double InverseCDF(double p)  // Quantile function (inverse CDF)
double CCDF(double x)        // P(X > x) = 1 - CDF(x)
double HF(double x)          // Hazard function
double LogPDF(double x)      // ln(PDF(x))
double LogCDF(double x)      // ln(CDF(x))
double LogCCDF(double x)     // ln(CCDF(x))
```

### Random Generation
```cs
double[] GenerateRandomValues(int sampleSize, int seed = -1)
```

## Creating Distributions

### Method 1: Direct Construction with Parameters

```cs
using Numerics.Distributions;

// Normal distribution: N(100, 15)
var normal = new Normal(mu: 100, sigma: 15);

// Generalized Extreme Value: GEV(1000, 200, -0.1)
var gev = new GeneralizedExtremeValue(xi: 1000, alpha: 200, kappa: -0.1);

// Log-Normal distribution
var lognormal = new LogNormal(mu: 4.5, sigma: 0.5);

// Gamma distribution
var gamma = new GammaDistribution(alpha: 5, beta: 2);
```

### Method 2: Using SetParameters

```cs
// Create with default parameters, then set
var weibull = new Weibull();
weibull.SetParameters(new double[] { 50, 2.5 }); // alpha=50, beta=2.5

// Or use named parameters
weibull.SetParameters(alpha: 50, beta: 2.5);
```

### Method 3: From Parameter Array

Useful when parameters are computed:

```cs
double[] gevParams = SomeEstimationFunction(data);
var gev = new GeneralizedExtremeValue();
gev.SetParameters(gevParams);
```

## Using Distributions

### Basic Probability Calculations

```cs
using Numerics.Distributions;

var normal = new Normal(100, 15);

// Probability density at x = 110
double pdf = normal.PDF(110);  // f(110)
Console.WriteLine($"PDF at 110: {pdf:F6}");

// Cumulative probability P(X ≤ 110)
double cdf = normal.CDF(110);
Console.WriteLine($"P(X ≤ 110) = {cdf:F4}");  // 0.7475

// Exceedance probability P(X > 110)
double ccdf = normal.CCDF(110);  // or 1 - cdf
Console.WriteLine($"P(X > 110) = {ccdf:F4}");  // 0.2525

// Find quantile: what value corresponds to 95th percentile?
double q95 = normal.InverseCDF(0.95);
Console.WriteLine($"95th percentile: {q95:F2}");  // 124.67
```

### Statistical Properties

```cs
var gev = new GeneralizedExtremeValue(xi: 1000, alpha: 200, kappa: -0.1);

Console.WriteLine($"Mean: {gev.Mean:F2}");
Console.WriteLine($"Std Dev: {gev.StandardDeviation:F2}");
Console.WriteLine($"Skewness: {gev.Skewness:F3}");
Console.WriteLine($"Kurtosis: {gev.Kurtosis:F3}");
Console.WriteLine($"Median: {gev.Median:F2}");
Console.WriteLine($"Mode: {gev.Mode:F2}");
```

### Hazard Function

The hazard function describes instantaneous failure rate:

```cs
var weibull = new Weibull(alpha: 100, beta: 2.5);

// Hazard at time t=50
double hazard = weibull.HF(50);
Console.WriteLine($"Hazard rate at t=50: {hazard:F6}");

// For Weibull, hazard increases with time when β > 1 (wear-out)
```

### Log-Space Calculations

For numerical stability with very small probabilities:

```cs
var normal = new Normal(0, 1);

// Regular CDF can underflow for extreme values
double x = -10;
double logCDF = normal.LogCDF(x);  // ln(CDF(x))
double cdf = Math.Exp(logCDF);

Console.WriteLine($"CDF(-10) = {cdf:E10}");
Console.WriteLine($"Log-CDF(-10) = {logCDF:F4}");
```

## Hydrological Distributions

### Log-Pearson Type III (LP3)

The LP3 distribution is the standard for USGS flood frequency analysis [[1]](#1):

```cs
using Numerics.Distributions;

double[] annualPeakFlows = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

// Fit LP3 using L-Moments (recommended for hydrologic data)
var lp3 = new LogPearsonTypeIII();
lp3.Estimate(annualPeakFlows, ParameterEstimationMethod.MethodOfLinearMoments);

// Or explicitly with ParametersFromLinearMoments
var lMoments = Statistics.LinearMoments(annualPeakFlows);
lp3.SetParameters(lp3.ParametersFromLinearMoments(lMoments));

Console.WriteLine($"LP3 Parameters:");
Console.WriteLine($"  μ: {lp3.Mu:F3}");
Console.WriteLine($"  σ: {lp3.Sigma:F3}");
Console.WriteLine($"  γ: {lp3.Gamma:F3}");

// Compute flood quantiles
double q100 = lp3.InverseCDF(0.99);  // 100-year flood (1% annual exceedance)
double q500 = lp3.InverseCDF(0.998); // 500-year flood
Console.WriteLine($"100-year flood: {q100:F0} cfs");
Console.WriteLine($"500-year flood: {q500:F0} cfs");
```

### Generalized Extreme Value (GEV)

GEV is widely used for extreme value analysis [[2]](#2):

```cs
// Annual maximum flood data
double[] annualMaxima = { 12500, 15300, 11200, 18700, 14100 };

var gev = new GeneralizedExtremeValue();
gev.Estimate(annualMaxima, ParameterEstimationMethod.MethodOfLinearMoments);

Console.WriteLine($"GEV Parameters:");
Console.WriteLine($"  Location (ξ): {gev.Xi:F2}");
Console.WriteLine($"  Scale (α): {gev.Alpha:F2}");
Console.WriteLine($"  Shape (κ): {gev.Kappa:F4}");

// Interpret shape parameter
if (gev.Kappa < 0)
    Console.WriteLine("  Type III (Weibull) - bounded upper tail");
else if (gev.Kappa > 0)
    Console.WriteLine("  Type II (Fréchet) - heavy upper tail");
else
    Console.WriteLine("  Type I (Gumbel) - exponential tail");
```

### Generalized Pareto Distribution (GPD)

For peaks-over-threshold analysis [[3]](#3):

```cs
// Values exceeding a threshold
double threshold = 10000;
var exceedances = annualPeakFlows.Where(x => x > threshold).Select(x => x - threshold).ToArray();

var gpd = new GeneralizedPareto();
gpd.Estimate(exceedances, ParameterEstimationMethod.MethodOfLinearMoments);

// Adjust location parameter for threshold
gpd.SetParameters(threshold, gpd.Alpha, gpd.Kappa);

Console.WriteLine($"GPD for exceedances over {threshold}:");
Console.WriteLine($"  ξ: {gpd.Xi:F2}");
Console.WriteLine($"  α: {gpd.Alpha:F2}");
Console.WriteLine($"  κ: {gpd.Kappa:F4}");
```

## Special Distribution Features

### Truncated Distributions

Create truncated versions of any distribution:

```cs
// Normal truncated to [0, 100]
var truncNormal = new TruncatedNormal(mu: 50, sigma: 15, lowerBound: 0, upperBound: 100);

// Or truncate any distribution
var normal = new Normal(50, 15);
var truncated = new TruncatedDistribution(normal, 0, 100);

double mean = truncated.Mean;  // Different from untruncated mean
```

### Mixture Distributions

Model multi-modal data with mixture distributions:

```cs
// Mixture of two normals (bimodal)
var component1 = new Normal(100, 10);
var component2 = new Normal(150, 15);
var weights = new double[] { 0.6, 0.4 };  // 60% from first, 40% from second

var mixture = new Mixture(new IUnivariateDistribution[] { component1, component2 }, weights);

// PDF will show two peaks
double pdf = mixture.PDF(125);  // Valley between modes
```

### Empirical Distribution

Non-parametric distribution from data:

```cs
double[] observations = { 12.5, 15.3, 11.2, 18.7, 14.1, 16.8, 13.4, 17.2 };

var empirical = new EmpiricalDistribution(observations);

// Uses linear interpolation for quantiles
double median = empirical.InverseCDF(0.5);
double q90 = empirical.InverseCDF(0.9);

Console.WriteLine($"Empirical median: {median:F2}");
Console.WriteLine($"Empirical 90th percentile: {q90:F2}");
```

### Kernel Density Estimation

Smooth non-parametric density estimation:

```cs
var kde = new KernelDensity(observations, bandwidth: 1.5);

// Smooth PDF
double density = kde.PDF(15.0);

// KDE-based CDF and quantiles
double cdf = kde.CDF(15.0);
double quantile = kde.InverseCDF(0.75);
```

### PERT Distributions

For expert judgment and project management:

```cs
// PERT from minimum, most likely, maximum
var pert = new Pert(min: 10, mode: 15, max: 25);

// PERT from percentile judgments
var pertPercentile = new PertPercentile(p10: 12, p50: 15, p90: 22);

// Use for duration or cost uncertainty
double expectedDuration = pert.Mean;
double variance = pert.Variance;
```

## Random Number Generation

All distributions can generate random samples:

```cs
var normal = new Normal(100, 15);

// Generate 1000 random values
double[] samples = normal.GenerateRandomValues(sampleSize: 1000, seed: 12345);

Console.WriteLine($"Sample mean: {samples.Average():F2}");
Console.WriteLine($"Sample std dev: {Statistics.StandardDeviation(samples):F2}");

// Use -1 or 0 for seed to use system clock
double[] randomSamples = normal.GenerateRandomValues(1000, seed: -1);
```

## Practical Examples

### Example 1: Computing Return Periods

```cs
// Fit distribution to annual maximum flood data
double[] annualMaxFlows = { 12500, 15300, 11200, 18700, 14100, 16800 };

var gev = new GeneralizedExtremeValue();
gev.Estimate(annualMaxFlows, ParameterEstimationMethod.MethodOfLinearMoments);

// Compute floods for different return periods
var returnPeriods = new int[] { 2, 5, 10, 25, 50, 100, 200, 500 };

Console.WriteLine("Return Period Analysis:");
Console.WriteLine("Return Period | Annual Exceedance Prob | Flood Magnitude");
Console.WriteLine("-----------------------------------------------------------");

foreach (var T in returnPeriods)
{
    double aep = 1.0 / T;  // Annual exceedance probability
    double nep = 1.0 - aep; // Non-exceedance probability
    double flood = gev.InverseCDF(nep);
    
    Console.WriteLine($"{T,13} | {aep,22:F6} | {flood,15:F0}");
}
```

### Example 2: Probability of Exceedance

```cs
var lp3 = new LogPearsonTypeIII(mu: 10.2, sigma: 0.3, gamma: 0.4);

// What's the probability a flood exceeds 50,000 cfs?
double threshold = 50000;
double exceedanceProb = lp3.CCDF(threshold);
double returnPeriod = 1.0 / exceedanceProb;

Console.WriteLine($"Probability of exceeding {threshold:N0} cfs: {exceedanceProb:F6}");
Console.WriteLine($"Equivalent return period: {returnPeriod:F1} years");
```

### Example 3: Comparing Distributions

```cs
double[] data = { 12.5, 15.3, 11.2, 18.7, 14.1, 16.8, 13.4, 17.2, 10.5, 19.3 };

// Fit multiple distributions
var normal = new Normal();
normal.Estimate(data, ParameterEstimationMethod.MethodOfMoments);

var lognormal = new LogNormal();
lognormal.Estimate(data, ParameterEstimationMethod.MethodOfMoments);

var gev = new GeneralizedExtremeValue();
gev.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);

// Compare at various quantiles
var probs = new double[] { 0.5, 0.9, 0.95, 0.99 };

Console.WriteLine("Quantile Comparison:");
Console.WriteLine("Probability | Normal | Log-Normal | GEV");
Console.WriteLine("----------------------------------------------");

foreach (var p in probs)
{
    Console.WriteLine($"{p,11:F2} | {normal.InverseCDF(p),6:F1} | {lognormal.InverseCDF(p),10:F1} | {gev.InverseCDF(p),3:F1}");
}
```

### Example 4: Reliability Analysis

```cs
// Component with Weibull failure time distribution
var weibull = new Weibull(alpha: 1000, beta: 2.5); // hours

// Reliability at time t (probability of survival)
double t = 500;  // hours
double reliability = weibull.CCDF(t);
double failureProb = weibull.CDF(t);

Console.WriteLine($"At t = {t} hours:");
Console.WriteLine($"  Reliability: {reliability:F4}");
Console.WriteLine($"  Failure probability: {failureProb:F4}");
Console.WriteLine($"  Hazard rate: {weibull.HF(t):E3}");

// Mean time to failure
Console.WriteLine($"  MTTF: {weibull.Mean:F1} hours");
```

### Example 5: Risk Assessment

```cs
// Annual probability of dam failure
var failureProb = new Beta(alpha: 2, beta: 1998); // ~0.001

// Generate scenarios
double[] scenarios = failureProb.GenerateRandomValues(10000, seed: 12345);

// Estimate risk metrics
Console.WriteLine($"Expected annual failure probability: {failureProb.Mean:E4}");
Console.WriteLine($"95th percentile: {failureProb.InverseCDF(0.95):E4}");
Console.WriteLine($"Scenarios > 0.002: {scenarios.Count(x => x > 0.002)} / 10000");
```

## Distribution Selection Guidelines

| Data Characteristics | Recommended Distribution(s) |
|---------------------|----------------------------|
| Symmetric, unbounded | Normal, Student's t (heavy tails) |
| Right-skewed, positive | Log-Normal, Gamma, Weibull |
| Left-skewed | Beta, Generalized Beta |
| Heavy tails | Student's t, Cauchy, Pareto |
| Bounded [a,b] | Uniform, Beta, Triangular, PERT |
| Extreme values (maxima) | GEV, Gumbel, Weibull |
| Extreme values (minima) | GEV (negative), Weibull (reversed) |
| Threshold exceedances | Generalized Pareto |
| Flood frequency | LP3, GEV, Pearson Type III |
| Failure/survival times | Weibull, Exponential, Gamma |
| Count data | Poisson, Binomial |
| Expert judgment | PERT, PERT Percentile, Triangular |
| Non-parametric | Empirical, Kernel Density |

## Parameter Bounds and Validation

All distributions validate parameters:

```cs
var gev = new GeneralizedExtremeValue();

// Check if parameters are valid
var params = new double[] { 1000, 200, 0.3 };
var exception = gev.ValidateParameters(params, throwException: false);

if (exception == null)
{
    gev.SetParameters(params);
    Console.WriteLine("Parameters are valid");
}
else
{
    Console.WriteLine($"Invalid parameters: {exception.Message}");
}

// Get parameter bounds
double[] minParams = gev.MinimumOfParameters;
double[] maxParams = gev.MaximumOfParameters;

Console.WriteLine("Parameter bounds:");
for (int i = 0; i < gev.NumberOfParameters; i++)
{
    Console.WriteLine($"  Param {i}: [{minParams[i]}, {maxParams[i]}]");
}
```

## Distribution Information

```cs
var normal = new Normal(100, 15);

// Display information
Console.WriteLine($"Distribution: {normal.DisplayName}");
Console.WriteLine($"Short name: {normal.ShortDisplayName}");
Console.WriteLine($"Type: {normal.Type}");
Console.WriteLine($"Parameters: {normal.DisplayLabel}");
Console.WriteLine($"Number of parameters: {normal.NumberOfParameters}");

// Parameter names
string[] paramNames = normal.ParameterNamesShortForm;
double[] paramValues = normal.GetParameters;

for (int i = 0; i < normal.NumberOfParameters; i++)
{
    Console.WriteLine($"  {paramNames[i]} = {paramValues[i]:F3}");
}
```

---

## References

<a id="1">[1]</a> Bulletin 17C: Guidelines for Determining Flood Flow Frequency. (2017). U.S. Geological Survey Techniques and Methods, Book 4, Chapter B5.

<a id="2">[2]</a> Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.

<a id="3">[3]</a> Hosking, J. R. M., & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press.

---

[← Back to Index](../index.md) | [Next: Parameter Estimation →](parameter-estimation.md)
