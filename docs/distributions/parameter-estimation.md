# Parameter Estimation

[← Previous: Univariate Distributions](univariate.md) | [Back to Index](../index.md) | [Next: Uncertainty Analysis →](uncertainty-analysis.md)

Parameter estimation is the process of fitting probability distributions to observed data. The ***Numerics*** library provides multiple estimation methods, each with different statistical properties and appropriate use cases. The library supports Method of Moments (MOM), L-Moments (L-MOM), Maximum Likelihood Estimation (MLE), and Method of Percentiles.

## Estimation Methods

The `ParameterEstimationMethod` enum defines the available methods:

```cs
public enum ParameterEstimationMethod
{
    MaximumLikelihood,        // Maximum likelihood estimation
    MethodOfMoments,          // Product moments
    MethodOfLinearMoments,    // L-moments
    MethodOfPercentiles       // Least squares / percentiles
}
```

### Comparison of Methods

| Method | Advantages | Disadvantages | Best For |
|--------|-----------|---------------|----------|
| **L-Moments** | Robust to outliers, unbiased for small samples, efficient | Computationally intensive | **Hydrological data, small samples** |
| **Maximum Likelihood** | Asymptotically efficient, optimal for large samples | Sensitive to outliers, can fail to converge | Large samples, well-behaved data |
| **Method of Moments** | Simple, fast | Inefficient, biased for small samples | Quick estimates, stable parameters |
| **Method of Percentiles** | Intuitive, robust | Less efficient | Expert judgment, special cases |

**Recommendation for Hydrological Applications:** L-Moments are recommended by USGS [[1]](#1) for flood frequency analysis due to superior performance with small samples and robustness to outliers.

## Using the Estimate() Method

The simplest way to fit a distribution is using the `Estimate()` method:

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

double[] annualPeakFlows = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

// Method 1: Using Estimate() with estimation method
var gev = new GeneralizedExtremeValue();
gev.Estimate(annualPeakFlows, ParameterEstimationMethod.MethodOfLinearMoments);

Console.WriteLine($"GEV Parameters (L-Moments):");
Console.WriteLine($"  Location (ξ): {gev.Xi:F2}");
Console.WriteLine($"  Scale (α): {gev.Alpha:F2}");
Console.WriteLine($"  Shape (κ): {gev.Kappa:F4}");

// Method 2: Using MLE
gev.Estimate(annualPeakFlows, ParameterEstimationMethod.MaximumLikelihood);

Console.WriteLine($"GEV Parameters (MLE):");
Console.WriteLine($"  Location (ξ): {gev.Xi:F2}");
Console.WriteLine($"  Scale (α): {gev.Alpha:F2}");
Console.WriteLine($"  Shape (κ): {gev.Kappa:F4}");
```

The `Estimate()` method automatically:
1. Computes the required moments or likelihood
2. Estimates parameters using the specified method
3. Sets the distribution parameters
4. Validates the parameters

## Manual Parameter Estimation

For more control, use the `ParametersFrom*` methods to compute parameters without setting them:

### L-Moments Approach

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

double[] data = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

// Step 1: Compute L-moments from data
double[] lMoments = Statistics.LinearMoments(data);

Console.WriteLine("Sample L-Moments:");
Console.WriteLine($"  λ₁ (mean): {lMoments[0]:F2}");
Console.WriteLine($"  λ₂ (L-scale): {lMoments[1]:F2}");
Console.WriteLine($"  τ₃ (L-skewness): {lMoments[2]:F4}");
Console.WriteLine($"  τ₄ (L-kurtosis): {lMoments[3]:F4}");

// Step 2: Estimate parameters from L-moments
var gev = new GeneralizedExtremeValue();
double[] parameters = gev.ParametersFromLinearMoments(lMoments);

Console.WriteLine($"\nEstimated Parameters:");
for (int i = 0; i < parameters.Length; i++)
{
    Console.WriteLine($"  Parameter {i}: {parameters[i]:F4}");
}

// Step 3: Set the parameters
gev.SetParameters(parameters);

// Verify: Check theoretical L-moments from fitted distribution
double[] theoreticalLMoments = gev.LinearMomentsFromParameters(parameters);

Console.WriteLine($"\nTheoretical L-Moments from fitted distribution:");
Console.WriteLine($"  λ₁: {theoreticalLMoments[0]:F2}");
Console.WriteLine($"  λ₂: {theoreticalLMoments[1]:F2}");
Console.WriteLine($"  τ₃: {theoreticalLMoments[2]:F4}");
Console.WriteLine($"  τ₄: {theoreticalLMoments[3]:F4}");
```

### Product Moments Approach

```cs
// Step 1: Compute product moments from data
double[] moments = Statistics.ProductMoments(data);

Console.WriteLine("Sample Product Moments:");
Console.WriteLine($"  Mean: {moments[0]:F2}");
Console.WriteLine($"  Std Dev: {moments[1]:F2}");
Console.WriteLine($"  Skewness: {moments[2]:F4}");
Console.WriteLine($"  Kurtosis: {moments[3]:F4}");

// Step 2: Estimate parameters from moments
var normal = new Normal();
double[] normParams = normal.ParametersFromMoments(moments);

Console.WriteLine($"\nNormal Parameters:");
Console.WriteLine($"  μ: {normParams[0]:F2}");
Console.WriteLine($"  σ: {normParams[1]:F2}");

// Step 3: Set parameters
normal.SetParameters(normParams);

// Verify: Check theoretical moments
double[] theoreticalMoments = normal.MomentsFromParameters(normParams);
Console.WriteLine($"\nTheoretical Moments:");
Console.WriteLine($"  Mean: {theoreticalMoments[0]:F2}");
Console.WriteLine($"  Std Dev: {theoreticalMoments[1]:F2}");
Console.WriteLine($"  Skewness: {theoreticalMoments[2]:F4}");
```

## L-Moments (Linear Moments)

L-moments are linear combinations of order statistics that provide robust alternatives to conventional moments [[2]](#2). They are especially valuable for:
- Small sample sizes (n < 50)
- Data with outliers
- Hydrological applications
- Extreme value analysis

### Properties of L-Moments

1. **More robust** than conventional moments - less influenced by outliers
2. **Less biased** for small samples
3. **More efficient** - smaller sampling variance
4. **Bounded** - L-moment ratios are bounded, unlike conventional moments
5. **Nearly unbiased** even for very small samples (n = 10)

### Computing L-Moments

```cs
using Numerics.Data.Statistics;

double[] sample = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9 };

// Compute L-moments
double[] lMoments = Statistics.LinearMoments(sample);

Console.WriteLine("L-Moments:");
Console.WriteLine($"  λ₁ (L-location/mean): {lMoments[0]:F3}");
Console.WriteLine($"  λ₂ (L-scale): {lMoments[1]:F3}");
Console.WriteLine($"  τ₃ (L-skewness): {lMoments[2]:F4}");
Console.WriteLine($"  τ₄ (L-kurtosis): {lMoments[3]:F4}");

// L-moment ratios
double tau3 = lMoments[2];  // L-skewness
double tau4 = lMoments[3];  // L-kurtosis

// Interpret L-skewness
if (Math.Abs(tau3) < 0.1)
    Console.WriteLine("Distribution is approximately symmetric");
else if (tau3 > 0)
    Console.WriteLine("Distribution is right-skewed");
else
    Console.WriteLine("Distribution is left-skewed");
```

### L-Moment Diagrams

L-moment diagrams plot L-skewness (τ₃) vs L-kurtosis (τ₄) to identify appropriate distributions [[2]](#2):

```cs
// Theoretical L-moment ratios for distributions
var distributions = new (string Name, IUnivariateDistribution Dist)[]
{
    ("GEV", new GeneralizedExtremeValue(1000, 200, -0.1)),
    ("Gumbel", new Gumbel(1000, 200)),
    ("Normal", new Normal(1000, 200)),
    ("LP3", new LogPearsonTypeIII(7.0, 0.2, 0.3))
};

Console.WriteLine("Distribution | τ₃ (L-skew) | τ₄ (L-kurt)");
Console.WriteLine("--------------------------------------------");

foreach (var (name, dist) in distributions)
{
    var lMom = dist.LinearMomentsFromParameters(dist.GetParameters);
    Console.WriteLine($"{name,-12} | {lMom[2],11:F4} | {lMom[3],11:F4}");
}

// Compare with sample
double[] sampleLM = Statistics.LinearMoments(sample);
Console.WriteLine($"{"Sample",-12} | {sampleLM[2],11:F4} | {sampleLM[3],11:F4}");
```

## Maximum Likelihood Estimation

MLE finds parameters that maximize the likelihood of observing the data [[3]](#3):

```cs
using Numerics.Distributions;

double[] observations = { 12.5, 15.3, 11.2, 18.7, 14.1, 16.8, 13.4, 17.2 };

// Fit using MLE
var weibull = new Weibull();
weibull.Estimate(observations, ParameterEstimationMethod.MaximumLikelihood);

Console.WriteLine($"Weibull Parameters (MLE):");
Console.WriteLine($"  Scale (λ): {weibull.Lambda:F3}");
Console.WriteLine($"  Shape (κ): {weibull.Kappa:F3}");

// Compute log-likelihood at fitted parameters
double logLikelihood = 0;
foreach (var x in observations)
{
    logLikelihood += weibull.LogPDF(x);
}

Console.WriteLine($"Log-likelihood: {logLikelihood:F4}");
```

### MLE Properties

**Advantages:**
- Asymptotically efficient (minimum variance for large n)
- Invariant under transformation
- Provides likelihood for model comparison (AIC, BIC)

**Disadvantages:**
- Can fail to converge for difficult distributions
- Sensitive to outliers
- Biased for small samples
- Computationally expensive (requires optimization)

### When MLE May Fail

```cs
// MLE can fail with difficult starting values or poor data
var gev = new GeneralizedExtremeValue();

try
{
    gev.Estimate(observations, ParameterEstimationMethod.MaximumLikelihood);
    Console.WriteLine("MLE converged successfully");
}
catch (Exception ex)
{
    Console.WriteLine($"MLE failed: {ex.Message}");
    
    // Fall back to L-moments
    Console.WriteLine("Falling back to L-moments...");
    gev.Estimate(observations, ParameterEstimationMethod.MethodOfLinearMoments);
    Console.WriteLine("L-moments estimation successful");
}
```

## Method of Moments

MOM matches sample moments with theoretical moments:

```cs
double[] data = { 100, 105, 98, 110, 95, 102, 108, 97, 103, 106 };

// Product moments
double[] moments = Statistics.ProductMoments(data);

// Fit Normal distribution
var normal = new Normal();
normal.Estimate(data, ParameterEstimationMethod.MethodOfMoments);

// Or manually
double[] params = normal.ParametersFromMoments(moments);
normal.SetParameters(params);

Console.WriteLine($"Normal Distribution (MOM):");
Console.WriteLine($"  μ = {normal.Mu:F2}");
Console.WriteLine($"  σ = {normal.Sigma:F2}");
Console.WriteLine($"  Sample mean = {moments[0]:F2}");
Console.WriteLine($"  Sample std dev = {moments[1]:F2}");
```

## Distribution-Specific Estimation

### Log-Pearson Type III (USGS Bulletin 17C)

For USGS flood frequency analysis [[1]](#1):

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

double[] annualPeaks = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200, 10500, 19300 };

// USGS recommends L-moments with Expected Moments Algorithm (EMA) adjustments
var lp3 = new LogPearsonTypeIII();
lp3.Estimate(annualPeaks, ParameterEstimationMethod.MethodOfLinearMoments);

Console.WriteLine("LP3 Parameters (Bulletin 17C method):");
Console.WriteLine($"  μ (log-space mean): {lp3.Mu:F4}");
Console.WriteLine($"  σ (log-space std dev): {lp3.Sigma:F4}");
Console.WriteLine($"  γ (log-space skew): {lp3.Gamma:F4}");

// Compute flood frequency curve
var returnPeriods = new int[] { 2, 5, 10, 25, 50, 100, 200, 500 };

Console.WriteLine("\nFlood Frequency Analysis:");
Console.WriteLine("T (years) | AEP     | Discharge");
foreach (var T in returnPeriods)
{
    double aep = 1.0 / T;
    double Q = lp3.InverseCDF(1 - aep);
    Console.WriteLine($"{T,9} | {aep,7:F5} | {Q,10:F0}");
}
```

### Generalized Extreme Value

```cs
// Block maxima approach
double[] blockMaxima = { 125, 153, 112, 187, 141, 168 };

var gev = new GeneralizedExtremeValue();

// Try MLE first
try
{
    gev.Estimate(blockMaxima, ParameterEstimationMethod.MaximumLikelihood);
    Console.WriteLine("Fitted using MLE");
}
catch
{
    // Fall back to L-moments
    gev.Estimate(blockMaxima, ParameterEstimationMethod.MethodOfLinearMoments);
    Console.WriteLine("Fitted using L-moments (MLE failed)");
}

Console.WriteLine($"ξ = {gev.Xi:F2}, α = {gev.Alpha:F2}, κ = {gev.Kappa:F4}");

// Classify GEV type
if (Math.Abs(gev.Kappa) < 0.01)
    Console.WriteLine("Approximately Gumbel (Type I)");
else if (gev.Kappa < 0)
    Console.WriteLine($"Weibull/Type III (bounded above at {gev.Xi - gev.Alpha / gev.Kappa:F1})");
else
    Console.WriteLine("Fréchet/Type II (heavy tail)");
```

### Two-Parameter Distributions

For simpler distributions:

```cs
double[] data = { 10.5, 12.3, 11.8, 15.2, 13.7, 14.1, 16.8, 12.9 };

// Exponential - two parameters (location + scale)
var exponential = new Exponential();
exponential.Estimate(data, ParameterEstimationMethod.MethodOfMoments);
Console.WriteLine($"Exponential ξ = {exponential.Xi:F4}, α = {exponential.Alpha:F4}");

// Log-Normal - two parameters
var lognormal = new LogNormal();
lognormal.Estimate(data, ParameterEstimationMethod.MethodOfMoments);
Console.WriteLine($"LogNormal μ = {lognormal.Mu:F4}, σ = {lognormal.Sigma:F4}");

// Weibull - two parameters (MLE only)
var weibull = new Weibull();
weibull.Estimate(data, ParameterEstimationMethod.MaximumLikelihood);
Console.WriteLine($"Weibull λ = {weibull.Lambda:F4}, κ = {weibull.Kappa:F4}");
```

## Tutorial: Complete Flood Frequency Analysis

This tutorial demonstrates a complete distribution fitting workflow using real streamflow data from the White River near Nora, Indiana. The data and expected results are drawn from published references [[4]](#4) and validated against the R `lmom` package [[2]](#2).

**Data source:** Rao, A. R. & Hamed, K. H. (2000). *Flood Frequency Analysis*. CRC Press, Table 7.1.2.
See also: [`example-data/white-river-nora-floods.csv`](../example-data/white-river-nora-floods.csv)

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

// White River near Nora, Indiana — 62 years of annual peak streamflow (cfs)
// Source: Rao & Hamed (2000), Table 7.1.2
double[] annualPeaks = {
    23200, 2950, 10300, 23200, 4540, 9960, 10800, 26900, 23300, 20400,
    8480, 3150, 9380, 32400, 20800, 11100, 7270, 9600, 14600, 14300,
    22500, 14700, 12700, 9740, 3050, 8830, 12000, 30400, 27000, 15200,
    8040, 11700, 20300, 22700, 30400, 9180, 4870, 14700, 12800, 13700,
    7960, 9830, 12500, 10700, 13200, 14700, 14300, 4050, 14600, 14400,
    19200, 7160, 12100, 8650, 10600, 24500, 14400, 6300, 9560, 15800,
    14300, 28700
};

Console.WriteLine($"Record length: {annualPeaks.Length} years");
Console.WriteLine($"Range: {annualPeaks.Min():F0} - {annualPeaks.Max():F0} cfs");

// Step 1: Compute sample L-moments
// L-moments are more robust than product moments for small to moderate samples.
// Validated against R lmom::samlmu()
double[] lMoments = Statistics.LinearMoments(annualPeaks);

Console.WriteLine($"\nSample L-moments:");
Console.WriteLine($"  λ₁ (L-location): {lMoments[0]:F1}");
Console.WriteLine($"  λ₂ (L-scale):    {lMoments[1]:F1}");
Console.WriteLine($"  τ₃ (L-skewness): {lMoments[2]:F4}");
Console.WriteLine($"  τ₄ (L-kurtosis): {lMoments[3]:F4}");

// Step 2: Fit candidate distributions using L-moments
var candidates = new List<(string Name, IUnivariateDistribution Dist)>
{
    ("LP3",    new LogPearsonTypeIII()),
    ("GEV",    new GeneralizedExtremeValue()),
    ("Gumbel", new Gumbel()),
    ("PIII",   new PearsonTypeIII())
};

foreach (var (name, dist) in candidates)
{
    dist.Estimate(annualPeaks, ParameterEstimationMethod.MethodOfLinearMoments);
    Console.WriteLine($"\n{name} fitted parameters:");
    var pNames = dist.ParameterNamesShortForm;
    var pValues = dist.GetParameters;
    for (int i = 0; i < dist.NumberOfParameters; i++)
        Console.WriteLine($"  {pNames[i]} = {pValues[i]:F4}");
}

// Step 3: Compare estimation methods for GEV
// Textbook (Rao & Hamed, Example 7.1.1, p. 218) provides MOM results for comparison.
Console.WriteLine("\nGEV: Comparing estimation methods:");
foreach (var method in new[] {
    ParameterEstimationMethod.MethodOfLinearMoments,
    ParameterEstimationMethod.MethodOfMoments,
    ParameterEstimationMethod.MaximumLikelihood })
{
    var gev = new GeneralizedExtremeValue();
    try
    {
        gev.Estimate(annualPeaks, method);
        Console.WriteLine($"  {method}: ξ={gev.Xi:F1}, α={gev.Alpha:F1}, κ={gev.Kappa:F4}");
    }
    catch (Exception ex)
    {
        Console.WriteLine($"  {method}: Failed — {ex.Message}");
    }
}

// Step 4: Compute flood frequency curve
Console.WriteLine("\nFlood Frequency Curve (LP3, L-moments):");
Console.WriteLine("T (years) |    AEP   | Discharge (cfs)");
Console.WriteLine("----------|----------|----------------");

var lp3 = (LogPearsonTypeIII)candidates[0].Dist;
foreach (int T in new[] { 2, 5, 10, 25, 50, 100, 200, 500 })
{
    double aep = 1.0 / T;
    double Q = lp3.InverseCDF(1 - aep);
    Console.WriteLine($"{T,9} | {aep,8:F5} | {Q,14:F0}");
}

// Step 5: Compare candidate distributions at key return periods
Console.WriteLine("\nQuantile Comparison (cfs):");
Console.Write("   AEP   ");
foreach (var c in candidates) Console.Write($"| {c.Name,8} ");
Console.WriteLine();
Console.WriteLine(new string('-', 55));

foreach (double p in new[] { 0.5, 0.9, 0.98, 0.99, 0.998 })
{
    Console.Write($"  {1 - p:F3}  ");
    foreach (var c in candidates)
        Console.Write($"| {c.Dist.InverseCDF(p),8:F0} ");
    Console.WriteLine();
}

// Step 6: Model selection using AIC
Console.WriteLine("\nModel Selection (AIC):");
foreach (var (name, dist) in candidates)
{
    double logLik = dist.LogLikelihood(annualPeaks);
    double aic = GoodnessOfFit.AIC(dist.NumberOfParameters, logLik);
    Console.WriteLine($"  {name,-8}: AIC = {aic:F1}");
}
```

## Estimation with Censored Data

For data with detection limits or censoring:

```cs
// Low flows below detection limit (left-censored)
double detectionLimit = 5.0;
var observed = data.Where(x => x >= detectionLimit).ToArray();
int nCensored = data.Length - observed.Length;

Console.WriteLine($"Observed: {observed.Length}, Censored: {nCensored}");

// Fit using only observed values
var lognormal = new LogNormal();
lognormal.Estimate(observed, ParameterEstimationMethod.MethodOfMoments);

// Note: This is a simple approach. For formal censored data analysis,
// use MLE with censored likelihood (requires custom implementation)
```

## Tips and Best Practices

### 1. Sample Size Requirements

```cs
// Rule of thumb: n > 10 * number of parameters for MLE
int minSamples = 10 * dist.NumberOfParameters;
if (data.Length < minSamples)
{
    Console.WriteLine($"Warning: Sample size ({data.Length}) below recommended minimum ({minSamples})");
    Console.WriteLine("Consider using L-moments for improved small-sample efficiency");
}
```

### 2. Checking Parameter Validity

```cs
var gev = new GeneralizedExtremeValue();
gev.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);

if (!gev.ParametersValid)
{
    Console.WriteLine("Warning: Invalid parameters estimated");
    Console.WriteLine("Try different estimation method or distribution");
}
```

### 3. Comparing Estimation Methods

```cs
var methods = new[]
{
    ParameterEstimationMethod.MethodOfLinearMoments,
    ParameterEstimationMethod.MethodOfMoments,
    ParameterEstimationMethod.MaximumLikelihood
};

foreach (var method in methods)
{
    var dist = new GeneralizedExtremeValue();
    try
    {
        dist.Estimate(data, method);
        Console.WriteLine($"{method}: κ = {dist.Kappa:F4}");
    }
    catch (Exception ex)
    {
        Console.WriteLine($"{method}: Failed - {ex.Message}");
    }
}
```

### 4. Outlier Detection

```cs
// Identify potential outliers before estimation
double[] sorted = data.OrderBy(x => x).ToArray();
double Q1 = Statistics.Quantile(sorted, 0.25);
double Q3 = Statistics.Quantile(sorted, 0.75);
double IQR = Q3 - Q1;

var outliers = data.Where(x => x < Q1 - 1.5 * IQR || x > Q3 + 1.5 * IQR).ToArray();

if (outliers.Length > 0)
{
    Console.WriteLine($"Potential outliers detected: {outliers.Length}");
    Console.WriteLine("Consider using L-moments (more robust to outliers)");
}
```

### 5. Historical Information

When historical data exists outside the systematic record:

```cs
// Combine systematic record with historical peaks
double[] systematicRecord = { 12.5, 15.3, 11.2, 18.7, 14.1 }; // Recent, complete
double[] historicalPeaks = { 22.3 };  // Known historical floods

// This is a simplified approach
// For formal analysis, use Expected Moments Algorithm (EMA)
var combined = systematicRecord.Concat(historicalPeaks).ToArray();

var lp3 = new LogPearsonTypeIII();
lp3.Estimate(combined, ParameterEstimationMethod.MethodOfLinearMoments);

Console.WriteLine("Fitted with historical information included");
```

## Common Pitfalls

1. **Using MOM for small samples** - Use L-moments instead
2. **Ignoring convergence failures in MLE** - Always check and have fallback
3. **Not checking parameter validity** - Always validate after estimation
4. **Wrong distribution family** - Use L-moment diagrams for selection
5. **Ignoring outliers** - L-moments are more robust than MOM or MLE
6. **Insufficient sample size** - Need at least 10× the number of parameters (e.g., 30 samples for a 3-parameter distribution)

---

## References

<a id="1">[1]</a> England, J. F., Jr., Cohn, T. A., Faber, B. A., Stedinger, J. R., Thomas, W. O., Jr., Veilleux, A. G., Kiang, J. E., & Mason, R. R., Jr. (2018). *Guidelines for Determining Flood Flow Frequency—Bulletin 17C* (ver. 1.1, May 2019): U.S. Geological Survey Techniques and Methods, book 4, chap. B5, 148 p.

<a id="2">[2]</a> Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society: Series B (Methodological)*, 52(1), 105-124.

<a id="3">[3]</a> Mood, A. M., Graybill, F. A., & Boes, D. C. (1974). *Introduction to the Theory of Statistics* (3rd ed.). McGraw-Hill.

<a id="4">[4]</a> Rao, A. R., & Hamed, K. H. (2000). *Flood Frequency Analysis*. CRC Press.

---

[← Previous: Univariate Distributions](univariate.md) | [Back to Index](../index.md) | [Next: Uncertainty Analysis →](uncertainty-analysis.md)
