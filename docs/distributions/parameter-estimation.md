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

**Recommendation for Hydrological Applications:** L-moments are recommended by USGS [[1]](#1) for flood frequency analysis due to superior performance with small samples and robustness to outliers.

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

Console.WriteLine("Sample L-moments:");
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

### Mathematical Formulation

L-moments are defined through probability-weighted moments (PWMs). For a random variable $X$ with CDF $F(x)$, the probability-weighted moments are:

```math
\beta_r = E\left[X \cdot F(X)^r\right] = \int_0^1 x(F) \cdot F^r \, dF, \quad r = 0, 1, 2, \ldots
```

The first four L-moments are linear combinations of the PWMs:

```math
\lambda_1 = \beta_0
```

```math
\lambda_2 = 2\beta_1 - \beta_0
```

```math
\lambda_3 = 6\beta_2 - 6\beta_1 + \beta_0
```

```math
\lambda_4 = 20\beta_3 - 30\beta_2 + 12\beta_1 - \beta_0
```

The L-moment ratios, which are dimensionless and bounded, are defined as:

```math
\tau = \frac{\lambda_2}{\lambda_1} \quad \text{(L-CV)}, \qquad \tau_3 = \frac{\lambda_3}{\lambda_2} \quad \text{(L-skewness)}, \qquad \tau_4 = \frac{\lambda_4}{\lambda_2} \quad \text{(L-kurtosis)}
```

L-skewness is bounded in $[-1, 1]$ and L-kurtosis in $[\frac{1}{4}(5\tau_3^2 - 1),\; 1]$, unlike conventional skewness and kurtosis which are unbounded. This boundedness makes L-moment ratios more interpretable and stable.

**Sample estimation.** Given a sorted sample $x_{1:n} \leq x_{2:n} \leq \cdots \leq x_{n:n}$, the unbiased sample PWM estimators are:

```math
b_r = \frac{1}{n}\sum_{j=r+1}^{n} \frac{\binom{j-1}{r}}{\binom{n-1}{r}} \, x_{j:n}, \quad r = 0, 1, 2, \ldots
```

The `Statistics.LinearMoments()` method computes these sample PWMs and returns the array $[\lambda_1,\; \lambda_2,\; \tau_3,\; \tau_4]$.

**Why L-moments are preferred for small samples.** Conventional moments involve powers of deviations from the mean, so a single extreme observation can dominate the skewness or kurtosis estimate. L-moments use only linear combinations of order statistics, which makes them far more robust to outliers and nearly unbiased even for samples as small as $n = 10$. For hydrological applications where sample sizes are often 30--60 years of annual data, this robustness is critical.

**L-moment ratio diagrams.** Plotting sample L-skewness ($\tau_3$) against L-kurtosis ($\tau_4$) and comparing to the theoretical curves of candidate distributions is a powerful tool for distribution identification. Each distribution family traces a distinct curve (or point) in L-moment ratio space, making visual comparison straightforward [[2]](#2).

### Properties of L-Moments

1. **More robust** than conventional moments -- less influenced by outliers
2. **Less biased** for small samples
3. **More efficient** -- smaller sampling variance
4. **Bounded** -- L-moment ratios are bounded, unlike conventional moments
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

Maximum Likelihood Estimation (MLE) finds the parameter values that make the observed data most probable under the assumed model [[3]](#3).

### Mathematical Formulation

Given independent observations $x_1, x_2, \ldots, x_n$ from a distribution with PDF $f(x|\boldsymbol{\theta})$, the likelihood function is the joint probability of the data viewed as a function of the parameters:

```math
L(\boldsymbol{\theta} \,|\, \mathbf{x}) = \prod_{i=1}^{n} f(x_i \,|\, \boldsymbol{\theta})
```

Because products are numerically unstable, optimization is performed on the log-likelihood:

```math
\ell(\boldsymbol{\theta}) = \sum_{i=1}^{n} \log f(x_i \,|\, \boldsymbol{\theta})
```

The MLE is the parameter vector that maximizes the log-likelihood:

```math
\hat{\boldsymbol{\theta}}_{\text{MLE}} = \underset{\boldsymbol{\theta}}{\text{argmax}} \; \ell(\boldsymbol{\theta})
```

For some distributions (e.g., Normal, Exponential), the MLE has a closed-form solution. For most distributions used in hydrology (GEV, LP3, Weibull), the optimization must be solved numerically. The library uses constrained optimization with initial values derived from L-moment estimates.

**Fisher Information and standard errors.** The Fisher Information matrix quantifies the curvature of the log-likelihood surface at the maximum:

```math
\mathcal{I}(\boldsymbol{\theta}) = -E\left[\frac{\partial^2 \ell}{\partial \boldsymbol{\theta} \, \partial \boldsymbol{\theta}^T}\right]
```

Under regularity conditions, the MLE is asymptotically normal [[5]](#5):

```math
\sqrt{n}\left(\hat{\boldsymbol{\theta}} - \boldsymbol{\theta}\right) \xrightarrow{d} N\left(\mathbf{0},\; \mathcal{I}(\boldsymbol{\theta})^{-1}\right) \quad \text{as } n \to \infty
```

This provides approximate standard errors for each parameter:

```math
\text{SE}(\hat{\theta}_j) \approx \frac{1}{\sqrt{\mathcal{I}(\hat{\boldsymbol{\theta}})_{jj}}}
```

**Strengths:** Asymptotically efficient (achieves the lowest possible variance among consistent estimators), asymptotically unbiased, invariant under reparameterization, provides a natural framework for model comparison via AIC and BIC.

**Weaknesses:** Requires numerical optimization that may fail to converge, sensitive to outliers, can be biased and inefficient for small samples, requires specification of the full probability model.

### Using MLE

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

The Method of Moments (MOM) is the oldest and simplest approach to parameter estimation. The core idea is to equate sample moments to the corresponding theoretical moments of the distribution and solve for the unknown parameters.

### Mathematical Formulation

Given a sample $x_1, x_2, \ldots, x_n$, the first four sample moments are the mean, standard deviation, skewness, and kurtosis:

```math
\bar{x} = \frac{1}{n}\sum_{i=1}^{n} x_i
```

```math
s = \sqrt{\frac{1}{n-1}\sum_{i=1}^{n}(x_i - \bar{x})^2}
```

```math
\hat{\gamma} = \frac{n}{(n-1)(n-2)} \sum_{i=1}^{n}\left(\frac{x_i - \bar{x}}{s}\right)^3
```

```math
\hat{\kappa} = \frac{n(n+1)}{(n-1)(n-2)(n-3)} \sum_{i=1}^{n}\left(\frac{x_i - \bar{x}}{s}\right)^4 - \frac{3(n-1)^2}{(n-2)(n-3)}
```

The `Statistics.ProductMoments()` method returns these four quantities as the array $[\bar{x},\; s,\; \hat{\gamma},\; \hat{\kappa}]$.

MOM estimation sets the theoretical moments equal to the sample moments and solves for the distribution parameters. For a two-parameter distribution, only the first two moments (mean and standard deviation) are needed. For three-parameter distributions, skewness is also required.

**Example: Normal distribution.** The Normal($\mu$, $\sigma$) has $E[X] = \mu$ and $\text{SD}[X] = \sigma$. Equating sample to theoretical moments yields:

```math
\hat{\mu} = \bar{x}, \quad \hat{\sigma} = s
```

**Example: Gamma distribution.** The Gamma($\kappa$, $\theta$) has $E[X] = \kappa\theta$ and $\text{Var}[X] = \kappa\theta^2$. Solving for the parameters:

```math
\hat{\kappa} = \frac{\bar{x}^2}{s^2}, \quad \hat{\theta} = \frac{s^2}{\bar{x}}
```

**Strengths:** Simple, closed-form solutions, always produces estimates, computationally fast.

**Weaknesses:** Not statistically efficient (higher variance than MLE), can produce invalid parameters for skewed distributions, estimates are sensitive to outliers because conventional moments give disproportionate weight to extreme values.

### Using Method of Moments

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

## Method of Percentiles

The Method of Percentiles (also called least-squares fitting or quantile matching) estimates parameters by matching theoretical quantiles of the distribution to empirical quantiles computed from the data.

### Mathematical Formulation

Given a sorted sample $x_{(1)} \leq x_{(2)} \leq \cdots \leq x_{(n)}$, each observation is assigned a plotting position $p_i$ that estimates $F(x_{(i)})$. A common choice is the Weibull plotting position:

```math
p_i = \frac{i}{n + 1}
```

The parameters $\boldsymbol{\theta}$ are then chosen so that the theoretical quantile function (inverse CDF) matches the observed data as closely as possible. For a distribution with quantile function $F^{-1}(p;\,\boldsymbol{\theta})$, the parameters minimize the sum of squared differences:

```math
\hat{\boldsymbol{\theta}} = \underset{\boldsymbol{\theta}}{\text{argmin}} \sum_{i=1}^{n} \left[x_{(i)} - F^{-1}(p_i;\,\boldsymbol{\theta})\right]^2
```

For a two-parameter distribution, it is sufficient to select two percentiles (e.g., the median and the 84th percentile) and solve the resulting system of two equations:

```math
F^{-1}(p_j;\,\boldsymbol{\theta}) = x_{(j)}, \quad j \in \{j_1,\, j_2\}
```

**Strengths:** Intuitive and easy to visualize, always produces estimates, moderately robust to outliers in the tails, useful when expert judgment suggests specific quantile targets.

**Weaknesses:** Uses only selected data points or gives equal weight to all quantiles (not statistically efficient), lower precision than MLE or L-moments for most distributions.

## Estimation Method Comparison

The choice of estimation method depends on sample size, data quality, and application requirements. The following table summarizes the key trade-offs:

| Method | Efficiency | Small Samples | Robustness | Complexity | Best For |
|--------|-----------|---------------|-----------|-----------|---------|
| **MOM** | Low | Fair | Low | Simple | Quick estimates, stable distributions |
| **L-Moments** | Moderate--High | Excellent | High | Moderate | Hydrological data, small samples |
| **MLE** | Highest (asymptotic) | Poor--Fair | Low | Complex | Large samples, model comparison |
| **Percentiles** | Low | Fair | Moderate | Simple | Visual fitting, expert judgment |

### Rules of Thumb

- **n < 50:** Prefer L-moments. With small samples, robustness matters more than asymptotic efficiency, and L-moment estimates are nearly unbiased.
- **n > 100:** MLE becomes competitive and provides standard errors via Fisher Information, enabling confidence intervals and hypothesis tests.
- **Skewed distributions:** L-moments substantially outperform MOM, because conventional skewness estimates are highly variable for small samples.
- **US flood frequency analysis:** L-moments are recommended by USGS Bulletin 17C [[1]](#1). The Expected Moments Algorithm (EMA) extends the framework to handle censored and historical data.
- **Model selection:** When comparing candidate distributions, MLE enables the use of information criteria (AIC, BIC) for objective model ranking.

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

The Numerics library provides two formal statistical tests for detecting outliers in flood frequency data, both based on the Grubbs-Beck framework.

#### Grubbs-Beck Test

The original Grubbs-Beck test ([Grubbs & Beck, 1972][6]) identifies high and low outlier thresholds at the 10% significance level. The test assumes the log-transformed data follow a normal distribution and uses a critical value $K_n$ based on sample size:

```math
X_{Hi} = \exp(\bar{y} + K_n \cdot s_y), \quad X_{Lo} = \exp(\bar{y} - K_n \cdot s_y)
```

where $\bar{y}$ and $s_y$ are the mean and standard deviation of the log-transformed sample.

```cs
using Numerics.Data.Statistics;

double[] annualPeaks = { 50, 180, 220, 310, 450, 520, 680, 720, 890, 1050, 1200, 5200 };

// Original Grubbs-Beck test at 10% significance
MultipleGrubbsBeckTest.GrubbsBeckTest(annualPeaks, out double XHi, out double XLo);

Console.WriteLine($"High outlier threshold: {XHi:F0} (values above are high outliers)");
Console.WriteLine($"Low outlier threshold:  {XLo:F0} (values below are low outliers)");

// Identify outliers
var highOutliers = annualPeaks.Where(x => x > XHi).ToArray();
var lowOutliers = annualPeaks.Where(x => x < XLo).ToArray();
Console.WriteLine($"High outliers: {highOutliers.Length}, Low outliers: {lowOutliers.Length}");
```

#### Multiple Grubbs-Beck Test

The Multiple Grubbs-Beck Test ([Cohn et al., 2013][7]) is a generalization that can detect multiple potentially influential low flows (PILFs). This is particularly important in Bulletin 17C flood frequency analysis, where low outliers can distort the fitted distribution and bias upper-tail quantile estimates.

```cs
using Numerics.Data.Statistics;

double[] annualPeaks = { 5, 12, 18, 180, 220, 310, 450, 520, 680, 720, 890, 1050 };

// Multiple Grubbs-Beck test — returns count of low outliers
int numLowOutliers = MultipleGrubbsBeckTest.Function(annualPeaks);

Console.WriteLine($"Number of low outliers detected: {numLowOutliers}");

if (numLowOutliers > 0)
{
    // Sort to identify which values are low outliers
    double[] sorted = annualPeaks.OrderBy(x => x).ToArray();
    Console.WriteLine("Low outlier values:");
    for (int i = 0; i < numLowOutliers; i++)
        Console.WriteLine($"  {sorted[i]}");
    Console.WriteLine("Consider censoring these values in LP3 frequency analysis");
}
```

> **Guidance:** Use the original Grubbs-Beck test for general-purpose outlier screening. Use the Multiple Grubbs-Beck test specifically for Bulletin 17C flood frequency analysis where multiple low outliers may influence the Log-Pearson Type III fit.

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

<a id="5">[5]</a> Casella, G. & Berger, R. L. (2002). *Statistical Inference* (2nd ed.). Duxbury/Thomson.

<a id="6">[6]</a> Grubbs, F. E., & Beck, G. (1972). Extension of sample sizes and percentage points for significance tests of outlying observations. *Technometrics*, 14(4), 847-854.

<a id="7">[7]</a> Cohn, T. A., England, J. F., Berenbrock, C. E., Mason, R. R., Stedinger, J. R., & Lamontagne, J. R. (2013). A generalized Grubbs-Beck test statistic for detecting multiple potentially influential low outliers in flood series. *Water Resources Research*, 49(8), 5047-5058.

---

[← Previous: Univariate Distributions](univariate.md) | [Back to Index](../index.md) | [Next: Uncertainty Analysis →](uncertainty-analysis.md)
