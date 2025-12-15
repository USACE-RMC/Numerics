# Parameter Estimation

This document describes the parameter estimation methods available in ***Numerics*** for fitting probability distributions to observed data.

## Overview

***Numerics*** supports three primary estimation methods:

| Method | Interface | Best For |
|--------|-----------|----------|
| Method of Moments (MOM) | `IMomentEstimation` | Quick estimates, symmetric distributions |
| L-Moments (LMOM) | `ILinearMomentEstimation` | Robust estimates, heavy-tailed distributions |
| Maximum Likelihood (MLE) | `IMaximumLikelihoodEstimation` | Efficient estimates, large samples |

## Method of Moments (MOM)

The Method of Moments equates sample moments to theoretical moments and solves for parameters [[1]](#ref1).

### Theory

For a distribution with $k$ parameters, MOM uses the first $k$ sample moments:

**Sample Moments:**
```math
m_r = \frac{1}{n}\sum_{i=1}^{n}x_i^r
```

**Central Moments:**
```math
\mu_r = \frac{1}{n}\sum_{i=1}^{n}(x_i - \bar{x})^r
```

The parameters are found by solving the system of equations where theoretical moments equal sample moments.

### Usage

```cs
using Numerics.Distributions;

double[] data = { 12.5, 15.2, 11.8, 18.9, 14.2, 16.5, 13.4, 17.8 };

// Normal distribution - uses mean and variance
var normal = new Normal();
if (normal is IMomentEstimation mom)
{
    double[] parameters = mom.ParametersFromMoments(data);
    normal.SetParameters(parameters);
    Console.WriteLine($"μ = {normal.Mean:F2}, σ = {normal.StandardDeviation:F2}");
}

// Get theoretical moments from fitted parameters
double[] theoreticalMoments = normal.MomentsFromParameters;
```

### Advantages and Limitations

**Advantages:**
- Simple, closed-form solutions for many distributions
- Computationally fast
- Intuitive interpretation

**Limitations:**
- Can be inefficient for small samples
- Sensitive to outliers (especially higher moments)
- May produce invalid parameters for some distributions

---

## L-Moments (LMOM)

L-Moments are linear combinations of order statistics, providing robust and efficient parameter estimates [[2]](#ref2).

### Theory

The $r$-th L-moment is defined as:

```math
\lambda_r = \frac{1}{r} \sum_{k=0}^{r-1} (-1)^k \binom{r-1}{k} E[X_{r-k:r}]
```

**Sample L-Moments** are computed from the probability-weighted moments:

```math
b_r = \frac{1}{n} \sum_{i=1}^{n} \frac{(i-1)(i-2)\cdots(i-r)}{(n-1)(n-2)\cdots(n-r)} x_{(i)}
```

**L-Moment Ratios** (dimensionless):
- L-CV: $\tau_2 = \lambda_2/\lambda_1$ (coefficient of L-variation)
- L-Skewness: $\tau_3 = \lambda_3/\lambda_2$
- L-Kurtosis: $\tau_4 = \lambda_4/\lambda_2$

### Usage

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

double[] data = { 12500, 15200, 11800, 18900, 14200, 16500, 13400, 17800 };

// Compute sample L-moments
double[] lmoms = Statistics.LinearMoments(data);
Console.WriteLine($"L1 (location): {lmoms[0]:F1}");
Console.WriteLine($"L2 (scale):    {lmoms[1]:F1}");
Console.WriteLine($"τ3 (L-skew):   {lmoms[2]:F3}");
Console.WriteLine($"τ4 (L-kurt):   {lmoms[3]:F3}");

// Fit GEV using L-Moments
var gev = new GeneralizedExtremeValue();
if (gev is ILinearMomentEstimation lmom)
{
    gev.SetParameters(lmom.ParametersFromLinearMoments(data));
    Console.WriteLine($"\nGEV Parameters:");
    Console.WriteLine($"  Location (μ): {gev.Xi:F1}");
    Console.WriteLine($"  Scale (σ):    {gev.Alpha:F1}");
    Console.WriteLine($"  Shape (ξ):    {gev.Kappa:F3}");
}
```

### L-Moment Ratio Diagram

L-moment ratios can help identify appropriate distribution families:

```cs
// Compare sample L-moments to theoretical values
double[] sample = Statistics.LinearMoments(data);
double sampleL3 = sample[2];  // L-skewness
double sampleL4 = sample[3];  // L-kurtosis

// Theoretical L-moment ratios for various distributions
// GEV: τ4 ≈ 0.1504 + 0.1052τ3²  
// Normal: τ3 = 0, τ4 = 0.1226
// Gumbel: τ3 = 0.1699, τ4 = 0.1504

Console.WriteLine($"Sample: τ3={sampleL3:F3}, τ4={sampleL4:F3}");
```

### Advantages and Limitations

**Advantages:**
- More robust to outliers than conventional moments
- More efficient for small samples
- Better characterization of distribution tails
- Unique L-moment ratios for each distribution family

**Limitations:**
- Requires ordered data
- Less familiar to some practitioners
- Fewer analytical results available

---

## Maximum Likelihood Estimation (MLE)

MLE finds parameters that maximize the probability of observing the data [[3]](#ref3).

### Theory

The likelihood function for independent observations is:

```math
L(\theta|x_1,\ldots,x_n) = \prod_{i=1}^{n} f(x_i|\theta)
```

The log-likelihood (more numerically stable):

```math
\ell(\theta) = \sum_{i=1}^{n} \ln f(x_i|\theta)
```

MLE finds $\hat{\theta} = \arg\max_\theta \ell(\theta)$.

### Usage

```cs
using Numerics.Distributions;

double[] data = { 12.5, 15.2, 11.8, 18.9, 14.2, 16.5, 13.4, 17.8 };

// Fit using MLE
var normal = new Normal();
if (normal is IMaximumLikelihoodEstimation mle)
{
    double[] parameters = mle.MLE(data);
    normal.SetParameters(parameters);
    
    // Get log-likelihood value
    double ll = mle.LogLikelihood(data);
    Console.WriteLine($"Log-likelihood: {ll:F2}");
    
    // Compute AIC for model comparison
    int k = normal.NumberOfParameters;
    double aic = -2 * ll + 2 * k;
    Console.WriteLine($"AIC: {aic:F2}");
}
```

### Constrained MLE

For distributions with parameter constraints:

```cs
// Log-Pearson III with regional skew constraint
var lp3 = new LogPearsonTypeIII();

// Unconstrained MLE
lp3.SetParameters(lp3.MLE(data));
double unconstrainedSkew = lp3.Gamma;

// Apply Bulletin 17C weighted skew
double stationSkew = unconstrainedSkew;
double regionalSkew = 0.0;  // From regional analysis
double mseStation = 0.3;    // MSE of station skew
double mseRegional = 0.15;  // MSE of regional skew

double weightedSkew = (mseRegional * stationSkew + mseStation * regionalSkew) 
                    / (mseStation + mseRegional);
```

### Advantages and Limitations

**Advantages:**
- Asymptotically efficient (minimum variance for large samples)
- Provides likelihood values for model comparison
- Well-developed theory for confidence intervals
- Handles censored and truncated data

**Limitations:**
- May require numerical optimization
- Can be sensitive to starting values
- May not have closed-form solutions
- Can be biased for small samples

---

## Comparing Estimation Methods

### Method Selection Guidelines

| Scenario | Recommended Method |
|----------|-------------------|
| Small sample (n < 30) | L-Moments |
| Heavy-tailed distribution | L-Moments |
| Symmetric distribution | MOM or MLE |
| Large sample (n > 100) | MLE |
| Outliers present | L-Moments |
| Model comparison needed | MLE (for AIC/BIC) |
| Quick exploratory analysis | MOM |

### Example Comparison

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

double[] data = { 12500, 15200, 11800, 18900, 14200, 16500, 
                  13400, 17800, 25000, 9800, 14800, 16200 };

Console.WriteLine("GEV Parameter Estimates by Method:");
Console.WriteLine("Method       Location    Scale      Shape");
Console.WriteLine("------       --------    -----      -----");

// Method of Moments
var gevMOM = new GeneralizedExtremeValue();
gevMOM.SetParameters(gevMOM.ParametersFromMoments(data));
Console.WriteLine($"MOM          {gevMOM.Xi,8:F1}    {gevMOM.Alpha,5:F1}    {gevMOM.Kappa,7:F4}");

// L-Moments
var gevLMOM = new GeneralizedExtremeValue();
gevLMOM.SetParameters(gevLMOM.ParametersFromLinearMoments(data));
Console.WriteLine($"L-Moments    {gevLMOM.Xi,8:F1}    {gevLMOM.Alpha,5:F1}    {gevLMOM.Kappa,7:F4}");

// MLE
var gevMLE = new GeneralizedExtremeValue();
gevMLE.SetParameters(gevMLE.MLE(data));
Console.WriteLine($"MLE          {gevMLE.Xi,8:F1}    {gevMLE.Alpha,5:F1}    {gevMLE.Kappa,7:F4}");

// Compare 100-year quantiles
Console.WriteLine($"\n100-year Flood Estimates:");
Console.WriteLine($"MOM:       {gevMOM.InverseCDF(0.99):N0} cfs");
Console.WriteLine($"L-Moments: {gevLMOM.InverseCDF(0.99):N0} cfs");
Console.WriteLine($"MLE:       {gevMLE.InverseCDF(0.99):N0} cfs");
```

---

## Distribution Selection

### Fitting Multiple Distributions

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

double[] annualMax = { /* your data */ };

// Candidate distributions
var candidates = new List<IUnivariateDistribution>
{
    new Normal(),
    new LogNormal(),
    new GeneralizedExtremeValue(),
    new Gumbel(),
    new LogPearsonTypeIII(),
    new PearsonTypeIII()
};

Console.WriteLine("Distribution         AIC       BIC      K-S p-value");
Console.WriteLine("------------         ---       ---      -----------");

foreach (var dist in candidates)
{
    string name = dist.GetType().Name;
    
    // Fit using L-Moments (most robust)
    if (dist is ILinearMomentEstimation lmom)
        dist.SetParameters(lmom.ParametersFromLinearMoments(annualMax));
    else if (dist is IMomentEstimation mom)
        dist.SetParameters(mom.ParametersFromMoments(annualMax));
    
    // Compute criteria
    int n = annualMax.Length;
    double ll = dist.LogLikelihood(annualMax);
    int k = dist.NumberOfParameters;
    double aic = GoodnessOfFit.AIC(k, ll);
    double bic = GoodnessOfFit.BIC(n, k, ll);
    
    // K-S test
    var ks = GoodnessOfFit.KolmogorovSmirnov(annualMax, dist);
    
    Console.WriteLine($"{name,-20} {aic,7:F1}   {bic,7:F1}   {ks.PValue,8:F4}");
}
```

---

## References

<a id="ref1">[1]</a> Johnson, N. L., Kotz, S., & Balakrishnan, N. (1994). *Continuous Univariate Distributions* (2nd ed.). Wiley.

<a id="ref2">[2]</a> Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society: Series B*, 52(1), 105-124.

<a id="ref3">[3]</a> Pawitan, Y. (2001). *In All Likelihood: Statistical Modelling and Inference Using Likelihood*. Oxford University Press.

<a id="ref4">[4]</a> Hosking, J. R. M., & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press.

<a id="ref5">[5]</a> England, J. F., et al. (2019). Guidelines for Determining Flood Flow Frequency—Bulletin 17C. *U.S. Geological Survey Techniques and Methods*, Book 4, Chapter B5.
