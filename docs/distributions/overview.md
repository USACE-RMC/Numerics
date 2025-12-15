# Probability Distributions Overview

The ***Numerics*** library provides a comprehensive collection of probability distributions for statistical modeling, risk analysis, and Monte Carlo simulation. This document describes the common interface and concepts shared across all distribution types.

## Distribution Hierarchy

```
IUnivariateDistribution (interface)
    └── UnivariateDistributionBase (abstract)
            ├── Normal
            ├── LogNormal
            ├── GeneralizedExtremeValue
            ├── Gumbel
            ├── LogPearsonTypeIII
            ├── ... (40+ distributions)
            └── Mixture

IMultivariateDistribution (interface)
    └── MultivariateDistribution (abstract)
            ├── MultivariateNormal
            └── BivariateEmpirical

IBivariateCopula (interface)
    └── BivariateCopula (abstract)
            ├── NormalCopula
            ├── ClaytonCopula
            ├── FrankCopula
            ├── GumbelCopula
            └── JoeCopula
```

## Core Interface: IUnivariateDistribution

All univariate distributions implement `IUnivariateDistribution`, which defines the standard probability distribution functions:

### Probability Functions

| Method | Description | Mathematical Definition |
|--------|-------------|------------------------|
| `PDF(x)` | Probability Density Function | $f(x)$ |
| `LogPDF(x)` | Log of the PDF (for numerical stability) | $\ln f(x)$ |
| `CDF(x)` | Cumulative Distribution Function | $F(x) = P(X \leq x)$ |
| `InverseCDF(p)` | Quantile Function | $F^{-1}(p) = x : F(x) = p$ |
| `CCDF(x)` | Complementary CDF | $\bar{F}(x) = 1 - F(x) = P(X > x)$ |

### Distribution Properties

| Property | Description |
|----------|-------------|
| `Mean` | Expected value $E[X]$ |
| `Median` | 50th percentile |
| `Mode` | Most likely value |
| `Variance` | $\text{Var}(X) = E[(X - \mu)^2]$ |
| `StandardDeviation` | $\sigma = \sqrt{\text{Var}(X)}$ |
| `Skew` | Third standardized moment |
| `Kurtosis` | Fourth standardized moment |
| `Minimum` | Lower bound of support |
| `Maximum` | Upper bound of support |

### Parameter Access

| Property/Method | Description |
|-----------------|-------------|
| `NumberOfParameters` | Count of distribution parameters |
| `GetParameters()` | Returns array of current parameters |
| `SetParameters(double[])` | Sets parameters from array |
| `ParametersValid` | Checks if current parameters are valid |

## Parameter Estimation Interfaces

Distributions support different estimation methods through specialized interfaces:

### IMomentEstimation

Method of Moments estimation using product moments:

```cs
public interface IMomentEstimation
{
    double[] ParametersFromMoments(IList<double> data);
    double[] MomentsFromParameters { get; }
}
```

### ILinearMomentEstimation

L-Moment estimation, preferred for heavy-tailed distributions [[1]](#ref1):

```cs
public interface ILinearMomentEstimation
{
    double[] ParametersFromLinearMoments(IList<double> data);
    double[] LinearMomentsFromParameters { get; }
}
```

### IMaximumLikelihoodEstimation

Maximum Likelihood Estimation for efficient parameter estimates:

```cs
public interface IMaximumLikelihoodEstimation
{
    double[] MLE(IList<double> data);
    double LogLikelihood(IList<double> data);
}
```

## Working with Distributions

### Creating Distributions

```cs
using Numerics.Distributions;

// Direct parameter specification
var normal = new Normal(100, 15);  // mean=100, std=15
var gev = new GeneralizedExtremeValue(1000, 200, -0.1);

// Default constructor (uses standard parameters)
var stdNormal = new Normal();  // N(0,1)

// From parameter array
var dist = new Normal();
dist.SetParameters(new double[] { 50, 10 });
```

### Computing Probabilities

```cs
var dist = new Normal(100, 15);

// What is the probability of X ≤ 120?
double p = dist.CDF(120);  // ≈ 0.9088

// What is the probability of X > 120?
double q = dist.CCDF(120);  // ≈ 0.0912

// What value has 95% probability of not being exceeded?
double x95 = dist.InverseCDF(0.95);  // ≈ 124.67

// What is the probability density at x = 100?
double density = dist.PDF(100);  // ≈ 0.0266
```

### Generating Random Samples

```cs
var dist = new LogNormal(3.0, 0.5);

// Generate 10,000 random samples
double[] samples = new double[10000];
dist.GenerateRandomValues(samples);

// With reproducible seed
dist.GenerateRandomValues(samples, seed: 42);

// Verify sample statistics match distribution
double sampleMean = Statistics.Mean(samples);
double trueMean = dist.Mean;
Console.WriteLine($"Sample mean: {sampleMean:F3}, True mean: {trueMean:F3}");
```

### Fitting to Data

```cs
double[] data = { 12.5, 15.2, 11.8, 18.9, 14.2, 16.5, 13.4 };

// L-Moments (robust, recommended for most cases)
var gev = new GeneralizedExtremeValue();
gev.SetParameters(gev.ParametersFromLinearMoments(data));

// Method of Moments
var normal = new Normal();
normal.SetParameters(normal.ParametersFromMoments(data));

// Maximum Likelihood
var lognormal = new LogNormal();
lognormal.SetParameters(lognormal.MLE(data));
```

## Available Distributions

### Continuous Distributions

| Distribution | Parameters | Common Use |
|-------------|------------|------------|
| Normal | μ (location), σ (scale) | General modeling |
| LogNormal | μ (log-mean), σ (log-std) | Positive skewed data |
| Exponential | λ (rate) | Time between events |
| Gamma | α (shape), β (rate) | Waiting times |
| Beta | α, β (shape parameters) | Proportions, probabilities |
| Uniform | a (min), b (max) | Equal likelihood |
| Triangular | min, mode, max | Expert elicitation |
| PERT | min, mode, max | Project management |
| Weibull | k (shape), λ (scale) | Reliability, survival |
| Gumbel | μ (location), σ (scale) | Extreme values (maxima) |
| GEV | μ, σ, ξ | Extreme values (general) |
| Generalized Pareto | μ, σ, ξ | Threshold exceedances |
| Pearson Type III | μ, σ, γ | Flood frequency |
| Log-Pearson Type III | μ, σ, γ | Flood frequency (US standard) |
| Rayleigh | σ (scale) | Wind speeds |
| Pareto | x_m (scale), α (shape) | Income, city sizes |

### Discrete Distributions

| Distribution | Parameters | Common Use |
|-------------|------------|------------|
| Binomial | n (trials), p (probability) | Success counts |
| Poisson | λ (rate) | Event counts |
| Geometric | p (probability) | Trials until success |
| Negative Binomial | r (successes), p (probability) | Overdispersed counts |

### Special Distributions

| Distribution | Description |
|-------------|-------------|
| Empirical | Non-parametric distribution from data |
| Mixture | Weighted combination of distributions |
| Truncated | Any distribution with bounds |
| Kernel Density | Smooth non-parametric estimate |

## Hydrologic Frequency Distributions

For flood frequency analysis, ***Numerics*** provides distributions commonly used in hydrology [[2]](#ref2) [[3]](#ref3):

```cs
// US standard: Log-Pearson Type III (Bulletin 17C)
var lp3 = new LogPearsonTypeIII();
lp3.SetParameters(lp3.ParametersFromLinearMoments(annualMaxFlows));

// UK standard: Generalized Logistic
var glo = new GeneralizedLogistic();
glo.SetParameters(glo.ParametersFromLinearMoments(annualMaxFlows));

// European: Generalized Extreme Value
var gev = new GeneralizedExtremeValue();
gev.SetParameters(gev.ParametersFromLinearMoments(annualMaxFlows));
```

## Thread Safety

Distribution objects are **not thread-safe**. For parallel applications, create separate distribution instances per thread:

```cs
Parallel.For(0, Environment.ProcessorCount, i =>
{
    // Each thread gets its own distribution instance
    var dist = new Normal(100, 15);
    double[] localSamples = new double[1000];
    dist.GenerateRandomValues(localSamples, seed: i);
    // Process samples...
});
```

---

## References

<a id="ref1">[1]</a> Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society: Series B*, 52(1), 105-124.

<a id="ref2">[2]</a> Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.

<a id="ref3">[3]</a> Stedinger, J. R., Vogel, R. M., & Foufoula-Georgiou, E. (1993). Frequency analysis of extreme events. In D. R. Maidment (Ed.), *Handbook of Hydrology* (Chapter 18). McGraw-Hill.
