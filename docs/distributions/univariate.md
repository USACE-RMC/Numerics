# Univariate Distributions

This reference documents the key univariate probability distributions available in ***Numerics***. Each distribution includes mathematical definitions, parameter descriptions, and usage examples.

## Table of Contents

1. [Normal (Gaussian)](#normal-gaussian)
2. [Log-Normal](#log-normal)
3. [Generalized Extreme Value (GEV)](#generalized-extreme-value-gev)
4. [Gumbel](#gumbel)
5. [Log-Pearson Type III](#log-pearson-type-iii)
6. [Pearson Type III (Gamma)](#pearson-type-iii-gamma)
7. [Exponential](#exponential)
8. [Weibull](#weibull)
9. [Beta](#beta)
10. [Uniform](#uniform)
11. [Triangular](#triangular)
12. [PERT](#pert)
13. [Generalized Pareto](#generalized-pareto)
14. [Empirical](#empirical)
15. [Mixture](#mixture)

---

## Normal (Gaussian)

The Normal distribution is the most widely used probability distribution, characterized by its symmetric bell-shaped curve [[1]](#ref1).

### Definition

```math
f(x|\mu,\sigma) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Mean | μ | $(-\infty, \infty)$ | Location (center) |
| Standard Deviation | σ | $(0, \infty)$ | Scale (spread) |

### Properties

| Property | Value |
|----------|-------|
| Mean | $\mu$ |
| Variance | $\sigma^2$ |
| Skewness | 0 |
| Kurtosis | 3 |
| Support | $(-\infty, \infty)$ |

### Usage

```cs
using Numerics.Distributions;

// Create with parameters
var normal = new Normal(100, 15);

// Standard normal
var stdNormal = new Normal(0, 1);

// Probability calculations
double p = normal.CDF(115);        // P(X ≤ 115)
double z = stdNormal.InverseCDF(0.975);  // ≈ 1.96

// Fit to data
double[] data = { 95, 102, 98, 105, 100 };
normal.SetParameters(normal.ParametersFromMoments(data));
```

---

## Log-Normal

A random variable whose logarithm follows a Normal distribution. Used for positive, right-skewed data [[1]](#ref1).

### Definition

```math
f(x|\mu,\sigma) = \frac{1}{x\sigma\sqrt{2\pi}} \exp\left(-\frac{(\ln x - \mu)^2}{2\sigma^2}\right), \quad x > 0
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Log-Mean | μ | $(-\infty, \infty)$ | Mean of ln(X) |
| Log-Std | σ | $(0, \infty)$ | Std dev of ln(X) |

### Properties

| Property | Value |
|----------|-------|
| Mean | $\exp(\mu + \sigma^2/2)$ |
| Variance | $[\exp(\sigma^2) - 1]\exp(2\mu + \sigma^2)$ |
| Mode | $\exp(\mu - \sigma^2)$ |
| Support | $(0, \infty)$ |

### Usage

```cs
var lognormal = new LogNormal(3.0, 0.5);

// Real-space mean is NOT μ
double mean = lognormal.Mean;  // = exp(3 + 0.25/2) ≈ 23.1

// Convert from real-space parameters
double realMean = 20;
double realStd = 10;
double sigma = Math.Sqrt(Math.Log(1 + (realStd/realMean)*(realStd/realMean)));
double mu = Math.Log(realMean) - sigma*sigma/2;
var dist = new LogNormal(mu, sigma);
```

---

## Generalized Extreme Value (GEV)

Models the distribution of block maxima (e.g., annual maximum floods) [[2]](#ref2).

### Definition

```math
F(x|\mu,\sigma,\xi) = \exp\left\{-\left[1 + \xi\left(\frac{x-\mu}{\sigma}\right)\right]^{-1/\xi}\right\}
```

For $\xi = 0$ (Gumbel case):
```math
F(x|\mu,\sigma,0) = \exp\left\{-\exp\left(-\frac{x-\mu}{\sigma}\right)\right\}
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Location | μ | $(-\infty, \infty)$ | Position |
| Scale | σ | $(0, \infty)$ | Spread |
| Shape | ξ | $(-\infty, \infty)$ | Tail behavior |

### Shape Parameter Interpretation

| ξ Value | Type | Tail | Example |
|---------|------|------|---------|
| ξ < 0 | Weibull | Bounded upper tail | Bounded phenomena |
| ξ = 0 | Gumbel | Light tail | Temperature extremes |
| ξ > 0 | Fréchet | Heavy tail | Flood peaks |

### Usage

```cs
var gev = new GeneralizedExtremeValue(1000, 200, -0.1);

// 100-year return level
double q100 = gev.InverseCDF(0.99);

// Fit to annual maxima
double[] annualMax = { 1200, 1500, 1100, 1800, 1350, 1600 };
gev.SetParameters(gev.ParametersFromLinearMoments(annualMax));
```

---

## Gumbel

Special case of GEV with ξ = 0, used for extreme values with light (exponential) tails [[2]](#ref2).

### Definition

```math
f(x|\mu,\sigma) = \frac{1}{\sigma}\exp\left(-\frac{x-\mu}{\sigma} - \exp\left(-\frac{x-\mu}{\sigma}\right)\right)
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Location | μ | $(-\infty, \infty)$ | Mode of the distribution |
| Scale | σ | $(0, \infty)$ | Spread |

### Properties

| Property | Value |
|----------|-------|
| Mean | $\mu + \gamma\sigma$ where $\gamma \approx 0.5772$ (Euler's constant) |
| Variance | $\pi^2\sigma^2/6$ |
| Mode | $\mu$ |
| Skewness | $\approx 1.14$ |

### Usage

```cs
var gumbel = new Gumbel(1000, 200);

// Relationship to GEV
var gevEquivalent = new GeneralizedExtremeValue(1000, 200, 0);

// Return periods
double T = 100;
double p = 1 - 1/T;  // Non-exceedance probability
double qT = gumbel.InverseCDF(p);
```

---

## Log-Pearson Type III

The standard distribution for flood frequency analysis in the United States (Bulletin 17C) [[3]](#ref3).

### Definition

The logarithms of the data follow a Pearson Type III (3-parameter Gamma) distribution.

If $Y = \ln(X)$ follows Pearson Type III with parameters $(\mu_Y, \sigma_Y, \gamma_Y)$:

```math
f_Y(y) = \frac{|\beta|^\alpha}{\Gamma(\alpha)}(y - \xi)^{\alpha-1}\exp(-\beta(y-\xi))
```

where $\alpha = 4/\gamma_Y^2$, $\beta = \text{sign}(\gamma_Y) \cdot 2/(\sigma_Y|\gamma_Y|)$, $\xi = \mu_Y - 2\sigma_Y/\gamma_Y$.

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Log-Mean | μ | $(-\infty, \infty)$ | Mean of ln(X) |
| Log-Std | σ | $(0, \infty)$ | Std dev of ln(X) |
| Log-Skew | γ | $(-\infty, \infty)$ | Skewness of ln(X) |

### Usage

```cs
// Create LP3 distribution
var lp3 = new LogPearsonTypeIII(3.0, 0.3, 0.5);

// Fit to annual maximum flows
double[] flows = { 12500, 15200, 11800, 18900, 14200, 16500 };
lp3.SetParameters(lp3.ParametersFromLinearMoments(flows));

// Compute flood quantiles
Console.WriteLine($"Q10:  {lp3.InverseCDF(0.90):N0} cfs");
Console.WriteLine($"Q100: {lp3.InverseCDF(0.99):N0} cfs");
Console.WriteLine($"Q500: {lp3.InverseCDF(0.998):N0} cfs");
```

---

## Pearson Type III (Gamma)

Three-parameter Gamma distribution, the parent of Log-Pearson Type III [[4]](#ref4).

### Definition

```math
f(x|\mu,\sigma,\gamma) = \frac{|\beta|^\alpha}{\Gamma(\alpha)}(x - \xi)^{\alpha-1}\exp(-\beta(x-\xi))
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Mean | μ | $(-\infty, \infty)$ | Location |
| Std Dev | σ | $(0, \infty)$ | Scale |
| Skewness | γ | $(-\infty, \infty)$ | Shape |

### Usage

```cs
var p3 = new PearsonTypeIII(100, 25, 1.5);

// When skewness = 0, equivalent to Normal
var normalEquiv = new PearsonTypeIII(100, 25, 0);

// Fit to data
p3.SetParameters(p3.ParametersFromLinearMoments(data));
```

---

## Exponential

Models time between events in a Poisson process [[1]](#ref1).

### Definition

```math
f(x|\lambda) = \lambda e^{-\lambda x}, \quad x \geq 0
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Rate | λ | $(0, \infty)$ | Events per unit time |

### Properties

| Property | Value |
|----------|-------|
| Mean | $1/\lambda$ |
| Variance | $1/\lambda^2$ |
| Skewness | 2 |
| Memoryless | Yes |

### Usage

```cs
// Average of 5 events per year → rate = 5
var exp = new Exponential(5);

// Probability of waiting more than 0.5 years
double p = exp.CCDF(0.5);  // = exp(-5 × 0.5) = exp(-2.5)

// Median waiting time
double median = exp.InverseCDF(0.5);  // = ln(2)/λ
```

---

## Weibull

Widely used in reliability engineering and survival analysis [[5]](#ref5).

### Definition

```math
f(x|k,\lambda) = \frac{k}{\lambda}\left(\frac{x}{\lambda}\right)^{k-1}\exp\left(-\left(\frac{x}{\lambda}\right)^k\right), \quad x \geq 0
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Shape | k | $(0, \infty)$ | Controls failure rate trend |
| Scale | λ | $(0, \infty)$ | Characteristic life |

### Shape Parameter Interpretation

| k Value | Failure Rate | Application |
|---------|--------------|-------------|
| k < 1 | Decreasing | Infant mortality |
| k = 1 | Constant | Random failures (Exponential) |
| k > 1 | Increasing | Wear-out failures |
| k ≈ 3.6 | Approximates Normal | |

### Usage

```cs
var weibull = new Weibull(2.5, 1000);

// Reliability at time t
double R = weibull.CCDF(500);  // Survival probability

// B10 life (10% failure probability)
double B10 = weibull.InverseCDF(0.10);
```

---

## Beta

Models random variables bounded between 0 and 1, useful for proportions and probabilities [[1]](#ref1).

### Definition

```math
f(x|\alpha,\beta) = \frac{x^{\alpha-1}(1-x)^{\beta-1}}{B(\alpha,\beta)}, \quad 0 \leq x \leq 1
```

where $B(\alpha,\beta) = \Gamma(\alpha)\Gamma(\beta)/\Gamma(\alpha+\beta)$.

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Alpha | α | $(0, \infty)$ | Shape parameter 1 |
| Beta | β | $(0, \infty)$ | Shape parameter 2 |

### Special Cases

| α | β | Shape |
|---|---|-------|
| 1 | 1 | Uniform |
| 0.5 | 0.5 | Arcsine |
| α = β | | Symmetric |
| α < β | | Left-skewed |
| α > β | | Right-skewed |

### Usage

```cs
var beta = new Beta(2, 5);

// Scale to [a, b] interval
double a = 10, b = 50;
double scaledValue = a + (b - a) * beta.InverseCDF(0.5);

// Fit to bounded data (normalize first)
double[] proportions = { 0.15, 0.22, 0.18, 0.25, 0.20 };
beta.SetParameters(beta.ParametersFromMoments(proportions));
```

---

## Uniform

Equal probability across a bounded interval [[1]](#ref1).

### Definition

```math
f(x|a,b) = \frac{1}{b-a}, \quad a \leq x \leq b
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Minimum | a | $(-\infty, b)$ | Lower bound |
| Maximum | b | $(a, \infty)$ | Upper bound |

### Properties

| Property | Value |
|----------|-------|
| Mean | $(a+b)/2$ |
| Variance | $(b-a)^2/12$ |
| Skewness | 0 |

### Usage

```cs
var uniform = new Uniform(10, 50);

// CDF is linear
double p = uniform.CDF(30);  // = (30-10)/(50-10) = 0.5

// Random sampling
double x = uniform.InverseCDF(new Random().NextDouble());
```

---

## Triangular

Simple three-parameter distribution defined by minimum, mode, and maximum [[6]](#ref6).

### Definition

```math
f(x|a,c,b) = \begin{cases}
\frac{2(x-a)}{(b-a)(c-a)} & a \leq x \leq c \\
\frac{2(b-x)}{(b-a)(b-c)} & c < x \leq b
\end{cases}
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Minimum | a | $(-\infty, c)$ | Lower bound |
| Mode | c | $[a, b]$ | Most likely value |
| Maximum | b | $(c, \infty)$ | Upper bound |

### Usage

```cs
// Expert estimate: min=10, most likely=25, max=50
var tri = new Triangular(10, 25, 50);

double mean = tri.Mean;  // = (10 + 25 + 50)/3 ≈ 28.3
```

---

## PERT

Modified Beta distribution commonly used in project management, smoother than Triangular [[6]](#ref6).

### Definition

PERT uses a Beta distribution scaled to $[a, b]$ with shape parameters derived from the mode:

```math
\alpha = 1 + \lambda\frac{c - a}{b - a}, \quad \beta = 1 + \lambda\frac{b - c}{b - a}
```

where $\lambda = 4$ is the default shape parameter.

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Minimum | a | $(-\infty, c)$ | Optimistic estimate |
| Mode | c | $[a, b]$ | Most likely estimate |
| Maximum | b | $(c, \infty)$ | Pessimistic estimate |

### Usage

```cs
// Project duration: optimistic=5, likely=8, pessimistic=15 days
var pert = new PERT(5, 8, 15);

// PERT mean weights the mode more heavily than Triangular
double pertMean = pert.Mean;  // ≈ (5 + 4×8 + 15)/6 = 8.67
```

---

## Generalized Pareto

Models exceedances over a threshold, fundamental to Peaks-Over-Threshold analysis [[2]](#ref2).

### Definition

```math
F(x|\mu,\sigma,\xi) = 1 - \left[1 + \xi\frac{x-\mu}{\sigma}\right]^{-1/\xi}
```

### Parameters

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Location | μ | $(-\infty, \infty)$ | Threshold |
| Scale | σ | $(0, \infty)$ | Spread |
| Shape | ξ | $(-\infty, \infty)$ | Tail behavior |

### Usage

```cs
// Model flood peaks above threshold of 1000 cfs
var gpd = new GeneralizedPareto(1000, 500, 0.2);

// Probability of exceeding 2500 cfs given > 1000 cfs
double pExceed = gpd.CCDF(2500);
```

---

## Empirical

Non-parametric distribution constructed directly from data using plotting positions [[7]](#ref7).

### Definition

The empirical CDF is defined by the sample order statistics and plotting positions:

```math
\hat{F}(x_{(i)}) = \frac{i - a}{n + 1 - 2a}
```

where $a$ is the plotting position parameter (Weibull: $a=0$, Cunnane: $a=0.4$).

### Usage

```cs
double[] data = { 10, 15, 12, 18, 14, 16, 13 };

// Create empirical distribution
var empirical = new Empirical(data);

// Interpolated quantiles
double median = empirical.InverseCDF(0.5);

// With specified plotting position
var emp2 = new Empirical(data, PlottingPosition.Cunnane);
```

---

## Mixture

Weighted combination of multiple distributions [[8]](#ref8).

### Definition

```math
f(x) = \sum_{i=1}^{k} w_i f_i(x), \quad \sum_{i=1}^{k} w_i = 1
```

### Usage

```cs
// Bimodal distribution: 70% N(0,1) + 30% N(5,1)
var components = new List<IUnivariateDistribution>
{
    new Normal(0, 1),
    new Normal(5, 1)
};
var weights = new double[] { 0.7, 0.3 };

var mixture = new Mixture(components, weights);

// Sample from mixture
double[] samples = new double[1000];
mixture.GenerateRandomValues(samples);
```

---

## Summary Table

| Distribution | Parameters | Support | Skewness | Common Use |
|-------------|------------|---------|----------|------------|
| Normal | μ, σ | $(-\infty, \infty)$ | 0 | General modeling |
| LogNormal | μ, σ | $(0, \infty)$ | >0 | Positive data |
| GEV | μ, σ, ξ | Varies | Varies | Block maxima |
| Gumbel | μ, σ | $(-\infty, \infty)$ | 1.14 | Extremes (light tail) |
| Log-Pearson III | μ, σ, γ | $(0, \infty)$ | Varies | Flood frequency |
| Pearson III | μ, σ, γ | Varies | γ | General |
| Exponential | λ | $[0, \infty)$ | 2 | Waiting times |
| Weibull | k, λ | $[0, \infty)$ | Varies | Reliability |
| Beta | α, β | $[0, 1]$ | Varies | Proportions |
| Uniform | a, b | $[a, b]$ | 0 | Equal likelihood |
| Triangular | a, c, b | $[a, b]$ | Varies | Expert judgment |
| PERT | a, c, b | $[a, b]$ | Varies | Project management |
| Gen. Pareto | μ, σ, ξ | $[\mu, \infty)$ | Varies | Threshold exceedances |
| Empirical | data | Data range | Data | Non-parametric |
| Mixture | components | Union | Varies | Multi-modal |

---

## References

<a id="ref1">[1]</a> Johnson, N. L., Kotz, S., & Balakrishnan, N. (1994). *Continuous Univariate Distributions* (2nd ed., Vols. 1-2). Wiley.

<a id="ref2">[2]</a> Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.

<a id="ref3">[3]</a> England, J. F., et al. (2019). Guidelines for Determining Flood Flow Frequency—Bulletin 17C. *U.S. Geological Survey Techniques and Methods*, Book 4, Chapter B5.

<a id="ref4">[4]</a> Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society: Series B*, 52(1), 105-124.

<a id="ref5">[5]</a> Weibull, W. (1951). A statistical distribution function of wide applicability. *Journal of Applied Mechanics*, 18(3), 293-297.

<a id="ref6">[6]</a> Vose, D. (2008). *Risk Analysis: A Quantitative Guide* (3rd ed.). Wiley.

<a id="ref7">[7]</a> Cunnane, C. (1978). Unbiased plotting positions—A review. *Journal of Hydrology*, 37(3-4), 205-222.

<a id="ref8">[8]</a> McLachlan, G., & Peel, D. (2000). *Finite Mixture Models*. Wiley.
