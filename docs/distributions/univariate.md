# Univariate Distributions

[← Previous: Hypothesis Tests](../statistics/hypothesis-tests.md) | [Back to Index](../index.md) | [Next: Parameter Estimation →](parameter-estimation.md)

The ***Numerics*** library provides 43 univariate probability distributions for statistical analysis, risk assessment, and uncertainty quantification. All distributions implement a common interface with consistent methods for computing probability density functions (PDF), cumulative distribution functions (CDF), quantiles, and statistical moments.

## Available Distributions

### Continuous Distributions

| Distribution | Parameters | Typical Applications |
|--------------|-----------|---------------------|
| **Normal** | μ (mean), σ (std dev) | General purpose, natural phenomena |
| **Log-Normal** | μ, σ | Right-skewed data, multiplicative processes |
| **Uniform** | a (min), b (max) | Maximum entropy, prior distributions |
| **Exponential** | ξ (location), α (scale) | Time between events, survival analysis |
| **Gamma** | θ (scale), κ (shape) | Waiting times, rainfall |
| **Beta** | α, β | Probabilities, proportions, [0,1] bounded |
| **Weibull** | λ (scale), κ (shape) | Failure times, wind speed |
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
| **Student's t** | μ (location), σ (scale), ν (degrees of freedom) | Heavy-tailed alternative to Normal |
| **Noncentral t** | ν, μ (noncentrality) | Power analysis, hypothesis testing |
| **Chi-Squared** | k (degrees of freedom) | Variance estimation, goodness-of-fit |
| **Inverse Gamma** | α, β | Bayesian priors for variance |
| **Inverse Chi-Squared** | ν, τ² (scale) | Bayesian inference |
| **Pareto** | xₘ (scale), α (shape) | Income distributions, city sizes |
| **PERT** | a, b, c | Project management, expert judgment |
| **PERT Percentile** | P₅, P₅₀, P₉₅ | Expert percentile elicitation |
| **PERT Percentile Z** | Similar to PERT Percentile | Alternative parametrization |
| **Truncated Normal** | μ, σ, a, b | Bounded normal distributions |
| **Truncated Distribution** | Any distribution + bounds | Bounded versions of distributions |
| **Mixture** | Multiple distributions | Multi-modal data |
| **Empirical** | Sample data | Non-parametric, data-driven |
| **Kernel Density** | Sample data, bandwidth | Smooth non-parametric estimation |
| **Deterministic** | Single value | Point estimates, constants |
| **Competing Risks** | Multiple distributions | Failure analysis with multiple causes |
| **Von Mises** | μ (mean direction), κ (concentration) | Circular data, flood seasonality |

### Discrete Distributions

| Distribution | Parameters | Typical Applications |
|--------------|-----------|---------------------|
| **Bernoulli** | p (success probability) | Binary outcomes |
| **Binomial** | n (trials), p (success prob) | Number of successes in n trials |
| **Poisson** | λ (rate) | Count data, rare events |
| **Geometric** | p | Number of failures before the first success (k = 0, 1, 2, ...) |
| **Uniform Discrete** | a, b | Discrete uniform outcomes |

The Binomial `Median` is the smallest k with CDF(k) ≥ 1/2. It always lies in [⌊np⌋, ⌈np⌉] and equals np when np is an integer [[5]](#5), so `Median` returns ⌊np⌋ when CDF(⌊np⌋) ≥ 1/2 and ⌈np⌉ otherwise. For p = 1/2 and odd n the distribution is symmetric, CDF((n − 1)/2) is exactly 1/2, and the median is (n − 1)/2.

The probability mass functions of the Binomial, Geometric, Poisson and Uniform Discrete distributions (returned by `PDF`) are zero away from the integers, and their CDFs floor a non-integer argument, so `CDF(2.5)` equals `CDF(2)`.

`Poisson` and `Binomial` also have static members that take one parameter set per element, for count models whose rate or success probability changes from element to element. They evaluate the same private kernels as the instance members, so every element is bitwise equal to the instance result, without creating a distribution per element. The Poisson quantile, instance and static alike, starts its search at the Cornish-Fisher estimate λ + √λ z + (z² − 1)/6 rather than at ⌊λ⌋, so it usually needs two CDF evaluations, and the static arrays of 4096 or more elements are evaluated in parallel ranges; neither changes a result:

```cs
double[] rates = { 0.8, 2.5, 12.0 };
double[] p = { 0.5, 0.9, 0.99 };
double[] counts = Poisson.InverseCDF(rates, p);                  // { 1, 5, 21 }
double[] logMass = Poisson.LogPDF(rates, new double[] { 0, 3, 12 });
double[] cdf = Poisson.CDF(rates, counts);

double[] binomialLogMass = Binomial.LogPDF(new[] { 0.3, 0.5 }, new[] { 10, 4 }, new double[] { 3, 2 });
double[] binomialCdf = Binomial.CDF(new[] { 0.3, 0.5 }, new[] { 10, 4 }, new double[] { 3, 2 });
```

The arrays must have equal lengths, and each rate must be positive and finite (each success probability in [0, 1] and each number of trials positive). The binomial log mass forms ln C(n, k) from log factorials, so its absolute error grows with ln n!: against mpmath it is 5.8E-13 at n = 500 and 1.3E-11 at n = 20000, and the binomial CDF is within 8E-13 relative at n = 20000.

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
double InverseCCDF(double q) // Survival quantile: x with CCDF(x) = q
double HF(double x)          // Hazard function
double LogPDF(double x)      // ln(PDF(x))
double LogCDF(double x)      // ln(CDF(x))
double LogCCDF(double x)     // ln(CCDF(x))
```

Each of `PDF`, `CDF`, `InverseCDF`, `CCDF`, `InverseCCDF`, `LogPDF`, `LogCDF` and `LogCCDF` also has a list form that takes an `IList<double>` and returns a `double[]`, each entry equal to the scalar function at that value. One call evaluates the whole list, so a caller crossing a language boundary (for example from Python) pays the crossing once:

```cs
var normal = new Normal(100, 15);
double[] exceedance = normal.CCDF(new[] { 100.0, 130.0, 160.0 });      // { 0.5, 0.02275, 3.167E-05 }
double[] levels = normal.InverseCCDF(new[] { 0.1, 0.01, 0.001 });      // { 119.22, 134.90, 146.35 }
double[] logDensity = normal.LogPDF(new[] { 85.0, 100.0, 115.0 });
```

`CCDF` equals `1 - CDF` in exact arithmetic. Families with a complementary special function (for example Beta, Cauchy, Chi-Squared, Exponential, Gamma, Inverse Chi-Squared and Inverse Gamma) evaluate the upper tail directly, so exceedance probabilities far below machine epsilon keep their relative accuracy instead of rounding to zero.

`InverseCCDF(q)` is the inverse in survival space, the value exceeded with probability q. By default it returns `InverseCDF(1 - q)` bit for bit, so forming 1 - q rounds away the low digits of a small q and every q below about 1.1E-16 maps to the upper support limit. `GeneralizedPareto` and `GeneralizedExtremeValue` override it with closed forms that never form 1 - q (see [Generalized Pareto Distribution](#generalized-pareto-distribution-gpd) and [Generalized Extreme Value](#generalized-extreme-value-gev)), and `EmpiricalDistribution` reads a descending (survival) table at q directly (see [Empirical Distribution](#empirical-distribution)). A list overload `InverseCCDF(IList<double>)` evaluates many probabilities in one call.

### Random Generation
```cs
double[] GenerateRandomValues(int sampleSize, int seed = -1)
```

## Creating Distributions

### Method 1: Direct Construction with Parameters

```cs
using Numerics.Distributions;

// Normal distribution: N(100, 15)
var normal = new Normal(mean: 100, standardDeviation: 15);

// Generalized Extreme Value: GEV(1000, 200, -0.1)
var gev = new GeneralizedExtremeValue(location: 1000, scale: 200, shape: -0.1);

// Log-Normal distribution
var lognormal = new LogNormal(meanOfLog: 4.5, standardDeviationOfLog: 0.5);

// Gamma distribution
var gamma = new GammaDistribution(scale: 5, shape: 2);
```

### Method 2: Using SetParameters

```cs
// Create with default parameters, then set
var weibull = new Weibull();
weibull.SetParameters(new double[] { 50, 2.5 }); // lambda=50, kappa=2.5

// Or use named parameters
weibull.SetParameters(scale: 50, shape: 2.5);
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
var gev = new GeneralizedExtremeValue(location: 1000, scale: 200, shape: -0.1);

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
var weibull = new Weibull(scale: 100, shape: 2.5);

// Hazard at time t=50
double hazard = weibull.HF(50);
Console.WriteLine($"Hazard rate at t=50: {hazard:F6}");

// For Weibull, hazard increases with time when κ > 1 (wear-out)
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

### Quantiles for a Parameter Ensemble

`ParameterSetQuantiles` evaluates the quantile function for every row of a parameter matrix, for example the posterior or bootstrap parameter sets of an uncertainty analysis, in one call. Each row is applied to a private clone, so the calling instance is not modified, and the results are bitwise identical to calling `SetParameters` and `InverseCDF` row by row. Rows that fail parameter validation are flagged in `valid` and their quantiles are `NaN`.

```cs
var gev = new GeneralizedExtremeValue(100, 20, -0.1);

// One parameter set (ξ, α, κ) per row
double[,] parameterSets =
{
    { 100, 20, -0.10 },
    { 105, 22, -0.05 },
    { 100, -1,  0.00 },   // invalid scale
};
double[] probabilities = { 0.5, 0.9, 0.99 };

double[,] quantiles = gev.ParameterSetQuantiles(parameterSets, probabilities, out bool[] valid);

for (int i = 0; i < parameterSets.GetLength(0); i++)
    Console.WriteLine($"Row {i}: valid = {valid[i]}, Q99 = {quantiles[i, 2]:F2}");
```

## Distribution Mathematical Definitions

This section provides the mathematical foundations for the most commonly used distributions in ***Numerics***. Understanding the underlying probability density functions (PDF), cumulative distribution functions (CDF), and moment structure is essential for correct application, particularly in safety-critical engineering analysis.

### Normal Distribution

The Normal (Gaussian) distribution is the most widely used continuous distribution. It arises naturally from the Central Limit Theorem, which states that the sum of many independent random variables tends toward a Normal distribution regardless of the underlying distributions. In ***Numerics***, the Normal distribution is parameterized by its mean $\mu$ and standard deviation $\sigma$.

#### Mathematical Definition

**Probability Density Function (PDF):**

```math
f(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\!\left(-\frac{(x - \mu)^2}{2\sigma^2}\right), \quad -\infty < x < \infty
```

**Cumulative Distribution Function (CDF):**

```math
F(x) = \Phi\!\left(\frac{x - \mu}{\sigma}\right) = \frac{1}{2}\left[1 + \text{erf}\!\left(\frac{x - \mu}{\sigma\sqrt{2}}\right)\right]
```

where $\Phi(\cdot)$ is the standard Normal CDF and $\text{erf}(\cdot)$ is the error function.

**Moments:**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \mu$ |
| Variance | $\text{Var}(X) = \sigma^2$ |
| Skewness | $\gamma_1 = 0$ |
| Kurtosis | $\kappa = 3$ |

**Parameters in Numerics:** `Normal(mean, standardDeviation)` where `mean` = $\mu$ and `standardDeviation` = $\sigma > 0$.

```cs
using Numerics.Distributions;

var normal = new Normal(mean: 100, standardDeviation: 15);

double pdf = normal.PDF(110);        // f(110)
double cdf = normal.CDF(110);        // P(X <= 110)
double q95 = normal.InverseCDF(0.95); // 95th percentile

Console.WriteLine($"Mean: {normal.Mean}");           // 100
Console.WriteLine($"Std Dev: {normal.StandardDeviation}"); // 15
Console.WriteLine($"Skewness: {normal.Skewness}");   // 0
```

**Estimation:** The maximum likelihood estimates have the closed form $\hat{\mu} = \bar{x}$ and $\hat{\sigma} = \sqrt{\tfrac{1}{n}\sum_{i}(x_i - \bar{x})^2}$, the population standard deviation with divisor $n$. The method of moments uses the sample standard deviation with divisor $n - 1$.

**Conditional moments:** `ConditionalMoments(a, b)` returns the conditional mean $E[X \mid a \le X \le b]$ followed by the second, third and fourth conditional moments about $\mu$, $E[(X-\mu)^k \mid a \le X \le b]$. They are closed-form truncated Normal moments in $\alpha = (a-\mu)/\sigma$ and $\beta = (b-\mu)/\sigma$. The boundary terms $\alpha^k\phi(\alpha)$ and $\beta^k\phi(\beta)$ are zero at infinite limits, and an interval lying entirely above or below the mean takes its mass from tail-probability differences, $\Phi(-\alpha) - \Phi(-\beta)$ or $\Phi(\beta) - \Phi(\alpha)$, so far-tail intervals do not lose the mass to cancellation.

```cs
double[] sample = { 9.8, 10.4, 10.1, 9.6, 10.9, 10.2 };

var normal = new Normal();
normal.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
Console.WriteLine($"μ = {normal.Mu:F4}, σ = {normal.Sigma:F4}"); // sample mean, divisor-n standard deviation

// Conditional mean of a standard Normal above 1: φ(1)/(1 − Φ(1))
double[] cm = new Normal(0, 1).ConditionalMoments(1.0, double.PositiveInfinity);
Console.WriteLine($"E[X | X ≥ 1] = {cm[0]:F6}"); // 1.525135
```

### Log-Normal Distribution

A random variable $X$ follows a Log-Normal distribution if $\log(X)$ follows a Normal distribution. The Log-Normal distribution is appropriate for strictly positive, right-skewed data arising from multiplicative processes. In ***Numerics***, the Log-Normal is parameterized by the mean and standard deviation of the log-transformed data, and uses base-10 logarithms by default. The `Base` property can be changed to use natural logarithms or any other base.

#### Mathematical Definition

For a general logarithmic base $b$, let $K = 1/\ln(b)$. The PDF and CDF are expressed in terms of $\log_b(x)$.

**Probability Density Function (PDF):**

```math
f(x) = \frac{K}{x \sigma \sqrt{2\pi}} \exp\!\left(-\frac{(\log_b x - \mu)^2}{2\sigma^2}\right), \quad x > 0
```

where $\mu$ and $\sigma$ are the mean and standard deviation of $\log_b(X)$.

**Cumulative Distribution Function (CDF):**

```math
F(x) = \frac{1}{2}\left[1 + \text{erf}\!\left(\frac{\log_b x - \mu}{\sigma\sqrt{2}}\right)\right]
```

**Moments (general base $b$, where $\beta = \ln(b)$):**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \exp\!\left[(\mu + \tfrac{1}{2}\sigma^2 \beta)\,\beta\right]$ |
| Mode | $\exp(\mu / K) = b^{\mu}$ |

**Key relationship:** If $X \sim \text{LogNormal}(\mu, \sigma)$, then $\log_b(X) \sim \text{Normal}(\mu, \sigma)$.

**Parameters in Numerics:** `LogNormal(meanOfLog, standardDeviationOfLog)` where `meanOfLog` = $\mu$ and `standardDeviationOfLog` = $\sigma > 0$. Default base is 10.

```cs
// Log-Normal with base-10 parameters
var lognormal = new LogNormal(meanOfLog: 4.5, standardDeviationOfLog: 0.5);

// Default base is 10
Console.WriteLine($"Base: {lognormal.Base}"); // 10

// Switch to natural log if needed
lognormal.Base = Math.E;

// Key relationship: log(X) ~ Normal
double x = lognormal.InverseCDF(0.5);
Console.WriteLine($"Median: {x:F2}");
Console.WriteLine($"log10(Median): {Math.Log10(x):F2}"); // equals Mu when Base=10
```

### Gamma Distribution

The Gamma distribution is a flexible two-parameter family for modeling positive-valued random variables. It generalizes the Exponential distribution and appears frequently in waiting-time problems, rainfall modeling, and Bayesian statistics. The ***Numerics*** library uses the shape/scale parameterization.

#### Mathematical Definition

**Probability Density Function (PDF):**

```math
f(x) = \frac{x^{\kappa-1}\, e^{-x/\theta}}{\theta^{\kappa}\,\Gamma(\kappa)}, \quad x > 0
```

where $\theta > 0$ is the scale parameter and $\kappa > 0$ is the shape parameter. $\Gamma(\cdot)$ is the gamma function.

**Cumulative Distribution Function (CDF):**

```math
F(x) = \frac{\gamma(\kappa,\, x/\theta)}{\Gamma(\kappa)} = P(\kappa,\, x/\theta)
```

where $\gamma(\kappa, z)$ is the lower incomplete gamma function and $P(\kappa, z)$ is the regularized lower incomplete gamma function.

`CDF` and `CCDF` evaluate $P(\kappa, x/\theta)$ and $Q(\kappa, x/\theta)$ with `Gamma.LowerIncomplete` and `Gamma.UpperIncomplete`, and `LogCDF` and `LogCCDF` use the same incomplete gamma core in logarithmic form (the factor $z^\kappa e^{-z}/\Gamma(\kappa)$ carried as its logarithm), so tails below the double range keep finite logarithms. `ChiSquared`, `InverseChiSquared`, `InverseGamma`, `PearsonTypeIII`, `LogPearsonTypeIII` and the large-rate `Poisson` tails share this log core, and `InverseCDF` solves on it. Against the mpmath oracle in `docs/mathematics/oracles/special-functions` (shapes from $10^{-3}$ to $10^6$) the relative error of `CDF` and `CCDF` is at most about $10^{-15}$ times $\max(1, |\ln F|)$; for $\kappa = 21544$ and $x = 1.25\kappa$, `CCDF` = 5.64E-254 is within $1.1 \times 10^{-13}$ of mpmath. Shapes up to $4 \times 10^{12}$ agree with R's `pgamma` to about $10^{-15}$ near the median.

**Moments:**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \kappa\theta$ |
| Variance | $\text{Var}(X) = \kappa\theta^2$ |
| Skewness | $\gamma_1 = 2/\sqrt{\kappa}$ |
| Kurtosis | $\kappa_4 = 3 + 6/\kappa$ |
| Mode | $(\kappa - 1)\theta$ for $\kappa \geq 1$, and 0 for $\kappa < 1$ |

**Special cases:** When $\kappa$ is a positive integer, the Gamma distribution is also known as the Erlang distribution.

**Parameters in Numerics:** `GammaDistribution(scale, shape)` where `scale` = $\theta > 0$ and `shape` = $\kappa > 0$.

```cs
var gamma = new GammaDistribution(scale: 5, shape: 2);

Console.WriteLine($"Mean: {gamma.Mean}");           // κθ = 10
Console.WriteLine($"Variance: {gamma.Variance}");   // κθ² = 50
Console.WriteLine($"Skewness: {gamma.Skewness:F4}"); // 2/√κ
Console.WriteLine($"Rate (1/θ): {gamma.Rate}");      // 0.2

double cdf = gamma.CDF(15.0);
double quantile = gamma.InverseCDF(0.95);
```

**L-moments:** `LinearMomentsFromParameters` gives $\lambda_1 = \kappa\theta$ and $\lambda_2 = \theta\,\Gamma(\kappa + \tfrac{1}{2})/(\sqrt{\pi}\,\Gamma(\kappa))$, with the gamma ratio formed from $\ln B(\kappa, \tfrac{1}{2})$ so that it stays finite for large shapes (small skews), where $\Gamma(\kappa)$ itself overflows beyond $\kappa \approx 171$. Against mpmath the relative error of $\lambda_2$ is below $10^{-15}$ for shapes from 0.01 to $10^6$. $\tau_3$ and $\tau_4$ use Hosking's rational-function approximations, accurate to about $10^{-6}$.

### Exponential Distribution

The Exponential distribution models the time between events in a Poisson process. Its defining property is *memorylessness*: the probability of an event occurring in the next $\Delta t$ time units is independent of how much time has already elapsed. In ***Numerics***, the Exponential distribution is a two-parameter (shifted) distribution with location $\xi$ and scale $\alpha$. Setting $\xi = 0$ yields the standard one-parameter Exponential.

#### Mathematical Definition

**Probability Density Function (PDF):**

```math
f(x) = \frac{1}{\alpha}\exp\!\left(-\frac{x - \xi}{\alpha}\right), \quad x \geq \xi
```

**Cumulative Distribution Function (CDF):**

```math
F(x) = 1 - \exp\!\left(-\frac{x - \xi}{\alpha}\right), \quad x \geq \xi
```

**Inverse CDF (Quantile Function):**

```math
Q(p) = \xi - \alpha \ln(1 - p)
```

**Moments:**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \xi + \alpha$ |
| Variance | $\text{Var}(X) = \alpha^2$ |
| Skewness | $\gamma_1 = 2$ |
| Kurtosis | $\kappa = 9$ |
| Mode | $\xi$ |

**Relationship to Gamma:** The Exponential($\xi$, $\alpha$) distribution is a special case of a shifted Gamma distribution with shape $\kappa = 1$.

**Parameters in Numerics:** `Exponential(location, scale)` where `location` = $\xi$ and `scale` = $\alpha > 0$. A single-parameter constructor `Exponential(scale)` sets $\xi = 0$.

```cs
// Two-parameter (shifted) exponential
var exp2 = new Exponential(location: 10, scale: 5);
Console.WriteLine($"Mean: {exp2.Mean}");   // ξ + α = 15
Console.WriteLine($"Mode: {exp2.Mode}");   // ξ = 10

// One-parameter (standard) exponential
var exp1 = new Exponential(scale: 5);
Console.WriteLine($"Mean: {exp1.Mean}");   // α = 5
Console.WriteLine($"P(X > 10): {exp1.CCDF(10):F4}"); // memoryless property
```

### Gumbel Distribution (Extreme Value Type I)

The Gumbel distribution is a special case of the Generalized Extreme Value (GEV) distribution with shape parameter $\kappa = 0$. It models the distribution of the maximum (or minimum) of a sample drawn from various distributions with exponentially decaying tails. It is widely used in hydrology for modeling annual maximum floods and in structural engineering for modeling extreme wind loads.

#### Mathematical Definition

**Probability Density Function (PDF):**

```math
f(x) = \frac{1}{\alpha}\exp\!\left[-(z + e^{-z})\right], \quad z = \frac{x - \xi}{\alpha}, \quad -\infty < x < \infty
```

**Cumulative Distribution Function (CDF):**

```math
F(x) = \exp\!\left(-e^{-z}\right), \quad z = \frac{x - \xi}{\alpha}
```

**Inverse CDF (Quantile Function):**

```math
Q(p) = \xi - \alpha \ln(-\ln p)
```

**Moments:**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \xi + \alpha\gamma_E$ where $\gamma_E \approx 0.5772$ is the Euler-Mascheroni constant |
| Variance | $\text{Var}(X) = \frac{\pi^2}{6}\alpha^2$ |
| Skewness | $\gamma_1 = 12\sqrt{6}\,\zeta(3)/\pi^3 \approx 1.13955$ |
| Kurtosis | $\kappa = 3 + 12/5 = 5.4$ |
| Mode | $\xi$ |
| Median | $\xi - \alpha\ln(\ln 2)$ |

**Parameters in Numerics:** `Gumbel(location, scale)` where `location` = $\xi$ and `scale` = $\alpha > 0$.

```cs
var gumbel = new Gumbel(location: 100, scale: 25);

Console.WriteLine($"Mode: {gumbel.Mode}");      // ξ = 100
Console.WriteLine($"Mean: {gumbel.Mean:F2}");    // ξ + αγ ≈ 114.43
Console.WriteLine($"Median: {gumbel.Median:F2}");

// 100-year return period quantile
double q100 = gumbel.InverseCDF(0.99);
Console.WriteLine($"1% AEP quantile: {q100:F2}");
```

### Uniform Distribution

The Uniform distribution assigns equal probability density to all values within a bounded interval $[a, b]$. It represents maximum uncertainty (maximum entropy) given only knowledge of the support bounds. It is commonly used as a non-informative prior in Bayesian analysis and for random number generation.

#### Mathematical Definition

**Probability Density Function (PDF):**

```math
f(x) = \frac{1}{b - a}, \quad a \leq x \leq b
```

**Cumulative Distribution Function (CDF):**

```math
F(x) = \frac{x - a}{b - a}, \quad a \leq x \leq b
```

**Inverse CDF (Quantile Function):**

```math
Q(p) = a + p(b - a)
```

**Moments:**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \frac{a + b}{2}$ |
| Variance | $\text{Var}(X) = \frac{(b - a)^2}{12}$ |
| Skewness | $\gamma_1 = 0$ |
| Kurtosis | $\kappa = 9/5 = 1.8$ |

**Parameters in Numerics:** `Uniform(min, max)` where `min` = $a$ and `max` = $b \geq a$.

```cs
var uniform = new Uniform(min: 0, max: 10);

Console.WriteLine($"Mean: {uniform.Mean}");     // 5
Console.WriteLine($"Std Dev: {uniform.StandardDeviation:F4}"); // 10/√12
Console.WriteLine($"PDF(5): {uniform.PDF(5)}"); // 0.1 everywhere in [0,10]

double median = uniform.InverseCDF(0.5); // 5
```

Setting `min` equal to `max` gives a point mass at that value: the CDF steps from 0 to 1 at `Min`, `InverseCDF` returns `Min` for every probability, and the PDF is 0 because a point mass has no density. `PDF(double.NaN)` returns NaN.

**Estimation:** Maximum likelihood gives the sample minimum and maximum. The method of moments gives $\bar{x} \mp \sqrt{3}\,s$, where $s$ is the sample standard deviation with divisor $n - 1$.

```cs
double[] obs = { 2.3, 7.9, 4.1, 5.6, 3.3, 6.8 };

var mle = new Uniform();
mle.Estimate(obs, ParameterEstimationMethod.MaximumLikelihood); // (min, max) = (2.3, 7.9)
var mom = new Uniform();
mom.Estimate(obs, ParameterEstimationMethod.MethodOfMoments);   // mean ∓ √3·s
Console.WriteLine($"MLE [{mle.Min}, {mle.Max}], MoM [{mom.Min:F4}, {mom.Max:F4}]");

var point = new Uniform(5, 5); // point mass at 5
Console.WriteLine($"CDF(4.9) = {point.CDF(4.9)}, CDF(5) = {point.CDF(5)}, Q(0.3) = {point.InverseCDF(0.3)}");
```

### Triangular Distribution

The Triangular distribution is defined by three parameters: minimum $a$, mode $c$, and maximum $b$. It provides a simple model for uncertainty when only the range and most likely value are known. It is frequently used in risk assessment, project management (PERT analysis), and expert elicitation.

#### Mathematical Definition

**Probability Density Function (PDF):**

```math
f(x) = \begin{cases}
\dfrac{2(x - a)}{(b - a)(c - a)} & a \leq x < c \\[6pt]
\dfrac{2}{b - a} & x = c \\[6pt]
\dfrac{2(b - x)}{(b - a)(b - c)} & c < x \leq b
\end{cases}
```

**Cumulative Distribution Function (CDF):**

```math
F(x) = \begin{cases}
\dfrac{(x - a)^2}{(b - a)(c - a)} & a \leq x \leq c \\[6pt]
1 - \dfrac{(b - x)^2}{(b - a)(b - c)} & c < x \leq b
\end{cases}
```

**Moments:**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \frac{a + b + c}{3}$ |
| Variance | $\text{Var}(X) = \frac{a^2 + b^2 + c^2 - ab - ac - bc}{18}$ |
| Mode | $c$ |

**Parameters in Numerics:** `Triangular(min, mode, max)` where `min` = $a$, `mode` = $c$, and `max` = $b$, with $a \leq c \leq b$.

```cs
var tri = new Triangular(min: 10, mode: 15, max: 25);

Console.WriteLine($"Mean: {tri.Mean:F2}");      // (10+25+15)/3 = 16.67
Console.WriteLine($"Mode: {tri.Mode}");          // 15
Console.WriteLine($"Median: {tri.Median:F2}");
Console.WriteLine($"Std Dev: {tri.StandardDeviation:F2}");
```

**Estimation:** Maximum likelihood bounds the minimum below by $-10^{\lceil \log_{10}|\min x| \rceil + 1}$ and the maximum above by $10^{\lceil \log_{10}|\max x| \rceil + 1}$. When that envelope cannot contain the starting values, for example when the sample minimum or maximum is exactly 0, the bounds are $\pm 10^{\lceil \log_{10} S \rceil + 1}$ with $S = \max(|\min x|, |\max x|, \max x - \min x)$, so such samples can still be fitted. The PERT distribution uses the same bounds.

### Beta Distribution

The Beta distribution is defined on the interval $[0, 1]$ and is parameterized by two positive shape parameters $\alpha$ and $\beta$. Its flexibility makes it ideal for modeling random proportions, probabilities, and percentages. In Bayesian statistics, it serves as the conjugate prior for the Bernoulli and Binomial distributions.

#### Mathematical Definition

**Probability Density Function (PDF):**

```math
f(x) = \frac{x^{\alpha-1}(1-x)^{\beta-1}}{B(\alpha, \beta)}, \quad 0 \leq x \leq 1
```

where $B(\alpha, \beta)$ is the Beta function:

```math
B(\alpha, \beta) = \frac{\Gamma(\alpha)\,\Gamma(\beta)}{\Gamma(\alpha + \beta)}
```

**Cumulative Distribution Function (CDF):**

```math
F(x) = I_x(\alpha, \beta)
```

where $I_x(\alpha, \beta)$ is the regularized incomplete Beta function.

**Moments:**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \frac{\alpha}{\alpha + \beta}$ |
| Variance | $\text{Var}(X) = \frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}$ |
| Mode | $\frac{\alpha - 1}{\alpha + \beta - 2}$ for $\alpha, \beta > 1$; when both shapes are below 1 the density is unbounded at 0 and 1, and `Mode` returns 0.5, where the density has its minimum |
| Skewness | $\frac{2(\beta - \alpha)\sqrt{\alpha + \beta + 1}}{(\alpha + \beta + 2)\sqrt{\alpha\beta}}$ |

**Special cases:** $\text{Beta}(1, 1) = \text{Uniform}(0, 1)$.

**Parameters in Numerics:** `BetaDistribution(alpha, beta)` where `alpha` = $\alpha > 0$ and `beta` = $\beta > 0$.

```cs
var beta = new BetaDistribution(alpha: 2, beta: 5);

Console.WriteLine($"Mean: {beta.Mean:F4}");      // α/(α+β) = 2/7
Console.WriteLine($"Mode: {beta.Mode:F4}");      // (α-1)/(α+β-2) = 1/5
Console.WriteLine($"Variance: {beta.Variance:F4}");

// Beta(1,1) is equivalent to Uniform(0,1)
var betaUniform = new BetaDistribution(alpha: 1, beta: 1);
Console.WriteLine($"Beta(1,1) PDF(0.5): {betaUniform.PDF(0.5)}"); // 1.0
```

### Student's t Distribution

The Student's t distribution in ***Numerics*** is the location-scale form with location $\mu$, scale $\sigma$ and degrees of freedom $\nu$. It has heavier tails than the Normal distribution and approaches it as $\nu \to \infty$.

#### Mathematical Definition

**PDF:**

```math
f(x) = \frac{\Gamma\!\left(\frac{\nu+1}{2}\right)}{\sigma\sqrt{\nu\pi}\,\Gamma\!\left(\frac{\nu}{2}\right)}\left(1 + \frac{z^2}{\nu}\right)^{-(\nu+1)/2}, \quad z = \frac{x - \mu}{\sigma}
```

**CDF:** With $I_y(a, b)$ the regularized incomplete beta function, the lower tail ($z < 0$, $z^2 \ge 1/2$) is evaluated in the cancellation-free form

```math
F(z) = I_{w}\!\left(\tfrac{\nu}{2}, \tfrac{\nu}{2}\right), \quad w = \frac{1}{2q(s + q)}, \quad s = \frac{|z|}{\sqrt{\nu}}, \quad q = \sqrt{1 + s^2}
```

near the center ($z^2 < 1/2$) as $F(z) = \tfrac{1}{2} - \tfrac{1}{2} I_{z^2/(\nu+z^2)}\!\left(\tfrac{1}{2}, \tfrac{\nu}{2}\right)$, and the upper tail by symmetry. `CCDF`, `LogCDF`, `LogCCDF` and `LogPDF` are evaluated directly, so the log functions stay finite after the probabilities underflow. For $\nu \ge 30$ and $\nu/(\nu + z^2) \ge 0.7$ the tail functions use the asymptotic expansion of the incomplete beta function in the large shape $\nu/2$ (routine BGRAT of TOMS algorithm 708), whose leading term is $\operatorname{erfc}\sqrt{(\nu/2 - 1/4)\ln(1 + z^2/\nu)}$, so they stay accurate for every $\nu$ without a Normal approximation. On the 97-point mpmath grid of the unit tests ($\nu$ from 1 to $10^{10}$, $|z|$ from $10^{-6}$ to $10^6$), the two-sided tail $2\,\text{CCDF}(|z|)$ is within 3.4E-15 relative where it is at least 0.01, 1.1E-14 down to $10^{-10}$, and 8E-13 in the far tail; further test points check $\nu$ up to $10^{30}$ to 1E-12. `InverseCDF` inverts the incomplete beta function of the two-sided centre probability $|1 - 2p|$ for $1/4 < p < 3/4$, which keeps the relative accuracy of $t$ near the median (for $\nu \ge 10^{16}$ the centre quantile is the Normal quantile to double precision); in the tails with $\nu \ge 30$ it solves $\ln F(t) = \ln p$ by Brent's method, and below $\nu = 30$ it inverts the incomplete beta function of $2\min(p, 1 - p)$, except where that probability is subnormal (the log CDF is solved) and for $\nu = 1$ below $10^{-150}$ (the Cauchy quantile $\mp 1/(\pi \min(p, 1 - p))$). For $\nu \ge 10^8$ the `PDF` uses the Normal distribution. `CDF(+∞)` = 1, `CDF(−∞)` = 0, and `InverseCDF` throws for a probability outside [0, 1], including NaN.

**Estimation:** Maximum likelihood fits $\mu$, $\sigma$ and $\nu$ jointly, with $\nu$ bounded to [1, 1000]; the method of moments is not supported. The class implements `IBootstrappable`.

**Parameters in Numerics:** `StudentT(location, scale, degreesOfFreedom)`.

```cs
var t = new StudentT(location: 0, scale: 1, degreesOfFreedom: 5);

Console.WriteLine($"CDF(-2) = {t.CDF(-2.0):F9}");     // 0.050969739
Console.WriteLine($"LogCDF(-1E10) = {t.LogCDF(-1E10):F4}"); // finite after the CDF underflows

double[] data = { -1.2, 0.4, 2.9, -0.3, 0.8, -4.1, 0.1, 1.5, -0.6, 0.2 };
var fitted = new StudentT();
fitted.Estimate(data, ParameterEstimationMethod.MaximumLikelihood);
Console.WriteLine($"μ = {fitted.Mu:F3}, σ = {fitted.Sigma:F3}, ν = {fitted.DegreesOfFreedom:F2}");
```

### Noncentral t Distribution

The noncentral t distribution is the distribution of $T = (Z + \mu)/\sqrt{V/\nu}$, where $Z$ is standard Normal and $V$ is an independent chi-squared variable with $\nu$ degrees of freedom. In ***Numerics*** the noncentrality parameter is named $\mu$ (property `Noncentrality`). It is used in power analysis and in confidence limits for Normal quantiles.

#### Mathematical Definition

**CDF:** Let $\lambda = \mu^2/2$, $p_j = e^{-\lambda}\lambda^j/j!$, $q_j = \mu e^{-\lambda}\lambda^j / (\sqrt{2}\,\Gamma(j + 3/2))$ and $y = x^2/(x^2 + \nu)$. For $x > 0$ the CDF is the Poisson-weighted incomplete beta series [[6]](#6), [[7]](#7)

```math
F(x) = \Phi(-\mu) + \frac{1}{2}\sum_{j=0}^{\infty}\left[p_j\, I_y\!\left(j + \tfrac{1}{2}, \tfrac{\nu}{2}\right) + q_j\, I_y\!\left(j + 1, \tfrac{\nu}{2}\right)\right]
```

The series starts at its modal term $j = \lfloor\lambda\rfloor$ with the weights evaluated in log space and is summed in both directions [[8]](#8), so it keeps its mass for any noncentrality. Summation stops when a bound on the remaining terms falls below machine precision relative to the sum. For $x < 0$ and $\mu \le 0$ the CDF is summed as the complementary series, all of whose terms are positive. For $x < 0$ and $\mu > 0$ the CDF is the positive integral

```math
F(x) = \int_0^\infty \phi(\mu + u)\, P\!\left(\tfrac{\nu}{2}, \tfrac{\nu u^2}{2x^2}\right) du
```

with $P$ the regularized lower incomplete gamma function, evaluated by adaptive Gauss-Kronrod quadrature. The CCDF uses the reflection $P(T > x) = F(-x; \nu, -\mu)$.

**PDF:** With $S = \sqrt{\chi^2_\nu/\nu}$, the density is the positive integral $f(x) = \int_0^\infty s\,\phi(xs - \mu) f_S(s)\, ds$, evaluated by adaptive Gauss-Kronrod quadrature. At $x = 0$ the closed form $\Gamma(\tfrac{\nu+1}{2}) e^{-\mu^2/2} / (\sqrt{\pi\nu}\,\Gamma(\tfrac{\nu}{2}))$ is used, through log-gamma functions when $\nu$ is large.

**Quantile:** `InverseCDF` brackets the quantile and refines it with Brent's method. The CDF, CCDF and PDF are accurate to about $10^{-12}$ relative.

**Estimation:** Maximum likelihood fits $\nu$ and $\mu$ jointly. The class implements `IBootstrappable`.

**Parameters in Numerics:** `NoncentralT(degreesOfFreedom, noncentrality)`.

```cs
var nct = new NoncentralT(degreesOfFreedom: 10, noncentrality: 2);

Console.WriteLine($"CDF(1.5)  = {nct.CDF(1.5):F10}");
Console.WriteLine($"CCDF(8)   = {nct.CCDF(8.0):E6}");   // upper tail without 1 − CDF cancellation
Console.WriteLine($"Q(0.95)   = {nct.InverseCDF(0.95):F6}");
```

### Pareto Distribution

The Pareto distribution with scale $x_m > 0$ and shape $\alpha > 0$ has support $[x_m, \infty)$ and a power-law upper tail.

#### Mathematical Definition

```math
f(x) = \frac{\alpha}{x}\left(\frac{x_m}{x}\right)^{\alpha}, \qquad F(x) = 1 - \left(\frac{x_m}{x}\right)^{\alpha}, \quad x \ge x_m
```

With $L = \ln(x/x_m)$, formed as $\text{log1p}((x - x_m)/x_m)$ near $x_m$, the CDF is evaluated as $-\text{expm1}(-\alpha L)$, the CCDF as $\exp(-\alpha L)$ and `LogCCDF` as $-\alpha L$, so both tails keep full relative precision.

**Estimation:** Maximum likelihood has the closed form $\hat{x}_m = \min x_i$ and $\hat{\alpha} = n / \sum_i \ln(x_i/\hat{x}_m)$. The class implements `IEstimation`, `IMaximumLikelihoodEstimation` and `IBootstrappable`.

**Parameters in Numerics:** `Pareto(scale, shape)`, with properties `Xm` and `Alpha`.

```cs
double[] incomes = { 1.2, 1.9, 1.4, 3.8, 1.1, 2.6, 1.7, 5.3 };

var pareto = new Pareto();
pareto.Estimate(incomes, ParameterEstimationMethod.MaximumLikelihood);
Console.WriteLine($"Xm = {pareto.Xm:F3}, α = {pareto.Alpha:F4}"); // Xm = min(x), α = n / Σ ln(x/Xm)
Console.WriteLine($"P(X > 10) = {pareto.CCDF(10.0):E4}");
```

### Poisson Distribution

The Poisson distribution with rate $\lambda > 0$ models counts of independent events.

#### Mathematical Definition

```math
P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}, \qquad F(k) = Q(\lfloor k \rfloor + 1, \lambda), \quad k = 0, 1, 2, \ldots
```

where $Q$ is the regularized upper incomplete gamma function. The mass is zero away from the nonnegative integers, the CDF floors $k$, and `CDF(+∞)` = 1. For $\lambda < 10^5$ the CDF uses `Gamma.UpperIncomplete`. For $\lambda \ge 10^5$ the mass is evaluated from the deviance form of the gamma kernel and the CDF from the uniform asymptotic (Temme) expansion of the incomplete gamma function [[10]](#10), which keep their relative accuracy at large rates. The CCDF and `LogCCDF` are evaluated from $P(\lfloor k \rfloor + 1, \lambda)$ in log space, so far upper-tail probabilities are not lost to $1 - F$. `InverseCDF` returns the smallest $k$ with $F(k) \ge p$, by a step-by-step walk below $\lambda = 10^5$ (from the Cornish-Fisher estimate for $p$ in $[10^{-300}, 1 - 10^{-10}]$, from $\lfloor\lambda\rfloor$ otherwise) and by bracketing and bisection above, and throws for a NaN probability.

**Estimation:** Maximum likelihood and the method of moments both give $\hat{\lambda} = \bar{x}$. The sample must consist of nonnegative integers with at least one positive count. The class implements `IEstimation`, `IMaximumLikelihoodEstimation`, `IMomentEstimation` and `IBootstrappable`.

**Parameters in Numerics:** `Poisson(rate)`, with property `Lambda`.

```cs
var poisson = new Poisson(rate: 4);

Console.WriteLine($"PMF(2)   = {poisson.PDF(2):F6}");   // 0.146525
Console.WriteLine($"PMF(2.5) = {poisson.PDF(2.5)}");    // 0, not an integer
Console.WriteLine($"CDF(2.5) = {poisson.CDF(2.5):F6}"); // equals CDF(2) = 0.238103
Console.WriteLine($"Median   = {poisson.InverseCDF(0.5)}"); // smallest k with CDF(k) ≥ 0.5

var large = new Poisson(rate: 1E6);
Console.WriteLine($"CDF(1E6) = {large.CDF(1E6):F10}");
```

### Rayleigh Distribution

The Rayleigh distribution with scale $\sigma > 0$ describes the magnitude of a two-dimensional vector with independent zero-mean Normal components.

#### Mathematical Definition

```math
f(x) = \frac{x}{\sigma^2}\, e^{-x^2/(2\sigma^2)}, \qquad F(x) = 1 - e^{-x^2/(2\sigma^2)}, \qquad Q(p) = \sigma\sqrt{-2\ln(1 - p)}, \quad x \ge 0
```

All functions are evaluated through $z = x/\sigma$, so $\sigma^2$ is never formed and large scales do not overflow. The CDF is evaluated as $-\text{expm1}(-z^2/2)$ and the quantile as $\sigma\sqrt{-2\,\text{log1p}(-p)}$; `CCDF`, `LogCCDF` and `LogPDF` are evaluated directly.

**Estimation:** The method of moments gives $\hat{\sigma} = \bar{x}/\sqrt{\pi/2}$. Maximum likelihood multiplies the closed-form estimate $\sqrt{\sum_i x_i^2/(2n)}$ by the bias correction $\Gamma(n)\sqrt{n}/\Gamma(n + \tfrac{1}{2})$, so the result is slightly larger than the plain maximum likelihood estimate and agrees with it as $n$ grows. Samples must be finite and nonnegative with at least one positive value. The class implements `IEstimation` and `IBootstrappable`.

```cs
var rayleigh = new Rayleigh(scale: 2.0);
Console.WriteLine($"CDF(3)  = {rayleigh.CDF(3.0):F6}");         // 1 − exp(−9/8) = 0.675348
Console.WriteLine($"Q(0.99) = {rayleigh.InverseCDF(0.99):F6}"); // 2·√(2 ln 100) = 6.069709

double[] heights = { 1.8, 2.6, 0.9, 3.1, 2.2, 1.4, 2.9, 1.7 };
var mom = new Rayleigh();
mom.Estimate(heights, ParameterEstimationMethod.MethodOfMoments);
var mle = new Rayleigh();
mle.Estimate(heights, ParameterEstimationMethod.MaximumLikelihood);
Console.WriteLine($"σ (MoM) = {mom.Sigma:F4}, σ (bias-corrected MLE) = {mle.Sigma:F4}");
```

### Uniform Discrete Distribution

The Uniform Discrete distribution assigns mass $1/N$ to each integer from `Min` through `Max`, with $N = \text{Max} - \text{Min} + 1$. Both bounds must be integers. `Min` = `Max` is a valid point mass, for which the kurtosis is undefined and returned as NaN.

The mass is zero away from the integers, the CDF $(\lfloor k \rfloor - \text{Min} + 1)/N$ floors a non-integer argument, and `InverseCDF` returns the smallest $k$ with $F(k) \ge p$ (`InverseCDF(0)` = `Min`).

**Estimation:** Maximum likelihood gives the sample minimum and maximum; every observation must be a finite integer. The class implements `IEstimation`, `IMaximumLikelihoodEstimation` and `IBootstrappable`.

```cs
var die = new UniformDiscrete(1, 6);

Console.WriteLine($"PMF(3)   = {die.PDF(3):F6}");     // 1/6
Console.WriteLine($"PMF(3.5) = {die.PDF(3.5)}");      // 0
Console.WriteLine($"CDF(3.5) = {die.CDF(3.5):F6}");   // CDF(3) = 0.5
Console.WriteLine($"Q(0.5)   = {die.InverseCDF(0.5)}"); // 3, the smallest k with CDF(k) ≥ 0.5
```

### Von Mises Distribution

The von Mises distribution is the circular analogue of the Normal distribution [[9]](#9), with mean direction $\mu \in [-\pi, \pi]$ and concentration $\kappa \ge 0$. It is used for directional data such as flood seasonality.

#### Mathematical Definition

```math
f(x) = \frac{e^{\kappa\cos(x - \mu)}}{2\pi I_0(\kappa)} = \frac{\exp\!\left(-2\kappa\sin^2\!\left(\frac{x - \mu}{2}\right)\right)}{2\pi\, I_{0e}(\kappa)}, \quad -\pi \le x \le \pi
```

where $I_{0e}(\kappa) = e^{-\kappa} I_0(\kappa)$ is `Bessel.I0e`. The half-angle form is finite for every $\kappa$. The CDF integrates the same kernel from $-\pi$ by adaptive Gauss-Kronrod quadrature, split at $\mu$ and $\mu \pm 8/\sqrt{\kappa}$ (for $\kappa \ge 1$) so the peak is resolved at large concentration; its absolute error is about $10^{-14}$.

**Summary statistics:** `Mean` and `Median` return the circular mean direction $\mu$. `StandardDeviation` returns $\sqrt{1 - A(\kappa)}$ with $A(\kappa) = I_1(\kappa)/I_0(\kappa)$ evaluated as `Bessel.I1e`/`Bessel.I0e`, the square root of the circular variance. It differs from the circular standard deviation $\sqrt{-2\ln \bar{R}}$ of Mardia and Jupp [[9]](#9) and from a linear standard deviation.

**Estimation:** Maximum likelihood sets $\hat{\mu}$ to the circular mean and solves $A(\hat{\kappa}) = \bar{R}$, the mean resultant length, by Newton-Raphson iterations. The sample needs at least two finite values; a sample with no angular dispersion ($\bar{R} = 1$) throws `ArgumentException`. The class implements `IEstimation`, `IMaximumLikelihoodEstimation` and `IBootstrappable`.

```cs
var vm = new VonMises(mu: 0, kappa: 1E4);

Console.WriteLine($"PDF(0)  = {vm.PDF(0):F6}");   // ≈ √(κ/2π) = 39.89
Console.WriteLine($"CDF(0)  = {vm.CDF(0):F12}");  // 0.5
Console.WriteLine($"SD      = {vm.StandardDeviation:F6}"); // √(1 − I₁e(κ)/I₀e(κ))

double[] angles = { 0.10, -0.25, 0.32, 0.05, -0.12, 0.41, -0.30, 0.18 };
var fitted = new VonMises();
fitted.Estimate(angles, ParameterEstimationMethod.MaximumLikelihood);
Console.WriteLine($"μ = {fitted.Mu:F4}, κ = {fitted.Kappa:F3}");
```

## Hydrological Distributions

### Log-Pearson Type III (LP3)

The LP3 distribution is the standard distribution for flood frequency analysis in the United States, as prescribed by Bulletin 17C [[1]](#1). It is constructed by applying the Pearson Type III distribution to the logarithms of the data. This means that if $X \sim \text{LP3}$, then $\log_b(X) \sim \text{Pearson Type III}$, where $b$ is the logarithmic base (default base 10 in ***Numerics***).

#### Mathematical Definition

The LP3 distribution is parameterized by the moments of the log-transformed data: $\mu$ (mean of log), $\sigma$ (standard deviation of log), and $\gamma$ (skewness of log). These parameters relate to the underlying Pearson Type III distribution through:

| LP3 Parameter | Pearson Type III Equivalent |
|---------------|---------------------------|
| Location $\xi$ | $\mu - 2\sigma/\gamma$ |
| Scale $\beta$ | $\sigma\gamma/2$ |
| Shape $\alpha$ | $4/\gamma^2$ |

The PDF and CDF of the LP3 are obtained by applying the Pearson Type III to the log-transformed variable. For a variable $Y = \log_b(X)$:

**Pearson Type III PDF (applied to Y):**

```math
f_Y(y) = \frac{|y - \xi|^{\alpha-1}\, \exp(-|y - \xi|/|\beta|)}{|\beta|^{\alpha}\,\Gamma(\alpha)}
```

where the sign conventions depend on the sign of $\gamma$: when $\gamma > 0$, $Y \geq \xi$; when $\gamma < 0$, $Y \leq \xi$.

**Quantile computation:** Quantiles follow the frequency factor form of Bulletin 17C:

```math
\log_b(X_p) = \mu + K_p \cdot \sigma
```

where $K_p$ is the Pearson Type III frequency factor, a function of the skewness $\gamma$ and probability $p$. `InverseCDF` evaluates $K_p$ exactly from the inverse of the gamma distribution. `WilsonHilfertyInverseCDF` uses the approximation `GammaDistribution.FrequencyFactorKp`: the Cornish-Fisher polynomial in $\gamma$ through $\gamma^6$ for $|\gamma| \le 2$, which is analytic through $\gamma = 0$ and equals the standard Normal quantile there, and the modified Wilson-Hilferty transformation above. Its fifth-order coefficient, $(12z^6 - 143z^4 - 923z^2 + 1472)/6531840$, is kept deliberately in place of the exact series term with $-243z^4$: against the actual gamma quantile its largest error is $1.4 \times 10^{-4}$, $1.2 \times 10^{-2}$ and $0.68$ for $|\gamma|$ up to 0.5, 1 and 2, versus $1.6 \times 10^{-4}$, $2.0 \times 10^{-2}$ and $0.93$ with the exact term, because it partly offsets the truncation after $\gamma^6$. Its published source has not been verified, so it is treated as empirically tuned. The Cornish-Fisher coefficients and the Wilson-Hilferty constants are attributed to Hoshi and Burges (1981); they are documented as implemented and have not been checked against Kirby (1972), Kite (1977), Hoshi and Burges (1981), PeakFQ, HEC-SSP or Bulletin 17B/17C reference code. `FrequencyFactorKp` is an approximation by design and does not return the gamma support limits at the probability endpoints: at $p = 1$ it is $+\infty$ for $|\gamma| < 10^{-4}$, not-a-number for $10^{-4} \le |\gamma| \le 2$, and for $|\gamma| > 2$ the Wilson-Hilferty value ($+\infty$ for positive skew, the clamped bound $-2/\gamma$ for negative skew); at $p = 0$ the Normal variate is finite (about $-38.4674$), so the approximation is evaluated there, and for $\gamma > 2$ the Wilson-Hilferty clamp gives $-2/\gamma$.

**Moments (of the log-transformed data):**

| Property | Formula |
|----------|---------|
| Mean of log | $\mu$ |
| Std dev of log | $\sigma$ |
| Skewness of log | $\gamma$ |

When $\gamma = 0$, the LP3 reduces to the Log-Normal distribution.

**Moments of the flows:** `Mean`, `StandardDeviation`, `Skewness` and `Kurtosis` describe $X$ itself. They follow from the moment-generating function of the Pearson Type III, $E[b^{rY}] = \exp\big(r\mu\ln b + \alpha[-\ln(1-u_r) - u_r]\big)$ with $u_r = r\sigma\gamma\ln(b)/2$, and moment $r$ exists only when $u_r < 1$. The closed-form central moments lose accuracy to cancellation as the log scale shrinks, so when $\sigma\ln(b)(1 + |\gamma|) < 0.15$ the skewness and kurtosis remove the cancelling terms analytically: with $d_r = K_r - rK_1$, where $K_r$ is the logarithm of the location-free moment, the normalized third central moment is $F(d_3) - 3F(d_2) + \Delta_3$ with $F(d) = e^d - 1 - d$ and the finite difference $\Delta_3 = K_3 - 3K_2 + 3K_1$ summed as a power series in $u = \sigma\gamma\ln(b)/2$ that starts at the first nonvanishing order, and the fourth is formed the same way. Against 60-digit references the skewness and kurtosis stay within about $1.1 \times 10^{-12}$ and $5 \times 10^{-13}$ relative, and the mean and standard deviation within about $2 \times 10^{-15}$.

**Parameters in Numerics:** `LogPearsonTypeIII(meanOfLog, standardDeviationOfLog, skewOfLog)`. Default logarithmic base is 10.

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

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

### Pearson Type III (P3)

The Pearson Type III distribution is a three-parameter gamma distribution parameterized by its mean $\mu$, standard deviation $\sigma$ and skewness $\gamma$; it is the distribution of the logarithms in the LP3 model above.

**Near-zero skew:** For $|\gamma| < 0.1$ the density and tails evaluate the gamma distribution at the unit coordinate $4/\gamma^2 + 2z/\gamma$, with $z = (x - \mu)/\sigma$, carrying the exact rounding residual of that coordinate and a first-order correction, so the results stay continuous as the skew approaches zero. Against 50-digit references the log density and log tails stay within about $5 \times 10^{-15}$ for probabilities from $10^{-6}$ to $1 - 10^{-6}$.

**Quantiles and their skew derivative:** `InverseCDF` inverts the gamma distribution as $\mu + \sigma(\gamma/2)(u - \alpha)$, with $u$ the unit gamma quantile and $\alpha = 4/\gamma^2$, which loses about $2\varepsilon/|\gamma|$ to cancellation for small skew, and the implicit skew derivative used by `QuantileGradient` loses about $4\varepsilon/\gamma^2$. For $|\gamma| \le 0.03$ with $|\gamma|(1 + |z|) \le 0.2$, where $z$ is the standard Normal quantile, both use instead the Cornish-Fisher series of the standardized gamma quantile through $\gamma^8$, $z + \sum_k \gamma^k P_k(z)$, with $P_1 = (z^2 - 1)/6$, $P_2 = z(z^2 - 7)/144$ and higher polynomials derived exactly, and its derivative. Against 60-digit references the quantile is within about $2 \times 10^{-15}$ in standardized units and the derivative within about $3 \times 10^{-14}$ relative on that region, and both are continuous with the gamma path at its boundary to about $10^{-14}$ and $3 \times 10^{-12}$.

**Maximum-likelihood covariance:** `ParameterCovariance` with maximum likelihood requires $|\gamma| < \sqrt{2}$. The skew information $a^3[\psi'(a) - 1/(a-1) + 1/(2(a-1)^2)]$, with $a = 4/\gamma^2$, is summed from its Bernoulli series through $B_{36}$ when $\gamma^2/4 < 0.15$ and evaluated with the trigamma function above; it agrees with 50-digit references to about $2 \times 10^{-14}$ relative.

**Conditional moments:** `ConditionalMoments(a, b)` returns the conditional mean followed by the second, third and fourth conditional moments about the unconditional mean $\mu$. For $|\gamma| < 0.6$ they are computed in standardized coordinates by the exact boundary recurrence

```math
J_{k+1} = k\,\frac{\gamma}{2}\,J_k + k\,J_{k-1} - \Big[z^k\left(1 + \tfrac{\gamma z}{2}\right) f(z)\Big]_{z_a}^{z_b}
```

where $J_k$ is the integral of $z^k f(z)$ over the standardized interval $[z_a, z_b]$ and $f$ is the standardized Pearson density. The recurrence is continuous through $\gamma = 0$, where it reduces to the truncated Normal recurrence. Larger skews use truncated gamma moments.

### Generalized Extreme Value (GEV)

The Generalized Extreme Value distribution unifies three classical extreme value distributions into a single three-parameter family [[2]](#2). It is the limiting distribution for block maxima (e.g., annual maximum floods, peak wind speeds) under very general conditions described by the Fisher-Tippett-Gnedenko theorem.

#### Mathematical Definition

The GEV is parameterized by location $\xi$, scale $\alpha > 0$, and shape $\kappa$. The ***Numerics*** library uses the Hosking parameterization where the sign convention for $\kappa$ follows L-moment theory.

**Probability Density Function (PDF):**

```math
f(x) = \frac{1}{\alpha} \exp\!\left[-(1-\kappa)y - e^{-y}\right]
```

where the reduced variate $y$ is defined as:

```math
y = \begin{cases}
-\dfrac{1}{\kappa}\ln\!\left(1 - \kappa\dfrac{x - \xi}{\alpha}\right) & \kappa \neq 0 \\[6pt]
\dfrac{x - \xi}{\alpha} & \kappa = 0
\end{cases}
```

**Cumulative Distribution Function (CDF):**

```math
F(x) = \exp(-e^{-y})
```

**Inverse CDF (Quantile Function):**

```math
Q(p) = \begin{cases}
\xi + \dfrac{\alpha}{\kappa}\left[1 - (-\ln p)^{\kappa}\right] & \kappa \neq 0 \\[6pt]
\xi - \alpha \ln(-\ln p) & \kappa = 0
\end{cases}
```

**Support:**

| Shape | Sub-type | Upper/Lower Bound |
|-------|----------|-------------------|
| $\kappa = 0$ | Type I (Gumbel) | $-\infty < x < \infty$ |
| $\kappa > 0$ | Type III (Weibull) | $x \leq \xi + \alpha/\kappa$ (bounded upper tail) |
| $\kappa < 0$ | Type II (Frechet) | $x \geq \xi + \alpha/\kappa$ (bounded lower tail, heavy upper tail) |

**Moments (exist when $|\kappa| < 1$ for the mean, $|\kappa| < 1/2$ for the variance):**

| Property | Formula ($\kappa \neq 0$) | Formula ($\kappa = 0$, Gumbel) |
|----------|--------------------------|-------------------------------|
| Mean | $\xi + \frac{\alpha}{\kappa}[1 - \Gamma(1+\kappa)]$ | $\xi + \alpha\gamma_E$ |
| Variance | $\frac{\alpha^2}{\kappa^2}[\Gamma(1+2\kappa) - \Gamma^2(1+\kappa)]$ | $\frac{\pi^2}{6}\alpha^2$ |
| Mode | $\xi + \frac{\alpha}{\kappa}[(1+\kappa)^{-\kappa} - 1]$ | $\xi$ |

**Parameters in Numerics:** `GeneralizedExtremeValue(location, scale, shape)` where `location` = $\xi$, `scale` = $\alpha > 0$, and `shape` = $\kappa$.

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
if (gev.Kappa > 0)
    Console.WriteLine("  Type III (Weibull) - bounded upper tail");
else if (gev.Kappa < 0)
    Console.WriteLine("  Type II (Fréchet) - heavy upper tail");
else
    Console.WriteLine("  Type I (Gumbel) - exponential tail");
```

`InverseCCDF(q)` returns the annual-maximum level exceeded with probability q. It evaluates the Gumbel coordinate as log(−log1p(−q)) and then the same quantile kernel as `InverseCDF`, so 1 − q is never formed:

```cs
var annualMax = new GeneralizedExtremeValue(10.0, 2.0, -0.2);   // heavy tail
double rare = annualMax.InverseCCDF(1E-20);                     // 100000.0000
double lost = annualMax.InverseCDF(1 - 1E-20);                  // +Infinity: 1 - 1e-20 rounds to 1
```

On a grid of nine shapes in [−0.5, 0.5] and survival probabilities from 0.89 down to 1E-300 (generator `docs/distributions/oracles/gev-survival-quantile.py`), the error against 400-bit references is within 3 + |κ g| ulp with g = log(−log1p(−q)), the amplification of the rounding of g in a heavy tail (65 ulp at most, at κ = −0.5 and q = 1E-300). `InverseCDF(1 - q)` is off by up to 1.3E11 ulp on the same grid for q down to 1E-12.

**L-moment estimation.** `ParametersFromLinearMoments` (and `Estimate` with `MethodOfLinearMoments`) is Hosking's routine PELGEV from the LMOMENTS library, version 3.04 [[12]](#12), followed by one Newton refinement. The shape solves $\tau_3 = 2(1 - 3^{-\kappa})/(1 - 2^{-\kappa}) - 3$. PELGEV gives the rational-function approximation of Donaldson [[13]](#13) for $-0.8 \le \tau_3 < 1$ and Newton-Raphson iteration below $-0.8$; one Newton step on the exact equation, written with `Expm1` so that it stays accurate as $\kappa \to 0$, then refines it. Against the exact root (mpmath, 50 digits) on $\tau_3 \in [-0.999, 0.999]$ the shape is within $3 \times 10^{-14}$ absolute for $\tau_3 \ge -0.8$ and $3 \times 10^{-14}$ relative below, from $3.2 \times 10^{-7}$ and $1.1 \times 10^{-11}$ for PELGEV alone. The scale is $\alpha = \lambda_2 \kappa / (\Gamma(1 + \kappa)(1 - 2^{-\kappa}))$ and the location $\xi = \lambda_1 - \alpha(1 - \Gamma(1 + \kappa))/\kappa$; for $|\kappa| \le 0.05$ both are evaluated from $\ln \Gamma(1 + \kappa)$ and the divided exponential, which pass continuously through the Gumbel values $\alpha = \lambda_2/\ln 2$ and $\xi = \lambda_1 - \gamma\alpha$. The method requires $\lambda_2 > 0$ and $|\tau_3| < 1$.

**L-moments from parameters.** `LinearMomentsFromParameters` evaluates $\lambda_1 = \xi + \alpha(1 - \Gamma(1 + \kappa))/\kappa$, $\lambda_2 = \alpha\Gamma(1 + \kappa)(1 - 2^{-\kappa})/\kappa$, $\tau_3 = 2(1 - 3^{-\kappa})/(1 - 2^{-\kappa}) - 3$ and $\tau_4 = (5(1 - 4^{-\kappa}) - 10(1 - 3^{-\kappa}) + 6(1 - 2^{-\kappa}))/(1 - 2^{-\kappa})$ [[12]](#12) with the same divided differences: $\ln \Gamma(1 + \kappa)$ with the divided exponential, and $(1 - j^{-\kappa})/\kappa = \ln j \cdot \mathrm{expm1}(-\kappa \ln j)/(-\kappa \ln j)$. The result is continuous through κ = 0, where it equals the Gumbel values $\lambda_1 = \xi + \gamma\alpha$, $\lambda_2 = \alpha \ln 2$, $\tau_3 = 2\ln 3/\ln 2 - 3 \approx 0.1699$ and $\tau_4 = 16 - 10\ln 3/\ln 2 \approx 0.1504$. Against mpmath (50 digits) for κ from $\pm 10^{-12}$ to $\pm 0.9$ the relative errors of $\lambda_1$ and $\lambda_2$ and the absolute errors of $\tau_3$ and $\tau_4$ are below $10^{-14}$, and `ParametersFromLinearMoments` recovers the parameters to about $10^{-14}$.

**Growth curves.** The list overload `InverseCDF(IList<double>)` of every distribution evaluates many probabilities under one parameter set, each bitwise equal to the scalar `InverseCDF`. In regional frequency analysis it gives the growth curve of a regional GEV or GPD fitted to L-moments scaled by the index flood:

```cs
var regional = new GeneralizedExtremeValue();
regional.SetParameters(regional.ParametersFromLinearMoments(new[] { 1.0, 0.2, 0.15 }));
double[] growth = regional.InverseCDF(new[] { 0.5, 0.9, 0.98, 0.99, 0.999 });
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

When the threshold is known, `MLE(sample, location)` fits the scale and shape by maximum likelihood with the location held fixed. The observations are the values themselves, not the excesses:

```cs
double[] peaks = { 10450, 11800, 10120, 13900, 10900, 12650, 15400, 10300 };
double[] fitted = new GeneralizedPareto().MLE(peaks, 10000);   // { ξ, α, κ }, with ξ = 10000
var gpdKnown = new GeneralizedPareto(fitted[0], fitted[1], fitted[2]);
Console.WriteLine($"100-yr exceedance level: {gpdKnown.InverseCDF(0.99):F0}");
```

The same fit is available through the `IFixedParameterEstimation` interface, the general form of an estimate at known parameters. `Estimate(sample, method, fixedParameters)` takes one entry per parameter in the order {ξ, α, κ}: a finite value holds that parameter fixed and NaN marks a parameter to estimate. With no parameter fixed it is `Estimate(sample, method)`. A fixed location with `MaximumLikelihood` sets the parameters of `MLE(sample, location)`; any other combination throws `NotSupportedException`.

```cs
var gpdFixed = new GeneralizedPareto();
gpdFixed.Estimate(peaks, ParameterEstimationMethod.MaximumLikelihood, new[] { 10000, double.NaN, double.NaN });
Console.WriteLine($"ξ = {gpdFixed.Xi}, α = {gpdFixed.Alpha:F1}, κ = {gpdFixed.Kappa:F4}");   // ξ = 10000, α = 2585.5, κ = 0.2991
```

**Method of moments:** `DirectMethodOfMoments` (and `Estimate` with `MethodOfMoments`) solves the skewness equation $2(1 - \kappa)\sqrt{1 + 2\kappa}/(1 + 3\kappa) = \gamma$ for κ by the Brent method, then sets $\alpha = \sigma(1 + \kappa)\sqrt{1 + 2\kappa}$ and $\xi = \mu - \alpha/(1 + \kappa)$. These are smooth in κ and pass continuously through the exponential values $\alpha = \sigma$ and $\xi = \mu - \sigma$ at κ = 0. The Brent solve uses an absolute tolerance of $10^{-15}$ with the minimum relative tolerance of `scipy.optimize.brentq`, so against mpmath κ is within $10^{-15}$ absolute and α and ξ within $10^{-15}$ relative for κ from −0.25 to 2, in about 11 skewness evaluations per solve.

The shape κ follows Hosking's sign convention throughout the GEV, GPD, generalized logistic, generalized normal and Kappa-Four classes: a positive κ gives a bounded upper tail. Libraries that use the (μ, σ, c) form of the GPD have c = −κ; `scipy.stats.genextreme` uses c = κ.

#### From Peaks over Threshold to the Annual Maximum

If the number of threshold exceedances in a year is Poisson with mean λ and each exceedance is GPD(ξ, α, κ), the annual maximum is exactly GEV with the same shape [[11]](#11):

```math
\xi_{\text{GEV}} = \xi + \frac{\alpha}{\kappa}\left(1 - \lambda^{-\kappa}\right), \qquad \alpha_{\text{GEV}} = \alpha\,\lambda^{-\kappa}, \qquad \kappa_{\text{GEV}} = \kappa
```

with the Gumbel limit ξ + α ln λ and scale α at κ = 0. `GeneralizedExtremeValue.FromPoissonGeneralizedPareto(rate, location, scale, shape)` returns that GEV:

```cs
// 3.7 exceedances per year of a GPD with threshold 10, scale 2 and a heavy tail (κ = -0.15)
var annual = GeneralizedExtremeValue.FromPoissonGeneralizedPareto(3.7, 10.0, 2.0, -0.15);
Console.WriteLine($"GEV: ξ = {annual.Xi:F4}, α = {annual.Alpha:F4}, κ = {annual.Kappa}");
// GEV: ξ = 12.8911, α = 2.4337, κ = -0.15
double q100 = annual.InverseCDF(0.99);   // 100-year annual-maximum level
```

For λ ≥ 1 the location equals `InverseCCDF(1 / λ)` of the GPD, the level exceeded once per year on average. It is evaluated as ξ + α expm1(−κ ln λ)/(−κ), which stays accurate as κ approaches zero, and the shape is kept exactly: a small nonzero κ is not rounded to zero. The relation holds above the threshold; below it the annual maximum has an atom exp(−λ) that the GEV extends smoothly. The rate and the scale must be positive and finite, and the location and the shape finite.

`GeneralizedPareto.FromPoissonGeneralizedExtremeValue(rate, location, scale, shape)` is the inverse: given the GEV(ξ, α, κ) of the annual maximum and the exceedance rate λ, it returns the GPD of the exceedances,

```math
\xi_{\text{GPD}} = \xi - \frac{\alpha}{\kappa}\left(\lambda^{\kappa} - 1\right), \qquad \alpha_{\text{GPD}} = \alpha\,\lambda^{\kappa}, \qquad \kappa_{\text{GPD}} = \kappa
```

with ξ − α ln λ and scale α at κ = 0. The location is evaluated as ξ − α expm1(κ ln λ)/κ, the shape is kept exactly, and the validation is that of the forward method. The returned threshold is the level whose GEV non-exceedance probability is exp(−λ).

```cs
var gpd = GeneralizedPareto.FromPoissonGeneralizedExtremeValue(3.7, annual.Xi, annual.Alpha, annual.Kappa);
Console.WriteLine($"GPD: ξ = {gpd.Xi:F4}, α = {gpd.Alpha:F4}, κ = {gpd.Kappa}");
// GPD: ξ = 10.0000, α = 2.0000, κ = -0.15
```

#### Survival-Space Quantiles

`GeneralizedPareto.InverseCCDF(q)` returns the level exceeded with probability q directly from q. With L = ln(1/q) it evaluates ξ + α L at κ = 0 and ξ + α expm1(−κ L)/(−κ) otherwise, so small exceedance probabilities keep their full relative precision:

```cs
var tail = new GeneralizedPareto(10.0, 2.0, -0.2);   // heavy tail
double x = tail.InverseCCDF(1E-20);                  // 100000.0000
double y = tail.InverseCDF(1 - 1E-20);               // +Infinity: 1 - 1e-20 rounds to 1
```

`InverseCCDF(0)` returns the upper support limit and `InverseCCDF(1)` the threshold. On 300 random cases (inputs from `numpy.random.default_rng(141)`, generator `docs/distributions/oracles/gpd-survival-quantile.py`) the override is within 2 + |κ| ln(1/q) ulp of a 200-bit reference (8 ulp at most, from the rounding of ln q in heavy tails), while `InverseCDF(1 - q)` is off by up to 6E10 ulp.

## Special Distribution Features

### Truncated Distributions

Create truncated versions of any distribution:

```cs
// Normal truncated to [0, 100]
var truncNormal = new TruncatedNormal(mean: 50, standardDeviation: 15, min: 0, max: 100);

// Or truncate any distribution
var normal = new Normal(50, 15);
var truncated = new TruncatedDistribution(normal, 0, 100);

double mean = truncated.Mean;  // Different from untruncated mean
```

**TruncatedNormal:** When the interval $[a, b]$ contains the parent mean, the mass, CDF, quantile and moments use the closed forms with $Z = \Phi(\beta) - \Phi(\alpha)$. When the interval lies wholly in one tail they are evaluated from $\log\Phi$ of the mirrored lower-tail interval, so intervals such as $[8\sigma, 9\sigma]$ or $[40\sigma, 41\sigma]$ above the mean stay accurate. In the far tails the closed-form skewness and kurtosis lose accuracy to cancellation and are approximate, as in SciPy's `truncnorm`. `Min` = `Max` is a point mass that follows the `Deterministic` conventions, and infinite bounds are rejected. The standardized bounds are cached when a parameter changes, so concurrent evaluations on one instance are safe.

The method of moments sets the bounds to the sample minimum and maximum and solves for the parent $(\mu, \sigma)$ such that the truncated mean and standard deviation equal the sample mean and sample standard deviation (divisor $n - 1$). $\sigma$ is limited to 5 times the sample range. When the sample standard deviation cannot be reached, which is always the case when it exceeds the uniform value $(b - a)/\sqrt{12}$, $\sigma$ is set to that limit and $\mu$ minimizes the squared differences between the fitted and sample moments. A constant sample gives the point mass. `TruncatedNormal` implements `IEstimation` and `IBootstrappable`.

```cs
double[] bounded = { 0.0, 0.45, 0.48, 0.50, 0.50, 0.51, 0.52, 0.55, 1.0 };

var tn = new TruncatedNormal();
tn.Estimate(bounded, ParameterEstimationMethod.MethodOfMoments);
Console.WriteLine($"µ = {tn.Mu:F4}, σ = {tn.Sigma:F4}, [{tn.Min}, {tn.Max}]");
Console.WriteLine($"Mean = {tn.Mean:F4}, SD = {tn.StandardDeviation:F4}"); // match the sample mean and SD (n − 1)

var far = new TruncatedNormal(0, 1, 8, 9); // interval 8σ to 9σ above the mean
Console.WriteLine($"CDF(8.5) = {far.CDF(8.5):F12}");
```

**TruncatedDistribution:** The base-distribution probabilities at the bounds are cached with a snapshot of the base parameters, so a change made to `BaseDistribution` is detected on the next evaluation. When the base CDF at the lower bound exceeds 1/2, the mass, CDF, CCDF and quantile use differences of the base survival function, which keep their relative accuracy in the upper tail. The moments are integrated by adaptive Gauss-Kronrod quadrature over the support (relative tolerance $10^{-12}$), `Mode` is the closed-form base mode, or the nearer bound when the base mode lies outside the bounds, for a unimodal base with a closed-form mode inside its support (Normal, Gamma, Weibull, Gumbel, GEV, Pearson III, Log-Pearson III and the like), and otherwise a global search (the density on 65 points of the central 99.8%, the highest grid peaks polished with Brent's method) that also compares the base distribution's mode, when it lies inside the bounds, and the density at the bounds, and an interval whose truncated mass is below $10^{-15}$ is rejected.

### Mixture Distributions

Model multi-modal data with mixture distributions:

```cs
// Mixture of two normals (bimodal)
var component1 = new Normal(100, 10);
var component2 = new Normal(150, 15);
var weights = new double[] { 0.6, 0.4 };  // 60% from first, 40% from second

var mixture = new Mixture(weights, new UnivariateDistributionBase[] { component1, component2 });

// PDF will show two peaks
double pdf = mixture.PDF(125);  // Valley between modes
```

A mixture density can have several local maxima, so `Mode` is a multi-start global search: Brent's method over the central 99.8% of the mixture and from a bracket about each component with positive weight, keeping the point with the largest density. For a zero-inflated mixture the search is restricted to positive values, and the point mass at zero is not compared with densities. `MinimumOfParameters` and `MaximumOfParameters` are aligned with `GetParameters` (component weights followed by component parameters); `ZeroWeight` is not part of the parameter vector and has no entry.

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

A table of exceedance probabilities is given with a descending probability order. `InverseCCDF(q)` then reads that survival table at q itself, so small exceedance probabilities keep their precision, including those below the 1E-16 far-tail floor of `InverseCDF`:

```cs
double[] levels = { 1, 2, 3, 4, 5, 6, 7, 8 };
double[] exceedance = { 0.5, 0.1, 1E-2, 1E-4, 1E-6, 1E-9, 1E-13, 1E-18 };
var survival = new EmpiricalDistribution(levels, exceedance, SortOrder.Ascending, SortOrder.Descending);

double x12 = survival.InverseCCDF(1E-12);        // between 6 and 7, exact in the table's normal-Z scale
double x17 = survival.InverseCCDF(1E-17);        // between 7 and 8
double lost = survival.InverseCDF(1 - 1E-17);    // 8: 1 - 1e-17 rounds to 1
```

Reading the table back at the returned level gives q to about 1E-14 relative at q = 1E-8, 1E-12, 1E-15 and 1E-17, whereas `InverseCDF(1 - q)` is off by 2E-5 relative at 1E-12 and 8E-4 at 1E-15. For q = 0, for q of at least 0.5, and for ascending tables the result equals `InverseCDF(1 - q)` bit for bit.

`Mode` is the largest point mass when the table has one (repeated x values whose probabilities differ, as a sample with ties has), and otherwise the midpoint of the highest bin of the density grid: the bins are the intervals between consecutive x values, each with density equal to its probability over its width, and ties go to the first. Under the Logarithmic and NormalZ transforms the density varies inside a bin, so the mode is the bin midpoint, not the argmax of `PDF`.

`StrictXValues` and `StrictProbabilities` report whether the table requires strictly monotonic x values and probabilities. A table built from arrays or a sample allows tied x values (a sample with ties) and requires strictly monotonic probabilities. A table built from `OrderedPairedData` takes the strictness of that data, for example tied probabilities for a stepped curve. Both settings are kept by `Clone` and by the XML round trip:

```cs
var levelsTable = new EmpiricalDistribution(new double[] { 1, 2, 3, 4, 5 }, new[] { 0.1, 0.3, 0.5, 0.8, 0.95 });
Console.WriteLine($"{levelsTable.StrictXValues} {levelsTable.StrictProbabilities}");   // False True

var steps = new OrderedPairedData(new double[] { 1, 2, 3, 4 }, new[] { 0.2, 0.5, 0.5, 0.9 },
    true, SortOrder.Ascending, false, SortOrder.Ascending);
var stepped = new EmpiricalDistribution(steps);
Console.WriteLine($"{stepped.StrictXValues} {stepped.StrictProbabilities}");           // True False

var restored = (EmpiricalDistribution)UnivariateDistributionFactory.CreateDistribution(stepped.ToXElement());
Console.WriteLine($"{restored.StrictXValues} {restored.StrictProbabilities}");         // True False
```

### Kernel Density Estimation

Smooth non-parametric density estimation:

```cs
var kde = new KernelDensity(observations, KernelDensity.KernelType.Gaussian, 1.5);

// Smooth PDF
double density = kde.PDF(15.0);

// KDE-based CDF and quantiles
double cdf = kde.CDF(15.0);
double quantile = kde.InverseCDF(0.75);
```

`Mode` is the global maximum of the estimate over its central 99.8%: the density is ranked on a grid with spacing at most a quarter of the bandwidth, using the linearly binned estimate (Wand and Jones, 1995), and the highest peaks are polished with Brent's method on the exact density. A bimodal estimate returns its taller peak.

### PERT Distributions

For expert judgment and project management:

```cs
// PERT from minimum, most likely, maximum
var pert = new Pert(min: 10, mode: 15, max: 25);

// PERT from percentile judgments
var pertPercentile = new PertPercentile(fifth: 12, fiftieth: 15, ninetyFifth: 22);

// Use for duration or cost uncertainty
double expectedDuration = pert.Mean;
double variance = pert.Variance;
```

**Pert:** `Median` is the exact median of the underlying four-parameter Beta distribution, `InverseCDF(0.5)`, rather than the approximation $(a + 6c + b)/8$. `Pert(a, a, a)` is the fixed value $a$. Maximum likelihood uses the scale-aware bounds described for the Triangular distribution, so samples whose minimum or maximum is 0 can be fitted, and `ParameterEstimationMethod.MethodOfPercentiles` with equal 5th and 95th sample percentiles gives the fixed-value `Pert(a, a, a)`.

```cs
var pert = new Pert(min: 10, mode: 15, max: 25);
Console.WriteLine($"Median = {pert.Median:F6}"); // exact Beta median, InverseCDF(0.5)

double[] fixedSample = { 7, 7, 7, 7, 7 };
var pertFixed = new Pert();
pertFixed.Estimate(fixedSample, ParameterEstimationMethod.MethodOfPercentiles); // Pert(7, 7, 7)
Console.WriteLine($"{pertFixed.Min}, {pertFixed.MostLikely}, {pertFixed.Max}");
```

**PertPercentile:** The PERT parameters are solved from the 5th, 50th and 95th percentiles on first use of a summary property, a probability function or `ToPert`; equal percentiles give the fixed-value PERT. `MinAllowableValue` and `MaxAllowableValue` censor the PERT: values beyond a bound are set to the bound, so each bound carries a point mass equal to the PERT probability beyond it. The CDF is 0 below the minimum and 1 at and above the maximum, the PDF is zero outside the bounds, the moments are those of the censored variable, the median is clamped to the bounds, and quantiles inside the bounds are unchanged. Finite bounds are written to and read from the XML representation.

```cs
// Expert judgment: fifth = 0, median = 0, ninety-fifth = 1, with values below zero impossible
var lifeLoss = new PertPercentile(fifth: 0, fiftieth: 0, ninetyFifth: 1);
lifeLoss.MinAllowableValue = 0;

Console.WriteLine($"CDF(0) = {lifeLoss.CDF(0):F4}");             // point mass at the bound
Console.WriteLine($"PERT CDF(0) = {lifeLoss.ToPert().CDF(0):F4}"); // the uncensored PERT probability below 0
Console.WriteLine($"CDF(-0.1) = {lifeLoss.CDF(-0.1)}");          // 0 below the bound
Console.WriteLine($"Mean = {lifeLoss.Mean:F4}");                  // mean of the censored variable
Console.WriteLine($"Q(0.95) = {lifeLoss.InverseCDF(0.95):F4}");   // quantiles above the bound are those of the PERT
```

**PertPercentileZ:** This variant describes a probability $X = \Phi(Z)$, where $Z$ follows a PERT distribution with minimum $a$ and maximum $b$ fitted to the z-scores of the three percentiles, which must lie strictly in (0, 1). The support is $[\Phi(a), \Phi(b)]$ inside (0, 1). By the change of variables, with $g$ and $G$ the PERT density and CDF in z-space,

```math
f(x) = \frac{g\left(\Phi^{-1}(x)\right)}{\phi\left(\Phi^{-1}(x)\right)}, \qquad F(x) = G\left(\Phi^{-1}(x)\right), \qquad Q(p) = \Phi\left(G^{-1}(p)\right)
```

The PDF is 0 and the CDF is 0 or 1 outside the support. The moments of $X$ are integrated by quadrature in z-space, the median is $\Phi$ of the median of $Z$, and the mode is the global maximum of $f$, which can be bimodal or peak at a bound when the z-range is wide.

```cs
// Percentiles of a probability, fitted as a PERT in standard Normal z-space
var pz = new PertPercentileZ(fifth: 0.02, fiftieth: 0.10, ninetyFifth: 0.40);

Console.WriteLine($"Support = [{pz.Minimum:F6}, {pz.Maximum:F6}]");
Console.WriteLine($"Median  = {pz.Median:F6}");       // Φ(median of Z) = 0.10
Console.WriteLine($"Mean    = {pz.Mean:F6}");         // E[Φ(Z)], not Φ(E[Z])
Console.WriteLine($"CDF(0.4) = {pz.CDF(0.40):F6}");   // 0.95
```

### Competing Risks

Model the distribution of the minimum (or maximum) of multiple independent or correlated random variables. This is commonly used in system reliability analysis where failure occurs when any component fails:

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

// Define failure modes for a levee system
var overtopping = new Normal(18.5, 2.0);   // Overtopping stage (ft)
var seepage = new LogNormal(2.85, 0.15);   // Seepage failure stage (ft)
var erosion = new GeneralizedExtremeValue(16.0, 2.5, -0.1);   // Erosion failure stage (ft)

// System fails at the MINIMUM failure stage
var system = new CompetingRisks(new UnivariateDistributionBase[] {
    overtopping, seepage, erosion
});

// Default: MinimumOfRandomVariables = true (first failure)
Console.WriteLine("System Failure Analysis (Independent Components):");
Console.WriteLine($"  Mean failure stage: {system.Mean:F2} ft");
Console.WriteLine($"  Median failure stage: {system.InverseCDF(0.5):F2} ft");
Console.WriteLine($"  P(failure ≤ 15 ft): {system.CDF(15.0):F4}");
Console.WriteLine($"  1% failure stage: {system.InverseCDF(0.01):F2} ft");
```

**Dependency options:**

```cs
// Independent components (default)
system.Dependency = Probability.DependencyType.Independent;

// Perfectly correlated: system CDF equals the weakest component
system.Dependency = Probability.DependencyType.PerfectlyPositive;

// Custom correlation structure
system.Dependency = Probability.DependencyType.CorrelationMatrix;
system.CorrelationMatrix = new double[,] {
    { 1.0, 0.6, 0.3 },
    { 0.6, 1.0, 0.4 },
    { 0.3, 0.4, 1.0 }
};

// Switch to maximum of random variables
system.MinimumOfRandomVariables = false;
```

Correlated dependency is evaluated with the product-of-conditional-marginals (PCM) approximation of the joint probabilities. `ToXElement` preserves component settings that are not parameters, such as the logarithm base of a `LogPearsonTypeIII` or `LogNormal` component, and nests table-backed components such as `EmpiricalDistribution`, so `UnivariateDistributionFactory.CreateDistribution` restores the same system.

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
var lp3 = new LogPearsonTypeIII(meanOfLog: 10.2, standardDeviationOfLog: 0.3, skewOfLog: 0.4);

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
var weibull = new Weibull(scale: 1000, shape: 2.5); // hours

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
var failureProb = new BetaDistribution(alpha: 2, beta: 1998); // ~0.001

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
| Circular/directional data | Von Mises |

## Parameter Bounds and Validation

All distributions validate parameters:

```cs
var gev = new GeneralizedExtremeValue();

// Check if parameters are valid
var parameters = new double[] { 1000, 200, 0.3 };
var exception = gev.ValidateParameters(parameters, throwException: false);

if (exception == null)
{
    gev.SetParameters(parameters);
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

## Parameter Interpretation Guide

Most continuous distributions can be understood through three fundamental types of parameters, each controlling a distinct aspect of the distribution's behavior:

### Location Parameters ($\xi$, $\mu$)

Location parameters shift the entire distribution left or right along the real line without changing its shape or spread. Changing the location parameter translates every quantile by the same amount.

- **Normal:** $\mu$ (mean) shifts the center of symmetry
- **Gumbel / GEV:** $\xi$ shifts the mode
- **Exponential:** $\xi$ shifts the lower bound of support

### Scale Parameters ($\alpha$, $\sigma$, $\theta$)

Scale parameters stretch or compress the distribution. Multiplying the scale by a constant $c$ multiplies the standard deviation (and all quantile deviations from the location) by $c$, without changing the shape.

- **Normal:** $\sigma$ (standard deviation) controls the spread
- **Gumbel / GEV / Exponential:** $\alpha$ controls the spread
- **Gamma:** $\theta$ scales the distribution; mean = $\kappa\theta$

### Shape Parameters ($\kappa$, $\alpha$, $\beta$, $\gamma$)

Shape parameters control the fundamental character of the distribution: tail weight, skewness, peakedness, and boundedness. Unlike location and scale, changing a shape parameter alters the qualitative behavior of the distribution.

- **GEV:** $\kappa$ determines the tail type (bounded vs. heavy-tailed vs. exponential)
- **Gamma:** $\kappa$ controls skewness ($2/\sqrt{\kappa}$); as $\kappa \to \infty$, the Gamma approaches the Normal
- **Beta:** $\alpha$ and $\beta$ jointly control the shape on $[0,1]$; equal values produce symmetry
- **LP3:** $\gamma$ (skewness of log) controls asymmetry; $\gamma = 0$ reduces LP3 to Log-Normal

## Distribution Relationships

Many distributions in the ***Numerics*** library are connected through limiting cases, transformations, or special parameterizations. Understanding these relationships helps in selecting appropriate models and verifying analytical results.

### Extreme Value Family

```
GEV(ξ, α, κ)
├── κ = 0  →  Gumbel(ξ, α)          [Type I: exponential tails]
├── κ > 0  →  Weibull-type           [Type III: bounded upper tail]
└── κ < 0  →  Fréchet-type           [Type II: heavy upper tail]
```

The Gumbel distribution `Gumbel(ξ, α)` is exactly equivalent to `GeneralizedExtremeValue(ξ, α, 0)`.

### Gamma Family

```
GammaDistribution(θ, κ)
├── κ = 1        →  Exponential(0, θ)     [memoryless special case]
├── κ = ν/2, θ=2 →  Chi-Squared(ν)        [sum of squared standard Normals]
└── κ → ∞       →  Normal(κθ, θ√κ)       [by Central Limit Theorem]
```

### Logarithmic Transforms

```
X ~ LogNormal(μ, σ)        ⟺  log(X) ~ Normal(μ, σ)
X ~ LogPearsonTypeIII(μ, σ, γ)  ⟺  log(X) ~ PearsonTypeIII(μ, σ, γ)
```

When the LP3 skewness parameter $\gamma = 0$, the LP3 reduces to the Log-Normal distribution.

### Beta and Uniform

```
Beta(1, 1) = Uniform(0, 1)
```

The Beta distribution with both shape parameters equal to 1 produces a uniform density on $[0, 1]$.

### Student's t and Normal

```
Student-t(ν) → Normal(0, 1)  as ν → ∞
```

The Student's t distribution approaches the standard Normal distribution as the degrees of freedom increase. For $\nu \ge 10^8$ the `StudentT` PDF uses the Normal distribution directly; the CDF, its complement and the quantile use the Student's t tail for every $\nu$.

---

## References

<a id="1">[1]</a> Interagency Advisory Committee on Water Data. (2019). *Guidelines for Determining Flood Flow Frequency, Bulletin 17C*. U.S. Geological Survey Techniques and Methods, Book 4, Chapter B5.

<a id="2">[2]</a> Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.

<a id="3">[3]</a> Hosking, J. R. M., & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press.

<a id="4">[4]</a> Johnson, N. L., Kotz, S., & Balakrishnan, N. (1994-1995). *Continuous Univariate Distributions*, Vols. 1-2 (2nd ed.). Wiley.

<a id="5">[5]</a> Kaas, R., & Buhrman, J. M. (1980). Mean, median and mode in binomial distributions. *Statistica Neerlandica*, 34(1), 13-18. https://doi.org/10.1111/j.1467-9574.1980.tb00681.x

<a id="6">[6]</a> Guenther, W. C. (1978). Evaluation of probabilities for the noncentral distributions and the difference of two T-variables with a desk calculator. *Journal of Statistical Computation and Simulation*, 6(3-4), 199-206. https://doi.org/10.1080/00949657808810188

<a id="7">[7]</a> Lenth, R. V. (1989). Algorithm AS 243: Cumulative distribution function of the non-central t distribution. *Applied Statistics*, 38(1), 185-189. https://doi.org/10.2307/2347693

<a id="8">[8]</a> Benton, D., & Krishnamoorthy, K. (2003). Computing discrete mixtures of continuous distributions: Noncentral chisquare, noncentral t and the distribution of the square of the sample multiple correlation coefficient. *Computational Statistics & Data Analysis*, 43(2), 249-267. https://doi.org/10.1016/S0167-9473(02)00283-9

<a id="9">[9]</a> Mardia, K. V., & Jupp, P. E. (2000). *Directional Statistics*. Wiley.

<a id="10">[10]</a> Olver, F. W. J., Olde Daalhuis, A. B., Lozier, D. W., Schneider, B. I., Boisvert, R. F., Clark, C. W., Miller, B. R., Saunders, B. V., Cohl, H. S., & McClain, M. A. (Eds.). *NIST Digital Library of Mathematical Functions*, Section 8.12. https://dlmf.nist.gov/8.12

<a id="11">[11]</a> Madsen, H., Rasmussen, P. F., & Rosbjerg, D. (1997). Comparison of annual maximum series and partial duration series methods for modeling extreme hydrologic events: 1. At-site modeling. *Water Resources Research*, 33(4), 747-757. https://doi.org/10.1029/96WR03848

<a id="12">[12]</a> Hosking, J. R. M. (1996). *Fortran routines for use with the method of L-moments, Version 3*. Research Report RC20525, IBM Research Division, Yorktown Heights, N.Y.

<a id="13">[13]</a> Donaldson, R. W. (1996). Calculating inverse CV, skew and PWM functions for Pearson-3, log-normal, extreme-value, and log-logistic distributions. *Communications in Statistics - Simulation and Computation*, 25(3), 741-747. https://doi.org/10.1080/03610919608813340

---

[← Previous: Hypothesis Tests](../statistics/hypothesis-tests.md) | [Back to Index](../index.md) | [Next: Parameter Estimation →](parameter-estimation.md)
