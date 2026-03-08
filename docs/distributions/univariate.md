# Univariate Distributions

[← Previous: Hypothesis Tests](../statistics/hypothesis-tests.md) | [Back to Index](../index.md) | [Next: Parameter Estimation →](parameter-estimation.md)

The ***Numerics*** library provides over 40 univariate probability distributions for statistical analysis, risk assessment, and uncertainty quantification. All distributions implement a common interface with consistent methods for computing probability density functions (PDF), cumulative distribution functions (CDF), quantiles, and statistical moments.

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
| **Student's t** | ν (degrees of freedom) | Heavy-tailed alternative to Normal |
| **Noncentral t** | ν, δ (noncentrality) | Power analysis, hypothesis testing |
| **Chi-Squared** | k (degrees of freedom) | Variance estimation, goodness-of-fit |
| **Inverse Gamma** | α, β | Bayesian priors for variance |
| **Inverse Chi-Squared** | ν | Bayesian inference |
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

**Moments:**

| Property | Formula |
|----------|---------|
| Mean | $E[X] = \kappa\theta$ |
| Variance | $\text{Var}(X) = \kappa\theta^2$ |
| Skewness | $\gamma_1 = 2/\sqrt{\kappa}$ |
| Kurtosis | $\kappa_4 = 3 + 6/\kappa$ |
| Mode | $(\kappa - 1)\theta$ for $\kappa \geq 1$ |

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
| Skewness | $\gamma_1 \approx 1.1396$ |
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
| Mode | $\frac{\alpha - 1}{\alpha + \beta - 2}$ for $\alpha, \beta > 1$ |
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

**Quantile computation:** In practice, LP3 quantiles are computed using the frequency factor approach from Bulletin 17C:

```math
\log_b(X_p) = \mu + K_p \cdot \sigma
```

where $K_p$ is the Pearson Type III frequency factor (a function of the skewness $\gamma$ and probability $p$), computed via the Cornish-Fisher or Wilson-Hilferty approximation.

**Moments (of the log-transformed data):**

| Property | Formula |
|----------|---------|
| Mean of log | $\mu$ |
| Std dev of log | $\sigma$ |
| Skewness of log | $\gamma$ |

When $\gamma = 0$, the LP3 reduces to the Log-Normal distribution.

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
var truncNormal = new TruncatedNormal(mean: 50, standardDeviation: 15, min: 0, max: 100);

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

var mixture = new Mixture(weights, new UnivariateDistributionBase[] { component1, component2 });

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
var kde = new KernelDensity(observations, KernelDensity.KernelType.Gaussian, 1.5);

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
var pertPercentile = new PertPercentile(fifth: 12, fiftieth: 15, ninetyFifth: 22);

// Use for duration or cost uncertainty
double expectedDuration = pert.Mean;
double variance = pert.Variance;
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

// Perfectly correlated — system CDF equals the weakest component
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

The Student's t distribution approaches the standard Normal distribution as the degrees of freedom increase.

---

## References

<a id="1">[1]</a> Interagency Advisory Committee on Water Data. (2019). *Guidelines for Determining Flood Flow Frequency, Bulletin 17C*. U.S. Geological Survey Techniques and Methods, Book 4, Chapter B5.

<a id="2">[2]</a> Coles, S. (2001). *An Introduction to Statistical Modeling of Extreme Values*. Springer.

<a id="3">[3]</a> Hosking, J. R. M., & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press.

<a id="4">[4]</a> Johnson, N. L., Kotz, S., & Balakrishnan, N. (1994-1995). *Continuous Univariate Distributions*, Vols. 1-2 (2nd ed.). Wiley.

---

[← Previous: Hypothesis Tests](../statistics/hypothesis-tests.md) | [Back to Index](../index.md) | [Next: Parameter Estimation →](parameter-estimation.md)
