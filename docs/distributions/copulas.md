# Copulas

[← Previous: Uncertainty Analysis](uncertainty-analysis.md) | [Back to Index](../index.md) | [Next: Multivariate Distributions →](multivariate.md)

Copulas separate the dependence structure of multivariate distributions from their marginal distributions. This separation is formalized by **Sklar's theorem** [[1]](#1): for any multivariate distribution $F$ with marginal CDFs $F_1, F_2, \ldots, F_n$, there exists a copula $C$ such that:

```math
F(x_1, x_2, \ldots, x_n) = C(F_1(x_1), F_2(x_2), \ldots, F_n(x_n))
```

A copula $C: [0,1]^n \rightarrow [0,1]$ is itself a multivariate CDF with uniform marginals. The power of this decomposition is that it allows us to:
1. Model marginal distributions independently (each can be any distribution)
2. Model dependence separately via the copula
3. Combine them to form any joint distribution

The ***Numerics*** library provides bivariate copula functions for modeling dependence between two random variables in risk assessment and multivariate analysis.

## Available Copulas

### Elliptical Copulas

Elliptical copulas are derived from elliptical distributions (Normal, Student's t). They produce symmetric dependence structures.

#### Normal (Gaussian) Copula

The Normal copula is derived from the bivariate normal distribution. Its density function is:

```math
c(u, v) = \frac{1}{\sqrt{1-\rho^2}} \exp\left(-\frac{\rho^2 s^2 + \rho^2 t^2 - 2\rho s t}{2(1-\rho^2)}\right)
```

where $s = \Phi^{-1}(u)$ and $t = \Phi^{-1}(v)$ are the standard normal quantiles, and $\rho \in [-1, 1]$ is the correlation parameter.

The Normal copula has **no tail dependence** ($\lambda_L = \lambda_U = 0$), meaning that extreme events in one variable do not increase the probability of extreme events in the other beyond what the overall correlation implies. This makes it unsuitable for modeling joint extremes.

The relationship between Kendall's tau and the copula parameter is:

```math
\tau = \frac{2}{\pi} \arcsin(\rho)
```

```cs
using Numerics.Distributions.Copulas;

// Normal (Gaussian) copula with correlation rho
double rho = 0.7;  // Correlation coefficient

var normalCopula = new NormalCopula(rho);

// Evaluate copula density
double u1 = 0.3, u2 = 0.7;
double density = normalCopula.PDF(u1, u2);

Console.WriteLine($"Normal copula density: {density:F4}");
```

#### Student's t Copula

The Student's t copula extends the Normal copula by adding **symmetric tail dependence**. Its density involves the bivariate Student's t distribution with $\nu$ degrees of freedom and correlation $\rho$. The density is computed in log-space for numerical stability:

```math
\log c(u,v) = \log \Gamma\!\left(\frac{\nu+2}{2}\right) + \log \Gamma\!\left(\frac{\nu}{2}\right) - 2\log \Gamma\!\left(\frac{\nu+1}{2}\right) - \frac{1}{2}\log(1-\rho^2) - \frac{\nu+2}{2}\log\!\left(1 + \frac{Q}{\nu(1-\rho^2)}\right) + \frac{\nu+1}{2}\left[\log\!\left(1 + \frac{x_1^2}{\nu}\right) + \log\!\left(1 + \frac{x_2^2}{\nu}\right)\right]
```

where $x_1 = t_\nu^{-1}(u)$, $x_2 = t_\nu^{-1}(v)$, and $Q = x_1^2 - 2\rho x_1 x_2 + x_2^2$.

The **tail dependence coefficient** is:

```math
\lambda = \lambda_L = \lambda_U = 2 \cdot t_{\nu+1}\!\left(-\sqrt{\frac{(\nu+1)(1-\rho)}{1+\rho}}\right)
```

As $\nu \to \infty$, the t copula converges to the Normal copula and tail dependence vanishes. As $\nu$ decreases, tail dependence increases.

```cs
// t-copula with 5 degrees of freedom
int nu = 5;
var tCopula = new StudentTCopula(rho, nu);

double density = tCopula.PDF(u1, u2);

Console.WriteLine($"t-copula density: {density:F4}");
Console.WriteLine($"Upper tail dependence: {tCopula.UpperTailDependence:F4}");
Console.WriteLine($"Lower tail dependence: {tCopula.LowerTailDependence:F4}");
```

### Archimedean Copulas

Archimedean copulas are defined through a **generator function** $\varphi(t)$, a continuous, strictly decreasing, convex function with $\varphi(1) = 0$. The copula CDF is expressed as:

```math
C(u, v) = \varphi^{-1}(\varphi(u) + \varphi(v))
```

The copula density is obtained by differentiating:

```math
c(u,v) = -\frac{\varphi''(\varphi^{-1}(\varphi(u) + \varphi(v)))}{[\varphi'(\varphi^{-1}(\varphi(u) + \varphi(v)))]^3} \cdot \varphi'(u) \cdot \varphi'(v)
```

Each Archimedean family is characterized by its generator, which determines the copula's dependence properties. The ***Numerics*** library implements five Archimedean families:

#### Clayton Copula

The Clayton copula exhibits **lower tail dependence**, making it suitable for modeling the joint occurrence of low values (e.g., concurrent droughts in multiple watersheds).

| Property | Value |
|----------|-------|
| Generator | $\varphi(t) = t^{-\theta} - 1$ |
| Inverse generator | $\varphi^{-1}(s) = (1+s)^{-1/\theta}$ |
| CDF | $C(u,v) = (u^{-\theta} + v^{-\theta} - 1)^{-1/\theta}$ |
| Parameter range | $\theta \in (0, \infty)$ |
| Kendall's tau | $\tau = \frac{\theta}{\theta + 2}$ |
| Lower tail dependence | $\lambda_L = 2^{-1/\theta}$ |
| Upper tail dependence | $\lambda_U = 0$ |

```cs
// Clayton copula (lower tail dependence), θ ∈ (0, ∞)
var claytonCopula = new ClaytonCopula(2.0);
```

#### Gumbel Copula

The Gumbel copula exhibits **upper tail dependence**, making it the natural choice for modeling joint flood events (e.g., concurrent high flows in tributaries).

| Property | Value |
|----------|-------|
| Generator | $\varphi(t) = (-\ln t)^{\theta}$ |
| Inverse generator | $\varphi^{-1}(s) = \exp(-s^{1/\theta})$ |
| Parameter range | $\theta \in [1, \infty)$ |
| Kendall's tau | $\tau = 1 - \frac{1}{\theta}$ |
| Lower tail dependence | $\lambda_L = 0$ |
| Upper tail dependence | $\lambda_U = 2 - 2^{1/\theta}$ |

```cs
// Gumbel copula (upper tail dependence), θ ∈ [1, ∞)
var gumbelCopula = new GumbelCopula(2.0);
```

#### Frank Copula

The Frank copula has **no tail dependence** and produces a symmetric dependence structure. It is useful for modeling moderate, general correlation without emphasizing joint extremes.

| Property | Value |
|----------|-------|
| Generator | $\varphi(t) = -\ln\!\left(\frac{e^{-\theta t} - 1}{e^{-\theta} - 1}\right)$ |
| CDF | $C(u,v) = -\frac{1}{\theta}\ln\!\left(1 + \frac{(e^{-\theta u}-1)(e^{-\theta v}-1)}{e^{-\theta}-1}\right)$ |
| Parameter range | $\theta \in (-\infty, \infty) \setminus \lbrace 0\rbrace$ |
| Tail dependence | $\lambda_L = \lambda_U = 0$ |

The Frank copula is the only Archimedean copula that allows both positive and negative dependence ($\theta > 0$ for positive, $\theta < 0$ for negative).

```cs
// Frank copula (no tail dependence), θ ∈ (-∞, ∞) \ {0}
var frankCopula = new FrankCopula(5.0);
```

#### Joe Copula

The Joe copula has **upper tail dependence**, similar to the Gumbel copula, but with a different dependence structure in the body of the distribution.

| Property | Value |
|----------|-------|
| Generator | $\varphi(t) = -\ln(1 - (1-t)^{\theta})$ |
| Parameter range | $\theta \in [1, \infty)$ |
| Lower tail dependence | $\lambda_L = 0$ |
| Upper tail dependence | $\lambda_U = 2 - 2^{1/\theta}$ |

```cs
// Joe copula (upper tail dependence), θ ∈ [1, ∞)
var joeCopula = new JoeCopula(2.0);
```

#### Ali-Mikhail-Haq (AMH) Copula

The AMH copula models **weak dependence structures** and has no tail dependence. Its parameter range is limited to $[-1, 1]$, which restricts it to Kendall's tau values in approximately $[-0.182, 0.333]$.

| Property | Value |
|----------|-------|
| Generator | $\varphi(t) = \ln\!\left(\frac{1-\theta(1-t)}{t}\right)$ |
| Parameter range | $\theta \in [-1, 1]$ |
| Kendall's tau range | $\tau \in \left[\frac{5 - 8\ln 2}{3}, \frac{1}{3}\right] \approx [-0.182, 0.333]$ |
| Tail dependence | $\lambda_L = \lambda_U = 0$ |

```cs
// Ali-Mikhail-Haq copula (weak dependence), θ ∈ [-1, 1]
var amhCopula = new AMHCopula(0.5);
```

### Copula Selection Guide

| Copula | Tail Dependence | Parameter Range | Best For |
|--------|----------------|-----------------|----------|
| Normal | None | $\rho \in [-1, 1]$ | General symmetric dependence |
| Student-t | Symmetric | $\rho \in [-1, 1]$, $\nu > 2$ | Heavy-tailed joint extremes |
| Clayton | Lower tail | $\theta \in (0, \infty)$ | Joint low extremes (droughts) |
| Gumbel | Upper tail | $\theta \in [1, \infty)$ | Joint high extremes (floods) |
| Frank | None | $\theta \in \mathbb{R} \setminus \lbrace 0\rbrace$ | Moderate symmetric dependence |
| Joe | Upper tail | $\theta \in [1, \infty)$ | Strong upper tail dependence |
| AMH | None | $\theta \in [-1, 1]$ | Weak dependence structures |

## Fitting Copulas to Data

The copula parameter can be estimated from data using the relationship between Kendall's tau and the copula parameter. The ***Numerics*** library provides the `SetThetaFromTau` method for Archimedean copulas:

```cs
using System.Linq;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

// Sample paired observations (e.g., peak flow and volume)
double[] x = { 1200, 1500, 1100, 1800, 1350, 1600, 1250, 1450, 1900, 1300 };
double[] y = { 45, 52, 42, 65, 48, 58, 44, 51, 68, 46 };

// Step 1: Fit marginal distributions
var gevX = new GeneralizedExtremeValue();
gevX.Estimate(x, ParameterEstimationMethod.MethodOfLinearMoments);

var gevY = new GeneralizedExtremeValue();
gevY.Estimate(y, ParameterEstimationMethod.MethodOfLinearMoments);

// Step 2: Transform to uniform margins using fitted CDFs
double[] u = x.Select(xi => gevX.CDF(xi)).ToArray();
double[] v = y.Select(yi => gevY.CDF(yi)).ToArray();

// Step 3: Estimate copula parameter using rank correlation
double tau = Correlation.KendallsTau(x, y);
Console.WriteLine($"Kendall's tau: {tau:F3}");

// The tau-to-parameter relationships:
// Clayton: θ = 2τ/(1-τ) for τ > 0
// Gumbel:  θ = 1/(1-τ) for τ > 0
// Normal:  ρ = sin(πτ/2)
```

Alternatively, copula parameters can be estimated by maximizing the pseudo-log-likelihood, which uses only the copula density evaluated at empirical probability-integral transforms:

```cs
// The BivariateCopula base class provides likelihood methods:
// PseudoLogLikelihood  — copula-only log-likelihood
// IFMLogLikelihood     — inference functions for margins (with pre-estimated marginals)
// LogLikelihood        — full log-likelihood (copula + marginals)
```

## Practical Example: Bivariate Distribution

Construct a bivariate distribution with arbitrary marginals and specified dependence:

```cs
using System.Linq;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

// Step 1: Define marginal distributions
var margin1 = new LogNormal(4.0, 0.5);  // Streamflow
var margin2 = new Gumbel(100, 20);      // Peak stage

// Step 2: Define dependence via copula
double rho = 0.8;  // Strong positive correlation
var copula = new NormalCopula(rho);

// Step 3: Sample from joint distribution
int n = 1000;
var samples = copula.GenerateRandomValues(n);

// Transform uniforms to actual distributions
double[] flow = new double[n];
double[] stage = new double[n];

for (int i = 0; i < n; i++)
{
    flow[i] = margin1.InverseCDF(samples[i, 0]);
    stage[i] = margin2.InverseCDF(samples[i, 1]);
}

Console.WriteLine($"Generated {n} correlated samples");
Console.WriteLine($"Flow range: [{flow.Min():F1}, {flow.Max():F1}]");
Console.WriteLine($"Stage range: [{stage.Min():F1}, {stage.Max():F1}]");

// Empirical correlation
var correlation = Correlation.Pearson(flow, stage);
Console.WriteLine($"Sample correlation: {correlation:F3}");
```

The sampling algorithm uses the **conditional method**: generate $u$ from $U(0,1)$, then generate $v$ from the conditional distribution $C(v|u) = \frac{\partial C(u,v)}{\partial u}$ by inverting this conditional CDF. The `GenerateRandomValues` method uses Latin Hypercube sampling for better coverage of the probability space.

## Applications in Risk Assessment

### Joint Probability Analysis

For dam safety, estimate probability of joint high flow and high stage:

```cs
// Critical thresholds
double flowThreshold = 10000;   // cfs
double stageThreshold = 150;    // feet

// Compute joint exceedance probability
double u1 = margin1.CDF(flowThreshold);
double u2 = margin2.CDF(stageThreshold);

// P(Flow > threshold AND Stage > threshold)
// = 1 - u1 - u2 + C(u1, u2)
double jointExceedance = copula.ANDJointExceedanceProbability(u1, u2);

Console.WriteLine($"Joint exceedance probability: {jointExceedance:E4}");
Console.WriteLine($"Return period: {1.0 / jointExceedance:F1} years");
```

The AND joint exceedance probability is computed as $P(X > x \text{ and } Y > y) = 1 - u - v + C(u, v)$, while the OR joint exceedance is $P(X > x \text{ or } Y > y) = 1 - C(u, v)$.

### Conditional Distributions

Given flow, what is the conditional distribution of stage? The conditional CDF can be computed numerically using the copula CDF via partial differentiation: $C(v|u) = \frac{\partial C(u,v)}{\partial u}$.

```cs
// Observed flow
double observedFlow = 12000;
double uFlow = margin1.CDF(observedFlow);

// Approximate the conditional CDF: dC(u,v)/du via finite difference
Func<double, double> conditionalCDF = (stage) =>
{
    double uStage = margin2.CDF(stage);
    double du = 1e-6;
    double uPlus = Math.Min(uFlow + du, 1.0);
    double uMinus = Math.Max(uFlow - du, 0.0);
    return (copula.CDF(uPlus, uStage) - copula.CDF(uMinus, uStage)) / (uPlus - uMinus);
};

// Conditional probability at specific values
Console.WriteLine($"Given flow = {observedFlow:F0} cfs:");
Console.WriteLine($"  P(Stage > 15 | Flow = {observedFlow}) = {1 - conditionalCDF(15):P1}");
```

## Tail Dependence

Tail dependence measures the probability of observing an extreme value in one variable given that the other is already extreme. The **upper tail dependence coefficient** is:

```math
\lambda_U = \lim_{u \to 1^-} P(Y > F_Y^{-1}(u) \mid X > F_X^{-1}(u)) = \lim_{u \to 1^-} \frac{1 - 2u + C(u,u)}{1-u}
```

The **lower tail dependence coefficient** is defined analogously as $u \to 0^+$:

```math
\lambda_L = \lim_{u \to 0^+} \frac{C(u,u)}{u}
```

A copula has tail dependence if $\lambda > 0$. This is a critical property for risk analysis: the Normal copula always has $\lambda = 0$, meaning it systematically underestimates joint extreme events. For flood analysis, copulas with upper tail dependence (Gumbel, Joe, Student-t) are preferred.

```cs
// Tail dependence comparison using UpperTailDependence / LowerTailDependence properties
Console.WriteLine("Tail Dependence Properties:");
Console.WriteLine($"  Student-t(ρ=0.7, ν=5):  λ_U = {new StudentTCopula(0.7, 5).UpperTailDependence:F4}");
Console.WriteLine($"  Student-t(ρ=0.7, ν=20): λ_U = {new StudentTCopula(0.7, 20).UpperTailDependence:F4}");
Console.WriteLine($"  Normal(ρ=0.7):           λ_U = {new NormalCopula(0.7).UpperTailDependence:F4}");
Console.WriteLine($"  Clayton(θ=2):            λ_L = {new ClaytonCopula(2.0).LowerTailDependence:F4}");
Console.WriteLine($"  Gumbel(θ=2):             λ_U = {new GumbelCopula(2.0).UpperTailDependence:F4}");
Console.WriteLine($"  Frank(θ=5):              λ_U = {new FrankCopula(5.0).UpperTailDependence:F4}");
```

## Higher-Dimensional Dependence

For problems with more than two variables, there are several approaches:

1. **Multivariate Normal distribution**: Use when Gaussian dependence is appropriate
2. **Nested Archimedean copulas**: Hierarchical structure for grouped variables
3. **Vine copulas**: Build from pairwise bivariate copulas (C-vine, D-vine, R-vine)

For hydrologic applications with 3+ correlated variables (e.g., peak flow, volume, duration), consider using the Multivariate Normal distribution for the initial analysis, which is available in the `Numerics.Distributions` namespace. See the [Multivariate Distributions](multivariate.md) documentation for details.

## Best Practices

1. **Check for dependence**: Use scatter plots and correlation tests before applying copulas. Kendall's tau is preferred over Pearson's $r$ because it is rank-based and invariant under monotonic transformations
2. **Choose copula family**: Match tail behavior to application — Gumbel for joint highs, Clayton for joint lows, Student-t for symmetric tail dependence
3. **Validate fit**: Compare empirical and theoretical joint exceedance curves. Use the pseudo-log-likelihood to compare copula families
4. **Sample size**: Need sufficient data ($n > 50$–$100$) for reliable parameter estimation via Kendall's tau
5. **Non-stationarity**: Check if dependence structure is time-varying — a copula fit to historical data may not represent future conditions

## Limitations

- Assumes marginals are correctly specified — copula inference is only valid when the marginal distributions are well-fitted
- Bivariate copulas may not capture complex nonlinear dependencies in higher dimensions
- Parameter estimation is challenging with limited data, especially for tail dependence
- Tail dependence is difficult to estimate from finite samples — the tail region by definition has few observations

---

## References

<a id="1">[1]</a> R. B. Nelsen, *An Introduction to Copulas*, 2nd ed., New York: Springer, 2006.

<a id="2">[2]</a> H. Joe, *Dependence Modeling with Copulas*, Boca Raton: CRC Press, 2014.

<a id="3">[3]</a> C. Genest and A.-C. Favre, "Everything you always wanted to know about copula modeling but were afraid to ask," *Journal of Hydrologic Engineering*, vol. 12, no. 4, pp. 347-368, 2007.

<a id="4">[4]</a> G. Salvadori, C. De Michele, N. T. Kottegoda and R. Rosso, *Extremes in Nature: An Approach Using Copulas*, Dordrecht: Springer, 2007.

---

[← Previous: Uncertainty Analysis](uncertainty-analysis.md) | [Back to Index](../index.md) | [Next: Multivariate Distributions →](multivariate.md)
