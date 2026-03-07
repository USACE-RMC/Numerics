# Copulas

[← Previous: Uncertainty Analysis](uncertainty-analysis.md) | [Back to Index](../index.md) | [Next: Multivariate Distributions →](multivariate.md)

Copulas separate the dependence structure of multivariate distributions from their marginal distributions. The ***Numerics*** library provides copula functions for modeling dependence between random variables in risk assessment and multivariate analysis [[1]](#1).

## Overview

A copula is a multivariate distribution with uniform marginals on [0,1]. For any multivariate distribution with marginals F₁, F₂, ..., Fₙ, there exists a copula C such that:

```
F(x₁, x₂, ..., xₙ) = C(F₁(x₁), F₂(x₂), ..., Fₙ(xₙ))
```

This allows us to:
1. Model marginal distributions independently
2. Model dependence separately via copula
3. Combine them to form joint distribution

## Available Copulas

The library provides common copula families:

### Normal (Gaussian) Copula

Models linear correlation with normal dependence structure:

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

### Student's t Copula

Similar to Gaussian but with tail dependence:

```cs
// t-copula with 5 degrees of freedom
int nu = 5;
var tCopula = new StudentTCopula(rho, nu);

double density = tCopula.PDF(u1, u2);

Console.WriteLine($"t-copula density: {density:F4}");
Console.WriteLine("t-copula has stronger tail dependence than Gaussian");
```

### Archimedean Copulas

Family of copulas defined through generator functions, each with distinct dependence structures:

```cs
// Clayton copula (lower tail dependence), θ ∈ (0, ∞)
var claytonCopula = new ClaytonCopula(2.0);

// Gumbel copula (upper tail dependence), θ ∈ [1, ∞)
var gumbelCopula = new GumbelCopula(2.0);

// Frank copula (no tail dependence), θ ∈ (-∞, ∞) \ {0}
var frankCopula = new FrankCopula(5.0);

// Joe copula (upper tail dependence), θ ∈ [1, ∞)
var joeCopula = new JoeCopula(2.0);

// Ali-Mikhail-Haq copula (weak dependence), θ ∈ [-1, 1]
var amhCopula = new AMHCopula(0.5);
```

**Copula selection guide:**

| Copula | Tail Dependence | Parameter Range | Best For |
|--------|----------------|-----------------|----------|
| Clayton | Lower tail | θ ∈ (0, ∞) | Joint low extremes (droughts) |
| Gumbel | Upper tail | θ ∈ [1, ∞) | Joint high extremes (floods) |
| Frank | None | θ ∈ (-∞, ∞) \ {0} | Moderate symmetric dependence |
| Joe | Upper tail | θ ∈ [1, ∞) | Strong upper tail dependence |
| AMH | None | θ ∈ [-1, 1] | Weak dependence structures |

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
double jointExceedance = copula.ANDJointExceedanceProbability(u1, u2);

Console.WriteLine($"Joint exceedance probability: {jointExceedance:E4}");
Console.WriteLine($"Return period: {1.0 / jointExceedance:F1} years");
```

### Conditional Distributions

Given flow, what is the conditional distribution of stage? The conditional CDF can be computed numerically using the copula CDF via partial differentiation: C(v|u) = dC(u,v)/du.

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

Different copulas have different tail dependence properties:

```cs
// Gaussian: No tail dependence (λ_L = λ_U = 0)
// t-copula: Symmetric tail dependence
// Clayton: Lower tail dependence only
// Gumbel: Upper tail dependence only
// Joe: Upper tail dependence only
// Frank: No tail dependence
// AMH: No tail dependence

Console.WriteLine("Tail Dependence Properties:");
Console.WriteLine("  Gaussian: No tail dependence");
Console.WriteLine("  Student-t: Symmetric tail dependence");
Console.WriteLine("  Clayton: Lower tail dependence (joint lows)");
Console.WriteLine("  Gumbel: Upper tail dependence (joint highs)");
Console.WriteLine("  Joe: Upper tail dependence (joint highs)");
Console.WriteLine("  Frank: No tail dependence");
Console.WriteLine("  AMH: No tail dependence");

// For flood analysis: Gumbel copula captures joint extremes
// For drought analysis: Clayton copula captures joint lows
```

## Fitting Copulas to Data

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

// Step 3: Estimate copula parameters using rank correlation
double tau = Correlation.KendallsTau(x, y);
Console.WriteLine($"Kendall's tau: {tau:F3}");

// Different copulas have different tau-to-parameter relationships
// Clayton: θ = 2τ/(1-τ) for τ > 0
// Gumbel: θ = 1/(1-τ) for τ > 0
```

## Higher-Dimensional Dependence

For problems with more than two variables, there are several approaches:

1. **Multivariate Normal distribution**: Use when Gaussian dependence is appropriate
2. **Nested Archimedean copulas**: Hierarchical structure for grouped variables
3. **Vine copulas**: Build from pairwise bivariate copulas (C-vine, D-vine, R-vine)

For hydrologic applications with 3+ correlated variables (e.g., peak flow, volume, duration), consider using the Multivariate Normal distribution for the initial analysis, which is available in the `Numerics.Distributions` namespace. See the [Multivariate Distributions](multivariate.md) documentation for details.

## Best Practices

1. **Check for dependence**: Use scatter plots and correlation tests before applying copulas
2. **Choose copula family**: Match tail behavior to application
   - Joint extremes → Gumbel
   - Joint lows → Clayton  
   - Moderate correlation → Gaussian
   - Heavy tails → Student-t
3. **Validate fit**: Check if copula captures observed dependence structure
4. **Sample size**: Need sufficient data (n > 50-100) for reliable fitting
5. **Non-stationarity**: Check if dependence structure is time-varying

## Limitations

- Assumes marginals are correctly specified
- May not capture complex nonlinear dependencies
- Parameter estimation challenging with limited data
- Tail dependence difficult to estimate from finite samples

---

## References

<a id="1">[1]</a> Nelsen, R. B. (2006). *An Introduction to Copulas* (2nd ed.). Springer.

<a id="2">[2]</a> Joe, H. (2014). *Dependence Modeling with Copulas*. CRC Press.

---

[← Previous: Uncertainty Analysis](uncertainty-analysis.md) | [Back to Index](../index.md) | [Next: Multivariate Distributions →](multivariate.md)
