# Copulas

[← Previous: Uncertainty Analysis](uncertainty-analysis.md) | [Back to Index](../index.md)

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

### Gaussian Copula

Models linear correlation with normal dependence structure:

```cs
using Numerics.Distributions.Copulas;

// Correlation matrix for 2D Gaussian copula
double rho = 0.7;  // Correlation coefficient
var corrMatrix = new double[,] { { 1.0, rho }, { rho, 1.0 } };

var gaussianCopula = new GaussianCopula(corrMatrix);

// Evaluate copula density
double u1 = 0.3, u2 = 0.7;
double density = gaussianCopula.PDF(new double[] { u1, u2 });

Console.WriteLine($"Gaussian copula density: {density:F4}");
```

### Student's t Copula

Similar to Gaussian but with tail dependence:

```cs
// t-copula with 5 degrees of freedom
int nu = 5;
var tCopula = new StudentTCopula(corrMatrix, nu);

double density = tCopula.PDF(new double[] { u1, u2 });

Console.WriteLine($"t-copula density: {density:F4}");
Console.WriteLine("t-copula has stronger tail dependence than Gaussian");
```

### Archimedean Copulas

Family of copulas with specific dependence structures:

```cs
// Clayton copula (lower tail dependence)
double theta = 2.0;  // Dependence parameter
var claytonCopula = new ClaytonCopula(theta);

// Gumbel copula (upper tail dependence)
var gumbelCopula = new GumbelCopula(theta);

// Frank copula (no tail dependence)
var frankCopula = new FrankCopula(theta);
```

## Practical Example: Bivariate Distribution

Construct a bivariate distribution with arbitrary marginals and specified dependence:

```cs
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

// Step 1: Define marginal distributions
var margin1 = new LogNormal(4.0, 0.5);  // Streamflow
var margin2 = new Gumbel(100, 20);      // Peak stage

// Step 2: Define dependence via copula
double rho = 0.8;  // Strong positive correlation
var corrMatrix = new double[,] { { 1.0, rho }, { rho, 1.0 } };
var copula = new GaussianCopula(corrMatrix);

// Step 3: Sample from joint distribution
int n = 1000;
var samples = copula.Sample(n);

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
double jointExceedance = copula.Survival(new double[] { u1, u2 });

Console.WriteLine($"Joint exceedance probability: {jointExceedance:E4}");
Console.WriteLine($"Return period: {1.0 / jointExceedance:F1} years");
```

### Conditional Distributions

Given flow, what is the conditional distribution of stage?

```cs
// Observed flow
double observedFlow = 12000;
double uFlow = margin1.CDF(observedFlow);

// Conditional CDF for stage | flow
Func<double, double> conditionalCDF = (stage) =>
{
    double uStage = margin2.CDF(stage);
    return copula.ConditionalCDF(uFlow, uStage, conditionIndex: 0);
};

// Find conditional quantiles
double[] condProbs = { 0.5, 0.9, 0.95 };
Console.WriteLine($"Given flow = {observedFlow:F0} cfs:");
foreach (var p in condProbs)
{
    // This would require numerical solution
    Console.WriteLine($"  {p:P0} quantile of stage");
}
```

## Tail Dependence

Different copulas have different tail dependence properties:

```cs
// Gaussian: No tail dependence (λ_L = λ_U = 0)
// t-copula: Symmetric tail dependence
// Clayton: Lower tail dependence only
// Gumbel: Upper tail dependence only
// Frank: No tail dependence

Console.WriteLine("Tail Dependence Properties:");
Console.WriteLine("  Gaussian: No tail dependence");
Console.WriteLine("  Student-t: Symmetric tail dependence");
Console.WriteLine("  Clayton: Lower tail dependence (joint lows)");
Console.WriteLine("  Gumbel: Upper tail dependence (joint highs)");
Console.WriteLine("  Frank: No tail dependence");

// For flood analysis: Gumbel copula captures joint extremes
// For drought analysis: Clayton copula captures joint lows
```

## Fitting Copulas to Data

```cs
double[] x = { /* observed data series 1 */ };
double[] y = { /* observed data series 2 */ };

// Step 1: Transform to uniform margins (pseudo-observations)
var u = x.Select(xi => (double)Array.FindIndex(x.OrderBy(v => v).ToArray(), v => v == xi) / x.Length);
var v = y.Select(yi => (double)Array.FindIndex(y.OrderBy(w => w).ToArray(), w => w == yi) / y.Length);

// Step 2: Fit copula to pseudo-observations
// Use maximum likelihood or rank correlation methods

// Step 3: Select best copula using AIC/BIC
var candidates = new[] { "Gaussian", "t", "Clayton", "Gumbel", "Frank" };

Console.WriteLine("Fit each candidate copula and select best by AIC");
```

## Vine Copulas

For higher dimensions, vine copulas decompose multivariate dependence:

```cs
// C-vine, D-vine, and R-vine structures available
// Allow flexible modeling of high-dimensional dependence
// Construct from pairwise bivariate copulas

Console.WriteLine("Vine copulas enable flexible high-dimensional modeling");
Console.WriteLine("Use for systems with > 2 correlated variables");
```

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

[← Previous: Uncertainty Analysis](uncertainty-analysis.md) | [Back to Index](../index.md)
