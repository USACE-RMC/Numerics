# Bivariate Copulas

Copulas provide a flexible framework for modeling dependence between random variables, independent of their marginal distributions [[1]](#ref1). This document describes the bivariate copula implementations in ***Numerics***.

## Introduction

### Sklar's Theorem

For any joint distribution function $H(x,y)$ with marginal distributions $F(x)$ and $G(y)$, there exists a copula $C$ such that:

```math
H(x,y) = C(F(x), G(y))
```

Conversely, for any copula $C$ and marginal distributions $F$ and $G$:

```math
C(u,v) = H(F^{-1}(u), G^{-1}(v))
```

This separation of marginals from dependence structure is the fundamental power of copulas.

### Applications

- **Hydrology**: Joint flood peaks and volumes
- **Risk Analysis**: Correlated failure modes
- **Finance**: Portfolio risk modeling
- **Reliability**: Dependent component failures

## Available Copulas

| Copula | Family | Tail Dependence | Parameter Range |
|--------|--------|-----------------|-----------------|
| Normal (Gaussian) | Elliptical | None | $-1 < \rho < 1$ |
| Clayton | Archimedean | Lower | $\theta > 0$ |
| Frank | Archimedean | None | $\theta \neq 0$ |
| Gumbel | Archimedean | Upper | $\theta \geq 1$ |
| Joe | Archimedean | Upper | $\theta \geq 1$ |

---

## Normal (Gaussian) Copula

The Normal copula is derived from the bivariate normal distribution [[2]](#ref2).

### Definition

```math
C_\rho^{Ga}(u,v) = \Phi_\rho(\Phi^{-1}(u), \Phi^{-1}(v))
```

where $\Phi_\rho$ is the bivariate standard normal CDF with correlation $\rho$, and $\Phi^{-1}$ is the standard normal quantile function.

### Density

```math
c_\rho^{Ga}(u,v) = \frac{1}{\sqrt{1-\rho^2}}\exp\left(-\frac{\rho^2(x^2+y^2) - 2\rho xy}{2(1-\rho^2)}\right)
```

where $x = \Phi^{-1}(u)$ and $y = \Phi^{-1}(v)$.

### Properties

| Property | Value |
|----------|-------|
| Parameter | $\rho \in (-1, 1)$ |
| Kendall's τ | $\frac{2}{\pi}\arcsin(\rho)$ |
| Upper tail dependence | 0 |
| Lower tail dependence | 0 |

### Usage

```cs
using Numerics.Distributions.Copulas;

// Create Normal copula with correlation 0.7
var normalCopula = new NormalCopula(0.7);

// Joint probability P(U ≤ 0.8, V ≤ 0.6)
double jointProb = normalCopula.CDF(0.8, 0.6);

// Copula density
double density = normalCopula.PDF(0.5, 0.5);

// Generate correlated uniform samples
double[,] samples = normalCopula.GenerateRandomValues(1000);
```

### Fit to Data

```cs
// Transform data to uniform margins
double[] u = TransformToUniform(xData, fittedDistX);
double[] v = TransformToUniform(yData, fittedDistY);

// Fit copula using MLE or rank correlation
var copula = new NormalCopula();
copula.FitFromData(u, v);
Console.WriteLine($"Estimated ρ: {copula.Rho:F3}");
```

---

## Clayton Copula

An Archimedean copula with lower tail dependence, suitable when extreme low values tend to occur together [[3]](#ref3).

### Definition

```math
C_\theta^{Cl}(u,v) = \left(u^{-\theta} + v^{-\theta} - 1\right)^{-1/\theta}
```

### Generator Function

```math
\phi(t) = \frac{1}{\theta}(t^{-\theta} - 1)
```

### Properties

| Property | Value |
|----------|-------|
| Parameter | $\theta > 0$ |
| Kendall's τ | $\frac{\theta}{\theta + 2}$ |
| Upper tail dependence | 0 |
| Lower tail dependence | $2^{-1/\theta}$ |

### Usage

```cs
// Clayton copula with θ = 2
var clayton = new ClaytonCopula(2.0);

// Lower tail dependence coefficient
double lambdaL = Math.Pow(2, -1.0 / clayton.Theta);
Console.WriteLine($"Lower tail dependence: {lambdaL:F3}");

// Conditional distribution: P(V ≤ v | U = u)
double condProb = clayton.ConditionalCDF(0.1, 0.15);
```

### When to Use

- Variables that tend to have joint extreme **low** values
- Drought analysis (joint low flows)
- Credit risk (joint defaults in economic downturns)

---

## Frank Copula

A symmetric Archimedean copula with no tail dependence [[3]](#ref3).

### Definition

```math
C_\theta^{Fr}(u,v) = -\frac{1}{\theta}\ln\left(1 + \frac{(e^{-\theta u}-1)(e^{-\theta v}-1)}{e^{-\theta}-1}\right)
```

### Properties

| Property | Value |
|----------|-------|
| Parameter | $\theta \neq 0$ |
| Kendall's τ | $1 - \frac{4}{\theta}\left(1 - D_1(\theta)\right)$ |
| Upper tail dependence | 0 |
| Lower tail dependence | 0 |
| Symmetry | Radially symmetric |

where $D_1$ is the first Debye function.

### Usage

```cs
var frank = new FrankCopula(5.0);

// Frank allows negative dependence
var frankNeg = new FrankCopula(-3.0);

// Kendall's tau from parameter
double tau = frank.KendallsTau;
```

### When to Use

- Moderate dependence without extreme co-movements
- When tail dependence is not expected
- General-purpose dependence modeling

---

## Gumbel Copula

An Archimedean copula with upper tail dependence, ideal for modeling joint extreme high values [[3]](#ref3).

### Definition

```math
C_\theta^{Gu}(u,v) = \exp\left(-\left[(-\ln u)^\theta + (-\ln v)^\theta\right]^{1/\theta}\right)
```

### Generator Function

```math
\phi(t) = (-\ln t)^\theta
```

### Properties

| Property | Value |
|----------|-------|
| Parameter | $\theta \geq 1$ |
| Kendall's τ | $1 - \frac{1}{\theta}$ |
| Upper tail dependence | $2 - 2^{1/\theta}$ |
| Lower tail dependence | 0 |

### Usage

```cs
// Gumbel copula with θ = 3
var gumbel = new GumbelCopula(3.0);

// Upper tail dependence
double lambdaU = 2 - Math.Pow(2, 1.0 / gumbel.Theta);
Console.WriteLine($"Upper tail dependence: {lambdaU:F3}");

// Joint probability of extreme values
double u = 0.99;  // 99th percentile of X
double v = 0.99;  // 99th percentile of Y
double jointExtreme = gumbel.CCDF(u, v);  // P(U > u, V > v)
```

### When to Use

- Variables with joint extreme **high** values
- Flood peaks and volumes
- Insurance losses in catastrophic events
- Maximum temperatures

---

## Joe Copula

Similar to Gumbel but with stronger upper tail dependence [[3]](#ref3).

### Definition

```math
C_\theta^{Jo}(u,v) = 1 - \left[(1-u)^\theta + (1-v)^\theta - (1-u)^\theta(1-v)^\theta\right]^{1/\theta}
```

### Properties

| Property | Value |
|----------|-------|
| Parameter | $\theta \geq 1$ |
| Kendall's τ | $1 + \frac{4}{\theta^2}\int_0^1 t\ln(t)(1-t)^{2(1-\theta)/\theta}dt$ |
| Upper tail dependence | $2 - 2^{1/\theta}$ |
| Lower tail dependence | 0 |

### Usage

```cs
var joe = new JoeCopula(2.5);

// Compare upper tail dependence with Gumbel
var gumbel = new GumbelCopula(2.5);
Console.WriteLine($"Joe λU:    {joe.UpperTailDependence:F3}");
Console.WriteLine($"Gumbel λU: {gumbel.UpperTailDependence:F3}");
```

---

## Working with Copulas

### Fitting Copulas to Data

```cs
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

// Original bivariate data
double[] xData = { /* ... */ };
double[] yData = { /* ... */ };

// Step 1: Fit marginal distributions
var distX = new GeneralizedExtremeValue();
distX.SetParameters(distX.ParametersFromLinearMoments(xData));

var distY = new GeneralizedExtremeValue();
distY.SetParameters(distY.ParametersFromLinearMoments(yData));

// Step 2: Transform to uniform margins (probability integral transform)
double[] u = xData.Select(x => distX.CDF(x)).ToArray();
double[] v = yData.Select(y => distY.CDF(y)).ToArray();

// Step 3: Fit copula
var copula = new GumbelCopula();
copula.FitFromData(u, v);

Console.WriteLine($"Estimated θ: {copula.Theta:F3}");
Console.WriteLine($"Kendall's τ: {copula.KendallsTau:F3}");
```

### Computing Joint Probabilities

```cs
// P(X ≤ x AND Y ≤ y)
double x = 1500;
double y = 2000;
double jointProb = copula.CDF(distX.CDF(x), distY.CDF(y));

// P(X > x AND Y > y) - joint exceedance
double jointExceed = copula.SurvivalFunction(distX.CDF(x), distY.CDF(y));

// P(X > x OR Y > y) - at least one exceeds
double eitherExceed = 1 - copula.CDF(distX.CDF(x), distY.CDF(y));

Console.WriteLine($"P(X ≤ {x} AND Y ≤ {y}) = {jointProb:F4}");
Console.WriteLine($"P(X > {x} AND Y > {y}) = {jointExceed:F4}");
Console.WriteLine($"P(X > {x} OR Y > {y})  = {eitherExceed:F4}");
```

### Conditional Distributions

```cs
// P(Y ≤ y | X = x)
double x = 1500;
double y = 2000;
double u = distX.CDF(x);
double v = distY.CDF(y);

double condProb = copula.ConditionalCDF_V_given_U(u, v);
Console.WriteLine($"P(Y ≤ {y} | X = {x}) = {condProb:F4}");

// Conditional quantile: y such that P(Y ≤ y | X = x) = p
double p = 0.95;
double vCond = copula.ConditionalInverseCDF_V_given_U(u, p);
double yCond = distY.InverseCDF(vCond);
Console.WriteLine($"95th percentile of Y given X = {x}: {yCond:F0}");
```

### Generating Correlated Samples

```cs
// Generate correlated samples
int nSamples = 10000;
double[,] uvSamples = copula.GenerateRandomValues(nSamples);

// Transform to original scale
double[] xSamples = new double[nSamples];
double[] ySamples = new double[nSamples];

for (int i = 0; i < nSamples; i++)
{
    xSamples[i] = distX.InverseCDF(uvSamples[i, 0]);
    ySamples[i] = distY.InverseCDF(uvSamples[i, 1]);
}
```

---

## Copula Selection

### Using Tail Dependence

| Data Characteristic | Recommended Copula |
|--------------------|-------------------|
| No tail dependence | Normal, Frank |
| Lower tail dependence | Clayton |
| Upper tail dependence | Gumbel, Joe |
| Symmetric dependence | Normal, Frank |

### Goodness-of-Fit

```cs
// Compare copulas using AIC
var copulas = new IBivariateCopula[]
{
    new NormalCopula(),
    new ClaytonCopula(),
    new FrankCopula(),
    new GumbelCopula(),
    new JoeCopula()
};

Console.WriteLine("Copula       θ        AIC");
Console.WriteLine("------       -        ---");

foreach (var cop in copulas)
{
    cop.FitFromData(u, v);
    double ll = cop.LogLikelihood(u, v);
    double aic = -2 * ll + 2 * 1;  // 1 parameter
    
    Console.WriteLine($"{cop.GetType().Name,-12} {cop.Theta:F3}   {aic:F1}");
}
```

---

## Example: Flood Peak and Volume Analysis

```cs
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

// Annual maximum flood peaks (cfs) and volumes (acre-feet)
double[] peaks = { 15200, 18900, 12500, 22100, 16800, 19500, 14200 };
double[] volumes = { 45000, 62000, 38000, 71000, 52000, 58000, 41000 };

// Fit marginal distributions
var peakDist = new GeneralizedExtremeValue();
peakDist.SetParameters(peakDist.ParametersFromLinearMoments(peaks));

var volDist = new GeneralizedExtremeValue();
volDist.SetParameters(volDist.ParametersFromLinearMoments(volumes));

// Transform to uniform margins
double[] u = peaks.Select(p => peakDist.CDF(p)).ToArray();
double[] v = volumes.Select(vol => volDist.CDF(vol)).ToArray();

// Fit Gumbel copula (appropriate for upper tail dependence)
var copula = new GumbelCopula();
copula.FitFromData(u, v);

// Design event: 100-year flood
double peak100 = peakDist.InverseCDF(0.99);
double vol100 = volDist.InverseCDF(0.99);

// Joint probability analysis
double uDesign = 0.99;
double vDesign = 0.99;

// P(Peak > peak100 AND Volume > vol100)
double jointExceed = 1 - uDesign - vDesign + copula.CDF(uDesign, vDesign);

// OR-return period: P(Peak > peak100 OR Volume > vol100)
double orExceed = 1 - copula.CDF(uDesign, vDesign);
double orReturnPeriod = 1 / orExceed;

// AND-return period
double andReturnPeriod = 1 / jointExceed;

Console.WriteLine($"100-year peak:   {peak100:N0} cfs");
Console.WriteLine($"100-year volume: {vol100:N0} acre-ft");
Console.WriteLine($"OR-return period:  {orReturnPeriod:F0} years");
Console.WriteLine($"AND-return period: {andReturnPeriod:F0} years");
```

---

## References

<a id="ref1">[1]</a> Nelsen, R. B. (2006). *An Introduction to Copulas* (2nd ed.). Springer.

<a id="ref2">[2]</a> Joe, H. (1997). *Multivariate Models and Dependence Concepts*. Chapman & Hall.

<a id="ref3">[3]</a> Genest, C., & Favre, A.-C. (2007). Everything you always wanted to know about copula modeling but were afraid to ask. *Journal of Hydrologic Engineering*, 12(4), 347-368.

<a id="ref4">[4]</a> Salvadori, G., De Michele, C., Kottegoda, N. T., & Rosso, R. (2007). *Extremes in Nature: An Approach Using Copulas*. Springer.

<a id="ref5">[5]</a> Salvadori, G., & De Michele, C. (2004). Frequency analysis via copulas: Theoretical aspects and applications to hydrological events. *Water Resources Research*, 40(12).
