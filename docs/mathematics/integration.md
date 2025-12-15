# Numerical Integration

The ***Numerics*** library provides comprehensive numerical integration methods for one-dimensional and multi-dimensional problems, ranging from simple quadrature rules to adaptive algorithms and Monte Carlo methods.

## Overview

| Method | Dimensions | Best For |
|--------|------------|----------|
| Trapezoidal Rule | 1D | Quick estimates, periodic functions |
| Simpson's Rule | 1D | Smooth functions |
| Adaptive Simpson's | 1D | General purpose, automatic accuracy |
| Gauss-Legendre | 1D | Smooth functions, high accuracy |
| Gauss-Kronrod | 1D | Automatic error estimation |
| Gauss-Lobatto | 1D | Functions with endpoint information |
| Monte Carlo | N-D | High dimensions, irregular domains |
| MISER | N-D | Stratified Monte Carlo |
| VEGAS | N-D | Adaptive importance sampling |

---

## One-Dimensional Integration

### Trapezoidal Rule

Approximates the integral using linear interpolation [[1]](#ref1):

```math
\int_a^b f(x)\,dx \approx \frac{h}{2}\left[f(a) + 2\sum_{i=1}^{n-1}f(x_i) + f(b)\right]
```

where $h = (b-a)/n$.

```cs
using Numerics.Mathematics.Integration;

Func<double, double> f = x => Math.Exp(-x * x);

var trap = new TrapezoidalRule(f, -3, 3, intervals: 100);
trap.Integrate();

Console.WriteLine($"Result: {trap.Result:F10}");
Console.WriteLine($"Function evaluations: {trap.FunctionEvaluations}");
```

**Error**: $O(h^2)$ - second-order accurate.

### Simpson's Rule

Uses quadratic interpolation for higher accuracy [[1]](#ref1):

```math
\int_a^b f(x)\,dx \approx \frac{h}{3}\left[f(a) + 4\sum_{i=1,3,5...}f(x_i) + 2\sum_{i=2,4,6...}f(x_i) + f(b)\right]
```

```cs
var simpson = new SimpsonsRule(f, -3, 3, intervals: 100);
simpson.Integrate();

Console.WriteLine($"Result: {simpson.Result:F10}");
```

**Error**: $O(h^4)$ - fourth-order accurate.

### Adaptive Simpson's Rule

Automatically subdivides intervals to achieve target accuracy:

```cs
var adaptive = new AdaptiveSimpsonsRule(f, -3, 3);
adaptive.RelativeTolerance = 1e-10;
adaptive.MaxIterations = 1000;
adaptive.Integrate();

Console.WriteLine($"Result: {adaptive.Result:F12}");
Console.WriteLine($"Estimated error: {adaptive.Error:E2}");
Console.WriteLine($"Function evaluations: {adaptive.FunctionEvaluations}");
```

**Algorithm**: Recursively bisects intervals where local error exceeds tolerance.

### Gauss-Legendre Quadrature

Optimal polynomial quadrature for smooth functions [[1]](#ref1):

```math
\int_{-1}^{1} f(x)\,dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

where $x_i$ are roots of Legendre polynomials and $w_i$ are corresponding weights.

```cs
var gl = new GaussLegendre(f, -3, 3, order: 20);
gl.Integrate();

Console.WriteLine($"Result: {gl.Result:F12}");
```

**Error**: Exact for polynomials of degree $2n-1$ or less.

### Adaptive Gauss-Kronrod (G10K21)

Pairs a 10-point Gauss rule with a 21-point Kronrod extension for error estimation [[2]](#ref2):

```cs
var gk = new AdaptiveGaussKronrod(f, -3, 3);
gk.RelativeTolerance = 1e-12;
gk.Integrate();

Console.WriteLine($"Result: {gk.Result:F14}");
Console.WriteLine($"Error estimate: {gk.Error:E2}");
```

**Advantages**: Efficient error estimation, reuses function evaluations.

### Adaptive Gauss-Lobatto

Includes endpoints in quadrature nodes, useful for endpoint singularities [[2]](#ref2):

```cs
var gl = new AdaptiveGaussLobatto(f, 0, 1);
gl.RelativeTolerance = 1e-10;
gl.Integrate();

Console.WriteLine($"Result: {gl.Result:F12}");
```

---

## Handling Difficult Integrands

### Singularities

For integrands with singularities, use transformations or specialized methods:

```cs
// Integral of 1/√x from 0 to 1 (singularity at x=0)
Func<double, double> singular = x => 1.0 / Math.Sqrt(x);

// Use substitution: let u = √x, then dx = 2u du
Func<double, double> transformed = u => 2.0;  // 1/√(u²) * 2u = 2

var integrator = new AdaptiveGaussKronrod(transformed, 0, 1);
integrator.Integrate();
Console.WriteLine($"∫₀¹ 1/√x dx = {integrator.Result:F10}");  // Should be 2
```

### Infinite Intervals

Transform infinite intervals to finite:

```cs
// Integral of e^(-x²) from 0 to ∞
// Use substitution: x = t/(1-t), dx = 1/(1-t)² dt

Func<double, double> infiniteIntegrand = x => Math.Exp(-x * x);

Func<double, double> transformed = t =>
{
    if (t >= 1) return 0;
    double x = t / (1 - t);
    double jacobian = 1.0 / ((1 - t) * (1 - t));
    return infiniteIntegrand(x) * jacobian;
};

var integrator = new AdaptiveGaussKronrod(transformed, 0, 1);
integrator.Integrate();
Console.WriteLine($"∫₀^∞ e^(-x²) dx = {integrator.Result:F10}");  // √π/2 ≈ 0.8862
```

### Oscillatory Integrands

For highly oscillatory functions, use more quadrature points or specialized methods:

```cs
// Integral of sin(100x) from 0 to π
Func<double, double> oscillatory = x => Math.Sin(100 * x);

var integrator = new AdaptiveGaussKronrod(oscillatory, 0, Math.PI);
integrator.RelativeTolerance = 1e-8;
integrator.MaxIterations = 10000;
integrator.Integrate();

Console.WriteLine($"Result: {integrator.Result:F10}");
```

---

## Two-Dimensional Integration

### Adaptive Simpson's 2D

Extends adaptive Simpson's rule to rectangular domains:

```cs
Func<double, double, double> f2d = (x, y) => Math.Exp(-(x*x + y*y));

var simpson2d = new AdaptiveSimpsonsRule2D(f2d, -3, 3, -3, 3);
simpson2d.RelativeTolerance = 1e-8;
simpson2d.Integrate();

Console.WriteLine($"Result: {simpson2d.Result:F10}");  // π ≈ 3.14159
```

### Product Rules

For separable integrands, use product of 1D rules:

```cs
// ∫∫ f(x)g(y) dx dy = (∫f(x)dx)(∫g(y)dy)
Func<double, double> fx = x => Math.Exp(-x * x);
Func<double, double> gy = y => Math.Exp(-y * y);

var intX = new AdaptiveGaussKronrod(fx, -3, 3);
var intY = new AdaptiveGaussKronrod(gy, -3, 3);

intX.Integrate();
intY.Integrate();

double result = intX.Result * intY.Result;
```

---

## Multi-Dimensional Integration

### Monte Carlo Integration

For high-dimensional integrals, Monte Carlo methods become essential [[3]](#ref3):

```math
\int_\Omega f(\mathbf{x})\,d\mathbf{x} \approx V \cdot \frac{1}{N}\sum_{i=1}^{N}f(\mathbf{x}_i)
```

where $V$ is the volume of the integration domain.

```cs
// 5-dimensional integral
Func<double[], double> f5d = x => 
{
    double sum = 0;
    for (int i = 0; i < x.Length; i++)
        sum += x[i] * x[i];
    return Math.Exp(-sum);
};

double[] lower = { -2, -2, -2, -2, -2 };
double[] upper = { 2, 2, 2, 2, 2 };

var mc = new MonteCarlo(f5d, lower, upper);
mc.Iterations = 100000;
mc.Integrate();

Console.WriteLine($"Result: {mc.Result:F6}");
Console.WriteLine($"Error estimate: {mc.Error:F6}");
```

**Error**: $O(1/\sqrt{N})$ - independent of dimension.

### MISER (Recursive Stratified Sampling)

Improves Monte Carlo by concentrating samples in high-variance regions [[3]](#ref3):

```cs
var miser = new MISER(f5d, lower, upper);
miser.Iterations = 100000;
miser.MinIterationsPerRegion = 100;
miser.Integrate();

Console.WriteLine($"Result: {miser.Result:F6}");
Console.WriteLine($"Error estimate: {miser.Error:F6}");
```

### VEGAS (Adaptive Importance Sampling)

Builds an adaptive importance sampling grid [[4]](#ref4):

```cs
var vegas = new Vegas(f5d, lower, upper);
vegas.Iterations = 10000;
vegas.WarmupIterations = 5000;
vegas.NumberOfBins = 50;
vegas.Integrate();

Console.WriteLine($"Result: {vegas.Result:F8}");
Console.WriteLine($"Error estimate: {vegas.Error:E2}");
Console.WriteLine($"Chi-squared/dof: {vegas.ChiSquaredPerDof:F2}");
```

**Advantages**: Excellent for peaked integrands, provides convergence diagnostics.

---

## Probability Space Integration

For risk analysis, integrate over probability distributions:

```cs
using Numerics.Distributions;
using Numerics.Mathematics.Integration;

// Expected value E[g(X)] where X ~ Normal(0,1)
var normal = new Normal(0, 1);

Func<double[], double> expectation = p =>
{
    // Transform from [0,1] to real line via inverse CDF
    double x = normal.InverseCDF(p[0]);
    
    // Function of interest
    double gx = x * x;  // E[X²] = 1 for standard normal
    
    return gx;  // Weight is implicitly 1 in probability space
};

var mc = new MonteCarlo(expectation, new[] { 0.0 }, new[] { 1.0 });
mc.Iterations = 100000;
mc.Integrate();

Console.WriteLine($"E[X²] = {mc.Result:F6}");  // Should be ≈ 1
```

### Failure Probability Integration

```cs
// P(g(X) < 0) where g is a limit state function
Func<double[], double> failureIndicator = p =>
{
    double x1 = new Normal(0, 1).InverseCDF(p[0]);
    double x2 = new Normal(0, 1).InverseCDF(p[1]);
    
    // Limit state function
    double g = 5 - x1 - x2;
    
    // Indicator function
    return g < 0 ? 1.0 : 0.0;
};

var vegas = new Vegas(failureIndicator, new[] { 0.0, 0.0 }, new[] { 1.0, 1.0 });
vegas.Iterations = 50000;
vegas.Integrate();

Console.WriteLine($"Failure probability: {vegas.Result:E4}");
```

---

## Choosing an Integration Method

| Scenario | Recommended Method |
|----------|-------------------|
| 1D, smooth function | Gauss-Kronrod |
| 1D, unknown smoothness | Adaptive Simpson's |
| 1D, endpoint singularity | Gauss-Lobatto |
| 1D, high accuracy needed | Gauss-Kronrod with tight tolerance |
| 2D, rectangular domain | Adaptive Simpson's 2D |
| 3D-6D, smooth | Product rules or VEGAS |
| High-D (>6), any | Monte Carlo or VEGAS |
| Peaked integrand | VEGAS |
| Probability integrals | Monte Carlo in probability space |

---

## Performance Comparison

```cs
Func<double, double> testFunc = x => Math.Exp(-x * x) * Math.Cos(x);
double trueValue = 1.3803884470431430;  // Known analytical result

var methods = new (string Name, Integrator Int)[]
{
    ("Trapezoidal (n=100)", new TrapezoidalRule(testFunc, -5, 5, 100)),
    ("Simpson's (n=100)", new SimpsonsRule(testFunc, -5, 5, 100)),
    ("Adaptive Simpson's", new AdaptiveSimpsonsRule(testFunc, -5, 5)),
    ("Gauss-Legendre (n=20)", new GaussLegendre(testFunc, -5, 5, 20)),
    ("Gauss-Kronrod", new AdaptiveGaussKronrod(testFunc, -5, 5))
};

Console.WriteLine("Method                  Result          Error       Evals");
Console.WriteLine("------                  ------          -----       -----");

foreach (var (name, integrator) in methods)
{
    integrator.Integrate();
    double error = Math.Abs(integrator.Result - trueValue);
    Console.WriteLine($"{name,-22}  {integrator.Result:F10}  {error:E2}  {integrator.FunctionEvaluations,5}");
}
```

---

## References

<a id="ref1">[1]</a> Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

<a id="ref2">[2]</a> Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W., & Kahaner, D. K. (1983). *QUADPACK: A Subroutine Package for Automatic Integration*. Springer.

<a id="ref3">[3]</a> Press, W. H., & Farrar, G. R. (1990). Recursive stratified sampling for multidimensional Monte Carlo integration. *Computers in Physics*, 4(2), 190-195.

<a id="ref4">[4]</a> Lepage, G. P. (1978). A new algorithm for adaptive multidimensional integration. *Journal of Computational Physics*, 27(2), 192-203.
