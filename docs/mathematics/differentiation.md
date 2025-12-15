# Numerical Differentiation

The `NumericalDerivative` class provides methods for approximating derivatives using finite difference formulas, including basic differences, Ridder's extrapolation method, and multivariate calculus operations.

## Overview

| Method | Accuracy | Best For |
|--------|----------|----------|
| Forward Difference | $O(h)$ | One-sided derivatives |
| Backward Difference | $O(h)$ | One-sided derivatives |
| Central Difference | $O(h^2)$ | General use |
| Ridder's Method | $O(h^{10})$ | High accuracy |

---

## Basic Finite Differences

### Forward Difference

First-order accurate approximation [[1]](#ref1):

```math
f'(x) \approx \frac{f(x+h) - f(x)}{h} + O(h)
```

```cs
using Numerics.Mathematics;

Func<double, double> f = x => Math.Sin(x);
double x = Math.PI / 4;
double h = 1e-6;

double forward = NumericalDerivative.ForwardDifference(f, x, h);
double exact = Math.Cos(x);

Console.WriteLine($"Forward difference: {forward:F10}");
Console.WriteLine($"Exact derivative:   {exact:F10}");
Console.WriteLine($"Error: {Math.Abs(forward - exact):E2}");
```

### Backward Difference

```math
f'(x) \approx \frac{f(x) - f(x-h)}{h} + O(h)
```

```cs
double backward = NumericalDerivative.BackwardDifference(f, x, h);
```

### Central Difference

Second-order accurate, generally preferred [[1]](#ref1):

```math
f'(x) \approx \frac{f(x+h) - f(x-h)}{2h} + O(h^2)
```

```cs
double central = NumericalDerivative.CentralDifference(f, x, h);
Console.WriteLine($"Central difference: {central:F10}");
```

### Optimal Step Size

The optimal step size balances truncation and roundoff errors:

```math
h_{opt} \approx \sqrt{\epsilon_m} \cdot |x|
```

where $\epsilon_m$ is machine epsilon (~2.2e-16 for double precision).

```cs
// Automatic step size selection
double derivative = NumericalDerivative.Derivative(f, x);
```

---

## Second Derivatives

### Central Difference for Second Derivative

```math
f''(x) \approx \frac{f(x+h) - 2f(x) + f(x-h)}{h^2} + O(h^2)
```

```cs
Func<double, double> f = x => Math.Exp(x);
double x = 1.0;

double second = NumericalDerivative.SecondDerivative(f, x);
double exact = Math.Exp(x);  // f''(x) = e^x

Console.WriteLine($"Numerical: {second:F10}");
Console.WriteLine($"Exact:     {exact:F10}");
```

---

## Ridder's Method

Ridder's extrapolation achieves high accuracy by combining estimates at different step sizes [[2]](#ref2):

```cs
Func<double, double> f = x => Math.Exp(x);
double x = 1.0;

var result = NumericalDerivative.Ridders(f, x, initialStepSize: 0.1);

Console.WriteLine($"Derivative: {result.Derivative:F12}");
Console.WriteLine($"Error estimate: {result.Error:E2}");
```

### Algorithm

1. Compute derivatives at step sizes $h, h/2, h/4, ...$
2. Use Neville's algorithm to extrapolate to $h \to 0$
3. Return when error estimate stops decreasing

### When to Use Ridder's Method

- When high accuracy is required
- When the function is smooth (infinitely differentiable)
- When function evaluations are not expensive

**Note**: Ridder's method is less effective for non-smooth functions or when $f$ contains numerical noise.

---

## Multivariate Derivatives

### Gradient

The gradient is the vector of partial derivatives [[1]](#ref1):

```math
\nabla f = \left[\frac{\partial f}{\partial x_1}, \frac{\partial f}{\partial x_2}, \ldots, \frac{\partial f}{\partial x_n}\right]
```

```cs
Func<double[], double> f = x => x[0]*x[0] + 2*x[1]*x[1] + x[0]*x[1];
double[] point = { 1.0, 2.0 };

double[] gradient = NumericalDerivative.Gradient(f, point);

Console.WriteLine($"∂f/∂x₁ = {gradient[0]:F6}");  // 2x₁ + x₂ = 4
Console.WriteLine($"∂f/∂x₂ = {gradient[1]:F6}");  // 4x₂ + x₁ = 9
```

### Gradient with Bounds

For optimization with bounds, use one-sided differences near boundaries:

```cs
double[] lowerBounds = { 0, 0 };
double[] upperBounds = { 10, 10 };

double[] gradient = NumericalDerivative.Gradient(f, point, lowerBounds, upperBounds);
```

### Jacobian Matrix

For vector-valued functions $\mathbf{f}: \mathbb{R}^n \to \mathbb{R}^m$:

```math
J_{ij} = \frac{\partial f_i}{\partial x_j}
```

```cs
Func<double[], double[]> vectorFunc = x => new double[]
{
    x[0] * x[0] + x[1],      // f₁
    x[0] * x[1] - x[1] * x[1] // f₂
};

double[] point = { 2.0, 3.0 };
double[,] jacobian = NumericalDerivative.Jacobian(vectorFunc, point, outputDim: 2);

Console.WriteLine("Jacobian:");
Console.WriteLine($"  [{jacobian[0,0]:F4}, {jacobian[0,1]:F4}]");
Console.WriteLine($"  [{jacobian[1,0]:F4}, {jacobian[1,1]:F4}]");
```

### Hessian Matrix

The Hessian is the matrix of second partial derivatives:

```math
H_{ij} = \frac{\partial^2 f}{\partial x_i \partial x_j}
```

```cs
Func<double[], double> f = x => x[0]*x[0] + 2*x[1]*x[1] + x[0]*x[1];
double[] point = { 1.0, 2.0 };

double[,] hessian = NumericalDerivative.Hessian(f, point);

Console.WriteLine("Hessian:");
Console.WriteLine($"  [{hessian[0,0]:F4}, {hessian[0,1]:F4}]");
Console.WriteLine($"  [{hessian[1,0]:F4}, {hessian[1,1]:F4}]");
// Expected: [[2, 1], [1, 4]]
```

---
## Error Analysis

### Truncation Error

For central differences with step size $h$:

```math
\text{Error} \approx \frac{h^2}{6}|f'''(\xi)|
```

### Roundoff Error

```math
\text{Error} \approx \frac{2\epsilon_m |f(x)|}{h}
```

### Total Error

The total error is minimized at:

```math
h_{opt} = \left(\frac{3\epsilon_m |f(x)|}{|f'''(x)|}\right)^{1/3}
```

```cs
// Demonstrate error vs step size
Func<double, double> f = x => Math.Sin(x);
double x = 1.0;
double exact = Math.Cos(x);

Console.WriteLine("Step Size     Error");
Console.WriteLine("---------     -----");

for (int i = 1; i <= 16; i++)
{
    double h = Math.Pow(10, -i);
    double approx = NumericalDerivative.CentralDifference(f, x, h);
    double error = Math.Abs(approx - exact);
    Console.WriteLine($"1e-{i,-2}        {error:E2}");
}
```

Typical output shows error decreasing until ~1e-8, then increasing due to roundoff.

---

## Applications

### Sensitivity Analysis

```cs
// Sensitivity of output to input parameters
Func<double[], double> model = parameters =>
{
    double a = parameters[0];
    double b = parameters[1];
    return a * Math.Exp(-b * 2.0);  // Model output at t=2
};

double[] nominalParams = { 1.0, 0.5 };
double[] gradient = NumericalDerivative.Gradient(model, nominalParams);

double output = model(nominalParams);
Console.WriteLine($"Output: {output:F6}");
Console.WriteLine($"Sensitivity to a: {gradient[0]:F6}");
Console.WriteLine($"Sensitivity to b: {gradient[1]:F6}");
```

### Newton-Raphson with Numerical Derivative

```cs
Func<double, double> f = x => x * x - 2;  // Find √2

double x = 1.0;
for (int i = 0; i < 10; i++)
{
    double fx = f(x);
    double fpx = NumericalDerivative.Derivative(f, x);
    x = x - fx / fpx;
    Console.WriteLine($"Iteration {i+1}: x = {x:F12}");
}
```

### Optimization Gradient Check

Verify analytical gradients against numerical:

```cs
Func<double[], double> f = x => x[0]*x[0] + x[1]*x[1];
Func<double[], double[]> analyticalGradient = x => new double[] { 2*x[0], 2*x[1] };

double[] point = { 3.0, 4.0 };

double[] numerical = NumericalDerivative.Gradient(f, point);
double[] analytical = analyticalGradient(point);

Console.WriteLine("Gradient Check:");
for (int i = 0; i < point.Length; i++)
{
    double relError = Math.Abs(numerical[i] - analytical[i]) / Math.Abs(analytical[i]);
    Console.WriteLine($"  Component {i}: Numerical={numerical[i]:F6}, Analytical={analytical[i]:F6}, RelError={relError:E2}");
}
```

---

## Best Practices

1. **Use central differences** as the default method
2. **Use automatic step size** unless you have specific requirements
3. **Use Ridder's method** when high accuracy is needed and the function is smooth
4. **Check derivatives numerically** when implementing analytical gradients
5. **Be cautious near discontinuities** - finite differences assume smoothness
6. **Scale variables** to similar magnitudes for multivariate problems

---

## References

<a id="ref1">[1]</a> Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

<a id="ref2">[2]</a> Ridders, C. J. F. (1982). Accurate computation of F'(x) and F'(x)F''(x). *Advances in Engineering Software*, 4(2), 75-76.

<a id="ref3">[3]</a> Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.). Springer.
