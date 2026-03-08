# Root Finding

[← Previous: Optimization](optimization.md) | [Back to Index](../index.md) | [Next: Linear Algebra →](linear-algebra.md)

Root finding is the process of determining the values of $x$ for which a function $f(x) = 0$. These are also called zeros, roots, or solutions of the equation. Root finding is fundamental to many numerical methods and applications, including solving nonlinear equations, finding equilibrium points, and numerical integration of differential equations.

## Overview

The ***Numerics*** library provides several robust algorithms for finding roots of univariate functions:

| Method | Requires | Convergence | Best For |
|--------|----------|-------------|----------|
| **Bisection** | Bracketing interval | Linear (slow but sure) | Robust, guaranteed convergence |
| **Secant** | Two initial points | Superlinear (~1.618) | When derivative unavailable |
| **Newton-Raphson** | Derivative | Quadratic (very fast) | Smooth functions, good initial guess |
| **Brent** | Bracketing interval | Superlinear | General purpose, best overall |

## Problem Formulation

Given a function $f: \mathbb{R} \rightarrow \mathbb{R}$, find $x^*$ such that:

```math
f(x^*) = 0
```

### Bracketing

Some methods require a **bracketing interval** $[a, b]$ where $f(a)$ and $f(b)$ have opposite signs. By the Intermediate Value Theorem, if $f$ is continuous, there must be at least one root in the interval.

## Bisection Method

The bisection method is the simplest root-finding algorithm. It repeatedly bisects an interval and selects the subinterval in which the root must lie [[1]](#1).

### Algorithm

1. Start with interval $[a, b]$ where $f(a) \cdot f(b) < 0$
2. Compute midpoint $c = (a + b) / 2$
3. If $f(c)$ is close enough to zero, return $c$
4. Otherwise, replace either $a$ or $b$ with $c$ based on the sign of $f(c)$
5. Repeat until convergence

### Convergence Rate

Bisection has **linear convergence** with a constant factor of $1/2$. After $n$ iterations, the error is bounded by:

```math
|e_n| \leq \frac{b - a}{2^n}
```

where $b - a$ is the initial interval width. This means each iteration gains approximately one binary digit of accuracy. To achieve a tolerance $\varepsilon$, the number of iterations required is:

```math
n \geq \frac{\log(b - a) - \log(\varepsilon)}{\log 2}
```

For example, starting with $[0, 3]$ and tolerance $10^{-10}$, bisection needs at most $\lceil \log_2(3 \times 10^{10}) \rceil = 35$ iterations. This predictability is one of bisection's strengths — you know exactly how many iterations you need before you start.

### Usage

```cs
using Numerics.Mathematics.RootFinding;

// Find root of f(x) = x² - 4 (roots at x = ±2)
Func<double, double> f = x => x * x - 4;

// Bisection requires a bracketing interval
double root = Bisection.Solve(f, 
    firstGuess: 1.5,      // Optional hint for initial interval
    lowerBound: 0,        // f(0) = -4 < 0
    upperBound: 3,        // f(3) = 5 > 0
    tolerance: 1e-8,
    maxIterations: 1000,
    reportFailure: true);

Console.WriteLine($"Root: {root:F10}");  // 2.0000000000
Console.WriteLine($"Verification: f({root}) = {f(root):E3}");
```

### Advantages and Disadvantages

**Advantages:**
- Always converges if initial interval brackets a root
- Very robust - works even for discontinuous functions
- Simple to implement and understand
- Guaranteed to find a root

**Disadvantages:**
- Slow convergence (linear, each iteration reduces error by half)
- Requires bracketing interval
- Cannot find roots where function doesn't change sign (e.g., $f(x) = x^2$)

**When to use:** When robustness is paramount, or for poorly behaved functions.

## Brent's Method

Brent's method combines the robustness of bisection with the speed of secant and inverse quadratic interpolation [[2]](#2). It's generally the best general-purpose root finder.

### Algorithm

The method maintains a bracketing interval and uses:
- **Inverse quadratic interpolation** when three points are available
- **Secant method** when two points are available
- **Bisection** as a fallback to guarantee convergence

### Usage

```cs
using Numerics.Mathematics.RootFinding;

// Find root of f(x) = cos(x) - x (root around x ≈ 0.739)
Func<double, double> f = x => Math.Cos(x) - x;

double root = Brent.Solve(f,
    lowerBound: 0,        // f(0) = 1 > 0
    upperBound: 1,        // f(1) = -0.46 < 0
    tolerance: 1e-10,
    maxIterations: 1000);

Console.WriteLine($"Root: {root:F12}");  // 0.739085133215
Console.WriteLine($"Verification: f({root}) = {f(root):E3}");
```

### Example: Finding where two functions intersect

To find where $f(x) = g(x)$, solve $h(x) = f(x) - g(x) = 0$:

```cs
Func<double, double> f = x => Math.Exp(-x);
Func<double, double> g = x => x * x;

// Find intersection: e^(-x) = x²
Func<double, double> h = x => f(x) - g(x);

double intersection = Brent.Solve(h, 0, 1);
Console.WriteLine($"Intersection at x = {intersection:F6}");
Console.WriteLine($"At this point: f({intersection}) = {f(intersection):F6}");
Console.WriteLine($"At this point: g({intersection}) = {g(intersection):F6}");
```

### Advantages and Disadvantages

**Advantages:**
- Very fast convergence (superlinear)
- Guaranteed to converge with bracketing interval
- Automatically switches between methods for optimal performance
- Widely considered the best general-purpose method

**Disadvantages:**
- Requires bracketing interval
- More complex than simpler methods

**When to use:** Default choice for most root-finding problems.

## Secant Method

The secant method is similar to Newton-Raphson but approximates the derivative using finite differences, eliminating the need for an analytical derivative [[1]](#1).

### Algorithm

The secant formula is:

```math
x_{n+1} = x_n - f(x_n) \frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}
```

This approximates the derivative as:

```math
f'(x_n) \approx \frac{f(x_n) - f(x_{n-1})}{x_n - x_{n-1}}
```

### Convergence Rate

The secant method has **superlinear convergence** of order $\varphi = (1 + \sqrt{5})/2 \approx 1.618$ (the golden ratio). Near a simple root $x^*$, the error satisfies:

```math
|e_{n+1}| \approx \left|\frac{f''(x^*)}{2f'(x^*)}\right|^{\varphi - 1} |e_n|^{\varphi}
```

This is faster than bisection's linear convergence but slower than Newton's quadratic convergence. However, since the secant method requires only **one function evaluation per iteration** (compared to Newton's two — one for $f$ and one for $f'$), it can be more efficient overall. In terms of function evaluations, the secant method's efficiency index is $\varphi^{1/1} \approx 1.618$, while Newton's is $2^{1/2} \approx 1.414$.

### Usage

```cs
using Numerics.Mathematics.RootFinding;

// Find root of f(x) = x³ - 2x - 5
Func<double, double> f = x => x * x * x - 2 * x - 5;

double root = Secant.Solve(f,
    lowerBound: 2,        // First point
    upperBound: 3,        // Second point
    tolerance: 1e-8,
    maxIterations: 1000);

Console.WriteLine($"Root: {root:F10}");  // 2.0945514815
Console.WriteLine($"Verification: f({root}) = {f(root):E3}");
```

### Advantages and Disadvantages

**Advantages:**
- Faster than bisection (superlinear convergence)
- Doesn't require derivative
- Only needs one function evaluation per iteration (vs. two for Newton-Raphson)

**Disadvantages:**
- Not guaranteed to converge
- May diverge for poor initial guesses
- Slower than Newton-Raphson when derivative is available

**When to use:** When derivative is unavailable or expensive to compute, and you have reasonable initial estimates.

## Newton-Raphson Method

Newton-Raphson uses the function's derivative to iteratively approach the root [[1]](#1). It has quadratic convergence, making it very fast when it works.

### Algorithm

The Newton-Raphson formula is:

```math
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
```

Geometrically, this finds where the tangent line at $(x_n, f(x_n))$ crosses the x-axis.

### Convergence Rate

Newton-Raphson has **quadratic convergence** near a simple root. If $e_n = x_n - x^*$ is the error at step $n$, then:

```math
e_{n+1} = -\frac{f''(x^*)}{2f'(x^*)} \cdot e_n^2 + O(e_n^3)
```

This means the number of correct digits roughly **doubles** each iteration. Starting from $|e_0| = 0.1$, after 4 iterations the error is approximately $10^{-16}$ — near machine precision. This is why Newton's method typically converges in 4–6 iterations.

The constant $C = |f''(x^*)|/(2|f'(x^*)|)$ determines how fast convergence kicks in. When $|f''|$ is large relative to $|f'|$ at the root (i.e., the function curves sharply while crossing zero slowly), the method needs a closer starting point for quadratic convergence to dominate.

For **repeated roots** where $f(x^*) = f'(x^*) = 0$, convergence degrades to linear. In such cases, modified Newton's method using $g(x) = f(x)/f'(x)$ restores quadratic convergence.

### Basin of Attraction

Newton's method is not globally convergent — its behavior depends critically on the initial guess. The **basin of attraction** for a root $x^*$ is the set of starting points $x_0$ from which Newton's method converges to $x^*$.

For polynomials with multiple roots, basins of attraction can have fractal boundaries. Even simple functions can produce surprising behavior:

- For $f(x) = x^3 - 1$ (with roots at $1$, $e^{2\pi i/3}$, $e^{4\pi i/3}$ in the complex plane), the boundary between basins of attraction is a fractal — the Newton fractal.

- For $f(x) = x^3 - 2x + 2$, starting at $x_0 = 0$ produces the cycle $0 \to 1 \to 1 \to \ldots$ that never converges.

- For $f(x) = |x|^a$ where $0 < a < 1/2$, Newton's method diverges from every starting point except the root itself.

A practical rule of thumb: Newton's method converges when $|f(x_0) \cdot f''(x_0)| < |f'(x_0)|^2$ (the Newton-Kantorovich condition). When unsure about the initial guess, use the `RobustSolve` method which falls back to bisection.

### Usage

```cs
using Numerics.Mathematics.RootFinding;

// Find root of f(x) = x² - 2 (square root of 2)
Func<double, double> f = x => x * x - 2;
Func<double, double> df = x => 2 * x;  // Derivative: f'(x) = 2x

double root = NewtonRaphson.Solve(f, df,
    firstGuess: 1.0,      // Initial guess
    tolerance: 1e-12,
    maxIterations: 100);

Console.WriteLine($"√2 = {root:F12}");  // 1.414213562373
Console.WriteLine($"Verification: {root}² = {root * root:F12}");
```

**Note:** If you don't have an analytical derivative, the ***Numerics*** library's `NumericalDerivative` class can compute it numerically, though this is less efficient:

```cs
using Numerics.Mathematics;
using Numerics.Mathematics.RootFinding;

Func<double, double> f = x => x * x - 2;
Func<double, double> df = x => NumericalDerivative.Derivative(f, x);

double root = NewtonRaphson.Solve(f, df, 1.0);
```

### Robust Newton-Raphson

The library also provides a robust version that combines Newton-Raphson with bisection to guarantee convergence:

```cs
// Robust version requires bracketing interval
double root = NewtonRaphson.RobustSolve(f, df,
    firstGuess: 1.0,
    lowerBound: 0.5,
    upperBound: 2.0,
    tolerance: 1e-12);
```

This method uses Newton-Raphson when it's making good progress, but falls back to bisection if the iteration goes outside the brackets or converges too slowly.

### Advantages and Disadvantages

**Advantages:**
- Very fast convergence (quadratic)
- Few iterations needed (typically 4-6)
- Can be extended to systems of equations

**Disadvantages:**
- Requires derivative (analytical or numerical)
- Not guaranteed to converge
- Can fail for poor initial guesses
- Can diverge or oscillate

**When to use:** When you have a good initial guess and can provide the derivative, or use the robust version for guaranteed convergence.

## Example: Solving for Implied Volatility

A practical example from financial mathematics - solving the Black-Scholes equation for implied volatility:

```cs
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

// Black-Scholes call option price
double BlackScholesCall(double S, double K, double T, double r, double sigma)
{
    double d1 = (Math.Log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * Math.Sqrt(T));
    double d2 = d1 - sigma * Math.Sqrt(T);

    // Standard normal CDF using the error function
    double N(double x) => 0.5 * (1.0 + Erf.Function(x / Math.Sqrt(2)));

    return S * N(d1) - K * Math.Exp(-r * T) * N(d2);
}

// Given market price, solve for implied volatility
double S = 100;   // Stock price
double K = 100;   // Strike price
double T = 1;     // Time to maturity (years)
double r = 0.05;  // Risk-free rate
double marketPrice = 12.34;  // Observed option price

// Define function: difference between model and market price
Func<double, double> f = sigma => BlackScholesCall(S, K, T, r, sigma) - marketPrice;

// Solve for implied volatility using Brent
double impliedVol = Brent.Solve(f, 0.01, 2.0);

Console.WriteLine($"Implied volatility: {impliedVol:P2}");
Console.WriteLine($"Verification: BS price = {BlackScholesCall(S, K, T, r, impliedVol):F2}");
```

## Finding Multiple Roots

To find multiple roots, solve in different intervals:

```cs
using Numerics.Mathematics.RootFinding;

// Function with multiple roots: f(x) = sin(x)
// Roots at x = 0, ±π, ±2π, ...
Func<double, double> f = x => Math.Sin(x);

// Find roots in different intervals
var roots = new List<double>();

// Find positive roots
for (int i = 0; i < 3; i++)
{
    double a = i * Math.PI + 0.1;
    double b = (i + 1) * Math.PI - 0.1;
    
    if (f(a) * f(b) < 0) // Bracket contains a root
    {
        double root = Brent.Solve(f, a, b);
        roots.Add(root);
    }
}

Console.WriteLine("Roots of sin(x):");
foreach (var root in roots)
{
    Console.WriteLine($"  x = {root:F10}  (≈ {root / Math.PI:F2}π)");
}
```

## Example: Critical Points via Root Finding

Find critical points of a function by solving $f'(x) = 0$:

```cs
using Numerics.Mathematics.RootFinding;

// Find critical points of f(x) = x³ - 3x² + 2
Func<double, double> f = x => x * x * x - 3 * x * x + 2;
Func<double, double> df = x => 3 * x * x - 6 * x;  // f'(x)

// f'(x) = 3x² - 6x = 3x(x - 2)
// Critical points at x = 0 and x = 2

double cp1 = Brent.Solve(df, -1, 1);    // Find critical point near 0
double cp2 = Brent.Solve(df, 1, 3);     // Find critical point near 2

Console.WriteLine($"Critical point 1: x = {cp1:F6}, f(x) = {f(cp1):F6}");
Console.WriteLine($"Critical point 2: x = {cp2:F6}, f(x) = {f(cp2):F6}");

// Determine if max or min using second derivative
Func<double, double> d2f = x => 6 * x - 6;  // f''(x)

Console.WriteLine($"At x={cp1}: f''(x)={d2f(cp1):F1} → " + 
                  (d2f(cp1) > 0 ? "Local minimum" : "Local maximum"));
Console.WriteLine($"At x={cp2}: f''(x)={d2f(cp2):F1} → " + 
                  (d2f(cp2) > 0 ? "Local minimum" : "Local maximum"));
```

## Choosing a Root Finding Method

| Scenario | Recommended Method | Notes |
|----------|-------------------|-------|
| General purpose | Brent | Best balance of speed and robustness |
| Have derivative | Newton-Raphson (robust version) | Fastest convergence |
| No derivative | Secant or Brent | Brent more robust, Secant faster |
| Difficult function | Bisection | Guaranteed convergence, but slow |
| Need absolute certainty | Bisection or Brent | Both guarantee convergence |
| Very smooth function | Newton-Raphson | Quadratic convergence |
| Poor initial guess | Brent or Robust Newton-Raphson | Fall back to bisection when needed |

## Convergence Criteria

All root-finding methods in ***Numerics*** use combined criteria:

1. **Function value tolerance**: $|f(x)| < \text{tolerance}$
2. **Parameter tolerance**: $|x_{n+1} - x_n| < \text{tolerance}$

Both must be satisfied for convergence. The default tolerance is $10^{-8}$.

## Best Practices

1. **Bracket First**: For robust methods (Bisection, Brent), ensure your initial interval brackets the root by checking that $f(a) \cdot f(b) < 0$.

2. **Plot the Function**: Before root finding, plot the function to understand its behavior and identify approximate root locations.

3. **Check Convergence**: Always verify the solution:
   ```cs
   double root = Brent.Solve(f, a, b);
   Console.WriteLine($"f({root}) = {f(root):E3}");  // Should be near zero
   ```

4. **Handle No Root Cases**: Check if a root exists in your interval before calling the solver:
   ```cs
   if (f(a) * f(b) >= 0)
   {
       Console.WriteLine("Function doesn't change sign - root may not exist");
   }
   ```

5. **Multiple Roots**: To find all roots, divide the search space and solve in each bracketing interval.

6. **Scale Appropriately**: Normalize your function if it has extreme magnitudes to avoid numerical issues.

7. **Set Reasonable Tolerances**: Tighter tolerances require more iterations:
   ```cs
   double root = Brent.Solve(f, a, b, tolerance: 1e-12);
   ```

8. **Use Robust Versions for Production**: When failure is not an option, use Brent or Robust Newton-Raphson.

## Common Pitfalls

1. **No Root in Interval**: Bisection and Brent will fail if $f(a)$ and $f(b)$ have the same sign.

2. **Multiple Roots**: These methods find one root at a time. Multiple roots in the same interval may cause issues.

3. **Flat Regions**: All methods struggle where $f'(x) \approx 0$ (Newton-Raphson especially).

4. **Discontinuities**: Newton-Raphson and Secant may fail near discontinuities. Brent and Bisection are more robust.

5. **Poor Initial Guess**: Newton-Raphson is sensitive to the starting point. Use Robust Newton-Raphson or Brent if initial guess quality is uncertain.

## Error Handling

```cs
using Numerics.Mathematics.RootFinding;

try
{
    double root = Brent.Solve(f, a, b, 
        tolerance: 1e-10,
        maxIterations: 1000,
        reportFailure: true);  // Throw exception on failure
    
    Console.WriteLine($"Root found: {root}");
}
catch (Exception ex)
{
    Console.WriteLine($"Root finding failed: {ex.Message}");
}

// Or suppress exceptions
double root2 = Brent.Solve(f, a, b, reportFailure: false);
if (double.IsNaN(root2))
{
    Console.WriteLine("Failed to converge");
}
```

## Performance Comparison

For the function $f(x) = x^3 - 2x - 5$ with root at $x \approx 2.0946$:

| Method | Initial Values | Iterations | Function Evals | Time (relative) |
|--------|---------------|-----------|----------------|-----------------|
| Bisection | [2, 3] | 27 | 27 | 1.0× |
| Secant | 2, 3 | 6 | 6 | 0.22× |
| Newton-Raphson | 2 | 4 | 8 | 0.30× |
| Brent | [2, 3] | 5 | 8 | 0.30× |

Newton-Raphson and Brent are typically fastest, while Bisection is slowest but most reliable.

---

## References

<a id="1">[1]</a> Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

<a id="2">[2]</a> Brent, R. P. (1973). *Algorithms for Minimization Without Derivatives*. Prentice-Hall, Englewood Cliffs, NJ.

<a id="3">[3]</a> Süli, E. & Mayers, D. (2003). *An Introduction to Numerical Analysis*. Cambridge University Press.

---

[← Previous: Optimization](optimization.md) | [Back to Index](../index.md) | [Next: Linear Algebra →](linear-algebra.md)
