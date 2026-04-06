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

Brent's method is the recommended general-purpose root-finding algorithm in the ***Numerics*** library. It combines the guaranteed convergence of bisection with the speed of the secant method and inverse quadratic interpolation, automatically switching between these strategies to achieve robust, fast convergence without requiring derivatives [[2]](#2) [[3]](#3). For the vast majority of root-finding problems, Brent's method should be the first method you reach for.

### Mathematical Description

Brent's method maintains a bracketing interval $[a, b]$ such that $f(a)$ and $f(b)$ have opposite signs, guaranteeing that a root exists within the interval by the Intermediate Value Theorem. At each iteration, it selects one of three strategies to propose the next approximation:

**Bisection.** The simplest strategy: take the midpoint of the current bracket.

```math
x_{\text{bisect}} = \frac{a + b}{2}
```

Bisection reduces the bracket by exactly half at each step, giving linear convergence with error bound $|e_n| \leq (b - a) / 2^n$. It is slow but absolutely reliable.

**Secant method.** When only two distinct points are available (i.e., the previous contrapoint $a$ equals the older contrapoint $c$), the algorithm uses the secant formula, which fits a line through the two most recent function values and finds where it crosses zero:

```math
x_{\text{secant}} = b - f(b) \cdot \frac{b - a}{f(b) - f(a)}
```

The secant method converges superlinearly with order approximately $\varphi \approx 1.618$ (the golden ratio) for well-behaved functions.

**Inverse quadratic interpolation (IQI).** When three distinct points $a$, $b$, $c$ with distinct function values are available, the algorithm fits an inverse quadratic (a parabola through the three points with $x$ as a function of $y$) and evaluates it at $y = 0$:

```math
x_{\text{IQI}} = \frac{f_b f_c \cdot a}{(f_a - f_b)(f_a - f_c)} + \frac{f_a f_c \cdot b}{(f_b - f_a)(f_b - f_c)} + \frac{f_a f_b \cdot c}{(f_c - f_a)(f_c - f_b)}
```

where $f_a = f(a)$, $f_b = f(b)$, $f_c = f(c)$. IQI can converge even faster than the secant method when the function is smooth and the iterates are close to the root.

### How It Works: The Decision Logic

The power of Brent's method lies in how it decides which strategy to use at each step. The algorithm follows a specific decision procedure that ensures it never loses the safety of bisection while taking faster steps whenever possible. Here is the logic as implemented in the ***Numerics*** library:

1. **Maintain the bracket.** The algorithm tracks three points: $b$ (the current best estimate, where $|f(b)|$ is smallest), $a$ (the previous iterate), and $c$ (the contrapoint such that $f(b)$ and $f(c)$ have opposite signs). If $f(b)$ and $f(c)$ have the same sign, the contrapoint is reset to $a$.

2. **Compute the midpoint and tolerance.** At each step, compute:

```math
x_m = \frac{c - b}{2}, \qquad \text{tol}_1 = 2 \varepsilon_{\text{mach}} |b| + \frac{\text{tol}}{2}
```

where $\varepsilon_{\text{mach}}$ is machine epsilon and $\text{tol}$ is the user-specified tolerance. If $|x_m| \leq \text{tol}_1$ or $f(b) = 0$, the root has been found.

3. **Try a fast step.** If the previous step $e$ was large enough ($|e| \geq \text{tol}_1$) and $b$ is improving ($|f(a)| > |f(b)|$), then attempt an open method:
   - If $a = c$ (only two distinct points available), use the **secant method**.
   - If $a \neq c$ (three distinct points available), use **inverse quadratic interpolation**.

4. **Accept or reject the fast step.** The proposed step $d = p/q$ is accepted only if it satisfies two conditions:
   - The step must be smaller than three-quarters of the distance to the midpoint: $2|p| < 3 |x_m \cdot q| - |\text{tol}_1 \cdot q|$
   - The step must be smaller than half the previous step: $2|p| < |e \cdot q|$

   These conditions ensure the method is making adequate progress toward the root. If either condition fails, the algorithm falls back to bisection.

5. **Update.** Apply the accepted step (or bisection fallback) to produce the new iterate $b$, evaluate $f(b)$, and repeat.

This design means that Brent's method takes fast steps when they are safe and productive, but always has bisection as a backstop, so it never diverges.

### Convergence Properties

Brent's method offers a rare combination of guaranteed convergence and fast practical performance:

- **Guaranteed convergence.** Like bisection, the method always maintains a valid bracket around the root. It will converge to a root for any continuous function where the initial interval satisfies $f(a) \cdot f(b) < 0$, regardless of how ill-conditioned the function is.

- **Superlinear convergence in practice.** For smooth, well-behaved functions, the algorithm typically achieves superlinear convergence by spending most iterations in secant or IQI mode. In the best case, IQI converges with order approximately 1.839.

- **Worst case is bisection.** When the function is poorly behaved (discontinuous derivatives, sharp curvature changes, or near-flat regions), the method gracefully degrades to bisection's linear convergence rate of $O(1/2^n)$. This is a floor, not a failure.

- **No pathological failure modes.** Unlike Newton-Raphson, which can cycle, diverge, or overshoot for bad initial guesses, and unlike the secant method, which can leave the bracket entirely, Brent's method has no failure modes beyond exceeding the maximum iteration count.

The convergence tolerance in the ***Numerics*** implementation uses a combined criterion: the root is accepted when the bracket half-width $|x_m|$ is within $\text{tol}_1 = 2\varepsilon_{\text{mach}}|b| + \text{tol}/2$, or when $f(b) = 0$ exactly. This accounts for both absolute and relative precision near the root.

### Why Brent's Method Is the Recommended Default

For solving $f(x) = 0$ in a single variable, Brent's method should be the default choice because:

- **No derivatives needed.** Unlike Newton-Raphson, Brent's method requires only function evaluations. There is no need to derive, implement, or numerically approximate $f'(x)$, which is a significant practical advantage.

- **Guaranteed to converge.** Given a valid bracketing interval where $f(a)$ and $f(b)$ have opposite signs, the method will always find a root. Newton-Raphson and the secant method have no such guarantee.

- **Typically faster than bisection.** While bisection requires approximately $\log_2((b-a)/\varepsilon)$ iterations, Brent's method usually converges in far fewer iterations by exploiting the smoothness of the function. For the test function $f(x) = x^3 - 2x - 5$ on $[2, 3]$, bisection needs 27 iterations while Brent needs only 5.

- **No pathological failure cases.** Newton-Raphson can diverge, oscillate, or cycle when the initial guess is poor, the derivative is near zero, or the function has inflection points. Brent's method avoids all of these failure modes.

- **Handles both well-behaved and ill-conditioned functions.** The adaptive switching between IQI, secant, and bisection means the algorithm performs well across a wide range of function behaviors without requiring the user to diagnose the function's properties in advance.

### When NOT to Use Brent's Method

While Brent's method is the best general-purpose choice, there are situations where a different method is more appropriate:

- **Derivatives are available and the function is smooth.** If you already have an analytical derivative $f'(x)$ and the function is smooth with a good initial guess, Newton-Raphson will converge faster (quadratically vs. superlinearly). For functions where the derivative is cheap to evaluate, the `NewtonRaphson.RobustSolve` method provides Newton's speed with bisection's safety net.

- **Systems of equations or multiple dimensions.** Brent's method is strictly a univariate solver. For systems of nonlinear equations $\mathbf{F}(\mathbf{x}) = \mathbf{0}$, use multidimensional optimization methods such as Nelder-Mead or Newton-based nonlinear solvers.

- **The root is not bracketed.** Brent's method requires an initial interval $[a, b]$ where $f(a)$ and $f(b)$ have opposite signs. If you cannot identify such a bracket, consider using the `Brent.Bracket` helper method (see below), the secant method, or Newton-Raphson which do not require bracketing.

- **The function does not change sign at the root.** For roots of even multiplicity (e.g., $f(x) = x^2$ has a root at $x = 0$ but does not change sign), bracketing methods cannot be used. In these cases, reformulate the problem or use Newton-Raphson on $g(x) = f(x)/f'(x)$.

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

The `Solve` method accepts the following parameters:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `f` | `Func<double, double>` | (required) | The function to solve |
| `lowerBound` | `double` | (required) | Lower bound $a$ of the bracketing interval |
| `upperBound` | `double` | (required) | Upper bound $b$ of the bracketing interval |
| `tolerance` | `double` | `1e-8` | Desired tolerance for the root |
| `maxIterations` | `int` | `1000` | Maximum number of iterations |
| `reportFailure` | `bool` | `true` | If true, throws an exception on failure; if false, returns the last iterate |

### Automatic Bracket Expansion

When you are unsure whether your initial interval brackets the root, the `Brent.Bracket` helper method can expand the interval outward until a sign change is found. It repeatedly expands the interval by a factor of 1.6, testing the function value at each new endpoint:

```cs
using Numerics.Mathematics.RootFinding;

// We think the root is near x = 5, but are not sure of the bracket
Func<double, double> f = x => x * x - 25;

double a = 4.0, b = 6.0;
double fa, fb;

bool found = Brent.Bracket(f, ref a, ref b, out fa, out fb, maxIterations: 10);

if (found)
{
    double root = Brent.Solve(f, a, b);
    Console.WriteLine($"Root: {root:F10}");  // 5.0000000000
}
else
{
    Console.WriteLine("Could not bracket the root.");
}
```

### Example: Finding Where Two Functions Intersect

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

### Example: Hydrologic Design with Manning's Equation

A common civil engineering problem: determine the water depth in a trapezoidal channel that produces a given discharge. Manning's equation relates flow rate to channel geometry:

```math
Q = \frac{1}{n} A R^{2/3} S^{1/2}
```

where $Q$ is discharge, $n$ is Manning's roughness coefficient, $A$ is cross-sectional area, $R = A/P$ is hydraulic radius, $P$ is wetted perimeter, and $S$ is channel slope. For a trapezoidal channel with bottom width $b_w$ and side slope $z$:

```math
A = y(b_w + z \cdot y), \qquad P = b_w + 2y\sqrt{1 + z^2}
```

Given a target discharge, solve for the normal depth $y$:

```cs
using Numerics.Mathematics.RootFinding;

// Channel parameters
double n = 0.030;    // Manning's roughness (natural channel)
double bw = 10.0;    // Bottom width (m)
double z = 2.0;      // Side slope (horizontal:vertical)
double S = 0.001;    // Channel slope
double Qtarget = 50; // Target discharge (m³/s)

// Manning's equation residual: f(y) = Q(y) - Qtarget
Func<double, double> f = y =>
{
    double A = y * (bw + z * y);
    double P = bw + 2 * y * Math.Sqrt(1 + z * z);
    double R = A / P;
    double Q = (1.0 / n) * A * Math.Pow(R, 2.0 / 3.0) * Math.Sqrt(S);
    return Q - Qtarget;
};

// Solve for normal depth
double depth = Brent.Solve(f,
    lowerBound: 0.01,    // Minimum physical depth
    upperBound: 20.0,    // Maximum reasonable depth
    tolerance: 1e-6);

Console.WriteLine($"Normal depth: {depth:F4} m");
Console.WriteLine($"Verification: Q = {f(depth) + Qtarget:F4} m³/s");
```

### Example: Controlling Failure Behavior

In production code, you may want to handle convergence failure gracefully rather than allowing exceptions:

```cs
using Numerics.Mathematics.RootFinding;

Func<double, double> f = x => Math.Cos(x) - x;

// Suppress exception on failure — returns last iterate instead
double root = Brent.Solve(f, 0, 1,
    tolerance: 1e-15,       // Very tight tolerance
    maxIterations: 5,       // Very few iterations allowed
    reportFailure: false);  // Don't throw on failure

// Always verify the result
double residual = Math.Abs(f(root));
if (residual < 1e-10)
    Console.WriteLine($"Root found: {root:F12}");
else
    Console.WriteLine($"Warning: residual {residual:E3} exceeds acceptable threshold");
```

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

### Method Comparison

The following table provides a comprehensive comparison of all root-finding methods available in the ***Numerics*** library:

| Feature | Bisection | Brent | Secant | Newton-Raphson |
|---------|-----------|-------|--------|----------------|
| **Convergence rate** | Linear $O(1/2^n)$ | Superlinear (adaptive) | ~1.618 (golden ratio) | Quadratic |
| **Requires derivative** | No | No | No | Yes |
| **Requires bracket** | Yes | Yes | No | No |
| **Guaranteed convergence** | Yes | Yes | No | No |
| **Function evals per step** | 1 | 1 | 1 | 2 ($f$ and $f'$) |
| **Typical iterations** | 30-50 | 5-10 | 5-10 | 4-6 |
| **Handles discontinuities** | Yes | Yes | Poorly | Poorly |
| **Risk of divergence** | None | None | Moderate | High |
| **Best for** | Proof of concept, robustness | General purpose (recommended) | When bracket unavailable | When $f'(x)$ is cheap and smooth |

### Scenario-Based Recommendations

| Scenario | Recommended Method | Notes |
|----------|-------------------|-------|
| General purpose | **Brent** | Best balance of speed and robustness |
| Have analytical derivative | Newton-Raphson (robust version) | Fastest convergence with safety net |
| No derivative available | Brent | More robust than Secant |
| Difficult or ill-conditioned function | Bisection or Brent | Guaranteed convergence |
| Need absolute certainty | Brent | Guaranteed convergence with fast speed |
| Very smooth function with good guess | Newton-Raphson | Quadratic convergence |
| Poor or uncertain initial guess | Brent or Robust Newton-Raphson | Fall back to bisection when needed |
| Cannot bracket the root | Secant or Newton-Raphson | Consider `Brent.Bracket` to find bracket first |

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

<a id="3">[3]</a> Sprott, J. C. (1991). *Numerical Recipes, Routines and Examples in Basic*. Cambridge University Press.

<a id="4">[4]</a> Süli, E. & Mayers, D. (2003). *An Introduction to Numerical Analysis*. Cambridge University Press.

---

[← Previous: Optimization](optimization.md) | [Back to Index](../index.md) | [Next: Linear Algebra →](linear-algebra.md)
