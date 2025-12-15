# Root Finding

The ***Numerics*** library provides algorithms for finding roots of univariate functions, from simple bracketing methods to Newton-type methods that use derivatives.

## Overview

| Method | Convergence | Requires | Best For |
|--------|-------------|----------|----------|
| Bisection | Linear | Bracket | Guaranteed convergence |
| Regula Falsi | Superlinear | Bracket | Simple problems |
| Brent's Method | Superlinear | Bracket | General purpose |
| Newton-Raphson | Quadratic | Derivative | Fast convergence |
| Secant | Superlinear | Two points | No derivative available |

---

## Bracketing Methods

These methods require an initial interval $[a, b]$ where $f(a)$ and $f(b)$ have opposite signs, guaranteeing a root exists by the Intermediate Value Theorem.

### Bisection Method

The simplest and most robust method [[1]](#ref1):

```cs
using Numerics.Mathematics.RootFinding;

Func<double, double> f = x => x * x - 2;  // Find √2

var bisect = new Bisection(f, 1, 2);
bisect.Tolerance = 1e-12;
bisect.MaxIterations = 100;

double root = bisect.Solve();

Console.WriteLine($"Root: {root:F12}");
Console.WriteLine($"f(root): {f(root):E2}");
Console.WriteLine($"Iterations: {bisect.Iterations}");
```

**Algorithm**: Repeatedly halves the interval, selecting the half containing the root.

**Convergence**: Linear—gains one bit of accuracy per iteration. After $n$ iterations, error ≤ $(b-a)/2^n$.

**Advantages**: Always converges, simple, robust.

**Disadvantages**: Slow compared to other methods.

### Regula Falsi (False Position)

Uses linear interpolation instead of midpoint [[1]](#ref1):

```cs
var regulaFalsi = new RegulaFalsi(f, 1, 2);
regulaFalsi.Tolerance = 1e-12;

double root = regulaFalsi.Solve();
Console.WriteLine($"Root: {root:F12}");
```

**Algorithm**: Draws a line between $(a, f(a))$ and $(b, f(b))$, uses the x-intercept as the new estimate.

**Convergence**: Superlinear for most functions, but can be slow if one endpoint remains fixed.

### Brent's Method

Combines bisection, secant, and inverse quadratic interpolation [[2]](#ref2):

```cs
var brent = new BrentSearch(f, 1, 2);
brent.Tolerance = 1e-12;

double root = brent.Solve();

Console.WriteLine($"Root: {root:F12}");
Console.WriteLine($"Iterations: {brent.Iterations}");
```

**Algorithm**: Uses inverse quadratic interpolation when safe, falls back to bisection otherwise.

**Convergence**: Superlinear (typically faster than secant), guaranteed to converge.

**Recommendation**: **Use Brent's method as the default bracketing method.**

---

## Open Methods

These methods do not require a bracket but may fail to converge if the initial guess is poor.

### Newton-Raphson Method

Uses the derivative for quadratic convergence [[1]](#ref1):

```math
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
```

```cs
Func<double, double> f = x => x * x - 2;
Func<double, double> df = x => 2 * x;  // Derivative

double x0 = 1.5;  // Initial guess

var newton = new NewtonRaphson(f, df, x0);
newton.Tolerance = 1e-12;
newton.MaxIterations = 50;

double root = newton.Solve();

Console.WriteLine($"Root: {root:F14}");
Console.WriteLine($"Iterations: {newton.Iterations}");
```

**Convergence**: Quadratic near the root—roughly doubles correct digits each iteration.

**Requirements**: Derivative must be provided, $f'(x) \neq 0$ near the root.

**Pitfalls**:
- May diverge if initial guess is far from root
- Fails at stationary points ($f'(x) = 0$)
- Cycles possible for some functions

### Newton-Raphson with Numerical Derivative

When the analytical derivative is not available:

```cs
Func<double, double> f = x => x * x - 2;

var newton = new NewtonRaphson(f, x0: 1.5);  // Uses numerical differentiation
double root = newton.Solve();
```

### Secant Method

Approximates the derivative using finite differences [[1]](#ref1):

```math
x_{n+1} = x_n - f(x_n)\frac{x_n - x_{n-1}}{f(x_n) - f(x_{n-1})}
```

```cs
Func<double, double> f = x => x * x - 2;

var secant = new SecantMethod(f, x0: 1, x1: 2);
secant.Tolerance = 1e-12;

double root = secant.Solve();

Console.WriteLine($"Root: {root:F14}");
Console.WriteLine($"Iterations: {secant.Iterations}");
```

**Convergence**: Superlinear (order ≈ 1.618, the golden ratio).

**Advantages**: No derivative needed, faster than bisection.

---

## Solving Equations

### Inverse CDF (Quantile) Computation

Find $x$ such that $F(x) = p$:

```cs
using Numerics.Distributions;

var dist = new Normal(100, 15);
double p = 0.95;

// Find x such that CDF(x) = 0.95
Func<double, double> f = x => dist.CDF(x) - p;

var brent = new BrentSearch(f, 50, 200);
double quantile = brent.Solve();

Console.WriteLine($"95th percentile: {quantile:F4}");
Console.WriteLine($"Verify: CDF({quantile:F4}) = {dist.CDF(quantile):F6}");
```

### Implicit Function Solutions

Find intersection points, equilibria, etc.:

```cs
// Find where sin(x) = x/2
Func<double, double> f = x => Math.Sin(x) - x / 2;

// Need to find brackets first
var brent = new BrentSearch(f, 0.1, 3);
double root = brent.Solve();

Console.WriteLine($"sin(x) = x/2 at x = {root:F10}");
```

### Internal Rate of Return (IRR)

Financial application—find rate where NPV = 0:

```cs
double[] cashflows = { -1000, 300, 400, 400, 300 };

Func<double, double> npv = r =>
{
    double sum = 0;
    for (int t = 0; t < cashflows.Length; t++)
        sum += cashflows[t] / Math.Pow(1 + r, t);
    return sum;
};

var brent = new BrentSearch(npv, 0, 1);
double irr = brent.Solve();

Console.WriteLine($"IRR: {irr * 100:F2}%");
```

---

## Finding Brackets

Before using bracketing methods, you need an interval containing the root.

### Sign Change Search

```cs
Func<double, double> f = x => Math.Exp(x) - 10;

double a = 0, b = 0;
bool found = false;

// Search outward from origin
for (double x = 0; x <= 10; x += 0.5)
{
    if (f(x) * f(x + 0.5) < 0)
    {
        a = x;
        b = x + 0.5;
        found = true;
        break;
    }
}

if (found)
{
    var brent = new BrentSearch(f, a, b);
    double root = brent.Solve();
    Console.WriteLine($"Root found: {root:F10}");
}
```

### Exponential Expansion

For unbounded search:

```cs
public static (double a, double b)? FindBracket(Func<double, double> f, double x0, double factor = 1.6, int maxIter = 50)
{
    double a = x0;
    double b = x0 + 1;
    
    for (int i = 0; i < maxIter; i++)
    {
        if (f(a) * f(b) < 0)
            return (Math.Min(a, b), Math.Max(a, b));
        
        if (Math.Abs(f(a)) < Math.Abs(f(b)))
            a += factor * (a - b);
        else
            b += factor * (b - a);
    }
    
    return null;  // No bracket found
}
```

---

## Multiple Roots

### Finding All Roots in an Interval

```cs
Func<double, double> f = x => Math.Sin(x);

double start = 0;
double end = 10;
double step = 0.5;
double tol = 1e-10;

var roots = new List<double>();

for (double x = start; x < end; x += step)
{
    if (f(x) * f(x + step) < 0)
    {
        var brent = new BrentSearch(f, x, x + step);
        brent.Tolerance = tol;
        double root = brent.Solve();
        
        // Avoid duplicates
        if (roots.Count == 0 || Math.Abs(root - roots.Last()) > step / 2)
            roots.Add(root);
    }
}

Console.WriteLine($"Found {roots.Count} roots:");
foreach (var root in roots)
    Console.WriteLine($"  x = {root:F10}, f(x) = {f(root):E2}");
```

### Polynomial Roots

For polynomials, use companion matrix eigenvalues:

```cs
// p(x) = x³ - 6x² + 11x - 6 = (x-1)(x-2)(x-3)
double[] coeffs = { -6, 11, -6, 1 };  // Constant term first

// Build companion matrix and find eigenvalues
// (Implementation depends on library support)
```

---

## Convergence Criteria

### Absolute Tolerance

Stop when $|f(x_n)| < \epsilon$:

```cs
solver.Tolerance = 1e-10;
```

### Relative Tolerance

Stop when $|x_{n+1} - x_n| < \epsilon |x_n|$:

```cs
solver.RelativeTolerance = 1e-10;
```

### Combined Criterion

```cs
// Stop when either criterion is met
bool converged = Math.Abs(f(x)) < absTol || 
                 Math.Abs(xNew - x) < relTol * Math.Abs(x);
```

---

## Troubleshooting

### Method Fails to Converge

1. **Check if root exists**: Verify $f(a)$ and $f(b)$ have opposite signs
2. **Widen the bracket**: The root may be outside the initial interval
3. **Check for discontinuities**: Function must be continuous
4. **Use more iterations**: Increase `MaxIterations`

### Slow Convergence

1. **Multiple root**: If $f(x) = (x-r)^m g(x)$ with $m > 1$, Newton's method slows to linear
2. **Poor initial guess**: Start closer to the root
3. **Nearly horizontal function**: $|f'(x)|$ is small near the root

### Newton's Method Diverges

1. **Bad initial guess**: Try a different starting point
2. **Stationary point nearby**: $f'(x) \approx 0$ causes large jumps
3. **Use damped Newton**: $x_{n+1} = x_n - \alpha \frac{f(x_n)}{f'(x_n)}$ with $\alpha < 1$

```cs
// Damped Newton-Raphson
double alpha = 0.5;  // Damping factor
x = x - alpha * f(x) / df(x);
```

---

## Comparison of Methods

| Method | Evaluations/Iter | Convergence Order | Robustness |
|--------|------------------|-------------------|------------|
| Bisection | 1 | 1.0 | Guaranteed |
| Regula Falsi | 1 | ~1.0-1.6 | Good |
| Brent | 1-2 | ~1.6 | Excellent |
| Newton | 2 (f and f') | 2.0 | Moderate |
| Secant | 1 | ~1.618 | Moderate |

**Recommendation**:
- **General use**: Brent's method
- **Speed critical, smooth function**: Newton-Raphson
- **No derivative, good initial guesses**: Secant method
- **Maximum robustness**: Bisection

---

## References

<a id="ref1">[1]</a> Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

<a id="ref2">[2]</a> Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*. Prentice-Hall.

<a id="ref3">[3]</a> Burden, R. L., & Faires, J. D. (2010). *Numerical Analysis* (9th ed.). Brooks/Cole.
