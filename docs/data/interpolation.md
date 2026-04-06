# Data and Interpolation

[← Previous: ODE Solvers](../mathematics/ode-solvers.md) | [Back to Index](../index.md) | [Next: Linear Regression →](regression.md)

Interpolation is the process of estimating values between known data points. Given a set of $n$ data points $(x_1, y_1), (x_2, y_2), \ldots, (x_n, y_n)$, interpolation constructs a function $p(x)$ that passes through all the data points, i.e., $p(x_i) = y_i$ for all $i$. This is distinct from regression, which fits a function that approximates the data while minimizing some error criterion.

The ***Numerics*** library provides interpolation methods for estimating values between known data points, essential for data analysis, curve fitting, and function approximation.

## Available Interpolation Methods

| Method | Class | Continuity | Error Order | Use Case |
|--------|-------|------------|-------------|----------|
| **Linear** | `Linear` | $C^0$ | $O(h^2)$ | Fast, simple, no overshooting |
| **Cubic Spline** | `CubicSpline` | $C^2$ | $O(h^4)$ | Smooth curves, physical phenomena |
| **Polynomial** | `Polynomial` | $C^\infty$ | Varies | Arbitrary order fitting |
| **Bilinear** | `Bilinear` | $C^0$ | $O(h^2)$ | 2D interpolation on grids |

## Linear Interpolation

Linear interpolation connects adjacent data points with straight line segments. For a query point $x$ in the interval $[x_i, x_{i+1}]$, the interpolated value is:

```math
p(x) = y_i + \frac{x - x_i}{x_{i+1} - x_i} \cdot (y_{i+1} - y_i)
```

This can be rewritten in terms of the normalized coordinate $t = \frac{x - x_i}{x_{i+1} - x_i}$ as $p(x) = (1-t) \cdot y_i + t \cdot y_{i+1}$, which is simply a weighted average of the two bracketing values.

**Error bound**: For a function $f(x)$ with bounded second derivative, the interpolation error satisfies:

```math
|f(x) - p(x)| \leq \frac{h^2}{8} \max |f''(\xi)|
```

where $h = x_{i+1} - x_i$ is the local spacing. The error is thus $O(h^2)$, meaning halving the data spacing reduces the error by a factor of four.

```cs
using Numerics.Data;

double[] xData = { 0, 1, 2, 3, 4, 5 };
double[] yData = { 1, 3, 2, 5, 4, 6 };

var linear = new Linear(xData, yData);

// Interpolate at new points
double y = linear.Interpolate(2.5);
Console.WriteLine($"y(2.5) = {y:F2}");

// Multiple points
double[] xNew = { 0.5, 1.5, 2.5, 3.5 };
double[] yNew = linear.Interpolate(xNew);

Console.WriteLine("Interpolated values:");
for (int i = 0; i < xNew.Length; i++)
{
    Console.WriteLine($"  y({xNew[i]}) = {yNew[i]:F2}");
}
```

**Properties:**
- Fast: $O(\log n)$ per evaluation (bisection search for interval)
- $C^0$ continuous (values continuous, derivatives not)
- No overshooting — interpolant stays within bracket values
- Good for piecewise linear trends or noisy data

The `Linear` class also supports coordinate transforms via `XTransform` and `YTransform` properties, enabling log-linear or probability-scale interpolation. Available transforms are `Transform.None` (default), `Transform.Logarithmic` ($\log_{10}$), and `Transform.NormalZ` (standard normal quantile).

## Cubic Spline Interpolation

A cubic spline constructs a piecewise cubic polynomial $S(x)$ that passes through all data points and has continuous first and second derivatives everywhere. This smoothness is what distinguishes splines from simple piecewise polynomial interpolation.

### Mathematical Foundation

On each subinterval $[x_i, x_{i+1}]$, the spline is a cubic polynomial. If we denote the second derivatives at the knot points as $M_i = S''(x_i)$, then the spline on $[x_i, x_{i+1}]$ can be written as:

```math
S(x) = \frac{M_i}{6h_i}(x_{i+1} - x)^3 + \frac{M_{i+1}}{6h_i}(x - x_i)^3 + \left(\frac{y_i}{h_i} - \frac{M_i h_i}{6}\right)(x_{i+1} - x) + \left(\frac{y_{i+1}}{h_i} - \frac{M_{i+1} h_i}{6}\right)(x - x_i)
```

where $h_i = x_{i+1} - x_i$. Requiring continuity of the first derivative $S'(x)$ at each interior knot point $x_i$ (for $i = 1, \ldots, n-2$) yields the tridiagonal system:

```math
h_{i-1} M_{i-1} + 2(h_{i-1} + h_i) M_i + h_i M_{i+1} = 6 \left( \frac{y_{i+1} - y_i}{h_i} - \frac{y_i - y_{i-1}}{h_{i-1}} \right)
```

This system of $n-2$ equations in $n$ unknowns requires two boundary conditions. The ***Numerics*** library uses **natural boundary conditions**, setting $M_0 = 0$ and $M_{n-1} = 0$ (zero second derivatives at the endpoints). The resulting tridiagonal system is solved efficiently using the Thomas algorithm (a specialized form of Gaussian elimination for tridiagonal matrices) in $O(n)$ time [[1]](#1).

**Error bound**: For a function $f(x)$ with bounded fourth derivative, the natural cubic spline error satisfies:

```math
|f(x) - S(x)| \leq \frac{5h^4}{384} \max |f^{(4)}(\xi)|
```

The $O(h^4)$ convergence rate means halving the data spacing reduces the error by a factor of sixteen — a significant improvement over linear interpolation.

```cs
using Numerics.Data;

double[] xData = { 0, 1, 2, 3, 4, 5 };
double[] yData = { 1, 3, 2, 5, 4, 6 };

// Natural cubic spline (zero second derivatives at endpoints)
var spline = new CubicSpline(xData, yData);

// Interpolate at a single point
double y = spline.Interpolate(2.5);
Console.WriteLine($"Spline y(2.5) = {y:F2}");

// Interpolate at multiple points
double[] xNew = { 0.5, 1.5, 2.5, 3.5 };
double[] yNew = spline.Interpolate(xNew);
for (int i = 0; i < xNew.Length; i++)
{
    Console.WriteLine($"  y({xNew[i]}) = {yNew[i]:F2}");
}
```

**Properties:**
- $C^2$ continuous (smooth second derivative)
- Unique solution through all points
- Natural boundary conditions ($S''=0$ at endpoints)
- May overshoot between data points, especially with oscillatory data
- Excellent for smooth physical phenomena

## Polynomial Interpolation

Given $n$ data points, there exists a unique polynomial of degree at most $n-1$ that passes through all of them. However, fitting a single high-degree polynomial to many data points is often a poor choice due to Runge's phenomenon (discussed below).

### Neville's Method

The ***Numerics*** library uses Neville's method [[1]](#1) to evaluate the interpolating polynomial. Rather than computing coefficients explicitly, Neville's method builds a tableau of progressively higher-degree polynomial approximations through recursive divided differences:

```math
P_{i,j}(x) = \frac{(x - x_j) P_{i,j-1}(x) - (x - x_i) P_{i+1,j}(x)}{x_i - x_j}
```

where $P_{i,i}(x) = y_i$ are the initial zeroth-degree polynomials. The method also provides an error estimate from the last correction term added to the approximation.

The `Polynomial` class fits a polynomial of specified order using a local window of `order + 1` points centered near the query point, which helps mitigate oscillation issues:

```cs
using Numerics.Data;

double[] xData = { 0, 1, 2, 3, 4 };
double[] yData = { 1, 3, 2, 5, 4 };

// Fit 3rd order polynomial (order is the first parameter)
var poly = new Polynomial(3, xData, yData);

// Interpolate
double y = poly.Interpolate(2.5);
Console.WriteLine($"Polynomial y(2.5) = {y:F2}");

// The error estimate from the most recent interpolation
Console.WriteLine($"Error estimate: {poly.Error:F6}");
```

### Runge's Phenomenon

A critical limitation of polynomial interpolation is Runge's phenomenon: as the polynomial degree increases, the interpolation error can grow dramatically near the edges of the interval, even for smooth functions. Consider the Runge function:

```math
f(x) = \frac{1}{1 + 25x^2}
```

With $n$ equally spaced points on $[-1, 1]$, the degree $n-1$ interpolating polynomial oscillates with increasing amplitude near $x = \pm 1$ as $n$ grows. For $n = 11$ (degree 10), the maximum error near the boundaries exceeds 1.0, even though $f(x)$ is perfectly smooth.

**Mitigation strategies:**
- Use **cubic splines** instead of high-degree polynomials — splines keep each piece low-degree while maintaining global smoothness
- Use **Chebyshev nodes** (non-uniform spacing clustered near endpoints) if you can control where data is collected
- Keep polynomial degree low ($\leq 5$) and use a local window, as the `Polynomial` class does

**Best practice:** Use splines instead of high-degree polynomials.

## Bilinear Interpolation

Bilinear interpolation extends linear interpolation to two-dimensional gridded data. For a query point $(x, y)$ in the cell bounded by grid points $(x_i, y_j)$, $(x_{i+1}, y_j)$, $(x_i, y_{j+1})$, $(x_{i+1}, y_{j+1})$, the interpolated value is computed by performing linear interpolation twice — first along one axis, then along the other:

```math
z(x, y) = (1-t)(1-u) \cdot z_{i,j} + t(1-u) \cdot z_{i+1,j} + t \cdot u \cdot z_{i+1,j+1} + (1-t) \cdot u \cdot z_{i,j+1}
```

where $t = \frac{x - x_i}{x_{i+1} - x_i}$ and $u = \frac{y - y_j}{y_{j+1} - y_j}$ are the normalized coordinates within the cell. The result is independent of the order in which the two linear interpolations are performed.

**Error**: For a function with bounded mixed partial derivative, $|f(x,y) - z(x,y)| = O(h^2)$, similar to 1D linear interpolation.

```cs
using Numerics.Data;

// Grid coordinates
double[] xGrid = { 0, 1, 2 };
double[] yGrid = { 0, 1, 2 };

// Grid values z[i,j] = z(xGrid[i], yGrid[j])
double[,] zGrid = {
    { 1, 2, 3 },
    { 2, 3, 4 },
    { 3, 4, 5 }
};

var bilinear = new Bilinear(xGrid, yGrid, zGrid);

// Interpolate at arbitrary point
double z = bilinear.Interpolate(0.5, 0.5);
Console.WriteLine($"z(0.5, 0.5) = {z:F2}");

// Multiple points (loop over individual point pairs)
double[] xNew = { 0.5, 1.5 };
double[] yNew = { 0.5, 1.5 };

for (int i = 0; i < xNew.Length; i++)
{
    double zi = bilinear.Interpolate(xNew[i], yNew[i]);
    Console.WriteLine($"z({xNew[i]}, {yNew[i]}) = {zi:F2}");
}
```

Like the `Linear` class, `Bilinear` supports coordinate transforms (`X1Transform`, `X2Transform`, `YTransform`) for log-linear or probability-scale interpolation in 2D.

**Applications:**
- Image resizing
- Terrain elevation maps
- Temperature/pressure fields
- Geographic data

## Practical Examples

### Example 1: Stage-Discharge Rating Curve

```cs
// Measured stage-discharge pairs
double[] stage = { 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0 };
double[] discharge = { 1000, 1500, 2200, 3100, 4200, 5500, 7000 };

// Create spline interpolator
var ratingCurve = new CubicSpline(stage, discharge);

// Interpolate discharge for observed stage
double observedStage = 6.3;
double estimatedQ = ratingCurve.Interpolate(observedStage);

Console.WriteLine($"Rating Curve Interpolation:");
Console.WriteLine($"Stage: {observedStage:F1} ft → Discharge: {estimatedQ:F0} cfs");

// Extrapolation warning
if (observedStage < stage.Min() || observedStage > stage.Max())
{
    Console.WriteLine("Warning: Extrapolating beyond data range");
}
```

### Example 2: Time Series Gap Filling

```cs
// Time series with gaps
var dates = new[] { 1.0, 2.0, 3.0, /* gap */ 6.0, 7.0, 8.0 };
var values = new[] { 10.5, 11.2, 10.8, /* gap */ 12.1, 12.5, 11.9 };

var interpolator = new CubicSpline(dates, values);

// Fill gaps
var missingDates = new[] { 4.0, 5.0 };
foreach (var t in missingDates)
{
    double filled = interpolator.Interpolate(t);
    Console.WriteLine($"Day {t}: {filled:F2} (interpolated)");
}
```

### Example 3: Comparing Methods

```cs
double[] x = { 0, 1, 2, 3, 4 };
double[] y = { 0, 1, 0, 1, 0 };  // Oscillating data

var linear = new Linear(x, y);
var spline = new CubicSpline(x, y);

double testPoint = 1.5;
double yLinear = linear.Interpolate(testPoint);
double ySpline = spline.Interpolate(testPoint);

Console.WriteLine($"At x = {testPoint}:");
Console.WriteLine($"  Linear: {yLinear:F3}");
Console.WriteLine($"  Spline: {ySpline:F3}");
Console.WriteLine("\nLinear connects with straight line");
Console.WriteLine("Spline creates smooth curve (may overshoot)");
```

Note the spline may produce values outside the range of the bracketing data points due to the smoothness constraint. For data that must remain monotonic, this overshoot can be problematic — in such cases, consider using linear interpolation instead.

### Example 4: 2D Surface Interpolation

```cs
// Elevation data at grid points
double[] eastings = { 0, 100, 200 };   // meters
double[] northings = { 0, 100, 200 };  // meters
double[,] elevations = {
    { 100, 105, 110 },
    { 102, 108, 115 },
    { 104, 112, 120 }
};

var terrain = new Bilinear(eastings, northings, elevations);

// Interpolate elevation at arbitrary location
double x = 150, y = 150;
double z = terrain.Interpolate(x, y);

Console.WriteLine($"Terrain elevation at ({x}, {y}): {z:F1} m");

// Sample elevations along a transect
Console.WriteLine("\nElevation transect from (2,2) to (8,8):");
for (double t = 0; t <= 1; t += 0.2)
{
    double xi = 2 + 6 * t;
    double yi = 2 + 6 * t;
    double zi = terrain.Interpolate(xi, yi);
    Console.WriteLine($"  ({xi:F1}, {yi:F1}): {zi:F1} m");
}
```

## Best Practices

1. **Data spacing**: Interpolation works best with reasonably uniform spacing. Highly irregular spacing can lead to large errors in some subintervals
2. **Extrapolation**: Avoid extrapolating beyond data range — splines and polynomials are especially unreliable outside the data bounds. The `Linear` class returns boundary values; the `Bilinear` class falls back to 1D interpolation at edges
3. **Smoothness**: Use splines for smooth physical phenomena ($C^2$ continuity), linear for piecewise trends or noisy data ($C^0$ continuity with no overshooting)
4. **Outliers**: Check for data errors before interpolating — splines will faithfully pass through outliers
5. **Monotonicity**: Cubic splines do not preserve monotonicity of the data. If your data should be monotonic (e.g., a CDF or rating curve), verify the interpolant doesn't violate this
6. **Periodic data**: Consider Fourier or trigonometric interpolation for periodic signals

## Choosing an Interpolation Method

| Data Characteristics | Recommended Method | Why |
|---------------------|-------------------|-----|
| Few points, simple trend | Linear | No overshooting, $O(h^2)$ error |
| Smooth physical process | Cubic Spline | $C^2$ smooth, $O(h^4)$ error |
| Need derivatives | Cubic Spline | Spline derivatives are well-defined |
| Noisy data | Linear | Splines amplify noise through curvature matching |
| 2D regular grid | Bilinear | Direct extension to 2D |
| Piecewise constant | Nearest neighbor | Preserves step structure |
| Exact polynomial | Polynomial (low degree) | Neville's method with error estimate |

## Common Pitfalls

1. **Runge's phenomenon**: High-degree polynomials ($>5$) oscillate wildly near interval boundaries — use splines instead
2. **Extrapolation**: Results outside the data range are unreliable for all methods
3. **Natural boundary conditions**: The zero-curvature constraint at endpoints can cause the spline to flatten near the boundaries, which may be physically inappropriate
4. **Monotonicity violation**: Cubic splines can introduce local extrema between data points, violating monotonicity
5. **Edge effects**: All methods lose accuracy near the boundaries of the data range

---

## References

<a id="1">[1]</a> W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery, *Numerical Recipes: The Art of Scientific Computing*, 3rd ed., Cambridge, UK: Cambridge University Press, 2007.

<a id="2">[2]</a> C. de Boor, *A Practical Guide to Splines*, Rev. ed., New York: Springer, 2001.

<a id="3">[3]</a> R. L. Burden and J. D. Faires, *Numerical Analysis*, 9th ed., Boston: Brooks/Cole, 2010.

---

[← Previous: ODE Solvers](../mathematics/ode-solvers.md) | [Back to Index](../index.md) | [Next: Linear Regression →](regression.md)
