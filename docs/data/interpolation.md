# Interpolation

The ***Numerics*** library provides interpolation methods for estimating values between known data points, from simple linear interpolation to smooth spline methods.

## Overview

| Method | Continuity | Best For |
|--------|------------|----------|
| Linear | C⁰ | Simple, fast |
| Cubic Spline | C² | Smooth curves |
| Akima Spline | C¹ | Data with outliers |
| Polynomial | C∞ | Small datasets |
| Bilinear | C⁰ | 2D gridded data |
| Bicubic | C¹ | Smooth 2D surfaces |

---

## One-Dimensional Interpolation

### Linear Interpolation

Connects points with straight line segments [[1]](#ref1):

```math
y = y_i + \frac{y_{i+1} - y_i}{x_{i+1} - x_i}(x - x_i)
```

```cs
using Numerics.Data.Interpolation;

double[] x = { 0, 1, 2, 3, 4, 5 };
double[] y = { 0, 1, 4, 9, 16, 25 };

var linear = new LinearInterpolation(x, y);

// Interpolate at a single point
double value = linear.Interpolate(2.5);
Console.WriteLine($"f(2.5) ≈ {value:F4}");

// Interpolate at multiple points
double[] xNew = { 0.5, 1.5, 2.5, 3.5, 4.5 };
double[] yNew = linear.Interpolate(xNew);
```

**Properties**:
- Continuous (C⁰) but not smooth (discontinuous first derivative)
- Never overshoots data
- Fast computation

### Cubic Spline Interpolation

Fits piecewise cubic polynomials with continuous second derivatives [[1]](#ref1):

```cs
var spline = new CubicSplineInterpolation(x, y);

double value = spline.Interpolate(2.5);
Console.WriteLine($"f(2.5) ≈ {value:F4}");

// Get derivative at a point
double derivative = spline.Derivative(2.5);
Console.WriteLine($"f'(2.5) ≈ {derivative:F4}");

// Get second derivative
double secondDeriv = spline.SecondDerivative(2.5);
```

**Boundary Conditions**:

```cs
// Natural spline (second derivative = 0 at endpoints)
var natural = new CubicSplineInterpolation(x, y, BoundaryCondition.Natural);

// Clamped spline (specified first derivatives at endpoints)
var clamped = new CubicSplineInterpolation(x, y, BoundaryCondition.Clamped, 
    leftDerivative: 0.0, rightDerivative: 10.0);

// Not-a-knot (default, continuous third derivative at second and second-to-last points)
var notAKnot = new CubicSplineInterpolation(x, y, BoundaryCondition.NotAKnot);
```

**Properties**:
- C² continuous (smooth first and second derivatives)
- Can overshoot/undershoot data
- Good for smooth underlying functions

### Akima Spline

Locally-determined spline that reduces overshoot [[2]](#ref2):

```cs
var akima = new AkimaSplineInterpolation(x, y);

double value = akima.Interpolate(2.5);
```

**Properties**:
- C¹ continuous (smooth first derivative only)
- Less oscillation than cubic spline
- Better for data with outliers or rapid changes

### Polynomial Interpolation

Fits a single polynomial through all points [[1]](#ref1):

```cs
var poly = new PolynomialInterpolation(x, y);

double value = poly.Interpolate(2.5);

// Get the polynomial coefficients
double[] coeffs = poly.Coefficients;
```

**Warning**: High-degree polynomials can oscillate wildly (Runge's phenomenon). Use only for small datasets (n ≤ 10).

### Monotonic Interpolation

Preserves monotonicity of the data (no spurious extrema):

```cs
var monotonic = new MonotonicInterpolation(x, y);

double value = monotonic.Interpolate(2.5);
```

**Properties**:
- Guarantees interpolant is monotonic if data is monotonic
- Useful for cumulative distribution functions
- C¹ continuous

---

## Extrapolation

By default, interpolation methods only work within the data range. To extrapolate:

```cs
var linear = new LinearInterpolation(x, y);
linear.AllowExtrapolation = true;

// Now works outside [x_min, x_max]
double extrapolated = linear.Interpolate(6.0);
```

**Caution**: Extrapolation can be unreliable, especially for splines.

---

## Two-Dimensional Interpolation

### Bilinear Interpolation

Linear interpolation on a 2D grid [[1]](#ref1):

```cs
double[] xGrid = { 0, 1, 2 };
double[] yGrid = { 0, 1, 2 };
double[,] zValues = {
    { 0, 1, 4 },
    { 1, 2, 5 },
    { 4, 5, 8 }
};

var bilinear = new BilinearInterpolation(xGrid, yGrid, zValues);

double value = bilinear.Interpolate(0.5, 0.5);
Console.WriteLine($"f(0.5, 0.5) ≈ {value:F4}");
```

**Algorithm**: 
1. Linear interpolation in x-direction at two y-values
2. Linear interpolation in y-direction between results

### Bicubic Interpolation

Cubic interpolation on a 2D grid for smoother surfaces:

```cs
var bicubic = new BicubicInterpolation(xGrid, yGrid, zValues);

double value = bicubic.Interpolate(0.5, 0.5);

// Partial derivatives
double dfdx = bicubic.PartialDerivativeX(0.5, 0.5);
double dfdy = bicubic.PartialDerivativeY(0.5, 0.5);
```

---

## Interpolation from Irregular Data

### Inverse Distance Weighting (IDW)

For scattered (non-gridded) 2D data:

```cs
double[] xPoints = { 0, 1, 2, 0, 1, 2 };
double[] yPoints = { 0, 0, 0, 1, 1, 1 };
double[] zPoints = { 0, 1, 4, 1, 2, 5 };

var idw = new InverseDistanceWeighting(xPoints, yPoints, zPoints);
idw.Power = 2;  // Inverse square weighting

double value = idw.Interpolate(0.5, 0.5);
```

**Formula**:
```math
z = \frac{\sum_{i=1}^{n} w_i z_i}{\sum_{i=1}^{n} w_i}, \quad w_i = \frac{1}{d_i^p}
```

where $d_i$ is the distance to point $i$ and $p$ is the power parameter.

---

## Applications

### Resampling Time Series

```cs
// Original irregular time series
double[] times = { 0, 1.5, 3.2, 4.0, 6.1, 8.0 };
double[] values = { 10, 15, 12, 18, 14, 20 };

var spline = new CubicSplineInterpolation(times, values);

// Resample to regular intervals
double dt = 0.5;
int n = (int)((times.Last() - times.First()) / dt) + 1;
double[] regularTimes = new double[n];
double[] regularValues = new double[n];

for (int i = 0; i < n; i++)
{
    regularTimes[i] = times.First() + i * dt;
    regularValues[i] = spline.Interpolate(regularTimes[i]);
}
```

### Rating Curve Interpolation

```cs
// Stage-discharge rating curve
double[] stage = { 0, 1, 2, 3, 4, 5 };      // feet
double[] discharge = { 0, 50, 200, 500, 1000, 1800 };  // cfs

// Use monotonic interpolation to ensure increasing discharge
var rating = new MonotonicInterpolation(stage, discharge);

// Convert stage reading to discharge
double observedStage = 2.5;
double estimatedQ = rating.Interpolate(observedStage);
Console.WriteLine($"Stage {observedStage} ft → {estimatedQ:F0} cfs");
```

### Probability Plotting

```cs
// Empirical CDF interpolation
double[] sortedData = data.OrderBy(x => x).ToArray();
int n = sortedData.Length;
double[] plottingPositions = new double[n];

for (int i = 0; i < n; i++)
    plottingPositions[i] = (i + 1 - 0.4) / (n + 0.2);  // Cunnane

// Interpolate to find quantiles
var empiricalCDF = new MonotonicInterpolation(sortedData, plottingPositions);
var inverseCDF = new MonotonicInterpolation(plottingPositions, sortedData);

// Find value at 95th percentile
double q95 = inverseCDF.Interpolate(0.95);
Console.WriteLine($"95th percentile: {q95:F2}");
```

### DEM Interpolation

```cs
// Digital Elevation Model (gridded)
double[] eastings = { 0, 100, 200, 300, 400 };
double[] northings = { 0, 100, 200, 300, 400 };
double[,] elevations = LoadDEM();  // 5x5 grid

var dem = new BicubicInterpolation(eastings, northings, elevations);

// Get elevation at arbitrary point
double x = 150, y = 175;
double elevation = dem.Interpolate(x, y);
Console.WriteLine($"Elevation at ({x}, {y}): {elevation:F1} m");

// Compute slope
double dzdx = dem.PartialDerivativeX(x, y);
double dzdy = dem.PartialDerivativeY(x, y);
double slope = Math.Atan(Math.Sqrt(dzdx * dzdx + dzdy * dzdy)) * 180 / Math.PI;
Console.WriteLine($"Slope: {slope:F1}°");
```

---

## Choosing an Interpolation Method

| Scenario | Recommended Method |
|----------|-------------------|
| Fast, simple | Linear |
| Smooth curves | Cubic Spline |
| Data with outliers | Akima Spline |
| Must preserve monotonicity | Monotonic |
| CDF interpolation | Monotonic |
| Small dataset, exact fit | Polynomial |
| 2D gridded data | Bilinear or Bicubic |
| 2D scattered data | IDW |

---

## Performance Considerations

1. **Precompute coefficients**: Spline setup is O(n), evaluation is O(log n)
2. **Binary search for interval**: All methods use binary search to find the containing interval
3. **Caching**: Reuse interpolation objects for multiple queries

```cs
// Efficient: Create once, use many times
var spline = new CubicSplineInterpolation(x, y);

for (int i = 0; i < 1000000; i++)
{
    double val = spline.Interpolate(queryPoints[i]);
}

// Inefficient: Recreating each time
for (int i = 0; i < 1000000; i++)
{
    var spline = new CubicSplineInterpolation(x, y);  // Don't do this!
    double val = spline.Interpolate(queryPoints[i]);
}
```

---

## Comparison of Methods

```cs
double[] x = { 0, 1, 2, 3, 4, 5 };
double[] y = { 0, 0.8, 0.9, 0.1, -0.8, -1.0 };

var linear = new LinearInterpolation(x, y);
var cubic = new CubicSplineInterpolation(x, y);
var akima = new AkimaSplineInterpolation(x, y);

Console.WriteLine("x      Linear   Cubic    Akima");
Console.WriteLine("---    ------   -----    -----");

for (double xi = 0; xi <= 5; xi += 0.5)
{
    Console.WriteLine($"{xi:F1}    {linear.Interpolate(xi),6:F3}   {cubic.Interpolate(xi),6:F3}   {akima.Interpolate(xi),6:F3}");
}
```

---

## References

<a id="ref1">[1]</a> Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

<a id="ref2">[2]</a> Akima, H. (1970). A new method of interpolation and smooth curve fitting based on local procedures. *Journal of the ACM*, 17(4), 589-602.

<a id="ref3">[3]</a> de Boor, C. (2001). *A Practical Guide to Splines* (Rev. ed.). Springer.
