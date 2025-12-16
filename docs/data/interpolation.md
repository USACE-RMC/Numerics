# Data and Interpolation

[← Back to Index](../index.md)

The ***Numerics*** library provides interpolation methods for estimating values between known data points, essential for data analysis, curve fitting, and function approximation.

## Available Interpolation Methods

| Method | Class | Use Case |
|--------|-------|----------|
| **Linear** | `Linear` | Fast, simple, C⁰ continuous |
| **Cubic Spline** | `CubicSpline` | Smooth curves, C² continuous |
| **Polynomial** | `Polynomial` | Arbitrary order fitting |
| **Bilinear** | `Bilinear` | 2D interpolation on grids |

## Linear Interpolation

Simplest method - connects points with straight lines:

```cs
using Numerics.Data.Interpolation;

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
- Fast: O(log n) per evaluation
- C⁰ continuous (values continuous, derivatives not)
- No overshooting
- Good for piecewise linear trends

## Cubic Spline Interpolation

Smooth curves passing through all data points:

```cs
using Numerics.Data.Interpolation;

double[] xData = { 0, 1, 2, 3, 4, 5 };
double[] yData = { 1, 3, 2, 5, 4, 6 };

// Natural cubic spline (zero second derivatives at endpoints)
var spline = new CubicSpline(xData, yData);

// Interpolate
double y = spline.Interpolate(2.5);
Console.WriteLine($"Spline y(2.5) = {y:F2}");

// Evaluate derivative
double dy = spline.Differentiate(2.5);
Console.WriteLine($"dy/dx(2.5) = {dy:F2}");

// Second derivative
double d2y = spline.Differentiate2(2.5);
Console.WriteLine($"d²y/dx²(2.5) = {d2y:F2}");
```

**Properties:**
- C² continuous (smooth second derivative)
- Unique solution through all points
- Natural boundary conditions
- May overshoot between points
- Excellent for smooth physical phenomena

### Boundary Conditions

```cs
// Natural spline (second derivative = 0 at endpoints)
var naturalSpline = new CubicSpline(xData, yData, 
    boundaryType: CubicSpline.BoundaryType.Natural);

// Clamped spline (specify first derivatives at endpoints)
double leftDerivative = 0.5;
double rightDerivative = 0.8;
var clampedSpline = new CubicSpline(xData, yData,
    boundaryType: CubicSpline.BoundaryType.Clamped,
    leftBoundaryValue: leftDerivative,
    rightBoundaryValue: rightDerivative);

// Not-a-knot spline (third derivative continuous at second and penultimate points)
var notAKnotSpline = new CubicSpline(xData, yData,
    boundaryType: CubicSpline.BoundaryType.NotAKnot);
```

## Polynomial Interpolation

Fits polynomial of specified degree:

```cs
using Numerics.Data.Interpolation;

double[] xData = { 0, 1, 2, 3, 4 };
double[] yData = { 1, 3, 2, 5, 4 };

// Fit 3rd degree polynomial
var poly = new Polynomial(xData, yData, degree: 3);

// Interpolate
double y = poly.Interpolate(2.5);
Console.WriteLine($"Polynomial y(2.5) = {y:F2}");

// Get polynomial coefficients
double[] coeffs = poly.Coefficients;
Console.WriteLine("Polynomial: y = " + 
    string.Join(" + ", coeffs.Select((c, i) => $"{c:F3}x^{i}")));
```

**Warning:** High-degree polynomials (> 5) can exhibit Runge's phenomenon (oscillations).

**Best practice:** Use splines instead of high-degree polynomials.

## Bilinear Interpolation

For 2D gridded data:

```cs
using Numerics.Data.Interpolation;

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

// Multiple points
double[] xNew = { 0.5, 1.5 };
double[] yNew = { 0.5, 1.5 };
double[] zNew = bilinear.Interpolate(xNew, yNew);

for (int i = 0; i < xNew.Length; i++)
{
    Console.WriteLine($"z({xNew[i]}, {yNew[i]}) = {zNew[i]:F2}");
}
```

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
Console.WriteLine("Spline creates smooth curve");
```

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

// Create contour at specific elevation
double contourElevation = 110;
Console.WriteLine($"\nFinding points at {contourElevation} m elevation...");
// Would need to scan grid and find where z = contourElevation
```

## Best Practices

1. **Data spacing**: Interpolation works best with reasonably uniform spacing
2. **Extrapolation**: Avoid extrapolating beyond data range (highly unreliable)
3. **Smoothness**: Use splines for smooth physical phenomena, linear for piecewise trends
4. **Outliers**: Check for data errors before interpolating
5. **Monotonicity**: If data should be monotonic, consider specialized methods
6. **Periodic data**: Consider Fourier or trigonometric interpolation

## Choosing an Interpolation Method

| Data Characteristics | Recommended Method |
|---------------------|-------------------|
| Few points, simple trend | Linear |
| Smooth physical process | Cubic Spline |
| Need derivatives | Cubic Spline |
| Noisy data | Linear or smoothing spline |
| 2D regular grid | Bilinear |
| Piecewise constant | Nearest neighbor |
| Exact polynomial | Polynomial (low degree) |

## Common Pitfalls

1. **Runge's phenomenon**: High-degree polynomials oscillate wildly
2. **Extrapolation**: Results outside data range are unreliable
3. **Unequal spacing**: Some methods assume uniform spacing
4. **Monotonicity**: Splines may violate monotonicity of data
5. **Edge effects**: Interpolation near boundaries less accurate

---

[← Back to Index](../index.md)
