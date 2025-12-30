# Linear Algebra

[← Back to Index](../index.md) | [Next: Special Functions →](special-functions.md)

The ***Numerics*** library provides `Matrix` and `Vector` classes for linear algebra operations. These classes support common operations needed for numerical computing, optimization, and statistical analysis.

## Matrix Class

### Creating Matrices

```cs
using Numerics.Mathematics.LinearAlgebra;

// Create matrix with dimensions
var m1 = new Matrix(3, 4);  // 3 rows, 4 columns

// Create square matrix
var m2 = new Matrix(3);     // 3x3 matrix

// Create from 2D array
double[,] data = {
    { 1, 2, 3 },
    { 4, 5, 6 },
    { 7, 8, 9 }
};
var m3 = new Matrix(data);

// Create from 1D array (column vector)
double[] columnData = { 1, 2, 3, 4 };
var m4 = new Matrix(columnData);  // 4x1 matrix

// Create from Vector
var vec = new Vector(new[] { 1.0, 2.0, 3.0 });
var m5 = new Matrix(vec);
```

### Special Matrices

```cs
// Identity matrix
var I = Matrix.Identity(3);  // 3x3 identity

// Zero matrix
var zeros = new Matrix(3, 3);  // Initialized to zeros by default

// Diagonal matrix
double[] diag = { 1, 2, 3 };
var D = Matrix.Diagonal(diag);

Console.WriteLine("Identity Matrix:");
Console.WriteLine(I.ToString());
```

### Accessing Elements

```cs
var m = new Matrix(new double[,] {
    { 1, 2, 3 },
    { 4, 5, 6 }
});

// Get/set elements
double value = m[0, 1];  // Get element at row 0, column 1 → 2
m[1, 2] = 10;            // Set element

// Properties
int rows = m.NumberOfRows;        // 2
int cols = m.NumberOfColumns;     // 3
bool isSquare = m.IsSquare;       // false

Console.WriteLine($"Matrix dimensions: {rows} x {cols}");
```

### Matrix Operations

#### Transpose

```cs
var A = new Matrix(new double[,] {
    { 1, 2, 3 },
    { 4, 5, 6 }
});

// Instance method
var AT = A.Transpose();  // 3x2 matrix

// Static method
var AT2 = Matrix.Transpose(A);

Console.WriteLine("A:");
Console.WriteLine(A.ToString());
Console.WriteLine("\nA^T:");
Console.WriteLine(AT.ToString());
```

#### Matrix Multiplication

```cs
var A = new Matrix(new double[,] {
    { 1, 2 },
    { 3, 4 }
});

var B = new Matrix(new double[,] {
    { 5, 6 },
    { 7, 8 }
});

// Matrix-matrix multiplication
var C = A.Multiply(B);  // A * B

// Operator overload
var C2 = A * B;

Console.WriteLine("A * B:");
Console.WriteLine(C.ToString());

// Scalar multiplication
var D = A.Multiply(2.0);  // or A * 2.0

Console.WriteLine("\n2 * A:");
Console.WriteLine(D.ToString());
```

#### Matrix-Vector Multiplication

```cs
var A = new Matrix(new double[,] {
    { 1, 2, 3 },
    { 4, 5, 6 }
});

var x = new Vector(new[] { 1.0, 2.0, 3.0 });

// Matrix-vector multiplication
var b = A.Multiply(x);  // A * x

Console.WriteLine("A * x:");
Console.WriteLine(b.ToString());
```

#### Matrix Addition and Subtraction

```cs
var A = new Matrix(new double[,] {
    { 1, 2 },
    { 3, 4 }
});

var B = new Matrix(new double[,] {
    { 5, 6 },
    { 7, 8 }
});

// Addition
var C = A + B;

// Subtraction
var D = A - B;

Console.WriteLine("A + B:");
Console.WriteLine(C.ToString());

Console.WriteLine("\nA - B:");
Console.WriteLine(D.ToString());
```

### Matrix Properties

#### Determinant

```cs
var A = new Matrix(new double[,] {
    { 4, 3 },
    { 6, 3 }
});

double det = A.Determinant();

Console.WriteLine($"Determinant: {det:F2}");

// For larger matrices
var B = new Matrix(new double[,] {
    { 1, 2, 3 },
    { 0, 1, 4 },
    { 5, 6, 0 }
});

double det2 = B.Determinant();
Console.WriteLine($"Determinant of 3x3: {det2:F2}");
```

#### Inverse

```cs
var A = new Matrix(new double[,] {
    { 4, 7 },
    { 2, 6 }
});

try
{
    var Ainv = A.Inverse();
    
    Console.WriteLine("A:");
    Console.WriteLine(A.ToString());
    Console.WriteLine("\nA^-1:");
    Console.WriteLine(Ainv.ToString());
    
    // Verify: A * A^-1 = I
    var I = A * Ainv;
    Console.WriteLine("\nA * A^-1 (should be I):");
    Console.WriteLine(I.ToString());
}
catch (InvalidOperationException ex)
{
    Console.WriteLine($"Matrix is singular: {ex.Message}");
}
```

### Row and Column Operations

```cs
var A = new Matrix(new double[,] {
    { 1, 2, 3 },
    { 4, 5, 6 },
    { 7, 8, 9 }
});

// Get row
Vector row1 = A.GetRow(1);  // Second row: [4, 5, 6]

// Get column
Vector col2 = A.GetColumn(2);  // Third column: [3, 6, 9]

// Set row
A.SetRow(0, new Vector(new[] { 10.0, 11.0, 12.0 }));

// Set column
A.SetColumn(1, new Vector(new[] { 20.0, 21.0, 22.0 }));

Console.WriteLine("Modified matrix:");
Console.WriteLine(A.ToString());
```

## Vector Class

### Creating Vectors

```cs
using Numerics.Mathematics.LinearAlgebra;

// Create from array
var v1 = new Vector(new[] { 1.0, 2.0, 3.0 });

// Create with size
var v2 = new Vector(5);  // Length 5, initialized to zeros

// Copy constructor
var v3 = new Vector(v1);

Console.WriteLine($"Vector v1: {v1.ToString()}");
Console.WriteLine($"Length: {v1.Length}");
```

### Vector Operations

#### Dot Product

```cs
var a = new Vector(new[] { 1.0, 2.0, 3.0 });
var b = new Vector(new[] { 4.0, 5.0, 6.0 });

double dot = a.DotProduct(b);  // 1*4 + 2*5 + 3*6 = 32

Console.WriteLine($"a · b = {dot}");
```

#### Norm (Magnitude)

```cs
var v = new Vector(new[] { 3.0, 4.0 });

double norm = v.Norm();  // √(3² + 4²) = 5

Console.WriteLine($"||v|| = {norm}");

// Unit vector
var u = v.Normalize();  // u = v / ||v||

Console.WriteLine($"Unit vector: {u.ToString()}");
Console.WriteLine($"||u|| = {u.Norm():F10}");  // Should be 1.0
```

#### Vector Addition and Scaling

```cs
var v1 = new Vector(new[] { 1.0, 2.0, 3.0 });
var v2 = new Vector(new[] { 4.0, 5.0, 6.0 });

// Addition
var v3 = v1 + v2;

// Subtraction
var v4 = v2 - v1;

// Scalar multiplication
var v5 = v1 * 2.0;  // or 2.0 * v1

Console.WriteLine($"v1 + v2 = {v3.ToString()}");
Console.WriteLine($"v2 - v1 = {v4.ToString()}");
Console.WriteLine($"2 * v1 = {v5.ToString()}");
```

### Accessing Elements

```cs
var v = new Vector(new[] { 10.0, 20.0, 30.0, 40.0 });

// Indexing
double x = v[0];  // 10.0
v[2] = 35.0;      // Set third element

Console.WriteLine($"Modified vector: {v.ToString()}");

// Convert to array
double[] array = v.ToArray();
```

## Practical Examples

### Example 1: Solving Linear Systems (Ax = b)

```cs
// Solve: 2x + 3y = 8
//        4x + y = 10

var A = new Matrix(new double[,] {
    { 2, 3 },
    { 4, 1 }
});

var b = new Vector(new[] { 8.0, 10.0 });

// Solve using matrix inversion: x = A^-1 * b
var Ainv = A.Inverse();
var x = Ainv.Multiply(b);

Console.WriteLine("Solution to Ax = b:");
Console.WriteLine($"x = {x[0]:F2}");
Console.WriteLine($"y = {x[1]:F2}");

// Verify solution
var check = A.Multiply(x);
Console.WriteLine($"\nVerification A*x = {check.ToString()}");
Console.WriteLine($"Expected b = {b.ToString()}");
```

### Example 2: Least Squares Regression

```cs
// Fit y = a + bx to data
double[] xData = { 1, 2, 3, 4, 5 };
double[] yData = { 2.1, 3.9, 6.2, 8.1, 9.8 };

int n = xData.Length;

// Build design matrix: X = [1, x]
var X = new Matrix(n, 2);
for (int i = 0; i < n; i++)
{
    X[i, 0] = 1.0;      // Intercept column
    X[i, 1] = xData[i]; // x column
}

var y = new Vector(yData);

// Normal equations: (X^T X) β = X^T y
var XTX = X.Transpose() * X;
var XTy = X.Transpose().Multiply(y);

// Solve for coefficients
var beta = XTX.Inverse().Multiply(XTy);

double intercept = beta[0];
double slope = beta[1];

Console.WriteLine($"Linear regression: y = {intercept:F3} + {slope:F3}x");

// Predictions
Console.WriteLine("\nPredictions:");
for (int i = 0; i < n; i++)
{
    double pred = intercept + slope * xData[i];
    Console.WriteLine($"x={xData[i]}: y_obs={yData[i]:F1}, y_pred={pred:F1}");
}
```

### Example 3: Covariance Matrix

```cs
// Compute covariance matrix of multivariate data
double[,] data = {
    { 1, 2, 3 },  // Observation 1
    { 4, 5, 6 },  // Observation 2
    { 7, 8, 9 },  // Observation 3
    { 2, 3, 4 }   // Observation 4
};

int n = data.GetLength(0);  // Number of observations
int p = data.GetLength(1);  // Number of variables

// Center data (subtract means)
double[] means = new double[p];
for (int j = 0; j < p; j++)
{
    for (int i = 0; i < n; i++)
        means[j] += data[i, j];
    means[j] /= n;
}

var X = new Matrix(n, p);
for (int i = 0; i < n; i++)
    for (int j = 0; j < p; j++)
        X[i, j] = data[i, j] - means[j];

// Covariance matrix: Σ = (1/(n-1)) X^T X
var XTX = X.Transpose() * X;
var Sigma = XTX * (1.0 / (n - 1));

Console.WriteLine("Covariance Matrix:");
Console.WriteLine(Sigma.ToString());

// Correlation matrix
var Corr = new Matrix(p, p);
for (int i = 0; i < p; i++)
{
    for (int j = 0; j < p; j++)
    {
        Corr[i, j] = Sigma[i, j] / Math.Sqrt(Sigma[i, i] * Sigma[j, j]);
    }
}

Console.WriteLine("\nCorrelation Matrix:");
Console.WriteLine(Corr.ToString());
```

### Example 4: Distance Calculations

```cs
var point1 = new Vector(new[] { 1.0, 2.0, 3.0 });
var point2 = new Vector(new[] { 4.0, 6.0, 8.0 });

// Euclidean distance
var diff = point2 - point1;
double distance = diff.Norm();

Console.WriteLine($"Distance between points: {distance:F3}");

// Manhattan distance
double manhattan = 0;
for (int i = 0; i < point1.Length; i++)
{
    manhattan += Math.Abs(point2[i] - point1[i]);
}

Console.WriteLine($"Manhattan distance: {manhattan:F1}");
```

## Performance Considerations

- Matrix operations create new objects - use in-place methods when available
- For large matrices, consider specialized libraries for decompositions
- Inverse is computationally expensive - avoid when possible
- For solving Ax=b, specialized solvers are more efficient than computing A^-1

## Common Operations Summary

| Operation | Method | Complexity |
|-----------|--------|------------|
| Matrix multiplication | `A.Multiply(B)` or `A * B` | O(n³) |
| Transpose | `A.Transpose()` | O(n²) |
| Inverse | `A.Inverse()` | O(n³) |
| Determinant | `A.Determinant()` | O(n³) |
| Vector norm | `v.Norm()` | O(n) |
| Dot product | `v.DotProduct(w)` | O(n) |

---

[← Back to Index](../index.md) | [Next: Special Functions →](special-functions.md)
