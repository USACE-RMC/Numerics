# Linear Algebra

The ***Numerics*** library provides matrix and vector classes with comprehensive linear algebra operations, including decompositions, solvers, and eigenvalue computations.

## Overview

| Class | Description |
|-------|-------------|
| `Vector` | 1D array with vector operations |
| `Matrix` | 2D array with matrix operations |
| `LUDecomposition` | LU factorization |
| `QRDecomposition` | QR factorization |
| `CholeskyDecomposition` | Cholesky factorization |
| `SingularValueDecomposition` | SVD |
| `EigenvalueDecomposition` | Eigenvalues and eigenvectors |

---

## Vector Operations

### Creating Vectors

```cs
using Numerics.Mathematics.LinearAlgebra;

// From array
double[] arr = { 1, 2, 3, 4, 5 };
var v1 = new Vector(arr);

// Specified length (zeros)
var v2 = new Vector(5);

// With initial value
var v3 = new Vector(5, 1.0);  // [1, 1, 1, 1, 1]

// Copy constructor
var v4 = new Vector(v1);
```

### Basic Operations

```cs
var a = new Vector(new[] { 1.0, 2.0, 3.0 });
var b = new Vector(new[] { 4.0, 5.0, 6.0 });

// Arithmetic
var sum = a + b;           // Element-wise addition
var diff = a - b;          // Element-wise subtraction
var scaled = 2.0 * a;      // Scalar multiplication
var scaled2 = a / 2.0;     // Scalar division

// Element-wise operations
var product = a * b;       // Element-wise multiplication
var quotient = a / b;      // Element-wise division

Console.WriteLine($"a + b = {sum}");
Console.WriteLine($"2 * a = {scaled}");
```

### Vector Products

```cs
// Dot product (inner product)
double dot = Vector.DotProduct(a, b);  // 1*4 + 2*5 + 3*6 = 32
double dot2 = a.Dot(b);                // Alternative syntax

// Cross product (3D only)
var cross = Vector.CrossProduct(a, b);

Console.WriteLine($"a · b = {dot}");
Console.WriteLine($"a × b = {cross}");
```

### Vector Norms

```cs
var v = new Vector(new[] { 3.0, 4.0 });

double norm1 = v.Norm(1);      // L1 norm: |3| + |4| = 7
double norm2 = v.Norm(2);      // L2 norm: √(9 + 16) = 5
double normInf = v.Norm(double.PositiveInfinity);  // Max norm: 4

// Euclidean norm (L2) is default
double euclidean = v.Norm();   // 5

// Normalize to unit vector
var unit = v.Normalize();
Console.WriteLine($"Unit vector: {unit}, Norm: {unit.Norm():F6}");
```

### Vector Properties

```cs
int length = v.Length;
double sum = v.Sum();
double min = v.Min();
double max = v.Max();
double mean = v.Mean();
```

---

## Matrix Operations

### Creating Matrices

```cs
// From 2D array
double[,] arr = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
var A = new Matrix(arr);

// Specified dimensions (zeros)
var B = new Matrix(3, 4);

// Identity matrix
var I = Matrix.Identity(3);

// Diagonal matrix
var D = Matrix.Diagonal(new[] { 1.0, 2.0, 3.0 });

// From vectors (as columns)
var col1 = new Vector(new[] { 1.0, 2.0, 3.0 });
var col2 = new Vector(new[] { 4.0, 5.0, 6.0 });
var M = Matrix.FromColumns(new[] { col1, col2 });
```

### Element Access

```cs
var A = new Matrix(new double[,] { { 1, 2 }, { 3, 4 } });

// Single element
double val = A[0, 1];  // Row 0, Column 1 = 2
A[1, 0] = 5;           // Set element

// Row and column vectors
Vector row0 = A.GetRow(0);
Vector col1 = A.GetColumn(1);

// Set row or column
A.SetRow(0, new Vector(new[] { 10.0, 20.0 }));
A.SetColumn(1, new Vector(new[] { 30.0, 40.0 }));
```

### Basic Arithmetic

```cs
var A = new Matrix(new double[,] { { 1, 2 }, { 3, 4 } });
var B = new Matrix(new double[,] { { 5, 6 }, { 7, 8 } });

// Addition and subtraction
var C = A + B;
var D = A - B;

// Scalar operations
var E = 2 * A;
var F = A / 2;

// Matrix multiplication
var G = A * B;

// Element-wise multiplication (Hadamard product)
var H = Matrix.ElementWiseMultiply(A, B);
```

### Matrix-Vector Multiplication

```cs
var A = new Matrix(new double[,] { { 1, 2 }, { 3, 4 } });
var x = new Vector(new[] { 1.0, 2.0 });

// A * x
Vector y = A * x;  // [1*1+2*2, 3*1+4*2] = [5, 11]

Console.WriteLine($"A * x = {y}");
```

### Transpose and Properties

```cs
var A = new Matrix(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });

Matrix AT = A.Transpose();

int rows = A.Rows;        // 2
int cols = A.Columns;     // 3
bool isSquare = A.IsSquare;  // false

// Trace (sum of diagonal elements)
var B = new Matrix(new double[,] { { 1, 2 }, { 3, 4 } });
double trace = B.Trace();  // 1 + 4 = 5
```

### Matrix Norms

```cs
var A = new Matrix(new double[,] { { 1, 2 }, { 3, 4 } });

double frobNorm = A.FrobeniusNorm();    // √(1+4+9+16) = √30
double norm1 = A.Norm1();               // Max column sum
double normInf = A.NormInfinity();      // Max row sum
```

---

## Matrix Decompositions

### LU Decomposition

Factors $A = LU$ (or $PA = LU$ with pivoting) [[1]](#ref1):

```cs
var A = new Matrix(new double[,] { { 2, 1, 1 }, { 4, 3, 3 }, { 8, 7, 9 } });

var lu = new LUDecomposition(A);

Matrix L = lu.L;  // Lower triangular
Matrix U = lu.U;  // Upper triangular
Matrix P = lu.P;  // Permutation matrix

// Solve Ax = b
var b = new Vector(new[] { 4.0, 10.0, 24.0 });
Vector x = lu.Solve(b);

Console.WriteLine($"Solution: {x}");

// Determinant
double det = lu.Determinant();
Console.WriteLine($"Determinant: {det}");

// Inverse
Matrix Ainv = lu.Inverse();
```

### QR Decomposition

Factors $A = QR$ where $Q$ is orthogonal and $R$ is upper triangular [[1]](#ref1):

```cs
var A = new Matrix(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } });

var qr = new QRDecomposition(A);

Matrix Q = qr.Q;  // Orthogonal (Q^T Q = I)
Matrix R = qr.R;  // Upper triangular

// Solve least squares: min ||Ax - b||
var b = new Vector(new[] { 1.0, 2.0, 3.0 });
Vector x = qr.Solve(b);

Console.WriteLine($"Least squares solution: {x}");
```

### Cholesky Decomposition

For symmetric positive definite matrices, $A = LL^T$ [[1]](#ref1):

```cs
// Symmetric positive definite matrix
var A = new Matrix(new double[,] { { 4, 2 }, { 2, 5 } });

var chol = new CholeskyDecomposition(A);

if (chol.IsPositiveDefinite)
{
    Matrix L = chol.L;
    
    // Solve Ax = b
    var b = new Vector(new[] { 8.0, 11.0 });
    Vector x = chol.Solve(b);
    
    Console.WriteLine($"Solution: {x}");
}
```

**Use for**: Covariance matrices, normal equations, simulation.

### Singular Value Decomposition (SVD)

Factors $A = U\Sigma V^T$ [[1]](#ref1):

```cs
var A = new Matrix(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } });

var svd = new SingularValueDecomposition(A);

Matrix U = svd.U;           // Left singular vectors
Matrix V = svd.V;           // Right singular vectors
double[] sigma = svd.S;     // Singular values

Console.WriteLine($"Singular values: {string.Join(", ", sigma.Select(s => s.ToString("F4")))}");

// Rank (number of non-zero singular values)
int rank = svd.Rank;

// Condition number
double cond = svd.ConditionNumber;
Console.WriteLine($"Condition number: {cond:F2}");

// Pseudoinverse (Moore-Penrose)
Matrix Apinv = svd.PseudoInverse();
```

### Eigenvalue Decomposition

For square matrices, finds eigenvalues and eigenvectors [[1]](#ref1):

```cs
var A = new Matrix(new double[,] { { 4, -2 }, { 1, 1 } });

var eig = new EigenvalueDecomposition(A);

double[] realEigenvalues = eig.RealEigenvalues;
double[] imagEigenvalues = eig.ImaginaryEigenvalues;
Matrix V = eig.Eigenvectors;

Console.WriteLine("Eigenvalues:");
for (int i = 0; i < realEigenvalues.Length; i++)
{
    if (Math.Abs(imagEigenvalues[i]) < 1e-10)
        Console.WriteLine($"  λ{i+1} = {realEigenvalues[i]:F4}");
    else
        Console.WriteLine($"  λ{i+1} = {realEigenvalues[i]:F4} + {imagEigenvalues[i]:F4}i");
}
```

---

## Solving Linear Systems

### Direct Solution

```cs
var A = new Matrix(new double[,] { { 3, 1 }, { 1, 2 } });
var b = new Vector(new[] { 9.0, 8.0 });

// Using LU decomposition (general)
var lu = new LUDecomposition(A);
Vector x = lu.Solve(b);

// Using Cholesky (if A is symmetric positive definite)
var chol = new CholeskyDecomposition(A);
if (chol.IsPositiveDefinite)
    x = chol.Solve(b);
```

### Multiple Right-Hand Sides

```cs
var A = new Matrix(new double[,] { { 3, 1 }, { 1, 2 } });
var B = new Matrix(new double[,] { { 9, 3 }, { 8, 4 } });  // Two RHS columns

var lu = new LUDecomposition(A);
Matrix X = lu.Solve(B);  // Solves AX = B
```

### Least Squares

For overdetermined systems ($m > n$):

```cs
// More equations than unknowns
var A = new Matrix(new double[,] { { 1, 1 }, { 1, 2 }, { 1, 3 } });
var b = new Vector(new[] { 1.0, 2.0, 2.0 });

var qr = new QRDecomposition(A);
Vector x = qr.Solve(b);  // Minimizes ||Ax - b||²

Console.WriteLine($"Least squares solution: {x}");

// Residual
Vector residual = A * x - b;
double residualNorm = residual.Norm();
Console.WriteLine($"Residual norm: {residualNorm:F6}");
```

---

## Special Matrix Operations

### Matrix Inverse

```cs
var A = new Matrix(new double[,] { { 4, 7 }, { 2, 6 } });

Matrix Ainv = A.Inverse();

// Verify: A * A^(-1) ≈ I
Matrix product = A * Ainv;
Console.WriteLine("A * A^(-1):");
Console.WriteLine(product);
```

### Determinant

```cs
double det = A.Determinant();
Console.WriteLine($"det(A) = {det}");
```

### Matrix Power

```cs
// A^n
Matrix A2 = A.Power(2);  // A * A
Matrix A3 = A.Power(3);  // A * A * A
```

### Matrix Exponential

For system dynamics $\dot{x} = Ax$, the solution is $x(t) = e^{At}x_0$:

```cs
Matrix expA = A.Exponential();
```

---

## Applications

### Covariance Matrix

```cs
// Data matrix (n observations × p variables)
var data = new Matrix(new double[,]
{
    { 1, 2, 3 },
    { 4, 5, 6 },
    { 7, 8, 9 },
    { 2, 3, 4 }
});

// Center the data
int n = data.Rows;
for (int j = 0; j < data.Columns; j++)
{
    double mean = data.GetColumn(j).Mean();
    for (int i = 0; i < n; i++)
        data[i, j] -= mean;
}

// Covariance matrix: (1/(n-1)) * X^T * X
Matrix cov = (1.0 / (n - 1)) * data.Transpose() * data;
Console.WriteLine("Covariance matrix:");
Console.WriteLine(cov);
```

### Principal Component Analysis (PCA)

```cs
// Eigendecomposition of covariance matrix
var eig = new EigenvalueDecomposition(cov);

// Principal components (eigenvectors)
Matrix pc = eig.Eigenvectors;

// Variance explained (eigenvalues)
double[] variances = eig.RealEigenvalues;
double totalVar = variances.Sum();

Console.WriteLine("Variance explained:");
for (int i = 0; i < variances.Length; i++)
{
    Console.WriteLine($"  PC{i+1}: {100 * variances[i] / totalVar:F1}%");
}
```

### Linear Regression

```cs
// Design matrix and response
var X = new Matrix(new double[,]
{
    { 1, 1.0 },
    { 1, 2.0 },
    { 1, 3.0 },
    { 1, 4.0 }
});
var y = new Vector(new[] { 2.1, 4.0, 5.9, 8.1 });

// Normal equations: β = (X^T X)^(-1) X^T y
Matrix XtX = X.Transpose() * X;
Vector Xty = X.Transpose() * y;

var chol = new CholeskyDecomposition(XtX);
Vector beta = chol.Solve(Xty);

Console.WriteLine($"Intercept: {beta[0]:F4}");
Console.WriteLine($"Slope: {beta[1]:F4}");
```

---

## Performance Tips

1. **Reuse decompositions**: Compute LU/Cholesky once, solve multiple times
2. **Use Cholesky for SPD matrices**: 2× faster than LU
3. **Avoid explicit inverse**: Use `Solve()` instead of computing inverse
4. **Check condition number**: Poor conditioning leads to inaccurate results

```cs
// Bad: Computing inverse explicitly
Matrix Ainv = A.Inverse();
Vector x = Ainv * b;

// Good: Using decomposition
var lu = new LUDecomposition(A);
Vector x = lu.Solve(b);
```

---

## References

<a id="ref1">[1]</a> Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.

<a id="ref2">[2]</a> Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

<a id="ref3">[3]</a> Trefethen, L. N., & Bau, D. (1997). *Numerical Linear Algebra*. SIAM.
