# Linear Algebra

[← Previous: Root Finding](root-finding.md) | [Back to Index](../index.md) | [Next: Special Functions →](special-functions.md)

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

// Diagonal matrix from a vector (static method creates a diagonal matrix)
var D = Matrix.Diagonal(new Vector(new[] { 1.0, 2.0, 3.0 }));

// Extract diagonal elements from an existing matrix (instance method)
double[] diagElements = D.Diagonal();  // Returns [1, 2, 3]

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

// Get row (returns double[])
double[] row1 = A.Row(1);  // Second row: [4, 5, 6]

// Get column (returns double[])
double[] col2 = A.Column(2);  // Third column: [3, 6, 9]

// Set row values using direct indexing
for (int j = 0; j < A.NumberOfColumns; j++)
    A[0, j] = new[] { 10.0, 11.0, 12.0 }[j];

// Set column values using direct indexing
for (int i = 0; i < A.NumberOfRows; i++)
    A[i, 1] = new[] { 20.0, 21.0, 22.0 }[i];

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

// Clone a vector
var v3 = v1.Clone();

Console.WriteLine($"Vector v1: {v1.ToString()}");
Console.WriteLine($"Length: {v1.Length}");
```

### Vector Operations

#### Dot Product

```cs
var a = new Vector(new[] { 1.0, 2.0, 3.0 });
var b = new Vector(new[] { 4.0, 5.0, 6.0 });

double dot = Vector.DotProduct(a, b);  // 1*4 + 2*5 + 3*6 = 32

Console.WriteLine($"a · b = {dot}");
```

#### Norm (Magnitude)

```cs
var v = new Vector(new[] { 3.0, 4.0 });

double norm = v.Norm();  // √(3² + 4²) = 5

Console.WriteLine($"||v|| = {norm}");

// Unit vector (normalize manually)
var u = v / v.Norm();  // u = v / ||v||

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
var v5 = v1 * 2.0;

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

## Matrix Decompositions

The library provides standard matrix decomposition algorithms for solving linear systems, computing determinants, and analyzing matrix structure. Matrix decompositions are the workhorses of numerical linear algebra -- rather than computing a matrix inverse directly (which is both expensive and numerically fragile), decompositions factor a matrix into structured components that can be used to solve systems, compute determinants, and analyze matrix properties efficiently and stably.

### LU Decomposition

#### Mathematical Background

LU decomposition factors a square matrix $A$ into the product of a lower triangular matrix $L$ and an upper triangular matrix $U$. In practice, partial pivoting is applied to improve numerical stability, yielding the factorization:

```math
PA = LU
```

where:

- $P$ is a permutation matrix that records row interchanges,
- $L$ is a lower triangular matrix with ones on the diagonal (unit lower triangular),
- $U$ is an upper triangular matrix.

**Why partial pivoting matters.** Without pivoting, small diagonal elements can appear during Gaussian elimination, causing division by near-zero values that amplify rounding errors catastrophically. Partial pivoting selects the largest available element in each column as the pivot, keeping the multipliers in $L$ bounded and ensuring numerical stability for most practical problems.

**Solving linear systems.** Once the factorization $PA = LU$ is computed, solving $Ax = b$ is reduced to two triangular solves:

```math
Ly = Pb \quad \text{(forward substitution)}
```

```math
Ux = y \quad \text{(back substitution)}
```

Each triangular solve costs only $O(n^2)$ operations. This makes LU decomposition especially efficient when solving multiple systems with the same coefficient matrix but different right-hand sides ($Ax = b_1, Ax = b_2, \ldots$), because the $O(n^3/3)$ factorization is performed only once.

**Computational cost:** $O(n^3/3)$ for the factorization, plus $O(n^2)$ per solve.

#### API Reference

The `LUDecomposition` class takes a square matrix in its constructor and immediately performs the factorization using outer-product Gaussian elimination with partial pivoting. The original matrix is not modified; an internal copy is used.

```cs
using Numerics.Mathematics.LinearAlgebra;

var A = new Matrix(new double[,] {
    { 2, 1, 1 },
    { 4, 3, 3 },
    { 8, 7, 9 }
});

var lu = new LUDecomposition(A);

// Solve Ax = b
var b = new Vector(new double[] { 1, 1, 1 });
Vector x = lu.Solve(b);

// Compute determinant
double det = lu.Determinant();
Console.WriteLine($"det(A) = {det:F4}");

// Compute inverse
Matrix Ainv = lu.InverseA();
```

**Key members:**

| Member | Type | Description |
|--------|------|-------------|
| `LU` | `Matrix` | The combined L and U factors stored in a single matrix |
| `A` | `Matrix` | A copy of the original input matrix |
| `Solve(Vector b)` | `Vector` | Solves $Ax = b$ using forward and back substitution |
| `Solve(Matrix B)` | `Matrix` | Solves $AX = B$ for multiple right-hand sides |
| `Determinant()` | `double` | Computes $\det(A)$ from the product of diagonal elements of $U$ |
| `InverseA()` | `Matrix` | Computes $A^{-1}$ by solving $AX = I$ |

### QR Decomposition

#### Mathematical Background

QR decomposition factors an $m \times n$ matrix $A$ (with $m \geq n$) into the product of an orthogonal matrix and an upper triangular matrix:

```math
A = QR
```

where:

- $Q$ is an $m \times m$ orthogonal matrix, meaning $Q^T Q = I$ (its columns are orthonormal),
- $R$ is an $m \times n$ upper triangular matrix.

The ***Numerics*** library computes the QR factorization using **Householder reflections**. At each step $k$, a Householder matrix $H_k = I - \beta v v^T$ is applied to zero out the subdiagonal entries of column $k$. The product $Q = H_1 H_2 \cdots H_n$ is accumulated explicitly as the full orthogonal factor.

**Application to least squares.** QR decomposition is the preferred method for solving overdetermined least squares problems:

```math
\min_x \| Ax - b \|_2
```

Because $Q$ is orthogonal, the 2-norm is preserved under multiplication by $Q^T$, and the problem reduces to solving the upper triangular system:

```math
Rx = Q^T b
```

This approach is significantly more numerically stable than forming and solving the normal equations $A^T A x = A^T b$ directly. The normal equations square the condition number of $A$ (i.e., $\kappa(A^T A) = \kappa(A)^2$), which can cause severe loss of accuracy for ill-conditioned problems. The QR approach avoids this squaring entirely.

**Computational cost:** $O(2mn^2 - 2n^3/3)$ for an $m \times n$ matrix.

#### API Reference

The `QRDecomposition` class takes a general real $m \times n$ matrix in its constructor and computes the factorization using Householder reflections.

```cs
var A = new Matrix(new double[,] {
    { 1, 1 },
    { 1, 2 },
    { 1, 3 }
});

var qr = new QRDecomposition(A);

// Access factors
Matrix Q = qr.Q;       // Orthogonal matrix
Matrix R = qr.RMatrix;  // Upper triangular

// Solve overdetermined system (least squares)
var b = new Vector(new double[] { 1, 2, 2 });
Vector x = qr.Solve(b);
```

**Key members:**

| Member | Type | Description |
|--------|------|-------------|
| `Q` | `Matrix` | The $m \times m$ orthogonal matrix |
| `RMatrix` | `Matrix` | The $m \times n$ upper triangular matrix |
| `Solve(Vector b)` | `Vector` | Solves $Ax = b$ in the least squares sense via back substitution on $Rx = Q^T b$ |
| `Solve(Matrix B)` | `Matrix` | Solves $AX = B$ for multiple right-hand sides |

### Cholesky Decomposition

#### Mathematical Background

Cholesky decomposition factors a **symmetric positive-definite** (SPD) matrix $A$ into the product of a lower triangular matrix and its transpose:

```math
A = LL^T
```

where $L$ is a lower triangular matrix with strictly positive diagonal entries.

A matrix $A$ is symmetric positive-definite if $A = A^T$ and $x^T A x > 0$ for all nonzero vectors $x$. Common examples of SPD matrices include covariance matrices, correlation matrices, Gram matrices ($X^T X$ for full-rank $X$), and stiffness matrices in structural analysis.

**Computational advantage.** Because the factorization exploits symmetry, Cholesky requires approximately half the work of LU decomposition:

```math
\text{Cost} \approx O(n^3/6)
```

This makes it the fastest general-purpose direct solver for SPD systems.

**Solving linear systems.** As with LU, the solution of $Ax = b$ proceeds by two triangular solves:

```math
Ly = b \quad \text{(forward substitution)}
```

```math
L^T x = y \quad \text{(back substitution)}
```

**Connection to multivariate Normal sampling.** Cholesky decomposition is essential for generating correlated random variables. If $\Sigma = LL^T$ is a covariance matrix and $z \sim N(0, I)$ is a vector of independent standard Normal samples, then:

```math
x = \mu + Lz \quad \implies \quad x \sim N(\mu, \Sigma)
```

This transformation is used extensively in Monte Carlo simulation, Bayesian inference, and risk analysis.

**Positive-definiteness diagnostic.** If the Cholesky decomposition fails (a diagonal element becomes zero or negative during factorization), the matrix is not positive-definite. This is a useful diagnostic for detecting ill-conditioned or indefinite covariance matrices.

#### API Reference

The `CholeskyDecomposition` class takes a symmetric positive-definite matrix in its constructor. If the matrix is not SPD, the constructor throws an exception.

```cs
// Symmetric positive-definite matrix (e.g., covariance matrix)
var A = new Matrix(new double[,] {
    { 4, 2 },
    { 2, 3 }
});

var chol = new CholeskyDecomposition(A);

// Check if decomposition succeeded
Console.WriteLine($"Positive definite: {chol.IsPositiveDefinite}");

// Lower triangular factor
Matrix L = chol.L;

// Solve Ax = b
var b = new Vector(new double[] { 1, 1 });
Vector x = chol.Solve(b);

// Compute log-determinant (numerically stable)
double logDet = chol.LogDeterminant();
Console.WriteLine($"log|A| = {logDet:F4}");
```

**Key members:**

| Member | Type | Description |
|--------|------|-------------|
| `L` | `Matrix` | The lower triangular Cholesky factor |
| `A` | `Matrix` | A copy of the original input matrix |
| `IsPositiveDefinite` | `bool` | `true` if the decomposition succeeded (matrix is SPD) |
| `Solve(Vector b)` | `Vector` | Solves $Ax = b$ via forward and back substitution on $LL^T$ |
| `Forward(Vector b)` | `Vector` | Solves the forward substitution step $Ly = b$ |
| `Backward(Vector y)` | `Vector` | Solves the back substitution step $L^T x = y$ |
| `InverseA()` | `Matrix` | Computes $A^{-1}$ using the stored decomposition |
| `Determinant()` | `double` | Computes $\det(A) = (\prod L_{ii})^2$ |
| `LogDeterminant()` | `double` | Computes $\log \det(A) = 2 \sum \log L_{ii}$ (numerically stable for large matrices) |

### Eigenvalue Decomposition

#### Mathematical Background

An eigenvalue $\lambda$ and corresponding eigenvector $v$ of a square matrix $A$ satisfy the fundamental equation:

```math
Av = \lambda v
```

That is, the matrix $A$ acts on the eigenvector $v$ by simply scaling it. The eigenvalues are the roots of the **characteristic polynomial**:

```math
\det(A - \lambda I) = 0
```

For a real **symmetric** matrix $A$, the eigenvalues are all real and the eigenvectors are orthogonal. This yields the **spectral decomposition**:

```math
A = Q \Lambda Q^T
```

where $Q$ is an orthogonal matrix whose columns are the eigenvectors and $\Lambda = \text{diag}(\lambda_1, \lambda_2, \ldots, \lambda_n)$ is a diagonal matrix of eigenvalues.

The ***Numerics*** library computes this decomposition using the **Jacobi rotation method**, which is an iterative algorithm that applies a sequence of plane rotations to systematically zero out off-diagonal elements. Each rotation is chosen to eliminate the largest remaining off-diagonal entry. The method converges when all off-diagonal elements are smaller than a tolerance ($10^{-12}$). The Jacobi method is robust and highly accurate for small to medium-sized symmetric matrices.

**Applications:**

- **Principal Component Analysis (PCA):** The eigenvectors of a covariance matrix define the principal directions of variation in data.
- **Stability analysis:** In dynamical systems, the eigenvalues of the system matrix determine whether the system is stable (all eigenvalues have negative real parts), neutrally stable, or unstable.
- **Modal analysis:** In structural engineering, the eigenvalues of the stiffness-mass system correspond to natural frequencies of vibration.

**Connection to condition number.** For symmetric matrices, the condition number can be computed directly from the eigenvalues:

```math
\kappa(A) = \left| \frac{\lambda_{\max}}{\lambda_{\min}} \right|
```

A large condition number indicates that the matrix is nearly singular and that solutions to linear systems involving $A$ may be unreliable.

#### API Reference

The `EigenValueDecomposition` class accepts a **symmetric** matrix and computes all eigenvalues and eigenvectors using the Jacobi rotation method. The constructor validates that the input matrix is both square and symmetric.

```cs
var A = new Matrix(new double[,] {
    { 2, 1 },
    { 1, 3 }
});

var eigen = new EigenValueDecomposition(A);

// Eigenvalues
Console.WriteLine("Eigenvalues:");
for (int i = 0; i < eigen.EigenValues.Length; i++)
    Console.WriteLine($"  λ{i} = {eigen.EigenValues[i]:F4}");

// Eigenvectors (columns of EigenVectors matrix)
Console.WriteLine("Eigenvectors:");
Matrix V = eigen.EigenVectors;
```

**Key members:**

| Member | Type | Description |
|--------|------|-------------|
| `EigenValues` | `Vector` | The eigenvalues $\lambda_1, \ldots, \lambda_n$ |
| `EigenVectors` | `Matrix` | Orthogonal matrix whose columns are the corresponding eigenvectors |
| `A` | `Matrix` | A copy of the original input matrix |
| `EffectiveSampleSize()` | `double` | Returns the effective sample size based on Dutilleul's method, computed as $(\sum \lambda_i)^2 / \sum \lambda_i^2$ |

### Singular Value Decomposition (SVD)

#### Mathematical Background

The Singular Value Decomposition is arguably the most important and versatile matrix factorization in numerical linear algebra. Every real $m \times n$ matrix $A$ (regardless of shape or rank) admits the decomposition:

```math
A = U \Sigma V^T
```

where:

- $U$ is an $m \times m$ orthogonal matrix whose columns are the **left singular vectors**,
- $\Sigma$ is an $m \times n$ diagonal matrix containing the **singular values** $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_{\min(m,n)} \geq 0$, arranged in non-increasing order,
- $V$ is an $n \times n$ orthogonal matrix whose columns are the **right singular vectors**.

The singular values are always non-negative real numbers. They represent the "stretching factors" along the principal axes of the linear transformation defined by $A$.

**Relationship to eigenvalues.** The singular values of $A$ are the square roots of the eigenvalues of $A^T A$ (or equivalently $A A^T$). The right singular vectors are the eigenvectors of $A^T A$, and the left singular vectors are the eigenvectors of $A A^T$.

**Pseudoinverse.** The SVD provides a natural way to compute the Moore-Penrose pseudoinverse $A^+$. If $A = U \Sigma V^T$, then:

```math
A^+ = V \Sigma^+ U^T
```

where $\Sigma^+$ is formed by taking the reciprocal of each nonzero singular value on the diagonal. Singular values below a threshold are treated as zero, providing a numerically stable inverse even for rank-deficient matrices.

**Numerical rank.** The number of singular values above a threshold determines the numerical rank of the matrix. This is far more reliable than computing rank via Gaussian elimination, which can be sensitive to rounding errors.

**Low-rank approximation.** The Eckart-Young theorem states that the best rank-$k$ approximation to $A$ (in both the 2-norm and Frobenius norm) is obtained by retaining only the $k$ largest singular values:

```math
A_k = \sum_{i=1}^{k} \sigma_i \, u_i \, v_i^T
```

**Condition number.** The condition number of $A$ in the 2-norm is the ratio of the largest to smallest singular value:

```math
\kappa(A) = \frac{\sigma_{\max}}{\sigma_{\min}}
```

**Applications:**

- **Least squares (rank-deficient):** SVD provides the minimum-norm least squares solution even when $A$ is rank-deficient, where QR decomposition would fail due to a singular $R$ factor.
- **Principal Component Analysis:** The right singular vectors define the principal components; the singular values quantify the variance explained by each component.
- **Signal/noise separation:** In data analysis, large singular values correspond to signal and small singular values correspond to noise.
- **Dimensionality reduction:** Truncated SVD retains only the most important components of a dataset.

#### API Reference

The `SingularValueDecomposition` class takes a general real $m \times n$ matrix and computes the full SVD. The implementation uses Golub-Kahan bidiagonalization followed by implicit QR iteration, and the singular values are automatically sorted in decreasing order.

```cs
using Numerics.Mathematics.LinearAlgebra;

var A = new Matrix(new double[,] {
    { 1, 2 },
    { 3, 4 },
    { 5, 6 }
});

var svd = new SingularValueDecomposition(A);

// Singular values (sorted in decreasing order)
Console.WriteLine("Singular values:");
for (int i = 0; i < svd.W.Length; i++)
    Console.WriteLine($"  σ{i} = {svd.W[i]:F6}");

// Left singular vectors (m x m)
Matrix U = svd.U;

// Right singular vectors (n x n)
Matrix V = svd.V;

// Condition number (reciprocal)
Console.WriteLine($"1/κ(A) = {svd.InverseCondition:E4}");

// Numerical rank
int rank = svd.Rank();
Console.WriteLine($"Rank: {rank}");

// Solve Ax = b using the pseudoinverse
var b = new Vector(new double[] { 1, 2, 3 });
Vector x = svd.Solve(b);
Console.WriteLine($"Least squares solution: {x.ToString()}");
```

**Key members:**

| Member | Type | Description |
|--------|------|-------------|
| `U` | `Matrix` | The $m \times n$ matrix of left singular vectors (column-orthogonal) |
| `V` | `Matrix` | The $n \times n$ orthogonal matrix of right singular vectors |
| `W` | `Vector` | The singular values $\sigma_1 \geq \sigma_2 \geq \cdots \geq 0$, stored as a vector |
| `A` | `Matrix` | A copy of the original input matrix |
| `Threshold` | `double` | Threshold below which singular values are treated as zero (default is based on machine precision) |
| `InverseCondition` | `double` | The reciprocal of the condition number: $\sigma_{\min} / \sigma_{\max}$ |
| `Rank(threshold)` | `int` | Number of singular values above the threshold |
| `Nullity(threshold)` | `int` | Number of singular values at or below the threshold |
| `Range(threshold)` | `Matrix` | Orthonormal basis for the column space of $A$ |
| `Nullspace(threshold)` | `Matrix` | Orthonormal basis for the null space of $A$ |
| `Solve(Vector b, threshold)` | `Vector` | Solves $Ax = b$ using the pseudoinverse |
| `Solve(Matrix B, threshold)` | `Matrix` | Solves $AX = B$ for multiple right-hand sides |
| `LogDeterminant()` | `double` | Computes $\sum \log \sigma_i$ |
| `LogPseudoDeterminant()` | `double` | Computes $\sum \log \sigma_i$ for nonzero $\sigma_i$ only |

## Condition Number and Numerical Stability

When solving linear systems or computing matrix operations, the **condition number** is the single most important diagnostic for assessing the reliability of your results. Understanding it is critical for any safety-critical application.

### Definition

The condition number of a matrix $A$ is defined as:

```math
\kappa(A) = \|A\| \cdot \|A^{-1}\|
```

Using the 2-norm (spectral norm), this simplifies to the ratio of the largest to smallest singular value:

```math
\kappa(A) = \frac{\sigma_{\max}}{\sigma_{\min}}
```

For symmetric matrices, this is equivalently the ratio of the largest to smallest eigenvalue in absolute value.

### Interpretation

The condition number quantifies how much errors in the input data (the matrix $A$ or the right-hand side $b$) are amplified in the solution. Specifically:

- A relative perturbation of size $\epsilon$ in $A$ or $b$ can produce a relative error of up to $\kappa(A) \cdot \epsilon$ in the solution $x$.
- In floating-point arithmetic with approximately 16 digits of precision (double), a condition number of $10^k$ means you may **lose up to $k$ digits of accuracy** in the computed solution.

### Practical Guidelines

| Condition Number | Assessment | Action |
|-----------------|------------|--------|
| $\kappa \approx 1$ | Perfectly conditioned | Results are reliable |
| $\kappa \approx 10^3$ | Mildly ill-conditioned | Results are likely fine for most applications |
| $\kappa \approx 10^6$ | Moderately ill-conditioned | Results should be verified independently |
| $\kappa > 10^{10}$ | Severely ill-conditioned | Results are suspect; consider reformulating the problem |
| $\kappa \approx 10^{16}$ or $\sigma_{\min} \approx 0$ | Numerically singular | The matrix is effectively singular in double precision |

### How to Check

The most reliable way to compute the condition number is through the SVD:

```cs
var A = new Matrix(new double[,] {
    { 1, 2 },
    { 3, 4 }
});

var svd = new SingularValueDecomposition(A);

// Reciprocal condition number (values near 0 indicate ill-conditioning)
double rcond = svd.InverseCondition;
Console.WriteLine($"1/κ(A) = {rcond:E4}");

// Full condition number
double cond = 1.0 / rcond;
Console.WriteLine($"κ(A) = {cond:F2}");

// Inspect individual singular values
Console.WriteLine("Singular values:");
for (int i = 0; i < svd.W.Length; i++)
    Console.WriteLine($"  σ{i} = {svd.W[i]:E6}");
```

**When to check.** Always check the condition number before trusting the results of a linear system solve, especially when:

- The matrix comes from measured or uncertain data,
- The problem involves interpolation or extrapolation,
- The matrix is large and its structure is not well understood,
- Results are used in safety-critical decisions.

## Choosing a Decomposition

Selecting the right decomposition depends on the structure of your matrix and the problem you are solving. The following table provides a practical decision guide:

| Problem | Matrix Properties | Recommended Decomposition | Rationale |
|---------|------------------|--------------------------|-----------|
| Solve $Ax = b$ (general) | Square, non-singular | **LU** | Fast $O(n^3/3)$ factorization; efficient for multiple right-hand sides |
| Solve $Ax = b$ (SPD) | Symmetric positive-definite | **Cholesky** | Half the cost of LU; exploits symmetry |
| Least squares $\min \|Ax - b\|_2$ | Overdetermined ($m > n$), full rank | **QR** | Numerically stable; avoids squaring the condition number |
| Least squares (rank-deficient) | Any shape, possibly rank-deficient | **SVD** | Provides minimum-norm solution even when $A$ is rank-deficient |
| Eigenvalues/eigenvectors | Symmetric | **Eigenvalue** | Spectral decomposition via Jacobi rotation |
| Condition number, numerical rank | Any | **SVD** | Singular values directly give $\kappa(A)$ and rank |
| Monte Carlo sampling | Covariance matrix (SPD) | **Cholesky** | $x = \mu + Lz$ generates correlated samples |
| Determinant | Square | **LU** or **Cholesky** | Product of diagonal elements (or their squares for Cholesky) |

### Rules of Thumb

1. **Start with LU** for general square systems. It is the standard workhorse.
2. **Use Cholesky** whenever you know the matrix is symmetric positive-definite. It is faster and the positive-definiteness check is a valuable diagnostic.
3. **Use QR** for least squares problems. Never solve normal equations with matrix inversion for production code -- it amplifies rounding errors.
4. **Use SVD** when you are unsure about the rank or conditioning of your matrix, or when you need the most robust solution possible.
5. **Use Eigenvalue decomposition** when you need the eigenvalues themselves (for spectral analysis, stability, PCA), not just to solve a linear system.

## Common Operations Summary

| Operation | Method | Complexity |
|-----------|--------|------------|
| Matrix multiplication | `A.Multiply(B)` or `A * B` | O(n³) |
| Transpose | `A.Transpose()` | O(n²) |
| Inverse | `A.Inverse()` | O(n³) |
| Determinant | `A.Determinant()` | O(n³) |
| Vector norm | `v.Norm()` | O(n) |
| Dot product | `Vector.DotProduct(v, w)` | O(n) |
| LU factorization | `new LUDecomposition(A)` | O(n³/3) |
| QR factorization | `new QRDecomposition(A)` | O(2mn² - 2n³/3) |
| Cholesky factorization | `new CholeskyDecomposition(A)` | O(n³/6) |
| SVD | `new SingularValueDecomposition(A)` | O(min(mn², m²n)) |
| Eigenvalue (symmetric) | `new EigenValueDecomposition(A)` | O(n³) iterative |

---

## References

<a id="1">[1]</a> Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.

<a id="2">[2]</a> Trefethen, L. N., & Bau, D. (1997). *Numerical Linear Algebra*. SIAM.

<a id="3">[3]</a> Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

---

[← Previous: Root Finding](root-finding.md) | [Back to Index](../index.md) | [Next: Special Functions →](special-functions.md)
