# Special functions, linear algebra and Fourier transforms

pynumerics exposes the mathematical building blocks of RMC.Numerics in three places:

- `pynumerics.special`: gamma, beta, error, Bessel and Debye functions, factorials and
  combinations, and polynomial evaluation (`Numerics.Mathematics.SpecialFunctions`).
- `pynumerics.linalg`: the `Matrix` and `Vector` types, the Cholesky, LU, QR, singular value
  and eigenvalue decompositions, and matrix regularization (`Numerics.Mathematics.LinearAlgebra`).
- `pynumerics.mathematics`: the fast Fourier transforms (`Numerics.Mathematics.Fourier`) and the
  numerical utilities of `Numerics.Tools`, such as `log_sum_exp`.

Every value is computed by Numerics. Arrays cross between NumPy and .NET in one bulk copy each
way, bit for bit.

## Special functions

The functions take scalars or arrays and broadcast like NumPy ufuncs. An array goes through the
Numerics array overload where there is one, and otherwise element by element through the scalar
member, which is bound once per process.

```python
import numpy as np
from pynumerics import special

print(special.gamma_function(5.0))                # 24.0
print(special.log_gamma(np.array([0.5, 10.0, 200.0])))  # ln|Gamma(x)|, as scipy.special.gammaln
print(special.digamma(1.0))                       # minus the Euler-Mascheroni constant
print(special.erf(0.5), special.inverse_erf(0.5))
print(special.beta_function(2.0, 3.0))            # 1/12
print(special.log_beta(0.5, 1e6))                  # ln B(a, b), as scipy.special.betaln
```

The regularized incomplete gamma functions take `(a, x)` in the order of
`scipy.special.gammainc`, and have inverses:

```python
from pynumerics import special

p = special.lower_incomplete_gamma(3.0, [0.5, 2.0, 8.0])   # P(a, x)
q = special.upper_incomplete_gamma(3.0, 2.0)               # Q(a, x) = 1 - P(a, x)
x = special.inverse_lower_incomplete_gamma(3.0, p)         # back to [0.5, 2.0, 8.0]
```

`incomplete_gamma(a, x)` is the AS 239 algorithm (`Gamma.Incomplete`); the docstrings give the
method and accuracy of each incomplete gamma function. The regularized incomplete beta function
and its inverse take `(a, b, x)` in the order of `scipy.special.betainc`, and `gamma_max()` is the
largest argument for which the gamma function is finite:

```python
from pynumerics import special

p = special.beta_incomplete(2.0, 3.0, [0.2, 0.4, 0.6])     # I_x(a, b)
x = special.beta_incomplete_inverse(2.0, 3.0, p)           # back to [0.2, 0.4, 0.6]
print(special.gamma_max())                                 # about 171.62
```

Bessel functions of integer order broadcast the order against the argument; the docstrings give
the method and accuracy of each.

```python
import numpy as np
from pynumerics import special

print(special.bessel_j0(2.4048255576957727))       # close to zero
print(special.bessel_jn(np.array([[2], [3]]), np.linspace(1.0, 5.0, 5)))  # shape (2, 5)
print(special.bessel_i0e(50.0))
```

Factorials, binomial coefficients and combinations:

```python
from pynumerics import special

print(special.factorial(20))                # 2432902008176640000.0
print(special.binomial_coefficient(52, 5))  # 2598960.0
print(special.find_combinations(2, 4))      # rows [0 1], [0 2], ..., [2 3]
print(special.all_combinations(3))          # 0/1 indicator rows of the 7 non-empty subsets
print(list(special.iter_all_combinations(2)))  # [array([0]), array([1]), array([0, 1])]
advanced, successor = special.next_combination([0, 3], 4)   # True, [1 2]
```

Polynomials, with the constant term first (`polynomial`) or last (`polynomial_reversed`, the
convention of `numpy.polyval`), and Chebyshev series in twice the Chebyshev variable (`chebyshev`,
the Cephes `chbevl` convention):

```python
from pynumerics import special

print(special.polynomial([1.0, 2.0, 3.0], 2.0))     # 1 + 2x + 3x^2 = 17
print(special.polynomial_reversed([1.0, 2.0, 3.0], 2.0))  # x^2 + 2x + 3 = 11
print(special.chebyshev([2.0, 3.0, -1.0, 0.5], 0.6))     # c0/2 + c1 T1(0.3) + c2 T2(0.3) + c3 T3(0.3) = 2.324
```

## Matrices and vectors

`Matrix` and `Vector` hold the .NET object as `.net`. Matrix-valued results stay in .NET as new
wrappers; `numpy.asarray` or `to_array()` copies the values out.

```python
import numpy as np
from pynumerics.linalg import Matrix, Vector

a = Matrix([[4.0, 1.0], [2.0, 3.0]])
b = a.multiply(a.transpose())          # A A^T, a Matrix
print(np.asarray(b))
print(a.determinant(), a.trace())
print(a.multiply([1.0, 1.0]).to_array())  # A x, a Vector
a[0, 1] = 5.0
print(Matrix.identity(3).to_array())
columns = Matrix([[1.0, 2.0], [3.0, 4.0]], by_column=True)   # columns [1, 2] and [3, 4]

v = Vector([3.0, 4.0])
print(v.norm())                        # 5.0
print(Vector.dot_product(v, [1.0, 2.0]))  # 11.0
```

`Matrix.apply` takes a Python callable or C# lambda source, which is compiled with Roslyn and runs
without a Python call per element:

```python
from pynumerics.linalg import Matrix

m = Matrix([[1.0, 4.0], [9.0, 16.0]])
m.apply(lambda value: value + 1.0)     # in place
m.apply("x => x * 2.0")                # C# source, compiled once
m.sqrt()                               # in place
print(m.to_array())
```

## Decompositions

The decompositions take a 2-D array or a `Matrix`, and return factors and solutions as NumPy
arrays.

```python
import numpy as np
from pynumerics.linalg import (CholeskyDecomposition, EigenValueDecomposition, LUDecomposition,
                               QRDecomposition, SingularValueDecomposition)

s = np.array([[4.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]])
b = np.array([1.0, 2.0, 3.0])

chol = CholeskyDecomposition(s)
print(chol.l, chol.solve(b), chol.log_determinant())

lu = LUDecomposition(s)
print(lu.solve(b), lu.solve(np.eye(3)), lu.determinant())

qr = QRDecomposition(np.random.default_rng(1).normal(size=(6, 3)))
print(qr.q, qr.r)

svd = SingularValueDecomposition(np.outer([1.0, 2.0, 3.0], [1.0, 1.0]))
print(svd.w, svd.rank(), svd.nullspace())

eig = EigenValueDecomposition(s)       # symmetric matrices, cyclic Jacobi
print(eig.eigenvalues, eig.eigenvectors, eig.is_converged)
ql = EigenValueDecomposition(s, method="ql")  # Householder + QL, ascending as numpy.linalg.eigh
print(ql.eigenvalues)
```

`CholeskyDecomposition.try_decompose` retries a matrix that is not numerically positive definite
with an escalating diagonal jitter:

```python
from pynumerics.linalg import CholeskyDecomposition

ok, decomposition, jitter = CholeskyDecomposition.try_decompose([[1.0, 1.0], [1.0, 1.0]])
print(ok, jitter)                      # True, a small positive jitter
```

`regularize` floors and caps the eigenvalues of a symmetric matrix, and
`make_symmetric_positive_definite` returns the symmetric part with the smallest diagonal ridge
that the Cholesky decomposition accepts:

```python
from pynumerics import linalg

print(linalg.regularize([[1.0, 0.9], [0.9, 1.0]]))
print(linalg.make_symmetric_positive_definite([[1.0, 0.5], [0.4, 1.0]]))  # [[1, 0.45], [0.45, 1]]
```

## Fourier transforms

Numerics uses the Numerical Recipes conventions: the forward transform is
H(k) = sum h(j) exp(+2 pi i j k / n), equal to `n * numpy.fft.ifft(x)`, and the inverse transform
is unnormalized, equal to `numpy.fft.fft(x)`.

```python
import numpy as np
from pynumerics import mathematics

x = np.arange(8.0) + 1j
print(np.allclose(mathematics.fft(x), 8 * np.fft.ifft(x)))   # True; the length must be a power of 2
print(np.allclose(mathematics.bluestein_fft(x[:6]), 6 * np.fft.ifft(x[:6]))) # any length
print(mathematics.fft_convolve([1.0, 2.0], [1.0, 1.0, 1.0])) # as numpy.convolve: [1, 3, 3, 2]
print(mathematics.fft_autocorrelation(np.cumsum(np.ones(20)), 3)) # columns: lag, autocorrelation
```

## Numerical utilities

```python
import numpy as np
from pynumerics import mathematics

print(mathematics.log_sum_exp([1000.0, 1000.0])) # 1000 + ln 2, without overflow
print(mathematics.log_sum_exp(np.log(0.25), np.log(0.5)))
print(mathematics.log_sum_exp([0.0, -30.0]))     # 9.357622968839737e-14, as scipy.special.logsumexp
print(mathematics.log_sum_exp(np.log([[0.25, 0.5], [0.75, 0.5]]), axis=0))  # about [0, 0], one value per column
print(mathematics.log_sum_exp([]), mathematics.log_sum_exp([1.0, np.inf]))  # -inf inf
print(mathematics.min_max([3.0, 1.0, 2.0]))   # (1.0, 3.0)
print(mathematics.sum_values([1.0, 2.0, 3.0], indicators=[1, 0, 1])) # 4.0
print(mathematics.sequence(0.0, 1.0, 0.25))
print(mathematics.constants.euler, mathematics.constants.golden_ratio)
```

## Compensated sums

`compensated_sum` is Neumaier's compensated summation (`Tools.CompensatedSum`): a running
correction keeps the exact rounding error of each addition, so ill-conditioned sums keep their
digits and a term larger than the running sum is handled, where plain Kahan summation fails.

```python
import math
from pynumerics import mathematics

x = [1e100, 1.0, -1e100]
print(mathematics.sum_values(x), mathematics.compensated_sum(x), math.fsum(x))   # 0.0 1.0 1.0
print(mathematics.compensated_sum([1.0, 1e100, 1.0, -1e100]))                   # 2.0 (Kahan gives 0)
print(mathematics.compensated_sum([1.0, 1e100, 1.0, -1e100], kahan=True))       # 0.0, classic Kahan bit for bit
print(mathematics.compensated_sum([1e100, 1.0, -1e100, 7.0], [1, 1, 1, 0]))     # 1.0, the values marked 1
print(mathematics.compensated_sum([1e100, 1.0, -1e100, 3.0], indices=[3, 0, 1, 2], start=1, count=3))  # 1.0
print(mathematics.compensated_cumulative_sum([1e100, 1.0, -1e100, 0.1]))        # [1e+100 1e+100 1. 1.1]
```

The indicator form selects values as `sum_values` does, and the gather form sums
`values[indices[i]]` for i in [`start`, `start + count`) without copying. NaN and infinities
give the IEEE sum: NaN for a NaN or for both infinities, where `math.fsum` raises for both
infinities, and the infinity for one infinity or an intermediate overflow, where `math.fsum`
raises for the overflow. As for `sum_values`, empty input gives NaN (`math.fsum` gives 0).

`Tools.Log` and `Tools.Log10` are guarded: `mathematics.log(0.0)` is about -36.84 rather than
-inf, because arguments in [0, 1e-16) are raised to 1e-16.
