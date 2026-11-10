# Special Functions

[← Previous: Linear Algebra](linear-algebra.md) | [Back to Index](../index.md) | [Next: ODE Solvers →](ode-solvers.md)

The ***Numerics*** library provides essential special functions commonly used in statistical distributions, numerical analysis, and scientific computing. These functions underpin many of the library's distribution and integration routines; for example, the Normal distribution CDF is computed using the error function, and the Chi-squared CDF uses the incomplete Gamma function.

## Gamma Function

The Gamma function is defined by the integral:

```math
\Gamma(x) = \int_0^{\infty} t^{x-1} e^{-t} \, dt, \quad x > 0
```

It extends the factorial function to real (and complex) numbers: $\Gamma(n) = (n-1)!$ for positive integers. The defining **recurrence relation** is:

```math
\Gamma(x+1) = x \cdot \Gamma(x)
```

Two other important identities are the **reflection formula** $\Gamma(x)\Gamma(1-x) = \frac{\pi}{\sin(\pi x)}$ and the **duplication formula** $\Gamma(x)\Gamma(x + \frac{1}{2}) = \frac{\sqrt{\pi}}{2^{2x-1}}\Gamma(2x)$. The special value $\Gamma(\frac{1}{2}) = \sqrt{\pi}$ follows directly from the reflection formula.

### Basic Gamma Function

```cs
using Numerics.Mathematics.SpecialFunctions;

// Gamma function
double g1 = Gamma.Function(5.0);    // Γ(5) = 4! = 24
double g2 = Gamma.Function(0.5);    // Γ(0.5) = √π ≈ 1.772
double g3 = Gamma.Function(3.5);    // Γ(3.5) ≈ 3.323

Console.WriteLine($"Γ(5) = {g1:F2}");
Console.WriteLine($"Γ(0.5) = {g2:F6}");
Console.WriteLine($"Γ(3.5) = {g3:F3}");

// Verify: Γ(n+1) = n·Γ(n)
double check = 3.5 * Gamma.Function(3.5);
Console.WriteLine($"3.5·Γ(3.5) = {check:F3}");
Console.WriteLine($"Γ(4.5) = {Gamma.Function(4.5):F3}");
```

### Log-Gamma Function

For large arguments, $\Gamma(x)$ overflows double-precision floating point. The log-gamma function $\ln \Gamma(x)$ grows much more slowly and is used internally throughout the library for distribution computations. For large $x$, Stirling's approximation gives:

```math
\ln \Gamma(x) \approx x \ln x - x + \frac{1}{2}\ln\!\left(\frac{2\pi}{x}\right)
```

Use log-gamma to avoid overflow:

```cs
// Regular gamma would overflow for large x
double x = 200.0;

// Use log-gamma
double logGamma = Gamma.LogGamma(x);
Console.WriteLine($"ln(Γ(200)) = {logGamma:F2}");

// Gamma itself would be huge: exp(logGamma)
// Don't compute directly - work in log space
```

`LogGamma` evaluates the rational minimax approximations of Cody and Hillstrom [[10]](#10), with the coefficients of the SPECFUN routine DLGAMA [[11]](#11), on $(0, 1.5]$, $(1.5, 4]$ and $(4, 12]$, and Stirling's formula with a minimax correction above 12. The approximations about 1 and 2 carry the factor $(x-1)$ or $(x-2)$, so the relative error stays below $10^{-15}$ next to the zeros $\ln\Gamma(1) = \ln\Gamma(2) = 0$ as well as across $(0, 10^{300})$.

For negative non-integer arguments `LogGamma` returns $\ln|\Gamma(x)|$ from the reflection formula $\ln|\Gamma(x)| = \ln\pi - \ln|\sin(\pi x)| - \ln\Gamma(1-x)$, as scipy `gammaln` does (for example `LogGamma(-0.5)` ≈ 1.26551212348465). Zero and the negative integers are poles and return positive infinity.

### Digamma and Trigamma

The digamma function $\psi(x) = \frac{d}{dx}\ln\Gamma(x) = \frac{\Gamma'(x)}{\Gamma(x)}$ and the trigamma function $\psi'(x) = \frac{d^2}{dx^2}\ln\Gamma(x)$ arise in maximum likelihood estimation for the Gamma, Beta, and Dirichlet distributions. The digamma function satisfies $\psi(x+1) = \psi(x) + \frac{1}{x}$ and for large $x$, $\psi(x) \approx \ln x - \frac{1}{2x}$.

```cs
// Digamma: ψ(x) = d/dx[ln(Γ(x))] = Γ'(x)/Γ(x)
double digamma = Gamma.Digamma(2.0);

// Trigamma: ψ'(x) = d²/dx²[ln(Γ(x))]
double trigamma = Gamma.Trigamma(2.0);

Console.WriteLine($"ψ(2) = {digamma:F6}");
Console.WriteLine($"ψ'(2) = {trigamma:F6}");

// Applications: moment matching in Gamma distribution MLE
```

`Trigamma` shifts its argument with $\psi'(x) = \psi'(x+1) + 1/x^2$ until it reaches 20 and then sums the asymptotic series with exact Bernoulli coefficients, so it is accurate to a few ulps for every $x > 0$ (including tiny $x$, where $\psi'(x) \approx 1/x^2 + \pi^2/6$). Non-positive arguments throw `ArgumentOutOfRangeException`.

### Incomplete Gamma Functions

The **lower incomplete gamma function** is defined as the integral with a finite upper limit:

```math
\gamma(a, x) = \int_0^{x} t^{a-1} e^{-t} \, dt
```

The library returns the **regularized** forms $P(a,x) = \gamma(a,x)/\Gamma(a)$ and $Q(a,x) = 1 - P(a,x)$, which range from 0 to 1. These are equivalent to the CDF and survival function of the Gamma distribution: the Chi-squared CDF with $k$ degrees of freedom is $P(k/2, x/2)$.

```cs
// Regularized lower incomplete gamma: P(a,x) = γ(a,x) / Γ(a)
double P = Gamma.LowerIncomplete(a: 2.0, x: 3.0);

// Regularized upper incomplete gamma: Q(a,x) = Γ(a,x) / Γ(a)
double Q = Gamma.UpperIncomplete(a: 2.0, x: 3.0);

// Verify: P(a,x) + Q(a,x) = 1
double sum = P + Q;

Console.WriteLine($"P(2, 3): {P:F6}");
Console.WriteLine($"Q(2, 3): {Q:F6}");
Console.WriteLine($"P + Q:   {sum:F6}");  // 1.000000
```

### Alternative: Gamma.Incomplete

An alternative method with different parameter ordering:

```cs
// Regularized incomplete gamma (CDF of Gamma distribution)
double P = Gamma.Incomplete(X: 3.0, alpha: 2.0);

Console.WriteLine($"P(2, 3) = {P:F6}");
Console.WriteLine("This equals the Gamma(2,1) CDF at x=3");

// Inverse: find x such that P(a,x) = p
double xInv = Gamma.InverseLowerIncomplete(a: 2.0, y: 0.9);
Console.WriteLine($"P(2, {xInv:F3}) = 0.9");
```

`InverseLowerIncomplete` and `InverseUpperIncomplete` invert the smaller tail directly (a probability above 0.5 is replaced by its exact complement and the other tail), so small probabilities in either tail keep full relative precision (for example `InverseLowerIncomplete(2.0, 1e-12)` returns 1.4142142290401938E-06, the correctly rounded root, where scipy `gammaincinv` is 4 ulps away). The root is found by a bracketed Halley iteration from a Wilson-Hilferty, small-$x$ series or asymptotic upper-tail start, so the inverses return a value for every $a > 0$ and $y \in [0, 1]$: `InverseLowerIncomplete` returns 0 at $y = 0$ and where the quantile is below the smallest positive double, $+\infty$ at $y = 1$, and NaN for invalid arguments. Edge values of `LowerIncomplete` and `UpperIncomplete` follow scipy `gammainc` and `gammaincc`: $P(0,x) = 1$ and $Q(0,x) = 0$ for $x > 0$, and a negative $a$ (or $a = x = 0$) returns NaN. For $x \le 0$ with $a > 0$ the library returns $P = 0$ and $Q = 1$ (the CDF below the support), where scipy returns NaN for negative $x$. An infinite shape gives $P = 0$ and $Q = 1$ for finite $x$. `Gamma.Incomplete(X, alpha)` keeps the argument order and validation of algorithm AS 239 (it throws for $\alpha \le 0$ or $X < 0$ and returns NaN for NaN arguments) and returns the same value as `LowerIncomplete(alpha, X)`.

`LowerIncomplete` and `UpperIncomplete` follow the organization of Cephes `igam` and `igamc` [[4]](#4), [[5]](#5) and always evaluate the smaller of $P$ and $Q$ directly. For $a > 20$ with $x$ close to $a$ ($|x - a|/a$ below 0.3 for $a < 200$, below $4.5/\sqrt{a}$ above, and also below 0.1 for $a \ge 10^4$, where the series would need of order $\sqrt{a}$ terms) they sum Temme's uniform asymptotic expansion [[6]](#6), [[7]](#7), [[3]](#3) (DLMF 8.12); elsewhere they use the power series of $P$, the continued fraction of $Q$, or the small-$x$ series of $Q$, with the factor $x^a e^{-x}/\Gamma(a)$ formed without cancellation. Against mpmath at 50 digits over $a \in [10^{-3}, 10^6]$ and $x$ from $10^{-4}a$ to $30a$ the relative error is at most about $2.4 \times 10^{-14}$ where the result exceeds $10^{-30}$; deeper in the tails it grows with $|\ln P|$ (about $7 |\ln P|$ units of roundoff at worst), the conditioning of the exponential factor. Results underflow only below the double range, not at the $e^{-180}$ floor of the AS 239 algorithm: `Gamma.Incomplete(0.004641588833612782, 46.41588833612782)` returns 1.8173331985826976E-167 (mpmath: 1.8173331985826257E-167). Against the 375 mpmath roots of the oracle grid (50 digits; $a$ from 0.01 to $10^5$, $y$ from $10^{-300}$ to $1 - 10^{-10}$, both tails) the relative error of the inverses is at most $1.8 \times 10^{-15}$ for $a \ge 1$ and $2.8 \times 10^{-14}$ for $0.01 \le a < 1$, and on the unit-test table at $a = 10^{-3}$ it is at most $2.5 \times 10^{-13}$; the condition number of the quantile is near $1/a$ for small quantiles. The oracle and its generators are in `docs/mathematics/oracles/special-functions`.

## Beta Function

The Beta function is defined by the integral:

```math
B(a,b) = \int_0^1 t^{a-1}(1-t)^{b-1}\,dt = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
```

The relationship to the Gamma function makes computation straightforward via $\ln B(a,b) = \ln\Gamma(a) + \ln\Gamma(b) - \ln\Gamma(a+b)$.

### Basic Beta Function

```cs
using Numerics.Mathematics.SpecialFunctions;

// Beta function
double beta = Beta.Function(a: 2.0, b: 3.0);

// Verify relation to Gamma
double gamma_a = Gamma.Function(2.0);
double gamma_b = Gamma.Function(3.0);
double gamma_ab = Gamma.Function(5.0);
double betaCheck = gamma_a * gamma_b / gamma_ab;

Console.WriteLine($"B(2,3) = {beta:F6}");
Console.WriteLine($"Γ(2)Γ(3)/Γ(5) = {betaCheck:F6}");
```

### Log-Beta Function

For large shapes $B(a,b)$ underflows, and forming $\ln\Gamma(a) + \ln\Gamma(b) - \ln\Gamma(a+b)$ subtracts large, nearly equal numbers: at $a = 1/2$, $b = 10^6$ the terms are near $1.3 \times 10^7$ and the difference is in error by about $10^{-9}$; at $b = 10^{15}$ it keeps no correct digit. `Beta.LogFunction` evaluates $\ln B(a,b)$ with the routines BETALN, ALGDIV and BCORR of algorithm 708 [[12]](#12), which combine the Stirling terms analytically and evaluate the difference of the Stirling corrections directly. Below 8 an argument above 2 is reduced to $(1, 2]$ by $\Gamma(x+1) = x\Gamma(x)$ before $\ln\Gamma$ is taken. `Beta.Function` is its exponential.

The unit tests check 20 random points per region against mpmath (400 digits), with $a = \min$ and $b = \max$ of the shapes from $10^{-10}$ to $10^{300}$; the error relative to $\max(1, |\ln B|)$ is at most:

| Region | Error bound |
|---|---|
| $a < 1$, $b \le 2$ | $2.8 \times 10^{-16}$ |
| $a < 1$, $2 < b < 8$ | $6.7 \times 10^{-16}$ |
| $1 \le a \le 2$, $b \le 2$ | $4.5 \times 10^{-16}$ |
| $1 \le a < 8$, $2 < b < 8$ | $4.1 \times 10^{-16}$ |
| $a < 8$, $b \ge 8$ | $1.2 \times 10^{-15}$ |
| $a \ge 8$ | $2.6 \times 10^{-16}$ |

Where $\ln B$ is close to zero the error is absolute, not relative: for example $\ln B(0.547, 2.642) \approx 1.84 \times 10^{-4}$ carries an error near $8 \times 10^{-17}$.

```cs
double logBeta = Beta.LogFunction(0.5, 1E6);   // -6.335390211057437
double logBetaEqual = Beta.LogFunction(1000, 1000);   // -1388.482601635902
```

### Incomplete Beta Function

The incomplete Beta function is defined as:

```math
B_x(a,b) = \int_0^x t^{a-1}(1-t)^{b-1}\,dt
```

The library returns the **regularized** form $I_x(a,b) = B_x(a,b) / B(a,b)$, which ranges from 0 to 1 and is the CDF of the Beta distribution. It is also used internally for the Student's t-test, the F-distribution CDF, and the binomial distribution CDF.

```cs
// Regularized incomplete beta: Iₓ(a,b) = Bₓ(a,b) / B(a,b)
double Ix = Beta.Incomplete(a: 2.0, b: 3.0, x: 0.4);

Console.WriteLine($"Iₓ(2,3,0.4) = {Ix:F6}");
Console.WriteLine("This is the CDF of Beta(2,3) at x=0.4");

// To recover the raw (non-regularized) incomplete beta:
double rawIncomplete = Ix * Beta.Function(2.0, 3.0);
Console.WriteLine($"Bₓ(2,3,0.4) = {rawIncomplete:F6}");
```

### Inverse Incomplete Beta

Quantile function for Beta distribution:

```cs
// Find x such that I(a,b,x) = p
double x = Beta.IncompleteInverse(aa: 2.0, bb: 3.0, yy0: 0.5);

Console.WriteLine($"50th percentile of Beta(2,3): {x:F4}");

// Verify (Beta.Incomplete is already regularized)
double check = Beta.Incomplete(2.0, 3.0, x);
Console.WriteLine($"Verification: I(2,3,{x:F4}) = {check:F6}");  // 0.500000
```

A small tail probability keeps its relative accuracy: where the search ends on the complement of a small tail, the root is refined by Newton's method on the logarithm of that tail. A root below $e^{-708}$ is the leading series term $(y\,a\,B(a,b))^{1/a}$, which returns a subnormal number, or 0 below the double range, instead of throwing. The tests check eleven mpmath roots with targets from $10^{-300}$ to $1 - 10^{-12}$ and shapes from 0.1 to $10^4$ to 1E-15 relative.

## Error Function

The error function is defined as the integral:

```math
\text{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2}\,dt
```

It is closely related to the Normal distribution CDF: $\Phi(x) = \frac{1}{2}\,\text{erfc}\!\left(-\frac{x}{\sqrt{2}}\right)$. The complementary error function $\text{erfc}(x) = 1 - \text{erf}(x)$ is computed directly, never by subtraction, as $\text{erfc}(x) = e^{-x^2}\,\text{erfcx}(x)$ for $x \ge 0$ and $2 - \text{erfc}(-x)$ below. The scaled function $\text{erfcx}(x) = e^{x^2}\text{erfc}(x)$ is the Chebyshev expansion of Schonfelder [[8]](#8) for $x < 3$ and, from 3 on, $g(1/x^2)/(x\sqrt{\pi})$ with $g$ a degree-21 polynomial fitted with mpmath (`docs/mathematics/oracles/special-functions/erfc_coefficients.py`, which also recomputes Schonfelder's coefficients). The factor $e^{-x^2}$ is formed with the splitting of Cody [[9]](#9): $x = h + l$ with $h$ truncated to 20 fraction bits, so $h^2$ is exact and $e^{-x^2} = e^{-h^2} e^{-(x-h)(x+h)}$; a rounded $x^2$ would cost about $x^2$ units in the last place. `Normal.StandardCDF` uses the same pieces, with $e^{-z^2/2}$ formed from $z$ itself, so the scaling by $\sqrt{2}$ never enters the exponential. Against mpmath the relative error of `Erf.Erfc` is at most $5 \times 10^{-16}$ (root mean square about $2 \times 10^{-16}$) for $x$ up to 26.5, where the result becomes subnormal, and that of the lower tail of `Normal.StandardCDF` is at most $6.3 \times 10^{-16}$ down to $z = -37.5$, the accuracy of R's `pnorm`; scipy's `erfc` errs by up to $5.7 \times 10^{-14}$ near $x = 26$.

### Error Function and Complement

```cs
using Numerics.Mathematics.SpecialFunctions;

// Error function: erf(x) = (2/√π) ∫₀ˣ e^(-t²) dt
double erf = Erf.Function(1.0);

// Complementary error function: erfc(x) = 1 - erf(x)
double erfc = Erf.Erfc(1.0);

Console.WriteLine($"erf(1) = {erf:F6}");
Console.WriteLine($"erfc(1) = {erfc:F6}");
Console.WriteLine($"Sum = {erf + erfc:F6}");  // Should be 1.0

// Relation to normal distribution
// Φ(x) = 0.5[1 + erf(x/√2)]
double x = 1.0;
double phi = 0.5 * (1 + Erf.Function(x / Math.Sqrt(2)));
Console.WriteLine($"Φ(1) via erf = {phi:F6}");
```

### Inverse Error Function

```cs
// Inverse erf: find x such that erf(x) = y
double y = 0.5;
double xInv = Erf.InverseErf(y);

Console.WriteLine($"erf⁻¹(0.5) = {xInv:F6}");
Console.WriteLine($"Verification: erf({xInv:F6}) = {Erf.Function(xInv):F6}");

// Inverse erfc
double xInvC = Erf.InverseErfc(0.5);
Console.WriteLine($"erfc⁻¹(0.5) = {xInvC:F6}");
```

`InverseErf` is odd in $y$ and keeps full relative precision across its range: for $|y| < 0.5$ a normal-quantile start is refined with two Halley steps on $\text{erf}(x) - |y|$, and for $|y| \ge 0.5$ it evaluates `InverseErfc(1 - |y|)`, where $1 - |y|$ is exact. `InverseErf(±1)` returns $\pm\infty$. `Erf.Function` returns $2x/\sqrt{\pi}$ for $|x| < 10^{-8}$, so tiny arguments such as $10^{-200}$ no longer underflow to zero.

## Factorial and Combinatorics

### Factorial

```cs
using Numerics.Mathematics.SpecialFunctions;

// Integer factorial
double fact5 = Factorial.Function(5);      // 5! = 120
double fact10 = Factorial.Function(10);    // 10! = 3,628,800

Console.WriteLine($"5! = {fact5:F0}");
Console.WriteLine($"10! = {fact10:N0}");

// Log factorial for large numbers
double logFact100 = Factorial.LogFactorial(100);
Console.WriteLine($"ln(100!) = {logFact100:F2}");

// Too large for double: use log space
double fact100 = Math.Exp(logFact100);
Console.WriteLine($"100! ≈ {fact100:E2}");
```

### Binomial Coefficients

```cs
// Binomial coefficient: C(n,k) = n!/(k!(n-k)!)
double c_10_3 = Factorial.BinomialCoefficient(n: 10, k: 3);

Console.WriteLine($"C(10,3) = {c_10_3:F0}");  // 120

// Verify: 10!/(3!·7!) = 3,628,800/(6·5040) = 120

// Pascal's triangle relation: C(n,k) = C(n-1,k-1) + C(n-1,k)
double check = Factorial.BinomialCoefficient(9, 2) + Factorial.BinomialCoefficient(9, 3);
Console.WriteLine($"C(9,2) + C(9,3) = {check:F0}");
Console.WriteLine($"C(10,3) = {c_10_3:F0}");
```

`BinomialCoefficient` uses the multiplicative form $r \leftarrow r\,(n-k+i)/i$ with $k \leftarrow \min(k, n-k)$, in exact 64-bit integer arithmetic while the result fits, so every result below $2^{53}$ is exact and every result below about $9.2 \times 10^{18}$ is the correctly rounded double. Larger results continue in floating point with a relative error of a few ulps per step and return $+\infty$ on overflow; invalid arguments ($n < 0$, $k < 0$ or $k > n$) return 0.

```cs
double c = Factorial.BinomialCoefficient(50, 25);
Console.WriteLine($"C(50,25) = {c:F0}");  // 126410606437752 (exact)
```

`Factorial.Function(long)` returns the same tabulated values as the `int` overload, so it is exact through 22!, and it returns $+\infty$ above 170! where the result exceeds the double range.

### Combinations

```cs
// Generate all combinations of m items from n
int m = 3;  // Choose 3
int n = 5;  // From 5

var combinations = Factorial.FindCombinations(m, n);

Console.WriteLine($"All ways to choose {m} items from {n}:");
foreach (var combo in combinations)
{
    Console.WriteLine($"  [{string.Join(", ", combo)}]");
}

// Total count should equal C(5,3) = 10
int count = combinations.Count();
Console.WriteLine($"Total: {count} combinations");
```

## Bessel Functions

Bessel functions arise in problems with cylindrical symmetry and in directional statistics (e.g., the Von Mises distribution). The library provides the Bessel functions of the first and second kind and the modified Bessel functions of the first and second kind.

### Modified Bessel Functions of the First Kind

```cs
using Numerics.Mathematics.SpecialFunctions;

// I₀(x) - Modified Bessel function of the first kind, order 0
double i0 = Bessel.I0(2.5);
Console.WriteLine($"I₀(2.5) = {i0:F6}");  // ≈ 3.289839

// I₁(x) - Modified Bessel function of the first kind, order 1
double i1 = Bessel.I1(2.5);
Console.WriteLine($"I₁(2.5) = {i1:F6}");  // ≈ 2.516716

// Iₙ(x) - Modified Bessel function of the first kind, integer order
double i_n = Bessel.In(2, 3.0);  // I₂(3.0)
Console.WriteLine($"I₂(3.0) = {i_n:F6}");
```

`I0` and `I1` are the Cephes expansions of `I0e` and `I1e` multiplied by $e^{|x|}$ [[4]](#4), [[5]](#5), accurate to about $10^{-15}$ relative, and split the exponential above 709, so they stay finite until the function itself exceeds the double range (`I0(711)` ≈ 9.087E306). `In` starts Miller's downward recurrence above both $n$ and $|x|$ and normalizes by `I0`.

### Exponentially Scaled Bessel Functions

`Bessel.I0e` and `Bessel.I1e` return the exponentially scaled functions

```math
I_{0e}(x) = e^{-|x|} I_0(x), \qquad I_{1e}(x) = e^{-|x|} I_1(x)
```

They sum the Chebyshev expansions of the Cephes library [[4]](#4), [[5]](#5) on $[0, 8]$ and $(8, \infty)$, are accurate to double precision, and are finite for every $x$. Prefer them to `I0` and `I1` whenever full precision matters or the argument is large: ratios such as $I_1(x)/I_0(x) = I_{1e}(x)/I_{0e}(x)$ and logarithms $\ln I_0(x) = \ln I_{0e}(x) + |x|$ are formed without overflow. The von Mises distribution uses them for its normalizing constant and mean resultant length.

```cs
// I₀e(x) = exp(−|x|) I₀(x) and I₁e(x) = exp(−|x|) I₁(x)
double i0e = Bessel.I0e(2.5);
double i1e = Bessel.I1e(2.5);
Console.WriteLine($"I₀e(2.5) = {i0e:F12}");  // ≈ 0.270046441612
Console.WriteLine($"I₁e(2.5) = {i1e:F12}");  // ≈ 0.206584649531

// Finite where I₀ and I₁ overflow
Console.WriteLine($"I₀e(1000) = {Bessel.I0e(1000.0):E12}");
```

### Log-Space Bessel Computations

For large arguments where the Bessel function overflows, compute the logarithm from the scaled function:

```cs
// ln I₀(x) = ln I₀e(x) + |x|, finite for any argument
double x = 5000.0;
double logI0 = Math.Log(Bessel.I0e(x)) + Math.Abs(x);
Console.WriteLine($"ln(I₀({x})) = {logI0:F4}");

// ln I₁(x) = ln I₁e(x) + |x| for x > 0
double logI1 = Math.Log(Bessel.I1e(x)) + x;
Console.WriteLine($"ln(I₁({x})) = {logI1:F4}");
```

### Bessel Function Ratios

Ratios of Bessel functions appear in maximum likelihood estimation for the Von Mises distribution:

```cs
// I₁(x)/I₀(x) ratio - used in Von Mises MLE; the scale factors cancel
double ratio = Bessel.I1e(5.0) / Bessel.I0e(5.0);
Console.WriteLine($"I₁(5)/I₀(5) = {ratio:F6}");

// For large kappa, this ratio approaches 1
double ratioLarge = Bessel.I1e(100.0) / Bessel.I0e(100.0);
Console.WriteLine($"I₁(100)/I₀(100) = {ratioLarge:F8}");
```

### Modified Bessel Functions of the Second Kind

```cs
// K₀(x) - Modified Bessel function of the second kind, order 0
double k0 = Bessel.K0(1.0);
Console.WriteLine($"K₀(1.0) = {k0:F6}");  // ≈ 0.421024

// K₁(x) - Order 1
double k1 = Bessel.K1(1.0);
Console.WriteLine($"K₁(1.0) = {k1:F6}");  // ≈ 0.601907

// Kₙ(x) - Integer order
double k_n = Bessel.Kn(2, 1.5);
Console.WriteLine($"K₂(1.5) = {k_n:F6}");
```

`K0` and `K1` follow Cephes `k0` and `k1` [[4]](#4): a Chebyshev expansion in $x^2$ plus the logarithmic term for $x \le 2$, and an expansion of $e^x \sqrt{x} K(x)$ in $8/x$ above, with the exponential split above 700 so the result underflows only where $K$ itself does. `Kn` uses the forward recurrence, which is stable for $K$.

### Bessel Functions of the First and Second Kind

```cs
double j0 = Bessel.J0(2.5);      // ≈ -0.048384
double y1 = Bessel.Y1(2.5);      // ≈  0.145918
double j5 = Bessel.Jn(5, 2.5);   // ≈  0.019502
Console.WriteLine($"J₀(2.5) = {j0:F6}, Y₁(2.5) = {y1:F6}, J₅(2.5) = {j5:F6}");
```

For $|x| \le 8$, `J0`, `J1`, `Y0` and `Y1` sum Chebyshev expansions from which the zeros in that range are factored out, as Cephes does for `j0` [[4]](#4), so the relative error stays near $10^{-15}$ even within $10^{-6}$ of a zero (scipy, like Cephes, is accurate there only in the absolute sense). Above 8 they use Hankel's $P$ and $Q$ [[1]](#1) (9.2.5 to 9.2.10), from Chebyshev expansions in $8/x$ up to 25 and from the asymptotic series above, with errors near $10^{-15}$ of the envelope $\sqrt{2/(\pi x)}$. `Yn` uses the forward recurrence; `Jn` uses the forward recurrence for $|x| > n$ and, for $|x| \le n$, the continued fraction for $J_n/J_{n-1}$ followed by the backward recurrence normalized by `J0` or `J1` [[2]](#2). The expansion coefficients were computed with mpmath by `docs/mathematics/oracles/special-functions/bessel_chebyshev_coefficients.py`.

| Function | Worst error against mpmath on the oracle grid |
|---|---|
| `I0`, `I1` | 8.2E-16 relative |
| `K0`, `K1` | 8.5E-16 relative |
| `In`, `Kn` (orders 2 to 50) | 5.2E-15 relative |
| `J0`, `J1`, `Y0`, `Y1` | 7.8E-16 of the envelope (1.1E-15 relative for $x \le 8$) |
| `Jn`, `Yn` (orders 2 to 50) | 6.1E-15 of the envelope |
| `Debye.Function` (order 3) | 4.4E-16 relative |
| `LowerIncomplete`, `UpperIncomplete` | 2.4E-14 relative above 1E-30 |
| `InverseLowerIncomplete`, `InverseUpperIncomplete` | 1.8E-15 relative for a ≥ 1, 2.8E-14 for 0.01 ≤ a < 1 (2.5E-13 at a = 1E-3, on the unit-test table) |

## Debye Functions

`Debye.Function` is the order-3 Debye function $D_3(x) = \frac{3}{x^3}\int_0^x \frac{t^3}{e^t - 1}\,dt$ and `Debye.FunctionOrderOne` the order-1 function $D_1$. $D_3$ sums its Maclaurin series in Bernoulli numbers for $x \le 3$ and, above, the exponential series written against the limit $\pi^4/15$ [[1]](#1) (27.1), accurate to about $10^{-15}$ relative.

```cs
double d3 = Debye.Function(1.0);          // ≈ 0.674416
double d1 = Debye.FunctionOrderOne(1.0);  // ≈ 0.777505
Console.WriteLine($"D₃(1) = {d3:F6}, D₁(1) = {d1:F6}");
```

## Polynomial and Chebyshev Evaluation

`Evaluate.Polynomial` sums $\sum c_i x^i$ with the constant term first, and `Evaluate.PolynomialRev` the same polynomial with the constant term last. `Evaluate.Chebyshev` sums the Chebyshev series $c_0/2 + \sum_{k \ge 1} c_k T_k(x/2)$ by the Clenshaw recurrence (Clenshaw, Mathematical Tables and Other Aids to Computation 9, 1955), with the argument already mapped to $[-2, 2]$ as in the Cephes routine `chbevl`. An optional term count truncates the series after its leading coefficients. The Bessel, scaled complementary error and `MultivariateNormal.MVNPHI` approximations are summed by this routine.

```cs
double[] c = { 2.0, 3.0, -1.0, 0.5 };          // c0, c1, c2, c3
double p = Evaluate.Polynomial(c, 0.3);        // 2 + 3x - x^2 + 0.5x^3
double s = Evaluate.Chebyshev(c, 2 * 0.3);     // c0/2 + c1 T1(0.3) + c2 T2(0.3) + c3 T3(0.3) = 2.324
double s2 = Evaluate.Chebyshev(c, 0.6, 2);     // the first two terms: 1 + 0.9 = 1.9
```

## Practical Applications

### Example 1: Gamma Distribution Moments

```cs
// Gamma distribution with shape α and rate β has:
// Mean = α/β
// Variance = α/β²

double alpha = 5.0;
double beta = 2.0;

// Can be computed using Gamma function
double mean = alpha / beta;
double variance = alpha / (beta * beta);

Console.WriteLine($"Gamma({alpha}, {beta}) distribution:");
Console.WriteLine($"  Mean = {mean:F2}");
Console.WriteLine($"  Variance = {variance:F2}");
Console.WriteLine($"  Std Dev = {Math.Sqrt(variance):F2}");

// Factorial moment: E[X(X-1)...(X-k+1)] = Γ(α+k)/(βᵏΓ(α))
int k = 2;
double factMoment = Gamma.Function(alpha + k) / (Math.Pow(beta, k) * Gamma.Function(alpha));
Console.WriteLine($"  E[X(X-1)] = {factMoment:F2}");
```

### Example 2: Normal Distribution CDF

```cs
// Compute Normal(0,1) CDF using error function
double x = 1.5;

// Φ(x) = 0.5[1 + erf(x/√2)]
double cdf = 0.5 * (1.0 + Erf.Function(x / Math.Sqrt(2.0)));

Console.WriteLine($"Φ({x}) = {cdf:F6}");
Console.WriteLine("Compare with Normal distribution class for verification");

// Tail probability
double tail = 1.0 - cdf;
Console.WriteLine($"P(X > {x}) = {tail:F6}");

// Using erfc for better precision in tails
double tailAlt = 0.5 * Erf.Erfc(x / Math.Sqrt(2.0));
Console.WriteLine($"P(X > {x}) via erfc = {tailAlt:F6}");
```

### Example 3: Chi-Squared CDF

```cs
// Chi-squared distribution with k degrees of freedom
// CDF = P(k/2, x/2) where P is regularized lower incomplete gamma

int k = 5;  // Degrees of freedom
double x = 8.0;

// CDF at x
double cdf = Gamma.Incomplete(X: x / 2.0, alpha: k / 2.0);

Console.WriteLine($"Chi-squared({k}) CDF at {x}:");
Console.WriteLine($"  P(X ≤ {x}) = {cdf:F6}");

// Critical value for α=0.05
double alpha_level = 0.05;
double critical = 2.0 * Gamma.InverseLowerIncomplete(a: k / 2.0, y: 1.0 - alpha_level);
Console.WriteLine($"  95th percentile = {critical:F3}");
```

### Example 4: Beta Distribution

```cs
// Beta(a,b) distribution CDF using incomplete beta
double a = 2.0;
double b = 3.0;
double x = 0.4;

// CDF = I(x;a,b) = Bₓ(a,b)/B(a,b); Beta.Incomplete is already regularized
double cdf = Beta.Incomplete(a, b, x);

Console.WriteLine($"Beta({a},{b}) CDF at x={x}:");
Console.WriteLine($"  P(X ≤ {x}) = {cdf:F6}");

// Median
double median = Beta.IncompleteInverse(a, b, 0.5);
Console.WriteLine($"  Median = {median:F4}");

// 90% confidence interval
double lower = Beta.IncompleteInverse(a, b, 0.05);
double upper = Beta.IncompleteInverse(a, b, 0.95);
Console.WriteLine($"  90% CI: [{lower:F4}, {upper:F4}]");
```

### Example 5: Stirling's Approximation

```cs
// Stirling's approximation for large factorials
// ln(n!) ≈ n·ln(n) - n + 0.5·ln(2πn)

int n = 50;

double exactLogFact = Factorial.LogFactorial(n);
double stirling = n * Math.Log(n) - n + 0.5 * Math.Log(2 * Math.PI * n);

Console.WriteLine($"ln({n}!):");
Console.WriteLine($"  Exact = {exactLogFact:F6}");
Console.WriteLine($"  Stirling = {stirling:F6}");
Console.WriteLine($"  Error = {Math.Abs(exactLogFact - stirling):F6}");

// Stirling is very accurate for large n
double relativeError = Math.Abs(exactLogFact - stirling) / exactLogFact;
Console.WriteLine($"  Relative error = {relativeError:P4}");
```

## Function Summary

| Function | Purpose | Key Methods |
|----------|---------|-------------|
| **Gamma** | Factorial extension | `Function()`, `LogGamma()`, `Digamma()` |
| **Incomplete Gamma** | Chi-squared, Gamma CDF | `LowerIncomplete()`, `UpperIncomplete()` |
| **Beta** | Beta distribution | `Function()`, `LogFunction()`, `Incomplete()` |
| **Error** | Normal distribution | `Function()`, `Erfc()`, `InverseErf()`, `InverseErfc()` |
| **Factorial** | Combinatorics | `Function()`, `BinomialCoefficient()` |
| **Bessel** | Cylindrical, directional stats | `I0()`, `I1()`, `I0e()`, `I1e()`, `In()`, `K0()`, `K1()`, `Kn()`, `J0()`, `J1()`, `Jn()`, `Y0()`, `Y1()`, `Yn()` |
| **Debye** | Heat capacity, thermal models | `Function()` (order 3), `FunctionOrderOne()` |

### List Overloads

Every scalar member of `Gamma`, `Erf`, `Bessel`, `Debye` and `Beta` also has a list form that evaluates a whole `IList<double>` in one call and returns a `double[]`. Functions of several arguments take one list per argument, of equal length (an `ArgumentException` otherwise); the integer-order Bessel functions take the order once and a list of arguments. Each value is the scalar result bit for bit, and an exception the scalar member throws for a value propagates unchanged. A caller crossing a language boundary, such as Python through pythonnet, pays the crossing once instead of once per value.

```cs
double[] p = Gamma.LowerIncomplete(new[] { 0.5, 2.0, 10.0 }, new[] { 1.0, 2.0, 12.0 });
double[] j = Bessel.Jn(3, new[] { 0.5, 5.0, 50.0 });
double[] d = Debye.Function(new List<double> { 0.1, 1.0, 10.0 });
```

## Implementation Notes

- All functions use high-precision polynomial or rational approximations
- Log-space variants prevent overflow for large arguments
- Inverse functions use Newton-Raphson iteration with appropriate starting values
- Special care for edge cases: $\Gamma(x)$ near poles at non-positive integers, $\text{erf}(x)$ for large $|x|$, incomplete gamma/beta for extreme parameter ratios

---

## References

<a id="1">[1]</a> M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions*, New York: Dover Publications, 1964.

<a id="2">[2]</a> W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery, *Numerical Recipes: The Art of Scientific Computing*, 3rd ed., Cambridge, UK: Cambridge University Press, 2007.

<a id="3">[3]</a> NIST Digital Library of Mathematical Functions, https://dlmf.nist.gov/.

<a id="4">[4]</a> S. L. Moshier, *Methods and Programs for Mathematical Functions*, Chichester, UK: Ellis Horwood, 1989.

<a id="5">[5]</a> Cephes Mathematical Library, https://www.netlib.org/cephes/.

<a id="6">[6]</a> N. M. Temme, "The asymptotic expansion of the incomplete gamma functions," *SIAM Journal on Mathematical Analysis*, vol. 10, no. 4, pp. 757-766, 1979.

<a id="7">[7]</a> A. R. DiDonato and A. H. Morris, "Computation of the incomplete gamma function ratios and their inverse," *ACM Transactions on Mathematical Software*, vol. 12, no. 4, pp. 377-393, 1986.

<a id="8">[8]</a> J. L. Schonfelder, "Chebyshev expansions for the error and related functions," *Mathematics of Computation*, vol. 32, no. 144, pp. 1232-1240, 1978.

<a id="9">[9]</a> W. J. Cody, "Rational Chebyshev approximations for the error function," *Mathematics of Computation*, vol. 23, no. 107, pp. 631-637, 1969.

<a id="10">[10]</a> W. J. Cody and K. E. Hillstrom, "Chebyshev approximations for the natural logarithm of the Gamma function," *Mathematics of Computation*, vol. 21, no. 98, pp. 198-203, 1967.

<a id="11">[11]</a> W. J. Cody, "Algorithm 715: SPECFUN, a portable FORTRAN package of special function routines and test drivers," *ACM Transactions on Mathematical Software*, vol. 19, no. 1, pp. 22-30, 1993.

<a id="12">[12]</a> A. R. DiDonato and A. H. Morris, "Algorithm 708: Significant digit computation of the incomplete beta function ratios," *ACM Transactions on Mathematical Software*, vol. 18, no. 3, pp. 360-373, 1992.

---

[← Previous: Linear Algebra](linear-algebra.md) | [Back to Index](../index.md) | [Next: ODE Solvers →](ode-solvers.md)
