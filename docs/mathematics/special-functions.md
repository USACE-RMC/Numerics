# Special Functions

[← Previous: Linear Algebra](linear-algebra.md) | [Back to Index](../index.md) | [Next: ODE Solvers →](ode-solvers.md)

The ***Numerics*** library provides essential special functions commonly used in statistical distributions, numerical analysis, and scientific computing. These include Gamma, Beta, Error functions, and combinatorial functions.

## Gamma Function

The Gamma function Γ(x) extends the factorial function to real and complex numbers: Γ(n) = (n-1)! for positive integers.

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
Console.WriteLine($"4.5·Γ(3.5) = {check:F3}");
Console.WriteLine($"Γ(4.5) = {Gamma.Function(4.5):F3}");
```

### Log-Gamma Function

For large arguments, use log-gamma to avoid overflow:

```cs
// Regular gamma would overflow for large x
double x = 200.0;

// Use log-gamma
double logGamma = Gamma.LogGamma(x);
Console.WriteLine($"ln(Γ(200)) = {logGamma:F2}");

// Gamma itself would be huge: exp(logGamma)
// Don't compute directly - work in log space
```

### Digamma and Trigamma

Derivatives of the log-gamma function:

```cs
// Digamma: ψ(x) = d/dx[ln(Γ(x))] = Γ'(x)/Γ(x)
double digamma = Gamma.Digamma(2.0);

// Trigamma: ψ'(x) = d²/dx²[ln(Γ(x))]
double trigamma = Gamma.Trigamma(2.0);

Console.WriteLine($"ψ(2) = {digamma:F6}");
Console.WriteLine($"ψ'(2) = {trigamma:F6}");

// Applications: moment matching in Gamma distribution MLE
```

### Incomplete Gamma Functions

Used in chi-squared and gamma distributions:

```cs
// Lower incomplete gamma: γ(a,x) = ∫₀ˣ t^(a-1)e^(-t) dt
double lowerIncomplete = Gamma.LowerIncomplete(a: 2.0, x: 3.0);

// Upper incomplete gamma: Γ(a,x) = ∫ₓ^∞ t^(a-1)e^(-t) dt
double upperIncomplete = Gamma.UpperIncomplete(a: 2.0, x: 3.0);

// Verify: γ(a,x) + Γ(a,x) = Γ(a)
double sum = lowerIncomplete + upperIncomplete;
double gamma = Gamma.Function(2.0);

Console.WriteLine($"Lower: {lowerIncomplete:F6}");
Console.WriteLine($"Upper: {upperIncomplete:F6}");
Console.WriteLine($"Sum: {sum:F6}, Γ(2): {gamma:F6}");
```

### Regularized Incomplete Gamma

Normalized version: P(a,x) = γ(a,x)/Γ(a)

```cs
// Regularized incomplete gamma (CDF of Gamma distribution)
double P = Gamma.Incomplete(x: 3.0, alpha: 2.0);

Console.WriteLine($"P(2, 3) = {P:F6}");
Console.WriteLine("This equals the Gamma(2,1) CDF at x=3");

// Inverse: find x such that P(a,x) = p
double xInv = Gamma.InverseLowerIncomplete(a: 2.0, y: 0.9);
Console.WriteLine($"P(2, {xInv:F3}) = 0.9");
```

## Beta Function

The Beta function relates to the Gamma function: B(a,b) = Γ(a)Γ(b)/Γ(a+b)

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

### Incomplete Beta Function

Used in Beta distribution and Student's t-test:

```cs
// Incomplete beta: Bₓ(a,b) = ∫₀ˣ t^(a-1)(1-t)^(b-1) dt
double incompleteBeta = Beta.Incomplete(a: 2.0, b: 3.0, x: 0.4);

Console.WriteLine($"Bₓ(2,3,0.4) = {incompleteBeta:F6}");
Console.WriteLine("This equals integral from 0 to 0.4");

// Regularized incomplete beta (CDF of Beta distribution)
double I = incompleteBeta / Beta.Function(2.0, 3.0);
Console.WriteLine($"I(2,3,0.4) = {I:F6}");
```

### Inverse Incomplete Beta

Quantile function for Beta distribution:

```cs
// Find x such that I(a,b,x) = p
double x = Beta.IncompleteInverse(aa: 2.0, bb: 3.0, yy0: 0.5);

Console.WriteLine($"50th percentile of Beta(2,3): {x:F4}");

// Verify
double check = Beta.Incomplete(2.0, 3.0, x) / Beta.Function(2.0, 3.0);
Console.WriteLine($"Verification: I(2,3,{x:F4}) = {check:F6}");
```

## Error Function

The error function is the integral of the Gaussian distribution:

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
double cdf = Gamma.Incomplete(x: x / 2.0, alpha: k / 2.0);

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

// CDF = I(x;a,b) = Bₓ(a,b)/B(a,b)
double incompleteBeta = Beta.Incomplete(a, b, x);
double betaFunc = Beta.Function(a, b);
double cdf = incompleteBeta / betaFunc;

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
| **Beta** | Beta distribution | `Function()`, `Incomplete()` |
| **Error** | Normal distribution | `Function()`, `Erfc()`, `InverseErf()` |
| **Factorial** | Combinatorics | `Function()`, `BinomialCoefficient()` |

## Implementation Notes

- All functions use high-precision approximations
- Log-space variants prevent overflow for large arguments
- Inverse functions use Newton-Raphson iteration
- Special care for edge cases and numerical stability

---

[← Previous: Linear Algebra](linear-algebra.md) | [Back to Index](../index.md) | [Next: ODE Solvers →](ode-solvers.md)
