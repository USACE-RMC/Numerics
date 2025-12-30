# Numerical Integration

[← Back to Index](../index.md) | [Next: Numerical Differentiation →](differentiation.md)

Numerical integration, also known as numerical quadrature, is a fundamental technique for approximating definite integrals. It has wide-ranging applications in various scientific and engineering fields. For example, in statistics, the expected value of a random variable is calculated using an integral, and numerical integration can be employed to approximate this expected value. Many problems in engineering and physics cannot be solved analytically and must rely on numerical methods to approximate solutions.

## Single Dimension Integration

The ***Numerics*** library provides several methods for performing numerical integration on single dimensional integrands. Each algorithm computes an approximation to a definite integral of the form:

```math
I = \int\limits_{a}^{b}f(x) \cdot dx
```

For the first example, let's consider a simple function with a single variable:

```math
f(x)=x^3
```

Integrating from $a=0$ to $b=1$ yields the exact solution:

```math
\int\limits_{0}^{1}f(x) \cdot dx = \frac{1}{4}x^4 \biggr|_0^1 = \frac{1}{4} \cdot 1^4 - 0 = 0.25
```

Definite integrals can be numerically solved using Riemann sums, such as the trapezoidal rule. This method works by approximating the region under the function $f(x)$ as a trapezoid and calculating its area:

```math
I =\int\limits_{0}^{1}f(x) \cdot  dx \approx \left(\frac{f(a) + f(b)}{2} \right)\cdot(b-a)
```

This approximation can be improved by partitioning (or binning) the integration interval $[a,b]$ and then applying the trapezoidal rule to each subinterval and summing the results:

```math
I =\int\limits_{0}^{1}f(x) \cdot dx \approx \sum_{i=1}^{N} \left(\frac{f(x_{i-1}) + f(x_i)}{2} \right)\cdot(x_i-x_{i-1})
```

Now, let's implement this in ***Numerics***. First, we need to reference the _Integration_ namespace:

```cs
using Numerics.Mathematics.Integration;
```

Next, create the test function:

```cs
/// <summary>
/// Test function: f(x) = x^3
/// </summary>
public double FX(double x)
{
    return Math.Pow(x, 3);
}
```

### Trapezoidal Rule

The _Integration_ class is a static class that contains the Midpoint Rule, Trapezoidal Rule, Simpson's Rule, and the 10-point Gauss-Legendre integration methods. Let's first compute the integral using the Trapezoidal Rule with 10 bins (or steps):

```cs
double result = Integration.TrapezoidalRule(FX, 0, 1, 10); // 0.25249999999999995
```

Increasing the number of steps will increase the accuracy. Let's compute it again using 1,000 steps:

```cs
double result = Integration.TrapezoidalRule(FX, 0, 1, 1000); // 0.25000025000000053
```

We can see that this is much more precise.

Alternatively, the ***Numerics*** library provides a _TrapezoidalRule_ class that extends the basic static method and provides additional functionality for computing integration error estimates:

```cs
var trap = new TrapezoidalRule(FX, 0, 1, intervals: 100);
trap.Integrate();
double result = trap.Result; // 0.25024999999999995
Console.WriteLine($"Function evaluations: {trap.FunctionEvaluations}");
```

### Simpson's Rule

Simpson's Rule provides more accurate approximations than the Trapezoidal Rule by using quadratic interpolation. The formula approximates the region under the function $f(x)$ as a weighted average of the trapezoidal and midpoint methods:

```math
I =\int\limits_{a}^{b}f(x) \cdot dx \approx \left[f(a) + 4 \cdot f \left(\frac{a+b}{2} \right) + f(b) \right]\cdot \left(\frac{b-a}{6} \right)
```

Similar to the Trapezoidal Rule, the accuracy is improved by partitioning the integration interval. Using the static method:

```cs
double result = Integration.SimpsonsRule(FX, 0, 1, 10); // 0.25
```

Or using the class-based approach:

```cs
var simpson = new SimpsonsRule(FX, 0, 1, intervals: 100);
simpson.Integrate();
double result = simpson.Result; // 0.25
```

**Error**: Simpson's Rule is fourth-order accurate, $O(h^4)$, making it significantly more accurate than the Trapezoidal Rule for smooth functions.

### Gauss-Legendre Quadrature

The Gauss-Legendre method uses optimal polynomial quadrature for smooth functions [[1]](#1). It evaluates the function at specific points (roots of Legendre polynomials) with corresponding weights to achieve high accuracy:

```math
\int_{-1}^{1} f(x)\,dx \approx \sum_{i=1}^{n} w_i f(x_i)
```

where $x_i$ are roots of Legendre polynomials and $w_i$ are corresponding weights. The ***Numerics*** library provides a 10-point Gauss-Legendre method:

```cs
double result = Integration.GaussLegendre(FX, 0, 1); // 0.25
```

**Error**: The 10-point Gauss-Legendre method is exact for polynomials of degree 19 or less.

### Midpoint Rule

The Midpoint Rule evaluates the function at the midpoint of each interval:

```cs
double result = Integration.Midpoint(FX, 0, 1, 100); // 0.24997500000000012
```

## Adaptive Integration

The challenge with static numerical integration methods, such as the trapezoidal rule mentioned above, is that the user must specify both the limits of integration and the number of integration bins. If the integrand function has subregions with high variance, this approach can lead to large approximation errors. Many real-world integrand functions have substantial weight concentrated in narrow subregions, resulting in wasted integration bins in areas that contribute little to the total weight.

Adaptive integration, a more refined numerical integration method, adjusts subintervals within the integration bounds based on the behavior of the function. These methods concentrate subintervals in regions that contribute the most to the integral, overcoming the limitations of static approaches.

The ***Numerics*** library provides three adaptive integration routines: the Adaptive Simpson's Rule, the Adaptive Gauss-Lobatto method, and the Adaptive Gauss-Kronrod method.

### Adaptive Simpson's Rule

The Adaptive Simpson's Rule (ASR) algorithm subdivides the integration interval recursively until a user-defined tolerance is achieved. In each subinterval, Simpson's Rule is used to approximate the region under the function. The criterion for determining when to stop subdividing an interval is:

```math
\frac{1}{15} \cdot \left| S(a, m) + S(m, b) - S(a,b) \right| \leq \epsilon + \epsilon \cdot \left| S(a, b) \right|
```

where $[a,b]$ is the integration interval, $m = \frac{a+b}{2}$, $S(\cdot)$ represents Simpson's Rule evaluated at those intervals, and $\epsilon$ is the absolute and relative error tolerance for the interval. Each subinterval is recursively subdivided and evaluated until the specified tolerance is met.

More details on the ASR method can be found in [[1]](#1).

To use the ASR method, follow these steps:

```cs
var asr = new AdaptiveSimpsonsRule(FX, 0, 1);
asr.Integrate();
double result = asr.Result; // 0.25
Console.WriteLine($"Function evaluations: {asr.FunctionEvaluations}");
```

For this simple test function, the ASR method requires only 5 function evaluations to converge with an absolute and relative tolerance of $1 \times 10^{-8}$. It should be noted that the ASR method gives exact results for 3rd degree (or less) polynomials.

You can customize the tolerance and iteration limits:

```cs
var asr = new AdaptiveSimpsonsRule(FX, 0, 1);
asr.RelativeTolerance = 1e-10;
asr.AbsoluteTolerance = 1e-10;
asr.MaxIterations = 1000;
asr.Integrate();
double result = asr.Result;
```

### Adaptive Gauss-Lobatto

The Adaptive Gauss-Lobatto (AGL) method includes endpoints in quadrature nodes, making it useful for integrands with endpoint singularities [[2]](#2). Alternatively, we can use the AGL method as follows:

```cs
var agl = new AdaptiveGaussLobatto(FX, 0, 1);
agl.Integrate();
double result = agl.Result; // 0.24999999999999997
Console.WriteLine($"Function evaluations: {agl.FunctionEvaluations}");
```

The AGL method requires 18 function evaluations to converge given an absolute and relative tolerance of $1 \times 10^{-8}$.

### Adaptive Gauss-Kronrod

The Adaptive Gauss-Kronrod method pairs a 10-point Gauss rule with a 21-point Kronrod extension for error estimation [[2]](#2). This method efficiently reuses function evaluations to provide accurate error estimates:

```cs
var gk = new AdaptiveGaussKronrod(FX, 0, 1);
gk.RelativeTolerance = 1e-12;
gk.Integrate();
double result = gk.Result;
Console.WriteLine($"Estimated error: {gk.Status}");
```

**Advantages**: Efficient error estimation, reuses function evaluations from the Gauss rule in the Kronrod extension.

### Example: Computing the Mean of a Gamma Distribution

For a more challenging test problem, let's compute the mean of a Gamma distribution with a scale of $\theta = 10$ and shape $\kappa = 5$. The true mean of the distribution is given by:

```math
\mu = \theta \cdot \kappa = 50
```

The probability density function (PDF) of the Gamma distribution is:

```math
f(x) = \frac{1}{\Gamma(\kappa)\theta^{\kappa}}x^{\kappa-1}e^{-\frac{x}{\theta}}
```

The mean of a continuous probability distribution is computed as:

```math
\mu = \mathbb{E} [X] = \int\limits_{-\infty}^{\infty} x \cdot f(x) \cdot  dx
```

Now, let's implement this in ***Numerics***. First, we need to reference the _Integration_ and _Distributions_ namespaces:

```cs
using Numerics.Mathematics.Integration;
using Numerics.Distributions;
```

Then, using the ASR method, follow these steps:

```cs
// Create the Gamma distribution and set the integration limits
var gamma = new GammaDistribution(10, 5);
double a = gamma.InverseCDF(1E-16); // Lower limit based on a very small cumulative probability
double b = gamma.InverseCDF(1 - 1E-16); // Upper limit based on a near-1 cumulative probability

// Define the integrand function
double I(double x)
{
    return x * gamma.PDF(x);
}

// Perform the integration
var asr = new AdaptiveSimpsonsRule(I, a, b);
asr.Integrate();
double result = asr.Result; // 50.000000004866415
Console.WriteLine($"Function evaluations: {asr.FunctionEvaluations}");
```

The ASR method requires 365 function evaluations to reach convergence.

## Two-Dimensional Integration

### Adaptive Simpson's 2D

The ***Numerics*** library extends adaptive Simpson's rule to rectangular domains in two dimensions:

```cs
using Numerics.Mathematics.Integration;

// Integrate f(x,y) = exp(-(x² + y²)) over [-3,3] × [-3,3]
Func<double, double, double> f2d = (x, y) => Math.Exp(-(x * x + y * y));

var simpson2d = new AdaptiveSimpsonsRule2D(f2d, -3, 3, -3, 3);
simpson2d.RelativeTolerance = 1e-8;
simpson2d.Integrate();

Console.WriteLine($"Result: {simpson2d.Result:F10}");  // Should be approximately π
Console.WriteLine($"Function evaluations: {simpson2d.FunctionEvaluations}");
```

## Multidimensional Integration

Multidimensional integration, also known as multiple or multivariate integration, involves evaluating integrals over functions of more than one variable. Instead of integrating over a single interval, as in one-dimensional integration, you integrate over a region in a multidimensional space. This is commonly used in fields like physics, engineering, and statistics where systems often depend on multiple variables.

Solving multidimensional integrals is computationally demanding. If traditional, nonadaptive numerical integration techniques were used, the solution would require $K^D$ iterations, where $K$ is the number of integration steps (or bins) and $D$ is the number of dimensions. If there were 100 integration steps and 5 dimensions, the solution would need 10 billion iterations.

To avoid these computation limitations, the ***Numerics*** library provides three multidimensional integration routines: Monte Carlo, Miser, and VEGAS. Each algorithm computes an approximation to a definite integral of the form:

```math
I = \int_{\Omega} f(\mathbf{x}) \, d\mathbf{x}
```

where $\Omega$ is the integration domain in $D$-dimensional space.

### Monte Carlo Integration

Monte Carlo integration uses random sampling to approximate integrals. For high-dimensional integrals, Monte Carlo methods become essential [[3]](#3):

```math
\int_\Omega f(\mathbf{x})\,d\mathbf{x} \approx V \cdot \frac{1}{N}\sum_{i=1}^{N}f(\mathbf{x}_i)
```

where $V$ is the volume of the integration domain, $N$ is the number of samples, and $\mathbf{x}_i$ are random points uniformly distributed in $\Omega$.

Let's use a simple 2D test problem to compute $\pi$. Consider the function:

```math
f(x, y) = 
\begin{cases}
1 & \text{if } x^2 + y^2 \leq 1 \\
0 & \text{otherwise}
\end{cases}
```

Integrating this function over the domain $[-1, 1] \times [-1, 1]$ gives the area of a unit circle, which equals $\pi$:

```math
\int\limits_{-1}^{1}\int\limits_{-1}^{1}f(x,y) \, dy \, dx = \pi
```

Now, let's implement this in ***Numerics***:

```cs
using Numerics.Mathematics.Integration;
using Numerics.Sampling;

// Define the integrand function
double PI(double[] x)
{
    if (x[0] * x[0] + x[1] * x[1] <= 1)
        return 1.0;
    return 0.0;
}

// Set integration bounds
var a = new double[] { -1, -1 };
var b = new double[] { 1, 1 };

// Create and configure the Monte Carlo integrator
var mc = new MonteCarloIntegration(PI, 2, a, b);
mc.Random = new MersenneTwister(12345); // Set the random number generator for repeatability
mc.MaxIterations = 100000;
mc.Integrate();

double result = mc.Result; // 3.13824
Console.WriteLine($"Result: {result:F6}");
Console.WriteLine($"Function evaluations: {mc.FunctionEvaluations}");
```

With 100,000 samples, we see that the result is close but still has a noticeable error. Now, let's run it again with the default setting, where the maximum iterations are $N=100,000,000$:

```cs
var mc = new MonteCarloIntegration(PI, 2, a, b);
mc.Random = new MersenneTwister(12345); // Set the random number generator for repeatability
mc.Integrate();
double result = mc.Result; // 3.1412028
Console.WriteLine($"Result: {result:F8}");
```

This result is much closer to the true value of $\pi$.

Unlike traditional methods, the complexity of Monte Carlo integration grows slowly with the number of dimensions, making it particularly useful for high-dimensional problems. The Monte Carlo approach is simple to implement in higher dimensions and can handle irregular domains and complex integrands. However, it converges slowly; the error decreases as $O \left( \frac{1}{\sqrt{N}} \right)$, meaning to halve the error, you need to quadruple the number of samples.

**Error**: $O(1/\sqrt{N})$ - independent of dimension.

### MISER (Recursive Stratified Sampling)

The Miser integration algorithm is a type of adaptive Monte Carlo method designed for efficient evaluation of multidimensional integrals. It is particularly well-suited for integrands that exhibit regions of high variance, as it allocates more samples to areas where the integrand contributes more to the total integral. The algorithm combines the flexibility of Monte Carlo integration with adaptive subdivision techniques to enhance accuracy and efficiency in complex, high-dimensional problems.

Key Concepts of the Miser Algorithm:

1. **Adaptive Subdivision:** Miser improves upon basic Monte Carlo integration by recursively subdividing the integration domain into smaller regions. The algorithm then allocates more samples to the subregions where the integrand has higher variance, focusing computational resources where they are most needed.

2. **Variance-Based Sampling:** The Miser algorithm estimates the variance of the integrand in different subregions. Subregions with higher variance are given a greater proportion of the total samples. This reduces the error by refining the integral in the parts of the domain that contribute the most to the integral's value.

For more details on the stratified sampling and the Miser algorithm, see [[2]](#2) and [[3]](#3).

Now, let's solve the $\pi$ test function using Miser with $N=100,000$ iterations:

```cs
var miser = new Miser(PI, 2, a, b);
miser.Random = new MersenneTwister(12345); // Set the random number generator for repeatability
miser.MaxIterations = 100000;
miser.Integrate();
double result = miser.Result; // 3.1420673978501474
Console.WriteLine($"Result: {result:F10}");
Console.WriteLine($"Function evaluations: {miser.FunctionEvaluations}");
```

With the same number of samples, Miser produces a more accurate result with smaller variance than basic Monte Carlo integration.

### VEGAS (Adaptive Importance Sampling)

The VEGAS integration method is a Monte Carlo-based numerical integration technique designed for efficiently evaluating high-dimensional integrals, particularly when dealing with functions that have significant variability in certain regions of the integration space [[4]](#4) [[5]](#5). It is widely used in computational physics and other fields requiring the evaluation of complex integrals.

Key Features of the VEGAS Algorithm:

1. **Importance Sampling:** VEGAS employs importance sampling to focus the integration effort on regions where the integrand contributes most significantly to the integral. This helps to improve the accuracy of the integral estimate while reducing variance.

2. **Adaptive Grid:** The algorithm adapts the sampling grid based on the characteristics of the integrand. It divides the integration domain into smaller subregions, and the sampling density is adjusted according to the estimated contribution of each region to the overall integral.

3. **Iterative Approach:** VEGAS works in iterations, refining the sampling strategy with each pass. In the first iteration, a uniform grid is typically used. After evaluating the integrand, the method estimates the probability distribution of the function values, allowing the grid to be adjusted in subsequent iterations to better capture areas with higher contributions.

For more details on the importance sampling and the VEGAS algorithm, see [[2]](#2) and [[3]](#3).

Now, let's solve the $\pi$ test function using VEGAS. The VEGAS method requires the integrand function to take a point $\mathbf{x}$ and an importance sampling weight $w$ as inputs. For this example, we will reuse the previous test function without utilizing the weight value:

```cs
var vegas = new Vegas((x, w) => { return PI(x); }, 2, a, b);
vegas.Random = new MersenneTwister(12345); // Set the random number generator for repeatability
vegas.Integrate();
double result = vegas.Result; // 3.1418009008273735
Console.WriteLine($"Result: {result:F10}");
Console.WriteLine($"Function evaluations: {vegas.FunctionEvaluations}");
```

In ***Numerics***, the VEGAS method iteratively adapts and refines the grid until convergence to a relative tolerance of $1 \times 10^{-3}$ is achieved. For this test problem, only 19,600 function evaluations are required to reach convergence.

**Advantages**: Excellent for peaked integrands, provides convergence diagnostics through chi-squared statistics.

### Example: Mean of Sum of Independent Normal Distributions

For a more challenging test, let's compute the mean of the sum of independent Normal distributions. The multidimensional integrand function can be written as:

```math
f(x_1, \cdots ,x_D) = \sum_{k=1}^{D} x_k \cdot \prod_{k=1}^{D} \phi (x_k | \mu_k, \sigma_k)
```

where $\phi(\cdot)$ is the PDF of the $k$-th Normal distribution with a mean $\mu_k$ and standard deviation $\sigma_k$. The exact solution for the mean of the sum of these random variables is:

```math
E[X] = \sum_{k=1}^{D} \mu_k
```

For this test, we use five Normal distributions with means $\mu = [10, 30, 17, 99, 68]$ and standard deviations $\sigma = [2, 15, 5, 14, 7]$. Therefore, the exact solution is:

```math
E[X] = 10+30+17+99+68=224
```

Here's how we can implement this in the ***Numerics*** library:

```cs
using Numerics.Mathematics.Integration;
using Numerics.Distributions;
using Numerics.Sampling;

// Create the Normal distributions and set the integration limits
var mu = new double[] { 10, 30, 17, 99, 68 };
var sigma = new double[] { 2, 15, 5, 14, 7 };
var dists = new Normal[5];
var min = new double[5];
var max = new double[5];
for (int i = 0; i < 5; i++)
{
    dists[i] = new Normal(mu[i], sigma[i]);
    min[i] = dists[i].InverseCDF(1E-16); // Lower limit based on a very small cumulative probability
    max[i] = dists[i].InverseCDF(1 - 1E-16); // Upper limit based on a near-1 cumulative probability
}

// Define the integrand function
double SumOfNormals(double[] x, double w)
{
    double sum = 0;
    double prod = 1;
    for (int i = 0; i < mu.Length; i++)
    {
        sum += x[i];
        prod *= dists[i].PDF(x[i]);
    }
    return sum * prod;
}

// Perform the integration
var vegas = new Vegas(SumOfNormals, 5, min, max);
vegas.Integrate();
double result = vegas.Result; // 224.07455771892427
Console.WriteLine($"Result: {result:F8}");
Console.WriteLine($"Function evaluations: {vegas.FunctionEvaluations}");
```

For this more complex test problem, 468,750 function evaluations are required to achieve convergence.

## Probability Space Integration

For risk analysis, it's often useful to integrate over probability distributions. The ***Numerics*** library allows you to compute expected values and failure probabilities by transforming the integration domain to probability space:

```cs
using Numerics.Distributions;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;

// Expected value E[g(X)] where X ~ Normal(0,1)
var normal = new Normal(0, 1);

double Expectation(double[] p)
{
    // Transform from [0,1] to real line via inverse CDF
    double x = normal.InverseCDF(p[0]);
    
    // Function of interest: g(x) = x²
    double gx = x * x;  // E[X²] = 1 for standard normal
    
    return gx;  // Weight is implicitly 1 in probability space
}

var mc = new MonteCarloIntegration(Expectation, 1, new[] { 0.0 }, new[] { 1.0 });
mc.MaxIterations = 100000;
mc.Integrate();

Console.WriteLine($"E[X²] = {mc.Result:F6}");  // Should be ≈ 1
```

### Failure Probability Integration

```cs
// P(g(X) < 0) where g is a limit state function
double FailureIndicator(double[] p)
{
    double x1 = new Normal(0, 1).InverseCDF(p[0]);
    double x2 = new Normal(0, 1).InverseCDF(p[1]);
    
    // Limit state function
    double g = 5 - x1 - x2;
    
    // Indicator function
    return g < 0 ? 1.0 : 0.0;
}

var vegas = new Vegas((x, w) => FailureIndicator(x), 2, new[] { 0.0, 0.0 }, new[] { 1.0, 1.0 });
vegas.MaxIterations = 50000;
vegas.Integrate();

Console.WriteLine($"Failure probability: {vegas.Result:E4}");
```

## Choosing an Integration Method

| Scenario | Recommended Method | Notes |
|----------|-------------------|-------|
| 1D, smooth function | Adaptive Gauss-Kronrod | Best accuracy with error estimates |
| 1D, unknown smoothness | Adaptive Simpson's | Reliable general-purpose method |
| 1D, endpoint singularity | Adaptive Gauss-Lobatto | Includes endpoints in quadrature |
| 1D, high accuracy needed | Adaptive Gauss-Kronrod | Tight tolerance, efficient reuse of evaluations |
| 1D, quick estimate | Static Simpson's or Gauss-Legendre | Fast but requires choosing number of points |
| 2D, rectangular domain | Adaptive Simpson's 2D | Direct extension of 1D method |
| 3D-6D, smooth | VEGAS | Adaptive importance sampling excels here |
| High-D (>6), any | VEGAS or Monte Carlo | Only practical methods for very high dimensions |
| Peaked integrand | VEGAS | Adaptive grid focuses on important regions |
| High variance regions | Miser | Stratified sampling concentrates effort |
| Probability integrals | Monte Carlo or VEGAS in [0,1] | Transform via inverse CDF |

## Common Properties for Integrators

All integration classes (except the static `Integration` methods) inherit from the `Integrator` base class and share these properties:

### Input Properties
- `MinIterations`: Minimum number of iterations allowed (default = 1)
- `MaxIterations`: Maximum number of iterations allowed (default = 10,000,000)
- `MinFunctionEvaluations`: Minimum function evaluations allowed (default = 1)
- `MaxFunctionEvaluations`: Maximum function evaluations allowed (default = 10,000,000)
- `AbsoluteTolerance`: Desired absolute tolerance (default = 1E-8)
- `RelativeTolerance`: Desired relative tolerance (default = 1E-8)
- `ReportFailure`: Whether to throw exception on convergence failure (default = true)

### Output Properties
- `Result`: The computed integral value
- `Iterations`: Number of iterations performed
- `FunctionEvaluations`: Number of function evaluations performed
- `Status`: Integration status (Success, Failure, MaxIterationsReached, etc.)

Example of using these properties:

```cs
var integrator = new AdaptiveSimpsonsRule(myFunction, 0, 1);
integrator.RelativeTolerance = 1e-12;
integrator.MaxIterations = 5000;
integrator.ReportFailure = false; // Don't throw exception on failure
integrator.Integrate();

if (integrator.Status == IntegrationStatus.Success)
{
    Console.WriteLine($"Result: {integrator.Result}");
    Console.WriteLine($"Function evaluations: {integrator.FunctionEvaluations}");
}
else
{
    Console.WriteLine($"Integration failed: {integrator.Status}");
}
```

## Additional VEGAS Features

The VEGAS integrator includes a special method for rare event analysis:

```cs
var vegas = new Vegas(myFunction, dimensions, min, max);
vegas.ConfigureForRareEvents(); // Optimizes settings for rare event detection
vegas.Integrate();
```

This method adjusts the internal parameters to better handle integrands with very small contributions over most of the domain, which is common in reliability analysis and rare event simulation.

---

## References

<a id="1">[1]</a> P. J. Davis and P. Rabinowitz, *Methods of Numerical Integration*, 2nd ed., Mineola, New York: Dover Publications, Inc., 2007.

<a id="2">[2]</a> W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery, *Numerical Recipes: The Art of Scientific Computing*, 3rd ed., Cambridge, UK: Cambridge University Press, 2007.

<a id="3">[3]</a> A. Ciric, *A Guide to Monte Carlo & Quantum Monte Carlo Methods*, Createspace Independent Publishing Platform, 2016.

<a id="4">[4]</a> G. Lepage, "A New Algorithm for Adaptive Multidimensional Integration," *Journal of Computational Physics*, vol. 27, no. 1, pp. 192-203, 1978.

<a id="5">[5]</a> G. Lepage, "VEGAS: An Adaptive Multidimensional Integration Program," Cornell University, 1980.

---

[← Back to Index](../index.md) | [Next: Numerical Differentiation →](differentiation.md)
