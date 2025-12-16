# Random Number Generation

[← Previous: Time Series](../data/time-series.md) | [Back to Index](../index.md)

The ***Numerics*** library provides multiple random number generation methods for different applications. These include high-quality pseudo-random generators, quasi-random sequences, and advanced sampling techniques.

## Pseudo-Random Number Generators

### Mersenne Twister

The Mersenne Twister is a high-quality pseudo-random number generator with excellent statistical properties [[1]](#1):

```cs
using Numerics.Sampling;

// Create with time-based seed
var rng = new MersenneTwister();

// Create with specific seed (for reproducibility)
var rng2 = new MersenneTwister(seed: 12345);

// Generate random integers
int randomInt = rng.Next();                    // [0, Int32.MaxValue)
int randomInRange = rng.Next(1, 100);          // [1, 100)

// Generate random doubles
double u1 = rng.NextDouble();                  // [0.0, 1.0)
double u2 = rng.GenRandReal1();                // [0, 1]
double u3 = rng.GenRandReal2();                // [0, 1)
double u4 = rng.GenRandReal3();                // (0, 1)
double u5 = rng.GenRandRes53();                // [0, 1) with 53-bit resolution

Console.WriteLine("Random numbers:");
Console.WriteLine($"Integer: {randomInt}");
Console.WriteLine($"Double [0,1): {u1:F6}");
Console.WriteLine($"Double (0,1): {u4:F6}");
```

**Properties:**
- Period: 2^19937 - 1 (extremely long)
- Excellent uniformity and independence
- Fast generation
- Reproducible with seeds

### Using with Distributions

```cs
using Numerics.Distributions;

var rng = new MersenneTwister(12345);

// Generate from distributions
var normal = new Normal(100, 15);
double[] samples = new double[1000];

for (int i = 0; i < samples.Length; i++)
{
    double u = rng.NextDouble();
    samples[i] = normal.InverseCDF(u);
}

Console.WriteLine($"Generated {samples.Length} Normal samples");
Console.WriteLine($"Mean: {samples.Average():F2}");
Console.WriteLine($"Std Dev: {Statistics.StandardDeviation(samples):F2}");
```

## Quasi-Random Sequences

### Sobol Sequence

Low-discrepancy quasi-random sequences for better coverage of parameter space [[2]](#2):

```cs
using Numerics.Sampling;

// Create 2D Sobol sequence
var sobol = new SobolSequence(dimension: 2);

Console.WriteLine("First 10 Sobol points:");
for (int i = 0; i < 10; i++)
{
    double[] point = sobol.NextDouble();
    Console.WriteLine($"Point {i}: ({point[0]:F4}, {point[1]:F4})");
}

// Skip to specific index
double[] pointAt100 = sobol.SkipTo(100);
Console.WriteLine($"\nPoint at index 100: ({pointAt100[0]:F4}, {pointAt100[1]:F4})");
```

**Properties:**
- Low discrepancy (better coverage than pseudo-random)
- Deterministic sequence
- Excellent for integration and optimization
- Converges faster than Monte Carlo

**When to use:**
- Numerical integration (better than Monte Carlo)
- Parameter space exploration
- Optimization initialization
- Sensitivity analysis

### Sobol vs. Pseudo-Random

```cs
int n = 100;
var random = new MersenneTwister(123);
var sobol = new SobolSequence(2);

// Generate pseudo-random points
var pseudoRandom = new List<(double, double)>();
for (int i = 0; i < n; i++)
{
    pseudoRandom.Add((random.NextDouble(), random.NextDouble()));
}

// Generate quasi-random points
var quasiRandom = new List<(double, double)>();
for (int i = 0; i < n; i++)
{
    double[] point = sobol.NextDouble();
    quasiRandom.Add((point[0], point[1]));
}

Console.WriteLine("Pseudo-random: Points may cluster");
Console.WriteLine("Quasi-random: Points evenly distributed");
Console.WriteLine("\nFor integration, quasi-random typically converges faster");
```

## Latin Hypercube Sampling

Stratified sampling for better parameter space coverage [[3]](#3):

```cs
using Numerics.Sampling;

// Generate Latin Hypercube sample
int sampleSize = 50;
int dimensions = 3;
int seed = 12345;

// Random LHS
double[,] lhsRandom = LatinHypercube.Random(sampleSize, dimensions, seed);

// Median LHS (centered in strata)
double[,] lhsMedian = LatinHypercube.Median(sampleSize, dimensions, seed);

Console.WriteLine("Latin Hypercube Sample (Random):");
Console.WriteLine("Sample | Dim 1  | Dim 2  | Dim 3");
Console.WriteLine("-------|--------|--------|--------");

for (int i = 0; i < Math.Min(10, sampleSize); i++)
{
    Console.WriteLine($"{i,6} | {lhsRandom[i, 0],6:F4} | {lhsRandom[i, 1],6:F4} | {lhsRandom[i, 2],6:F4}");
}

// Transform to actual distributions
var normal = new Normal(100, 15);
var lognormal = new LogNormal(4, 0.5);
var uniform = new Uniform(0, 10);

double[,] transformed = new double[sampleSize, 3];
for (int i = 0; i < sampleSize; i++)
{
    transformed[i, 0] = normal.InverseCDF(lhsRandom[i, 0]);
    transformed[i, 1] = lognormal.InverseCDF(lhsRandom[i, 1]);
    transformed[i, 2] = uniform.InverseCDF(lhsRandom[i, 2]);
}

Console.WriteLine("\nTransformed to distributions:");
Console.WriteLine("Sample | Normal | LogNormal | Uniform");
for (int i = 0; i < Math.Min(5, sampleSize); i++)
{
    Console.WriteLine($"{i,6} | {transformed[i, 0],6:F1} | {transformed[i, 1],9:F2} | {transformed[i, 2],7:F2}");
}
```

**Properties:**
- Stratified sampling (one sample per stratum)
- Better coverage than simple random sampling
- Reduced variance in estimates
- Efficient for small sample sizes

**When to use:**
- Monte Carlo simulation with limited budget
- Sensitivity analysis
- Calibration with expensive models
- Risk assessment studies

## Practical Examples

### Example 1: Monte Carlo Integration

```cs
// Integrate f(x) = x² from 0 to 1 using different methods

Func<double, double> f = x => x * x;
int n = 1000;

// Simple Monte Carlo (pseudo-random)
var rng = new MersenneTwister(123);
double mcSum = 0;
for (int i = 0; i < n; i++)
{
    mcSum += f(rng.NextDouble());
}
double mcEstimate = mcSum / n;  // Approximates ∫₀¹ x² dx = 1/3

// Quasi-Monte Carlo (Sobol)
var sobol = new SobolSequence(1);
double qmcSum = 0;
for (int i = 0; i < n; i++)
{
    qmcSum += f(sobol.NextDouble()[0]);
}
double qmcEstimate = qmcSum / n;

double exact = 1.0 / 3.0;

Console.WriteLine("Monte Carlo Integration of x² from 0 to 1:");
Console.WriteLine($"Exact value: {exact:F6}");
Console.WriteLine($"MC estimate: {mcEstimate:F6} (error: {Math.Abs(mcEstimate - exact):E4})");
Console.WriteLine($"QMC estimate: {qmcEstimate:F6} (error: {Math.Abs(qmcEstimate - exact):E4})");
Console.WriteLine("\nQMC typically has smaller error for same sample size");
```

### Example 2: Uncertainty Propagation

```cs
// Model: y = a*x + b*x² where a, b are uncertain

var normal_a = new Normal(2.0, 0.3);
var normal_b = new Normal(1.0, 0.1);
double x = 5.0;

// Monte Carlo with pseudo-random
var rng = new MersenneTwister(123);
int nSamples = 10000;

double[] y_mc = new double[nSamples];
for (int i = 0; i < nSamples; i++)
{
    double a = normal_a.InverseCDF(rng.NextDouble());
    double b = normal_b.InverseCDF(rng.NextDouble());
    y_mc[i] = a * x + b * x * x;
}

// Latin Hypercube Sampling
var lhs = LatinHypercube.Random(nSamples, 2, seed: 123);

double[] y_lhs = new double[nSamples];
for (int i = 0; i < nSamples; i++)
{
    double a = normal_a.InverseCDF(lhs[i, 0]);
    double b = normal_b.InverseCDF(lhs[i, 1]);
    y_lhs[i] = a * x + b * x * x;
}

Console.WriteLine("Uncertainty Propagation:");
Console.WriteLine($"MC  - Mean: {y_mc.Average():F2}, Std: {Statistics.StandardDeviation(y_mc):F2}");
Console.WriteLine($"LHS - Mean: {y_lhs.Average():F2}, Std: {Statistics.StandardDeviation(y_lhs):F2}");
Console.WriteLine("\nLHS typically more stable with fewer samples");
```

### Example 3: Global Optimization Initialization

```cs
// Initialize population for global optimization

int popSize = 20;
int dimensions = 3;

// Parameter bounds
double[] lowerBounds = { -10, -5, 0 };
double[] upperBounds = { 10, 5, 100 };

// Generate initial population with LHS
var lhs = LatinHypercube.Random(popSize, dimensions, seed: 123);

// Scale to actual bounds
double[,] population = new double[popSize, dimensions];
for (int i = 0; i < popSize; i++)
{
    for (int d = 0; d < dimensions; d++)
    {
        population[i, d] = lowerBounds[d] + lhs[i, d] * (upperBounds[d] - lowerBounds[d]);
    }
}

Console.WriteLine("Initial Population for Optimization:");
Console.WriteLine("Individual | Param 1 | Param 2 | Param 3");
Console.WriteLine("-----------|---------|---------|----------");

for (int i = 0; i < Math.Min(10, popSize); i++)
{
    Console.WriteLine($"{i,10} | {population[i, 0],7:F2} | {population[i, 1],7:F2} | {population[i, 2],8:F2}");
}

Console.WriteLine("\nLHS ensures good coverage of parameter space");
```

### Example 4: Sensitivity Analysis

```cs
// Compute Sobol sensitivity indices using quasi-random sampling

Func<double[], double> model = x => 
    x[0] + 2 * x[1] + 3 * x[2] + 4 * x[1] * x[2];

int n = 1000;
int dim = 3;

// Generate two independent LHS samples
var A = LatinHypercube.Random(n, dim, seed: 123);
var B = LatinHypercube.Random(n, dim, seed: 456);

// Evaluate model
double[] yA = new double[n];
double[] yB = new double[n];

for (int i = 0; i < n; i++)
{
    yA[i] = model(new[] { A[i, 0], A[i, 1], A[i, 2] });
    yB[i] = model(new[] { B[i, 0], B[i, 1], B[i, 2] });
}

// First-order sensitivity indices
double varY = Statistics.Variance(yA);

Console.WriteLine("Sensitivity Analysis:");
Console.WriteLine("Parameter | First-Order Index");
Console.WriteLine("----------|------------------");

for (int j = 0; j < dim; j++)
{
    double[] yABj = new double[n];
    
    for (int i = 0; i < n; i++)
    {
        double[] x = new double[dim];
        for (int k = 0; k < dim; k++)
        {
            x[k] = (k == j) ? B[i, k] : A[i, k];
        }
        yABj[i] = model(x);
    }
    
    double S_j = (yA.Zip(yABj, (ya, yabj) => ya * yabj).Average() - 
                  yA.Average() * yA.Average()) / varY;
    
    Console.WriteLine($"Param {j + 1}  | {S_j,17:F4}");
}
```

## Choosing a Random Number Generator

| Method | Use Case | Pros | Cons |
|--------|----------|------|------|
| **Mersenne Twister** | General purpose | Fast, long period | Clusters in high dimensions |
| **Sobol** | Integration, optimization | Low discrepancy, deterministic | Not random, dimension limit |
| **Latin Hypercube** | Small sample studies | Stratified, efficient | Requires planning |

### Decision Guide

**Use Mersenne Twister when:**
- Need standard random numbers
- Simulating stochastic processes
- Large sample sizes available
- Distribution sampling

**Use Sobol when:**
- Numerical integration
- Parameter space exploration
- Deterministic sequence needed
- Integration convergence critical

**Use Latin Hypercube when:**
- Limited computational budget
- Sensitivity analysis
- Need efficient stratification
- Small to medium samples (10-1000)

## Reproducibility

### Setting Seeds

```cs
// Pseudo-random - use same seed for reproducibility
var rng1 = new MersenneTwister(12345);
var rng2 = new MersenneTwister(12345);

// Generate same sequence
for (int i = 0; i < 5; i++)
{
    double r1 = rng1.NextDouble();
    double r2 = rng2.NextDouble();
    Console.WriteLine($"RNG1: {r1:F6}, RNG2: {r2:F6}, Same: {r1 == r2}");
}

// Quasi-random - deterministic by design
var sobol1 = new SobolSequence(2);
var sobol2 = new SobolSequence(2);

for (int i = 0; i < 3; i++)
{
    var point1 = sobol1.NextDouble();
    var point2 = sobol2.NextDouble();
    Console.WriteLine($"Sobol sequences identical: {point1[0] == point2[0] && point1[1] == point2[1]}");
}
```

## Best Practices

1. **Always set seeds** for reproducible research
2. **Document RNG choices** in methods section
3. **Use appropriate method** for application
4. **Check sample size** requirements
5. **Validate distributions** with statistical tests
6. **Consider quasi-random** for integration
7. **Use LHS** for expensive models

## Performance Considerations

- **Mersenne Twister**: Very fast, suitable for large samples
- **Sobol**: Slightly slower, excellent for moderate samples
- **Latin Hypercube**: Overhead for stratification, best for small samples

---

## References

<a id="1">[1]</a> Matsumoto, M., & Nishimura, T. (1998). Mersenne twister: a 623-dimensionally equidistributed uniform pseudo-random number generator. *ACM Transactions on Modeling and Computer Simulation*, 8(1), 3-30.

<a id="2">[2]</a> Sobol, I. M. (1967). On the distribution of points in a cube and the approximate evaluation of integrals. *USSR Computational Mathematics and Mathematical Physics*, 7(4), 86-112.

<a id="3">[3]</a> McKay, M. D., Beckman, R. J., & Conover, W. J. (1979). A comparison of three methods for selecting values of input variables in the analysis of output from a computer code. *Technometrics*, 21(2), 239-245.

---

[← Previous: Time Series](../data/time-series.md) | [Back to Index](../index.md)
