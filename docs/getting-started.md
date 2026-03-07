# Getting Started

This guide will help you get up and running with the ***Numerics*** library quickly.

## Installation

### NuGet Package

The easiest way to install ***Numerics*** is via NuGet:

```bash
dotnet add package RMC.Numerics
```

Or using the Package Manager Console:

```
Install-Package RMC.Numerics
```

### Manual Installation

Download the compiled DLL from the releases page and add a reference to your project:

```xml
<Reference Include="Numerics">
  <HintPath>Numerics.dll</HintPath>
</Reference>
```

## Required Namespaces

Import the namespaces you need at the top of your C# files:

```cs
// Core distributions
using Numerics.Distributions;

// Statistical functions
using Numerics.Data.Statistics;

// Numerical methods
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.RootFinding;

// Sampling and MCMC
using Numerics.Sampling;
using Numerics.Sampling.MCMC;

// Interpolation
using Numerics.Data;
```

## Working with Distributions

### Creating a Distribution

All univariate distributions implement a consistent interface. You can create distributions by specifying parameters directly:

```cs
using Numerics.Distributions;

// Normal distribution with mean=100, std=15
var normal = new Normal(100, 15);

// Generalized Extreme Value with location=1000, scale=200, shape=-0.1
var gev = new GeneralizedExtremeValue(1000, 200, -0.1);

// Log-Pearson Type III (commonly used in hydrology)
var lp3 = new LogPearsonTypeIII(3.0, 0.5, 0.2);
```

### Probability Functions

Every distribution provides standard probability functions:

```cs
var dist = new Normal(0, 1);  // Standard normal

// Probability Density Function (PDF)
double density = dist.PDF(1.5);

// Cumulative Distribution Function (CDF)
double probability = dist.CDF(1.96);  // ≈ 0.975

// Inverse CDF (Quantile Function)
double quantile = dist.InverseCDF(0.975);  // ≈ 1.96

// Complementary CDF (Survival Function)
double exceedance = dist.CCDF(1.96);  // ≈ 0.025
```

### Distribution Properties

Access statistical properties directly:

```cs
var dist = new Normal(100, 15);

Console.WriteLine($"Mean:     {dist.Mean}");
Console.WriteLine($"Median:   {dist.Median}");
Console.WriteLine($"Mode:     {dist.Mode}");
Console.WriteLine($"Variance: {dist.Variance}");
Console.WriteLine($"Std Dev:  {dist.StandardDeviation}");
Console.WriteLine($"Skewness: {dist.Skewness}");
Console.WriteLine($"Kurtosis: {dist.Kurtosis}");
```

### Random Number Generation

Generate random samples from any distribution:

```cs
var dist = new Normal(100, 15);

// Single random value
double x = dist.InverseCDF(new Random().NextDouble());

// Multiple random values (more efficient)
double[] samples = dist.GenerateRandomValues(1000);

// With specific seed for reproducibility
double[] seededSamples = dist.GenerateRandomValues(1000, seed: 12345);
```

## Fitting Distributions to Data

### Parameter Estimation Methods

***Numerics*** supports multiple estimation methods. The simplest approach is `Estimate()`:

```cs
using Numerics.Distributions;

double[] data = { 10.2, 15.1, 12.3, 18.7, 14.2, 16.8, 13.1, 17.5 };

// L-Moments (recommended for hydrological data)
var gev = new GeneralizedExtremeValue();
gev.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);

// Maximum Likelihood Estimation
var normal = new Normal();
normal.Estimate(data, ParameterEstimationMethod.MaximumLikelihood);

// Method of Moments
var lognormal = new LogNormal();
lognormal.Estimate(data, ParameterEstimationMethod.MethodOfMoments);
```

### Hydrologic Frequency Analysis Example

This example uses annual peak streamflow from the Tippecanoe River near Delphi, Indiana (Rao & Hamed, 2000, Table 5.1.1). See [`example-data/tippecanoe-river-streamflow.csv`](example-data/tippecanoe-river-streamflow.csv).

```cs
using Numerics.Distributions;
using Numerics.Data.Statistics;

// Tippecanoe River near Delphi, IN — 48 years of annual peak streamflow (cfs)
// Source: Rao & Hamed (2000), Table 5.1.1
double[] annualPeaks = {
    6290, 2700, 13100, 16900, 14600, 9600, 7740, 8490, 8130, 12000,
    17200, 15000, 12400, 6960, 6500, 5840, 10400, 18800, 21400, 22600,
    14200, 11000, 12800, 15700, 4740, 6950, 11800, 12100, 20600, 14600,
    14600, 8900, 10600, 14200, 14100, 14100, 12500, 7530, 13400, 17600,
    13400, 19200, 16900, 15500, 14500, 21900, 10400, 7460
};

// Fit Log-Pearson Type III using L-Moments (USGS Bulletin 17C method)
var lp3 = new LogPearsonTypeIII();
lp3.Estimate(annualPeaks, ParameterEstimationMethod.MethodOfLinearMoments);

// Compute flood quantiles for key return periods
Console.WriteLine("Return Period Analysis:");
Console.WriteLine($"  10-year flood (10% AEP): {lp3.InverseCDF(0.90):N0} cfs");
Console.WriteLine($"  50-year flood (2% AEP):  {lp3.InverseCDF(0.98):N0} cfs");
Console.WriteLine($" 100-year flood (1% AEP):  {lp3.InverseCDF(0.99):N0} cfs");
Console.WriteLine($" 500-year flood (0.2% AEP): {lp3.InverseCDF(0.998):N0} cfs");

// Assess goodness-of-fit
double aic = GoodnessOfFit.AIC(lp3.NumberOfParameters, lp3.LogLikelihood(annualPeaks));
Console.WriteLine($"\nAIC: {aic:F2}");
```

## Numerical Integration

### One-Dimensional Integration

```cs
using Numerics.Mathematics.Integration;

// Define the function to integrate
Func<double, double> f = x => Math.Exp(-x * x);

// Adaptive Simpson's Rule (good general-purpose method)
var simpson = new AdaptiveSimpsonsRule(f, -5, 5);
simpson.Integrate();
Console.WriteLine($"Simpson's: {simpson.Result:F10}");

// Gauss-Kronrod (higher accuracy for smooth functions)
var gk = new AdaptiveGaussKronrod(f, -5, 5);
gk.Integrate();
Console.WriteLine($"Gauss-Kronrod: {gk.Result:F10}");
```

### Multi-Dimensional Integration

```cs
using Numerics.Mathematics.Integration;

// 3D function: f(x,y,z) = x*y*z over [0,1]³
Func<double[], double> f3d = x => x[0] * x[1] * x[2];
double[] lower = { 0, 0, 0 };
double[] upper = { 1, 1, 1 };

// Monte Carlo integration
var mc = new MonteCarloIntegration(f3d, 3, lower, upper);
mc.MaxIterations = 100000;
mc.Integrate();
Console.WriteLine($"Monte Carlo: {mc.Result:F6} ± {mc.StandardError:F6}");

// VEGAS adaptive importance sampling
Func<double[], double, double> f3dVegas = (x, w) => x[0] * x[1] * x[2];
var vegas = new Vegas(f3dVegas, 3, lower, upper);
vegas.MaxIterations = 10000;
vegas.Integrate();
Console.WriteLine($"VEGAS: {vegas.Result:F6} ± {vegas.StandardError:F6}");
```

## Optimization

### Local Optimization

```cs
using Numerics.Mathematics.Optimization;

// Rosenbrock function (classic test function)
Func<double[], double> rosenbrock = x =>
{
    double a = 1 - x[0];
    double b = x[1] - x[0] * x[0];
    return a * a + 100 * b * b;
};

// BFGS (quasi-Newton method)
var bfgs = new BFGS(rosenbrock, 2, new double[] { -1, -1 }, new double[] { -5, -5 }, new double[] { 5, 5 });
bfgs.Minimize();
Console.WriteLine($"BFGS minimum at: ({bfgs.BestParameterSet.Values[0]:F6}, {bfgs.BestParameterSet.Values[1]:F6})");
Console.WriteLine($"Function value: {bfgs.BestParameterSet.Fitness:E6}");
```

### Global Optimization

```cs
using Numerics.Mathematics.Optimization;

// Multi-modal function with many local minima
Func<double[], double> rastrigin = x =>
{
    double sum = 10 * x.Length;
    for (int i = 0; i < x.Length; i++)
        sum += x[i] * x[i] - 10 * Math.Cos(2 * Math.PI * x[i]);
    return sum;
};

double[] lower = { -5.12, -5.12 };
double[] upper = { 5.12, 5.12 };

// Differential Evolution
var de = new DifferentialEvolution(rastrigin, 2, lower, upper);
de.Minimize();
Console.WriteLine($"DE minimum at: ({de.BestParameterSet.Values[0]:F6}, {de.BestParameterSet.Values[1]:F6})");
```

## MCMC Sampling

### Basic Bayesian Inference

```cs
using Numerics.Sampling.MCMC;
using Numerics.Distributions;

// Observed data
double[] observations = { 5.2, 4.8, 5.1, 5.5, 4.9, 5.3, 5.0, 5.2 };

// Prior distributions for mean and standard deviation
var priors = new List<IUnivariateDistribution>
{
    new Normal(5, 2),       // Prior for mean: N(5, 2)
    new Uniform(0.1, 5)     // Prior for std: U(0.1, 5)
};

// Log-likelihood function
double ComputeLogLikelihood(double[] theta)
{
    double mu = theta[0];
    double sigma = theta[1];
    if (sigma <= 0) return double.NegativeInfinity;

    var model = new Normal(mu, sigma);
    return model.LogLikelihood(observations);
}

// Run DE-MCz sampler
var sampler = new DEMCz(priors, ComputeLogLikelihood);
sampler.Iterations = 20000;
sampler.WarmupIterations = 5000;
sampler.Sample();

// Analyze results - flatten chains and extract parameter values
var allSamples = sampler.Output.SelectMany(chain => chain).ToList();

for (int i = 0; i < 2; i++)
{
    var values = allSamples.Select(ps => ps.Values[i]).ToArray();
    double mean = Statistics.Mean(values);
    double lower = Statistics.Percentile(values, 0.025);
    double upper = Statistics.Percentile(values, 0.975);
    Console.WriteLine($"Parameter {i}: {mean:F3} [{lower:F3}, {upper:F3}]");
}
```

## Performance Tips

1. **Use parallelization** for large-scale Monte Carlo simulations (enabled by default; verify it is enabled):
   ```cs
   sampler.ParallelizeChains = true; // default is true
   ```

2. **Generate random values efficiently**:
   ```cs
   double[] values = dist.GenerateRandomValues(100000);
   ```

3. **Choose appropriate integration method** based on function smoothness:
   - Smooth functions → Gauss-Kronrod or Adaptive Simpson's
   - Functions with discontinuities → Adaptive methods with smaller tolerance
   - High dimensions → Monte Carlo methods

4. **Set reasonable tolerances** for iterative methods to balance accuracy and speed.

## Next Steps

- Explore the [Univariate Distributions](distributions/univariate.md) reference
- Learn about [Parameter Estimation](distributions/parameter-estimation.md) methods
- Understand [Goodness-of-Fit](statistics/goodness-of-fit.md) metrics
- Dive into [MCMC Methods](sampling/mcmc.md) for Bayesian inference

---

## References

[1] Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.
