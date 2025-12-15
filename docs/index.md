# Numerics Library Documentation

**Version:** 1.0  
**Author:** U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC)  
**License:** BSD-3-Clause

## Overview

***Numerics*** is a comprehensive numerical computing library for .NET, developed by the U.S. Army Corps of Engineers Risk Management Center for quantitative risk assessment applications. The library provides robust implementations of probability distributions, statistical analysis, numerical methods, optimization algorithms, and Markov Chain Monte Carlo (MCMC) sampling techniques.

The library is designed for engineers, scientists, and researchers who need reliable numerical computing capabilities, particularly in the domains of:

- Hydrological frequency analysis
- Infrastructure risk assessment
- Monte Carlo simulation
- Bayesian parameter estimation
- Statistical model validation

## Key Features

### Probability Distributions
- 40+ univariate probability distributions with PDF, CDF, and inverse CDF
- Multiple parameter estimation methods (Method of Moments, L-Moments, Maximum Likelihood)
- Uncertainty analysis via bootstrap resampling
- Bivariate copulas for dependency modeling
- Multivariate normal distribution

### Statistical Analysis
- Comprehensive goodness-of-fit metrics (NSE, KGE, RMSE, PBIAS, AIC/BIC)
- Descriptive statistics and hypothesis tests
- Autocorrelation and time series analysis
- Outlier detection (Multiple Grubbs-Beck Test)

### Numerical Methods
- Adaptive numerical integration (Simpson's Rule, Gauss-Lobatto, Gauss-Kronrod)
- Multidimensional integration (Monte Carlo, MISER, VEGAS)
- Numerical differentiation (two-point, Ridder's method)
- Root finding algorithms (Bisection, Brent, Newton-Raphson)
- ODE solvers (Runge-Kutta methods)

### Optimization
- Local optimization (BFGS, Nelder-Mead, Powell, Golden Section)
- Global optimization (Differential Evolution, SCE-UA, Simulated Annealing, Particle Swarm)
- Constrained optimization (Augmented Lagrangian)

### MCMC Sampling
- Random Walk Metropolis-Hastings (RWMH)
- Adaptive Random Walk Metropolis-Hastings (ARWMH)
- Differential Evolution MCMC (DE-MCz, DE-MCzs)
- Hamiltonian Monte Carlo (HMC)
- Gibbs sampling
- Convergence diagnostics (Gelman-Rubin, Effective Sample Size)

## Quick Start

### Creating a Distribution

```cs
using Numerics.Distributions;

// Create a Normal distribution
var normal = new Normal(100, 15);

// Compute probability functions
double pdf = normal.PDF(110);           // f(110)
double cdf = normal.CDF(110);           // P(X ≤ 110)
double quantile = normal.InverseCDF(0.95); // x such that P(X ≤ x) = 0.95

Console.WriteLine($"PDF(110) = {pdf:F6}");
Console.WriteLine($"CDF(110) = {cdf:F4}");
Console.WriteLine($"95th percentile = {quantile:F2}");
```

### Fitting a Distribution to Data

```cs
using Numerics.Distributions;

double[] annualMaxFlows = { 1200, 1500, 1100, 1800, 1350, 1600, 1250, 1450 };

// Fit using L-Moments (recommended for hydrologic data)
var gev = new GeneralizedExtremeValue();
gev.SetParameters(gev.ParametersFromLinearMoments(annualMaxFlows));

// Compute the 100-year flood (1% annual exceedance probability)
double q100 = gev.InverseCDF(0.99);
Console.WriteLine($"100-year flood estimate: {q100:F0} cfs");
```

### Numerical Integration

```cs
using Numerics.Mathematics.Integration;

// Integrate f(x) = x² from 0 to 1
Func<double, double> f = x => x * x;
var asr = new AdaptiveSimpsonsRule(f, 0, 1);
asr.Integrate();
double adaptiveResult = asr.Result;
```

### MCMC Sampling

```cs
using Numerics.Sampling.MCMC;
using Numerics.Distributions;

// Define prior distributions
var priors = new List<IUnivariateDistribution>
{
    new Normal(0, 10),    // Prior for parameter 1
    new Uniform(0, 100)   // Prior for parameter 2
};

// Define log-likelihood function
double LogLikelihood(double[] parameters)
{
    // Your likelihood calculation here
    return -0.5 * Math.Pow(parameters[0] - 5, 2);
}

// Create and run sampler
var sampler = new DEMCz(priors, LogLikelihood);
sampler.Iterations = 10000;
sampler.Sample();

// Get results
var results = sampler.Output;
```

## Documentation Structure

| Document | Description |
|----------|-------------|
| [Getting Started](getting-started.md) | Installation and basic usage patterns |
| **Distributions** | |
| [Overview](distributions/overview.md) | Distribution concepts and interfaces |
| [Univariate Distributions](distributions/univariate.md) | Complete reference for univariate distributions |
| [Parameter Estimation](distributions/parameter-estimation.md) | Fitting distributions to data |
| [Uncertainty Analysis](distributions/uncertainty-analysis.md) | Bootstrap and confidence intervals |
| [Bivariate Copulas](distributions/copulas.md) | Dependency modeling with copulas |
| **Statistics** | |
| [Descriptive Statistics](statistics/descriptive.md) | Summary statistics functions |
| [Goodness-of-Fit](statistics/goodness-of-fit.md) | Model evaluation metrics |
| **Mathematics** | |
| [Numerical Integration](mathematics/integration.md) | Quadrature methods |
| [Numerical Differentiation](mathematics/differentiation.md) | Derivative approximation |
| [Optimization](mathematics/optimization.md) | Local and global optimization |
| [Linear Algebra](mathematics/linear-algebra.md) | Matrix and vector operations |
| [Root Finding](mathematics/root-finding.md) | Equation solving methods |
| **Data** | |
| [Interpolation](data/interpolation.md) | Interpolation methods |
| [Time Series](data/time-series.md) | Time series operations |
| **Sampling** | |
| [Random Generation](sampling/random-generation.md) | RNG and quasi-random sequences |
| [MCMC Methods](sampling/mcmc.md) | Markov Chain Monte Carlo samplers |
| [References](references.md) | Complete bibliography |

## Namespaces

| Namespace | Description |
|-----------|-------------|
| `Numerics.Distributions` | Probability distributions and copulas |
| `Numerics.Data.Statistics` | Statistical functions and tests |
| `Numerics.Data.Interpolation` | Interpolation methods |
| `Numerics.Data.TimeSeries` | Time series data structures |
| `Numerics.Mathematics.Integration` | Numerical integration |
| `Numerics.Mathematics.Differentiation` | Numerical differentiation |
| `Numerics.Mathematics.Optimization` | Optimization algorithms |
| `Numerics.Mathematics.LinearAlgebra` | Matrix and vector operations |
| `Numerics.Mathematics.RootFinding` | Root finding algorithms |
| `Numerics.Mathematics.SpecialFunctions` | Gamma, Beta, Error functions |
| `Numerics.Sampling` | Random sampling and stratification |
| `Numerics.Sampling.MCMC` | MCMC samplers and diagnostics |

## Support and Contributing

This library is developed and maintained by the U.S. Army Corps of Engineers Risk Management Center. For questions, bug reports, or feature requests, please contact the development team.

## License

This software is provided under a BSD-3-Clause license. See the LICENSE file for complete terms.

---

## References

[1] Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

[2] Hosking, J. R. M. (1990). L-moments: Analysis and estimation of distributions using linear combinations of order statistics. *Journal of the Royal Statistical Society: Series B*, 52(1), 105-124.

[3] ter Braak, C. J. F., & Vrugt, J. A. (2008). Differential Evolution Markov Chain with snooker updater and fewer chains. *Statistics and Computing*, 18(4), 435-446.

[4] Moriasi, D. N., et al. (2007). Model evaluation guidelines for systematic quantification of accuracy in watershed simulations. *Transactions of the ASABE*, 50(3), 885-900.
