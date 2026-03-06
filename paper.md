---
title: 'Numerics: A .NET Library for Numerical Computing, Statistical Analysis, and Risk Assessment'
tags:
  - C#
  - .NET
  - numerical methods
  - statistics
  - probability distributions
  - Bayesian inference
  - hydrology
  - risk assessment
authors:
  - name: Haden Smith
    orcid: 0000-0000-0000-0000
    affiliation: 1
    corresponding: true
  - name: Woodrow Fields
    affiliation: 1
  - name: Brennan Beam
    affiliation: 1
  - name: Julian Gonzalez
    affiliation: 1
affiliations:
  - name: U.S. Army Corps of Engineers, Risk Management Center, Lakewood, Colorado, USA
    index: 1
date: 31 January 2026
bibliography: paper.bib
---

# Summary

Numerics is a free and open-source numerical computing library for .NET, developed by the U.S. Army Corps of Engineers Risk Management Center (USACE-RMC). The library provides comprehensive implementations of probability distributions, statistical analysis methods, numerical algorithms, optimization routines, and Markov Chain Monte Carlo (MCMC) sampling techniques. Designed for engineers, scientists, and researchers requiring reliable numerical computing capabilities, Numerics is particularly focused on hydrological frequency analysis, infrastructure risk assessment, Monte Carlo simulation, and Bayesian parameter estimation.

The library implements over 40 univariate probability distributions with multiple parameter estimation methods including L-moments [@hosking1990], maximum likelihood, and method of moments. It provides robust numerical methods for integration, differentiation, root finding, and optimization, alongside six MCMC samplers with convergence diagnostics. Numerics serves as the computational foundation for several USACE-RMC software applications used for dam and levee safety risk analysis across the United States.

# Statement of Need

Infrastructure risk assessment—particularly for dams and levees—requires sophisticated probabilistic analysis that integrates extreme value statistics, uncertainty quantification, and Monte Carlo simulation. Analysts in this domain face unique challenges: small sample sizes from historical records, the need for robust outlier handling, regulatory requirements for specific statistical methods (such as USGS Bulletin 17C for flood frequency analysis), and the integration of expert judgment with observational data through Bayesian inference.

While the Python ecosystem offers mature numerical libraries like NumPy [@harris2020] and SciPy [@virtanen2020], and the R ecosystem provides specialized hydrology packages, the .NET ecosystem—which dominates enterprise and government software development—lacks equivalent comprehensive libraries tailored for engineering risk analysis. Existing .NET numerical libraries such as Math.NET Numerics provide general-purpose capabilities but do not include hydrology-specific distributions (e.g., Log-Pearson Type III), L-moment parameter estimation, specialized extreme value distributions (Kappa-Four, Generalized Logistic), or the suite of MCMC samplers necessary for Bayesian uncertainty quantification.

Numerics addresses this gap by providing:

1. **Domain-specific probability distributions**: The Log-Pearson Type III distribution (mandated by USGS Bulletin 17C for U.S. flood frequency analysis [@england2019]), Generalized Extreme Value, Pearson Type III, and other distributions commonly used in hydrological and reliability analysis.

2. **L-moment parameter estimation**: L-moments provide more robust and efficient parameter estimates than conventional moments for small samples typical in hydrological records [@hosking1997], yet few libraries implement them comprehensively.

3. **Integrated uncertainty quantification**: Bootstrap resampling methods combined with analytical standard errors enable practitioners to report confidence intervals on design estimates—a requirement in dam safety risk assessments.

4. **Bayesian inference capabilities**: Multiple MCMC samplers (Random Walk Metropolis-Hastings, Adaptive RWMH, Differential Evolution MCMC, Hamiltonian Monte Carlo, and Gibbs sampling) with convergence diagnostics allow analysts to incorporate prior knowledge and quantify parametric uncertainty.

5. **Global optimization algorithms**: Infrastructure risk models often require calibrating complex, multi-modal objective functions. Numerics provides Differential Evolution [@storn1997], Shuffled Complex Evolution (SCE-UA) [@duan1994], Particle Swarm Optimization, and other global optimizers commonly used in hydrological model calibration.

Numerics serves as the computational engine for multiple USACE-RMC software applications:

- **RMC-BestFit**: Distribution fitting and frequency analysis with uncertainty quantification
- **RMC-RFA**: Regional frequency analysis for ungauged sites
- **RMC-TotalRisk**: Comprehensive dam and levee risk quantification
- **LifeSim**: Life loss estimation for dam and levee failures
- **Levee Screening Tool (LST)**: Portfolio-level levee risk screening
- **Dam Screening Tool (DST)**: Portfolio-level dam risk screening

These tools are used by USACE and other federal agencies for infrastructure safety decisions affecting millions of Americans. By open-sourcing Numerics, USACE-RMC enables independent verification, community contributions, and adoption by the broader water resources engineering community.

# Key Features

## Probability Distributions

Numerics implements 40+ univariate probability distributions organized into continuous and discrete categories. Each distribution provides:

- Probability density function (PDF) and log-PDF
- Cumulative distribution function (CDF) and survival function
- Inverse CDF (quantile function) for computing design values
- Random variate generation
- Parameter estimation via multiple methods
- Theoretical moments and L-moments

Distributions particularly relevant for extreme value analysis include:

| Distribution | Parameters | Primary Application |
|-------------|-----------|---------------------|
| Log-Pearson Type III | μ, σ, γ | USGS flood frequency (Bulletin 17C) |
| Generalized Extreme Value | ξ, α, κ | Block maxima analysis |
| Generalized Pareto | ξ, α, κ | Peaks-over-threshold |
| Kappa-Four | ξ, α, κ, h | Flexible extreme value |
| Pearson Type III | μ, σ, γ | Precipitation, flood flows |

The library also provides bivariate copulas (Clayton, Frank, Gumbel, Joe, Normal) for modeling dependence structures between variables, essential for multivariate risk analysis [@nelsen2006; @salvadori2007].

## Parameter Estimation

Numerics supports four parameter estimation methods:

```csharp
// L-moments: recommended for small hydrological samples
var gev = new GeneralizedExtremeValue();
gev.Estimate(annualPeaks, ParameterEstimationMethod.MethodOfLinearMoments);

// Maximum likelihood for large samples
gev.Estimate(annualPeaks, ParameterEstimationMethod.MaximumLikelihood);
```

L-moment estimation follows the algorithms of @hosking1990 and provides superior performance for small samples (n < 50) typical in flood frequency analysis. The Multiple Grubbs-Beck Test [@cohn2013] is implemented for detecting potentially influential low outliers in flood series.

## MCMC Sampling

Six MCMC samplers enable Bayesian inference:

| Sampler | Description | Use Case |
|---------|-------------|----------|
| RWMH | Random Walk Metropolis-Hastings | Baseline sampler |
| ARWMH | Adaptive RWMH [@haario2001] | Self-tuning proposal |
| DE-MCz | Differential Evolution MCMC [@terbraak2008] | High-dimensional problems |
| DE-MCzs | DE-MC with snooker update | Multimodal posteriors |
| HMC | Hamiltonian Monte Carlo [@neal2011] | Efficient gradient-based |
| Gibbs | Gibbs sampling | Conditional conjugacy |

Convergence diagnostics include the Gelman-Rubin statistic (R-hat) [@vehtari2021] and effective sample size.

## Numerical Methods

The library provides production-quality implementations of:

- **Integration**: Adaptive Simpson's rule, Gauss-Kronrod, Gauss-Lobatto, Monte Carlo (including VEGAS and MISER for multidimensional integration)
- **Optimization**: BFGS, Nelder-Mead, Differential Evolution, SCE-UA, Particle Swarm, Simulated Annealing
- **Root finding**: Bisection, Brent's method, Newton-Raphson, Secant
- **Linear algebra**: LU, QR, Cholesky, SVD, and eigenvalue decompositions
- **ODE solvers**: Runge-Kutta 4th order and adaptive RK4/5 (Dormand-Prince)

## Data Analysis

Time series capabilities include support for irregular and regular intervals, multiple temporal scales, moving averages, and direct data download from USGS water services. Interpolation methods (linear, cubic spline, polynomial, bilinear) support typical engineering workflows. Comprehensive goodness-of-fit metrics include Nash-Sutcliffe Efficiency [@nash1970], Kling-Gupta Efficiency [@gupta2009], and information criteria (AIC, BIC) for model selection [@burnham2002].

# Example Usage

A typical flood frequency analysis workflow:

```csharp
using Numerics.Distributions;
using Numerics.Data.Statistics;

// Annual peak flow data (cfs)
double[] annualPeaks = { 12500, 15300, 11200, 18700, 14100,
                          16800, 13400, 17200, 10500, 19300 };

// Fit Log-Pearson Type III using L-moments
var lp3 = new LogPearsonTypeIII();
lp3.Estimate(annualPeaks, ParameterEstimationMethod.MethodOfLinearMoments);

// Compute design floods with uncertainty
int[] returnPeriods = { 10, 50, 100, 500 };
foreach (int T in returnPeriods)
{
    double aep = 1.0 / T;
    double designFlood = lp3.InverseCDF(1 - aep);
    Console.WriteLine($"{T}-year flood: {designFlood:N0} cfs");
}

// Bootstrap confidence intervals
var bootstrap = lp3.Bootstrap(annualPeaks, 1000);
double q100_lower = bootstrap.Quantile(0.99, 0.05);  // 5th percentile
double q100_upper = bootstrap.Quantile(0.99, 0.95); // 95th percentile
```

# Testing and Quality Assurance

Numerics includes over 1,000 unit tests organized across 145 test classes, validating numerical accuracy against published references including @press2007, @rao2000, and official USGS test cases from Bulletin 17C. Tests cover:

- Distribution PDF, CDF, and quantile function accuracy
- Parameter estimation correctness
- Numerical method convergence and error bounds
- MCMC sampler mixing and convergence
- Edge cases and numerical stability

The library supports .NET 9.0, .NET 8.0, and .NET Framework 4.8.1, ensuring compatibility with both modern and legacy enterprise systems. Continuous integration via GitHub Actions runs the test suite on each pull request.

# Documentation

Comprehensive documentation includes 24 markdown files covering all major features with mathematical explanations and code examples. The documentation follows the library's namespace structure:

- Mathematics: Integration, differentiation, optimization, linear algebra, root finding, ODE solvers
- Distributions: 40+ distributions, parameter estimation, uncertainty analysis, copulas
- Statistics: Descriptive statistics, goodness-of-fit, hypothesis tests
- Sampling: MCMC methods, convergence diagnostics, quasi-random sequences
- Data: Time series, interpolation

A consolidated bibliography provides 64 references to foundational literature in numerical analysis, statistics, hydrology, and probability theory.

# Acknowledgments

The authors acknowledge the U.S. Army Corps of Engineers Risk Management Center for supporting the development and open-source release of this software. The library benefits from foundational algorithms described in @press2007 and statistical methods developed by @hosking1990 and @hosking1997.

# References
