# Numerics

[![CI](https://github.com/USACE-RMC/Numerics/actions/workflows/Integration.yml/badge.svg)](https://github.com/USACE-RMC/Numerics/actions/workflows/Integration.yml)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.10540/status.svg)](https://doi.org/10.21105/joss.10540)
[![DOI](https://zenodo.org/badge/697446868.svg)](https://zenodo.org/badge/latestdoi/697446868)
[![NuGet](https://img.shields.io/nuget/v/RMC.Numerics)](https://www.nuget.org/packages/RMC.Numerics/)
[![PyPI](https://img.shields.io/pypi/v/pynumerics)](https://pypi.org/project/pynumerics/)
[![License: 0BSD](https://img.shields.io/badge/License-0BSD-blue.svg)](LICENSE)

***Numerics*** is a free and open-source numerical computing library for .NET developed by the U.S. Army Corps of Engineers Risk Management Center (USACE-RMC). It provides methods and algorithms for probability distributions, statistical analysis, numerical methods, optimization, machine learning, and Bayesian MCMC sampling, with a focus on hydrological and risk assessment applications.

## Supported Frameworks

| Framework | Version |
|-----------|---------|
| .NET | 10.0 |
| .NET Framework | 4.8.1 |

The .NET 10 build runs on Windows and Linux. Version 2.3.0 targets `net10.0` and `net481`; applications on .NET 8 or .NET 9 can stay on RMC.Numerics 2.2.0.

Install via NuGet:
```bash
dotnet add package RMC.Numerics --version 2.3.0
```
Or search for [RMC.Numerics](https://www.nuget.org/packages/RMC.Numerics/) in the NuGet Package Manager.

### Python

[pynumerics](python/pynumerics/README.md) ships the same `Numerics.dll` with typed Python wrappers over every Numerics area. It requires CPython 3.11 to 3.14 on Linux or Windows and the .NET 10 runtime:
```bash
pip install pynumerics
```

## What's New in 2.3.0

- Hosking-Wallis regional frequency analysis: discordancy, heterogeneity and goodness-of-fit measures, and regional growth curves.
- Machine learning: Gaussian process regression, kriging, elastic net regression, and principal component analysis.
- Optimization: L-BFGS-B, CMA-ES, dynamically dimensioned search, Pareto fronts, dynamic programming, and knapsack solvers.
- Numerical methods: a QL eigensolver, an opt-in QUADPACK QAGS integrator, multi-source shortest paths, column statistics, and compensated sums.
- Accuracy fixes in the Student-t distribution, the incomplete beta and inverse incomplete gamma functions, the Archimedean copula densities, and the GEV, Gamma and Pearson Type III L-moments.
- The pynumerics Python package, published to PyPI with each release.

See the [v2.3.0 release notes](https://github.com/USACE-RMC/Numerics/releases/tag/v2.3.0) for the full list, including breaking changes and the results that move.

## Documentation

**[User Guide and API Documentation](docs/index.md)** - Comprehensive documentation with code examples and mathematical explanations.

| Section | Topics |
|---------|--------|
| [Mathematics](docs/mathematics/integration.md) | Adaptive one- and two-dimensional integration, differentiation, optimization, root finding, linear algebra, ODE solvers, special functions, Fourier methods |
| [Data](docs/data/interpolation.md) | Interpolation, linear regression, time series analysis |
| [Statistics](docs/statistics/descriptive.md) | Descriptive and weighted statistics, goodness-of-fit metrics, hypothesis tests, global sensitivity analysis |
| [Distributions](docs/distributions/univariate.md) | 43 univariate distributions, parameter estimation, uncertainty analysis, copulas, multivariate distributions, regional frequency analysis |
| [Functions](docs/functions/index.md) | Uncertain univariate function forms, composition, posterior ensembles, link functions |
| [Machine Learning](docs/machine-learning/machine-learning.md) | GLM, elastic net, Gaussian process regression, kriging, decision trees, random forests, KNN, naive Bayes, k-means, GMM, PCA |
| [Sampling](docs/sampling/mcmc.md) | MCMC (RWMH, ARWMH, DE-MCz, DE-MCzs, HMC, NUTS, Gibbs), random generation, scrambled quasi-random sequences, convergence diagnostics |
| [References](docs/references.md) | Consolidated bibliography |

## Prerequisites
- .NET 10 runtime on Windows or Linux, or .NET Framework 4.8.1 on Windows. Install the [.NET SDK](https://dotnet.microsoft.com/download) if you don't already have it.
  - The Microsoft .NET SDK is available on the App Portal for Corps users.

## Support

USACE-RMC is committed to maintaining and supporting the library with regular updates, bug fixes, and enhancements.

The repository includes a unit testing library with over 4,000 tests that also serve as usage examples for the classes and methods in the library, and the pynumerics package adds over 5,000 Python tests.

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on bug reports, feature requests, and pull requests.

## AI Disclosure

Frontier generative AI models and coding agents are used to assist with software development, documentation, testing, and verification workflows. Project maintainers remain responsible for reviewing contributions, validating scientific methods and results, and deciding what is accepted or released.

## License

See [LICENSE](LICENSE) for details.
