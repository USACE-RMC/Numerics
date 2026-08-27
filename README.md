# Numerics

[![CI](https://github.com/USACE-RMC/Numerics/actions/workflows/Integration.yml/badge.svg)](https://github.com/USACE-RMC/Numerics/actions/workflows/Integration.yml)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.10540/status.svg)](https://doi.org/10.21105/joss.10540)
[![DOI](https://zenodo.org/badge/697446868.svg)](https://zenodo.org/badge/latestdoi/697446868)
[![NuGet](https://img.shields.io/nuget/v/RMC.Numerics)](https://www.nuget.org/packages/RMC.Numerics/)
[![License: 0BSD](https://img.shields.io/badge/License-0BSD-blue.svg)](LICENSE)

***Numerics*** is a free and open-source numerical computing library for .NET developed by the U.S. Army Corps of Engineers Risk Management Center (USACE-RMC). It provides methods and algorithms for probability distributions, statistical analysis, numerical methods, optimization, machine learning, and Bayesian MCMC sampling - with a focus on hydrological and risk assessment applications.

## Supported Frameworks

| Framework | Version |
|-----------|---------|
| .NET | 10.0, 9.0, 8.0 |
| .NET Framework | 4.8.1 |

Install via NuGet:
```
dotnet add package RMC.Numerics
```
Or search for [RMC.Numerics](https://www.nuget.org/packages/RMC.Numerics/) in the NuGet Package Manager.

## Documentation

**[User Guide and API Documentation](docs/index.md)** - Comprehensive documentation with code examples and mathematical explanations.

| Section | Topics |
|---------|--------|
| [Mathematics](docs/mathematics/integration.md) | Adaptive one- and two-dimensional integration, differentiation, optimization, root finding, linear algebra, ODE solvers, special functions |
| [Data](docs/data/interpolation.md) | Interpolation, linear regression, time series analysis |
| [Statistics](docs/statistics/descriptive.md) | Descriptive and weighted statistics, goodness-of-fit metrics, hypothesis tests, global sensitivity analysis |
| [Distributions](docs/distributions/univariate.md) | 43 univariate distributions, parameter estimation, uncertainty analysis, copulas, multivariate distributions |
| [Machine Learning](docs/machine-learning/machine-learning.md) | GLM, decision trees, random forests, KNN, naive Bayes, k-means, GMM |
| [Sampling](docs/sampling/mcmc.md) | MCMC (RWMH, ARWMH, DE-MCz, HMC, NUTS, Gibbs), random generation, scrambled quasi-random sequences, convergence diagnostics |
| [References](docs/references.md) | Consolidated bibliography |

## Prerequisites
- .NET 8+ runtime (or .NET Framework 4.8.1 on Windows). Install the [.NET SDK](https://dotnet.microsoft.com/download) if you don't already have it.
  - The Microsoft .NET SDK is available on the App Portal for Corps users.

## Support

USACE-RMC is committed to maintaining and supporting the library with regular updates, bug fixes, and enhancements.

The repository includes a unit testing library with over 2,300 tests that also serve as usage examples for the classes and methods in the library.

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on bug reports, feature requests, and pull requests.

## License

See [LICENSE](LICENSE) for details.
