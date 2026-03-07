# Numerics

***Numerics*** is a free and open-source numerical computing library for .NET developed by the U.S. Army Corps of Engineers Risk Management Center (USACE-RMC). It provides methods and algorithms for probability distributions, statistical analysis, numerical methods, optimization, machine learning, and Bayesian MCMC sampling — with a focus on hydrological and risk assessment applications.

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

**[User Guide and API Documentation](docs/index.md)** — Comprehensive documentation with code examples and mathematical explanations.

| Section | Topics |
|---------|--------|
| [Distributions](docs/distributions/univariate.md) | 40+ univariate distributions, parameter estimation, uncertainty analysis, copulas, multivariate distributions |
| [Statistics](docs/statistics/descriptive.md) | Descriptive statistics, goodness-of-fit metrics, hypothesis tests |
| [Data](docs/data/interpolation.md) | Interpolation, linear regression, time series analysis |
| [Machine Learning](docs/machine-learning.md) | GLM, decision trees, random forests, KNN, naive Bayes, k-means, GMM |
| [Sampling](docs/sampling/mcmc.md) | MCMC (RWMH, ARWMH, DE-MCz, HMC, NUTS, Gibbs), random generation, convergence diagnostics |
| [Mathematics](docs/mathematics/integration.md) | Integration, differentiation, optimization, root finding, linear algebra, ODE solvers, special functions |
| [References](docs/references.md) | Consolidated bibliography |

## Support

USACE-RMC is committed to maintaining and supporting the library with regular updates, bug fixes, and enhancements.

The repository includes a unit testing library with over 1,000 tests that also serve as usage examples for the classes and methods in the library.

## License

See [LICENSE](LICENSE) for details.
