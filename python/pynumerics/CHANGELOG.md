# Changelog

All notable changes to pynumerics. The pynumerics version equals the RMC.Numerics version of the
bundled `Numerics.dll`; the RMC.Numerics changes of each version are in its
[GitHub release notes](https://github.com/USACE-RMC/Numerics/releases).

## 2.3.0 (2026-11-10)

The first release of pynumerics on PyPI. pynumerics replaces the `rmc-numerics` package that was
used with the RMC.Numerics 2.3 pre-release builds and was never published. Code written against
`rmc-numerics` imports `pynumerics` instead; there is no alias module (see "Migrating from
rmc-numerics" in the getting-started guide).

### What the package contains

- `Numerics.dll` from `lib/net10.0` of the RMC.Numerics 2.3.0 NuGet package, byte for byte. The
  wheel records the package version, the assembly version and the SHA-256 of the file;
  `numerics_identity()` reports the loaded assembly and `assert_numerics_identity()` checks it.
- One CoreCLR runtime and one `Numerics` assembly per process (`runtime`). Importing the package
  does not start .NET; the first wrapper call does, or `initialize()`. `PYNUMERICS_DLL_DIR`
  selects a different `Numerics.dll`, and `add_assembly()` loads dependent assemblies.
- Bulk, bit-exact copies between NumPy and .NET arrays (`marshal`).
- Python exceptions that carry the Numerics messages, raised by every wrapper (`errors`).
- Callbacks from Python callables or from C# source compiled with the bundled Roslyn 4.14.0
  compiler (`callbacks.compile_csharp`).
- Typed wrappers with NumPy-style docstrings for every Numerics area, and a user guide with
  runnable examples for each area.

### Requirements

- CPython 3.11 to 3.14 on Linux or Windows.
- The .NET 10 runtime (`Microsoft.NETCore.App` 10.0). Set `DOTNET_ROOT` when .NET is not installed
  in the default location.
- pythonnet 3.1 and NumPy 1.26 or later, installed as dependencies.

### Module layout

Each name lives only in the module that mirrors its Numerics namespace.

| Module | Contents |
|---|---|
| `distributions`, `multivariate`, `copulas` | Univariate and multivariate distributions, bivariate copulas, and the per-element Poisson and Binomial array functions |
| `regional` | Hosking-Wallis regional frequency analysis |
| `statistics`, `hypothesis_tests`, `goodness_of_fit`, `probability` | Descriptive and column statistics, hypothesis tests, goodness-of-fit measures, and probability combination rules |
| `sampling`, `mcmc`, `bootstrap` | Random generation and sampling designs, the MCMC samplers and convergence diagnostics, and the bootstrap |
| `integration`, `differentiation`, `roots`, `ode` | Numerical integration, derivatives, root finding, and ODE solvers |
| `optimization`, `graphs` | Local, global and multi-objective optimizers, dynamic programming and knapsack, and network shortest paths |
| `special`, `linalg`, `mathematics` | Special functions, matrices and decompositions, and the `Numerics.Tools` utilities with Fourier transforms |
| `functions`, `data`, `timeseries`, `download` | Univariate and link functions, paired data and interpolation, time series, and USGS, GHCN-Daily, Water Survey of Canada and BOM downloads with offline parsers |
| `machine_learning` | Regression, classification, clustering, and dimension reduction |
| `enums` | The Numerics enums the wrappers accept |
| `runtime`, `marshal`, `identity`, `errors`, `callbacks` | The loader, array marshalling, assembly identity, exception translation, and delegates |

### Notable APIs

- `regional.RegionalData`, `RegionalHomogeneityTest` and `RegionalGrowthCurve`.
- `machine_learning.GaussianProcessRegression`, `Kriging`, `ElasticNetRegression` and
  `PrincipalComponentAnalysis`, alongside `GeneralizedLinearModel`, `RandomForest`, `KMeans`,
  `GaussianMixtureModel` and the other learners.
- `optimization.LBFGSB`, `CMAES`, `DynamicallyDimensionedSearch`, `ParetoFront`,
  `DynamicProgramming` and `knapsack`; every optimizer takes a progress callback and a
  cancellation token.
- `mcmc.NUTS`, `HMC`, `DEMCzs` and the other samplers, with `gelman_rubin` and
  `effective_sample_size`.
- `linalg.EigenValueDecomposition` with `EigenSolverMethod.QL`, and
  `integration.integrate_adaptive_gauss_kronrod` with `AdaptiveGaussKronrodMethod.QAGS`.
- `graphs.shortest_distances` for multi-source shortest distances on a sparse graph.
- `statistics.column_means`, `column_variances`, `column_percentiles` and
  `column_linear_moments`, and `goodness_of_fit.column_kolmogorov_smirnov_uniform` and
  `column_anderson_darling_uniform`.
- `distributions.poisson_inverse_cdf`, `poisson_log_pmf`, `poisson_cdf`, `binomial_log_pmf` and
  `binomial_cdf` for per-element arrays.
- `mathematics.log_sum_exp` with an `axis` argument.
- `download.from_usgs` (USGS Water Data API) and `parse_usgs`, `parse_ghcn`, `parse_chmn` and
  `parse_abom` for saved responses.

### Release process

Each release is built from the RMC.Numerics package on nuget.org and published by
`.github/workflows/PythonPublish.yml` after the NuGet release, with trusted publishing. The wheel
is tested on Linux and Windows before it is uploaded.
