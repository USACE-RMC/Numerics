# Numerics Port Review — Fix Summary

Thank you for the detailed C++ port notes. We independently checked the reported issues against the C# source and corrected the items that could affect numerical accuracy, reliability, or API behavior.

## What we fixed

- Corrected Student-t density scaling and removed false saturation in extreme inverse-CDF tails.
- Stabilized Pearson III, Log-Pearson III, Gamma-gradient, and Generalized Logistic calculations near zero skew/shape and at large parameters.
- Corrected Beta-family boundary modes while preserving degenerate PERT behavior.
- Repaired the Genz multivariate-normal singular-covariance branch and clarified the original Fortran `MVNDNT` status initialization without changing its mathematics.
- Fixed descending bisection/hunt searches and made Brent bracketing use its geometric expansion factor with bounded failure behavior.
- Completed univariate factory mappings and added explicit exceptions for unsupported distribution types.
- Fixed stale or mutated state in Competing Risks, Bivariate Empirical interpolation, and related MVN caching.
- Restored Histogram adaptive-range behavior, Bilinear guarded logarithms, and the HPCM extreme-tail underflow guard.
- Corrected MCMC MAP fitness sign handling and ensured NUTS uses custom gradients during initialization as well as trajectory integration.
- Fixed the guarded `Tools.Log10` implementation. We also verified `Tools.ParallelAdd` is race-safe, documented its floating-point ordering limitation, and added concurrent and `NaN` regression tests.

Items such as the noncentral-t AGK integration note were reviewed and retained where the C# behavior was intentional. No public signatures were removed or changed.

## Validation

- Added focused regression coverage for every behavior change, including analytical and extreme-value cases.
- Release build completes with **zero errors and zero warnings**.
- **1,875/1,875 tests pass** on each of .NET Framework 4.8.1, .NET 8, .NET 9, and .NET 10.
- Brent bracketing for a distant minimum improved from **10,002 to 16 objective evaluations**; the measured Powell median improved from **7 ms to 5 ms**.

The remediation is contained in commits `33dc1af` and `651035e`.
