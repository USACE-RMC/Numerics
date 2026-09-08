# Independent generalized-family Fisher information oracle

Generated with R version 4.4.3 (2025-02-28 ucrt) on x86_64-w64-mingw32 .

The generator uses base R only. It does not load Numerics or any production numerical helper.
All CSV matrices use location=0, scale=1, sample_size=1 and parameter order xi, alpha, kappa[, hondo].
GNO denotes Hosking's transformed normal, not the symmetric exponential-power distribution.
GLO fixes Kappa Four hondo=-1; GEV fixes hondo=0. A fixed-shape family inverts its own information block.

## Method and acceptance checks

R integrate (QUADPACK): relative tolerance 2e-12; absolute tolerance 2e-13; at most 1500 subdivisions.
GNO uses direct integration over the entire normal latent line, with an independent probability-coordinate cross-check.
GLO, GEV and K4 pair p=u^a/2 and 1-p=u^a/2 over the complete u interval [0,1]. Powers 16/24 are cross-checked; boundary-approach cases use 128/192.
Log probabilities and log survival probabilities are retained separately. Below exp(-36), log(t)=log(survival) is used to binary64 relative precision; no probability tail is discarded.
Analytical scores are checked against five-point fixed-observation derivatives of an independently expressed log density at probabilities .1, .5 and .9, using two step sizes.
Finite matrices, near-zero mean scores, positive eigenvalues, inverse residuals and agreement between integration coordinates are required before writing the fixture.
The covariance error column is a conservative matrix-norm perturbation estimate from the numerical integration error estimates; it is not a rigorous interval bound.
The alternate-coordinate difference is independent numerical corroboration, not a rigorous bound. Values are binary64 R references, not arbitrary-precision goldens.
The exact GNO kappa=0 covariance is [[7/6,0,1/3],[0,1/2,0],[1/3,0,2/3]] with all three parameters estimated.

The additional moment rows use row=1 mean, row=2 standard deviation, row=3 skewness, row=4 non-excess kurtosis. GNO uses analytical lognormal moments; GLO uses trigonometric moments, with factored coefficient-array series for abs(kappa)<=.01.
Moment absolute_error_estimate is an arithmetic comparison allowance of 5e-13*max(1,abs(value)), not a quadrature error estimate. The GLO near-zero series is cross-checked against independent 70-digit Decimal values at kappa=.00010001 and .0002.

## Mathematical domain

GNO information is finite for every finite kappa. GLO uses abs(kappa)<1/2. K4 uses kappa<1/2, hondo<1/2 and kappa*hondo<1/2. GEV uses kappa<1/2.
These are local asymptotic-information conditions, not guarantees of global-MLE existence, estimator convergence, or finite-sample coverage.
MLE covariance is not L-moment or product-moment covariance. The fixture does not authorize substitution between estimators.

## Primary-source mapping and references

- [Hosking lmom R defining quantiles](https://raw.githubusercontent.com/cran/lmom/master/R/lmom.r): quagno(p,c(xi,alpha,k)), quaglo(p,c(xi,alpha,k)), quakap(p,c(xi,alpha,k,h)), quagev(p,c(xi,alpha,k)).
- [Park and Kim (2007), Fisher information matrix for a four-parameter kappa distribution](https://doi.org/10.1016/j.spl.2007.03.002). The present scores/domain are independently derived; inaccessible full-paper formulae were not treated as verified values.
- [Wang and Flournoy (2015), local likelihood estimation for the three-parameter lognormal](https://doi.org/10.1016/j.spl.2015.05.021): finite Fisher information does not imply bounded global likelihood.

## Numerical summary

| Family | kappa | hondo | Minimum eigenvalue | Condition number | Max mean score | Information error norm | Coordinate delta norm | Covariance error estimate | Fixed-x derivative relative error |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GNO | 0 | NA | 0.75 | 2.66667 | 2.88e-17 | 2.05e-12 | 1.48e-15 | 3.65e-12 | 2.84e-12 |
| GNO | -1e-06 | NA | 0.75 | 2.66667 | 9.17e-17 | 2.05e-12 | 1.3e-15 | 3.65e-12 | 1.3e-11 |
| GNO | 1e-06 | NA | 0.75 | 2.66667 | 9.17e-17 | 2.05e-12 | 1.38e-15 | 3.65e-12 | 3.87e-12 |
| GNO | -0.2 | NA | 0.742867 | 3.98025 | 1.34e-16 | 3.56e-12 | 1.09e-15 | 6.45e-12 | 4.91e-12 |
| GNO | 0.2 | NA | 0.742867 | 3.98025 | 1.34e-16 | 3.56e-12 | 1.09e-15 | 6.45e-12 | 3.5e-12 |
| GNO | -1 | NA | 0.422625 | 79.9633 | 5.42e-17 | 4.83e-11 | 2.93e-14 | 2.71e-10 | 5.08e-12 |
| GNO | 1 | NA | 0.422625 | 79.9633 | 5.42e-17 | 4.83e-11 | 2.93e-14 | 2.71e-10 | 2.47e-12 |
| GNO | -2 | NA | 0.177334 | 110247 | 1.11e-16 | 1.54e-08 | 3.14e-11 | 4.89e-07 | 3.47e-12 |
| GNO | 2 | NA | 0.177334 | 110247 | 1.11e-16 | 1.54e-08 | 3.14e-11 | 4.89e-07 | 3.5e-12 |
| GLO | 0 | -1 | 0.320008 | 11.8793 | 1.11e-16 | 5.58e-14 | 9.32e-16 | 5.45e-13 | 7.35e-12 |
| GLO | -1e-06 | -1 | 0.320008 | 11.8793 | 2.01e-16 | 5.58e-14 | 1.13e-15 | 5.45e-13 | 8.1e-12 |
| GLO | 1e-06 | -1 | 0.320008 | 11.8793 | 2.78e-17 | 5.58e-14 | 9.35e-16 | 5.45e-13 | 6.36e-12 |
| GLO | -0.2 | -1 | 0.309535 | 18.5564 | 8.33e-17 | 6.61e-13 | 2.37e-15 | 6.9e-12 | 7.92e-12 |
| GLO | 0.2 | -1 | 0.309535 | 18.5564 | 4.37e-16 | 6.61e-13 | 2.49e-15 | 6.9e-12 | 7.31e-12 |
| GLO | -0.4 | -1 | 0.276988 | 78.6219 | 6.61e-16 | 3.37e-11 | 4.96e-14 | 4.39e-10 | 7.97e-12 |
| GLO | 0.4 | -1 | 0.276988 | 78.6219 | 7.98e-17 | 3.37e-11 | 4.84e-14 | 4.39e-10 | 5.35e-12 |
| GLO | -0.45 | -1 | 0.266216 | 178.618 | 9.71e-17 | 9.57e-12 | 1.67e-14 | 1.35e-10 | 5.44e-12 |
| GLO | 0.45 | -1 | 0.266216 | 178.618 | 4.65e-16 | 9.58e-12 | 1.52e-14 | 1.35e-10 | 5.37e-12 |
| GLO | -0.49 | -1 | 0.257239 | 999.511 | 3.07e-16 | 2.07e-10 | 3.97e-13 | 3.13e-09 | 4.25e-12 |
| GLO | 0.49 | -1 | 0.257239 | 999.511 | 1.39e-16 | 2.07e-10 | 3.04e-13 | 3.13e-09 | 6.35e-12 |
| GEV | 0 | 0 | 0.672476 | 3.87267 | 1.11e-16 | 1.1e-12 | 4.78e-16 | 2.43e-12 | 7.93e-12 |
| GEV | -1e-06 | 0 | 0.672476 | 3.87266 | 1.67e-16 | 1.16e-12 | 1.01e-15 | 2.56e-12 | 5.94e-12 |
| GEV | 1e-06 | 0 | 0.672477 | 3.87268 | 7.29e-17 | 1.16e-12 | 5.04e-16 | 2.56e-12 | 5.56e-12 |
| GEV | -0.2 | 0 | 0.58145 | 4.69117 | 1.11e-16 | 7.71e-13 | 2.01e-15 | 2.28e-12 | 6.77e-12 |
| GEV | 0.2 | 0 | 0.779251 | 6.89922 | 6.38e-16 | 3.32e-12 | 2.63e-16 | 5.47e-12 | 6.48e-12 |
| GEV | -1 | 0 | 0.292005 | 48.1418 | 3.47e-16 | 1.46e-12 | 2.59e-15 | 1.71e-11 | 1.24e-11 |
| GEV | 0.45 | 0 | 0.822612 | 57.7779 | 1.39e-16 | 9.02e-12 | 2.49e-14 | 1.33e-11 | 4.14e-12 |
| GEV | 0.49 | 0 | 0.804731 | 319.524 | 2.22e-16 | 1.83e-10 | 2.81e-13 | 2.82e-10 | 5.09e-12 |
| K4 | 0 | 0 | 0.0696442 | 56.6234 | 1.25e-16 | 1.1e-12 | 7.19e-16 | 2.27e-10 | 7.93e-12 |
| K4 | 1e-06 | 0 | 0.0696443 | 56.6233 | 1.25e-16 | 1.16e-12 | 5.51e-15 | 2.39e-10 | 5.56e-12 |
| K4 | -1e-06 | 0 | 0.069644 | 56.6236 | 1.67e-16 | 1.16e-12 | 1.02e-15 | 2.39e-10 | 5.94e-12 |
| K4 | 0 | 1e-06 | 0.0696443 | 56.6234 | 3.05e-16 | 1.1e-12 | 1.22e-15 | 2.27e-10 | 4.28e-12 |
| K4 | 0 | -1e-06 | 0.069644 | 56.6235 | 3.33e-16 | 1.1e-12 | 2.3e-15 | 2.27e-10 | 8.2e-12 |
| K4 | 0 | -1 | 0.00507523 | 770.598 | 1.11e-16 | 1.45e-13 | 9.33e-16 | 5.63e-09 | 7.35e-12 |
| K4 | 0.2 | -1 | 0.0155506 | 370.485 | 4.37e-16 | 6.7e-13 | 2.58e-15 | 2.77e-09 | 8.15e-12 |
| K4 | -0.2 | -1 | 0.00824959 | 719.157 | 8.33e-17 | 6.78e-13 | 2.46e-15 | 9.96e-09 | 8.89e-12 |
| K4 | 0.2 | 0.2 | 0.139053 | 45.1301 | 2.91e-16 | 2.75e-12 | 4.88e-15 | 1.42e-10 | 7.47e-12 |
| K4 | -0.2 | -0.2 | 0.0212752 | 171.195 | 1.53e-16 | 3.11e-13 | 2.05e-15 | 6.88e-10 | 1.07e-11 |
| K4 | 0.4 | -1 | 0.0328533 | 662.914 | 7.98e-17 | 3.37e-11 | 4.84e-14 | 3.12e-08 | 5.35e-12 |
| K4 | -0.4 | -1 | 0.0336709 | 650.134 | 6.61e-16 | 3.37e-11 | 4.96e-14 | 2.97e-08 | 8.89e-12 |
| K4 | -1 | -0.2 | 0.0174094 | 1205.62 | 7.22e-16 | 3.92e-12 | 4.99e-15 | 1.29e-08 | 6.07e-12 |
| K4 | -0.2 | -2 | 0.0160321 | 3193.35 | 1.85e-15 | 8.67e-11 | 7.51e-14 | 3.38e-07 | 1.11e-11 |
| K4 | 0.49 | 0.2 | 0.18099 | 1420.86 | 1.67e-16 | 1.91e-10 | 3.15e-13 | 5.83e-09 | 4.44e-12 |
| K4 | 0.2 | 0.49 | 0.242565 | 1060.47 | 8.05e-16 | 1.68e-10 | 5e-13 | 2.85e-09 | 4.31e-12 |
| K4 | -0.49 | -1 | 0.0537668 | 4782.28 | 3.07e-16 | 2.07e-10 | 3.97e-13 | 7.18e-08 | 4.25e-12 |
| K4 | -0.2 | -2.45 | 0.01518 | 56179.1 | 1.48e-15 | 6.17e-10 | 2.48e-12 | 2.68e-06 | 1.38e-11 |

1333 CSV rows from 46 parameter cases.
