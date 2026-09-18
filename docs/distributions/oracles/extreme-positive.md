# Extreme and positive distribution oracle fixtures

Generated with **R version 4.4.3 (2025-02-28 ucrt)**, using only base R and `stats`.
The generator does not load Numerics, invoke an optimizer, run a test suite, or
change production files. Its probability, moment, and derivative calculations
are independent of the implementation being repaired.

Run from the repository root:

```powershell
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla docs/distributions/oracles/extreme-positive.R
```

An optional first argument changes the CSV output path. The generator requires
R 4.4.3 so a reference change cannot silently substitute another R version.
Regeneration is deterministic, with LF line endings and 17-digit numeric text.

## Contents and parameter mapping

| CSV family | Numerics parameters | Independent mapping | Rows |
| --- | --- | --- | ---: |
| Exponential | xi, scale=Alpha | R `dexp/pexp/qexp` in unit-scale coordinates, then affine scale | 97 |
| Gamma | scale=Theta, shape=Kappa; xi=0 | R `dgamma/pgamma/qgamma(shape=Kappa, scale=Theta)` | 701 |
| GEV | xi, scale=Alpha, shape=Kappa | Hosking GEV defining formulas; SciPy `genextreme(c=Kappa)` | 560 |
| GPA | xi, scale=Alpha, shape=Kappa | Hosking GPA defining formulas; SciPy `genpareto(c=-Kappa)` | 513 |
| Gumbel | xi, scale=Alpha | Maximum Gumbel defining formulas; SciPy `gumbel_r` | 47 |
| Weibull | scale=Lambda, shape=Kappa; xi=0 | R `dweibull/pweibull/qweibull(shape=Kappa, scale=Lambda)` | 579 |

The CSV has **2,497 rows**. All kurtoses are ordinary kurtosis, not excess.
Xi is zero except for selected tail-overflow cases; affine location/scale
metamorphic checks supplement the external goldens. Shapes include both signs,
exact zero, and nonzero values as small as +/-1e-12. Scales include 1e-200 and
1e200, plus selected 1e-320 log-density reproductions.

The added GEV shape 50, 100, and 200 moment rows use R `lgamma` and normalized
central-moment ratios. Scale is combined with the moment logarithm before
exponentiation, so finite scaled moments remain referenceable when a raw Gamma
function or the corresponding unit-scale moment overflows. These 36 rows have
relative tolerance 2e-10 and exact finite/infinite classification.

Eighty-four additional Weibull moment rows use the defining transform
`X = scale * T^(1/kappa)`, where T has the unit exponential distribution.
Shapes .005, .01, and .02 use normalized R `lgamma` central moments at scales
1 and 1e±200, with relative tolerance 2e-10. Shapes 1e3, 1e6, 1e12, and 1e200
use full-probability integration of the centered, divided transform
`expm1(c*log(T))/c`, with `c=1/kappa`, at the same scales. Integrating this
transform preserves its central moments as c tends to zero; scale and c are
restored through logarithms. These large-shape rows use relative tolerance
2e-8. No production GEV moment routine participates in either reference.

The final review added 56 rows without changing any of the earlier 2,441 rows:

- Thirty-six GEV/GPA density and tail logarithms evaluate the defining Hosking
  support expression from physical-coordinate logarithms when the difference,
  standardized observation, or shape product exceeds the floating-point range.
- Twelve Gamma covariance entries use R integration of the scaled Fisher
  residual kernel for shapes 1e16 and 1e155, with the same Fisher inverse and
  physical scale/count factors. In the coordinate `u=kappa*t`, the integrated
  kernel is `exp(-u)*(u/2+u^2/(12*kappa))`; the first omitted integrated
  Bernoulli term has magnitude `1/(30*kappa^3)`.
- Two Gamma median variance rows use the defining concentrated-shape limit
  `kappa/n` for MLE and MoM at kappa=1e16. The median expansion
  `q=kappa-1/3+O(1/kappa)` gives a relative variance correction of order
  `1/kappa`, below the rows' 2e-10 relative acceptance. These rows protect
  the positive mean-direction term from cancellation in a rounded matrix.
- One Weibull variance reference analytically contracts its retained rounded
  covariance coefficients at lambda=kappa=1e100. One median reference combines
  the defining quantile logarithm with lambda=1e200 before exponentiation.
- Four GPA skewness/kurtosis entries evaluate the same rational moments in
  reciprocal-shape coordinates at shapes 1e200 and 1e308.

The GEV/GPA definitions are evaluated with `expm1` divided differences, without
zeroing small nonzero shapes. Mathematical support and one-sided endpoint
densities are explicit. Infinite density at an integrable endpoint is retained.
These endpoints must not be replaced blindly with a package's endpoint policy.

R 4.4.3 `dweibull(...,log=TRUE)` and `pweibull(...,log.p=TRUE)` can return
negative infinity after an internal positive power underflows, even when the
mathematical logarithm is finite. The selected `k=100,z=1e-20` rows therefore
use the defining log-density and log-CDF limit, identified by the oracle label
`Weibull-defining-log-formula-R-underflow`. The density uses
`log(k)+(k-1)*log(z)-exp(k*log(z))-log(scale)`; the CDF logarithm equals
`k*log(z)` to binary64 precision after this underflow. Other native R rows
remain unchanged.

## Schema and tolerance interpretation

Coordinates: `family,case,xi,scale,shape,x,p,n`. `NaN` in an unused coordinate
means not applicable, not an invalid-parameter test. `quantity` names the result;
`value` is its independent expected value. Readers must map `Inf`, `-Inf`, and
`NaN` to their corresponding floating-point classifications.

For finite rows, use `absolute_tolerance + relative_tolerance * abs(value)`;
classification must be checked separately for infinite/undefined rows. A
relative-only tolerance intentionally protects tiny nonzero values from passing
as zero. `estimated_absolute_error` records an oracle convergence indicator,
not a rigorous proof or an additional allowed test tolerance.

`Gradient1/2/3` follow the distribution's public parameter order. Gamma instead
uses `GradientScale` and `GradientShape` for clarity. `Covarianceij` is the
upper triangle of the public-order covariance matrix, already divided by `n`.
`CovarianceDefined=0` rows denote unsupported uncertainty domains; their
`case` distinguishes MLE from MoM. They are not invalid distribution parameters.
The standalone `QuantileVariance` and `Median` rows are exercised by focused
regressions in addition to the family-wide probability/moment/covariance readers.

Most native R/analytic rows use relative tolerance 2e-11. Centered moment
integration uses 2e-8 relative plus 2e-9 absolute, protecting zero centered
moments from meaningless relative checks. Covariance uses 2e-9 relative plus
1e-12 absolute. Weibull gamma-function moment cancellation receives 2e-8
relative tolerance. Exact classification still applies at moment-existence
boundaries.

## Gamma actual-quantile derivatives

The **49 finite shape derivatives** differentiate actual R `qgamma` results;
they do not differentiate a Wilson-Hilferty or Cornish-Fisher approximation.
The shape grid is `.001,.01,.1,.5,1,2,10,100,10000`, and probabilities are
`1e-12,1e-6,.01,.5,.99,1-1e-12`.

The generator differentiates `log(qgamma(p, shape=k*exp(u)))` with respect to
`u`, then multiplies by `q/k`. A five-point central derivative is extrapolated
to sixth order, and refinement selects the most stable pair of estimates.
Using log-quantiles retains the shape=.001 median near 5.24e-302. Upper-tail
quantiles are requested with `lower.tail=FALSE, log.p=TRUE` to avoid derivative
noise from an unnecessary lower-tail inversion.

Every finite derivative is independently cross-checked by differentiating
`pgamma` at the fixed quantile and applying implicit differentiation. Lower or
upper log-probability is selected according to the requested tail. The maximum
observed relative disagreement in the frozen grid is below **5e-10**; fixture
acceptance is a conservative **5e-8 relative**. The generator refuses to write
the CSV if this cross-check exceeds that acceptance threshold.

These are convergence-checked **binary64** references, not arbitrary-precision
calculations. Five lower-tail quantiles underflow in R; their shape derivatives
have status `unresolved-quantile-underflow` and value NaN. They must not become
zero-derivative expected values. Quantiles/scale derivatives that underflow are
separately marked `underflowed`. A broader extreme-tail contract would require
an independently validated arbitrary-precision oracle.

There are no legacy named quantile-approximation API rows in this file. Those
APIs need their own separately named fixtures if their existing contracts are
changed; they must not replace these actual-quantile references. The two Gamma
large-shape variance-limit rows explicitly identify their asymptotic bound.

## Moment and covariance provenance

GEV/GPA centered moments are integrated over the **full probability interval**
with paired tails `p=u^8/2`, relative tolerance 1e-10, absolute tolerance 1e-12,
and a 1,000-subdivision ceiling. Existence is checked first: moment order r
requires `kappa > -1/r`. No probability-tail truncation is used. Undefined
moments follow the planned Numerics NaN convention. Other moments use analytic
identities and R's gamma function with scale applied after standardization.

GEV expected Fisher information is independently integrated from analytic
scores in the latent Gumbel coordinate; its positive-definite inverse is frozen
for kappa `-.2,-1e-6,0,1e-6,.2,.49`. No Numerics information implementation or
matrix routine participates. GPA MLE/MoM covariance uses Hosking-Wallis formulas,
with separate rejected-domain rows. Gamma covariance uses R `trigamma` for MLE
and an independent sample-mean/sample-variance Jacobian for MoM. Gumbel/Weibull
covariance uses Euler/pi constants, rather than the production rounded decimals.

Exponential MLE covariance deliberately corresponds to the **actual MLE**
`location=min(sample), scale=mean(sample)-min(sample)`: diagonal entries are
`scale^2/n^2` and `scale^2*(n-1)/n^2`, with zero covariance. This differs from
bias-corrected estimators and is explicitly named `actual-MLE-order-statistics`.
The n=50000 case also protects against integer multiplication overflow.

## Sources

- [R pgamma implementation](https://raw.githubusercontent.com/wch/r-source/trunk/src/nmath/pgamma.c)
- [R qgamma implementation](https://raw.githubusercontent.com/wch/r-source/trunk/src/nmath/qgamma.c)
- [SciPy 1.16.2 distribution source, for parameterization cross-checks](https://raw.githubusercontent.com/scipy/scipy/v1.16.2/scipy/stats/_continuous_distns.py)
- [Hosking and Wallis, Parameter and Quantile Estimation for the Generalized Pareto Distribution, 1987](https://stat.cmu.edu/technometrics/80-89/VOL-29-03/v2903339.pdf)

The current CSV SHA-256 is
`E5824F8D5E09B2D582F42DB962F5F49A24D87641F8B5E3DFC1913BCE0BC72BCA`.
An independent regeneration on 2026-09-08 produced the same hash and all
2,497 coordinate/quantity combinations are unique. The task's
`extreme-positive-report.md` records current implementation and validation
evidence; `extreme-positive-oracles-report.md` retains the earlier 2,321-row
reference-generation checkpoint.
