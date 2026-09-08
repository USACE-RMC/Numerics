# Distribution robustness and uncertainty evidence

This record accompanies the [approved repair plan](bestfit-robustness-plan.md). The implementation starts from Numerics `c0d67b9` and BestFit `b732703`, in isolated `codex/distribution-robustness` worktrees. The original checkouts, including unrelated BestFit changes, are preserved. No package release, push, or BestFit Verification execution is part of this work.

## Probability and likelihood contracts

The 15 BestFit families now provide direct logarithmic lower and upper tails. Normal tails use a continued-fraction Mills ratio in the remote tail. Gamma tails use direct lower series, upper continued fractions, a small-shape complementary series and a large-shape transition expansion; actual gamma quantiles use those tails and implicit shape differentiation. Named approximate frequency-factor methods remain separate APIs. GEV, GPA, GLO and GNO evaluate the actual nonzero shape with continuous divided differences, retaining its finite support endpoint.

`LogPDF` preserves mathematical endpoint limits, including positive infinity. The established aggregate sample log-likelihood convention still maps nonfinite totals to negative infinity. Zero censoring counts contribute zero. Interval likelihood means `(lower, upper]`, uses the better-conditioned difference of logarithmic tails, and handles mixture atoms explicitly. A collapsed continuous interval has a density-quadrature fallback, limited to explicitly recognized continuous families.

Independent CompetingRisks combines logarithmic survival/CDF/density terms and resolves endpoint product limits. Its minimum support is the minimum of component lower endpoints through the minimum of upper endpoints; maximum support uses the corresponding maxima. Dependent probability-combination rules are retained, with bounded numerical density differentiation and explicit failure for unresolved derivatives. Ordinary mixture moments combine component central moments about a scaled mixture center. Hurdle conditioning uses logarithmic positive mass, with the zero atom included explicitly in moments and intervals. Mutable component, weight, base, transform and dependence settings invalidate derived caches.

Representative regressions include:

| Case | Independent expected result |
|---|---:|
| Standard Normal `LogCDF(-40)` | -804.6084420137538 |
| Standard Normal log probability of `(9,10]` | -43.628216632280818 |
| Standard Normal log probability of `(0,1e-20]` | -46.970640393085588 |
| Weibull(scale 2, shape 2) median gradient | [0.832554611157698, 0.152571011039570] |
| GNO(0,1,2) mean | -3.1945280494653252 |
| Mixture of identical Normal(1e10,1) components, SD | 1 |
| Weights .9/.1 at Normal means -1e308/+1e308, SD | 6e307 |
| Minimum of two unit exponentials, `LogPDF(1000)` | -1999.3068528194401 |
| Maximum of two Gamma(scale 1, shape .5), density at zero | 4/pi |

BestFit's aggregate and pointwise univariate/point-process quantile priors consume the additive `IStandardError.LogAbsQuantileJacobian` extension. The existing interface and `QuantileJacobian` signatures are unchanged. Scaled elimination sums log pivots; an exact binary-rational fallback distinguishes an unresolved small pivot from true singularity without installing an artificial pivot. Repeated probabilities produce negative-infinite log absolute determinant. Finite log determinants remain available when the raw determinant overflows or underflows.

BestFit mixture EM keeps exact, censored, interval, positive-conditional and measurement-error observation calculations logarithmic through responsibility normalization. Component log observations are centered before adding log weights, so identical components retain weights .4/.6 even at a common log density near -5e199. Existing GL20 measurement-error nodes, domain and retained-mass normalization are preserved. Failed automatic initialization is reported through model validation while retaining editable parameter state; it does not fabricate observations or prior bounds.

## Public parameter coordinates and uncertainty

Every covariance row/column and quantile-gradient entry follows the public parameter order below. `n` is the number of independent observations; it is not a posterior draw count. Probabilities supplied to uncertainty methods must be finite and strictly between zero and one, and `n` must be positive.

| Family | Public coordinates | Implemented uncertainty |
|---|---|---|
| Exponential | Xi, Alpha (location, scale) | Existing MoM; exact covariance of implemented minimum/mean-minus-minimum MLE |
| GammaDistribution | Theta, Kappa (scale, shape) | Existing MoM and MLE, with actual-quantile gradients |
| GEV | Xi, Alpha, Kappa | Full three-parameter local MLE, Kappa < 1/2 |
| GLO | Xi, Alpha, Kappa | New local MLE, abs(Kappa) < 1/2 |
| GNO | Xi, Alpha, Kappa | New closed-form local MLE for finite shape and positive finite scale |
| GPA | Xi, Alpha, Kappa | Existing MLE for Kappa < 1/2; MoM for Kappa > -1/4; location covariance also requires n + 2*Kappa > 0 |
| Gumbel | Xi, Alpha | Existing MLE covariance and delta method |
| KappaFour | Xi, Alpha, Kappa, Hondo | New local MLE for Kappa < 1/2, Hondo < 1/2 and Kappa*Hondo < 1/2 |
| LnNormal | Mean, StandardDeviation in observation space | Existing indirect log-moment/MLE estimators, transformed into physical coordinates |
| Logistic | Xi, Alpha | Existing MoM/MLE covariance and delta method |
| LogNormal | Mu, Sigma of log-base observations | Existing indirect log-moment/MLE estimators; Base is fixed, not an estimated coordinate |
| LogPearsonTypeIII | Mu, Sigma, Gamma of log-base observations | Existing MoM and regular full-family MLE; fixed Base |
| Normal | Mu, Sigma | Existing MoM/MLE leading asymptotic covariance |
| PearsonTypeIII | Mu, Sigma, Gamma of observations | Existing MoM and regular full-family MLE |
| Weibull | Lambda, Kappa (scale, shape) | Existing MLE covariance with repaired actual-quantile gradient |

Pearson and Log-Pearson regular MLE covariance requires `abs(Gamma) < sqrt(2)`, equivalent to shifted-gamma shape greater than two. The smooth zero-skew covariance is `diag(Sigma^2, Sigma^2/2, 6)/n`, retaining skew-estimation uncertainty. GPA's scalar quantile variance continues to omit its order-`n^-2` location term, as documented by the existing method. These uncertainty domains do not restrict distribution validity or fitting bounds.

For Exponential MLE the covariance is `Alpha^2 * diag(1/n^2, (n-1)/n^2)`. For LnNormal, source inspection resolved an ambiguity in the initial review: `Estimate(MethodOfMoments)` fits moments of the log observations. A direct physical-moment estimator would have a different covariance and is not substituted under that enum. At physical mean/SD `(1,1)`, median variance for the existing indirect estimator and `n=100` is .0034657359027997265; the initially proposed .00875 belongs to a different estimator. Supplied physical LnNormal parameters are preserved exactly through serialization and cloning until log-coordinate setters change them.

GNO uses the approved inverse expected information with `v=k^2`, `c=-expm1(-v/2)/v`, `R=((1+v)*exp(v)-1-2*v)/v^2` and `S=diag(Alpha,Alpha,1)`:

```
Cov = S * C * S / n
C = [[1+c*c/R, -k, c/R],
     [-k, v+1/2, k/2],
     [c/R, k/2, v/2+1/R]]
```

At zero shape, `c=1/2`, `R=3/2`, and `C=[[7/6,0,1/3],[0,1/2,0],[1/3,0,2/3]]`. Series and logarithmic scale restoration avoid cancellation and intermediate overflow. Its scalar delta variance is evaluated through an exact sum of nonnegative squares, preserving finite uncertainty after endpoint-gradient cancellation. Other scalar delta methods similarly contract in normalized coordinates before restoring physical scale. This does not regularize a covariance or replace the estimator.

GLO, Kappa and GEV integrate analytical fixed-observation score products over both complete probability tails. GLO fixes Hondo=-1 and inverts the three-parameter principal **information** block; GEV fixes Hondo=0. The implementation checks relative error `1e-10`, absolute error `1e-12`, score means, positive-definite information, finite solves and the inverse residual. It adds no jitter, clipping or pseudoinverse. Integration/solve failures throw `InvalidOperationException`; mathematical regularity violations throw `ArgumentOutOfRangeException`. These matrices describe local asymptotic uncertainty, not global MLE existence, optimization convergence or finite-sample coverage.

## Moments, support and initialization

GNO uses analytical transformed-normal moments. GLO stabilizes its existing reciprocal-sinc formulas. GEV and GPA require `k > -1/r` for moment order `r`; GLO requires `r*abs(k) < 1`. Log-Pearson guards each required gamma moment-generating-function value before evaluating its central moments. Lognormal standardized moments depend only on shape and fixed base. Weibull uses log-Gamma moments and the exact centered-power relation to GEV to preserve very small and very large shape limits. Ordinary mixture moments avoid subtraction of large raw moments, including when the component locations themselves span nearly the full floating-point range.

The inherited `CentralMoments(int)` remains an approximate bin calculation, now accumulating central powers in scaled coordinates. `CentralMoments(double)` retains its documented numerical approximation. These overloads are distinct from analytical scalar properties. Competing and positive-truncated mixture moments use checked full-support central integration where general component formulas are unavailable; they do not silently renormalize estimated mass or clip divergent moments.

Initialization validates the supplied sample and candidate vector, including finite values, required positivity, minimum length and nonconstant data. Signed location bounds handle negative and zero-centered samples; positive scale bounds are finite and ordered with feasible starts. Existing fitting algorithms, estimator enums, random seeds, convergence settings, public parameter coordinates, log bases, serialization conventions and existing minimum-scale policies are retained.

Corrected legacy expectations include GEV/GPA moment existence, median/mode values, the GLO shape-one endpoint mode, and the Weibull standard-error test that instantiated Gumbel. Correcting a Kappa endpoint allowed an old optimizer fixture to converge under unchanged settings; its requirement that that particular sample exhaust iterations was replaced by a deterministic invalid-data rejection/unchanged-state regression. Endpoint evaluation, not optimizer tuning, changed that trajectory.

## Independent oracles and reproducibility

The new oracle regression tests embed CSV values and do not load R, Python, a package service or a network source. Companion generators and provenance are in [oracles](oracles). R 4.4.3 `stats` supplies ordinary-family tails/quantiles; standalone defining formulas cover transformations and R underflow limitations. Python 3.12 Decimal/Fraction scripts use exact input doubles and up to 400-digit arithmetic for adjacent Kappa endpoints, scalar variance identities and gamma expansion coefficients.

| Evidence | Coverage |
|---|---|
| [extreme-positive.csv](oracles/extreme-positive.csv), [.R](oracles/extreme-positive.R), [notes](oracles/extreme-positive.md) | 2,497 values for six families; 49 checked actual gamma derivatives, log tails, covariance, moments, scale extremes and existence boundaries |
| [generalized-fisher.csv](oracles/generalized-fisher.csv), [.R](oracles/generalized-fisher.R), [notes](oracles/generalized-fisher.md) | 1,333 rows, including 46 information/covariance cases and 88 moment values; score means, alternative tail transformations and independent fixed-observation score checks |
| [generalized-adjacent-boundaries.csv](oracles/generalized-adjacent-boundaries.csv), [.py](oracles/generalized-adjacent-boundaries.py) | First through fourth interior doubles at 45 finite Kappa endpoints, 180 rows |
| [generalized-scalar-variance.py](oracles/generalized-scalar-variance.py) | High-precision scalar contractions and the exact GNO positive quadratic at extreme scales/shapes |
| [generalized-affine-overflow.R](oracles/generalized-affine-overflow.R), [notes](oracles/generalized-affine-overflow.md) | Both GNO/GLO shape signs with finite log tails after the standardized observation overflows |
| [normal-pearson.csv](oracles/normal-pearson.csv), [.R](oracles/normal-pearson.R), [notes](oracles/normal-pearson-evidence.txt) | Reflected/transformed tails, covariance coordinates and standardized-moment/uncertainty scale regressions |
| [generate-gamma-temme-coefficients.py](oracles/generate-gamma-temme-coefficients.py) | Exact rational derivation of the gamma transition-expansion coefficients |

Each oracle records its mapping, precision/tolerance and limitations. The [manifest](oracles/manifest.json) records row counts and checksums with portable LF normalization. GNO is Hosking's transformed normal, not SciPy `gennorm`; GLO is Hosking's generalized logistic, not SciPy `genlogistic`. The R Weibull routines can internally underflow a positive power before computing a logarithm; selected such rows explicitly use the independent defining log formula. Five gamma derivative cases whose reference unit quantile underflows are identified as unresolved, not converted into zero expected derivatives. Published rounded Gumbel/Weibull covariance coefficients retain their documented precision.

Primary references: [R normal log-tail interface](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Normal.html), [R gamma interface](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/GammaDist.html), [Hosking-compatible lmom definitions](https://raw.githubusercontent.com/cran/lmom/master/R/lmom.r), [DLMF gamma asymptotic expansions](https://dlmf.nist.gov/8.12), and [Hosking and Wallis GPA estimation](https://stat.cmu.edu/technometrics/80-89/VOL-29-03/v2903339.pdf).

## Retained limitations

- New GNO/GLO/Kappa covariance is MLE only. L-moment and product-moment covariance still require estimator-specific derivations. Other previously unsupported combinations remain unsupported; composites do not gain generic component-estimator covariance APIs.
- A mathematically regular covariance can remain numerically unresolved near a regularity boundary or outside double range. Checked Fisher integration reports failure instead of modifying the information. Scalar variance can sometimes remain representable even when a full physical covariance is not.
- Dependent CompetingRisks retains its existing ordinary-probability dependence backend; logarithmic combinations are direct for the independent case. This change does not introduce a different dependence model or an arbitrary-precision multivariate tail calculation.
- Subtracting two component log values cannot recover corrections already lost by rounding those values. For example, a positive-conditioned Normal with mean -1e100 and scale 1 requires a specialized conditional-ratio calculation to resolve its density around x=1e-100; generic double-valued log-tail subtraction is insufficient there. This extreme conditional-ratio domain remains a documented limitation.
- The tiny continuous-interval density fallback is a fixed quadrature rule, not a rigorous arbitrary-precision error enclosure. Generic moment quadrature can explicitly fail for unresolved or divergent integrals. Existing full uncertain-observation likelihood integration outside the repaired EM path retains its established arithmetic and quadrature policy.
- The separate, pre-existing `GammaDistribution.MLE_NR` method uses the instance `Kappa` rather than the newly solved shape when computing its scale. BestFit uses the main MLE; that alternate-estimator defect was recorded during review and left outside this repair inventory.

## Validation and review

Targeted red runs preceded repairs, including the final scale/cancellation cases. Independent code review covered family formulas, shared probability/derivative helpers, composite semantics and BestFit integration. The [machine-readable results](bestfit-robustness-validation.json) record final gate counts, report checksums and the matching Numerics dependency hash in all three BestFit test outputs.

| Gate | Result |
|---|---|
| Numerics Release build, XML enforcement, net481/net8.0/net9.0/net10.0 | 0 warnings, 0 errors |
| Numerics complete Release tests, all four frameworks | 2,675/2,675 passed on each; 10,700 passes total, no failures or skips |
| BestFit core fast project, explicit Numerics worktree | 3,406/3,406 passed; strict documentation; one test worker |
| BestFit UI fast project, explicit Numerics worktree | 593/593 passed; strict documentation |
| BestFit App fast project, explicit Numerics worktree | 443/443 passed; strict documentation |
| BestFit Verification | Not executed; outside authorization |

Numerics commands are `dotnet build -c Release -p:EnforceXmlDocumentation=true` and `dotnet test -c Release --no-build`. BestFit project commands set `UseLocalRmcNumerics=true` and `RmcNumericsProjectPath` to the intended isolated Numerics project's absolute path. The BestFit MSTest.Sdk 3.6.4 host uses Microsoft.Testing.Platform; targeted filters are application arguments after `--`, and every reported result count was checked in the test output/TRX.

The initial complete gates exposed an unintended change from Normal's established decade bounds and shape-one subnormal quantile round-trips on .NET Framework. Both were corrected; the seeded RWMH test now reproduces every existing literal unchanged, and all four final framework runs pass. One initial live BOM download returned HTTP 500; all four final runs pass that unchanged network test.

BestFit's unrelated `RunAsync_MultipleAnalyses_Parallel` stopwatch test passed in isolation but failed under default inter-test contention (369 ms and 673 ms against its unchanged 250 ms limit). The final core gate uses the optional `docs/validation/distribution-robustness.runsettings` in the paired BestFit worktree to schedule test cases with one worker. All 3,406 cases execute; the test still runs its three mock analyses concurrently and uses its original assertion. UI/App use their normal scheduling. No runner default, test threshold, application concurrency or production numerical setting was changed. The configuration mechanism is documented by [MSTest execution control](https://learn.microsoft.com/en-us/dotnet/core/testing/unit-testing-mstest-writing-tests-controlling-execution).
