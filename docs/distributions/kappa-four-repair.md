# Kappa Four numerical repair and regression evidence

Review baseline: `94d1713e03a4c416dd472d6fbef0ce907881e9d1`. Approved repair implemented September 8, 2026.

The user approved the five numerical repair slices after reviewing the audit. This repair preserves public signatures, parameter order, serialization property names, constructor defaults, optimizer selection, fitting bounds, random-number behavior, and the July zero-kappa formula corrections. It does not modify Differential Evolution, restrict the distribution family, penalize the likelihood, or change convergence tolerances.

A declaration comparison against the baseline found no removed or changed public declaration. The four additions are overrides of the existing inherited `LogPDF`, `LogCDF`, `CCDF`, and `LogCCDF` methods.

## Implemented slices

| Slice | Result |
| --- | --- |
| Distribution evaluation and support | Logarithmic PDF/CDF/tail calculations; `Log1p`/`Expm1` divided differences; exact zero-shape limits; stable finite support; mathematical endpoint densities. A private compensated-arithmetic helper retains support residuals when ordinary double logarithms nearly cancel. Affine overflow is avoided when the final answer remains representable. |
| L-moments | Log-gamma probability-weighted moments, exact Gumbel and zero-hondo limits, integration of the actual quantile near zero shapes, strict finite-input and moment-existence checks. Hosking's iteration retains tolerance `1e-6`, at most 20 iterations, ten step reductions, initial shapes and estimation region. Exhaustion throws; final location/scale use moments at the accepted shapes. |
| Fitting | Accept initializers only with valid parameters, unchanged bounds and finite sample likelihood. Use the existing GEV initializer within established K4 bounds when needed. Require successful Nelder-Mead status and a valid finite-likelihood final fit. Estimate and bootstrap propagate failure without installing or returning failed estimates. |
| Moment and mode properties | Standardized quantile quadrature over the full probability interval, with central moments about the computed mean. Relative tolerance `1e-10`, absolute tolerance `1e-12`; return `NaN` for nonexistent or unresolved moments. Modes include endpoints and the existing mathematical stationary point; nonunique modes return `NaN`. Cache invalidation remains in the parameter setters. |
| Derivatives and documentation | Continuous zero-shape quantile gradients and stable nearby evaluations. Jacobian layout and LU determinant retained. Numerical domains, estimator failure and the unimplemented uncertainty methods documented. |

The moment integral pairs the two transformed tails before estimating error, so a near-zero mean is checked against the complete integral rather than two separately large contributions. Each half uses `p=u^8/2` or its survival counterpart, including its Jacobian. No tail truncation or actual-shape substitution is used.

`ParameterCovariance` and `QuantileVariance` remain explicitly unimplemented. Floating-point calculations and estimated quadrature errors are not a proof of absolute accuracy for all real inputs.

## Reported fitting regression

The exact 38 observations are frozen in `Test_KappaFourRegression.FittingSample`. The following controlled comparison was run against the audit baseline before source changes. These historical statuses describe that baseline; repaired probability arithmetic can change subsequent optimizer trajectories.

| Baseline configuration | Status | Iterations | Evaluations | Log likelihood |
| --- | --- | ---: | ---: | ---: |
| Current DE, seed 12345, population 40 | MaximumIterationsReached | 10,000 | 400,000 | -23.472600041817795 |
| Only restore pre-September-1 DE boundary repair | Success | 1,231 | 49,280 | -23.522258636620737 |
| Current DE with pre-July zero-kappa PDF behavior | Same as current DE | 10,000 | 400,000 | -23.472600041817795 |

No DE objective evaluation used exactly zero kappa in these comparisons. The July correction therefore did not explain the reported DE convergence change.

Current-DE best baseline parameters were `[-0.02304299227317367, 1.3774672309514362, 0.376420487117973, 1.264141217131441]`, with lower endpoint `0.28599211573600747`. Old-repair parameters were `[-0.022253095483755644, 1.3752818361277244, 0.34855569679255555, 1.2628508861048067]`, with lower endpoint `0.2859921157360075`. Both approach the sample minimum `0.2859921157360077`.

The old L-moment initializer excluded two observations: its lower support was approximately `0.3818457477`. Repaired initialization uses the existing GEV candidate `[0.750064308812039, 0.464304902648505, -0.176396212686591, 0]` within the original K4 bounds: lower `[-10, DoubleMachineEpsilon, -10, -2]`, upper `[10, 100, 10, 2]`. It has finite sample likelihood. The built-in MLE uses Nelder-Mead; unsuccessful termination is now reported as an exception and leaves the distribution unchanged.

### Why unrestricted MLE remains unresolved

Let `m` be the minimum observation. Set `alpha=1`, `kappa=0`, `hondo=1.5`, and `xi=m-log(1.5)-epsilon`, with positive epsilon tending to zero. These parameters remain within the existing bounds. The minimum observation approaches the lower endpoint from inside support, where density diverges because `hondo>1`; all other observation contributions remain finite. Thus the likelihood is unbounded and there is no finite global MLE for this problem.

The deterministic regression tests require increasing finite likelihood along this sequence. They do not require DE to find a finite global maximum. Selecting an MLE restriction or replacement estimator remains a separate scientific-policy decision. Infinite endpoint log-density contributions retain the existing negative-infinity likelihood convention.

## Relevant history

| Commit | Date | Audit finding |
| --- | --- | --- |
| `44c70af` | 2025-04-28 | PDF/CDF/quantile rewrite introduced the incorrect zero-kappa PDF and inverse expression. |
| `7329914` | 2026-05-24 | Inherited invalid-likelihood handling changed to negative infinity. |
| `b9af61c` | 2026-07-24 | Corrected the zero-kappa PDF and inverse expression; retained. |
| `48c99a2` | 2026-07-24 | Added finite-shape tests; did not change the class. |
| `2b57771` | 2026-09-01 | DE boundary repair changed; controlled comparison altered convergence status. |

Several L-moment defects predate these changes and exist in the earliest available September 2023 class history.

## Independent evidence and validation

The frozen oracle fixtures and their generator live in `Test_Numerics/Distributions/Univariate/Fixtures`. Python's standard-library Decimal evaluates the defining formulas at 90 digits using exact conversions of the input doubles. Shape derivatives are computed independently by high-precision symmetric differences. The 45 endpoint cases use 400 digits to retain subnormal distances above zero. Tests require no Python, R, network, or external numerical package at runtime.

Regeneration was checked byte for byte. SHA-256: probabilities `483358B5F45D909C6EDAB850FAF1C502B783691E408437557CA48BDBC2A3D454`; boundaries `B057D04D9D32F3FF45A73F6A2C7E4204EC3FC241DFAF01172EB3E494105E0909`.

- 1,546 frozen cases cover 15 values of each shape (both signs, exact zeros and near zeros), central probabilities and tails. Checks include quantile, PDF, CDF and both shape gradients.
- 45 additional frozen cases evaluate log probabilities within four representable steps of finite support boundaries. A one-step reverse-exponential boundary case is also checked.
- When a quantile rounds onto the public represented support endpoint, its PDF/CDF follow the documented one-sided limit; interior oracle comparisons are restricted to represented interior arguments.
- Analytical Gumbel, exponential, uniform, logistic, generalized Pareto and zero-hondo reductions, normalization, moment existence, endpoint modes, invalid data and the reported fitting sample supplement the formula grid.
- The heavy-tail mean counterexample `(0,1,-0.5,0)` is checked against `2*(sqrt(pi)-1) = 1.544907701811032`.
- Exact exponential/Gumbel inverse L-moments, iteration exhaustion and gamma-overflow examples have dedicated regressions.

Formula/source cross-checks from the audit: [SciPy kappa4](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kappa4.html), [Hosking's lmom Fortran](https://raw.githubusercontent.com/cran/lmom/master/src/lmoments.f), and [lmomco parkap](https://rdrr.io/cran/lmomco/src/R/parkap.R). The referenced [nsRFA implementation](https://rdrr.io/cran/nsRFA/src/R/KAPPA.R) was not accepted as an unquestioned oracle; its zero-shape substitution is numerically unreliable.

### Completed validation

`dotnet build -c Release` completed with **zero warnings and zero errors** on all four target frameworks. XML documentation enforcement was enabled.

The final full run (`dotnet test -c Release --no-build --logger 'trx;LogFilePrefix=kappa-four-repair'`) produced:

| Framework | Passed | Failed | K4 methods passed |
| --- | ---: | ---: | ---: |
| net481 | 2,524 | 1 | 32 |
| net8.0 | 2,524 | 1 | 32 |
| net9.0 | 2,524 | 1 | 32 |
| net10.0 | 2,524 | 1 | 32 |

The sole failure on each framework was `Data.TimeSeriesAnalysis.Test_TimeSeriesDownload.BOM_FullPor_Goodradigbee_Discharge`: the live BOM API returned HTTP 500, `DatasourceError`, `Error connecting to WDP.` The same failure persisted in a focused rerun. An earlier full run also saw intermittent Cotter River and Murray River API failures; these passed in the final run. No download code or tests were changed or suppressed.

All 32 K4 methods (12 original and 20 new) passed on each framework, including rejection of the reported sample's iteration-exhausted MLE and propagation of a seeded bootstrap estimation failure. The full-run TRXs were inspected to verify those exact K4 outcomes.

Final evidence files in `Test_Numerics/TestResults`:

- `kappa-four-repair_net481_20260908112347.trx`
- `kappa-four-repair_net8.0_20260908112339.trx`
- `kappa-four-repair_net9.0_20260908112337.trx`
- `kappa-four-repair_net10.0_20260908112340.trx`

Haden explicitly authorized committing these reviewed changes with the documented external BOM service failure on September 8, 2026. This exception applies to this repair's commit gate. The MLE estimator-policy decision remains separate, as approved in the repair plan.
