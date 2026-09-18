# BestFit distribution robustness and uncertainty repair

Approved implementation plan, September 8, 2026.

## Baselines and authorization

- Numerics: `c0d67b9c52d62dc49bdbfb744ec47b9a205cc86b`.
- BestFit: `b7327036ecd25770201597dc163945353370347e`.
- Haden explicitly approved implementation of the review and repair plan, including the stated formula, reference, support, derivative, moment, covariance and focused BestFit changes. Missing MLE uncertainty for GNO/GLO/KappaFour is included. Other estimator combinations remain unsupported.
- Work in isolated worktrees; preserve unrelated original-checkout changes. Do not push/publish. Preserve optimizer selection, convergence settings, seeds, dependence rules, public signatures and serialization coordinates. Additive log-Jacobian extension is approved.
- No BestFit Verification runs without separately named authorization. All numerical oracle regression tests belong in Numerics; BestFit fast tests cover integration contracts.

## Review evidence and required behavior

All 15 BestFit distributions plus CompetingRisks/Mixture were reviewed. Only KappaFour had direct lower/upper log tails. Baseline net8 build had zero warnings/errors and all 344 selected distribution tests passed despite reproduced defects.

| Family | Approved repairs |
|---|---|
| Exponential | Stable tiny probabilities and log tails; covariance integer overflow and implemented-MLE estimator alignment; validation/constraints. |
| GammaDistribution | Shape-one endpoint; direct log gamma tails; exact-quantile derivatives; scale-stable moments; frequency-factor reflection/cap corrections. |
| GEV | Exact-zero shape/support, gradients/Fisher limits, correct moment existence, median/mode/endpoints, feasible constraints. |
| GLO | Exact-zero shape/support, gradients, stable existing analytical moments, mode/endpoints; MLE covariance. |
| GNO | Exact-zero shape/support, analytical moments and gradients, symmetric L-moment initializer, mode/endpoints; closed-form MLE covariance. |
| GPA | Exact-zero shape/support, gradients, moment existence, median/mode/endpoints, covariance domains and validation. |
| Gumbel | Direct log tails; affine and scale overflow avoidance; input/constraint validation. |
| KappaFour | Retain recent repairs; MLE uncertainty; singular Jacobians and first-through-fourth adjacent boundary regressions. |
| LnNormal | Public real-space mean/SD covariance and gradients; correct MoM variance; conversion and negative-mean validation. |
| Logistic | Stable log density/tails; signed/zero-center initialization bounds; uncertainty validation. |
| LogNormal | Base/moment conversion, mode, base validation, shape-only moments, actual-quantile gradient with complete covariance transformations. |
| LogPearsonTypeIII | Moment existence, reflected tails, endpoints, transformed mode, base preservation, correct uncertainty coordinates and exact quantile derivatives. |
| Normal | Stable normal log tails and log scale evaluation; uncertainty/constraints validation. |
| PearsonTypeIII | Reflected tails, endpoints, mode/constraints, full public-coordinate covariance and regularity, exact derivatives and full zero-skew limit. |
| Weibull | Direct log tails, actual-quantile derivatives, full-rank Jacobian, stable analytical moments and scale-homogeneous uncertainty; correct the test that constructed Gumbel. |
| CompetingRisks | No invented out-of-support density; independent log combinations, min/max support, supplied-vector validation, cache invalidation; audit dependent differentiation without changing copula rules. |
| Mixture | Inactive/singular components; conditional log tails/support/caches; central component-moment combination and checked positive-conditional integration. |

Representative independent regressions: Normal log-CDF(-40)=-804.6084420137538; Normal interval (9,10] log probability=-43.628216632280818; Weibull(2,2) median gradient=[.832554611157698,.152571011039570]; GNO(0,1,2) mean=-3.1945280494653252; identical Normal(1e10,1) mixture SD=1; minimum of two unit exponentials log-PDF(1000)=-1999.3068528194401. Existing GEV/GPA moment/median/mode goldens and a Weibull test incorrectly targeting Gumbel require corresponding fixes.

## Implementation and interfaces

1. Freeze oracle fixtures before each repair and run their failing regressions. Record versions/mappings/precision/tolerances. Tests need no external numerical runtime/network.
2. Internal probability helpers: stable log complements/differences, normal log tails, gamma log P/Q and implicit gamma-quantile shape differentiation. Preserve represented support endpoints and true one-sided density limits. Exact nonzero shapes must not be flattened to zero.
3. Shared likelihood: stable differences of log tails for (lower,upper] intervals, including atoms; zero censoring count contributes zero. Keep aggregate infinite-log-density handling.
4. Existing uncertainty: all gradients/Jacobians/covariances use public parameter order/coordinates; exact quantile derivatives, strict finite interior probabilities, appropriate positive sample sizes, unchanged unsupported estimator behavior. GEV/GPA MLE k<.5; GPA MoM k>-.25 and its location-moment sample-size condition; Pearson regular MLE |gamma|<sqrt(2), including full smooth zero-skew covariance diag(sigma^2,sigma^2/2,6)/n. Restrictions apply only to uncertainty, not fitting families.
5. New MLE uncertainty: GNO closed form, GLO analytical Fisher for |k|<.5, Kappa analytical Fisher for k<.5,h<.5,kh<.5. Integrate full probability domain in standardized coordinates, paired tails, relative tolerance 1e-10 and absolute 1e-12. Scale weighted scores before outer products. Verify integration status/errors, positive definiteness and solve residuals. No jitter, clipping or pseudoinverse. Report unresolved calculations explicitly. Compute quantile variance g^T Cov g. These are local asymptotic uncertainties, not guarantees of a global MLE.
6. Analytical moments: stable GNO formulas, factored GLO series; GEV/GPA rth moment exists iff k>-1/r. Guard LP3 moments and stabilize lognormal standardized moments. Ordinary mixtures combine central component moments; positive-conditional cases retain checked integration. Document generic approximate numerical CentralMoments separately.
7. Support/validation/initializers: validate candidate vectors, correct bounds and inactive components, refresh derived caches when dependencies change; reject invalid/degenerate data; preserve signed log coordinates, bases and serialization.
8. Add LogAbsQuantileJacobian extension without modifying IStandardError. Use scaled pivoting and log pivots; exact singularity=-Infinity without artificial pivots. Update BestFit aggregate/pointwise Univariate and PointProcess priors, and preserve logarithmic mixture EM observation/responsibility calculations.

## Oracle definitions

- R 4.4.3 stats: Normal(mu,sigma); logistic(xi,alpha); exponential(xi,scale alpha); gamma(shape k,scale theta); Weibull(shape k,scale lambda).
- LnNormal oracle uses meanlog/sdlog internal conversion from public physical moments. LogNormal uses meanlog=Mu*log(Base), sdlog=Sigma*log(Base).
- PIII a=4/gamma^2, signed beta=sigma*gamma/2, xi=mu-2*sigma/gamma; shifted/reflected gamma with direct complementary tails. LP3 is Base^PIII with density Jacobian 1/(x log(Base)).
- Hosking GEV uses repository shape k; scipy genpareto uses shape -k. GLO/GNO use Hosking definitions, not scipy genlogistic/gennorm. scipy kappa4 shape order is (h,k).
- Independent sources: R stats numerical functions; Hosking lmom source; 90+ digit defining formulas and independent fixed-observation log-density differences for Fisher matrices. Do not call production helpers to generate goldens.

## Execution ledger

- [x] Original baselines and dirty state recorded; isolated worktrees created.
- [x] Shared probability/uncertainty helpers and likelihood regressions.
- [x] Exponential/Gamma/GEV/GPA/Gumbel/Weibull repairs.
- [x] Normal/logistic/lognormal/Pearson repairs.
- [x] GLO/GNO/KappaFour repairs and new uncertainty.
- [x] CompetingRisks/Mixture repairs.
- [x] BestFit integration repairs.
- [x] Independent code review and corrective follow-up.
- [x] Full required Release/XML/four-framework Numerics checks and all BestFit fast projects.
- [x] Reviewed implementation committed locally; current destination branches incorporated and combined Release gates passed. Local integration targets are `bug-fixes-and-enhancements` and `documentation-verification-updates`; no push.

## Completion gates

Cover signs, zeros/near-zeros, adjacent boundaries, scale/translation extremes, invalid/sample/constraint cases, existence/regularity thresholds, independent actual-quantile derivatives, covariance and determinant properties, and new method numerical failure. Before commits run dotnet build -c Release and dotnet test -c Release --no-build on all supported targets. BestFit tests reference the intended Numerics worktree explicitly. New unapproved scientific changes require a decision; failures do not authorize adjusting optimizers, seeds, likelihood policy or goldens beyond this plan.

Implementation details, parameter coordinates, evidence and retained limitations are recorded in [bestfit-robustness-evidence.md](bestfit-robustness-evidence.md).
