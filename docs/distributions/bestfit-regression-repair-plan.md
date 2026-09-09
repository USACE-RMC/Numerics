# Iterative BestFit distribution regression repair

Approved by Haden Smith on 2026-09-08. This is an execution plan; scoped regression repairs and each exact method in the companion allowlist are authorized.

## Objective and frozen contracts

Eliminate numerical and reproducible runtime regressions attributable to the recent Numerics distribution hardening and BestFit integration. Restore the previous `GetParameterConstraints()` logic, initialization, family-specific magnitude rounding, and valid prior envelopes. Exceptional-input handling must remain separate from previously working cases. Apply complete covariance transformations into the verification tests' existing coordinates, retaining cross-terms and recovery rules.

Do not change tolerances, acceptance thresholds, confidence levels, reference values, seeds, sample sizes, sampler settings, optimizer settings, or default scientific policies. Do not skip failing tests. A conflict requiring a new scientific policy goes to Haden for a specific decision.

## Source versions and isolation

- Numerics comparison: `d80bfa8621c48a78cf4ebad7f01326841fac37aa`; repair start: `a3bd2afd64286e1b8df3ada8be711d55180cbb11`.
- BestFit comparison: `fbe0989a87a0bd38a297f8a42e8d979585d0f4ad`; repair start: `3e0371a995b2c9d7993f8b4ddd42b48238b9bb7a`.
- Earlier Kappa comparison if needed: `94d1713e03a4c416dd472d6fbef0ce907881e9d1` to `c0d67b9c52d62dc49bdbfb744ec47b9a205cc86b`.
- Both repositories have isolated repair and detached baseline worktrees under `artifacts/worktrees/bestfit-regression-{repair,baseline}`. Original dirty checkouts are preserved.
- Resolve local Numerics explicitly through `UseLocalRmcNumerics=true` and `RmcNumericsProjectPath`; record source and assembly fingerprints.

## Allowlist and iteration

The companion inventory freezes 293 exact catalog identities: 50 distribution-fitting, 40 estimation/diagnostic, 117 univariate, 25 bivariate, 16 rating-curve, 22 time-series, and 23 spatial. It records 35 unrelated exclusions and 56 source methods absent in the original dirty BestFit checkout. Reappearing coverage sources in a comparison worktree do not expand authorization.

1. Run one allowlisted fully qualified method through BestFit's guarded runner, serially. Independently verify one TRX, one result, and the exact requested class/method.
2. Record result, test duration separately from build/host time, source versions/fingerprints, dependencies, and exact TRX path.
3. On failure or material slowdown, immediately compare the unchanged method with the relevant baseline; isolate Numerics versus integration effects and inspect inputs, initialization, bounds, likelihoods, estimates, diagnostics, and computational hot paths.
4. Fix attributable regressions immediately, add a focused regression check, rerun the method, and mark earlier dependent passes stale. Continue independent work if a case needs a scientific decision.
5. Begin with Pearson and Log-Pearson covariance failures and confirmed constraint regressions; then order by dependencies and observed failures.

## Tasks and validation

1. Freeze the allowlist, dirty-file fingerprints, and resumable ledger; reproduce the two covariance failures against the pinned baseline.
2. Restore legacy constraints and initialization for previously valid inputs, preserving readable family-specific rounding; retain robust exceptional-input fallback. Add baseline-derived checks and verify resulting BestFit Uniform priors and nonstationary bounds.
3. Correct full covariance coordinate transformations for Pearson, Log-Pearson, and LnNormal verification helpers. Keep existing numerical acceptance unchanged and validate cross-term propagation independently.
4. Execute the iterative 293-method loop; investigate all affected density/tail/interval/quantile/support/moment/derivative/cache/composite paths as evidence requires.
5. After repairs stabilize, obtain a fresh complete 293-method pass. Run Numerics Release and framework gates, BestFit fast and XML documentation gates, inspect prohibited-setting and unrelated-work diffs, and commit validated changes locally. Do not push.

The ledger is authoritative for progress. A partial run is not completion. Any change to production dependencies invalidates earlier affected evidence until rerun.
