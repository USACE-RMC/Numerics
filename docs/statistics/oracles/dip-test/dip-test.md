# Hartigan dip test oracle fixtures

Generated with **R version 4.4.3 (2025-02-28 ucrt)** and the CRAN package **diptest 0.77-2**
(Maechler and Ringach, 2025). diptest is GPL-licensed; the generator calls it as a black box to
produce reference values and nothing from it is copied into Numerics. The Numerics implementation
follows Hartigan (1985), Algorithm AS 217, and simulates p-values from its own uniform samples.

Run from the repository root:

```powershell
& 'C:/Program Files/R/R-4.4.3/bin/Rscript.exe' --vanilla docs/statistics/oracles/dip-test/dip-test.R
```

An optional first argument changes the CSV output path. An optional second argument also writes
the 400-sample randomized validation set (seed 20261105; n from 1 to 2,000; normal, two-component,
rounded, exponential, uniform and rounded-mixture samples) as one line per sample, the dip followed
by the values. Numerics matches every dip in that set to within 1.4e-16 relative. The set is not
committed because of its size. The generator checks the R and diptest versions, and regeneration is
deterministic, with LF line endings and 17-digit numeric text.

## Contents

| Case | Quantity | Oracle | Test |
|---|---|---|---|
| empty, single, two, four, constant5 | dip | `diptest::dip` | `Test_HartiganDip_MatchesR` |
| ties10, ties15, halfsteps20 | dip of tied samples | `diptest::dip` | `Test_HartiganDip_MatchesR` |
| mixture20, bimodal30 | dip | `diptest::dip` | `Test_HartiganDip_MatchesR` |
| bimodal30 | p-value from the interpolated table, and simulated with B = 100,000 | `diptest::dip.test` | `Test_HartiganDipTest_MatchesR` |
| quantiles40 | minimum dip 1/(2n) and p = 1 | `diptest::dip`, `dip.test` | `Test_HartiganDipTest_MatchesR` |
| three | p = 1 for n = 3 | `diptest::dip.test` | `Test_HartiganDipTest_MatchesR` |

Numerics p-values are Monte Carlo estimates (2,000 realizations by default), so the tests compare
them with the R p-values within Monte Carlo error rather than exactly.

## References

Hartigan, J. A., & Hartigan, P. M. (1985). The dip test of unimodality. *The Annals of Statistics*, 13(1), 70-84.

Hartigan, P. M. (1985). Algorithm AS 217: Computation of the dip statistic to test for unimodality. *Journal of the Royal Statistical Society: Series C (Applied Statistics)*, 34, 320-325.

Maechler, M., & Ringach, D. (2025). diptest: Hartigan's Dip Test Statistic for Unimodality - Corrected. R package version 0.77-2. https://CRAN.R-project.org/package=diptest
