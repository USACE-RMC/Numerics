"""Hypothesis-test wrappers against scipy, statsmodels, R and published flood-frequency examples.

Each wrapper calls one Numerics member; the oracle computes the same p-value independently, so the
tolerances cover only the different operation order or a published value's rounding.
"""

from __future__ import annotations

import math
import time
from fractions import Fraction

import numpy as np
import pytest

from pynumerics import hypothesis_tests as ht

stats = pytest.importorskip("scipy.stats")
sm = pytest.importorskip("statsmodels.api")
diagnostic = pytest.importorskip("statsmodels.stats.diagnostic")

# Table 1.2 of Bobee and Ashkar (1991), The Gamma Family and Derived Distributions Applied in
# Hydrology: annual maximum discharge (m3/s), Harricana River at Amos, Quebec.
HARRICANA = np.array([
    122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235,
    216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117,
    192, 337, 125, 166, 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187,
    180, 173, 174,
])


def _samples(seed: int = 3) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    return rng.normal(10.0, 2.0, 25), rng.normal(11.5, 3.0, 31)


def test_t_tests_equal_scipy() -> None:
    # Arrange
    a, b = _samples()
    c = a + np.random.default_rng(4).normal(0.4, 0.5, a.size)
    # Act / Assert
    assert ht.one_sample_t_test(a, 9.5) == pytest.approx(stats.ttest_1samp(a, 9.5).pvalue, rel=1e-10)
    assert ht.one_sample_t_test(a - a.mean()) == pytest.approx(1.0, abs=1e-12)
    assert ht.equal_variance_t_test(a, b) == pytest.approx(stats.ttest_ind(a, b).pvalue, rel=1e-10)
    assert ht.unequal_variance_t_test(a, b) == pytest.approx(stats.ttest_ind(a, b, equal_var=False).pvalue, rel=1e-10)
    assert ht.paired_t_test(a, c) == pytest.approx(stats.ttest_rel(a, c).pvalue, rel=1e-10)


def test_f_tests_equal_the_f_distribution_and_statsmodels() -> None:
    # Arrange
    a, b = _samples()
    f = np.var(b, ddof=1) / np.var(a, ddof=1)
    two_sided = 2.0 * min(stats.f.cdf(f, b.size - 1, a.size - 1), stats.f.sf(f, b.size - 1, a.size - 1))
    x = np.arange(30.0)
    y = 1.0 + 0.3 * x + 0.02 * x**2 + np.random.default_rng(5).normal(0.0, 1.0, x.size)
    full = sm.OLS(y, sm.add_constant(np.column_stack([x, x**2]))).fit()
    restricted = sm.OLS(y, sm.add_constant(x)).fit()
    f_ref, p_ref, _ = full.compare_f_test(restricted)
    # Act
    f_stat, p_value = ht.f_test_models(restricted.ssr, full.ssr, int(restricted.df_resid), int(full.df_resid))
    # Assert
    assert ht.f_test(a, b) == pytest.approx(two_sided, rel=1e-10)
    assert f_stat == pytest.approx(f_ref, rel=1e-12)
    assert p_value == pytest.approx(p_ref, rel=1e-9)
    assert ht.f_test_models(10.0, 12.0, 20, 18) == (pytest.approx(-1.5), 1.0)
    with pytest.raises(ValueError):
        ht.f_test_models(40.0, 100.0, 7, 9)
    with pytest.raises(ValueError):
        ht.f_test([1.0], [1.0, 2.0])


def test_normality_and_autocorrelation_tests_equal_scipy_and_statsmodels() -> None:
    # Arrange
    rng = np.random.default_rng(8)
    sample = rng.gamma(2.0, 1.0, 120)
    series = np.cumsum(rng.normal(size=80)) * 0.1 + rng.normal(size=80)
    lb = diagnostic.acorr_ljungbox(series, lags=[6])
    # Act / Assert
    assert ht.jarque_bera_test(sample) == pytest.approx(stats.jarque_bera(sample).pvalue, rel=1e-10)
    assert ht.ljung_box_test(series, 6) == pytest.approx(float(lb["lb_pvalue"].iloc[0]), rel=1e-10)
    default_lag = math.floor(min(10 * math.log10(series.size), series.size - 1))
    lb_default = diagnostic.acorr_ljungbox(series, lags=[default_lag])
    assert ht.ljung_box_test(series) == pytest.approx(float(lb_default["lb_pvalue"].iloc[0]), rel=1e-10)


def test_wald_wolfowitz_mann_whitney_and_mann_kendall_match_published_values() -> None:
    # Arrange: Bobee and Ashkar (1991) report U = 1.167 (Wald-Wolfowitz) and |z| = 0.54
    # (Mann-Whitney, the last 19 years against the first 50); Mann-Kendall p = 0.7757.
    first, second = HARRICANA[50:], HARRICANA[:50]
    # Act
    ww = ht.wald_wolfowitz_test(HARRICANA)
    mw = ht.mann_whitney_test(first, second)
    mk = ht.mann_kendall_test(HARRICANA)
    # Assert
    assert ww == pytest.approx(2.0 * stats.norm.sf(1.167), abs=1e-3)
    assert mw == pytest.approx(2.0 * stats.norm.sf(0.54), abs=1e-2)
    assert mw == pytest.approx(stats.mannwhitneyu(first, second, method="asymptotic",
                                                  use_continuity=False).pvalue, rel=1e-10)
    assert mk == pytest.approx(0.7757, abs=1e-4)
    with pytest.raises(ValueError):
        ht.mann_whitney_test(second, first)
    with pytest.raises(ValueError):
        ht.mann_kendall_test(HARRICANA[:9])


def test_mann_kendall_equals_the_textbook_formula_with_ties() -> None:
    # Arrange
    x = np.array([3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0, 5.0, 8.0, 9.0, 7.0, 9.0])
    n = x.size
    s = sum(np.sign(x[j] - x[i]) for i in range(n - 1) for j in range(i + 1, n))
    _, counts = np.unique(x, return_counts=True)
    var = (n * (n - 1) * (2 * n + 5) - np.sum(counts * (counts - 1) * (2 * counts + 5))) / 18.0
    z = (s - np.sign(s)) / math.sqrt(var)
    # Act / Assert
    assert ht.mann_kendall_test(x) == pytest.approx(2.0 * stats.norm.sf(abs(z)), rel=1e-12)


def test_linear_trend_test_equals_statsmodels_ols() -> None:
    # Arrange
    t = np.arange(1.0, 41.0)
    y = 0.05 * t + np.random.default_rng(9).normal(size=t.size)
    fit = sm.OLS(y, sm.add_constant(t)).fit()
    # Act / Assert
    assert ht.linear_trend_test(t, y) == pytest.approx(fit.pvalues[1], rel=1e-9)
    with pytest.raises(ValueError):
        ht.linear_trend_test(t[:-1], y)


def test_dip_equals_r_diptest_and_unimodality_matches_mclust() -> None:
    # Arrange: R 4.4.3, diptest::dip on the same values
    bimodal = [1.2, 2.3, 1.9, 2.8, 3.1, 7.9, 8.4, 9.1, 8.8, 7.5, 2.2, 8.1]
    unimodal = [-1.5, -0.8, -0.3, 0.1, 0.4, 0.9, 1.3, 2.0, -0.1, 0.6]
    # mclust likelihood-ratio p-values, from the Numerics reference tests
    flat = [4.5, 5.2, 5.1, 4.9, 5.0, 5.3, 5.4, 4.8, 4.7, 5.2, 5.1, 4.6, 5.0, 5.3, 5.1, 4.9, 5.2, 5.0, 4.8, 5.3,
            4.9, 5.1, 5.2, 4.8, 5.0, 5.1, 5.2, 4.7, 5.3, 5.0]
    # Act / Assert
    assert ht.hartigan_dip(bimodal) == pytest.approx(0.17460317460317457, abs=1e-15)
    assert ht.hartigan_dip(unimodal) == pytest.approx(0.060000000000000012, abs=1e-15)
    p = ht.hartigan_dip_test(bimodal)
    assert p < 0.01  # R dip.test interpolates its table: 0.00044
    assert ht.hartigan_dip_test(bimodal, 500, 7) == ht.hartigan_dip_test(bimodal, 500, 7)
    assert ht.hartigan_dip_test([1.0, 2.0, 3.0]) == 1.0
    assert ht.unimodality_test(flat) == pytest.approx(0.4142441, abs=1e-4)
    boot = ht.unimodality_test(flat, 200, 3)
    assert 0.0 <= boot <= 1.0 and boot == ht.unimodality_test(flat, 200, 3)
    with pytest.raises(ValueError):
        ht.hartigan_dip([1.0, np.nan])
    with pytest.raises(ValueError):
        ht.unimodality_test(flat[:9])


def test_grubbs_beck_thresholds_match_the_published_example_and_the_formula() -> None:
    # Arrange: Bobee and Ashkar (1991), page 9: high 378.2, low 91.2
    logs = np.log(HARRICANA)
    n = HARRICANA.size
    k = -3.62201 + 6.28446 * n**0.25 - 2.49835 * n**0.5 + 0.491436 * n**0.75 - 0.037911 * n
    # Act
    high, low = ht.grubbs_beck_test(HARRICANA)
    # Assert
    assert high == pytest.approx(378.2, abs=0.1) and low == pytest.approx(91.2, abs=0.1)
    assert high == pytest.approx(math.exp(logs.mean() + k * logs.std(ddof=1)), rel=1e-13)
    with pytest.raises(ValueError):
        ht.grubbs_beck_test([1.0, 0.0, 3.0])


def test_multiple_grubbs_beck_matches_hec_ssp() -> None:
    # Arrange: HEC-SSP Bulletin 17C analyses (HCK: one low outlier; LOP: one), and the Numerics
    # reference sample with none.
    hck = [23520, 16220, 15290, 14750, 13980, 13740, 13030, 12930, 12550, 12140, 11470, 10610, 10220, 10060, 10010,
           9880, 9090, 8780, 8700, 8420, 8410, 8280, 8140, 7990, 7770, 7540, 7460, 7400, 7400, 6770, 6590, 6500, 6420,
           6200, 6190, 6160, 6100, 5900, 5800, 5730, 5690, 5670, 5630, 5620, 5590, 5580, 5560, 5440, 5410, 5340, 5330,
           5310, 5150, 5080, 4690, 4680, 4620, 4600, 4580, 4490, 4380, 4110, 4080, 4050, 4020, 4010, 3980, 3960, 3900,
           3840, 3770, 3700, 3670, 3620, 3300, 3220, 2820, 2810, 2670, 2440, 2370, 2030, 1750, 1570, 920]
    rng = np.random.default_rng(2)
    lognormal = np.exp(rng.normal(6.0, 0.4, 60))
    with_lows = np.concatenate([lognormal, [5.0, 8.0, 10.0]])
    # Act / Assert
    assert ht.multiple_grubbs_beck_test(hck) == 1
    assert ht.multiple_grubbs_beck_test(hck[::-1]) == 1
    assert ht.multiple_grubbs_beck_test(with_lows) >= 3
    assert isinstance(ht.multiple_grubbs_beck_test(lognormal), int)


def test_hot_path_tests_on_large_samples_equal_scipy() -> None:
    # Arrange
    rng = np.random.default_rng(1)
    a, b = rng.normal(size=200_000), rng.normal(0.01, 1.0, 200_000)
    # Act
    p = ht.equal_variance_t_test(a, b)
    # Assert
    assert p == pytest.approx(stats.ttest_ind(a, b).pvalue, rel=1e-8)


@pytest.mark.benchmark
def test_hot_path_tests_cross_once_for_large_samples() -> None:
    # Arrange
    rng = np.random.default_rng(1)
    a, b = rng.normal(size=200_000), rng.normal(0.01, 1.0, 200_000)
    ht.equal_variance_t_test(a[:10], b[:10])
    # Act
    start = time.perf_counter()
    ht.equal_variance_t_test(a, b)
    elapsed = time.perf_counter() - start
    # Assert: one bulk copy each, not one crossing per value (a per-element path takes seconds)
    assert elapsed < 0.5


def test_two_sample_kolmogorov_smirnov_equals_scipy() -> None:
    # Arrange
    rng = np.random.default_rng(8)
    a, b = rng.normal(size=70), rng.normal(0.3, 1.2, 55)
    ties_a, ties_b = np.round(a), np.round(b)
    # Act / Assert
    assert ht.two_sample_kolmogorov_smirnov(a, b) == pytest.approx(stats.ks_2samp(a, b).statistic, abs=1e-15)
    assert ht.two_sample_kolmogorov_smirnov(ties_a, ties_b) == pytest.approx(
        stats.ks_2samp(ties_a, ties_b).statistic, abs=1e-15)
    assert ht.two_sample_kolmogorov_smirnov(np.sort(a), np.sort(b), data_is_sorted=True) == \
        ht.two_sample_kolmogorov_smirnov(a, b)
    with pytest.raises(ValueError, match="sorted"):
        ht.two_sample_kolmogorov_smirnov([2.0, 1.0], [1.0], data_is_sorted=True)
    with pytest.raises(ValueError):
        ht.two_sample_kolmogorov_smirnov([], [1.0])


def _exact_ks(a: list[float], b: list[float]) -> Fraction:
    """Return max |F_a - F_b| over the pooled values as an exact fraction, by counting."""
    return max(abs(Fraction(sum(v <= t for v in a), len(a)) - Fraction(sum(v <= t for v in b), len(b)))
               for t in set(a) | set(b))


def test_two_sample_kolmogorov_smirnov_is_the_correctly_rounded_rational() -> None:
    # Arrange: tie-heavy quarter-step samples of unequal sizes, plus the documented small cases
    rng = np.random.default_rng(20261107)
    cases = [([0.0, 1.0, 1.0], [1.0, 1.0, 1.0]), ([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]),
             ([1.0, 1.0, 2.0], [1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0]), ([2.0, 2.0], [5.0, 5.0, 5.0]),
             ([-1.0, 0.0, 0.0, 3.0], [-1.0, 0.0, 0.0, 3.0])]
    for _ in range(300):
        levels = int(rng.integers(1, 13))
        cases.append(((rng.integers(0, levels + 1, int(rng.integers(1, 61))) / 4.0).tolist(),
                      (rng.integers(0, levels + 1, int(rng.integers(1, 61))) / 4.0).tolist()))
    for a, b in cases:
        # Act
        value = ht.two_sample_kolmogorov_smirnov(a, b)
        # Assert: the nearest double to the exact rational, and scipy to an ulp or so (scipy versions
        # differ on whether the statistic is rounded from the integer count)
        assert value == float(_exact_ks(a, b))
        assert value == pytest.approx(float(stats.ks_2samp(a, b).statistic), rel=1e-14, abs=1e-15)
    assert ht.two_sample_kolmogorov_smirnov([0.0, 1.0, 1.0], [1.0, 1.0, 1.0]) == 1.0 / 3.0
    assert ht.two_sample_kolmogorov_smirnov([2.0, 2.0], [5.0, 5.0, 5.0]) == 1.0
    assert ht.two_sample_kolmogorov_smirnov([-1.0, 0.0, 3.0], [-1.0, 0.0, 3.0]) == 0.0
