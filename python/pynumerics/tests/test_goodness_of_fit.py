"""Goodness-of-fit wrappers against scipy, scikit-learn and the textbook formulas.

Each wrapper calls one Numerics member; the oracle computes the same quantity independently.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest

from pynumerics import distributions, statistics
from pynumerics import goodness_of_fit as gof

stats = pytest.importorskip("scipy.stats")
metrics = pytest.importorskip("sklearn.metrics")


def _pair(n: int = 40, seed: int = 6) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    observed = rng.gamma(3.0, 10.0, n) + 5.0
    return observed, observed * 0.9 + rng.normal(0.0, 4.0, n) + 3.0


def test_information_criteria_and_weights_equal_the_formulas() -> None:
    # Arrange
    ll, k, n = -123.4, 3, 50
    aics = np.array([250.0, 252.5, 260.0])
    rmses = np.array([1.5, 2.0, 4.0])
    # Act / Assert
    assert gof.aic(k, ll) == pytest.approx(2 * k - 2 * ll, rel=1e-15)
    assert gof.aicc(n, k, ll) == pytest.approx(2 * k - 2 * ll + 2 * k * (k + 1) / (n - k - 1), rel=1e-15)
    assert gof.bic(n, k, ll) == pytest.approx(k * math.log(n) - 2 * ll, rel=1e-15)
    w = np.exp(-0.5 * (aics - aics.min()))
    np.testing.assert_allclose(gof.aic_weights(aics), w / w.sum(), rtol=1e-15)
    np.testing.assert_allclose(gof.rmse_weights(rmses), rmses**-2 / np.sum(rmses**-2), rtol=1e-15)


def test_error_metrics_equal_scikit_learn() -> None:
    # Arrange
    o, m = _pair()
    # Act / Assert
    assert gof.mse(o, m) == pytest.approx(metrics.mean_squared_error(o, m), rel=1e-13)
    assert gof.mae(o, m) == pytest.approx(metrics.mean_absolute_error(o, m), rel=1e-13)
    assert gof.mape(o, m) == pytest.approx(100.0 * metrics.mean_absolute_percentage_error(o, m), rel=1e-13)
    assert gof.smape(o, m) == pytest.approx(100.0 * np.mean(2 * np.abs(o - m) / (np.abs(o) + np.abs(m))), rel=1e-13)
    assert gof.rmse(o, m) == pytest.approx(math.sqrt(metrics.mean_squared_error(o, m)), rel=1e-13)
    assert gof.rmse(o, m, k=2) == pytest.approx(math.sqrt(np.sum((o - m) ** 2) / (o.size - 2)), rel=1e-13)
    with pytest.raises(ValueError):
        gof.rmse(o, m, k=o.size)
    with pytest.raises(ValueError):
        gof.mse(o, m[:-1])


def test_rmse_against_a_distribution_uses_weibull_quantiles_and_its_parameter_count() -> None:
    # Arrange
    o, _ = _pair()
    model = distributions.Normal(float(o.mean()), float(o.std(ddof=1)))
    x = np.sort(o)
    pp = np.arange(1, o.size + 1) / (o.size + 1.0)
    q = stats.norm.ppf(pp, o.mean(), o.std(ddof=1))
    # Act / Assert
    assert gof.rmse(o, model) == pytest.approx(math.sqrt(np.sum((q - x) ** 2) / (o.size - 2)), rel=1e-12)
    assert gof.rmse(x, model, plotting_positions=pp) == pytest.approx(gof.rmse(o, model), rel=1e-15)
    with pytest.raises(ValueError):
        gof.rmse(o, o, plotting_positions=pp)


def test_hydrologic_efficiencies_equal_their_formulas() -> None:
    # Arrange
    o, m = _pair()
    om = o.mean()
    r = np.corrcoef(o, m)[0, 1]
    alpha, beta = m.std(ddof=1) / o.std(ddof=1), m.mean() / om
    gamma = (m.std(ddof=1) / m.mean()) / (o.std(ddof=1) / om)
    sum_abs = np.sum(np.abs(m - o))
    eps = om / 100.0
    lo, lm = np.log(o + eps), np.log(m + eps)
    # Act / Assert
    assert gof.nash_sutcliffe_efficiency(o, m) == pytest.approx(metrics.r2_score(o, m), rel=1e-12)
    assert gof.log_nash_sutcliffe_efficiency(o, m) == pytest.approx(metrics.r2_score(lo, lm), rel=1e-12)
    assert gof.log_nash_sutcliffe_efficiency(o, m, 0.0) == pytest.approx(
        metrics.r2_score(np.log(o), np.log(m)), rel=1e-12)
    assert gof.kling_gupta_efficiency(o, m) == pytest.approx(
        1 - math.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2), rel=1e-12)
    assert gof.kling_gupta_efficiency_mod(o, m) == pytest.approx(
        1 - math.sqrt((r - 1) ** 2 + (gamma - 1) ** 2 + (beta - 1) ** 2), rel=1e-12)
    assert gof.kling_gupta_efficiency(o, np.full(o.size, 3.0)) == -10.0
    assert gof.pbias(o, m) == pytest.approx(100.0 * np.sum(o - m) / np.sum(o), rel=1e-12)  # Moriasi (2007)
    assert gof.pbias([1.0, 2.0, 3.0, 4.0], [0.8, 1.6, 2.4, 3.2]) == pytest.approx(20.0, rel=1e-12)
    assert gof.pbias([1.0, 2.0, 3.0, 4.0], [2.0, 2.0, 4.0, 4.0]) == pytest.approx(-20.0, rel=1e-12)
    assert gof.rsr(o, m) == pytest.approx(math.sqrt(np.sum((o - m) ** 2) / np.sum((o - om) ** 2)), rel=1e-12)
    assert gof.r_squared(o, m) == pytest.approx(stats.pearsonr(o, m).statistic ** 2, rel=1e-12)
    assert gof.index_of_agreement(o, m) == pytest.approx(
        1 - np.sum((m - o) ** 2) / np.sum((np.abs(m - om) + np.abs(o - om)) ** 2), rel=1e-12)
    assert gof.modified_index_of_agreement(o, m) == pytest.approx(
        1 - sum_abs / np.sum(np.abs(m - om) + np.abs(o - om)), rel=1e-12)
    two_dev = 2 * np.sum(np.abs(o - om))
    refined = 1 - sum_abs / two_dev if sum_abs <= two_dev else two_dev / sum_abs - 1
    assert gof.refined_index_of_agreement(o, m) == pytest.approx(refined, rel=1e-12)
    assert gof.volumetric_efficiency(o, m) == pytest.approx(1 - sum_abs / np.sum(np.abs(o)), rel=1e-12)


def test_classification_scores_equal_scikit_learn() -> None:
    # Arrange
    rng = np.random.default_rng(12)
    truth = (rng.random(200) < 0.3).astype(float)
    predicted = np.where(rng.random(200) < 0.8, truth, 1.0 - truth)
    tn, fp, fn, tp = metrics.confusion_matrix(truth, predicted).ravel()
    # Act / Assert
    assert gof.accuracy(truth, predicted) == pytest.approx(100.0 * metrics.accuracy_score(truth, predicted))
    assert gof.precision(truth, predicted) == pytest.approx(metrics.precision_score(truth, predicted))
    assert gof.recall(truth, predicted) == pytest.approx(metrics.recall_score(truth, predicted))
    assert gof.f1_score(truth, predicted) == pytest.approx(metrics.f1_score(truth, predicted))
    assert gof.specificity(truth, predicted) == pytest.approx(tn / (tn + fp))
    assert gof.balanced_accuracy(truth, predicted) == pytest.approx(metrics.balanced_accuracy_score(truth, predicted))
    assert tp + fn > 0


def test_distribution_statistics_equal_scipy() -> None:
    # Arrange
    rng = np.random.default_rng(21)
    sample = rng.gumbel(100.0, 20.0, 75)
    model = distributions.Gumbel(100.0, 20.0)
    frozen = stats.gumbel_r(100.0, 20.0)
    x = np.sort(sample)
    n = x.size
    i = np.arange(1, n + 1)
    a2 = -n - np.sum((2 * i - 1) * (frozen.logcdf(x) + frozen.logsf(x[::-1]))) / n
    lower, upper, counts = statistics.Histogram(x).bins()
    expected = n * (frozen.cdf(upper) - frozen.cdf(lower))
    # Act
    ks = gof.kolmogorov_smirnov(sample, model)
    # Assert
    assert ks == pytest.approx(stats.kstest(sample, frozen.cdf).statistic, rel=1e-12)
    assert gof.kolmogorov_smirnov(x, model, data_is_sorted=True) == ks
    assert gof.anderson_darling(sample, model) == pytest.approx(a2, rel=1e-11)
    assert gof.chi_squared(sample, model) == pytest.approx(np.sum((counts - expected) ** 2 / expected), rel=1e-11)
    with pytest.raises(TypeError):
        gof.kolmogorov_smirnov(sample, frozen)
    with pytest.raises(ValueError):
        gof.kolmogorov_smirnov([], model)


def test_paired_metrics_on_large_series_equal_sklearn() -> None:
    # Arrange
    o, m = _pair(500_000, 3)
    # Act
    nse = gof.nash_sutcliffe_efficiency(o, m)
    # Assert
    assert nse == pytest.approx(metrics.r2_score(o, m), rel=1e-10)


@pytest.mark.benchmark
def test_paired_metrics_cross_once_for_large_series() -> None:
    # Arrange
    o, m = _pair(500_000, 3)
    gof.nash_sutcliffe_efficiency(o[:10], m[:10])
    # Act
    start = time.perf_counter()
    gof.nash_sutcliffe_efficiency(o, m)
    elapsed = time.perf_counter() - start
    # Assert: bulk copies, not one crossing per value
    assert elapsed < 0.5


def test_column_uniform_statistics_equal_scipy_and_the_definition() -> None:
    # Arrange
    rng = np.random.default_rng(3)
    u = np.sort(rng.random((25, 12)), axis=0)
    # Act
    ks = gof.column_kolmogorov_smirnov_uniform(u)
    ad = gof.column_anderson_darling_uniform(u)
    # Assert
    uniform = distributions.Uniform(0.0, 1.0)
    for j in range(u.shape[1]):
        n = u.shape[0]
        i = np.arange(1, n + 1)
        expected = -n - np.sum((2 * i - 1) * (np.log(u[:, j]) + np.log1p(-u[::-1, j]))) / n
        assert ks[j] == pytest.approx(stats.kstest(u[:, j], "uniform").statistic, abs=1e-15)
        assert ad[j] == pytest.approx(expected, rel=1e-13)
        assert ks[j] == gof.kolmogorov_smirnov(u[:, j], uniform, data_is_sorted=True)
        assert ad[j] == gof.anderson_darling(u[:, j], uniform, data_is_sorted=True)
    for bad in (u[::-1], np.where(u > 0.9, np.nan, u), u + 0.5):
        with pytest.raises(ValueError):
            gof.column_kolmogorov_smirnov_uniform(bad)
