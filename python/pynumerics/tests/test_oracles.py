"""Correlation, ranks, regression and kernel density wrappers against scipy, statsmodels and scikit-learn.

Each wrapper calls one Numerics member; the oracle computes the same quantity independently, so
the tolerances cover only the different operation order.
"""

from __future__ import annotations

import numpy as np
import pytest

from pynumerics import data, distributions, enums, statistics
from pynumerics.errors import net_errors

stats = pytest.importorskip("scipy.stats")
sm = pytest.importorskip("statsmodels.api")
SklearnKernelDensity = pytest.importorskip("sklearn.neighbors").KernelDensity


def _paired_sample(n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    return x, 0.6 * x + rng.normal(size=n)


def test_ranks_equal_scipy_average_ranks_with_ties_and_infinities() -> None:
    # Arrange
    values = np.array([3.0, 1.0, 2.0, 2.0, np.inf, np.inf, 5.0, -np.inf])
    original = values.copy()
    # Act
    ranks = statistics.ranks(values)
    same, ties = statistics.ranks(values, return_ties=True)
    # Assert
    np.testing.assert_array_equal(ranks, stats.rankdata(values, method="average"))
    np.testing.assert_array_equal(same, ranks)
    np.testing.assert_array_equal(ties, [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
    np.testing.assert_array_equal(values, original)


@pytest.mark.parametrize("n", [20, 300])
def test_correlations_equal_scipy(n: int) -> None:
    # Arrange: n = 300 takes the O(n log n) Kendall count, n = 20 the pair loop
    x, y = _paired_sample(n, seed=n)
    y[::7] = y[0]  # ties in one sample
    # Act
    pearson = statistics.pearson(x, y)
    spearman = statistics.spearman(x, y)
    tau = statistics.kendalls_tau(x, y)
    # Assert
    assert pearson == pytest.approx(stats.pearsonr(x, y).statistic, abs=1e-14)
    assert spearman == pytest.approx(stats.spearmanr(x, y).statistic, abs=1e-14)
    assert tau == pytest.approx(stats.kendalltau(x, y).statistic, abs=1e-14)


def test_correlation_matrices_equal_numpy_and_scipy() -> None:
    # Arrange
    rng = np.random.default_rng(11)
    samples = rng.normal(size=(60, 4))
    samples[:, 1] += samples[:, 0]
    # Act
    pearson = statistics.pearson(samples)
    spearman = statistics.spearman(samples)
    # Assert
    np.testing.assert_allclose(pearson, np.corrcoef(samples, rowvar=False), atol=1e-14)
    np.testing.assert_allclose(spearman, stats.spearmanr(samples).statistic, atol=1e-14)


def test_kendalls_tau_is_nan_for_nan_or_constant_input_and_refuses_unequal_lengths() -> None:
    # Act / Assert
    assert np.isnan(statistics.kendalls_tau([1.0, np.nan, 3.0], [1.0, 2.0, 3.0]))
    assert np.isnan(statistics.kendalls_tau([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]))
    with pytest.raises(ValueError), net_errors():
        statistics.kendalls_tau([1.0, 2.0], [1.0, 2.0, 3.0])


def test_linear_regression_equals_statsmodels_ols() -> None:
    # Arrange
    rng = np.random.default_rng(3)
    x = rng.normal(size=(40, 2))
    y = 1.0 + x @ [2.0, -1.0] + rng.normal(size=40)
    new = rng.normal(size=(3, 2))
    ols = sm.OLS(y, sm.add_constant(x)).fit()
    frame = ols.get_prediction(sm.add_constant(new, has_constant="add")).summary_frame(alpha=0.1)
    # Act
    fit = data.LinearRegression(x, y)
    intervals = fit.prediction_intervals(new, alpha=0.1)
    # Assert
    np.testing.assert_allclose(fit.parameters, ols.params, atol=1e-13)
    np.testing.assert_allclose(fit.parameter_standard_errors, ols.bse, atol=1e-13)
    np.testing.assert_allclose(fit.parameter_t_stats, ols.tvalues, atol=1e-12)
    np.testing.assert_allclose(fit.covariance, ols.normalized_cov_params, atol=1e-14)
    np.testing.assert_allclose(fit.residuals, ols.resid, atol=1e-13)
    assert fit.standard_error == pytest.approx(np.sqrt(ols.scale), abs=1e-14)
    assert fit.r_squared == pytest.approx(ols.rsquared, abs=1e-14)
    assert fit.adjusted_r_squared == pytest.approx(ols.rsquared_adj, abs=1e-14)
    assert (fit.sample_size, fit.degrees_of_freedom, fit.has_intercept) == (40, 37, True)
    assert fit.parameter_names[0] == "Intercept" and len(fit.parameter_names) == 3
    np.testing.assert_allclose(fit.predict(new), frame["mean"], atol=1e-13)
    np.testing.assert_allclose(intervals[:, 0], frame["obs_ci_lower"], atol=1e-12)
    np.testing.assert_allclose(intervals[:, 1], frame["obs_ci_upper"], atol=1e-12)
    np.testing.assert_array_equal(intervals[:, 2], fit.predict(new))
    assert any("Intercept" in line for line in fit.summary())


def test_linear_regression_takes_a_single_predictor_as_a_1d_array_and_no_intercept() -> None:
    # Arrange
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y = np.array([2.1, 3.9, 6.2, 7.8, 10.1])
    # Act
    with_intercept = data.LinearRegression(x, y)
    through_origin = data.LinearRegression(x, y, has_intercept=False)
    # Assert
    np.testing.assert_allclose(with_intercept.parameters, np.polyfit(x, y, 1)[::-1], atol=1e-13)
    assert through_origin.parameters[0] == pytest.approx(x @ y / (x @ x), abs=1e-14)
    assert not through_origin.has_intercept
    with pytest.raises(ArithmeticError), net_errors():
        data.LinearRegression([1.0, 2.0], [1.0, 2.0])


def test_gaussian_kernel_density_equals_scipy_gaussian_kde() -> None:
    # Arrange
    sample = np.random.default_rng(5).normal(size=50)
    x = np.array([-1.0, 0.0, 0.5, 2.5])
    # Act
    kde = distributions.KernelDensity(sample)
    oracle = stats.gaussian_kde(sample, bw_method=kde.bandwidth / np.std(sample, ddof=1))
    # Assert
    assert kde.bandwidth == pytest.approx((4.0 / (3.0 * 50)) ** 0.2 * np.std(sample, ddof=1), rel=1e-15)
    assert kde.bandwidth == kde.bandwidth_rule(sample)
    np.testing.assert_allclose(kde.pdf(x), oracle(x), atol=1e-15)
    assert (kde.kernel, kde.sample_size, kde.bounded_by_data) == ("Gaussian", 50, True)
    np.testing.assert_array_equal(kde.sample_data, sample)


@pytest.mark.parametrize(("kernel", "sklearn_kernel"), [
    ("epanechnikov", "epanechnikov"), ("triangular", "linear"), ("uniform", "tophat"),
])
def test_compact_kernel_densities_equal_scikit_learn(kernel: str, sklearn_kernel: str) -> None:
    # Arrange
    sample = np.random.default_rng(6).normal(size=40)
    x = np.array([-1.3, 0.0, 0.4, 1.9])
    oracle = SklearnKernelDensity(kernel=sklearn_kernel, bandwidth=0.7).fit(sample[:, None])
    # Act
    kde = distributions.KernelDensity(sample, kernel, bandwidth=0.7)
    # Assert
    np.testing.assert_allclose(kde.pdf(x), np.exp(oracle.score_samples(x[:, None])), atol=1e-15)


def test_weighted_kernel_density_equals_the_weighted_gaussian_mixture() -> None:
    # Arrange
    sample = np.array([0.0, 1.0, 3.0])
    weights = np.array([1.0, 2.0, 1.0])
    x = np.array([-0.5, 1.0, 2.0])
    expected = (weights * stats.norm.pdf((x[:, None] - sample) / 0.8)).sum(axis=1) / (weights.sum() * 0.8)
    # Act
    kde = distributions.KernelDensity(sample, bandwidth=0.8, weights=weights)
    default = distributions.KernelDensity(sample, weights=weights)
    # Assert
    np.testing.assert_allclose(kde.pdf(x), expected, atol=1e-15)
    mean = weights @ sample / weights.sum()
    sd = np.sqrt(weights @ (sample - mean) ** 2 / weights.sum())
    assert default.bandwidth == pytest.approx((4.0 / 9.0) ** 0.2 * sd, rel=1e-15)


def test_kernel_density_cdf_table_is_close_to_the_exact_mixture_cdf_and_round_trips() -> None:
    # Arrange: the CDF interpolates a 1000-bin table, so it agrees with the exact CDF to the table resolution
    sample = np.random.default_rng(8).normal(size=30)
    kde = distributions.KernelDensity(sample, bandwidth=0.5)
    kde.bounded_by_data = False
    x = np.array([-1.0, 0.0, 1.0])
    # Act
    cdf = kde.cdf(x)
    # Assert
    np.testing.assert_allclose(cdf, stats.norm.cdf((x[:, None] - sample) / 0.5).mean(axis=1), atol=5e-3)
    np.testing.assert_allclose(kde.cdf(kde.inverse_cdf(cdf)), cdf, atol=1e-9)
    assert kde.minimum == pytest.approx(sample.min() - 1.5, abs=1e-15)


def test_kernel_density_setters_clone_and_invalid_bandwidth() -> None:
    # Arrange
    kde = distributions.KernelDensity([1.0, 2.0, 4.0], "gaussian", 0.5)
    # Act
    kde.kernel = enums.KernelType.EPANECHNIKOV
    kde.bandwidth = 2.0
    kde.x_transform = "none"
    kde.probability_transform = "none"
    copy = kde.clone()
    # Assert
    assert isinstance(copy, distributions.KernelDensity)
    assert (copy.kernel, copy.bandwidth, copy.probability_transform) == ("Epanechnikov", 2.0, "None")
    assert isinstance(distributions.wrap(kde.net), distributions.KernelDensity)
    with pytest.raises(ValueError), net_errors():
        kde.bandwidth = 0.0
    with pytest.raises(ValueError), net_errors():
        distributions.KernelDensity([1.0, np.nan])
