"""Descriptive statistics, transformations, histograms, streaming covariance and given-data sensitivity.

Wrappers are checked against NumPy and scipy (``scipy.stats``), and against the analytic value of a
sensitivity index where the estimator's finite-sample bias is known to be small.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest

from pynumerics import callbacks, distributions, statistics

stats = pytest.importorskip("scipy.stats")


def _sample(n: int = 60, seed: int = 5) -> np.ndarray:
    return np.random.default_rng(seed).gamma(2.5, 4.0, n) + 1.0


def test_descriptive_statistics_equal_numpy_and_scipy() -> None:
    # Arrange
    x = _sample()
    y = x * 0.5 + np.random.default_rng(6).normal(size=x.size)
    # Act / Assert
    assert statistics.minimum(x) == x.min() and statistics.maximum(x) == x.max()
    assert statistics.sum(x) == pytest.approx(x.sum(), rel=1e-14)
    assert statistics.parallel_mean(x) == pytest.approx(x.mean(), rel=1e-14)
    assert statistics.geometric_mean(x) == pytest.approx(stats.gmean(x), rel=1e-13)
    assert statistics.harmonic_mean(x) == pytest.approx(stats.hmean(x), rel=1e-13)
    assert statistics.population_variance(x) == pytest.approx(np.var(x), rel=1e-13)
    assert statistics.population_standard_deviation(x) == pytest.approx(np.std(x), rel=1e-13)
    mean, var = statistics.mean_variance(x)
    assert mean == pytest.approx(x.mean(), rel=1e-14) and var == pytest.approx(np.var(x, ddof=1), rel=1e-13)
    mean, sd = statistics.mean_standard_deviation(x)
    assert sd == pytest.approx(np.std(x, ddof=1), rel=1e-13)
    assert statistics.coefficient_of_variation(x) == pytest.approx(np.std(x, ddof=1) / x.mean(), rel=1e-13)
    assert statistics.skewness(x) == pytest.approx(stats.skew(x, bias=False), rel=1e-12)
    assert statistics.kurtosis(x) == pytest.approx(stats.kurtosis(x, bias=False), rel=1e-12)
    assert statistics.skewness(x, np.ones(x.size)) == pytest.approx(statistics.skewness(x), rel=1e-12)
    assert statistics.kurtosis(x, np.ones(x.size), "reliability") == pytest.approx(statistics.kurtosis(x), rel=1e-10)
    assert statistics.covariance(x, y) == pytest.approx(np.cov(x, y)[0, 1], rel=1e-13)
    assert statistics.population_covariance(x, y) == pytest.approx(np.cov(x, y, ddof=0)[0, 1], rel=1e-13)
    moments = statistics.product_moments(x)
    np.testing.assert_allclose(moments, [x.mean(), np.std(x, ddof=1), stats.skew(x, bias=False),
                                         stats.kurtosis(x, bias=False)], rtol=1e-12)
    assert np.isnan(statistics.product_moments([1.0, 2.0, 3.0])).all()
    np.testing.assert_allclose(statistics.five_number_summary(x), np.percentile(x, [0, 25, 50, 75, 100]), rtol=1e-14)
    np.testing.assert_allclose(statistics.seven_number_summary(x),
                               np.percentile(x, [0, 5, 25, 50, 75, 95, 100]), rtol=1e-14)
    assert np.isnan(statistics.minimum([1.0, np.nan]))


def test_jackknife_and_entropy_with_python_and_csharp_callbacks() -> None:
    # Arrange
    x = _sample(30)
    leave_one_out = np.array([np.delete(x, i).mean() for i in range(x.size)])
    se = math.sqrt((x.size - 1) / x.size * np.sum((leave_one_out - leave_one_out.mean()) ** 2))
    normal = distributions.Normal(0.0, 1.0)
    points = np.linspace(-2.0, 2.0, 9)
    density = stats.norm.pdf(points)
    # Act / Assert
    assert statistics.jackknife_standard_error(x, np.mean) == pytest.approx(se, rel=1e-12)
    np.testing.assert_allclose(statistics.jackknife_sample(x, np.mean), leave_one_out, rtol=1e-13)
    assert statistics.jackknife_sample([], np.mean) is None
    if callbacks.roslyn_available():
        compiled = statistics.jackknife_standard_error(x, "x => x.Average()")
        assert compiled == pytest.approx(se, rel=1e-12)
    assert statistics.entropy(points, normal.pdf) == pytest.approx(-np.sum(density * np.log(density)), rel=1e-12)


def test_effective_sample_size_is_the_sequential_kish_ratio() -> None:
    # Arrange: [1, 1, 2] gives 8/3 bitwise
    w = np.random.default_rng(3).random(50)
    total, squares = 0.0, 0.0
    for value in w:
        total += value
        squares += value * value
    # Act / Assert
    assert statistics.effective_sample_size([1.0, 1.0, 2.0]) == 2.6666666666666665
    assert statistics.effective_sample_size(w) == total * total / squares
    assert statistics.effective_sample_size(np.ones(7)) == 7.0
    assert np.isnan(statistics.effective_sample_size([]))
    with pytest.raises(ValueError):
        statistics.effective_sample_size([1.0, -1.0])


def test_box_cox_and_yeo_johnson_equal_scipy() -> None:
    # Arrange
    x = _sample()
    signed = np.random.default_rng(9).normal(1.0, 3.0, 80)
    lam = 0.37
    # Act / Assert: scipy's llf omits the constant -n/2 (1 + ln 2 pi), so differences in lambda agree
    assert statistics.box_cox_log_jacobian(x, lam) == pytest.approx((lam - 1) * np.sum(np.log(x)), rel=1e-13)
    diff = statistics.box_cox_log_likelihood(x, lam) - statistics.box_cox_log_likelihood(x, 1.2)
    assert diff == pytest.approx(stats.boxcox_llf(lam, x) - stats.boxcox_llf(1.2, x), rel=1e-10)
    success, fitted = statistics.box_cox_fit_lambda(x)
    assert success and fitted == pytest.approx(stats.boxcox_normmax(x, method="mle"), abs=1e-4)
    np.testing.assert_allclose(statistics.yeo_johnson_transform(signed, lam), stats.yeojohnson(signed, lam),
                               rtol=1e-13)
    assert statistics.yeo_johnson_transform(-2.0, lam) == pytest.approx(stats.yeojohnson(np.array([-2.0]), lam)[0])
    back = statistics.yeo_johnson_inverse_transform(statistics.yeo_johnson_transform(signed, lam), lam)
    np.testing.assert_allclose(back, signed, rtol=1e-12, atol=1e-12)
    h = 1e-6
    numeric = (stats.yeojohnson(signed + h, lam) - stats.yeojohnson(signed - h, lam)) / (2 * h)
    np.testing.assert_allclose(statistics.yeo_johnson_derivative(signed, lam), numeric, rtol=1e-6)
    jacobian = (lam - 1) * np.sum(np.sign(signed) * np.log1p(np.abs(signed)))
    assert statistics.yeo_johnson_log_jacobian(signed, lam) == pytest.approx(jacobian, rel=1e-12)
    diff = statistics.yeo_johnson_log_likelihood(signed, lam) - statistics.yeo_johnson_log_likelihood(signed, 1.1)
    assert diff == pytest.approx(stats.yeojohnson_llf(lam, signed) - stats.yeojohnson_llf(1.1, signed), rel=1e-9)
    assert statistics.yeo_johnson_fit_lambda(signed) == pytest.approx(stats.yeojohnson_normmax(signed), abs=1e-4)
    with pytest.raises(ValueError):
        statistics.yeo_johnson_fit_lambda([1.0])


def test_running_statistics_population_moments_and_clone() -> None:
    # Arrange
    x = _sample(40)
    rs = statistics.RunningStatistics(x)
    # Act
    copy = rs.clone()
    copy.push(1e6)
    # Assert
    assert rs.population_variance == pytest.approx(np.var(x), rel=1e-12)
    assert rs.population_standard_deviation == pytest.approx(np.std(x), rel=1e-12)
    assert rs.coefficient_of_variation == pytest.approx(np.std(x, ddof=1) / x.mean(), rel=1e-12)
    assert rs.population_skewness == pytest.approx(stats.skew(x), rel=1e-10)
    assert rs.population_kurtosis == pytest.approx(stats.kurtosis(x), rel=1e-10)
    assert rs.count == x.size and copy.count == x.size + 1


def test_autocorrelation_validates_its_arguments_and_gives_the_band() -> None:
    # Arrange
    series = np.sin(np.arange(20.0))
    # Act / Assert
    for bad, match in [(np.ones((4, 2)), "data must be one-dimensional"), ([1.0], "data needs at least two"),
                       ([1.0, np.nan, 2.0], "data must be finite"), ([1.0, np.inf], "data must be finite")]:
        with pytest.raises(ValueError, match=match):
            statistics.autocorrelation(bad)
    for lag in (0, 20, -1, 2.0, True):
        with pytest.raises(ValueError, match="max_lag"):
            statistics.autocorrelation(series, lag)
    lags, _ = statistics.autocorrelation(series, np.int64(19))
    assert lags[-1] == 19
    lower, upper = statistics.autocorrelation_confidence_interval(100)
    z = stats.norm.ppf(0.975) / 10.0
    assert lower == pytest.approx(-z, rel=1e-9) and upper == pytest.approx(z, rel=1e-9)


def test_uncentered_autocorrelation_equals_numpy_and_median_propagates_nan() -> None:
    # Arrange
    x = np.array([0.3, -1.2, 0.8, 0.5, -0.7, 1.1, -0.4, -0.9, 0.6, 0.2])
    expected = np.array([np.dot(x[k:], x[:x.size - k]) for k in range(4)])
    # Act
    _, cov = statistics.autocorrelation(x, 3, "covariance", centered=False)
    lags, acf = statistics.autocorrelation(x, 3, centered=False)
    _, centered = statistics.autocorrelation(x, 3, centered=True)
    # Assert
    assert list(lags) == [0, 1, 2, 3]
    np.testing.assert_allclose(cov, expected / x.size, rtol=0, atol=1e-15)
    np.testing.assert_allclose(acf, expected / expected[0], rtol=0, atol=1e-15)
    assert np.array_equal(centered, statistics.autocorrelation(x, 3)[1])
    assert math.isnan(statistics.median([1.0, np.nan, 3.0]))
    assert statistics.median([1e308, 1e308]) == 1e308


def test_histogram_counts_equal_numpy_and_bins_are_values() -> None:
    # Arrange
    x = _sample(200)
    hist = statistics.Histogram(x, 8)
    edges = np.linspace(x.min(), x.max(), 9)
    # Act
    lower, upper, counts = hist.bins()
    hist.add_data(x.max() + 5.0)
    hist.add_data(np.array([x.min() - 1.0, x.mean()]))
    # Assert
    np.testing.assert_array_equal(counts, np.histogram(x, edges)[0])
    np.testing.assert_allclose(lower, edges[:-1], rtol=1e-12)
    np.testing.assert_allclose(upper, edges[1:], rtol=1e-12)
    assert len(hist) == hist.number_of_bins == 8
    assert hist.data_count == 203 and hist.upper_bound == x.max() + 5.0 and hist.lower_bound == x.min() - 1.0
    assert hist.bin_width == pytest.approx((x.max() - x.min()) / 8)
    mids = (lower + upper) / 2
    plain = statistics.Histogram(x, 8)
    assert plain.mean == pytest.approx(np.sum(mids * counts) / counts.sum(), rel=1e-12)
    assert plain.standard_deviation == pytest.approx(
        math.sqrt(np.sum(counts * (mids - plain.mean) ** 2) / counts.sum()), rel=1e-12)
    assert plain.mode == mids[np.argmax(counts)]
    assert plain.median == mids[np.searchsorted(np.cumsum(counts), counts.sum() / 2)]
    assert plain.bin_index(edges[3] + 1e-9) == 3
    assert plain[-1].upper_bound == x.max()
    with pytest.raises(IndexError):
        plain[8]
    with pytest.raises(ValueError):
        plain.bin_index(x.max() + 100.0)
    assert statistics.Histogram(x).number_of_bins == math.ceil(2 * x.size ** (1 / 3)) + 1


def test_histogram_bins_compare_clone_and_mutate() -> None:
    # Arrange
    first, second = statistics.HistogramBin(0.0, 1.0, 3), statistics.HistogramBin(1.0, 2.0)
    # Act
    copy = first.clone()
    copy.frequency = 7
    copy.upper_bound = 1.5
    copy.lower_bound = 0.5
    # Assert
    assert first.frequency == 3 and copy.frequency == 7 and copy.midpoint == 1.0
    assert first.compare_to(second) == -1 and second.compare_to(first) == 1 and first.compare_to(None) == 1
    assert first == statistics.HistogramBin(0.0, 1.0, 3) and first != second
    with pytest.raises(TypeError):
        hash(first)
    hist = statistics.Histogram([0.2, 0.4, 1.6], 2)
    hist.add_bin(statistics.HistogramBin(5.0, 6.0))
    assert len(hist) == 2 and hist.bins()[0].size == 2  # number_of_bins is the construction count


def test_running_covariance_matrix_equals_numpy_plus_the_identity_prior() -> None:
    # Arrange
    x = np.random.default_rng(4).multivariate_normal([1.0, -2.0, 0.5], [[2, 0.6, 0], [0.6, 1, 0.3], [0, 0.3, 3]], 400)
    rcm = statistics.RunningCovarianceMatrix(3)
    n = x.shape[0]
    m2 = (x - x.mean(0)).T @ (x - x.mean(0))
    # Act
    rcm.push(x[:10])
    rcm.push(x[10])
    rcm.push(x[11:])
    # Assert
    assert rcm.n == n
    np.testing.assert_allclose(rcm.mean, x.mean(0), rtol=1e-12)
    np.testing.assert_allclose(rcm.covariance, m2 + np.eye(3), rtol=1e-10)
    np.testing.assert_allclose(rcm.sample_covariance, (m2 + np.eye(3)) / (n - 1), rtol=1e-10)
    np.testing.assert_allclose(rcm.population_covariance, (m2 + np.eye(3)) / n, rtol=1e-10)
    d = np.sqrt(np.diag(rcm.sample_covariance))
    np.testing.assert_allclose(rcm.sample_correlation, rcm.sample_covariance / np.outer(d, d), rtol=1e-10)
    d = np.sqrt(np.diag(rcm.population_covariance))
    np.testing.assert_allclose(rcm.population_correlation, rcm.population_covariance / np.outer(d, d), rtol=1e-10)
    with pytest.raises(ValueError):
        rcm.push([1.0, 2.0])


def test_given_data_sensitivity_indices_match_scipy_and_analytic_values() -> None:
    # Arrange: y = x1 + x2 + small noise with independent uniform inputs has first-order index 1/2
    rng = np.random.default_rng(10)
    n, bins = 20_000, 20
    x1, x2 = rng.random(n), rng.random(n)
    y = x1 + x2 + 0.01 * rng.normal(size=n)
    order = np.lexsort((np.arange(n), x1))
    edges = [b * n // bins for b in range(bins + 1)]
    reference = [stats.ks_2samp(y[order[edges[b]:edges[b + 1]]], y).statistic for b in range(bins)]
    # Act
    per_bin = statistics.pawn(x1, y, bins)
    # Assert
    assert statistics.first_order_sobol(x1, y) == pytest.approx(0.5, abs=0.01)
    np.testing.assert_allclose(per_bin, reference, atol=1e-15)
    assert statistics.pawn_median(x1, y, bins) == pytest.approx(np.median(per_bin), abs=1e-15)
    assert statistics.borgonovo_delta(x1, x1) == pytest.approx(1 - 1 / 20, abs=1e-12)
    assert statistics.borgonovo_delta(x1, rng.random(n)) < 0.1  # finite-sample bias ~ sqrt(classes / bin size)
    assert statistics.borgonovo_delta(x1, y) == statistics.borgonovo_delta(x1, np.exp(y))  # rank invariance
    with pytest.raises(ValueError):
        statistics.first_order_sobol(x1[:10], y[:10])


def test_bulk_statistics_on_large_arrays_equal_numpy_and_scipy() -> None:
    # Arrange
    x = np.random.default_rng(1).normal(size=500_000)
    # Act
    summary = statistics.seven_number_summary(x)
    skew = statistics.skewness(x)
    ess = statistics.effective_sample_size(np.abs(x))
    # Assert
    assert summary[3] == pytest.approx(np.median(x), rel=1e-12)
    assert skew == pytest.approx(stats.skew(x, bias=False), rel=1e-9) and ess > 0


@pytest.mark.benchmark
def test_bulk_statistics_cross_once_for_large_arrays() -> None:
    # Arrange
    x = np.random.default_rng(1).normal(size=500_000)
    statistics.skewness(x[:10])
    # Act
    start = time.perf_counter()
    statistics.seven_number_summary(x)
    statistics.skewness(x)
    statistics.effective_sample_size(np.abs(x))
    elapsed = time.perf_counter() - start
    # Assert: bulk copies, not one crossing per value
    assert elapsed < 1.0


def test_column_dispersion_and_linear_moments_equal_numpy_and_one_sample_members() -> None:
    # Arrange
    rng = np.random.default_rng(11)
    data = 1e6 + 50.0 * rng.random((60, 30))
    data[rng.random(data.shape) < 0.05] = np.nan
    data[3, 4] = np.inf
    finite = np.where(np.isfinite(data), data, np.nan)
    # Act
    variances, excluded = statistics.column_variances(data)
    deviations, _ = statistics.column_standard_deviations(data)
    moments, lexcluded = statistics.column_linear_moments(data)
    # Assert
    np.testing.assert_allclose(variances, np.nanvar(finite, axis=0, ddof=1), rtol=1e-12)
    np.testing.assert_allclose(deviations, np.nanstd(finite, axis=0, ddof=1), rtol=1e-12)
    assert np.array_equal(excluded, (~np.isfinite(data)).sum(axis=0)) and np.array_equal(excluded, lexcluded)
    assert moments.shape == (4, 30)
    for j in range(30):
        column = data[np.isfinite(data[:, j]), j]
        assert variances[j] == statistics.variance(column)
        assert deviations[j] == statistics.standard_deviation(column)
        assert np.array_equal(moments[:, j], statistics.linear_moments(column))
        assert moments[1, j] / moments[0, j] == statistics.linear_coefficient_of_variation(column)


def test_weighted_hazen_percentile_matches_the_definition() -> None:
    # Arrange
    rng = np.random.default_rng(14)
    x = np.sort(rng.gamma(2.0, 3.0, 40))
    w = rng.random(40)
    w[::7] = 0.0
    total = 0.0
    positions = []
    for weight in w:
        positions.append((total + 0.5 * weight) / sum(w.tolist()))
        total += weight
    # Act
    got = statistics.plotting_position_percentile(x, [0.1, 0.5, 0.9], "hazen", data_is_sorted=True, weights=w)
    # Assert
    np.testing.assert_allclose(got, np.interp([0.1, 0.5, 0.9], positions, x), rtol=1e-12)
    unsorted = statistics.plotting_position_percentile(x[::-1], 0.5, "hazen", weights=w[::-1])
    assert unsorted == pytest.approx(got[1], rel=1e-14)
    shuffled = rng.permutation(40)
    batch = statistics.plotting_position_percentile(x[shuffled], [0.1, 0.5, 0.9], "hazen", weights=w[shuffled])
    singles = [statistics.plotting_position_percentile(x[shuffled], q, "hazen", weights=w[shuffled])
               for q in (0.1, 0.5, 0.9)]
    assert batch.tolist() == singles


def _exact_linear_moments(x: np.ndarray, nmom: int) -> list:
    """Exact rational unbiased L-moments (Hosking 1990) of the double data, as Fractions."""
    from fractions import Fraction  # noqa: PLC0415

    xs = sorted(Fraction(float(v)) for v in x)
    n = len(xs)
    b = [sum(xs[i - 1] * math.comb(i - 1, r) for i in range(1, n + 1)) / (n * math.comb(n - 1, r)) for r in range(nmom)]
    lam = [sum(Fraction((-1) ** (r - k) * math.factorial(r + k), math.factorial(k) ** 2 * math.factorial(r - k)) * b[k]
               for k in range(r + 1)) for r in range(nmom)]
    return lam[:2] + [lam[r] / lam[1] for r in range(2, nmom)]


def test_linear_moments_with_nmom_match_exact_rationals_and_the_four_value_member() -> None:
    # Arrange
    rng = np.random.default_rng(23)
    worst = [0.0] * 6
    for _ in range(40):
        x = rng.lognormal(3.0, 1.0, int(rng.integers(6, 150)))
        # Act
        got = statistics.linear_moments(x, nmom=6)
        exact = _exact_linear_moments(x, 6)
        # Assert
        assert np.array_equal(got[:4], statistics.linear_moments(x))
        for j in range(6):
            error = abs(float(got[j]) - float(exact[j]))
            worst[j] = max(worst[j], error / abs(float(exact[j])) if j < 2 else error)
    assert worst[0] < 1e-15 and worst[1] < 1e-14
    assert max(worst[2:4]) < 1e-13 and worst[4] < 1e-12 and worst[5] < 1e-11


def test_linear_moments_with_nmom_edges() -> None:
    assert statistics.linear_moments([1.0, 3.0], nmom=2).tolist() == [2.0, 1.0]
    assert np.isnan(statistics.linear_moments([1.0, 2.0, 3.0, 4.0], nmom=5)).all()
    assert statistics.linear_moments([1.0, 2.0, 3.0, 4.0], nmom=5).shape == (5,)
    with pytest.raises(ValueError):
        statistics.linear_moments([1.0, 2.0, 3.0], nmom=0)
    with pytest.raises(TypeError):
        statistics.linear_moments([1.0, 2.0, 3.0, 4.0, 5.0], nmom=5.0)
    assert statistics.linear_moments([1.0, 3.0], nmom=np.int64(2)).tolist() == [2.0, 1.0]


def test_weighted_percentile_list_overload_is_bitwise_the_scalar_member() -> None:
    # Arrange
    rng = np.random.default_rng(41)
    x = np.round(rng.gamma(2.0, 3.0, 60), 1)
    w = rng.random(60)
    w[::9] = 0.0
    k = np.array([0.9, 0.1, 0.5, 0.5, 0.0, 1.0, 0.37, 0.02])
    for kind in ("weibull", "median", "blom", "cunnane", "gringorten", "hazen", "kaplan_meier"):
        # Act
        batch = statistics.plotting_position_percentile(x, k, kind, weights=w)
        singles = [statistics.plotting_position_percentile(x, float(q), kind, weights=w) for q in k]
        # Assert
        assert batch.tolist() == singles
    grid = statistics.plotting_position_percentile(x, k.reshape(2, 4), "hazen", weights=w)
    assert grid.shape == (2, 4)
    with pytest.raises(ValueError):
        statistics.plotting_position_percentile(x, [0.5, 1.5], "hazen", weights=w)


def test_blom_and_hazen_percentiles_are_closer_to_exact_rationals() -> None:
    # Arrange: the exact rational interpolation at the return-period probabilities 1 - 1/T.
    from fractions import Fraction  # noqa: PLC0415

    rng = np.random.default_rng(12)
    for kind, alpha in (("blom", Fraction(3, 8)), ("hazen", Fraction(1, 2))):
        new_error = old_error = 0.0
        for _ in range(300):
            n = int(rng.integers(5, 121))
            x = np.sort(np.exp(3.0 + 2.0 * (rng.random(n) - 0.5)))
            k = 1.0 - 1.0 / int(rng.integers(2, 101))
            h_double = (n + 1 - 2.0 * float(alpha)) * k + float(alpha)
            h = (n + 1 - 2 * alpha) * Fraction(k) + alpha
            r = math.floor(h)
            if not 1 <= r < n or math.floor(h_double) != r:
                continue
            exact = Fraction(float(x[r - 1])) + (h - r) * (Fraction(float(x[r])) - Fraction(float(x[r - 1])))
            # Act
            got = statistics.plotting_position_percentile(x, k, kind, data_is_sorted=True)
            rank_weight = float(x[r - 1]) + (h_double - r) * (float(x[r]) - float(x[r - 1]))
            # Assert (accumulated in ulps of the exact value)
            unit = math.ulp(float(exact))
            new_error += abs(float((Fraction(got) - exact) / Fraction(unit)))
            old_error += abs(float((Fraction(rank_weight) - exact) / Fraction(unit)))
        assert new_error < old_error, kind
