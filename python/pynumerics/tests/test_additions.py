"""Wrappers over the Numerics statistics, time series and distribution members they call directly."""

from __future__ import annotations

import importlib

import numpy as np
import pytest

from pynumerics import distributions, statistics
from pynumerics.marshal import to_net_doubles


def test_linear_coefficient_of_variation_is_bitwise_the_l_moment_ratio() -> None:
    # Arrange
    sample = np.exp(np.random.default_rng(4).normal(size=60))
    moments = statistics.linear_moments(sample)
    # Act
    lcv = statistics.linear_coefficient_of_variation(sample)
    # Assert
    assert lcv == moments[1] / moments[0]
    assert np.isnan(statistics.linear_coefficient_of_variation([1.0, 2.0, 3.0]))


def test_box_cox_transform_derivative_and_fit_equal_scipy() -> None:
    # Arrange
    from scipy import stats  # noqa: PLC0415

    sample = np.exp(np.random.default_rng(9).normal(size=80))
    # Act
    transformed = statistics.box_cox_transform(sample, 0.3)
    back = statistics.box_cox_inverse_transform(transformed, 0.3)
    derivative = statistics.box_cox_derivative(sample, 0.3)
    success, lam = statistics.box_cox_fit_lambda(sample)
    failed, nan = statistics.box_cox_fit_lambda([2.0, 2.0, 2.0])
    # Assert
    np.testing.assert_allclose(transformed, stats.boxcox(sample, 0.3), rtol=1e-14, atol=1e-15)
    assert statistics.box_cox_transform(float(sample[0]), 0.3) == transformed[0]
    np.testing.assert_allclose(back, sample, rtol=1e-13)
    np.testing.assert_allclose(derivative, sample ** -0.7, rtol=1e-14)
    assert statistics.box_cox_derivative(2.0, 0.0) == 0.5
    assert np.isnan(statistics.box_cox_derivative(-1.0, 0.3))
    assert success and lam == pytest.approx(stats.boxcox(sample)[1], abs=1e-6)
    assert not failed and np.isnan(nan)


def test_running_statistics_clear_matches_a_new_accumulator_bitwise() -> None:
    # Arrange
    values = [1.5, -2.0, 4.25, 8.0]
    used = statistics.RunningStatistics([100.0, 200.0])
    # Act
    used.clear()
    empty_count = used.count
    used.push(values)
    fresh = statistics.RunningStatistics(values)
    # Assert
    assert empty_count == 0
    for name in ("count", "minimum", "maximum", "mean", "variance", "skewness", "kurtosis"):
        assert getattr(used, name) == getattr(fresh, name), name


def test_time_series_ticks_and_values_equal_the_ordinates() -> None:
    # Arrange
    import System  # noqa: PLC0415

    from pynumerics import data  # noqa: PLC0415

    net_data = importlib.import_module("Numerics.Data")
    start = System.DateTime(2001, 2, 28, 6, 0, 0)
    values = np.array([1.5, np.nan, -0.0, 4.0])
    series = net_data.TimeSeries(net_data.TimeInterval.OneDay, start, to_net_doubles(values))
    # Act
    ticks = data.time_series_index_ticks(series)
    read = data.time_series_values(series)
    # Assert
    assert ticks.dtype == np.int64
    assert ticks.tolist() == [int(series[i].Index.Ticks) for i in range(series.Count)]
    assert ticks[1] - ticks[0] == 864_000_000_000  # one day in 100 ns ticks
    dates = ((ticks - 621_355_968_000_000_000) // 10).astype("datetime64[us]")
    assert str(dates[0]) == "2001-02-28T06:00:00.000000"
    np.testing.assert_array_equal(read.view(np.uint64), values.view(np.uint64))


def test_column_percentiles_equal_the_per_column_percentile_bitwise() -> None:
    # Arrange
    rng = np.random.default_rng(7)
    data = 100.0 * np.exp(rng.standard_normal((300, 25)))
    data[rng.random(data.shape) < 0.05] = np.nan
    data[3, 4] = np.inf
    data[:, 7] = np.nan
    k = np.array([0.05, 0.25, 0.5, 0.75, 0.95])
    # Act
    table, excluded = statistics.column_percentiles(data, k)
    median, _ = statistics.column_percentiles(data, 0.5)
    # Assert
    assert table.shape == (5, 25) and excluded.dtype == np.int32
    for column in range(data.shape[1]):
        values = data[:, column]
        finite = values[np.isfinite(values)]
        assert excluded[column] == values.size - finite.size
        if finite.size == 0:
            assert np.all(np.isnan(table[:, column]))
            continue
        reference = statistics.percentile(finite, k)
        np.testing.assert_array_equal(table[:, column].view(np.uint64), reference.view(np.uint64))
    np.testing.assert_array_equal(median, table[2])


def test_generalized_pareto_fit_known_location_matches_the_reference_maximum() -> None:
    # Arrange
    excesses = np.array([1.2, 0.4, 3.9, 2.2, 0.15, 5.6, 1.05, 0.8, 7.3, 2.9, 0.33, 1.7, 4.4, 0.62, 2.5, 11.2, 0.09,
                         1.33, 3.1, 0.95])
    # Act
    fitted = distributions.GeneralizedPareto.fit_known_location(excesses + 10.0, 10.0)
    # Assert (the reference maximizes the same likelihood with SciPy at xatol 1e-12)
    location, scale, shape = fitted.parameters
    assert location == 10.0
    assert scale == pytest.approx(2.388710825197821, rel=1e-3)
    assert shape == pytest.approx(-0.07721168825622311, abs=1e-3)
    assert fitted.log_likelihood(excesses + 10.0) == pytest.approx(-38.959310020261206, abs=1e-6)


def test_weighted_plotting_positions_reduce_to_the_unweighted_formula() -> None:
    # Arrange
    weights = np.array([0.5, 2.0, 0.0, 1.25, 3.0])
    # Act
    positions = statistics.weighted_plotting_positions(weights, 0.44)
    unit = statistics.weighted_plotting_positions(np.ones(40), 0.44)
    # Assert
    cumulative = np.cumsum(weights)
    total = cumulative[-1]
    expected = (cumulative - 0.44 * weights) / (total + (1.0 - 2.0 * 0.44) * (total / weights.size))
    np.testing.assert_allclose(positions, expected, rtol=1e-15, atol=0.0)
    np.testing.assert_array_equal(unit.view(np.uint64), statistics.plotting_positions(40, 0.44).view(np.uint64))


def test_wilson_score_interval_known_values_and_operation_order() -> None:
    # Arrange
    import math

    # Act
    lower, upper = statistics.wilson_score_interval(5, 20, 0.9)
    # Assert (SciPy binomtest(5, 20).proportion_ci(0.9, method="wilson"))
    assert lower == pytest.approx(0.12737716852035774, abs=1e-12)
    assert upper == pytest.approx(0.4322017501601245, abs=1e-12)
    z = float(distributions.Normal(0.0, 1.0).inverse_cdf(0.5 + 0.8 / 2.0))
    p_hat = 12.5 / 37.25
    z2n = z * z / 37.25
    denominator = 1.0 + z2n
    center = (p_hat + z2n / 2.0) / denominator
    half = z * math.sqrt(p_hat * (1.0 - p_hat) / 37.25 + z2n / (4.0 * 37.25)) / denominator
    assert statistics.wilson_score_interval(12.5, 37.25, 0.8) == (max(0.0, center - half), min(1.0, center + half))
    assert statistics.wilson_score_interval(0, 10, 0.95)[0] == pytest.approx(0.0, abs=1e-15)


def test_weibull_plotting_position_percentile_agrees_with_numpy_type_6() -> None:
    # Arrange
    rng = np.random.default_rng(11)
    data = rng.gamma(2.0, 30.0, 47)
    k = np.linspace(1 / 48, 47 / 48, 200)
    # Act
    weibull = statistics.plotting_position_percentile(data, k, "weibull")
    scalar = statistics.plotting_position_percentile(data, 0.9, "weibull")
    hazen = statistics.plotting_position_percentile([3, 1, 4, 1, 5, 9, 2, 6], 0.9, "hazen")
    # Assert: NumPy forms the rank as n q + q - 1, Numerics as (n + 1) q, so the two agree to rounding.
    np.testing.assert_allclose(weibull, np.quantile(data, k, method="weibull"), rtol=1e-13, atol=0.0)
    assert scalar == statistics.plotting_position_percentile(data, [0.9], "weibull")[0]
    assert scalar == pytest.approx(np.quantile(data, 0.9, method="weibull"), rel=1e-13)
    assert hazen == pytest.approx(8.1, abs=1e-12)


def test_weighted_mean_curve_reproduces_the_unweighted_curve_at_equal_weights() -> None:
    # Arrange
    from pynumerics.bootstrap import UncertaintyAnalysisResults

    members = [distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1),
               distributions.GeneralizedExtremeValue(110.0, 25.0, 0.0),
               distributions.GeneralizedExtremeValue(95.0, 15.0, 0.1)]
    p = [0.01, 0.5, 0.99]
    # Act
    plain = UncertaintyAnalysisResults().process_mean_curve(members, p)
    equal = UncertaintyAnalysisResults().process_mean_curve(members, p, weights=[0.2, 0.2, 0.2])
    weighted = UncertaintyAnalysisResults().process_mean_curve(members, p, weights=[0.6, 0.3, 0.1])
    # Assert
    np.testing.assert_array_equal(equal.view(np.uint64), plain.view(np.uint64))
    assert not np.array_equal(weighted, plain)


def test_the_added_wrappers_refuse_invalid_input_as_value_errors() -> None:
    # Arrange
    from pynumerics.bootstrap import UncertaintyAnalysisResults
    from pynumerics.errors import net_errors

    members = [distributions.Normal(0.0, 1.0), distributions.Normal(1.0, 1.0)]
    calls = [
        lambda: statistics.column_percentiles(np.ones((3, 2)), [1.5]),
        lambda: distributions.GeneralizedPareto.fit_known_location([1.0, 2.0, -1.0], 0.0),
        lambda: statistics.weighted_plotting_positions([1.0, -1.0], 0.5),
        lambda: statistics.weighted_plotting_positions([1.0, 1.0], 0.75),
        lambda: statistics.wilson_score_interval(11, 10, 0.9),
        lambda: statistics.plotting_position_percentile([1.0, 2.0], 1.5, "weibull"),
        lambda: UncertaintyAnalysisResults().process_mean_curve(members, [0.5], weights=[1.0, -1.0]),
    ]
    # Act / Assert
    for call in calls:
        with pytest.raises(ValueError), net_errors():
            call()
