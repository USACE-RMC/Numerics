"""TimeSeries wrapper: construction, arithmetic, smoothing, block series, resampling and statistics.

pandas and NumPy are the oracles (``rolling``, ``resample``, ``groupby``, ``cumsum``, ``diff``,
``interpolate``, ``DateOffset`` and ``numpy.percentile``); every value comes from Numerics.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
import pytest

from pynumerics import timeseries
from pynumerics.timeseries import BlockFunction, TimeInterval, TimeSeries


def _daily(values, start="2020-01-01") -> TimeSeries:
    return TimeSeries(np.asarray(values, dtype=np.float64), interval="one_day", start=start)


def _gamma(n: int, seed: int = 7) -> np.ndarray:
    return np.random.default_rng(seed).gamma(2.0, 3.0, n)


# Construction and transfer --------------------------------------------------------------------


def test_regular_series_dates_match_pandas_date_range() -> None:
    # Arrange
    values = np.linspace(1.0, 10.0, 40)
    # Act
    series = _daily(values)
    # Assert
    np.testing.assert_array_equal(series.values, values)
    np.testing.assert_array_equal(series.dates, pd.date_range("2020-01-01", periods=40, freq="D").to_numpy())
    assert series.interval is TimeInterval.ONE_DAY
    assert series.start_date == np.datetime64("2020-01-01") and series.end_date == np.datetime64("2020-02-09")
    assert len(series) == series.count == 40
    assert series[2] == (np.datetime64("2020-01-03T00:00:00.000000"), values[2])


_CALENDAR_STEPS = [
    (TimeInterval.ONE_MONTH, "2020-01-31", {"months": 1}),
    (TimeInterval.ONE_QUARTER, "2019-11-30", {"months": 3}),
    (TimeInterval.ONE_YEAR, "2020-02-29", {"years": 1}),
]


def _chained(start: str, offset: dict, n: int) -> np.ndarray:
    """Return pandas ``DateOffset`` applied to the previous date, n - 1 times from the start."""
    dates = [pd.Timestamp(start)]
    for _ in range(n - 1):
        dates.append(dates[-1] + pd.DateOffset(**offset))
    return np.array(dates, dtype="datetime64[us]")


@pytest.mark.parametrize(("interval", "start", "offset"), _CALENDAR_STEPS)
def test_calendar_steps_chain_from_the_previous_date_like_pandas_date_offset(interval, start, offset) -> None:
    # Act
    series = TimeSeries(np.arange(5.0), interval=interval, start=start)
    ranged = TimeSeries(interval=interval, start=start, end=series.end_date)
    # Assert (each date is the previous one plus a step, so a clipped day of month carries on)
    expected = _chained(start, offset, 5)
    np.testing.assert_array_equal(series.dates, expected)
    np.testing.assert_array_equal(ranged.dates, expected)


def test_calendar_steps_chain_in_shift_fill_and_resample() -> None:
    # Arrange: 31 Jan, 29 Feb, 29 Mar, 29 Apr, 29 May (2020)
    chained = _chained("2020-01-31", {"months": 1}, 5)
    series = TimeSeries(np.arange(5.0), interval="one_month", start="2019-01-01")
    sparse = TimeSeries([7.0], interval="one_month", dates=[chained[2]])
    source = TimeSeries(np.sin(np.arange(24.0)), interval="one_month", start=chained[0])
    # Act
    series.shift_all_dates(chained[0])
    filled = sparse.fill_missing_dates(chained[0], chained[-1], -1.0)
    # Assert
    np.testing.assert_array_equal(series.dates, chained)
    np.testing.assert_array_equal(filled.dates, chained)
    np.testing.assert_array_equal(filled.values, [-1.0, -1.0, 7.0, -1.0, -1.0])
    np.testing.assert_array_equal(source.resample_with_knn(5, 3, seed=7).dates, chained)
    np.testing.assert_array_equal(source.resample_with_block_bootstrap(5, 3, seed=7).dates, chained)


def test_start_and_end_series_has_no_date_after_the_end() -> None:
    # Act: chained monthly dates from 31 January fall on the 29th, so 30 June lies between two of them
    series = TimeSeries(interval="one_month", start="2020-01-31", end="2020-06-30")
    # Assert
    np.testing.assert_array_equal(series.dates, _chained("2020-01-31", {"months": 1}, 6))
    assert len(series) == 6 and series.end_date <= np.datetime64("2020-06-30")


@pytest.mark.parametrize(("interval", "start", "offset"), _CALENDAR_STEPS)
def test_anchored_dates_count_from_the_start_like_pandas_date_offset(interval, start, offset) -> None:
    # Act
    dates = TimeSeries.anchored_dates(start, 5, interval)
    series = TimeSeries(np.arange(5.0), interval=interval, dates=dates)
    # Assert (pandas DateOffset applied k times to the start: 31 Jan, 29 Feb, 31 Mar, 30 Apr, 31 May)
    step = pd.DateOffset(**offset)
    expected = np.array([pd.Timestamp(start) + k * step for k in range(5)], dtype="datetime64[us]")
    np.testing.assert_array_equal(dates, expected)
    np.testing.assert_array_equal(series.dates, expected)


def test_anchored_dates_match_a_regular_series_where_they_agree() -> None:
    # Fixed steps, and calendar steps from day 28 or earlier, give the dates of a regular series
    for interval, start in [("six_hour", "2020-01-31T18:30"), ("seven_day", "2020-02-29"), ("one_month", "2019-01-28"),
                            ("one_year", "2019-03-15")]:
        regular = TimeSeries(np.zeros(30), interval=interval, start=start)
        np.testing.assert_array_equal(TimeSeries.anchored_dates(start, 30, interval), regular.dates)
    assert TimeSeries.anchored_dates("2020-01-31", 0, "one_month").size == 0
    with pytest.raises(ValueError):
        TimeSeries.anchored_dates("2020-01-31", 3, "irregular")
    with pytest.raises(ValueError):
        TimeSeries.anchored_dates("2020-01-31", -1, "one_month")


def test_explicit_dates_and_values_keep_every_tick_and_bit() -> None:
    # Arrange
    dates = np.array(["1850-06-01T00:00:00.1234567", "2020-02-29T23:59:59", "2100-01-01", "2100-01-02"],
                     dtype="datetime64[ns]")
    values = np.array([np.uint64(0x7FF8000000001234).view(np.float64), -0.0, np.inf, np.nan])
    # Act
    series = TimeSeries(values, interval="irregular", dates=dates)
    # Assert (the ticks constructor stores every value as given, a NaN payload included)
    np.testing.assert_array_equal(series.ticks, timeseries.to_ticks(dates))
    assert series.ticks[0] % 10 == 7
    np.testing.assert_array_equal(series.values.view(np.uint64), values.view(np.uint64))
    assert series.interval is TimeInterval.IRREGULAR


def test_start_end_constructor_fills_and_empty_series_appends() -> None:
    # Act
    filled = TimeSeries(interval="one_hour", start="2020-01-01", end="2020-01-01T05:00", fill=2.5)
    empty = TimeSeries(interval="one_year")
    empty.append("2015-01-01", 1200.0)
    empty.append("2016-01-01", 1250.0)
    # Assert
    np.testing.assert_array_equal(filled.values, np.full(6, 2.5))
    np.testing.assert_array_equal(empty.values, [1200.0, 1250.0])
    with pytest.raises(ValueError):
        TimeSeries([1.0, 2.0])


def test_pandas_and_xml_round_trips_are_exact() -> None:
    # Arrange
    frame = pd.Series(_gamma(30), index=pd.date_range("2001-03-04 06:00", periods=30, freq="6h"))
    # Act
    series = TimeSeries.from_pandas(frame, interval="six_hour")
    back = TimeSeries.from_xml(series.to_xml())
    # Assert
    pd.testing.assert_series_equal(series.to_pandas(), frame, check_freq=False, check_index_type=False)
    np.testing.assert_array_equal(back.ticks, series.ticks)
    np.testing.assert_array_equal(back.values, series.values)
    assert back.interval is TimeInterval.SIX_HOUR


def test_collection_members_insert_remove_and_index() -> None:
    # Arrange
    series = _daily([1.0, 2.0, 3.0])
    # Act
    series.insert(0, "2019-12-31", 0.5)
    series[1] = ("2020-01-01", 1.5)
    position = series.index_of("2020-01-02", 2.0)
    present = series.contains("2020-01-03", 3.0)
    removed = series.remove("2020-01-03", 3.0)
    series.remove_at(0)
    # Assert
    assert position == 2 and present and removed
    np.testing.assert_array_equal(series.values, [1.5, 2.0])
    assert list(series)[0] == (np.datetime64("2020-01-01"), 1.5)
    series.clear()
    assert len(series) == 0


def test_index_of_date_matches_pandas_get_indexer() -> None:
    # Arrange
    series = _daily(_gamma(400), start="1999-12-30")
    index = pd.DatetimeIndex(series.dates)
    queries = ["1999-12-30", "2000-02-29", "2001-02-02", "1999-12-29", np.datetime64("2000-03-01T00:00:00.0000001"),
               np.datetime64("2001-02-01T00:00:00"), pd.Timestamp("2000-06-15")]
    # Act
    found = [series.index_of_date(query) for query in queries]
    # Assert
    expected = index.get_indexer(pd.DatetimeIndex([pd.Timestamp(q) for q in queries])).tolist()
    assert found == expected
    assert found[0] == 0 and found[1] == 61 and found[3] == -1 and found[4] == -1 and found[5] == 399


def test_index_of_date_on_an_unsorted_series_returns_the_first_match() -> None:
    # Arrange
    dates = np.array(["2020-03-01", "2020-01-01", "2020-02-01", "2020-01-01"], dtype="datetime64[us]")
    series = TimeSeries([1.0, 2.0, 3.0, 4.0], interval="irregular", dates=dates)
    empty = TimeSeries(interval="one_day")
    # Act / Assert
    assert series.index_of_date("2020-01-01") == 1
    assert series.index_of_date("2020-03-01") == 0
    assert series.index_of_date("2020-04-01") == -1
    assert empty.index_of_date("2020-01-01") == -1


def test_index_of_date_crosses_the_dates_once() -> None:
    # Arrange: 100,000 ordinates; one lookup copies the ticks in one transfer (about 0.3 ms measured)
    series = _daily(np.zeros(100_000), start="1800-01-01")
    last = series.dates[-1]
    series.index_of_date(last)
    # Act
    begin = time.perf_counter()
    position = series.index_of_date(last)
    seconds = time.perf_counter() - begin
    # Assert
    assert position == 99_999
    assert seconds < 0.02, seconds


def test_wrapping_a_net_series_shares_it() -> None:
    # Arrange
    series = _daily([1.0, 2.0])
    # Act
    shared = TimeSeries.from_net(series.net)
    shared.add(1.0)
    copied = series.copy()
    copied.add(1.0)
    # Assert
    np.testing.assert_array_equal(series.values, [2.0, 3.0])
    np.testing.assert_array_equal(copied.values, [3.0, 4.0])


# Arithmetic and missing data ------------------------------------------------------------------


def test_in_place_arithmetic_matches_numpy() -> None:
    # Arrange
    values = np.array([10.0, 15.0, 20.0, 25.0, 30.0])
    series = _daily(values)
    expected = values.copy()
    # Act
    series.add(5.0)
    series.subtract(3.0)
    series.multiply(2.0)
    series.divide(4.0)
    series.add(10.0, [0, 2])
    series.multiply(1.5, [1])
    series.exponentiate(2.0)
    series.log_transform(10.0)
    expected = (expected + 5.0 - 3.0) * 2.0 / 4.0
    expected[[0, 2]] += 10.0
    expected[1] *= 1.5
    expected = np.log10(expected**2)
    # Assert
    np.testing.assert_allclose(series.values, expected, rtol=1e-15)
    series.inverse([0])
    series.subtract(10.0)
    series.absolute_value()
    expected[0] = 1.0 / expected[0]
    np.testing.assert_allclose(series.values, np.abs(expected - 10.0), rtol=1e-15)
    series.log_transform(np.e, [1])
    expected = np.abs(expected - 10.0)
    expected[1] = np.log(expected[1])
    np.testing.assert_allclose(series.values, expected, rtol=1e-14)


def test_standardize_matches_numpy_z_scores() -> None:
    # Arrange
    values = _gamma(50)
    series = _daily(values)
    # Act
    series.standardize()
    # Assert
    np.testing.assert_allclose(series.values, (values - values.mean()) / values.std(ddof=1), rtol=1e-12, atol=1e-14)


def test_missing_values_are_counted_replaced_and_interpolated_like_pandas() -> None:
    # Arrange
    values = _gamma(30)
    values[[4, 10, 11, 12, 20, 21, 22, 23, 24]] = np.nan
    series = _daily(values)
    # Act
    missing = series.number_of_missing_values()
    series.interpolate_missing_data(3)
    expected = pd.Series(values).interpolate(limit_area="inside")
    for start, stop in ((20, 25),):
        expected[start:stop] = np.nan
    # Assert
    assert missing == 9 and series.has_missing_values
    np.testing.assert_allclose(series.values, expected.to_numpy(), rtol=1e-14)
    series.replace_missing_data(-1.0)
    assert not series.has_missing_values
    np.testing.assert_array_equal(series.values[20:25], -1.0)


def test_replace_missing_data_at_indexes_only() -> None:
    # Arrange
    series = _daily([np.nan, 1.0, np.nan])
    # Act
    series.replace_missing_data(0.0, [2])
    # Assert
    assert np.isnan(series.values[0]) and series.values[2] == 0.0


def test_sorting_by_value_and_time() -> None:
    # Arrange
    series = _daily([3.0, 1.0, 2.0])
    # Act
    series.sort_by_value()
    by_value = series.values
    series.sort_by_value(descending=True)
    descending = series.values
    series.sort_by_time()
    # Assert
    np.testing.assert_array_equal(by_value, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(descending, [3.0, 2.0, 1.0])
    np.testing.assert_array_equal(series.values, [3.0, 1.0, 2.0])


# Smoothing and derived series -----------------------------------------------------------------


def test_cumulative_sum_and_differences_match_numpy_and_pandas() -> None:
    # Arrange
    values = _gamma(25)
    series = _daily(values)
    # Act
    cumulative = series.cumulative_sum()
    first = series.difference()
    second = series.difference(1, 2)
    seasonal = series.difference(lag=7)
    # Assert
    np.testing.assert_allclose(cumulative.values, np.cumsum(values), rtol=1e-14)
    np.testing.assert_allclose(first.values, np.diff(values), rtol=1e-13)
    np.testing.assert_allclose(second.values, np.diff(values, 2), rtol=1e-12, atol=1e-13)
    np.testing.assert_allclose(seasonal.values, pd.Series(values).diff(7).dropna().to_numpy(), rtol=1e-13)
    # By default each difference is dated at its earlier value.
    np.testing.assert_array_equal(first.dates, series.dates[:-1])
    np.testing.assert_array_equal(second.dates, series.dates[:-2])
    np.testing.assert_array_equal(seasonal.dates, series.dates[:-7])


def test_difference_dated_at_the_later_value_matches_pandas_diff() -> None:
    # Arrange
    values = _gamma(25)
    series = _daily(values)
    frame = pd.Series(values, index=pd.DatetimeIndex(series.dates))
    # Act
    first = series.difference(label_at_later_date=True)
    second = series.difference(1, 2, label_at_later_date=True)
    seasonal = series.difference(lag=7, label_at_later_date=True)
    # Assert: pandas Series.diff keeps the index of x[t]
    for result, expected in [(first, frame.diff()), (second, frame.diff().diff()), (seasonal, frame.diff(7))]:
        expected = expected.dropna()
        np.testing.assert_array_equal(result.dates, expected.index.to_numpy())
        np.testing.assert_allclose(result.values, expected.to_numpy(), rtol=1e-12, atol=1e-13)


@pytest.mark.parametrize("interval", ["one_day", "irregular"])
def test_difference_keeps_the_real_dates_of_a_gapped_or_irregular_series(interval: str) -> None:
    # Arrange
    dates = np.array(["2020-01-01", "2020-01-02", "2020-01-05", "2020-01-06", "2020-01-20"], dtype="datetime64[us]")
    values = np.array([1.0, 4.0, 9.0, 16.0, 25.0])
    series = TimeSeries(values, interval=interval, dates=dates)
    frame = pd.Series(values, index=pd.DatetimeIndex(dates))
    # Act
    earlier = series.difference()
    later = series.difference(label_at_later_date=True)
    # Assert: pandas diff is positional, so the labels are the series' own dates
    np.testing.assert_array_equal(earlier.dates, dates[:-1])
    np.testing.assert_array_equal(earlier.values, np.diff(values))
    np.testing.assert_array_equal(later.dates, frame.diff().dropna().index.to_numpy())
    np.testing.assert_array_equal(later.values, frame.diff().dropna().to_numpy())
    assert earlier.interval is TimeInterval(interval)


@pytest.mark.parametrize("min_valid", [None, 2])
def test_moving_average_and_sum_match_pandas_rolling(min_valid: int | None) -> None:
    # Arrange
    values = _gamma(60)
    values[[5, 30, 31]] = np.nan
    series = _daily(values)
    rolling = pd.Series(values).rolling(3, min_periods=3 if min_valid is None else min_valid)
    # Act
    average = series.moving_average(3, min_valid)
    total = series.moving_sum(3, min_valid)
    # Assert
    np.testing.assert_allclose(average.values, rolling.mean().to_numpy()[2:], rtol=1e-12)
    np.testing.assert_allclose(total.values, rolling.sum().to_numpy()[2:], rtol=1e-12)
    np.testing.assert_array_equal(average.dates, series.dates[2:])


def test_smoothed_series_dispatches_to_the_smoothing_members() -> None:
    # Arrange
    series = _daily(_gamma(20))
    # Act / Assert
    np.testing.assert_array_equal(series.smoothed_series("moving_average", 4).values, series.moving_average(4).values)
    np.testing.assert_array_equal(series.smoothed_series("moving_sum", 4).values, series.moving_sum(4).values)
    np.testing.assert_array_equal(series.smoothed_series("difference", 2).values, series.difference(2).values)
    np.testing.assert_array_equal(series.smoothed_series("none").values, series.values)


def test_date_shifts_clip_and_fill_match_pandas_offsets() -> None:
    # Arrange
    series = _daily([1.0, 2.0, 3.0, 4.0], start="2020-01-31")
    index = pd.DatetimeIndex(series.dates)
    # Act
    by_day = series.shift_dates_by_day(-2)
    by_month = series.shift_dates_by_month(1)
    by_year = series.shift_dates_by_year(1)
    clipped = series.clip("2020-02-01", "2020-02-02")
    gappy = TimeSeries([1.0, 3.0], interval="one_day", dates=["2020-01-01", "2020-01-03"])
    filled = gappy.fill_missing_dates("2020-01-01", "2020-01-04", 9.0)
    # Assert
    np.testing.assert_array_equal(by_day.dates, (index - pd.Timedelta(days=2)).to_numpy())
    np.testing.assert_array_equal(by_month.dates, np.array([d + pd.DateOffset(months=1) for d in index],
                                                           dtype="datetime64[us]"))
    np.testing.assert_array_equal(by_year.dates, np.array([d + pd.DateOffset(years=1) for d in index],
                                                          dtype="datetime64[us]"))
    np.testing.assert_array_equal(clipped.values, [2.0, 3.0])
    np.testing.assert_array_equal(filled.values, [1.0, 9.0, 3.0, 9.0])
    with pytest.raises(ValueError):
        series.clip("2019-01-01", "2020-02-01")


def test_shift_all_dates_restarts_the_regular_steps() -> None:
    # Arrange
    series = _daily([1.0, 2.0, 3.0])
    # Act
    series.shift_all_dates("2021-06-01")
    # Assert
    np.testing.assert_array_equal(series.dates, pd.date_range("2021-06-01", periods=3).to_numpy())


def test_interval_arithmetic_matches_pandas() -> None:
    # Act / Assert
    assert TimeSeries.add_time_interval("2020-01-31", "one_month") == np.datetime64("2020-02-29")
    assert TimeSeries.subtract_time_interval("2020-03-01", "one_day") == np.datetime64("2020-02-29")
    assert TimeSeries.add_time_interval("2020-01-01T23:45", "fifteen_minute") == np.datetime64("2020-01-02")
    assert TimeSeries.time_interval_in_hours("seven_day") == 168.0
    assert np.isnan(TimeSeries.time_interval_in_hours("one_year"))


# Interval conversion and aggregation ----------------------------------------------------------


def test_daily_to_hourly_conversion_interpolates_linearly_like_numpy() -> None:
    # Arrange
    values = _gamma(10)
    series = _daily(values)
    # Act
    hourly = series.convert_time_interval("one_hour")
    # Assert
    hours = np.arange(len(hourly), dtype=np.float64)
    np.testing.assert_allclose(hourly.values, np.interp(hours, 24.0 * np.arange(10), values), rtol=1e-13)
    np.testing.assert_array_equal(hourly.dates, pd.date_range("2020-01-01", periods=len(hourly), freq="h").to_numpy())


def test_hourly_to_daily_conversion_averages_or_sums_blocks() -> None:
    # Arrange
    values = _gamma(72)
    series = TimeSeries(values, interval="one_hour", start="2020-01-01")
    # Act
    daily_mean = series.convert_time_interval("one_day")
    daily_sum = series.convert_time_interval("one_day", average=False)
    # Assert
    # pandas resample("D") of 72 hourly values: three days and no empty block after the data.
    frame = pd.Series(values, index=pd.date_range("2020-01-01", periods=72, freq="h"))
    np.testing.assert_allclose(daily_mean.values, frame.resample("D").mean().to_numpy(), rtol=1e-13)
    np.testing.assert_allclose(daily_sum.values, frame.resample("D").sum().to_numpy(), rtol=1e-13)
    np.testing.assert_array_equal(daily_mean.dates, frame.resample("D").mean().index.to_numpy())
    partial = TimeSeries(values[:50], interval="one_hour", start="2020-01-01").convert_time_interval("one_day")
    np.testing.assert_allclose(partial.values, frame[:50].resample("D").mean().to_numpy(), rtol=1e-13)


@pytest.mark.parametrize(("block", "how"), [("average", "mean"), ("sum", "sum"), ("minimum", "min"),
                                            ("maximum", "max")])
def test_aggregate_by_interval_matches_pandas_resample(block: str, how: str) -> None:
    # Arrange
    values = _gamma(120)
    values[[3, 40, 41]] = np.nan
    series = _daily(values, start="2020-01-15")
    frame = pd.Series(values, index=pd.DatetimeIndex(series.dates))
    # Act
    monthly = series.aggregate_by_interval("one_month", block)
    # Assert
    expected = getattr(frame.resample("MS"), how)()
    np.testing.assert_allclose(monthly.values, expected.to_numpy(), rtol=1e-13)
    np.testing.assert_array_equal(monthly.dates, expected.index.to_numpy())


def test_aggregate_by_interval_min_valid_count_matches_pandas_min_count() -> None:
    # Arrange
    values = np.arange(48.0)
    values[0:20] = np.nan
    series = TimeSeries(values, interval="one_hour", start="2020-01-01")
    frame = pd.Series(values, index=pd.DatetimeIndex(series.dates))
    # Act
    daily = series.aggregate_by_interval("one_day", BlockFunction.SUM, min_valid_count=5)
    # Assert
    np.testing.assert_allclose(daily.values, frame.resample("D").sum(min_count=5).to_numpy())


# Block series and peaks -----------------------------------------------------------------------


def _three_years() -> TimeSeries:
    return _daily(_gamma(3 * 365 + 1, seed=11), start="2019-01-01")


def test_calendar_year_maxima_match_pandas_groupby() -> None:
    # Arrange
    series = _three_years()
    frame = series.to_pandas()
    # Act
    annual = series.calendar_year_series()
    annual_mean = series.calendar_year_series("average")
    # Assert
    grouped = frame.groupby(frame.index.year)
    np.testing.assert_array_equal(annual.values, grouped.max().to_numpy())
    np.testing.assert_array_equal(annual.dates, grouped.idxmax().to_numpy())
    np.testing.assert_allclose(annual_mean.values, grouped.mean().to_numpy(), rtol=1e-13)
    assert annual.interval is TimeInterval.IRREGULAR


def test_water_year_maxima_match_pandas_groupby_on_the_ending_year() -> None:
    # Arrange
    series = _three_years()
    frame = series.to_pandas()
    # Act
    water = series.custom_year_series(10)
    # Assert
    water_year = frame.index.year + (frame.index.month >= 10)
    np.testing.assert_array_equal(water.values, frame.groupby(water_year).max().to_numpy())


def test_seasonal_window_and_monthly_and_quarterly_maxima_match_pandas() -> None:
    # Arrange
    series = _three_years()
    frame = series.to_pandas()
    # Act
    summer = series.custom_year_series(6, 8, "maximum")
    monthly = series.monthly_series()
    quarterly = series.quarterly_series("minimum")
    # Assert
    in_summer = frame[(frame.index.month >= 6) & (frame.index.month <= 8)]
    np.testing.assert_array_equal(summer.values, in_summer.groupby(in_summer.index.year).max().to_numpy())
    np.testing.assert_array_equal(monthly.values, frame.resample("MS").max().to_numpy())
    np.testing.assert_array_equal(quarterly.values, frame.resample("QS").min().to_numpy())


def test_block_series_of_a_smoothed_series_equal_blocks_of_the_moving_average() -> None:
    # Arrange
    series = _three_years()
    # Act
    smoothed = series.calendar_year_series("maximum", "moving_average", 3)
    # Assert
    expected = series.moving_average(3).calendar_year_series("maximum")
    np.testing.assert_array_equal(smoothed.values, expected.values)


def test_peaks_over_threshold_returns_one_peak_per_separated_cluster() -> None:
    # Arrange
    values = np.array([0, 5, 7, 6, 0, 0, 0, 9, 0, 8, 0, 0, 0, 0, 3.0])
    series = _daily(values)
    # Act
    peaks = series.peaks_over_threshold(4.0, min_steps_between_events=2)
    # Assert
    np.testing.assert_array_equal(peaks.values, [7.0, 9.0])
    np.testing.assert_array_equal(peaks.dates, series.dates[[2, 7]])


# Resampling and decomposition -----------------------------------------------------------------


def test_block_bootstrap_is_reproducible_and_made_of_original_blocks() -> None:
    # Arrange
    values = _gamma(60)
    series = _daily(values)
    # Act
    first = series.resample_with_block_bootstrap(100, 10, seed=42)
    second = series.resample_with_block_bootstrap(100, 10, seed=42)
    # Assert
    np.testing.assert_array_equal(first.values, second.values)
    windows = {tuple(values[i:i + 10]) for i in range(51)}
    for start in range(0, 100, 10):
        assert tuple(first.values[start:start + 10]) in windows
    np.testing.assert_array_equal(first.dates, pd.date_range("2020-01-01", periods=100).to_numpy())


def test_knn_resampling_is_reproducible_and_draws_observed_values() -> None:
    # Arrange
    values = _gamma(60)
    series = _daily(values)
    # Act
    first = series.resample_with_knn(200, 5, seed=41)
    second = series.resample_with_knn(200, 5, seed=41)
    # Assert
    np.testing.assert_array_equal(first.values, second.values)
    assert np.isin(first.values, values).all() and len(first) == 200


def test_seasonal_decomposition_adds_back_to_the_series() -> None:
    # Arrange
    t = np.arange(84, dtype=np.float64)
    values = 40.0 + 0.05 * t + 6.0 * np.sin(2.0 * np.pi * t / 7.0) + np.random.default_rng(3).normal(0, 1, 84)
    series = _daily(values)
    # Act
    trend, seasonal, residual = series.seasonal_decompose(7)
    # Assert
    np.testing.assert_array_equal(trend.values, series.moving_average(7).values)
    np.testing.assert_allclose(residual.values, values[6:] - trend.values - seasonal[6:], rtol=1e-12, atol=1e-12)
    assert seasonal.shape == (84,)
    assert abs(np.corrcoef(seasonal, 6.0 * np.sin(2.0 * np.pi * t / 7.0))[0, 1]) > 0.9


# Statistics -----------------------------------------------------------------------------------


def test_summary_values_match_numpy_with_missing_values_skipped() -> None:
    # Arrange
    values = _gamma(101)
    values[[7, 70]] = np.nan
    series = _daily(values)
    # Act
    summary = series.summary_statistics()
    # Assert
    assert series.min_value() == np.nanmin(values) and series.max_value() == np.nanmax(values)
    assert series.mean_value() == pytest.approx(np.nanmean(values), rel=1e-14)
    assert series.standard_deviation() == pytest.approx(np.nanstd(values, ddof=1), rel=1e-13)
    assert summary["Record Length"] == 101 and summary["Missing Values"] == 2
    assert summary["Mean"] == pytest.approx(np.nanmean(values), rel=1e-14)
    assert summary["50%"] == pytest.approx(np.nanpercentile(values, 50), rel=1e-14)


def test_percentiles_are_type_7_like_numpy() -> None:
    # Arrange
    values = _gamma(57)
    values[3] = np.nan
    series = _daily(values)
    k = np.array([0.01, 0.1, 0.5, 0.9, 0.99])
    # Act / Assert
    np.testing.assert_allclose(series.percentiles(k), np.nanpercentile(values, 100 * k), rtol=1e-13)
    np.testing.assert_allclose(series.summary_percentiles(),
                               np.nanpercentile(values, [5, 25, 50, 75, 95]), rtol=1e-13)


def test_duration_curve_uses_weibull_plotting_positions() -> None:
    # Arrange
    values = _gamma(20)
    series = _daily(values)
    # Act
    curve = series.duration()
    # Assert
    np.testing.assert_allclose(curve[:, 0], 100.0 * np.arange(1, 21) / 21.0, rtol=1e-14)
    np.testing.assert_array_equal(curve[:, 1], np.sort(values)[::-1])


def test_monthly_statistics_match_pandas_groupby_month() -> None:
    # Arrange
    series = _three_years()
    frame = series.to_pandas()
    by_month = frame.groupby(frame.index.month)
    # Act
    frequency = series.monthly_frequency()
    summary = series.monthly_summary_statistics()
    percentiles = series.monthly_percentiles([0.1, 0.5])
    # Assert
    np.testing.assert_array_equal(frequency, by_month.size().to_numpy())
    np.testing.assert_array_equal(summary[:, 0], by_month.min().to_numpy())
    np.testing.assert_array_equal(summary[:, 6], by_month.max().to_numpy())
    np.testing.assert_allclose(summary[:, 7], by_month.mean().to_numpy(), rtol=1e-13)
    np.testing.assert_allclose(summary[:, 3], by_month.median().to_numpy(), rtol=1e-14)
    np.testing.assert_allclose(percentiles[:, 0], by_month.quantile(0.1).to_numpy(), rtol=1e-13)


def test_summary_hypothesis_tests_report_seven_p_values() -> None:
    # Arrange
    series = _daily(_gamma(80))
    # Act
    tests = series.summary_hypothesis_test()
    # Assert
    assert len(tests) == 7
    assert all(0.0 <= p <= 1.0 for p in tests.values())


def test_invalid_arguments_raise_value_error() -> None:
    # Arrange
    series = _daily([1.0, 2.0, 3.0])
    # Act / Assert
    with pytest.raises(ValueError, match="period"):
        series.moving_average(5)
    with pytest.raises(ValueError):
        series.custom_year_series(13)
    with pytest.raises(ValueError):
        TimeSeries([1.0], interval="irregular", start="2020-01-01")
    with pytest.raises(ValueError, match="one value per date"):
        TimeSeries([1.0, 2.0], dates=["2020-01-01"])


# Speed ----------------------------------------------------------------------------------------


def test_values_and_dates_cross_in_bulk_faster_than_per_ordinate() -> None:
    # Arrange
    series = _daily(_gamma(20000))
    net = series.net
    # Act
    start = time.perf_counter()
    values, dates = series.values, series.dates
    bulk = time.perf_counter() - start
    start = time.perf_counter()
    slow = np.array([float(o.Value) for o in net])
    per_element = time.perf_counter() - start
    # Assert
    np.testing.assert_array_equal(values, slow)
    assert len(dates) == 20000
    assert bulk * 10.0 < per_element, (bulk, per_element)


def test_explicit_dates_build_faster_than_appending_each_ordinate() -> None:
    # Arrange
    n = 5000
    dates = np.datetime64("2000-01-01") + np.sort(np.random.default_rng(1).choice(20000, n, replace=False))
    values = _gamma(n)
    # Act
    start = time.perf_counter()
    bulk = TimeSeries(values, interval="irregular", dates=dates)
    bulk_time = time.perf_counter() - start
    start = time.perf_counter()
    slow = TimeSeries(interval="irregular")
    for d, v in zip(dates, values, strict=True):
        slow.append(d, v)
    append_time = time.perf_counter() - start
    # Assert
    np.testing.assert_array_equal(bulk.ticks, slow.ticks)
    np.testing.assert_array_equal(bulk.values, slow.values)
    assert bulk_time < append_time, (bulk_time, append_time)
