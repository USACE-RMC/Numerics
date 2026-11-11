"""Interpolators and searches from ``Numerics.Data``.

``scipy.interpolate.CubicSpline`` (natural), ``RegularGridInterpolator``, ``BarycentricInterpolator``,
``numpy.polyfit`` and ``numpy.searchsorted`` are the oracles.
"""

from __future__ import annotations

import copy
import time

import numpy as np
import pytest
from scipy import interpolate, stats

from pynumerics import data

KNOTS = np.array([0.0, 0.7, 1.5, 2.0, 3.2, 4.0, 5.5])
VALUES = np.sin(KNOTS) + 0.1 * KNOTS


def test_cubic_spline_matches_scipy_natural_spline_inside_and_beyond_the_knots() -> None:
    # Arrange
    spline = data.CubicSpline(KNOTS, VALUES)
    x = np.linspace(-0.5, 6.0, 101)
    # Act
    y = spline.interpolate(x)
    # Assert
    np.testing.assert_allclose(y, interpolate.CubicSpline(KNOTS, VALUES, bc_type="natural")(x), rtol=1e-12,
                               atol=1e-13)
    assert spline(1.0) == pytest.approx(float(interpolate.CubicSpline(KNOTS, VALUES, bc_type="natural")(1.0)),
                                        rel=1e-12)


def test_cubic_spline_on_descending_knots_matches_scipy() -> None:
    # Arrange
    spline = data.CubicSpline(KNOTS[::-1], VALUES[::-1], "descending")
    x = np.linspace(0.1, 5.4, 37)
    # Act
    spline.set_second_derivatives()
    y = spline(x)
    # Assert
    np.testing.assert_allclose(y, interpolate.CubicSpline(KNOTS, VALUES, bc_type="natural")(x), rtol=1e-12)
    assert spline.sort_order.value == "Descending"


def test_polynomial_reproduces_a_cubic_like_numpy_polyfit() -> None:
    # Arrange
    coefficients = np.polyfit(KNOTS, 0.5 * KNOTS**3 - KNOTS + 2.0, 3)
    y = np.polyval(coefficients, KNOTS)
    polynomial = data.Polynomial(3, KNOTS, y)
    x = np.linspace(0.0, 5.5, 23)
    # Act
    values = polynomial(x)
    # Assert
    np.testing.assert_allclose(values, np.polyval(coefficients, x), rtol=1e-11, atol=1e-12)
    assert polynomial.order == 3 and np.isfinite(polynomial.error)


def test_polynomial_uses_the_knots_from_the_bracketing_index() -> None:
    # Arrange
    polynomial = data.Polynomial(2, KNOTS, VALUES)
    x = 1.8
    # Act
    value = polynomial(x)
    # Assert
    start = int(np.searchsorted(KNOTS, x) - 1)
    window = slice(start, start + 3)
    expected = interpolate.BarycentricInterpolator(KNOTS[window], VALUES[window])(x)
    assert value == pytest.approx(float(expected), rel=1e-13)
    polynomial.order = 1
    assert polynomial(x) == pytest.approx(float(np.interp(x, KNOTS, VALUES)), rel=1e-13)
    with pytest.raises(ValueError, match="order"):
        data.Polynomial(7, KNOTS, VALUES)


def test_bilinear_matches_regular_grid_interpolator_and_holds_the_edges() -> None:
    # Arrange
    x1, x2 = np.array([0.0, 1.0, 2.5, 4.0]), np.array([10.0, 20.0, 35.0])
    y = np.add.outer(x1**2, np.sqrt(x2))
    grid = data.Bilinear(x1, x2, y)
    points = np.random.default_rng(2).uniform([-1.0, 5.0], [5.0, 40.0], (200, 2))
    # Act
    values = grid(points[:, 0], points[:, 1])
    # Assert
    clipped = np.column_stack([np.clip(points[:, 0], 0, 4), np.clip(points[:, 1], 10, 35)])
    np.testing.assert_allclose(values, interpolate.RegularGridInterpolator((x1, x2), y)(clipped), rtol=1e-13)
    np.testing.assert_array_equal(grid.y_values, y)
    np.testing.assert_array_equal(grid.x1_values, x1)
    np.testing.assert_array_equal(grid.x2_values, x2)
    scalar = grid.interpolate(2.0, 20.0)
    assert isinstance(scalar, float)
    assert scalar == pytest.approx(float(interpolate.RegularGridInterpolator((x1, x2), y)([2.0, 20.0])[0]), rel=1e-14)


def test_bilinear_in_log_space_and_its_settings() -> None:
    # Arrange
    x1, x2 = np.array([1.0, 10.0, 100.0]), np.array([1.0, 10.0])
    y = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    grid = data.Bilinear(x1, x2, y, x1_transform="logarithmic")
    # Act
    value = grid(np.sqrt(10.0), 1.0)
    grid.x2_transform = "logarithmic"
    grid.y_transform = "none"
    grid.use_smart_search = False
    # Assert
    assert value == pytest.approx(2.0, rel=1e-14)
    assert grid(np.sqrt(10.0), np.sqrt(10.0)) == pytest.approx(2.5, rel=1e-14)
    assert grid.x1_transform.value == "Logarithmic" and grid.x2_transform.value == "Logarithmic"
    assert grid.y_transform.value == "None" and not grid.use_smart_search and grid.sort_order.value == "Ascending"
    with pytest.raises(ValueError):
        data.Bilinear(x1, x2, y.T)


def test_linear_transforms_match_numpy_interp_in_the_transformed_spaces() -> None:
    # Arrange
    p = np.array([0.01, 0.1, 0.5, 0.9, 0.99])
    q = np.array([1.0, 10.0, 50.0, 200.0, 1000.0])
    line = data.Linear(p, q, x_transform="normal_z", y_transform="logarithmic")
    x = np.array([0.05, 0.3, 0.95])
    # Act
    y = line(x)
    # Assert
    expected = 10 ** np.interp(stats.norm.ppf(x), stats.norm.ppf(p), np.log10(q))
    np.testing.assert_allclose(y, expected, rtol=1e-12)
    assert line.x_transform.value == "NormalZ" and line.y_transform.value == "Logarithmic"
    line.y_transform = "none"
    np.testing.assert_allclose(line(x), np.interp(stats.norm.ppf(x), stats.norm.ppf(p), q), rtol=1e-12)


def test_linear_clone_and_copy_share_knots_but_not_the_cursor() -> None:
    # Arrange
    line = data.Linear(KNOTS, VALUES)
    # Act
    clone = line.shallow_clone()
    deep = copy.deepcopy(line)
    clone.search_start = 4
    # Assert
    assert clone.net is not line.net and line.search_start == 0
    np.testing.assert_array_equal(clone(KNOTS), line(KNOTS))
    np.testing.assert_array_equal(deep(KNOTS), line(KNOTS))
    np.testing.assert_array_equal(line.x_values, KNOTS)
    np.testing.assert_array_equal(line.y_values, VALUES)
    assert len(line) == line.count == len(KNOTS)


def test_interpolater_searches_match_numpy_searchsorted() -> None:
    # Arrange
    line = data.Linear(KNOTS, VALUES)
    x = np.linspace(0.05, 5.45, 40)
    expected = np.searchsorted(KNOTS, x) - 1
    # Act
    bisection = [line.bisection_search(v) for v in x]
    hunt = [line.hunt_search(v) for v in x]
    sequential = [line.sequential_search(v) for v in x]
    smart = [line.search(v) for v in x]
    line.use_smart_search = False
    # Assert
    np.testing.assert_array_equal(bisection, expected)
    np.testing.assert_array_equal(hunt, expected)
    np.testing.assert_array_equal(sequential, expected)
    np.testing.assert_array_equal(smart, expected)
    assert not line.use_smart_search
    assert line.base_interpolate(1.0, 1) == pytest.approx(float(np.interp(1.0, KNOTS, VALUES)), rel=1e-14)


@pytest.mark.parametrize("method", list(data.SearchMethod))
def test_search_functions_match_numpy_searchsorted(method: data.SearchMethod) -> None:
    # Arrange
    x = np.linspace(0.05, 5.45, 25)
    pairs = data.OrderedPairedData(KNOTS, np.arange(7.0))
    # Act
    on_values = [data.search(v, KNOTS, method) for v in x]
    on_pairs = [data.search(v, pairs, method.value) for v in x]
    descending = [data.search(v, KNOTS[::-1], method, order="descending") for v in x]
    # Assert
    expected = np.searchsorted(KNOTS, x) - 1
    np.testing.assert_array_equal(on_values, expected)
    np.testing.assert_array_equal(on_pairs, expected)
    np.testing.assert_array_equal(descending, len(KNOTS) - 2 - expected)


def test_array_interpolation_gives_the_scalar_call_bits() -> None:
    # Arrange
    spline = data.CubicSpline(KNOTS, VALUES)
    x = np.random.default_rng(4).uniform(0.0, 5.5, 20000)
    # Act
    bulk = spline(x)
    scalar = np.array([spline(v) for v in x])
    # Assert
    np.testing.assert_array_equal(bulk.view(np.uint64), scalar.view(np.uint64))


@pytest.mark.benchmark
def test_array_interpolation_crosses_once_faster_than_scalar_calls() -> None:
    # Arrange
    spline = data.CubicSpline(KNOTS, VALUES)
    x = np.random.default_rng(4).uniform(0.0, 5.5, 20000)
    # Act
    start = time.perf_counter()
    spline(x)
    bulk_time = time.perf_counter() - start
    start = time.perf_counter()
    _ = [spline(v) for v in x]
    scalar_time = time.perf_counter() - start
    # Assert
    assert bulk_time * 10.0 < scalar_time, (bulk_time, scalar_time)


# List overloads -------------------------------------------------------------------------------


def _bilinear_grids() -> tuple[data.Bilinear, data.Bilinear, np.ndarray]:
    """Return two identical logarithmic-x1 bilinear grids and 3,000 query points."""
    x1, x2 = np.array([1.0, 2.0, 4.0, 8.0]), np.array([10.0, 20.0, 35.0])
    y = np.add.outer(np.log(x1), np.sqrt(x2))
    points = np.random.default_rng(4).uniform([0.5, 5.0], [9.0, 40.0], (3000, 2))
    grid = data.Bilinear(x1, x2, y, x1_transform="logarithmic")
    reference = data.Bilinear(x1, x2, y, x1_transform="logarithmic")
    return grid, reference, points


def test_bilinear_arrays_give_the_scalar_bits() -> None:
    # Arrange
    grid, reference, points = _bilinear_grids()
    # Act
    values = grid(points[:, 0], points[:, 1])
    expected = np.array([reference.net.Interpolate(float(a), float(b)) for a, b in points])
    # Assert
    np.testing.assert_array_equal(values.view(np.uint64), expected.view(np.uint64))


@pytest.mark.benchmark
def test_bilinear_arrays_cross_once_faster_than_scalar_calls() -> None:
    # Arrange
    grid, reference, points = _bilinear_grids()
    # Act
    begin = time.perf_counter()
    grid(points[:, 0], points[:, 1])
    list_seconds = time.perf_counter() - begin
    begin = time.perf_counter()
    _ = [reference.net.Interpolate(float(a), float(b)) for a, b in points]
    loop_seconds = time.perf_counter() - begin
    # Assert (one crossing takes a fraction of the time of the scalar calls)
    assert list_seconds < loop_seconds / 3


def test_linear_extrapolate_arrays_cross_once_with_the_scalar_bits() -> None:
    # Arrange: R Hmisc approxExtrap(x, y, xout) for the straight end segments
    x, y = np.array([50.0, 100.0, 150.0, 200.0, 250.0]), np.array([100.0, 200.0, 300.0, 400.0, 500.0])
    out = np.array([300.0, 310.0, 350.0, 400.0, 440.0, 0.0, -10.0, 75.0])
    line = data.Linear(x, y)
    logged = data.Linear(x, y, x_transform="logarithmic")
    reference = data.Linear(x, y, x_transform="logarithmic")
    points = np.random.default_rng(3).uniform(10.0, 400.0, 500)
    # Act
    values = line.extrapolate(out)
    log_values = logged.extrapolate(points)
    # Assert
    np.testing.assert_allclose(values, [600.0, 620.0, 700.0, 800.0, 880.0, 0.0, -20.0, 150.0], rtol=1e-14, atol=1e-12)
    expected = np.array([reference.net.Extrapolate(float(v)) for v in points])
    np.testing.assert_array_equal(log_values.view(np.uint64), expected.view(np.uint64))
    assert isinstance(line.extrapolate(300.0), float)
