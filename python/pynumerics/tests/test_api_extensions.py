"""Wrappers over the generalized members: the plotting-position alpha table and its weighted typed overload."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import optimize
from scipy import stats as sps

from pynumerics import distributions, mathematics, statistics


def _bits(values: object) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64)


def test_plotting_position_alpha_table_and_the_weighted_named_formulas() -> None:
    # Arrange
    rng = np.random.default_rng(23)
    weights = rng.random(80) * 3.0
    table = {"weibull": 0.0, "median": 0.3175, "blom": 0.375, "cunnane": 0.4, "gringorten": 0.44, "hazen": 0.5}
    # Act / Assert
    for name, alpha in table.items():
        assert statistics.plotting_position_alpha(name) == alpha
        typed = statistics.weighted_plotting_positions(weights, name)
        assert _bits(typed).tolist() == _bits(statistics.weighted_plotting_positions(weights, alpha)).tolist()
        assert _bits(statistics.weighted_plotting_positions(np.ones(25), name)).tolist() == _bits(
            statistics.plotting_positions(25, name)).tolist()
    # Kaplan-Meier: the cumulative weight share, by NumPy (exact for these weights)
    km = statistics.weighted_plotting_positions([1.0, 3.0, 0.0, 4.0], "kaplan_meier")
    assert km.tolist() == (np.cumsum([1.0, 3.0, 0.0, 4.0]) / 8.0).tolist()
    assert _bits(statistics.weighted_plotting_positions(np.ones(30), "kaplan_meier")).tolist() == _bits(
        statistics.plotting_positions(30, "kaplan_meier")).tolist()
    with pytest.raises(ValueError, match="not a member"):
        statistics.weighted_plotting_positions([1.0], "nope")


def test_wilson_score_interval_returns_a_tuple_of_floats() -> None:
    # Arrange / Act: Numerics returns a (Lower, Upper) value tuple; the wrapper returns a tuple of floats
    bounds = statistics.wilson_score_interval(3.0, 17.0, 0.9)
    # Assert: the wrapper's shape is unchanged, and the values agree with the closed form by SciPy
    assert isinstance(bounds, tuple) and len(bounds) == 2 and all(type(b) is float for b in bounds)
    z = sps.norm.ppf(0.95)
    p, n = 3.0 / 17.0, 17.0
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    np.testing.assert_allclose(bounds, (center - half, center + half), rtol=1e-12)


def test_fixed_location_estimate_equals_the_known_location_fit_and_scipy() -> None:
    # Arrange
    excesses = np.array([1.2, 0.4, 3.9, 2.2, 0.15, 5.6, 1.05, 0.8, 7.3, 2.9, 0.33, 1.7, 4.4, 0.62, 2.5, 11.2, 0.09,
                         1.33, 3.1, 0.95])
    fit = distributions.create_distribution("GeneralizedPareto")
    # Act
    fit.estimate(excesses, "maximum_likelihood", [0.0, np.nan, np.nan])
    known = distributions.GeneralizedPareto.fit_known_location(excesses, 0.0)
    # Assert: bitwise the known-location fit, and close to SciPy's maximum (Nelder-Mead tolerance)
    assert _bits(fit.parameters).tolist() == _bits(known.parameters).tolist()
    best = optimize.minimize(lambda v: -sps.genpareto.logpdf(excesses, -v[1], 0.0, v[0]).sum(), [2.0, 0.0],
                             method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20000})
    np.testing.assert_allclose(fit.parameters[1:], best.x, rtol=1e-3, atol=1e-3)
    unconstrained = distributions.create_distribution("GeneralizedPareto")
    unconstrained.estimate(excesses, "method_of_linear_moments", [np.nan, np.nan, np.nan])
    reference = distributions.create_distribution("GeneralizedPareto")
    reference.estimate(excesses, "method_of_linear_moments")
    assert _bits(unconstrained.parameters).tolist() == _bits(reference.parameters).tolist()
    with pytest.raises(TypeError, match="fixed parameters"):
        distributions.create_distribution("Normal").estimate(excesses, "maximum_likelihood", [0.0, np.nan])


def test_column_means_equal_the_per_column_mean_and_numpy() -> None:
    # Arrange
    rng = np.random.default_rng(29)
    data = 1e3 * np.exp(rng.standard_normal((500, 40)))
    data[rng.random(data.shape) < 0.03] = np.nan
    data[7, 3] = -np.inf
    data[:, 11] = np.nan
    # Act
    means, excluded = statistics.column_means(data)
    # Assert: bitwise the sequential mean of the finite entries; NumPy (pairwise sums) to 1e-14
    for column in range(data.shape[1]):
        finite = data[np.isfinite(data[:, column]), column]
        assert excluded[column] == data.shape[0] - finite.size
        if finite.size == 0:
            assert np.isnan(means[column])
            continue
        assert _bits(means[column]).tolist() == _bits(statistics.mean(finite)).tolist()
        np.testing.assert_allclose(means[column], finite.mean(), rtol=1e-14)


def test_log1p_agrees_with_python_math_log1p_within_one_unit_in_the_last_place() -> None:
    # Arrange
    rng = np.random.default_rng(37)
    x = np.concatenate([-rng.random(20000) * 0.999999, 10.0 ** rng.uniform(-16, 3, 20000),
                        -(10.0 ** rng.uniform(-16, -0.0001, 20000)), 10.0 ** rng.uniform(3, 300, 500),
                        np.linspace(1e-4, 1e-2, 2001), -np.linspace(1e-4, 1e-2, 2001)])
    # Act
    values = mathematics.log1p(x)
    reference = np.array([math.log1p(v) for v in x.tolist()])
    # Assert: at most one unit in the last place (0 on Windows, where both use the same C runtime log)
    ulps = np.abs(values.view(np.int64) - reference.view(np.int64))
    assert int(ulps.max()) <= 1
    assert math.isinf(mathematics.log1p(-1.0)) and mathematics.log1p(-1.0) < 0
    assert math.isnan(mathematics.log1p(-1.5)) and math.isnan(mathematics.log1p(float("nan")))
    assert mathematics.log1p(float("inf")) == float("inf")
