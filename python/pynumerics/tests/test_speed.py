"""Bound overloads, cached enum members and the list forms that remove per-element crossings.

Every faster path must return the bits of the path it replaces (the scalar member element by
element, or the unbound method group). Where SciPy or NumPy compute the same quantity, the values
are also checked against them at a tolerance that covers the different operation order.
"""

from __future__ import annotations

import importlib
import json
import math

import numpy as np
import pytest
from scipy import special, stats

from pynumerics import _util, copulas, data, distributions, mathematics, statistics
from pynumerics import special as pn_special
from pynumerics.marshal import to_net_doubles, to_numpy_doubles
from pynumerics.runtime import initialize


def _bits(values: object) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64)


def _has_overload(namespace: str, type_name: str, member: str, *parameters: str) -> bool:
    initialize()
    import clr  # noqa: PLC0415
    import System  # noqa: PLC0415

    owner = clr.GetClrType(getattr(importlib.import_module(namespace), type_name))
    types = System.Array[System.Type]([System.Type.GetType(name) for name in parameters])
    return owner.GetMethod(member, types) is not None


_LIST = "System.Collections.Generic.IList`1[System.Double]"


def test_net_enum_resolves_names_once_and_still_refuses_unknown_names() -> None:
    # Arrange
    stats_ns = statistics._stats()
    kind = stats_ns.PlottingPositions.PlottingPostionType
    # Act
    first = _util.net_enum(kind, "kaplan_meier")
    second = _util.net_enum(kind, "KaplanMeier")
    # Assert
    assert first == kind.KaplanMeier and second == kind.KaplanMeier
    assert (kind, "kaplanmeier") in _util._ENUM_MEMBERS
    with pytest.raises(ValueError, match="is not a member of PlottingPostionType"):
        _util.net_enum(kind, "not_a_formula")


def test_bound_statistics_overloads_return_the_bits_of_the_method_group() -> None:
    # Arrange
    rng = np.random.default_rng(11)
    sample = rng.standard_normal(301)
    sample[:4] = [0.0, -0.0, 0.0, -0.0]
    weights = rng.random(301)
    k = np.linspace(0.0, 1.0, 41)
    group = statistics._stats().Statistics
    net_sample, net_k, net_w = to_net_doubles(sample), to_net_doubles(k), to_net_doubles(weights)
    # Act / Assert
    assert _bits(statistics.percentile(sample, k)).tolist() == _bits(
        to_numpy_doubles(group.Percentile(net_sample, net_k, False))).tolist()
    assert _bits(statistics.percentile(sample, k, weights=weights)).tolist() == _bits(
        to_numpy_doubles(group.Percentile(net_sample, net_k, net_w, False))).tolist()
    assert statistics.percentile(sample, 0.37) == float(group.Percentile(net_sample, 0.37, False))
    assert statistics.mean(sample) == float(group.Mean(net_sample))
    assert statistics.mean(sample, weights) == float(group.Mean(net_sample, net_w))
    assert statistics.variance(sample) == float(group.Variance(net_sample))
    assert statistics.standard_deviation(sample, weights, "reliability") == float(
        group.StandardDeviation(net_sample, net_w, statistics._stats().WeightType.Reliability))


@pytest.mark.parametrize(("family", "parameters", "frozen", "points"), [
    ("GeneralizedExtremeValue", [100.0, 20.0, 0.1], stats.genextreme(0.1, 100.0, 20.0), [50.0, 100.0, 150.0, 250.0]),
    ("GeneralizedPareto", [0.0, 10.0, -0.1], stats.genpareto(0.1, 0.0, 10.0), [0.5, 5.0, 20.0, 60.0]),
    ("Normal", [100.0, 15.0], stats.norm(100.0, 15.0), [40.0, 100.0, 130.0, 170.0]),
    ("Gumbel", [100.0, 20.0], stats.gumbel_r(100.0, 20.0), [60.0, 100.0, 150.0, 250.0]),
])
def test_ccdf_and_log_members_on_arrays_equal_the_scalar_members_and_scipy(
        family: str, parameters: list[float], frozen: object, points: list[float]) -> None:
    # Arrange
    distribution = distributions.create_distribution(family, parameters)
    x = np.concatenate([points, np.linspace(-50.0, 400.0, 91), [-0.0, 0.0, np.inf, -np.inf, np.nan]])
    # Act / Assert: the array path is the scalar member, bit for bit
    for name, member in (("ccdf", "CCDF"), ("log_pdf", "LogPDF"), ("log_cdf", "LogCDF"), ("log_ccdf", "LogCCDF")):
        array = getattr(distribution, name)(x)
        scalar = np.array([float(getattr(distribution.net, member)(float(v))) for v in x])
        assert _bits(array).tolist() == _bits(scalar).tolist(), name
    # SciPy oracle at interior points (different operation order: 1e-12 relative)
    p = np.asarray(points)
    np.testing.assert_allclose(distribution.ccdf(p), frozen.sf(p), rtol=1e-12)
    np.testing.assert_allclose(distribution.log_pdf(p), frozen.logpdf(p), rtol=1e-12)
    np.testing.assert_allclose(distribution.log_cdf(p), frozen.logcdf(p), rtol=1e-12)
    np.testing.assert_allclose(distribution.log_ccdf(p), frozen.logsf(p), rtol=1e-12)


def test_ccdf_arrays_use_the_list_overload() -> None:
    # Arrange
    distribution = distributions.create_distribution("Normal", [0.0, 1.0])
    # Act
    distribution.ccdf(np.zeros(3))
    # Assert
    assert distribution._optional_lists["CCDF"] is not None


def test_expm1_and_log1p_arrays_equal_the_scalar_members_and_numpy() -> None:
    # Arrange
    rng = np.random.default_rng(3)
    x = np.concatenate([[-800.0, -1.0, -1.0001e-4, -1e-4, -1e-12, -0.0, 0.0, 1e-12, 1e-4, 0.7, 709.0, 710.0, np.inf,
                         np.nan], rng.standard_normal(200) * 1e-3])
    tools = importlib.import_module("Numerics").Tools
    # Act
    e = mathematics.expm1(x)
    log_x = x[~(x <= -1.0)]
    lg = mathematics.log1p(log_x)
    # Assert
    assert _bits(e).tolist() == _bits([float(tools.Expm1(float(v))) for v in x]).tolist()
    assert _bits(lg).tolist() == _bits([float(tools.Log1p(float(v))) for v in log_x]).tolist()
    assert mathematics.expm1(1e-7) == float(tools.Expm1(1e-7))
    finite = np.isfinite(x) & (np.abs(x) < 700)
    np.testing.assert_allclose(e[finite], np.expm1(x[finite]), rtol=1e-15)
    # Below |x| = 1e-4 Tools.Log1p agrees with the C library to 1e-15; test_api_extensions holds it within one unit
    # in the last place at every argument.
    compensated = np.isfinite(log_x) & (np.abs(log_x) <= 1e-4)
    np.testing.assert_allclose(lg[compensated], np.log1p(log_x[compensated]), rtol=1e-15)
    np.testing.assert_allclose(lg[np.isfinite(log_x)], np.log1p(log_x[np.isfinite(log_x)]), rtol=1e-11)


def test_beta_incomplete_bound_once_matches_scipy() -> None:
    # Arrange
    x = np.linspace(0.0, 1.0, 21)
    # Act
    values = pn_special.beta_incomplete(2.0, 3.0, x)
    inverse = pn_special.beta_incomplete_inverse(2.0, 3.0, values)
    # Assert
    np.testing.assert_allclose(values, special.betainc(2.0, 3.0, x), rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(inverse, x, atol=1e-10)


@pytest.mark.parametrize("copula", [
    lambda: copulas.NormalCopula(0.6), lambda: copulas.ClaytonCopula(2.0), lambda: copulas.FrankCopula(5.0),
    lambda: copulas.GumbelCopula(2.0), lambda: copulas.JoeCopula(2.0), lambda: copulas.AMHCopula(0.5),
])
def test_copula_arrays_equal_the_scalar_pairs(copula: object) -> None:
    # Arrange
    c = copula()
    rng = np.random.default_rng(5)
    u = np.concatenate([[0.0, 1.0, 0.5, 1e-12], rng.random(100)])
    v = np.concatenate([[0.5, 0.5, 0.0, 0.3], rng.random(100)])
    # Act / Assert
    for name, member in (("pdf", "PDF"), ("cdf", "CDF"), ("conditional_cdf", "ConditionalCDF"),
                         ("inverse_conditional_cdf", "InverseConditionalCDF")):
        array = getattr(c, name)(u, v)
        scalar = [float(getattr(c.net, member)(float(a), float(b))) for a, b in zip(u, v, strict=True)]
        assert _bits(array).tolist() == _bits(scalar).tolist(), name
    second = [float(c.inverse_cdf(a, b)[1]) for a, b in zip(u, v, strict=True)]
    assert _bits(c.inverse_conditional_cdf(u, v)).tolist() == _bits(second).tolist()


def test_normal_copula_conditional_inverse_matches_the_closed_form() -> None:
    # Arrange
    rng = np.random.default_rng(9)
    u, t, rho = rng.random(500), rng.random(500), 0.6
    # Act
    values = copulas.NormalCopula(rho).inverse_conditional_cdf(u, t)
    # Assert: Phi(rho z_u + sqrt(1 - rho^2) z_t), by SciPy
    expected = special.ndtr(rho * special.ndtri(u) + math.sqrt(1 - rho * rho) * special.ndtri(t))
    np.testing.assert_allclose(values, expected, rtol=1e-12, atol=1e-15)


def test_column_plotting_position_percentiles_equal_the_per_column_member_and_numpy() -> None:
    # Arrange
    rng = np.random.default_rng(13)
    cells = rng.gamma(3.0, 30.0, size=(40, 120))
    cells[2, 5] = np.nan
    cells[:, 9] = np.nan
    k = np.array([0.5, 0.9, 0.96, 0.975])
    # Act
    table, excluded = statistics.column_percentiles(cells, k, "weibull")
    # Assert
    assert table.shape == (4, 120) and excluded[5] == 1 and excluded[9] == 40
    for column in range(cells.shape[1]):
        finite = cells[np.isfinite(cells[:, column]), column]
        if finite.size == 0:
            assert np.all(np.isnan(table[:, column]))
            continue
        per_column = statistics.plotting_position_percentile(finite, k, "weibull")
        assert _bits(table[:, column]).tolist() == _bits(per_column).tolist()
        np.testing.assert_allclose(table[:, column], np.quantile(finite, k, method="weibull"), rtol=1e-12)


@pytest.mark.parametrize(("family", "parameters", "spread"), [
    ("GeneralizedExtremeValue", [100.0, 20.0, -0.1], [10.0, 15.0, 0.6]),
    ("Normal", [0.0, 1.0], [1.0, 1.0]),
    ("PearsonTypeIII", [100.0, 20.0, 0.5], [10.0, 15.0, 3.0]),
    ("LogPearsonTypeIII", [2.0, 0.2, 0.3], [1.0, 0.3, 3.0]),
    ("GeneralizedPareto", [0.0, 10.0, -0.1], [1.0, 8.0, 1.0]),
    ("KappaFour", [100.0, 20.0, -0.1, -0.2], [10.0, 15.0, 1.0, 1.0]),
])
def test_parameter_set_validity_equals_the_per_row_validation(family: str, parameters: list[float],
                                                              spread: list[float]) -> None:
    # Arrange
    distribution = distributions.create_distribution(family, parameters)
    rng = np.random.default_rng(17)
    rows = np.asarray(parameters) + rng.standard_normal((300, len(parameters))) * np.asarray(spread)
    rows[0, 1] = np.nan
    # Act
    valid = distribution.parameter_set_validity(rows)
    # Assert
    expected = [distribution.validate_parameters(row) is None for row in rows]
    assert valid.tolist() == expected
    assert 0 < valid.sum() < rows.shape[0]


def test_ordered_paired_data_bound_lookups_return_the_bits_of_the_method_group() -> None:
    # Arrange
    rng = np.random.default_rng(19)
    x = np.sort(rng.random(60)) * 100.0
    y = np.sort(rng.random(60))[::-1]
    table = data.OrderedPairedData(x, y, strict_x=False, y_order="descending", strict_y=False)
    ns = data._data()
    xt, yt, ext = ns.Transform.Logarithmic, ns.Transform.NormalZ, getattr(ns.ExtrapolationSides, "None")
    q = np.linspace(-5.0, 110.0, 50)
    # Act
    array = table.get_y_from_x(q, "logarithmic", "normal_z")
    scalar = table.get_y_from_x(42.0, "logarithmic", "normal_z")
    inverse = table.get_x_from_y(np.linspace(0.05, 0.95, 19), "logarithmic", "normal_z")
    # Assert
    expected = to_numpy_doubles(table.net.GetYFromX(to_net_doubles(q), xt, yt, ext))
    assert _bits(array).tolist() == _bits(expected).tolist()
    assert scalar == float(table.net.GetYFromX(42.0, xt, yt, ext))
    assert _bits(inverse).tolist() == _bits(to_numpy_doubles(
        table.net.GetXFromY(to_net_doubles(np.linspace(0.05, 0.95, 19)), xt, yt, ext))).tolist()


def test_runtime_config_compiles_loops_optimized_without_profile_guided_tiering() -> None:
    # Arrange
    from pynumerics import identity, runtime  # noqa: PLC0415

    initialize()
    import System  # noqa: PLC0415

    packaged = json.loads((runtime._PACKAGED_BIN / runtime.RUNTIME_CONFIG_FILE).read_text(encoding="utf-8"))
    loaded_file = runtime._resolve_runtime_config(identity.numerics_identity().location.parent)
    loaded = json.loads(loaded_file.read_text(encoding="utf-8"))["runtimeOptions"]["configProperties"]
    # Act / Assert: the packaged settings, and the running runtime reports the settings it was started with
    properties = packaged["runtimeOptions"]["configProperties"]
    assert properties["System.Runtime.TieredCompilation.QuickJitForLoops"] is False
    assert properties["System.Runtime.TieredPGO"] is False
    for name in ("System.Runtime.TieredCompilation.QuickJitForLoops", "System.Runtime.TieredPGO"):
        if name in loaded:
            assert str(System.AppContext.GetData(name)).lower() == str(loaded[name]).lower()
