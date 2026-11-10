"""Paired data, uncertain paired data, line simplification and probability ordinates.

Shapely's Douglas-Peucker ``simplify``, ``numpy.trapezoid``, ``numpy.searchsorted`` and
``scipy.stats.norm`` are the oracles; Visvalingam-Whyatt and Lang are checked on small
hand-worked curves (Visvalingam and Whyatt, 1993, Cartogr. J. 30(1):46-51).
"""

from __future__ import annotations

import enum
import importlib
import inspect

import numpy as np
import pytest
from numpydoc import validate as numpydoc_validate
from scipy import stats

from pynumerics import data, distributions, errors

X = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
Y = np.array([0.0, 0.3, -0.1, 2.0, 2.2, 1.9, 5.0, 5.1])


def test_ordinate_values_transforms_and_errors() -> None:
    # Arrange
    point = data.Ordinate(100.0, 0.975)
    # Act
    moved = point.transform("logarithmic", "normal_z")
    restored = data.Ordinate.from_xml(point.to_xml())
    bad = data.Ordinate(float("nan"), 1.0)
    # Assert
    assert (moved.x, moved.y) == (2.0, pytest.approx(stats.norm.ppf(0.975), rel=1e-14))
    assert (restored.x, restored.y) == (100.0, 0.975) and point.is_valid
    assert not bad.is_valid and bad.ordinate_errors()
    assert point.ordinate_valid(data.Ordinate(200.0, 0.99))
    assert not point.ordinate_valid(data.Ordinate(50.0, 0.99))
    assert point.ordinate_errors(data.Ordinate(50.0, 0.99))
    assert repr(point) == "Ordinate(100.0, 0.975)"


def test_paired_data_members_match_numpy() -> None:
    # Arrange
    pairs = data.OrderedPairedData(X, Y, strict_y=False, y_order="none")
    # Act
    area_y = pairs.trapezoidal_area_under_y()
    increasing = np.cumsum(np.abs(Y) + 0.1)
    area_x = data.OrderedPairedData(X, increasing).trapezoidal_area_under_x()
    inverted = pairs.invert()
    # Assert
    assert area_y == pytest.approx(np.trapezoid(Y, X), rel=1e-14)
    assert area_x == pytest.approx(np.trapezoid(X, increasing), rel=1e-14)
    with pytest.raises(RuntimeError, match="sort order"):
        pairs.trapezoidal_area_under_x()
    np.testing.assert_array_equal(inverted.x, Y)
    np.testing.assert_array_equal(inverted.y, X)
    assert pairs[3].y == 2.0 and len(pairs) == pairs.count == 8
    assert [o.x for o in pairs] == X.tolist()


def test_paired_data_searches_match_numpy_searchsorted() -> None:
    # Arrange
    y = np.cumsum(np.arange(1.0, 9.0))
    pairs = data.OrderedPairedData(X, y)
    probes = np.linspace(0.1, 6.9, 15)
    expected = np.searchsorted(X, probes) - 1
    y_probes = np.interp(probes, X, y)
    # Act / Assert
    for member in ("search_x", "sequential_search_x", "bisection_search_x", "hunt_search_x"):
        np.testing.assert_array_equal([getattr(pairs, member)(v) for v in probes], expected, err_msg=member)
    for member in ("search_y", "sequential_search_y", "bisection_search_y", "hunt_search_y"):
        np.testing.assert_array_equal([getattr(pairs, member)(v) for v in y_probes], expected, err_msg=member)
    assert pairs.binary_search_x(3.0) == 3 and pairs.binary_search_x(3.5) == ~4
    assert pairs.binary_search_y(y[5]) == 5
    pairs.use_smart_search = False
    pairs.x_search_start = 2
    pairs.y_search_start = 1
    assert (pairs.use_smart_search, pairs.x_search_start, pairs.y_search_start) == (False, 2, 1)


def test_paired_data_collection_members_and_validity() -> None:
    # Arrange
    pairs = data.OrderedPairedData([1.0, 2.0], [10.0, 20.0])
    # Act
    pairs.append(3.0, 30.0)
    pairs.insert(0, 0.0, 0.0)
    pairs[1] = data.Ordinate(1.0, 11.0)
    found = pairs.index_of(2.0, 20.0)
    contained = data.Ordinate(3.0, 30.0) in pairs
    pairs.append(2.5, 5.0)
    invalid, errors_found = pairs.is_valid, pairs.errors
    removed = pairs.remove(2.5, 5.0)
    # Assert
    assert found == 2 and contained and removed
    assert not invalid and errors_found and pairs.is_valid
    np.testing.assert_array_equal(pairs.y, [0.0, 11.0, 20.0, 30.0])
    pairs.remove_at(0)
    pairs.remove_range(0, 2)
    np.testing.assert_array_equal(pairs.x, [3.0])
    pairs.clear()
    assert len(pairs) == 0


def test_paired_data_order_settings_copies_and_xml() -> None:
    # Arrange
    pairs = data.OrderedPairedData([3.0, 2.0, 1.0], [1.0, 2.0, 3.0], x_order="descending")
    # Act
    copied = pairs.copy()
    unvalidated = pairs.copy_without_validation()
    restored = data.OrderedPairedData.from_xml(pairs.to_xml())
    built = data.OrderedPairedData.from_ordinates([data.Ordinate(1.0, 1.0), data.Ordinate(2.0, 4.0)])
    pairs.strict_x = False
    pairs.order_y = "descending"
    # Assert
    assert pairs.order_x.value == "Descending" and not pairs.strict_x and pairs.strict_y
    assert pairs.order_y.value == "Descending" and not pairs.is_valid
    for other in (copied, unvalidated, restored):
        np.testing.assert_array_equal(other.x, [3.0, 2.0, 1.0])
        assert other.is_valid and other.order_y.value == "Ascending"
    np.testing.assert_array_equal(built.y, [1.0, 4.0])
    pairs.strict_y = False
    pairs.order_x = "ascending"
    assert not pairs.strict_y and pairs.order_x.value == "Ascending"


def test_douglas_peucker_matches_shapely() -> None:
    # Arrange
    shapely = pytest.importorskip("shapely")
    pairs = data.OrderedPairedData(X, Y, strict_y=False, y_order="none")
    line = shapely.LineString(np.column_stack([X, Y]))
    for tolerance in (0.05, 0.3, 1.0):
        # Act
        kept = pairs.douglas_peucker_simplify(tolerance)
        rx, ry = data.ramer_douglas_peucker(X, Y, tolerance)
        # Assert
        expected = np.asarray(shapely.simplify(line, tolerance, preserve_topology=False).coords)
        np.testing.assert_array_equal(np.column_stack([kept.x, kept.y]), expected, err_msg=str(tolerance))
        np.testing.assert_array_equal(np.column_stack([rx, ry]), expected, err_msg=str(tolerance))


def test_visvalingam_whyatt_drops_the_smallest_effective_areas_first() -> None:
    # Arrange (the triangle areas at x = 1, 2, 3 are 0.05, 1.0 and 0.75)
    pairs = data.OrderedPairedData([0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 0.1, 1.0, 0.0, 0.0], strict_y=False,
                                   y_order="none")
    # Act
    four = pairs.visvalingam_whyatt_simplify(4)
    three = pairs.visvalingam_whyatt_simplify(3)
    # Assert
    np.testing.assert_array_equal(four.x, [0.0, 2.0, 3.0, 4.0])
    np.testing.assert_array_equal(three.x, [0.0, 2.0, 4.0])


def test_lang_keeps_the_points_that_leave_the_chord() -> None:
    # Arrange
    x = np.arange(9.0)
    y = np.array([0.0, 0.01, 0.0, 0.01, 3.0, 0.01, 0.0, 0.01, 0.0])
    pairs = data.OrderedPairedData(x, y, strict_y=False, y_order="none")
    # Act
    kept = pairs.lang_simplify(0.1, 4)
    same = pairs.lang_simplify(0.0, 4)
    # Assert
    assert 4.0 in kept.x.tolist() and kept.x[0] == 0.0 and kept.x[-1] == 8.0 and len(kept) < 9
    np.testing.assert_array_equal(same.y, y)


def test_perpendicular_distance_matches_the_cross_product() -> None:
    # Arrange
    point, start, end = np.array([2.0, 3.0]), np.array([-1.0, 0.5]), np.array([4.0, -1.0])
    # Act
    distance = data.perpendicular_distance(point, start, end)
    on_ordinates = data.perpendicular_distance(data.Ordinate(*point), data.Ordinate(*start), data.Ordinate(*end))
    # Assert
    d, v = end - start, point - start
    expected = abs(d[0] * v[1] - d[1] * v[0]) / np.hypot(*d)
    assert distance == pytest.approx(expected, rel=1e-14) and on_ordinates == distance
    with pytest.raises(ValueError, match="Not enough points"):
        data.ramer_douglas_peucker([1.0], [1.0], 0.1)


def _normals() -> list:
    return [distributions.Normal(10.0, 1.0), distributions.Normal(20.0, 2.0), distributions.Normal(40.0, 3.0)]


def test_uncertain_curve_samples_are_the_quantiles_of_each_distribution() -> None:
    # Arrange
    table = data.UncertainOrderedPairedData([1.0, 2.0, 4.0], _normals())
    # Act
    upper = table.curve_sample(0.95)
    mean = table.curve_sample()
    # Assert
    np.testing.assert_allclose(upper.y, stats.norm.ppf(0.95, [10, 20, 40], [1, 2, 3]), rtol=1e-14)
    np.testing.assert_array_equal(mean.y, [10.0, 20.0, 40.0])
    np.testing.assert_array_equal(table.x, [1.0, 2.0, 4.0])
    assert table.distribution == "Normal" and table.is_valid and not table.errors
    assert [d.mean for d in table.distributions] == [10.0, 20.0, 40.0]


def test_uncertain_paired_data_validates_the_quantile_order() -> None:
    # Arrange
    crossing = [distributions.Normal(10.0, 1.0), distributions.Normal(10.5, 5.0)]
    # Act
    table = data.UncertainOrderedPairedData([1.0, 2.0], crossing)
    # Assert
    assert not table.is_valid and table.errors
    table.validation_tail_probability = 0.4
    assert table.validation_tail_probability == 0.4
    table.allow_different_distribution_types = True
    assert table.allow_different_distribution_types


def test_uncertain_paired_data_collection_members() -> None:
    # Arrange
    table = data.UncertainOrderedPairedData([1.0, 2.0, 4.0], _normals())
    # Act
    table.append(5.0, distributions.Normal(50.0, 3.0))
    table.extend([6.0, 7.0], [distributions.Normal(60.0, 3.0), distributions.Normal(70.0, 3.0)])
    table.insert(0, 0.5, distributions.Normal(5.0, 1.0))
    table.insert_range(1, [0.7, 0.8], [distributions.Normal(6.0, 1.0), distributions.Normal(7.0, 1.0)])
    first = table[0]
    position = table.index_of(first)
    contained = first in table
    table[0] = data.UncertainOrdinate(0.4, distributions.Normal(4.0, 1.0))
    removed = table.remove(table[0])
    table.remove_at(0)
    table.remove_range(0, 1)
    table.remove_range([0, 2])
    # Assert
    assert position == 0 and contained and removed
    np.testing.assert_array_equal(table.x, [2.0, 5.0, 6.0, 7.0])
    copied = table.copy()
    table.clear()
    assert len(table) == 0 and len(copied) == copied.count == 4


def test_uncertain_paired_data_settings_and_xml() -> None:
    # Arrange
    table = data.UncertainOrderedPairedData([1.0, 2.0, 4.0], _normals(), distribution_type="Normal")
    # Act
    restored = data.UncertainOrderedPairedData.from_xml(table.to_xml())
    table.strict_x = False
    table.strict_y = False
    table.order_x = "ascending"
    table.order_y = "ascending"
    table.validate()
    # Assert
    np.testing.assert_array_equal(restored.curve_sample(0.9).y, table.curve_sample(0.9).y)
    assert not table.strict_x and not table.strict_y and table.order_x.value == table.order_y.value == "Ascending"
    assert [o.x for o in table] == [1.0, 2.0, 4.0]


def test_uncertain_ordinate_members() -> None:
    # Arrange
    low = data.UncertainOrdinate(1.0, distributions.Normal(10.0, 1.0))
    high = data.UncertainOrdinate(2.0, distributions.Normal(20.0, 2.0))
    # Act
    restored = data.UncertainOrdinate.from_xml(low.to_xml())
    legacy = data.UncertainOrdinate.from_xml('<UncertainOrdinate X="3" Mu="30" Sigma="2" />', "Normal")
    # Assert
    assert low.get_ordinate().y == 10.0
    assert low.get_ordinate(0.9).y == pytest.approx(stats.norm.ppf(0.9, 10, 1), rel=1e-14)
    assert low.ordinate_valid(high) and not high.ordinate_valid(low)
    assert low.ordinate_valid(high, tail_probability=0.01)
    assert high.ordinate_errors(low) and not low.ordinate_errors() and low.is_valid
    assert restored.x == 1.0 and restored.y.mean == 10.0 and legacy.x == 3.0
    assert (legacy.y.mean, legacy.y.standard_deviation) == (30.0, 2.0)


def test_probability_ordinates_defaults_parse_and_validation() -> None:
    # Arrange
    ordinates = data.ProbabilityOrdinates()
    # Act
    text = ordinates.to_delimited_string()
    parsed = data.ProbabilityOrdinates.parse(text)
    custom = data.ProbabilityOrdinates([0.5, 0.1])
    # Assert
    assert len(ordinates) == 25 and ordinates[0] == 1e-6 and ordinates[24] == 0.99
    np.testing.assert_array_equal(parsed.values, ordinates.values)
    assert ordinates.validate() == (True, [])
    valid, messages = custom.validate()
    assert not valid and messages
    custom.from_delimited_string("0.2;0.4", ";")
    np.testing.assert_array_equal(custom.values, [0.2, 0.4])
    custom.reset_to_defaults()
    assert len(custom) == 25 and data.ProbabilityOrdinates.DEFAULT_DELIMITER == "|"


def test_probability_ordinates_list_members() -> None:
    # Arrange
    ordinates = data.ProbabilityOrdinates([0.1, 0.2])
    # Act
    ordinates.append(0.5)
    ordinates.extend([0.8, 0.9])
    ordinates.insert(0, 0.01)
    ordinates.insert_range(1, [0.02, 0.05])
    ordinates[0] = 0.001
    removed = ordinates.remove(0.5)
    ordinates.remove_at(0)
    ordinates.remove_range(0, 2)
    # Assert
    assert removed
    np.testing.assert_array_equal(ordinates.values, [0.1, 0.2, 0.8, 0.9])
    ordinates.clear()
    assert len(ordinates) == 0 and ordinates.to_delimited_string(",") == ""


def test_linear_regression_exposes_the_design_and_responses() -> None:
    # Arrange
    x = np.arange(10.0)
    y = 1.0 + 2.0 * x + np.random.default_rng(8).normal(0, 0.1, 10)
    # Act
    model = data.LinearRegression(x, y)
    # Assert
    np.testing.assert_array_equal(model.x, np.column_stack([np.ones(10), x]))
    np.testing.assert_array_equal(model.y, y)


# Docstrings and translation of the data, function and time series modules -----------------------

MODULES = ("pynumerics.data", "pynumerics.functions", "pynumerics.timeseries")
CHECKS = frozenset({"GL02", "GL03", "GL05", "GL06", "GL07", "GL08", "GL09", "GL10", "SS01", "SS02", "SS03", "SS04",
                    "SS05", "SS06", "PR01", "PR02", "PR03", "PR04", "PR05", "PR06", "PR07", "PR08", "PR09", "PR10",
                    "RT01", "RT02", "RT03", "RT04", "RT05", "YD01"})


def _public(module_name: str) -> list[tuple[str, object]]:
    module = importlib.import_module(module_name)
    names = []
    for name in module.__all__:
        obj = getattr(module, name)
        if getattr(obj, "__module__", None) != module_name or not (inspect.isfunction(obj) or inspect.isclass(obj)):
            continue
        if inspect.isclass(obj) and issubclass(obj, enum.Enum):
            continue
        names.append((f"{module_name}.{name}", obj))
        if inspect.isclass(obj):
            for member, value in vars(obj).items():
                if (not member.startswith("_") or member == "__init__") and (
                        inspect.isfunction(value) or isinstance(value, (property, classmethod, staticmethod))):
                    names.append((f"{module_name}.{name}.{member}", value))
    return names


PUBLIC = [entry for module in MODULES for entry in _public(module)]


@pytest.mark.parametrize("name", [n for n, _ in PUBLIC if not n.endswith(".__init__")])
def test_public_docstrings_pass_numpydoc(name: str) -> None:
    # Act
    report = numpydoc_validate.validate(name)
    # Assert
    found = [f"{code}: {message}" for code, message in report["errors"] if code in CHECKS]
    assert not found, f"{name}\n" + "\n".join(found)


@pytest.mark.parametrize(("name", "obj"), [e for e in PUBLIC if not inspect.isclass(e[1])], ids=lambda v: str(v))
def test_public_members_translate_net_exceptions(name: str, obj: object) -> None:
    # Assert
    assert errors.is_translated(obj), name
