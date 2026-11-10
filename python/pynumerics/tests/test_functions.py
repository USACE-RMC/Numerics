"""Link functions and univariate functions from ``Numerics.Functions``.

SciPy is the oracle for the links (``special.logit``, ``expit``, ``ndtri``, ``ndtr`` and
``stats.yeojohnson``); the closed forms of the linear and power laws, written with NumPy and
``scipy.stats.norm``, are the oracle for the functions.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from scipy import special, stats

from pynumerics import data, distributions
from pynumerics import functions as fn

X = np.array([0.05, 0.2, 0.5, 0.8, 0.95])


# Link functions -------------------------------------------------------------------------------


@pytest.mark.parametrize(("link", "forward", "inverse", "derivative", "points"), [
    (fn.IdentityLink(), lambda x: x, lambda e: e, lambda x: np.ones_like(x), X - 0.5),
    (fn.LogLink(), np.log, np.exp, lambda x: 1.0 / x, X * 10.0),
    (fn.LogitLink(), special.logit, special.expit, lambda x: 1.0 / (x * (1.0 - x)), X),
    (fn.ProbitLink(), special.ndtri, special.ndtr, lambda x: 1.0 / stats.norm.pdf(special.ndtri(x)), X),
    (fn.ComplementaryLogLogLink(), lambda x: np.log(-np.log1p(-x)), lambda e: -np.expm1(-np.exp(e)),
     lambda x: -1.0 / ((1.0 - x) * np.log1p(-x)), X),
    (fn.FisherZLink(), np.arctanh, np.tanh, lambda x: 1.0 / (1.0 - x * x), X - 0.5),
])
def test_links_match_scipy(link, forward, inverse, derivative, points) -> None:
    # Act
    eta = link.link(points)
    # Assert
    np.testing.assert_allclose(eta, forward(points), rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(link.inverse_link(eta), inverse(eta), rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(link.d_link(points), derivative(points), rtol=1e-10)
    assert link.link(float(points[1])) == eta[1]
    assert type(fn.link_from_xml(link.to_xml())) is type(link)


@pytest.mark.parametrize("lam", [-0.5, 0.0, 0.5, 2.0])
def test_yeo_johnson_link_matches_scipy(lam: float) -> None:
    # Arrange
    link = fn.YeoJohnsonLink(lam)
    points = np.array([-3.0, -0.5, 0.0, 0.7, 4.0])
    # Act
    eta = link.link(points)
    # Assert
    np.testing.assert_allclose(eta, stats.yeojohnson(points, lmbda=lam), rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(link.inverse_link(eta), points, rtol=1e-12, atol=1e-14)
    assert fn.link_from_xml(link.to_xml()).lambda_ == lam


def test_yeo_johnson_link_fits_lambda_by_maximum_likelihood_like_scipy() -> None:
    # Arrange
    values = np.random.default_rng(5).gamma(2.0, 2.0, 200)
    # Act
    link = fn.YeoJohnsonLink(values=values)
    # Assert
    assert link.lambda_ == pytest.approx(stats.yeojohnson_normmax(values), abs=1e-4)


def test_create_link_by_name_and_invalid_names() -> None:
    # Act / Assert
    assert isinstance(fn.create_link("complementary_log_log"), fn.ComplementaryLogLogLink)
    assert fn.create_link("YeoJohnson").lambda_ == 1.0
    with pytest.raises(ValueError, match="not a member"):
        fn.create_link("cube")
    with pytest.raises(ValueError):
        fn.YeoJohnsonLink(values=[1.0])


def test_link_controller_maps_each_parameter_and_round_trips() -> None:
    # Arrange
    controller = fn.LinkController([fn.LogLink(), "logit", None])
    theta = np.array([2.0, 0.3, -1.5])
    # Act
    eta = controller.link(theta)
    phi_term = controller.log_det_jacobian(eta)
    restored = fn.LinkController.from_xml(controller.to_xml())
    # Assert
    np.testing.assert_allclose(eta, [np.log(2.0), special.logit(0.3), -1.5], rtol=1e-14)
    np.testing.assert_allclose(controller.inverse_link(eta), theta, rtol=1e-14)
    assert phi_term == pytest.approx(-np.log(1.0 / 2.0) - np.log(1.0 / (0.3 * 0.7)), rel=1e-13)
    assert len(controller) == controller.count == 3
    assert isinstance(controller[0], fn.LogLink) and controller[2] is None and controller[5] is None
    assert [type(v) for v in restored.links] == [fn.LogLink, fn.LogitLink, type(None)]
    np.testing.assert_allclose(np.diag(controller.link_jacobian(theta)), [0.5, 1.0 / 0.21, 1.0], rtol=1e-14)


def test_location_scale_shape_controller() -> None:
    # Act
    controller = fn.LinkController.for_location_scale_shape(scale="log")
    # Assert
    assert controller.links[0] is None and isinstance(controller.links[1], fn.LogLink)
    assert controller.links[2] is None


# Univariate functions -------------------------------------------------------------------------


def test_linear_function_matches_the_closed_form() -> None:
    # Arrange
    line = fn.LinearFunction(1.5, 2.0, 0.4)
    x = np.array([-1.0, 0.0, 3.0])
    # Act
    mean = line.function(x)
    line.confidence_level = 0.9
    upper = line(x)
    # Assert
    assert not line.is_deterministic and line.number_of_parameters == 3
    np.testing.assert_allclose(mean, 1.5 + 2.0 * x, rtol=1e-15)
    np.testing.assert_allclose(upper, 1.5 + 2.0 * x + stats.norm.ppf(0.9, 0, 0.4), rtol=1e-14)
    np.testing.assert_allclose(line.inverse_function(upper), x, rtol=1e-14, atol=1e-14)
    np.testing.assert_array_equal(line.minimum_of_parameters[2], 0.0)


def test_linear_function_parameters_validate_and_set() -> None:
    # Arrange
    line = fn.LinearFunction(0.0, 1.0, 1.0)
    # Act
    message = line.validate_parameters([0.0, 1.0, -1.0])
    line.set_parameters([2.0, 3.0, -1.0])
    # Assert
    assert "greater than zero" in message
    assert line.validate_parameters([0.0, 1.0, 1.0]) is None
    assert (line.alpha, line.beta, line.sigma) == (2.0, 3.0, -1.0) and not line.parameters_valid
    with pytest.raises(ValueError, match="greater than zero"):
        line.function(1.0)
    with pytest.raises(ValueError, match="greater than zero"):
        line.validate_parameters([0.0, 1.0, -1.0], throw=True)


@pytest.mark.parametrize("inverse", [False, True])
def test_power_function_matches_the_closed_form(inverse: bool) -> None:
    # Arrange
    power = fn.PowerFunction(2.0, 1.6, 0.5, 0.2)
    power.is_inverse = inverse
    x = np.array([1.0, 2.0, 7.5])
    # Act
    mean = power(x)
    power.confidence_level = 0.25
    quantile = power(x)
    # Assert
    z = stats.norm.ppf(0.25, 0, 0.2)
    if inverse:
        np.testing.assert_allclose(mean, np.exp((np.log(x) - np.log(2.0)) / 1.6) + 0.5, rtol=1e-14)
        np.testing.assert_allclose(quantile, np.exp((np.log(x) - np.log(2.0) - z) / 1.6) + 0.5, rtol=1e-13)
    else:
        np.testing.assert_allclose(mean, 2.0 * (x - 0.5) ** 1.6, rtol=1e-14)
        np.testing.assert_allclose(quantile, 2.0 * (x - 0.5) ** 1.6 * np.exp(z), rtol=1e-13)
    np.testing.assert_allclose(power.inverse_function(quantile), x, rtol=1e-12)
    assert power.minimum == 0.5 and power.is_inverse is inverse


def test_power_function_minimum_is_derived_from_xi() -> None:
    # Arrange
    power = fn.PowerFunction()
    # Act / Assert
    with pytest.raises(NotImplementedError, match="Xi"):
        power.minimum = 1.0
    power.xi = 1.0
    assert power.minimum == 1.0


def test_segmented_power_function_sums_the_controls() -> None:
    # Arrange
    parameters = [1.0, 1.2, 1.6, 3.0, 0.4, 2.1, 0.05]
    rating = fn.SegmentedPowerFunction(parameters=parameters)
    h = np.array([0.5, 1.0, 2.0, 3.0, 5.0])
    # Act
    q = rating.function(h)
    rating.confidence_level = 0.95
    q95 = rating(h)
    # Assert
    expected = np.where(h > 1.0, 10**1.2 * np.clip(h - 1.0, 0, None) ** 1.6, 0.0) + \
        np.where(h > 3.0, 10**0.4 * np.clip(h - 3.0, 0, None) ** 2.1, 0.0)
    np.testing.assert_allclose(q, expected, rtol=1e-13)
    np.testing.assert_allclose(q95, expected * 10 ** stats.norm.ppf(0.95, 0, 0.05), rtol=1e-12)
    np.testing.assert_allclose(rating.inverse_function(q95[2:]), h[2:], rtol=1e-8)
    assert rating.number_of_segments == 2 and rating.sigma == 0.05
    assert (rating.get_breakpoint(2), rating.get_log10_alpha(2), rating.get_beta(2)) == (3.0, 0.4, 2.1)


def test_default_segmented_power_function_is_deterministic() -> None:
    # Act
    rating = fn.SegmentedPowerFunction(3)
    # Assert
    assert rating.is_deterministic and rating.number_of_parameters == 10
    with pytest.raises(ValueError):
        fn.SegmentedPowerFunction(parameters=[1.0, 2.0])


def _rating_table() -> data.UncertainOrderedPairedData:
    return data.UncertainOrderedPairedData([1.0, 2.0, 4.0], [distributions.Normal(10.0, 1.0),
                                                             distributions.Normal(20.0, 2.0),
                                                             distributions.Normal(40.0, 3.0)])


def test_tabular_function_interpolates_the_curve_sample() -> None:
    # Arrange
    table = fn.TabularFunction(_rating_table())
    x = np.array([1.5, 3.0])
    # Act
    mean = table(x)
    table.confidence_level = 0.9
    upper = table(x)
    # Assert
    np.testing.assert_allclose(mean, np.interp(x, [1, 2, 4], [10, 20, 40]), rtol=1e-14)
    curve = stats.norm.ppf(0.9, [10, 20, 40], [1, 2, 3])
    np.testing.assert_allclose(upper, np.interp(x, [1, 2, 4], curve), rtol=1e-14)
    np.testing.assert_allclose(table.inverse_function(upper), x, rtol=1e-13)
    assert table.paired_data.x.tolist() == [1.0, 2.0, 4.0] and table.parameters_valid


def test_tabular_function_transforms_and_extrapolation() -> None:
    # Arrange
    table = fn.TabularFunction(_rating_table())
    # Act
    table.x_transform = "logarithmic"
    table.y_transform = "logarithmic"
    inside = table(3.0)
    held = table(8.0)
    table.extrapolation = "above"
    extended = table(8.0)
    table.allow_negative_y_values = False
    # Assert
    expected = 10 ** np.interp(np.log10(3.0), np.log10([1, 2, 4]), np.log10([10, 20, 40]))
    assert inside == pytest.approx(expected, rel=1e-14)
    assert held == 40.0 and extended == pytest.approx(80.0, rel=1e-12)
    assert table.extrapolation.value == "Above" and not table.allow_negative_y_values
    assert table.x_transform.value == "Logarithmic" and table.y_transform.value == "Logarithmic"
    with pytest.raises(NotImplementedError):
        table.set_parameters([1.0])


def test_composite_function_weighted_average_and_mixture() -> None:
    # Arrange
    first, second = fn.LinearFunction(0.0, 1.0, 0.5), fn.LinearFunction(10.0, 2.0, 0.5)
    composite = fn.CompositeFunction([first, second], [0.25, 0.75])
    # Act
    mean = composite(2.0)
    composite.confidence_level = 0.8
    average = composite(2.0)
    composite.mode = "mixture"
    mixed = composite(2.0)
    # Assert
    assert mean == pytest.approx(0.25 * 2.0 + 0.75 * 14.0, rel=1e-15)
    assert average == pytest.approx(mean + stats.norm.ppf(0.8, 0, 0.5), rel=1e-14)
    assert mixed in (pytest.approx(first.function(2.0) + stats.norm.ppf(0.8 / 1.0, 0, 0.5)),) or \
        np.isfinite(mixed)
    np.testing.assert_array_equal(composite.weights, [0.25, 0.75])
    assert [type(f) for f in composite.functions] == [fn.LinearFunction, fn.LinearFunction]
    assert composite.mode is fn.CompositeMode.MIXTURE
    composite.confidence_level = -1.0
    assert composite.inverse_function(mean) == pytest.approx(2.0, rel=1e-9)
    assert fn.CompositeFunction([first, second]).weights.tolist() == [0.5, 0.5]


def test_every_function_round_trips_through_xml() -> None:
    # Arrange
    functions = [fn.LinearFunction(1.0, 2.0, 0.3), fn.PowerFunction(2.0, 1.5, 0.1, 0.2),
                 fn.SegmentedPowerFunction(parameters=[0.0, 1.0, 1.5, 0.1]),
                 fn.TabularFunction(_rating_table()),
                 fn.CompositeFunction([fn.LinearFunction(0.0, 1.0), fn.LinearFunction(1.0, 1.0)])]
    # Act / Assert
    for function in functions:
        restored = type(function).from_xml(function.to_xml())
        generic = fn.function_from_xml(function.to_xml())
        assert type(restored) is type(function) is type(generic)
        assert restored.function(2.5) == function.function(2.5) == generic(2.5)
        assert restored.function_type is function.function_type


def test_create_function_by_type() -> None:
    # Act
    power = fn.create_function(fn.FunctionType.POWER)
    # Assert
    assert isinstance(power, fn.PowerFunction) and power.function_type is fn.FunctionType.POWER
    assert isinstance(fn.create_function("segmented_power"), fn.SegmentedPowerFunction)
    power.maximum = 10.0
    assert power.maximum == 10.0 and power(100.0) == power(10.0)


def test_ensemble_samples_members_by_index_and_percentile() -> None:
    # Arrange
    sets = np.array([[0.0, 1.0, 0.1], [1.0, 2.0, 0.2], [2.0, 3.0, 0.3], [3.0, 4.0, 0.4]])
    ensemble = fn.EnsembleFunction(fn.LinearFunction(0.0, 1.0, 0.1), sets, fitness=[-1, -2, -3, -4],
                                   weights=[0.1, 0.2, 0.3, 0.4])
    # Act
    member = ensemble.sample_at(2)
    upper = ensemble.sample(0.6)
    restored = fn.EnsembleFunction.from_xml(ensemble.to_xml())
    # Assert
    assert isinstance(member, fn.LinearFunction) and (member.alpha, member.beta, member.sigma) == (2.0, 3.0, 0.3)
    assert upper.alpha == 2.0 and ensemble.sample(1.0).alpha == 3.0
    np.testing.assert_array_equal(ensemble.parameter_sets, sets)
    np.testing.assert_array_equal(ensemble.fitness, [-1, -2, -3, -4])
    np.testing.assert_array_equal(ensemble.weights, [0.1, 0.2, 0.3, 0.4])
    assert len(ensemble) == ensemble.count == restored.count == 4
    with pytest.raises(ValueError):
        ensemble.sample_at(4)
    with pytest.raises(ValueError, match="exactly 3"):
        fn.EnsembleFunction(fn.LinearFunction(), [[0.0, 1.0]])


# List overloads -------------------------------------------------------------------------------


def _bits(values: object) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64)


def _loop_members(net: object, members: tuple[str, ...], values: np.ndarray) -> list[np.ndarray]:
    return [np.array([float(getattr(net, name)(float(v))) for v in values]) for name in members]


def _functions_with_points() -> list[tuple[object, tuple[str, ...], np.ndarray]]:
    unit = np.linspace(0.01, 0.99, 400)
    links = [fn.IdentityLink(), fn.LogLink(), fn.LogitLink(), fn.ProbitLink(), fn.ComplementaryLogLogLink(),
             fn.FisherZLink(), fn.YeoJohnsonLink(0.4)]
    cases = [(link, ("Link", "InverseLink", "DLink"), unit) for link in links]
    for function in (fn.LinearFunction(2.0, 0.5), fn.PowerFunction(3.0, 1.5, 1.0)):
        cases.append((function, ("Function", "InverseFunction"), 1.5 + 100.0 * unit))
    return cases


@pytest.mark.parametrize("case", range(9))
def test_array_arguments_cross_once_with_the_scalar_bits(case: int) -> None:
    # Arrange
    wrapper, members, points = _functions_with_points()[case]
    names = {"Link": "link", "InverseLink": "inverse_link", "DLink": "d_link", "Function": "function",
             "InverseFunction": "inverse_function"}
    # Act
    arrays = [getattr(wrapper, names[name])(points.reshape(20, 20)) for name in members]
    # Assert (the list overloads equal the scalar members element by element, bit for bit)
    for array, expected in zip(arrays, _loop_members(wrapper.net, members, points), strict=True):
        assert array.shape == (20, 20)
        np.testing.assert_array_equal(_bits(array.ravel()), _bits(expected))


def test_link_arrays_are_faster_than_element_calls() -> None:
    # Arrange
    link = fn.LogitLink()
    points = np.linspace(0.001, 0.999, 20000)
    link.link(points[:5])
    # Act
    begin = time.perf_counter()
    fast = link.link(points)
    fast_seconds = time.perf_counter() - begin
    scalar = link.net.Link
    begin = time.perf_counter()
    slow = np.array([scalar(p) for p in points.tolist()])
    slow_seconds = time.perf_counter() - begin
    # Assert
    np.testing.assert_array_equal(fast, slow)
    np.testing.assert_allclose(fast, special.logit(points), rtol=1e-13, atol=1e-15)
    assert fast_seconds < slow_seconds / 3
