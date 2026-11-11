"""Integration, differentiation, root finding and ODE wrappers against SciPy and closed forms.

Every wrapper takes a Python callable, C# source or a .NET delegate; the compiled and Python forms
of the same function must give the same bits, since Numerics runs the same arithmetic on the same
values either way.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
from scipy import integrate, optimize

from pynumerics import callbacks, differentiation, integration, ode, roots, sampling
from pynumerics._util import namespace
from pynumerics.integration import IntegrationStatus

needs_roslyn = pytest.mark.skipif(not callbacks.roslyn_available(), reason="Roslyn assemblies not staged")


def _gauss(x: float) -> float:
    return math.exp(-x * x)


def _bits(values: object) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64)


# Fixed rules of the static Integration class.

def test_fixed_rules_match_scipy_on_the_same_nodes() -> None:
    # Arrange
    a, b, n = -1.0, 2.0, 200
    grid = np.linspace(a, b, n + 1)
    mid = a + (np.arange(n) + 0.5) * (b - a) / n
    fine = np.linspace(a, b, 2 * n + 1)
    # Act
    gl10 = integration.gauss_legendre(_gauss, a, b)
    gl20 = integration.gauss_legendre_20(_gauss, a, b)
    trap = integration.trapezoidal_rule(_gauss, a, b, n)
    simp = integration.simpsons_rule(_gauss, a, b, n)
    midp = integration.midpoint_rule(_gauss, a, b, n)
    # Assert
    assert gl10 == pytest.approx(integrate.fixed_quad(np.vectorize(_gauss), a, b, n=10)[0], rel=1e-14)
    assert gl20 == pytest.approx(integrate.fixed_quad(np.vectorize(_gauss), a, b, n=20)[0], rel=1e-14)
    assert trap == pytest.approx(integrate.trapezoid(np.exp(-grid ** 2), grid), rel=1e-13)
    assert simp == pytest.approx(integrate.simpson(np.exp(-fine ** 2), x=fine), rel=1e-13)
    assert midp == pytest.approx(float(np.sum(np.exp(-mid ** 2)) * (b - a) / n), rel=1e-13)


def test_gauss_legendre_is_exact_for_polynomials_of_its_degree() -> None:
    # Act / Assert
    assert integration.gauss_legendre(lambda x: x ** 19 + x ** 18, 0.0, 1.0) == pytest.approx(1 / 20 + 1 / 19, 1e-14)
    assert integration.gauss_legendre_20(lambda x: x ** 39, 0.0, 1.0) == pytest.approx(1 / 40, rel=1e-13)


# Integrator subclasses.

ONE_D = [
    integration.integrate_trapezoidal, integration.integrate_simpson, integration.integrate_adaptive_simpson,
    integration.integrate_adaptive_gauss_kronrod, integration.integrate_adaptive_gauss_lobatto,
]


@pytest.mark.parametrize("integrate_fn", ONE_D, ids=lambda fn: fn.__name__)
def test_one_dimensional_integrators_match_scipy_quad(integrate_fn) -> None:
    # Arrange
    expected = integrate.quad(_gauss, -3.0, 3.0, epsabs=1e-13, epsrel=1e-13)[0]
    # Act
    result = integrate_fn(_gauss, -3.0, 3.0, absolute_tolerance=1e-10, relative_tolerance=1e-10)
    # Assert
    assert result.value == pytest.approx(expected, rel=1e-9)
    assert result.status is IntegrationStatus.SUCCESS
    assert result.function_evaluations > 0


def test_standard_error_is_reported_only_by_integrators_that_estimate_one() -> None:
    # Act
    kronrod = integration.integrate_adaptive_gauss_kronrod(math.sin, 0.0, 3.0)
    lobatto = integration.integrate_adaptive_gauss_lobatto(math.sin, 0.0, 3.0)
    # Assert
    assert math.isfinite(kronrod.standard_error) and kronrod.standard_error >= 0.0
    assert math.isnan(lobatto.standard_error)


def test_gauss_kronrod_qags_matches_scipy_quad() -> None:
    # Arrange: endpoint singularities, where QAGS extrapolates; scipy.integrate.quad runs QUADPACK DQAGS.
    cases = [(lambda x: 1.0 / math.sqrt(x), 2.0), (math.log, -1.0), (lambda x: math.log(x) / math.sqrt(x), -4.0)]
    for f, exact in cases:
        # Act
        qags = integration.integrate_adaptive_gauss_kronrod(f, 0.0, 1.0, method="qags")
        local = integration.integrate_adaptive_gauss_kronrod(f, 0.0, 1.0)
        value, _, info = integrate.quad(f, 0.0, 1.0, epsabs=1e-8, epsrel=1e-8, full_output=1)[:3]
        # Assert
        assert qags.status is IntegrationStatus.SUCCESS
        assert qags.value == pytest.approx(exact, abs=1e-12)
        assert qags.value == pytest.approx(value, abs=1e-12)
        assert qags.function_evaluations == info["neval"]
        assert abs(qags.value - exact) < abs(local.value - exact)
        assert local.function_evaluations > qags.function_evaluations


def test_gauss_kronrod_method_default_and_names() -> None:
    # Arrange
    f = lambda x: math.exp(-x * x)  # noqa: E731
    # Act
    default = integration.integrate_adaptive_gauss_kronrod(f, -3.0, 3.0)
    local = integration.integrate_adaptive_gauss_kronrod(
        f, -3.0, 3.0, method=integration.AdaptiveGaussKronrodMethod.LOCAL_BISECTION)
    qags = integration.integrate_adaptive_gauss_kronrod(f, -3.0, 3.0, method="Qags", relative_tolerance=0.0,
                                                        absolute_tolerance=1e-12)
    # Assert
    assert default == local
    assert qags.value == pytest.approx(math.sqrt(math.pi) * math.erf(3.0), abs=1e-12)
    with pytest.raises(ValueError):
        integration.integrate_adaptive_gauss_kronrod(f, -3.0, 3.0, method="romberg")
    with pytest.raises(ValueError):
        integration.integrate_adaptive_gauss_kronrod(f, -3.0, 3.0, relative_tolerance=0.0, absolute_tolerance=1e-12)


def test_gauss_kronrod_qags_recorder_partitions_the_range() -> None:
    # Arrange
    nodes = []
    # Act
    result = integration.integrate_adaptive_gauss_kronrod(_gauss, -3.0, 3.0, method="qags",
                                                          recorder=lambda x, w, fx: nodes.append((x, w, fx)))
    # Assert
    weights = np.array([w for _, w, _ in nodes])
    assert len(nodes) % 21 == 0
    assert weights.sum() == pytest.approx(6.0, rel=1e-14)
    assert sum(w * fx for _, w, fx in nodes) == pytest.approx(result.value, rel=1e-14)


def test_gauss_kronrod_recorder_reproduces_the_result() -> None:
    # Arrange
    nodes = []
    # Act
    result = integration.integrate_adaptive_gauss_kronrod(_gauss, -3.0, 3.0,
                                                          recorder=lambda x, w, fx: nodes.append((x, w, fx)))
    # Assert
    weights = np.array([w for _, w, _ in nodes])
    assert len(nodes) % 21 == 0
    assert weights.sum() == pytest.approx(6.0, rel=1e-14)
    assert sum(w * fx for _, w, fx in nodes) == pytest.approx(result.value, rel=1e-14)


@pytest.mark.parametrize("integrate_fn", [integration.integrate_adaptive_gauss_kronrod,
                                          integration.integrate_adaptive_simpson], ids=lambda fn: fn.__name__)
def test_bins_integrate_over_stratification_bins(integrate_fn) -> None:
    # Arrange
    bins = sampling.stratify_x_values(-3.0, 3.0, 4)
    expected = integrate.quad(_gauss, -3.0, 3.0, epsabs=1e-13, epsrel=1e-13)[0]
    # Act
    result = integrate_fn(_gauss, -3.0, 3.0, bins=bins)
    # Assert
    assert result.value == pytest.approx(expected, rel=1e-8)


@pytest.mark.parametrize("integrate_fn", [integration.integrate_adaptive_simpson_2d,
                                          integration.integrate_adaptive_gauss_kronrod_2d],
                         ids=lambda fn: fn.__name__)
def test_two_dimensional_integrators_match_scipy_dblquad(integrate_fn) -> None:
    # Arrange
    def f(x: float, y: float) -> float:
        return math.exp(-x * x - 0.5 * y * y) * math.cos(x * y)

    expected = integrate.dblquad(lambda y, x: f(x, y), -1.0, 2.0, -1.5, 1.0, epsabs=1e-13, epsrel=1e-13)[0]
    # Act
    result = integrate_fn(f, -1.0, 2.0, -1.5, 1.0, absolute_tolerance=1e-10, relative_tolerance=1e-10)
    # Assert
    assert result.value == pytest.approx(expected, rel=1e-8)
    assert result.status is IntegrationStatus.SUCCESS


def test_gauss_kronrod_2d_recorder_weights_cover_the_rectangle() -> None:
    # Arrange
    total = []
    # Act
    result = integration.integrate_adaptive_gauss_kronrod_2d(
        lambda x, y: x * x + y, 0.0, 1.0, 0.0, 2.0, recorder=lambda x, y, w, fxy: total.append((w, fxy)))
    # Assert
    assert sum(w for w, _ in total) == pytest.approx(2.0, rel=1e-13)
    assert sum(w * v for w, v in total) == pytest.approx(result.value, rel=1e-13)
    assert result.value == pytest.approx(2.0 / 3.0 + 2.0, rel=1e-13)


def _sum_of_squares(x: np.ndarray) -> float:
    return float(np.sum(x * x))


@pytest.mark.parametrize(("integrate_fn", "options"), [
    (integration.integrate_monte_carlo, {"seed": 7, "max_iterations": 40000}),
    (integration.integrate_monte_carlo, {"use_sobol_sequence": True, "max_iterations": 40000}),
    (integration.integrate_miser, {"seed": 7, "max_function_evaluations": 40000}),
    (integration.integrate_miser, {"use_sobol_sequence": False, "seed": 3, "fraction": 0.2, "dither": 0.05,
                                   "minimum_number_of_subregion_points": 20, "minimum_number_of_bisections": 80,
                                   "max_function_evaluations": 40000}),
], ids=["monte_carlo", "monte_carlo_sobol", "miser", "miser_options"])
def test_sampling_integrators_match_the_closed_form(integrate_fn, options) -> None:
    # Arrange: the integral of |x|^2 over the unit cube in three dimensions is 1.
    lower, upper = [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]
    # Act
    first = integrate_fn(_sum_of_squares, lower, upper, **options)
    second = integrate_fn(_sum_of_squares, lower, upper, **options)
    # Assert
    assert first.value == pytest.approx(1.0, abs=max(5.0 * first.standard_error, 1e-3))
    assert first == second


def test_vegas_matches_the_closed_form_and_reports_its_grid() -> None:
    # Arrange
    lower, upper = [0.0, 0.0], [1.0, 2.0]
    # Act
    result = integration.integrate_vegas(lambda x, w: float(np.exp(-x[0]) * x[1]), lower, upper)
    seeded = integration.integrate_vegas(lambda x, w: float(np.exp(-x[0]) * x[1]), lower, upper, seed=5,
                                         use_sobol_sequence=False, function_calls=5000, number_of_bins=60,
                                         alpha=1.4, independent_evaluations=500, check_convergence=True,
                                         relative_tolerance=1e-3)
    # Assert
    expected = (1.0 - math.exp(-1.0)) * 2.0
    assert result.value == pytest.approx(expected, abs=5.0 * result.standard_error + 1e-6)
    assert seeded.value == pytest.approx(expected, abs=5.0 * seeded.standard_error + 1e-6)
    assert result.grid.shape == (2, 50) and seeded.grid.shape == (2, 60)
    assert math.isfinite(result.chi_squared)


def test_vegas_rare_event_configuration_estimates_a_tail_probability() -> None:
    # Arrange: P(U1 + U2 > 1.98) = 0.0002 for independent uniforms.
    def tail(x: np.ndarray, w: float) -> float:
        return 1.0 if x[0] + x[1] > 1.98 else 0.0

    # Act
    result = integration.integrate_vegas(tail, [0.0, 0.0], [1.0, 1.0], rare_event_probability=2e-4,
                                         sobol_seed=11, tail_focus_parameter=None)
    # Assert
    assert result.value == pytest.approx(2e-4, rel=0.25)


def test_settings_reach_numerics_and_unknown_settings_raise() -> None:
    # Act
    limited = integration.integrate_trapezoidal(_gauss, -3.0, 3.0, max_function_evaluations=20)
    capped = integration.integrate_simpson(_gauss, -3.0, 3.0, max_iterations=3, min_iterations=1,
                                           min_function_evaluations=1)
    # Assert
    assert limited.status is IntegrationStatus.MAXIMUM_FUNCTION_EVALUATIONS_REACHED
    assert capped.status is IntegrationStatus.MAXIMUM_ITERATIONS_REACHED and capped.iterations == 3
    with pytest.raises(TypeError, match="Unknown integrator settings: tolerance"):
        integration.integrate_trapezoidal(_gauss, 0.0, 1.0, tolerance=1e-6)
    with pytest.raises(ValueError, match="relative tolerance"):
        integration.integrate_adaptive_gauss_kronrod(_gauss, 0.0, 1.0, relative_tolerance=2.0)
    with pytest.raises(ValueError, match="maximum value"):
        integration.integrate_adaptive_gauss_kronrod(_gauss, 1.0, 0.0)


def test_failures_follow_report_failure() -> None:
    # Arrange
    def broken(x: float) -> float:
        raise ZeroDivisionError("integrand failed")

    # Act
    quiet = integration.integrate_trapezoidal(broken, 0.0, 1.0, report_failure=False)
    # Assert
    assert quiet.status is IntegrationStatus.FAILURE and math.isnan(quiet.value)
    with pytest.raises(Exception, match="integrand failed"):
        integration.integrate_trapezoidal(broken, 0.0, 1.0)


def test_integrators_apply_the_parallel_policy(monkeypatch: pytest.MonkeyPatch) -> None:
    # Arrange
    seen = []
    original = callbacks.apply_parallel_policy
    monkeypatch.setattr(callbacks, "apply_parallel_policy",
                        lambda net, *fns: seen.append((type(net).__name__, len(fns))) or original(net, *fns))
    # Act
    integration.integrate_vegas(lambda x, w: 1.0, [0.0], [1.0], function_calls=200)
    integration.integrate_miser(lambda x: 1.0, [0.0], [1.0], max_function_evaluations=200)
    integration.integrate_adaptive_gauss_kronrod(_gauss, 0.0, 1.0, recorder=lambda x, w, fx: None)
    # Assert
    assert seen == [("Vegas", 1), ("Miser", 1), ("AdaptiveGaussKronrod", 2)]


# Python and compiled C# delegates give identical bits.

@needs_roslyn
def test_python_and_compiled_csharp_integrands_give_identical_bits() -> None:
    # Arrange
    def py2(x: float, y: float) -> float:
        return math.exp(-x * x) * math.cos(y)

    cs2 = "(x, y) => Math.Exp(-x * x) * Math.Cos(y)"
    # Act
    pairs = [
        (integration.integrate_adaptive_gauss_kronrod(_gauss, -3.0, 3.0),
         integration.integrate_adaptive_gauss_kronrod("x => Math.Exp(-x * x)", -3.0, 3.0)),
        (integration.integrate_adaptive_simpson(_gauss, -3.0, 3.0),
         integration.integrate_adaptive_simpson(callbacks.compile_csharp("x => Math.Exp(-x * x)"), -3.0, 3.0)),
        (integration.integrate_adaptive_gauss_kronrod_2d(py2, 0.0, 1.0, 0.0, 1.0),
         integration.integrate_adaptive_gauss_kronrod_2d(cs2, 0.0, 1.0, 0.0, 1.0)),
        (integration.integrate_miser(lambda x: float(x[0] * x[1]), [0.0, 0.0], [1.0, 1.0],
                                     max_function_evaluations=5000),
         integration.integrate_miser("x => x[0] * x[1]", [0.0, 0.0], [1.0, 1.0], max_function_evaluations=5000)),
    ]
    vegas_py = integration.integrate_vegas(lambda x, w: float(x[0] * x[1]), [0.0, 0.0], [1.0, 1.0])
    vegas_cs = integration.integrate_vegas("(x, w) => x[0] * x[1]", [0.0, 0.0], [1.0, 1.0])
    # Assert
    for py, cs in pairs:
        assert _bits(py.value) == _bits(cs.value) and py == cs
    assert _bits([vegas_py.value, vegas_py.standard_error]).tolist() == \
        _bits([vegas_cs.value, vegas_cs.standard_error]).tolist()
    assert np.array_equal(_bits(vegas_py.grid), _bits(vegas_cs.grid))


@needs_roslyn
def test_python_and_compiled_csharp_roots_derivatives_and_odes_give_identical_bits() -> None:
    # Arrange
    def cube(x: float) -> float:
        return x * x * x - 2.0

    # Act
    py = [
        roots.bisection_solve(cube, 1.0, 0.0, 2.0),
        roots.secant_solve(cube, 0.0, 2.0),
        roots.newton_raphson_solve(cube, lambda x: 3.0 * x * x, 1.0),
        differentiation.ridders(math.exp, 0.5)[0],
        differentiation.hessian(lambda t: float(t[0] * t[0] * t[1]), [1.0, 2.0])[0, 1],
        ode.runge_kutta_fehlberg(lambda t, y: -2.0 * t * y, 1.0, 0.0, 1.0, 1e-5, 1e-9),
    ]
    cs = [
        roots.bisection_solve("x => x * x * x - 2.0", 1.0, 0.0, 2.0),
        roots.secant_solve("x => x * x * x - 2.0", 0.0, 2.0),
        roots.newton_raphson_solve("x => x * x * x - 2.0", "x => 3.0 * x * x", 1.0),
        differentiation.ridders("x => Math.Exp(x)", 0.5)[0],
        differentiation.hessian("t => t[0] * t[0] * t[1]", [1.0, 2.0])[0, 1],
        ode.runge_kutta_fehlberg("(t, y) => -2.0 * t * y", 1.0, 0.0, 1.0, 1e-5, 1e-9),
    ]
    # Assert
    assert _bits(py).tolist() == _bits(cs).tolist()


def test_dotnet_delegates_pass_through_unchanged() -> None:
    # Arrange
    namespace("Numerics")
    import System  # noqa: PLC0415

    delegate = System.Func[System.Double, System.Double](_gauss)
    # Act
    result = integration.integrate_adaptive_gauss_kronrod(delegate, -3.0, 3.0)
    # Assert
    assert result == integration.integrate_adaptive_gauss_kronrod(_gauss, -3.0, 3.0)


# Root finding.

ROOT_CASES = [
    (lambda x: x * x - 2.0, lambda x: 2.0 * x, 0.0, 4.0, 1.0),
    (lambda x: 2.0 * math.sin(x) - 3.0 * math.cos(x) - 0.5, lambda x: 2.0 * math.cos(x) + 3.0 * math.sin(x),
     0.0, math.pi, 0.5),
    (lambda x: math.exp(-x) - x, lambda x: -math.exp(-x) - 1.0, -2.0, 2.0, 1.0),
]


@pytest.mark.parametrize(("f", "df", "lower", "upper", "guess"), ROOT_CASES, ids=["sqrt2", "trig", "omega"])
def test_root_finders_match_scipy(f, df, lower, upper, guess) -> None:
    # Arrange
    expected = optimize.brentq(f, lower, upper, xtol=1e-15, rtol=4 * np.finfo(float).eps)
    # Act
    found = [
        roots.bisection_solve(f, guess, lower, upper),
        roots.newton_raphson_solve(f, df, guess),
        roots.newton_raphson_robust_solve(f, df, guess, lower, upper),
    ]
    secant = roots.secant_solve(f, lower, upper)
    # Assert: the secant method is not bracketed, so it may converge to another root of f.
    assert found == pytest.approx([expected] * 3, abs=1e-8)
    assert f(secant) == pytest.approx(0.0, abs=1e-12)
    assert roots.newton_raphson_solve(f, df, guess, tolerance=1e-14) == pytest.approx(
        optimize.newton(f, guess, fprime=df, tol=1e-14), abs=1e-13)


BRENTQ_CASES = [
    (lambda x: x * x - 2.0, 0.0, 2.0),
    (lambda x: x * x * x - x - 1.0, 1.0, 2.0),
    (lambda x: (x - 1.0) * (x - 1.0) * (x - 1.0), 0.0, 3.0),
    (lambda x: (x * x - 3.0) / (x + 1.0), 0.0, 5.0),
    (lambda x: math.exp(-x) - x, 0.0, 1.0),
]


@pytest.mark.parametrize(("f", "lower", "upper"), BRENTQ_CASES, ids=["sqrt2", "cubic", "triple", "rational", "omega"])
@pytest.mark.parametrize(("xtol", "rtol"), [(2e-12, None), (1e-12, 1e-10), (1e-300, 1e-14)])
def test_roots_brent_solve_is_bitwise_scipy_brentq(f, lower, upper, xtol, rtol) -> None:
    # Arrange
    kwargs = {} if rtol is None else {"rtol": rtol}
    expected = optimize.brentq(f, lower, upper, xtol=xtol, maxiter=1000, **kwargs)
    # Act
    found = roots.brent_solve(f, lower, upper, xtol, 1000,
                              **({} if rtol is None else {"relative_tolerance": rtol}))
    # Assert
    assert np.float64(found).view(np.uint64) == np.float64(expected).view(np.uint64)


def test_roots_brent_relative_tolerance_constant_and_validation() -> None:
    # Arrange
    net = namespace("Numerics.Mathematics.RootFinding").Brent.MinimumRelativeTolerance
    # Assert
    assert roots.BRENT_MINIMUM_RELATIVE_TOLERANCE == 4 * np.finfo(float).eps == float(net)
    with pytest.raises(ValueError):
        roots.brent_solve(lambda x: x - 0.5, 0.0, 1.0, relative_tolerance=1e-16)


@pytest.mark.parametrize(("f", "df", "lower", "upper", "guess"), ROOT_CASES, ids=["sqrt2", "trig", "omega"])
def test_roots_brent_solve_matches_scipy_brentq(f, df, lower, upper, guess) -> None:
    # Arrange
    expected = optimize.brentq(f, lower, upper, xtol=1e-15, rtol=4 * np.finfo(float).eps)
    # Act
    found = roots.brent_solve(f, lower, upper, tolerance=1e-12)
    # Assert
    assert "brent_solve" in roots.__all__ and "brent_bracket" in roots.__all__
    assert found == pytest.approx(expected, abs=1e-11)
    assert np.float64(found).view(np.uint64) == np.float64(roots.brent_solve(f, lower, upper, 1e-12)).view(np.uint64)


def test_newton_raphson_system_matches_scipy_fsolve() -> None:
    # Arrange
    def f(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] ** 2 + x[1] ** 2 - 4.0, math.exp(x[0]) + x[1] - 1.0])

    def jac(x: np.ndarray) -> np.ndarray:
        return np.array([[2.0 * x[0], 2.0 * x[1]], [math.exp(x[0]), 1.0]])

    expected = optimize.fsolve(f, [1.0, -1.0], fprime=jac, xtol=1e-14)
    # Act
    solution = roots.newton_raphson_solve_system(f, jac, [1.0, -1.0])
    # Assert
    assert solution == pytest.approx(expected, abs=1e-10)


def test_root_finders_raise_on_bad_brackets_and_honour_report_failure() -> None:
    # Act / Assert
    with pytest.raises(ValueError):
        roots.bisection_solve(lambda x: x * x + 1.0, 0.5, 0.0, 1.0)
    with pytest.raises(ValueError, match="failed to find root"):
        roots.newton_raphson_solve(lambda x: x * x + 1.0, lambda x: 2.0 * x, 3.0, max_iterations=20)
    last = roots.newton_raphson_solve(lambda x: x * x + 1.0, lambda x: 2.0 * x, 3.0, max_iterations=20,
                                      report_failure=False)
    assert math.isfinite(last)


# Differentiation.

def test_single_variable_derivatives_match_the_closed_form() -> None:
    # Arrange
    x = 0.7
    exact = math.exp(x)
    # Act
    first = [differentiation.forward_difference(math.exp, x), differentiation.backward_difference(math.exp, x),
             differentiation.central_difference(math.exp, x), differentiation.derivative(math.exp, x)]
    value, error = differentiation.ridders(math.exp, x)
    second = [differentiation.second_derivative(math.exp, x),
              differentiation.second_derivative_forward(math.exp, x),
              differentiation.second_derivative_backward(math.exp, x)]
    # Assert
    assert first[:2] == pytest.approx([exact] * 2, rel=1e-7)
    assert first[2:] == pytest.approx([exact] * 2, rel=1e-9)
    # The default large step 0.1 (1 + |x|) gives about 1e-15; a small step such as calculate_step_size leaves 1e-8
    assert value == pytest.approx(exact, rel=2e-14) and error > 0.0
    small = differentiation.ridders(math.exp, x, step_size=differentiation.calculate_step_size(x))[0]
    assert abs(small - exact) > 100 * abs(value - exact)
    assert differentiation.ridders(math.exp, x, step_size=0.1)[0] == pytest.approx(exact, rel=1e-12)
    assert second == pytest.approx([exact] * 3, rel=1e-4)
    assert differentiation.central_difference(math.exp, x, step_size=1e-3) == pytest.approx(
        (math.exp(x + 1e-3) - math.exp(x - 1e-3)) / 2e-3, rel=1e-15)


def test_step_size_follows_the_documented_formula() -> None:
    # Act / Assert
    eps = float(namespace("Numerics").Tools.DoubleMachineEpsilon)
    assert eps == pytest.approx(np.finfo(float).eps / 2.0, rel=1e-13)
    assert differentiation.calculate_step_size(2.0) == pytest.approx(math.sqrt(eps) * 3.0, rel=1e-15)
    assert differentiation.calculate_step_size(-1.0, 2) == pytest.approx(eps ** (1 / 3) * 2.0, rel=1e-15)


def test_gradient_jacobian_and_hessian_match_scipy_and_the_closed_form() -> None:
    # Arrange
    rosen = optimize.rosen
    theta = np.array([0.8, 1.3, -0.4])
    # Act
    grad = differentiation.gradient(rosen, theta)
    hess = differentiation.hessian(rosen, theta)
    jac = differentiation.jacobian(lambda t: np.array([t[0] * t[1], math.sin(t[2]), t[0] ** 2]), theta)
    at = differentiation.jacobian_at(lambda x, t: t[0] * math.exp(t[1] * x), [0.0, 1.0, 2.0], [1.5, 0.5])
    bounded = differentiation.gradient(lambda t: float(t[0] ** 3), [1.0], lower_bounds=[0.0], upper_bounds=[1.0])
    # Assert
    assert grad == pytest.approx(optimize.rosen_der(theta), rel=1e-7)
    assert hess == pytest.approx(optimize.rosen_hess(theta), rel=1e-5, abs=1e-5)
    assert np.array_equal(hess, hess.T)
    assert jac == pytest.approx(np.array([[1.3, 0.8, 0.0], [0.0, 0.0, math.cos(-0.4)], [1.6, 0.0, 0.0]]),
                                abs=1e-8)
    assert at == pytest.approx(np.array([[1.0, 0.0], [math.exp(0.5), 1.5 * math.exp(0.5)],
                                         [math.exp(1.0), 3.0 * math.exp(1.0)]]), rel=1e-7)
    assert bounded == pytest.approx([3.0], rel=1e-4)


# ODE solvers.

def test_runge_kutta_matches_scipy_solve_ivp_and_the_closed_form() -> None:
    # Arrange: dy/dt = -2 t y, y(0) = 1, so y = exp(-t^2).
    def f(t: float, y: float) -> float:
        return -2.0 * t * y

    times = np.linspace(0.0, 2.0, 201)
    reference = integrate.solve_ivp(lambda t, y: -2.0 * t * y, (0.0, 2.0), [1.0], t_eval=times,
                                    rtol=1e-12, atol=1e-14).y[0]
    # Act
    rk2 = ode.runge_kutta_second_order(f, 1.0, 0.0, 2.0, 201)
    rk4 = ode.runge_kutta_fourth_order(f, 1.0, 0.0, 2.0, 201)
    step = ode.runge_kutta_fourth_order_step(f, 1.0, 0.0, 0.01)
    fehlberg = ode.runge_kutta_fehlberg(f, 1.0, 0.0, 2.0, 1e-6, 1e-10)
    cash_karp = ode.runge_kutta_cash_karp(f, 1.0, 0.0, 2.0, 1e-6, 1e-10)
    # Assert
    assert rk2 == pytest.approx(reference, abs=1e-4)
    assert rk4 == pytest.approx(reference, abs=1e-9)
    assert rk4 == pytest.approx(np.exp(-times ** 2), abs=1e-9)
    assert step == rk4[1]
    assert [fehlberg, cash_karp] == pytest.approx([math.exp(-4.0)] * 2, abs=1e-8)


# Defaults are the Numerics defaults.

@pytest.mark.parametrize(("fn", "type_name", "member", "count"), [
    (roots.bisection_solve, "Numerics.Mathematics.RootFinding.Bisection", "Solve", 7),
    (roots.secant_solve, "Numerics.Mathematics.RootFinding.Secant", "Solve", 6),
    (roots.newton_raphson_solve, "Numerics.Mathematics.RootFinding.NewtonRaphson", "Solve", 6),
    (roots.newton_raphson_robust_solve, "Numerics.Mathematics.RootFinding.NewtonRaphson", "RobustSolve", 8),
    (ode.runge_kutta_fehlberg, "Numerics.Mathematics.ODESolvers.RungeKutta", "Fehlberg", 6),
    (ode.runge_kutta_cash_karp, "Numerics.Mathematics.ODESolvers.RungeKutta", "CashKarp", 6),
    (differentiation.gradient, "Numerics.Mathematics.NumericalDerivative", "Gradient", 7),
    (differentiation.hessian, "Numerics.Mathematics.NumericalDerivative", "Hessian", 7),
    (differentiation.jacobian, "Numerics.Mathematics.NumericalDerivative", "Jacobian", 7),
    (differentiation.calculate_step_size, "Numerics.Mathematics.NumericalDerivative", "CalculateStepSize", 2),
], ids=lambda value: value if isinstance(value, str) else None)
def test_python_defaults_are_the_numerics_defaults(fn, type_name, member, count) -> None:
    # Arrange
    import inspect  # noqa: PLC0415

    namespace("Numerics")
    import System  # noqa: PLC0415

    assembly = next(a for a in System.AppDomain.CurrentDomain.GetAssemblies() if a.GetName().Name == "Numerics")
    method = next(m for m in assembly.GetType(type_name).GetMethods()
                  if m.Name == member and len(m.GetParameters()) == count)
    net_defaults = [p.DefaultValue for p in method.GetParameters() if p.HasDefaultValue]
    # Act
    py_defaults = [p.default for p in inspect.signature(fn).parameters.values()
                   if p.default is not inspect.Parameter.empty]
    # Assert
    assert len(py_defaults) == len(net_defaults)
    for py, net in zip(py_defaults, net_defaults, strict=True):
        assert (py is None and net is None) or py == net


# Speed: a compiled integrand avoids the per-evaluation crossing into Python.

@pytest.mark.benchmark
@needs_roslyn
def test_compiled_integrand_is_faster_than_a_python_callback() -> None:
    # Arrange
    def run(f: object) -> float:
        start = time.perf_counter()
        for _ in range(5):
            integration.integrate_miser(f, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], max_function_evaluations=20000)
        return time.perf_counter() - start

    compiled = callbacks.compile_csharp("x => x[0] * x[1] * x[2]", "Func<double[],double>")
    run(compiled)
    # Act
    python_time = run(lambda x: float(x[0] * x[1] * x[2]))
    compiled_time = run(compiled)
    # Assert
    assert compiled_time < python_time


def test_ridders_default_step_outside_the_domain_is_halved() -> None:
    # Arrange: the default step 0.1 * 1.05 exceeds 0.05, so log(0.05 - h) is NaN for the first tableau
    def log(x: float) -> float:
        return math.log(x) if x > 0 else math.nan

    def sqrt(x: float) -> float:
        return math.sqrt(x) if x >= 0 else math.nan

    # Act
    value, error = differentiation.ridders(log, 0.05)
    # Assert
    assert value == pytest.approx(20.0, rel=1e-12)
    assert error < 1e-10
    assert differentiation.ridders(sqrt, 0.01)[0] == pytest.approx(5.0, rel=1e-12)
    assert differentiation.ridders("x => Math.Log(x)", 1e-6)[0] == pytest.approx(1e6, rel=1e-12)
    # A given step is used as it is
    assert differentiation.ridders(log, 0.05, step_size=0.1)[0] == pytest.approx(20.0, rel=1e-8)
