"""Optimizer wrappers against scipy.optimize and textbook optima, the shared result, callbacks and speed.

Each optimizer runs the Numerics search; the oracle is ``scipy.optimize`` (``minimize``,
``minimize_scalar``, ``differential_evolution``, ``linprog``) or the known optimum of the test
function, at a tolerance that covers the different stopping rules.
"""

from __future__ import annotations

import enum
import importlib
import inspect
import time

import numpy as np
import pytest
from numpydoc import validate as numpydoc_validate
from test_docstrings import ENFORCED_CHECKS

from pynumerics import callbacks, errors, optimization
from pynumerics.optimization import (
    ADAM,
    BFGS,
    CMAES,
    LBFGSB,
    MLSL,
    AugmentedLagrange,
    BrentSearch,
    Constraint,
    ConstraintType,
    DifferentialEvolution,
    DynamicallyDimensionedSearch,
    GoldenSection,
    GradientDescent,
    LocalMethod,
    MultiStart,
    NelderMead,
    OptimizationProgress,
    OptimizationStatus,
    ParameterSet,
    ParetoFront,
    ParticleSwarm,
    Powell,
    ShuffledComplexEvolution,
    SimulatedAnnealing,
)

opt = pytest.importorskip("scipy.optimize")
needs_roslyn = pytest.mark.skipif(not callbacks.roslyn_available(), reason="Roslyn assemblies not staged")

MODULES = ("pynumerics.optimization", "pynumerics.graphs")


def rosenbrock(x: np.ndarray) -> float:
    return float(100.0 * (x[1] - x[0] ** 2) ** 2 + (1.0 - x[0]) ** 2)


def rosenbrock_gradient(x: np.ndarray) -> np.ndarray:
    return np.array([-400.0 * x[0] * (x[1] - x[0] ** 2) - 2.0 * (1.0 - x[0]), 200.0 * (x[1] - x[0] ** 2)])


def mccormick(x: np.ndarray) -> float:
    return float(np.sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 1.5 * x[0] + 2.5 * x[1] + 1.0)


def ackley(x: np.ndarray) -> float:
    return float(-20.0 * np.exp(-0.2 * np.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2)))
                 - np.exp(0.5 * (np.cos(2 * np.pi * x[0]) + np.cos(2 * np.pi * x[1]))) + np.e + 20.0)


ACKLEY_CS = ("x => -20.0 * Math.Exp(-0.2 * Math.Sqrt(0.5 * (x[0] * x[0] + x[1] * x[1])))"
             " - Math.Exp(0.5 * (Math.Cos(2 * Math.PI * x[0]) + Math.Cos(2 * Math.PI * x[1]))) + Math.E + 20.0")

LOWER, UPPER = [-2.048, -2.048], [2.048, 2.048]


def _public_objects() -> list[tuple[str, str]]:
    """Return ``(module, dotted path)`` of every public function, class and class member of the modules."""
    names = []
    for module_name in MODULES:
        module = importlib.import_module(module_name)
        for name in module.__all__:
            obj = getattr(module, name)
            if not (inspect.isfunction(obj) or inspect.isclass(obj)) or (
                    inspect.isclass(obj) and issubclass(obj, enum.Enum)):
                continue
            names.append((module_name, name))
            if inspect.isclass(obj):
                names += [(module_name, f"{name}.{m}") for m, v in vars(obj).items() if not m.startswith("_")
                          and (inspect.isfunction(v) or isinstance(v, (property, classmethod, staticmethod)))]
    return names


@pytest.mark.parametrize(("module_name", "path"), _public_objects())
def test_public_docstrings_and_exception_translation(module_name: str, path: str) -> None:
    # Arrange
    owner, _, member = path.partition(".")
    obj = getattr(importlib.import_module(module_name), owner)
    target = vars(obj)[member] if member else obj
    # Act
    report = numpydoc_validate.validate(f"{module_name}.{path}")
    # Assert
    problems = [f"{code}: {message}" for code, message in report["errors"] if code in ENFORCED_CHECKS]
    assert not problems, f"{path}\n" + "\n".join(problems)
    if not inspect.isclass(target):
        assert errors.is_translated(target), f"{path} does not translate .NET exceptions"


@pytest.mark.parametrize("factory", [NelderMead, Powell, BFGS, LBFGSB])
def test_local_optimizers_match_scipy_minimize_on_rosenbrock(factory) -> None:
    # Arrange
    reference = opt.minimize(rosenbrock, [0.0, 0.0], jac=rosenbrock_gradient, method="L-BFGS-B",
                             bounds=list(zip(LOWER, UPPER, strict=True)), tol=1e-14)
    # Act
    result = factory(rosenbrock, [0.0, 0.0], LOWER, UPPER).minimize()
    # Assert
    assert result.status is OptimizationStatus.SUCCESS and result.success
    np.testing.assert_allclose(result.parameters, reference.x, atol=1e-6)
    assert result.value == pytest.approx(rosenbrock(result.parameters), abs=0.0)
    assert result.function_evaluations > 0 and result.iterations > 0
    assert result.hessian.shape == (2, 2)
    np.testing.assert_allclose(result.hessian, [[802.0, -400.0], [-400.0, 200.0]], rtol=1e-3)


@pytest.mark.parametrize("factory", [NelderMead, Powell, BFGS])
def test_mccormick_minimum_matches_scipy(factory) -> None:
    # Arrange
    reference = opt.minimize(mccormick, [0.0, 0.0], method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-14})
    # Act
    result = factory(mccormick, [0.0, 0.0], [-3.0, -3.0], [4.0, 4.0]).minimize()
    # Assert
    np.testing.assert_allclose(result.parameters, reference.x, atol=2e-4)
    assert result.value == pytest.approx(reference.fun, abs=1e-8)
    assert result.value == pytest.approx(-1.9132229549810362, abs=1e-8)


@pytest.mark.parametrize("factory", [BFGS, LBFGSB])
def test_analytic_gradient_is_used_and_gives_scipy_optimum(factory) -> None:
    # Arrange
    calls = []

    def gradient(x: np.ndarray) -> np.ndarray:
        calls.append(1)
        return rosenbrock_gradient(x)

    # Act
    optimizer = factory(rosenbrock, [-1.2, 1.0], LOWER, UPPER, gradient=gradient)
    result = optimizer.minimize()
    # Assert
    assert calls and optimizer.gradient is not None
    np.testing.assert_allclose(result.parameters, [1.0, 1.0], atol=1e-6)


def test_bound_constrained_minimum_matches_scipy_lbfgsb() -> None:
    # Arrange
    def f(x: np.ndarray) -> float:
        return float((x[0] - 3.0) ** 2 + (x[1] + 1.0) ** 2 + x[0] * x[1])

    bounds = [(-1.0, 2.0), (0.0, 5.0)]
    reference = opt.minimize(f, [0.5, 0.5], method="L-BFGS-B", bounds=bounds, tol=1e-14)
    # Act
    results = [factory(f, [0.5, 0.5], [-1.0, 0.0], [2.0, 5.0]).minimize() for factory in (LBFGSB, BFGS, NelderMead)]
    # Assert
    for result in results:
        np.testing.assert_allclose(result.parameters, reference.x, atol=1e-6)


@pytest.mark.parametrize("factory", [GradientDescent, ADAM])
def test_first_order_methods_reach_the_quadratic_minimum(factory) -> None:
    # Arrange
    centre = np.array([1.5, -0.5])

    def f(x: np.ndarray) -> float:
        return float(np.sum((x - centre) ** 2))

    # Act
    optimizer = factory(f, [0.0, 0.0], [-5.0, -5.0], [5.0, 5.0], alpha=0.05, gradient=lambda x: 2.0 * (x - centre))
    result = optimizer.minimize()
    # Assert
    assert optimizer.alpha == 0.05
    np.testing.assert_allclose(result.parameters, opt.minimize(f, [0.0, 0.0]).x, atol=1e-4)


def test_adam_moment_settings_round_trip() -> None:
    # Act
    adam = ADAM(lambda x: float(x[0] ** 2), [1.0], [-2.0], [2.0], beta1=0.8, beta2=0.99)
    # Assert
    assert (adam.beta1, adam.beta2, adam.alpha) == (0.8, 0.99, 0.001)


@pytest.mark.parametrize("factory", [GoldenSection, BrentSearch])
def test_scalar_searches_match_minimize_scalar(factory) -> None:
    # Arrange
    def f(x: float) -> float:
        return (x - 2.0) ** 4 + np.cos(3.0 * x)

    reference = opt.minimize_scalar(f, bounds=(0.0, 4.0), method="bounded", options={"xatol": 1e-12})
    # Act
    search = factory(f, 0.0, 4.0)
    result = search.minimize()
    # Assert
    assert (search.lower_bound, search.upper_bound) == (0.0, 4.0)
    assert result.parameters[0] == pytest.approx(reference.x, abs=1e-6)
    assert result.value == pytest.approx(reference.fun, abs=1e-12)


def test_brent_bracket_moves_the_interval_around_the_minimum() -> None:
    # Arrange
    search = BrentSearch(lambda x: (x - 10.0) ** 2, 0.0, 1.0)
    # Act
    lower, upper = search.bracket()
    result = search.minimize()
    # Assert
    assert lower < 10.0 < upper
    assert result.parameters[0] == pytest.approx(10.0, abs=1e-7)


@pytest.mark.parametrize(("f", "bounds"), [
    (lambda x: (x + 3.0) * ((x - 1.0) * (x - 1.0)), (-3.0, 3.0)),
    (lambda x: x * x * x * x - 3.0 * x * x * x + 2.0, (0.0, 5.0)),
    (lambda x: abs(x - 0.3), (-1.0, 1.0)),
    (lambda x: (x - 2.0) ** 4 + np.cos(3.0 * x), (0.0, 4.0)),
], ids=["cubic", "quartic", "kink", "wavy"])
def test_brent_search_is_bitwise_minimize_scalar_bounded_at_scipy_settings(f, bounds) -> None:
    # Arrange
    reference = opt.minimize_scalar(f, bounds=bounds, method="bounded")
    # Act
    search = BrentSearch(f, *bounds, relative_tolerance=np.sqrt(2.2e-16), absolute_tolerance=1e-5,
                         max_function_evaluations=500)
    result = search.minimize()
    # Assert
    assert np.float64(result.parameters[0]).view(np.uint64) == np.float64(reference.x).view(np.uint64)
    assert search.function_evaluations == reference.nfev


def test_differential_evolution_matches_scipy_on_ackley() -> None:
    # Arrange
    bounds = [(-5.0, 5.0), (-5.0, 5.0)]
    reference = opt.differential_evolution(ackley, bounds, seed=1, tol=1e-12, polish=True)
    # Act
    result = DifferentialEvolution(ackley, [-5.0, -5.0], [5.0, 5.0], population_size=40).minimize()
    # Assert
    np.testing.assert_allclose(result.parameters, reference.x, atol=1e-6)
    assert result.value == pytest.approx(reference.fun, abs=1e-6)


@pytest.mark.parametrize(("factory", "settings"), [
    (ParticleSwarm, {"population_size": 50}),
    (ShuffledComplexEvolution, {"complexes": 5}),
    (SimulatedAnnealing, {}),
    (CMAES, {}),
    (DifferentialEvolution, {"population_size": 50, "polish": True, "polish_method": "bfgs"}),
])
def test_global_optimizers_find_the_rosenbrock_minimum(factory, settings) -> None:
    # Arrange
    reference = opt.differential_evolution(rosenbrock, list(zip(LOWER, UPPER, strict=True)), seed=3, tol=1e-12)
    # Act
    optimizer = factory(rosenbrock, LOWER, UPPER, **settings)
    result = optimizer.minimize()
    # Assert
    np.testing.assert_allclose(result.parameters, reference.x, atol=1e-3)
    for name, value in settings.items():
        expected = LocalMethod(value) if name == "polish_method" else value
        assert getattr(optimizer, name) == expected


def test_dynamically_dimensioned_search_reaches_the_sphere_minimum() -> None:
    # Act
    dds = DynamicallyDimensionedSearch(lambda x: float(np.sum(x ** 2)), [-5.0] * 4, [5.0] * 4,
                                       initial_values=[3.0, -2.0, 1.0, 4.0], evaluation_budget=20000)
    result = dds.minimize()
    # Assert
    np.testing.assert_allclose(dds.initial_values, [3.0, -2.0, 1.0, 4.0])
    assert dds.evaluation_budget == 20000
    assert result.value < 1e-2


@pytest.mark.parametrize("factory", [MultiStart, MLSL])
def test_multi_start_methods_find_the_three_hump_camel_minimum(factory) -> None:
    # Arrange
    def camel(x: np.ndarray) -> float:
        return float(2 * x[0] ** 2 - 1.05 * x[0] ** 4 + x[0] ** 6 / 6 + x[0] * x[1] + x[1] ** 2)

    reference = opt.differential_evolution(camel, [(-5.0, 5.0), (-5.0, 5.0)], seed=5, tol=1e-12)
    # Act
    optimizer = factory(camel, [3.0, 3.0], [-5.0, -5.0], [5.0, 5.0], method=LocalMethod.NELDER_MEAD, max_iterations=30)
    result = optimizer.minimize()
    # Assert
    assert optimizer.method is LocalMethod.NELDER_MEAD
    np.testing.assert_allclose(result.parameters, reference.x, atol=1e-5)
    if factory is MLSL:
        assert optimizer.sampled_points and optimizer.local_minimums
        assert any(p.minimized for p in optimizer.sampled_points)


def test_maximize_reports_the_objective_with_its_sign() -> None:
    # Arrange
    def bump(x: np.ndarray) -> float:
        return float(5.0 - (x[0] - 1.0) ** 2 - (x[1] + 0.5) ** 2)

    optimizer = NelderMead(bump, [0.0, 0.0], [-3.0, -3.0], [3.0, 3.0])
    # Act
    result = optimizer.maximize()
    # Assert
    assert result.maximized and result.value == pytest.approx(5.0, abs=1e-7)
    assert optimizer.best_parameter_set.fitness == -result.value
    np.testing.assert_allclose(result.parameters, [1.0, -0.5], atol=1e-3)
    assert result.trace_values[-1] == result.value and np.all(np.diff(result.trace_values) >= 0.0)


def test_trace_rows_are_the_running_best_of_every_evaluation() -> None:
    # Act
    optimizer = NelderMead(rosenbrock, [0.0, 0.0], LOWER, UPPER)
    result = optimizer.minimize()
    trace = optimizer.parameter_set_trace
    # Assert
    assert result.trace_parameters.shape == (result.function_evaluations, 2)
    assert np.all(np.diff(result.trace_values) <= 0.0)
    np.testing.assert_array_equal(result.trace_parameters[-1], result.parameters)
    assert [p.fitness for p in trace] == list(result.trace_values)


@needs_roslyn
def test_bulk_trace_copy_gives_the_bits_of_the_entry_by_entry_copy() -> None:
    # Arrange
    optimizer = DifferentialEvolution(ackley, [-5.0, -5.0], [5.0, 5.0])
    optimizer.minimize()
    trace = optimizer.net.ParameterSetTrace
    # Act
    bulk = optimization._trace_arrays(trace, 1.0, bulk=True)
    single = optimization._trace_arrays(trace, 1.0, bulk=False)
    # Assert
    for a, b in zip(bulk, single, strict=True):
        np.testing.assert_array_equal(a.view(np.uint64), b.view(np.uint64))


def test_record_traces_off_gives_an_empty_trace() -> None:
    # Act
    result = NelderMead(rosenbrock, [0.0, 0.0], LOWER, UPPER, record_traces=False, compute_hessian=False).minimize()
    # Assert
    assert result.trace_values.size == 0 and result.hessian is None


def test_settings_round_trip_and_unknown_settings_raise() -> None:
    # Act
    nm = NelderMead(rosenbrock, [0.0, 0.0], LOWER, UPPER, max_iterations=500, relative_tolerance=1e-10,
                    absolute_tolerance=1e-10, simplex_scale=0.1, enable_start_point_probe=True,
                    start_point_probe_count=4, start_point_probe_radius=0.01, progress_interval=7,
                    max_function_evaluations=100000, report_failure=False)
    # Assert
    assert (nm.max_iterations, nm.relative_tolerance, nm.simplex_scale) == (500, 1e-10, 0.1)
    assert nm.start_point_probe_count == 4
    assert nm.enable_start_point_probe and not nm.report_failure and nm.progress_interval == 7
    assert nm.number_of_parameters == 2 and nm.max_function_evaluations == 100000
    np.testing.assert_array_equal(nm.lower_bounds, LOWER)
    np.testing.assert_array_equal(nm.initial_values, [0.0, 0.0])
    with pytest.raises(TypeError, match="no setting 'bogus'"):
        NelderMead(rosenbrock, [0.0, 0.0], LOWER, UPPER, bogus=1)
    with pytest.raises(TypeError, match="no setting 'number_of_parameters'"):
        NelderMead(rosenbrock, [0.0, 0.0], LOWER, UPPER, number_of_parameters=3)


def test_iteration_limit_sets_the_status_and_invalid_settings_raise() -> None:
    # Act
    limited = NelderMead(rosenbrock, [0.0, 0.0], LOWER, UPPER, max_iterations=10).minimize()
    # Assert
    assert limited.status is OptimizationStatus.MAXIMUM_ITERATIONS_REACHED and not limited.success
    with pytest.raises(ValueError, match="greater than 10"):
        NelderMead(rosenbrock, [0.0, 0.0], LOWER, UPPER, max_iterations=5).minimize()


def test_objective_errors_propagate_and_clear_results_resets() -> None:
    # Arrange
    def broken(x: np.ndarray) -> float:
        raise ZeroDivisionError("objective failed")

    optimizer = NelderMead(rosenbrock, [0.0, 0.0], LOWER, UPPER)
    optimizer.minimize()
    # Act
    optimizer.clear_results()
    # Assert
    assert optimizer.status is OptimizationStatus.NONE and optimizer.function_evaluations == 0
    optimizer.objective_function = broken
    with pytest.raises(Exception, match="objective failed"):
        optimizer.minimize()


def test_progress_reports_reach_the_callable_and_cancellation_stops_the_search() -> None:
    # Arrange
    import System  # noqa: PLC0415 - the optimizer started the runtime

    reports: list[OptimizationProgress] = []
    source = System.Threading.CancellationTokenSource()

    def on_progress(report: OptimizationProgress) -> None:
        reports.append(report)
        if len(reports) == 3:
            source.Cancel()

    optimizer = DifferentialEvolution(ackley, [-5.0, -5.0], [5.0, 5.0], progress_interval=50)
    optimizer.progress = on_progress
    optimizer.cancellation_token = source.Token
    # Act
    result = optimizer.minimize()
    # Assert
    assert optimizer.progress is on_progress and optimizer.cancellation_token.IsCancellationRequested
    assert result.status is OptimizationStatus.CANCELLED
    assert [r.function_evaluations for r in reports[:3]] == [50, 100, 150]
    assert reports[0].best_values.shape == (2,)
    round_trip = OptimizationProgress.from_net(reports[0].to_net())
    assert round_trip.best_fitness == reports[0].best_fitness


def test_parameter_set_xml_and_clone_round_trip() -> None:
    # Arrange
    original = ParameterSet(np.array([0.1, -2.5, 1e-300]), 3.25, 0.5)
    # Act
    xml = original.to_xml()
    parsed = ParameterSet.from_xml(xml)
    copy = original.clone()
    # Assert
    assert xml.startswith("<ParameterSet")
    np.testing.assert_array_equal(parsed.values, original.values)
    assert (parsed.fitness, parsed.weight) == (3.25, 0.5)
    assert copy.values is not original.values and np.array_equal(copy.values, original.values)
    assert ParameterSet.from_net(original.to_net()).fitness == 3.25


@pytest.mark.parametrize(("relation", "rhs", "scipy_constraint"), [
    (ConstraintType.LESSER_THAN_OR_EQUAL_TO, 2.0, {"type": "ineq", "fun": lambda x: 2.0 - x[0] ** 2 - x[1] ** 2}),
    (ConstraintType.EQUAL_TO, 1.0, {"type": "eq", "fun": lambda x: x[0] ** 2 + x[1] ** 2 - 1.0}),
    ("greater_than_or_equal_to", 10.0, {"type": "ineq", "fun": lambda x: x[0] ** 2 + x[1] ** 2 - 10.0}),
])
def test_augmented_lagrange_matches_scipy_slsqp(relation, rhs, scipy_constraint) -> None:
    # Arrange
    def objective(x: np.ndarray) -> float:
        return float((x[0] - 2.0) ** 2 + (x[1] - 1.0) ** 2)

    reference = opt.minimize(objective, [0.5, 0.5], method="SLSQP", constraints=[scipy_constraint],
                             bounds=[(-5.0, 5.0)] * 2, options={"ftol": 1e-14})
    inner = BFGS(objective, [0.5, 0.5], [-5.0, -5.0], [5.0, 5.0])
    constraint = Constraint(lambda x: float(x[0] ** 2 + x[1] ** 2), 2, rhs, relation)
    solver = AugmentedLagrange(objective, inner, [constraint])
    # Act
    result = solver.minimize()
    # Assert
    np.testing.assert_allclose(result.parameters, reference.x, atol=1e-4)
    assert result.value == pytest.approx(reference.fun, abs=1e-4)
    assert solver.optimizer is inner and solver.constraints == (constraint,)
    assert constraint.type is ConstraintType(relation) and constraint.value == rhs
    assert constraint.number_of_parameters == 2 and constraint.tolerance == 1e-8 and constraint.function is not None
    multipliers = {ConstraintType.EQUAL_TO: solver.lambda_, ConstraintType.LESSER_THAN_OR_EQUAL_TO: solver.mu,
                   ConstraintType.GREATER_THAN_OR_EQUAL_TO: solver.nu}
    assert multipliers[ConstraintType(relation)].shape == (1,)


def test_augmented_lagrange_solves_a_linear_program_like_linprog() -> None:
    # Arrange: minimize -x - 2y subject to x + y <= 4, x - y >= -2, 0 <= x, y <= 3.
    reference = opt.linprog([-1.0, -2.0], A_ub=[[1.0, 1.0], [-1.0, 1.0]], b_ub=[4.0, 2.0], bounds=[(0, 3), (0, 3)])

    def objective(x: np.ndarray) -> float:
        return float(-x[0] - 2.0 * x[1])

    constraints = [Constraint(lambda x: float(x[0] + x[1]), 2, 4.0, "LesserThanOrEqualTo"),
                   Constraint(lambda x: float(x[0] - x[1]), 2, -2.0, "GreaterThanOrEqualTo")]
    inner = NelderMead(objective, [0.5, 0.5], [0.0, 0.0], [3.0, 3.0])
    # Act
    result = AugmentedLagrange(objective, inner, constraints).minimize()
    # Assert
    np.testing.assert_allclose(result.parameters, reference.x, atol=1e-3)
    assert result.value == pytest.approx(reference.fun, abs=1e-3)


def test_augmented_lagrange_needs_a_constraint() -> None:
    # Act / Assert
    inner = BFGS(rosenbrock, [0.0, 0.0], LOWER, UPPER)
    with pytest.raises(ValueError, match="at least one constraint"):
        AugmentedLagrange(rosenbrock, inner, [])


@needs_roslyn
def test_parallel_policy_keeps_python_objectives_serial_and_compiled_ones_parallel() -> None:
    # Act
    python = DifferentialEvolution(ackley, [-5.0, -5.0], [5.0, 5.0], use_parallel=True)
    compiled = DifferentialEvolution(ACKLEY_CS, [-5.0, -5.0], [5.0, 5.0], use_parallel=True,
                                     max_degree_of_parallelism=4)
    swarm = ParticleSwarm(ackley, [-5.0, -5.0], [5.0, 5.0])
    swarm.use_parallel = True
    # Assert
    assert not python.use_parallel and not swarm.use_parallel
    assert compiled.use_parallel and compiled.max_degree_of_parallelism == 4
    parallel = compiled.minimize()
    np.testing.assert_allclose(parallel.parameters, [0.0, 0.0], atol=1e-6)


@needs_roslyn
def test_augmented_lagrange_applies_the_policy_to_python_constraints() -> None:
    # Arrange
    inner = DifferentialEvolution("x => (x[0] - 2) * (x[0] - 2) + (x[1] - 1) * (x[1] - 1)", [-5.0, -5.0],
                                  [5.0, 5.0], use_parallel=True)
    constraint = Constraint(lambda x: float(x[0] + x[1]), 2, 1.0, "EqualTo")
    # Act
    AugmentedLagrange("x => (x[0] - 2) * (x[0] - 2) + (x[1] - 1) * (x[1] - 1)", inner, [constraint])
    # Assert
    assert not inner.use_parallel


def _run_ackley_differential_evolution(objective: object) -> tuple[np.ndarray, float, float]:
    """Minimize Ackley with DifferentialEvolution three times; return the parameters, value and best seconds."""
    best = np.inf
    for _ in range(3):
        optimizer = DifferentialEvolution(objective, [-5.0, -5.0], [5.0, 5.0], population_size=40,
                                          record_traces=False, compute_hessian=False)
        start = time.perf_counter()
        result = optimizer.minimize()
        best = min(best, time.perf_counter() - start)
    return result.parameters, result.value, best


def _ackley_math(x: np.ndarray) -> float:
    """Return the two-dimensional Ackley function with the math module, as ACKLEY_CS computes it."""
    import math  # noqa: PLC0415

    return (-20.0 * math.exp(-0.2 * math.sqrt(0.5 * (x[0] * x[0] + x[1] * x[1])))
            - math.exp(0.5 * (math.cos(2 * math.pi * x[0]) + math.cos(2 * math.pi * x[1]))) + math.e + 20.0)


@needs_roslyn
def test_differential_evolution_compiled_objective_gives_the_python_bits() -> None:
    # Act
    python_x, python_value, _ = _run_ackley_differential_evolution(_ackley_math)
    compiled_x, compiled_value, _ = _run_ackley_differential_evolution(ACKLEY_CS)
    # Assert
    np.testing.assert_allclose(compiled_x, python_x, atol=1e-12)
    assert compiled_value == pytest.approx(python_value, abs=1e-14)


@pytest.mark.benchmark
@needs_roslyn
def test_differential_evolution_compiled_objective_is_faster_than_python() -> None:
    # Act
    _, _, python_seconds = _run_ackley_differential_evolution(_ackley_math)
    _, _, compiled_seconds = _run_ackley_differential_evolution(ACKLEY_CS)
    # Assert
    assert compiled_seconds * 3 < python_seconds, (compiled_seconds, python_seconds)


def test_differential_evolution_initial_population_and_integer_parameters() -> None:
    # Arrange
    population = np.random.default_rng(2).uniform(-5.0, 5.0, size=(30, 2))

    def f(x: np.ndarray) -> float:
        return float((x[0] - 2.3) ** 2 + (x[1] + 1.0) ** 2)

    # Act
    de = DifferentialEvolution(f, [-5.0, -5.0], [5.0, 5.0], population_size=40, mutation=0.6, dither_rate=0.5,
                               crossover_probability=0.8)
    de.initial_population = population
    de.integer_parameters = [True, False]
    result = de.minimize()
    # Assert
    np.testing.assert_array_equal(de.initial_population, population)
    np.testing.assert_array_equal(de.integer_parameters, [True, False])
    assert (de.mutation, de.dither_rate, de.crossover_probability) == (0.6, 0.5, 0.8)
    np.testing.assert_allclose(result.parameters, [2.0, -1.0], atol=1e-6)
    de.initial_population = None
    de.integer_parameters = None
    assert de.initial_population is None and de.integer_parameters is None


def test_global_settings_round_trip() -> None:
    # Act
    sa = SimulatedAnnealing(rosenbrock, LOWER, UPPER, initial_temperature=5.0, min_temperature=0.2, cooling_rate=0.9,
                            update_cycles=3, temperature_cycles=8, tolerance_steps=10, prng_seed=7)
    sce = ShuffledComplexEvolution(rosenbrock, LOWER, UPPER, cce_iterations=4, tolerance_steps=9)
    cmaes = CMAES(rosenbrock, LOWER, UPPER, initial_values=[0.5, 0.5], population_size=12, initial_step_size=0.2,
                  restarts=2)
    mlsl = MLSL(rosenbrock, [0.0, 0.0], LOWER, UPPER, sample_size=10, gamma=0.1, sigma=3.0, min_no_improvement=2,
                max_no_improvement=4, local_absolute_tolerance=1e-7, local_relative_tolerance=1e-7, polish=False,
                prng_seed=3)
    bfgs = BFGS(rosenbrock, [0.0, 0.0], LOWER, UPPER, scale_initial_inverse_hessian=True)
    lbfgsb = LBFGSB(rosenbrock, [0.0, 0.0], LOWER, UPPER, memory_size=5, max_line_search_steps=30)
    cmaes.minimize()
    # Assert
    assert (sa.initial_temperature, sa.min_temperature, sa.cooling_rate) == (5.0, 0.2, 0.9)
    assert (sa.update_cycles, sa.temperature_cycles, sa.tolerance_steps, sa.prng_seed) == (3, 8, 10, 7)
    assert (sce.cce_iterations, sce.tolerance_steps) == (4, 9)
    assert (cmaes.population_size, cmaes.initial_step_size, cmaes.restarts) == (12, 0.2, 2)
    assert 0 <= cmaes.restarts_performed <= 2
    np.testing.assert_array_equal(cmaes.initial_values, [0.5, 0.5])
    assert (mlsl.sample_size, mlsl.gamma, mlsl.sigma, mlsl.min_no_improvement, mlsl.max_no_improvement) == (
        10, 0.1, 3.0, 2, 4)
    assert (mlsl.local_absolute_tolerance, mlsl.local_relative_tolerance, mlsl.polish) == (1e-7, 1e-7, False)
    assert mlsl.prng_seed == 3 and mlsl.sampled_points == [] and mlsl.local_minimums == []
    assert bfgs.scale_initial_inverse_hessian and bfgs.gradient is None
    assert (lbfgsb.memory_size, lbfgsb.max_line_search_steps) == (5, 30)


def _dominates(a: np.ndarray, b: np.ndarray) -> bool:
    return bool(np.all(a <= b) and np.any(a < b))


def _fronts(points: np.ndarray) -> list[list[int]]:
    """Textbook front numbering: one plus the largest front among the dominators."""
    order = np.lexsort(points.T[::-1])
    rank = np.full(len(points), -1)
    for i in order:
        rank[i] = 1 + max([rank[j] for j in range(len(points)) if rank[j] >= 0 and _dominates(points[j], points[i])],
                          default=-1)
    return [sorted(np.flatnonzero(rank == k).tolist()) for k in range(rank.max() + 1)]


def _hypervolume(points: np.ndarray, reference: np.ndarray) -> float:
    """Exact hypervolume by summing the cells of the coordinate grid that some point dominates."""
    inside = points[np.all(points < reference, axis=1)]
    grids = [np.unique(np.append(inside[:, k], reference[k])) for k in range(points.shape[1])]
    mesh = np.meshgrid(*[g[:-1] for g in grids], indexing="ij")
    widths = np.meshgrid(*[np.diff(g) for g in grids], indexing="ij")
    corners = np.stack([m.ravel() for m in mesh], axis=1)
    covered = np.any(np.all(inside[None, :, :] <= corners[:, None, :], axis=2), axis=1)
    return float(np.sum(np.prod(np.stack([w.ravel() for w in widths], axis=1), axis=1)[covered]))


def _crowding(front: np.ndarray) -> np.ndarray:
    n, m = front.shape
    distance = np.zeros(n)
    if n <= 2:
        return np.full(n, np.inf)
    for k in range(m):
        order = np.argsort(front[:, k], kind="stable")
        span = front[order[-1], k] - front[order[0], k]
        if span == 0:
            continue
        distance[order[0]] = distance[order[-1]] = np.inf
        distance[order[1:-1]] += (front[order[2:], k] - front[order[:-2], k]) / span
    return distance


@pytest.mark.parametrize("objectives", [2, 3, 4])
def test_pareto_sorting_and_hypervolume_match_textbook_definitions(objectives) -> None:
    # Arrange
    points = np.round(np.random.default_rng(objectives).random((40, objectives)), 2)
    reference = np.full(objectives, 1.1)
    # Act
    fronts = ParetoFront.non_dominated_sort(points)
    first = ParetoFront.non_dominated_indices(points)
    volume = ParetoFront.hypervolume(points, reference)
    # Assert
    assert [f.tolist() for f in fronts] == _fronts(points)
    assert first.tolist() == _fronts(points)[0]
    assert volume == pytest.approx(_hypervolume(points, reference), rel=1e-12)


def test_crowding_and_generational_distances_match_numpy() -> None:
    # Arrange
    rng = np.random.default_rng(9)
    population = rng.random((30, 2))
    front = ParetoFront.non_dominated_sort(population)[0]
    reference = np.column_stack([np.linspace(0, 1, 50), 1 - np.linspace(0, 1, 50)])
    spatial = pytest.importorskip("scipy.spatial")
    distances = spatial.distance.cdist(population, reference)
    # Act / Assert
    np.testing.assert_allclose(ParetoFront.crowding_distance(population[front]), _crowding(population[front]))
    np.testing.assert_allclose(ParetoFront.crowding_distance(population, front), _crowding(population[front]))
    assert ParetoFront.generational_distance(population, reference) == pytest.approx(distances.min(1).mean(), rel=1e-14)
    assert ParetoFront.inverted_generational_distance(population, reference) == pytest.approx(
        distances.min(0).mean(), rel=1e-14)
    assert ParetoFront.dominates([1.0, 2.0], [1.0, 3.0]) and not ParetoFront.dominates([1.0, 2.0], [1.0, 2.0])
    with pytest.raises(ValueError):
        ParetoFront.dominates([1.0, np.nan], [1.0, 2.0])
