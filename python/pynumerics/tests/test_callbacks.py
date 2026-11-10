"""Python callables and Roslyn-compiled C# as Numerics delegates, and the parallel rule."""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
from scipy import integrate

from pynumerics import callbacks, roots
from pynumerics._util import namespace
from pynumerics.marshal import to_net_doubles, to_net_ints

needs_roslyn = pytest.mark.skipif(not callbacks.roslyn_available(), reason="Roslyn assemblies not staged")


def _integrate(delegate: object) -> tuple[float, int, float]:
    """Integrate exp(-x^2) on [-10, 10] with AdaptiveGaussKronrod; return result, evaluations, seconds each."""
    integration = namespace("Numerics.Mathematics.Integration")
    best, result, count = math.inf, math.nan, 0
    for _ in range(15):
        integrator = integration.AdaptiveGaussKronrod(delegate, -10.0, 10.0)
        start = time.perf_counter()
        integrator.Integrate()
        elapsed = time.perf_counter() - start
        result, count = float(integrator.Result), int(integrator.FunctionEvaluations)
        best = min(best, elapsed / count)
    return result, count, best


def test_scalar_python_callback_integrates_like_scipy() -> None:
    # Act
    result, count, _ = _integrate(callbacks.scalar_function(lambda x: math.exp(-x * x)))
    # Assert
    expected, _ = integrate.quad(lambda x: math.exp(-x * x), -10.0, 10.0, epsabs=1e-14)
    assert count > 0
    assert result == pytest.approx(math.sqrt(math.pi), rel=1e-12)
    assert result == pytest.approx(expected, rel=1e-12)


@needs_roslyn
def test_compiled_lambda_gives_the_python_callback_bits_and_runs_faster() -> None:
    # Arrange
    python = callbacks.scalar_function(lambda x: math.exp(-x * x))
    compiled = callbacks.compile_csharp("x => Math.Exp(-x * x)")
    # Act
    python_result, python_count, python_cost = _integrate(python)
    compiled_result, compiled_count, compiled_cost = _integrate(compiled)
    # Assert
    assert compiled_result == python_result and compiled_count == python_count
    assert compiled_cost * 5.0 < python_cost, (compiled_cost, python_cost)


@needs_roslyn
def test_compile_csharp_caches_by_source_and_type() -> None:
    # Act
    first = callbacks.compile_csharp("x => x * x + 1.0")
    second = callbacks.compile_csharp("x => x * x + 1.0", "Func<double, double>")
    other = callbacks.compile_csharp("(x, y) => x * x + 1.0", "Func<double,double,double>")
    # Assert
    assert first is second and other is not first
    assert first(3.0) == 10.0 and other(3.0, 0.0) == 10.0


@needs_roslyn
def test_compile_csharp_accepts_static_methods_and_numerics_types() -> None:
    # Act
    total = callbacks.compile_csharp("static double Total(double[] x) { return x.Sum(); }", "Func<double[],double>")
    likelihood = callbacks.compile_csharp("p => new Normal(0, 1).LogPDF(p[0])", "LogLikelihood")
    statistic = callbacks.compile_csharp("x => Statistics.Mean(x)", "Func<double[],double>",
                                         usings=("Numerics.Data.Statistics",))
    # Assert
    assert total(to_net_doubles([1.0, 2.0, 3.5])) == 6.5
    assert likelihood(to_net_doubles([0.0])) == pytest.approx(-0.5 * math.log(2.0 * math.pi), rel=1e-15)
    assert statistic(to_net_doubles([1.0, 2.0, 6.0])) == 3.0


@needs_roslyn
def test_compile_errors_raise_value_error_with_roslyn_diagnostics() -> None:
    # Act / Assert
    with pytest.raises(ValueError, match=r"\(1,6\): error CS0103: The name 'y' does not exist"):
        callbacks.compile_csharp("x => y * x")
    with pytest.raises(ValueError, match="error CS"):
        callbacks.compile_csharp("x => x[0]")
    with pytest.raises(ValueError, match="exactly one static method"):
        callbacks.compile_csharp("static double A(double x) => x; static double B(double x) => -x;")


@needs_roslyn
def test_wrappers_accept_csharp_source_wherever_they_take_a_callable() -> None:
    # Act
    compiled = roots.brent_solve("x => x * x - 2.0", 0.0, 2.0)
    python = roots.brent_solve(lambda x: x * x - 2.0, 0.0, 2.0)
    # Assert
    assert compiled == python == pytest.approx(math.sqrt(2.0), abs=1e-8)


def test_to_delegate_converts_arrays_vectors_and_tuples() -> None:
    # Arrange
    linear = namespace("Numerics.Mathematics.LinearAlgebra")
    seen: list[object] = []

    def objective(x: np.ndarray) -> float:
        seen.append(type(x))
        return float(np.sum(x * x))

    # Act
    vector_value = callbacks.vector_function(objective)(to_net_doubles([1.0, 2.0]))
    mapped = callbacks.to_delegate(lambda x: x * 2.0, "Func<double[],double[]>")(to_net_doubles([1.0, 2.0]))
    vector = callbacks.to_delegate(lambda v: v[::-1], "Func<Vector,Vector>")(linear.Vector(to_net_doubles([1.0, 3.0])))
    pair = callbacks.to_delegate(lambda p: (-0.5 * float(p @ p), -p), "NUTS.LogLikelihoodAndGradient")(
        to_net_doubles([1.0, 2.0]))
    two = callbacks.bivariate_function(lambda t, y: t * y)(2.0, 3.5)
    likelihood = callbacks.log_likelihood(lambda p: -float(p[0]) ** 2)(to_net_doubles([3.0]))
    # Assert
    assert vector_value == 5.0 and seen == [np.ndarray]
    assert list(mapped) == [2.0, 4.0]
    assert list(vector.Array) == [3.0, 1.0]
    assert pair.Item1 == -2.5 and list(pair.Item2.Array) == [-1.0, -2.0]
    assert two == 7.0 and likelihood == -9.0


def test_transition_distribution_callbacks_fill_the_net_arrays() -> None:
    # Arrange
    def transition(stage: int, state: int, decision: int, next_states: object, probabilities: object) -> int:
        next_states[0] = state + decision
        probabilities[0] = 1.0
        return 1

    delegate = callbacks.to_delegate(transition, "TransitionDistribution")
    states, probabilities = to_net_ints([0, 0]), to_net_doubles([0.0, 0.0])
    # Act
    count = delegate(0, 2, 3, states, probabilities)
    # Assert
    assert count == 1 and states[0] == 5 and probabilities[0] == 1.0


def test_delegates_pass_through_and_types_are_checked() -> None:
    # Arrange
    delegate = callbacks.scalar_function(math.sin)
    # Act / Assert
    assert callbacks.scalar_function(delegate) is delegate
    with pytest.raises(TypeError, match="is not a"):
        callbacks.vector_function(delegate)
    with pytest.raises(TypeError, match="Expected a callable"):
        callbacks.scalar_function(3.0)
    with pytest.raises(ValueError, match="Unknown delegate type"):
        callbacks.delegate_type("Func<string>")
    with pytest.raises(TypeError, match="not a .NET delegate type"):
        callbacks.delegate_type(namespace("System").Double)
    assert "LogLikelihood" in callbacks.delegate_types()
    for name in callbacks.delegate_types():
        assert callbacks.delegate_type(name) is callbacks.delegate_type(name.replace(",", ", "))


@needs_roslyn
def test_parallel_policy_turns_parallel_off_only_for_python_callbacks() -> None:
    # Arrange
    optimization = namespace("Numerics.Mathematics.Optimization")
    compiled = callbacks.compile_csharp("x => x[0] * x[0]", "Func<double[],double>")
    python = callbacks.vector_function(lambda x: float(x[0] ** 2))
    bounds = to_net_doubles([-1.0]), to_net_doubles([1.0])
    native = optimization.DifferentialEvolution(compiled, 1, *bounds)
    hosted = optimization.DifferentialEvolution(python, 1, *bounds)
    native.UseParallel = hosted.UseParallel = True
    # Act
    native_changed = callbacks.apply_parallel_policy(native, compiled)
    hosted_changed = callbacks.apply_parallel_policy(hosted, python)
    # Assert
    assert callbacks.uses_python_callback(lambda x: x) and callbacks.uses_python_callback(python)
    assert not callbacks.uses_python_callback(compiled, "x => x", None)
    assert not native_changed and bool(native.UseParallel)
    assert hosted_changed and not bool(hosted.UseParallel)


_PATTERNS = np.array([0x7FF8000000000001, 0xFFF8DEADBEEF0001, 0x8000000000000000, 0x0000000000000001,
                      0x7FF0000000000000, 0xFFF0000000000000, 0x3FB999999999999A], dtype=np.uint64)


def _received(delegate_type: str, wrap: object, lengths: tuple[int, ...]) -> list[tuple[np.ndarray, np.ndarray]]:
    """Call a recording callback with arrays of each length; return (sent bits, received array) pairs."""
    seen: list[np.ndarray] = []
    delegate = callbacks.to_delegate(lambda x: seen.append(x) or 0.0, delegate_type)
    sent = []
    for length in lengths:
        bits = np.resize(_PATTERNS, length)
        sent.append(bits)
        delegate(wrap(to_net_doubles(bits.view(np.float64))))
    return list(zip(sent, seen, strict=True))


@pytest.mark.parametrize("delegate_type, wrap", [
    ("Func<double[],double>", lambda net: net),
    ("LogLikelihood", lambda net: net),
    ("Func<IList<double>,double>", lambda net: net),
    ("Func<IList<double>,double>", lambda net: namespace("System.Collections.Generic").List[float](net)),
    ("Func<Vector,Vector>", lambda net: namespace("Numerics.Mathematics.LinearAlgebra").Vector(net)),
])
def test_array_arguments_reach_python_bit_exact_as_owned_arrays(delegate_type: str, wrap: object) -> None:
    # Arrange
    limit = callbacks.BUFFERED_ARRAY_LENGTH
    lengths = (0, 1, 7, limit, limit + 1, 3 * limit, 2)
    # Act
    pairs = _received(delegate_type, wrap, lengths)
    # Assert
    for bits, received in pairs:
        assert isinstance(received, np.ndarray) and received.dtype == np.float64
        assert received.flags.owndata and received.flags.writeable and received.flags.c_contiguous
        np.testing.assert_array_equal(received.view(np.uint64), bits)
    kept = pairs[2][1]
    np.testing.assert_array_equal(kept.view(np.uint64), np.resize(_PATTERNS, 7))


def test_buffered_callbacks_are_reentrant_and_thread_safe() -> None:
    # Arrange
    import threading  # noqa: PLC0415

    outer_seen: list[np.ndarray] = []
    holder: dict[str, object] = {}

    def nested(x: np.ndarray) -> float:
        if x.size == 3:
            outer_seen.append(x)
            holder["delegate"](to_net_doubles([9.0, 8.0]))
            outer_seen.append(x.copy())
        return float(x.sum())

    holder["delegate"] = callbacks.vector_function(nested)
    summed = callbacks.vector_function(lambda x: float(x.sum()))
    errors: list[str] = []

    def worker(seed: int) -> None:
        rng = np.random.default_rng(seed)
        for _ in range(300):
            values = rng.standard_normal(int(rng.integers(1, 40)))
            if summed(to_net_doubles(values)) != float(values.sum()):
                errors.append(f"thread {seed} received another thread's values")
                return

    threads = [threading.Thread(target=worker, args=(seed,)) for seed in range(6)]
    # Act
    nested_value = holder["delegate"](to_net_doubles([1.0, 2.0, 3.0]))
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    # Assert
    assert nested_value == 6.0
    assert outer_seen[0].tolist() == outer_seen[1].tolist() == [1.0, 2.0, 3.0]
    assert errors == []


def test_buffered_callbacks_raise_python_errors_and_count_as_python() -> None:
    # Arrange
    def failing(x: np.ndarray) -> float:
        raise ZeroDivisionError(f"length {x.size}")

    delegate = callbacks.vector_function(failing)
    # Act / Assert
    with pytest.raises(ZeroDivisionError, match="length 2"):
        delegate(to_net_doubles([1.0, 2.0]))
    with pytest.raises(ZeroDivisionError, match=f"length {callbacks.BUFFERED_ARRAY_LENGTH + 1}"):
        delegate(to_net_doubles(np.zeros(callbacks.BUFFERED_ARRAY_LENGTH + 1)))
    assert callbacks.uses_python_callback(delegate)
    assert callbacks.uses_python_callback(callbacks.to_delegate(lambda x: x, "Func<Vector,Vector>"))


def _per_call(delegate: object, argument: object, calls: int = 20000) -> float:
    """Return the best per-call time over several rounds, in seconds."""
    delegate(argument)
    best = math.inf
    for _ in range(5):
        start = time.perf_counter()
        for _ in range(calls):
            delegate(argument)
        best = min(best, (time.perf_counter() - start) / calls)
    return best


def test_array_callback_costs_close_to_a_raw_pythonnet_delegate() -> None:
    # Arrange: measured on a desktop, the overhead over a raw pythonnet delegate is about 0.4 us at
    # 2 to 50 values and 1 us at 1,000; the bound leaves a wide margin for slower machines.
    system = namespace("System")
    raw_type = system.Func[system.Array[system.Double], system.Double]
    constant = lambda x: 0.0  # noqa: E731
    wrapped = callbacks.vector_function(constant)
    raw = raw_type(constant)
    # Act
    overheads = {}
    for length in (2, 5, 50, callbacks.BUFFERED_ARRAY_LENGTH):
        argument = to_net_doubles(np.linspace(0.0, 1.0, length))
        overheads[length] = _per_call(wrapped.Invoke, argument) - _per_call(raw.Invoke, argument)
    # Assert
    assert all(overhead < 2.5e-6 for overhead in overheads.values()), overheads
