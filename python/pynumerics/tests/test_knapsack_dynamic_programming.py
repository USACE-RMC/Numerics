"""Knapsack and discrete dynamic programming against scipy, brute force and textbook values.

The knapsack is checked against ``scipy.optimize.milp`` and the dynamic programs against brute-force
enumeration and NumPy value iteration.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from pynumerics import callbacks, optimization
from pynumerics.optimization import DynamicProgramming

opt = pytest.importorskip("scipy.optimize")
needs_roslyn = pytest.mark.skipif(not callbacks.roslyn_available(), reason="Roslyn assemblies not staged")


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_knapsack_matches_scipy_milp(seed) -> None:
    # Arrange
    rng = np.random.default_rng(seed)
    values, weights = rng.integers(1, 50, size=15).astype(float), rng.integers(1, 30, size=15)
    capacity = int(weights.sum() // 3)
    reference = opt.milp(-values, constraints=opt.LinearConstraint(weights[None, :], 0, capacity),
                         integrality=np.ones(15), bounds=opt.Bounds(0, 1))
    # Act
    total, selected = optimization.knapsack(values, weights, capacity)
    # Assert
    assert total == pytest.approx(-reference.fun, abs=1e-9)
    assert values[selected].sum() == total and weights[selected].sum() <= capacity
    assert selected.dtype == np.int32 and np.all(np.diff(selected) > 0)


def test_multi_constraint_knapsack_matches_scipy_milp() -> None:
    # Arrange
    rng = np.random.default_rng(11)
    values, weights = rng.integers(1, 40, size=10).astype(float), rng.integers(0, 8, size=(10, 2))
    capacities = [15, 12]
    reference = opt.milp(-values, constraints=opt.LinearConstraint(weights.T, 0, capacities), integrality=np.ones(10),
                         bounds=opt.Bounds(0, 1))
    # Act
    total, selected = optimization.knapsack(values, weights, capacities)
    # Assert
    assert total == pytest.approx(-reference.fun, abs=1e-9)
    assert np.all(weights[selected].sum(axis=0) <= capacities)


def test_textbook_knapsack() -> None:
    # Act / Assert: the classic instance (values 60, 100, 120; weights 10, 20, 30; capacity 50).
    total, selected = optimization.knapsack([60.0, 100.0, 120.0], [10, 20, 30], 50)
    assert total == 220.0 and selected.tolist() == [1, 2]
    with pytest.raises(ValueError):
        optimization.knapsack([1.0, 2.0], [1, -1], 3)


@pytest.mark.parametrize(("weights", "capacity"), [
    ([1.5, 2.0], 3), ([1, 2], 3.5), (np.array([[1.0, 0.5], [2.0, 1.0]]), [3, 2]), ([1, 2], [3.0, np.nan])])
def test_knapsack_refuses_fractional_weights_and_capacities(weights, capacity) -> None:
    # Act / Assert: truncating 1.5 to 1 would accept items whose real total weight exceeds the capacity.
    with pytest.raises(ValueError, match="must be integers"):
        optimization.knapsack([1.0, 2.0], weights, capacity)


def test_knapsack_accepts_integral_floats() -> None:
    # Act
    total, selected = optimization.knapsack([60.0, 100.0, 120.0], [10.0, 20.0, 30.0], 50.0)
    # Assert
    assert total == 220.0 and selected.tolist() == [1, 2]


def _inventory() -> tuple[int, int, int, callable, callable]:
    """A deterministic inventory problem: order a units at cost 2a + 1 (if a > 0), hold s at cost 0.5 s, demand 1."""
    stages, states, decisions = 4, 4, 3

    def cost(t: int, s: int, a: int) -> float:
        return (2.0 * a + 1.0 if a > 0 else 0.0) + 0.5 * s

    def step(t: int, s: int, a: int) -> int:
        nxt = s + a - 1
        return nxt if 0 <= nxt < states else -1

    return stages, states, decisions, cost, step


def test_deterministic_dynamic_program_matches_brute_force() -> None:
    # Arrange
    stages, states, decisions, cost, step = _inventory()
    dp = DynamicProgramming(stages, states, decisions, cost, step)
    # Act
    dp.solve()
    path = dp.simulate(1)
    # Assert: enumerate every decision sequence from every state.
    for s0 in range(states):
        best = np.inf
        for plan in itertools.product(range(decisions), repeat=stages):
            s, total = s0, 0.0
            for t, a in enumerate(plan):
                nxt = step(t, s, a)
                if nxt < 0:
                    total = np.inf
                    break
                total, s = total + cost(t, s, a), nxt
            best = min(best, total)
        assert dp.value[0, s0] == pytest.approx(best, abs=1e-12)
    assert dp.is_solved and not dp.is_stochastic and dp.max_outcomes == 1
    assert (dp.stages, dp.states, dp.decisions) == (4, 4, 3)
    assert dp.value.shape == (5, 4) and dp.policy.shape == (4, 4) and dp.policy.dtype == np.int32
    assert path.total_value == dp.value[0, 1] and path.states.size == 5 and path.decisions.size == 4
    np.testing.assert_array_equal(path.decisions, [dp.policy[t, s] for t, s in enumerate(path.states[:-1])])
    assert path.stage_values.sum() == pytest.approx(path.total_value)


def test_maximize_discount_and_terminal_value() -> None:
    # Arrange
    dp = DynamicProgramming(2, 2, 2, lambda t, s, a: float(s + a), lambda t, s, a: a)
    dp.maximize = True
    dp.discount_factor = 0.5
    dp.terminal_value = [0.0, 10.0]
    # Act
    dp.solve()
    # Assert: V1(s) = s + 1 + 0.5 * 10; V0(s) = s + 1 + 0.5 * V1(1).
    np.testing.assert_allclose(dp.value[1], [6.0, 7.0])
    np.testing.assert_allclose(dp.value[0], [1.0 + 3.5, 2.0 + 3.5])
    assert dp.maximize and dp.discount_factor == 0.5 and dp.terminal_value.tolist() == [0.0, 10.0]
    with pytest.raises(ValueError):
        dp.discount_factor = 1.5


def _mdp() -> tuple[np.ndarray, np.ndarray]:
    """A three-state, two-action MDP: rewards r[s, a] and transition probabilities p[s, a, s']."""
    rewards = np.array([[1.0, 0.0], [0.0, 2.0], [-1.0, 0.5]])
    p = np.array([[[0.7, 0.3, 0.0], [0.0, 0.5, 0.5]],
                  [[0.4, 0.6, 0.0], [0.0, 0.2, 0.8]],
                  [[0.5, 0.0, 0.5], [0.1, 0.1, 0.8]]])
    return rewards, p


def _transition(p: np.ndarray):
    def write(t: int, s: int, a: int, next_states, probabilities) -> int:
        count = 0
        for nxt in range(p.shape[2]):
            if p[s, a, nxt] > 0:
                next_states[count] = nxt
                probabilities[count] = float(p[s, a, nxt])
                count += 1
        return count

    return write


def test_stochastic_dynamic_program_matches_numpy_value_iteration() -> None:
    # Arrange
    rewards, p = _mdp()
    gamma = 0.9
    dp = DynamicProgramming(6, 3, 2, lambda t, s, a: float(rewards[s, a]), _transition(p), 3)
    dp.maximize = True
    dp.discount_factor = gamma
    value = np.zeros(3)
    finite = [value]
    for _ in range(6):
        value = np.max(rewards + gamma * p @ value, axis=1)
        finite.append(value)
    # Act
    dp.solve()
    converged = dp.solve_steady_state(1e-12, 100000)
    stationary = np.zeros(3)
    for _ in range(2000):
        stationary = np.max(rewards + gamma * p @ stationary, axis=1)
    # Assert
    assert dp.is_stochastic and dp.max_outcomes == 3
    np.testing.assert_allclose(dp.value[::-1], np.array(finite), atol=1e-12)
    assert converged and dp.steady_state_converged and dp.steady_state_iterations > 10
    np.testing.assert_allclose(dp.steady_state_value, stationary, atol=1e-10)
    np.testing.assert_array_equal(dp.steady_state_policy, np.argmax(rewards + gamma * p @ stationary, axis=1))
    sampled = dp.simulate(0, seed=5)
    replayed = dp.simulate(0, outcomes=[0] * 6)
    assert sampled.states.size == 7 and replayed.states[1] == 0
    np.testing.assert_array_equal(dp.simulate(0, seed=5).states, sampled.states)
    with pytest.raises(ValueError, match="not both"):
        dp.simulate(0, outcomes=[0] * 6, seed=1)


@needs_roslyn
def test_dynamic_program_parallel_policy_and_compiled_callbacks() -> None:
    # Arrange
    stages, states, decisions, cost, step = _inventory()
    python = DynamicProgramming(stages, states, decisions, cost, step)
    python.use_parallel = True
    compiled = DynamicProgramming(stages, states, decisions, "(t, s, a) => (a > 0 ? 2.0 * a + 1.0 : 0.0) + 0.5 * s",
                                  "(t, s, a) => s + a - 1 >= 0 && s + a - 1 < 4 ? s + a - 1 : -1")
    compiled.use_parallel = True
    # Act
    python.solve()
    compiled.solve()
    # Assert
    assert not python.use_parallel and compiled.use_parallel
    np.testing.assert_array_equal(compiled.value, python.value)
    np.testing.assert_array_equal(compiled.policy, python.policy)
