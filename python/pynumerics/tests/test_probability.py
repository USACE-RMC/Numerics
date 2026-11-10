"""Probability-combination wrappers against scipy (bivariate and multivariate normal, quadrature).

Exact rules (independence, the Frechet-Hoeffding bounds) are checked to rounding; the PCM method is
an approximation of the multivariate normal, checked against scipy at its accuracy.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from pynumerics import multivariate
from pynumerics import probability as pr

stats = pytest.importorskip("scipy.stats")
integrate = pytest.importorskip("scipy.integrate")

P = np.array([0.1, 0.25, 0.4])
R = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])


def _mvn_joint(p: np.ndarray, corr: np.ndarray, subset: tuple[int, ...]) -> float:
    """Joint probability that every event in ``subset`` occurs, by scipy's multivariate normal."""
    idx = list(subset)
    if len(idx) == 1:
        return float(p[idx[0]])
    z = stats.norm.ppf(p[idx])
    return float(stats.multivariate_normal(np.zeros(len(idx)), corr[np.ix_(idx, idx)], seed=0).cdf(z))


def _mvn_union(p: np.ndarray, corr: np.ndarray) -> float:
    total = 0.0
    for k in range(1, p.size + 1):
        for subset in itertools.combinations(range(p.size), k):
            total += (-1) ** (k + 1) * _mvn_joint(p, corr, subset)
    return total


def _combinations(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Subset counts and indicator rows in subset-size order (lexicographic within a size)."""
    rows, counts = [], []
    for k in range(1, n + 1):
        subsets = list(itertools.combinations(range(n), k))
        counts.append(len(subsets))
        for subset in subsets:
            row = np.zeros(n, dtype=np.int32)
            row[list(subset)] = 1
            rows.append(row)
    return np.array(counts, dtype=np.int32), np.array(rows)


def test_two_event_rules_equal_the_bivariate_normal() -> None:
    # Arrange
    a, b, rho = 0.15, 0.62, 0.45
    za, zb = stats.norm.ppf([a, b])
    joint = stats.multivariate_normal([0.0, 0.0], [[1.0, rho], [rho, 1.0]]).cdf([za, zb])
    # Act / Assert
    assert pr.a_and_b(a, b, rho) == pytest.approx(joint, abs=1e-7)
    assert pr.a_or_b(a, b, rho) == pytest.approx(a + b - joint, abs=1e-7)
    assert pr.a_not_b(a, b, rho) == pytest.approx(a - joint, abs=1e-7)
    assert pr.b_not_a(a, b, rho) == pytest.approx(b - joint, abs=1e-7)
    assert pr.a_given_b(a, b, rho) == pytest.approx(joint / b, abs=1e-6)
    assert pr.b_given_a(a, b, rho) == pytest.approx(joint / a, abs=1e-6)
    assert pr.a_and_b(a, b) == pytest.approx(a * b, rel=1e-12)


def test_joint_probabilities_under_the_closed_form_rules() -> None:
    # Act / Assert
    assert pr.independent_joint_probability(P) == pytest.approx(np.prod(P), rel=1e-15)
    assert pr.independent_joint_probability(P, [1, 0, 1]) == pytest.approx(0.1 * 0.4, rel=1e-15)
    assert pr.positive_joint_probability(P) == 0.1
    assert pr.positive_joint_probability(P, [0, 1, 1]) == 0.25
    assert pr.negative_joint_probability([0.8, 0.9]) == pytest.approx(0.7, abs=1e-15)
    assert pr.negative_joint_probability([0.8, 0.9, 0.5], [1, 1, 0]) == pytest.approx(0.7, abs=1e-15)
    assert pr.joint_probability(P) == pr.independent_joint_probability(P)
    assert pr.joint_probability(P, dependency="perfectly_positive") == 0.1
    assert pr.joint_probability(P, [1, 1, 1], dependency=pr.DependencyType.INDEPENDENT) == pytest.approx(np.prod(P))
    assert math.isnan(pr.joint_probability(P, [1, 1, 1]))  # correlation-matrix default needs a matrix
    with pytest.raises(ValueError):
        pr.joint_probability(P, correlation_matrix=R)


def test_pcm_joint_probabilities_approximate_the_multivariate_normal() -> None:
    # Arrange
    exact = _mvn_joint(P, R, (0, 1, 2))
    # Act
    hpcm, conditional = pr.joint_probability_hpcm(P, [1, 1, 1], R, return_conditional=True)
    pcm = pr.joint_probability_pcm(P, [1, 1, 1], R)
    # Assert
    assert hpcm == pytest.approx(exact, rel=0.05)
    assert pcm == pytest.approx(exact, rel=0.1)
    assert pr.joint_probability(P, [1, 1, 1], R) == hpcm
    assert conditional[0] == pytest.approx(P[0], rel=1e-12)
    assert pr.joint_probability_hpcm(P, [1, 1, 1], np.eye(3)) == pytest.approx(np.prod(P), rel=1e-10)
    counts, rows = _combinations(3)
    batch = pr.joint_probabilities_pcm(P, rows, R)
    np.testing.assert_allclose(batch[:3], P, rtol=1e-15)  # single-event rows are the marginals
    assert batch[-1] == pytest.approx(pr.joint_probability_pcm(P, [1, 1, 1], R), rel=1e-12)
    assert counts.tolist() == [3, 3, 1]


def test_mvn_joint_and_union_equal_scipy() -> None:
    # Arrange
    mvn = multivariate.MultivariateNormal(np.zeros(2), R[:2, :2])
    p2 = P[:2]
    counts, rows = _combinations(2)
    # Act / Assert
    assert pr.joint_probability_mvn(p2, [1, 1], mvn) == pytest.approx(_mvn_joint(p2, R, (0, 1)), abs=1e-7)
    batch = pr.joint_probabilities_mvn(p2, rows, mvn)
    assert batch[2] == pytest.approx(_mvn_joint(p2, R, (0, 1)), abs=1e-7)
    assert pr.union_mvn(p2, mvn) == pytest.approx(_mvn_union(p2, R[:2, :2]), abs=1e-7)
    assert pr.union_mvn(p2, mvn.net, counts, rows) == pytest.approx(pr.union_mvn(p2, mvn), abs=1e-12)
    exclusive = pr.exclusive_mvn(p2, mvn)
    assert exclusive.sum() == pytest.approx(pr.union_mvn(p2, mvn), abs=1e-12)
    events, indicators = pr.exclusive_mvn(p2, mvn, counts, rows, return_events=True)
    assert events.shape[0] == indicators.shape[0] and indicators.dtype == np.int32


def test_unions_equal_the_closed_forms_scipy_and_quadrature() -> None:
    # Arrange
    rho = 0.35
    b = stats.norm.ppf(P)

    def integrand(z: float) -> float:
        conditional = stats.norm.cdf((b - math.sqrt(rho) * z) / math.sqrt(1 - rho))
        return stats.norm.pdf(z) * (1 - np.prod(1 - conditional))

    single_factor = integrate.quad(integrand, -10, 10, epsabs=1e-13, epsrel=1e-12)[0]
    equi = np.full((3, 3), rho) + (1 - rho) * np.eye(3)
    # Act / Assert
    assert pr.independent_union(P) == pytest.approx(1 - np.prod(1 - P), rel=1e-15)
    assert pr.union(P) == pr.independent_union(P)
    assert pr.positively_dependent_union(P) == 0.4
    assert pr.union(P, "perfectly_positive") == 0.4
    assert pr.negatively_dependent_union(P) == pytest.approx(min(1.0, P.sum()), rel=1e-15)
    assert pr.union_single_factor(P, rho) == pytest.approx(single_factor, rel=1e-7)
    assert pr.union_single_factor(P, rho) == pytest.approx(_mvn_union(P, equi), abs=2e-5)  # scipy QMC cdf
    conditional = pr.single_factor_conditional_probabilities(b, rho, 0.7)
    np.testing.assert_allclose(conditional, stats.norm.cdf((b - math.sqrt(rho) * 0.7) / math.sqrt(1 - rho)),
                               rtol=1e-12)


def test_pcm_unions_and_exclusive_probabilities_are_consistent() -> None:
    # Arrange
    counts, rows = _combinations(3)
    exact = _mvn_union(P, R)
    # Act
    union = pr.union_pcm(P, R)
    lazy, status = pr.union_pcm_lazy(P, R)
    explicit = pr.union_pcm(P, R, counts, rows)
    value, events, indicators = pr.union_pcm(P, R, counts, rows, return_events=True)
    exclusive = pr.exclusive_pcm(P, R)
    lazy_events, lazy_rows, lazy_status = pr.exclusive_pcm_lazy(P, R)
    # Assert
    assert union == pytest.approx(exact, abs=5e-3)
    assert lazy == pytest.approx(union, abs=1e-12) and status is pr.ExclusiveEnumerationStatus.COMPLETE
    assert explicit == pytest.approx(union, abs=1e-12) and value == pytest.approx(union, abs=1e-12)
    assert events.size == indicators.shape[0]
    assert exclusive.sum() == pytest.approx(exact, abs=5e-3)
    np.testing.assert_allclose(pr.exclusive_pcm(P, R, counts, rows), exclusive, rtol=1e-12)
    ev, ind = pr.exclusive_pcm(P, R, counts, rows, return_events=True)
    assert ind.shape == (ev.size, 3)
    assert lazy_rows.shape == (lazy_events.size, 3) and lazy_status is pr.ExclusiveEnumerationStatus.COMPLETE
    with pytest.raises(ValueError):
        pr.union_pcm(P, R, return_events=True)


def test_independent_and_positive_exclusive_probabilities_equal_the_closed_forms() -> None:
    # Arrange
    counts, rows = _combinations(3)
    expected = np.array([np.prod(np.where(r == 1, P, 1 - P)) for r in rows])
    # Act
    full = pr.independent_exclusive(P)
    one = pr.independent_exclusive(P, rows[4])
    batch = pr.independent_exclusive(P, rows)
    events, indicators = pr.independent_exclusive(P, rows, counts)
    lazy_events, lazy_rows, status = pr.independent_exclusive_lazy(P, include_no_event_row=True)
    _, _, capped = pr.independent_exclusive_lazy(P, max_emitted_combinations=2)
    positive = pr.positively_dependent_exclusive(P)
    pos_events, pos_rows = pr.positively_dependent_exclusive(P, rows, counts)
    pos_lazy, _, pos_status = pr.positively_dependent_exclusive_lazy(P)
    # Assert
    np.testing.assert_allclose(full, expected, rtol=1e-14)
    assert one == pytest.approx(expected[4], rel=1e-14)
    np.testing.assert_allclose(batch, expected, rtol=1e-14)
    assert events.sum() == pytest.approx(1 - np.prod(1 - P), rel=1e-12)
    assert indicators.shape == (events.size, 3)
    assert lazy_events[0] == pytest.approx(np.prod(1 - P), rel=1e-14) and lazy_rows[0].tolist() == [0, 0, 0]
    assert status is pr.ExclusiveEnumerationStatus.COMPLETE and capped is pr.ExclusiveEnumerationStatus.CAPPED
    assert positive.sum() == pytest.approx(0.4, rel=1e-14)
    assert pos_events.sum() == pytest.approx(0.4, rel=1e-12) and pos_rows.shape[1] == 3
    assert pos_lazy.sum() == pytest.approx(0.4, rel=1e-12) and pos_status is pr.ExclusiveEnumerationStatus.COMPLETE


def test_adjustment_factors() -> None:
    # Act / Assert
    assert pr.common_cause_adjustment(P) == pytest.approx((1 - np.prod(1 - P)) / P.sum(), rel=1e-12)
    assert pr.common_cause_adjustment(P, np.eye(3)) == pytest.approx(pr.common_cause_adjustment(P), rel=1e-10)
    assert pr.common_cause_adjustment([0.3]) == 1.0
    assert 0.0 < pr.mutually_exclusive_adjustment(P) <= 1.0
    assert pr.DependencyType("CorrelationMatrix") is pr.DependencyType.CORRELATION_MATRIX


def test_indicator_matrices_cross_in_one_block_copy() -> None:
    # Arrange
    rng = np.random.default_rng(0)
    rows = (rng.random((50_000, 6)) < 0.5).astype(np.int64)
    p = np.full(6, 0.2)
    # Act
    result = pr.independent_exclusive(p, rows)
    # Assert
    expected = np.prod(np.where(rows == 1, p, 1 - p), axis=1)
    np.testing.assert_allclose(result, expected, rtol=1e-13)
    with pytest.raises(ValueError):
        pr._int_matrix(np.zeros(3))
