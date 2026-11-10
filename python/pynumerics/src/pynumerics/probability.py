"""Probability combination for risk analysis from ``Numerics.Data.Statistics.Probability``.

Joint, union and mutually exclusive probabilities of events under independence, perfect positive
or negative dependence, a correlation matrix through the product of conditional marginals (PCM),
a multivariate normal, or an equicorrelated single-factor Gaussian structure. Each function calls
one Numerics member.

Indicator rows mark which events occur (1) or do not (0). ``binomial_combinations`` and
``indicators`` are the subset counts and indicator rows of the inclusion-exclusion expansion in
subset-size order; the forms without them enumerate every combination. Event lists come back as
a probability vector and an int32 indicator matrix with one row per exclusive event.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics._util import namespace, net_enum
from pynumerics.enums import NetEnum
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _int32_array,
    _net_rectangular,
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_ints,
)

_TOLERANCE = 1e-4
"""Numerics' default absolute and relative inclusion-exclusion tolerances."""


class DependencyType(NetEnum):
    """Dependence between events (``Probability.DependencyType``)."""

    INDEPENDENT = "Independent"
    PERFECTLY_POSITIVE = "PerfectlyPositive"
    PERFECTLY_NEGATIVE = "PerfectlyNegative"
    CORRELATION_MATRIX = "CorrelationMatrix"


class ExclusiveEnumerationStatus(NetEnum):
    """Outcome of a lazily enumerated exclusive-probability expansion (``Probability.ExclusiveEnumerationStatus``)."""

    COMPLETE = "Complete"
    CONVERGED = "Converged"
    CAPPED = "Capped"


def _probability() -> Any:
    """Return ``Numerics.Data.Statistics.Probability`` after starting the runtime."""
    return namespace("Numerics.Data.Statistics").Probability


def _dependency(value: DependencyType | str) -> Any:
    """Map a dependency to the .NET enum member."""
    return net_enum(_probability().DependencyType, value)


def _status(value: Any) -> ExclusiveEnumerationStatus:
    """Map a .NET enumeration status to the Python enum."""
    return ExclusiveEnumerationStatus(str(value))


def _int_matrix(values: Any) -> Any:
    """Copy a 2-D integer array into a new ``Int32[,]`` in one copy, refusing values outside the int32 range."""
    array = np.asarray(values)
    if array.ndim != 2:
        raise ValueError(f"indicators must be two-dimensional; got shape {array.shape}.")
    return _net_rectangular(_int32_array(array))


def _rows(indicators: Any) -> tuple[Any, bool]:
    """Return the indicators as ``Int32[]`` or ``Int32[,]`` and whether they are one row."""
    array = np.asarray(indicators)
    if array.ndim == 1:
        return to_net_ints(array), True
    return _int_matrix(array), False


def _mvn(multivariate_normal: Any) -> Any:
    """Return the .NET ``MultivariateNormal`` behind a wrapper (or the object itself)."""
    return getattr(multivariate_normal, "net", multivariate_normal)


def _event_lists() -> tuple[Any, Any]:
    """Return new empty .NET ``List<double>`` and ``List<int[]>`` output lists."""
    import System  # noqa: PLC0415 - the runtime is running once _probability returns
    from System.Collections.Generic import List  # noqa: PLC0415

    return List[System.Double](), List[System.Array[System.Int32]]()


def _events(probabilities: Any, indicators: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return .NET event lists as a probability vector and an int32 indicator matrix."""
    values = to_numpy_doubles(probabilities)
    rows = [to_numpy_ints(row) for row in indicators]
    if not rows:
        return values, np.empty((0, 0), dtype=np.int32)
    return values, np.vstack(rows)


def _pair_member(name: str, a: float, b: float, rho: float) -> float:
    return float(getattr(_probability(), name)(float(a), float(b), float(rho)))


def a_and_b(a: float, b: float, rho: float = 0.0) -> float:
    """Return P(A and B) of two events with correlated Gaussian variates (``Probability.AAndB``).

    Parameters
    ----------
    a, b : float
        The marginal probabilities.
    rho : float, optional
        The Pearson correlation of the underlying standard normal variates; default 0.

    Returns
    -------
    float
        The joint probability, ``a b`` at zero correlation.
    """
    return _pair_member("AAndB", a, b, rho)


def a_or_b(a: float, b: float, rho: float = 0.0) -> float:
    """Return P(A or B) of two events (``Probability.AOrB``).

    Parameters
    ----------
    a, b : float
        The marginal probabilities.
    rho : float, optional
        The correlation of the underlying standard normal variates; default 0.

    Returns
    -------
    float
        ``a + b - P(A and B)``.
    """
    return _pair_member("AOrB", a, b, rho)


def a_not_b(a: float, b: float, rho: float = 0.0) -> float:
    """Return P(A and not B) of two events (``Probability.ANotB``).

    Parameters
    ----------
    a, b : float
        The marginal probabilities.
    rho : float, optional
        The correlation of the underlying standard normal variates; default 0.

    Returns
    -------
    float
        ``a - P(A and B)``.
    """
    return _pair_member("ANotB", a, b, rho)


def b_not_a(a: float, b: float, rho: float = 0.0) -> float:
    """Return P(B and not A) of two events (``Probability.BNotA``).

    Parameters
    ----------
    a, b : float
        The marginal probabilities.
    rho : float, optional
        The correlation of the underlying standard normal variates; default 0.

    Returns
    -------
    float
        ``b - P(A and B)``.
    """
    return _pair_member("BNotA", a, b, rho)


def a_given_b(a: float, b: float, rho: float = 0.0) -> float:
    """Return P(A | B) of two events (``Probability.AGivenB``).

    Parameters
    ----------
    a, b : float
        The marginal probabilities.
    rho : float, optional
        The correlation of the underlying standard normal variates; default 0.

    Returns
    -------
    float
        ``P(A and B) / b``.
    """
    return _pair_member("AGivenB", a, b, rho)


def b_given_a(a: float, b: float, rho: float = 0.0) -> float:
    """Return P(B | A) of two events (``Probability.BGivenA``).

    Parameters
    ----------
    a, b : float
        The marginal probabilities.
    rho : float, optional
        The correlation of the underlying standard normal variates; default 0.

    Returns
    -------
    float
        ``P(A and B) / a``.
    """
    return _pair_member("BGivenA", a, b, rho)


def joint_probability(probabilities: Any, indicators: Any = None, correlation_matrix: Any = None,
                      dependency: DependencyType | str | None = None) -> float:
    """Return the joint probability of events under a dependence model (``Probability.JointProbability``).

    Without ``indicators`` every event occurs and ``dependency`` defaults to independence. With
    ``indicators`` the indicated events occur and the others do not; ``dependency`` then defaults to
    the correlation matrix (the PCM method of :func:`joint_probability_hpcm`), and NaN results when
    no correlation matrix is given.

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like, optional
        One 0/1 indicator per event.
    correlation_matrix : array_like, optional
        The correlation matrix, used with indicators and the correlation-matrix dependency.
    dependency : DependencyType or str, optional
        The dependence model.

    Returns
    -------
    float
        The joint probability.
    """
    p = _probability()
    values = to_net_doubles(probabilities)
    if indicators is None:
        if correlation_matrix is not None:
            raise ValueError("correlation_matrix requires indicators.")
        return float(p.JointProbability(values, _dependency(dependency or DependencyType.INDEPENDENT)))
    matrix = None if correlation_matrix is None else to_net_matrix(correlation_matrix)
    return float(p.JointProbability(values, to_net_ints(indicators), matrix,
                                    _dependency(dependency or DependencyType.CORRELATION_MATRIX)))


def _indicated(name: str, probabilities: Any, indicators: Any) -> float:
    member = getattr(_probability(), name)
    values = to_net_doubles(probabilities)
    return float(member(values) if indicators is None else member(values, to_net_ints(indicators)))


def independent_joint_probability(probabilities: Any, indicators: Any = None) -> float:
    """Return the joint probability under independence (``Probability.IndependentJointProbability``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like, optional
        One 0/1 indicator per event; only indicated events take part. Default all.

    Returns
    -------
    float
        The product of the participating probabilities.
    """
    return _indicated("IndependentJointProbability", probabilities, indicators)


def positive_joint_probability(probabilities: Any, indicators: Any = None) -> float:
    """Return the joint probability under perfect positive dependence (``Probability.PositiveJointProbability``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like, optional
        One 0/1 indicator per event; only indicated events take part. Default all.

    Returns
    -------
    float
        The smallest participating probability (the Frechet-Hoeffding upper bound).
    """
    return _indicated("PositiveJointProbability", probabilities, indicators)


def negative_joint_probability(probabilities: Any, indicators: Any = None) -> float:
    """Return the Frechet-Hoeffding lower bound of the joint probability (``Probability.NegativeJointProbability``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like, optional
        One 0/1 indicator per event; only indicated events take part. Default all.

    Returns
    -------
    float
        ``max(0, sum p - (k - 1))`` over the k participating events; 1 when none participate.
    """
    return _indicated("NegativeJointProbability", probabilities, indicators)


def _pcm(name: str, probabilities: Any, indicators: Any, correlation_matrix: Any,
         return_conditional: bool) -> float | tuple[float, np.ndarray]:
    import System  # noqa: PLC0415 - the runtime is running once _probability returns

    member = getattr(_probability(), name)
    values = np.asarray(probabilities, dtype=np.float64)
    conditional = System.Array[System.Double](values.size) if return_conditional else None
    result = float(member(to_net_doubles(values), to_net_ints(indicators), to_net_matrix(correlation_matrix),
                          conditional))
    return (result, to_numpy_doubles(conditional)) if return_conditional else result


def joint_probability_hpcm(probabilities: Any, indicators: Any, correlation_matrix: Any,
                           return_conditional: bool = False) -> float | tuple[float, np.ndarray]:
    """Return the joint probability by the modified PCM method (``Probability.JointProbabilityHPCM``).

    Haden Smith's modification of Pandey's product of conditional marginals.

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like
        One 0/1 indicator per event.
    correlation_matrix : array_like
        The symmetric correlation matrix of the underlying normal variates.
    return_conditional : bool, optional
        Also return the conditional probability of each event; default False.

    Returns
    -------
    probability : float
        The joint probability.
    conditional : numpy.ndarray
        Only when ``return_conditional`` is True: the conditional probabilities.
    """
    return _pcm("JointProbabilityHPCM", probabilities, indicators, correlation_matrix, return_conditional)


def joint_probability_pcm(probabilities: Any, indicators: Any, correlation_matrix: Any,
                          return_conditional: bool = False) -> float | tuple[float, np.ndarray]:
    """Return the joint probability by Pandey's PCM method (``Probability.JointProbabilityPCM``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like
        One 0/1 indicator per event.
    correlation_matrix : array_like
        The symmetric correlation matrix of the underlying normal variates.
    return_conditional : bool, optional
        Also return the conditional probability of each event; default False.

    Returns
    -------
    probability : float
        The joint probability.
    conditional : numpy.ndarray
        Only when ``return_conditional`` is True: the conditional probabilities.
    """
    return _pcm("JointProbabilityPCM", probabilities, indicators, correlation_matrix, return_conditional)


def joint_probabilities_pcm(probabilities: Any, indicators: Any, correlation_matrix: Any) -> np.ndarray:
    """Return the PCM joint probability of each indicator row (``Probability.JointProbabilitiesPCM``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like
        A 2-D 0/1 array, one row per combination, in subset-size order. Precondition: the first
        ``len(probabilities)`` rows must be the singleton rows, row ``i`` flagging event ``i`` alone.
        Numerics returns ``probabilities[i]`` for those rows without reading them and does not check
        the layout.
    correlation_matrix : array_like
        The correlation matrix.

    Returns
    -------
    numpy.ndarray
        One joint probability per row.
    """
    p = _probability()
    return to_numpy_doubles(p.JointProbabilitiesPCM(to_net_doubles(probabilities), _int_matrix(indicators),
                                                    to_net_matrix(correlation_matrix)))


def joint_probability_mvn(probabilities: Any, indicators: Any, multivariate_normal: Any) -> float:
    """Return the joint probability from a multivariate normal (``Probability.JointProbabilityMVN``).

    Dimensions above two use the distribution's randomized lattice, which advances with each call.

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like
        One 0/1 indicator per event.
    multivariate_normal : MultivariateNormal
        The standard multivariate normal of the underlying variates (a wrapper or the .NET object).

    Returns
    -------
    float
        The joint probability.
    """
    p = _probability()
    return float(p.JointProbabilityMVN(to_net_doubles(probabilities), to_net_ints(indicators),
                                       _mvn(multivariate_normal)))


def joint_probabilities_mvn(probabilities: Any, indicators: Any, multivariate_normal: Any) -> np.ndarray:
    """Return the multivariate normal joint probability of each indicator row (``Probability.JointProbabilitiesMVN``).

    Rows are evaluated serially in order, so a freshly seeded distribution gives a reproducible batch.

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like
        A 2-D 0/1 array, one row per combination, in subset-size order. Precondition: the first
        ``len(probabilities)`` rows must be the singleton rows, row ``i`` flagging event ``i`` alone.
        Numerics returns ``probabilities[i]`` for those rows without reading them and does not check
        the layout.
    multivariate_normal : MultivariateNormal
        The standard multivariate normal of the underlying variates.

    Returns
    -------
    numpy.ndarray
        One joint probability per row.
    """
    p = _probability()
    return to_numpy_doubles(p.JointProbabilitiesMVN(to_net_doubles(probabilities), _int_matrix(indicators),
                                                    _mvn(multivariate_normal)))


def union(probabilities: Any, dependency: DependencyType | str = DependencyType.INDEPENDENT) -> float:
    """Return the probability of the union of events (``Probability.Union``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    dependency : DependencyType or str, optional
        Independence (default), perfect positive or perfect negative dependence.

    Returns
    -------
    float
        The union probability.
    """
    return float(_probability().Union(to_net_doubles(probabilities), _dependency(dependency)))


def independent_union(probabilities: Any) -> float:
    """Return the union probability under independence, ``1 - prod(1 - p)`` (``Probability.IndependentUnion``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.

    Returns
    -------
    float
        The union probability.
    """
    return float(_probability().IndependentUnion(to_net_doubles(probabilities)))


def positively_dependent_union(probabilities: Any) -> float:
    """Return the union probability under perfect positive dependence, ``max p`` (``PositivelyDependentUnion``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.

    Returns
    -------
    float
        The unimodal lower bound of the union.
    """
    return float(_probability().PositivelyDependentUnion(to_net_doubles(probabilities)))


def negatively_dependent_union(probabilities: Any) -> float:
    """Return the union probability under perfect negative dependence, ``min(1, sum p)`` (``NegativelyDependentUnion``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.

    Returns
    -------
    float
        The unimodal upper bound of the union.
    """
    return float(_probability().NegativelyDependentUnion(to_net_doubles(probabilities)))


def union_single_factor(probabilities: Any, rho: float, relative_tolerance: float = 1e-8) -> float:
    """Return the union probability under an equicorrelated single-factor Gaussian (``UnionSingleFactor``).

    ``P = int phi(z) [1 - prod(1 - Phi((b_i - sqrt(rho) z) / sqrt(1 - rho)))] dz`` with
    ``b_i = Phi^-1(p_i)``: the events are independent given one shared standard normal factor.

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    rho : float
        The common correlation, in [0, 1].
    relative_tolerance : float, optional
        The quadrature relative tolerance; default 1e-8.

    Returns
    -------
    float
        The union probability.
    """
    return float(_probability().UnionSingleFactor(to_net_doubles(probabilities), float(rho),
                                                  float(relative_tolerance)))


def single_factor_conditional_probabilities(normal_thresholds: Any, rho: float, z: float) -> np.ndarray:
    """Return ``Phi((b_i - sqrt(rho) z) / sqrt(1 - rho))`` per threshold (``SingleFactorConditionalProbabilities``).

    Parameters
    ----------
    normal_thresholds : array_like
        The standard normal thresholds ``b_i = Phi^-1(p_i)``.
    rho : float
        The common correlation, in [0, 1).
    z : float
        The shared factor value.

    Returns
    -------
    numpy.ndarray
        The conditional event probabilities.
    """
    import System  # noqa: PLC0415 - the runtime is running once _probability returns

    p = _probability()
    thresholds = to_net_doubles(normal_thresholds)
    out = System.Array[System.Double](int(thresholds.Length))
    p.SingleFactorConditionalProbabilities(thresholds, float(rho), float(z), out)
    return to_numpy_doubles(out)


def union_pcm(probabilities: Any, correlation_matrix: Any, binomial_combinations: Any = None, indicators: Any = None,
              absolute_tolerance: float = _TOLERANCE, relative_tolerance: float = _TOLERANCE,
              return_events: bool = False) -> float | tuple[float, np.ndarray, np.ndarray]:
    """Return the union probability by inclusion-exclusion with PCM joint probabilities (``Probability.UnionPCM``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    correlation_matrix : array_like
        The correlation matrix.
    binomial_combinations : array_like, optional
        The number of combinations of each subset size; given together with ``indicators``.
    indicators : array_like, optional
        The 2-D indicator rows in subset-size order.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances of the expansion; default 1e-4.
    return_events : bool, optional
        Also return the exclusive event probabilities and indicators; needs
        ``binomial_combinations`` and ``indicators``. Default False.

    Returns
    -------
    probability : float
        The union probability.
    event_probabilities : numpy.ndarray
        Only when ``return_events`` is True.
    event_indicators : numpy.ndarray
        Only when ``return_events`` is True: one int32 row per event.
    """
    p = _probability()
    values, matrix = to_net_doubles(probabilities), to_net_matrix(correlation_matrix)
    atol, rtol = float(absolute_tolerance), float(relative_tolerance)
    if binomial_combinations is None and indicators is None:
        if return_events:
            raise ValueError("return_events needs binomial_combinations and indicators.")
        return float(p.UnionPCM(values, matrix, atol, rtol))
    combinations, rows = to_net_ints(binomial_combinations), _int_matrix(indicators)
    if not return_events:
        return float(p.UnionPCM(values, combinations, rows, matrix, atol, rtol))
    result, events, event_rows = p.UnionPCM(values, combinations, rows, matrix, None, None, atol, rtol)
    return (float(result), *_events(events, event_rows))


def union_pcm_lazy(probabilities: Any, correlation_matrix: Any, absolute_tolerance: float = _TOLERANCE,
                   relative_tolerance: float = _TOLERANCE) -> tuple[float, ExclusiveEnumerationStatus]:
    """Return the PCM union probability with combinations generated lazily (``Probability.UnionPCMLazy``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    correlation_matrix : array_like
        The correlation matrix.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances; default 1e-4.

    Returns
    -------
    probability : float
        The union probability.
    status : ExclusiveEnumerationStatus
        Whether the expansion completed or converged early.
    """
    p = _probability()
    result, status = p.UnionPCMLazy(to_net_doubles(probabilities), to_net_matrix(correlation_matrix),
                                    p.ExclusiveEnumerationStatus.Complete, float(absolute_tolerance),
                                    float(relative_tolerance))
    return float(result), _status(status)


def union_mvn(probabilities: Any, multivariate_normal: Any, binomial_combinations: Any = None,
              indicators: Any = None) -> float:
    """Return the union probability by inclusion-exclusion with multivariate normal joints (``Probability.UnionMVN``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    multivariate_normal : MultivariateNormal
        The standard multivariate normal of the underlying variates.
    binomial_combinations : array_like, optional
        The number of combinations of each subset size; given together with ``indicators``.
    indicators : array_like, optional
        The 2-D indicator rows in subset-size order.

    Returns
    -------
    float
        The union probability.
    """
    p = _probability()
    values, mvn = to_net_doubles(probabilities), _mvn(multivariate_normal)
    if binomial_combinations is None and indicators is None:
        return float(p.UnionMVN(values, mvn))
    return float(p.UnionMVN(values, to_net_ints(binomial_combinations), _int_matrix(indicators), mvn))


def _exclusive(name: str, probabilities: Any, indicators: Any, binomial_combinations: Any,
               absolute_tolerance: float, relative_tolerance: float) -> Any:
    """Dispatch the shared overload set of ``IndependentExclusive`` and ``PositivelyDependentExclusive``."""
    member = getattr(_probability(), name)
    values = to_net_doubles(probabilities)
    if binomial_combinations is not None:
        if indicators is None:
            raise ValueError("binomial_combinations needs indicators.")
        events, rows = member(values, to_net_ints(binomial_combinations), _int_matrix(indicators), None, None,
                              float(absolute_tolerance), float(relative_tolerance))
        return _events(events, rows)
    if indicators is None:
        return to_numpy_doubles(member(values))
    net_rows, single = _rows(indicators)
    result = member(values, net_rows)
    return float(result) if single else to_numpy_doubles(result)


def independent_exclusive(probabilities: Any, indicators: Any = None, binomial_combinations: Any = None,
                          absolute_tolerance: float = _TOLERANCE, relative_tolerance: float = _TOLERANCE) -> Any:
    """Return mutually exclusive event probabilities under independence (``Probability.IndependentExclusive``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like, optional
        One 0/1 row (one exclusive probability) or a 2-D array of rows; default every combination.
    binomial_combinations : array_like, optional
        The number of combinations of each subset size; with 2-D ``indicators`` this runs the
        convergent inclusion-exclusion expansion and returns the event lists.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances of the expansion; default 1e-4.

    Returns
    -------
    float, numpy.ndarray or tuple of numpy.ndarray
        The exclusive probability of one row, one per row (or per combination), or, with
        ``binomial_combinations``, the event probabilities and their int32 indicator rows.
    """
    return _exclusive("IndependentExclusive", probabilities, indicators, binomial_combinations, absolute_tolerance,
                      relative_tolerance)


def positively_dependent_exclusive(probabilities: Any, indicators: Any = None, binomial_combinations: Any = None,
                                   absolute_tolerance: float = _TOLERANCE,
                                   relative_tolerance: float = _TOLERANCE) -> Any:
    """Return exclusive probabilities under perfect positive dependence (``PositivelyDependentExclusive``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    indicators : array_like, optional
        One 0/1 row or a 2-D array of rows; default every combination.
    binomial_combinations : array_like, optional
        The number of combinations of each subset size; with 2-D ``indicators`` this runs the
        convergent expansion and returns the event lists.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances of the expansion; default 1e-4.

    Returns
    -------
    float, numpy.ndarray or tuple of numpy.ndarray
        As for :func:`independent_exclusive`.
    """
    return _exclusive("PositivelyDependentExclusive", probabilities, indicators, binomial_combinations,
                      absolute_tolerance, relative_tolerance)


def independent_exclusive_lazy(probabilities: Any, include_no_event_row: bool = False,
                               max_emitted_combinations: int = 0, absolute_tolerance: float = _TOLERANCE,
                               relative_tolerance: float = _TOLERANCE,
                               ) -> tuple[np.ndarray, np.ndarray, ExclusiveEnumerationStatus]:
    """Enumerate exclusive event probabilities under independence lazily (``Probability.IndependentExclusiveLazy``).

    Combinations are generated in subset-size order instead of read from a materialized indicator
    matrix.

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    include_no_event_row : bool, optional
        Emit the all-zero row, carrying ``prod(1 - p)``, first; default False.
    max_emitted_combinations : int, optional
        The cap on emitted combinations; 0 (default) for none.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances; default 1e-4.

    Returns
    -------
    event_probabilities : numpy.ndarray
        The exclusive event probabilities.
    event_indicators : numpy.ndarray
        One int32 indicator row per event.
    status : ExclusiveEnumerationStatus
        Complete, converged or capped.
    """
    events, rows = _event_lists()
    status = _probability().IndependentExclusiveLazy(to_net_doubles(probabilities), events, rows,
                                                     bool(include_no_event_row), int(max_emitted_combinations),
                                                     float(absolute_tolerance), float(relative_tolerance))
    return (*_events(events, rows), _status(status))


def positively_dependent_exclusive_lazy(probabilities: Any, absolute_tolerance: float = _TOLERANCE,
                                        relative_tolerance: float = _TOLERANCE,
                                        ) -> tuple[np.ndarray, np.ndarray, ExclusiveEnumerationStatus]:
    """Enumerate exclusive probabilities under perfect positive dependence (``PositivelyDependentExclusiveLazy``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances; default 1e-4.

    Returns
    -------
    event_probabilities : numpy.ndarray
        The exclusive event probabilities.
    event_indicators : numpy.ndarray
        One int32 indicator row per event.
    status : ExclusiveEnumerationStatus
        Complete or converged.
    """
    events, rows = _event_lists()
    status = _probability().PositivelyDependentExclusiveLazy(to_net_doubles(probabilities), events, rows,
                                                             float(absolute_tolerance), float(relative_tolerance))
    return (*_events(events, rows), _status(status))


def exclusive_pcm(probabilities: Any, correlation_matrix: Any, binomial_combinations: Any = None,
                  indicators: Any = None, absolute_tolerance: float = _TOLERANCE,
                  relative_tolerance: float = _TOLERANCE, return_events: bool = False) -> Any:
    """Return mutually exclusive event probabilities with PCM dependence (``Probability.ExclusivePCM``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    correlation_matrix : array_like
        The correlation matrix.
    binomial_combinations : array_like, optional
        The number of combinations of each subset size; given together with ``indicators``.
    indicators : array_like, optional
        The 2-D indicator rows in subset-size order.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances of the expansion (``return_events`` only); default 1e-4.
    return_events : bool, optional
        Run the convergent expansion and return the event lists; needs ``binomial_combinations``
        and ``indicators``. Default False.

    Returns
    -------
    numpy.ndarray or tuple of numpy.ndarray
        One exclusive probability per combination, or the event probabilities and their int32
        indicator rows.
    """
    p = _probability()
    values, matrix = to_net_doubles(probabilities), to_net_matrix(correlation_matrix)
    if binomial_combinations is None and indicators is None:
        if return_events:
            raise ValueError("return_events needs binomial_combinations and indicators.")
        return to_numpy_doubles(p.ExclusivePCM(values, matrix))
    combinations, rows = to_net_ints(binomial_combinations), _int_matrix(indicators)
    if not return_events:
        return to_numpy_doubles(p.ExclusivePCM(values, combinations, rows, matrix))
    events, event_rows = p.ExclusivePCM(values, combinations, rows, matrix, None, None, float(absolute_tolerance),
                                        float(relative_tolerance))
    return _events(events, event_rows)


def exclusive_pcm_lazy(probabilities: Any, correlation_matrix: Any, absolute_tolerance: float = _TOLERANCE,
                       relative_tolerance: float = _TOLERANCE,
                       ) -> tuple[np.ndarray, np.ndarray, ExclusiveEnumerationStatus]:
    """Enumerate exclusive event probabilities with PCM dependence lazily (``Probability.ExclusivePCMLazy``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    correlation_matrix : array_like
        The correlation matrix.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances; default 1e-4.

    Returns
    -------
    event_probabilities : numpy.ndarray
        The exclusive event probabilities.
    event_indicators : numpy.ndarray
        One int32 indicator row per event.
    status : ExclusiveEnumerationStatus
        Complete or converged.
    """
    events, rows = _event_lists()
    status = _probability().ExclusivePCMLazy(to_net_doubles(probabilities), to_net_matrix(correlation_matrix),
                                             events, rows, float(absolute_tolerance), float(relative_tolerance))
    return (*_events(events, rows), _status(status))


def exclusive_mvn(probabilities: Any, multivariate_normal: Any, binomial_combinations: Any = None,
                  indicators: Any = None, absolute_tolerance: float = _TOLERANCE,
                  relative_tolerance: float = _TOLERANCE, return_events: bool = False) -> Any:
    """Return mutually exclusive event probabilities with multivariate normal dependence (``Probability.ExclusiveMVN``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    multivariate_normal : MultivariateNormal
        The standard multivariate normal of the underlying variates.
    binomial_combinations : array_like, optional
        The number of combinations of each subset size; given together with ``indicators``.
    indicators : array_like, optional
        The 2-D indicator rows in subset-size order.
    absolute_tolerance, relative_tolerance : float, optional
        The convergence tolerances of the expansion (``return_events`` only); default 1e-4.
    return_events : bool, optional
        Run the convergent expansion and return the event lists; default False.

    Returns
    -------
    numpy.ndarray or tuple of numpy.ndarray
        As for :func:`exclusive_pcm`.
    """
    p = _probability()
    values, mvn = to_net_doubles(probabilities), _mvn(multivariate_normal)
    if binomial_combinations is None and indicators is None:
        if return_events:
            raise ValueError("return_events needs binomial_combinations and indicators.")
        return to_numpy_doubles(p.ExclusiveMVN(values, mvn))
    combinations, rows = to_net_ints(binomial_combinations), _int_matrix(indicators)
    if not return_events:
        return to_numpy_doubles(p.ExclusiveMVN(values, combinations, rows, mvn))
    events, event_rows = p.ExclusiveMVN(values, combinations, rows, mvn, None, None, float(absolute_tolerance),
                                        float(relative_tolerance))
    return _events(events, event_rows)


def common_cause_adjustment(probabilities: Any, correlation_matrix: Any = None,
                            dependency: DependencyType | str = DependencyType.CORRELATION_MATRIX) -> float:
    """Return the common cause adjustment factor (``Probability.CommonCauseAdjustment``).

    ``(1 - P(no event)) / sum p``: the union of the events under the dependence model over the sum
    of their probabilities. Without a correlation matrix the events are independent.

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.
    correlation_matrix : array_like, optional
        The correlation matrix for the correlation-matrix dependence.
    dependency : DependencyType or str, optional
        The dependence model with a correlation matrix; default the correlation matrix.

    Returns
    -------
    float
        The adjustment factor; 1 for a single event or zero total.
    """
    p = _probability()
    values = to_net_doubles(probabilities)
    if correlation_matrix is None:
        return float(p.CommonCauseAdjustment(values))
    return float(p.CommonCauseAdjustment(values, to_net_matrix(correlation_matrix), _dependency(dependency)))


def mutually_exclusive_adjustment(probabilities: Any) -> float:
    """Return the mutually exclusive adjustment factor (``Probability.MutuallyExclusiveAdjustment``).

    Parameters
    ----------
    probabilities : array_like
        The marginal probabilities.

    Returns
    -------
    float
        The adjustment factor.
    """
    return float(_probability().MutuallyExclusiveAdjustment(to_net_doubles(probabilities)))


__all__ = [
    "DependencyType", "ExclusiveEnumerationStatus", "a_and_b", "a_given_b", "a_not_b", "a_or_b", "b_given_a",
    "b_not_a", "common_cause_adjustment", "exclusive_mvn", "exclusive_pcm", "exclusive_pcm_lazy",
    "independent_exclusive", "independent_exclusive_lazy", "independent_joint_probability", "independent_union",
    "joint_probabilities_mvn", "joint_probabilities_pcm", "joint_probability", "joint_probability_hpcm",
    "joint_probability_mvn", "joint_probability_pcm", "mutually_exclusive_adjustment", "negative_joint_probability",
    "negatively_dependent_union", "positive_joint_probability", "positively_dependent_exclusive",
    "positively_dependent_exclusive_lazy", "positively_dependent_union", "single_factor_conditional_probabilities",
    "union", "union_mvn", "union_pcm", "union_pcm_lazy", "union_single_factor",
]

translate_public(globals())
