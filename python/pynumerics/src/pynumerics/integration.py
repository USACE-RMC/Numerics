"""Numerical integration from ``Numerics.Mathematics.Integration``.

The fixed rules of the static ``Integration`` class (:func:`gauss_legendre`, :func:`trapezoidal_rule`,
...) return a float. Every ``Integrator`` subclass has an ``integrate_*`` function that builds the
Numerics integrator, applies the settings given, calls ``Integrate`` and returns an
:class:`IntegrationResult` (a :class:`VegasResult` for Vegas). Settings left as None keep the
Numerics defaults.

Integrands are Python callables, C# source compiled by :func:`pynumerics.callbacks.compile_csharp`
or .NET delegates (see :mod:`pynumerics.callbacks`). A compiled integrand runs without crossing into
Python per evaluation.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import namespace, net_enum
from pynumerics.enums import NetEnum
from pynumerics.errors import translate_public
from pynumerics.marshal import to_net_doubles, to_numpy_matrix


class IntegrationStatus(NetEnum):
    """Integration outcomes (``Numerics.Mathematics.IntegrationStatus``)."""

    NONE = "None"
    SUCCESS = "Success"
    MAXIMUM_ITERATIONS_REACHED = "MaximumIterationsReached"
    MAXIMUM_FUNCTION_EVALUATIONS_REACHED = "MaximumFunctionEvaluationsReached"
    FAILURE = "Failure"


class AdaptiveGaussKronrodMethod(NetEnum):
    """Adaptive strategies of ``AdaptiveGaussKronrod``.

    Wraps ``Numerics.Mathematics.Integration.AdaptiveGaussKronrodMethod``.

    ``LOCAL_BISECTION`` is the default recursive local bisection. ``QAGS`` is the QUADPACK QAGS algorithm
    behind ``scipy.integrate.quad``: global bisection of the largest error with epsilon-algorithm
    extrapolation, suited to endpoint singularities.
    """

    LOCAL_BISECTION = "LocalBisection"
    QAGS = "Qags"


INTEGRATOR_SETTINGS: Mapping[str, str] = MappingProxyType({
    "min_iterations": "MinIterations",
    "max_iterations": "MaxIterations",
    "min_function_evaluations": "MinFunctionEvaluations",
    "max_function_evaluations": "MaxFunctionEvaluations",
    "absolute_tolerance": "AbsoluteTolerance",
    "relative_tolerance": "RelativeTolerance",
    "report_failure": "ReportFailure",
})
"""The ``**settings`` keywords every ``integrate_*`` function takes, and the ``Integrator`` property each sets."""


@dataclass(frozen=True, slots=True)
class IntegrationResult:
    """The outcome of a Numerics integrator.

    Parameters
    ----------
    value : float
        The integral estimate (``Result``).
    standard_error : float
        The error estimate (``StandardError``); NaN for integrators without one.
    iterations : int
        Iterations used (``Iterations``).
    function_evaluations : int
        Integrand evaluations used (``FunctionEvaluations``).
    status : IntegrationStatus
        How the integration ended (``Status``).
    """

    value: float
    standard_error: float
    iterations: int
    function_evaluations: int
    status: IntegrationStatus


@dataclass(frozen=True, slots=True)
class VegasResult(IntegrationResult):
    """The outcome of a Vegas integration.

    Parameters
    ----------
    value : float
        The integral estimate (``Result``).
    standard_error : float
        The error estimate (``StandardError``).
    iterations : int
        Iterations used (``Iterations``).
    function_evaluations : int
        Integrand evaluations used (``FunctionEvaluations``).
    status : IntegrationStatus
        How the integration ended (``Status``).
    chi_squared : float
        The chi-squared per degree of freedom of the iteration estimates (``ChiSquared``).
    grid : numpy.ndarray
        The adapted importance-sampling grid, one row per dimension (``Grid``).
    """

    chi_squared: float
    grid: np.ndarray


def _integration() -> Any:
    """Return the ``Numerics.Mathematics.Integration`` namespace after starting the runtime."""
    return namespace("Numerics.Mathematics.Integration")


def _set(net: Any, **values: Any) -> None:
    """Set each .NET property named in ``values`` whose value is not None."""
    for name, value in values.items():
        if value is not None:
            setattr(net, name, value)


def _optional(kind: Callable[[Any], Any], value: Any) -> Any:
    """Return ``kind(value)``, or None when ``value`` is None."""
    return None if value is None else kind(value)


def _run(net: Any, settings: Mapping[str, Any], *functions: Any, bins: Any = None,
         standard_error: bool = True) -> IntegrationResult:
    """Apply the ``Integrator`` settings and the parallel policy, integrate and read the result.

    ``bins`` (a ``StratificationBins`` record or a .NET ``List<StratificationBin>``) selects the
    ``Integrate(List<StratificationBin>)`` overload.
    """
    unknown = sorted(set(settings) - set(INTEGRATOR_SETTINGS))
    if unknown:
        raise TypeError(f"Unknown integrator settings: {', '.join(unknown)}; known: {', '.join(INTEGRATOR_SETTINGS)}.")
    for key, value in settings.items():
        kind = bool if key == "report_failure" else float if key.endswith("tolerance") else int
        _set(net, **{INTEGRATOR_SETTINGS[key]: _optional(kind, value)})
    callbacks.apply_parallel_policy(net, *functions)
    if bins is None:
        net.Integrate()
    else:
        net.Integrate(getattr(bins, "net", bins))
    return IntegrationResult(float(net.Result), float(net.StandardError) if standard_error else float("nan"),
                             int(net.Iterations), int(net.FunctionEvaluations), IntegrationStatus(str(net.Status)))


def _random(seed: int | None) -> Any:
    """Return a seeded ``MersenneTwister``, or None to keep the integrator's own generator."""
    return None if seed is None else namespace("Numerics.Sampling").MersenneTwister(int(seed))


def _bounds(lower: Sequence[float], upper: Sequence[float]) -> tuple[int, Any, Any]:
    """Return the dimension count and the bounds as .NET arrays."""
    low, high = to_net_doubles(lower), to_net_doubles(upper)
    return int(low.Length), low, high


def gauss_legendre(f: Any, a: float, b: float) -> float:
    """Integrate ``f`` on [a, b] with the 10-point Gauss-Legendre rule (``Integration.GaussLegendre``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float (see :mod:`pynumerics.callbacks`).
    a, b : float
        The limits.

    Returns
    -------
    float
        The estimate, exact for polynomials up to degree 19.
    """
    return float(_integration().Integration.GaussLegendre(callbacks.scalar_function(f), float(a), float(b)))


def gauss_legendre_20(f: Any, a: float, b: float) -> float:
    """Integrate ``f`` on [a, b] with the 20-point Gauss-Legendre rule (``Integration.GaussLegendre20``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        The limits.

    Returns
    -------
    float
        The estimate, exact for polynomials up to degree 39.
    """
    return float(_integration().Integration.GaussLegendre20(callbacks.scalar_function(f), float(a), float(b)))


def trapezoidal_rule(f: Any, a: float, b: float, steps: int = 2) -> float:
    """Integrate ``f`` on [a, b] with the composite trapezoidal rule (``Integration.TrapezoidalRule``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        The limits.
    steps : int, optional
        Number of equal panels; default 2.

    Returns
    -------
    float
        The estimate.
    """
    return float(_integration().Integration.TrapezoidalRule(callbacks.scalar_function(f), float(a), float(b),
                                                            int(steps)))


def simpsons_rule(f: Any, a: float, b: float, steps: int = 2) -> float:
    """Integrate ``f`` on [a, b] with the composite Simpson's rule (``Integration.SimpsonsRule``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        The limits.
    steps : int, optional
        Number of equal panels, each evaluated at its ends and midpoint; default 2.

    Returns
    -------
    float
        The estimate.
    """
    return float(_integration().Integration.SimpsonsRule(callbacks.scalar_function(f), float(a), float(b),
                                                         int(steps)))


def midpoint_rule(f: Any, a: float, b: float, steps: int = 2) -> float:
    """Integrate ``f`` on [a, b] with the composite midpoint rule (``Integration.Midpoint``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        The limits.
    steps : int, optional
        Number of equal panels; default 2.

    Returns
    -------
    float
        The estimate.
    """
    return float(_integration().Integration.Midpoint(callbacks.scalar_function(f), float(a), float(b), int(steps)))


def integrate_trapezoidal(f: Any, a: float, b: float, **settings: Any) -> IntegrationResult:
    """Integrate ``f`` on [a, b] by trapezoidal refinement to tolerance (``TrapezoidalRule``).

    Each iteration doubles the panels until successive estimates agree within the tolerances.

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        Finite limits with ``a < b``.
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, counts and status; ``standard_error`` is NaN.

    Raises
    ------
    ValueError
        If a limit or setting is out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.scalar_function(f)
    return _run(_integration().TrapezoidalRule(fn, float(a), float(b)), settings, fn, standard_error=False)


def integrate_simpson(f: Any, a: float, b: float, **settings: Any) -> IntegrationResult:
    """Integrate ``f`` on [a, b] by Simpson refinement to tolerance (``SimpsonsRule``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        Finite limits with ``a < b``.
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, counts and status; ``standard_error`` is NaN.

    Raises
    ------
    ValueError
        If a limit or setting is out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.scalar_function(f)
    return _run(_integration().SimpsonsRule(fn, float(a), float(b)), settings, fn, standard_error=False)


def integrate_adaptive_simpson(f: Any, a: float, b: float, *, min_depth: int | None = None,
                               max_depth: int | None = None, bins: Any = None,
                               **settings: Any) -> IntegrationResult:
    """Integrate ``f`` on [a, b] with the adaptive Simpson's rule (``AdaptiveSimpsonsRule``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        Finite limits with ``a < b``.
    min_depth, max_depth : int, optional
        Recursion depth limits (``MinDepth``, ``MaxDepth``).
    bins : pynumerics.sampling.StratificationBins, optional
        Integrate over these bins instead of [a, b] (``Integrate(List<StratificationBin>)``).
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, standard error, counts and status.

    Raises
    ------
    ValueError
        If a limit or setting is out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.scalar_function(f)
    net = _integration().AdaptiveSimpsonsRule(fn, float(a), float(b))
    _set(net, MinDepth=_optional(int, min_depth), MaxDepth=_optional(int, max_depth))
    return _run(net, settings, fn, bins=bins)


def integrate_adaptive_gauss_kronrod(f: Any, a: float, b: float, *, min_depth: int | None = None,
                                     max_depth: int | None = None, recorder: Any = None,
                                     bins: Any = None, method: AdaptiveGaussKronrodMethod | str | None = None,
                                     **settings: Any) -> IntegrationResult:
    """Integrate ``f`` on [a, b] with the adaptive Gauss-Kronrod G10K21 rule (``AdaptiveGaussKronrod``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        Finite limits with ``a < b``.
    min_depth, max_depth : int, optional
        Recursion depth limits (``MinDepth``, ``MaxDepth``).
    recorder : callable, str or System.Delegate, optional
        Called as ``recorder(x, weight, fx)`` for each node of the accepted intervals (``Recorder``);
        the weights sum to the domain width and ``sum(weight * fx)`` reproduces the result.
    bins : pynumerics.sampling.StratificationBins, optional
        Integrate over these bins instead of [a, b] (``Integrate(List<StratificationBin>)``).
    method : AdaptiveGaussKronrodMethod or str, optional
        The adaptive strategy (``Method``); left at the Numerics default, local bisection, when None.
        ``"qags"`` selects QUADPACK QAGS, which suits endpoint singularities; there ``recorder``
        receives the final partition and ``relative_tolerance=0`` is legal with a positive
        ``absolute_tolerance``.
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, standard error, counts and status.

    Raises
    ------
    ValueError
        If a limit or setting is out of range (raised by Numerics), or ``method`` is not a method name.
    ArithmeticError
        If an accepted interval is not finite and ``report_failure`` is on (raised by Numerics).
    TypeError
        If a setting name is unknown.

    Examples
    --------
    >>> from pynumerics.integration import integrate_adaptive_gauss_kronrod
    >>> integrate_adaptive_gauss_kronrod("x => Math.Exp(-x * x)", -3.0, 3.0).value  # doctest: +SKIP
    1.772414696519...
    """
    fn = callbacks.scalar_function(f)
    record = None if recorder is None else callbacks.to_delegate(recorder, "Action<double,double,double>")
    integration = _integration()
    net = integration.AdaptiveGaussKronrod(fn, float(a), float(b))
    if method is not None:
        net.Method = net_enum(integration.AdaptiveGaussKronrodMethod, AdaptiveGaussKronrodMethod(method))
    _set(net, MinDepth=_optional(int, min_depth), MaxDepth=_optional(int, max_depth), Recorder=record)
    return _run(net, settings, fn, record, bins=bins)


def integrate_adaptive_gauss_lobatto(f: Any, a: float, b: float, **settings: Any) -> IntegrationResult:
    """Integrate ``f`` on [a, b] with the adaptive Gauss-Lobatto rule (``AdaptiveGaussLobatto``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of one float.
    a, b : float
        Finite limits with ``a < b``.
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, counts and status; ``standard_error`` is NaN.

    Raises
    ------
    ValueError
        If a limit or setting is out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.scalar_function(f)
    return _run(_integration().AdaptiveGaussLobatto(fn, float(a), float(b)), settings, fn, standard_error=False)


def integrate_adaptive_simpson_2d(f: Any, x_min: float, x_max: float, y_min: float, y_max: float, *,
                                  min_depth: int | None = None, max_depth: int | None = None,
                                  **settings: Any) -> IntegrationResult:
    """Integrate ``f(x, y)`` over a rectangle with the adaptive 2-D Simpson's rule (``AdaptiveSimpsonsRule2D``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of two floats.
    x_min, x_max, y_min, y_max : float
        The rectangle.
    min_depth, max_depth : int, optional
        Recursion depth limits (``MinDepth``, ``MaxDepth``).
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, standard error, counts and status.

    Raises
    ------
    ValueError
        If a limit or setting is out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.bivariate_function(f)
    net = _integration().AdaptiveSimpsonsRule2D(fn, float(x_min), float(x_max), float(y_min), float(y_max))
    _set(net, MinDepth=_optional(int, min_depth), MaxDepth=_optional(int, max_depth))
    return _run(net, settings, fn)


def integrate_adaptive_gauss_kronrod_2d(f: Any, x_min: float, x_max: float, y_min: float, y_max: float, *,
                                        min_depth: int | None = None, max_depth: int | None = None,
                                        recorder: Any = None, **settings: Any) -> IntegrationResult:
    """Integrate ``f(x, y)`` over a rectangle with the adaptive 2-D Gauss-Kronrod rule (``AdaptiveGaussKronrod2D``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of two floats.
    x_min, x_max, y_min, y_max : float
        The rectangle.
    min_depth, max_depth : int, optional
        Recursion depth limits (``MinDepth``, ``MaxDepth``).
    recorder : callable, str or System.Delegate, optional
        Called as ``recorder(x, y, weight, fxy)`` for each node of the accepted cells (``Recorder``).
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, standard error, counts and status.

    Raises
    ------
    ValueError
        If a limit or setting is out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.bivariate_function(f)
    record = None if recorder is None else callbacks.to_delegate(recorder, "Action<double,double,double,double>")
    net = _integration().AdaptiveGaussKronrod2D(fn, float(x_min), float(x_max), float(y_min), float(y_max))
    _set(net, MinDepth=_optional(int, min_depth), MaxDepth=_optional(int, max_depth), Recorder=record)
    return _run(net, settings, fn, record)


def integrate_monte_carlo(f: Any, lower: Sequence[float], upper: Sequence[float], *, seed: int | None = None,
                          use_sobol_sequence: bool | None = None, **settings: Any) -> IntegrationResult:
    """Integrate ``f(x)`` over a box by plain Monte Carlo sampling (``MonteCarloIntegration``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of a float64 array of length ``len(lower)``.
    lower, upper : array_like
        The box bounds, one per dimension.
    seed : int, optional
        Seed of the ``MersenneTwister`` set as ``Random``; None keeps the Numerics generator.
    use_sobol_sequence : bool, optional
        Sample a Sobol sequence instead of ``Random`` (``UseSobolSequence``; Numerics default False).
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, standard error, counts and status.

    Raises
    ------
    ValueError
        If the bounds or a setting are out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.vector_function(f)
    dimensions, low, high = _bounds(lower, upper)
    net = _integration().MonteCarloIntegration(fn, dimensions, low, high)
    _set(net, Random=_random(seed), UseSobolSequence=_optional(bool, use_sobol_sequence))
    return _run(net, settings, fn)


def integrate_miser(f: Any, lower: Sequence[float], upper: Sequence[float], *, seed: int | None = None,
                    use_sobol_sequence: bool | None = None, fraction: float | None = None,
                    minimum_number_of_subregion_points: int | None = None,
                    minimum_number_of_bisections: int | None = None, dither: float | None = None,
                    **settings: Any) -> IntegrationResult:
    """Integrate ``f(x)`` over a box by recursive stratified sampling (``Miser``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of a float64 array of length ``len(lower)``.
    lower, upper : array_like
        The box bounds, one per dimension.
    seed : int, optional
        Seed of the ``MersenneTwister`` set as ``Random``; None keeps the Numerics generator.
    use_sobol_sequence : bool, optional
        Sample a Sobol sequence instead of ``Random`` (``UseSobolSequence``; Numerics default True).
    fraction : float, optional
        Fraction of the points used to explore variance (``Fraction``).
    minimum_number_of_subregion_points : int, optional
        Points below which a region is sampled directly (``MinimumNumberOfSubregionPoints``).
    minimum_number_of_bisections : int, optional
        Bisection threshold (``MinimumNumberOfBisections``).
    dither : float, optional
        Bisection dither (``Dither``).
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS`.

    Returns
    -------
    IntegrationResult
        The estimate, standard error, counts and status.

    Raises
    ------
    ValueError
        If the bounds or a setting are out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.vector_function(f)
    dimensions, low, high = _bounds(lower, upper)
    net = _integration().Miser(fn, dimensions, low, high)
    _set(net, Random=_random(seed), UseSobolSequence=_optional(bool, use_sobol_sequence),
         Fraction=_optional(float, fraction),
         MinimumNumberOfSubregionPoints=_optional(int, minimum_number_of_subregion_points),
         MinimumNumberOfBisections=_optional(int, minimum_number_of_bisections), Dither=_optional(float, dither))
    return _run(net, settings, fn)


def integrate_vegas(f: Any, lower: Sequence[float], upper: Sequence[float], *, seed: int | None = None,
                    use_sobol_sequence: bool | None = None, sobol_seed: int | None = None,
                    check_convergence: bool | None = None, independent_evaluations: int | None = None,
                    function_calls: int | None = None, alpha: float | None = None,
                    number_of_bins: int | None = None, tail_focus_parameter: float | None = None,
                    rare_event_probability: float | None = None, **settings: Any) -> VegasResult:
    """Integrate ``f(x, weight)`` over a box by adaptive importance sampling (``Vegas``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Integrand of a float64 array of length ``len(lower)`` and the sample weight
        (``Func<double[], double, double>``).
    lower, upper : array_like
        The box bounds, one per dimension.
    seed : int, optional
        Seed of the ``MersenneTwister`` set as ``Random``; None keeps the Numerics generator.
    use_sobol_sequence : bool, optional
        Sample a Sobol sequence instead of ``Random`` (``UseSobolSequence``; Numerics default True).
    sobol_seed : int, optional
        Seed of the Sobol sequence (``SobolSeed``).
    check_convergence : bool, optional
        Stop once the tolerances are met (``CheckConvergence``).
    independent_evaluations : int, optional
        Evaluations of the final independent estimate (``IndependentEvaluations``).
    function_calls : int, optional
        Evaluations per iteration (``FunctionCalls``).
    alpha : float, optional
        Grid adaptation rate (``Alpha``).
    number_of_bins : int, optional
        Grid bins per dimension (``NumberOfBins``).
    tail_focus_parameter : float, optional
        Power-transform exponent that concentrates samples near the upper bounds (``TailFocusParameter``).
    rare_event_probability : float, optional
        Tune the tail focus, bins and alpha for an event of this probability, applied after the
        other options (``ConfigureForRareEvents``).
    **settings
        ``Integrator`` settings, the keys of :data:`INTEGRATOR_SETTINGS` (Vegas sets
        ``RelativeTolerance`` to 1e-3 by default).

    Returns
    -------
    VegasResult
        The estimate, standard error, counts, status, chi-squared and adapted grid.

    Raises
    ------
    ValueError
        If the bounds or a setting are out of range (raised by Numerics).
    TypeError
        If a setting name is unknown.
    """
    fn = callbacks.to_delegate(f, "Func<double[],double,double>")
    dimensions, low, high = _bounds(lower, upper)
    net = _integration().Vegas(fn, dimensions, low, high)
    _set(net, Random=_random(seed), UseSobolSequence=_optional(bool, use_sobol_sequence),
         SobolSeed=_optional(int, sobol_seed), CheckConvergence=_optional(bool, check_convergence),
         IndependentEvaluations=_optional(int, independent_evaluations),
         FunctionCalls=_optional(int, function_calls), Alpha=_optional(float, alpha),
         NumberOfBins=_optional(int, number_of_bins), TailFocusParameter=_optional(float, tail_focus_parameter))
    if rare_event_probability is not None:
        net.ConfigureForRareEvents(float(rare_event_probability))
    base = _run(net, settings, fn)
    return VegasResult(base.value, base.standard_error, base.iterations, base.function_evaluations, base.status,
                       float(net.ChiSquared), to_numpy_matrix(net.Grid))


__all__ = [
    "INTEGRATOR_SETTINGS", "AdaptiveGaussKronrodMethod", "IntegrationResult", "IntegrationStatus", "VegasResult",
    "gauss_legendre",
    "gauss_legendre_20", "integrate_adaptive_gauss_kronrod", "integrate_adaptive_gauss_kronrod_2d",
    "integrate_adaptive_gauss_lobatto", "integrate_adaptive_simpson", "integrate_adaptive_simpson_2d",
    "integrate_miser", "integrate_monte_carlo", "integrate_simpson", "integrate_trapezoidal", "integrate_vegas",
    "midpoint_rule", "simpsons_rule", "trapezoidal_rule",
]

translate_public(globals())
