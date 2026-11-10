"""Finite-difference derivatives from ``Numerics.Mathematics.NumericalDerivative``.

Functions are Python callables, C# source or .NET delegates (see :mod:`pynumerics.callbacks`).
A ``step_size`` of None lets Numerics choose the step (:func:`calculate_step_size`); the other
defaults are the Numerics defaults.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import namespace
from pynumerics.errors import translate_public
from pynumerics.marshal import _optional_doubles, to_net_doubles, to_numpy_doubles, to_numpy_matrix


def _derivative() -> Any:
    """Return ``Numerics.Mathematics.NumericalDerivative`` after starting the runtime."""
    return namespace("Numerics.Mathematics").NumericalDerivative


def _step(step_size: float | None) -> float:
    """Return the Numerics step argument: -1 (choose the step) for None."""
    return -1.0 if step_size is None else float(step_size)


def _scalar(name: str, f: Any, point: float, step_size: float | None) -> float:
    """Call the single-variable member ``name`` of ``NumericalDerivative``."""
    return float(getattr(_derivative(), name)(callbacks.scalar_function(f), float(point), _step(step_size)))


def forward_difference(f: Any, point: float, step_size: float | None = None) -> float:
    """Return f'(point) by the forward difference [f(x + h) − f(x)] / h (``ForwardDifference``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    point : float
        Where to differentiate.
    step_size : float, optional
        The step h; None lets Numerics choose it.

    Returns
    -------
    float
        The derivative estimate.
    """
    return _scalar("ForwardDifference", f, point, step_size)


def backward_difference(f: Any, point: float, step_size: float | None = None) -> float:
    """Return f'(point) by the backward difference [f(x) − f(x − h)] / h (``BackwardDifference``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    point : float
        Where to differentiate.
    step_size : float, optional
        The step h; None lets Numerics choose it.

    Returns
    -------
    float
        The derivative estimate.
    """
    return _scalar("BackwardDifference", f, point, step_size)


def central_difference(f: Any, point: float, step_size: float | None = None) -> float:
    """Return f'(point) by the central difference [f(x + h) − f(x − h)] / 2h (``CentralDifference``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    point : float
        Where to differentiate.
    step_size : float, optional
        The step h; None lets Numerics choose it.

    Returns
    -------
    float
        The derivative estimate.
    """
    return _scalar("CentralDifference", f, point, step_size)


def derivative(f: Any, point: float, step_size: float | None = None) -> float:
    """Return f'(point) by the Numerics default rule, the central difference (``Derivative``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    point : float
        Where to differentiate.
    step_size : float, optional
        The step h; None lets Numerics choose it.

    Returns
    -------
    float
        The derivative estimate.
    """
    return _scalar("Derivative", f, point, step_size)


def ridders(f: Any, point: float, step_size: float | None = None) -> tuple[float, float]:
    """Return f'(point) and its error estimate by Ridders' extrapolation (``RiddersMethod``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    point : float
        Where to differentiate.
    step_size : float, optional
        The initial step; None uses the large step 0.1 * (1 + abs(point)) that Numerical Recipes
        recommends, which the extrapolation refines to about 1e-15 relative for a smooth function.
        Pass a smaller step when ``f`` is not smooth over that scale (near a singularity or the
        edge of its domain). With the default step, a tableau that meets a NaN or infinite value
        (a logarithm at 0.05 that returns NaN for a negative argument, say) is rebuilt from half
        the step, up to 30 times, until it is finite.

    Returns
    -------
    tuple of float
        ``(derivative, error)``.
    """
    value, error = _derivative().RiddersMethod(callbacks.scalar_function(f), float(point), 0.0, _step(step_size))
    return float(value), float(error)


def second_derivative(f: Any, point: float, step_size: float | None = None) -> float:
    """Return f''(point) by the central second difference (``SecondDerivative``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    point : float
        Where to differentiate.
    step_size : float, optional
        The step h; None lets Numerics choose it.

    Returns
    -------
    float
        The second-derivative estimate.
    """
    return _scalar("SecondDerivative", f, point, step_size)


def second_derivative_forward(f: Any, point: float, step_size: float | None = None) -> float:
    """Return f''(point) by the forward second difference (``SecondDerivativeForward``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    point : float
        Where to differentiate.
    step_size : float, optional
        The step h; None lets Numerics choose it.

    Returns
    -------
    float
        The second-derivative estimate.
    """
    return _scalar("SecondDerivativeForward", f, point, step_size)


def second_derivative_backward(f: Any, point: float, step_size: float | None = None) -> float:
    """Return f''(point) by the backward second difference (``SecondDerivativeBackward``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    point : float
        Where to differentiate.
    step_size : float, optional
        The step h; None lets Numerics choose it.

    Returns
    -------
    float
        The second-derivative estimate.
    """
    return _scalar("SecondDerivativeBackward", f, point, step_size)


def jacobian_at(f: Any, x: Any, point: Any, step_size: float | None = None) -> np.ndarray:
    """Return J[i, j] = ∂f(x_i, θ)/∂θ_j at θ = point by central differences (``Jacobian(f, x, point)``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function ``f(x_i, theta)`` of a float and a float64 parameter array.
    x : array_like
        The fixed locations, one Jacobian row each.
    point : array_like
        The parameters θ, one Jacobian column each.
    step_size : float, optional
        The step h; None lets Numerics choose one per parameter.

    Returns
    -------
    numpy.ndarray
        The ``len(x)`` by ``len(point)`` Jacobian.
    """
    return to_numpy_matrix(_derivative().Jacobian(callbacks.to_delegate(f, "Func<double,double[],double>"),
                                                  to_net_doubles(x), to_net_doubles(point), _step(step_size)))


def jacobian(g: Any, theta: Any, lower_bounds: Any = None, upper_bounds: Any = None, rel_step: float = 1e-5,
             abs_step: float = 1e-7, max_backtrack: int = 5) -> np.ndarray:
    """Return the Jacobian of a vector function with bound-aware adaptive steps (``Jacobian(g, theta, ...)``).

    Central differences are used where both sides fit inside the bounds and one-sided differences
    near a bound; the step shrinks when an evaluation is not finite.

    Parameters
    ----------
    g : callable, str or System.Delegate
        Function of a float64 array returning a float64 array (``Func<double[], double[]>``).
    theta : array_like
        Where to differentiate.
    lower_bounds, upper_bounds : array_like, optional
        Parameter bounds; None for unbounded.
    rel_step : float, optional
        Relative step; default 1e-5.
    abs_step : float, optional
        Minimum absolute step; default 1e-7.
    max_backtrack : int, optional
        Step reductions allowed after a failed evaluation; default 5.

    Returns
    -------
    numpy.ndarray
        The ``len(g(theta))`` by ``len(theta)`` Jacobian.
    """
    return to_numpy_matrix(_derivative().Jacobian(
        callbacks.to_delegate(g, "Func<double[],double[]>"), to_net_doubles(theta), _optional_doubles(lower_bounds),
        _optional_doubles(upper_bounds), float(rel_step), float(abs_step), int(max_backtrack)))


def gradient(f: Any, theta: Any, lower_bounds: Any = None, upper_bounds: Any = None, rel_step: float = 1e-5,
             abs_step: float = 1e-7, max_backtrack: int = 5) -> np.ndarray:
    """Return the gradient of a scalar function with bound-aware adaptive steps (``Gradient``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of a float64 array returning a float.
    theta : array_like
        Where to differentiate.
    lower_bounds, upper_bounds : array_like, optional
        Parameter bounds; None for unbounded.
    rel_step : float, optional
        Relative step; default 1e-5.
    abs_step : float, optional
        Minimum absolute step; default 1e-7.
    max_backtrack : int, optional
        Step reductions allowed after a failed evaluation; default 5.

    Returns
    -------
    numpy.ndarray
        The gradient.

    Raises
    ------
    ArithmeticError
        If f(theta) is not finite (raised by Numerics).
    """
    return to_numpy_doubles(_derivative().Gradient(
        callbacks.vector_function(f), to_net_doubles(theta), _optional_doubles(lower_bounds),
        _optional_doubles(upper_bounds), float(rel_step), float(abs_step), int(max_backtrack)))


def hessian(f: Any, theta: Any, lower_bounds: Any = None, upper_bounds: Any = None, rel_step: float = 1e-4,
            abs_step: float = 1e-6, max_backtrack: int = 6) -> np.ndarray:
    """Return the symmetric Hessian of a scalar function with bound-aware adaptive steps (``Hessian``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of a float64 array returning a float.
    theta : array_like
        Where to differentiate.
    lower_bounds, upper_bounds : array_like, optional
        Parameter bounds; None for unbounded.
    rel_step : float, optional
        Relative step; default 1e-4.
    abs_step : float, optional
        Minimum absolute step; default 1e-6.
    max_backtrack : int, optional
        Step reductions allowed after a failed evaluation; default 6.

    Returns
    -------
    numpy.ndarray
        The square Hessian.

    Raises
    ------
    ArithmeticError
        If f(theta) is not finite (raised by Numerics).
    """
    return to_numpy_matrix(_derivative().Hessian(
        callbacks.vector_function(f), to_net_doubles(theta), _optional_doubles(lower_bounds),
        _optional_doubles(upper_bounds), float(rel_step), float(abs_step), int(max_backtrack)))


def calculate_step_size(x: float, order: int = 1) -> float:
    """Return the finite-difference step ε^(1/(1 + order))·(1 + abs(x)) Numerics uses (``CalculateStepSize``).

    Parameters
    ----------
    x : float
        The point.
    order : int, optional
        The derivative order; default 1.

    Returns
    -------
    float
        The step.
    """
    return float(_derivative().CalculateStepSize(float(x), int(order)))


__all__ = [
    "backward_difference", "calculate_step_size", "central_difference", "derivative", "forward_difference",
    "gradient", "hessian", "jacobian", "jacobian_at", "ridders", "second_derivative", "second_derivative_backward",
    "second_derivative_forward",
]

translate_public(globals())
