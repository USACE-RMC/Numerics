"""Scalar ODE solvers from ``Numerics.Mathematics.ODESolvers.RungeKutta``.

Each solver integrates dy/dt = f(t, y) for one state variable. ``f`` is a Python callable
``f(t, y)``, C# source or a ``Func<double, double, double>`` (see :mod:`pynumerics.callbacks`).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import namespace
from pynumerics.errors import translate_public
from pynumerics.marshal import to_numpy_doubles


def _runge_kutta() -> Any:
    """Return ``Numerics.Mathematics.ODESolvers.RungeKutta`` after starting the runtime."""
    return namespace("Numerics.Mathematics.ODESolvers").RungeKutta


def runge_kutta_second_order(f: Any, initial_value: float, start_time: float, end_time: float,
                             time_steps: int) -> np.ndarray:
    """Solve dy/dt = f(t, y) on an even grid with Heun's second-order method (``RungeKutta.SecondOrder``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        The right-hand side ``f(t, y)``.
    initial_value : float
        The value of y at ``start_time``.
    start_time, end_time : float
        The time span.
    time_steps : int
        Number of grid points, ends included.

    Returns
    -------
    numpy.ndarray
        The value of y at the ``time_steps`` equally spaced times from ``start_time`` to ``end_time``.
    """
    return to_numpy_doubles(_runge_kutta().SecondOrder(callbacks.bivariate_function(f), float(initial_value),
                                                       float(start_time), float(end_time), int(time_steps)))


def runge_kutta_fourth_order(f: Any, initial_value: float, start_time: float, end_time: float,
                             time_steps: int) -> np.ndarray:
    """Solve dy/dt = f(t, y) on an even grid with the classical fourth-order method (``RungeKutta.FourthOrder``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        The right-hand side ``f(t, y)``.
    initial_value : float
        The value of y at ``start_time``.
    start_time, end_time : float
        The time span.
    time_steps : int
        Number of grid points, ends included.

    Returns
    -------
    numpy.ndarray
        The value of y at the ``time_steps`` equally spaced times from ``start_time`` to ``end_time``.
    """
    return to_numpy_doubles(_runge_kutta().FourthOrder(callbacks.bivariate_function(f), float(initial_value),
                                                       float(start_time), float(end_time), int(time_steps)))


def runge_kutta_fourth_order_step(f: Any, initial_value: float, start_time: float, dt: float) -> float:
    """Advance dy/dt = f(t, y) one classical fourth-order step (``RungeKutta.FourthOrder(f, y0, t0, dt)``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        The right-hand side ``f(t, y)``.
    initial_value : float
        The value of y at ``start_time``.
    start_time : float
        The start time.
    dt : float
        The step.

    Returns
    -------
    float
        The value of y at ``start_time + dt``.
    """
    return float(_runge_kutta().FourthOrder(callbacks.bivariate_function(f), float(initial_value),
                                            float(start_time), float(dt)))


def runge_kutta_fehlberg(f: Any, initial_value: float, start_time: float, dt: float, dt_min: float,
                         tolerance: float = 1e-3) -> float:
    """Advance dy/dt = f(t, y) over ``dt`` with adaptive Runge-Kutta-Fehlberg 4(5) steps (``RungeKutta.Fehlberg``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        The right-hand side ``f(t, y)``.
    initial_value : float
        The value of y at ``start_time``.
    start_time : float
        The start time.
    dt : float
        The span to advance, which is also the largest step.
    dt_min : float
        The smallest step.
    tolerance : float, optional
        Absolute error tolerance per step; default 1e-3.

    Returns
    -------
    float
        The value of y at ``start_time + dt``.
    """
    return float(_runge_kutta().Fehlberg(callbacks.bivariate_function(f), float(initial_value), float(start_time),
                                         float(dt), float(dt_min), float(tolerance)))


def runge_kutta_cash_karp(f: Any, initial_value: float, start_time: float, dt: float, dt_min: float,
                          tolerance: float = 1e-3) -> float:
    """Advance dy/dt = f(t, y) over ``dt`` with adaptive Cash-Karp 4(5) steps (``RungeKutta.CashKarp``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        The right-hand side ``f(t, y)``.
    initial_value : float
        The value of y at ``start_time``.
    start_time : float
        The start time.
    dt : float
        The span to advance, which is also the largest step.
    dt_min : float
        The smallest step.
    tolerance : float, optional
        Absolute error tolerance per step; default 1e-3.

    Returns
    -------
    float
        The value of y at ``start_time + dt``.
    """
    return float(_runge_kutta().CashKarp(callbacks.bivariate_function(f), float(initial_value), float(start_time),
                                         float(dt), float(dt_min), float(tolerance)))


__all__ = [
    "runge_kutta_cash_karp", "runge_kutta_fehlberg", "runge_kutta_fourth_order", "runge_kutta_fourth_order_step",
    "runge_kutta_second_order",
]

translate_public(globals())
