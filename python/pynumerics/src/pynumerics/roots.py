"""Root finding from ``Numerics.Mathematics.RootFinding``: bisection, Brent, Newton-Raphson and secant.

Each function calls the static Numerics solver with the function given as a Python callable,
C# source or a .NET delegate (see :mod:`pynumerics.callbacks`). The defaults are the Numerics
defaults.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import namespace
from pynumerics.errors import translate_public
from pynumerics.marshal import to_net_doubles, to_numpy_doubles


def _root_finding() -> Any:
    """Return the ``Numerics.Mathematics.RootFinding`` namespace after starting the runtime."""
    return namespace("Numerics.Mathematics.RootFinding")


def bisection_solve(f: Any, first_guess: float, lower: float, upper: float, tolerance: float = 1e-8,
                    max_iterations: int = 1000, report_failure: bool = True) -> float:
    """Find a root of ``f`` in a sign-changing bracket by bisection (``Bisection.Solve``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    first_guess : float
        A point inside the bracket; it is validated but does not steer the search.
    lower, upper : float
        Bracket with f(lower)·f(upper) ≤ 0.
    tolerance : float, optional
        Tolerance on the step between successive estimates; default 1e-8.
    max_iterations : int, optional
        Iteration limit; default 1000.
    report_failure : bool, optional
        Raise when the solver fails to converge; default True.

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        If the bracket does not change sign or the guess lies outside it, or if the solver does not
        converge and ``report_failure`` is True (raised by Numerics).
    """
    return float(_root_finding().Bisection.Solve(callbacks.scalar_function(f), float(first_guess), float(lower),
                                                 float(upper), float(tolerance), int(max_iterations),
                                                 bool(report_failure)))


BRENT_MINIMUM_RELATIVE_TOLERANCE = 8.881784197001252e-16
"""``Brent.MinimumRelativeTolerance``: 4 * 2**-52, the default and smallest ``relative_tolerance`` of
:func:`brent_solve` (the default ``rtol`` of ``scipy.optimize.brentq``)."""


def brent_solve(f: Any, lower: float, upper: float, tolerance: float = 1e-8, max_iterations: int = 1000,
                report_failure: bool = True, *,
                relative_tolerance: float = BRENT_MINIMUM_RELATIVE_TOLERANCE) -> float:
    """Find a root of ``f`` in a sign-changing bracket by Brent's method (``Brent.Solve``).

    The solver is a literal port of SciPy's ``brentq.c``, so the root is bitwise equal to
    ``scipy.optimize.brentq(f, lower, upper, xtol=tolerance, rtol=relative_tolerance, maxiter=max_iterations)``.

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    lower, upper : float
        Bracket with f(lower)·f(upper) ≤ 0.
    tolerance : float, optional
        Absolute tolerance (``xtol``), non-negative; default 1e-8.
    max_iterations : int, optional
        Iteration limit; default 1000.
    report_failure : bool, optional
        Raise when the solver fails (the Numerics default); default True.
    relative_tolerance : float, optional
        Relative tolerance (``rtol``), at least :data:`BRENT_MINIMUM_RELATIVE_TOLERANCE`, the default.

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        If the bracket does not change sign, or a tolerance is out of range (raised by Numerics).
    """
    return float(_root_finding().Brent.Solve(callbacks.scalar_function(f), float(lower), float(upper),
                                             float(tolerance), float(relative_tolerance), int(max_iterations),
                                             bool(report_failure)))


def brent_bracket(f: Any, lower: float, upper: float,
                  max_iterations: int = 10) -> tuple[bool, float, float, float, float]:
    """Expand a range outward until it brackets a sign change (``Brent.Bracket``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    lower, upper : float
        Starting range.
    max_iterations : int, optional
        Expansion limit; default 10.

    Returns
    -------
    tuple
        ``(found, lower, upper, f_lower, f_upper)``.
    """
    found, lo, hi, f1, f2 = _root_finding().Brent.Bracket(callbacks.scalar_function(f), float(lower), float(upper),
                                                          0.0, 0.0, int(max_iterations))
    return bool(found), float(lo), float(hi), float(f1), float(f2)


def secant_solve(f: Any, lower: float, upper: float, tolerance: float = 1e-8, max_iterations: int = 1000,
                 report_failure: bool = True) -> float:
    """Find a root of ``f`` by the secant method started from two points (``Secant.Solve``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    lower, upper : float
        The two starting points, which bracket the root.
    tolerance : float, optional
        Tolerance on the step between successive estimates; default 1e-8.
    max_iterations : int, optional
        Iteration limit; default 1000.
    report_failure : bool, optional
        Raise when the solver fails to converge; default True.

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        If the solver does not converge and ``report_failure`` is True (raised by Numerics).
    """
    return float(_root_finding().Secant.Solve(callbacks.scalar_function(f), float(lower), float(upper),
                                              float(tolerance), int(max_iterations), bool(report_failure)))


def newton_raphson_solve(f: Any, df: Any, first_guess: float, tolerance: float = 1e-8, max_iterations: int = 1000,
                         report_failure: bool = True) -> float:
    """Find a root of ``f`` by Newton-Raphson iteration (``NewtonRaphson.Solve``).

    The iteration is unbounded; :func:`newton_raphson_robust_solve` keeps it inside a bracket.

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    df : callable, str or System.Delegate
        The derivative of ``f``.
    first_guess : float
        Starting point.
    tolerance : float, optional
        Tolerance on the step between successive estimates; default 1e-8.
    max_iterations : int, optional
        Iteration limit; default 1000.
    report_failure : bool, optional
        Raise when the solver fails to converge; default True.

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        If the solver does not converge and ``report_failure`` is True (raised by Numerics).
    """
    return float(_root_finding().NewtonRaphson.Solve(callbacks.scalar_function(f), callbacks.scalar_function(df),
                                                     float(first_guess), float(tolerance), int(max_iterations),
                                                     bool(report_failure)))


def newton_raphson_robust_solve(f: Any, df: Any, first_guess: float, lower: float, upper: float,
                                tolerance: float = 1e-8, max_iterations: int = 1000,
                                report_failure: bool = True) -> float:
    """Find a root of ``f`` by Newton-Raphson with bisection fallback in a bracket (``NewtonRaphson.RobustSolve``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of one float.
    df : callable, str or System.Delegate
        The derivative of ``f``.
    first_guess : float
        Starting point inside the bracket.
    lower, upper : float
        Bracket with f(lower)·f(upper) ≤ 0.
    tolerance : float, optional
        Tolerance on the step between successive estimates; default 1e-8.
    max_iterations : int, optional
        Iteration limit; default 1000.
    report_failure : bool, optional
        Raise when the solver fails to converge; default True.

    Returns
    -------
    float
        The root.

    Raises
    ------
    ValueError
        If the bracket does not change sign, or if the solver does not converge and ``report_failure``
        is True (raised by Numerics).
    """
    return float(_root_finding().NewtonRaphson.RobustSolve(
        callbacks.scalar_function(f), callbacks.scalar_function(df), float(first_guess), float(lower), float(upper),
        float(tolerance), int(max_iterations), bool(report_failure)))


def newton_raphson_solve_system(f: Any, jacobian: Any, first_guess: Any, tolerance: float = 1e-8,
                                max_iterations: int = 1000, report_failure: bool = True) -> np.ndarray:
    """Solve the system f(x) = 0 by multivariate Newton-Raphson (``NewtonRaphson.Solve(Func<Vector, Vector>, ...)``).

    Parameters
    ----------
    f : callable, str or System.Delegate
        Function of a float64 array returning an array of the same length (``Func<Vector, Vector>``).
    jacobian : callable, str or System.Delegate
        Function of a float64 array returning the square Jacobian matrix (``Func<Vector, Matrix>``).
    first_guess : array_like
        Starting point.
    tolerance : float, optional
        Tolerance on the step norm; default 1e-8.
    max_iterations : int, optional
        Iteration limit; default 1000.
    report_failure : bool, optional
        Raise when the solver fails to converge; default True.

    Returns
    -------
    numpy.ndarray
        The solution.

    Raises
    ------
    ValueError
        If the Jacobian is singular, or if the solver does not converge and ``report_failure`` is True
        (raised by Numerics).
    """
    linear_algebra = namespace("Numerics.Mathematics.LinearAlgebra")
    solution = _root_finding().NewtonRaphson.Solve(
        callbacks.to_delegate(f, "Func<Vector,Vector>"), callbacks.to_delegate(jacobian, "Func<Vector,Matrix>"),
        linear_algebra.Vector(to_net_doubles(first_guess)), float(tolerance), int(max_iterations),
        bool(report_failure))
    return to_numpy_doubles(solution.Array)


__all__ = [
    "BRENT_MINIMUM_RELATIVE_TOLERANCE", "bisection_solve", "brent_bracket", "brent_solve",
    "newton_raphson_robust_solve", "newton_raphson_solve", "newton_raphson_solve_system", "secant_solve",
]

translate_public(globals())
