"""Special functions from ``Numerics.Mathematics.SpecialFunctions``.

Gamma, beta, error, Bessel and Debye functions, factorials and combinations, and polynomial
evaluation. Every value is computed by the Numerics routine, and arguments broadcast like NumPy
ufuncs. A scalar argument calls the scalar member, bound once per process. An array argument of a
gamma, beta, error, Bessel or Debye function crosses into .NET once, through the ``IList<double>``
overload of the member (one call per distinct order for the integer-order Bessel functions), and
every value is the bits of the scalar member.

The incomplete gamma functions take ``(a, x)`` in the order of ``scipy.special.gammainc``, and the
regularized incomplete beta function and its inverse take ``(a, b, x)`` in the order of
``scipy.special.betainc``.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import numpy as np

from pynumerics._util import bind, broadcast_loop, cached_overload, loop, namespace, overload
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _flat_doubles,
    _numpy_rectangular,
    to_net_doubles,
    to_net_ints,
    to_numpy_doubles,
    to_numpy_ints,
)


def _special() -> Any:
    """Return the ``Numerics.Mathematics.SpecialFunctions`` namespace after starting the runtime."""
    return namespace("Numerics.Mathematics.SpecialFunctions")


def _member(type_name: str, name: str, arity: int = 1) -> Any:
    """Return ``<type_name>.<name>`` taking ``arity`` doubles, bound once per process."""
    return cached_overload(f"{type_name}.{name}({arity} double)",
                           lambda: bind(getattr(getattr(_special(), type_name), name), arity))


def _typed_member(type_name: str, name: str, *types: str) -> Any:
    """Return the overload of ``<type_name>.<name>`` taking the named ``System`` types, bound once per process."""
    return cached_overload(f"{type_name}.{name}({', '.join(types)})",
                           lambda: overload(getattr(getattr(_special(), type_name), name), *types))


def _list_member(type_name: str, name: str, *types: str) -> Any:
    """Return the list overload of ``<type_name>.<name>`` taking ``types``, resolved once per process."""
    return cached_overload(f"{type_name}.{name}({', '.join(types)})",
                           lambda: overload(getattr(getattr(_special(), type_name), name), *types))


_APPLY_MEMBERS: dict[tuple[str, str, int], tuple[Any, Any]] = {}
"""The (scalar, list) overload pair of each ``_apply`` member, keyed by type, member name and arity."""


def _apply(type_name: str, name: str, *args: Any, order: tuple[int, ...] | None = None) -> Any:
    """Evaluate the double member ``<type_name>.<name>`` over the broadcast of ``args``.

    ``order`` gives the Numerics parameter order as indices into ``args`` (for example ``(1, 0)``
    when the member takes the arguments in the reverse order). Scalars call the bound scalar member
    and give a Python float. Arrays are broadcast, flattened and passed as one ``System.Double[]``
    per parameter to the ``IList<double>`` overload, so they cross once; the result has the
    broadcast shape and the bits of the scalar member at every element. A single argument needs no
    broadcasting, and each broadcast argument crosses in one strided copy, with no contiguous copy
    first.
    """
    members = _APPLY_MEMBERS.get((type_name, name, len(args)))
    if members is None:
        members = _APPLY_MEMBERS[(type_name, name, len(args))] = (
            _member(type_name, name, len(args)), _list_member(type_name, name, *(("IList[Double]",) * len(args))))
    scalar, vector = members
    order = order or range(len(args))
    if len(args) == 1:
        arrays = (np.asarray(args[0], dtype=np.float64),)
    else:
        arrays = np.broadcast_arrays(*(np.asarray(a, dtype=np.float64) for a in args))
    if arrays[0].ndim == 0:
        return float(scalar(*(float(arrays[i]) for i in order)))
    result = to_numpy_doubles(vector(*(_flat_doubles(arrays[i]) for i in order)))
    return result if arrays[0].ndim == 1 else result.reshape(arrays[0].shape)


def gamma_function(x: Any) -> Any:
    """Return the gamma function (``Gamma.Function``).

    Parameters
    ----------
    x : float or array_like
        Arguments; large positive arguments give +inf.

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``x``.

    Raises
    ------
    ArithmeticError
        At zero and at the negative integers (poles; raised by Numerics).
    """
    return _apply("Gamma", "Function", x)


def log_gamma(x: Any) -> Any:
    """Return ln|Gamma(x)| (``Gamma.LogGamma``), the quantity of ``scipy.special.gammaln``.

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``x``.
    """
    return _apply("Gamma", "LogGamma", x)


def digamma(x: Any) -> Any:
    """Return the digamma function, the derivative of ln Gamma (``Gamma.Digamma``, algorithm AS 103).

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``x``.
    """
    return _apply("Gamma", "Digamma", x)


def trigamma(x: Any) -> Any:
    """Return the trigamma function, the second derivative of ln Gamma (``Gamma.Trigamma``).

    Parameters
    ----------
    x : float or array_like
        Positive arguments.

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``x``.
    """
    return _apply("Gamma", "Trigamma", x)


def gamma_max() -> float:
    """Return ``Gamma.GammaMax``, the largest argument for which the gamma function is finite.

    Returns
    -------
    float
        The constant the loaded assembly uses.
    """
    return float(_special().Gamma.GammaMax)


def stirling_gamma(x: Any) -> Any:
    """Return Stirling's approximation of the gamma function (``Gamma.Stirling``).

    An asymptotic series, accurate for large arguments only.

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The approximations, shaped like ``x``.
    """
    return _apply("Gamma", "Stirling", x)


def lanczos_gamma(x: Any) -> Any:
    """Return the Lanczos approximation of the gamma function (``Gamma.Lanczos``).

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The approximations, shaped like ``x``.
    """
    return _apply("Gamma", "Lanczos", x)


def incomplete_gamma(a: Any, x: Any) -> Any:
    """Return the regularized lower incomplete gamma integral by algorithm AS 239 (``Gamma.Incomplete``).

    ``Gamma.Incomplete(X, alpha)`` is called with ``X = x`` and ``alpha = a``. It keeps the argument
    order and validation of AS 239 and returns the value of :func:`lower_incomplete_gamma`, accurate
    to about 1e-14 relative and nonzero wherever the result is representable.

    Parameters
    ----------
    a : float or array_like
        Positive shape parameters.
    x : float or array_like
        Non-negative upper limits; arrays broadcast with ``a``.

    Returns
    -------
    float or numpy.ndarray
        P(a, x), in the broadcast shape.

    Raises
    ------
    ValueError
        If ``a <= 0`` or ``x < 0`` (raised by Numerics).
    """
    return _apply("Gamma", "Incomplete", a, x, order=(1, 0))


def lower_incomplete_gamma(a: Any, x: Any) -> Any:
    """Return the regularized lower incomplete gamma function P(a, x) (``Gamma.LowerIncomplete``).

    The quantity of ``scipy.special.gammainc``; P(a, x) = 0 for x <= 0. Series, continued fraction
    and Temme's uniform expansion, accurate to about 1e-14 relative where P exceeds 1e-30.

    Parameters
    ----------
    a : float or array_like
        Shape parameters.
    x : float or array_like
        Upper limits; arrays broadcast with ``a``.

    Returns
    -------
    float or numpy.ndarray
        P(a, x), in the broadcast shape.
    """
    return _apply("Gamma", "LowerIncomplete", a, x)


def upper_incomplete_gamma(a: Any, x: Any) -> Any:
    """Return the regularized upper incomplete gamma function Q(a, x) (``Gamma.UpperIncomplete``).

    The quantity of ``scipy.special.gammaincc``; Q(a, x) = 1 for x <= 0. Accurate to about 1e-14
    relative where Q exceeds 1e-30.

    Parameters
    ----------
    a : float or array_like
        Shape parameters.
    x : float or array_like
        Lower limits; arrays broadcast with ``a``.

    Returns
    -------
    float or numpy.ndarray
        Q(a, x), in the broadcast shape.
    """
    return _apply("Gamma", "UpperIncomplete", a, x)


def inverse_lower_incomplete_gamma(a: Any, y: Any) -> Any:
    """Return x with P(a, x) = y (``Gamma.InverseLowerIncomplete``), as ``scipy.special.gammaincinv``.

    Parameters
    ----------
    a : float or array_like
        Positive shape parameters.
    y : float or array_like
        Probabilities in [0, 1]; arrays broadcast with ``a``.

    Returns
    -------
    float or numpy.ndarray
        The quantiles (0 at y = 0, +inf at y = 1), in the broadcast shape.
    """
    return _apply("Gamma", "InverseLowerIncomplete", a, y)


def inverse_upper_incomplete_gamma(a: Any, y: Any) -> Any:
    """Return x with Q(a, x) = y (``Gamma.InverseUpperIncomplete``), as ``scipy.special.gammainccinv``.

    Parameters
    ----------
    a : float or array_like
        Positive shape parameters.
    y : float or array_like
        Probabilities in [0, 1]; arrays broadcast with ``a``.

    Returns
    -------
    float or numpy.ndarray
        The quantiles, in the broadcast shape.
    """
    return _apply("Gamma", "InverseUpperIncomplete", a, y)


def beta_function(a: Any, b: Any) -> Any:
    """Return the beta function B(a, b) = Gamma(a) Gamma(b) / Gamma(a + b) (``Beta.Function``).

    Parameters
    ----------
    a, b : float or array_like
        Positive arguments; arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The values, in the broadcast shape.
    """
    return _apply("Beta", "Function", a, b)


def log_beta(a: Any, b: Any) -> Any:
    """Return the natural logarithm of the beta function, ln B(a, b) (``Beta.LogFunction``).

    The value comes from the routines BETALN, ALGDIV and BCORR of algorithm 708 (DiDonato and
    Morris, 1992), so it keeps its accuracy at large or very unequal shapes, where
    ``log_gamma(a) + log_gamma(b) - log_gamma(a + b)`` subtracts large, nearly equal numbers. It
    agrees with ``scipy.special.betaln``, including ln|B(a, b)| for non-positive arguments.

    Parameters
    ----------
    a, b : float or array_like
        Arguments; arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The values, in the broadcast shape.
    """
    return _apply("Beta", "LogFunction", a, b)


def beta_incomplete(a: Any, b: Any, x: Any) -> Any:
    """Return the regularized incomplete beta function I_x(a, b) (``Beta.Incomplete``).

    Parameters
    ----------
    a, b : float or array_like
        Positive shape parameters; arrays broadcast.
    x : float or array_like
        Arguments in [0, 1]; arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The values, in the broadcast shape.
    """
    return _apply("Beta", "Incomplete", a, b, x)


def beta_incomplete_inverse(a: Any, b: Any, probability: Any) -> Any:
    """Return the inverse of the regularized incomplete beta function (``Beta.IncompleteInverse``).

    Parameters
    ----------
    a, b : float or array_like
        Positive shape parameters; arrays broadcast.
    probability : float or array_like
        Probabilities in [0, 1]; arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The x with I_x(a, b) equal to each probability, in the broadcast shape.
    """
    return _apply("Beta", "IncompleteInverse", a, b, probability)


def beta_continued_fraction1(a: Any, b: Any, x: Any) -> Any:
    """Return continued fraction expansion 1 of the incomplete beta integral (``Beta.Incbcf``).

    One of the expansions ``beta_incomplete`` combines (Cephes ``incbcf``).

    Parameters
    ----------
    a, b : float or array_like
        Shape parameters.
    x : float or array_like
        Arguments in (0, 1); arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The values, in the broadcast shape.
    """
    return _apply("Beta", "Incbcf", a, b, x)


def beta_continued_fraction2(a: Any, b: Any, x: Any) -> Any:
    """Return continued fraction expansion 2 of the incomplete beta integral (``Beta.Incbd``).

    One of the expansions ``beta_incomplete`` combines (Cephes ``incbd``).

    Parameters
    ----------
    a, b : float or array_like
        Shape parameters.
    x : float or array_like
        Arguments in (0, 1); arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The values, in the broadcast shape.
    """
    return _apply("Beta", "Incbd", a, b, x)


def beta_power_series(a: Any, b: Any, x: Any) -> Any:
    """Return the power series of the regularized incomplete beta integral (``Beta.PowerSeries``).

    Accurate when b x is small and x is not close to 1.

    Parameters
    ----------
    a, b : float or array_like
        Shape parameters.
    x : float or array_like
        Arguments in [0, 1); arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The values, in the broadcast shape.
    """
    return _apply("Beta", "PowerSeries", a, b, x)


def beta_incomplete_ratio(a: Any, b: Any, x: Any, log_beta: Any) -> Any:
    """Return the incomplete beta function ratio I_x(a, b) by algorithm AS 63 (``Beta.IncompleteRatio``).

    ``Beta.IncompleteRatio(x, p, q, beta)`` is called with ``p = a``, ``q = b`` and
    ``beta = log_beta``.

    Parameters
    ----------
    a, b : float or array_like
        Positive shape parameters.
    x : float or array_like
        Arguments in [0, 1].
    log_beta : float or array_like
        The log beta function ln B(a, b), for example ``scipy.special.betaln(a, b)``; arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        I_x(a, b), in the broadcast shape.
    """
    return _apply("Beta", "IncompleteRatio", a, b, x, log_beta, order=(2, 0, 1, 3))


def erf(x: Any) -> Any:
    """Return the error function (``Erf.Function``).

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``x``.
    """
    return _apply("Erf", "Function", x)


def erfc(x: Any) -> Any:
    """Return the complementary error function 1 - erf(x) (``Erf.Erfc``).

    Evaluated directly, never as 1 - erf(x), from the scaled function exp(x**2) erfc(x) with the factor
    exp(-x**2) formed without rounding x**2; within about 5e-16 relative of mpmath up to x = 26.5, where the
    result becomes subnormal (scipy's erfc errs by up to about 6e-14 there).

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``x``.
    """
    return _apply("Erf", "Erfc", x)


def inverse_erf(y: Any) -> Any:
    """Return the inverse error function (``Erf.InverseErf``), as ``scipy.special.erfinv``.

    Parameters
    ----------
    y : float or array_like
        Values in [-1, 1]; y = +-1 gives +-inf.

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``y``.

    Raises
    ------
    ValueError
        If abs(y) > 1 (raised by Numerics).
    """
    return _apply("Erf", "InverseErf", y)


def inverse_erfc(y: Any) -> Any:
    """Return the inverse complementary error function (``Erf.InverseErfc``), as ``scipy.special.erfcinv``.

    Parameters
    ----------
    y : float or array_like
        Values in [0, 2].

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``y``.
    """
    return _apply("Erf", "InverseErfc", y)


def _bessel_n(name: str, n: Any, x: Any) -> Any:
    """Evaluate the integer-order Bessel member ``name`` over the broadcast of ``n`` and ``x``.

    Scalars call the bound scalar member. Arrays cross once per distinct order into the
    ``(int, IList<double>)`` overload, and every value is the scalar member's.
    """
    scalar = _typed_member("Bessel", name, "Int32", "Double")
    vector = _list_member("Bessel", name, "Int32", "IList[Double]")
    orders, values = np.broadcast_arrays(np.asarray(n, dtype=np.int64), np.asarray(x, dtype=np.float64))
    if values.ndim == 0:
        return float(scalar(int(orders), float(values)))
    out = np.empty(values.shape, dtype=np.float64)
    for order in np.unique(orders):
        mask = orders == order
        out[mask] = to_numpy_doubles(vector(int(order), to_net_doubles(values[mask])))
    return out


def bessel_i0(x: Any) -> Any:
    """Return the modified Bessel function of the first kind of order 0 (``Bessel.I0``).

    The Cephes Chebyshev expansion of exp(-abs(x)) I0(x) times exp(abs(x)), accurate to about 1e-15 relative.

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        I0(x), shaped like ``x``.
    """
    return _apply("Bessel", "I0", x)


def bessel_i1(x: Any) -> Any:
    """Return the modified Bessel function of the first kind of order 1 (``Bessel.I1``).

    The Cephes Chebyshev expansion of exp(-abs(x)) I1(x) times exp(abs(x)), accurate to about 1e-15 relative.

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        I1(x), shaped like ``x``.
    """
    return _apply("Bessel", "I1", x)


def bessel_i0e(x: Any) -> Any:
    """Return the exponentially scaled exp(-abs(x)) I0(x) (``Bessel.I0e``), as ``scipy.special.i0e``.

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The scaled values, shaped like ``x``.
    """
    return _apply("Bessel", "I0e", x)


def bessel_i1e(x: Any) -> Any:
    """Return the exponentially scaled exp(-abs(x)) I1(x) (``Bessel.I1e``), as ``scipy.special.i1e``.

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The scaled values, shaped like ``x``.
    """
    return _apply("Bessel", "I1e", x)


def bessel_in(n: Any, x: Any) -> Any:
    """Return the modified Bessel function of the first kind of integer order n (``Bessel.In``).

    Miller's downward recurrence normalized by I0.

    Parameters
    ----------
    n : int or array_like of int
        Orders.
    x : float or array_like
        Arguments; arrays broadcast with ``n``.

    Returns
    -------
    float or numpy.ndarray
        In(x), in the broadcast shape.
    """
    return _bessel_n("In", n, x)


def bessel_k0(x: Any) -> Any:
    """Return the modified Bessel function of the second kind of order 0 (``Bessel.K0``).

    Parameters
    ----------
    x : float or array_like
        Positive arguments.

    Returns
    -------
    float or numpy.ndarray
        K0(x), shaped like ``x``.
    """
    return _apply("Bessel", "K0", x)


def bessel_k1(x: Any) -> Any:
    """Return the modified Bessel function of the second kind of order 1 (``Bessel.K1``).

    Parameters
    ----------
    x : float or array_like
        Positive arguments.

    Returns
    -------
    float or numpy.ndarray
        K1(x), shaped like ``x``.
    """
    return _apply("Bessel", "K1", x)


def bessel_kn(n: Any, x: Any) -> Any:
    """Return the modified Bessel function of the second kind of integer order n (``Bessel.Kn``).

    Forward recurrence from K0 and K1.

    Parameters
    ----------
    n : int or array_like of int
        Orders.
    x : float or array_like
        Positive arguments; arrays broadcast with ``n``.

    Returns
    -------
    float or numpy.ndarray
        Kn(x), in the broadcast shape.
    """
    return _bessel_n("Kn", n, x)


def bessel_j0(x: Any) -> Any:
    """Return the Bessel function of the first kind of order 0 (``Bessel.J0``).

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        J0(x), shaped like ``x``.
    """
    return _apply("Bessel", "J0", x)


def bessel_j1(x: Any) -> Any:
    """Return the Bessel function of the first kind of order 1 (``Bessel.J1``).

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        J1(x), shaped like ``x``.
    """
    return _apply("Bessel", "J1", x)


def bessel_jn(n: Any, x: Any) -> Any:
    """Return the Bessel function of the first kind of integer order n (``Bessel.Jn``).

    Forward recurrence for abs(x) > n; otherwise the continued fraction for Jn/Jn-1 and the backward
    recurrence, normalized by J0 or J1.

    Parameters
    ----------
    n : int or array_like of int
        Orders.
    x : float or array_like
        Arguments; arrays broadcast with ``n``.

    Returns
    -------
    float or numpy.ndarray
        Jn(x), in the broadcast shape.
    """
    return _bessel_n("Jn", n, x)


def bessel_y0(x: Any) -> Any:
    """Return the Bessel function of the second kind of order 0 (``Bessel.Y0``).

    Parameters
    ----------
    x : float or array_like
        Positive arguments.

    Returns
    -------
    float or numpy.ndarray
        Y0(x), shaped like ``x``.
    """
    return _apply("Bessel", "Y0", x)


def bessel_y1(x: Any) -> Any:
    """Return the Bessel function of the second kind of order 1 (``Bessel.Y1``).

    Parameters
    ----------
    x : float or array_like
        Positive arguments.

    Returns
    -------
    float or numpy.ndarray
        Y1(x), shaped like ``x``.
    """
    return _apply("Bessel", "Y1", x)


def bessel_yn(n: Any, x: Any) -> Any:
    """Return the Bessel function of the second kind of integer order n (``Bessel.Yn``).

    Forward recurrence from Y0 and Y1.

    Parameters
    ----------
    n : int or array_like of int
        Orders.
    x : float or array_like
        Positive arguments; arrays broadcast with ``n``.

    Returns
    -------
    float or numpy.ndarray
        Yn(x), in the broadcast shape.
    """
    return _bessel_n("Yn", n, x)


def debye3(x: Any) -> Any:
    """Return the order-3 Debye function D3(x) = 3 / x**3 * int_0^x t**3 / (e**t - 1) dt (``Debye.Function``).

    Parameters
    ----------
    x : float or array_like
        Positive arguments.

    Returns
    -------
    float or numpy.ndarray
        D3(x), shaped like ``x``.
    """
    return _apply("Debye", "Function", x)


def debye1(x: Any) -> Any:
    """Return the order-1 Debye function D1(x) = 1 / x * int_0^x t / (e**t - 1) dt (``Debye.FunctionOrderOne``).

    Parameters
    ----------
    x : float or array_like
        Positive arguments.

    Returns
    -------
    float or numpy.ndarray
        D1(x), shaped like ``x``.
    """
    return _apply("Debye", "FunctionOrderOne", x)


def factorial(n: Any) -> Any:
    """Return n! as a double (``Factorial.Function``, the ``long`` overload).

    Exact through 22!; above 170 the result is +inf.

    Parameters
    ----------
    n : int or array_like of int
        Non-negative integers.

    Returns
    -------
    float or numpy.ndarray
        The factorials, shaped like ``n``.

    Raises
    ------
    ValueError
        If ``n`` is negative (raised by Numerics).
    """
    return broadcast_loop(_typed_member("Factorial", "Function", "Int64"), n, kinds=(int,))


def log_factorial(n: Any) -> Any:
    """Return ln(n!) (``Factorial.LogFactorial``).

    Parameters
    ----------
    n : int or array_like of int
        Non-negative integers.

    Returns
    -------
    float or numpy.ndarray
        The values, shaped like ``n``.
    """
    return broadcast_loop(_typed_member("Factorial", "LogFactorial", "Int32"), n, kinds=(int,))


def binomial_coefficient(n: Any, k: Any) -> Any:
    """Return the binomial coefficient n choose k (``Factorial.BinomialCoefficient``).

    Correctly rounded, and exact below 2**53; invalid arguments (n < 0, k < 0 or k > n) give 0.

    Parameters
    ----------
    n, k : int or array_like of int
        The counts; arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The coefficients, in the broadcast shape.
    """
    member = _typed_member("Factorial", "BinomialCoefficient", "Int32", "Int32")
    return broadcast_loop(member, n, k, kinds=(int, int))


def find_combinations(m: int, n: int) -> np.ndarray:
    """Return every increasing index combination of size m from n items (``Factorial.FindCombinations``).

    Parameters
    ----------
    m : int
        The combination size.
    n : int
        The number of items.

    Returns
    -------
    numpy.ndarray
        A ``(n choose m, m)`` int32 array, one combination per row in lexicographic order.
    """
    rows = [to_numpy_ints(row) for row in _special().Factorial.FindCombinations(int(m), int(n))]
    return np.vstack(rows) if rows else np.empty((0, int(m)), dtype=np.int32)


def all_combinations(n: int) -> np.ndarray:
    """Return every non-empty subset of n items as a 0/1 indicator row (``Factorial.AllCombinations``).

    Parameters
    ----------
    n : int
        The number of items.

    Returns
    -------
    numpy.ndarray
        A ``(2**n - 1, n)`` int32 array; rows ordered by subset size, then lexicographically.
    """
    table = _special().Factorial.AllCombinations(int(n))
    return _numpy_rectangular(table, np.int32)


def iter_all_combinations(n: int) -> Iterator[np.ndarray]:
    """Yield the indices of every non-empty subset of n items lazily (``Factorial.AllCombinationsLazy``).

    The order is that of :func:`all_combinations`, without forming the indicator table.

    Parameters
    ----------
    n : int
        The number of items.

    Yields
    ------
    numpy.ndarray
        The int32 indices of one subset.
    """
    for combination in _special().Factorial.AllCombinationsLazy(int(n)):
        yield to_numpy_ints(combination)


def next_combination(combination: Any, n: int) -> tuple[bool, np.ndarray]:
    """Advance a k-combination of [0, n) to its lexicographic successor (``Factorial.NextCombination``).

    Parameters
    ----------
    combination : array_like of int
        The current strictly increasing index tuple; it is not modified.
    n : int
        The number of items.

    Returns
    -------
    advanced : bool
        False when ``combination`` was the last one.
    successor : numpy.ndarray
        The int32 indices Numerics left in the array (the successor when ``advanced``).
    """
    net = to_net_ints(combination)
    advanced = bool(_special().Factorial.NextCombination(net, int(n)))
    return advanced, to_numpy_ints(net)


def _polynomial(name: str, coefficients: Any, x: Any, *extra: int) -> Any:
    """Evaluate ``Evaluate.<name>`` with the coefficients copied once and ``x`` element by element."""
    member = getattr(_special().Evaluate, name)
    net = to_net_doubles(coefficients)
    return loop(x, lambda v: member(net, v, *extra))


def polynomial(coefficients: Any, x: Any) -> Any:
    """Evaluate sum c[i] x**i with ``coefficients[0]`` the constant term (``Evaluate.Polynomial``).

    Parameters
    ----------
    coefficients : array_like
        The coefficients, constant term first.
    x : float or array_like
        Points.

    Returns
    -------
    float or numpy.ndarray
        The polynomial values, shaped like ``x``.
    """
    return _polynomial("Polynomial", coefficients, x)


def polynomial_reversed(coefficients: Any, x: Any, n: int | None = None) -> Any:
    """Evaluate a polynomial with the constant term last by Horner's rule (``Evaluate.PolynomialRev``).

    ``coefficients[0] * x**n + ... + coefficients[n]``, the convention of ``numpy.polyval``.

    Parameters
    ----------
    coefficients : array_like
        The coefficients, highest power first.
    x : float or array_like
        Points.
    n : int, optional
        Use only the first n + 1 coefficients; None (the Numerics default, -1) uses them all.

    Returns
    -------
    float or numpy.ndarray
        The polynomial values, shaped like ``x``.

    Raises
    ------
    ValueError
        If ``n`` exceeds the number of coefficients (raised by Numerics).
    """
    return _polynomial("PolynomialRev", coefficients, x, -1 if n is None else int(n))


def polynomial_reversed_monic(coefficients: Any, x: Any) -> Any:
    """Evaluate the monic polynomial x**N + c[0] x**(N-1) + ... + c[N-1] (``Evaluate.PolynomialRev_1``).

    Parameters
    ----------
    coefficients : array_like
        The N coefficients after the implicit leading 1, highest power first.
    x : float or array_like
        Points.

    Returns
    -------
    float or numpy.ndarray
        The polynomial values, shaped like ``x``.
    """
    return _polynomial("PolynomialRev_1", coefficients, x)


def chebyshev(coefficients: Any, x: Any, terms: int | None = None) -> Any:
    """Sum the Chebyshev series c[0]/2 + sum c[k] T_k(x/2) by the Clenshaw recurrence (``Evaluate.Chebyshev``).

    The argument is twice the Chebyshev variable, the convention of the Cephes routine ``chbevl``,
    so ``chebyshev(c, x)`` equals ``numpy.polynomial.chebyshev.chebval(x / 2, c)`` with ``c[0]`` halved.

    Parameters
    ----------
    coefficients : array_like
        The coefficients, order-zero term first.
    x : float or array_like
        Points in [-2, 2].
    terms : int, optional
        Sum only the first ``terms`` coefficients; None sums them all.

    Returns
    -------
    float or numpy.ndarray
        The series values, shaped like ``x``.

    Raises
    ------
    ValueError
        If ``terms`` is below one or exceeds the number of coefficients (raised by Numerics).
    """
    if terms is None:
        return _polynomial("Chebyshev", coefficients, x)
    return _polynomial("Chebyshev", coefficients, x, int(terms))


__all__ = [
    "all_combinations", "bessel_i0", "bessel_i0e", "bessel_i1", "bessel_i1e", "bessel_in", "bessel_j0", "bessel_j1",
    "bessel_jn", "bessel_k0", "bessel_k1", "bessel_kn", "bessel_y0", "bessel_y1", "bessel_yn",
    "beta_continued_fraction1", "beta_continued_fraction2", "beta_function", "beta_incomplete",
    "beta_incomplete_inverse", "beta_incomplete_ratio", "beta_power_series", "binomial_coefficient", "chebyshev",
    "debye1", "debye3",
    "digamma", "erf", "erfc", "factorial", "find_combinations", "gamma_function", "gamma_max", "incomplete_gamma",
    "inverse_erf", "inverse_erfc", "inverse_lower_incomplete_gamma", "inverse_upper_incomplete_gamma",
    "iter_all_combinations", "lanczos_gamma", "log_beta", "log_factorial", "log_gamma", "lower_incomplete_gamma",
    "next_combination", "polynomial", "polynomial_reversed", "polynomial_reversed_monic", "stirling_gamma", "trigamma",
    "upper_incomplete_gamma",
]

translate_public(globals())
