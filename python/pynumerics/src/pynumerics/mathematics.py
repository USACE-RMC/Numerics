"""Numerical utilities of ``Numerics.Tools`` and fast Fourier transforms of ``Numerics.Mathematics.Fourier``.

``expm1`` and ``log1p`` call ``Tools.Expm1`` and ``Tools.Log1p``; an array crosses once through
their ``IList<double>`` overloads. The fast Fourier transforms (``fft``, ``real_fft``,
``bluestein_fft``, ``fft_convolve``, ``fft_correlation`` and ``fft_autocorrelation``) call
``Numerics.Mathematics.Fourier``, and the other utilities of ``Numerics.Tools`` (``log_sum_exp``,
distances, scaling, masked reductions, sequences, compression and the :data:`constants`) are
wrapped one to one. Special functions are in :mod:`pynumerics.special`, linear algebra in
:mod:`pynumerics.linalg`, root finding in :mod:`pynumerics.roots`, optimization in
:mod:`pynumerics.optimization` and link functions in :mod:`pynumerics.functions`.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics._util import (
    bind,
    bind_list,
    broadcast_loop,
    cached_overload,
    evaluate,
    loop,
    namespace,
    overload,
)
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _bytes_from_net,
    _bytes_to_net,
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_ints,
    to_numpy_matrix,
)


def _tools(name: str, x: Any) -> Any:
    """Evaluate ``Tools.Expm1`` or ``Tools.Log1p`` on a scalar or an array.

    The scalar overload is bound once per process; an array crosses once into the
    ``IList<double>`` overload. Every value is the scalar member's.
    """
    scalar = cached_overload(f"Tools.{name}(double)", lambda: bind(getattr(namespace("Numerics").Tools, name)))
    vector = cached_overload(f"Tools.{name}(IList<double>)",
                             lambda: bind_list(getattr(namespace("Numerics").Tools, name)))
    return evaluate(x, scalar, vector)


def expm1(x: Any) -> Any:
    """Return exp(x) − 1 accurately near zero (``Tools.Expm1``).

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The values exp(x) − 1, shaped like ``x``.
    """
    return _tools("Expm1", x)


def log1p(x: Any) -> Any:
    """Return log(1 + x) accurately near zero (``Tools.Log1p``).

    Parameters
    ----------
    x : float or array_like
        Arguments.

    Returns
    -------
    float or numpy.ndarray
        The values log(1 + x), shaped like ``x``.
    """
    return _tools("Log1p", x)


def _tools_member(name: str, *types: str) -> Any:
    """Return the ``Numerics.Tools`` overload of ``name`` taking ``types``, bound once per process."""
    return cached_overload(f"Tools.{name}({', '.join(types)})",
                           lambda: overload(getattr(namespace("Numerics").Tools, name), *types))


class NumericsConstants:
    """The numerical constants of ``Numerics.Tools``.

    Each property reads the .NET constant, so the values are the ones the loaded assembly uses.
    Use the module instance ``constants``.
    """

    @property
    def double_machine_epsilon(self) -> float:
        """``Tools.DoubleMachineEpsilon``, the unit roundoff 2**-53 as Numerics writes it."""
        return float(namespace("Numerics").Tools.DoubleMachineEpsilon)

    @property
    def single_machine_epsilon(self) -> float:
        """``Tools.SingleMachineEpsilon``, the single-precision unit roundoff 2**-24 widened from float32."""
        return float(namespace("Numerics").Tools.SingleMachineEpsilon)

    @property
    def euler(self) -> float:
        """``Tools.Euler``, the Euler-Mascheroni constant."""
        return float(namespace("Numerics").Tools.Euler)

    @property
    def golden_ratio(self) -> float:
        """``Tools.GoldenRatio``, (1 + sqrt(5)) / 2."""
        return float(namespace("Numerics").Tools.GoldenRatio)

    @property
    def log_sqrt_2pi(self) -> float:
        """``Tools.LogSqrt2PI``, log(sqrt(2 pi))."""
        return float(namespace("Numerics").Tools.LogSqrt2PI)

    @property
    def sqrt_2pi(self) -> float:
        """``Tools.Sqrt2PI``, sqrt(2 pi)."""
        return float(namespace("Numerics").Tools.Sqrt2PI)

    @property
    def sqrt2(self) -> float:
        """``Tools.Sqrt2``, sqrt(2)."""
        return float(namespace("Numerics").Tools.Sqrt2)

    @property
    def log2(self) -> float:
        """``Tools.Log2``, log(2)."""
        return float(namespace("Numerics").Tools.Log2)


constants = NumericsConstants()


def copy_sign(x: Any, y: Any) -> Any:
    """Return the magnitude of ``x`` with the sign of ``y`` (``Tools.CopySign``).

    Unlike ``math.copysign``, the sign test is ``y >= 0``, so ``y = -0.0`` gives a positive result.

    Parameters
    ----------
    x, y : float or array_like
        Magnitudes and sign sources; arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The values, in the broadcast shape.
    """
    return broadcast_loop(_tools_member("CopySign", "Double", "Double"), x, y)


def is_finite(x: float) -> bool:
    """Return whether ``x`` is neither NaN nor infinite (``Tools.IsFinite``).

    Parameters
    ----------
    x : float
        The value.

    Returns
    -------
    bool
        True for a finite value.
    """
    return bool(_tools_member("IsFinite", "Double")(float(x)))


def is_power_of_two(n: int) -> bool:
    """Return whether ``n`` is a positive power of 2 (``Tools.IsPowerOfTwo``).

    Parameters
    ----------
    n : int
        A 32-bit integer.

    Returns
    -------
    bool
        True when ``n`` is 1, 2, 4, ...
    """
    return bool(_tools_member("IsPowerOfTwo", "Int32")(int(n)))


def next_power_of_two(n: int) -> int:
    """Return the smallest power of 2 at least ``n``, or 1 for ``n < 1`` (``Tools.NextPowerOfTwo``).

    Parameters
    ----------
    n : int
        A 32-bit integer.

    Returns
    -------
    int
        The power of 2.
    """
    return int(_tools_member("NextPowerOfTwo", "Int32")(int(n)))


def clamp(x: Any, minimum: float, maximum: float) -> Any:
    """Clamp values to the inclusive range [minimum, maximum] (``Tools.Clamp``).

    Parameters
    ----------
    x : float or array_like
        The values.
    minimum, maximum : float
        The bounds.

    Returns
    -------
    float or numpy.ndarray
        The clamped values, shaped like ``x``; NaN passes through.
    """
    member = _tools_member("Clamp", "Double", "Double", "Double")
    lo, hi = float(minimum), float(maximum)
    return loop(x, lambda v: member(v, lo, hi))


def hypot(a: Any, b: Any) -> Any:
    """Return sqrt(a**2 + b**2) without destructive underflow or overflow (``Tools.Hypot``).

    Parameters
    ----------
    a, b : float or array_like
        The two values; arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The Euclidean norms, in the broadcast shape; within 2 ulps of the exact value. Special values follow
        IEEE 754 and ``numpy.hypot``: an infinite argument gives ``inf`` even when the other is NaN.
    """
    return broadcast_loop(_tools_member("Hypot", "Double", "Double"), a, b)


def sign(a: Any, b: Any) -> Any:
    """Return the magnitude of ``a`` with the sign of ``b`` (``Tools.Sign``, the Fortran ``SIGN``).

    Parameters
    ----------
    a, b : float or array_like
        Magnitudes and sign sources (``b >= 0`` counts as positive); arrays broadcast.

    Returns
    -------
    float or numpy.ndarray
        The values, in the broadcast shape.
    """
    return broadcast_loop(_tools_member("Sign", "Double", "Double"), a, b)


def sqr(a: Any) -> Any:
    """Return ``a * a`` (``Tools.Sqr``).

    Parameters
    ----------
    a : float or array_like
        The values.

    Returns
    -------
    float or numpy.ndarray
        The squares, shaped like ``a``.
    """
    return loop(a, _tools_member("Sqr", "Double"))


def power(a: Any, b: Any) -> Any:
    """Raise ``a`` to the integer power ``b`` by repeated squaring (``Tools.Pow``).

    Parameters
    ----------
    a : float or array_like
        Bases.
    b : int or array_like of int
        Integer exponents; arrays broadcast with ``a``.

    Returns
    -------
    float or numpy.ndarray
        The powers; ``0 ** -k`` is +inf.
    """
    return broadcast_loop(_tools_member("Pow", "Double", "Int32"), a, b, kinds=(float, int))


def log10(x: Any) -> Any:
    """Return the guarded base-10 logarithm (``Tools.Log10``).

    Values in [0, 1e-16) are raised to 1e-16, so ``log10(0)`` is -16 rather than -inf.

    Parameters
    ----------
    x : float or array_like
        The values.

    Returns
    -------
    float or numpy.ndarray
        The logarithms, shaped like ``x``.
    """
    return loop(x, _tools_member("Log10", "Double"))


def log(x: Any) -> Any:
    """Return the guarded natural logarithm (``Tools.Log``).

    Values in [0, 1e-16), including -0.0, are raised to 1e-16, so ``log(0)`` is about -36.84
    rather than -inf.

    Parameters
    ----------
    x : float or array_like
        The values.

    Returns
    -------
    float or numpy.ndarray
        The logarithms, shaped like ``x``.
    """
    return loop(x, _tools_member("Log", "Double"))


def point_distance(x1: float, y1: float, x2: float, y2: float) -> float:
    """Return the Euclidean distance between two points of the plane (``Tools.Distance`` of four doubles).

    Parameters
    ----------
    x1, y1 : float
        The first point.
    x2, y2 : float
        The second point.

    Returns
    -------
    float
        The distance sqrt((x2 - x1)**2 + (y2 - y1)**2).
    """
    member = _tools_member("Distance", "Double", "Double", "Double", "Double")
    return float(member(float(x1), float(y1), float(x2), float(y2)))


def distance(x: Any, y: Any) -> float:
    """Return the Euclidean distance ||x - y|| (``Tools.Distance`` of two lists).

    Parameters
    ----------
    x, y : array_like
        The points, one coordinate per element.

    Returns
    -------
    float
        The distance.
    """
    member = _tools_member("Distance", "IList[Double]", "IList[Double]")
    return float(member(to_net_doubles(x), to_net_doubles(y)))


def normalized_distance(x: Any, y: Any, lower: Any, upper: Any) -> float:
    """Return the Euclidean distance after min-max scaling each coordinate (``Tools.NormalizedDistance``).

    A coordinate with ``upper - lower <= 0`` (or NaN) contributes nothing.

    Parameters
    ----------
    x, y : array_like
        The points.
    lower, upper : array_like
        The bounds of each coordinate.

    Returns
    -------
    float
        The normalized distance.
    """
    tools = namespace("Numerics").Tools
    return float(tools.NormalizedDistance(to_net_doubles(x), to_net_doubles(y), to_net_doubles(lower),
                                          to_net_doubles(upper)))


def weighted_distance(x: Any, y: Any, weights: Any) -> float:
    """Return the weighted Euclidean distance sqrt(sum w (x - y)**2) (``Tools.WeightedDistance``).

    Parameters
    ----------
    x, y : array_like
        The points.
    weights : array_like
        One weight per coordinate.

    Returns
    -------
    float
        The weighted distance.
    """
    tools = namespace("Numerics").Tools
    return float(tools.WeightedDistance(to_net_doubles(x), to_net_doubles(y), to_net_doubles(weights)))


def normalize(values: Any) -> np.ndarray:
    """Scale values to [0, 1] by their minimum and maximum (``Tools.Normalize``).

    Empty input, input holding NaN, or constant input gives zeros.

    Parameters
    ----------
    values : array_like
        The values.

    Returns
    -------
    numpy.ndarray
        (values - min) / (max - min).
    """
    return to_numpy_doubles(namespace("Numerics").Tools.Normalize(to_net_doubles(values)))


def denormalize(values: Any, minimum: float, maximum: float) -> np.ndarray:
    """Invert :func:`normalize` with the original minimum and maximum (``Tools.Denormalize``).

    Parameters
    ----------
    values : array_like
        Normalized values.
    minimum, maximum : float
        The minimum and maximum of the original data.

    Returns
    -------
    numpy.ndarray
        The values times (maximum - minimum) plus minimum.
    """
    return to_numpy_doubles(namespace("Numerics").Tools.Denormalize(to_net_doubles(values), float(minimum),
                                                                     float(maximum)))


def standardize(values: Any) -> np.ndarray:
    """Return the z-scores (values - mean) / sd with the sample standard deviation (``Tools.Standardize``).

    A zero or NaN standard deviation gives zeros.

    Parameters
    ----------
    values : array_like
        The values.

    Returns
    -------
    numpy.ndarray
        The standardized values.
    """
    return to_numpy_doubles(namespace("Numerics").Tools.Standardize(to_net_doubles(values)))


def destandardize(values: Any, mean: float, standard_deviation: float) -> np.ndarray:
    """Invert :func:`standardize` with the original mean and standard deviation (``Tools.Destandardize``).

    Parameters
    ----------
    values : array_like
        Standardized values.
    mean, standard_deviation : float
        The moments of the original data.

    Returns
    -------
    numpy.ndarray
        The values times standard_deviation plus mean.
    """
    return to_numpy_doubles(namespace("Numerics").Tools.Destandardize(to_net_doubles(values), float(mean),
                                                                       float(standard_deviation)))


def _masked(name: str, values: Any, indicators: Any, use_complement: bool) -> float:
    """Call ``Tools.<name>`` on doubles, or on doubles with 0/1 indicators and the complement flag."""
    if indicators is None:
        return float(_tools_member(name, "IList[Double]")(to_net_doubles(values)))
    member = _tools_member(name, "IList[Double]", "IList[Int32]", "Boolean")
    return float(member(to_net_doubles(values), to_net_ints(indicators), bool(use_complement)))


def sum_values(values: Any, indicators: Any = None, use_complement: bool = False) -> float:
    """Return the sum of the values, optionally of those an indicator selects (``Tools.Sum``).

    Integer input without indicators is summed by the ``IList<int>`` overload.

    Parameters
    ----------
    values : array_like
        The values.
    indicators : array_like of int, optional
        0 or 1 per value; with ``use_complement`` False the values marked 1 are summed, with it True
        those marked 0. None sums every value.
    use_complement : bool, optional
        Select the values marked 0; default False.

    Returns
    -------
    float
        The sum.
    """
    if indicators is None and np.asarray(values).dtype.kind in "iub":
        return float(_tools_member("Sum", "IList[Int32]")(to_net_ints(values)))
    return _masked("Sum", values, indicators, use_complement)


def _compensated_sum_member(*types: str) -> Any:
    """Return the ``Tools.CompensatedSum`` overload taking ``types``, bound once per process."""
    return cached_overload(
        f"Tools.CompensatedSum({', '.join(types)})",
        lambda: overload(namespace("Numerics").Tools.CompensatedSum, *types))


def compensated_sum(values: Any, indicators: Any = None, use_complement: bool = False, *, indices: Any = None,
                    start: int = 0, count: Any = None, kahan: bool = False) -> float:
    """Return the compensated (Neumaier) sum of the values (``Tools.CompensatedSum``).

    Optionally of those an indicator selects, or of the gathered values ``values[indices[i]]``
    for i in [``start``, ``start + count``).

    Parameters
    ----------
    values : array_like
        The values.
    indicators : array_like of int, optional
        0 or 1 per value, as for :func:`sum_values`. None sums every value.
    use_complement : bool, optional
        Select the values marked 0; default False.
    indices : array_like of int, optional
        Zero-based positions in ``values`` to gather and sum, in this order; a position may repeat.
        Not combined with ``indicators``.
    start : int, optional
        The first entry of ``indices`` to use; default 0.
    count : int, optional
        The number of entries of ``indices`` to use; default all from ``start`` on.
    kahan : bool, optional
        Use classic Kahan summation instead of Neumaier's method; default False. Kahan's
        ``a = v - c; f = t + a; c = (f - t) - a; t = f`` is reproduced bit for bit, for reproducing
        results of the classic Kahan loop; Neumaier's method is the more accurate. A final infinite
        sum is returned as is, and any value added after the sum becomes infinite gives NaN.

    Returns
    -------
    float
        The compensated sum. As for :func:`sum_values`, an empty input (or an empty range of
        ``indices``) gives NaN, where :func:`math.fsum` gives 0, and indicators of a different
        length from the values give NaN. A NaN, or both infinities, give NaN and a single
        infinity is returned, the IEEE sum; :func:`math.fsum` raises ``ValueError`` for both
        infinities and ``OverflowError`` for an intermediate overflow, which returns the infinity here.

    Raises
    ------
    ValueError
        If both ``indicators`` and ``indices`` are given, or the range or a gathered index is out of bounds.

    Notes
    -----
    Neumaier's improvement of Kahan summation keeps a running correction of the exact rounding
    error of each addition, so the error is bounded by 2u|s| + O(n u²) Σ|x_i| and stays small when
    a term exceeds the running sum. The input crosses into Numerics once.
    """
    extra = ("Boolean",) if kahan else ()
    flag = (True,) if kahan else ()
    if indices is not None:
        if indicators is not None:
            raise ValueError("Give either indicators or indices, not both.")
        net_indices = to_net_ints(indices)
        take = int(net_indices.Length) - int(start) if count is None else int(count)
        member = _compensated_sum_member("IList[Double]", "IList[Int32]", "Int32", "Int32", *extra)
        return float(member(to_net_doubles(values), net_indices, int(start), take, *flag))
    if indicators is None:
        return float(_compensated_sum_member("IList[Double]", *extra)(to_net_doubles(values), *flag))
    member = _compensated_sum_member("IList[Double]", "IList[Int32]", "Boolean", *extra)
    return float(member(to_net_doubles(values), to_net_ints(indicators), bool(use_complement), *flag))


def compensated_cumulative_sum(values: Any, kahan: bool = False) -> np.ndarray:
    """Return the compensated (Neumaier) running sums of the values (``Tools.CompensatedCumulativeSum``).

    Parameters
    ----------
    values : array_like
        The values.
    kahan : bool, optional
        Use classic Kahan summation instead of Neumaier's method, as in :func:`compensated_sum`;
        entry i is then the Kahan running sum t after value i. Default False.

    Returns
    -------
    numpy.ndarray
        Entry i is :func:`compensated_sum` of the first i + 1 values, bit for bit; empty for empty input.
    """
    types = ("IList[Double]", "Boolean") if kahan else ("IList[Double]",)
    member = cached_overload(
        f"Tools.CompensatedCumulativeSum({', '.join(types)})",
        lambda: overload(namespace("Numerics").Tools.CompensatedCumulativeSum, *types))
    return to_numpy_doubles(member(to_net_doubles(values), True) if kahan else member(to_net_doubles(values)))


def sum_product(values1: Any, values2: Any) -> float:
    """Return sum(values1 * values2) (``Tools.SumProduct``).

    Parameters
    ----------
    values1, values2 : array_like
        Equal-length values.

    Returns
    -------
    float
        The sum of the products.
    """
    return float(namespace("Numerics").Tools.SumProduct(to_net_doubles(values1), to_net_doubles(values2)))


def mean_values(values: Any, indicators: Any = None, use_complement: bool = False) -> float:
    """Return the arithmetic mean, optionally of the values an indicator selects (``Tools.Mean``).

    Parameters
    ----------
    values : array_like
        The values.
    indicators : array_like of int, optional
        0 or 1 per value, as for :func:`sum_values`.
    use_complement : bool, optional
        Select the values marked 0; default False.

    Returns
    -------
    float
        The mean.
    """
    return _masked("Mean", values, indicators, use_complement)


def product(values: Any, indicators: Any = None, use_complement: bool = False) -> float:
    """Return the product, optionally of the values an indicator selects (``Tools.Product``).

    Parameters
    ----------
    values : array_like
        The values.
    indicators : array_like of int, optional
        0 or 1 per value, as for :func:`sum_values`.
    use_complement : bool, optional
        Select the values marked 0; default False.

    Returns
    -------
    float
        The product.
    """
    return _masked("Product", values, indicators, use_complement)


def min_value(values: Any, indicators: Any = None, use_complement: bool = False) -> float:
    """Return the smallest value, optionally of those an indicator selects (``Tools.Min``).

    Empty input or input holding NaN gives NaN.

    Parameters
    ----------
    values : array_like
        The values.
    indicators : array_like of int, optional
        0 or 1 per value, as for :func:`sum_values`.
    use_complement : bool, optional
        Select the values marked 0; default False.

    Returns
    -------
    float
        The minimum.
    """
    return _masked("Min", values, indicators, use_complement)


def max_value(values: Any, indicators: Any = None, use_complement: bool = False) -> float:
    """Return the largest value, optionally of those an indicator selects (``Tools.Max``).

    Empty input or input holding NaN gives NaN.

    Parameters
    ----------
    values : array_like
        The values.
    indicators : array_like of int, optional
        0 or 1 per value, as for :func:`sum_values`.
    use_complement : bool, optional
        Select the values marked 0; default False.

    Returns
    -------
    float
        The maximum.
    """
    return _masked("Max", values, indicators, use_complement)


def min_max(values: Any) -> tuple[float, float]:
    """Return the minimum and the maximum in one pass (``Tools.MinMax``).

    Parameters
    ----------
    values : array_like
        The values.

    Returns
    -------
    tuple of float
        ``(minimum, maximum)``; both NaN for empty input or input holding NaN.
    """
    low, high = namespace("Numerics").Tools.MinMax(to_net_doubles(values), 0.0, 0.0)
    return float(low), float(high)


def arg_min(values: Any) -> int:
    """Return the index of the smallest value, ignoring NaN (``Tools.ArgMin``).

    Parameters
    ----------
    values : array_like
        The values.

    Returns
    -------
    int
        The zero-based index of the first minimum, or -1 when every value is NaN or there are none.
    """
    return int(namespace("Numerics").Tools.ArgMin(to_net_doubles(values)))


def arg_max(values: Any) -> int:
    """Return the index of the largest value, ignoring NaN (``Tools.ArgMax``).

    Parameters
    ----------
    values : array_like
        The values.

    Returns
    -------
    int
        The zero-based index of the first maximum, or -1 when every value is NaN or there are none.
    """
    return int(namespace("Numerics").Tools.ArgMax(to_net_doubles(values)))


def log_sum_exp(values: Any, v: Any = None, *, axis: int | None = None) -> Any:
    """Return log(sum(exp(values))) without overflow (``Tools.LogSumExp``).

    With one argument the whole array is reduced in one crossing by the ``IList<double>``
    overload; with ``axis`` a two-dimensional array is reduced along that axis in one crossing by
    the ``(double[,], int)`` overload; with two arguments, log(exp(values) + exp(v)) is evaluated
    pairwise by the ``(double, double)`` overload, broadcasting arrays. The evaluation and the
    special values are those of ``scipy.special.logsumexp``.

    Parameters
    ----------
    values : array_like or float
        Log values to reduce, or the first log value of each pair.
    v : float or array_like, optional
        The second log value of each pair.
    axis : int, optional
        0 or 1 (or -2 or -1): reduce a two-dimensional ``values`` over its rows, giving one value
        per column, or over its columns, giving one value per row. None, the default, reduces
        every value.

    Returns
    -------
    float or numpy.ndarray
        The reduction, the reductions along ``axis``, or the pairwise values. An empty reduction
        and a reduction of -inf values give -inf, any +inf gives +inf, and NaN propagates.

    Raises
    ------
    ValueError
        If ``axis`` is given with ``v``, or ``values`` is not two-dimensional, or ``axis`` is not
        one of 0, 1, -2 and -1.
    TypeError
        If ``axis`` is a bool or not an integer.
    """
    if axis is not None:
        if v is not None:
            raise ValueError("Give either v or axis, not both.")
        if isinstance(axis, bool) or not isinstance(axis, (int, np.integer)):
            raise TypeError(f"axis must be an integer; got {type(axis).__name__}.")
        if axis not in (0, 1, -2, -1):
            raise ValueError(f"axis must be 0, 1, -2 or -1; got {axis}.")
        array = np.asarray(values, dtype=np.float64)
        if array.ndim != 2:
            raise ValueError(f"axis needs a 2-D array; got shape {array.shape}.")
        return to_numpy_doubles(namespace("Numerics").Tools.LogSumExp(to_net_matrix(array), int(axis) % 2))
    if v is None:
        return float(_tools_member("LogSumExp", "IList[Double]")(to_net_doubles(values)))
    return broadcast_loop(_tools_member("LogSumExp", "Double", "Double"), values, v)


def sequence(start: float, end: float, step: float = 1) -> np.ndarray:
    """Return start, start + step, ... up to and including ``end`` (``Tools.Sequence``).

    Integer arguments use the ``int`` overload and give int32 values; any float argument uses
    the ``double`` overload, which accumulates ``v += step``. A zero step gives an empty array.

    Parameters
    ----------
    start, end : int or float
        The first value and the inclusive limit.
    step : int or float, optional
        The increment (negative counts down); default 1.

    Returns
    -------
    numpy.ndarray
        The sequence, int32 or float64.
    """
    if all(isinstance(a, (int, np.integer)) and not isinstance(a, bool) for a in (start, end, step)):
        member = _tools_member("Sequence", "Int32", "Int32", "Int32")
        return to_numpy_ints(member(int(start), int(end), int(step)))
    member = _tools_member("Sequence", "Double", "Double", "Double")
    return to_numpy_doubles(member(float(start), float(end), float(step)))


def compress(data: bytes) -> bytes:
    """Compress bytes with the .NET raw deflate stream at the optimal level (``Tools.Compress``).

    The output is a raw deflate stream, which ``zlib.decompress(out, -15)`` also reads.

    Parameters
    ----------
    data : bytes
        The bytes to compress.

    Returns
    -------
    bytes
        The compressed bytes.
    """
    source = _bytes_to_net(np.frombuffer(bytes(data), dtype=np.uint8), "bytes")
    return _bytes_from_net(namespace("Numerics").Tools.Compress(source), np.uint8).tobytes()


def decompress(data: bytes) -> bytes:
    """Inflate a raw deflate stream (``Tools.Decompress``).

    Parameters
    ----------
    data : bytes
        Bytes from :func:`compress` or any raw deflate stream.

    Returns
    -------
    bytes
        The decompressed bytes.
    """
    source = _bytes_to_net(np.frombuffer(bytes(data), dtype=np.uint8), "bytes")
    return _bytes_from_net(namespace("Numerics").Tools.Decompress(source), np.uint8).tobytes()


def _complex_to_net(values: Any) -> Any:
    """Copy a 1-D complex array into a ``System.Double[]`` of interleaved real and imaginary parts."""
    array = np.ascontiguousarray(values, dtype=np.complex128)
    if array.ndim != 1:
        raise ValueError(f"Expected a 1-D array; got shape {array.shape}.")
    return to_net_doubles(array.view(np.float64))


def _fourier() -> Any:
    """Return ``Numerics.Mathematics.Fourier`` after starting the runtime."""
    return namespace("Numerics.Mathematics").Fourier


def fft(values: Any, inverse: bool = False) -> np.ndarray:
    """Return the discrete Fourier transform of a power-of-2 complex sequence (``Fourier.FFT``).

    Numerics uses the Numerical Recipes conventions: the forward transform is
    H(k) = sum h(j) exp(+2 pi i j k / n), which equals ``n * numpy.fft.ifft(values)``, and the
    inverse transform returns n times the inverse DFT, which equals ``numpy.fft.fft(values)``.
    The array crosses once each way as interleaved doubles.

    Parameters
    ----------
    values : array_like of complex
        The sequence; its length must be a power of 2.
    inverse : bool, optional
        Compute the (unnormalized) inverse transform; default False.

    Returns
    -------
    numpy.ndarray
        The complex transform.

    Raises
    ------
    ValueError
        If the length is not a power of 2 (raised by Numerics).
    """
    data = _complex_to_net(values)
    _fourier().FFT(data, bool(inverse))
    return to_numpy_doubles(data).view(np.complex128)


def bluestein_fft(values: Any, inverse: bool = False) -> np.ndarray:
    """Return the discrete Fourier transform of a complex sequence of any length (``Fourier.BluesteinFFT``).

    The conventions are those of :func:`fft`; a length that is not a power of 2 is transformed
    with Bluestein's chirp-z convolution.

    Parameters
    ----------
    values : array_like of complex
        The sequence, of any positive length.
    inverse : bool, optional
        Compute the (unnormalized) inverse transform; default False.

    Returns
    -------
    numpy.ndarray
        The complex transform.

    Raises
    ------
    ValueError
        If the sequence is empty (raised by Numerics).
    """
    member = _fourier().BluesteinFFT
    data = _complex_to_net(values)
    member(data, bool(inverse))
    return to_numpy_doubles(data).view(np.complex128)


def real_fft(values: Any, inverse: bool = False) -> np.ndarray:
    """Return the packed transform of a power-of-2 real sequence (``Fourier.RealFFT``).

    The forward result holds the positive-frequency half of the :func:`fft` of the real input:
    element 0 is H(0), element 1 is H(n/2) (both real), and elements 2k, 2k + 1 are the real and
    imaginary parts of H(k) for k = 1 .. n/2 - 1. With ``inverse`` a packed array is turned back
    into the real sequence times n/2.

    Parameters
    ----------
    values : array_like
        The real sequence (or packed transform); its length must be a power of 2.
    inverse : bool, optional
        Compute the inverse transform; default False.

    Returns
    -------
    numpy.ndarray
        The packed transform, or the unscaled inverse.

    Raises
    ------
    ValueError
        If the length is not a power of 2 (raised by Numerics).
    """
    data = to_net_doubles(values)
    _fourier().RealFFT(data, bool(inverse))
    return to_numpy_doubles(data)


def fft_convolve(values1: Any, values2: Any) -> np.ndarray:
    """Return the full linear convolution of two real sequences by FFT (``Fourier.Convolve``).

    The same quantity as ``numpy.convolve(values1, values2, mode="full")``, with an absolute
    round-off of about machine epsilon times log2(n) times the largest output.

    Parameters
    ----------
    values1, values2 : array_like
        The sequences, each of length at least 1.

    Returns
    -------
    numpy.ndarray
        The convolution, of length ``len(values1) + len(values2) - 1``.

    Raises
    ------
    ValueError
        If either sequence is empty (raised by Numerics).
    """
    member = _fourier().Convolve
    return to_numpy_doubles(member(to_net_doubles(values1), to_net_doubles(values2)))


def fft_correlation(values1: Any, values2: Any) -> np.ndarray:
    """Return the circular cross-correlation of two real power-of-2 sequences (``Fourier.Correlation``).

    Element k is the sum over j of values1[(j + k) mod n] * values2[j], in wraparound order: lags
    0 .. n/2 - 1 first, then the negative lags -n/2 .. -1. This is not a normalized correlation
    coefficient.

    Parameters
    ----------
    values1, values2 : array_like
        Equal-length real sequences whose length is a power of 2.

    Returns
    -------
    numpy.ndarray
        The correlation at each lag.

    Raises
    ------
    ValueError
        If the lengths are not equal powers of 2 (raised by Numerics).
    """
    return to_numpy_doubles(_fourier().Correlation(to_net_doubles(values1), to_net_doubles(values2)))


def fft_autocorrelation(series: Any, lag_max: int) -> np.ndarray:
    """Return the sample autocorrelation function by FFT (``Fourier.Autocorrelation``).

    r(k) = sum (x_t - mean)(x_(t+k) - mean) / sum (x_t - mean)**2, the linear (not circular)
    estimate that ``statsmodels.tsa.stattools.acf`` also computes.

    Parameters
    ----------
    series : array_like
        The series.
    lag_max : int
        The largest lag.

    Returns
    -------
    numpy.ndarray
        A ``(lag_max + 1, 2)`` array: the lags 0 .. lag_max in column 0 and r(k) in column 1.
    """
    return to_numpy_matrix(_fourier().Autocorrelation(to_net_doubles(series), int(lag_max)))


__all__ = [
    "NumericsConstants", "arg_max", "arg_min", "bluestein_fft", "clamp", "compensated_cumulative_sum",
    "compensated_sum", "compress", "constants", "copy_sign",
    "decompress", "denormalize", "destandardize", "distance", "expm1", "fft", "fft_autocorrelation", "fft_convolve",
    "fft_correlation", "hypot", "is_finite", "is_power_of_two", "log", "log10", "log1p", "log_sum_exp", "max_value",
    "mean_values", "min_max", "min_value", "next_power_of_two", "normalize", "normalized_distance", "point_distance",
    "power", "product", "real_fft", "sequence", "sign", "sqr", "standardize", "sum_product", "sum_values",
    "weighted_distance",
]

translate_public(globals())
