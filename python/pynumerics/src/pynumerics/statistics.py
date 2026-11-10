"""Descriptive statistics from ``Numerics.Data.Statistics``.

Thin calls into ``Statistics`` (sequential sums, Type-7 and weighted percentiles, ranks and
L-moments), ``Correlation`` (Pearson, Spearman and Kendall's tau-b), ``BoxCox``,
``RunningStatistics`` (streaming moments), ``PlottingPositions`` and ``Autocorrelation``.
"""

from __future__ import annotations

import operator
from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import bind, bind_list, cached_overload, evaluate, loop, namespace, net_enum
from pynumerics.enums import AutocorrelationType, PlottingPosition, WeightType
from pynumerics.errors import translate_public
from pynumerics.marshal import _bound, to_net_doubles, to_net_matrix, to_numpy_doubles, to_numpy_ints, to_numpy_matrix


def _stats() -> Any:
    """Return the ``Numerics.Data.Statistics`` namespace after starting the runtime."""
    return namespace("Numerics.Data.Statistics")


def _statistics_member(name: str, *signature: str) -> Any:
    """Return one overload of a ``Statistics`` method group, bound once per process.

    ``signature`` names the parameter types: ``"list"`` (``IList<double>``), ``"double"``,
    ``"bool"``, ``"int"``, ``"kind"`` (``PlottingPositions.PlottingPostionType``) and ``"weight"``
    (``WeightType``). The bound overload is the member pythonnet would resolve for the same
    arguments; binding it once removes the per-call resolution (see ``_util.cached_overload``).
    """
    def factory() -> Any:
        stats = _stats()
        import System  # noqa: PLC0415 - the runtime is running
        from System.Collections.Generic import IList  # noqa: PLC0415

        types = {"list": IList[System.Double], "double": System.Double, "bool": System.Boolean, "int": System.Int32,
                 "kind": stats.PlottingPositions.PlottingPostionType, "weight": stats.WeightType}
        return getattr(stats.Statistics, name).Overloads[tuple(types[item] for item in signature)]

    return cached_overload(f"Statistics.{name}({', '.join(signature)})", factory)


def mean(data: Any, weights: Any = None) -> float:
    """Return the mean through ``Statistics.Mean`` (sequential summation).

    Parameters
    ----------
    data : array_like
        Values.
    weights : array_like, optional
        Weights for ``Statistics.Mean(IList, IList)``.

    Returns
    -------
    float
        The (weighted) mean.
    """
    if weights is None:
        return float(_statistics_member("Mean", "list")(to_net_doubles(data)))
    return float(_statistics_member("Mean", "list", "list")(to_net_doubles(data), to_net_doubles(weights)))


def variance(data: Any, weights: Any = None, weight_type: WeightType | str = WeightType.FREQUENCY) -> float:
    """Return the sample variance through ``Statistics.Variance``.

    Parameters
    ----------
    data : array_like
        Values.
    weights : array_like, optional
        Weights; when given ``weight_type`` selects frequency or reliability weights.
    weight_type : WeightType or str, optional
        Weight interpretation; default frequency.

    Returns
    -------
    float
        The variance.
    """
    stats = _stats()
    if weights is None:
        return float(_statistics_member("Variance", "list")(to_net_doubles(data)))
    member = _statistics_member("Variance", "list", "list", "weight")
    return float(member(to_net_doubles(data), to_net_doubles(weights), net_enum(stats.WeightType, weight_type)))


def standard_deviation(data: Any, weights: Any = None, weight_type: WeightType | str = WeightType.FREQUENCY) -> float:
    """Return the sample standard deviation through ``Statistics.StandardDeviation``.

    Parameters
    ----------
    data : array_like
        Values.
    weights : array_like, optional
        Weights; when given ``weight_type`` selects frequency or reliability weights.
    weight_type : WeightType or str, optional
        Weight interpretation; default frequency.

    Returns
    -------
    float
        The standard deviation.
    """
    stats = _stats()
    if weights is None:
        return float(_statistics_member("StandardDeviation", "list")(to_net_doubles(data)))
    member = _statistics_member("StandardDeviation", "list", "list", "weight")
    return float(member(to_net_doubles(data), to_net_doubles(weights), net_enum(stats.WeightType, weight_type)))


def median(data: Any, data_is_sorted: bool = False) -> float:
    """Return the median through ``Statistics.Median``.

    For an even count the result is the correctly rounded midpoint ``0.5 * (lo + hi)`` of the two
    middle values, as ``numpy.median`` computes it; ``percentile(data, 0.5)`` (Type 7) can differ
    from it in the last bit.

    Parameters
    ----------
    data : array_like
        Values.
    data_is_sorted : bool, optional
        Whether ``data`` is already ascending; default False.

    Returns
    -------
    float
        The median.
    """
    return float(_statistics_member("Median", "list", "bool")(to_net_doubles(data), bool(data_is_sorted)))


def percentile(data: Any, k: Any, weights: Any = None, data_is_sorted: bool = False) -> Any:
    """Return percentiles through ``Statistics.Percentile``.

    Unweighted percentiles are Type 7 (position (n − 1)k). Weighted percentiles use
    ``Statistics.Percentile(IList, k, IList weights, bool)``, which reduces to Type 7 at
    unit weights.

    Parameters
    ----------
    data : array_like
        Values.
    k : float or array_like
        Percentiles as fractions in [0, 1].
    weights : array_like, optional
        Weights, one per value.
    data_is_sorted : bool, optional
        Whether ``data`` is already ascending; default False.

    Returns
    -------
    float or numpy.ndarray
        Percentile values, shaped like ``k``.
    """
    values = to_net_doubles(data)
    karr = np.asarray(k, dtype=np.float64)
    scalar = karr.ndim == 0
    ks = float(karr) if scalar else to_net_doubles(karr.ravel())
    k_type = "double" if scalar else "list"
    if weights is None:
        result = _statistics_member("Percentile", "list", k_type, "bool")(values, ks, bool(data_is_sorted))
    else:
        result = _statistics_member("Percentile", "list", k_type, "list", "bool")(values, ks, to_net_doubles(weights),
                                                                                 bool(data_is_sorted))
    return float(result) if scalar else to_numpy_doubles(result).reshape(karr.shape)


def plotting_position_percentile(data: Any, k: Any,
                                 kind: PlottingPosition | str = PlottingPosition.WEIBULL,
                                 data_is_sorted: bool = False, weights: Any = None) -> Any:
    """Return percentiles interpolated between order statistics at plotting positions.

    Calls ``Statistics.Percentile(IList, k, PlottingPostionType, bool)``: the i-th order
    statistic sits at its plotting position and the percentile interpolates linearly between
    them, held at the sample extremes outside the first and last positions. Weibull is
    Hyndman-Fan type 6 (NumPy ``method="weibull"``), Hazen type 5 and Kaplan-Meier type 4. For
    Weibull, Blom and Hazen the interpolation weight is formed from ``k`` and the two bracketing
    positions, which is closer to exact rational evaluation than the fractional part of the rank.

    With ``weights`` it calls ``Statistics.Percentile(IList, k, IList weights, PlottingPostionType,
    bool)``, with an array ``k`` the list overload, which validates and sorts the sample once and
    serves every probability in one walk, each value bitwise the scalar call. Point i sits at the
    weighted position ``(c_(i-1) + (1 - alpha) w_i) / (W + (1 - 2 alpha) w_mean)``
    (``(c_(i-1) + w_i / 2) / W`` for Hazen, ``c_i / W`` for Kaplan-Meier), zero weights included.

    Parameters
    ----------
    data : array_like
        Values.
    k : float or array_like
        Percentiles as fractions in [0, 1].
    kind : PlottingPosition or str, optional
        The plotting-position formula; default Weibull.
    data_is_sorted : bool, optional
        Whether ``data`` is already ascending; default False.
    weights : array_like, optional
        The non-negative weight of each value, aligned with ``data``; default unweighted.

    Returns
    -------
    float or numpy.ndarray
        Percentile values, shaped like ``k``.
    """
    kind_net = net_enum(_stats().PlottingPositions.PlottingPostionType, kind)
    values = to_net_doubles(data)
    karr = np.asarray(k, dtype=np.float64)
    if weights is not None:
        k_type = "double" if karr.ndim == 0 else "list"
        member = _statistics_member("Percentile", "list", k_type, "list", "kind", "bool")
        net_weights = to_net_doubles(weights)
        if karr.ndim == 0:
            return float(member(values, float(karr), net_weights, kind_net, bool(data_is_sorted)))
        result = member(values, to_net_doubles(karr.ravel()), net_weights, kind_net, bool(data_is_sorted))
        return to_numpy_doubles(result).reshape(karr.shape)
    if karr.ndim == 0:
        return float(_statistics_member("Percentile", "list", "double", "kind", "bool")(
            values, float(karr), kind_net, bool(data_is_sorted)))
    result = _statistics_member("Percentile", "list", "list", "kind", "bool")(
        values, to_net_doubles(karr.ravel()), kind_net, bool(data_is_sorted))
    return to_numpy_doubles(result).reshape(karr.shape)


def column_percentiles(data: Any, k: Any,
                       kind: PlottingPosition | str | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return the percentiles of each column over its finite entries (``Statistics.ColumnPercentiles``).

    Each column's finite entries go through ``Statistics.Percentile`` in one call for the
    whole matrix, so every value is bitwise the per-column ``percentile`` (Type 7, the
    default) or ``plotting_position_percentile`` (when ``kind`` is given) of the finite
    entries; NaN and infinite entries are excluded rather than sorted in.

    Parameters
    ----------
    data : array_like
        Two-dimensional samples, one row per member (for example a realization) and one
        column per variable (for example a curve node).
    k : float or array_like
        Percentiles as fractions in [0, 1].
    kind : PlottingPosition or str, optional
        A plotting-position formula (``"weibull"`` is Hyndman-Fan type 6) for
        ``ColumnPercentiles(double[,], IList, PlottingPostionType, out int[])``; None (the
        default) gives the Type 7 percentiles.

    Returns
    -------
    percentiles : numpy.ndarray
        Shape ``(len(k), columns)`` (``(columns,)`` for a scalar ``k``); NaN for a column
        without a finite entry.
    excluded : numpy.ndarray
        The int32 number of non-finite entries excluded from each column.
    """
    stats = _stats()
    karr = np.asarray(k, dtype=np.float64)
    matrix, ks = to_net_matrix(data), to_net_doubles(karr.ravel())
    if kind is None:
        table, excluded = stats.Statistics.ColumnPercentiles(matrix, ks, None)
    else:
        kind_net = net_enum(stats.PlottingPositions.PlottingPostionType, kind)
        table, excluded = stats.Statistics.ColumnPercentiles(matrix, ks, kind_net, None)
    result = to_numpy_matrix(table)
    if karr.ndim == 0:
        result = result[0]
    return result, to_numpy_ints(excluded)


def column_means(data: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return the mean of each column over its finite entries (``Statistics.ColumnMeans``).

    Each column's finite entries go through ``Statistics.Mean`` (a sequential sum) in one call for
    the whole matrix, so every value is bitwise :func:`mean` of the column's finite entries; NaN and
    infinite entries are excluded, as in :func:`column_percentiles`.

    Parameters
    ----------
    data : array_like
        Two-dimensional samples, one row per member and one column per variable.

    Returns
    -------
    means : numpy.ndarray
        One mean per column; NaN for a column without a finite entry.
    excluded : numpy.ndarray
        The int32 number of non-finite entries excluded from each column.
    """
    means, excluded = _column_statistic("ColumnMeans", data)
    return to_numpy_doubles(means), excluded


def _column_statistic(member: str, data: Any) -> tuple[Any, np.ndarray]:
    """Call a ``Statistics`` column member through its ``(double[,], out int[])`` overload, bound once per process.

    The matrix crosses in one copy; the overload is resolved on first use only.
    """
    def factory() -> Any:
        group = getattr(_stats().Statistics, member)
        import clr  # noqa: PLC0415 - the runtime is running
        import System  # noqa: PLC0415

        return group.Overloads[_bound()["doubles_matrix"], clr.GetClrType(System.Array[System.Int32]).MakeByRefType()]

    values, excluded = cached_overload(f"Statistics.{member}(double[,], out int[])", factory)(to_net_matrix(data), None)
    return values, to_numpy_ints(excluded)


def column_variances(data: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return the sample variance of each column over its finite entries (``Statistics.ColumnVariances``).

    Every value is bitwise :func:`variance` (n - 1 divisor, corrected two-pass) of the column's
    finite entries, computed for the whole matrix in one call; NaN and infinite entries are excluded,
    as in :func:`column_means`. Equals ``numpy.nanvar(data, axis=0, ddof=1)`` to rounding when every
    entry is finite or NaN.

    Parameters
    ----------
    data : array_like
        Two-dimensional samples, one row per member and one column per variable.

    Returns
    -------
    variances : numpy.ndarray
        One variance per column; NaN for a column with fewer than two finite entries.
    excluded : numpy.ndarray
        The int32 number of non-finite entries excluded from each column.
    """
    values, excluded = _column_statistic("ColumnVariances", data)
    return to_numpy_doubles(values), excluded


def column_standard_deviations(data: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return the sample standard deviation of each column over its finite entries.

    Wraps ``Statistics.ColumnStandardDeviations``.

    Every value is bitwise :func:`standard_deviation` of the column's finite entries, the square
    root of the :func:`column_variances` value.

    Parameters
    ----------
    data : array_like
        Two-dimensional samples, one row per member and one column per variable.

    Returns
    -------
    standard_deviations : numpy.ndarray
        One standard deviation per column; NaN for a column with fewer than two finite entries.
    excluded : numpy.ndarray
        The int32 number of non-finite entries excluded from each column.
    """
    values, excluded = _column_statistic("ColumnStandardDeviations", data)
    return to_numpy_doubles(values), excluded


def column_linear_moments(data: Any) -> tuple[np.ndarray, np.ndarray]:
    """Return the L-moments of each column over its finite entries (``Statistics.ColumnLinearMoments``).

    Every column is bitwise :func:`linear_moments` of its finite entries. The L-CV of a column is
    ``moments[1] / moments[0]``, bitwise :func:`linear_coefficient_of_variation`.

    Parameters
    ----------
    data : array_like
        Two-dimensional samples, one row per member and one column per variable.

    Returns
    -------
    moments : numpy.ndarray
        A (4, p) array with rows lambda1, lambda2, tau3 and tau4; NaN for a column with fewer than
        four finite entries.
    excluded : numpy.ndarray
        The int32 number of non-finite entries excluded from each column.
    """
    values, excluded = _column_statistic("ColumnLinearMoments", data)
    return to_numpy_matrix(values), excluded


def wilson_score_interval(successes: float, trials: float, confidence_level: float) -> tuple[float, float]:
    """Return the Wilson score interval of a binomial proportion (``Statistics.WilsonScoreInterval``).

    Parameters
    ----------
    successes : float
        The number of successes, in [0, ``trials``]; may be fractional.
    trials : float
        The number of trials (positive).
    confidence_level : float
        The two-sided level in (0, 1).

    Returns
    -------
    tuple of float
        The lower and upper bounds, clamped to [0, 1].
    """
    bounds = _stats().Statistics.WilsonScoreInterval(float(successes), float(trials), float(confidence_level))
    return float(bounds.Item1), float(bounds.Item2)


def linear_moments(data: Any, nmom: int | None = None) -> np.ndarray:
    """Return the sample L-moments through ``Statistics.LinearMoments``.

    With ``nmom`` it calls ``Statistics.LinearMoments(IList, int)``, the unbiased estimator of
    Hosking's SAMLMR through probability weighted moments, whose first four values are bitwise those
    without ``nmom``; ``nmom=5`` adds the ratio ``t5`` used in regional frequency analysis.

    Parameters
    ----------
    data : array_like
        Values (sorted by Numerics).
    nmom : int, optional
        The number of L-moments, at least 1; default the four-value member.

    Returns
    -------
    numpy.ndarray
        ``[l1, l2, t3, t4]``: the first two L-moments and the L-skewness and L-kurtosis
        ratios, as Numerics returns them; with ``nmom``, ``[l1, l2, t3, ..., t_nmom]``. All NaN
        when there are fewer than four values (fewer than ``nmom`` with ``nmom``).

    Raises
    ------
    TypeError
        If ``nmom`` is not an integer, for example a float.
    ValueError
        If ``nmom`` is less than 1.
    """
    if nmom is None:
        return to_numpy_doubles(_statistics_member("LinearMoments", "list")(to_net_doubles(data)))
    member = _statistics_member("LinearMoments", "list", "int")
    return to_numpy_doubles(member(to_net_doubles(data), operator.index(nmom)))


def linear_coefficient_of_variation(data: Any) -> float:
    """Return the sample L-coefficient of variation through ``Statistics.LinearCoefficientOfVariation``.

    Parameters
    ----------
    data : array_like
        Values (sorted by Numerics).

    Returns
    -------
    float
        ``l2 / l1``, bitwise the ratio of the first two entries of :func:`linear_moments`; NaN for
        fewer than four values or any NaN.
    """
    return float(_stats().Statistics.LinearCoefficientOfVariation(to_net_doubles(data)))


def ranks(data: Any, return_ties: bool = False) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
    """Return the one-based average ranks of the values through ``Statistics.RanksInPlace``.

    Values that compare equal, infinities included, share their average rank, as
    ``scipy.stats.rankdata(method="average")`` assigns them. The input is not modified.

    Parameters
    ----------
    data : array_like
        Values, in any order.
    return_ties : bool, optional
        Also return the tie-run lengths of ``RanksInPlace(double[], out double[])``; default False.

    Returns
    -------
    ranks : numpy.ndarray
        One rank per value, in input order.
    ties : numpy.ndarray
        Only when ``return_ties`` is True: in sorted order, the last position of each run of k equal
        values holds k - 1 and every other entry is zero.
    """
    s = _stats().Statistics
    values = to_net_doubles(data)
    if not return_ties:
        return to_numpy_doubles(s.RanksInPlace(values))
    result, ties = s.RanksInPlace(values, None)
    return to_numpy_doubles(result), to_numpy_doubles(ties)


def _correlation(name: str, sample1: Any, sample2: Any) -> Any:
    """Call a ``Correlation`` member on two samples, or on the columns of one 2-D array."""
    correlation = _stats().Correlation
    if sample2 is None:
        return to_numpy_matrix(getattr(correlation, name)(to_net_matrix(sample1)))
    return float(getattr(correlation, name)(to_net_doubles(sample1), to_net_doubles(sample2)))


def pearson(sample1: Any, sample2: Any = None) -> Any:
    """Return the Pearson correlation through ``Correlation.Pearson``.

    Parameters
    ----------
    sample1 : array_like
        The first sample, or, when ``sample2`` is omitted, a 2-D array with one row per observation
        and one column per variable.
    sample2 : array_like, optional
        The second sample, of the same length.

    Returns
    -------
    float or numpy.ndarray
        The coefficient, or the correlation matrix of the columns of ``sample1``.
    """
    return _correlation("Pearson", sample1, sample2)


def spearman(sample1: Any, sample2: Any = None) -> Any:
    """Return the Spearman rank correlation through ``Correlation.Spearman``.

    The Pearson correlation of the average ranks of :func:`ranks`.

    Parameters
    ----------
    sample1 : array_like
        The first sample, or, when ``sample2`` is omitted, a 2-D array with one row per observation
        and one column per variable.
    sample2 : array_like, optional
        The second sample, of the same length.

    Returns
    -------
    float or numpy.ndarray
        The coefficient, or the rank correlation matrix of the columns of ``sample1``.
    """
    return _correlation("Spearman", sample1, sample2)


def kendalls_tau(sample1: Any, sample2: Any) -> float:
    """Return Kendall's tau-b through ``Correlation.KendallsTau``.

    ``(n_c - n_d) / sqrt(n_1 n_2)``, with ties counted as ``scipy.stats.kendalltau`` counts them.
    Pairs are counted in O(n log n) for samples of 48 or more values (Knight, 1966).

    Parameters
    ----------
    sample1, sample2 : array_like
        Samples of the same length.

    Returns
    -------
    float
        Tau-b; NaN when either sample contains NaN or is constant.
    """
    return float(_stats().Correlation.KendallsTau(to_net_doubles(sample1), to_net_doubles(sample2)))


def box_cox_transform(values: Any, lam: float) -> Any:
    """Return the Box-Cox transformation through ``BoxCox.Transform``.

    ``(x^lam - 1) / lam``, or ``ln x`` for ``|lam| < 1e-8``.

    Parameters
    ----------
    values : float or array_like
        Positive values.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float or numpy.ndarray
        The transformed values, shaped like ``values``; NaN for a non-positive value or ``|lam| > 5``.
    """
    box_cox = _stats().BoxCox
    exponent = float(lam)
    scalar = cached_overload("BoxCox.Transform(double, double)", lambda: bind(box_cox.Transform, 2))
    return evaluate(values, lambda x: scalar(x, exponent), lambda v: box_cox.Transform(v, exponent))


def box_cox_inverse_transform(values: Any, lam: float) -> Any:
    """Return the inverse Box-Cox transformation through ``BoxCox.InverseTransform``.

    Parameters
    ----------
    values : float or array_like
        Transformed values.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float or numpy.ndarray
        The original-scale values, shaped like ``values``.
    """
    box_cox = _stats().BoxCox
    exponent = float(lam)
    scalar = cached_overload("BoxCox.InverseTransform(double, double)", lambda: bind(box_cox.InverseTransform, 2))
    return evaluate(values, lambda y: scalar(y, exponent), lambda v: box_cox.InverseTransform(v, exponent))


def box_cox_derivative(values: Any, lam: float) -> Any:
    """Return the derivative of the Box-Cox transformation through ``BoxCox.Derivative``.

    ``x^(lam - 1)``, or ``1 / x`` where the transformation is the logarithm. Numerics has no array
    overload, so an array is evaluated element by element with the scalar member.

    Parameters
    ----------
    values : float or array_like
        Positive values.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float or numpy.ndarray
        The derivatives, shaped like ``values``; NaN for a non-positive value or ``|lam| > 5``.
    """
    member, exponent = _stats().BoxCox.Derivative, float(lam)
    return loop(values, lambda x: member(x, exponent))


def box_cox_fit_lambda(values: Any) -> tuple[bool, float]:
    """Fit the Box-Cox exponent by maximum likelihood through ``BoxCox.TryFitLambda``.

    Parameters
    ----------
    values : array_like
        Positive values, not all equal.

    Returns
    -------
    success : bool
        Whether a finite exponent was found.
    lam : float
        The exponent in [-5, 5], or NaN when the sample cannot support a fit or the search fails.
    """
    success, lam = _stats().BoxCox.TryFitLambda(to_net_doubles(values), 0.0)
    return bool(success), float(lam)


def plotting_positions(n: int, kind: PlottingPosition | str | float = PlottingPosition.WEIBULL) -> np.ndarray:
    """Return plotting-position probabilities through ``PlottingPositions.Function``.

    Parameters
    ----------
    n : int
        Sample size.
    kind : PlottingPosition, str or float, optional
        A named formula, or the plotting-position parameter α of (i − α)/(n + 1 − 2α).

    Returns
    -------
    numpy.ndarray
        The n ascending probabilities.
    """
    pp = _stats().PlottingPositions
    if isinstance(kind, (int, float)) and not isinstance(kind, bool):
        return to_numpy_doubles(pp.Function(int(n), float(kind)))
    return to_numpy_doubles(pp.Function(int(n), net_enum(pp.PlottingPostionType, kind)))


def weighted_plotting_positions(weights: Any, alpha: float | PlottingPosition | str) -> np.ndarray:
    """Return weighted plotting positions through ``PlottingPositions.Function(IList, ...)``.

    ``F_i = (c_i - alpha w_i) / (W + (1 - 2 alpha) W / n)`` with ``c_i`` the sequential cumulative
    weight and ``W`` the total; unit weights give bitwise the positions of
    :func:`plotting_positions`. A named formula selects ``Function(IList, PlottingPostionType)``:
    the alpha formulas at :func:`plotting_position_alpha`, and Kaplan-Meier as ``c_i / W``.

    Parameters
    ----------
    weights : array_like
        Non-negative weights of the observations in ascending order of their values.
    alpha : float, PlottingPosition or str
        The plotting-position parameter in [0, 0.5] (0 Weibull, 0.44 Gringorten, 0.5 Hazen; above
        one half a small last weight would give a position above one), or a named formula.

    Returns
    -------
    numpy.ndarray
        The non-decreasing positions, one per weight.
    """
    pp = _stats().PlottingPositions
    if isinstance(alpha, (int, float)) and not isinstance(alpha, bool):
        return to_numpy_doubles(pp.Function(to_net_doubles(weights), float(alpha)))
    return to_numpy_doubles(pp.Function(to_net_doubles(weights), net_enum(pp.PlottingPostionType, alpha)))


def plotting_position_alpha(kind: PlottingPosition | str) -> float:
    """Return the alpha coefficient of a named plotting-position formula (``PlottingPositions.Alpha``).

    Parameters
    ----------
    kind : PlottingPosition or str
        An alpha formula (Weibull, Median, Blom, Cunnane, Gringorten or Hazen).

    Returns
    -------
    float
        The coefficient of (i - alpha) / (n + 1 - 2 alpha).
    """
    pp = _stats().PlottingPositions
    return float(pp.Alpha(net_enum(pp.PlottingPostionType, kind)))


class RunningStatistics:
    """Streaming moments (``Numerics.Data.Statistics.RunningStatistics``).

    Parameters
    ----------
    values : array_like, optional
        Initial values to push.

    Attributes
    ----------
    net : Numerics.Data.Statistics.RunningStatistics
        The wrapped .NET object.
    """

    def __init__(self, values: Any = None) -> None:
        stats = _stats()
        self.net = stats.RunningStatistics() if values is None else stats.RunningStatistics(to_net_doubles(values))

    def push(self, values: Any) -> None:
        """Push one value or an array (``Push(double)`` / ``Push(IList<double>)``).

        Parameters
        ----------
        values : float or array_like
            Values to add.
        """
        array = np.asarray(values, dtype=np.float64)
        if array.ndim == 0:
            self.net.Push(float(array))
        else:
            self.net.Push(to_net_doubles(array.ravel()))

    def clear(self) -> None:
        """Reset to the empty state in place (``RunningStatistics.Clear``).

        Pushing the same values afterwards gives bitwise the statistics of a new accumulator.
        """
        self.net.Clear()

    @classmethod
    def combine(cls, a: RunningStatistics, b: RunningStatistics) -> RunningStatistics:
        """Merge two accumulators (``RunningStatistics.Combine``).

        Parameters
        ----------
        a, b : RunningStatistics
            Accumulators.

        Returns
        -------
        RunningStatistics
            The combined accumulator.
        """
        out = cls.__new__(cls)
        out.net = _stats().RunningStatistics.Combine(a.net, b.net)
        return out

    @property
    def count(self) -> int:
        """Number of values pushed."""
        return int(self.net.Count)

    @property
    def minimum(self) -> float:
        """Smallest value."""
        return float(self.net.Minimum)

    @property
    def maximum(self) -> float:
        """Largest value."""
        return float(self.net.Maximum)

    @property
    def mean(self) -> float:
        """Mean."""
        return float(self.net.Mean)

    @property
    def variance(self) -> float:
        """Sample variance."""
        return float(self.net.Variance)

    @property
    def standard_deviation(self) -> float:
        """Sample standard deviation."""
        return float(self.net.StandardDeviation)

    @property
    def skewness(self) -> float:
        """Sample skewness."""
        return float(self.net.Skewness)

    @property
    def kurtosis(self) -> float:
        """Sample kurtosis."""
        return float(self.net.Kurtosis)

    @property
    def population_variance(self) -> float:
        """Population variance (divisor n); 0 for a single value."""
        return float(self.net.PopulationVariance)

    @property
    def population_standard_deviation(self) -> float:
        """Population standard deviation (divisor n); 0 for a single value."""
        return float(self.net.PopulationStandardDeviation)

    @property
    def coefficient_of_variation(self) -> float:
        """Coefficient of variation, the sample standard deviation over the mean."""
        return float(self.net.CoefficientOfVariation)

    @property
    def population_skewness(self) -> float:
        """Population skewness (type 1, no bias correction)."""
        return float(self.net.PopulationSkewness)

    @property
    def population_kurtosis(self) -> float:
        """Population kurtosis (type 1, no bias correction)."""
        return float(self.net.PopulationKurtosis)

    def clone(self) -> RunningStatistics:
        """Return an independent copy of the accumulator (``RunningStatistics.Clone``).

        Returns
        -------
        RunningStatistics
            A new accumulator with the same state.
        """
        out = type(self).__new__(type(self))
        out.net = self.net.Clone()
        return out


def autocorrelation(data: Any, max_lag: int | None = None,
                    kind: AutocorrelationType | str = AutocorrelationType.CORRELATION,
                    centered: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Return a sample autocorrelation function through ``Autocorrelation.Function``.

    ``COVARIANCE`` is ``R(k) = n^-1 sum_{t=k}^{n-1} (x_t - xbar)(x_{t-k} - xbar)`` (divisor ``n`` at
    every lag, the mean by ``Statistics.Mean``), ``CORRELATION`` is ``R(k) / R(0)`` and ``PARTIAL``
    is the Durbin-Levinson partial autocorrelation from lag 1. With ``centered=False`` the mean is
    taken as zero, for a process whose mean is known to be zero: ``R(k) = n^-1 sum x_t x_{t-k}`` and
    the correlation is ``sum x_t x_{t-k} / sum x_t^2``.

    Parameters
    ----------
    data : array_like
        The series, in time order: one-dimensional, finite, with at least two values.
    max_lag : int, optional
        The largest lag, an integer in [1, n - 1]; default Numerics' ``floor(min(10 log10(n), n - 1))``.
    kind : AutocorrelationType or str, optional
        ``CORRELATION`` (default), ``COVARIANCE`` or ``PARTIAL``.
    centered : bool, optional
        Remove the sample mean (default True); False evaluates about zero.

    Returns
    -------
    lags : numpy.ndarray
        The lags (``0..max_lag``; ``1..max_lag`` for ``PARTIAL``), as integers.
    values : numpy.ndarray
        The function at those lags.

    Raises
    ------
    ValueError
        If ``data`` is not one-dimensional, holds a non-finite value or has fewer than two values;
        if ``max_lag`` is not an integer in [1, n - 1]; or if Numerics returns no function (a
        zero-variance series for ``CORRELATION``). Each message names the argument at fault.
    """
    series = np.asarray(data, dtype=np.float64)
    if series.ndim != 1:
        raise ValueError(f"data must be one-dimensional; got shape {series.shape}.")
    if series.size < 2:
        raise ValueError(f"data needs at least two values; got {series.size}.")
    if not np.all(np.isfinite(series)):
        raise ValueError("data must be finite; it holds NaN or infinite values.")
    if max_lag is None:
        lag = -1
    else:
        if isinstance(max_lag, (bool, np.bool_)) or not isinstance(max_lag, (int, np.integer)):
            raise ValueError(f"max_lag must be an integer; got {max_lag!r}.")
        lag = int(max_lag)
        if not 1 <= lag <= series.size - 1:
            raise ValueError(f"max_lag must lie in [1, {series.size - 1}] for {series.size} values; got {lag}.")
    stats = _stats()
    net_kind = net_enum(stats.Autocorrelation.Type, kind)
    if centered:
        table = stats.Autocorrelation.Function(to_net_doubles(series), lag, net_kind)
    else:
        table = stats.Autocorrelation.Function(to_net_doubles(series), lag, net_kind, False)
    if table is None:
        raise ValueError("Autocorrelation.Function returned no function: the series needs at least two values, "
                         "a maximum lag of at least 1 and, for the correlation, a nonzero variance.")
    result = to_numpy_matrix(table)
    return result[:, 0].astype(np.int64), np.ascontiguousarray(result[:, 1])


def _one(name: str, data: Any) -> float:
    """Call a one-sample ``Statistics`` member through its bound ``IList<double>`` overload."""
    return float(_statistics_member(name, "list")(to_net_doubles(data)))


def minimum(data: Any) -> float:
    """Return the smallest value through ``Statistics.Minimum``.

    Parameters
    ----------
    data : array_like
        Values, in any order.

    Returns
    -------
    float
        The minimum; NaN for an empty sample or any NaN.
    """
    return _one("Minimum", data)


def maximum(data: Any) -> float:
    """Return the largest value through ``Statistics.Maximum``.

    Parameters
    ----------
    data : array_like
        Values, in any order.

    Returns
    -------
    float
        The maximum; NaN for an empty sample or any NaN.
    """
    return _one("Maximum", data)


def sum(data: Any) -> float:
    """Return the sequential sum through ``Statistics.Sum``.

    Parameters
    ----------
    data : array_like
        Values, in any order.

    Returns
    -------
    float
        The sum in index order; NaN for an empty sample or any NaN.
    """
    return _one("Sum", data)


def parallel_mean(data: Any) -> float:
    """Return the mean through ``Statistics.ParallelMean`` (deterministic chunked reduction).

    The chunk boundaries depend only on the sample length, so the result is reproducible across
    runs and machines; it can differ from :func:`mean` in the last bits.

    Parameters
    ----------
    data : array_like
        Values, in any order.

    Returns
    -------
    float
        The mean; NaN for an empty sample or any NaN.
    """
    return _one("ParallelMean", data)


def geometric_mean(data: Any) -> float:
    """Return the geometric mean through ``Statistics.GeometricMean``.

    Parameters
    ----------
    data : array_like
        Positive values.

    Returns
    -------
    float
        The geometric mean; NaN for an empty sample or any NaN.
    """
    return _one("GeometricMean", data)


def harmonic_mean(data: Any) -> float:
    """Return the harmonic mean through ``Statistics.HarmonicMean``.

    Parameters
    ----------
    data : array_like
        Nonzero values.

    Returns
    -------
    float
        The harmonic mean; NaN for an empty sample or any NaN.
    """
    return _one("HarmonicMean", data)


def population_variance(data: Any) -> float:
    """Return the population variance (divisor n) through ``Statistics.PopulationVariance``.

    Parameters
    ----------
    data : array_like
        Values.

    Returns
    -------
    float
        The variance with divisor n; NaN for an empty sample or any NaN.
    """
    return _one("PopulationVariance", data)


def population_standard_deviation(data: Any) -> float:
    """Return the population standard deviation (divisor n) through ``Statistics.PopulationStandardDeviation``.

    Parameters
    ----------
    data : array_like
        Values.

    Returns
    -------
    float
        The standard deviation with divisor n; NaN for an empty sample or any NaN.
    """
    return _one("PopulationStandardDeviation", data)


def _pair(result: Any) -> tuple[float, float]:
    """Return a .NET ``ValueTuple<double, double>`` as a Python pair."""
    return float(result.Item1), float(result.Item2)


def mean_variance(data: Any) -> tuple[float, float]:
    """Return the mean and sample variance in one pass through ``Statistics.MeanVariance``.

    Parameters
    ----------
    data : array_like
        Values.

    Returns
    -------
    tuple of float
        The mean and the sample variance (divisor n - 1).
    """
    return _pair(_stats().Statistics.MeanVariance(to_net_doubles(data)))


def mean_standard_deviation(data: Any) -> tuple[float, float]:
    """Return the mean and sample standard deviation through ``Statistics.MeanStandardDeviation``.

    Parameters
    ----------
    data : array_like
        Values.

    Returns
    -------
    tuple of float
        The mean and the sample standard deviation (divisor n - 1).
    """
    return _pair(_stats().Statistics.MeanStandardDeviation(to_net_doubles(data)))


def coefficient_of_variation(data: Any) -> float:
    """Return the coefficient of variation through ``Statistics.CoefficientOfVariation``.

    Parameters
    ----------
    data : array_like
        Values.

    Returns
    -------
    float
        The sample standard deviation over the mean.
    """
    return _one("CoefficientOfVariation", data)


def _weighted_moment(name: str, data: Any, weights: Any, weight_type: WeightType | str) -> float:
    """Call ``Statistics.Skewness`` or ``Kurtosis``, unweighted or with weights."""
    if weights is None:
        return _one(name, data)
    member = _statistics_member(name, "list", "list", "weight")
    return float(member(to_net_doubles(data), to_net_doubles(weights), net_enum(_stats().WeightType, weight_type)))


def skewness(data: Any, weights: Any = None, weight_type: WeightType | str = WeightType.FREQUENCY) -> float:
    """Return the bias-corrected sample skewness through ``Statistics.Skewness``.

    The adjusted Fisher-Pearson coefficient ``G1``, as ``scipy.stats.skew(bias=False)``. With
    weights the correction uses the effective sample size of ``weight_type``.

    Parameters
    ----------
    data : array_like
        Values.
    weights : array_like, optional
        Non-negative, finite weights, one per value.
    weight_type : WeightType or str, optional
        Weight interpretation; default frequency.

    Returns
    -------
    float
        The skewness; NaN for an empty sample, any NaN, or an effective sample size of two or less.
    """
    return _weighted_moment("Skewness", data, weights, weight_type)


def kurtosis(data: Any, weights: Any = None, weight_type: WeightType | str = WeightType.FREQUENCY) -> float:
    """Return the bias-corrected sample excess kurtosis through ``Statistics.Kurtosis``.

    The small-sample estimator ``G2``, as ``scipy.stats.kurtosis(bias=False)``. With weights the
    correction uses the effective sample size of ``weight_type``.

    Parameters
    ----------
    data : array_like
        Values.
    weights : array_like, optional
        Non-negative, finite weights, one per value.
    weight_type : WeightType or str, optional
        Weight interpretation; default frequency.

    Returns
    -------
    float
        The excess kurtosis; NaN for an empty sample, any NaN, or an effective sample size of three
        or less.
    """
    return _weighted_moment("Kurtosis", data, weights, weight_type)


def covariance(data1: Any, data2: Any) -> float:
    """Return the sample covariance (divisor n - 1) through ``Statistics.Covariance``.

    Parameters
    ----------
    data1, data2 : array_like
        Samples of the same length.

    Returns
    -------
    float
        The covariance; NaN for fewer than two pairs or any NaN.
    """
    return float(_stats().Statistics.Covariance(to_net_doubles(data1), to_net_doubles(data2)))


def population_covariance(data1: Any, data2: Any) -> float:
    """Return the population covariance (divisor n) through ``Statistics.PopulationCovariance``.

    Parameters
    ----------
    data1, data2 : array_like
        Samples of the same length.

    Returns
    -------
    float
        The covariance; NaN for an empty sample or any NaN.
    """
    return float(_stats().Statistics.PopulationCovariance(to_net_doubles(data1), to_net_doubles(data2)))


def product_moments(data: Any) -> np.ndarray:
    """Return the first four product moments through ``Statistics.ProductMoments``.

    Parameters
    ----------
    data : array_like
        Values.

    Returns
    -------
    numpy.ndarray
        ``[mean, standard deviation, skewness, excess kurtosis]`` with the bias-corrected skewness
        and kurtosis; four NaN values for fewer than four values or any NaN.
    """
    return to_numpy_doubles(_stats().Statistics.ProductMoments(to_net_doubles(data)))


def five_number_summary(data: Any) -> np.ndarray:
    """Return the five-number summary through ``Statistics.FiveNumberSummary``.

    Parameters
    ----------
    data : array_like
        Values, in any order.

    Returns
    -------
    numpy.ndarray
        ``[min, 25th, 50th, 75th percentile, max]``, the percentiles Type 7.
    """
    return to_numpy_doubles(_stats().Statistics.FiveNumberSummary(to_net_doubles(data)))


def seven_number_summary(data: Any) -> np.ndarray:
    """Return the seven-number summary through ``Statistics.SevenNumberSummary``.

    Parameters
    ----------
    data : array_like
        Values, in any order.

    Returns
    -------
    numpy.ndarray
        ``[min, 5th, 25th, 50th, 75th, 95th percentile, max]``, the percentiles Type 7.
    """
    return to_numpy_doubles(_stats().Statistics.SevenNumberSummary(to_net_doubles(data)))


def jackknife_standard_error(data: Any, statistic: Any) -> float:
    """Return the jackknife standard error of a statistic through ``Statistics.JackKnifeStandardError``.

    Numerics evaluates the leave-one-out samples in parallel. A Python callable serializes on the
    GIL; pass C# source (compiled with :func:`pynumerics.callbacks.compile_csharp`) for a
    statistic that runs without crossing into Python.

    Parameters
    ----------
    data : array_like
        Values.
    statistic : callable, str or System.Delegate
        The statistic of one sample, a ``Func<IList<double>, double>``: a Python callable taking a
        NumPy array, C# source, or a .NET delegate.

    Returns
    -------
    float
        The standard error; 0 for a single value.
    """
    fn = callbacks.to_delegate(statistic, "Func<IList<double>,double>")
    return float(_stats().Statistics.JackKnifeStandardError(to_net_doubles(data), fn))


def jackknife_sample(data: Any, statistic: Any) -> np.ndarray | None:
    """Return the statistic of every leave-one-out sample through ``Statistics.JackKnifeSample``.

    Parameters
    ----------
    data : array_like
        Values.
    statistic : callable, str or System.Delegate
        The statistic of one sample, as for :func:`jackknife_standard_error`.

    Returns
    -------
    numpy.ndarray or None
        One value per left-out index, or None for an empty sample.
    """
    fn = callbacks.to_delegate(statistic, "Func<IList<double>,double>")
    result = _stats().Statistics.JackKnifeSample(to_net_doubles(data), fn)
    return None if result is None else to_numpy_doubles(result)


def entropy(data: Any, pdf: Any) -> float:
    """Return ``-sum p log p`` over the density values at the data through ``Statistics.Entropy``.

    Values where the density is zero or negative are skipped.

    Parameters
    ----------
    data : array_like
        Points at which the density is evaluated.
    pdf : callable, str or System.Delegate
        The density, a ``Func<double, double>``: a Python callable (for example a distribution's
        ``pdf``), C# source, or a .NET delegate.

    Returns
    -------
    float
        The entropy sum.
    """
    return float(_stats().Statistics.Entropy(to_net_doubles(data), callbacks.scalar_function(pdf)))


def effective_sample_size(weights: Any) -> float:
    """Return the Kish effective sample size ``(sum w)^2 / sum w^2`` through ``Statistics.EffectiveSampleSize``.

    The total and the sum of squares are accumulated sequentially in index order, the order the
    reliability-weighted statistics use. This is the effective size of a weighted sample, not of
    a correlated series: for the effective size of MCMC draws use
    :func:`pynumerics.mcmc.effective_sample_size`.

    Parameters
    ----------
    weights : array_like
        Non-negative, finite weights.

    Returns
    -------
    float
        The effective sample size, between one and the number of positive weights; NaN for an
        empty input.

    Raises
    ------
    ValueError
        If a weight is negative or not finite.
    """
    return float(_stats().Statistics.EffectiveSampleSize(to_net_doubles(weights)))


def box_cox_log_likelihood(values: Any, lam: float) -> float:
    """Return the Box-Cox profile log-likelihood through ``BoxCox.LogLikelihood``.

    Parameters
    ----------
    values : array_like
        Positive values.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float
        The normal log-likelihood of the transformed values plus the log-Jacobian.
    """
    return float(_stats().BoxCox.LogLikelihood(to_net_doubles(values), float(lam)))


def box_cox_log_jacobian(values: Any, lam: float) -> float:
    """Return the log-Jacobian ``(lam - 1) sum ln x`` of the Box-Cox transformation (``BoxCox.LogJacobian``).

    Parameters
    ----------
    values : array_like
        Positive values.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float
        The log-Jacobian.
    """
    return float(_stats().BoxCox.LogJacobian(to_net_doubles(values), float(lam)))


def yeo_johnson_transform(values: Any, lam: float) -> Any:
    """Return the Yeo-Johnson transformation through ``YeoJohnson.Transform``.

    Parameters
    ----------
    values : float or array_like
        Values of any sign.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float or numpy.ndarray
        The transformed values, shaped like ``values``.
    """
    yeo_johnson = _stats().YeoJohnson
    exponent = float(lam)
    scalar = cached_overload("YeoJohnson.Transform(double, double)", lambda: bind(yeo_johnson.Transform, 2))
    return evaluate(values, lambda x: scalar(x, exponent), lambda v: yeo_johnson.Transform(v, exponent))


def yeo_johnson_inverse_transform(values: Any, lam: float) -> Any:
    """Return the inverse Yeo-Johnson transformation through ``YeoJohnson.InverseTransform``.

    Parameters
    ----------
    values : float or array_like
        Transformed values.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float or numpy.ndarray
        The original-scale values, shaped like ``values``.
    """
    yeo_johnson = _stats().YeoJohnson
    exponent = float(lam)
    scalar = cached_overload("YeoJohnson.InverseTransform(double, double)",
                             lambda: bind(yeo_johnson.InverseTransform, 2))
    return evaluate(values, lambda y: scalar(y, exponent), lambda v: yeo_johnson.InverseTransform(v, exponent))


def yeo_johnson_derivative(values: Any, lam: float) -> Any:
    """Return the derivative of the Yeo-Johnson transformation through ``YeoJohnson.Derivative``.

    Numerics has no array overload, so an array is evaluated element by element.

    Parameters
    ----------
    values : float or array_like
        Values of any sign.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float or numpy.ndarray
        The derivatives, shaped like ``values``.
    """
    member, exponent = _stats().YeoJohnson.Derivative, float(lam)
    return loop(values, lambda x: member(x, exponent))


def yeo_johnson_fit_lambda(values: Any) -> float:
    """Fit the Yeo-Johnson exponent by maximum likelihood through ``YeoJohnson.FitLambda``.

    Parameters
    ----------
    values : array_like
        At least two values.

    Returns
    -------
    float
        The exponent in [-5, 5], or NaN when the search fails.

    Raises
    ------
    ValueError
        If fewer than two values are given.
    """
    member = cached_overload("YeoJohnson.FitLambda(IList)", lambda: bind_list(_stats().YeoJohnson.FitLambda))
    return float(member(to_net_doubles(values)))


def yeo_johnson_log_likelihood(values: Any, lam: float) -> float:
    """Return the Yeo-Johnson profile log-likelihood through ``YeoJohnson.LogLikelihood``.

    Parameters
    ----------
    values : array_like
        Values of any sign.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float
        The normal log-likelihood of the transformed values plus the log-Jacobian.
    """
    return float(_stats().YeoJohnson.LogLikelihood(to_net_doubles(values), float(lam)))


def yeo_johnson_log_jacobian(values: Any, lam: float) -> float:
    """Return the log-Jacobian of the Yeo-Johnson transformation through ``YeoJohnson.LogJacobian``.

    Parameters
    ----------
    values : array_like
        Values of any sign.
    lam : float
        The exponent, in [-5, 5].

    Returns
    -------
    float
        ``(lam - 1) sum sign(x) ln(1 + |x|)``.
    """
    return float(_stats().YeoJohnson.LogJacobian(to_net_doubles(values), float(lam)))


def autocorrelation_confidence_interval(sample_size: int, interval: float = 0.95) -> tuple[float, float]:
    """Return the white-noise band of the sample autocorrelation (``Autocorrelation.CorrelationConfidenceInterval``).

    Parameters
    ----------
    sample_size : int
        The series length.
    interval : float, optional
        The two-sided level; default 0.95.

    Returns
    -------
    tuple of float
        The lower and upper bounds ``-+ z / sqrt(n)``.
    """
    bounds = _stats().Autocorrelation.CorrelationConfidenceInterval(int(sample_size), float(interval))
    return float(bounds[0]), float(bounds[1])


def _sorted(values: Any, data_is_sorted: bool) -> Any:
    """Return the values as a ``System.Double[]``, sorted ascending unless already sorted."""
    array = np.asarray(values, dtype=np.float64)
    return to_net_doubles(array if data_is_sorted else np.sort(array))


def first_order_sobol(x: Any, y: Any, bins: int = 20) -> float:
    """Estimate the first-order Sobol index of an input from given data (``GlobalSensitivity.FirstOrderSobol``).

    The variance of the conditional output mean over equal-frequency input bins, divided by the
    output variance (Plischke, Borgonovo and Smith, 2013).

    Parameters
    ----------
    x : array_like
        The input sample.
    y : array_like
        The output sample, aligned with ``x``.
    bins : int, optional
        The number of equal-frequency input bins; default 20.

    Returns
    -------
    float
        The index, clamped to [0, 1]; 0 for a constant output.

    Raises
    ------
    ValueError
        If the samples differ in length, are shorter than ``bins``, ``bins`` is below two, or a
        value is not finite.
    """
    return float(_stats().GlobalSensitivity.FirstOrderSobol(to_net_doubles(x), to_net_doubles(y), int(bins)))


def pawn(x: Any, y: Any, bins: int = 20) -> np.ndarray:
    """Return the per-bin PAWN statistics of an input from given data (``GlobalSensitivity.Pawn``).

    For each equal-frequency input bin, the two-sample Kolmogorov-Smirnov statistic between the
    conditional and the unconditional output (Pianosi and Wagener, 2015).

    Parameters
    ----------
    x : array_like
        The input sample.
    y : array_like
        The output sample, aligned with ``x``.
    bins : int, optional
        The number of equal-frequency input bins; default 20.

    Returns
    -------
    numpy.ndarray
        One statistic per bin, in bin order.

    Raises
    ------
    ValueError
        As for :func:`first_order_sobol`.
    """
    return to_numpy_doubles(_stats().GlobalSensitivity.Pawn(to_net_doubles(x), to_net_doubles(y), int(bins)))


def pawn_median(x: Any, y: Any, bins: int = 20) -> float:
    """Return the PAWN index, the median of :func:`pawn` (``GlobalSensitivity.PawnMedian``).

    Parameters
    ----------
    x : array_like
        The input sample.
    y : array_like
        The output sample, aligned with ``x``.
    bins : int, optional
        The number of equal-frequency input bins; default 20.

    Returns
    -------
    float
        The median per-bin statistic.

    Raises
    ------
    ValueError
        As for :func:`first_order_sobol`.
    """
    return float(_stats().GlobalSensitivity.PawnMedian(to_net_doubles(x), to_net_doubles(y), int(bins)))


def borgonovo_delta(x: Any, y: Any, x_bins: int = 20, y_bins: int = 20) -> float:
    """Estimate the Borgonovo moment-independent delta from given data (``GlobalSensitivity.BorgonovoDelta``).

    Half the bin-weighted total variation between conditional and unconditional frequencies of
    equal-frequency output classes; invariant to increasing transforms of the output.

    Parameters
    ----------
    x : array_like
        The input sample.
    y : array_like
        The output sample, aligned with ``x``.
    x_bins : int, optional
        The number of equal-frequency input bins; default 20.
    y_bins : int, optional
        The number of equal-frequency output classes; default 20.

    Returns
    -------
    float
        The delta, in [0, 1 - 1 / y_bins].

    Raises
    ------
    ValueError
        As for :func:`first_order_sobol`, or if ``y_bins`` is below two or above the sample size.
    """
    return float(_stats().GlobalSensitivity.BorgonovoDelta(to_net_doubles(x), to_net_doubles(y), int(x_bins),
                                                           int(y_bins)))


class HistogramBin:
    """One histogram bin (``Numerics.Data.Statistics.Histogram.Bin``).

    Parameters
    ----------
    lower_bound : float
        The lower bound.
    upper_bound : float
        The upper bound.
    frequency : int, optional
        The count; default 0.

    Attributes
    ----------
    net : Numerics.Data.Statistics.Histogram.Bin
        The wrapped .NET object.
    """

    __hash__ = None  # mutable: compared by value, never hashed

    def __init__(self, lower_bound: float, upper_bound: float, frequency: int = 0) -> None:
        self.net = _stats().Histogram.Bin(float(lower_bound), float(upper_bound), int(frequency))

    @classmethod
    def _from_net(cls, net: Any) -> HistogramBin:
        out = cls.__new__(cls)
        out.net = net
        return out

    @property
    def lower_bound(self) -> float:
        """Lower bound."""
        return float(self.net.LowerBound)

    @lower_bound.setter
    def lower_bound(self, value: float) -> None:
        self.net.LowerBound = float(value)

    @property
    def upper_bound(self) -> float:
        """Upper bound."""
        return float(self.net.UpperBound)

    @upper_bound.setter
    def upper_bound(self, value: float) -> None:
        self.net.UpperBound = float(value)

    @property
    def midpoint(self) -> float:
        """Midpoint of the bounds."""
        return float(self.net.Midpoint)

    @property
    def frequency(self) -> int:
        """Count of values in the bin."""
        return int(self.net.Frequency)

    @frequency.setter
    def frequency(self, value: int) -> None:
        self.net.Frequency = int(value)

    def clone(self) -> HistogramBin:
        """Return an independent copy (``Bin.Clone``).

        Returns
        -------
        HistogramBin
            The copy.
        """
        return HistogramBin._from_net(self.net.Clone())

    def compare_to(self, other: HistogramBin | None) -> int:
        """Order two non-overlapping bins (``Bin.CompareTo``).

        Parameters
        ----------
        other : HistogramBin or None
            The bin to compare with.

        Returns
        -------
        int
            0 when both bounds are bitwise equal, +1 when this bin is higher (or ``other`` is
            None), -1 otherwise.
        """
        return int(self.net.CompareTo(None if other is None else other.net))

    def __eq__(self, other: object) -> bool:
        return isinstance(other, HistogramBin) and bool(self.net.Equals(other.net))

    def __repr__(self) -> str:
        return f"HistogramBin({self.lower_bound!r}, {self.upper_bound!r}, {self.frequency!r})"


class Histogram:
    """An adaptive histogram (``Numerics.Data.Statistics.Histogram``).

    Equal-width bins span the data range; values added later outside the range widen the first or
    last bin. A value on an interior edge counts in the bin above it.

    Parameters
    ----------
    data : array_like
        The initial values.
    number_of_bins : int, optional
        The bin count; default Numerics' ``ceil(2 n^(1/3)) + 1``.

    Attributes
    ----------
    net : Numerics.Data.Statistics.Histogram
        The wrapped .NET object.
    """

    def __init__(self, data: Any, number_of_bins: int | None = None) -> None:
        values = to_net_doubles(data)
        h = _stats().Histogram
        self.net = h(values) if number_of_bins is None else h(values, int(number_of_bins))

    @property
    def lower_bound(self) -> float:
        """Lower edge of the first bin."""
        return float(self.net.LowerBound)

    @property
    def upper_bound(self) -> float:
        """Upper edge of the last bin."""
        return float(self.net.UpperBound)

    @property
    def number_of_bins(self) -> int:
        """Number of bins."""
        return int(self.net.NumberOfBins)

    @property
    def bin_width(self) -> float:
        """Initial bin width."""
        return float(self.net.BinWidth)

    @property
    def data_count(self) -> int:
        """Total count over the bins."""
        return int(self.net.DataCount)

    @property
    def mean(self) -> float:
        """Frequency-weighted mean of the bin midpoints; 0 when empty."""
        return float(self.net.Mean)

    @property
    def median(self) -> float:
        """Midpoint of the bin holding the median count; NaN when empty."""
        return float(self.net.Median)

    @property
    def mode(self) -> float:
        """Midpoint of the most frequent bin."""
        return float(self.net.Mode)

    @property
    def standard_deviation(self) -> float:
        """Frequency-weighted standard deviation of the bin midpoints (divisor n); 0 when empty."""
        return float(self.net.StandardDeviation)

    def __len__(self) -> int:
        return int(self.net.NumberOfBins)

    def __getitem__(self, index: int) -> HistogramBin:
        """Return a copy of the bin at an index, in ascending order (``Histogram[index]``).

        Parameters
        ----------
        index : int
            The bin index; negative indices count from the end.

        Returns
        -------
        HistogramBin
            A copy of the bin.

        Raises
        ------
        IndexError
            If the index is out of range.
        """
        count = len(self)
        position = int(index) + count if int(index) < 0 else int(index)
        if not 0 <= position < count:
            raise IndexError(f"Bin index {index} is out of range for {count} bins.")
        return HistogramBin._from_net(self.net[position])

    def bins(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the lower bounds, upper bounds and frequencies of every bin, in ascending order.

        Returns
        -------
        lower : numpy.ndarray
            The lower bounds.
        upper : numpy.ndarray
            The upper bounds.
        frequency : numpy.ndarray
            The int64 counts.
        """
        items = [self.net[i] for i in range(len(self))]
        return (np.array([float(b.LowerBound) for b in items]), np.array([float(b.UpperBound) for b in items]),
                np.array([int(b.Frequency) for b in items], dtype=np.int64))

    def add_bin(self, bin: HistogramBin) -> None:
        """Append a bin (``Histogram.AddBin``).

        Parameters
        ----------
        bin : HistogramBin
            The bin; it must not overlap an existing one.
        """
        self.net.AddBin(bin.net)

    def add_data(self, values: Any) -> None:
        """Count one value or an array in one call (``Histogram.AddData``).

        Parameters
        ----------
        values : float or array_like
            The values; one outside the range widens the first or last bin.
        """
        array = np.asarray(values, dtype=np.float64)
        if array.ndim == 0:
            self.net.AddData(float(array))
        else:
            self.net.AddData(to_net_doubles(array.ravel()))

    def bin_index(self, value: float) -> int:
        """Return the index of the bin containing a value (``Histogram.GetBinIndexOf``).

        Parameters
        ----------
        value : float
            A value inside the histogram limits.

        Returns
        -------
        int
            The bin index.

        Raises
        ------
        ValueError
            If the value lies outside the histogram limits.
        """
        return int(self.net.GetBinIndexOf(float(value)))


class RunningCovarianceMatrix:
    """A streaming, identity-regularized covariance matrix (``Numerics.Data.Statistics.RunningCovarianceMatrix``).

    The co-moment matrix starts at the identity, so after N vectors ``covariance`` is ``M2 + I``,
    ``sample_covariance`` is ``(M2 + I) / (N - 1)`` and ``population_covariance`` is
    ``(M2 + I) / N``: the sample statistics shrunk toward the identity with a weight that fades as
    1 / N. Subtract the identity from ``covariance`` for the plain co-moment matrix.

    Parameters
    ----------
    size : int
        The vector length.

    Attributes
    ----------
    net : Numerics.Data.Statistics.RunningCovarianceMatrix
        The wrapped .NET object.
    """

    def __init__(self, size: int) -> None:
        self.net = _stats().RunningCovarianceMatrix(int(size))

    def push(self, values: Any) -> None:
        """Add one vector, or each row of a 2-D array, with a Welford update (``Push``).

        Numerics takes one vector per call, so a 2-D array crosses once per row.

        Parameters
        ----------
        values : array_like
            A vector of length ``size``, or an array of shape ``(m, size)``.

        Raises
        ------
        ValueError
            If a vector length differs from ``size``.
        """
        array = np.asarray(values, dtype=np.float64)
        push = self.net.Push
        for row in (array,) if array.ndim == 1 else array:
            push(to_net_doubles(row))

    @property
    def n(self) -> int:
        """Number of vectors pushed."""
        return int(self.net.N)

    @property
    def mean(self) -> np.ndarray:
        """Mean vector."""
        return to_numpy_matrix(self.net.Mean.ToArray()).ravel()

    @property
    def covariance(self) -> np.ndarray:
        """Co-moment matrix ``M2 + I``, not divided by the sample size."""
        return to_numpy_matrix(self.net.Covariance.ToArray())

    @property
    def sample_covariance(self) -> np.ndarray:
        """Regularized sample covariance ``(M2 + I) / (N - 1)``."""
        return to_numpy_matrix(self.net.SampleCovariance.ToArray())

    @property
    def population_covariance(self) -> np.ndarray:
        """Regularized population covariance ``(M2 + I) / N``."""
        return to_numpy_matrix(self.net.PopulationCovariance.ToArray())

    @property
    def sample_correlation(self) -> np.ndarray:
        """Correlation matrix of :attr:`sample_covariance`."""
        return to_numpy_matrix(self.net.SampleCorrelation.ToArray())

    @property
    def population_correlation(self) -> np.ndarray:
        """Correlation matrix of :attr:`population_covariance`."""
        return to_numpy_matrix(self.net.PopulationCorrelation.ToArray())


__all__ = [
    "Histogram", "HistogramBin", "RunningCovarianceMatrix", "RunningStatistics", "autocorrelation",
    "autocorrelation_confidence_interval", "borgonovo_delta", "box_cox_derivative", "box_cox_fit_lambda",
    "box_cox_inverse_transform", "box_cox_log_jacobian", "box_cox_log_likelihood", "box_cox_transform",
    "coefficient_of_variation", "column_linear_moments", "column_means", "column_percentiles",
    "column_standard_deviations", "column_variances", "covariance", "effective_sample_size",
    "entropy", "first_order_sobol", "five_number_summary", "geometric_mean", "harmonic_mean", "jackknife_sample",
    "jackknife_standard_error", "kendalls_tau", "kurtosis", "linear_coefficient_of_variation", "linear_moments",
    "maximum", "mean", "mean_standard_deviation", "mean_variance", "median", "minimum", "parallel_mean", "pawn",
    "pawn_median",
    "pearson", "percentile", "plotting_position_alpha", "plotting_position_percentile", "plotting_positions",
    "population_covariance", "population_standard_deviation", "population_variance", "product_moments", "ranks",
    "seven_number_summary", "skewness", "spearman", "standard_deviation", "sum",
    "variance", "weighted_plotting_positions", "wilson_score_interval", "yeo_johnson_derivative",
    "yeo_johnson_fit_lambda", "yeo_johnson_inverse_transform", "yeo_johnson_log_jacobian",
    "yeo_johnson_log_likelihood", "yeo_johnson_transform",
]

translate_public(globals())
