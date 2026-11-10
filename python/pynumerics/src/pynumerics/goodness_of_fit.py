"""Goodness-of-fit measures from ``Numerics.Data.Statistics.GoodnessOfFit``.

Information criteria (AIC, AICc, BIC and their weights), error metrics (RMSE, MSE, MAE, MAPE,
sMAPE), hydrologic efficiencies (Nash-Sutcliffe, Kling-Gupta, PBIAS, RSR, the Willmott indices
of agreement, volumetric efficiency), binary classification scores, and the
Kolmogorov-Smirnov, Anderson-Darling and chi-squared statistics of a sample against a
distribution. Each function calls one Numerics member. The paired metrics take the observed
values first and the modeled values second, as Numerics does.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics._util import cached_overload, namespace
from pynumerics.errors import translate_public
from pynumerics.marshal import to_net_doubles, to_net_matrix, to_numpy_doubles
from pynumerics.statistics import _sorted


def _gof() -> Any:
    """Return ``Numerics.Data.Statistics.GoodnessOfFit`` after starting the runtime."""
    return namespace("Numerics.Data.Statistics").GoodnessOfFit


def _distribution(model: Any) -> Any:
    """Return the .NET ``UnivariateDistributionBase`` behind a distribution, or None for anything else."""
    import clr  # noqa: PLC0415 - the runtime is running once _gof returns

    net = getattr(model, "net", model)
    base = clr.GetClrType(namespace("Numerics.Distributions").UnivariateDistributionBase)
    try:
        return net if base.IsInstanceOfType(net) else None
    except TypeError:
        return None


def _required_distribution(model: Any) -> Any:
    """Return the .NET distribution behind ``model``, or raise TypeError."""
    _gof()
    net = _distribution(model)
    if net is None:
        raise TypeError(f"Expected a univariate distribution; got {type(model).__name__}.")
    return net


def _paired(name: str, observed: Any, modeled: Any) -> float:
    """Call a paired ``GoodnessOfFit`` metric on observed and modeled values."""
    return float(getattr(_gof(), name)(to_net_doubles(observed), to_net_doubles(modeled)))


def aic(number_of_parameters: int, log_likelihood: float) -> float:
    """Return the Akaike information criterion ``2k - 2 ln L`` (``GoodnessOfFit.AIC``).

    Parameters
    ----------
    number_of_parameters : int
        The number of fitted parameters k.
    log_likelihood : float
        The maximized log-likelihood.

    Returns
    -------
    float
        The AIC; lower is preferred.
    """
    return float(_gof().AIC(int(number_of_parameters), float(log_likelihood)))


def aicc(sample_size: int, number_of_parameters: int, log_likelihood: float) -> float:
    """Return the small-sample corrected AIC ``AIC + 2k(k + 1) / (n - k - 1)`` (``GoodnessOfFit.AICc``).

    Parameters
    ----------
    sample_size : int
        The sample size n.
    number_of_parameters : int
        The number of fitted parameters k.
    log_likelihood : float
        The maximized log-likelihood.

    Returns
    -------
    float
        The AICc; lower is preferred.
    """
    return float(_gof().AICc(int(sample_size), int(number_of_parameters), float(log_likelihood)))


def bic(sample_size: int, number_of_parameters: int, log_likelihood: float) -> float:
    """Return the Bayesian information criterion ``k ln n - 2 ln L`` (``GoodnessOfFit.BIC``).

    Parameters
    ----------
    sample_size : int
        The sample size n.
    number_of_parameters : int
        The number of fitted parameters k.
    log_likelihood : float
        The maximized log-likelihood.

    Returns
    -------
    float
        The BIC; lower is preferred.
    """
    return float(_gof().BIC(int(sample_size), int(number_of_parameters), float(log_likelihood)))


def aic_weights(aic_values: Any) -> np.ndarray:
    """Return the Akaike weights of a set of models (``GoodnessOfFit.AICWeights``).

    Parameters
    ----------
    aic_values : array_like
        One AIC per model.

    Returns
    -------
    numpy.ndarray
        ``exp(-(AIC_i - AIC_min) / 2)`` normalized to sum to one.
    """
    return to_numpy_doubles(_gof().AICWeights(to_net_doubles(aic_values)))


def rmse(observed: Any, model: Any, k: int = 0, plotting_positions: Any = None) -> float:
    """Return the root mean square error (``GoodnessOfFit.RMSE``).

    With modeled values, ``sqrt(sum (m - o)^2 / (n - k))``. With a distribution, the observed
    values are compared with its quantiles at ``plotting_positions`` (Weibull positions of the
    sorted observations by default), and k is the distribution's parameter count.

    Parameters
    ----------
    observed : array_like
        The observed values.
    model : array_like or UnivariateDistribution
        The modeled values, one per observation, or a fitted distribution.
    k : int, optional
        The number of fitted parameters for modeled values; default 0. Ignored for a distribution.
    plotting_positions : array_like, optional
        The non-exceedance probabilities of the observations, in their order, for a distribution.

    Returns
    -------
    float
        The RMSE.

    Raises
    ------
    ValueError
        If the lengths differ or ``k`` is not in [0, n).
    """
    gof = _gof()
    net = _distribution(model)
    if net is None:
        if plotting_positions is not None:
            raise ValueError("plotting_positions applies only when model is a distribution.")
        return float(gof.RMSE(to_net_doubles(observed), to_net_doubles(model), int(k)))
    if plotting_positions is None:
        return float(gof.RMSE(to_net_doubles(observed), net))
    return float(gof.RMSE(to_net_doubles(observed), to_net_doubles(plotting_positions), net))


def rmse_weights(rmse_values: Any) -> np.ndarray:
    """Return inverse-MSE model weights (``GoodnessOfFit.RMSEWeights``).

    Parameters
    ----------
    rmse_values : array_like
        One RMSE per model.

    Returns
    -------
    numpy.ndarray
        ``RMSE_i^-2`` normalized to sum to one.
    """
    return to_numpy_doubles(_gof().RMSEWeights(to_net_doubles(rmse_values)))


def mse(observed: Any, modeled: Any) -> float:
    """Return the mean squared error (``GoodnessOfFit.MSE``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The MSE.
    """
    return _paired("MSE", observed, modeled)


def mae(observed: Any, modeled: Any) -> float:
    """Return the mean absolute error (``GoodnessOfFit.MAE``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The MAE.
    """
    return _paired("MAE", observed, modeled)


def mape(observed: Any, modeled: Any) -> float:
    """Return the mean absolute percentage error ``100/n sum |o - m| / |o|`` (``GoodnessOfFit.MAPE``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The MAPE in percent.
    """
    return _paired("MAPE", observed, modeled)


def smape(observed: Any, modeled: Any) -> float:
    """Return the symmetric MAPE ``100/n sum 2|o - m| / (|o| + |m|)`` (``GoodnessOfFit.sMAPE``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The sMAPE in percent, in [0, 200].
    """
    return _paired("sMAPE", observed, modeled)


def nash_sutcliffe_efficiency(observed: Any, modeled: Any) -> float:
    """Return the Nash-Sutcliffe efficiency ``1 - sum (o - m)^2 / sum (o - mean o)^2`` (``NashSutcliffeEfficiency``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The NSE, at most one.
    """
    return _paired("NashSutcliffeEfficiency", observed, modeled)


def log_nash_sutcliffe_efficiency(observed: Any, modeled: Any, epsilon: float | None = None) -> float:
    """Return the Nash-Sutcliffe efficiency of ``ln(x + epsilon)`` (``GoodnessOfFit.LogNashSutcliffeEfficiency``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.
    epsilon : float, optional
        The non-negative offset; default the observed mean over 100 (Pushpalatha et al., 2012).

    Returns
    -------
    float
        The log NSE.
    """
    eps = -1.0 if epsilon is None else float(epsilon)
    return float(_gof().LogNashSutcliffeEfficiency(to_net_doubles(observed), to_net_doubles(modeled), eps))


def kling_gupta_efficiency(observed: Any, modeled: Any) -> float:
    """Return the Kling-Gupta efficiency of Gupta et al. (2009) (``GoodnessOfFit.KlingGuptaEfficiency``).

    ``1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)`` with r the Pearson correlation,
    alpha the ratio of standard deviations and beta the ratio of means (modeled over observed).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The KGE; -10 when either series is constant, NaN for a single pair.
    """
    return _paired("KlingGuptaEfficiency", observed, modeled)


def kling_gupta_efficiency_mod(observed: Any, modeled: Any) -> float:
    """Return the modified Kling-Gupta efficiency of Kling et al. (2012) (``GoodnessOfFit.KlingGuptaEfficiencyMod``).

    As :func:`kling_gupta_efficiency` with alpha replaced by the ratio of coefficients of variation.

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The KGE'; -10 when either series is constant, NaN for a single pair.
    """
    return _paired("KlingGuptaEfficiencyMod", observed, modeled)


def pbias(observed: Any, modeled: Any) -> float:
    """Return the percent bias ``100 sum (o - m) / sum o`` of Moriasi et al. (2007) (``GoodnessOfFit.PBIAS``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The PBIAS in percent; positive when the model underestimates the observed total and negative
        when it overestimates. R's hydroGOF ``pbias`` reports the opposite sign.
    """
    return _paired("PBIAS", observed, modeled)


def rsr(observed: Any, modeled: Any) -> float:
    """Return ``sqrt(sum (o - m)^2 / sum (o - mean o)^2)``, the RMSE over the observed spread (``GoodnessOfFit.RSR``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The RSR, zero for a perfect model.
    """
    return _paired("RSR", observed, modeled)


def r_squared(observed: Any, modeled: Any) -> float:
    """Return the coefficient of determination, the squared Pearson correlation (``GoodnessOfFit.RSquared``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        R squared, in [0, 1].
    """
    return _paired("RSquared", observed, modeled)


def index_of_agreement(observed: Any, modeled: Any) -> float:
    """Return Willmott's (1981) index of agreement d (``GoodnessOfFit.IndexOfAgreement``).

    ``1 - sum (m - o)^2 / sum (|m - mean o| + |o - mean o|)^2``.

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The index d, in [0, 1].
    """
    return _paired("IndexOfAgreement", observed, modeled)


def modified_index_of_agreement(observed: Any, modeled: Any) -> float:
    """Return the modified index of agreement d1 of Willmott et al. (1985) (``ModifiedIndexOfAgreement``).

    ``1 - sum |m - o| / sum (|m - mean o| + |o - mean o|)``.

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The index d1, in [0, 1].
    """
    return _paired("ModifiedIndexOfAgreement", observed, modeled)


def refined_index_of_agreement(observed: Any, modeled: Any) -> float:
    """Return the refined index of agreement dr of Willmott et al. (2012) (``RefinedIndexOfAgreement``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        The index dr, in [-1, 1].
    """
    return _paired("RefinedIndexOfAgreement", observed, modeled)


def volumetric_efficiency(observed: Any, modeled: Any) -> float:
    """Return the volumetric efficiency ``1 - sum |m - o| / sum |o|`` of Criss and Winston (2008).

    Calls ``GoodnessOfFit.VolumetricEfficiency``.

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled values of the same length.

    Returns
    -------
    float
        VE, at most one.
    """
    return _paired("VolumetricEfficiency", observed, modeled)


def accuracy(observed: Any, modeled: Any) -> float:
    """Return the percentage of exact matches (``GoodnessOfFit.Accuracy``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled labels of the same length.

    Returns
    -------
    float
        The accuracy in percent.
    """
    return _paired("Accuracy", observed, modeled)


def precision(observed: Any, modeled: Any) -> float:
    """Return the precision ``TP / (TP + FP)`` of binary labels (``GoodnessOfFit.Precision``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled 0/1 labels of the same length.

    Returns
    -------
    float
        The precision, in [0, 1].
    """
    return _paired("Precision", observed, modeled)


def recall(observed: Any, modeled: Any) -> float:
    """Return the recall ``TP / (TP + FN)`` of binary labels (``GoodnessOfFit.Recall``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled 0/1 labels of the same length.

    Returns
    -------
    float
        The recall, in [0, 1].
    """
    return _paired("Recall", observed, modeled)


def f1_score(observed: Any, modeled: Any) -> float:
    """Return the F1 score ``2TP / (2TP + FP + FN)`` of binary labels (``GoodnessOfFit.F1Score``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled 0/1 labels of the same length.

    Returns
    -------
    float
        The F1 score, in [0, 1].
    """
    return _paired("F1Score", observed, modeled)


def specificity(observed: Any, modeled: Any) -> float:
    """Return the specificity ``TN / (TN + FP)`` of binary labels (``GoodnessOfFit.Specificity``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled 0/1 labels of the same length.

    Returns
    -------
    float
        The specificity, in [0, 1].
    """
    return _paired("Specificity", observed, modeled)


def balanced_accuracy(observed: Any, modeled: Any) -> float:
    """Return the mean of recall and specificity of binary labels (``GoodnessOfFit.BalancedAccuracy``).

    Parameters
    ----------
    observed, modeled : array_like
        Observed and modeled 0/1 labels of the same length.

    Returns
    -------
    float
        The balanced accuracy, in [0, 1].
    """
    return _paired("BalancedAccuracy", observed, modeled)


def _column_uniform(member: str, sorted_samples: Any) -> np.ndarray:
    """Call a batched ``GoodnessOfFit`` statistic against Uniform(0, 1) on one sample per column.

    The matrix crosses in one copy, and the member is looked up on first use only.
    """
    bound = cached_overload(f"GoodnessOfFit.{member}(double[,])", lambda: getattr(_gof(), member))
    return to_numpy_doubles(bound(to_net_matrix(sorted_samples)))


def column_kolmogorov_smirnov_uniform(sorted_samples: Any) -> np.ndarray:
    """Return the Kolmogorov-Smirnov D of each column against Uniform(0, 1).

    Wraps ``GoodnessOfFit.ColumnKolmogorovSmirnovUniform``.

    Each value is bitwise :func:`kolmogorov_smirnov` of the column against ``Uniform(0, 1)``, as
    ``scipy.stats.kstest(u, "uniform").statistic``, for all columns in one call.

    Parameters
    ----------
    sorted_samples : array_like
        An (n, m) array holding m samples of n values, one per column, each ascending, in [0, 1].

    Returns
    -------
    numpy.ndarray
        One statistic per column.

    Raises
    ------
    ValueError
        If there are no rows, or a value is NaN or outside [0, 1], or a column is not ascending.
    """
    return _column_uniform("ColumnKolmogorovSmirnovUniform", sorted_samples)


def column_anderson_darling_uniform(sorted_samples: Any) -> np.ndarray:
    """Return the Anderson-Darling A^2 of each column against Uniform(0, 1).

    Wraps ``GoodnessOfFit.ColumnAndersonDarlingUniform``.

    ``A^2 = -n - n^-1 sum_i (2i - 1)(ln u_(i) + ln(1 - u_(n+1-i)))`` for each column, bitwise
    :func:`anderson_darling` against ``Uniform(0, 1)``, for all columns in one call. A column with a
    value at 0 or 1 gives infinity.

    Parameters
    ----------
    sorted_samples : array_like
        An (n, m) array holding m samples of n values, one per column, each ascending, in [0, 1].

    Returns
    -------
    numpy.ndarray
        One statistic per column.

    Raises
    ------
    ValueError
        If there are no rows, or a value is NaN or outside [0, 1], or a column is not ascending.
    """
    return _column_uniform("ColumnAndersonDarlingUniform", sorted_samples)


def kolmogorov_smirnov(sample: Any, distribution: Any, data_is_sorted: bool = False) -> float:
    """Return the one-sample Kolmogorov-Smirnov statistic against a distribution (``GoodnessOfFit.KolmogorovSmirnov``).

    ``D = max_i max(F(x_(i)) - (i - 1)/n, i/n - F(x_(i)))``, as ``scipy.stats.kstest(...).statistic``.
    Unless ``data_is_sorted`` is True the sample is sorted first.

    Parameters
    ----------
    sample : array_like
        The observations.
    distribution : UnivariateDistribution
        The hypothesized distribution.
    data_is_sorted : bool, optional
        Whether ``sample`` is already ascending, which skips the sort; default False.

    Returns
    -------
    float
        The statistic D.

    Raises
    ------
    TypeError
        If ``distribution`` is not a univariate distribution.
    ValueError
        If the sample is empty.
    """
    net = _required_distribution(distribution)
    return float(_gof().KolmogorovSmirnov(_sorted(sample, data_is_sorted), net, True))


def anderson_darling(sample: Any, distribution: Any, data_is_sorted: bool = False) -> float:
    """Return the Anderson-Darling statistic against a distribution (``GoodnessOfFit.AndersonDarling``).

    ``A^2 = -n - (1/n) sum (2i - 1) (ln F(x_(i)) + ln(1 - F(x_(n+1-i))))`` with the distribution's
    log CDF and log survival function. Unless ``data_is_sorted`` is True the sample is sorted first.

    Parameters
    ----------
    sample : array_like
        The observations.
    distribution : UnivariateDistribution
        The hypothesized distribution.
    data_is_sorted : bool, optional
        Whether ``sample`` is already ascending, which skips the sort; default False.

    Returns
    -------
    float
        The statistic A squared.

    Raises
    ------
    TypeError
        If ``distribution`` is not a univariate distribution.
    ValueError
        If the sample is empty.
    """
    net = _required_distribution(distribution)
    return float(_gof().AndersonDarling(_sorted(sample, data_is_sorted), net, True))


def chi_squared(sample: Any, distribution: Any, data_is_sorted: bool = False) -> float:
    """Return the chi-squared statistic against a distribution (``GoodnessOfFit.ChiSquared``).

    The sample is binned with the default :class:`pynumerics.statistics.Histogram`, and the
    statistic is ``sum (O_j - E_j)^2 / E_j`` with ``E_j = n (F(upper_j) - F(lower_j))``. Unless
    ``data_is_sorted`` is True the sample is sorted first.

    Parameters
    ----------
    sample : array_like
        The observations.
    distribution : UnivariateDistribution
        The hypothesized distribution.
    data_is_sorted : bool, optional
        Whether ``sample`` is already ascending, which skips the sort; default False.

    Returns
    -------
    float
        The statistic.

    Raises
    ------
    TypeError
        If ``distribution`` is not a univariate distribution.
    ValueError
        If the sample is empty.
    """
    net = _required_distribution(distribution)
    return float(_gof().ChiSquared(_sorted(sample, data_is_sorted), net))


__all__ = [
    "accuracy", "aic", "aic_weights", "aicc", "anderson_darling", "balanced_accuracy", "bic", "chi_squared",
    "column_anderson_darling_uniform", "column_kolmogorov_smirnov_uniform",
    "f1_score", "index_of_agreement", "kling_gupta_efficiency", "kling_gupta_efficiency_mod", "kolmogorov_smirnov",
    "log_nash_sutcliffe_efficiency", "mae", "mape", "modified_index_of_agreement", "mse", "nash_sutcliffe_efficiency",
    "pbias", "precision", "r_squared", "recall", "refined_index_of_agreement", "rmse", "rmse_weights", "rsr",
    "smape", "specificity", "volumetric_efficiency",
]

translate_public(globals())
