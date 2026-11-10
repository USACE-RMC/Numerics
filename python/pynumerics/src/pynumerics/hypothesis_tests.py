"""Hypothesis tests from ``Numerics.Data.Statistics.HypothesisTests`` and ``MultipleGrubbsBeckTest``.

Each function calls one Numerics member and returns its p-value (or statistic). The t-tests,
the Mann-Whitney, Mann-Kendall and Wald-Wolfowitz tests are two-sided; the F-test of two
samples is two-sided and the nested-model F-test, Jarque-Bera and Ljung-Box tests are upper-tail.
``two_sample_kolmogorov_smirnov`` returns the exact two-sample statistic without a p-value.
The Grubbs-Beck tests detect low outliers in flood series (Bulletin 17C).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics._util import namespace
from pynumerics.errors import translate_public
from pynumerics.marshal import to_net_doubles
from pynumerics.statistics import _sorted


def _tests() -> Any:
    """Return ``Numerics.Data.Statistics.HypothesisTests`` after starting the runtime."""
    return namespace("Numerics.Data.Statistics").HypothesisTests


def one_sample_t_test(sample: Any, population_mean: float = 0.0) -> float:
    """Return the two-sided p-value of the one-sample t-test (``HypothesisTests.OneSampleTtest``).

    Parameters
    ----------
    sample : array_like
        The sample.
    population_mean : float, optional
        The hypothesized mean; default 0.

    Returns
    -------
    float
        The p-value, as ``scipy.stats.ttest_1samp``.
    """
    return float(_tests().OneSampleTtest(to_net_doubles(sample), float(population_mean)))


def equal_variance_t_test(sample1: Any, sample2: Any) -> float:
    """Return the two-sided p-value of the pooled two-sample t-test (``HypothesisTests.EqualVarianceTtest``).

    Parameters
    ----------
    sample1, sample2 : array_like
        The samples.

    Returns
    -------
    float
        The p-value, as ``scipy.stats.ttest_ind(equal_var=True)``.
    """
    return float(_tests().EqualVarianceTtest(to_net_doubles(sample1), to_net_doubles(sample2)))


def unequal_variance_t_test(sample1: Any, sample2: Any) -> float:
    """Return the two-sided p-value of Welch's t-test (``HypothesisTests.UnequalVarianceTtest``).

    Parameters
    ----------
    sample1, sample2 : array_like
        The samples.

    Returns
    -------
    float
        The p-value, as ``scipy.stats.ttest_ind(equal_var=False)``.
    """
    return float(_tests().UnequalVarianceTtest(to_net_doubles(sample1), to_net_doubles(sample2)))


def paired_t_test(sample1: Any, sample2: Any) -> float:
    """Return the two-sided p-value of the paired t-test (``HypothesisTests.PairedTtest``).

    Parameters
    ----------
    sample1, sample2 : array_like
        Paired samples of the same length.

    Returns
    -------
    float
        The p-value, as ``scipy.stats.ttest_rel``.
    """
    return float(_tests().PairedTtest(to_net_doubles(sample1), to_net_doubles(sample2)))


def f_test(sample1: Any, sample2: Any) -> float:
    """Return the two-sided p-value of the F-test for equal variances (``HypothesisTests.Ftest``).

    Parameters
    ----------
    sample1, sample2 : array_like
        Samples of at least two values each.

    Returns
    -------
    float
        The p-value, as R's ``var.test``; 1 when both variances are zero.

    Raises
    ------
    ValueError
        If a sample has fewer than two values.
    """
    return float(_tests().Ftest(to_net_doubles(sample1), to_net_doubles(sample2)))


def f_test_models(sse_restricted: float, sse_full: float, df_restricted: int, df_full: int) -> tuple[float, float]:
    """Return the F statistic and p-value comparing nested models (``HypothesisTests.FtestModels``).

    ``F = ((SSE_R - SSE_F) / (df_R - df_F)) / (SSE_F / df_F)`` with the upper-tail p-value, as R's
    ``anova`` for nested linear models.

    Parameters
    ----------
    sse_restricted : float
        The sum of squared errors of the restricted model.
    sse_full : float
        The sum of squared errors of the full model.
    df_restricted : int
        The residual degrees of freedom of the restricted model.
    df_full : int
        The residual degrees of freedom of the full model, positive and fewer than ``df_restricted``.

    Returns
    -------
    f_statistic : float
        The F statistic.
    p_value : float
        The p-value; 1 when the statistic is zero or negative.

    Raises
    ------
    ValueError
        If the degrees of freedom are equal, ``df_full`` is not positive, or ``df_restricted`` is
        smaller than ``df_full``.
    """
    f_stat, p_value = _tests().FtestModels(float(sse_restricted), float(sse_full), int(df_restricted), int(df_full),
                                           0.0, 0.0)
    return float(f_stat), float(p_value)


def jarque_bera_test(sample: Any) -> float:
    """Return the p-value of the Jarque-Bera normality test (``HypothesisTests.JarqueBeraTest``).

    Parameters
    ----------
    sample : array_like
        The sample.

    Returns
    -------
    float
        The chi-squared(2) upper-tail p-value, as ``scipy.stats.jarque_bera``.
    """
    return float(_tests().JarqueBeraTest(to_net_doubles(sample)))


def wald_wolfowitz_test(sample: Any) -> float:
    """Return the two-sided p-value of the Wald-Wolfowitz test for independence (``HypothesisTests.WaldWolfowitzTest``).

    The circular lag-one serial product statistic with its permutation moments (Wald and
    Wolfowitz, 1943), referred to the standard normal.

    Parameters
    ----------
    sample : array_like
        The series, in time order.

    Returns
    -------
    float
        The p-value.
    """
    return float(_tests().WaldWolfowitzTest(to_net_doubles(sample)))


def ljung_box_test(sample: Any, max_lag: int | None = None) -> float:
    """Return the p-value of the Ljung-Box test of autocorrelation (``HypothesisTests.LjungBoxTest``).

    Parameters
    ----------
    sample : array_like
        The series, in time order.
    max_lag : int, optional
        The number of lags; default Numerics' ``floor(min(10 log10(n), n - 1))``.

    Returns
    -------
    float
        The chi-squared(``max_lag``) upper-tail p-value, as ``statsmodels`` ``acorr_ljungbox``;
        NaN when the autocorrelation is undefined.
    """
    return float(_tests().LjungBoxTest(to_net_doubles(sample), -1 if max_lag is None else int(max_lag)))


def mann_whitney_test(sample1: Any, sample2: Any) -> float:
    """Return the two-sided p-value of the Mann-Whitney test (``HypothesisTests.MannWhitneyTest``).

    The normal approximation with the tie correction and no continuity correction, as
    ``scipy.stats.mannwhitneyu(method="asymptotic", use_continuity=False)``.

    Parameters
    ----------
    sample1 : array_like
        The first sample, no longer than ``sample2`` and of more than 3 values.
    sample2 : array_like
        The second sample, of more than 3 values; together more than 20.

    Returns
    -------
    float
        The p-value.

    Raises
    ------
    ValueError
        If ``sample1`` is longer than ``sample2``, a sample has 3 or fewer values, or the two
        have 20 or fewer together.
    """
    return float(_tests().MannWhitneyTest(to_net_doubles(sample1), to_net_doubles(sample2)))


def mann_kendall_test(sample: Any) -> float:
    """Return the two-sided p-value of the Mann-Kendall trend test (``HypothesisTests.MannKendallTest``).

    The normal approximation of S with the tie-corrected variance and the continuity correction.

    Parameters
    ----------
    sample : array_like
        The series, in time order, of at least 10 values.

    Returns
    -------
    float
        The p-value.

    Raises
    ------
    ValueError
        If the series has fewer than 10 values.
    """
    return float(_tests().MannKendallTest(to_net_doubles(sample)))


def linear_trend_test(indices: Any, sample: Any) -> float:
    """Return the two-sided p-value of the slope of a linear regression (``HypothesisTests.LinearTrendTest``).

    Parameters
    ----------
    indices : array_like
        The time indices.
    sample : array_like
        The series, aligned with ``indices``.

    Returns
    -------
    float
        The p-value of the slope t statistic, as an ordinary least squares fit with intercept.

    Raises
    ------
    ValueError
        If the lengths differ.
    """
    return float(_tests().LinearTrendTest(to_net_doubles(indices), to_net_doubles(sample)))


def hartigan_dip(sample: Any) -> float:
    """Return Hartigan's dip statistic (``HypothesisTests.HartiganDip``, algorithm AS 217).

    Parameters
    ----------
    sample : array_like
        Finite values, in any order.

    Returns
    -------
    float
        The dip, as R's ``diptest::dip``; between 1/(2n) and 1/4, 0 for an empty sample.

    Raises
    ------
    ValueError
        If a value is not finite.
    """
    return float(_tests().HartiganDip(to_net_doubles(sample)))


def hartigan_dip_test(sample: Any, realizations: int = 2000, seed: int = 12345) -> float:
    """Return the simulated p-value of Hartigan's dip test of unimodality (``HypothesisTests.HartiganDipTest``).

    Parameters
    ----------
    sample : array_like
        Finite values.
    realizations : int, optional
        The number of simulated uniform samples; default 2000.
    seed : int, optional
        The simulation seed; default 12345.

    Returns
    -------
    float
        The share of uniform samples with a dip at least the sample's; 1 for three or fewer values.

    Raises
    ------
    ValueError
        If a value is not finite or ``realizations`` is below one.
    """
    return float(_tests().HartiganDipTest(to_net_doubles(sample), int(realizations), int(seed)))


def unimodality_test(sample: Any, realizations: int | None = None, seed: int = 12345) -> float:
    """Return the p-value of the Gaussian mixture likelihood-ratio unimodality test (``UnimodalityTest``).

    One Gaussian against a two-component mixture: without ``realizations`` the asymptotic
    chi-squared p-value, otherwise the parametric bootstrap p-value.

    Parameters
    ----------
    sample : array_like
        At least 10 values.
    realizations : int, optional
        The number of bootstrap samples; None (default) for the asymptotic p-value.
    seed : int, optional
        The bootstrap seed; default 12345.

    Returns
    -------
    float
        The p-value, or NaN when the mixture fit fails.

    Raises
    ------
    ValueError
        If the sample has fewer than 10 values or ``realizations`` is below one.
    """
    import System  # noqa: PLC0415 - the runtime is running once _tests returns
    from System.Collections.Generic import IList  # noqa: PLC0415

    group = _tests().UnimodalityTest
    values = to_net_doubles(sample)
    if realizations is None:
        return float(group.Overloads[IList[System.Double]](values))
    member = group.Overloads[IList[System.Double], System.Int32, System.Int32]
    return float(member(values, int(realizations), int(seed)))


def multiple_grubbs_beck_test(flows: Any) -> int:
    """Return the number of low outliers by the multiple Grubbs-Beck test (``MultipleGrubbsBeckTest.Function``).

    The generalized Grubbs-Beck test of Cohn et al. (2013) used by Bulletin 17C, on the base-10
    logarithms of the flows, converted from the PeakfqSA source.

    Parameters
    ----------
    flows : array_like
        Annual peak flows in real space; a value below 1e-88 counts as 1e-88.

    Returns
    -------
    int
        The number of potentially influential low flows (the smallest values to censor).
    """
    stats = namespace("Numerics.Data.Statistics")
    return int(stats.MultipleGrubbsBeckTest.Function(to_net_doubles(flows)))


def grubbs_beck_test(sample: Any) -> tuple[float, float]:
    """Return the outlier thresholds of the original Grubbs-Beck test (``MultipleGrubbsBeckTest.GrubbsBeckTest``).

    The test at the 10 percent level on the natural logarithms: ``exp(mean +- K_N s)`` with the
    Pilon et al. (1985) approximation ``K_N = -3.62201 + 6.28446 N^(1/4) - 2.49835 N^(1/2)
    + 0.491436 N^(3/4) - 0.037911 N``.

    Parameters
    ----------
    sample : array_like
        Positive flows.

    Returns
    -------
    high : float
        Values above it are high outliers.
    low : float
        Values below it are low outliers.

    Raises
    ------
    ValueError
        If a value is not positive.
    """
    stats = namespace("Numerics.Data.Statistics")
    high, low = stats.MultipleGrubbsBeckTest.GrubbsBeckTest(to_net_doubles(np.asarray(sample, dtype=np.float64)),
                                                             0.0, 0.0)
    return float(high), float(low)


def two_sample_kolmogorov_smirnov(first: Any, second: Any, data_is_sorted: bool = False) -> float:
    """Return the exact two-sample Kolmogorov-Smirnov statistic (``HypothesisTests.TwoSampleKolmogorovSmirnov``).

    ``sup |F1(t) - F2(t)|`` over the merged values, ties counted as single jumps. The statistic is
    ``max |i*n2 - j*n1| / (n1*n2)`` formed from integer counts, so it is the correctly rounded
    rational; ``scipy.stats.ks_2samp(...).statistic`` agrees to an ulp or so, and bit for bit
    where scipy rounds it from the integer count. Numerics requires ascending inputs and checks the
    order; unless ``data_is_sorted`` is True the samples are sorted first.

    Parameters
    ----------
    first, second : array_like
        Nonempty, finite samples.
    data_is_sorted : bool, optional
        Whether both samples are already ascending, which skips the sort; default False.

    Returns
    -------
    float
        The statistic, in [0, 1].

    Raises
    ------
    ValueError
        If a sample is empty, holds a non-finite value, or is not ascending when ``data_is_sorted``
        is True.
    """
    return float(_tests().TwoSampleKolmogorovSmirnov(_sorted(first, data_is_sorted), _sorted(second, data_is_sorted)))


__all__ = [
    "equal_variance_t_test", "f_test", "f_test_models", "grubbs_beck_test", "hartigan_dip", "hartigan_dip_test",
    "jarque_bera_test", "linear_trend_test", "ljung_box_test", "mann_kendall_test", "mann_whitney_test",
    "multiple_grubbs_beck_test", "one_sample_t_test", "paired_t_test", "two_sample_kolmogorov_smirnov",
    "unequal_variance_t_test", "unimodality_test", "wald_wolfowitz_test",
]

translate_public(globals())
