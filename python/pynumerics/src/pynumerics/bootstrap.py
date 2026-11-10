"""Uncertainty analysis from ``Numerics.Distributions``: ensembles and the parametric bootstrap.

:class:`UncertaintyAnalysisResults` wraps the ensemble summaries (Type-7 percentile
bands, the expected-probability curve, the mode curve); :class:`BootstrapAnalysis`
wraps the parametric resample-and-refit bootstrap of one distribution.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from pynumerics._util import _nd, net_enum
from pynumerics.distributions import UnivariateDistribution, _net_property, _parse_xml, wrap
from pynumerics.enums import EstimationMethod
from pynumerics.errors import translate_public
from pynumerics.marshal import _parameter_set_arrays, to_net_doubles, to_numpy_doubles, to_numpy_matrix


def _members(distributions: Any) -> Any:
    """Build a ``UnivariateDistributionBase[]``; None entries stay null (Numerics treats them as missing).

    A .NET array (for example from ``UnivariateDistribution.parameter_set_distributions``) is passed
    through unchanged.
    """
    nd = _nd()
    import clr  # noqa: PLC0415
    import System  # noqa: PLC0415

    if isinstance(distributions, System.Array) and clr.GetClrType(nd.UnivariateDistributionBase).IsAssignableFrom(
            distributions.GetType().GetElementType()):
        return distributions
    return System.Array[nd.UnivariateDistributionBase]([None if d is None else d.net for d in distributions])


class UncertaintyAnalysisResults:
    """Ensemble summaries (``Numerics.Distributions.UncertaintyAnalysisResults``).

    Construct empty and call the ``process_*`` members, or use :meth:`from_ensemble`
    for the full constructor. Probabilities are non-exceedance probabilities.

    Parameters
    ----------
    net : Numerics.Distributions.UncertaintyAnalysisResults, optional
        An existing .NET results object to wrap; None, the default, creates an empty one.

    Attributes
    ----------
    net : Numerics.Distributions.UncertaintyAnalysisResults
        The wrapped .NET object.
    """

    def __init__(self, net: Any = None) -> None:
        self.net = _nd().UncertaintyAnalysisResults() if net is None else net

    @classmethod
    def from_ensemble(cls, parent: UnivariateDistribution, sampled: Sequence[UnivariateDistribution | None],
                      probabilities: Any, alpha: float = 0.1, min_probability: float = 0.001,
                      max_probability: float = 1 - 1e-9,
                      record_parameter_sets: bool = False) -> UncertaintyAnalysisResults:
        """Run the full constructor: mode curve, confidence intervals and mean curve.

        Parameters
        ----------
        parent : UnivariateDistribution
            The point-estimate distribution (mode curve).
        sampled : sequence of UnivariateDistribution or System.Array
            Ensemble members; None (null) members are missing.
        probabilities : array_like
            Non-exceedance probabilities at which to summarize.
        alpha : float, optional
            Two-sided band level; the band is [α/2, 1 − α/2]; default 0.1.
        min_probability, max_probability : float, optional
            Range of the mean-curve value grid.
        record_parameter_sets : bool, optional
            Store the members' parameter sets; default False.

        Returns
        -------
        UncertaintyAnalysisResults
            The populated results.
        """
        net = _nd().UncertaintyAnalysisResults(parent.net, _members(sampled), to_net_doubles(probabilities),
                                               float(alpha), float(min_probability), float(max_probability),
                                               bool(record_parameter_sets))
        return cls(net)

    def process_mode_curve(self, parent: UnivariateDistribution, probabilities: Any) -> np.ndarray:
        """Compute the mode curve ``parent.InverseCDF(probabilities)`` (``ProcessModeCurve``).

        Parameters
        ----------
        parent : UnivariateDistribution
            The point-estimate distribution.
        probabilities : array_like
            Non-exceedance probabilities.

        Returns
        -------
        numpy.ndarray
            The mode curve.
        """
        self.net.ProcessModeCurve(parent.net, to_net_doubles(probabilities))
        return self.mode_curve

    def process_confidence_intervals(self, sampled: Sequence[UnivariateDistribution | None], probabilities: Any,
                                     alpha: float = 0.1) -> np.ndarray:
        """Compute the Type-7 percentile band of member quantiles (``ProcessConfidenceIntervals``).

        Non-finite member quantiles are dropped before the percentiles.

        Parameters
        ----------
        sampled : sequence of UnivariateDistribution or System.Array
            Ensemble members; None (null) members are missing.
        probabilities : array_like
            Non-exceedance probabilities.
        alpha : float, optional
            Two-sided band level; default 0.1.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), 2)``: lower and upper limits.
        """
        self.net.ProcessConfidenceIntervals(_members(sampled), to_net_doubles(probabilities), float(alpha))
        return self.confidence_intervals

    def process_mean_curve(self, sampled: Sequence[UnivariateDistribution | None], probabilities: Any,
                           min_probability: float = 0.001, max_probability: float = 1 - 1e-9,
                           weights: Any = None) -> np.ndarray:
        """Compute the expected-probability curve E[F] (``ProcessMeanCurve``).

        Parameters
        ----------
        sampled : sequence of UnivariateDistribution or System.Array
            Ensemble members; None (null) members are missing.
        probabilities : array_like
            Non-exceedance probabilities at which the curve is read.
        min_probability, max_probability : float, optional
            Range of the value grid.
        weights : array_like, optional
            Non-negative member weights (for example logic-tree branch weights), one per
            member; the weighted overload averages ``sum(w F) / sum(w)``. Equal weights give
            bitwise the unweighted curve. Default None, the unweighted mean.

        Returns
        -------
        numpy.ndarray
            The mean curve.
        """
        if weights is None:
            self.net.ProcessMeanCurve(_members(sampled), to_net_doubles(probabilities), float(min_probability),
                                      float(max_probability))
        else:
            self.net.ProcessMeanCurve(_members(sampled), to_net_doubles(probabilities), to_net_doubles(weights),
                                      float(min_probability), float(max_probability))
        return self.mean_curve

    @property
    def confidence_intervals(self) -> np.ndarray | None:
        """The band, shape ``(n, 2)``."""
        ci = self.net.ConfidenceIntervals
        return None if ci is None else to_numpy_matrix(ci)

    @property
    def mean_curve(self) -> np.ndarray | None:
        """The expected-probability curve."""
        curve = self.net.MeanCurve
        return None if curve is None else to_numpy_doubles(curve)

    @property
    def mode_curve(self) -> np.ndarray | None:
        """The mode curve."""
        curve = self.net.ModeCurve
        return None if curve is None else to_numpy_doubles(curve)

    aic = _net_property("AIC", "The Akaike information criterion of the parent fit (``AIC``); settable.")
    bic = _net_property("BIC", "The Bayesian information criterion of the parent fit (``BIC``); settable.")
    dic = _net_property("DIC", "The deviance information criterion (``DIC``); settable.")
    rmse = _net_property("RMSE", "The root mean square error of the parent fit (``RMSE``); settable.")
    erl = _net_property("ERL", "The effective record length (``ERL``); settable.")

    @property
    def parent_distribution(self) -> UnivariateDistribution | None:
        """The point-estimate distribution, wrapped (``ParentDistribution``); settable."""
        parent = self.net.ParentDistribution
        return None if parent is None else wrap(parent)

    @parent_distribution.setter
    def parent_distribution(self, value: UnivariateDistribution | None) -> None:
        self.net.ParentDistribution = None if value is None else value.net

    @property
    def parameter_sets(self) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
        """The recorded member parameter sets as ``(values, fitness, weight)`` (``ParameterSets``), or None.

        ``values`` has one row per member (NaN for a missing member); ``fitness`` and ``weight`` hold
        one entry per member.
        """
        sets = self.net.ParameterSets
        return None if sets is None else _parameter_set_arrays(sets)

    def process_parameter_sets(self, sampled: Sequence[UnivariateDistribution | None]) -> None:
        """Record the members' parameter sets (``ProcessParameterSets``).

        Parameters
        ----------
        sampled : sequence of UnivariateDistribution or System.Array
            Ensemble members; None (null) members are recorded as NaN rows.
        """
        self.net.ProcessParameterSets(_members(sampled))

    def to_bytes(self) -> bytes:
        """Serialize the results as Numerics does (``ToByteArray``, UTF-8 JSON).

        Returns
        -------
        bytes
            The serialized results; :meth:`from_bytes` reads them back.
        """
        return bytes(_nd().UncertaintyAnalysisResults.ToByteArray(self.net))

    @classmethod
    def from_bytes(cls, data: bytes) -> UncertaintyAnalysisResults:
        """Read results serialized by :meth:`to_bytes` (``FromByteArray``).

        Parameters
        ----------
        data : bytes
            The serialized results.

        Returns
        -------
        UncertaintyAnalysisResults
            The results.
        """
        nd = _nd()
        import System  # noqa: PLC0415 - the runtime is running (_nd started it)

        return cls(nd.UncertaintyAnalysisResults.FromByteArray(System.Array[System.Byte](bytes(data))))

    def to_xml(self) -> str:
        """Return the results as XML text (``ToXElement``).

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.ToXElement().ToString())

    @classmethod
    def from_xml(cls, xml: str) -> UncertaintyAnalysisResults:
        """Read results written by :meth:`to_xml` (``FromXElement``).

        Parameters
        ----------
        xml : str
            The element as text.

        Returns
        -------
        UncertaintyAnalysisResults
            The results.
        """
        return cls(_nd().UncertaintyAnalysisResults.FromXElement(_parse_xml(xml)))


class BootstrapAnalysis:
    """Parametric bootstrap of one distribution (``Numerics.Distributions.BootstrapAnalysis``).

    Parameters
    ----------
    distribution : UnivariateDistribution
        The parent distribution; it must support bootstrapping (``IBootstrappable``).
    method : EstimationMethod or str
        Refit method for each replicate.
    sample_size : int
        Size of each resampled data set.
    replications : int, optional
        Number of replicates; default 10,000.
    seed : int, optional
        Seed; default 12345.

    Attributes
    ----------
    net : Numerics.Distributions.BootstrapAnalysis
        The wrapped .NET object.
    """

    def __init__(self, distribution: UnivariateDistribution, method: EstimationMethod | str, sample_size: int,
                 replications: int = 10000, seed: int = 12345) -> None:
        nd = _nd()
        self.net = nd.BootstrapAnalysis(distribution.net, net_enum(nd.ParameterEstimationMethod, method),
                                        int(sample_size), int(replications), int(seed))
        self._distributions: Any = None

    @property
    def failed_replications(self) -> int:
        """Replicates whose refit failed (set after :meth:`distributions`)."""
        return int(self.net.FailedReplications)

    def _net_distributions(self) -> Any:
        """Return the cached .NET replicate array, computing it once."""
        if self._distributions is None:
            self._distributions = self.net.Distributions()
        return self._distributions

    def members(self) -> Any:
        """Run (once) and return the replicates as a ``UnivariateDistributionBase[]`` with nulls for failures.

        The array stays on the .NET side, ready for the members that take one (for example
        ``UncertaintyAnalysisResults.from_ensemble``), so a large ensemble is never wrapped
        replicate by replicate in Python.

        Returns
        -------
        System.Array[Numerics.Distributions.UnivariateDistributionBase]
            One fitted replicate per replication; null where the refit failed.
        """
        nd = _nd()
        import System  # noqa: PLC0415

        replicates = self._net_distributions()
        members = System.Array.CreateInstance(nd.UnivariateDistributionBase, len(replicates))
        for index, replicate in enumerate(replicates):
            if replicate is not None:
                members[index] = getattr(replicate, "__implementation__", replicate)
        return members

    def parameters(self) -> np.ndarray:
        """Run (once) and return the replicates' parameters (``Parameters``), NaN rows for failures.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(replications, number_of_parameters)``.
        """
        return to_numpy_matrix(self.net.Parameters(self._net_distributions()))

    def distributions(self) -> list[UnivariateDistribution | None]:
        """Run (once) and return the bootstrap replicate distributions (``Distributions``).

        Returns
        -------
        list of UnivariateDistribution or None
            One fitted distribution per replicate; None where the refit failed.
        """
        return [None if d is None else wrap(d) for d in self._net_distributions()]

    def quantiles(self, probabilities: Any) -> np.ndarray:
        """Return replicate quantiles (``Quantiles``).

        Parameters
        ----------
        probabilities : array_like
            Non-exceedance probabilities.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(replications, len(probabilities))``.
        """
        return to_numpy_matrix(self.net.Quantiles(to_net_doubles(probabilities), self._net_distributions()))

    def percentile_quantile_ci(self, probabilities: Any, alpha: float = 0.1) -> np.ndarray:
        """Return percentile confidence intervals of quantiles (``PercentileQuantileCI``).

        Parameters
        ----------
        probabilities : array_like
            Non-exceedance probabilities.
        alpha : float, optional
            Two-sided level; default 0.1.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), 2)``.
        """
        return to_numpy_matrix(self.net.PercentileQuantileCI(to_net_doubles(probabilities), float(alpha),
                                                             self._net_distributions()))

    def expected_probabilities(self, quantiles: Any) -> np.ndarray:
        """Return the mean replicate CDF at given values (``ExpectedProbabilities``).

        Parameters
        ----------
        quantiles : array_like
            Values at which the replicate CDFs are averaged.

        Returns
        -------
        numpy.ndarray
            Expected non-exceedance probabilities.
        """
        return to_numpy_doubles(self.net.ExpectedProbabilities(to_net_doubles(quantiles), self._net_distributions()))

    @property
    def distribution(self) -> UnivariateDistribution:
        """The parent distribution, wrapped (``Distribution``)."""
        return wrap(self.net.Distribution)

    @property
    def estimation_method(self) -> str:
        """The refit method name (``EstimationMethod``)."""
        return str(self.net.EstimationMethod)

    @property
    def sample_size(self) -> int:
        """The size of each resampled data set (``SampleSize``)."""
        return int(self.net.SampleSize)

    @property
    def replications(self) -> int:
        """The number of replicates (``Replications``)."""
        return int(self.net.Replications)

    @property
    def prng_seed(self) -> int:
        """The seed of the replicate seeds (``PRNGSeed``)."""
        return int(self.net.PRNGSeed)

    def parameter_sets(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Run (once) and return the replicate parameter sets as ``(values, fitness, weight)`` (``ParameterSets``).

        Returns
        -------
        values : numpy.ndarray
            Array of shape ``(replications, number_of_parameters)``; NaN rows for failed refits.
        fitness : numpy.ndarray
            One entry per replicate (NaN: the bootstrap records none).
        weight : numpy.ndarray
            One entry per replicate.
        """
        return _parameter_set_arrays(self.net.ParameterSets(self._net_distributions()))

    def product_moments(self) -> np.ndarray:
        """Return the product moments of every resampled data set (``ProductMoments``).

        Returns
        -------
        numpy.ndarray
            Array of shape ``(replications, 4)``: mean, standard deviation, skewness, kurtosis.
        """
        return to_numpy_matrix(self.net.ProductMoments())

    def linear_moments(self) -> np.ndarray:
        """Return the L-moments of every resampled data set (``LinearMoments``).

        Returns
        -------
        numpy.ndarray
            Array of shape ``(replications, 4)``: λ1, λ2, τ3, τ4.
        """
        return to_numpy_matrix(self.net.LinearMoments())

    def probabilities(self, quantiles: Any) -> np.ndarray:
        """Return every replicate CDF at given values (``Probabilities``).

        Parameters
        ----------
        quantiles : array_like
            Values at which the replicate CDFs are evaluated.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(replications, len(quantiles))``.
        """
        return to_numpy_matrix(self.net.Probabilities(to_net_doubles(quantiles), self._net_distributions()))

    def compute_min_max_quantiles(self, min_probability: float, max_probability: float) -> tuple[float, float]:
        """Return the smallest and largest replicate quantiles at two probabilities (``ComputeMinMaxQuantiles``).

        Parameters
        ----------
        min_probability : float
            The probability of the lower quantile.
        max_probability : float
            The probability of the upper quantile.

        Returns
        -------
        tuple of float
            The minimum over replicates of ``InverseCDF(min_probability)`` and the maximum of
            ``InverseCDF(max_probability)``.
        """
        values = self.net.ComputeMinMaxQuantiles(float(min_probability), float(max_probability),
                                                 self._net_distributions())
        return float(values[0]), float(values[1])

    def bias_corrected_quantile_ci(self, probabilities: Any, alpha: float = 0.1) -> np.ndarray:
        """Return bias-corrected percentile confidence intervals of quantiles (``BiasCorrectedQuantileCI``).

        Parameters
        ----------
        probabilities : array_like
            Non-exceedance probabilities.
        alpha : float, optional
            Two-sided level; default 0.1.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), 2)``.
        """
        return to_numpy_matrix(self.net.BiasCorrectedQuantileCI(to_net_doubles(probabilities), float(alpha),
                                                                self._net_distributions()))

    def normal_quantile_ci(self, probabilities: Any, alpha: float = 0.1) -> np.ndarray:
        """Return normal-approximation confidence intervals of quantiles (``NormalQuantileCI``).

        Parameters
        ----------
        probabilities : array_like
            Non-exceedance probabilities.
        alpha : float, optional
            Two-sided level; default 0.1.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), 2)``.
        """
        return to_numpy_matrix(self.net.NormalQuantileCI(to_net_doubles(probabilities), float(alpha),
                                                         self._net_distributions()))

    def bca_quantile_ci(self, sample: Any, probabilities: Any, alpha: float = 0.1) -> np.ndarray:
        """Return bias-corrected and accelerated confidence intervals of quantiles (``BCaQuantileCI``).

        Parameters
        ----------
        sample : array_like
            The observed data (for the jackknife acceleration).
        probabilities : array_like
            Non-exceedance probabilities.
        alpha : float, optional
            Two-sided level; default 0.1.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), 2)``.
        """
        return to_numpy_matrix(self.net.BCaQuantileCI(to_net_doubles(sample), to_net_doubles(probabilities),
                                                      float(alpha)))

    def bootstrap_t_quantile_ci(self, probabilities: Any, alpha: float = 0.1) -> np.ndarray:
        """Return bootstrap-t confidence intervals of quantiles (``BootstrapTQuantileCI``).

        Parameters
        ----------
        probabilities : array_like
            Non-exceedance probabilities.
        alpha : float, optional
            Two-sided level; default 0.1.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), 2)``.
        """
        return to_numpy_matrix(self.net.BootstrapTQuantileCI(to_net_doubles(probabilities), float(alpha)))

    def estimate(self, probabilities: Any, alpha: float = 0.1,
                 record_parameter_sets: bool = True) -> UncertaintyAnalysisResults:
        """Summarize the bootstrap ensemble (``Estimate``).

        Parameters
        ----------
        probabilities : array_like
            Non-exceedance probabilities.
        alpha : float, optional
            Two-sided level; default 0.1.
        record_parameter_sets : bool, optional
            Store the replicate parameter sets; default True.

        Returns
        -------
        UncertaintyAnalysisResults
            Mode curve, band and mean curve.
        """
        return UncertaintyAnalysisResults(self.net.Estimate(to_net_doubles(probabilities), float(alpha),
                                                            self._net_distributions(), bool(record_parameter_sets)))


__all__ = ["BootstrapAnalysis", "UncertaintyAnalysisResults"]

translate_public(globals())
