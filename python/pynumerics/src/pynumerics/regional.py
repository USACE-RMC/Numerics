"""Hosking-Wallis regional frequency analysis by L-moments from ``Numerics.Distributions.Regional``.

:class:`RegionalData` holds a region's sites (name, record length, mean and L-moment ratios), built from raw samples
or from precomputed L-moments, with the record-length weighted regional average (Hosking and Wallis, 1997, equation
4.3) and the discordancy measure ``D_i`` (section 3.2.3). :class:`RegionalHomogeneityTest` runs the Monte Carlo
simulations of a homogeneous kappa region and gives the heterogeneity measures ``H1``, ``H2`` and ``H3`` (section
4.3.3) and the goodness-of-fit measures ``Z`` of GLO, GEV, GNO, PE3 and GPA (section 5.2.3), as Hosking's ``REGTST``
and R's ``lmomRFA::regtst`` do. :class:`RegionalGrowthCurve` fits a distribution to the regional L-moments with mean 1
and gives growth factors and site quantiles (index flood times growth factor).

The simulations are seeded per simulated region from the master seed, so a seed gives the same results whatever the
number of threads. The sites are simulated as independent.

Distributions are named by their ``UnivariateDistributionType`` member (``"GeneralizedExtremeValue"``) or by
Hosking's abbreviation (``"GEV"``), case-insensitively.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics._util import _nd, namespace, net_enum
from pynumerics.distributions import UnivariateDistribution, wrap
from pynumerics.errors import translate_public
from pynumerics.marshal import to_net_doubles, to_net_ints, to_numpy_doubles, to_numpy_matrix

ABBREVIATIONS: dict[str, str] = {
    "GLO": "GeneralizedLogistic",
    "GEV": "GeneralizedExtremeValue",
    "GNO": "GeneralizedNormal",
    "PE3": "PearsonTypeIII",
    "GPA": "GeneralizedPareto",
    "KAP": "KappaFour",
}
"""Hosking's three-letter abbreviations, mapped to the ``UnivariateDistributionType`` member names."""


def _regional() -> Any:
    """Return the ``Numerics.Distributions.Regional`` namespace after starting the runtime."""
    return namespace("Numerics.Distributions.Regional")


def _type(distribution: Any) -> Any:
    """Return the .NET ``UnivariateDistributionType`` member for a name or abbreviation."""
    if isinstance(distribution, str):
        distribution = ABBREVIATIONS.get(distribution.upper(), distribution)
    return net_enum(_nd().UnivariateDistributionType, distribution)


def _strings(values: Any) -> Any:
    """Copy names into a new ``System.String[]``."""
    import System  # noqa: PLC0415 - the runtime is running

    items = [str(v) for v in values]
    array = System.Array.CreateInstance(System.String, len(items))
    for index, item in enumerate(items):
        array[index] = item
    return array


class RegionalData:
    """The sites of a region: names, record lengths, means and L-moment ratios (``RegionalData``).

    Parameters
    ----------
    names : sequence of str
        The site names.
    record_lengths : array_like of int
        The record length of each site, at least 1.
    mean : array_like
        The sample mean ``l1`` of each site (its index flood).
    l_cv : array_like
        The sample L-CV ``t = l2 / l1`` of each site.
    tau3 : array_like
        The sample L-skewness of each site.
    tau4 : array_like
        The sample L-kurtosis of each site.
    tau5 : array_like, optional
        The sample ``t5`` of each site; None leaves every ``t5`` NaN.

    Raises
    ------
    ValueError
        If the lists differ in length, there are no sites, a record length is below 1, a ratio is not finite, an
        L-CV is not positive or a ``t5`` is infinite.
    """

    def __init__(self, names: Any, record_lengths: Any, mean: Any, l_cv: Any, tau3: Any, tau4: Any,
                 tau5: Any = None) -> None:
        self._net = _regional().RegionalData(
            _strings(names), to_net_ints(record_lengths), to_net_doubles(mean), to_net_doubles(l_cv),
            to_net_doubles(tau3), to_net_doubles(tau4), None if tau5 is None else to_net_doubles(tau5))

    @classmethod
    def from_net(cls, net: Any) -> RegionalData:
        """Wrap an existing .NET ``RegionalData``.

        Parameters
        ----------
        net : Numerics.Distributions.Regional.RegionalData
            The .NET object.

        Returns
        -------
        RegionalData
            The wrapper.
        """
        instance = cls.__new__(cls)
        instance._net = net
        return instance

    @classmethod
    def from_samples(cls, samples: Any, names: Any = None) -> RegionalData:
        """Build regional data from each site's raw sample (``FromSamples``).

        Parameters
        ----------
        samples : sequence of array_like, or 2-D array_like
            The sample of each site; a 2-D array has one site per row. NaN marks a missing year and is left out.
        names : sequence of str, optional
            The site names; None names the sites ``"1"``, ``"2"``, and so on.

        Returns
        -------
        RegionalData
            The sites, with the unbiased sample L-moments of ``Statistics.LinearMoments`` up to ``t5``.

        Raises
        ------
        ValueError
            If there are no samples, or a sample holds an infinity, has fewer than four values, a nonpositive mean or
            no spread.
        """
        regional = _regional()  # starts the runtime, which the System imports below need
        from System import Double  # noqa: PLC0415
        from System.Collections.Generic import IList, List  # noqa: PLC0415

        rows = List[IList[Double]]()
        for sample in samples:
            rows.Add(to_net_doubles(np.asarray(sample, dtype=np.float64).ravel()))
        return cls.from_net(regional.RegionalData.FromSamples(rows, None if names is None else _strings(names)))

    @property
    def net(self) -> Any:
        """The wrapped .NET ``RegionalData``."""
        return self._net

    @property
    def site_count(self) -> int:
        """The number of sites (``SiteCount``)."""
        return int(self._net.SiteCount)

    @property
    def names(self) -> list[str]:
        """The site names (``Names``)."""
        return [str(name) for name in self._net.Names]

    @property
    def record_lengths(self) -> np.ndarray:
        """The record length of each site, int32 (``RecordLengths``)."""
        return np.asarray([int(n) for n in self._net.RecordLengths], dtype=np.int32)

    @property
    def mean(self) -> np.ndarray:
        """The sample mean ``l1`` of each site (``Mean``)."""
        return to_numpy_doubles(self._net.Mean)

    @property
    def l_cv(self) -> np.ndarray:
        """The sample L-CV of each site (``LCV``)."""
        return to_numpy_doubles(self._net.LCV)

    @property
    def tau3(self) -> np.ndarray:
        """The sample L-skewness of each site (``Tau3``)."""
        return to_numpy_doubles(self._net.Tau3)

    @property
    def tau4(self) -> np.ndarray:
        """The sample L-kurtosis of each site (``Tau4``)."""
        return to_numpy_doubles(self._net.Tau4)

    @property
    def tau5(self) -> np.ndarray:
        """The sample ``t5`` of each site, NaN where unavailable (``Tau5``)."""
        return to_numpy_doubles(self._net.Tau5)

    def site_linear_moments(self, site: int) -> np.ndarray:
        """Return one site's L-moments (``SiteLinearMoments``).

        Parameters
        ----------
        site : int
            The zero-based site index.

        Returns
        -------
        numpy.ndarray
            ``(l1, l2, t3, t4, t5)``.
        """
        return to_numpy_doubles(self._net.SiteLinearMoments(int(site)))

    def regional_average(self, weights: Any = None) -> np.ndarray:
        """Return the regional average L-moments (``RegionalAverage``).

        Parameters
        ----------
        weights : array_like, optional
            One nonnegative weight per site; None weights by record length (equation 4.3).

        Returns
        -------
        numpy.ndarray
            ``(1, t^R, t3^R, t4^R, t5^R)``, the regional L-moments scaled to ``l1 = 1``.
        """
        if weights is None:
            return to_numpy_doubles(self._net.RegionalAverage())
        return to_numpy_doubles(self._net.RegionalAverage(to_net_doubles(weights)))

    def discordancy(self) -> np.ndarray:
        """Return the discordancy measure of each site (``Discordancy``).

        Returns
        -------
        numpy.ndarray
            ``D_i`` by site; all 1 for three or fewer sites, all NaN when the cross-product matrix is singular.
        """
        return to_numpy_doubles(self._net.Discordancy())

    @staticmethod
    def discordancy_critical_value(site_count: int, significance: float = 0.10) -> float:
        """Return the critical discordancy for a region (``DiscordancyCriticalValue``).

        Parameters
        ----------
        site_count : int
            The number of sites, at least 1.
        significance : float, optional
            0.10 (Hosking and Wallis's Table 3.1, 3 from 15 sites) or 0.05 (4 from 19 sites).

        Returns
        -------
        float
            The critical value.
        """
        return float(_regional().RegionalData.DiscordancyCriticalValue(int(site_count), float(significance)))


class RegionalHomogeneityTest:
    """Heterogeneity ``H`` and goodness-of-fit ``Z`` of a region by kappa simulation (``RegionalHomogeneityTest``).

    Parameters
    ----------
    data : RegionalData
        The region, at least two sites with record lengths of at least 4.
    simulations : int, optional
        The number of simulated regions, at least 2 (Hosking and Wallis use 500).
    seed : int, optional
        The master seed; each simulated region takes its own seed drawn from it.
    max_degree_of_parallelism : int, optional
        The largest number of simulations run at once; -1 lets the scheduler choose. The results do not depend
        on it.

    Raises
    ------
    ValueError
        If the region has fewer than two sites or a record shorter than 4, or ``simulations`` is below 2.
    """

    def __init__(self, data: RegionalData, simulations: int = 500, seed: int = 12345,
                 max_degree_of_parallelism: int = -1) -> None:
        self._data = data
        self._net = _regional().RegionalHomogeneityTest(data.net, int(simulations), int(seed),
                                                        int(max_degree_of_parallelism))

    @property
    def net(self) -> Any:
        """The wrapped .NET ``RegionalHomogeneityTest``."""
        return self._net

    @property
    def data(self) -> RegionalData:
        """The regional data tested (``Data``)."""
        return self._data

    @property
    def simulations(self) -> int:
        """The number of simulated regions (``Simulations``)."""
        return int(self._net.Simulations)

    @property
    def seed(self) -> int:
        """The master seed (``Seed``)."""
        return int(self._net.Seed)

    @property
    def regional_linear_moments(self) -> np.ndarray:
        """The regional average ``(1, t^R, t3^R, t4^R, t5^R)`` (``RegionalLinearMoments``)."""
        return to_numpy_doubles(self._net.RegionalLinearMoments)

    @property
    def discordancy(self) -> np.ndarray:
        """The discordancy ``D_i`` of each site (``Discordancy``)."""
        return to_numpy_doubles(self._net.Discordancy)

    @property
    def discordancy_critical_value(self) -> float:
        """The 10% critical discordancy for the number of sites (``DiscordancyCriticalValue``)."""
        return float(self._net.DiscordancyCriticalValue)

    @property
    def kappa_parameters(self) -> np.ndarray:
        """The simulated kappa parent ``(xi, alpha, k, h)``; ``h = -1`` for the GLO stand-in (``KappaParameters``)."""
        return to_numpy_doubles(self._net.KappaParameters)

    @property
    def is_generalized_logistic_parent(self) -> bool:
        """Whether the kappa fit failed and the GLO is the parent (``IsGeneralizedLogisticParent``)."""
        return bool(self._net.IsGeneralizedLogisticParent)

    @property
    def observed_dispersion(self) -> np.ndarray:
        """The observed ``(V1, V2, V3)`` (``ObservedDispersion``)."""
        return to_numpy_doubles(self._net.ObservedDispersion)

    @property
    def simulated_dispersion_mean(self) -> np.ndarray:
        """The mean of each V over the simulations (``SimulatedDispersionMean``)."""
        return to_numpy_doubles(self._net.SimulatedDispersionMean)

    @property
    def simulated_dispersion_standard_deviation(self) -> np.ndarray:
        """The sample standard deviation of each V (``SimulatedDispersionStandardDeviation``)."""
        return to_numpy_doubles(self._net.SimulatedDispersionStandardDeviation)

    @property
    def heterogeneity(self) -> np.ndarray:
        """The heterogeneity measures ``(H1, H2, H3)`` (``Heterogeneity``)."""
        return to_numpy_doubles(self._net.Heterogeneity)

    @property
    def nonpositive_mean_simulations(self) -> int:
        """Simulated regions with a site of zero or negative mean, which corrupts H (``NonpositiveMeanSimulations``)."""
        return int(self._net.NonpositiveMeanSimulations)

    @property
    def tau4_bias(self) -> float:
        """The bias ``B4`` of the regional L-kurtosis (``Tau4Bias``)."""
        return float(self._net.Tau4Bias)

    @property
    def tau4_standard_deviation(self) -> float:
        """``sigma4``, the standard deviation of the simulated regional L-kurtosis (``Tau4StandardDeviation``)."""
        return float(self._net.Tau4StandardDeviation)

    @property
    def distributions(self) -> list[str]:
        """The candidate distributions, GLO, GEV, GNO, PE3 and GPA, by type name (``Distributions``)."""
        return [str(d) for d in self._net.Distributions]

    @property
    def distribution_tau4(self) -> np.ndarray:
        """The L-kurtosis of each fitted candidate (``DistributionTau4``)."""
        return to_numpy_doubles(self._net.DistributionTau4)

    @property
    def z(self) -> dict[str, float]:
        """``Z^DIST`` of each candidate, keyed by Hosking's abbreviation (``Z``)."""
        values = to_numpy_doubles(self._net.Z)
        return dict(zip(("GLO", "GEV", "GNO", "PE3", "GPA"), (float(v) for v in values), strict=True))

    def get_z(self, distribution: str) -> float:
        """Return ``Z^DIST`` of one candidate (``GetZ``).

        Parameters
        ----------
        distribution : str
            A candidate, by abbreviation or type name.

        Returns
        -------
        float
            The goodness-of-fit measure.
        """
        return float(self._net.GetZ(_type(distribution)))

    def acceptable_distributions(self, critical_value: float = 1.64) -> list[str]:
        """Return the candidates with ``|Z|`` at most the critical value, best first (``AcceptableDistributions``).

        Parameters
        ----------
        critical_value : float, optional
            The largest acceptable ``|Z|``; 1.64 is the 10% level.

        Returns
        -------
        list of str
            The accepted candidates' type names.
        """
        return [str(d) for d in self._net.AcceptableDistributions(float(critical_value))]


class RegionalGrowthCurve:
    """The regional growth curve: a distribution with mean 1 fitted to regional L-moments (``RegionalGrowthCurve``).

    Parameters
    ----------
    distribution : str
        ``"GLO"``, ``"GEV"``, ``"GNO"``, ``"PE3"``, ``"GPA"`` or ``"KAP"``, or the type name.
    data : RegionalData
        The region; the curve is fitted to its record-length weighted regional average.

    Raises
    ------
    ValueError
        If the distribution is not supported or the fit fails.
    """

    def __init__(self, distribution: str, data: RegionalData) -> None:
        self._data: RegionalData | None = data
        self._net = _regional().RegionalGrowthCurve(_type(distribution), data.net)

    @classmethod
    def from_ratios(cls, distribution: str, l_cv: float, tau3: float,
                    tau4: float = float("nan")) -> RegionalGrowthCurve:
        """Fit the growth curve to regional L-moment ratios.

        Parameters
        ----------
        distribution : str
            The family, as for the constructor.
        l_cv : float
            The regional L-CV, positive.
        tau3 : float
            The regional L-skewness.
        tau4 : float, optional
            The regional L-kurtosis; required only by the kappa.

        Returns
        -------
        RegionalGrowthCurve
            The fitted curve.
        """
        instance = cls.__new__(cls)
        instance._data = None
        instance._net = _regional().RegionalGrowthCurve(_type(distribution), float(l_cv), float(tau3), float(tau4))
        return instance

    @property
    def net(self) -> Any:
        """The wrapped .NET ``RegionalGrowthCurve``."""
        return self._net

    @property
    def data(self) -> RegionalData | None:
        """The regional data the curve was fitted to, or None for ratios (``Data``)."""
        return self._data

    @property
    def distribution_type(self) -> str:
        """The fitted family's type name (``DistributionType``)."""
        return str(self._net.DistributionType)

    @property
    def parameters(self) -> np.ndarray:
        """The growth curve parameters in the family's order (``Parameters``)."""
        return to_numpy_doubles(self._net.Parameters)

    def get_distribution(self) -> UnivariateDistribution:
        """Return a copy of the growth curve distribution, whose mean is 1 (``GetDistribution``).

        Returns
        -------
        UnivariateDistribution
            The typed distribution wrapper.
        """
        return wrap(self._net.GetDistribution())

    def growth_factors(self, probabilities: Any) -> Any:
        """Return the growth factors at non-exceedance probabilities (``GrowthFactor``, ``GrowthFactors``).

        Parameters
        ----------
        probabilities : float or array_like
            Non-exceedance probabilities in [0, 1].

        Returns
        -------
        float or numpy.ndarray
            A float for a scalar, otherwise an array shaped like ``probabilities``.
        """
        p = np.asarray(probabilities, dtype=np.float64)
        if p.ndim == 0:
            return float(self._net.GrowthFactor(float(p)))
        return to_numpy_doubles(self._net.GrowthFactors(to_net_doubles(p.ravel()))).reshape(p.shape)

    def site_quantiles(self, index_flood: Any, probabilities: Any) -> np.ndarray:
        """Return site quantiles, index flood times growth factor (``SiteQuantiles``).

        Parameters
        ----------
        index_flood : float or RegionalData
            One site's index flood, or a region whose site means are the index floods.
        probabilities : array_like
            Non-exceedance probabilities in [0, 1].

        Returns
        -------
        numpy.ndarray
            One quantile per probability, or a ``(sites, probabilities)`` array for a region.
        """
        p = to_net_doubles(np.asarray(probabilities, dtype=np.float64).ravel())
        if isinstance(index_flood, RegionalData):
            return to_numpy_matrix(self._net.SiteQuantiles(index_flood.net, p))
        return to_numpy_doubles(self._net.SiteQuantiles(float(index_flood), p))


__all__ = ["ABBREVIATIONS", "RegionalData", "RegionalGrowthCurve", "RegionalHomogeneityTest"]

translate_public(globals())
