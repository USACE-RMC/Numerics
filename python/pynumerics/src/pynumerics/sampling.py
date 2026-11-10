"""Sampling designs and resampling from ``Numerics.Sampling``.

Stratification, Latin hypercube and Sobol designs, the Mersenne Twister, and the general
bootstrap (:class:`Bootstrap`, ``Numerics.Sampling.Bootstrap<TData>``) with its confidence
intervals. The MCMC samplers are in :mod:`pynumerics.mcmc`.

Stratification bins come back as a :class:`StratificationBins` record of NumPy arrays
holding the .NET bin list (``.net``), which can be passed back to Numerics members that
take ``List<StratificationBin>``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from pynumerics import callbacks
from pynumerics._util import namespace, net_enum
from pynumerics.enums import ImportanceDistribution, NetEnum
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _numerics_matrix,
    _parameter_set_arrays,
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_ints,
    to_numpy_matrix,
)


def _sampling() -> Any:
    """Return the ``Numerics.Sampling`` namespace after starting the runtime."""
    return namespace("Numerics.Sampling")


def _xelement(text: str) -> Any:
    """Parse XML text into a ``System.Xml.Linq.XElement``."""
    _sampling()
    from System.Xml.Linq import XElement  # noqa: PLC0415 - the runtime is running

    return XElement.Parse(str(text))


@dataclass(frozen=True)
class StratificationBin:
    """One stratification bin as Numerics holds it (``StratificationBin``).

    Parameters
    ----------
    lower_bound : float
        Lower edge.
    upper_bound : float
        Upper edge.
    weight : float
        Weight (probability mass); Numerics sets the bin width for value bins built
        without one.
    midpoint : float
        ``Midpoint`` as Numerics computes it, (lower + upper) / 2.
    """

    lower_bound: float
    upper_bound: float
    weight: float
    midpoint: float

    def contains(self, x: float) -> bool:
        """Return whether ``lower_bound <= x < upper_bound``.

        Parameters
        ----------
        x : float
            The value.

        Returns
        -------
        bool
            True when the value lies in the half-open bin.
        """
        return self.lower_bound <= x < self.upper_bound

    def to_net(self) -> Any:
        """Return the bin as a new .NET ``StratificationBin``.

        Returns
        -------
        Numerics.Sampling.StratificationBin
            A bin with the same bounds and weight.
        """
        return _sampling().StratificationBin(float(self.lower_bound), float(self.upper_bound), float(self.weight))

    def compare_to(self, other: StratificationBin) -> int:
        """Order two non-overlapping bins (``CompareTo``).

        Parameters
        ----------
        other : StratificationBin
            The other bin.

        Returns
        -------
        int
            1 when this bin lies above ``other``, -1 when below.

        Raises
        ------
        ValueError
            If the bins overlap; Numerics counts a bin as overlapping itself.
        """
        return int(self.to_net().CompareTo(other.to_net()))


@dataclass(frozen=True)
class StratificationBins:
    """Stratification bins as arrays.

    Parameters
    ----------
    lower : numpy.ndarray
        Lower bound of each bin.
    upper : numpy.ndarray
        Upper bound of each bin.
    weight : numpy.ndarray
        Weight (probability mass) of each bin.
    net : System.Collections.Generic.List[Numerics.Sampling.StratificationBin]
        The .NET bin list the arrays were read from.
    """

    lower: np.ndarray
    upper: np.ndarray
    weight: np.ndarray
    net: Any = field(repr=False, compare=False)

    @property
    def midpoint(self) -> np.ndarray:
        """Bin midpoints, (lower + upper) / 2 as Numerics computes them."""
        return np.array([float(b.Midpoint) for b in self.net], dtype=np.float64)

    def __len__(self) -> int:
        return int(self.lower.size)

    @property
    def bins(self) -> tuple[StratificationBin, ...]:
        """The bins as :class:`StratificationBin` records (``LowerBound``, ``UpperBound``, ``Weight``, ``Midpoint``)."""
        return tuple(StratificationBin(float(b.LowerBound), float(b.UpperBound), float(b.Weight), float(b.Midpoint))
                     for b in self.net)

    @classmethod
    def from_net(cls, bins: Any) -> StratificationBins:
        """Read a .NET ``List<StratificationBin>``.

        Parameters
        ----------
        bins : System.Collections.Generic.List[Numerics.Sampling.StratificationBin]
            The .NET bins.

        Returns
        -------
        StratificationBins
            The arrays and the original list.
        """
        lower = np.array([float(b.LowerBound) for b in bins], dtype=np.float64)
        upper = np.array([float(b.UpperBound) for b in bins], dtype=np.float64)
        weight = np.array([float(b.Weight) for b in bins], dtype=np.float64)
        return cls(lower, upper, weight, bins)

    def to_xml(self) -> str:
        """Serialize the bins (``StratificationBin.SaveToXElement(IList<StratificationBin>)``).

        Returns
        -------
        str
            The ``StratificationBinList`` XML element, bounds and weights in G17.
        """
        return str(_sampling().StratificationBin.SaveToXElement(self.net).ToString())

    @classmethod
    def from_xml(cls, text: str) -> StratificationBins:
        """Read bins written by :meth:`to_xml` (``StratificationBin.XElementToStratificationBinList``).

        Parameters
        ----------
        text : str
            A ``StratificationBinList`` XML element.

        Returns
        -------
        StratificationBins
            The bins.
        """
        return cls.from_net(_sampling().StratificationBin.XElementToStratificationBinList(_xelement(text)))


class StratificationOptions:
    """The range and bin count of one stratification (``Numerics.Sampling.StratificationOptions``).

    Parameters
    ----------
    lower_bound : float
        Lower end of the range.
    upper_bound : float
        Upper end of the range.
    number_of_bins : int
        Number of bins.
    is_probability : bool, optional
        Whether the range is a probability range; default False.

    Attributes
    ----------
    net : Numerics.Sampling.StratificationOptions
        The wrapped .NET object.
    """

    def __init__(self, lower_bound: float, upper_bound: float, number_of_bins: int,
                 is_probability: bool = False) -> None:
        self.net = _sampling().StratificationOptions(float(lower_bound), float(upper_bound), int(number_of_bins),
                                                     bool(is_probability))

    @classmethod
    def _wrap(cls, net: Any) -> StratificationOptions:
        instance = cls.__new__(cls)
        instance.net = net
        return instance

    def __repr__(self) -> str:
        return (f"StratificationOptions({self.lower_bound!r}, {self.upper_bound!r}, {self.number_of_bins!r}, "
                f"is_probability={self.is_probability!r})")

    @property
    def lower_bound(self) -> float:
        """Lower end of the range (``LowerBound``)."""
        return float(self.net.LowerBound)

    @property
    def upper_bound(self) -> float:
        """Upper end of the range (``UpperBound``)."""
        return float(self.net.UpperBound)

    @property
    def number_of_bins(self) -> int:
        """Number of bins (``NumberOfBins``)."""
        return int(self.net.NumberOfBins)

    @property
    def is_probability(self) -> bool:
        """Whether the range is a probability range (``IsProbability``)."""
        return bool(self.net.IsProbability)

    @property
    def is_valid(self) -> bool:
        """Whether the options are valid (``IsValid``)."""
        return bool(self.net.IsValid)

    def get_errors(self) -> list[str]:
        """Return the validation messages (``GetErrors``).

        Returns
        -------
        list of str
            Empty when the options are valid.
        """
        return [str(e) for e in self.net.GetErrors()]

    def clone(self) -> StratificationOptions:
        """Return a copy (``Clone``).

        Returns
        -------
        StratificationOptions
            The copy.
        """
        return StratificationOptions._wrap(self.net.Clone())

    def to_xml(self) -> str:
        """Serialize the options (``SaveToXElement``).

        Returns
        -------
        str
            The ``StratificationOptions`` XML element.
        """
        return str(self.net.SaveToXElement().ToString())

    @classmethod
    def from_xml(cls, text: str) -> StratificationOptions:
        """Read options written by :meth:`to_xml` (``StratificationOptions(XElement)``).

        Parameters
        ----------
        text : str
            A ``StratificationOptions`` XML element.

        Returns
        -------
        StratificationOptions
            The options.
        """
        return cls._wrap(_sampling().StratificationOptions(_xelement(text)))

    @staticmethod
    def list_to_xml(options: Sequence[StratificationOptions]) -> str:
        """Serialize several options (``StratificationOptions.SaveToXElement(IList<StratificationOptions>)``).

        Parameters
        ----------
        options : sequence of StratificationOptions
            The options.

        Returns
        -------
        str
            The XML element holding every option.
        """
        return str(_sampling().StratificationOptions.SaveToXElement(_options_list(options)).ToString())

    @staticmethod
    def list_from_xml(text: str) -> list[StratificationOptions]:
        """Read options written by :meth:`list_to_xml` (``XElementToStratificationOptions``).

        Parameters
        ----------
        text : str
            The XML element.

        Returns
        -------
        list of StratificationOptions
            The options.
        """
        items = _sampling().StratificationOptions.XElementToStratificationOptions(_xelement(text))
        return [StratificationOptions._wrap(o) for o in items]


def _options(lower: float, upper: float, number_of_bins: int, is_probability: bool) -> Any:
    """Build ``StratificationOptions``."""
    return _sampling().StratificationOptions(float(lower), float(upper), int(number_of_bins), bool(is_probability))


def _options_list(options: Any, is_probability: bool = False) -> Any:
    """Build ``List<StratificationOptions>`` from :class:`StratificationOptions` or ``(lower, upper, bins)`` triples."""
    s = _sampling()
    from System.Collections.Generic import List  # noqa: PLC0415 - the runtime is started by _sampling

    items = List[s.StratificationOptions]()
    for option in options:
        if isinstance(option, StratificationOptions):
            items.Add(option.net)
        else:
            lower, upper, number_of_bins = option
            items.Add(_options(lower, upper, number_of_bins, is_probability))
    return items


def stratify_x_values(lower: float, upper: float, number_of_bins: int,
                      is_logarithmic: bool = False) -> StratificationBins:
    """Stratify a value range into equal-width bins (``Stratify.XValues``).

    Parameters
    ----------
    lower, upper : float
        Range to stratify.
    number_of_bins : int
        Number of bins.
    is_logarithmic : bool, optional
        Equal widths in log10 space; default False.

    Returns
    -------
    StratificationBins
        The bins; Numerics leaves their weights unset (−1) for value bins.
    """
    s = _sampling()
    return StratificationBins.from_net(s.Stratify.XValues(_options(lower, upper, number_of_bins, False),
                                                          bool(is_logarithmic)))


def stratify_x_values_many(options: Any, is_logarithmic: bool = False) -> StratificationBins:
    """Stratify several value ranges into one bin list (``Stratify.XValues(List<StratificationOptions>)``).

    Parameters
    ----------
    options : sequence of (float, float, int) or StratificationOptions
        One ``(lower, upper, number_of_bins)`` triple or options object per range.
    is_logarithmic : bool, optional
        Equal widths in log10 space; default False.

    Returns
    -------
    StratificationBins
        The bins of every range, in the order Numerics returns them.
    """
    s = _sampling()
    return StratificationBins.from_net(s.Stratify.XValues(_options_list(options), bool(is_logarithmic)))


def stratify_probabilities(lower: float, upper: float, number_of_bins: int,
                           importance: ImportanceDistribution | str = ImportanceDistribution.UNIFORM,
                           is_exhaustive: bool = True) -> StratificationBins:
    """Stratify a non-exceedance probability range (``Stratify.Probabilities``).

    Parameters
    ----------
    lower, upper : float
        Probability range in [0, 1].
    number_of_bins : int
        Number of bins.
    importance : ImportanceDistribution or str, optional
        Space in which the bins are equally wide; default uniform.
    is_exhaustive : bool, optional
        Whether the end bins absorb the remaining probability; default True.

    Returns
    -------
    StratificationBins
        The bins.
    """
    s = _sampling()
    kind = net_enum(s.Stratify.ImportanceDistribution, importance)
    return StratificationBins.from_net(s.Stratify.Probabilities(_options(lower, upper, number_of_bins, True), kind,
                                                                bool(is_exhaustive)))


def stratify_probabilities_many(options: Any, importance: ImportanceDistribution | str = ImportanceDistribution.UNIFORM,
                                is_exhaustive: bool = True) -> StratificationBins:
    """Stratify several probability ranges into one bin list (``Stratify.Probabilities(List<...>)``).

    Parameters
    ----------
    options : sequence of (float, float, int) or StratificationOptions
        One ``(lower, upper, number_of_bins)`` probability range or options object per range.
    importance : ImportanceDistribution or str, optional
        Space in which the bins are equally wide; default uniform.
    is_exhaustive : bool, optional
        Whether the end bins absorb the remaining probability; default True.

    Returns
    -------
    StratificationBins
        The bins of every range.
    """
    s = _sampling()
    kind = net_enum(s.Stratify.ImportanceDistribution, importance)
    return StratificationBins.from_net(s.Stratify.Probabilities(_options_list(options, True), kind,
                                                                bool(is_exhaustive)))


def stratify_exceedance_probabilities(lower: float, upper: float, number_of_bins: int,
                                      importance: ImportanceDistribution | str = ImportanceDistribution.UNIFORM,
                                      is_exhaustive: bool = True) -> StratificationBins:
    """Stratify an exceedance probability range (``Stratify.ExceedanceProbabilities``).

    Parameters
    ----------
    lower, upper : float
        Exceedance probability range in [0, 1].
    number_of_bins : int
        Number of bins.
    importance : ImportanceDistribution or str, optional
        Space in which the bins are equally wide; default uniform.
    is_exhaustive : bool, optional
        Whether the end bins absorb the remaining probability; default True.

    Returns
    -------
    StratificationBins
        The bins.
    """
    s = _sampling()
    kind = net_enum(s.Stratify.ImportanceDistribution, importance)
    return StratificationBins.from_net(s.Stratify.ExceedanceProbabilities(
        _options(lower, upper, number_of_bins, True), kind, bool(is_exhaustive)))


def stratify_exceedance_probabilities_many(options: Any,
                                           importance: ImportanceDistribution | str = ImportanceDistribution.UNIFORM,
                                           is_exhaustive: bool = True) -> StratificationBins:
    """Stratify several exceedance probability ranges (``Stratify.ExceedanceProbabilities(List<...>)``).

    Parameters
    ----------
    options : sequence of (float, float, int) or StratificationOptions
        One ``(lower, upper, number_of_bins)`` range or options object per range.
    importance : ImportanceDistribution or str, optional
        Space in which the bins are equally wide; default uniform.
    is_exhaustive : bool, optional
        Whether the end bins absorb the remaining probability; default True.

    Returns
    -------
    StratificationBins
        The bins of every range.
    """
    s = _sampling()
    kind = net_enum(s.Stratify.ImportanceDistribution, importance)
    return StratificationBins.from_net(s.Stratify.ExceedanceProbabilities(_options_list(options, True), kind,
                                                                          bool(is_exhaustive)))


def stratify_multivariate_probabilities(lower: float, upper: float, number_of_bins: int,
                                        importance: ImportanceDistribution | str = ImportanceDistribution.UNIFORM,
                                        is_exhaustive: bool = True, dimension: int = 1, seed: int = -1,
                                        correlation: Any = None) -> list[StratificationBins]:
    """Stratify a probability range in several dimensions (``Stratify.MultivariateProbabilities``).

    Parameters
    ----------
    lower, upper : float
        Probability range in [0, 1].
    number_of_bins : int
        Number of bins per dimension.
    importance : ImportanceDistribution or str, optional
        Space in which the bins are equally wide; default uniform.
    is_exhaustive : bool, optional
        Whether the end bins absorb the remaining probability; default True.
    dimension : int, optional
        Number of dimensions; default 1.
    seed : int, optional
        Seed of the random pairing; zero or negative seeds from the clock (not reproducible).
    correlation : array_like, optional
        Correlation matrix of shape ``(dimension, dimension)``; default independence.

    Returns
    -------
    list of StratificationBins
        The bins of each dimension.
    """
    s = _sampling()
    kind = net_enum(s.Stratify.ImportanceDistribution, importance)
    matrix = None if correlation is None else to_net_matrix(correlation)
    result = s.Stratify.MultivariateProbabilities(_options(lower, upper, number_of_bins, True), kind,
                                                  bool(is_exhaustive), int(dimension), int(seed), matrix)
    return [StratificationBins.from_net(bins) for bins in result]


def _bins(bins: Any) -> Any:
    """Return the .NET bin list of :class:`StratificationBins` or pass a .NET list through."""
    return getattr(bins, "net", bins)


def x_to_probability(x_bins: Any, transform: Any, is_exhaustive: bool = True) -> StratificationBins:
    """Map value bins to non-exceedance probability bins (``Stratify.XToProbability``).

    Parameters
    ----------
    x_bins : StratificationBins
        Value bins.
    transform : callable, str or System.Delegate
        The CDF ``x -> p`` (``Func<double, double>``), for example ``Normal(0, 1).cdf``.
    is_exhaustive : bool, optional
        Whether the end bins absorb the remaining probability; default True.

    Returns
    -------
    StratificationBins
        The probability bins.
    """
    s = _sampling()
    return StratificationBins.from_net(s.Stratify.XToProbability(_bins(x_bins), callbacks.scalar_function(transform),
                                                                 bool(is_exhaustive)))


def x_to_exceedance_probability(x_bins: Any, transform: Any, is_exhaustive: bool = True) -> StratificationBins:
    """Map value bins to exceedance probability bins (``Stratify.XToExceedanceProbability``).

    Parameters
    ----------
    x_bins : StratificationBins
        Value bins.
    transform : callable, str or System.Delegate
        The survival function ``x -> 1 - F(x)`` (``Func<double, double>``).
    is_exhaustive : bool, optional
        Whether the end bins absorb the remaining probability; default True.

    Returns
    -------
    StratificationBins
        The exceedance probability bins.
    """
    s = _sampling()
    return StratificationBins.from_net(s.Stratify.XToExceedanceProbability(
        _bins(x_bins), callbacks.scalar_function(transform), bool(is_exhaustive)))


def probability_to_x(probability_bins: Any, transform: Any) -> StratificationBins:
    """Map non-exceedance probability bins to value bins (``Stratify.ProbabilityToX``).

    Parameters
    ----------
    probability_bins : StratificationBins
        Probability bins.
    transform : callable, str or System.Delegate
        The quantile function ``p -> x`` (``Func<double, double>``).

    Returns
    -------
    StratificationBins
        The value bins; Numerics sets each weight to the bin width, as for value bins.
    """
    s = _sampling()
    return StratificationBins.from_net(s.Stratify.ProbabilityToX(_bins(probability_bins),
                                                                 callbacks.scalar_function(transform)))


def exceedance_probability_to_x(probability_bins: Any, transform: Any) -> StratificationBins:
    """Map exceedance probability bins to value bins (``Stratify.ExceedanceProbabilityToX``).

    Parameters
    ----------
    probability_bins : StratificationBins
        Exceedance probability bins.
    transform : callable, str or System.Delegate
        The function ``q -> x`` (``Func<double, double>``), for example ``q -> F^-1(1 - q)``.

    Returns
    -------
    StratificationBins
        The value bins, ascending in x; Numerics sets each weight to the bin width.
    """
    s = _sampling()
    return StratificationBins.from_net(s.Stratify.ExceedanceProbabilityToX(_bins(probability_bins),
                                                                           callbacks.scalar_function(transform)))


def latin_hypercube(sample_size: int, dimension: int = 1, seed: Any = -1, median: bool = False) -> np.ndarray:
    """Draw a Latin hypercube design on the unit cube (``LatinHypercube.Random``/``Median``).

    Parameters
    ----------
    sample_size : int
        Number of rows.
    dimension : int, optional
        Number of columns; default 1.
    seed : int or MersenneTwister, optional
        Seed, or a generator to draw from (the ``Random`` overloads); zero or negative seeds
        from the clock (not reproducible).
    median : bool, optional
        Use bin midpoints instead of uniform positions within bins; default False.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(sample_size, dimension)``.
    """
    lhs = _sampling().LatinHypercube
    method = lhs.Median if median else lhs.Random
    generator = seed.net if isinstance(seed, MersenneTwister) else seed
    if not isinstance(generator, int):
        return to_numpy_matrix(method(int(sample_size), int(dimension), generator))
    return to_numpy_matrix(method(int(sample_size), int(dimension), int(generator)))


class MersenneTwister:
    """The Numerics Mersenne Twister (``Numerics.Sampling.MersenneTwister``).

    Exposed for cross-checks of NumPy re-implementations against the .NET stream, and as the
    generator passed to Numerics members that take a ``System.Random``.

    Parameters
    ----------
    seed : int or sequence of int, optional
        Seed, or an array of seeds (``init_by_array``); default the clock (not reproducible).

    Attributes
    ----------
    net : Numerics.Sampling.MersenneTwister
        The wrapped .NET generator.
    """

    def __init__(self, seed: Any = None) -> None:
        s = _sampling()
        if seed is None:
            self.net = s.MersenneTwister()
        elif isinstance(seed, (int, np.integer)):
            self.net = s.MersenneTwister(int(seed))
        else:
            self.net = s.MersenneTwister(to_net_ints(seed))

    def initialize(self, seed: int) -> None:
        """Reseed the generator with an unsigned 32-bit seed (``Initialize(uint)``).

        Parameters
        ----------
        seed : int
            Seed in [0, 2**32).
        """
        import System  # noqa: PLC0415 - the runtime is running

        self.net.Initialize(System.UInt32(int(seed)))

    def next_double(self) -> float:
        """Return the next double in [0, 1) (``NextDouble``).

        Returns
        -------
        float
            The draw.
        """
        return float(self.net.NextDouble())

    def next_doubles(self, length: int) -> np.ndarray:
        """Return the next doubles (``ExtensionMethods.NextDoubles``).

        Parameters
        ----------
        length : int
            Number of draws.

        Returns
        -------
        numpy.ndarray
            The draws.
        """
        root = namespace("Numerics")
        return to_numpy_doubles(root.ExtensionMethods.NextDoubles(self.net, int(length)))

    def next_integers(self, length: int) -> np.ndarray:
        """Return the next non-negative integers (``ExtensionMethods.NextIntegers``).

        Parameters
        ----------
        length : int
            Number of draws.

        Returns
        -------
        numpy.ndarray
            The draws as int32.
        """
        root = namespace("Numerics")
        return to_numpy_ints(root.ExtensionMethods.NextIntegers(self.net, int(length)))

    def gen_rand_int32(self) -> int:
        """Return the next raw 32-bit output (``GenRandInt32``).

        Returns
        -------
        int
            An unsigned 32-bit integer.
        """
        return int(self.net.GenRandInt32())

    def gen_rand_int31(self) -> int:
        """Return the next 31-bit output (``GenRandInt31``).

        Returns
        -------
        int
            An integer in [0, 2**31).
        """
        return int(self.net.GenRandInt31())

    def gen_rand_real1(self) -> float:
        """Return the next double in [0, 1] (``GenRandReal1``).

        Returns
        -------
        float
            The draw.
        """
        return float(self.net.GenRandReal1())

    def gen_rand_real2(self) -> float:
        """Return the next double in [0, 1) (``GenRandReal2``).

        Returns
        -------
        float
            The draw.
        """
        return float(self.net.GenRandReal2())

    def gen_rand_real3(self) -> float:
        """Return the next double in (0, 1) (``GenRandReal3``).

        Returns
        -------
        float
            The draw.
        """
        return float(self.net.GenRandReal3())

    def gen_rand_res53(self) -> float:
        """Return the next double in [0, 1) with 53-bit resolution (``GenRandRes53``).

        Returns
        -------
        float
            The draw.
        """
        return float(self.net.GenRandRes53())


class SobolSequence:
    """The Sobol low-discrepancy sequence (``Numerics.Sampling.SobolSequence``).

    Parameters
    ----------
    dimension : int, optional
        Space dimension; default 1.
    seed : int, optional
        Seed of the linear matrix scramble and digital shift; default None, the original
        unrandomized sequence.

    Attributes
    ----------
    net : Numerics.Sampling.SobolSequence
        The wrapped .NET sequence.
    """

    def __init__(self, dimension: int = 1, seed: int | None = None) -> None:
        s = _sampling()
        self.net = s.SobolSequence(int(dimension)) if seed is None else s.SobolSequence(int(dimension), int(seed))

    @property
    def dimension(self) -> int:
        """The space dimension (``Dimension``)."""
        return int(self.net.Dimension)

    @property
    def seed(self) -> int | None:
        """The scramble seed, or None for the unrandomized sequence (``Seed``)."""
        value = self.net.Seed
        return None if value is None else int(value)

    def next_double(self) -> np.ndarray:
        """Return the next point of the sequence (``NextDouble``).

        Returns
        -------
        numpy.ndarray
            A point in [0, 1)^dimension.
        """
        return to_numpy_doubles(self.net.NextDouble())

    def next_doubles(self, length: int) -> np.ndarray:
        """Return the next points, one ``NextDouble`` each.

        Parameters
        ----------
        length : int
            Number of points.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(length, dimension)``.
        """
        out = np.empty((int(length), self.dimension), dtype=np.float64)
        for row in range(int(length)):
            out[row] = to_numpy_doubles(self.net.NextDouble())
        return out

    def skip_to(self, index: int) -> np.ndarray:
        """Move to a position of the sequence and return the point there (``SkipTo``).

        Indices count the origin as point 0, which Numerics never returns: ``skip_to(i)`` gives
        point ``i`` of ``scipy.stats.qmc.Sobol(d, scramble=False)`` for ``i >= 1``, ``skip_to(0)``
        gives point 1, and :meth:`next_double` continues from there.

        Parameters
        ----------
        index : int
            Position in the sequence.

        Returns
        -------
        numpy.ndarray
            The point at ``index``.
        """
        return to_numpy_doubles(self.net.SkipTo(int(index)))


class BootstrapCIMethod(NetEnum):
    """Bootstrap confidence interval methods (``Numerics.Sampling.BootstrapCIMethod``)."""

    PERCENTILE = "Percentile"
    BIAS_CORRECTED = "BiasCorrected"
    BCA = "BCa"
    NORMAL = "Normal"
    BOOTSTRAP_T = "BootstrapT"


class PivotalBootstrapInvalidDrawPolicy(NetEnum):
    """Handling of invalid pivotal draws (``Numerics.Sampling.PivotalBootstrapInvalidDrawPolicy``)."""

    DROP = "Drop"
    USE_RAW = "UseRaw"
    USE_PARENT = "UseParent"


def _parameter_set(values: Any, fitness: float = float("nan")) -> Any:
    """Build a .NET ``ParameterSet`` from values, or pass a .NET one through."""
    ps = namespace("Numerics.Mathematics.Optimization").ParameterSet
    if isinstance(values, ps):
        return values
    return ps(to_net_doubles(np.ravel(np.asarray(values, dtype=np.float64))), float(fitness))


class BootstrapFit:
    """A fitted parameter set with its covariance (``Numerics.Sampling.BootstrapFit``).

    Parameters
    ----------
    parameters : array_like or Numerics.Sampling.BootstrapFit
        The fitted parameters, or a .NET fit to wrap.
    covariance : array_like, optional
        The parameter covariance, shape ``(p, p)``; required unless wrapping.

    Attributes
    ----------
    net : Numerics.Sampling.BootstrapFit
        The wrapped .NET object.
    """

    def __init__(self, parameters: Any, covariance: Any = None) -> None:
        s = _sampling()
        if isinstance(parameters, s.BootstrapFit):
            self.net = parameters
            return
        if covariance is None:
            raise ValueError("BootstrapFit needs a covariance matrix.")
        matrix = _numerics_matrix(np.atleast_2d(np.asarray(covariance, dtype=np.float64)))
        values = np.ravel(np.asarray(parameters, dtype=np.float64))
        self.net = s.BootstrapFit(to_net_doubles(values), matrix)

    def __repr__(self) -> str:
        return f"BootstrapFit(parameters={self.parameters!r})"

    @property
    def parameters(self) -> np.ndarray:
        """The fitted parameter values (``Parameters.Values``)."""
        return to_numpy_doubles(self.net.Parameters.Values)

    @property
    def covariance(self) -> np.ndarray:
        """The parameter covariance (``Covariance``)."""
        return to_numpy_matrix(self.net.Covariance.ToArray())

    @property
    def parameter_count(self) -> int:
        """The number of parameters (``ParameterCount``)."""
        return int(self.net.ParameterCount)


class PivotalBootstrapContext:
    """The parent fit and raw fits a pivotal link factory sees (``Numerics.Sampling.PivotalBootstrapContext``).

    Parameters
    ----------
    parent_fit : BootstrapFit or Numerics.Sampling.PivotalBootstrapContext
        The parent fit, or a .NET context to wrap.
    raw_bootstrap_fits : sequence of BootstrapFit, optional
        The accepted raw bootstrap fits; required unless wrapping.

    Attributes
    ----------
    net : Numerics.Sampling.PivotalBootstrapContext
        The wrapped .NET object.
    """

    def __init__(self, parent_fit: Any, raw_bootstrap_fits: Sequence[BootstrapFit] | None = None) -> None:
        s = _sampling()
        if isinstance(parent_fit, s.PivotalBootstrapContext):
            self.net = parent_fit
            return
        import System  # noqa: PLC0415 - the runtime is running

        fits = System.Array[s.BootstrapFit]([f.net for f in (raw_bootstrap_fits or ())])
        self.net = s.PivotalBootstrapContext(parent_fit.net, fits)

    @property
    def parent_fit(self) -> BootstrapFit:
        """The parent fit (``ParentFit``)."""
        return BootstrapFit(self.net.ParentFit)

    @property
    def raw_bootstrap_fits(self) -> list[BootstrapFit]:
        """The accepted raw bootstrap fits (``RawBootstrapFits``)."""
        return [BootstrapFit(f) for f in self.net.RawBootstrapFits]

    @property
    def parameter_count(self) -> int:
        """The number of parameters (``ParameterCount``)."""
        return int(self.net.ParameterCount)

    def get_raw_parameter_values(self, parameter_index: int) -> np.ndarray:
        """Return one parameter across the raw fits (``GetRawParameterValues``).

        Parameters
        ----------
        parameter_index : int
            Zero-based parameter index.

        Returns
        -------
        numpy.ndarray
            The values, one per raw fit.
        """
        return to_numpy_doubles(self.net.GetRawParameterValues(int(parameter_index)))


class PivotalBootstrapDiagnostics:
    """Counts and timings of a pivotal bootstrap run (``Numerics.Sampling.PivotalBootstrapDiagnostics``).

    Parameters
    ----------
    net : Numerics.Sampling.PivotalBootstrapDiagnostics, optional
        The .NET object to wrap; default a new empty one.

    Attributes
    ----------
    net : Numerics.Sampling.PivotalBootstrapDiagnostics
        The wrapped .NET object.
    """

    def __init__(self, net: Any = None) -> None:
        self.net = _sampling().PivotalBootstrapDiagnostics() if net is None else net

    @property
    def requested_replicates(self) -> int:
        """Replicates requested (``RequestedReplicates``)."""
        return int(self.net.RequestedReplicates)

    @property
    def rejected_raw_replicates(self) -> int:
        """Raw fits rejected by the filter (``RejectedRawReplicates``)."""
        return int(self.net.RejectedRawReplicates)

    @property
    def failed_raw_replicates(self) -> int:
        """Raw fits that failed (``FailedRawReplicates``)."""
        return int(self.net.FailedRawReplicates)

    @property
    def accepted_raw_replicates(self) -> int:
        """Raw fits accepted (``AcceptedRawReplicates``)."""
        return int(self.net.AcceptedRawReplicates)

    @property
    def invalid_pivotal_replicates(self) -> int:
        """Pivotal draws found invalid (``InvalidPivotalReplicates``)."""
        return int(self.net.InvalidPivotalReplicates)

    @property
    def retained_pivotal_replicates(self) -> int:
        """Pivotal draws retained (``RetainedPivotalReplicates``)."""
        return int(self.net.RetainedPivotalReplicates)

    @property
    def resampling_time(self) -> float:
        """Seconds spent resampling and fitting (``ResamplingTime``)."""
        return float(self.net.ResamplingTime.TotalSeconds)

    @property
    def transformation_time(self) -> float:
        """Seconds spent on the pivotal transformation (``TransformationTime``)."""
        return float(self.net.TransformationTime.TotalSeconds)


class BootstrapStatisticResult:
    """The interval of one statistic or parameter (``Numerics.Sampling.BootstrapStatisticResult``).

    Parameters
    ----------
    net : Numerics.Sampling.BootstrapStatisticResult, optional
        The .NET object to wrap; default a new empty one.

    Attributes
    ----------
    net : Numerics.Sampling.BootstrapStatisticResult
        The wrapped .NET object.
    """

    def __init__(self, net: Any = None) -> None:
        self.net = _sampling().BootstrapStatisticResult() if net is None else net

    def __repr__(self) -> str:
        return (f"BootstrapStatisticResult(estimate={self.population_estimate!r}, lower_ci={self.lower_ci!r}, "
                f"upper_ci={self.upper_ci!r}, standard_error={self.standard_error!r})")

    @property
    def population_estimate(self) -> float:
        """The estimate from the original data (``PopulationEstimate``)."""
        return float(self.net.PopulationEstimate)

    @property
    def lower_ci(self) -> float:
        """Lower confidence bound (``LowerCI``)."""
        return float(self.net.LowerCI)

    @property
    def upper_ci(self) -> float:
        """Upper confidence bound (``UpperCI``)."""
        return float(self.net.UpperCI)

    @property
    def valid_count(self) -> int:
        """Replicates with a finite value (``ValidCount``)."""
        return int(self.net.ValidCount)

    @property
    def total_count(self) -> int:
        """Replicates requested (``TotalCount``)."""
        return int(self.net.TotalCount)

    @property
    def standard_error(self) -> float:
        """Bootstrap standard error (``StandardError``)."""
        return float(self.net.StandardError)

    @property
    def mean(self) -> float:
        """Mean of the bootstrap replicates (``Mean``)."""
        return float(self.net.Mean)


class BootstrapResults:
    """Confidence intervals of a bootstrap run (``Numerics.Sampling.BootstrapResults``).

    Parameters
    ----------
    net : Numerics.Sampling.BootstrapResults, optional
        The .NET object to wrap; default a new empty one.

    Attributes
    ----------
    net : Numerics.Sampling.BootstrapResults
        The wrapped .NET object.
    """

    def __init__(self, net: Any = None) -> None:
        self.net = _sampling().BootstrapResults() if net is None else net

    @property
    def method(self) -> BootstrapCIMethod:
        """The interval method (``Method``)."""
        return BootstrapCIMethod(str(self.net.Method))

    @property
    def alpha(self) -> float:
        """The two-sided level (``Alpha``)."""
        return float(self.net.Alpha)

    @property
    def statistic_results(self) -> list[BootstrapStatisticResult]:
        """The interval of each statistic by the requested method (``StatisticResults``)."""
        return [BootstrapStatisticResult(r) for r in self.net.StatisticResults]

    @property
    def parameter_results(self) -> list[BootstrapStatisticResult]:
        """The percentile interval of each parameter, whatever the method (``ParameterResults``)."""
        return [BootstrapStatisticResult(r) for r in self.net.ParameterResults]

    @property
    def failed_replicates(self) -> int:
        """Replicates that failed after every retry (``FailedReplicates``)."""
        return int(self.net.FailedReplicates)


_PROGRESS: list[Any] = []


def _progress_type() -> Any:
    """Return a .NET ``IProgress<int>`` implementation that forwards ``Report`` to a Python callable, built once."""
    if not _PROGRESS:
        _sampling()
        import System  # noqa: PLC0415 - the runtime is running

        class CallableProgress(System.IProgress[System.Int32]):
            __namespace__ = "Pynumerics.Sampling"

            def __init__(self, fn: Any) -> None:
                self.fn = fn

            def Report(self, value: Any) -> None:  # noqa: N802 - the .NET interface member
                self.fn(int(value))

        _PROGRESS.append(CallableProgress)
    return _PROGRESS[0]


def _is_python(fn: Any) -> bool:
    """Return whether ``fn`` is a Python callable (not C# source or a .NET delegate)."""
    return not isinstance(fn, str) and callbacks.uses_python_callback(fn) and not _is_net_delegate(fn)


def _is_net_delegate(fn: Any) -> bool:
    """Return whether ``fn`` is a .NET delegate."""
    import System  # noqa: PLC0415 - the runtime is running

    return isinstance(fn, System.Delegate)


class Bootstrap:
    """The general bootstrap (``Numerics.Sampling.Bootstrap<TData>``).

    The analysis is defined by functions: ``resample`` draws a data set, ``fit`` fits the model
    and ``statistic`` derives the quantities of interest, replicated ``replicates`` times. Each
    function is a Python callable, C# source compiled with Roslyn or a .NET delegate of the
    ``Func`` type named below. Python callables see NumPy arrays: a ``ParameterSet`` argument
    arrives as its values, and ``fit`` may return parameter values (built into a ``ParameterSet``
    with NaN fitness). C# source and delegates use the .NET types directly.

    Runs are reproducible: replicate i draws from its own Mersenne Twister seeded from ``seed``
    (``PRNGSeed``), whatever the thread order. Numerics runs the replicates with ``Parallel.For``
    while :attr:`use_parallel` is True. A Python callable sets it to False (the callback policy of
    :mod:`pynumerics.callbacks`), so the replicates run in order on the calling thread with the same
    bits; C# source and .NET delegates keep the replicates parallel. The leave-one-out jackknife of
    the BCa intervals follows :attr:`use_parallel` too, with the same bits either way. An exception
    in a function fails the attempt, which is retried up to ``max_retries`` times (see
    :attr:`failed_attempts`).

    Parameters
    ----------
    data : array_like or .NET object
        The original data. A NumPy-convertible 1-D array makes ``TData`` ``double[]``; a .NET
        object makes ``TData`` its type (or ``data_type``).
    original_parameters : array_like or BootstrapFit
        The fit to the original data, or a :class:`BootstrapFit` with covariance for the pivotal
        bootstrap.
    resample : callable, str or System.Delegate, optional
        ``Func<TData, ParameterSet, Random, TData>``: ``resample(data, parameters, prng)``.
    fit : callable, str or System.Delegate, optional
        ``Func<TData, ParameterSet>``: ``fit(data)``.
    statistic : callable, str or System.Delegate, optional
        ``Func<ParameterSet, double[]>``: ``statistic(parameters)``.
    fit_with_covariance : callable, str or System.Delegate, optional
        ``Func<TData, BootstrapFit>`` for the pivotal bootstrap; a Python callable returns a
        :class:`BootstrapFit` or ``(parameters, covariance)``.
    jackknife : callable, str or System.Delegate, optional
        ``Func<TData, int, TData>``: the data without observation i (for BCa).
    sample_size : callable, str or System.Delegate, optional
        ``Func<TData, int>``: the number of observations (for BCa).
    replicates : int, optional
        Number of replicates; default 10000.
    seed : int, optional
        ``PRNGSeed``; default 12345.
    data_type : type, optional
        The .NET type ``TData`` for .NET data.

    Attributes
    ----------
    net : Numerics.Sampling.Bootstrap
        The wrapped .NET bootstrap.
    """

    def __init__(self, data: Any, original_parameters: Any, *, resample: Any = None, fit: Any = None,
                 statistic: Any = None, fit_with_covariance: Any = None, jackknife: Any = None,
                 sample_size: Any = None, replicates: int | None = None, seed: int | None = None,
                 data_type: Any = None) -> None:
        s = _sampling()
        import System  # noqa: PLC0415 - the runtime is running

        if data_type is None and isinstance(data, System.Object) and not isinstance(data, System.Array[System.Double]):
            data_type = type(data)
        if data_type is None:
            self._data_type = System.Array[System.Double]
            net_data = data if isinstance(data, System.Array[System.Double]) else to_net_doubles(
                np.ravel(np.asarray(data, dtype=np.float64)))
        else:
            self._data_type, net_data = data_type, data
        self._is_array = self._data_type is System.Array[System.Double]
        boot = s.Bootstrap[self._data_type]
        if isinstance(original_parameters, BootstrapFit):
            self.net = boot(net_data, original_parameters.net)
        else:
            self.net = boot(net_data, _parameter_set(original_parameters))
        for name, value in (("resample_function", resample), ("fit_function", fit), ("statistic_function", statistic),
                            ("fit_with_covariance_function", fit_with_covariance),
                            ("jackknife_function", jackknife), ("sample_size_function", sample_size)):
            if value is not None:
                setattr(self, name, value)
        if replicates is not None:
            self.replicates = replicates
        if seed is not None:
            self.seed = seed

    def _data_in(self, value: Any) -> Any:
        """Return a data argument as Python sees it."""
        return to_numpy_doubles(value) if self._is_array else value

    def _data_out(self, value: Any) -> Any:
        """Return a Python data result as .NET ``TData``."""
        import System  # noqa: PLC0415 - the runtime is running

        if self._is_array and not isinstance(value, System.Array[System.Double]):
            return to_net_doubles(np.ravel(np.asarray(value, dtype=np.float64)))
        return value

    def _func(self, *types: Any) -> Any:
        """Return the ``System.Func`` delegate type of the given argument and result types."""
        import System  # noqa: PLC0415 - the runtime is running

        return System.Func[types]

    def _types(self) -> tuple[Any, Any, Any, Any]:
        """Return ``TData``, ``ParameterSet``, ``BootstrapFit`` and ``System.Random``."""
        import System  # noqa: PLC0415 - the runtime is running

        ps = namespace("Numerics.Mathematics.Optimization").ParameterSet
        return self._data_type, ps, _sampling().BootstrapFit, System.Random

    def _set(self, prop: str, fn: Any, dtype: Any, adapter: Any) -> None:
        """Set a delegate property from a Python callable (through ``adapter``), C# source or delegate."""
        if fn is None:
            setattr(self.net, prop, None)
            return
        setattr(self.net, prop, dtype(adapter(fn)) if _is_python(fn) else callbacks.to_delegate(fn, dtype))
        callbacks.apply_parallel_policy(self.net, fn)

    @property
    def resample_function(self) -> Any:
        """``resample(data, parameters, prng) -> data`` (``ResampleFunction``)."""
        return self.net.ResampleFunction

    @resample_function.setter
    def resample_function(self, fn: Any) -> None:
        data, ps, _, rnd = self._types()
        self._set("ResampleFunction", fn, self._func(data, ps, rnd, data),
                  lambda f: lambda d, p, r: self._data_out(f(self._data_in(d), to_numpy_doubles(p.Values), r)))

    @property
    def fit_function(self) -> Any:
        """``fit(data) -> parameters`` (``FitFunction``)."""
        return self.net.FitFunction

    @fit_function.setter
    def fit_function(self, fn: Any) -> None:
        data, ps, _, _ = self._types()
        self._set("FitFunction", fn, self._func(data, ps), lambda f: lambda d: _parameter_set(f(self._data_in(d))))

    @property
    def fit_with_covariance_function(self) -> Any:
        """``fit_with_covariance(data) -> BootstrapFit`` (``FitWithCovarianceFunction``)."""
        return self.net.FitWithCovarianceFunction

    @fit_with_covariance_function.setter
    def fit_with_covariance_function(self, fn: Any) -> None:
        data, _, bf, _ = self._types()

        def adapt(f: Any) -> Any:
            def call(d: Any) -> Any:
                value = f(self._data_in(d))
                return value.net if isinstance(value, BootstrapFit) else BootstrapFit(*value).net
            return call

        self._set("FitWithCovarianceFunction", fn, self._func(data, bf), adapt)

    @property
    def statistic_function(self) -> Any:
        """``statistic(parameters) -> values`` (``StatisticFunction``)."""
        return self.net.StatisticFunction

    @statistic_function.setter
    def statistic_function(self, fn: Any) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        _, ps, _, _ = self._types()
        self._set("StatisticFunction", fn, self._func(ps, System.Array[System.Double]),
                  lambda f: lambda p: to_net_doubles(np.ravel(np.asarray(f(to_numpy_doubles(p.Values)),
                                                                         dtype=np.float64))))

    @property
    def jackknife_function(self) -> Any:
        """``jackknife(data, index) -> data`` without observation ``index`` (``JackknifeFunction``)."""
        return self.net.JackknifeFunction

    @jackknife_function.setter
    def jackknife_function(self, fn: Any) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        data = self._data_type
        self._set("JackknifeFunction", fn, self._func(data, System.Int32, data),
                  lambda f: lambda d, i: self._data_out(f(self._data_in(d), int(i))))

    @property
    def sample_size_function(self) -> Any:
        """``sample_size(data) -> int`` (``SampleSizeFunction``)."""
        return self.net.SampleSizeFunction

    @sample_size_function.setter
    def sample_size_function(self, fn: Any) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        self._set("SampleSizeFunction", fn, self._func(self._data_type, System.Int32),
                  lambda f: lambda d: int(f(self._data_in(d))))

    @property
    def transform(self) -> Any:
        """Transform applied before the Normal and Bootstrap-t intervals (``Transform``); default cube root."""
        return self.net.Transform

    @transform.setter
    def transform(self, fn: Any) -> None:
        self.net.Transform = callbacks.scalar_function(fn)

    @property
    def inverse_transform(self) -> Any:
        """Inverse of :attr:`transform` (``InverseTransform``); default cube."""
        return self.net.InverseTransform

    @inverse_transform.setter
    def inverse_transform(self, fn: Any) -> None:
        self.net.InverseTransform = callbacks.scalar_function(fn)

    @property
    def pivotal_link_factory(self) -> Any:
        """``factory(context) -> links``, one ``ILinkFunction`` or None per parameter (``PivotalLinkFactory``)."""
        return self.net.PivotalLinkFactory

    @pivotal_link_factory.setter
    def pivotal_link_factory(self, fn: Any) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        link = namespace("Numerics.Functions").ILinkFunction
        s = _sampling()

        def adapt(f: Any) -> Any:
            return lambda context: System.Array[link](list(f(PivotalBootstrapContext(context))))

        self._set("PivotalLinkFactory", fn, self._func(s.PivotalBootstrapContext, System.Array[link]), adapt)

    @property
    def pivotal_replicate_filter(self) -> Any:
        """``filter(fit) -> bool`` applied to valid raw fits (``PivotalReplicateFilter``)."""
        return self.net.PivotalReplicateFilter

    @pivotal_replicate_filter.setter
    def pivotal_replicate_filter(self, fn: Any) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        self._set("PivotalReplicateFilter", fn, self._func(_sampling().BootstrapFit, System.Boolean),
                  lambda f: lambda fit: bool(f(BootstrapFit(fit))))

    @property
    def pivotal_parameter_validator(self) -> Any:
        """``validator(values) -> bool`` applied to pivotal draws (``PivotalParameterValidator``)."""
        return self.net.PivotalParameterValidator

    @pivotal_parameter_validator.setter
    def pivotal_parameter_validator(self, fn: Any) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        self._set("PivotalParameterValidator", fn, self._func(System.Array[System.Double], System.Boolean),
                  lambda f: lambda v: bool(f(to_numpy_doubles(v))))

    @property
    def replicates(self) -> int:
        """Number of replicates (``Replicates``)."""
        return int(self.net.Replicates)

    @replicates.setter
    def replicates(self, value: int) -> None:
        self.net.Replicates = int(value)

    @property
    def seed(self) -> int:
        """The PRNG seed (``PRNGSeed``)."""
        return int(self.net.PRNGSeed)

    @seed.setter
    def seed(self, value: int) -> None:
        self.net.PRNGSeed = int(value)

    @property
    def max_retries(self) -> int:
        """Attempts per replicate before it counts as failed (``MaxRetries``)."""
        return int(self.net.MaxRetries)

    @max_retries.setter
    def max_retries(self, value: int) -> None:
        self.net.MaxRetries = int(value)

    @property
    def inner_replicates(self) -> int:
        """Inner replicates of the Bootstrap-t standard errors (``InnerReplicates``)."""
        return int(self.net.InnerReplicates)

    @inner_replicates.setter
    def inner_replicates(self, value: int) -> None:
        self.net.InnerReplicates = int(value)

    @property
    def use_parallel(self) -> bool:
        """Whether the outer replicates run on several threads (``UseParallel``).

        Numerics defaults to True; giving a Python callable turns it off. False runs the replicates
        in order on the calling thread, with the same per-replicate seeds and bit-identical results.
        """
        return bool(self.net.UseParallel)

    @use_parallel.setter
    def use_parallel(self, value: bool) -> None:
        self.net.UseParallel = bool(value)

    @property
    def progress(self) -> Any:
        """The ``IProgress<int>`` receiving completed outer replicates, or None (``Progress``).

        A Python callable ``progress(count)`` is wrapped in an ``IProgress<int>`` whose ``Report``
        calls it synchronously from the worker threads, as Numerics reports (counts can arrive out
        of order); an exception it raises fails the run.
        """
        return self.net.Progress

    @progress.setter
    def progress(self, value: Any) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        if value is not None and callable(value) and not isinstance(value, System.Object):
            value = _progress_type()(value)
        self.net.Progress = value

    @property
    def cancellation_token(self) -> Any:
        """A ``System.Threading.CancellationToken`` that stops the run (``CancellationToken``)."""
        return self.net.CancellationToken

    @cancellation_token.setter
    def cancellation_token(self, value: Any) -> None:
        self.net.CancellationToken = value

    @property
    def original_covariance(self) -> np.ndarray | None:
        """The covariance of the original fit for the pivotal bootstrap (``OriginalCovariance``)."""
        value = self.net.OriginalCovariance
        return None if value is None else to_numpy_matrix(value.ToArray())

    @original_covariance.setter
    def original_covariance(self, value: Any) -> None:
        self.net.OriginalCovariance = None if value is None else _numerics_matrix(
            np.atleast_2d(np.asarray(value, dtype=np.float64)))

    @property
    def pivotal_invalid_draw_policy(self) -> PivotalBootstrapInvalidDrawPolicy:
        """Handling of invalid pivotal draws (``PivotalInvalidDrawPolicy``)."""
        return PivotalBootstrapInvalidDrawPolicy(str(self.net.PivotalInvalidDrawPolicy))

    @pivotal_invalid_draw_policy.setter
    def pivotal_invalid_draw_policy(self, value: PivotalBootstrapInvalidDrawPolicy | str) -> None:
        self.net.PivotalInvalidDrawPolicy = net_enum(_sampling().PivotalBootstrapInvalidDrawPolicy, value)

    @property
    def regularize_pivotal_covariances(self) -> bool:
        """Whether pivotal covariances are regularized before Cholesky (``RegularizePivotalCovariances``)."""
        return bool(self.net.RegularizePivotalCovariances)

    @regularize_pivotal_covariances.setter
    def regularize_pivotal_covariances(self, value: bool) -> None:
        self.net.RegularizePivotalCovariances = bool(value)

    @property
    def pivotal_z_limit(self) -> float | None:
        """Absolute limit of the standardized pivotal components, or None (``PivotalZLimit``)."""
        value = self.net.PivotalZLimit
        return None if value is None else float(value)

    @pivotal_z_limit.setter
    def pivotal_z_limit(self, value: float | None) -> None:
        import System  # noqa: PLC0415 - the runtime is running

        self.net.PivotalZLimit = None if value is None else System.Nullable[System.Double](float(value))

    @property
    def add_pivotal_jitter(self) -> bool:
        """Whether Gaussian jitter is added to pivotal vectors (``AddPivotalJitter``)."""
        return bool(self.net.AddPivotalJitter)

    @add_pivotal_jitter.setter
    def add_pivotal_jitter(self, value: bool) -> None:
        self.net.AddPivotalJitter = bool(value)

    @property
    def pivotal_jitter_scale(self) -> float:
        """Base jitter standard deviation, divided by sqrt(p) when applied (``PivotalJitterScale``)."""
        return float(self.net.PivotalJitterScale)

    @pivotal_jitter_scale.setter
    def pivotal_jitter_scale(self, value: float) -> None:
        self.net.PivotalJitterScale = float(value)

    @property
    def bootstrap_parameter_sets(self) -> np.ndarray:
        """Replicate parameters, shape ``(replicates, parameters)`` (``BootstrapParameterSets``).

        There are as many columns as the longest parameter set. A shorter set is padded with NaN, and a
        failed replicate, which Numerics fills with NaN values of the full length, reads as NaN throughout.
        """
        return _parameter_set_arrays(self.net.BootstrapParameterSets)[0]

    @property
    def bootstrap_statistics(self) -> np.ndarray:
        """The replicate statistics, shape ``(replicates, statistics)`` (``BootstrapStatistics``)."""
        return to_numpy_matrix(self.net.BootstrapStatistics)

    @property
    def raw_bootstrap_fits(self) -> list[BootstrapFit]:
        """The accepted raw fits of the last pivotal run (``RawBootstrapFits``)."""
        return [BootstrapFit(f) for f in self.net.RawBootstrapFits]

    @property
    def raw_bootstrap_parameter_sets(self) -> np.ndarray:
        """The accepted raw parameters of the last pivotal run (``RawBootstrapParameterSets``).

        There are as many columns as the longest parameter set; a shorter set is padded with NaN.
        """
        return _parameter_set_arrays(self.net.RawBootstrapParameterSets)[0]

    @property
    def raw_bootstrap_statistics(self) -> np.ndarray | None:
        """The raw statistics of the last pivotal run, or None (``RawBootstrapStatistics``)."""
        value = self.net.RawBootstrapStatistics
        return None if value is None else to_numpy_matrix(value)

    @property
    def pivotal_links(self) -> list[Any]:
        """The link functions of the last pivotal run (``PivotalLinks``)."""
        return list(self.net.PivotalLinks)

    @property
    def pivotal_diagnostics(self) -> PivotalBootstrapDiagnostics | None:
        """Diagnostics of the last pivotal run, or None (``PivotalDiagnostics``)."""
        value = self.net.PivotalDiagnostics
        return None if value is None else PivotalBootstrapDiagnostics(value)

    @property
    def failed_replicates(self) -> int:
        """Replicates that failed after every retry (``FailedReplicates``)."""
        return int(self.net.FailedReplicates)

    @property
    def failed_attempts(self) -> int:
        """Attempts that raised in the last run (``FailedAttempts``)."""
        return int(self.net.FailedAttempts)

    @property
    def failed_jackknife_replicates(self) -> int:
        """Leave-one-out replicates that failed in the last BCa request (``FailedJackknifeReplicates``)."""
        return int(self.net.FailedJackknifeReplicates)

    def run(self) -> Bootstrap:
        """Run the regular bootstrap (``Run``).

        Returns
        -------
        Bootstrap
            This object, for chaining :meth:`get_confidence_intervals`.
        """
        self.net.Run()
        return self

    def run_double_bootstrap(self, inner_replicates: int = 300) -> Bootstrap:
        """Run the double bootstrap with bias correction (``RunDoubleBootstrap``).

        Parameters
        ----------
        inner_replicates : int, optional
            Inner replicates per outer replicate; default 300.

        Returns
        -------
        Bootstrap
            This object.
        """
        self.net.RunDoubleBootstrap(int(inner_replicates))
        return self

    def run_with_studentized_bootstrap(self) -> Bootstrap:
        """Run the bootstrap with inner standard errors for Bootstrap-t (``RunWithStudentizedBootstrap``).

        Returns
        -------
        Bootstrap
            This object.
        """
        self.net.RunWithStudentizedBootstrap()
        return self

    def run_pivotal_bootstrap(self) -> Bootstrap:
        """Run the covariance-aware pivotal bootstrap (``RunPivotalBootstrap``).

        Returns
        -------
        Bootstrap
            This object.
        """
        self.net.RunPivotalBootstrap()
        return self

    def transform_pivotal_bootstrap(self, raw_fits: Sequence[BootstrapFit]) -> Bootstrap:
        """Apply the pivotal transformation to given raw fits (``TransformPivotalBootstrap``).

        Parameters
        ----------
        raw_fits : sequence of BootstrapFit
            Raw covariance-aware bootstrap fits.

        Returns
        -------
        Bootstrap
            This object.
        """
        import System  # noqa: PLC0415 - the runtime is running

        fits = System.Array[_sampling().BootstrapFit]([getattr(f, "net", f) for f in raw_fits])
        self.net.TransformPivotalBootstrap(fits)
        return self

    def get_confidence_intervals(self, method: BootstrapCIMethod | str = BootstrapCIMethod.PERCENTILE,
                                 alpha: float = 0.1) -> BootstrapResults:
        """Return confidence intervals of the last run (``GetConfidenceIntervals``).

        Parameters
        ----------
        method : BootstrapCIMethod or str, optional
            The interval method; default percentile.
        alpha : float, optional
            Two-sided level; default 0.1 (90 percent intervals).

        Returns
        -------
        BootstrapResults
            Intervals of every parameter and statistic.
        """
        kind = net_enum(_sampling().BootstrapCIMethod, method)
        return BootstrapResults(self.net.GetConfidenceIntervals(kind, float(alpha)))

    def get_raw_pivotal_confidence_intervals(self, alpha: float = 0.1) -> BootstrapResults:
        """Return percentile intervals of the raw fits of the last pivotal run (``GetRawPivotalConfidenceIntervals``).

        Parameters
        ----------
        alpha : float, optional
            Two-sided level; default 0.1.

        Returns
        -------
        BootstrapResults
            The intervals.
        """
        return BootstrapResults(self.net.GetRawPivotalConfidenceIntervals(float(alpha)))


__all__ = [
    "Bootstrap", "BootstrapCIMethod", "BootstrapFit", "BootstrapResults", "BootstrapStatisticResult",
    "MersenneTwister", "PivotalBootstrapContext", "PivotalBootstrapDiagnostics", "PivotalBootstrapInvalidDrawPolicy",
    "SobolSequence", "StratificationBin", "StratificationBins", "StratificationOptions",
    "exceedance_probability_to_x", "latin_hypercube", "probability_to_x", "stratify_exceedance_probabilities",
    "stratify_exceedance_probabilities_many", "stratify_multivariate_probabilities", "stratify_probabilities",
    "stratify_probabilities_many", "stratify_x_values", "stratify_x_values_many", "x_to_exceedance_probability",
    "x_to_probability",
]

translate_public(globals())
