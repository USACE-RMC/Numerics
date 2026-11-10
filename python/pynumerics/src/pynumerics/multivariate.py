"""Multivariate distributions from ``Numerics.Distributions``.

Every class wraps one .NET distribution object (``.net``) and derives from
:class:`MultivariateDistribution`, which carries the members every Numerics multivariate
distribution declares (``PDF``, ``CDF``, ``CCDF`` and their logarithms, ``Dimension``,
``Clone``). Numerics evaluates those members at one point (``double[]``); a two-dimensional
array is evaluated row by row with the bound member, one crossing per row. Samples cross back
as one ``double[,]`` block copy.

Random generators: members that take a ``System.Random`` accept a
:class:`pynumerics.sampling.MersenneTwister` or any .NET ``System.Random``; the ``seed`` argument
of ``generate_random_values`` (and of ``MultivariateNormal.latin_hypercube_random_values``)
accepts an integer seed or such a generator.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics._util import _nd, cached_overload, namespace, net_enum, overload, unwrap
from pynumerics.enums import NetEnum, Transform
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _flat_doubles,
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_ints,
    to_numpy_matrix,
)
from pynumerics.sampling import _bins


class MultivariateDistributionType(NetEnum):
    """Multivariate distribution types (``Numerics.Distributions.MultivariateDistributionType``)."""

    BIVARIATE_EMPIRICAL_DISTRIBUTION = "BivariateEmpiricalDistribution"
    MULTIVARIATE_NORMAL = "MultivariateNormal"
    DIRICHLET = "Dirichlet"
    MULTINOMIAL = "Multinomial"
    MULTIVARIATE_STUDENT_T = "MultivariateStudentT"


def _categorical_indices() -> Any:
    """Return ``Multinomial.CategoricalIndices(double[], double[], bool)``.

    The two-argument overload calls this one with ``inclusive`` false, so binding the
    three-argument overload once serves both forms without a per-call overload resolution.
    """
    return overload(_nd().Multinomial.CategoricalIndices, "Double[]", "Double[]", "Boolean")


def _generator(random: Any) -> Any:
    """Return the .NET ``System.Random`` behind a :class:`pynumerics.sampling.MersenneTwister` or a .NET generator."""
    return getattr(random, "net", random)


def _wrap_generator(net: Any) -> Any:
    """Wrap a .NET generator as :class:`pynumerics.sampling.MersenneTwister` when it is one, else return it."""
    from pynumerics.sampling import MersenneTwister  # noqa: PLC0415 - imported on use

    concrete = unwrap(net)
    if str(concrete.GetType().FullName) != "Numerics.Sampling.MersenneTwister":
        return concrete
    instance = MersenneTwister.__new__(MersenneTwister)
    instance.net = concrete
    return instance


def _seeded_generator(value: Any) -> Any:
    """Return the .NET generator for a generator argument, or a Numerics ``MersenneTwister`` for an integer seed."""
    if isinstance(value, (int, np.integer)):
        return namespace("Numerics.Sampling").MersenneTwister(int(value))
    return _generator(value)


def _draw(member: Any, sample_size: int, seed: Any) -> np.ndarray:
    """Call a ``(int, int seed)`` or ``(int, Random)`` sampling overload and return the ``double[,]`` as NumPy.

    An integer ``seed`` calls the seed overload; any other value is taken as a generator for the
    ``System.Random`` overload.
    """
    import System  # noqa: PLC0415 - the runtime is running (the owner exists)

    if isinstance(seed, (int, np.integer)):
        return to_numpy_matrix(member.Overloads[System.Int32, System.Int32](int(sample_size), int(seed)))
    return to_numpy_matrix(member.Overloads[System.Int32, System.Random](int(sample_size), _generator(seed)))


def _message(error: Any) -> str | None:
    """Return the message of an exception Numerics built without throwing, or None."""
    return None if error is None else str(error.Message)


class MultivariateDistribution:
    """A wrapped ``Numerics.Distributions.MultivariateDistribution``.

    Parameters
    ----------
    net : Numerics.Distributions.MultivariateDistribution
        The .NET distribution; interface and base-class wrappers are unwrapped.

    Attributes
    ----------
    net : Numerics.Distributions.MultivariateDistribution
        The wrapped .NET object.

    Notes
    -----
    The point members accept one point of length ``dimension`` and return a float, or a
    two-dimensional array of points (one per row) and return one value per row.
    """

    def __init__(self, net: Any) -> None:
        self.net = unwrap(net)
        self._bind()

    def _bind(self) -> None:
        """Bind members once after construction; subclasses override it."""

    @classmethod
    def from_net(cls, net: Any) -> MultivariateDistribution:
        """Wrap an existing .NET distribution without copying it.

        Parameters
        ----------
        net : Numerics.Distributions.MultivariateDistribution
            The .NET object; interface and base-class wrappers are unwrapped.

        Returns
        -------
        MultivariateDistribution
            A wrapper of class ``cls`` sharing ``net``.
        """
        instance = cls.__new__(cls)
        MultivariateDistribution.__init__(instance, net)
        return instance

    def _point(self, name: str, x: Any) -> Any:
        """Evaluate a ``(double[])`` member at one point, or row by row for a 2-D array."""
        cache = self.__dict__.setdefault("_points", {})
        member = cache.get(name)
        if member is None:
            import System  # noqa: PLC0415 - the runtime is running (self.net exists)

            group = getattr(self.net, name)
            try:
                member = group.Overloads[System.Array[System.Double]]
            except (KeyError, TypeError):
                member = group
            cache[name] = member
        array = np.asarray(x, dtype=np.float64)
        if array.ndim <= 1:
            return float(member(to_net_doubles(array)))
        if array.ndim != 2:
            raise ValueError(f"Expected one point or a 2-D array of points; got shape {array.shape}.")
        return np.array([float(member(to_net_doubles(row))) for row in array], dtype=np.float64)

    @property
    def dimension(self) -> int:
        """The number of variables."""
        return int(self.net.Dimension)

    @property
    def type(self) -> MultivariateDistributionType:
        """The distribution type (``Type``)."""
        return MultivariateDistributionType(str(self.net.Type))

    @property
    def display_name(self) -> str:
        """The Numerics display name."""
        return str(self.net.DisplayName)

    @property
    def short_display_name(self) -> str:
        """The Numerics short display name."""
        return str(self.net.ShortDisplayName)

    @property
    def parameters_valid(self) -> bool:
        """Whether the current parameters are valid (``ParametersValid``)."""
        return bool(self.net.ParametersValid)

    def pdf(self, x: Any) -> Any:
        """Return the density (``PDF(double[])``).

        Parameters
        ----------
        x : array_like
            One point of length d, or an array of shape (n, d).

        Returns
        -------
        float or numpy.ndarray
            The density at the point, or one per row.
        """
        return self._point("PDF", x)

    def log_pdf(self, x: Any) -> Any:
        """Return the log-density (``LogPDF(double[])``).

        Parameters
        ----------
        x : array_like
            One point of length d, or an array of shape (n, d).

        Returns
        -------
        float or numpy.ndarray
            The log-density at the point, or one per row.
        """
        return self._point("LogPDF", x)

    def cdf(self, x: Any) -> Any:
        """Return the joint CDF P(X ≤ x) (``CDF(double[])``).

        Parameters
        ----------
        x : array_like
            One point of length d, or an array of shape (n, d).

        Returns
        -------
        float or numpy.ndarray
            The joint non-exceedance probability at the point, or one per row.
        """
        return self._point("CDF", x)

    def log_cdf(self, x: Any) -> Any:
        """Return the log of the joint CDF (``LogCDF(double[])``).

        Parameters
        ----------
        x : array_like
            One point of length d, or an array of shape (n, d).

        Returns
        -------
        float or numpy.ndarray
            The log joint non-exceedance probability.
        """
        return self._point("LogCDF", x)

    def ccdf(self, x: Any) -> Any:
        """Return the complement of the joint CDF, 1 − P(X ≤ x) (``CCDF(double[])``).

        Parameters
        ----------
        x : array_like
            One point of length d, or an array of shape (n, d).

        Returns
        -------
        float or numpy.ndarray
            One minus the joint non-exceedance probability.
        """
        return self._point("CCDF", x)

    def log_ccdf(self, x: Any) -> Any:
        """Return the log of the complement of the joint CDF (``LogCCDF(double[])``).

        Parameters
        ----------
        x : array_like
            One point of length d, or an array of shape (n, d).

        Returns
        -------
        float or numpy.ndarray
            The log of one minus the joint non-exceedance probability.
        """
        return self._point("LogCCDF", x)

    def clone(self) -> MultivariateDistribution:
        """Return a deep copy through ``Clone``.

        Returns
        -------
        MultivariateDistribution
            A wrapper of the cloned .NET object, of the matching Python class.
        """
        return wrap(self.net.Clone())


class MultivariateNormal(MultivariateDistribution):
    """Multivariate normal distribution (``Numerics.Distributions.MultivariateNormal``).

    Parameters
    ----------
    mean : array_like or int
        Mean vector μ of length d, or the dimension d for a zero mean.
    covariance : array_like, optional
        Covariance matrix Σ of shape (d, d); default the identity.
    decomposition : str, optional
        Factorization of Σ (``DecompositionMethod`` member name, for example ``"Cholesky"``);
        default the Numerics default.

    Attributes
    ----------
    net : Numerics.Distributions.MultivariateNormal
        The wrapped .NET object.
    """

    def __init__(self, mean: Any, covariance: Any = None, decomposition: str | None = None) -> None:
        nd = _nd()
        is_dimension = isinstance(mean, (int, np.integer))
        args: list[Any] = [int(mean) if is_dimension else to_net_doubles(mean)]
        if covariance is not None:
            if is_dimension:
                raise ValueError("A covariance matrix needs a mean vector, not a dimension.")
            args.append(to_net_matrix(covariance))
        if decomposition is not None:
            linear_algebra = namespace("Numerics.Mathematics.LinearAlgebra")
            args.append(net_enum(linear_algebra.DecompositionMethod, decomposition))
        MultivariateDistribution.__init__(self, nd.MultivariateNormal(*args))

    @classmethod
    def univariate(cls, mu: float, sigma: float) -> MultivariateNormal:
        """Build the one-dimensional normal N(μ, σ²) (``MultivariateNormal.Univariate``).

        Parameters
        ----------
        mu : float
            Mean.
        sigma : float
            Standard deviation.

        Returns
        -------
        MultivariateNormal
            The distribution.
        """
        return cls.from_net(_nd().MultivariateNormal.Univariate(float(mu), float(sigma)))

    @classmethod
    def bivariate(cls, mu1: float, mu2: float, sigma1: float, sigma2: float, rho: float) -> MultivariateNormal:
        """Build a bivariate normal from means, standard deviations and correlation (``Bivariate``).

        Parameters
        ----------
        mu1 : float
            Mean of the first variable.
        mu2 : float
            Mean of the second variable.
        sigma1 : float
            Standard deviation of the first variable.
        sigma2 : float
            Standard deviation of the second variable.
        rho : float
            Correlation.

        Returns
        -------
        MultivariateNormal
            The distribution.
        """
        return cls.from_net(_nd().MultivariateNormal.Bivariate(float(mu1), float(mu2), float(sigma1), float(sigma2),
                                                               float(rho)))

    @staticmethod
    def default_mvnuni_seed() -> int:
        """Return the seed of the default CDF generator (``DefaultMVNUNISeed``).

        Returns
        -------
        int
            The seed every new instance's ``mvnuni`` starts from.
        """
        return int(_nd().MultivariateNormal.DefaultMVNUNISeed)

    @property
    def mvnuni(self) -> Any:
        """The generator of the randomized lattice rule behind the CDF above two dimensions (``MVNUNI``).

        Reads as a :class:`pynumerics.sampling.MersenneTwister` when it is one. Accepts a generator or
        an integer seed (a new Numerics ``MersenneTwister``).
        """
        return _wrap_generator(self.net.MVNUNI)

    @mvnuni.setter
    def mvnuni(self, value: Any) -> None:
        self.net.MVNUNI = _seeded_generator(value)

    @property
    def max_evaluations(self) -> int:
        """The maximum number of integrand evaluations of the CDF above two dimensions; settable."""
        return int(self.net.MaxEvaluations)

    @max_evaluations.setter
    def max_evaluations(self, value: int) -> None:
        self.net.MaxEvaluations = int(value)

    @property
    def absolute_error(self) -> float:
        """The absolute error tolerance of the CDF above two dimensions; settable."""
        return float(self.net.AbsoluteError)

    @absolute_error.setter
    def absolute_error(self, value: float) -> None:
        self.net.AbsoluteError = float(value)

    @property
    def relative_error(self) -> float:
        """The relative error tolerance of the CDF above two dimensions; settable."""
        return float(self.net.RelativeError)

    @relative_error.setter
    def relative_error(self, value: float) -> None:
        self.net.RelativeError = float(value)

    @property
    def mean(self) -> np.ndarray:
        """The mean vector."""
        return to_numpy_doubles(self.net.Mean)

    @property
    def median(self) -> np.ndarray:
        """The median vector (equal to the mean)."""
        return to_numpy_doubles(self.net.Median)

    @property
    def mode(self) -> np.ndarray:
        """The mode vector (equal to the mean)."""
        return to_numpy_doubles(self.net.Mode)

    @property
    def variance(self) -> np.ndarray:
        """The marginal variances (the diagonal of Σ)."""
        return to_numpy_doubles(self.net.Variance)

    @property
    def standard_deviation(self) -> np.ndarray:
        """The marginal standard deviations."""
        return to_numpy_doubles(self.net.StandardDeviation)

    @property
    def covariance(self) -> np.ndarray:
        """The covariance matrix."""
        return to_numpy_matrix(self.net.Covariance)

    @property
    def decomposition(self) -> str:
        """The factorization of Σ in use (``DecompositionMethod`` member name)."""
        return str(self.net.Decomposition)

    @property
    def is_positive_definite(self) -> bool:
        """Whether Σ is positive definite."""
        return bool(self.net.IsPositiveDefinite)

    @property
    def is_density_valid(self) -> bool:
        """Whether the density is defined (Σ is non-singular)."""
        return bool(self.net.IsDensityValid)

    def set_parameters(self, mean: Any, covariance: Any) -> None:
        """Replace the mean and covariance (``SetParameters``).

        Parameters
        ----------
        mean : array_like
            Mean vector of length d.
        covariance : array_like
            Covariance matrix of shape (d, d).
        """
        self.net.SetParameters(to_net_doubles(mean), to_net_matrix(covariance))

    def set_mean(self, mean: Any) -> None:
        """Replace the mean and keep the covariance and its factorization (``SetMean``).

        Parameters
        ----------
        mean : array_like
            Mean vector of length d.
        """
        self.net.SetMean(to_net_doubles(mean))

    def validate_parameters(self, mean: Any, covariance: Any) -> str | None:
        """Return why a mean and covariance are invalid, or None when they are valid.

        Parameters
        ----------
        mean : array_like
            Mean vector.
        covariance : array_like
            Covariance matrix.

        Returns
        -------
        str or None
            The Numerics message of ``ValidateParameters(mean, covariance, false)``, or None.
        """
        return _message(self.net.ValidateParameters(to_net_doubles(mean), to_net_matrix(covariance), False))

    def try_set_parameters(self, mean: Any, covariance: Any) -> bool:
        """Set the parameters when they are valid (``TrySetParameters``).

        Parameters
        ----------
        mean : array_like
            Mean vector.
        covariance : array_like
            Covariance matrix.

        Returns
        -------
        bool
            True when the parameters were set; the distribution is unchanged otherwise.
        """
        return bool(self.net.TrySetParameters(to_net_doubles(mean), to_net_matrix(covariance)))

    def try_set_covariance(self, covariance: Any) -> bool:
        """Set the covariance when it is valid, keeping the mean (``TrySetCovariance``).

        Parameters
        ----------
        covariance : array_like
            Covariance matrix of shape (d, d).

        Returns
        -------
        bool
            True when the covariance was set.
        """
        return bool(self.net.TrySetCovariance(to_net_matrix(covariance)))

    def marginal(self, *indices: int) -> MultivariateNormal:
        """Return the marginal distribution of a subset of the variables (``Marginal``).

        Parameters
        ----------
        *indices : int
            Zero-based variable indices to keep, in the order wanted.

        Returns
        -------
        MultivariateNormal
            The marginal normal.
        """
        return MultivariateNormal.from_net(self.net.Marginal(to_net_ints(list(indices))))

    def conditional(self, observed_indices: Any, observed_values: Any) -> MultivariateNormal:
        """Return the distribution of the other variables given observed values (``Conditional``).

        Parameters
        ----------
        observed_indices : array_like of int
            Zero-based indices of the observed variables.
        observed_values : array_like
            Their observed values.

        Returns
        -------
        MultivariateNormal
            The conditional normal of the unobserved variables, in index order.
        """
        return MultivariateNormal.from_net(self.net.Conditional(to_net_ints(observed_indices),
                                                                to_net_doubles(observed_values)))

    def independent_normal(self, index: int) -> Any:
        """Return the univariate normal marginal of one variable (``IndependentNormal``).

        Parameters
        ----------
        index : int
            Zero-based variable index.

        Returns
        -------
        pynumerics.distributions.Normal
            The marginal.
        """
        from pynumerics.distributions import wrap as wrap_univariate  # noqa: PLC0415 - imported on use

        return wrap_univariate(self.net.IndependentNormal(int(index)))

    def mahalanobis(self, x: Any) -> float:
        """Return the squared Mahalanobis distance (x − μ)ᵀ Σ⁻¹ (x − μ) (``Mahalanobis``).

        Parameters
        ----------
        x : array_like
            Point of length d.

        Returns
        -------
        float
            The squared distance.
        """
        return float(self.net.Mahalanobis(to_net_doubles(x)))

    def interval(self, lower: Any, upper: Any) -> float:
        """Return the probability of the box lower ≤ X ≤ upper (``Interval``).

        Parameters
        ----------
        lower : array_like
            Lower limits of length d.
        upper : array_like
            Upper limits of length d.

        Returns
        -------
        float
            The box probability.
        """
        return float(self.net.Interval(to_net_doubles(lower), to_net_doubles(upper)))

    def inverse_cdf(self, probabilities: Any) -> np.ndarray:
        """Map d marginal probabilities to a correlated point (``InverseCDF(double[])``).

        Parameters
        ----------
        probabilities : array_like
            Probabilities of length d; each maps to a standard normal deviate, and the
            deviates are correlated as A·z + μ.

        Returns
        -------
        numpy.ndarray
            The point of length d.
        """
        return to_numpy_doubles(self.net.InverseCDF(to_net_doubles(probabilities)))

    def mvndst(self, lower: Any, upper: Any, infin: Any, correl: Any, max_points: int, absolute_error: float,
               relative_error: float) -> tuple[float, float, int]:
        """Integrate a standard multivariate normal over a box with Genz's MVNDST (``MVNDST``).

        Parameters
        ----------
        lower : array_like
            Lower limits of length N.
        upper : array_like
            Upper limits of length N.
        infin : array_like of int
            Limit flags: negative for (−∞, ∞), 0 for (−∞, upper], 1 for [lower, ∞), 2 for [lower, upper].
        correl : array_like
            Strict lower triangle of the correlation matrix by rows, N(N − 1)/2 values.
        max_points : int
            Maximum number of integrand evaluations.
        absolute_error : float
            Absolute error tolerance.
        relative_error : float
            Relative error tolerance.

        Returns
        -------
        error : float
            Estimated absolute error (99% confidence).
        value : float
            Estimated probability.
        inform : int
            0 for normal completion, 1 when ``max_points`` was reached first, 2 for N outside 1 to 500.

        Notes
        -----
        The lattice shifts come from this instance's :attr:`mvnuni`, which the call advances.
        """
        low = to_net_doubles(lower)
        error, value, inform = self.net.MVNDST(len(low), low, to_net_doubles(upper), to_net_ints(infin),
                                               to_net_doubles(correl), int(max_points), float(absolute_error),
                                               float(relative_error), 0.0, 0.0, 0)
        return float(error), float(value), int(inform)

    @staticmethod
    def bivariate_cdf(z1: float, z2: float, r: float) -> float:
        """Return the standard bivariate normal upper orthant P(X > z1, Y > z2) (``BivariateCDF``).

        Parameters
        ----------
        z1 : float
            Lower limit of the first variable.
        z2 : float
            Lower limit of the second variable.
        r : float
            Correlation.

        Returns
        -------
        float
            The probability (Drezner and Wesolowsky, as modified by Genz).

        Notes
        -----
        Numerics follows Genz's BVN convention, an upper orthant despite the name; the joint CDF
        P(X ≤ a, Y ≤ b) is ``bivariate_cdf(-a, -b, r)``, which is how Numerics calls it.
        """
        return float(_nd().MultivariateNormal.BivariateCDF(float(z1), float(z2), float(r)))

    @staticmethod
    def bivnor(ah: float, ak: float, r: float) -> float:
        """Return the standard bivariate normal upper orthant P(X ≥ ah, Y ≥ ak) (``bivnor``).

        Parameters
        ----------
        ah : float
            Lower limit of the first variable.
        ak : float
            Lower limit of the second variable.
        r : float
            Correlation.

        Returns
        -------
        float
            The probability (Donnelly, ACM Algorithm 462).
        """
        return float(_nd().MultivariateNormal.bivnor(float(ah), float(ak), float(r)))

    @staticmethod
    def bvu(sh: float, sk: float, r: float) -> float:
        """Return the standard bivariate normal upper orthant P(X > sh, Y > sk) (``BVU``).

        Parameters
        ----------
        sh : float
            Lower limit of the first variable.
        sk : float
            Lower limit of the second variable.
        r : float
            Correlation.

        Returns
        -------
        float
            The probability (Genz).
        """
        return float(_nd().MultivariateNormal.BVU(float(sh), float(sk), float(r)))

    @staticmethod
    def mvnphi(z: float) -> float:
        """Return the standard normal CDF used by MVNDST (``MVNPHI``).

        Parameters
        ----------
        z : float
            The z-score.

        Returns
        -------
        float
            Φ(z).
        """
        return float(_nd().MultivariateNormal.MVNPHI(float(z)))

    @staticmethod
    def phinvs(p: float) -> float:
        """Return the standard normal quantile used by MVNDST (``PHINVS``, algorithm AS 241).

        Parameters
        ----------
        p : float
            Lower tail probability.

        Returns
        -------
        float
            Φ⁻¹(p).
        """
        return float(_nd().MultivariateNormal.PHINVS(float(p)))

    def generate_random_values(self, sample_size: int, seed: Any = -1) -> np.ndarray:
        """Draw a Monte Carlo sample (``GenerateRandomValues``).

        Parameters
        ----------
        sample_size : int
            Number of points.
        seed : int or MersenneTwister, optional
            Seed (zero or negative seeds from the clock, not reproducible) or a generator.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(sample_size, d)``.
        """
        return _draw(self.net.GenerateRandomValues, sample_size, seed)

    def latin_hypercube_random_values(self, sample_size: int, seed: Any) -> np.ndarray:
        """Draw a Latin hypercube sample (``LatinHypercubeRandomValues``).

        Parameters
        ----------
        sample_size : int
            Number of points.
        seed : int or MersenneTwister
            Seed or a generator.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(sample_size, d)``.
        """
        return _draw(self.net.LatinHypercubeRandomValues, sample_size, seed)

    def stratified_random_values(self, stratification_bins: Any, seed: int) -> np.ndarray:
        """Draw a sample stratified in the first variable (``StratifiedRandomValues``).

        Parameters
        ----------
        stratification_bins : pynumerics.sampling.StratificationBins
            Probability bins; one point is drawn per bin at its midpoint.
        seed : int
            Seed of the other variables.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(stratification_bins), d)``.
        """
        return to_numpy_matrix(self.net.StratifiedRandomValues(_bins(stratification_bins), int(seed)))


class MultivariateStudentT(MultivariateDistribution):
    """Multivariate Student-t distribution (``Numerics.Distributions.MultivariateStudentT``).

    Parameters
    ----------
    degrees_of_freedom : float
        Degrees of freedom ν > 0.
    location : array_like or int
        Location vector μ of length d, or the dimension d for a zero location.
    scale_matrix : array_like, optional
        Scale matrix Σ of shape (d, d), positive definite; default the identity. The covariance
        is ν/(ν − 2)·Σ for ν > 2.

    Attributes
    ----------
    net : Numerics.Distributions.MultivariateStudentT
        The wrapped .NET object.
    """

    def __init__(self, degrees_of_freedom: float, location: Any, scale_matrix: Any = None) -> None:
        nd = _nd()
        if isinstance(location, (int, np.integer)):
            if scale_matrix is not None:
                raise ValueError("A scale matrix needs a location vector, not a dimension.")
            net = nd.MultivariateStudentT(int(location), float(degrees_of_freedom))
        elif scale_matrix is None:
            net = nd.MultivariateStudentT(float(degrees_of_freedom), to_net_doubles(location))
        else:
            net = nd.MultivariateStudentT(float(degrees_of_freedom), to_net_doubles(location),
                                          to_net_matrix(scale_matrix))
        MultivariateDistribution.__init__(self, net)

    @property
    def mvnuni(self) -> Any:
        """The generator of the Genz-Bretz lattice rule of :meth:`cdf` above two dimensions (``MVNUNI``).

        Reads as a :class:`pynumerics.sampling.MersenneTwister` when it is one. Accepts a generator or
        an integer seed (a new Numerics ``MersenneTwister``).
        """
        return _wrap_generator(self.net.MVNUNI)

    @mvnuni.setter
    def mvnuni(self, value: Any) -> None:
        self.net.MVNUNI = _seeded_generator(value)

    @property
    def max_evaluations(self) -> int:
        """The maximum number of lattice-rule evaluations of the CDF above two dimensions; default 25,000."""
        return int(self.net.MaxEvaluations)

    @max_evaluations.setter
    def max_evaluations(self, value: int) -> None:
        self.net.MaxEvaluations = int(value)

    @property
    def absolute_error(self) -> float:
        """The absolute error tolerance of the CDF above two dimensions; default 1e-5."""
        return float(self.net.AbsoluteError)

    @absolute_error.setter
    def absolute_error(self, value: float) -> None:
        self.net.AbsoluteError = float(value)

    @property
    def relative_error(self) -> float:
        """The relative error tolerance of the CDF above two dimensions; default 0."""
        return float(self.net.RelativeError)

    @relative_error.setter
    def relative_error(self, value: float) -> None:
        self.net.RelativeError = float(value)

    @property
    def degrees_of_freedom(self) -> float:
        """The degrees of freedom ν."""
        return float(self.net.DegreesOfFreedom)

    @property
    def location(self) -> np.ndarray:
        """The location vector μ."""
        return to_numpy_doubles(self.net.Location)

    @property
    def mean(self) -> np.ndarray:
        """The mean vector (μ for ν > 1)."""
        return to_numpy_doubles(self.net.Mean)

    @property
    def median(self) -> np.ndarray:
        """The median vector μ."""
        return to_numpy_doubles(self.net.Median)

    @property
    def mode(self) -> np.ndarray:
        """The mode vector μ."""
        return to_numpy_doubles(self.net.Mode)

    @property
    def variance(self) -> np.ndarray:
        """The marginal variances ν/(ν − 2)·Σᵢᵢ for ν > 2."""
        return to_numpy_doubles(self.net.Variance)

    @property
    def standard_deviation(self) -> np.ndarray:
        """The marginal standard deviations."""
        return to_numpy_doubles(self.net.StandardDeviation)

    @property
    def scale_matrix(self) -> np.ndarray:
        """The scale matrix Σ."""
        return to_numpy_matrix(self.net.ScaleMatrix)

    @property
    def covariance(self) -> np.ndarray:
        """The covariance matrix ν/(ν − 2)·Σ for ν > 2."""
        return to_numpy_matrix(self.net.Covariance)

    @property
    def is_positive_definite(self) -> bool:
        """Whether Σ is positive definite."""
        return bool(self.net.IsPositiveDefinite)

    def set_parameters(self, degrees_of_freedom: float, location: Any, scale_matrix: Any) -> None:
        """Replace the parameters (``SetParameters``).

        Parameters
        ----------
        degrees_of_freedom : float
            Degrees of freedom ν.
        location : array_like
            Location vector μ.
        scale_matrix : array_like
            Scale matrix Σ.
        """
        self.net.SetParameters(float(degrees_of_freedom), to_net_doubles(location), to_net_matrix(scale_matrix))

    def validate_parameters(self, degrees_of_freedom: float, location: Any, scale_matrix: Any) -> str | None:
        """Return why a parameter set is invalid, or None when it is valid.

        Parameters
        ----------
        degrees_of_freedom : float
            Degrees of freedom ν.
        location : array_like
            Location vector μ.
        scale_matrix : array_like
            Scale matrix Σ.

        Returns
        -------
        str or None
            The Numerics message of ``ValidateParameters(..., false)``, or None.
        """
        return _message(self.net.ValidateParameters(float(degrees_of_freedom), to_net_doubles(location),
                                                    to_net_matrix(scale_matrix), False))

    def mahalanobis(self, x: Any) -> float:
        """Return the squared Mahalanobis distance (x − μ)ᵀ Σ⁻¹ (x − μ) (``Mahalanobis``).

        Parameters
        ----------
        x : array_like
            Point of length d.

        Returns
        -------
        float
            The squared distance.
        """
        return float(self.net.Mahalanobis(to_net_doubles(x)))

    def inverse_cdf(self, probabilities: Any) -> np.ndarray:
        """Map d + 1 probabilities to a point (``InverseCDF(double[])``).

        Parameters
        ----------
        probabilities : array_like
            Length d + 1: the first d map to standard normal deviates, the last to the χ²(ν) mixing
            variate W, and the point is μ + L·z·√(ν/W).

        Returns
        -------
        numpy.ndarray
            The point of length d.
        """
        return to_numpy_doubles(self.net.InverseCDF(to_net_doubles(probabilities)))

    def generate_random_values(self, sample_size: int, seed: Any = -1) -> np.ndarray:
        """Draw a Monte Carlo sample (``GenerateRandomValues``).

        Parameters
        ----------
        sample_size : int
            Number of points.
        seed : int or MersenneTwister, optional
            Seed (zero or negative seeds from the clock, not reproducible) or a generator.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(sample_size, d)``.
        """
        return _draw(self.net.GenerateRandomValues, sample_size, seed)

    def latin_hypercube_random_values(self, sample_size: int, seed: int) -> np.ndarray:
        """Draw a Latin hypercube sample of the normal part (``LatinHypercubeRandomValues``).

        Parameters
        ----------
        sample_size : int
            Number of points.
        seed : int
            Seed; the χ² mixing variates use ``seed + 1``.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(sample_size, d)``.
        """
        return to_numpy_matrix(self.net.LatinHypercubeRandomValues(int(sample_size), int(seed)))

    def stratified_random_values(self, stratification_bins: Any, seed: int) -> np.ndarray:
        """Draw a sample stratified in the first variable (``StratifiedRandomValues``).

        Parameters
        ----------
        stratification_bins : pynumerics.sampling.StratificationBins
            Probability bins; one point is drawn per bin at its midpoint.
        seed : int
            Seed of the other variables; the χ² mixing variates use ``seed + 1``.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(stratification_bins), d)``.
        """
        return to_numpy_matrix(self.net.StratifiedRandomValues(_bins(stratification_bins), int(seed)))


class Multinomial(MultivariateDistribution):
    """Multinomial distribution of category counts (``Numerics.Distributions.Multinomial``).

    Parameters
    ----------
    number_of_trials : int
        Number of trials N ≥ 1.
    probabilities : array_like
        Category probabilities p of length k ≥ 2, each in [0, 1], summing to one.

    Attributes
    ----------
    net : Numerics.Distributions.Multinomial
        The wrapped .NET object.

    Notes
    -----
    Points are count vectors of length k; ``pdf`` is the probability mass function.
    """

    def __init__(self, number_of_trials: int, probabilities: Any) -> None:
        MultivariateDistribution.__init__(self, _nd().Multinomial(int(number_of_trials), to_net_doubles(probabilities)))

    @property
    def number_of_trials(self) -> int:
        """The number of trials N."""
        return int(self.net.NumberOfTrials)

    @property
    def probabilities(self) -> np.ndarray:
        """The category probabilities p."""
        return to_numpy_doubles(self.net.Probabilities)

    @property
    def mean(self) -> np.ndarray:
        """The mean counts N·pᵢ."""
        return to_numpy_doubles(self.net.Mean)

    @property
    def variance(self) -> np.ndarray:
        """The count variances N·pᵢ(1 − pᵢ)."""
        return to_numpy_doubles(self.net.Variance)

    def covariance(self, i: int, j: int) -> float:
        """Return one covariance entry (``Covariance(int, int)``).

        Parameters
        ----------
        i : int
            Row category index.
        j : int
            Column category index.

        Returns
        -------
        float
            Cov(Xᵢ, Xⱼ): N·pᵢ(1 − pᵢ) on the diagonal, −N·pᵢ·pⱼ off it.
        """
        return float(self.net.Covariance(int(i), int(j)))

    def log_pmf(self, x: Any) -> Any:
        """Return the log probability mass (``LogPMF``).

        Parameters
        ----------
        x : array_like
            One count vector of length k, or an array of shape (n, k).

        Returns
        -------
        float or numpy.ndarray
            The log mass; negative infinity off the support.
        """
        return self._point("LogPMF", x)

    def generate_random_values(self, sample_size: int, seed: int = -1) -> np.ndarray:
        """Draw count vectors by sequential binomial sampling (``GenerateRandomValues``).

        Parameters
        ----------
        sample_size : int
            Number of draws.
        seed : int, optional
            Seed; zero or negative seeds from the clock (not reproducible).

        Returns
        -------
        numpy.ndarray
            Array of shape ``(sample_size, k)``; each row sums to N.
        """
        return to_numpy_matrix(self.net.GenerateRandomValues(int(sample_size), int(seed)))

    @staticmethod
    def sample(weights: Any, random: Any) -> int:
        """Draw one category index with probability proportional to its weight (``Multinomial.Sample``).

        Parameters
        ----------
        weights : array_like
            Non-negative weights, at least one positive; normalization is not required.
        random : MersenneTwister
            The generator (a :class:`pynumerics.sampling.MersenneTwister` or a .NET ``System.Random``);
            one draw is consumed.

        Returns
        -------
        int
            The zero-based category index; a draw on a boundary selects the lower category.
        """
        return int(_nd().Multinomial.Sample(to_net_doubles(weights), _generator(random)))

    @staticmethod
    def categorical_indices(cumulative_weights: Any, uniforms: Any, inclusive: bool = False) -> Any:
        """Map uniforms to category indices through a cumulative weight table (``CategoricalIndices``).

        Parameters
        ----------
        cumulative_weights : array_like
            Non-decreasing finite cumulative weights c of length n ≥ 1; normalization is not required.
        uniforms : array_like
            Variates on the scale of ``cumulative_weights``, any shape, no NaN.
        inclusive : bool, optional
            False (the default) selects min{k : c[k] > u}; True selects min{k : c[k] ≥ u}, the
            inverse-transform count min{k : F(k) ≥ u} of a discrete CDF table F. Both clip to n - 1.

        Returns
        -------
        int or numpy.ndarray
            One zero-based index per variate (int32, the shape of ``uniforms``), equal to
            ``numpy.clip(numpy.searchsorted(c, u, side="right"), 0, n - 1)``, or with ``side="left"``
            when ``inclusive`` is True.

        Notes
        -----
        The two forms differ only for a variate equal to a table entry. The whole array crosses into
        Numerics once, and the indices cross back once.
        """
        member = cached_overload("Multinomial.CategoricalIndices(double[], double[], bool)", _categorical_indices)
        array = np.asarray(uniforms, dtype=np.float64)
        indices = to_numpy_ints(member(to_net_doubles(cumulative_weights), _flat_doubles(array), bool(inclusive)))
        return int(indices[0]) if array.ndim == 0 else indices.reshape(array.shape)


class BivariateEmpirical(MultivariateDistribution):
    """Bivariate distribution interpolated from a table of joint non-exceedance probabilities.

    Wraps ``Numerics.Distributions.BivariateEmpirical``: the CDF is bilinear interpolation of
    ``probability_values`` on the grid ``x1_values`` by ``x2_values``, in the chosen transforms.

    Parameters
    ----------
    x1_values : array_like, optional
        Ascending grid of the first variable. Omit all three tables for an empty distribution
        to fill later with :meth:`set_parameters`.
    x2_values : array_like, optional
        Ascending grid of the second variable.
    probability_values : array_like, optional
        Joint probabilities of shape ``(len(x1_values), len(x2_values))``.
    x1_transform : Transform or str, optional
        Interpolation transform of the first variable; default none.
    x2_transform : Transform or str, optional
        Interpolation transform of the second variable; default none.
    probability_transform : Transform or str, optional
        Interpolation transform of the probabilities; default none.

    Attributes
    ----------
    net : Numerics.Distributions.BivariateEmpirical
        The wrapped .NET object.
    """

    def __init__(self, x1_values: Any = None, x2_values: Any = None, probability_values: Any = None,
                 x1_transform: Transform | str = Transform.NONE, x2_transform: Transform | str = Transform.NONE,
                 probability_transform: Transform | str = Transform.NONE) -> None:
        nd = _nd()
        kind = namespace("Numerics.Data").Transform
        transforms = (net_enum(kind, x1_transform), net_enum(kind, x2_transform), net_enum(kind, probability_transform))
        if x1_values is None and x2_values is None and probability_values is None:
            net = nd.BivariateEmpirical(*transforms)
        else:
            net = nd.BivariateEmpirical(to_net_doubles(x1_values), to_net_doubles(x2_values),
                                        to_net_matrix(probability_values), *transforms)
        MultivariateDistribution.__init__(self, net)

    def _bind(self) -> None:
        """Bind the ``(double, double)`` overloads of ``PDF`` and ``CDF`` once."""
        import System  # noqa: PLC0415 - the runtime is running (self.net exists)

        self._pair = {name: getattr(self.net, name).Overloads[System.Double, System.Double] for name in ("PDF", "CDF")}

    def _pairs(self, name: str, x1: Any, x2: Any) -> Any:
        """Evaluate a point member at a point, or at broadcast (x1, x2) pairs with the scalar overload."""
        if x2 is None:
            return self._point(name, x1)
        member = self._pair[name]
        a, b = np.broadcast_arrays(np.asarray(x1, dtype=np.float64), np.asarray(x2, dtype=np.float64))
        if a.ndim == 0:
            return float(member(float(a), float(b)))
        out = np.array([float(member(float(p), float(q))) for p, q in zip(a.ravel(), b.ravel(), strict=True)])
        return out.reshape(a.shape)

    @property
    def x1_values(self) -> np.ndarray:
        """The grid of the first variable."""
        return to_numpy_doubles(self.net.X1Values)

    @property
    def x2_values(self) -> np.ndarray:
        """The grid of the second variable."""
        return to_numpy_doubles(self.net.X2Values)

    @property
    def probability_values(self) -> np.ndarray:
        """The table of joint non-exceedance probabilities."""
        return to_numpy_matrix(self.net.ProbabilityValues)

    @property
    def x1_transform(self) -> Transform:
        """The interpolation transform of the first variable; settable."""
        return Transform(str(self.net.X1Transform))

    @x1_transform.setter
    def x1_transform(self, value: Transform | str) -> None:
        self.net.X1Transform = net_enum(namespace("Numerics.Data").Transform, value)

    @property
    def x2_transform(self) -> Transform:
        """The interpolation transform of the second variable; settable."""
        return Transform(str(self.net.X2Transform))

    @x2_transform.setter
    def x2_transform(self, value: Transform | str) -> None:
        self.net.X2Transform = net_enum(namespace("Numerics.Data").Transform, value)

    @property
    def probability_transform(self) -> Transform:
        """The interpolation transform of the probabilities; settable."""
        return Transform(str(self.net.ProbabilityTransform))

    @probability_transform.setter
    def probability_transform(self, value: Transform | str) -> None:
        self.net.ProbabilityTransform = net_enum(namespace("Numerics.Data").Transform, value)

    def set_parameters(self, x1_values: Any, x2_values: Any, probability_values: Any) -> None:
        """Replace the grids and the probability table (``SetParameters``).

        Parameters
        ----------
        x1_values : array_like
            Grid of the first variable.
        x2_values : array_like
            Grid of the second variable.
        probability_values : array_like
            Joint probabilities of shape ``(len(x1_values), len(x2_values))``.
        """
        self.net.SetParameters(to_net_doubles(x1_values), to_net_doubles(x2_values), to_net_matrix(probability_values))

    def validate_parameters(self, x1_values: Any, x2_values: Any, probability_values: Any) -> str | None:
        """Return why a grid and table are invalid, or None when they are valid.

        Parameters
        ----------
        x1_values : array_like
            Grid of the first variable.
        x2_values : array_like
            Grid of the second variable.
        probability_values : array_like
            Joint probabilities.

        Returns
        -------
        str or None
            The Numerics message of ``ValidateParameters(..., false)``, or None.
        """
        return _message(self.net.ValidateParameters(to_net_doubles(x1_values), to_net_doubles(x2_values),
                                                    to_net_matrix(probability_values), False))

    def pdf(self, x1: Any, x2: Any = None) -> Any:
        """Return the density, which Numerics does not define for the interpolated surface (NaN).

        Parameters
        ----------
        x1 : float or array_like
            The first variable, or one point (x1, x2) when ``x2`` is omitted.
        x2 : float or array_like, optional
            The second variable, broadcast with ``x1``.

        Returns
        -------
        float or numpy.ndarray
            NaN for every point.
        """
        return self._pairs("PDF", x1, x2)

    def cdf(self, x1: Any, x2: Any = None) -> Any:
        """Return the interpolated joint CDF P(X1 ≤ x1, X2 ≤ x2) (``CDF(double, double)``).

        Parameters
        ----------
        x1 : float or array_like
            The first variable, or one point (x1, x2) when ``x2`` is omitted.
        x2 : float or array_like, optional
            The second variable, broadcast with ``x1``.

        Returns
        -------
        float or numpy.ndarray
            Joint non-exceedance probabilities.
        """
        return self._pairs("CDF", x1, x2)


class Dirichlet(MultivariateDistribution):
    """Dirichlet distribution on the probability simplex (``Numerics.Distributions.Dirichlet``).

    Parameters
    ----------
    alpha : array_like
        Concentration parameters α of length K ≥ 2, each positive and finite. Numerics refuses
        anything else with ``ArgumentOutOfRangeException`` (``ValueError`` inside
        :func:`pynumerics.errors.net_errors`).

    Attributes
    ----------
    net : Numerics.Distributions.Dirichlet
        The wrapped .NET object.

    Notes
    -----
    Every value comes from Numerics; nothing is computed in Python. Draws use the gamma
    representation inside ``Dirichlet.GenerateRandomValues``: K independent
    Gamma(αᵢ, 1) inverse-CDF variates from one seeded Mersenne Twister, normalized to sum
    to one. The sample crosses back as one ``double[,]`` block copy, so a draw is bit for bit
    the ``double[,]`` the .NET call returns. The density is zero (log-density negative
    infinity) off the simplex.
    """

    def __init__(self, alpha: Any) -> None:
        MultivariateDistribution.__init__(self, _nd().Dirichlet(to_net_doubles(alpha)))

    @classmethod
    def symmetric(cls, dimension: int, alpha: float) -> Dirichlet:
        """Build the symmetric Dirichlet with every αᵢ = ``alpha`` (``Dirichlet(int, double)``).

        Parameters
        ----------
        dimension : int
            Number of components K ≥ 2.
        alpha : float
            Common positive concentration parameter.

        Returns
        -------
        Dirichlet
            The wrapped distribution.
        """
        return cls.from_net(_nd().Dirichlet(int(dimension), float(alpha)))

    def _bind(self) -> None:
        """Bind ``GenerateRandomValues(int, int)`` once, so calls skip overload resolution."""
        import System  # noqa: PLC0415 - the runtime is started by the constructor

        self._generate = self.net.GenerateRandomValues.Overloads[System.Int32, System.Int32]

    @property
    def alpha(self) -> np.ndarray:
        """The concentration parameters α (a copy)."""
        return to_numpy_doubles(self.net.Alpha)

    @property
    def alpha_sum(self) -> float:
        """The concentration sum α₀ = Σαᵢ."""
        return float(self.net.AlphaSum)

    @property
    def mean(self) -> np.ndarray:
        """The mean vector αᵢ/α₀."""
        return to_numpy_doubles(self.net.Mean)

    @property
    def variance(self) -> np.ndarray:
        """The marginal variances αᵢ(α₀ − αᵢ)/(α₀²(α₀ + 1))."""
        return to_numpy_doubles(self.net.Variance)

    @property
    def mode(self) -> np.ndarray:
        """The mode (αᵢ − 1)/(α₀ − K); Numerics refuses it unless every αᵢ > 1."""
        return to_numpy_doubles(self.net.Mode)

    def covariance(self, i: int, j: int) -> float:
        """Return one covariance entry (``Covariance(int, int)``).

        Parameters
        ----------
        i : int
            Row component index.
        j : int
            Column component index.

        Returns
        -------
        float
            Cov(Xᵢ, Xⱼ).
        """
        return float(self.net.Covariance(int(i), int(j)))

    def covariance_matrix(self) -> np.ndarray:
        """Return the full covariance matrix (``CovarianceMatrix()``).

        Returns
        -------
        numpy.ndarray
            Array of shape ``(K, K)``.
        """
        return to_numpy_matrix(self.net.CovarianceMatrix())

    def generate_random_values(self, sample_size: int, seed: int = -1) -> np.ndarray:
        """Draw a Monte Carlo sample (``GenerateRandomValues(int, int)``).

        Parameters
        ----------
        sample_size : int
            Number of points; zero returns an empty ``(0, K)`` array.
        seed : int, optional
            Mersenne Twister seed; zero or negative seeds from the clock (not reproducible).

        Returns
        -------
        numpy.ndarray
            Array of shape ``(sample_size, K)``; each row lies on the simplex.
        """
        return to_numpy_matrix(self._generate(int(sample_size), int(seed)))

    @staticmethod
    def log_multivariate_beta(alpha: Any) -> float:
        """Return log B(α) = Σ log Γ(αᵢ) − log Γ(Σαᵢ) (``Dirichlet.LogMultivariateBeta``).

        Parameters
        ----------
        alpha : array_like
            Concentration parameters.

        Returns
        -------
        float
            The log of the multivariate beta function.
        """
        return float(_nd().Dirichlet.LogMultivariateBeta(to_net_doubles(alpha)))


_CLASSES: dict[str, type[MultivariateDistribution]] = {
    "MultivariateNormal": MultivariateNormal,
    "MultivariateStudentT": MultivariateStudentT,
    "Multinomial": Multinomial,
    "BivariateEmpirical": BivariateEmpirical,
    "Dirichlet": Dirichlet,
}


def wrap(net: Any) -> MultivariateDistribution:
    """Wrap a .NET multivariate distribution in the matching Python class.

    Parameters
    ----------
    net : Numerics.Distributions.MultivariateDistribution or IMultivariateDistribution
        The .NET object.

    Returns
    -------
    MultivariateDistribution
        The specific wrapper class when one exists, otherwise the generic wrapper; the .NET object
        is shared, not copied.
    """
    concrete = unwrap(net)
    cls = _CLASSES.get(str(concrete.GetType().Name), MultivariateDistribution)
    return cls.from_net(concrete)


__all__ = [
    "BivariateEmpirical", "Dirichlet", "Multinomial", "MultivariateDistribution", "MultivariateDistributionType",
    "MultivariateNormal", "MultivariateStudentT", "wrap",
]

translate_public(globals())
