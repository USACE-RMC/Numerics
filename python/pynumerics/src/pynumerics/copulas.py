"""Bivariate copulas from ``Numerics.Distributions.Copulas``.

Every class wraps one .NET copula (``.net``). Paired (u, v) arrays cross once into the
``(IList<double>, IList<double>)`` overloads of ``BivariateCopula`` where Numerics has them, and
are evaluated pair by pair with the scalar member otherwise. Parameter estimation
(``BivariateCopulaEstimation.Estimate``), the factory (``CopulaFactory``) and the XML round trip
(``ToXElement`` and ``CopulaFactory.CreateCopula(XElement)``) are Numerics calls.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from pynumerics._util import bind, loop, namespace, net_enum, unwrap
from pynumerics.distributions import FrozenDistribution, UnivariateDistribution
from pynumerics.distributions import wrap as wrap_univariate
from pynumerics.enums import NetEnum
from pynumerics.errors import translate_public
from pynumerics.marshal import to_net_doubles, to_numpy_doubles, to_numpy_matrix
from pynumerics.multivariate import _draw, _message


def _nc() -> Any:
    """Return the ``Numerics.Distributions.Copulas`` namespace after starting the runtime."""
    return namespace("Numerics.Distributions.Copulas")


class CopulaType(NetEnum):
    """Copula families (``Numerics.Distributions.Copulas.CopulaType``)."""

    ALI_MIKHAIL_HAQ = "AliMikhailHaq"
    CLAYTON = "Clayton"
    FRANK = "Frank"
    GUMBEL = "Gumbel"
    JOE = "Joe"
    NORMAL = "Normal"
    STUDENT_T = "StudentT"
    INDEPENDENCE = "Independence"


class CopulaEstimationMethod(NetEnum):
    """Copula estimation methods (``Numerics.Distributions.Copulas.CopulaEstimationMethod``)."""

    FULL_LIKELIHOOD = "FullLikelihood"
    PSEUDO_LIKELIHOOD = "PseudoLikelihood"
    INFERENCE_FROM_MARGINS = "InferenceFromMargins"


def _pairwise(u: Any, v: Any, copula: Any, name: str) -> Any:
    """Evaluate a (u, v) copula member over broadcast arrays.

    A scalar pair calls the ``(double, double)`` overload. Arrays cross once into the
    ``(IList<double>, IList<double>)`` overload of ``BivariateCopula`` when the member has one, or
    pair by pair otherwise; every value is the scalar member's. Both overloads are bound once per
    copula object.
    """
    import System  # noqa: PLC0415 - the runtime is running (the copula exists)
    from System.Collections.Generic import IList  # noqa: PLC0415

    cache = copula.__dict__.setdefault("_bound", {})
    entry = cache.get(name)
    if entry is None or entry[0] is not copula.net:
        group = getattr(copula.net, name)
        scalar = group.Overloads[System.Double, System.Double]
        try:
            vector = group.Overloads[IList[System.Double], IList[System.Double]]
        except (KeyError, TypeError):
            vector = None
        entry = cache[name] = (copula.net, scalar, vector)
    _, scalar, vector = entry
    uu, vv = np.broadcast_arrays(np.asarray(u, dtype=np.float64), np.asarray(v, dtype=np.float64))
    if uu.ndim == 0:
        return float(scalar(float(uu), float(vv)))
    if vector is not None:
        return to_numpy_doubles(vector(to_net_doubles(uu.ravel()), to_net_doubles(vv.ravel()))).reshape(uu.shape)
    out = np.empty(uu.size, dtype=np.float64)
    for index, (a, b) in enumerate(zip(uu.ravel(), vv.ravel(), strict=True)):
        out[index] = scalar(float(a), float(b))
    return out.reshape(uu.shape)


def _marginal_net(distribution: UnivariateDistribution | None) -> Any:
    """Return the .NET marginal to attach: None, the distribution's object, or a clone of a frozen one.

    A frozen distribution is read-only, and estimation with the full likelihood refits attached
    marginals in place, so a frozen marginal is attached as a clone.
    """
    if distribution is None:
        return None
    if isinstance(distribution, FrozenDistribution):
        return distribution.net.Clone()
    return distribution.net


class BivariateCopula:
    """A wrapped ``Numerics.Distributions.Copulas.BivariateCopula``.

    Parameters
    ----------
    net : Numerics.Distributions.Copulas.BivariateCopula
        The .NET copula.

    Attributes
    ----------
    net : Numerics.Distributions.Copulas.BivariateCopula
        The wrapped .NET object.
    """

    def __init__(self, net: Any) -> None:
        self.net = unwrap(net)

    def _attach(self, marginal_x: UnivariateDistribution | None, marginal_y: UnivariateDistribution | None) -> None:
        """Attach the marginals given to a constructor, leaving unset ones untouched."""
        if marginal_x is not None or marginal_y is not None:
            self.set_marginals(marginal_x, marginal_y)

    @property
    def type(self) -> CopulaType:
        """The copula family (``Type``)."""
        return CopulaType(str(self.net.Type))

    @property
    def display_name(self) -> str:
        """The Numerics display name."""
        return str(self.net.DisplayName)

    @property
    def short_display_name(self) -> str:
        """The Numerics short display name."""
        return str(self.net.ShortDisplayName)

    @property
    def theta(self) -> float:
        """The dependence parameter θ (ρ for the normal and Student-t copulas); settable."""
        return float(self.net.Theta)

    @theta.setter
    def theta(self, value: float) -> None:
        self.net.Theta = float(value)

    @property
    def theta_minimum(self) -> float:
        """The smallest admissible θ."""
        return float(self.net.ThetaMinimum)

    @property
    def theta_maximum(self) -> float:
        """The largest admissible θ."""
        return float(self.net.ThetaMaximum)

    @property
    def parameter_to_string(self) -> tuple[tuple[str, ...], ...]:
        """The parameter names and values as Numerics formats them, one (name, value) row per parameter."""
        table = self.net.ParameterToString
        return tuple(tuple(str(table[i, j]) for j in range(table.GetLength(1))) for i in range(table.GetLength(0)))

    @property
    def parameter_name_short_form(self) -> str:
        """The short parameter name (for example θ or ρ)."""
        return str(self.net.ParameterNameShortForm)

    @property
    def parameter_property_name(self) -> str:
        """The name of the parameter property (``GetParameterPropertyName``)."""
        return str(self.net.GetParameterPropertyName)

    @property
    def parameters_valid(self) -> bool:
        """Whether the current parameters are valid (``ParametersValid``)."""
        return bool(self.net.ParametersValid)

    @property
    def upper_tail_dependence(self) -> float:
        """The upper tail dependence coefficient λ_U."""
        return float(self.net.UpperTailDependence)

    @property
    def lower_tail_dependence(self) -> float:
        """The lower tail dependence coefficient λ_L."""
        return float(self.net.LowerTailDependence)

    @property
    def number_of_copula_parameters(self) -> int:
        """The number of copula parameters (0, 1, or 2 for the Student-t copula)."""
        return int(self.net.NumberOfCopulaParameters)

    @property
    def copula_parameters(self) -> np.ndarray:
        """The copula parameters in Numerics order (``GetCopulaParameters``)."""
        return to_numpy_doubles(self.net.GetCopulaParameters)

    def set_copula_parameters(self, parameters: Any) -> None:
        """Set the copula parameters in Numerics order (``SetCopulaParameters``).

        Parameters
        ----------
        parameters : array_like
            One value per copula parameter.
        """
        self.net.SetCopulaParameters(to_net_doubles(parameters))

    def validate_parameter(self, parameter: float) -> str | None:
        """Return why a dependence parameter is invalid, or None when it is valid.

        Parameters
        ----------
        parameter : float
            Candidate θ.

        Returns
        -------
        str or None
            The Numerics message of ``ValidateParameter(parameter, false)``, or None.
        """
        return _message(self.net.ValidateParameter(float(parameter), False))

    def parameter_constraints(self, x: Any, y: Any) -> np.ndarray:
        """Return the search bounds of each copula parameter for a sample (``ParameterConstraints``).

        Parameters
        ----------
        x : array_like
            Sample of the first variable.
        y : array_like
            Sample of the second variable.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(number_of_copula_parameters, 2)`` of lower and upper bounds.
        """
        return to_numpy_matrix(self.net.ParameterConstraints(to_net_doubles(x), to_net_doubles(y)))

    @property
    def marginal_distribution_x(self) -> UnivariateDistribution | None:
        """The marginal of the first variable, or None; settable (a frozen marginal is attached as a clone)."""
        net = self.net.MarginalDistributionX
        return None if net is None else wrap_univariate(net)

    @marginal_distribution_x.setter
    def marginal_distribution_x(self, value: UnivariateDistribution | None) -> None:
        self.net.MarginalDistributionX = _marginal_net(value)

    @property
    def marginal_distribution_y(self) -> UnivariateDistribution | None:
        """The marginal of the second variable, or None; settable (a frozen marginal is attached as a clone)."""
        net = self.net.MarginalDistributionY
        return None if net is None else wrap_univariate(net)

    @marginal_distribution_y.setter
    def marginal_distribution_y(self, value: UnivariateDistribution | None) -> None:
        self.net.MarginalDistributionY = _marginal_net(value)

    def set_marginals(self, x: UnivariateDistribution | None, y: UnivariateDistribution | None) -> None:
        """Set the marginal distributions used by sampling and the likelihoods.

        Parameters
        ----------
        x : UnivariateDistribution or None
            Marginal of the first variable.
        y : UnivariateDistribution or None
            Marginal of the second variable.

        Notes
        -----
        The copula holds the distribution's own .NET object, so estimation with the full
        likelihood refits it in place; a frozen distribution is attached as a clone.
        """
        self.net.MarginalDistributionX = _marginal_net(x)
        self.net.MarginalDistributionY = _marginal_net(y)

    def pdf(self, u: Any, v: Any) -> Any:
        """Return the copula density c(u, v).

        Parameters
        ----------
        u : float or array_like
            First uniform margin.
        v : float or array_like
            Second uniform margin, broadcast with ``u``.

        Returns
        -------
        float or numpy.ndarray
            Densities.
        """
        return _pairwise(u, v, self, "PDF")

    def log_pdf(self, u: Any, v: Any) -> Any:
        """Return the log copula density log c(u, v) (``LogPDF``).

        Parameters
        ----------
        u : float or array_like
            First uniform margin.
        v : float or array_like
            Second uniform margin, broadcast with ``u``.

        Returns
        -------
        float or numpy.ndarray
            Log-densities.
        """
        return _pairwise(u, v, self, "LogPDF")

    def cdf(self, u: Any, v: Any) -> Any:
        """Return the copula C(u, v).

        Parameters
        ----------
        u : float or array_like
            First uniform margin.
        v : float or array_like
            Second uniform margin, broadcast with ``u``.

        Returns
        -------
        float or numpy.ndarray
            Joint non-exceedance probabilities.
        """
        return _pairwise(u, v, self, "CDF")

    def or_joint_exceedance_probability(self, u: Any, v: Any) -> Any:
        """Return P(U > u or V > v) = 1 − C(u, v) (``ORJointExceedanceProbability``).

        Parameters
        ----------
        u : float or array_like
            First uniform margin.
        v : float or array_like
            Second uniform margin, broadcast with ``u``.

        Returns
        -------
        float or numpy.ndarray
            Probabilities that either variable exceeds its threshold.
        """
        return _pairwise(u, v, self, "ORJointExceedanceProbability")

    def and_joint_exceedance_probability(self, u: Any, v: Any) -> Any:
        """Return P(U > u and V > v) = 1 − u − v + C(u, v) (``ANDJointExceedanceProbability``).

        Parameters
        ----------
        u : float or array_like
            First uniform margin.
        v : float or array_like
            Second uniform margin, broadcast with ``u``.

        Returns
        -------
        float or numpy.ndarray
            Probabilities that both variables exceed their thresholds.
        """
        return _pairwise(u, v, self, "ANDJointExceedanceProbability")

    def conditional_cdf(self, u: Any, v: Any) -> Any:
        """Return the conditional CDF ``ConditionalCDF(u, v)``.

        Parameters
        ----------
        u : float or array_like
            First uniform margin.
        v : float or array_like
            Second uniform margin, broadcast with ``u``.

        Returns
        -------
        float or numpy.ndarray
            Conditional probabilities.
        """
        return _pairwise(u, v, self, "ConditionalCDF")

    def inverse_conditional_cdf(self, u: Any, t: Any) -> Any:
        """Return the inverse conditional CDF ``InverseConditionalCDF(u, t)``.

        Parameters
        ----------
        u : float or array_like
            Conditioning margin.
        t : float or array_like
            Conditional probability.

        Returns
        -------
        float or numpy.ndarray
            The second margin v.
        """
        return _pairwise(u, t, self, "InverseConditionalCDF")

    def inverse_cdf(self, u: float, v: float) -> np.ndarray:
        """Map two independent uniforms to a dependent pair (``InverseCDF(u, v)``).

        Parameters
        ----------
        u : float
            First independent uniform.
        v : float
            Second independent uniform.

        Returns
        -------
        numpy.ndarray
            The dependent uniform pair.
        """
        return to_numpy_doubles(self.net.InverseCDF(float(u), float(v)))

    def generate_random_values(self, sample_size: int, seed: Any = -1) -> np.ndarray:
        """Draw a sample (``GenerateRandomValues``; marginal values when both marginals are set).

        Parameters
        ----------
        sample_size : int
            Number of pairs.
        seed : int or MersenneTwister, optional
            Seed (zero or negative seeds from the clock, not reproducible) or a generator.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(sample_size, 2)``.
        """
        return _draw(self.net.GenerateRandomValues, sample_size, seed)

    def pseudo_log_likelihood(self, x: Any, y: Any) -> float:
        """Return the pseudo log-likelihood Σ log c(xᵢ, yᵢ) (``PseudoLogLikelihood``).

        Parameters
        ----------
        x : array_like
            Plotting positions of the first variable.
        y : array_like
            Plotting positions of the second variable.

        Returns
        -------
        float
            The log-likelihood; negative infinity when it is not finite.
        """
        return float(self.net.PseudoLogLikelihood(to_net_doubles(x), to_net_doubles(y)))

    def ifm_log_likelihood(self, x: Any, y: Any) -> float:
        """Return the inference-from-margins log-likelihood Σ log c(F(xᵢ), G(yᵢ)) (``IFMLogLikelihood``).

        Parameters
        ----------
        x : array_like
            Sample of the first variable.
        y : array_like
            Sample of the second variable.

        Returns
        -------
        float
            The log-likelihood with the attached marginals; negative infinity when it is not finite.
        """
        return float(self.net.IFMLogLikelihood(to_net_doubles(x), to_net_doubles(y)))

    def log_likelihood(self, x: Any, y: Any) -> float:
        """Return the full log-likelihood of copula and marginals (``LogLikelihood``).

        Parameters
        ----------
        x : array_like
            Sample of the first variable.
        y : array_like
            Sample of the second variable.

        Returns
        -------
        float
            Σ [log c(F(xᵢ), G(yᵢ)) + log f(xᵢ) + log g(yᵢ)]; negative infinity when it is not finite.
        """
        return float(self.net.LogLikelihood(to_net_doubles(x), to_net_doubles(y)))

    def estimate(self, x: Any, y: Any,
                 method: CopulaEstimationMethod | str = CopulaEstimationMethod.PSEUDO_LIKELIHOOD) -> None:
        """Fit the copula in place (``BivariateCopulaEstimation.Estimate``).

        Parameters
        ----------
        x : array_like
            Sample of the first variable (plotting positions for the pseudo likelihood).
        y : array_like
            Sample of the second variable (plotting positions for the pseudo likelihood).
        method : CopulaEstimationMethod or str, optional
            Default the pseudo likelihood. Inference from margins needs both marginals attached;
            the full likelihood also refits them, so they must support maximum likelihood.

        Notes
        -----
        One-parameter copulas are fitted with Brent's method over the parameter constraints and
        the Student-t copula with Nelder-Mead; the full likelihood fits the copula and both
        marginals jointly with Nelder-Mead, restarted once. A zero-parameter copula is unchanged.
        """
        nc = _nc()
        self.net = unwrap(nc.BivariateCopulaEstimation.Estimate(
            self.net, to_net_doubles(x), to_net_doubles(y), net_enum(nc.CopulaEstimationMethod, method)))

    def clone(self) -> BivariateCopula:
        """Return a deep copy, marginals included (``Clone``).

        Returns
        -------
        BivariateCopula
            A wrapper of the cloned .NET copula, of the matching Python class.
        """
        return wrap(self.net.Clone())

    def to_xml(self, include_marginals: bool = False) -> str:
        """Return the copula's ``ToXElement`` as XML text.

        Parameters
        ----------
        include_marginals : bool, optional
            Also write the attached marginals (``ToXElement(true)``).

        Returns
        -------
        str
            The element as .NET writes it; :func:`copula_from_xml` reads it back.
        """
        if not include_marginals:
            return str(self.net.ToXElement().ToString())
        import System  # noqa: PLC0415 - the runtime is running (self.net exists)

        return str(self.net.ToXElement.Overloads[System.Boolean](True).ToString())


class ArchimedeanCopula(BivariateCopula):
    """A wrapped ``Numerics.Distributions.Copulas.ArchimedeanCopula`` with its generator φ.

    Parameters
    ----------
    net : Numerics.Distributions.Copulas.ArchimedeanCopula
        The .NET copula.

    Notes
    -----
    C(u, v) = φ⁻¹(φ(u) + φ(v)). The generator members accept a scalar or an array, evaluated
    element by element with the bound scalar member.
    """

    def _generator_member(self, name: str, t: Any) -> Any:
        """Evaluate a scalar generator member, bound once per copula object."""
        cache = self.__dict__.setdefault("_generators", {})
        entry = cache.get(name)
        if entry is None or entry[0] is not self.net:
            entry = cache[name] = (self.net, bind(getattr(self.net, name)))
        return loop(t, entry[1])

    def generator(self, t: Any) -> Any:
        """Return the generator φ(t) (``Generator``).

        Parameters
        ----------
        t : float or array_like
            Arguments in (0, 1].

        Returns
        -------
        float or numpy.ndarray
            The value φ(t).
        """
        return self._generator_member("Generator", t)

    def generator_inverse(self, t: Any) -> Any:
        """Return the inverse generator φ⁻¹(t) (``GeneratorInverse``).

        Parameters
        ----------
        t : float or array_like
            Arguments in [0, ∞).

        Returns
        -------
        float or numpy.ndarray
            The value φ⁻¹(t).
        """
        return self._generator_member("GeneratorInverse", t)

    def generator_prime(self, t: Any) -> Any:
        """Return the first derivative φ′(t) (``GeneratorPrime``).

        Parameters
        ----------
        t : float or array_like
            Arguments in (0, 1].

        Returns
        -------
        float or numpy.ndarray
            The value φ′(t).
        """
        return self._generator_member("GeneratorPrime", t)

    def generator_prime2(self, t: Any) -> Any:
        """Return the second derivative φ″(t) (``GeneratorPrime2``).

        Parameters
        ----------
        t : float or array_like
            Arguments in (0, 1].

        Returns
        -------
        float or numpy.ndarray
            The value φ″(t).
        """
        return self._generator_member("GeneratorPrime2", t)

    def generator_prime_inverse(self, t: Any) -> Any:
        """Return the inverse of the first derivative, (φ′)⁻¹(t) (``GeneratorPrimeInverse``).

        Parameters
        ----------
        t : float or array_like
            Values of φ′.

        Returns
        -------
        float or numpy.ndarray
            The argument whose φ′ is ``t``.
        """
        return self._generator_member("GeneratorPrimeInverse", t)

    def set_theta_from_tau(self, x: Any, y: Any) -> None:
        """Set θ by inverting Kendall's tau of a sample (``SetThetaFromTau``).

        Parameters
        ----------
        x : array_like
            Sample of the first variable.
        y : array_like
            Sample of the second variable.
        """
        self.net.SetThetaFromTau(to_net_doubles(x), to_net_doubles(y))


def _build(name: str, parameters: tuple[Any, ...]) -> Any:
    """Construct a .NET copula with its parameters, or with the Numerics defaults when none are given."""
    net_type = getattr(_nc(), name)
    return net_type(*(float(p) for p in parameters)) if parameters else net_type()


class NormalCopula(BivariateCopula):
    """Gaussian copula.

    Parameters
    ----------
    rho : float, optional
        Correlation ρ; default the Numerics default.
    marginal_x : UnivariateDistribution, optional
        Marginal of the first variable.
    marginal_y : UnivariateDistribution, optional
        Marginal of the second variable.
    """

    def __init__(self, rho: float | None = None, marginal_x: UnivariateDistribution | None = None,
                 marginal_y: UnivariateDistribution | None = None) -> None:
        super().__init__(_build("NormalCopula", () if rho is None else (rho,)))
        self._attach(marginal_x, marginal_y)


class StudentTCopula(BivariateCopula):
    """Student-t copula.

    Parameters
    ----------
    rho : float, optional
        Correlation ρ. Give both ``rho`` and ``degrees_of_freedom``, or neither for the Numerics defaults.
    degrees_of_freedom : float, optional
        Degrees of freedom ν.
    marginal_x : UnivariateDistribution, optional
        Marginal of the first variable.
    marginal_y : UnivariateDistribution, optional
        Marginal of the second variable.
    """

    def __init__(self, rho: float | None = None, degrees_of_freedom: float | None = None,
                 marginal_x: UnivariateDistribution | None = None,
                 marginal_y: UnivariateDistribution | None = None) -> None:
        if (rho is None) != (degrees_of_freedom is None):
            raise TypeError("Give both rho and degrees_of_freedom, or neither.")
        parameters = () if rho is None else (rho, degrees_of_freedom)
        super().__init__(_build("StudentTCopula", parameters))
        self._attach(marginal_x, marginal_y)

    @property
    def degrees_of_freedom(self) -> float:
        """The degrees of freedom ν; settable."""
        return float(self.net.DegreesOfFreedom)

    @degrees_of_freedom.setter
    def degrees_of_freedom(self, value: float) -> None:
        self.net.DegreesOfFreedom = float(value)

    def validate_parameters(self, rho: float, degrees_of_freedom: float) -> str | None:
        """Return why a correlation and degrees of freedom are invalid, or None when they are valid.

        Parameters
        ----------
        rho : float
            Candidate correlation.
        degrees_of_freedom : float
            Candidate degrees of freedom.

        Returns
        -------
        str or None
            The Numerics message of ``ValidateParameters(rho, nu, false)``, or None.
        """
        return _message(self.net.ValidateParameters(float(rho), float(degrees_of_freedom), False))


class ClaytonCopula(ArchimedeanCopula):
    """Clayton copula.

    Parameters
    ----------
    theta : float, optional
        Dependence parameter θ; default the Numerics default.
    marginal_x : UnivariateDistribution, optional
        Marginal of the first variable.
    marginal_y : UnivariateDistribution, optional
        Marginal of the second variable.
    """

    def __init__(self, theta: float | None = None, marginal_x: UnivariateDistribution | None = None,
                 marginal_y: UnivariateDistribution | None = None) -> None:
        super().__init__(_build("ClaytonCopula", () if theta is None else (theta,)))
        self._attach(marginal_x, marginal_y)


class GumbelCopula(ArchimedeanCopula):
    """Gumbel copula.

    Parameters
    ----------
    theta : float, optional
        Dependence parameter θ ≥ 1; default the Numerics default.
    marginal_x : UnivariateDistribution, optional
        Marginal of the first variable.
    marginal_y : UnivariateDistribution, optional
        Marginal of the second variable.
    """

    def __init__(self, theta: float | None = None, marginal_x: UnivariateDistribution | None = None,
                 marginal_y: UnivariateDistribution | None = None) -> None:
        super().__init__(_build("GumbelCopula", () if theta is None else (theta,)))
        self._attach(marginal_x, marginal_y)


class FrankCopula(ArchimedeanCopula):
    """Frank copula.

    Parameters
    ----------
    theta : float, optional
        Dependence parameter θ ≠ 0; default the Numerics default.
    marginal_x : UnivariateDistribution, optional
        Marginal of the first variable.
    marginal_y : UnivariateDistribution, optional
        Marginal of the second variable.
    """

    def __init__(self, theta: float | None = None, marginal_x: UnivariateDistribution | None = None,
                 marginal_y: UnivariateDistribution | None = None) -> None:
        super().__init__(_build("FrankCopula", () if theta is None else (theta,)))
        self._attach(marginal_x, marginal_y)


class JoeCopula(ArchimedeanCopula):
    """Joe copula.

    Parameters
    ----------
    theta : float, optional
        Dependence parameter θ ≥ 1; default the Numerics default.
    marginal_x : UnivariateDistribution, optional
        Marginal of the first variable.
    marginal_y : UnivariateDistribution, optional
        Marginal of the second variable.
    """

    def __init__(self, theta: float | None = None, marginal_x: UnivariateDistribution | None = None,
                 marginal_y: UnivariateDistribution | None = None) -> None:
        super().__init__(_build("JoeCopula", () if theta is None else (theta,)))
        self._attach(marginal_x, marginal_y)


class AMHCopula(ArchimedeanCopula):
    """Ali-Mikhail-Haq copula.

    Parameters
    ----------
    theta : float, optional
        Dependence parameter θ in [−1, 1); default the Numerics default.
    marginal_x : UnivariateDistribution, optional
        Marginal of the first variable.
    marginal_y : UnivariateDistribution, optional
        Marginal of the second variable.
    """

    def __init__(self, theta: float | None = None, marginal_x: UnivariateDistribution | None = None,
                 marginal_y: UnivariateDistribution | None = None) -> None:
        super().__init__(_build("AMHCopula", () if theta is None else (theta,)))
        self._attach(marginal_x, marginal_y)


class IndependenceCopula(BivariateCopula):
    """Independence copula C(u, v) = u·v.

    Parameters
    ----------
    marginal_x : UnivariateDistribution, optional
        Marginal of the first variable.
    marginal_y : UnivariateDistribution, optional
        Marginal of the second variable.
    """

    def __init__(self, marginal_x: UnivariateDistribution | None = None,
                 marginal_y: UnivariateDistribution | None = None) -> None:
        super().__init__(_build("IndependenceCopula", ()))
        self._attach(marginal_x, marginal_y)


_CLASSES: dict[str, type[BivariateCopula]] = {
    "NormalCopula": NormalCopula, "StudentTCopula": StudentTCopula, "ClaytonCopula": ClaytonCopula,
    "GumbelCopula": GumbelCopula, "FrankCopula": FrankCopula, "JoeCopula": JoeCopula, "AMHCopula": AMHCopula,
    "IndependenceCopula": IndependenceCopula,
}


def wrap(net: Any) -> BivariateCopula:
    """Wrap a .NET bivariate copula in the matching Python class.

    Parameters
    ----------
    net : Numerics.Distributions.Copulas.BivariateCopula or IBivariateCopula
        The .NET object.

    Returns
    -------
    BivariateCopula
        The specific wrapper class when one exists, otherwise the generic wrapper; the .NET object
        is shared, not copied.
    """
    concrete = unwrap(net)
    cls = _CLASSES.get(str(concrete.GetType().Name), BivariateCopula)
    instance = cls.__new__(cls)
    BivariateCopula.__init__(instance, concrete)
    return instance


def create_copula(copula_type: CopulaType | str) -> BivariateCopula:
    """Create a copula of a family with its default parameters (``CopulaFactory.CreateCopula``).

    Parameters
    ----------
    copula_type : CopulaType or str
        The family.

    Returns
    -------
    BivariateCopula
        The typed copula.
    """
    nc = _nc()
    return wrap(nc.CopulaFactory.CreateCopula(net_enum(nc.CopulaType, copula_type)))


def copula_from_xml(xml: str, include_marginals: bool = False) -> BivariateCopula:
    """Rebuild a copula from :meth:`BivariateCopula.to_xml` text (``CopulaFactory.CreateCopula(XElement)``).

    Parameters
    ----------
    xml : str
        The serialized ``Copula`` element.
    include_marginals : bool, optional
        Also restore the serialized marginals (``CreateCopula(XElement, bool)``).

    Returns
    -------
    BivariateCopula
        The validated typed copula.

    Raises
    ------
    ValueError
        If the element is malformed or its parameters do not define a valid copula.
    """
    nc = _nc()
    element = namespace("System.Xml.Linq").XElement.Parse(str(xml))
    if not include_marginals:
        return wrap(nc.CopulaFactory.CreateCopula(element))
    import System  # noqa: PLC0415 - the runtime is running

    return wrap(nc.CopulaFactory.CreateCopula.Overloads[type(element), System.Boolean](element, True))


__all__ = [
    "AMHCopula", "ArchimedeanCopula", "BivariateCopula", "ClaytonCopula", "CopulaEstimationMethod", "CopulaType",
    "FrankCopula", "GumbelCopula", "IndependenceCopula", "JoeCopula", "NormalCopula", "StudentTCopula",
    "copula_from_xml", "create_copula", "wrap",
]

translate_public(globals())
