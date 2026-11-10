"""Univariate functions and link functions from ``Numerics.Functions``.

:class:`LinearFunction`, :class:`PowerFunction`, :class:`SegmentedPowerFunction`,
:class:`TabularFunction` and :class:`CompositeFunction` wrap the ``IUnivariateFunction``
implementations (rating curves and other deterministic or stochastic transfer functions);
:class:`EnsembleFunction` holds an ensemble of parameter sets for one template.
:class:`IdentityLink`, :class:`LogLink`, :class:`LogitLink`, :class:`ProbitLink`,
:class:`ComplementaryLogLogLink`, :class:`FisherZLink` and :class:`YeoJohnsonLink` wrap the
``ILinkFunction`` implementations, and :class:`LinkController` holds one link per parameter.

An array argument crosses once, through the list overloads ``UnivariateFunctionExtensions`` and
``LinkFunctionExtensions``, which give the scalar members' values bit for bit.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np

from pynumerics._util import evaluate, namespace, net_enum, unwrap
from pynumerics.data import UncertainOrderedPairedData, parse_xml
from pynumerics.enums import Extrapolation, NetEnum, Transform
from pynumerics.errors import translate_public
from pynumerics.marshal import to_net_doubles, to_numpy_doubles, to_numpy_matrix


def _functions() -> Any:
    """Return the ``Numerics.Functions`` namespace after starting the runtime."""
    return namespace("Numerics.Functions")


_LIST_OVERLOADS: dict[tuple[str, str], Any] = {}


def _apply(net: Any, member: str, extensions: str, values: Any) -> Any:
    """Evaluate a scalar member of ``net`` on a scalar or an array, the array through its list overload.

    The list overload is the static ``extensions.member(net, IList<double>)``, which calls the scalar
    member once per element, so both give the scalar member's values bit for bit.
    """
    key = (extensions, member)
    vector = _LIST_OVERLOADS.get(key)
    if vector is None:
        vector = _LIST_OVERLOADS[key] = getattr(getattr(_functions(), extensions), member)
    return evaluate(values, getattr(net, member), lambda array: vector(net, array))


class FunctionType(NetEnum):
    """Univariate function families (``Numerics.Functions.UnivariateFunctionType``)."""

    LINEAR = "Linear"
    POWER = "Power"
    TABULAR = "Tabular"
    SEGMENTED_POWER = "SegmentedPower"
    COMPOSITE = "Composite"


class CompositeMode(NetEnum):
    """How a :class:`CompositeFunction` combines its children (``Numerics.Functions.CompositeFunctionMode``)."""

    WEIGHTED_AVERAGE = "WeightedAverage"
    MIXTURE = "Mixture"


class LinkFunction:
    """A link function ``eta = g(x)`` (``Numerics.Functions.ILinkFunction``).

    Parameters
    ----------
    net : Numerics.Functions.ILinkFunction
        The .NET link to wrap; the subclasses construct their own.

    Attributes
    ----------
    net : Numerics.Functions.ILinkFunction
        The wrapped .NET object.
    """

    def __init__(self, net: Any) -> None:
        self.net = unwrap(net)

    def link(self, x: Any) -> Any:
        """Map values to the link space (``Link``).

        Parameters
        ----------
        x : float or array_like
            Values in the natural space.

        Returns
        -------
        float or numpy.ndarray
            ``g(x)``, shaped like ``x``.
        """
        return _apply(self.net, "Link", "LinkFunctionExtensions", x)

    def inverse_link(self, eta: Any) -> Any:
        """Map link-space values back to the natural space (``InverseLink``).

        Parameters
        ----------
        eta : float or array_like
            Values in the link space.

        Returns
        -------
        float or numpy.ndarray
            ``g^-1(eta)``, shaped like ``eta``.
        """
        return _apply(self.net, "InverseLink", "LinkFunctionExtensions", eta)

    def d_link(self, x: Any) -> Any:
        """Return the derivative ``dg/dx`` (``DLink``).

        Parameters
        ----------
        x : float or array_like
            Values in the natural space.

        Returns
        -------
        float or numpy.ndarray
            The derivative, shaped like ``x``.
        """
        return _apply(self.net, "DLink", "LinkFunctionExtensions", x)

    def to_xml(self) -> str:
        """Return the link's ``ToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.ToXElement().ToString())


class IdentityLink(LinkFunction):
    """The identity link ``g(x) = x`` (``Numerics.Functions.IdentityLink``)."""

    def __init__(self) -> None:
        super().__init__(_functions().IdentityLink())


class LogLink(LinkFunction):
    """The log link ``g(x) = ln x`` for positive parameters (``Numerics.Functions.LogLink``)."""

    def __init__(self) -> None:
        super().__init__(_functions().LogLink())


class LogitLink(LinkFunction):
    """The logit link ``g(x) = ln(x / (1 - x))`` for parameters in (0, 1) (``Numerics.Functions.LogitLink``)."""

    def __init__(self) -> None:
        super().__init__(_functions().LogitLink())


class ProbitLink(LinkFunction):
    """The probit link ``g(x) = Phi^-1(x)`` for parameters in (0, 1) (``Numerics.Functions.ProbitLink``)."""

    def __init__(self) -> None:
        super().__init__(_functions().ProbitLink())


class ComplementaryLogLogLink(LinkFunction):
    """The complementary log-log link ``g(x) = ln(-ln(1 - x))`` (``ComplementaryLogLogLink``)."""

    def __init__(self) -> None:
        super().__init__(_functions().ComplementaryLogLogLink())


class FisherZLink(LinkFunction):
    """The Fisher z link ``g(x) = atanh(x)`` for correlations in (-1, 1) (``Numerics.Functions.FisherZLink``)."""

    def __init__(self) -> None:
        super().__init__(_functions().FisherZLink())


class YeoJohnsonLink(LinkFunction):
    """The Yeo-Johnson power link (``Numerics.Functions.YeoJohnsonLink``).

    Parameters
    ----------
    lambda_ : float, optional
        The transformation parameter; default 1 (the identity).
    values : array_like, optional
        Representative values (at least two, finite) to fit lambda to by maximum likelihood
        (``YeoJohnson.FitLambda``) instead of giving ``lambda_``.
    """

    def __init__(self, lambda_: float = 1.0, values: Any = None) -> None:
        functions = _functions()
        if values is not None:
            super().__init__(functions.YeoJohnsonLink(to_net_doubles(values)))
        else:
            super().__init__(functions.YeoJohnsonLink(float(lambda_)))

    @property
    def lambda_(self) -> float:
        """The transformation parameter."""
        return float(self.net.Lambda)


_LINK_CLASSES: dict[str, type[LinkFunction]] = {
    cls.__name__: cls for cls in (IdentityLink, LogLink, LogitLink, ProbitLink, ComplementaryLogLogLink,
                                  FisherZLink, YeoJohnsonLink)}


def wrap_link(net: Any) -> LinkFunction | None:
    """Wrap a .NET link function in the matching Python class.

    Parameters
    ----------
    net : Numerics.Functions.ILinkFunction or None
        The .NET link.

    Returns
    -------
    LinkFunction or None
        The specific wrapper class (the generic :class:`LinkFunction` for an unknown type), or
        None for a null link.
    """
    if net is None:
        return None
    concrete = unwrap(net)
    cls = _LINK_CLASSES.get(str(concrete.GetType().Name), LinkFunction)
    instance = cls.__new__(cls)
    LinkFunction.__init__(instance, concrete)
    return instance


def create_link(link_type: str) -> LinkFunction:
    """Create a link function by type name (``LinkFunctionFactory.Create``).

    Parameters
    ----------
    link_type : str
        A ``LinkFunctionType`` member name: ``"identity"``, ``"log"``, ``"logit"``, ``"probit"``,
        ``"complementary_log_log"``, ``"yeo_johnson"`` or ``"fisher_z"``.

    Returns
    -------
    LinkFunction
        The wrapped link (a Yeo-Johnson link has lambda 1).
    """
    functions = _functions()
    return wrap_link(functions.LinkFunctionFactory.Create(net_enum(functions.LinkFunctionType, link_type)))


def link_from_xml(xml: str) -> LinkFunction:
    """Restore a link function from its XML (``LinkFunctionFactory.CreateFromXElement``).

    Parameters
    ----------
    xml : str
        The text :meth:`LinkFunction.to_xml` returns.

    Returns
    -------
    LinkFunction
        The wrapped link.
    """
    return wrap_link(_functions().LinkFunctionFactory.CreateFromXElement(parse_xml(xml)))


def _net_link(link: Any) -> Any:
    """Return the .NET link for a wrapper, a ``LinkFunctionType`` name, a .NET link or None."""
    if link is None:
        return None
    if isinstance(link, LinkFunction):
        return link.net
    if isinstance(link, str):
        return create_link(link).net
    return link


class LinkController:
    """Per-parameter link functions (``Numerics.Functions.LinkController``).

    Maps parameters to and from the link space, with the Jacobian, the log-determinant of the
    inverse-link Jacobian, indexing and the XML round trip.

    Parameters
    ----------
    links : sequence of LinkFunction, str or None
        One link per parameter: a :class:`LinkFunction`, a ``LinkFunctionType`` member name, or
        None for a parameter without a link.

    Attributes
    ----------
    net : Numerics.Functions.LinkController
        The wrapped .NET object.
    """

    def __init__(self, links: Sequence[Any]) -> None:
        functions = _functions()
        import System  # noqa: PLC0415 - the runtime is running

        self.net = functions.LinkController(System.Array[functions.ILinkFunction]([_net_link(v) for v in links]))

    @classmethod
    def _from_net(cls, net: Any) -> LinkController:
        """Wrap an existing .NET controller."""
        instance = cls.__new__(cls)
        instance.net = net
        return instance

    @classmethod
    def for_location_scale_shape(cls, location: Any = None, scale: Any = None, shape: Any = None) -> LinkController:
        """Return a three-parameter controller (``ForLocationScaleShape``).

        Parameters
        ----------
        location, scale, shape : LinkFunction, str or None, optional
            The link of each parameter; None leaves the parameter unlinked.

        Returns
        -------
        LinkController
            The controller.
        """
        return cls._from_net(_functions().LinkController.ForLocationScaleShape(
            _net_link(location), _net_link(scale), _net_link(shape)))

    @classmethod
    def from_xml(cls, xml: str) -> LinkController:
        """Restore a controller from its XML (the ``XElement`` constructor).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        LinkController
            The controller.
        """
        return cls._from_net(_functions().LinkController(parse_xml(xml)))

    def __len__(self) -> int:
        return int(self.net.Count)

    @property
    def count(self) -> int:
        """The number of parameter slots."""
        return int(self.net.Count)

    def __getitem__(self, index: int) -> LinkFunction | None:
        return wrap_link(self.net[int(index)])

    @property
    def links(self) -> list[LinkFunction | None]:
        """The link of each parameter (None where unlinked)."""
        return [wrap_link(link) for link in self.net.Links]

    def link(self, x: Any) -> np.ndarray:
        """Map parameters to the link space (``Link``).

        Parameters
        ----------
        x : array_like
            Parameters, one per link.

        Returns
        -------
        numpy.ndarray
            The linked values.
        """
        return to_numpy_doubles(self.net.Link(to_net_doubles(x)))

    def inverse_link(self, eta: Any) -> np.ndarray:
        """Map link-space values back to parameters (``InverseLink``).

        Parameters
        ----------
        eta : array_like
            Link-space values, one per link.

        Returns
        -------
        numpy.ndarray
            The parameters.
        """
        return to_numpy_doubles(self.net.InverseLink(to_net_doubles(eta)))

    def link_jacobian(self, x: Any) -> np.ndarray:
        """Return the diagonal Jacobian d(eta)/d(x) at the parameters (``LinkJacobian``).

        Parameters
        ----------
        x : array_like
            Parameters, one per link.

        Returns
        -------
        numpy.ndarray
            The square Jacobian matrix.
        """
        return to_numpy_matrix(self.net.LinkJacobian(to_net_doubles(x)).ToArray())

    def log_det_jacobian(self, phi: Any) -> float:
        """Return ``log |d theta / d phi|`` of the inverse link at link-space values (``LogDetJacobian``).

        Parameters
        ----------
        phi : array_like
            Link-space values, one per parameter.

        Returns
        -------
        float
            The log-determinant, the change-of-variables term of a prior in the link space.
        """
        return float(self.net.LogDetJacobian(to_net_doubles(phi)))

    def to_xml(self) -> str:
        """Return the controller's ``ToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.ToXElement().ToString())


class UnivariateFunction:
    """A deterministic or stochastic function ``y = f(x)`` (``Numerics.Functions.IUnivariateFunction``).

    A stochastic function evaluates at :attr:`confidence_level` of its residual; with the level
    outside [0, 1] or :attr:`is_deterministic` set, it evaluates the mean curve.

    Parameters
    ----------
    net : Numerics.Functions.IUnivariateFunction
        The .NET function to wrap; the subclasses construct their own.

    Attributes
    ----------
    net : Numerics.Functions.IUnivariateFunction
        The wrapped .NET object.
    """

    def __init__(self, net: Any) -> None:
        self.net = unwrap(net)

    @property
    def number_of_parameters(self) -> int:
        """The length of the parameter vector."""
        return int(self.net.NumberOfParameters)

    @property
    def parameters_valid(self) -> bool:
        """Whether the current parameters are valid."""
        return bool(self.net.ParametersValid)

    @property
    def minimum(self) -> float:
        """The lower bound of the support (settable where the function allows it)."""
        return float(self.net.Minimum)

    @minimum.setter
    def minimum(self, value: float) -> None:
        self.net.Minimum = float(value)

    @property
    def maximum(self) -> float:
        """The upper bound of the support."""
        return float(self.net.Maximum)

    @maximum.setter
    def maximum(self, value: float) -> None:
        self.net.Maximum = float(value)

    @property
    def minimum_of_parameters(self) -> np.ndarray:
        """The lower bound of each parameter."""
        return to_numpy_doubles(self.net.MinimumOfParameters)

    @property
    def maximum_of_parameters(self) -> np.ndarray:
        """The upper bound of each parameter."""
        return to_numpy_doubles(self.net.MaximumOfParameters)

    @property
    def is_deterministic(self) -> bool:
        """Whether the function ignores its residual."""
        return bool(self.net.IsDeterministic)

    @is_deterministic.setter
    def is_deterministic(self, value: bool) -> None:
        self.net.IsDeterministic = bool(value)

    @property
    def confidence_level(self) -> float:
        """The residual quantile evaluated; outside [0, 1] (default -1) the mean curve."""
        return float(self.net.ConfidenceLevel)

    @confidence_level.setter
    def confidence_level(self, value: float) -> None:
        self.net.ConfidenceLevel = float(value)

    def set_parameters(self, parameters: Any) -> None:
        """Set the parameter vector (``SetParameters``); invalid values are recorded, not refused.

        Parameters
        ----------
        parameters : array_like
            The parameters in Numerics order.
        """
        self.net.SetParameters(to_net_doubles(parameters))

    def validate_parameters(self, parameters: Any, throw: bool = False) -> str | None:
        """Check a parameter vector (``ValidateParameters``).

        Parameters
        ----------
        parameters : array_like
            The parameters in Numerics order.
        throw : bool, optional
            Whether to raise instead of returning the message; default False.

        Returns
        -------
        str or None
            The validation message, or None when the parameters are valid.

        Raises
        ------
        ValueError
            If ``throw`` is True and the parameters are invalid.
        """
        error = self.net.ValidateParameters(to_net_doubles(parameters), bool(throw))
        return None if error is None else str(error.Message)

    def function(self, x: Any) -> Any:
        """Evaluate ``y = f(x)`` (``Function``).

        Parameters
        ----------
        x : float or array_like
            Arguments.

        Returns
        -------
        float or numpy.ndarray
            The function values, shaped like ``x``.
        """
        return _apply(self.net, "Function", "UnivariateFunctionExtensions", x)

    def inverse_function(self, y: Any) -> Any:
        """Evaluate ``x = f^-1(y)`` (``InverseFunction``).

        Parameters
        ----------
        y : float or array_like
            Function values.

        Returns
        -------
        float or numpy.ndarray
            The arguments, shaped like ``y``.
        """
        return _apply(self.net, "InverseFunction", "UnivariateFunctionExtensions", y)

    def __call__(self, x: Any) -> Any:
        return self.function(x)

    def to_xml(self) -> str:
        """Return the function's ``ToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.ToXElement().ToString())

    @property
    def function_type(self) -> FunctionType:
        """The family (``UnivariateFunctionFactory.GetFunctionType``)."""
        return FunctionType(str(_functions().UnivariateFunctionFactory.GetFunctionType(self.net)))


class LinearFunction(UnivariateFunction):
    """The linear function ``y = alpha + beta x`` with a Normal residual (``Numerics.Functions.LinearFunction``).

    Parameters
    ----------
    alpha, beta : float, optional
        Intercept and slope; default 0 and 1.
    sigma : float, optional
        Residual standard error; when given the function is stochastic, otherwise deterministic.
    """

    def __init__(self, alpha: float = 0.0, beta: float = 1.0, sigma: float | None = None) -> None:
        functions = _functions()
        net = (functions.LinearFunction(float(alpha), float(beta)) if sigma is None
               else functions.LinearFunction(float(alpha), float(beta), float(sigma)))
        super().__init__(net)

    @classmethod
    def from_xml(cls, xml: str) -> LinearFunction:
        """Restore the function from its XML (``FromXElement``).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        LinearFunction
            The function.
        """
        return wrap_function(_functions().LinearFunction.FromXElement(parse_xml(xml)))

    @property
    def alpha(self) -> float:
        """The intercept."""
        return float(self.net.Alpha)

    @alpha.setter
    def alpha(self, value: float) -> None:
        self.net.Alpha = float(value)

    @property
    def beta(self) -> float:
        """The slope."""
        return float(self.net.Beta)

    @beta.setter
    def beta(self, value: float) -> None:
        self.net.Beta = float(value)

    @property
    def sigma(self) -> float:
        """The residual standard error."""
        return float(self.net.Sigma)

    @sigma.setter
    def sigma(self, value: float) -> None:
        self.net.Sigma = float(value)


class PowerFunction(UnivariateFunction):
    """The power function ``y = alpha (x - xi)^beta`` with a log-space Normal residual (``PowerFunction``).

    Parameters
    ----------
    alpha, beta : float, optional
        Coefficient and exponent; default 1 and 1.5.
    xi : float, optional
        Location (the lower bound of the support); default 0.
    sigma : float, optional
        Residual standard error in natural-log space; when given the function is stochastic.
    """

    def __init__(self, alpha: float = 1.0, beta: float = 1.5, xi: float = 0.0, sigma: float | None = None) -> None:
        functions = _functions()
        net = (functions.PowerFunction(float(alpha), float(beta), float(xi)) if sigma is None
               else functions.PowerFunction(float(alpha), float(beta), float(xi), float(sigma)))
        super().__init__(net)

    @classmethod
    def from_xml(cls, xml: str) -> PowerFunction:
        """Restore the function from its XML (``FromXElement``).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        PowerFunction
            The function.
        """
        return wrap_function(_functions().PowerFunction.FromXElement(parse_xml(xml)))

    @property
    def alpha(self) -> float:
        """The coefficient."""
        return float(self.net.Alpha)

    @alpha.setter
    def alpha(self, value: float) -> None:
        self.net.Alpha = float(value)

    @property
    def beta(self) -> float:
        """The exponent."""
        return float(self.net.Beta)

    @beta.setter
    def beta(self, value: float) -> None:
        self.net.Beta = float(value)

    @property
    def xi(self) -> float:
        """The location (lower bound of the support)."""
        return float(self.net.Xi)

    @xi.setter
    def xi(self, value: float) -> None:
        self.net.Xi = float(value)

    @property
    def sigma(self) -> float:
        """The residual standard error in natural-log space."""
        return float(self.net.Sigma)

    @sigma.setter
    def sigma(self, value: float) -> None:
        self.net.Sigma = float(value)

    @property
    def is_inverse(self) -> bool:
        """Whether :meth:`function` evaluates the inverse power law (x and y swapped)."""
        return bool(self.net.IsInverse)

    @is_inverse.setter
    def is_inverse(self, value: bool) -> None:
        self.net.IsInverse = bool(value)


class SegmentedPowerFunction(UnivariateFunction):
    """A sum of power-law controls ``sum_k alpha_k (x - h_k)^beta_k [x > h_k]`` (``SegmentedPowerFunction``).

    The parameter layout is ``[h_1, log10 alpha_1, beta_1, ..., h_K, log10 alpha_K, beta_K, sigma]``
    with a Gaussian residual in log10 space (Le Coz et al., 2014, J. Hydrol. 509:573-587).

    Parameters
    ----------
    number_of_segments : int, optional
        Segments of a deterministic function at the default parameters; default 1.
    parameters : array_like, optional
        A full parameter vector (length ``3 K + 1``) instead; the function is then stochastic.
    """

    def __init__(self, number_of_segments: int = 1, parameters: Any = None) -> None:
        functions = _functions()
        net = (functions.SegmentedPowerFunction(int(number_of_segments)) if parameters is None
               else functions.SegmentedPowerFunction(to_net_doubles(parameters)))
        super().__init__(net)

    @classmethod
    def from_xml(cls, xml: str) -> SegmentedPowerFunction:
        """Restore the function from its XML (``FromXElement``).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        SegmentedPowerFunction
            The function.
        """
        return wrap_function(_functions().SegmentedPowerFunction.FromXElement(parse_xml(xml)))

    @property
    def number_of_segments(self) -> int:
        """The number of controls."""
        return int(self.net.NumberOfSegments)

    @property
    def sigma(self) -> float:
        """The residual standard error in log10 space."""
        return float(self.net.Sigma)

    def get_breakpoint(self, segment: int) -> float:
        """Return the activation stage ``h_k`` of a segment (``GetBreakpoint``).

        Parameters
        ----------
        segment : int
            The one-based segment number.

        Returns
        -------
        float
            The breakpoint.
        """
        return float(self.net.GetBreakpoint(int(segment)))

    def get_log10_alpha(self, segment: int) -> float:
        """Return ``log10 alpha_k`` of a segment (``GetLog10Alpha``).

        Parameters
        ----------
        segment : int
            The one-based segment number.

        Returns
        -------
        float
            The log10 coefficient.
        """
        return float(self.net.GetLog10Alpha(int(segment)))

    def get_beta(self, segment: int) -> float:
        """Return the exponent ``beta_k`` of a segment (``GetBeta``).

        Parameters
        ----------
        segment : int
            The one-based segment number.

        Returns
        -------
        float
            The exponent.
        """
        return float(self.net.GetBeta(int(segment)))


class TabularFunction(UnivariateFunction):
    """A function tabulated by uncertain paired data (``Numerics.Functions.TabularFunction``).

    The curve at :attr:`confidence_level` is ``paired_data.curve_sample(level)`` (the means
    when the level is negative), interpolated in the transformed spaces.

    Parameters
    ----------
    paired_data : UncertainOrderedPairedData
        The x values and the distribution of y at each.
    """

    def __init__(self, paired_data: UncertainOrderedPairedData) -> None:
        super().__init__(_functions().TabularFunction(paired_data.net))

    @classmethod
    def from_xml(cls, xml: str) -> TabularFunction:
        """Restore the function from its XML (``FromXElement``).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        TabularFunction
            The function.
        """
        return wrap_function(_functions().TabularFunction.FromXElement(parse_xml(xml)))

    @property
    def paired_data(self) -> UncertainOrderedPairedData:
        """The uncertain paired data."""
        return UncertainOrderedPairedData._from_net(self.net.PairedData)

    @property
    def x_transform(self) -> Transform:
        """The interpolation space of x."""
        return Transform(str(self.net.XTransform))

    @x_transform.setter
    def x_transform(self, value: Transform | str) -> None:
        self.net.XTransform = net_enum(namespace("Numerics.Data").Transform, value)

    @property
    def y_transform(self) -> Transform:
        """The interpolation space of y."""
        return Transform(str(self.net.YTransform))

    @y_transform.setter
    def y_transform(self, value: Transform | str) -> None:
        self.net.YTransform = net_enum(namespace("Numerics.Data").Transform, value)

    @property
    def extrapolation(self) -> Extrapolation:
        """The sides on which the table extends linearly; default none (end values held)."""
        return Extrapolation(str(self.net.Extrapolation))

    @extrapolation.setter
    def extrapolation(self, value: Extrapolation | str) -> None:
        self.net.Extrapolation = net_enum(namespace("Numerics.Data").ExtrapolationSides, value)

    @property
    def allow_negative_y_values(self) -> bool:
        """Whether negative or NaN values are kept; when False they become 0."""
        return bool(self.net.AllowNegativeYValues)

    @allow_negative_y_values.setter
    def allow_negative_y_values(self, value: bool) -> None:
        self.net.AllowNegativeYValues = bool(value)


class CompositeFunction(UnivariateFunction):
    """A weighted combination of univariate functions (``Numerics.Functions.CompositeFunction``).

    With :attr:`mode` weighted average, a stochastic evaluation draws every child at the same
    confidence level; with mixture, the level selects one child by its weight.

    Parameters
    ----------
    functions : sequence of UnivariateFunction
        The children (at least one).
    weights : array_like, optional
        Non-negative weights summing to one, one per child; default equal weights.
    """

    def __init__(self, functions: Sequence[UnivariateFunction], weights: Any = None) -> None:
        net_functions = _functions()
        import System  # noqa: PLC0415 - the runtime is running

        array = System.Array[net_functions.IUnivariateFunction]([f.net for f in functions])
        net = (net_functions.CompositeFunction(array) if weights is None
               else net_functions.CompositeFunction(array, to_net_doubles(weights)))
        super().__init__(net)

    @classmethod
    def from_xml(cls, xml: str) -> CompositeFunction:
        """Restore the function from its XML (``FromXElement``).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        CompositeFunction
            The function.
        """
        return wrap_function(_functions().CompositeFunction.FromXElement(parse_xml(xml)))

    @property
    def mode(self) -> CompositeMode:
        """How the children combine."""
        return CompositeMode(str(self.net.Mode))

    @mode.setter
    def mode(self, value: CompositeMode | str) -> None:
        self.net.Mode = net_enum(_functions().CompositeFunctionMode, value)

    @property
    def functions(self) -> list[UnivariateFunction]:
        """The children."""
        return [wrap_function(f) for f in self.net.Functions]

    @property
    def weights(self) -> np.ndarray:
        """The child weights (the parameter vector)."""
        return to_numpy_doubles(unwrap(self.net.Weights))


_FUNCTION_CLASSES: dict[str, type[UnivariateFunction]] = {
    cls.__name__: cls for cls in (LinearFunction, PowerFunction, SegmentedPowerFunction, TabularFunction,
                                  CompositeFunction)}


def wrap_function(net: Any) -> UnivariateFunction:
    """Wrap a .NET univariate function in the matching Python class.

    Parameters
    ----------
    net : Numerics.Functions.IUnivariateFunction
        The .NET function.

    Returns
    -------
    UnivariateFunction
        The specific wrapper class (the generic :class:`UnivariateFunction` for an unknown type).
    """
    concrete = unwrap(net)
    cls = _FUNCTION_CLASSES.get(str(concrete.GetType().Name), UnivariateFunction)
    instance = cls.__new__(cls)
    UnivariateFunction.__init__(instance, concrete)
    return instance


def create_function(function_type: FunctionType | str) -> UnivariateFunction:
    """Create a function at its default parameters (``UnivariateFunctionFactory.CreateFunction``).

    Parameters
    ----------
    function_type : FunctionType or str
        The family: ``"linear"``, ``"power"`` or ``"segmented_power"`` (a tabular or composite
        function needs its data, so the factory refuses them).

    Returns
    -------
    UnivariateFunction
        The wrapped function.
    """
    functions = _functions()
    return wrap_function(functions.UnivariateFunctionFactory.CreateFunction(
        net_enum(functions.UnivariateFunctionType, function_type)))


def function_from_xml(xml: str) -> UnivariateFunction:
    """Restore any univariate function from its XML (``UnivariateFunctionFactory.CreateFromXElement``).

    Parameters
    ----------
    xml : str
        The text :meth:`UnivariateFunction.to_xml` returns.

    Returns
    -------
    UnivariateFunction
        The wrapped function.
    """
    return wrap_function(_functions().UnivariateFunctionFactory.CreateFromXElement(parse_xml(xml)))


class EnsembleFunction:
    """An ensemble of parameter sets for one function template (``Numerics.Functions.EnsembleFunction``).

    Each member is the template with one parameter set (for example a posterior sample of a
    rating curve); :meth:`sample` picks the member at a percentile of the ensemble order.

    Parameters
    ----------
    template : UnivariateFunction
        The function whose configuration every member copies (not a :class:`TabularFunction`).
    parameter_sets : array_like
        Shape ``(members, template.number_of_parameters)``: one finite, valid parameter vector per row.
    fitness : array_like, optional
        One fitness value per member (for example the log-posterior); default 0.
    weights : array_like, optional
        One weight per member; default 0.

    Attributes
    ----------
    net : Numerics.Functions.EnsembleFunction
        The wrapped .NET object.
    """

    def __init__(self, template: UnivariateFunction, parameter_sets: Any, fitness: Any = None,
                 weights: Any = None) -> None:
        rows = np.atleast_2d(np.asarray(parameter_sets, dtype=np.float64))
        count = rows.shape[0]
        fit = np.zeros(count) if fitness is None else np.asarray(fitness, dtype=np.float64).reshape(count)
        weight = np.zeros(count) if weights is None else np.asarray(weights, dtype=np.float64).reshape(count)
        optimization = namespace("Numerics.Mathematics.Optimization")
        import System  # noqa: PLC0415 - the runtime is running

        sets = System.Array[optimization.ParameterSet]([
            optimization.ParameterSet(to_net_doubles(row), float(f), float(w))
            for row, f, w in zip(rows, fit, weight, strict=True)])
        self.net = _functions().EnsembleFunction(template.net, sets)

    @classmethod
    def from_xml(cls, xml: str) -> EnsembleFunction:
        """Restore the ensemble from its XML (``FromXElement``).

        Parameters
        ----------
        xml : str
            The text :meth:`to_xml` returns.

        Returns
        -------
        EnsembleFunction
            The ensemble.
        """
        instance = cls.__new__(cls)
        instance.net = _functions().EnsembleFunction.FromXElement(parse_xml(xml))
        return instance

    def __len__(self) -> int:
        return int(self.net.Count)

    @property
    def count(self) -> int:
        """The number of members."""
        return int(self.net.Count)

    @property
    def parameter_sets(self) -> np.ndarray:
        """The parameter vectors, shape ``(members, parameters)`` (copies)."""
        return np.array([to_numpy_doubles(s.Values) for s in self.net.ParameterSets], dtype=np.float64)

    @property
    def fitness(self) -> np.ndarray:
        """The fitness of each member."""
        return np.array([float(s.Fitness) for s in self.net.ParameterSets], dtype=np.float64)

    @property
    def weights(self) -> np.ndarray:
        """The weight of each member."""
        return np.array([float(s.Weight) for s in self.net.ParameterSets], dtype=np.float64)

    def sample_at(self, index: int) -> UnivariateFunction:
        """Return the member at an index as an independent function (``SampleAt``).

        Parameters
        ----------
        index : int
            The member, in [0, count).

        Returns
        -------
        UnivariateFunction
            A copy of the template with that member's parameters.
        """
        return wrap_function(self.net.SampleAt(int(index)))

    def sample(self, percentile: float) -> UnivariateFunction:
        """Return the member at a percentile of the ensemble order (``Sample``).

        Parameters
        ----------
        percentile : float
            In [0, 1]; the member index is ``floor(percentile * count)``, capped at the last.

        Returns
        -------
        UnivariateFunction
            A copy of the template with that member's parameters.
        """
        return wrap_function(self.net.Sample(float(percentile)))

    def to_xml(self) -> str:
        """Return the ensemble's ``ToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it.
        """
        return str(self.net.ToXElement().ToString())


__all__ = [
    "ComplementaryLogLogLink", "CompositeFunction", "CompositeMode", "EnsembleFunction", "FisherZLink",
    "FunctionType", "IdentityLink", "LinearFunction", "LinkController", "LinkFunction", "LogLink", "LogitLink",
    "PowerFunction", "ProbitLink", "SegmentedPowerFunction", "TabularFunction", "UnivariateFunction",
    "YeoJohnsonLink", "create_function", "create_link", "function_from_xml", "link_from_xml", "wrap_function",
    "wrap_link",
]

translate_public(globals())
