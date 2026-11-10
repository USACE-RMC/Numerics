"""Univariate distributions from ``Numerics.Distributions``.

Every class wraps one .NET distribution object (``.net``). Evaluation members accept a
scalar or any array: a scalar calls the scalar member, an array crosses once as a
``System.Double[]`` into the ``UnivariateDistributionBase`` array overloads
(``PDF(IList<double>)``, ``CDF(IList<double>)``, ``InverseCDF(IList<double>)``), and
random generation uses ``GenerateRandomValues(n, seed)``, so no numerical work happens
in Python. Parameters follow the Numerics conventions (for example the generalized
extreme value shape is Hosking's kappa, and ``LogNormal`` is base 10 while ``LnNormal``
is natural-log with real-space mean and standard deviation).

Distributions compare by value: two wrappers are equal when they wrap the same .NET family with
bitwise-equal parameters (and the same ``Base`` where the family has one). A mutable distribution
is not hashable, because its value can change; :meth:`UnivariateDistribution.frozen` returns a
read-only :class:`FrozenDistribution` that is hashable by the same value. :data:`FAMILIES` maps each
Numerics type name to its typed class.
"""

from __future__ import annotations

import types
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np

from pynumerics import enums
from pynumerics._util import (
    _data,
    _nd,
    bind,
    bind_list,
    cached_overload,
    evaluate,
    list_overload,
    loop,
    namespace,
    net_enum,
    unwrap,
)
from pynumerics.errors import net_errors, translate, translate_class, translate_public, translated
from pynumerics.marshal import (
    to_net_doubles,
    to_net_ints,
    to_net_matrix,
    to_numpy_bools,
    to_numpy_doubles,
    to_numpy_matrix,
)

_PROTOTYPES: dict[str, Any] = {}


def _prototype(type_name: str) -> Any:
    """Return a default-parameter .NET instance of a family, created once per process.

    Used read-only, for the class-level ``parameter_names`` and ``validate_parameters``. A family
    the factory cannot build without parameters (``Mixture`` needs its components) has none: the
    factory's ``NotSupportedException`` is raised as ``NotImplementedError``.
    """
    net = _PROTOTYPES.get(type_name)
    if net is None:
        nd = _nd()
        with net_errors():
            net = _PROTOTYPES[type_name] = unwrap(nd.UnivariateDistributionFactory.CreateDistribution(
                net_enum(nd.UnivariateDistributionType, type_name)))
    return net


def _fixed_arity(net: Any) -> bool:
    """Return whether the family's parameter vector describes it (``GetParameters`` has ``NumberOfParameters`` values).

    False for the table families (``EmpiricalDistribution``, ``KernelDensity``), whose state is not
    a parameter vector; their arity is not checked.
    """
    return len(net.GetParameters) == int(net.NumberOfParameters)


def _check_arity(net: Any, parameters: Any) -> np.ndarray:
    """Return ``parameters`` as float64, or raise ValueError when their count is not the family's."""
    values = np.atleast_1d(np.asarray(parameters, dtype=np.float64))
    if values.ndim != 1:
        raise ValueError(f"Expected a one-dimensional parameter vector; got shape {values.shape}.")
    expected = int(net.NumberOfParameters)
    if _fixed_arity(net) and values.size != expected:
        names = ", ".join(str(name) for name in net.GetParameterPropertyNames)
        raise ValueError(f"{net.Type} takes {expected} parameters ({names}); got {values.size}.")
    return values


def _net_of(owner: Any) -> Any:
    """Return the .NET object a family method runs on: the instance's, or the typed class's default instance.

    Raises
    ------
    TypeError
        If ``owner`` is the generic :class:`UnivariateDistribution` class.
    """
    if isinstance(owner, type):
        type_name = _class_type_name(owner)
        if type_name is None:
            raise TypeError("This method needs a typed family class or an instance.")
        return _prototype(type_name)
    return owner.net


def _interface(net: Any, name: str) -> Any:
    """Return ``net`` viewed through the ``Numerics.Distributions`` interface ``name``.

    Raises
    ------
    NotImplementedError
        If the family does not implement the interface.
    """
    import clr  # noqa: PLC0415 - the runtime is running (net exists)

    interface = getattr(_nd(), name)
    if not clr.GetClrType(interface).IsInstanceOfType(net):
        raise NotImplementedError(f"{net.DisplayName} does not implement {name}.")
    return interface(net)


def _net_property(name: str, doc: str, kind: Callable[[Any], Any] = float, settable: bool = True) -> property:
    """Return a property that reads (and, when ``settable``, writes) the .NET property ``name`` of ``self.net``.

    The setter calls the owner's ``_ensure_mutable`` first when it has one (distributions do) and
    passes the value through ``kind``; Numerics validates it as its own setter does.
    """

    def fget(self: Any) -> Any:
        return kind(getattr(self.net, name))

    if not settable:
        return property(fget, doc=doc)

    def fset(self: Any, value: Any) -> None:
        ensure_mutable = getattr(self, "_ensure_mutable", None)
        if ensure_mutable is not None:
            ensure_mutable()
        setattr(self.net, name, kind(value))

    return property(fget, fset, doc=doc)


def _enum_property(name: str, net_type: Callable[[], Any], doc: str) -> property:
    """Return a property that reads the .NET enum property ``name`` as its member name and sets it by name."""

    def fget(self: Any) -> str:
        return str(getattr(self.net, name))

    def fset(self: Any, value: Any) -> None:
        self._ensure_mutable()
        setattr(self.net, name, net_enum(net_type(), value))

    return property(fget, fset, doc=doc)


def _bins_net(bins: Any) -> Any:
    """Return the .NET ``List<StratificationBin>`` of a :class:`pynumerics.sampling.StratificationBins`.

    A sequence of them becomes a ``List<StratificationBin>[]`` (the many-range overloads).
    """
    if hasattr(bins, "net"):
        return bins.net
    import System  # noqa: PLC0415 - the runtime is running
    from System.Collections.Generic import List  # noqa: PLC0415

    sampling = namespace("Numerics.Sampling")
    return System.Array[List[sampling.StratificationBin]]([b.net for b in bins])


def _members_array(distributions: Sequence[UnivariateDistribution]) -> Any:
    """Return a ``UnivariateDistributionBase[]`` of the wrapped .NET objects."""
    import System  # noqa: PLC0415 - the runtime is running

    return System.Array[_nd().UnivariateDistributionBase]([d.net for d in distributions])


def _parse_xml(xml: str) -> Any:
    """Parse XML text into a ``System.Xml.Linq.XElement``."""
    namespace("Numerics")
    import clr  # noqa: PLC0415

    clr.AddReference("System.Xml.Linq")
    from System.Xml.Linq import XElement  # noqa: PLC0415

    return XElement.Parse(str(xml))


_STATE: dict[str, tuple[bool, bool]] = {}


def _state_info(net: Any) -> tuple[bool, bool]:
    """Return (has ``Base``, needs the XML) for the .NET type of ``net``, cached per type.

    The XML joins the value key where the parameter vector is not the whole state: families that
    write their own ``ToXElement`` (mixtures, tables, competing risks) or whose ``GetParameters`` is
    not a full parameter vector.
    """
    net_type = net.GetType()
    key = str(net_type.FullName)
    info = _STATE.get(key)
    if info is None:
        declaring = net_type.GetMethod("ToXElement", _empty_types()).DeclaringType
        own_xml = str(declaring.FullName) != "Numerics.Distributions.UnivariateDistributionBase"
        info = _STATE[key] = (net_type.GetProperty("Base") is not None, own_xml or not _fixed_arity(net))
    return info


def _empty_types() -> Any:
    """Return an empty ``System.Type[]`` (the signature of a parameterless method lookup)."""
    import System  # noqa: PLC0415 - the runtime is running

    return System.Type.EmptyTypes


class _FamilyAttribute:
    """A read-only attribute with a class-level and an instance-level value.

    On a typed class (``GeneralizedExtremeValue.type_name``) it gives the family's value without an
    instance; on an instance it reads the wrapped .NET object, so the generic wrapper of any
    family answers too.
    """

    def __init__(self, on_class: Callable[[type], Any], on_instance: Callable[[Any], Any], doc: str) -> None:
        self._on_class = on_class
        self._on_instance = on_instance
        self.__doc__ = doc

    def __get__(self, obj: Any, owner: type) -> Any:
        if obj is None:
            return self._on_class(owner)
        return self._on_instance(obj)


class _FamilyMethod:
    """A method callable on a typed class (on the family's default instance) or on an instance."""

    def __init__(self, func: Callable[..., Any]) -> None:
        self._func = translated(func)
        self.__doc__ = func.__doc__
        self.__wrapped__ = func

    def __get__(self, obj: Any, owner: type) -> Any:
        return types.MethodType(self._func, owner if obj is None else obj)


def _class_type_name(cls: type) -> str | None:
    return getattr(cls, "_type_name", None)


def _class_parameter_names(cls: type) -> tuple[str, ...] | None:
    type_name = _class_type_name(cls)
    if type_name is None:
        return None
    try:
        prototype = _prototype(type_name)
    except NotImplementedError:
        return None  # the family has no default instance, so its names come from an instance only
    return tuple(str(name) for name in prototype.GetParameterPropertyNames)


class UnivariateDistribution:
    """A wrapped ``Numerics.Distributions.UnivariateDistributionBase``.

    Parameters
    ----------
    net : Numerics.Distributions.UnivariateDistributionBase
        The .NET distribution object; interface-typed values are unwrapped to their
        concrete type.

    Attributes
    ----------
    net : Numerics.Distributions.UnivariateDistributionBase
        The wrapped .NET object, for passing back into Numerics members.
    """

    _type_name: str | None = None

    def __init__(self, net: Any) -> None:
        self.net = unwrap(net)

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        translate_class(cls)

    def __deepcopy__(self, memo: dict) -> UnivariateDistribution:
        """Return a wrapper of the distribution's ``Clone``, so a deep copy owns its own .NET object."""
        return wrap(self.net.Clone())

    def _ensure_mutable(self) -> None:
        """Do nothing; :class:`FrozenDistribution` raises here. Members that change the distribution call it first."""

    def _value_key(self) -> tuple:
        """Return the value compared by ``==``: .NET type, parameter bits, ``Base`` bits and, where needed, the XML."""
        net = self.net
        has_base, needs_xml = _state_info(net)
        key: tuple = (str(net.GetType().FullName), to_numpy_doubles(net.GetParameters).tobytes())
        if has_base:
            key += (np.float64(float(net.Base)).tobytes(),)
        if needs_xml:
            key += (str(net.ToXElement().ToString()),)
        return key

    def __eq__(self, other: object) -> bool:
        """Return whether ``other`` is a distribution of the same family with bitwise-equal parameters.

        The family is the wrapped .NET type; parameters compare by their bits (so ``0.0`` and ``-0.0``
        differ and a NaN equals the same NaN), as does ``Base`` where the family has one. Families whose
        parameter vector is not their whole state (mixtures, empirical tables, kernel densities,
        competing risks) also compare their ``ToXElement`` text. Frozen and mutable wrappers of equal
        values are equal.
        """
        if not isinstance(other, UnivariateDistribution):
            return NotImplemented
        return self is other or self._value_key() == other._value_key()

    __hash__ = None  # type: ignore[assignment]  # mutable values are unhashable; see frozen()

    def frozen(self) -> FrozenDistribution:
        """Return a read-only, hashable copy of this distribution.

        The copy wraps its own ``Clone`` of the .NET object, so later changes to this distribution do
        not reach it. It is an instance of :class:`FrozenDistribution` and of this distribution's
        class, so the family's members stay available.

        Returns
        -------
        FrozenDistribution
            The frozen copy; ``set_parameters``, ``estimate`` and attribute assignment raise, and
            ``copy.deepcopy`` returns it unchanged.
        """
        cls = _frozen_class(type(self))
        instance = cls.__new__(cls)
        instance.__dict__["net"] = unwrap(self.net.Clone())
        return instance

    def _scalar(self, name: str) -> Any:
        """Return the ``double`` overload of a scalar member, bound once per instance (see ``_util.bind``)."""
        cache = self.__dict__.setdefault("_bound", {})
        member = cache.get(name)
        if member is None:
            member = cache[name] = bind(getattr(self.net, name))
        return member

    def _vector(self, name: str) -> Any:
        """Return the ``IList<double>`` overload of an array member, bound once per instance (``_util.bind_list``)."""
        cache = self.__dict__.setdefault("_bound_lists", {})
        member = cache.get(name)
        if member is None:
            member = cache[name] = bind_list(getattr(self.net, name))
        return member

    def _evaluate(self, name: str, values: Any) -> Any:
        """Evaluate a scalar member on a scalar or an array, through its array overload when Numerics has one.

        A member without an ``IList<double>`` overload (``HF``) evaluates an array element by element
        with the bound scalar overload. Both give the distribution's own scalar values, bit for bit.
        """
        cache = self.__dict__.setdefault("_optional_lists", {})
        if name not in cache:
            cache[name] = list_overload(getattr(self.net, name))
        vector = cache[name]
        if vector is None:
            return loop(values, self._scalar(name))
        return evaluate(values, self._scalar(name), vector)

    def __repr__(self) -> str:
        return f"{type(self).__name__}({', '.join(repr(float(p)) for p in self.parameters)})"

    @property
    def display_name(self) -> str:
        """The Numerics display name."""
        return str(self.net.DisplayName)

    @property
    def parameters(self) -> np.ndarray:
        """The parameter values in Numerics order (``GetParameters``)."""
        return to_numpy_doubles(self.net.GetParameters)

    @property
    def parameters_valid(self) -> bool:
        """Whether the current parameters are valid (``ParametersValid``)."""
        return bool(self.net.ParametersValid)

    @property
    def mean(self) -> float:
        """The distribution mean."""
        return float(self.net.Mean)

    @property
    def standard_deviation(self) -> float:
        """The distribution standard deviation."""
        return float(self.net.StandardDeviation)

    @property
    def skewness(self) -> float:
        """The distribution skewness."""
        return float(self.net.Skewness)

    @property
    def kurtosis(self) -> float:
        """The distribution kurtosis."""
        return float(self.net.Kurtosis)

    @property
    def minimum(self) -> float:
        """The lower end of the support."""
        return float(self.net.Minimum)

    @property
    def maximum(self) -> float:
        """The upper end of the support."""
        return float(self.net.Maximum)

    def pdf(self, x: Any) -> Any:
        """Evaluate the probability density (or mass) function.

        Parameters
        ----------
        x : float or array_like
            Values at which to evaluate.

        Returns
        -------
        float or numpy.ndarray
            Densities, shaped like ``x``.
        """
        return evaluate(x, self._scalar("PDF"), self._vector("PDF"))

    def cdf(self, x: Any) -> Any:
        """Evaluate the cumulative distribution function (non-exceedance probability).

        Parameters
        ----------
        x : float or array_like
            Values at which to evaluate.

        Returns
        -------
        float or numpy.ndarray
            Non-exceedance probabilities, shaped like ``x``.
        """
        return evaluate(x, self._scalar("CDF"), self._vector("CDF"))

    def ccdf(self, x: Any) -> Any:
        """Evaluate the complementary CDF (exceedance probability) through ``CCDF``.

        An array crosses once into ``CCDF(IList<double>)``; every value is the distribution's own ``CCDF``.

        Parameters
        ----------
        x : float or array_like
            Values at which to evaluate.

        Returns
        -------
        float or numpy.ndarray
            Exceedance probabilities, shaped like ``x``.
        """
        return self._evaluate("CCDF", x)

    def inverse_cdf(self, probability: Any) -> Any:
        """Evaluate the quantile function at non-exceedance probabilities.

        Parameters
        ----------
        probability : float or array_like
            Non-exceedance probabilities in [0, 1].

        Returns
        -------
        float or numpy.ndarray
            Quantiles, shaped like ``probability``.
        """
        return evaluate(probability, self._scalar("InverseCDF"), self._vector("InverseCDF"))

    def log_pdf(self, x: Any) -> Any:
        """Evaluate the natural logarithm of the density through ``LogPDF``.

        An array crosses once into ``LogPDF(IList<double>)``; every value is the distribution's own ``LogPDF``.

        Parameters
        ----------
        x : float or array_like
            Values at which to evaluate.

        Returns
        -------
        float or numpy.ndarray
            Log densities, shaped like ``x``.
        """
        return self._evaluate("LogPDF", x)

    def log_cdf(self, x: Any) -> Any:
        """Evaluate the natural logarithm of the CDF through ``LogCDF``.

        An array crosses once into ``LogCDF(IList<double>)``.

        Parameters
        ----------
        x : float or array_like
            Values at which to evaluate.

        Returns
        -------
        float or numpy.ndarray
            Log non-exceedance probabilities, shaped like ``x``.
        """
        return self._evaluate("LogCDF", x)

    def log_ccdf(self, x: Any) -> Any:
        """Evaluate the natural logarithm of the complementary CDF through ``LogCCDF``.

        An array crosses once into ``LogCCDF(IList<double>)``.

        Parameters
        ----------
        x : float or array_like
            Values at which to evaluate.

        Returns
        -------
        float or numpy.ndarray
            Log exceedance probabilities, shaped like ``x``.
        """
        return self._evaluate("LogCCDF", x)

    type_name = _FamilyAttribute(
        _class_type_name, lambda self: str(self.net.Type),
        """The ``UnivariateDistributionType`` member name, for example ``"GeneralizedExtremeValue"``.

        Available on the typed classes without an instance (``Gumbel.type_name``); None on the
        generic :class:`UnivariateDistribution` class.
        """)

    parameter_names = _FamilyAttribute(
        _class_parameter_names, lambda self: tuple(str(name) for name in self.net.GetParameterPropertyNames),
        """The parameter property names in ``GetParameters`` order (``GetParameterPropertyNames``).

        Available on the typed classes without an instance, from the family's default instance;
        None on the generic :class:`UnivariateDistribution` class and on a family Numerics cannot
        build without parameters (``Mixture``).
        """)

    @property
    def number_of_parameters(self) -> int:
        """The number of parameters (``NumberOfParameters``)."""
        return int(self.net.NumberOfParameters)

    def set_parameters(self, parameters: Any) -> None:
        """Set the parameters in place through ``SetParameters(IList<double>)``.

        Invalid values are accepted and recorded, as in Numerics: ``parameters_valid``
        becomes False and the evaluation members raise.

        Parameters
        ----------
        parameters : array_like
            The parameters in Numerics order.

        Raises
        ------
        TypeError
            If the distribution is frozen.
        """
        self._ensure_mutable()
        self.net.SetParameters(to_net_doubles(parameters))

    @_FamilyMethod
    def validate_parameters(self, parameters: Any) -> str | None:
        """Return why a parameter vector is invalid for this family, or None when it is valid.

        Callable on a typed class without an instance (``GeneralizedExtremeValue.validate_parameters(p)``,
        checked on the family's default instance) or on an instance. The count of values is checked
        first; then ``ValidateParameters(IList<double>, false)`` runs and the message of the exception
        Numerics builds, without throwing it, is returned.

        Parameters
        ----------
        parameters : array_like
            A parameter vector in Numerics order.

        Returns
        -------
        str or None
            The Numerics message, or None for a valid vector.

        Raises
        ------
        ValueError
            If the vector does not have the family's number of parameters (families described by
            a parameter vector only).
        TypeError
            If called on the generic :class:`UnivariateDistribution` class.
        NotImplementedError
            If called on the class of a family Numerics cannot build without parameters (``Mixture``).
        """
        if isinstance(self, type):
            type_name = _class_type_name(self)
            if type_name is None:
                raise TypeError("validate_parameters needs a typed family class or an instance.")
            net = _prototype(type_name)
        else:
            net = self.net
        values = _check_arity(net, parameters)
        error = net.ValidateParameters(to_net_doubles(values), False)
        return None if error is None else str(error.Message)

    @property
    def has_parameter_covariance(self) -> bool:
        """Whether the family implements ``IStandardError`` (an asymptotic parameter covariance)."""
        nd = _nd()
        import clr  # noqa: PLC0415

        return bool(clr.GetClrType(nd.IStandardError).IsInstanceOfType(self.net))

    def parameter_covariance(self, sample_size: int,
                             method: enums.EstimationMethod | str = enums.EstimationMethod.MAXIMUM_LIKELIHOOD,
                             ) -> np.ndarray:
        """Return the asymptotic parameter covariance through ``IStandardError.ParameterCovariance``.

        Parameters
        ----------
        sample_size : int
            The sample size n; the covariance scales as 1/n.
        method : EstimationMethod or str, optional
            The estimation method whose covariance is returned; default maximum likelihood.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(number_of_parameters, number_of_parameters)``.

        Raises
        ------
        NotImplementedError
            If the distribution does not implement ``IStandardError``.
        """
        nd = _nd()
        if not self.has_parameter_covariance:
            raise NotImplementedError(f"{self.display_name} exposes no parameter covariance.")
        covariance = nd.IStandardError(self.net).ParameterCovariance(int(sample_size),
                                                                     net_enum(nd.ParameterEstimationMethod, method))
        return to_numpy_matrix(covariance)

    def parameter_set_quantiles(self, parameter_sets: Any, probability: float, *, validate: bool = False) -> np.ndarray:
        """Return one quantile per parameter set, each from this family at that set.

        A clone of this distribution receives every row through ``SetParameters`` and
        returns ``InverseCDF(probability)``: the probe pattern a .NET caller uses for many
        parameter sets, with the scalar overloads bound once. The rows cross one at a time
        and are not validated; :meth:`parameter_set_quantile_table` passes the whole matrix
        to ``ParameterSetQuantiles`` in one call.

        Parameters
        ----------
        parameter_sets : array_like
            Array of shape ``(count, number_of_parameters)``.
        probability : float
            The non-exceedance probability.
        validate : bool, optional
            When True, each row is checked with ``ValidateParameters(row, false)`` before it is applied,
            and the first invalid row raises its own Numerics message with its index. Default False:
            rows are not validated and an invalid row raises from ``InverseCDF`` without an index.

        Returns
        -------
        numpy.ndarray
            ``count`` quantiles.

        Raises
        ------
        ValueError
            If the rows do not have ``number_of_parameters`` columns, or ``validate`` is True and a
            row is invalid (``"Parameter set <index>: <Numerics message>"``).
        """
        rows = np.atleast_2d(np.asarray(parameter_sets, dtype=np.float64))
        arity = self.number_of_parameters
        if rows.shape[1] != arity:
            raise ValueError(f"Expected {arity} parameters per row; got {rows.shape[1]}.")
        if validate:
            self._raise_first_invalid(rows, self.parameter_set_validity(rows))
        probe = self.net.Clone()
        quantile = bind(probe.InverseCDF)
        p = float(probability)
        out = np.empty(rows.shape[0], dtype=np.float64)
        setter = _typed_member(probe, "SetParameters", arity) if arity > 1 else None
        for index, row in enumerate(rows.tolist()):
            if setter is not None:
                setter(probe, *row)
            else:
                probe.SetParameters(to_net_doubles(row))
            out[index] = quantile(p)
        return out

    def parameter_set_quantile_table(self, parameter_sets: Any, probabilities: Any, *,
                                     validate: bool = False) -> tuple[np.ndarray, np.ndarray]:
        """Return the quantiles of many probabilities under many parameter sets, in one call into Numerics.

        Calls ``UnivariateDistributionBase.ParameterSetQuantiles(double[,], IList<double>, out bool[])``:
        the parameter matrix crosses into .NET once, every row is validated with
        ``ValidateParameters(row, false)``, a valid row is applied to one private clone with
        ``SetParameters(row)`` and evaluated with ``InverseCDF`` at each probability, and the
        quantile matrix and the validity flags cross back once. The values are bitwise those
        of the per-row loop; the rows of invalid parameter sets are NaN.

        Parameters
        ----------
        parameter_sets : array_like
            Array of shape ``(count, number_of_parameters)``.
        probabilities : array_like
            The non-exceedance probabilities, one-dimensional.
        validate : bool, optional
            When True, an invalid row raises its own Numerics message with its index instead of
            returning a NaN row. Default False.

        Returns
        -------
        quantiles : numpy.ndarray
            Array of shape ``(count, len(probabilities))``; NaN rows where the set is invalid.
        valid : numpy.ndarray
            Boolean array of length ``count``, True where ``ValidateParameters(row, false)`` is null.

        Raises
        ------
        ValueError
            If the rows do not have ``number_of_parameters`` columns, ``probabilities`` is not
            one-dimensional, or ``validate`` is True and a row is invalid
            (``"Parameter set <index>: <Numerics message>"``, for the first invalid row).
        """
        rows = np.atleast_2d(np.asarray(parameter_sets, dtype=np.float64))
        arity = self.number_of_parameters
        if rows.shape[1] != arity:
            raise ValueError(f"Expected {arity} parameters per row; got {rows.shape[1]}.")
        p = np.asarray(probabilities, dtype=np.float64)
        if p.ndim != 1:
            raise ValueError(f"Expected one-dimensional probabilities; got shape {p.shape}.")
        quantiles, valid = self.net.ParameterSetQuantiles(to_net_matrix(rows), to_net_doubles(p), None)
        valid = to_numpy_bools(valid)
        if validate:
            self._raise_first_invalid(rows, valid)
        return to_numpy_matrix(quantiles), valid

    def _raise_first_invalid(self, rows: np.ndarray, valid: np.ndarray) -> None:
        """Raise ValueError for the first row marked invalid, with its index and its own Numerics message."""
        invalid = np.flatnonzero(~valid)
        if invalid.size == 0:
            return
        index = int(invalid[0])
        error = self.net.ValidateParameters(to_net_doubles(rows[index]), False)
        message = "invalid parameters" if error is None else str(error.Message)
        raise ValueError(f"Parameter set {index}: {message}") from (None if error is None else translate(error))

    def parameter_set_validity(self, parameter_sets: Any) -> np.ndarray:
        """Return which parameter sets this family accepts (``ValidateParameters(row, false)`` is null).

        The matrix crosses once, into ``ParameterSetQuantiles`` with no probabilities, whose validity
        flags are this test row by row.

        Parameters
        ----------
        parameter_sets : array_like
            Array of shape ``(count, number_of_parameters)``.

        Returns
        -------
        numpy.ndarray
            Boolean array of length ``count``, True where the row is valid.

        Raises
        ------
        ValueError
            If the rows do not have ``number_of_parameters`` columns.
        """
        rows = np.atleast_2d(np.asarray(parameter_sets, dtype=np.float64))
        arity = self.number_of_parameters
        if rows.shape[1] != arity:
            raise ValueError(f"Expected {arity} parameters per row; got {rows.shape[1]}.")
        _, valid = self.net.ParameterSetQuantiles(to_net_matrix(rows), to_net_doubles(np.empty(0)), None)
        return to_numpy_bools(valid)

    def parameter_set_distributions(self, parameter_sets: Any, valid: Any = None) -> Any:
        """Build one distribution of this family per parameter set, as a .NET array.

        Each member is a ``Clone`` of this distribution with ``SetParameters(row)``, so fixed
        settings that are not parameters (a ``LogNormal`` base) carry over. The result stays
        on the .NET side for the members that take ``UnivariateDistributionBase[]``
        (``UncertaintyAnalysisResults``, ``BootstrapAnalysis.Quantiles``); rows marked False
        in ``valid`` are null entries, which those members treat as missing.

        Parameters
        ----------
        parameter_sets : array_like
            Array of shape ``(count, number_of_parameters)``.
        valid : array_like of bool, optional
            Which rows to build; None builds every row.

        Returns
        -------
        System.Array[Numerics.Distributions.UnivariateDistributionBase]
            ``count`` members, null where ``valid`` is False.

        Raises
        ------
        ValueError
            If the rows do not have ``number_of_parameters`` columns or ``valid`` has another length.
        """
        rows = np.atleast_2d(np.asarray(parameter_sets, dtype=np.float64))
        arity = self.number_of_parameters
        if rows.shape[1] != arity:
            raise ValueError(f"Expected {arity} parameters per row; got {rows.shape[1]}.")
        keep = np.ones(rows.shape[0], dtype=bool) if valid is None else np.asarray(valid, dtype=bool)
        if keep.shape != (rows.shape[0],):
            raise ValueError("valid must hold one flag per row.")
        nd = _nd()
        import System  # noqa: PLC0415 - the runtime is running (self.net exists)

        members = System.Array.CreateInstance(nd.UnivariateDistributionBase, rows.shape[0])
        clone = self.net.Clone
        setter = _typed_member(self.net, "SetParameters", arity) if arity > 1 else None
        for index, row in enumerate(rows.tolist()):
            if not keep[index]:
                continue
            member = clone()
            if setter is not None:
                setter(member, *row)
            else:
                member.SetParameters(to_net_doubles(row))
            members[index] = member
        return members

    def to_xml(self) -> str:
        """Return the distribution's ``ToXElement()`` as XML text.

        Returns
        -------
        str
            The element as .NET writes it (``XElement.ToString``).
        """
        return str(self.net.ToXElement().ToString())

    def generate_random_values(self, sample_size: int, seed: int = -1) -> np.ndarray:
        """Draw a Monte Carlo sample through ``GenerateRandomValues``.

        Parameters
        ----------
        sample_size : int
            Number of values.
        seed : int, optional
            Mersenne Twister seed; zero or negative seeds from the clock (not
            reproducible).

        Returns
        -------
        numpy.ndarray
            The sample.
        """
        return to_numpy_doubles(self.net.GenerateRandomValues(int(sample_size), int(seed)))

    def log_likelihood(self, sample: Any) -> float:
        """Return the log-likelihood of a sample (``LogLikelihood(IList<double>)``).

        Parameters
        ----------
        sample : array_like
            Observations.

        Returns
        -------
        float
            The log-likelihood.
        """
        return float(self.net.LogLikelihood(to_net_doubles(sample)))

    def estimate(self, sample: Any, method: enums.EstimationMethod | str = enums.EstimationMethod.MAXIMUM_LIKELIHOOD,
                 fixed_parameters: Any = None) -> None:
        """Fit the parameters in place through ``IEstimation.Estimate``.

        With ``fixed_parameters`` the fit holds those parameters at their values through
        ``IFixedParameterEstimation.Estimate(sample, method, fixed)``, for example the generalized
        Pareto scale and shape by maximum likelihood at a known location.

        Parameters
        ----------
        sample : array_like
            Observations.
        method : EstimationMethod or str, optional
            Estimation method; default maximum likelihood.
        fixed_parameters : array_like, optional
            One entry per parameter in Numerics order: a finite value holds that parameter fixed and
            NaN marks a parameter to estimate. None (the default) estimates every parameter.

        Raises
        ------
        TypeError
            If the distribution does not implement ``IEstimation``, or ``fixed_parameters`` is given and
            it does not implement ``IFixedParameterEstimation``.

        Notes
        -----
        A frozen distribution raises ``TypeError``.
        """
        self._ensure_mutable()
        nd = _nd()
        import clr  # noqa: PLC0415

        if not clr.GetClrType(nd.IEstimation).IsInstanceOfType(self.net):
            raise TypeError(f"{self.display_name} does not support parameter estimation.")
        kind = net_enum(nd.ParameterEstimationMethod, method)
        if fixed_parameters is None:
            nd.IEstimation(self.net).Estimate(to_net_doubles(sample), kind)
            return
        fixed = nd.IFixedParameterEstimation
        if not clr.GetClrType(fixed).IsInstanceOfType(self.net):
            raise TypeError(f"{self.display_name} does not support estimation with fixed parameters.")
        fixed(self.net).Estimate(to_net_doubles(sample), kind, to_net_doubles(fixed_parameters))

    def clone(self) -> UnivariateDistribution:
        """Return a deep copy through ``Clone``.

        Returns
        -------
        UnivariateDistribution
            A wrapper of the cloned .NET object, of the same Python class when possible.
        """
        return wrap(self.net.Clone())

    @property
    def median(self) -> float:
        """The distribution median."""
        return float(self.net.Median)

    @property
    def mode(self) -> float:
        """The distribution mode."""
        return float(self.net.Mode)

    @property
    def variance(self) -> float:
        """The distribution variance (``Variance``, the square of ``StandardDeviation``)."""
        return float(self.net.Variance)

    @property
    def coefficient_of_variation(self) -> float:
        """The coefficient of variation, standard deviation over mean (``CoefficientOfVariation``)."""
        return float(self.net.CoefficientOfVariation)

    @property
    def short_display_name(self) -> str:
        """The Numerics short display name, for example ``"GEV"``."""
        return str(self.net.ShortDisplayName)

    @property
    def display_label(self) -> str:
        """The Numerics label: the short display name followed by the parameter values."""
        return str(self.net.DisplayLabel)

    @property
    def parameter_display_names(self) -> tuple[str, ...]:
        """The parameter names for display, for example ``"Location (ξ)"`` (``ParameterNames``)."""
        return tuple(str(name) for name in self.net.ParameterNames)

    @property
    def parameter_names_short_form(self) -> tuple[str, ...]:
        """The short parameter symbols, for example ``"ξ"`` (``ParameterNamesShortForm``)."""
        return tuple(str(name) for name in self.net.ParameterNamesShortForm)

    @property
    def parameters_to_string(self) -> list[tuple[str, str]]:
        """The (name, value) text rows Numerics shows for the parameters (``ParametersToString``)."""
        table = self.net.ParametersToString
        return [(str(table[i, 0]), str(table[i, 1])) for i in range(int(table.GetLength(0)))]

    @property
    def minimum_of_parameters(self) -> np.ndarray:
        """The lower bound of each parameter, in Numerics order (``MinimumOfParameters``)."""
        return to_numpy_doubles(self.net.MinimumOfParameters)

    @property
    def maximum_of_parameters(self) -> np.ndarray:
        """The upper bound of each parameter, in Numerics order (``MaximumOfParameters``)."""
        return to_numpy_doubles(self.net.MaximumOfParameters)

    def hf(self, x: Any) -> Any:
        """Evaluate the hazard function, the density over the survival probability (``HF``).

        Parameters
        ----------
        x : float or array_like
            Values at which to evaluate.

        Returns
        -------
        float or numpy.ndarray
            Hazard rates, shaped like ``x``; Numerics has no array overload, so an array is evaluated
            element by element with the bound scalar member.
        """
        return self._evaluate("HF", x)

    def inverse_ccdf(self, probability: Any) -> Any:
        """Evaluate the survival quantile: x with ``ccdf(x) = probability`` (``InverseCCDF``).

        Families with a closed form in the exceedance probability (``GeneralizedPareto``) never form
        ``1 - probability``, which keeps the precision of small exceedance probabilities; the others
        return ``InverseCDF(1 - probability)``. An array crosses once into ``InverseCCDF(IList<double>)``.

        Parameters
        ----------
        probability : float or array_like
            Exceedance probabilities in [0, 1].

        Returns
        -------
        float or numpy.ndarray
            Quantiles, shaped like ``probability``.
        """
        return self._evaluate("InverseCCDF", probability)

    def central_moments(self, tolerance: float = 1e-8, steps: int | None = None) -> np.ndarray:
        """Return the mean, standard deviation, skewness and kurtosis by numerical integration (``CentralMoments``).

        With ``steps`` None the adaptive Gauss-Kronrod overload ``CentralMoments(double tolerance)``
        runs over ``[InverseCDF(1e-16), InverseCDF(1 - 1e-16)]``; with ``steps`` the fixed-step
        bin-expectation overload ``CentralMoments(int steps)`` runs instead and ``tolerance`` is unused.

        Parameters
        ----------
        tolerance : float, optional
            Relative tolerance of the adaptive integration; default 1e-8.
        steps : int, optional
            Number of bins of the fixed-step overload (Numerics' default is 300).

        Returns
        -------
        numpy.ndarray
            ``[mean, standard deviation, skewness, kurtosis]``; a moment whose integration fails is NaN.
        """
        import System  # noqa: PLC0415 - the runtime is running

        member = self.net.CentralMoments
        if steps is None:
            return to_numpy_doubles(member.Overloads[System.Double](float(tolerance)))
        return to_numpy_doubles(member.Overloads[System.Int32](int(steps)))

    def conditional_moments(self, a: float, b: float) -> np.ndarray:
        """Return the first four moments of X restricted to [a, b] by 300-bin integration (``ConditionalMoments``).

        Parameters
        ----------
        a : float
            Lower limit.
        b : float
            Upper limit.

        Returns
        -------
        numpy.ndarray
            ``[E[X | a < X < b], m2, m3, m4]``, where ``mk`` is the conditional expectation of
            ``(X - mean)**k`` about the unconditional distribution mean, as Numerics computes them;
            NaN when ``a >= b`` or the interval has no probability. They are not the central moments
            of the truncated distribution: its variance is ``m2 - (m1 - mean)**2``.
        """
        return to_numpy_doubles(self.net.ConditionalMoments(float(a), float(b)))

    def conditional_mean(self, a: float, b: float, tolerance: float = 1e-8) -> float:
        """Return E[X | a < X < b] by adaptive Gauss-Kronrod integration (``ConditionalMean``).

        Parameters
        ----------
        a : float
            Lower limit.
        b : float
            Upper limit.
        tolerance : float, optional
            Relative integration tolerance; default 1e-8.

        Returns
        -------
        float
            The integral of ``x f(x)`` over [a, b] divided by ``F(b) - F(a)``; NaN if the integration fails.
        """
        return float(self.net.ConditionalMean(float(a), float(b), float(tolerance)))

    def conditional_variance(self, a: float, b: float, mean: float = 0.0, tolerance: float = 1e-8) -> float:
        """Return E[(X - mean)^2 | a < X < b] by adaptive integration (``ConditionalVariance``).

        Parameters
        ----------
        a : float
            Lower limit.
        b : float
            Upper limit.
        mean : float, optional
            The centre of the moment; default 0 (pass :meth:`conditional_mean` for the central moment).
        tolerance : float, optional
            Relative integration tolerance; default 1e-8.

        Returns
        -------
        float
            The conditional second moment about ``mean``.
        """
        return float(self.net.ConditionalVariance(float(a), float(b), float(mean), float(tolerance)))

    def conditional_skewness(self, a: float, b: float, mean: float = 0.0, tolerance: float = 1e-8) -> float:
        """Return E[(X - mean)^3 | a < X < b] by adaptive integration (``ConditionalSkewness``).

        Parameters
        ----------
        a : float
            Lower limit.
        b : float
            Upper limit.
        mean : float, optional
            The centre of the moment; default 0.
        tolerance : float, optional
            Relative integration tolerance; default 1e-8.

        Returns
        -------
        float
            The conditional third moment about ``mean`` (not standardized).
        """
        return float(self.net.ConditionalSkewness(float(a), float(b), float(mean), float(tolerance)))

    def conditional_kurtosis(self, a: float, b: float, mean: float = 0.0, tolerance: float = 1e-8) -> float:
        """Return E[(X - mean)^4 | a < X < b] by adaptive integration (``ConditionalKurtosis``).

        Parameters
        ----------
        a : float
            Lower limit.
        b : float
            Upper limit.
        mean : float, optional
            The centre of the moment; default 0.
        tolerance : float, optional
            Relative integration tolerance; default 1e-8.

        Returns
        -------
        float
            The conditional fourth moment about ``mean`` (not standardized).
        """
        return float(self.net.ConditionalKurtosis(float(a), float(b), float(mean), float(tolerance)))

    def conditional_expected_value(self, alpha: float, tolerance: float = 1e-8) -> float:
        """Return the mean of the upper tail beyond the ``alpha`` quantile (``ConditionalExpectedValue``).

        Parameters
        ----------
        alpha : float
            The non-exceedance level; the tail is ``[InverseCDF(alpha), InverseCDF(1 - 1e-16)]``.
        tolerance : float, optional
            Relative integration tolerance; default 1e-8.

        Returns
        -------
        float
            The integral of ``x f(x)`` over the tail divided by ``1 - alpha`` (the conditional value at risk).
        """
        return float(self.net.ConditionalExpectedValue(float(alpha), float(tolerance)))

    def create_pdf_graph(self, bins: Any = None) -> np.ndarray:
        """Return a (x, density) table for plotting (``CreatePDFGraph``).

        Parameters
        ----------
        bins : pynumerics.sampling.StratificationBins or sequence of them, optional
            X-value bins whose midpoints are the table points; None, the default, lets Numerics
            choose the bounds.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(n, 2)``: x and the density.
        """
        member = self.net.CreatePDFGraph
        return to_numpy_matrix(member() if bins is None else member(_bins_net(bins)))

    def create_cdf_graph(self, bins: Any = None) -> np.ndarray:
        """Return a (x, non-exceedance probability) table for plotting (``CreateCDFGraph``).

        Parameters
        ----------
        bins : pynumerics.sampling.StratificationBins or sequence of them, optional
            X-value bins; None, the default, lets Numerics choose the bounds.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(n, 2)``: x and the CDF.
        """
        member = self.net.CreateCDFGraph
        return to_numpy_matrix(member() if bins is None else member(_bins_net(bins)))

    def create_inverse_cdf_graph(self, bins: Any = None) -> np.ndarray:
        """Return a (probability, quantile) table for plotting (``CreateInverseCDFGraph``).

        Parameters
        ----------
        bins : pynumerics.sampling.StratificationBins or sequence of them, optional
            Probability bins; None, the default, lets Numerics choose the probabilities.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(n, 2)``: the probability and its quantile.
        """
        member = self.net.CreateInverseCDFGraph
        return to_numpy_matrix(member() if bins is None else member(_bins_net(bins)))

    def log_likelihood_left_censored(self, threshold: float, number_below: int) -> float:
        """Return the log-likelihood of ``number_below`` values known only to lie below ``threshold``.

        Parameters
        ----------
        threshold : float
            The censoring threshold.
        number_below : int
            The count of censored values.

        Returns
        -------
        float
            ``number_below * log F(threshold)`` as Numerics computes it (``LogLikelihood_LeftCensored``).
        """
        return float(self.net.LogLikelihood_LeftCensored(float(threshold), int(number_below)))

    def log_likelihood_right_censored(self, threshold: float, number_above: int) -> float:
        """Return the log-likelihood of ``number_above`` values known only to exceed ``threshold``.

        Parameters
        ----------
        threshold : float
            The censoring threshold.
        number_above : int
            The count of censored values.

        Returns
        -------
        float
            ``number_above * log(1 - F(threshold))`` as Numerics computes it (``LogLikelihood_RightCensored``).
        """
        return float(self.net.LogLikelihood_RightCensored(float(threshold), int(number_above)))

    def log_likelihood_intervals(self, lower: float, upper: float) -> float:
        """Return the log-likelihood of one value known to lie in [lower, upper] (``LogLikelihood_Intervals``).

        Parameters
        ----------
        lower : float
            Lower limit of the interval.
        upper : float
            Upper limit of the interval.

        Returns
        -------
        float
            ``log(F(upper) - F(lower))`` as Numerics computes it.
        """
        return float(self.net.LogLikelihood_Intervals(float(lower), float(upper)))

    def compare_to(self, other: UnivariateDistribution | None) -> int:
        """Order two distributions by their means (``CompareTo``).

        Parameters
        ----------
        other : UnivariateDistribution or None
            The distribution to compare with; None orders below every distribution.

        Returns
        -------
        int
            Negative, zero or positive as this mean is below, equal to or above the other's.
        """
        return int(self.net.CompareTo(None if other is None else other.net))

    def bootstrap(self, method: enums.EstimationMethod | str, sample_size: int,
                  seed: int = -1) -> UnivariateDistribution:
        """Draw a sample from this distribution and refit the family to it (``IBootstrappable.Bootstrap``).

        Parameters
        ----------
        method : EstimationMethod or str
            The estimation method of the refit.
        sample_size : int
            The size of the drawn sample.
        seed : int, optional
            Mersenne Twister seed; zero or negative seeds from the clock (not reproducible).

        Returns
        -------
        UnivariateDistribution
            The refitted distribution (a new object; this one is unchanged).

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IBootstrappable``.
        """
        nd = _nd()
        replicate = _interface(self.net, "IBootstrappable").Bootstrap(
            net_enum(nd.ParameterEstimationMethod, method), int(sample_size), int(seed))
        return wrap(replicate)

    def mle(self, sample: Any) -> np.ndarray:
        """Return the maximum likelihood parameters of a sample, without changing this distribution.

        Parameters
        ----------
        sample : array_like
            Observations.

        Returns
        -------
        numpy.ndarray
            The parameters in Numerics order (``IMaximumLikelihoodEstimation.MLE``).

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IMaximumLikelihoodEstimation``.
        """
        return to_numpy_doubles(_interface(self.net, "IMaximumLikelihoodEstimation").MLE(to_net_doubles(sample)))

    def parameter_constraints(self, sample: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the starting values and bounds the maximum likelihood search uses for a sample.

        Parameters
        ----------
        sample : array_like
            Observations.

        Returns
        -------
        initial : numpy.ndarray
            The starting parameters (``IMaximumLikelihoodEstimation.GetParameterConstraints``).
        lower : numpy.ndarray
            The lower bound of each parameter.
        upper : numpy.ndarray
            The upper bound of each parameter.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IMaximumLikelihoodEstimation``.
        """
        result = _interface(self.net, "IMaximumLikelihoodEstimation").GetParameterConstraints(to_net_doubles(sample))
        return to_numpy_doubles(result.Item1), to_numpy_doubles(result.Item2), to_numpy_doubles(result.Item3)

    @_FamilyMethod
    def parameters_from_moments(self, moments: Any) -> np.ndarray:
        """Return the parameters matching product moments (``IMomentEstimation.ParametersFromMoments``).

        Callable on a typed class without an instance (on the family's default instance) or on an
        instance; families with a setting that is not a parameter (a ``LogNormal`` base) use the
        instance's setting.

        Parameters
        ----------
        moments : array_like
            ``[mean, standard deviation, skewness, kurtosis]`` (the leading entries the family needs).

        Returns
        -------
        numpy.ndarray
            The parameters in Numerics order.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IMomentEstimation``.
        """
        net = _interface(_net_of(self), "IMomentEstimation")
        return to_numpy_doubles(net.ParametersFromMoments(to_net_doubles(moments)))

    @_FamilyMethod
    def moments_from_parameters(self, parameters: Any) -> np.ndarray:
        """Return the product moments of a parameter vector (``IMomentEstimation.MomentsFromParameters``).

        Callable on a typed class without an instance or on an instance.

        Parameters
        ----------
        parameters : array_like
            The parameters in Numerics order.

        Returns
        -------
        numpy.ndarray
            ``[mean, standard deviation, skewness, kurtosis]``.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IMomentEstimation``.
        """
        net = _interface(_net_of(self), "IMomentEstimation")
        return to_numpy_doubles(net.MomentsFromParameters(to_net_doubles(parameters)))

    @_FamilyMethod
    def parameters_from_linear_moments(self, moments: Any) -> np.ndarray:
        """Return the parameters matching L-moments (``ILinearMomentEstimation.ParametersFromLinearMoments``).

        Callable on a typed class without an instance
        (``GeneralizedExtremeValue.parameters_from_linear_moments(l)``) or on an instance.

        Parameters
        ----------
        moments : array_like
            ``[λ1, λ2, τ3, τ4]``: the first two L-moments and the L-skewness and L-kurtosis ratios,
            as :func:`pynumerics.statistics.linear_moments` returns them.

        Returns
        -------
        numpy.ndarray
            The parameters in Numerics order.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``ILinearMomentEstimation``.
        """
        net = _interface(_net_of(self), "ILinearMomentEstimation")
        return to_numpy_doubles(net.ParametersFromLinearMoments(to_net_doubles(moments)))

    @_FamilyMethod
    def linear_moments_from_parameters(self, parameters: Any) -> np.ndarray:
        """Return the L-moments of a parameter vector (``ILinearMomentEstimation.LinearMomentsFromParameters``).

        Callable on a typed class without an instance or on an instance.

        Parameters
        ----------
        parameters : array_like
            The parameters in Numerics order.

        Returns
        -------
        numpy.ndarray
            ``[λ1, λ2, τ3, τ4]``.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``ILinearMomentEstimation``.
        """
        net = _interface(_net_of(self), "ILinearMomentEstimation")
        return to_numpy_doubles(net.LinearMomentsFromParameters(to_net_doubles(parameters)))

    def quantile_variance(self, probability: float, sample_size: int,
                          method: enums.EstimationMethod | str = enums.EstimationMethod.MAXIMUM_LIKELIHOOD) -> float:
        """Return the asymptotic variance of the quantile estimate (``IStandardError.QuantileVariance``).

        Parameters
        ----------
        probability : float
            The non-exceedance probability.
        sample_size : int
            The sample size n.
        method : EstimationMethod or str, optional
            The estimation method; default maximum likelihood.

        Returns
        -------
        float
            The variance of the estimated ``inverse_cdf(probability)``.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IStandardError``.
        """
        net = _interface(self.net, "IStandardError")
        return float(net.QuantileVariance(float(probability), int(sample_size),
                                          net_enum(_nd().ParameterEstimationMethod, method)))

    def quantile_gradient(self, probability: float) -> np.ndarray:
        """Return the gradient of the quantile with respect to the parameters (``IStandardError.QuantileGradient``).

        Parameters
        ----------
        probability : float
            The non-exceedance probability.

        Returns
        -------
        numpy.ndarray
            The partial derivatives of ``inverse_cdf(probability)``, in parameter order.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IStandardError``.
        """
        return to_numpy_doubles(_interface(self.net, "IStandardError").QuantileGradient(float(probability)))

    def quantile_jacobian(self, probabilities: Any) -> tuple[np.ndarray, float]:
        """Return the Jacobian of quantiles with respect to the parameters (``IStandardError.QuantileJacobian``).

        Parameters
        ----------
        probabilities : array_like
            Non-exceedance probabilities, one per parameter.

        Returns
        -------
        jacobian : numpy.ndarray
            Square matrix; row i is the gradient of the quantile at ``probabilities[i]``.
        determinant : float
            Its determinant.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IStandardError``.
        """
        jacobian, determinant = _interface(self.net, "IStandardError").QuantileJacobian(
            to_net_doubles(probabilities), 0.0)
        return to_numpy_matrix(jacobian), float(determinant)

    def log_abs_quantile_jacobian(self, probabilities: Any) -> float:
        """Return log abs(det J) of the quantile Jacobian (``StandardErrorExtensions.LogAbsQuantileJacobian``).

        Parameters
        ----------
        probabilities : array_like
            Non-exceedance probabilities, one per parameter.

        Returns
        -------
        float
            The natural logarithm of the absolute Jacobian determinant.

        Raises
        ------
        NotImplementedError
            If the family does not implement ``IStandardError``.
        """
        net = _interface(self.net, "IStandardError")
        return float(_nd().StandardErrorExtensions.LogAbsQuantileJacobian(net, to_net_doubles(probabilities)))


def _typed_member(net: Any, name: str, arity: int) -> Any:
    """Return the unbound ``name(double x arity)`` overload of ``net``'s concrete type, or None.

    Called as ``member(instance, *values)``. Binding the overload once on the type, rather than
    resolving the method group on every call, is what makes the per-row loops of
    ``parameter_set_*`` cheap; families without such an overload use the ``IList<double>`` one.
    """
    import System  # noqa: PLC0415 - the runtime is running (net exists)

    types = (System.Double,) * arity
    try:
        return getattr(type(net), name).Overloads[types if len(types) > 1 else types[0]]
    except (AttributeError, KeyError, TypeError):
        return None


class FrozenDistribution(UnivariateDistribution):
    """A read-only, hashable distribution, made by :meth:`UnivariateDistribution.frozen`.

    Each frozen object is an instance of ``FrozenDistribution`` and of the family's class (for
    example ``FrozenGeneralizedExtremeValue`` derives from both), and wraps its own ``Clone`` of the
    .NET object. ``set_parameters``, ``estimate`` and every attribute assignment (including the
    settable properties) raise; ``copy.copy`` and ``copy.deepcopy`` return the object itself;
    :meth:`clone` returns a mutable copy. It hashes by the value ``==`` compares. Changing the
    ``.net`` object through .NET directly bypasses the protection and must not be done.

    Parameters
    ----------
    net : Numerics.Distributions.UnivariateDistributionBase
        Unused; frozen distributions are made with :meth:`UnivariateDistribution.frozen`.
    """

    def __setattr__(self, name: str, value: Any) -> None:
        raise AttributeError(f"{type(self).__name__} is frozen; use clone() for a mutable copy.")

    def __delattr__(self, name: str) -> None:
        raise AttributeError(f"{type(self).__name__} is frozen; use clone() for a mutable copy.")

    def _ensure_mutable(self) -> None:
        raise TypeError(f"{type(self).__name__} is frozen; use clone() for a mutable copy.")

    def __hash__(self) -> int:
        key = self.__dict__.get("_hash")
        if key is None:
            key = self.__dict__["_hash"] = hash(self._value_key())
        return key

    def __copy__(self) -> FrozenDistribution:
        return self

    def __deepcopy__(self, memo: dict) -> FrozenDistribution:
        return self

    def frozen(self) -> FrozenDistribution:
        """Return this distribution, which is already frozen.

        Returns
        -------
        FrozenDistribution
            ``self``.
        """
        return self


_FROZEN_CLASSES: dict[type, type] = {}


def _frozen_class(cls: type) -> type:
    """Return the frozen class of a wrapper class, created once: ``Frozen<Name>(FrozenDistribution, cls)``."""
    if issubclass(cls, FrozenDistribution):
        return cls
    frozen = _FROZEN_CLASSES.get(cls)
    if frozen is None:
        if cls is UnivariateDistribution:
            frozen = FrozenDistribution
        else:
            frozen = types.new_class(f"Frozen{cls.__name__}", (FrozenDistribution, cls), exec_body=lambda ns: ns.update(
                {"__module__": cls.__module__, "__doc__": f"Frozen {cls.__name__}; see FrozenDistribution."}))
        _FROZEN_CLASSES[cls] = frozen
    return frozen


_CLASS_BY_NET_NAME: dict[str, type[UnivariateDistribution]] = {}
_FAMILIES: dict[str, type[UnivariateDistribution]] = {}

FAMILIES: Mapping[str, type[UnivariateDistribution]] = types.MappingProxyType(_FAMILIES)
"""The typed distribution classes, keyed by their Numerics type name (``UnivariateDistributionType``).

Read-only. ``FAMILIES[name].type_name == name`` for every entry, so a persisted type name finds its
class without an instance; families without a typed class are created by :func:`create_distribution`.
"""


def wrap(net: Any) -> UnivariateDistribution:
    """Wrap a .NET univariate distribution in the matching Python class.

    Parameters
    ----------
    net : Numerics.Distributions.UnivariateDistributionBase or IUnivariateDistribution
        The .NET object.

    Returns
    -------
    UnivariateDistribution
        The specific wrapper class when one exists, otherwise the generic wrapper.
    """
    concrete = unwrap(net)
    cls = _CLASS_BY_NET_NAME.get(str(concrete.GetType().Name), UnivariateDistribution)
    instance = cls.__new__(cls)
    UnivariateDistribution.__init__(instance, concrete)
    return instance


def _register(net_name: str, type_name: str | None = None):
    """Register a wrapper class for a .NET class name (for :func:`wrap`) and its type name (for :data:`FAMILIES`).

    ``type_name`` is the ``UnivariateDistributionType`` member, which is the .NET class name unless
    given (``BetaDistribution`` is ``Beta``).
    """

    def decorator(cls: type[UnivariateDistribution]) -> type[UnivariateDistribution]:
        name = type_name or net_name
        cls._type_name = name
        _CLASS_BY_NET_NAME[net_name] = cls
        _FAMILIES[name] = cls
        return cls

    return decorator


_LOCATION = "The location parameter ξ (``Xi``); settable."
_SCALE = "The scale parameter α > 0 (``Alpha``); settable."
_SHAPE = "The shape parameter κ, Hosking sign (``Kappa``); settable."


class _KappaMoments:
    """Members of the three-parameter Hosking families with a skew-to-shape solution (GEV, GPA, GLO)."""

    def solve_for_kappa(self, skew: float) -> float:
        """Return the shape κ whose product-moment skewness is ``skew`` (``SolveForKappa``).

        Parameters
        ----------
        skew : float
            The skewness coefficient.

        Returns
        -------
        float
            The shape κ, Hosking sign.
        """
        return float(self.net.SolveForKappa(float(skew)))

    def direct_method_of_moments(self, moments: Any) -> np.ndarray:
        """Return the parameters of the direct method of moments (``DirectMethodOfMoments``).

        Parameters
        ----------
        moments : array_like
            ``[mean, standard deviation, skewness]`` (further entries are ignored).

        Returns
        -------
        numpy.ndarray
            ``[location, scale, shape]``.
        """
        return to_numpy_doubles(self.net.DirectMethodOfMoments(to_net_doubles(moments)))


class _LogMoments:
    """Indirect (log-space) moment fits of the log-transformed families (LnNormal, LogNormal, LogPearsonTypeIII)."""

    def indirect_method_of_moments(self, sample: Any) -> np.ndarray:
        """Return the product moments of the log-transformed sample (``IndirectMethodOfMoments``).

        The indirect method of moments (U.S. Water Resources Council, 1967) fits the log-space
        parameters to these moments.

        Parameters
        ----------
        sample : array_like
            Positive observations.

        Returns
        -------
        numpy.ndarray
            ``[mean, standard deviation, skewness, kurtosis]`` of log_base(sample), in this
            distribution's base.
        """
        return to_numpy_doubles(self.net.IndirectMethodOfMoments(to_net_doubles(sample)))

    def indirect_method_of_linear_moments(self, sample: Any) -> np.ndarray:
        """Return the L-moments of the log-transformed sample (``IndirectMethodOfLinearMoments``).

        Parameters
        ----------
        sample : array_like
            Positive observations.

        Returns
        -------
        numpy.ndarray
            ``[λ1, λ2, τ3, τ4]`` of log_base(sample), in this distribution's base.
        """
        return to_numpy_doubles(self.net.IndirectMethodOfLinearMoments(to_net_doubles(sample)))


class _WilsonHilferty:
    """The Wilson-Hilferty quantile approximation of the gamma-type families."""

    def wilson_hilferty_inverse_cdf(self, probability: float) -> float:
        """Return the Wilson-Hilferty approximation of the quantile (``WilsonHilfertyInverseCDF``).

        Parameters
        ----------
        probability : float
            The non-exceedance probability.

        Returns
        -------
        float
            The approximate quantile; :meth:`UnivariateDistribution.inverse_cdf` gives the exact one.
        """
        return float(self.net.WilsonHilfertyInverseCDF(float(probability)))


class _PearsonForm(_WilsonHilferty):
    """The three-parameter gamma form (ξ, β, α) of the Pearson type III families."""

    support_bound = _net_property(
        "Xi", "The location ξ of the gamma form: the lower support bound for positive skew, the upper for "
              "negative (``Xi``; in log space for LogPearsonTypeIII).", settable=False)
    scale = _net_property("Beta", "The scale β of the gamma form, negative for negative skew (``Beta``).",
                          settable=False)
    shape = _net_property("Alpha", "The shape α of the gamma form, 4 / skew² (``Alpha``).", settable=False)

    def quantile_gradient_for_moments(self, probability: float) -> np.ndarray:
        """Return the quantile gradient in (mean, standard deviation, skew) (``QuantileGradientForMoments``).

        Parameters
        ----------
        probability : float
            The non-exceedance probability.

        Returns
        -------
        numpy.ndarray
            The three partial derivatives.
        """
        return to_numpy_doubles(self.net.QuantileGradientForMoments(float(probability)))


class _MonteCarloIntervals:
    """Monte Carlo quantile confidence intervals of the normal families."""

    def monte_carlo_confidence_intervals(self, sample_size: int, realizations: int, probabilities: Any,
                                         percentiles: Any) -> np.ndarray:
        """Return Monte Carlo confidence limits of quantiles (``MonteCarloConfidenceIntervals``).

        Parameters
        ----------
        sample_size : int
            The sample size of each realization.
        realizations : int
            The number of simulated samples.
        probabilities : array_like
            Non-exceedance probabilities of the quantiles.
        percentiles : array_like
            Confidence percentiles, for example ``[0.05, 0.95]``.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), len(percentiles))``.
        """
        return to_numpy_matrix(self.net.MonteCarloConfidenceIntervals(int(sample_size), int(realizations),
                                                                      to_net_doubles(probabilities),
                                                                      to_net_doubles(percentiles)))


for _mixin in (_KappaMoments, _LogMoments, _WilsonHilferty, _PearsonForm, _MonteCarloIntervals):
    translate_class(_mixin)


@_register("GeneralizedExtremeValue")
class GeneralizedExtremeValue(_KappaMoments, UnivariateDistribution):
    """Generalized extreme value distribution (Hosking parameterization).

    Parameters
    ----------
    location : float, optional
        Location ξ; default 100.
    scale : float, optional
        Scale α > 0; default 10.
    shape : float, optional
        Shape κ (Hosking sign: κ < 0 is heavy-tailed); default 0.
    """

    def __init__(self, location: float = 100.0, scale: float = 10.0, shape: float = 0.0) -> None:
        super().__init__(_nd().GeneralizedExtremeValue(float(location), float(scale), float(shape)))

    location = _net_property("Xi", _LOCATION)
    scale = _net_property("Alpha", _SCALE)
    shape = _net_property("Kappa", _SHAPE)

    @classmethod
    def from_poisson_generalized_pareto(cls, rate: float, location: float, scale: float,
                                        shape: float) -> GeneralizedExtremeValue:
        """Return the annual-maximum GEV of a Poisson count of generalized Pareto exceedances.

        Calls ``GeneralizedExtremeValue.FromPoissonGeneralizedPareto``: with N ~ Poisson(``rate``)
        exceedances a year of GPD(``location``, ``scale``, ``shape``), the maximum is GEV with
        ``scale * rate**(-shape)`` and the same shape (Hosking sign), the location formed through
        ``expm1`` so it stays accurate near ``shape = 0``.

        Parameters
        ----------
        rate : float
            The mean number of exceedances per year (> 0, finite).
        location : float
            The GPD location (the threshold).
        scale : float
            The GPD scale (> 0).
        shape : float
            The GPD shape κ (Hosking sign).

        Returns
        -------
        GeneralizedExtremeValue
            The annual-maximum distribution.

        Raises
        ------
        ValueError
            If ``rate`` or ``scale`` is not positive and finite.
        """
        member = _nd().GeneralizedExtremeValue.FromPoissonGeneralizedPareto
        instance = cls.__new__(cls)
        UnivariateDistribution.__init__(instance, member(float(rate), float(location), float(scale), float(shape)))
        return instance

    def expected_information_matrix(self, sample_size: int) -> np.ndarray:
        """Return the expected (Fisher) information matrix of a sample (``ExpectedInformationMatrix``).

        Parameters
        ----------
        sample_size : int
            The sample size n.

        Returns
        -------
        numpy.ndarray
            The 3 x 3 information matrix in parameter order.
        """
        return to_numpy_matrix(self.net.ExpectedInformationMatrix(int(sample_size)).ToArray())


@_register("GeneralizedPareto")
class GeneralizedPareto(_KappaMoments, UnivariateDistribution):
    """Generalized Pareto distribution (Hosking parameterization).

    Parameters
    ----------
    location : float, optional
        Location ξ; default 100.
    scale : float, optional
        Scale α > 0; default 10.
    shape : float, optional
        Shape κ (Hosking sign); default 0.
    """

    def __init__(self, location: float = 100.0, scale: float = 10.0, shape: float = 0.0) -> None:
        super().__init__(_nd().GeneralizedPareto(float(location), float(scale), float(shape)))

    location = _net_property("Xi", _LOCATION)
    scale = _net_property("Alpha", _SCALE)
    shape = _net_property("Kappa", _SHAPE)
    peaks_per_block = _net_property("Lambda", "The average number of peaks per block (``Lambda``); settable.")

    @classmethod
    def from_poisson_generalized_extreme_value(cls, rate: float, location: float, scale: float,
                                               shape: float) -> GeneralizedPareto:
        """Return the exceedance GPD whose Poisson annual count gives a stated GEV annual maximum.

        Calls ``GeneralizedPareto.FromPoissonGeneralizedExtremeValue``, the inverse of
        :meth:`GeneralizedExtremeValue.from_poisson_generalized_pareto`: for an annual-maximum
        GEV(``location``, ``scale``, ``shape``) and N ~ Poisson(``rate``) exceedances a year, the
        exceedances are GPD with scale ``scale * rate**shape`` and the same shape (Hosking sign), the
        location ``location - scale * expm1(shape * log(rate)) / shape`` formed through ``expm1`` so it
        stays accurate near ``shape = 0`` (``location - scale * log(rate)`` at ``shape = 0``).

        Parameters
        ----------
        rate : float
            The mean number of exceedances per year (> 0, finite).
        location : float
            The GEV location.
        scale : float
            The GEV scale (> 0).
        shape : float
            The GEV shape κ (Hosking sign).

        Returns
        -------
        GeneralizedPareto
            The exceedance distribution; its location is the threshold.

        Raises
        ------
        ValueError
            If ``rate`` or ``scale`` is not positive and finite, ``location`` or ``shape`` is not finite,
            or the mapped parameters overflow.
        """
        member = _nd().GeneralizedPareto.FromPoissonGeneralizedExtremeValue
        instance = cls.__new__(cls)
        UnivariateDistribution.__init__(instance, member(float(rate), float(location), float(scale), float(shape)))
        return instance

    @classmethod
    def fit_known_location(cls, sample: Any, location: float) -> GeneralizedPareto:
        """Fit the scale and shape by maximum likelihood at a known location (``MLE(sample, location)``).

        Parameters
        ----------
        sample : array_like
            Observations, each at or above ``location`` (for example the peaks over a threshold).
        location : float
            The known location ξ (for example the threshold).

        Returns
        -------
        GeneralizedPareto
            The fitted distribution, with location ``location``.
        """
        parameters = to_numpy_doubles(_nd().GeneralizedPareto().MLE(to_net_doubles(sample), float(location)))
        return cls(float(parameters[0]), float(parameters[1]), float(parameters[2]))

    def mle(self, sample: Any, location: float | None = None) -> np.ndarray:
        """Return the maximum likelihood parameters of a sample, without changing this distribution.

        Parameters
        ----------
        sample : array_like
            Observations.
        location : float, optional
            A known location; when given only the scale and shape are fitted (``MLE(sample, location)``).

        Returns
        -------
        numpy.ndarray
            ``[location, scale, shape]``.
        """
        if location is None:
            return super().mle(sample)
        return to_numpy_doubles(self.net.MLE(to_net_doubles(sample), float(location)))

    def modified_method_of_moments(self, sample: Any) -> np.ndarray:
        """Return the parameters of the modified method of moments (``ModifiedMethodOfMoments``).

        Parameters
        ----------
        sample : array_like
            Observations.

        Returns
        -------
        numpy.ndarray
            ``[location, scale, shape]``.
        """
        return to_numpy_doubles(self.net.ModifiedMethodOfMoments(to_net_doubles(sample)))


@_register("Gumbel")
class Gumbel(UnivariateDistribution):
    """Gumbel (extreme value type I) distribution.

    Parameters
    ----------
    location : float, optional
        Location ξ; default 100.
    scale : float, optional
        Scale α > 0; default 10.
    """

    def __init__(self, location: float = 100.0, scale: float = 10.0) -> None:
        super().__init__(_nd().Gumbel(float(location), float(scale)))

    location = _net_property("Xi", _LOCATION)
    scale = _net_property("Alpha", _SCALE)

    def set_parameters_from_mle(self, sample: Any) -> None:
        """Fit the parameters in place by maximum likelihood (``SetParametersFromMLE``).

        Parameters
        ----------
        sample : array_like
            Observations.

        Raises
        ------
        TypeError
            If the distribution is frozen.
        """
        self._ensure_mutable()
        self.net.SetParametersFromMLE(to_net_doubles(sample))


@_register("Normal")
class Normal(_MonteCarloIntervals, UnivariateDistribution):
    """Normal distribution.

    Parameters
    ----------
    mean : float, optional
        Mean μ; default 0.
    standard_deviation : float, optional
        Standard deviation σ > 0; default 1.
    """

    def __init__(self, mean: float = 0.0, standard_deviation: float = 1.0) -> None:
        super().__init__(_nd().Normal(float(mean), float(standard_deviation)))

    mean = _net_property("Mu", "The mean μ (``Mu``); settable.")
    standard_deviation = _net_property("Sigma", "The standard deviation σ > 0 (``Sigma``); settable.")

    def normal_confidence_intervals(self, sample_size: int, probabilities: Any, percentiles: Any) -> np.ndarray:
        """Return normal-approximation confidence limits of quantiles (``NormalConfidenceIntervals``).

        Parameters
        ----------
        sample_size : int
            The sample size.
        probabilities : array_like
            Non-exceedance probabilities of the quantiles.
        percentiles : array_like
            Confidence percentiles, for example ``[0.05, 0.95]``.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), len(percentiles))``.
        """
        return to_numpy_matrix(self.net.NormalConfidenceIntervals(int(sample_size), to_net_doubles(probabilities),
                                                                  to_net_doubles(percentiles)))

    def noncentral_t_confidence_intervals(self, sample_size: int, probabilities: Any,
                                          percentiles: Any) -> np.ndarray:
        """Return exact noncentral-t confidence limits of quantiles (``NoncentralTConfidenceIntervals``).

        Parameters
        ----------
        sample_size : int
            The sample size.
        probabilities : array_like
            Non-exceedance probabilities of the quantiles.
        percentiles : array_like
            Confidence percentiles, for example ``[0.05, 0.95]``.

        Returns
        -------
        numpy.ndarray
            Array of shape ``(len(probabilities), len(percentiles))``.
        """
        return to_numpy_matrix(self.net.NoncentralTConfidenceIntervals(
            int(sample_size), to_net_doubles(probabilities), to_net_doubles(percentiles)))

    def expected_probability(self, sample_size: int, probability: float) -> float:
        """Return the expected non-exceedance probability of the estimated quantile (``ExpectedProbability``).

        Parameters
        ----------
        sample_size : int
            The sample size.
        probability : float
            The nominal non-exceedance probability.

        Returns
        -------
        float
            The expected probability, which accounts for parameter uncertainty.
        """
        return float(self.net.ExpectedProbability(int(sample_size), float(probability)))


@_register("LogNormal")
class LogNormal(_LogMoments, _MonteCarloIntervals, UnivariateDistribution):
    """Log-normal distribution with log-space parameters (``Numerics.Distributions.LogNormal``).

    Parameters
    ----------
    mean_of_log : float, optional
        Mean of log_base(X); default 3.
    standard_deviation_of_log : float, optional
        Standard deviation of log_base(X); default 0.5.
    base : float, optional
        The logarithm base, finite and greater than one; default 10.
    """

    def __init__(self, mean_of_log: float = 3.0, standard_deviation_of_log: float = 0.5, base: float = 10.0) -> None:
        net = _nd().LogNormal(float(mean_of_log), float(standard_deviation_of_log))
        net.Base = float(base)
        super().__init__(net)

    mean_of_log = _net_property("Mu", "The mean of log_base(X) (``Mu``); settable.")
    standard_deviation_of_log = _net_property("Sigma", "The standard deviation of log_base(X) (``Sigma``); settable.")
    base = _net_property("Base", "The logarithm base, finite and greater than one (``Base``); settable.")

    def direct_method_of_moments(self, mean: float, standard_deviation: float) -> tuple[float, float]:
        """Return (mean_of_log, standard_deviation_of_log) for real-space moments (``DirectMethodOfMoments``).

        Parameters
        ----------
        mean : float
            The mean of X.
        standard_deviation : float
            The standard deviation of X.

        Returns
        -------
        tuple of float
            The log-space parameters in this distribution's base.
        """
        values = self.net.DirectMethodOfMoments(float(mean), float(standard_deviation))
        return float(values[0]), float(values[1])


@_register("LnNormal")
class LnNormal(_LogMoments, UnivariateDistribution):
    """Natural-log normal distribution with real-space moments (``LnNormal``).

    Parameters
    ----------
    mean : float, optional
        Real-space mean of X; default 10.
    standard_deviation : float, optional
        Real-space standard deviation of X; default 10.
    """

    def __init__(self, mean: float = 10.0, standard_deviation: float = 10.0) -> None:
        super().__init__(_nd().LnNormal(float(mean), float(standard_deviation)))

    mean_of_log = _net_property("Mu", "The mean μ of ln X (``Mu``); settable.")
    standard_deviation_of_log = _net_property("Sigma", "The standard deviation σ of ln X (``Sigma``); settable.")

    @staticmethod
    def direct_method_of_moments(mean: float, standard_deviation: float) -> tuple[float, float]:
        """Return (mu, sigma) of ln X for real-space moments (``LnNormal.DirectMethodOfMoments``).

        Parameters
        ----------
        mean : float
            The mean of X.
        standard_deviation : float
            The standard deviation of X.

        Returns
        -------
        tuple of float
            ``(mu, sigma)`` as Numerics computes them, before the ``Sigma`` setter's clamp.
        """
        values = _nd().LnNormal.DirectMethodOfMoments(float(mean), float(standard_deviation))
        return float(values[0]), float(values[1])

    @classmethod
    def from_log_parameters(cls, mu: float, sigma: float) -> LnNormal:
        """Create the distribution from mu and sigma of ln X (the ``Mu`` and ``Sigma`` setters).

        Parameters
        ----------
        mu : float
            The mean of ln X.
        sigma : float
            The standard deviation of ln X.

        Returns
        -------
        LnNormal
            The distribution; Numerics derives the real-space moments from mu and sigma.
        """
        net = _nd().LnNormal()
        net.Mu = float(mu)
        net.Sigma = float(sigma)
        instance = cls.__new__(cls)
        UnivariateDistribution.__init__(instance, net)
        return instance

    def indirect_method_of_moments(self, sample: Any) -> np.ndarray:
        """Return the product moments of ln(sample) (``LnNormal.IndirectMethodOfMoments``).

        Parameters
        ----------
        sample : array_like
            Positive observations.

        Returns
        -------
        numpy.ndarray
            ``[mean, standard deviation, skewness, kurtosis]`` of the natural logarithms.
        """
        return to_numpy_doubles(_nd().LnNormal.IndirectMethodOfMoments(to_net_doubles(sample)))


@_register("PearsonTypeIII")
class PearsonTypeIII(_PearsonForm, UnivariateDistribution):
    """Pearson type III distribution.

    Parameters
    ----------
    mean : float, optional
        Mean μ; default 100.
    standard_deviation : float, optional
        Standard deviation σ > 0; default 10.
    skew : float, optional
        Skewness γ; default 0.
    """

    def __init__(self, mean: float = 100.0, standard_deviation: float = 10.0, skew: float = 0.0) -> None:
        super().__init__(_nd().PearsonTypeIII(float(mean), float(standard_deviation), float(skew)))

    mean = _net_property("Mu", "The mean μ (``Mu``); settable.")
    standard_deviation = _net_property("Sigma", "The standard deviation σ > 0 (``Sigma``); settable.")
    skew = _net_property("Gamma", "The skewness γ (``Gamma``); settable.")


@_register("LogPearsonTypeIII")
class LogPearsonTypeIII(_LogMoments, _PearsonForm, UnivariateDistribution):
    """Log-Pearson type III distribution: Pearson type III in log_base(X).

    Parameters
    ----------
    mean_of_log : float, optional
        Mean of log_base(X); default 3.
    standard_deviation_of_log : float, optional
        Standard deviation of log_base(X); default 0.5.
    skew_of_log : float, optional
        Skewness of log_base(X); default 0.
    base : float, optional
        The logarithm base, finite and greater than one; default 10.
    """

    def __init__(self, mean_of_log: float = 3.0, standard_deviation_of_log: float = 0.5, skew_of_log: float = 0.0,
                 base: float = 10.0) -> None:
        net = _nd().LogPearsonTypeIII(float(mean_of_log), float(standard_deviation_of_log), float(skew_of_log))
        net.Base = float(base)
        super().__init__(net)

    mean_of_log = _net_property("Mu", "The mean of log_base(X) (``Mu``); settable.")
    standard_deviation_of_log = _net_property("Sigma", "The standard deviation of log_base(X) (``Sigma``); settable.")
    skew_of_log = _net_property("Gamma", "The skewness of log_base(X) (``Gamma``); settable.")
    base = _net_property("Base", "The logarithm base, finite and greater than one (``Base``); settable.")

    def log_space(self) -> PearsonTypeIII:
        """Return the Pearson type III distribution of log_base(X).

        Returns
        -------
        PearsonTypeIII
            ``PearsonTypeIII(mean_of_log, standard_deviation_of_log, skew_of_log)``.
        """
        return PearsonTypeIII(self.mean_of_log, self.standard_deviation_of_log, self.skew_of_log)


@_register("KappaFour")
class KappaFour(UnivariateDistribution):
    """Four-parameter kappa distribution (Hosking).

    Parameters
    ----------
    location : float, optional
        Location ξ; default 100.
    scale : float, optional
        Scale α > 0; default 10.
    shape : float, optional
        Shape κ; default 0.
    shape2 : float, optional
        Second shape h; default 0.
    """

    def __init__(self, location: float = 100.0, scale: float = 10.0, shape: float = 0.0, shape2: float = 0.0) -> None:
        super().__init__(_nd().KappaFour(float(location), float(scale), float(shape), float(shape2)))

    location = _net_property("Xi", _LOCATION)
    scale = _net_property("Alpha", _SCALE)
    shape = _net_property("Kappa", _SHAPE)
    shape2 = _net_property("Hondo", "The second shape parameter h (``Hondo``); settable.")


@_register("Uniform")
class Uniform(UnivariateDistribution):
    """Continuous uniform distribution on [minimum, maximum] (``Uniform``).

    Parameters
    ----------
    minimum : float, optional
        Lower end; default 0.
    maximum : float, optional
        Upper end; default 1.
    """

    def __init__(self, minimum: float = 0.0, maximum: float = 1.0) -> None:
        super().__init__(_nd().Uniform(float(minimum), float(maximum)))

    minimum = _net_property("Min", "The lower end of the support (``Min``); settable.")
    maximum = _net_property("Max", "The upper end of the support (``Max``); settable.")


@_register("BetaDistribution", "Beta")
class Beta(UnivariateDistribution):
    """Beta distribution on [0, 1] (``BetaDistribution``).

    Parameters
    ----------
    alpha : float, optional
        First shape α > 0; default 2.
    beta : float, optional
        Second shape β > 0; default 2.
    """

    def __init__(self, alpha: float = 2.0, beta: float = 2.0) -> None:
        super().__init__(_nd().BetaDistribution(float(alpha), float(beta)))

    alpha = _net_property("Alpha", "The first shape α > 0 (``Alpha``); settable.")
    beta = _net_property("Beta", "The second shape β > 0 (``Beta``); settable.")


@_register("Poisson")
class Poisson(UnivariateDistribution):
    """Poisson distribution (discrete; ``pdf`` is the probability mass).

    Parameters
    ----------
    rate : float, optional
        Rate λ > 0; default 1.
    """

    def __init__(self, rate: float = 1.0) -> None:
        super().__init__(_nd().Poisson(float(rate)))

    rate = _net_property("Lambda", "The rate λ > 0 (``Lambda``); settable.")


@_register("Bernoulli")
class Bernoulli(UnivariateDistribution):
    """Bernoulli distribution on {0, 1} (discrete; ``pdf`` is the probability mass).

    Parameters
    ----------
    probability : float, optional
        Probability of a one, in [0, 1]; default 0.5.
    """

    def __init__(self, probability: float = 0.5) -> None:
        super().__init__(_nd().Bernoulli(float(probability)))

    probability = _net_property("Probability", "The probability p of a one (``Probability``); settable.")
    complement = _net_property("Complement", "The probability 1 - p of a zero (``Complement``).", settable=False)


@_register("Binomial")
class Binomial(UnivariateDistribution):
    """Binomial distribution: successes in ``number_of_trials`` trials (discrete).

    Parameters
    ----------
    probability : float, optional
        Probability of success per trial, in [0, 1]; default 0.5.
    number_of_trials : int, optional
        The number of trials n >= 0; default 10.
    """

    def __init__(self, probability: float = 0.5, number_of_trials: int = 10) -> None:
        super().__init__(_nd().Binomial(float(probability), int(number_of_trials)))

    probability = _net_property("ProbabilityOfSuccess", "The probability of success p (``ProbabilityOfSuccess``); "
                                                        "settable.")
    number_of_trials = _net_property("NumberOfTrials", "The number of trials n (``NumberOfTrials``); settable.",
                                     kind=int)
    complement = _net_property("Complement", "The probability of failure 1 - p (``Complement``).", settable=False)


@_register("Cauchy")
class Cauchy(UnivariateDistribution):
    """Cauchy distribution.

    Parameters
    ----------
    location : float, optional
        Location x0; default 0.
    scale : float, optional
        Scale γ > 0; default 1.
    """

    def __init__(self, location: float = 0.0, scale: float = 1.0) -> None:
        super().__init__(_nd().Cauchy(float(location), float(scale)))

    location = _net_property("X0", "The location x0 (``X0``); settable.")
    scale = _net_property("Gamma", "The scale γ > 0 (``Gamma``); settable.")


@_register("ChiSquared")
class ChiSquared(UnivariateDistribution):
    """Chi-squared distribution.

    Parameters
    ----------
    degrees_of_freedom : int, optional
        Degrees of freedom ν >= 1; default 10.
    """

    def __init__(self, degrees_of_freedom: int = 10) -> None:
        super().__init__(_nd().ChiSquared(int(degrees_of_freedom)))

    degrees_of_freedom = _net_property("DegreesOfFreedom", "The degrees of freedom ν (``DegreesOfFreedom``); "
                                                           "settable.", kind=int)


@_register("Deterministic")
class Deterministic(UnivariateDistribution):
    """A point mass at one value (``Deterministic``).

    Parameters
    ----------
    value : float, optional
        The value; default 0.5.
    """

    def __init__(self, value: float = 0.5) -> None:
        super().__init__(_nd().Deterministic(float(value)))

    value = _net_property("Value", "The value of the point mass (``Value``); settable.")


@_register("Exponential")
class Exponential(UnivariateDistribution):
    """Shifted exponential distribution: ``F(x) = 1 - exp(-(x - ξ) / α)`` for x >= ξ.

    Parameters
    ----------
    location : float, optional
        Location ξ; default 100.
    scale : float, optional
        Scale α > 0 (the mean excess); default 10.
    """

    def __init__(self, location: float = 100.0, scale: float = 10.0) -> None:
        super().__init__(_nd().Exponential(float(location), float(scale)))

    location = _net_property("Xi", _LOCATION)
    scale = _net_property("Alpha", _SCALE)


@_register("GammaDistribution")
class GammaDistribution(_WilsonHilferty, UnivariateDistribution):
    """Gamma distribution with scale θ and shape κ (``GammaDistribution``).

    Parameters
    ----------
    scale : float, optional
        Scale θ > 0; default 10.
    shape : float, optional
        Shape κ > 0; default 2.
    """

    def __init__(self, scale: float = 10.0, shape: float = 2.0) -> None:
        super().__init__(_nd().GammaDistribution(float(scale), float(shape)))

    scale = _net_property("Theta", "The scale θ > 0 (``Theta``); settable.")
    shape = _net_property("Kappa", "The shape κ > 0 (``Kappa``); settable.")
    rate = _net_property("Rate", "The rate 1 / θ (``Rate``).", settable=False)

    def estimate_mle_newton_raphson(self, sample: Any) -> None:
        """Fit the parameters in place by maximum likelihood with Newton-Raphson iterations (``MLE_NR``).

        Parameters
        ----------
        sample : array_like
            Positive observations.

        Raises
        ------
        TypeError
            If the distribution is frozen.
        """
        self._ensure_mutable()
        self.net.MLE_NR(to_net_doubles(sample))

    def estimate_mle_bobee(self, sample: Any) -> None:
        """Fit the parameters in place by the maximum likelihood approximation of Bobee (``MLE_Bobee``).

        Parameters
        ----------
        sample : array_like
            Positive observations.

        Raises
        ------
        TypeError
            If the distribution is frozen.
        """
        self._ensure_mutable()
        self.net.MLE_Bobee(to_net_doubles(sample))

    @staticmethod
    def frequency_factor_kp(skewness: float, probability: float) -> float:
        """Return the approximate Pearson type III frequency factor K(γ, p) (``GammaDistribution.FrequencyFactorKp``).

        The Cornish-Fisher expansion for abs(γ) <= 2 and the modified Wilson-Hilferty transformation
        (Kirby, 1972) beyond; it approximates, and is not, the exact standardized quantile.

        Parameters
        ----------
        skewness : float
            The skewness γ.
        probability : float
            The non-exceedance probability.

        Returns
        -------
        float
            The standardized quantile, so the quantile is ``mean + K * standard_deviation``.
        """
        return float(_nd().GammaDistribution.FrequencyFactorKp(float(skewness), float(probability)))

    @staticmethod
    def partial_kp(skewness: float, probability: float) -> float:
        """Return the derivative of the frequency factor with respect to the skewness (``PartialKp``).

        Parameters
        ----------
        skewness : float
            The skewness γ.
        probability : float
            The non-exceedance probability.

        Returns
        -------
        float
            The derivative dK/dγ at ``(skewness, probability)``.
        """
        return float(_nd().GammaDistribution.PartialKp(float(skewness), float(probability)))


@_register("GeneralizedBeta")
class GeneralizedBeta(UnivariateDistribution):
    """Beta distribution rescaled to [minimum, maximum] (``GeneralizedBeta``).

    Parameters
    ----------
    alpha : float, optional
        First shape α > 0; default 2.
    beta : float, optional
        Second shape β > 0; default 2.
    minimum : float, optional
        Lower end; default 0.
    maximum : float, optional
        Upper end; default 1.
    """

    def __init__(self, alpha: float = 2.0, beta: float = 2.0, minimum: float = 0.0, maximum: float = 1.0) -> None:
        super().__init__(_nd().GeneralizedBeta(float(alpha), float(beta), float(minimum), float(maximum)))

    alpha = _net_property("Alpha", "The first shape α > 0 (``Alpha``); settable.")
    beta = _net_property("Beta", "The second shape β > 0 (``Beta``); settable.")
    minimum = _net_property("Min", "The lower end of the support (``Min``); settable.")
    maximum = _net_property("Max", "The upper end of the support (``Max``); settable.")

    @classmethod
    def pert(cls, minimum: float, mode: float, maximum: float, scale: float | None = None) -> GeneralizedBeta:
        """Return the PERT beta distribution through a mode (``GeneralizedBeta.PERT``).

        Parameters
        ----------
        minimum : float
            Lower end.
        mode : float
            The most likely value.
        maximum : float
            Upper end.
        scale : float, optional
            The PERT weight λ of the mode; None, the default, is Numerics' ``PERT(min, mode, max)``
            (λ = 4).

        Returns
        -------
        GeneralizedBeta
            ``α = 1 + λ (mode - min) / (max - min)``, ``β = 1 + λ (max - mode) / (max - min)``.
        """
        pert = _nd().GeneralizedBeta.PERT
        net = (pert(float(minimum), float(mode), float(maximum)) if scale is None
               else pert(float(minimum), float(mode), float(maximum), float(scale)))
        instance = cls.__new__(cls)
        UnivariateDistribution.__init__(instance, net)
        return instance

    def set_parameters_from_moments(self, mean: float, standard_deviation: float, minimum: float,
                                    maximum: float) -> None:
        """Set the shapes in place from a mean and standard deviation on [minimum, maximum].

        Parameters
        ----------
        mean : float
            The mean.
        standard_deviation : float
            The standard deviation.
        minimum : float
            Lower end.
        maximum : float
            Upper end.

        Raises
        ------
        TypeError
            If the distribution is frozen.

        Notes
        -----
        Calls ``SetParametersFromMoments``.
        """
        self._ensure_mutable()
        self.net.SetParametersFromMoments(float(mean), float(standard_deviation), float(minimum), float(maximum))


@_register("GeneralizedLogistic")
class GeneralizedLogistic(_KappaMoments, UnivariateDistribution):
    """Generalized logistic distribution (Hosking parameterization).

    ``F(x) = 1 / (1 + exp(-y))`` with ``y = -ln(1 - κ (x - ξ) / α) / κ`` (``y = (x - ξ) / α`` at κ = 0).

    Parameters
    ----------
    location : float, optional
        Location ξ; default 100.
    scale : float, optional
        Scale α > 0; default 10.
    shape : float, optional
        Shape κ (Hosking sign); default 0.
    """

    def __init__(self, location: float = 100.0, scale: float = 10.0, shape: float = 0.0) -> None:
        super().__init__(_nd().GeneralizedLogistic(float(location), float(scale), float(shape)))

    location = _net_property("Xi", _LOCATION)
    scale = _net_property("Alpha", _SCALE)
    shape = _net_property("Kappa", _SHAPE)

    @_FamilyMethod
    def moments_from_parameters(self, parameters: Any) -> np.ndarray:
        """Return the product moments of a parameter vector (``MomentsFromParameters``).

        Callable on the class without an instance or on an instance.

        Parameters
        ----------
        parameters : array_like
            ``[location, scale, shape]``.

        Returns
        -------
        numpy.ndarray
            ``[mean, standard deviation, skewness, kurtosis]``.
        """
        return to_numpy_doubles(_net_of(self).MomentsFromParameters(to_net_doubles(parameters)))


@_register("GeneralizedNormal")
class GeneralizedNormal(UnivariateDistribution):
    """Generalized normal distribution (Hosking parameterization).

    ``F(x) = Φ(y)`` with ``y = -ln(1 - κ (x - ξ) / α) / κ`` (``y = (x - ξ) / α`` at κ = 0).

    Parameters
    ----------
    location : float, optional
        Location ξ; default 100.
    scale : float, optional
        Scale α > 0; default 10.
    shape : float, optional
        Shape κ (Hosking sign); default 0.
    """

    def __init__(self, location: float = 100.0, scale: float = 10.0, shape: float = 0.0) -> None:
        super().__init__(_nd().GeneralizedNormal(float(location), float(scale), float(shape)))

    location = _net_property("Xi", _LOCATION)
    scale = _net_property("Alpha", _SCALE)
    shape = _net_property("Kappa", _SHAPE)


@_register("Geometric")
class Geometric(UnivariateDistribution):
    """Geometric distribution: failures before the first success (discrete, support 0, 1, 2, ...).

    Parameters
    ----------
    probability : float, optional
        Probability of success per trial, in (0, 1]; default 0.5.
    """

    def __init__(self, probability: float = 0.5) -> None:
        super().__init__(_nd().Geometric(float(probability)))

    probability = _net_property("ProbabilityOfSuccess", "The probability of success p (``ProbabilityOfSuccess``); "
                                                        "settable.")


@_register("InverseChiSquared")
class InverseChiSquared(UnivariateDistribution):
    """Scaled inverse chi-squared distribution (``InverseChiSquared``).

    The inverse gamma distribution with shape ν/2 and scale ν τ²/2; the mean is ν τ²/(ν - 2) for ν > 2.

    Parameters
    ----------
    degrees_of_freedom : int, optional
        Degrees of freedom ν >= 1; default 10.
    scale : float, optional
        The scale τ² > 0 (Numerics ``Sigma``); default 1.
    """

    def __init__(self, degrees_of_freedom: int = 10, scale: float = 1.0) -> None:
        super().__init__(_nd().InverseChiSquared(int(degrees_of_freedom), float(scale)))

    degrees_of_freedom = _net_property("DegreesOfFreedom", "The degrees of freedom ν (``DegreesOfFreedom``); "
                                                           "settable.", kind=int)
    scale = _net_property("Sigma", "The scale τ² > 0 (``Sigma``); settable.")
    precision = _net_property("Precision", "The value 1 / τ⁴, the reciprocal of the squared scale (``Precision``).",
                              settable=False)


@_register("InverseGamma")
class InverseGamma(UnivariateDistribution):
    """Inverse gamma distribution with scale β and shape α (``InverseGamma``).

    Parameters
    ----------
    scale : float, optional
        Scale β > 0; default 0.5.
    shape : float, optional
        Shape α > 0; default 2.
    """

    def __init__(self, scale: float = 0.5, shape: float = 2.0) -> None:
        super().__init__(_nd().InverseGamma(float(scale), float(shape)))

    scale = _net_property("Beta", "The scale β > 0 (``Beta``); settable.")
    shape = _net_property("Alpha", "The shape α > 0 (``Alpha``); settable.")


@_register("Logistic")
class Logistic(UnivariateDistribution):
    """Logistic distribution.

    Parameters
    ----------
    location : float, optional
        Location ξ; default 0.
    scale : float, optional
        Scale α > 0; default 0.1.
    """

    def __init__(self, location: float = 0.0, scale: float = 0.1) -> None:
        super().__init__(_nd().Logistic(float(location), float(scale)))

    location = _net_property("Xi", _LOCATION)
    scale = _net_property("Alpha", _SCALE)


@_register("NoncentralT")
class NoncentralT(UnivariateDistribution):
    """Noncentral Student's t distribution.

    Parameters
    ----------
    degrees_of_freedom : float, optional
        Degrees of freedom ν > 0; default 10.
    noncentrality : float, optional
        Noncentrality μ; default 0.
    """

    def __init__(self, degrees_of_freedom: float = 10.0, noncentrality: float = 0.0) -> None:
        super().__init__(_nd().NoncentralT(float(degrees_of_freedom), float(noncentrality)))

    degrees_of_freedom = _net_property("DegreesOfFreedom", "The degrees of freedom ν > 0 (``DegreesOfFreedom``); "
                                                           "settable.")
    noncentrality = _net_property("Noncentrality", "The noncentrality μ (``Noncentrality``); settable.")


@_register("Pareto")
class Pareto(UnivariateDistribution):
    """Pareto (type I) distribution: ``F(x) = 1 - (xm / x)^α`` for x >= xm.

    Parameters
    ----------
    scale : float, optional
        Scale xm > 0, the lower end of the support; default 1.
    shape : float, optional
        Shape α > 0; default 10.
    """

    def __init__(self, scale: float = 1.0, shape: float = 10.0) -> None:
        super().__init__(_nd().Pareto(float(scale), float(shape)))

    scale = _net_property("Xm", "The scale xm > 0 (``Xm``); settable.")
    shape = _net_property("Alpha", "The shape α > 0 (``Alpha``); settable.")


@_register("Pert")
class Pert(UnivariateDistribution):
    """PERT distribution: a beta on [minimum, maximum] weighted towards the most likely value.

    Parameters
    ----------
    minimum : float, optional
        Lower end a; default 0.
    most_likely : float, optional
        The mode c; default 0.5.
    maximum : float, optional
        Upper end b; default 1.
    """

    def __init__(self, minimum: float = 0.0, most_likely: float = 0.5, maximum: float = 1.0) -> None:
        super().__init__(_nd().Pert(float(minimum), float(most_likely), float(maximum)))

    minimum = _net_property("Min", "The lower end a (``Min``); settable.")
    most_likely = _net_property("MostLikely", "The most likely value c (``MostLikely``); settable.")
    maximum = _net_property("Max", "The upper end b (``Max``); settable.")

    @property
    def beta(self) -> GeneralizedBeta:
        """The underlying beta distribution on [a, b], wrapped (``Beta``)."""
        return wrap(self.net.Beta)


class _PertPercentiles:
    """The percentile parameters shared by the percentile PERT families."""

    percentile_5th = _net_property("Percentile5th", "The 5th percentile (``Percentile5th``); settable.")
    percentile_50th = _net_property("Percentile50th", "The 50th percentile (``Percentile50th``); settable.")
    percentile_95th = _net_property("Percentile95th", "The 95th percentile (``Percentile95th``); settable.")

    def solve_parameters(self) -> None:
        """Solve the PERT parameters that reproduce the three percentiles (``SolveParameters``).

        Raises
        ------
        TypeError
            If the distribution is frozen.
        """
        self._ensure_mutable()
        self.net.SolveParameters()

    def to_pert(self) -> Pert:
        """Return the solved PERT distribution (``ToPert``).

        Returns
        -------
        Pert
            The PERT distribution whose percentiles match.
        """
        return wrap(self.net.ToPert())


translate_class(_PertPercentiles)


@_register("PertPercentile")
class PertPercentile(_PertPercentiles, UnivariateDistribution):
    """PERT distribution given by its 5th, 50th and 95th percentiles (``PertPercentile``).

    Parameters
    ----------
    percentile_5th : float, optional
        The 5th percentile; default 0.05.
    percentile_50th : float, optional
        The median; default 0.5.
    percentile_95th : float, optional
        The 95th percentile; default 0.95.
    """

    def __init__(self, percentile_5th: float = 0.05, percentile_50th: float = 0.5,
                 percentile_95th: float = 0.95) -> None:
        super().__init__(_nd().PertPercentile(float(percentile_5th), float(percentile_50th), float(percentile_95th)))

    min_allowable_value = _net_property("MinAllowableValue", "The smallest value the solved support may take "
                                                             "(``MinAllowableValue``); settable.")
    max_allowable_value = _net_property("MaxAllowableValue", "The largest value the solved support may take "
                                                             "(``MaxAllowableValue``); settable.")


@_register("PertPercentileZ")
class PertPercentileZ(_PertPercentiles, UnivariateDistribution):
    """PERT distribution given by three percentiles, solved in a transformed space (``PertPercentileZ``).

    Parameters
    ----------
    percentile_5th : float, optional
        The 5th percentile, in (0, 1); default 0.05.
    percentile_50th : float, optional
        The median, in (0, 1); default 0.5.
    percentile_95th : float, optional
        The 95th percentile, in (0, 1); default 0.95.
    """

    def __init__(self, percentile_5th: float = 0.05, percentile_50th: float = 0.5,
                 percentile_95th: float = 0.95) -> None:
        super().__init__(_nd().PertPercentileZ(float(percentile_5th), float(percentile_50th),
                                               float(percentile_95th)))


@_register("Rayleigh")
class Rayleigh(UnivariateDistribution):
    """Rayleigh distribution.

    Parameters
    ----------
    scale : float, optional
        Scale σ > 0; default 10.
    """

    def __init__(self, scale: float = 10.0) -> None:
        super().__init__(_nd().Rayleigh(float(scale)))

    scale = _net_property("Sigma", "The scale σ > 0 (``Sigma``); settable.")


@_register("StudentT")
class StudentT(UnivariateDistribution):
    """Location-scale Student's t distribution.

    Parameters
    ----------
    location : float, optional
        Location μ; default 0.
    scale : float, optional
        Scale σ > 0; default 1.
    degrees_of_freedom : float, optional
        Degrees of freedom ν > 0; default 10.
    """

    def __init__(self, location: float = 0.0, scale: float = 1.0, degrees_of_freedom: float = 10.0) -> None:
        super().__init__(_nd().StudentT(float(location), float(scale), float(degrees_of_freedom)))

    location = _net_property("Mu", "The location μ (``Mu``); settable.")
    scale = _net_property("Sigma", "The scale σ > 0 (``Sigma``); settable.")
    degrees_of_freedom = _net_property("DegreesOfFreedom", "The degrees of freedom ν > 0 (``DegreesOfFreedom``); "
                                                           "settable.")


@_register("Triangular")
class Triangular(UnivariateDistribution):
    """Triangular distribution on [minimum, maximum] with mode ``most_likely``.

    Parameters
    ----------
    minimum : float, optional
        Lower end a; default 0.
    most_likely : float, optional
        The mode c; default 0.5.
    maximum : float, optional
        Upper end b; default 1.
    """

    def __init__(self, minimum: float = 0.0, most_likely: float = 0.5, maximum: float = 1.0) -> None:
        super().__init__(_nd().Triangular(float(minimum), float(most_likely), float(maximum)))

    minimum = _net_property("Min", "The lower end a (``Min``); settable.")
    most_likely = _net_property("MostLikely", "The mode c (``MostLikely``); settable.")
    maximum = _net_property("Max", "The upper end b (``Max``); settable.")


class TruncatedDistribution(UnivariateDistribution):
    """Any univariate distribution truncated to [lower_bound, upper_bound] (``TruncatedDistribution``).

    The density is ``f(x) / (F(upper) - F(lower))`` on the interval. ``type_name`` on an instance is
    the base distribution's (Numerics reports the base type), so this class is not in :data:`FAMILIES`.

    Parameters
    ----------
    distribution : UnivariateDistribution
        The base distribution; Numerics keeps its own copy.
    lower_bound : float
        Lower truncation point.
    upper_bound : float
        Upper truncation point, above ``lower_bound``.
    """

    def __init__(self, distribution: UnivariateDistribution, lower_bound: float, upper_bound: float) -> None:
        super().__init__(_nd().TruncatedDistribution(distribution.net, float(lower_bound), float(upper_bound)))

    lower_bound = _net_property("Min", "The lower truncation point (``Min``).", settable=False)
    upper_bound = _net_property("Max", "The upper truncation point (``Max``).", settable=False)

    @property
    def base_distribution(self) -> UnivariateDistribution:
        """The base distribution, wrapped (``BaseDistribution``)."""
        return wrap(self.net.BaseDistribution)

    def _value_key(self) -> tuple:
        """Return the truncated value key extended with the base distribution's family and state."""
        return super()._value_key() + self.base_distribution._value_key()


_CLASS_BY_NET_NAME["TruncatedDistribution"] = TruncatedDistribution


@_register("TruncatedNormal")
class TruncatedNormal(UnivariateDistribution):
    """Normal distribution truncated to [minimum, maximum] (``TruncatedNormal``).

    ``mu`` and ``sigma`` are the parameters of the parent normal; the moment properties
    (``mean``, ``standard_deviation``) are those of the truncated distribution.

    Parameters
    ----------
    mu : float, optional
        Mean of the parent normal; default 0.5.
    sigma : float, optional
        Standard deviation of the parent normal, > 0; default 0.2.
    minimum : float, optional
        Lower truncation point; default 0.
    maximum : float, optional
        Upper truncation point; default 1.
    """

    def __init__(self, mu: float = 0.5, sigma: float = 0.2, minimum: float = 0.0, maximum: float = 1.0) -> None:
        super().__init__(_nd().TruncatedNormal(float(mu), float(sigma), float(minimum), float(maximum)))

    mu = _net_property("Mu", "The mean of the parent normal (``Mu``); settable.")
    sigma = _net_property("Sigma", "The standard deviation of the parent normal (``Sigma``); settable.")
    minimum = _net_property("Min", "The lower truncation point (``Min``); settable.")
    maximum = _net_property("Max", "The upper truncation point (``Max``); settable.")


@_register("UniformDiscrete")
class UniformDiscrete(UnivariateDistribution):
    """Discrete uniform distribution on the integers minimum, ..., maximum.

    Parameters
    ----------
    minimum : float, optional
        Lowest value; default 0.
    maximum : float, optional
        Highest value; default 1.
    """

    def __init__(self, minimum: float = 0.0, maximum: float = 1.0) -> None:
        super().__init__(_nd().UniformDiscrete(float(minimum), float(maximum)))

    minimum = _net_property("Min", "The lowest value (``Min``); settable.")
    maximum = _net_property("Max", "The highest value (``Max``); settable.")
    count = _net_property("N", "The number of values, maximum - minimum + 1 (``N``).", settable=False)


@_register("VonMises")
class VonMises(UnivariateDistribution):
    """Von Mises (circular normal) distribution on [-π, π]; the CDF accumulates from -π.

    Parameters
    ----------
    mu : float, optional
        Mean direction μ in [-π, π]; default 0.
    kappa : float, optional
        Concentration κ >= 0; default 1.
    """

    def __init__(self, mu: float = 0.0, kappa: float = 1.0) -> None:
        super().__init__(_nd().VonMises(float(mu), float(kappa)))

    mu = _net_property("Mu", "The mean direction μ (``Mu``); settable.")
    kappa = _net_property("Kappa", "The concentration κ (``Kappa``); settable.")


@_register("Weibull")
class Weibull(UnivariateDistribution):
    """Two-parameter Weibull distribution: ``F(x) = 1 - exp(-(x / λ)^κ)``.

    Parameters
    ----------
    scale : float, optional
        Scale λ > 0; default 10.
    shape : float, optional
        Shape κ > 0; default 2.
    """

    def __init__(self, scale: float = 10.0, shape: float = 2.0) -> None:
        super().__init__(_nd().Weibull(float(scale), float(shape)))

    scale = _net_property("Lambda", "The scale λ > 0 (``Lambda``); settable.")
    shape = _net_property("Kappa", "The shape κ > 0 (``Kappa``); settable.")

    def solve_mle(self, sample: Any) -> np.ndarray:
        """Return the maximum likelihood ``[scale, shape]`` from the profile-likelihood root (``SolveMLE``).

        Parameters
        ----------
        sample : array_like
            Positive observations.

        Returns
        -------
        numpy.ndarray
            ``[scale, shape]``; this distribution is unchanged.
        """
        return to_numpy_doubles(self.net.SolveMLE(to_net_doubles(sample)))


def _dependency_type() -> Any:
    """Return the ``Probability.DependencyType`` enum type."""
    return namespace("Numerics.Data.Statistics").Probability.DependencyType


@_register("CompetingRisks")
class CompetingRisks(UnivariateDistribution):
    """The minimum (or maximum) of several random variables, the competing-risks failure time.

    With independent components the CDF of the minimum is ``1 - prod(1 - F_i(x))``; ``dependency``
    and ``correlation_matrix`` select a dependent combination. The CDF is tabulated on a grid and
    interpolated (``x_transform``, ``probability_transform``).

    Parameters
    ----------
    distributions : sequence of UnivariateDistribution
        The component distributions.
    """

    def __init__(self, distributions: Sequence[UnivariateDistribution]) -> None:
        super().__init__(_nd().CompetingRisks(_members_array(distributions)))

    def __repr__(self) -> str:
        return f"CompetingRisks({self.distributions!r})"

    @property
    def distributions(self) -> list[UnivariateDistribution]:
        """The component distributions, wrapped (``Distributions``)."""
        return [wrap(component) for component in self.net.Distributions]

    prng_seed = _net_property("PRNGSeed", "The seed of the random draws of dependent combinations (``PRNGSeed``); "
                                          "settable.", kind=int)
    minimum_of_random_variables = _net_property(
        "MinimumOfRandomVariables", "True for the minimum of the components (default), False for the maximum "
                                    "(``MinimumOfRandomVariables``); settable.", kind=bool)
    x_transform = _enum_property("XTransform", lambda: _data().Transform,
                                 "The X interpolation transform of the CDF table; settable.")
    probability_transform = _enum_property("ProbabilityTransform", lambda: _data().Transform,
                                           "The probability interpolation transform of the CDF table; settable.")
    dependency = _enum_property("Dependency", _dependency_type,
                                "The dependency between components (``Probability.DependencyType`` name); settable.")

    @property
    def correlation_matrix(self) -> np.ndarray | None:
        """The correlation matrix of the components for a correlated dependency; settable."""
        matrix = self.net.CorrelationMatrix
        return None if matrix is None else to_numpy_matrix(matrix)

    @correlation_matrix.setter
    def correlation_matrix(self, value: Any) -> None:
        self._ensure_mutable()
        self.net.CorrelationMatrix = None if value is None else to_net_matrix(value)

    def set_parameters(self, parameters: Any) -> None:
        """Replace the components (``SetParameters(UnivariateDistributionBase[])``) or set a parameter vector.

        Parameters
        ----------
        parameters : sequence of UnivariateDistribution or array_like
            New component distributions, or the parameter vector in Numerics order.

        Raises
        ------
        TypeError
            If the distribution is frozen.
        """
        self._ensure_mutable()
        items = list(parameters) if not isinstance(parameters, np.ndarray) else []
        if items and all(isinstance(item, UnivariateDistribution) for item in items):
            self.net.SetParameters(_members_array(items))
        else:
            super().set_parameters(parameters)

    def create_empirical_cdf(self) -> None:
        """Build the tabulated CDF now rather than on first evaluation (``CreateEmpiricalCDF``)."""
        self.net.CreateEmpiricalCDF()

    def cumulative_incidence_functions(self, bins: Any = None) -> list[Empirical]:
        """Return the cumulative incidence function of each component (``CumulativeIncidenceFunctions``).

        Parameters
        ----------
        bins : pynumerics.sampling.StratificationBins, optional
            The X bins of the tables; None, the default, lets Numerics choose them.

        Returns
        -------
        list of Empirical
            One table per component: the probability that it is the first to occur by x.
        """
        functions = self.net.CumulativeIncidenceFunctions(None if bins is None else _bins_net(bins))
        return [wrap(function) for function in functions]

    def generate_random_values_with_dependency(self, sample_size: int, seed: int = -1) -> np.ndarray:
        """Draw the minimum (or maximum) of dependent component draws (``GenerateRandomValuesWithDependency``).

        Parameters
        ----------
        sample_size : int
            Number of values.
        seed : int, optional
            Mersenne Twister seed; zero or negative seeds from the clock.

        Returns
        -------
        numpy.ndarray
            The sample.
        """
        return to_numpy_doubles(self.net.GenerateRandomValuesWithDependency(int(sample_size), int(seed)))


@_register("Mixture")
class Mixture(UnivariateDistribution):
    """Finite mixture of univariate distributions.

    Parameters
    ----------
    weights : array_like
        Component weights.
    distributions : sequence of UnivariateDistribution
        Components, one per weight.
    """

    def __init__(self, weights: Any, distributions: Sequence[UnivariateDistribution]) -> None:
        super().__init__(_nd().Mixture(to_net_doubles(weights), _members_array(distributions)))

    @property
    def weights(self) -> np.ndarray:
        """The component weights (``Weights``)."""
        return to_numpy_doubles(self.net.Weights)

    @property
    def components(self) -> list[UnivariateDistribution]:
        """The component distributions, wrapped (``Distributions``)."""
        return [wrap(component) for component in self.net.Distributions]

    is_zero_inflated = _net_property("IsZeroInflated", "Whether a point mass at zero with weight ``zero_weight`` "
                                                       "joins the components (``IsZeroInflated``); settable.",
                                     kind=bool)
    zero_weight = _net_property("ZeroWeight", "The weight of the point mass at zero (``ZeroWeight``); settable.")
    max_iterations = _net_property("MaxIterations", "The iteration limit of the maximum likelihood fit "
                                                    "(``MaxIterations``); settable.", kind=int)
    tolerance = _net_property("Tolerance", "The convergence tolerance of the maximum likelihood fit (``Tolerance``); "
                                           "settable.")
    iterations = _net_property("Iterations", "The iterations the last fit used (``Iterations``).", kind=int,
                               settable=False)
    x_transform = _enum_property("XTransform", lambda: _data().Transform,
                                 "The X interpolation transform of the inverse-CDF table; settable.")
    probability_transform = _enum_property("ProbabilityTransform", lambda: _data().Transform,
                                           "The probability interpolation transform of the inverse-CDF table; "
                                           "settable.")

    def create_empirical_cdf(self) -> None:
        """Build the tabulated CDF used by ``inverse_cdf`` now rather than on first use (``CreateEmpiricalCDF``)."""
        self.net.CreateEmpiricalCDF()


@_register("EmpiricalDistribution", "Empirical")
class Empirical(UnivariateDistribution):
    """Empirical distribution defined by value and probability tables.

    With the default orders the probabilities are non-exceedance probabilities
    (strictly ascending). With ``probability_order="descending"`` the table is a
    survival function: the probabilities are exceedance probabilities S(x) (strictly
    descending), ``cdf`` returns 1 − S(x) and ``inverse_cdf(q)`` reads the exceedance
    axis at 1 − q. Interpolation is in (X, Φ⁻¹(P)) space by default
    (``probability_transform="normal_z"``) and holds the endpoints.

    Parameters
    ----------
    values : array_like
        Ascending X values (ties allowed).
    probabilities : array_like
        Probability values in [0, 1], strictly monotone in ``probability_order``.
    x_order : SortOrder or str, optional
        Sort order of ``values``; must be ascending (the default).
    probability_order : SortOrder or str, optional
        ``"ascending"`` (non-exceedance, default) or ``"descending"`` (survival).

    Notes
    -----
    Uses ``EmpiricalDistribution(IList, IList, SortOrder, SortOrder)``.
    """

    def __init__(self, values: Any, probabilities: Any, x_order: enums.SortOrder | str = enums.SortOrder.ASCENDING,
                 probability_order: enums.SortOrder | str = enums.SortOrder.ASCENDING) -> None:
        nd = _nd()
        data = _data()
        super().__init__(nd.EmpiricalDistribution(to_net_doubles(values), to_net_doubles(probabilities),
                                                  net_enum(data.SortOrder, x_order),
                                                  net_enum(data.SortOrder, probability_order)))

    @classmethod
    def survival(cls, values: Any, exceedance: Any) -> Empirical:
        """Build a survival-orientation table from exceedance probabilities.

        Parameters
        ----------
        values : array_like
            Ascending X values.
        exceedance : array_like
            Strictly descending exceedance probabilities S(x).

        Returns
        -------
        Empirical
            The distribution with ``probability_order`` descending.
        """
        return cls(values, exceedance, enums.SortOrder.ASCENDING, enums.SortOrder.DESCENDING)

    @classmethod
    def from_sample(cls, sample: Any,
                    plotting_position: enums.PlottingPosition | str = enums.PlottingPosition.WEIBULL) -> Empirical:
        """Build the distribution of a sample with plotting-position probabilities.

        Parameters
        ----------
        sample : array_like
            Observations; sorted by Numerics.
        plotting_position : PlottingPosition or str, optional
            Plotting-position formula; default Weibull.

        Returns
        -------
        Empirical
            ``EmpiricalDistribution(IList, PlottingPostionType)``.
        """
        nd = _nd()
        stats = namespace("Numerics.Data.Statistics")
        kind = net_enum(stats.PlottingPositions.PlottingPostionType, plotting_position)
        instance = cls.__new__(cls)
        UnivariateDistribution.__init__(instance, nd.EmpiricalDistribution(to_net_doubles(sample), kind))
        return instance

    @classmethod
    def from_paired_data(cls, paired_data: Any) -> Empirical:
        """Build the distribution from ordered paired data, keeping its orders and strictness.

        Parameters
        ----------
        paired_data : pynumerics.data.OrderedPairedData or Numerics.Data.OrderedPairedData
            X ascending; Y the probabilities.

        Returns
        -------
        Empirical
            ``EmpiricalDistribution(OrderedPairedData)``.
        """
        net = getattr(paired_data, "net", paired_data)
        instance = cls.__new__(cls)
        UnivariateDistribution.__init__(instance, _nd().EmpiricalDistribution(net))
        return instance

    @property
    def x_values(self) -> np.ndarray:
        """The X table."""
        return to_numpy_doubles(self.net.XValues)

    @property
    def probability_values(self) -> np.ndarray:
        """The probability table as stored (exceedance for a survival table)."""
        return to_numpy_doubles(self.net.ProbabilityValues)

    @property
    def probability_order(self) -> str:
        """The probability sort order: ``"Ascending"``, ``"Descending"`` or ``"None"``."""
        return str(self.net.ProbabilityOrder)

    @property
    def x_transform(self) -> str:
        """The X interpolation transform; settable with a :class:`~pynumerics.enums.Transform` or name."""
        return str(self.net.XTransform)

    @x_transform.setter
    def x_transform(self, value: enums.Transform | str) -> None:
        self.net.XTransform = net_enum(_data().Transform, value)

    @property
    def probability_transform(self) -> str:
        """The probability interpolation transform (default ``NormalZ``); settable."""
        return str(self.net.ProbabilityTransform)

    @probability_transform.setter
    def probability_transform(self, value: enums.Transform | str) -> None:
        self.net.ProbabilityTransform = net_enum(_data().Transform, value)

    @property
    def extrapolation(self) -> str:
        """The X-space extrapolation sides (default ``None``: endpoint hold); settable."""
        return str(self.net.Extrapolation)

    @extrapolation.setter
    def extrapolation(self, value: enums.Extrapolation | str) -> None:
        self.net.Extrapolation = net_enum(_data().ExtrapolationSides, value)

    @property
    def x_value_order(self) -> str:
        """The X sort order: ``"Ascending"`` (``XValueOrder``)."""
        return str(self.net.XValueOrder)

    @property
    def strict_x_values(self) -> bool:
        """Whether the X values are strictly monotone, with no ties (``StrictXValues``)."""
        return bool(self.net.StrictXValues)

    @property
    def strict_probabilities(self) -> bool:
        """Whether the probabilities are strictly monotone (``StrictProbabilities``)."""
        return bool(self.net.StrictProbabilities)

    @classmethod
    def convolve(cls, distributions: Sequence[Empirical], number_of_points: int = 1024,
                 log_spaced_output: bool = False) -> Empirical:
        """Return the distribution of the sum of independent empirical variables (``Convolve``).

        Parameters
        ----------
        distributions : sequence of Empirical
            Two or more distributions; ``log_spaced_output`` needs exactly two.
        number_of_points : int, optional
            Points of the result table; default 1024.
        log_spaced_output : bool, optional
            Space the output X values logarithmically (the two-distribution overload); default False.

        Returns
        -------
        Empirical
            The table of the sum.
        """
        empirical = _nd().EmpiricalDistribution
        items = list(distributions)
        if log_spaced_output:
            if len(items) != 2:
                raise ValueError("log_spaced_output needs exactly two distributions.")
            result = empirical.Convolve(items[0].net, items[1].net, int(number_of_points), True)
        elif len(items) == 2:
            result = empirical.Convolve(items[0].net, items[1].net, int(number_of_points))
        else:
            from System.Collections.Generic import List  # noqa: PLC0415 - the runtime is running

            members = List[empirical]()
            for item in items:
                members.Add(item.net)
            result = empirical.Convolve(members, int(number_of_points))
        return wrap(result)

    @staticmethod
    def convolve_discrete(values1: Any, masses1: Any, values2: Any, masses2: Any,
                          lattice_points: int) -> tuple[np.ndarray, np.ndarray]:
        """Convolve two discrete distributions on a common lattice (``ConvolveDiscrete``).

        Parameters
        ----------
        values1, masses1 : array_like
            Support points and probability masses of the first variable.
        values2, masses2 : array_like
            Support points and probability masses of the second variable.
        lattice_points : int
            Number of lattice points of the result.

        Returns
        -------
        values : numpy.ndarray
            The lattice of the sum.
        masses : numpy.ndarray
            The probability masses on the lattice.
        """
        values, masses = _nd().EmpiricalDistribution.ConvolveDiscrete(
            to_net_doubles(values1), to_net_doubles(masses1), to_net_doubles(values2), to_net_doubles(masses2),
            int(lattice_points), None, None)
        return to_numpy_doubles(values), to_numpy_doubles(masses)


@_register("KernelDensity")
class KernelDensity(UnivariateDistribution):
    """Kernel density estimate of a sample (``Numerics.Distributions.KernelDensity``).

    The density is ``f(x) = sum_i w_i K((x - x_i) / h) / (h sum_i w_i)`` (unit weights when none
    are given). The CDF and inverse CDF interpolate a 1000-bin cumulative table over the support,
    which is the sample range when ``bounded_by_data`` is True (the default) and the range widened
    by ``3 h`` on each side otherwise.

    Parameters
    ----------
    sample : array_like
        Finite sample values, in any order.
    kernel : KernelType or str, optional
        The kernel; default Gaussian.
    bandwidth : float, optional
        The bandwidth ``h`` (positive). None, the default, applies the Numerics rule of thumb
        ``(4 / (3 n))^(1/5) s``: ``s`` is the sample standard deviation (divisor ``n - 1``), or
        with weights the weighted standard deviation with divisor ``sum_i w_i``.
    weights : array_like, optional
        Non-negative weights, one per value, with a positive total.

    Notes
    -----
    Uses ``KernelDensity(IList, KernelType[, double])``, or
    ``KernelDensity(IList, IList, KernelType, double?)`` with weights.
    """

    def __init__(self, sample: Any, kernel: enums.KernelType | str = enums.KernelType.GAUSSIAN,
                 bandwidth: float | None = None, weights: Any = None) -> None:
        kd = _nd().KernelDensity
        kind, values = net_enum(kd.KernelType, kernel), to_net_doubles(sample)
        if weights is not None:
            net = kd(values, to_net_doubles(weights), kind, None if bandwidth is None else float(bandwidth))
        elif bandwidth is None:
            net = kd(values, kind)
        else:
            net = kd(values, kind, float(bandwidth))
        super().__init__(net)

    def __repr__(self) -> str:
        return f"KernelDensity(kernel={self.kernel!r}, bandwidth={self.bandwidth!r}, sample_size={self.sample_size})"

    @property
    def sample_data(self) -> np.ndarray:
        """The sample values, in input order."""
        return to_numpy_doubles(self.net.SampleData)

    @property
    def sample_size(self) -> int:
        """The number of sample values."""
        return int(self.net.SampleSize)

    @property
    def kernel(self) -> str:
        """The kernel name (``KernelDistribution``); settable with a :class:`~pynumerics.enums.KernelType` or name."""
        return str(self.net.KernelDistribution)

    @kernel.setter
    def kernel(self, value: enums.KernelType | str) -> None:
        self.net.KernelDistribution = net_enum(_nd().KernelDensity.KernelType, value)

    @property
    def bandwidth(self) -> float:
        """The bandwidth ``h``; settable with a positive finite value."""
        return float(self.net.Bandwidth)

    @bandwidth.setter
    def bandwidth(self, value: float) -> None:
        self.net.Bandwidth = float(value)

    @property
    def bounded_by_data(self) -> bool:
        """Whether the support is the sample range (default) rather than the range widened by ``3 h``; settable."""
        return bool(self.net.BoundedByData)

    @bounded_by_data.setter
    def bounded_by_data(self, value: bool) -> None:
        self.net.BoundedByData = bool(value)

    @property
    def x_transform(self) -> str:
        """The X interpolation transform of the CDF table (default ``None``); settable."""
        return str(self.net.XTransform)

    @x_transform.setter
    def x_transform(self, value: enums.Transform | str) -> None:
        self.net.XTransform = net_enum(_data().Transform, value)

    @property
    def probability_transform(self) -> str:
        """The probability interpolation transform of the CDF table (default ``NormalZ``); settable."""
        return str(self.net.ProbabilityTransform)

    @probability_transform.setter
    def probability_transform(self, value: enums.Transform | str) -> None:
        self.net.ProbabilityTransform = net_enum(_data().Transform, value)

    def bandwidth_rule(self, sample: Any, weights: Any = None) -> float:
        """Return the rule-of-thumb bandwidth of a sample (``BandwidthRule``).

        Parameters
        ----------
        sample : array_like
            Finite sample values.
        weights : array_like, optional
            Weights, one per value; the weighted standard deviation is used.

        Returns
        -------
        float
            ``(4 / (3 n))^(1/5)`` times the standard deviation (as for the ``bandwidth`` default),
            kept finite and positive for a constant sample.
        """
        values = to_net_doubles(sample)
        if weights is None:
            return float(self.net.BandwidthRule(values))
        return float(self.net.BandwidthRule(values, to_net_doubles(weights)))

    def set_sample_data(self, sample: Any, weights: Any = None) -> None:
        """Replace the sample (and weights) in place, keeping the kernel and bandwidth (``SetSampleData``).

        Parameters
        ----------
        sample : array_like
            Finite sample values.
        weights : array_like, optional
            Non-negative weights, one per value; None, the default, gives unit weights.

        Raises
        ------
        TypeError
            If the distribution is frozen.
        """
        self._ensure_mutable()
        if weights is None:
            self.net.SetSampleData(to_net_doubles(sample))
        else:
            self.net.SetSampleData(to_net_doubles(sample), to_net_doubles(weights))

    degenerate_relative_bandwidth = _FamilyAttribute(
        lambda cls: float(_nd().KernelDensity.DegenerateRelativeBandwidth),
        lambda self: float(_nd().KernelDensity.DegenerateRelativeBandwidth),
        """The bandwidth of a constant sample relative to its magnitude (``DegenerateRelativeBandwidth``).""")

    degenerate_absolute_bandwidth = _FamilyAttribute(
        lambda cls: float(_nd().KernelDensity.DegenerateAbsoluteBandwidth),
        lambda self: float(_nd().KernelDensity.DegenerateAbsoluteBandwidth),
        """The bandwidth of a constant sample at zero (``DegenerateAbsoluteBandwidth``).""")


def create_distribution(type_name: str, parameters: Any = None) -> UnivariateDistribution:
    """Create a distribution by type name through ``UnivariateDistributionFactory``.

    ``CreateDistribution(UnivariateDistributionType)`` builds the family at its default
    parameters; ``parameters``, when given, are then set through ``SetParameters``.
    Invalid parameters are recorded, not refused (``parameters_valid`` is False).

    Parameters
    ----------
    type_name : str
        A ``UnivariateDistributionType`` member name, for example ``"Gumbel"``; matched
        case-insensitively with underscores ignored.
    parameters : array_like, optional
        The parameters in Numerics order.

    Returns
    -------
    UnivariateDistribution
        The wrapped distribution (its specific wrapper class when one exists).

    Raises
    ------
    ValueError
        If ``type_name`` is not a member of ``UnivariateDistributionType``, or ``parameters`` does
        not hold the family's number of parameters (the message names the expected count).
    """
    nd = _nd()
    net = unwrap(nd.UnivariateDistributionFactory.CreateDistribution(
        net_enum(nd.UnivariateDistributionType, type_name)))
    if parameters is not None:
        net.SetParameters(to_net_doubles(_check_arity(net, parameters)))
    return wrap(net)


def standard_z(probability: Any) -> Any:
    """Return standard normal quantiles through ``Normal.StandardZ``.

    Parameters
    ----------
    probability : float or array_like
        Non-exceedance probabilities; 0 and 1 map to ±``double.Epsilon`` probabilities
        as in Numerics.

    Returns
    -------
    float or numpy.ndarray
        Standard normal deviates, shaped like ``probability``.
    """
    n = _nd().Normal
    return evaluate(probability, _static("StandardZ", n), _static_list("StandardZ", n))


def standard_cdf(z: Any) -> Any:
    """Return the standard normal CDF through ``Normal.StandardCDF``.

    Parameters
    ----------
    z : float or array_like
        Standard normal deviates.

    Returns
    -------
    float or numpy.ndarray
        Non-exceedance probabilities, shaped like ``z``.
    """
    n = _nd().Normal
    return evaluate(z, _static("StandardCDF", n), _static_list("StandardCDF", n))


_STATIC_BOUND: dict[str, Any] = {}


def _static(name: str, owner: Any) -> Any:
    """Return a static member's ``double`` overload, bound once per process (see ``_util.bind``)."""
    member = _STATIC_BOUND.get(name)
    if member is None:
        member = _STATIC_BOUND[name] = bind(getattr(owner, name))
    return member


_STATIC_LIST_BOUND: dict[str, Any] = {}


def _static_list(name: str, owner: Any) -> Any:
    """Return a static member's ``IList<double>`` overload, bound once per process (see ``_util.bind_list``)."""
    member = _STATIC_LIST_BOUND.get(name)
    if member is None:
        member = _STATIC_LIST_BOUND[name] = bind_list(getattr(owner, name))
    return member


def standard_pdf(z: Any) -> Any:
    """Return the standard normal density through ``Normal.StandardPDF``.

    Parameters
    ----------
    z : float or array_like
        Standard normal deviates.

    Returns
    -------
    float or numpy.ndarray
        Densities, shaped like ``z``.
    """
    n = _nd().Normal
    return evaluate(z, _static("StandardPDF", n), _static_list("StandardPDF", n))


def try_create_distribution(type_name: str) -> UnivariateDistribution | None:
    """Create a family at its default parameters, or return None when Numerics cannot (``TryCreateDistribution``).

    Parameters
    ----------
    type_name : str
        A ``UnivariateDistributionType`` member name, matched case-insensitively with underscores ignored.

    Returns
    -------
    UnivariateDistribution or None
        The wrapped distribution; None for the families that need components (``Mixture``,
        ``CompetingRisks``) and ``UserDefined``.

    Raises
    ------
    ValueError
        If ``type_name`` is not a member of ``UnivariateDistributionType``.
    """
    nd = _nd()
    created, net = nd.UnivariateDistributionFactory.TryCreateDistribution(
        net_enum(nd.UnivariateDistributionType, type_name), None)
    return wrap(net) if created else None


def from_xml(xml: str) -> UnivariateDistribution:
    """Rebuild a distribution from the XML :meth:`UnivariateDistribution.to_xml` writes.

    Calls ``UnivariateDistributionFactory.CreateDistribution(XElement)``, which dispatches to the
    family's own ``FromXElement`` for mixtures, competing risks, percentile PERT, empirical and
    kernel density tables.

    Parameters
    ----------
    xml : str
        The ``Distribution`` element as text.

    Returns
    -------
    UnivariateDistribution
        The wrapped distribution.

    Raises
    ------
    ValueError
        If the type or a parameter is missing or invalid.
    """
    return wrap(_nd().UnivariateDistributionFactory.CreateDistribution(_parse_xml(xml)))


def standard_log_pdf(z: Any) -> Any:
    """Return the log of the standard normal density through ``Normal.StandardLogPDF``.

    Parameters
    ----------
    z : float or array_like
        Standard normal deviates.

    Returns
    -------
    float or numpy.ndarray
        Log densities, shaped like ``z``; Numerics has no array overload, so an array is evaluated
        element by element with the bound scalar member.
    """
    return loop(z, _static("StandardLogPDF", _nd().Normal))


def _elementwise(*arrays: Any) -> tuple[tuple[int, ...], list[np.ndarray]]:
    """Broadcast the per-element arguments of a static array member and flatten them."""
    broadcast = np.broadcast_arrays(*(np.asarray(a) for a in arrays))
    shape = broadcast[0].shape
    return shape, [np.ascontiguousarray(b).ravel() for b in broadcast]


def _shaped(values: Any, shape: tuple[int, ...]) -> Any:
    """Return a .NET result array as a float for a scalar shape, otherwise as an array of ``shape``."""
    out = to_numpy_doubles(values)
    return float(out[0]) if shape == () else out.reshape(shape)


def _poisson_static(name: str, rate: Any, values: Any) -> Any:
    """Call the per-element static ``Poisson.<name>(double[], double[])`` on broadcast arguments."""
    shape, (rates, points) = _elementwise(np.asarray(rate, dtype=np.float64), np.asarray(values, dtype=np.float64))

    def factory() -> Any:
        poisson = _nd().Poisson  # starts the runtime before System is imported
        import System  # noqa: PLC0415 - the runtime is running

        doubles = System.Array[System.Double]
        return getattr(poisson, name).Overloads[doubles, doubles]

    member = cached_overload(f"Poisson.{name}(double[], double[])", factory)
    return _shaped(member(to_net_doubles(rates), to_net_doubles(points)), shape)


def poisson_inverse_cdf(rate: Any, probability: Any) -> Any:
    """Return Poisson quantiles with one rate per probability (``Poisson.InverseCDF(double[], double[])``).

    Parameters
    ----------
    rate : float or array_like
        Rates λ > 0, broadcast against ``probability``.
    probability : float or array_like
        Non-exceedance probabilities in [0, 1].

    Returns
    -------
    float or numpy.ndarray
        The smallest integer k with P(X <= k) >= p for each element, shaped like the broadcast
        arguments; probability 0 gives 0 and 1 gives infinity. Each value is bitwise equal to
        ``Poisson(rate).inverse_cdf(probability)``.

    Raises
    ------
    ValueError
        If a rate is not positive and finite, a probability is outside [0, 1], or the arguments
        do not broadcast.
    """
    return _poisson_static("InverseCDF", rate, probability)


def poisson_log_pmf(rate: Any, k: Any) -> Any:
    """Return Poisson log probability masses with one rate per count (``Poisson.LogPDF(double[], double[])``).

    Parameters
    ----------
    rate : float or array_like
        Rates λ > 0, broadcast against ``k``.
    k : float or array_like
        Counts; a value that is not a nonnegative integer has log mass ``-inf``.

    Returns
    -------
    float or numpy.ndarray
        The log masses, shaped like the broadcast arguments, each bitwise equal to
        ``Poisson(rate).log_pdf(k)``.

    Raises
    ------
    ValueError
        If a rate is not positive and finite, or the arguments do not broadcast.
    """
    return _poisson_static("LogPDF", rate, k)


def poisson_cdf(rate: Any, k: Any) -> Any:
    """Return Poisson cumulative probabilities with one rate per count (``Poisson.CDF(double[], double[])``).

    Parameters
    ----------
    rate : float or array_like
        Rates λ > 0, broadcast against ``k``.
    k : float or array_like
        Counts, floored when not integers.

    Returns
    -------
    float or numpy.ndarray
        P(X <= k), shaped like the broadcast arguments, each bitwise equal to ``Poisson(rate).cdf(k)``.

    Raises
    ------
    ValueError
        If a rate is not positive and finite, or the arguments do not broadcast.
    """
    return _poisson_static("CDF", rate, k)


def _binomial_static(name: str, probability: Any, number_of_trials: Any, k: Any) -> Any:
    """Call the per-element static ``Binomial.<name>(double[], int[], double[])`` on broadcast arguments."""
    trials = np.asarray(number_of_trials)
    if trials.dtype.kind not in "iuf" or (trials.dtype.kind == "f" and not np.array_equal(trials, np.trunc(trials))):
        raise ValueError("The number of trials must be an integer.")
    if trials.size and (trials.min() < np.iinfo(np.int32).min or trials.max() > np.iinfo(np.int32).max):
        raise ValueError("The number of trials must fit in a 32-bit integer.")
    shape, (p, n, points) = _elementwise(np.asarray(probability, dtype=np.float64), trials,
                                         np.asarray(k, dtype=np.float64))

    def factory() -> Any:
        binomial = _nd().Binomial  # starts the runtime before System is imported
        import System  # noqa: PLC0415 - the runtime is running

        doubles = System.Array[System.Double]
        return getattr(binomial, name).Overloads[doubles, System.Array[System.Int32], doubles]

    member = cached_overload(f"Binomial.{name}(double[], int[], double[])", factory)
    return _shaped(member(to_net_doubles(p), to_net_ints(n), to_net_doubles(points)), shape)


def binomial_log_pmf(probability: Any, number_of_trials: Any, k: Any) -> Any:
    """Return binomial log probability masses with one parameter set per count (``Binomial.LogPDF``).

    Calls ``Binomial.LogPDF(double[], int[], double[])``.

    Parameters
    ----------
    probability : float or array_like
        Success probabilities in [0, 1], broadcast against the other arguments.
    number_of_trials : int or array_like
        Numbers of trials n >= 1.
    k : float or array_like
        Numbers of successes; a value that is not an integer in [0, n] has log mass ``-inf``.

    Returns
    -------
    float or numpy.ndarray
        The log masses, shaped like the broadcast arguments, each bitwise equal to
        ``Binomial(probability, number_of_trials).log_pdf(k)``.

    Raises
    ------
    ValueError
        If a probability is outside [0, 1], a number of trials is not a positive integer that fits
        in 32 bits (booleans, text and fractional values are refused), or the arguments do not broadcast.
    """
    return _binomial_static("LogPDF", probability, number_of_trials, k)


def binomial_cdf(probability: Any, number_of_trials: Any, k: Any) -> Any:
    """Return binomial cumulative probabilities with one parameter set per count (``Binomial.CDF``).

    Calls ``Binomial.CDF(double[], int[], double[])``.

    Parameters
    ----------
    probability : float or array_like
        Success probabilities in [0, 1], broadcast against the other arguments.
    number_of_trials : int or array_like
        Numbers of trials n >= 1.
    k : float or array_like
        Numbers of successes, floored when not integers.

    Returns
    -------
    float or numpy.ndarray
        P(X <= k), shaped like the broadcast arguments, each bitwise equal to
        ``Binomial(probability, number_of_trials).cdf(k)``.

    Raises
    ------
    ValueError
        If a probability is outside [0, 1], a number of trials is not a positive integer that fits
        in 32 bits (booleans, text and fractional values are refused), or the arguments do not broadcast.
    """
    return _binomial_static("CDF", probability, number_of_trials, k)


__all__ = [
    "FAMILIES", "FrozenDistribution", "Bernoulli", "Beta", "Binomial", "Cauchy", "ChiSquared", "CompetingRisks",
    "Deterministic", "Empirical", "Exponential", "GammaDistribution", "GeneralizedBeta", "GeneralizedExtremeValue",
    "GeneralizedLogistic", "GeneralizedNormal", "GeneralizedPareto", "Geometric", "Gumbel", "InverseChiSquared",
    "InverseGamma", "KappaFour", "KernelDensity", "LnNormal", "LogNormal", "LogPearsonTypeIII", "Logistic", "Mixture",
    "NoncentralT", "Normal", "Pareto", "PearsonTypeIII", "Pert", "PertPercentile", "PertPercentileZ", "Poisson",
    "Rayleigh", "StudentT", "Triangular", "TruncatedDistribution", "TruncatedNormal", "Uniform", "UniformDiscrete",
    "UnivariateDistribution", "VonMises", "Weibull", "binomial_cdf", "binomial_log_pmf", "create_distribution",
    "from_xml", "poisson_cdf", "poisson_inverse_cdf", "poisson_log_pmf", "standard_cdf", "standard_log_pdf",
    "standard_pdf", "standard_z", "try_create_distribution", "wrap",
]

translate_public(globals())
