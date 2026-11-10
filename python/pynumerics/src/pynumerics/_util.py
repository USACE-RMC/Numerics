"""Shared helpers for the typed wrappers: enum mapping, array dispatch and unwrapping."""

from __future__ import annotations

import enum
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

from pynumerics.errors import translated
from pynumerics.marshal import _flat_doubles, to_numpy_doubles
from pynumerics.runtime import initialize


def unwrap(obj: Any) -> Any:
    """Return the concrete .NET object behind a pythonnet interface or base-class wrapper.

    pythonnet 3 wraps return values at their declared type; ``__implementation__``
    exposes the runtime type so members of the concrete class are reachable.
    """
    return getattr(obj, "__implementation__", obj)


_ENUM_MEMBERS: dict[tuple[Any, str], Any] = {}


def fold_name(name: str) -> str:
    """Return the folded form of an enum name: underscores removed, lower case.

    ``"normal_z"``, ``"NORMAL_Z"`` and ``"NormalZ"`` all fold to ``"normalz"``. :func:`net_enum` and
    the ``_missing_`` lookup of the pynumerics enums match names by this rule.
    """
    return name.replace("_", "").lower()


def net_enum(net_type: Any, value: Any) -> Any:
    """Map a Python enum member or a name to a member of a .NET enum.

    Names match case-insensitively and ignore underscores, so ``"normal_z"``,
    ``"NormalZ"`` and ``Transform.NORMAL_Z`` all select ``Transform.NormalZ``.
    A resolved name is remembered per enum type: resolving it again through
    ``System.Enum.GetNames`` costs about fifteen microseconds a call, more than many of
    the members the value is passed to.

    Raises
    ------
    ValueError
        If the name matches no member of ``net_type``.
    """
    if isinstance(value, enum.Enum):
        value = value.value
    if not isinstance(value, str):
        return value
    key = fold_name(value)
    member = _ENUM_MEMBERS.get((net_type, key))
    if member is not None:
        return member
    import clr  # noqa: PLC0415
    import System  # noqa: PLC0415 - the runtime is started by the caller

    clr_type = clr.GetClrType(net_type)
    for name in System.Enum.GetNames(clr_type):
        if str(name).lower() == key:
            member = _ENUM_MEMBERS[(net_type, key)] = getattr(net_type, str(name))
            return member
    names = ", ".join(str(name) for name in System.Enum.GetNames(clr_type))
    raise ValueError(f"'{value}' is not a member of {clr_type.Name} ({names}).")


@translated
def evaluate(values: Any, scalar: Callable[[float], float], vector: Callable[[Any], Any]) -> Any:
    """Evaluate a Numerics member on a scalar or on an array through its array overload.

    A scalar input calls ``scalar`` and returns a Python float. Any other input is
    flattened to float64, passed to ``vector`` as one ``System.Double[]`` (the
    ``IList<double>`` overload), and the result is reshaped to the input's shape.
    Callers pass ``vector`` already bound with :func:`bind_list` where the member is a
    method group, so the overload is not resolved again on every call.
    """
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 0:
        return float(scalar(float(array)))
    result = to_numpy_doubles(vector(_flat_doubles(array)))
    return result if array.ndim == 1 else result.reshape(array.shape)


@translated
def loop(values: Any, scalar: Callable[[float], float]) -> Any:
    """Evaluate a scalar-only Numerics member element by element.

    Used only where Numerics has no array overload; the arithmetic stays in Numerics.
    """
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 0:
        return float(scalar(float(array)))
    out = np.empty(array.size, dtype=np.float64)
    for index, value in enumerate(array.ravel()):
        out[index] = scalar(float(value))
    return out.reshape(array.shape)


def bind(member: Any, arity: int = 1) -> Any:
    """Return the overload of a .NET method group that takes ``arity`` doubles.

    pythonnet resolves an overloaded method group on every call, which costs several times
    the call itself for the cheap scalar members (``InverseCDF(double)`` beside
    ``InverseCDF(IList<double>)``); binding the overload once makes repeated scalar calls
    about six times faster. The arithmetic is unchanged: the bound overload is the member
    the group would have resolved to. A group without that overload is returned as it is.
    """
    import System  # noqa: PLC0415 - the runtime is started by the caller

    try:
        return member.Overloads[(System.Double,) * arity] if arity > 1 else member.Overloads[System.Double]
    except (KeyError, TypeError, AttributeError):
        return member


def bind_list(member: Any) -> Any:
    """Return the ``IList<double>`` overload of a .NET method group.

    pythonnet resolves an overloaded method group on every call. For the array members
    (``InverseCDF(IList<double>)`` beside ``InverseCDF(double)``) the resolution costs about
    30 microseconds a call, several times the evaluation of a short array; binding the
    overload once removes it. The arithmetic is unchanged: the bound overload is the member
    the group resolves to for a ``System.Double[]`` argument. A group without that overload
    is returned as it is.
    """
    import System  # noqa: PLC0415 - the runtime is started by the caller
    from System.Collections.Generic import IList  # noqa: PLC0415

    try:
        return member.Overloads[IList[System.Double]]
    except (KeyError, TypeError, AttributeError):
        return member


def list_overload(member: Any) -> Any:
    """Return the ``IList<double>`` overload of a .NET method group, or None when it has none.

    Like :func:`bind_list`, but a group without the array overload gives None, so a wrapper
    can evaluate an array element by element with the scalar member instead.
    """
    bound = bind_list(member)
    return None if bound is member else bound


_OVERLOADS: dict[str, Any] = {}


def cached_overload(key: str, factory: Callable[[], Any]) -> Any:
    """Return a bound .NET overload, resolving it with ``factory`` on first use only.

    pythonnet resolves an overloaded method group on every call; for a group with many
    overloads (``Statistics.Percentile`` has eight) the resolution costs 40 to 70
    microseconds, several times the call itself. Binding the exact overload once per
    process, under a key naming the member and its signature, removes that cost without
    changing which member runs.
    """
    member = _OVERLOADS.get(key)
    if member is None:
        member = _OVERLOADS[key] = factory()
    return member


_NAMESPACES: dict[str, Any] = {}


def namespace(name: str) -> Any:
    """Start the runtime and import a .NET namespace by dotted name, once per process.

    Wrappers call this on every evaluation; returning the imported module from a dictionary
    skips the import machinery, which costs a few microseconds a call.
    """
    module = _NAMESPACES.get(name)
    if module is None:
        initialize()
        import importlib  # noqa: PLC0415

        module = _NAMESPACES[name] = importlib.import_module(name)
    return module


def _nd() -> Any:
    """Return the ``Numerics.Distributions`` namespace after starting the runtime."""
    return namespace("Numerics.Distributions")


def _data() -> Any:
    """Return the ``Numerics.Data`` namespace after starting the runtime."""
    return namespace("Numerics.Data")


def _opt() -> Any:
    """Return the ``Numerics.Mathematics.Optimization`` namespace after starting the runtime."""
    return namespace("Numerics.Mathematics.Optimization")


def overload(member: Any, *types: str) -> Any:
    """Return the overload of a .NET method group taking the named types.

    Each name is a ``System`` type (``"Double"``, ``"Int32"``, ``"Boolean"``), an array of one
    (``"Double[]"``) or a list of one (``"IList[Double]"``).
    """
    import System  # noqa: PLC0415 - the caller started the runtime
    from System.Collections.Generic import IList  # noqa: PLC0415

    resolved = []
    for name in types:
        if name.startswith("IList["):
            resolved.append(IList[getattr(System, name[6:-1])])
        elif name.endswith("[]"):
            resolved.append(System.Array[getattr(System, name[:-2])])
        else:
            resolved.append(getattr(System, name))
    return member.Overloads[tuple(resolved)] if len(resolved) > 1 else member.Overloads[resolved[0]]


@translated
def broadcast_loop(member: Callable[..., Any], *args: Any, kinds: Sequence[type] | None = None) -> Any:
    """Evaluate a scalar-only Numerics member of several arguments over their broadcast shape.

    Each argument is converted with its entry of ``kinds`` (``float`` by default, ``int`` for
    integer parameters such as a Bessel order). Scalars give a Python float; arrays give a float64
    array of the broadcast shape. Used only where Numerics has no array overload; every value is
    the scalar member's.
    """
    kinds = kinds or (float,) * len(args)
    arrays = np.broadcast_arrays(*(np.asarray(a, dtype=np.int64 if k is int else np.float64)
                                   for a, k in zip(args, kinds, strict=True)))
    if arrays[0].ndim == 0:
        return float(member(*(k(a) for a, k in zip(arrays, kinds, strict=True))))
    out = np.empty(arrays[0].size, dtype=np.float64)
    for index, values in enumerate(zip(*(a.ravel().tolist() for a in arrays), strict=True)):
        out[index] = member(*values)
    return out.reshape(arrays[0].shape)
