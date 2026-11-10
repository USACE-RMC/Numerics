"""Map the exceptions Numerics throws onto the Python exceptions callers expect.

pythonnet raises a .NET exception as a Python exception of the .NET type, which ``except ValueError``
does not catch. Every public wrapper in this package translates them, keeping the .NET message
text exactly and chaining the .NET exception as the cause:

========================================================  ========================
.NET exception                                            Python exception
========================================================  ========================
``ArgumentException`` (``ArgumentOutOfRangeException``)   ``ValueError``
``OverflowException``                                     ``OverflowError``
``DivideByZeroException``                                 ``ZeroDivisionError``
``ArithmeticException`` (any other)                       ``ArithmeticError``
``NotImplementedException``, ``NotSupportedException``    ``NotImplementedError``
any other ``System.Exception``                            ``RuntimeError``
========================================================  ========================

A ``TargetInvocationException``, or an ``AggregateException`` with one inner exception, is
translated by its inner exception. A Python exception, including one raised inside a Python
callback that Numerics invoked, propagates unchanged.

:func:`net_errors` applies the translation to a block of code, :func:`translated` to a function,
and :func:`translate_public` to every public function and class of a module.
"""

from __future__ import annotations

import enum
import functools
import inspect
from collections.abc import Callable, Iterator, MutableMapping
from contextlib import contextmanager
from typing import Any, TypeVar

from pynumerics import runtime

F = TypeVar("F", bound=Callable[..., Any])

_MARKER = "__net_errors__"

_WRAPPED_DUNDERS = frozenset({"__init__", "__call__", "__getitem__", "__setitem__", "__len__", "__contains__",
                              "__iter__", "__eq__", "__hash__"})
"""Special methods that call into .NET in the wrappers and so are translated as well."""


def _system() -> Any:
    """Return the ``System`` namespace, or None while the runtime is not running."""
    if not runtime.is_initialized():
        return None
    import System  # noqa: PLC0415 - the runtime is running

    return System


def translate(error: BaseException) -> BaseException | None:
    """Return the Python exception for a .NET exception.

    Parameters
    ----------
    error : BaseException
        The raised exception.

    Returns
    -------
    BaseException or None
        A new ``ValueError``, ``OverflowError``, ``ZeroDivisionError``, ``ArithmeticError``,
        ``NotImplementedError`` or ``RuntimeError`` carrying the .NET message, or None when
        ``error`` is not a .NET exception.
    """
    system = _system()
    if system is None or not isinstance(error, system.Exception):
        return None
    inner = error
    while True:
        if isinstance(inner, system.Reflection.TargetInvocationException) and inner.InnerException is not None:
            inner = inner.InnerException
        elif isinstance(inner, system.AggregateException) and inner.InnerExceptions.Count == 1:
            inner = inner.InnerExceptions[0]
        else:
            break
    message = str(inner.Message)
    if isinstance(inner, system.ArgumentException):
        return ValueError(message)
    if isinstance(inner, system.OverflowException):
        return OverflowError(message)
    if isinstance(inner, system.DivideByZeroException):
        return ZeroDivisionError(message)
    if isinstance(inner, system.ArithmeticException):
        return ArithmeticError(message)
    if isinstance(inner, (system.NotImplementedException, system.NotSupportedException)):
        return NotImplementedError(message)
    return RuntimeError(message)


@contextmanager
def net_errors() -> Iterator[None]:
    """Re-raise .NET exceptions raised in the block as their Python counterparts.

    Yields
    ------
    None
        Nothing; the body runs inside the translation.

    Examples
    --------
    >>> from pynumerics.distributions import Normal
    >>> from pynumerics.errors import net_errors
    >>> try:
    ...     with net_errors():
    ...         Normal().net.InverseCDF(1.5)
    ... except ValueError as error:
    ...     print("refused")
    refused
    """
    try:
        yield
    except Exception as error:
        mapped = translate(error)
        if mapped is None:
            raise
        raise mapped from error


def translated(func: F) -> F:
    """Return ``func`` wrapped so that .NET exceptions it raises surface as Python exceptions.

    The wrapper keeps the name, docstring and signature of ``func`` (``functools.wraps``) and adds
    one Python call to a call that does not raise. Wrapping an already translated function
    returns it unchanged.

    Parameters
    ----------
    func : callable
        The function to wrap.

    Returns
    -------
    callable
        The wrapped function.
    """
    if getattr(func, _MARKER, False):
        return func

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        try:
            return func(*args, **kwargs)
        except Exception as error:
            mapped = translate(error)
            if mapped is None:
                raise
            raise mapped from error

    setattr(wrapper, _MARKER, True)
    return wrapper  # type: ignore[return-value]


def is_translated(obj: Any) -> bool:
    """Return whether a function, property or class member translates .NET exceptions.

    Parameters
    ----------
    obj : object
        A function, ``property``, ``classmethod`` or ``staticmethod``.

    Returns
    -------
    bool
        True when every function the object runs is wrapped by :func:`translated`.
    """
    if isinstance(obj, property):
        return all(f is None or getattr(f, _MARKER, False) for f in (obj.fget, obj.fset, obj.fdel))
    if isinstance(obj, (classmethod, staticmethod)):
        return bool(getattr(obj.__func__, _MARKER, False))
    return bool(getattr(obj, _MARKER, False))


def _translate_member(value: Any) -> Any:
    """Return a class-dictionary entry with its functions wrapped by :func:`translated`."""
    if isinstance(value, property):
        if is_translated(value):
            return value
        return property(*(None if f is None else translated(f) for f in (value.fget, value.fset, value.fdel)),
                        doc=value.__doc__)
    if isinstance(value, (classmethod, staticmethod)):
        return value if is_translated(value) else type(value)(translated(value.__func__))
    if inspect.isfunction(value):
        return translated(value)
    return value


def translate_class(klass: type) -> type:
    """Wrap the methods and properties a class defines so .NET exceptions surface as Python exceptions.

    Covers the functions, properties, class methods and static methods in the class's own
    dictionary whose names are public, plus ``__init__``, ``__call__``, ``__getitem__``,
    ``__setitem__``, ``__len__``, ``__contains__``, ``__iter__``, ``__eq__`` and ``__hash__``.
    Enums and exception classes are returned unchanged.

    Parameters
    ----------
    klass : type
        The class to update in place.

    Returns
    -------
    type
        ``klass``.
    """
    if issubclass(klass, (enum.Enum, BaseException)):
        return klass
    for name, value in list(vars(klass).items()):
        if (name.startswith("_") and name not in _WRAPPED_DUNDERS) or value is None:
            continue
        member = _translate_member(value)
        if member is not value:
            setattr(klass, name, member)
    return klass


def translate_public(namespace: MutableMapping[str, Any]) -> None:
    """Apply :func:`translated` and :func:`translate_class` to the public objects a module defines.

    Called at the end of a wrapper module as ``translate_public(globals())``. The objects are
    the module's ``__all__`` (or every public name when it has none) whose ``__module__`` is the
    module itself; re-exported objects are left to their own module.

    Parameters
    ----------
    namespace : MutableMapping
        The module's ``globals()``.
    """
    module = namespace.get("__name__")
    names = namespace.get("__all__") or [n for n in list(namespace) if not n.startswith("_")]
    for name in names:
        obj = namespace.get(name)
        if getattr(obj, "__module__", None) != module:
            continue
        if inspect.isclass(obj):
            translate_class(obj)
        elif inspect.isfunction(obj):
            namespace[name] = translated(obj)


__all__ = ["is_translated", "net_errors", "translate", "translate_class", "translate_public", "translated"]
