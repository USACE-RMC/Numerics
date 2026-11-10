"""Python callables and C# source as the .NET delegates Numerics members take.

Numerics takes functions as delegates: ``Func<double, double>`` for integrands and root-finding
functions, ``Func<double[], double>`` for objectives, ``LogLikelihood`` for the MCMC samplers and
others. Every wrapper that takes one accepts any of:

- a Python callable, wrapped by :func:`to_delegate`: arrays reach it as NumPy float64 copies and
  its return value is converted to the delegate's return type. Each evaluation crosses from
  .NET into Python and holds the GIL. Measured on a desktop, the crossing of a raw pythonnet
  delegate costs about 2 microseconds; the wrapper adds about 0.2 microseconds for scalar
  arguments, about 0.4 microseconds for an array of up to :data:`BUFFERED_ARRAY_LENGTH` values,
  and one copy at memory speed beyond that (about 30 microseconds for 100,000 values);
- a string of C# source, compiled by :func:`compile_csharp` into a delegate that runs entirely in
  .NET, with no Python crossing per evaluation;
- a .NET delegate of the right type, passed through unchanged.

Parallel rule: a Python callback serializes on the GIL, so a Numerics object that would call it
from several threads runs serially instead. :func:`uses_python_callback` tells which case applies,
and :func:`apply_parallel_policy` turns the object's parallel switch off when any of its delegates
calls Python; compiled and .NET delegates leave it as it was.

:func:`compile_csharp` uses the Roslyn C# compiler (``Microsoft.CodeAnalysis.CSharp``) bundled in
the wheel's ``_bin`` directory, pinned by version and SHA-256 in ``pynumerics._version``. The
first compilation in a process takes one to two seconds while Roslyn starts; later ones take tens
of milliseconds, and a source compiled once is reused from a cache.
"""

from __future__ import annotations

import hashlib
import threading
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import numpy as np

from pynumerics import _version, runtime
from pynumerics._util import namespace
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    to_net_doubles,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_ints,
    to_numpy_matrix,
)

ROSLYN_ASSEMBLIES: tuple[tuple[str, str], ...] = (
    ("Microsoft.CodeAnalysis.dll", _version.ROSLYN_CODEANALYSIS_SHA256),
    ("Microsoft.CodeAnalysis.CSharp.dll", _version.ROSLYN_CSHARP_SHA256),
)
"""The bundled Roslyn assemblies and their pinned SHA-256 digests, in load order."""

DEFAULT_USINGS: tuple[str, ...] = (
    "System", "System.Collections.Generic", "System.Linq", "Numerics", "Numerics.Distributions",
    "Numerics.Mathematics", "Numerics.Mathematics.LinearAlgebra", "Numerics.Mathematics.SpecialFunctions",
)
"""The namespaces C# callback source can use without qualification."""

PARALLEL_SWITCHES: tuple[str, ...] = ("ParallelizeChains", "UseParallel", "IsParallel", "Parallel")
"""The Numerics boolean properties that run work on several threads, in the order they are looked for."""

_LOCK = threading.RLock()
_DELEGATE_TYPES: dict[str, Any] = {}
_COMPILED: dict[tuple[str, str, tuple[str, ...]], Any] = {}
_REFERENCES: dict[str, Any] = {}
_ROSLYN: dict[str, Any] = {}
_SIGNATURES: dict[str, tuple[tuple[str, ...], str]] = {}
_STAGING: dict[str, tuple[Any, Any] | None] = {}

_RAW_DELEGATES = frozenset({"Numerics.Mathematics.Optimization.TransitionDistribution"})
"""Delegates whose array arguments are outputs the callback fills; they reach Python as the .NET arrays."""

_KEYWORDS = {
    "System.Double": "double", "System.Single": "float", "System.Int32": "int", "System.Int64": "long",
    "System.Boolean": "bool", "System.String": "string", "System.Object": "object", "System.Void": "void",
}

_IN_DOUBLE_LIST = "System.Collections.Generic.IList`1[[System.Double"


def _delegate_table() -> dict[str, Any]:
    """Return the named delegate types, keyed by their C# spelling without spaces; built once."""
    if not _DELEGATE_TYPES:
        namespace("Numerics")
        import System  # noqa: PLC0415 - the runtime is running
        from System.Collections.Generic import IList  # noqa: PLC0415

        linear = namespace("Numerics.Mathematics.LinearAlgebra")
        optimization = namespace("Numerics.Mathematics.Optimization")
        mcmc = namespace("Numerics.Sampling.MCMC")
        d, darr, i, b = System.Double, System.Array[System.Double], System.Int32, System.Boolean
        dlist = IList[System.Double]
        func, action = System.Func, System.Action
        _DELEGATE_TYPES.update({
            "Func<double,double>": func[d, d],
            "Func<double,double,double>": func[d, d, d],
            "Func<double,double,double,double>": func[d, d, d, d],
            "Func<double[],double>": func[darr, d],
            "Func<double[],double[]>": func[darr, darr],
            "Func<double[],double,double>": func[darr, d, d],
            "Func<double,double[],double>": func[d, darr, d],
            "Func<double,double[],double[]>": func[d, darr, darr],
            "Func<double[],double[],double>": func[darr, darr, d],
            "Func<double[],bool>": func[darr, b],
            "Func<IList<double>,double>": func[dlist, d],
            "Func<IList<double>,double[]>": func[dlist, darr],
            "Func<int,int,int,double>": func[i, i, i, d],
            "Func<int,int,int,int>": func[i, i, i, i],
            "Func<Vector,Vector>": func[linear.Vector, linear.Vector],
            "Func<Vector,Matrix>": func[linear.Vector, linear.Matrix],
            "Func<ParameterSet,double[]>": func[optimization.ParameterSet, darr],
            "Action<double,double,double>": action[d, d, d],
            "Action<double,double,double,double>": action[d, d, d, d],
            "LogLikelihood": mcmc.LogLikelihood,
            "HMC.Gradient": mcmc.HMC.Gradient,
            "Gibbs.Proposal": mcmc.Gibbs.Proposal,
            "MCMCSampler.ProgressChangedEventHandler": mcmc.MCMCSampler.ProgressChangedEventHandler,
            "NUTS.LogLikelihoodAndGradient": mcmc.NUTS.LogLikelihoodAndGradient,
            "TransitionDistribution": optimization.TransitionDistribution,
        })
    return _DELEGATE_TYPES


def delegate_types() -> tuple[str, ...]:
    """Return the delegate type names :func:`delegate_type` resolves.

    Returns
    -------
    tuple of str
        C# spellings such as ``"Func<double,double>"`` and Numerics delegate names such as
        ``"LogLikelihood"`` and ``"HMC.Gradient"``.
    """
    return tuple(_delegate_table())


def delegate_type(spec: Any) -> Any:
    """Resolve a delegate type from its C# spelling, or return a .NET delegate type unchanged.

    Parameters
    ----------
    spec : str or type
        A name from :func:`delegate_types` (spaces are ignored, so ``"Func<double, double>"``
        works), or a pythonnet delegate type such as ``System.Func[System.Double, System.Double]``.

    Returns
    -------
    type
        The pythonnet delegate type.

    Raises
    ------
    ValueError
        If ``spec`` names no known delegate type.
    TypeError
        If ``spec`` is a type that is not a .NET delegate type.
    """
    table = _delegate_table()
    if isinstance(spec, str):
        key = "".join(spec.split())
        if key not in table:
            raise ValueError(f"Unknown delegate type {spec!r}; known: {', '.join(table)}.")
        return table[key]
    import clr  # noqa: PLC0415 - the runtime is running
    import System  # noqa: PLC0415

    try:
        clr_type = clr.GetClrType(spec)
    except Exception as error:
        raise TypeError(f"{spec!r} is not a .NET delegate type.") from error
    if not clr_type.IsSubclassOf(clr.GetClrType(System.MulticastDelegate)):
        raise TypeError(f"{clr_type.FullName} is not a .NET delegate type.")
    return spec


_resolve_delegate_type = delegate_type  # compile_csharp's parameter of the same name shadows it


def _signature(dtype: Any) -> tuple[tuple[str, ...], str]:
    """Return the full names of a delegate's ``Invoke`` parameter types and return type, cached."""
    import clr  # noqa: PLC0415 - the runtime is running

    clr_type = clr.GetClrType(dtype)
    key = str(clr_type.AssemblyQualifiedName)
    signature = _SIGNATURES.get(key)
    if signature is None:
        invoke = clr_type.GetMethod("Invoke")
        signature = _SIGNATURES[key] = (tuple(str(p.ParameterType.FullName) for p in invoke.GetParameters()),
                                        str(invoke.ReturnType.FullName))
    return signature


def _argument(type_name: str) -> Callable[[Any], Any] | None:
    """Return the converter of a .NET argument for a Python callback, or None to pass it as it is."""
    if type_name == "System.Double[]" or type_name.startswith(_IN_DOUBLE_LIST):
        return to_numpy_doubles
    if type_name == "System.Int32[]":
        return to_numpy_ints
    if type_name == "System.Double[,]":
        return to_numpy_matrix
    if type_name == "Numerics.Mathematics.LinearAlgebra.Vector":
        return lambda vector: to_numpy_doubles(vector.Array)
    if type_name == "Numerics.Mathematics.LinearAlgebra.Matrix":
        return lambda matrix: to_numpy_matrix(matrix.ToArray())
    return None


def _result(type_name: str) -> Callable[[Any], Any] | None:
    """Return the converter of a Python callback's result to the .NET return type, or None to pass it."""
    import System  # noqa: PLC0415 - the runtime is running

    def net(convert: Callable[[Any], Any]) -> Callable[[Any], Any]:
        return lambda value: value if isinstance(value, System.Object) else convert(value)

    if type_name == "System.Double":
        return float
    if type_name in ("System.Int32", "System.Int64"):
        return int
    if type_name == "System.Boolean":
        return bool
    if type_name == "System.Void":
        return lambda value: None
    if type_name == "System.Double[]":
        return net(to_net_doubles)
    if type_name == "System.Double[,]":
        return net(to_net_matrix)
    linear = namespace("Numerics.Mathematics.LinearAlgebra")
    if type_name == "Numerics.Mathematics.LinearAlgebra.Vector":
        return net(lambda value: linear.Vector(to_net_doubles(value)))
    if type_name == "Numerics.Mathematics.LinearAlgebra.Matrix":
        return net(lambda value: linear.Matrix(to_net_matrix(value)))
    if type_name.startswith("System.ValueTuple`2[[System.Double") and "LinearAlgebra.Vector" in type_name:
        pair = System.ValueTuple[System.Double, linear.Vector]
        vector = net(lambda value: linear.Vector(to_net_doubles(value)))
        return net(lambda value: pair(float(value[0]), vector(value[1])))
    return None


def to_delegate(fn: Any, delegate_type_spec: Any) -> Any:
    """Return ``fn`` as a .NET delegate of the given type.

    A Python callable is wrapped: array, ``IList<double>``, ``Vector`` and ``Matrix`` arguments
    reach it as new NumPy arrays (``TransitionDistribution`` receives its output arrays as the .NET
    arrays it fills), scalars as Python numbers and anything else as the .NET object; the result is
    converted to the return type (a float for ``double``, a ``System.Double[]`` for ``double[]``, a
    ``Vector`` for ``Vector``, and so on). A delegate whose arguments and result are all ``double``
    calls ``float(fn(...))`` with no other conversion. A delegate whose one argument is a
    ``double[]``, ``IList<double>`` or ``Vector`` of at most :data:`BUFFERED_ARRAY_LENGTH` values
    stages it in a per-thread buffer, so an evaluation costs about 0.4 microseconds more than a
    raw pythonnet delegate; longer arrays cross in one copy. The values are bit for bit the ones
    .NET passed either way. A string is compiled by :func:`compile_csharp`. A .NET delegate of the
    type is returned unchanged.

    Parameters
    ----------
    fn : callable, str or System.Delegate
        The function.
    delegate_type_spec : str or type
        The delegate type, as for :func:`delegate_type`.

    Returns
    -------
    System.Delegate
        A delegate of the requested type.

    Raises
    ------
    TypeError
        If ``fn`` is neither callable, a string nor a delegate of the type.
    ValueError
        If the delegate type is unknown, or C# source does not compile.
    """
    dtype = delegate_type(delegate_type_spec)
    if isinstance(fn, str):
        return compile_csharp(fn, dtype)
    import clr  # noqa: PLC0415 - the runtime is running
    import System  # noqa: PLC0415

    if isinstance(fn, System.Delegate):
        if not clr.GetClrType(dtype).IsInstanceOfType(fn):
            raise TypeError(f"{fn.GetType().FullName} is not a {clr.GetClrType(dtype).FullName}.")
        return fn
    if not callable(fn):
        raise TypeError(f"Expected a callable, C# source or a .NET delegate; got {type(fn).__name__}.")
    arguments, result = _signature(dtype)
    if result == "System.Double" and all(a == "System.Double" for a in arguments):
        return dtype(_double_wrapper(fn, len(arguments)))
    raw = str(clr.GetClrType(dtype).FullName) in _RAW_DELEGATES
    converters = [None if raw else _argument(a) for a in arguments]
    back = _result(result)
    wrapper = dtype(_converting_wrapper(fn, converters, back))
    if not raw and len(arguments) == 1 and result != "System.Void":
        return _buffered_delegate(fn, dtype, arguments[0], back, wrapper)
    return wrapper


BUFFERED_ARRAY_LENGTH = 1024
"""Longest array a one-argument Python callback receives through its per-thread staging buffer.

Below this length the .NET side copies the argument into a staging buffer that NumPy already
views, and Python copies it out with one slice, which avoids acquiring a buffer view of a new .NET
array on every evaluation. Longer arrays use the direct one-copy path, which is faster there.
"""


def _buffered_delegate(fn: Callable[..., Any], dtype: Any, argument: str, back: Callable[[Any], Any] | None,
                       direct: Any) -> Any:
    """Return a delegate that hands ``fn`` a ``double[]``, ``IList<double>`` or ``Vector`` argument cheaply.

    The delegate, built by the factory of :func:`_staging_factory`, copies an argument of at most
    :data:`BUFFERED_ARRAY_LENGTH` values into a staging ``double[]`` owned by the calling thread,
    then calls Python with only the length; Python returns ``fn`` of a copy of the first ``length``
    staged values, bit for bit the values .NET passed. A null or longer argument goes to ``direct``,
    the plain converting delegate. Each thread has its own staging buffer, registered with Python
    the first time that thread calls, so concurrent evaluations never share one; a nested call on
    the same thread is safe because Python copies the staged values before calling ``fn``. Other
    argument types return ``direct`` unchanged.
    """
    factory = _staging_factory(dtype, argument)
    if factory is None:
        return direct
    import System  # noqa: PLC0415 - the runtime is running
    from System.Threading import ThreadLocal  # noqa: PLC0415

    views: dict[int, np.ndarray] = {}
    ident = threading.get_ident

    def register(buffer: Any) -> None:
        views[ident()] = np.frombuffer(buffer, dtype=np.float64)

    if back is None:
        def staged(length: int) -> Any:
            return fn(views[ident()][:length].copy())
    else:
        def staged(length: int) -> Any:
            return back(fn(views[ident()][:length].copy()))

    doubles = System.Array[System.Double]
    python, _ = factory
    return factory[1](direct, python(staged), System.Action[doubles](register), ThreadLocal[doubles]())


def _staging_factory(dtype: Any, argument: str) -> tuple[Any, Any] | None:
    """Return ``(python delegate type, factory)`` for the staging delegates of ``dtype``, compiled once.

    The factory is a compiled ``System.Linq.Expressions`` lambda taking the direct delegate, the
    Python ``Func<int, R>`` that reads the staged values, the ``Action<double[]>`` that registers a
    new staging buffer, and the ``ThreadLocal<double[]>`` holding each thread's buffer; it returns
    the staging delegate of ``dtype`` closed over them. Compiling once per delegate type keeps each
    :func:`to_delegate` call in the microseconds. Returns None for an argument type it cannot stage.
    """
    import clr  # noqa: PLC0415 - the runtime is running
    import System  # noqa: PLC0415

    clr_type = clr.GetClrType(dtype)
    key = str(clr_type.AssemblyQualifiedName)
    with _LOCK:
        if key in _STAGING:
            return _STAGING[key]
        if argument == "System.Double[]":
            source = "array"
        elif argument.startswith(_IN_DOUBLE_LIST):
            source = "list"
        elif argument == "Numerics.Mathematics.LinearAlgebra.Vector":
            source = "vector"
        else:
            _STAGING[key] = None
            return None
        from System.Collections.Generic import ICollection, List  # noqa: PLC0415
        from System.Linq.Expressions import Expression, ParameterExpression  # noqa: PLC0415
        from System.Threading import ThreadLocal  # noqa: PLC0415

        doubles = System.Array[System.Double]
        double_array, integer = clr.GetClrType(doubles), clr.GetClrType(System.Int32)
        invoke = clr_type.GetMethod("Invoke")
        python_type = System.Func[System.Int32, invoke.ReturnType]
        register_type, local_type = System.Action[doubles], ThreadLocal[doubles]
        direct = Expression.Parameter(clr_type, "direct")
        python = Expression.Parameter(clr.GetClrType(python_type), "python")
        register = Expression.Parameter(clr.GetClrType(register_type), "register")
        local = Expression.Parameter(clr.GetClrType(local_type), "local")
        capacity = Expression.Constant(System.Int32(BUFFERED_ARRAY_LENGTH))
        x = Expression.Parameter(invoke.GetParameters()[0].ParameterType, "x")
        buffer = Expression.Variable(double_array, "buffer")
        length = Expression.Variable(integer, "length")
        array = Expression.Variable(double_array, "array")
        no_array = Expression.Constant(None, double_array)
        if source == "list":
            collection = Expression.Convert(x, clr.GetClrType(ICollection[System.Double]))
            measure = Expression.Assign(length, Expression.Property(collection, "Count"))
            copy = Expression.Call(collection, "CopyTo", None, buffer, Expression.Constant(System.Int32(0)))
        else:
            values = x if source == "array" else array
            measure = Expression.Assign(length, Expression.ArrayLength(values))
            copy = Expression.Call(clr.GetClrType(System.Array), "Copy", None, values, buffer, length)
        stored = Expression.Property(local, "Value")
        fresh = Expression.Block(
            Expression.Assign(buffer, Expression.NewArrayBounds(clr.GetClrType(System.Double), capacity)),
            Expression.Assign(stored, buffer),
            Expression.Invoke(register, buffer))
        fast = Expression.Block(
            Expression.Assign(buffer, stored),
            Expression.IfThen(Expression.Equal(buffer, no_array), fresh),
            copy,
            Expression.Invoke(python, length))
        slow = Expression.Invoke(direct, x)
        steps = []
        missing = Expression.Equal(x, Expression.Constant(None, x.Type))
        if source == "vector":
            vector_array = Expression.Condition(missing, no_array, Expression.Property(x, "Array"))
            steps.append(Expression.Assign(array, vector_array))
            missing = Expression.Equal(array, no_array)
        steps.append(Expression.Condition(
            missing, slow,
            Expression.Block(measure, Expression.Condition(Expression.GreaterThan(length, capacity), slow, fast))))
        variables = List[ParameterExpression]()
        for variable in (buffer, length, array):
            variables.Add(variable)
        body = Expression.Block(invoke.ReturnType, variables, System.Array[Expression](steps))
        inner = Expression.Lambda(clr_type, body, System.Array[ParameterExpression]([x]))
        factory_type = System.Func[dtype, python_type, register_type, local_type, dtype]
        outer = Expression.Lambda(clr.GetClrType(factory_type), inner,
                                  System.Array[ParameterExpression]([direct, python, register, local]))
        _STAGING[key] = (python_type, outer.Compile())
        return _STAGING[key]


def _converting_wrapper(fn: Callable[..., Any], converters: list[Callable[[Any], Any] | None],
                        back: Callable[[Any], Any] | None) -> Callable[..., Any]:
    """Return a wrapper that converts the arguments, calls ``fn`` and converts its result.

    A one-argument delegate (objectives, log-likelihoods, gradients) gets a closure with its
    converters bound directly, so an evaluation makes no per-argument loop or tuple.
    """
    if len(converters) == 1:
        (convert,) = converters
        if convert is None:
            return (lambda x: fn(x)) if back is None else (lambda x: back(fn(x)))
        return (lambda x: fn(convert(x))) if back is None else (lambda x: back(fn(convert(x))))

    def wrapper(*values: Any) -> Any:
        values = tuple(v if c is None else c(v) for c, v in zip(converters, values, strict=True))
        value = fn(*values)
        return value if back is None else back(value)

    return wrapper


def _double_wrapper(fn: Callable[..., Any], arity: int) -> Callable[..., float]:
    """Return a wrapper of an all-``double`` callback that only converts the result with ``float``."""
    if arity == 1:
        return lambda x: float(fn(float(x)))
    if arity == 2:
        return lambda x, y: float(fn(x, y))
    return lambda *values: float(fn(*values))


def scalar_function(fn: Any) -> Any:
    """Return ``fn`` as a ``Func<double, double>`` (integrands, root-finding and 1-D optimization).

    Parameters
    ----------
    fn : callable, str or System.Delegate
        A function of one float, C# source or a delegate, as for :func:`to_delegate`.

    Returns
    -------
    System.Func[Double, Double]
        The delegate.
    """
    return to_delegate(fn, "Func<double,double>")


def vector_function(fn: Any) -> Any:
    """Return ``fn`` as a ``Func<double[], double>`` (objective functions and multidimensional integrands).

    Parameters
    ----------
    fn : callable, str or System.Delegate
        A function of a float64 array returning a float, C# source or a delegate.

    Returns
    -------
    System.Func[Double[], Double]
        The delegate.
    """
    return to_delegate(fn, "Func<double[],double>")


def bivariate_function(fn: Any) -> Any:
    """Return ``fn`` as a ``Func<double, double, double>`` (2-D integrands and ODE right-hand sides).

    Parameters
    ----------
    fn : callable, str or System.Delegate
        A function of two floats, C# source or a delegate.

    Returns
    -------
    System.Func[Double, Double, Double]
        The delegate.
    """
    return to_delegate(fn, "Func<double,double,double>")


def log_likelihood(fn: Any) -> Any:
    """Return ``fn`` as a ``Numerics.Sampling.MCMC.LogLikelihood`` (the MCMC target density).

    Parameters
    ----------
    fn : callable, str or System.Delegate
        A function of the parameter array returning the log-likelihood, C# source or a delegate.

    Returns
    -------
    Numerics.Sampling.MCMC.LogLikelihood
        The delegate.
    """
    return to_delegate(fn, "LogLikelihood")


def _is_python_delegate(value: Any) -> bool:
    """Return whether a .NET delegate calls into Python (pythonnet builds those in a generated assembly).

    A compiled expression lambda, such as the staging delegate of :func:`to_delegate`, calls Python
    when any delegate among its closure constants or captured variables does.
    """
    import System  # noqa: PLC0415 - the runtime is running

    for item in value.GetInvocationList():
        target = item.Target
        if target is None:
            continue
        target_type = target.GetType()
        if str(target_type.Assembly.GetName().Name).startswith("__Python_Runtime"):
            return True
        if str(target_type.FullName) == "System.Runtime.CompilerServices.Closure":
            captured = [*(target.Constants or ()), *(getattr(box, "Value", None) for box in (target.Locals or ()))]
            if any(isinstance(c, System.Delegate) and _is_python_delegate(c) for c in captured):
                return True
    return False


def uses_python_callback(*functions: Any) -> bool:
    """Return whether any of the functions calls into Python when Numerics invokes it.

    Parameters
    ----------
    *functions : callable, str, System.Delegate or None
        The functions a wrapper passes to Numerics. A Python callable, or a .NET delegate built
        from one, calls Python; C# source and compiled or native .NET delegates do not. None is
        ignored.

    Returns
    -------
    bool
        True when at least one function calls Python, so the Numerics object must not run it
        from several threads.
    """
    for fn in functions:
        if fn is None or isinstance(fn, str):
            continue
        if runtime.is_initialized():
            import System  # noqa: PLC0415 - the runtime is running

            if isinstance(fn, System.Delegate):
                if _is_python_delegate(fn):
                    return True
                continue
        if callable(fn):
            return True
    return False


def apply_parallel_policy(net: Any, *functions: Any) -> bool:
    """Turn off a Numerics object's parallel switch when any of its functions calls Python.

    The first property of :data:`PARALLEL_SWITCHES` the object has is set to False when
    :func:`uses_python_callback` is True, and left as it is otherwise.

    Parameters
    ----------
    net : object
        The .NET object (for example an MCMC sampler or ``DifferentialEvolution``).
    *functions : callable, str, System.Delegate or None
        The functions given to the object.

    Returns
    -------
    bool
        True when the switch was turned off.
    """
    if not uses_python_callback(*functions):
        return False
    for name in PARALLEL_SWITCHES:
        if hasattr(net, name):
            setattr(net, name, False)
            return True
    return False


def _roslyn_directory() -> Path:
    """Return the directory holding the Roslyn assemblies: the assembly directory, else the packaged ``_bin``."""
    first = ROSLYN_ASSEMBLIES[0][0]
    for directory in (runtime.runtime_dir(), Path(runtime.__file__).resolve().parent / "_bin"):
        if (directory / first).is_file():
            return directory
    raise runtime.RuntimeLoadError(
        f"{first} is not beside Numerics.dll or in the packaged _bin directory; install the pynumerics wheel or "
        "stage it with python/pynumerics/scripts/build_wheel.py --roslyn-only.")


def roslyn_available() -> bool:
    """Return whether the bundled Roslyn assemblies are present, so :func:`compile_csharp` can run.

    Returns
    -------
    bool
        True when both assemblies are found (their digests are checked when they load).
    """
    try:
        directory = _roslyn_directory()
    except runtime.RuntimeLoadError:
        return False
    return all((directory / name).is_file() for name, _ in ROSLYN_ASSEMBLIES)


def _load_roslyn() -> dict[str, Any]:
    """Verify and load the Roslyn assemblies once, and return the Roslyn types the compiler uses."""
    if _ROSLYN:
        return _ROSLYN
    directory = _roslyn_directory()
    for name, pinned in ROSLYN_ASSEMBLIES:
        file = directory / name
        digest = hashlib.sha256(file.read_bytes()).hexdigest()
        if digest != pinned:
            raise runtime.RuntimeLoadError(f"{file} has SHA-256 {digest}; pynumerics pins {pinned}.")
        runtime.add_assembly(file)
    code_analysis = namespace("Microsoft.CodeAnalysis")
    csharp = namespace("Microsoft.CodeAnalysis.CSharp")
    _ROSLYN.update({"ca": code_analysis, "cs": csharp})
    return _ROSLYN


def _references() -> Any:
    """Return the compilation references: the .NET shared framework, Numerics and every added assembly."""
    import System  # noqa: PLC0415 - the runtime is running
    from System.Collections.Generic import List  # noqa: PLC0415

    roslyn = _load_roslyn()
    added = tuple(sorted(str(p) for p in runtime.loaded_assemblies().values()))
    key = "|".join(added)
    references = _REFERENCES.get(key)
    if references is None:
        metadata = roslyn["ca"].MetadataReference
        references = List[metadata]()
        platform = str(System.AppContext.GetData("TRUSTED_PLATFORM_ASSEMBLIES")).split(System.IO.Path.PathSeparator)
        for path in [*platform, str(runtime.runtime_dir() / runtime.ASSEMBLY_FILE), *added]:
            if path:
                references.Add(metadata.CreateFromFile(path))
        _REFERENCES.clear()
        _REFERENCES[key] = references
    return references


def _csharp_name(clr_type: Any) -> str:
    """Return the fully qualified C# spelling of a .NET type (``global::System.Func<double, double>``)."""
    if clr_type.IsArray:
        return f"{_csharp_name(clr_type.GetElementType())}[{',' * (clr_type.GetArrayRank() - 1)}]"
    full = str(clr_type.FullName) if clr_type.FullName is not None else str(clr_type.Name)
    if full in _KEYWORDS:
        return _KEYWORDS[full]
    name = str(clr_type.Name).split("`")[0]
    if clr_type.IsGenericType:
        arguments = [a for a in clr_type.GetGenericArguments()]
        if clr_type.IsNested:
            arguments = arguments[len(clr_type.DeclaringType.GetGenericArguments()):]
        name += "<" + ", ".join(_csharp_name(a) for a in arguments) + ">"
    if clr_type.IsNested:
        return f"{_csharp_name(clr_type.DeclaringType)}.{name}"
    return f"global::{clr_type.Namespace}.{name}" if clr_type.Namespace else f"global::{name}"


_CLASS = "PynumericsCallback"


def _emit(code: str, assembly_name: str) -> Any:
    """Compile C# code into an in-memory assembly and load it; raise ValueError with the errors."""
    import System  # noqa: PLC0415 - the runtime is running

    roslyn = _load_roslyn()
    ca, cs = roslyn["ca"], roslyn["cs"]
    tree = cs.CSharpSyntaxTree.ParseText(code)
    options = cs.CSharpCompilationOptions(ca.OutputKind.DynamicallyLinkedLibrary).WithOptimizationLevel(
        ca.OptimizationLevel.Release)
    compilation = cs.CSharpCompilation.Create(assembly_name, System.Array[ca.SyntaxTree]([tree]), _references(),
                                              options)
    stream = System.IO.MemoryStream()
    emitted = compilation.Emit(stream)
    if not emitted.Success:
        messages = [str(d) for d in emitted.Diagnostics if d.Severity == ca.DiagnosticSeverity.Error]
        raise ValueError("C# callback source does not compile:\n" + "\n".join(messages))
    return System.Reflection.Assembly.Load(stream.ToArray())


def compile_csharp(source: str, delegate_type: Any = "Func<double,double>", *,
                   usings: Iterable[str] = ()) -> Any:
    """Compile C# source into a .NET delegate that runs without calling back into Python.

    ``source`` is either an expression convertible to the delegate type, usually a lambda
    (``"x => Math.Exp(-x * x)"``), or one or more static method declarations of which exactly one
    matches the delegate's signature (``"static double F(double x) { return x * x; }"``). The code
    is compiled with Roslyn against the .NET shared framework, ``Numerics.dll`` and every assembly
    added with :func:`pynumerics.add_assembly`, so it can use Numerics types; the namespaces in
    :data:`DEFAULT_USINGS` are imported. Compiled delegates are cached by source, delegate type and
    usings, so compiling the same source again returns the same delegate.

    Parameters
    ----------
    source : str
        The C# lambda, expression or static methods.
    delegate_type : str or type, optional
        The delegate type, as for :func:`delegate_type`; default ``Func<double,double>``.
    usings : iterable of str, optional
        More namespaces to import, for example ``("Numerics.Data.Statistics",)``.

    Returns
    -------
    System.Delegate
        The compiled delegate.

    Raises
    ------
    ValueError
        If the source does not compile (the message lists Roslyn's errors with line and column
        numbers counted in ``source``), if no or several static methods match the delegate, or if
        the delegate type is unknown.
    pynumerics.RuntimeLoadError
        If the Roslyn assemblies are missing or do not match their pinned SHA-256.

    Examples
    --------
    >>> from pynumerics.callbacks import compile_csharp
    >>> f = compile_csharp("x => Math.Exp(-x * x)")  # doctest: +SKIP
    >>> f(0.0)  # doctest: +SKIP
    1.0
    """
    dtype = _resolve_delegate_type(delegate_type)
    import clr  # noqa: PLC0415 - the runtime is running
    import System  # noqa: PLC0415

    clr_type = clr.GetClrType(dtype)
    extra = tuple(usings)
    key = (source, str(clr_type.AssemblyQualifiedName), extra)
    with _LOCK:
        compiled = _COMPILED.get(key)
        if compiled is not None:
            return compiled
        roslyn = _load_roslyn()
        header = "".join(f"using {name};\n" for name in (*DEFAULT_USINGS, *extra))
        assembly_name = f"{_CLASS}{len(_COMPILED) + 1}"
        type_name = _csharp_name(clr_type)
        expression = roslyn["cs"].SyntaxFactory.ParseExpression(source)
        if not expression.ContainsDiagnostics:
            code = (f"{header}public static class {_CLASS}\n{{\n    public static {type_name} Create()\n    {{\n"
                    f"        return\n#line 1\n{source}\n#line default\n        ;\n    }}\n}}\n")
            assembly = _emit(code, assembly_name)
            compiled = assembly.GetType(_CLASS).GetMethod("Create").Invoke(None, None)
        else:
            code = f"{header}public static class {_CLASS}\n{{\n#line 1\n{source}\n#line default\n}}\n"
            assembly = _emit(code, assembly_name)
            flags = (System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.NonPublic
                     | System.Reflection.BindingFlags.Static | System.Reflection.BindingFlags.DeclaredOnly)
            matches = [m for m in assembly.GetType(_CLASS).GetMethods(flags)
                       if System.Delegate.CreateDelegate(clr_type, m, False) is not None]
            if len(matches) != 1:
                raise ValueError(f"Expected exactly one static method matching {type_name}; found {len(matches)}.")
            compiled = System.Delegate.CreateDelegate(clr_type, matches[0])
        _COMPILED[key] = compiled
        return compiled


__all__ = [
    "BUFFERED_ARRAY_LENGTH", "DEFAULT_USINGS", "PARALLEL_SWITCHES", "ROSLYN_ASSEMBLIES", "apply_parallel_policy",
    "bivariate_function",
    "compile_csharp", "delegate_type", "delegate_types", "log_likelihood", "roslyn_available", "scalar_function",
    "to_delegate", "uses_python_callback", "vector_function",
]

translate_public(globals())
