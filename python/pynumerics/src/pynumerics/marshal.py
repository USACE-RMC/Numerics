"""Copy arrays between NumPy and .NET in bulk, preserving every bit.

Each non-empty transfer is one bulk copy through the Python buffer protocol, never
element-by-element indexing, so values keep their exact bit patterns (NaN payloads, signed
zeros, infinities and subnormals included) and large arrays copy at memory speed.

pythonnet exports a .NET array of a primitive type (``Double``, ``Single``, ``Int32``,
``Int64``, ``Byte``, ``SByte``, of any rank) through the buffer protocol. Acquiring the buffer
pins the array with a ``GCHandle`` and releasing it frees the handle; every copy here holds the
buffer in a ``with memoryview(...)`` block, so the handle is freed when the block ends, also
when the copy raises. A write allocates the .NET array and assigns the NumPy values into a
``numpy.frombuffer`` view of its storage; a read copies such a view into a new NumPy array.
Both are one ``memcpy`` for contiguous data, the same bits as
``System.Runtime.InteropServices.Marshal.Copy`` at less than half of its call cost. Callbacks
rely on this, since every evaluation copies its array arguments.

Two-dimensional arrays (``double[,]``, which Numerics returns from samplers such as
``LatinHypercube.Random`` and takes as covariance matrices, and the ``int[,]`` and ``float[,]``
tables of the graph, optimization and probability wrappers) copy directly into or out of the
rectangular .NET array: a .NET rectangular array is stored in row major order, the same layout
as a C-order NumPy array, so element ``[i, j]`` maps exactly.

Space-time blocks cross as a flat ``System.Single[]`` in C order together with explicit
``(steps, rows, columns)`` dimensions: ``flat[(t * rows + r) * columns + c]`` equals
``block[t, r, c]``.

pythonnet describes a ``System.Boolean[]`` buffer with the four-byte marshalled size of
``Boolean`` rather than its one-byte storage, so ``Boolean[]`` arrays cross as one buffer copy
into a ``System.Byte[]`` followed by one ``System.Buffer.BlockCopy``, which reinterprets the
bytes without conversion (a .NET ``Boolean`` occupies one byte, as a NumPy ``bool`` does).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from pynumerics.runtime import initialize

_INT32_MIN = -(2**31)
_INT32_MAX = 2**31 - 1
_INT64_MIN = -(2**63)
_INT64_MAX = 2**63 - 1
_INT8_MIN = -(2**7)
_INT8_MAX = 2**7 - 1


def _clr() -> tuple[Any, Any]:
    """Return the (System, Marshal) handles after guaranteeing runtime initialization."""
    initialize()
    import System  # noqa: PLC0415 - resolvable only after runtime initialization
    from System.Runtime.InteropServices import Marshal  # noqa: PLC0415

    return System, Marshal


_BOUND: dict[str, Any] = {}


def _bound() -> dict[str, Any]:
    """Return the array constructors and the ``Marshal.Copy`` read overloads, resolved once per process.

    pythonnet resolves an overloaded method group on every call; ``Marshal.Copy`` has more than a
    dozen overloads and ``Array.CreateInstance`` six, and resolving them costs many times the copy
    of a short array. Binding the exact overloads once keeps a small transfer in the microseconds.
    The ``Marshal.Copy`` overloads serve only reads of an argument that is not an array of the
    expected element type, which they convert or refuse as they always have.
    """
    if not _BOUND:
        System, Marshal = _clr()
        import clr  # noqa: PLC0415 - resolvable only after runtime initialization

        clr_type = clr.GetClrType
        pointer, length = System.IntPtr, System.Int32
        doubles, ints, singles = System.Array[System.Double], System.Array[System.Int32], System.Array[System.Single]
        longs, raw = System.Array[System.Int64], System.Array[System.Byte]
        create_matrix = System.Array.CreateInstance.Overloads[System.Type, length, length]
        _BOUND.update({
            "doubles": doubles, "ints": ints, "singles": singles, "pointer": pointer,
            "longs": longs, "bytes": raw, "sbytes": System.Array[System.SByte], "bools": System.Array[System.Boolean],
            "block_copy": System.Buffer.BlockCopy, "double_type": clr_type(System.Double),
            "create_matrix": create_matrix, "objects": System.Array[System.Object],
            "element_types": {np.dtype(np.float64): clr_type(System.Double), np.dtype(np.int32): clr_type(System.Int32),
                              np.dtype(np.float32): clr_type(System.Single)},
            "doubles_matrix": type(create_matrix(clr_type(System.Double), 0, 0)),
            "doubles_out": Marshal.Copy.Overloads[doubles, length, pointer, length],
            "ints_out": Marshal.Copy.Overloads[ints, length, pointer, length],
            "singles_out": Marshal.Copy.Overloads[singles, length, pointer, length],
            "longs_out": Marshal.Copy.Overloads[longs, length, pointer, length],
        })
    return _BOUND


def _fill(target: Any, array: np.ndarray) -> Any:
    """Copy ``array`` into the storage of ``target``, a new .NET primitive array of the same size, in C order.

    ``array`` must already have the native NumPy dtype of the element type of ``target`` (float64
    for ``Double``, int32 for ``Int32`` and so on), and ``target`` must hold ``array.size``
    elements, at any rank. Acquiring the buffer pins ``target``; the ``with`` block releases it,
    also when the copy raises. A same-dtype NumPy assignment copies the bits unchanged: one
    ``memcpy`` for contiguous input, a C-order strided copy otherwise. Returns ``target``.
    """
    if array.size:
        with memoryview(target) as view:
            if array.ndim == 1:
                np.frombuffer(view, dtype=array.dtype)[:] = array
            else:
                np.frombuffer(view, dtype=array.dtype).reshape(array.shape)[...] = array
    return target


def _read(values: Any, dtype: Any) -> np.ndarray:
    """Copy a .NET primitive array whose element type matches ``dtype`` into a new NumPy array of its shape.

    One ``memcpy`` of the pinned .NET storage; the ``with`` block releases the buffer. The caller
    guarantees the element type, since the buffer bytes are read as ``dtype``.
    """
    with memoryview(values) as view:
        shape = view.shape
        if len(shape) == 1:
            return np.frombuffer(view, dtype=dtype).copy()
        return np.frombuffer(view, dtype=dtype).reshape(shape).copy()


def _flat_doubles(array: np.ndarray) -> Any:
    """Copy a float64 array of any shape (contiguous, strided or broadcast) into a new ``System.Double[]`` in C order.

    The values cross in one assignment, with no intermediate contiguous copy; the elements are
    those of ``to_net_doubles(array.ravel())``, bit for bit.
    """
    return _fill(_bound()["doubles"](array.size), array)


def _net_rectangular(array: np.ndarray) -> Any:
    """Copy a two-dimensional float64, int32 or float32 array into a new .NET rectangular array of that element type."""
    bound = _bound()
    rows, cols = (int(dimension) for dimension in array.shape)
    return _fill(bound["create_matrix"](bound["element_types"][array.dtype], rows, cols), array)


_BUFFER_FORMATS = {np.dtype(np.float64): "d", np.dtype(np.int32): "i", np.dtype(np.float32): "f"}
"""The buffer-protocol format pythonnet reports for a .NET ``Double``, ``Int32`` and ``Single`` array."""


def _numpy_rectangular(values: Any, dtype: Any) -> np.ndarray:
    """Copy a .NET rectangular array (``double[,]``, ``int[,]`` or ``float[,]``) into a new 2-D array of ``dtype``.

    Raises
    ------
    TypeError
        If ``values`` is not a rank-2 array whose element type matches ``dtype``.
    """
    with memoryview(values) as view:
        if view.ndim != 2 or view.format != _BUFFER_FORMATS[np.dtype(dtype)]:
            raise TypeError(f"Expected a rank-2 .NET array of {np.dtype(dtype)}; got format {view.format!r} "
                            f"and rank {view.ndim}.")
        return np.frombuffer(view, dtype=dtype).reshape(view.shape).copy()


def to_net_doubles(values: Any) -> Any:
    """Copy a one-dimensional array into a new ``System.Double[]``, preserving every bit.

    The input is converted with ``numpy.ascontiguousarray(values, dtype=numpy.float64)``
    first. That conversion is exact for float32 values and for integers up to 2**53 in
    magnitude, and it turns a scalar into a one-element array. Starts the .NET runtime if
    it is not already running.

    Parameters
    ----------
    values : array_like
        Values to copy, one-dimensional after conversion.

    Returns
    -------
    System.Double[]
        A new .NET array holding the float64 values.

    Raises
    ------
    ValueError
        If ``values`` is not one-dimensional after conversion.
    """
    bound = _bound()
    array = np.ascontiguousarray(values, dtype=np.float64)
    if array.ndim == 0:
        array = array.reshape(1)
    if array.ndim != 1:
        raise ValueError(f"Expected a 1-D array; got shape {array.shape}.")
    return _fill(bound["doubles"](array.size), array)


def _optional_doubles(values: Any) -> Any:
    """Copy values into a ``System.Double[]``, passing None through as a null reference."""
    return None if values is None else to_net_doubles(values)


_LINEAR_ALGEBRA: list[Any] = []


def _linear_algebra() -> Any:
    """Return the ``Numerics.Mathematics.LinearAlgebra`` namespace after starting the runtime, imported once."""
    if not _LINEAR_ALGEBRA:
        initialize()
        import importlib  # noqa: PLC0415

        _LINEAR_ALGEBRA.append(importlib.import_module("Numerics.Mathematics.LinearAlgebra"))
    return _LINEAR_ALGEBRA[0]


_CONSTRUCTORS: dict[str, Any] = {}


def _construct(name: str, argument: Any) -> Any:
    """Construct a Numerics ``Matrix`` from a ``double[,]`` or a ``Vector`` from a ``double[]`` (``name``).

    pythonnet resolves a constructor overload on every call, about 30 microseconds for the
    ``Matrix`` constructors; invoking the one ``ConstructorInfo``, found once per process, runs
    the same constructor at about a quarter of that cost. Both constructors throw only for a null
    argument, and the callers pass a new array, so no exception reaches the reflection wrapper.
    """
    invoke = _CONSTRUCTORS.get(name)
    if invoke is None:
        bound = _bound()  # starts the runtime before the System import
        import clr  # noqa: PLC0415
        import System  # noqa: PLC0415

        owner = getattr(_linear_algebra(), name)
        argument_type = clr.GetClrType(bound["doubles_matrix"] if name == "Matrix" else bound["doubles"])
        constructor = clr.GetClrType(owner).GetConstructor(System.Array[System.Type]([argument_type]))
        invoke = _CONSTRUCTORS[name] = constructor.Invoke.Overloads[bound["objects"]]
    return invoke(_bound()["objects"]([argument]))


def _numerics_vector(values: Any, *, unwrap: bool = False, pass_through: bool = False) -> Any:
    """Return a Numerics ``Vector`` for ``values``, by default a new one copied from a 1-D array.

    With ``unwrap``, a wrapper whose ``net`` attribute holds a Numerics ``Vector`` (such as
    :class:`pynumerics.linalg.Vector`) gives that vector, not a copy. With ``pass_through``, a
    Numerics ``Vector`` is returned as it is, None gives None (a null reference), and any other
    input is flattened before the copy. Otherwise the values are copied with :func:`to_net_doubles`.
    """
    linear_algebra = _linear_algebra()
    if unwrap and isinstance(getattr(values, "net", None), linear_algebra.Vector):
        return values.net
    if pass_through:
        if values is None or isinstance(values, linear_algebra.Vector):
            return values
        values = np.ravel(np.asarray(values, dtype=np.float64))
    return _construct("Vector", to_net_doubles(values))


def to_numpy_doubles(values: Any) -> np.ndarray:
    """Copy a ``System.Double[]`` into a new float64 NumPy array, preserving every bit.

    The result is always a copy, so it never aliases .NET state. A ``System.Double[]`` is read
    with one copy through the buffer protocol. A .NET collection of doubles that is not an array
    (``List<double>``, ``ReadOnlyCollection<double>``) is first copied into a ``System.Double[]``
    with one ``ICollection<double>.CopyTo`` call, never element by element.

    Parameters
    ----------
    values : System.Double[] or System.Collections.Generic.ICollection[System.Double]
        .NET array or collection to copy.

    Returns
    -------
    numpy.ndarray
        A new one-dimensional float64 array with one element per value.
    """
    bound = _bound()
    doubles = bound["doubles"]
    if type(values) is not doubles:
        # An interface wrapper (IList<double>) exposes the object behind it.
        values = getattr(values, "__implementation__", values)
    if type(values) is doubles:
        with memoryview(values) as view:
            return np.frombuffer(view, dtype=np.float64).copy()
    if not hasattr(values, "Length"):
        array = doubles(int(values.Count))
        values.CopyTo(array, 0)
        return _read(array, np.float64)
    length = int(values.Length)
    out = np.empty(length, dtype=np.float64)
    if length:
        bound["doubles_out"](values, 0, bound["pointer"](out.ctypes.data), length)
    return out


def to_net_matrix(values: Any) -> Any:
    """Copy a two-dimensional array into a new ``System.Double[,]``, preserving every bit.

    Parameters
    ----------
    values : array_like
        Values to copy, two-dimensional after conversion to C-order float64.

    Returns
    -------
    System.Double[,]
        A new rectangular .NET array with the same shape.

    Raises
    ------
    ValueError
        If ``values`` is not two-dimensional after conversion.
    """
    array = np.ascontiguousarray(values, dtype=np.float64)
    if array.ndim != 2:
        raise ValueError(f"Expected a 2-D array; got shape {array.shape}.")
    return _net_rectangular(array)


def _numerics_matrix(values: Any, *, unwrap: bool = False, pass_through: bool = False) -> Any:
    """Return a Numerics ``Matrix`` for ``values``, by default a new one copied from a 2-D array.

    With ``unwrap``, a wrapper whose ``net`` attribute holds a Numerics ``Matrix`` (such as
    :class:`pynumerics.linalg.Matrix`) gives that matrix, not a copy. With ``pass_through``, a
    Numerics ``Matrix`` is returned as it is and any other input is promoted to at least two
    dimensions before the copy, so a scalar gives a 1 x 1 matrix. Otherwise the values are copied
    with :func:`to_net_matrix`.
    """
    linear_algebra = _linear_algebra()
    if unwrap and isinstance(getattr(values, "net", None), linear_algebra.Matrix):
        return values.net
    if pass_through:
        if isinstance(values, linear_algebra.Matrix):
            return values
        values = np.atleast_2d(np.asarray(values, dtype=np.float64))
    return _construct("Matrix", to_net_matrix(values))


def _parameter_set_arrays(sets: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a .NET ``ParameterSet[]`` into ``(values, fitness, weight)`` arrays, one row or entry per set.

    ``values`` has one row per set and as many columns as the longest set. A set without values
    (a null ``Values`` array) gives a row of NaN, and a shorter set is padded with NaN, so every
    value Numerics returned is kept. Numerics fills a failed replicate with NaN values of the full
    length, so such a row reads as NaN throughout.
    """
    count = len(sets)
    width = max((len(item.Values) for item in sets if item.Values is not None), default=0)
    values = np.full((count, width), np.nan)
    fitness, weight = np.empty(count), np.empty(count)
    for index, item in enumerate(sets):
        if item.Values is not None:
            row = to_numpy_doubles(item.Values)
            values[index, :row.size] = row
        fitness[index], weight[index] = float(item.Fitness), float(item.Weight)
    return values, fitness, weight


def to_numpy_matrix(values: Any) -> np.ndarray:
    """Copy a ``System.Double[,]`` into a new two-dimensional float64 array, preserving every bit.

    Parameters
    ----------
    values : System.Double[,]
        Rectangular .NET array to copy.

    Returns
    -------
    numpy.ndarray
        A new C-order float64 array of shape ``(values.GetLength(0), values.GetLength(1))``.

    Raises
    ------
    ValueError
        If ``values`` is not a rank-2 array.
    """
    bound = _bound()
    if type(values) is bound["doubles_matrix"]:
        return _read(values, np.float64)
    System, _ = _clr()
    if int(values.Rank) != 2:
        raise ValueError(f"Expected a rank-2 .NET array; got rank {int(values.Rank)}.")
    rows, cols = int(values.GetLength(0)), int(values.GetLength(1))
    size = rows * cols
    if size == 0:
        return np.empty((rows, cols), dtype=np.float64)
    if values.GetType().GetElementType() == bound["double_type"]:
        return _read(values, np.float64)
    flat = System.Array.CreateInstance(System.Double, size)
    System.Buffer.BlockCopy(values, 0, flat, 0, size * 8)
    return to_numpy_doubles(flat).reshape((rows, cols))


def _range_error(label: str, minimum: int, maximum: int, low: object, high: object) -> OverflowError:
    """Return the error for values from ``low`` to ``high`` that leave the ``label`` range."""
    return OverflowError(f"Values must lie in the {label} range [{minimum}, {maximum}]; got values from {low} to "
                         f"{high}.")


def _checked_integers(values: Any, dtype: type, label: str, minimum: int, maximum: int) -> np.ndarray:
    """Convert values to a contiguous integer array of ``dtype``, refusing any value outside its range.

    Booleans cross as 0 and 1 and floating-point values truncate toward zero; a value outside
    ``[minimum, maximum]``, or a floating-point value that is not finite, raises ``OverflowError``
    instead of wrapping. Complex, date, time, text and other values raise ``TypeError``.
    """
    source = np.asarray(values)
    kind = source.dtype.kind
    if kind == "O":
        items = source.ravel().tolist()
        if not all(isinstance(item, (int, np.integer)) for item in items):
            raise TypeError(f"Values must be booleans, integers, or floating-point numbers to cross as {label}; got "
                            "an object array holding another type.")
        if items and (min(items) < minimum or max(items) > maximum):
            raise _range_error(label, minimum, maximum, min(items), max(items))
        return np.ascontiguousarray(np.array(items, dtype=np.int64).reshape(source.shape), dtype=dtype)
    if kind not in "biuf":
        raise TypeError(f"Values must be booleans, integers, or floating-point numbers to cross as {label}; got dtype "
                        f"{source.dtype}.")
    if kind in "iuf" and source.size:
        if kind == "f":
            truncated = np.trunc(source)
            low, high = float(truncated.min()), float(truncated.max())
            if not (math.isfinite(low) and math.isfinite(high)):
                raise OverflowError(f"Values must lie in the {label} range [{minimum}, {maximum}]; "
                                    "got a value that is not finite.")
        else:
            low, high = int(source.min()), int(source.max())
        if low < minimum or high > maximum:
            raise _range_error(label, minimum, maximum, low, high)
    return np.ascontiguousarray(source, dtype=dtype)


def _int32_array(values: Any) -> np.ndarray:
    """Convert values to a contiguous int32 array, refusing any value outside the int32 range."""
    return _checked_integers(values, np.int32, "Int32", _INT32_MIN, _INT32_MAX)


def to_net_ints(values: Any) -> Any:
    """Copy a one-dimensional array into a new ``System.Int32[]``, preserving every bit.

    A value outside the int32 range raises ``OverflowError`` instead of wrapping,
    whether it arrives as a Python integer, a wider NumPy integer, or a floating-point
    value. Floating-point values inside the range truncate toward zero. Booleans cross
    as 0 and 1; complex, date, time, text and other values are refused.

    Parameters
    ----------
    values : array_like
        Integer values to copy, one-dimensional after conversion.

    Returns
    -------
    System.Int32[]
        A new .NET array holding the int32 values.

    Raises
    ------
    OverflowError
        If a value lies outside the int32 range or is not finite; raised before the
        runtime starts.
    TypeError
        If the values are not booleans, integers, or floating-point numbers.
    ValueError
        If ``values`` is not one-dimensional after conversion.
    """
    array = _int32_array(values)
    bound = _bound()
    if array.ndim != 1:
        raise ValueError(f"Expected a 1-D array; got shape {array.shape}.")
    return _fill(bound["ints"](array.size), array)


def to_numpy_ints(values: Any) -> np.ndarray:
    """Copy a ``System.Int32[]`` into a new int32 NumPy array, preserving every bit.

    Parameters
    ----------
    values : System.Int32[]
        .NET array to copy.

    Returns
    -------
    numpy.ndarray
        A new one-dimensional int32 array with ``values.Length`` elements.
    """
    bound = _bound()
    if type(values) is bound["ints"]:
        return _read(values, np.int32)
    length = int(values.Length)
    out = np.empty(length, dtype=np.int32)
    if length:
        bound["ints_out"](values, 0, bound["pointer"](out.ctypes.data), length)
    return out


def to_net_longs(values: Any) -> Any:
    """Copy a one-dimensional array into a new ``System.Int64[]``, preserving every bit.

    A value outside the int64 range raises ``OverflowError`` instead of wrapping, whether
    it arrives as a Python integer, an unsigned 64-bit NumPy integer, or a floating-point
    value. Floating-point values inside the range truncate toward zero. Booleans cross as
    0 and 1; complex, date, time, text and other values are refused.

    Parameters
    ----------
    values : array_like
        Integer values to copy, one-dimensional after conversion.

    Returns
    -------
    System.Int64[]
        A new .NET array holding the int64 values.

    Raises
    ------
    OverflowError
        If a value lies outside the int64 range or is not finite; raised before the
        runtime starts.
    TypeError
        If the values are not booleans, integers, or floating-point numbers.
    ValueError
        If ``values`` is not one-dimensional after conversion.
    """
    array = _checked_integers(values, np.int64, "Int64", _INT64_MIN, _INT64_MAX)
    bound = _bound()
    if array.ndim != 1:
        raise ValueError(f"Expected a 1-D array; got shape {array.shape}.")
    return _fill(bound["longs"](array.size), array)


def to_numpy_longs(values: Any) -> np.ndarray:
    """Copy a ``System.Int64[]`` into a new int64 NumPy array, preserving every bit.

    Parameters
    ----------
    values : System.Int64[]
        .NET array to copy, for example the ticks of a ``DateTime`` series.

    Returns
    -------
    numpy.ndarray
        A new one-dimensional int64 array with ``values.Length`` elements.
    """
    bound = _bound()
    if type(values) is bound["longs"]:
        return _read(values, np.int64)
    length = int(values.Length)
    out = np.empty(length, dtype=np.int64)
    if length:
        bound["longs_out"](values, 0, bound["pointer"](out.ctypes.data), length)
    return out


def _bytes_to_net(array: np.ndarray, kind: str) -> Any:
    """Copy a contiguous one-byte-per-element array into a new .NET array of ``kind`` through one byte copy."""
    bound = _bound()
    target = bound[kind](array.size)
    if kind != "bools":
        return _fill(target, array)
    if array.size:
        raw = _fill(bound["bytes"](array.size), array.view(np.uint8))
        bound["block_copy"](raw, 0, target, 0, array.size)
    return target


def _bytes_from_net(values: Any, dtype: type) -> np.ndarray:
    """Copy a one-byte-per-element .NET array into a new NumPy array of ``dtype`` through one byte copy."""
    bound = _bound()
    if type(values) is bound["bytes"] or type(values) is bound["sbytes"]:
        return _read(values, dtype)
    length = int(values.Length)
    raw = bound["bytes"](length)
    if length:
        bound["block_copy"](values, 0, raw, 0, length)
    return _read(raw, dtype)


def to_net_sbytes(values: Any) -> Any:
    """Copy a one-dimensional array into a new ``System.SByte[]``, preserving every bit.

    A value outside the int8 range raises ``OverflowError`` instead of wrapping.
    Floating-point values inside the range truncate toward zero. Booleans cross as 0 and 1;
    complex, date, time, text and other values are refused.

    Parameters
    ----------
    values : array_like
        Integer values to copy, one-dimensional after conversion.

    Returns
    -------
    System.SByte[]
        A new .NET array holding the int8 values.

    Raises
    ------
    OverflowError
        If a value lies outside the int8 range or is not finite; raised before the
        runtime starts.
    TypeError
        If the values are not booleans, integers, or floating-point numbers.
    ValueError
        If ``values`` is not one-dimensional after conversion.
    """
    array = _checked_integers(values, np.int8, "SByte", _INT8_MIN, _INT8_MAX)
    if array.ndim != 1:
        raise ValueError(f"Expected a 1-D array; got shape {array.shape}.")
    return _bytes_to_net(array, "sbytes")


def to_numpy_sbytes(values: Any) -> np.ndarray:
    """Copy a ``System.SByte[]`` into a new int8 NumPy array, preserving every bit.

    Parameters
    ----------
    values : System.SByte[]
        .NET array to copy.

    Returns
    -------
    numpy.ndarray
        A new one-dimensional int8 array with ``values.Length`` elements.
    """
    return _bytes_from_net(values, np.int8)


def to_net_bools(values: Any) -> Any:
    """Copy a one-dimensional boolean array into a new ``System.Boolean[]``.

    Only boolean input is accepted: integers or floating-point values are refused rather
    than interpreted as truth values. An empty input of any type gives an empty array.

    Parameters
    ----------
    values : array_like
        Boolean values to copy, one-dimensional after conversion.

    Returns
    -------
    System.Boolean[]
        A new .NET array holding the values.

    Raises
    ------
    TypeError
        If the values are not booleans; raised before the runtime starts.
    ValueError
        If ``values`` is not one-dimensional after conversion.
    """
    source = np.asarray(values)
    if source.dtype.kind != "b" and source.size:
        raise TypeError(f"Values must be booleans to cross as Boolean; got dtype {source.dtype}.")
    array = np.ascontiguousarray(source, dtype=np.bool_)
    if array.ndim != 1:
        raise ValueError(f"Expected a 1-D array; got shape {array.shape}.")
    return _bytes_to_net(array, "bools")


def to_numpy_bools(values: Any) -> np.ndarray:
    """Copy a ``System.Boolean[]`` into a new boolean NumPy array.

    Parameters
    ----------
    values : System.Boolean[]
        .NET array to copy.

    Returns
    -------
    numpy.ndarray
        A new one-dimensional bool array with ``values.Length`` elements; any nonzero
        .NET byte reads as True.
    """
    return _bytes_from_net(values, np.uint8) != 0


def to_net_singles(block: Any) -> tuple[Any, tuple[int, int, int]]:
    """Flatten a ``(steps, rows, columns)`` block into a ``System.Single[]`` and its shape.

    The block is converted with ``numpy.ascontiguousarray(block, dtype=numpy.float32)``
    first, which rounds float64 values to the nearest float32, and is then copied in C
    order, so ``flat[(t * rows + r) * columns + c]`` equals ``block[t, r, c]``.

    Parameters
    ----------
    block : array_like
        Block with shape ``(steps, rows, columns)``.

    Returns
    -------
    flat : System.Single[]
        A new .NET array holding the ``steps * rows * columns`` float32 values.
    dims : tuple of int
        The block shape ``(steps, rows, columns)``.

    Raises
    ------
    ValueError
        If ``block`` is not three-dimensional after conversion.
    """
    array = np.ascontiguousarray(block, dtype=np.float32)
    if array.ndim != 3:
        raise ValueError(f"Expected a (steps, rows, cols) block; got shape {array.shape}.")
    steps, rows, cols = (int(dimension) for dimension in array.shape)
    return _fill(_bound()["singles"](array.size), array), (steps, rows, cols)


def to_numpy_singles(values: Any, dims: tuple[int, int, int]) -> np.ndarray:
    """Copy a flat ``System.Single[]`` back into a ``(steps, rows, columns)`` float32 block.

    Parameters
    ----------
    values : System.Single[]
        Flat C-order .NET array, as produced by :func:`to_net_singles`.
    dims : tuple of int
        Block shape ``(steps, rows, columns)``.

    Returns
    -------
    numpy.ndarray
        A new float32 array with shape ``dims``, preserving every bit.

    Raises
    ------
    ValueError
        If ``values.Length`` differs from ``steps * rows * columns``.
    """
    bound = _bound()
    steps, rows, cols = (int(dimension) for dimension in dims)
    length = int(values.Length)
    if length != steps * rows * cols:
        raise ValueError(f"Flat length {length} does not match dims {dims}.")
    if type(values) is bound["singles"]:
        return _read(values, np.float32).reshape((steps, rows, cols))
    out = np.empty(length, dtype=np.float32)
    if length:
        bound["singles_out"](values, 0, bound["pointer"](out.ctypes.data), length)
    return out.reshape((steps, rows, cols))
