"""Buffer-protocol copies give the same .NET and NumPy arrays, bit for bit, as Marshal.Copy and element-wise copies.

The reference functions below are the ``Marshal.Copy`` transfers pynumerics used before it
copied through the buffer protocol, kept here only as the oracle; every comparison is on raw
bits (``view`` as unsigned integers), so NaN payloads, signalling NaNs and signed zeros count.
"""

from __future__ import annotations

import gc
from typing import Any

import numpy as np
import pytest

from pynumerics import marshal

_DTYPES = {np.dtype(np.float64): ("Double", np.uint64), np.dtype(np.float32): ("Single", np.uint32),
           np.dtype(np.int32): ("Int32", np.uint32), np.dtype(np.int64): ("Int64", np.uint64)}

_EDGE_DOUBLES = np.array([0x7FF8000000000001, 0xFFF8DEADBEEF0001, 0x7FF0000000000001, 0x7FF4000000000000,
                          0x8000000000000000, 0x0000000000000000, 0x0000000000000001, 0x000FFFFFFFFFFFFF,
                          0x7FF0000000000000, 0xFFF0000000000000, 0x7FEFFFFFFFFFFFFF, 0x0010000000000000,
                          0x3FB999999999999A], dtype=np.uint64).view(np.float64)
_EDGE_SINGLES = np.array([0x7FC00001, 0xFFC0BEEF, 0x7F800001, 0x80000000, 0x00000000, 0x00000001, 0x007FFFFF,
                          0x7F800000, 0xFF800000, 0x7F7FFFFF, 0x3DCCCCCD], dtype=np.uint32).view(np.float32)


def _system() -> Any:
    marshal._bound()
    import System  # noqa: PLC0415 - resolvable once the runtime has started

    return System


def _reference_to_net(values: np.ndarray) -> Any:
    """Copy a contiguous 1-D array into a new .NET array with one ``Marshal.Copy`` (the previous transfer)."""
    System = _system()
    from System.Runtime.InteropServices import Marshal  # noqa: PLC0415

    array = np.ascontiguousarray(values)
    element = getattr(System, _DTYPES[array.dtype][0])
    kind = System.Array[element]
    target = kind(array.size)
    if array.size:
        copy = Marshal.Copy.Overloads[System.IntPtr, kind, System.Int32, System.Int32]
        copy(System.IntPtr(array.ctypes.data), target, 0, array.size)
    return target


def _reference_to_numpy(values: Any, dtype: Any) -> np.ndarray:
    """Copy a 1-D .NET array into NumPy with one ``Marshal.Copy`` (the previous transfer)."""
    System = _system()
    from System.Runtime.InteropServices import Marshal  # noqa: PLC0415

    kind = System.Array[getattr(System, _DTYPES[np.dtype(dtype)][0])]
    out = np.empty(int(values.Length), dtype=dtype)
    if out.size:
        copy = Marshal.Copy.Overloads[kind, System.Int32, System.IntPtr, System.Int32]
        copy(values, 0, System.IntPtr(out.ctypes.data), out.size)
    return out


def _reference_flatten(values: Any, dtype: Any) -> np.ndarray:
    """Read a rectangular .NET array through a flat copy and ``Buffer.BlockCopy`` (the previous 2-D transfer)."""
    System = _system()
    rows, cols = int(values.GetLength(0)), int(values.GetLength(1))
    flat = System.Array[getattr(System, _DTYPES[np.dtype(dtype)][0])](rows * cols)
    if rows * cols:
        System.Buffer.BlockCopy(values, 0, flat, 0, rows * cols * np.dtype(dtype).itemsize)
    return _reference_to_numpy(flat, dtype).reshape(rows, cols)


def _bits(values: np.ndarray) -> np.ndarray:
    """Return the raw bits of an array as unsigned integers of the same width."""
    return np.ascontiguousarray(values).view(_DTYPES[values.dtype][1])


def _assert_same_bits(actual: np.ndarray, expected: np.ndarray) -> None:
    assert actual.dtype == expected.dtype and actual.shape == expected.shape
    np.testing.assert_array_equal(_bits(actual), _bits(expected))


def _inputs(dtype: Any) -> list[np.ndarray]:
    """Return edge, random, empty, strided and large inputs of ``dtype``."""
    rng = np.random.default_rng(20261108)
    if np.dtype(dtype) == np.float64:
        base = np.concatenate([_EDGE_DOUBLES, rng.standard_normal(97) * 1e300, rng.standard_normal(50) * 1e-310])
    elif np.dtype(dtype) == np.float32:
        base = np.concatenate([_EDGE_SINGLES, rng.standard_normal(97).astype(np.float32)])
    else:
        info = np.iinfo(dtype)
        base = np.concatenate([np.array([info.min, info.min + 1, -1, 0, 1, info.max - 1, info.max], dtype=dtype),
                               rng.integers(info.min, info.max, 97, dtype=dtype, endpoint=True)])
    large = np.resize(base, 100_003)
    return [base, base[::-3], base[1:2], base[:0], large, large[::7]]


@pytest.mark.parametrize("dtype", [np.float64, np.int32, np.int64])
def test_one_dimensional_writes_and_reads_match_marshal_copy_bit_for_bit(dtype: Any) -> None:
    # Arrange
    write = {np.float64: marshal.to_net_doubles, np.int32: marshal.to_net_ints, np.int64: marshal.to_net_longs}[dtype]
    read = {np.float64: marshal.to_numpy_doubles, np.int32: marshal.to_numpy_ints,
            np.int64: marshal.to_numpy_longs}[dtype]
    for values in _inputs(dtype):
        # Act
        net = write(values)
        reference = _reference_to_net(values)
        # Assert: the .NET arrays hold the same bits, and both read paths agree
        assert int(net.Length) == values.size
        _assert_same_bits(_reference_to_numpy(net, dtype), _reference_to_numpy(reference, dtype))
        _assert_same_bits(read(reference), _reference_to_numpy(reference, dtype))
        _assert_same_bits(read(net), np.ascontiguousarray(values))


def test_double_writes_match_element_by_element_construction() -> None:
    # Arrange: integers and ordinary doubles convert exactly through pythonnet indexing
    System = _system()
    values = np.concatenate([np.array([-0.0, 0.0, np.inf, -np.inf, 5e-324, 1.7976931348623157e308]),
                             np.random.default_rng(3).standard_normal(64)])
    ints = np.random.default_rng(4).integers(-(2**31), 2**31, 64).astype(np.int32)
    longs = np.random.default_rng(5).integers(-(2**63), 2**63 - 1, 64, dtype=np.int64)
    # Act
    expected_doubles = System.Array[System.Double](values.size)
    expected_ints = System.Array[System.Int32](ints.size)
    expected_longs = System.Array[System.Int64](longs.size)
    for index in range(values.size):
        expected_doubles[index] = float(values[index])
    for index in range(ints.size):
        expected_ints[index] = int(ints[index])
        expected_longs[index] = int(longs[index])
    # Assert
    _assert_same_bits(marshal.to_numpy_doubles(marshal.to_net_doubles(values)),
                      _reference_to_numpy(expected_doubles, np.float64))
    _assert_same_bits(_reference_to_numpy(marshal.to_net_ints(ints), np.int32),
                      _reference_to_numpy(expected_ints, np.int32))
    _assert_same_bits(_reference_to_numpy(marshal.to_net_longs(longs), np.int64),
                      _reference_to_numpy(expected_longs, np.int64))


def test_conversions_before_the_copy_are_unchanged() -> None:
    # Arrange: float32, integer, boolean and scalar inputs convert to float64 or int before crossing
    singles = _EDGE_SINGLES[np.isfinite(_EDGE_SINGLES)]
    cases = [singles, np.arange(-5, 5, dtype=np.int16), np.array([True, False]), 2.5, [1, 2.5, -3]]
    for values in cases:
        # Act
        net = marshal.to_net_doubles(values)
        # Assert
        _assert_same_bits(_reference_to_numpy(net, np.float64),
                          _reference_to_numpy(_reference_to_net(np.ascontiguousarray(values, dtype=np.float64)),
                                              np.float64))
    _assert_same_bits(marshal.to_numpy_ints(marshal.to_net_ints([2.9, -2.9, True])), np.array([2, -2, 1], np.int32))
    with pytest.raises(OverflowError):
        marshal.to_net_ints(np.array([2**31], dtype=np.int64))
    with pytest.raises(OverflowError, match="not finite"):
        marshal.to_net_ints([np.inf])
    with pytest.raises(TypeError):
        marshal.to_net_ints(np.array([1 + 1j]))


def test_single_blocks_match_marshal_copy_bit_for_bit() -> None:
    # Arrange
    rng = np.random.default_rng(6)
    blocks = [np.resize(_EDGE_SINGLES, (3, 4, 5)), rng.standard_normal((2, 3, 4)).astype(np.float32)[:, ::-1, ::2],
              np.asfortranarray(rng.standard_normal((4, 3, 2))), np.zeros((0, 2, 3), dtype=np.float32)]
    for block in blocks:
        # Act
        flat, dims = marshal.to_net_singles(block)
        back = marshal.to_numpy_singles(flat, dims)
        expected = np.ascontiguousarray(block, dtype=np.float32)
        # Assert
        assert dims == expected.shape
        _assert_same_bits(_reference_to_numpy(flat, np.float32), expected.ravel())
        _assert_same_bits(back, expected)
        _assert_same_bits(marshal.to_numpy_singles(_reference_to_net(expected.ravel()), dims), expected)


@pytest.mark.parametrize("shape", [(5, 7), (1, 1), (0, 3), (3, 0), (316, 316)])
def test_matrices_match_block_copy_bit_for_bit_in_every_memory_order(shape: tuple[int, int]) -> None:
    # Arrange
    rng = np.random.default_rng(7)
    values = np.resize(np.concatenate([_EDGE_DOUBLES, rng.standard_normal(200)]), shape)
    wide = np.zeros((shape[0], 2 * shape[1]))
    wide[:, ::2] = values
    inputs = [values, np.asfortranarray(values), values.T.copy().T, wide[:, ::2]]
    for array in inputs:
        # Act
        net = marshal.to_net_matrix(array)
        back = marshal.to_numpy_matrix(net)
        # Assert: element [i, j] maps exactly in both directions
        assert (int(net.GetLength(0)), int(net.GetLength(1))) == shape
        _assert_same_bits(_reference_flatten(net, np.float64), np.ascontiguousarray(values))
        _assert_same_bits(back, np.ascontiguousarray(values))
        assert back.flags.owndata and back.flags.c_contiguous


def test_integer_and_single_tables_match_block_copy_bit_for_bit() -> None:
    # Arrange
    from pynumerics import graphs, optimization, probability, special  # noqa: PLC0415

    rng = np.random.default_rng(8)
    ints = rng.integers(-(2**31), 2**31, (6, 4)).astype(np.int32)
    table = np.resize(_EDGE_SINGLES, (7, 3))
    # Act
    net_ints = optimization._ints_matrix_to_net(np.asfortranarray(ints))
    net_indicators = probability._int_matrix(ints[:, ::-1])
    net_table = graphs._table_to_net(table)
    combinations = special.all_combinations(4)
    # Assert
    _assert_same_bits(_reference_flatten(net_ints, np.int32), ints)
    _assert_same_bits(optimization._ints_matrix_to_numpy(net_ints), ints)
    _assert_same_bits(_reference_flatten(net_indicators, np.int32), np.ascontiguousarray(ints[:, ::-1]))
    _assert_same_bits(graphs._table_to_numpy(net_table), table)
    _assert_same_bits(_reference_flatten(net_table, np.float32), table)
    assert combinations.shape == (15, 4) and combinations.dtype == np.int32
    with pytest.raises(OverflowError):
        optimization._ints_matrix_to_net(np.array([[2**31]]))
    with pytest.raises(TypeError, match="rank-2"):
        marshal._numpy_rectangular(net_ints, np.float64)
    with pytest.raises(TypeError, match="rank-2"):
        marshal._numpy_rectangular(marshal.to_net_ints([1, 2]), np.int32)


def test_byte_sized_arrays_match_block_copy() -> None:
    # Arrange
    System = _system()
    sbytes = np.arange(-128, 128, dtype=np.int8)
    flags = np.random.default_rng(9).integers(0, 2, 37).astype(bool)
    # Act
    net_sbytes = marshal.to_net_sbytes(sbytes[::-1])
    net_flags = marshal.to_net_bools(flags)
    raw = System.Array[System.Byte](sbytes.size)
    System.Buffer.BlockCopy(net_sbytes, 0, raw, 0, sbytes.size)
    # Assert
    assert [int(raw[i]) for i in range(3)] == [127, 126, 125]
    np.testing.assert_array_equal(marshal.to_numpy_sbytes(net_sbytes), sbytes[::-1])
    assert [bool(net_flags[i]) for i in range(flags.size)] == flags.tolist()
    np.testing.assert_array_equal(marshal.to_numpy_bools(net_flags), flags)
    assert marshal._bytes_from_net(marshal._bytes_to_net(np.frombuffer(b"\x00\xffab", np.uint8), "bytes"),
                                   np.uint8).tobytes() == b"\x00\xffab"


def test_reads_refuse_other_element_types_as_before() -> None:
    # Arrange
    doubles = marshal.to_net_doubles([1.0, 2.0])
    ints = marshal.to_net_ints([1, 2])
    # Act / Assert: an array of another element type is converted or refused by Marshal.Copy, never reinterpreted
    with pytest.raises(TypeError):
        marshal.to_numpy_ints(doubles)
    with pytest.raises(TypeError):
        marshal.to_numpy_longs(ints)
    with pytest.raises(TypeError):
        marshal.to_numpy_doubles(ints)
    with pytest.raises(ValueError, match="rank-2"):
        marshal.to_numpy_matrix(doubles)


def test_buffers_are_released_after_every_copy_and_after_a_failed_copy() -> None:
    # Arrange
    System = _system()
    target = System.Array[System.Double](3)
    # Act: a copy whose shapes disagree raises inside the pinned block
    with pytest.raises(ValueError):
        marshal._fill(target, np.arange(4.0))
    for _ in range(1000):
        marshal.to_numpy_doubles(marshal.to_net_doubles(np.arange(8.0)))
    System.GC.Collect()
    gc.collect()
    # Assert: the array can be exported, released and written again
    with memoryview(target) as view:
        assert view.nbytes == 24
    marshal._fill(target, np.array([1.0, -0.0, np.nan]))
    np.testing.assert_array_equal(_bits(marshal.to_numpy_doubles(target)),
                                  _bits(np.array([1.0, -0.0, np.nan])))


def test_dispatch_helpers_keep_their_results() -> None:
    # Arrange
    from pynumerics import linalg, special  # noqa: PLC0415
    from pynumerics.distributions import Normal  # noqa: PLC0415
    from pynumerics.multivariate import Multinomial  # noqa: PLC0415

    x = np.linspace(-3.0, 3.0, 12).reshape(3, 4)
    a, b = np.array([0.5, 1.0, 2.5]), np.array([[1.0], [2.0]])
    normal = Normal(0.0, 1.0)
    cumulative = np.array([1.0, 3.0, 3.0, 6.0])
    uniforms = np.array([[0.0, 1.0, 2.9], [3.0, 5.99, 6.0]])
    # Act / Assert: arrays give the scalar member's bits at every element, in the input's shape
    assert normal.cdf(x).shape == (3, 4)
    _assert_same_bits(normal.cdf(x), np.array([[normal.cdf(float(v)) for v in row] for row in x]))
    positive = (x + 3.25)[:, ::2]
    _assert_same_bits(special.gamma_function(positive),
                      np.array([[special.gamma_function(float(v)) for v in row] for row in positive]))
    broadcast = special.beta_function(a, b)
    assert broadcast.shape == (2, 3)
    _assert_same_bits(broadcast, np.array([[special.beta_function(float(i), float(j[0])) for i in a] for j in b]))
    exclusive = Multinomial.categorical_indices(cumulative, uniforms)
    inclusive = Multinomial.categorical_indices(cumulative, uniforms, inclusive=True)
    np.testing.assert_array_equal(exclusive, np.clip(np.searchsorted(cumulative, uniforms, "right"), 0, 3))
    np.testing.assert_array_equal(inclusive, np.clip(np.searchsorted(cumulative, uniforms, "left"), 0, 3))
    assert Multinomial.categorical_indices(cumulative, 3.0) == 3
    matrix = np.array([[2.0, -1.0], [0.5, 4.0]])
    np.testing.assert_array_equal(np.asarray(linalg.Matrix(matrix).multiply([1.0, -2.0])), matrix @ [1.0, -2.0])
    np.testing.assert_array_equal(np.asarray(linalg.Matrix(matrix).multiply(matrix)), matrix @ matrix)
