"""Bulk marshaling keeps every bit, including NaN payloads, signed zeros and subnormals."""

from __future__ import annotations

import numpy as np
import pytest

from pynumerics import marshal

PATTERNS = np.array(
    [
        0x7FF8000000000001,  # quiet NaN with a payload
        0xFFF8DEADBEEF0001,  # negative quiet NaN with a payload
        0x7FF0000000000001,  # signalling NaN pattern
        0x8000000000000000,  # -0.0
        0x0000000000000000,  # +0.0
        0x0000000000000001,  # smallest subnormal
        0x000FFFFFFFFFFFFF,  # largest subnormal
        0x7FF0000000000000,  # +inf
        0xFFF0000000000000,  # -inf
        0x7FEFFFFFFFFFFFFF,  # largest finite
        0x3FB999999999999A,  # 0.1
    ],
    dtype=np.uint64,
)


def test_doubles_round_trip_bit_exact() -> None:
    # Arrange
    values = PATTERNS.view(np.float64)
    # Act
    back = marshal.to_numpy_doubles(marshal.to_net_doubles(values))
    # Assert
    np.testing.assert_array_equal(back.view(np.uint64), PATTERNS)


def test_doubles_are_bit_exact_inside_net() -> None:
    # Arrange
    values = PATTERNS.view(np.float64)
    # Act
    net = marshal.to_net_doubles(values)
    import System  # resolvable once the runtime has started

    bits = [int(System.BitConverter.DoubleToInt64Bits(net[i])) & 0xFFFFFFFFFFFFFFFF for i in range(net.Length)]
    # Assert
    assert bits == [int(p) for p in PATTERNS]


def test_net_collections_copy_out_bit_exact() -> None:
    # Arrange
    values = PATTERNS.view(np.float64)
    net = marshal.to_net_doubles(values)
    import System  # resolvable once the runtime has started
    from System.Collections.Generic import List  # noqa: PLC0415
    from System.Collections.ObjectModel import ReadOnlyCollection  # noqa: PLC0415

    as_list = List[System.Double](net)
    # Act
    from_list = marshal.to_numpy_doubles(as_list)
    from_read_only = marshal.to_numpy_doubles(ReadOnlyCollection[System.Double](as_list))
    empty = marshal.to_numpy_doubles(List[System.Double]())
    # Assert
    np.testing.assert_array_equal(from_list.view(np.uint64), PATTERNS)
    np.testing.assert_array_equal(from_read_only.view(np.uint64), PATTERNS)
    assert empty.shape == (0,) and empty.dtype == np.float64


def test_buffer_reads_are_owned_copies_and_refuse_other_element_types() -> None:
    # Arrange
    values = PATTERNS.view(np.float64)
    net = marshal.to_net_doubles(values)
    import System  # resolvable once the runtime has started
    from System.Collections.Generic import IList  # noqa: PLC0415

    as_interface = IList[System.Double](net)
    # Act
    first = marshal.to_numpy_doubles(net)
    through_interface = marshal.to_numpy_doubles(as_interface)
    first[0] = 1.0
    again = marshal.to_numpy_doubles(net)
    # Assert
    assert type(as_interface).__name__.startswith("IList")
    assert first.flags.owndata and first.flags.writeable and first.flags.c_contiguous
    np.testing.assert_array_equal(through_interface.view(np.uint64), PATTERNS)
    np.testing.assert_array_equal(again.view(np.uint64), PATTERNS)
    assert net[0] != 1.0
    with pytest.raises(TypeError):
        marshal.to_numpy_doubles(marshal.to_net_ints([1, 2]))


def test_matrix_reads_are_owned_row_major_copies() -> None:
    # Arrange
    values = np.arange(6.0).reshape(2, 3)
    net = marshal.to_net_matrix(values)
    # Act
    copied = marshal.to_numpy_matrix(net)
    copied[0, 0] = -1.0
    # Assert
    np.testing.assert_array_equal(marshal.to_numpy_matrix(net), values)
    assert copied.flags.owndata and copied.flags.writeable and copied.flags.c_contiguous


def test_empty_and_scalar_inputs() -> None:
    # Act / Assert
    assert marshal.to_numpy_doubles(marshal.to_net_doubles([])).size == 0
    assert marshal.to_numpy_doubles(marshal.to_net_doubles(2.5)).tolist() == [2.5]
    with pytest.raises(ValueError, match="1-D"):
        marshal.to_net_doubles(np.zeros((2, 2)))


def test_matrix_round_trip_bit_exact_and_row_major() -> None:
    # Arrange
    values = PATTERNS[:10].view(np.float64).reshape(2, 5)
    # Act
    net = marshal.to_net_matrix(values)
    back = marshal.to_numpy_matrix(net)
    # Assert
    assert (int(net.GetLength(0)), int(net.GetLength(1))) == (2, 5)
    assert np.float64(net[1, 3]).view(np.uint64) == values[1, 3].view(np.uint64)
    np.testing.assert_array_equal(back.view(np.uint64), values.view(np.uint64))


def test_ints_round_trip_and_range_refusal() -> None:
    # Arrange
    values = np.array([-(2**31), -1, 0, 1, 2**31 - 1], dtype=np.int64)
    # Act / Assert
    np.testing.assert_array_equal(marshal.to_numpy_ints(marshal.to_net_ints(values)), values.astype(np.int32))
    with pytest.raises(OverflowError):
        marshal.to_net_ints([2**31])
    with pytest.raises(TypeError):
        marshal.to_net_ints(np.array(["a"]))


def test_singles_block_round_trip_c_order() -> None:
    # Arrange
    block = np.arange(24, dtype=np.float32).reshape(2, 3, 4)
    block[1, 2, 3] = np.float32(np.nan)
    # Act
    flat, dims = marshal.to_net_singles(block)
    back = marshal.to_numpy_singles(flat, dims)
    # Assert
    assert dims == (2, 3, 4)
    assert float(flat[(1 * 3 + 0) * 4 + 2]) == float(block[1, 0, 2])
    np.testing.assert_array_equal(back.view(np.uint32), block.view(np.uint32))


def test_longs_round_trip_bit_exact_and_range_refusal() -> None:
    # Arrange
    values = np.array([-(2**63), -1, 0, 1, 637_134_336_000_000_000, 2**63 - 1], dtype=np.int64)
    # Act
    net = marshal.to_net_longs(values)
    back = marshal.to_numpy_longs(net)
    # Assert
    assert int(net[5]) == 2**63 - 1 and int(net[0]) == -(2**63)
    np.testing.assert_array_equal(back, values)
    assert back.dtype == np.int64
    assert marshal.to_numpy_longs(marshal.to_net_longs([])).size == 0
    np.testing.assert_array_equal(marshal.to_numpy_longs(marshal.to_net_longs([2**40, -3.9])), [2**40, -3])
    with pytest.raises(OverflowError, match="Int64"):
        marshal.to_net_longs([2**63])
    with pytest.raises(OverflowError, match="Int64"):
        marshal.to_net_longs(np.array([2**64 - 1], dtype=np.uint64))
    with pytest.raises(OverflowError, match="not finite"):
        marshal.to_net_longs([np.nan])
    with pytest.raises(TypeError):
        marshal.to_net_longs(np.array(["a"]))
    with pytest.raises(ValueError, match="1-D"):
        marshal.to_net_longs(np.zeros((2, 2), dtype=np.int64))


def test_sbytes_round_trip_bit_exact_and_range_refusal() -> None:
    # Arrange
    values = np.array([-128, -1, 0, 1, 127], dtype=np.int8)
    # Act
    net = marshal.to_net_sbytes(values)
    back = marshal.to_numpy_sbytes(net)
    # Assert
    assert [int(net[i]) for i in range(net.Length)] == [-128, -1, 0, 1, 127]
    np.testing.assert_array_equal(back, values)
    assert back.dtype == np.int8
    assert marshal.to_numpy_sbytes(marshal.to_net_sbytes([])).size == 0
    with pytest.raises(OverflowError, match="SByte"):
        marshal.to_net_sbytes([128])
    with pytest.raises(OverflowError, match="SByte"):
        marshal.to_net_sbytes(np.array([-129], dtype=np.int16))


def test_bools_round_trip_and_type_refusal() -> None:
    # Arrange
    values = np.array([True, False, False, True, True], dtype=bool)
    # Act
    net = marshal.to_net_bools(values)
    back = marshal.to_numpy_bools(net)
    # Assert
    assert [bool(net[i]) for i in range(net.Length)] == values.tolist()
    np.testing.assert_array_equal(back, values)
    assert back.dtype == np.bool_
    assert marshal.to_numpy_bools(marshal.to_net_bools([])).size == 0
    with pytest.raises(TypeError, match="booleans"):
        marshal.to_net_bools([1, 0])
    with pytest.raises(ValueError, match="1-D"):
        marshal.to_net_bools(np.zeros((2, 2), dtype=bool))


def test_numerics_matrix_and_vector_copy_unwrap_or_pass_through() -> None:
    # Arrange
    from pynumerics import linalg  # noqa: PLC0415

    values = np.array([[4.0, 1.0], [1.0, 3.0]])
    wrapped, wrapped_vector = linalg.Matrix(values), linalg.Vector([1.0, 2.0])
    net = marshal._numerics_matrix(values)
    net_vector = marshal._numerics_vector([1.0, 2.0])
    # Act / Assert: by default every input is copied into a new object
    assert marshal._numerics_matrix(values) is not net
    np.testing.assert_array_equal(marshal.to_numpy_matrix(net.ToArray()), values)
    with pytest.raises(ValueError, match="2-D"):
        marshal._numerics_matrix([1.0, 2.0])
    # unwrap reuses the Numerics object a wrapper holds, and copies anything else
    assert marshal._numerics_matrix(wrapped, unwrap=True) is wrapped.net
    assert marshal._numerics_vector(wrapped_vector, unwrap=True) is wrapped_vector.net
    np.testing.assert_array_equal(marshal.to_numpy_matrix(marshal._numerics_matrix(values, unwrap=True).ToArray()),
                                  values)
    # pass_through reuses a Numerics object, maps None to null and promotes or flattens arrays
    assert marshal._numerics_matrix(net, pass_through=True) is net
    assert marshal._numerics_vector(net_vector, pass_through=True) is net_vector
    assert marshal._numerics_vector(None, pass_through=True) is None
    assert marshal.to_numpy_matrix(marshal._numerics_matrix(0.25, pass_through=True).ToArray()).tolist() == [[0.25]]
    flat = marshal._numerics_vector([[1.0], [2.0]], pass_through=True)
    assert marshal.to_numpy_doubles(flat.Array).tolist() == [1.0, 2.0]


def test_parameter_set_arrays_keep_every_value_and_pad_with_nan() -> None:
    # Arrange
    from pynumerics._util import _opt  # noqa: PLC0415

    ps = _opt().ParameterSet
    sets = [ps(marshal.to_net_doubles([1.0, 2.0, 3.0]), 0.5, 2.0), ps(marshal.to_net_doubles([4.0, 5.0]), 1.5, 1.0),
            ps()]
    # Act
    values, fitness, weight = marshal._parameter_set_arrays(sets)
    empty = marshal._parameter_set_arrays([])
    # Assert: a shorter set is padded, a set without values is NaN throughout
    np.testing.assert_array_equal(values, [[1.0, 2.0, 3.0], [4.0, 5.0, np.nan], [np.nan] * 3])
    np.testing.assert_array_equal(fitness[:2], [0.5, 1.5])
    np.testing.assert_array_equal(weight[:2], [2.0, 1.0])
    assert empty[0].shape == (0, 0) and empty[1].size == 0
