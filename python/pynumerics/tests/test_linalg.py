"""Matrices, vectors and decompositions against numpy.linalg.

Elementwise operations that numpy performs with the same arithmetic (sums, differences, scaling,
transposes, triangles) must match bit for bit; products and factorizations match to 1e-14 relative
to the matrix scale, since BLAS and LAPACK order their sums differently.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from pynumerics import linalg
from pynumerics._util import namespace
from pynumerics.marshal import to_net_matrix


def _bits(values: object) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64)


def _general(n: int, seed: int) -> np.ndarray:
    return np.random.default_rng(seed).normal(size=(n, n))


def _spd(n: int, seed: int) -> np.ndarray:
    a = _general(n, seed)
    return a @ a.T + n * np.eye(n)


def _close(actual: object, expected: object, scale: float = 1.0, rtol: float = 1e-14) -> None:
    np.testing.assert_allclose(actual, expected, rtol=0, atol=rtol * scale)


def test_matrix_constructors_and_accessors() -> None:
    # Arrange
    a = _general(4, 1)
    # Act
    m = linalg.Matrix(a)
    # Assert
    np.testing.assert_array_equal(_bits(m.to_array()), _bits(a))
    np.testing.assert_array_equal(np.asarray(m), a)
    np.testing.assert_array_equal(m.array, a)
    assert m.shape == (4, 4) and m.number_of_rows == 4 and m.number_of_columns == 4 and m.is_square
    assert m[1, 2] == a[1, 2]
    m[1, 2] = 7.5
    assert m[1, 2] == 7.5 and a[1, 2] != 7.5
    assert linalg.Matrix(3).shape == (3, 3) and linalg.Matrix(2, 5).shape == (2, 5)
    np.testing.assert_array_equal(linalg.Matrix([1.0, 2.0]).to_array(), [[1.0], [2.0]])
    np.testing.assert_array_equal(linalg.Matrix(linalg.Vector([1.0, 2.0])).to_array(), [[1.0], [2.0]])
    np.testing.assert_array_equal(linalg.Matrix([[1.0, 2.0], [3.0, 4.0]], by_column=True).to_array(),
                                  [[1.0, 3.0], [2.0, 4.0]])
    np.testing.assert_array_equal(linalg.Matrix([[1.0, 2.0], [3.0, 4.0]], by_column=False).to_array(),
                                  [[1.0, 2.0], [3.0, 4.0]])
    np.testing.assert_array_equal(m.row(0), a[0])
    np.testing.assert_array_equal(m.column(3), m.to_array()[:, 3])
    np.testing.assert_array_equal(m.diagonal(), np.diag(m.to_array()))


def test_matrix_xml_header_clone_and_symmetry() -> None:
    # Arrange
    s = _spd(3, 2)
    m = linalg.Matrix(s)
    # Act
    m.header = ["a", "b", "c"]
    copy = m.clone()
    copy[0, 0] = -1.0
    restored = linalg.Matrix.from_xml(m.to_xml())
    # Assert
    assert m.header == ["a", "b", "c"]
    assert m.is_symmetric and not linalg.Matrix(_general(3, 3)).is_symmetric
    assert m[0, 0] == s[0, 0]
    np.testing.assert_array_equal(_bits(restored.to_array()), _bits(s))
    assert linalg.Matrix.from_net(m.net).net is m.net


def test_matrix_elementwise_operations_match_numpy_bit_for_bit() -> None:
    # Arrange
    a, b = _general(5, 4), _general(5, 5)
    m = linalg.Matrix(a)
    # Act / Assert
    np.testing.assert_array_equal(_bits(np.asarray(m.add(b))), _bits(a + b))
    np.testing.assert_array_equal(_bits(np.asarray(m.subtract(linalg.Matrix(b)))), _bits(a - b))
    np.testing.assert_array_equal(_bits(np.asarray(m.multiply(2.5))), _bits(a * 2.5))
    np.testing.assert_array_equal(_bits(np.asarray(m.divide(3.0))), _bits(a / 3.0))
    np.testing.assert_array_equal(_bits(np.asarray(m.transpose())), _bits(a.T))
    np.testing.assert_array_equal(np.asarray(m.upper_triangle()), np.triu(a))
    np.testing.assert_array_equal(np.asarray(m.lower_triangle()), np.tril(a))
    np.testing.assert_array_equal(np.asarray(linalg.Matrix.identity(4)), np.eye(4))
    np.testing.assert_array_equal(np.asarray(linalg.Matrix.diagonal_matrix([1.0, 2.0])), np.diag([1.0, 2.0]))
    np.testing.assert_array_equal(np.asarray(linalg.Matrix.diagonal_matrix(a)), np.diag(np.diag(a)))
    np.testing.assert_array_equal(np.asarray(linalg.Matrix.outer([1.0, 2.0], [3.0, 4.0, 5.0])),
                                  np.outer([1.0, 2.0], [3.0, 4.0, 5.0]))
    np.testing.assert_allclose(m.sum(), a.sum(), rtol=1e-14)
    np.testing.assert_allclose(m.trace(), np.trace(a), rtol=1e-14)
    np.testing.assert_allclose(m.column_means(), a.mean(axis=0), rtol=1e-14, atol=1e-16)


def test_matrix_pointwise_functions_match_numpy() -> None:
    # Arrange
    a = np.abs(_general(4, 6)) + 0.1
    expected = {"sqrt": np.sqrt(a), "sqr": a * a, "log": np.log(a), "exp": np.exp(a)}
    for name, values in expected.items():
        m = linalg.Matrix(a)
        # Act
        getattr(m, name)()
        # Assert
        np.testing.assert_allclose(m.to_array(), values, rtol=1e-15)
    m = linalg.Matrix(a)
    m.apply(lambda v: 2.0 * v + 1.0)
    np.testing.assert_array_equal(m.to_array(), 2.0 * a + 1.0)


def test_matrix_products_match_numpy() -> None:
    # Arrange
    a, b = _general(6, 7), _general(6, 8)
    x = np.arange(6.0)
    m = linalg.Matrix(a)
    # Act
    product = m.multiply(b)
    vector = m.multiply(x)
    gram = m.transpose_multiply(m)
    # Assert
    assert isinstance(product, linalg.Matrix) and isinstance(vector, linalg.Vector)
    _close(np.asarray(product), a @ b, 10.0)
    _close(np.asarray(vector), a @ x, 30.0)
    _close(np.asarray(gram), a.T @ a, 10.0)
    np.testing.assert_array_equal(np.asarray(gram), np.asarray(gram).T)


def test_matrix_determinant_and_inverse_match_numpy() -> None:
    # Arrange
    a = _general(5, 9)
    m = linalg.Matrix(a)
    # Act / Assert
    np.testing.assert_allclose(m.determinant(), np.linalg.det(a), rtol=1e-13)
    _close(np.asarray(m.inverse()), np.linalg.inv(a), np.abs(np.linalg.inv(a)).max(), 1e-13)


def test_vector_members_match_numpy() -> None:
    # Arrange
    x, y = np.array([1.0, -2.0, 3.5]), np.array([0.5, 4.0, -1.0])
    v = linalg.Vector(x)
    a = _general(3, 10)
    # Act / Assert
    assert len(v) == 3 and v.length == 3 and v[1] == -2.0
    np.testing.assert_array_equal(v.to_array(), x)
    np.testing.assert_array_equal(np.asarray(v), x)
    np.testing.assert_array_equal(v.array, x)
    assert v.to_list() == x.tolist()
    np.testing.assert_array_equal(v.add(y).to_array(), x + y)
    np.testing.assert_array_equal(v.add(linalg.Vector(y)).to_array(), x + y)
    np.testing.assert_array_equal(v.subtract(y).to_array(), x - y)
    np.testing.assert_array_equal(v.subtract(linalg.Vector(y)).to_array(), x - y)
    np.testing.assert_array_equal(v.multiply(y).to_array(), x * y)
    np.testing.assert_array_equal(v.multiply(linalg.Vector(y)).to_array(), x * y)
    np.testing.assert_array_equal(v.multiply(2.0).to_array(), 2.0 * x)
    np.testing.assert_array_equal(v.divide(4.0).to_array(), x / 4.0)
    _close(v.multiply(linalg.Matrix(a)).to_array(), a @ x, 10.0)
    np.testing.assert_allclose(v.sum(), x.sum(), rtol=1e-15)
    np.testing.assert_allclose(v.norm(), np.linalg.norm(x), rtol=1e-15)
    np.testing.assert_allclose(v.norm_squared(), x @ x, rtol=1e-15)
    np.testing.assert_allclose(linalg.Vector.distance(x, y), np.linalg.norm(x - y), rtol=1e-15)
    np.testing.assert_allclose(linalg.Vector.dot_product(x, y), x @ y, rtol=1e-15)
    np.testing.assert_allclose(linalg.Vector.project(x, y).to_array(), (x @ y) / (y @ y) * y, rtol=1e-15)


def test_vector_construction_mutation_and_header() -> None:
    # Arrange
    v = linalg.Vector(4, 2.5)
    w = linalg.Vector(3)
    # Act
    copy = v.clone()
    v[0] = 1.0
    w.copy_from([7.0, 8.0, 9.0])
    v.header = "flows"
    # Assert
    np.testing.assert_array_equal(copy.to_array(), [2.5] * 4)
    np.testing.assert_array_equal(v.to_array(), [1.0, 2.5, 2.5, 2.5])
    np.testing.assert_array_equal(w.to_array(), [7.0, 8.0, 9.0])
    assert v.header == "flows"
    w.clear()
    np.testing.assert_array_equal(w.to_array(), np.zeros(3))
    assert linalg.Vector.from_net(v.net).net is v.net


def test_cholesky_matches_numpy() -> None:
    # Arrange
    s = _spd(6, 11)
    b = np.arange(1.0, 7.0)
    # Act
    chol = linalg.CholeskyDecomposition(linalg.Matrix(s))
    # Assert
    _close(chol.l, np.linalg.cholesky(s), np.sqrt(np.abs(s).max()))
    np.testing.assert_array_equal(chol.a, s)
    assert chol.is_positive_definite
    assert chol.relative_tolerance == linalg.CholeskyDecomposition.default_relative_tolerance(6)
    np.testing.assert_allclose(chol.log_determinant(), np.linalg.slogdet(s)[1], rtol=1e-14)
    np.testing.assert_allclose(chol.determinant(), np.linalg.det(s), rtol=1e-13)
    _close(chol.solve(b), np.linalg.solve(s, b), np.abs(np.linalg.solve(s, b)).max())
    y = chol.forward(b)
    _close(y, np.linalg.solve(np.linalg.cholesky(s), b), np.abs(y).max())
    _close(chol.backward(y), np.linalg.solve(s, b), np.abs(np.linalg.solve(s, b)).max())
    _close(chol.inverse(), np.linalg.inv(s), np.abs(np.linalg.inv(s)).max())


def test_cholesky_refuses_an_indefinite_matrix_and_try_decompose_adds_jitter() -> None:
    # Arrange
    singular = np.array([[1.0, 1.0], [1.0, 1.0]])
    # Act
    ok, decomposition, jitter = linalg.CholeskyDecomposition.try_decompose(singular)
    plain_ok, _, plain_jitter = linalg.CholeskyDecomposition.try_decompose(np.eye(2))
    failed, none, nan_jitter = linalg.CholeskyDecomposition.try_decompose([[1.0, 2.0], [2.0, 1.0]], 1e-10, 1e-8)
    # Assert
    with pytest.raises(RuntimeError, match="not positive-definite"):
        linalg.CholeskyDecomposition([[1.0, 2.0], [2.0, 1.0]])
    assert ok and jitter > 0.0
    np.testing.assert_allclose(decomposition.a, singular + jitter * np.eye(2), rtol=0, atol=0)
    assert plain_ok and plain_jitter == 0.0
    assert not failed and none is None and np.isnan(nan_jitter)


def test_lu_matches_numpy() -> None:
    # Arrange
    a = _general(6, 12)
    b = np.arange(1.0, 7.0)
    rhs = _general(6, 13)[:, :3]
    # Act
    lu = linalg.LUDecomposition(a)
    # Assert
    scale = np.abs(np.linalg.inv(a)).max() * 10
    _close(lu.solve(b), np.linalg.solve(a, b), scale)
    _close(lu.solve(rhs), np.linalg.solve(a, rhs), scale)
    _close(lu.solve(linalg.Vector(b)), np.linalg.solve(a, b), scale)
    _close(lu.inverse(), np.linalg.inv(a), scale)
    np.testing.assert_allclose(lu.determinant(), np.linalg.det(a), rtol=1e-13)
    np.testing.assert_array_equal(lu.a, a)
    assert lu.lu.shape == (6, 6)
    with pytest.raises(ValueError):
        linalg.LUDecomposition(np.ones((2, 3)))


def test_qr_matches_numpy() -> None:
    # Arrange
    a = np.random.default_rng(14).normal(size=(7, 4))
    b = np.arange(7.0)
    # Act
    qr = linalg.QRDecomposition(a)
    # Assert
    _close(qr.q @ qr.r, a, 10.0)
    _close(qr.q.T @ qr.q, np.eye(qr.q.shape[1]), 10.0)
    _close(np.abs(qr.r[:4]), np.abs(np.linalg.qr(a)[1]), 10.0)
    _close(qr.solve(b), np.linalg.lstsq(a, b, rcond=None)[0], 10.0)
    with pytest.raises(ValueError, match="rank-deficient"):
        linalg.QRDecomposition(np.ones((3, 3))).solve(np.ones(3))


def test_svd_matches_numpy() -> None:
    # Arrange
    a = np.random.default_rng(15).normal(size=(6, 4))
    deficient = np.outer(np.arange(1.0, 5.0), np.arange(1.0, 4.0))
    b = np.arange(6.0)
    # Act
    svd = linalg.SingularValueDecomposition(a)
    rank1 = linalg.SingularValueDecomposition(deficient)
    # Assert
    np.testing.assert_allclose(svd.w, np.linalg.svd(a, compute_uv=False), rtol=1e-14)
    _close(svd.u @ np.diag(svd.w) @ svd.v.T, a, 10.0)
    np.testing.assert_array_equal(svd.a, a)
    np.testing.assert_allclose(svd.inverse_condition, 1.0 / np.linalg.cond(a), rtol=1e-13)
    _close(svd.solve(b), np.linalg.lstsq(a, b, rcond=None)[0], 10.0)
    square = _general(4, 20)
    _close(linalg.SingularValueDecomposition(square).solve(np.eye(4)), np.linalg.inv(square),
           10.0 * np.abs(np.linalg.inv(square)).max())
    np.testing.assert_allclose(svd.log_determinant(), np.sum(np.log(np.linalg.svd(a, compute_uv=False))), rtol=1e-14)
    assert rank1.rank() == np.linalg.matrix_rank(deficient) == 1
    assert rank1.nullity() == 2
    assert rank1.range().shape == (4, 1) and rank1.nullspace().shape == (3, 2)
    _close(deficient @ rank1.nullspace(), np.zeros((4, 2)), 100.0)
    assert rank1.threshold > 0.0
    np.testing.assert_allclose(rank1.log_pseudo_determinant(-1.0), np.log(np.linalg.svd(deficient)[1][0]),
                               rtol=1e-14)
    assert rank1.log_pseudo_determinant() <= rank1.log_pseudo_determinant(-1.0)


def test_eigen_decomposition_matches_numpy() -> None:
    # Arrange
    s = _spd(5, 16)
    # Act
    eig = linalg.EigenValueDecomposition(s)
    # Assert
    order = np.argsort(eig.eigenvalues)
    np.testing.assert_allclose(eig.eigenvalues[order], np.linalg.eigvalsh(s), rtol=1e-14)
    v = eig.eigenvectors
    # the Jacobi sweeps stop once every off-diagonal element is below 1e-12 of the largest
    _close(v @ np.diag(eig.eigenvalues) @ v.T, s, np.abs(s).max(), 1e-12)
    _close(v.T @ v, np.eye(5), 1.0, 1e-14)
    assert eig.is_converged
    lam = eig.eigenvalues / np.abs(eig.eigenvalues).max()
    np.testing.assert_allclose(eig.effective_sample_size(), lam.sum() ** 2 / (lam**2).sum(), rtol=1e-14)
    np.testing.assert_array_equal(eig.a, s)
    with pytest.raises(ValueError):
        linalg.EigenValueDecomposition(_general(3, 17))


def test_gauss_jordan_matches_numpy_and_leaves_inputs_alone() -> None:
    # Arrange
    a = _general(4, 18)
    b = np.arange(4.0)
    m = linalg.Matrix(a)
    # Act
    inverse, x = linalg.gauss_jordan_solve(m, b)
    _, xs = linalg.gauss_jordan_solve(a, np.eye(4))
    # Assert
    scale = np.abs(np.linalg.inv(a)).max() * 10
    _close(inverse, np.linalg.inv(a), scale)
    _close(x, np.linalg.solve(a, b), scale)
    _close(xs, np.linalg.inv(a), scale)
    np.testing.assert_array_equal(m.to_array(), a)


def test_regularize_floors_and_caps_eigenvalues() -> None:
    # Arrange
    s = np.diag([1e-12, 1.0, 2.0, 1000.0])
    # Act
    out = linalg.regularize(s)
    net = namespace("Numerics.Mathematics.LinearAlgebra").MatrixRegularization.Regularize(to_net_matrix(s), 1e-6, 50.0)
    # Assert
    eigenvalues = np.linalg.eigvalsh(out)
    assert eigenvalues.min() > 1e-12 and eigenvalues.max() <= 50.0 * 1.5 * (1 + 1e-12)
    np.testing.assert_array_equal(_bits(out), _bits(np.asarray([[net[i, j] for j in range(4)] for i in range(4)])))
    np.testing.assert_array_equal(linalg.regularize(linalg.Matrix(s)), out)


def test_make_symmetric_positive_definite_returns_the_cholesky_accepted_symmetric_part() -> None:
    # Arrange: an asymmetric positive-definite matrix and a symmetric indefinite one
    a = np.array([[1.0, 0.5], [0.4, 1.0]])
    indefinite = np.array([[1.0, 2.0], [2.0, 1.0]])
    # Act
    sym = linalg.make_symmetric_positive_definite(a)
    ridged = linalg.make_symmetric_positive_definite(linalg.Matrix(indefinite))
    net = namespace("Numerics.Mathematics.LinearAlgebra").MatrixRegularization.MakeSymmetricPositiveDefinite(
        linalg.Matrix(a).net)
    # Assert: the symmetric part when it factors, otherwise the symmetric part plus a diagonal ridge
    np.testing.assert_array_equal(sym, 0.5 * (a + a.T))
    np.testing.assert_array_equal(_bits(sym), _bits(np.asarray([[net[i, j] for j in range(2)] for i in range(2)])))
    np.testing.assert_array_equal(ridged - np.diag(np.diag(ridged)), indefinite - np.diag(np.diag(indefinite)))
    assert np.all(np.diag(ridged) > 1.0) and np.linalg.eigvalsh(ridged).min() > 0.0
    linalg.CholeskyDecomposition(ridged)
    with pytest.raises(ValueError, match="square"):
        linalg.make_symmetric_positive_definite(np.ones((2, 3)))


def test_decomposition_method_enum_resolves_names() -> None:
    # Act / Assert
    assert linalg.DecompositionMethod("singular_value") is linalg.DecompositionMethod.SINGULAR_VALUE
    assert linalg.DecompositionMethod("Cholesky") is linalg.DecompositionMethod.CHOLESKY


def test_large_matrices_cross_in_one_copy() -> None:
    # Arrange
    a = np.random.default_rng(19).normal(size=(500, 500))
    m = linalg.Matrix(a)
    # Act
    start = time.perf_counter()
    for _ in range(5):
        values = m.to_array()
        linalg.Matrix(a)
    elapsed = time.perf_counter() - start
    # Assert: 250,000 doubles element by element would take seconds
    np.testing.assert_array_equal(_bits(values), _bits(a))
    assert elapsed < 0.5


def test_ql_eigen_decomposition_equals_numpy_eigh() -> None:
    # Arrange
    rng = np.random.default_rng(8)
    b = rng.random((40, 40)) - 0.5
    a = b @ b.T
    expected = np.linalg.eigh(a)[0]
    # Act
    ql = linalg.EigenValueDecomposition(a, method="ql")
    jacobi = linalg.EigenValueDecomposition(a)
    # Assert
    values, vectors = ql.eigenvalues, ql.eigenvectors
    np.testing.assert_allclose(values, expected, rtol=0, atol=1e-13 * np.abs(a).max())
    assert np.all(np.diff(values) >= 0) and ql.is_converged
    np.testing.assert_allclose(a @ vectors, vectors * values, rtol=0, atol=1e-12 * np.abs(a).max())
    np.testing.assert_allclose(vectors.T @ vectors, np.eye(40), rtol=0, atol=1e-13)
    by_enum = linalg.EigenValueDecomposition(a, linalg.EigenSolverMethod.JACOBI)
    assert np.array_equal(by_enum.eigenvalues, jacobi.eigenvalues)
    with pytest.raises(ValueError):
        linalg.EigenValueDecomposition(a, method="power")
