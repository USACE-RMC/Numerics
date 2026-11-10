"""Matrices, vectors and decompositions from ``Numerics.Mathematics.LinearAlgebra``.

:class:`Matrix` and :class:`Vector` hold a .NET ``Matrix`` or ``Vector`` (``.net``) and expose
its members; matrix-valued results stay in .NET as new wrappers, and ``numpy.asarray`` copies
the values out. The decompositions (:class:`CholeskyDecomposition`, :class:`LUDecomposition`,
:class:`QRDecomposition`, :class:`SingularValueDecomposition` and
:class:`EigenValueDecomposition`) take any 2-D array or a :class:`Matrix` and return their
factors and solutions as NumPy arrays. :func:`regularize` and :func:`make_symmetric_positive_definite`
wrap ``MatrixRegularization``. Every array crosses in one bulk copy each way, bit for bit.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

from pynumerics._util import cached_overload, net_enum, overload
from pynumerics.callbacks import scalar_function
from pynumerics.enums import NetEnum
from pynumerics.errors import translate_public
from pynumerics.marshal import (
    _construct,
    _linear_algebra,
    _numerics_matrix,
    _numerics_vector,
    to_net_doubles,
    to_net_matrix,
    to_numpy_doubles,
    to_numpy_matrix,
)


def _numpy_matrix(net: Any) -> np.ndarray:
    """Copy a .NET ``Matrix`` into a new 2-D float64 array."""
    return to_numpy_matrix(net.Array)


def _numpy_vector(net: Any) -> np.ndarray:
    """Copy a .NET ``Vector`` into a new 1-D float64 array."""
    return to_numpy_doubles(net.Array)


def _is_matrix_operand(values: Any) -> bool:
    """Return whether an operand is a :class:`Matrix` or a 2-D array (otherwise a vector)."""
    return isinstance(values, Matrix) or (not isinstance(values, Vector) and np.ndim(values) == 2)


def _solve(net: Any, b: Any, *extra: float) -> np.ndarray:
    """Call ``Solve`` with a vector or a matrix right-hand side and copy the solution out."""
    if _is_matrix_operand(b):
        return _numpy_matrix(net.Solve(_numerics_matrix(b, unwrap=True), *extra))
    return _numpy_vector(net.Solve(_numerics_vector(b, unwrap=True), *extra))


class DecompositionMethod(NetEnum):
    """Matrix decomposition choices (``Numerics.Mathematics.LinearAlgebra.DecompositionMethod``)."""

    CHOLESKY = "Cholesky"
    SINGULAR_VALUE = "SingularValue"


class EigenSolverMethod(NetEnum):
    """Symmetric eigen-decomposition algorithms (``Numerics.Mathematics.LinearAlgebra.EigenSolverMethod``)."""

    JACOBI = "Jacobi"
    QL = "QL"


class Matrix:
    """A dense matrix of doubles (``Numerics.Mathematics.LinearAlgebra.Matrix``).

    Parameters
    ----------
    values : array_like or int or Vector
        A 2-D array; a 1-D array or a :class:`Vector` (a single column); or the number of rows
        of a zero matrix (square unless ``columns`` is given). With ``by_column``, a sequence of
        1-D arrays.
    columns : int, optional
        The number of columns of a zero matrix.
    by_column : bool, optional
        Build from a sequence of 1-D arrays, as columns (True) or rows (False).

    Attributes
    ----------
    net : Numerics.Mathematics.LinearAlgebra.Matrix
        The wrapped .NET matrix.
    """

    def __init__(self, values: Any, columns: int | None = None, *, by_column: bool | None = None) -> None:
        la = _linear_algebra()
        if isinstance(values, Vector):
            self.net = la.Matrix(values.net)
        elif isinstance(values, (int, np.integer)) and not isinstance(values, bool):
            self.net = la.Matrix(int(values)) if columns is None else la.Matrix(int(values), int(columns))
        elif by_column is not None:
            import System  # noqa: PLC0415 - started by _linear_algebra()
            from System.Collections.Generic import List  # noqa: PLC0415

            arrays = List[System.Array[System.Double]]()
            for item in values:
                arrays.Add(to_net_doubles(item))
            self.net = la.Matrix(arrays, bool(by_column))
        elif np.ndim(values) == 1:
            self.net = la.Matrix(to_net_doubles(values))
        else:
            self.net = _construct("Matrix", to_net_matrix(values))

    @classmethod
    def from_net(cls, net: Any) -> Matrix:
        """Wrap an existing .NET matrix without copying it.

        Parameters
        ----------
        net : Numerics.Mathematics.LinearAlgebra.Matrix
            The .NET matrix.

        Returns
        -------
        Matrix
            The wrapper.
        """
        matrix = cls.__new__(cls)
        matrix.net = net
        return matrix

    @classmethod
    def from_xml(cls, xml: str) -> Matrix:
        """Read a matrix written by :meth:`to_xml` (the ``Matrix(XElement)`` constructor).

        Parameters
        ----------
        xml : str
            The ``<Matrix>`` element text.

        Returns
        -------
        Matrix
            The matrix.
        """
        from System.Xml.Linq import XElement  # noqa: PLC0415

        _linear_algebra()
        return cls.from_net(_linear_algebra().Matrix(XElement.Parse(str(xml))))

    def __array__(self, dtype: Any = None, copy: Any = None) -> np.ndarray:
        values = _numpy_matrix(self.net)
        return values if dtype is None else values.astype(dtype)

    def __getitem__(self, index: tuple[int, int]) -> float:
        row, column = index
        return float(self.net[int(row), int(column)])

    def __setitem__(self, index: tuple[int, int], value: float) -> None:
        row, column = index
        self.net[int(row), int(column)] = float(value)

    def __repr__(self) -> str:
        return f"Matrix({self.to_array()!r})"

    @property
    def number_of_rows(self) -> int:
        """The number of rows (``NumberOfRows``)."""
        return int(self.net.NumberOfRows)

    @property
    def number_of_columns(self) -> int:
        """The number of columns (``NumberOfColumns``)."""
        return int(self.net.NumberOfColumns)

    @property
    def shape(self) -> tuple[int, int]:
        """``(number_of_rows, number_of_columns)``."""
        return self.number_of_rows, self.number_of_columns

    @property
    def header(self) -> list[str] | None:
        """The column header text (``Header``)."""
        header = self.net.Header
        return None if header is None else [str(h) for h in header]

    @header.setter
    def header(self, value: Sequence[str] | None) -> None:
        import System  # noqa: PLC0415 - the matrix exists, so the runtime is running

        self.net.Header = None if value is None else System.Array[System.String]([str(v) for v in value])

    @property
    def is_symmetric(self) -> bool:
        """Whether every pair agrees to 1e-12 times the largest magnitude (``IsSymmetric``)."""
        return bool(self.net.IsSymmetric)

    @property
    def is_square(self) -> bool:
        """Whether the matrix is square (``IsSquare``)."""
        return bool(self.net.IsSquare)

    @property
    def array(self) -> np.ndarray:
        """A copy of the underlying ``double[,]`` (``Array``)."""
        return _numpy_matrix(self.net)

    def to_array(self) -> np.ndarray:
        """Return the values as a new 2-D array (``ToArray``).

        Returns
        -------
        numpy.ndarray
            The values.
        """
        return to_numpy_matrix(self.net.ToArray())

    def to_xml(self) -> str:
        """Return the matrix as XML text (``ToXElement``).

        Returns
        -------
        str
            The ``<Matrix>`` element, values written with round-trip precision.
        """
        return str(self.net.ToXElement().ToString())

    def clone(self) -> Matrix:
        """Return a copy (``Clone``).

        Returns
        -------
        Matrix
            The copy.
        """
        return Matrix.from_net(self.net.Clone())

    def row(self, index: int) -> np.ndarray:
        """Return a copy of a row (``Row``).

        Parameters
        ----------
        index : int
            The zero-based row index.

        Returns
        -------
        numpy.ndarray
            The row.
        """
        return to_numpy_doubles(self.net.Row(int(index)))

    def column(self, index: int) -> np.ndarray:
        """Return a copy of a column (``Column``).

        Parameters
        ----------
        index : int
            The zero-based column index.

        Returns
        -------
        numpy.ndarray
            The column.
        """
        return to_numpy_doubles(self.net.Column(int(index)))

    def diagonal(self) -> np.ndarray:
        """Return the min(rows, columns) diagonal elements (``Diagonal()``).

        Returns
        -------
        numpy.ndarray
            The diagonal.
        """
        return to_numpy_doubles(self.net.Diagonal())

    @staticmethod
    def diagonal_matrix(values: Any) -> Matrix:
        """Return a diagonal matrix (the static ``Matrix.Diagonal``).

        Parameters
        ----------
        values : Matrix or Vector or array_like
            A square matrix (its diagonal is kept) or a 1-D array of diagonal entries.

        Returns
        -------
        Matrix
            The diagonal matrix.
        """
        la = _linear_algebra()
        if _is_matrix_operand(values):
            return Matrix.from_net(la.Matrix.Diagonal.Overloads[la.Matrix](_numerics_matrix(values, unwrap=True)))
        return Matrix.from_net(la.Matrix.Diagonal.Overloads[la.Vector](_numerics_vector(values, unwrap=True)))

    def upper_triangle(self) -> Matrix:
        """Return the upper triangle, with zeros below the diagonal (``UpperTriangle``).

        Returns
        -------
        Matrix
            The upper-triangular matrix.
        """
        return Matrix.from_net(self.net.UpperTriangle())

    def lower_triangle(self) -> Matrix:
        """Return the lower triangle, with zeros above the diagonal (``LowerTriangle``).

        Returns
        -------
        Matrix
            The lower-triangular matrix.
        """
        return Matrix.from_net(self.net.LowerTriangle())

    def determinant(self) -> float:
        """Return the determinant, by LU decomposition (``Determinant``).

        Returns
        -------
        float
            The determinant det(A).
        """
        return float(self.net.Determinant())

    def inverse(self) -> Matrix:
        """Return the inverse, by LU decomposition (``Inverse``).

        Returns
        -------
        Matrix
            A^-1.
        """
        return Matrix.from_net(self.net.Inverse())

    def transpose(self) -> Matrix:
        """Return the transpose (``Transpose``).

        Returns
        -------
        Matrix
            A^T.
        """
        return Matrix.from_net(self.net.Transpose())

    def trace(self) -> float:
        """Return the sum of the diagonal (``Trace``).

        Returns
        -------
        float
            The trace tr(A).
        """
        return float(self.net.Trace())

    @staticmethod
    def identity(size: int) -> Matrix:
        """Return the identity matrix (``Identity``).

        Parameters
        ----------
        size : int
            The number of rows and columns.

        Returns
        -------
        Matrix
            I.
        """
        return Matrix.from_net(_linear_algebra().Matrix.Identity(int(size)))

    @staticmethod
    def outer(x: Any, y: Any) -> Matrix:
        """Return the outer product x y^T (``Outer``).

        Parameters
        ----------
        x, y : Vector or array_like
            1-D operands.

        Returns
        -------
        Matrix
            The ``len(x)`` by ``len(y)`` product.
        """
        return Matrix.from_net(_linear_algebra().Matrix.Outer(_numerics_vector(x, unwrap=True),
                                                              _numerics_vector(y, unwrap=True)))

    def column_means(self) -> np.ndarray:
        """Return the mean of each column (``ColumnMeans``).

        Returns
        -------
        numpy.ndarray
            One mean per column.
        """
        return _numpy_vector(self.net.ColumnMeans())

    def apply(self, func: Callable[[float], float] | str) -> None:
        """Replace every element ``a`` by ``func(a)`` in place (``Apply``).

        Parameters
        ----------
        func : callable or str
            A Python callable (one crossing per element) or C# lambda source compiled with
            Roslyn (no crossing per element), as :func:`pynumerics.callbacks.scalar_function`
            accepts.
        """
        self.net.Apply(scalar_function(func))

    def sqrt(self) -> None:
        """Take the square root of every element in place (``Sqrt``)."""
        self.net.Sqrt()

    def sqr(self) -> None:
        """Square every element in place (``Sqr``)."""
        self.net.Sqr()

    def log(self) -> None:
        """Take the natural logarithm of every element in place (``Log``)."""
        self.net.Log()

    def exp(self) -> None:
        """Exponentiate every element in place (``Exp``)."""
        self.net.Exp()

    def sum(self) -> float:
        """Return the sum of every element (``Sum``).

        Returns
        -------
        float
            The sum.
        """
        return float(self.net.Sum())

    def multiply(self, other: Any) -> Matrix | Vector:
        """Return the product with a matrix, a vector or a scalar (``Multiply``).

        Parameters
        ----------
        other : Matrix or Vector or array_like or float
            A matrix (or 2-D array), a vector (or 1-D array) or a scalar.

        Returns
        -------
        Matrix or Vector
            A B or A s (a :class:`Matrix`), or A x (a :class:`Vector`).
        """
        if _is_matrix_operand(other):
            member = cached_overload("Matrix.Multiply(Matrix)",
                                     lambda: _linear_algebra().Matrix.Multiply.Overloads[_linear_algebra().Matrix])
            return Matrix.from_net(member(self.net, _numerics_matrix(other, unwrap=True)))
        if isinstance(other, Vector) or np.ndim(other) == 1:
            member = cached_overload("Matrix.Multiply(Vector)",
                                     lambda: _linear_algebra().Matrix.Multiply.Overloads[_linear_algebra().Vector])
            return Vector.from_net(member(self.net, _numerics_vector(other, unwrap=True)))
        return Matrix.from_net(overload(self.net.Multiply, "Double")(float(other)))

    def transpose_multiply(self, other: Any) -> Matrix:
        """Return A^T B without forming the transpose (``TransposeMultiply``).

        Parameters
        ----------
        other : Matrix or array_like
            B, with as many rows as this matrix.

        Returns
        -------
        Matrix
            A^T B; exactly symmetric when B is A.
        """
        return Matrix.from_net(self.net.TransposeMultiply(_numerics_matrix(other, unwrap=True)))

    def divide(self, scalar: float) -> Matrix:
        """Return the matrix divided by a scalar (``Divide``).

        Parameters
        ----------
        scalar : float
            The divisor.

        Returns
        -------
        Matrix
            A / s.
        """
        return Matrix.from_net(self.net.Divide(float(scalar)))

    def add(self, other: Any) -> Matrix:
        """Return the elementwise sum (``Add``).

        Parameters
        ----------
        other : Matrix or array_like
            A matrix of the same shape.

        Returns
        -------
        Matrix
            A + B.
        """
        return Matrix.from_net(self.net.Add(_numerics_matrix(other, unwrap=True)))

    def subtract(self, other: Any) -> Matrix:
        """Return the elementwise difference (``Subtract``).

        Parameters
        ----------
        other : Matrix or array_like
            A matrix of the same shape.

        Returns
        -------
        Matrix
            A - B.
        """
        return Matrix.from_net(self.net.Subtract(_numerics_matrix(other, unwrap=True)))


class Vector:
    """A dense vector of doubles (``Numerics.Mathematics.LinearAlgebra.Vector``).

    Parameters
    ----------
    values : array_like or int
        The 1-D values, or the length of a new vector.
    fill : float, optional
        The value of every element of a new vector of length ``values``; None gives zeros.

    Attributes
    ----------
    net : Numerics.Mathematics.LinearAlgebra.Vector
        The wrapped .NET vector.
    """

    def __init__(self, values: Any, fill: float | None = None) -> None:
        la = _linear_algebra()
        if isinstance(values, (int, np.integer)) and not isinstance(values, bool):
            self.net = la.Vector(int(values)) if fill is None else la.Vector(int(values), float(fill))
        else:
            self.net = _construct("Vector", to_net_doubles(values))

    @classmethod
    def from_net(cls, net: Any) -> Vector:
        """Wrap an existing .NET vector without copying it.

        Parameters
        ----------
        net : Numerics.Mathematics.LinearAlgebra.Vector
            The .NET vector.

        Returns
        -------
        Vector
            The wrapper.
        """
        vector = cls.__new__(cls)
        vector.net = net
        return vector

    def __array__(self, dtype: Any = None, copy: Any = None) -> np.ndarray:
        values = _numpy_vector(self.net)
        return values if dtype is None else values.astype(dtype)

    def __len__(self) -> int:
        return int(self.net.Length)

    def __getitem__(self, index: int) -> float:
        return float(self.net[int(index)])

    def __setitem__(self, index: int, value: float) -> None:
        self.net[int(index)] = float(value)

    def __repr__(self) -> str:
        return f"Vector({self.to_array()!r})"

    @property
    def length(self) -> int:
        """The number of elements (``Length``)."""
        return int(self.net.Length)

    @property
    def header(self) -> str | None:
        """The header text (``Header``)."""
        header = self.net.Header
        return None if header is None else str(header)

    @header.setter
    def header(self, value: str | None) -> None:
        self.net.Header = None if value is None else str(value)

    @property
    def array(self) -> np.ndarray:
        """A copy of the underlying ``double[]`` (``Array``)."""
        return _numpy_vector(self.net)

    def to_array(self) -> np.ndarray:
        """Return the values as a new 1-D array (``ToArray``).

        Returns
        -------
        numpy.ndarray
            The values.
        """
        return to_numpy_doubles(self.net.ToArray())

    def to_list(self) -> list[float]:
        """Return the values as a Python list (``ToList``).

        Returns
        -------
        list of float
            The values.
        """
        return to_numpy_doubles(self.net.ToList()).tolist()

    def clone(self) -> Vector:
        """Return a copy (``Clone``).

        Returns
        -------
        Vector
            The copy.
        """
        return Vector.from_net(self.net.Clone())

    def clear(self) -> None:
        """Set every element to zero in place (``Clear``)."""
        self.net.Clear()

    def copy_from(self, other: Any) -> None:
        """Overwrite the elements with those of another vector of the same length (``CopyFrom``).

        Parameters
        ----------
        other : Vector or array_like
            The source.
        """
        self.net.CopyFrom(_numerics_vector(other, unwrap=True))

    def sum(self) -> float:
        """Return the sum of the elements (``Sum``).

        Returns
        -------
        float
            The sum.
        """
        return float(self.net.Sum())

    def norm(self) -> float:
        """Return the Euclidean norm (``Norm``).

        Returns
        -------
        float
            ||v||.
        """
        return float(self.net.Norm())

    def norm_squared(self) -> float:
        """Return the squared Euclidean norm (``NormSquared``).

        Returns
        -------
        float
            ||v||**2.
        """
        return float(self.net.NormSquared())

    @staticmethod
    def distance(a: Any, b: Any) -> float:
        """Return the Euclidean distance ||a - b|| (``Distance``).

        Parameters
        ----------
        a, b : Vector or array_like
            Equal-length vectors.

        Returns
        -------
        float
            The distance.
        """
        return float(_linear_algebra().Vector.Distance(_numerics_vector(a, unwrap=True),
                                                       _numerics_vector(b, unwrap=True)))

    @staticmethod
    def dot_product(a: Any, b: Any) -> float:
        """Return the dot product a . b (``DotProduct``).

        Parameters
        ----------
        a, b : Vector or array_like
            Equal-length vectors.

        Returns
        -------
        float
            The dot product.
        """
        return float(_linear_algebra().Vector.DotProduct(_numerics_vector(a, unwrap=True),
                                                         _numerics_vector(b, unwrap=True)))

    @staticmethod
    def project(a: Any, b: Any) -> Vector:
        """Return the projection of a onto b, (a . b / b . b) b (``Project``).

        Parameters
        ----------
        a, b : Vector or array_like
            Equal-length vectors.

        Returns
        -------
        Vector
            The projection.
        """
        return Vector.from_net(_linear_algebra().Vector.Project(_numerics_vector(a, unwrap=True),
                                                                _numerics_vector(b, unwrap=True)))

    def add(self, other: Any) -> Vector:
        """Return the elementwise sum (``Add``).

        Parameters
        ----------
        other : Vector or array_like
            An equal-length vector.

        Returns
        -------
        Vector
            The sum v + w.
        """
        if isinstance(other, Vector):
            return Vector.from_net(self.net.Add.Overloads[_linear_algebra().Vector](other.net))
        return Vector.from_net(overload(self.net.Add, "Double[]")(to_net_doubles(other)))

    def subtract(self, other: Any) -> Vector:
        """Return the elementwise difference (``Subtract``).

        Parameters
        ----------
        other : Vector or array_like
            An equal-length vector.

        Returns
        -------
        Vector
            The difference v - w.
        """
        if isinstance(other, Vector):
            return Vector.from_net(self.net.Subtract.Overloads[_linear_algebra().Vector](other.net))
        return Vector.from_net(overload(self.net.Subtract, "Double[]")(to_net_doubles(other)))

    def multiply(self, other: Any) -> Vector:
        """Return a matrix times this vector, the elementwise product, or the scaled vector (``Multiply``).

        Parameters
        ----------
        other : Matrix or Vector or array_like or float
            A matrix M (or 2-D array) whose column count is the vector length, for M v; a vector
            (or 1-D array) for the elementwise product; or a scalar.

        Returns
        -------
        Vector
            The product.
        """
        la = _linear_algebra()
        if _is_matrix_operand(other):
            return Vector.from_net(self.net.Multiply.Overloads[la.Matrix](_numerics_matrix(other, unwrap=True)))
        if isinstance(other, Vector):
            return Vector.from_net(self.net.Multiply.Overloads[la.Vector](other.net))
        if np.ndim(other) == 1:
            return Vector.from_net(overload(self.net.Multiply, "Double[]")(to_net_doubles(other)))
        return Vector.from_net(overload(self.net.Multiply, "Double")(float(other)))

    def divide(self, scalar: float) -> Vector:
        """Return the vector divided by a scalar (``Divide``).

        Parameters
        ----------
        scalar : float
            The divisor.

        Returns
        -------
        Vector
            The quotient v / s.
        """
        return Vector.from_net(self.net.Divide(float(scalar)))


class CholeskyDecomposition:
    """The Cholesky factorization A = L L^T of a symmetric positive-definite matrix (``CholeskyDecomposition``).

    Only the upper triangle and the diagonal of A are read. A pivot is rejected when it is at or
    below ``relative_tolerance`` times the corresponding diagonal entry.

    Parameters
    ----------
    matrix : Matrix or array_like
        The symmetric positive-definite matrix.
    relative_tolerance : float, optional
        The scale-relative pivot tolerance in [0, 1); None uses
        :meth:`default_relative_tolerance` of the dimension.

    Attributes
    ----------
    net : Numerics.Mathematics.LinearAlgebra.CholeskyDecomposition
        The wrapped .NET decomposition.

    Raises
    ------
    RuntimeError
        If the matrix is not positive definite (raised by Numerics).
    ValueError
        If the matrix is not square or the tolerance is outside [0, 1) (raised by Numerics).
    """

    def __init__(self, matrix: Any, relative_tolerance: float | None = None) -> None:
        la = _linear_algebra()
        self.net = (la.CholeskyDecomposition(_numerics_matrix(matrix, unwrap=True)) if relative_tolerance is None
                    else la.CholeskyDecomposition(_numerics_matrix(matrix, unwrap=True), float(relative_tolerance)))

    @classmethod
    def try_decompose(cls, matrix: Any, minimum_relative_jitter: float | None = None,
                      maximum_relative_jitter: float | None = None) -> tuple[bool, CholeskyDecomposition | None,
                                                                             float]:
        """Factor A, retrying on A + jitter I with an escalating jitter (``TryCholesky``).

        The jitters are s * minimum * 10**k up to s * maximum, with s the mean diagonal of A
        (or 1); the defaults are the Numerics ones, 1e-10 and 1e-2.

        Parameters
        ----------
        matrix : Matrix or array_like
            The symmetric matrix.
        minimum_relative_jitter, maximum_relative_jitter : float, optional
            The jitter range relative to s; give both or neither.

        Returns
        -------
        success : bool
            Whether A or a jittered A was factored.
        decomposition : CholeskyDecomposition or None
            The factorization of A + jitter I on success.
        jitter : float
            The absolute jitter added: 0 when A factors as given, NaN on failure.
        """
        member = _linear_algebra().CholeskyDecomposition.TryCholesky
        if minimum_relative_jitter is None and maximum_relative_jitter is None:
            ok, net, jitter = member(_numerics_matrix(matrix, unwrap=True), None, 0.0)
        else:
            ok, net, jitter = member(_numerics_matrix(matrix, unwrap=True), float(minimum_relative_jitter),
                                     float(maximum_relative_jitter), None, 0.0)
        if not ok:
            return False, None, float(jitter)
        decomposition = cls.__new__(cls)
        decomposition.net = net
        return True, decomposition, float(jitter)

    @staticmethod
    def default_relative_tolerance(dimension: int) -> float:
        """Return the default pivot tolerance for a dimension (``DefaultRelativeTolerance``).

        Parameters
        ----------
        dimension : int
            The matrix order n.

        Returns
        -------
        float
            The tolerance (2 n eps).
        """
        return float(_linear_algebra().CholeskyDecomposition.DefaultRelativeTolerance(int(dimension)))

    @property
    def relative_tolerance(self) -> float:
        """The pivot tolerance used (``RelativeTolerance``)."""
        return float(self.net.RelativeTolerance)

    @property
    def l(self) -> np.ndarray:  # noqa: E743 - the factor's conventional name
        """The lower-triangular factor L (``L``)."""
        return _numpy_matrix(self.net.L)

    @property
    def a(self) -> np.ndarray:
        """The decomposed matrix (``A``; A + jitter I after :meth:`try_decompose`)."""
        return _numpy_matrix(self.net.A)

    @property
    def is_positive_definite(self) -> bool:
        """Whether the factorization succeeded (``IsPositiveDefinite``)."""
        return bool(self.net.IsPositiveDefinite)

    def solve(self, b: Any) -> np.ndarray:
        """Solve A x = b (``Solve``).

        Parameters
        ----------
        b : Vector or array_like
            The right-hand side.

        Returns
        -------
        numpy.ndarray
            The solution x.
        """
        return _numpy_vector(self.net.Solve(_numerics_vector(b, unwrap=True)))

    def forward(self, b: Any) -> np.ndarray:
        """Solve L y = b by forward substitution (``Forward``).

        Parameters
        ----------
        b : Vector or array_like
            The right-hand side.

        Returns
        -------
        numpy.ndarray
            The solution y.
        """
        return _numpy_vector(self.net.Forward(_numerics_vector(b, unwrap=True)))

    def backward(self, y: Any) -> np.ndarray:
        """Solve L^T x = y by back substitution (``Backward``).

        Parameters
        ----------
        y : Vector or array_like
            The right-hand side.

        Returns
        -------
        numpy.ndarray
            The solution x.
        """
        return _numpy_vector(self.net.Backward(_numerics_vector(y, unwrap=True)))

    def inverse(self) -> np.ndarray:
        """Return A^-1 from the factor (``InverseA``).

        Returns
        -------
        numpy.ndarray
            The inverse.
        """
        return _numpy_matrix(self.net.InverseA())

    def determinant(self) -> float:
        """Return det(A) from the factor (``Determinant``).

        Returns
        -------
        float
            The determinant.
        """
        return float(self.net.Determinant())

    def log_determinant(self) -> float:
        """Return ln det(A) = 2 sum ln L[i, i] (``LogDeterminant``).

        Returns
        -------
        float
            The log determinant.
        """
        return float(self.net.LogDeterminant())


class LUDecomposition:
    """The LU decomposition with partial pivoting of a square matrix (``LUDecomposition``).

    Parameters
    ----------
    matrix : Matrix or array_like
        The square matrix; it is copied, not modified.

    Attributes
    ----------
    net : Numerics.Mathematics.LinearAlgebra.LUDecomposition
        The wrapped .NET decomposition.

    Raises
    ------
    ValueError
        If the matrix is not square (raised by Numerics).
    """

    def __init__(self, matrix: Any) -> None:
        self.net = _linear_algebra().LUDecomposition(_numerics_matrix(matrix, unwrap=True))

    @property
    def lu(self) -> np.ndarray:
        """The packed factors of the row-permuted matrix (``LU``; unit diagonal of L implied)."""
        return _numpy_matrix(self.net.LU)

    @property
    def a(self) -> np.ndarray:
        """The decomposed matrix (``A``)."""
        return _numpy_matrix(self.net.A)

    def solve(self, b: Any) -> np.ndarray:
        """Solve A x = b for a vector or a matrix of right-hand sides (``Solve``).

        Parameters
        ----------
        b : Vector or Matrix or array_like
            A 1-D right-hand side or a 2-D array of them (one per column).

        Returns
        -------
        numpy.ndarray
            The solution x, shaped like ``b``.
        """
        return _solve(self.net, b)

    def inverse(self) -> np.ndarray:
        """Return A^-1 (``InverseA``).

        Returns
        -------
        numpy.ndarray
            The inverse.
        """
        return _numpy_matrix(self.net.InverseA())

    def determinant(self) -> float:
        """Return det(A) (``Determinant``).

        Returns
        -------
        float
            The determinant.
        """
        return float(self.net.Determinant())


class QRDecomposition:
    """The Householder QR decomposition A = Q R of an m-by-n matrix (``QRDecomposition``).

    Parameters
    ----------
    matrix : Matrix or array_like
        The matrix.

    Attributes
    ----------
    net : Numerics.Mathematics.LinearAlgebra.QRDecomposition
        The wrapped .NET decomposition.
    """

    def __init__(self, matrix: Any) -> None:
        self.net = _linear_algebra().QRDecomposition(_numerics_matrix(matrix, unwrap=True))

    @property
    def q(self) -> np.ndarray:
        """The orthogonal factor Q (``Q``)."""
        return _numpy_matrix(self.net.Q)

    @property
    def r(self) -> np.ndarray:
        """The upper-triangular factor R (``RMatrix``)."""
        return _numpy_matrix(self.net.RMatrix)

    def solve(self, b: Any) -> np.ndarray:
        """Solve A x = b (least squares when m > n) for a vector or a matrix (``Solve``).

        Parameters
        ----------
        b : Vector or Matrix or array_like
            A 1-D right-hand side or a 2-D array of them.

        Returns
        -------
        numpy.ndarray
            The solution x.

        Raises
        ------
        ValueError
            If ``b`` does not match A or A is numerically rank-deficient (raised by Numerics).
        """
        return _solve(self.net, b)


class SingularValueDecomposition:
    """The singular value decomposition A = U diag(W) V^T (``SingularValueDecomposition``).

    A ``threshold`` argument below zero, the Numerics default, means the roundoff-based
    :attr:`threshold`; singular values below the threshold count as zero.

    Parameters
    ----------
    matrix : Matrix or array_like
        The m-by-n matrix.

    Attributes
    ----------
    net : Numerics.Mathematics.LinearAlgebra.SingularValueDecomposition
        The wrapped .NET decomposition.
    """

    def __init__(self, matrix: Any) -> None:
        self.net = _linear_algebra().SingularValueDecomposition(_numerics_matrix(matrix, unwrap=True))

    @property
    def a(self) -> np.ndarray:
        """The decomposed matrix (``A``)."""
        return _numpy_matrix(self.net.A)

    @property
    def u(self) -> np.ndarray:
        """The m-by-n column-orthogonal factor U (``U``)."""
        return _numpy_matrix(self.net.U)

    @property
    def v(self) -> np.ndarray:
        """The n-by-n orthogonal factor V, so A = U diag(w) V^T (``V``)."""
        return _numpy_matrix(self.net.V)

    @property
    def w(self) -> np.ndarray:
        """The singular values, in decreasing order (``W``)."""
        return _numpy_vector(self.net.W)

    @property
    def threshold(self) -> float:
        """The default zero threshold 0.5 sqrt(m + n + 1) w[0] eps (``Threshold``)."""
        return float(self.net.Threshold)

    @property
    def inverse_condition(self) -> float:
        """The reciprocal condition number w[-1] / w[0] (``InverseCondition``)."""
        return float(self.net.InverseCondition)

    def rank(self, threshold: float = -1.0) -> int:
        """Return the number of singular values above the threshold (``Rank``).

        Parameters
        ----------
        threshold : float, optional
            The zero threshold; negative (default) uses :attr:`threshold`.

        Returns
        -------
        int
            The rank.
        """
        return int(self.net.Rank(float(threshold)))

    def nullity(self, threshold: float = -1.0) -> int:
        """Return the number of singular values at or below the threshold (``Nullity``).

        Parameters
        ----------
        threshold : float, optional
            The zero threshold; negative (default) uses :attr:`threshold`.

        Returns
        -------
        int
            The nullity.
        """
        return int(self.net.Nullity(float(threshold)))

    def range(self, threshold: float = -1.0) -> np.ndarray:
        """Return an orthonormal basis of the range of A, one vector per column (``Range``).

        Parameters
        ----------
        threshold : float, optional
            The zero threshold; negative (default) uses :attr:`threshold`.

        Returns
        -------
        numpy.ndarray
            The m-by-rank basis.
        """
        return _numpy_matrix(self.net.Range(float(threshold)))

    def nullspace(self, threshold: float = -1.0) -> np.ndarray:
        """Return an orthonormal basis of the null space of A, one vector per column (``Nullspace``).

        Parameters
        ----------
        threshold : float, optional
            The zero threshold; negative (default) uses :attr:`threshold`.

        Returns
        -------
        numpy.ndarray
            The n-by-nullity basis.
        """
        return _numpy_matrix(self.net.Nullspace(float(threshold)))

    def solve(self, b: Any, threshold: float = -1.0) -> np.ndarray:
        """Return the pseudo-inverse solution of A x = b for a vector or a matrix (``Solve``).

        The matrix overload of Numerics requires ``b`` to have as many rows as A has columns and
        reads its first m rows, so a matrix right-hand side is meaningful for a square A only.

        Parameters
        ----------
        b : Vector or Matrix or array_like
            A 1-D right-hand side of length m, or a 2-D array of them (square A).
        threshold : float, optional
            The zero threshold; negative (default) uses :attr:`threshold`.

        Returns
        -------
        numpy.ndarray
            The solution x.
        """
        return _solve(self.net, b, float(threshold))

    def log_determinant(self) -> float:
        """Return the sum of the logs of the singular values, ln abs(det A) (``LogDeterminant``).

        Returns
        -------
        float
            The log absolute determinant.
        """
        return float(self.net.LogDeterminant())

    def log_pseudo_determinant(self, threshold: float | None = None) -> float:
        """Return the sum of the logs of the nonzero singular values (``LogPseudoDeterminant``).

        Parameters
        ----------
        threshold : float, optional
            None (the no-argument overload) drops only singular values exactly equal to zero; a
            value drops those below it, and a negative value uses :attr:`threshold`, matching
            :meth:`rank`.

        Returns
        -------
        float
            The log pseudo-determinant.
        """
        if threshold is None:
            return float(self.net.LogPseudoDeterminant())
        return float(self.net.LogPseudoDeterminant(float(threshold)))


class EigenValueDecomposition:
    """The eigen-decomposition A = V diag(lambda) V^T of a symmetric matrix (``EigenValueDecomposition``).

    ``method="jacobi"`` (the default) uses cyclic Jacobi rotations and leaves the eigenvalues in
    diagonal order, not sorted. ``method="ql"`` uses Householder tridiagonalization and implicit QL
    iterations, returns the eigenvalues in ascending order as ``numpy.linalg.eigh`` does, and is
    many times faster for large matrices.

    Parameters
    ----------
    matrix : Matrix or array_like
        The symmetric matrix.
    method : EigenSolverMethod or str, optional
        ``"jacobi"`` (default) or ``"ql"``.

    Attributes
    ----------
    net : Numerics.Mathematics.LinearAlgebra.EigenValueDecomposition
        The wrapped .NET decomposition.

    Raises
    ------
    ValueError
        If the matrix is not square, not symmetric or not finite (raised by Numerics), or
        ``method`` names no method.
    """

    def __init__(self, matrix: Any, method: EigenSolverMethod | str = EigenSolverMethod.JACOBI) -> None:
        la = _linear_algebra()
        kind = EigenSolverMethod(method)
        a = _numerics_matrix(matrix, unwrap=True)
        if kind is EigenSolverMethod.JACOBI:
            self.net = la.EigenValueDecomposition(a)
            return
        self.net = la.EigenValueDecomposition(a, net_enum(la.EigenSolverMethod, kind))

    @property
    def a(self) -> np.ndarray:
        """The decomposed matrix (``A``)."""
        return _numpy_matrix(self.net.A)

    @property
    def eigenvalues(self) -> np.ndarray:
        """The eigenvalues (``EigenValues``)."""
        return _numpy_vector(self.net.EigenValues)

    @property
    def eigenvectors(self) -> np.ndarray:
        """The orthonormal eigenvectors, one per column (``EigenVectors``)."""
        return _numpy_matrix(self.net.EigenVectors)

    @property
    def is_converged(self) -> bool:
        """Whether every off-diagonal element fell to 1e-12 of the largest (``IsConverged``)."""
        return bool(self.net.IsConverged)

    def effective_sample_size(self) -> float:
        """Return Dutilleul's effective sample size (sum lambda)**2 / sum lambda**2 (``EffectiveSampleSize``).

        Returns
        -------
        float
            The effective sample size, or 0 for an all-zero spectrum.
        """
        return float(self.net.EffectiveSampleSize())


def gauss_jordan_solve(a: Any, b: Any) -> tuple[np.ndarray, np.ndarray]:
    """Solve A X = B by Gauss-Jordan elimination with full pivoting (``GaussJordanElimination.Solve``).

    Parameters
    ----------
    a : Matrix or array_like
        The square matrix A; a :class:`Matrix` is not modified (a copy is eliminated).
    b : Matrix or array_like
        The right-hand sides, one per column (2-D) or a single 1-D one.

    Returns
    -------
    inverse : numpy.ndarray
        A^-1.
    solution : numpy.ndarray
        X, shaped like ``b``.

    Raises
    ------
    ValueError
        If A is not square, B has a different number of rows, or A is not finite (raised by Numerics).
    """
    rhs = np.asarray(b, dtype=np.float64)
    column = rhs.ndim == 1
    net_a = _numerics_matrix(a, unwrap=True).Clone()
    net_b = _numerics_matrix(rhs.reshape(-1, 1) if column else rhs) if not isinstance(b, Matrix) else b.net.Clone()
    inverse, solution = _linear_algebra().GaussJordanElimination.Solve(net_a, net_b)
    out = _numpy_matrix(solution)
    return _numpy_matrix(inverse), out.ravel() if column else out


def regularize(matrix: Any, eps: float = 1e-6, cap_multiplier: float = 50.0) -> np.ndarray:
    """Return a positive semi-definite version of a symmetric matrix (``MatrixRegularization.Regularize``).

    The matrix is symmetrized and eigen-decomposed, its eigenvalues are floored at
    eps * trace / p (eps when the trace is not positive) and capped at cap_multiplier times
    their median, and it is rebuilt as Q diag(lambda) Q^T.

    Parameters
    ----------
    matrix : array_like
        The (approximately) symmetric matrix.
    eps : float, optional
        The relative eigenvalue floor; default 1e-6.
    cap_multiplier : float, optional
        The eigenvalue cap as a multiple of the median; default 50.

    Returns
    -------
    numpy.ndarray
        The regularized matrix.
    """
    values = np.asarray(matrix.to_array() if isinstance(matrix, Matrix) else matrix, dtype=np.float64)
    return to_numpy_matrix(_linear_algebra().MatrixRegularization.Regularize(to_net_matrix(values), float(eps),
                                                                             float(cap_multiplier)))


def make_symmetric_positive_definite(matrix: Any) -> np.ndarray:
    """Return the nearest Cholesky-accepted symmetric matrix (``MatrixRegularization.MakeSymmetricPositiveDefinite``).

    The symmetric part when its Cholesky factor exists, otherwise the symmetric part plus
    the smallest trace-scaled decade ridge whose factor exists.

    Parameters
    ----------
    matrix : Matrix or array_like
        A finite square matrix.

    Returns
    -------
    numpy.ndarray
        The regularized matrix.

    Raises
    ------
    ValueError
        If the matrix is not square or holds non-finite entries (raised by Numerics).
    RuntimeError
        If no ridge on the ladder yields a factor (raised by Numerics).
    """
    regularized = _linear_algebra().MatrixRegularization.MakeSymmetricPositiveDefinite(
        _numerics_matrix(matrix, unwrap=True))
    return _numpy_matrix(regularized)


__all__ = [
    "CholeskyDecomposition", "DecompositionMethod", "EigenSolverMethod", "EigenValueDecomposition", "LUDecomposition",
    "Matrix", "QRDecomposition", "SingularValueDecomposition", "Vector", "gauss_jordan_solve",
    "make_symmetric_positive_definite", "regularize",
]

translate_public(globals())
