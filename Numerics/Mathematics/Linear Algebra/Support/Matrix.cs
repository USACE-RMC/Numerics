using System;
using System.Collections.Generic;
using System.Globalization;
using System.Xml.Linq;

namespace Numerics.Mathematics.LinearAlgebra
{

    /// <summary>
    /// A simple class for performing Matrix operations.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// This class is contains the basis for all Matrix operations used by Numerics.
    /// </para>
    /// </remarks>
    [Serializable]
    public class Matrix
    {

        #region Construction

        /// <summary>
        /// Construct a new matrix with specified number of rows and columns.
        /// </summary>
        /// <param name="numberOfRows">The number of rows.</param>
        /// <param name="numberOfColumns">The number of columns.</param>
        public Matrix(int numberOfRows, int numberOfColumns)
        {
            _matrix = new double[numberOfRows, numberOfColumns];
        }

        /// <summary>
        /// Constructs a new square matrix.
        /// </summary>
        /// <param name="size">The number of rows or columns (rows = columns).</param>
        public Matrix(int size)
        {
            _matrix = new double[size, size];
        }

        /// <summary>
        /// Construct a new matrix based on an initial array.
        /// </summary>
        /// <param name="initialArray">Initializing array.</param>
        public Matrix(double[,] initialArray)
        {
            if (initialArray == null) throw new ArgumentNullException(nameof(initialArray));
            int rows = initialArray.GetLength(0);
            int cols = initialArray.GetLength(1);
            _matrix = new double[rows, cols];
            // Deep copy of the contiguous row-major storage.
            Buffer.BlockCopy(initialArray, 0, _matrix, 0, rows * cols * sizeof(double));
        }

        /// <summary>
        /// Constructs a new matrix based on a single column array.
        /// </summary>
        /// <param name="initialArray">Initializing array.</param>
        public Matrix(double[] initialArray)
        {
            _matrix = new double[initialArray.Length, 1];
            Buffer.BlockCopy(initialArray, 0, _matrix, 0, initialArray.Length * sizeof(double));
        }

        /// <summary>
        /// Constructs a new matrix based on a single column array.
        /// </summary>
        /// <param name="initialArray">Initializing array.</param>
        public Matrix(float[] initialArray)
        {
            _matrix = new double[initialArray.Length, 1];
            for (int i = 0; i < _matrix.Length; i++)
                _matrix[i, 0] = initialArray[i];
        }

        /// <summary>
        /// Constructs a new matrix based on a vector.
        /// </summary>
        /// <param name="initialVector">Initializing vector.</param>
        public Matrix(Vector initialVector)
        {
            _matrix = new double[initialVector.Length, 1];
            for (int i = 0; i < _matrix.Length; i++) _matrix[i, 0] = initialVector[i];
        }

        /// <summary>
        /// Constructs a new matrix based on a list of single column arrays.
        /// </summary>
        /// <param name="listOfArrays">List of initializing arrays.</param>
        /// <param name="byColumn">Determines if the list should be added as columns or rows. Default = columns.</param>
        public Matrix(List<double[]> listOfArrays, bool byColumn = true)
        {
            if (listOfArrays == null || listOfArrays.Count == 0)
                throw new ArgumentException("List of arrays must not be empty.");

            if (byColumn)
            {
                int rowLength = listOfArrays[0].Length;
                foreach (var array in listOfArrays)
                {
                    if (array.Length != rowLength)
                        throw new ArgumentException("All arrays must have the same length.");
                }

                _matrix = new double[rowLength, listOfArrays.Count];
                for (int i = 0; i < rowLength; i++)
                    for (int j = 0; j < listOfArrays.Count; j++)
                        _matrix[i, j] = listOfArrays[j][i];
            }
            else
            {
                int colLength = listOfArrays[0].Length;
                foreach (var array in listOfArrays)
                {
                    if (array.Length != colLength)
                        throw new ArgumentException("All arrays must have the same length.");
                }

                _matrix = new double[listOfArrays.Count, colLength];
                for (int i = 0; i < listOfArrays.Count; i++)
                    for (int j = 0; j < colLength; j++)
                        _matrix[i, j] = listOfArrays[i][j];
            }

        }

        /// <summary>
        /// Construct a new matrix from XElement.
        /// </summary>
        /// <param name="xElement">The XElement to deserialize.</param>
        public Matrix(XElement xElement)
        {
            if (xElement.Name != nameof(Matrix))
            {
                return;
            }
                
          
            int nrow = 0, ncol = 0;
            var nrowAttr = xElement.Attribute(nameof(NumberOfRows));
            if (nrowAttr != null) int.TryParse(nrowAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out nrow);
            var ncolAttr = xElement.Attribute(nameof(NumberOfColumns));
            if (ncolAttr != null) int.TryParse(ncolAttr.Value, NumberStyles.Any, CultureInfo.InvariantCulture, out ncol);

            _matrix = new double[nrow, ncol];
            int rowCount = 0;
            foreach (var row in xElement.Elements("row"))
            {
                if (rowCount >= nrow)
                    break;

                var parts = row.Value.Split('|');
                int maxCols = Math.Min(parts.Length, ncol);

                for (int j = 0; j < maxCols; j++)
                {
                    if (double.TryParse(parts[j], NumberStyles.Any, CultureInfo.InvariantCulture, out var val))
                    {
                        _matrix[rowCount, j] = val;
                    }
                    else
                    {
                        _matrix[rowCount, j] = double.NaN;
                    }
                }

                rowCount++;
            }          
        }

        #endregion

        #region Members

        private double[,] _matrix = new double[0, 0];

        /// <summary>
        /// Below this number of result columns the products use plain scalar loops.
        /// </summary>
        private const int VectorizedWidth = 8;

        /// <summary>
        /// Gets the number of rows.
        /// </summary>
        public int NumberOfRows
        {
            get { return _matrix.GetLength(0); }
        }

        /// <summary>
        /// Gets the number of columns.
        /// </summary>
        public int NumberOfColumns
        {
            get { return _matrix.GetLength(1); }
        }

        /// <summary>
        /// Get and sets the element at the specific row and column index.
        /// </summary>
        /// <param name="rowIndex">The zero-based row index of the element to get or set.</param>
        /// <param name="columnIndex">The zero-based column index of the element to get or set.</param>
        public double this[int rowIndex, int columnIndex]
        {
            get { return _matrix[rowIndex, columnIndex]; }
            set { _matrix[rowIndex, columnIndex] = value; }
        }

        /// <summary>
        /// The matrix column header text. 
        /// </summary>
        public string[] Header { get; set; } = null!;

        /// <summary>
        /// Evaluates whether this matrix is symmetric.
        /// </summary>
        /// <remarks>
        /// Each pair must satisfy |A[i,j] - A[j,i]| &lt;= 1E-12 * max |A[k,l]|, so the test is invariant to the
        /// scale of the matrix. A pair involving NaN is not symmetric.
        /// </remarks>
        public bool IsSymmetric
        {
            get
            {
                if (!IsSquare) return false;
                double scale = 0d;
                for (int i = 0; i < NumberOfRows; i++)
                    for (int j = 0; j < NumberOfColumns; j++)
                        scale = Math.Max(scale, Math.Abs(_matrix[i, j]));
                double tolerance = 1e-12 * scale;
                for (int i = 0; i < NumberOfRows; i++)
                    for (int j = i + 1; j < NumberOfColumns; j++)
                        if (!(Math.Abs(this[i, j] - this[j, i]) <= tolerance))
                            return false;
                return true;
            }
        }

        /// <summary>
        /// Determines whether this matrix is square.
        /// </summary>
        public bool IsSquare => NumberOfRows == NumberOfColumns;

        /// <summary>
        /// Returns the underlying array as-is.
        /// </summary>
        public double[,] Array => _matrix;

        #endregion

        #region Methods

        /// <summary>
        /// Returns a copy of this matrix.  
        /// </summary>
        public Matrix Clone()
        {
            return new Matrix(_matrix);
        }

        /// <summary>
        /// Returns the matrix as XElement.
        /// </summary>
        public XElement ToXElement()
        {         
            var result = new XElement(nameof(Matrix));
            if (_matrix == null) return result;
            result.SetAttributeValue(nameof(NumberOfRows), NumberOfRows.ToString(CultureInfo.InvariantCulture));
            result.SetAttributeValue(nameof(NumberOfColumns), NumberOfColumns.ToString(CultureInfo.InvariantCulture));
            for (int i = 0; i < NumberOfRows; i++)
            {
                var rowElement = new XElement("row");
                var row = new double[NumberOfColumns];
                for (int j = 0; j < NumberOfColumns; j++)
                    row[j] = _matrix[i, j];
                
                rowElement.Value = string.Join("|", row.Select(v => v.ToString("G17", CultureInfo.InvariantCulture)));
                result.Add(rowElement);
            }
            return result;
        }

        /// <summary>
        /// Returns the matrix as an array.
        /// </summary>
        public double[,] ToArray()
        {
            return (double[,])_matrix.Clone();
        }

        /// <summary>
        /// Returns a copy of a specific row.
        /// </summary>
        /// <param name="index">The zero-based row index</param>
        public double[] Row(int index)
        {
            var vector = new double[NumberOfColumns];
            for (int j = 0; j < NumberOfColumns; j++)
                vector[j] = _matrix[index, j];
            return vector;
        }

        /// <summary>
        /// Returns a copy of a specific column.
        /// </summary>
        /// <param name="index">The zero-based column index</param>
        public double[] Column(int index)
        {
            var vector = new double[NumberOfRows];
            for (int i = 0; i < NumberOfRows; i++)
                vector[i] = _matrix[i, index];
            return vector;
        }

        #region Operations

        /// <summary>
        /// Returns a new matrix containing the upper triangle of this matrix.
        /// </summary>
        public Matrix UpperTriangle()
        {
            var result = new Matrix(NumberOfRows, NumberOfColumns);
            for (int i = 0; i < NumberOfRows; i++)
            {
                for (int j = i; j < NumberOfColumns; j++)
                    result[i, j] = this[i, j];
            }
            return result;
        }

        /// <summary>
        /// Returns a new matrix containing the lower triangle of this matrix.
        /// </summary>
        public Matrix LowerTriangle()
        {
            var result = new Matrix(NumberOfRows, NumberOfColumns);
            for (int i = 0; i < NumberOfRows; i++)
            {
                for (int j = 0; j < NumberOfColumns; j++)
                {
                    if (j > i) break;
                    result[i, j] = this[i, j];
                }
            }
            return result;
        }

        /// <summary>
        /// Returns the elements of the diagonal.
        /// </summary>
        /// <remarks>
        /// For non-square matrices, the method returns Min(Rows, Columns) elements where i = j (i is the row index, and j is the column index).
        /// </remarks>
        public double[] Diagonal()
        {
            int min = Math.Min(NumberOfRows, NumberOfColumns);
            var result = new double[min];
            for (int i = 0; i < min; i++)
                result[i] = this[i, i];
            return result;
        }

        /// <summary>
        /// Returns the matrix determinant.
        /// </summary>
        public double Determinant()
        {
            if (IsSquare == false)
                throw new ArgumentException("The matrix must be square.");

            var lU = new LUDecomposition(this);
            return lU.Determinant();
        }

        /// <summary>
        /// Returns the matrix inverse A^-1.
        /// </summary>
        public Matrix Inverse()
        {
            if (IsSquare == false)
                throw new ArgumentException("The matrix must be square.");
            var lU = new LUDecomposition(this);
            return lU.InverseA();
        }

        /// <summary>
        /// Returns the transpose the matrix.
        /// </summary>
        public Matrix Transpose()
        {
            return Transpose(this);
        }

        /// <summary>
        /// Returns the sum of the diagonal elements of a square matrix.
        /// </summary>
        public double Trace()
        {
            if (IsSquare == false)
                throw new ArgumentException("The matrix must be square.");

            double sum = 0;
            for (int i = 0; i < NumberOfRows; i++)
            {
                sum += this[i, i];
            }
            return sum;
        }

        #endregion

        #region Static Operations

        /// <summary>
        /// Returns the diagonal matrix.
        /// </summary>
        /// <param name="A">The square matrix.</param>
        public static Matrix Diagonal(Matrix A)
        {
            if (A.IsSquare == false)
                throw new ArgumentException("The matrix must be square.");
            var D = new Matrix(A.NumberOfColumns);
            for (int i = 0; i < A.NumberOfRows; i++)
            {
                for (int j = 0; j < A.NumberOfColumns; j++)
                {
                    D[i, j] = i == j ? A[i, j] : 0;
                }
            }
            return D;
        }

        /// <summary>
        /// Returns the diagonal matrix.
        /// </summary>
        /// <param name="A">The diagonal vector.</param>
        public static Matrix Diagonal(Vector A)
        {
            var D = new Matrix(A.Length);
            for (int i = 0; i < D.NumberOfRows; i++)
            {
                for (int j = 0; j < D.NumberOfColumns; j++)
                {
                    D[i, j] = i == j ? A[i] : 0;
                }
            }
            return D;
        }

        /// <summary>
        /// Returns the identity matrix of size n.
        /// </summary>
        /// <param name="size">The number of rows or columns (rows = columns).</param>
        /// <para>
        /// <see href = "https://en.wikipedia.org/wiki/Identity_matrix" />
        /// </para>
        public static Matrix Identity(int size)
        {
            var I = new Matrix(size);
            for (int j = 0; j < size; j++)
                I[j, j] = 1d;
            return I;
        }

        /// <summary>
        /// Returns the transposed matrix.
        /// </summary>
        /// <param name="A">The matrix to transpose.</param>
        /// <returns>The transpose of A.</returns>
        /// <remarks>
        /// The copy proceeds in square tiles so that both the source rows and the destination rows stay in cache.
        /// </remarks>
        public static Matrix Transpose(Matrix A)
        {
            int rows = A.NumberOfRows, cols = A.NumberOfColumns;
            var t = new Matrix(cols, rows);
            var a = A._matrix;
            var b = t._matrix;
            const int tile = 32;
            for (int i0 = 0; i0 < rows; i0 += tile)
            {
                int i1 = Math.Min(rows, i0 + tile);
                for (int j0 = 0; j0 < cols; j0 += tile)
                {
                    int j1 = Math.Min(cols, j0 + tile);
                    for (int i = i0; i < i1; i++)
                        for (int j = j0; j < j1; j++)
                            b[j, i] = a[i, j];
                }
            }
            return t;
        }

        #endregion

        #region Mathematics

        /// <summary>
        /// Computes the mean of each column.
        /// </summary>
        /// <returns>A vector where each element is the mean of the corresponding column in the matrix.</returns>
        public Vector ColumnMeans()
        {
            var means = new Vector(NumberOfColumns);
            for (int j = 0; j < NumberOfColumns; j++)
            {
                double sum = 0.0;
                for (int i = 0; i < NumberOfRows; i++)
                {
                    sum += _matrix[i, j];
                }
                means[j] = sum / NumberOfRows;
            }
            return means;
        }

        /// <summary>
        /// Applies a function to each element in the matrix.
        /// </summary>
        /// <param name="func">The function to apply to each element.</param>
        public void Apply(Func<double, double> func)
        {
            for (int i = 0; i < NumberOfRows; i++)
                for (int j = 0; j < NumberOfColumns; j++)
                    _matrix[i, j] = func(_matrix[i, j]);
        }

        /// <summary>
        /// Computes the square root of the matrix point wise. 
        /// </summary>
        public void Sqrt() => Apply(Math.Sqrt);

        /// <summary>
        /// Computes the square of the matrix point wise. 
        /// </summary>
        public void Sqr() => Apply(Tools.Sqr);

        /// <summary>
        /// Computes the log of the matrix point wise. 
        /// </summary>
        public void Log() => Apply(Math.Log);

        /// <summary>
        /// Computes the exponential of the matrix point wise. 
        /// </summary>
        public void Exp() => Apply(Math.Exp);

        /// <summary>
        /// Returns the sum of every point in the matrix.
        /// </summary>
        public double Sum()
        {
            double result = 0;
            for (int i = 0; i < NumberOfRows; i++)
                for (int j = 0; j < NumberOfColumns; j++)
                    result += this[i, j];
            return result;
        }

        /// <summary>
        /// Returns the outer product A = x yᵀ (size: x.Length × y.Length).
        /// </summary>
        public static Matrix Outer(Vector x, Vector y)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));

            int m = x.Length;
            int n = y.Length;
            var A = new Matrix(m, n);

            for (int i = 0; i < m; i++)
            {
                double xi = x[i];
                for (int j = 0; j < n; j++)
                    A[i, j] = xi * y[j];
            }
            return A;
        }

        /// <summary>
        /// Multiply by another matrix.
        /// </summary>
        /// <param name="matrix">The right-side matrix.</param>
        /// <returns>The product of this matrix and the right-side matrix.</returns>
        /// <exception cref="ArgumentException">Thrown when the number of columns of this matrix differs from the number of
        /// rows of the right-side matrix.</exception>
        /// <remarks>
        /// Each entry is accumulated from 0 over k in increasing order, C[i, j] = ((0 + A[i, 0]B[0, j]) + A[i, 1]B[1, j]) + ....
        /// For wider results each row of C is built as a sum of scaled rows of B, vectorized across j, which performs the
        /// same additions in the same order for every entry, so the result is identical to the plain triple loop.
        /// </remarks>
        public Matrix Multiply(Matrix matrix)
        {
            if (NumberOfColumns != matrix.NumberOfRows)
                throw new ArgumentException("The number of rows in the right-hand matrix must be equal to the number of columns in this matrix.");

            int rows = NumberOfRows, inner = NumberOfColumns, cols = matrix.NumberOfColumns;
            var result = new Matrix(rows, cols);
            var right = matrix._matrix;
            if (cols < VectorizedWidth)
            {
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        double sum = 0.0d;
                        for (int k = 0; k < inner; k++)
                            sum += _matrix[i, k] * right[k, j];
                        result._matrix[i, j] = sum;
                    }
                }
                return result;
            }
            var b = new double[inner * cols];
            Buffer.BlockCopy(right, 0, b, 0, b.Length * sizeof(double));
            var row = new double[cols];
            var c = result._matrix;
            for (int i = 0; i < rows; i++)
            {
                System.Array.Clear(row, 0, cols);
                for (int k = 0; k < inner; k++)
                    CrossProductKernel.Axpy(_matrix[i, k], b, k * cols, row, 0, cols);
                Buffer.BlockCopy(row, 0, c, i * cols * sizeof(double), cols * sizeof(double));
            }
            return result;
        }

        /// <summary>
        /// Multiply the transpose of this matrix by another matrix, Aᵀ·B, without forming the transpose.
        /// </summary>
        /// <param name="matrix">The right-side matrix, with the same number of rows as this matrix.</param>
        /// <returns>The product Aᵀ·B.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="matrix"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the numbers of rows differ.</exception>
        /// <remarks>
        /// <para>
        /// For a different right-hand matrix, each entry is accumulated from 0 over the rows k in increasing order, as in
        /// <see cref="Multiply(Matrix)"/>, so the result is identical to <c>Transpose().Multiply(matrix)</c>.
        /// </para>
        /// <para>
        /// When <paramref name="matrix"/> is this matrix, the Gram matrix AᵀA is accumulated by a SIMD rank-8 update over
        /// fixed blocks of rows, in parallel for large matrices, and mirrored from its upper triangle. The result is
        /// exactly symmetric and reproducible on a given machine. Its entries agree with the row-by-row sum to rounding
        /// (a few units in the last place), and because the SIMD register width sets the order of the partial sums, the
        /// last digits can differ between processor architectures.
        /// </para>
        /// </remarks>
        public Matrix TransposeMultiply(Matrix matrix)
        {
            if (matrix == null) throw new ArgumentNullException(nameof(matrix));
            if (NumberOfRows != matrix.NumberOfRows)
                throw new ArgumentException("The number of rows in the right-hand matrix must be equal to the number of rows in this matrix.");

            int inner = NumberOfRows, rows = NumberOfColumns, cols = matrix.NumberOfColumns;
            var result = new Matrix(rows, cols);
            var c = result._matrix;
            if (ReferenceEquals(matrix, this))
            {
                var gram = new double[rows * rows];
                CrossProductKernel.RowMajorGram(_matrix, new double[rows], new double[rows], gram);
                Buffer.BlockCopy(gram, 0, c, 0, gram.Length * sizeof(double));
                return result;
            }
            var a = new double[inner * rows];
            Buffer.BlockCopy(_matrix, 0, a, 0, a.Length * sizeof(double));
            var b = new double[inner * cols];
            Buffer.BlockCopy(matrix._matrix, 0, b, 0, b.Length * sizeof(double));
            var row = new double[cols];
            for (int i = 0; i < rows; i++)
            {
                System.Array.Clear(row, 0, cols);
                for (int k = 0; k < inner; k++)
                    CrossProductKernel.Axpy(a[k * rows + i], b, k * cols, row, 0, cols);
                Buffer.BlockCopy(row, 0, c, i * cols * sizeof(double), cols * sizeof(double));
            }
            return result;
        }

        /// <summary>
        /// Multiply by a vector.
        /// </summary>
        /// <param name="vector">The right-side vector.</param>
        public Vector Multiply(Vector vector)
        {
            if (NumberOfColumns != vector.Length)
                throw new ArgumentException("The number of rows in vector must be equal to the number of columns in the matrix.");
            var result = new double[NumberOfRows];
            for (int i = 0; i < NumberOfRows; i++)
            {
                double sum = 0.0d;
                for (int j = 0; j < NumberOfColumns; j++)
                    sum += _matrix[i, j] * vector[j];
                result[i] = sum;
            }
            return new Vector(result);
        }

        /// <summary>
        /// Multiply by a scalar.
        /// </summary>
        /// <param name="scalar">The scalar to multiply by.</param>
        public Matrix Multiply(double scalar)
        {
            var result = new Matrix(NumberOfRows, NumberOfColumns);
            for (int i = 0; i < NumberOfRows; i++)
                for (int j = 0; j < NumberOfColumns; j++)
                    result[i, j] = _matrix[i, j] * scalar;
            return result;
        }

        /// <summary>
        /// Divide by a scalar.
        /// </summary>
        /// <param name="scalar">The scalar value to divide by.</param>
        public Matrix Divide(double scalar)
        {
            var result = new Matrix(NumberOfRows, NumberOfColumns);
            for (int i = 0; i < NumberOfRows; i++)
                for (int j = 0; j < NumberOfColumns; j++)
                    result[i, j] = _matrix[i, j] / scalar;
            return result;
        }

        /// <summary>
        /// Add the matrix.
        /// </summary>
        /// <param name="matrix">The right-side matrix.</param>
        public Matrix Add(Matrix matrix)
        {
            if (NumberOfColumns != matrix.NumberOfColumns)
                throw new ArgumentException("The matrix must have the same number of columns.");
            if (NumberOfRows != matrix.NumberOfRows)
                throw new ArgumentException("The matrix must have the same number of rows.");

            var result = new Matrix(NumberOfRows, NumberOfColumns);
            for (int i = 0; i < NumberOfRows; i++)
                for (int j = 0; j < NumberOfColumns; j++)
                    result[i, j] = _matrix[i, j] + matrix[i, j];
            return result;
        }

        /// <summary>
        /// Subtract the matrix.
        /// </summary>
        /// <param name="matrix">The right-side matrix.</param>
        public Matrix Subtract(Matrix matrix)
        {
            if (NumberOfColumns != matrix.NumberOfColumns)
                throw new ArgumentException("The matrix must have the same number of columns.");
            if (NumberOfRows != matrix.NumberOfRows)
                throw new ArgumentException("The matrix must have the same number of rows.");

            var result = new Matrix(NumberOfRows, NumberOfColumns);
            for (int i = 0; i < NumberOfRows; i++)
                for (int j = 0; j < NumberOfColumns; j++)
                    result[i, j] = _matrix[i, j] - matrix[i, j];
            return result;
        }

        #endregion

        #region Operators

        /// <summary>
        /// Multiplies matrix A and B and returns the results as a matrix.
        /// </summary>
        /// <param name="A">Left-side matrix.</param>
        /// <param name="B">Right-side matrix.</param>
        public static Matrix operator *(Matrix A, Matrix B)
        {
            return A.Multiply(B);
        }

        /// <summary>
        /// Multiplies matrix A with a vector B and returns the results as an array vector.
        /// </summary>
        /// <param name="A">The matrix to multiply.</param>
        /// <param name="vector">The vector to multiply with.</param>
        public static double[] operator *(Matrix A, double[] vector)
        {
            return A.Multiply(new Vector(vector)).Array;
        }

        /// <summary>
        /// Multiplies matrix A with a vector B and returns the results as an array vector.
        /// </summary>
        /// <param name="A">The matrix to multiply.</param>
        /// <param name="B">The vector to multiply with.</param>
        public static Vector operator *(Matrix A, Vector B)
        {
            return A.Multiply(B);
        }

        /// <summary>
        /// Multiplies a matrix A with a scalar and returns the results as a matrix.
        /// </summary>
        /// <param name="matrix">The matrix to multiply.</param>
        /// <param name="scalar">The scalar to multiply by.</param>
        public static Matrix operator *(Matrix matrix, double scalar)
        {
            return matrix.Multiply(scalar);
        }

        /// <summary>
        /// Multiplies a matrix A with a scalar and returns the results as a matrix.
        /// </summary>
        /// <param name="scalar">The scalar to multiply by.</param>
        /// <param name="matrix">The matrix to multiply.</param>
        public static Matrix operator *(double scalar, Matrix matrix)
        {
            return matrix.Multiply(scalar);
        }

        /// <summary>
        /// Divides a matrix A with a scalar and returns the results as a matrix.
        /// </summary>
        /// <param name="matrix">The matrix to divide.</param>
        /// <param name="scalar">The scalar to divide by.</param>
        public static Matrix operator /(Matrix matrix, double scalar)
        {
            return matrix.Divide(scalar);
        }

        /// <summary>
        /// Add matrix A and B and returns the results as a matrix.
        /// </summary>
        /// <param name="A">Left-side matrix.</param>
        /// <param name="B">Right-side matrix.</param>
        public static Matrix operator +(Matrix A, Matrix B)
        {
            return A.Add(B);
        }

        /// <summary>
        /// Subtract matrix A and B and returns the results as a matrix.
        /// </summary>
        /// <param name="A">Left-side matrix.</param>
        /// <param name="B">Right-side matrix.</param>
        public static Matrix operator -(Matrix A, Matrix B)
        {
            return A.Subtract(B);
        }

        /// <summary>
        /// Returns the transposed matrix.
        /// </summary>
        /// <param name="matrix">The matrix to transpose.</param>
        public static Matrix operator ~(Matrix matrix)
        {
            return matrix.Transpose();
        }

        /// <summary>
        /// Returns the inverse of the matrix.
        /// </summary>
        /// <param name="matrix">The matrix to inverse.</param>
        public static Matrix operator !(Matrix matrix)
        {
            return matrix.Inverse();
        }

        #endregion

        #endregion

    }
}