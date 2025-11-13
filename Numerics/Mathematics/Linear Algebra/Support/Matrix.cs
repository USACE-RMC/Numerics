/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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

            // Deep copy - compatible with .NET Framework 4.8.1 and later
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    _matrix[i, j] = initialArray[i, j];
                }
            }
        }

        /// <summary>
        /// Constructs a new matrix based on a single column array.
        /// </summary>
        /// <param name="initialArray">Initializing array.</param>
        public Matrix(double[] initialArray)
        {
            _matrix = new double[initialArray.Length, 1];
            for (int i = 0; i < _matrix.Length; i++)
                _matrix[i, 0] = initialArray[i];
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
            if (xElement.Attribute(nameof(NumberOfRows)) != null) int.TryParse(xElement.Attribute(nameof(NumberOfRows)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out nrow);
            if (xElement.Attribute(nameof(NumberOfColumns)) != null) int.TryParse(xElement.Attribute(nameof(NumberOfColumns)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out ncol);

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

        private double[,] _matrix;

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
        public string[] Header { get; set; }

        /// <summary>
        /// Evaluates whether this matrix is symmetric.
        /// </summary>
        public bool IsSymmetric
        {
            get
            {
                if (!IsSquare) return false;
                const double epsilon = 1e-12;
                for (int i = 0; i < NumberOfRows; i++)
                    for (int j = i + 1; j < NumberOfColumns; j++)
                        if (Math.Abs(this[i, j] - this[j, i]) > epsilon)
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
            return new Matrix(ToArray());
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
            var t = new Matrix(NumberOfColumns, NumberOfRows);
            for (int i = 0; i < NumberOfRows; i++)
            {
                for (int j = 0; j < NumberOfColumns; j++)
                    t[j, i] = this[i, j];
            }
            return t;
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
        public static Matrix Transpose(Matrix A)
        {
            var t = new Matrix(A.NumberOfColumns, A.NumberOfRows);
            for (int i = 0; i < A.NumberOfRows; i++)
            {
                for (int j = 0; j < A.NumberOfColumns; j++)
                    t[j, i] = A[i, j];
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
        public Matrix Multiply(Matrix matrix)
        {
            if (NumberOfColumns != matrix.NumberOfRows)
                throw new ArgumentException("The number of rows in the right-hand matrix must be equal to the number of columns in this matrix.");

            var result = new Matrix(NumberOfRows, matrix.NumberOfColumns);
            for (int i = 0; i < NumberOfRows; i++)
            {
                for (int j = 0; j < matrix.NumberOfColumns; j++)
                {
                    double sum = 0.0d;
                    for (int k = 0; k < NumberOfColumns; k++)
                        sum += this[i, k] * matrix[k, j];
                    result[i, j] = sum;
                }
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