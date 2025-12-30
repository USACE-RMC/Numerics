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

using Numerics.Sampling;
using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;

namespace Numerics
{
    /// <summary>
    /// A class for extension methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    ///     </list>
    /// </para>
    /// </remarks>
    public static class ExtensionMethods
    {

        #region Enum

        /// <summary>
        /// Gets an attribute on an enum field value
        /// </summary>
        /// <typeparam name="T">The type of the attribute you want to retrieve</typeparam>
        /// <param name="enumValue">The enum value</param>
        /// <returns>The attribute of type T that exists on the enum value</returns>
        public static T GetAttributeOfType<T>(this Enum enumValue) where T : Attribute
        {
            var type = enumValue.GetType();
            var memInfo = type.GetMember(enumValue.ToString());
            var attributes = memInfo[0].GetCustomAttributes(typeof(T), false);
            return (attributes.Length > 0) ? (T)attributes[0] : null!;
        }

        #endregion

        #region Double

        /// <summary>
        /// Determines if the values are almost equal to each other.
        /// </summary>
        /// <param name="a">The left-side value.</param>
        /// <param name="b">The right-side value.</param>
        /// <param name="epsilon">The absolute tolerance level. Default = 1E-15.</param>
        public static bool AlmostEquals(this double a, double b, double epsilon = 1E-15)
        {
            return Math.Abs(a - b) <= epsilon;
        }

        #endregion

        #region Random

        /// <summary>
        /// Returns an array of random integers.
        /// </summary>
        /// <param name="random">A random number generator.</param>
        /// <param name="length">The number of samples to return.</param>
        public static int[] NextIntegers(this Random random, int length)
        {
            var values = new int[length];
            for (int i = 0; i < length; i++)
                values[i] = random.Next();
            return values;
        }

        /// <summary>
        /// Returns an array of random integers between a min and max value. 
        /// </summary>
        /// <param name="random">A random number generator.</param>
        /// <param name="minValue">The minimum value to sample between.</param>
        /// <param name="maxValue">The maximum value to sample between.</param>
        /// <param name="length">The number of samples to return.</param>
        public static int[] NextIntegers(this Random random, int minValue, int maxValue, int length)
        {
            var values = new int[length];
            for (int i = 0; i < length; i++)
                values[i] = random.Next(minValue, maxValue);
            return values;
        }

        /// <summary>
        /// Returns an array of random integers between a min and max value. 
        /// </summary>
        /// <param name="random">A random number generator.</param>
        /// <param name="minValue">The minimum value to sample between.</param>
        /// <param name="maxValue">The maximum value to sample between.</param>
        /// <param name="length">The number of samples to return.</param>
        /// <param name="replace">Determines whether or not to sample with replacement. Default = true.</param>
        public static int[] NextIntegers(this Random random, int minValue, int maxValue, int length, bool replace = true)
        {
            if (replace == true)
            {
                var values = new int[length];
                for (int i = 0; i < length; i++)
                    values[i] = random.Next(minValue, maxValue);
                return values;
            }
            else
            {       
                if (length > maxValue - minValue)
                    throw new ArgumentException("When sampling without replacement, the length must be less than or equal to the range of values.");

                var bins = new List<int>();
                for (int i = minValue; i < maxValue; i++)
                    bins.Add(i);

                // Sample random bin without replacement
                var values = new int[length];
                for (int i = 0; i < length; i++)
                {
                    int r = random.Next(0, bins.Count);
                    values[i] = bins[r];
                    bins.RemoveAt(r);
                }
                return values;
            }

        }

        /// <summary>
        /// Returns an array of random doubles.
        /// </summary>
        /// <param name="random">A random number generator.</param>
        /// <param name="length">The number of samples to return.</param>
        public static double[] NextDoubles(this Random random, int length)
        {
            var values = new double[length];
            for (int i = 0; i < length; i++)
                values[i] = random.NextDouble();
            return values;
        }

        /// <summary>
        /// Returns a 2-D array of random doubles.
        /// </summary>
        /// <param name="random">A random number generator.</param>
        /// <param name="length">The number of samples to return.</param>
        /// <param name="dimension">The spatial dimension</param>
        public static double[,] NextDoubles(this Random random, int length, int dimension)
        {
            var values = new double[length, dimension];
            var prngs = new Random[dimension];
            for (int i = 0; i < dimension; i++)
            {
                prngs[i] = random.GetType() == typeof(MersenneTwister) ? new MersenneTwister(random.Next()) : new Random(random.Next());
                for (int j = 0; j < length; j++)
                {
                    values[j, i] = prngs[i].NextDouble();
                }
            }
            return values;
        }

        #endregion

        #region 1-D Array

        /// <summary>
        /// Applies a function to each element in the array.
        /// </summary>
        /// <param name="array">The array.</param>
        /// <param name="func">The function to apply to each element.</param>
        public static void Apply(this double[] array, Func<double, double> func)
        {
            for (int i = 0; i < array.Length; i++)
                array[i] = func(array[i]);
        }

        /// <summary>
        /// Applies a function to each element in the array
        /// </summary>
        /// <param name="array">The array.</param>
        /// <param name="func">The function to apply to each element.</param>
        public static double[] Map(this double[] array, Func<double, double> func)
        {
            var result = new double[array.Length];
            for (int i = 0; i < array.Length; i++)
                result[i] = func(array[i]);
            return result;
        }

        /// <summary>
        /// Adds corresponding elements of arrays.
        /// </summary>
        /// <param name="array">The array.</param>
        /// <param name="values">The values to add.</param>
        /// <returns>The array after addition.</returns>
        public static double[] Add(this double[] array, double[] values)
        {
            if (array.Length != values.Length) throw new ArgumentException(nameof(array), "The arrays must be the same length.");
            var result = new double[array.Length];
            for (int i = 0; i < array.Length; i++)
                result[i] = array[i] + values[i];
            return result;
        }

        /// <summary>
        /// Subtracts corresponding elements of arrays.
        /// </summary>
        /// <param name="array">The array</param>
        /// <param name="values">The values to subtract.</param>
        /// <returns>The array after subtraction.</returns>
        public static double[] Subtract(this double[] array, double[] values)
        {
            if (array.Length != values.Length) throw new ArgumentException(nameof(array), "The arrays must be the same length.");
            var result = new double[array.Length];
            for (int i = 0; i < array.Length; i++)
                result[i] = array[i] - values[i];
            return result;
        }

        /// <summary>
        /// Multiplies an array with a scalar.
        /// </summary>
        /// <param name="array">The array</param>
        /// <param name="scalar">The scalar to multiply by.</param>
        /// <returns>The array after multiplication.</returns>
        public static double[] Multiply(this double[] array, double scalar)
        {
            var result = new double[array.Length];
            for (int i = 0; i < array.Length; i++)
                result[i] = array[i] * scalar;
            return result;
        }

        /// <summary>
        /// Divides an array with a scalar.
        /// </summary>
        /// <param name="array">The array</param>
        /// <param name="scalar">The scalar to divide by.</param>
        /// <returns>The array after division.</returns>
        public static double[] Divide(this double[] array, double scalar)
        {
            var result = new double[array.Length];
            for (int i = 0; i < array.Length; i++)
                result[i] = array[i] / scalar;
            return result;
        }

        /// <summary>
        /// Returns the dot product of two arrays of values.
        /// </summary>
        /// <param name="array">The array</param>
        /// <param name="values">The values to multiply.</param>
        /// <returns>The array after performing the dot product.</returns>
        public static double DotProduct(this double[] array, double[] values)
        {
            if (array.Length != values.Length) throw new ArgumentException(nameof(array), "The arrays must be the same length.");
            double result = 0.0;
            for (int i = 0; i < array.Length; i++)
                result += array[i] * values[i];
            return result;
        }

        /// <summary>
        /// Returns a subset of the array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="startIndex">The start index.</param>
        public static T[] Subset<T>(this T[] array, int startIndex)
        {
            int t = 0;
            var result = new T[array.Length - startIndex];
            for (int i = startIndex; i < array.Length; i++)
            {
                result[t] = array[i];
                t++;
            }
            return result;
        }

        /// <summary>
        /// Returns a subset of the array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="startIndex">The start index.</param>
        /// <param name="endIndex">The end index.</param>
        public static T[] Subset<T>(this T[] array, int startIndex, int endIndex)
        {
            int t = 0;
            var result = new T[endIndex - startIndex + 1];
            for (int i = startIndex; i <= endIndex; i++)
            {
                result[t] = array[i];
                t++;
            }
            return result;
        }

        /// <summary>
        /// Returns a random subset of the array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="length">The number of samples to return.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="replace">Optional. Determines whether or not to sample with replacement. Default = false.</param>
        public static T[] RandomSubset<T>(this T[] array, int length, int seed = -1, bool replace = false)
        {
            if (length > array.Length)
                throw new ArgumentException("The subset length must be less than or equal to the array length.");

            Random random = seed > 0 ? new Random(seed) : new Random();
            var indexes = random.NextIntegers(0, array.Length, length, replace);
            var result = new T[indexes.Length];
            for (int i = 0; i < indexes.Length; i++)
                result[i] = array[indexes[i]];
            return result;
        }

        /// <summary>
        /// Fills an array with the specified value.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="value">The fill value.</param>
        public static void Fill<T>(this T[] array, T value)
        {
            for (int i = 0; i < array.Length; i++)
                array[i] = value;
        }

        #endregion

        #region 2-D Array

        /// <summary>
        /// Applies a function to each element in the array.
        /// </summary>
        /// <param name="array">The array.</param>
        /// <param name="func">The function to apply to each element.</param>
        public static void Apply(this double[,] array, Func<double, double> func)
        {
            for (int i = 0; i < array.GetLength(0); i++)
                for (int j = 0; j < array.GetLength(1); j++)
                    array[i, j] = func(array[i, j]);
        }

        /// <summary>
        /// Gets a specific column from a 2-D array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="index">Zero-based index of the column.</param>
        public static T[] GetColumn<T>(this T[,] array, int index)
        {
            var col = new T[array.GetLength(0)];
            for (int i = 0; i < array.GetLength(0); i++)
                col[i] = array[i, index];
            return col;
        }

        /// <summary>
        /// Gets a specific row from a 2-D array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="index">Zero-based index of the row.</param>
        public static T[] GetRow<T>(this T[,] array, int index)
        {
            var row = new T[array.GetLength(1)];
            for (int i = 0; i < array.GetLength(1); i++)
                row[i] = array[index, i];
            return row;
        }

        /// <summary>
        /// Sets a specific row in a 2-D array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="index">Zero-based index of the row.</param>
        /// <param name="values">The new values.</param>
        public static void SetRow<T>(this T[,] array, int index, T[] values)
        {
            for (int i = 0; i < array.GetLength(1); i++)
                array[index, i] = values[i];
        }

        /// <summary>
        /// Sets a specific column in a 2-D array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="index">Zero-based index of the column.</param>
        /// <param name="values">The new values.</param>
        public static void SetColumn<T>(this T[,] array, int index, T[] values)
        {
            for (int i = 0; i < array.GetLength(0); i++)
                array[i, index] = values[i];
        }


        /// <summary>
        /// Returns a subset of the array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="startIndex">The start index.</param>
        public static T[,] Subset<T>(this T[,] array, int startIndex)
        {
            int t = 0;
            var result = new T[array.GetLength(0) - startIndex, array.GetLength(1)];
            for (int i = startIndex; i < array.GetLength(0); i++)
            {
                for (int j = 0; j < array.GetLength(1); j++)
                {
                    result[t, j] = array[i, j];
                }
                t++;
            }
            return result;
        }

        /// <summary>
        /// Returns a subset of the array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="startIndex">The start index.</param>
        /// <param name="endIndex">The end index.</param>
        public static T[,] Subset<T>(this T[,] array, int startIndex, int endIndex)
        {
            int t = 0;
            var result = new T[endIndex - startIndex + 1, array.GetLength(1)];
            for (int i = startIndex; i <= endIndex; i++)
            {
                for (int j = 0; j < array.GetLength(1); j++)
                {
                    result[t, j] = array[i, j];
                }
                t++;
            }
            return result;
        }

        /// <summary>
        /// Returns a random subset of the array.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="length">The number of samples to return.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="replace">Optional. Determines whether or not to sample with replacement. Default = false.</param>
        public static T[,] RandomSubset<T>(this T[,] array, int length, int seed = -1, bool replace = false)
        {
            if (length > array.GetLength(0))
                throw new ArgumentException("The subset length must be less than or equal to the array length.");

            Random random = seed > 0 ? new Random(seed) : new Random();
            var indexes = random.NextIntegers(0, array.GetLength(0), length, replace);
            var result = new T[indexes.Length, array.GetLength(1)];
            for (int i = 0; i < indexes.Length; i++)
            {
                for (int j = 0; j < array.GetLength(1); j++)
                {
                    result[i, j] = array[indexes[i], j];
                }
            }
            return result;
        }

        /// <summary>
        /// Fills a 2-D array with the specified value.
        /// </summary>
        /// <typeparam name="T">The array value type.</typeparam>
        /// <param name="array">The array.</param>
        /// <param name="value">The fill value.</param>
        public static void Fill<T>(this T[,] array, T value)
        {
            for (int i = 0; i < array.GetLength(0); i++)
                for (int j = 0; j < array.GetLength(1); j++)
                    array[i, j] = value;
        }

        #endregion

        #region Vector

        /// <summary>
        /// Returns a subset of the vector.
        /// </summary>
        /// <param name="vector">The vector.</param>
        /// <param name="startIndex">The start index.</param>
        public static Vector Subset(this Vector vector, int startIndex)
        {
            int t = 0;
            var result = new Vector(vector.Length - startIndex);
            for (int i = startIndex; i < vector.Length; i++)
            {
                result[t] = vector[i];
                t++;
            }
            return result;
        }

        /// <summary>
        /// Returns a subset of the vector.
        /// </summary>
        /// <param name="vector">The vector.</param>
        /// <param name="startIndex">The start index.</param>
        /// <param name="endIndex">The end index.</param>
        public static Vector Subset(this Vector vector, int startIndex, int endIndex)
        {
            int t = 0;
            var result = new Vector(endIndex - startIndex + 1);
            for (int i = startIndex; i <= endIndex; i++)
            {
                result[t] = vector[i];
                t++;
            }
            return result;
        }

        /// <summary>
        /// Returns a random subset of the vector.
        /// </summary>
        /// <param name="vector">The vector.</param>
        /// <param name="length">The number of samples to return.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="replace">Optional. Determines whether or not to sample with replacement. Default = false.</param>
        public static Vector RandomSubset(this Vector vector, int length, int seed = -1, bool replace = false)
        {
            if (length > vector.Length)
                throw new ArgumentException("The subset length must be less than or equal to the vector length.");

            Random random = seed > 0 ? new Random(seed) : new Random();
            var indexes = random.NextIntegers(0, vector.Length, length, replace);
            var result = new Vector(indexes.Length);
            for (int i = 0; i < indexes.Length; i++)
                result[i] = vector[indexes[i]];
            return result;
        }

        /// <summary>
        /// Fills a vector with the specified value.
        /// </summary>
        /// <param name="vector">The vector.</param>
        /// <param name="value">The fill value.</param>
        public static void Fill(this Vector vector, double value)
        {
            for (int i = 0; i < vector.Length; i++)
                vector[i] = value;
        }

        #endregion

        #region Matrix

        /// <summary>
        /// Gets a specific column from a matrix.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="index">Zero-based index of the column.</param>
        public static double[] GetColumn(this Matrix matrix, int index)
        {
            var col = new double[matrix.NumberOfRows];
            for (int i = 0; i < matrix.NumberOfRows; i++)
                col[i] = matrix[i, index];
            return col;
        }

        /// <summary>
        /// Gets a specific row from a matrix.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="index">Zero-based index of the row.</param>
        public static double[] GetRow(this Matrix matrix, int index)
        {
            var row = new double[matrix.NumberOfColumns];
            for (int i = 0; i < matrix.NumberOfColumns; i++)
                row[i] = matrix[index, i];
            return row;
        }

        /// <summary>
        /// Sets a specific row in a matrix.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="index">Zero-based index of the row.</param>
        /// <param name="values">The new values.</param>
        public static void SetRow(this Matrix matrix, int index, double[] values)
        {
            for (int i = 0; i < matrix.NumberOfColumns; i++)
                matrix[index, i] = values[i];
        }

        /// <summary>
        /// Sets a specific column in a matrix.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="index">Zero-based index of the column.</param>
        /// <param name="values">The new values.</param>
        public static void SetColumn(this Matrix matrix, int index, double[] values)
        {
            for (int i = 0; i < matrix.NumberOfRows; i++)
                matrix[i, index] = values[i];
        }

        /// <summary>
        /// Returns a subset of the matrix.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="startIndex">The start index.</param>
        public static Matrix Subset(this Matrix matrix, int startIndex)
        {
            int t = 0;
            var result = new Matrix(matrix.NumberOfRows - startIndex, matrix.NumberOfColumns);
            for (int i = startIndex; i < matrix.NumberOfRows; i++)
            {
                for (int j = 0; j < matrix.NumberOfColumns; j++)
                {
                    result[t, j] = matrix[i, j];
                }
                t++;
            }
            return result;
        }

        /// <summary>
        /// Returns a subset of the matrix.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="startIndex">The start index.</param>
        /// <param name="endIndex">The end index.</param>
        public static Matrix Subset(this Matrix matrix, int startIndex, int endIndex)
        {
            int t = 0;
            var result = new Matrix(endIndex - startIndex + 1, matrix.NumberOfColumns);
            for (int i = startIndex; i <= endIndex; i++)
            {
                for (int j = 0; j < matrix.NumberOfColumns; j++)
                {
                    result[t, j] = matrix[i, j];
                }
                t++;
            }
            return result;
        }

        /// <summary>
        /// Returns a random subset of the matrix.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="length">The number of samples to return.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <param name="replace">Optional. Determines whether or not to sample with replacement. Default = false.</param>
        public static Matrix RandomSubset(this Matrix matrix, int length, int seed = -1, bool replace = false)
        {
            if (length > matrix.NumberOfRows)
                throw new ArgumentException("The subset length must be less than or equal to the number of rows in the matrix.");

            Random random = seed > 0 ? new Random(seed) : new Random();
            var indexes = random.NextIntegers(0, matrix.NumberOfRows, length, replace);
            var result = new Matrix(indexes.Length, matrix.NumberOfColumns);
            for (int i = 0; i < indexes.Length; i++)
            {
                for (int j = 0; j < matrix.NumberOfColumns; j++)
                {
                    result[i, j] = matrix[indexes[i], j];
                }
            }
            return result;
        }

        /// <summary>
        /// Fills a matrix with the specified value.
        /// </summary>
        /// <param name="matrix">The matrix.</param>
        /// <param name="value">The fill value.</param>
        public static void Fill(this Matrix matrix, double value)
        {
            for (int i = 0; i < matrix.NumberOfRows; i++)
                for (int j = 0; j < matrix.NumberOfColumns; j++)
                    matrix[i, j] = value;
        }

        #endregion

    }
}
