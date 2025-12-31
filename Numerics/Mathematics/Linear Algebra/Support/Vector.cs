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

using Numerics.Data.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.LinearAlgebra
{

    /// <summary>
    /// A simple vector class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///    <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class Vector
    {
     
        /// <summary>
        /// Construct a new vector with specified length.
        /// </summary>
        /// <param name="length">The number of elements in the vector.</param>
        public Vector(int length)
        {
            _vector = new double[length];
        }

        /// <summary>
        /// Construct a new vector with specified length and fill it with a constant value.
        /// </summary>
        /// <param name="length"></param>
        /// <param name="fill">Fill the vector with a constant fill value.</param>
        public Vector(int length, double fill)
        {
            _vector = new double[length];
            for(int i = 0; i < length; i++)
            {
                _vector[i] = fill;
            }
        }

        /// <summary>
        /// Construct a new vector based an initial array.
        /// </summary>
        /// <param name="initialArray">Initializing array.</param>
        public Vector(double[] initialArray)
        {
            _vector = initialArray;
        }

        private double[] _vector;

        /// <summary>
        /// Returns the underlying array as-is.
        /// </summary>
        public double[] Array => _vector;

        /// <summary>
        /// The length of the vector.
        /// </summary>
        public int Length {
            get { return _vector.Length; }
        }

        /// <summary>
        /// Get the element at the specific index.
        /// </summary>
        /// <param name="index">The zero-based row index of the element to get or set.</param>
        public double this[int index]
        {
            get { return _vector[index]; }
            set { _vector[index] = value; }
        }

        /// <summary>
        /// The vector header text. 
        /// </summary>
        public string Header { get; set; } = null!;


        /// <summary>
        /// Returns the vector as an array.
        /// </summary>
        public double[] ToArray()
        {
            return (double[])_vector.Clone();
        }

        /// <summary>
        /// Returns the vector as a list.
        /// </summary>
        public List<double> ToList()
        {
            return _vector.ToList();
        }

        /// <summary>
        /// Returns a clone of the vector.
        /// </summary>
        public Vector Clone()
        {
            return new Vector(ToArray());
        }

        /// <summary>
        /// Clear the vector.
        /// </summary>
        public void Clear() => System.Array.Clear(_vector, 0, _vector.Length);


        /// <summary>
        /// Copy from another vector.
        /// </summary>
        /// <param name="other">The vector to copy from.</param>
        public void CopyFrom(Vector other)
        {
            if (other.Length != Length)
                throw new ArgumentException("Vectors must be the same length.");
            System.Array.Copy(other._vector, _vector, Length);
        }

        /// <summary>
        /// Returns the sum of the vector.
        /// </summary>
        public double Sum()
        {
            return _vector.Sum();
        }

        /// <summary>
        /// Returns the vector Norm.
        /// </summary>
        public double Norm()
        {
            double d = 0;
            for (int i = 0; i < Length; i++)
                d += _vector[i] * _vector[i];
            return Math.Sqrt(d);
        }

        /// <summary>
        /// Returns the vector Norm squared.
        /// </summary>
        public double NormSquared()
        {
            double d = 0;
            for (int i = 0; i < Length; i++)
                d += _vector[i] * _vector[i];
            return d;
        }

        /// <summary>
        /// Returns the Euclidean distance between two vectors ||x - y||.
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static double Distance(Vector A, Vector B)
        {
            if (A.Length != B.Length) throw new ArgumentException(nameof(A.Length), "The vectors must be the same length.");
            double d = 0;
            for (int i = 0; i < A.Length; i++)
            {
                double dx = A[i] - B[i];
                d += dx * dx;
            }
            return Math.Sqrt(d);
        }

        /// <summary>
        /// Returns the dot product of two vectors.
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static double DotProduct(Vector A, Vector B)
        {
            if (A.Length != B.Length) throw new ArgumentException(nameof(A.Length), "The vectors must be the same length.");
            double sum = 0;
            for (int i = 0; i < A.Length; i++)
                sum += A[i] * B[i];
            return sum;
        }

        /// <summary>
        /// Project vector A onto B.
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector Project(Vector A, Vector B)
        {
            var ab = DotProduct(A, B);
            var bb = B.NormSquared();
            return B * (ab / bb);
        }

        /// <summary>
        /// Add the vector.
        /// </summary>
        /// <param name="vector">The right-side vector.</param>
        public Vector Add(Vector vector)
        {
            if (Length != vector.Length)
                throw new ArgumentException("Vectors must be the same length.");
            var result = new Vector(Length);
            for (int i = 0; i < Length; i++)
                result[i] = _vector[i] + vector[i];
            return result;
        }

        /// <summary>
        /// Add the array.
        /// </summary>
        /// <param name="vector">The right-side array.</param>
        public Vector Add(double[] vector)
        {
            if (Length != vector.Length)
                throw new ArgumentException("Vectors must be the same length.");
            var result = new Vector(Length);
            for (int i = 0; i < Length; i++)
                result[i] = _vector[i] + vector[i];
            return result;
        }

        /// <summary>
        /// Subtract the vector.
        /// </summary>
        /// <param name="vector">The right-side vector.</param>
        public Vector Subtract(Vector vector)
        {
            if (Length != vector.Length)
                throw new ArgumentException("Vectors must be the same length.");
            var result = new Vector(Length);
            for (int i = 0; i < Length; i++)
                result[i] = _vector[i] - vector[i];
            return result;
        }

        /// <summary>
        /// Subtract the array.
        /// </summary>
        /// <param name="vector">The right-side array.</param>
        public Vector Subtract(double[] vector)
        {
            if (Length != vector.Length)
                throw new ArgumentException("Vectors must be the same length.");
            var result = new Vector(Length);
            for (int i = 0; i < Length; i++)
                result[i] = _vector[i] - vector[i];
            return result;
        }

        /// <summary>
        /// Multiply by a matrix.
        /// </summary>
        /// <param name="matrix">The right-side matrix.</param>
        public Vector Multiply(Matrix matrix)
        {
            if (matrix.NumberOfColumns != Length)
                throw new ArgumentException("The number of rows in vector must be equal to the number of columns in the matrix.");
            var result = new double[matrix.NumberOfRows];
            for (int i = 0; i < matrix.NumberOfRows; i++)
            {
                double sum = 0.0d;
                for (int j = 0; j < matrix.NumberOfColumns; j++)
                    sum += matrix[i, j] * _vector[j];
                result[i] = sum;
            }
            return new Vector(result);
        }

        /// <summary>
        /// Multiply by a vector.
        /// </summary>
        /// <param name="vector">The right-side vector.</param>
        public Vector Multiply(Vector vector)
        {
            if (Length != vector.Length) throw new ArgumentException(nameof(Length), "The vectors must be the same length.");
            var result = new Vector(Length);
            for (int i = 0; i < Length; i++)
                result[i] = _vector[i] * vector[i];
            return result;
        }

        /// <summary>
        /// Multiply by an array.
        /// </summary>
        /// <param name="vector">The right-side array.</param>
        public Vector Multiply(double[] vector)
        {
            if (Length != vector.Length) throw new ArgumentException(nameof(Length), "The vectors must be the same length.");
            var result = new Vector(Length);
            for (int i = 0; i < Length; i++)
                result[i] = _vector[i] * vector[i];
            return result;
        }

        /// <summary>
        /// Multiply by a scalar.
        /// </summary>
        /// <param name="scalar">The scalar to multiply by.</param>
        public Vector Multiply(double scalar)
        {
            var result = new Vector(Length);
            for (int i = 0; i < Length; i++)
                result[i] = _vector[i] * scalar;
            return result;
        }

        /// <summary>
        /// Divide by a scalar.
        /// </summary>
        /// <param name="scalar">The scalar to divide by.</param>
        public Vector Divide(double scalar)
        {
            var result = new Vector(Length);
            for (int i = 0; i < Length; i++)
                result[i] = _vector[i] / scalar;
            return result;
        }

        /// <summary>
        /// Multiplies vectors A and matrix B. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side matrix.</param>
        public static Vector operator *(Vector A, Matrix B)
        {
            return A.Multiply(B);
        }

        /// <summary>
        /// Multiplies vectors A and B by multiplying corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator *(Vector A, Vector B)
        {
            return A.Multiply(B);
        }

        /// <summary>
        /// Multiplies vectors A and B by multiplying corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator *(Vector A, double[] B)
        {
            return A.Multiply(B);
        }

        /// <summary>
        /// Multiplies a vector A with a scalar.
        /// </summary>
        /// <param name="A">The vector to multiply.</param>
        /// <param name="scalar">The scalar to multiply by.</param>
        public static Vector operator *(Vector A, double scalar)
        {
            return A.Multiply(scalar);
        }

        /// <summary>
        /// Divides a vector A by a scalar.
        /// </summary>
        /// <param name="A">The vector to multiply.</param>
        /// <param name="scalar">The scalar to divide by.</param>
        public static Vector operator /(Vector A, double scalar)
        {
            return A.Divide(scalar);
        }

        /// <summary>
        /// Raises each element of the vector A by the exponent b.
        /// </summary>
        /// <param name="A">The vector to multiply.</param>
        /// <param name="b">The exponent b..</param>
        public static Vector operator ^(Vector A, double b)
        {          
            var v = new Vector(A.Length);
            for (int i = 0; i < A.Length; i++)
                v[i] = Math.Pow(A[i], b);
            return v;
        }

        /// <summary>
        /// Adds vectors A and B by summing corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator +(Vector A, Vector B)
        {
            return A.Add(B);
        }

        /// <summary>
        /// Adds vectors A and B by summing corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator +(Vector A, double[] B)
        {
            return A.Add(B);
        }

        /// <summary>
        /// Subtracts vectors A and B by subtracting corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator -(Vector A, Vector B)
        {
            return A.Subtract(B);
        }

        /// <summary>
        /// Subtracts vectors A and B by subtracting corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator -(Vector A, double[] B)
        {
            return A.Subtract(B);
        }

    }
}