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
            _array = new double[length];
        }

        /// <summary>
        /// Construct a new vector with specified length and fill it with a constant value.
        /// </summary>
        /// <param name="length"></param>
        /// <param name="fill">Fill the vector with a constant fill value.</param>
        public Vector(int length, double fill)
        {
            _array = new double[length];
            for(int i = 0; i < length; i++)
            {
                _array[i] = fill;
            }
        }

        /// <summary>
        /// Construct a new vector based an initial array.
        /// </summary>
        /// <param name="initialArray">Initializing array.</param>
        public Vector(double[] initialArray)
        {
            _array = initialArray;
        }

        private double[] _array;

        /// <summary>
        /// Returns the underlying array as-is.
        /// </summary>
        public double[] Array => _array;

        /// <summary>
        /// The length of the vector.
        /// </summary>
        public int Length {
            get { return _array.Length; }
        }

        /// <summary>
        /// Get the element at the specific index.
        /// </summary>
        /// <param name="index">The zero-based row index of the element to get or set.</param>
        public double this[int index]
        {
            get { return _array[index]; }
            set { _array[index] = value; }
        }

         /// <summary>
         /// The vector header text. 
         /// </summary>
        public string Header { get; set; }


        /// <summary>
        /// Returns the vector as an array.
        /// </summary>
        public double[] ToArray()
        {
            return (double[])_array.Clone();
        }

        /// <summary>
        /// Returns the vector as a list.
        /// </summary>
        public List<double> ToList()
        {
            return _array.ToList();
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
        public void Clear() => System.Array.Clear(_array, 0, _array.Length);


        /// <summary>
        /// Copy from another vector.
        /// </summary>
        /// <param name="other">The vector to copy from.</param>
        public void CopyFrom(Vector other)
        {
            if (other.Length != Length)
                throw new ArgumentException("Vectors must be the same length.");
            System.Array.Copy(other._array, _array, Length);
        }

        /// <summary>
        /// Returns the vector Norm.
        /// </summary>
        public double Norm()
        {
            double d = 0;
            for (int i = 0; i < Length; i++)
                d += _array[i] * _array[i];
            return Math.Sqrt(d);
        }

        /// <summary>
        /// Returns the vector Norm squared.
        /// </summary>
        public double NormSquared()
        {
            double d = 0;
            for (int i = 0; i < Length; i++)
                d += _array[i] * _array[i];
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
        /// Multiplies a vector A with a scalar.
        /// </summary>
        /// <param name="A">The vector to multiply.</param>
        /// <param name="scalar">The scalar to multiply by.</param>
        public static Vector operator *(Vector A, double scalar)
        {
            var v = new Vector(A.Length);
            for (int i = 0; i < A.Length; i++)
                v[i] = A[i] * scalar;
            return v;
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
        /// Multiplies vectors A and B by multiplying corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator *(Vector A, Vector B)
        {
            if (A.Length != B.Length) throw new ArgumentException(nameof(A.Length), "The vectors must be the same length.");
            var v = new Vector(A.Length);
            for (int i = 0; i < A.Length; i++)
                v[i] = A[i] * B[i];
            return v;
        }

        /// <summary>
        /// Adds vectors A and B by summing corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator +(Vector A, Vector B)
        {
            if (A.Length != B.Length) throw new ArgumentException(nameof(A.Length), "The vectors must be the same length.");
            var v = new Vector(A.Length);
            for (int i = 0; i < A.Length; i++)
                v[i] = A[i] + B[i];
            return v;
        }

        /// <summary>
        /// Subtracts vectors A and B by subtracting corresponding elements. Vectors A and B must the same size. 
        /// </summary>
        /// <param name="A">Left-side vector.</param>
        /// <param name="B">Right-side vector.</param>
        public static Vector operator -(Vector A, Vector B)
        {
            if (A.Length != B.Length) throw new ArgumentException(nameof(A.Length), "The vectors must be the same length.");
            var v = new Vector(A.Length);
            for (int i = 0; i < A.Length; i++)
                v[i] = A[i] - B[i];
            return v;
        }

        /// <summary>
        /// Performs in-place addition: result = A + B.
        /// </summary>
        public static void Add(Vector A, Vector B, Vector result)
        {
            if (A.Length != B.Length || result.Length != A.Length)
                throw new ArgumentException("Vectors must be the same length.");
            for (int i = 0; i < A.Length; i++)
                result[i] = A[i] + B[i];
        }

        /// <summary>
        /// Performs in-place subtraction: result = A - B.
        /// </summary>
        public static void Subtract(Vector A, Vector B, Vector result)
        {
            if (A.Length != B.Length || result.Length != A.Length)
                throw new ArgumentException("Vectors must be the same length.");
            for (int i = 0; i < A.Length; i++)
                result[i] = A[i] - B[i];
        }

        /// <summary>
        /// Performs in-place scalar multiplication: result = A * scalar.
        /// </summary>
        public static void Multiply(Vector A, double scalar, Vector result)
        {
            if (result.Length != A.Length)
                throw new ArgumentException("Vectors must be the same length.");
            for (int i = 0; i < A.Length; i++)
                result[i] = A[i] * scalar;
        }

    }
}