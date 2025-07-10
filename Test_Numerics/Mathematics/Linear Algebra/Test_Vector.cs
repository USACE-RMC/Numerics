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

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;

namespace Mathematics.LinearAlgebra
{
    /// <summary>
    /// A class testing all vector operations.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Vector
    {
        /// <summary>
        /// Testing length of a vector
        /// </summary>
        [TestMethod()]
        public void Test_Length()
        {
            var a = new Vector(new[] { 1d, 2d, 3d });
            var length = a.Length;
            Assert.AreEqual(3, length);
        }

        /// <summary>
        /// Testing value at specified index
        /// </summary>
        [TestMethod()]
        public void Test_Index()
        {
            var a = new Vector(new[] { 1d, 2, 3 });
            var index = a[0];
            Assert.AreEqual(1, index);
        }

        /// <summary>
        /// Testing Vector to Array
        /// </summary>
        [TestMethod()]
        public void Test_ToArray()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var array = a.ToArray();
            var result = new [] { 1d, 2, 3d };
            for (int i = 0; i < array.Length; i++)
            {
                Assert.AreEqual(result[i], array[i]);
            }
        }

        /// <summary>
        /// Testing Vector to List
        /// </summary>
        [TestMethod()]
        public void Test_ToList()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var list = a.ToList();
            List<double> result = new List<double> { 1, 2, 3d };
            for (int i = 0; i<list.Count; i++)
            {
                Assert.AreEqual(list[i], result[i]);
            }
        }
        
        /// <summary>
        /// Cloning a vector
        /// </summary>
        [TestMethod()]
        public void Test_Clone()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var result = a.Clone();
            for (int i = 0;i<a.Length;i++)
            {
                Assert.AreEqual(a[i], result[i]);
            }
        }

        /// <summary>
        /// Testing norm of a vector
        /// </summary>
        [TestMethod()]
        public void Test_Norm()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var result = a.Norm();

            var true_result = Math.Sqrt(14);
            Assert.AreEqual(true_result, result);
        }
        
        /// <summary>
        /// Distance between two vectors
        /// </summary>
        [TestMethod()]
        public void Test_Distance()
        {
            //https://www.varsitytutors.com/calculus_3-help/distance-between-vectors
            var a = new Vector(new[] { 1, 0, 5d });
            var b = new Vector(new[] { 0, 2d, 4 });
            var result = Vector.Distance(a, b);
            var true_result = Math.Sqrt(6);
            Assert.AreEqual(true_result, result);
        }

        /// <summary>
        /// Dot product between two vectors
        /// </summary>
        [TestMethod()]
        public void Test_DotProduct()
        {
            var a = new Vector(new[] { 6, 2, -1d });
            var b = new Vector(new[] { 5, -8, 2d });
            var result = Vector.DotProduct(a, b);
            Assert.AreEqual(12, result);
        }


        /// <summary>
        /// Testing vector projection
        /// </summary>
        [TestMethod]
        public void Test_Projection()
        {
            // https://math.libretexts.org/Bookshelves/Applied_Mathematics/Mathematics_for_Game_Developers_(Burzynski)/02%3A_Vectors_In_Two_Dimensions/2.06%3A_The_Vector_Projection_of_One_Vector_onto_Another
            var u = new Vector(new[] { 4d, 3d });
            var v = new Vector(new[] { 2d, 8d });
            var vp = Vector.Project(u, v);
            Assert.AreEqual(16d / 17d, vp[0], 1E-6);
            Assert.AreEqual(64d / 17d, vp[1], 1E-6);
        }

        /// <summary>
        /// Testing multiplication operator between a scalar and a vector
        /// </summary>
        [TestMethod]
        public void Test_MultiplicationScalar()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var mult = a * 2;
            var result = new Vector(new[] { 2, 4, 6d });
            for (int i = 0; i < mult.Length; i++)
            {
                Assert.AreEqual(mult[i], result[i]);
            }
        }

        /// <summary>
        /// Testing multiplication operator between a vector and a vector
        /// </summary>
        [TestMethod]
        public void Test_MultiplicationVector()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var b = new Vector(new[] { 4, 5, 6d });
            var result = a * b;
            var true_result = new Vector(new[] { 4, 10, 18d });
            for (int i = 0; i < result.Length; i++)
            {
                Assert.AreEqual(result[i], true_result[i]);
            }
        }

        /// <summary>
        /// Raises vector to a power chosen by user
        /// </summary>
        [TestMethod]
        public  void Test_PowerOperator()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var pow = a ^ 2;
            var result = new Vector(new[] { 1, 4, 9d });
            for (int i = 0; i < result.Length; i++)
            {
                Assert.AreEqual(result[i], pow[i]);
            }
        }

        /// <summary>
        /// Adding two vectors
        /// </summary>
        [TestMethod]
        public void Test_AddOperator()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var b = new Vector(new[] { 4, 5, 6d });
            var result = a + b;
            var true_result = new Vector(new[] { 5, 7, 9d });
            for (int i = 0; i < result.Length; i++)
            {
                Assert.AreEqual(result[i], true_result[i]);
            }
        }

        /// <summary>
        /// Subtracting two vectors
        /// </summary>
        [TestMethod]
        public void Test_SubtractOperator()
        {
            var a = new Vector(new[] { 1, 2, 3d });
            var b = new Vector(new[] { 4, 5, 6d });
            var result = b - a;
            var true_result = new Vector(new[] { 3, 3, 3d });
            for (int i = 0; i < result.Length; i++)
            {
                Assert.AreEqual(result[i], true_result[i]);
            }
        }
    }
}
