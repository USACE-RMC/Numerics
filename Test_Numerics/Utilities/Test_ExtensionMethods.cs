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
using Numerics;
using System;
using System.Linq;
using Numerics.Data;
using Numerics.Mathematics.LinearAlgebra;

namespace Utilities
{
    /// <summary>
    /// Testing the ExtensionMethods class
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
    [TestClass]
    public class Test_ExtensionMethods
    {

        /// <summary>
        /// Testing get attribute with the Description Attribute example.
        /// </summary>
        [TestMethod]
        public void Test_GetAttribute()
        {
            var expectedDescription = "1-Hr";
            var attribute = TimeInterval.OneHour.GetAttributeOfType<System.ComponentModel.DescriptionAttribute>();
            Assert.AreEqual(expectedDescription, attribute.Description);
        }

        /// <summary>
        /// Testing almost equals with a false return.
        /// </summary>
        [TestMethod]
        public void Test_AlmostEquals()
        {
            var a = 2;
            var b = 1;
            var result = ExtensionMethods.AlmostEquals(a, b);
            Assert.IsFalse(result);
        }

        /// <summary>
        /// Testing a random integer array generator. Should output 5 random integers. 
        /// </summary>
        [TestMethod]
        public void Test_NextIntegers()
        {
            var random = new Random();
            var result = random.NextIntegers(5);
            Assert.AreEqual(5, result.Length);
        }

        /// <summary>
        /// Should output 5 random integers between 0 and 10.
        /// </summary>
        [TestMethod]
        public void Test_NextIntegersMinMax()
        {
            var random = new Random();
            var result = random.NextIntegers(0, 10, 5);
            Assert.AreEqual(5, result.Length);
            for (int i = 0; i < result.Length; i++)
                Assert.IsTrue(result[i] >= 0 && result[i] < 10);
        }

        /// <summary>
        /// Test next with replacement.
        /// </summary>
        [TestMethod]
        public void Test_NextIntegers_WithReplacement()
        {
            var random = new Random(12345);
            var result = random.NextIntegers(0, 1000, 1000, true);
            Assert.AreNotEqual(1000, result.Distinct().Count());
            for (int i = 0; i < result.Length; i++)
                Assert.IsTrue(result[i] >= 0 && result[i] < 1000);
        }

        /// <summary>
        /// Test next integers without replacement.
        /// </summary>
        [TestMethod]
        public void Test_NextIntegers_WithoutReplacement()
        {
            var random = new Random(12345);
            var result = random.NextIntegers(0, 1000, 1000, false);
            Assert.AreEqual(1000, result.Distinct().Count());
            for (int i = 0; i < result.Length; i++)
                Assert.IsTrue(result[i] >= 0 && result[i] < 1000);
        }

        /// <summary>
        /// Test next doubles. 
        /// </summary>
        [TestMethod]
        public void Test_NextDoubles()
        {
            var random = new Random();
            var result = random.NextDoubles(5);
            Assert.AreEqual(5, result.Length);
            for (int i = 0; i < result.Length; i++)
                Assert.IsTrue(result[i] >= 0 && result[i] < 1);
        }

        /// <summary>
        /// Test next doubles in multiple dimensions 
        /// </summary>
        [TestMethod]
        public void Test_NextDoubles2D()
        {
            var random = new Random();
            var result = random.NextDoubles(5, 3);
            for (int i = 0; i < result.GetLength(0); i++)
            {
                for (int j = 0; j < result.GetLength(1); j++)
                {
                    Assert.IsTrue(result[i, j] >= 0 && result[i, j] < 1);
                }
            }
        }


        /// <summary>
        /// Test addition of two arrays.
        /// </summary>
        [TestMethod]
        public void Test_Add()
        {
            double[] array = new double[] { 1, 2, 3 };
            double[] values = new double[] { 1, 2, 3 };
            var trueResult = new double[] { 2, 4, 6 };
            var result = array.Add(values);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Test edge case of different sized arrays throwing exception.
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Test_AddEdge()
        {
            double[] array = new double[] { 1, 2, 3 };
            double[] values = new double[] { 1, 2, 3, 4};
            var result = array.Add(values);
        }

        /// <summary>
        /// Test subtraction of two arrays.
        /// </summary>
        [TestMethod]
        public void Test_Subtract()
        {
            var array = new double[] { 3, 2, 1 };
            var values = new double[] { 1, 2, 3 };
            var trueResult = new double[] { 2, 0, -2 };
            var result = array.Subtract(values);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Test edge case of different sized arrays throwing exception.
        /// </summary>
        [TestMethod]
        [ExpectedException(typeof(ArgumentException))]
        public void Test_SubtractEdge()
        {
            var array = new double[] { 3, 2, 1, 0 };
            var values = new double[] { 1, 2, 3 };
            var result = array.Subtract(values);       
        }

        /// <summary>
        /// Test multiplying an array by a scalar.
        /// </summary>
        [TestMethod]
        public void Test_Multiply()
        {
            var scalar = 2;
            var array = new double[] { 4, 5, 6, 7 };
            var trueResult = new double[] { 8, 10, 12, 14 };
            var result = array.Multiply(scalar);
            CollectionAssert.AreEqual(trueResult, result);
        }
        
        /// <summary>
        /// Test dividing an array by a scalar. 
        /// </summary>
        [TestMethod]
        public void Test_Divide()
        {
            var scalar = 2;
            var array = new double[] { 4, -5, 6, 7 };
            var trueResult = new double[] { 2, -2.5, 3, 3.5 };
            var result = array.Divide(scalar);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a subset of an array. 
        /// </summary>
        [TestMethod]
        public void Test_Subset_Start()
        {
            var startIndex = 3;
            var array = new double[] { 1, 2, 3, 4, 5, 6 };
            var trueResult = new double[] { 4, 5, 6 };
            var result = array.Subset(startIndex);
            CollectionAssert.IsSubsetOf(result, array);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a subset of an array. 
        /// </summary>
        [TestMethod]
        public void Test_Subset_Start_End()
        {
            var startIndex = 0;
            var endIndex = 2;
            var array = new double[] { 1, 2, 3, 4, 5, 6 };
            var trueResult = new double[] { 1, 2, 3 };
            var result = array.Subset(startIndex, endIndex);
            CollectionAssert.IsSubsetOf(result, array);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a random subset of an array. 
        /// </summary>
        [TestMethod]
        public void Test_Random_Subset()
        {
            var array = new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
            var result = array.RandomSubset(5);          
            CollectionAssert.IsSubsetOf(result, array);
            Assert.AreEqual(5, result.Length);
        }

        /// <summary>
        /// Testing the fill function to see if the array was populated with the specified value.
        /// </summary>
        [TestMethod]
        public void Test_Fill()
        {
            var array = new double[3];
            var fillValue = 1.2;
            array.Fill(fillValue);
            for (int i = 0; i < array.Length; i++)
                Assert.AreEqual(fillValue, array[i]);
        }

        /// <summary>
        /// Test GetColumn from 2D array.
        /// </summary>
        [TestMethod]
        public void Test_GetColumn()
        {
            int[,] array = new int[,]
            {
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 }
            };
            int expectedColumnIndex = 1;
            int[] expectedColumn = new int[] { 2, 5, 8 };
            int[] result = array.GetColumn(expectedColumnIndex);
            CollectionAssert.AreEqual(expectedColumn, result);
        }

        /// <summary>
        /// Test GetRow from 2D array.
        /// </summary>
        [TestMethod]
        public void Test_GetRow()
        {
            int[,] array = new int[,]
            {
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 }
            };
            int expectedRowIndex = 1;
            int[] expectedColumn = new int[] { 4, 5, 6 };
            int[] result = array.GetRow(expectedRowIndex);
            CollectionAssert.AreEqual(expectedColumn, result);
        }

        /// <summary>
        /// Testing a subset of a 2D array. 
        /// </summary>
        [TestMethod]
        public void Test_Subset_Start_2D()
        {
            var startIndex = 2;
            var array = new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 }, { 7, 8 }, { 9, 10 } };
            var trueResult = new double[,] { { 5, 6 }, { 7, 8 }, { 9, 10 } };
            var result = array.Subset(startIndex);
            CollectionAssert.IsSubsetOf(result, array);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a subset of a 2D array. 
        /// </summary>
        [TestMethod]
        public void Test_Subset_Start_End_2D()
        {
            var startIndex = 2;
            var endIndex = 3;
            var array = new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 }, { 7, 8 }, { 9, 10 } };
            var trueResult = new double[,] { { 5, 6 }, { 7, 8 } };
            var result = array.Subset(startIndex, endIndex);
            CollectionAssert.IsSubsetOf(result, array);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a random subset of a 2D array. 
        /// </summary>
        [TestMethod]
        public void Test_Random_Subset_2D()
        {
            var array = new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 }, { 7, 8 }, { 9, 10 } };
            var result = array.RandomSubset(3);
            CollectionAssert.IsSubsetOf(result, array);
            Assert.AreEqual(3, result.GetLength(0));
        }

        /// <summary>
        /// Testing the fill function to see if the 2D array was populated with the specified value.
        /// </summary>
        [TestMethod]
        public void Test_Fill_2D()
        {
            var array = new double[3, 3];
            var fillValue = 1.2;
            array.Fill(fillValue);
            for (int i = 0; i < array.GetLength(0); i++)
                for (int j = 0; j < array.GetLength(1); j++)
                    Assert.AreEqual(fillValue, array[i, j]);
        }

        /// <summary>
        /// Testing a subset of a vector. 
        /// </summary>
        [TestMethod]
        public void Test_Vector_Subset_Start()
        {
            var startIndex = 3;
            var vector = new Vector(new double[] { 1, 2, 3, 4, 5, 6 });
            var trueResult = new double[] { 4, 5, 6 };
            var result = vector.Subset(startIndex).Array;
            CollectionAssert.IsSubsetOf(result, vector.Array);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a subset of a vector. 
        /// </summary>
        [TestMethod]
        public void Test_Vector_Subset_Start_End()
        {
            var startIndex = 0;
            var endIndex = 2;
            var vector = new Vector(new double[] { 1, 2, 3, 4, 5, 6 });
            var trueResult = new double[] { 1, 2, 3 };
            var result = vector.Subset(startIndex, endIndex).Array;
            CollectionAssert.IsSubsetOf(result, vector.Array);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a random subset of a vector. 
        /// </summary>
        [TestMethod]
        public void Test_Vector_Random_Subset()
        {
            var vector = new Vector(new double[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 });
            var result = vector.RandomSubset(5).Array;
            CollectionAssert.IsSubsetOf(result, vector.Array);
            Assert.AreEqual(5, result.Length);
        }

        /// <summary>
        /// Testing the fill function to see if the vector was populated with the specified value.
        /// </summary>
        [TestMethod]
        public void Test_Vector_Fill()
        {
            var vector = new Vector(3);
            var fillValue = 1.2;
            vector.Fill(fillValue);
            for (int i = 0; i < vector.Length; i++)
                Assert.AreEqual(fillValue, vector[i]);
        }

        /// <summary>
        /// Test GetColumn from matrix.
        /// </summary>
        [TestMethod]
        public void Test_Matrix_GetColumn()
        {
            var matrix = new Matrix(new double[,]
            {
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 }
            });
            int expectedColumnIndex = 1;
            double[] expectedColumn = new double[] { 2, 5, 8 };
            double[] result = matrix.GetColumn(expectedColumnIndex);
            CollectionAssert.AreEqual(expectedColumn, result);
        }

        /// <summary>
        /// Test GetRow from matrix.
        /// </summary>
        [TestMethod]
        public void Test_Matrix_GetRow()
        {
            var matrix = new Matrix(new double[,]
                        {
                { 1, 2, 3 },
                { 4, 5, 6 },
                { 7, 8, 9 }
                        });
            int expectedRowIndex = 1;
            double[] expectedColumn = new double[] { 4, 5, 6 };
            double[] result = matrix.GetRow(expectedRowIndex);
            CollectionAssert.AreEqual(expectedColumn, result);
        }

        /// <summary>
        /// Testing a subset of a matrix. 
        /// </summary>
        [TestMethod]
        public void Test_Subset_Start_Matrix()
        {
            var startIndex = 2;
            var matrix = new Matrix(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 }, { 7, 8 }, { 9, 10 } });
            var trueResult = new double[,] { { 5, 6 }, { 7, 8 }, { 9, 10 } };
            var result = matrix.Subset(startIndex).Array;
            CollectionAssert.IsSubsetOf(result, matrix.Array);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a subset of a matrix. 
        /// </summary>
        [TestMethod]
        public void Test_Subset_Start_End_Matrix()
        {
            var startIndex = 2;
            var endIndex = 3;
            var matrix = new Matrix(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 }, { 7, 8 }, { 9, 10 } });
            var trueResult = new double[,] { { 5, 6 }, { 7, 8 } };
            var result = matrix.Subset(startIndex, endIndex).Array;
            CollectionAssert.IsSubsetOf(result, matrix.Array);
            CollectionAssert.AreEqual(trueResult, result);
        }

        /// <summary>
        /// Testing a random subset of a matrix. 
        /// </summary>
        [TestMethod]
        public void Test_Random_Subset_Matrix()
        {
            var matrix = new Matrix(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 }, { 7, 8 }, { 9, 10 } });
            var result = matrix.RandomSubset(3).Array;
            CollectionAssert.IsSubsetOf(result, matrix.Array);
            Assert.AreEqual(3, result.GetLength(0));
        }

        /// <summary>
        /// Testing the fill function to see if the matrix was populated with the specified value.
        /// </summary>
        [TestMethod]
        public void Test_Fill_Matrix()
        {
            var matrix = new Matrix(new double[3, 3]);
            var fillValue = 1.2;
            matrix.Fill(fillValue);
            for (int i = 0; i < matrix.NumberOfRows; i++)
                for (int j = 0; j < matrix.NumberOfColumns; j++)
                    Assert.AreEqual(fillValue, matrix[i, j]);
        }

    }
}
