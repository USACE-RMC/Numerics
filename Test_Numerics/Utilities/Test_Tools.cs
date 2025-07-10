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
using System.Collections.Generic;

namespace Utilities
{
    /// <summary>
    /// Testing the Tools class
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b> Authors: </b>
    /// <list type="bullet"> 
    /// <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Tools
    {
        /// <summary>
        /// Test Sign function with varying inputs.
        /// </summary>
        [TestMethod]
        public void Test_Sign()
        {
            double a = 1;
            double b = -2;
            var result = Tools.Sign(a, b);

            Assert.AreEqual(-1, result);

            double c = -1;
            double d = -2;
            var result2 = Tools.Sign(c, d);
            Assert.AreEqual(-1, result2);

            double e = -1;
            double f = 2;
            var result3 = Tools.Sign(e, f);
            Assert.AreEqual(1, result3);
        }

        /// <summary>
        /// Testing the squared value function.
        /// </summary>
        [TestMethod]
        public void Test_Sqr()
        {
            double a = 9;
            var result = Tools.Sqr(a);
            Assert.AreEqual(81, result);
        }

        /// <summary>
        /// Testing Pow function.
        /// </summary>
        [TestMethod]
        public void Test_Pow()
        {
            double a = 2;
            int b = 3;
            var result = Tools.Pow(a, b);
            Assert.AreEqual(8, result);

            double c = 123;
            int d = 0;
            var result2 = Tools.Pow(c, d);
            Assert.AreEqual(1, result2);
        }

        /// <summary>
        /// Testing Log10 function with different inputs.
        /// </summary>
        [TestMethod]
        public void Test_Log10()
        {
            double x = 100;
            var result = Tools.Log10(x);
            Assert.AreEqual(2, result);

            double x2 = 1;
            var result2 = Tools.Log10(x2);
            Assert.AreEqual(0, result2);
        }

        /// <summary>
        /// Testing natural log with different inputs.
        /// </summary>
        [TestMethod]
        public void Test_Log()
        {
            double x = 1;
            var result = Tools.Log(x);
            Assert.AreEqual(0, result);

            double x2 = 2.9;
            var result2 = Tools.Log(x2);
            Assert.AreEqual(1.0647, result2, 1E-04);
        }

        /// <summary>
        /// Testing Euclidean Distance between two points.
        /// </summary>
        [TestMethod]
        public void Test_Distance()
        {
            double x1 = 3;
            double y1 = 0;
            double x2 = 0;
            double y2 = 4;

            var result = Tools.Distance(x1, y1, x2, y2);
            Assert.AreEqual(5, result);

            List<double> point1 = new List<double> { 3d, 0 };
            List<double> point2 = new List<double> { 0, 4d };
            var result2 = Tools.Distance(point1, point2);
            Assert.AreEqual(5, result2);
        }

        /// <summary>
        /// Completing Min-Max normalization and denormalization on data
        /// </summary>
        [TestMethod]
        public void Test_Normalize_Denormalize()
        {
            var trueData = new double[] { 0, 250, 500, 750, 1000 };
            var trueNorm = new double[] { 0, 0.25, 0.5, 0.75, 1.0 };
            // Test normalization
            var dataNorm = Tools.Normalize(trueData);
            for (int i = 0; i < trueData.Length; i++)
            {
                Assert.AreEqual(trueNorm[i], dataNorm[i]);
            }
            // Test de-normalization
            var dataDenorm = Tools.Denormalize(dataNorm, Tools.Min(trueData), Tools.Max(trueData));
            for (int i = 0; i < trueData.Length; i++)
            {
                Assert.AreEqual(trueData[i], dataDenorm[i]);
            }
        }

        /// <summary>
        /// Standardizing (Z-score normalization) each value in a list.
        /// </summary>
        [TestMethod]
        public void Test_Standardize()
        {
            List<double> values = new List<double> { 3, 3, 4, 4, 6 };
            var result = Tools.Standardize(values);
            var true_result = new double[] { -0.8164, -0.8164, 0, 0, 1.63299 };
            for(int i = 0;i < values.Count; i++)
            {
                Assert.AreEqual(result[i], true_result[i], 1E-04);
            }
            
        }

        /// <summary>
        /// Destandardizing each value in a list.
        /// </summary>
        [TestMethod]
        public void Test_Destandardize()
        {
            List<double> values = new List<double> { -0.8164, -0.8164, 0, 0, 1.63299 };
            var result = Tools.Destandardize(values, 4, 1.224745);
            var true_result = new double[] { 3, 3, 4, 4, 6 };
            for (int i = 0; i < values.Count; i++)
            {
                Assert.AreEqual(result[i], true_result[i], 1E-03);
            }
        }

        /// <summary>
        /// Summing the int values in a list.
        /// </summary>
        [TestMethod]
        public void Test_SumInt()
        {
            List<int> values = new List<int> { 1, 2, 3 };
            var result = Tools.Sum(values);
            Assert.AreEqual(6, result);
        }

        /// <summary>
        /// Summing the double values in a list.
        /// </summary>
        [TestMethod]
        public void Test_SumDouble()
        {
            List<double> values = new List<double> { 1.4d, 2.3, 3.3 };
            var result = Tools.Sum(values);
            Assert.AreEqual(7, result);
        }

        /// <summary>
        /// Testing the sum of a list of values with indicator variables.
        /// </summary>
        [TestMethod]
        public void Test_SumIndicator()
        {
            List<double> values = new List<double> { 1.4, 2.3, 3.2, 3.3d };
            List<int> predictors = new List<int> { 1, 1, 1, 0 };
            var result = Tools.Sum(values,predictors);
            Assert.AreEqual(result, 6.9);
        }

        /// <summary>
        /// Testing the sum of the product of two lists' values.
        /// </summary>
        [TestMethod]
        public void Test_SumProduct()
        {
            List<double> list1 = new List<double> { 1, 2, 3 };
            List<double> list2 = new List<double> { 4, 5, 6 };

            var result = Tools.SumProduct(list1, list2);
            Assert.AreEqual(result, 32);
        }

        /// <summary>
        /// Testing the average of values in a list.
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.Mean(values);
            Assert.AreEqual(result, 2);
        }

        /// <summary>
        /// Testing average of values in a list with indicator variables.
        /// </summary>
        [TestMethod]
        public void Test_MeanIndicator()
        {
            List<double> values = new List<double> { 1, 2, 3,4 };
            List<int> indicators = new List<int> { 1, 1, 1, 0 };
            var result = Tools.Mean(values,indicators);
            Assert.AreEqual(result, 2);
        }

        /// <summary>
        /// Testing product of a list of values.
        /// </summary>
        [TestMethod]
        public void Test_Product()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.Product(values);
            Assert.AreEqual(result, 6);
        }

        /// <summary>
        /// Testing product of a list of values with indicator variables
        /// </summary>
        [TestMethod]
        public void Test_ProductIndicator()
        {
            List<double> values = new List<double> { 1, 2, 3,4 };
            List<int> indicators = new List<int> { 1,1,1,0 };
            var result = Tools.Product(values,indicators);
            Assert.AreEqual(result, 6);
        }

        /// <summary>
        /// Testing minimum value in a list
        /// </summary>
        [TestMethod]
        public void Test_Min()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.Min(values);
            Assert.AreEqual(result, 1);
        }

        /// <summary>
        /// Testing arg min value in a list
        /// </summary>
        [TestMethod]
        public void Test_ArgMin()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.ArgMin(values);
            Assert.AreEqual(result, 0);
        }

        /// <summary>
        /// Testing max value in a list
        /// </summary>
        [TestMethod]
        public void Test_Max()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.Max(values);
            Assert.AreEqual(result, 3);
        }


        /// <summary>
        /// Testing arg max value in a list
        /// </summary>
        [TestMethod]
        public void Test_ArgMax()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.ArgMax(values);
            Assert.AreEqual(result, 2);
        }

        /// <summary>
        /// Testing minimum value in a list with indicator variables
        /// </summary>
        [TestMethod]
        public void Test_MinIndicator()
        {
            List<double> values = new List<double> { 1, 2, 3, 4};
            List<int> indicators = new List<int> { 0, 1, 1, 1 };
            var result = Tools.Min(values,indicators);
            Assert.AreEqual(result, 2);
        }

        /// <summary>
        /// Testing maximum value in a list with indicator variables
        /// </summary>
        [TestMethod]
        public void Test_MaxIndicator()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            List<int> indicators = new List<int> { 1, 0, 0 };
            var result = Tools.Max(values, indicators);
            Assert.AreEqual(result, 1);
        }

        /// <summary>
        /// Testing the Log-sum-exponential function.
        /// <para> 
        /// <see href="https://www.rdocumentation.org/packages/matrixStats/versions/1.3.0/topics/logSumExp"/>
        /// </para>
        /// </summary>
        [TestMethod]
        public void Test_LogSumExp()
        {
            double u = 1000.01;
            double v = 1000.02;
            List<double> values = new List<double> { 1000.01, 1000.02 };

            var result = Tools.LogSumExp(u, v);
            var result2 = Tools.LogSumExp(values);
            Assert.AreEqual(result, 1000.70815, 1E-04);
            Assert.AreEqual(result2, 1000.70815, 1E-04);
        }

        /// <summary>
        /// Testing integer sequence.
        /// </summary>
        [TestMethod]
        public void Test_IntegerSequence()
        {
            int start = 0;
            int end = 3;
            var result = Tools.Sequence(start, end);
            var true_result = new int[] { 0, 1, 2, 3 };
            for (int i = 0;  i < true_result.Length; i++)
            {
                Assert.AreEqual(result[i], true_result[i]);
            }
        }

        /// <summary>
        /// Testing compress into byte array.
        /// </summary>
        [TestMethod]
        public void Test_Compress()
        {
            byte[] data = new byte[3];
            data[0] = 0;
            data[1] = 128;
            data[2] = 255;
            var result = Tools.Compress(data);
            Assert.IsTrue(data.Length <= result.Length);
        }

        /// <summary>
        /// Testing decompress a byte array.
        /// </summary>
        [TestMethod]
        public void Test_Decompress()
        {
            byte[] data = new byte[3];
            data[0] = 0;
            data[1] = 128;
            data[2] = 255;
            var result = Tools.Decompress(data);
            Assert.IsTrue(data.Length >= result.Length);
        }
    }
    
}
