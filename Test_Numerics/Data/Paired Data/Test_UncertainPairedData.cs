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
using Numerics.Data;
using Numerics.Distributions;
using System.Collections.Generic;
using System.Linq;
using System.Xml.Linq;

namespace Data.PairedData
{
    /// <summary>
    /// Unit tests for the UncertainOrderedPairedData class
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    ///     Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// </remarks>
    [TestClass]
    public class Test_UncertainPairedData
    {
        private readonly UncertainOrderedPairedData _dataset1;
        private readonly UncertainOrderedPairedData _dataset2;
        private readonly UncertainOrderedPairedData _dataset3;
        private readonly UncertainOrderedPairedData _dataset4;

        /// <summary>
        /// Creating a set of UncertainOrderedPairObjects to run tests on
        /// </summary>
        public Test_UncertainPairedData()
        {
            double[] xVals = new double[] { 1, 2, 3, 5 };
            UnivariateDistributionBase[] yVals = new UnivariateDistributionBase[] { new Triangular(1, 2, 3), new Triangular(2, 4, 5), new Triangular(6, 8, 12), new Triangular(13, 19, 20) };

            _dataset1 = new UncertainOrderedPairedData(xVals, yVals, true, SortOrder.Ascending, true, SortOrder.Ascending, UnivariateDistributionType.Triangular);
            _dataset2 = new UncertainOrderedPairedData(xVals.Reverse().ToArray(), yVals, true, SortOrder.Descending, true, SortOrder.Ascending, UnivariateDistributionType.Triangular);
            _dataset3 = new UncertainOrderedPairedData(xVals, yVals.Reverse().ToArray(), true, SortOrder.Ascending, true, SortOrder.Descending, UnivariateDistributionType.Triangular);
            _dataset4 = new UncertainOrderedPairedData(xVals.Reverse().ToArray(), yVals.Reverse().ToArray(), true, SortOrder.Descending, true, SortOrder.Descending, UnivariateDistributionType.Triangular);
        }

        /// <summary>
        /// Test method the SaveToXElement(), which saves the UncertainOrderedPairedData object as an XElement, 
        /// and test the conversion back to an OrderPairedData object from an XElement object
        /// </summary>
        [TestMethod]
        public void Test_ReadWriteXElement()
        {
            var el1 = _dataset1.SaveToXElement();
            var el2 = _dataset2.SaveToXElement();
            var el3 = _dataset3.SaveToXElement();
            var el4 = _dataset4.SaveToXElement();

            UncertainOrderedPairedData newDataset1 = new UncertainOrderedPairedData(el1);
            UncertainOrderedPairedData newDataset2 = new UncertainOrderedPairedData(el2);
            UncertainOrderedPairedData newDataset3 = new UncertainOrderedPairedData(el3);
            UncertainOrderedPairedData newDataset4 = new UncertainOrderedPairedData(el4);

            Assert.IsTrue(_dataset1 == newDataset1);
            Assert.IsTrue(_dataset2 == newDataset2);
            Assert.IsTrue(_dataset3 == newDataset3);
            Assert.IsTrue(_dataset4 == newDataset4);
        }

        /// <summary>
        /// Test the CurveSample() method
        /// </summary>
        [TestMethod]
        public void Test_CurveSample()
        {
            OrderedPairedData data1 = _dataset1.CurveSample();
            OrderedPairedData data2 = _dataset2.CurveSample();
            OrderedPairedData data3 = _dataset3.CurveSample();
            OrderedPairedData data4 = _dataset4.CurveSample();

            double[] xVals = new double[] { 1, 2, 3, 5 };
            // mean = (min + max + mode) / 3
            double[] yMeanVals = new double[] { 2, 3.66667, 8.66667, 17.33333 };

            var data1Expected = new OrderedPairedData(xVals, yMeanVals, true, SortOrder.Ascending, true, SortOrder.Ascending);
            var data2Expected = new OrderedPairedData(xVals.Reverse().ToArray(), yMeanVals, true, SortOrder.Descending, true, SortOrder.Ascending);
            var data3Expected = new OrderedPairedData(xVals, yMeanVals.Reverse().ToArray(), true, SortOrder.Ascending, true, SortOrder.Descending);
            var data4Expected = new OrderedPairedData(xVals.Reverse().ToArray(), yMeanVals.Reverse().ToArray(), true, SortOrder.Descending, true, SortOrder.Descending);

            for (int i = 0; i < data1.Count; i++)
            {
                Assert.AreEqual(data1Expected[i].X, data1[i].X, 1E-5);
                Assert.AreEqual(data1Expected[i].Y, data1[i].Y, 1E-5);

                Assert.AreEqual(data2Expected[i].X, data2[i].X, 1E-5);
                Assert.AreEqual(data2Expected[i].Y, data2[i].Y, 1E-5);

                Assert.AreEqual(data3Expected[i].X, data3[i].X, 1E-5);
                Assert.AreEqual(data3Expected[i].Y, data3[i].Y, 1E-5);

                Assert.AreEqual(data4Expected[i].X, data4[i].X, 1E-5);
                Assert.AreEqual(data4Expected[i].Y, data4[i].Y, 1E-5);
            }
        }

        /// <summary>
        /// Test the CurveSample() method that takes a probability sample. The output values were tested against R's "qtri()" function from the 
        /// "EnvStats" package which calculates the inverse cdf of the triangular distribution.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Millard SP (2013). EnvStats: An R Package for Environmental Statistics. Springer, New York.
        /// </remarks>
        [TestMethod]
        public void Test_Curve_Sample_Probability()
        {
            OrderedPairedData data1 = _dataset1.CurveSample(0.5);
            OrderedPairedData data2 = _dataset2.CurveSample(0.5);
            OrderedPairedData data3 = _dataset3.CurveSample(0.5);
            OrderedPairedData data4 = _dataset4.CurveSample(0.5);

            double[] xVals = new double[] { 1, 2, 3, 5 };
            // inverse cdf at the given probability (0.5)
            double[] yInverseVals = new double[] { 2, 3.732051, 8.535898, 17.58258 };

            var data1Expected = new OrderedPairedData(xVals, yInverseVals, true, SortOrder.Ascending, true, SortOrder.Ascending);
            var data2Expected = new OrderedPairedData(xVals.Reverse().ToArray(), yInverseVals, true, SortOrder.Descending, true, SortOrder.Ascending);
            var data3Expected = new OrderedPairedData(xVals, yInverseVals.Reverse().ToArray(), true, SortOrder.Ascending, true, SortOrder.Descending);
            var data4Expected = new OrderedPairedData(xVals.Reverse().ToArray(), yInverseVals.Reverse().ToArray(), true, SortOrder.Descending, true, SortOrder.Descending);

            for (int i = 0; i < data1.Count; i++)
            {
                Assert.AreEqual(data1Expected[i].X, data1[i].X, 1E-5);
                Assert.AreEqual(data1Expected[i].Y, data1[i].Y, 1E-5);

                Assert.AreEqual(data2Expected[i].X, data2[i].X, 1E-5);
                Assert.AreEqual(data2Expected[i].Y, data2[i].Y, 1E-5);

                Assert.AreEqual(data3Expected[i].X, data3[i].X, 1E-5);
                Assert.AreEqual(data3Expected[i].Y, data3[i].Y, 1E-5);

                Assert.AreEqual(data4Expected[i].X, data4[i].X, 1E-5);
                Assert.AreEqual(data4Expected[i].Y, data4[i].Y, 1E-5);
            }
        }

        /// <summary>
        /// Test the various UncertainOrderedPairedData IList methods
        /// </summary>
        [TestMethod]
        public void Test_IList()
        {
            var pairedData = _dataset1.Clone();
            var ordinate = pairedData[2];

            // Test IndexOf
            int test1 = pairedData.IndexOf(ordinate);
            Assert.AreEqual(2, test1);

            // Test Remove and Contains
            pairedData.Remove(ordinate);
            bool test2 = pairedData.Contains(ordinate);
            Assert.AreEqual(false, test2);

            // Test RemoveAt and Contains
            pairedData.RemoveAt(2);
            bool test3 = pairedData.Contains(ordinate);
            Assert.AreEqual(false, test3);

            // Test Insert and IndexOf
            pairedData.Insert(2, ordinate);
            int test4 = pairedData.IndexOf(ordinate);
            Assert.AreEqual(2, test4);

            // Test Add and IndexOf
            var newOrdinate = new UncertainOrdinate(7, new Triangular(16, 22, 28));
            pairedData.Add(newOrdinate);
            int test5 = pairedData.IndexOf(newOrdinate);
            Assert.AreEqual(pairedData.Count - 1, test5);

            // Test RemoveRange
            int prevCount = pairedData.Count();
            ordinate = pairedData[2];
            pairedData.RemoveRange(0, 2);
            UncertainOrdinate test6 = pairedData[0];
            int currCount = pairedData.Count();
            Assert.AreEqual(prevCount - 2, currCount);
            Assert.AreEqual(ordinate, test6);

            // Test CopyTo
            var array = new UncertainOrdinate[2];
            pairedData.CopyTo(array, 0);
            for (int i = 0; i < array.Length; i++)
            {
                Assert.AreEqual(pairedData[i], array[i]);
            }

            // Test InsertRange
            var toInsert = new List<UncertainOrdinate>() { new UncertainOrdinate(1, new Triangular(1, 2, 3)), new UncertainOrdinate(2, new Triangular(2, 4, 5)) };
            pairedData.InsertRange(0, toInsert);
            for (int j = 0; j < toInsert.Count; j++)
            {
                Assert.IsTrue(pairedData.Contains(toInsert[j]));
            }

            // Test AddRange
            var toAdd = new List<UncertainOrdinate>() { new UncertainOrdinate(3, new Triangular(6, 8, 12)), new UncertainOrdinate(5, new Triangular(13, 19, 20)), new UncertainOrdinate(7, new Triangular(16, 22, 28)) };
            pairedData.AddRange(toAdd);
            var test7 = pairedData.Count();
            Assert.AreEqual(7, test7);
            for (int k = 0; k < toAdd.Count(); k++)
            {
                Assert.IsTrue(pairedData.Contains(toAdd[k]));
            }

            // Test Clear
            pairedData.Clear();
            Assert.AreEqual(0, pairedData.Count());

        }

        /// <summary>
        /// Test the overloaded equality operators for the UncertainOrderedPairedData object
        /// </summary>
        [TestMethod]
        public void Test_Equality()
        {
            var dataset5 = _dataset1.Clone();
            bool test1 = (_dataset1 == dataset5);
            Assert.IsTrue(test1);

            bool test2 = (_dataset1 == _dataset2);
            Assert.IsFalse(test2);

            bool test3 = (_dataset2 != _dataset3);
            Assert.IsTrue(test3);

            bool test4 = (_dataset1 != dataset5);
            Assert.IsFalse(test4);
        }
    }
}
