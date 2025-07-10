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
using System.Xml.Linq;
using System.Globalization;
using System.Collections.Generic;


namespace Data.PairedData
{
    /// <summary>
    /// Unit tests for the Ordinate class.
    /// </summary>
    ///  <remarks>
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
    public class Test_Ordinate
    {
        /// <summary>
        /// Test the construction of the ordinate class
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var ordinate1 = new Ordinate(2, 4);

            var xElement = new XElement("xElement");
            xElement.SetAttributeValue("X", 2);
            xElement.SetAttributeValue("Y", 4);
            var ordinate2 = new Ordinate(xElement);

            var ordinate3 = new Ordinate(2, double.PositiveInfinity);
            var ordinate4 = new Ordinate(double.NaN, 4);

            Assert.AreEqual(ordinate1, ordinate2);
            Assert.AreEqual(ordinate1.X, 2);
            Assert.AreEqual(ordinate1.Y, 4);
            Assert.AreEqual(ordinate2.X, 2);
            Assert.AreEqual(ordinate2.Y, 4);
            
            Assert.AreEqual(ordinate1.IsValid, true);
            Assert.AreEqual(ordinate3.IsValid, false);
            Assert.AreEqual(ordinate4.IsValid, false);

            Assert.AreNotEqual(ordinate1, ordinate3);
            Assert.AreNotEqual(ordinate1, ordinate4);
        }

        /// <summary>
        /// Test the OrdinateValid() method that indicates if the ordinate being compare is valid given
        /// the criteria
        /// </summary>
        [TestMethod]
        public void Test_OrdinateValid()
        {
            var ordinate = new Ordinate(2, 4);
            var compareOrdinate1 = new Ordinate(3, 3);
            var compareOrdinate2 = new Ordinate(3, 5);
            var compareOrdinate3 = new Ordinate(1, 3);
            var compareOrdinate4 = new Ordinate(1, 5);

            var test1 = new bool[4];
            test1[0] = ordinate.OrdinateValid(compareOrdinate1, true, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test1[1] = ordinate.OrdinateValid(compareOrdinate2, true, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test1[2] = ordinate.OrdinateValid(compareOrdinate3, true, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test1[3] = ordinate.OrdinateValid(compareOrdinate4, true, true, SortOrder.Ascending, SortOrder.Ascending, true);

            var test1Expected = new bool[] { false, true, false, false };
            CollectionAssert.AreEqual(test1Expected, test1);

            var test2 = new bool[4];
            test2[0] = ordinate.OrdinateValid(compareOrdinate1, true, true, SortOrder.Ascending, SortOrder.Ascending, false);
            test2[1] = ordinate.OrdinateValid(compareOrdinate2, true, true, SortOrder.Ascending, SortOrder.Ascending, false);
            test2[2] = ordinate.OrdinateValid(compareOrdinate3, true, true, SortOrder.Ascending, SortOrder.Ascending, false);
            test2[3] = ordinate.OrdinateValid(compareOrdinate4, true, true, SortOrder.Ascending, SortOrder.Ascending, false);

            var test2Expected = new bool[] { false, false, true, false };
            CollectionAssert.AreEqual(test2Expected, test2);

            var test3 = new bool[4];
            test3[0] = ordinate.OrdinateValid(compareOrdinate1, false, true, SortOrder.None, SortOrder.Ascending, true);
            test3[1] = ordinate.OrdinateValid(compareOrdinate2, false, true, SortOrder.None, SortOrder.Ascending, true);
            test3[2] = ordinate.OrdinateValid(compareOrdinate3, false, true, SortOrder.None, SortOrder.Ascending, true);
            test3[3] = ordinate.OrdinateValid(compareOrdinate4, false, true, SortOrder.None, SortOrder.Ascending, true);

            var test3Expected = new bool[] { false, true, false, true };
            CollectionAssert.AreEqual(test3Expected, test3);

            var test4 = new bool[4];
            test4[0] = ordinate.OrdinateValid(compareOrdinate1, true, false, SortOrder.Ascending, SortOrder.None, true);
            test4[1] = ordinate.OrdinateValid(compareOrdinate2, true, false, SortOrder.Ascending, SortOrder.None, true);
            test4[2] = ordinate.OrdinateValid(compareOrdinate3, true, false, SortOrder.Ascending, SortOrder.None, true);
            test4[3] = ordinate.OrdinateValid(compareOrdinate4, true, false, SortOrder.Ascending, SortOrder.None, true);

            var test4Expected = new bool[] { true, true, false, false };
            CollectionAssert.AreEqual(test4Expected, test4);
            
            var test5 = new bool[4];
            test5[0] = ordinate.OrdinateValid(compareOrdinate1, false, false, SortOrder.Descending, SortOrder.Ascending, true);
            test5[1] = ordinate.OrdinateValid(compareOrdinate2, false, false, SortOrder.Descending, SortOrder.Ascending, true);
            test5[2] = ordinate.OrdinateValid(compareOrdinate3, false, false, SortOrder.Descending, SortOrder.Ascending, true);
            test5[3] = ordinate.OrdinateValid(compareOrdinate4, false, false, SortOrder.Descending, SortOrder.Ascending, true);

            var test5Expected = new bool[] { false, false, false, true };
            CollectionAssert.AreEqual(test5Expected, test5);

            var test6 = new bool[4];
            test6[0] = ordinate.OrdinateValid(compareOrdinate1, false, false, SortOrder.Ascending, SortOrder.Descending, true);
            test6[1] = ordinate.OrdinateValid(compareOrdinate2, false, false, SortOrder.Ascending, SortOrder.Descending, true);
            test6[2] = ordinate.OrdinateValid(compareOrdinate3, false, false, SortOrder.Ascending, SortOrder.Descending, true);
            test6[3] = ordinate.OrdinateValid(compareOrdinate4, false, false, SortOrder.Ascending, SortOrder.Descending, true);

            var test6Expected = new bool[] { true, false, false, false };
            CollectionAssert.AreEqual(test6Expected, test6);
        }

        /// <summary>
        /// Test the OrdinateErrors() method that gives error messages about the ordinate being compared
        /// given the criteria
        /// </summary>
        [TestMethod]
        public void Test_OrdinateErrors()
        {
            var ordinate = new Ordinate(2, 4);
            var compareOrdinate2 = new Ordinate(3, 5);
            var compareOrdinate3 = new Ordinate(1, 3);
            var compareOrdinate5 = new Ordinate(2, 4);

            var test1 = new List<string>();
            test1 = ordinate.OrdinateErrors(compareOrdinate5, true, true, SortOrder.Ascending, SortOrder.Descending, true);
            var test1Expected = new List<string>() { "Y values must be strictly decreasing.", "X values must be strictly increasing." };
            CollectionAssert.AreEqual(test1Expected, test1);

            var test2 = ordinate.OrdinateErrors(compareOrdinate2, false, false, SortOrder.Descending, SortOrder.Descending, true);
            var test2Expected = new List<string>() { "Y values must decrease.", "X values must decrease." };
            CollectionAssert.AreEqual(test2Expected, test2);

            var test3 = ordinate.OrdinateErrors(compareOrdinate3, false, false, SortOrder.Ascending, SortOrder.Ascending, true);
            var test3Expected = new List<string>() { "Y values must increase.", "X values must increase." };
            CollectionAssert.AreEqual(test3Expected, test3);

            var test4 = ordinate.OrdinateErrors(compareOrdinate3, false, false, SortOrder.Descending, SortOrder.Descending, false);
            var test4Expected = new List<string>() { "Y values must decrease.", "X values must decrease." };
            CollectionAssert.AreEqual(test4Expected, test4);

            var test5 = ordinate.OrdinateErrors(compareOrdinate2, false, false, SortOrder.Ascending, SortOrder.Ascending, false);
            var test5Expected = new List<string>() { "Y values must increase.", "X values must increase." };
            CollectionAssert.AreEqual(test5Expected, test5);
        }

        /// <summary>
        /// Test the transform methods in the class
        /// </summary>
        [TestMethod]
        public void Test_Transform()
        {
            var original = new Ordinate(50, 100);
            var logX = original.Transform(Transform.Logarithmic, Transform.None);
            var logY = original.Transform(Transform.None, Transform.Logarithmic);
            var logXY = original.Transform(Transform.Logarithmic, Transform.Logarithmic);

            Assert.AreEqual(1.69897, logX.X, 1E-6);
            Assert.AreEqual(2, logY.Y, 1E-6);
            Assert.AreEqual(1.69897, logXY.X, 1E-6);
            Assert.AreEqual(2, logXY.Y, 1E-6);

            var originalZ = new Ordinate(0.3, 0.8);
            var zX = originalZ.Transform(Transform.NormalZ, Transform.None);
            var zY = originalZ.Transform(Transform.None, Transform.NormalZ);
            var zXY = originalZ.Transform(Transform.NormalZ, Transform.NormalZ);

            Assert.AreEqual(-0.5244005, zX.X, 1E-6);
            Assert.AreEqual(0.8416212, zY.Y, 1E-6);
            Assert.AreEqual(-0.5244005, zXY.X, 1E-6);
            Assert.AreEqual(0.8416212, zXY.Y, 1E-6);
        }

        /// <summary>
        /// Test the overloaded equality (==) and inequality (!=) operators
        /// </summary>
        [TestMethod]
        public void Test_Equality()
        {
            var ordinate1 = new Ordinate(7.325, 6.389);
            var ordinate2 = new Ordinate(7.325, 6.389);
            var ordinate3 = new Ordinate(8.36, 25.99);

            var test = new bool[4];
            test[0] = (ordinate1 == ordinate2);
            test[1] = (ordinate1 == ordinate3);
            test[2] = (ordinate2 != ordinate3);
            test[3] = (ordinate2 != ordinate1);

            var expected = new bool[] { true, false, true, false };
            CollectionAssert.AreEqual(expected, test);
        }

        /// <summary>
        /// Test the method that converts an Ordinate object to an XElement object
        /// </summary>
        [TestMethod]
        public void Test_ToXElement()
        {
            double X = 4.2;
            double Y = 3.869;
            var original = new Ordinate(X, Y);
            XElement xElement = original.ToXElement();
            double x, y;
            double.TryParse(xElement.Attribute(nameof(X)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out x);
            double.TryParse(xElement.Attribute(nameof(Y)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out y);

            Assert.AreEqual(X, x);
            Assert.AreEqual(Y, y);
        }
    }
}
