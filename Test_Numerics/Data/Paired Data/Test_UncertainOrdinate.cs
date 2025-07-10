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
using System.Xml.Linq;
using System.Globalization;
using System.Collections.Generic;

namespace Data.PairedData
{
    /// <summary>
    /// Unit tests for the UncertainOrdinate class
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
    public class Test_UncertainOrdinate
    {
        /// <summary>
        /// Test the construction of the uncertain ordinate class
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var distribution = new Normal();
            var unordinate1 = new UncertainOrdinate(2, distribution);

            var xElement = new XElement("xElement");
            xElement.SetAttributeValue("X", 2);
            xElement.SetAttributeValue("Mu", 0);
            xElement.SetAttributeValue("Sigma", 1);

            var unordinate2 = new UncertainOrdinate(xElement, UnivariateDistributionType.Normal);

            var xElement2 = new XElement("xElement2");
            xElement2.SetElementValue("Distribution", UnivariateDistributionType.Normal);
            xElement2.SetAttributeValue("X", 2);
            xElement2.Element("Distribution").SetAttributeValue("Type", UnivariateDistributionType.Normal);
            xElement2.Element("Distribution").SetAttributeValue("Mu", 0);
            xElement2.Element("Distribution").SetAttributeValue("Sigma", 1);

            var unordinate3 = new UncertainOrdinate(xElement2);

            var unordinate4 = new UncertainOrdinate(double.NegativeInfinity, distribution);
            var unordinate5 = new UncertainOrdinate(3, null);

            // Also testing overloaded equality operators
            Assert.IsTrue(unordinate1 == unordinate3);
            Assert.IsTrue(unordinate1 == unordinate2);
            Assert.AreEqual(unordinate1.X, 2);
            Assert.AreEqual(unordinate1.Y, distribution);
            Assert.AreEqual(unordinate3.X, 2);
            Assert.IsTrue(unordinate3.Y == distribution);

            Assert.AreEqual(unordinate1.IsValid, true);
            Assert.AreEqual(unordinate4.IsValid, false);
            Assert.AreEqual(unordinate5.IsValid, false);

            Assert.IsTrue(unordinate1 != unordinate4);
            Assert.IsTrue(unordinate1 != unordinate5);
        }

        /// <summary>
        /// Test the OrdinateValid() method that indicates if the ordinate being compare is valid given
        /// the criteria
        /// </summary>
        [TestMethod]
        public void Test_OrdinateValid()
        {
            var unordinate = new UncertainOrdinate(3, new Triangular(6, 8, 12));
            var compareUnordinate1 = new UncertainOrdinate(0, new Triangular(1, 2, 3));
            var compareUnordinate2 = new UncertainOrdinate(2, new Triangular(2, 4, 5));
            var compareUnordinate3 = new UncertainOrdinate(5, new Triangular(13, 19, 20));
            var compareUnordinate4 = new UncertainOrdinate(6, new Uniform(7,14));

            var test1 = new bool[4];
            test1[0] = unordinate.OrdinateValid(compareUnordinate1, false, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test1[1] = unordinate.OrdinateValid(compareUnordinate2, false, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test1[2] = unordinate.OrdinateValid(compareUnordinate3, false, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test1[3] = unordinate.OrdinateValid(compareUnordinate4, false, true, SortOrder.Ascending, SortOrder.Ascending, true);

            var test1Expected = new bool[] { false, false, true, false };
            CollectionAssert.AreEqual(test1Expected, test1);

            var test2 = new bool[4];
            test2[0] = unordinate.OrdinateValid(compareUnordinate1, false, true, SortOrder.Ascending, SortOrder.Ascending, false);
            test2[1] = unordinate.OrdinateValid(compareUnordinate2, false, true, SortOrder.Ascending, SortOrder.Ascending, false);
            test2[2] = unordinate.OrdinateValid(compareUnordinate3, false, true, SortOrder.Ascending, SortOrder.Ascending, false);
            test2[3] = unordinate.OrdinateValid(compareUnordinate4, false, true, SortOrder.Ascending, SortOrder.Ascending, false);

            var test2Expected = new bool[] { true, true, false, false };
            CollectionAssert.AreEqual(test2Expected, test2);

            var test3 = new bool[4];
            test3[0] = unordinate.OrdinateValid(compareUnordinate1, true, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test3[1] = unordinate.OrdinateValid(compareUnordinate2, true, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test3[2] = unordinate.OrdinateValid(compareUnordinate3, true, true, SortOrder.Ascending, SortOrder.Ascending, true);
            test3[3] = unordinate.OrdinateValid(compareUnordinate4, true, true, SortOrder.Ascending, SortOrder.Ascending, true);

            var test3Expected = new bool[] { false, false, true, false };
            CollectionAssert.AreEqual(test3Expected, test3);

            var test4 = new bool[4];
            test4[0] = unordinate.OrdinateValid(compareUnordinate1, false, false, SortOrder.Ascending, SortOrder.Ascending, true);
            test4[1] = unordinate.OrdinateValid(compareUnordinate2, false, false, SortOrder.Ascending, SortOrder.Ascending, true);
            test4[2] = unordinate.OrdinateValid(compareUnordinate3, false, false, SortOrder.Ascending, SortOrder.Ascending, true);
            test4[3] = unordinate.OrdinateValid(compareUnordinate4, false, false, SortOrder.Ascending, SortOrder.Ascending, true);

            var test4Expected = new bool[] { false, false, true, false };
            CollectionAssert.AreEqual(test4Expected, test4);

            var test5 = new bool[4];
            test5[0] = unordinate.OrdinateValid(compareUnordinate1, false, true, SortOrder.Descending, SortOrder.Ascending, true);
            test5[1] = unordinate.OrdinateValid(compareUnordinate2, false, true, SortOrder.Descending, SortOrder.Ascending, true);
            test5[2] = unordinate.OrdinateValid(compareUnordinate3, false, true, SortOrder.Descending, SortOrder.Ascending, true);
            test5[3] = unordinate.OrdinateValid(compareUnordinate4, false, true, SortOrder.Descending, SortOrder.Ascending, true);

            var test5Expected = new bool[] { false, false, false, false };
            CollectionAssert.AreEqual(test5Expected, test5);

            var test6 = new bool[4];
            test6[0] = unordinate.OrdinateValid(compareUnordinate1, false, true, SortOrder.Ascending, SortOrder.Descending, true); 
            test6[1] = unordinate.OrdinateValid(compareUnordinate2, false, true, SortOrder.Ascending, SortOrder.Descending, true);
            test6[2] = unordinate.OrdinateValid(compareUnordinate3, false, true, SortOrder.Ascending, SortOrder.Descending, true);
            test6[3] = unordinate.OrdinateValid(compareUnordinate4, false, true, SortOrder.Ascending, SortOrder.Descending, true);

            var test6Expected = new bool[] { false, false, false, false };
            CollectionAssert.AreEqual(test6Expected, test6);

            var test7 = new bool[4];
            test7[0] = unordinate.OrdinateValid(compareUnordinate1, false, true, SortOrder.Descending, SortOrder.Descending, true);
            test7[1] = unordinate.OrdinateValid(compareUnordinate2, false, true, SortOrder.Descending, SortOrder.Descending, true);
            test7[2] = unordinate.OrdinateValid(compareUnordinate3, false, true, SortOrder.Descending, SortOrder.Descending, true);
            test7[3] = unordinate.OrdinateValid(compareUnordinate4, false, true, SortOrder.Descending, SortOrder.Descending, true);

            var test7Expected = new bool[] { true, true, false, false };
            CollectionAssert.AreEqual(test7Expected, test7);

            var test8 = new bool[4];
            test8[0] = unordinate.OrdinateValid(compareUnordinate1, false, true, SortOrder.Ascending, SortOrder.Ascending, true, true);
            test8[1] = unordinate.OrdinateValid(compareUnordinate2, false, true, SortOrder.Ascending, SortOrder.Ascending, true, true);
            test8[2] = unordinate.OrdinateValid(compareUnordinate3, false, true, SortOrder.Ascending, SortOrder.Ascending, true, true);
            test8[3] = unordinate.OrdinateValid(compareUnordinate4, false, true, SortOrder.Ascending, SortOrder.Ascending, true, true);

            var test8Expected = new bool[] { false, false, true, true };
            CollectionAssert.AreEqual(test8Expected, test8);
        }

        /// <summary>
        /// Test the OrdinateErrors() method that gives error messages about the ordinate being compared
        /// given the criteria
        /// </summary>
        [TestMethod]
        public void Test_OrdinateErrors()
        {

            var unordinate = new UncertainOrdinate(3, new Triangular(6, 8, 12));
            var compareUnordinate1 = new UncertainOrdinate(2, new Triangular(2, 4, 5));
            var compareUnordinate2= new UncertainOrdinate(5, new Triangular(13, 19, 20));
            var compareUnordinate3 = new UncertainOrdinate(double.PositiveInfinity, new Uniform(7, 14));

            var test1 = unordinate.OrdinateErrors( compareUnordinate3, true, true, SortOrder.Ascending, SortOrder.Ascending, true);
            var test1Expected = new List<string>() { "Ordinate X value can not be infinity.", "Can't compare two ordinates with different distribution types." };
            CollectionAssert.AreEqual (test1Expected, test1);

            var test2 = unordinate.OrdinateErrors(compareUnordinate1, true, true, SortOrder.Ascending, SortOrder.Ascending, true);
            // Each of the three tests (lower bound, central tendency, upper bound) will give the same error
            var test2Expected = new List<string>() { "Y values must increase.", "X values must increase.", "Y values must increase.", "X values must increase.", "Y values must increase.", "X values must increase." };
            CollectionAssert.AreEqual(test2Expected, test2);

            var test3 = unordinate.OrdinateErrors(compareUnordinate2, true, true, SortOrder.Descending, SortOrder.Descending, true);
            // Each of the three tests (lower bound, central tendency, upper bound) will give the same error
            var test3Expected = new List<string>() { "Y values must decrease.", "X values must decrease.", "Y values must decrease.", "X values must decrease.", "Y values must decrease.", "X values must decrease." };
            CollectionAssert.AreEqual(test3Expected, test3);
                
        }

        /// <summary>
        /// Test the method that converts an UncertainOrdinate object to an XElement object
        /// </summary>
        [TestMethod]
        public void Test_ToXElement()
        {
            double X = 2;
            var distribution = new Normal();
            var orginal = new UncertainOrdinate(X, distribution);
            XElement xElement = orginal.ToXElement();
            double x = 0;
            UnivariateDistributionBase dist = null;
            if (xElement.Attribute(nameof(X)) != null) double.TryParse(xElement.Attribute(nameof(X)).Value, NumberStyles.Any, CultureInfo.InvariantCulture, out x);
            if (xElement.Element("Distribution") != null) { dist = UnivariateDistributionFactory.CreateDistribution(xElement.Element("Distribution")); }

            Assert.AreEqual(X, x);
            Assert.IsTrue(distribution == dist);
        }
    }
}
