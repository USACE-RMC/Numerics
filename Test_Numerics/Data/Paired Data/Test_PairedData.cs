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

using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;

namespace Data.PairedData
{
    /// <summary>
    /// Unit tests for the OrderedPairedData class
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
    public class Test_PairedData
    {
        private readonly OrderedPairedData _dataset1;
        private readonly OrderedPairedData _dataset2;
        private readonly OrderedPairedData _dataset3;
        private readonly OrderedPairedData _dataset4;

        /// <summary>
        /// Creating a set of OrderedPairObjects to run tests on
        /// </summary>
        public Test_PairedData()
        {
            double[] xVals = new double[] { 230408, 288010, 345611, 403213, 460815, 518417, 576019, 633612, 691223, 748825, 806427, 864029, 921631, 1036834, 1152038 };
            double[] yVals = new double[] { 1519.7, 1520.5, 1520.9, 1521.7, 1523.5, 1525.9, 1528.4, 1530.9, 1533.2, 1534.7, 1535.9, 1538, 1541.3, 1547.7, 1552.7 };
            _dataset1 = new OrderedPairedData(xVals, yVals, true, SortOrder.Ascending, true, SortOrder.Ascending);
            _dataset2 = new OrderedPairedData(xVals.Reverse().ToArray(), yVals, true, SortOrder.Descending, true, SortOrder.Ascending);
            _dataset3 = new OrderedPairedData(xVals, yVals.Reverse().ToArray(), true, SortOrder.Ascending, true, SortOrder.Descending);
            _dataset4 = new OrderedPairedData(xVals.Reverse().ToArray(), yVals.Reverse().ToArray(), true, SortOrder.Descending, true, SortOrder.Descending);
        }

        /// <summary>
        /// Helper function to test the GetYfromX and GetYValues function. This method calls the GetYfromX on each of the given
        /// samples, and the GetYValues on the array of samples, and compares the results to the expected values.
        /// </summary>
        /// <param name="data"> The OrderedPairedData object being tested </param>
        /// <param name="samples"> The array of values we want to test (i.e. pass to the GetYfromX method) </param>
        /// <param name="expectedResults"> The array of expected Y's for the given X samples </param>
        private void SampleDatasetY(OrderedPairedData data, double[] samples, double[] expectedResults)
        {
            // Test passing a single value
            for (int i = 0; i < samples.Count(); i++)
                Assert.AreEqual(expectedResults[i], data.GetYFromX(samples[i]), 1E-5);

            // Test passing a list
            double[] xFromY = data.GetYFromX(samples);
            for (int i = 0; i < samples.Count(); i++)
            {
                Assert.AreEqual(expectedResults[i], xFromY[i], 1E-5);
            }
        }

        /// <summary>
        /// Helper function to test the GetYfromX function. This method calls the GetXfromY on each of the given
        /// samples and compares the result to the expected value.
        /// </summary>
        /// <param name="data"> The OrderedPairedData object being tested </param>
        /// <param name="samples"> The array of values we want to test (i.e. pass to the GetXfromY method) </param>
        /// <param name="expectedResults"> The array of expected X's for the given Y samples </param>
        private void SampleDatasetX(OrderedPairedData data, double[] samples, double[] expectedResults)
        {
            // Test passing a single value
            for (int i = 0; i < samples.Count(); i++)
                Assert.AreEqual(expectedResults[i], data.GetXFromY(samples[i]), 1E-5);

            // Test passing a list
            double[] xFromY = data.GetXFromY(samples);
            for (int i = 0; i < samples.Count(); i++)
            {
                Assert.AreEqual(expectedResults[i], xFromY[i], 1E-5);
            }
        }

        /// <summary>
        /// Test the methods that get Y from given X values
        /// </summary>
        [TestMethod()]
        public void Test_GetY()
        {
            double[] xVals = new double[] { 1018627, 742619, 770076, 350167, 260164, 502421, 1034555, 810438, 655158, 253951, 1149424, 973525, 397450, 1128872, 640330, 303494, 668286, 373493, 731302, 518190, 553336, 897626, 682656, 580591, 915612, 949596, 411371, 660301, 793095, 651386, 289081, 958856, 346709, 523205, 930560, 500896, 1056651, 581860, 653012, 1000323, 592979, 1042845, 418936, 494984, 447825, 653758, 250286, 233175, 532844, 567222, 508867, 941194, 607007, 275898, 565432, 654433, 1138135, 465458, 433568, 983439, 765617, 810132, 570295, 1135444, 1039781, 346336, 851548, 607131, 312623, 836800, 348436, 756059, 275085, 903306, 955898, 374842, 597348, 741079 };
            double[] yVals = new double[] { 1546.68852634046, 1534.53839102809, 1535.14271379466, 1520.96327558071, 1520.11326342835, 1525.23352314156, 1547.57339218597, 1536.04622929759, 1531.76017947961, 1520.02697475782, 1552.58654907816, 1544.18292492383, 1521.61996111246, 1551.6945661609, 1531.16820225304, 1520.60752591101, 1532.28428772283, 1521.28723655429, 1534.24368771918, 1525.89054199507, 1527.41552897469, 1539.92476129301, 1532.85798024683, 1528.5984616186, 1540.95517343148, 1542.85357065354, 1521.95492864831, 1531.96550311572, 1535.62225964376, 1531.60959018243, 1520.50743737088, 1543.36800170134, 1520.9152494705, 1526.10780528454, 1541.79604263778, 1525.16998368112, 1548.56008298323, 1528.65354643794, 1531.67450486886, 1545.67166393236, 1529.13620057993, 1547.96088503871, 1522.19132669005, 1524.92365890073, 1523.09407659456, 1531.70428737567, 1519.97607374744, 1519.73842922121, 1526.52615013368, 1528.01819902087, 1525.5020971494, 1542.38680503112, 1529.74512874829, 1520.33178361862, 1527.94051074615, 1531.73123535436, 1552.09659213222, 1523.69345161626, 1522.64856081386, 1544.73368835881, 1535.04982118676, 1536.03507343495, 1528.151571126, 1551.97979931252, 1547.82790354502, 1520.91006909482, 1537.544979341, 1529.75051134687, 1520.67092064374, 1537.00731051005, 1520.93923474879, 1534.85070310059, 1520.32049234402, 1540.25016666088, 1543.20367264741, 1521.30597201486, 1529.32585036376, 1534.49828825388 };
            SampleDatasetY(_dataset1, xVals, yVals);
            SampleDatasetY(_dataset4, xVals, yVals);
            yVals = new double[] { 1520.56321710372, 1526.16934828652, 1525.01457241068, 1541.03898822958, 1550.11710357279, 1535.03323842922, 1520.50791298838, 1523.37466060206, 1529.96502230477, 1550.65640776362, 1519.71815214749, 1520.71981719226, 1538.3301604111, 1519.86086941426, 1530.60847581191, 1545.97958542386, 1529.39533943171, 1539.70264921357, 1526.66052046804, 1534.70472900247, 1533.79068261519, 1521.23339120169, 1528.77176060127, 1533.01741531089, 1520.98359431964, 1520.80290183415, 1537.70258324364, 1529.74184443943, 1524.05548071248, 1530.12870632345, 1547.58100206594, 1520.77074989367, 1541.23709593417, 1534.57531682928, 1520.86899733514, 1535.06500815944, 1520.36238672268, 1532.9667372771, 1530.05814688167, 1520.62677100423, 1532.52269546646, 1520.45825839381, 1537.42678552828, 1535.18817054963, 1536.37357730634, 1530.02577459166, 1550.9745390785, 1552.45981736745, 1534.32430991979, 1533.42908058748, 1534.8989514253, 1520.83207468556, 1531.96248155158, 1548.7513523836, 1533.47569355231, 1529.99648331048, 1519.79654525885, 1535.80327419187, 1536.89334571716, 1520.68539447757, 1525.20035762647, 1523.38422277004, 1533.3490573244, 1519.81523211, 1520.4795354328, 1541.25846498385, 1522.09001770772, 1531.95752956088, 1544.96526970018, 1522.55087670567, 1541.13815666123, 1525.59859379883, 1548.8219228499, 1521.15450505191, 1520.78102045954, 1539.6253654387, 1532.34821766534, 1526.23618624353 };
            SampleDatasetY(_dataset2, xVals, yVals);
            SampleDatasetY(_dataset3, xVals, yVals);

        }

        /// <summary>
        /// Test the method that get X from given Y values
        /// </summary>
        [TestMethod]
        public void Test_GetX()
        {
            double[] xVals = new double[] { 1018627, 742619, 770076, 350167, 260164, 502421, 1034555, 810438, 655158, 253951, 1149424, 973525, 397450, 1128872, 640330, 303494, 668286, 373493, 731302, 518190, 553336, 897626, 682656, 580591, 915612, 949596, 411371, 660301, 793095, 651386, 289081, 958856, 346709, 523205, 930560, 500896, 1056651, 581860, 653012, 1000323, 592979, 1042845, 418936, 494984, 447825, 653758, 250286, 233175, 532844, 567222, 508867, 941194, 607007, 275898, 565432, 654433, 1138135, 465458, 433568, 983439, 765617, 810132, 570295, 1135444, 1039781, 346336, 851548, 607131, 312623, 836800, 348436, 756059, 275085, 903306, 955898, 374842, 597348, 741079 };
            double[] yVals = new double[] { 1546.68852634046, 1534.53839102809, 1535.14271379466, 1520.96327558071, 1520.11326342835, 1525.23352314156, 1547.57339218597, 1536.04622929759, 1531.76017947961, 1520.02697475782, 1552.58654907816, 1544.18292492383, 1521.61996111246, 1551.6945661609, 1531.16820225304, 1520.60752591101, 1532.28428772283, 1521.28723655429, 1534.24368771918, 1525.89054199507, 1527.41552897469, 1539.92476129301, 1532.85798024683, 1528.5984616186, 1540.95517343148, 1542.85357065354, 1521.95492864831, 1531.96550311572, 1535.62225964376, 1531.60959018243, 1520.50743737088, 1543.36800170134, 1520.9152494705, 1526.10780528454, 1541.79604263778, 1525.16998368112, 1548.56008298323, 1528.65354643794, 1531.67450486886, 1545.67166393236, 1529.13620057993, 1547.96088503871, 1522.19132669005, 1524.92365890073, 1523.09407659456, 1531.70428737567, 1519.97607374744, 1519.73842922121, 1526.52615013368, 1528.01819902087, 1525.5020971494, 1542.38680503112, 1529.74512874829, 1520.33178361862, 1527.94051074615, 1531.73123535436, 1552.09659213222, 1523.69345161626, 1522.64856081386, 1544.73368835881, 1535.04982118676, 1536.03507343495, 1528.151571126, 1551.97979931252, 1547.82790354502, 1520.91006909482, 1537.544979341, 1529.75051134687, 1520.67092064374, 1537.00731051005, 1520.93923474879, 1534.85070310059, 1520.32049234402, 1540.25016666088, 1543.20367264741, 1521.30597201486, 1529.32585036376, 1534.49828825388 };
            SampleDatasetX(_dataset1, yVals, xVals);
            SampleDatasetX(_dataset4, yVals, xVals);
            yVals = new double[] { 1520.56321710372, 1526.16934828652, 1525.01457241068, 1541.03898822958, 1550.11710357279, 1535.03323842922, 1520.50791298838, 1523.37466060206, 1529.96502230477, 1550.65640776362, 1519.71815214749, 1520.71981719226, 1538.3301604111, 1519.86086941426, 1530.60847581191, 1545.97958542386, 1529.39533943171, 1539.70264921357, 1526.66052046804, 1534.70472900247, 1533.79068261519, 1521.23339120169, 1528.77176060127, 1533.01741531089, 1520.98359431964, 1520.80290183415, 1537.70258324364, 1529.74184443943, 1524.05548071248, 1530.12870632345, 1547.58100206594, 1520.77074989367, 1541.23709593417, 1534.57531682928, 1520.86899733514, 1535.06500815944, 1520.36238672268, 1532.9667372771, 1530.05814688167, 1520.62677100423, 1532.52269546646, 1520.45825839381, 1537.42678552828, 1535.18817054963, 1536.37357730634, 1530.02577459166, 1550.9745390785, 1552.45981736745, 1534.32430991979, 1533.42908058748, 1534.8989514253, 1520.83207468556, 1531.96248155158, 1548.7513523836, 1533.47569355231, 1529.99648331048, 1519.79654525885, 1535.80327419187, 1536.89334571716, 1520.68539447757, 1525.20035762647, 1523.38422277004, 1533.3490573244, 1519.81523211, 1520.4795354328, 1541.25846498385, 1522.09001770772, 1531.95752956088, 1544.96526970018, 1522.55087670567, 1541.13815666123, 1525.59859379883, 1548.8219228499, 1521.15450505191, 1520.78102045954, 1539.6253654387, 1532.34821766534, 1526.23618624353 };
            SampleDatasetX(_dataset2, yVals, xVals);
            SampleDatasetX(_dataset3, yVals, xVals);
        }

        /// <summary>
        /// Test method the SaveToXElement(), which saves the OrderedPairedData object as an XElement, 
        /// and test the conversion back to an OrderPairedData object from an XElement object
        /// </summary>
        [TestMethod()]
        public void Test_ReadWriteXElement()
        {
            var el1 = _dataset1.SaveToXElement();
            var el2 = _dataset2.SaveToXElement();
            var el3 = _dataset3.SaveToXElement();
            var el4 = _dataset4.SaveToXElement();

            OrderedPairedData newDataset1 = new OrderedPairedData(el1);
            OrderedPairedData newDataset2 = new OrderedPairedData(el2);
            OrderedPairedData newDataset3 = new OrderedPairedData(el3);
            OrderedPairedData newDataset4 = new OrderedPairedData(el4);

            Assert.IsTrue(_dataset1 == newDataset1);
            Assert.IsTrue(_dataset2 == newDataset2);
            Assert.IsTrue(_dataset3 == newDataset3);
            Assert.IsTrue(_dataset4 == newDataset4);
        }

        /// <summary>
        /// Test the various OrderedPairedData object indexing and manipulation methods
        /// </summary>
        [TestMethod]
        public void Test_Indexing()
        {
            var dataset11 = _dataset1.Clone();
            Ordinate ordinate = new Ordinate(460815, 1523.5);

            Ordinate test1 = dataset11[4];
            Assert.AreEqual(ordinate, test1);

            int test2 = dataset11.IndexOf(ordinate);
            int test3 = dataset11.IndexOf(460815, 1523.5);
            Assert.AreEqual(4, test2);
            Assert.AreEqual(4, test3);

            dataset11.Remove(ordinate);
            bool test4 = dataset11.Contains(ordinate);
            Assert.AreEqual(false, test4);

            Ordinate newOrdinate = dataset11[4];
            dataset11.RemoveAt(4);
            bool test5 = dataset11.Contains(newOrdinate);
            Assert.AreEqual(false, test5);

            Ordinate newOrdinate2 = new Ordinate(1243177, 1563.8);
            dataset11.Add(newOrdinate2);
            int test6 = dataset11.IndexOf(newOrdinate2);
            Assert.AreEqual(dataset11.Count - 1, test6);

            dataset11.Insert(4, ordinate);
            int test7 = dataset11.IndexOf(ordinate);
            Assert.AreEqual(4, test7);
            
            int prevCount = dataset11.Count();
            Ordinate newOrdinate3 = dataset11[3];
            dataset11.RemoveRange(0, 3);
            Ordinate test8 = dataset11[0];
            int currCount = dataset11.Count();
            Assert.AreEqual(prevCount-3, currCount);
            Assert.AreEqual(newOrdinate3, test8);

            var array = new Ordinate[12];
            dataset11.CopyTo(array, 0);
            for( int i = 0; i < array.Length; i++ )
            {
                Assert.AreEqual(dataset11[i], array[i]);
            }

            var inDataset11 = dataset11.Invert();
            for(int i = 0; i < dataset11.Count(); i++ )
            {
                Assert.AreEqual(dataset11[i].X, inDataset11[i].Y);
                Assert.AreEqual(dataset11[i].Y, inDataset11[i].X);
            }

            dataset11.Clear();
            Assert.AreEqual(0, dataset11.Count());
        }

        /// <summary>
        /// Test the overloaded equality operators for the OrderedPairedData object
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

        /// <summary>
        /// Test the methods of finding area under a curve (TrapezodialAreaUnderY and TrapezodialAreaUnderX) against
        /// R's "AUC()" function from "DescTools" package
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Andri Signorell et mult. al. (2017). DescTools: Tools for descriptive statistics. R package version 0.99.23.
        /// </remarks>
        [TestMethod]
        public void Test_TrapezoidalArea()
        {
            var areaY1 = _dataset1.TrapezoidalAreaUnderY();
            Assert.AreEqual(1413175623, areaY1, 1);

            var areaX1 = _dataset1.TrapezoidalAreaUnderX();
            Assert.AreEqual(25442742, areaX1, 1);

            var areaY2 = _dataset2.TrapezoidalAreaUnderY();
            Assert.AreEqual(1410070832, areaY2, 1);

            var areaX2 = _dataset2.TrapezoidalAreaUnderX();
            Assert.AreEqual(17073185, areaX2, 1);

            var areaY3 = _dataset3.TrapezoidalAreaUnderY();
            Assert.AreEqual(1410070832, areaY3, 1);

            var areaX3 = _dataset3.TrapezoidalAreaUnderX();
            Assert.AreEqual(17073185, areaX3, 1);

            var areaY4 = _dataset4.TrapezoidalAreaUnderY();
            Assert.AreEqual(1413175623, areaY4, 1);

            var areaX4 = _dataset4.TrapezoidalAreaUnderX();
            Assert.AreEqual(25442742, areaX4, 1);
        }

    }
}
