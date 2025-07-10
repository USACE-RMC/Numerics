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
using Numerics.Data.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the Correlation class. These methods were tested against values attained from R's "cor()" function from the "stats" package.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, 
    /// Vienna, Austria. ISBN 3-900051-07-0, URL <see href=" http://www.R-project.org/."/>
    /// </remarks>
    [TestClass]
    public class Test_Correlation
    {
        /// <summary>
        /// Test the Pearson correlation coefficient on a small dataset
        /// </summary>
        [TestMethod()]
        public void Test_Pearson()
        {
            var XArray = new double[] { 14d, 8d, 32d, 7d, 3d, 15d };
            var YArray = new double[] { 10d, 5d, 7d, 4d, 3d, 8d };
            double R = Correlation.Pearson(XArray, YArray);
            double trueR = 0.54502739907793d;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Test the Pearson correlation coefficient on a larger dataset
        /// </summary>
        [TestMethod]
        public void Test_Pearson_Big()
        {
            var XArray = new double[] { 230408, 288010, 345611, 403213, 460815, 518417, 576019, 633612, 691223, 748825, 806427, 864029, 921631, 1036834, 1152038 };
            var YArray = new double[] { 1519.7, 1520.5, 1520.9, 1521.7, 1523.5, 1525.9, 1528.4, 1530.9, 1533.2, 1534.7, 1535.9, 1538, 1541.3, 1547.7, 1552.7 };
            double R = Correlation.Pearson(XArray, YArray);
            double trueR = 0.988054377242161;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Test the Spearman ranked correlation coefficient on a small dataset
        /// </summary>
        [TestMethod()]
        public void Test_Spearman()
        {
            var XArray = new double[] { 14d, 8d, 32d, 7d, 3d, 15d };
            var YArray = new double[] { 10d, 5d, 7d, 4d, 3d, 8d };
            double R = Correlation.Spearman(XArray, YArray);
            double trueR = 0.771428571428571d;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Test the Spearman ranked correlation coefficient on a larger dataset
        /// </summary>
        [TestMethod]
        public void Test_Spearman_Big()
        {
            var XArray = new double[] { 230408, 288010, 345611, 403213, 460815, 518417, 576019, 633612, 691223, 748825, 806427, 864029, 921631, 1036834, 1152038 };
            var YArray = new double[] { 1519.7, 1520.5, 1520.9, 1521.7, 1523.5, 1525.9, 1528.4, 1530.9, 1533.2, 1534.7, 1535.9, 1538, 1541.3, 1547.7, 1552.7 };
            double R = Correlation.Spearman(XArray, YArray);
            double trueR = 1;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Test the Kendall's Tau ranked correlation coefficient on a small dataset
        /// </summary>
        [TestMethod()]
        public void Test_KendallsTau()
        {
            var XArray = new double[] { 14d, 8d, 32d, 7d, 3d, 15d };
            var YArray = new double[] { 10d, 5d, 7d, 4d, 3d, 8d };
            double R = Correlation.KendallsTau(XArray, YArray);
            double trueR = 0.6d;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }

        /// <summary>
        /// Test the Kendall's Tau ranked correlation coefficient on a larger dataset
        /// </summary>
        [TestMethod]
        public void Test_KendallsTau_Big()
        {
            var XArray = new double[] { 230408, 288010, 345611, 403213, 460815, 518417, 576019, 633612, 691223, 748825, 806427, 864029, 921631, 1036834, 1152038 };
            var YArray = new double[] { 1519.7, 1520.5, 1520.9, 1521.7, 1523.5, 1525.9, 1528.4, 1530.9, 1533.2, 1534.7, 1535.9, 1538, 1541.3, 1547.7, 1552.7 };
            double R = Correlation.KendallsTau(XArray, YArray);
            double trueR = 1;
            Assert.AreEqual(R, trueR, 0.0000000001d);
        }
    }
}
