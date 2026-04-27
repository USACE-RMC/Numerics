using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using System.Collections.Generic;

namespace Data.PairedData
{
    /// <summary>
    /// Test the OrderedPairData's line simplification methods. The results were test against values obtained from
    /// R's "RamerDouglasPeucker( )" method from the "RDP" package.
    /// </summary>
    /// <remarks>
    /// <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    ///     Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <b> References: </b>
    /// Robert Dahl Jacobsen (2023). RDP: The Ramer-Douglas-Peucker Algorithm. https://github.com/robertdj/RDP
    /// </remarks>
    [TestClass]
    public class Test_PairedDataLineSimplification
    {
        /// <summary>
        /// Test the Douglas-Peucker simplification algorithm
        /// </summary>
        [TestMethod]
        public void Test_DouglasPeuckerSimplify()
        {
            var data = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(3.14 / 2, 1), new Ordinate(3.14, 0), new Ordinate(3 * 3.14 / 2, -1), new Ordinate(2 * 3.14, 0) };
            var orderedPair = new OrderedPairedData(data, true, SortOrder.Ascending, false, SortOrder.None);
            var lineSimp = new List<Ordinate>();
            LineSimplification.RamerDouglasPeucker(data, 0.01, ref lineSimp);

            var test = orderedPair.DouglasPeuckerSimplify(0.01);
            var valid = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(1.57, 1), new Ordinate(4.71, -1), new Ordinate(6.28, 0) };

            for (int i = 0; i < test.Count; i++)
            {
                Assert.AreEqual(lineSimp[i].X, test[i].X);
                Assert.AreEqual(lineSimp[i].Y, test[i].Y);
                Assert.AreEqual(valid[i].X, test[i].X);
                Assert.AreEqual(valid[i].Y, test[i].Y);
            }
        }

        /// <summary>
        /// Test the Visvaligam-Whyatt simplification algorithm
        /// </summary>
        [TestMethod]
        public void Test_VisvaligamWhyattSimplify()
        {
            var data = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(3.14 / 2, 1), new Ordinate(3.14, 0), new Ordinate(3 * 3.14 / 2, -1), new Ordinate(2 * 3.14, 0) };
            var orderedPair = new OrderedPairedData(data, true, SortOrder.Ascending, false, SortOrder.None);

            var test = orderedPair.VisvaligamWhyattSimplify(4);
            var valid = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(1.57, 1), new Ordinate(4.71, -1), new Ordinate(6.28, 0) };
            for (int i = 0; i < test.Count; i++)
            {
                Assert.AreEqual(valid[i].X, test[i].X);
                Assert.AreEqual(valid[i].Y, test[i].Y);
            }
        }

        /// <summary>
        /// Test the Lang simplification algorithm
        /// </summary>
        [TestMethod]
        public void Test_LangSimplify()
        {
            var data = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(3.14 / 2, 1), new Ordinate(3.14, 0), new Ordinate(3 * 3.14 / 2, -1), new Ordinate(2 * 3.14, 0) };
            var orderedPair = new OrderedPairedData(data, true, SortOrder.Ascending, false, SortOrder.None);

            var test = orderedPair.LangSimplify(0.01, 2);
            var valid = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(1.57, 1), new Ordinate(4.71, -1), new Ordinate(6.28, 0) };
            for (int i = 0; i < test.Count; i++)
            {
                Assert.AreEqual(valid[i].X, test[i].X);
                Assert.AreEqual(valid[i].Y, test[i].Y);
            }
        }
    }
}
