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
        /// Test the Lang simplification algorithm. The count is asserted before the contents: an
        /// earlier defect dropped the curve's final point, and a loop bounded by the short result's
        /// own length compared only the surviving points and passed.
        /// </summary>
        [TestMethod]
        public void Test_LangSimplify()
        {
            var data = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(3.14 / 2, 1), new Ordinate(3.14, 0), new Ordinate(3 * 3.14 / 2, -1), new Ordinate(2 * 3.14, 0) };
            var orderedPair = new OrderedPairedData(data, true, SortOrder.Ascending, false, SortOrder.None);

            var test = orderedPair.LangSimplify(0.01, 2);
            var valid = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(1.57, 1), new Ordinate(4.71, -1), new Ordinate(6.28, 0) };
            Assert.AreEqual(valid.Count, test.Count);
            for (int i = 0; i < test.Count; i++)
            {
                Assert.AreEqual(valid[i].X, test[i].X);
                Assert.AreEqual(valid[i].Y, test[i].Y);
            }
        }

        /// <summary>
        /// The guarded early return of LangSimplify hands back a distinct object, matching the
        /// other simplifiers' contract; it previously returned the receiver itself, so mutating the
        /// "simplified" result silently mutated the original curve.
        /// </summary>
        [TestMethod]
        public void Test_LangSimplify_GuardReturnsADistinctObject()
        {
            var data = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(1, 1), new Ordinate(2, 0) };
            var orderedPair = new OrderedPairedData(data, true, SortOrder.Ascending, false, SortOrder.None);

            var guarded = orderedPair.LangSimplify(0d, 2);
            Assert.AreNotSame(orderedPair, guarded);
            Assert.AreEqual(orderedPair.Count, guarded.Count);

            guarded.Add(new Ordinate(3, 5));
            Assert.AreEqual(3, orderedPair.Count);
        }

        /// <summary>
        /// RamerDouglasPeucker's output parameter holds the simplified curve alone on every path.
        /// The recursion branch previously appended to whatever the caller's list already held,
        /// while the endpoint branch replaced it, so a pre-populated list's final content depended
        /// on which branch the curve happened to take.
        /// </summary>
        [TestMethod]
        public void Test_RamerDouglasPeucker_ClearsThePrePopulatedOutput()
        {
            var data = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(3.14 / 2, 1), new Ordinate(3.14, 0), new Ordinate(3 * 3.14 / 2, -1), new Ordinate(2 * 3.14, 0) };
            var prePopulated = new List<Ordinate>() { new Ordinate(-99, -99) };
            LineSimplification.RamerDouglasPeucker(data, 0.01, ref prePopulated);
            Assert.HasCount(4, prePopulated);
            Assert.AreEqual(0d, prePopulated[0].X);

            var alsoPrePopulated = new List<Ordinate>() { new Ordinate(-99, -99) };
            LineSimplification.RamerDouglasPeucker(new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(1, 0) }, 0.01, ref alsoPrePopulated);
            Assert.HasCount(2, alsoPrePopulated);
            Assert.AreEqual(0d, alsoPrePopulated[0].X);
        }

        /// <summary>
        /// A curve whose first and last points coincide has a zero-length base segment, where the
        /// perpendicular distance degenerates to the point-to-endpoint distance. Unguarded, the
        /// 0/0 = NaN distance failed every comparison and the interior spike was silently dropped.
        /// </summary>
        [TestMethod]
        public void Test_DouglasPeuckerSimplify_CoincidentEndpoints()
        {
            var data = new List<Ordinate>() { new Ordinate(1, 5), new Ordinate(1, 7), new Ordinate(1, 5) };
            var orderedPair = new OrderedPairedData(data, false, SortOrder.Ascending, false, SortOrder.None);
            var test = orderedPair.DouglasPeuckerSimplify(0.5);
            Assert.AreEqual(3, test.Count);
            Assert.AreEqual(7d, test[1].Y);
        }
    }
}
