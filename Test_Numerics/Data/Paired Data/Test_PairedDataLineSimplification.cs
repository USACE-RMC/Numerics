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
