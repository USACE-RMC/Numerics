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
    /// Unit tests for the Ramer-Douglas-Peucker line simplification algorithm. The algorithm was tested against R's
    /// "RamerDouglasPeucker( )" method from the "RDP" package. The other edge case tests are modeled after the validation testing of the
    /// same package.
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
    /// Robert Dahl Jacobsen (2023). RDP: The Ramer-Douglas-Peucker Algorithm. <see href="https://github.com/robertdj/RDP"/>
    /// </remarks>
    [TestClass]
    public class Test_LineSimplification
    {
        /// <summary>
        /// Test the RDP algorithm
        /// </summary>
        [TestMethod]
        public void Test_RamerDouglasPeucker()
        {
            // coordinates of the sin(x) function
            var inputs = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(3.14 / 2, 1), new Ordinate(3.14, 0), new Ordinate(3 * 3.14 / 2, -1), new Ordinate(2 * 3.14, 0) };
            var outputs = new List<Ordinate>();
            double epsilon = 0.1;

            LineSimplification.RamerDouglasPeucker(inputs, epsilon, ref outputs);
            var valid = new List<Ordinate>() { new Ordinate(0, 0), new Ordinate(1.57, 1), new Ordinate(4.71, -1), new Ordinate(6.28, 0) };

            for(int i = 0; i < valid.Count; i++)
            {
                Assert.AreEqual(valid[i].X, outputs[i].X, 1E-6);
                Assert.AreEqual(valid[i].Y, outputs[i].Y, 1E-6);
            }
        }

        /// <summary>
        /// Test when epsilon is zero, all points should be kept.
        /// </summary>
        [TestMethod]
        public void Test_ZeroEpsilon()
        {
            var inputs = new List<Ordinate>() { new Ordinate(3.6, 8.5), new Ordinate(19.66, 0.33), new Ordinate(88.17, 64.9), new Ordinate(-5.63, 93.2), new Ordinate(-22.35, -7.5), new Ordinate(-2, -2)};
            var outputs = new List<Ordinate>();
            double epsilon = 0;

            LineSimplification.RamerDouglasPeucker(inputs, epsilon, ref outputs);
            CollectionAssert.AreEquivalent(inputs, outputs);

        }

        /// <summary>
        /// Test the algorithm on a linear line, only the endpoints should be kept.
        /// </summary>
        [TestMethod]
        public void Test_Line()
        {
            var inputs = new List<Ordinate>() { new Ordinate(1, 2), new Ordinate(2, 4), new Ordinate(3, 6), new Ordinate(4, 8), new Ordinate(5, 10), new Ordinate(6, 12) };
            var outputs = new List<Ordinate>();
            double epsilon = 0.1;

            LineSimplification.RamerDouglasPeucker(inputs, epsilon, ref outputs);
            var expected = new List<Ordinate>() { inputs[0], inputs[inputs.Count - 1] };

            CollectionAssert.AreEquivalent(expected, outputs);
        }

        /// <summary>
        /// Test the algorithm when all input points are equal, only the endpoints should be kept.
        /// </summary>
        [TestMethod]
        public void Test_EqualPoints()
        {
            var inputs = new List<Ordinate> { new Ordinate(1,30), new Ordinate(1, 30), new Ordinate(1, 30), new Ordinate(1, 30), new Ordinate(1, 30), new Ordinate(1, 30) };
            var outputs = new List<Ordinate>();
            double epsilon = 0.1;

            LineSimplification.RamerDouglasPeucker(inputs, epsilon, ref outputs);
            var expected = new List<Ordinate>() { inputs[0], inputs[inputs.Count - 1] };

            CollectionAssert.AreEquivalent(expected, outputs);
        }
    }
}
