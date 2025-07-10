/**
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
* **/

using Microsoft.VisualStudio.TestTools.UnitTesting;
using System.Collections.Generic;
using Numerics.Mathematics.Optimization;
using System.Diagnostics;
using System;
using Numerics.Distributions;
namespace Mathematics.Optimization
{
    [TestClass]
    public class ShortestPathTesting
    {

        [TestMethod]
        public void SimpleNetworkRouting()
        {

            //Simple Network Node Setup
            // 0 - 1 - 2 - 3 - 4
            // |   | \ |   |   |
            // 5 - 6 - 7 - 8 - 9


            List<Edge> edges = new List<Edge>();
            edges.Add(new Edge(0, 5, 1, 0));
            edges.Add(new Edge(0, 1, 30, 1));

            edges.Add(new Edge(1, 0, 30, 1));
            edges.Add(new Edge(1, 2, 1, 2));
            edges.Add(new Edge(1, 6, 15, 3));
            edges.Add(new Edge(1, 7, 2, 4));

            edges.Add(new Edge(2, 1, 1, 2));
            edges.Add(new Edge(2, 3, 5, 5));
            edges.Add(new Edge(2, 7, 5, 6));

            edges.Add(new Edge(3, 2, 5, 5));
            edges.Add(new Edge(3, 8, 2, 7));
            edges.Add(new Edge(3, 4, 1, 8));

            edges.Add(new Edge(4, 3, 1, 8));
            edges.Add(new Edge(4, 9, 30, 9));

            edges.Add(new Edge(5, 0, 1, 0));
            edges.Add(new Edge(5, 6, 3, 10));

            edges.Add(new Edge(6, 5, 3, 10));
            edges.Add(new Edge(6, 1, 15, 3));
            edges.Add(new Edge(6, 7, 1, 11));

            edges.Add(new Edge(7, 6, 1, 11));
            edges.Add(new Edge(7, 1, 2, 4));
            edges.Add(new Edge(7, 2, 5, 6));
            edges.Add(new Edge(7, 8, 1, 12));

            edges.Add(new Edge(8, 7, 1, 12));
            edges.Add(new Edge(8, 3, 2, 7));
            edges.Add(new Edge(8, 9, 2, 13));

            edges.Add(new Edge(9, 8, 2, 13));
            edges.Add(new Edge(9, 4, 30, 9));


            float[,] result = Dijkstra.Solve(edges.ToArray(), 9);

            Assert.AreEqual(result[0, 0] == 5 && result[0, 2] == 8, true);
            Assert.AreEqual(result[1, 0] == 7 && result[1, 2] == 5, true);
            Assert.AreEqual(result[2, 0] == 1 && result[2, 2] == 6, true);
            Assert.AreEqual(result[3, 0] == 8 && result[3, 2] == 4, true);
            Assert.AreEqual(result[4, 0] == 3 && result[4, 2] == 5, true);
            Assert.AreEqual(result[5, 0] == 6 && result[5, 2] == 7, true);
            Assert.AreEqual(result[6, 0] == 7 && result[6, 2] == 4, true);
            Assert.AreEqual(result[7, 0] == 8 && result[7, 2] == 3, true);
            Assert.AreEqual(result[8, 0] == 9 && result[8, 2] == 2, true);
            Assert.AreEqual(result[9, 0] == 9 && result[9, 2] == 0, true);

            //for (int i = 0; i < result.GetLength(0); i++)
            //{
            //    Debug.Print($"{i}, {result[i, 0]}, {result[i, 2]}");
            //}

        }

        [TestMethod]
        public void BidirectionalRouting()
        {
            List<Edge> edges = new List<Edge>();
            edges.Add(new Edge(0, 1, 6, 0));
            edges.Add(new Edge(0, 3, 1, 1));

            edges.Add(new Edge(1, 0, 6, 0));
            edges.Add(new Edge(1, 2, 5, 2));
            edges.Add(new Edge(1, 3, 2, 3));
            edges.Add(new Edge(1, 4, 2, 4));

            edges.Add(new Edge(2, 1, 5, 2));
            edges.Add(new Edge(2, 4, 5, 5));

            edges.Add(new Edge(3, 0, 1, 1));
            edges.Add(new Edge(3, 1, 2, 3));
            edges.Add(new Edge(3, 4, 1, 6));

            edges.Add(new Edge(4, 1, 2, 4));
            edges.Add(new Edge(4, 2, 5, 5));
            edges.Add(new Edge(4, 3, 1, 6));

            float[,] result = Dijkstra.Solve(edges.ToArray(), 4);
            Assert.AreEqual(result[0, 0] == 3 && result[0, 2] == 2, true);
            Assert.AreEqual(result[1, 0] == 4 && result[1, 2] == 2, true);
            Assert.AreEqual(result[2, 0] == 4 && result[2, 2] == 5, true);
            Assert.AreEqual(result[3, 0] == 4 && result[3, 2] == 1, true);
            Assert.AreEqual(result[4, 0] == 4 && result[4, 2] == 0, true);

            //    for (int i = 0; i < result.GetLength(0); i++)
            //    {
            //        Debug.Print($"{i}, {result[i, 0]}, {result[i, 2]}");
            //    }

        }
    }
}
