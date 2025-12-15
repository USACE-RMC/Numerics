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


using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    [TestClass]
    public class ShortestPathTesting
    {
        /// <summary>
        /// Testing a something that cost doesn't really matter.
        /// </summary>
        [TestMethod]
        public void SimpleEdgeGraphCost()
        {
            List<Edge> edges = new List<Edge>();
            edges.Add(new Edge(0, 1, 2, 0));
            edges.Add(new Edge(0, 2, 4, 2));
            edges.Add(new Edge(1, 2, 1, 2));
            edges.Add(new Edge(1, 3, 7, 3));
            edges.Add(new Edge(2, 3, 3, 4));
            edges.Add(new Edge(4, 0, 1, 5));

            float[,] result = Dijkstra.Solve(edges, 3,6);

            Assert.AreEqual(result[3, 2], 0f);
            Assert.AreEqual(result[2, 2], 3f);
            Assert.AreEqual(result[1, 2], 4f);
            Assert.AreEqual(result[0, 2], 6f);
            Assert.AreEqual(result[4, 2], 7f);
            Assert.IsTrue(float.IsPositiveInfinity(result[5,2]));

        }

        /// <summary>
        /// Simple network run, testing to see if algorithm chooses 
        /// the lowest cost path as it should.
        /// </summary>
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


            float[,] result = Dijkstra.Solve(edges,9);

            Assert.AreEqual(result[0, 0], 5f); //Algorithm is choosing the next node that yields the shortest paths
            Assert.AreEqual(result[0, 2], 8f);

            Assert.AreEqual(result[1, 0], 7);
            Assert.AreEqual(result[1, 2], 5);

            Assert.AreEqual(1, result[2, 0]);
            Assert.AreEqual(6, result[2, 2]);

            Assert.AreEqual(8, result[3, 0]);
            Assert.AreEqual(4, result[3, 2]);

            Assert.AreEqual(3, result[4, 0]);
            Assert.AreEqual(5, result[4, 2]);

            Assert.AreEqual(6, result[5, 0]);
            Assert.AreEqual(7, result[5, 2]);

            Assert.AreEqual(7, result[6, 0]);
            Assert.AreEqual(4, result[6, 2]);

            Assert.AreEqual(8, result[7, 0]);
            Assert.AreEqual(3, result[7, 2]);

            Assert.AreEqual(9, result[8, 0]);
            Assert.AreEqual(2, result[8, 2]);

            Assert.AreEqual(9, result[9, 0]);
            Assert.AreEqual(0, result[9, 2]);

        }

        /// <summary>
        /// Testing edges with with bidirectionality.
        /// </summary>
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

            Assert.AreEqual(3, result[0, 0]);
            Assert.AreEqual(2, result[0, 2]);

            Assert.AreEqual(4, result[1, 0]);
            Assert.AreEqual(2, result[1, 2]);

            Assert.AreEqual(4, result[2, 0]);
            Assert.AreEqual(5, result[2, 2]);

            Assert.AreEqual(4, result[3, 0]);
            Assert.AreEqual(1, result[3, 2]);

            Assert.AreEqual(4, result[4, 0]);
            Assert.AreEqual(0, result[4, 2]);

        }

        /// <summary>
        /// Testing that a disconnected node returns a positive infinity.
        /// </summary>
        [TestMethod]
        public void DisconnectedNodesTest()
        {
            var edges = new List<Edge>
        {
            new Edge(0, 1, 1, 0),
             new Edge(1, 2, 1, 1),
             // Node 3 is disconnected
        };

            var result = Dijkstra.Solve(edges, 2,4);

            Assert.AreEqual(2, result[1, 0]);
            Assert.AreEqual(1, result[1, 2]);

            Assert.AreEqual(1, result[0, 0]);
            Assert.AreEqual(2, result[0, 2]);

            // Unreachable node 3 should remain with default values
            Assert.IsTrue(float.IsPositiveInfinity(result[3, 2]));
        }

        /// <summary>
        /// Simple multiple dest path
        /// </summary>
        [TestMethod]
        public void MultipleDestSharedPath()
        {
            // Graph:
            // 0 - 1 - 2
            //     |
            //     3

            var edges = new List<Edge>
            {
                new Edge(0,1,1,0),
                new Edge(1,0,3,1),
                new Edge(1,2,1,2),
                new Edge(2,1,2,3),
                new Edge(1,3,3,4)
            };

            var result = Dijkstra.Solve(edges, [0,3],4);

            Assert.AreEqual(result[1, 0], 0);
            Assert.AreEqual(result[1, 2], 3);
            Assert.AreEqual(result[2, 0], 1);
            Assert.AreEqual(result[2, 2], 5);
        }

        /// <summary>
        /// Checking paths are indeed disconnected.
        /// </summary>
        [TestMethod]
        public void DisconnectedComponent()
        {
            // Graph: 
            // 0 - 1      2 - 3
            var edges = new List<Edge>
            {
                new Edge(0,1,1,0),
                new Edge(1,0,3,1),
                new Edge(2,3,1,2)
            };
            var result = Dijkstra.Solve(edges, [0, 3], 4);
            Assert.AreEqual(result[1,0],0);
            Assert.AreEqual(result[1, 2], 3);
            Assert.AreEqual(result[2, 0], 3);
            Assert.AreEqual(result[2, 2], 1);
        }

        /// <summary>
        /// Checking disconnected with 1 destination.
        /// </summary>
        [TestMethod]
        public void DisconnectedComponent2()
        {
            var edges = new List<Edge>
            {
                new Edge(0,1,1,0),
                new Edge(1,0,3,1),
                new Edge(2,3,1,2),
            };

            var result = Dijkstra.Solve(edges, [0], 4);
            Assert.IsTrue(float.IsPositiveInfinity(result[2, 2]));
        }

        /// <summary>
        /// Testing two destinations when they are connected by an edge.
        /// </summary>
        [TestMethod]
        public void TrianglePath()
        {
            // Graph:
            // 0 <-> 1
            // 1 -> 2
            // 2 -> 0
            var edges = new List<Edge>
            {
                new Edge(0,1,1,0),
                new Edge(1,0,4,1),
                new Edge(1,2,1,2),
                new Edge(2,0,10,3)
            };

            var result = Dijkstra.Solve(edges, [0, 2], 3);
            Assert.AreEqual(result[0, 0], 0);
            Assert.AreEqual(result[0, 2], 0);
            Assert.AreEqual(result[1, 0], 2);
            Assert.AreEqual(result[1, 2], 1);
            Assert.AreEqual(result[2, 0], 2);
        }
    }
}
