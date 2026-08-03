
using System.Diagnostics;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Tests shortest-path routing behavior for Dijkstra networks.
    /// </summary>
    [TestClass]
    public class Test_ShortestPath
    {
        /// <summary>
        /// Tests predecessor routing and cumulative costs on a simple edge graph.
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

            Assert.AreEqual(0f, result[3, 2]);
            Assert.AreEqual(3f, result[2, 2]);
            Assert.AreEqual(4f, result[1, 2]);
            Assert.AreEqual(6f, result[0, 2]);
            Assert.AreEqual(7f, result[4, 2]);
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

            Assert.AreEqual(5f, result[0, 0]); //Algorithm is choosing the next node that yields the shortest paths
            Assert.AreEqual(8f, result[0, 2]);

            Assert.AreEqual(7, result[1, 0]);
            Assert.AreEqual(5, result[1, 2]);

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

            Assert.AreEqual(0, result[1, 0]);
            Assert.AreEqual(3, result[1, 2]);
            Assert.AreEqual(1, result[2, 0]);
            Assert.AreEqual(5, result[2, 2]);
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
            Assert.AreEqual(0, result[1, 0]);
            Assert.AreEqual(3, result[1, 2]);
            Assert.AreEqual(3, result[2, 0]);
            Assert.AreEqual(1, result[2, 2]);
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
            Assert.AreEqual(0, result[0, 0]);
            Assert.AreEqual(0, result[0, 2]);
            Assert.AreEqual(2, result[1, 0]);
            Assert.AreEqual(1, result[1, 2]);
            Assert.AreEqual(2, result[2, 0]);
        }

        /// <summary>
        /// Pins the edge-index column: every routed node records the index of the edge it takes,
        /// the destination and unreachable nodes record -1.
        /// </summary>
        [TestMethod]
        public void ResultTableEdgeIndexColumn()
        {
            var edges = new List<Edge>
            {
                new Edge(0, 1, 2, 0),
                new Edge(0, 2, 4, 2),
                new Edge(1, 2, 1, 2),
                new Edge(1, 3, 7, 3),
                new Edge(2, 3, 3, 4),
                new Edge(4, 0, 1, 5),
            };

            var result = Dijkstra.Solve(edges, 3, 6);

            Assert.AreEqual(0f, result[0, 1]);   // 0 departs on edge 0 toward node 1
            Assert.AreEqual(2f, result[1, 1]);   // 1 departs on edge 2 toward node 2
            Assert.AreEqual(4f, result[2, 1]);   // 2 departs on edge 4 toward node 3
            Assert.AreEqual(-1f, result[3, 1]);  // the destination takes no edge
            Assert.AreEqual(5f, result[4, 1]);   // 4 departs on edge 5 toward node 0
            Assert.AreEqual(-1f, result[5, 1]);  // unreachable
            Assert.AreEqual(1f, result[0, 0]);
            Assert.AreEqual(2f, result[1, 0]);
            Assert.AreEqual(3f, result[2, 0]);
            Assert.AreEqual(0f, result[4, 0]);
        }

        /// <summary>
        /// The destination row is exactly (itself, -1, 0) for single- and multi-destination solves.
        /// </summary>
        [TestMethod]
        public void DestinationRowContract()
        {
            var edges = new List<Edge> { new Edge(0, 1, 1, 0), new Edge(1, 2, 1, 1) };

            var single = Dijkstra.Solve(edges, 2, 3);
            Assert.AreEqual(2f, single[2, 0]);
            Assert.AreEqual(-1f, single[2, 1]);
            Assert.AreEqual(0f, single[2, 2]);

            var multi = Dijkstra.Solve(edges, [0, 2], 3);
            Assert.AreEqual(0f, multi[0, 0]);
            Assert.AreEqual(-1f, multi[0, 1]);
            Assert.AreEqual(0f, multi[0, 2]);
            Assert.AreEqual(2f, multi[2, 0]);
            Assert.AreEqual(-1f, multi[2, 1]);
            Assert.AreEqual(0f, multi[2, 2]);
        }

        /// <summary>
        /// An unreachable node's row is exactly (-1, -1, positive infinity) in all three columns.
        /// </summary>
        [TestMethod]
        public void UnreachableRowContract()
        {
            var edges = new List<Edge> { new Edge(0, 1, 1, 0) };

            var result = Dijkstra.Solve(edges, 1, 3);

            Assert.AreEqual(-1f, result[2, 0]);
            Assert.AreEqual(-1f, result[2, 1]);
            Assert.IsTrue(float.IsPositiveInfinity(result[2, 2]));
        }

        /// <summary>
        /// Solving the same inputs twice produces element-for-element identical tables.
        /// </summary>
        [TestMethod]
        public void RepeatedSolveIsBitIdentical()
        {
            var edges = BuildRandomGraph(new MersenneTwister(45678), 300, 1200);

            var first = Dijkstra.Solve(edges, 7, 300);
            var second = Dijkstra.Solve(edges, 7, 300);

            for (int i = 0; i < 300; i++)
            {
                for (int c = 0; c < 3; c++)
                {
                    Assert.AreEqual(first[i, c], second[i, c], 0f, $"Cell [{i},{c}] differs between runs.");
                }
            }
        }

        /// <summary>
        /// On an exact cost tie between destinations, the earlier destination in the array wins;
        /// reversing the array flips the winner.
        /// </summary>
        [TestMethod]
        public void MultiDestinationTieBreaksByArrayOrder()
        {
            var edges = new List<Edge>
            {
                new Edge(0, 1, 5, 0),
                new Edge(0, 2, 5, 1),
            };

            var forward = Dijkstra.Solve(edges, [1, 2], 3);
            Assert.AreEqual(1f, forward[0, 0]);
            Assert.AreEqual(0f, forward[0, 1]);
            Assert.AreEqual(5f, forward[0, 2]);

            var reversed = Dijkstra.Solve(edges, [2, 1], 3);
            Assert.AreEqual(2f, reversed[0, 0]);
            Assert.AreEqual(1f, reversed[0, 1]);
            Assert.AreEqual(5f, reversed[0, 2]);
        }

        /// <summary>
        /// The cost column equals an independent Bellman-Ford re-derivation on fixed-seed random
        /// graphs of several sizes and densities (integer weights make the comparison exact).
        /// </summary>
        [TestMethod]
        public void SingleDestinationMatchesBellmanFordOnRandomGraphs()
        {
            (int nodes, int edges, int seed)[] cases = { (200, 800, 101), (500, 2500, 102), (1000, 6000, 103) };
            foreach (var (nodeCount, edgeCount, seed) in cases)
            {
                var edges = BuildRandomGraph(new MersenneTwister(seed), nodeCount, edgeCount);
                int destination = nodeCount / 2;

                var result = Dijkstra.Solve(edges, destination, nodeCount);
                float[] oracle = BellmanFordToDestination(edges, nodeCount, destination);

                for (int i = 0; i < nodeCount; i++)
                {
                    Assert.AreEqual(oracle[i], result[i, 2], 0f, $"Cost mismatch at node {i} (n={nodeCount}).");
                }
            }
        }

        /// <summary>
        /// Multi-destination merged costs equal the per-destination oracle minima.
        /// </summary>
        [TestMethod]
        public void MultiDestinationMatchesPerDestinationOracleMin()
        {
            var edges = BuildRandomGraph(new MersenneTwister(104), 400, 1600);
            int[] destinations = { 3, 200, 397 };

            var result = Dijkstra.Solve(edges, destinations, 400);

            var oracles = new float[destinations.Length][];
            for (int d = 0; d < destinations.Length; d++)
            {
                oracles[d] = BellmanFordToDestination(edges, 400, destinations[d]);
            }
            for (int i = 0; i < 400; i++)
            {
                float expected = float.PositiveInfinity;
                for (int d = 0; d < destinations.Length; d++)
                {
                    if (oracles[d][i] < expected) expected = oracles[d][i];
                }
                Assert.AreEqual(expected, result[i, 2], 0f, $"Merged cost mismatch at node {i}.");
            }
        }

        /// <summary>
        /// Structural self-consistency on random graphs: every reachable non-destination row
        /// names a real edge from the node to its recorded next node, and the costs telescope
        /// exactly along that edge.
        /// </summary>
        [TestMethod]
        public void ResultTableIsSelfConsistentOnRandomGraphs()
        {
            var edges = BuildRandomGraph(new MersenneTwister(105), 300, 1500);
            int destination = 17;

            var result = Dijkstra.Solve(edges, destination, 300);

            for (int i = 0; i < 300; i++)
            {
                if (i == destination || float.IsPositiveInfinity(result[i, 2])) continue;
                int next = (int)result[i, 0];
                int edgeIndex = (int)result[i, 1];

                bool found = false;
                for (int j = 0; j < edges.Count; j++)
                {
                    var edge = edges[j];
                    if (edge.FromIndex == i && edge.ToIndex == next && edge.Index == edgeIndex
                        && result[i, 2] == result[next, 2] + edge.Weight)
                    {
                        found = true;
                        break;
                    }
                }
                Assert.IsTrue(found, $"Node {i} routes over edge index {edgeIndex} to {next}, but no such edge telescopes the costs.");
            }
        }

        /// <summary>
        /// A complete directed graph keeps every node in flight at once; the solver's
        /// exact-node-count heap capacity must never overflow.
        /// </summary>
        [TestMethod]
        public void DenseGraphNeverOverflowsExactHeapCapacity()
        {
            const int n = 200;
            var edges = new List<Edge>(n * (n - 1));
            int index = 0;
            for (int a = 0; a < n; a++)
            {
                for (int b = 0; b < n; b++)
                {
                    if (a == b) continue;
                    edges.Add(new Edge(a, b, 1f + ((a + b) % 7), index++));
                }
            }

            var result = Dijkstra.Solve(edges, n - 1, n);

            for (int i = 0; i < n; i++)
            {
                Assert.IsTrue(Dijkstra.PathExists(result, i), $"Node {i} must reach the destination in a complete graph.");
            }
        }

        /// <summary>
        /// Malformed inputs fail with clear argument exceptions instead of raw index or LINQ
        /// crashes.
        /// </summary>
        [TestMethod]
        public void SolveValidationThrows()
        {
            var edges = new List<Edge> { new Edge(0, 1, 1, 0) };

            Assert.Throws<ArgumentNullException>(() => Dijkstra.Solve(null!, 0, 2));
            Assert.Throws<ArgumentNullException>(() => Dijkstra.Solve(edges, null!, 2));
            Assert.Throws<ArgumentException>(() => Dijkstra.Solve(new List<Edge>(), 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, 0, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, -1, 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, 2, 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, [0, 5], 2));
            Assert.Throws<ArgumentException>(() => Dijkstra.Solve(new List<Edge> { new Edge(0, 9, 1, 0) }, 0, 3));

            Assert.Throws<ArgumentNullException>(() => Dijkstra.PathExists(null!, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Dijkstra.PathExists(new float[2, 3], 2));
            Assert.Throws<ArgumentException>(() => Dijkstra.PathExists(new float[2, 2], 0));
        }

        /// <summary>
        /// Negative, NaN, and positive-infinity weights pass through without rejection: a
        /// negative weight computes, a NaN weight never relaxes its edge, and an infinite
        /// weight is impassable.
        /// </summary>
        [TestMethod]
        public void NegativeNaNAndInfinityWeightsAreAccepted()
        {
            var negative = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, -5, 0) }, 1, 2);
            Assert.AreEqual(-5f, negative[0, 2]);

            var nan = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, float.NaN, 0) }, 1, 2);
            Assert.IsTrue(float.IsPositiveInfinity(nan[0, 2]), "A NaN weight never relaxes, severing its edge.");

            var infinite = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, float.PositiveInfinity, 0) }, 1, 2);
            Assert.IsTrue(float.IsPositiveInfinity(infinite[0, 2]), "An infinite weight is impassable.");
        }

        /// <summary>
        /// An empty edge list with an explicit node count yields the destination row and
        /// all-unreachable rows, and an empty destination array yields an all-unreachable table.
        /// </summary>
        [TestMethod]
        public void EmptyInputsPreserveDocumentedSemantics()
        {
            var noEdges = Dijkstra.Solve(new List<Edge>(), 0, 3);
            Assert.AreEqual(0f, noEdges[0, 0]);
            Assert.AreEqual(0f, noEdges[0, 2]);
            Assert.IsTrue(float.IsPositiveInfinity(noEdges[1, 2]));
            Assert.IsTrue(float.IsPositiveInfinity(noEdges[2, 2]));

            var noDestinations = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, 1, 0) }, new int[0], 2);
            Assert.IsTrue(float.IsPositiveInfinity(noDestinations[0, 2]));
            Assert.IsTrue(float.IsPositiveInfinity(noDestinations[1, 2]));
        }

        /// <summary>
        /// The path walker returns the ordered edge indices to the destination, and the try
        /// variant returns the total cost read from the table.
        /// </summary>
        [TestMethod]
        public void GetPathWalksEdgeIndices()
        {
            var edges = new List<Edge>
            {
                new Edge(0, 1, 2, 0),
                new Edge(0, 2, 4, 2),
                new Edge(1, 2, 1, 2),
                new Edge(1, 3, 7, 3),
                new Edge(2, 3, 3, 4),
                new Edge(4, 0, 1, 5),
            };
            var result = Dijkstra.Solve(edges, 3, 6);

            var path = Dijkstra.GetPath(result, 0);
            Assert.IsNotNull(path);
            CollectionAssert.AreEqual(new List<int> { 0, 2, 4 }, path);

            Assert.IsTrue(Dijkstra.TryGetPath(result, 4, out var fromFour, out float cost));
            CollectionAssert.AreEqual(new List<int> { 5, 0, 2, 4 }, fromFour);
            Assert.AreEqual(7f, cost, 0f);
        }

        /// <summary>
        /// Starting the walk at the destination returns an empty path with zero cost, not null.
        /// </summary>
        [TestMethod]
        public void GetPathStartAtDestinationReturnsEmpty()
        {
            var result = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, 1, 0) }, 1, 2);

            var path = Dijkstra.GetPath(result, 1);
            Assert.IsNotNull(path);
            Assert.AreEqual(0, path!.Count);

            Assert.IsTrue(Dijkstra.TryGetPath(result, 1, out var tryPath, out float cost));
            Assert.AreEqual(0, tryPath.Count);
            Assert.AreEqual(0f, cost, 0f);
        }

        /// <summary>
        /// Walking from an unreachable node returns null (or false with an infinite cost).
        /// </summary>
        [TestMethod]
        public void GetPathUnreachableReturnsNull()
        {
            var result = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, 1, 0) }, 1, 3);

            Assert.IsNull(Dijkstra.GetPath(result, 2));
            Assert.IsFalse(Dijkstra.TryGetPath(result, 2, out var path, out float cost));
            Assert.AreEqual(0, path.Count);
            Assert.IsTrue(float.IsPositiveInfinity(cost));
        }

        /// <summary>
        /// A cyclic or non-converging table fails the walk loudly instead of hanging.
        /// </summary>
        [TestMethod]
        public void GetPathThrowsOnInconsistentTable()
        {
            var cyclic = new float[2, 3];
            cyclic[0, 0] = 1; cyclic[0, 1] = 0; cyclic[0, 2] = 1;
            cyclic[1, 0] = 0; cyclic[1, 1] = 1; cyclic[1, 2] = 1;

            Assert.Throws<ArgumentException>(() => Dijkstra.GetPath(cyclic, 0));
        }

        /// <summary>
        /// The single-pass nearest-destination solve matches the merged multi-destination solve:
        /// costs are identical on random graphs, and every column is identical on a tie-free
        /// graph.
        /// </summary>
        [TestMethod]
        public void SolveNearestMatchesMergedSolve()
        {
            var edges = BuildRandomGraph(new MersenneTwister(106), 400, 1600);
            int[] destinations = { 11, 222, 333 };
            var merged = Dijkstra.Solve(edges, destinations, 400);
            var nearest = Dijkstra.SolveNearest(edges, destinations, 400);
            for (int i = 0; i < 400; i++)
            {
                Assert.AreEqual(merged[i, 2], nearest[i, 2], 0f, $"Cost mismatch at node {i}.");
            }

            // Tie-free line graph: 0 -1- 1 -1- 2 -10- 3 -1- 4; destinations 0 and 4.
            var line = new List<Edge>
            {
                new Edge(1, 0, 1, 0),
                new Edge(2, 1, 1, 1),
                new Edge(3, 2, 10, 2),
                new Edge(3, 4, 1, 3),
                new Edge(2, 3, 10, 4),
                new Edge(1, 2, 1, 5),
                new Edge(0, 1, 1, 6),
                new Edge(4, 3, 1, 7),
            };
            var mergedLine = Dijkstra.Solve(line, [0, 4], 5);
            var nearestLine = Dijkstra.SolveNearest(line, [0, 4], 5);
            for (int i = 0; i < 5; i++)
            {
                for (int c = 0; c < 3; c++)
                {
                    Assert.AreEqual(mergedLine[i, c], nearestLine[i, c], 0f, $"Cell [{i},{c}] differs on the tie-free graph.");
                }
            }
        }

        /// <summary>
        /// The nearest-destination solve is deterministic across repeated calls and tolerates
        /// duplicate destination indices.
        /// </summary>
        [TestMethod]
        public void SolveNearestIsDeterministic()
        {
            var edges = BuildRandomGraph(new MersenneTwister(107), 250, 1000);

            var first = Dijkstra.SolveNearest(edges, [5, 5, 100], 250);
            var second = Dijkstra.SolveNearest(edges, [5, 100], 250);
            for (int i = 0; i < 250; i++)
            {
                for (int c = 0; c < 3; c++)
                {
                    Assert.AreEqual(first[i, c], second[i, c], 0f, $"Cell [{i},{c}] differs with duplicate destinations.");
                }
            }
        }

        /// <summary>
        /// The nearest-destination solve rejects null and empty destination arrays and
        /// out-of-range destinations.
        /// </summary>
        [TestMethod]
        public void SolveNearestValidationThrows()
        {
            var edges = new List<Edge> { new Edge(0, 1, 1, 0) };

            Assert.Throws<ArgumentNullException>(() => Dijkstra.SolveNearest(null!, [0], 2));
            Assert.Throws<ArgumentNullException>(() => Dijkstra.SolveNearest(edges, null!, 2));
            Assert.Throws<ArgumentException>(() => Dijkstra.SolveNearest(edges, new int[0], 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => Dijkstra.SolveNearest(edges, [7], 2));
        }

        /// <summary>
        /// Performance smoke: a 316 by 316 four-neighbor grid (99,856 nodes) solves to a corner
        /// destination without heap overflow, and the far corner's cost equals the Manhattan
        /// distance.
        /// </summary>
        [TestMethod]
        [TestCategory("Performance")]
        public void P1_LargeGridSingleDestination()
        {
            const int side = 316;
            var edges = BuildGrid(side);

            var stopwatch = Stopwatch.StartNew();
            var result = Dijkstra.Solve(edges, 0, side * side);
            stopwatch.Stop();
            Console.WriteLine($"P1 grid {side}x{side}: {stopwatch.Elapsed.TotalSeconds:F3} s");

            Assert.AreEqual(2f * (side - 1), result[side * side - 1, 2], 0f);
        }

        /// <summary>
        /// Performance smoke: the same grid solved to its four corners through the
        /// multi-destination overload; the center's merged cost equals the nearest corner's
        /// Manhattan distance.
        /// </summary>
        [TestMethod]
        [TestCategory("Performance")]
        public void P2_LargeGridMultiDestination()
        {
            const int side = 316;
            var edges = BuildGrid(side);
            int[] corners = { 0, side - 1, side * (side - 1), side * side - 1 };

            var stopwatch = Stopwatch.StartNew();
            var result = Dijkstra.Solve(edges, corners, side * side);
            stopwatch.Stop();
            Console.WriteLine($"P2 grid {side}x{side} x4 destinations: {stopwatch.Elapsed.TotalSeconds:F3} s");

            stopwatch.Restart();
            var nearest = Dijkstra.SolveNearest(edges, corners, side * side);
            stopwatch.Stop();
            Console.WriteLine($"P2 grid {side}x{side} SolveNearest: {stopwatch.Elapsed.TotalSeconds:F3} s");

            int center = (side / 2) * side + (side / 2);
            // The nearest corner from (158, 158) is (315, 315): 157 + 157 steps.
            float expected = 2f * (side - 1 - side / 2);
            Assert.AreEqual(expected, result[center, 2], 0f);
            Assert.AreEqual(expected, nearest[center, 2], 0f);
        }

        /// <summary>
        /// Builds a random directed graph with integer-valued float weights (exact float sums).
        /// </summary>
        /// <param name="randy">The seeded generator.</param>
        /// <param name="nodeCount">The node count.</param>
        /// <param name="edgeCount">The edge count.</param>
        /// <returns>The edge list.</returns>
        private static List<Edge> BuildRandomGraph(Random randy, int nodeCount, int edgeCount)
        {
            var edges = new List<Edge>(edgeCount);
            for (int j = 0; j < edgeCount; j++)
            {
                edges.Add(new Edge(randy.Next(0, nodeCount), randy.Next(0, nodeCount), randy.Next(1, 20), j));
            }
            return edges;
        }

        /// <summary>
        /// Builds a bidirectional four-neighbor grid with unit weights.
        /// </summary>
        /// <param name="side">The grid side length.</param>
        /// <returns>The edge list.</returns>
        private static List<Edge> BuildGrid(int side)
        {
            var edges = new List<Edge>(4 * side * (side - 1));
            int index = 0;
            for (int r = 0; r < side; r++)
            {
                for (int c = 0; c < side; c++)
                {
                    int node = r * side + c;
                    if (c + 1 < side)
                    {
                        edges.Add(new Edge(node, node + 1, 1, index++));
                        edges.Add(new Edge(node + 1, node, 1, index++));
                    }
                    if (r + 1 < side)
                    {
                        edges.Add(new Edge(node, node + side, 1, index++));
                        edges.Add(new Edge(node + side, node, 1, index++));
                    }
                }
            }
            return edges;
        }

        /// <summary>
        /// Independent Bellman-Ford re-derivation of every node's cost to the destination,
        /// following edges backward exactly like the solver under test.
        /// </summary>
        /// <param name="edges">The network edges.</param>
        /// <param name="nodeCount">The node count.</param>
        /// <param name="destination">The destination node.</param>
        /// <returns>The per-node costs to the destination.</returns>
        private static float[] BellmanFordToDestination(List<Edge> edges, int nodeCount, int destination)
        {
            var dist = new float[nodeCount];
            for (int i = 0; i < nodeCount; i++) dist[i] = float.PositiveInfinity;
            dist[destination] = 0f;

            for (int pass = 0; pass < nodeCount; pass++)
            {
                bool changed = false;
                for (int j = 0; j < edges.Count; j++)
                {
                    var edge = edges[j];
                    float candidate = dist[edge.ToIndex] + edge.Weight;
                    if (candidate < dist[edge.FromIndex])
                    {
                        dist[edge.FromIndex] = candidate;
                        changed = true;
                    }
                }
                if (!changed) break;
            }
            return dist;
        }
    }
}
