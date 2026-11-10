
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
            Assert.Throws<ArgumentException>(() => Dijkstra.Solve(new List<Edge> { new Edge(0, 1, 1, -1) }, 1, 2));

            var incomingEdges = new List<Edge>[]
            {
                new(),
                new() { new Edge(0, 1, 1, -1) }
            };
            Assert.Throws<ArgumentException>(() => Dijkstra.Solve(edges, 1, 2, incomingEdges));

            Assert.Throws<ArgumentNullException>(() => Dijkstra.PathExists(null!, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Dijkstra.PathExists(new float[2, 3], 2));
            Assert.Throws<ArgumentException>(() => Dijkstra.PathExists(new float[2, 2], 0));
        }

        /// <summary>
        /// NaN and positive-infinity weights pass through without rejection: a NaN weight never relaxes its edge, and
        /// an infinite weight is impassable. A negative zero weight is not negative and is accepted.
        /// </summary>
        [TestMethod]
        public void NaNAndInfinityWeightsAreAccepted()
        {
            var negativeZero = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, -0f, 0) }, 1, 2);
            Assert.AreEqual(0f, negativeZero[0, 2]);

            var nan = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, float.NaN, 0) }, 1, 2);
            Assert.IsTrue(float.IsPositiveInfinity(nan[0, 2]), "A NaN weight never relaxes, severing its edge.");

            var infinite = Dijkstra.Solve(new List<Edge> { new Edge(0, 1, float.PositiveInfinity, 0) }, 1, 2);
            Assert.IsTrue(float.IsPositiveInfinity(infinite[0, 2]), "An infinite weight is impassable.");
        }

        /// <summary>
        /// A single negative edge weight makes every edge-list solver throw at entry, naming the edge's position,
        /// including when the caller supplies the per-node incoming-edge lists.
        /// </summary>
        [TestMethod]
        public void NegativeWeightThrowsInEveryEdgeListSolver()
        {
            var edges = new List<Edge> { new Edge(0, 1, 2f, 0), new Edge(1, 2, -5f, 1) };
            var ex = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, 2, 3));
            StringAssert.Contains(ex.Message, "The weight of the edge at position 1 must be non-negative.");
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, new[] { 2 }, 3));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.SolveNearest(edges, new[] { 2 }, 3));

            var incoming = new List<Edge>[] { new List<Edge>(), new List<Edge> { edges[0] }, new List<Edge> { edges[1] } };
            ex = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, 2, 3, incoming));
            StringAssert.Contains(ex.Message, "The weight of an incoming edge of node 2 must be non-negative.");
        }

        /// <summary>
        /// A graph with a negative cycle (0 -> 1 -> 0 with total weight -1) throws at entry instead of relaxing
        /// around the cycle, in the static solvers and in every network method that would search with the weights.
        /// The network itself is still constructed, since its methods with custom weights do not use the constructor's.
        /// </summary>
        [TestMethod]
        public void NegativeCycleThrowsQuickly()
        {
            var edges = new[] { new Edge(0, 1, 1f, 0), new Edge(1, 0, -2f, 1), new Edge(1, 2, 1f, 2) };
            var watch = System.Diagnostics.Stopwatch.StartNew();
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, 2, 3));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.Solve(edges, new[] { 0, 2 }, 3));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.SolveNearest(edges, new[] { 2 }, 3));

            var network = new Network(edges, new[] { 2 });
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.Solve(2));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.Solve(new[] { 2 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.SolveNearest());
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.GetPath(new int[0], 0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.GetPath(new int[0], 0, new float[3, 3]));

            var negativeWeights = new[] { 1f, -2f, 1f };
            var table = new float[3, 3];
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.Solve(negativeWeights));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.Solve(negativeWeights, table));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.SolveNearest(negativeWeights, table));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => network.GetPath(new int[0], 0, table, negativeWeights));
            Assert.IsLessThan(5000L, watch.ElapsedMilliseconds, "Every call must throw at entry.");

            // Custom non-negative weights still solve on the same network.
            network.SolveNearest(new[] { 1f, 2f, 1f }, table);
            Assert.AreEqual(1f, table[1, 2]);
            Assert.AreEqual(2f, table[0, 2]);
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
            Assert.IsEmpty(path!);

            Assert.IsTrue(Dijkstra.TryGetPath(result, 1, out var tryPath, out float cost));
            Assert.IsEmpty(tryPath);
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
            Assert.IsEmpty(path);
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

        /// <summary>The decimal-weight network of the scipy.sparse.csgraph.dijkstra oracle cases (weights rounded to float32).</summary>
        private static List<Edge> DecimalWeightEdges()
        {
            var e = new (int From, int To, float W)[]
            {
                (0, 1, 0.1f), (1, 2, 0.2f), (0, 2, 0.3f), (2, 3, 0.7f), (1, 3, 0.9f), (3, 4, 0.1f), (2, 4, 0.8f),
                (4, 5, 0.3f), (0, 5, 1.4f), (5, 6, 0.6f), (6, 7, 0.1f), (3, 7, 1.0f), (7, 0, 2.5f), (8, 7, 0.2f),
            };
            var edges = new List<Edge>();
            for (int i = 0; i < e.Length; i++) edges.Add(new Edge(e[i].From, e[i].To, e[i].W, i));
            return edges;
        }

        /// <summary>
        /// The single-precision tables are float32 of the double-precision distances of
        /// scipy.sparse.csgraph.dijkstra (SciPy 1.18.1) on the reversed graph, bitwise, for one and for
        /// two destinations, through every float signature.
        /// </summary>
        [TestMethod]
        public void FloatTablesMatchScipyDoubleDistances()
        {
            var edges = DecimalWeightEdges();
            float[] toSeven = { (float)2.0, (float)1.899999976158142, (float)1.7000000476837158, (float)1.0, (float)1.0,
                (float)0.7000000476837158, (float)0.10000000149011612, 0f, (float)0.20000000298023224, float.PositiveInfinity };
            float[] toSevenOrFour = { (float)1.100000023841858, (float)1.0, (float)0.800000011920929, (float)0.10000000149011612, 0f,
                (float)0.7000000476837158, (float)0.10000000149011612, 0f, (float)0.20000000298023224, float.PositiveInfinity };

            AssertCosts(toSeven, Dijkstra.Solve(edges, 7, 10));
            AssertCosts(toSeven, new Network(edges.ToArray(), new[] { 7 }).Solve(7));
            AssertCosts(toSevenOrFour, Dijkstra.Solve(edges, new[] { 7, 4 }, 10));
            AssertCosts(toSevenOrFour, Dijkstra.SolveNearest(edges, new[] { 7, 4 }, 10));
            var network = new Network(edges.ToArray(), new[] { 7, 4 });
            AssertCosts(toSevenOrFour, network.SolveNearest());
            var weights = new float[edges.Count];
            for (int i = 0; i < weights.Length; i++) weights[i] = edges[i].Weight;
            AssertCosts(toSevenOrFour, network.Solve(weights));
            var table = new float[network.NodeCount, 3];
            network.SolveNearest(weights, table);
            AssertCosts(toSevenOrFour, table);
        }

        private static void AssertCosts(float[] expected, float[,] table)
        {
            // A network derives its node count from the edges, so its table omits the isolated node 9.
            for (int i = 0; i < table.GetLength(0); i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i]), BitConverter.DoubleToInt64Bits(table[i, 2]), $"node {i}: expected {expected[i]:R}, actual {table[i, 2]:R}");
        }

        /// <summary>
        /// A finite double cost beyond single-precision range is written as an unreachable row, (-1, -1, +inf), so the
        /// table keeps its contract; a cost within range is kept.
        /// </summary>
        [TestMethod]
        public void CostBeyondFloatRangeIsAnUnreachableRow()
        {
            var edges = new List<Edge> { new Edge(0, 1, 3E38f, 0), new Edge(1, 2, 3E38f, 1) };
            float[,] table = Dijkstra.Solve(edges, 2, 3);
            Assert.AreEqual(2f, table[1, 0]);
            Assert.AreEqual(1f, table[1, 1]);
            Assert.AreEqual(3E38f, table[1, 2]);
            Assert.AreEqual(-1f, table[0, 0]);
            Assert.AreEqual(-1f, table[0, 1]);
            Assert.IsTrue(float.IsPositiveInfinity(table[0, 2]));
            Assert.IsFalse(Dijkstra.PathExists(table, 0));
            Assert.IsNull(Dijkstra.GetPath(table, 0));
            float[,] nearest = new Network(edges.ToArray(), new[] { 2 }).SolveNearest();
            Assert.AreEqual(-1f, nearest[0, 0]);
            Assert.IsTrue(float.IsPositiveInfinity(nearest[0, 2]));
        }

        /// <summary>The CSR graph of the scipy.sparse.csgraph.dijkstra oracle cases.</summary>
        private static readonly int[] CsrRowStart = { 0, 2, 4, 5, 7, 7, 8 };
        private static readonly int[] CsrToNode = { 1, 2, 2, 3, 4, 4, 5, 0 };
        private static readonly double[] CsrWeight = { 0.1, 0.3, 0.2, 0.9, 0.7, 0.1, 2.0, 0.4 };

        /// <summary>
        /// Distances, predecessors and nearest sources equal scipy.sparse.csgraph.dijkstra(min_only=True,
        /// return_predecessors=True) (SciPy 1.18.1) bitwise, with -1 where SciPy reports -9999, for one and two sources.
        /// </summary>
        [TestMethod]
        public void ShortestDistancesMatchesScipy()
        {
            var predecessors = new int[6];
            var nearest = new int[6];
            double[] d = Dijkstra.ShortestDistances(6, CsrRowStart, CsrToNode, CsrWeight, new[] { 0 }, double.PositiveInfinity, predecessors, nearest);
            CollectionAssert.AreEqual(new[] { 0.0, 0.1, 0.3, 1.0, 1.0, 3.0 }, d);
            CollectionAssert.AreEqual(new[] { -1, 0, 0, 1, 2, 3 }, predecessors);
            CollectionAssert.AreEqual(new[] { 0, 0, 0, 0, 0, 0 }, nearest);

            d = Dijkstra.ShortestDistances(6, CsrRowStart, CsrToNode, CsrWeight, new List<int> { 0, 3, 3 }, double.PositiveInfinity, predecessors, nearest);
            CollectionAssert.AreEqual(new[] { 0.0, 0.1, 0.3, 0.0, 0.1, 2.0 }, d);
            CollectionAssert.AreEqual(new[] { -1, 0, 0, -1, 3, 3 }, predecessors);
            CollectionAssert.AreEqual(new[] { 0, 0, 0, 3, 3, 3 }, nearest);

            // Without the optional outputs the distances are the same.
            CollectionAssert.AreEqual(d, Dijkstra.ShortestDistances(6, CsrRowStart, CsrToNode, CsrWeight, new[] { 0, 3 }));
        }

        /// <summary>
        /// Tie routing regression on a 4 x 4 grid with unit weights in both directions and sources 0 and 15, where
        /// most nodes have several equal-cost predecessors and the two sources tie along the anti-diagonal. The
        /// predecessors and nearest sources pin the heap's settle order among equal distances.
        /// </summary>
        [TestMethod]
        public void ShortestDistancesTieRoutingIsPinned()
        {
            int[] rowStart = { 0, 2, 5, 8, 10, 13, 17, 21, 24, 27, 31, 35, 38, 40, 43, 46, 48 };
            int[] toNode = { 1, 4, 2, 5, 0, 3, 6, 1, 7, 2, 5, 8, 0, 6, 9, 4, 1, 7, 10, 5, 2, 11, 6, 3,
                9, 12, 4, 10, 13, 8, 5, 11, 14, 9, 6, 15, 10, 7, 13, 8, 14, 12, 9, 15, 13, 10, 14, 11 };
            var weight = new double[toNode.Length];
            for (int k = 0; k < weight.Length; k++) weight[k] = 1d;
            var predecessors = new int[16];
            var nearest = new int[16];
            double[] d = Dijkstra.ShortestDistances(16, rowStart, toNode, weight, new[] { 0, 15 }, double.PositiveInfinity, predecessors, nearest);
            CollectionAssert.AreEqual(new[] { 0.0, 1, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 2, 1, 0 }, d);
            CollectionAssert.AreEqual(new[] { -1, 0, 1, 2, 0, 4, 2, 11, 4, 13, 11, 15, 13, 14, 15, -1 }, predecessors);
            CollectionAssert.AreEqual(new[] { 0, 0, 0, 0, 0, 0, 0, 15, 0, 15, 15, 15, 15, 15, 15, 15 }, nearest);
        }

        /// <summary>
        /// A limit leaves every node beyond it unreached, as scipy.sparse.csgraph.dijkstra(limit=0.35) does, and
        /// an empty source list leaves every node unreached.
        /// </summary>
        [TestMethod]
        public void ShortestDistancesLimitAndEmptySources()
        {
            var predecessors = new int[6];
            var nearest = new int[6];
            double inf = double.PositiveInfinity;
            double[] d = Dijkstra.ShortestDistances(6, CsrRowStart, CsrToNode, CsrWeight, new[] { 0 }, 0.35, predecessors, nearest);
            CollectionAssert.AreEqual(new[] { 0.0, 0.1, 0.3, inf, inf, inf }, d);
            CollectionAssert.AreEqual(new[] { -1, 0, 0, -1, -1, -1 }, predecessors);
            CollectionAssert.AreEqual(new[] { 0, 0, 0, -1, -1, -1 }, nearest);

            d = Dijkstra.ShortestDistances(6, CsrRowStart, CsrToNode, CsrWeight, new int[0], inf, predecessors, nearest);
            CollectionAssert.AreEqual(new[] { inf, inf, inf, inf, inf, inf }, d);
            CollectionAssert.AreEqual(new[] { -1, -1, -1, -1, -1, -1 }, nearest);

            // An infinite weight is impassable.
            d = Dijkstra.ShortestDistances(2, new[] { 0, 1, 1 }, new[] { 1 }, new[] { inf }, new[] { 0 });
            CollectionAssert.AreEqual(new[] { 0.0, inf }, d);
        }

        /// <summary>
        /// On random graphs every finite distance is reached through its predecessor with the exact
        /// double-precision sum, and agrees with an independent Bellman-Ford relaxation to 1E-14 relative.
        /// </summary>
        [TestMethod]
        public void ShortestDistancesIsSelfConsistentOnRandomGraphs()
        {
            var random = new Random(20261108);
            for (int trial = 0; trial < 20; trial++)
            {
                int n = 50, m = 300;
                var from = new int[m];
                var to = new int[m];
                var w = new double[m];
                for (int k = 0; k < m; k++) { from[k] = random.Next(n); to[k] = random.Next(n); w[k] = random.NextDouble(); }
                // Counting sort into CSR by start node.
                var rowStart = new int[n + 1];
                for (int k = 0; k < m; k++) rowStart[from[k] + 1]++;
                for (int i = 0; i < n; i++) rowStart[i + 1] += rowStart[i];
                var cursor = (int[])rowStart.Clone();
                var toNode = new int[m];
                var weight = new double[m];
                for (int k = 0; k < m; k++) { int slot = cursor[from[k]]++; toNode[slot] = to[k]; weight[slot] = w[k]; }

                var sources = new[] { random.Next(n), random.Next(n) };
                var predecessors = new int[n];
                var nearest = new int[n];
                double[] d = Dijkstra.ShortestDistances(n, rowStart, toNode, weight, sources, double.PositiveInfinity, predecessors, nearest);

                var bellman = new double[n];
                for (int i = 0; i < n; i++) bellman[i] = double.PositiveInfinity;
                foreach (int s in sources) bellman[s] = 0d;
                for (int pass = 0; pass < n; pass++)
                    for (int k = 0; k < m; k++)
                        if (bellman[from[k]] + w[k] < bellman[to[k]]) bellman[to[k]] = bellman[from[k]] + w[k];

                for (int v = 0; v < n; v++)
                {
                    if (double.IsPositiveInfinity(bellman[v])) { Assert.IsTrue(double.IsPositiveInfinity(d[v])); Assert.AreEqual(-1, nearest[v]); continue; }
                    Assert.AreEqual(bellman[v], d[v], 1E-14 * Math.Max(1d, bellman[v]));
                    if (predecessors[v] < 0) { Assert.AreEqual(0d, d[v]); Assert.AreEqual(v, nearest[v]); continue; }
                    int u = predecessors[v];
                    double best = double.PositiveInfinity;
                    for (int k = rowStart[u]; k < rowStart[u + 1]; k++) if (toNode[k] == v) best = Math.Min(best, d[u] + weight[k]);
                    Assert.AreEqual(d[v], best, 0d);
                    Assert.AreEqual(nearest[u], nearest[v]);
                }
            }
        }

        /// <summary>Invalid CSR arrays, sources, weights, limits and output arrays throw.</summary>
        [TestMethod]
        public void ShortestDistancesValidationThrows()
        {
            int[] r = CsrRowStart, t = CsrToNode;
            double[] w = CsrWeight;
            var s = new[] { 0 };
            Assert.ThrowsExactly<ArgumentNullException>(() => Dijkstra.ShortestDistances(6, null!, t, w, s));
            Assert.ThrowsExactly<ArgumentNullException>(() => Dijkstra.ShortestDistances(6, r, null!, w, s));
            Assert.ThrowsExactly<ArgumentNullException>(() => Dijkstra.ShortestDistances(6, r, t, null!, s));
            Assert.ThrowsExactly<ArgumentNullException>(() => Dijkstra.ShortestDistances(6, r, t, w, null!));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.ShortestDistances(0, new[] { 0 }, new int[0], new double[0], new int[0]));
            Assert.ThrowsExactly<ArgumentException>(() => Dijkstra.ShortestDistances(5, r, t, w, s));
            Assert.ThrowsExactly<ArgumentException>(() => Dijkstra.ShortestDistances(6, r, t, new double[7], s));
            Assert.ThrowsExactly<ArgumentException>(() => Dijkstra.ShortestDistances(6, new[] { 0, 2, 1, 5, 7, 7, 8 }, t, w, s));
            Assert.ThrowsExactly<ArgumentException>(() => Dijkstra.ShortestDistances(6, new[] { 1, 2, 4, 5, 7, 7, 8 }, t, w, s));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.ShortestDistances(6, r, new[] { 1, 2, 2, 3, 4, 4, 6, 0 }, w, s));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.ShortestDistances(6, r, t, new[] { 0.1, 0.3, 0.2, -0.9, 0.7, 0.1, 2.0, 0.4 }, s));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.ShortestDistances(6, r, t, new[] { 0.1, 0.3, 0.2, double.NaN, 0.7, 0.1, 2.0, 0.4 }, s));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.ShortestDistances(6, r, t, w, new[] { 6 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.ShortestDistances(6, r, t, w, new[] { -1 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => Dijkstra.ShortestDistances(6, r, t, w, s, double.NaN));
            Assert.ThrowsExactly<ArgumentException>(() => Dijkstra.ShortestDistances(6, r, t, w, s, double.PositiveInfinity, new int[5]));
            Assert.ThrowsExactly<ArgumentException>(() => Dijkstra.ShortestDistances(6, r, t, w, s, double.PositiveInfinity, null, new int[7]));
        }
    }
}
