using System;
using System.Collections.Generic;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Tests the compiled-topology network wrapper over the Dijkstra solvers.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_Network
    {
        /// <summary>
        /// Builds the two-way five-node test network used throughout: a short chain 0-1-2 with a
        /// slow bypass 0-3-2 and a spur 2-4.
        /// </summary>
        /// <returns>The edges and the network destinations (node 4).</returns>
        private static (Edge[] Edges, int[] Destinations) BuildFixture()
        {
            var edges = new[]
            {
                new Edge(0, 1, 1, 0),
                new Edge(1, 0, 1, 0),
                new Edge(1, 2, 1, 1),
                new Edge(2, 1, 1, 1),
                new Edge(0, 3, 5, 2),
                new Edge(3, 0, 5, 2),
                new Edge(3, 2, 5, 3),
                new Edge(2, 3, 5, 3),
                new Edge(2, 4, 1, 4),
                new Edge(4, 2, 1, 4),
            };
            return (edges, new[] { 4 });
        }

        /// <summary>
        /// Construction succeeds for a non-empty edge set, groups the adjacency correctly with
        /// empty lists for edge-free directions, and reports the node count.
        /// </summary>
        [TestMethod]
        public void CtorBuildsAdjacencyForNonEmptyEdges()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            Assert.AreEqual(5, network.NodeCount);
            Assert.AreEqual(5, network.IncomingEdges.Length);
            Assert.AreEqual(5, network.OutgoingEdges.Length);
            Assert.AreEqual(2, network.OutgoingEdges[0].Count);   // 0->1 and 0->3
            Assert.AreEqual(1, network.IncomingEdges[4].Count);
            Assert.AreEqual(2, network.IncomingEdges[4][0].FromIndex);
            CollectionAssert.AreEqual(new[] { 4 }, network.DestinationIndices);

            // An isolated node (index 5 present only via the node count) is impossible here;
            // instead verify a node with no incoming edges gets an empty list, not null.
            var oneWay = new Network(new[] { new Edge(0, 1, 1, 0) }, new[] { 1 });
            Assert.IsNotNull(oneWay.IncomingEdges[0]);
            Assert.AreEqual(0, oneWay.IncomingEdges[0].Count);
            Assert.IsNotNull(oneWay.OutgoingEdges[1]);
            Assert.AreEqual(0, oneWay.OutgoingEdges[1].Count);
        }

        /// <summary>
        /// The constructor copies both input arrays: mutating them afterward changes nothing.
        /// </summary>
        [TestMethod]
        public void CtorCopiesInputs()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);
            var before = network.Solve(4);

            edges[0] = new Edge(0, 1, 100, 0);
            destinations[0] = 0;

            var after = network.Solve(4);
            Assert.AreEqual(before[0, 2], after[0, 2], 0f);
            CollectionAssert.AreEqual(new[] { 4 }, network.DestinationIndices);
        }

        /// <summary>
        /// Malformed construction inputs fail with clear argument exceptions.
        /// </summary>
        [TestMethod]
        public void CtorValidationThrows()
        {
            var edges = new[] { new Edge(0, 1, 1, 0) };

            Assert.Throws<ArgumentNullException>(() => new Network(null!, [0]));
            Assert.Throws<ArgumentNullException>(() => new Network(edges, null!));
            Assert.Throws<ArgumentException>(() => new Network(new Edge[0], [0]));
            Assert.Throws<ArgumentException>(() => new Network(edges, new int[0]));
            Assert.Throws<ArgumentOutOfRangeException>(() => new Network(edges, [5]));
            Assert.Throws<ArgumentException>(() => new Network(new[] { new Edge(-1, 1, 1, 0) }, [0]));
        }

        /// <summary>
        /// The network solves equal the static Dijkstra solves on the same inputs, cell for
        /// cell, for single and multiple destinations.
        /// </summary>
        [TestMethod]
        public void SolveMatchesDijkstra()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            var networkSingle = network.Solve(4);
            var directSingle = Dijkstra.Solve(edges, 4, network.NodeCount);
            var networkMulti = network.Solve(new[] { 0, 4 });
            var directMulti = Dijkstra.Solve(edges, new[] { 0, 4 }, network.NodeCount);

            for (int i = 0; i < network.NodeCount; i++)
            {
                for (int c = 0; c < 3; c++)
                {
                    Assert.AreEqual(directSingle[i, c], networkSingle[i, c], 0f, $"Single cell [{i},{c}].");
                    Assert.AreEqual(directMulti[i, c], networkMulti[i, c], 0f, $"Multi cell [{i},{c}].");
                }
            }

            Assert.Throws<ArgumentOutOfRangeException>(() => network.Solve(9));
            Assert.Throws<ArgumentNullException>(() => network.Solve((int[])null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => network.Solve(new[] { 9 }));
        }

        /// <summary>
        /// Custom edge weights are actually applied: raising the short chain's weights reroutes
        /// traffic over the bypass.
        /// </summary>
        [TestMethod]
        public void SolveCustomWeightsAreApplied()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            // Base network: node 0 routes 0-1-2-4 at cost 3.
            var baseTable = network.Solve(4);
            Assert.AreEqual(3f, baseTable[0, 2], 0f);
            Assert.AreEqual(1f, baseTable[0, 0], 0f);

            // Flood the 0-1 and 1-2 links: node 0 must reroute 0-3-2-4 at cost 11.
            var weights = new float[] { 100, 100, 100, 100, 5, 5, 5, 5, 1, 1 };
            var flooded = network.Solve(weights);
            Assert.AreEqual(11f, flooded[0, 2], 0f);
            Assert.AreEqual(3f, flooded[0, 0], 0f);

            Assert.Throws<ArgumentNullException>(() => network.Solve((float[])null!));
            Assert.Throws<ArgumentException>(() => network.Solve(new float[3]));
        }

        /// <summary>
        /// An infinite custom weight severs its edge: with every route to the destination
        /// flooded to infinity, the start node becomes unreachable.
        /// </summary>
        [TestMethod]
        public void SolveCustomWeightsInfinityBlocksEdge()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);

            float inf = float.PositiveInfinity;
            // Sever 1->2 and 3->2 (positions 2 and 6): node 0 cannot reach node 4.
            var weights = new float[] { 1, 1, inf, 1, 5, 5, inf, 5, 1, 1 };
            var severed = network.Solve(weights);

            Assert.IsTrue(float.IsPositiveInfinity(severed[0, 2]));
            Assert.IsTrue(float.IsPositiveInfinity(severed[1, 2]));
            Assert.AreEqual(1f, severed[2, 2], 0f);
        }

        /// <summary>
        /// The table-reuse overloads produce cell-identical results to the allocating solves,
        /// repeated reuse into the same buffer stays identical, and wrong dimensions throw.
        /// </summary>
        [TestMethod]
        public void SolveReuseOverloadIsBitIdentical()
        {
            var (edges, destinations) = BuildFixture();
            var network = new Network(edges, destinations);
            var weights = new float[] { 2, 2, 3, 3, 4, 4, 6, 6, 1, 1 };

            var allocated = network.Solve(weights);
            var reused = new float[network.NodeCount, 3];
            network.Solve(weights, reused);
            network.Solve(weights, reused); // steady-state second pass into the same buffer

            var nearestAllocated = network.SolveNearest();
            var nearestReused = new float[network.NodeCount, 3];
            var baseWeights = new float[] { 1, 1, 1, 1, 5, 5, 5, 5, 1, 1 };
            network.SolveNearest(baseWeights, nearestReused);

            for (int i = 0; i < network.NodeCount; i++)
            {
                for (int c = 0; c < 3; c++)
                {
                    Assert.AreEqual(allocated[i, c], reused[i, c], 0f, $"Reuse cell [{i},{c}].");
                    Assert.AreEqual(nearestAllocated[i, c], nearestReused[i, c], 0f, $"Nearest reuse cell [{i},{c}].");
                }
            }

            Assert.Throws<ArgumentNullException>(() => network.Solve(weights, null!));
            Assert.Throws<ArgumentException>(() => network.Solve(weights, new float[2, 3]));
            Assert.Throws<ArgumentException>(() => network.SolveNearest(weights, new float[network.NodeCount, 2]));
        }

        /// <summary>
        /// The single-pass nearest solve matches the merged solve costs over the network's
        /// destinations.
        /// </summary>
        [TestMethod]
        public void SolveNearestMatchesMergedNetworkSolve()
        {
            var (edges, _) = BuildFixture();
            var network = new Network(edges, new[] { 0, 4 });

            var merged = network.Solve(new[] { 0, 4 });
            var nearest = network.SolveNearest();

            for (int i = 0; i < network.NodeCount; i++)
            {
                Assert.AreEqual(merged[i, 2], nearest[i, 2], 0f, $"Cost mismatch at node {i}.");
            }
        }

        /// <summary>
        /// Performance smoke: one thousand custom-weight re-solves through the table-reuse
        /// overload on a 100 by 100 grid.
        /// </summary>
        [TestMethod]
        [TestCategory("Performance")]
        public void P3_RepeatedWeightResolves()
        {
            const int side = 100;
            var edges = new List<Edge>();
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
            var network = new Network(edges.ToArray(), new[] { 0, side * side - 1 });
            var weights = new float[edges.Count];
            for (int j = 0; j < weights.Length; j++) weights[j] = 1f;
            var table = new float[network.NodeCount, 3];

            var stopwatch = Stopwatch.StartNew();
            for (int pass = 0; pass < 1000; pass++)
            {
                weights[pass % weights.Length] = 1f + (pass % 5);
                network.SolveNearest(weights, table);
            }
            stopwatch.Stop();
            Console.WriteLine($"P3 1000 re-solves on {side}x{side}: {stopwatch.Elapsed.TotalSeconds:F3} s ({stopwatch.Elapsed.TotalMilliseconds / 1000.0:F3} ms per solve)");

            Assert.AreEqual(0f, table[0, 2], 0f);
        }
    }
}
